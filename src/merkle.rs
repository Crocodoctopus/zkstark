use sha2::Digest;
use sha2::Sha256;
use std::ops::Index;

#[derive(Clone, Debug)]
pub struct Merkle(Box<[[u8; 32]]>);

impl Merkle {
    pub fn new(size: usize, data: impl Iterator<Item = u32>) -> Self {
        // Calculate size
        let mut i = size;
        let mut size = size;
        while i != 0 {
            i /= 2;
            size += i;
        }

        // Create output vec
        let mut out = vec![[0u8; 32]; size];

        // First round of hashing
        let mut offset = out.len() / 2;
        out.splice(
            offset..,
            data.map(|v| {
                let mut hasher = Sha256::new();
                hasher.update(v.to_be_bytes());
                hasher.finalize().into()
            }),
        );

        // The rest
        while offset > 0 {
            offset /= 2;
            for it in 0..offset + 1 {
                let index = offset + it;
                let mut hasher = Sha256::new();
                hasher.update(out[2 * index + 1]);
                hasher.update(out[2 * index + 2]);
                out[index] = hasher.finalize().into();
            }
        }

        // Return
        Self(out.into_boxed_slice())
    }

    pub fn trace(&self, mut i: usize) -> Box<[[u8; 32]]> {
        let mut v = vec![];
        i += self.0.len() / 2;

        while i != 0 {
            // If right
            if i % 2 == 0 {
                v.push(self[i - 1]);
                i -= 2;
            } else {
                v.push(self[i + 1]);
                i -= 1;
            }
            i >>= 1;
        }

        v.into_boxed_slice()
    }
}

impl Index<usize> for Merkle {
    type Output = [u8; 32];
    fn index(&self, i: usize) -> &Self::Output {
        &self.0[i]
    }
}

pub fn compute_root_from_path(element: u32, mut index: usize, path: &Box<[[u8; 32]]>) -> [u8; 32] {
    // BITS
    index += (1 << path.len()) - 1;

    // Generate current hash
    let mut hasher = Sha256::new();
    hasher.update(element.to_be_bytes());
    let mut current = hasher.finalize().into();

    // Step through the path
    for hash in path.iter() {
        // If index is a right node
        if index % 2 == 0 {
            let mut hasher = Sha256::new();
            hasher.update(hash);
            hasher.update(current);
            current = hasher.finalize().into();
            index -= 2;
        } else {
            let mut hasher = Sha256::new();
            hasher.update(current);
            hasher.update(hash);
            current = hasher.finalize().into();
            index -= 1;
        }
        index >>= 1;
    }

    // Return final hash
    return current;
}

#[test]
fn merkle_test() {
    let merkle = Merkle::new(4, [0x01, 0x02, 0x03, 0x04].into_iter());

    // Level 1:
    let i3 = [
        0xb4, 0x07, 0x11, 0xa8, 0x8c, 0x70, 0x39, 0x75, 0x6f, 0xb8, 0xa7, 0x38, 0x27, 0xea, 0xbe,
        0x2c, 0x0f, 0xe5, 0xa0, 0x34, 0x6c, 0xa7, 0xe0, 0xa1, 0x04, 0xad, 0xc0, 0xfc, 0x76, 0x4f,
        0x52, 0x8d,
    ];
    let i4 = [
        0x43, 0x3e, 0xbf, 0x5b, 0xc0, 0x3d, 0xff, 0xa3, 0x85, 0x36, 0x67, 0x32, 0x07, 0xa2, 0x12,
        0x81, 0x61, 0x2c, 0xef, 0x5f, 0xaa, 0x9b, 0xc7, 0xa4, 0xd5, 0xb9, 0xbe, 0x2f, 0xdb, 0x12,
        0xcf, 0x1a,
    ];
    let i5 = [
        0x88, 0x18, 0x5d, 0x12, 0x8d, 0x99, 0x22, 0xe0, 0xe6, 0xbc, 0xd3, 0x2b, 0x07, 0xb6, 0xc7,
        0xf2, 0x0f, 0x27, 0x96, 0x8e, 0xab, 0x44, 0x7a, 0x1d, 0x8d, 0x1c, 0xdf, 0x25, 0x0f, 0x79,
        0xf7, 0xd3,
    ];
    let i6 = [
        0x1b, 0xc5, 0xd0, 0xe3, 0xdf, 0x0e, 0xa1, 0x2c, 0x4d, 0x00, 0x78, 0x66, 0x8d, 0x14, 0x92,
        0x4f, 0x95, 0x10, 0x6b, 0xbe, 0x17, 0x3e, 0x19, 0x6d, 0xe5, 0x0f, 0xe1, 0x3a, 0x90, 0x0b,
        0x09, 0x37,
    ];
    // Level 2:
    let i1 = [
        0xbe, 0x8d, 0xc3, 0x57, 0xde, 0xcb, 0x6e, 0x09, 0xc8, 0xe5, 0xad, 0x87, 0x4d, 0x3c, 0x4f,
        0xa7, 0xfc, 0x09, 0x73, 0x0b, 0xbb, 0x5e, 0x90, 0xf4, 0x2c, 0x97, 0xda, 0xd2, 0x0e, 0x00,
        0x12, 0xd4,
    ];
    let i2 = [
        0x6b, 0xed, 0x5b, 0x6d, 0x7a, 0xe0, 0x93, 0xd1, 0x81, 0x2a, 0xb9, 0xbe, 0x5c, 0xbf, 0xa1,
        0xce, 0x78, 0x78, 0x12, 0xa0, 0x03, 0xd9, 0x5c, 0x11, 0x44, 0x87, 0x20, 0xa4, 0x07, 0xb6,
        0x17, 0x27,
    ];
    // Level 3:
    let i0 = [
        0x32, 0x7c, 0xf2, 0x13, 0xe1, 0x73, 0x8d, 0xe4, 0x20, 0x6b, 0xfd, 0x14, 0x29, 0x7c, 0x26,
        0xc6, 0x82, 0x96, 0x17, 0x50, 0xcb, 0x56, 0x89, 0x7e, 0xd5, 0xe8, 0xf5, 0x19, 0xb0, 0x54,
        0x8f, 0xf2,
    ];

    // Assert merkle elements
    assert_eq!(merkle[0], i0);
    assert_eq!(merkle[1], i1);
    assert_eq!(merkle[2], i2);
    assert_eq!(merkle[3], i3);
    assert_eq!(merkle[4], i4);
    assert_eq!(merkle[5], i5);
    assert_eq!(merkle[6], i6);

    // Trace test on merkle
    let trace0 = merkle.trace(0);
    let trace1 = merkle.trace(1);
    let trace2 = merkle.trace(2);
    let trace3 = merkle.trace(3);

    // Assert traces
    assert_eq!(trace0[0], i4);
    assert_eq!(trace0[1], i2);
    assert_eq!(trace1[0], i3);
    assert_eq!(trace1[1], i2);
    assert_eq!(trace2[0], i6);
    assert_eq!(trace2[1], i1);
    assert_eq!(trace3[0], i5);
    assert_eq!(trace3[1], i1);

    // Assert compute
    assert_eq!(compute_root_from_path(0x01, 0, &trace0), merkle[0]);
}
