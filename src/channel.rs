use crate::merkle::Hash;
use crate::proof::Proof;
use serde::ser::Serialize;
use sha2::{Digest, Sha256};

pub struct Channel {
    state: Hash,
    data: Vec<u8>,
}

impl Channel {
    pub fn new() -> Self {
        Channel {
            state: [0; 32],
            data: Vec::new(),
        }
    }

    pub fn commit(&mut self, data: impl Serialize) {
        let mut v = bincode::serialize(&data).unwrap();
        let mut hasher = Sha256::new();
        hasher.update(&self.state);
        hasher.update(&v);
        self.state = hasher.finalize().into();
        self.data.append(&mut v);
    }

    pub fn get_u32(&mut self) -> u32 {
        let f = u32::from_be_bytes([self.state[0], self.state[1], self.state[2], self.state[3]]);
        self.commit(f);
        f
    }

    pub fn finalize(self) -> Proof {
        Proof::new(self.state, self.data.into_boxed_slice())
    }
}
