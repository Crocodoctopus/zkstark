use crate::merkle::{self, AuthPath, Hash};
use crate::F;
use num_traits::Pow;

pub struct Proof {
    state: Hash,
    data: Box<[u8]>,
}

impl Proof {
    pub fn new(state: Hash, data: Box<[u8]>) -> Self {
        Self { state, data }
    }

    pub fn verify(&self) {
        use bincode::deserialize_from;
        let mut reader = std::io::BufReader::new(&self.data[..]);

        // Pull elements out of proof
        let f_eval_merkle_root: Hash = deserialize_from(&mut reader).unwrap();

        let alpha0: u32 = deserialize_from(&mut reader).unwrap();
        let alpha1: u32 = deserialize_from(&mut reader).unwrap();
        let alpha2: u32 = deserialize_from(&mut reader).unwrap();
        let cp_eval_merkle_root: Hash = deserialize_from(&mut reader).unwrap();

        let mut betas = vec![0];
        let mut cp_eval_merkle_roots = vec![cp_eval_merkle_root];
        for _ in 0..10 {
            let beta: u32 = deserialize_from(&mut reader).unwrap();
            let fri_eval_merkle_root: Hash = deserialize_from(&mut reader).unwrap();
            betas.push(beta);
            cp_eval_merkle_roots.push(fri_eval_merkle_root);
        }
        let fri_free_term: u32 = deserialize_from(&mut reader).unwrap();

        let test_point: u32 = deserialize_from(&mut reader).unwrap();
        let f_x: (u32, AuthPath) = deserialize_from(&mut reader).unwrap();
        let f_gx: (u32, AuthPath) = deserialize_from(&mut reader).unwrap();
        let f_ggx: (u32, AuthPath) = deserialize_from(&mut reader).unwrap();
        let cp0_x: (u32, AuthPath) = deserialize_from(&mut reader).unwrap();
        let mut fri_layers = vec![];
        for _ in 0..10 {
            let fri_layer: (u32, u32, AuthPath, AuthPath) = deserialize_from(&mut reader).unwrap();
            fri_layers.push(fri_layer);
        }

        // Protocol consts
        let primitive_root = F::generator();
        let generator_g = primitive_root.pow(3145728);
        let generator_h = primitive_root.pow(393216);
        let g: Vec<F> = (0..1024).map(|n| generator_g.pow(n)).collect();
        let h: Vec<F> = (0..8192).map(|n| generator_h.pow(n)).collect();
        let f_domain: Vec<F> = h.iter().map(|n| primitive_root * *n).collect();

        ///////////////////
        // Prove trace

        // Wrap test point to range 0..8176
        let test_point = test_point as usize % (8192 - 16);

        // Verify computation
        {
            let x = f_domain[test_point];
            let f_x = F::from(f_x.0);
            let f_gx = F::from(f_gx.0);
            let f_ggx = F::from(f_ggx.0);

            let p0 = (f_x - 1) / (x - g[0]);
            let p1 = (f_x - 2338775057) / (x - g[1022]);
            let p2 = (f_ggx - f_gx.pow(2) - f_x.pow(2))
                / ((x.pow(1024) - 1) / ((x - g[1021]) * (x - g[1022]) * (x - g[1023])));
            let cp0 = alpha0 * p0 + alpha1 * p1 + alpha2 * p2;

            //
            assert_eq!(cp0.residue(), cp0_x.0);
        }

        // Verify trace
        assert_eq!(
            merkle::compute_root_from_path(f_x.0, test_point, &f_x.1),
            f_eval_merkle_root
        );
        assert_eq!(
            merkle::compute_root_from_path(f_gx.0, test_point + 8, &f_gx.1),
            f_eval_merkle_root
        );
        assert_eq!(
            merkle::compute_root_from_path(f_ggx.0, test_point + 16, &f_ggx.1),
            f_eval_merkle_root
        );
        assert_eq!(
            merkle::compute_root_from_path(cp0_x.0, test_point, &cp0_x.1),
            cp_eval_merkle_root
        );

        ///////////////////
        // Prove FRI layers

        // Verify computation for first 9 layers
        for n in 0..9 {
            // Get cp(x) and cp(-x) for layer n, and cp(x^2) for layer n + 1
            let (cp0_x, cp0_nx, _, _) = fri_layers[n];
            let (cp1_xx, _, _, _) = fri_layers[n + 1];
            let x = f_domain[test_point].pow(2u32.pow(n as u32));

            // NOTE: the tutorial video got this part wrong!!
            // The numerator of g(x^2) is NOT cp(x) - cp(-x), it is cp(x) + cp(-x)
            // cp(x) - cp(-x) will yield something closer to h(x^2) and give bad results
            let g_xx = (F::from(cp0_x) + F::from(cp0_nx)) / F::from(2);
            let h_xx = (F::from(cp0_x) - F::from(cp0_nx)) / (x * 2);
            let calc_cp1_xx = g_xx + betas[n + 1] * h_xx;
            assert_eq!(cp1_xx, calc_cp1_xx.residue());
        }

        // Verify computation for layer 9 and free term
        {
            let (cp9_x, cp9_nx, _, _) = fri_layers[9];
            let cp10_xx = fri_free_term;
            let x = f_domain[test_point].pow(2u32.pow(9 as u32));

            let g_xx = (F::from(cp9_x) + F::from(cp9_nx)) / F::from(2);
            let h_xx = (F::from(cp9_x) - F::from(cp9_nx)) / (x * 2);
            let calc_cp10_xx = g_xx + betas[10] * h_xx;
            assert_eq!(cp10_xx, calc_cp10_xx.residue());
        }

        // Verify auth paths for FRI layers
        for n in 0..10 {
            // Get data for FRI layer n
            let size = 8192 >> n;
            let (cp0_x, cp0_nx, cp0_x_auth_path, cp0_nx_auth_path) = &fri_layers[n];
            let merkle_root = cp_eval_merkle_roots[n];

            //
            assert_eq!(
                merkle::compute_root_from_path(*cp0_x, test_point % size, cp0_x_auth_path),
                merkle_root,
            );
            assert_eq!(
                merkle::compute_root_from_path(
                    *cp0_nx,
                    (test_point + size / 2) % size,
                    cp0_nx_auth_path
                ),
                merkle_root,
            );
        }
    }

    pub fn size(&self) -> usize {
        use std::mem::size_of;
        size_of::<Self>() + self.data.len()
    }
}

/*pub struct Proof {
    // State
    pub init_state: Hash,
    pub final_state: Hash,

    // Commits (in order)
    pub f_eval_merkle_root: Hash,
    pub alpha0: u32,
    pub alpha1: u32,
    pub alpha2: u32,
    pub cp_eval_merkle_root: Hash,
    pub betas: [u32; 10],
    pub fri_eval_merkle_roots: [Hash; 10],
    pub fri_free_term: u32,

    // Decommits (in order)
    pub test_point: u32,
    pub f_x: (u32, AuthPath),
    pub f_gx: (u32, AuthPath),
    pub f_ggx: (u32, AuthPath),
    pub cp0_x: (u32, AuthPath),
    pub fri_layers: Box<[(u32, AuthPath, u32, AuthPath)]>,
}

impl Proof {
    // Verify that a proof is correct. The video series doesn't go into this part,
    // but we should still be able to get something good just based on our understanding
    // of how the data was generated in the first place.
    pub fn verify(&self) {
    }

    pub fn size(&self) -> usize {
        use std::mem::size_of;

        let mut size = 0;

        size += size_of::<Self>();
        size += self.f_x.1.len() * size_of::<Hash>();
        size += self.f_gx.1.len() * size_of::<Hash>();
        size += self.f_ggx.1.len() * size_of::<Hash>();
        size += self.cp0_x.1.len() * size_of::<Hash>();
        size += self.fri_layers.len() * size_of::<(u32, AuthPath, u32, AuthPath)>();
        for data in self.fri_layers.iter() {
            size += data.1.len() * size_of::<Hash>();
            size += data.3.len() * size_of::<Hash>();
        }

        return size;
    }
}
*/
