use crate::merkle;
use crate::F;
use num_traits::Pow;

#[derive(PartialEq, Debug)]
pub enum ProofErr {
    //
    CpMismatch(u32),
    AuthPathErr(&'static str),
}

pub struct Proof {
    // State
    pub init_state: [u8; 32],
    pub final_state: [u8; 32],

    // Commits (in order)
    pub f_eval_merkle_root: [u8; 32],
    pub alpha0: u32,
    pub alpha1: u32,
    pub alpha2: u32,
    pub cp_eval_merkle_root: [u8; 32],
    pub betas: [u32; 10],
    pub fri_eval_merkle_roots: [[u8; 32]; 10],
    pub fri_free_term: u32,

    // Decommits (in order)
    pub test_point: u32,
    pub f_x: (u32, Box<[[u8; 32]]>),
    pub f_gx: (u32, Box<[[u8; 32]]>),
    pub f_ggx: (u32, Box<[[u8; 32]]>),
    pub cp0_x: (u32, Box<[[u8; 32]]>),
    pub fri_layers: Box<[(u32, Box<[[u8; 32]]>, u32, Box<[[u8; 32]]>)]>,
}

impl Proof {
    pub fn verify(&self) -> Result<(), ProofErr> {
        // Protocol consts
        let primitive_root = F::generator();
        let generator_g = primitive_root.pow(3145728);
        let generator_h = primitive_root.pow(393216);
        let g: Vec<F> = (0..1024).map(|n| generator_g.pow(n)).collect();
        let h: Vec<F> = (0..8192).map(|n| generator_h.pow(n)).collect();
        let f_domain: Vec<F> = h.iter().map(|n| primitive_root * *n).collect();

        ///////////////////
        // Prove trace

        // Verify computation
        let x = f_domain[self.test_point as usize];
        let f_x = F::from(self.f_x.0);
        let f_gx = F::from(self.f_gx.0);
        let f_ggx = F::from(self.f_ggx.0);

        let p0 = (f_x - 1) / (x - g[0]);
        let p1 = (f_x - 2338775057) / (x - g[1022]);
        let p2 = (f_ggx - f_gx.pow(2) - f_x.pow(2))
            / ((x.pow(1024) - 1) / ((x - g[1021]) * (x - g[1022]) * (x - g[1023])));
        let cp0 = self.alpha0 * p0 + self.alpha1 * p1 + self.alpha2 * p2;

        //
        if cp0.residue() != self.cp0_x.0 {
            return Err(ProofErr::CpMismatch(cp0.residue()));
        }

        // Verify trace
        if merkle::compute_root_from_path(self.f_x.0, self.test_point as _, &self.f_x.1)
            != self.f_eval_merkle_root
        {
            return Err(ProofErr::AuthPathErr(
                "Auth path for f(x) doesn't evaluate to merkle root.",
            ));
        }
        if merkle::compute_root_from_path(self.f_gx.0, self.test_point as usize + 8, &self.f_gx.1)
            != self.f_eval_merkle_root
        {
            return Err(ProofErr::AuthPathErr(
                "Auth path for f(gx) doesn't evalutate to merkle root.",
            ));
        }
        if merkle::compute_root_from_path(
            self.f_ggx.0,
            self.test_point as usize + 16,
            &self.f_ggx.1,
        ) != self.f_eval_merkle_root
        {
            return Err(ProofErr::AuthPathErr(
                "Auth path for f(ggx) doesn't evalutate to merkle root.",
            ));
        }
        if merkle::compute_root_from_path(self.cp0_x.0, self.test_point as _, &self.cp0_x.1)
            != self.cp_eval_merkle_root
        {
            return Err(ProofErr::AuthPathErr(
                "Auth path for cp0(x) doesn't evalutate to merkle root.",
            ));
        }

        ///////////////////
        // Prove FRI layers

        // Verify computation for first 9 layers
        for n in 0..9 {
            // Get cp(x) and cp(-x) for layer n, and cp(x^2) for layer n + 1
            let (cp0_x, _, cp0_nx, _) = self.fri_layers[n];
            let (cp1_xx, _, _, _) = self.fri_layers[n + 1];
            let x = f_domain[self.test_point as usize].pow(2u32.pow(n as u32));

            let g_xx = (F::from(cp0_x) + F::from(cp0_nx)) / F::from(2);
            let h_xx = (F::from(cp0_x) - F::from(cp0_nx)) / (x * 2);
            let calc_cp1_xx = g_xx + self.betas[n] * h_xx;
            assert_eq!(cp1_xx, calc_cp1_xx.residue());
        }

        // Verify computation for layer 9 and free term
        let (cp9_x, _, cp9_nx, _) = self.fri_layers[9];
        let cp10_xx = self.fri_free_term;
        let x = f_domain[self.test_point as usize].pow(2u32.pow(9 as u32));

        let g_xx = (F::from(cp9_x) + F::from(cp9_nx)) / F::from(2);
        let h_xx = (F::from(cp9_x) - F::from(cp9_nx)) / (x * 2);
        let calc_cp10_xx = g_xx + self.betas[9] * h_xx;
        assert_eq!(cp10_xx, calc_cp10_xx.residue());

        // Verify auth path for first layer
        let (cp0_x, cp0_x_auth_path, cp0_nx, cp0_nx_auth_path) = &self.fri_layers[0];
        assert_eq!(
            merkle::compute_root_from_path(*cp0_x, self.test_point as usize, cp0_x_auth_path),
            self.cp_eval_merkle_root,
        );

        // Verify auth paths for last 9 layers
        let mut size = 8192;
        for n in 1..10 {
            let (cp0_x, cp0_x_auth_path, cp0_nx, cp0_nx_auth_path) = &self.fri_layers[n];
            size >>= 1;
            assert_eq!(
                merkle::compute_root_from_path(
                    *cp0_x,
                    self.test_point as usize % size,
                    cp0_x_auth_path
                ),
                self.fri_eval_merkle_roots[n - 1],
            );
        }

        // Proof success
        Ok(())
    }

    pub fn size(&self) -> usize {
        use std::mem::size_of;

        let mut size = 0;

        size += size_of::<Self>();
        size += self.f_x.1.len() * size_of::<[u8; 32]>();
        size += self.f_gx.1.len() * size_of::<[u8; 32]>();
        size += self.f_ggx.1.len() * size_of::<[u8; 32]>();
        size += self.cp0_x.1.len() * size_of::<[u8; 32]>();
        size += self.fri_layers.len() * size_of::<(u32, Box<[[u8; 32]]>, u32, Box<[[u8; 32]]>)>();
        for data in self.fri_layers.iter() {
            size += data.1.len() * size_of::<[u8; 32]>();
            size += data.3.len() * size_of::<[u8; 32]>();
        }

        return size;
    }
}
