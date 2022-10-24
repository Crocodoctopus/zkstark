use crate::F;
use num_traits::Pow;

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
    pub fx: (u32, Box<[[u8; 32]]>),
    pub fgx: (u32, Box<[[u8; 32]]>),
    pub fggx: (u32, Box<[[u8; 32]]>),
    pub cp0x: (u32, Box<[[u8; 32]]>),
}

impl Proof {
    pub fn verify(&self) -> bool {
        // Protocol consts
        let primitive_root = F::generator();
        let generator_g = primitive_root.pow(3145728);
        let generator_h = primitive_root.pow(393216);
        let g: Vec<F> = (0..1024).map(|n| generator_g.pow(n)).collect();
        let h: Vec<F> = (0..8192).map(|n| generator_h.pow(n)).collect();
        let f_domain: Vec<F> = h.iter().map(|n| primitive_root * *n).collect();

        // Prove trace
        let x = f_domain[self.test_point as usize];
        let fx = self.fx.0;
        let fgx = self.fgx.0;
        let fggx = self.fggx.0;

        let p0 = (F::from(fx) - 1) / (x - g[0]);
        let p1 = (F::from(fx) - 2338775057) / (x - g[1022]);
        let p2 = (F::from(fggx) - F::from(fgx).pow(2) - F::from(fx).pow(2))
            / ((x.pow(1024) - 1) / ((x - g[1021]) * (x - g[1022]) * (x - g[1023])));
        let cp0 = self.alpha0 * p0 + self.alpha1 * p1 + self.alpha2 * p2;

        //
        if cp0.residue() != self.cp0x.0 {
            return false;
        }

        true
    }
}
