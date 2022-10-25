use crate::field::Gf;
use crate::merkle::{AuthPath, Hash};
use crate::proof::Proof;

pub struct Channel {
    state: Hash,

    // Commits
    f_eval_merkle_root: Option<Hash>,
    pub alpha0: Option<u32>,
    pub alpha1: Option<u32>,
    pub alpha2: Option<u32>,
    cp_eval_merkle_root: Option<Hash>,
    pub betas: [Option<u32>; 10],
    fri_eval_merkle_roots: [Option<Hash>; 10],
    fri_free_term: Option<u32>,

    // Decommit
    test_point: Option<u32>,
    f_x: Option<(u32, AuthPath)>,
    f_gx: Option<(u32, AuthPath)>,
    f_ggx: Option<(u32, AuthPath)>,
    cp0_x: Option<(u32, AuthPath)>,
    fri_layers: Vec<(u32, AuthPath, u32, AuthPath)>,
}

impl Channel {
    pub fn new() -> Self {
        Channel {
            state: [0; 32],

            alpha0: None,
            alpha1: None,
            alpha2: None,
            f_eval_merkle_root: None,
            cp_eval_merkle_root: None,
            betas: [None; 10],
            fri_eval_merkle_roots: [None; 10],
            fri_free_term: None,

            test_point: None,
            f_x: None,
            f_gx: None,
            f_ggx: None,
            cp0_x: None,
            fri_layers: Vec::with_capacity(10),
        }
    }

    pub fn into_proof(self) -> Proof {
        let t0 = self.betas;
        let t1 = self.fri_eval_merkle_roots;
        Proof {
            init_state: [0; 32],
            final_state: self.state,

            alpha0: self.alpha0.unwrap(),
            alpha1: self.alpha1.unwrap(),
            alpha2: self.alpha2.unwrap(),
            f_eval_merkle_root: self.f_eval_merkle_root.unwrap(),
            cp_eval_merkle_root: self.cp_eval_merkle_root.unwrap(),
            betas: [
                t0[0].unwrap(),
                t0[1].unwrap(),
                t0[2].unwrap(),
                t0[3].unwrap(),
                t0[4].unwrap(),
                t0[5].unwrap(),
                t0[6].unwrap(),
                t0[7].unwrap(),
                t0[8].unwrap(),
                t0[9].unwrap(),
            ],
            fri_eval_merkle_roots: [
                t1[0].unwrap(),
                t1[1].unwrap(),
                t1[2].unwrap(),
                t1[3].unwrap(),
                t1[4].unwrap(),
                t1[5].unwrap(),
                t1[6].unwrap(),
                t1[7].unwrap(),
                t1[8].unwrap(),
                t1[9].unwrap(),
            ],
            fri_free_term: self.fri_free_term.unwrap(),

            test_point: self.test_point.unwrap(),
            f_x: self.f_x.unwrap(),
            f_gx: self.f_gx.unwrap(),
            f_ggx: self.f_ggx.unwrap(),
            cp0_x: self.cp0_x.unwrap(),
            fri_layers: self.fri_layers.into_boxed_slice(),
        }
    }

    pub fn commit_f_eval_merkle_root(&mut self, merkle: Hash) {
        assert_eq!(self.f_eval_merkle_root, None);
        self.f_eval_merkle_root = Some(merkle);
    }

    pub fn get_alpha0(&mut self) -> Gf<3221225473> {
        assert_eq!(self.alpha0, None);
        self.alpha0 = Some(0);
        Gf::from(self.alpha0.unwrap())
    }

    pub fn get_alpha1(&mut self) -> Gf<3221225473> {
        assert_eq!(self.alpha1, None);
        self.alpha1 = Some(787618507);
        Gf::from(self.alpha1.unwrap())
    }

    pub fn get_alpha2(&mut self) -> Gf<3221225473> {
        assert_eq!(self.alpha2, None);
        self.alpha2 = Some(2154038926);
        Gf::from(self.alpha2.unwrap())
    }

    pub fn commit_cp_eval_merkle_root(&mut self, merkle: Hash) {
        assert_eq!(self.cp_eval_merkle_root, None);
        self.cp_eval_merkle_root = Some(merkle);
    }

    pub fn get_beta(&mut self, i: usize) -> Gf<3221225473> {
        assert_eq!(self.betas[i], None);
        self.betas[i] = Some(3);
        Gf::from(self.betas[i].unwrap())
    }

    pub fn commit_fri_eval_merkle_root(&mut self, i: usize, merkle: Hash) {
        assert_eq!(self.fri_eval_merkle_roots[i], None);
        self.fri_eval_merkle_roots[i] = Some(merkle);
    }

    pub fn commit_fri_free_term(&mut self, term: u32) {
        assert_eq!(self.fri_free_term, None);
        self.fri_free_term = Some(term);
    }

    pub fn get_test_point(&mut self) -> u32 {
        assert_eq!(self.test_point, None);
        self.test_point = Some(3);
        3
    }

    pub fn decommit_trace_f_x(&mut self, f_x: u32, f_x_auth_path: AuthPath) {
        assert_eq!(self.f_x, None);
        self.f_x = Some((f_x, f_x_auth_path));
    }

    pub fn decommit_trace_f_gx(&mut self, f_gx: u32, f_gx_auth_path: AuthPath) {
        assert_eq!(self.f_gx, None);
        self.f_gx = Some((f_gx, f_gx_auth_path));
    }

    pub fn decommit_trace_f_ggx(&mut self, f_ggx: u32, f_ggx_auth_path: AuthPath) {
        assert_eq!(self.f_ggx, None);
        self.f_ggx = Some((f_ggx, f_ggx_auth_path));
    }

    pub fn decommit_trace_cp0_x(&mut self, cp0_x: u32, cp0_x_auth_path: AuthPath) {
        assert_eq!(self.cp0_x, None);
        self.cp0_x = Some((cp0_x, cp0_x_auth_path));
    }

    pub fn decommit_fri_layer(
        &mut self,
        cp_x: u32,
        cp_x_auth_path: AuthPath,
        cp_nx: u32,
        cp_nx_auth_path: AuthPath,
    ) {
        self.fri_layers
            .push((cp_x, cp_x_auth_path, cp_nx, cp_nx_auth_path));
    }
}
