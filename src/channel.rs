use crate::field::Gf;

pub struct Channel {
    state: [u8; 32],

    // Commits
    f_eval_merkle_root: Option<[u8; 32]>,
    alpha0: Option<Gf<3221225473>>,
    alpha1: Option<Gf<3221225473>>,
    alpha2: Option<Gf<3221225473>>,
    cp_eval_merkle_root: Option<[u8; 32]>,
    betas: [Option<Gf<3221225473>>; 10],
    fri_eval_merkle_roots: [Option<[u8; 32]>; 10],
    fri_free_term: Option<u32>,
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
        }
    }

    pub fn commit_f_eval_merkle_root(&mut self, merkle: [u8; 32]) {
        assert_eq!(self.f_eval_merkle_root, None);
        self.f_eval_merkle_root = Some(merkle);
    }

    pub fn get_alpha0(&mut self) -> Gf<3221225473> {
        assert_eq!(self.alpha0, None);
        self.alpha0 = Some(Gf::from(0));
        self.alpha0.unwrap()
    }

    pub fn get_alpha1(&mut self) -> Gf<3221225473> {
        assert_eq!(self.alpha1, None);
        self.alpha1 = Some(Gf::from(787618507));
        self.alpha1.unwrap()
    }

    pub fn get_alpha2(&mut self) -> Gf<3221225473> {
        assert_eq!(self.alpha2, None);
        self.alpha2 = Some(Gf::from(-1067186547));
        self.alpha2.unwrap()
    }

    pub fn commit_cp_eval_merkle_root(&mut self, merkle: [u8; 32]) {
        assert_eq!(self.cp_eval_merkle_root, None);
        self.cp_eval_merkle_root = Some(merkle);
    }

    pub fn get_beta(&mut self, i: usize) -> Gf<3221225473> {
        assert_eq!(self.betas[i], None);
        self.betas[i] = Some(Gf::from(3));
        self.betas[i].unwrap()
    }

    pub fn commit_fri_eval_merkle_root(&mut self, i: usize, merkle: [u8; 32]) {
        assert_eq!(self.fri_eval_merkle_roots[i], None);
        self.fri_eval_merkle_roots[i] = Some(merkle);
    }

    pub fn commit_fri_free_term(&mut self, term: u32) {
        assert_eq!(self.fri_free_term, None);
        self.fri_free_term = Some(term);
    }
}