use crate::field;

pub struct Channel {
    state: [u8; 32],

    // Commits
    f_eval_merkle_root: Option<[u8; 32]>,
    cp_eval_merkle_root: Option<[u8; 32]>,
    fri_eval_merkle_roots: [Option<[u8; 32]>; 10],
    fri_free_term: Option<u32>,
}

impl Channel {
    pub fn new() -> Self {
        Channel {
            state: [0; 32],

            f_eval_merkle_root: None,
            cp_eval_merkle_root: None,
            fri_eval_merkle_roots: [None; 10],
            fri_free_term: None,
        }
    }

    pub fn get_random_element(&mut self) -> field::Gf<3221225473> {
        // Such random wow
        <_>::from(3)
    }

    pub fn commit_f_eval_merkle_root(&mut self, merkle: [u8; 32]) {
        assert_eq!(self.f_eval_merkle_root, None);
        self.f_eval_merkle_root = Some(merkle);
    }

    pub fn commit_cp_eval_merkle_root(&mut self, merkle: [u8; 32]) {
        assert_eq!(self.cp_eval_merkle_root, None);
        self.cp_eval_merkle_root = Some(merkle);
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