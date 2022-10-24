use crate::field::Gf;

pub struct Channel {
    state: [u8; 32],

    // Commits
    f_eval_merkle_root: Option<[u8; 32]>,
    alpha0: Option<u32>,
    alpha1: Option<u32>,
    alpha2: Option<u32>,
    cp_eval_merkle_root: Option<[u8; 32]>,
    betas: [Option<u32>; 10],
    fri_eval_merkle_roots: [Option<[u8; 32]>; 10],
    fri_free_term: Option<u32>,

    // Decommit
    test_point: Option<usize>,
    fx: Option<(u32, Box<[[u8; 32]]>)>,
    fgx: Option<(u32, Box<[[u8; 32]]>)>,
    fggx: Option<(u32, Box<[[u8; 32]]>)>,
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
            fx: None,
            fgx: None,
            fggx: None,
        }
    }

    pub fn print(&self) {
        println!("Channel output:");
        println!("  // Commit");
        println!(
            "  f_eval_merkle_root: {:02X?}",
            self.f_eval_merkle_root.unwrap()
        );
        println!("  alpha0: {:?}", self.alpha0.unwrap());
        println!("  alpha1: {:?}", self.alpha1.unwrap());
        println!("  alpha2: {:?}", self.alpha2.unwrap());
        println!(
            "  cp_eval_merkle_root: {:02X?}",
            self.cp_eval_merkle_root.unwrap()
        );
        println!("  FRI:");
        for i in 0..10 {
            println!(
                "    [cp{}] beta: {:?}; fri_eval_merkle_roots: {:02X?}",
                i + 1,
                self.betas[i].unwrap(),
                self.fri_eval_merkle_roots[i].unwrap()
            );
        }
        println!("  fri_free_term: {:?}", self.fri_free_term.unwrap());
        println!("  // Decommit");
        println!("  fx: {:?}", self.fx.as_ref().unwrap());
    }

    pub fn commit_f_eval_merkle_root(&mut self, merkle: [u8; 32]) {
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

    pub fn commit_cp_eval_merkle_root(&mut self, merkle: [u8; 32]) {
        assert_eq!(self.cp_eval_merkle_root, None);
        self.cp_eval_merkle_root = Some(merkle);
    }

    pub fn get_beta(&mut self, i: usize) -> Gf<3221225473> {
        assert_eq!(self.betas[i], None);
        self.betas[i] = Some(3);
        Gf::from(self.betas[i].unwrap())
    }

    pub fn commit_fri_eval_merkle_root(&mut self, i: usize, merkle: [u8; 32]) {
        assert_eq!(self.fri_eval_merkle_roots[i], None);
        self.fri_eval_merkle_roots[i] = Some(merkle);
    }

    pub fn commit_fri_free_term(&mut self, term: u32) {
        assert_eq!(self.fri_free_term, None);
        self.fri_free_term = Some(term);
    }

    pub fn get_test_point(&mut self) -> usize {
        assert_eq!(self.test_point, None);
        self.test_point = Some(5);
        5
    }

    pub fn commit_fx(&mut self, fx: u32, fx_auth_path: Box<[[u8; 32]]>) {
        assert_eq!(self.fx, None);
        self.fx = Some((fx, fx_auth_path));
    }

    pub fn commit_fgx(&mut self, fgx: u32, fgx_auth_path: Box<[[u8; 32]]>) {
        assert_eq!(self.fgx, None);
        self.fgx = Some((fgx, fgx_auth_path));
    }

    pub fn commit_fggx(&mut self, fggx: u32, fggx_auth_path: Box<[[u8; 32]]>) {
        assert_eq!(self.fggx, None);
        self.fggx = Some((fggx, fggx_auth_path));
    }
}
