mod channel;
mod field;
mod merkle;
mod polynomial;
mod proof;

use merkle::Merkle;
use num_traits::Pow;
use polynomial::{lagrange, Polynomial};

type F = field::Gf<3221225473>;
use num_traits::{One, Zero};

fn main() {
    let mut channel = channel::Channel::new();

    ///////////////////
    // Part 1

    // The trace sequence
    let mut a = [F::zero(); 1023];
    a[0] = F::one();
    a[1] = F::from(3141592); // The secret
    for i in 2..1023 {
        let t0 = a[i - 2].pow(2);
        let t1 = a[i - 1].pow(2);
        a[i] = t0 + t1;
    }

    // Assert the trace is correct
    assert_eq!(a[1022].residue(), 2338775057);

    // Generate a primitive root of F_3221225473
    let primitive_root = F::generator();

    // Create generators for cyclic groups of sizes 1024 and 8192
    let generator_g = primitive_root.pow(3145728);
    let generator_h = primitive_root.pow(393216);

    // Assert generators are of the correct order
    assert_eq!(generator_g.order(), 1024);
    assert_eq!(generator_h.order(), 8192);

    // Generate respective cyclic groups
    let g: Vec<F> = (0..1024).map(|n| generator_g.pow(n)).collect();
    let h: Vec<F> = (0..8192).map(|n| generator_h.pow(n)).collect();

    // Generate lagrange polynomial going through points (g[i], a[i]) for i <= 1022
    let points: Vec<(F, F)> = std::iter::zip(&g, &a).map(|(&x, &y)| (x, y)).collect();
    let f_poly = lagrange::<F>(&points);

    // Assert that the polynomial has the correct solutions
    for (x, y) in std::iter::zip(&g, &a) {
        assert_eq!(f_poly.solve(*x), *y);
    }

    // Solve polynomial over h, shifted by the primitive root
    let f_domain: Vec<F> = h.iter().map(|n| primitive_root * *n).collect();
    let f_eval: Vec<F> = f_domain.iter().map(|&n| f_poly.solve(n)).collect();

    // Assert a few elements of eval are correct
    assert_eq!(f_eval[0].residue(), 576067152);
    assert_eq!(f_eval[1].residue(), 3100214617);
    assert_eq!(f_eval[2].residue(), 2091264768);
    assert_eq!(f_eval[8189].residue(), 800520420);
    assert_eq!(f_eval[8190].residue(), 1199720174);
    assert_eq!(f_eval[8191].residue(), 1076821037);

    // Generate merkle tree from f_eval
    let f_eval_merkle = Merkle::new(1024, f_eval.iter().map(|f| f.residue()));
    let f_eval_merkle_root = f_eval_merkle[0];

    // Commit f_eval merkle root
    channel.commit_f_eval_merkle_root(f_eval_merkle_root);

    ///////////////////
    // Part 2
    use polynomial::x;

    // Constraint 0:
    // f(x) - a[0]
    // -----------
    //  x - g[0]
    let numerator = &f_poly - &x(a[0], 0);
    let denominator = Polynomial::from([F::one(), -g[0]]);
    let (c0, c0r) = Polynomial::<F>::div(numerator, denominator);

    // Constraint 1:
    // f(x) - a[1022]
    // --------------
    //  x - g[1022]
    let numerator = &f_poly - &x(a[1022], 0);
    let denominator = Polynomial::from([F::one(), -g[1022]]);
    let (c1, c1r) = Polynomial::<F>::div(numerator, denominator);

    // Constraint 2:
    //              f(g^2 x) - f(g x)^2 - f(x)^2
    // ------------------------------------------------------
    // (x^1024 - 1)/(x - g[1021])/(x - g[1022])/(x - g[1023])
    let t0 = f_poly.clone().apply_const(g[2]);
    let t1 = f_poly.clone().apply_const(g[1]);
    let t1 = &t1 * &t1;
    let t2 = &f_poly * &f_poly;
    let numerator = t0 - t1 - t2;

    let denominator = x(F::one(), 1024) - x(F::one(), 0);
    let tp0 = Polynomial::from([F::one(), -g[1021]]);
    let tp1 = Polynomial::from([F::one(), -g[1022]]);
    let tp2 = Polynomial::from([F::one(), -g[1023]]);
    let (denominator, t2r) = Polynomial::<F>::div(denominator, &tp2 * &tp0 * tp1);
    let (c2, c2r) = Polynomial::<F>::div(numerator, denominator);

    // Assert constraints have no remainders
    assert_eq!(c0r.degree(), None);
    assert_eq!(c1r.degree(), None);
    assert_eq!(t2r.degree(), None);
    assert_eq!(c2r.degree(), None);

    // Assert constraints resolve correctly
    assert_eq!(c2.degree(), Some(1023));
    assert_eq!(c0.solve(F::from(2718)).residue(), 2509888982);
    assert_eq!(c1.solve(F::from(5772)).residue(), 232961446);
    assert_eq!(c2.solve(F::from(31415)).residue(), 2090051528);

    // Generate composition polynomial (normally, these would be random)
    let a0 = channel.get_alpha0(); //F::from(0);
    let a1 = channel.get_alpha1(); //F::from(787618507);
    let a2 = channel.get_alpha2(); //F::from(-1067186547);
    let cp_poly = c0 * x(a0, 0) + c1 * x(a1, 0) + c2 * x(a2, 0);

    // Assert composition polynomial resolves correctly
    assert_eq!(cp_poly.degree(), Some(1023));
    //assert_eq!(cp.solve(F::from(2439804)).residue(), 838767343);

    // Evaluate cp over f_domain
    let cp_domain = f_domain;
    let cp_eval: Vec<F> = cp_domain.iter().map(|&n| cp_poly.solve(n)).collect();

    // Assert a few elements of cp_eval are correct
    assert_eq!(cp_eval[0].residue(), 551740506);
    assert_eq!(cp_eval[1].residue(), 716458408);
    assert_eq!(cp_eval[2].residue(), 2091260387);
    assert_eq!(cp_eval[8189].residue(), 412406999);
    assert_eq!(cp_eval[8190].residue(), 782538909);
    assert_eq!(cp_eval[8191].residue(), 811632985);

    // Generate merkle tree from cp_eval
    let cp_eval_merkle = Merkle::new(8192, cp_eval.iter().map(|f| f.residue()));
    let cp_eval_merkle_root = cp_eval_merkle[0];

    // Commit cp_eval merkle root
    channel.commit_cp_eval_merkle_root(cp_eval_merkle_root);

    ///////////////////
    // Part 3

    let mut cp_polys: Vec<Polynomial<F>> = vec![cp_poly];
    let mut cp_domains: Vec<Vec<F>> = vec![cp_domain];
    let mut cp_evals: Vec<Vec<F>> = vec![cp_eval];
    let mut cp_eval_merkles: Vec<Merkle> = vec![cp_eval_merkle];

    // Perform FRI operation
    for i in 0..10 {
        // Get new fri poly
        let beta = channel.get_beta(i);
        let fri_poly = polynomial::fri::<F>(cp_polys.last().unwrap(), beta);

        // Get new fri domain
        let mut fri_domain = cp_domains.last().unwrap().clone();
        fri_domain.truncate(fri_domain.len() / 2);
        for e in &mut fri_domain {
            *e = e.pow(2);
        }

        // Solve over new domain
        let fri_eval: Vec<_> = fri_domain.iter().map(|&n| fri_poly.solve(n)).collect();

        // Generate merkle tree from fri_eval
        let fri_eval_merkle = Merkle::new(fri_domain.len(), fri_eval.iter().map(|f| f.residue()));
        let fri_eval_merkle_root = fri_eval_merkle[0];

        // Push
        cp_polys.push(fri_poly);
        cp_domains.push(fri_domain);
        cp_evals.push(fri_eval);
        cp_eval_merkles.push(fri_eval_merkle);

        // Commit fri_eval merkle root
        channel.commit_fri_eval_merkle_root(i, fri_eval_merkle_root);
    }

    // Assert the degree of the FRI polynomials
    assert_eq!(cp_polys[0].degree(), Some(1023));
    assert_eq!(cp_polys[1].degree(), Some(511));
    assert_eq!(cp_polys[2].degree(), Some(255));
    assert_eq!(cp_polys[3].degree(), Some(127));
    assert_eq!(cp_polys[4].degree(), Some(63));
    assert_eq!(cp_polys[5].degree(), Some(31));
    assert_eq!(cp_polys[6].degree(), Some(15));
    assert_eq!(cp_polys[7].degree(), Some(7));
    assert_eq!(cp_polys[8].degree(), Some(3));
    assert_eq!(cp_polys[9].degree(), Some(1));
    assert_eq!(cp_polys[10].degree(), Some(0));

    // Assert fri_eval sizes
    assert_eq!(cp_evals[0].len(), 8192);
    assert_eq!(cp_evals[1].len(), 4096);
    assert_eq!(cp_evals[2].len(), 2048);
    assert_eq!(cp_evals[3].len(), 1024);
    assert_eq!(cp_evals[4].len(), 512);
    assert_eq!(cp_evals[5].len(), 256);
    assert_eq!(cp_evals[6].len(), 128);
    assert_eq!(cp_evals[7].len(), 64);
    assert_eq!(cp_evals[8].len(), 32);
    assert_eq!(cp_evals[9].len(), 16);
    assert_eq!(cp_evals[10].len(), 8);

    // Commit free term of the final polynomial
    channel.commit_fri_free_term(cp_polys[10][0].residue());

    ///////////////////
    // Part 4

    // Get test point
    let x = channel.get_test_point() as usize;

    // Decommit on trace
    let fx = f_eval[x].residue();
    let fx_auth_path = f_eval_merkle.trace(x);
    channel.commit_fx(fx, fx_auth_path);

    let fgx = f_eval[x + 8].residue();
    let fgx_auth_path = f_eval_merkle.trace(x + 8);
    channel.commit_fgx(fgx, fgx_auth_path);

    let fggx = f_eval[x + 16].residue();
    let fggx_auth_path = f_eval_merkle.trace(x + 16);
    channel.commit_fggx(fggx, fggx_auth_path);

    let cp0x = cp_evals[0][x].residue();
    let cp0x_auth_path = cp_eval_merkles[0].trace(x);
    channel.commit_cp0x(cp0x, cp0x_auth_path);

    // VERIFY:
    let proof = channel.into_proof();
    assert_eq!(proof.verify(), true);
}
