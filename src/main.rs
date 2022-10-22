mod channel;
mod field;
mod polynomial;

use num_traits::Pow;
use polynomial::{lagrange, Polynomial};

type F = field::Gf<3221225473>;

fn main() {
    let mut channel = channel::Channel::new();

    ///////////////////
    // Part 1:

    // The trace sequence
    let mut a = [F::from(0); 1023];
    a[0] = F::from(1);
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

    // Create generators cyclic groups of sizes 1024 and 8192
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
    let f = lagrange::<F>(&points);

    // Assert that the polynomial has the correct solutions
    for (x, y) in std::iter::zip(&g, &a) {
        assert_eq!(f.solve(*x), *y);
    }

    // Solve polynomial over h, shifted by the primitive root
    let eval_domain: Vec<F> = h.iter().map(|n| primitive_root * *n).collect();
    let f_eval: Vec<F> = eval_domain.iter().map(|&n| f.solve(n)).collect();

    // Assert a few elements of eval are correct
    assert_eq!(f_eval[0].residue(), 576067152);
    assert_eq!(f_eval[1].residue(), 3100214617);
    assert_eq!(f_eval[2].residue(), 2091264768);
    assert_eq!(f_eval[8189].residue(), 800520420);
    assert_eq!(f_eval[8190].residue(), 1199720174);
    assert_eq!(f_eval[8191].residue(), 1076821037);

    // Generate merkle tree from f_eval
    let f_eval_merkle_root = [1; 32]; // [TODO]

    // Commit f_eval merkle root
    channel.commit_f_eval_merkle_root(f_eval_merkle_root);

    ///////////////////
    // Part 2
    use polynomial::x;

    // Constraint 0:
    // f(x) - a[0]
    // -----------
    //  x - g[0]
    let numerator = &f - &x(a[0], 0);
    let denominator = [F::from(1), -g[0]].into();
    let (c0, c0r) = Polynomial::<F>::div(numerator, denominator);

    // Constraint 1:
    // f(x) - a[1022]
    // --------------
    //  x - g[1022]
    let numerator = &f - &x(a[1022], 0);
    let denominator = [F::from(1), -g[1022]].into();
    let (c1, c1r) = Polynomial::<F>::div(numerator, denominator);

    // Constraint 2:
    //              f(g^2 x) - f(g x)^2 - f(x)^2
    // ------------------------------------------------------
    // (x^1024 - 1)/(x - g[1021])/(x - g[1022])/(x - g[1023])
    let t0 = f.clone().apply_const(g[2]);
    let t1 = f.clone().apply_const(g[1]);
    let t1 = &t1 * &t1;
    let t2 = &f * &f;
    let numerator = t0 - t1 - t2;

    let denominator = x(F::from(1), 1024) - x(F::from(1), 0);
    let (denominator, r0) = Polynomial::<F>::div(denominator, [F::from(1), -g[1021]].into());
    let (denominator, r1) = Polynomial::<F>::div(denominator, [F::from(1), -g[1022]].into());
    let (denominator, r2) = Polynomial::<F>::div(denominator, [F::from(1), -g[1023]].into());
    let (c2, c2r) = Polynomial::<F>::div(numerator, denominator);

    // Assert constraints have no remainders
    assert_eq!(c0r.degree(), None);
    assert_eq!(c1r.degree(), None);
    assert_eq!(r0.degree(), None);
    assert_eq!(r1.degree(), None);
    assert_eq!(r2.degree(), None);
    assert_eq!(c2r.degree(), None);

    // Assert constraints resolve correctly
    assert_eq!(c2.degree(), Some(1023));
    assert_eq!(c0.solve(F::from(2718)).residue(), 2509888982);
    assert_eq!(c1.solve(F::from(5772)).residue(), 232961446);
    assert_eq!(c2.solve(F::from(31415)).residue(), 2090051528);

    // Generate composition polynomial (normally, these would be random)
    let a0 = channel.get_random_element(); //F::from(0);
    let a1 = channel.get_random_element(); //F::from(787618507);
    let a2 = channel.get_random_element(); //F::from(-1067186547);
    let cp = c0 * x(a0, 0) + c1 * x(a1, 0) + c2 * x(a2, 0);

    // Assert composition polynomial resolves correctly
    assert_eq!(cp.degree(), Some(1023));
    //assert_eq!(cp.solve(F::from(2439804)).residue(), 838767343);

    // Evaluate cp over eval_domain
    let cp_eval: Vec<F> = eval_domain.iter().map(|&n| cp.solve(n)).collect();

    // Assert a few elements of cp_eval are correct
    /*assert_eq!(cp_eval[0].residue(), 551740506);
    assert_eq!(cp_eval[1].residue(), 716458408);
    assert_eq!(cp_eval[2].residue(), 2091260387);
    assert_eq!(cp_eval[8189].residue(), 412406999);
    assert_eq!(cp_eval[8190].residue(), 782538909);
    assert_eq!(cp_eval[8191].residue(), 811632985);*/

    // Generate merkle tree from cp_eval
    let cp_eval_merkle_root = [2u8; 32];

    // Commit cp_eval merkle root
    channel.commit_cp_eval_merkle_root(cp_eval_merkle_root);

    ///////////////////
    // Part 3

    let mut cp_evals: Vec<Vec<F>> = vec![cp_eval];
    let mut cp_polys: Vec<Polynomial<F>> = vec![cp];

    // Perform FRI operation
    let mut fri_domain = eval_domain;
    for i in 0..10 {
        // Get new fri domain
        fri_domain.truncate(fri_domain.len() / 2);
        for e in &mut fri_domain {
            *e = e.pow(2);
        }

        // Get new fri poly
        let beta = channel.get_random_element();
        let fri_poly = polynomial::fri::<F>(cp_polys.last().unwrap(), beta);

        // Solve over new domain
        let fri_eval = fri_domain.iter().map(|&n| fri_poly.solve(n)).collect();

        // Push
        cp_polys.push(fri_poly);
        cp_evals.push(fri_eval);

        // Generate merkle tree from fri_eval
        let fri_eval_merkle_root = [i + 3; 32];

        // Commit fri_eval merkle root
        channel.commit_fri_eval_merkle_root(i as usize, fri_eval_merkle_root);
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
}
