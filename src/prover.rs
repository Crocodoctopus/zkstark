use crate::channel::Channel;
use crate::merkle::Merkle;
use crate::polynomial::{x, fri, lagrange, Polynomial};
use crate::proof::Proof;
use crate::F;
use num_traits::Pow;
use num_traits::{One, Zero};

pub fn generate_proof(mut channel: Channel) -> Proof {
    // I'll do my best to explain things, at least how I understand them thus far.
    //
    // The proof is divided into 4 parts:
    // 1) Generate a trace sequence, 
    // 2) Create constraint polynomial to prove the trace evaluates correctly
    // 3) Perform FRI operator on constraint polynomial to prove is "almost" low degree
    // 4) Probe polynomials with test values to ensure all the math works out
    //
    // The way I understand the whole process is, at a high level, we're just creating a
    // polynomial that will be low degree (and evaluate correctly at all points) if and only
    // if it was constructed faithfully from valid trace data. The rest of the proof is just
    // showing that it is both a low degree poly, and evaluates correctly (at least at a few
    // test points).

    ///////////////////
    // Part 1:
    //   In this part, we generate a trace sequence, a lagrange polynomial for that sequence,
    // and then evaluate said polynomial over an extended domain. I won't go too far into
    // this part, I found it relatively straight forward from the video guide this project
    // is based on.

    // Generate the trace sequence of 1023 elements.
    let mut a = [F::zero(); 1023];
    a[0] = F::one();
    a[1] = F::from(3141592); // The secret
    for i in 2..1023 {
        let t0 = a[i - 2].pow(2);
        let t1 = a[i - 1].pow(2);
        a[i] = t0 + t1;
    }

    // Assert the trace is correct at the last value
    assert_eq!(a[1022].residue(), 2338775057);

    // Generate a primitive root of F_3221225473 (this ends up being 5 in the python codebase)
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

    // Assert that the polynomial has the correct solutions at each domain input
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
    let f_eval_merkle = Merkle::new(8192, f_eval.iter().map(|f| f.residue()));
    let f_eval_merkle_root = f_eval_merkle[0];

    // Commit f_eval merkle root
    channel.commit_f_eval_merkle_root(f_eval_merkle_root);

    ///////////////////
    // Part 2:
    //  This part is a bit mathy. Basically, we want to create polynomials that both
    // mathematically tie each element of the group together, and are low degree if
    // and only if the math that created the sequence was correct.

    // Constraint 0:
    // f(x) - a[0]
    // -----------
    //  x - g[0]
    // So, f(x) is a degree 1022, and at g[0] is evaluates to a[0], by definition
    // (we used lagrange precisely for this property). Therefor, f(x) - a[0] = 0 
    // at g[0]. Therefor, g[0] is a root, and the first constraint divides evenly to
    // produce a polynomial of degree 1021.
    let numerator = &f_poly - &x(a[0], 0);
    let denominator = Polynomial::from([F::one(), -g[0]]);
    let (c0, c0r) = Polynomial::<F>::div(numerator, denominator);

    // Constraint 1:
    // f(x) - a[1022]
    // --------------
    //  x - g[1022]
    // Following the logic from constraint 0, f(x) at g[1022] evaluates to a[1022].
    // Once again, we can use this fact to produce a degree 1021 polynomial.
    let numerator = &f_poly - &x(a[1022], 0);
    let denominator = Polynomial::from([F::one(), -g[1022]]);
    let (c1, c1r) = Polynomial::<F>::div(numerator, denominator);

    // Constraint 2:
    //              f(g^2 x) - f(g x)^2 - f(x)^2
    // ------------------------------------------------------
    // (x^1024 - 1)/(x - g[1021])/(x - g[1022])/(x - g[1023])
    // We have constraints for our publically known trace values (0, 1) and (1022, 2338775057),
    // but a constraint for the other 1021 values is a bit more involved. The idea is the same,
    // but with a twist. Notice g[n] * g == g[n + 1]. In particular, f(g[n] * g) evaluates
    // the poly at g[n + 1]. We can multiply the input by g to "slide it forward". Our original
    // trace equation was a[n + 2] = a[n + 1]^2 + a[n]^2. Because f(g[n]) == a[n], we can
    // transform the equation to f(g[n + 2]) = f(g[n + 1])^2 + f(g[n])^2, and using our sliding
    // rule we get f(g*g*g[n]) = f(g*g[n]) + f(g[n])^2. Move some terms to the side and we get
    // f(g*g*x) - f(g*x) - f(x)^2 = 0 where x <- g[n] for each n. Just like the other 2
    // contraints, we know that this equation has roots = { g[n] | 0 < n < 1022 }. Thus, we can
    // divide by (x - g[1])(x - g[2])(x - g[3])...(x - g[1021]) to produce a polynomial of 
    // degree 1023.
    //
    // Note, the denom is not (x - g[1])(x - g[2])... like I said. This multiplcation is very
    // expensive. Instead, start with precomputed (x^1024 - 1) and divide out the unwanted
    // terms instead.
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
    assert_eq!(c0.degree(), Some(1021));
    assert_eq!(c1.degree(), Some(1021));
    assert_eq!(c2.degree(), Some(1023));
    assert_eq!(c0.solve(F::from(2718)).residue(), 2509888982);
    assert_eq!(c1.solve(F::from(5772)).residue(), 232961446);
    assert_eq!(c2.solve(F::from(31415)).residue(), 2090051528);

    // Generate composition polynomial, with a degree of at most 1023 (this is true because poly 
    // addition can't produce higher power terms)
    let a0 = channel.get_alpha0();
    let a1 = channel.get_alpha1();
    let a2 = channel.get_alpha2();
    let cp_poly = c0 * x(a0, 0) + c1 * x(a1, 0) + c2 * x(a2, 0);

    // Assert composition polynomial resolves correctly
    assert_eq!(cp_poly.degree(), Some(1023));

    // Evaluate cp over f_domain
    let cp_domain = f_domain;
    let cp_eval: Vec<F> = cp_domain.iter().map(|&n| cp_poly.solve(n)).collect();

    // Generate merkle tree from cp_eval
    let cp_eval_merkle = Merkle::new(8192, cp_eval.iter().map(|f| f.residue()));
    let cp_eval_merkle_root = cp_eval_merkle[0];

    // Commit cp_eval merkle root
    channel.commit_cp_eval_merkle_root(cp_eval_merkle_root);

    ///////////////////
    // Part 3:
    //   Instead of proving that composite polynomial from earlier is degree 1023, we're
    // going to prove that it is "close". The video series doesn't go too into depth on
    // how this work, this is mostly just going through the motions. Basically, we perform
    // an iterative process of our composite polynomial that cuts the degree in half.
    // We know the degree is probably 1023, so after 10 iterations, we should have a polynomial
    // that is degree 0 (just a constant). If it isn't, then our constraint polynomial roots
    // didn't cancel out well, and the computation wasn't faithful. At least how I understand it.

    let mut cp_polys: Vec<Polynomial<F>> = vec![cp_poly];
    let mut cp_domains: Vec<Vec<F>> = vec![cp_domain];
    let mut cp_evals: Vec<Vec<F>> = vec![cp_eval];
    let mut cp_eval_merkles: Vec<Merkle> = vec![cp_eval_merkle];

    // Perform FRI operation
    for i in 0..10 {
        // Get new fri poly
        let beta = channel.get_beta(i);
        let fri_poly = fri::<F>(cp_polys.last().unwrap(), beta);

        // Get new fri domain
        let mut fri_domain = cp_domains.last().unwrap().clone();
        fri_domain.truncate(fri_domain.len() / 2);
        for e in &mut fri_domain {
            *e = e.pow(2);
        }

        // Solve over new domain
        let fri_eval: Vec<_> = fri_domain.iter().map(|&n| fri_poly.solve(n)).collect();

        // Generate merkle tree from fri_eval
        let fri_eval_merkle = Merkle::new(fri_eval.len(), fri_eval.iter().map(|f| f.residue()));
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
    //   This is the "decommit" phase. The above 3 steps generate all the data we need, now
    // we just receive a test point and evaluate that point through each stage, showing that
    // the math between each stage follows. Again, this was mostly just going through the motions.

    // Get test point
    let x = channel.get_test_point() as usize % (8192 - 16);

    // Decommit on trace
    let f_x = f_eval[x].residue();
    let f_x_auth_path = f_eval_merkle.trace(x);
    let f_gx = f_eval[x + 8].residue();
    let f_gx_auth_path = f_eval_merkle.trace(x + 8);
    let f_ggx = f_eval[x + 16].residue();
    let f_ggx_auth_path = f_eval_merkle.trace(x + 16);
    let cp0_x = cp_evals[0][x].residue();
    let cp0_x_auth_path = cp_eval_merkles[0].trace(x);
    channel.decommit_trace_f_x(f_x, f_x_auth_path);
    channel.decommit_trace_f_gx(f_gx, f_gx_auth_path);
    channel.decommit_trace_f_ggx(f_ggx, f_ggx_auth_path);
    channel.decommit_trace_cp0_x(cp0_x, cp0_x_auth_path);

    // Decommit on FRI
    for i in 0..10 {
        let len = cp_domains[i].len();
        let x = x % len;
        let nx = (x + len / 2) % len;
        let cp_x = cp_evals[i][x].residue();
        let cp_x_auth_path = cp_eval_merkles[i].trace(x);
        let cp_nx = cp_evals[i][nx].residue();
        let cp_nx_auth_path = cp_eval_merkles[i].trace(nx);
        channel.decommit_fri_layer(cp_x, cp_x_auth_path, cp_nx, cp_nx_auth_path);
    }

    // Done
    return channel.into_proof();
}
