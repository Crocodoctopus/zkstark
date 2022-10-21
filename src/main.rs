mod field;
mod polynomial;

use num_traits::Pow;
use polynomial::{lagrange, Polynomial};

type F = field::Gf<3221225473>;

fn main() {
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
    let poly = lagrange::<F>(&points);

    // Assert that the polynomial has the correct solutions
    for (x, y) in std::iter::zip(&g, &a) {
        assert_eq!(poly.solve(*x), *y);
    }

    // Solve polynomial over h, shifted by the primitive root
    let eval_domain: Vec<F> = h.iter().map(|n| primitive_root * *n).collect();
    let eval: Vec<F> = eval_domain.iter().map(|&n| poly.solve(n)).collect();

    // Assert a few elements of eval are correct
    assert_eq!(eval[0].residue(), 576067152);
    assert_eq!(eval[1].residue(), 3100214617);
    assert_eq!(eval[2].residue(), 2091264768);
    assert_eq!(eval[8189].residue(), 800520420);
    assert_eq!(eval[8190].residue(), 1199720174);
    assert_eq!(eval[8191].residue(), 1076821037);

    ///////////////////
    // Part 2
    use polynomial::x;

    // Constraint 0:
    // f(x) - a[0]
    // -----------
    //  x - g[0]
    let numerator = &poly - &x(a[0], 0);
    let denominator = [F::from(1), -g[0]].into();
    let (c0, c0r) = Polynomial::<F>::div(numerator, denominator);

    // Constraint 1:
    // f(x) - a[1022]
    // --------------
    //  x - g[1022]
    let numerator = &poly - &x(a[1022], 0);
    let denominator = [F::from(1), -g[1022]].into();
    let (c1, c1r) = Polynomial::<F>::div(numerator, denominator);

    // Constraint 2:
    //              f(g^2 x) - f(g x)^2 - f(x)^2
    // ------------------------------------------------------
    // (x^1024 - 1)/(x - g[1021])/(x - g[1022])/(x - g[1023])
    let t0 = poly.clone().apply_const(g[2]);
    let t1 = poly.clone().apply_const(g[1]);
    let t1 = &t1 * &t1;
    let t2 = &poly * &poly;
    let numerator = t0 - t1 - t2;

    let denominator = x(F::from(1), 1024) - x(F::from(1), 0);
    let (denominator, r0) = Polynomial::<F>::div(denominator, [F::from(1), -g[1021]].into());
    let (denominator, r1) = Polynomial::<F>::div(denominator, [F::from(1), -g[1022]].into());
    let (denominator, r2) = Polynomial::<F>::div(denominator, [F::from(1), -g[1023]].into());
    let (c2, c2r) = Polynomial::<F>::div(numerator, denominator);

    // Assert constraints have no remainders
    assert!(c0r.degree().is_none());
    assert!(c1r.degree().is_none());
    assert!(r0.degree().is_none());
    assert!(r1.degree().is_none());
    assert!(r2.degree().is_none());
    assert!(c2r.degree().is_none());

    // Assert constraints resolve correctly
    assert_eq!(c2.degree(), Some(1023));
    assert_eq!(c0.solve(F::from(2718)).residue(), 2509888982);
    assert_eq!(c1.solve(F::from(5772)).residue(), 232961446);
    assert_eq!(c2.solve(F::from(31415)).residue(), 2090051528);

    // Generate composition polynomial (normally, these would be random)
    let a0 = F::from(0);
    let a1 = F::from(787618507);
    let a2 = F::from(-1067186547);
    let cp = c0 * x(a0, 0) + c1 * x(a1, 0) + c2 * x(a2, 0);

    // Assert composition polynomial resolves correctly
    assert_eq!(cp.degree(), Some(1023));
    assert_eq!(cp.solve(F::from(2439804)).residue(), 838767343);

    // Evaluate cp over eval_domain
    let cp_eval: Vec<F> = eval_domain.iter().map(|&n| cp.solve(n)).collect();

    // Assert a few elements of cp_eval are correct
    assert_eq!(cp_eval[0].residue(), 551740506);
    assert_eq!(cp_eval[1].residue(), 716458408);
    assert_eq!(cp_eval[2].residue(), 2091260387);
    assert_eq!(cp_eval[8189].residue(), 412406999);
    assert_eq!(cp_eval[8190].residue(), 782538909);
    assert_eq!(cp_eval[8191].residue(), 811632985);

    ///////////////////
    // Part 3

    let mut fri_evals: Vec<Vec<F>> = vec![cp_eval];
    let mut fri_domain = eval_domain;
    let mut fri_poly = cp;
    while fri_poly.degree().unwrap() > 0 {
        // Get new fri domain
        fri_domain.truncate(fri_domain.len() / 2);
        for e in &mut fri_domain {
            *e = e.pow(2);
        }

        // Get new fri poly
        fri_poly = polynomial::fri::<F>(&fri_poly, F::from(3));

        // Solve over new domain
        let fri_eval = fri_domain.iter().map(|&n| fri_poly.solve(n)).collect();
        fri_evals.push(fri_eval);
    }

    // Assert that it took 10 (plus) iterations to reach degree 0
    assert_eq!(fri_evals.len(), 11);

    // Assert fri_eval sizes
    assert_eq!(fri_evals.get(0).map(|v| v.len()), Some(8192));
    assert_eq!(fri_evals.get(1).map(|v| v.len()), Some(4096));
    assert_eq!(fri_evals.get(2).map(|v| v.len()), Some(2048));
    assert_eq!(fri_evals.get(3).map(|v| v.len()), Some(1024));
    assert_eq!(fri_evals.get(4).map(|v| v.len()), Some(512));
    assert_eq!(fri_evals.get(5).map(|v| v.len()), Some(256));
    assert_eq!(fri_evals.get(6).map(|v| v.len()), Some(128));
    assert_eq!(fri_evals.get(7).map(|v| v.len()), Some(64));
    assert_eq!(fri_evals.get(8).map(|v| v.len()), Some(32));
    assert_eq!(fri_evals.get(9).map(|v| v.len()), Some(16));
    assert_eq!(fri_evals.get(10).map(|v| v.len()), Some(8));
}
