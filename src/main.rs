mod field;
mod merkle;
mod polynomial;

use num_traits::Pow;
use polynomial::{lagrange, Polynomial};

type F = field::Gf<3221225473>;

fn main() {
    ///////////////////
    // Part 1

    // The value of x_1
    let super_secret = 3141592;

    // Trace polynomial, containing x_0, x_1, all the way to x_1022
    let mut trace: Vec<F> = vec![F::from(1), F::from(super_secret)];
    for i in 2..1023 {
        let a = trace[i - 2] * trace[i - 2];
        let b = trace[i - 1] * trace[i - 1];
        trace.push(a + b);
    }

    // Assert the trace polynomial is correct
    assert!(trace.len() == 1023);
    assert!(trace[0].residue() == 1);
    assert!(trace[1].residue() == super_secret);
    assert!(trace[1022].residue() == 2338775057);

    // Generate a primitive root of F_3221225473
    let primitive_root = F::generator();

    // Create generators for 1024 and 8192 element cyclic groups
    let generator_g = primitive_root.pow(3145728);
    let generator_h = primitive_root.pow(393216);

    // Assert generators are of the correct order
    assert_eq!(generator_g.order(), 1024);
    assert_eq!(generator_h.order(), 8192);

    // Generate groups
    let g: Vec<F> = (0..1024).map(|n| generator_g.pow(n)).collect();
    let h: Vec<F> = (0..8192).map(|n| generator_h.pow(n)).collect();

    // Generate lagrange polynomial
    let poly = lagrange(&std::iter::zip(&g, &trace).map(|(&x, &y)| (x, y)).collect());

    // Check that the generated polynomial passes through each point in (G[i], trace[i]) | i < 1023
    for (x, y) in std::iter::zip(&g, &trace) {
        assert_eq!(poly.solve(*x), *y);
    }

    // Solve polynomial over H, shifted by the primitive root
    let eval: Vec<F> = h.iter().map(|n| poly.solve(primitive_root * *n)).collect();

    // Assert a few elements of eval
    assert_eq!(eval[0].residue(), 576067152);
    assert_eq!(eval[1].residue(), 3100214617);
    assert_eq!(eval[2].residue(), 2091264768);
    assert_eq!(eval[8189].residue(), 800520420);
    assert_eq!(eval[8190].residue(), 1199720174);
    assert_eq!(eval[8191].residue(), 1076821037);

    ///////////////////
    // Part 2

    // Constraint 0:
    // (f(x) - trace[0]) / (x - g[0])
    let poly0 = poly.clone() - Polynomial::from([trace[0]]);
    let poly1 = Polynomial::from([F::from(1), -g[0]]);
    let (c0, c0r) = Polynomial::rdiv(poly0, poly1);

    // Constraint 1:
    // (f(x) - trace[1022]) / (x - g[1022])
    let poly0 = poly.clone() - Polynomial::from([trace[1022]]);
    let poly1 = Polynomial::from([F::from(1), -g[1022]]);
    let (c1, c1r) = Polynomial::rdiv(poly0, poly1);

    // Check that constraints have no remainders
    assert!(c0r.degree().is_none());
    assert!(c1r.degree().is_none());
}
