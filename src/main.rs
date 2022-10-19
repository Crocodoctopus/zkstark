mod field;
mod polynomial;

use num_traits::Pow;

type F = field::Gf<3221225473>;

fn main() {
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

    // Domains
    let g_domain = (0..1023).map(|n| generator_g.pow(n));
    let h_domain = (0..8191).map(|n| generator_h.pow(n));

    // Generate lagrange polynomial
    let poly = polynomial::lagrange(
        &g_domain.clone().zip(&trace).map(|(x, y)| (x, *y)).collect()
    );

    // Check that the generated polynomial passes through each point in (G[i], trace[i]) | i < 1023
    g_domain.clone().zip(&trace).for_each(|(x, y)| assert_eq!(poly.solve(x), *y));

    // Solve polynomial over larger domain
    let eval: Vec<F> = h_domain
        .map(|n| generator_h * n)
        .map(|n| poly.solve(n))
        .collect();
}
