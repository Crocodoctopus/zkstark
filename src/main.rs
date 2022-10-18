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

    // Generator for a multiplicative group of 3221225472 elements
    let generator = F::generator();

    // Generator for a multiplicative group of 1024 elements
    let generator = generator.pow(3145728);

    // Ensure the generator is of the correct order
    assert!(generator.order() == 1024);

    // Generate 1024 cycle group group
    let group: Vec<(F, F)> = (0..1024).map(|n| (F::from(n), generator.pow(n))).collect();

    // Create interpolating polynomial from the 1024 cyclic group
    let poly = polynomial::lagrange(&group);

    // Assert that the lagrange poly hits every points in 1024 cyclic group
    for (x, y) in group {
        assert_eq!(poly.solve(x), y);
    }
}
