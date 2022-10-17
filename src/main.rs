mod field;

type F = field::Gf<3221225473>;

fn main() {
    let super_secret = 3141592;

    let mut trace: Vec<F> = vec![F::new(1), F::new(super_secret)];
    for i in 2..1023 {
        let a = trace[i - 2] * trace[i - 2];
        let b = trace[i - 1] * trace[i - 1];
        trace.push(a + b);
    }

    assert!(trace.len() == 1023);
    assert!(trace[0].residue() == 1);
    assert!(trace[1].residue() == super_secret);
    assert!(trace[1022].residue() == 2338775057);

    //
    let generator = F::generator();
    let group: Vec<F> = (0..1024).map(|n| generator.pow(n)).collect();
}
