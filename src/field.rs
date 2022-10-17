use num_modular::ModularInteger;
use num_modular::MontgomeryInt;
use num_traits::pow::Pow;
use num_traits::Inv;

#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Gf<const P: u32>(MontgomeryInt<u32>);

impl<const P: u32> Gf<P> {
    pub fn new(v: u32) -> Self {
        Self(MontgomeryInt::new(v, &P))
    }

    pub fn residue(self) -> u32 {
        self.0.residue()
    }

    pub fn pow(self, rhs: u32) -> Self {
        Self(self.0.pow(rhs))
    }

    // Finds first multiplicative primitive element over F_P 
    pub fn generator() -> Self {
        // Collect unqiue prime factors
        let mut prime_factors = vec![];
        let mut p = P - 1;
        let mut it = 2;
        while p != 1 {
            if p % it == 0 {
                prime_factors.push(it);
            }
            while p % it == 0 {
                p /= it;
            }
            it += 1;
        }

        // Convert factors into exponents
        let exps: Vec<u32> = prime_factors
            .into_iter()
            .map(|factor| {
                (MontgomeryInt::new(P - 1, &P) * MontgomeryInt::new(factor, &P).inv()).residue()
            })
            .collect();

        // Test for all x^exp != 1 (mod p) for 2 < x < P, exp <- exps
        let mut found = 0;
        'l1: for x in 2..P {
            let x = MontgomeryInt::new(x, &P);
            for exp in &exps {
                if x.pow(*exp).residue() == 1 {
                    continue 'l1;
                }
            }
            found += 1;
            //println!("{}", x.residue());
            return Self(x);
        }
        unreachable!();
    }
}

impl<const P: u32> std::ops::Add for Gf<P> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl<const P: u32> std::ops::Mul for Gf<P> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}

#[test]
fn generator_test() {
    let g = Gf::<4391>::generator();
    let mut elems = [false; 4391];

    // Touch all elems once, except for elems[0]
    for i in 0..4390 {
        let e = g.pow(i).residue() as usize;
        assert_eq!(elems[e], false);
        elems[e] = true;
    }
    assert_eq!(elems[0], false);
}
