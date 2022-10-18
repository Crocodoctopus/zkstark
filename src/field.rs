use num_modular::ModularInteger;
use num_modular::MontgomeryInt;
use num_traits::pow::Pow;
use num_traits::Inv;
use std::ops::Neg;

#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Gf<const P: u32>(pub MontgomeryInt<u32>);

impl<const P: u32> From<i32> for Gf<P> {
    fn from(f: i32) -> Self {
        Self(if f < 0 {
            MontgomeryInt::new(f.abs() as u32, &P).neg()
        } else {
            MontgomeryInt::new(f as u32, &P)
        })
    }
}

impl<const P: u32> From<u32> for Gf<P> {
    fn from(f: u32) -> Self {
        Self(MontgomeryInt::new(f, &P))
    }
}

impl<const P: u32> Pow<u32> for Gf<P> {
    type Output = Self;
    fn pow(self, rhs: u32) -> Self::Output {
        Self(self.0.pow(rhs))
    }
}

impl<const P: u32> Gf<P> {
    pub fn residue(self) -> u32 {
        self.0.residue()
    }

    pub fn order(self) -> u32 {
        (1..)
            .find_map(|it| (self.pow(it).residue() == 1).then(|| it))
            .unwrap()
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
        'l1: for x in 2..P {
            let x = MontgomeryInt::new(x, &P);
            for exp in &exps {
                if x.pow(*exp).residue() == 1 {
                    continue 'l1;
                }
            }
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

impl<const P: u32> num_traits::identities::Zero for Gf<P> {
    fn zero() -> Self {
        Self(MontgomeryInt::new(0, &P))
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<const P: u32> num_traits::identities::One for Gf<P> {
    fn one() -> Self {
        Self(MontgomeryInt::new(1, &P))
    }

    fn is_one(&self) -> bool {
        self.0.is_zero()
    }
}

impl<const P: u32> std::ops::Neg for Gf<P> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self(self.0.neg())
    }
}

impl<const P: u32> std::ops::Sub for Gf<P> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self(self.0 - rhs.0)
    }
}

impl<const P: u32> Inv for Gf<P> {
    type Output = Self;
    fn inv(self) -> Self::Output {
        Self(self.0.inv())
    }
}

#[test]
fn generator_test() {
    let g = Gf::<4391>::generator();
    assert_eq!(g.order(), 4390);

    // Touch all elems once, except for elems[0]
    let mut elems = [false; 4391];
    for i in 0..4390 {
        let e = g.pow(i).residue() as usize;
        assert_eq!(elems[e], false);
        elems[e] = true;
    }
    assert_eq!(elems[0], false);
}
