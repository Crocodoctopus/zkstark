use num_traits::identities::{One, Zero};
use num_traits::Inv;
use std::ops::{Add, Mul, MulAssign, Neg, Sub};

#[derive(Clone, Debug)]
pub struct Polynomial<T>(pub Vec<T>);

impl<T> Polynomial<T> {
    pub fn map<I>(self, f: impl Fn(T) -> I) -> Polynomial<I> {
        Polynomial(self.0.into_iter().map(|t| f(t)).collect())
    }

    pub fn new_degree1(c0: T) -> Self {
        Self(vec![c0])
    }

    pub fn new_degree2(c0: T, c1: T) -> Self {
        Self(vec![c1, c0])
    }
}

impl<T> Polynomial<T> {
    fn len(&self) -> usize {
        self.0.len()
    }
}

impl<T> Polynomial<T>
where
    T: Zero + PartialEq,
{
    pub fn degree(&self) -> usize {
        self.0
            .iter()
            .enumerate()
            .rev()
            .find_map(|(degree, coeffs)| (*coeffs == T::zero()).then(|| degree))
            .unwrap_or(0)
    }
}

impl<T> Mul<Self> for Polynomial<T>
where
    T: Mul<T, Output = T> + Add<T, Output = T> + Copy + Zero,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut out = vec![T::zero(); (self.0.len() - 1) + (rhs.0.len() - 1) + 1];
        for (degree0, coeff0) in self.0.iter().enumerate() {
            for (degree1, coeff1) in rhs.0.iter().enumerate() {
                out[degree0 + degree1] = out[degree0 + degree1] + *coeff0 * *coeff1;
            }
        }
        Polynomial(out)
    }
}

impl<T> Mul<T> for Polynomial<T>
where
    T: Mul<T, Output = T> + Copy,
{
    type Output = Self;
    fn mul(mut self, rhs: T) -> Self::Output {
        for coeff in self.0.iter_mut() {
            *coeff = *coeff * rhs;
        }
        self
    }
}

impl<T> MulAssign<T> for Polynomial<T>
where
    T: Mul<T, Output = T> + Copy,
{
    fn mul_assign(&mut self, rhs: T) {
        for coeff in self.0.iter_mut() {
            *coeff = *coeff * rhs;
        }
    }
}

impl<T> Add<Self> for Polynomial<T>
where
    T: Add<T, Output = T> + Copy,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let (short, mut long) = match (self.len(), rhs.len()) {
            (a, b) if a < b => (self, rhs),
            (a, b) => (rhs, self),
        };

        for i in 0..short.len() {
            long.0[i] = long.0[i] + short.0[i];
        }
        long
    }
}

pub fn lagrange<T>(points: Vec<(T, T)>) -> Vec<Polynomial<T>>
where
    T: Inv<Output = T>
        + Neg<Output = T>
        + One
        + Zero
        + Copy
        + Mul<T, Output = T>
        + Sub<T, Output = T>
        + PartialEq,
{
    // Generate non-normalized basis polynomials
    let mut bases: Vec<Polynomial<T>> = {
        // Generate left polynomial expansions
        let mut ll = vec![Polynomial::new_degree1(T::one()); points.len()];
        for i in 1..points.len() {
            ll[i] = ll[i - 1].clone() * Polynomial::new_degree2(T::one(), -points[i - 1].0);
        }

        // Generate right polynomial expansions
        let mut lr = vec![Polynomial::new_degree1(T::one()); points.len()];
        for i in (0..points.len() - 1).rev() {
            lr[i] = lr[i + 1].clone() * Polynomial::new_degree2(T::one(), -points[i + 1].0);
        }

        // Combine
        std::iter::zip(ll, lr).map(|(l, r)| l * r).collect()
    };

    bases
    //println!("{:?}", bases);

    // Collect normalizing factors in 'mod p'
    /*let nf: Vec<T> = points
        .iter()
        .map(|(x, _)| {
            points
                .iter()
                .filter(|(ix, _)| *ix != *x)
                .map(|(ix, _)| *x - *ix)
                .reduce(|acc, n| acc * n)
                .unwrap()
        })
        .map(|x: T| x.inv())
        .collect();

    // Apply factors to "normalize" bases
    for (i, poly) in bases.iter_mut().enumerate() {
        let y = points[i].1;
        *poly *= nf[i] * y;
    }

    // Add all bases together
    bases.into_iter().reduce(|acc, poly| acc + poly).unwrap()*/
}
