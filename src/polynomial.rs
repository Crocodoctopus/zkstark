use num_traits::identities::{One, Zero};
use num_traits::{Inv, Pow};
use std::ops::{Add, Mul, MulAssign, Neg, Sub};

#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial<T>(Vec<T>);

impl<T> Polynomial<T> {
    pub fn map<I>(self, f: impl Fn(T) -> I) -> Polynomial<I> {
        Polynomial(self.0.into_iter().map(|t| f(t)).collect())
    }

    pub fn coeff(&self, i: usize) -> &T {
        &self.0[i]
    }

    pub fn coeff_mut(&mut self, i: usize) -> &mut T {
        &mut self.0[i]
    }

    fn len(&self) -> usize {
        self.0.len()
    }
}

impl<T, II> From<II> for Polynomial<T>
where
    II: IntoIterator<Item = T>,
    II::IntoIter: DoubleEndedIterator,
{
    fn from(f: II) -> Self {
        Self(f.into_iter().rev().collect())
    }
}

impl<T> Polynomial<T>
where
    T: Add<T, Output = T> + Mul<T, Output = T> + Pow<u32, Output = T> + Copy + Zero,
{
    pub fn solve(&self, t: T) -> T {
        self.0
            .iter()
            .enumerate()
            .fold(T::zero(), |acc, (degree, coeff)| {
                acc + *coeff * t.pow(degree as u32)
            })
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
            _ => (rhs, self),
        };

        for i in 0..short.len() {
            long.0[i] = long.0[i] + short.0[i];
        }
        long
    }
}

pub fn lagrange<T>(points: &Vec<(T, T)>) -> Polynomial<T>
where
    T: Inv<Output = T>
        + Neg<Output = T>
        + One
        + Zero
        + Copy
        + Mul<T, Output = T>
        + Sub<T, Output = T>
        + Pow<u32, Output = T>
        + PartialEq
{
    // Generate non-normalized basis polynomials
    let mut bases: Vec<Polynomial<T>> = {
        // Generate left polynomial expansions
        let mut ll = vec![Polynomial::from([T::one()]); points.len()];
        for i in 1..points.len() {
            ll[i] = ll[i - 1].clone() * Polynomial::from([T::one(), -points[i - 1].0]);
        }

        // Generate right polynomial expansions
        let mut lr = vec![Polynomial::from([T::one()]); points.len()];
        for i in (0..points.len() - 1).rev() {
            lr[i] = lr[i + 1].clone() * Polynomial::from([T::one(), -points[i + 1].0]);
        }

        // Combine
        std::iter::zip(ll, lr).map(|(l, r)| l * r).collect()
    };

    // Normalize
    for (poly, &(tx, ty)) in bases.iter_mut().zip(points) {
        *poly *= poly.solve(tx).inv() * ty;
    }

    // Add all bases together
    bases.into_iter().reduce(|acc, poly| acc + poly).unwrap()
}

#[test]
fn lagrange_test() {
    type F = crate::field::Gf<7>;

    // Pick 4 points
    let p0 = (F::from(0), F::from(2));
    let p1 = (F::from(1), F::from(4));
    let p2 = (F::from(2), F::from(3));
    let p3 = (F::from(3), F::from(1));

    // Generate lagrange for those 4 points
    let poly = lagrange(vec![p0, p1, p2, p3]);

    // Solve for the remaining points in GF7
    let p4x = F::from(4);
    let p4 = (p4x, poly.solve(p4x));
    let p5x = F::from(5);
    let p5 = (p5x, poly.solve(p5x));
    let p6x = F::from(6);
    let p6 = (p6x, poly.solve(p6x));

    // Assert the lagrange for combinations of points are equal
    assert_eq!(poly, lagrange(vec![p0, p3, p5, p6]));
    assert_eq!(poly, lagrange(vec![p1, p6, p3, p2]));
    assert_eq!(poly, lagrange(vec![p3, p2, p1, p0]));
    assert_eq!(poly, lagrange(vec![p6, p5, p4, p3]));
}
