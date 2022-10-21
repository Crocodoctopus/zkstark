use num_traits::identities::{One, Zero};
use num_traits::{Inv, Pow};
use std::ops::{Add, Div, Mul, MulAssign, Neg, Sub};
use std::ops::{Index, IndexMut};

pub fn x<T>(t: T, e: usize) -> Polynomial<T>
where
    T: Zero + One + Clone,
{
    let mut temp = vec![T::zero(); e + 1];
    temp[0] = t;
    Polynomial::from(temp)
}

fn reduce<T>(mut v: Vec<T>) -> Vec<T>
where
    T: Zero + PartialEq,
{
    let l = v
        .iter()
        .enumerate()
        .rev()
        .find_map(|(degree, coeffs)| (*coeffs != T::zero()).then(|| degree))
        .map(|v| v + 1)
        .unwrap_or(0);
    v.truncate(l);
    return v;
}

#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial<T>(Box<[T]>);

impl<T, II> From<II> for Polynomial<T>
where
    II: IntoIterator<Item = T>,
    II::IntoIter: DoubleEndedIterator,
{
    fn from(f: II) -> Self {
        Self(f.into_iter().rev().collect())
    }
}

// Hmm, not sure about this one
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

impl<T> Polynomial<T> {
    pub fn degree(&self) -> Option<usize> {
        if self.0.len() == 0 {
            None
        } else {
            Some(self.0.len() - 1)
        }
    }
}

impl<T> Polynomial<T>
where
    for<'a> &'a T: Pow<u32, Output = T>,
    T: MulAssign<T>,
{
    pub fn apply_const(mut self, t: T) -> Self {
        for (degree, coeff) in self.0.iter_mut().enumerate() {
            *coeff *= (&t).pow(degree as u32)
        }
        self
    }
}

impl<T> Index<usize> for Polynomial<T> {
    type Output = T;
    fn index(&self, i: usize) -> &Self::Output {
        let degree = self.degree().unwrap();
        &self.0[degree - i]
    }
}

impl<T> IndexMut<usize> for Polynomial<T> {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        let degree = self.degree().unwrap();
        &mut self.0[degree - i]
    }
}

impl<T> Add for &Polynomial<T>
where
    for<'a> &'a T: Add<Output = T>,
    T: Zero + Clone + PartialEq,
{
    type Output = Polynomial<T>;
    fn add(self, rhs: Self) -> Self::Output {
        let degree0 = self.degree();
        let degree1 = rhs.degree();

        // If either is None, return the other
        let (degree0, degree1) = match (degree0, degree1) {
            (Some(d0), Some(d1)) => (d0, d1),
            (None, None) => return Polynomial(Box::new([])),
            (None, _) => return rhs.clone(),
            (_, None) => return self.clone(),
        };

        // Construct poly
        let mut poly = vec![T::zero(); usize::max(degree0, degree1) + 1];

        // Perform addition
        for (i, t) in self.0.iter().enumerate() {
            poly[i] = &poly[i] + t;
        }
        for (i, t) in rhs.0.iter().enumerate() {
            poly[i] = &poly[i] + t;
        }

        Polynomial(reduce(poly).into_boxed_slice())
    }
}

impl<T> Add for Polynomial<T>
where
    for<'a> &'a T: Add<Output = T> + Sub<Output = T>,
    T: Zero + Clone + PartialEq,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl<T> Sub for &Polynomial<T>
where
    for<'a> &'a T: Add<Output = T> + Sub<Output = T>,
    T: Zero + Clone + PartialEq,
{
    type Output = Polynomial<T>;
    fn sub(self, rhs: Self) -> Self::Output {
        let degree0 = self.degree();
        let degree1 = rhs.degree();

        // If either is None, return the other
        let (degree0, degree1) = match (degree0, degree1) {
            (Some(d0), Some(d1)) => (d0, d1),
            (None, None) => return Polynomial(Box::new([])),
            (None, _) => return rhs.clone(),
            (_, None) => return self.clone(),
        };

        // Construct poly
        let mut poly = vec![T::zero(); usize::max(degree0, degree1) + 1];

        // Perform addition
        for (i, t) in self.0.iter().enumerate() {
            poly[i] = &poly[i] + t;
        }
        for (i, t) in rhs.0.iter().enumerate() {
            poly[i] = &poly[i] - t;
        }

        Polynomial(reduce(poly).into_boxed_slice())
    }
}

impl<T> Sub for Polynomial<T>
where
    for<'a> &'a T: Add<Output = T> + Sub<Output = T>,
    T: Zero + Clone + PartialEq,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        <&Self as Sub<&Self>>::sub(&self, &rhs)
    }
}

impl<T> Mul for &Polynomial<T>
where
    for<'a> &'a T: Mul<Output = T> + Add<Output = T>,
    T: Zero + Clone,
{
    type Output = Polynomial<T>;
    fn mul(self, rhs: Self) -> Self::Output {
        let degree0 = self.degree();
        let degree1 = rhs.degree();

        // If either degree is None, return "None"
        let (degree0, degree1) = match (degree0, degree1) {
            (Some(d0), Some(d1)) => (d0, d1),
            (_, _) => return Polynomial(Box::new([])),
        };

        // Allocate space for new poly
        let new_degree = degree0 + degree1;
        let mut poly = vec![T::zero(); new_degree + 1];

        // Perform multiplication
        for (degree0, coeff0) in self.0.iter().enumerate() {
            for (degree1, coeff1) in rhs.0.iter().enumerate() {
                poly[degree0 + degree1] = &poly[degree0 + degree1] + &(coeff0 * coeff1);
            }
        }

        Polynomial(poly.into_boxed_slice())
    }
}

impl<T> Mul for Polynomial<T>
where
    for<'a> &'a T: Mul<Output = T> + Add<Output = T>,
    T: Zero + Clone,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
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

impl<T> Polynomial<T>
where
    for<'a> &'a T: Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + Add<Output = T>,
    T: Zero + PartialEq + Clone,
{
    pub fn rdiv(lhs: Self, rhs: Self) -> (Self, Self) {
        // Get degree of each poly
        let lhs_degree = lhs.degree().unwrap_or(0);
        let rhs_degree = rhs.degree().unwrap_or(0);

        // Return early if division is undoable
        if lhs_degree < rhs_degree {
            return (Polynomial::from([]), lhs);
        }

        // Get leading coeff
        let lhs_lead = &lhs[0];
        let rhs_lead = &rhs[0];

        // Construct division poly
        let diff = lhs_degree - rhs_degree;
        let mut div = Polynomial::from(vec![T::zero(); diff + 1]);
        div[0] = lhs_lead / rhs_lead;

        // Calculate remainder
        let r = lhs - &div * &rhs;

        // Reapply division on remainder
        let (div2, r) = Polynomial::<T>::rdiv(r, rhs);

        // Return
        return (div + div2, r);
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
        + PartialEq,
    for<'a> &'a T: Mul<Output = T> + Add<Output = T>,
{
    // Generate non-normalized basis polynomials
    let mut bases: Vec<Polynomial<T>> = {
        // Generate left polynomial expansions
        let mut ll = vec![Polynomial::from([T::one()]); points.len()];
        for i in 1..points.len() {
            ll[i] = &ll[i - 1] * &Polynomial::from([T::one(), -points[i - 1].0]);
        }

        // Generate right polynomial expansions
        let mut lr = vec![Polynomial::from([T::one()]); points.len()];
        for i in (0..points.len() - 1).rev() {
            lr[i] = &lr[i + 1] * &Polynomial::from([T::one(), -points[i + 1].0]);
        }

        // Combine
        std::iter::zip(ll, lr).map(|(l, r)| &l * &r).collect()
    };

    // Normalize
    for (poly, &(tx, ty)) in bases.iter_mut().zip(points) {
        *poly *= poly.solve(tx).inv() * ty;
    }

    // Add all bases together
    bases.into_iter().reduce(|acc, poly| &acc + &poly).unwrap()
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
    let poly = lagrange::<F>(&vec![p0, p1, p2, p3]);

    // Solve for the remaining points in GF7
    let p4x = F::from(4);
    let p4 = (p4x, poly.solve(p4x));
    let p5x = F::from(5);
    let p5 = (p5x, poly.solve(p5x));
    let p6x = F::from(6);
    let p6 = (p6x, poly.solve(p6x));

    // Assert the lagrange for combinations of points are equal
    assert_eq!(poly, lagrange::<F>(&vec![p0, p3, p5, p6]));
    assert_eq!(poly, lagrange::<F>(&vec![p1, p6, p3, p2]));
    assert_eq!(poly, lagrange::<F>(&vec![p3, p2, p1, p0]));
    assert_eq!(poly, lagrange::<F>(&vec![p6, p5, p4, p3]));
}

#[test]
fn rdiv_test() {
    // Pick two polynomials
    let p0 = Polynomial::from([1, -3, -10]); // x² -3x -10
    let p1 = Polynomial::from([1, 2]); // x +2

    // Perform rdiv
    let (d, r) = Polynomial::<i32>::rdiv(p0, p1);

    // Assert
    assert_eq!(d, Polynomial::from([1, -5])); // x -5
    assert_eq!(r, Polynomial::from([])); // 0

    // Pick two more polynomials
    let p0 = Polynomial::from([2, -5, -1]); // 2x² -5x -1
    let p1 = Polynomial::from([1, -3]); // x -3

    // Perform rdiv
    let (d, r) = Polynomial::<i32>::rdiv(p0, p1);

    // Assert
    assert_eq!(d, Polynomial::from([2, 1])); // 2x +1
    assert_eq!(r, Polynomial::from([2])); // 2

    // Last two
    let p0 = Polynomial::from([1, 0, 2, 0, 0, 6, -9]); // x^6 +2x^4 +6x -9
    let p1 = Polynomial::from([1, 0, 0, 3]); // x^3 +3

    // Perform rdiv
    let (d, r) = Polynomial::<i32>::rdiv(p0, p1);

    // Assert
    assert_eq!(d, Polynomial::from([1, 0, 2, -3])); // x^3 +2x -3
    assert_eq!(r, Polynomial::from([])); // 0
}
