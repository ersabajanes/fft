use num_traits::{Num, Zero};

#[derive(Clone, Copy, Debug)]
pub struct Complex<T: Num + Copy> {
    re: T,
    im: T,
}

impl<T: Num + Copy> Complex<T> {
    pub fn new(re: T, im: T) -> Self {
        Self { re: re, im: im }
    }
}

impl<T: Num + Copy> Zero for Complex<T> {
    fn zero() -> Self {
        Self { re: T::zero(), im: T::zero() }
    }

    fn is_zero(&self) -> bool {
        T::is_zero(&self.re) && T::is_zero(&self.im)
    }
}

impl<T: Num + Copy> std::fmt::Display for Complex<T>
where
    T: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.3} + {:.3}i", &self.re, &self.im)
    }
}

impl<T: Num + Copy> std::ops::Add for Complex<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl<T: Num + Copy> std::ops::AddAssign for Complex<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<T: Num + Copy> std::ops::Sub for Complex<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl<T: Num + Copy> std::ops::SubAssign for Complex<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<T: Num + Copy> std::ops::Mul for Complex<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl<T: Num + Copy> std::ops::MulAssign for Complex<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<T: Num + Copy> std::ops::Div for Complex<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        let dot = rhs.re * rhs.re + rhs.im * rhs.im;
        Self {
            re: (self.re * rhs.re + self.im * rhs.im) / dot,
            im: (self.im * rhs.re - self.re * rhs.im) / dot,
        }
    }
}

impl<T: Num + Copy> std::ops::DivAssign for Complex<T> {
    #[inline]
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}
