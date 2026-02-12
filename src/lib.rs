pub mod complex;

use complex::Complex;
use core::f32;
use num_traits::Zero;

pub fn fft(xt: &[Complex<f32>], radix: usize) -> Vec<Complex<f32>> {
    let ns = xt.len();
    let mut xf = vec![Complex::zero(); ns];

    if ns == 1 {
        xf[0] += xt[0];
        return xf;
    }

    for i in 0..radix {
        let xt0: Vec<Complex<f32>> = Vec::from_iter(xt.iter().skip(i).step_by(radix).cloned());
        let xf0 = fft(&xt0, radix);

        for m in 0..ns {
            let m2 = m % (ns / radix);
            let th = -f32::consts::TAU * m as f32 / ns as f32 * i as f32;
            xf[m] += Complex::new(f32::cos(th), f32::sin(th)) * xf0[m2];
        }
    }

    xf
}

pub fn ifft(xf: &[Complex<f32>], radix: usize) -> Vec<Complex<f32>> {
    let ns = xf.len();
    let mut xt = vec![Complex::zero(); ns];

    if ns == 1 {
        xt[0] += xf[0];
        return xt;
    }

    for i in 0..radix {
        let xf0: Vec<Complex<f32>> = Vec::from_iter(xf.iter().skip(i).step_by(radix).cloned());
        let xt0 = ifft(&xf0, radix);

        for n in 0..ns {
            let m2 = n % (ns / radix);
            let th = f32::consts::TAU * n as f32 / ns as f32 * i as f32;
            xt[n] +=
                Complex::new(f32::cos(th) / radix as f32, f32::sin(th) / radix as f32) * xt0[m2];
        }
    }

    xt
}
