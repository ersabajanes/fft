pub mod complex;

use complex::Complex;
use core::f32;
use num_traits::Zero;

fn next_pow(num: usize, radix: usize) -> usize {
    let mut foo = num - 1;
    let mut bar = 1;

    while foo > 0 {
        foo /= radix;
        bar *= radix;
    }

    bar
}

pub fn fft(xt: &[Complex<f32>], radix: usize) -> Vec<Complex<f32>> {
    let ns = xt.len();
    if ns == 1 {
        return xt.to_vec();
    }

    let mut xf = vec![Complex::zero(); ns];
    for i in 0..radix {
        let xt0 = Vec::from_iter(xt.iter().skip(i).step_by(radix).cloned());
        let xf0 = fft(&xt0, radix);

        for m in 0..ns {
            let m2 = m % (ns / radix);
            let th = -f32::consts::TAU * m as f32 / ns as f32 * i as f32;
            xf[m] += Complex::new(f32::cos(th), f32::sin(th)) * xf0[m2];
        }
    }

    xf
}

pub fn fft_it(xt: &[Complex<f32>], radix: usize) -> Vec<Complex<f32>> {
    let ns = xt.len();
    if ns == 1 {
        return xt.to_vec();
    }

    let mut xfa = xt.to_vec();
    let mut xfb = vec![Complex::zero(); ns];

    let mut foo;
    let mut bar = ns;
    while bar > 1 {
        foo = bar;
        bar /= radix;

        for n in 0..ns {
            xfb[(n / foo) * foo + (n % radix) * bar + (n / radix) % bar] = xfa[n];
        }

        (xfa, xfb) = (xfb, xfa);
    }

    let mut foo;
    let mut bar = 1;
    while bar < ns {
        foo = bar;
        bar *= radix;

        for m in 0..ns {
            xfb[m] = Complex::zero();
            for i in 0..radix {
                let m2 = (m / bar) * bar + (m % foo) + i * foo;
                let th = -f32::consts::TAU * (m % bar) as f32 / bar as f32 * i as f32;
                xfb[m] += Complex::new(f32::cos(th), f32::sin(th)) * xfa[m2];
            }
        }

        (xfa, xfb) = (xfb, xfa);
    }

    xfa
}

pub fn ifft(xf: &[Complex<f32>], radix: usize) -> Vec<Complex<f32>> {
    let ns = xf.len();
    let mut xt = vec![Complex::zero(); ns];

    if ns == 1 {
        xt[0] += xf[0];
        return xt;
    }

    for i in 0..radix {
        let xf0 = Vec::from_iter(xf.iter().skip(i).step_by(radix).cloned());
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

pub fn fft2(xt: &[Complex<f32>], w: usize, h: usize, r: usize) -> Vec<Complex<f32>> {
    let w2 = next_pow(w, r);
    let h2 = next_pow(h, r);

    let mut xf = vec![Complex::zero(); w2 * h2];

    for i in 0..h {
        let mut xti = xt[i * w..i * w + w].to_vec();
        xti.resize(w2, Complex::zero());

        let xfi = fft(&xti, r);
        for j in 0..w2 {
            xf[i * w2 + j] = xfi[j];
        }
    }

    for j in 0..w2 {
        let xtj = Vec::from_iter(xf.iter().skip(j).step_by(w2).cloned());

        let xfj = fft(&xtj, r);
        for i in 0..h2 {
            xf[i * w2 + j] = xfj[i];
        }
    }

    xf
}

pub fn ifft2(xt: &[Complex<f32>], w: usize, h: usize, r: usize) -> Vec<Complex<f32>> {
    let mut xf = vec![Complex::zero(); xt.len()];

    for i in 0..h {
        let xi = ifft(&xt[i * w..i * w + w], r);
        for j in 0..w {
            xf[i * w + j] = xi[j];
        }
    }

    for j in 0..w {
        let xj = Vec::from_iter(xf.iter().skip(j).step_by(w).cloned());
        let xj = ifft(&xj, r);
        for i in 0..h {
            xf[i * w + j] = xj[i];
        }
    }

    xf
}
