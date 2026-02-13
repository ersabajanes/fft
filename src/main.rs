use fft::{complex::Complex};
use num_traits::Zero;

fn next_pow2(num: usize) -> usize {    
    let mut foo = num - 1;
    let mut bar = 1;

    while foo > 0 {
        foo >>= 1;
        bar <<= 1;
    }

    bar
}

fn main() {
    // Image
    let mut w0: u32 = 0;
    let mut h0: u32 = 0;
    let mut d0: u8  = 0;
    let mut c0: u8  = 0;
    let img = qoi::load("IMG-20260127-WA0001.qoi", &mut w0, &mut h0, &mut d0, &mut c0).unwrap();

    // Kernel
    let mut w1: u32 = 0;
    let mut h1: u32 = 0;
    let mut d1: u8  = 0;
    let mut c1: u8  = 0;
    let ker = qoi::load("flower.qoi", &mut w1, &mut h1, &mut d1, &mut c1).unwrap();

    // FFT size
    let w2: usize = w0 as usize + w1 as usize;
    let h2: usize = h0 as usize + h1 as usize;

    // Image to Complex
    let mut rt0 = vec![Complex::zero(); w2 * h2];
    let mut gt0 = vec![Complex::zero(); w2 * h2];
    let mut bt0 = vec![Complex::zero(); w2 * h2];
    for i in 0..h0 as usize {
        for j in 0..w0 as usize {
            let r = img[(i * w0 as usize + j) * d0 as usize + 0] as f32 / 255.0;
            let g = img[(i * w0 as usize + j) * d0 as usize + 1] as f32 / 255.0;
            let b = img[(i * w0 as usize + j) * d0 as usize + 2] as f32 / 255.0;
            let a = if d0 > 3 { img[(i * w0 as usize + j) * d0 as usize + 3] as f32 / 255.0 } else { 1.0 };

            rt0[i * w2 + j] = Complex::new((1.0 - r) * a, 0.0);
            gt0[i * w2 + j] = Complex::new((1.0 - g) * a, 0.0);
            bt0[i * w2 + j] = Complex::new((1.0 - b) * a, 0.0);
        }
    }

    // Kernel to complex
    let mut rt1 = vec![Complex::zero(); w2 * h2];
    let mut gt1 = vec![Complex::zero(); w2 * h2];
    let mut bt1 = vec![Complex::zero(); w2 * h2];
    for i in 0..h1 as usize {
        for j in 0..w1 as usize {
            rt1[i * w2 + j] = Complex::new(ker[(i * w1 as usize + j) * d1 as usize + 0] as f32 / 255.0, 0.0);
            gt1[i * w2 + j] = Complex::new(ker[(i * w1 as usize + j) * d1 as usize + 1] as f32 / 255.0, 0.0);
            bt1[i * w2 + j] = Complex::new(ker[(i * w1 as usize + j) * d1 as usize + 2] as f32 / 255.0, 0.0);
        }
    }

    // Image FFT
    let rf0 = fft::fft2(&rt0, w2, h2, 2);
    let gf0 = fft::fft2(&gt0, w2, h2, 2);
    let bf0 = fft::fft2(&bt0, w2, h2, 2);

    // Kernel FFT
    let rf1 = fft::fft2(&rt1, w2, h2, 2);
    let gf1 = fft::fft2(&gt1, w2, h2, 2);
    let bf1 = fft::fft2(&bt1, w2, h2, 2);

    let w2: usize = next_pow2(w2 as usize);
    let h2: usize = next_pow2(h2 as usize);

    let mut rf2 = vec![Complex::zero(); w2 * h2];
    let mut gf2 = vec![Complex::zero(); w2 * h2];
    let mut bf2 = vec![Complex::zero(); w2 * h2];

    // Convolution
    for i in 0..h2 as usize {
        for j in 0..w2 as usize {
            rf2[i * w2 + j] = rf0[i * w2 + j] * rf1[i * w2 + j];
            gf2[i * w2 + j] = gf0[i * w2 + j] * gf1[i * w2 + j];
            bf2[i * w2 + j] = bf0[i * w2 + j] * bf1[i * w2 + j];
        }
    }

    // IFFT
    let rt2 = fft::ifft2(&rf2, w2, h2, 2);
    let gt2 = fft::ifft2(&gf2, w2, h2, 2);
    let bt2 = fft::ifft2(&bf2, w2, h2, 2);

    let rts = rt1.iter().cloned().reduce(|a, b| a + b).unwrap();
    let gts = gt1.iter().cloned().reduce(|a, b| a + b).unwrap();
    let bts = bt1.iter().cloned().reduce(|a, b| a + b).unwrap();

    // let w0: u32 = w2 as u32;
    // let h0: u32 = h2 as u32;

    // Result
    let mut pix = vec![0u8; w0 as usize * h0 as usize * 3];
    for i in 0..h0 as usize {
        for j in 0..w0 as usize {
            let r = rt2[(i + h1 as usize / 2) * w2 + (j + w1 as usize / 2)];
            let g = gt2[(i + h1 as usize / 2) * w2 + (j + w1 as usize / 2)];
            let b = bt2[(i + h1 as usize / 2) * w2 + (j + w1 as usize / 2)];

            let r = r.re / rts.re;
            let g = g.re / gts.re;
            let b = b.re / bts.re;

            let r = r * 255.0;
            let g = g * 255.0;
            let b = b * 255.0;

            // let r = r.pow(1.0 / 2.2);
            // let g = g.pow(1.0 / 2.2);
            // let b = b.pow(1.0 / 2.2);

            // let r = r.re * r.re + r.im * r.im;
            // let g = g.re * g.re + g.im * g.im;
            // let b = b.re * b.re + b.im * b.im;

            pix[(i * w0 as usize + j) * 3 + 0] = r as u8;
            pix[(i * w0 as usize + j) * 3 + 1] = g as u8;
            pix[(i * w0 as usize + j) * 3 + 2] = b as u8;
        }
    }

    qoi::save("IMG-20260127-WA0001_fft.qoi", w0, h0, 3, 0, &pix).unwrap();
}
