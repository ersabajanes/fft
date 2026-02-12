use std::f32::{INFINITY, NEG_INFINITY};

// use core::f32;
use fft::{complex::Complex};
use num_traits::Zero;

fn main() {
    let mut w: u32 = 0;
    let mut h: u32 = 0;
    let mut d: u8  = 0;
    let mut c: u8  = 0;
    let mut pix = qoi::load("noise_blue.qoi", &mut w, &mut h, &mut d, &mut c).unwrap();

    let mut xt = vec![Complex::zero(); w as usize * h as usize];
    for i in 0..h as usize {
        for j in 0..w as usize {
            xt[i * w as usize + j] = Complex::new(pix[(i * w as usize + j) * 3 + 0] as f32 / 255.0 * 2.0 - 1.0, 0.0);
        }
    }

    let xf = fft::fft2(&xt, w as usize, h as usize, 2);


    let mut min_a = INFINITY;
    let mut max_a = NEG_INFINITY;
    for i in 0..h as usize {
        for j in 0..w as usize {
            let f = xf[i * w as usize + j];
            let a = f32::sqrt(f.re * f.re + f.im * f.im);
            min_a = f32::min(a, min_a);
            max_a = f32::max(a, max_a);
        }
    }

    for i in 0..h as usize {
        for j in 0..w as usize {
            let iw = (i + (h as usize + 1) / 2) % h as usize;
            let jw = (j + (w as usize + 1) / 2) % w as usize;
            let f = xf[iw * w as usize + jw];

            // pix[(i * w as usize + j) * 3 + 0] = ((f.re - min_a) / (max_a - min_a) * 255.0) as u8;
            // pix[(i * w as usize + j) * 3 + 1] = ((f.im - min_a) / (max_a - min_a) * 255.0) as u8;
            // pix[(i * w as usize + j) * 3 + 2] = 0;

            let a = (f32::sqrt(f.re * f.re + f.im * f.im) - min_a) / (max_a - min_a) * 255.0;
            pix[(i * w as usize + j) * 3 + 0] = a as u8;
            pix[(i * w as usize + j) * 3 + 1] = a as u8;
            pix[(i * w as usize + j) * 3 + 2] = a as u8;
        }
    }

    qoi::save("fft2.qoi", w, h, d, c, &pix).unwrap();
}
