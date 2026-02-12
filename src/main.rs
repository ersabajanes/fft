use fft::{fft, ifft, complex::Complex};
use core::f32;
use num_traits::Zero;

fn main() {
    let fs = 8;
    let ns = fs;

    let fr = 1.0;
    let am = 1.0;
    let ph = 0.0;
    let mut x = vec![Complex::zero(); ns];
    for n in 0..ns {
        x[n] = Complex::new(am * f32::cos(2.0 * f32::consts::PI * n as f32 / fs as f32 * fr + ph), 0.0);
    }

    let y =  fft(&x[..], 2);
    let z = ifft(&y[..], 2);

    for n in 0..ns {
        println!("{} | ({}) - ({}) = {}", n, x[n], z[n], x[n] - z[n]);
    }

    println!();
    for n in 0..ns {
        let m = (n + (ns + 1) / 2) % ns;
        println!("{} | {}", m, y[m]);
    }
}
