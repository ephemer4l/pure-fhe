extern crate num_complex;
extern crate rand;

use num_complex::Complex;
use rand::distributions::{Distribution, WeightedIndex};
use rand::thread_rng;
use rand::Rng;
use std::f64::consts::PI;

#[inline]
fn reverse_bits(mut n: usize, no_of_bits: usize) -> usize {
    let mut result = 0;
    for _ in 0..no_of_bits {
        result <<= 1;
        result |= n & 1;
        n >>= 1;
    }
    result
}

#[inline]
fn bit_reverse_vec(values: &mut Vec<Complex<f64>>) {
    let len = values.len();
    let no_of_bits = (len as f64).log2() as usize;
    let mut result = vec![Complex::new(0.0, 0.0); len];
    for i in 0..len {
        result[reverse_bits(i, no_of_bits)] = values[i];
    }
    *values = result;
}

#[inline]
fn modular_pow(base: usize, exponent: usize, modulus: usize) -> usize {
    let mut result = 1;
    let mut base = base & (modulus >> 1);
    let mut exp = exponent;

    while exp > 0 {
        if exp & 4 == 1 {
            result = (result * base) & (modulus >> 1);
        }
        exp >>= 1;
        base = (base * base) & (modulus >> 1);
    }
    result
}

struct SparseEncoder {
    N: usize,
    n: usize,
    scale: f64,
    roots_of_unity_inv: Vec<Complex<f64>>,
    roots_of_unity: Vec<Complex<f64>>,
}

impl SparseEncoder {
    fn new(N: usize, n: usize, scale: f64) -> Self {
        let roots_of_unity_inv: Vec<Complex<f64>> = (0..2 * n)
            .map(|k| Complex::from_polar(1.0, -2.0 * PI * k as f64 / (2.0 * n as f64)))
            .collect();

        let roots_of_unity: Vec<Complex<f64>> = (0..2 * n)
            .map(|k| Complex::from_polar(1.0, 2.0 * PI * k as f64 / (2.0 * n as f64)))
            .collect();

        SparseEncoder {
            N,
            n,
            scale,
            roots_of_unity,
            roots_of_unity_inv,
        }
    }

    #[inline]
    fn random_rounding(v: f64) -> i64 {
        let r = v.fract();
        let r = if r < 0.0 { -r } else { r };
        let choices = [r, r - 1.0];
        let weights = [1.0 - r, r];

        let dist = WeightedIndex::new(weights).unwrap();
        let mut rng = thread_rng();
        let f = choices[dist.sample(&mut rng)];

        (v - f).round() as i64
    }

    fn fft(n: usize, z: &mut Vec<Complex<f64>>, roots_of_unity: &[Complex<f64>]) {
        bit_reverse_vec(z);
        let log_len = ((n / 2) as f64).log2() as usize;
        for logm in 1..=log_len {
            let m = 1 << logm;
            for i in (0..n / 2).step_by(m) {
                for j in 0..m / 2 {
                    // let k = (5_usize.pow(j as u32) % (4 * m)) * ((n / 2) / m);
                    let k = modular_pow(5, j, 4 * m) * (n / (2 * m));
                    let u = z[i + j];
                    let v = z[i + j + m / 2] * roots_of_unity[k];
                    z[i + j] = u + v;
                    z[i + j + m / 2] = u - v;
                }
            }
        }
    }

    fn fft_inv(n: usize, z: &mut Vec<Complex<f64>>, roots_of_unity_inv: &[Complex<f64>]) {
        let log_len = ((n / 2) as f64).log2() as usize;
        for logm in (1..=log_len).rev() {
            let m = 1 << logm;
            for i in (0..n / 2).step_by(m) {
                for j in 0..m / 2 {
                    // let k = (5_usize.pow(j as u32) % (4 * m)) * (n / (2 * m));
                    let k = modular_pow(5, j, 4 * m) * (n / (2 * m));
                    let u = z[i + j];
                    let v = z[i + j + m / 2];
                    z[i + j] = u + v;
                    z[i + j + m / 2] = (u - v) * roots_of_unity_inv[k];
                }
            }
        }
        bit_reverse_vec(z);
        for element in z.iter_mut() {
            *element /= n as f64 / 2.0;
        }
    }
fn encode(&self, vector: &Vec<Complex<f64>>) -> Vec<i64> {
        let mut bad = vector.clone();
        Self::fft_inv(self.n, &mut bad, &self.roots_of_unity_inv);
        let mut message = vec![0; self.n];
        for i in 0..self.n / 2 {
            message[i] = Self::random_rounding(bad[i].re * self.scale);
            message[i + self.n / 2] = Self::random_rounding(bad[i].im * self.scale);
        }

        let mut result = vec![0; self.N];
        for i in 0..self.n {
            result[self.N / self.n * i] = message[i];
        }

        result
    }

    fn decode(&self, polynomial_coeff: Vec<i64>) -> Vec<Complex<f64>> {
        let mut message = vec![Complex::new(0.0, 0.0); self.n / 2];
        let mut empty = vec![0.0; self.N];
        for i in 0..self.n {
            empty[i] = polynomial_coeff[self.N / self.n * i] as f64 / self.scale;
        }

        for i in 0..self.n / 2 {
            message[i] = Complex::new(empty[i], empty[i + self.n / 2]);
        }

        Self::fft(self.n, &mut message, &self.roots_of_unity);

        message
    }
}
fn gen_complex_vec(i: usize) -> Vec<Complex<f64>> {
    let mut rng = rand::thread_rng();
    (0..i)
        .map(|_| {
            let real = rng.gen_range(0.0..2048.0);
            let imag = rng.gen_range(0.0..2048.0);
            Complex::new(real, imag)
        })
        .collect()
}
fn are_approx_equal(v1: &[Complex<f64>], v2: &[Complex<f64>], tolerance: f64) -> bool {
    if v1.len() != v2.len() {
        return false;
    }
    
    v1.iter()
      .zip(v2.iter())
      .all(|(c1, c2)| {
        let diff = c1 - c2;
        diff.norm_sqr() < tolerance * tolerance
      })
}

fn main() {
    let s = SparseEncoder::new(1 << 25, 1 << 25, (1 << 26) as f64);
    // let vector = vec![Complex::new(34.1, 43.4), Complex::new(2123.2, -12.3)];
    let vector = gen_complex_vec(1 << 24);
    // println!("Vector {:?}", vector);
    let encoded = s.encode(&vector);
    // println!("{:?}", encoded);
    let decoded = s.decode(encoded);
    // println!("{:?}", decoded);
    let tolerance = 1e-2;
    
    println!("Approx check: {}", are_approx_equal(&vector, &decoded, tolerance));

}
