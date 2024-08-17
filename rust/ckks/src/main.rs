extern crate num_complex;
extern crate rand;

use num_complex::Complex;
use rand::distributions::{Distribution, WeightedIndex};
use rand::thread_rng;
use rand::Rng;
use std::f64::consts::PI;

#[inline]
fn reverse_bits(n: usize, no_of_bits: usize) -> usize {
    n.reverse_bits() >> (usize::BITS as usize - no_of_bits)
}

#[inline]
fn bit_reverse_vec(values: &mut [Complex<f64>]) {
    let len = values.len();
    let no_of_bits = (len as f64).log2() as usize;
    for i in 0..len {
        let j = reverse_bits(i, no_of_bits);
        if i < j {
            values.swap(i, j);
        }
    }
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
        let r = v.fract().abs();
        let choices = [r, r - 1.0];
        let weights = [1.0 - r, r];

        let dist = WeightedIndex::new(&weights).unwrap();
        let mut rng = thread_rng();
        let f = choices[dist.sample(&mut rng)];

        (v - f).round() as i64
    }

    fn fft(&self, z: &mut [Complex<f64>], roots_of_unity: &[Complex<f64>]) {
        bit_reverse_vec(z);
        let n_half = self.n / 2;
        let log_len = (n_half as f64).log2() as usize;
        for logm in 1..=log_len {
            let m = 1 << logm;
            for i in (0..n_half).step_by(m) {
                for j in 0..m / 2 {
                    let k = modular_pow(5, j, 4 * m) * (n_half / m);
                    let u = z[i + j];
                    let v = z[i + j + m / 2] * roots_of_unity[k];
                    z[i + j] = u + v;
                    z[i + j + m / 2] = u - v;
                }
            }
        }
    }

    fn fft_inv(&self, z: &mut [Complex<f64>], roots_of_unity_inv: &[Complex<f64>]) {
        let n_half = self.n / 2;
        let log_len = (n_half as f64).log2() as usize;
        for logm in (1..=log_len).rev() {
            let m = 1 << logm;
            for i in (0..n_half).step_by(m) {
                for j in 0..m / 2 {
                    let k = modular_pow(5, j, 4 * m) * (n_half / m);
                    let u = z[i + j];
                    let v = z[i + j + m / 2];
                    z[i + j] = u + v;
                    z[i + j + m / 2] = (u - v) * roots_of_unity_inv[k];
                }
            }
        }
        bit_reverse_vec(z);
        let n_recip = 1.0 / (n_half as f64);
        for element in z.iter_mut() {
            *element *= n_recip;
        }
    }
    fn encode(&self, vector: &mut [Complex<f64>]) -> Vec<i64> {
        self.fft_inv(vector, &self.roots_of_unity_inv);
        let mut result = vec![0; self.N];
        let n_half = self.n / 2;
        let factor = self.N / self.n;
        for i in 0..n_half {
            result[factor * i] = Self::random_rounding(vector[i].re * self.scale);
            result[factor * (i + n_half)] = Self::random_rounding(vector[i].im * self.scale);
        }
        result
    }

    fn decode(&self, polynomial_coeff: &[i64]) -> Vec<Complex<f64>> {
        let n_half = self.n / 2;
        let factor = self.N / self.n;
        let mut message = Vec::with_capacity(n_half);
        for i in 0..n_half {
            let real_part = polynomial_coeff[factor * i] as f64 / self.scale;
            let imag_part = polynomial_coeff[factor * (i + n_half)] as f64 / self.scale;
            message.push(Complex::new(real_part, imag_part));
        }

        self.fft(&mut message, &self.roots_of_unity);

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
    let s = SparseEncoder::new(1 << 20, 1 << 20, (1 << 26) as f64);
    // let mut vector = vec![Complex::new(34.1, 43.4), Complex::new(2123.2, -12.3)];
    let mut vector = gen_complex_vec(1 << 19);
    let checkvec = vector.clone();
    // println!("Vector {:?}", vector);
    let encoded = s.encode(&mut vector);
    // println!("{:?}", encoded);
    let decoded = s.decode(&encoded);
    // println!("{:?}", decoded);
    let tolerance = 1e-2;
    
    println!("Approx check: {}", are_approx_equal(&checkvec, &decoded, tolerance));
}
