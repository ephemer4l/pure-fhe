"""
References:
    https://eprint.iacr.org/2018/153.pdf
    https://eprint.iacr.org/2018/1043.pdf
"""

from math import log, pi
from cmath import exp
import random

from vector import Vector
from polynomial import Polynomial


def reverse_bits(n, no_of_bits):
    """Reverse bits of an integer with bitwise shifting

    :n: Number
    :no_of_bits: Number of bits in the binary representation of n
    """
    result = 0
    for i in range(no_of_bits):
        result <<= 1
        result |= n & 1
        n >>= 1
    return result


def bit_reverse_vec(values):
    """Reverses list by reversing the bits of the indices.

    Reverse indices of the given list.
    For example, reversing the list [0, 1, 2, 3, 4, 5, 6, 7] would become
    [0, 4, 2, 6, 1, 5, 3, 7], since 1 = 0b001 reversed is 0b100 = 4,
    3 = 0b011 reversed is 0b110 = 6.

    Args:
        values (list): List of values to be reversed. Length of list must be a power of two.

    Returns:
        The reversed list based on indices.
    """
    result = [0] * len(values)
    for i in range(len(values)):
        result[i] = values[reverse_bits(i, int(log(len(values), 2)))]
    return result


class SparseEncoder:
    """Sparse encoding for CKKS algorithm using power of two cyclotomic
    polynomials. We consider a subring Z[x^(N/n)]/(x^N+1)=Z[Y]/(Y^n+1) of the
    usual polynomial ring Z[x]/(x^N+1). Encoding and decoding is done as usual
    on the subring and then converted to the 'notation of' the full ring (i.e.
    y is mapped to x^(N/n) and vice versa). The usual decoding algorithm
    applied to encoded plaintexts just give the concatenation

    Attributes
    ----------
    :N: The power of two cyclotomic polynomial degree
    :n: A divisor of n determining the dimension of the plaintext space C^n
    :scale: The scale denoted by delta in the original paper for preserving
    significant bits
    """

    def __init__(self, N: int, n: int, scale: int):
        self.N = N
        self.n = n
        self.scale = scale

        # Precompute roots of unity for FFT
        self.roots_of_unity = [
            exp(2 * pi * 1j * k / (2 * self.n)) for k in range(0, 2 * self.n)
        ]

        # Precompute inverses of roots of unity for FFT
        self.roots_of_unity_inv = [
            exp(-2 * pi * 1j * k / (2 * self.n)) for k in range(0, 2 * self.n)
        ]

    @staticmethod
    def random_rounding(v: float):
        """Standard implementation of coordinate-wise random rounding as
        described in section 2.4.2 of https://eprint.iacr.org/2013/293.pdf

        :v: A float
        :returns: A rounded integer with the expected value equal to input
        vector
        """
        r = v % 1
        f = random.choices([r, r - 1], weights=[1 - r, r])[0]
        rounded_coord = v - f
        return int(rounded_coord)

    @staticmethod
    def fft(n, z, roots_of_unity):
        """Implementation based on https://eprint.iacr.org/2018/1043.pdf
        Algorithm 1

        :n: Twice the ring size
        :z: Complex vector of length n//2
        :roots_of_unity: 2*n roots of unity generated by some primitive 2n-th root
        """
        w = z
        w = bit_reverse_vec(w)
        log_len = int(log(len(z), 2))
        for logm in range(1, log_len + 1, 1):
            m = 1 << logm
            for i in range(0, n // 2, m):
                for j in range(0, m // 2, 1):
                    k = (5**j % (4 * m)) * ((n // 2) // m)
                    U = w[i + j]
                    V = w[i + j + m // 2]
                    V = V * roots_of_unity[k]
                    w[i + j] = U + V
                    w[i + j + m // 2] = U - V
        return w

    @staticmethod
    def fft_inv(n, z, roots_of_unity_inv):
        """Implementation based on https://eprint.iacr.org/2018/1043.pdf
        Algorithm 1

        :n: Twice the ring dimension
        :z: Complex vector
        :roots_of_unity_inv: 2*n many 2n-th roots of unity starting from the
        conjugate of a primitive 2n-th root (depending on the roots of unity
        used for the fft)
        """
        w = z.copy()
        log_len = int(log(len(z), 2))
        for logm in range(log_len, 0, -1):
            m = 1 << logm
            for i in range(0, n // 2, m):
                for j in range(0, m // 2, 1):
                    k = (5**j % (4 * m)) * (n // (2 * m))
                    U = w[i + j]
                    V = w[i + j + m // 2]
                    w[i + j] = U + V
                    w[i + j + m // 2] = (U - V) * roots_of_unity_inv[k]
        w = bit_reverse_vec(w)
        for i in range(n // 2):
            w[i] /= n // 2
        return w

    def encode(self, vector):
        """Find the polynomial representation (in the ring Z[x^(N/n)]/(x^N+1))
        of a vector in C^(n/2).

        :vector: A vector in C^(n/2)
        :returns: Polynomial representation of vector in Z[x^(N/n)]/(x^N+1)"""

        # Apply FFT and fix coefficient order
        bad = __class__.fft_inv(self.n, vector, self.roots_of_unity_inv)
        message = [0] * ((self.n // 2) << 1)
        for i in range(self.n // 2):
            message[i] = __class__.random_rounding(bad[i].real * self.scale)
            message[i + self.n // 2] = __class__.random_rounding(
                bad[i].imag * self.scale
            )

        # View as an element of the ring Z[x^(N/n)]/(x^N+1) by shifting coefficients
        empty = (self.N) * [0]
        for i in range(self.n):
            empty[self.N // self.n * i] = message[i]

        return Polynomial(empty)

    def decode(self, polynomial_coeff: list):
        """Function for the canonical embedding of Z[Y]/(Y^n+1) into C^(n/2) by
        evaluation at roots of unity

        :polynomial_coeff: coefficients of polynomial viewed as an element of Z[Y]/(Y^n+1)
        :returns: the de-scaled image of the polynomial under the canonical
        embedding
        """
        # View the polynomial in Z[x^(N/n)]/(x^N+1) as an element of Z[y]/(y^n+1) and
        # pack coefficients into a complex vector after rounding
        message = []
        for i in range(self.n // 2):
            message.append(
                complex(
                    polynomial_coeff[(self.N // self.n) * i] / self.scale,
                    polynomial_coeff[(self.N // self.n) * (i + (self.n // 2))]
                    / self.scale,
                )
            )

        # Apply FFT
        ev = __class__.fft(self.n, message, self.roots_of_unity)

        return Vector(ev)


if __name__ == "__main__":
    s = SparseEncoder(N=4, n=4, scale=32)
    p = s.encode([38.1 + 83.4j, 9123.2 - 12.3j])
    print(p)
    print(s.decode(p))
