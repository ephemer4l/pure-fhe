"""
References:
    https://eprint.iacr.org/2018/153.pdf
    https://eprint.iacr.org/2018/1043.pdf
"""

import numpy as np
import random


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
        root = np.exp((-1j * np.pi) / self.n)

        # We can generate roots of unity by using powers of 5 since the order
        # of 5 in the multiplicative group of integers modulo n is 2^(n-2)
        # (giving all elements congruent to 1 mod 4). The rest are elements
        # congruent to 3 mod 4 and they give the conjugates of our roots.
        self.roots_of_unity = [(root) ** (5**j % self.N) for j in range(0, self.n // 2)]
        self.vandermonde = np.vander(self.roots_of_unity, self.n, increasing=True)

    def decode(self, polynomial: np.polynomial.polynomial.Polynomial):
        """Function for the canonical embedding of Z[Y]/(Y^n+1) into C^(n/2) by
        evaluation at roots of unity

        :polynomial: a polynomial viewed as an element of Z[Y]/(Y^n+1)
        :returns: the de-scaled image of the polynomial under the canonical
        embedding
        """

        # View the polynomial in Z[x^(N/n)]/(x^N+1) as an element of Z[y]/(y^n+1)
        empty = (self.N) * [0]
        for i in range(self.n):
            empty[i] = polynomial.coef[(self.N // self.n) * i]
        polynomial = np.polynomial.Polynomial(empty)

        # Evaluate the polynomial at roots of unity
        ev = polynomial(self.roots_of_unity)

        # Scale for precision
        scaled_ev = ev / self.scale
        return scaled_ev

    @staticmethod
    def pi(vector: np.ndarray):
        """The natural isomorphism between H^k and C^(k/2)

        :vector: Vector in C^(k/2) for some k
        :returns: Vector in H^k constructed from the original vector where the
        (-i)-th coordinate is the conjugate of the (i)-th coordinate for all i
        in k/2
        """
        vector_conjugate = [np.conjugate(x) for x in vector[::-1]]
        return np.concatenate([vector, vector_conjugate])

    @staticmethod
    def coordinate_wise_random_rounding(vector: np.ndarray):
        """Standard implementation of coordinate-wise random rounding as
        described in section 2.4.2 of https://eprint.iacr.org/2013/293.pdf

        :vector: A vector in C^k for some k
        :returns: An integer coefficient vector in C^k with the expected value
        of each coordinate equal to the corresponding coordinate of the input
        vector
        """
        r = vector - np.floor(vector)
        f = np.array(
            [np.random.choice([c, c - 1], 1, p=[1 - c, c]) for c in r]
        ).reshape(-1)
        rounded_coordinates = vector - f
        rounded_coordinates = [int(coeff) for coeff in rounded_coordinates]
        return rounded_coordinates

    def proj(self, vector):
        """Take the orthogonal projection of a vector onto the image of the
        canonical embedding. This is the closest vector in the image to the
        original vector.

        :vector: A vector in C^(n)
        :returns: A vector in the image of the canonical embedding which is
        closest to the original vector
        """

        # Basis for the image of the canonical embedding is just the evaluation
        # of the basis of the polynomial ring at roots of unity (including the
        # conjugates). This basis is orthogonal only because N and so n is a
        # power of 2.
        basis = np.vstack(
            (self.vandermonde, np.flip(self.vandermonde.conj(), axis=0))
        ).T

        # Find the coordinates of the projection in terms of basis vectors
        coordinates = np.array(
            [np.real((np.vdot(vector, u) / np.vdot(u, u))) for u in basis]
        )

        # Round coordinates to actually fall into the image of the integer
        # polynomial ring
        rounded_coordinates = __class__.coordinate_wise_random_rounding(coordinates)

        # Calculate the projection vector using the coordinates
        return np.matmul(basis.T, rounded_coordinates)

    def encode(self, vector):
        """Find the polynomial representation (in the ring Z[x^(N/n)]/(x^N+1))
        of a vector in C^(n/2).

        :vector: A vector in C^(n/2)
        :returns: Polynomial representation of vector in Z[x^(N/n)]/(x^N+1)"""

        # Convert possible list to np.array
        vector = np.array(vector)

        # View as an element of H^n instead of C^(n/2) using the natural
        # isomorphism pi
        expansion = __class__.pi(vector)

        # Scale up for better precision
        scaled_expansion = self.scale * expansion

        # Project to orthogonal basis of im(embedding) and snip the conjugate
        # part
        projection = self.proj(scaled_expansion)[: self.n // 2]

        # Get coefficients of the polynomial using CRT
        coefficients = (
            np.real(
                np.dot(self.vandermonde.conj().T, projection)
                + np.dot(self.vandermonde.T, projection.conj())
            )
            / self.n
        )

        # Round numpy's numerical imprecision
        rounded_coeff = np.round(np.real(coefficients)).astype(int)

        # View as an element of the ring Z[x^(N/n)]/(x^N+1) by shifting coefficients
        empty = (self.N) * [0]
        for i in range(self.n):
            empty[self.N // self.n * i] = rounded_coeff[i]

        return np.polynomial.Polynomial(empty)

np.set_printoptions(precision=3, suppress=True, linewidth=np.nan)

s = SparseEncoder(N=512, n=256, scale=32)
p = s.encode(64 * [3+4j,2-1j])

print(p)
print(s.decode(p))
