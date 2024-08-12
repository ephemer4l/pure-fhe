"""Tests for sparse_encoding_full.py."""
import unittest
import random

from sparse_encoding_lite import SparseEncoder

def check_complex_vector_approx_eq(vec1, vec2, error=0.00001,
                                   error_message="Error: vectors are not approximately equal."):
    """Checks whether two vectors are approximately equal.

    Checks that each entry of two vectors of complex numbers are within the given error.

    Args:
        vec1 (list (complex)): Vector 1.
        vec2 (list (complex)): Vector 2.
        error (float): How close each entry of the two vectors must be to be approximately equal.
        error_message (String): Message to output if not equal.

    Returns:
        A Ciphertext which encrypts the same message under a different key.
    """

    assert(len(vec1) == len(vec2)), 'Length of v1 = %d, Length of v2 = %d' % (len(vec1), len(vec2))
    for i in range(len(vec1)):
        if (abs(vec1[i].real - vec2[i].real) > error) or (abs(vec1[i].imag - vec2[i].imag) > error):
            print("-------- VALUES DO NOT MATCH AT INDEX %d --------" % (i))
            print(str(vec1[i]) + " != " + str(vec2[i]))
            print("v1: " + str(vec1[:10]))
            print("v2: " + str(vec2[:10]))
            raise ValueError(error_message)

def sample_random_complex_vector(length):
    """Samples a random complex vector,

    Samples a vector with elements of the form a + bi where a and b
    are chosen uniformly at random from the set [-2048, 2048).

    Args:
        length (int): Length of vector.

    Returns:
        A list of randomly sampled complex values.
    """
    sample = [0] * length
    for i in range(length):   
        a = random.choice([-1,1]) * random.random() * 2048
        b = random.choice([-1,1]) * random.random() * 2048
        sample[i] = a + b * 1j
    return sample
class TestEncoder(unittest.TestCase):
    def setUp(self):
        self.ciph_modulus = 1 << 40
        self.big_modulus = 1 << 1200
        self.scaling_factor = 1 << 30
        self.degree = 1 << 12
    def run_test_encode_decode(self, vec, n):
        """Checks that encode and decode are inverses.

        Encodes the input vector, decodes the result, and checks that
        they match.

        Args:
            vec (list (complex)): Vector of complex numbers to encode.

        Raises:
            ValueError: An error if test fails.
        """
        self.encoder = SparseEncoder(N=1 << 12, n=1 << n, scale = self.scaling_factor)
        plain = self.encoder.encode(vec)
        value = self.encoder.decode(plain)
        check_complex_vector_approx_eq(vec, value, error=0.1)
    def test_encode_decode_01(self):
        for n in range(1, 13):
            vec = sample_random_complex_vector(1 << (n-1))
            print(vec)
            self.run_test_encode_decode(vec, n)
    if __name__ == '__main__':
        res = unittest.main(verbosity=3, exit=False)
