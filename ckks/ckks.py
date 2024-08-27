from encoding import SparseEncoder
from probability import sample_discrete_gaussian, sample_zo, sample_hwt, sample_uniform
from vector import Vector
from ciphertext import Ciphertext
class CKKS():

    """ CKKS Scheme """

    def __init__(self, big_degree: int, small_degree: int, scale, q0: int, levels, hwt: int, variance: int):
        self.big_degree = big_degree
        self.small_degree = small_degree
        # Decimal precision
        self.scale = scale
        self.hwt = hwt
        # Integer precision
        self.q0 = q0
        self.variance = variance
        self.levels = levels
        self.modulus = self.scale**self.levels * self.q0
        self.encoder = SparseEncoder(N=self.big_degree, n=self.small_degree, scale=self.scale)

    def KSGen(self, s_prime):
        #TODO
        a_prime = Polynomial(sample_uniform)

    def KeyGen(self):
        self.sk = Vector([1, sample_hwt(self.hwt, N=self.big_degree)])
        a = Polynomial(sample_uniform(self.big_degree))
        e = sample_discrete_gaussian(self.variance)
        b = -1 * a * s + e
        self.pk = Vector([b, a])

    def encrypt(self, plaintext):
        v = sample_zo(0.5)
        errors = [Polynomial(sample_discrete_gaussian(self.variance, self.big_degree)) for i in range(1)]
        return v * self.pk + Vector([plaintext, 0]) + errors

    def decrypt(self, ciphertext):
        return ciphertext @ self.sk

    def encode(self, vector):
        return self.encoder.encode(vector)

    def decode(self, plaintext):
        return self.encoder.decode(plaintext)

    def multiply_scalar(self, ciphertext, scalar):
        return Vector([scalar * ciphertext[0], scalar * ciphertext[1]])

    def rescale(self, ciphertext):
        # Only for decreasing the level by 1
        c_0, c_1 = ciphertext
        c_0 = (1/self.q0) * c_0
        new_c_0 = Polynomial(c_0, c_0.modulus - 1)
        c_1 = (1/self.q0) * c_1
        new_c_1 = Polynomial(c_1, c_1.modulus - 1)
        return Vector([new_c_0, new_c_1])

    def multiply_ciphertext(self, ciphertext1, ciphertext2):
        b_1, a_1 = ciphertext1
        b_2, a_2 = ciphertext2
        d_0, d_1, d_2 = Vector([b_1 * b_2, a_1 * b_2 + a_2 * b_1, a_1 * a_2])
        return Vector([d_0, d_1]) + Vector([round((1/self.P) * d_2 * self.evk)])

    def add_ciphertext(self, ciphertext1, ciphertext2):
        return ciphertext1 + ciphertext2

    def rotate(self, ciphertext, idx, rot_key):
        c_0, c_1 = ciphertext


    def linear_transformation(self, ciphertext, rot_keys, encoder, matrix):
        # Compute two factors of matrix_len (a power of two), both near its square root.
        matrix_len = len(matrix)
        matrix_len_factor1 = int(sqrt(matrix_len))
        if matrix_len != matrix_len_factor1 * matrix_len_factor1:
            matrix_len_factor1 = int(sqrt(2 * matrix_len))
        matrix_len_factor2 = matrix_len // matrix_len_factor1

        # Compute rotations.
        ciph_rots = [0] * matrix_len_factor1
        ciph_rots[0] = ciph
        for i in range(1, matrix_len_factor1):
            ciph_rots[i] = self.rotate(ciph, i, rot_keys[i])

        # Compute sum.
        outer_sum = None
        for j in range(matrix_len_factor2):
            inner_sum = None
            shift = matrix_len_factor1 * j
            for i in range(matrix_len_factor1):
                diagonal = util.matrix_operations.diagonal(matrix, shift + i)
                diagonal = util.matrix_operations.rotate(diagonal, -shift)
                diagonal_plain = encoder.encode(diagonal, self.scaling_factor)
                dot_prod = self.multiply_plain(ciph_rots[i], diagonal_plain)
                if inner_sum:
                    inner_sum = self.add(inner_sum, dot_prod)
                else:
                    inner_sum = dot_prod

            rotated_sum = self.rotate(inner_sum, shift, rot_keys[shift])
            if outer_sum:
                outer_sum = self.add(outer_sum, rotated_sum)
            else:
                outer_sum = rotated_sum

        outer_sum = self.rescale(outer_sum, self.scaling_factor)
        return outer_sum

    def bootstrap(self, ciphertext: Ciphertext):
        ciphertext.bootstrapper.mod_raise(ciphertext)
        slotted = ciphertext.bootstrapper.coeff_to_slot(raised_modulus)
        q_reduced = ciphertext.bootstrapper.img_ext(eval_exp(slotted))
        coeffs = ciphertext.bootstrapper.slot_to_coeff(q_reduced)
        ciphertext.value = coeffs
        ciphertext.modulus = ciphertext.scheme_instance.scale
        return 0
