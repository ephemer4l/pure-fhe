class Polynomial(list):

    def __init__(self, coeffs, modulus = 0):
        self.modulus = modulus
        if self.modulus > 0:
            coeffs = list(map(lambda x: x % self.modulus, coeffs))
        super().__init__(coeffs)
        self.coeffs = coeffs
        self.ring_degree = len(coeffs)

    def rotate(self, r):
        """Rotates plaintext coefficients by r.

        Rotates all the plaintext coefficients to the left such that the x^r
        coefficient is now the coefficient for x^0. We do so by applying the
        transformation m(X) -> m(X^k), where k = 5^r in the ciphertext
        polynomial.

        Returns:
            A rotated Polynomial.
        """
        k = 5 ** r
        new_coeffs = [0] * self.ring_degree
        for i in range(self.ring_degree):
            index = (i * k) % (2 * self.ring_degree)
            if index < self.ring_degree:
                new_coeffs[index] = self.coeffs[i]
            else:
                new_coeffs[index - self.ring_degree] = -self.coeffs[i]
        return Polynomial(self.ring_degree, new_coeffs)

    def evaluate(self, inp):
        """Evaluates the polynomial at the given input value.

        Evaluates the polynomial using Horner's method.

        Args:
            inp (int): Value to evaluate polynomial at.

        Returns:
            Evaluation of polynomial at input.
        """
        result = self.coeffs[-1]

        for i in range(self.ring_degree - 2, -1, -1):
            result = result * inp + self.coeffs[i]

        return result


    def __add__(self, poly, coeff_modulus=self.modulus):
        """Adds two polynomials in the ring.

        Adds the current polynomial to poly inside the ring R_a.

        Args:
            poly (Polynomial): Polynomial to be added to the current
                polynomial.
            coeff_modulus (int): Modulus a of coefficients of polynomial
                ring R_a.

        Returns:
            A Polynomial which is the sum of the two polynomials.
        """
        assert isinstance(poly, Polynomial)

        poly_sum = Polynomial(self.ring_degree, [0] * self.ring_degree)

        poly_sum.coeffs = [self.coeffs[i] + poly.coeffs[i] for i in range(self.ring_degree)]
        if coeff_modulus:
            poly_sum = poly_sum.mod(coeff_modulus)
        return poly_sum

    def __mod__(self, modulus):
        return Polynomial([coeff % modulus for coeff in self.coeffs])

    def __mul__(self, scalar, coeff_modulus=self.modulus):
        return Polynomial([scalar * self[i] for i in self.ring_degree])
        
    __rmul__ = __mul__

    def __matmul__(self, poly, coeff_modulus=self.modulus):
        """Multiplies two polynomials in the ring in O(n^2).

        Multiplies the current polynomial to poly inside the ring R_a
        naively in O(n^2) time.

        Args:
            poly (Polynomial): Polynomial to be multiplied to the current
                polynomial.
            coeff_modulus (int): Modulus a of coefficients of polynomial
                ring R_a.

        Returns:
            A Polynomial which is the product of the two polynomials.
        """
        assert isinstance(poly, Polynomial)

        poly_prod = Polynomial(self.ring_degree,
                               [0] * self.ring_degree)

        for d in range(2 * self.ring_degree - 1):
            # Since x^d = -1, the degree is taken mod d, and the sign
            # changes when the exponent is > d.
            index = d % self.ring_degree
            sign = int(d < self.ring_degree) * 2 - 1

            # Perform a convolution to compute the coefficient for x^d.
            coeff = 0
            for i in range(self.ring_degree):
                if 0 <= d - i < self.ring_degree:
                    coeff += self.coeffs[i] * poly.coeffs[d - i]
            poly_prod.coeffs[index] += sign * coeff
            
            if coeff_modulus:
                poly_prod.coeffs[index] %= coeff_modulus

        return poly_prod

    def __str__(self):
        """Represents polynomial as a readable string.

        Returns:
            A string which represents the Polynomial.
        """
        s = ''
        for i in range(self.ring_degree - 1, -1, -1):
            if self.coeffs[i] != 0:
                if s != '':
                    s += ' + '
                if i == 0 or self.coeffs[i] != 1:
                    s += str(int(self.coeffs[i]))
                if i != 0:
                    s += 'x'
                if i > 1:
                    s += '^' + str(i)
        return s
