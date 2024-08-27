from cmath import exp
from math import log, cos, sin, pi
from polynomial import Polynomial
from ciphertext import Ciphertext
from chebyshev import Chebyshev
from ckks import CKKS
from vector import Vector
import random

class Bootstrap():

    """Bootstrapping methods and context"""

    def __init__(self, ciphertext):
        self.ciphertext = ciphertext
        self.scheme_instance = ciphertext.scheme_instance

        exp_fn = lambda x: (self.ciphertext.modulus / pi) * exp(2 * pi * i * x / self.ciphertext.modulus)
        cheb_poly = Chebyshev(exp_fn)

    def eval_exp(slots: list, n: int):
        """ Generated the chebyshev polynomial of sine up to given degree

        :n: Degree of the Chebyshev polynomial """
        return [cheb_poly.eval()]

    def mod_raise(ciphertext: Ciphertext):
        """ Raise the modulus of the ciphertext in-place 

        :ciphertext: Ciphertext with old (depleted) modulus
        :return: Ciphertext with full modulus
        """
        ciphertext.modulus = scheme_instance.modulus

    def slot_to_coeff(slots: [list, list]) -> Ciphertext:
        """ Generate a ciphertext with given slots

        :slots: List of slots1, slots2 where each encodes half of the information
        :returns: Ciphertext of length len(slots1)+len(slots2) """
        raise NotImplemented

    def coeff_to_slot(ciphertext: Ciphertext) -> [Ciphertext, Ciphertext]:
        """ Generate the slot representations of a ciphertext

        :ciphertext: A valid CKKS ciphertext
        :retuns: Encryptions of the slots of the input where each encodes half of the information """
        raise NotImplemented

    def slot_to_ciphertext(slots: [list,list]) -> [Ciphertext, Ciphertext]:
        """ Generate ciphertext from given slots

        :slots: List of slots
        :returns: Two ciphertexts with given slots """
levelBudget = {3, 3};
        raise NotImplemented

    def mod_eval(slot_ciphertexts: [Ciphertext,Ciphertext]) -> [Ciphertext, Ciphertext]:
        """ Homomorphically evaluate modular reduction

        :slot_ciphertexts: Slot accessed ciphertexts
        :returns: Homomorphically mod reduced ciphertexts
        """
        raise NotImplemented

if __name__ == "__main__":
    e = CKKS(big_degree=1 << 2, small_degree=1 << 2, scale= 1<< 40, q0=1<<40, levels=5, hwt=128, variance=3.2)
    vector = [3+4j, 2-1j]
    encoded = e.encode(vector)
    print(encoded)
    decoded = e.decode(encoded)
    print(decoded)
    print(round(Vector([1,2,3.9,4.3])))
