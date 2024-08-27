class Ciphertext():

    """A CKKS ciphertext

    Attributes:
    -----------
    ciphertext: The value of the ciphertext
    modulus: Current modulus
    """

    def __init__(self, ciphertext, scale, schemeObject):
        self.ciphertext = ciphertext
        self.modulus = scale
        self.scheme_instance = schemeObject
        self.bootstrapper = Bootstrap(self, self.scheme_instance)

    def mod_down(self, new_modulus):
        """ Reduce remaining modulus after multiplication and addition """
        self.modulus = new_modulus
        
    def __add__(self, other):
        return self.scheme_instance.add_ciphertext(self, other)

    def __mult__(self, other):
        if type(other) == Ciphertext:
            return self.scheme_instance.multiply_ciphertext(self, other)
        if type(other) == Plaintext:
            return self.scheme_instance.multiply_scalar(self, other)

    def bootstrap(self):
        """ Bootstrap the ciphertext """
        self.scheme_instance.bootstrap(ciphertext)
        return 0

    __rmult__ = __mult__
