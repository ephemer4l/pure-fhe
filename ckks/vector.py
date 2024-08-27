class Vector(list):

    def __init__(self, components):
        super().__init__(components)
        self.components = components

    def __matmul__(self, other):
        if len(self) != len(other):
            raise ValueError("Vectors are of different size")
        return Vector([self[i]*other[i] for i in range(len(self))])

    def __add__(self, other):
        if len(self) != len(other):
            raise ValueError("Vectors are of different size")
        return Vector([self[i]+other[i] for i in range(len(self))])
    
    def __mul__(self, other):
        return Vector([other * self[i] for i in range(len(self))])

    __rmul__ = __mul__

    def __mod__(self, modulus):
        return Vector([self[i] % modulus for i in range(len(self))])

    def __round__(self):
        return Vector([round(self[i]) for i in range(len(self))])
