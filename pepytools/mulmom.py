from __future__ import print_function

""" Multipole moments

"""

class MulMom(list):
    """ A multipole moment """
    def __init__(self, *values):
        """ Creates a multipole moment

            The input values determine the order of the multipole moment
            so for example if you initialize with

            >>> q = MulMom(0.5)

            you are effectively creating a charge whereas a dipole is
            constructed via

            >>> d = MulMom(1.0, 0.0, 0.0)

            Only symmetry-unique elements have to be provided, so for a
            quadrupole you must provide the elements in the canonical
            order, i.e. xx, xy, xz, yy, yz, zz

            >>> O = MulMom(-4.0, -0.1, 0.0, -4.6, 0.0, -5.0)

            Arguments:
            value -- values of the multipole moment
        """
        moment_orders = {1: 0, 3: 1, 6: 2} # orders
        if not len(values) in moment_orders:
            raise ValueError("Input data not understood.")
        self.order = moment_orders[len(values)]
        self.verbose = False

        tensor_sizes = {0: 1, 1: 3, 2: 6}
        list.__init__(self, [0.0]*tensor_sizes[self.order])

        for i, v in enumerate(values):
            self[i] = v

    def __add__(self, other):
        if self.order != other.order:
            raise ValueError("Order must be the same.")

        result = [0.0]*len(self)
        for i, (v1, v2) in enumerate(zip(self, other)):
            result[i] = v1 + v2

        return MulMom(*result)

    def __mul__(self, other):
        if not isinstance(other, MulMom):
            return other * self

        if self.verbose:
            print("M.mul:", repr(self), "*", repr(other))
        if self.order != other.order:
            raise ValueError("Order must be the same.")

        if self.order == 0: # and other.order == 0:
            return MulMom(self[0] * other[0])

        if self.order == 1: # and other.order == 0:
            result = sum([self[i] * other[i] for i in range(3)])
            return MulMom(result)

        if self.order == 2 and other.order == 2:
            result = sum([self[i] * other[i] for i in [0, 1, 2, 3, 4, 5, 1, 2, 4]])
            return MulMom(result)


    #def __rmul__(self, other):
    #    #print "MulMom.rmul:", type(self), type(other), " repr(other)", repr(other)
    #    if isinstance(other, float):
    #        for i in range(len(self)):
    #            self[i] *= other
    #        return self

    #    return NotImplementedError

    #def __str__(self):
    #    s_num = "{0:9.4f} "
    #   s = ""
    #    for v in self:
    #        s += s_num.format(v)
    #    return s[:-1]

    #def __repr__(self):
    #    s_num = "{0:0.8f},"
    #    s = ""
    #    for v in self:
    #        s += s_num.format(v)

    #    return "MulMom({})".format(s[:-1])

    def is_traceless(self):
        if self.order < 2:
            raise NotImplementedError("Tracelessness is not defined for moments with order below 2.")

        trace = self[0] + self[3] + self[5]
        return trace < 1.0e-8

    def make_traceless(self):
        if self.order < 2:
            raise NotImplementedError("Tracelessness is not defined for moments with order below 2.")

        trace = self[0] + self[3] + self[5]
        xxp = 3*self[0] - trace
        yyp = 3*self[3] - trace
        zzp = 3*self[5] - trace
        self[0] = xxp
        self[3] = yyp
        self[5] = zzp
