import math
import pytest
import numpy

from .mulmom import MulMom

class T(list):
    """ Implements an interaction tensor T according Stones

        Elements are stored in canonical order for orders > 0, i.e.
        the order that one would expect: x, y, z for order = 1 or
        xx, xy, xz, yy, yz, zz for order = 2
    """
    def __init__(self, order, Rc, **kwargs):
        """ Creates an interaction tensor

            The interaction tensor is defined for a separation distance Rc
            to a specific order both given as arguments

            Arguments:
            order -- the order of the interaction tensor
            Rc -- the separation (or interaction) distance

            Keyword Arguments:
            verbose -- specifies verbosity. Accepted values are True or False
        """
        self.verbose = kwargs.get('verbose', False)
        if not isinstance(self.verbose, bool):
            raise ValueError("Keyword Argument 'verbose' has wrong type.")

        self.order = order
        tensor_sizes = {0: 1, 1: 3, 2: 6, 3: 10, 4: 15}
        if order < 0:
            raise ValueError("negative order not supported")
        if order > 5:
            raise NotImplementedError("order {} not implemented for interaction tensor".format(order))

        list.__init__(self, [0.0]*tensor_sizes[order])
        x = Rc[0]
        y = Rc[1]
        z = Rc[2]
        r = math.sqrt(x*x + y*y + z*z)

        if r < 1.0e-9:
            raise ValueError("Distance R is close to zero.")

        self._x = x
        self._y = y
        self._z = z

        if order == 0:
            # eq 3.1.4
            self[0] = 1 / r
            self._np = list(self)
        elif order == 1:
            # eq 3.1.5
            R3 = r**3
            self[0] = -x / R3
            self[1] = -y / R3
            self[2] = -z / R3
            self._np = list(self)
        elif order == 2:
            # eq 3.1.6
            R5 = r**5
            self[0] = (3*x*x - r**2) / R5 # xx
            self[1] = (3*x*y       ) / R5 # xy
            self[2] = (3*x*z       ) / R5 # xz
            self[3] = (3*y*y - r**2) / R5 # yy
            self[4] = (3*y*z       ) / R5 # xz
            self[5] = (3*z*z - r**2) / R5 # zz
            self._np = [[self[0], self[1], self[2]],
                        [self[1], self[3], self[4]],
                        [self[2], self[4], self[5]]]
        elif order == 3:
            # eq 3.1.7
            R7 = r**7
            self[0] = -(15*x*x*x - 3*r*r * (x + x + x))/ R7 # xxx
            self[1] = -(15*x*x*y - 3*r*r * y) / R7          # xxy
            self[2] = -(15*x*x*z - 3*r*r * z) / R7          # xxz
            self[3] = -(15*x*y*y - 3*r*r * x) / R7          # xyy
            self[4] = -(15*x*y*z) / R7                      # xyz
            self[5] = -(15*x*z*z - 3*r*r * x) / R7          # xzz
            self[6] = -(15*y*y*y - 3*r*r * (y + y + y))/ R7 # yyy
            self[7] = -(15*y*y*z - 3*r*r * z) / R7          # yyz
            self[8] = -(15*y*z*z - 3*r*r * y) / R7          # yzz
            self[9] = -(15*z*z*z - 3*r*r * (z + z + z))/ R7 # zzz
            self._np = None
        elif order == 4:
            # eq 3.1.8
            R9 = r**9
            R2 = r**2
            R4 = R2 * R2
            self[ 0] = (105 *x*x*x*x - 15* R2 * (x*x + x*x + x*x + x*x + x*x + x*x) + 3* R4 * (1 + 1 + 1)) / R9  # xxxx
            self[ 1] = (105 *x*x*x*y - 15* R2 * (  0 +   0 + x*y +   0 + x*y + x*y) + 3* R4 * (0 + 0 + 0)) / R9  # xxxy
            self[ 2] = (105 *x*x*x*z - 15* R2 * (  0 +   0 + x*z +   0 + x*z + x*z) + 3* R4 * (0 + 0 + 0)) / R9  # xxxz
            self[ 3] = (105 *x*x*y*y - 15* R2 * (x*x +   0 +   0 +   0 +   0 + y*y) + 3* R4 * (1 + 0 + 0)) / R9  # xxyy
            self[ 4] = (105 *x*x*y*z - 15* R2 * (  0 +   0 +   0 +   0 +   0 + y*z) + 3* R4 * (0 + 0 + 0)) / R9  # xxyz
            self[ 5] = (105 *x*x*z*z - 15* R2 * (x*x +   0 +   0 +   0 +   0 + z*z) + 3* R4 * (1 + 0 + 0)) / R9  # xxzz
            self[ 6] = (105 *x*y*y*y - 15* R2 * (x*y + x*y + x*y +   0 +   0 +   0) + 3* R4 * (0 + 0 + 0)) / R9  # xyyy
            self[ 7] = (105 *x*y*y*z - 15* R2 * (  0 +   0 + x*z +   0 +   0 +   0) + 3* R4 * (0 + 0 + 0)) / R9  # xyyz
            self[ 8] = (105 *x*y*z*z - 15* R2 * (x*y +   0 +   0 +   0 +   0 +   0) + 3* R4 * (0 + 0 + 0)) / R9  # xyzz
            self[ 9] = (105 *x*z*z*z - 15* R2 * (x*z + x*z + x*z +   0 +   0 +   0) + 3* R4 * (0 + 0 + 0)) / R9  # xzzz
            self[10] = (105 *y*y*y*y - 15* R2 * (y*y + y*y + y*y + y*y + y*y + y*y) + 3* R4 * (1 + 1 + 1)) / R9  # yyyy
            self[11] = (105 *y*y*y*z - 15* R2 * (  0 +   0 + y*z +   0 + y*z + y*z) + 3* R4 * (0 + 0 + 0)) / R9  # yyyz
            self[12] = (105 *y*y*z*z - 15* R2 * (y*y +   0 +   0 +   0 +   0 + z*z) + 3* R4 * (1 + 0 + 0)) / R9  # yyzz
            self[13] = (105 *y*z*z*z - 15* R2 * (y*z + y*z + y*z +   0 +   0 +   0) + 3* R4 * (0 + 0 + 0)) / R9  # yzzz
            self[14] = (105 *z*z*z*z - 15* R2 * (z*z + z*z + z*z + z*z + z*z + z*z) + 3* R4 * (1 + 1 + 1)) / R9  # zzzz
        else:
            raise NotImplementedError

    def __array__(self):
        if self._np is None:
            raise NotImplementedError
        return numpy.array(self._np)

    def __repr__(self):
        return "T({}, [{}, {}, {}])".format(self.order, self._x, self._y, self._z)

    def __mul__(self, other):
        if self.verbose:
            print("T.mul:", repr(self), "*", repr(other))
        if isinstance(other, MulMom):

            #print "  orders:", self.order, other.order
            # potential of charge
            if self.order == 0 and other.order == 0:
                return MulMom(self[0] * other[0])

            # electric field of charge
            if self.order == 1 and other.order == 0:
                result = [0.0, 0.0, 0.0]
                try:
                    for i in range(3):
                        result[i] = self[i] * other[0]
                except TypeError:
                    print("other:", other[0], " type:", type(other[0]))
                    print("self :", self[i], " type: ", type(self[i]))
                return MulMom(*result)

            # electric field gradient of charge
            if self.order == 2 and other.order == 0:
                result = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                for i in range(6):
                    result[i] = self[i] * other[0]
                return MulMom(*result)

            # T(0, R) * dipole is undefined
            if self.order == 0 and other.order == 1:
                raise ValueError("T(0, R) * dipole is undefined.")

            # potential of dipole
            if self.order == 1 and other.order == 1:
                result = sum([self[i] * other[i] for i in range(3)])
                return MulMom(result)

            # electric field of dipole
            if self.order == 2 and other.order == 1:
                result = [0.0, 0.0, 0.0]
                result[0] = self[0] * other[0] + self[1] * other[1] + self[2] * other[2] # xx * x + xy * y + xz * z
                result[1] = self[1] * other[0] + self[3] * other[1] + self[4] * other[2] # yx * x + yy * y + yz * z
                result[2] = self[2] * other[0] + self[4] * other[1] + self[5] * other[2] # zx * x + zy * y + zz * z
                return MulMom(*result)

            # T(0, R) * quadrupole is undefined
            if self.order == 0 and other.order == 2:
                raise ValueError("T(0, R) * quadrupole is undefined.")

            # T(1, R) * quadropole is undefined
            if self.order == 1 and other.order == 2:
                raise ValueError("T(1, R) * quadrupole is undefined.")

            # potential of quadrupole
            if self.order == 2 and other.order == 2:
                # first add diagonal parts xx, yy, zz
                result = self[0]*other[0] + self[3]*other[3] + self[5]*other[5]
                # then add off-diagonal parts (multiplied by two to include double counting) xy, xz, yz
                result += 2.0*(self[1]*other[1] + self[1]*other[1] + self[4]*other[4])
                return MulMom(result)

            # electric field of quadrupole
            if self.order == 3 and other.order == 2:
                result = [0.0, 0.0, 0.0]
                result[0] = (self[0] * other[0] + self[1] * other[1] + self[2] * other[2]
                          +  self[1] * other[1] + self[3] * other[3] + self[4] * other[4]
                          +  self[2] * other[2] + self[4] * other[4] + self[0] * other[5])
                result[1] = (self[1] * other[0] + self[3] * other[1] + self[4] * other[2]
                          +  self[3] * other[1] + self[6] * other[3] + self[7] * other[4]
                          +  self[4] * other[2] + self[7] * other[4] + self[8] * other[5])
                result[2] = (self[2] * other[0] + self[4] * other[1] + self[5] * other[2]
                          +  self[4] * other[1] + self[7] * other[3] + self[8] * other[4]
                          +  self[5] * other[2] + self[8] * other[4] + self[9] * other[5])
                return MulMom(*result)

            raise NotImplementedError("Your request of T({}, R) * multipole order {} is not implemented.".format(self.order, other.order))

    #        # T(1, R) * dipole is the potential of a dipole
    #        #if self.order == 1 and other.order == 1:
    #        #    return other * self
    #        if self.order == 1 and other.order == 1:
    #            value = 0.0
    #            for i, (v1, v2) in enumerate(zip(self, other)):
    #                value += v1 * v2
    #            return MulMom(value)
    #        raise TypeError("LOL")
        else:
            return TypeError

    #def __rmul__(self, other):
    #    #print "T.rmul:", type(self), type(other), " repr(other)", repr(other)
    #    if isinstance(other, MulMom):
    #        return self*other

    #    if isinstance(other, float):
    #        for i in range(len(self)):
    #            self[i] *= other
    #        return self

    #    raise TypeError("multiplication of T-tensor with {0:s} not implemented.".format(type(other)))


@pytest.fixture
def setup_vectors():
    import random

    R = [0.0, 0.0, 2.0]
    randR = [(random.random()-1) for i in range(3)]
    r = math.sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2])
    return R, randR, r


def test_order0():
    """ tests correct scaling of order 0 """
    R, _, r = setup_vectors()
    tensor = T(0, R)
    assert 1.0/r == pytest.approx(tensor[0])


def test_order1():
    """ tests correct scaling of order 1 """
    R, _, r = setup_vectors()
    tensor = T(1, R)
    assert -1.0/r**2 == pytest.approx(tensor[2])


def test_order2():
    """ tests tracelessness of order 2 interaction tensor

        The condition that an interaction tensor is traceless
        i.e sum_alpha T_{alpha,alpha} = 0 should be True
        for all input vectors.
    """
    _, randR, _ = setup_vectors()
    tensor = T(2, randR)
    assert 0.0 == pytest.approx(tensor[0] + tensor[3] + tensor[5])


def test_order3():
    """ tests tracelessness of order 3 interaction tensor

        The condition that an interaction tensor is traceless
        i.e sum_alpha T_{alpha,alpha,...} = 0 should be True
        for all input vectors.
    """
    _, randR, _ = setup_vectors()
    tensor = T(3, randR)
    # tracelessness of interaction tensor is defined on
    # pp 45 in Stones' Theory of Intermolecular Forces (2nd ed)
    trace_indices = [0, 1, 2, 3, 5, 6, 7, 8, 9]
    assert 0.0 == pytest.approx( sum( [tensor[i] for i in trace_indices] ) )


def test_order4():
    """ tests tracelessness of order 3 interaction tensor

        The condition that an interaction tensor is traceless
        i.e sum_alpha T_{alpha,alpha,...} = 0 should be True
        for all input vectors.
    """
    _, randR, _ = setup_vectors()
    tensor = T(4, randR)
    # tracelessness of interaction tensor is defined on
    # pp 45 in Stones' Theory of Intermolecular Forces (2nd ed)
    trace_indices = list(range(15))  + [3, 5, 12]
    for i in trace_indices:
        print(i, tensor[i])
    assert 0.0 == pytest.approx( sum( [tensor[i] for i in trace_indices] ) )


if __name__ == '__main__':
    #unittest.main()

    #R = [-2, 0.2, 5.0]
    R = [0, 0,  -5.0]
    print("T0:", list(T(0, R)))
    print("T1:", list(T(1, R)))
    print("T2:", list(T(2, R)))
    dd = MulMom(1.0)
    T1 = T(0, R)
    print(T1 * dd)
    #print tt._np
    #q = MulMom(1.0)
    #print T(1, R) * dd
