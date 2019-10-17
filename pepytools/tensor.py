import math
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
                mux = other[0]
                muy = other[1]
                muz = other[2]
                #print("oo", list(other))
                #print("mu", mux, muy, muz)
                result[0] = self[0] * mux + self[1] * muy + self[2] * muz # xx * x + xy * y + xz * z
                result[1] = self[1] * mux + self[3] * muy + self[4] * muz # yx * x + yy * y + yz * z
                result[2] = self[2] * mux + self[4] * muy + self[5] * muz # zx * x + zy * y + zz * z
                return MulMom(*result)

            # electric field gradient of dipole
            if self.order == 3 and other.order == 1:
                result = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                mux = other[0]
                muy = other[1]
                muz = other[2]
                result[0] = self[0] * mux + self[1] * muy + self[2] * muz
                result[1] = self[1] * mux + self[3] * muy + self[4] * muz
                result[2] = self[2] * mux + self[4] * muy + self[5] * muz
                result[3] = self[3] * mux + self[6] * muy + self[7] * muz
                result[4] = self[4] * mux + self[7] * muy + self[8] * muz
                result[5] = self[5] * mux + self[8] * muy + self[9] * muz
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
                Oxx = other[0]
                Oxy = other[1]
                Oxz = other[2]
                Oyy = other[3]
                Oyz = other[4]
                Ozz = other[5]
                result[0] = self[0] * Oxx + self[1] * Oxy + self[2] * Oxz + self[ 3] * Oyy + self[ 4] * Oyz + self[ 5] * Ozz
                result[1] = self[1] * Oxx + self[3] * Oxy + self[4] * Oxz + self[ 6] * Oyy + self[ 7] * Oyz + self[ 8] * Ozz
                result[2] = self[2] * Oxx + self[4] * Oxy + self[5] * Oxz + self[ 7] * Oyy + self[ 8] * Oyz + self[ 9] * Ozz
                return MulMom(*result)

            # electric field gradient of quadrupole
            if self.order == 4 and other.order == 2:
                result = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                Oxx = other[0]
                Oxy = other[1]
                Oxz = other[2]
                Oyy = other[3]
                Oyz = other[4]
                Ozz = other[5]
                result[0] = self[0] * Oxx + self[1] * Oxy + self[2] * Oxz + self[ 3] * Oyy + self[ 4] * Oyz + self[ 5] * Ozz
                result[1] = self[1] * Oxx + self[3] * Oxy + self[4] * Oxz + self[ 6] * Oyy + self[ 7] * Oyz + self[ 8] * Ozz
                result[2] = self[2] * Oxx + self[4] * Oxy + self[5] * Oxz + self[ 7] * Oyy + self[ 8] * Oyz + self[ 9] * Ozz
                result[3] = self[3] * Oxx + self[6] * Oxy + self[7] * Oxz + self[10] * Oyy + self[11] * Oyz + self[12] * Ozz
                result[4] = self[4] * Oxx + self[7] * Oxy + self[8] * Oxz + self[11] * Oyy + self[12] * Oyz + self[13] * Ozz
                result[5] = self[5] * Oxx + self[8] * Oxy + self[9] * Oxz + self[12] * Oyy + self[13] * Oyz + self[14] * Ozz
                return MulMom(*result)

            raise NotImplementedError("Your request of T({}, R) * multipole order {} is not implemented.".format(self.order, other.order))
        else:
            return TypeError("Expected {}".format(type(other)))
