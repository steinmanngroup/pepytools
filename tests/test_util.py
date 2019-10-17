import unittest

import numpy
from numpy.linalg import norm

import util
from pepytools.potential import Potential


class TestUtil(unittest.TestCase):
    """ Testing Util Library """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_center_of_mass(self):
        """ center of mass """
        c = numpy.asarray([[0.0, 0.0, 0.0]])
        m = numpy.asarray([1.0])
        self.assertEqual(norm(util.center_of_mass(c, m)), 0.0)

        c = numpy.asarray([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        m = numpy.asarray([0.5, 0.5])
        self.assertEqual(norm(util.center_of_mass(c, m)), 1.0)

    def test_get_pol_matrix(self):
        """ polarization matrix """
        p = Potential.from_file('nomul.pot')
        pol = p.polarizabilities[0]
        mat = util.get_polarization_matrix(p)
        self.assertEqual(pol[0], mat[0,0])
        self.assertEqual(pol[3], mat[1,1])
        self.assertEqual(pol[5], mat[2,2])
        


def suite():
    s = unittest.TestSuite()
    s.addTest(unittest.makeSuite(TestRegressions))
    return s

if __name__ == '__main__':
    unittest.main()
