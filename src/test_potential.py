from potential import Potential
import unittest


class TestPotentialModule(unittest.TestCase):

    def setUp(self):
        c1 = [[0.0, 0.0, 0.0]]
        q1 = [[1.0]]
        self.p1 = Potential.from_multipoles( c1, q1 )
        c2 = [[5.4, 0.0, 0.0]]
        d2 = [[1.0, 1.0, 0.0]]
        self.p2 = Potential.from_multipoles( c2, d2 )

        # a potential from a file
        self.pfile = Potential.from_file('pehf_iter.pot')

    def tearDown(self):
        pass

    def test_basic_properties(self):
        self.assertEqual(self.p1.nsites, 1)
        self.assertEqual(self.p1.npols, 0)

    def test_add_potentials_no_polarization(self):
        p3 = self.p1 + self.p2
        self.assertEqual(p3.nsites, 2)
        self.assertEqual(p3.npols, 0)

    def test_add_potentials_with_polarization(self):
        alpha1 = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        alpha2 = [[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]]
        self.p1.polarizabilities = alpha1
        self.p2.polarizabilities = alpha2

        self.assertEqual(self.p1.npols, 1)
        p3 = self.p1 + self.p2
        self.assertEqual(p3.npols, 2)

    def test_add_potentials_with_polarization_mixed(self):
        """ Adds a potential with polarization to a potential
            without
        """
        alpha1 = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        self.p1.polarizabilities = alpha1

        p3 = self.p1 + self.p2
        self.assertEqual(p3.npols, 1)
        self.assertEqual(p3.nsites, 2)

    def test_basic_potential_from_file(self):
        self.assertEqual(self.pfile.nsites, 6)
        self.assertEqual(self.pfile.npols, 6)


def suite():
    s = unittest.TestSuite()
    s.addTest(unittest.makeSuite(TestPotentialModule))
    return s

if __name__ == '__main__':
    unittest.main()
