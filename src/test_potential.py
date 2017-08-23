import os
import unittest

from potential import Potential

def delete_file(filename):
    try:
        f = open(filename)
    except IOError:
        return
    finally:
        f.close()
        os.remove(filename)

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
        """ Potential basic properties """
        self.assertEqual(self.p1.nsites, 1)
        self.assertEqual(self.p1.npols, 0)

    def test_add_potentials_no_polarization(self):
        """ Potential addition without polarization """
        p3 = self.p1 + self.p2
        self.assertEqual(p3.nsites, 2)
        self.assertEqual(p3.npols, 0)

    def test_add_potentials_with_polarization(self):
        """ Potential addition with polarization """
        alpha1 = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        alpha2 = [[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]]
        self.p1.polarizabilities = alpha1
        self.p2.polarizabilities = alpha2

        self.assertEqual(self.p1.npols, 1)
        p3 = self.p1 + self.p2
        self.assertEqual(p3.npols, 2)

    def test_add_potentials_with_polarization_mixed(self):
        """ Adds a potential with polarization to a potential without
        """
        alpha1 = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        self.p1.polarizabilities = alpha1

        p3 = self.p1 + self.p2
        self.assertEqual(p3.npols, 1)
        self.assertEqual(p3.nsites, 2)

        p4 = self.p2 + self.p1
        self.assertEqual(p4.npols, 1)
        self.assertEqual(p4.nsites, 2)

    def test_add_empty_potentials(self):
        """ Potential addition with empty potential """
        # first we test construction of empty potential
        p = Potential()
        self.assertEqual(p.labels, [])
        self.assertEqual(p.nsites, 0)
        self.assertEqual(p.npols, 0)

        p2 = p + p # maybe we should just raise a value error instead?
        self.assertEqual(p2.nsites, 0)
        self.assertEqual(p2.npols, 0)

    def test_add_potentials_with_only_polarization(self):
        """ Potential addition with only polarization """
        p = Potential.from_file('nomul.pot')
        p2 = p + p
        self.assertEqual(p2.nsites, 2)
        self.assertEqual(p2.npols, 2)

    def test_basic_potential_from_file(self):
        """ Load potential from file """
        self.assertEqual(self.pfile.nsites, 6)
        self.assertEqual(self.pfile.npols, 6)

    def test_potential_charge(self):
        """ tests the total charge """
        p = Potential.from_file('m0.pot')
        self.assertEqual(p.charge, 0)

        p = Potential.from_file('nomul.pot')
        self.assertEqual(p.charge, 0)

    def test_potential_equal(self):
        """ tests equality of potentials """
        p = Potential.from_file('m0.pot')
        self.assertEqual(p == p, True)

    def test_printed_potential_match(self):
        p = Potential.from_file('m2p2.pot')
        with open('printed_potential.pot', 'w') as f:
            f.write(str(p))
        p2 = Potential.from_file('printed_potential.pot')
        self.assertEqual(p == p2, True)
        delete_file('printed_potential.pot')
        
        

def suite():
    s = unittest.TestSuite()
    s.addTest(unittest.makeSuite(TestPotentialModule))
    return s

if __name__ == '__main__':
    unittest.main()
