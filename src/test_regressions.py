from potential import Potential
import unittest


class TestRegressions(unittest.TestCase):
    """ Errors that were found during applying pepytools """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_add_potentials_with_multipoles_mixed(self):
        """ Adds potentials with mixed multipole moments """
        m2p2 = Potential.from_file('m2p2.pot')
        m0 = Potential.from_file('m0.pot')

        # check we can add both ways
        p2 = m0 + m2p2
        self.assertEqual(p2.nsites, 15)
        self.assertEqual(len(p2.polarizabilities), 15)
        self.assertEqual(max(p2.multipoles.keys()), 2)

        p1 = m2p2 + m0
        self.assertEqual(p1.nsites, 15)
        self.assertEqual(len(p1.polarizabilities), 15)
        self.assertEqual(max(p1.multipoles.keys()), 2)

def suite():
    s = unittest.TestSuite()
    s.addTest(unittest.makeSuite(TestRegressions))
    return s

if __name__ == '__main__':
    unittest.main()
