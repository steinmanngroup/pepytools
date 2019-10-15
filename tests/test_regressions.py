import os

import numpy

from pepytools import Potential

def test_commutative_addition():
    m2p2 = Potential.from_file('tests/m2p2.pot')
    m0 = Potential.from_file('tests/m0.pot')

    # check we can add both ways
    p2 = m0 + m2p2
    assert p2.nsites == 15
    assert len(p2.polarizabilities)== 15
    assert max(p2.multipoles.keys())== 2

    p1 = m2p2 + m0
    assert p1.nsites == 15
    assert len(p1.polarizabilities)== 15
    assert max(p1.multipoles.keys())== 2

def test_addition_nomul_info():
    p1 = Potential()
    p1.setCoordinates(numpy.asarray([[0.0, 0.0, 0.0]]))
    p1.setLabels(['He'])

    p2 = Potential()
    p2.setCoordinates(numpy.asarray([[1.0, 0.0, 0.0]]))
    p2.setLabels(['He'])

    p = p1 + p2
