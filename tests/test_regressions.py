import os

from pepytools import Potential

def test_commutative_addition():
    m2p2 = Potential.from_file('m2p2.pot')
    m0 = Potential.from_file('m0.pot')

    # check we can add both ways
    p2 = m0 + m2p2
    assert p2.nsites == 15
    assert len(p2.polarizabilities)== 15
    assert max(p2.multipoles.keys())== 2

    p1 = m2p2 + m0
    assert p1.nsites == 15
    assert len(p1.polarizabilities)== 15
    assert max(p1.multipoles.keys())== 2
