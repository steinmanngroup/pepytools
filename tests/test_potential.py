import os
from pepytools import Potential

def test_m0_from_file():
    p = Potential.from_file('tests/m0.pot')
    assert p.nsites == 9  # 3 water molecules
    assert p.npols == 0
    print(p)

def test_m2p2_from_file():
    p = Potential.from_file('tests/m2p2.pot')
    assert p.nsites == 6
    assert p.npols == 6

def test_equality():
    p1 = Potential.from_file('tests/m0.pot')
    p2 = Potential.from_file('tests/m0.pot')
    print(p1, p2)
    assert p1 == p2

def test_from_multipoles_quad():
    c1 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    q1 = [[1.0],[2.0]]
    mu1 = [[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]]
    theta1 = [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
              [-1.0, 2.0, -3.0, 4.0, -5.0, 6.0]]
    p1 = Potential.from_multipoles( c1, [q1, mu1, theta1], max_k=2)
    p1.save("quad_test.pot")

    p2 = Potential.from_file("tests/quad_test.pot")
    assert p1 == p2

def test_from_multipoles():
    c1 = [[0.0, 0.0, 0.0]]
    q1 = [[1.0]]
    p1 = Potential.from_multipoles( c1, q1 )

    assert p1.nsites == 1
    assert p1.multipoles[0] == [[1.0]]
    assert p1.npols == 0

    c2 = [[5.4, 0.0, 0.0]]
    d2 = [[1.0, 1.0, 0.0]]
    p2 = Potential.from_multipoles( c2, d2 )
    assert p2.nsites == 1
    assert p2.multipoles[0] == [[0.0]]
    assert p2.multipoles[1] == [[1.0, 1.0, 0.0]]
    assert p2.npols == 0

def test_addition():
    """ Tests various kinds of potential addition """

    # first test addition of empty potentials
    p0 = Potential()
    p = p0 + p0
    assert p.nsites == 2*p0.nsites
    assert p.npols == 2*p0.npols

    # do simple potentials and make sure they add both a + b and b + a
    c1 = [[0.0, 0.0, 0.0]]
    q1 = [[1.0]]
    p1 = Potential.from_multipoles( c1, q1 )

    c2 = [[5.4, 0.0, 0.0]]
    d2 = [[1.0, 1.0, 0.0]]
    p2 = Potential.from_multipoles( c2, d2 )

    p = p1 + p2
    assert p.nsites == p1.nsites + p2.nsites
    assert p.npols == p1.npols + p2.npols

    p = p2 + p1
    assert p.nsites == p1.nsites + p2.nsites
    assert p.npols == p1.npols + p2.npols

    # now add polarizabilities to one
    p1.polarizabilities = [[1.0, 0.0, 0.0, 1.0, 0.0, 1.0]]
    assert p1.npols == 1
    p = p1 + p2
    assert p.nsites == p1.nsites + p2.nsites
    assert p.npols == p1.npols + p2.npols
    p = p2 + p1
    assert p.nsites == p1.nsites + p2.nsites
    assert p.npols == p1.npols + p2.npols

    p2.polarizabilities = [[0.0, 2.0, 2.0, 0.0, 0.0, 2.0]]
    assert p2.npols == 1
    p = p1 + p2
    assert p.nsites == p1.nsites + p2.nsites
    assert p.npols == p1.npols + p2.npols

    # finally add two potentials with only polarization
    p3 = Potential.from_file('tests/nomul.pot')
    print(p3)
    assert p3.npols == 1
    p = p3 + p3
    assert p.npols == 2
