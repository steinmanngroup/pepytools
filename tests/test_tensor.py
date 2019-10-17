import math
import numpy
import pytest

from pepytools.tensor import T
from pepytools.mulmom import MulMom

@pytest.fixture
def setup_vectors():
    import random

    R = [0.0, 0.0, 2.0]
    randR = [(random.random()-1) for i in range(3)]
    r = math.sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2])
    return R, randR, r


def test_order0(setup_vectors):
    """ tests correct scaling of order 0 """
    R = setup_vectors[0]
    r = setup_vectors[2]
    tensor = T(0, R)
    assert 1.0/r == pytest.approx(tensor[0])


def test_order1(setup_vectors):
    """ tests correct scaling of order 1 """
    R = setup_vectors[0]
    r = setup_vectors[2]
    tensor = T(1, R)
    assert -1.0/r**2 == pytest.approx(tensor[2])


def test_order2(setup_vectors):
    """ tests tracelessness of order 2 interaction tensor

        The condition that an interaction tensor is traceless
        i.e sum_alpha T_{alpha,alpha} = 0 should be True
        for all input vectors.
    """
    randR = setup_vectors[1]
    tensor = T(2, randR)
    assert 0.0 == pytest.approx(tensor[0] + tensor[3] + tensor[5])


def test_order3(setup_vectors):
    """ tests tracelessness of order 3 interaction tensor

        The condition that an interaction tensor is traceless
        i.e sum_alpha T_{alpha,alpha,...} = 0 should be True
        for all input vectors.
    """
    randR = setup_vectors[1]
    tensor = T(3, randR)
    # tracelessness of interaction tensor is defined on
    # pp 45 in Stones' Theory of Intermolecular Forces (2nd ed)
    trace_indices = [0, 1, 2, 3, 5, 6, 7, 8, 9]
    assert 0.0 == pytest.approx( sum( [tensor[i] for i in trace_indices] ) )


def test_order4(setup_vectors):
    """ tests tracelessness of order 3 interaction tensor

        The condition that an interaction tensor is traceless
        i.e sum_alpha T_{alpha,alpha,...} = 0 should be True
        for all input vectors.
    """
    randR = setup_vectors[1]
    tensor = T(4, randR)
    # tracelessness of interaction tensor is defined on
    # pp 45 in Stones' Theory of Intermolecular Forces (2nd ed)
    trace_indices = list(range(15))  + [3, 5, 12]
    for i in trace_indices:
        print(i, tensor[i])
    assert 0.0 == pytest.approx( sum( [tensor[i] for i in trace_indices] ) )
