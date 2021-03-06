import pytest
import numpy

def intersect(sc, oc):
    """ Computes the intersection between two potential coordinates

        :param numpy.ndarray sc: coordinates of self
        :param numpy.ndarray oc: coordinates of other
    """
    try:
        from pepytools.fintersect import fintersect
    except ModuleNotFoundError:
        pytest.xfail("Cannot test equality because of missing fortran module.")
    else:
        nmax = max(len(sc), len(oc))
        return fintersect(nmax, sc, oc)
