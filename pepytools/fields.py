""" obtain the static field from a set of charges and
    dipoles at polarizable points.
"""
import numpy
from .tensor import T
from .mulmom import MulMom as M

def get_static_field_from_file(potential, filename):
    f = open(filename, 'r')
    field = []
    for i, line in enumerate(f):
        d = line.split()
        if i == 0:
            continue
        else:
            field.append(list(map(float, d)))
    f.close()

    return numpy.ravel(field)


def get_static_field(potential, **kwargs):
    """
    """

    verbose = kwargs.get('verbose', False)
    filename = kwargs.pop('filename', None)

    F_static = numpy.zeros(3 * potential.npols)

    if filename is not None:
        if verbose:
            print("Loading static field from file '{0}'".format(filename))

        F_static = get_static_field_from_file(potential, filename)
    else:

        try:
            from .field import static_field
            ex = numpy.array([potential.exclusion_list[k] for k in range(len(potential.exclusion_list))])
            q = numpy.zeros(potential.nsites)
            d = numpy.zeros((potential.nsites,3))

            multipoles = potential.multipoles
            if 0 in multipoles.keys():
                q = numpy.array([q[0] for q in multipoles[0]])
            if 1 in multipoles.keys():
                d = numpy.array([d for d in multipoles[1]])

            F_static = static_field(potential.npols, potential.coordinates, potential.has_alpha, ex, q, d)
        except ModuleNotFoundError:
            if verbose:
                print("INFO: static field calculated using (slow) python version.")
            offset = 0
            for isite in range(potential.nsites):
                itensor = potential.has_alpha[isite]
                is_polarizable_point = (itensor > -1)
                Ri = potential.coordinates[isite]

                if is_polarizable_point:
                    iexclusion_list = potential.exclusion_list[isite]

                    for jsite in range(potential.nsites):
                        if jsite == isite:
                            continue

                        if jsite in iexclusion_list:
                            continue

                        jtensor = potential.has_alpha[jsite]
                        js_polarizable_point = (jtensor > -1)
                        Rj = potential.coordinates[jsite]
                        dRij = Rj - Ri

                        T1 = T(1, dRij)
                        try:
                            M0 = M(*potential.multipoles[0][jsite])
                        except KeyError:
                            F0 = numpy.zeros(3)
                        else:
                            F0 = numpy.array(M0 * T1).ravel()
                        finally:
                            F_static[offset:offset + 3] -= F0

                        T2 = T(2, dRij)
                        try:
                            M1 = M(*potential.multipoles[1][jsite])
                        except KeyError:
                            F1 = numpy.zeros(3)
                        else:
                            F1 = numpy.array(M1 * T2)
                        finally:
                            F_static[offset:offset + 3] += F1

                        T3 = T(3, dRij)
                        try:
                            M2 = M(*potential.multipoles[2][jsite])
                        except KeyError:
                            F2 = numpy.zeros(3)
                        else:
                            F2 = numpy.array(M2 * T3)
                        finally:
                            F_static[offset:offset + 3] += F2
                    offset += 3

    return F_static
