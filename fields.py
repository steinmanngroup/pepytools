import numpy


def get_static_field(potential, **kwargs):
    """
    """

    verbose = kwargs.get('verbose', False)

    F_static = numpy.zeros(3 * potential.npols)
    try:
        from ffield import fstatic_field
        ex = numpy.array([potential.exclusion_list[k] for k in range(len(potential.exclusion_list))])
        q = numpy.array([q[0] for q in potential.multipoles[0]])
        d = numpy.array([d for d in potential.multipoles[1]])
        F_static = fstatic_field(potential.npols, potential.coordinates, potential.hasalpha, ex, q, d)
    except ImportError:
        if verbose:
            print("INFO: static field calculated using (slow) python version.")
        offset = 0
        for isite in range(potential.nsites):
            itensor = potential.hasalpha[isite]
            is_polarizable_point = (itensor > -1)
            Ri = potential.coordinates[isite]

            if is_polarizable_point:
                iexclusion_list = potential.exclusion_list[isite]

                for jsite in range(potential.nsites):
                    if jsite == isite:
                        continue

                    if jsite in iexclusion_list:
                        continue

                    jtensor = potential.hasalpha[jsite]
                    js_polarizable_point = (jtensor > -1)
                    Rj = potential.coordinates[jsite]
                    Qj = potential.multipoles[0][jsite]
                    Pj = potential.multipoles[1][jsite]

                    dRij = Rj - Ri
                    Rij = numpy.sqrt(numpy.dot(dRij, dRij))
                    R1i = 1.0 / Rij
                    R2i = R1i * R1i
                    R3i = R2i * R1i
                    R5i = R3i * R2i
                    uRij = dRij * R1i
                    M0 = Qj * R2i * uRij
                    M1 = Pj * R3i
                    M1 -= 3.0 * dRij * R5i * numpy.dot(dRij, Pj)
                    F_static[offset:offset + 3] -= M0 + M1
                offset += 3

    return F_static
