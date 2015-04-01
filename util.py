import numpy
import numpy.linalg

def get_polarization_matrix(potential, **kwargs):

    Aij = numpy.zeros((3 * potential.npols, 3 * potential.npols))
    for isite in range(potential.nsites):
        itensor = potential.has_alpha[isite]
        if itensor == -1:
            continue

        tensor = numpy.zeros((3, 3))
        tensor[0] = potential.polarizabilities[isite][0:3]
        tensor[1][1:] = potential.polarizabilities[isite][3:5]
        tensor[2][2:] = potential.polarizabilities[isite][5]
        tensor[1][0] = tensor[0][1]
        tensor[2][0] = tensor[0][2]
        tensor[2][1] = tensor[1][2]

        ii = 3 * itensor
        jj = ii + 3

        Aij[ii:jj][0][ii:jj] = tensor[0]
        Aij[ii:jj][1][ii:jj] = tensor[1]
        Aij[ii:jj][2][ii:jj] = tensor[2]

    return Aij


def get_interaction_matrix(potential, **kwargs):
    verbose = kwargs.get('verbose', False)
    TT = numpy.zeros((3 * potential.npols, 3 * potential.npols))
    try:
        from field import interaction_matrix
        ex = numpy.array([potential.exclusion_list[k] for k in range(len(potential.exclusion_list))])
        q = numpy.array([q[0] for q in potential.multipoles[0]])
        d = numpy.array([d for d in potential.multipoles[1]])
        return interaction_matrix(potential.npols, potential.coordinates, potential.has_alpha, ex)
    except:
        if verbose:
            print("INFO: interaction matrix calculated using (slow) python version.")
        for isite in range(potential.nsites):
            itensor = potential.has_alpha[isite]
            is_polarizable_point = (itensor > -1)
            Ri = potential.coordinates[isite]

            if is_polarizable_point:
                iexclusion_list = potential.exclusion_list[isite]

                for jsite in range(potential.nsites):
                    if jsite in iexclusion_list:
                        continue

                    if jsite == isite:
                        continue

                    jtensor = potential.has_alpha[jsite]
                    js_polarizable_point = (jtensor > -1)
                    Rj = potential.coordinates[jsite]

                    if js_polarizable_point:
                        Rij = Rj - Ri
                        R = numpy.sqrt(numpy.dot(Rij, Rij))
                        R1i = 1.0 / R
                        R3i = R1i * R1i * R1i
                        R5i = R3i * R1i * R1i
                        for ii in range(3):
                            for jj in range(3):
                                iii = 3 * itensor + ii
                                jjj = 3 * jtensor + jj
                                TT[iii][jjj] = TT[iii][jjj] + 3.0 * Rij[ii] * Rij[jj] * R5i
                                if ii == jj:
                                    TT[iii][jjj] -= R3i

        return TT
