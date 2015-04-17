import numpy
import numpy.linalg

def center_of_mass( coordinates, masses=None ):
    total_mass = 0.0
    c = numpy.zeros(3)
    for cin,mass in zip(coordinates, masses):
        c += cin * mass
        total_mass += mass

    return c / total_mass

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
    damping = kwargs.get('induced_damping', False)
    damping_factor = kwargs.get('induced_damping_factor', 2.1304)

    from field import interaction_matrix
    ex = numpy.array([potential.exclusion_list[k] for k in range(len(potential.exclusion_list))])
    q = numpy.array([q[0] for q in potential.multipoles[0]])
    d = numpy.array([d for d in potential.multipoles[1]])
    a = []
    for aa in potential.polarizabilities:
        if numpy.sum(numpy.abs(aa)) > 0.0001:
            a.append((aa[0] + aa[3] + aa[5])/3.0)
    return interaction_matrix(potential.coordinates, potential.has_alpha, ex, a, damping, damping_factor)
