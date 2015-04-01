import numpy
import numpy.linalg

from fields import get_static_field


class BaseSolver(object):
    """ Baseclass for solvers.
    """
    def __init__(self, potential, **kwargs):
        self.potential = potential
        self.verbose = kwargs.get('verbose', False)

    def Solve(self):
        raise NotImplementedError

    def Step(self):
        raise NotImplementedError


class IterativeSolver(BaseSolver):
    """ Basic iterative solver using a simple forward power-iteration
        approach to obtain the induced moments from a potential.

        Used together with a potential like:

        pot = Potential.fromFile( 'mypotential.pot' )
        solver = IterativeSolver( pot )
        muind = solver.Solve()

        NB: The static field and other computationally expensive elements
            are evaluated during instantiation.

    """
    def __init__(self, potential, max_iter=30, threshold=1.0e-7, **kwargs):
        """ Create an IterativeSolver from a potential

            Arguments:
            potential -- the polarizable embedding (PE) potential to evaluate

            Keyword arguments:
            max_iter -- maximum number of iterations to use (Default 30)
            threshold -- convergence threshold for subsequent induced moments (Default 1.0e-7)

        """
        super(IterativeSolver, self).__init__(potential, **kwargs)
        self.max_iterations = max_iter
        self.threshold = threshold

        self.static_field = get_static_field(self.potential, **kwargs)
        self.polarization_matrix = get_polarization_matrix(self.potential, **kwargs)
        self.interaction_matrix = get_interaction_matrix(self.potential, **kwargs)

    def Solve(self):
        """ Iteratively solves for the induced moments by evaluating the induced
            moments from the total field and updating the total field from the
            newly obtain induced moments.

            muind = A * ( F_static + F_ind )
            F_ind = T * muind

            where A is the polarization matrix and T is the interaction matrix
        """
        #muind = numpy.dot(self.polarization_matrix, self.static_field)
        muind = self.polarization_matrix.dot( self.static_field )

        for k in range(1, self.max_iterations):
            # obtain the field from induced dipoles and update total field
            self.induced_field = numpy.dot(self.interaction_matrix, muind)
            total_field = self.static_field + self.induced_field

            # Calculate updated induced dipoles based on the field
            (muind_new, muind_diff) = self.Step(muind, total_field, k)

            # absolute error of updated induced dipoles
            abs_err = numpy.sqrt(numpy.dot(muind_diff, muind_diff))

            # extremely numpy-ish way of calculating the total induced dipole
            n, = numpy.shape(muind)
            muind_temp = numpy.reshape(muind, (n / 3, 3))
            muind_tot = numpy.sum(muind_temp, axis=0)
            D = numpy.sqrt(numpy.dot(muind_tot, muind_tot))

            e_pol = -0.5 * numpy.dot(self.static_field, muind)
            if self.verbose:
                print("iter = {0:3d}, energy = {1:12.8f}, err = {3:12.7f}, |D| = {2:9.5f}".format(k, e_pol, D, abs_err))
            muind = muind_new[:]
            if abs_err < self.threshold:
                break

        # unless we really want to save these two matrices, they
        # are destroyed in order to conserve memory.
        self.polarization_matrix = None
        self.interaction_matrix = None
        return muind

    def Step(self, muind_old, total_field, k):
        ''' Generate induced dipoles for the next step (i+1) using regular
            forward power iteration, i.e.

            dmu^(i+1) = A * (F_static + F_ind(mu^i))
        '''
        muind_new = numpy.dot(self.polarization_matrix, total_field)
        muind_diff = muind_new - muind_old
        return muind_new, muind_diff


class IterativeDIISSolver(IterativeSolver):
    """ Iterative DIIS solver using a simple forward power-iteration
        approach with DIIS accelleration to obtain the induced moments
        from a potential.

        Used together with a potential like:

        pot = Potential.fromFile( 'mypotential.pot' )
        solver = IterativeDIISSolver( pot )
        muind = solver.Solve()

        NB: The static field and other computationally expensive elements
            are evaluated during instantiation.

    """
    def __init__(self, potential, max_iter=30, threshold=1.0e-7, **kwargs):
        """ Create an IterativeSolver from a potential

            Arguments:
            potential -- the polarizable embedding (PE) potential to evaluate

            Keyword arguments:
            max_iter -- maximum number of iterations to use (Default 30)
            threshold -- convergence threshold for subsequent induced moments (Default 1.0e-7)

        """
        super(IterativeDIISSolver, self).__init__(potential, max_iter, threshold, **kwargs)

        self.diis_start_from_niter = kwargs.get('diis_start_from_niter', 1)
        self.diis_vector_threshold = kwargs.get('diis_vector_threshold', 1.0e-4)

        self.muind_storage = []
        self.muind_diff_storage = []

        self.started = False

    def Step(self, muind_old, total_field, k):
        ''' Generate a better guess at the induced dipoles using DIIS

            muind_old  :   old induced dipoles
            total_field:   total field from static multipoles and old induced dipoles
            k          :   current step

            we generate induced dipoles for the next step (i+1) using regular
            forward power iteration, i.e.

            mu^(i+1) = A * (F_static + F_ind(mu^i))

            however, a better guess at the new dipoles can be obtained
            by solving a set of linear equation of the error (difference)

            Bx = a

            where elements of the matrix B, B_ij, is

            B_ij = < dmu^(i) | dmu^(j) >,

            where

            dmu^(i) = mu^(i) - mu^(i-1)

            and the last row and column having -1 for elements except for the
            last element having a zero. x is a vector containing coefficients c_i
            and a lambda and a is the solution vector which is equal to zero, except
            for the last item which is -1.

            solving this gives a set of coefficients c_i for which we can write the
            induced dipoles at the next step as a linear combination of all the
            previous ones, i.e.

            u^(i+1) = \sum_i c_(i) * mu^(i)
        '''
        muind_new = numpy.dot(self.polarization_matrix, total_field)
        muind_diff = muind_new - muind_old

        # store the updated induced dipoles and differences
        if k >= self.diis_start_from_niter:
            if not self.started:
                if self.verbose:
                    print("  *** DIIS STARTED ***")
                self.started = True
            self.muind_storage.append(muind_new)
            self.muind_diff_storage.append(muind_diff)

        if len(self.muind_diff_storage) > 1:
            n = len(self.muind_diff_storage)
            rhs = numpy.zeros(n + 1)
            rhs[n] = -1

            # setup system of linear equations
            B = numpy.zeros((n + 1, n + 1))
            B[:][n] = -1
            B = B.transpose()
            B[:][n] = -1
            B[n][n] = 0

            for i, vi in enumerate(self.muind_diff_storage):
                for j, vj in enumerate(self.muind_diff_storage):
                    B[i, j] = numpy.dot(vi, vj)

            # solve the system of linear equations and return only the
            # n elements we need
            c = numpy.linalg.solve(B, rhs)[:n]

            # make a new guess at the induced dipoles
            muind_new = numpy.zeros(numpy.shape(muind_diff))
            for i in range(n):
                muind_new += c[i] * self.muind_storage[i]

            # remove, if needed, items from the DIIS space
            sel = list(numpy.where(numpy.abs(c) < self.diis_vector_threshold)[0])
            sel.reverse()
            for i in sel:
                self.muind_storage.pop(i)
                self.muind_diff_storage.pop(i)

        return muind_new, muind_diff


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
        from ffield import interaction_matrix
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

if __name__ == '__main__':
    import sys
    from potential import Potential
    filename = sys.argv[1]
    potential = Potential.from_file(filename)
    s = IterativeDIISSolver(potential, verbose=True)
    s.Solve()
