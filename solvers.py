import numpy


class BaseSolver(object):
    """ Baseclass for solvers.
    """
    def __init__(self, **kwargs):
        self.verbose = kwargs.get('verbose', False)

    def Solve(self):
        raise NotImplementedError

    def Step(self):
        raise NotImplementedError


class IterativeSolver(BaseSolver):
    """ Basic iterative solver using a simple forward power-iteration approach

        Used together with a potential like:

        pot = Potential.from_file( 'mypotential.pot' )
        solver = IterativeSolver( pot )
        muind = solver.Solve()

        NB: The static field and other computationally expensive elements
            are evaluated during instantiation.

    """
    def __init__(self, polmat, intmat, field, **kwargs):
        super(IterativeSolver, self).__init__(**kwargs)
        self.max_iterations = kwargs.get('max_iterations', 50)
        self.threshold = kwargs.get('threshold', 1.0e-9)

        self.polarization_matrix = polmat
        self.interaction_matrix = intmat
        self.F = field

    def Solve(self):
        """ Iteratively solves for the quantity S which can be obtained by
            the following functional form

            >> S = A * (F_s + F_n)
            >> F_n = T * S

            A is known as the polarization matrix and T is the interaction
            matrix. F is the field or potential that through the polarization
            matrix generates the property of interest, S.
        """
        if self.verbose:
            print("{0:>6s}{1:>16s}{2:>16s}{3:>16s}".format("iter", "energy", "rms error", "max error"))
        guess = self.polarization_matrix.dot(self.F)

        for k in range(1, self.max_iterations):
            F = self.F + self.interaction_matrix.dot(guess)

            (guess, difference) = self.Step(guess, F, k)

            # error is RMS in the difference
            rms_err = numpy.std(difference)
            max_err = numpy.max(difference)

            energy = -0.5 * self.F.dot(guess)
            if self.verbose:
                print("{0:6d}{1:16.8f}{2:16.8f}{3:16.8f}".format(k, energy, rms_err, max_err))

            # according to:
            #   Scalmani G et al. Theo. Chem. Account. 2004, (111), 90 - 100, DOI: 10.1007/s00214-003-0527-2
            #
            # we use a convergence threshold on both the RMS and MAX of the residual
            if rms_err < self.threshold and max_err < self.threshold:
                break

        # destroy matrices to conserve memory
        self.polarization_matrix = None
        self.interaction_matrix = None
        return guess

    def Step(self, old_guess, F, k):
        new_guess = self.polarization_matrix.dot(F)
        difference = new_guess - old_guess
        return new_guess, difference


class IterativeDIISSolver(IterativeSolver):
    """ Iterative solver with DIIS accelleration

        Used together with a potential like:

        pot = Potential.from_file( 'mypotential.pot' )
        solver = IterativeDIISSolver( pot )
        muind = solver.Solve()
    """
    def __init__(self, polmat, intmat, field, **kwargs):
        super(IterativeDIISSolver, self).__init__(polmat, intmat, field, **kwargs)
        self.diis_start_from_iter = kwargs.get('diis_start_from_iter', 1)
        self.diis_vector_threshold = kwargs.get('diis_vector_threshold', 1.0e-4)
        self.diis_started = False

        self.guess_storage = []
        self.diff_storage = []

    def Step(self, old_guess, F, k):
        new_guess = self.polarization_matrix.dot(F)
        difference = new_guess - old_guess

        if k >= self.diis_start_from_iter:

            new_guess = self.diis_step(new_guess, difference)

        return new_guess, difference

    def diis_step(self, guess, difference):
        self.guess_storage.append(guess)
        self.diff_storage.append(difference)

        new_guess = guess[:]

        n = len(self.guess_storage)
        if n > 1:
            if not self.diis_started:
                if self.verbose:
                    print("{0:^60s}".format("--------- DIIS STARTED ---------"))
                self.diis_started = True

            # right hand side of system of equations first
            rhs = numpy.zeros(n + 1)
            rhs[n] = -1

            # prepare setup of system of linear equations
            B = numpy.zeros((n + 1, n + 1))
            B[:][n] = -1
            B = B.transpose()
            B[:][n] = -1
            B[n][n] = 0

            # iterate over storage vectors and take inner product of
            # differences
            for i, idiff in enumerate(self.diff_storage):
                for j, jdiff in enumerate(self.diff_storage):
                    B[i, j] = idiff.dot(jdiff)

            # solve the system of linear equations and only
            # extract the n solutions we want
            coeff = numpy.linalg.solve(B, rhs)[:n]

            # now make a better guess at the solution by taking a
            # linear combination of all the stored guesses so far
            new_guess = numpy.zeros(numpy.shape(self.guess_storage[0]))
            for i in range(n):
                new_guess += coeff[i] * self.guess_storage[i]

            # remove unwanted (non-contributing) solutions from the
            # DIIS space
            selection = list(numpy.where(numpy.abs(coeff) < self.diis_vector_threshold)[0])
            selection.reverse()
            for i in selection:
                self.guess_storage.pop(i)
                self.diff_storage.pop(i)

        return new_guess

if __name__ == '__main__':
    import sys

    from fields import get_static_field
    from potential import Potential
    from util import get_interaction_matrix, get_polarization_matrix

    filename = sys.argv[1]
    p = Potential.from_file(filename)

    polmat = get_polarization_matrix(p)
    intmat = get_interaction_matrix(p)
    field = get_static_field(p)
    s = IterativeDIISSolver(polmat, intmat, field, **{'verbose': True, 'threshold': 1.0e-5})
    result = s.Solve()
