import numpy

from reader import read_potential
from constants import BOHRTOAA


class Potential(object):
    """ Instantiated using a Polarizable Embedding (PE) potential file.

        pot = Potential.fromFile( "myfile.pot" )
    """

    def __init__(self):
        """ Create a potential without content
        """
        pass

    @classmethod
    def from_file(cls, filename):
        """ Creates a potential from a file
            Arguments:
            filename -- the filename to instantiate the potential from
        """
        a = cls()
        a._filename = filename
        (c, l, m, p, e) = read_potential(filename)
        a.coordinates = c
        a.labels = l
        a.multipoles = m
        a.polarizabilities = p
        a.exclusion_list = e
        return a

    def getCoordinates(self):
        """ Returns a nsites x 3 array with the coordinates of all
            sites in the potential.
        """
        return self._coordinates

    def setCoordinates(self, coordinates):
        """ Sets the coordinates of system and evaluates the number of sites.

            NB: The coordinates does also include coordinates of non-atomic
                sites such as polarizable sites.

            Arguments:
            coordinates -- the coordinates of the potential
        """
        self._coordinates = coordinates
        self._nsites, n = numpy.shape(coordinates)

    coordinates = property(getCoordinates, setCoordinates,
                           doc='Gets or sets the coordinates of the potential.')

    def getLabels(self):
        return self._labels

    def setLabels(self, labels):
        self._labels = labels

    labels = property(getLabels, setLabels,
                      doc='Gets or stes atomic labels of the potential.')

    def getMultipoles(self):
        return self._multipoles

    def setMultipoles(self, multipoles):
        self._multipoles = multipoles

    multipoles = property(getMultipoles, setMultipoles, doc='Gets or sets the multipoles of the potential.')

    def getPolarizabilities(self):
        """ Returns the polarizabilities of the system.
        """
        return self._polarizabilities

    def setPolarizabilities(self, polarizabilities):
        """ Sets the polarizabilities of the system.

            Arguments:
            polarizabilities -- The polarizability tensors to add.

            NOTE: This will evaluate whether a site is polarizable or not based
                  on the size of the tensor. Currently, it checks if the absolute
                  value of the maximum element is larger than zero.
        """
        self._polarizabilities = polarizabilities
        self._hasalpha = [-1 for p in polarizabilities]
        itensor = 0
        for i, tensor in enumerate(polarizabilities):
            if numpy.abs(numpy.max(tensor)) > 0.0:
                self._hasalpha[i] = itensor
                itensor += 1
        self._npols = itensor

    polarizabilities = property(getPolarizabilities, setPolarizabilities, doc='Gets or sets the polarizabilities of the potential.')

    def getExclusionList(self):
        return self._exclusion_list

    def setExclusionList(self, exclusion_list):
        self._exclusion_list = exclusion_list

    exclusion_list = property(getExclusionList, setExclusionList, doc='Gets or sets the exclusion list of the potential.')

    def getNSites(self):
        return self._nsites

    nsites = property(getNSites, doc='Gets the number of classical sites.')

    def getNPols(self):
        return self._npols

    npols = property(getNPols, doc='Gets the number of polarizable points.')

    def getHasAlpha(self):
        return self._hasalpha

    hasalpha = property(getHasAlpha, doc='List relating coordinate indices to polarizable points.')

    def saveToFile(self, filename):
        """ Save the potential to a file.

            Arguments:
            filename -- The filename to save to.
        """
        f = open(filename, 'w')
        f.write(str(self))
        f.close()

    def __str__(self):
        """ Converts the potential to a string readable format
        """
        sc = "@COORDINATES\n{0}\nAA\n".format(self.nsites)
        for label, coord in zip(self.labels, self.coordinates):
            cc = coord * BOHRTOAA
            sc += "{0:2s}{1:14.8f}{2:14.8f}{3:14.8f}\n".format(label, cc[0], cc[1], cc[2])

        sm = "@MULTIPOLES\n"
        for order in self.multipoles.keys():
            sm += "ORDER {0}\n{1}\n".format(order, self.nsites)
            for i, m in enumerate(self.multipoles[order]):
                sm += "{0:3d}".format(i + 1)
                for v in m:
                    sm += "{0:14.8f}".format(v)
                sm += "\n"

        sp = "@POLARIZABILITIES\nORDER 1 1\n{0}\n".format(self.nsites)
        for i, poltensor in enumerate(self.polarizabilities):
            sp += "{0:3d}".format(i + 1)
            for v in poltensor:
                sp += "{0:14.8f}".format(v)
            sp += "\n"

        se = "EXCLISTS\n{0:d} {1:d}\n".format(self.nsites, len(self.exclusion_list[0]) + 1)
        for i in self.exclusion_list.keys():
            se += "{0:>5d}".format(i + 1)
            excl = self.exclusion_list[i] + 1
            for v in excl:
                se += "{0:5d}".format(v)
            se += "\n"

        return sc + sm + sp + se[:-1]

    def __add__(self, other):
        """ Adds two potential files together

            The strategy is to make an empty potential and copy over
            everything from "self" and "other"
        """
        p = Potential()
        c1 = list(self.coordinates)[:]
        c2 = list(other.coordinates)[:]
        c1.extend(c2)
        p.coordinates = numpy.array(c1)

        l1 = self.labels[:]
        l2 = other.labels[:]
        l1.extend(l2)
        p.labels = l1

        # multipoles are stored as a dictionary starting from 0,
        # here we just copy down everything
        m = dict()
        for key in self.multipoles:
            m[key] = self.multipoles[key][:]
            m[key].extend(other.multipoles[key][:])
        p.multipoles = m

        pol1 = self.polarizabilities[:]
        pol2 = other.polarizabilities[:]
        pol1.extend(pol2)
        p.polarizabilities = pol1

        # exclusionlists have to be updated so that the id's in the
        # list reflect correct atoms. Offset items in the "other" by
        # the number of polarizable sites. Items with a "-1" should
        # not be updated
        n1 = self.npols
        n2 = other.npols
        e1 = self.exclusion_list.copy()
        e2 = other.exclusion_list.copy()

        # find out which one has the more items per key
        ne1 = len(e1[0])
        ne2 = len(e2[0])

        # if e2 has the more items, update elements in e1
        # with the proper length
        if (ne2 > ne1):
            for k in e1.keys():
                items = -1 * numpy.ones(ne2, dtype=int)
                items[:ne1] = e1[k]
                e1[k] = items

        # append e2 to e1
        for k in e2.keys():
            items = -1 * numpy.ones(ne2, dtype=int)
            if (ne1 > ne2):
                items = -1 * numpy.ones(ne1, dtype=int)
            items[:ne2] = e2[k]
            selection = numpy.where(items != -1)
            items[selection] += n1
            e1[k + n1] = items

        p.exclusion_list = e1
        return p

    def makeIsotropicPolarizabilities(self):
        pol1 = self.polarizabilities[:]
        pol2 = []
        for pol in pol1:
            value = (pol[0] + pol[3] + pol[5]) / 3.0
            isopol = numpy.zeros(6)
            isopol[0] = value
            isopol[3] = value
            isopol[5] = value
            pol2.append(isopol)
        self.polarizabilities = pol2[:]

    def __and__(self, other):
        """ Intersection can be computed for two potentials
            by finding appropriate overlapping atoms
        """

        satoms = []
        oatoms = []

        slen = self.nsites
        olen = other.nsites
        sc = numpy.array( self.coordinates )
        oc = numpy.array( other.coordinates )

        try:
            from intersect import intersect
            nmax = max(len(sc), len(oc))
            F,n = intersect(nmax, sc, oc)
            satoms = F[:n,0]
            oatoms = F[:n,1]
        except:

            eps = 1.0e-2
            eps2 = eps*eps

            # first pass is linear in time by comparing atoms directly as we loop over them
            # to avoid expensive pairwise calculations. now the lists might not be equal
            # in length so take the shorter one
            nsites = slen
            if slen > olen:
                nsites = olen

            for ic in range(nsites):
                dr = sc[ic] - oc[ic]
                R2 = dr.dot(dr)
                if R2 < eps2:
                    satoms.append(ic)
                    oatoms.append(ic)

            # single atoms could have been removed, so let's just try again with a spray of offsets
            for ic in range(nsites):
                if ic in satoms:
                    continue

                for offset in [-4, -3, -2, -1, 1, 2, 3, 4]:
                    dr = sc[ic] - oc[ic+offset]
                    R2 = dr.dot(dr)
                    if R2 < eps2:
                        satoms.append(ic)
                        oatoms.append(ic+offset)
                        break

            # next loop is pairwise, looking for all the other ones. do not attempt to use any
            # atoms already in the {s,o}atoms lists.
            for ic in range(nsites):
                if ic in satoms:
                    continue

                for jc in range(nsites):
                    if jc in oatoms:
                        continue

                    dr = sc[ic] - oc[jc]
                    R2 = dr.dot(dr)
                    if R2 < eps2:
                        satoms.append(ic)
                        oatoms.append(jc)
                        break

        # transfer stuff
        p = Potential()
        c = []
        l = []
        pol = []
        m = dict()
        e = dict()

        for i, ic in enumerate(satoms):
            c.append( self.coordinates[ic] )
            l.append( self.labels[ic] )
            pol.append( self.polarizabilities[ic] )

            ex_unfixed = self.exclusion_list[ic]
            ex_fixed = self.fix_exclusion_list( ex_unfixed, satoms )
            e[i] = numpy.array(ex_fixed[:])

        for key in self.multipoles:
            m[key] = []
            for ic in satoms:
                m[key].append(self.multipoles[key][ic])

        p.coordinates = numpy.array(c)
        p.labels = numpy.array(l)
        p.multipoles = m
        p.polarizabilities = pol[:]
        p.exclusion_list = e

        return p

    def fix_exclusion_list( self, exlist, ids ):
        new_list = [-1 for i in exlist]
        offset = 0
        for i, value in enumerate(exlist):
            if value == -1:
                break
            try:
                new_list[i+offset] = ids.index( value )
            except ValueError:
                # offset index to correct the list for values that are not part of the system
                offset -= 1
                #print("index {} was not found in atom list.".format(value))
        return new_list

class TransitionPotential(Potential):

    def __init__(self):
        super(TransitionPotential, self).__init__()

    def removePolarizablePoints( self, coordinates, distance ):
        """ The PE library does not remove anything, but it merely sets the polarizabilites
            (and multipoles) to zero.
        """
        d2 = distance*distance

        Cp = list(self.coordinates)
        nCp = range(len(Cp))
        nCp.reverse()

        coordinates_to_remove = []

        for iCp, p in enumerate(Cp):
            for coord in coordinates:
                dr = p - coord
                R2 = dr.dot(dr)
                #print("R = {0:12.6f}, D = {1:12.6f}".format(R2, d2))
                if R2 < d2:
                    coordinates_to_remove.append( iCp )

        # make the list unique
        coordinates_to_remove = [x for x in set(coordinates_to_remove)]
        coordinates_to_remove.sort()
        #coordinates_to_remove.reverse()

        Ct = list(self.coordinates)
        Cl = list(self.labels)
        M0 = list(self._multipoles[0])
        M1 = list(self._multipoles[1])
        P2 = list(self.polarizabilities)
        EX = self.exclusion_list

        for i in coordinates_to_remove:
            #print "removing point", i
            #Ct.pop(i)
            #Cl.pop(i)
            M0[i] = [0.0]
            M1[i] = [0.0, 0.0, 0.0]
            P2[i] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            #EX.pop(i)

        #self.coordinates = numpy.array(Ct)
        #self.labels = Cl
        self._multipoles[0] = numpy.array(M0)
        self._multipoles[1] = numpy.array(M1)
        #self.exclusion_list = EX
        self.polarizabilities = numpy.array(P2)

        self.coordinates_to_remove = coordinates_to_remove[:]
        #print coordinates_to_remove

    def makeTransitionPotentialFromCharges(self, Ci, Qi):
        """ Sets the static part of the potential to zero

            Arguments:
            Ci -- Coordinates of the transition density charges
            Qi -- Transition density charges
        """
        self.makeTransitionPotential()

        Ct = list(self.coordinates)
        Ct.extend(Ci)
        self.coordinates = numpy.array(Ct)

        Cl = list(self.labels)
        Cl.extend(['Z' for c in Ci])
        self.labels = Cl

        # adds new charges
        M0 = list(self._multipoles[0])
        M0.extend([[q] for q in Qi])
        self._multipoles[0] = numpy.array(M0)

        # dipoles
        M1 = list(self._multipoles[1])
        M1.extend([[0.0, 0.0, 0.0] for q in Qi])
        self._multipoles[1] = numpy.array(M1)

        # adds new zero polarizabilities
        # adds new exclusion list
        P2 = list(self.polarizabilities)
        ex = self.exclusion_list
        maxex = max(ex.keys()) + 1
        lenex = len(ex[0])
        for k in range(maxex, maxex + len(Ci)):
            ex[k] = numpy.array([-1 for t in range(lenex)])
            P2.append(numpy.array([0.0 for t in range(6)]))

        self.exclusion_list = ex
        self.polarizabilities = numpy.array(P2)

    def makeTransitionPotentialFromDipole(self, Ci, Mui):
        """ Sets the static part of the potential to zero

            Arguments:
            Ci  -- coordinate of the transition dipole
            Mui -- Transition dipole
        """
        self.makeTransitionPotential()

        Ct = list(self.coordinates)
        Ct.extend(Ci)
        self.coordinates = numpy.array(Ct)

        Cl = list(self.labels)
        Cl.extend(['Z' for c in Ci])
        self.labels = Cl

        # adds new charges
        M0 = list(self._multipoles[0])
        M0.extend([[0.0] for q in [0]])
        self._multipoles[0] = numpy.array(M0)

        # dipoles
        M1 = list(self._multipoles[1])
        M1.extend([Mui])
        self._multipoles[1] = numpy.array(M1)

        # adds new zero polarizabilities
        # adds new exclusion list
        P2 = list(self.polarizabilities)
        ex = self.exclusion_list
        maxex = max(ex.keys()) + 1
        lenex = len(ex[0])
        for k in range(maxex, maxex + len(Ci)):
            ex[k] = numpy.array([-1 for t in range(lenex)])
            P2.append(numpy.array([0.0 for t in range(6)]))

        self.exclusion_list = ex
        self.polarizabilities = numpy.array(P2)

    def makeTransitionPotential(self):
        for key in self._multipoles:
            self.multipoles[key] = numpy.zeros(numpy.shape(self._multipoles[key]))

if __name__ == '__main__':
    import sys
    from solvers import IterativeSolver, IterativeDIISSolver

    p1 = Potential.from_file(sys.argv[1])
    solver = IterativeDIISSolver(p1, threshold=1.0e-5, max_iter=30, verbose=True, diis_start_from_niter=1)

    muind = solver.Solve()
