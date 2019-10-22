"""
https://gitlab.com/cstein/pepytools

The polarizable embedding python tools gives pythonic access
to polarizable embedding potential files.

Modules:

  - solvers: Solvers for induced charges and moments

  - version: holds the current API version

Classes:

  - Potential: Class to manipulate polarizable embedding potentials

  - IterativeDIISSolver: DIIS solver for induced moments

"""
from .solvers import IterativeDIISSolver

from .potential import Potential
from .intersect import intersect

from .version import __version__

from .tensor import T
from .mulmom import MulMom

__all__ = ['Potential', 'solvers', 'IterativeDIISSolver']
__author__ = "Casper Steinmann"
__copyright__ = "Copyright 2015-2019, Casper Steinmann"
__license__ = "BSD 2-clause"
__email__ = "casper.steinmann@gmail.com"
__status__ = "Development"
__description__ = "Pythonic access to polarizable embedding potentials"
__url__ = 'https://gitlab.com/cstein/pepytools'
