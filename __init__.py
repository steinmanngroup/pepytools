"""
https://gitlab.com/cstein/pepylib

The Polarizable Embedding Python Library gives pythonic access
to PE potential files.

Modules:

- potential: Class to manipulate PE potentials

- solvers: Solvers for induced charges and moments

- version: holds the current API version

"""

__all__ = ['potential', 'solvers']
__author__ = "Casper Steinmann"
__copyright__ = "Copyright 2015, Casper Steinmann"
__license__ = "BSD 2-clause"
__maintainer__ = "Casper Steinmann"
__email__ = "casper.steinmann@gmail.com"
__status__ = "Development"

from . version import __version__
import potential
import solvers

