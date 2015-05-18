# Polarizable Embedding Python Library (PEPYlib)

The Polarizable Embedding Python Library (PEPYlib) is a [Python][] library designed to manipulate [polarizable embedding][] (PE) potentials for use in embedded QM/MM calculations in the [DALTON][] quantum chemistry code.

[Python]: http://www.python.org
[polarizable embedding]: https://gitlab.com/pe-software/pelib-public
[DALTON]: http://daltonprogram.org/

The main class of the library, the `Potential` class, allows you to open potentials from files, add multiple potentials, construct potentials from multipole moments, polarizabilities and coordinates and write those potentials back to to a file for later use.

The PE potentials have the option of carrying anisotropic polarizabilities at the expansion points and this library also comes with solvers for induced moments.
