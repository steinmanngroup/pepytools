# Polarizable Embedding Python Tools (pepytools)

The polarizable embedding python tools (pepytools) is a set of executables and a [Python][] API designed to manipulate [polarizable embedding][] (PE) potentials for use in embedded QM/MM calculations in the [DALTON][] quantum chemistry code.

[Python]: http://www.python.org
[polarizable embedding]: https://gitlab.com/pe-software/pelib-public
[DALTON]: http://daltonprogram.org/

The potentials are typically made from a static charge distribution described by charges, dipoles and quadrupoles placed on atomic sites - for example the atoms of a water molecule - and by a classical, dynamic charge distribution described by (an)isotropic electric dipole-dipole polarizabilities - also placed on atomic sites.

## Example
The most typical action to use pepytools is to add individual potentials.
This is done either through the program supplied with pepytools called `pepy_add`

### shell
```sh
pepy_add p1.pot p2.pot > p_new.pot
```

### API
There is of course also access to the API should you choose to make a program that uses potential files in some form.
The equivalent of the above example would be something like:

```python
from pepytools import Potential

pot1 = Potential.from_file("p1.pot")
pot2 = Potential.from_file("p2.pot")

pot3 = pot1 + pot2

print(pot3)
```
The main class of the library, the `Potential` class, allows you to open potentials from files, add multiple potentials as shown above or even construct potentials from multipole moments, polarizabilities and coordinates and write those potentials back to to a file for later use in your QM/MM calculations.

The polarizable embedding potentials also carry anisotropic polarizabilities at the expansion points and this library also comes with solvers for induced moments.
