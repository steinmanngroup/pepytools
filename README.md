# Polarizable Embedding Python Tools (pepytools)

The polarizable embedding python tools (pepytools) is a [Python][] library designed to manipulate [polarizable embedding][] (PE) potentials for use in embedded QM/MM calculations in the [DALTON][] quantum chemistry code.

[Python]: http://www.python.org
[polarizable embedding]: https://gitlab.com/pe-software/pelib-public
[DALTON]: http://daltonprogram.org/

The main class of the library, the `Potential` class, allows you to open potentials from files, add multiple potentials, construct potentials from multipole moments, polarizabilities and coordinates and write those potentials back to to a file for later use.

The polarizable embedding potentials have the option of carrying anisotropic polarizabilities at the expansion points and this library also comes with solvers for induced moments.

## Example
The most typical action one performs would be to add individual potentials.
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
