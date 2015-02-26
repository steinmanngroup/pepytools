ffield.so: ffield.f90
	f2py -c --fcompiler=gnu95 --f90flags="-fopenmp" -lgomp -m ffield ffield.f90

clean:
	rm -f ffield.so
	rm -f *.pyc
