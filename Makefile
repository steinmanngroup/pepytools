all: field.so intersect.so

field.so: field.f90
	f2py -c --fcompiler=gnu95 --f90flags="-fopenmp" -lgomp -m field field.f90

intersect.so: intersect.f90
	f2py -c --fcompiler=gnu95 -m intersect intersect.f90

clean:
	rm -f field.so
	rm -f intersect.so
	rm -f *.pyc
