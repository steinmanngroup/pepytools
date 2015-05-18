all: field.so intersect.so qm_fields.so

field.so: field.f90
	f2py -c --fcompiler=gnu95 --f90flags="-fopenmp" -lgomp -m field $^

intersect.so: intersect.f90
	f2py -c --fcompiler=gnu95 -m intersect $^

qm_fields.so: qm_fields.f90
	f2py -c --fcompiler=gnu95 -m qm_fields $^

clean:
	rm -f field.so
	rm -f intersect.so
	rm -f *.pyc
