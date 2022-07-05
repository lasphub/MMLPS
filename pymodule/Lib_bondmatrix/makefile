f = arrayprocess.f90 structure.f90 bondmatrix.f90 

FC      = ifort
FCFLAGS = -shared -fpic -O2 -xHost #-openmp 

bondmatrix.so: ${f}
	${FC} ${FCFLAGS}  -o bondmatrix.so ${f}

.PHONY: clean

clean:
	-rm -rf *.o *.so *.mod 
