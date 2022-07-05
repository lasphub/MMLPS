f = fakebondmatrix.f90 simplefun.f90 

FC      = ifort
FCFLAGS = -shared -fpic -O2 -xHost #-openmp 

bondmatrix.so: ${f}
	${FC} ${FCFLAGS}  -o checkminbond.so ${f}

.PHONY: clean

clean:
	-rm -rf *.o *.so *.mod 
