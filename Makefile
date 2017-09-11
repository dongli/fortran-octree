FC = gfortran
FFLAGS = -O3 #-fopenmp

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

test.exe: octree.o test.o
	$(FC) $(FFLAGS) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o *.mod gmon.out *.exe
