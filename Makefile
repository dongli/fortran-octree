FC = gfortran

%.o: %.F90
	$(FC) -c $<

test.exe: octree.o test.o
	$(FC) -o $@ $^

.PHONY: clean

clean:
	rm -f *.o *.mod *.exe
