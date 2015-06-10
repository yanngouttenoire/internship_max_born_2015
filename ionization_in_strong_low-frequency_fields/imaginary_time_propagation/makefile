# A simple hand-made makefile for a package including applications
# built from Fortran 90 sources, taking into account the usual
# dependency cases.

# This makefile works with the GNU make command, the one find on
# GNU/Linux systems and often called gmake on non-GNU systems, if you
# are using an old style make command, please see the file
# Makefile_oldstyle provided with the package.

# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = gfortran
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS = -I ~/Documents/fftw-3.3.4
LDFLAGS =  -lfftw3
EXEC=exec.out
SOURCE= $(wildcard *.f90)
OBJET= $(SOURCE:.cpp=.o)
HEADER= $(wildcard *.h)

run: exec.out
	@ ./exec.out

exec.out: $(OBJET)
	@ $(FC) -o $@ $^ $(FCFLAGS) $(LDFLAGS)

main.o: $(HEADER)

%.o: %.c
	@ $(FC) -o $@ -c $< $(FCFLAGS) $(LDFLAGS)

.PHONY: clean mrproper run 

clean:
	@ rm -rf *.o

mrproper: clean
	@ rm -rf $(EXEC)

		

