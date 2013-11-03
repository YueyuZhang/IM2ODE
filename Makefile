.SUFFIXES: .F90 .o

FC = gfortran 
C = gcc
LINKER = $(FC)
FLAGS = -g 

MODULES = kinds.o constants.o parameters.o spacegroup.o \
lattice_mod.o Elec_Str.o tools.o init_mod.o \
init_struct.o sort_results.o differencial_evolution_1.o \
run_lammps.o test.o

ALL = $(MODULES)

default: DE

DE: $(ALL) 
	$(LINKER) $(FLAGS) -o de.x $(ALL)

$(MODULES):%.o:%.F90
	$(FC) -c $(FLAGS) $*.F90
clean:
	-rm -f de.x *.o *.mod

