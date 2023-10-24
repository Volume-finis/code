F90 = gfortran
FLAGS = -O0 -g -ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all -Wall
EXE = executable
all: $(EXE)

# ce qui est en bleu est une cible
# le run est juste un appel

$(EXE): mod_precision.o mod_maillage.o mod_sortie.o Principal_Program.o
	$(F90) $(FLAGS) -o executable mod_precision.o mod_maillage.o mod_sortie.o Principal_Program.o
mod_precision.o : mod_precision.f90
	$(F90) $(FLAGS) -c mod_precision.f90
mod_maillage.o : mod_precision.o mod_maillage.f90
	$(F90) $(FLAGS) -c mod_maillage.f90
mod_sortie.o : mod_precision.o mod_maillage.o mod_sortie.f90
	$(F90) $(FLAGS) -c mod_sortie.f90
Principal_Program.o : mod_precision.o mod_maillage.o mod_sortie.o Principal_Program.f90
	$(F90) $(FLAGS) -c Principal_Program.f90

clean :.
	rm *.o *.mod *.vtk
