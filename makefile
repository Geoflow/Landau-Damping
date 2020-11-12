#
## Compilateur
#
FC = gfortran
#
## Options de compilation
#
FFLAGS = -O3 -g -fbacktrace -ffixed-line-length-none 



#
## Editeur de liens
#
LD = $(FC)
#
## Options d'édition des liens
#
LDFLAGS = -O2

#
#
# Fichiers objets
#
#
OBJS = numerics.o utils.o  main.o godunov.o semilag.o
MODS = numerics.mod utils.mod godunov.mod semilag.mod
#
#
## Programme
#
EXEC = exe

#
## Règles de compilation





all: $(EXEC)

$(EXEC): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXEC)


numerics.o: numerics.f90 
	$(FC) -c $(FFLAGS) numerics.f90



utils.o:numerics.o utils.f90 
	$(FC) -c $(FFLAGS) utils.f90


godunov.o: numerics.o utils.o godunov.f90
	$(FC) -c $(FFLAGS)  godunov.f90

semilag.o: semilag.f90  numerics.o utils.o 
	$(FC) -c $(FFLAGS)  semilag.f90 


main.o: godunov.o numerics.o utils.o main.f90
	$(FC) -c $(FFLAGS)  main.f90

#
#

clean:
	\rm -f $(OBJS) $(MODS)
#
#

