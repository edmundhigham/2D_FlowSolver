# makefile for FORTRAN 2D Flow Slover
# ME4 Final Year Project 2013/14
# E. Higham

# Compiler type
COMP = ifort

# Fortran Standard
STD =  -stand f03

# Enable FPP Pre-Processing
PREPROC = -fpp
# Generate Intermediates
#PREPROC += -save-temps bin/temps.f90

# Turn on optimisations
OPTIM = -O3
# Generate optimisation report
#OPTIMREP = -opt-report 2

# Enable Parallel Processing (generattion of multi-thread loops)
#PARAL = -parallel

# Allocate Arrays on Heap
ARRAY = -heap-arrays 1024

# Compile In Debug Mode
#DBUG = -g -warn all -e03 -fstack-security-check -traceback -check bounds

# Path to store .mod files
MODPATH = -module bin/

# COMBINE FLAGS
FFLAGS = $(STD) $(PREPROC) $(OPTIM) $(ARRAY) $(DBUG) $(PARAL) $(OPTIMREP)

# Name of Executable
NAME = Cavity65xy

bin/Precision.o : src/Precision.f90
	mkdir -p bin
	$(COMP) $(FFLAGS) $(MODPATH) $ -c $ $^ $ $ -o $ $@

bin/MeshGenmod.o : src/MeshGenmod.f90
	mkdir -p bin
	$(COMP) $(FFLAGS) $(MODPATH) $ -c $ $^ $ $ -o $ $@ 

bin/DerivativeMod.o : src/DerivativeMod.f90
	mkdir -p bin
	$(COMP) $(FFLAGS) $(MODPATH) $ -c $ $^ $ $ -o $ $@

bin/VortexRebound.o : src/VortexRebound.f90
	mkdir -p bin
	$(COMP) $(FFLAGS) $(MODPATH) $ -c $ $^ $ $ -o $ $@

bin/CavityFlow.o : src/CavityFlow.f90
	mkdir -p bin
	$(COMP) $(FFLAGS) $(MODPATH) $ -c $ $^ $ $ -o $ $@

bin/Poissonmod.o : src/Poissonmod.f90
	mkdir -p bin
	$(COMP) $(FFLAGS) $(MODPATH) $ -c $ $^ $ $ -o $ $@

bin/RKMod.o : src/RKMod.f90
	mkdir -p bin
	$(COMP) $(FFLAGS) $(MODPATH) $ -c $ $^ $ $ -o $ $@

bin/main.o : src/main.f90
	mkdir -p bin
	$(COMP) $(FFLAGS) $(MODPATH) $ -c $ $^ $ $ -o $ $@

all : bin/Precision.o bin/MeshGenmod.o bin/DerivativeMod.o bin/VortexRebound.o bin/CavityFlow.o bin/Poissonmod.o bin/RKMod.o bin/main.o
	$(COMP) $(FFLAGS) $(MODPATH) bin/main.o bin/RKMod.o bin/Poissonmod.o bin/CavityFlow.o bin/VortexRebound.o bin/DerivativeMod.o bin/MeshGenmod.o bin/Precision.o $ -o $ bin/$(NAME)


clean : 
	rm bin/*