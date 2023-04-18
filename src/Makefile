# makefile for program GranPrix

OBJ = breaks.o \
defect.o \
findforces.o \
force.o \
generate.o \
init.o \
inputReader.o \
linklist.o \
mycommons.o \
postScriptTool.o \
principalstress.o \
pssubs.o \
rotate.o \
rw_conf.o \
smoothNonExtended.o \
updatePressure.o \
outputWriter.o \
vel_ver.o

SRC = mycommons.f90 \
inputReader.f90 \
outputWriter.f90 \
findforces.f90 \
generate.f90 \
linklist.f90 \
principalstress.f90 \
breaks.f90 \
defect.f90 \
force.f90 \
init.f90 \
postScriptTool.f90 \
rotate.f90 \
pssubs.f90 \
rw_conf.f90 \
smoothNonExtended.f90 \
updatePressure.f90 \
vel_ver.f90 

FC     =  gfortran
CFLAGS = -O3 -std=f2008 -Wall 
CFLAGS2 = -O3
EXE    = GranFrixrm

all: clean $(OBJ) GranLayer.f90
	$(FC) $(CFLAGS) -o $(EXE) GranLayer.f90 $(OBJ)

intel: clean GranLayer.f90
	ifort $(CFLAGS2) -c $(SRC)
	ifort $(CFLAGS2) -o $(EXE) GranLayer.f90 $(OBJ)

$(OBJ):
	$(FC) $(CFLAGS) -c $(SRC)

clean-mod:
	rm -f *.mod
	
clean-obj:
	rm -f *.o

clean: clean-mod clean-obj 
	rm -f GranFrixrm
