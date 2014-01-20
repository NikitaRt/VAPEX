#####################################################################
# Makefile for compiling the Vapex program
#####################################################################
# $Id: Makefile 5 2014-01-13 11:44:50Z Sergey $
#####################################################################

#-------------- Define the compiler and its options -----------------
FC = ifort
FFLAGS =  -I. -I.. -80 -w -cm -C #-mpentiumpro -malign-double 
F90FLAGS = -I. -I.. -CU -CB -CV

OPT = -O0 -g

TARGET = WaterTest

LIBDIR = . #../lib/
LIBRARY =  

LIBS = -llapack -lm -lg2c 


MODS= 	../MOD_Globals.o\
	../MOD_Hydrogen.o\
	../MOD_WaterProps.o\
	../MOD_Disp.o\
	../MOD_Correl.o\
	../MOD_Evapor.o\
	../MOD_EqCoefs.o\
	$(NULL)

OBJS = WaterTest.o\
	$(NULL)


.SUFFIXES: .o .f .f90 .for .FOR

.f.o:
	$(FC) $(FFLAGS) $(OPT) -c $(*F).f
.f90.o:
	$(FC) $(F90FLAGS) $(OPT) -c $(*F).f90
.for.o:
	$(FC) $(FFLAGS) $(OPT) -c $(*F).for
.FOR.o:
	$(FC) $(FFLAGS) $(OPT) -c $(*F).FOR


$(TARGET): $(MODS) $(OBJS) #$(LIBDIR)lib$(LIBRARY).a
#	$(FC) $(CFLAGS) $(OPT) $(OBJS) -L$(LIBDIR) -o $(TARGET) -lm
	$(FC) $(CFLAGS) $(OPT) $(MODS) $(OBJS) -L$(LIBDIR) $(LIBS) -o $(TARGET) 

#lib$(LIBRARY).a:
#	cd ..;make

#--------------------------- use optimization -------------------------------
opt:
	@make -f Makefile OPT=-O3	
#--------------------------- remove the library -----------------------------
clean:
	@rm -f $(TARGET) *.dat *.f~ *.c *.c~ *.i~ *.o *.f90~ *.for~ Makefile~ core *.mod ifc* *.inp~
realclean:
	@make clean
