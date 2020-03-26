#
# setting for various compilers
# files should be standard f90 + getarg / iargc
#

##FC=xlf90
##FCFLAGS=-O2 -qsuffix=f=f90 -qmaxmem=128000 -bmaxdata:1000000000 -bmaxstack:128000000

##FC=ifc
##FCFLAGS=-O2 -Vaxlib

##FC=ifort
##FCFLAGS=-O2

##FC=pgf90
##FCFLAGS=-O2

##FC=g95
##FCFLAGS=-O2

FC=gfortran-mp-5
FCFLAGS=-O2

EXE=-o cubecruncher.x

#
# should not require modification
#
cubecruncher.x: utils.f90 periodic_table.f90 cubecruncher.f90 Makefile
	$(FC) $(FCFLAGS) $(EXE) utils.f90 periodic_table.f90 cubecruncher.f90

clean:
	rm -f *.MOD *.mod *.o core

distclean:
	rm -f *.MOD *.mod *.o core *.x
