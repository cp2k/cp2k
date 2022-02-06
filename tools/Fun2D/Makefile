#
# requires a sufficiently recent gfortran and mpfr ....
#
# unfortunately, compiling can take ~1-2h
#
# if the openmp version is to be used, mpfr needs to be thread safe
# this requires an mpfr configured with './configure --enable-thread-safe ...'
#
YUKAWA   = -D__YUKAWA
T_C_G0   = -D__T_C_G0
DFLAGS   = $(T_C_G0)
CPP      = cpp
CPPFLAGS = -C -traditional $(DFLAGS)
FC       = gfortran-4.3 -ffree-line-length-512 -g -fopenmp -ffree-form $(DFLAGS)
# mpfr needs to be thread safe if used with OMP
LIBS     = -lmpfr -lgmp
# these flags can greatly affect the performance of the generated code
# but should only be used to compile the tester program, compilation is otherwise way too slow
FFLAGS   = -O3 -march=native -funroll-loops

.SUFFIXES: .f90 .x .o

%.o: %.f90
ifneq ($(CPP),)
	$(CPP) $(CPPFLAGS) $< > $*.f
	$(FC) -c $*.f
else
	$(FC) -c $<
endif

drvrec_mpfr.o : drvrec_mpfr.f90 mpfr.o mpfr.o functions.o padua2_mpfr.o
evalf.o : evalf.f90 mpfr_cutoff_gamma.o
functions.o : functions.f90 function_types.o mpfr_yukawa.o mpfr_cutoff_gamma.o mpfr.o mpfr.o
function_types.o : function_types.f90
kinds.o : kinds.f90
mpfr_cutoff_gamma.o : mpfr_cutoff_gamma.f90 mpfr.o mpfr.o
mpfr.o : mpfr.f90
mpfr_yukawa.o : mpfr_yukawa.f90 mpfr.o mpfr.o
padua2_mpfr.o : padua2_mpfr.f90 mpfr.o mpfr.o
recurse.o : recurse.f90 function_types.o
t_c_g0.o : t_c_g0.f90 kinds.o
yukawa.o : yukawa.f90 kinds.o

.PHONY: recurse tester test
recurse: recurse.x
tester: tester.x

#
# the test program depends on the choice of -DFLAGS
#
ifeq ($(DFLAGS),$(T_C_G0))
TESTEROBJ=kinds.o t_c_g0.o
else
TESTEROBJ=kinds.o yukawa.o
endif
tester.o: tester.f90 $(TESTEROBJ)
tester.x: tester.o $(TESTEROBJ)
	$(FC) $(FFLAGS) -o tester.x tester.o $(TESTEROBJ)

#
# drvrec_mpfr.x
#
DRVRECOBJ=mpfr.o functions.o padua2_mpfr.o function_types.o mpfr_yukawa.o mpfr_cutoff_gamma.o
drvrec_mpfr.x: drvrec_mpfr.o $(DRVRECOBJ)
	$(FC) -o drvrec_mpfr.x drvrec_mpfr.o $(DRVRECOBJ) $(LIBS)

#
# recurse.x
#
RECURSEOBJ=function_types.o
recurse.x: recurse.o $(RECURSEOBJ)
	$(FC) -o recurse.x recurse.o $(RECURSEOBJ)

#
# evalf is a quick hack
#
evalf.x: evalf.o mpfr_cutoff_gamma.o mpfr.o
	$(FC) -o drvrec_mpfr.x evalf.o mpfr_cutoff_gamma.o mpfr.o $(LIBS)

# generated files need special dependencies
t_c_g0.f90: recurse.x drvrec_mpfr.x
	./recurse.x > t_c_g0.f90
yukawa.f90: recurse.x drvrec_mpfr.x
	./recurse.x > yukawa.f90

test: tester.x
	./tester.x

clean:
	rm -f *.o *.x *.mod *.MOD *.f

realclean: clean
	rm -f *.dat *.in *.out t_c_g0.f90 fort.* yukawa.f90
