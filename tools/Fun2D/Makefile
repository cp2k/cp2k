#
# requires a sufficiently recent gfortran and mpfr ....
#
# if the openmp version is to be used, mpfr needs to be thread safe
# this requires an mpfr configured with './configure --enable-thread-safe ...'
#
FC       = gfortran-4.3 -ffree-line-length-512 -g -fopenmp
FFLAGS   = -O3 -march=native -funroll-loops
LIBSMPFR = -lmpfr -lgmp

mpfr_cutoff_gamma.o: mpfr_cutoff_gamma.f90
	$(FC) -c mpfr.f90 mpfr_cutoff_gamma.f90

evalf: mpfr_cutoff_gamma.o
	$(FC) -o evalf.x mpfr.f90 mpfr_cutoff_gamma.o evalf.f90 $(LIBSMPFR)

drvrec_mpfr: mpfr_cutoff_gamma.o
	$(FC) -o drvrec_mpfr.x mpfr.f90 padua2_mpfr.f90 mpfr_cutoff_gamma.o drvrec_mpfr.f90  $(LIBSMPFR)

recurse: drvrec_mpfr evalf
	$(FC) -o recurse.x recurse.f90

t_c_g0: recurse 
	./recurse.x > t_c_g0.f90

tester: t_c_g0
	$(FC) $(FFLAGS) -o tester.x kinds.f90 t_c_g0.f90 tester.f90

test: tester
	./tester.x

clean: 
	rm -f *.o *.x *.mod *.MOD

realclean: clean 
	rm -f *.dat *.in *.out t_c_g0.f90 fort.*
