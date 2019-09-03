# -*- Autoconf -*-
# LICENCE : GPL
# Author : Mathieu Taillefumier
# Institution : ETH Zurich
#
# this m4 macro that check if intel mkl is installed on the system. by default
# it is on but it can be deactivated easily with the option --without-mkl. In
# that case it will search for others blas and lapack libraries using ACX_BLAS.
# The macro checks which compiler is used for compilation and set the different
# MKL_* variables accordingly. It also check for the presence of mpi and
# automatically setup the mpi dependent part of mkl accordingly.
#
# if mkl is found the function will set these variables for the 32 integers interface
#
# MKL_PARALLEL_LIBS : mkl libraries with omp support
# MKL_SEQUENTIAL_LIBS : mkl sequential libraries
# MKL_LDFLAGS : linker flags
# MKL_CFLAGS : compilation parameters for C/C++ compilers (with icc we can just use -mkl=sequential or parallel)
# MKL_FFLAGS : compilations parameters for fortran compilers
# MKL_CFLAGS_I64 : compilation parameters for C/C++ compilers (64 bits integers for indices)
# MKL_FFLAGS_I64 : compilations parameters for fortran compilers (64 bits integers for indices)
#
# and for 64 integers interface
#
# MKL_PARALLEL_LIBS_I64
# MKL_SEQUENTIAL_LIBS_I64
# MKL_CFLAGS_I64
# MKL_FFLAGS_I64
#
# all these variables are empty if mkl is not found so there is no arm putting
# @MKL_CFLAGS@, etc in your Makefile.am. For the 64 bit integers, simply add
# @MKL_CFLAGS_I64@ to the compiler flags. The libraries on the other end are
# specific.
#
# ============================================================================
#
# NB : the macro also defines BLAS_CFLAGS, BLAS_LDFLAGS, and BLAS_LIBS
# accordingly so that makefile can use these variables instead of the MKL_*
# variables for the 32 bits interface (or 64 if the --enable-mkl-int64 option on
# activated). Note that they by default using the sequential interface. This
# behavior can be changed with the --with-mkl-mode=parallel
#
#=============================================================================
#
#
# MKL also has scalapack, fftw, dft, etc integrated to it. it mostly add one
# parameter to the CFLAGS and FFLAGS and eventually add the needed libraries to
# MKL_PARALLEL_LIBS
#

AU_ALIAS([ACX_MKL], [AX_MKL])
AC_DEFUN([ACX_MKL], [
MKL_SCALAPACK_CFLAGS=""
MKL_PARALLEL_LIBS=""
MKL_SEQUENTIAL_LIBS=""
MKL_SCALAPACK_LIBS=""
MKL_PARALLEL_LIBS_I64=""
MKL_SEQUENTIAL_LIBS_I64=""
MKL_PARALLEL_FORTRAN_LIBS=""
MKL_SEQUENTIAL_FORTRAN_LIBS=""
MKL_CFLAGS_I64=""
MKL_FFLAGS_I64=""
MKL_CFLAGS=""
MKL_FFLAGS=""
have_mkl_scalapack=no
mkl_compiler_failure=no
mkl_lib_found=no
mklint64=no
AC_ARG_WITH(mkl,
            [AC_HELP_STRING([--with-mkl=yes/no/PATH],
                            [compile using intel mkl for blas/lapack, fftw])],
           [mkl_prefix=$withval],
           [mkl_prefix="/opt/intel/mkl"])

SCLAPACK_FOUND="no"
AS_IF([test "$with_mkl" != "no"],
      [

      dnl check for different environment variables such as MKLROOT, MKL_ROOT, etc
      AS_IF([test "$withval" == "yes"], [mkl_prefix="/opt/intel/mkl"])
      AS_IF([test "${MKLROOT+yes}" != "yes"], [mkl_prefix="$MKLROOT"])
      AS_IF([test "${MKLROOT+yes}" != "yes"], [mkl_prefix="$MKL_ROOT"])

      AC_ARG_WITH(mkl-mode,
                  [AS_HELP_STRING([--with-mkl-mode=sequential/parallel/both],
                                  [use multithreaded or single threaded version of mkl])],
                  [mklmode=$withval],
                  [mklmode="sequential"])

      AC_ARG_WITH(mkl-interface,
                  [AS_HELP_STRING([--with-mkl-interface=gcc/icc/pgi],
                  [indicate which interface to use for linking since the libraries are compiler dependent (the openmp part). By default it is given by the value of the CC variable])],
                  [mklinterface=$withval],
                  [mklinterface=${CC}])

      dnl check if need scalapack (will activate other dependencies)
      AC_ARG_ENABLE(mkl-mkl-fftw,
                  [AS_HELP_STRING([--enable-mkl-fftw],
                  [use mkl fftw3 wrapper])],
                  [mkl_fftw_wrapper="yes"],
                  [mkl_fftw_wrapper="no"])

    dnl enable 64 bit integer mode
    dnl check if need scalapack (will activate other dependencies)
    AC_ARG_ENABLE(mkl-int64,
                  [AS_HELP_STRING([ --enable-mkl-int64],
                                  [use mkl 64 bits integer interface])],
                  [mklint64="yes"],
                  [mklint64="no"])
      dnl check if need scalapack (will activate other dependencies)
      AC_ARG_ENABLE(mkl-scalapack,
                    [AS_HELP_STRING([ --enable-mkl-scalapack],
                                    [use mkl scalapack])],
                                    [AS_IF([test "${MPICC+yes}" == "yes"],
                                    [mklscalapack="yes"],
                                    [AC_MSG_ERROR([Could not find mpi. scalapack extension will be disabled])
                                    mklscalapack="no"])],
                                    [mklscalapack="no"])

      AS_IF([test $mklscalapack == "yes"],
            [
             AC_MSG_CHECKING([which version of mpi we use: ])

             temp=`mpirun --version 2>&1`
             AS_IF([echo $temp | grep -s  "HYDRA"], [mkl_mpi=mpich])
             AS_IF([echo $temp | grep -s  "Open MPI"], [mkl_mpi=openmpi])
             AS_IF([echo $temp | grep -s  "Intel(R) MPI"], [mkl_mpi=intelmpi])

             AS_IF([test ${mkl_mpi+yes} == "yes"],
                   [AC_MSG_RESULT([none. scalapack will be activated])],
                   [AC_MSG_RESULT($mkl_mpi)
                    SCALAPACK_LIBS_i32="-lmkl_scalapack_lp64"
                    SCALAPACK_LIBS_i64="-lmkl_scalapack_ilp64"
                    BLACS_LIBS_i32="-lmkl_blacs_${mkl_mpi}_lp64"
                    BLACS_LIBS_i64="-lmkl_blacs_${mkl_mpi}_ilp64"
                    SCALAPACK_LIBS=""
                    SCALAPACK_LDFLAGS=""
                    SCALAPACK_FFLAGS=""
                    SCALAPACK_FOUND="yes"
                    have_mkl_scalapack="yes"])
             ])



      dnl check that the library mkl_core and include files are there.
      AC_CHECK_LIB([mkl_sequential],
                   [mkl_free],
                   [MKL_SERIAL_CFLAGS="-m64 -I${mkl_prefix}/include "
                   MKL_PARALLEL_CFLAGS="-m64 -I${mkl_prefix}/include "
                   mkl_lib_found=yes],
                   [AC_MSG_NOTICE(["could not find the mkl libraries"])],
                   [-L$mkl_prefix/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lpthread -lm -ldl])

      AS_IF([test "$mkl_lib_found" == "yes"],
                  [AS_CASE("$mklinterface",
                           ["icc"], [
                                     cp_vendor=intel
                                     OMP_LIBS="-liomp5"
                                     MKL_CFLAGS_I64="-DMKL_ILP64 "
                                     MKL_CFLAGS="-m64 -I${mkl_prefix}/include"
                                     MKL_FFLAGS="-m64 -I${mkl_prefix}/include "
                                     MKL_FFLAGS_I64="-i8 "],
                           ["gcc"], [
                                     cp_vendor=gnu
                                     OMP_LIBS="-lgomp"
                                     MKL_CFLAGS=" -m64 -I${mkl_prefix}/include"
                                     MKL_CFLAGS_I64="-DMKL_ILP64 "
                                     MKL_FFLAGS="-m64 -I${mkl_prefix}/include"
                                     MKL_FFLAGS_I64="-fdefault-integer-8 -m64 -I${mkl_prefix}/include "],
                           ["pgcc"], [
                                      cp_vendor=pgi
                                      OMP_LIBS=" -pgf90libs -mp"
                                      MKL_SERIAL_CFLAGS+="-pgf90libs "
                                      MKL_PARALLEL_CFLAGS="-pgf90libs "
                                      MKL_CFLAGS_I64="-i8 "],
                           [AC_MSG_ERROR([the mkl parameters for the compiler "$CC" are not known.])])
                           MKL_LDFLAGS="-L$mkl_prefix/lib/intel64"
                           MKL_SEQUENTIAL_LIBS_I64="-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core"
                           MKL_PARALLEL_LIBS_I64="$SCALAPACK_LIBS_i64 -lmkl_intel_ilp64 -lmkl_${cp_vendor}_thread -lmkl_core $BLACS_LIBS_i64 $OMP_LIBS"
                           MKL_SEQUENTIAL_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
                           MKL_PARALLEL_LIBS="$SCALAPACK_LIBS_i32  -lmkl_intel_lp64 -lmkl_${cp_vendor}_thread -lmkl_core  $BLACS_LIBS_i32 $OMP_LIBS"
                           BLAS_CFLAGS=$MKL_CFLAGS
                           BLAS_LDFLAGS=$MKL_LDFLAGS
                           BLAS_FFLAGS=$MKL_FFLAGS
                           AS_IF([test "$mklmode" = "parallel"],
                                 [AS_IF([test "$mklint64" != "yes"], [BLAS_LIBS=$MKL_PARALLEL_LIBS], [BLAS_LIBS=$MKL_PARALLEL_LIBS_I64])],
                                 [
                                 AS_IF([test "$mklint64" != "yes"], [BLAS_LIBS=$MKL_SEQUENTIAL_LIBS], [BLAS_LIBS=$MKL_SEQUENTIAL_LIBS_I64])])
                          AS_IF([test "$mklint64" = "yes"], [BLAS_CFLAGS="${BLAS_CFLAGS} ${MKL_CFLAGS_I64}"
                                                             BLAS_FFLAGS="${MKL_FFLAGS_I64}"])
                 ])
      ])

      AC_SUBST([MKL_PARALLEL_LIBS])
      AC_SUBST([MKL_SEQUENTIAL_LIBS])
      AC_SUBST([MKL_PARALLEL_FORTRAN_LIBS])
      AC_SUBST([MKL_SEQUENTIAL_FORTRAN_LIBS])
      AC_SUBST([MKL_CFLAGS])
      AC_SUBST([MKL_FFLAGS])
      AC_SUBST([MKL_PARALLEL_LIBS_I64])
      AC_SUBST([MKL_SEQUENTIAL_LIBS_I64])
      AC_SUBST([MKL_PARALLEL_FORTRAN_LIBS_I64])
      AC_SUBST([MKL_SEQUENTIAL_FORTRAN_LIBS_I64])
      AC_SUBST([MKL_CFLAGS_I64])
      AC_SUBST([MKL_FFLAGS_I64])
      AC_SUBST([MKL_LDFLAGS])
      AC_SUBST([BLAS_CFLAGS])
      AC_SUBST([BLAS_LDFLAGS])
      AC_SUBST([BLAS_LIBS])
      AC_SUBST([BLAS_FFLAGS])
      AC_SUBST([OMP_LIBS])

      AS_IF([test "${mkl_fftw_wrapper+yes}" = "yes"],
            [AC_DEFINE([HAVE_MKL_FFTW], 1, [use the mkl fftw wrapper])])

      AS_IF([test "$mkl_lib_found" == "yes"],
            [AC_DEFINE([HAVE_MKL], 1, [mkl library found])
             AC_DEFINE([HAVE_CBLAS], 1, [mkl cblas library found])
            ])

      AS_IF([test "$mkl_lib_found" == "yes"],
            [$1], [$2])
])
