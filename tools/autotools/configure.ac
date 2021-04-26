dnl #                                               -*- Autoconf -*-
dnl # Process this file with autoconf to produce a configure script.

AC_INIT([CP2K], [6.1.0], [mathieu.taillefumier@free.fr])
AC_PREREQ([2.52])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADERS([config.h])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_ISC_POSIX
AM_INIT_AUTOMAKE([1.6 dist-bzip2 silent-rules subdir-objects])

dnl Small macro to print out a nice, pretty section title.
define(SECTION_TITLE,
[
echo >& AS_MESSAGE_FD
echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
echo '   $1' >& AS_MESSAGE_FD
echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
])

LT_INIT([disable_shared])
AC_DISABLE_SHARED
PKG_PROG_PKG_CONFIG([0.10])

dnl list of all possible modules
all_options="mpi cuda openmp libint2 blas mkl smm xsmm spglib sirius libxc libvdwxc elpa scalapack superlu_dist ptscotch ParMETIS fftw3 gsl pexsi"

dnl # Checks for programs.
AC_PROG_CXX(g++ icc clang pgcc xlc++)
AC_PROG_CC(gcc icc clang pgcc xlc)
AC_PROG_CC_STDC
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_PROG_CC_C99
dnl # Checks for library functions.
AC_FUNC_MALLOC

SECTION_TITLE([Fortran])
AC_PROG_FC([gfortran ifort xlf pgf90])

AC_FC_FREEFORM([], [])

AC_F77_LIBRARY_LDFLAGS
SECTION_TITLE([multithreading])
dnl # check for pthreads. It is needed by almost all libraries
AX_PTHREAD([options="$options thread"
CC="$PTHREAD_CC"], [])

SECTION_TITLE([OpenMP support])
dnl check openmp support
AX_OPENMP([options="$options openmp"
DFLAGS_SMP=""], [])
AC_SUBST([OPENMP_CFLAGS])
AC_SUBST([OPENMP_CXXFLAGS])
AC_SUBST([OPENMP_FFLAGS])
AC_SUBST([OPENMP_LIBS])
AC_SUBST([OPENMP_LDFLAGS])
dnl optional modules to be compiled against
DFLAGS=""

SECTION_TITLE([Fortran compiler warnings options])
AC_LANG_PUSH([Fortran])
AX_CHECK_COMPILE_FLAG(-m64)
AX_CHECK_COMPILE_FLAG("-Werror=aliasing", [WFLAGS+="-Werror=aliasing "])
AX_CHECK_COMPILE_FLAG("-Werror=ampersand", [WFLAGS+="-Werror=ampersand "])
AX_CHECK_COMPILE_FLAG("-Werror=c-binding-type", [WFLAGS+="-Werror=c-binding-type "])
AX_CHECK_COMPILE_FLAG("-Werror=intrinsic-shadow", [WFLAGS+="-Werror=intrinsic-shadow "])
AX_CHECK_COMPILE_FLAG("-Werror=intrinsics-std", [WFLAGS+="-Werror=intrinsics-std "])
AX_CHECK_COMPILE_FLAG("-Werror=line-truncation", [WFLAGS+="-Werror=line-truncation "])
AX_CHECK_COMPILE_FLAG("-Werror=tabs", [WFLAGS+="-Werror=tabs "])
AX_CHECK_COMPILE_FLAG("-Werror=target-lifetime", [WFLAGS="-Werror=target-lifetime "])
AX_CHECK_COMPILE_FLAG("-Werror=underflow", [WFLAGS+="-Werror=underflow "])
AX_CHECK_COMPILE_FLAG("-Werror=unused-but-set-variable", [WFLAGS+="-Werror=unused-but-set-variable "])
AX_CHECK_COMPILE_FLAG("-Werror=unused-variable", [WFLAGS+="-Werror=unused-variable "])
AX_CHECK_COMPILE_FLAG("-Werror=unused-dummy-argument", [WFLAGS+="-Werror=unused-dummy-argument "])
AX_CHECK_COMPILE_FLAG("-Werror=conversion", [WFLAGS+="-Werror=conversion "])
AX_CHECK_COMPILE_FLAG("-Werror=zerotrip", [WFLAGS+="-Werror=zerotrip "])
AX_CHECK_COMPILE_FLAG("-Werror=uninitialized", [WFLAGS+="-Werror=uninitialized "])
AX_CHECK_COMPILE_FLAG("-Wno-maybe-uninitialized", [WFLAGS+="-Wno-maybe-uninitialized "])
AX_CHECK_COMPILE_FLAG("-Wuse-without-only", [WFLAGS+="-Wuse-without-only "])
echo $WFLAGS
AC_LANG_POP([Fortran])
AC_SUBST(WFLAGS)

SECTION_TITLE(Compilation parameters)
PSMP_CFLAGS="-O2 -mtune=native"
SOPT_FLAGS="-g -O2"
PSMP_FFLAGS="-O2 -mtune=native"
AC_SUBST([PSMP_FFLAGS])
AC_SUBST([PSMP_CFLAGS])
AC_SUBST([SOPT_FLAGS])

SECTION_TITLE([mpi support])

AC_LANG_PUSH([Fortran])
LX_FIND_MPI([options="$options mpi"], [AC_MSG_ERROR([cp2k needs mpi])])
AC_LANG_POP([Fortran])

AC_LANG_PUSH([C++])
LX_FIND_MPI([], [AC_MSG_ERROR([cp2k needs C++ mpi support])])
AC_LANG_POP([C++])

AC_LANG_PUSH([C])
LX_FIND_MPI([], [AC_MSG_ERROR([cp2k needs C mpi support])])
AC_LANG_POP([C])
MPI_DFLAGS="-D__parallel"

AC_SUBST([MPICC])
AC_SUBST([MPICXX])
AC_SUBST([MPIFC])
AC_SUBST([MPI_DFLAGS])

SECTION_TITLE([GPU support])
dnl check for cuda support
CUDA_DFLAGS=""
ACX_CUDA([options="$options cuda"
          CUDA_DFLAGS="-D__DBCSR_ACC -D__PW_CUDA"], [])
AC_SUBST(GPUVER)
AC_SUBST(CUDA_DFLAGS)

SECTION_TITLE([BLAS, LAPACK, MKL, etc...])
dnl check for mkl. Note that blas is checked only if mkl is not found

MKL_DFLAGS=""
ACX_MKL([AC_DEFINE([HAVE_BLAS], [1], [use blas])
         options="$options mkl" MKL_DFLAGS="-D__MKL"], [ACX_BLAS([options="$options blas"], [options="$options blas"], [AC_MSG_ERROR([a library implementing blas is required])])])
AC_SUBST(BLAS_CFLAGS)
AC_SUBST(BLAS_LDFLAGS)
AC_SUBST(BLAS_LIBS)
AC_SUBST(BLAS_FFLAGS)
AC_SUBST(MKL_DFLAGS)
SECTION_TITLE([small matrix multiplication library])
dnl check for libsmm
ACX_LIBXSMM([options="$options xsmm"
             LIBXSMM_DFLAGS="-D__LIBXSMM"], [])

ACX_LIBSMM([options="$options smm"
             LIBSMM_DFLAGS="-D__SMM"], [])
AC_SUBST(LIBSMM_DFLAGS)
AC_SUBST(LIBXSMM_DFLAGS)

SECTION_TITLE([SCALAPACK])
dnl check for scalapack
ACX_SCALAPACK([options="$options scalapack"; SCALAPACK_DFLAGS="-D__SCALAPACK"], [])
AC_SUBST(SCALAPACK_DFLAGS)

SECTION_TITLE([elpa])
dnl check for elpa
ACX_ELPA([options="$options elpa"; ELPA_DFLAGS="-D__ELPA=$elpa_version"], [])
AC_SUBST(ELPA_DFLAGS)


SECTION_TITLE([gsl library])
# detect the gsl if asked for
AX_PATH_GSL([1.12], [GSL_DFLAGS="-D__GSL"
             options="$options gsl"])
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_SUBST(GSL_DFLAGS)


SECTION_TITLE([FFTW3])
dnl check for fftw
AC_ARG_ENABLE(fftw-mkl-wrapper,
              [AC_HELP_STRING([--enable-fftw-mkl-wrapper],
                              [use the fftw wrappoer included in mkl])
              ],
              [],
              [])

dnl if enabled check that we have MKL detected, otherwise fall back to searching fftw
AS_IF([test "$have_mkl" = "no"], [enable_fftw_mkl_wrapper=no])

dnl fall back to fftw
AS_IF([test $enable_fftw_mkl_wrapper=="no"], [ACX_FFTW([options="$options fftw3"], [])])
FFTW3_DFLAGS="-D__FFTW3"
AC_SUBST(FFTW3_DFLAGS)

SECTION_TITLE([HDF5])
dnl Check for serial HDF5 support
AX_LIB_HDF5([serial], [])
AS_IF([test "${HDF5_CC}yes" != "yes"], [HDF5_DFLAGS="-D__HDF5"])
dnl Check for parallel HDF5 support
dnl AX_LIB_HDF5([parallel])
AC_SUBST(HDF5_DFLAGS)

SECTION_TITLE([spglib])
dnl check for spglib
ACX_SPGLIB([options="$options spglib"; SPGLIB_DFLAGS="-D__SPGLIB"], [])
AC_SUBST(SPGLIB_DFLAGS)

SECTION_TITLE([libint2])
dnl check for libint
ACX_LIBINT([options="$options libint2"; LIBINT2_DFLAGS="-D__LIBINT"], [])
AC_SUBST(LIBINT2_DFLAGS)

SECTION_TITLE([libxc])
dnl check for libxc
ACX_LIBXC([options="$options libxc"; LIBXC_DFLAGS="-D__LIBXC"], [])
AC_SUBST(LIBXC_DFLAGS)

SECTION_TITLE([libvdwxc (Van der Walls functionals)])
dnl check for libvdwxc
ACX_LIBVDWXC([options="$options libvdwxc"; LIBVDWXC_DFLAGS="-D__LIBVDWXC"], [])
AC_SUBST(LIBVDWXC_DFLAGS)

SECTION_TITLE([ParMETIS, PTSCOTCH and SuperLU])
dnl check for parmetis
ACX_PARMETIS([options="$options ParMETIS"; PARMETIS_DFLAGS="-D__PARMETIS"], [])
AC_SUBST(PARMETIS_DFLAGS)

dnl check for PTSCOTCH
ACX_PTSCOTCH([options="$options ptscotch"; PTSCOTCH_DFLAGS="-D__PTSCOTCH"], [])
AC_SUBST(PTSCOTCH_DFLAGS)

dnl check for superlu
ACX_SUPERLU([options="$options superlu_dist"; SUPERLU_DIST_DFLAGS="-D__SUPERLU_DIST"], [])
AC_SUBST(SUPERLU_DIST_DFLAGS)

SECTION_TITLE([PEXSI])
dnl check for pexsi
ACX_PEXSI([options="$options pexsi"; PEXSI_DFLAGS="-D__LIBPEXSI"], [])
AC_SUBST(PEXSI_DFLAGS)

SECTION_TITLE([SIRIUS])
dnl check for sirius
ACX_SIRIUS([options="$options sirius"; SIRIUS_DFLAGS="-D__SIRIUS"], [])
AC_SUBST(SIRIUS_DFLAGS)

SECTION_TITLE([Configuration results])

AC_SUBST(SOPT_CFLAGS)
AC_SUBST(SOPT_FFLAGS)

echo "Options:" >& AS_MESSAGE_FD
echo "$options" >& AS_MESSAGE_FD
echo >& AS_MESSAGE_FD

echo "Missing options:" >& AS_MESSAGE_FD
echo "$options" "$all_options" | tr ' ' '\n' | sort | uniq -u | tr '\n' ' ' >& AS_MESSAGE_FD


echo -e "\n\n" >& AS_MESSAGE_FD
echo "CXXFLAGS=$CXXFLAGS" >& AS_MESSAGE_FD
echo "LDFLAGS=$LDFLAGS" >& AS_MESSAGE_FD
echo "LIBS=$LIBS" >& AS_MESSAGE_FD

AS_IF([test "$ac_have_cuda_" = "yes"], [AC_CONFIG_FILES([local_cuda.sopt local_cuda.psmp])])
AC_CONFIG_FILES([local.sopt local.psmp local.ssmp])
AC_OUTPUT()

#           src/Makefile
#           src/bin/Makefile
#           src/lib/Makefile],
#           []
#           )
