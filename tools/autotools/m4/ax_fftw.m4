# -*- Autoconf -*-
#
# Copyright (C) 2019 Mathieu Taillefumier
#
# This file is part of the fftw software package. For license information,
# please see the COPYING file in the top-level directory of the fftw3 source
# distribution.
#

AU_ALIAS([ACX_FFTW], [AX_FFTW])
AC_DEFUN([ACX_FFTW],[
  dnl Init
  have_fftw3=no
  have_fftw3_mpi=no

  AC_ARG_WITH([fftw3],
              [AC_HELP_STRING([--with-fftw3=yes/no/PATH],
                              [compile with fftw3 support])],
              [fftw3_prefix=$withval],
              [fftw3_prefix="/usr/local"])
  AC_ARG_ENABLE([fftw3-mpi],
                [AC_HELP_STRING([--enable-fftw3-mpi], [check wether the fftw mpi module is present])])

  AS_IF([test "$with_fftw3" != "no"],
        [
        AS_IF([test "$fftw3_prefix" == "yes"], [fftw3_prefix="/usr/local"])

        ac_pkg_config_backup=$PKG_CONFIG_PATH
        PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$fftw3_prefix/lib/pkgconfig
        pc_fftw3_found="no"
        fftw3_ok="unknown"
        fftw3_has_hdrs="unknown"
        fftw3_has_libs="unknown"
        AS_IF([test "${PKG_CONFIG}yes" != "yes"], [echo "test"
                                                   PKG_CHECK_EXISTS([fftw3],
                                                                    [PKG_CHECK_MODULES(FFTW3,
                                                                                      [fftw3 >= 3.3.0],
                                                                                      [pc_fftw3_found="yes"
                                                                                      fftw3_has_hdrs="yes"
                                                                                      fftw3_has_libs="yes"],
                                                                                      [pc_fftw3_found="no"])],
                                                                    [pc_fftw3_found="no"])])
        PKG_CONFIG_PATH=$ac_pkg_config_backup

        dnl old fashion way if pkg-config is not found
        AS_IF([test "$pc_fftw3_found" == "no"],
              [
              dnl Backup environment
              saved_CPPFLAGS="${CPPFLAGS}"
              saved_LIBS="${LIBS}"

              dnl Prepare build parameters
              __ac_cpp_flags=$CPPFLAGS
              CPPFLAGS="${CPPFLAGS} -I${fftw3_prefix}/include"
              dnl Look for C includes
              AC_LANG_PUSH([C])
              AC_CHECK_HEADERS([fftw3.h],
                               [fftw3_has_hdrs="yes"],
                               [fftw3_has_hdrs="no"])
              dnl Look for C libraries and routines
              AC_CHECK_LIB(fftw3, [fftw_execute_dft], [fftw3_has_libs="yes"], [fftw3_has_libs="no"], [$PTHREAD_LIBS -lm])
              AC_LANG_POP([C])
              AS_IF([test "$fftw3_has_libs" == "yes" -a "$fftw3_has_hdrs" != "no"],
                    [FFTW3_CFLAGS=" -I$fftw3_prefix/include"
                    FFTW3_LIBS="-lfftw3"
                    FFTW3_THREAD_LIBS="-lfftw3_thread"
                    FFTW3_OMP_LIBS="-lfftw3_omp"])
              CPPFLAGS=$__ac_cpp_flags
              ])

        dnl Take final decision
        AC_MSG_CHECKING([whether we have basic FFTW3 support])
        AS_IF([test "${fftw3_has_hdrs}" = "yes" -a "${fftw3_has_libs}" = "yes"],
              [fftw3_ok="yes"
              AC_DEFINE([HAVE_FFTW3], 1, [fftw3 is found])],
              [fftw3_ok="no"])
        AC_MSG_RESULT([${fftw3_ok}])

        dnl check for mpi support
        AS_IF([test "${fftw3_has_libs}" = "yes" -a "$enable_fftw3_mpi" == "yes"],
              [
              AS_IF([test "$have_C_mpi" = "yes"],
                    [
                    AC_LANG_PUSH([C])
                    __ac_fftw3_cc=$CC
                    __ac_cpp_flags=$CPPFLAGS
                    CPPFLAGS="${CPPFLAGS} -I${fftw3_prefix}/include $FFTW3_CFLAGS"
                    CC=mpicc
                    AC_CHECK_HEADERS([fftw3-mpi.h],
                                     [fftw3_has_mpi_hdrs="yes"],
                                     [fftw3_has_mpi_hdrs="no"])
                    dnl Look for C libraries and routines
                    AC_CHECK_LIB(fftw3_mpi,
                                 [fftw_mpi_init],
                                 [fftw3_has_mpi_libs="yes"],
                                 [fftw3_has_mpi_libs="no"],
                                 [-lfftw3 $PTHREAD_LIBS -lm])
                    AC_LANG_POP([C])
                    CC=$__ac_fftw3_cc
                    CPPFLAGS=$__ac_cpp_flags])

            AS_IF([test "$enable_fftw3_mpi" == "yes"],
                  [
                   AC_MSG_CHECKING([whether the MPI FFTW3 libraries are installed])
                   AS_IF([test "$fftw3_has_mpi_libs" == "yes" -a "$fftw3_has_mpi_hdrs" == "yes"],
                         [fftw3_mpi_ok=yes; FFTW3_MPI_LIBS="-lfftw3_mpi"; AC_DEFINE([HAVE_FFTW3_MPI], 1, [fftw3_mpi is found])],
                         [fftw3_mpi_ok=no])
                   AC_MSG_RESULT([${fftw3_mpi_ok}])

                   ])
            ])
      ])
      have_fftw3=$fftw3_ok
      have_fftw3_mpi=$fftw3_mpi_ok
      AS_IF([test "$have_fftw3" == "yes"], [$1], [$2])
AC_SUBST([FFTW3_CFLAGS])
AC_SUBST([FFTW3_LIBS])
AC_SUBST([FFTW3_LDFLAGS])
AC_SUBST([FFTW3_THREAD_LIBS])
AC_SUBST([FFTW3_OMP_LIBS])
AC_SUBST([FFTW3_MPI_LIBS])
])
#
