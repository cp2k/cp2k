# -*- Autoconf -*-
AU_ALIAS([ACX_LIBVDWXC], [AX_LIBVDWXC])
AC_DEFUN([ACX_LIBVDWXC],
        [
        AC_REQUIRE([LX_FIND_MPI()])
        have_libvdwxc=no

        AC_ARG_WITH(libvdwxc,
                    [AC_HELP_STRING([--with-libvdxxc=yes/no/prefix],
                                   ["use the libvdwxc library for van der waals correction"])],
                    [libvdwxc_prefix="$withval"],
                    [libvdwxc_prefix="/usr/local"])

       AS_IF([test "$with_libvdwxc" != "no"],
             [
             __ac_cppflags_=$CPPFLAGS
             AS_IF([test "${libvdwxc_prefix}" == "yes"], [libvdwxc_prefix="/usr/local"])

             CPPFLAGS+=" -I${libvdwxc_prefix}/include "
             AC_LANG_PUSH([C])
             AC_CHECK_HEADER(vdwxc.h,
                             [libvdwxc_headers_found=yes],
                             [libvdwxc_headers_found=no])

             dnl libvdwxc has no mpi support so search for
             AC_CHECK_LIB(vdwxc,
                          [vdwxc_new],
                          [libvdwxc_found=yes],
                          [libvdwxc_found=no],
                          [-L$libvdwxc_prefix/lib -lvdwxc $FFTW3_LDFLAGS ${FFTW_LIBS} -lm])

             AS_IF([test "$have_C_mpi" = "yes"],
                   [
                    dnl do a self-consistency check. first the headers
                    AC_CHECK_HEADER(vdwxc_mpi.h,
                                   [libvdwxc_mpi_headers_found=yes],
                                   [libvdwxc_mpi_headers_found=no])

                    CPPFLAGS+=" ${MPI_CFLAGS}"
                    dnl the library with mpi
                    AC_CHECK_LIB(vdwxc,
                                 [vdwxc_mpi_init],
                                 [libvdwxc_mpi_found=yes; libvdwxc_found=yes],
                                 [libvdwxc_mpi_found=no],
                                 [-L$libint_prefix/lib -lvdwxc $FFTW3_LDFLAGS ${FFTW3_MPI_LIBS} ${MPI_LDFLAGS} ${FFTW3_LIBS} -lm])
                    dnl libvdwxc has no mpi support so search for
                    AC_LANG_POP([C])
                   ])
             CPPFLAGS=$__ac_cppflags_
             AC_MSG_CHECKING([whether the vdwxc library is installed and usable])
             AS_IF([test "${libvdwxc_found}" == "yes" -a "${libvdwxc_headers_found}" == "yes"],
                   [
                    AC_MSG_RESULT(yes)
                    AC_MSG_CHECKING([whether the vdwxc library is compiled with mpi support])

                    AS_IF([test "$libvdw_mpi_found" == "yes" -a "$libvdwxc_mpi_headers_found" == "yes"],
                          [AC_MSG_RESULT(yes)
                           AC_DEFINE([HAVE_LIBVDWXC_MPI], 1, [libvdwxc supports mpi])],
                          [AC_MSG_RESULT(no)])
                    LIBVDWXC_CFLAGS="-I${libvdwxc_prefix}/include"
                    LIBVDWXC_LIBS="-lvdwxc"
                    LIBVDWXC_LDFLAGS="-L${libvdwxc_prefix}/lib"
                    AC_DEFINE([HAVE_LIBVDWXC], 1, [libvdwxc is found])],
                   [AC_MSG_RESULT(no)])
               ])
        AC_SUBST(LIBVDWXC_CFLAGS)
        AC_SUBST(LIBVDWXC_LIBS)
        AC_SUBST(LIBVDWXC_LDFLAGS)
        AS_IF([test "$libvdwxc_found" == "yes" -a "$libvdwxc_headers_found" == "yes"], [$1], [$2])
])
#
