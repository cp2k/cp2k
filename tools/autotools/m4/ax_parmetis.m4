# -*- Autoconf -*-

AU_ALIAS([ACX_PARMETIS], [AX_PARMETIS])
AC_DEFUN([ACX_PARMETIS], [

AC_REQUIRE([LX_FIND_MPI])

have_parmetis=no
AC_ARG_WITH(parmetis,
            [AC_HELP_STRING([--with-parmetis=yes/no/prefix], [use the parmetis library (default yes) ])],
            [libparmetis_prefix="$withval"],
            [libparmetis_prefix="/usr/local"])

AS_IF([test "$with_libparmetis" != "no" -a $have_C_mpi == "yes"],
      [
            AS_IF([test "$libparmetis_prefix" == "yes"], [libparmetis_prefix="/usr/local"])
            __ac_cppflags_=$CPPFLAGS
            __ac_cc_compiler=$CC
            CC=${MPICC}
            CPPFLAGS+=" -I$libparmetis_prefix/include $MPI_CFLAGS"

            AC_LANG_PUSH([C])
            AC_CHECK_HEADERS([parmetis.h],
                             [parmetis_has_hdrs="yes"],
                             [parmetis_has_hdrs="no"])
            AC_CHECK_LIB(parmetis,
                         [ParMETIS_V3_PartKway],
                         [parmetis_found=yes],
                         [parmetis_found=no],
                         [-L$libparmetis_prefix/lib $MPI_CLDFLAGS $PTHREAD_LIBS -lm])
            AC_LANG_POP([C])

            CC=$__ac_cc_compiler
            CPPFLAGS=$__ac_cppflags_
            AS_IF([test "$parmetis_found" == "yes"],
                  [LIBPARMETIS_CFLAGS="-I${libparmetis_prefix}/include"
                   LIBPARMETIS_LIBS="-lparmetis"
                   LIBMETIS_CFLAGS="-I${libparmetis_prefix/include}"
                   LIBMETIS_LIBS="-lmetis"
                   LIBPARMETIS_LDFLAGS="-L${libparmetis_prefix}/lib"
                   AC_DEFINE([HAVE_PARMETIS], 1, ["parmetis is found"])
                   AC_DEFINE([HAVE_METIS], 1, ["metis is found"])
                   have_parmetis=yes])
      ])

AC_SUBST(PARMETIS_CFLAGS)
AC_SUBST(PARMETIS_LIBS)
AC_SUBST(PARMETIS_LDFLAGS)
AC_SUBST(METIS_CFLAGS)
AC_SUBST(METIS_LIBS)

AS_IF([test "$libparmetis_found" == "yes"], [$1], [$2])
])
