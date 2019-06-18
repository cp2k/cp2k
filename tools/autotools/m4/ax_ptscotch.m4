# -*- Autoconf -*-

AU_ALIAS([ACX_PTSCOTCH], [AX_PTSCOTCH])
AC_DEFUN([ACX_PTSCOTCH],
         [
         AC_REQUIRE([LX_FIND_MPI])
         AC_REQUIRE([ACX_PARMETIS])
         AC_ARG_WITH(ptscotch,
                     [AC_HELP_STRING([--with-ptscotch=yes/no/prefix], [enable ptscotsch support (default yes)])],
                     [ptscotch_prefix="$withval"],
                     [ptscotch_prefix="/usr/local"])
         ptscotch_found=no
         ptscotch_headers_found=no
         AS_IF([test "$with_ptscotch" != "no" -a "$have_C_mpi" != "yes"],
               [
                             __ac_cppflags_=$CPPFLAGS
                             CPPFLAGS+=" -I${ptscotch_prefix}/include $MPI_CFLAGS"
                             __ac_cc_backup=${CC}
                             CC=${MPICC}
                             AC_CHECK_HEADER(ptscotch.h,
                                             [ptscotch_headers_found=yes],
                                             [ptscotch_headers_found=no])
                             dnl ptscotch depends theoretically on parmetis but
                             dnl some distribution do not link against it for licencing reasons.
                             dnl so check if ptscotch can be linked independently or not.
                             AC_CHECK_LIB(ptscotch,
                                          [SCOTCH_dgraphInit],
                                          [ptscotch_found=yes],
                                          [ptscotch_found=no],
                                          [-L${ptscotch_prefix}/lib -lptscotcherr -lscotchmetis -lscotch  -lscotcherr  -lptscotchparmetis $MPI_CLDFLAGS -lm])

                             AC_MSG_CHECKING([whether ptscotsch depends explicitly on parmetis])
                             AS_IF([test $ptscotch_found == "yes"], [AC_MSG_RESULT(no)], [AC_MSG_RESULT(yes)])
                             CPPFLAGS=$__ac_cppflags_
                             CC=$__ac_cc_backup
                             PTSCOTCH_CFLAGS="-I${ptscotch_prefix}/include"
                             PTSCOTCH_LDFLAGS="-L${ptscotch_prefix}/lib"
                             PTSCOTCH_LIBS="-lptscotch -lptscotcherr -lscotchmetis -lscotch  -lscotcherr  -lptscotchparmetis"
                             AC_DEFINE([HAVE_PTSCOTCH], 1, ["ptscotch is found"])
                         ])
               ])

         AC_SUBST(PTSCOTCH_CFLAGS)
         AC_SUBST(PTSCOTCH_LIBS)
         AC_SUBST(PTSCOTCH_LDFLAGS)

         AS_IF([test "$ptscotch_found" == "yes" -a "$ptscotch_headers_found" == "yes"],
               [$1],
               [$2])
])
