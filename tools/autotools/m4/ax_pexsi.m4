# -*- Autoconf -*-

AU_ALIAS([ACX_PEXSI], [AX_PEXSI])
AC_DEFUN([ACX_PEXSI],
         [
         have_pexsi="yes"
         AC_REQUIRE([LX_FIND_MPI])
         AC_ARG_WITH(pexsi,
                     [AC_HELP_STRING([--with-pexsi=yes/no/prefix], [enable pexsi support (default yes)])],
                     [pexsi_prefix="$withval"],
                     [pexsi_prefix="/usr/local"])
         pexsi_found=no
         pexsi_headers_found=no
         AS_IF([test "$with_pexsi" != "no" -a "$have_C_mpi" != "yes"],
               [
                             __ac_cppflags_=$CPPFLAGS
                             CPPFLAGS+=" -I${pexsi_prefix}/include $MPI_CFLAGS"
                             __ac_cc_backup=${CC}
                             CC=${MPICC}
                             AC_CHECK_HEADER(c_pexsi_interface.h,
                                             [pexsi_headers_found=yes],
                                             [pexsi_headers_found=no])
                             dnl ptscotch depends theoretically on parmetis but
                             dnl some distribution do not link against it for licencing reasons.
                             dnl so check if ptscotch can be linked independently or not.
                             AC_CHECK_LIB(pexsi,
                                          [PPEXSIPlanInitialize],
                                          [pexsi_found=yes],
                                          [pexsi_found=no],
                                          [-L${pexsi_prefix}/lib \
                                           -lpexsi \
                                           $SUPERLU_DIST_LDFLAGS \
                                           $SUPERLU_DIST_LIBS \
                                           $PARMETIS_LDFLAGS \
                                           $PARMETIS_LIBS \
                                           $MPI_CLDFLAGS -lm])

                             CPPFLAGS=$__ac_cppflags_
                             CC=$__ac_cc_backup
                             AS_IF([test "$pexsi_found" = "yes"], [
                                         PEXSI_CFLAGS="-I${pexsi_prefix}/include"
                                         PEXSI_LDFLAGS="-L${pexsi_prefix}/lib"
                                         PEXSI_LIBS="-lpexsi"
                                         have_pexsi=yes]);

                             AC_DEFINE([HAVE_PEXSI], 1, ["pexsi is found"])
                         ])
               ])

         AC_SUBST(PEXSI_CFLAGS)
         AC_SUBST(PEXSI_LIBS)
         AC_SUBST(PEXSI_LDFLAGS)

         AS_IF([test "$pexsi_found" == "yes" -a "$pexsi_headers_found" == "yes"],
               [$1],
               [$2])
])
