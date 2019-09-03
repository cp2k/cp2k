# -*- Autoconf -*-

AU_ALIAS([ACX_SUPERLU], [AX_SUPERLU])
AC_DEFUN([ACX_SUPERLU], [

AC_REQUIRE([LX_FIND_MPI()])
AC_REQUIRE([ACX_PARMETIS])

dnl it is not exported outside the configure script
have_superlu=no

AC_ARG_WITH(superlu,
            [AC_HELP_STRING([--with-superlu=yes/no/prefix], [use the superlu library for LU decomposition (default yes)])],
            [libsuperlu_prefix="$withval"],
            [libsuperlu_prefix="/usr/local"])

            AS_IF([test "$with_libsuperlu" != "no" -a $have_C_mpi == "yes"],
                  [
                   AS_IF([test "$libsuperlu_prefix" == "yes"], [libsuperlu_prefix="/usr/local"])

                   __ac_cppflags_=$CPPFLAGS
                   __ac_cc_compiler=$CC
                   CC=${MPICC}
                   CPPFLAGS+="-I${libsuperlu_prefix}/include -I${libsuperlu_prefix}/include/superlu_dist ${MPI_CFLAGS}"

                   AC_LANG_PUSH([C])
                   AC_CHECK_HEADERS([superlu_ddefs.h],
                                    [superlu_has_hdrs="yes"],
                                    [superlu_has_hdrs="no"],
                                    [])
                   AC_CHECK_HEADERS([superlu_dist/superlu_ddefs.h],
                                    [superlu_has_hdrs="yes"],
                                    [superlu_has_hdrs="no"],
                                    [])

                   dnl check if superlu_dist needs parmetis. Fedora for instance
                   dnl does not provide parmetis and compiles superlu_dist
                   dnl without. So we check without parmetis and if it does not
                   dnl work check again with

                   AC_CHECK_LIB(superlu_dist,
                                [dCreate_CompCol_Matrix_dist],
                                [superlu_found=yes],
                                [superlu_found=no],
                                [-L${libsuperlu_prefix}/lib -lsuperlu_dist $MPI_CLDFLAGS -lpthread -lrt -lm])
                   AC_MSG_CHECKING([whether superlu_dist can be linked without parmetis])

                   AS_IF([test "$superlu_found" = "yes"],
                         [AC_MSG_RESULT(yes)],
                         [AC_MSG_RESULT(no)
                          AC_CHECK_LIB(superlu_dist,
                                       [dCreate_CompCol_Matrix_dist],
                                       [superlu_found=yes],
                                       [superlu_found=no],
                                       [-L${libsuperlu_prefix}/lib -lsuperlu_dist $PARMETIS_LDFLAGS $PARMETIS_LIBS $MPI_CLDFLAGS -lpthread -lrt -lm])
                         ])

                   AC_LANG_POP([C])

                   CC=$__ac_cc_compiler
                   CPPFLAGS=$__ac_cppflags_
                   AS_IF([test "$superlu_found" == "yes"],
                         [LIBSUPERLU_DIST_CFLAGS="-I${libsuperlu_prefix}/include"
                          LIBSUPERLU_DIST_LIBS="-lsuperlu_dist"
                          LIBSUPERLU_DIST_LDFLAGS="-L${libsuperlu_prefix}/lib"
                          have_superlu=yes
                          AC_DEFINE([HAVE_SUPERLU_DIST], 1, ["superlu is found"])
                          ], [have_superlu=no])
                ])

          AC_SUBST(SUPERLU_DIST_CFLAGS)
          AC_SUBST(SUPERLU_DIST_LIBS)
          AC_SUBST(SUPERLU_DIST_LDFLAGS)

          AS_IF([test "$superlu_found" == "yes"], [$1], [$2])
])
