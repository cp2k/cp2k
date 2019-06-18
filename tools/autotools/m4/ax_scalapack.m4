# -*- Autoconf -*-

AU_ALIAS([ACX_SCALAPACK], [AX_SCALAPACK])
AC_DEFUN([ACX_SCALAPACK], [
AC_ARG_WITH(scalapack,
            [AC_HELP_STRING([--with-scalapack=yes/no/prefix], [use scalapack])],
            [scalapack_prefix="$withval"],
            [scalapack_prefix="/usr/local"])

have_scalapack="no"

AS_IF([test "${have_mkl_scalapack}" = "yes"], [AC_DEFINE([HAVE_SCALAPACK], 1, ["scalapack is found"])
                                               have_scalapack="yes"
                                               AC_MSG_NOTICE([Using mkl scalapack])])

AS_IF([test "$with_scalapack" != "no" -a "${have_mkl_scalapack}" = "no"],
      [
      AS_IF([test "$have_C_mpi" == "yes"],
              [
                __ac_cppflags_=$CPPFLAGS

                dnl setup the default directory if necessary
                AS_IF([test "${scalapack_prefix}" == "yes"],
                      [scalapack_prefix="/usr/local"])

                CPPFLAGS+=" -I${scalapack_prefix}/include "
                AC_LANG_PUSH([Fortran])
                __ac_backup_fc=$FC
                FC=${MPIFC}

                AC_CHECK_LIB(scalapack,
                             [pcgetrf],
                             [scalapack_found=yes],
                             [scalapack_found=no],
                             [-L$scalapack_prefix/lib -lm])
               __ac_backup_fc=$FC
               MPIFC=$__ac_backup_fc
               AC_LANG_POP([Fortran])

               AS_IF([test "$scalapack_found" == "yes"],
                     [SCALAPACK_CFLAGS=" -I${scalapack_prefix}/include"
                     SCALAPACK_LIBS="-lscalapack"
                     SCALAPACK_LDFLAGS="-L$scalapack_prefix/lib"])
               dnl restore cpp flags
               CPPFLAGS=$__ac_cppflags_
               ])
       AC_DEFINE([HAVE_SCALAPACK], 1, ["scalapack is found"])
])


AC_SUBST(SCALAPACK_CFLAGS)
AC_SUBST(SCALAPACK_LIBS)
AC_SUBST(SCALAPACK_LDFLAGS)

AS_IF([test "$scalapack_found" == "yes"], [$1], [$2])
])
