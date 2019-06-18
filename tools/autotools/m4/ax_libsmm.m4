# -*- Autoconf -*-

AU_ALIAS([ACX_LIBSMM], [AX_LIBSMM])
AC_DEFUN([ACX_LIBSMM], [
AC_ARG_WITH(libsmm,
            [AC_HELP_STRING([--with-libsmm=yes/no/prefix], [use libsmm library (small matrix multiplication)])],
            [libsmm_prefix="$withval"],
            [libsmm_prefix="/usr/local"])

AS_IF([test "$with_libsmm" != "no"],
      [
        AS_IF([test "$libsmm_prefix" == "yes"],
              [libsmm_prefix="/usr/local"])
        __ac_cppflags_=$CPPFLAGS
        CPPFLAGS+="-I$libsmm_prefix/include"

        AC_CHECK_LIB(smm_dnn,
                     [xc_func_end],
                     [smm_found=yes],
                     [smm_found=no],
                     [-L$libsnn_prefix/lib -lm])

       AS_IF([test "$libsmm_found" == "yes"],
             [LIBSMM_CFLAGS="-I${libxc_prefix}/include"
             LIBSMM_FFLAGS="-I${libxc_prefix}/include"
             LIBSMM_LIBS="-I${libxc_prefix}/lib -lsmm_dnn"])
             AC_DEFINE([HAVE_LIBSMM], 1, ["libsmm is found"])
             ]
)

AC_SUBST(LIBSMM_FFLAGS)
AC_SUBST(LIBSMM_CFLAGS)
AC_SUBST(LIBSMM_LIBS)
AC_SUBST(LIBSMM_LDFLAGS)

AS_IF([test "$libsmm_found" == "yes"], [$1], [$2])
])
