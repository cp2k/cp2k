# -*- Autoconf -*-

AU_ALIAS([ACX_LIBXSMM], [AX_LIBXSMM])
AC_DEFUN([ACX_LIBXSMM], [
AC_ARG_WITH(libxsmm,
            [AC_HELP_STRING([--with-libxsmm=yes/no/prefix], [use libsmm library (small matrix multiplication)])],
            [libxsmm_prefix="$withval"],
            [libxsmm_prefix="/usr/local"])

AS_IF([test "$with_libxsmm" != "no"],
      [
       AS_IF([test "$libxsmm_prefix" = "yes"],
             [libxsmm_prefix="/usr/local"])
       __ac_cppflags_=$CPPFLAGS
       CPPFLAGS+="-I$libxsmm_prefix/include"

       AC_CHECK_LIB(xsmm,
                    [libxsmm_blas_xgemm],
                    [xsmm_found=yes],
                    [xsmm_found=no],
                    [-L$libxsmm_prefix/lib $BLAS_LIBS -lpthread -lrt -lm -ldl])
       AC_CHECK_LIB(xsmmext,
                    [libxsmm_mmbatch_omp],
                    [xsmmext_found=yes],
                    [xsmmext_found=no],
                    [-L$libxsmm_prefix/lib -lxsmm $BLAS_LIBS $OMP_LIBS -lpthread -lrt -lm -ldl])

       AS_IF([test "$xsmm_found" = "yes"],
             [LIBXSMM_CFLAGS="-I${libxsmm_prefix}/include"
             LIBXSMM_LIBS="-lxsmm"
             LIBXSMM_FORTRAN_LIBS="-lxsmmf"
             LIBXSMM_LDFLAGS="-L${libxsmm_prefix}/lib"
             AC_DEFINE([HAVE_LIBXSMM], 1, ["libxsmm is found"])])

      AS_IF([test "$xsmmext_found" = "yes"],
            [LIBXSMMEXT_CFLAGS="-I${libxsmm_prefix}/include"
            LIBXSMMEXT_LIBS="-lxsmmext"])
      ])


AC_SUBST(LIBXSMM_CFLAGS)
AC_SUBST(LIBXSMM_LIBS)
AC_SUBST(LIBXSMM_LDFLAGS)
AC_SUBST(LIBXSMM_FORTRAN_LIBS)
AC_SUBST(LIBXSMMEXT_LIBS)
AC_SUBST(LIBXSMM_LDFLAGS)

AS_IF([test "$xsmm_found" == "yes"], [$1], [$2])
])
