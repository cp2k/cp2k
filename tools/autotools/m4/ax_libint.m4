# -*- Autoconf -*-

AU_ALIAS([ACX_LIBINT], [AX_LIBINT])
AC_DEFUN([ACX_LIBINT], [
AC_ARG_WITH(libint2,
            [AC_HELP_STRING([--with-libint=yes/no/prefix], [use libint2 for calculating gaussian integrals])],
            [libint_prefix="$withval"],
            [libint_prefix="/usr/local"])

AS_IF([test "$with_libint" != "no"],
      [
        ac_pkg_config_backup=$PKG_CONFIG_PATH
        PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$fftw3_prefix/lib/pkgconfig
        AS_IF([test "${PKG_CONFIG}yes" != "yes"],
              [PKG_CHECK_EXISTS([libint2],
                                [PKG_CHECK_MODULES(LIBINT2,
                                                   [libint2 >= 2.1.0],
                                                   [pc_libint2_found="yes"
                                                   libint2_has_hdrs="yes"
                                                   libint2_has_libs="yes"
                                                   libint2_found="yes"],
                                                   [pc_libint2_found="no"])],
                                                   [pc_libint2_found="no"])])
        PKG_CONFIG_PATH=$ac_pkg_config_backup

        dnl old fashion way if pkg-config is not found
        AS_IF([test "$pc_libint2_found" == "no"],
               [__ac_cppflags_=$CPPFLAGS
                AS_IF([test "${libint_prefix}" == "yes"],
                      [libint_prefix="/usr/local"])

                CPPFLAGS+=" -I${libint_prefix}/include "

                AC_CHECK_LIB(int2,
                            [libint2_init_default],
                            [libint2_found=yes],
                            [libint2_found=no],
                            [-L$libint_prefix/lib $libint_prefix/lib64 -lint2 -lm])
                AC_CHECK_HEADER(libint.h,
                                [libint2_headers_found=yes],
                                [libint2_headers_found=no])

                AS_IF([test "$libint2_found" == "yes" -a "$libint2_headers_found" == "yes"],
                      [LIBINT_CFLAGS=" -I${libint_prefix}/include -I${libint_prefix}/include/libint2"
                      LIBINT_LIBS="-L${libint_prefix}/lib64 $libint_prefix/lib -lint2"
                      AC_DEFINE([HAVE_LIBINT2], 1, ["libint2 is found"])
                      libint2_found="yes"])
                      dnl restore cpp flags
                CPPFLAGS=$__ac_cppflags_
                ],
                [libint2_found="yes"])
       ])

       AC_SUBST(LIBINT2_CFLAGS)
       AC_SUBST(LIBINT2_LIBS)
       AC_SUBST(LIBINT2_FFLAGS)
       AS_IF([test "$libint2_found" == "yes"], [$1], [$2])
])
