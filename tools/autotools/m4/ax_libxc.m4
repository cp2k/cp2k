# -*- Autoconf -*-

AU_ALIAS([ACX_LIBXC], [AX_LIBXC])
AC_DEFUN([ACX_LIBXC], [
AC_ARG_WITH(libxc,
            [AC_HELP_STRING([--with-libxc=yes/no/prefix], [use libxc for exchange functionals])],
            [libxc_prefix="$withval"],
            [libxc_prefix="/usr/local"])

AS_IF([test "$with_libxc" != "no"],
      [
        __ac_cppflags_=$CPPFLAGS
        CPPFLAGS+=" -I$libxc_prefix/include "

        AC_CHECK_LIB(xc,
                     [xc_func_end],
                     [libxc_found=yes],
                     [libxc_found=no],
                     [-L$libxc_prefix/lib -lm])
       AC_CHECK_HEADER(xc.h,
                      [xc_headers_found=yes],
                      [xc_headers_found=no])

       AC_LANG_PUSH([Fortran])
       AC_CHECK_LIB(xcf90, [xc_f90_version],
                     [xc_f90_found=yes; LIBXC_F90_LIBS="-lxcf90"; LIBXC_F03_LIBS="-lxcf03"],
                     [xc_f90_found=no],
                     [-L${libxc_prefix/lib} -lxc -lm])

       AS_IF([test "$xc_f90_found" = "yes"], [LIBXC_FFLAGS="-I$libxc_prefix/include"])
       AC_LANG_POP([Fortran])

       AS_IF([test "$libxc_found" == "yes" -a "$xc_headers_found" == "yes"],
              [LIBXC_CFLAGS="-I${libxc_prefix}/include"; LIBXC_LIBS="-lxc"; LIBXC_LDFLAGS="-L${libxc_prefix}/lib"
               AC_DEFINE([HAVE_LIBXC], 1, [libxc is found])])
       dnl restore cpp flags
       CPPFLAGS=$__ac_cppflags_
       ]
)

AC_SUBST(LIBXC_LDFLAGS)
AC_SUBST(LIBXC_CFLAGS)
AC_SUBST(LIBXC_LIBS)
AC_SUBST(LIBXC_FFLAGS)
AC_SUBST(LIBXC_F90_LIBS)
AC_SUBST(LIBXC_F03_LIBS)

AS_IF([test "$libxc_found" == "yes" -a "$xc_headers_found" == "yes"], [$1], [$2])
])
