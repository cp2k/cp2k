# -*- Autoconf -*-

AU_ALIAS([ACX_SPGLIB], [AX_SPGLIB])
AC_DEFUN([ACX_SPGLIB],
         [
         have_spglib=no
         AC_ARG_WITH(spglib,
                     [AC_HELP_STRING([--with-spglib=yes/no/PATH],
                                      [check for spglib support (default yes)])],
                                      [spglib_prefix=$withval],
                                      [spglib_prefix=/usr/local])

         AS_IF([test "$with_spglib" != "no"],
               [
               __ac_cppflags_=$CPPFLAGS
               CPPFLAGS="-I$spglib_prefix/include"
               #check for the header file
               AC_CHECK_HEADER([spglib.h],
                               [SPGLIB_CFLAGS=-I$spglib_prefix/include
                                spglib_header_found=yes],
                                [AC_CHECK_HEADER([spglib/spglib.h],
                                                 [SPGLIB_CFLAGS="-I$spglib_prefix/include -I$spglib_prefix/include/spglib"
                                                  spglib_header_found=yes],
                                                 [spglib_header_found=no])
                                ])
              CPPFLAGS=$__ac_cppflags_
              AS_IF([test "$spglib_header_found" == "yes"], [
              # strickly speaking we should also check for the library, but if the header file
              # is not present no point for checking the library
              SPGLIB_LIBS="-lsymspg"
              SPGLIB_LDFLAGS="-L$spglib_prefix/lib"
              have_spglib="yes"])
       ])
       AC_SUBST(SPGLIB_CFLAGS)
       AC_SUBST(SPGLIB_LDFLAGS)
       AC_SUBST(SPGLIB_LIBS)

       AS_IF([test "$spglib_header_found" = "yes"],
             [AC_DEFINE(HAVE_SPGLIB, 1, [spglib is present])])
       AS_IF([test "$have_spglib" = "yes"], [$1], [$2])
      ])
