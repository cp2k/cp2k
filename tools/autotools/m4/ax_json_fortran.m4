# -*- Autoconf -*-
AU_ALIAS([ACX_JSON_FORTRAN], [AX_JSON_FORTRAN])
AC_DEFUN([ACX_JSON_FORTRAN],
        [
        _ac_have_json_fortran=no

        AC_ARG_WITH(json-fortran,
                    [AC_HELP_STRING([--with-json-fortran=yes/no/prefix],
                                   ["use the json fortran library (SIRIUS dependency)"])],
                    [json_fortran_prefix="$withval"],
                    [json_fortran_prefix="/usr/local"])
        pkg_config_path_save=$PKG_CONFIG_PATH
        PKG_CONFIG_PATH=${json_fortran_prefix}/lib/pkgconfig:$PKG_CONFIG_PATH
       AS_IF([test "$with_json_fortran" != "no"],
             [PKG_CHECK_EXISTS([json-fortran], 
                               [PKG_CHECK_MODULES([JSON_FORTRAN], [json-fortran], [_ac_have_json_fortran=yes], [AC_MSG_WARN([json-fortran not found])])], 
                               [])
             ])
        PKG_CONFIG_PATH=$pkg_config_path_save 
        AC_SUBST(JSON_FORTRAN_CFLAGS)
        AC_SUBST(JSON_FORTRAN_LIBS)
        AC_SUBST(JSON_FORTRAN_LDFLAGS)
        AS_IF([test "_ac_have_json_fortran" == "yes"], [$1], [$2])
])
#
