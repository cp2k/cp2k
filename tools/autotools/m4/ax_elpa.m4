# -*- Autoconf -*-
AU_ALIAS([ACX_ELPA], [AX_ELPA])
AC_DEFUN([ACX_ELPA],
         [
         have_elpa=no
         AC_ARG_WITH(elpa,
                     [AC_HELP_STRING([--with-elpa=yes/no/prefix],
                                     ["use the elpa library for diagonalization"])],
                                     [elpa_prefix="$withval"],
                                     [elpa_prefix="/usr"])
        AC_ARG_WITH(elpa-version,
                    [AC_HELP_STRING([--with-elpa-version=2017.05.003],
                                     ["indicate the library version. Might be useful if the script does not detect it"])],
                                     [elpa_version="$withval"],
                                     [elpa_version="2017.05.003"])

         AS_IF([test "$with_elpa" != "no"],
               [
               AS_IF([test "${elpa_prefix}" == "yes"], [elpa_prefix="/usr"])
               dnl check for pc files first
               pushd "$elpa_prefix/lib64/pkgconfig" >/dev/null
               elpa_pc=`ls elpa-*.pc 2>/dev/null`
               popd >/dev/null
               AS_IF([test "${elpa_pc}yes" != "yes"], [elpa_version=${elpa_pc:5:11}
               AC_MSG_CHECKING(which elpa version is found)
               AC_MSG_RESULT($elpa_version)
               ])

               AS_IF([test "${PKG_CONFIG}yes" != "yes"],
                     [PKG_CHECK_EXISTS([elpa-${elpa_version}],
                                       [PKG_CHECK_MODULES(ELPA,
                                       [elpa-${elpa_version}],
                                       [elpa_found="yes"
                                       elpa_has_hdrs="yes"
                                       elpa_has_libs="yes"
                                       have_elpa=yes],
                                       [elpa_found="no"])],
                                       [elpa_found="no"])
                      PKG_CHECK_EXISTS([elpa_openmp-${elpa_version}],
                                       [PKG_CHECK_MODULES(ELPA_OPENMP,
                                                         [elpa_openmp-${elpa_version}],
                                                         [elpa_openmp_found="yes"
                                                          elpa_openmp_has_hdrs="yes"
                                                          elpa_openmp_has_libs="yes"
                                                          have_elpa_openmp=yes],
                                                         [elpa_openmp_found="no"])],
                                       [elpa_openmp_found="no"])
                     ])
              ])
        AC_SUBST(ELPA_CFLAGS)
        AC_SUBST(ELPA_LDFLAGS)
        AC_SUBST(ELPA_LIBS)
        AC_SUBST(ELPA_OPENMP_CFLAGS)
        AC_SUBST(ELPA_OPENMP_LDFLAGS)
        AC_SUBST(ELPA_OPENMP_LIBS)
        AS_IF([test "$elpa_found" == "yes"], [$1], [$2])
])
