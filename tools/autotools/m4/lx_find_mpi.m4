# -*- Autoconf -*-

#################################################################################################
# Copyright (c) 2010, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory
# Written by Todd Gamblin, tgamblin@llnl.gov.
# LLNL-CODE-417602
# All rights reserved.
#
# This file is part of Libra. For details, see http://github.com/tgamblin/libra.
# Please also read the LICENSE file for further information.
#
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice, this list of
#    conditions and the disclaimer below.
#  * Redistributions in binary form must reproduce the above copyright notice, this list of
#    conditions and the disclaimer (as noted below) in the documentation and/or other materials
#    provided with the distribution.
#  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse
#    or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################################

#
# LX_FIND_MPI()
#  ------------------------------------------------------------------------
# This macro finds a MPI compiler and extracts includes and libraries from
# it for use in automake projects.  The script exports the following variables:
#
# AC_DEFINE variables:
#     HAVE_MPI         AC_DEFINE'd to 1 if we found MPI
#
# AC_SUBST variables:
#     MPICC            Name of MPI compiler
#     MPI_CFLAGS       Includes and defines for MPI C compilation
#     MPI_CLDFLAGS     Libraries and library paths for linking MPI C programs
#
#     MPICXX           Name of MPI C++ compiler
#     MPI_CXXFLAGS     Includes and defines for MPI C++ compilation
#     MPI_CXXLDFLAGS   Libraries and library paths for linking MPI C++ programs
#
#     MPIF77           Name of MPI Fortran 77 compiler
#     MPI_F77FLAGS     Includes and defines for MPI Fortran 77 compilation
#     MPI_F77LDFLAGS   Libraries and library paths for linking MPI Fortran 77 programs
#
#     MPIFC            Name of MPI Fortran compiler
#     MPI_FFLAGS       Includes and defines for MPI Fortran compilation
#     MPI_FLDFLAGS     Libraries and library paths for linking MPI Fortran programs
#
# Shell variables output by this macro:
#     have_C_mpi       'yes' if we found MPI for C, 'no' otherwise
#     have_CXX_mpi     'yes' if we found MPI for C++, 'no' otherwise
#     have_F77_mpi     'yes' if we found MPI for F77, 'no' otherwise
#     have_F_mpi       'yes' if we found MPI for Fortran, 'no' otherwise
#
AC_DEFUN([LX_FIND_MPI],
[
     AC_LANG_CASE(
     [C], [
         AC_REQUIRE([AC_PROG_CC])
         AS_IF([test "${MPICC}yes" != "yes"], [LX_QUERY_MPI_COMPILER(MPICC, [$MPICC], C)], [LX_QUERY_MPI_COMPILER(MPICC, [mpicc mpiicc mpixlc mpipgcc], C)])
          ],
     [C++], [
         AC_REQUIRE([AC_PROG_CXX])
         if [[ ! -z "$MPICXX" ]]; then
             LX_QUERY_MPI_COMPILER(MPICXX, [$MPICXX], CXX)
         else
             LX_QUERY_MPI_COMPILER(MPICXX, [mpicxx mpiCC mpic++ mpig++ mpiicpc mpipgCC mpixlC], CXX)
         fi
     ],
     [F77], [
         AC_REQUIRE([AC_PROG_F77])
         if [[ ! -z "$MPIF77" ]]; then
             LX_QUERY_MPI_COMPILER(MPIF77, [$MPIF77], F77)
         else
             LX_QUERY_MPI_COMPILER(MPIF77, [mpif77 mpiifort mpixlf77 mpixlf77_r], F77)
         fi
     ],
     [Fortran], [
         AC_REQUIRE([AC_PROG_FC])
         if [[ ! -z "$MPIFC" ]]; then
             LX_QUERY_MPI_COMPILER(MPIFC, [$MPIFC], F)
         else
             mpi_default_fc="mpif95 mpif90 mpigfortran mpif2003"
             mpi_intel_fc="mpiifort"
             mpi_xl_fc="mpixlf95 mpixlf95_r mpixlf90 mpixlf90_r mpixlf2003 mpixlf2003_r"
             mpi_pg_fc="mpipgf95 mpipgf90"
             LX_QUERY_MPI_COMPILER(MPIFC, [$mpi_default_fc $mpi_intel_fc $mpi_xl_fc $mpi_pg_fc], F)
         fi
     ])

     AS_IF([test "$have_C_mpi" == "yes" -o "$have_CXX_mpi" == "yes" -o "$have_F77_mpi" == "yes" -o "$have_F_mpi" == "yes"], [$1], [$2])
])


#
# LX_QUERY_MPI_COMPILER([compiler-var-name], [compiler-names], [output-var-prefix])
#  ------------------------------------------------------------------------
# AC_SUBST variables:
#     MPI_<prefix>FLAGS       Includes and defines for MPI compilation
#     MPI_<prefix>LDFLAGS     Libraries and library paths for linking MPI C programs
#
# Shell variables output by this macro:
#     found_mpi_flags         'yes' if we were able to get flags, 'no' otherwise
#
AC_DEFUN([LX_QUERY_MPI_COMPILER],
[
     # Try to find a working MPI compiler from the supplied names
     AC_PATH_PROGS($1, [$2], [not-found])

     # Figure out what the compiler responds to to get it to show us the compile
     # and link lines.  After this part of the macro, we'll have a valid
     # lx_mpi_command_line
     echo -n "Checking whether $$1 responds to '-showme:compile'... "
     lx_mpi_compile_line=`$$1 -showme:compile 2>/dev/null`
     if [[ "$?" -eq 0 ]]; then
         echo yes
         lx_mpi_link_line=`$$1 -showme:link 2>/dev/null`
     else
         echo no
         echo -n "Checking whether $$1 responds to '-showme'... "
         lx_mpi_command_line=`$$1 -showme 2>/dev/null`
         if [[ "$?" -ne 0 ]]; then
             echo no
             echo -n "Checking whether $$1 responds to '-compile-info'... "
             lx_mpi_compile_line=`$$1 -compile-info 2>/dev/null`
             if [[ "$?" -eq 0 ]]; then
                 echo yes
                 lx_mpi_link_line=`$$1 -link-info 2>/dev/null`
             else
                 echo no
                 echo -n "Checking whether $$1 responds to '-show'... "
                 lx_mpi_command_line=`$$1 -show 2>/dev/null`
                 if [[ "$?" -eq 0 ]]; then
                     echo yes
                 else
                     echo no
                 fi
             fi
         else
             echo yes
         fi
     fi

     if [[ ! -z "$lx_mpi_compile_line" -a ! -z "$lx_mpi_link_line" ]]; then
         lx_mpi_command_line="$lx_mpi_compile_line $lx_mpi_link_line"
     fi

     if [[ ! -z "$lx_mpi_command_line" ]]; then
         # Now extract the different parts of the MPI command line.  Do these separately in case we need to
         # parse them all out in future versions of this macro.
         lx_mpi_defines=`    echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-D\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         lx_mpi_includes=`   echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-I\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         lx_mpi_link_paths=` echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-L\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         lx_mpi_libs=`       echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-l\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`
         lx_mpi_link_args=`  echo "$lx_mpi_command_line" | grep -o -- '\(^\| \)-Wl,\([[^\"[:space:]]]\+\|\"[[^\"[:space:]]]\+\"\)'`

         # Create variables and clean up newlines and multiple spaces
         MPI_$3FLAGS="$lx_mpi_defines $lx_mpi_includes"
         MPI_$3LDFLAGS="$lx_mpi_link_paths $lx_mpi_libs $lx_mpi_link_args"
         MPI_$3FLAGS=`  echo "$MPI_$3FLAGS"   | tr '\n' ' ' | sed 's/^[[ \t]]*//;s/[[ \t]]*$//' | sed 's/  +/ /g'`
         MPI_$3LDFLAGS=`echo "$MPI_$3LDFLAGS" | tr '\n' ' ' | sed 's/^[[ \t]]*//;s/[[ \t]]*$//' | sed 's/  +/ /g'`

         OLD_CPPFLAGS=$CPPFLAGS
         OLD_LIBS=$LIBS
         CPPFLAGS=$MPI_$3FLAGS
         LIBS=$MPI_$3LDFLAGS

         AS_IF([test "$3" = "C" -o "$3" = "CXX"], 
               [
                AC_TRY_LINK([#include <mpi.h>],
                             [int rank, size;
                              MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                              MPI_Comm_size(MPI_COMM_WORLD, &size);],
                              [# Add a define for testing at compile time.
                              AC_DEFINE([HAVE_MPI], [1], [Define to 1 if you have MPI libs and headers.])
                              have_$3_mpi='yes'],
                              [# zero out mpi flags so we don't link against the faulty library.
                               MPI_$3FLAGS=""
                               MPI_$3LDFLAGS=""
                               have_$3_mpi='no'])], 
                             [AC_DEFINE([HAVE_MPI], [1], [Define to 1 if you have MPI libs and headers.])
                              have_$3_mpi='yes'])

         # AC_SUBST everything.
         AC_SUBST($1)
         AC_SUBST(MPI_$3FLAGS)
         AC_SUBST(MPI_$3LDFLAGS)

         LIBS=$OLD_LIBS
         CPPFLAGS=$OLD_CPPFLAGS
     else
         echo Unable to find suitable MPI Compiler. Try setting $1.
         have_$3_mpi='no'
     fi
])
