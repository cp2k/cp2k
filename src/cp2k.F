!!****** cp2k/cp2k [1.0] *
!!
!!   NAME
!!     cp2k, The main program
!!
!!   COPYRIGHT
!!I---------------------------------------------------------------------------I
!!I                                                                           I
!!I CP2K: A general program to perform molecular dynamics simulations         I
!!I Copyright (C) 2000  CP2K developers group                                 I
!!I                                                                           I
!!I This program is free software; you can redistribute it and/or modify      I
!!I it under the terms of the GNU General Public License as published by      I
!!I the Free Software Foundation; either version 2 of the License, or         I
!!I (at your option) any later version.                                       I
!!I                                                                           I
!!I This program is distributed in the hope that it will be useful,           I
!!I but WITHOUT ANY WARRANTY; without even the implied warranty of            I
!!I MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             I
!!I GNU General Public License for more details.                              I
!!I                                                                           I
!!I You should have received a copy of the GNU General Public License         I
!!I along with this program; if not, write to the Free Software               I
!!I Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!
!!   FUNCTION
!!      Main Program
!!
!!   AUTHOR
!!      CJM & JGH
!!
!!   MODIFICATION HISTORY
!!      none
!!
!!   NOTES
!!     Main Program
!!
!!   SOURCE
!------------------------------------------------------------------------------

PROGRAM cp2k

  USE atom_module,                     ONLY: atomic_code
  USE environment,                     ONLY: cp2k_finalize,&
                                             cp2k_init
  USE fist_main,                       ONLY: fist
  USE global_types,                    ONLY: global_environment_type
  USE kg_main,                         ONLY: kg
  USE library_tests,                   ONLY: lib_test
  USE qs_main,                         ONLY: quickstep
  USE timings,                         ONLY: trace_debug

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: routine = "PROGRAM cp2k"

  TYPE(global_environment_type), TARGET :: globenv
  INTEGER :: handle

! -----------------------------------------------------------------------------

  CALL cp2k_init(globenv)

  CALL trace_debug("START",handle)

  SELECT CASE (globenv%program_name)
  CASE ("ATOM")
    CALL atomic_code(globenv)
  CASE ("FIST")
    CALL fist(globenv)
  CASE ("KG")
    CALL kg(globenv)
  CASE ("QS")
    CALL quickstep(globenv)
  CASE ("TEST")
    CALL lib_test(globenv)
  END SELECT

  CALL cp2k_finalize(globenv)

END PROGRAM cp2k
!!*****
