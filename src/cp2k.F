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

  USE atom_module,ONLY : atomic_code
  USE control_module, ONLY : control
  USE cp2k_input, ONLY : read_cp2k_section
  USE environment, ONLY : initialisation, trailer
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl, print_kind_info
  USE library_tests,ONLY : lib_test
  USE parallel, ONLY : start_parallel, end_parallel
  USE physcon, ONLY : write_physcon
  USE tbmd_module, ONLY : tbmd
  USE wave_module, ONLY : wave
  USE fit_module, ONLY : fit
  USE qs_main, ONLY : quickstep
  USE termination,  ONLY: stop_program
  USE timings, ONLY: trace_debug

  IMPLICIT NONE

  TYPE ( global_environment_type ) :: globenv

!------------------------------------------------------------------------------
  print*
  !CALL trace_debug("START")
  CALL start_parallel ( globenv )

  CALL initialisation ( globenv )

  CALL read_cp2k_section ( globenv )

  IF (globenv%ionode.AND.globenv%print%basic_data_types) THEN
    CALL print_kind_info(globenv%scr)
  END IF
  IF (globenv%ionode.AND.globenv%print%physcon) THEN
    CALL write_physcon(globenv%scr)
  END IF

  SELECT CASE ( globenv % program_name )

  CASE DEFAULT

     CALL stop_program ( "cp2k", "No valid program name given "&
                         //globenv % program_name )

  CASE ( "FIST", "POL" )

     CALL control ( globenv )

  CASE ( "FIT" )

     CALL fit ( globenv )

  CASE ( "WAVE" )

     CALL wave ( globenv )

  CASE ( "TBMD" )

     CALL tbmd ( globenv )

  CASE ( "QUICKSTEP", "QS" )

     CALL quickstep ( globenv )

  CASE ( "ATOM" )

     CALL atomic_code ( globenv )

  CASE ( "TEST" )

     CALL lib_test ( globenv )

  END SELECT

  CALL end_parallel ( )

  CALL trailer ( globenv )

!------------------------------------------------------------------------------

END PROGRAM cp2k

!------------------------------------------------------------------------------
!!*****
