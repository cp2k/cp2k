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
!!      10.2003 modified to use force_env [fawzi]
!!
!!   NOTES
!!     Main Program
!!
!!   SOURCE
!------------------------------------------------------------------------------
#include "cp_prep_globals.h"

PROGRAM cp2k

  USE atom_module,                     ONLY: atomic_code
  USE environment,                     ONLY: cp2k_finalize,&
       cp2k_init
!  USE fist_main,                       ONLY: fist
  USE global_types,                    ONLY: global_environment_type
  USE kg_main,                         ONLY: kg_create_force_env
  USE library_tests,                   ONLY: lib_test
  USE qs_main,                         ONLY: quickstep_create_force_env
  !  USE qs_ep_methods,                   ONLY: ep_create_force_env
  USE timings,                         ONLY: trace_debug
  USE cp_blacs_env, ONLY: start_blacs, end_blacs
  use cp_para_types, ONLY: cp_para_env_type, cp_blacs_env_type
  USE force_env_types, ONLY: force_env_type, force_env_calculate_energy, force_env_release
  USE cp_log_handling, ONLY: cp_logger_type, cp_failure_level, cp_logger_create, cp_logger_release, cp_to_string
  USE cp_error_handling, ONLY: cp_error_type, cp_a_l, cp_assert, cp_error_init, cp_error_dealloc_ref, cp_assertion_failed
  USE geo_opt, only: qs_geo_opt
  USE md_run, only: qs_mol_dyn 
  USE machine, ONLY: print_memory
  USE qs_mc,                           ONLY: qs_mon_car
  USE trajana_methods, ONLY: trajana

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: routineN = "cp2k",&
       moduleN="cp2k", routineP=moduleN//':'//routineN

  TYPE(global_environment_type), pointer :: globenv
  TYPE(cp_error_type) :: error
  TYPE(force_env_type), pointer :: force_env
  TYPE(cp_blacs_env_type), pointer :: blacs_env
  INTEGER :: handle, stat, i
  logical :: failure

  ! -----------------------------------------------------------------------------

  failure=.false.
  NULLIFY(globenv, force_env, blacs_env)
  ALLOCATE(globenv,stat=stat)
  IF (stat/=0) STOP "cp2k:cp2k, glob_env allocation"

  CALL cp2k_init(globenv)
  call cp_error_init(error)

  ! initialize force_env
  SELECT CASE (globenv%program_name)
  CASE ("ATOM")
     CALL atomic_code(globenv)
  CASE ("TRAJANA")
     CALL trajana(globenv)
  CASE ("FIST")
!     CALL fist(globenv)
  CASE ("KG")
     CALL kg_create_force_env(force_env, globenv, error=error)
  CASE ("QS")
     !   *** Initialize the BLACS environment ***
     CALL start_blacs(blacs_env,globenv)
!FM     DO i=1,100
        CALL quickstep_create_force_env(force_env, blacs_env, globenv, error=error)
!FM        CALL print_memory("qs_created"//cp_to_string(i))
!FM        CALL force_env_release(force_env, error=error)
!FM        CALL print_memory("qs_removed"//cp_to_string(i))
!FM     END DO
  CASE ("EP")
     !   *** Initialize the BLACS environment ***
     CALL start_blacs(blacs_env,globenv)

     !     CALL ep_create_force_env(force_env, blacs_env, globenv, error=error)
  CASE ("TEST")
     CALL lib_test(globenv)
  CASE default
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "Invalid run TYPE <"//TRIM(globenv%run_type)//"> was specified, "//&
          CPSourceFileRef,&
          error,failure)
  END SELECT

  ! perform requested calculation
  IF (globenv%run_type/="NONE") THEN 
     !FM this test should be reduced to the first one

     IF (globenv%program_name=="QS".OR.globenv%program_name=="EP") THEN
       SELECT CASE (globenv%run_type)
       CASE ("GEOMETRY OPTIMIZATION")
          CALL qs_geo_opt(force_env,globenv,error=error)
       CASE ("MONTE CARLO")
         CALL qs_mon_car ( force_env, blacs_env, globenv, error=error )
       CASE ("MOLECULAR DYNAMICS")
          CALL qs_mol_dyn ( force_env, globenv, error=error )
       CASE ("WAVEFUNCTION OPTIMIZATION","ENERGY CALCULATION")
          CALL force_env_calculate_energy(force_env,error=error)
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Invalid run TYPE <"//TRIM(globenv%run_type)//"> was specified, "//&
               CPSourceFileRef,&
               error,failure)
       END SELECT
       CALL force_env_release(force_env, error=error)
#if defined(__HMD)
       CALL _dump_allocated(%VAL(0))
#endif
     ELSE IF(globenv%program_name=="KG") THEN
       SELECT CASE (globenv%run_type)
       CASE ("MOLECULAR DYNAMICS")
          CALL qs_mol_dyn ( force_env, globenv, error=error )
!      CASE ("WAVEFUNCTION OPTIMIZATION","ENERGY CALCULATION")
!         CALL force_env_calculate_energy(force_env,error=error)
!      CASE ("DEBUG")
!         CALL kg_debug(force_env,error=error)
       END SELECT
       CALL force_env_release(force_env, error=error)
     END IF
     ! method specific cleanup
     SELECT CASE (globenv%program_name)
     CASE ("ATOM")
     CASE ("FIST")
     CASE ("KG")
     CASE ("QS")
        !   *** Release BLACS environment ***
        CALL end_blacs(blacs_env,globenv)
     CASE ("EP")
        !   *** Release BLACS environment ***
        CALL end_blacs(blacs_env,globenv)
     CASE ("TEST")
     END SELECT
  END IF

  CALL cp_error_dealloc_ref(error)

  CALL cp2k_finalize(globenv)

END PROGRAM cp2k
!!*****
