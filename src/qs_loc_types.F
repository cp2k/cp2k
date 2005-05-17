!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_loc_types [1.0] *
!!
!!   NAME
!!     qs_loc_types
!!
!!   FUNCTION
!!     New version of the module for the localization of the molecular orbitals
!!     This should be able to use different definition of the spread functional
!!     It should also calculate the integrals analytically so that it can be
!!     used irrespective of the pw_env and the collocation of wfn on the grids
!!     It should also work with a selected set of states, instead than all of them,
!!     in this case one should check that the selected states have the same occupation number
!!     The spread functional can be only estimated, or also optimized by minimization
!!     and in principle also maximization should be available.
!!     This operations can be required irrespective of the printing requirements
!!     It would be highly desirable to do all this along a MD run every N steps,
!!     and have a trajectory of the centeroids of the localized wfn
!!     In addition these functions can be used for properties calculations 
!!     like NMR and XAS. Therefore it is necessary that the rotated wfn are then copied
!!     in the mos fm matrix to be available for next use.
!!
!!   AUTHORS
!!     MI (05-2005)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_loc_types

! *****************************************************************************

  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pool_create_fm
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE particle_types,                  ONLY: particle_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_release,&
                                             localized_wfn_control_retain,&
                                             localized_wfn_control_type
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: module_name = "qs_loc_types"


! *** Public ***
  PUBLIC :: qs_loc_env_create, qs_loc_env_destroy, &
            qs_loc_env_release, qs_loc_env_retain, &
            get_qs_loc_env, set_qs_loc_env
  PUBLIC :: qs_loc_env_new_type

!****t* qs_loc_types/qs_loc_env_new_type [1.0] *
!!
!!   NAME
!!     qs_loc_env_new_type
!!
!!   FUNCTION
!!     contains all the info needed by quickstep to calculate
!!     the spread of a selectedset of orbitals and if required
!!     to minimize or mazimize the spread by rotation of the orbitals
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     is_initialized : logical control variable
!!     para_env: info for the distribution of the calculations
!!     mo_coeff: full matrix containing only the selected subset of orbitals
!!     local_molecules: molecules distributed 
!!     cell : box that contains the system
!!     localized_wfn_control: variables and parameter that define the spread 
!!                            functional and the optimization algorithm
!!     particle_set: position, type, ao_indexes etc for each atom
!!     op_sm_set: set of sparse matrices used to define the spread operator
!!                 when the functional is defined by the use operator acting on the
!!                 basis functions, e.g. the Berry phase definition
!!                 The matrix element of the type <a|O|b> are computed in initialization
!!                 of qs_loc_env
!!     op_fm_set: set of full matrices used to define the spread operator
!!                when the functional has to be defined directly using the products of MOS
!!                as in the case of the Pipek-Mezek definition.
!!     weights: for a spread defined as extension of the orbitral in the box, these
!!              factors renormilize with respect to the box size 
!!get_qs_loc_env
!!   NOTE
!!     this type should replace the previous set up for the localization of the wfn
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     04-05 created
!!
!!*** *********************************************************************


  TYPE qs_loc_env_new_type

       INTEGER :: is_initialized
       INTEGER :: ref_count
       LOGICAL :: molecular_states
       TYPE ( cp_para_env_type ), POINTER       :: para_env
       TYPE ( cp_fm_p_type), DIMENSION(:),  &
         POINTER                                :: moloc_coeff
       TYPE ( cp_fm_p_type), DIMENSION(:,:),&
         POINTER                                :: op_fm_set
       TYPE(distribution_1d_type), POINTER      :: local_molecules
       TYPE ( cell_type ), POINTER              :: cell
       TYPE ( localized_wfn_control_type ), &
         POINTER                                :: localized_wfn_control
       TYPE (particle_type), DIMENSION(:),&
         POINTER                                :: particle_set 
       TYPE (real_matrix_p_type), DIMENSION(:,:),&
         POINTER                                :: op_sm_set

       REAL (KIND = dp)                         :: weights ( 6 )
  END TYPE qs_loc_env_new_type
  

! *****************************************************************************

CONTAINS

!******************************************************************************
!!****f* qs_loc_types/qs_loc_env_create[1.0] *

  SUBROUTINE qs_loc_env_create(qs_loc_env,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                 :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'qs_loc_env_create', &
      routineP = module_name//':'//routine_name

    INTEGER                                    :: istat
    LOGICAL                                    :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(qs_loc_env),cp_failure_level,routineP,error,failure)

    IF(.NOT. failure) THEN
       ALLOCATE(qs_loc_env,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       qs_loc_env%ref_count = 1

       NULLIFY(qs_loc_env%para_env)
       NULLIFY(qs_loc_env%cell)
       NULLIFY(qs_loc_env%op_sm_set)
       NULLIFY(qs_loc_env%op_fm_set)
       NULLIFY(qs_loc_env%local_molecules)
       NULLIFY(qs_loc_env%moloc_coeff)
       NULLIFY(qs_loc_env%particle_set)
       NULLIFY(qs_loc_env%localized_wfn_control)
       qs_loc_env%weights = 0.0_dp

    END IF

  END SUBROUTINE qs_loc_env_create

!******************************************************************************
!!****f* qs_loc_types/qs_loc_env_destroy [1.0] *

  SUBROUTINE qs_loc_env_destroy(qs_loc_env,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'qs_loc_env_destroy', &
      routineP = module_name//':'//routine_name

    INTEGER                                    :: i, istat, j
    LOGICAL                                    :: failure

    failure =.FALSE.
    CPPrecondition(ASSOCIATED(qs_loc_env),cp_warning_level,routineP,error,failure)

    IF(.NOT. failure) THEN
       IF(ASSOCIATED(qs_loc_env%cell)) CALL cell_release(qs_loc_env%cell,error=error)
       IF(ASSOCIATED(qs_loc_env%local_molecules)) &
          CALL distribution_1d_release(qs_loc_env%local_molecules,error=error)
       IF (ASSOCIATED(qs_loc_env%localized_wfn_control)) THEN
         CALL localized_wfn_control_release(qs_loc_env%localized_wfn_control,&
                                            error=error)
       END IF
       IF(ASSOCIATED(qs_loc_env%para_env)) CALL cp_para_env_release(qs_loc_env%para_env)
       IF(ASSOCIATED(qs_loc_env%particle_set))  NULLIFY(qs_loc_env%particle_set)

       IF(ASSOCIATED(qs_loc_env%moloc_coeff)) THEN
         DO i=1,SIZE ( qs_loc_env % moloc_coeff,1) 
           CALL cp_fm_release(qs_loc_env%moloc_coeff(i)%matrix,error=error)
         END DO
         DEALLOCATE(qs_loc_env%moloc_coeff,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF(ASSOCIATED(qs_loc_env%op_fm_set)) THEN
         DO i=1,SIZE ( qs_loc_env % op_fm_set,2) 
           DO j=1,SIZE ( qs_loc_env % op_fm_set,1)
             CALL cp_fm_release(qs_loc_env%op_fm_set(j,i)%matrix,error=error)
           END DO
         END DO
         DEALLOCATE(qs_loc_env%op_fm_set,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF(ASSOCIATED(qs_loc_env%op_sm_set)) THEN
         DO i=1,SIZE ( qs_loc_env % op_sm_set, 2 )
           DO j=1,SIZE ( qs_loc_env % op_sm_set, 1 )
             CALL deallocate_matrix(qs_loc_env%op_sm_set(j,i)%matrix)
           ENDDO
         END DO
         DEALLOCATE(qs_loc_env%op_sm_set,STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
       
    END IF 

  END SUBROUTINE qs_loc_env_destroy

  SUBROUTINE qs_loc_env_release(qs_loc_env,error)

    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routine_name='qs_loc_env_release',&
      routineP = module_name//':'//routine_name

    INTEGER :: istat
    LOGICAL :: failure
    failure=.FALSE.

    IF (ASSOCIATED(qs_loc_env)) THEN
      CPPreconditionNoFail(qs_loc_env%ref_count>0,cp_failure_level,routineP,error)
      qs_loc_env%ref_count = qs_loc_env%ref_count -1
      IF (qs_loc_env%ref_count==0) THEN
        CALL qs_loc_env_destroy(qs_loc_env,error)
      END IF 
    END IF
  END SUBROUTINE qs_loc_env_release

  SUBROUTINE qs_loc_env_retain(qs_loc_env,error)

    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routine_name='qs_loc_env_retain',&
      routineP = module_name//':'//routine_name

    INTEGER :: istat
    LOGICAL :: failure
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(qs_loc_env),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
      CPPreconditionNoFail(qs_loc_env%ref_count>0,cp_failure_level,routineP,error)
      qs_loc_env%ref_count = qs_loc_env%ref_count +1
    END IF
  END SUBROUTINE qs_loc_env_retain

!******************************************************************************
!!****f* qs_loc_types/get_qs_loc_env [1.0] *

  SUBROUTINE get_qs_loc_env(qs_loc_env,cell,local_molecules,&
             localized_wfn_control,moloc_coeff,&
             op_sm_set,op_fm_set,para_env,particle_set,weights,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(cell_type), OPTIONAL, POINTER         :: cell
    TYPE(cp_fm_p_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: moloc_coeff
    TYPE(cp_fm_p_type), DIMENSION(:,:),&
      POINTER, OPTIONAL                        :: op_fm_set
   TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                  :: local_molecules
    TYPE(localized_wfn_control_type), &
      POINTER, OPTIONAL                        :: localized_wfn_control
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      OPTIONAL, POINTER                        :: op_sm_set
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                  :: para_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: particle_set
    REAL(dp), DIMENSION(6), OPTIONAL           :: weights
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error

    CHARACTER(len=*), PARAMETER :: routine_name = 'get_qs_loc_env', &
      routineP = module_name//':'//routine_name

    INTEGER                                    :: istat
    LOGICAL                                    :: failure

     failure =.FALSE.
     CPPrecondition(ASSOCIATED(qs_loc_env),cp_failure_level,routineP,error,failure)
 
     IF(.NOT.failure) THEN
       IF (PRESENT(cell)) cell =>  qs_loc_env%cell
       IF (PRESENT(moloc_coeff)) moloc_coeff => qs_loc_env%moloc_coeff
       IF (PRESENT(local_molecules)) local_molecules => qs_loc_env%local_molecules
       IF (PRESENT(localized_wfn_control)) &
          localized_wfn_control => qs_loc_env%localized_wfn_control
       IF (PRESENT(op_sm_set)) op_sm_set =>  qs_loc_env%op_sm_set
       IF (PRESENT(op_fm_set)) op_fm_set =>  qs_loc_env%op_fm_set
       IF (PRESENT(para_env)) para_env =>  qs_loc_env%para_env
       IF (PRESENT(particle_set)) particle_set =>  qs_loc_env%particle_set
       IF (PRESENT(weights)) weights(1:6)=  qs_loc_env%weights(1:6)
     END IF
 
  END SUBROUTINE get_qs_loc_env

!******************************************************************************
!!****f* qs_loc_types/set_qs_loc_env [1.0] *

  SUBROUTINE set_qs_loc_env(qs_loc_env,cell,local_molecules,&
             localized_wfn_control,moloc_coeff,&
             op_sm_set,op_fm_set,para_env,particle_set,weights,error)

    TYPE(qs_loc_env_new_type), POINTER         :: qs_loc_env
    TYPE(cell_type), OPTIONAL, POINTER         :: cell
    TYPE(cp_fm_p_type), DIMENSION(:),&
      POINTER, OPTIONAL                        :: moloc_coeff
    TYPE(cp_fm_p_type), DIMENSION(:,:),&
      POINTER, OPTIONAL                        :: op_fm_set
   TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                  :: local_molecules
    TYPE(localized_wfn_control_type), &
      POINTER, OPTIONAL                        :: localized_wfn_control
    TYPE (real_matrix_p_type), DIMENSION(:,:),&
      OPTIONAL, POINTER                        :: op_sm_set
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                  :: para_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                        :: particle_set
    REAL(dp), DIMENSION(6), OPTIONAL           :: weights
    TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                                :: error


    CHARACTER(len=*), PARAMETER :: routine_name = 'set_qs_loc_env', &
      routineP = module_name//':'//routine_name

    INTEGER                                    :: istat
    LOGICAL                                    :: failure

     failure =.FALSE.
     CPPrecondition(ASSOCIATED(qs_loc_env),cp_failure_level,routineP,error,failure)

     IF(.NOT.failure) THEN
     IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(qs_loc_env%cell,error=error)
       qs_loc_env%cell => cell
     END IF

     IF (PRESENT(local_molecules)) THEN
        CALL distribution_1d_retain(local_molecules,error=error)
        IF(ASSOCIATED(qs_loc_env%local_molecules)) &
             CALL distribution_1d_release(qs_loc_env%local_molecules,error=error)
        qs_loc_env%local_molecules => local_molecules
     END IF

     IF(PRESENT(localized_wfn_control)) THEN
       CALL localized_wfn_control_retain(localized_wfn_control,error=error)
       CALL localized_wfn_control_release(qs_loc_env%localized_wfn_control,error=error)
       qs_loc_env % localized_wfn_control => localized_wfn_control
     END IF
     IF(PRESENT(para_env)) THEN
      CALL cp_para_env_retain(para_env)
      CALL cp_para_env_release(qs_loc_env%para_env)
      qs_loc_env%para_env => para_env
     END IF
     IF (PRESENT(particle_set)) qs_loc_env%particle_set => particle_set
     IF(PRESENT(moloc_coeff)) THEN
       qs_loc_env%moloc_coeff => moloc_coeff
     END IF
     IF(PRESENT(op_sm_set)) THEN
       qs_loc_env%op_sm_set => op_sm_set
     END IF
     IF(PRESENT(op_fm_set)) THEN
       qs_loc_env%op_fm_set => op_fm_set
     END IF
     IF(PRESENT(weights)) THEN
       qs_loc_env%weights = weights
     END IF
     END IF
   
  END SUBROUTINE set_qs_loc_env
!******************************************************************************

END MODULE qs_loc_types

