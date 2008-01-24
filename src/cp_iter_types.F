!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/cp_iter_types [1.0] *
!!
!!   NAME
!!     cp_iter_types
!!
!!   FUNCTION
!!     Collection of routines to handle the iteration info
!!
!!   NOTES
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE cp_iter_types
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_abort

  IMPLICIT NONE
  PRIVATE

  ! iteration_info
  PUBLIC :: cp_iteration_info_type, cp_iteration_info_create, cp_iteration_info_retain, &
            cp_iteration_info_release

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_iter_types'
  LOGICAL, PRIVATE, PARAMETER          :: debug_this_module=.FALSE.
  INTEGER, SAVE, PRIVATE               :: last_it_info_id=0

  ! When adding a new iteration level PLEASE update the following list with the proper name!
  CHARACTER(LEN=default_path_length), PARAMETER, PUBLIC, DIMENSION(12) :: each_possible_labels= (/&
       "__ROOT__          ",&
       "JUST_ENERGY       ",&
       "QS_SCF            ",&
       "XAS_SCF           ",&
       "MD                ",&
       "METADYNAMICS      ",&
       "GEO_OPT           ",&
       "ROT_OPT           ",&
       "BAND              ",&
       "EP_LIN_SOLVER     ",&
       "SPLINE_FIND_COEFFS",&
       "REPLICA_EVAL      "/)

  CHARACTER(LEN=default_path_length), PARAMETER, PUBLIC, DIMENSION(12) ::  each_desc_labels= (/&
       "Iteration level for __ROOT__ (fictitious iteration level)                      ",&
       "Iteration level for an ENERGY/ENERGY_FORCE calculation.                        ",&
       "Iteration level for the SCF Steps.                                             ",&
       "Iteration level for the X-Ray Absorption Spectroscopy (XAS) SCF Steps.         ",&
       "Iteration level for the MD steps.                                              ",&
       "Iteration level for the METADYNAMICS steps (number of hills added).            ",&
       "Iteration level for the geometry optimization steps.                           ",&
       "Iteration level for the rotational optimization steps in the Dimer Calculation.",&
       "Iteration level for the Band Calculation Steps                                 ",&
       "Iteration level for the Energy Perturbation (EP) linear solver                 ",&
       "Iteration level for the solution of the coefficients of the splines            ",&
       "Iteration level for the evaluation of the Replica Environment                  "/)

!!****s* cp_iter_types/cp_iteration_info_type *
!!
!!   NAME
!!     cp_iteration_info
!!
!!   FUNCTION
!!     contains the information about the current state of the program
!!     to be able to decide if output is necessary
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_iteration_info_type
     INTEGER                              :: ref_count, id_nr
     INTEGER                              :: print_level, n_rlevel
     INTEGER, DIMENSION(:), POINTER       :: iteration
     LOGICAL, DIMENSION(:), POINTER       :: last_iter
     CHARACTER(len=default_string_length) :: project_name
     CHARACTER(LEN=default_string_length),&
          DIMENSION(:), POINTER           :: level_name
  END TYPE cp_iteration_info_type

CONTAINS

!!****f* cp_iter_types/cp_iteration_info_create *
!!
!!   NAME
!!     cp_iteration_info_create
!!
!!   FUNCTION
!!     creates an output info object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - iteration_info: the object to create
!!     - project_name: name of the project, used to create the filenames
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE cp_iteration_info_create(iteration_info,project_name)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    CHARACTER(len=*), INTENT(in)             :: project_name

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iteration_info_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    ALLOCATE(iteration_info,stat=stat)
    IF (stat/=0) THEN
       CALL mp_abort(routineP//" could not allocate iteration_info")
    END IF
    last_it_info_id = last_it_info_id+1
    iteration_info%id_nr = last_it_info_id
    iteration_info%ref_count  = 1
    iteration_info%print_level= 2
    iteration_info%n_rlevel   = 1
    iteration_info%project_name = project_name
    NULLIFY(iteration_info%iteration)
    NULLIFY(iteration_info%level_name)
    NULLIFY(iteration_info%last_iter)
    ALLOCATE(iteration_info%iteration(iteration_info%n_rlevel),stat=stat)
    IF (stat/=0) THEN
       CALL mp_abort(routineP//" iteration_info%iteration allocation")
    END IF
    ALLOCATE(iteration_info%level_name(iteration_info%n_rlevel),stat=stat)
    IF (stat/=0) THEN
       CALL mp_abort(routineP//" iteration_info%level_name allocation")
    END IF
    ALLOCATE(iteration_info%last_iter(iteration_info%n_rlevel),stat=stat)
    IF (stat/=0) THEN
       CALL mp_abort(routineP//" iteration_info%last_iter allocation")
    END IF
    iteration_info%iteration(iteration_info%n_rlevel)  = 1
    iteration_info%level_name(iteration_info%n_rlevel) = "__ROOT__"
    iteration_info%last_iter(iteration_info%n_rlevel)  = .FALSE.

  END SUBROUTINE cp_iteration_info_create
!***************************************************************************

!!****f* cp_iter_types/cp_iteration_info_retain *
!!
!!   NAME
!!     cp_iteration_info_retain
!!
!!   FUNCTION
!!     retains the iteration_info (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - iteration_info: the iteration_info to retain
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE cp_iteration_info_retain(iteration_info)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iteration_info_retain', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(iteration_info)) THEN
       CALL mp_abort(routineP//" iteration_info not associated")
    END IF
    IF (iteration_info%ref_count<=0) THEN
       CALL mp_abort(routineP//" iteration_info%ref_counf<=0")
    END IF
    iteration_info%ref_count=iteration_info%ref_count+1
  END SUBROUTINE cp_iteration_info_retain
!***************************************************************************

!!****f* cp_iter_types/cp_iteration_info_release *
!!
!!   NAME
!!     cp_iteration_info_release
!!
!!   FUNCTION
!!     releases the iteration_info (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - iteration_info: the iteration_info to release
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE cp_iteration_info_release(iteration_info)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iteration_info_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    IF (ASSOCIATED(iteration_info)) THEN
       IF (iteration_info%ref_count<=0) THEN
          CALL mp_abort(routineP//" iteration_info%ref_counf<=0")
       END IF
       iteration_info%ref_count=iteration_info%ref_count-1
       IF (iteration_info%ref_count==0) THEN
          IF (ASSOCIATED(iteration_info%iteration)) THEN
             DEALLOCATE(iteration_info%iteration,stat=stat)
             IF (stat/=0) THEN
                CALL mp_abort(routineP//" iteration_info%iteration deallocation stat ")
             END IF
          END IF
          IF (ASSOCIATED(iteration_info%last_iter)) THEN
             DEALLOCATE(iteration_info%last_iter,stat=stat)
             IF (stat/=0) THEN
                CALL mp_abort(routineP//" iteration_info%last_iter deallocation stat ")
             END IF
          END IF
          IF (ASSOCIATED(iteration_info%level_name)) THEN
             DEALLOCATE(iteration_info%level_name,stat=stat)
             IF (stat/=0) THEN
                CALL mp_abort(routineP//" iteration_info%level_name deallocation stat ")
             END IF
          END IF
          DEALLOCATE(iteration_info,stat=stat)
          IF (stat/=0) THEN
             CALL mp_abort(routineP//" iteration_info deallocation stat ")
          END IF
       END IF
    END IF
  END SUBROUTINE cp_iteration_info_release
!***************************************************************************

END MODULE cp_iter_types

!******************************************************************************

