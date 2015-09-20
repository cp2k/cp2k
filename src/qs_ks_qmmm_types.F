!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      05.2004 [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qs_ks_qmmm_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_release,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_qmmm_types'

  PUBLIC :: qs_ks_qmmm_env_type
  PUBLIC :: qs_ks_qmmm_release, qs_ks_qmmm_retain

! *****************************************************************************
!> \brief calculation environement to calculate the ks_qmmm matrix,
!>      holds the QM/MM potential and all the needed variables to
!>      compute the QM/MM electrostatic 1-electron ks matrix
!>      assumes that the core hamiltonian and energy are up to date.
!>      v_metal_rspace is the potential at the metal sites within the image
!>      charge approach
!> \par History
!>      05.2004 created [tlaino]
!>      01.2012 added v_metal_rspace [dgolze]
!> \author Teodoro Laino
! *****************************************************************************
  TYPE qs_ks_qmmm_env_type
     INTEGER :: n_evals, &
                id_nr, ref_count
     REAL(KIND=dp)                               :: pc_ener
     TYPE(pw_env_type), POINTER                  :: pw_env
     TYPE(pw_p_type)                             :: v_qmmm_rspace
     TYPE(pw_p_type)                             :: v_metal_rspace
     TYPE(cube_info_type),DIMENSION(:), POINTER  :: cube_info
     TYPE(cp_dbcsr_p_type), DIMENSION(:), &
          POINTER                                :: matrix_h
  END TYPE qs_ks_qmmm_env_type

! *****************************************************************************
!> \brief type to build arrays of pointers
!> \param ks_env the ks_env pointer
!> \par History
!>      05.2004 [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  TYPE qs_ks_qmmm_env_p_type
     TYPE(qs_ks_qmmm_env_type), POINTER :: ks_env
  END TYPE qs_ks_qmmm_env_p_type
CONTAINS

! *****************************************************************************
!> \brief releases the ks_qmmm_env (see doc/ReferenceCounting.html)
!> \param ks_qmmm_env the ks_qmmm_env to be released
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qs_ks_qmmm_release(ks_qmmm_env)
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(pw_pool_type), POINTER              :: pool

    failure=.FALSE.

    IF (ASSOCIATED(ks_qmmm_env)) THEN
       CPASSERT(ks_qmmm_env%ref_count>0)
       ks_qmmm_env%ref_count=ks_qmmm_env%ref_count-1

       IF (ks_qmmm_env%ref_count<1) THEN
          CALL pw_env_get(ks_qmmm_env%pw_env,auxbas_pw_pool=pool)
          CALL pw_pool_give_back_pw(pool,ks_qmmm_env%v_qmmm_rspace%pw)
          CALL pw_env_release(ks_qmmm_env%pw_env)
          IF (ASSOCIATED(ks_qmmm_env%cube_info))THEN
             DO i=1,SIZE(ks_qmmm_env%cube_info)
                CALL destroy_cube_info(ks_qmmm_env%cube_info(i))
             END DO
             DEALLOCATE(ks_qmmm_env%cube_info, stat=stat)
             CPASSERT(stat==0)
          END IF
          IF (ASSOCIATED(ks_qmmm_env%matrix_h)) THEN
             CALL cp_dbcsr_deallocate_matrix_set(ks_qmmm_env%matrix_h)
          END IF
          DEALLOCATE(ks_qmmm_env, stat=stat)
          CPASSERT(stat==0)
       END IF
    END IF
    NULLIFY(ks_qmmm_env)
  END SUBROUTINE qs_ks_qmmm_release

! *****************************************************************************
!> \brief retains the given ks_environment
!> \param ks_qmmm_env the KohnSham QM/MM environment to retain
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
SUBROUTINE qs_ks_qmmm_retain(ks_qmmm_env)
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(ks_qmmm_env))
  CPASSERT(ks_qmmm_env%ref_count>0)
  ks_qmmm_env%ref_count=ks_qmmm_env%ref_count+1
END SUBROUTINE qs_ks_qmmm_retain

END MODULE qs_ks_qmmm_types
