!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types needed for a linear scaling quickstep SCF run based on the density
!>        matrix
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_release
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: ls_cluster_atomic,&
                                             ls_cluster_molecular
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_types'

  PUBLIC :: ls_scf_env_type, ls_mstruct_type, ls_cluster_atomic, ls_cluster_molecular

  TYPE ls_mstruct_type
     INTEGER :: cluster_type
     LOGICAL :: single_precision
     INTEGER, DIMENSION(:), ALLOCATABLE :: atom_to_molecule
  END TYPE

  TYPE ls_mat_history_type
     INTEGER :: istore,nstore
     TYPE(cp_dbcsr_type), DIMENSION(:,:), ALLOCATABLE :: matrix
  END TYPE

  TYPE ls_scf_env_type
     INTEGER               :: nspins, natoms
     INTEGER               :: nelectron_total
     INTEGER, DIMENSION(2) :: nelectron_spin
     REAL(KIND=dp), DIMENSION(2) ::  mu_spin
     REAL(KIND=dp), DIMENSION(2) ::  homo_spin
     REAL(KIND=dp), DIMENSION(2) ::  lumo_spin

     TYPE(ls_mat_history_type) :: scf_history
     INTEGER :: extrapolation_order

     LOGICAL :: has_unit_metric

     INTEGER :: s_preconditioner_type
     INTEGER :: s_inversion_type
     INTEGER :: purification_method
     INTEGER :: sign_sqrt_order

     LOGICAL               :: needs_s_inv, has_s_preconditioner, fixed_mu, dynamic_threshold
     REAL(KIND=dp)         :: eps_filter, eps_scf

     REAL(KIND=dp)         :: mixing_fraction
     INTEGER               :: max_scf

     TYPE(cp_dbcsr_type)   :: matrix_s_inv
     TYPE(cp_dbcsr_type)   :: matrix_s
     TYPE(cp_dbcsr_type)   :: matrix_bs_sqrt, matrix_bs_sqrt_inv
     TYPE(cp_dbcsr_type)   :: matrix_s_sqrt, matrix_s_sqrt_inv
     TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks
     TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p

     LOGICAL  :: report_all_sparsities, perform_mu_scan, use_s_sqrt

     TYPE(ls_mstruct_type) :: ls_mstruct

     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE ls_scf_env_type

  PUBLIC :: ls_scf_release

CONTAINS

! *****************************************************************************
!> \brief release the LS type.
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE ls_scf_release(ls_scf_env,error)
    TYPE(ls_scf_env_type), POINTER           :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ls_scf_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, istat, istore
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CALL cp_para_env_release(ls_scf_env%para_env,error)

    DEALLOCATE(ls_scf_env%ls_mstruct%atom_to_molecule,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! set up the buffer for the history of matrices
    DO istore=1,MIN(ls_scf_env%scf_history%istore,ls_scf_env%scf_history%nstore)
       DO ispin=1,SIZE(ls_scf_env%scf_history%matrix,1)
          CALL cp_dbcsr_release(ls_scf_env%scf_history%matrix(ispin,istore),error=error)
       ENDDO
    ENDDO
    DEALLOCATE(ls_scf_env%scf_history%matrix,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)


    DEALLOCATE(ls_scf_env,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE ls_scf_release

END MODULE dm_ls_scf_types


