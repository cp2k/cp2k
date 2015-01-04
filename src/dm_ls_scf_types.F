!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types needed for a linear scaling quickstep SCF run based on the density
!>        matrix
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE dm_ls_scf_types
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_release,&
                                             cp_dbcsr_type
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: ls_cluster_atomic,&
                                             ls_cluster_molecular
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_types'

  PUBLIC :: ls_scf_env_type, ls_mstruct_type, ls_cluster_atomic, ls_cluster_molecular,&
            ls_scf_curvy_type

  TYPE ls_mstruct_type
     INTEGER :: cluster_type
     LOGICAL :: single_precision
     INTEGER, DIMENSION(:), ALLOCATABLE :: atom_to_molecule
  END TYPE

  TYPE ls_mat_history_type
     INTEGER :: istore,nstore
     TYPE(cp_dbcsr_type), DIMENSION(:,:), ALLOCATABLE :: matrix
  END TYPE

  TYPE ls_scf_curvy_type
      TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_dp
      TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p
      TYPE(cp_dbcsr_type), DIMENSION(:,:), ALLOCATABLE :: matrix_psave
      TYPE(cp_dbcsr_type), DIMENSION(:,:), ALLOCATABLE :: matrix_BCH
      REAL(KIND=dp)                                    :: step_size(2)
      REAL(KIND=dp)                                    :: shift(2)
      REAL(KIND=dp)                                    :: cg_denom(2)
      REAL(KIND=dp)                                    :: cg_numer(2)
      REAL(KIND=dp), DIMENSION(6)                      :: energies
      INTEGER                                          :: line_search_step
      INTEGER                                          :: BCH_saved(2)
      LOGICAL                                          :: double_step_size
      LOGICAL                                          :: fix_shift(2)

      INTEGER                                          :: line_search_type
      INTEGER                                          :: n_bch_hist
      REAL(KIND=dp)                                    :: scale_filter
      REAL(KIND=dp)                                    :: filter_factor
      REAL(KIND=dp)                                    :: min_shift
      REAL(KIND=dp)                                    :: min_filter
  END TYPE


  TYPE chebyshev_type
     LOGICAL :: compute_chebyshev
     INTEGER :: n_chebyshev
     INTEGER :: n_gridpoint_dos
     REAL(KIND=dp), DIMENSION(:), POINTER :: min_energy=>NULL()
     REAL(KIND=dp), DIMENSION(:), POINTER :: max_energy=>NULL()
     TYPE(section_vals_type), POINTER :: print_key_dos=>NULL()
     TYPE(section_vals_type), POINTER :: print_key_cube=>NULL()
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

     LOGICAL :: curvy_steps
     INTEGER :: s_preconditioner_type
     INTEGER :: s_inversion_type
     INTEGER :: purification_method
     INTEGER :: sign_sqrt_order

     LOGICAL               :: needs_s_inv, has_s_preconditioner, fixed_mu, dynamic_threshold
     LOGICAL               :: restart_read, restart_write, non_monotonic
     REAL(KIND=dp)         :: eps_filter, eps_scf

     REAL(KIND=dp)         :: eps_lanczos
     INTEGER               :: max_iter_lanczos

     REAL(KIND=dp)         :: mixing_fraction
     INTEGER               :: max_scf
     LOGICAL               :: ls_diis
     INTEGER               :: iter_ini_diis
     INTEGER               :: nmixing, max_diis
     REAL(KIND=dp)         :: eps_diis
     REAL(KIND=dp)         :: energy_init

     TYPE(cp_dbcsr_type)   :: matrix_s_inv
     TYPE(cp_dbcsr_type)   :: matrix_s
     TYPE(cp_dbcsr_type)   :: matrix_bs_sqrt, matrix_bs_sqrt_inv
     TYPE(cp_dbcsr_type)   :: matrix_s_sqrt, matrix_s_sqrt_inv
     TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks
     TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p

     LOGICAL  :: report_all_sparsities, perform_mu_scan, use_s_sqrt

     TYPE(ls_mstruct_type) :: ls_mstruct
     TYPE(ls_scf_curvy_type) :: curvy_data

     TYPE(chebyshev_type) :: chebyshev

     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE ls_scf_env_type

  PUBLIC :: ls_scf_release

CONTAINS

! *****************************************************************************
!> \brief release the LS type.
!> \param ls_scf_env ...
!> \param error ...
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

    IF (ASSOCIATED(ls_scf_env%chebyshev%print_key_dos)) &
        CALL section_vals_release(ls_scf_env%chebyshev%print_key_dos,error=error)
    IF (ASSOCIATED(ls_scf_env%chebyshev%print_key_cube)) &
        CALL section_vals_release(ls_scf_env%chebyshev%print_key_cube,error=error)
    IF (ASSOCIATED(ls_scf_env%chebyshev%min_energy)) THEN
       DEALLOCATE(ls_scf_env%chebyshev%min_energy,stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(ls_scf_env%chebyshev%max_energy)) THEN
       DEALLOCATE(ls_scf_env%chebyshev%max_energy,stat=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    DEALLOCATE(ls_scf_env,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE ls_scf_release

END MODULE dm_ls_scf_types


