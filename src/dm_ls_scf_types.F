!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for a linear scaling quickstep SCF run based on the density
!>        matrix
!> \par History
!>       2010.10 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE dm_ls_scf_types
   USE dbcsr_api,                       ONLY: dbcsr_release,&
                                              dbcsr_type
   USE input_constants,                 ONLY: ls_cluster_atomic,&
                                              ls_cluster_molecular
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE pao_types,                       ONLY: pao_env_type,&
                                              pao_finalize
   USE pexsi_types,                     ONLY: lib_pexsi_env,&
                                              lib_pexsi_finalize
   USE qs_density_mixing_types,         ONLY: mixing_storage_release,&
                                              mixing_storage_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dm_ls_scf_types'

   PUBLIC :: ls_scf_env_type, ls_mstruct_type, ls_cluster_atomic, ls_cluster_molecular, &
             ls_scf_curvy_type

   TYPE ls_mstruct_type
      INTEGER :: cluster_type = -1
      LOGICAL :: single_precision = .FALSE.
      LOGICAL :: do_pao = .FALSE.
      INTEGER, DIMENSION(:), ALLOCATABLE :: atom_to_molecule
      TYPE(dbcsr_type)                :: matrix_A
      TYPE(dbcsr_type)                :: matrix_B
   END TYPE

   TYPE ls_mat_history_type
      INTEGER :: istore = 0, nstore = 0
      TYPE(dbcsr_type), DIMENSION(:, :), ALLOCATABLE :: matrix
   END TYPE

   TYPE ls_scf_curvy_type
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_dp
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p
      TYPE(dbcsr_type), DIMENSION(:, :), ALLOCATABLE :: matrix_psave
      TYPE(dbcsr_type), DIMENSION(:, :), ALLOCATABLE :: matrix_BCH
      REAL(KIND=dp), DIMENSION(2)                      :: step_size = 0.0_dp
      REAL(KIND=dp), DIMENSION(2)                      :: shift = 0.0_dp
      REAL(KIND=dp), DIMENSION(2)                      :: cg_denom = 0.0_dp
      REAL(KIND=dp), DIMENSION(2)                      :: cg_numer = 0.0_dp
      REAL(KIND=dp), DIMENSION(6)                      :: energies = 0.0_dp
      INTEGER                                          :: line_search_step = 0
      INTEGER, DIMENSION(2)                            :: BCH_saved = 0
      LOGICAL                                          :: double_step_size = .FALSE.
      LOGICAL, DIMENSION(2)                            :: fix_shift = .FALSE.

      INTEGER                                          :: line_search_type = 0
      INTEGER                                          :: n_bch_hist = 0
      REAL(KIND=dp)                                    :: scale_filter = 0.0_dp
      REAL(KIND=dp)                                    :: filter_factor = 0.0_dp
      REAL(KIND=dp)                                    :: min_shift = 0.0_dp
      REAL(KIND=dp)                                    :: min_filter = 0.0_dp
   END TYPE

   TYPE chebyshev_type
      LOGICAL :: compute_chebyshev = .FALSE.
      INTEGER :: n_chebyshev = 0
      INTEGER :: n_gridpoint_dos = 0
      REAL(KIND=dp), DIMENSION(:), POINTER :: min_energy => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER :: max_energy => NULL()
      TYPE(section_vals_type), POINTER :: print_key_dos => NULL()
      TYPE(section_vals_type), POINTER :: print_key_cube => NULL()
   END TYPE

   TYPE ls_scf_env_type
      INTEGER               :: nspins = 0, natoms = 0
      INTEGER               :: nelectron_total = 0
      INTEGER, DIMENSION(2) :: nelectron_spin = 0
      REAL(KIND=dp), DIMENSION(2) ::  mu_spin = 0.0_dp
      REAL(KIND=dp), DIMENSION(2) ::  homo_spin = 0.0_dp
      REAL(KIND=dp), DIMENSION(2) ::  lumo_spin = 0.0_dp

      ! This prevents a bug in GCC up to version 9.x (9.5 seems to work though)
      TYPE(ls_mat_history_type) :: scf_history = ls_mat_history_type(matrix=null())
      INTEGER :: extrapolation_order = -1

      LOGICAL :: has_unit_metric = .FALSE.

      LOGICAL :: curvy_steps = .FALSE.
      INTEGER :: s_preconditioner_type = 0
      INTEGER :: s_inversion_type = 0
      INTEGER :: purification_method = 0
      INTEGER :: sign_method = 0
      INTEGER :: sign_order = 0
      LOGICAL :: sign_symmetric = .FALSE.
      INTEGER :: submatrix_sign_method = -1
      INTEGER :: s_sqrt_method = 0
      INTEGER :: s_sqrt_order = 0

      LOGICAL               :: needs_s_inv = .FALSE., has_s_preconditioner = .FALSE., fixed_mu = .FALSE., &
                               dynamic_threshold = .FALSE., check_s_inv = .FALSE.
      LOGICAL               :: restart_read = .FALSE., restart_write = .FALSE., non_monotonic = .FALSE.
      REAL(KIND=dp)         :: eps_filter = 0.0_dp, eps_scf = 0.0_dp

      REAL(KIND=dp)         :: eps_lanczos = 0.0_dp
      INTEGER               :: max_iter_lanczos = 0

      REAL(KIND=dp)         :: mixing_fraction = 0.0_dp
      INTEGER               :: max_scf = 0
      LOGICAL               :: ls_diis = .FALSE.
      INTEGER               :: iter_ini_diis = 0
      INTEGER               :: nmixing = 0, max_diis = 0
      REAL(KIND=dp)         :: eps_diis = 0.0_dp
      REAL(KIND=dp)         :: energy_init = 0.0_dp

      TYPE(dbcsr_type)   :: matrix_s_inv
      TYPE(dbcsr_type)   :: matrix_s
      TYPE(dbcsr_type)   :: matrix_bs_sqrt, matrix_bs_sqrt_inv
      TYPE(dbcsr_type)   :: matrix_s_sqrt, matrix_s_sqrt_inv
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p

      LOGICAL  :: report_all_sparsities = .FALSE., perform_mu_scan = .FALSE., use_s_sqrt = .FALSE.

      TYPE(ls_mstruct_type) :: ls_mstruct
      TYPE(ls_scf_curvy_type) :: curvy_data = ls_scf_curvy_type(matrix_dp=null(), matrix_p=null(), &
                                                                matrix_psave=null(), matrix_bch=null())

      TYPE(chebyshev_type) :: chebyshev = chebyshev_type()

      LOGICAL :: do_rho_mixing = .FALSE.
      INTEGER :: density_mixing_method = 0
      TYPE(mixing_storage_type), POINTER :: mixing_store => NULL()

      LOGICAL :: do_transport = .FALSE.
      LOGICAL :: do_pexsi = .FALSE.

      LOGICAL :: calculate_forces = .FALSE.

      TYPE(lib_pexsi_env) :: pexsi

      TYPE(mp_para_env_type), POINTER :: para_env => NULL()
      LOGICAL                 :: do_pao = .FALSE.
      TYPE(pao_env_type)      :: pao_env
   END TYPE ls_scf_env_type

   PUBLIC :: ls_scf_release

CONTAINS

! **************************************************************************************************
!> \brief release the LS type.
!> \param ls_scf_env ...
!> \par History
!>       2012.11 created [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE ls_scf_release(ls_scf_env)
      TYPE(ls_scf_env_type), POINTER                     :: ls_scf_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ls_scf_release'

      INTEGER                                            :: handle, ispin, istore

      CALL timeset(routineN, handle)

      CALL mp_para_env_release(ls_scf_env%para_env)

      DEALLOCATE (ls_scf_env%ls_mstruct%atom_to_molecule)

      ! set up the buffer for the history of matrices
      DO istore = 1, MIN(ls_scf_env%scf_history%istore, ls_scf_env%scf_history%nstore)
         DO ispin = 1, SIZE(ls_scf_env%scf_history%matrix, 1)
            CALL dbcsr_release(ls_scf_env%scf_history%matrix(ispin, istore))
         END DO
      END DO
      DEALLOCATE (ls_scf_env%scf_history%matrix)

      IF (ALLOCATED(ls_scf_env%matrix_p)) THEN
         DO ispin = 1, SIZE(ls_scf_env%matrix_p)
            CALL dbcsr_release(ls_scf_env%matrix_p(ispin))
         END DO
         DEALLOCATE (ls_scf_env%matrix_p)
      END IF

      IF (ASSOCIATED(ls_scf_env%chebyshev%print_key_dos)) &
         CALL section_vals_release(ls_scf_env%chebyshev%print_key_dos)
      IF (ASSOCIATED(ls_scf_env%chebyshev%print_key_cube)) &
         CALL section_vals_release(ls_scf_env%chebyshev%print_key_cube)
      IF (ASSOCIATED(ls_scf_env%chebyshev%min_energy)) THEN
         DEALLOCATE (ls_scf_env%chebyshev%min_energy)
      END IF
      IF (ASSOCIATED(ls_scf_env%chebyshev%max_energy)) THEN
         DEALLOCATE (ls_scf_env%chebyshev%max_energy)
      END IF

      IF (ASSOCIATED(ls_scf_env%mixing_store)) THEN
         CALL mixing_storage_release(ls_scf_env%mixing_store)
         DEALLOCATE (ls_scf_env%mixing_store)
      END IF

      IF (ls_scf_env%do_pexsi) THEN
         CALL lib_pexsi_finalize(ls_scf_env%pexsi)
      END IF

      IF (ls_scf_env%do_pao) &
         CALL pao_finalize(ls_scf_env%pao_env)

      DEALLOCATE (ls_scf_env)

      CALL timestop(handle)

   END SUBROUTINE ls_scf_release

END MODULE dm_ls_scf_types
