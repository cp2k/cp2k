!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types of  preconditioners
!> \par History
!>      Separate types from construction and application
!> \author Joost VandeVondele (09.2002)
! **************************************************************************************************
MODULE preconditioner_types
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_release,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_release_p,&
                                              dbcsr_type
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE input_constants,                 ONLY: cholesky_reduce,&
                                              ot_lattice_fft_off,&
                                              ot_low_rank_base_overlap,&
                                              ot_precond_solver_default
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner_types'

   PUBLIC  :: preconditioner_type, preconditioner_p_type
   PUBLIC  :: lattice_rhs_distribution_type
   PUBLIC  :: init_preconditioner
   PUBLIC  :: preconditioner_in_use
   PUBLIC  :: destroy_preconditioner
   PUBLIC  :: release_lattice_rhs_distribution

! **************************************************************************************************
   TYPE lattice_rhs_distribution_type
      INTEGER :: nrow_global = 0, ncol_global = 0
      INTEGER :: nrow_local = 0, ncol_local = 0, num_pe = 0
      INTEGER, ALLOCATABLE, DIMENSION(:) :: global_columns
      INTEGER, ALLOCATABLE, DIMENSION(:) :: local_rows, local_columns
      INTEGER, ALLOCATABLE, DIMENSION(:) :: send_counts, send_displacements
      INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_counts, recv_displacements
      INTEGER, ALLOCATABLE, DIMENSION(:) :: send_local_rows, send_local_columns
      INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_global_rows, recv_global_columns
   END TYPE lattice_rhs_distribution_type

! **************************************************************************************************
   TYPE preconditioner_type
!    PRIVATE
      TYPE(dbcsr_type), POINTER :: sparse_matrix => NULL()
      TYPE(cp_fm_type), POINTER :: fm => NULL()
      TYPE(cp_cfm_type), POINTER :: complex_fm => NULL()
      LOGICAL :: owns_complex_fm = .TRUE.
      TYPE(cp_fm_type), POINTER :: base_fm => NULL()
      TYPE(dbcsr_type), POINTER           :: dbcsr_matrix => NULL()
      TYPE(dbcsr_type), POINTER           :: base_dbcsr_matrix => NULL()
      TYPE(dbcsr_type), POINTER           :: occ_rotation_dbcsr => NULL()
      TYPE(dbcsr_type), POINTER           :: max_ev_vector => NULL()
      TYPE(dbcsr_type), POINTER           :: min_ev_vector => NULL()
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:) :: inverse_history => NULL()
      TYPE(mp_para_env_type), POINTER   :: para_env => NULL()
      TYPE(cp_blacs_env_type), POINTER   :: ctxt => NULL()
      TYPE(cp_fm_type), POINTER          :: occ_rotation => NULL()
      INTEGER :: in_use = -1, solver = -1, ihistory = -1, cholesky_use = -1
      REAL(KIND=dp), DIMENSION(:), POINTER :: occ_evals => NULL(), full_evals => NULL()
      REAL(KIND=dp) :: energy_gap = -1.0_dp
      REAL(KIND=dp) :: condition_num = -1.0_dp
      INTEGER :: polynomial_degree = 8
      REAL(KIND=dp) :: polynomial_min = -1.0_dp, polynomial_max = -1.0_dp
      INTEGER :: spectral_rank = 0
      INTEGER :: spectral_max_rank = 48
      REAL(KIND=dp) :: spectral_reference = 0.0_dp, spectral_window = 1.0_dp, &
                       spectral_base_scale = 1.0_dp
      INTEGER :: low_rank_base = ot_low_rank_base_overlap
      INTEGER :: lattice_fft_mode = ot_lattice_fft_off
      LOGICAL :: lattice_fft_active = .FALSE.
      INTEGER, DIMENSION(3) :: lattice_dims = [1, 1, 1]
      INTEGER :: lattice_block_size = 0
      INTEGER :: lattice_local_cells = 0
      REAL(KIND=dp) :: lattice_projection_error = 0.0_dp, &
                       lattice_corrected_projection_error = 0.0_dp, &
                       lattice_cost_ratio = 0.0_dp, &
                       lattice_storage_ratio = 0.0_dp
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: lattice_inverse_k
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: lattice_state_inverse_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: lattice_occ_vectors, lattice_occ_dual
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lattice_local_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: lattice_local_operator, &
                                                     lattice_local_coarse_inverse
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: lattice_state_local_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: lattice_state_local_operator, &
                                                        lattice_state_local_coarse_inverse
      TYPE(lattice_rhs_distribution_type) :: lattice_rhs_distribution
   END TYPE preconditioner_type

! **************************************************************************************************
   TYPE preconditioner_p_type
      TYPE(preconditioner_type), POINTER :: preconditioner => NULL()
   END TYPE preconditioner_p_type

CONTAINS

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param preconditioner ...
!> \return ...
! **************************************************************************************************
   FUNCTION preconditioner_in_use(preconditioner)
      TYPE(preconditioner_type)                          :: preconditioner
      LOGICAL                                            :: preconditioner_in_use

      preconditioner_in_use = .NOT. (preconditioner%in_use == 0)
   END FUNCTION preconditioner_in_use

! **************************************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param para_env ...
!> \param blacs_env ...
! **************************************************************************************************
   SUBROUTINE init_preconditioner(preconditioner_env, para_env, blacs_env)

      TYPE(preconditioner_type)                          :: preconditioner_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      NULLIFY (preconditioner_env%sparse_matrix)
      NULLIFY (preconditioner_env%fm)
      NULLIFY (preconditioner_env%complex_fm)
      preconditioner_env%owns_complex_fm = .TRUE.
      NULLIFY (preconditioner_env%base_fm)
      NULLIFY (preconditioner_env%dbcsr_matrix)
      NULLIFY (preconditioner_env%base_dbcsr_matrix)
      NULLIFY (preconditioner_env%occ_rotation_dbcsr)
      NULLIFY (preconditioner_env%occ_rotation)
      NULLIFY (preconditioner_env%occ_evals)
      NULLIFY (preconditioner_env%full_evals)
      NULLIFY (preconditioner_env%inverse_history)
      NULLIFY (preconditioner_env%max_ev_vector)
      NULLIFY (preconditioner_env%min_ev_vector)
      preconditioner_env%solver = ot_precond_solver_default
      preconditioner_env%para_env => para_env
      preconditioner_env%ctxt => blacs_env
      !inverse is used for filtering in update set it to something huge to
      ! avoid filtering if the information is not available
      preconditioner_env%condition_num = -1.0_dp
      preconditioner_env%polynomial_degree = 8
      preconditioner_env%polynomial_min = -1.0_dp
      preconditioner_env%polynomial_max = -1.0_dp
      preconditioner_env%spectral_rank = 0
      preconditioner_env%spectral_max_rank = 48
      preconditioner_env%spectral_reference = 0.0_dp
      preconditioner_env%spectral_window = 1.0_dp
      preconditioner_env%spectral_base_scale = 1.0_dp
      preconditioner_env%low_rank_base = ot_low_rank_base_overlap
      preconditioner_env%lattice_fft_mode = ot_lattice_fft_off
      preconditioner_env%lattice_fft_active = .FALSE.
      preconditioner_env%lattice_dims = [1, 1, 1]
      preconditioner_env%lattice_block_size = 0
      preconditioner_env%lattice_local_cells = 0
      preconditioner_env%lattice_projection_error = 0.0_dp
      preconditioner_env%lattice_corrected_projection_error = 0.0_dp
      preconditioner_env%lattice_cost_ratio = 0.0_dp
      preconditioner_env%lattice_storage_ratio = 0.0_dp
      preconditioner_env%ihistory = 0

      CALL preconditioner_env%para_env%retain()
      CALL preconditioner_env%ctxt%retain()

   END SUBROUTINE init_preconditioner

! **************************************************************************************************
!> \brief ...
!> \param preconditioner_env ...
! **************************************************************************************************
   SUBROUTINE destroy_preconditioner(preconditioner_env)

      TYPE(preconditioner_type)                          :: preconditioner_env

      CHARACTER(len=*), PARAMETER :: routineN = 'destroy_preconditioner'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(preconditioner_env%sparse_matrix)) THEN
         CALL dbcsr_deallocate_matrix(preconditioner_env%sparse_matrix)
         NULLIFY (preconditioner_env%sparse_matrix)
      END IF

      IF (ASSOCIATED(preconditioner_env%fm)) THEN
         CALL cp_fm_release(preconditioner_env%fm)
         DEALLOCATE (preconditioner_env%fm)
         NULLIFY (preconditioner_env%fm)
      END IF
      IF (preconditioner_env%owns_complex_fm) THEN
         IF (ASSOCIATED(preconditioner_env%complex_fm)) THEN
            CALL cp_cfm_release(preconditioner_env%complex_fm)
            DEALLOCATE (preconditioner_env%complex_fm)
         END IF
      END IF
      NULLIFY (preconditioner_env%complex_fm)
      preconditioner_env%owns_complex_fm = .TRUE.
      IF (ASSOCIATED(preconditioner_env%base_fm)) THEN
         CALL cp_fm_release(preconditioner_env%base_fm)
         DEALLOCATE (preconditioner_env%base_fm)
         NULLIFY (preconditioner_env%base_fm)
      END IF
      IF (ASSOCIATED(preconditioner_env%dbcsr_matrix)) THEN
         CALL dbcsr_release_p(preconditioner_env%dbcsr_matrix)
      END IF
      IF (ASSOCIATED(preconditioner_env%base_dbcsr_matrix)) THEN
         CALL dbcsr_release_p(preconditioner_env%base_dbcsr_matrix)
      END IF
      IF (ASSOCIATED(preconditioner_env%occ_rotation_dbcsr)) THEN
         CALL dbcsr_release_p(preconditioner_env%occ_rotation_dbcsr)
      END IF
      IF (ASSOCIATED(preconditioner_env%occ_rotation)) THEN
         CALL cp_fm_release(preconditioner_env%occ_rotation)
         DEALLOCATE (preconditioner_env%occ_rotation)
         NULLIFY (preconditioner_env%occ_rotation)
      END IF
      IF (ASSOCIATED(preconditioner_env%max_ev_vector)) THEN
         CALL dbcsr_release_p(preconditioner_env%max_ev_vector)
      END IF
      IF (ASSOCIATED(preconditioner_env%min_ev_vector)) THEN
         CALL dbcsr_release_p(preconditioner_env%min_ev_vector)
      END IF
      IF (ASSOCIATED(preconditioner_env%occ_evals)) THEN
         DEALLOCATE (preconditioner_env%occ_evals)
      END IF
      IF (ASSOCIATED(preconditioner_env%full_evals)) THEN
         DEALLOCATE (preconditioner_env%full_evals)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_inverse_k)) THEN
         DEALLOCATE (preconditioner_env%lattice_inverse_k)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_state_inverse_k)) THEN
         DEALLOCATE (preconditioner_env%lattice_state_inverse_k)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_occ_vectors)) THEN
         DEALLOCATE (preconditioner_env%lattice_occ_vectors)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_occ_dual)) THEN
         DEALLOCATE (preconditioner_env%lattice_occ_dual)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_local_indices)) THEN
         DEALLOCATE (preconditioner_env%lattice_local_indices)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_local_operator)) THEN
         DEALLOCATE (preconditioner_env%lattice_local_operator)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_local_coarse_inverse)) THEN
         DEALLOCATE (preconditioner_env%lattice_local_coarse_inverse)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_state_local_indices)) THEN
         DEALLOCATE (preconditioner_env%lattice_state_local_indices)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_state_local_operator)) THEN
         DEALLOCATE (preconditioner_env%lattice_state_local_operator)
      END IF
      IF (ALLOCATED(preconditioner_env%lattice_state_local_coarse_inverse)) THEN
         DEALLOCATE (preconditioner_env%lattice_state_local_coarse_inverse)
      END IF
      CALL release_lattice_rhs_distribution(preconditioner_env%lattice_rhs_distribution)
      IF (ASSOCIATED(preconditioner_env%inverse_history)) THEN
         DO i = 1, SIZE(preconditioner_env%inverse_history)
            CALL dbcsr_release_p(preconditioner_env%inverse_history(i)%matrix)
         END DO
         DEALLOCATE (preconditioner_env%inverse_history)
      END IF
      CALL mp_para_env_release(preconditioner_env%para_env)
      CALL cp_blacs_env_release(preconditioner_env%ctxt)

      preconditioner_env%in_use = 0
      preconditioner_env%cholesky_use = cholesky_reduce

      CALL timestop(handle)

   END SUBROUTINE destroy_preconditioner

! **************************************************************************************************
!> \brief Release a cached lattice right-hand-side communication map.
!> \param distribution ...
! **************************************************************************************************
   SUBROUTINE release_lattice_rhs_distribution(distribution)

      TYPE(lattice_rhs_distribution_type), INTENT(INOUT) :: distribution

      IF (ALLOCATED(distribution%global_columns)) DEALLOCATE (distribution%global_columns)
      IF (ALLOCATED(distribution%local_rows)) DEALLOCATE (distribution%local_rows)
      IF (ALLOCATED(distribution%local_columns)) DEALLOCATE (distribution%local_columns)
      IF (ALLOCATED(distribution%send_counts)) DEALLOCATE (distribution%send_counts)
      IF (ALLOCATED(distribution%send_displacements)) DEALLOCATE (distribution%send_displacements)
      IF (ALLOCATED(distribution%recv_counts)) DEALLOCATE (distribution%recv_counts)
      IF (ALLOCATED(distribution%recv_displacements)) DEALLOCATE (distribution%recv_displacements)
      IF (ALLOCATED(distribution%send_local_rows)) DEALLOCATE (distribution%send_local_rows)
      IF (ALLOCATED(distribution%send_local_columns)) DEALLOCATE (distribution%send_local_columns)
      IF (ALLOCATED(distribution%recv_global_rows)) DEALLOCATE (distribution%recv_global_rows)
      IF (ALLOCATED(distribution%recv_global_columns)) DEALLOCATE (distribution%recv_global_columns)
      distribution%nrow_global = 0
      distribution%ncol_global = 0
      distribution%nrow_local = 0
      distribution%ncol_local = 0
      distribution%num_pe = 0

   END SUBROUTINE release_lattice_rhs_distribution

END MODULE preconditioner_types
