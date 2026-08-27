!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE optbas_frontier_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_uplo_to_full
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_element,&
                                              cp_fm_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_condnum,                      ONLY: overlap_condnum
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_frontier_utils'

   TYPE, PUBLIC :: frontier_result_type
      REAL(KIND=dp) :: loss_occupied = 0.0_dp
      REAL(KIND=dp) :: loss_virtual = 0.0_dp
      REAL(KIND=dp) :: loss_gap = 0.0_dp
      REAL(KIND=dp) :: loss_coefficients = 0.0_dp
      REAL(KIND=dp) :: occupied_overlap_percent = 0.0_dp
      REAL(KIND=dp) :: virtual_overlap_percent = 0.0_dp
      REAL(KIND=dp) :: gap_reference = 0.0_dp
      REAL(KIND=dp) :: gap_candidate = 0.0_dp
      REAL(KIND=dp) :: condition_number = 0.0_dp
   END TYPE frontier_result_type

   TYPE, PUBLIC :: frontier_reference_type
      INTEGER :: nocc = 0
      INTEGER :: nstates = 0
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: eigenvalues
      TYPE(cp_fm_type) :: eigenvectors
   END TYPE frontier_reference_type

   PUBLIC :: closed_shell_nocc, compute_frontier_losses, evaluate_frontier_objective, &
             frontier_reference_init, frontier_reference_release, smooth_frontier_weight

CONTAINS

! **************************************************************************************************
!> \brief Convert a restricted-channel electron count to occupied spatial orbitals.
!> \param nelectron even closed-shell electron count
!> \return number of doubly occupied spatial orbitals
! **************************************************************************************************
   PURE FUNCTION closed_shell_nocc(nelectron) RESULT(nocc)
      INTEGER, INTENT(IN)                                :: nelectron
      INTEGER                                            :: nocc

      nocc = nelectron/2

   END FUNCTION closed_shell_nocc

! **************************************************************************************************
!> \brief Numerically stable Fermi-like weight for the virtual frontier window.
!> \param energy_from_lumo orbital energy relative to the LUMO
!> \param cutoff outer energy of the target window
!> \param smoothing smoothing width at the outer boundary
!> \return frontier weight in [0, 1]
! **************************************************************************************************
   PURE FUNCTION smooth_frontier_weight(energy_from_lumo, cutoff, smoothing) RESULT(weight)
      REAL(KIND=dp), INTENT(IN)                          :: energy_from_lumo, cutoff, smoothing
      REAL(KIND=dp)                                      :: weight

      REAL(KIND=dp)                                      :: argument

      argument = (energy_from_lumo - cutoff)/smoothing
      IF (argument >= 40.0_dp) THEN
         weight = 0.0_dp
      ELSE IF (argument <= -40.0_dp) THEN
         weight = 1.0_dp
      ELSE
         weight = 1.0_dp/(1.0_dp + EXP(argument))
      END IF

   END FUNCTION smooth_frontier_weight

! **************************************************************************************************
!> \brief Evaluate the normalized occupied, virtual, gap, and coefficient losses.
!> \param overlap C_R^T S_RC C_C for all reference and candidate eigenstates
!> \param eigenvalues_reference frozen-reference eigenvalues in ascending order
!> \param eigenvalues_candidate candidate eigenvalues in ascending order
!> \param nocc number of occupied orbitals
!> \param virtual_cutoff virtual-window cutoff in atomic units
!> \param virtual_smoothing virtual-window smoothing width in atomic units
!> \param gap_scale gap normalization energy in atomic units
!> \param coefficients current ordered optimization coefficient vector
!> \param coefficients_initial initial ordered optimization coefficient vector
!> \param result individual normalized losses and reporting quantities
! **************************************************************************************************
   PURE SUBROUTINE compute_frontier_losses(overlap, eigenvalues_reference, eigenvalues_candidate, &
                                           nocc, virtual_cutoff, virtual_smoothing, gap_scale, &
                                           coefficients, coefficients_initial, RESULT)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: overlap
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues_reference, &
                                                            eigenvalues_candidate
      INTEGER, INTENT(IN)                                :: nocc
      REAL(KIND=dp), INTENT(IN)                          :: virtual_cutoff, virtual_smoothing, &
                                                            gap_scale
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coefficients, coefficients_initial
      TYPE(frontier_result_type), INTENT(OUT)            :: result

      INTEGER                                            :: icand, iref
      REAL(KIND=dp)                                      :: candidate_norm, common_content, &
                                                            occupied_content, reference_norm, &
                                                            w_candidate, w_reference

      RESULT = frontier_result_type()

      occupied_content = SUM(overlap(1:nocc, 1:nocc)**2)
      result%loss_occupied = MAX(0.0_dp, 1.0_dp - occupied_content/REAL(nocc, KIND=dp))
      result%occupied_overlap_percent = 100.0_dp*occupied_content/REAL(nocc, KIND=dp)

      reference_norm = 0.0_dp
      DO iref = nocc + 1, SIZE(eigenvalues_reference)
         w_reference = smooth_frontier_weight(eigenvalues_reference(iref) - &
                                              eigenvalues_reference(nocc + 1), &
                                              virtual_cutoff, virtual_smoothing)
         reference_norm = reference_norm + w_reference**2
      END DO

      candidate_norm = 0.0_dp
      DO icand = nocc + 1, SIZE(eigenvalues_candidate)
         w_candidate = smooth_frontier_weight(eigenvalues_candidate(icand) - &
                                              eigenvalues_candidate(nocc + 1), &
                                              virtual_cutoff, virtual_smoothing)
         candidate_norm = candidate_norm + w_candidate**2
      END DO

      common_content = 0.0_dp
      DO icand = nocc + 1, SIZE(eigenvalues_candidate)
         w_candidate = smooth_frontier_weight(eigenvalues_candidate(icand) - &
                                              eigenvalues_candidate(nocc + 1), &
                                              virtual_cutoff, virtual_smoothing)
         DO iref = nocc + 1, SIZE(eigenvalues_reference)
            w_reference = smooth_frontier_weight(eigenvalues_reference(iref) - &
                                                 eigenvalues_reference(nocc + 1), &
                                                 virtual_cutoff, virtual_smoothing)
            common_content = common_content + w_reference*w_candidate*overlap(iref, icand)**2
         END DO
      END DO

      IF (reference_norm > TINY(1.0_dp)) THEN
         result%loss_virtual = MAX(0.0_dp, &
                                   (reference_norm + candidate_norm - 2.0_dp*common_content)/(2.0_dp*reference_norm))
         result%virtual_overlap_percent = 100.0_dp*common_content/reference_norm
      END IF

      result%gap_reference = eigenvalues_reference(nocc + 1) - eigenvalues_reference(nocc)
      result%gap_candidate = eigenvalues_candidate(nocc + 1) - eigenvalues_candidate(nocc)
      result%loss_gap = ((result%gap_candidate - result%gap_reference)/gap_scale)**2

      IF (SIZE(coefficients_initial) > 0 .AND. SUM(coefficients_initial**2) > TINY(1.0_dp)) THEN
         result%loss_coefficients = SUM((coefficients - coefficients_initial)**2)/ &
                                    SUM(coefficients_initial**2)
      END IF

   END SUBROUTINE compute_frontier_losses

! **************************************************************************************************
!> \brief Diagonalize and cache the large frozen-reference eigenproblem.
!> \param qs_env initialized reference Quickstep environment
!> \param reference cached reference eigenvalues and eigenvectors
! **************************************************************************************************
   SUBROUTINE frontier_reference_init(qs_env, reference)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(frontier_reference_type), INTENT(OUT)         :: reference

      INTEGER                                            :: nspins
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      LOGICAL                                            :: do_kpoints
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: hamiltonian, overlap, work
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (blacs_env, dft_control, fm_struct, matrix_ks, matrix_s)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, dft_control=dft_control, do_kpoints=do_kpoints, &
                      matrix_ks=matrix_ks, matrix_s=matrix_s, nelectron_spin=nelectron_spin)

      nspins = SIZE(matrix_ks)
      IF (nspins /= 1 .OR. dft_control%multiplicity /= 1) THEN
         CPABORT("FRONTIER_STATES currently supports only closed-shell calculations")
      END IF
      IF (do_kpoints) CPABORT("FRONTIER_STATES currently supports only Gamma-point calculations")
      IF (ASSOCIATED(qs_env%x_data)) THEN
         CPABORT("FRONTIER_STATES does not yet support hybrid functionals")
      END IF

      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=reference%nstates)
      ! In the single restricted-spin channel each spatial orbital holds two electrons.
      reference%nocc = closed_shell_nocc(nelectron_spin(1))
      IF (2*reference%nocc /= nelectron_spin(1)) THEN
         CPABORT("FRONTIER_STATES requires an even closed-shell electron count")
      END IF
      IF (reference%nocc < 1 .OR. reference%nocc >= reference%nstates) THEN
         CPABORT("FRONTIER_STATES requires both occupied and unoccupied reference states")
      END IF

      ALLOCATE (reference%eigenvalues(reference%nstates))
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, &
                               nrow_global=reference%nstates, ncol_global=reference%nstates)
      CALL cp_fm_create(hamiltonian, fm_struct, name="frontier_reference_h")
      CALL cp_fm_create(overlap, fm_struct, name="frontier_reference_s")
      CALL cp_fm_create(reference%eigenvectors, fm_struct, name="frontier_reference_c")
      CALL cp_fm_create(work, fm_struct, name="frontier_reference_work")
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix_ks(1)%matrix, hamiltonian)
      CALL cp_fm_uplo_to_full(hamiltonian, work)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, overlap)
      CALL cp_fm_uplo_to_full(overlap, work)
      CALL cp_fm_geeig(hamiltonian, overlap, reference%eigenvectors, reference%eigenvalues, work)

      CALL cp_fm_release(hamiltonian)
      CALL cp_fm_release(overlap)
      CALL cp_fm_release(work)

   END SUBROUTINE frontier_reference_init

! **************************************************************************************************
!> \brief Release a frozen-reference cache.
!> \param reference cached reference data
! **************************************************************************************************
   SUBROUTINE frontier_reference_release(reference)
      TYPE(frontier_reference_type), INTENT(INOUT)       :: reference

      IF (ALLOCATED(reference%eigenvalues)) DEALLOCATE (reference%eigenvalues)
      IF (ASSOCIATED(reference%eigenvectors%matrix_struct)) CALL cp_fm_release(reference%eigenvectors)
      reference%nocc = 0
      reference%nstates = 0

   END SUBROUTINE frontier_reference_release

! **************************************************************************************************
!> \brief Build and diagonalize the projected frozen candidate Hamiltonian, then evaluate losses.
!> \param reference cached full reference eigensystem
!> \param matrix_s_candidate candidate AO overlap matrix
!> \param matrix_s_candidate_reference candidate/reference AO cross-overlap matrix
!> \param virtual_cutoff virtual window cutoff
!> \param virtual_smoothing virtual window smoothing width
!> \param gap_scale gap normalization energy
!> \param coefficients current ordered coefficient vector
!> \param coefficients_initial initial ordered coefficient vector
!> \param result objective components and reporting values
! **************************************************************************************************
   SUBROUTINE evaluate_frontier_objective(reference, matrix_s_candidate, matrix_s_candidate_reference, &
                                          virtual_cutoff, virtual_smoothing, gap_scale, coefficients, &
                                          coefficients_initial, RESULT)
      TYPE(frontier_reference_type), INTENT(IN)          :: reference
      TYPE(dbcsr_type), POINTER                          :: matrix_s_candidate, &
                                                            matrix_s_candidate_reference
      REAL(KIND=dp), INTENT(IN)                          :: virtual_cutoff, virtual_smoothing, &
                                                            gap_scale
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coefficients, coefficients_initial
      TYPE(frontier_result_type), INTENT(OUT)            :: result

      INTEGER                                            :: istate, naux, nref
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_candidate
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap_dense
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: ao_cross_struct, aux_struct, &
                                                            cross_struct, ref_struct
      TYPE(cp_fm_type) :: candidate_eigenvectors, candidate_hamiltonian, candidate_overlap, &
         candidate_work, cross_ao, cross_eigen, eigenvalue_matrix, overlap_eigen
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: smat

      NULLIFY (ao_cross_struct, aux_struct, blacs_env, cross_struct, ref_struct, smat)
      nref = reference%nstates
      CALL dbcsr_get_info(matrix_s_candidate, nfullrows_total=naux)
      IF (naux <= reference%nocc) THEN
         CPABORT("FRONTIER_STATES candidate basis has no unoccupied states")
      END IF

      blacs_env => reference%eigenvectors%matrix_struct%context
      ALLOCATE (smat(1, 1))
      smat(1, 1)%matrix => matrix_s_candidate
      CALL overlap_condnum(smat, condnum, -1, .FALSE., .TRUE., .FALSE., blacs_env)
      result%condition_number = condnum(2)
      DEALLOCATE (smat)
      IF (.NOT. ieee_is_finite(result%condition_number) .OR. &
          result%condition_number > 1.0E14_dp .OR. result%condition_number <= 0.0_dp) THEN
         result%loss_occupied = 1.0E6_dp
         result%loss_virtual = 1.0E6_dp
         result%loss_gap = 1.0E6_dp
         IF (SIZE(coefficients_initial) > 0 .AND. SUM(coefficients_initial**2) > TINY(1.0_dp)) THEN
            result%loss_coefficients = SUM((coefficients - coefficients_initial)**2)/ &
                                       SUM(coefficients_initial**2)
         END IF
         RETURN
      END IF

      CALL cp_fm_struct_create(aux_struct, context=blacs_env, nrow_global=naux, ncol_global=naux)
      CALL cp_fm_struct_create(ao_cross_struct, context=blacs_env, nrow_global=naux, ncol_global=nref)
      CALL cp_fm_struct_create(cross_struct, context=blacs_env, nrow_global=nref, ncol_global=naux)
      CALL cp_fm_struct_create(ref_struct, context=blacs_env, nrow_global=nref, ncol_global=nref)

      CALL cp_fm_create(candidate_hamiltonian, aux_struct, name="frontier_candidate_h")
      CALL cp_fm_create(candidate_overlap, aux_struct, name="frontier_candidate_s")
      CALL cp_fm_create(candidate_eigenvectors, aux_struct, name="frontier_candidate_c")
      CALL cp_fm_create(candidate_work, aux_struct, name="frontier_candidate_work")
      CALL cp_fm_create(cross_ao, matrix_struct=ao_cross_struct, name="frontier_cross_ao")
      CALL cp_fm_create(cross_eigen, matrix_struct=cross_struct, name="frontier_cross_eigen")
      CALL cp_fm_create(overlap_eigen, matrix_struct=cross_struct, name="frontier_overlap_eigen")
      CALL cp_fm_create(eigenvalue_matrix, matrix_struct=ref_struct, name="frontier_reference_e")
      CALL cp_fm_struct_release(ao_cross_struct)
      CALL cp_fm_struct_release(aux_struct)
      CALL cp_fm_struct_release(cross_struct)
      CALL cp_fm_struct_release(ref_struct)

      ! S_CR is stored candidate-by-reference. Transform it to reference eigenstate-by-candidate AO.
      CALL copy_dbcsr_to_fm(matrix_s_candidate_reference, cross_ao)
      CALL parallel_gemm('T', 'T', nref, naux, nref, 1.0_dp, reference%eigenvectors, cross_ao, &
                         0.0_dp, cross_eigen)

      ! Galerkin representation H_C = (C_R^T S_RC)^T diag(e_R) (C_R^T S_RC).
      CALL cp_fm_set_all(eigenvalue_matrix, 0.0_dp)
      DO istate = 1, nref
         CALL cp_fm_set_element(eigenvalue_matrix, istate, istate, reference%eigenvalues(istate))
      END DO
      CALL parallel_gemm('N', 'N', nref, naux, nref, 1.0_dp, eigenvalue_matrix, cross_eigen, &
                         0.0_dp, overlap_eigen)
      CALL parallel_gemm('T', 'N', naux, naux, nref, 1.0_dp, cross_eigen, overlap_eigen, &
                         0.0_dp, candidate_hamiltonian)

      CALL copy_dbcsr_to_fm(matrix_s_candidate, candidate_overlap)
      CALL cp_fm_uplo_to_full(candidate_overlap, candidate_work)
      ALLOCATE (eigenvalues_candidate(naux))
      CALL cp_fm_geeig(candidate_hamiltonian, candidate_overlap, candidate_eigenvectors, &
                       eigenvalues_candidate, candidate_work)

      CALL parallel_gemm('N', 'N', nref, naux, naux, 1.0_dp, cross_eigen, candidate_eigenvectors, &
                         0.0_dp, overlap_eigen)
      ALLOCATE (overlap_dense(nref, naux))
      CALL cp_fm_get_submatrix(overlap_eigen, overlap_dense)
      CALL compute_frontier_losses(overlap_dense, reference%eigenvalues, eigenvalues_candidate, &
                                   reference%nocc, virtual_cutoff, virtual_smoothing, gap_scale, &
                                   coefficients, coefficients_initial, RESULT)
      result%condition_number = condnum(2)

      DEALLOCATE (eigenvalues_candidate, overlap_dense)
      CALL cp_fm_release(candidate_eigenvectors)
      CALL cp_fm_release(candidate_hamiltonian)
      CALL cp_fm_release(candidate_overlap)
      CALL cp_fm_release(candidate_work)
      CALL cp_fm_release(cross_ao)
      CALL cp_fm_release(cross_eigen)
      CALL cp_fm_release(eigenvalue_matrix)
      CALL cp_fm_release(overlap_eigen)

   END SUBROUTINE evaluate_frontier_objective

END MODULE optbas_frontier_utils
