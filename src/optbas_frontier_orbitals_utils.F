!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Utilities for frontier-orbital basis optimization
!> \par History
!>      28.08.2026 created [Jan Wilhelm]
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE optbas_frontier_orbitals_utils
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
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
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_condnum,                      ONLY: overlap_condnum
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_frontier_orbitals_utils'

   TYPE, PUBLIC :: frontier_orbitals_result_type
      REAL(KIND=dp) :: loss_occupied_subspace = 0.0_dp
      REAL(KIND=dp) :: loss_virtual_subspace = 0.0_dp
      REAL(KIND=dp) :: loss_empty_subspace = 0.0_dp
      REAL(KIND=dp) :: loss_gap = 0.0_dp
      REAL(KIND=dp) :: occupied_subspace_overlap = 0.0_dp
      REAL(KIND=dp) :: virtual_subspace_similarity = 0.0_dp
      REAL(KIND=dp) :: empty_subspace_overlap = 0.0_dp
      REAL(KIND=dp) :: gap_reference = 0.0_dp
      REAL(KIND=dp) :: gap_candidate = 0.0_dp
      REAL(KIND=dp) :: condition_number = 0.0_dp
      REAL(KIND=dp) :: reference_condition_number = 0.0_dp
      REAL(KIND=dp) :: calculation_time = 0.0_dp
      INTEGER :: number_candidate_orbitals = 0
      INTEGER :: number_reference_orbitals = 0
   END TYPE frontier_orbitals_result_type

   TYPE, PUBLIC :: frontier_orbitals_reference_type
      INTEGER :: number_occupied_orbitals = 0
      INTEGER :: number_reference_orbitals = 0
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: eigenvalues
      TYPE(cp_fm_type) :: eigenvectors
      REAL(KIND=dp) :: condition_number = 0.0_dp
      INTEGER :: number_kinds = 0
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: element_symbol
      CHARACTER(LEN=default_string_length), DIMENSION(:), ALLOCATABLE :: reference_basis_name
      CHARACTER(LEN=default_string_length), DIMENSION(:), ALLOCATABLE :: candidate_basis_name
   END TYPE frontier_orbitals_reference_type

   PUBLIC :: compute_frontier_orbital_loss_function, &
             evaluate_frontier_orbitals_objective, &
             frontier_orbitals_reference_init, frontier_orbitals_reference_release, &
             print_frontier_orbital_basis_sets, &
             smooth_frontier_orbital_weight

CONTAINS

! **************************************************************************************************
!> \brief Print reference and initial basis-set names for every fitted atom kind.
!> \param unit_nr output unit
!> \param element_symbols fitted atom-kind symbols
!> \param reference_basis_names reference basis-set names
!> \param initial_basis_names initial small-basis names
! **************************************************************************************************
   SUBROUTINE print_frontier_orbital_basis_sets(unit_nr, element_symbols, reference_basis_names, &
                                                initial_basis_names)
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)         :: element_symbols, reference_basis_names, &
                                                            initial_basis_names

      INTEGER                                            :: ikind

      CPASSERT(SIZE(element_symbols) == SIZE(reference_basis_names))
      CPASSERT(SIZE(element_symbols) == SIZE(initial_basis_names))
      WRITE (unit_nr, '(1X,A,T14,A,T51,A)') "Atom kind", "Ref. basis set", "Initial basis set"
      DO ikind = 1, SIZE(element_symbols)
         WRITE (unit_nr, '(1X,A,T14,A,T51,A)') TRIM(element_symbols(ikind)), &
            TRIM(reference_basis_names(ikind)), TRIM(initial_basis_names(ikind))
      END DO

   END SUBROUTINE print_frontier_orbital_basis_sets

! **************************************************************************************************
!> \brief Diagonalize and store the frozen reference Hamiltonian.
!> \param qs_env initialized reference Quickstep environment
!> \param reference reference eigenvalues and eigenvectors
! **************************************************************************************************
   SUBROUTINE frontier_orbitals_reference_init(qs_env, reference)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(frontier_orbitals_reference_type), &
         INTENT(OUT)                                     :: reference

      INTEGER                                            :: nspins
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      LOGICAL                                            :: do_kpoints
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: reference_hamiltonian, &
                                                            reference_overlap, work
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: smat
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (blacs_env, dft_control, fm_struct, matrix_ks, matrix_s, smat)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, dft_control=dft_control, do_kpoints=do_kpoints, &
                      matrix_ks=matrix_ks, matrix_s=matrix_s, nelectron_spin=nelectron_spin)

      nspins = SIZE(matrix_ks)
      IF (nspins /= 1 .OR. dft_control%multiplicity /= 1) THEN
         CPABORT("FRONTIER_ORBITALS currently supports only closed-shell calculations")
      END IF
      IF (do_kpoints) THEN
         CPABORT("FRONTIER_ORBITALS currently supports only Γ-point calculations")
      END IF
      IF (ASSOCIATED(qs_env%x_data)) THEN
         CPABORT("FRONTIER_ORBITALS does not yet support hybrid functionals")
      END IF

      CALL dbcsr_get_info(matrix_s(1)%matrix, &
                          nfullrows_total=reference%number_reference_orbitals)
      ALLOCATE (smat(1, 1))
      smat(1, 1)%matrix => matrix_s(1)%matrix
      CALL overlap_condnum(smat, condnum, -1, .FALSE., .TRUE., .FALSE., blacs_env)
      reference%condition_number = condnum(2)
      DEALLOCATE (smat)

      CALL store_basis_metadata(qs_env, reference)

      ! A restricted spatial orbital contains two electrons.
      reference%number_occupied_orbitals = nelectron_spin(1)/2
      IF (2*reference%number_occupied_orbitals /= nelectron_spin(1)) THEN
         CPABORT("FRONTIER_ORBITALS requires an even closed-shell electron count")
      END IF
      IF (reference%number_occupied_orbitals < 1 .OR. &
          reference%number_occupied_orbitals >= reference%number_reference_orbitals) THEN
         CPABORT("FRONTIER_ORBITALS requires occupied and virtual reference orbitals")
      END IF

      ALLOCATE (reference%eigenvalues(reference%number_reference_orbitals))
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, &
                               nrow_global=reference%number_reference_orbitals, &
                               ncol_global=reference%number_reference_orbitals)
      CALL cp_fm_create(reference_hamiltonian, fm_struct, &
                        name="frontier_orbitals_reference_hamiltonian")
      CALL cp_fm_create(reference_overlap, fm_struct, &
                        name="frontier_orbitals_reference_overlap")
      CALL cp_fm_create(reference%eigenvectors, fm_struct, &
                        name="frontier_orbitals_reference_eigenvectors")
      CALL cp_fm_create(work, fm_struct, name="frontier_orbitals_reference_work")
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix_ks(1)%matrix, reference_hamiltonian)
      CALL cp_fm_uplo_to_full(reference_hamiltonian, work)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, reference_overlap)
      CALL cp_fm_uplo_to_full(reference_overlap, work)

      ! H_R C_R = S_R C_R diag(ε_R), with C_R^T S_R C_R = I.
      CALL cp_fm_geeig(reference_hamiltonian, reference_overlap, reference%eigenvectors, &
                       reference%eigenvalues, work)

      CALL cp_fm_release(reference_hamiltonian)
      CALL cp_fm_release(reference_overlap)
      CALL cp_fm_release(work)

   END SUBROUTINE frontier_orbitals_reference_init

! **************************************************************************************************
!> \brief Release the stored reference eigenvalues and eigenvectors.
!> \param reference reference eigenvalues and eigenvectors
! **************************************************************************************************
   SUBROUTINE frontier_orbitals_reference_release(reference)
      TYPE(frontier_orbitals_reference_type), &
         INTENT(INOUT)                                   :: reference

      IF (ALLOCATED(reference%eigenvalues)) DEALLOCATE (reference%eigenvalues)
      IF (ALLOCATED(reference%element_symbol)) DEALLOCATE (reference%element_symbol)
      IF (ALLOCATED(reference%reference_basis_name)) DEALLOCATE (reference%reference_basis_name)
      IF (ALLOCATED(reference%candidate_basis_name)) DEALLOCATE (reference%candidate_basis_name)
      IF (ASSOCIATED(reference%eigenvectors%matrix_struct)) THEN
         CALL cp_fm_release(reference%eigenvectors)
      END IF
      reference%number_occupied_orbitals = 0
      reference%number_reference_orbitals = 0
      reference%number_kinds = 0
      reference%condition_number = 0.0_dp

   END SUBROUTINE frontier_orbitals_reference_release

! **************************************************************************************************
!> \brief Store the reference and candidate basis assignment for every atomic kind.
!> \param qs_env initialized Quickstep environment
!> \param reference reference data receiving the metadata
! **************************************************************************************************
   SUBROUTINE store_basis_metadata(qs_env, reference)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(frontier_orbitals_reference_type), &
         INTENT(INOUT)                                   :: reference

      INTEGER                                            :: ikind
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis_set, qs_kind_set)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      reference%number_kinds = SIZE(qs_kind_set)
      ALLOCATE (reference%element_symbol(reference%number_kinds))
      ALLOCATE (reference%reference_basis_name(reference%number_kinds))
      ALLOCATE (reference%candidate_basis_name(reference%number_kinds))
      reference%reference_basis_name = ""
      reference%candidate_basis_name = ""

      DO ikind = 1, reference%number_kinds
         CALL get_qs_kind(qs_kind_set(ikind), element_symbol=reference%element_symbol(ikind))
         NULLIFY (basis_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, basis_type="ORB")
         IF (ASSOCIATED(basis_set)) THEN
            CALL get_gto_basis_set(basis_set, name=reference%reference_basis_name(ikind))
         END IF
         NULLIFY (basis_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_set, basis_type="AUX_OPT")
         IF (ASSOCIATED(basis_set)) THEN
            CALL get_gto_basis_set(basis_set, name=reference%candidate_basis_name(ikind))
         END IF
      END DO

   END SUBROUTINE store_basis_metadata

! **************************************************************************************************
!> \brief Construct and diagonalize the candidate Hamiltonian and evaluate the loss function.
!> \param reference reference eigenvalues and eigenvectors C_R
!> \param matrix_s_candidate candidate AO overlap matrix S_C
!> \param matrix_s_candidate_reference candidate/reference AO overlap matrix S_CR
!> \param virtual_cutoff virtual-orbital energy cutoff
!> \param virtual_smoothing virtual-orbital energy smoothing width
!> \param gap_scale normalization energy for the HOMO-LUMO gap loss
!> \param objective_result loss-function terms and result information
!> \note The reference Hamiltonian is transformed to the candidate AO basis using
!> B = C_R^T S_RC
!> H_C = B^T diag(ε_R) B
!> The candidate orbitals are obtained from
!> H_C C_C = S_C C_C diag(ε_C)
!> C_C^T S_C C_C = I. The reference/candidate orbital overlap matrix is
!> M = C_R^T S_RC C_C = B C_C.
! **************************************************************************************************
   SUBROUTINE evaluate_frontier_orbitals_objective(reference, matrix_s_candidate, &
                                                   matrix_s_candidate_reference, virtual_cutoff, &
                                                   virtual_smoothing, gap_scale, objective_result)
      TYPE(frontier_orbitals_reference_type), INTENT(IN) :: reference
      TYPE(dbcsr_type), POINTER                          :: matrix_s_candidate, &
                                                            matrix_s_candidate_reference
      REAL(KIND=dp), INTENT(IN)                          :: virtual_cutoff, virtual_smoothing, &
                                                            gap_scale
      TYPE(frontier_orbitals_result_type), INTENT(OUT)   :: objective_result

      INTEGER                                            :: istate, naux, nref
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_candidate
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: orbital_overlap_dense
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: ao_cross_struct, aux_struct, &
                                                            cross_struct, ref_struct
      TYPE(cp_fm_type) :: candidate_eigenvectors, candidate_hamiltonian, candidate_overlap, &
         candidate_work, reference_candidate_ao_overlap, reference_candidate_orbital_overlap, &
         reference_eigenvalue_matrix, reference_orbital_candidate_ao_overlap
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: smat

      NULLIFY (ao_cross_struct, aux_struct, blacs_env, cross_struct, ref_struct, smat)
      nref = reference%number_reference_orbitals
      CALL dbcsr_get_info(matrix_s_candidate, nfullrows_total=naux)
      IF (naux <= reference%number_occupied_orbitals) THEN
         CPABORT("FRONTIER_ORBITALS candidate basis has no virtual orbitals")
      END IF

      blacs_env => reference%eigenvectors%matrix_struct%context
      ALLOCATE (smat(1, 1))
      smat(1, 1)%matrix => matrix_s_candidate
      CALL overlap_condnum(smat, condnum, -1, .FALSE., .TRUE., .FALSE., blacs_env)
      DEALLOCATE (smat)
      objective_result%condition_number = condnum(2)
      objective_result%number_candidate_orbitals = naux
      IF (.NOT. ieee_is_finite(objective_result%condition_number) .OR. &
          objective_result%condition_number > 1.0E14_dp .OR. &
          objective_result%condition_number <= 0.0_dp) THEN
         objective_result%loss_occupied_subspace = 1.0E6_dp
         objective_result%loss_virtual_subspace = 1.0E6_dp
         objective_result%loss_empty_subspace = 1.0E6_dp
         objective_result%loss_gap = 1.0E6_dp
         RETURN
      END IF

      CALL cp_fm_struct_create(aux_struct, context=blacs_env, nrow_global=naux, ncol_global=naux)
      CALL cp_fm_struct_create(ao_cross_struct, context=blacs_env, &
                               nrow_global=naux, ncol_global=nref)
      CALL cp_fm_struct_create(cross_struct, context=blacs_env, nrow_global=nref, ncol_global=naux)
      CALL cp_fm_struct_create(ref_struct, context=blacs_env, nrow_global=nref, ncol_global=nref)

      CALL cp_fm_create(candidate_hamiltonian, aux_struct, &
                        name="frontier_orbitals_candidate_hamiltonian")
      CALL cp_fm_create(candidate_overlap, aux_struct, &
                        name="frontier_orbitals_candidate_overlap")
      CALL cp_fm_create(candidate_eigenvectors, aux_struct, &
                        name="frontier_orbitals_candidate_eigenvectors")
      CALL cp_fm_create(candidate_work, aux_struct, name="frontier_orbitals_candidate_work")
      CALL cp_fm_create(reference_candidate_ao_overlap, matrix_struct=ao_cross_struct, &
                        name="frontier_orbitals_reference_candidate_ao_overlap")
      CALL cp_fm_create(reference_orbital_candidate_ao_overlap, matrix_struct=cross_struct, &
                        name="frontier_orbitals_reference_orbital_candidate_ao_overlap")
      CALL cp_fm_create(reference_candidate_orbital_overlap, matrix_struct=cross_struct, &
                        name="frontier_orbitals_reference_candidate_orbital_overlap")
      CALL cp_fm_create(reference_eigenvalue_matrix, matrix_struct=ref_struct, &
                        name="frontier_orbitals_reference_eigenvalues")
      CALL cp_fm_struct_release(ao_cross_struct)
      CALL cp_fm_struct_release(aux_struct)
      CALL cp_fm_struct_release(cross_struct)
      CALL cp_fm_struct_release(ref_struct)

      ! B = C_R^T S_RC. The supplied cross-overlap matrix is stored as S_CR.
      CALL copy_dbcsr_to_fm(matrix_s_candidate_reference, reference_candidate_ao_overlap)
      CALL parallel_gemm('T', 'T', nref, naux, nref, 1.0_dp, reference%eigenvectors, &
                         reference_candidate_ao_overlap, 0.0_dp, &
                         reference_orbital_candidate_ao_overlap)

      ! H_C = B^T diag(ε_R) B.
      CALL cp_fm_set_all(reference_eigenvalue_matrix, 0.0_dp)
      DO istate = 1, nref
         CALL cp_fm_set_element(reference_eigenvalue_matrix, istate, istate, &
                                reference%eigenvalues(istate))
      END DO
      CALL parallel_gemm('N', 'N', nref, naux, nref, 1.0_dp, reference_eigenvalue_matrix, &
                         reference_orbital_candidate_ao_overlap, 0.0_dp, &
                         reference_candidate_orbital_overlap)
      CALL parallel_gemm('T', 'N', naux, naux, nref, 1.0_dp, &
                         reference_orbital_candidate_ao_overlap, &
                         reference_candidate_orbital_overlap, 0.0_dp, candidate_hamiltonian)

      CALL copy_dbcsr_to_fm(matrix_s_candidate, candidate_overlap)
      CALL cp_fm_uplo_to_full(candidate_overlap, candidate_work)
      ALLOCATE (eigenvalues_candidate(naux))

      ! H_C C_C = S_C C_C diag(ε_C).
      CALL cp_fm_geeig(candidate_hamiltonian, candidate_overlap, candidate_eigenvectors, &
                       eigenvalues_candidate, candidate_work)

      ! M = B C_C, where M_ij = <ψ_i^R | ψ_j^C>.
      CALL parallel_gemm('N', 'N', nref, naux, naux, 1.0_dp, &
                         reference_orbital_candidate_ao_overlap, candidate_eigenvectors, &
                         0.0_dp, reference_candidate_orbital_overlap)
      ALLOCATE (orbital_overlap_dense(nref, naux))
      CALL cp_fm_get_submatrix(reference_candidate_orbital_overlap, orbital_overlap_dense)
      CALL compute_frontier_orbital_loss_function( &
         orbital_overlap_dense, reference%eigenvalues, eigenvalues_candidate, &
         reference%number_occupied_orbitals, virtual_cutoff, virtual_smoothing, gap_scale, objective_result)
      objective_result%condition_number = condnum(2)
      objective_result%number_candidate_orbitals = naux

      DEALLOCATE (eigenvalues_candidate, orbital_overlap_dense)
      CALL cp_fm_release(candidate_eigenvectors)
      CALL cp_fm_release(candidate_hamiltonian)
      CALL cp_fm_release(candidate_overlap)
      CALL cp_fm_release(candidate_work)
      CALL cp_fm_release(reference_candidate_ao_overlap)
      CALL cp_fm_release(reference_candidate_orbital_overlap)
      CALL cp_fm_release(reference_eigenvalue_matrix)
      CALL cp_fm_release(reference_orbital_candidate_ao_overlap)

   END SUBROUTINE evaluate_frontier_orbitals_objective

! **************************************************************************************************
!> \brief Compute the frontier-orbital loss function.
!> \param orbital_overlap reference/candidate orbital overlaps
!> \param eigenvalues_reference reference eigenvalues
!> \param eigenvalues_candidate candidate eigenvalues
!> \param number_occupied_orbitals number of occupied spatial orbitals
!> \param virtual_cutoff virtual-orbital energy cutoff E_cut
!> \param virtual_smoothing virtual-orbital energy smoothing width ΔE
!> \param gap_scale normalization energy E_scale for the gap loss
!> \param objective_result loss-function terms and result information
!> \note Let i,j denote occupied orbitals, a,b virtual orbitals, and N_occ the number of occupied
!> orbitals. R and C label reference and candidate quantities, with
!> M_ij = <ψ_i^R | ψ_j^C> and M_ab = <ψ_a^R | ψ_b^C>. The occupied-subspace loss is
!> L_occ = 1 - sum_(i,j in occupied) |M_ij|^2 / N_occ.
!>
!> For virtual orbitals, define the smooth energy-window weights and their norms as
!> w_a^X = 1 / {1 + exp[(ε_a^X - ε_LUMO^X - E_cut) / ΔE]}, for X = R or C,
!> N_X = sum_(a in virtual) (w_a^X)^2.
!> The virtual-subspace loss is the normalized squared Frobenius distance between the
!> energy-weighted virtual-space operators. Its expanded form is
!> L_vir = [N_R + N_C - 2 sum_(a,b in virtual) w_a^R w_b^C |M_ab|^2] / (2 N_R).
!>
!> The empty-subspace loss is
!> L_empty = 1 - sum_(a,b in virtual) (w_a^R)^2 |M_ab|^2 / N_R.
!>
!> With E_gap^X = ε_LUMO^X - ε_HOMO^X, the gap loss is
!> L_gap = [(E_gap^C - E_gap^R) / E_scale]^2.
! **************************************************************************************************
   PURE SUBROUTINE compute_frontier_orbital_loss_function( &
      orbital_overlap, eigenvalues_reference, eigenvalues_candidate, number_occupied_orbitals, &
      virtual_cutoff, virtual_smoothing, gap_scale, objective_result)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: orbital_overlap
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues_reference, &
                                                            eigenvalues_candidate
      INTEGER, INTENT(IN)                                :: number_occupied_orbitals
      REAL(KIND=dp), INTENT(IN)                          :: virtual_cutoff, virtual_smoothing, &
                                                            gap_scale
      TYPE(frontier_orbitals_result_type), INTENT(OUT)   :: objective_result

      INTEGER                                            :: candidate_orbital, reference_orbital
      REAL(KIND=dp) :: candidate_virtual_weight_norm, occupied_subspace_overlap, &
         reference_virtual_subspace_overlap, reference_virtual_weight_norm, &
         virtual_subspace_overlap, weight_candidate, weight_reference

      ! O_occ = sum_{i,j in occupied} |M_ij|^2 / N_occ.
      occupied_subspace_overlap = &
         SUM(orbital_overlap(1:number_occupied_orbitals, &
                             1:number_occupied_orbitals)**2)/ &
         REAL(number_occupied_orbitals, KIND=dp)
      objective_result%loss_occupied_subspace = MAX(0.0_dp, 1.0_dp - occupied_subspace_overlap)
      objective_result%occupied_subspace_overlap = occupied_subspace_overlap

      ! N_R = sum_{a in virtual} (w_a^R)^2.
      reference_virtual_weight_norm = 0.0_dp
      DO reference_orbital = number_occupied_orbitals + 1, SIZE(eigenvalues_reference)
         weight_reference = smooth_frontier_orbital_weight( &
                            eigenvalues_reference(reference_orbital) - &
                            eigenvalues_reference(number_occupied_orbitals + 1), &
                            virtual_cutoff, virtual_smoothing)
         reference_virtual_weight_norm = reference_virtual_weight_norm + weight_reference**2
      END DO

      ! N_C = sum_{b in virtual} (w_b^C)^2.
      candidate_virtual_weight_norm = 0.0_dp
      DO candidate_orbital = number_occupied_orbitals + 1, SIZE(eigenvalues_candidate)
         weight_candidate = smooth_frontier_orbital_weight( &
                            eigenvalues_candidate(candidate_orbital) - &
                            eigenvalues_candidate(number_occupied_orbitals + 1), &
                            virtual_cutoff, virtual_smoothing)
         candidate_virtual_weight_norm = candidate_virtual_weight_norm + weight_candidate**2
      END DO

      ! O_vir = sum_{a,b in virtual} w_a^R w_b^C |M_ab|^2.
      virtual_subspace_overlap = 0.0_dp
      reference_virtual_subspace_overlap = 0.0_dp
      DO candidate_orbital = number_occupied_orbitals + 1, SIZE(eigenvalues_candidate)
         weight_candidate = smooth_frontier_orbital_weight( &
                            eigenvalues_candidate(candidate_orbital) - &
                            eigenvalues_candidate(number_occupied_orbitals + 1), &
                            virtual_cutoff, virtual_smoothing)
         DO reference_orbital = number_occupied_orbitals + 1, SIZE(eigenvalues_reference)
            weight_reference = smooth_frontier_orbital_weight( &
                               eigenvalues_reference(reference_orbital) - &
                               eigenvalues_reference(number_occupied_orbitals + 1), &
                               virtual_cutoff, virtual_smoothing)
            virtual_subspace_overlap = virtual_subspace_overlap + &
                                       weight_reference*weight_candidate* &
                                       orbital_overlap(reference_orbital, candidate_orbital)**2
            reference_virtual_subspace_overlap = reference_virtual_subspace_overlap + &
                                                 weight_reference**2* &
                                                 orbital_overlap(reference_orbital, candidate_orbital)**2
         END DO
      END DO

      IF (reference_virtual_weight_norm > TINY(1.0_dp)) THEN
         ! L_vir = (N_R + N_C - 2 O_vir)/(2 N_R).
         objective_result%loss_virtual_subspace = &
            MAX(0.0_dp, (reference_virtual_weight_norm + candidate_virtual_weight_norm - &
                         2.0_dp*virtual_subspace_overlap)/(2.0_dp*reference_virtual_weight_norm))
         ! Bound the computed overlap against small roundoff excursions from its exact range [0, 1].
         objective_result%empty_subspace_overlap = &
            MIN(1.0_dp, MAX(0.0_dp, reference_virtual_subspace_overlap/reference_virtual_weight_norm))
         objective_result%loss_empty_subspace = 1.0_dp - objective_result%empty_subspace_overlap
         IF (candidate_virtual_weight_norm > TINY(1.0_dp)) THEN
            objective_result%virtual_subspace_similarity = &
               MIN(1.0_dp, MAX(0.0_dp, &
                               virtual_subspace_overlap/ &
                               SQRT(reference_virtual_weight_norm*candidate_virtual_weight_norm)))
         END IF
      END IF

      ! E_gap^X = ε_LUMO^X - ε_HOMO^X.
      objective_result%gap_reference = eigenvalues_reference(number_occupied_orbitals + 1) - &
                                       eigenvalues_reference(number_occupied_orbitals)
      objective_result%gap_candidate = eigenvalues_candidate(number_occupied_orbitals + 1) - &
                                       eigenvalues_candidate(number_occupied_orbitals)

      ! L_gap = [(E_gap^C - E_gap^R)/E_scale]^2.
      objective_result%loss_gap = &
         ((objective_result%gap_candidate - objective_result%gap_reference)/gap_scale)**2

   END SUBROUTINE compute_frontier_orbital_loss_function

! **************************************************************************************************
!> \brief Return the smooth energy-window weight for a virtual orbital.
!> \param energy_from_lumo orbital energy relative to the LUMO
!> \param cutoff outer energy of the virtual-orbital window
!> \param smoothing smoothing width at the outer boundary
!> \return energy-window weight in the interval [0, 1]
! **************************************************************************************************
   PURE FUNCTION smooth_frontier_orbital_weight(energy_from_lumo, cutoff, smoothing) RESULT(weight)
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

   END FUNCTION smooth_frontier_orbital_weight

END MODULE optbas_frontier_orbitals_utils
