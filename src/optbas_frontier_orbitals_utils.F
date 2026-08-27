!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE optbas_frontier_orbitals_utils
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optbas_frontier_orbitals_utils'

   TYPE, PUBLIC :: frontier_orbitals_result_type
      REAL(KIND=dp) :: loss_occupied_subspace = 0.0_dp
      REAL(KIND=dp) :: loss_virtual_subspace = 0.0_dp
      REAL(KIND=dp) :: loss_empty_subspace = 0.0_dp
      REAL(KIND=dp) :: loss_gap = 0.0_dp
      REAL(KIND=dp) :: occupied_subspace_overlap = 0.0_dp
      REAL(KIND=dp) :: empty_subspace_overlap = 0.0_dp
      REAL(KIND=dp) :: gap_reference = 0.0_dp
      REAL(KIND=dp) :: gap_candidate = 0.0_dp
      REAL(KIND=dp) :: condition_number = 0.0_dp
   END TYPE frontier_orbitals_result_type

   TYPE, PUBLIC :: frontier_orbitals_reference_type
      INTEGER :: number_occupied_orbitals = 0
      INTEGER :: number_reference_orbitals = 0
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: eigenvalues
      TYPE(cp_fm_type) :: eigenvectors
   END TYPE frontier_orbitals_reference_type

   PUBLIC :: compute_frontier_orbital_loss_function, evaluate_frontier_orbitals_objective, &
             frontier_orbitals_reference_init, frontier_orbitals_reference_release, &
             smooth_frontier_orbital_weight

CONTAINS

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
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: reference_hamiltonian, &
                                                            reference_overlap, work
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (blacs_env, dft_control, fm_struct, matrix_ks, matrix_s)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, dft_control=dft_control, do_kpoints=do_kpoints, &
                      matrix_ks=matrix_ks, matrix_s=matrix_s, nelectron_spin=nelectron_spin)

      nspins = SIZE(matrix_ks)
      IF (nspins /= 1 .OR. dft_control%multiplicity /= 1) THEN
         CPABORT("FRONTIER_ORBITALS currently supports only closed-shell calculations")
      END IF
      IF (do_kpoints) THEN
         CPABORT("FRONTIER_ORBITALS currently supports only Gamma-point calculations")
      END IF
      IF (ASSOCIATED(qs_env%x_data)) THEN
         CPABORT("FRONTIER_ORBITALS does not yet support hybrid functionals")
      END IF

      CALL dbcsr_get_info(matrix_s(1)%matrix, &
                          nfullrows_total=reference%number_reference_orbitals)

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

      ! H_R C_R = S_R C_R diag(epsilon_R), with C_R^T S_R C_R = I.
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
      IF (ASSOCIATED(reference%eigenvectors%matrix_struct)) THEN
         CALL cp_fm_release(reference%eigenvectors)
      END IF
      reference%number_occupied_orbitals = 0
      reference%number_reference_orbitals = 0

   END SUBROUTINE frontier_orbitals_reference_release

! **************************************************************************************************
!> \brief Construct and diagonalize the candidate Hamiltonian and evaluate the loss function.
!> \param reference reference eigenvalues and eigenvectors \f$C_R\f$
!> \param matrix_s_candidate candidate AO overlap matrix \f$S_C\f$
!> \param matrix_s_candidate_reference candidate/reference AO overlap matrix \f$S_{CR}\f$
!> \param virtual_cutoff virtual-orbital energy cutoff
!> \param virtual_smoothing virtual-orbital energy smoothing width
!> \param gap_scale normalization energy for the HOMO-LUMO gap loss
!> \param result loss-function terms and diagnostic quantities
!>
!> The reference Hamiltonian is transformed to the candidate AO basis using
!> \f[
!> B=C_R^\mathrm{T}S_{RC}, \qquad H_C=B^\mathrm{T}
!>       \operatorname{diag}(\epsilon_R)B .
!> \f]
!> The candidate orbitals are obtained from
!> \f[
!> H_C C_C=S_C C_C\operatorname{diag}(\epsilon_C),
!> \qquad C_C^\mathrm{T}S_C C_C=I ,
!> \f]
!> and the reference/candidate orbital overlap matrix is
!> \f[
!> M=C_R^\mathrm{T}S_{RC}C_C=BC_C .
!> \f]
! **************************************************************************************************
   SUBROUTINE evaluate_frontier_orbitals_objective(reference, matrix_s_candidate, &
                                                   matrix_s_candidate_reference, virtual_cutoff, &
                                                   virtual_smoothing, gap_scale, RESULT)
      TYPE(frontier_orbitals_reference_type), INTENT(IN) :: reference
      TYPE(dbcsr_type), POINTER                          :: matrix_s_candidate, &
                                                            matrix_s_candidate_reference
      REAL(KIND=dp), INTENT(IN)                          :: virtual_cutoff, virtual_smoothing, &
                                                            gap_scale
      TYPE(frontier_orbitals_result_type), INTENT(OUT)   :: result

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
      result%condition_number = condnum(2)
      IF (.NOT. ieee_is_finite(result%condition_number) .OR. &
          result%condition_number > 1.0E14_dp .OR. result%condition_number <= 0.0_dp) THEN
         result%loss_occupied_subspace = 1.0E6_dp
         result%loss_virtual_subspace = 1.0E6_dp
         result%loss_empty_subspace = 1.0E6_dp
         result%loss_gap = 1.0E6_dp
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

      ! H_C = B^T diag(epsilon_R) B.
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

      ! H_C C_C = S_C C_C diag(epsilon_C).
      CALL cp_fm_geeig(candidate_hamiltonian, candidate_overlap, candidate_eigenvectors, &
                       eigenvalues_candidate, candidate_work)

      ! M = B C_C, where M_ij = <psi_i^R | psi_j^C>.
      CALL parallel_gemm('N', 'N', nref, naux, naux, 1.0_dp, &
                         reference_orbital_candidate_ao_overlap, candidate_eigenvectors, &
                         0.0_dp, reference_candidate_orbital_overlap)
      ALLOCATE (orbital_overlap_dense(nref, naux))
      CALL cp_fm_get_submatrix(reference_candidate_orbital_overlap, orbital_overlap_dense)
      CALL compute_frontier_orbital_loss_function( &
         orbital_overlap_dense, reference%eigenvalues, eigenvalues_candidate, &
         reference%number_occupied_orbitals, virtual_cutoff, virtual_smoothing, gap_scale, RESULT)
      result%condition_number = condnum(2)

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
!> \param orbital_overlap \f$M_{ij}=\langle\psi_i^R|\psi_j^C\rangle\f$
!> \param eigenvalues_reference reference eigenvalues \f$\epsilon_i^R\f$
!> \param eigenvalues_candidate candidate eigenvalues \f$\epsilon_j^C\f$
!> \param number_occupied_orbitals number of occupied spatial orbitals
!> \param virtual_cutoff virtual-orbital energy cutoff \f$E_\mathrm{cut}\f$
!> \param virtual_smoothing virtual-orbital energy smoothing width \f$\Delta E\f$
!> \param gap_scale normalization energy \f$E_\mathrm{scale}\f$ for the gap loss
!> \param result loss-function terms and diagnostic quantities
!>
!> Let \f$i\f$ denote a reference orbital, \f$j\f$ a candidate orbital, and
!> \f$N_\mathrm{occ}\f$ the number of occupied orbitals. The occupied-subspace loss is
!> \f[
!> \mathcal{L}_\mathrm{occ}
!> =1-\frac{1}{N_\mathrm{occ}}
!>   \sum_{i=1}^{N_\mathrm{occ}}\sum_{j=1}^{N_\mathrm{occ}}|M_{ij}|^2 .
!> \f]
!> For virtual orbitals \f$i,j>N_\mathrm{occ}\f$, define
!> \f[
!> w_i^X=\left[1+\exp\left(
!> \frac{\epsilon_i^X-\epsilon_\mathrm{LUMO}^X-E_\mathrm{cut}}{\Delta E}
!> \right)\right]^{-1},\qquad
!> N_X=\sum_{i>N_\mathrm{occ}}(w_i^X)^2,\quad X\in\{R,C\}.
!> \f]
!> The virtual-subspace loss compares the energy-weighted reference and candidate subspaces:
!> \f[
!> \mathcal{L}_\mathrm{vir}
!> =\frac{N_R+N_C
!> -2\sum_{i>N_\mathrm{occ}}\sum_{j>N_\mathrm{occ}}
!> w_i^R w_j^C |M_{ij}|^2}{2N_R}.
!> \f]
!> With \f$E_g^X=\epsilon_{N_\mathrm{occ}+1}^X-\epsilon_{N_\mathrm{occ}}^X\f$,
!> the gap loss is
!> \f[
!> \mathcal{L}_g=\left(\frac{E_g^C-E_g^R}{E_\mathrm{scale}}\right)^2 .
!> \f]
! **************************************************************************************************
   PURE SUBROUTINE compute_frontier_orbital_loss_function( &
      orbital_overlap, eigenvalues_reference, eigenvalues_candidate, number_occupied_orbitals, &
      virtual_cutoff, virtual_smoothing, gap_scale, RESULT)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: orbital_overlap
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues_reference, &
                                                            eigenvalues_candidate
      INTEGER, INTENT(IN)                                :: number_occupied_orbitals
      REAL(KIND=dp), INTENT(IN)                          :: virtual_cutoff, virtual_smoothing, &
                                                            gap_scale
      TYPE(frontier_orbitals_result_type), INTENT(OUT)   :: result

      INTEGER                                            :: candidate_orbital, reference_orbital
      REAL(KIND=dp) :: candidate_virtual_weight_norm, occupied_subspace_overlap, &
         reference_virtual_weight_norm, reference_virtual_subspace_overlap, &
         virtual_subspace_overlap, weight_candidate, weight_reference

      ! O_occ = sum_{i,j in occupied} |M_ij|^2 / N_occ.
      occupied_subspace_overlap = &
         SUM(orbital_overlap(1:number_occupied_orbitals, &
                             1:number_occupied_orbitals)**2)/ &
         REAL(number_occupied_orbitals, KIND=dp)
      result%loss_occupied_subspace = MAX(0.0_dp, 1.0_dp - occupied_subspace_overlap)
      result%occupied_subspace_overlap = occupied_subspace_overlap

      ! N_R = sum_{i in virtual} (w_i^R)^2.
      reference_virtual_weight_norm = 0.0_dp
      DO reference_orbital = number_occupied_orbitals + 1, SIZE(eigenvalues_reference)
         weight_reference = smooth_frontier_orbital_weight( &
                            eigenvalues_reference(reference_orbital) - &
                            eigenvalues_reference(number_occupied_orbitals + 1), &
                            virtual_cutoff, virtual_smoothing)
         reference_virtual_weight_norm = reference_virtual_weight_norm + weight_reference**2
      END DO

      ! N_C = sum_{j in virtual} (w_j^C)^2.
      candidate_virtual_weight_norm = 0.0_dp
      DO candidate_orbital = number_occupied_orbitals + 1, SIZE(eigenvalues_candidate)
         weight_candidate = smooth_frontier_orbital_weight( &
                            eigenvalues_candidate(candidate_orbital) - &
                            eigenvalues_candidate(number_occupied_orbitals + 1), &
                            virtual_cutoff, virtual_smoothing)
         candidate_virtual_weight_norm = candidate_virtual_weight_norm + weight_candidate**2
      END DO

      ! O_vir = sum_{i,j in virtual} w_i^R w_j^C |M_ij|^2.
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
         result%loss_virtual_subspace = &
            MAX(0.0_dp, (reference_virtual_weight_norm + candidate_virtual_weight_norm - &
                         2.0_dp*virtual_subspace_overlap)/(2.0_dp*reference_virtual_weight_norm))
         result%empty_subspace_overlap = &
            MIN(1.0_dp, MAX(0.0_dp, reference_virtual_subspace_overlap/reference_virtual_weight_norm))
         result%loss_empty_subspace = 1.0_dp - result%empty_subspace_overlap
      END IF

      ! E_g^X = epsilon_LUMO^X - epsilon_HOMO^X.
      result%gap_reference = eigenvalues_reference(number_occupied_orbitals + 1) - &
                             eigenvalues_reference(number_occupied_orbitals)
      result%gap_candidate = eigenvalues_candidate(number_occupied_orbitals + 1) - &
                             eigenvalues_candidate(number_occupied_orbitals)

      ! L_g = [(E_g^C - E_g^R)/E_scale]^2.
      result%loss_gap = ((result%gap_candidate - result%gap_reference)/gap_scale)**2

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
