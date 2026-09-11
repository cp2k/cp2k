!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Automatic RI basis set optimization for molecular GW.
!> \par History
!>      09.2026 created [Jan Wilhelm]
! **************************************************************************************************
MODULE gw_auto_ri_generate_RI_basis
   USE ai_contraction_sphi,             ONLY: ab_contract,&
                                              abc_contract
   USE atomic_kind_types,               ONLY: get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE constants_operator,              ONLY: operator_coulomb,&
                                              operator_truncated
   USE coulomb_integral_interface,      ONLY: compute_coulomb_2c,&
                                              compute_coulomb_3c,&
                                              coulomb_integral_cleanup,&
                                              coulomb_integral_context_type,&
                                              coulomb_integral_init,&
                                              coulomb_operator_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_put_block, dbcsr_release, dbcsr_type, &
        dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_files,                        ONLY: close_file,&
                                              discover_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_syrk,&
                                              cp_fm_uplo_to_full
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE generic_os_integrals,            ONLY: int_operators_r12_ab_os
   USE gw_auto_ri_types,                ONLY: auto_ri_type
   USE gw_integrals,                    ONLY: build_3c_integral_block_ctx,&
                                              gw_3c_ctx_create,&
                                              gw_3c_ctx_release,&
                                              gw_3c_ctx_type,&
                                              gw_3c_ws_create,&
                                              gw_3c_ws_release,&
                                              gw_3c_ws_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_truncated,&
                                              small_cell_full_kp
   USE integral_library_types,          ONLY: active_integral_library,&
                                              library_libint
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_HOST,&
                                              local_gemm_ctxt_type
   USE machine,                         ONLY: m_walltime
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              init_t_c_g0 => init
   USE util,                            ONLY: sort
#include "base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_auto_ri_generate_RI_basis'

   PUBLIC :: generate_auto_ri_basis

   TYPE fit_group_list_type
      INTEGER, ALLOCATABLE                               :: group(:)
      INTEGER, ALLOCATABLE                               :: first_group_column(:)
      INTEGER, ALLOCATABLE                               :: first_block_column(:)
      INTEGER, ALLOCATABLE                               :: ncolumn(:)
   END TYPE fit_group_list_type

   TYPE optimized_atom_transform_type
      INTEGER, ALLOCATABLE                               :: reference_atom(:)
      INTEGER, ALLOCATABLE                               :: first_row(:)
      REAL(KIND=dp), ALLOCATABLE                         :: coefficient(:, :)
   END TYPE optimized_atom_transform_type

CONTAINS

!**************************************************************************************************!
! The full derivation of the AUTO_RI algorithm will be published separately.                       !
! The equations below define the notation used throughout this module.                             !
!                                                                                                  !
! The RI scalar product uses the truncated Coulomb operator:                                       !
!                                                                                                  !
! (f|g) = ∫dr∫dr' f(r) V_rc(r-r') g(r').                                                   (1)     !
!                                                                                                  !
! For the RI basis-set optimization, the orbital-product density formed by                         !
! AOs on atoms A and B is                                                                          !
!                                                                                                  !
! ρ_mn^AB(r) = Σ_{μ∈A} Σ_{ν∈B} C_μm C_νn Φ_μ(r) Φ_ν(r).                                    (2)     !
!                                                                                                  !
! The atom-local correlation target mimics the density response χ, whereas                         !
! the exchange target describes the exchange contribution to the self-energy:                      !
!                                                                                                  !
! B_PQ^{A,C} = (1/N_C^A) Σ_ia^oe w_ia (φ_P^A|ρ_ia^AA)(ρ_ia^AA|φ_Q^A).                      (3)     !
!                                                                                                  !
! B_PQ^{A,X} = (1/N_X^A) Σ_{n∈W_Δ} Σ_i^occ (φ_P^A|ρ_ni^AA)(ρ_ni^AA|φ_Q^A).                 (4)     !
!                                                                                                  !
! Here w_ia gives larger weight to frontier-orbital products close to the HOMO                     !
! and LUMO:                                                                                        !
!                                                                                                  !
! w_ia = 1/(|ε_i-E_F| |ε_a-E_F|),   E_F = (ε_HOMO+ε_LUMO)/2.                               (5)     !
!                                                                                                  !
! Diagonalizing B^A = B^{A,C}+B^{A,X} identifies the most important                                !
! atom-local RI functions. The contraction coefficients solve                                      !
!                                                                                                  !
! B^A U^A = M^A U^A λ^A,   Σ_Q B_PQ^A U_Qp^A = λ_p^A Σ_Q M_PQ^A U_Qp^A.                    (6)     !
!                                                                                                  !
! The eigenvectors with the largest eigenvalues are the contracted RI functions used for the       !
! optimized RI basis set:                                                                          !
!                                                                                                  !
! φ_p^A(r) = Σ_P U_Pp^A φ_P^A(r).                                                          (7)     !
!                                                                                                  !
! We also use superpositions of Gaussians on neighboring atoms in the RI basis. Before optimizing  !
! the expansion coefficients, the AA from Eq. (7) used in AB contractions are projected out:       !
!                                                                                                  !
! span{φ_P^{AB,⊥}} = span{φ_P^A,φ_Q^B} \ span{φ_p^A,φ_q^B}.                                (8)     !
!                                                                                                  !
! For obtaining the AB contraction coefficients, we use the symmetric off-diagonal density:        !
!                                                                                                  !
! ρ̅_mn^AB(r) = ρ_mn^AB(r)+ρ_mn^BA(r).                                                      (9)     !
!                                                                                                  !
! To identify the most important functions for a neighboring pair A-B, define                      !
! correlation and exchange targets in this projected basis:                                        !
!                                                                                                  !
! B_PQ^{AB,C} = (1/N_C^{AB}) Σ_ia^oe w_ia (φ_P^{AB,⊥}|ρ̅_ia^AB)(ρ̅_ia^AB|φ_Q^{AB,⊥}).       (10)     !
!                                                                                                  !
! B_PQ^{AB,X} = (1/N_X^{AB}) Σ_{n∈W_Δ} Σ_i^occ (φ_P^{AB,⊥}|ρ̅_ni^AB)(ρ̅_ni^AB|φ_Q^{AB,⊥}).  (11)     !
!                                                                                                  !
! Diagonalizing their sum gives the AB-neighbor-pair contraction coefficients:                     !
!                                                                                                  !
! B^AB U^AB = U^AB λ^AB,   B^AB = B^{AB,C}+B^{AB,X}.                                      (12)     !
!                                                                                                  !
! The eigenvectors with the largest eigenvalues expand the pair functions:                         !
!                                                                                                  !
! φ_p^AB(r) = Σ_P U_Pp^AB φ_P^{AB,⊥}(r).                                                  (13)     !
!                                                                                                  !
! AUTO_RI uses the same number of AA and AB contractions in the final optimized RI basis:          !
!                                                                                                  !
! Σ_A N_A = Σ_{A>B} N_AB = N_RI/2.                                                        (14)     !
!                                                                                                  !
!**************************************************************************************************!

! **************************************************************************************************
!> \brief Executes the AUTO_RI algorithm defined by Eqs. (1)-(14).
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE generate_auto_ri_basis(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'generate_auto_ri_basis'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: time_start

      CALL timeset(routineN, handle)

      time_start = m_walltime()

      ! Compute the AA target matrices and contractions in Eqs. (3)-(7).
      CALL compute_AA_contractions(bs_env)

      ! Compute the AB residual target matrices and contractions in Eqs. (8)-(13).
      CALL compute_AB_contractions(qs_env, bs_env)

      ! Transform and invert the optimized RI metric matrix M_PQ = (P|Q).
      CALL compute_optimized_RI_metric_matrix(qs_env, bs_env)

      ! Activate the optimized per-atom RI block sizes for RI-RS.
      CALL set_optimized_RI_basis_sizes(bs_env)

      ! Print the optimized RI basis dimension and total optimization time.
      CALL print_auto_ri_summary(bs_env, time_start)

      CALL timestop(handle)

   END SUBROUTINE generate_auto_ri_basis

! **************************************************************************************************
!> \brief Computes the atom-local contraction coefficients U^A from Eqs. (3)-(7):
!>
!>          B_PQ^{A,C} = (1/N_C^A) Σ_ia^oe w_ia (φ_P^A|ρ_ia^AA)(ρ_ia^AA|φ_Q^A),         (3)
!>
!>          B_PQ^{A,X} = (1/N_X^A) Σ_{n∈W_Δ} Σ_i^occ (φ_P^A|ρ_ni^AA)(ρ_ni^AA|φ_Q^A),    (4)
!>
!>          w_ia = 1/(|ε_i-E_F| |ε_a-E_F|),  E_F = (ε_HOMO+ε_LUMO)/2,                   (5)
!>
!>          B^A U^A = M^A U^A λ^A,                                                      (6)
!>
!>          Contracted RI basis functions for optimized RI basis:
!>
!>          φ_p^A(r) = Σ_P U_Pp^A φ_P^A(r).                                             (7)
!>
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_AA_contractions(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_AA_contractions'

      INTEGER                                            :: handle, iatom, ikind, max_ao, &
                                                            max_reference, n_AA_requested, &
                                                            nao_kind, natom, nkind, root_rank
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: counts, kind_of
      REAL(KIND=dp)                                      :: normalization
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: B_PQ, coefficients, metric, scores
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: fitting_metric_by_kind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: atom_local_RI_integrals
      TYPE(libint_potential_type)                        :: fitting_potential
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      ! Determine reference dimensions and the number of AA RI functions (i.e. number of
      ! eigenvectors U_Pp^A used for RI-basis contractions).
      CALL initialize_AA_optimization(bs_env, para_env, kind_of, natom, nkind, &
                                      max_reference, n_AA_requested)
      ALLOCATE (scores(max_reference, natom), source=0.0_dp)
      ALLOCATE (counts(natom))

      CALL initialize_auto_ri_integrals(bs_env)
      fitting_potential = bs_env%ri_metric
      max_ao = MAXVAL(bs_env%sizes_AO)
      ALLOCATE (bs_env%auto_ri%transform(max_reference, max_reference, natom), &
                fitting_metric_by_kind(max_reference, max_reference, nkind), source=0.0_dp)
      ALLOCATE (atom_local_RI_integrals(max_ao, max_ao, max_reference, nkind), source=0.0_dp)

      ! The few atom-kind integral blocks are inexpensive; compute them on every rank.
      DO ikind = 1, nkind
         nao_kind = bs_env%basis_set_AO(ikind)%gto_basis_set%nsgf
         CALL compute_atom_local_RI_integrals( &
            ao_basis=bs_env%basis_set_AO(ikind)%gto_basis_set, &
            ri_basis=bs_env%basis_set_RI(ikind)%gto_basis_set, &
            potential=fitting_potential, metric=metric, &
            products=atom_local_RI_integrals( &
            1:nao_kind, 1:nao_kind, 1:bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf, ikind))
         CALL add_RI_metric_regularization(bs_env%basis_set_RI(ikind)%gto_basis_set, &
                                           bs_env%regularization_RI, metric)
         fitting_metric_by_kind(1:SIZE(metric, 1), 1:SIZE(metric, 2), ikind) = metric
         DEALLOCATE (metric)
      END DO

      bs_env%auto_ri%optimized_size = 0
      DO iatom = 1, natom
         root_rank = MODULO(iatom - 1, para_env%num_pe)
         ikind = kind_of(iatom)
         nao_kind = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
         ALLOCATE (metric(nao_kind, nao_kind))
         metric(:, :) = fitting_metric_by_kind(1:nao_kind, 1:nao_kind, ikind)
         ! Eqs. (3)-(5): form B^{A,C}, B^{A,X}, and their weighted sum B^A.
         CALL compute_B_PQ_AA( &
            bs_env, iatom, root_rank, &
            atom_local_RI_integrals(1:bs_env%sizes_AO(iatom), &
                                    1:bs_env%sizes_AO(iatom), 1:nao_kind, ikind), &
            metric, B_PQ, normalization)
         DEALLOCATE (metric)
         IF (root_rank /= para_env%mepos) THEN
            DEALLOCATE (B_PQ)
            CYCLE
         END IF
         ! The fitting metric, target and product norm must use the same operator.
         nao_kind = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
         ! Eq. (6): solve B^A U^A = M^A U^A λ^A.
         CALL solve_AA_eigenproblem(fitting_metric_by_kind(1:nao_kind, 1:nao_kind, ikind), &
                                    B_PQ, normalization, bs_env%auto_ri%min_score_ratio, &
                                    bs_env%auto_ri%metric_null_threshold, &
                                    bs_env%auto_ri%psd_tolerance, coefficients, eigenvalues)
         scores(1:SIZE(coefficients, 2), iatom) = &
            eigenvalues(1:SIZE(coefficients, 2))
         bs_env%auto_ri%optimized_size(iatom) = SIZE(coefficients, 2)
         bs_env%auto_ri%transform(1:SIZE(coefficients, 1), &
                                  1:SIZE(coefficients, 2), iatom) = coefficients
         DEALLOCATE (coefficients, eigenvalues, B_PQ)
      END DO
      CALL para_env%sum(bs_env%auto_ri%optimized_size)
      CALL para_env%sum(bs_env%auto_ri%transform)
      CALL para_env%sum(scores)
      CALL select_RI_counts(scores, n_AA_requested, 1, counts)
      bs_env%auto_ri%optimized_size(:) = counts
      DO iatom = 1, natom
         bs_env%auto_ri%transform(:, counts(iatom) + 1:, iatom) = 0.0_dp
      END DO
      DEALLOCATE (atom_local_RI_integrals, fitting_metric_by_kind)

      CALL timestop(handle)

   END SUBROUTINE compute_AA_contractions

! **************************************************************************************************
!> \brief Initializes dimensions and atom-kind mappings for the AA optimization.
!>
!>        This routine does not evaluate an equation of the objective. It maps every atom
!>        to its reference RI basis, allocates the persistent per-atom arrays, and determines
!>        the AA part of the requested basis size from
!>
!>          Σ_A N_A = N_RI/2.                                             (14)
!>
!>        At least one AA contraction is assigned to every atom.
!> \param bs_env ...
!> \param para_env ...
!> \param kind_of ...
!> \param natom ...
!> \param nkind ...
!> \param max_reference ...
!> \param n_AA_requested ...
! **************************************************************************************************
   SUBROUTINE initialize_AA_optimization(bs_env, para_env, kind_of, &
                                         natom, nkind, max_reference, &
                                         n_AA_requested)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: kind_of(:)
      INTEGER, INTENT(OUT)                               :: natom, nkind, max_reference, &
                                                            n_AA_requested

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_AA_optimization'

      INTEGER                                            :: handle, ikind, n_RI_requested
      INTEGER, ALLOCATABLE                               :: reference_size(:)

      CALL timeset(routineN, handle)

      para_env => bs_env%para_env
      natom = bs_env%n_atom
      nkind = SIZE(bs_env%basis_set_RI)
      CPASSERT(ASSOCIATED(bs_env%ri_rs%atomic_kind_set))
      CPASSERT(ASSOCIATED(bs_env%ri_rs%particle_set))
      CPASSERT(ASSOCIATED(para_env) .AND. natom > 0 .AND. nkind > 0)
      CPASSERT(ALLOCATED(bs_env%basis_set_AO))
      CPASSERT(ALLOCATED(bs_env%basis_set_RI))
      CPASSERT(SIZE(bs_env%basis_set_AO) == nkind)
      CPASSERT(SIZE(bs_env%basis_set_RI) == nkind)
      CPASSERT(bs_env%do_gw_ri_rs)
      CPASSERT(.NOT. bs_env%do_periodic)
      CPASSERT(bs_env%small_cell_full_kp_or_large_cell_Gamma /= small_cell_full_kp)

      CALL get_atomic_kind_set(atomic_kind_set=bs_env%ri_rs%atomic_kind_set, &
                               kind_of=kind_of)
      ALLOCATE (bs_env%auto_ri%optimized_size(natom))
      ALLOCATE (reference_size(nkind))
      DO ikind = 1, nkind
         CPASSERT(ASSOCIATED(bs_env%basis_set_AO(ikind)%gto_basis_set))
         CPASSERT(ASSOCIATED(bs_env%basis_set_RI(ikind)%gto_basis_set))
         reference_size(ikind) = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
      END DO
      max_reference = MAXVAL(reference_size)
      n_RI_requested = NINT( &
                       MIN(bs_env%auto_ri%ri_ao_ratio, &
                           REAL(SUM(reference_size(kind_of)), dp)/ &
                           REAL(SUM(bs_env%sizes_AO), dp))*REAL(SUM(bs_env%sizes_AO), dp))
      n_AA_requested = MAX(natom, NINT(0.5_dp*REAL(n_RI_requested, dp)))
      DEALLOCATE (reference_size)

      CALL timestop(handle)

   END SUBROUTINE initialize_AA_optimization

! **************************************************************************************************
!> \brief Sets the RI-RS block sizes to the optimized RI basis dimensions.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE set_optimized_RI_basis_sizes(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_optimized_RI_basis_sizes'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      bs_env%sizes_RI(:) = bs_env%auto_ri%optimized_size
      bs_env%auto_ri%ready = .TRUE.

      CALL timestop(handle)

   END SUBROUTINE set_optimized_RI_basis_sizes

! **************************************************************************************************
!> \brief Prints the optimized RI basis size and total optimization time.
!> \param bs_env ...
!> \param time_start ...
! **************************************************************************************************
   SUBROUTINE print_auto_ri_summary(bs_env, time_start)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: time_start

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_auto_ri_summary'

      INTEGER                                            :: handle, n_optimized, unit_nr
      REAL(KIND=dp)                                      :: elapsed

      CALL timeset(routineN, handle)

      n_optimized = SUM(bs_env%auto_ri%optimized_size)
      elapsed = m_walltime() - time_start
      CALL bs_env%para_env%max(elapsed)
      unit_nr = bs_env%unit_nr
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(A)') ' '
         WRITE (unit_nr, '(T2,A,T74,I9)') &
            'AUTO_RI| Number of automatic RI functions for χ, ε, W:', n_optimized
         WRITE (unit_nr, '(T2,A,T71,F8.1,A)') &
            'AUTO_RI| Total RI basis optimization time:', elapsed, ' s'
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_auto_ri_summary

! **************************************************************************************************
!> \brief Adds REGULARIZATION_RI to a two-center RI metric matrix (P|Q).
!> \param basis ...
!> \param regularization ...
!> \param metric ...
!> \param offset ...
! **************************************************************************************************
   SUBROUTINE add_RI_metric_regularization(basis, regularization, metric, offset)
      USE basis_set_types, ONLY: gto_basis_set_type
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(KIND=dp), INTENT(IN)                          :: regularization
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: metric
      INTEGER, INTENT(IN), OPTIONAL                      :: offset

      CHARACTER(LEN=*), PARAMETER :: routineN = 'add_RI_metric_regularization'

      INTEGER                                            :: first, handle, i, iset, nset, start
      INTEGER, DIMENSION(:), POINTER                     :: npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(KIND=dp)                                      :: shift
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet

      CALL timeset(routineN, handle)

      start = 0
      IF (PRESENT(offset)) start = offset
      CALL get_gto_basis_set(basis, first_sgf=first_sgf, npgf=npgf, nset=nset, &
                             nsgf_set=nsgf_set, zet=zet)
      CPASSERT(start >= 0 .AND. start + basis%nsgf <= SIZE(metric, 1))
      CPASSERT(start + basis%nsgf <= SIZE(metric, 2))
      DO iset = 1, nset
         shift = regularization*MAX(1.0_dp, 1.0_dp/MINVAL(zet(1:npgf(iset), iset)))
         first = start + first_sgf(1, iset)
         DO i = first, first + nsgf_set(iset) - 1
            metric(i, i) = metric(i, i) + shift
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE add_RI_metric_regularization

! **************************************************************************************************
!> \brief Copies atom-local diagonal blocks M^AA from a distributed metric M.
!> \param local_metric locally stored entries of M
!> \param AA_metric locally stored entries of the block-diagonal matrix
!> \param row_indices global indices of the locally owned rows
!> \param col_indices global indices of the locally owned columns
!> \param block_sizes optimized RI functions in each atomic matrix block
! **************************************************************************************************
   SUBROUTINE copy_AA_metric_blocks(local_metric, AA_metric, row_indices, &
                                    col_indices, block_sizes)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: local_metric
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: AA_metric
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_indices, col_indices, block_sizes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_AA_metric_blocks'

      INTEGER                                            :: handle, i, j, last_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: block_ends, col_atom, row_atom

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(local_metric, 1) == SIZE(row_indices))
      CPASSERT(SIZE(local_metric, 2) == SIZE(col_indices))
      CPASSERT(ALL(SHAPE(AA_metric) == SHAPE(local_metric)))
      CPASSERT(SIZE(block_sizes) > 0 .AND. ALL(block_sizes > 0))

      ALLOCATE (block_ends(SIZE(block_sizes)))
      block_ends(1) = block_sizes(1)
      DO i = 2, SIZE(block_sizes)
         block_ends(i) = block_ends(i - 1) + block_sizes(i)
      END DO
      last_index = block_ends(SIZE(block_ends))
      CPASSERT(ALL(row_indices >= 1 .AND. row_indices <= last_index))
      CPASSERT(ALL(col_indices >= 1 .AND. col_indices <= last_index))

      ALLOCATE (row_atom(SIZE(row_indices)), col_atom(SIZE(col_indices)))
      DO i = 1, SIZE(row_indices)
         row_atom(i) = metric_index_atom(row_indices(i), block_ends)
      END DO
      DO j = 1, SIZE(col_indices)
         col_atom(j) = metric_index_atom(col_indices(j), block_ends)
      END DO

      AA_metric = 0.0_dp
      DO j = 1, SIZE(col_indices)
         DO i = 1, SIZE(row_indices)
            IF (row_atom(i) == col_atom(j)) THEN
               AA_metric(i, j) = local_metric(i, j)
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE copy_AA_metric_blocks

! **************************************************************************************************
!> \brief Forms f_i=1/sqrt(λ_i) for eigenvalues of the RI metric matrix M_PQ.
!> \param eigenvalues generalized metric eigenvalues
!> \param relative_threshold smallest retained fraction of the largest eigenvalue
!> \param factors inverse-square-root weights, or zero for discarded modes
! **************************************************************************************************
   SUBROUTINE inverse_sqrt_metric_eigenvalues(eigenvalues, relative_threshold, factors)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenvalues
      REAL(KIND=dp), INTENT(IN)                          :: relative_threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: factors

      CHARACTER(LEN=*), PARAMETER :: routineN = 'inverse_sqrt_metric_eigenvalues'

      INTEGER                                            :: handle, i, n_dependent
      REAL(KIND=dp)                                      :: largest, threshold

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(eigenvalues) > 0)
      CPASSERT(relative_threshold >= 0.0_dp .AND. relative_threshold <= 1.0_dp)

      largest = MAXVAL(eigenvalues)
      CPASSERT(largest > 0.0_dp)
      threshold = relative_threshold*largest
      ALLOCATE (factors(SIZE(eigenvalues)))
      factors = 0.0_dp
      n_dependent = 0
      DO i = 1, SIZE(eigenvalues)
         IF (eigenvalues(i) > 0.0_dp .AND. eigenvalues(i) >= threshold) THEN
            factors(i) = 1.0_dp/SQRT(eigenvalues(i))
         ELSE
            n_dependent = n_dependent + 1
         END IF
      END DO
      CPASSERT(n_dependent < SIZE(eigenvalues))

      CALL timestop(handle)

   END SUBROUTINE inverse_sqrt_metric_eigenvalues

! **************************************************************************************************
!> \brief Finds the atom owning a global RI metric index.
!> \param global_index global metric index
!> \param block_ends inclusive end index of each consecutive atom block
!> \return atom index
! **************************************************************************************************
   PURE FUNCTION metric_index_atom(global_index, block_ends) RESULT(atom_index)
      INTEGER, INTENT(IN)                                :: global_index
      INTEGER, DIMENSION(:), INTENT(IN)                  :: block_ends
      INTEGER                                            :: atom_index

      INTEGER                                            :: high, low, middle

      low = 1
      high = SIZE(block_ends)
      DO WHILE (low < high)
         middle = (low + high)/2
         IF (global_index <= block_ends(middle)) THEN
            high = middle
         ELSE
            low = middle + 1
         END IF
      END DO
      atom_index = low

   END FUNCTION metric_index_atom

! **************************************************************************************************
!> \brief Computes X=U diag(1/sqrt(λ)) with X^T M X=I.
!>
!>        The columns of X span the non-null eigenspace of the symmetric
!>        two-center RI metric matrix M=U diag(λ) U^T.
!> \param matrix ...
!> \param factor ...
!> \param rank ...
!> \param null_threshold Relative eigenvalue threshold for the null space.
!> \param psd_threshold Relative tolerance for negative eigenvalues.
! **************************************************************************************************
   SUBROUTINE compute_metric_inverse_factor(matrix, factor, rank, null_threshold, psd_threshold)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: factor
      INTEGER, INTENT(OUT)                               :: rank
      REAL(KIND=dp), INTENT(IN)                          :: null_threshold, psd_threshold

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_metric_inverse_factor'

      INTEGER                                            :: column, handle, i, n
      REAL(KIND=dp)                                      :: largest
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvectors

      CALL timeset(routineN, handle)

      n = SIZE(matrix, 1)
      CPASSERT(n > 0 .AND. SIZE(matrix, 2) == n)
      ALLOCATE (eigenvectors, source=0.5_dp*(matrix + TRANSPOSE(matrix)))
      ALLOCATE (eigenvalues(n))
      CALL diamat_all(eigenvectors, eigenvalues, dac=.TRUE.)
      largest = MAXVAL(eigenvalues)
      IF (largest <= 0.0_dp) CPABORT("AUTO_RI metric has no positive non-null range")
      IF (MINVAL(eigenvalues) < -psd_threshold*largest) THEN
         CPABORT("AUTO_RI metric has no positive non-null range")
      END IF
      rank = COUNT(eigenvalues > null_threshold*largest)
      IF (rank == 0) CPABORT("AUTO_RI metric has no positive non-null range")
      ALLOCATE (factor(n, rank))
      column = 0
      DO i = 1, n
         IF (eigenvalues(i) <= null_threshold*largest) CYCLE
         column = column + 1
         factor(:, column) = eigenvectors(:, i)/SQRT(eigenvalues(i))
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_metric_inverse_factor

! **************************************************************************************************
!> \brief Solves eigenvalue problem B^A U^A=M^A U^A λ^A in Eq. (6).
!>
!>        Columns are ordered from the largest to the smallest λ^A.
!>        Columns are retained when λ_i/λ_max >= min_eigenvalue_kept.
!> \param metric ...
!> \param B_PQ normalized atom-local target matrix B^A
!> \param target_normalization ...
!> \param min_eigenvalue_kept ...
!> \param metric_null_threshold Relative eigenvalue threshold for the metric
!>        null space.
!> \param psd_tolerance Relative tolerance for negative metric eigenvalues.
!> \param coefficients ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE solve_AA_eigenproblem(metric, B_PQ, target_normalization, &
                                    min_eigenvalue_kept, metric_null_threshold, &
                                    psd_tolerance, coefficients, eigenvalues)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: metric, B_PQ
      REAL(KIND=dp), INTENT(IN)                          :: target_normalization, &
                                                            min_eigenvalue_kept, &
                                                            metric_null_threshold, psd_tolerance
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: coefficients
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: eigenvalues

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_AA_eigenproblem'

      INTEGER                                            :: handle, i, j, metric_rank, n, nret
      REAL(KIND=dp)                                      :: largest_score, relative_fitting_defect
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: overlap_values, score_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: check, normalized_target, rotated, &
                                                            score_vectors, whitener

      CALL timeset(routineN, handle)

      n = SIZE(metric, 1)
      CPASSERT(n > 0 .AND. SIZE(metric, 2) == n .AND. ALL(SHAPE(B_PQ) == [n, n]))
      CPASSERT(target_normalization > 0.0_dp)
      CPASSERT(min_eigenvalue_kept > 0.0_dp .AND. min_eigenvalue_kept <= 1.0_dp)
      CALL compute_metric_inverse_factor(metric, whitener, metric_rank, &
                                         metric_null_threshold, psd_tolerance)

      ALLOCATE (normalized_target(metric_rank, metric_rank))
      normalized_target(:, :) = MATMUL(TRANSPOSE(whitener), &
                                       MATMUL(0.5_dp*(B_PQ + TRANSPOSE(B_PQ)), &
                                              whitener))/target_normalization
      ALLOCATE (score_vectors, source=normalized_target)
      ALLOCATE (score_values(metric_rank))
      CALL diamat_all(score_vectors, score_values, dac=.TRUE.)
      largest_score = MAXVAL(score_values)
      IF (largest_score <= 0.0_dp) THEN
         CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
      END IF
      IF (MINVAL(score_values) < -psd_tolerance*largest_score) THEN
         CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
      END IF
      nret = COUNT(score_values/largest_score >= min_eigenvalue_kept)
      ALLOCATE (coefficients(n, nret), eigenvalues(metric_rank))
      DO i = 1, metric_rank
         eigenvalues(i) = score_values(metric_rank - i + 1)
      END DO
      DO i = 1, nret
         j = metric_rank - i + 1
         coefficients(:, i) = MATMUL(whitener, score_vectors(:, j))
      END DO

      relative_fitting_defect = 1.0_dp
      DO i = 1, nret
         j = metric_rank - i + 1
         relative_fitting_defect = relative_fitting_defect - &
                                   DOT_PRODUCT(score_vectors(:, j), &
                                               MATMUL(normalized_target, score_vectors(:, j)))
      END DO
      IF (ABS(relative_fitting_defect) < 100.0_dp*EPSILON(1.0_dp)) THEN
         relative_fitting_defect = 0.0_dp
      END IF
      ! A projection cannot capture more than the exact product norm.
      IF (relative_fitting_defect < -1.0E-8_dp .OR. &
          relative_fitting_defect > 1.0_dp + 1.0E-8_dp) THEN
         CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
      END IF
      IF (nret > 0) THEN
         ALLOCATE (check(nret, nret), overlap_values(nret), rotated(n, nret))
         check(:, :) = MATMUL(TRANSPOSE(coefficients), &
                              MATMUL(0.5_dp*(metric + TRANSPOSE(metric)), coefficients))
         CALL diamat_all(check, overlap_values, dac=.TRUE.)
         IF (MINVAL(overlap_values) <= metric_null_threshold) THEN
            CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
         END IF
         rotated(:, :) = MATMUL(coefficients, check)
         DO i = 1, nret
            rotated(:, i) = rotated(:, i)/SQRT(overlap_values(i))
         END DO
         coefficients(:, :) = MATMUL(rotated, TRANSPOSE(check))
         check(:, :) = MATMUL(TRANSPOSE(coefficients), &
                              MATMUL(0.5_dp*(metric + TRANSPOSE(metric)), coefficients))
         DO i = 1, nret
            check(i, i) = check(i, i) - 1.0_dp
         END DO
         IF (MAXVAL(ABS(check)) >= 1.0E-7_dp) THEN
            CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE solve_AA_eigenproblem

! **************************************************************************************************
!> \brief Solves eigenvalue problem B^AB U^AB=U^AB λ^AB of Eq. (8).
!>
!>        Columns are ordered from the largest to the smallest λ^AB.
!> \param metric combined AB reference-basis metric
!> \param B_PQ normalized two-center target matrix B^AB
!> \param atomic retained AA columns embedded in the combined AB reference space
!> \param relative_cutoff Minimum score relative to the largest residual score.
!> \param absolute_cutoff Minimum score on the unprojected normalization scale.
!> \param metric_null_threshold Relative eigenvalue threshold for the metric
!>        null space.
!> \param psd_tolerance Relative tolerance for negative metric eigenvalues.
!> \param coefficients Selected physical pair contractions.
!> \param eigenvalues normalized scores in descending order
! **************************************************************************************************
   SUBROUTINE solve_AB_eigenproblem(metric, B_PQ, atomic, relative_cutoff, &
                                    absolute_cutoff, metric_null_threshold, psd_tolerance, &
                                    coefficients, eigenvalues)
      REAL(KIND=dp), INTENT(IN) :: metric(:, :), B_PQ(:, :), atomic(:, :), relative_cutoff, &
         absolute_cutoff, metric_null_threshold, psd_tolerance
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: coefficients(:, :), eigenvalues(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_AB_eigenproblem'

      INTEGER                                            :: handle, i, k, n, nkeep, rank
      REAL(KIND=dp)                                      :: threshold
      REAL(KIND=dp), ALLOCATABLE                         :: eig(:), gram(:, :), projection(:, :), &
                                                            q(:, :), residual(:, :), score(:, :), &
                                                            u(:, :), vectors(:, :), x(:, :)

      CALL timeset(routineN, handle)

      n = SIZE(metric, 1)
      CALL compute_metric_inverse_factor(metric, x, rank, metric_null_threshold, psd_tolerance)

      ! Xᵀ M X=I; Q spans the retained atomic columns in these orthonormal coordinates.
      ALLOCATE (q(rank, SIZE(atomic, 2)))
      q(:, :) = MATMUL(TRANSPOSE(x), MATMUL(metric, atomic))
      ALLOCATE (projection(rank, rank), source=0.0_dp)
      DO i = 1, rank
         projection(i, i) = 1.0_dp
      END DO
      IF (SIZE(atomic, 2) > 0) THEN
         ALLOCATE (gram(SIZE(atomic, 2), SIZE(atomic, 2)))
         gram(:, :) = MATMUL(TRANSPOSE(q), q)
         CALL compute_metric_inverse_factor(gram, u, k, metric_null_threshold, psd_tolerance)
         ALLOCATE (vectors(rank, SIZE(u, 2)))
         vectors(:, :) = MATMUL(q, u)
         projection(:, :) = projection - MATMUL(vectors, TRANSPOSE(vectors))
         DEALLOCATE (vectors)
      END IF
      ALLOCATE (score(rank, rank))
      score(:, :) = MATMUL(TRANSPOSE(x), MATMUL(B_PQ, x))
      ALLOCATE (residual(rank, rank))
      residual(:, :) = MATMUL(projection, MATMUL(score, projection))
      ALLOCATE (eig(rank))
      CALL diamat_all(residual, eig, dac=.TRUE.)
      threshold = MAX(absolute_cutoff, relative_cutoff*MAXVAL(eig), &
                      100.0_dp*EPSILON(1.0_dp)*MAXVAL(ABS(score)))
      nkeep = COUNT(eig > threshold)
      ALLOCATE (coefficients(n, nkeep), eigenvalues(nkeep))
      k = 0
      DO i = rank, 1, -1
         IF (eig(i) <= threshold) CYCLE
         k = k + 1
         coefficients(:, k) = MATMUL(x, MATMUL(projection, residual(:, i)))
         eigenvalues(k) = eig(i)
      END DO

      CALL timestop(handle)

   END SUBROUTINE solve_AB_eigenproblem

! **************************************************************************************************
!> \brief Selects the requested number of largest positive eigenvalues.
!> \param scores descending eigenvalues; zero pads unavailable modes
!> \param n_requested requested total number of functions
!> \param minimum Minimum retained count in every nonempty block.
!> \param counts number retained from each block
! **************************************************************************************************
   SUBROUTINE select_RI_counts(scores, n_requested, minimum, counts)
      REAL(KIND=dp), INTENT(IN)                          :: scores(:, :)
      INTEGER, INTENT(IN)                                :: n_requested, minimum
      INTEGER, INTENT(OUT)                               :: counts(:)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'select_RI_counts'

      INTEGER                                            :: block, handle, i, n, nmode, remaining
      INTEGER, ALLOCATABLE                               :: indices(:)
      REAL(KIND=dp), ALLOCATABLE                         :: values(:)

      CALL timeset(routineN, handle)

      nmode = SIZE(scores, 1)
      n = SIZE(scores)
      CPASSERT(SIZE(counts) == SIZE(scores, 2) .AND. minimum >= 0)
      counts(:) = 0
      IF (n > 0) THEN
         ALLOCATE (values(n), indices(n))
         values(:) = -RESHAPE(scores, [n])
         DO block = 1, SIZE(counts)
            counts(block) = MIN(minimum, COUNT(scores(:, block) > 0.0_dp))
            values((block - 1)*nmode + 1:(block - 1)*nmode + counts(block)) = 0.0_dp
         END DO
         CALL sort(values, n, indices)
         remaining = MIN(MAX(0, n_requested - SUM(counts)), COUNT(values < 0.0_dp))
         DO i = 1, remaining
            block = (indices(i) - 1)/nmode + 1
            counts(block) = counts(block) + 1
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE select_RI_counts

! **************************************************************************************************
!> \brief Computes orbital factors for the occupied-empty weight in Eq. (5).
!>
!>        weights(n) = 1/|ε_n-E_F|, so weights(i)*weights(a) = w_ia.
!> \param energies molecular eigenvalues, occupied first
!> \param homo Last occupied state in this spin channel.
!> \param weights factors whose occupied-empty products give w_ia
! **************************************************************************************************
   SUBROUTINE compute_orbital_weights(energies, homo, weights)
      REAL(KIND=dp), INTENT(IN)                          :: energies(:)
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(OUT)                         :: weights(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_orbital_weights'

      INTEGER                                            :: handle, n
      REAL(KIND=dp)                                      :: gap, middle

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(weights) == SIZE(energies))
      IF (homo < 1 .OR. homo >= SIZE(energies)) THEN
         CPABORT("AUTO_RI orbital weights require occupied and empty states")
      END IF
      IF (.NOT. ALL(ieee_is_finite(energies))) THEN
         CPABORT("AUTO_RI orbital weights require finite eigenvalues")
      END IF
      gap = energies(homo + 1) - energies(homo)
      IF (.NOT. ieee_is_finite(gap) .OR. gap <= 0.0_dp) THEN
         CPABORT("AUTO_RI orbital weights require a finite positive gap")
      END IF
      middle = 0.5_dp*(energies(homo) + energies(homo + 1))
      DO n = 1, SIZE(energies)
         weights(n) = 1.0_dp/ABS(energies(n) - middle)
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_orbital_weights

! **************************************************************************************************
!> \brief Computes the occupied-empty contribution B_PQ^C in Eqs. (3) or (10).
!>
!>        B_PQ^C = Σ_ia w_ia (P|ρ_ia)(ρ_ia|Q).
!>
!> \param products Atom-local or atom-pair three-center integrals (μν|P).
!> \param coefficients One atom's AO rows of the molecular coefficient matrix.
!> \param weights Per-orbital weights for this spin channel.
!> \param homo Last occupied state.
!> \param B_PQ unnormalized transition matrix for this spin channel
! **************************************************************************************************
   SUBROUTINE compute_B_PQ_transition(products, coefficients, weights, homo, B_PQ)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: products(:, :, :), coefficients(:, :), &
                                                            weights(:)
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(OUT)                         :: B_PQ(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_B_PQ_transition'

      INTEGER                                            :: b, handle, n, naux, nstates
      REAL(KIND=dp), ALLOCATABLE                         :: empty(:, :), occupied(:, :), scaled(:, :)
      TYPE(local_gemm_ctxt_type)                         :: gemm_context

      CALL timeset(routineN, handle)

      CALL gemm_context%create(LOCAL_GEMM_PU_HOST)
      b = SIZE(coefficients, 1)
      nstates = SIZE(coefficients, 2)
      naux = SIZE(products, 3)
      CPASSERT(SIZE(products, 1) == b .AND. SIZE(products, 2) == b)
      CPASSERT(SIZE(B_PQ, 1) == naux .AND. SIZE(B_PQ, 2) == naux)
      CPASSERT(SIZE(weights) == nstates .AND. homo > 0 .AND. homo < nstates)
      ALLOCATE (scaled(b, nstates), occupied(b, b), empty(b, b))
      DO n = 1, nstates
         scaled(:, n) = coefficients(:, n)*SQRT(weights(n))
      END DO
      CALL gemm_context%gemm( &
         'N', 'T', b, b, homo, 1.0_dp, scaled, b, scaled, b, &
         0.0_dp, occupied, b)
      CALL gemm_context%gemm('N', 'T', b, b, nstates - homo, 1.0_dp, scaled(:, homo + 1:), b, &
                             scaled(:, homo + 1:), b, 0.0_dp, empty, b)
      CALL contract_B_PQ(products, occupied, empty, B_PQ)
      CALL gemm_context%destroy()

      CALL timestop(handle)

   END SUBROUTINE compute_B_PQ_transition

! **************************************************************************************************
!> \brief Computes the frontier-exchange contribution B_PQ^X in Eqs. (4) or (11).
!>
!>        B_PQ^X = Σ_{n∈W_Δ} Σ_i^occ (P|ρ_ni)(ρ_ni|Q).
!>
!> \param products Atom-local or atom-pair three-center integrals (μν|P).
!> \param coefficients Molecular coefficients restricted to this atom's AO rows.
!> \param energies Molecular eigenvalues, occupied first, in Hartree.
!> \param homo Last occupied state in this spin channel.
!> \param delta frontier-window extension in Hartree
!> \param B_PQ unnormalized exchange matrix
! **************************************************************************************************
   SUBROUTINE compute_B_PQ_exchange(products, coefficients, energies, homo, delta, B_PQ)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: products(:, :, :), coefficients(:, :), &
                                                            energies(:)
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(IN)                          :: delta
      REAL(KIND=dp), INTENT(OUT)                         :: B_PQ(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_B_PQ_exchange'

      INTEGER                                            :: b, handle, n, nstates
      REAL(KIND=dp), ALLOCATABLE                         :: frontier(:, :), occupied(:, :), &
                                                            selected(:, :)
      TYPE(local_gemm_ctxt_type)                         :: gemm_context

      CALL timeset(routineN, handle)

      B_PQ(:, :) = 0.0_dp
      nstates = SIZE(energies)
      IF (homo < 1 .OR. homo >= nstates) THEN
         CPABORT("AUTO_RI exchange target requires occupied and empty states")
      END IF
      IF (.NOT. ieee_is_finite(delta) .OR. &
          .NOT. ALL(ieee_is_finite(energies))) THEN
         CPABORT("AUTO_RI exchange target requires finite input")
      END IF
      IF (delta < 0.0_dp .OR. energies(homo + 1) <= energies(homo)) THEN
         CPABORT("AUTO_RI exchange target requires a valid frontier window")
      END IF
      b = SIZE(coefficients, 1)
      CPASSERT(SIZE(coefficients, 2) == nstates)
      ALLOCATE (occupied(b, b), frontier(b, b), selected(b, nstates), source=0.0_dp)
      DO n = 1, nstates
         IF (energies(n) < energies(homo) - delta .OR. &
             energies(n) > energies(homo + 1) + delta) CYCLE
         selected(:, n) = coefficients(:, n)
      END DO
      CALL gemm_context%create(LOCAL_GEMM_PU_HOST)
      ! D_μν = Σ_i^occ C_μi C_νi.
      ! F_μν = Σ_{n: ε_n ∈ [ε_H−δ,ε_L+δ]} C_μn C_νn.
      CALL gemm_context%gemm('N', 'T', b, b, homo, 1.0_dp, coefficients, b, coefficients, b, &
                             0.0_dp, occupied, b)
      CALL gemm_context%gemm('N', 'T', b, b, nstates, 1.0_dp, selected, b, selected, b, &
                             0.0_dp, frontier, b)
      CALL contract_B_PQ(products, occupied, frontier, B_PQ)
      CALL gemm_context%destroy()

      CALL timestop(handle)

   END SUBROUTINE compute_B_PQ_exchange

! **************************************************************************************************
!> \brief Contracts AO density matrices and (μν|P) to form B_PQ.
!>
!>        B_PQ = Σ_μν (μν|P) [left (..|Q) right]_μν.
!>
!> \param products Atom-local or atom-pair three-center integrals (μν|P).
!> \param left First density block, for occupied partners.
!> \param right Second density block, for empty or hard-window target states.
!> \param B_PQ symmetric RI-basis product matrix
! **************************************************************************************************
   SUBROUTINE contract_B_PQ(products, left, right, B_PQ)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: products(:, :, :), left(:, :), &
                                                            right(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: B_PQ(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_B_PQ'

      INTEGER                                            :: b, handle, mu, naux, nu, p, q
      REAL(KIND=dp), ALLOCATABLE                         :: weighted(:, :), work(:, :)
      TYPE(local_gemm_ctxt_type)                         :: gemm_context

      CALL timeset(routineN, handle)

      b = SIZE(products, 1)
      naux = SIZE(products, 3)
      CPASSERT(SIZE(products, 2) == b)
      CPASSERT(SIZE(left, 1) == b .AND. SIZE(left, 2) == b)
      CPASSERT(SIZE(right, 1) == b .AND. SIZE(right, 2) == b)
      CPASSERT(SIZE(B_PQ, 1) == naux .AND. SIZE(B_PQ, 2) == naux)
      ALLOCATE (work(b, b), weighted(b, b))
      CALL gemm_context%create(LOCAL_GEMM_PU_HOST)
      DO q = 1, naux
         ! B_PQ = Σ_μν (μν|P) [left (..|Q) right]_μν.
         CALL gemm_context%gemm( &
            'N', 'N', b, b, b, 1.0_dp, left, b, products(:, :, q), b, &
            0.0_dp, work, b)
         CALL gemm_context%gemm('N', 'N', b, b, b, 1.0_dp, work, b, right, b, 0.0_dp, weighted, b)
         DO p = 1, q
            B_PQ(p, q) = 0.0_dp
            DO nu = 1, b
               DO mu = 1, b
                  B_PQ(p, q) = B_PQ(p, q) + &
                               products(mu, nu, p)*weighted(mu, nu)
               END DO
            END DO
            B_PQ(q, p) = B_PQ(p, q)
         END DO
      END DO
      CALL gemm_context%destroy()

      CALL timestop(handle)

   END SUBROUTINE contract_B_PQ

! **************************************************************************************************
!> \brief Adds separately normalized transition and exchange matrices.
!> \param factor Inverse metric factor, factor factor^T = M^+.
!> \param transition transition matrix, replaced by the combined target
!> \param exchange Unnormalized hard-window exchange covariance.
!> \param normalization number of nonzero normalized terms
! **************************************************************************************************
   SUBROUTINE combine_transition_and_exchange_targets(factor, transition, exchange, normalization)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: factor(:, :)
      REAL(KIND=dp), INTENT(INOUT)                       :: transition(:, :)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: exchange(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: normalization

      CHARACTER(LEN=*), PARAMETER :: routineN = 'combine_transition_and_exchange_targets'

      INTEGER                                            :: handle, k
      REAL(KIND=dp)                                      :: exchange_norm, transition_norm
      REAL(KIND=dp), ALLOCATABLE                         :: work(:, :)

      CALL timeset(routineN, handle)

      ALLOCATE (work(SIZE(factor, 1), SIZE(factor, 2)))
      work(:, :) = MATMUL(transition, factor)
      transition_norm = 0.0_dp
      DO k = 1, SIZE(factor, 2)
         transition_norm = transition_norm + DOT_PRODUCT(factor(:, k), work(:, k))
      END DO
      work(:, :) = MATMUL(exchange, factor)
      exchange_norm = 0.0_dp
      DO k = 1, SIZE(factor, 2)
         exchange_norm = exchange_norm + DOT_PRODUCT(factor(:, k), work(:, k))
      END DO
      normalization = 0.0_dp
      IF (.NOT. ieee_is_finite(transition_norm) .OR. &
          .NOT. ieee_is_finite(exchange_norm) .OR. &
          transition_norm < 0.0_dp .OR. exchange_norm < 0.0_dp) THEN
         CPABORT("AUTO_RI target normalization is invalid")
      END IF
      ! B = B_ia / Tr(M^+ B_ia) + B_ni / Tr(M^+ B_ni), with zero-norm terms omitted.
      IF (transition_norm > 0.0_dp) THEN
         transition(:, :) = transition/transition_norm
         normalization = normalization + 1.0_dp
      ELSE
         transition(:, :) = 0.0_dp
      END IF
      IF (exchange_norm > 0.0_dp) THEN
         transition(:, :) = transition + exchange/exchange_norm
         normalization = normalization + 1.0_dp
      END IF
      IF (normalization <= 0.0_dp .OR. .NOT. ALL(ieee_is_finite(transition))) THEN
         CPABORT("AUTO_RI target normalization is invalid")
      END IF

      CALL timestop(handle)

   END SUBROUTINE combine_transition_and_exchange_targets

! **************************************************************************************************
!> \brief Computes the atom-local matrix B_PQ^A in Eqs. (3)-(6).
!>
!>        B_PQ^A = B_PQ^{A,C}/N_C^A + B_PQ^{A,X}/N_X^A.
!>
!>        The molecular-orbital coefficient rows belonging to one atom are
!>        collected on one MPI rank. No global AO-pair matrix is formed.
!> \param bs_env molecular spectrum and distributed MO coefficients
!> \param iatom Atom whose AO rows are needed.
!> \param root_rank MPI rank performing the atomic contractions.
!> \param products Atom-local three-center integrals (μν|P).
!> \param metric Atom-local two-center RI metric matrix (P|Q).
!> \param B_PQ normalized atom-local matrix B_PQ^A
!> \param normalization Number of nonzero normalized contributions
! **************************************************************************************************
   SUBROUTINE compute_B_PQ_AA(bs_env, iatom, root_rank, products, metric, &
                              B_PQ, normalization)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom, root_rank
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: products(:, :, :), metric(:, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: B_PQ(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: normalization

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_B_PQ_AA'

      INTEGER                                            :: b, column, first, handle, i, j, naux, &
                                                            nstates, rank, row, spin
      INTEGER, POINTER                                   :: cols(:), rows(:)
      REAL(KIND=dp), ALLOCATABLE                         :: exchange(:, :), factor(:, :), &
                                                            slab(:, :), spin_target(:, :), &
                                                            weights(:)

      CALL timeset(routineN, handle)
      NULLIFY (rows, cols)
      b = SIZE(products, 1)
      naux = SIZE(products, 3)
      nstates = bs_env%n_mo_retained
      CPASSERT(ALLOCATED(bs_env%eigenval_scf_Gamma))
      first = 1 + SUM(bs_env%sizes_AO(:iatom - 1))
      ALLOCATE (B_PQ(naux, naux), source=0.0_dp)
      ALLOCATE (slab(b, nstates), weights(nstates), spin_target(naux, naux))
      normalization = 0.0_dp
      ALLOCATE (exchange(naux, naux), source=0.0_dp)
      DO spin = 1, bs_env%n_spin
         CALL compute_orbital_weights(bs_env%eigenval_scf_Gamma(1:nstates, spin), &
                                      bs_env%n_occ(spin), weights)
         CALL cp_fm_get_info(bs_env%fm_mo_coeff_Gamma(spin), row_indices=rows, col_indices=cols)
         slab(:, :) = 0.0_dp
         DO i = 1, SIZE(rows)
            row = rows(i) - first + 1
            IF (row < 1 .OR. row > b) CYCLE
            DO j = 1, SIZE(cols)
               column = cols(j)
               IF (column > nstates) CYCLE
               slab(row, column) = bs_env%fm_mo_coeff_Gamma(spin)%local_data(i, j)
            END DO
         END DO
         CALL bs_env%para_env%sum(slab, root_rank)
         IF (bs_env%para_env%mepos /= root_rank) CYCLE
         CALL compute_B_PQ_transition(products, slab, weights, bs_env%n_occ(spin), spin_target)
         B_PQ(:, :) = B_PQ + spin_target
         CALL compute_B_PQ_exchange(products, slab, &
                                    bs_env%eigenval_scf_Gamma(1:nstates, spin), &
                                    bs_env%n_occ(spin), bs_env%auto_ri%delta, spin_target)
         exchange(:, :) = exchange + spin_target
      END DO
      IF (bs_env%para_env%mepos /= root_rank) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL compute_metric_inverse_factor( &
         metric, factor, rank, bs_env%auto_ri%metric_null_threshold, &
         bs_env%auto_ri%psd_tolerance)
      CALL combine_transition_and_exchange_targets(factor, B_PQ, exchange, &
                                                   normalization)
      CALL timestop(handle)

   END SUBROUTINE compute_B_PQ_AA

! **************************************************************************************************
!> \brief Normalizes each contraction coefficient vector U_Pp.
!>
!>        U_Pp <- U_Pp/sqrt(Σ_P |U_Pp|^2).
!> \param coefficient ...
! **************************************************************************************************
   SUBROUTINE normalize_RI_coefficients(coefficient)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: coefficient

      CHARACTER(LEN=*), PARAMETER :: routineN = 'normalize_RI_coefficients'

      INTEGER                                            :: handle, j
      REAL(KIND=dp)                                      :: column_norm

      CALL timeset(routineN, handle)

      DO j = 1, SIZE(coefficient, 2)
         column_norm = NORM2(coefficient(:, j))
         IF (column_norm <= TINY(1.0_dp)) THEN
            CPABORT("AUTO_RI produced a numerically zero localized column")
         END IF
         coefficient(:, j) = coefficient(:, j)/column_norm
      END DO

      CALL timestop(handle)

   END SUBROUTINE normalize_RI_coefficients

! **************************************************************************************************
!> \brief Creates a matrix in the existing molecular BLACS context.
!> \param qs_env ...
!> \param matrix ...
!> \param nr ...
!> \param nc ...
! **************************************************************************************************
   SUBROUTINE create_distributed_matrix(qs_env, matrix, nr, nc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(OUT)                      :: matrix
      INTEGER, INTENT(IN)                                :: nr, nc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_distributed_matrix'

      INTEGER                                            :: handle
      TYPE(cp_blacs_env_type), POINTER                   :: blacs
      TYPE(cp_fm_struct_type), POINTER                   :: layout
      TYPE(mp_para_env_type), POINTER                    :: para

      CALL timeset(routineN, handle)

      NULLIFY (blacs, para, layout)
      CALL get_qs_env(qs_env, blacs_env=blacs, para_env=para)
      CALL cp_fm_struct_create(layout, context=blacs, nrow_global=nr, ncol_global=nc, para_env=para)
      CALL cp_fm_create(matrix, layout, set_zero=.TRUE.)
      CALL cp_fm_struct_release(layout)

      CALL timestop(handle)

   END SUBROUTINE create_distributed_matrix

! **************************************************************************************************
!> \brief Computes and stores the two-center contractions U^AB in Eqs. (8)-(13).
!>
!>        The selected AA space is projected out according to Eq. (8). The symmetric
!>        pair density in Eq. (9) defines B^{AB,C} and B^{AB,X} in Eqs. (10) and (11).
!>        Diagonalizing their sum gives U^AB in Eq. (12), which expands the contracted
!>        pair functions in Eq. (13).
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_AB_contractions(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_AB_contractions'

      INTEGER :: atomic_columns, center_index, current_atom, expected_rank, first_pair, &
         fit_coefficient_size, handle, handle_pair_batches, iatom, jatom, k, last_pair, n, &
         n_RI_requested, na, nao, nb, nlocal, npair, nret, pair_columns, pair_storage_size, rank, &
         root_rank, row, spin
      INTEGER, ALLOCATABLE :: aa_sizes(:), block_next(:), counts(:), pair_atom_a(:), &
         pair_atom_b(:), pair_coefficient_offset(:), pair_n_to_a(:), sizes(:)
      REAL(KIND=dp), ALLOCATABLE :: atomic(:, :), B_PQ(:, :), c(:, :), eigenvalues(:), &
         metric(:, :), pair_coefficient(:), pair_mos(:, :, :), scores(:, :), slab(:, :)
      TYPE(gw_3c_ctx_type)                               :: ctx
      TYPE(gw_3c_ws_type)                                :: ws
      TYPE(mp_para_env_type), POINTER                    :: para

      CALL timeset(routineN, handle)

      NULLIFY (para)
      CALL get_qs_env(qs_env, para_env=para)
      ALLOCATE (sizes(bs_env%n_atom), aa_sizes(bs_env%n_atom))
      DO iatom = 1, bs_env%n_atom
         sizes(iatom) = reference_RI_size(bs_env, iatom)
      END DO
      aa_sizes(:) = bs_env%auto_ri%optimized_size
      n = SUM(sizes)
      atomic_columns = SUM(aa_sizes)

      npair = 0
      DO iatom = 1, bs_env%n_atom - 1
         DO jatom = iatom + 1, bs_env%n_atom
            IF (NORM2(bs_env%ri_rs%particle_set(iatom)%r - bs_env%ri_rs%particle_set(jatom)%r) <= &
                bs_env%auto_ri%neighbor_radius) THEN
               npair = npair + 1
            END IF
         END DO
      END DO
      ALLOCATE (pair_atom_a(npair), pair_atom_b(npair), pair_n_to_a(npair), &
                pair_coefficient_offset(npair), counts(npair))
      pair_n_to_a = 0
      pair_coefficient_offset = 0
      ALLOCATE (scores(2*MAXVAL(sizes), npair), source=0.0_dp)

      CALL gw_3c_ctx_create(ctx, qs_env, bs_env%ri_metric, bs_env%basis_set_AO, &
                            bs_env%basis_set_AO, bs_env%basis_set_RI)
      CALL gw_3c_ws_create(ws, ctx)
      k = 0
      DO iatom = 1, bs_env%n_atom - 1
         DO jatom = iatom + 1, bs_env%n_atom
            IF (NORM2(bs_env%ri_rs%particle_set(iatom)%r - bs_env%ri_rs%particle_set(jatom)%r) > &
                bs_env%auto_ri%neighbor_radius) CYCLE
            k = k + 1
            pair_atom_a(k) = iatom
            pair_atom_b(k) = jatom
         END DO
      END DO
      pair_storage_size = 0
      DO k = 1, npair
         root_rank = MODULO(k - 1, para%num_pe)
         IF (para%mepos /= root_rank) CYCLE
         nlocal = sizes(pair_atom_a(k)) + sizes(pair_atom_b(k))
         pair_coefficient_offset(k) = pair_storage_size + 1
         pair_storage_size = pair_storage_size + nlocal*nlocal
      END DO
      ALLOCATE (pair_coefficient(pair_storage_size))
      CALL timeset(routineN//'_pair_batches', handle_pair_batches)
      DO first_pair = 1, npair, para%num_pe
         last_pair = MIN(npair, first_pair + para%num_pe - 1)
         ! Gather only one pair's MO rows per rank, then solve all batch pairs concurrently.
         DO k = first_pair, last_pair
            root_rank = k - first_pair
            iatom = pair_atom_a(k)
            jatom = pair_atom_b(k)
            IF (para%mepos == root_rank) THEN
               ALLOCATE (pair_mos( &
                         bs_env%sizes_AO(iatom) + bs_env%sizes_AO(jatom), &
                         bs_env%n_mo_retained, bs_env%n_spin))
            END IF
            DO spin = 1, bs_env%n_spin
               row = 0
               DO center_index = 1, 2
                  current_atom = iatom
                  IF (center_index == 2) current_atom = jatom
                  nao = bs_env%sizes_AO(current_atom)
                  ALLOCATE (slab(nao, bs_env%n_mo_retained))
                  CALL cp_fm_get_submatrix(bs_env%fm_mo_coeff_Gamma(spin), slab, &
                                           start_row=1 + SUM(bs_env%sizes_AO(:current_atom - 1)))
                  IF (para%mepos == root_rank) pair_mos(row + 1:row + nao, :, spin) = slab
                  DEALLOCATE (slab)
                  row = row + nao
               END DO
            END DO
         END DO
         k = first_pair + para%mepos
         IF (k <= last_pair) THEN
            iatom = pair_atom_a(k)
            jatom = pair_atom_b(k)
            na = sizes(iatom)
            nb = sizes(jatom)
            nlocal = na + nb
            ! Eqs. (9)-(11): form the symmetric pair density and B^{AB,C}+B^{AB,X}.
            CALL compute_B_PQ_AB(bs_env, ctx, ws, iatom, jatom, pair_mos, metric, B_PQ)
            ALLOCATE (atomic(nlocal, aa_sizes(iatom) + aa_sizes(jatom)), source=0.0_dp)
            atomic(1:na, 1:aa_sizes(iatom)) = &
               bs_env%auto_ri%transform( &
               1:na, 1:aa_sizes(iatom), iatom)
            atomic(na + 1:, aa_sizes(iatom) + 1:) = &
               bs_env%auto_ri%transform( &
               1:nb, 1:aa_sizes(jatom), jatom)
            ! Eqs. (8) and (12): project out the AA space and diagonalize the remainder.
            CALL solve_AB_eigenproblem( &
               metric, B_PQ, atomic, bs_env%auto_ri%min_score_ratio, &
               bs_env%auto_ri%pair_absolute_threshold, &
               bs_env%auto_ri%metric_null_threshold, bs_env%auto_ri%psd_tolerance, &
               c, eigenvalues)
            scores(1:SIZE(eigenvalues), k) = eigenvalues
            DEALLOCATE (eigenvalues)
            pair_coefficient( &
               pair_coefficient_offset(k):pair_coefficient_offset(k) + SIZE(c) - 1) = &
               RESHAPE(c, [SIZE(c)])
            DEALLOCATE (c)
            DEALLOCATE (atomic, B_PQ, metric, pair_mos)
         END IF
      END DO
      CALL timestop(handle_pair_batches)
      k = npair
      CALL gw_3c_ws_release(ws)
      CALL gw_3c_ctx_release(ctx)
      n_RI_requested = MAX(bs_env%n_atom, NINT(MIN(bs_env%auto_ri%ri_ao_ratio, &
                                                   REAL(n, dp)/REAL(SUM(bs_env%sizes_AO), dp))* &
                                               REAL(SUM(bs_env%sizes_AO), dp)))
      CALL para%sum(scores)
      CALL select_RI_counts(scores, n_RI_requested - SUM(aa_sizes), 0, counts)
      pair_columns = SUM(counts)
      bs_env%auto_ri%optimized_size(:) = aa_sizes
      DO k = 1, npair
         iatom = pair_atom_a(k)
         jatom = pair_atom_b(k)
         nret = counts(k)
         pair_n_to_a(k) = balance_AB_columns(bs_env%auto_ri%optimized_size(iatom), &
                                             bs_env%auto_ri%optimized_size(jatom), nret)
         bs_env%auto_ri%optimized_size(iatom) = &
            bs_env%auto_ri%optimized_size(iatom) + pair_n_to_a(k)
         bs_env%auto_ri%optimized_size(jatom) = &
            bs_env%auto_ri%optimized_size(jatom) + nret - pair_n_to_a(k)
      END DO

      rank = SUM(bs_env%auto_ri%optimized_size)
      expected_rank = atomic_columns + pair_columns
      CPASSERT(rank == expected_rank)
      CPASSERT(ALL(bs_env%auto_ri%optimized_size > 0))
      IF (rank > n) THEN
         CPABORT("AUTO_RI request exceeds independent reference RI directions")
      END IF
      fit_coefficient_size = 0
      DO iatom = 1, bs_env%n_atom
         fit_coefficient_size = fit_coefficient_size + sizes(iatom)*aa_sizes(iatom)
      END DO
      DO k = 1, npair
         fit_coefficient_size = fit_coefficient_size + &
                                (sizes(pair_atom_a(k)) + sizes(pair_atom_b(k)))*counts(k)
      END DO
      ALLOCATE (bs_env%auto_ri%fit_group_atom_a(bs_env%n_atom + npair), source=0)
      ALLOCATE (bs_env%auto_ri%fit_group_atom_b(bs_env%n_atom + npair), source=0)
      ALLOCATE (bs_env%auto_ri%fit_group_first_a(bs_env%n_atom + npair), source=0)
      ALLOCATE (bs_env%auto_ri%fit_group_first_b(bs_env%n_atom + npair), source=0)
      ALLOCATE (bs_env%auto_ri%fit_group_n_to_a(bs_env%n_atom + npair), source=0)
      ALLOCATE (bs_env%auto_ri%fit_group_nrow(bs_env%n_atom + npair), source=0)
      ALLOCATE (bs_env%auto_ri%fit_group_ncolumn(bs_env%n_atom + npair), source=0)
      ALLOCATE (bs_env%auto_ri%fit_group_offset(bs_env%n_atom + npair), source=0)
      ALLOCATE (bs_env%auto_ri%fit_group_coefficient(fit_coefficient_size), source=0.0_dp)
      bs_env%auto_ri%fit_group_count = 0

      ALLOCATE (block_next(bs_env%n_atom))
      block_next(1) = 1
      DO iatom = 2, bs_env%n_atom
         block_next(iatom) = block_next(iatom - 1) + bs_env%auto_ri%optimized_size(iatom - 1)
      END DO
      DO iatom = 1, bs_env%n_atom
         IF (aa_sizes(iatom) > 0) bs_env%auto_ri%fit_group_count = &
            bs_env%auto_ri%fit_group_count + 1
         CALL store_AA_columns(bs_env, iatom, aa_sizes(iatom), block_next)
      END DO

      DO k = 1, npair
         iatom = pair_atom_a(k)
         jatom = pair_atom_b(k)
         nret = counts(k)
         IF (nret == 0) CYCLE
         root_rank = MODULO(k - 1, para%num_pe)
         nlocal = sizes(iatom) + sizes(jatom)
         ALLOCATE (c(nlocal, nret))
         IF (para%mepos == root_rank) THEN
            c(:, :) = RESHAPE( &
                      pair_coefficient(pair_coefficient_offset(k): &
                                       pair_coefficient_offset(k) + nlocal*nret - 1), &
                      [nlocal, nret])
         END IF
         CALL para%bcast(c, root_rank)
         bs_env%auto_ri%fit_group_count = bs_env%auto_ri%fit_group_count + 1
         ! Eq. (13): store the selected neighbor-pair contractions in the global transform.
         CALL store_AB_columns(bs_env, c, iatom, jatom, pair_n_to_a(k), block_next)
         DEALLOCATE (c)
      END DO

      DO iatom = 1, bs_env%n_atom
         expected_rank = 1 + SUM(bs_env%auto_ri%optimized_size(:iatom))
         CPASSERT(block_next(iatom) == expected_rank)
      END DO
      DEALLOCATE (block_next, pair_atom_a, pair_atom_b, pair_n_to_a, &
                  pair_coefficient_offset, pair_coefficient)

      CALL timestop(handle)

   END SUBROUTINE compute_AB_contractions

! **************************************************************************************************
!> \brief Divides AB columns between atoms A and B for load balance.
!> \param load_a current number of columns owned by A
!> \param load_b current number of columns owned by B
!> \param ncolumn number of new pair columns
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION balance_AB_columns(load_a, load_b, ncolumn) RESULT(n_to_a)
      INTEGER, INTENT(IN)                                :: load_a, load_b, ncolumn

      n_to_a = (load_b + ncolumn - load_a + 1)/2
      n_to_a = MAX(0, MIN(ncolumn, n_to_a))

   END FUNCTION balance_AB_columns

! **************************************************************************************************
!> \brief Stores the AA contraction columns for one atom.
!> \param bs_env ...
!> \param atom ...
!> \param ncolumn ...
!> \param block_next ...
! **************************************************************************************************
   SUBROUTINE store_AA_columns(bs_env, atom, ncolumn, block_next)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom, ncolumn
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: block_next

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'store_AA_columns'

      INTEGER                                            :: first, group, handle, n_reference
      REAL(KIND=dp), ALLOCATABLE                         :: coefficient(:, :)

      CALL timeset(routineN, handle)

      first = block_next(atom)
      group = bs_env%auto_ri%fit_group_count
      n_reference = reference_RI_size(bs_env, atom)
      ALLOCATE (coefficient(n_reference, ncolumn))
      IF (ncolumn > 0) THEN
         coefficient(:, :) = bs_env%auto_ri%transform(1:n_reference, 1:ncolumn, atom)
         CALL normalize_RI_coefficients(coefficient)
         bs_env%auto_ri%fit_group_atom_a(group) = atom
         bs_env%auto_ri%fit_group_first_a(group) = &
            first - SUM(bs_env%auto_ri%optimized_size(:atom - 1))
         bs_env%auto_ri%fit_group_n_to_a(group) = ncolumn
         CALL store_fit_group_coefficients(bs_env%auto_ri, group, coefficient)
      END IF
      DEALLOCATE (coefficient)
      block_next(atom) = first + ncolumn

      CALL timestop(handle)

   END SUBROUTINE store_AA_columns

! **************************************************************************************************
!> \brief Stores an AB contraction in the blocks assigned to atoms A and B.
!> \param bs_env ...
!> \param coefficient ...
!> \param atom_a ...
!> \param atom_b ...
!> \param n_to_a ...
!> \param block_next ...
! **************************************************************************************************
   SUBROUTINE store_AB_columns(bs_env, coefficient, atom_a, atom_b, n_to_a, block_next)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: coefficient
      INTEGER, INTENT(IN)                                :: atom_a, atom_b, n_to_a
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: block_next

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'store_AB_columns'

      INTEGER                                            :: first_a, first_b, group, handle, ncolumn

      CALL timeset(routineN, handle)

      ncolumn = SIZE(coefficient, 2)
      group = bs_env%auto_ri%fit_group_count
      first_a = block_next(atom_a)
      first_b = block_next(atom_b)
      CALL normalize_RI_coefficients(coefficient)

      bs_env%auto_ri%fit_group_atom_a(group) = atom_a
      bs_env%auto_ri%fit_group_atom_b(group) = atom_b
      bs_env%auto_ri%fit_group_first_a(group) = &
         first_a - SUM(bs_env%auto_ri%optimized_size(:atom_a - 1))
      bs_env%auto_ri%fit_group_first_b(group) = &
         first_b - SUM(bs_env%auto_ri%optimized_size(:atom_b - 1))
      bs_env%auto_ri%fit_group_n_to_a(group) = n_to_a
      CALL store_fit_group_coefficients(bs_env%auto_ri, group, coefficient)
      block_next(atom_a) = first_a + n_to_a
      block_next(atom_b) = first_b + ncolumn - n_to_a

      CALL timestop(handle)

   END SUBROUTINE store_AB_columns

! **************************************************************************************************
!> \brief Packs one AA or AB coefficient matrix into the persistent AUTO_RI data.
!> \param auto_ri ...
!> \param group ...
!> \param coefficient ...
! **************************************************************************************************
   SUBROUTINE store_fit_group_coefficients(auto_ri, group, coefficient)
      TYPE(auto_ri_type), INTENT(INOUT)                  :: auto_ri
      INTEGER, INTENT(IN)                                :: group
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: coefficient

      INTEGER                                            :: first, last

      first = 1
      IF (group > 1) THEN
         first = auto_ri%fit_group_offset(group - 1) + &
                 auto_ri%fit_group_nrow(group - 1)*auto_ri%fit_group_ncolumn(group - 1)
      END IF
      last = first + SIZE(coefficient) - 1
      CPASSERT(last <= SIZE(auto_ri%fit_group_coefficient))
      auto_ri%fit_group_offset(group) = first
      auto_ri%fit_group_nrow(group) = SIZE(coefficient, 1)
      auto_ri%fit_group_ncolumn(group) = SIZE(coefficient, 2)
      auto_ri%fit_group_coefficient(first:last) = RESHAPE(coefficient, [SIZE(coefficient)])

   END SUBROUTINE store_fit_group_coefficients

! **************************************************************************************************
!> \brief Computes the two-center matrix B_PQ^AB in Eqs. (9)-(12).
!>
!>        Its numerator uses the symmetric off-diagonal density
!>        ρ̅_mn^AB=ρ_mn^AB+ρ_mn^BA after projection from the retained AA
!>        space. Its normalization uses the complete two-atom density.
!> \param bs_env ...
!> \param ctx ...
!> \param ws ...
!> \param iatom ...
!> \param jatom ...
!> \param c ...
!> \param metric ...
!> \param B_PQ normalized two-center matrix B_PQ^AB
! **************************************************************************************************
   SUBROUTINE compute_B_PQ_AB(bs_env, ctx, ws, iatom, jatom, c, metric, B_PQ)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      TYPE(gw_3c_ws_type), INTENT(INOUT)                 :: ws
      INTEGER, INTENT(IN)                                :: iatom, jatom
      REAL(KIND=dp), INTENT(IN)                          :: c(:, :, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: metric(:, :), B_PQ(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_B_PQ_AB'

      INTEGER                                            :: ao_off(2), ao_size(2), atoms(2), handle, &
                                                            homo, i, j, k, kind_a, kind_b, n, nr, &
                                                            rank, ri_off(2), ri_size(2), spin
      LOGICAL                                            :: screened
      REAL(KIND=dp)                                      :: nex, nrsp
      REAL(KIND=dp), ALLOCATABLE :: block(:, :), exchange(:, :), full_exchange(:, :), &
         full_transition(:, :), pair_products(:, :, :), products(:, :, :), transition(:, :), &
         weights(:), work(:, :), x(:, :)

      CALL timeset(routineN, handle)

      atoms(1) = iatom
      atoms(2) = jatom
      ao_size = bs_env%sizes_AO(atoms)
      DO i = 1, 2
         ri_size(i) = reference_RI_size(bs_env, atoms(i))
      END DO
      ao_off(1) = 0
      ao_off(2) = ao_size(1)
      ri_off(1) = 0
      ri_off(2) = ri_size(1)
      n = SUM(ao_size)
      nr = SUM(ri_size)
      ALLOCATE (B_PQ(nr, nr), metric(nr, nr), source=0.0_dp)
      ALLOCATE (products(n, n, nr), source=0.0_dp)
      DO i = 1, 2
         DO j = 1, 2
            CALL compute_AB_RI_metric_matrix(bs_env, atoms(i), atoms(j), block, bs_env%ri_metric)
            metric(ri_off(i) + 1:ri_off(i) + ri_size(i), &
                   ri_off(j) + 1:ri_off(j) + ri_size(j)) = block
            DEALLOCATE (block)
            DO k = 1, 2
               CALL build_3c_integral_block_ctx( &
                  products, ctx, ws, atom_j=atoms(j), atom_k=atoms(k), &
                  atom_i=atoms(i), j_offset=ao_off(j), k_offset=ao_off(k), &
                  i_offset=ri_off(i), screened=screened)
            END DO
         END DO
      END DO
      kind_a = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
      kind_b = bs_env%ri_rs%particle_set(jatom)%atomic_kind%kind_number
      CALL add_RI_metric_regularization(bs_env%basis_set_RI(kind_a)%gto_basis_set, &
                                        bs_env%regularization_RI, metric)
      CALL add_RI_metric_regularization(bs_env%basis_set_RI(kind_b)%gto_basis_set, &
                                        bs_env%regularization_RI, metric, offset=ri_size(1))
      ALLOCATE (pair_products, source=products)
      pair_products(1:ao_size(1), 1:ao_size(1), :) = 0.0_dp
      pair_products(ao_size(1) + 1:, ao_size(1) + 1:, :) = 0.0_dp
      ALLOCATE (weights(bs_env%n_mo_retained), work(nr, nr))
      ALLOCATE (transition(nr, nr), exchange(nr, nr), &
                full_transition(nr, nr), full_exchange(nr, nr), &
                source=0.0_dp)
      DO spin = 1, bs_env%n_spin
         homo = bs_env%n_occ(spin)
         CALL compute_orbital_weights(energies= &
                                      bs_env%eigenval_scf_Gamma(1:bs_env%n_mo_retained, spin), &
                                      homo=homo, weights=weights)
         CALL compute_B_PQ_transition(products, c(:, :, spin), weights, homo, work)
         full_transition(:, :) = full_transition + work
         CALL compute_B_PQ_transition(pair_products, c(:, :, spin), weights, homo, work)
         transition(:, :) = transition + work
         CALL compute_B_PQ_exchange(products, c(:, :, spin), &
                                    bs_env%eigenval_scf_Gamma(1:bs_env%n_mo_retained, spin), &
                                    homo, bs_env%auto_ri%delta, work)
         full_exchange(:, :) = full_exchange + work
         CALL compute_B_PQ_exchange(pair_products, c(:, :, spin), &
                                    bs_env%eigenval_scf_Gamma(1:bs_env%n_mo_retained, spin), &
                                    homo, bs_env%auto_ri%delta, work)
         exchange(:, :) = exchange + work
      END DO
      CALL compute_metric_inverse_factor( &
         metric, x, rank, bs_env%auto_ri%metric_null_threshold, &
         bs_env%auto_ri%psd_tolerance)
      DEALLOCATE (work)
      ALLOCATE (work(nr, rank))
      work(:, :) = MATMUL(full_transition, x)
      nrsp = SUM(x*work)
      work(:, :) = MATMUL(full_exchange, x)
      nex = SUM(x*work)
      ! Eqs. (10)-(12): B_PQ^AB = (B_PQ^{AB,C}/N_C^{AB}+B_PQ^{AB,X}/N_X^{AB})/2.
      IF (nrsp > 0.0_dp) B_PQ(:, :) = B_PQ + 0.5_dp*transition/nrsp
      IF (nex > 0.0_dp) B_PQ(:, :) = B_PQ + 0.5_dp*exchange/nex

      CALL timestop(handle)

   END SUBROUTINE compute_B_PQ_AB

! **************************************************************************************************
!> \brief Builds the lists of fit groups whose columns belong to each optimized atom block.
!> \param bs_env ...
!> \param lists ...
!> \param max_reference_recomputation upper bound for recomputing one reference atom-pair block
! **************************************************************************************************
   SUBROUTINE build_fit_group_lists(bs_env, lists, max_reference_recomputation)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(fit_group_list_type), ALLOCATABLE, &
         INTENT(OUT)                                     :: lists(:)
      INTEGER, INTENT(OUT)                               :: max_reference_recomputation

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_fit_group_lists'

      INTEGER                                            :: atom_a, atom_b, group, handle, n_to_a, &
                                                            ncolumn, position
      INTEGER, ALLOCATABLE                               :: counts(:), reference_support_count(:)

      CALL timeset(routineN, handle)

      ALLOCATE (lists(bs_env%n_atom))
      ALLOCATE (counts(bs_env%n_atom), source=0)
      ALLOCATE (reference_support_count(bs_env%n_atom), source=1)
      DO group = 1, bs_env%auto_ri%fit_group_count
         atom_a = bs_env%auto_ri%fit_group_atom_a(group)
         atom_b = bs_env%auto_ri%fit_group_atom_b(group)
         ncolumn = bs_env%auto_ri%fit_group_ncolumn(group)
         n_to_a = bs_env%auto_ri%fit_group_n_to_a(group)
         IF (n_to_a > 0) counts(atom_a) = counts(atom_a) + 1
         IF (atom_b > 0 .AND. n_to_a < ncolumn) counts(atom_b) = counts(atom_b) + 1
         IF (atom_b > 0 .AND. n_to_a > 0) THEN
            reference_support_count(atom_b) = reference_support_count(atom_b) + 1
         END IF
         IF (atom_b > 0 .AND. n_to_a < ncolumn) THEN
            reference_support_count(atom_a) = reference_support_count(atom_a) + 1
         END IF
      END DO
      max_reference_recomputation = MAXVAL(reference_support_count)**2

      DO atom_a = 1, bs_env%n_atom
         ALLOCATE (lists(atom_a)%group(counts(atom_a)))
         ALLOCATE (lists(atom_a)%first_group_column(counts(atom_a)))
         ALLOCATE (lists(atom_a)%first_block_column(counts(atom_a)))
         ALLOCATE (lists(atom_a)%ncolumn(counts(atom_a)))
      END DO
      counts = 0
      DO group = 1, bs_env%auto_ri%fit_group_count
         atom_a = bs_env%auto_ri%fit_group_atom_a(group)
         atom_b = bs_env%auto_ri%fit_group_atom_b(group)
         ncolumn = bs_env%auto_ri%fit_group_ncolumn(group)
         n_to_a = bs_env%auto_ri%fit_group_n_to_a(group)
         IF (n_to_a > 0) THEN
            position = counts(atom_a) + 1
            counts(atom_a) = position
            lists(atom_a)%group(position) = group
            lists(atom_a)%first_group_column(position) = 1
            lists(atom_a)%first_block_column(position) = &
               bs_env%auto_ri%fit_group_first_a(group)
            lists(atom_a)%ncolumn(position) = n_to_a
         END IF
         IF (atom_b > 0 .AND. n_to_a < ncolumn) THEN
            position = counts(atom_b) + 1
            counts(atom_b) = position
            lists(atom_b)%group(position) = group
            lists(atom_b)%first_group_column(position) = n_to_a + 1
            lists(atom_b)%first_block_column(position) = &
               bs_env%auto_ri%fit_group_first_b(group)
            lists(atom_b)%ncolumn(position) = ncolumn - n_to_a
         END IF
      END DO
      DO atom_a = 1, bs_env%n_atom
         IF (SUM(lists(atom_a)%ncolumn) /= bs_env%auto_ri%optimized_size(atom_a)) THEN
            CPABORT("AUTO_RI fit-group columns do not fill the optimized atom block")
         END IF
      END DO
      DEALLOCATE (counts, reference_support_count)

      CALL timestop(handle)

   END SUBROUTINE build_fit_group_lists

! **************************************************************************************************
!> \brief Returns the location of atom in atoms(:count), or zero if it is absent.
!> \param atom ...
!> \param atoms ...
!> \param count ...
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION find_atom_position(atom, atoms, count) RESULT(position)
      INTEGER, INTENT(IN)                                :: atom
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atoms
      INTEGER, INTENT(IN)                                :: count

      INTEGER                                            :: i

      position = 0
      DO i = 1, count
         IF (atoms(i) == atom) THEN
            position = i
            RETURN
         END IF
      END DO

   END FUNCTION find_atom_position

! **************************************************************************************************
!> \brief Collects the sparse reference coefficients for one optimized atom block.
!> \param bs_env ...
!> \param atom ...
!> \param list ...
!> \param transform ...
! **************************************************************************************************
   SUBROUTINE build_optimized_atom_transform(bs_env, atom, list, transform)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom
      TYPE(fit_group_list_type), INTENT(IN)              :: list
      TYPE(optimized_atom_transform_type), INTENT(OUT)   :: transform

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_optimized_atom_transform'

      INTEGER :: atom_a, atom_b, first, first_block, first_group, first_row, group, handle, i, &
         last, ncolumn, nrow, nrow_a, nsupport, position
      INTEGER, ALLOCATABLE                               :: support_atoms(:)
      REAL(KIND=dp), ALLOCATABLE                         :: group_coefficient(:, :)

      CALL timeset(routineN, handle)

      ALLOCATE (support_atoms(2*SIZE(list%group)))
      nsupport = 0
      DO i = 1, SIZE(list%group)
         group = list%group(i)
         atom_a = bs_env%auto_ri%fit_group_atom_a(group)
         atom_b = bs_env%auto_ri%fit_group_atom_b(group)
         position = find_atom_position(atom_a, support_atoms, nsupport)
         IF (position == 0) THEN
            nsupport = nsupport + 1
            support_atoms(nsupport) = atom_a
         END IF
         IF (atom_b > 0) THEN
            position = find_atom_position(atom_b, support_atoms, nsupport)
            IF (position == 0) THEN
               nsupport = nsupport + 1
               support_atoms(nsupport) = atom_b
            END IF
         END IF
      END DO
      CPASSERT(nsupport > 0)
      ALLOCATE (transform%reference_atom(nsupport), transform%first_row(nsupport))
      transform%reference_atom(:) = support_atoms(:nsupport)
      first_row = 1
      DO i = 1, nsupport
         transform%first_row(i) = first_row
         first_row = first_row + reference_RI_size(bs_env, transform%reference_atom(i))
      END DO
      ALLOCATE (transform%coefficient(first_row - 1, bs_env%auto_ri%optimized_size(atom)), &
                source=0.0_dp)

      DO i = 1, SIZE(list%group)
         group = list%group(i)
         atom_a = bs_env%auto_ri%fit_group_atom_a(group)
         atom_b = bs_env%auto_ri%fit_group_atom_b(group)
         nrow = bs_env%auto_ri%fit_group_nrow(group)
         ncolumn = bs_env%auto_ri%fit_group_ncolumn(group)
         first = bs_env%auto_ri%fit_group_offset(group)
         last = first + nrow*ncolumn - 1
         ALLOCATE (group_coefficient(nrow, ncolumn))
         group_coefficient(:, :) = RESHAPE( &
                                   bs_env%auto_ri%fit_group_coefficient(first:last), &
                                   [nrow, ncolumn])
         first_group = list%first_group_column(i)
         first_block = list%first_block_column(i)
         ncolumn = list%ncolumn(i)
         IF (first_group < 1 .OR. &
             first_group + ncolumn - 1 > SIZE(group_coefficient, 2)) THEN
            CPABORT("AUTO_RI fit-group column range is invalid")
         END IF
         IF (first_block < 1 .OR. &
             first_block + ncolumn - 1 > bs_env%auto_ri%optimized_size(atom)) THEN
            CPABORT("AUTO_RI optimized atom-block column range is invalid")
         END IF

         nrow_a = reference_RI_size(bs_env, atom_a)
         position = find_atom_position(atom_a, transform%reference_atom, nsupport)
         CPASSERT(position > 0)
         first_row = transform%first_row(position)
         transform%coefficient(first_row:first_row + nrow_a - 1, &
                               first_block:first_block + ncolumn - 1) = &
            group_coefficient(1:nrow_a, first_group:first_group + ncolumn - 1)
         IF (atom_b > 0) THEN
            position = find_atom_position(atom_b, transform%reference_atom, nsupport)
            CPASSERT(position > 0)
            first_row = transform%first_row(position)
            transform%coefficient(first_row:first_row + reference_RI_size(bs_env, atom_b) - 1, &
                                  first_block:first_block + ncolumn - 1) = &
               group_coefficient(nrow_a + 1:, first_group:first_group + ncolumn - 1)
         END IF
         DEALLOCATE (group_coefficient)
      END DO
      DEALLOCATE (support_atoms)

      CALL timestop(handle)

   END SUBROUTINE build_optimized_atom_transform

! **************************************************************************************************
!> \brief Adds REGULARIZATION_RI to one atom-local reference metric block.
!> \param bs_env ...
!> \param atom ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE regularize_reference_RI_metric_block(bs_env, atom, block)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: block

      CHARACTER(LEN=*), PARAMETER :: routineN = 'regularize_reference_RI_metric_block'

      INTEGER                                            :: handle, i, iset, kind, nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(block, 1) == SIZE(block, 2))
      kind = bs_env%ri_rs%particle_set(atom)%atomic_kind%kind_number
      CALL get_gto_basis_set(bs_env%basis_set_RI(kind)%gto_basis_set, &
                             first_sgf=first_sgf, npgf=npgf, nset=nset, &
                             nsgf_set=nsgf_set, zet=zet)
      DO iset = 1, nset
         DO i = first_sgf(1, iset), first_sgf(1, iset) + nsgf_set(iset) - 1
            block(i, i) = block(i, i) + &
                          bs_env%regularization_RI* &
                          MAX(1.0_dp, 1.0_dp/MINVAL(zet(1:npgf(iset), iset)))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE regularize_reference_RI_metric_block

! **************************************************************************************************
!> \brief Directly contracts one optimized atom-pair metric block.
!>
!>        M_pq^{AB} = Σ_IJ (U_Ip^A)^T M_IJ U_Jq^B. Each small reference
!>        atom-pair block M_IJ is discarded immediately after this contraction.
!> \param bs_env ...
!> \param transform_a ...
!> \param transform_b ...
!> \param coulomb true for Coulomb; false for the GW fitting operator
!> \param symmetric true when both transforms describe the same optimized block
!> \param block ...
!> \param reference_block_count ...
! **************************************************************************************************
   SUBROUTINE compute_optimized_RI_metric_block( &
      bs_env, transform_a, transform_b, coulomb, symmetric, block, reference_block_count)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(optimized_atom_transform_type), INTENT(IN)    :: transform_a, transform_b
      LOGICAL, INTENT(IN)                                :: coulomb, symmetric
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: block
      INTEGER, INTENT(INOUT)                             :: reference_block_count

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_optimized_RI_metric_block'

      INTEGER                                            :: atom_a, atom_b, first_a, first_b, &
                                                            handle, i, j, last_a, last_b, nopt_a, &
                                                            nopt_b
      REAL(KIND=dp), ALLOCATABLE                         :: contribution(:, :), &
                                                            reference_block(:, :), work(:, :)

      CALL timeset(routineN, handle)

      nopt_a = SIZE(transform_a%coefficient, 2)
      nopt_b = SIZE(transform_b%coefficient, 2)
      IF (symmetric) THEN
         CPASSERT(nopt_a == nopt_b)
      END IF
      ALLOCATE (block(nopt_a, nopt_b), source=0.0_dp)
      DO i = 1, SIZE(transform_a%reference_atom)
         atom_a = transform_a%reference_atom(i)
         first_a = transform_a%first_row(i)
         last_a = first_a + reference_RI_size(bs_env, atom_a) - 1
         DO j = 1, SIZE(transform_b%reference_atom)
            IF (symmetric .AND. j > i) CYCLE
            atom_b = transform_b%reference_atom(j)
            first_b = transform_b%first_row(j)
            last_b = first_b + reference_RI_size(bs_env, atom_b) - 1
            IF (coulomb) THEN
               CALL compute_AB_RI_metric_matrix(bs_env, atom_a, atom_b, reference_block)
            ELSE
               CALL compute_AB_RI_metric_matrix( &
                  bs_env, atom_a, atom_b, reference_block, bs_env%ri_metric)
               IF (atom_a == atom_b .AND. bs_env%regularization_RI > 0.0_dp) THEN
                  CALL regularize_reference_RI_metric_block(bs_env, atom_a, reference_block)
               END IF
            END IF
            ALLOCATE (work(SIZE(reference_block, 1), nopt_b))
            work(:, :) = MATMUL(reference_block, &
                                transform_b%coefficient(first_b:last_b, :))
            ALLOCATE (contribution(nopt_a, nopt_b))
            contribution(:, :) = MATMUL( &
                                 TRANSPOSE(transform_a%coefficient(first_a:last_a, :)), work)
            block(:, :) = block + contribution
            IF (symmetric .AND. i /= j) block(:, :) = block + TRANSPOSE(contribution)
            reference_block_count = reference_block_count + 1
            DEALLOCATE (contribution, reference_block, work)
         END DO
      END DO
      IF (.NOT. ALL(ieee_is_finite(block))) THEN
         CPABORT("AUTO_RI optimized metric block contains invalid values")
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_optimized_RI_metric_block

! **************************************************************************************************
!> \brief Builds M_opt or V_opt directly from contracted atom-pair blocks.
!> \param qs_env ...
!> \param bs_env ...
!> \param lists ...
!> \param metric ...
!> \param coulomb true for Coulomb; false for the GW fitting operator
!> \param reference_block_count ...
! **************************************************************************************************
   SUBROUTINE compute_contracted_RI_metric_matrix( &
      qs_env, bs_env, lists, metric, coulomb, reference_block_count)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(fit_group_list_type), DIMENSION(:), &
         INTENT(IN)                                      :: lists
      TYPE(cp_fm_type), INTENT(OUT)                      :: metric
      LOGICAL, INTENT(IN)                                :: coulomb
      INTEGER, INTENT(OUT)                               :: reference_block_count

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_contracted_RI_metric_matrix'

      INTEGER                                            :: handle, iatom, jatom, n, npcol, nprow, &
                                                            pair_index
      INTEGER, POINTER                                   :: col_dist(:), row_dist(:), sizes(:)
      LOGICAL                                            :: transform_a_ready
      REAL(KIND=dp), ALLOCATABLE                         :: block(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_distribution_type), POINTER             :: template
      TYPE(dbcsr_type)                                   :: blocks
      TYPE(optimized_atom_transform_type)                :: transform_a, transform_b

      CALL timeset(routineN, handle)

      ALLOCATE (sizes(bs_env%n_atom), row_dist(bs_env%n_atom), col_dist(bs_env%n_atom))
      sizes(:) = bs_env%auto_ri%optimized_size
      CALL get_qs_env(qs_env, dbcsr_dist=template)
      CALL dbcsr_distribution_get(template, nprows=nprow, npcols=npcol)
      DO iatom = 1, SIZE(sizes)
         row_dist(iatom) = MOD(iatom - 1, nprow)
         col_dist(iatom) = MOD(iatom - 1, npcol)
      END DO
      n = SUM(sizes)
      CALL create_distributed_matrix(qs_env, metric, n, n)
      CALL dbcsr_distribution_new(dist, template=template, row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(blocks, name='AUTO_RI contracted metric', dist=dist, &
                        matrix_type=dbcsr_type_symmetric, row_blk_size=sizes, col_blk_size=sizes)

      reference_block_count = 0
      pair_index = 0
      DO iatom = 1, SIZE(sizes)
         transform_a_ready = .FALSE.
         DO jatom = 1, iatom
            pair_index = pair_index + 1
            IF (MOD(pair_index - 1, bs_env%para_env%num_pe) /= bs_env%para_env%mepos) CYCLE
            IF (.NOT. transform_a_ready) THEN
               CALL build_optimized_atom_transform(bs_env, iatom, lists(iatom), transform_a)
               transform_a_ready = .TRUE.
            END IF
            IF (iatom == jatom) THEN
               CALL compute_optimized_RI_metric_block( &
                  bs_env, transform_a, transform_a, coulomb, .TRUE., block, reference_block_count)
               block(:, :) = 0.5_dp*(block + TRANSPOSE(block))
            ELSE
               CALL build_optimized_atom_transform(bs_env, jatom, lists(jatom), transform_b)
               CALL compute_optimized_RI_metric_block( &
                  bs_env, transform_a, transform_b, coulomb, .FALSE., block, reference_block_count)
            END IF
            CALL dbcsr_put_block(blocks, iatom, jatom, block)
            DEALLOCATE (block)
         END DO
      END DO
      CALL bs_env%para_env%sum(reference_block_count)
      CALL dbcsr_finalize(blocks)
      CALL copy_dbcsr_to_fm(blocks, metric)
      CALL dbcsr_release(blocks)
      CALL dbcsr_distribution_release(dist)
      DEALLOCATE (sizes, row_dist, col_dist)

      CALL timestop(handle)

   END SUBROUTINE compute_contracted_RI_metric_matrix

! **************************************************************************************************
!> \brief Computes the pseudoinverse of the optimized fitting metric M_opt.
!>
!>        The generalized problem M_opt X=M_AA X λ separates linearly
!>        dependent directions. The result is M_opt^+=X diag(1/λ) X^T.
!> \param qs_env ...
!> \param bs_env ...
!> \param ri_metric ...
! **************************************************************************************************
   SUBROUTINE compute_optimized_RI_metric_inverse(qs_env, bs_env, ri_metric)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: ri_metric

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_optimized_RI_metric_inverse'

      INTEGER                                            :: handle, n, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), ALLOCATABLE                         :: eigenvalues(:), factors(:)
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_metric, local_scaling
      TYPE(cp_fm_type)                                   :: eigenvectors, scaling, work

      CALL timeset(routineN, handle)

      n = SUM(bs_env%auto_ri%optimized_size)
      CALL create_distributed_matrix(qs_env, bs_env%auto_ri%inverse_ri_metric, n, n)
      CALL create_distributed_matrix(qs_env, eigenvectors, n, n)
      CALL create_distributed_matrix(qs_env, scaling, n, n)
      CALL create_distributed_matrix(qs_env, work, n, n)
      CALL cp_fm_get_info(ri_metric, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, &
                          local_data=local_metric)
      CALL cp_fm_get_info(scaling, local_data=local_scaling)
      local_scaling = 0.0_dp
      CALL copy_AA_metric_blocks(local_metric(:nrow_local, :ncol_local), &
                                 local_scaling(:nrow_local, :ncol_local), &
                                 row_indices(:nrow_local), col_indices(:ncol_local), &
                                 bs_env%auto_ri%optimized_size)

      ALLOCATE (eigenvalues(n))
      CALL cp_fm_geeig(ri_metric, scaling, eigenvectors, eigenvalues, work)
      CALL inverse_sqrt_metric_eigenvalues( &
         eigenvalues, bs_env%auto_ri%metric_rank_threshold, factors)
      CALL cp_fm_column_scale(eigenvectors, factors)
      CALL cp_fm_syrk("U", "N", n, 1.0_dp, eigenvectors, 1, 1, 0.0_dp, &
                      bs_env%auto_ri%inverse_ri_metric)
      CALL cp_fm_uplo_to_full(bs_env%auto_ri%inverse_ri_metric, work)

      DEALLOCATE (eigenvalues, factors)
      CALL cp_fm_release(eigenvectors)
      CALL cp_fm_release(scaling)
      CALL cp_fm_release(work)

      CALL timestop(handle)

   END SUBROUTINE compute_optimized_RI_metric_inverse

! **************************************************************************************************
!> \brief Computes M_opt and V_opt without constructing M_ref or V_ref.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_optimized_RI_metric_matrix(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_optimized_RI_metric_matrix'

      INTEGER                                            :: coulomb_block_count, &
                                                            fitting_block_count, handle, &
                                                            max_reference_recomputation, unit_nr
      REAL(KIND=dp)                                      :: average_recomputation
      TYPE(cp_fm_type)                                   :: ri_metric
      TYPE(fit_group_list_type), ALLOCATABLE             :: lists(:)

      CALL timeset(routineN, handle)

      CALL build_fit_group_lists(bs_env, lists, max_reference_recomputation)
      CALL initialize_auto_ri_integrals(bs_env)
      CALL compute_contracted_RI_metric_matrix( &
         qs_env, bs_env, lists, ri_metric, .FALSE., fitting_block_count)
      CALL compute_optimized_RI_metric_inverse(qs_env, bs_env, ri_metric)
      CALL cp_fm_release(ri_metric)
      CALL compute_contracted_RI_metric_matrix( &
         qs_env, bs_env, lists, bs_env%auto_ri%coulomb_matrix, .TRUE., coulomb_block_count)
      CPASSERT(fitting_block_count == coulomb_block_count)

      average_recomputation = REAL(fitting_block_count, dp)/ &
                              REAL(bs_env%n_atom*(bs_env%n_atom + 1)/2, dp)
      unit_nr = bs_env%unit_nr
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,T71,F10.2)') &
            'AUTO_RI| Average reference-block recomputation:', average_recomputation
         WRITE (unit_nr, '(T2,A,T72,I9)') &
            'AUTO_RI| Maximum reference-block recomputation bound:', &
            max_reference_recomputation
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_optimized_RI_metric_matrix

! **************************************************************************************************
!> \brief Computes the atom-local two-center RI metric (P|Q) and the
!>        three-center integrals (μν|P) used in Eqs. (3) and (4).
!> \param ao_basis ...
!> \param ri_basis ...
!> \param potential fitting operator used for both integral types
!> \param metric ...
!> \param products ...
! **************************************************************************************************
   SUBROUTINE compute_atom_local_RI_integrals(ao_basis, ri_basis, potential, metric, products)
      TYPE(gto_basis_set_type), INTENT(IN)               :: ao_basis, ri_basis
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: metric
      REAL(KIND=dp), INTENT(OUT)                         :: products(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_atom_local_RI_integrals'

      INTEGER                                            :: handle, iset, jset, kset, max_ao_l, &
                                                            max_l, nao, ncoa, ncob, ncoc, nri, &
                                                            nsgfa, nsgfb, nsgfc, sgfa, sgfb, sgfc
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_c, lmin_a, lmin_c, npgf_a, &
                                                            npgf_c, nsgf_set_a, nsgf_set_c
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_a, first_sgf_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: block_2c, raw_2c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_3c, raw_3c, three_center
      REAL(KIND=dp), DIMENSION(3)                        :: center
      REAL(KIND=dp), DIMENSION(:), POINTER               :: radius_a, radius_b, radius_c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pgf_radius_a, pgf_radius_c, sphi_a, &
                                                            sphi_c, zeta_a, zeta_c
      TYPE(coulomb_integral_context_type)                :: context

      CALL timeset(routineN, handle)

      CALL get_gto_basis_set(ao_basis, nsgf=nao, lmax=lmax_a, lmin=lmin_a, &
                             npgf=npgf_a, nsgf_set=nsgf_set_a, first_sgf=first_sgf_a, &
                             pgf_radius=pgf_radius_a, sphi=sphi_a, zet=zeta_a)
      CALL get_gto_basis_set(ri_basis, nsgf=nri, lmax=lmax_c, lmin=lmin_c, &
                             npgf=npgf_c, nsgf_set=nsgf_set_c, first_sgf=first_sgf_c, &
                             pgf_radius=pgf_radius_c, sphi=sphi_c, zet=zeta_c)
      CPASSERT(nao > 0 .AND. nri > 0)
      max_ao_l = MAXVAL(lmax_a)
      max_l = MAX(max_ao_l, MAXVAL(lmax_c))
      CALL init_orbital_pointers(MAX(max_l, 2*max_ao_l))
      CALL init_spherical_harmonics(max_l, -1)
      IF (potential%potential_type == do_potential_truncated) THEN
         IF (active_integral_library%coulomb2_library /= library_libint .OR. &
             active_integral_library%coulomb3_library /= library_libint) THEN
            CALL cp_abort(__LOCATION__, &
                          "AUTO_RI TRUNCATED fitting requires LIBINT for "// &
                          "COULOMB_2C and COULOMB_3C")
         END IF
         CALL coulomb_integral_init(context, max_am_2c=MAX(max_l, 2*max_ao_l), max_am_3c=max_l)
      ELSE
         CALL coulomb_integral_init(context, max_am_2c=max_l, max_am_3c=max_l)
      END IF
      center = 0.0_dp

      ALLOCATE (metric(nri, nri), three_center(nao, nao, nri))
      metric = 0.0_dp
      three_center = 0.0_dp

      DO iset = 1, ri_basis%nset
         ncoa = npgf_c(iset)*ncoset(lmax_c(iset))
         sgfa = first_sgf_c(1, iset)
         nsgfa = nsgf_set_c(iset)
         radius_a => pgf_radius_c(1:npgf_c(iset), iset)
         DO jset = 1, ri_basis%nset
            ncob = npgf_c(jset)*ncoset(lmax_c(jset))
            sgfb = first_sgf_c(1, jset)
            nsgfb = nsgf_set_c(jset)
            radius_b => pgf_radius_c(1:npgf_c(jset), jset)
            ALLOCATE (raw_2c(ncoa, ncob), block_2c(nsgfa, nsgfb))
            raw_2c = 0.0_dp
            CALL compute_coulomb_2c(context, lmin_c(iset), lmax_c(iset), &
                                    lmin_c(jset), lmax_c(jset), &
                                    npgf_c(iset), npgf_c(jset), zeta_c(1:npgf_c(iset), iset), &
                                    zeta_c(1:npgf_c(jset), jset), &
                                    radius_a, radius_b, center, center, raw_2c, potential)
            CALL ab_contract(block_2c, raw_2c, sphi_c(:, sgfa:), sphi_c(:, sgfb:), &
                             ncoa, ncob, nsgfa, nsgfb)
            metric(sgfa:sgfa + nsgfa - 1, sgfb:sgfb + nsgfb - 1) = block_2c
            DEALLOCATE (block_2c, raw_2c)
         END DO
      END DO

      DO iset = 1, ao_basis%nset
         ncoa = npgf_a(iset)*ncoset(lmax_a(iset))
         sgfa = first_sgf_a(1, iset)
         nsgfa = nsgf_set_a(iset)
         radius_a => pgf_radius_a(1:npgf_a(iset), iset)
         DO jset = 1, ao_basis%nset
            ncob = npgf_a(jset)*ncoset(lmax_a(jset))
            sgfb = first_sgf_a(1, jset)
            nsgfb = nsgf_set_a(jset)
            radius_b => pgf_radius_a(1:npgf_a(jset), jset)
            DO kset = 1, ri_basis%nset
               ncoc = npgf_c(kset)*ncoset(lmax_c(kset))
               sgfc = first_sgf_c(1, kset)
               nsgfc = nsgf_set_c(kset)
               radius_c => pgf_radius_c(1:npgf_c(kset), kset)
               ALLOCATE (raw_3c(ncoa, ncob, ncoc), block_3c(nsgfa, nsgfb, nsgfc))
               raw_3c = 0.0_dp
               CALL compute_coulomb_3c(context, lmin_a(iset), lmax_a(iset), &
                                       lmin_a(jset), lmax_a(jset), &
                                       lmin_c(kset), lmax_c(kset), npgf_a(iset), npgf_a(jset), &
                                       npgf_c(kset), zeta_a(1:npgf_a(iset), iset), &
                                       zeta_a(1:npgf_a(jset), jset), zeta_c(1:npgf_c(kset), kset), &
                                       radius_a, radius_b, radius_c, center, center, center, &
                                       raw_3c, potential)
               CALL abc_contract(block_3c, raw_3c, sphi_a(:, sgfa:), sphi_a(:, sgfb:), &
                                 sphi_c(:, sgfc:), ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc)
               three_center(sgfa:sgfa + nsgfa - 1, sgfb:sgfb + nsgfb - 1, &
                            sgfc:sgfc + nsgfc - 1) = block_3c
               DEALLOCATE (block_3c, raw_3c)
            END DO
         END DO
      END DO

      products(:, :, :) = three_center
      CALL coulomb_integral_cleanup(context)
      IF (.NOT. ALL(ieee_is_finite(metric))) THEN
         CPABORT("AUTO_RI metric contains invalid values")
      END IF
      IF (.NOT. ALL(ieee_is_finite(products))) THEN
         CPABORT("AUTO_RI three-center integrals contain invalid values")
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_atom_local_RI_integrals

! **************************************************************************************************
!> \brief Initializes the integral kernels used by the AUTO_RI optimization.
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE initialize_auto_ri_integrals(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_auto_ri_integrals'

      CHARACTER(LEN=default_path_length)                 :: table_file
      INTEGER                                            :: handle, ikind, max_l, table_lmax, unit_id

      CALL timeset(routineN, handle)

      max_l = 0
      DO ikind = 1, SIZE(bs_env%basis_set_RI)
         CPASSERT(ASSOCIATED(bs_env%basis_set_RI(ikind)%gto_basis_set))
         max_l = MAX(max_l, MAXVAL(bs_env%basis_set_RI(ikind)%gto_basis_set%lmax))
      END DO
      ! The OS two-center kernel accesses l+1 even without force evaluation.
      CALL init_orbital_pointers(max_l + 1)
      CALL init_spherical_harmonics(max_l, -1)

      IF (bs_env%ri_metric%potential_type == do_potential_truncated) THEN
         ! The truncated-Coulomb kernel needs auxiliary functions through l_A+l_B+1.
         table_lmax = 2*max_l + 1
         IF (get_lmax_init() < table_lmax) THEN
            unit_id = -1
            IF (bs_env%para_env%is_source()) THEN
               table_file = discover_file(bs_env%ri_metric%filename)
               CALL open_file(unit_number=unit_id, file_name=TRIM(table_file))
            END IF
            CALL init_t_c_g0(table_lmax, unit_id, bs_env%para_env%mepos, bs_env%para_env)
            IF (bs_env%para_env%is_source()) CALL close_file(unit_id)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE initialize_auto_ri_integrals

! **************************************************************************************************
!> \brief Computes the atom-pair two-center RI metric matrix (P_A|Q_B).
!> \param bs_env ...
!> \param atom_a ...
!> \param atom_b ...
!> \param block ...
!> \param ri_potential ...
! **************************************************************************************************
   SUBROUTINE compute_AB_RI_metric_matrix(bs_env, atom_a, atom_b, block, ri_potential)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom_a, atom_b
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: block
      TYPE(libint_potential_type), INTENT(IN), OPTIONAL  :: ri_potential

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_AB_RI_metric_matrix'

      INTEGER                                            :: atom_1, atom_2, handle, kind_1, kind_2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: canonical_block
      TYPE(libint_potential_type)                        :: potential

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(bs_env%ri_rs%particle_set))
      CPASSERT(atom_a >= 1 .AND. atom_a <= bs_env%n_atom)
      CPASSERT(atom_b >= 1 .AND. atom_b <= bs_env%n_atom)
      IF (atom_a >= atom_b) THEN
         atom_1 = atom_a
         atom_2 = atom_b
      ELSE
         atom_1 = atom_b
         atom_2 = atom_a
      END IF
      kind_1 = bs_env%ri_rs%particle_set(atom_1)%atomic_kind%kind_number
      kind_2 = bs_env%ri_rs%particle_set(atom_2)%atomic_kind%kind_number

      IF (PRESENT(ri_potential)) THEN
         CPASSERT(ri_potential%potential_type == do_potential_truncated)
         potential = ri_potential
      ELSE
         potential%potential_type = do_potential_coulomb
      END IF
      CALL compute_2c_RI_metric( &
         basis_a=bs_env%basis_set_RI(kind_1)%gto_basis_set, &
         basis_b=bs_env%basis_set_RI(kind_2)%gto_basis_set, &
         center_a=bs_env%ri_rs%particle_set(atom_1)%r, &
         center_b=bs_env%ri_rs%particle_set(atom_2)%r, &
         potential=potential, block=canonical_block)

      IF (atom_a >= atom_b) THEN
         CALL MOVE_ALLOC(canonical_block, block)
      ELSE
         ALLOCATE (block(SIZE(canonical_block, 2), SIZE(canonical_block, 1)))
         block(:, :) = TRANSPOSE(canonical_block)
         DEALLOCATE (canonical_block)
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_AB_RI_metric_matrix

! **************************************************************************************************
!> \brief Computes one contracted two-center RI metric block (P_A|Q_B).
!> \param basis_a ...
!> \param basis_b ...
!> \param center_a ...
!> \param center_b ...
!> \param potential ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE compute_2c_RI_metric(basis_a, basis_b, center_a, center_b, potential, block)
      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: center_a, center_b
      TYPE(libint_potential_type), INTENT(IN)            :: potential
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: block

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_2c_RI_metric'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (block(basis_a%nsgf, basis_b%nsgf), source=0.0_dp)
      SELECT CASE (potential%potential_type)
      CASE (do_potential_truncated)
         CALL int_operators_r12_ab_os(operator_truncated, block, rab=center_b - center_a, &
                                      fba=basis_a, fbb=basis_b, &
                                      r_cutoff=potential%cutoff_radius, calculate_forces=.FALSE.)
      CASE (do_potential_coulomb)
         CALL int_operators_r12_ab_os(operator_coulomb, block, rab=center_b - center_a, &
                                      fba=basis_a, fbb=basis_b, calculate_forces=.FALSE.)
      CASE DEFAULT
         CPABORT("AUTO_RI two-center metric requires COULOMB or TRUNCATED")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE compute_2c_RI_metric

! **************************************************************************************************
!> \brief Returns the reference RI basis size of one atom.
!> \param bs_env ...
!> \param iatom ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION reference_RI_size(bs_env, iatom) RESULT(n)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom

      INTEGER                                            :: ikind

      ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
      n = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf

   END FUNCTION reference_RI_size

END MODULE gw_auto_ri_generate_RI_basis
