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

   ! Numerical cutoffs for the reference M_PQ eigenspace and final M_pq pseudoinverse.
   REAL(KIND=dp), PARAMETER :: reference_M_eigenval_threshold = 1.0E-10_dp
   REAL(KIND=dp), PARAMETER :: optimized_M_eigenval_threshold = 1.0E-6_dp
   ! Discard null target directions and reject non-positive-semidefinite M_PQ or B_PQ spectra.
   REAL(KIND=dp), PARAMETER :: relative_target_eigenval_threshold = 1.0E-12_dp
   REAL(KIND=dp), PARAMETER :: AB_absolute_target_eigenval_threshold = 1.0E-12_dp
   REAL(KIND=dp), PARAMETER :: negative_eigenval_tolerance = 1.0E-10_dp

   PUBLIC :: generate_auto_ri_basis

   TYPE AB_block_list_type
      INTEGER, ALLOCATABLE                               :: block(:)
      INTEGER, ALLOCATABLE                               :: first_p_AB(:)
      INTEGER, ALLOCATABLE                               :: first_p_atom(:)
      INTEGER, ALLOCATABLE                               :: size_opt_RI(:)
   END TYPE AB_block_list_type

   TYPE optimized_atom_U_type
      INTEGER, ALLOCATABLE                               :: ref_atom(:)
      INTEGER, ALLOCATABLE                               :: first_P(:)
      REAL(KIND=dp), ALLOCATABLE                         :: U_Pp(:, :)
   END TYPE optimized_atom_U_type

CONTAINS

!**************************************************************************************************!
! The full derivation of the AUTO_RI algorithm will be published separately.                       !
!                                                                                                  !
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
! The atom-local correlation target B_PQ^{A,C} mimics contributions to the density response χ,     !
! and the exchange target B_PQ^{A,X} mimics contributions to the exchange self-energy:             !
!                                                                                                  !
! B_PQ^{A,C} = (1/N_C^A) Σ_ia^oe w_ia (φ_P^A|ρ_ia^AA)(ρ_ia^AA|φ_Q^A).                      (3)     !
!                                                                                                  !
! B_PQ^{A,X} = (1/N_X^A) Σ_{n∈W_occ} Σ_i^occ (φ_P^A|ρ_ni^AA)(ρ_ni^AA|φ_Q^A).                 (4)     !
!                                                                                                  !
! W_occ = [ε_HOMO-E_win, ε_LUMO+E_win], where E_win is OCC_EMPTY_FRONTIER_ORBITAL_WINDOW.          !
!                                                                                                  !
! Here w_ia gives larger weight to frontier-orbital products close to the HOMO and LUMO:           !
!                                                                                                  !
! w_ia = 1/(|ε_i-E_F| |ε_a-E_F|),   E_F = (ε_HOMO+ε_LUMO)/2.                               (5)     !
!                                                                                                  !
! Diagonalizing B^A = B^{A,C}+B^{A,X} identifies the most important atom-local RI functions        !
! that we use in the optimized RI basis set. The contraction coefficients U_Pp^A [see Eq. (7)]     !
! used for optimized RI basis functions are obtained from the eigenvalue problem                   !
!                                                                                                  !
! B^A U^A = M^A U^A λ^A,   Σ_Q B_PQ^A U_Qp^A = λ_p^A Σ_Q M_PQ^A U_Qp^A.                    (6)     !
!                                                                                                  !
! The eigenvectors with the largest eigenvalues are the contracted RI functions φ_p^A(r)           !
! used for the optimized RI basis set:                                                             !
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
! B_PQ^{AB,X} = (1/N_X^{AB}) Σ_{n∈W_occ} Σ_i^occ (φ_P^{AB,⊥}|ρ̅_ni^AB)(ρ̅_ni^AB|φ_Q^{AB,⊥}).  (11)     !
!                                                                                                  !
! Diagonalizing their sum gives the AB-neighbor-pair contraction coefficients:                     !
!                                                                                                  !
! B^AB U^AB = U^AB λ^AB,   B^AB = B^{AB,C}+B^{AB,X}.                                      (12)     !
!                                                                                                  !
! The eigenvectors with the largest eigenvalues are used as AB contractions in the optimized       !
! RI basis set:                                                                                    !
!                                                                                                  !
! φ_p^AB(r) = Σ_P U_Pp^AB φ_P^{AB,⊥}(r).                                                  (13)     !
!                                                                                                  !
! Eq. (13) means that the RI basis is not an atom-centered basis any more, but can be delocalized  !
! over two neighboring atoms AB.                                                                   !
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
      REAL(KIND=dp), ALLOCATABLE                         :: U_Pp_A(:, :, :)

      CALL timeset(routineN, handle)

      time_start = m_walltime()

      ! Compute AA target matrices B^{AA,C}, B^{AA,X} and contraction U_Pp^A in Eqs. (3)-(7).
      CALL compute_AA_optimized_RI_basis(bs_env, U_Pp_A)

      ! Compute AB target matrices B^{AB,C}, B^{AB,X} and contraction U_Pp^AB in Eqs. (8)-(13).
      CALL compute_AB_optimized_RI_basis(qs_env, bs_env, U_Pp_A)
      DEALLOCATE (U_Pp_A)

      ! Compute and invert M_pq = (φ_p|φ_q) of the optimized RI basis in Eqs. (7) and (13).
      CALL compute_M_pq_inv_and_V_pq(qs_env, bs_env)

      ! Print the optimized RI basis dimension and total optimization time.
      CALL print_auto_ri_summary(bs_env, time_start)

      CALL timestop(handle)

   END SUBROUTINE generate_auto_ri_basis

! **************************************************************************************************
!> \brief Computes the atom-local contraction coefficients U^A from Eqs. (3)-(7):
!>
!>          B_PQ^{A,C} = (1/N_C^A) Σ_ia^oe w_ia (φ_P^A|ρ_ia^AA)(ρ_ia^AA|φ_Q^A),         (3)
!>
!>          B_PQ^{A,X} = (1/N_X^A) Σ_{n∈W_occ} Σ_i^occ (φ_P^A|ρ_ni^AA)(ρ_ni^AA|φ_Q^A),    (4)
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
!> \param U_Pp_A temporary atom-local contraction coefficients used to construct the AB space
! **************************************************************************************************
   SUBROUTINE compute_AA_optimized_RI_basis(bs_env, U_Pp_A)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: U_Pp_A(:, :, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_AA_optimized_RI_basis'

      INTEGER :: handle, iatom, ikind, max_nAO_kind, max_nRI_kind, nAO_kind, natom, nkind, &
         nRI_kind, RI_size_opt_AA_requested, root_rank
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of_atom, sizes_opt_RI_AA
      REAL(KIND=dp)                                      :: target_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: B_PQ, eigenval_AA, M_PQ, U_Pp_A_atom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: M_PQ_by_kind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: Int_3c
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      ! Determine reference dimensions and the number of AA RI functions (i.e. number of
      ! eigenvectors U_Pp^A used for RI-basis contractions).
      CALL initialize_AA_optimization(bs_env, para_env, kind_of_atom, natom, nkind, &
                                      max_nRI_kind, RI_size_opt_AA_requested)
      ALLOCATE (eigenval_AA(max_nRI_kind, natom), source=0.0_dp)
      ALLOCATE (sizes_opt_RI_AA(natom))

      CALL initialize_auto_ri_integrals(bs_env)
      max_nAO_kind = MAXVAL(bs_env%sizes_AO)
      ALLOCATE (U_Pp_A(max_nRI_kind, max_nRI_kind, natom), &
                M_PQ_by_kind(max_nRI_kind, max_nRI_kind, nkind), source=0.0_dp)
      ALLOCATE (Int_3c(max_nAO_kind, max_nAO_kind, max_nRI_kind, nkind), source=0.0_dp)

      ! The few atom-kind integral blocks are inexpensive; compute them on every rank.
      DO ikind = 1, nkind
         nAO_kind = bs_env%basis_set_AO(ikind)%gto_basis_set%nsgf
         nRI_kind = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
         CALL compute_int_3c( &
            ao_basis=bs_env%basis_set_AO(ikind)%gto_basis_set, &
            ri_basis=bs_env%basis_set_RI(ikind)%gto_basis_set, &
            potential=bs_env%ri_metric, M_PQ=M_PQ, &
            Int_3c=Int_3c(1:nAO_kind, 1:nAO_kind, 1:nRI_kind, ikind))
         CALL add_RI_metric_regularization(bs_env%basis_set_RI(ikind)%gto_basis_set, &
                                           bs_env%regularization_RI, M_PQ)
         M_PQ_by_kind(1:nRI_kind, 1:nRI_kind, ikind) = M_PQ
         DEALLOCATE (M_PQ)
      END DO

      bs_env%auto_ri%sizes_opt_RI = 0
      DO iatom = 1, natom
         root_rank = MODULO(iatom - 1, para_env%num_pe)
         ikind = kind_of_atom(iatom)
         nRI_kind = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
         nAO_kind = bs_env%sizes_AO(iatom)
         ALLOCATE (M_PQ(nRI_kind, nRI_kind))
         M_PQ(:, :) = M_PQ_by_kind(1:nRI_kind, 1:nRI_kind, ikind)
         ! Eqs. (3)-(5): form B^{A,C}, B^{A,X}, and their weighted sum B^A.
         CALL compute_B_PQ_AA(bs_env, iatom, root_rank, &
                              Int_3c(1:nAO_kind, 1:nAO_kind, 1:nRI_kind, ikind), &
                              M_PQ, B_PQ, target_norm)
         DEALLOCATE (M_PQ)
         IF (root_rank /= para_env%mepos) THEN
            DEALLOCATE (B_PQ)
            CYCLE
         END IF
         ! Eq. (6): solve B^A U^A = M^A U^A λ^A.
         CALL solve_AA_eigenproblem(M_PQ_by_kind(1:nRI_kind, 1:nRI_kind, ikind), &
                                    B_PQ, target_norm, U_Pp_A_atom, eigenval)
         eigenval_AA(1:SIZE(U_Pp_A_atom, 2), iatom) = eigenval(1:SIZE(U_Pp_A_atom, 2))
         bs_env%auto_ri%sizes_opt_RI(iatom) = SIZE(U_Pp_A_atom, 2)
         U_Pp_A(1:SIZE(U_Pp_A_atom, 1), 1:SIZE(U_Pp_A_atom, 2), iatom) = U_Pp_A_atom
         DEALLOCATE (U_Pp_A_atom, eigenval, B_PQ)
      END DO
      CALL para_env%sum(bs_env%auto_ri%sizes_opt_RI)
      CALL para_env%sum(U_Pp_A)
      CALL para_env%sum(eigenval_AA)
      CALL select_RI_sizes(eigenval_AA, RI_size_opt_AA_requested, 1, sizes_opt_RI_AA)
      bs_env%auto_ri%sizes_opt_RI(:) = sizes_opt_RI_AA
      DO iatom = 1, natom
         U_Pp_A(:, sizes_opt_RI_AA(iatom) + 1:, iatom) = 0.0_dp
      END DO
      DEALLOCATE (Int_3c, M_PQ_by_kind)

      CALL timestop(handle)

   END SUBROUTINE compute_AA_optimized_RI_basis

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
!> \param kind_of_atom ...
!> \param natom ...
!> \param nkind ...
!> \param max_nRI_kind ...
!> \param RI_size_opt_AA_requested ...
! **************************************************************************************************
   SUBROUTINE initialize_AA_optimization(bs_env, para_env, kind_of_atom, &
                                         natom, nkind, max_nRI_kind, &
                                         RI_size_opt_AA_requested)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: kind_of_atom(:)
      INTEGER, INTENT(OUT)                               :: natom, nkind, max_nRI_kind, &
                                                            RI_size_opt_AA_requested

      CHARACTER(LEN=*), PARAMETER :: routineN = 'initialize_AA_optimization'

      INTEGER                                            :: handle, ikind, RI_size_requested
      INTEGER, ALLOCATABLE                               :: sizes_ref_RI_by_kind(:)

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
                               kind_of=kind_of_atom)
      ALLOCATE (bs_env%auto_ri%sizes_opt_RI(natom))
      ALLOCATE (sizes_ref_RI_by_kind(nkind))
      DO ikind = 1, nkind
         CPASSERT(ASSOCIATED(bs_env%basis_set_AO(ikind)%gto_basis_set))
         CPASSERT(ASSOCIATED(bs_env%basis_set_RI(ikind)%gto_basis_set))
         sizes_ref_RI_by_kind(ikind) = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf
      END DO
      max_nRI_kind = MAXVAL(sizes_ref_RI_by_kind)
      RI_size_requested = NINT( &
                          MIN(bs_env%auto_ri%ri_ao_ratio, &
                              REAL(SUM(sizes_ref_RI_by_kind(kind_of_atom)), dp)/ &
                              REAL(SUM(bs_env%sizes_AO), dp))*REAL(SUM(bs_env%sizes_AO), dp))
      RI_size_opt_AA_requested = MAX(natom, NINT(0.5_dp*REAL(RI_size_requested, dp)))
      DEALLOCATE (sizes_ref_RI_by_kind)

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

      bs_env%sizes_RI(:) = bs_env%auto_ri%sizes_opt_RI
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

      INTEGER                                            :: handle, RI_size_optimized, unit_nr
      REAL(KIND=dp)                                      :: elapsed

      CALL timeset(routineN, handle)

      RI_size_optimized = SUM(bs_env%auto_ri%sizes_opt_RI)
      elapsed = m_walltime() - time_start
      CALL bs_env%para_env%max(elapsed)
      unit_nr = bs_env%unit_nr
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(A)') ' '
         WRITE (unit_nr, '(T2,A,T74,I9)') &
            'AUTO_RI| Number of automatic RI functions for χ, ε, W:', RI_size_optimized
         WRITE (unit_nr, '(T2,A,T71,F8.1,A)') &
            'AUTO_RI| Total RI basis optimization time:', elapsed, ' s'
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_auto_ri_summary

! **************************************************************************************************
!> \brief Adds REGULARIZATION_RI to a two-center RI metric matrix (P|Q).
!> \param basis ...
!> \param regularization ...
!> \param M_PQ reference RI metric M_PQ
!> \param offset ...
! **************************************************************************************************
   SUBROUTINE add_RI_metric_regularization(basis, regularization, M_PQ, offset)
      USE basis_set_types, ONLY: gto_basis_set_type
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      REAL(KIND=dp), INTENT(IN)                          :: regularization
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: M_PQ
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
      CPASSERT(start >= 0 .AND. start + basis%nsgf <= SIZE(M_PQ, 1))
      CPASSERT(start + basis%nsgf <= SIZE(M_PQ, 2))
      DO iset = 1, nset
         shift = regularization*MAX(1.0_dp, 1.0_dp/MINVAL(zet(1:npgf(iset), iset)))
         first = start + first_sgf(1, iset)
         DO i = first, first + nsgf_set(iset) - 1
            M_PQ(i, i) = M_PQ(i, i) + shift
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE add_RI_metric_regularization

! **************************************************************************************************
!> \brief Copies atom-local diagonal blocks M^AA from a distributed metric M.
!> \param M_pq_local locally stored entries of the optimized RI metric M_pq
!> \param M_pq_AA_local locally stored entries of the atom-block-diagonal metric
!> \param row_indices global indices of the locally owned rows
!> \param col_indices global indices of the locally owned columns
!> \param sizes_RI optimized RI functions in each atomic matrix block
! **************************************************************************************************
   SUBROUTINE copy_M_pq_AA_blocks(M_pq_local, M_pq_AA_local, row_indices, &
                                  col_indices, sizes_RI)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: M_pq_local
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: M_pq_AA_local
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_indices, col_indices, sizes_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'copy_M_pq_AA_blocks'

      INTEGER                                            :: handle, i, j, last_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: block_ends, col_atom, row_atom

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(M_pq_local, 1) == SIZE(row_indices))
      CPASSERT(SIZE(M_pq_local, 2) == SIZE(col_indices))
      CPASSERT(ALL(SHAPE(M_pq_AA_local) == SHAPE(M_pq_local)))
      CPASSERT(SIZE(sizes_RI) > 0 .AND. ALL(sizes_RI > 0))

      ALLOCATE (block_ends(SIZE(sizes_RI)))
      block_ends(1) = sizes_RI(1)
      DO i = 2, SIZE(sizes_RI)
         block_ends(i) = block_ends(i - 1) + sizes_RI(i)
      END DO
      last_index = block_ends(SIZE(block_ends))
      CPASSERT(ALL(row_indices >= 1 .AND. row_indices <= last_index))
      CPASSERT(ALL(col_indices >= 1 .AND. col_indices <= last_index))

      ALLOCATE (row_atom(SIZE(row_indices)), col_atom(SIZE(col_indices)))
      DO i = 1, SIZE(row_indices)
         row_atom(i) = M_pq_index_atom(row_indices(i), block_ends)
      END DO
      DO j = 1, SIZE(col_indices)
         col_atom(j) = M_pq_index_atom(col_indices(j), block_ends)
      END DO

      M_pq_AA_local = 0.0_dp
      DO j = 1, SIZE(col_indices)
         DO i = 1, SIZE(row_indices)
            IF (row_atom(i) == col_atom(j)) THEN
               M_pq_AA_local(i, j) = M_pq_local(i, j)
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE copy_M_pq_AA_blocks

! **************************************************************************************************
!> \brief Forms f_i=1/sqrt(λ_i) for eigenvalues of the RI metric matrix M_PQ.
!> \param eigenval generalized metric eigenvalues
!> \param factors inverse-square-root weights, or zero for discarded modes
! **************************************************************************************************
   SUBROUTINE compute_optimized_M_pq_inv_factors(eigenval, factors)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: factors

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_optimized_M_pq_inv_factors'

      INTEGER                                            :: handle, i, n_dependent
      REAL(KIND=dp)                                      :: largest, threshold

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(eigenval) > 0)

      largest = MAXVAL(eigenval)
      CPASSERT(largest > 0.0_dp)
      threshold = optimized_M_eigenval_threshold*largest
      ALLOCATE (factors(SIZE(eigenval)))
      factors = 0.0_dp
      n_dependent = 0
      DO i = 1, SIZE(eigenval)
         IF (eigenval(i) > 0.0_dp .AND. eigenval(i) >= threshold) THEN
            factors(i) = 1.0_dp/SQRT(eigenval(i))
         ELSE
            n_dependent = n_dependent + 1
         END IF
      END DO
      CPASSERT(n_dependent < SIZE(eigenval))

      CALL timestop(handle)

   END SUBROUTINE compute_optimized_M_pq_inv_factors

! **************************************************************************************************
!> \brief Finds the atom owning a global RI metric index.
!> \param global_index global metric index
!> \param block_ends inclusive end index of each consecutive atom block
!> \return atom index
! **************************************************************************************************
   PURE FUNCTION M_pq_index_atom(global_index, block_ends) RESULT(atom_index)
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

   END FUNCTION M_pq_index_atom

! **************************************************************************************************
!> \brief Computes X=U diag(1/sqrt(λ)) with X^T M X=I.
!>
!>        The columns of X span the non-null eigenspace of the symmetric
!>        two-center RI metric matrix M=U diag(λ) U^T.
!> \param M_PQ ...
!> \param X_Pi ...
!> \param rank ...
! **************************************************************************************************
   SUBROUTINE compute_reference_M_PQ_inv_factor(M_PQ, X_Pi, rank)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: M_PQ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: X_Pi
      INTEGER, INTENT(OUT)                               :: rank

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_reference_M_PQ_inv_factor'

      INTEGER                                            :: column, handle, i, n
      REAL(KIND=dp)                                      :: largest
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvec

      CALL timeset(routineN, handle)

      n = SIZE(M_PQ, 1)
      CPASSERT(n > 0 .AND. SIZE(M_PQ, 2) == n)
      ALLOCATE (eigenvec, source=0.5_dp*(M_PQ + TRANSPOSE(M_PQ)))
      ALLOCATE (eigenval(n))
      CALL diamat_all(eigenvec, eigenval, dac=.TRUE.)
      largest = MAXVAL(eigenval)
      IF (largest <= 0.0_dp) CPABORT("AUTO_RI metric has no positive non-null range")
      IF (MINVAL(eigenval) < -negative_eigenval_tolerance*largest) THEN
         CPABORT("AUTO_RI metric has no positive non-null range")
      END IF
      rank = COUNT(eigenval > reference_M_eigenval_threshold*largest)
      IF (rank == 0) CPABORT("AUTO_RI metric has no positive non-null range")
      ALLOCATE (X_Pi(n, rank))
      column = 0
      DO i = 1, n
         IF (eigenval(i) <= reference_M_eigenval_threshold*largest) CYCLE
         column = column + 1
         X_Pi(:, column) = eigenvec(:, i)/SQRT(eigenval(i))
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_reference_M_PQ_inv_factor

! **************************************************************************************************
!> \brief Solves eigenvalue problem B^A U^A=M^A U^A λ^A in Eq. (6).
!>
!>        Columns are ordered from the largest to the smallest λ^A.
!>        Columns are retained when λ_i/λ_max >= eigenval_threshold.
!> \param M_PQ metric of the large reference RI basis
!> \param B_PQ normalized atom-local target matrix B^A
!> \param target_normalization ...
!> \param U_Pp_A selected AA contractions from the large RI basis P to optimized functions p
!> \param eigenval normalized target eigenvalues in descending order
! **************************************************************************************************
   SUBROUTINE solve_AA_eigenproblem(M_PQ, B_PQ, target_normalization, U_Pp_A, eigenval)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: M_PQ, B_PQ
      REAL(KIND=dp), INTENT(IN)                          :: target_normalization
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: U_Pp_A
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: eigenval

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_AA_eigenproblem'

      INTEGER                                            :: handle, i, j, M_PQ_rank, n, nret
      REAL(KIND=dp)                                      :: largest_score, relative_fitting_defect
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval_B, eigenval_M
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: B_ij, M_pp_check, U_ij, U_Pp_A_rotated, &
                                                            X_Pi

      CALL timeset(routineN, handle)

      n = SIZE(M_PQ, 1)
      CPASSERT(n > 0)
      CPASSERT(SIZE(M_PQ, 2) == n)
      CPASSERT(ALL(SHAPE(B_PQ) == [n, n]))
      CPASSERT(target_normalization > 0.0_dp)
      CALL compute_reference_M_PQ_inv_factor(M_PQ, X_Pi, M_PQ_rank)

      ALLOCATE (B_ij(M_PQ_rank, M_PQ_rank))
      B_ij(:, :) = MATMUL(TRANSPOSE(X_Pi), &
                          MATMUL(0.5_dp*(B_PQ + TRANSPOSE(B_PQ)), X_Pi))/target_normalization
      ALLOCATE (U_ij, source=B_ij)
      ALLOCATE (eigenval_B(M_PQ_rank))
      CALL diamat_all(U_ij, eigenval_B, dac=.TRUE.)
      largest_score = MAXVAL(eigenval_B)
      IF (largest_score <= 0.0_dp) THEN
         CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
      END IF
      IF (MINVAL(eigenval_B) < -negative_eigenval_tolerance*largest_score) THEN
         CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
      END IF
      nret = COUNT(eigenval_B/largest_score >= relative_target_eigenval_threshold)
      ALLOCATE (U_Pp_A(n, nret), eigenval(M_PQ_rank))
      DO i = 1, M_PQ_rank
         eigenval(i) = eigenval_B(M_PQ_rank - i + 1)
      END DO
      DO i = 1, nret
         j = M_PQ_rank - i + 1
         U_Pp_A(:, i) = MATMUL(X_Pi, U_ij(:, j))
      END DO

      relative_fitting_defect = 1.0_dp
      DO i = 1, nret
         j = M_PQ_rank - i + 1
         relative_fitting_defect = relative_fitting_defect - &
                                   DOT_PRODUCT(U_ij(:, j), MATMUL(B_ij, U_ij(:, j)))
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
         ALLOCATE (M_pp_check(nret, nret), eigenval_M(nret), U_Pp_A_rotated(n, nret))
         M_pp_check(:, :) = MATMUL(TRANSPOSE(U_Pp_A), &
                                   MATMUL(0.5_dp*(M_PQ + TRANSPOSE(M_PQ)), U_Pp_A))
         CALL diamat_all(M_pp_check, eigenval_M, dac=.TRUE.)
         IF (MINVAL(eigenval_M) <= reference_M_eigenval_threshold) THEN
            CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
         END IF
         U_Pp_A_rotated(:, :) = MATMUL(U_Pp_A, M_pp_check)
         DO i = 1, nret
            U_Pp_A_rotated(:, i) = U_Pp_A_rotated(:, i)/SQRT(eigenval_M(i))
         END DO
         U_Pp_A(:, :) = MATMUL(U_Pp_A_rotated, TRANSPOSE(M_pp_check))
         M_pp_check(:, :) = MATMUL(TRANSPOSE(U_Pp_A), &
                                   MATMUL(0.5_dp*(M_PQ + TRANSPOSE(M_PQ)), U_Pp_A))
         DO i = 1, nret
            M_pp_check(i, i) = M_pp_check(i, i) - 1.0_dp
         END DO
         IF (MAXVAL(ABS(M_pp_check)) >= 1.0E-7_dp) THEN
            CPABORT("AUTO_RI conditioned atom-local eigensolver failed")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE solve_AA_eigenproblem

! **************************************************************************************************
!> \brief Solves eigenvalue problem B^AB U^AB=U^AB λ^AB of Eq. (8).
!>
!>        Columns are ordered from the largest to the smallest λ^AB.
!> \param M_PQ combined AB reference-basis metric
!> \param B_PQ normalized two-center target matrix B^AB
!> \param U_Pp_AA retained AA columns embedded in the combined AB reference space
!> \param U_Pp_AB selected AB contractions from the large RI basis P to optimized functions p
!> \param eigenval normalized scores in descending order
! **************************************************************************************************
   SUBROUTINE solve_AB_eigenproblem(M_PQ, B_PQ, U_Pp_AA, U_Pp_AB, eigenval)
      REAL(KIND=dp), INTENT(IN)                          :: M_PQ(:, :), B_PQ(:, :), U_Pp_AA(:, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: U_Pp_AB(:, :), eigenval(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_AB_eigenproblem'

      INTEGER                                            :: handle, i, k, n, nkeep, rank
      REAL(KIND=dp)                                      :: threshold
      REAL(KIND=dp), ALLOCATABLE :: B_ij(:, :), B_ij_residual(:, :), eigenval_residual(:), &
         M_AA(:, :), P_ij(:, :), Q_ia(:, :), U_Aa(:, :), U_ia(:, :), X_Pi(:, :)

      CALL timeset(routineN, handle)

      n = SIZE(M_PQ, 1)
      CALL compute_reference_M_PQ_inv_factor(M_PQ, X_Pi, rank)

      ! Xᵀ M X=I; Q spans the retained atomic columns in these orthonormal coordinates.
      ALLOCATE (Q_ia(rank, SIZE(U_Pp_AA, 2)))
      Q_ia(:, :) = MATMUL(TRANSPOSE(X_Pi), MATMUL(M_PQ, U_Pp_AA))
      ALLOCATE (P_ij(rank, rank), source=0.0_dp)
      DO i = 1, rank
         P_ij(i, i) = 1.0_dp
      END DO
      IF (SIZE(U_Pp_AA, 2) > 0) THEN
         ALLOCATE (M_AA(SIZE(U_Pp_AA, 2), SIZE(U_Pp_AA, 2)))
         M_AA(:, :) = MATMUL(TRANSPOSE(Q_ia), Q_ia)
         CALL compute_reference_M_PQ_inv_factor(M_AA, U_Aa, k)
         ALLOCATE (U_ia(rank, SIZE(U_Aa, 2)))
         U_ia(:, :) = MATMUL(Q_ia, U_Aa)
         P_ij(:, :) = P_ij - MATMUL(U_ia, TRANSPOSE(U_ia))
         DEALLOCATE (U_ia)
      END IF
      ALLOCATE (B_ij(rank, rank))
      B_ij(:, :) = MATMUL(TRANSPOSE(X_Pi), MATMUL(B_PQ, X_Pi))
      ALLOCATE (B_ij_residual(rank, rank))
      B_ij_residual(:, :) = MATMUL(P_ij, MATMUL(B_ij, P_ij))
      ALLOCATE (eigenval_residual(rank))
      CALL diamat_all(B_ij_residual, eigenval_residual, dac=.TRUE.)
      threshold = MAX(AB_absolute_target_eigenval_threshold, &
                      relative_target_eigenval_threshold*MAXVAL(eigenval_residual), &
                      100.0_dp*EPSILON(1.0_dp)*MAXVAL(ABS(B_ij)))
      nkeep = COUNT(eigenval_residual > threshold)
      ALLOCATE (U_Pp_AB(n, nkeep), eigenval(nkeep))
      k = 0
      DO i = rank, 1, -1
         IF (eigenval_residual(i) <= threshold) CYCLE
         k = k + 1
         U_Pp_AB(:, k) = MATMUL(X_Pi, MATMUL(P_ij, B_ij_residual(:, i)))
         eigenval(k) = eigenval_residual(i)
      END DO

      CALL timestop(handle)

   END SUBROUTINE solve_AB_eigenproblem

! **************************************************************************************************
!> \brief Selects the requested RI sizes from the largest positive eigenvalues.
!> \param eigenval_by_block descending eigenvalues; zero pads unavailable modes
!> \param RI_size_requested requested total number of functions
!> \param RI_size_minimum minimum retained size in every nonempty block
!> \param sizes_RI_selected number retained from each block
! **************************************************************************************************
   SUBROUTINE select_RI_sizes(eigenval_by_block, RI_size_requested, RI_size_minimum, &
                              sizes_RI_selected)
      REAL(KIND=dp), INTENT(IN)                          :: eigenval_by_block(:, :)
      INTEGER, INTENT(IN)                                :: RI_size_requested, RI_size_minimum
      INTEGER, INTENT(OUT)                               :: sizes_RI_selected(:)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'select_RI_sizes'

      INTEGER                                            :: block, handle, i, n, nmode, remaining
      INTEGER, ALLOCATABLE                               :: indices(:)
      REAL(KIND=dp), ALLOCATABLE                         :: values(:)

      CALL timeset(routineN, handle)

      nmode = SIZE(eigenval_by_block, 1)
      n = SIZE(eigenval_by_block)
      CPASSERT(SIZE(sizes_RI_selected) == SIZE(eigenval_by_block, 2))
      CPASSERT(RI_size_minimum >= 0)
      sizes_RI_selected(:) = 0
      IF (n > 0) THEN
         ALLOCATE (values(n), indices(n))
         values(:) = -RESHAPE(eigenval_by_block, [n])
         DO block = 1, SIZE(sizes_RI_selected)
            sizes_RI_selected(block) = &
               MIN(RI_size_minimum, COUNT(eigenval_by_block(:, block) > 0.0_dp))
            values((block - 1)*nmode + 1: &
                   (block - 1)*nmode + sizes_RI_selected(block)) = 0.0_dp
         END DO
         CALL sort(values, n, indices)
         remaining = MIN(MAX(0, RI_size_requested - SUM(sizes_RI_selected)), COUNT(values < 0.0_dp))
         DO i = 1, remaining
            block = (indices(i) - 1)/nmode + 1
            sizes_RI_selected(block) = sizes_RI_selected(block) + 1
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE select_RI_sizes

! **************************************************************************************************
!> \brief Computes orbital factors for the occupied-empty weight in Eq. (5).
!>
!>        weights(n) = 1/|ε_n-E_F|, so weights(i)*weights(a) = w_ia.
!> \param eigenval_MO molecular-orbital eigenvalues, occupied first
!> \param homo Last occupied state in this spin channel.
!> \param weights factors whose occupied-empty products give w_ia
! **************************************************************************************************
   SUBROUTINE compute_orbital_weights(eigenval_MO, homo, weights)
      REAL(KIND=dp), INTENT(IN)                          :: eigenval_MO(:)
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(OUT)                         :: weights(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_orbital_weights'

      INTEGER                                            :: handle, n
      REAL(KIND=dp)                                      :: gap, middle

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(weights) == SIZE(eigenval_MO))
      IF (homo < 1 .OR. homo >= SIZE(eigenval_MO)) THEN
         CPABORT("AUTO_RI orbital weights require occupied and empty states")
      END IF
      IF (.NOT. ALL(ieee_is_finite(eigenval_MO))) THEN
         CPABORT("AUTO_RI orbital weights require finite eigenvalues")
      END IF
      gap = eigenval_MO(homo + 1) - eigenval_MO(homo)
      IF (.NOT. ieee_is_finite(gap) .OR. gap <= 0.0_dp) THEN
         CPABORT("AUTO_RI orbital weights require a finite positive gap")
      END IF
      middle = 0.5_dp*(eigenval_MO(homo) + eigenval_MO(homo + 1))
      DO n = 1, SIZE(eigenval_MO)
         weights(n) = 1.0_dp/ABS(eigenval_MO(n) - middle)
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_orbital_weights

! **************************************************************************************************
!> \brief Computes the occupied-empty contribution B_PQ^C in Eqs. (3) or (10).
!>
!>        B_PQ^C = Σ_ia w_ia (P|ρ_ia)(ρ_ia|Q).
!>
!> \param Int_3c atom-local or atom-pair three-center integrals (μν|P)
!> \param C_mu_n molecular-orbital coefficients restricted to the relevant AO rows
!> \param weights Per-orbital weights for this spin channel.
!> \param homo Last occupied state.
!> \param B_PQ unnormalized transition matrix for this spin channel
! **************************************************************************************************
   SUBROUTINE compute_B_PQ_transition(Int_3c, C_mu_n, weights, homo, B_PQ)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: Int_3c(:, :, :), C_mu_n(:, :), weights(:)
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(OUT)                         :: B_PQ(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_B_PQ_transition'

      INTEGER                                            :: handle, n, nAO, nRI, nstates
      REAL(KIND=dp), ALLOCATABLE                         :: empty(:, :), occupied(:, :), scaled(:, :)
      TYPE(local_gemm_ctxt_type)                         :: gemm_context

      CALL timeset(routineN, handle)

      CALL gemm_context%create(LOCAL_GEMM_PU_HOST)
      nAO = SIZE(C_mu_n, 1)
      nstates = SIZE(C_mu_n, 2)
      nRI = SIZE(Int_3c, 3)
      CPASSERT(SIZE(Int_3c, 1) == nAO .AND. SIZE(Int_3c, 2) == nAO)
      CPASSERT(SIZE(B_PQ, 1) == nRI .AND. SIZE(B_PQ, 2) == nRI)
      CPASSERT(SIZE(weights) == nstates .AND. homo > 0 .AND. homo < nstates)
      ALLOCATE (scaled(nAO, nstates), occupied(nAO, nAO), empty(nAO, nAO))
      DO n = 1, nstates
         scaled(:, n) = C_mu_n(:, n)*SQRT(weights(n))
      END DO
      CALL gemm_context%gemm( &
         'N', 'T', nAO, nAO, homo, 1.0_dp, scaled, nAO, scaled, nAO, &
         0.0_dp, occupied, nAO)
      CALL gemm_context%gemm('N', 'T', nAO, nAO, nstates - homo, 1.0_dp, &
                             scaled(:, homo + 1:), nAO, scaled(:, homo + 1:), nAO, &
                             0.0_dp, empty, nAO)
      CALL contract_B_PQ(Int_3c, occupied, empty, B_PQ)
      CALL gemm_context%destroy()

      CALL timestop(handle)

   END SUBROUTINE compute_B_PQ_transition

! **************************************************************************************************
!> \brief Computes the frontier-exchange contribution B_PQ^X in Eqs. (4) or (11).
!>
!>        B_PQ^X = Σ_{n∈W_occ} Σ_i^occ (P|ρ_ni)(ρ_ni|Q).
!>
!> \param Int_3c atom-local or atom-pair three-center integrals (μν|P)
!> \param C_mu_n molecular-orbital coefficients restricted to the relevant AO rows
!> \param eigenval_MO molecular-orbital eigenvalues, occupied first, in Hartree
!> \param homo Last occupied state in this spin channel.
!> \param occ_energy_window frontier-window extension in Hartree
!> \param B_PQ unnormalized exchange matrix
! **************************************************************************************************
   SUBROUTINE compute_B_PQ_exchange(Int_3c, C_mu_n, eigenval_MO, homo, occ_energy_window, B_PQ)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: Int_3c(:, :, :), C_mu_n(:, :), &
                                                            eigenval_MO(:)
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(IN)                          :: occ_energy_window
      REAL(KIND=dp), INTENT(OUT)                         :: B_PQ(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_B_PQ_exchange'

      INTEGER                                            :: handle, n, nAO, nstates
      REAL(KIND=dp), ALLOCATABLE                         :: frontier(:, :), occupied(:, :), &
                                                            selected(:, :)
      TYPE(local_gemm_ctxt_type)                         :: gemm_context

      CALL timeset(routineN, handle)

      B_PQ(:, :) = 0.0_dp
      nstates = SIZE(eigenval_MO)
      IF (homo < 1 .OR. homo >= nstates) THEN
         CPABORT("AUTO_RI exchange target requires occupied and empty states")
      END IF
      IF (.NOT. ieee_is_finite(occ_energy_window) .OR. &
          .NOT. ALL(ieee_is_finite(eigenval_MO))) THEN
         CPABORT("AUTO_RI exchange target requires finite input")
      END IF
      IF (occ_energy_window < 0.0_dp .OR. eigenval_MO(homo + 1) <= eigenval_MO(homo)) THEN
         CPABORT("AUTO_RI exchange target requires a valid frontier window")
      END IF
      nAO = SIZE(C_mu_n, 1)
      CPASSERT(SIZE(C_mu_n, 2) == nstates)
      ALLOCATE (occupied(nAO, nAO), frontier(nAO, nAO), selected(nAO, nstates), source=0.0_dp)
      DO n = 1, nstates
         IF (eigenval_MO(n) < eigenval_MO(homo) - occ_energy_window .OR. &
             eigenval_MO(n) > eigenval_MO(homo + 1) + occ_energy_window) CYCLE
         selected(:, n) = C_mu_n(:, n)
      END DO
      CALL gemm_context%create(LOCAL_GEMM_PU_HOST)
      ! D_μν = Σ_i^occ C_μi C_νi.
      ! F_μν = Σ_{n: ε_n ∈ [ε_H−E_win,ε_L+E_win]} C_μn C_νn.
      CALL gemm_context%gemm('N', 'T', nAO, nAO, homo, 1.0_dp, C_mu_n, nAO, C_mu_n, nAO, &
                             0.0_dp, occupied, nAO)
      CALL gemm_context%gemm('N', 'T', nAO, nAO, nstates, 1.0_dp, selected, nAO, selected, nAO, &
                             0.0_dp, frontier, nAO)
      CALL contract_B_PQ(Int_3c, occupied, frontier, B_PQ)
      CALL gemm_context%destroy()

      CALL timestop(handle)

   END SUBROUTINE compute_B_PQ_exchange

! **************************************************************************************************
!> \brief Contracts AO density matrices and (μν|P) to form B_PQ.
!>
!>        B_PQ = Σ_μν (μν|P) [left (..|Q) right]_μν.
!>
!> \param Int_3c atom-local or atom-pair three-center integrals (μν|P)
!> \param left First density block, for occupied partners.
!> \param right Second density block, for empty or hard-window target states.
!> \param B_PQ symmetric RI-basis product matrix
! **************************************************************************************************
   SUBROUTINE contract_B_PQ(Int_3c, left, right, B_PQ)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: Int_3c(:, :, :), left(:, :), right(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: B_PQ(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_B_PQ'

      INTEGER                                            :: b, handle, mu, naux, nu, p, q
      REAL(KIND=dp), ALLOCATABLE                         :: weighted(:, :), work(:, :)
      TYPE(local_gemm_ctxt_type)                         :: gemm_context

      CALL timeset(routineN, handle)

      b = SIZE(Int_3c, 1)
      naux = SIZE(Int_3c, 3)
      CPASSERT(SIZE(Int_3c, 2) == b)
      CPASSERT(SIZE(left, 1) == b .AND. SIZE(left, 2) == b)
      CPASSERT(SIZE(right, 1) == b .AND. SIZE(right, 2) == b)
      CPASSERT(SIZE(B_PQ, 1) == naux .AND. SIZE(B_PQ, 2) == naux)
      ALLOCATE (work(b, b), weighted(b, b))
      CALL gemm_context%create(LOCAL_GEMM_PU_HOST)
      DO q = 1, naux
         ! B_PQ = Σ_μν (μν|P) [left (..|Q) right]_μν.
         CALL gemm_context%gemm('N', 'N', b, b, b, 1.0_dp, left, b, Int_3c(:, :, q), b, &
                                0.0_dp, work, b)
         CALL gemm_context%gemm('N', 'N', b, b, b, 1.0_dp, work, b, right, b, 0.0_dp, weighted, b)
         DO p = 1, q
            B_PQ(p, q) = 0.0_dp
            DO nu = 1, b
               DO mu = 1, b
                  B_PQ(p, q) = B_PQ(p, q) + Int_3c(mu, nu, p)*weighted(mu, nu)
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
!> \param X_Pi inverse metric factor, X_Pi X_Qi = M_PQ^+
!> \param B_PQ_transition transition matrix, replaced by the combined target
!> \param B_PQ_exchange unnormalized hard-window exchange covariance
!> \param target_norm number of nonzero normalized terms
! **************************************************************************************************
   SUBROUTINE combine_transition_and_exchange_targets( &
      X_Pi, B_PQ_transition, B_PQ_exchange, target_norm)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: X_Pi(:, :)
      REAL(KIND=dp), INTENT(INOUT)                       :: B_PQ_transition(:, :)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: B_PQ_exchange(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: target_norm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'combine_transition_and_exchange_targets'

      INTEGER                                            :: handle, k
      REAL(KIND=dp)                                      :: exchange_norm, transition_norm
      REAL(KIND=dp), ALLOCATABLE                         :: work(:, :)

      CALL timeset(routineN, handle)

      ALLOCATE (work(SIZE(X_Pi, 1), SIZE(X_Pi, 2)))
      work(:, :) = MATMUL(B_PQ_transition, X_Pi)
      transition_norm = 0.0_dp
      DO k = 1, SIZE(X_Pi, 2)
         transition_norm = transition_norm + DOT_PRODUCT(X_Pi(:, k), work(:, k))
      END DO
      work(:, :) = MATMUL(B_PQ_exchange, X_Pi)
      exchange_norm = 0.0_dp
      DO k = 1, SIZE(X_Pi, 2)
         exchange_norm = exchange_norm + DOT_PRODUCT(X_Pi(:, k), work(:, k))
      END DO
      target_norm = 0.0_dp
      IF (.NOT. ieee_is_finite(transition_norm) .OR. &
          .NOT. ieee_is_finite(exchange_norm) .OR. &
          transition_norm < 0.0_dp .OR. exchange_norm < 0.0_dp) THEN
         CPABORT("AUTO_RI target normalization is invalid")
      END IF
      ! B = B_ia / Tr(M^+ B_ia) + B_ni / Tr(M^+ B_ni), with zero-norm terms omitted.
      IF (transition_norm > 0.0_dp) THEN
         B_PQ_transition(:, :) = B_PQ_transition/transition_norm
         target_norm = target_norm + 1.0_dp
      ELSE
         B_PQ_transition(:, :) = 0.0_dp
      END IF
      IF (exchange_norm > 0.0_dp) THEN
         B_PQ_transition(:, :) = B_PQ_transition + B_PQ_exchange/exchange_norm
         target_norm = target_norm + 1.0_dp
      END IF
      IF (target_norm <= 0.0_dp .OR. .NOT. ALL(ieee_is_finite(B_PQ_transition))) THEN
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
!> \param Int_3c atom-local three-center integrals (μν|P)
!> \param M_PQ atom-local metric of the large reference RI basis
!> \param B_PQ normalized atom-local matrix B_PQ^A
!> \param target_norm number of nonzero normalized contributions
! **************************************************************************************************
   SUBROUTINE compute_B_PQ_AA(bs_env, iatom, root_rank, Int_3c, M_PQ, &
                              B_PQ, target_norm)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom, root_rank
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: Int_3c(:, :, :), M_PQ(:, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: B_PQ(:, :)
      REAL(KIND=dp), INTENT(OUT)                         :: target_norm

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_B_PQ_AA'

      INTEGER                                            :: column, first, handle, i, j, nAO, nRI, &
                                                            nstates, rank, row, spin
      INTEGER, POINTER                                   :: cols(:), rows(:)
      REAL(KIND=dp), ALLOCATABLE                         :: B_PQ_exchange(:, :), B_PQ_spin(:, :), &
                                                            C_mu_n(:, :), weights(:), X_Pi(:, :)

      CALL timeset(routineN, handle)
      NULLIFY (rows, cols)
      nAO = SIZE(Int_3c, 1)
      nRI = SIZE(Int_3c, 3)
      nstates = bs_env%n_mo_retained
      CPASSERT(ALLOCATED(bs_env%eigenval_scf_Gamma))
      first = 1 + SUM(bs_env%sizes_AO(:iatom - 1))
      ALLOCATE (B_PQ(nRI, nRI), source=0.0_dp)
      ALLOCATE (C_mu_n(nAO, nstates), weights(nstates), B_PQ_spin(nRI, nRI))
      target_norm = 0.0_dp
      ALLOCATE (B_PQ_exchange(nRI, nRI), source=0.0_dp)
      DO spin = 1, bs_env%n_spin
         CALL compute_orbital_weights(bs_env%eigenval_scf_Gamma(1:nstates, spin), &
                                      bs_env%n_occ(spin), weights)
         CALL cp_fm_get_info(bs_env%fm_mo_coeff_Gamma(spin), row_indices=rows, col_indices=cols)
         C_mu_n(:, :) = 0.0_dp
         DO i = 1, SIZE(rows)
            row = rows(i) - first + 1
            IF (row < 1 .OR. row > nAO) CYCLE
            DO j = 1, SIZE(cols)
               column = cols(j)
               IF (column > nstates) CYCLE
               C_mu_n(row, column) = bs_env%fm_mo_coeff_Gamma(spin)%local_data(i, j)
            END DO
         END DO
         CALL bs_env%para_env%sum(C_mu_n, root_rank)
         IF (bs_env%para_env%mepos /= root_rank) CYCLE
         CALL compute_B_PQ_transition(Int_3c, C_mu_n, weights, bs_env%n_occ(spin), B_PQ_spin)
         B_PQ(:, :) = B_PQ + B_PQ_spin
         CALL compute_B_PQ_exchange(Int_3c, C_mu_n, &
                                    bs_env%eigenval_scf_Gamma(1:nstates, spin), &
                                    bs_env%n_occ(spin), bs_env%auto_ri%occ_energy_window, B_PQ_spin)
         B_PQ_exchange(:, :) = B_PQ_exchange + B_PQ_spin
      END DO
      IF (bs_env%para_env%mepos /= root_rank) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL compute_reference_M_PQ_inv_factor(M_PQ, X_Pi, rank)
      CALL combine_transition_and_exchange_targets(X_Pi, B_PQ, B_PQ_exchange, target_norm)
      CALL timestop(handle)

   END SUBROUTINE compute_B_PQ_AA

! **************************************************************************************************
!> \brief Normalizes each contraction coefficient vector U_Pp.
!>
!>        U_Pp <- U_Pp/sqrt(Σ_P |U_Pp|^2).
!> \param U_Pp contractions from the large RI basis P to optimized functions p
! **************************************************************************************************
   SUBROUTINE normalize_RI_coefficients(U_Pp)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: U_Pp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'normalize_RI_coefficients'

      INTEGER                                            :: handle, j
      REAL(KIND=dp)                                      :: column_norm

      CALL timeset(routineN, handle)

      DO j = 1, SIZE(U_Pp, 2)
         column_norm = NORM2(U_Pp(:, j))
         IF (column_norm <= TINY(1.0_dp)) THEN
            CPABORT("AUTO_RI produced a numerically zero localized column")
         END IF
         U_Pp(:, j) = U_Pp(:, j)/column_norm
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
!> \param U_Pp_A temporary atom-local contractions used to project the AA space
! **************************************************************************************************
   SUBROUTINE compute_AB_optimized_RI_basis(qs_env, bs_env, U_Pp_A)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: U_Pp_A

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_AB_optimized_RI_basis'

      INTEGER :: AB_pair_count, center_index, current_atom, expected_RI_size, first_AB_pair, &
         handle, handle_pair_batches, iatom, jatom, k, last_AB_pair, nAO_kind, nRI_A, nRI_AB_ref, &
         nRI_B, RI_size_opt_AA, RI_size_opt_AB, RI_size_optimized, RI_size_ref, RI_size_requested, &
         root_rank, row, size_opt_RI_AB, size_U_Pp_AB, size_U_Pp_AB_storage, spin
      INTEGER, ALLOCATABLE :: AB_atom_A(:), AB_atom_B(:), AB_size_opt_RI_to_A(:), block_next(:), &
         sizes_opt_RI_AA(:), sizes_opt_RI_AB(:), sizes_ref_RI(:), U_Pp_AB_offset(:)
      REAL(KIND=dp), ALLOCATABLE :: B_PQ(:, :), C_mu_n_AB(:, :, :), C_mu_n_atom(:, :), &
         eigenval(:), eigenval_AB(:, :), M_PQ(:, :), U_Pp_AA(:, :), U_Pp_AB(:, :), &
         U_Pp_AB_packed(:)
      TYPE(gw_3c_ctx_type)                               :: ctx
      TYPE(gw_3c_ws_type)                                :: ws
      TYPE(mp_para_env_type), POINTER                    :: para

      CALL timeset(routineN, handle)

      NULLIFY (para)
      CALL get_qs_env(qs_env, para_env=para)
      ALLOCATE (sizes_ref_RI(bs_env%n_atom), sizes_opt_RI_AA(bs_env%n_atom))
      DO iatom = 1, bs_env%n_atom
         sizes_ref_RI(iatom) = get_ref_RI_size(bs_env, iatom)
      END DO
      sizes_opt_RI_AA(:) = bs_env%auto_ri%sizes_opt_RI
      RI_size_ref = SUM(sizes_ref_RI)
      RI_size_opt_AA = SUM(sizes_opt_RI_AA)

      AB_pair_count = 0
      DO iatom = 1, bs_env%n_atom - 1
         DO jatom = iatom + 1, bs_env%n_atom
            IF (NORM2(bs_env%ri_rs%particle_set(iatom)%r - bs_env%ri_rs%particle_set(jatom)%r) <= &
                bs_env%auto_ri%neighbor_radius) THEN
               AB_pair_count = AB_pair_count + 1
            END IF
         END DO
      END DO
      ALLOCATE (AB_atom_A(AB_pair_count), AB_atom_B(AB_pair_count), &
                AB_size_opt_RI_to_A(AB_pair_count), sizes_opt_RI_AB(AB_pair_count), &
                U_Pp_AB_offset(AB_pair_count))
      AB_size_opt_RI_to_A = 0
      U_Pp_AB_offset = 0
      ALLOCATE (eigenval_AB(2*MAXVAL(sizes_ref_RI), AB_pair_count), source=0.0_dp)

      CALL gw_3c_ctx_create(ctx, qs_env, bs_env%ri_metric, bs_env%basis_set_AO, &
                            bs_env%basis_set_AO, bs_env%basis_set_RI)
      CALL gw_3c_ws_create(ws, ctx)
      k = 0
      DO iatom = 1, bs_env%n_atom - 1
         DO jatom = iatom + 1, bs_env%n_atom
            IF (NORM2(bs_env%ri_rs%particle_set(iatom)%r - bs_env%ri_rs%particle_set(jatom)%r) > &
                bs_env%auto_ri%neighbor_radius) CYCLE
            k = k + 1
            AB_atom_A(k) = iatom
            AB_atom_B(k) = jatom
         END DO
      END DO
      size_U_Pp_AB_storage = 0
      DO k = 1, AB_pair_count
         root_rank = MODULO(k - 1, para%num_pe)
         IF (para%mepos /= root_rank) CYCLE
         nRI_AB_ref = sizes_ref_RI(AB_atom_A(k)) + sizes_ref_RI(AB_atom_B(k))
         U_Pp_AB_offset(k) = size_U_Pp_AB_storage + 1
         size_U_Pp_AB_storage = size_U_Pp_AB_storage + nRI_AB_ref*nRI_AB_ref
      END DO
      ALLOCATE (U_Pp_AB_packed(size_U_Pp_AB_storage))
      CALL timeset(routineN//'_pair_batches', handle_pair_batches)
      DO first_AB_pair = 1, AB_pair_count, para%num_pe
         last_AB_pair = MIN(AB_pair_count, first_AB_pair + para%num_pe - 1)
         ! Gather only one pair's MO rows per rank, then solve all batch pairs concurrently.
         DO k = first_AB_pair, last_AB_pair
            root_rank = k - first_AB_pair
            iatom = AB_atom_A(k)
            jatom = AB_atom_B(k)
            IF (para%mepos == root_rank) THEN
               ALLOCATE (C_mu_n_AB( &
                         bs_env%sizes_AO(iatom) + bs_env%sizes_AO(jatom), &
                         bs_env%n_mo_retained, bs_env%n_spin))
            END IF
            DO spin = 1, bs_env%n_spin
               row = 0
               DO center_index = 1, 2
                  current_atom = iatom
                  IF (center_index == 2) current_atom = jatom
                  nAO_kind = bs_env%sizes_AO(current_atom)
                  ALLOCATE (C_mu_n_atom(nAO_kind, bs_env%n_mo_retained))
                  CALL cp_fm_get_submatrix(bs_env%fm_mo_coeff_Gamma(spin), C_mu_n_atom, &
                                           start_row=1 + SUM(bs_env%sizes_AO(:current_atom - 1)))
                  IF (para%mepos == root_rank) THEN
                     C_mu_n_AB(row + 1:row + nAO_kind, :, spin) = C_mu_n_atom
                  END IF
                  DEALLOCATE (C_mu_n_atom)
                  row = row + nAO_kind
               END DO
            END DO
         END DO
         k = first_AB_pair + para%mepos
         IF (k <= last_AB_pair) THEN
            iatom = AB_atom_A(k)
            jatom = AB_atom_B(k)
            nRI_A = sizes_ref_RI(iatom)
            nRI_B = sizes_ref_RI(jatom)
            nRI_AB_ref = nRI_A + nRI_B
            ! Eqs. (9)-(11): form the symmetric pair density and B^{AB,C}+B^{AB,X}.
            CALL compute_B_PQ_AB(bs_env, ctx, ws, iatom, jatom, C_mu_n_AB, M_PQ, B_PQ)
            ALLOCATE (U_Pp_AA(nRI_AB_ref, sizes_opt_RI_AA(iatom) + sizes_opt_RI_AA(jatom)), &
                      source=0.0_dp)
            U_Pp_AA(1:nRI_A, 1:sizes_opt_RI_AA(iatom)) = &
               U_Pp_A(1:nRI_A, 1:sizes_opt_RI_AA(iatom), iatom)
            U_Pp_AA(nRI_A + 1:, sizes_opt_RI_AA(iatom) + 1:) = &
               U_Pp_A(1:nRI_B, 1:sizes_opt_RI_AA(jatom), jatom)
            ! Eqs. (8) and (12): project out the AA space and diagonalize the remainder.
            CALL solve_AB_eigenproblem(M_PQ, B_PQ, U_Pp_AA, U_Pp_AB, eigenval)
            eigenval_AB(1:SIZE(eigenval), k) = eigenval
            DEALLOCATE (eigenval)
            U_Pp_AB_packed(U_Pp_AB_offset(k):U_Pp_AB_offset(k) + SIZE(U_Pp_AB) - 1) = &
               RESHAPE(U_Pp_AB, [SIZE(U_Pp_AB)])
            DEALLOCATE (U_Pp_AB)
            DEALLOCATE (U_Pp_AA, B_PQ, M_PQ, C_mu_n_AB)
         END IF
      END DO
      CALL timestop(handle_pair_batches)
      k = AB_pair_count
      CALL gw_3c_ws_release(ws)
      CALL gw_3c_ctx_release(ctx)
      RI_size_requested = NINT(MIN(bs_env%auto_ri%ri_ao_ratio, &
                                   REAL(RI_size_ref, dp)/REAL(SUM(bs_env%sizes_AO), dp))* &
                               REAL(SUM(bs_env%sizes_AO), dp))
      RI_size_requested = MAX(RI_size_requested, bs_env%n_atom)

      CALL para%sum(eigenval_AB)
      CALL select_RI_sizes(eigenval_AB, RI_size_requested - SUM(sizes_opt_RI_AA), 0, &
                           sizes_opt_RI_AB)
      RI_size_opt_AB = SUM(sizes_opt_RI_AB)
      bs_env%auto_ri%sizes_opt_RI(:) = sizes_opt_RI_AA
      DO k = 1, AB_pair_count
         iatom = AB_atom_A(k)
         jatom = AB_atom_B(k)
         size_opt_RI_AB = sizes_opt_RI_AB(k)
         AB_size_opt_RI_to_A(k) = balance_AB_columns(bs_env%auto_ri%sizes_opt_RI(iatom), &
                                                     bs_env%auto_ri%sizes_opt_RI(jatom), &
                                                     size_opt_RI_AB)
         bs_env%auto_ri%sizes_opt_RI(iatom) = &
            bs_env%auto_ri%sizes_opt_RI(iatom) + AB_size_opt_RI_to_A(k)
         bs_env%auto_ri%sizes_opt_RI(jatom) = &
            bs_env%auto_ri%sizes_opt_RI(jatom) + size_opt_RI_AB - AB_size_opt_RI_to_A(k)
      END DO

      RI_size_optimized = SUM(bs_env%auto_ri%sizes_opt_RI)
      expected_RI_size = RI_size_opt_AA + RI_size_opt_AB
      CPASSERT(RI_size_optimized == expected_RI_size)
      CPASSERT(ALL(bs_env%auto_ri%sizes_opt_RI > 0))
      IF (RI_size_optimized > RI_size_ref) THEN
         CPABORT("AUTO_RI request exceeds independent reference RI directions")
      END IF
      size_U_Pp_AB = 0
      DO iatom = 1, bs_env%n_atom
         size_U_Pp_AB = size_U_Pp_AB + sizes_ref_RI(iatom)*sizes_opt_RI_AA(iatom)
      END DO
      DO k = 1, AB_pair_count
         size_U_Pp_AB = size_U_Pp_AB + &
                        (sizes_ref_RI(AB_atom_A(k)) + &
                         sizes_ref_RI(AB_atom_B(k)))*sizes_opt_RI_AB(k)
      END DO
      ALLOCATE (bs_env%auto_ri%AB_atom_A(bs_env%n_atom + AB_pair_count), source=0)
      ALLOCATE (bs_env%auto_ri%AB_atom_B(bs_env%n_atom + AB_pair_count), source=0)
      ALLOCATE (bs_env%auto_ri%AB_first_p_A(bs_env%n_atom + AB_pair_count), source=0)
      ALLOCATE (bs_env%auto_ri%AB_first_p_B(bs_env%n_atom + AB_pair_count), source=0)
      ALLOCATE (bs_env%auto_ri%AB_size_opt_RI_to_A(bs_env%n_atom + AB_pair_count), source=0)
      ALLOCATE (bs_env%auto_ri%AB_size_ref_RI(bs_env%n_atom + AB_pair_count), source=0)
      ALLOCATE (bs_env%auto_ri%AB_size_opt_RI(bs_env%n_atom + AB_pair_count), source=0)
      ALLOCATE (bs_env%auto_ri%U_Pp_AB_offset(bs_env%n_atom + AB_pair_count), source=0)
      ALLOCATE (bs_env%auto_ri%U_Pp_AB(size_U_Pp_AB), source=0.0_dp)
      bs_env%auto_ri%AB_block_count = 0

      ALLOCATE (block_next(bs_env%n_atom))
      block_next(1) = 1
      DO iatom = 2, bs_env%n_atom
         block_next(iatom) = block_next(iatom - 1) + bs_env%auto_ri%sizes_opt_RI(iatom - 1)
      END DO
      DO iatom = 1, bs_env%n_atom
         IF (sizes_opt_RI_AA(iatom) > 0) bs_env%auto_ri%AB_block_count = &
            bs_env%auto_ri%AB_block_count + 1
         CALL store_AA_columns(bs_env, iatom, &
                               U_Pp_A(1:sizes_ref_RI(iatom), 1:sizes_opt_RI_AA(iatom), iatom), &
                               block_next)
      END DO

      DO k = 1, AB_pair_count
         iatom = AB_atom_A(k)
         jatom = AB_atom_B(k)
         size_opt_RI_AB = sizes_opt_RI_AB(k)
         IF (size_opt_RI_AB == 0) CYCLE
         root_rank = MODULO(k - 1, para%num_pe)
         nRI_AB_ref = sizes_ref_RI(iatom) + sizes_ref_RI(jatom)
         ALLOCATE (U_Pp_AB(nRI_AB_ref, size_opt_RI_AB))
         IF (para%mepos == root_rank) THEN
            U_Pp_AB(:, :) = RESHAPE( &
                            U_Pp_AB_packed(U_Pp_AB_offset(k): &
                                           U_Pp_AB_offset(k) + &
                                           nRI_AB_ref*size_opt_RI_AB - 1), &
                            [nRI_AB_ref, size_opt_RI_AB])
         END IF
         CALL para%bcast(U_Pp_AB, root_rank)
         bs_env%auto_ri%AB_block_count = bs_env%auto_ri%AB_block_count + 1
         ! Eq. (13): store the selected neighbor-pair contractions in the global transform.
         CALL store_AB_columns(bs_env, U_Pp_AB, iatom, jatom, AB_size_opt_RI_to_A(k), block_next)
         DEALLOCATE (U_Pp_AB)
      END DO

      DO iatom = 1, bs_env%n_atom
         expected_RI_size = 1 + SUM(bs_env%auto_ri%sizes_opt_RI(:iatom))
         CPASSERT(block_next(iatom) == expected_RI_size)
      END DO
      DEALLOCATE (block_next, AB_atom_A, AB_atom_B, AB_size_opt_RI_to_A, &
                  U_Pp_AB_offset, U_Pp_AB_packed)

      CALL timestop(handle)

   END SUBROUTINE compute_AB_optimized_RI_basis

! **************************************************************************************************
!> \brief Divides AB columns between atoms A and B for load balance.
!> \param load_a current number of columns owned by A
!> \param load_b current number of columns owned by B
!> \param RI_size_opt_AB number of new pair RI functions
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION balance_AB_columns(load_a, load_b, RI_size_opt_AB) RESULT(RI_size_to_a)
      INTEGER, INTENT(IN)                                :: load_a, load_b, RI_size_opt_AB

      RI_size_to_a = (load_b + RI_size_opt_AB - load_a + 1)/2
      RI_size_to_a = MAX(0, MIN(RI_size_opt_AB, RI_size_to_a))

   END FUNCTION balance_AB_columns

! **************************************************************************************************
!> \brief Stores the AA contraction columns for one atom.
!> \param bs_env ...
!> \param atom ...
!> \param U_Pp_A atom-local contraction coefficients
!> \param block_next ...
! **************************************************************************************************
   SUBROUTINE store_AA_columns(bs_env, atom, U_Pp_A, block_next)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: U_Pp_A
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: block_next

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'store_AA_columns'

      INTEGER                                            :: AB_block, first, handle, RI_size_opt_AA
      REAL(KIND=dp), ALLOCATABLE                         :: U_Pp_A_normalized(:, :)

      CALL timeset(routineN, handle)

      first = block_next(atom)
      AB_block = bs_env%auto_ri%AB_block_count
      RI_size_opt_AA = SIZE(U_Pp_A, 2)
      ALLOCATE (U_Pp_A_normalized(SIZE(U_Pp_A, 1), RI_size_opt_AA))
      IF (RI_size_opt_AA > 0) THEN
         U_Pp_A_normalized(:, :) = U_Pp_A
         CALL normalize_RI_coefficients(U_Pp_A_normalized)
         bs_env%auto_ri%AB_atom_A(AB_block) = atom
         bs_env%auto_ri%AB_atom_B(AB_block) = atom
         bs_env%auto_ri%AB_first_p_A(AB_block) = &
            first - SUM(bs_env%auto_ri%sizes_opt_RI(:atom - 1))
         bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block) = RI_size_opt_AA
         CALL store_U_Pp_AB(bs_env%auto_ri, AB_block, U_Pp_A_normalized)
      END IF
      DEALLOCATE (U_Pp_A_normalized)
      block_next(atom) = first + RI_size_opt_AA

      CALL timestop(handle)

   END SUBROUTINE store_AA_columns

! **************************************************************************************************
!> \brief Stores an AB contraction in the blocks assigned to atoms A and B.
!> \param bs_env ...
!> \param U_Pp_AB neighbor-pair contractions
!> \param atom_A ...
!> \param atom_B ...
!> \param RI_size_to_a number of optimized pair RI functions assigned to atom A
!> \param block_next ...
! **************************************************************************************************
   SUBROUTINE store_AB_columns(bs_env, U_Pp_AB, atom_A, atom_B, RI_size_to_a, block_next)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: U_Pp_AB
      INTEGER, INTENT(IN)                                :: atom_A, atom_B, RI_size_to_a
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: block_next

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'store_AB_columns'

      INTEGER                                            :: AB_block, first_a, first_b, handle, &
                                                            RI_size_opt_AB

      CALL timeset(routineN, handle)

      RI_size_opt_AB = SIZE(U_Pp_AB, 2)
      AB_block = bs_env%auto_ri%AB_block_count
      first_a = block_next(atom_A)
      first_b = block_next(atom_B)
      CALL normalize_RI_coefficients(U_Pp_AB)

      bs_env%auto_ri%AB_atom_A(AB_block) = atom_A
      bs_env%auto_ri%AB_atom_B(AB_block) = atom_B
      bs_env%auto_ri%AB_first_p_A(AB_block) = first_a - SUM(bs_env%auto_ri%sizes_opt_RI(:atom_A - 1))
      bs_env%auto_ri%AB_first_p_B(AB_block) = first_b - SUM(bs_env%auto_ri%sizes_opt_RI(:atom_B - 1))
      bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block) = RI_size_to_a
      CALL store_U_Pp_AB(bs_env%auto_ri, AB_block, U_Pp_AB)
      block_next(atom_A) = first_a + RI_size_to_a
      block_next(atom_B) = first_b + RI_size_opt_AB - RI_size_to_a

      CALL timestop(handle)

   END SUBROUTINE store_AB_columns

! **************************************************************************************************
!> \brief Packs one AB U_Pp matrix into the persistent AUTO_RI data; B=A denotes an AA block.
!> \param auto_ri ...
!> \param AB_block ...
!> \param U_Pp ...
! **************************************************************************************************
   SUBROUTINE store_U_Pp_AB(auto_ri, AB_block, U_Pp)
      TYPE(auto_ri_type), INTENT(INOUT)                  :: auto_ri
      INTEGER, INTENT(IN)                                :: AB_block
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: U_Pp

      INTEGER                                            :: first, last

      first = 1
      IF (AB_block > 1) THEN
         first = auto_ri%U_Pp_AB_offset(AB_block - 1) + &
                 auto_ri%AB_size_ref_RI(AB_block - 1)* &
                 auto_ri%AB_size_opt_RI(AB_block - 1)
      END IF
      last = first + SIZE(U_Pp) - 1
      CPASSERT(last <= SIZE(auto_ri%U_Pp_AB))
      auto_ri%U_Pp_AB_offset(AB_block) = first
      auto_ri%AB_size_ref_RI(AB_block) = SIZE(U_Pp, 1)
      auto_ri%AB_size_opt_RI(AB_block) = SIZE(U_Pp, 2)
      auto_ri%U_Pp_AB(first:last) = RESHAPE(U_Pp, [SIZE(U_Pp)])

   END SUBROUTINE store_U_Pp_AB

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
!> \param C_mu_n_AB molecular-orbital coefficients on the AO rows of atoms A and B
!> \param M_PQ combined metric of the large reference RI bases on A and B
!> \param B_PQ normalized two-center matrix B_PQ^AB
! **************************************************************************************************
   SUBROUTINE compute_B_PQ_AB(bs_env, ctx, ws, iatom, jatom, C_mu_n_AB, M_PQ, B_PQ)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      TYPE(gw_3c_ws_type), INTENT(INOUT)                 :: ws
      INTEGER, INTENT(IN)                                :: iatom, jatom
      REAL(KIND=dp), INTENT(IN)                          :: C_mu_n_AB(:, :, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: M_PQ(:, :), B_PQ(:, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_B_PQ_AB'

      INTEGER :: atoms(2), handle, homo, i, j, k, kind_a, kind_b, nAO_AB, nRI_AB_ref, &
         offset_AO(2), offset_ref_RI(2), rank, sizes_AO(2), sizes_ref_RI(2), spin
      LOGICAL                                            :: screened
      REAL(KIND=dp)                                      :: norm_exchange, norm_transition
      REAL(KIND=dp), ALLOCATABLE :: B_Pi(:, :), B_PQ_exchange_AB(:, :), B_PQ_exchange_full(:, :), &
         B_PQ_spin(:, :), B_PQ_transition_AB(:, :), B_PQ_transition_full(:, :), &
         Int_3c_AB(:, :, :), Int_3c_all(:, :, :), M_PQ_block(:, :), weights(:), X_Pi(:, :)

      CALL timeset(routineN, handle)

      atoms(1) = iatom
      atoms(2) = jatom
      sizes_AO = bs_env%sizes_AO(atoms)
      DO i = 1, 2
         sizes_ref_RI(i) = get_ref_RI_size(bs_env, atoms(i))
      END DO
      offset_AO(1) = 0
      offset_AO(2) = sizes_AO(1)
      offset_ref_RI(1) = 0
      offset_ref_RI(2) = sizes_ref_RI(1)
      nAO_AB = SUM(sizes_AO)
      nRI_AB_ref = SUM(sizes_ref_RI)
      ALLOCATE (B_PQ(nRI_AB_ref, nRI_AB_ref), &
                M_PQ(nRI_AB_ref, nRI_AB_ref), source=0.0_dp)
      ALLOCATE (Int_3c_all(nAO_AB, nAO_AB, nRI_AB_ref), source=0.0_dp)
      DO i = 1, 2
         DO j = 1, 2
            CALL compute_M_PQ_AB( &
               bs_env, atoms(i), atoms(j), M_PQ_block, bs_env%ri_metric)
            M_PQ(offset_ref_RI(i) + 1:offset_ref_RI(i) + sizes_ref_RI(i), &
                 offset_ref_RI(j) + 1:offset_ref_RI(j) + sizes_ref_RI(j)) = M_PQ_block
            DEALLOCATE (M_PQ_block)
            DO k = 1, 2
               CALL build_3c_integral_block_ctx( &
                  Int_3c_all, ctx, ws, atom_j=atoms(j), atom_k=atoms(k), &
                  atom_i=atoms(i), j_offset=offset_AO(j), k_offset=offset_AO(k), &
                  i_offset=offset_ref_RI(i), screened=screened)
            END DO
         END DO
      END DO
      kind_a = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
      kind_b = bs_env%ri_rs%particle_set(jatom)%atomic_kind%kind_number
      CALL add_RI_metric_regularization(bs_env%basis_set_RI(kind_a)%gto_basis_set, &
                                        bs_env%regularization_RI, M_PQ)
      CALL add_RI_metric_regularization(bs_env%basis_set_RI(kind_b)%gto_basis_set, &
                                        bs_env%regularization_RI, M_PQ, &
                                        offset=sizes_ref_RI(1))
      ALLOCATE (Int_3c_AB, source=Int_3c_all)
      Int_3c_AB(1:sizes_AO(1), 1:sizes_AO(1), :) = 0.0_dp
      Int_3c_AB(sizes_AO(1) + 1:, sizes_AO(1) + 1:, :) = 0.0_dp
      ALLOCATE (weights(bs_env%n_mo_retained), B_PQ_spin(nRI_AB_ref, nRI_AB_ref))
      ALLOCATE (B_PQ_transition_AB(nRI_AB_ref, nRI_AB_ref), &
                B_PQ_exchange_AB(nRI_AB_ref, nRI_AB_ref), &
                B_PQ_transition_full(nRI_AB_ref, nRI_AB_ref), &
                B_PQ_exchange_full(nRI_AB_ref, nRI_AB_ref), &
                source=0.0_dp)
      DO spin = 1, bs_env%n_spin
         homo = bs_env%n_occ(spin)
         CALL compute_orbital_weights( &
            eigenval_MO=bs_env%eigenval_scf_Gamma(1:bs_env%n_mo_retained, spin), &
            homo=homo, weights=weights)
         CALL compute_B_PQ_transition(Int_3c_all, C_mu_n_AB(:, :, spin), weights, homo, B_PQ_spin)
         B_PQ_transition_full(:, :) = B_PQ_transition_full + B_PQ_spin
         CALL compute_B_PQ_transition(Int_3c_AB, C_mu_n_AB(:, :, spin), weights, homo, B_PQ_spin)
         B_PQ_transition_AB(:, :) = B_PQ_transition_AB + B_PQ_spin
         CALL compute_B_PQ_exchange(Int_3c_all, C_mu_n_AB(:, :, spin), &
                                    bs_env%eigenval_scf_Gamma(1:bs_env%n_mo_retained, spin), &
                                    homo, bs_env%auto_ri%occ_energy_window, B_PQ_spin)
         B_PQ_exchange_full(:, :) = B_PQ_exchange_full + B_PQ_spin
         CALL compute_B_PQ_exchange(Int_3c_AB, C_mu_n_AB(:, :, spin), &
                                    bs_env%eigenval_scf_Gamma(1:bs_env%n_mo_retained, spin), &
                                    homo, bs_env%auto_ri%occ_energy_window, B_PQ_spin)
         B_PQ_exchange_AB(:, :) = B_PQ_exchange_AB + B_PQ_spin
      END DO
      CALL compute_reference_M_PQ_inv_factor(M_PQ, X_Pi, rank)
      DEALLOCATE (B_PQ_spin)
      ALLOCATE (B_Pi(nRI_AB_ref, rank))
      B_Pi(:, :) = MATMUL(B_PQ_transition_full, X_Pi)
      norm_transition = SUM(X_Pi*B_Pi)
      B_Pi(:, :) = MATMUL(B_PQ_exchange_full, X_Pi)
      norm_exchange = SUM(X_Pi*B_Pi)
      ! Eqs. (10)-(12): B_PQ^AB = (B_PQ^{AB,C}/N_C^{AB}+B_PQ^{AB,X}/N_X^{AB})/2.
      IF (norm_transition > 0.0_dp) THEN
         B_PQ(:, :) = B_PQ + 0.5_dp*B_PQ_transition_AB/norm_transition
      END IF
      IF (norm_exchange > 0.0_dp) THEN
         B_PQ(:, :) = B_PQ + 0.5_dp*B_PQ_exchange_AB/norm_exchange
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_B_PQ_AB

! **************************************************************************************************
!> \brief Builds lists of AB blocks whose columns belong to each atom; B=A denotes an AA block.
!> \param bs_env ...
!> \param lists ...
! **************************************************************************************************
   SUBROUTINE build_AB_block_lists(bs_env, lists)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(AB_block_list_type), ALLOCATABLE, INTENT(OUT) :: lists(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_AB_block_lists'

      INTEGER                                            :: AB_block, atom_A, atom_B, handle, &
                                                            position, RI_size_to_a, size_opt_RI_AB
      INTEGER, ALLOCATABLE                               :: n_AB_blocks_by_atom(:)

      CALL timeset(routineN, handle)

      ALLOCATE (lists(bs_env%n_atom))
      ALLOCATE (n_AB_blocks_by_atom(bs_env%n_atom), source=0)
      DO AB_block = 1, bs_env%auto_ri%AB_block_count
         atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
         atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
         size_opt_RI_AB = bs_env%auto_ri%AB_size_opt_RI(AB_block)
         RI_size_to_a = bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
         IF (RI_size_to_a > 0) n_AB_blocks_by_atom(atom_A) = n_AB_blocks_by_atom(atom_A) + 1
         IF (atom_B /= atom_A .AND. RI_size_to_a < size_opt_RI_AB) THEN
            n_AB_blocks_by_atom(atom_B) = n_AB_blocks_by_atom(atom_B) + 1
         END IF
      END DO

      DO atom_A = 1, bs_env%n_atom
         ALLOCATE (lists(atom_A)%block(n_AB_blocks_by_atom(atom_A)))
         ALLOCATE (lists(atom_A)%first_p_AB(n_AB_blocks_by_atom(atom_A)))
         ALLOCATE (lists(atom_A)%first_p_atom(n_AB_blocks_by_atom(atom_A)))
         ALLOCATE (lists(atom_A)%size_opt_RI(n_AB_blocks_by_atom(atom_A)))
      END DO
      n_AB_blocks_by_atom = 0
      DO AB_block = 1, bs_env%auto_ri%AB_block_count
         atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
         atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
         size_opt_RI_AB = bs_env%auto_ri%AB_size_opt_RI(AB_block)
         RI_size_to_a = bs_env%auto_ri%AB_size_opt_RI_to_A(AB_block)
         IF (RI_size_to_a > 0) THEN
            position = n_AB_blocks_by_atom(atom_A) + 1
            n_AB_blocks_by_atom(atom_A) = position
            lists(atom_A)%block(position) = AB_block
            lists(atom_A)%first_p_AB(position) = 1
            lists(atom_A)%first_p_atom(position) = bs_env%auto_ri%AB_first_p_A(AB_block)
            lists(atom_A)%size_opt_RI(position) = RI_size_to_a
         END IF
         IF (atom_B /= atom_A .AND. RI_size_to_a < size_opt_RI_AB) THEN
            position = n_AB_blocks_by_atom(atom_B) + 1
            n_AB_blocks_by_atom(atom_B) = position
            lists(atom_B)%block(position) = AB_block
            lists(atom_B)%first_p_AB(position) = RI_size_to_a + 1
            lists(atom_B)%first_p_atom(position) = &
               bs_env%auto_ri%AB_first_p_B(AB_block)
            lists(atom_B)%size_opt_RI(position) = size_opt_RI_AB - RI_size_to_a
         END IF
      END DO
      DO atom_A = 1, bs_env%n_atom
         IF (SUM(lists(atom_A)%size_opt_RI) /= bs_env%auto_ri%sizes_opt_RI(atom_A)) THEN
            CPABORT("AUTO_RI AB-block columns do not fill the optimized atom block")
         END IF
      END DO
      DEALLOCATE (n_AB_blocks_by_atom)

      CALL timestop(handle)

   END SUBROUTINE build_AB_block_lists

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
!> \brief Collects U_Pp from the large reference basis for one optimized atom block.
!> \param bs_env ...
!> \param atom ...
!> \param list ...
!> \param atom_U ...
! **************************************************************************************************
   SUBROUTINE build_optimized_atom_U(bs_env, atom, list, atom_U)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom
      TYPE(AB_block_list_type), INTENT(IN)               :: list
      TYPE(optimized_atom_U_type), INTENT(OUT)           :: atom_U

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_optimized_atom_U'

      INTEGER :: AB_block, atom_A, atom_B, first, first_p_AB, first_p_atom, first_row, handle, i, &
         last, nRI_ref_A, nRI_ref_block, nsupport, position, size_opt_RI_AB
      INTEGER, ALLOCATABLE                               :: support_atoms(:)
      REAL(KIND=dp), ALLOCATABLE                         :: U_Pp_AB(:, :)

      CALL timeset(routineN, handle)

      ALLOCATE (support_atoms(2*SIZE(list%block)))
      nsupport = 0
      DO i = 1, SIZE(list%block)
         AB_block = list%block(i)
         atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
         atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
         position = find_atom_position(atom_A, support_atoms, nsupport)
         IF (position == 0) THEN
            nsupport = nsupport + 1
            support_atoms(nsupport) = atom_A
         END IF
         IF (atom_B /= atom_A) THEN
            position = find_atom_position(atom_B, support_atoms, nsupport)
            IF (position == 0) THEN
               nsupport = nsupport + 1
               support_atoms(nsupport) = atom_B
            END IF
         END IF
      END DO
      CPASSERT(nsupport > 0)
      ALLOCATE (atom_U%ref_atom(nsupport), atom_U%first_P(nsupport))
      atom_U%ref_atom(:) = support_atoms(:nsupport)
      first_row = 1
      DO i = 1, nsupport
         atom_U%first_P(i) = first_row
         first_row = first_row + get_ref_RI_size(bs_env, atom_U%ref_atom(i))
      END DO
      ALLOCATE (atom_U%U_Pp(first_row - 1, bs_env%auto_ri%sizes_opt_RI(atom)), &
                source=0.0_dp)

      DO i = 1, SIZE(list%block)
         AB_block = list%block(i)
         atom_A = bs_env%auto_ri%AB_atom_A(AB_block)
         atom_B = bs_env%auto_ri%AB_atom_B(AB_block)
         nRI_ref_block = bs_env%auto_ri%AB_size_ref_RI(AB_block)
         size_opt_RI_AB = bs_env%auto_ri%AB_size_opt_RI(AB_block)
         first = bs_env%auto_ri%U_Pp_AB_offset(AB_block)
         last = first + nRI_ref_block*size_opt_RI_AB - 1
         ALLOCATE (U_Pp_AB(nRI_ref_block, size_opt_RI_AB))
         U_Pp_AB(:, :) = RESHAPE(bs_env%auto_ri%U_Pp_AB(first:last), &
                                 [nRI_ref_block, size_opt_RI_AB])
         first_p_AB = list%first_p_AB(i)
         first_p_atom = list%first_p_atom(i)
         size_opt_RI_AB = list%size_opt_RI(i)
         IF (first_p_AB < 1 .OR. &
             first_p_AB + size_opt_RI_AB - 1 > SIZE(U_Pp_AB, 2)) THEN
            CPABORT("AUTO_RI AB-block column range is invalid")
         END IF
         IF (first_p_atom < 1 .OR. &
             first_p_atom + size_opt_RI_AB - 1 > bs_env%auto_ri%sizes_opt_RI(atom)) THEN
            CPABORT("AUTO_RI optimized atom-block column range is invalid")
         END IF

         nRI_ref_A = get_ref_RI_size(bs_env, atom_A)
         position = find_atom_position(atom_A, atom_U%ref_atom, nsupport)
         CPASSERT(position > 0)
         first_row = atom_U%first_P(position)
         atom_U%U_Pp(first_row:first_row + nRI_ref_A - 1, &
                     first_p_atom:first_p_atom + size_opt_RI_AB - 1) = &
            U_Pp_AB(1:nRI_ref_A, first_p_AB:first_p_AB + size_opt_RI_AB - 1)
         IF (atom_B /= atom_A) THEN
            position = find_atom_position(atom_B, atom_U%ref_atom, nsupport)
            CPASSERT(position > 0)
            first_row = atom_U%first_P(position)
            atom_U%U_Pp(first_row:first_row + get_ref_RI_size(bs_env, atom_B) - 1, &
                        first_p_atom:first_p_atom + size_opt_RI_AB - 1) = &
               U_Pp_AB(nRI_ref_A + 1:, first_p_AB:first_p_AB + size_opt_RI_AB - 1)
         END IF
         DEALLOCATE (U_Pp_AB)
      END DO
      DEALLOCATE (support_atoms)

      CALL timestop(handle)

   END SUBROUTINE build_optimized_atom_U

! **************************************************************************************************
!> \brief Adds REGULARIZATION_RI to one atom-local reference metric block.
!> \param bs_env ...
!> \param atom ...
!> \param M_PQ atom-local metric block in the large reference RI basis
! **************************************************************************************************
   SUBROUTINE regularize_M_PQ_block(bs_env, atom, M_PQ)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: M_PQ

      CHARACTER(LEN=*), PARAMETER :: routineN = 'regularize_M_PQ_block'

      INTEGER                                            :: handle, i, iset, kind, nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zet

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(M_PQ, 1) == SIZE(M_PQ, 2))
      kind = bs_env%ri_rs%particle_set(atom)%atomic_kind%kind_number
      CALL get_gto_basis_set(bs_env%basis_set_RI(kind)%gto_basis_set, &
                             first_sgf=first_sgf, npgf=npgf, nset=nset, &
                             nsgf_set=nsgf_set, zet=zet)
      DO iset = 1, nset
         DO i = first_sgf(1, iset), first_sgf(1, iset) + nsgf_set(iset) - 1
            M_PQ(i, i) = M_PQ(i, i) + &
                         bs_env%regularization_RI* &
                         MAX(1.0_dp, 1.0_dp/MINVAL(zet(1:npgf(iset), iset)))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE regularize_M_PQ_block

! **************************************************************************************************
!> \brief Directly contracts one optimized atom-pair metric block.
!>
!>        M_pq^{AB} = Σ_IJ (U_Ip^A)^T M_IJ U_Jq^B. Each reference
!>        atom-pair block M_IJ is discarded immediately after this contraction.
!> \param bs_env ...
!> \param atom_U_a contractions U_Ip for optimized block A
!> \param atom_U_b contractions U_Jq for optimized block B
!> \param coulomb true for Coulomb; false for the GW fitting operator
!> \param symmetric true when both transforms describe the same optimized block
!> \param matrix_optimized contracted M_pq or V_pq block in the optimized RI basis
! **************************************************************************************************
   SUBROUTINE compute_optimized_RI_matrix_block( &
      bs_env, atom_U_a, atom_U_b, coulomb, symmetric, matrix_optimized)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(optimized_atom_U_type), INTENT(IN)            :: atom_U_a, atom_U_b
      LOGICAL, INTENT(IN)                                :: coulomb, symmetric
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: matrix_optimized

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_optimized_RI_matrix_block'

      INTEGER                                            :: atom_A, atom_B, first_a, first_b, &
                                                            handle, i, j, last_a, last_b, &
                                                            nRI_small_A, nRI_small_B
      REAL(KIND=dp), ALLOCATABLE :: matrix_pq_contribution(:, :), matrix_reference(:, :), &
         matrix_reference_times_U_b(:, :)

      CALL timeset(routineN, handle)

      nRI_small_A = SIZE(atom_U_a%U_Pp, 2)
      nRI_small_B = SIZE(atom_U_b%U_Pp, 2)
      IF (symmetric) THEN
         CPASSERT(nRI_small_A == nRI_small_B)
      END IF
      ALLOCATE (matrix_optimized(nRI_small_A, nRI_small_B), source=0.0_dp)
      DO i = 1, SIZE(atom_U_a%ref_atom)
         atom_A = atom_U_a%ref_atom(i)
         first_a = atom_U_a%first_P(i)
         last_a = first_a + get_ref_RI_size(bs_env, atom_A) - 1
         DO j = 1, SIZE(atom_U_b%ref_atom)
            IF (symmetric .AND. j > i) CYCLE
            atom_B = atom_U_b%ref_atom(j)
            first_b = atom_U_b%first_P(j)
            last_b = first_b + get_ref_RI_size(bs_env, atom_B) - 1
            IF (coulomb) THEN
               CALL compute_M_PQ_AB(bs_env, atom_A, atom_B, matrix_reference)
            ELSE
               CALL compute_M_PQ_AB( &
                  bs_env, atom_A, atom_B, matrix_reference, bs_env%ri_metric)
               IF (atom_A == atom_B .AND. bs_env%regularization_RI > 0.0_dp) THEN
                  CALL regularize_M_PQ_block(bs_env, atom_A, matrix_reference)
               END IF
            END IF
            ALLOCATE (matrix_reference_times_U_b(SIZE(matrix_reference, 1), nRI_small_B))
            matrix_reference_times_U_b(:, :) = &
               MATMUL(matrix_reference, atom_U_b%U_Pp(first_b:last_b, :))
            ALLOCATE (matrix_pq_contribution(nRI_small_A, nRI_small_B))
            matrix_pq_contribution(:, :) = &
               MATMUL(TRANSPOSE(atom_U_a%U_Pp(first_a:last_a, :)), &
                      matrix_reference_times_U_b)
            matrix_optimized(:, :) = matrix_optimized + matrix_pq_contribution
            IF (symmetric .AND. i /= j) THEN
               matrix_optimized(:, :) = matrix_optimized + TRANSPOSE(matrix_pq_contribution)
            END IF
            DEALLOCATE (matrix_pq_contribution, matrix_reference, matrix_reference_times_U_b)
         END DO
      END DO
      IF (.NOT. ALL(ieee_is_finite(matrix_optimized))) THEN
         CPABORT("AUTO_RI optimized metric block contains invalid values")
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_optimized_RI_matrix_block

! **************************************************************************************************
!> \brief Builds M_opt or V_opt directly from contracted atom-pair blocks.
!> \param qs_env ...
!> \param bs_env ...
!> \param lists ...
!> \param matrix_pq contracted M_pq or V_pq in the optimized RI basis
!> \param coulomb true for Coulomb; false for the GW fitting operator
! **************************************************************************************************
   SUBROUTINE compute_contracted_RI_matrix(qs_env, bs_env, lists, matrix_pq, coulomb)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(AB_block_list_type), DIMENSION(:), INTENT(IN) :: lists
      TYPE(cp_fm_type), INTENT(OUT)                      :: matrix_pq
      LOGICAL, INTENT(IN)                                :: coulomb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_contracted_RI_matrix'

      INTEGER                                            :: handle, iatom, jatom, npcol, nprow, &
                                                            pair_index, RI_size_opt
      INTEGER, POINTER                                   :: col_dist(:), row_dist(:), sizes_opt_RI(:)
      LOGICAL                                            :: atom_U_a_ready
      REAL(KIND=dp), ALLOCATABLE                         :: matrix_pq_block(:, :)
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_distribution_type), POINTER             :: template
      TYPE(dbcsr_type)                                   :: blocks
      TYPE(optimized_atom_U_type)                        :: atom_U_a, atom_U_b

      CALL timeset(routineN, handle)

      ALLOCATE (sizes_opt_RI(bs_env%n_atom), row_dist(bs_env%n_atom), col_dist(bs_env%n_atom))
      sizes_opt_RI(:) = bs_env%auto_ri%sizes_opt_RI
      CALL get_qs_env(qs_env, dbcsr_dist=template)
      CALL dbcsr_distribution_get(template, nprows=nprow, npcols=npcol)
      DO iatom = 1, SIZE(sizes_opt_RI)
         row_dist(iatom) = MOD(iatom - 1, nprow)
         col_dist(iatom) = MOD(iatom - 1, npcol)
      END DO
      RI_size_opt = SUM(sizes_opt_RI)
      CALL create_distributed_matrix(qs_env, matrix_pq, RI_size_opt, RI_size_opt)
      CALL dbcsr_distribution_new(dist, template=template, row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(blocks, name='AUTO_RI contracted metric', dist=dist, &
                        matrix_type=dbcsr_type_symmetric, &
                        row_blk_size=sizes_opt_RI, col_blk_size=sizes_opt_RI)

      pair_index = 0
      DO iatom = 1, SIZE(sizes_opt_RI)
         atom_U_a_ready = .FALSE.
         DO jatom = 1, iatom
            pair_index = pair_index + 1
            IF (MOD(pair_index - 1, bs_env%para_env%num_pe) /= bs_env%para_env%mepos) CYCLE
            IF (.NOT. atom_U_a_ready) THEN
               CALL build_optimized_atom_U(bs_env, iatom, lists(iatom), atom_U_a)
               atom_U_a_ready = .TRUE.
            END IF
            IF (iatom == jatom) THEN
               CALL compute_optimized_RI_matrix_block( &
                  bs_env, atom_U_a, atom_U_a, coulomb, .TRUE., &
                  matrix_pq_block)
               matrix_pq_block(:, :) = 0.5_dp*(matrix_pq_block + TRANSPOSE(matrix_pq_block))
            ELSE
               CALL build_optimized_atom_U(bs_env, jatom, lists(jatom), atom_U_b)
               CALL compute_optimized_RI_matrix_block( &
                  bs_env, atom_U_a, atom_U_b, coulomb, .FALSE., &
                  matrix_pq_block)
            END IF
            CALL dbcsr_put_block(blocks, iatom, jatom, matrix_pq_block)
            DEALLOCATE (matrix_pq_block)
         END DO
      END DO
      CALL dbcsr_finalize(blocks)
      CALL copy_dbcsr_to_fm(blocks, matrix_pq)
      CALL dbcsr_release(blocks)
      CALL dbcsr_distribution_release(dist)
      DEALLOCATE (sizes_opt_RI, row_dist, col_dist)

      CALL timestop(handle)

   END SUBROUTINE compute_contracted_RI_matrix

! **************************************************************************************************
!> \brief Computes the pseudoinverse of the optimized fitting metric M_opt.
!>
!>        The generalized problem M_opt X=M_AA X λ separates linearly
!>        dependent directions. The result is M_opt^+=X diag(1/λ) X^T.
!> \param qs_env ...
!> \param bs_env ...
!> \param M_pq optimized RI fitting metric
! **************************************************************************************************
   SUBROUTINE compute_M_pq_inv(qs_env, bs_env, M_pq)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: M_pq

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_M_pq_inv'

      INTEGER                                            :: handle, ncol_local, nrow_local, &
                                                            RI_size_opt
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), ALLOCATABLE                         :: eigenval(:), factors(:)
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: M_pq_AA_local, M_pq_local
      TYPE(cp_fm_type)                                   :: M_pq_AA, work, X_pq

      CALL timeset(routineN, handle)

      RI_size_opt = SUM(bs_env%auto_ri%sizes_opt_RI)
      CALL create_distributed_matrix(qs_env, bs_env%auto_ri%M_pq_inv, &
                                     RI_size_opt, RI_size_opt)
      CALL create_distributed_matrix(qs_env, X_pq, RI_size_opt, RI_size_opt)
      CALL create_distributed_matrix(qs_env, M_pq_AA, RI_size_opt, RI_size_opt)
      CALL create_distributed_matrix(qs_env, work, RI_size_opt, RI_size_opt)
      CALL cp_fm_get_info(M_pq, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices, &
                          local_data=M_pq_local)
      CALL cp_fm_get_info(M_pq_AA, local_data=M_pq_AA_local)
      M_pq_AA_local = 0.0_dp
      CALL copy_M_pq_AA_blocks(M_pq_local(:nrow_local, :ncol_local), &
                               M_pq_AA_local(:nrow_local, :ncol_local), &
                               row_indices(:nrow_local), col_indices(:ncol_local), &
                               bs_env%auto_ri%sizes_opt_RI)

      ALLOCATE (eigenval(RI_size_opt))
      CALL cp_fm_geeig(M_pq, M_pq_AA, X_pq, eigenval, work)
      CALL compute_optimized_M_pq_inv_factors(eigenval, factors)
      CALL cp_fm_column_scale(X_pq, factors)
      CALL cp_fm_syrk("U", "N", RI_size_opt, 1.0_dp, X_pq, 1, 1, 0.0_dp, &
                      bs_env%auto_ri%M_pq_inv)
      CALL cp_fm_uplo_to_full(bs_env%auto_ri%M_pq_inv, work)

      DEALLOCATE (eigenval, factors)
      CALL cp_fm_release(X_pq)
      CALL cp_fm_release(M_pq_AA)
      CALL cp_fm_release(work)

      CALL timestop(handle)

   END SUBROUTINE compute_M_pq_inv

! **************************************************************************************************
!> \brief Computes M_opt and V_opt without constructing M_ref or V_ref.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_M_pq_inv_and_V_pq(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_M_pq_inv_and_V_pq'

      INTEGER                                            :: handle
      TYPE(AB_block_list_type), ALLOCATABLE              :: lists(:)
      TYPE(cp_fm_type)                                   :: M_pq

      CALL timeset(routineN, handle)

      CALL build_AB_block_lists(bs_env, lists)
      CALL initialize_auto_ri_integrals(bs_env)
      CALL compute_contracted_RI_matrix(qs_env, bs_env, lists, M_pq, .FALSE.)
      CALL compute_M_pq_inv(qs_env, bs_env, M_pq)
      CALL cp_fm_release(M_pq)
      CALL compute_contracted_RI_matrix(qs_env, bs_env, lists, bs_env%auto_ri%V_pq, .TRUE.)

      ! Activate the optimized per-atom RI block sizes for RI-RS.
      CALL set_optimized_RI_basis_sizes(bs_env)

      CALL timestop(handle)

   END SUBROUTINE compute_M_pq_inv_and_V_pq

! **************************************************************************************************
!> \brief Computes the atom-local two-center RI metric (P|Q) and the
!>        three-center integrals (μν|P) used in Eqs. (3) and (4).
!> \param ao_basis ...
!> \param ri_basis ...
!> \param potential fitting operator used for both integral types
!> \param M_PQ two-center metric of the large reference RI basis
!> \param Int_3c three-center integrals (μν|P)
! **************************************************************************************************
   SUBROUTINE compute_int_3c(ao_basis, ri_basis, potential, M_PQ, Int_3c)
      TYPE(gto_basis_set_type), INTENT(IN)               :: ao_basis, ri_basis
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: M_PQ
      REAL(KIND=dp), INTENT(OUT)                         :: Int_3c(:, :, :)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_int_3c'

      INTEGER                                            :: handle, iset, jset, kset, max_ao_l, &
                                                            max_l, nAO, ncoa, ncob, ncoc, nRI, &
                                                            nsgfa, nsgfb, nsgfc, sgfa, sgfb, sgfc
      INTEGER, DIMENSION(:), POINTER                     :: lmax_a, lmax_c, lmin_a, lmin_c, npgf_a, &
                                                            npgf_c, nsgf_set_a, nsgf_set_c
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_a, first_sgf_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Int_2c_cart, Int_2c_spherical
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Int_3c_cart, Int_3c_spherical_all, &
                                                            Int_3c_spherical_block
      REAL(KIND=dp), DIMENSION(3)                        :: center
      REAL(KIND=dp), DIMENSION(:), POINTER               :: radius_a, radius_b, radius_c
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pgf_radius_a, pgf_radius_c, sphi_a, &
                                                            sphi_c, zeta_a, zeta_c
      TYPE(coulomb_integral_context_type)                :: context

      CALL timeset(routineN, handle)

      CALL get_gto_basis_set(ao_basis, nsgf=nAO, lmax=lmax_a, lmin=lmin_a, &
                             npgf=npgf_a, nsgf_set=nsgf_set_a, first_sgf=first_sgf_a, &
                             pgf_radius=pgf_radius_a, sphi=sphi_a, zet=zeta_a)
      CALL get_gto_basis_set(ri_basis, nsgf=nRI, lmax=lmax_c, lmin=lmin_c, &
                             npgf=npgf_c, nsgf_set=nsgf_set_c, first_sgf=first_sgf_c, &
                             pgf_radius=pgf_radius_c, sphi=sphi_c, zet=zeta_c)
      CPASSERT(nAO > 0 .AND. nRI > 0)
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

      ALLOCATE (M_PQ(nRI, nRI), Int_3c_spherical_all(nAO, nAO, nRI))
      M_PQ = 0.0_dp
      Int_3c_spherical_all = 0.0_dp

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
            ALLOCATE (Int_2c_cart(ncoa, ncob), Int_2c_spherical(nsgfa, nsgfb))
            Int_2c_cart = 0.0_dp
            CALL compute_coulomb_2c(context, lmin_c(iset), lmax_c(iset), &
                                    lmin_c(jset), lmax_c(jset), &
                                    npgf_c(iset), npgf_c(jset), zeta_c(1:npgf_c(iset), iset), &
                                    zeta_c(1:npgf_c(jset), jset), &
                                    radius_a, radius_b, center, center, Int_2c_cart, potential)
            CALL ab_contract(Int_2c_spherical, Int_2c_cart, sphi_c(:, sgfa:), sphi_c(:, sgfb:), &
                             ncoa, ncob, nsgfa, nsgfb)
            M_PQ(sgfa:sgfa + nsgfa - 1, sgfb:sgfb + nsgfb - 1) = Int_2c_spherical
            DEALLOCATE (Int_2c_spherical, Int_2c_cart)
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
               ALLOCATE (Int_3c_cart(ncoa, ncob, ncoc), &
                         Int_3c_spherical_block(nsgfa, nsgfb, nsgfc))
               Int_3c_cart = 0.0_dp
               CALL compute_coulomb_3c(context, lmin_a(iset), lmax_a(iset), &
                                       lmin_a(jset), lmax_a(jset), &
                                       lmin_c(kset), lmax_c(kset), npgf_a(iset), npgf_a(jset), &
                                       npgf_c(kset), zeta_a(1:npgf_a(iset), iset), &
                                       zeta_a(1:npgf_a(jset), jset), zeta_c(1:npgf_c(kset), kset), &
                                       radius_a, radius_b, radius_c, center, center, center, &
                                       Int_3c_cart, potential)
               CALL abc_contract(Int_3c_spherical_block, Int_3c_cart, &
                                 sphi_a(:, sgfa:), sphi_a(:, sgfb:), &
                                 sphi_c(:, sgfc:), ncoa, ncob, ncoc, nsgfa, nsgfb, nsgfc)
               Int_3c_spherical_all(sgfa:sgfa + nsgfa - 1, sgfb:sgfb + nsgfb - 1, &
                                    sgfc:sgfc + nsgfc - 1) = Int_3c_spherical_block
               DEALLOCATE (Int_3c_spherical_block, Int_3c_cart)
            END DO
         END DO
      END DO

      Int_3c(:, :, :) = Int_3c_spherical_all
      CALL coulomb_integral_cleanup(context)
      IF (.NOT. ALL(ieee_is_finite(M_PQ))) THEN
         CPABORT("AUTO_RI metric contains invalid values")
      END IF
      IF (.NOT. ALL(ieee_is_finite(Int_3c))) THEN
         CPABORT("AUTO_RI three-center integrals contain invalid values")
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_int_3c

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
!> \param atom_A ...
!> \param atom_B ...
!> \param M_PQ metric block in the large reference RI basis
!> \param ri_potential ...
! **************************************************************************************************
   SUBROUTINE compute_M_PQ_AB(bs_env, atom_A, atom_B, M_PQ, ri_potential)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: atom_A, atom_B
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: M_PQ
      TYPE(libint_potential_type), INTENT(IN), OPTIONAL  :: ri_potential

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_M_PQ_AB'

      INTEGER                                            :: atom_1, atom_2, handle, kind_1, kind_2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: M_PQ_canonical
      TYPE(libint_potential_type)                        :: potential

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(bs_env%ri_rs%particle_set))
      CPASSERT(atom_A >= 1 .AND. atom_A <= bs_env%n_atom)
      CPASSERT(atom_B >= 1 .AND. atom_B <= bs_env%n_atom)
      IF (atom_A >= atom_B) THEN
         atom_1 = atom_A
         atom_2 = atom_B
      ELSE
         atom_1 = atom_B
         atom_2 = atom_A
      END IF
      kind_1 = bs_env%ri_rs%particle_set(atom_1)%atomic_kind%kind_number
      kind_2 = bs_env%ri_rs%particle_set(atom_2)%atomic_kind%kind_number

      IF (PRESENT(ri_potential)) THEN
         CPASSERT(ri_potential%potential_type == do_potential_truncated)
         potential = ri_potential
      ELSE
         potential%potential_type = do_potential_coulomb
      END IF
      CALL compute_M_PQ( &
         basis_a=bs_env%basis_set_RI(kind_1)%gto_basis_set, &
         basis_b=bs_env%basis_set_RI(kind_2)%gto_basis_set, &
         center_a=bs_env%ri_rs%particle_set(atom_1)%r, &
         center_b=bs_env%ri_rs%particle_set(atom_2)%r, &
         potential=potential, M_PQ=M_PQ_canonical)

      IF (atom_A >= atom_B) THEN
         CALL MOVE_ALLOC(M_PQ_canonical, M_PQ)
      ELSE
         ALLOCATE (M_PQ(SIZE(M_PQ_canonical, 2), SIZE(M_PQ_canonical, 1)))
         M_PQ(:, :) = TRANSPOSE(M_PQ_canonical)
         DEALLOCATE (M_PQ_canonical)
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_M_PQ_AB

! **************************************************************************************************
!> \brief Computes one contracted two-center RI metric block (P_A|Q_B).
!> \param basis_a ...
!> \param basis_b ...
!> \param center_a ...
!> \param center_b ...
!> \param potential ...
!> \param M_PQ metric block in the large reference RI basis
! **************************************************************************************************
   SUBROUTINE compute_M_PQ(basis_a, basis_b, center_a, center_b, potential, M_PQ)
      TYPE(gto_basis_set_type), POINTER                  :: basis_a, basis_b
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: center_a, center_b
      TYPE(libint_potential_type), INTENT(IN)            :: potential
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: M_PQ

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_M_PQ'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (M_PQ(basis_a%nsgf, basis_b%nsgf), source=0.0_dp)
      SELECT CASE (potential%potential_type)
      CASE (do_potential_truncated)
         CALL int_operators_r12_ab_os(operator_truncated, M_PQ, rab=center_b - center_a, &
                                      fba=basis_a, fbb=basis_b, &
                                      r_cutoff=potential%cutoff_radius, calculate_forces=.FALSE.)
      CASE (do_potential_coulomb)
         CALL int_operators_r12_ab_os(operator_coulomb, M_PQ, rab=center_b - center_a, &
                                      fba=basis_a, fbb=basis_b, calculate_forces=.FALSE.)
      CASE DEFAULT
         CPABORT("AUTO_RI two-center metric requires COULOMB or TRUNCATED")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE compute_M_PQ

! **************************************************************************************************
!> \brief Returns the reference RI basis size of one atom.
!> \param bs_env ...
!> \param iatom ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION get_ref_RI_size(bs_env, iatom) RESULT(nRI_ref)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: iatom

      INTEGER                                            :: ikind

      ikind = bs_env%ri_rs%particle_set(iatom)%atomic_kind%kind_number
      nRI_ref = bs_env%basis_set_RI(ikind)%gto_basis_set%nsgf

   END FUNCTION get_ref_RI_size

END MODULE gw_auto_ri_generate_RI_basis
