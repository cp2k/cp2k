!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate MAO's and analyze wavefunctions
!> \par History
!>      03.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE mao_methods
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: create_primitive_basis_set,&
                                              get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type,&
                                              write_gto_basis_set
   USE bibliography,                    ONLY: Ehrhardt1985,&
                                              Heinzmann1976,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_restore
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, cp_dbcsr_create, &
        cp_dbcsr_deallocate_matrix_set, cp_dbcsr_desymmetrize, cp_dbcsr_distribution, &
        cp_dbcsr_get_block_diag, cp_dbcsr_get_block_p, cp_dbcsr_get_info, cp_dbcsr_init, &
        cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
        cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_norm, &
        cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_replicate_all, cp_dbcsr_reserve_diag_blocks, &
        cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_type, dbcsr_distribution_mp, &
        dbcsr_distribution_obj, dbcsr_mp_group, dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm, &
        dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_plus_fm_fm_t
   USE cp_fm_diag,                      ONLY: cp_fm_geeig
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: mao_basis_ext,&
                                              mao_basis_orb,&
                                              mao_basis_prim
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE lapack,                          ONLY: lapack_ssyev,&
                                              lapack_ssygv
   USE mathlib,                         ONLY: invmat_symm
   USE message_passing,                 ONLY: mp_sum
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_overlap,                      ONLY: build_overlap_matrix_simple
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE block_type
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE  :: mat
   END TYPE block_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mao_methods'

   PUBLIC ::  mao_analysis, mao_basis_simple

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input_section ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE mao_analysis(qs_env, input_section, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_section
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_analysis', routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol, esa, esb, esc
      INTEGER :: fall, handle, ia, iab, iabc, iatom, ib, ic, icol, ikind, irow, ispin, iter, &
         jatom, mao_basis, max_iter, me, na, nab, nabc, natom, nb, nc, nimages, nspin, ssize
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, mao_blk, mao_blk_sizes, &
                                                            orb_blk, row_blk_sizes
      LOGICAL                                            :: analyze_ua, explicit, fo, for, fos, &
                                                            found, neglect_abc, print_basis
      REAL(KIND=dp) :: a1, a2, alpha, an, beta, deltaq, electra(2), eps_ab, eps_abc, eps_filter, &
         eps_fun, eps_grad, fa1, fa2, fnnew, fnold, fval, grad_norm, senabc, senmax, threshold, &
         total_charge, total_spin, ua_charge(2), zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: occnumA, occnumABC, qab, qmatab, qmatac, &
                                                            qmatbc, raq, sab, selnABC, sinv, &
                                                            smatab, smatac, smatbc, uaq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: occnumAB, selnAB
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block, cmao, diag, qblka, qblkb, qblkc, &
                                                            rblkl, rblku, sblk, sblka, sblkb, sblkc
      TYPE(block_type), ALLOCATABLE, DIMENSION(:)        :: rowblock
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_coef, mao_dmat, mao_grad, mao_qmat, &
                                                            mao_smat, matrix_q, matrix_smm, &
                                                            matrix_smo
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: matrix_ks, matrix_p, matrix_s
      TYPE(cp_dbcsr_type)                                :: amat, axmat, cgmat, cholmat, crumat, &
                                                            qmat, qmat_diag, rumat, smat_diag, &
                                                            sumat, tmat
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: mao_basis_set_list, orb_basis_set_list
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb, smm_list, smo_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

! only do MAO analysis if explicitely requested

      CALL section_vals_get(input_section, explicit=explicit)
      IF (.NOT. explicit) RETURN

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
         WRITE (UNIT=unit_nr, FMT="(T36,A)") "MAO ANALYSIS"
         WRITE (UNIT=unit_nr, FMT="(T12,A)") "Claus Ehrhardt and Reinhart Ahlrichs, TCA 68:231-245 (1985)"
         WRITE (unit_nr, '(T2,A)') '!-----------------------------------------------------------------------------!'
      END IF
      CALL cite_reference(Heinzmann1976)
      CALL cite_reference(Ehrhardt1985)

      ! input options
      CALL section_vals_val_get(input_section, "REFERENCE_BASIS", i_val=mao_basis)
      CALL section_vals_val_get(input_section, "EPS_FILTER", r_val=eps_filter)
      CALL section_vals_val_get(input_section, "EPS_FUNCTION", r_val=eps_fun)
      CALL section_vals_val_get(input_section, "EPS_GRAD", r_val=eps_grad)
      CALL section_vals_val_get(input_section, "MAX_ITER", i_val=max_iter)
      CALL section_vals_val_get(input_section, "PRINT_BASIS", l_val=print_basis)
      CALL section_vals_val_get(input_section, "NEGLECT_ABC", l_val=neglect_abc)
      CALL section_vals_val_get(input_section, "AB_THRESHOLD", r_val=eps_ab)
      CALL section_vals_val_get(input_section, "ABC_THRESHOLD", r_val=eps_abc)
      CALL section_vals_val_get(input_section, "ANALYZE_UNASSIGNED_CHARGE", l_val=analyze_ua)

      ! k-points?
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nimages = dft_control%nimages
      IF (nimages > 1) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A)") &
               "K-Points: MAO's determined and analyzed using Gamma-Point only."
         END IF
      END IF

      ! Reference basis set
      NULLIFY (mao_basis_set_list, orb_basis_set_list)
      CALL mao_reference_basis(qs_env, mao_basis, mao_basis_set_list, orb_basis_set_list, &
                               unit_nr, print_basis)

      ! neighbor lists
      NULLIFY (smm_list, smo_list)
      CALL setup_neighbor_list(smm_list, mao_basis_set_list, qs_env=qs_env)
      CALL setup_neighbor_list(smo_list, mao_basis_set_list, orb_basis_set_list, qs_env=qs_env)

      ! overlap matrices
      NULLIFY (matrix_smm, matrix_smo)
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL build_overlap_matrix_simple(ks_env, matrix_smm, &
                                       mao_basis_set_list, mao_basis_set_list, smm_list)
      CALL build_overlap_matrix_simple(ks_env, matrix_smo, &
                                       mao_basis_set_list, orb_basis_set_list, smo_list)

      ! get reference density matrix and overlap matrix
      CALL get_qs_env(qs_env, rho=rho, matrix_s_kp=matrix_s)
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      nspin = SIZE(matrix_p, 1)
      !
      ! Q matrix
      IF (nimages == 1) THEN
         CALL mao_build_q(matrix_q, matrix_p, matrix_s, matrix_smm, matrix_smo, smm_list, electra, eps_filter)
      ELSE
         CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, kpoints=kpoints)
         CALL mao_build_q(matrix_q, matrix_p, matrix_s, matrix_smm, matrix_smo, smm_list, electra, eps_filter, &
                          nimages=nimages, kpoints=kpoints, matrix_ks=matrix_ks, sab_orb=sab_orb)
      END IF

      ! check for extended basis sets
      fall = 0
      CALL neighbor_list_iterator_create(nl_iterator, smm_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom)
         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF
         CALL cp_dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
                                   row=irow, col=icol, block=block, found=found)
         IF (.NOT. found) fall = fall+1
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL get_qs_env(qs_env=qs_env, para_env=para_env)
      CALL mp_sum(fall, para_env%group)
      IF (unit_nr > 0 .AND. fall > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(/,T2,A,/,T2,A,/)") &
            "Warning: Extended MAO basis used with original basis filtered density matrix", &
            "Warning: Possible errors can be controled with EPS_PGF_ORB"
      END IF

      ! MAO matrices
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, natom=natom)
      CALL get_ks_env(ks_env=ks_env, particle_set=particle_set, dbcsr_dist=dbcsr_dist)
      NULLIFY (mao_coef, mao_grad)
      CALL cp_dbcsr_allocate_matrix_set(mao_coef, nspin)
      CALL cp_dbcsr_allocate_matrix_set(mao_grad, nspin)
      ALLOCATE (row_blk_sizes(natom), col_blk_sizes(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, &
                            basis=mao_basis_set_list)
      CALL get_particle_set(particle_set, qs_kind_set, nmao=col_blk_sizes)
      ! check if MAOs have been specified
      DO iab = 1, natom
         IF (col_blk_sizes(iab) < 0) &
            CPABORT("Number of MAOs has to be specified in KIND section for all elements")
      END DO
      DO ispin = 1, nspin
         ! coeficients
         ALLOCATE (mao_coef(ispin)%matrix)
         CALL cp_dbcsr_init(mao_coef(ispin)%matrix)
         CALL cp_dbcsr_create(matrix=mao_coef(ispin)%matrix, &
                              name="MAO_COEF", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL cp_dbcsr_reserve_diag_blocks(matrix=mao_coef(ispin)%matrix)
         ! gradients
         ALLOCATE (mao_grad(ispin)%matrix)
         CALL cp_dbcsr_init(mao_grad(ispin)%matrix)
         CALL cp_dbcsr_create(matrix=mao_grad(ispin)%matrix, &
                              name="MAO_GRAD", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL cp_dbcsr_reserve_diag_blocks(matrix=mao_grad(ispin)%matrix)
      END DO
      DEALLOCATE (row_blk_sizes, col_blk_sizes)

      ! initialize MAO coeficients from diagonal blocks of the Q matrix
      DO ispin = 1, nspin
         CALL mao_initialization(mao_coef(ispin)%matrix, &
                                 matrix_q(ispin)%matrix, matrix_smm(1)%matrix)
      END DO

      IF (max_iter < 1) THEN
         ! projection only
         CALL cp_dbcsr_get_info(mao_coef(1)%matrix, col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
         CALL cp_dbcsr_init(tmat)
         CALL cp_dbcsr_create(tmat, name="Binv", dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                              row_blk_size=col_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         IF (unit_nr > 0) WRITE (unit_nr, *)
         DO ispin = 1, nspin
            CALL mao_function(mao_coef(ispin)%matrix, fval, matrix_q(ispin)%matrix, &
                              matrix_smm(1)%matrix, tmat, .FALSE.)
            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(T2,A,T20,A,I2,T56,A,F12.8)") &
                  "MAO Projection", "Spin =", ispin, "Completness =", fval/electra(ispin)
            END IF
         END DO
         CALL cp_dbcsr_release(tmat)
      ELSE
         ! optimize MAOs
         alpha = 0.25_dp
         beta = 0.0_dp
         CALL cp_dbcsr_get_info(mao_coef(1)%matrix, col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
         CALL cp_dbcsr_init(tmat)
         CALL cp_dbcsr_create(tmat, name="Binv", dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                              row_blk_size=col_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL cp_dbcsr_init(cgmat)
         CALL cp_dbcsr_create(cgmat, template=mao_grad(1)%matrix)
         CALL cp_dbcsr_init(amat)
         CALL cp_dbcsr_create(amat, template=mao_coef(1)%matrix)
         DO ispin = 1, nspin
            CALL mao_function_gradient(mao_coef(ispin)%matrix, fval, mao_grad(ispin)%matrix, &
                                       matrix_q(ispin)%matrix, matrix_smm(1)%matrix, tmat, .FALSE.)
            CALL cp_dbcsr_copy(cgmat, mao_grad(ispin)%matrix)
            CALL cp_dbcsr_norm(mao_grad(ispin)%matrix, dbcsr_norm_maxabsnorm, norm_scalar=grad_norm)
            CALL cp_dbcsr_norm(mao_grad(ispin)%matrix, dbcsr_norm_frobenius, norm_scalar=fnold)
            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(/,T2,A,T73,A,I2)") "MAO OPTIMIZATION", "Spin =", ispin
               WRITE (UNIT=unit_nr, FMT="(T2,A,T24,A,F11.8,T48,A,F11.8,T69,A,F6.3)") &
                  "Initialization", "fval =", fval/electra(ispin), "grad =", grad_norm, "step =", alpha
            END IF
            DO iter = 1, max_iter
               IF (grad_norm < eps_grad) EXIT
               IF ((1.0_dp-fval/electra(ispin)) < eps_fun) EXIT
               CALL cp_dbcsr_add(mao_coef(ispin)%matrix, cgmat, 1.0_dp, alpha)
               CALL mao_orthogonalization(mao_coef(ispin)%matrix, matrix_smm(1)%matrix)
               CALL mao_function_gradient(mao_coef(ispin)%matrix, fval, mao_grad(ispin)%matrix, &
                                          matrix_q(ispin)%matrix, matrix_smm(1)%matrix, tmat, .TRUE.)
               CALL cp_dbcsr_norm(mao_grad(ispin)%matrix, dbcsr_norm_maxabsnorm, norm_scalar=grad_norm)
               CALL cp_dbcsr_norm(mao_grad(ispin)%matrix, dbcsr_norm_frobenius, norm_scalar=fnnew)
               IF (unit_nr > 0) THEN
                  WRITE (UNIT=unit_nr, FMT="(T2,A,i8,T24,A,F11.8,T48,A,F11.8,T69,A,F6.3)") &
                     "iter=", iter, "fval =", fval/electra(ispin), "grad =", grad_norm, "step =", alpha
               END IF
               beta = fnnew/fnold
               CALL cp_dbcsr_add(cgmat, mao_grad(ispin)%matrix, beta, 1.0_dp)
               fnold = fnnew
               ! line search, update alpha
               CALL cp_dbcsr_copy(amat, mao_coef(ispin)%matrix)
               CALL cp_dbcsr_add(amat, cgmat, 1.0_dp, 0.5_dp*alpha)
               CALL mao_orthogonalization(amat, matrix_smm(1)%matrix)
               CALL mao_function(amat, fa1, matrix_q(ispin)%matrix, matrix_smm(1)%matrix, tmat, .TRUE.)
               CALL cp_dbcsr_copy(amat, mao_coef(ispin)%matrix)
               CALL cp_dbcsr_add(amat, cgmat, 1.0_dp, alpha)
               CALL mao_orthogonalization(amat, matrix_smm(1)%matrix)
               CALL mao_function(amat, fa2, matrix_q(ispin)%matrix, matrix_smm(1)%matrix, tmat, .TRUE.)
               a2 = (4._dp*fa1-fa2-3._dp*fval)/alpha
               a1 = (fa2-fval-a2*alpha)/(alpha*alpha)
               an = -a2/(2._dp*a1)
               an = MIN(an, 2.0_dp*alpha)
               IF (an < 0.05_dp .OR. a1 > 0.0_dp) THEN
                  CALL cp_dbcsr_copy(cgmat, mao_grad(ispin)%matrix)
                  alpha = 0.25_dp
               ELSE
                  alpha = an
               END IF
            END DO
         END DO
         CALL cp_dbcsr_release(tmat)
         CALL cp_dbcsr_release(cgmat)
         CALL cp_dbcsr_release(amat)
      END IF

      ! Analyze the MAO basis
      CALL mao_basis_analysis(mao_coef, matrix_smm, mao_basis_set_list, particle_set, &
                              qs_kind_set, unit_nr, para_env)

      ! Calculate the overlap and density matrix in the new MAO basis
      NULLIFY (mao_dmat, mao_smat, mao_qmat)
      CALL cp_dbcsr_allocate_matrix_set(mao_qmat, nspin)
      CALL cp_dbcsr_allocate_matrix_set(mao_dmat, nspin)
      CALL cp_dbcsr_allocate_matrix_set(mao_smat, nspin)
      CALL cp_dbcsr_get_info(mao_coef(1)%matrix, col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
      DO ispin = 1, nspin
         ALLOCATE (mao_dmat(ispin)%matrix)
         CALL cp_dbcsr_init(mao_dmat(ispin)%matrix)
         CALL cp_dbcsr_create(mao_dmat(ispin)%matrix, name="MAO density", dist=dbcsr_dist, &
                              matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                              col_blk_size=col_blk_sizes, nze=0)
         ALLOCATE (mao_smat(ispin)%matrix)
         CALL cp_dbcsr_init(mao_smat(ispin)%matrix)
         CALL cp_dbcsr_create(mao_smat(ispin)%matrix, name="MAO overlap", dist=dbcsr_dist, &
                              matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                              col_blk_size=col_blk_sizes, nze=0)
         ALLOCATE (mao_qmat(ispin)%matrix)
         CALL cp_dbcsr_init(mao_qmat(ispin)%matrix)
         CALL cp_dbcsr_create(mao_qmat(ispin)%matrix, name="MAO covar density", dist=dbcsr_dist, &
                              matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                              col_blk_size=col_blk_sizes, nze=0)
      END DO
      CALL cp_dbcsr_init(amat)
      CALL cp_dbcsr_create(amat, name="MAO overlap", template=mao_dmat(1)%matrix)
      CALL cp_dbcsr_init(tmat)
      CALL cp_dbcsr_create(tmat, name="MAO Overlap Inverse", template=amat)
      CALL cp_dbcsr_init(qmat)
      CALL cp_dbcsr_create(qmat, name="MAO covar density", template=amat)
      CALL cp_dbcsr_init(cgmat)
      CALL cp_dbcsr_create(cgmat, name="TEMP matrix", template=mao_coef(1)%matrix)
      CALL cp_dbcsr_init(axmat)
      CALL cp_dbcsr_create(axmat, name="TEMP", template=amat, matrix_type=dbcsr_type_no_symmetry)
      DO ispin = 1, nspin
         ! calculate MAO overlap matrix
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_smm(1)%matrix, mao_coef(ispin)%matrix, &
                                0.0_dp, cgmat)
         CALL cp_dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, amat)
         ! calculate inverse of MAO overlap
         threshold = 1.e-8_dp
         CALL invert_Hotelling(tmat, amat, threshold, norm_convergence=1.e-4_dp, silent=.TRUE.)
         CALL cp_dbcsr_copy(mao_smat(ispin)%matrix, amat)
         ! calculate q-matrix q = C*Q*C
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_q(ispin)%matrix, mao_coef(ispin)%matrix, &
                                0.0_dp, cgmat, filter_eps=eps_filter)
         CALL cp_dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, &
                                0.0_dp, qmat, filter_eps=eps_filter)
         CALL cp_dbcsr_copy(mao_qmat(ispin)%matrix, qmat)
         ! calculate density matrix
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, qmat, tmat, 0.0_dp, axmat, filter_eps=eps_filter)
         CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmat, axmat, 0.0_dp, mao_dmat(ispin)%matrix, &
                                filter_eps=eps_filter)
      END DO
      CALL cp_dbcsr_release(amat)
      CALL cp_dbcsr_release(tmat)
      CALL cp_dbcsr_release(qmat)
      CALL cp_dbcsr_release(cgmat)
      CALL cp_dbcsr_release(axmat)

      ! calculate unassigned charge : n - Tr PS
      DO ispin = 1, nspin
         CALL cp_dbcsr_trace(mao_dmat(ispin)%matrix, mao_smat(ispin)%matrix, ua_charge(ispin))
         ua_charge(ispin) = electra(ispin)-ua_charge(ispin)
      END DO
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         DO ispin = 1, nspin
            WRITE (UNIT=unit_nr, FMT="(T2,A,T32,A,i2,T55,A,F12.8)") &
               "Unassigned charge", "Spin ", ispin, "delta charge =", ua_charge(ispin)
         END DO
      END IF

      ! occupation numbers: single atoms
      ! We use S_A = 1
      ! At the gamma point we use an effective MIC
      CALL get_qs_env(qs_env, natom=natom)
      ALLOCATE (occnumA(natom, nspin))
      occnumA = 0.0_dp
      DO ispin = 1, nspin
         DO iatom = 1, natom
            CALL cp_dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, &
                                      row=iatom, col=iatom, block=block, found=found)
            IF (found) THEN
               DO iab = 1, SIZE(block, 1)
                  occnumA(iatom, ispin) = occnumA(iatom, ispin)+block(iab, iab)
               END DO
            END IF
         END DO
      END DO
      CALL mp_sum(occnumA, para_env%group)

      ! occupation numbers: atom pairs
      ALLOCATE (occnumAB(natom, natom, nspin))
      occnumAB = 0.0_dp
      DO ispin = 1, nspin
         CALL cp_dbcsr_init(qmat_diag)
         CALL cp_dbcsr_create(qmat_diag, name="MAO diagonal density", template=mao_dmat(1)%matrix)
         CALL cp_dbcsr_init(smat_diag)
         CALL cp_dbcsr_create(smat_diag, name="MAO diagonal overlap", template=mao_dmat(1)%matrix)
         ! replicate the diagonal blocks of the density and overlap matrices
         CALL cp_dbcsr_get_block_diag(mao_qmat(ispin)%matrix, qmat_diag)
         CALL cp_dbcsr_replicate_all(qmat_diag)
         CALL cp_dbcsr_get_block_diag(mao_smat(ispin)%matrix, smat_diag)
         CALL cp_dbcsr_replicate_all(smat_diag)
         DO ia = 1, natom
            DO ib = ia+1, natom
               iab = 0
               CALL cp_dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, &
                                         row=ia, col=ib, block=block, found=found)
               IF (found) iab = 1
               CALL mp_sum(iab, para_env%group)
               CPASSERT(iab <= 1)
               IF (iab == 0 .AND. para_env%ionode) THEN
                  ! AB block is not available N_AB = N_A + N_B
                  ! Do this only on the "source" processor
                  occnumAB(ia, ib, ispin) = occnumA(ia, ispin)+occnumA(ib, ispin)
                  occnumAB(ib, ia, ispin) = occnumA(ia, ispin)+occnumA(ib, ispin)
               ELSE IF (found) THEN
                  ! owner of AB block performs calculation
                  na = SIZE(block, 1)
                  nb = SIZE(block, 2)
                  nab = na+nb
                  ALLOCATE (sab(nab, nab), qab(nab, nab), sinv(nab, nab))
                  ! qmat
                  qab(1:na, na+1:nab) = block(1:na, 1:nb)
                  qab(na+1:nab, 1:na) = TRANSPOSE(block(1:na, 1:nb))
                  CALL cp_dbcsr_get_block_p(matrix=qmat_diag, row=ia, col=ia, block=diag, found=fo)
                  CPASSERT(fo)
                  qab(1:na, 1:na) = diag(1:na, 1:na)
                  CALL cp_dbcsr_get_block_p(matrix=qmat_diag, row=ib, col=ib, block=diag, found=fo)
                  CPASSERT(fo)
                  qab(na+1:nab, na+1:nab) = diag(1:nb, 1:nb)
                  ! smat
                  CALL cp_dbcsr_get_block_p(matrix=mao_smat(ispin)%matrix, &
                                            row=ia, col=ib, block=block, found=fo)
                  CPASSERT(fo)
                  sab(1:na, na+1:nab) = block(1:na, 1:nb)
                  sab(na+1:nab, 1:na) = TRANSPOSE(block(1:na, 1:nb))
                  CALL cp_dbcsr_get_block_p(matrix=smat_diag, row=ia, col=ia, block=diag, found=fo)
                  CPASSERT(fo)
                  sab(1:na, 1:na) = diag(1:na, 1:na)
                  CALL cp_dbcsr_get_block_p(matrix=smat_diag, row=ib, col=ib, block=diag, found=fo)
                  CPASSERT(fo)
                  sab(na+1:nab, na+1:nab) = diag(1:nb, 1:nb)
                  ! inv smat
                  sinv(1:nab, 1:nab) = sab(1:nab, 1:nab)
                  CALL invmat_symm(sinv)
                  ! Tr(Q*Sinv)
                  occnumAB(ia, ib, ispin) = SUM(qab*sinv)
                  occnumAB(ib, ia, ispin) = occnumAB(ia, ib, ispin)
                  !
                  DEALLOCATE (sab, qab, sinv)
               END IF
            END DO
         END DO
         CALL cp_dbcsr_release(qmat_diag)
         CALL cp_dbcsr_release(smat_diag)
      END DO
      CALL mp_sum(occnumAB, para_env%group)

      ! calculate shared electron numbers (AB)
      ALLOCATE (selnAB(natom, natom, nspin))
      selnAB = 0.0_dp
      DO ispin = 1, nspin
         DO ia = 1, natom
            DO ib = ia+1, natom
               selnAB(ia, ib, ispin) = occnumA(ia, ispin)+occnumA(ib, ispin)-occnumAB(ia, ib, ispin)
               selnAB(ib, ia, ispin) = selnAB(ia, ib, ispin)
            END DO
         END DO
      END DO

      IF (.NOT. neglect_abc) THEN
         ! calculate N_ABC
         nabc = (natom*(natom-1)*(natom-2))/6
         ALLOCATE (occnumABC(nabc, nspin))
         occnumABC = -1.0_dp
         DO ispin = 1, nspin
            CALL cp_dbcsr_init(qmat_diag)
            CALL cp_dbcsr_create(qmat_diag, name="MAO diagonal density", template=mao_dmat(1)%matrix)
            CALL cp_dbcsr_init(smat_diag)
            CALL cp_dbcsr_create(smat_diag, name="MAO diagonal overlap", template=mao_dmat(1)%matrix)
            ! replicate the diagonal blocks of the density and overlap matrices
            CALL cp_dbcsr_get_block_diag(mao_qmat(ispin)%matrix, qmat_diag)
            CALL cp_dbcsr_replicate_all(qmat_diag)
            CALL cp_dbcsr_get_block_diag(mao_smat(ispin)%matrix, smat_diag)
            CALL cp_dbcsr_replicate_all(smat_diag)
            iabc = 0
            DO ia = 1, natom
               CALL cp_dbcsr_get_block_p(matrix=qmat_diag, row=ia, col=ia, block=qblka, found=fo)
               CPASSERT(fo)
               CALL cp_dbcsr_get_block_p(matrix=smat_diag, row=ia, col=ia, block=sblka, found=fo)
               CPASSERT(fo)
               na = SIZE(qblka, 1)
               DO ib = ia+1, natom
                  ! screen with SEN(AB)
                  IF (selnAB(ia, ib, ispin) < eps_abc) THEN
                     iabc = iabc+(natom-ib)
                     CYCLE
                  END IF
                  CALL cp_dbcsr_get_block_p(matrix=qmat_diag, row=ib, col=ib, block=qblkb, found=fo)
                  CPASSERT(fo)
                  CALL cp_dbcsr_get_block_p(matrix=smat_diag, row=ib, col=ib, block=sblkb, found=fo)
                  CPASSERT(fo)
                  nb = SIZE(qblkb, 1)
                  nab = na+nb
                  ALLOCATE (qmatab(na, nb), smatab(na, nb))
                  CALL cp_dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, row=ia, col=ib, &
                                            block=block, found=found)
                  qmatab = 0.0_dp
                  IF (found) qmatab(1:na, 1:nb) = block(1:na, 1:nb)
                  CALL mp_sum(qmatab, para_env%group)
                  CALL cp_dbcsr_get_block_p(matrix=mao_smat(ispin)%matrix, row=ia, col=ib, &
                                            block=block, found=found)
                  smatab = 0.0_dp
                  IF (found) smatab(1:na, 1:nb) = block(1:na, 1:nb)
                  CALL mp_sum(smatab, para_env%group)
                  DO ic = ib+1, natom
                     ! screen with SEN(AB)
                     IF ((selnAB(ia, ic, ispin) < eps_abc) .OR. (selnAB(ib, ic, ispin) < eps_abc)) THEN
                        iabc = iabc+1
                        CYCLE
                     END IF
                     CALL cp_dbcsr_get_block_p(matrix=qmat_diag, row=ic, col=ic, block=qblkc, found=fo)
                     CPASSERT(fo)
                     CALL cp_dbcsr_get_block_p(matrix=smat_diag, row=ic, col=ic, block=sblkc, found=fo)
                     CPASSERT(fo)
                     nc = SIZE(qblkc, 1)
                     ALLOCATE (qmatac(na, nc), smatac(na, nc))
                     CALL cp_dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, row=ia, col=ic, &
                                               block=block, found=found)
                     qmatac = 0.0_dp
                     IF (found) qmatac(1:na, 1:nc) = block(1:na, 1:nc)
                     CALL mp_sum(qmatac, para_env%group)
                     CALL cp_dbcsr_get_block_p(matrix=mao_smat(ispin)%matrix, row=ia, col=ic, &
                                               block=block, found=found)
                     smatac = 0.0_dp
                     IF (found) smatac(1:na, 1:nc) = block(1:na, 1:nc)
                     CALL mp_sum(smatac, para_env%group)
                     ALLOCATE (qmatbc(nb, nc), smatbc(nb, nc))
                     CALL cp_dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, row=ib, col=ic, &
                                               block=block, found=found)
                     qmatbc = 0.0_dp
                     IF (found) qmatbc(1:nb, 1:nc) = block(1:nb, 1:nc)
                     CALL mp_sum(qmatbc, para_env%group)
                     CALL cp_dbcsr_get_block_p(matrix=mao_smat(ispin)%matrix, row=ib, col=ic, &
                                               block=block, found=found)
                     smatbc = 0.0_dp
                     IF (found) smatbc(1:nb, 1:nc) = block(1:nb, 1:nc)
                     CALL mp_sum(smatbc, para_env%group)
                     !
                     nabc = na+nb+nc
                     ALLOCATE (sab(nabc, nabc), sinv(nabc, nabc), qab(nabc, nabc))
                     !
                     qab(1:na, 1:na) = qblka(1:na, 1:na)
                     qab(na+1:nab, na+1:nab) = qblkb(1:nb, 1:nb)
                     qab(nab+1:nabc, nab+1:nabc) = qblkc(1:nc, 1:nc)
                     qab(1:na, na+1:nab) = qmatab(1:na, 1:nb)
                     qab(na+1:nab, 1:na) = TRANSPOSE(qmatab(1:na, 1:nb))
                     qab(1:na, nab+1:nabc) = qmatac(1:na, 1:nc)
                     qab(nab+1:nabc, 1:na) = TRANSPOSE(qmatac(1:na, 1:nc))
                     qab(na+1:nab, nab+1:nabc) = qmatbc(1:nb, 1:nc)
                     qab(nab+1:nabc, na+1:nab) = TRANSPOSE(qmatbc(1:nb, 1:nc))
                     !
                     sab(1:na, 1:na) = sblka(1:na, 1:na)
                     sab(na+1:nab, na+1:nab) = sblkb(1:nb, 1:nb)
                     sab(nab+1:nabc, nab+1:nabc) = sblkc(1:nc, 1:nc)
                     sab(1:na, na+1:nab) = smatab(1:na, 1:nb)
                     sab(na+1:nab, 1:na) = TRANSPOSE(smatab(1:na, 1:nb))
                     sab(1:na, nab+1:nabc) = smatac(1:na, 1:nc)
                     sab(nab+1:nabc, 1:na) = TRANSPOSE(smatac(1:na, 1:nc))
                     sab(na+1:nab, nab+1:nabc) = smatbc(1:nb, 1:nc)
                     sab(nab+1:nabc, na+1:nab) = TRANSPOSE(smatbc(1:nb, 1:nc))
                     ! inv smat
                     sinv(1:nabc, 1:nabc) = sab(1:nabc, 1:nabc)
                     CALL invmat_symm(sinv)
                     ! Tr(Q*Sinv)
                     iabc = iabc+1
                     me = MOD(iabc, para_env%num_pe)
                     IF (me == para_env%mepos) THEN
                        occnumABC(iabc, ispin) = SUM(qab*sinv)
                     ELSE
                        occnumABC(iabc, ispin) = 0.0_dp
                     END IF
                     !
                     DEALLOCATE (sab, sinv, qab)
                     DEALLOCATE (qmatac, smatac)
                     DEALLOCATE (qmatbc, smatbc)
                  END DO
                  DEALLOCATE (qmatab, smatab)
               END DO
            END DO
            CALL cp_dbcsr_release(qmat_diag)
            CALL cp_dbcsr_release(smat_diag)
         END DO
         CALL mp_sum(occnumABC, para_env%group)
      END IF

      IF (.NOT. neglect_abc) THEN
         ! calculate shared electron numbers (ABC)
         nabc = (natom*(natom-1)*(natom-2))/6
         ALLOCATE (selnABC(nabc, nspin))
         selnABC = 0.0_dp
         DO ispin = 1, nspin
            iabc = 0
            DO ia = 1, natom
               DO ib = ia+1, natom
                  DO ic = ib+1, natom
                     iabc = iabc+1
                     IF (occnumABC(iabc, ispin) >= 0.0_dp) THEN
                        selnABC(iabc, ispin) = occnumA(ia, ispin)+occnumA(ib, ispin)+occnumA(ic, ispin)- &
                                               occnumAB(ia, ib, ispin)-occnumAB(ia, ic, ispin)-occnumAB(ib, ic, ispin)+ &
                                               occnumABC(iabc, ispin)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END IF

      ! calculate atomic charge
      ALLOCATE (raq(natom, nspin))
      raq = 0.0_dp
      DO ispin = 1, nspin
         DO ia = 1, natom
            raq(ia, ispin) = occnumA(ia, ispin)
            DO ib = 1, natom
               raq(ia, ispin) = raq(ia, ispin)-0.5_dp*selnAB(ia, ib, ispin)
            END DO
         END DO
         IF (.NOT. neglect_abc) THEN
            iabc = 0
            DO ia = 1, natom
               DO ib = ia+1, natom
                  DO ic = ib+1, natom
                     iabc = iabc+1
                     raq(ia, ispin) = raq(ia, ispin)+selnABC(iabc, ispin)/3._dp
                     raq(ib, ispin) = raq(ib, ispin)+selnABC(iabc, ispin)/3._dp
                     raq(ic, ispin) = raq(ic, ispin)+selnABC(iabc, ispin)/3._dp
                  END DO
               END DO
            END DO
         END IF
      END DO

      ! calculate unassigned charge (from sum over atomic charges)
      DO ispin = 1, nspin
         deltaq = (electra(ispin)-SUM(raq(1:natom, ispin)))-ua_charge(ispin)
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A,T32,A,i2,T55,A,F12.8)") &
               "Cutoff error on charge", "Spin ", ispin, "error charge =", deltaq
         END IF
      END DO

      ! analyze unassigned charge
      ALLOCATE (uaq(natom, nspin))
      uaq = 0.0_dp
      IF (analyze_ua) THEN
         CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
         CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb, sab_all=sab_all)
         CALL cp_dbcsr_get_info(mao_coef(1)%matrix, row_blk_size=mao_blk_sizes, &
                                col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
         CALL cp_dbcsr_get_info(matrix_s(1, 1)%matrix, row_blk_size=row_blk_sizes)
         CALL cp_dbcsr_init(amat)
         CALL cp_dbcsr_create(amat, name="temp", template=matrix_s(1, 1)%matrix)
         CALL cp_dbcsr_init(tmat)
         CALL cp_dbcsr_create(tmat, name="temp", template=mao_coef(1)%matrix)
         ! replicate diagonal of smm matrix
         CALL cp_dbcsr_init(smat_diag)
         CALL cp_dbcsr_get_block_diag(matrix_smm(1)%matrix, smat_diag)
         CALL cp_dbcsr_replicate_all(smat_diag)

         ALLOCATE (orb_blk(natom), mao_blk(natom))
         DO ia = 1, natom
            orb_blk = row_blk_sizes
            mao_blk = row_blk_sizes
            mao_blk(ia) = col_blk_sizes(ia)
            CALL cp_dbcsr_init(sumat)
            CALL cp_dbcsr_create(sumat, name="Smat", dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                                 row_blk_size=mao_blk, col_blk_size=mao_blk, nze=0)
            CALL cp_dbcsr_alloc_block_from_nbl(sumat, sab_orb)
            CALL cp_dbcsr_init(cholmat)
            CALL cp_dbcsr_create(cholmat, name="Cholesky matrix", dist=dbcsr_dist, &
                                 matrix_type=dbcsr_type_no_symmetry, row_blk_size=mao_blk, col_blk_size=mao_blk, nze=0)
            CALL cp_dbcsr_init(rumat)
            CALL cp_dbcsr_create(rumat, name="Rmat", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                                 row_blk_size=orb_blk, col_blk_size=mao_blk, nze=0)
            CALL cp_dbcsr_alloc_block_from_nbl(rumat, sab_orb, .TRUE.)
            CALL cp_dbcsr_init(crumat)
            CALL cp_dbcsr_create(crumat, name="Rmat*Umat", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                                 row_blk_size=orb_blk, col_blk_size=mao_blk, nze=0)
            ! replicate row and col of smo matrix
            ALLOCATE (rowblock(natom))
            DO ib = 1, natom
               na = mao_blk_sizes(ia)
               nb = row_blk_sizes(ib)
               ALLOCATE (rowblock(ib)%mat(na, nb))
               rowblock(ib)%mat = 0.0_dp
               CALL cp_dbcsr_get_block_p(matrix=matrix_smo(1)%matrix, row=ia, col=ib, &
                                         block=block, found=found)
               IF (found) rowblock(ib)%mat(1:na, 1:nb) = block(1:na, 1:nb)
               CALL mp_sum(rowblock(ib)%mat, para_env%group)
            END DO
            !
            DO ispin = 1, nspin
               CALL cp_dbcsr_copy(tmat, mao_coef(ispin)%matrix)
               CALL cp_dbcsr_replicate_all(tmat)
               CALL cp_dbcsr_iterator_start(dbcsr_iter, matrix_s(1, 1)%matrix)
               DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
                  CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, block)
                  CALL cp_dbcsr_get_block_p(matrix=sumat, row=iatom, col=jatom, block=sblk, found=fos)
                  CPASSERT(fos)
                  CALL cp_dbcsr_get_block_p(matrix=rumat, row=iatom, col=jatom, block=rblku, found=for)
                  CPASSERT(for)
                  CALL cp_dbcsr_get_block_p(matrix=rumat, row=jatom, col=iatom, block=rblkl, found=for)
                  CPASSERT(for)
                  CALL cp_dbcsr_get_block_p(matrix=tmat, row=ia, col=ia, block=cmao, found=found)
                  CPASSERT(found)
                  IF (iatom /= ia .AND. jatom /= ia) THEN
                     ! copy original overlap matrix
                     sblk = block
                     rblku = block
                     rblkl = TRANSPOSE(block)
                  ELSE IF (iatom /= ia) THEN
                     rblkl = TRANSPOSE(block)
                     sblk = MATMUL(TRANSPOSE(rowblock(iatom)%mat), cmao)
                     rblku = sblk
                  ELSE IF (jatom /= ia) THEN
                     rblku = block
                     sblk = MATMUL(TRANSPOSE(cmao), rowblock(jatom)%mat)
                     rblkl = TRANSPOSE(sblk)
                  ELSE
                     CALL cp_dbcsr_get_block_p(matrix=smat_diag, row=ia, col=ia, block=block, found=found)
                     CPASSERT(found)
                     sblk = MATMUL(TRANSPOSE(cmao), MATMUL(block, cmao))
                     rblku = MATMUL(TRANSPOSE(rowblock(ia)%mat), cmao)
                  END IF
               END DO
               CALL cp_dbcsr_iterator_stop(dbcsr_iter)
               ! Cholesky decomposition of SUMAT = U'U
               CALL cp_dbcsr_desymmetrize(sumat, cholmat)
               CALL cp_dbcsr_cholesky_decompose(cholmat, para_env=para_env, blacs_env=blacs_env)
               ! T = R*inv(U)
               ssize = SUM(mao_blk)
               CALL cp_dbcsr_cholesky_restore(rumat, ssize, cholmat, crumat, op="SOLVE", pos="RIGHT", &
                                              transa="N", para_env=para_env, blacs_env=blacs_env)
               ! A = T*transpose(T)
               CALL cp_dbcsr_multiply("N", "T", 1.0_dp, crumat, crumat, 0.0_dp, amat, &
                                      filter_eps=eps_filter)
               ! Tr(P*A)
               CALL cp_dbcsr_trace(matrix_p(ispin, 1)%matrix, amat, uaq(ia, ispin))
               uaq(ia, ispin) = uaq(ia, ispin)-electra(ispin)
            END DO
            !
            CALL cp_dbcsr_release(sumat)
            CALL cp_dbcsr_release(cholmat)
            CALL cp_dbcsr_release(rumat)
            CALL cp_dbcsr_release(crumat)
            !
            DO ib = 1, natom
               DEALLOCATE (rowblock(ib)%mat)
            END DO
            DEALLOCATE (rowblock)
         END DO
         CALL cp_dbcsr_release(smat_diag)
         CALL cp_dbcsr_release(amat)
         CALL cp_dbcsr_release(tmat)
         DEALLOCATE (orb_blk, mao_blk)
      END IF
      !
      raq(1:natom, 1:nspin) = raq(1:natom, 1:nspin)-uaq(1:natom, 1:nspin)
      DO ispin = 1, nspin
         deltaq = electra(ispin)-SUM(raq(1:natom, ispin))
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A,T32,A,i2,T55,A,F12.8)") &
               "Charge/Atom redistributed", "Spin ", ispin, "delta charge =", &
               (deltaq+ua_charge(ispin))/REAL(natom, KIND=dp)
         END IF
      END DO

      ! output charges
      IF (unit_nr > 0) THEN
         IF (nspin == 1) THEN
            WRITE (unit_nr, "(/,T2,A,T40,A,T75,A)") "MAO atomic charges ", "Atom", "Charge"
         ELSE
            WRITE (unit_nr, "(/,T2,A,T40,A,T55,A,T70,A)") "MAO atomic charges ", "Atom", "Charge", "Spin Charge"
         END IF
         DO ispin = 1, nspin
            deltaq = electra(ispin)-SUM(raq(1:natom, ispin))
            raq(:, ispin) = raq(:, ispin)+deltaq/REAL(natom, KIND=dp)
         END DO
         total_charge = 0.0_dp
         total_spin = 0.0_dp
         DO iatom = 1, natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 element_symbol=element_symbol, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
            IF (nspin == 1) THEN
               WRITE (unit_nr, "(T30,I6,T42,A2,T69,F12.6)") iatom, element_symbol, zeff-raq(iatom, 1)
               total_charge = total_charge+(zeff-raq(iatom, 1))
            ELSE
               WRITE (unit_nr, "(T30,I6,T42,A2,T48,F12.6,T69,F12.6)") iatom, element_symbol, &
                  zeff-raq(iatom, 1)-raq(iatom, 2), raq(iatom, 1)-raq(iatom, 2)
               total_charge = total_charge+(zeff-raq(iatom, 1)-raq(iatom, 2))
               total_spin = total_spin+(raq(iatom, 1)-raq(iatom, 2))
            END IF
         END DO
         IF (nspin == 1) THEN
            WRITE (unit_nr, "(T2,A,T69,F12.6)") "Total Charge", total_charge
         ELSE
            WRITE (unit_nr, "(T2,A,T49,F12.6,T69,F12.6)") "Total Charge", total_charge, total_spin
         END IF
      END IF

      IF (analyze_ua) THEN
         ! output unassigned charges
         IF (unit_nr > 0) THEN
            IF (nspin == 1) THEN
               WRITE (unit_nr, "(/,T2,A,T40,A,T75,A)") "MAO hypervalent charges ", "Atom", "Charge"
            ELSE
               WRITE (unit_nr, "(/,T2,A,T40,A,T55,A,T70,A)") "MAO hypervalent charges ", "Atom", &
                  "Charge", "Spin Charge"
            END IF
            total_charge = 0.0_dp
            total_spin = 0.0_dp
            DO iatom = 1, natom
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                    element_symbol=element_symbol)
               IF (nspin == 1) THEN
                  WRITE (unit_nr, "(T30,I6,T42,A2,T69,F12.6)") iatom, element_symbol, uaq(iatom, 1)
                  total_charge = total_charge+uaq(iatom, 1)
               ELSE
                  WRITE (unit_nr, "(T30,I6,T42,A2,T48,F12.6,T69,F12.6)") iatom, element_symbol, &
                     uaq(iatom, 1)+uaq(iatom, 2), uaq(iatom, 1)-uaq(iatom, 2)
                  total_charge = total_charge+uaq(iatom, 1)+uaq(iatom, 2)
                  total_spin = total_spin+uaq(iatom, 1)-uaq(iatom, 2)
               END IF
            END DO
            IF (nspin == 1) THEN
               WRITE (unit_nr, "(T2,A,T69,F12.6)") "Total Charge", total_charge
            ELSE
               WRITE (unit_nr, "(T2,A,T49,F12.6,T69,F12.6)") "Total Charge", total_charge, total_spin
            END IF
         END IF
      END IF

      ! output shared electron numbers AB
      IF (unit_nr > 0) THEN
         IF (nspin == 1) THEN
            WRITE (unit_nr, "(/,T2,A,T31,A,T40,A,T78,A)") "Shared electron numbers ", "Atom", "Atom", "SEN"
         ELSE
            WRITE (unit_nr, "(/,T2,A,T31,A,T40,A,T51,A,T63,A,T71,A)") "Shared electron numbers ", "Atom", "Atom", &
               "SEN(1)", "SEN(2)", "SEN(total)"
         END IF
         DO ia = 1, natom
            DO ib = ia+1, natom
               CALL get_atomic_kind(atomic_kind=particle_set(ia)%atomic_kind, element_symbol=esa)
               CALL get_atomic_kind(atomic_kind=particle_set(ib)%atomic_kind, element_symbol=esb)
               IF (nspin == 1) THEN
                  IF (selnAB(ia, ib, 1) > eps_ab) THEN
                     WRITE (unit_nr, "(T26,I6,' ',A2,T35,I6,' ',A2,T69,F12.6)") ia, esa, ib, esb, selnAB(ia, ib, 1)
                  END IF
               ELSE
                  IF ((selnAB(ia, ib, 1)+selnAB(ia, ib, 2)) > eps_ab) THEN
                     WRITE (unit_nr, "(T26,I6,' ',A2,T35,I6,' ',A2,T45,3F12.6)") ia, esa, ib, esb, &
                        selnAB(ia, ib, 1), selnAB(ia, ib, 2), (selnAB(ia, ib, 1)+selnAB(ia, ib, 2))
                  END IF
               END IF
            END DO
         END DO
      END IF

      IF (.NOT. neglect_abc) THEN
         ! output shared electron numbers ABC
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, "(/,T2,A,T40,A,T49,A,T58,A,T78,A)") "Shared electron numbers ABC", &
               "Atom", "Atom", "Atom", "SEN"
            senmax = 0.0_dp
            iabc = 0
            DO ia = 1, natom
               DO ib = ia+1, natom
                  DO ic = ib+1, natom
                     iabc = iabc+1
                     senabc = SUM(selnABC(iabc, :))
                     senmax = MAX(senmax, senabc)
                     IF (senabc > eps_abc) THEN
                        CALL get_atomic_kind(atomic_kind=particle_set(ia)%atomic_kind, element_symbol=esa)
                        CALL get_atomic_kind(atomic_kind=particle_set(ib)%atomic_kind, element_symbol=esb)
                        CALL get_atomic_kind(atomic_kind=particle_set(ic)%atomic_kind, element_symbol=esc)
                        WRITE (unit_nr, "(T35,I6,' ',A2,T44,I6,' ',A2,T53,I6,' ',A2,T69,F12.6)") &
                           ia, esa, ib, esb, ic, esc, senabc
                     END IF
                  END DO
               END DO
            END DO
            WRITE (unit_nr, "(T2,A,T69,F12.6)") "Maximum SEN value calculated", senmax
         END IF
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') &
            '!---------------------------END OF MAO ANALYSIS-------------------------------!'
      END IF

      ! Deallocate temporary arrays
      DEALLOCATE (occnumA, occnumAB, selnAB, raq, uaq)
      IF (.NOT. neglect_abc) THEN
         DEALLOCATE (occnumABC, selnABC)
      END IF

      ! Deallocate the neighbor list structure
      IF (ASSOCIATED(smm_list)) THEN
         DO iab = 1, SIZE(smm_list)
            CALL deallocate_neighbor_list_set(smm_list(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (smm_list)
      END IF
      IF (ASSOCIATED(smo_list)) THEN
         DO iab = 1, SIZE(smo_list)
            CALL deallocate_neighbor_list_set(smo_list(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (smo_list)
      END IF

      DEALLOCATE (mao_basis_set_list, orb_basis_set_list)

      IF (ASSOCIATED(matrix_smm)) CALL cp_dbcsr_deallocate_matrix_set(matrix_smm)
      IF (ASSOCIATED(matrix_smo)) CALL cp_dbcsr_deallocate_matrix_set(matrix_smo)
      IF (ASSOCIATED(matrix_q)) CALL cp_dbcsr_deallocate_matrix_set(matrix_q)

      IF (ASSOCIATED(mao_coef)) CALL cp_dbcsr_deallocate_matrix_set(mao_coef)
      IF (ASSOCIATED(mao_grad)) CALL cp_dbcsr_deallocate_matrix_set(mao_grad)
      IF (ASSOCIATED(mao_dmat)) CALL cp_dbcsr_deallocate_matrix_set(mao_dmat)
      IF (ASSOCIATED(mao_smat)) CALL cp_dbcsr_deallocate_matrix_set(mao_smat)
      IF (ASSOCIATED(mao_qmat)) CALL cp_dbcsr_deallocate_matrix_set(mao_qmat)

      CALL timestop(handle)

   END SUBROUTINE mao_analysis

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mao_coef ...
!> \param ref_basis_set ...
!> \param pmat_external ...
!> \param smat_external ...
!> \param molecular ...
!> \param max_iter ...
!> \param eps_grad ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE mao_basis_simple(qs_env, mao_coef, ref_basis_set, pmat_external, smat_external, &
                               molecular, max_iter, eps_grad, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_coef
      CHARACTER(len=*), OPTIONAL                         :: ref_basis_set
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: pmat_external, smat_external
      LOGICAL, INTENT(IN), OPTIONAL                      :: molecular
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iter
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_grad
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_basis_simple', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=10)                                  :: mao_basis_set
      INTEGER                                            :: handle, iab, ikind, ispin, iter, iw, &
                                                            mao_max_iter, natom, nbas, nimages, &
                                                            nkind, nmao, nspin
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: molecule
      REAL(KIND=dp)                                      :: a1, a2, alpha, an, beta, electra(2), &
                                                            eps_filter, eps_fun, fa1, fa2, fnnew, &
                                                            fnold, fval, grad_norm, mao_eps_grad
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_grad, matrix_q, matrix_smm, &
                                                            matrix_smo
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: matrix_p, matrix_s
      TYPE(cp_dbcsr_type)                                :: amat, cgmat, tmat
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: mao_basis_set_list, orb_basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: smm_list, smo_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
         WRITE (UNIT=unit_nr, FMT="(T37,A)") "MAO BASIS"
         WRITE (unit_nr, '(T2,A)') '!-----------------------------------------------------------------------------!'
      END IF

      ! k-points?
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nimages = dft_control%nimages
      CPASSERT(nimages == 1)

      ! output
      IF (PRESENT(unit_nr)) THEN
         iw = unit_nr
      ELSE
         iw = -1
      END IF

      ! molecules
      IF (PRESENT(molecular)) THEN
         molecule = molecular
      ELSE
         molecule = .FALSE.
      END IF

      ! iterations
      IF (PRESENT(max_iter)) THEN
         mao_max_iter = max_iter
      ELSE
         mao_max_iter = 0
      END IF

      ! threshold
      IF (PRESENT(eps_grad)) THEN
         mao_eps_grad = eps_grad
      ELSE
         mao_eps_grad = 0.00001_dp
      END IF
      eps_fun = 10._dp*mao_eps_grad

      ! Reference basis set
      IF (PRESENT(ref_basis_set)) THEN
         mao_basis_set = ref_basis_set
      ELSE
         mao_basis_set = "ORB"
      END IF
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (mao_basis_set_list(nkind), orb_basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         NULLIFY (mao_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (orb_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (basis_set_a, basis_set_b)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type="ORB")
         IF (ASSOCIATED(basis_set_a)) orb_basis_set_list(ikind)%gto_basis_set => basis_set_a
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_b, basis_type=mao_basis_set)
         IF (ASSOCIATED(basis_set_b)) mao_basis_set_list(ikind)%gto_basis_set => basis_set_b
      END DO
      IF (iw > 0) THEN
         DO ikind = 1, nkind
            IF (.NOT. ASSOCIATED(mao_basis_set_list(ikind)%gto_basis_set)) THEN
               WRITE (UNIT=iw, FMT="(T2,A,I4)") &
                  "WARNING: No MAO basis set associated with Kind ", ikind
            ELSE
               CALL get_qs_kind(qs_kind_set(ikind), mao=nmao)
               nbas = mao_basis_set_list(ikind)%gto_basis_set%nsgf
               WRITE (UNIT=iw, FMT="(T2,A,I4,T30,A,I2,T56,A,I10)") &
                  "MAO basis set Kind ", ikind, " Number of MAO:", nmao, " Number of BSF:", nbas
            END IF
         END DO
      END IF

      ! neighbor lists
      NULLIFY (smm_list, smo_list)
      CALL setup_neighbor_list(smm_list, mao_basis_set_list, molecular=molecule, qs_env=qs_env)
      CALL setup_neighbor_list(smo_list, mao_basis_set_list, orb_basis_set_list, &
                               molecular=molecule, qs_env=qs_env)

      ! overlap matrices
      NULLIFY (matrix_smm, matrix_smo)
      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL build_overlap_matrix_simple(ks_env, matrix_smm, &
                                       mao_basis_set_list, mao_basis_set_list, smm_list)
      CALL build_overlap_matrix_simple(ks_env, matrix_smo, &
                                       mao_basis_set_list, orb_basis_set_list, smo_list)

      ! get reference density matrix and overlap matrix
      IF (PRESENT(pmat_external)) THEN
         matrix_p => pmat_external
      ELSE
         CALL get_qs_env(qs_env, rho=rho)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      END IF
      IF (PRESENT(smat_external)) THEN
         matrix_s => smat_external
      ELSE
         CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
      END IF

      nspin = SIZE(matrix_p, 1)
      eps_filter = 0.0_dp
      ! Q matrix
      CALL mao_build_q(matrix_q, matrix_p, matrix_s, matrix_smm, matrix_smo, smm_list, electra, eps_filter)

      ! MAO matrices
      CALL get_qs_env(qs_env=qs_env, natom=natom)
      CALL get_ks_env(ks_env=ks_env, particle_set=particle_set, dbcsr_dist=dbcsr_dist)
      NULLIFY (mao_coef)
      CALL cp_dbcsr_allocate_matrix_set(mao_coef, nspin)
      ALLOCATE (row_blk_sizes(natom), col_blk_sizes(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, &
                            basis=mao_basis_set_list)
      CALL get_particle_set(particle_set, qs_kind_set, nmao=col_blk_sizes)
      ! check if MAOs have been specified
      DO iab = 1, natom
         IF (col_blk_sizes(iab) < 0) &
            CPABORT("Number of MAOs has to be specified in KIND section for all elements")
      END DO
      DO ispin = 1, nspin
         ! coeficients
         ALLOCATE (mao_coef(ispin)%matrix)
         CALL cp_dbcsr_init(mao_coef(ispin)%matrix)
         CALL cp_dbcsr_create(matrix=mao_coef(ispin)%matrix, &
                              name="MAO_COEF", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL cp_dbcsr_reserve_diag_blocks(matrix=mao_coef(ispin)%matrix)
      END DO
      DEALLOCATE (row_blk_sizes, col_blk_sizes)

      ! initialize MAO coeficients from diagonal blocks of the Q matrix
      DO ispin = 1, nspin
         CALL mao_initialization(mao_coef(ispin)%matrix, &
                                 matrix_q(ispin)%matrix, matrix_smm(1)%matrix)
      END DO

      IF (mao_max_iter < 1) THEN
         ! projection only
         CALL cp_dbcsr_get_info(mao_coef(1)%matrix, col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
         CALL cp_dbcsr_init(tmat)
         CALL cp_dbcsr_create(tmat, name="Binv", dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                              row_blk_size=col_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         DO ispin = 1, nspin
            CALL mao_function(mao_coef(ispin)%matrix, fval, matrix_q(ispin)%matrix, &
                              matrix_smm(1)%matrix, tmat, .FALSE.)
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(T2,A,T31,A,I2,T57,A,F12.8)") &
                  "MAO Projection", "Spin:", ispin, "Completness:", fval/electra(ispin)
            END IF
         END DO
         CALL cp_dbcsr_release(tmat)
      ELSE
         ! optimize MAOs
         NULLIFY (mao_grad)
         CALL cp_dbcsr_allocate_matrix_set(mao_grad, nspin)
         DO ispin = 1, nspin
            ALLOCATE (mao_grad(ispin)%matrix)
            CALL cp_dbcsr_init(mao_grad(ispin)%matrix)
            CALL cp_dbcsr_create(matrix=mao_grad(ispin)%matrix, name="MAO_GRAD", template=mao_coef(1)%matrix)
            CALL cp_dbcsr_reserve_diag_blocks(matrix=mao_grad(ispin)%matrix)
         END DO
         alpha = 0.25_dp
         beta = 0.0_dp
         CALL cp_dbcsr_get_info(mao_coef(1)%matrix, col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
         CALL cp_dbcsr_init(tmat)
         CALL cp_dbcsr_create(tmat, name="Binv", dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                              row_blk_size=col_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL cp_dbcsr_init(cgmat)
         CALL cp_dbcsr_create(cgmat, template=mao_grad(1)%matrix)
         CALL cp_dbcsr_init(amat)
         CALL cp_dbcsr_create(amat, template=mao_coef(1)%matrix)
         DO ispin = 1, nspin
            CALL mao_function_gradient(mao_coef(ispin)%matrix, fval, mao_grad(ispin)%matrix, &
                                       matrix_q(ispin)%matrix, matrix_smm(1)%matrix, tmat, .FALSE.)
            CALL cp_dbcsr_copy(cgmat, mao_grad(ispin)%matrix)
            CALL cp_dbcsr_norm(mao_grad(ispin)%matrix, dbcsr_norm_maxabsnorm, norm_scalar=grad_norm)
            fnold = mao_scalar_product(mao_grad(ispin)%matrix, mao_grad(ispin)%matrix)
            IF (iw > 0) THEN
               WRITE (UNIT=iw, FMT="(/,T2,A,T73,A,I2)") "MAO OPTIMIZATION", "Spin =", ispin
               WRITE (UNIT=iw, FMT="(T2,A,T24,A,F11.8,T48,A,F11.8,T69,A,F6.3)") &
                  "Initialization", "fval =", fval/electra(ispin), "grad =", grad_norm, "step =", alpha
            END IF
            DO iter = 1, mao_max_iter
               IF (grad_norm < mao_eps_grad) EXIT
               IF ((1.0_dp-fval/electra(ispin)) < eps_fun) EXIT
               CALL cp_dbcsr_add(mao_coef(ispin)%matrix, cgmat, 1.0_dp, alpha)
               CALL mao_orthogonalization(mao_coef(ispin)%matrix, matrix_smm(1)%matrix)
               CALL mao_function_gradient(mao_coef(ispin)%matrix, fval, mao_grad(ispin)%matrix, &
                                          matrix_q(ispin)%matrix, matrix_smm(1)%matrix, tmat, .TRUE.)
               CALL cp_dbcsr_norm(mao_grad(ispin)%matrix, dbcsr_norm_maxabsnorm, norm_scalar=grad_norm)
               IF (iw > 0) THEN
                  WRITE (UNIT=iw, FMT="(T2,A,i8,T24,A,F11.8,T48,A,F11.8,T69,A,F6.3)") &
                     "iter=", iter, "fval =", fval/electra(ispin), "grad =", grad_norm, "step =", alpha
               END IF
               fnnew = mao_scalar_product(mao_grad(ispin)%matrix, mao_grad(ispin)%matrix)
               beta = fnnew/fnold
               CALL cp_dbcsr_add(cgmat, mao_grad(ispin)%matrix, beta, 1.0_dp)
               fnold = fnnew
               ! line search, update alpha
               CALL cp_dbcsr_copy(amat, mao_coef(ispin)%matrix)
               CALL cp_dbcsr_add(amat, cgmat, 1.0_dp, 0.5_dp*alpha)
               CALL mao_orthogonalization(amat, matrix_smm(1)%matrix)
               CALL mao_function(amat, fa1, matrix_q(ispin)%matrix, matrix_smm(1)%matrix, tmat, .TRUE.)
               CALL cp_dbcsr_copy(amat, mao_coef(ispin)%matrix)
               CALL cp_dbcsr_add(amat, cgmat, 1.0_dp, alpha)
               CALL mao_orthogonalization(amat, matrix_smm(1)%matrix)
               CALL mao_function(amat, fa2, matrix_q(ispin)%matrix, matrix_smm(1)%matrix, tmat, .TRUE.)
               a2 = (4._dp*fa1-fa2-3._dp*fval)/alpha
               a1 = (fa2-fval-a2*alpha)/(alpha*alpha)
               an = -a2/(2._dp*a1)
               an = MIN(an, 2.0_dp*alpha)
               IF (an < 0.05_dp .OR. a1 > 0.0_dp) THEN
                  CALL cp_dbcsr_copy(cgmat, mao_grad(ispin)%matrix)
                  alpha = 0.25_dp
               ELSE
                  alpha = an
               END IF
            END DO
         END DO
         CALL cp_dbcsr_release(tmat)
         CALL cp_dbcsr_release(cgmat)
         CALL cp_dbcsr_release(amat)
         IF (ASSOCIATED(mao_grad)) CALL cp_dbcsr_deallocate_matrix_set(mao_grad)
      END IF

      ! Deallocate the neighbor list structure
      IF (ASSOCIATED(smm_list)) THEN
         DO iab = 1, SIZE(smm_list)
            CALL deallocate_neighbor_list_set(smm_list(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (smm_list)
      END IF
      IF (ASSOCIATED(smo_list)) THEN
         DO iab = 1, SIZE(smo_list)
            CALL deallocate_neighbor_list_set(smo_list(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (smo_list)
      END IF

      DEALLOCATE (mao_basis_set_list, orb_basis_set_list)

      IF (ASSOCIATED(matrix_smm)) CALL cp_dbcsr_deallocate_matrix_set(matrix_smm)
      IF (ASSOCIATED(matrix_smo)) CALL cp_dbcsr_deallocate_matrix_set(matrix_smo)
      IF (ASSOCIATED(matrix_q)) CALL cp_dbcsr_deallocate_matrix_set(matrix_q)

      CALL timestop(handle)

   END SUBROUTINE mao_basis_simple

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param pmat ...
!> \param smat ...
! **************************************************************************************************
   SUBROUTINE mao_initialization(mao_coef, pmat, smat)
      TYPE(cp_dbcsr_type)                                :: mao_coef, pmat, smat

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_initialization', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, info, jatom, lwork, m, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, bmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cblock, pblock, sblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, mao_coef)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, cblock)
         CPASSERT(iatom == jatom)
         m = SIZE(cblock, 2)
         NULLIFY (pblock, sblock)
         CALL cp_dbcsr_get_block_p(matrix=pmat, row=iatom, col=jatom, block=pblock, found=found)
         CPASSERT(found)
         CALL cp_dbcsr_get_block_p(matrix=smat, row=iatom, col=jatom, block=sblock, found=found)
         CPASSERT(found)
         n = SIZE(sblock, 1)
         lwork = MAX(n*n, 100)
         ALLOCATE (amat(n, n), bmat(n, n), w(n), work(lwork))
         amat(1:n, 1:n) = pblock(1:n, 1:n)
         bmat(1:n, 1:n) = sblock(1:n, 1:n)
         info = 0
         CALL lapack_ssygv(1, "V", "U", n, amat, n, bmat, n, w, work, lwork, info)
         CPASSERT(info == 0)
         cblock(1:n, 1:m) = amat(1:n, n:n-m+1:-1)
         DEALLOCATE (amat, bmat, w, work)
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE mao_initialization

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param fval ...
!> \param qmat ...
!> \param smat ...
!> \param binv ...
!> \param reuse ...
! **************************************************************************************************
   SUBROUTINE mao_function(mao_coef, fval, qmat, smat, binv, reuse)
      TYPE(cp_dbcsr_type)                                :: mao_coef
      REAL(KIND=dp), INTENT(OUT)                         :: fval
      TYPE(cp_dbcsr_type)                                :: qmat, smat, binv
      LOGICAL, INTENT(IN)                                :: reuse

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_function', routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: convergence, threshold
      TYPE(cp_dbcsr_type)                                :: bmat, scmat, tmat

      threshold = 1.e-8_dp
      convergence = 1.e-6_dp
      ! temp matrices
      CALL cp_dbcsr_init(scmat)
      CALL cp_dbcsr_create(scmat, template=mao_coef)
      CALL cp_dbcsr_init(bmat)
      CALL cp_dbcsr_create(bmat, template=binv)
      CALL cp_dbcsr_init(tmat)
      CALL cp_dbcsr_create(tmat, template=qmat)
      ! calculate B=C(T)*S*C matrix, S=(MAO,MAO) overlap
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, smat, mao_coef, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, mao_coef, scmat, 0.0_dp, bmat)
      ! calculate inverse of B
      CALL invert_Hotelling(binv, bmat, threshold, use_inv_as_guess=reuse, &
                            norm_convergence=convergence, silent=.TRUE.)
      ! calculate Binv*C and T=C(T)*Binv*C
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mao_coef, binv, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, scmat, mao_coef, 0.0_dp, tmat)
      ! function = Tr(Q*T)
      CALL cp_dbcsr_trace(qmat, tmat, fval)
      ! free temp matrices
      CALL cp_dbcsr_release(scmat)
      CALL cp_dbcsr_release(bmat)
      CALL cp_dbcsr_release(tmat)

   END SUBROUTINE mao_function

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param fval ...
!> \param mao_grad ...
!> \param qmat ...
!> \param smat ...
!> \param binv ...
!> \param reuse ...
! **************************************************************************************************
   SUBROUTINE mao_function_gradient(mao_coef, fval, mao_grad, qmat, smat, binv, reuse)
      TYPE(cp_dbcsr_type)                                :: mao_coef
      REAL(KIND=dp), INTENT(OUT)                         :: fval
      TYPE(cp_dbcsr_type)                                :: mao_grad, qmat, smat, binv
      LOGICAL, INTENT(IN)                                :: reuse

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_function_gradient', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: convergence, threshold
      TYPE(cp_dbcsr_type)                                :: bmat, scmat, t2mat, tmat

      threshold = 1.e-8_dp
      convergence = 1.e-6_dp
      ! temp matrices
      CALL cp_dbcsr_init(scmat)
      CALL cp_dbcsr_create(scmat, template=mao_coef)
      CALL cp_dbcsr_init(bmat)
      CALL cp_dbcsr_create(bmat, template=binv)
      CALL cp_dbcsr_init(tmat)
      CALL cp_dbcsr_create(tmat, template=qmat)
      CALL cp_dbcsr_init(t2mat)
      CALL cp_dbcsr_create(t2mat, template=scmat)
      ! calculate B=C(T)*S*C matrix, S=(MAO,MAO) overlap
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, smat, mao_coef, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("T", "N", 1.0_dp, mao_coef, scmat, 0.0_dp, bmat)
      ! calculate inverse of B
      CALL invert_Hotelling(binv, bmat, threshold, use_inv_as_guess=reuse, &
                            norm_convergence=convergence, silent=.TRUE.)
      ! calculate R=C*Binv and T=C*Binv*C(T)=R*C(T)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mao_coef, binv, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("N", "T", 1.0_dp, scmat, mao_coef, 0.0_dp, tmat)
      ! function = Tr(Q*T)
      CALL cp_dbcsr_trace(qmat, tmat, fval)
      ! Gradient part 1: g = 2*Q*C*Binv = 2*Q*R
      CALL cp_dbcsr_multiply("N", "N", 2.0_dp, qmat, scmat, 0.0_dp, mao_grad, &
                             retain_sparsity=.TRUE.)
      ! Gradient part 2: g = -2*S*T*X; X = Q*R
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, qmat, scmat, 0.0_dp, t2mat)
      CALL cp_dbcsr_multiply("N", "N", 1.0_dp, tmat, t2mat, 0.0_dp, scmat)
      CALL cp_dbcsr_multiply("N", "N", -2.0_dp, smat, scmat, 1.0_dp, mao_grad, &
                             retain_sparsity=.TRUE.)
      ! free temp matrices
      CALL cp_dbcsr_release(scmat)
      CALL cp_dbcsr_release(bmat)
      CALL cp_dbcsr_release(tmat)
      CALL cp_dbcsr_release(t2mat)

      CALL mao_project_gradient(mao_coef, mao_grad, smat)

   END SUBROUTINE mao_function_gradient

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param smat ...
! **************************************************************************************************
   SUBROUTINE mao_orthogonalization(mao_coef, smat)
      TYPE(cp_dbcsr_type)                                :: mao_coef, smat

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_orthogonalization', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iatom, info, jatom, lwork, m, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, bmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cblock, sblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, mao_coef)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, cblock)
         CPASSERT(iatom == jatom)
         m = SIZE(cblock, 2)
         n = SIZE(cblock, 1)
         NULLIFY (sblock)
         CALL cp_dbcsr_get_block_p(matrix=smat, row=iatom, col=jatom, block=sblock, found=found)
         CPASSERT(found)
         lwork = MAX(n*n, 100)
         ALLOCATE (amat(n, m), bmat(m, m), w(m), work(lwork))
         amat(1:n, 1:m) = MATMUL(sblock(1:n, 1:n), cblock(1:n, 1:m))
         bmat(1:m, 1:m) = MATMUL(TRANSPOSE(cblock(1:n, 1:m)), amat(1:n, 1:m))
         info = 0
         CALL lapack_ssyev("V", "U", m, bmat, m, w, work, lwork, info)
         CPASSERT(info == 0)
         CPASSERT(ALL(w > 0.0_dp))
         w = 1.0_dp/SQRT(w)
         DO i = 1, m
            amat(1:m, i) = bmat(1:m, i)*w(i)
         END DO
         bmat(1:m, 1:m) = MATMUL(amat(1:m, 1:m), TRANSPOSE(bmat(1:m, 1:m)))
         cblock(1:n, 1:m) = MATMUL(cblock(1:n, 1:m), bmat(1:m, 1:m))
         DEALLOCATE (amat, bmat, w, work)
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE mao_orthogonalization

! **************************************************************************************************
!> \brief ...
!> \param mao_coef ...
!> \param mao_grad ...
!> \param smat ...
! **************************************************************************************************
   SUBROUTINE mao_project_gradient(mao_coef, mao_grad, smat)
      TYPE(cp_dbcsr_type)                                :: mao_coef, mao_grad, smat

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_project_gradient', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, jatom, m, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cblock, gblock, sblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, mao_coef)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, cblock)
         CPASSERT(iatom == jatom)
         m = SIZE(cblock, 2)
         n = SIZE(cblock, 1)
         NULLIFY (sblock)
         CALL cp_dbcsr_get_block_p(matrix=smat, row=iatom, col=jatom, block=sblock, found=found)
         CPASSERT(found)
         NULLIFY (gblock)
         CALL cp_dbcsr_get_block_p(matrix=mao_grad, row=iatom, col=jatom, block=gblock, found=found)
         CPASSERT(found)
         ALLOCATE (amat(m, m))
         amat(1:m, 1:m) = MATMUL(TRANSPOSE(cblock(1:n, 1:m)), MATMUL(sblock(1:n, 1:n), gblock(1:n, 1:m)))
         gblock(1:n, 1:m) = gblock(1:n, 1:m)-MATMUL(cblock(1:n, 1:m), amat(1:m, 1:m))
         DEALLOCATE (amat)
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE mao_project_gradient

! **************************************************************************************************
!> \brief ...
!> \param fmat1 ...
!> \param fmat2 ...
!> \retval spro ...
! **************************************************************************************************
   FUNCTION mao_scalar_product(fmat1, fmat2) RESULT(spro)
      TYPE(cp_dbcsr_type)                                :: fmat1, fmat2
      REAL(KIND=dp)                                      :: spro

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_scalar_product', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, iatom, jatom, m, n
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ablock, bblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      spro = 0.0_dp

      CALL cp_dbcsr_iterator_start(dbcsr_iter, fmat1)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, ablock)
         CPASSERT(iatom == jatom)
         m = SIZE(ablock, 2)
         n = SIZE(ablock, 1)
         CALL cp_dbcsr_get_block_p(matrix=fmat2, row=iatom, col=jatom, block=bblock, found=found)
         CPASSERT(found)
         spro = spro+SUM(ablock(1:n, 1:m)*bblock(1:n, 1:m))
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

      group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(fmat1)))
      CALL mp_sum(spro, group)

   END FUNCTION mao_scalar_product

! **************************************************************************************************
!> \brief Calculate the density matrix at the Gamma point
!> \param pmat ...
!> \param ksmat ...
!> \param smat ...
!> \param kpoints      Kpoint environment
!> \param nmos         Number of occupied orbitals
!> \param occ          Maximum occupation per orbital
!> \par History
!>      04.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE calculate_p_gamma(pmat, ksmat, smat, kpoints, nmos, occ)

      TYPE(cp_dbcsr_type)                                :: pmat, ksmat, smat
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: nmos
      REAL(KIND=dp), INTENT(IN)                          :: occ

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_p_gamma', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: norb
      REAL(KIND=dp)                                      :: de
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_dbcsr_type)                                :: tempmat
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: fmksmat, fmsmat, fmvec, fmwork

      ! FM matrices

      CALL cp_dbcsr_get_info(smat, nfullrows_total=norb)
      CALL cp_fm_struct_create(fmstruct=matrix_struct, context=kpoints%blacs_env_all, &
                               nrow_global=norb, ncol_global=norb)
      CALL cp_fm_create(fmksmat, matrix_struct)
      CALL cp_fm_create(fmsmat, matrix_struct)
      CALL cp_fm_create(fmvec, matrix_struct)
      CALL cp_fm_create(fmwork, matrix_struct)
      ALLOCATE (eigenvalues(norb))

      ! DBCSR matrix
      CALL cp_dbcsr_init(tempmat)
      CALL cp_dbcsr_create(tempmat, template=smat, matrix_type=dbcsr_type_no_symmetry)

      ! transfer to FM
      CALL cp_dbcsr_desymmetrize(smat, tempmat)
      CALL copy_dbcsr_to_fm(tempmat, fmsmat)
      CALL cp_dbcsr_desymmetrize(ksmat, tempmat)
      CALL copy_dbcsr_to_fm(tempmat, fmksmat)

      ! diagonalize
      CALL cp_fm_geeig(fmksmat, fmsmat, fmvec, eigenvalues, fmwork)
      de = eigenvalues(nmos+1)-eigenvalues(nmos)
      IF (de < 0.001_dp) THEN
         CALL cp_warn(__LOCATION__, "MAO: No band gap at "// &
                      "Gamma point. MAO analysis not reliable.")
      END IF
      ! density matrix
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=pmat, matrix_v=fmvec, ncol=nmos, alpha=occ)

      DEALLOCATE (eigenvalues)
      CALL cp_dbcsr_release(tempmat)
      CALL cp_fm_release(fmksmat)
      CALL cp_fm_release(fmsmat)
      CALL cp_fm_release(fmvec)
      CALL cp_fm_release(fmwork)
      CALL cp_fm_struct_release(matrix_struct)

   END SUBROUTINE calculate_p_gamma

! **************************************************************************************************
!> \brief Define the MAO reference basis set
!> \param qs_env ...
!> \param mao_basis ...
!> \param mao_basis_set_list ...
!> \param orb_basis_set_list ...
!> \param iunit ...
!> \param print_basis ...
!> \par History
!>      07.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE mao_reference_basis(qs_env, mao_basis, mao_basis_set_list, orb_basis_set_list, &
                                  iunit, print_basis)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: mao_basis
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: mao_basis_set_list, orb_basis_set_list
      INTEGER, INTENT(IN), OPTIONAL                      :: iunit
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_basis

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_reference_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, nbas, nkind, unit_nr
      REAL(KIND=dp)                                      :: eps_pgf_orb
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, pbasis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      ! Reference basis set
      CPASSERT(.NOT. ASSOCIATED(mao_basis_set_list))
      CPASSERT(.NOT. ASSOCIATED(orb_basis_set_list))

      ! options
      IF (PRESENT(iunit)) THEN
         unit_nr = iunit
      ELSE
         unit_nr = -1
      END IF

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (mao_basis_set_list(nkind), orb_basis_set_list(nkind))
      DO ikind = 1, nkind
         NULLIFY (mao_basis_set_list(ikind)%gto_basis_set)
         NULLIFY (orb_basis_set_list(ikind)%gto_basis_set)
      END DO
      !
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
         IF (ASSOCIATED(basis_set)) orb_basis_set_list(ikind)%gto_basis_set => basis_set
      END DO
      !
      SELECT CASE (mao_basis)
      CASE (mao_basis_orb)
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
            IF (ASSOCIATED(basis_set)) mao_basis_set_list(ikind)%gto_basis_set => basis_set
         END DO
      CASE (mao_basis_prim)
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
            NULLIFY (pbasis)
            IF (ASSOCIATED(basis_set)) THEN
               CALL create_primitive_basis_set(basis_set, pbasis)
               CALL get_qs_env(qs_env, dft_control=dft_control)
               eps_pgf_orb = dft_control%qs_control%eps_pgf_orb
               CALL init_interaction_radii_orb_basis(pbasis, eps_pgf_orb)
               pbasis%kind_radius = basis_set%kind_radius
               mao_basis_set_list(ikind)%gto_basis_set => pbasis
               CALL add_basis_set_to_container(qs_kind%basis_sets, pbasis, "MAO")
            END IF
         END DO
      CASE (mao_basis_ext)
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="MAO")
            IF (ASSOCIATED(basis_set)) THEN
               basis_set%kind_radius = orb_basis_set_list(ikind)%gto_basis_set%kind_radius
               mao_basis_set_list(ikind)%gto_basis_set => basis_set
            END IF
         END DO
      CASE DEFAULT
         CPABORT("Unknown option for MAO basis")
      END SELECT
      IF (unit_nr > 0) THEN
         DO ikind = 1, nkind
            IF (.NOT. ASSOCIATED(mao_basis_set_list(ikind)%gto_basis_set)) THEN
               WRITE (UNIT=unit_nr, FMT="(T2,A,I4)") &
                  "WARNING: No MAO basis set associated with Kind ", ikind
            ELSE
               nbas = mao_basis_set_list(ikind)%gto_basis_set%nsgf
               WRITE (UNIT=unit_nr, FMT="(T2,A,I4,T56,A,I10)") &
                  "MAO basis set Kind ", ikind, " Number of BSF:", nbas
            END IF
         END DO
      END IF

      IF (PRESENT(print_basis)) THEN
         IF (print_basis) THEN
            DO ikind = 1, nkind
               basis_set => mao_basis_set_list(ikind)%gto_basis_set
               CALL write_gto_basis_set(basis_set, unit_nr, "MAO REFERENCE BASIS")
            END DO
         END IF
      END IF

   END SUBROUTINE mao_reference_basis

! **************************************************************************************************
!> \brief Analyze the MAO basis, projection on angular functions
!> \param mao_coef ...
!> \param matrix_smm ...
!> \param mao_basis_set_list ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param unit_nr ...
!> \param para_env ...
!> \par History
!>      07.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE mao_basis_analysis(mao_coef, matrix_smm, mao_basis_set_list, particle_set, &
                                 qs_kind_set, unit_nr, para_env)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_coef, matrix_smm
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: mao_basis_set_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_basis_analysis', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol
      INTEGER                                            :: ia, iab, iatom, ikind, iset, ishell, &
                                                            ispin, l, lmax, lshell, m, ma, na, &
                                                            natom, nspin
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cmask, vec1, vec2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weight
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block, cmao
      TYPE(gto_basis_set_type), POINTER                  :: basis_set

      ! Analyze the MAO basis
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(/,A)") " Analyze angular momentum character of MAOs "
         WRITE (unit_nr, "(T7,A,T15,A,T20,A,T40,A,T50,A,T60,A,T70,A,T80,A)") &
            "ATOM", "Spin", "MAO", "S", "P", "D", "F", "G"
      END IF
      lmax = 4 ! analyze up to g-functions
      natom = SIZE(particle_set)
      nspin = SIZE(mao_coef)
      DO iatom = 1, natom
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                              element_symbol=element_symbol, kind_number=ikind)
         basis_set => mao_basis_set_list(ikind)%gto_basis_set
         CALL get_qs_kind(qs_kind_set(ikind), mao=na)
         CALL get_gto_basis_set(basis_set, nsgf=ma)
         ALLOCATE (cmask(ma), vec1(ma), vec2(ma), weight(0:lmax, na))
         weight = 0.0_dp
         CALL cp_dbcsr_get_block_p(matrix=matrix_smm(1)%matrix, row=iatom, col=iatom, &
                                   block=block, found=found)
         DO ispin = 1, nspin
            CALL cp_dbcsr_get_block_p(matrix=mao_coef(ispin)%matrix, row=iatom, col=iatom, &
                                      block=cmao, found=found)
            IF (found) THEN
               DO l = 0, lmax
                  cmask = 0.0_dp
                  iab = 0
                  DO iset = 1, basis_set%nset
                     DO ishell = 1, basis_set%nshell(iset)
                        lshell = basis_set%l(ishell, iset)
                        DO m = -lshell, lshell
                           iab = iab+1
                           IF (l == lshell) cmask(iab) = 1.0_dp
                        END DO
                     END DO
                  END DO
                  DO ia = 1, na
                     vec1(1:ma) = cmask*cmao(1:ma, ia)
                     vec2(1:ma) = MATMUL(block, vec1)
                     weight(l, ia) = SUM(vec1(1:ma)*vec2(1:ma))
                  END DO
               END DO
            END IF
            CALL mp_sum(weight, para_env%group)
            IF (unit_nr > 0) THEN
               DO ia = 1, na
                  IF (ispin == 1 .AND. ia == 1) THEN
                     WRITE (unit_nr, "(i6,T9,A2,T17,i2,T20,i3,T31,5F10.4)") &
                        iatom, element_symbol, ispin, ia, weight(0:lmax, ia)
                  ELSE
                     WRITE (unit_nr, "(T17,i2,T20,i3,T31,5F10.4)") ispin, ia, weight(0:lmax, ia)
                  END IF
               END DO
            END IF
         END DO
         DEALLOCATE (cmask, weight, vec1, vec2)
      END DO
   END SUBROUTINE mao_basis_analysis

! **************************************************************************************************
!> \brief Calculte the Q=APA(T) matrix, A=(MAO,ORB) overlap
!> \param matrix_q ...
!> \param matrix_p ...
!> \param matrix_s ...
!> \param matrix_smm ...
!> \param matrix_smo ...
!> \param smm_list ...
!> \param electra ...
!> \param eps_filter ...
!> \param nimages ...
!> \param kpoints ...
!> \param matrix_ks ...
!> \param sab_orb ...
!> \par History
!>      08.2016 created [JGH]
! **************************************************************************************************
   SUBROUTINE mao_build_q(matrix_q, matrix_p, matrix_s, matrix_smm, matrix_smo, smm_list, &
                          electra, eps_filter, nimages, kpoints, matrix_ks, sab_orb)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_q
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: matrix_p, matrix_s
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_smm, matrix_smo
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: smm_list
      REAL(KIND=dp), DIMENSION(2), INTENT(OUT)           :: electra
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN), OPTIONAL                      :: nimages
      TYPE(kpoint_type), OPTIONAL, POINTER               :: kpoints
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_ks
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb

      CHARACTER(len=*), PARAMETER :: routineN = 'mao_build_q', routineP = moduleN//':'//routineN

      INTEGER                                            :: im, ispin, nim, nocc, norb, nspin
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp)                                      :: elex, xkp(3)
      TYPE(cp_dbcsr_type)                                :: ksmat, pmat, smat, tmat

      nim = 1
      IF (PRESENT(nimages)) nim = nimages
      IF (nim > 1) THEN
         CPASSERT(PRESENT(kpoints))
         CPASSERT(PRESENT(matrix_ks))
         CPASSERT(PRESENT(sab_orb))
      END IF

      ! Reference
      nspin = SIZE(matrix_p, 1)
      DO ispin = 1, nspin
         electra(ispin) = 0.0_dp
         DO im = 1, nim
            CALL cp_dbcsr_trace(matrix_p(ispin, im)%matrix, matrix_s(1, im)%matrix, elex)
            electra(ispin) = electra(ispin)+elex
         END DO
      END DO

      ! Q matrix
      NULLIFY (matrix_q)
      CALL cp_dbcsr_allocate_matrix_set(matrix_q, nspin)
      DO ispin = 1, nspin
         ALLOCATE (matrix_q(ispin)%matrix)
         CALL cp_dbcsr_init(matrix_q(ispin)%matrix)
         CALL cp_dbcsr_create(matrix_q(ispin)%matrix, template=matrix_smm(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_q(ispin)%matrix, smm_list)
      END DO
      ! temp matrix
      CALL cp_dbcsr_init(tmat)
      CALL cp_dbcsr_create(tmat, template=matrix_smo(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      ! Q=APA(T)
      DO ispin = 1, nspin
         IF (nim == 1) THEN
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_smo(1)%matrix, matrix_p(ispin, 1)%matrix, &
                                   0.0_dp, tmat, filter_eps=eps_filter)
            CALL cp_dbcsr_multiply("N", "T", 1.0_dp, tmat, matrix_smo(1)%matrix, &
                                   0.0_dp, matrix_q(ispin)%matrix, filter_eps=eps_filter)
         ELSE
            ! k-points
            CALL cp_dbcsr_init(pmat)
            CALL cp_dbcsr_create(pmat, template=matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_init(smat)
            CALL cp_dbcsr_create(smat, template=matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_init(ksmat)
            CALL cp_dbcsr_create(ksmat, template=matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(pmat, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(smat, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(ksmat, sab_orb)
            NULLIFY (cell_to_index)
            CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
            ! calculate density matrix at gamma point
            xkp = 0.0_dp
            ! transform KS and S matrices to the gamma point
            CALL cp_dbcsr_set(ksmat, 0.0_dp)
            CALL rskp_transform(rmatrix=ksmat, rsmat=matrix_ks, ispin=ispin, &
                                xkp=xkp, cell_to_index=cell_to_index, sab_nl=sab_orb)
            CALL cp_dbcsr_set(smat, 0.0_dp)
            CALL rskp_transform(rmatrix=smat, rsmat=matrix_s, ispin=1, &
                                xkp=xkp, cell_to_index=cell_to_index, sab_nl=sab_orb)
            norb = NINT(electra(ispin))
            nocc = MOD(2, nspin)+1
            CALL calculate_p_gamma(pmat, ksmat, smat, kpoints, norb, REAL(nocc, KIND=dp))
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_smo(1)%matrix, pmat, &
                                   0.0_dp, tmat, filter_eps=eps_filter)
            CALL cp_dbcsr_multiply("N", "T", 1.0_dp, tmat, matrix_smo(1)%matrix, &
                                   0.0_dp, matrix_q(ispin)%matrix, filter_eps=eps_filter)
            CALL cp_dbcsr_release(pmat)
            CALL cp_dbcsr_release(smat)
            CALL cp_dbcsr_release(ksmat)
         END IF
      END DO
      ! free temp matrix
      CALL cp_dbcsr_release(tmat)

   END SUBROUTINE mao_build_q

END MODULE mao_methods
