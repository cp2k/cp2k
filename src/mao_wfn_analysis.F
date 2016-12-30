!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate MAO's and analyze wavefunctions
!> \par History
!>      03.2016 created [JGH]
!>      12.2016 split into four modules [JGH]
!> \author JGH
! **************************************************************************************************
MODULE mao_wfn_analysis
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE bibliography,                    ONLY: Ehrhardt1985,&
                                              Heinzmann1976,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_restore
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_allocate_matrix_set, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix_set, &
        dbcsr_desymmetrize, dbcsr_distribution_obj, dbcsr_get_block_diag, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_init, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_replicate_all, dbcsr_reserve_diag_blocks, dbcsr_trace, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE mao_methods,                     ONLY: mao_basis_analysis,&
                                              mao_build_q,&
                                              mao_reference_basis
   USE mao_optimizer,                   ONLY: mao_optimize
   USE mathlib,                         ONLY: invmat_symm
   USE message_passing,                 ONLY: mp_sum
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mao_wfn_analysis'

   PUBLIC ::  mao_analysis

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
      INTEGER :: fall, handle, ia, iab, iabc, iatom, ib, ic, icol, ikind, irow, ispin, jatom, &
         mao_basis, max_iter, me, na, nab, nabc, natom, nb, nc, nimages, nspin, ssize
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, mao_blk, mao_blk_sizes, &
                                                            orb_blk, row_blk_sizes
      LOGICAL                                            :: analyze_ua, explicit, fo, for, fos, &
                                                            found, neglect_abc, print_basis
      REAL(KIND=dp) :: deltaq, electra(2), eps_ab, eps_abc, eps_filter, eps_fun, eps_grad, senabc, &
         senmax, threshold, total_charge, total_spin, ua_charge(2), zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: occnumA, occnumABC, qab, qmatab, qmatac, &
                                                            qmatbc, raq, sab, selnABC, sinv, &
                                                            smatab, smatac, smatbc, uaq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: occnumAB, selnAB
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block, cmao, diag, qblka, qblkb, qblkc, &
                                                            rblkl, rblku, sblk, sblka, sblkb, sblkc
      TYPE(block_type), ALLOCATABLE, DIMENSION(:)        :: rowblock
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist
      TYPE(dbcsr_iterator_type)                          :: dbcsr_iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef, mao_dmat, mao_qmat, mao_smat, &
                                                            matrix_q, matrix_smm, matrix_smo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_p, matrix_s
      TYPE(dbcsr_type)                                   :: amat, axmat, cgmat, cholmat, crumat, &
                                                            qmat, qmat_diag, rumat, smat_diag, &
                                                            sumat, tmat
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
         CALL dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
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
      NULLIFY (mao_coef)
      CALL dbcsr_allocate_matrix_set(mao_coef, nspin)
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
         CALL dbcsr_init(mao_coef(ispin)%matrix)
         CALL dbcsr_create(matrix=mao_coef(ispin)%matrix, &
                           name="MAO_COEF", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL dbcsr_reserve_diag_blocks(matrix=mao_coef(ispin)%matrix)
      END DO
      DEALLOCATE (row_blk_sizes, col_blk_sizes)

      ! optimze MAOs
      CALL mao_optimize(mao_coef, matrix_q, matrix_smm, electra, max_iter, eps_grad, unit_nr)

      ! Analyze the MAO basis
      CALL mao_basis_analysis(mao_coef, matrix_smm, mao_basis_set_list, particle_set, &
                              qs_kind_set, unit_nr, para_env)

      ! Calculate the overlap and density matrix in the new MAO basis
      NULLIFY (mao_dmat, mao_smat, mao_qmat)
      CALL dbcsr_allocate_matrix_set(mao_qmat, nspin)
      CALL dbcsr_allocate_matrix_set(mao_dmat, nspin)
      CALL dbcsr_allocate_matrix_set(mao_smat, nspin)
      CALL dbcsr_get_info(mao_coef(1)%matrix, col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
      DO ispin = 1, nspin
         ALLOCATE (mao_dmat(ispin)%matrix)
         CALL dbcsr_init(mao_dmat(ispin)%matrix)
         CALL dbcsr_create(mao_dmat(ispin)%matrix, name="MAO density", dist=dbcsr_dist, &
                           matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                           col_blk_size=col_blk_sizes, nze=0)
         ALLOCATE (mao_smat(ispin)%matrix)
         CALL dbcsr_init(mao_smat(ispin)%matrix)
         CALL dbcsr_create(mao_smat(ispin)%matrix, name="MAO overlap", dist=dbcsr_dist, &
                           matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                           col_blk_size=col_blk_sizes, nze=0)
         ALLOCATE (mao_qmat(ispin)%matrix)
         CALL dbcsr_init(mao_qmat(ispin)%matrix)
         CALL dbcsr_create(mao_qmat(ispin)%matrix, name="MAO covar density", dist=dbcsr_dist, &
                           matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                           col_blk_size=col_blk_sizes, nze=0)
      END DO
      CALL dbcsr_init(amat)
      CALL dbcsr_create(amat, name="MAO overlap", template=mao_dmat(1)%matrix)
      CALL dbcsr_init(tmat)
      CALL dbcsr_create(tmat, name="MAO Overlap Inverse", template=amat)
      CALL dbcsr_init(qmat)
      CALL dbcsr_create(qmat, name="MAO covar density", template=amat)
      CALL dbcsr_init(cgmat)
      CALL dbcsr_create(cgmat, name="TEMP matrix", template=mao_coef(1)%matrix)
      CALL dbcsr_init(axmat)
      CALL dbcsr_create(axmat, name="TEMP", template=amat, matrix_type=dbcsr_type_no_symmetry)
      DO ispin = 1, nspin
         ! calculate MAO overlap matrix
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_smm(1)%matrix, mao_coef(ispin)%matrix, &
                             0.0_dp, cgmat)
         CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, amat)
         ! calculate inverse of MAO overlap
         threshold = 1.e-8_dp
         CALL invert_Hotelling(tmat, amat, threshold, norm_convergence=1.e-4_dp, silent=.TRUE.)
         CALL dbcsr_copy(mao_smat(ispin)%matrix, amat)
         ! calculate q-matrix q = C*Q*C
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_q(ispin)%matrix, mao_coef(ispin)%matrix, &
                             0.0_dp, cgmat, filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, &
                             0.0_dp, qmat, filter_eps=eps_filter)
         CALL dbcsr_copy(mao_qmat(ispin)%matrix, qmat)
         ! calculate density matrix
         CALL dbcsr_multiply("N", "N", 1.0_dp, qmat, tmat, 0.0_dp, axmat, filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, tmat, axmat, 0.0_dp, mao_dmat(ispin)%matrix, &
                             filter_eps=eps_filter)
      END DO
      CALL dbcsr_release(amat)
      CALL dbcsr_release(tmat)
      CALL dbcsr_release(qmat)
      CALL dbcsr_release(cgmat)
      CALL dbcsr_release(axmat)

      ! calculate unassigned charge : n - Tr PS
      DO ispin = 1, nspin
         CALL dbcsr_trace(mao_dmat(ispin)%matrix, mao_smat(ispin)%matrix, ua_charge(ispin))
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
            CALL dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, &
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
         CALL dbcsr_init(qmat_diag)
         CALL dbcsr_create(qmat_diag, name="MAO diagonal density", template=mao_dmat(1)%matrix)
         CALL dbcsr_init(smat_diag)
         CALL dbcsr_create(smat_diag, name="MAO diagonal overlap", template=mao_dmat(1)%matrix)
         ! replicate the diagonal blocks of the density and overlap matrices
         CALL dbcsr_get_block_diag(mao_qmat(ispin)%matrix, qmat_diag)
         CALL dbcsr_replicate_all(qmat_diag)
         CALL dbcsr_get_block_diag(mao_smat(ispin)%matrix, smat_diag)
         CALL dbcsr_replicate_all(smat_diag)
         DO ia = 1, natom
            DO ib = ia+1, natom
               iab = 0
               CALL dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, &
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
                  CALL dbcsr_get_block_p(matrix=qmat_diag, row=ia, col=ia, block=diag, found=fo)
                  CPASSERT(fo)
                  qab(1:na, 1:na) = diag(1:na, 1:na)
                  CALL dbcsr_get_block_p(matrix=qmat_diag, row=ib, col=ib, block=diag, found=fo)
                  CPASSERT(fo)
                  qab(na+1:nab, na+1:nab) = diag(1:nb, 1:nb)
                  ! smat
                  CALL dbcsr_get_block_p(matrix=mao_smat(ispin)%matrix, &
                                         row=ia, col=ib, block=block, found=fo)
                  CPASSERT(fo)
                  sab(1:na, na+1:nab) = block(1:na, 1:nb)
                  sab(na+1:nab, 1:na) = TRANSPOSE(block(1:na, 1:nb))
                  CALL dbcsr_get_block_p(matrix=smat_diag, row=ia, col=ia, block=diag, found=fo)
                  CPASSERT(fo)
                  sab(1:na, 1:na) = diag(1:na, 1:na)
                  CALL dbcsr_get_block_p(matrix=smat_diag, row=ib, col=ib, block=diag, found=fo)
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
         CALL dbcsr_release(qmat_diag)
         CALL dbcsr_release(smat_diag)
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
            CALL dbcsr_init(qmat_diag)
            CALL dbcsr_create(qmat_diag, name="MAO diagonal density", template=mao_dmat(1)%matrix)
            CALL dbcsr_init(smat_diag)
            CALL dbcsr_create(smat_diag, name="MAO diagonal overlap", template=mao_dmat(1)%matrix)
            ! replicate the diagonal blocks of the density and overlap matrices
            CALL dbcsr_get_block_diag(mao_qmat(ispin)%matrix, qmat_diag)
            CALL dbcsr_replicate_all(qmat_diag)
            CALL dbcsr_get_block_diag(mao_smat(ispin)%matrix, smat_diag)
            CALL dbcsr_replicate_all(smat_diag)
            iabc = 0
            DO ia = 1, natom
               CALL dbcsr_get_block_p(matrix=qmat_diag, row=ia, col=ia, block=qblka, found=fo)
               CPASSERT(fo)
               CALL dbcsr_get_block_p(matrix=smat_diag, row=ia, col=ia, block=sblka, found=fo)
               CPASSERT(fo)
               na = SIZE(qblka, 1)
               DO ib = ia+1, natom
                  ! screen with SEN(AB)
                  IF (selnAB(ia, ib, ispin) < eps_abc) THEN
                     iabc = iabc+(natom-ib)
                     CYCLE
                  END IF
                  CALL dbcsr_get_block_p(matrix=qmat_diag, row=ib, col=ib, block=qblkb, found=fo)
                  CPASSERT(fo)
                  CALL dbcsr_get_block_p(matrix=smat_diag, row=ib, col=ib, block=sblkb, found=fo)
                  CPASSERT(fo)
                  nb = SIZE(qblkb, 1)
                  nab = na+nb
                  ALLOCATE (qmatab(na, nb), smatab(na, nb))
                  CALL dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, row=ia, col=ib, &
                                         block=block, found=found)
                  qmatab = 0.0_dp
                  IF (found) qmatab(1:na, 1:nb) = block(1:na, 1:nb)
                  CALL mp_sum(qmatab, para_env%group)
                  CALL dbcsr_get_block_p(matrix=mao_smat(ispin)%matrix, row=ia, col=ib, &
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
                     CALL dbcsr_get_block_p(matrix=qmat_diag, row=ic, col=ic, block=qblkc, found=fo)
                     CPASSERT(fo)
                     CALL dbcsr_get_block_p(matrix=smat_diag, row=ic, col=ic, block=sblkc, found=fo)
                     CPASSERT(fo)
                     nc = SIZE(qblkc, 1)
                     ALLOCATE (qmatac(na, nc), smatac(na, nc))
                     CALL dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, row=ia, col=ic, &
                                            block=block, found=found)
                     qmatac = 0.0_dp
                     IF (found) qmatac(1:na, 1:nc) = block(1:na, 1:nc)
                     CALL mp_sum(qmatac, para_env%group)
                     CALL dbcsr_get_block_p(matrix=mao_smat(ispin)%matrix, row=ia, col=ic, &
                                            block=block, found=found)
                     smatac = 0.0_dp
                     IF (found) smatac(1:na, 1:nc) = block(1:na, 1:nc)
                     CALL mp_sum(smatac, para_env%group)
                     ALLOCATE (qmatbc(nb, nc), smatbc(nb, nc))
                     CALL dbcsr_get_block_p(matrix=mao_qmat(ispin)%matrix, row=ib, col=ic, &
                                            block=block, found=found)
                     qmatbc = 0.0_dp
                     IF (found) qmatbc(1:nb, 1:nc) = block(1:nb, 1:nc)
                     CALL mp_sum(qmatbc, para_env%group)
                     CALL dbcsr_get_block_p(matrix=mao_smat(ispin)%matrix, row=ib, col=ic, &
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
            CALL dbcsr_release(qmat_diag)
            CALL dbcsr_release(smat_diag)
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
         CALL dbcsr_get_info(mao_coef(1)%matrix, row_blk_size=mao_blk_sizes, &
                             col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
         CALL dbcsr_get_info(matrix_s(1, 1)%matrix, row_blk_size=row_blk_sizes)
         CALL dbcsr_init(amat)
         CALL dbcsr_create(amat, name="temp", template=matrix_s(1, 1)%matrix)
         CALL dbcsr_init(tmat)
         CALL dbcsr_create(tmat, name="temp", template=mao_coef(1)%matrix)
         ! replicate diagonal of smm matrix
         CALL dbcsr_init(smat_diag)
         CALL dbcsr_get_block_diag(matrix_smm(1)%matrix, smat_diag)
         CALL dbcsr_replicate_all(smat_diag)

         ALLOCATE (orb_blk(natom), mao_blk(natom))
         DO ia = 1, natom
            orb_blk = row_blk_sizes
            mao_blk = row_blk_sizes
            mao_blk(ia) = col_blk_sizes(ia)
            CALL dbcsr_init(sumat)
            CALL dbcsr_create(sumat, name="Smat", dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                              row_blk_size=mao_blk, col_blk_size=mao_blk, nze=0)
            CALL cp_dbcsr_alloc_block_from_nbl(sumat, sab_orb)
            CALL dbcsr_init(cholmat)
            CALL dbcsr_create(cholmat, name="Cholesky matrix", dist=dbcsr_dist, &
                              matrix_type=dbcsr_type_no_symmetry, row_blk_size=mao_blk, col_blk_size=mao_blk, nze=0)
            CALL dbcsr_init(rumat)
            CALL dbcsr_create(rumat, name="Rmat", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=orb_blk, col_blk_size=mao_blk, nze=0)
            CALL cp_dbcsr_alloc_block_from_nbl(rumat, sab_orb, .TRUE.)
            CALL dbcsr_init(crumat)
            CALL dbcsr_create(crumat, name="Rmat*Umat", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=orb_blk, col_blk_size=mao_blk, nze=0)
            ! replicate row and col of smo matrix
            ALLOCATE (rowblock(natom))
            DO ib = 1, natom
               na = mao_blk_sizes(ia)
               nb = row_blk_sizes(ib)
               ALLOCATE (rowblock(ib)%mat(na, nb))
               rowblock(ib)%mat = 0.0_dp
               CALL dbcsr_get_block_p(matrix=matrix_smo(1)%matrix, row=ia, col=ib, &
                                      block=block, found=found)
               IF (found) rowblock(ib)%mat(1:na, 1:nb) = block(1:na, 1:nb)
               CALL mp_sum(rowblock(ib)%mat, para_env%group)
            END DO
            !
            DO ispin = 1, nspin
               CALL dbcsr_copy(tmat, mao_coef(ispin)%matrix)
               CALL dbcsr_replicate_all(tmat)
               CALL dbcsr_iterator_start(dbcsr_iter, matrix_s(1, 1)%matrix)
               DO WHILE (dbcsr_iterator_blocks_left(dbcsr_iter))
                  CALL dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, block)
                  CALL dbcsr_get_block_p(matrix=sumat, row=iatom, col=jatom, block=sblk, found=fos)
                  CPASSERT(fos)
                  CALL dbcsr_get_block_p(matrix=rumat, row=iatom, col=jatom, block=rblku, found=for)
                  CPASSERT(for)
                  CALL dbcsr_get_block_p(matrix=rumat, row=jatom, col=iatom, block=rblkl, found=for)
                  CPASSERT(for)
                  CALL dbcsr_get_block_p(matrix=tmat, row=ia, col=ia, block=cmao, found=found)
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
                     CALL dbcsr_get_block_p(matrix=smat_diag, row=ia, col=ia, block=block, found=found)
                     CPASSERT(found)
                     sblk = MATMUL(TRANSPOSE(cmao), MATMUL(block, cmao))
                     rblku = MATMUL(TRANSPOSE(rowblock(ia)%mat), cmao)
                  END IF
               END DO
               CALL dbcsr_iterator_stop(dbcsr_iter)
               ! Cholesky decomposition of SUMAT = U'U
               CALL dbcsr_desymmetrize(sumat, cholmat)
               CALL cp_dbcsr_cholesky_decompose(cholmat, para_env=para_env, blacs_env=blacs_env)
               ! T = R*inv(U)
               ssize = SUM(mao_blk)
               CALL cp_dbcsr_cholesky_restore(rumat, ssize, cholmat, crumat, op="SOLVE", pos="RIGHT", &
                                              transa="N", para_env=para_env, blacs_env=blacs_env)
               ! A = T*transpose(T)
               CALL dbcsr_multiply("N", "T", 1.0_dp, crumat, crumat, 0.0_dp, amat, &
                                   filter_eps=eps_filter)
               ! Tr(P*A)
               CALL dbcsr_trace(matrix_p(ispin, 1)%matrix, amat, uaq(ia, ispin))
               uaq(ia, ispin) = uaq(ia, ispin)-electra(ispin)
            END DO
            !
            CALL dbcsr_release(sumat)
            CALL dbcsr_release(cholmat)
            CALL dbcsr_release(rumat)
            CALL dbcsr_release(crumat)
            !
            DO ib = 1, natom
               DEALLOCATE (rowblock(ib)%mat)
            END DO
            DEALLOCATE (rowblock)
         END DO
         CALL dbcsr_release(smat_diag)
         CALL dbcsr_release(amat)
         CALL dbcsr_release(tmat)
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

      IF (ASSOCIATED(matrix_smm)) CALL dbcsr_deallocate_matrix_set(matrix_smm)
      IF (ASSOCIATED(matrix_smo)) CALL dbcsr_deallocate_matrix_set(matrix_smo)
      IF (ASSOCIATED(matrix_q)) CALL dbcsr_deallocate_matrix_set(matrix_q)

      IF (ASSOCIATED(mao_coef)) CALL dbcsr_deallocate_matrix_set(mao_coef)
      IF (ASSOCIATED(mao_dmat)) CALL dbcsr_deallocate_matrix_set(mao_dmat)
      IF (ASSOCIATED(mao_smat)) CALL dbcsr_deallocate_matrix_set(mao_smat)
      IF (ASSOCIATED(mao_qmat)) CALL dbcsr_deallocate_matrix_set(mao_qmat)

      CALL timestop(handle)

   END SUBROUTINE mao_analysis

END MODULE mao_wfn_analysis
