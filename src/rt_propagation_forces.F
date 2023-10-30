!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

MODULE rt_propagation_forces
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_get_block_p, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: one,&
                                              zero
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: add_qs_force,&
                                              qs_force_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_force
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: calc_c_mat_force, &
             rt_admm_force

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_forces'

CONTAINS

! **************************************************************************************************
!> \brief calculates the three additional force contributions needed in EMD
!>        P_imag*C , P_imag*B*S^-1*S_der , P*S^-1*H*S_der
!> \param qs_env ...
!> \par History
!>      02.2014 switched to dbcsr matrices [Samuel Andermatt]
!>      10.2023 merge MO-based and all-atom into one routine [Guillaume Le Breton]
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE calc_c_mat_force(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_c_mat_force'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: handle, i, im, ispin, re
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: C_mat, rho_ao, rho_ao_im, rho_new, &
                                                            S_der, SinvB, SinvH, SinvH_imag
      TYPE(dbcsr_type), POINTER                          :: S_inv, tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      CALL timeset(routineN, handle)
      NULLIFY (rtp, particle_set, atomic_kind_set, dft_control)

      CALL get_qs_env(qs_env, &
                      rtp=rtp, &
                      rho=rho, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      force=force, &
                      dft_control=dft_control)

      rtp_control => dft_control%rtp_control
      CALL get_rtp(rtp=rtp, C_mat=C_mat, S_der=S_der, S_inv=S_inv, &
                   SinvH=SinvH, SinvB=SinvB)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)

      NULLIFY (tmp)
      ALLOCATE (tmp)
      CALL dbcsr_create(tmp, template=SinvB(1)%matrix)

      IF (rtp%linear_scaling) THEN
         CALL get_rtp(rtp=rtp, rho_new=rho_new)
      ELSE
         CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao, rho_ao_im=rho_ao_im)
      END IF

      ! If SinvH has an imaginary part (the minus sign is already in SinvH_imag)
      IF (rtp%propagate_complex_ks) CALL get_rtp(rtp=rtp, SinvH_imag=SinvH_imag)

      DO ispin = 1, SIZE(SinvH)
         re = 2*ispin - 1
         im = 2*ispin
         IF (rtp%linear_scaling) THEN
            CALL dbcsr_multiply("N", "N", one, SinvH(ispin)%matrix, rho_new(re)%matrix, zero, tmp, &
                                filter_eps=rtp%filter_eps)
            IF (rtp%propagate_complex_ks) &
               CALL dbcsr_multiply("N", "N", one, SinvH_imag(ispin)%matrix, rho_new(im)%matrix, one, tmp, &
                                   filter_eps=rtp%filter_eps)
            CALL dbcsr_multiply("N", "N", one, SinvB(ispin)%matrix, rho_new(im)%matrix, one, tmp, &
                                filter_eps=rtp%filter_eps)
            CALL compute_forces(force, tmp, S_der, rho_new(im)%matrix, C_mat, kind_of, atom_of_kind)
         ELSE
            CALL dbcsr_multiply("N", "N", one, SinvH(ispin)%matrix, rho_ao(ispin)%matrix, zero, tmp)
            IF (rtp%propagate_complex_ks) &
               CALL dbcsr_multiply("N", "N", one, SinvH_imag(ispin)%matrix, rho_ao_im(ispin)%matrix, one, tmp)
            CALL dbcsr_multiply("N", "N", one, SinvB(ispin)%matrix, rho_ao_im(ispin)%matrix, one, tmp)
            CALL compute_forces(force, tmp, S_der, rho_ao_im(ispin)%matrix, C_mat, kind_of, atom_of_kind)
         END IF
      END DO

      ! recall QS forces, at this point have the other sign.
      DO i = 1, SIZE(force)
         force(i)%ehrenfest(:, :) = -force(i)%ehrenfest(:, :)
      END DO

      CALL dbcsr_deallocate_matrix(tmp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param force ...
!> \param tmp ...
!> \param S_der ...
!> \param rho_im ...
!> \param C_mat ...
!> \param kind_of ...
!> \param atom_of_kind ...
! **************************************************************************************************
   SUBROUTINE compute_forces(force, tmp, S_der, rho_im, C_mat, kind_of, atom_of_kind)
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbcsr_type), POINTER                          :: tmp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: S_der
      TYPE(dbcsr_type), POINTER                          :: rho_im
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: C_mat
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, atom_of_kind

      INTEGER                                            :: col_atom, i, ikind, kind_atom, row_atom
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:), POINTER                    :: block_values, block_values2
      TYPE(dbcsr_iterator_type)                          :: iter

      DO i = 1, 3
         !Calculate the sum over the hadmard product
         !S_der part

         CALL dbcsr_iterator_start(iter, tmp)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
            CALL dbcsr_get_block_p(S_der(i)%matrix, row_atom, col_atom, block_values2, found=found)
            IF (found) THEN
               ikind = kind_of(col_atom)
               kind_atom = atom_of_kind(col_atom)
               !The block_values are in a vector format,
               ! so the dot_product is the sum over all elements of the hamand product, that I need
               force(ikind)%ehrenfest(i, kind_atom) = force(ikind)%ehrenfest(i, kind_atom) + &
                                                      2.0_dp*DOT_PRODUCT(block_values, block_values2)
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)

         !C_mat part

         CALL dbcsr_iterator_start(iter, rho_im)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
            CALL dbcsr_get_block_p(C_mat(i)%matrix, row_atom, col_atom, block_values2, found=found)
            IF (found) THEN
               ikind = kind_of(col_atom)
               kind_atom = atom_of_kind(col_atom)
               !The block_values are in a vector format, so the dot_product is
               ! the sum over all elements of the hamand product, that I need
               force(ikind)%ehrenfest(i, kind_atom) = force(ikind)%ehrenfest(i, kind_atom) + &
                                                      2.0_dp*DOT_PRODUCT(block_values, block_values2)
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO

   END SUBROUTINE compute_forces

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rt_admm_force(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos, mos_admm
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: KS_aux_im, KS_aux_re, matrix_s_aux_fit, &
                                                            matrix_s_aux_fit_vs_orb
      TYPE(rt_prop_type), POINTER                        :: rtp

      CALL get_qs_env(qs_env, &
                      admm_env=admm_env, &
                      rtp=rtp)
      CALL get_admm_env(admm_env, matrix_ks_aux_fit=KS_aux_re, &
                        matrix_ks_aux_fit_im=KS_aux_im, &
                        matrix_s_aux_fit=matrix_s_aux_fit, &
                        matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb)

      CALL get_rtp(rtp=rtp, mos_new=mos, admm_mos=mos_admm)

      ! currently only none option
      CALL rt_admm_forces_none(qs_env, admm_env, KS_aux_re, KS_aux_im, &
                               matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, mos_admm, mos)

   END SUBROUTINE rt_admm_force

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param admm_env ...
!> \param KS_aux_re ...
!> \param KS_aux_im ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param mos_admm ...
!> \param mos ...
! **************************************************************************************************
   SUBROUTINE rt_admm_forces_none(qs_env, admm_env, KS_aux_re, KS_aux_im, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, mos_admm, mos)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: KS_aux_re, KS_aux_im, matrix_s_aux_fit, &
                                                            matrix_s_aux_fit_vs_orb
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_admm, mos

      INTEGER                                            :: im, ispin, jspin, nao, natom, naux, nmo, &
                                                            re
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: admm_force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_struct_type), POINTER                   :: mstruct
      TYPE(cp_fm_type), DIMENSION(2)                     :: tmp_aux_aux, tmp_aux_mo, tmp_aux_mo1, &
                                                            tmp_aux_nao
      TYPE(dbcsr_type), POINTER                          :: matrix_w_q, matrix_w_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit_asymm, sab_aux_fit_vs_orb
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      NULLIFY (sab_aux_fit_asymm, sab_aux_fit_vs_orb, ks_env)

      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL get_admm_env(admm_env, sab_aux_fit_asymm=sab_aux_fit_asymm, &
                        sab_aux_fit_vs_orb=sab_aux_fit_vs_orb)

      ALLOCATE (matrix_w_s)
      CALL dbcsr_create(matrix_w_s, template=matrix_s_aux_fit(1)%matrix, &
                        name='W MATRIX AUX S', matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s, sab_aux_fit_asymm)

      ALLOCATE (matrix_w_q)
      CALL dbcsr_copy(matrix_w_q, matrix_s_aux_fit_vs_orb(1)%matrix, &
                      "W MATRIX AUX Q")

      DO jspin = 1, 2
         CALL cp_fm_create(tmp_aux_aux(jspin), admm_env%work_aux_aux%matrix_struct, name="taa")
         CALL cp_fm_create(tmp_aux_nao(jspin), admm_env%work_aux_orb%matrix_struct, name="tao")
      END DO

      DO ispin = 1, SIZE(KS_aux_re)
         re = 2*ispin - 1; im = 2*ispin
         naux = admm_env%nao_aux_fit; nmo = admm_env%nmo(ispin); nao = admm_env%nao_orb

         mstruct => admm_env%work_aux_nmo(ispin)%matrix_struct
         DO jspin = 1, 2
            CALL cp_fm_create(tmp_aux_mo(jspin), mstruct, name="tam")
            CALL cp_fm_create(tmp_aux_mo1(jspin), mstruct, name="tam")
         END DO

! First calculate H=KS_aux*C~, real part ends on work_aux_aux2, imaginary part ends at work_aux_aux3
         CALL cp_dbcsr_sm_fm_multiply(KS_aux_re(ispin)%matrix, mos_admm(re), tmp_aux_mo(re), nmo, 4.0_dp, 0.0_dp)
         CALL cp_dbcsr_sm_fm_multiply(KS_aux_re(ispin)%matrix, mos_admm(im), tmp_aux_mo(im), nmo, 4.0_dp, 0.0_dp)
         CALL cp_dbcsr_sm_fm_multiply(KS_aux_im(ispin)%matrix, mos_admm(im), tmp_aux_mo(re), nmo, -4.0_dp, 1.0_dp)
         CALL cp_dbcsr_sm_fm_multiply(KS_aux_im(ispin)%matrix, mos_admm(re), tmp_aux_mo(im), nmo, 4.0_dp, 1.0_dp)

! Next step compute S-1*H
         CALL parallel_gemm('N', 'N', naux, nmo, naux, 1.0_dp, admm_env%S_inv, tmp_aux_mo(re), 0.0_dp, tmp_aux_mo1(re))
         CALL parallel_gemm('N', 'N', naux, nmo, naux, 1.0_dp, admm_env%S_inv, tmp_aux_mo(im), 0.0_dp, tmp_aux_mo1(im))

! Here we go on with Ws=S-1*H * C^H (take care of sign of the imaginary part!!!)

         CALL parallel_gemm("N", "T", naux, nao, nmo, -1.0_dp, tmp_aux_mo1(re), mos(re), 0.0_dp, &
                            tmp_aux_nao(re))
         CALL parallel_gemm("N", "T", naux, nao, nmo, -1.0_dp, tmp_aux_mo1(im), mos(im), 1.0_dp, &
                            tmp_aux_nao(re))
         CALL parallel_gemm("N", "T", naux, nao, nmo, 1.0_dp, tmp_aux_mo1(re), mos(im), 0.0_dp, &
                            tmp_aux_nao(im))
         CALL parallel_gemm("N", "T", naux, nao, nmo, -1.0_dp, tmp_aux_mo1(im), mos(re), 1.0_dp, &
                            tmp_aux_nao(im))

! Let's do the final bit  Wq=S-1*H * C^H * A^T
         CALL parallel_gemm('N', 'T', naux, naux, nao, -1.0_dp, tmp_aux_nao(re), admm_env%A, 0.0_dp, tmp_aux_aux(re))
         CALL parallel_gemm('N', 'T', naux, naux, nao, -1.0_dp, tmp_aux_nao(im), admm_env%A, 0.0_dp, tmp_aux_aux(im))

         ! *** copy to sparse matrix
         CALL copy_fm_to_dbcsr(tmp_aux_nao(re), matrix_w_q, keep_sparsity=.TRUE.)

         ! *** copy to sparse matrix
         CALL copy_fm_to_dbcsr(tmp_aux_aux(re), matrix_w_s, keep_sparsity=.TRUE.)

         DO jspin = 1, 2
            CALL cp_fm_release(tmp_aux_mo(jspin))
            CALL cp_fm_release(tmp_aux_mo1(jspin))
         END DO

! *** This can be done in one call w_total = w_alpha + w_beta
         ! allocate force vector
         CALL get_qs_env(qs_env=qs_env, natom=natom)
         ALLOCATE (admm_force(3, natom))
         admm_force = 0.0_dp
         CALL build_overlap_force(ks_env, admm_force, &
                                  basis_type_a="AUX_FIT", basis_type_b="AUX_FIT", &
                                  sab_nl=sab_aux_fit_asymm, matrix_p=matrix_w_s)
         CALL build_overlap_force(ks_env, admm_force, &
                                  basis_type_a="AUX_FIT", basis_type_b="ORB", &
                                  sab_nl=sab_aux_fit_vs_orb, matrix_p=matrix_w_q)
         ! add forces
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                         force=force)
         CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set)
         DEALLOCATE (admm_force)

         ! *** Deallocated weighted density matrices
         CALL dbcsr_deallocate_matrix(matrix_w_s)
         CALL dbcsr_deallocate_matrix(matrix_w_q)
      END DO

      DO jspin = 1, 2
         CALL cp_fm_release(tmp_aux_aux(jspin))
         CALL cp_fm_release(tmp_aux_nao(jspin))
      END DO

   END SUBROUTINE rt_admm_forces_none

END MODULE rt_propagation_forces
