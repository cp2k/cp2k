!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common routines for PAO parametrizations.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_methods
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_complete_redistribute, dbcsr_create, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_diag_blocks, dbcsr_scale, dbcsr_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE dm_ls_scf_qs,                    ONLY: matrix_decluster
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                              ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_methods'

   PUBLIC :: pao_calc_grad_lnv_wrt_U, pao_calc_AB_from_U, pao_calc_grad_lnv_wrt_AB

CONTAINS

! **************************************************************************************************
!> \brief Helper routine, calculates partial derivative dE/dU
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_M_diag the derivate wrt U, matrix uses pao%diag_distribution
! **************************************************************************************************
   SUBROUTINE pao_calc_grad_lnv_wrt_U(qs_env, ls_scf_env, matrix_M_diag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      TYPE(dbcsr_type)                                   :: matrix_M_diag

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_grad_lnv_wrt_U'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: filter_eps
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_M, matrix_Ma, matrix_Mb, matrix_NM
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao

      CALL timeset(routineN, handle)

      ls_mstruct => ls_scf_env%ls_mstruct
      pao => ls_scf_env%pao_env
      filter_eps = ls_scf_env%eps_filter
      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      CALL pao_calc_grad_lnv_wrt_AB(qs_env, ls_scf_env, matrix_Ma, matrix_Mb)

      ! Calculation uses distr. of matrix_s, afterwards we redistribute to pao%diag_distribution.
      CALL dbcsr_create(matrix_M, template=matrix_s(1)%matrix, matrix_type="N")
      CALL dbcsr_reserve_diag_blocks(matrix_M)

      CALL dbcsr_create(matrix_NM, template=ls_mstruct%matrix_A, matrix_type="N")

      CALL dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N_inv, matrix_Ma, &
                          1.0_dp, matrix_NM, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, matrix_Mb, &
                          1.0_dp, matrix_NM, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_NM, pao%matrix_Y, &
                          1.0_dp, matrix_M, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! redistribute using pao%diag_distribution
      CALL dbcsr_create(matrix_M_diag, &
                        name="PAO matrix_M", &
                        matrix_type="N", &
                        dist=pao%diag_distribution, &
                        template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_M_diag)
      CALL dbcsr_complete_redistribute(matrix_M, matrix_M_diag)

      !---------------------------------------------------------------------------
      ! cleanup:
      CALL dbcsr_release(matrix_M)
      CALL dbcsr_release(matrix_Ma)
      CALL dbcsr_release(matrix_Mb)
      CALL dbcsr_release(matrix_NM)

      CALL timestop(handle)
   END SUBROUTINE pao_calc_grad_lnv_wrt_U

! **************************************************************************************************
!> \brief Takes current matrix_X and calculates the matrices A and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_U_diag ...
! **************************************************************************************************
   SUBROUTINE pao_calc_AB_from_U(pao, qs_env, ls_scf_env, matrix_U_diag)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      TYPE(dbcsr_type)                                   :: matrix_U_diag

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_AB_from_U'

      INTEGER                                            :: acol, arow, handle, iatom
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_A, block_B, block_N, block_N_inv, &
                                                            block_U, block_Y
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_U
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      ls_mstruct => ls_scf_env%ls_mstruct

      ! --------------------------------------------------------------------------------------------
      ! sanity check matrix U
      CALL pao_assert_unitary(pao, matrix_U_diag)

      ! --------------------------------------------------------------------------------------------
      ! redistribute matrix_U_diag from diag_distribution to distribution of matrix_s
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_create(matrix_U, matrix_type="N", template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_U)
      CALL dbcsr_complete_redistribute(matrix_U_diag, matrix_U)

      ! --------------------------------------------------------------------------------------------
      ! calculate matrix A and B from matrix U
      ! Multiplying diagonal matrices is a local operation.
      ! To take advantage of this we're using an iterator instead of calling dbcsr_multiply().
!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,ls_mstruct,matrix_U) &
!$OMP PRIVATE(iter,arow,acol,iatom,block_U,block_Y,block_A,block_B,block_N,block_N_inv,found)
      CALL dbcsr_iterator_start(iter, matrix_U)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_U)
         iatom = arow; CPASSERT(arow == acol)

         CALL dbcsr_get_block_p(matrix=pao%matrix_Y, row=iatom, col=iatom, block=block_Y, found=found)
         CPASSERT(ASSOCIATED(block_Y))

         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_A, row=iatom, col=iatom, block=block_A, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N_inv, row=iatom, col=iatom, block=block_N_inv, found=found)
         CPASSERT(ASSOCIATED(block_A) .AND. ASSOCIATED(block_N_inv))

         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_B, row=iatom, col=iatom, block=block_B, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N, row=iatom, col=iatom, block=block_N, found=found)
         CPASSERT(ASSOCIATED(block_B) .AND. ASSOCIATED(block_N))

         block_A = MATMUL(MATMUL(block_N_inv, block_U), block_Y)
         block_B = MATMUL(MATMUL(block_N, block_U), block_Y)

      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_release(matrix_U)

      CALL timestop(handle)
   END SUBROUTINE pao_calc_AB_from_U

! **************************************************************************************************
!> \brief Debugging routine, check unitaryness of U
!> \param pao ...
!> \param matrix_U ...
! **************************************************************************************************
   SUBROUTINE pao_assert_unitary(pao, matrix_U)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(dbcsr_type)                                   :: matrix_U

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_assert_unitary'

      INTEGER                                            :: acol, arow, group_handle, handle, i, &
                                                            iatom, M, N
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri
      REAL(dp)                                           :: delta_max
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_test, tmp1, tmp2
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_comm_type)                                 :: group

      IF (pao%check_unitary_tol < 0.0_dp) RETURN ! no checking

      CALL timeset(routineN, handle)
      delta_max = 0.0_dp

      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri, col_blk_size=blk_sizes_pao)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,matrix_U,blk_sizes_pri,blk_sizes_pao,delta_max) &
!$OMP PRIVATE(iter,arow,acol,iatom,N,M,block_test,tmp1,tmp2)
      CALL dbcsr_iterator_start(iter, matrix_U)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_test)
         iatom = arow; CPASSERT(arow == acol)
         N = blk_sizes_pri(iatom) ! size of primary basis
         M = blk_sizes_pao(iatom) ! size of pao basis

         ! we only need the upper left "PAO-corner" to be unitary
         ALLOCATE (tmp1(N, M), tmp2(M, M))
         tmp1 = block_test(:, 1:M)
         tmp2 = MATMUL(TRANSPOSE(tmp1), tmp1)
         DO i = 1, M
            tmp2(i, i) = tmp2(i, i) - 1.0_dp
         END DO

!$OMP ATOMIC
         delta_max = MAX(delta_max, MAXVAL(ABS(tmp2)))

         DEALLOCATE (tmp1, tmp2)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_get_info(matrix_U, group=group_handle)
      CALL group%set_handle(group_handle)

      CALL group%max(delta_max)
      IF (pao%iw > 0) WRITE (pao%iw, *) 'PAO| checked unitaryness, max delta:', delta_max
      IF (delta_max > pao%check_unitary_tol) &
         CPABORT("Found bad unitaryness:"//cp_to_string(delta_max))

      CALL timestop(handle)
   END SUBROUTINE pao_assert_unitary

! **************************************************************************************************
!> \brief Helper routine, calculates partial derivative dE/dA and dE/dB.
!>        As energy functional serves the definition by LNV (Li, Nunes, Vanderbilt).
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_Ma the derivate wrt A, matrix uses s_matrix-distribution.
!> \param matrix_Mb the derivate wrt B, matrix uses s_matrix-distribution.
! **************************************************************************************************
   SUBROUTINE pao_calc_grad_lnv_wrt_AB(qs_env, ls_scf_env, matrix_Ma, matrix_Mb)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      TYPE(dbcsr_type)                                   :: matrix_Ma, matrix_Mb

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_grad_lnv_wrt_AB'

      INTEGER                                            :: handle, nspin
      INTEGER, DIMENSION(:), POINTER                     :: pao_blk_sizes
      REAL(KIND=dp)                                      :: filter_eps
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, rho_ao
      TYPE(dbcsr_type) :: matrix_HB, matrix_HPS, matrix_M, matrix_M1, matrix_M1_dc, matrix_M2, &
         matrix_M2_dc, matrix_M3, matrix_M3_dc, matrix_PA, matrix_PH, matrix_PHP, matrix_PSP, &
         matrix_SB, matrix_SP
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      ls_mstruct => ls_scf_env%ls_mstruct
      pao => ls_scf_env%pao_env

      CALL get_qs_env(qs_env, &
                      rho=rho, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control)
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      nspin = dft_control%nspins
      filter_eps = ls_scf_env%eps_filter

      CALL dbcsr_get_info(ls_mstruct%matrix_A, col_blk_size=pao_blk_sizes)

      IF (nspin /= 1) CPABORT("open shell not yet implemented")
      !TODO: handle openshell case properly

      ! Notation according to equation (4.6) on page 50 from:
      ! https://dx.doi.org/10.3929%2Fethz-a-010819495

      !---------------------------------------------------------------------------
      ! calculate need products in pao basis
      CALL dbcsr_create(matrix_PH, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), ls_scf_env%matrix_ks(1), &
                          0.0_dp, matrix_PH, filter_eps=filter_eps)

      CALL dbcsr_create(matrix_PHP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_PH, ls_scf_env%matrix_p(1), &
                          0.0_dp, matrix_PHP, filter_eps=filter_eps)

      CALL dbcsr_create(matrix_SP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s, ls_scf_env%matrix_p(1), &
                          0.0_dp, matrix_SP, filter_eps=filter_eps)

      IF (nspin == 1) CALL dbcsr_scale(matrix_SP, 0.5_dp)

      CALL dbcsr_create(matrix_HPS, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "T", 1.0_dp, ls_scf_env%matrix_ks(1), matrix_SP, &
                          0.0_dp, matrix_HPS, filter_eps=filter_eps)

      CALL dbcsr_create(matrix_PSP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), matrix_SP, &
                          0.0_dp, matrix_PSP, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! M1 = dE_lnv / dP_pao
      CALL dbcsr_create(matrix_M1, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL dbcsr_multiply("N", "T", 3.0_dp, ls_scf_env%matrix_ks(1), matrix_SP, &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", 3.0_dp, matrix_SP, ls_scf_env%matrix_ks(1), &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "T", -2.0_dp, matrix_HPS, matrix_SP, &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_SP, matrix_HPS, &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "T", -2.0_dp, matrix_SP, matrix_HPS, &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL dbcsr_create(matrix_M1_dc, &
                        template=matrix_s(1)%matrix, &
                        row_blk_size=pao_blk_sizes, &
                        col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M1_dc, matrix_M1, ls_mstruct)

      !---------------------------------------------------------------------------
      ! M2 = dE_lnv / dH
      CALL dbcsr_create(matrix_M2, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL dbcsr_add(matrix_M2, matrix_PSP, 1.0_dp, 3.0_dp)

      CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_PSP, matrix_SP, &
                          1.0_dp, matrix_M2, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL dbcsr_create(matrix_M2_dc, &
                        template=matrix_s(1)%matrix, &
                        row_blk_size=pao_blk_sizes, &
                        col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M2_dc, matrix_M2, ls_mstruct)

      !---------------------------------------------------------------------------
      ! M3 = dE_lnv / dS
      CALL dbcsr_create(matrix_M3, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL dbcsr_add(matrix_M3, matrix_PHP, 1.0_dp, 3.0_dp)

      CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_PHP, matrix_SP, &
                          1.0_dp, matrix_M3, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "T", -2.0_dp, matrix_PSP, matrix_PH, &
                          1.0_dp, matrix_M3, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL dbcsr_create(matrix_M3_dc, &
                        template=matrix_s(1)%matrix, &
                        row_blk_size=pao_blk_sizes, &
                        col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M3_dc, matrix_M3, ls_mstruct)

      !---------------------------------------------------------------------------
      ! assemble Ma and Mb
      ! matrix_Ma = dE_lnv / dA = P * A * M1
      ! matrix_Mb = dE_lnv / dB = H * B * M2  +  S * B * M3
      CALL dbcsr_create(matrix_Ma, template=ls_mstruct%matrix_A, matrix_type="N")
      CALL dbcsr_reserve_diag_blocks(matrix_Ma)
      CALL dbcsr_create(matrix_Mb, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL dbcsr_reserve_diag_blocks(matrix_Mb)

      !---------------------------------------------------------------------------
      ! combine M1 with matrices from primary basis
      CALL dbcsr_create(matrix_PA, template=ls_mstruct%matrix_A, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, rho_ao(1)%matrix, ls_mstruct%matrix_A, &
                          0.0_dp, matrix_PA, filter_eps=filter_eps)

      ! matrix_Ma = P * A * M1
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_PA, matrix_M1_dc, &
                          0.0_dp, matrix_Ma, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! combine M2 with matrices from primary basis
      CALL dbcsr_create(matrix_HB, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_ks(1)%matrix, ls_mstruct%matrix_B, &
                          0.0_dp, matrix_HB, filter_eps=filter_eps)

      ! matrix_Mb = H * B * M2
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_HB, matrix_M2_dc, &
                          0.0_dp, matrix_Mb, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! combine M3 with matrices from primary basis
      CALL dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B, &
                          0.0_dp, matrix_SB, filter_eps=filter_eps)

      IF (nspin == 1) CALL dbcsr_scale(matrix_SB, 0.5_dp)

      ! matrix_Mb += S * B * M3
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_SB, matrix_M3_dc, &
                          1.0_dp, matrix_Mb, filter_eps=filter_eps)

      IF (nspin == 1) CALL dbcsr_scale(matrix_Ma, 2.0_dp)
      IF (nspin == 1) CALL dbcsr_scale(matrix_Mb, 2.0_dp)

      !---------------------------------------------------------------------------
      ! cleanup: TODO release matrices as early as possible
      CALL dbcsr_release(matrix_PH)
      CALL dbcsr_release(matrix_PHP)
      CALL dbcsr_release(matrix_SP)
      CALL dbcsr_release(matrix_HPS)
      CALL dbcsr_release(matrix_PSP)
      CALL dbcsr_release(matrix_M)
      CALL dbcsr_release(matrix_M1)
      CALL dbcsr_release(matrix_M2)
      CALL dbcsr_release(matrix_M3)
      CALL dbcsr_release(matrix_M1_dc)
      CALL dbcsr_release(matrix_M2_dc)
      CALL dbcsr_release(matrix_M3_dc)
      CALL dbcsr_release(matrix_PA)
      CALL dbcsr_release(matrix_HB)
      CALL dbcsr_release(matrix_SB)

      CALL timestop(handle)
   END SUBROUTINE pao_calc_grad_lnv_wrt_AB

END MODULE pao_param_methods
