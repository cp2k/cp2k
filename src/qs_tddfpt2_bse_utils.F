!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_bse_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE exstates_types,                  ONLY: excited_energy_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_bse_utils'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC:: tddfpt_apply_bse
   PUBLIC:: zeroth_order_gw

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param Aop_evects ...
!> \param evects ...
!> \param S_evects ...
!> \param gs_mos ...
!> \param matrix_s ...
!> \param matrix_ks ...
! **************************************************************************************************
   SUBROUTINE zeroth_order_gw(qs_env, Aop_evects, evects, S_evects, gs_mos, matrix_s, matrix_ks)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(INOUT)   :: Aop_evects
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: evects, S_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_type), INTENT(in), POINTER              :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_ks

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'zeroth_order_gw'

      INTEGER                                            :: handle, i, ispin, ivect, nactive, nao, &
                                                            nmo, nspins, nvects, occ, virt
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: gw_occ, gw_virt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, matrix_struct
      TYPE(cp_fm_type)                                   :: fms, hevec, matrixtmp, matrixtmp2, &
                                                            matrixtmp3, Sweighted_vect, &
                                                            weighted_vect
      TYPE(dbcsr_type)                                   :: matrixf
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      NULLIFY (ex_env, sab_orb)
      CALL get_qs_env(qs_env, exstate_env=ex_env, sab_orb=sab_orb)

      nspins = SIZE(matrix_ks, 1)
      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      DO ispin = 1, nspins

         CALL dbcsr_create(matrixf, template=matrix_s)
         nmo = SIZE(ex_env%gw_eigen)
         CALL cp_fm_get_info(matrix=evects(ispin, 1), matrix_struct=matrix_struct, &
                             nrow_global=nao, ncol_global=nactive)
         NULLIFY (blacs_env, para_env)
         CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

         occ = SIZE(gs_mos(ispin)%evals_occ)
         nmo = SIZE(ex_env%gw_eigen)
         virt = SIZE(gs_mos(ispin)%evals_virt)
         NULLIFY (fmstruct)
         CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, &
                                  context=blacs_env, nrow_global=virt, ncol_global=virt)
         CALL cp_fm_create(matrixtmp, fmstruct)
         CALL cp_fm_struct_release(fmstruct)

         NULLIFY (fmstruct)
         CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, &
                                  context=blacs_env, nrow_global=virt, ncol_global=nao)
         CALL cp_fm_create(matrixtmp2, fmstruct)
         CALL cp_fm_struct_release(fmstruct)

         NULLIFY (fmstruct)
         CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, &
                                  context=blacs_env, nrow_global=nao, ncol_global=nao)
         CALL cp_fm_create(matrixtmp3, fmstruct)
         CALL cp_fm_create(fms, fmstruct)
         CALL cp_fm_struct_release(fmstruct)
         CALL cp_dbcsr_alloc_block_from_nbl(matrixf, sab_orb)

!--add virt eigenvalues
         CALL dbcsr_set(matrixf, 0.0_dp)
         CALL cp_fm_create(weighted_vect, gs_mos(ispin)%mos_virt%matrix_struct)
         CALL cp_fm_create(Sweighted_vect, gs_mos(ispin)%mos_virt%matrix_struct)
         CALL cp_fm_to_fm(gs_mos(ispin)%mos_virt, weighted_vect)
         CALL copy_dbcsr_to_fm(matrix_s, fms)

         ALLOCATE (gw_virt(virt))
         ALLOCATE (gw_occ(occ))
         gw_virt(1:virt) = ex_env%gw_eigen(occ + 1:nmo)
         gw_occ(1:occ) = ex_env%gw_eigen(1:occ)

         CALL cp_fm_set_all(matrixtmp, 0.0_dp)
         DO i = 1, virt
            CALL cp_fm_set_element(matrixtmp, i, i, gw_virt(i))
         END DO
         DEALLOCATE (gw_virt)
         CALL parallel_gemm('N', 'N', nao, virt, nao, 1.0_dp, fms, weighted_vect, 0.0_dp, Sweighted_vect)
         CALL parallel_gemm('N', 'T', virt, nao, virt, 1.0_dp, matrixtmp, Sweighted_vect, 0.0_dp, matrixtmp2)
         CALL parallel_gemm('N', 'N', nao, nao, virt, 1.0_dp, Sweighted_vect, matrixtmp2, 0.0_dp, matrixtmp3)
         CALL copy_fm_to_dbcsr(matrixtmp3, matrixf)

         CALL cp_fm_release(weighted_vect)
         CALL cp_fm_release(Sweighted_vect)
         CALL cp_fm_release(fmS)
!--add occ eigenvalues
         CALL cp_fm_get_info(matrix=evects(ispin, 1), matrix_struct=matrix_struct, &
                             nrow_global=nao, ncol_global=nactive)
         CALL cp_fm_create(hevec, matrix_struct)

         DO ivect = 1, nvects
            CALL cp_dbcsr_sm_fm_multiply(matrixf, evects(ispin, ivect), &
                                         Aop_evects(ispin, ivect), ncol=nactive, &
                                         alpha=1.0_dp, beta=1.0_dp)

            CALL cp_fm_to_fm(S_evects(ispin, ivect), hevec)
            CALL cp_fm_column_scale(hevec, gw_occ)

            CALL cp_fm_scale_and_add(1.0_dp, Aop_evects(ispin, ivect), -1.0_dp, hevec)
         END DO !ivect
         DEALLOCATE (gw_occ)

         CALL cp_fm_release(matrixtmp)
         CALL cp_fm_release(matrixtmp2)
         CALL cp_fm_release(matrixtmp3)

         CALL dbcsr_release(matrixf)
         CALL cp_fm_release(hevec)
      END DO !ispin

      virt = SIZE(Aop_evects, 2)
      CALL timestop(handle)

   END SUBROUTINE zeroth_order_gw

! **************************************************************************************************
!> \brief Update action of TDDFPT operator on trial vectors by adding BSE W term.
!> \param Aop_evects ...
!> \param evects ...
!> \param gs_mos ...
!> \param qs_env ...
!> \note Based on the subroutine tddfpt_apply_hfx() which was originally created by
!>       Mohamed Fawzi on 10.2002.
! **************************************************************************************************
   SUBROUTINE tddfpt_apply_bse(Aop_evects, evects, gs_mos, qs_env)

      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(INOUT)   :: Aop_evects
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(IN)      :: evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddfpt_apply_bse'

      INTEGER :: a_nao_col, a_virt_col, b_nao_col, c_virt_col, handle, i_occ_row, i_row_global, &
         ii, ispin, ivect, j_col_global, j_occ_row, jj, k_occ_col, mu_col_global, nao, ncol_block, &
         ncol_local, nrow_block, nrow_local, nspins, nvects, nvirt
      INTEGER, DIMENSION(2)                              :: nactive
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, matrix_struct
      TYPE(cp_fm_type)                                   :: CSvirt, fms, WXaoao, WXmat2, WXvirtao
      TYPE(cp_fm_type), POINTER                          :: bse_w_matrix_MO
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(excited_energy_type), POINTER                 :: ex_env
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)
      IF (nspins > 1) THEN
         alpha = 1.0_dp
      ELSE
         alpha = 2.0_dp
      END IF
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)
      DO ispin = 1, nspins
         CALL cp_fm_get_info(evects(ispin, 1), ncol_global=nactive(ispin))
      END DO

      NULLIFY (ex_env, para_env, blacs_env, matrix_s)
      CALL get_qs_env(qs_env, exstate_env=ex_env, para_env=para_env, blacs_env=blacs_env, &
                      matrix_s=matrix_s)

      CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, &
                               context=blacs_env, nrow_global=nao, ncol_global=nao)
      CALL cp_fm_create(fms, fmstruct)
      CALL cp_fm_struct_release(fmstruct)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, fms)

      NULLIFY (bse_w_matrix_MO)
      bse_w_matrix_MO => ex_env%bse_w_matrix_MO(1, 1)

      DO ivect = 1, nvects
         DO ispin = 1, nspins
            NULLIFY (matrix_struct, fmstruct)
            CALL cp_fm_get_info(matrix=evects(ispin, 1), matrix_struct=matrix_struct, &
                                nrow_global=nao, ncol_global=nactive(ispin))
            nvirt = SIZE(gs_mos(ispin)%evals_virt)

            CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, &
                                     context=blacs_env, nrow_global=nvirt, ncol_global=nao)
            CALL cp_fm_create(CSvirt, fmstruct)
            CALL cp_fm_struct_release(fmstruct)

            CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, &
                                     context=blacs_env, nrow_global=nactive(ispin)*nactive(ispin), &
                                     ncol_global=nvirt*nao)
            CALL cp_fm_create(WXvirtao, fmstruct)
            CALL cp_fm_struct_release(fmstruct)
            CALL cp_fm_set_all(WXvirtao, 0.0_dp)

            CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, &
                                     context=blacs_env, nrow_global=nactive(ispin)*nactive(ispin), &
                                     ncol_global=nao*nao)
            CALL cp_fm_create(WXaoao, fmstruct)
            CALL cp_fm_struct_release(fmstruct)
            CALL cp_fm_set_all(WXaoao, 0.0_dp)

            CALL parallel_gemm('T', 'N', nvirt, nao, nao, 1.0_dp, gs_mos(ispin)%mos_virt, fms, 0.0_dp, CSvirt)
            NULLIFY (row_indices, col_indices)
            CALL cp_fm_get_info(matrix=WXvirtao, nrow_local=nrow_local, ncol_local=ncol_local, &
                                row_indices=row_indices, col_indices=col_indices, &
                                nrow_block=nrow_block, ncol_block=ncol_block)

            CALL cp_fm_set_all(WXvirtao, 0.0_dp)
            DO ii = 1, nrow_local
               i_row_global = row_indices(ii)
               DO jj = 1, ncol_local
                  j_col_global = col_indices(jj)

                  i_occ_row = (i_row_global - 1)/nactive(ispin) + 1
                  j_occ_row = MOD(i_row_global - 1, nactive(ispin)) + 1

                  a_virt_col = (j_col_global - 1)/nao + 1
                  b_nao_col = MOD(j_col_global - 1, nao) + 1

                  DO c_virt_col = 1, nvirt
                     mu_col_global = (a_virt_col - 1)*nvirt + c_virt_col

                     WXvirtao%local_data(i_row_global, j_col_global) = WXvirtao%local_data(i_row_global, j_col_global) + &
                                    bse_w_matrix_MO%local_data(i_row_global, mu_col_global)*CSvirt%local_data(c_virt_col, b_nao_col)

                  END DO
               END DO
            END DO

            NULLIFY (row_indices, col_indices) ! redefine indices
            CALL cp_fm_get_info(matrix=WXaoao, nrow_local=nrow_local, ncol_local=ncol_local, &
                                row_indices=row_indices, col_indices=col_indices, &
                                nrow_block=nrow_block, ncol_block=ncol_block)

            CALL cp_fm_set_all(WXaoao, 0.0_dp)
            DO ii = 1, nrow_local
               i_row_global = row_indices(ii)
               DO jj = 1, ncol_local
                  j_col_global = col_indices(jj)

                  i_occ_row = (i_row_global - 1)/nactive(ispin) + 1
                  j_occ_row = MOD(i_row_global - 1, nactive(ispin)) + 1

                  a_nao_col = (j_col_global - 1)/nao + 1
                  b_nao_col = MOD(j_col_global - 1, nao) + 1

                  DO k_occ_col = 1, nvirt
                     mu_col_global = (k_occ_col - 1)*nao + a_nao_col

                     WXaoao%local_data(i_row_global, j_col_global) = WXaoao%local_data(i_row_global, j_col_global) + &
                                            WXvirtao%local_data(i_row_global, mu_col_global)*CSvirt%local_data(k_occ_col, b_nao_col)

                  END DO
               END DO
            END DO

            CALL cp_fm_release(WXvirtao)
            CALL cp_fm_release(CSvirt)

            CALL cp_fm_struct_create(fmstruct=fmstruct, para_env=para_env, &
                                     context=blacs_env, nrow_global=nao, ncol_global=nactive(ispin))
            CALL cp_fm_create(WXmat2, fmstruct)
            CALL cp_fm_struct_release(fmstruct)
            CALL cp_fm_set_all(WXmat2, 0.0_dp)

            DO ii = 1, nrow_local
               i_row_global = row_indices(ii)
               DO jj = 1, ncol_local
                  j_col_global = col_indices(jj)

                  i_occ_row = (i_row_global - 1)/nactive(ispin) + 1
                  j_occ_row = MOD(i_row_global - 1, nactive(ispin)) + 1
                  a_nao_col = (j_col_global - 1)/nao + 1
                  b_nao_col = MOD(j_col_global - 1, nao) + 1
                  IF (a_nao_col .EQ. b_nao_col) THEN
                     WXmat2%local_data(b_nao_col, i_occ_row) = WXmat2%local_data(b_nao_col, i_occ_row) + &
                                 WXaoao%local_data(i_row_global, j_col_global)*evects(ispin, ivect)%local_data(a_nao_col, j_occ_row)
                  END IF
                  IF (a_nao_col .NE. b_nao_col) THEN
                     WXmat2%local_data(a_nao_col, i_occ_row) = WXmat2%local_data(a_nao_col, i_occ_row) + &
                                 WXaoao%local_data(i_row_global, j_col_global)*evects(ispin, ivect)%local_data(b_nao_col, j_occ_row)
                  END IF
               END DO
            END DO

            CALL cp_fm_release(WXaoao)

            CALL cp_fm_scale_and_add(1.0_dp, Aop_evects(ispin, ivect), -1.0_dp, WXmat2)

            CALL cp_fm_release(WXmat2)

         END DO! ispin
      END DO !ivect

      CALL cp_fm_release(fms)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_apply_bse
END MODULE qs_tddfpt2_bse_utils
