!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_fhxc
   USE cp_control_types,                ONLY: stda_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_get_info, &
        dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_symmetric
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_kernel_types,                 ONLY: full_kernel_env_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density
   USE qs_tddfpt2_operators,            ONLY: tddfpt_apply_coulomb,&
                                              tddfpt_apply_xc
   USE qs_tddfpt2_stda_types,           ONLY: stda_env_type
   USE qs_tddfpt2_stda_utils,           ONLY: stda_calculate_kernel
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_work_matrices
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_fhxc'

   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: fhxc_kernel, stda_kernel

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Compute action matrix-vector products with the FHxc Kernel
!> \param Aop_evects            action of TDDFPT operator on trial vectors (modified on exit)
!> \param evects                TDDFPT trial vectors
!> \param is_rks_triplets       indicates that a triplet excited states calculation using
!>                              spin-unpolarised molecular orbitals has been requested
!> \param do_hfx                flag that activates computation of exact-exchange terms
!> \param do_admm ...
!> \param qs_env                Quickstep environment
!> \param kernel_env            kernel environment
!> \param kernel_env_admm_aux   kernel environment for ADMM correction
!> \param sub_env               parallel (sub)group environment
!> \param work_matrices         collection of work matrices (modified on exit)
!> \param admm_symm             use symmetric definition of ADMM kernel correction
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
!>    * 03.2017 refactored [Sergey Chulkov]
!>    * 04.2019 refactored [JHU]
! **************************************************************************************************
   SUBROUTINE fhxc_kernel(Aop_evects, evects, is_rks_triplets, &
                          do_hfx, do_admm, qs_env, kernel_env, kernel_env_admm_aux, &
                          sub_env, work_matrices, admm_symm)
      TYPE(cp_fm_p_type), DIMENSION(:, :)                :: Aop_evects, evects
      LOGICAL, INTENT(in)                                :: is_rks_triplets, do_hfx, do_admm
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(full_kernel_env_type), POINTER                :: kernel_env, kernel_env_admm_aux
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices
      LOGICAL, INTENT(in)                                :: admm_symm

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fhxc_kernel'

      INTEGER                                            :: handle, ispin, ivect, nao, nao_aux, &
                                                            nspins, nvects
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes
      INTEGER, DIMENSION(maxspins)                       :: nactive
      TYPE(cp_fm_type), POINTER                          :: work_aux_orb, work_orb_orb
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: A_xc_munu_sub, rho_ia_ao, &
                                                            rho_ia_ao_aux_fit
      TYPE(dbcsr_type), POINTER                          :: dbwork
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_g_aux_fit, rho_ia_r, &
                                                            rho_ia_r_aux_fit, tau_ia_r, &
                                                            tau_ia_r_aux_fit, V_rspace_sub
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)
      IF (do_admm) THEN
         CPASSERT(do_hfx)
         CPASSERT(ASSOCIATED(sub_env%admm_A))
      END IF

      CALL cp_fm_get_info(evects(1, 1)%matrix, nrow_global=nao)
      DO ispin = 1, nspins
         CALL cp_fm_get_info(evects(ispin, 1)%matrix, ncol_global=nactive(ispin))
      END DO

      CALL qs_rho_get(work_matrices%rho_orb_struct_sub, rho_ao=rho_ia_ao, &
                      rho_g=rho_ia_g, rho_r=rho_ia_r, tau_r=tau_ia_r)
      IF (do_hfx .AND. do_admm) THEN
         CALL qs_rho_get(work_matrices%rho_aux_fit_struct_sub, &
                         rho_ao=rho_ia_ao_aux_fit, rho_g=rho_ia_g_aux_fit, &
                         rho_r=rho_ia_r_aux_fit, tau_r=tau_ia_r_aux_fit)
      END IF

      DO ivect = 1, nvects
         IF (ALLOCATED(work_matrices%evects_sub)) THEN
            IF (ASSOCIATED(work_matrices%evects_sub(1, ivect)%matrix)) THEN
               DO ispin = 1, nspins
                  ! Can this be done using cp_dbcsr_plus_fm_fm_t ?
                  CALL cp_gemm('N', 'T', nao, nao, nactive(ispin), &
                               0.5_dp, sub_env%mos_occ(ispin)%matrix, &
                               work_matrices%evects_sub(ispin, ivect)%matrix, &
                               0.0_dp, work_matrices%rho_ao_orb_fm_sub)
                  CALL cp_gemm('N', 'T', nao, nao, nactive(ispin), &
                               0.5_dp, work_matrices%evects_sub(ispin, ivect)%matrix, &
                               sub_env%mos_occ(ispin)%matrix, &
                               1.0_dp, work_matrices%rho_ao_orb_fm_sub)

                  CALL copy_fm_to_dbcsr(work_matrices%rho_ao_orb_fm_sub, &
                                        rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
               END DO
            ELSE
               ! skip trial vectors which are assigned to different parallel groups
               CYCLE
            END IF
         ELSE
            DO ispin = 1, nspins
               ! Can this be done using cp_dbcsr_plus_fm_fm_t ?
               CALL cp_gemm('N', 'T', nao, nao, nactive(ispin), 0.5_dp, sub_env%mos_occ(ispin)%matrix, &
                            evects(ispin, ivect)%matrix, 0.0_dp, work_matrices%rho_ao_orb_fm_sub)
               CALL cp_gemm('N', 'T', nao, nao, nactive(ispin), 0.5_dp, evects(ispin, ivect)%matrix, &
                            sub_env%mos_occ(ispin)%matrix, 1.0_dp, work_matrices%rho_ao_orb_fm_sub)

               CALL copy_fm_to_dbcsr(work_matrices%rho_ao_orb_fm_sub, &
                                     rho_ia_ao(ispin)%matrix, keep_sparsity=.TRUE.)
            END DO
         END IF

         CALL qs_rho_update_rho(work_matrices%rho_orb_struct_sub, qs_env, &
                                pw_env_external=sub_env%pw_env, task_list_external=sub_env%task_list_orb)

         DO ispin = 1, nspins
            CALL dbcsr_set(work_matrices%A_ia_munu_sub(ispin)%matrix, 0.0_dp)
         END DO

         ! electron-hole exchange-correlation interaction
         DO ispin = 1, nspins
            CALL pw_zero(work_matrices%A_ia_rspace_sub(ispin)%pw)
         END DO

         ! C_x d^{2}E_{x}^{DFT}[\rho] / d\rho^2
         ! + C_{HF} d^{2}E_{x, ADMM}^{DFT}[\rho] / d\rho^2 in case of ADMM calculation
         CALL tddfpt_apply_xc(A_ia_rspace=work_matrices%A_ia_rspace_sub, kernel_env=kernel_env, &
                              rho_ia_struct=work_matrices%rho_orb_struct_sub, &
                              is_rks_triplets=is_rks_triplets, pw_env=sub_env%pw_env, &
                              work_v_xc=work_matrices%wpw_rspace_sub)

         ! ADMM correction
         IF (do_admm) THEN
            CALL tddfpt_construct_aux_fit_density(rho_orb_struct=work_matrices%rho_orb_struct_sub, &
                                                  rho_aux_fit_struct=work_matrices%rho_aux_fit_struct_sub, &
                                                  qs_env=qs_env, sub_env=sub_env, &
                                                  wfm_rho_orb=work_matrices%rho_ao_orb_fm_sub, &
                                                  wfm_rho_aux_fit=work_matrices%rho_ao_aux_fit_fm_sub, &
                                                  wfm_aux_orb=work_matrices%wfm_aux_orb_sub)

            ! - C_{HF} d^{2}E_{x, ADMM}^{DFT}[\hat{\rho}] / d\hat{\rho}^2
            IF (admm_symm) THEN
               CALL dbcsr_get_info(rho_ia_ao_aux_fit(1)%matrix, row_blk_size=blk_sizes)
               ALLOCATE (A_xc_munu_sub(nspins))
               DO ispin = 1, nspins
                  ALLOCATE (A_xc_munu_sub(ispin)%matrix)
                  CALL dbcsr_create(matrix=A_xc_munu_sub(ispin)%matrix, name="ADMM_XC", &
                                    dist=sub_env%dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                                    row_blk_size=blk_sizes, col_blk_size=blk_sizes, nze=0)
                  CALL cp_dbcsr_alloc_block_from_nbl(A_xc_munu_sub(ispin)%matrix, sub_env%sab_aux_fit)
                  CALL dbcsr_set(A_xc_munu_sub(ispin)%matrix, 0.0_dp)
               END DO

               CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
               ALLOCATE (V_rspace_sub(nspins))
               DO ispin = 1, nspins
                  NULLIFY (V_rspace_sub(ispin)%pw)
                  CALL pw_pool_create_pw(auxbas_pw_pool, V_rspace_sub(ispin)%pw, &
                                         use_data=REALDATA3D, in_space=REALSPACE)
                  CALL pw_zero(V_rspace_sub(ispin)%pw)
               END DO

               CALL tddfpt_apply_xc(A_ia_rspace=V_rspace_sub, &
                                    kernel_env=kernel_env_admm_aux, &
                                    rho_ia_struct=work_matrices%rho_aux_fit_struct_sub, &
                                    is_rks_triplets=is_rks_triplets, pw_env=sub_env%pw_env, &
                                    work_v_xc=work_matrices%wpw_rspace_sub)
               DO ispin = 1, nspins
                  CALL pw_scale(V_rspace_sub(ispin)%pw, V_rspace_sub(ispin)%pw%pw_grid%dvol)
                  CALL integrate_v_rspace(v_rspace=V_rspace_sub(ispin), &
                                          hmat=A_xc_munu_sub(ispin), &
                                          qs_env=qs_env, calculate_forces=.FALSE., &
                                          pw_env_external=sub_env%pw_env, &
                                          basis_type="AUX_FIT", &
                                          task_list_external=sub_env%task_list_aux_fit)
               END DO
               ALLOCATE (dbwork)
               CALL dbcsr_create(dbwork, template=work_matrices%A_ia_munu_sub(1)%matrix)
               CALL cp_fm_create(work_aux_orb, &
                                 matrix_struct=work_matrices%wfm_aux_orb_sub%matrix_struct)
               CALL cp_fm_create(work_orb_orb, &
                                 matrix_struct=work_matrices%rho_ao_orb_fm_sub%matrix_struct)
               CALL cp_fm_get_info(work_aux_orb, nrow_global=nao_aux, ncol_global=nao)
               DO ispin = 1, nspins
                  CALL cp_dbcsr_sm_fm_multiply(A_xc_munu_sub(ispin)%matrix, sub_env%admm_A, &
                                               work_aux_orb, nao)
                  CALL cp_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, sub_env%admm_A, &
                               work_aux_orb, 0.0_dp, work_orb_orb)
                  CALL dbcsr_copy(dbwork, work_matrices%A_ia_munu_sub(1)%matrix)
                  CALL dbcsr_set(dbwork, 0.0_dp)
                  CALL copy_fm_to_dbcsr(work_orb_orb, dbwork, keep_sparsity=.TRUE.)
                  CALL dbcsr_add(work_matrices%A_ia_munu_sub(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
               END DO
               CALL dbcsr_release(dbwork)
               DEALLOCATE (dbwork)
               DO ispin = 1, nspins
                  CALL pw_pool_give_back_pw(auxbas_pw_pool, V_rspace_sub(ispin)%pw)
               END DO
               DEALLOCATE (V_rspace_sub)
               CALL cp_fm_release(work_aux_orb)
               CALL cp_fm_release(work_orb_orb)
               DO ispin = 1, nspins
                  CALL dbcsr_deallocate_matrix(A_xc_munu_sub(ispin)%matrix)
               END DO
               DEALLOCATE (A_xc_munu_sub)
            ELSE
               CALL tddfpt_apply_xc(A_ia_rspace=work_matrices%A_ia_rspace_sub, &
                                    kernel_env=kernel_env_admm_aux, &
                                    rho_ia_struct=work_matrices%rho_aux_fit_struct_sub, &
                                    is_rks_triplets=is_rks_triplets, pw_env=sub_env%pw_env, &
                                    work_v_xc=work_matrices%wpw_rspace_sub)
            END IF
         END IF

         ! electron-hole Coulomb interaction
         IF (.NOT. is_rks_triplets) THEN
            ! a sum J_i{alpha}a{alpha}_munu + J_i{beta}a{beta}_munu can be computed by solving
            ! the Poisson equation for combined density (rho_{ia,alpha} + rho_{ia,beta}) .
            ! The following action will destroy reciprocal-space grid in spin-unrestricted case.
            DO ispin = 2, nspins
               CALL pw_axpy(rho_ia_g(ispin)%pw, rho_ia_g(1)%pw)
            END DO

            CALL tddfpt_apply_coulomb(A_ia_rspace=work_matrices%A_ia_rspace_sub, &
                                      rho_ia_g=rho_ia_g(1)%pw, pw_env=sub_env%pw_env, &
                                      work_v_gspace=work_matrices%wpw_gspace_sub(1), &
                                      work_v_rspace=work_matrices%wpw_rspace_sub(1))
         END IF

         ! convert from the plane-wave representation into the Gaussian basis set representation
         DO ispin = 1, nspins
            CALL pw_scale(work_matrices%A_ia_rspace_sub(ispin)%pw, &
                          work_matrices%A_ia_rspace_sub(ispin)%pw%pw_grid%dvol)

            CALL integrate_v_rspace(v_rspace=work_matrices%A_ia_rspace_sub(ispin), &
                                    hmat=work_matrices%A_ia_munu_sub(ispin), &
                                    qs_env=qs_env, calculate_forces=.FALSE., gapw=.FALSE., &
                                    pw_env_external=sub_env%pw_env, task_list_external=sub_env%task_list_orb)
         END DO

         IF (ALLOCATED(work_matrices%evects_sub)) THEN
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(work_matrices%A_ia_munu_sub(ispin)%matrix, &
                                            sub_env%mos_occ(ispin)%matrix, &
                                            work_matrices%Aop_evects_sub(ispin, ivect)%matrix, &
                                            ncol=nactive(ispin), alpha=1.0_dp, beta=0.0_dp)
            END DO
         ELSE
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(work_matrices%A_ia_munu_sub(ispin)%matrix, &
                                            sub_env%mos_occ(ispin)%matrix, &
                                            Aop_evects(ispin, ivect)%matrix, &
                                            ncol=nactive(ispin), alpha=1.0_dp, beta=0.0_dp)
            END DO
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE fhxc_kernel

! **************************************************************************************************
!> \brief Compute action matrix-vector products with the sTDA Kernel
!> \param Aop_evects            action of TDDFPT operator on trial vectors (modified on exit)
!> \param evects                TDDFPT trial vectors
!> \param is_rks_triplets       indicates that a triplet excited states calculation using
!>                              spin-unpolarised molecular orbitals has been requested
!> \param qs_env                Quickstep environment
!> \param stda_control          control parameters for sTDA kernel
!> \param stda_env ...
!> \param sub_env               parallel (sub)group environment
!> \param work_matrices         collection of work matrices (modified on exit)
!> \par History
!>    * 04.2019 initial version [JHU]
! **************************************************************************************************
   SUBROUTINE stda_kernel(Aop_evects, evects, is_rks_triplets, &
                          qs_env, stda_control, stda_env, &
                          sub_env, work_matrices)

      TYPE(cp_fm_p_type), DIMENSION(:, :)                :: Aop_evects, evects
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_control_type)                            :: stda_control
      TYPE(stda_env_type)                                :: stda_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'stda_kernel'

      INTEGER                                            :: handle, ivect, nvects

      CALL timeset(routineN, handle)

      nvects = SIZE(evects, 2)

      DO ivect = 1, nvects
         IF (ALLOCATED(work_matrices%evects_sub)) THEN
            IF (ASSOCIATED(work_matrices%evects_sub(1, ivect)%matrix)) THEN
               CALL stda_calculate_kernel(qs_env, stda_control, stda_env, sub_env, work_matrices, &
                                          is_rks_triplets, work_matrices%evects_sub(:, ivect), &
                                          work_matrices%Aop_evects_sub(:, ivect))
            ELSE
               ! skip trial vectors which are assigned to different parallel groups
               CYCLE
            END IF
         ELSE
            CALL stda_calculate_kernel(qs_env, stda_control, stda_env, sub_env, work_matrices, &
                                       is_rks_triplets, evects(:, ivect), Aop_evects(:, ivect))
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE stda_kernel

! **************************************************************************************************

END MODULE qs_tddfpt2_fhxc
