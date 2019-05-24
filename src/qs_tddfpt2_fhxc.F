!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_fhxc
   USE cp_control_types,                ONLY: stda_control_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_p_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE kinds,                           ONLY: dp
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_zero
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density
   USE qs_tddfpt2_operators,            ONLY: tddfpt_apply_coulomb,&
                                              tddfpt_apply_xc
   USE qs_tddfpt2_stda_types,           ONLY: stda_env_type
   USE qs_tddfpt2_stda_utils,           ONLY: stda_calculate_kernel
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: full_kernel_env_type,&
                                              tddfpt_work_matrices
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
!> \param qs_env                Quickstep environment
!> \param kernel_env            kernel environment
!> \param kernel_env_admm_aux   kernel environment for ADMM correction
!> \param sub_env               parallel (sub)group environment
!> \param work_matrices         collection of work matrices (modified on exit)
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
!>    * 03.2017 refactored [Sergey Chulkov]
!>    * 04.2019 refactored [JHU]
! **************************************************************************************************
   SUBROUTINE fhxc_kernel(Aop_evects, evects, is_rks_triplets, &
                          do_hfx, qs_env, kernel_env, kernel_env_admm_aux, &
                          sub_env, work_matrices)
      TYPE(cp_fm_p_type), DIMENSION(:, :)                :: Aop_evects, evects
      LOGICAL, INTENT(in)                                :: is_rks_triplets, do_hfx
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(full_kernel_env_type), POINTER                :: kernel_env, kernel_env_admm_aux
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fhxc_kernel', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nao, nspins, nvects
      INTEGER, DIMENSION(maxspins)                       :: nactive
      LOGICAL                                            :: do_admm
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ia_ao, rho_ia_ao_aux_fit
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ia_g, rho_ia_g_aux_fit, rho_ia_r, &
                                                            rho_ia_r_aux_fit, tau_ia_r, &
                                                            tau_ia_r_aux_fit

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)
      do_admm = ASSOCIATED(sub_env%admm_A)

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
            CALL tddfpt_apply_xc(A_ia_rspace=work_matrices%A_ia_rspace_sub, &
                                 kernel_env=kernel_env_admm_aux, &
                                 rho_ia_struct=work_matrices%rho_aux_fit_struct_sub, &
                                 is_rks_triplets=is_rks_triplets, pw_env=sub_env%pw_env, &
                                 work_v_xc=work_matrices%wpw_rspace_sub)
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'stda_kernel', routineP = moduleN//':'//routineN

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
