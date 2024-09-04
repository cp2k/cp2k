!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_densities
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_scale
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE pw_env_types,                    ONLY: pw_env_get
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
   USE qs_rho_methods,                  ONLY: qs_rho_copy,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE task_list_types,                 ONLY: task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_densities'

   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt_construct_ground_state_orb_density, tddfpt_construct_aux_fit_density

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Compute the ground-state charge density expressed in primary basis set.
!> \param rho_orb_struct      ground-state density in primary basis set
!> \param rho_xc_struct       ground-state density in primary basis set for GAPW_XC
!> \param is_rks_triplets     indicates that the triplet excited states calculation using
!>                            spin-unpolarised molecular orbitals has been requested
!> \param qs_env              Quickstep environment
!> \param sub_env             parallel (sub)group environment
!> \param wfm_rho_orb         work dense matrix with shape [nao x nao] distributed among
!>                            processors of the given parallel group (modified on exit)
!> \par History
!>    * 06.2018 created by splitting the subroutine tddfpt_apply_admm_correction() in two
!>              subroutines tddfpt_construct_ground_state_orb_density() and
!>              tddfpt_construct_aux_fit_density [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_construct_ground_state_orb_density(rho_orb_struct, rho_xc_struct, is_rks_triplets, &
                                                        qs_env, sub_env, wfm_rho_orb)
      TYPE(qs_rho_type), POINTER                         :: rho_orb_struct, rho_xc_struct
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(cp_fm_type), INTENT(IN)                       :: wfm_rho_orb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_construct_ground_state_orb_density'

      INTEGER                                            :: handle, ispin, nao, nspins
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ij_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      nspins = SIZE(sub_env%mos_occ)
      DO ispin = 1, nspins
         CALL cp_fm_get_info(sub_env%mos_occ(ispin), nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      CALL qs_rho_get(rho_orb_struct, rho_ao=rho_ij_ao)
      DO ispin = 1, nspins
         CALL parallel_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, &
                            sub_env%mos_occ(ispin), sub_env%mos_occ(ispin), &
                            0.0_dp, wfm_rho_orb)

         CALL copy_fm_to_dbcsr(wfm_rho_orb, rho_ij_ao(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO

      ! take into account that all MOs are doubly occupied in spin-restricted case
      IF (nspins == 1 .AND. (.NOT. is_rks_triplets)) &
         CALL dbcsr_scale(rho_ij_ao(1)%matrix, 2.0_dp)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (dft_control%qs_control%gapw) THEN
         CALL qs_rho_update_rho(rho_orb_struct, qs_env, &
                                local_rho_set=sub_env%local_rho_set, &
                                pw_env_external=sub_env%pw_env, &
                                task_list_external=sub_env%task_list_orb_soft, &
                                para_env_external=sub_env%para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set=sub_env%local_rho_set)
      ELSEIF (dft_control%qs_control%gapw_xc) THEN
         CALL qs_rho_update_rho(rho_orb_struct, qs_env, &
                                rho_xc_external=rho_xc_struct, &
                                local_rho_set=sub_env%local_rho_set, &
                                pw_env_external=sub_env%pw_env, &
                                task_list_external=sub_env%task_list_orb, &
                                task_list_external_soft=sub_env%task_list_orb_soft, &
                                para_env_external=sub_env%para_env)
         CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL qs_rho_copy(rho_xc_struct, rho_orb_struct, auxbas_pw_pool, nspins)
         CALL prepare_gapw_den(qs_env, local_rho_set=sub_env%local_rho_set, do_rho0=.FALSE.)
      ELSE
         CALL qs_rho_update_rho(rho_orb_struct, qs_env, &
                                pw_env_external=sub_env%pw_env, &
                                task_list_external=sub_env%task_list_orb, &
                                para_env_external=sub_env%para_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_construct_ground_state_orb_density

! **************************************************************************************************
!> \brief Project a charge density expressed in primary basis set into the auxiliary basis set.
!> \param rho_orb_struct      response density in primary basis set
!> \param rho_aux_fit_struct  response density in auxiliary basis set (modified on exit)
!> \param local_rho_set       GAPW density of auxiliary basis set density
!> \param qs_env              Quickstep environment
!> \param sub_env             parallel (sub)group environment
!> \param wfm_rho_orb         work dense matrix with shape [nao x nao] distributed among
!>                            processors of the given parallel group (modified on exit)
!> \param wfm_rho_aux_fit     work dense matrix with shape [nao_aux x nao_aux] distributed among
!>                            processors of the given parallel group (modified on exit)
!> \param wfm_aux_orb         work dense matrix with shape [nao_aux x nao] distributed among
!>                            processors of the given parallel group (modified on exit)
!> \par History
!>    * 03.2017 the subroutine tddfpt_apply_admm_correction() was originally created by splitting
!>              the subroutine tddfpt_apply_hfx() in two parts [Sergey Chulkov]
!>    * 06.2018 created by splitting the subroutine tddfpt_apply_admm_correction() in two subroutines
!>              tddfpt_construct_ground_state_orb_density() and tddfpt_construct_aux_fit_density()
!>              in order to avoid code duplication [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_construct_aux_fit_density(rho_orb_struct, rho_aux_fit_struct, local_rho_set, &
                                               qs_env, sub_env, &
                                               wfm_rho_orb, wfm_rho_aux_fit, wfm_aux_orb)
      TYPE(qs_rho_type), POINTER                         :: rho_orb_struct, rho_aux_fit_struct
      TYPE(local_rho_type), POINTER                      :: local_rho_set
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(cp_fm_type), INTENT(IN)                       :: wfm_rho_orb, wfm_rho_aux_fit, wfm_aux_orb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_construct_aux_fit_density'

      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER                                            :: handle, ispin, nao, nao_aux, nspins
      REAL(kind=dp), DIMENSION(:), POINTER               :: tot_rho_aux_fit_r
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao_aux_fit, rho_ao_orb
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_aux_fit_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_aux_fit_r
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(sub_env%admm_A))

      CALL get_qs_env(qs_env, ks_env=ks_env, admm_env=admm_env)
      CALL qs_rho_get(rho_orb_struct, rho_ao=rho_ao_orb)
      CALL qs_rho_get(rho_aux_fit_struct, rho_ao=rho_ao_aux_fit, rho_g=rho_aux_fit_g, &
                      rho_r=rho_aux_fit_r, tot_rho_r=tot_rho_aux_fit_r)

      nspins = SIZE(rho_ao_orb)

      IF (admm_env%do_gapw) THEN
         basis_type = "AUX_FIT_SOFT"
         task_list => sub_env%task_list_aux_fit_soft
      ELSE
         basis_type = "AUX_FIT"
         task_list => sub_env%task_list_aux_fit
      END IF

      CALL cp_fm_get_info(sub_env%admm_A, nrow_global=nao_aux, ncol_global=nao)
      DO ispin = 1, nspins
         ! TO DO: consider sub_env%admm_A to be a DBCSR matrix
         CALL copy_dbcsr_to_fm(rho_ao_orb(ispin)%matrix, wfm_rho_orb)
         CALL parallel_gemm('N', 'N', nao_aux, nao, nao, 1.0_dp, sub_env%admm_A, &
                            wfm_rho_orb, 0.0_dp, wfm_aux_orb)
         CALL parallel_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, sub_env%admm_A, wfm_aux_orb, &
                            0.0_dp, wfm_rho_aux_fit)
         CALL copy_fm_to_dbcsr(wfm_rho_aux_fit, rho_ao_aux_fit(ispin)%matrix, keep_sparsity=.TRUE.)

         CALL calculate_rho_elec(matrix_p=rho_ao_aux_fit(ispin)%matrix, &
                                 rho=rho_aux_fit_r(ispin), rho_gspace=rho_aux_fit_g(ispin), &
                                 total_rho=tot_rho_aux_fit_r(ispin), ks_env=ks_env, &
                                 soft_valid=.FALSE., basis_type=basis_type, &
                                 pw_env_external=sub_env%pw_env, task_list_external=task_list)
      END DO
      IF (admm_env%do_gapw) THEN
         CALL get_admm_env(qs_env%admm_env, sab_aux_fit=sab_aux_fit)
         CALL calculate_rho_atom_coeff(qs_env, rho_ao_aux_fit, &
                                       rho_atom_set=local_rho_set%rho_atom_set, &
                                       qs_kind_set=admm_env%admm_gapw_env%admm_kind_set, &
                                       oce=admm_env%admm_gapw_env%oce, sab=sab_aux_fit, para_env=sub_env%para_env)
         CALL prepare_gapw_den(qs_env, local_rho_set=local_rho_set, &
                               do_rho0=.FALSE., kind_set_external=admm_env%admm_gapw_env%admm_kind_set)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_construct_aux_fit_density

END MODULE qs_tddfpt2_densities
