!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_kernel_methods
   USE cp_dbcsr_api,                    ONLY: dbcsr_distribution_type
   USE input_section_types,             ONLY: section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fxc,                          ONLY: qs_fxc_prep
   USE qs_kernel_types,                 ONLY: full_kernel_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_copy
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE xc_fxc_kernel,                   ONLY: calc_fxc_kernel
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kernel_methods'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   PUBLIC :: create_kernel_env, create_fxc_kernel

CONTAINS

! **************************************************************************************************
!> \brief Create kernel environment.
!> \param kernel_env       kernel environment (allocated and initialised on exit)
!> \param xc_section       input section which defines an exchange-correlation functional
!> \param is_rks_triplets  indicates that the triplet excited states calculation using
!>                         spin-unpolarised molecular orbitals has been requested
!> \param rho_struct_sub   ground state charge density, if not associated on input, it will be associated on output
!> \param sub_env          parallel group environment
!> \param qs_env ...
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
!>    * 06.2018 the charge density needs to be provided via a dummy argument [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE create_kernel_env(kernel_env, xc_section, is_rks_triplets, rho_struct_sub, &
                                sub_env, qs_env)
      TYPE(full_kernel_env_type), INTENT(inout)          :: kernel_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(qs_rho_type), POINTER                         :: rho_struct_sub
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(qs_environment_type), INTENT(in), POINTER     :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_kernel_env'

      INTEGER                                            :: handle, nspins
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc, weights

      CALL timeset(routineN, handle)

      pw_env => sub_env%pw_env
      dbcsr_dist => sub_env%dbcsr_dist
      sab_orb => sub_env%sab_orb

      nspins = SIZE(sub_env%mos_occ)

      NULLIFY (weights, rho_nlcc)
      weights => sub_env%xcint_weights
      IF (sub_env%is_split .AND. ASSOCIATED(weights)) THEN
         CPABORT("Split communicators and integration weights")
      END IF
      CALL get_qs_env(qs_env, rho_nlcc=rho_nlcc)
      IF (sub_env%is_split .AND. ASSOCIATED(rho_nlcc)) THEN
         CPABORT("Split communicators and NLCC not implemented")
      END IF

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      IF (.NOT. ASSOCIATED(kernel_env%rho0_struct)) THEN
         ALLOCATE (kernel_env%rho0_struct)
         CALL qs_rho_create(kernel_env%rho0_struct)
      END IF
      CALL qs_rho_copy(rho_struct_sub, kernel_env%rho0_struct, auxbas_pw_pool, nspins)

      ! ++ allocate structure for response density
      kernel_env%xc_section => xc_section
      ALLOCATE (kernel_env%xc_rho_set)
      CALL qs_fxc_prep(qs_env, rho_struct_sub, &
                       kernel_env%xc_rho_set, kernel_env%xc_deriv_set, &
                       xc_section, pw_env, is_rks_triplets)

      kernel_env%alpha = 1.0_dp
      kernel_env%beta = 0.0_dp
      ! kernel_env%beta is taken into account in spin-restricted case only
      IF (nspins == 1) THEN
         IF (is_rks_triplets) THEN
            ! K_{triplets} = K_{alpha,alpha} - K_{alpha,beta}
            kernel_env%beta = -1.0_dp
         ELSE
            !                                                 alpha                 beta
            ! K_{singlets} = K_{alpha,alpha} + K_{alpha,beta} = 2 * K_{alpha,alpha} + 0 * K_{alpha,beta},
            ! due to the following relation : K_{alpha,alpha,singlets} == K_{alpha,beta,singlets}
            kernel_env%alpha = 2.0_dp
         END IF
      END IF

      ! finite differences
      kernel_env%deriv2_analytic = section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")
      kernel_env%deriv3_analytic = section_get_lval(xc_section, "3RD_DERIV_ANALYTICAL")

      CALL timestop(handle)

   END SUBROUTINE create_kernel_env

! **************************************************************************************************
!> \brief Create the xc kernel potential for the approximate Fxc kernel model
!> \param rho_struct ...
!> \param fxc_rspace ...
!> \param xc_section ...
!> \param is_rks_triplets ...
!> \param sub_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE create_fxc_kernel(rho_struct, fxc_rspace, xc_section, is_rks_triplets, sub_env, qs_env)
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: fxc_rspace
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL, INTENT(IN)                                :: is_rks_triplets
      TYPE(tddfpt_subgroup_env_type), INTENT(IN)         :: sub_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_fxc_kernel'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: rho_g_valid, tau_r_valid
      REAL(KIND=dp)                                      :: factor
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, tau_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc
      TYPE(section_vals_type), POINTER                   :: xc_kernel

      CALL timeset(routineN, handle)

      pw_env => sub_env%pw_env
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      NULLIFY (rho_r, rho_g, tau_r)
      CALL qs_rho_get(rho_struct, &
                      tau_r_valid=tau_r_valid, &
                      rho_g_valid=rho_g_valid, &
                      rho_r=rho_r, &
                      rho_g=rho_g, &
                      tau_r=tau_r)

      IF (.NOT. tau_r_valid) NULLIFY (tau_r)
      IF (.NOT. rho_g_valid) NULLIFY (rho_g)

      nspins = SIZE(rho_r)

      NULLIFY (rho_nlcc, rho_nlcc_g)
      CALL get_qs_env(qs_env, &
                      rho_nlcc=rho_nlcc, &
                      rho_nlcc_g=rho_nlcc_g)
      ! add the nlcc densities
      IF (ASSOCIATED(rho_nlcc)) THEN
         factor = 1.0_dp
         DO ispin = 1, nspins
            CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
            CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
         END DO
      END IF

      DO ispin = 1, SIZE(fxc_rspace)
         CALL pw_zero(fxc_rspace(ispin))
      END DO

      xc_kernel => section_vals_get_subs_vals(xc_section, "XC_KERNEL")
      CALL calc_fxc_kernel(fxc_rspace, rho_r, rho_g, tau_r, &
                           xc_kernel, is_rks_triplets, auxbas_pw_pool)

      ! remove the nlcc densities (keep stuff in original state)
      IF (ASSOCIATED(rho_nlcc)) THEN
         factor = -1.0_dp
         DO ispin = 1, nspins
            CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
            CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE create_fxc_kernel

END MODULE qs_kernel_methods
