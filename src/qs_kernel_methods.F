!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_kernel_methods
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_distribution_type,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_get_rval,&
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
   USE qs_kernel_types,                 ONLY: full_kernel_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_dbcsr_create_by_dist,&
                                              tddfpt_subgroup_env_type
   USE xc,                              ONLY: xc_prep_2nd_deriv
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_fxc_kernel,                   ONLY: calc_fxc_kernel
   USE xc_rho_set_types,                ONLY: xc_rho_set_create
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
!> \param lsd_singlets ...
!> \param do_excitations ...
!> \param sub_env          parallel group environment
!> \param qs_env ...
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
!>    * 06.2018 the charge density needs to be provided via a dummy argument [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE create_kernel_env(kernel_env, xc_section, is_rks_triplets, rho_struct_sub, &
                                lsd_singlets, do_excitations, sub_env, qs_env)
      TYPE(full_kernel_env_type), INTENT(inout)          :: kernel_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(qs_rho_type), POINTER                         :: rho_struct_sub
      LOGICAL, INTENT(in), OPTIONAL                      :: lsd_singlets, do_excitations
      TYPE(tddfpt_subgroup_env_type), INTENT(in), &
         OPTIONAL                                        :: sub_env
      TYPE(qs_environment_type), INTENT(in), OPTIONAL, &
         POINTER                                         :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_kernel_env'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: lsd, my_excitations, my_singlets
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ia_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_ij_r, rho_ij_r2, tau_ij_r, tau_ij_r2
      TYPE(section_vals_type), POINTER                   :: xc_fun_section

      CALL timeset(routineN, handle)

      IF (PRESENT(sub_env)) THEN
         pw_env => sub_env%pw_env
         dbcsr_dist => sub_env%dbcsr_dist
         sab_orb => sub_env%sab_orb

         nspins = SIZE(sub_env%mos_occ)
      ELSE
         CPASSERT(PRESENT(qs_env))

         CALL get_qs_env(qs_env=qs_env, &
                         pw_env=pw_env, &
                         dbcsr_dist=dbcsr_dist, &
                         sab_orb=sab_orb, &
                         dft_control=dft_control)

         nspins = dft_control%nspins

         IF (.NOT. ASSOCIATED(rho_struct_sub)) THEN
            ! Build rho_set
            NULLIFY (rho_ia_ao)
            CALL dbcsr_allocate_matrix_set(rho_ia_ao, nspins)
            DO ispin = 1, nspins
               CALL dbcsr_init_p(rho_ia_ao(ispin)%matrix)
               CALL tddfpt_dbcsr_create_by_dist(rho_ia_ao(ispin)%matrix, template=matrix_s(1)%matrix, &
                                                dbcsr_dist=dbcsr_dist, sab=sab_orb)
            END DO

            ALLOCATE (rho_struct_sub)
            CALL qs_rho_create(rho_struct_sub)
            CALL qs_rho_set(rho_struct_sub, rho_ao=rho_ia_ao)
            CALL qs_rho_rebuild(rho_struct_sub, qs_env, rebuild_ao=.FALSE., &
                                rebuild_grids=.TRUE., pw_env_external=pw_env)
         END IF
      END IF

      lsd = (nspins > 1) .OR. is_rks_triplets

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL qs_rho_get(rho_struct_sub, rho_r=rho_ij_r, tau_r=tau_ij_r)

      NULLIFY (kernel_env%xc_rho_set, kernel_env%xc_rho1_set)

      ALLOCATE (kernel_env%xc_rho_set)
      IF (is_rks_triplets) THEN
         ! we are about to compute triplet states using spin-restricted reference MOs;
         ! we still need the beta-spin density component in order to compute the TDDFT kernel
         ALLOCATE (rho_ij_r2(2))
         rho_ij_r2(1) = rho_ij_r(1)
         rho_ij_r2(2) = rho_ij_r(1)

         IF (ASSOCIATED(tau_ij_r)) THEN
            ALLOCATE (tau_ij_r2(2))
            tau_ij_r2(1) = tau_ij_r(1)
            tau_ij_r2(2) = tau_ij_r(1)
         END IF

         CALL xc_prep_2nd_deriv(kernel_env%xc_deriv_set, kernel_env%xc_rho_set, rho_ij_r2, &
                                auxbas_pw_pool, xc_section=xc_section, tau_r=tau_ij_r2)

         IF (ASSOCIATED(tau_ij_r)) DEALLOCATE (tau_ij_r2)

         DEALLOCATE (rho_ij_r2)
      ELSE
         CALL xc_prep_2nd_deriv(kernel_env%xc_deriv_set, kernel_env%xc_rho_set, rho_ij_r, &
                                auxbas_pw_pool, xc_section=xc_section, tau_r=tau_ij_r)
      END IF

      ! ++ allocate structure for response density
      kernel_env%xc_section => xc_section
      kernel_env%deriv_method_id = section_get_ival(xc_section, "XC_GRID%XC_DERIV")
      kernel_env%rho_smooth_id = section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO")

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      kernel_env%xc_rho1_cflags = xc_functionals_get_needs(functionals=xc_fun_section, lsd=lsd, &
                                                           calc_potential=.TRUE.)

      IF (.NOT. ASSOCIATED(kernel_env%xc_rho1_set)) THEN
         NULLIFY (kernel_env%xc_rho1_set)
         ALLOCATE (kernel_env%xc_rho1_set)
      END IF
      CALL xc_rho_set_create(kernel_env%xc_rho1_set, auxbas_pw_pool%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      my_excitations = .TRUE.
      IF (PRESENT(do_excitations)) my_excitations = do_excitations

      my_singlets = .FALSE.
      IF (PRESENT(lsd_singlets)) my_singlets = lsd_singlets

      kernel_env%alpha = 1.0_dp
      kernel_env%beta = 0.0_dp

      ! kernel_env%beta is taken into account in spin-restricted case only
      IF (nspins == 1 .AND. my_excitations) THEN
         IF (is_rks_triplets) THEN
            ! K_{triplets} = K_{alpha,alpha} - K_{alpha,beta}
            kernel_env%beta = -1.0_dp
         ELSE
            !                                                 alpha                 beta
            ! K_{singlets} = K_{alpha,alpha} + K_{alpha,beta} = 2 * K_{alpha,alpha} + 0 * K_{alpha,beta},
            ! due to the following relation : K_{alpha,alpha,singlets} == K_{alpha,beta,singlets}
            kernel_env%alpha = 2.0_dp

            IF (my_singlets) THEN
               kernel_env%beta = 1.0_dp
            END IF
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
