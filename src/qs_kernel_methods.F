!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_kernel_methods
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_distribution_type,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_lval,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_p_type
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
   USE xc_rho_set_types,                ONLY: xc_rho_set_create
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kernel_methods'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   PUBLIC :: create_kernel_env

CONTAINS

! **************************************************************************************************
!> \brief Create kernel environment.
!> \param kernel_env       kernel environment (allocated and initialised on exit)
!> \param xc_section       input section which defines an exchange-correlation functional
!> \param is_rks_triplets  indicates that the triplet excited states calculation using
!>                         spin-unpolarised molecular orbitals has been requested
!> \param rho_struct_sub   ground state charge densityi, if not associated on input, it will be associated on output
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
      TYPE(full_kernel_env_type), INTENT(out)            :: kernel_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(qs_rho_type), POINTER                         :: rho_struct_sub
      LOGICAL, INTENT(in), OPTIONAL                      :: lsd_singlets, do_excitations
      TYPE(tddfpt_subgroup_env_type), INTENT(in), &
         OPTIONAL                                        :: sub_env
      TYPE(qs_environment_type), INTENT(in), OPTIONAL, &
         POINTER                                         :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_kernel_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: lsd, my_excitations, my_singlets
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ia_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ij_r, rho_ij_r2, tau_ij_r, tau_ij_r2
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
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

            CALL qs_rho_create(rho_struct_sub)
            CALL qs_rho_set(rho_struct_sub, rho_ao=rho_ia_ao)
            CALL qs_rho_rebuild(rho_struct_sub, qs_env, rebuild_ao=.FALSE., &
                                rebuild_grids=.TRUE., pw_env_external=pw_env)
         END IF
      END IF

      lsd = (nspins > 1) .OR. is_rks_triplets

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL qs_rho_get(rho_struct_sub, rho_r=rho_ij_r, tau_r=tau_ij_r)

      NULLIFY (kernel_env%xc_rho_set, kernel_env%xc_rho1_set, kernel_env%xc_deriv_set)

      IF (is_rks_triplets) THEN
         ! we are about to compute triplet states using spin-restricted reference MOs;
         ! we still need the beta-spin density component in order to compute the TDDFT kernel
         ALLOCATE (rho_ij_r2(2))
         rho_ij_r2(1)%pw => rho_ij_r(1)%pw
         rho_ij_r2(2)%pw => rho_ij_r(1)%pw

         IF (ASSOCIATED(tau_ij_r)) THEN
            ALLOCATE (tau_ij_r2(2))
            tau_ij_r2(1)%pw => tau_ij_r(1)%pw
            tau_ij_r2(2)%pw => tau_ij_r(1)%pw
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
                                                           add_basic_components=.TRUE.)

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

END MODULE qs_kernel_methods
