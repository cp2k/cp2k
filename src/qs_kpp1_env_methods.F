!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief module that builds the second order perturbation kernel
!>      kpp1 = delta_rho|_P delta_rho|_P E drho(P1) drho
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE qs_kpp1_env_methods
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_scale,&
                                              dbcsr_set
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              do_method_gapw,&
                                              do_method_gapw_xc,&
                                              tddfpt_excitations,&
                                              tddfpt_triplet
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type,&
                                              lri_kind_type
   USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix
   USE message_passing,                 ONLY: mp_sum
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_create,&
                                              pw_p_type,&
                                              pw_release,&
                                              pw_retain
   USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                              deallocate_qs_energy,&
                                              qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_integrate_potential,          ONLY: integrate_v_rspace,&
                                              integrate_v_rspace_diagonal,&
                                              integrate_v_rspace_one_center
   USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc_atom,                     ONLY: calculate_xc_2nd_deriv_atom
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_dset_release
   USE xc_rho_set_types,                ONLY: xc_rho_set_release
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpp1_env_methods'

   PUBLIC :: kpp1_create, &
             kpp1_calc_k_p_p1, &
             kpp1_calc_k_p_p1_fdiff, &
             kpp1_did_change, &
             kpp1_check_i_alloc, &
             calc_kpp1

CONTAINS

! **************************************************************************************************
!> \brief allocates and initializes a kpp1_env
!> \param kpp1_env the environment to initialize
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE kpp1_create(kpp1_env)
      TYPE(qs_kpp1_env_type)                             :: kpp1_env

      NULLIFY (kpp1_env%v_ao, kpp1_env%rho_set, kpp1_env%deriv_set, &
               kpp1_env%rho_set_admm, kpp1_env%deriv_set_admm)
   END SUBROUTINE kpp1_create

! **************************************************************************************************
!> \brief calculates the k_p_p1 kernel of the perturbation theory
!> \param p_env perturbation environment containing kpp1 kernel and kpp1_env
!> \param qs_env kpp1's qs_env
!> \param rho1 the density that represent the first direction along which
!>        you should evaluate the derivatives
!> \param rho1_xc ...
! **************************************************************************************************
   SUBROUTINE kpp1_calc_k_p_p1(p_env, qs_env, rho1, rho1_xc)

      TYPE(qs_p_env_type)                                :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho1
      TYPE(qs_rho_type), OPTIONAL, POINTER               :: rho1_xc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'kpp1_calc_k_p_p1'

      INTEGER                                            :: excitations, handle, res_etype
      LOGICAL                                            :: do_excitations, do_triplet, explicit, &
                                                            lsd_singlets
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routineN, handle)

      NULLIFY (input)

      CPASSERT(ASSOCIATED(rho1))

      CALL get_qs_env(qs_env=qs_env, &
                      input=input)

      CALL section_vals_val_get(input, "DFT%EXCITATIONS", &
                                i_val=excitations)
      CALL section_vals_val_get(input, "DFT%TDDFPT%LSD_SINGLETS", &
                                l_val=lsd_singlets)
      CALL section_vals_val_get(input, "DFT%TDDFPT%RES_ETYPE", &
                                i_val=res_etype)

      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      IF (excitations == tddfpt_excitations) THEN
         xc_section => section_vals_get_subs_vals(input, "DFT%TDDFPT%XC")
         CALL section_vals_get(xc_section, explicit=explicit)
         IF (.NOT. explicit) THEN
            xc_section => section_vals_get_subs_vals(input, "DFT%XC")
         END IF
      END IF

      do_excitations = (excitations == tddfpt_excitations)
      do_triplet = (res_etype == tddfpt_triplet)

      CALL calc_kpp1(rho1_xc, rho1, xc_section, .TRUE., &
                     lsd_singlets, .FALSE., do_excitations, &
                     do_triplet, qs_env, p_env)

      CALL timestop(handle)
   END SUBROUTINE kpp1_calc_k_p_p1

! **************************************************************************************************
!> \brief ...
!> \param rho1_xc ...
!> \param rho1 ...
!> \param xc_section ...
!> \param do_tddft ...
!> \param lsd_singlets ...
!> \param lrigpw ...
!> \param do_excitations ...
!> \param do_triplet ...
!> \param qs_env ...
!> \param p_env ...
!> \param calc_forces ...
!> \param calc_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calc_kpp1(rho1_xc, rho1, xc_section, do_tddft, lsd_singlets, lrigpw, &
                        do_excitations, do_triplet, qs_env, p_env, calc_forces, &
                        calc_virial, virial)

      TYPE(qs_rho_type), POINTER                         :: rho1_xc, rho1
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: do_tddft, lsd_singlets, lrigpw, &
                                                            do_excitations, do_triplet
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: calc_forces, calc_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_kpp1', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, ispin, nkind, ns, nspins, &
                                                            output_unit
      LOGICAL                                            :: gapw, gapw_xc, lsd, my_calc_forces
      REAL(KIND=dp)                                      :: alpha, energy_hartree, energy_hartree_1c
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: vxg
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: k1mat, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, psmat
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho1_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_g, rho1_g_pw, rho1_r, rho1_r_pw, &
                                                            tau1_r, tau1_r_pw, v_rspace_new, v_xc, &
                                                            v_xc_tau
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, scf_section

      CALL timeset(routineN, handle)

      NULLIFY (v_xc, rho1_g, pw_env, rho1_g_pw, tau1_r_pw)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(p_env%kpp1))
      CPASSERT(ASSOCIATED(p_env%kpp1_env))
      CPASSERT(ASSOCIATED(rho1))

      nspins = SIZE(p_env%kpp1)
      lsd = (nspins == 2)

      my_calc_forces = .FALSE.
      IF (PRESENT(calc_forces)) my_calc_forces = calc_forces

      CALL get_qs_env(qs_env, &
                      pw_env=pw_env, &
                      input=input, &
                      para_env=para_env, &
                      rho=rho)

      CPASSERT(ASSOCIATED(rho1))

      IF (lrigpw) THEN
         CALL get_qs_env(qs_env, &
                         lri_env=lri_env, &
                         lri_density=lri_density, &
                         atomic_kind_set=atomic_kind_set)
      END IF

      gapw = (section_get_ival(input, "DFT%QS%METHOD") == do_method_gapw)
      gapw_xc = (section_get_ival(input, "DFT%QS%METHOD") == do_method_gapw_xc)
      IF (gapw_xc) THEN
         CPASSERT(ASSOCIATED(rho1_xc))
      END IF

      CALL kpp1_check_i_alloc(p_env%kpp1_env, qs_env, do_excitations, lsd_singlets, do_triplet)

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL qs_rho_get(rho1, rho_g=rho1_g)

! gets the tmp grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      ALLOCATE (v_rspace_new(nspins))
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      IF (gapw .OR. gapw_xc) &
         CALL prepare_gapw_den(qs_env, p_env%local_rho_set, do_rho0=(.NOT. gapw_xc))

! *** calculate the hartree potential on the total density ***

      CALL pw_pool_create_pw(auxbas_pw_pool, rho1_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL pw_copy(rho1_g(1)%pw, rho1_tot_gspace%pw)
      DO ispin = 2, nspins
         CALL pw_axpy(rho1_g(ispin)%pw, rho1_tot_gspace%pw)
      END DO
      IF (gapw) &
         CALL pw_axpy(p_env%local_rho_set%rho0_mpole%rho0_s_gs%pw, rho1_tot_gspace%pw)

      scf_section => section_vals_get_subs_vals(input, "DFT%SCF")
      IF (cp_print_key_should_output(logger%iter_info, scf_section, "PRINT%TOTAL_DENSITIES") &
          /= 0) THEN
         output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%TOTAL_DENSITIES", &
                                            extension=".scfLog")
         CALL print_densities(rho1, rho1_tot_gspace, output_unit)
         CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                           "PRINT%TOTAL_DENSITIES")
      END IF

      IF (.NOT. (nspins == 1 .AND. do_excitations .AND. do_triplet)) THEN
         CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_poisson_solve(poisson_env, rho1_tot_gspace%pw, &
                               energy_hartree, &
                               v_hartree_gspace%pw)
         CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
         CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho1_tot_gspace%pw)

! *** calculate the xc potential ***
      IF (gapw_xc) THEN
         CALL qs_rho_get(rho1_xc, rho_r=rho1_r, tau_r=tau1_r)
      ELSE
         CALL qs_rho_get(rho1, rho_r=rho1_r, tau_r=tau1_r)
      END IF

      IF (nspins == 1 .AND. do_excitations .AND. &
          (lsd_singlets .OR. do_triplet)) THEN

         lsd = .TRUE.
         ALLOCATE (rho1_r_pw(2))
         DO ispin = 1, 2
            NULLIFY (rho1_r_pw(ispin)%pw)
            CALL pw_create(rho1_r_pw(ispin)%pw, rho1_r(1)%pw%pw_grid, &
                           rho1_r(1)%pw%in_use, rho1_r(1)%pw%in_space)
            CALL pw_transfer(rho1_r(1)%pw, rho1_r_pw(ispin)%pw)
         END DO

         IF (ASSOCIATED(tau1_r)) THEN
            ALLOCATE (tau1_r_pw(2))
            DO ispin = 1, 2
               NULLIFY (tau1_r_pw(ispin)%pw)
               CALL pw_create(tau1_r_pw(ispin)%pw, tau1_r(1)%pw%pw_grid, &
                              tau1_r(1)%pw%in_use, tau1_r(1)%pw%in_space)
               CALL pw_transfer(tau1_r(1)%pw, tau1_r_pw(ispin)%pw)
            END DO
         END IF

      ELSE

         ALLOCATE (rho1_r_pw(nspins))
         DO ispin = 1, nspins
            rho1_r_pw(ispin)%pw => rho1_r(ispin)%pw
            CALL pw_retain(rho1_r_pw(ispin)%pw)
         END DO

         IF (ASSOCIATED(tau1_r)) THEN
            ALLOCATE (tau1_r_pw(nspins))
            DO ispin = 1, nspins
               tau1_r_pw(ispin)%pw => tau1_r(ispin)%pw
               CALL pw_retain(tau1_r_pw(ispin)%pw)
            END DO
         END IF

      END IF

      CALL xc_calc_2nd_deriv(v_xc, v_xc_tau, p_env%kpp1_env%deriv_set, p_env%kpp1_env%rho_set, &
                             rho1_r_pw, rho1_g_pw, tau1_r_pw, auxbas_pw_pool, xc_section, .FALSE., &
                             NULL(vxg), lsd_singlets, do_excitations, do_triplet, do_tddft, &
                             compute_virial=calc_virial, virial_xc=virial)

      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         v_rspace_new(ispin)%pw => v_xc(ispin)%pw
      END DO
      IF (SIZE(v_xc) /= nspins) CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(2)%pw)
      DEALLOCATE (v_xc)
      IF (ASSOCIATED(v_xc_tau)) THEN
      DO ispin = 1, nspins
         CALL pw_scale(v_xc_tau(ispin)%pw, v_xc_tau(ispin)%pw%pw_grid%dvol)
      END DO
      IF (SIZE(v_xc_tau) /= nspins) CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc_tau(2)%pw)
      END IF

      IF (gapw) CALL calculate_xc_2nd_deriv_atom(p_env%local_rho_set, qs_env, xc_section, para_env, &
                                                 do_tddft=do_tddft, do_triplet=do_triplet)

      DO ispin = 1, SIZE(rho1_r_pw)
         CALL pw_release(rho1_r_pw(ispin)%pw)
      END DO
      DEALLOCATE (rho1_r_pw)
      IF (ASSOCIATED(tau1_r_pw)) THEN
      DO ispin = 1, SIZE(tau1_r_pw)
         IF (ASSOCIATED(tau1_r_pw(ispin)%pw)) CALL pw_release(tau1_r_pw(ispin)%pw)
      END DO
      DEALLOCATE (tau1_r_pw)
      END IF

      alpha = 1.0_dp
      IF (do_excitations .AND. nspins == 1) alpha = 2.0_dp

      !-------------------------------!
      ! Add both hartree and xc terms !
      !-------------------------------!
      DO ispin = 1, nspins
         CALL dbcsr_set(p_env%kpp1_env%v_ao(ispin)%matrix, 0.0_dp)

         ! XC and Hartree are integrated separatedly
         ! XC uses the soft basis set only
         IF (gapw_xc) THEN

            IF (do_excitations .AND. nspins == 1) THEN
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=my_calc_forces, gapw=gapw_xc)

               IF (ASSOCIATED(v_xc_tau)) THEN
                  CALL integrate_v_rspace(v_rspace=v_xc_tau(ispin), &
                                          pmat=rho_ao(ispin), &
                                          hmat=p_env%kpp1_env%v_ao(ispin), &
                                          qs_env=qs_env, &
                                          compute_tau=.TRUE., &
                                          calculate_forces=my_calc_forces, gapw=gapw_xc)
               END IF

               ! add hartree only for SINGLETS
               IF (.NOT. do_triplet) THEN
                  v_rspace_new(1)%pw%cr3d = v_hartree_rspace%pw%cr3d

                  CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                          pmat=rho_ao(ispin), &
                                          hmat=p_env%kpp1_env%v_ao(ispin), &
                                          qs_env=qs_env, &
                                          calculate_forces=my_calc_forces, gapw=gapw)
               END IF
            ELSE
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=my_calc_forces, gapw=gapw_xc)

               IF (ASSOCIATED(v_xc_tau)) THEN
                  CALL integrate_v_rspace(v_rspace=v_xc_tau(ispin), &
                                          pmat=rho_ao(ispin), &
                                          hmat=p_env%kpp1_env%v_ao(ispin), &
                                          qs_env=qs_env, &
                                          compute_tau=.TRUE., &
                                          calculate_forces=my_calc_forces, gapw=gapw_xc)
               END IF

               v_rspace_new(ispin)%pw%cr3d = v_hartree_rspace%pw%cr3d
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=my_calc_forces, gapw=gapw)
            END IF

         ELSE

            IF (do_excitations .AND. nspins == 1) THEN

               ! add hartree only for SINGLETS
               IF (.NOT. do_triplet) THEN
                  v_rspace_new(1)%pw%cr3d = v_rspace_new(1)%pw%cr3d + &
                                            v_hartree_rspace%pw%cr3d
               END IF
            ELSE
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d + &
                                             v_hartree_rspace%pw%cr3d
            END IF

            IF (lrigpw) THEN
               IF (ASSOCIATED(v_xc_tau)) CPABORT("Meta-GGA functionals not supported with LRI!")

               lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
               CALL get_qs_env(qs_env, nkind=nkind)
               DO ikind = 1, nkind
                  lri_v_int(ikind)%v_int = 0.0_dp
               END DO
               CALL integrate_v_rspace_one_center(v_rspace_new(ispin), qs_env, &
                                                  lri_v_int, .FALSE., "LRI_AUX")
               DO ikind = 1, nkind
                  CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
               END DO
               ALLOCATE (k1mat(1))
               k1mat(1)%matrix => p_env%kpp1_env%v_ao(ispin)%matrix
               IF (lri_env%exact_1c_terms) THEN
                  CALL integrate_v_rspace_diagonal(v_rspace_new(ispin), k1mat(1)%matrix, &
                                                   rho_ao(ispin)%matrix, qs_env, my_calc_forces, "ORB")
               END IF
               CALL calculate_lri_ks_matrix(lri_env, lri_v_int, k1mat, atomic_kind_set)
               DEALLOCATE (k1mat)
            ELSE
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=my_calc_forces, gapw=gapw)

               IF (ASSOCIATED(v_xc_tau)) THEN
                  CALL integrate_v_rspace(v_rspace=v_xc_tau(ispin), &
                                          pmat=rho_ao(ispin), &
                                          hmat=p_env%kpp1_env%v_ao(ispin), &
                                          qs_env=qs_env, &
                                          compute_tau=.TRUE., &
                                          calculate_forces=my_calc_forces, gapw=gapw)
               END IF
            END IF
         END IF

         CALL dbcsr_add(p_env%kpp1(ispin)%matrix, p_env%kpp1_env%v_ao(ispin)%matrix, 1.0_dp, alpha)
      END DO

      IF (gapw) THEN
         IF (.NOT. (do_excitations .AND. nspins == 1 .AND. do_triplet)) THEN
            CALL Vh_1c_gg_integrals(qs_env, energy_hartree_1c, &
                                    p_env%hartree_local%ecoul_1c, &
                                    p_env%local_rho_set, &
                                    para_env, tddft=.TRUE.)
            CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, para_env, &
                                       calculate_forces=my_calc_forces, &
                                       local_rho_set=p_env%local_rho_set)
         END IF
         !  ***  Add single atom contributions to the KS matrix ***
         ! remap pointer
         ns = SIZE(p_env%kpp1)
         ksmat(1:ns, 1:1) => p_env%kpp1(1:ns)
         ns = SIZE(rho_ao)
         psmat(1:ns, 1:1) => rho_ao(1:ns)
         CALL update_ks_atom(qs_env, ksmat, psmat, forces=my_calc_forces, tddft=.TRUE., &
                             rho_atom_external=p_env%local_rho_set%rho_atom_set)
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      DO ispin = 1, SIZE(v_rspace_new)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_new(ispin)%pw)
      END DO
      DEALLOCATE (v_rspace_new)
      IF (ASSOCIATED(v_xc_tau)) THEN
      DO ispin = 1, SIZE(v_xc_tau)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc_tau(ispin)%pw)
      END DO
      DEALLOCATE (v_xc_tau)
      END IF

      CALL timestop(handle)
   END SUBROUTINE calc_kpp1

! **************************************************************************************************
!> \brief calcualtes the k_p_p1 kernel of the perturbation theory with finite
!>      differences
!> \param qs_env kpp1's qs_env
!> \param k_p_p1 the sparse matrix that will contain the kernel k_p_p1
!> \param rho the density where to evaluate the derivatives (i.e. p along
!>        with with its grid representations, that must be valid)
!> \param rho1 the density that represent the first direction along which
!>        you should evaluate the derivatives
!> \param diff the amount of the finite difference step
!> \par History
!>      01.2003 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      useful for testing purposes.
!>      rescale my_diff depending on the norm of rho1?
! **************************************************************************************************
   SUBROUTINE kpp1_calc_k_p_p1_fdiff(qs_env, k_p_p1, rho, rho1, &
                                     diff)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: k_p_p1
      TYPE(qs_rho_type), POINTER                         :: rho, rho1
      REAL(KIND=dp), INTENT(in), OPTIONAL                :: diff

      INTEGER                                            :: ispin, nspins
      REAL(KIND=dp)                                      :: my_diff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_2, matrix_s, rho1_ao, rho_ao
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_g, rho1_r, rho_g, rho_r
      TYPE(qs_energy_type), POINTER                      :: qs_energy

      NULLIFY (ks_2, matrix_s, qs_energy, rho_ao, rho1_ao, rho_r, rho1_r, rho_g, rho1_g)
      nspins = SIZE(k_p_p1)
      my_diff = 1.0e-6_dp
      IF (PRESENT(diff)) my_diff = diff
      CALL allocate_qs_energy(qs_energy)

      CALL qs_rho_get(rho, rho_ao=rho_ao, rho_r=rho_r, rho_g=rho_g)
      CALL qs_rho_get(rho1, rho_ao=rho1_ao, rho_r=rho1_r, rho_g=rho1_g)
      CALL get_qs_env(qs_env, matrix_s=matrix_s)

! rho = rho0+h/2*rho1
      my_diff = my_diff/2.0_dp
      DO ispin = 1, SIZE(k_p_p1)
         CALL dbcsr_add(rho_ao(ispin)%matrix, rho1_ao(ispin)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=my_diff)
         rho_r(ispin)%pw%cr3d = rho_r(ispin)%pw%cr3d + my_diff*rho1_r(ispin)%pw%cr3d
         rho_g(ispin)%pw%cc = rho_g(ispin)%pw%cc + my_diff*rho1_g(ispin)%pw%cc
      END DO

      CALL qs_ks_build_kohn_sham_matrix(qs_env, &
                                        ext_ks_matrix=k_p_p1, &
                                        calculate_forces=.FALSE., &
                                        just_energy=.FALSE.)

      CALL dbcsr_allocate_matrix_set(ks_2, nspins)
      DO ispin = 1, nspins
         ALLOCATE (ks_2(ispin)%matrix)
         CALL dbcsr_copy(ks_2(ispin)%matrix, matrix_s(1)%matrix, &
                         name="tmp_ks2-"//ADJUSTL(cp_to_string(ispin)))
      END DO

! rho = rho0-h/2*rho1
      my_diff = -2.0_dp*my_diff
      DO ispin = 1, nspins
         CALL dbcsr_add(rho_ao(ispin)%matrix, rho1_ao(ispin)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=my_diff)
         rho_r(ispin)%pw%cr3d = rho_r(ispin)%pw%cr3d + my_diff*rho1_r(ispin)%pw%cr3d
         rho_g(ispin)%pw%cc = rho_g(ispin)%pw%cc + my_diff*rho1_g(ispin)%pw%cc
      END DO

      CALL qs_ks_build_kohn_sham_matrix(qs_env, &
                                        ext_ks_matrix=ks_2, &
                                        calculate_forces=.FALSE., &
                                        just_energy=.FALSE.)

! rho = rho0
      my_diff = -0.5_dp*my_diff
      DO ispin = 1, nspins
         CALL dbcsr_add(rho_ao(ispin)%matrix, rho1_ao(ispin)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=my_diff)
         rho_r(ispin)%pw%cr3d = rho_r(ispin)%pw%cr3d + my_diff*rho1_r(ispin)%pw%cr3d
         rho_g(ispin)%pw%cc = rho_g(ispin)%pw%cc + my_diff*rho1_g(ispin)%pw%cc
      END DO

! k_p_p1=(H(rho0+h/2 rho1)-H(rho0-h/2 rho1))/h
      DO ispin = 1, nspins
         CALL dbcsr_add(k_p_p1(ispin)%matrix, ks_2(ispin)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         CALL dbcsr_scale(k_p_p1(ispin)%matrix, alpha_scalar=0.5_dp/my_diff)
      END DO

      CALL dbcsr_deallocate_matrix_set(ks_2)
      CALL deallocate_qs_energy(qs_energy)
   END SUBROUTINE kpp1_calc_k_p_p1_fdiff

! **************************************************************************************************
!> \brief checks that the intenal storage is allocated, and allocs it if needed
!> \param kpp1_env the environment to check
!> \param qs_env the qs environment this kpp1_env lives in
!> \param do_excitations ...
!> \param lsd_singlets ...
!> \param do_triplet ...
!> \author Fawzi Mohamed
!> \note
!>      private routine
! **************************************************************************************************
   SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, do_excitations, lsd_singlets, do_triplet)

      TYPE(qs_kpp1_env_type)                             :: kpp1_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      LOGICAL, INTENT(IN)                                :: do_excitations, lsd_singlets, do_triplet

      CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_check_i_alloc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: my_rho_r, my_tau_r, rho_r, tau_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: admm_xc_section, input, xc_section

! ------------------------------------------------------------------

      NULLIFY (pw_env, auxbas_pw_pool, matrix_s, rho, rho_r, admm_env, dft_control, my_rho_r, my_tau_r)

      CALL get_qs_env(qs_env, pw_env=pw_env, &
                      matrix_s=matrix_s, rho=rho, input=input, &
                      admm_env=admm_env, dft_control=dft_control)

      CALL qs_rho_get(rho, rho_r=rho_r, tau_r=tau_r)
      nspins = SIZE(rho_r)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      IF (.NOT. ASSOCIATED(kpp1_env%v_ao)) THEN
         CALL dbcsr_allocate_matrix_set(kpp1_env%v_ao, nspins)
         DO ispin = 1, nspins
            ALLOCATE (kpp1_env%v_ao(ispin)%matrix)
            CALL dbcsr_copy(kpp1_env%v_ao(ispin)%matrix, matrix_s(1)%matrix, &
                            name="kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)))
         END DO
      END IF

      IF (.NOT. ASSOCIATED(kpp1_env%deriv_set)) THEN

         IF (nspins == 1 .AND. (do_excitations .AND. &
                                (lsd_singlets .OR. do_triplet))) THEN
            ALLOCATE (my_rho_r(2))
            DO ispin = 1, 2
               CALL pw_pool_create_pw(auxbas_pw_pool, my_rho_r(ispin)%pw, &
                                      use_data=rho_r(1)%pw%in_use, in_space=rho_r(1)%pw%in_space)
               my_rho_r(ispin)%pw%cr3d = 0.5_dp*rho_r(1)%pw%cr3d
            END DO
            IF (dft_control%use_kinetic_energy_density) THEN
               ALLOCATE (my_tau_r(2))
               DO ispin = 1, 2
                  CALL pw_pool_create_pw(auxbas_pw_pool, my_tau_r(ispin)%pw, &
                                         use_data=tau_r(1)%pw%in_use, in_space=tau_r(1)%pw%in_space)
                  my_tau_r(ispin)%pw%cr3d = 0.5_dp*tau_r(1)%pw%cr3d
               END DO
            END IF
         ELSE
            ALLOCATE (my_rho_r(SIZE(rho_r)))
            DO ispin = 1, SIZE(rho_r)
               my_rho_r(ispin)%pw => rho_r(ispin)%pw
               CALL pw_retain(my_rho_r(ispin)%pw)
            END DO
            IF (dft_control%use_kinetic_energy_density) THEN
               ALLOCATE (my_tau_r(SIZE(tau_r)))
               DO ispin = 1, SIZE(tau_r)
                  my_tau_r(ispin)%pw => tau_r(ispin)%pw
                  CALL pw_retain(my_tau_r(ispin)%pw)
               END DO
            END IF
         END IF

         IF (dft_control%do_admm) THEN
            xc_section => admm_env%xc_section_primary
         ELSE
            xc_section => section_vals_get_subs_vals(input, "DFT%XC")
         END IF

         ALLOCATE (kpp1_env%deriv_set)
         CALL xc_prep_2nd_deriv(kpp1_env%deriv_set, kpp1_env%rho_set, &
                                my_rho_r, auxbas_pw_pool, &
                                xc_section=xc_section, tau_r=my_tau_r)

         DO ispin = 1, SIZE(my_rho_r)
            CALL pw_release(my_rho_r(ispin)%pw)
         END DO
         DEALLOCATE (my_rho_r)
         IF (dft_control%use_kinetic_energy_density) THEN
            DO ispin = 1, SIZE(my_tau_r)
               CALL pw_release(my_tau_r(ispin)%pw)
            END DO
            DEALLOCATE (my_tau_r)
         END IF
      END IF

      ! ADMM Correction
      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func /= do_admm_aux_exch_func_none) THEN
            IF (.NOT. ASSOCIATED(kpp1_env%deriv_set_admm)) THEN
               CPASSERT(.NOT. do_triplet)
               admm_xc_section => admm_env%xc_section_aux
               CALL get_qs_env(qs_env=qs_env, rho_aux_fit=rho)
               CALL qs_rho_get(rho, rho_r=rho_r)
               ALLOCATE (kpp1_env%deriv_set_admm)
               CALL xc_prep_2nd_deriv(kpp1_env%deriv_set_admm, kpp1_env%rho_set_admm, &
                                      rho_r, auxbas_pw_pool, &
                                      xc_section=admm_xc_section)
            END IF
         END IF
      END IF

   END SUBROUTINE kpp1_check_i_alloc

! **************************************************************************************************
!> \brief function to advise of changes either in the grids
!> \param kpp1_env the kpp1_env
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE kpp1_did_change(kpp1_env)
      TYPE(qs_kpp1_env_type)                             :: kpp1_env

      IF (ASSOCIATED(kpp1_env%deriv_set)) THEN
         CALL xc_dset_release(kpp1_env%deriv_set)
         DEALLOCATE (kpp1_env%deriv_set)
         NULLIFY (kpp1_env%deriv_set)
      END IF
      CALL xc_rho_set_release(kpp1_env%rho_set)

   END SUBROUTINE kpp1_did_change

! **************************************************************************************************
!> \brief ...
!> \param rho1 ...
!> \param rho1_tot_gspace ...
!> \param out_unit ...
! **************************************************************************************************
   SUBROUTINE print_densities(rho1, rho1_tot_gspace, out_unit)

      TYPE(qs_rho_type), POINTER                         :: rho1
      TYPE(pw_p_type), INTENT(IN)                        :: rho1_tot_gspace
      INTEGER                                            :: out_unit

      REAL(KIND=dp)                                      :: total_rho_gspace
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho1_r

      NULLIFY (tot_rho1_r)

      total_rho_gspace = pw_integrate_function(rho1_tot_gspace%pw, isign=-1)
      IF (out_unit > 0) THEN
         CALL qs_rho_get(rho1, tot_rho_r=tot_rho1_r)
         WRITE (UNIT=out_unit, FMT="(T3,A,T60,F20.10)") &
            "KPP1 total charge density (r-space):", &
            accurate_sum(tot_rho1_r), &
            "KPP1 total charge density (g-space):", &
            total_rho_gspace
      END IF

   END SUBROUTINE print_densities

END MODULE qs_kpp1_env_methods
