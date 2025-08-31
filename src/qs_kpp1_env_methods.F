!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              do_method_gapw,&
                                              do_method_gapw_xc
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type,&
                                              lri_kind_type
   USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gapw_densities,               ONLY: prepare_gapw_den
   USE qs_integrate_potential,          ONLY: integrate_v_rspace,&
                                              integrate_v_rspace_diagonal,&
                                              integrate_v_rspace_one_center
   USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
   USE qs_ks_atom,                      ONLY: update_ks_atom
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
   USE qs_rho_atom_types,               ONLY: rho_atom_type
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
             kpp1_did_change, &
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
!> \brief ...
!> \param rho1_xc ...
!> \param rho1 ...
!> \param xc_section ...
!> \param lrigpw ...
!> \param do_triplet ...
!> \param qs_env ...
!> \param p_env ...
!> \param calc_forces ...
!> \param calc_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE calc_kpp1(rho1_xc, rho1, xc_section, lrigpw, do_triplet, qs_env, p_env, &
                        calc_forces, calc_virial, virial)

      TYPE(qs_rho_type), POINTER                         :: rho1_xc, rho1
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: lrigpw, do_triplet
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type)                                :: p_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: calc_forces, calc_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_kpp1'

      INTEGER                                            :: handle, ikind, ispin, nkind, ns, nspins, &
                                                            output_unit
      LOGICAL                                            :: gapw, gapw_xc, lsd, my_calc_forces
      REAL(KIND=dp)                                      :: alpha, energy_hartree, energy_hartree_1c
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: k1mat, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, psmat
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type)                               :: rho1_tot_gspace
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho1_g, rho1_g_pw
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: v_hartree_rspace
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho1_r, rho1_r_pw, tau1_r, tau1_r_pw, &
                                                            v_rspace_new, v_xc, v_xc_tau
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho1_atom_set, rho_atom_set
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

      CALL kpp1_check_i_alloc(p_env%kpp1_env, qs_env, do_triplet)

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL qs_rho_get(rho1, rho_g=rho1_g)

      ! gets the tmp grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL auxbas_pw_pool%create_pw(v_hartree_rspace)

      IF (gapw .OR. gapw_xc) &
         CALL prepare_gapw_den(qs_env, p_env%local_rho_set, do_rho0=(.NOT. gapw_xc))

      ! *** calculate the hartree potential on the total density ***
      CALL auxbas_pw_pool%create_pw(rho1_tot_gspace)

      CALL pw_copy(rho1_g(1), rho1_tot_gspace)
      DO ispin = 2, nspins
         CALL pw_axpy(rho1_g(ispin), rho1_tot_gspace)
      END DO
      IF (gapw) THEN
         CALL pw_axpy(p_env%local_rho_set%rho0_mpole%rho0_s_gs, rho1_tot_gspace)
         IF (ASSOCIATED(p_env%local_rho_set%rho0_mpole%rhoz_cneo_s_gs)) THEN
            CALL pw_axpy(p_env%local_rho_set%rho0_mpole%rhoz_cneo_s_gs, rho1_tot_gspace)
         END IF
      END IF

      scf_section => section_vals_get_subs_vals(input, "DFT%SCF")
      IF (cp_print_key_should_output(logger%iter_info, scf_section, "PRINT%TOTAL_DENSITIES") &
          /= 0) THEN
         output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%TOTAL_DENSITIES", &
                                            extension=".scfLog")
         CALL print_densities(rho1, rho1_tot_gspace, output_unit)
         CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                           "PRINT%TOTAL_DENSITIES")
      END IF

      IF (.NOT. (nspins == 1 .AND. do_triplet)) THEN
         BLOCK
            TYPE(pw_c1d_gs_type) :: v_hartree_gspace
            CALL auxbas_pw_pool%create_pw(v_hartree_gspace)
            CALL pw_poisson_solve(poisson_env, rho1_tot_gspace, &
                                  energy_hartree, &
                                  v_hartree_gspace)
            CALL pw_transfer(v_hartree_gspace, v_hartree_rspace)
            CALL auxbas_pw_pool%give_back_pw(v_hartree_gspace)
         END BLOCK
         CALL pw_scale(v_hartree_rspace, v_hartree_rspace%pw_grid%dvol)
      END IF

      CALL auxbas_pw_pool%give_back_pw(rho1_tot_gspace)

      ! *** calculate the xc potential ***
      IF (gapw_xc) THEN
         CALL qs_rho_get(rho1_xc, rho_r=rho1_r, tau_r=tau1_r)
      ELSE
         CALL qs_rho_get(rho1, rho_r=rho1_r, tau_r=tau1_r)
      END IF

      IF (nspins == 1 .AND. do_triplet) THEN

         lsd = .TRUE.
         ALLOCATE (rho1_r_pw(2))
         DO ispin = 1, 2
            CALL rho1_r_pw(ispin)%create(rho1_r(1)%pw_grid)
            CALL pw_transfer(rho1_r(1), rho1_r_pw(ispin))
         END DO

         IF (ASSOCIATED(tau1_r)) THEN
            ALLOCATE (tau1_r_pw(2))
            DO ispin = 1, 2
               CALL tau1_r_pw(ispin)%create(tau1_r(1)%pw_grid)
               CALL pw_transfer(tau1_r(1), tau1_r_pw(ispin))
            END DO
         END IF

      ELSE

         rho1_r_pw => rho1_r

         tau1_r_pw => tau1_r

      END IF

      CALL xc_calc_2nd_deriv(v_xc, v_xc_tau, p_env%kpp1_env%deriv_set, p_env%kpp1_env%rho_set, &
                             rho1_r_pw, rho1_g_pw, tau1_r_pw, auxbas_pw_pool, xc_section, .FALSE., &
                             do_excitations=.TRUE., do_triplet=do_triplet, &
                             compute_virial=calc_virial, virial_xc=virial)

      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
      END DO
      v_rspace_new => v_xc
      IF (SIZE(v_xc) /= nspins) THEN
         CALL auxbas_pw_pool%give_back_pw(v_xc(2))
      END IF
      NULLIFY (v_xc)
      IF (ASSOCIATED(v_xc_tau)) THEN
      DO ispin = 1, nspins
         CALL pw_scale(v_xc_tau(ispin), v_xc_tau(ispin)%pw_grid%dvol)
      END DO
      IF (SIZE(v_xc_tau) /= nspins) THEN
         CALL auxbas_pw_pool%give_back_pw(v_xc_tau(2))
      END IF
      END IF

      IF (gapw .OR. gapw_xc) THEN
         CALL get_qs_env(qs_env, rho_atom_set=rho_atom_set)
         rho1_atom_set => p_env%local_rho_set%rho_atom_set
         CALL calculate_xc_2nd_deriv_atom(rho_atom_set, rho1_atom_set, qs_env, xc_section, para_env, &
                                          do_triplet=do_triplet)
      END IF

      IF (nspins == 1 .AND. do_triplet) THEN
         DO ispin = 1, SIZE(rho1_r_pw)
            CALL rho1_r_pw(ispin)%release()
         END DO
         DEALLOCATE (rho1_r_pw)
         IF (ASSOCIATED(tau1_r_pw)) THEN
         DO ispin = 1, SIZE(tau1_r_pw)
            CALL tau1_r_pw(ispin)%release()
         END DO
         DEALLOCATE (tau1_r_pw)
         END IF
      END IF

      alpha = 1.0_dp
      IF (nspins == 1) alpha = 2.0_dp

      !-------------------------------!
      ! Add both hartree and xc terms !
      !-------------------------------!
      DO ispin = 1, nspins
         CALL dbcsr_set(p_env%kpp1_env%v_ao(ispin)%matrix, 0.0_dp)

         ! XC and Hartree are integrated separatedly
         ! XC uses the soft basis set only
         IF (gapw_xc) THEN

            IF (nspins == 1) THEN
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
                  CALL pw_copy(v_hartree_rspace, v_rspace_new(1))

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

               CALL pw_copy(v_hartree_rspace, v_rspace_new(ispin))
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat=rho_ao(ispin), &
                                       hmat=p_env%kpp1_env%v_ao(ispin), &
                                       qs_env=qs_env, &
                                       calculate_forces=my_calc_forces, gapw=gapw)
            END IF

         ELSE

            IF (nspins == 1) THEN

               ! add hartree only for SINGLETS
               IF (.NOT. do_triplet) THEN
                  CALL pw_axpy(v_hartree_rspace, v_rspace_new(1))
               END IF
            ELSE
               CALL pw_axpy(v_hartree_rspace, v_rspace_new(ispin))
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
                  CALL para_env%sum(lri_v_int(ikind)%v_int)
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
         IF (.NOT. (nspins == 1 .AND. do_triplet)) THEN
            CALL Vh_1c_gg_integrals(qs_env, energy_hartree_1c, &
                                    p_env%hartree_local%ecoul_1c, &
                                    p_env%local_rho_set, &
                                    para_env, tddft=.TRUE., core_2nd=.TRUE.)
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
      ELSEIF (gapw_xc) THEN
         ns = SIZE(p_env%kpp1)
         ksmat(1:ns, 1:1) => p_env%kpp1(1:ns)
         ns = SIZE(rho_ao)
         psmat(1:ns, 1:1) => rho_ao(1:ns)
         CALL update_ks_atom(qs_env, ksmat, psmat, forces=my_calc_forces, tddft=.TRUE., &
                             rho_atom_external=p_env%local_rho_set%rho_atom_set)
      END IF

      CALL auxbas_pw_pool%give_back_pw(v_hartree_rspace)
      DO ispin = 1, SIZE(v_rspace_new)
         CALL auxbas_pw_pool%give_back_pw(v_rspace_new(ispin))
      END DO
      DEALLOCATE (v_rspace_new)
      IF (ASSOCIATED(v_xc_tau)) THEN
      DO ispin = 1, SIZE(v_xc_tau)
         CALL auxbas_pw_pool%give_back_pw(v_xc_tau(ispin))
      END DO
      DEALLOCATE (v_xc_tau)
      END IF

      CALL timestop(handle)
   END SUBROUTINE calc_kpp1

! **************************************************************************************************
!> \brief checks that the intenal storage is allocated, and allocs it if needed
!> \param kpp1_env the environment to check
!> \param qs_env the qs environment this kpp1_env lives in
!> \param do_triplet ...
!> \author Fawzi Mohamed
!> \note
!>      private routine
! **************************************************************************************************
   SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, do_triplet)

      TYPE(qs_kpp1_env_type)                             :: kpp1_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      LOGICAL, INTENT(IN)                                :: do_triplet

      INTEGER                                            :: ispin, nspins
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: my_rho_r, my_tau_r, rho_r, tau_r
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

         IF (nspins == 1 .AND. do_triplet) THEN
            ALLOCATE (my_rho_r(2))
            DO ispin = 1, 2
               CALL auxbas_pw_pool%create_pw(my_rho_r(ispin))
               CALL pw_axpy(rho_r(1), my_rho_r(ispin), 0.5_dp, 0.0_dp)
            END DO
            IF (dft_control%use_kinetic_energy_density) THEN
               ALLOCATE (my_tau_r(2))
               DO ispin = 1, 2
                  CALL auxbas_pw_pool%create_pw(my_tau_r(ispin))
                  CALL pw_axpy(tau_r(1), my_tau_r(ispin), 0.5_dp, 0.0_dp)
               END DO
            END IF
         ELSE
            my_rho_r => rho_r
            IF (dft_control%use_kinetic_energy_density) THEN
               my_tau_r => tau_r
            END IF
         END IF

         IF (dft_control%do_admm) THEN
            xc_section => admm_env%xc_section_primary
         ELSE
            xc_section => section_vals_get_subs_vals(input, "DFT%XC")
         END IF

         ALLOCATE (kpp1_env%deriv_set, kpp1_env%rho_set)
         CALL xc_prep_2nd_deriv(kpp1_env%deriv_set, kpp1_env%rho_set, &
                                my_rho_r, auxbas_pw_pool, &
                                xc_section=xc_section, tau_r=my_tau_r)

         IF (nspins == 1 .AND. do_triplet) THEN
            DO ispin = 1, SIZE(my_rho_r)
               CALL my_rho_r(ispin)%release()
            END DO
            DEALLOCATE (my_rho_r)
            IF (ASSOCIATED(my_tau_r)) THEN
               DO ispin = 1, SIZE(my_tau_r)
                  CALL my_tau_r(ispin)%release()
               END DO
               DEALLOCATE (my_tau_r)
            END IF
         END IF
      END IF

      ! ADMM Correction
      IF (dft_control%do_admm) THEN
         IF (admm_env%aux_exch_func /= do_admm_aux_exch_func_none) THEN
            IF (.NOT. ASSOCIATED(kpp1_env%deriv_set_admm)) THEN
               CPASSERT(.NOT. do_triplet)
               admm_xc_section => admm_env%xc_section_aux
               CALL get_admm_env(qs_env%admm_env, rho_aux_fit=rho)
               CALL qs_rho_get(rho, rho_r=rho_r)
               ALLOCATE (kpp1_env%deriv_set_admm, kpp1_env%rho_set_admm)
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
      IF (ASSOCIATED(kpp1_env%rho_set)) THEN
         CALL xc_rho_set_release(kpp1_env%rho_set)
         DEALLOCATE (kpp1_env%rho_set)
      END IF

   END SUBROUTINE kpp1_did_change

! **************************************************************************************************
!> \brief ...
!> \param rho1 ...
!> \param rho1_tot_gspace ...
!> \param out_unit ...
! **************************************************************************************************
   SUBROUTINE print_densities(rho1, rho1_tot_gspace, out_unit)

      TYPE(qs_rho_type), POINTER                         :: rho1
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: rho1_tot_gspace
      INTEGER                                            :: out_unit

      REAL(KIND=dp)                                      :: total_rho_gspace
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho1_r

      NULLIFY (tot_rho1_r)

      total_rho_gspace = pw_integrate_function(rho1_tot_gspace, isign=-1)
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
