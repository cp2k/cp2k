!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!> \author MI
! **************************************************************************************************
MODULE qs_vxc

   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE gapw_gspace_reconstruction,      ONLY: calculate_rhotot_elec_gspace
   USE input_constants,                 ONLY: sic_ad,&
                                              sic_eo,&
                                              sic_mauri_spz,&
                                              sic_mauri_us,&
                                              sic_none,&
                                              xc_none,&
                                              xc_vdw_fun_nonloc
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: get_pw_grid_info,&
                                              pw_grid_compare
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integral_ab,&
                                              pw_integrate_function,&
                                              pw_multiply,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_dispersion_nonloc,            ONLY: calculate_dispersion_nonloc
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE skala_gpw_functional,            ONLY: get_gauxc_section,&
                                              native_skala_gapw_composite_direct_ao,&
                                              native_skala_gapw_composite_reference,&
                                              skala_gpw_eval,&
                                              xc_section_uses_native_skala_evaluator
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: calc_xc_density,&
                                              xc_exc_calc,&
                                              xc_vxc_pw_create
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: qs_vxc_create, qs_xc_density

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vxc'

CONTAINS

! **************************************************************************************************
!> \brief calculates and allocates the xc potential, already reducing it to
!>      the dependence on rho and the one on tau
!> \param ks_env to get all the needed things
!> \param rho_struct density for which v_xc is calculated
!> \param xc_section ...
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic tau part of v_xc
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param exc ...
!> \param just_energy if true calculates just the energy, and does not
!>        allocate v_*_rspace
!> \param edisp ...
!> \param dispersion_env ...
!> \param adiabatic_rescale_factor ...
!> \param pw_env_external    external plane wave environment
!> \param native_skala_atom_force ...
!> \param qs_env_external parent Quickstep environment for GAPW reconstruction
!> \param native_gapw_composite_override ...
!> \param native_skala_defer_to_atom_composite prepare VXC storage while deferring the
!>        native SKALA evaluation to the atom-centered GAPW composite path
!> \par History
!>      - 05.2002 modified to use the mp_allgather function each pe
!>        computes only part of the grid and this is broadcasted to all
!>        instead of summed.
!>        This scales significantly better (e.g. factor 3 on 12 cpus
!>        32 H2O) [Joost VdV]
!>      - moved to qs_ks_methods [fawzi]
!>      - sic alterations [Joost VandeVondele]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_vxc_create(ks_env, rho_struct, xc_section, vxc_rho, vxc_tau, exc, &
                            just_energy, edisp, dispersion_env, adiabatic_rescale_factor, &
                            pw_env_external, native_skala_atom_force, qs_env_external, &
                            native_gapw_composite_override, native_skala_defer_to_atom_composite)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: vxc_rho, vxc_tau
      REAL(KIND=dp), INTENT(out)                         :: exc
      LOGICAL, INTENT(in), OPTIONAL                      :: just_energy
      REAL(KIND=dp), INTENT(out), OPTIONAL               :: edisp
      TYPE(qs_dispersion_type), OPTIONAL, POINTER        :: dispersion_env
      REAL(KIND=dp), INTENT(in), OPTIONAL                :: adiabatic_rescale_factor
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: native_skala_atom_force
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env_external
      LOGICAL, INTENT(in), OPTIONAL :: native_gapw_composite_override, &
         native_skala_defer_to_atom_composite

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_vxc_create'

      INTEGER                                            :: handle, i, ispin, j, k, mspin, myfun, &
                                                            nelec_spin(2), output_unit, vdw
      LOGICAL :: compute_virial, defer_native_skala_to_atom_composite, do_adiabatic_rescaling, &
         my_just_energy, native_gapw_composite_direct_ao, native_gapw_composite_reference, &
         native_grid_diagnostics, native_skala_grid, rho_g_valid, sic_scaling_b_zero, tau_g_valid, &
         tau_r_valid, uf_grid, vdW_nl
      REAL(KIND=dp) :: composite_hard_integral, composite_rho_max, composite_rho_min, &
         composite_rho_r_integral, composite_soft_integral, composite_tau_max, composite_tau_min, &
         composite_tau_r_integral, delta, direct_rho_integral, direct_tau_integral, exc_m, factor, &
         my_adiabatic_rescale_factor, my_scaling, nelec_s_inv, q_max, rho_composite_hard, &
         rho_composite_soft, rho_diff_l2, rho_diff_max, rho_ref_l2, tau_composite_hard, &
         tau_composite_soft, tau_diff_l2, tau_diff_max, tau_ref_l2, volume
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc_tmp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: tmp_g
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER :: rho_direct_g, rho_g, rho_hard_g, &
         rho_m_gspace, rho_smooth_g, rho_smooth_model_g, rho_struct_g, tau_direct_g, tau_hard_g, &
         tau_smooth_g, tau_smooth_model_g, tau_struct_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g, rho_nlcc_g_use, rho_nlcc_g_xc
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, diagnostics_pw_pool, &
                                                            vdw_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: my_vxc_rho, my_vxc_tau, rho_direct_r, &
         rho_hard_r, rho_m_rspace, rho_r, rho_smooth_model_r, rho_smooth_r, rho_struct_r, tau, &
         tau_direct_r, tau_hard_r, tau_smooth_model_r, tau_smooth_r, tau_struct_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc, rho_nlcc_use, rho_nlcc_xc, &
                                                            tmp_pw, weights, weights_use, &
                                                            weights_xc
      TYPE(section_vals_type), POINTER                   :: gauxc_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(vxc_rho))
      CPASSERT(.NOT. ASSOCIATED(vxc_tau))
      NULLIFY (dft_control, pw_env, auxbas_pw_pool, diagnostics_pw_pool, xc_pw_pool, vdw_pw_pool, &
               cell, my_vxc_rho, &
               tmp_pw, my_vxc_tau, rho_g, rho_r, tau, rho_m_rspace, &
               rho_m_gspace, rho_nlcc, rho_nlcc_g, rho_nlcc_g_use, rho_nlcc_g_xc, &
               rho_nlcc_use, rho_nlcc_xc, rho_struct_r, rho_struct_g, tau_struct_g, tau_struct_r, &
               weights_use, weights_xc, particle_set, rho_ao_kp, rho_direct_g, rho_direct_r, &
               rho_hard_g, rho_hard_r, rho_smooth_g, rho_smooth_model_g, rho_smooth_model_r, &
               rho_smooth_r, tau_direct_g, tau_direct_r, tau_hard_g, tau_hard_r, tau_smooth_g, &
               tau_smooth_model_g, tau_smooth_model_r, tau_smooth_r, gauxc_section)

      exc = 0.0_dp
      my_just_energy = .FALSE.
      IF (PRESENT(just_energy)) my_just_energy = just_energy
      my_adiabatic_rescale_factor = 1.0_dp
      do_adiabatic_rescaling = .FALSE.
      IF (PRESENT(adiabatic_rescale_factor)) THEN
         my_adiabatic_rescale_factor = adiabatic_rescale_factor
         do_adiabatic_rescaling = .TRUE.
      END IF

      CALL get_ks_env(ks_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env, &
                      cell=cell, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      xcint_weights=weights, &
                      virial=virial, &
                      rho_nlcc=rho_nlcc, &
                      rho_nlcc_g=rho_nlcc_g)
      rho_nlcc_use => rho_nlcc
      rho_nlcc_g_use => rho_nlcc_g
      weights_use => weights

      CALL qs_rho_get(rho_struct, &
                      tau_r_valid=tau_r_valid, &
                      tau_g_valid=tau_g_valid, &
                      rho_g_valid=rho_g_valid, &
                      rho_r=rho_struct_r, &
                      rho_g=rho_struct_g, &
                      tau_g=tau_struct_g, &
                      tau_r=tau_struct_r)
      rho_smooth_r => rho_struct_r
      rho_smooth_g => rho_struct_g
      tau_smooth_r => tau_struct_r
      tau_smooth_g => tau_struct_g

      compute_virial = virial%pv_calculate .AND. (.NOT. virial%pv_numer)
      IF (compute_virial) THEN
         virial%pv_xc = 0.0_dp
      END IF

      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)
      CALL section_vals_val_get(xc_section, "VDW_POTENTIAL%POTENTIAL_TYPE", &
                                i_val=vdw)

      vdW_nl = (vdw == xc_vdw_fun_nonloc)
      ! this combination has not been investigated
      CPASSERT(.NOT. (do_adiabatic_rescaling .AND. vdW_nl))
      ! are the necessary inputs available
      IF (.NOT. (PRESENT(dispersion_env) .AND. PRESENT(edisp))) THEN
         vdW_nl = .FALSE.
      END IF
      IF (PRESENT(edisp)) edisp = 0.0_dp
      native_skala_grid = xc_section_uses_native_skala_evaluator(xc_section)
      defer_native_skala_to_atom_composite = .FALSE.
      IF (PRESENT(native_skala_defer_to_atom_composite)) THEN
         defer_native_skala_to_atom_composite = native_skala_defer_to_atom_composite
      END IF
      CPASSERT(.NOT. defer_native_skala_to_atom_composite .OR. native_skala_grid)
      native_gapw_composite_reference = native_skala_grid .AND. &
                                        native_skala_gapw_composite_reference(xc_section) .AND. &
                                        (dft_control%qs_control%gapw .OR. &
                                         dft_control%qs_control%gapw_xc)
      IF (PRESENT(native_gapw_composite_override)) THEN
         native_gapw_composite_reference = native_skala_grid .AND. &
                                           native_gapw_composite_override .AND. &
                                           (dft_control%qs_control%gapw .OR. &
                                            dft_control%qs_control%gapw_xc)
      END IF
      native_gapw_composite_direct_ao = native_gapw_composite_reference .AND. &
                                        native_skala_gapw_composite_direct_ao(xc_section)
      native_grid_diagnostics = .FALSE.
      IF (native_gapw_composite_reference) THEN
         gauxc_section => get_gauxc_section(xc_section)
         CPASSERT(ASSOCIATED(gauxc_section))
         CALL section_vals_val_get(gauxc_section, "NATIVE_GRID_DIAGNOSTICS", &
                                   l_val=native_grid_diagnostics)
      END IF
      IF (native_gapw_composite_reference) NULLIFY (weights_use)

      IF (myfun /= xc_none .OR. vdW_nl) THEN

         ! test if the real space density is available
         CPASSERT(ASSOCIATED(rho_struct))
         IF (dft_control%nspins /= 1 .AND. dft_control%nspins /= 2) THEN
            CPABORT("nspins must be 1 or 2")
         END IF
         mspin = SIZE(rho_struct_r)
         IF (dft_control%nspins == 2 .AND. mspin == 1) THEN
            CPABORT("Spin count mismatch")
         END IF

         ! there are some options related to SIC here.
         ! Normal DFT computes E(rho_alpha,rho_beta) (or its variant E(2*rho_alpha) for non-LSD)
         ! SIC can             E(rho_alpha,rho_beta)-b*(E(rho_alpha,rho_beta)-E(rho_beta,rho_beta))
         ! or compute          E(rho_alpha,rho_beta)-b*E(rho_alpha-rho_beta,0)

         ! my_scaling is the scaling needed of the standard E(rho_alpha,rho_beta) term
         my_scaling = 1.0_dp
         SELECT CASE (dft_control%sic_method_id)
         CASE (sic_none)
            ! all fine
         CASE (sic_mauri_spz, sic_ad)
            ! no idea yet what to do here in that case
            CPASSERT(.NOT. tau_r_valid)
         CASE (sic_mauri_us)
            my_scaling = 1.0_dp - dft_control%sic_scaling_b
            ! no idea yet what to do here in that case
            CPASSERT(.NOT. tau_r_valid)
         CASE (sic_eo)
            ! NOTHING TO BE DONE
         CASE DEFAULT
            ! this case has not yet been treated here
            CPABORT("NYI")
         END SELECT

         IF (dft_control%sic_scaling_b == 0.0_dp) THEN
            sic_scaling_b_zero = .TRUE.
         ELSE
            sic_scaling_b_zero = .FALSE.
         END IF

         IF (PRESENT(pw_env_external)) THEN
            pw_env => pw_env_external
         END IF
         CALL pw_env_get(pw_env, xc_pw_pool=xc_pw_pool, auxbas_pw_pool=auxbas_pw_pool)

         IF (native_gapw_composite_reference) THEN
            CPASSERT(tau_r_valid)
            CPASSERT(PRESENT(qs_env_external))
            ALLOCATE (rho_hard_r(mspin), rho_hard_g(mspin), tau_hard_r(mspin), tau_hard_g(mspin))
            DO ispin = 1, mspin
               CALL auxbas_pw_pool%create_pw(rho_hard_r(ispin))
               CALL auxbas_pw_pool%create_pw(rho_hard_g(ispin))
               CALL auxbas_pw_pool%create_pw(tau_hard_r(ispin))
               CALL auxbas_pw_pool%create_pw(tau_hard_g(ispin))
            END DO
            IF (native_gapw_composite_direct_ao) THEN
               CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
               CPASSERT(ASSOCIATED(rho_ao_kp))
               CPASSERT(SIZE(rho_ao_kp, 2) == 1)
               DO ispin = 1, mspin
                  CALL calculate_rho_elec( &
                     matrix_p=rho_ao_kp(ispin, 1)%matrix, rho=rho_hard_r(ispin), &
                     rho_gspace=rho_hard_g(ispin), ks_env=ks_env)
                  CALL calculate_rho_elec( &
                     matrix_p=rho_ao_kp(ispin, 1)%matrix, rho=tau_hard_r(ispin), &
                     rho_gspace=tau_hard_g(ispin), ks_env=ks_env, compute_tau=.TRUE.)
               END DO
            ELSE
               CALL get_pw_grid_info(rho_hard_g(1)%pw_grid, vol=volume)
               q_max = SQRT(MAXVAL(rho_hard_g(1)%pw_grid%gsq))
               CALL calculate_rhotot_elec_gspace( &
                  qs_env=qs_env_external, auxbas_pw_pool=auxbas_pw_pool, &
                  rhotot_elec_gspace=rho_hard_g(1), q_max=q_max, &
                  rho_hard=composite_hard_integral, rho_soft=composite_soft_integral, &
                  rho_source=rho_struct, allow_nonorthorhombic=.TRUE.)
               rho_composite_hard = composite_hard_integral
               rho_composite_soft = composite_soft_integral
               CALL calculate_rhotot_elec_gspace( &
                  qs_env=qs_env_external, auxbas_pw_pool=auxbas_pw_pool, &
                  rhotot_elec_gspace=tau_hard_g(1), q_max=q_max, &
                  rho_hard=composite_hard_integral, rho_soft=composite_soft_integral, &
                  compute_tau=.TRUE., rho_source=rho_struct, allow_nonorthorhombic=.TRUE.)
               tau_composite_hard = composite_hard_integral
               tau_composite_soft = composite_soft_integral
               IF (para_env%mepos == 0) THEN
                  output_unit = cp_logger_get_default_io_unit()
                  IF (output_unit > 0) THEN
                     WRITE (UNIT=output_unit, FMT="(/,T2,A,2(1X,ES19.11))") &
                        "SKALA_GPW| Composite rho hard and soft integrals", &
                        rho_composite_hard, rho_composite_soft
                     WRITE (UNIT=output_unit, FMT="(T2,A,2(1X,ES19.11))") &
                        "SKALA_GPW| Composite tau hard and soft integrals", &
                        tau_composite_hard, tau_composite_soft
                  END IF
               END IF
               IF (mspin == 1) THEN
                  CALL pw_scale(rho_hard_g(1), 1.0_dp/volume)
                  CALL pw_scale(tau_hard_g(1), 1.0_dp/volume)
               ELSE
                  CALL calculate_rhotot_elec_gspace( &
                     qs_env=qs_env_external, auxbas_pw_pool=auxbas_pw_pool, &
                     rhotot_elec_gspace=rho_hard_g(2), q_max=q_max, &
                     rho_hard=composite_hard_integral, rho_soft=composite_soft_integral, &
                     fsign=-1.0_dp, rho_source=rho_struct, allow_nonorthorhombic=.TRUE.)
                  CALL calculate_rhotot_elec_gspace( &
                     qs_env=qs_env_external, auxbas_pw_pool=auxbas_pw_pool, &
                     rhotot_elec_gspace=tau_hard_g(2), q_max=q_max, &
                     rho_hard=composite_hard_integral, rho_soft=composite_soft_integral, &
                     fsign=-1.0_dp, compute_tau=.TRUE., rho_source=rho_struct, &
                     allow_nonorthorhombic=.TRUE.)
                  CALL pw_scale(rho_hard_g(1), 0.5_dp/volume)
                  CALL pw_scale(rho_hard_g(2), 0.5_dp/volume)
                  CALL auxbas_pw_pool%create_pw(tmp_g)
                  CALL pw_copy(rho_hard_g(1), tmp_g)
                  CALL pw_axpy(rho_hard_g(2), rho_hard_g(1), 1.0_dp)
                  CALL pw_axpy(rho_hard_g(2), tmp_g, -1.0_dp)
                  CALL pw_copy(tmp_g, rho_hard_g(2))
                  CALL pw_scale(tau_hard_g(1), 0.5_dp/volume)
                  CALL pw_scale(tau_hard_g(2), 0.5_dp/volume)
                  CALL pw_copy(tau_hard_g(1), tmp_g)
                  CALL pw_axpy(tau_hard_g(2), tau_hard_g(1), 1.0_dp)
                  CALL pw_axpy(tau_hard_g(2), tmp_g, -1.0_dp)
                  CALL pw_copy(tmp_g, tau_hard_g(2))
                  CALL auxbas_pw_pool%give_back_pw(tmp_g)
               END IF
            END IF
            DO ispin = 1, mspin
               CALL pw_transfer(rho_hard_g(ispin), rho_hard_r(ispin))
               CALL pw_transfer(tau_hard_g(ispin), tau_hard_r(ispin))
            END DO
            IF (native_grid_diagnostics .AND. .NOT. native_gapw_composite_direct_ao) THEN
               CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
               CPASSERT(ASSOCIATED(rho_ao_kp))
               CPASSERT(SIZE(rho_ao_kp, 2) == 1)
               ALLOCATE (rho_direct_r(mspin), rho_direct_g(mspin), &
                         tau_direct_r(mspin), tau_direct_g(mspin))
               rho_diff_l2 = 0.0_dp
               rho_diff_max = 0.0_dp
               rho_ref_l2 = 0.0_dp
               tau_diff_l2 = 0.0_dp
               tau_diff_max = 0.0_dp
               tau_ref_l2 = 0.0_dp
               direct_rho_integral = 0.0_dp
               direct_tau_integral = 0.0_dp
               DO ispin = 1, mspin
                  CALL auxbas_pw_pool%create_pw(rho_direct_r(ispin))
                  CALL auxbas_pw_pool%create_pw(rho_direct_g(ispin))
                  CALL auxbas_pw_pool%create_pw(tau_direct_r(ispin))
                  CALL auxbas_pw_pool%create_pw(tau_direct_g(ispin))
                  CALL calculate_rho_elec( &
                     matrix_p=rho_ao_kp(ispin, 1)%matrix, rho=rho_direct_r(ispin), &
                     rho_gspace=rho_direct_g(ispin), ks_env=ks_env)
                  CALL calculate_rho_elec( &
                     matrix_p=rho_ao_kp(ispin, 1)%matrix, rho=tau_direct_r(ispin), &
                     rho_gspace=tau_direct_g(ispin), ks_env=ks_env, compute_tau=.TRUE.)
                  direct_rho_integral = direct_rho_integral + &
                                        pw_integrate_function(rho_direct_r(ispin))
                  direct_tau_integral = direct_tau_integral + &
                                        pw_integrate_function(tau_direct_r(ispin))
                  DO k = LBOUND(rho_direct_r(ispin)%array, 3), UBOUND(rho_direct_r(ispin)%array, 3)
                     DO j = LBOUND(rho_direct_r(ispin)%array, 2), UBOUND(rho_direct_r(ispin)%array, 2)
                        DO i = LBOUND(rho_direct_r(ispin)%array, 1), UBOUND(rho_direct_r(ispin)%array, 1)
                           delta = rho_hard_r(ispin)%array(i, j, k) - &
                                   rho_direct_r(ispin)%array(i, j, k)
                           rho_diff_l2 = rho_diff_l2 + delta*delta
                           rho_diff_max = MAX(rho_diff_max, ABS(delta))
                           rho_ref_l2 = rho_ref_l2 + rho_direct_r(ispin)%array(i, j, k)**2
                           delta = tau_hard_r(ispin)%array(i, j, k) - &
                                   tau_direct_r(ispin)%array(i, j, k)
                           tau_diff_l2 = tau_diff_l2 + delta*delta
                           tau_diff_max = MAX(tau_diff_max, ABS(delta))
                           tau_ref_l2 = tau_ref_l2 + tau_direct_r(ispin)%array(i, j, k)**2
                        END DO
                     END DO
                  END DO
               END DO
               CALL para_env%sum(rho_diff_l2)
               CALL para_env%sum(rho_ref_l2)
               CALL para_env%sum(tau_diff_l2)
               CALL para_env%sum(tau_ref_l2)
               CALL para_env%max(rho_diff_max)
               CALL para_env%max(tau_diff_max)
               rho_diff_l2 = SQRT(rho_diff_l2/MAX(rho_ref_l2, TINY(1.0_dp)))
               tau_diff_l2 = SQRT(tau_diff_l2/MAX(tau_ref_l2, TINY(1.0_dp)))
               IF (para_env%mepos == 0) THEN
                  output_unit = cp_logger_get_default_io_unit()
                  IF (output_unit > 0) THEN
                     WRITE (UNIT=output_unit, FMT="(T2,A,2(1X,ES19.11))") &
                        "SKALA_GPW| Direct AO rho and tau integrals", &
                        direct_rho_integral, direct_tau_integral
                     WRITE (UNIT=output_unit, FMT="(T2,A,2(1X,ES19.11))") &
                        "SKALA_GPW| Composite/direct relative L2 rho and tau", &
                        rho_diff_l2, tau_diff_l2
                     WRITE (UNIT=output_unit, FMT="(T2,A,2(1X,ES19.11))") &
                        "SKALA_GPW| Composite/direct Linf rho and tau", &
                        rho_diff_max, tau_diff_max
                  END IF
               END IF
               DO ispin = 1, mspin
                  CALL auxbas_pw_pool%give_back_pw(rho_direct_r(ispin))
                  CALL auxbas_pw_pool%give_back_pw(rho_direct_g(ispin))
                  CALL auxbas_pw_pool%give_back_pw(tau_direct_r(ispin))
                  CALL auxbas_pw_pool%give_back_pw(tau_direct_g(ispin))
               END DO
               DEALLOCATE (rho_direct_r, rho_direct_g, tau_direct_r, tau_direct_g)
            END IF
            composite_rho_r_integral = 0.0_dp
            composite_tau_r_integral = 0.0_dp
            composite_rho_min = HUGE(1.0_dp)
            composite_rho_max = -HUGE(1.0_dp)
            composite_tau_min = HUGE(1.0_dp)
            composite_tau_max = -HUGE(1.0_dp)
            DO ispin = 1, mspin
               composite_rho_r_integral = composite_rho_r_integral + &
                                          pw_integrate_function(rho_hard_r(ispin))
               composite_tau_r_integral = composite_tau_r_integral + &
                                          pw_integrate_function(tau_hard_r(ispin))
               composite_rho_min = MIN(composite_rho_min, MINVAL(rho_hard_r(ispin)%array))
               composite_rho_max = MAX(composite_rho_max, MAXVAL(rho_hard_r(ispin)%array))
               composite_tau_min = MIN(composite_tau_min, MINVAL(tau_hard_r(ispin)%array))
               composite_tau_max = MAX(composite_tau_max, MAXVAL(tau_hard_r(ispin)%array))
            END DO
            CALL para_env%min(composite_rho_min)
            CALL para_env%max(composite_rho_max)
            CALL para_env%min(composite_tau_min)
            CALL para_env%max(composite_tau_max)
            IF (para_env%mepos == 0) THEN
               output_unit = cp_logger_get_default_io_unit()
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T2,A,1X,ES19.11)") &
                     "SKALA_GPW| Composite real-grid rho integral", composite_rho_r_integral
                  WRITE (UNIT=output_unit, FMT="(T2,A,1X,ES19.11)") &
                     "SKALA_GPW| Composite real-grid tau integral", composite_tau_r_integral
                  WRITE (UNIT=output_unit, FMT="(T2,A,2(1X,ES19.11))") &
                     "SKALA_GPW| Composite real-grid rho min and max", &
                     composite_rho_min, composite_rho_max
                  WRITE (UNIT=output_unit, FMT="(T2,A,2(1X,ES19.11))") &
                     "SKALA_GPW| Composite real-grid tau min and max", &
                     composite_tau_min, composite_tau_max
               END IF
            END IF
            rho_struct_r => rho_hard_r
            rho_struct_g => rho_hard_g
            tau_struct_r => tau_hard_r
            tau_struct_g => tau_hard_g
            rho_g_valid = .TRUE.
            tau_g_valid = .TRUE.
         END IF

         uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)

         IF (.NOT. uf_grid) THEN
            rho_r => rho_struct_r

            IF (tau_r_valid) THEN
               tau => tau_struct_r
            END IF

            ! for gradient corrected functional the density in g space might
            ! be useful so if we have it, we pass it in
            IF (rho_g_valid) THEN
               rho_g => rho_struct_g
            END IF
         ELSE
            CPASSERT(rho_g_valid)
            ALLOCATE (rho_r(mspin))
            ALLOCATE (rho_g(mspin))
            DO ispin = 1, mspin
               CALL xc_pw_pool%create_pw(rho_g(ispin))
               CALL pw_transfer(rho_struct_g(ispin), rho_g(ispin))
            END DO
            DO ispin = 1, mspin
               CALL xc_pw_pool%create_pw(rho_r(ispin))
               CALL pw_transfer(rho_g(ispin), rho_r(ispin))
            END DO
            IF (tau_r_valid) THEN
               ALLOCATE (tau(mspin))
               DO ispin = 1, mspin
                  CALL xc_pw_pool%create_pw(tau(ispin))
                  BLOCK
                     TYPE(pw_c1d_gs_type) :: tau_g_aux, tau_g_xc
                     CALL xc_pw_pool%create_pw(tau_g_xc)
                     IF (tau_g_valid) THEN
                        CALL pw_transfer(tau_struct_g(ispin), tau_g_xc)
                     ELSE
                        CALL auxbas_pw_pool%create_pw(tau_g_aux)
                        CALL pw_transfer(tau_struct_r(ispin), tau_g_aux)
                        CALL pw_transfer(tau_g_aux, tau_g_xc)
                        CALL auxbas_pw_pool%give_back_pw(tau_g_aux)
                     END IF
                     CALL pw_transfer(tau_g_xc, tau(ispin))
                     CALL xc_pw_pool%give_back_pw(tau_g_xc)
                  END BLOCK
               END DO
            END IF
            IF (ASSOCIATED(weights) .AND. .NOT. native_gapw_composite_reference) THEN
               ALLOCATE (weights_xc)
               CALL xc_pw_pool%create_pw(weights_xc)
               BLOCK
                  TYPE(pw_c1d_gs_type) :: weights_g_aux, weights_g_xc
                  CALL auxbas_pw_pool%create_pw(weights_g_aux)
                  CALL xc_pw_pool%create_pw(weights_g_xc)
                  CALL pw_transfer(weights, weights_g_aux)
                  CALL pw_transfer(weights_g_aux, weights_g_xc)
                  CALL pw_transfer(weights_g_xc, weights_xc)
                  CALL xc_pw_pool%give_back_pw(weights_g_xc)
                  CALL auxbas_pw_pool%give_back_pw(weights_g_aux)
               END BLOCK
               weights_use => weights_xc
            END IF
            IF (ASSOCIATED(rho_nlcc)) THEN
               CPASSERT(ASSOCIATED(rho_nlcc_g))
               ALLOCATE (rho_nlcc_g_xc, rho_nlcc_xc)
               CALL xc_pw_pool%create_pw(rho_nlcc_g_xc)
               CALL xc_pw_pool%create_pw(rho_nlcc_xc)
               CALL pw_transfer(rho_nlcc_g, rho_nlcc_g_xc)
               CALL pw_transfer(rho_nlcc_g_xc, rho_nlcc_xc)
               rho_nlcc_use => rho_nlcc_xc
               rho_nlcc_g_use => rho_nlcc_g_xc
            END IF
         END IF

         IF (native_gapw_composite_reference) THEN
            composite_rho_r_integral = 0.0_dp
            composite_tau_r_integral = 0.0_dp
            DO ispin = 1, mspin
               composite_rho_r_integral = composite_rho_r_integral + &
                                          pw_integrate_function(rho_r(ispin))
               composite_tau_r_integral = composite_tau_r_integral + &
                                          pw_integrate_function(tau(ispin))
            END DO
            IF (para_env%mepos == 0) THEN
               output_unit = cp_logger_get_default_io_unit()
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T2,A,1X,ES19.11)") &
                     "SKALA_GPW| Composite XC-grid rho integral", composite_rho_r_integral
                  WRITE (UNIT=output_unit, FMT="(T2,A,1X,ES19.11)") &
                     "SKALA_GPW| Composite XC-grid tau integral", composite_tau_r_integral
               END IF
            END IF
         END IF

         ! NLCC augments rho and grad(rho); tau remains valence-only.
         IF (ASSOCIATED(rho_nlcc_use)) THEN
            factor = 1.0_dp
            DO ispin = 1, mspin
               CALL pw_axpy(rho_nlcc_use, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g_use, rho_g(ispin), factor)
            END DO
         END IF

         !
         ! here the rho_r, rho_g, tau is what it should be
         ! we get back the right my_vxc_rho and my_vxc_tau as required
         !
         IF (defer_native_skala_to_atom_composite) THEN
            ! The atom-centered composite path supplies the complete SKALA energy and
            ! VXC below qs_vxc_create. Keep only storage for that VXC (and any additive
            ! nonlocal-dispersion potential) instead of evaluating the smooth grid here.
            ALLOCATE (my_vxc_rho(mspin), my_vxc_tau(mspin))
            DO ispin = 1, mspin
               CALL xc_pw_pool%create_pw(my_vxc_rho(ispin))
               CALL xc_pw_pool%create_pw(my_vxc_tau(ispin))
               CALL pw_zero(my_vxc_rho(ispin))
               CALL pw_zero(my_vxc_tau(ispin))
            END DO
         ELSE IF (native_skala_grid) THEN
            CALL skala_gpw_eval(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, exc=exc, &
                                rho_r=rho_r, rho_g=rho_g, tau=tau, xc_section=xc_section, &
                                weights=weights_use, pw_pool=xc_pw_pool, &
                                particle_set=particle_set, cell=cell, &
                                compute_virial=compute_virial, virial_xc=virial%pv_xc, &
                                just_energy=my_just_energy, atom_force=native_skala_atom_force)
         ELSE IF (my_just_energy) THEN
            exc = xc_exc_calc(rho_r=rho_r, tau=tau, &
                              rho_g=rho_g, xc_section=xc_section, &
                              weights=weights_use, pw_pool=xc_pw_pool)

         ELSE
            CALL xc_vxc_pw_create(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, rho_r=rho_r, &
                                  rho_g=rho_g, tau=tau, exc=exc, &
                                  xc_section=xc_section, &
                                  weights=weights_use, pw_pool=xc_pw_pool, &
                                  compute_virial=compute_virial, &
                                  virial_xc=virial%pv_xc)
         END IF

         ! The composite density is assembled through an R-to-G-to-R projection. Apply the
         ! discrete adjoint of that projection before integrating VXC with the smooth AO density.
         IF (native_gapw_composite_reference .AND. .NOT. my_just_energy) THEN
            IF (uf_grid) THEN
               diagnostics_pw_pool => xc_pw_pool
            ELSE
               diagnostics_pw_pool => auxbas_pw_pool
            END IF
            CALL diagnostics_pw_pool%create_pw(tmp_g)
            DO ispin = 1, mspin
               CALL pw_transfer(my_vxc_rho(ispin), tmp_g)
               CALL pw_transfer(tmp_g, my_vxc_rho(ispin))
               CALL pw_transfer(my_vxc_tau(ispin), tmp_g)
               CALL pw_transfer(tmp_g, my_vxc_tau(ispin))
            END DO
            CALL diagnostics_pw_pool%give_back_pw(tmp_g)
            NULLIFY (diagnostics_pw_pool)
         END IF

         IF (native_gapw_composite_reference .AND. native_grid_diagnostics .AND. &
             .NOT. my_just_energy) THEN
            IF (uf_grid) THEN
               diagnostics_pw_pool => xc_pw_pool
            ELSE
               diagnostics_pw_pool => auxbas_pw_pool
            END IF
            ALLOCATE (rho_smooth_model_g(mspin), rho_smooth_model_r(mspin), &
                      tau_smooth_model_g(mspin), tau_smooth_model_r(mspin))
            DO ispin = 1, mspin
               CALL diagnostics_pw_pool%create_pw(rho_smooth_model_g(ispin))
               CALL diagnostics_pw_pool%create_pw(rho_smooth_model_r(ispin))
               CALL diagnostics_pw_pool%create_pw(tau_smooth_model_g(ispin))
               CALL diagnostics_pw_pool%create_pw(tau_smooth_model_r(ispin))
               IF (uf_grid) THEN
                  IF (rho_g_valid) THEN
                     CALL pw_transfer(rho_smooth_g(ispin), rho_smooth_model_g(ispin))
                  ELSE
                     CALL auxbas_pw_pool%create_pw(tmp_g)
                     CALL pw_transfer(rho_smooth_r(ispin), tmp_g)
                     CALL pw_transfer(tmp_g, rho_smooth_model_g(ispin))
                     CALL auxbas_pw_pool%give_back_pw(tmp_g)
                  END IF
                  IF (tau_g_valid) THEN
                     CALL pw_transfer(tau_smooth_g(ispin), tau_smooth_model_g(ispin))
                  ELSE
                     CALL auxbas_pw_pool%create_pw(tmp_g)
                     CALL pw_transfer(tau_smooth_r(ispin), tmp_g)
                     CALL pw_transfer(tmp_g, tau_smooth_model_g(ispin))
                     CALL auxbas_pw_pool%give_back_pw(tmp_g)
                  END IF
               ELSE
                  IF (rho_g_valid) THEN
                     CALL pw_copy(rho_smooth_g(ispin), rho_smooth_model_g(ispin))
                  ELSE
                     CALL pw_transfer(rho_smooth_r(ispin), rho_smooth_model_g(ispin))
                  END IF
                  IF (tau_g_valid) THEN
                     CALL pw_copy(tau_smooth_g(ispin), tau_smooth_model_g(ispin))
                  ELSE
                     CALL pw_transfer(tau_smooth_r(ispin), tau_smooth_model_g(ispin))
                  END IF
               END IF
               CALL pw_transfer(rho_smooth_model_g(ispin), rho_smooth_model_r(ispin))
               CALL pw_transfer(tau_smooth_model_g(ispin), tau_smooth_model_r(ispin))
               IF (ASSOCIATED(rho_nlcc_use)) THEN
                  CALL pw_axpy(rho_nlcc_use, rho_smooth_model_r(ispin), 1.0_dp)
                  CALL pw_axpy(rho_nlcc_g_use, rho_smooth_model_g(ispin), 1.0_dp)
               END IF
            END DO
            CALL diagnose_gapw_composite_direction( &
               rho_r, rho_g, tau, rho_smooth_model_r, rho_smooth_model_g, &
               tau_smooth_model_r, my_vxc_rho, my_vxc_tau, xc_section, weights_use, &
               diagnostics_pw_pool, particle_set, cell, para_env)
            DO ispin = 1, mspin
               CALL diagnostics_pw_pool%give_back_pw(rho_smooth_model_g(ispin))
               CALL diagnostics_pw_pool%give_back_pw(rho_smooth_model_r(ispin))
               CALL diagnostics_pw_pool%give_back_pw(tau_smooth_model_g(ispin))
               CALL diagnostics_pw_pool%give_back_pw(tau_smooth_model_r(ispin))
            END DO
            DEALLOCATE (rho_smooth_model_g, rho_smooth_model_r, &
                        tau_smooth_model_g, tau_smooth_model_r)
            NULLIFY (diagnostics_pw_pool)
         END IF

         ! remove the nlcc densities (keep stuff in original state)
         IF (ASSOCIATED(rho_nlcc_use)) THEN
            factor = -1.0_dp
            DO ispin = 1, mspin
               CALL pw_axpy(rho_nlcc_use, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g_use, rho_g(ispin), factor)
            END DO
         END IF

         ! calclulate non-local vdW functional
         ! only if this XC_SECTION has it
         ! if yes, we use the dispersion_env from ks_env
         ! this is dangerous, as it assumes a special connection xc_section -> qs_env
         IF (vdW_nl) THEN
            CALL get_ks_env(ks_env=ks_env, para_env=para_env)
            ! no SIC functionals allowed
            CPASSERT(dft_control%sic_method_id == sic_none)
            !
            CALL pw_env_get(pw_env, vdw_pw_pool=vdw_pw_pool)
            IF (my_just_energy) THEN
               CALL calculate_dispersion_nonloc(my_vxc_rho, rho_r, rho_g, edisp, dispersion_env, &
                                                my_just_energy, vdw_pw_pool, xc_pw_pool, para_env)
            ELSE
               CALL calculate_dispersion_nonloc(my_vxc_rho, rho_r, rho_g, edisp, dispersion_env, &
                                                my_just_energy, vdw_pw_pool, xc_pw_pool, para_env, virial=virial)
            END IF
         END IF

         !! Apply rescaling to the potential if requested
         IF (.NOT. my_just_energy) THEN
            IF (do_adiabatic_rescaling) THEN
               IF (ASSOCIATED(my_vxc_rho)) THEN
                  DO ispin = 1, SIZE(my_vxc_rho)
                     CALL pw_scale(my_vxc_rho(ispin), my_adiabatic_rescale_factor)
                  END DO
               END IF
            END IF
         END IF

         IF (my_scaling /= 1.0_dp) THEN
            exc = exc*my_scaling
            IF (ASSOCIATED(my_vxc_rho)) THEN
               DO ispin = 1, SIZE(my_vxc_rho)
                  CALL pw_scale(my_vxc_rho(ispin), my_scaling)
               END DO
            END IF
            IF (ASSOCIATED(my_vxc_tau)) THEN
               DO ispin = 1, SIZE(my_vxc_tau)
                  CALL pw_scale(my_vxc_tau(ispin), my_scaling)
               END DO
            END IF
         END IF

         ! we have pw data for the xc, qs_ks requests coeff structure, here we transfer
         ! pw -> coeff
         IF (ASSOCIATED(my_vxc_rho)) THEN
            vxc_rho => my_vxc_rho
            NULLIFY (my_vxc_rho)
         END IF
         IF (ASSOCIATED(my_vxc_tau)) THEN
            vxc_tau => my_vxc_tau
            NULLIFY (my_vxc_tau)
         END IF
         IF (uf_grid) THEN
            DO ispin = 1, SIZE(rho_r)
               CALL xc_pw_pool%give_back_pw(rho_r(ispin))
            END DO
            DEALLOCATE (rho_r)
            IF (ASSOCIATED(rho_g)) THEN
               DO ispin = 1, SIZE(rho_g)
                  CALL xc_pw_pool%give_back_pw(rho_g(ispin))
               END DO
               DEALLOCATE (rho_g)
            END IF
         END IF

         ! compute again the xc but now for Exc(m,o) and the opposite sign
         IF (dft_control%sic_method_id == sic_mauri_spz .AND. .NOT. sic_scaling_b_zero) THEN
            ALLOCATE (rho_m_rspace(2), rho_m_gspace(2))
            CALL xc_pw_pool%create_pw(rho_m_gspace(1))
            CALL xc_pw_pool%create_pw(rho_m_rspace(1))
            CALL pw_copy(rho_struct_r(1), rho_m_rspace(1))
            CALL pw_axpy(rho_struct_r(2), rho_m_rspace(1), alpha=-1._dp)
            CALL pw_copy(rho_struct_g(1), rho_m_gspace(1))
            CALL pw_axpy(rho_struct_g(2), rho_m_gspace(1), alpha=-1._dp)
            ! bit sad, these will be just zero...
            CALL xc_pw_pool%create_pw(rho_m_gspace(2))
            CALL xc_pw_pool%create_pw(rho_m_rspace(2))
            CALL pw_zero(rho_m_rspace(2))
            CALL pw_zero(rho_m_gspace(2))

            IF (my_just_energy) THEN
               exc_m = xc_exc_calc(rho_r=rho_m_rspace, tau=tau, &
                                   rho_g=rho_m_gspace, xc_section=xc_section, &
                                   weights=weights_use, pw_pool=xc_pw_pool)
            ELSE
               ! virial untested
               CPASSERT(.NOT. compute_virial)
               CALL xc_vxc_pw_create(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, rho_r=rho_m_rspace, &
                                     rho_g=rho_m_gspace, tau=tau, exc=exc_m, &
                                     xc_section=xc_section, &
                                     weights=weights_use, pw_pool=xc_pw_pool, &
                                     compute_virial=.FALSE., &
                                     virial_xc=virial_xc_tmp)
            END IF

            exc = exc - dft_control%sic_scaling_b*exc_m

            ! and take care of the potential only vxc_rho is taken into account
            IF (.NOT. my_just_energy) THEN
               CALL pw_axpy(my_vxc_rho(1), vxc_rho(1), -dft_control%sic_scaling_b)
               CALL pw_axpy(my_vxc_rho(1), vxc_rho(2), dft_control%sic_scaling_b)
               CALL my_vxc_rho(1)%release()
               CALL my_vxc_rho(2)%release()
               DEALLOCATE (my_vxc_rho)
            END IF

            DO ispin = 1, 2
               CALL xc_pw_pool%give_back_pw(rho_m_rspace(ispin))
               CALL xc_pw_pool%give_back_pw(rho_m_gspace(ispin))
            END DO
            DEALLOCATE (rho_m_rspace)
            DEALLOCATE (rho_m_gspace)

         END IF

         ! now we have - sum_s N_s * Exc(rho_s/N_s,0)
         IF (dft_control%sic_method_id == sic_ad .AND. .NOT. sic_scaling_b_zero) THEN

            ! find out how many elecs we have
            CALL get_ks_env(ks_env, nelectron_spin=nelec_spin)

            ALLOCATE (rho_m_rspace(2), rho_m_gspace(2))
            DO ispin = 1, 2
               CALL xc_pw_pool%create_pw(rho_m_gspace(ispin))
               CALL xc_pw_pool%create_pw(rho_m_rspace(ispin))
            END DO

            DO ispin = 1, 2
               IF (nelec_spin(ispin) > 0.0_dp) THEN
                  nelec_s_inv = 1.0_dp/nelec_spin(ispin)
               ELSE
                  ! does it matter if there are no electrons with this spin (H) ?
                  nelec_s_inv = 0.0_dp
               END IF
               CALL pw_copy(rho_struct_r(ispin), rho_m_rspace(1))
               CALL pw_copy(rho_struct_g(ispin), rho_m_gspace(1))
               CALL pw_scale(rho_m_rspace(1), nelec_s_inv)
               CALL pw_scale(rho_m_gspace(1), nelec_s_inv)
               CALL pw_zero(rho_m_rspace(2))
               CALL pw_zero(rho_m_gspace(2))

               IF (my_just_energy) THEN
                  exc_m = xc_exc_calc(rho_r=rho_m_rspace, tau=tau, &
                                      rho_g=rho_m_gspace, xc_section=xc_section, &
                                      weights=weights_use, pw_pool=xc_pw_pool)
               ELSE
                  ! virial untested
                  CPASSERT(.NOT. compute_virial)
                  CALL xc_vxc_pw_create(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, rho_r=rho_m_rspace, &
                                        rho_g=rho_m_gspace, tau=tau, exc=exc_m, &
                                        xc_section=xc_section, &
                                        weights=weights_use, pw_pool=xc_pw_pool, &
                                        compute_virial=.FALSE., &
                                        virial_xc=virial_xc_tmp)
               END IF

               exc = exc - dft_control%sic_scaling_b*nelec_spin(ispin)*exc_m

               ! and take care of the potential only vxc_rho is taken into account
               IF (.NOT. my_just_energy) THEN
                  CALL pw_axpy(my_vxc_rho(1), vxc_rho(ispin), -dft_control%sic_scaling_b)
                  CALL my_vxc_rho(1)%release()
                  CALL my_vxc_rho(2)%release()
                  DEALLOCATE (my_vxc_rho)
               END IF
            END DO

            DO ispin = 1, 2
               CALL xc_pw_pool%give_back_pw(rho_m_rspace(ispin))
               CALL xc_pw_pool%give_back_pw(rho_m_gspace(ispin))
            END DO
            DEALLOCATE (rho_m_rspace)
            DEALLOCATE (rho_m_gspace)

         END IF

         ! compute again the xc but now for Exc(n_down,n_down)
         IF (dft_control%sic_method_id == sic_mauri_us .AND. .NOT. sic_scaling_b_zero) THEN
            ALLOCATE (rho_r(2))
            rho_r(1) = rho_struct_r(2)
            rho_r(2) = rho_struct_r(2)
            IF (rho_g_valid) THEN
               ALLOCATE (rho_g(2))
               rho_g(1) = rho_struct_g(2)
               rho_g(2) = rho_struct_g(2)
            END IF

            IF (my_just_energy) THEN
               exc_m = xc_exc_calc(rho_r=rho_r, tau=tau, &
                                   rho_g=rho_g, xc_section=xc_section, &
                                   weights=weights_use, pw_pool=xc_pw_pool)
            ELSE
               ! virial untested
               CPASSERT(.NOT. compute_virial)
               CALL xc_vxc_pw_create(vxc_rho=my_vxc_rho, vxc_tau=my_vxc_tau, rho_r=rho_r, &
                                     rho_g=rho_g, tau=tau, exc=exc_m, &
                                     xc_section=xc_section, &
                                     weights=weights_use, pw_pool=xc_pw_pool, &
                                     compute_virial=.FALSE., &
                                     virial_xc=virial_xc_tmp)
            END IF

            exc = exc + dft_control%sic_scaling_b*exc_m

            ! and take care of the potential
            IF (.NOT. my_just_energy) THEN
               ! both go to minority spin
               CALL pw_axpy(my_vxc_rho(1), vxc_rho(2), 2.0_dp*dft_control%sic_scaling_b)
               CALL my_vxc_rho(1)%release()
               CALL my_vxc_rho(2)%release()
               DEALLOCATE (my_vxc_rho)
            END IF
            DEALLOCATE (rho_r, rho_g)

         END IF

         !
         ! cleanups
         !
         IF (uf_grid .AND. (ASSOCIATED(vxc_rho) .OR. ASSOCIATED(vxc_tau))) THEN
            BLOCK
               TYPE(pw_r3d_rs_type) :: tmp_pw
               TYPE(pw_c1d_gs_type) :: tmp_g, tmp_g2
               CALL xc_pw_pool%create_pw(tmp_g)
               CALL auxbas_pw_pool%create_pw(tmp_g2)
               IF (ASSOCIATED(vxc_rho)) THEN
                  DO ispin = 1, SIZE(vxc_rho)
                     CALL auxbas_pw_pool%create_pw(tmp_pw)
                     CALL pw_transfer(vxc_rho(ispin), tmp_g)
                     CALL pw_transfer(tmp_g, tmp_g2)
                     CALL pw_transfer(tmp_g2, tmp_pw)
                     CALL xc_pw_pool%give_back_pw(vxc_rho(ispin))
                     vxc_rho(ispin) = tmp_pw
                  END DO
               END IF
               IF (ASSOCIATED(vxc_tau)) THEN
                  DO ispin = 1, SIZE(vxc_tau)
                     CALL auxbas_pw_pool%create_pw(tmp_pw)
                     CALL pw_transfer(vxc_tau(ispin), tmp_g)
                     CALL pw_transfer(tmp_g, tmp_g2)
                     CALL pw_transfer(tmp_g2, tmp_pw)
                     CALL xc_pw_pool%give_back_pw(vxc_tau(ispin))
                     vxc_tau(ispin) = tmp_pw
                  END DO
               END IF
               CALL auxbas_pw_pool%give_back_pw(tmp_g2)
               CALL xc_pw_pool%give_back_pw(tmp_g)
            END BLOCK
         END IF
         IF (ASSOCIATED(tau) .AND. uf_grid) THEN
            DO ispin = 1, SIZE(tau)
               CALL xc_pw_pool%give_back_pw(tau(ispin))
            END DO
            DEALLOCATE (tau)
         END IF
         IF (ASSOCIATED(weights_xc)) THEN
            CALL xc_pw_pool%give_back_pw(weights_xc)
            DEALLOCATE (weights_xc)
         END IF
         IF (ASSOCIATED(rho_nlcc_xc)) THEN
            CALL xc_pw_pool%give_back_pw(rho_nlcc_xc)
            DEALLOCATE (rho_nlcc_xc)
         END IF
         IF (ASSOCIATED(rho_nlcc_g_xc)) THEN
            CALL xc_pw_pool%give_back_pw(rho_nlcc_g_xc)
            DEALLOCATE (rho_nlcc_g_xc)
         END IF

         IF (ASSOCIATED(rho_hard_r)) THEN
            DO ispin = 1, SIZE(rho_hard_r)
               CALL auxbas_pw_pool%give_back_pw(rho_hard_r(ispin))
               CALL auxbas_pw_pool%give_back_pw(rho_hard_g(ispin))
               CALL auxbas_pw_pool%give_back_pw(tau_hard_r(ispin))
               CALL auxbas_pw_pool%give_back_pw(tau_hard_g(ispin))
            END DO
            DEALLOCATE (rho_hard_r, rho_hard_g, tau_hard_r, tau_hard_g)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_vxc_create

! **************************************************************************************************
!> \brief Check the SKALA energy derivative along the reconstructed GAPW one-center correction.
!> \param rho_combined_r combined smooth plus hard-minus-soft density
!> \param rho_combined_g reciprocal representation of rho_combined_r
!> \param tau_combined combined kinetic-energy density
!> \param rho_smooth_r smooth density on the model grid
!> \param rho_smooth_g reciprocal representation of rho_smooth_r
!> \param tau_smooth smooth kinetic-energy density on the model grid
!> \param vxc_rho derivative of the SKALA energy with respect to rho
!> \param vxc_tau derivative of the SKALA energy with respect to tau
!> \param xc_section ...
!> \param weights ...
!> \param pw_pool ...
!> \param particle_set ...
!> \param cell ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE diagnose_gapw_composite_direction( &
      rho_combined_r, rho_combined_g, tau_combined, rho_smooth_r, rho_smooth_g, tau_smooth, &
      vxc_rho, vxc_tau, xc_section, weights, pw_pool, particle_set, cell, para_env)

      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_combined_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_combined_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau_combined, rho_smooth_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_smooth_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau_smooth, vxc_rho, vxc_tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_r3d_rs_type), POINTER                      :: weights
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env

      REAL(KIND=dp), DIMENSION(3), PARAMETER :: fd_steps = [1.0E-2_dp, 3.0E-3_dp, 1.0E-3_dp]

      INTEGER                                            :: ispin, istep, nspin, output_unit
      REAL(KIND=dp) :: exc_minus, exc_plus, fd_derivative, fd_step, rho_contraction, &
         rho_g_contraction, tau_contraction, tau_g_contraction, total_contraction
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_dummy
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_test_g, tau_test_g, vxc_rho_g, &
                                                            vxc_tau_g
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_test_r, tau_test, vxc_dummy_rho, &
                                                            vxc_dummy_tau

      NULLIFY (rho_test_g, rho_test_r, tau_test, tau_test_g, vxc_dummy_rho, vxc_dummy_tau, &
               vxc_rho_g, vxc_tau_g)
      nspin = SIZE(rho_combined_r)
      CPASSERT(SIZE(rho_combined_g) == nspin)
      CPASSERT(SIZE(tau_combined) == nspin)
      CPASSERT(SIZE(rho_smooth_r) == nspin)
      CPASSERT(SIZE(rho_smooth_g) == nspin)
      CPASSERT(SIZE(tau_smooth) == nspin)
      CPASSERT(SIZE(vxc_rho) == nspin)
      CPASSERT(SIZE(vxc_tau) == nspin)

      ALLOCATE (rho_test_g(nspin), rho_test_r(nspin), tau_test(nspin), tau_test_g(nspin), &
                vxc_rho_g(nspin), vxc_tau_g(nspin))
      DO ispin = 1, nspin
         CALL pw_pool%create_pw(rho_test_g(ispin))
         CALL pw_pool%create_pw(rho_test_r(ispin))
         CALL pw_pool%create_pw(tau_test(ispin))
         CALL pw_pool%create_pw(tau_test_g(ispin))
         CALL pw_pool%create_pw(vxc_rho_g(ispin))
         CALL pw_pool%create_pw(vxc_tau_g(ispin))
         CALL pw_copy(rho_combined_g(ispin), rho_test_g(ispin))
         CALL pw_axpy(rho_smooth_g(ispin), rho_test_g(ispin), -1.0_dp)
         CALL pw_copy(rho_combined_r(ispin), rho_test_r(ispin))
         CALL pw_axpy(rho_smooth_r(ispin), rho_test_r(ispin), -1.0_dp)
         CALL pw_copy(tau_combined(ispin), tau_test(ispin))
         CALL pw_axpy(tau_smooth(ispin), tau_test(ispin), -1.0_dp)
         CALL pw_transfer(tau_test(ispin), tau_test_g(ispin))
         CALL pw_transfer(vxc_rho(ispin), vxc_rho_g(ispin))
         CALL pw_transfer(vxc_tau(ispin), vxc_tau_g(ispin))
      END DO
      rho_contraction = 0.0_dp
      rho_g_contraction = 0.0_dp
      tau_contraction = 0.0_dp
      tau_g_contraction = 0.0_dp
      DO ispin = 1, nspin
         rho_contraction = rho_contraction + pw_integral_ab(rho_test_r(ispin), vxc_rho(ispin))
         rho_g_contraction = rho_g_contraction + pw_integral_ab(rho_test_g(ispin), vxc_rho_g(ispin))
         tau_contraction = tau_contraction + pw_integral_ab(tau_test(ispin), vxc_tau(ispin))
         tau_g_contraction = tau_g_contraction + pw_integral_ab(tau_test_g(ispin), vxc_tau_g(ispin))
      END DO

      total_contraction = rho_contraction + tau_contraction
      DO istep = 1, SIZE(fd_steps)
         fd_step = fd_steps(istep)
         DO ispin = 1, nspin
            CALL pw_copy(rho_combined_g(ispin), rho_test_g(ispin))
            CALL pw_axpy(rho_combined_g(ispin), rho_test_g(ispin), fd_step)
            CALL pw_axpy(rho_smooth_g(ispin), rho_test_g(ispin), -fd_step)
            CALL pw_copy(rho_combined_r(ispin), rho_test_r(ispin))
            CALL pw_axpy(rho_combined_r(ispin), rho_test_r(ispin), fd_step)
            CALL pw_axpy(rho_smooth_r(ispin), rho_test_r(ispin), -fd_step)
            CALL pw_copy(tau_combined(ispin), tau_test(ispin))
            CALL pw_axpy(tau_combined(ispin), tau_test(ispin), fd_step)
            CALL pw_axpy(tau_smooth(ispin), tau_test(ispin), -fd_step)
         END DO
         CALL skala_gpw_eval( &
            vxc_rho=vxc_dummy_rho, vxc_tau=vxc_dummy_tau, exc=exc_plus, rho_r=rho_test_r, &
            rho_g=rho_test_g, tau=tau_test, xc_section=xc_section, weights=weights, pw_pool=pw_pool, &
            particle_set=particle_set, cell=cell, compute_virial=.FALSE., virial_xc=virial_dummy, &
            just_energy=.TRUE.)
         CPASSERT(.NOT. ASSOCIATED(vxc_dummy_rho))
         CPASSERT(.NOT. ASSOCIATED(vxc_dummy_tau))

         DO ispin = 1, nspin
            CALL pw_copy(rho_combined_g(ispin), rho_test_g(ispin))
            CALL pw_axpy(rho_combined_g(ispin), rho_test_g(ispin), -fd_step)
            CALL pw_axpy(rho_smooth_g(ispin), rho_test_g(ispin), fd_step)
            CALL pw_copy(rho_combined_r(ispin), rho_test_r(ispin))
            CALL pw_axpy(rho_combined_r(ispin), rho_test_r(ispin), -fd_step)
            CALL pw_axpy(rho_smooth_r(ispin), rho_test_r(ispin), fd_step)
            CALL pw_copy(tau_combined(ispin), tau_test(ispin))
            CALL pw_axpy(tau_combined(ispin), tau_test(ispin), -fd_step)
            CALL pw_axpy(tau_smooth(ispin), tau_test(ispin), fd_step)
         END DO
         CALL skala_gpw_eval( &
            vxc_rho=vxc_dummy_rho, vxc_tau=vxc_dummy_tau, exc=exc_minus, rho_r=rho_test_r, &
            rho_g=rho_test_g, tau=tau_test, xc_section=xc_section, weights=weights, pw_pool=pw_pool, &
            particle_set=particle_set, cell=cell, compute_virial=.FALSE., virial_xc=virial_dummy, &
            just_energy=.TRUE.)
         CPASSERT(.NOT. ASSOCIATED(vxc_dummy_rho))
         CPASSERT(.NOT. ASSOCIATED(vxc_dummy_tau))

         fd_derivative = (exc_plus - exc_minus)/(2.0_dp*fd_step)
         IF (para_env%mepos == 0) THEN
            output_unit = cp_logger_get_default_io_unit()
            IF (output_unit > 0) THEN
               WRITE (output_unit, "(/,T2,A,1X,ES12.4,2(1X,ES20.10))") &
                  "SKALA_GPW| Composite-direction FD step and energies", &
                  fd_step, exc_plus, exc_minus
               WRITE (output_unit, "(T2,A,1X,ES20.10)") &
                  "SKALA_GPW| Composite-direction energy FD", fd_derivative
               WRITE (output_unit, "(T2,A,1X,ES20.10)") &
                  "SKALA_GPW| Composite-direction FD/VXC difference", &
                  fd_derivative - total_contraction
            END IF
         END IF
      END DO
      IF (para_env%mepos == 0) THEN
         output_unit = cp_logger_get_default_io_unit()
         IF (output_unit > 0) THEN
            WRITE (output_unit, "(T2,A,1X,ES20.10)") &
               "SKALA_GPW| Composite-direction rho VXC contraction", rho_contraction
            WRITE (output_unit, "(T2,A,1X,ES20.10)") &
               "SKALA_GPW| Composite-direction rho G-space contraction", rho_g_contraction
            WRITE (output_unit, "(T2,A,1X,ES20.10)") &
               "SKALA_GPW| Composite-direction tau VXC contraction", tau_contraction
            WRITE (output_unit, "(T2,A,1X,ES20.10)") &
               "SKALA_GPW| Composite-direction tau G-space contraction", tau_g_contraction
            WRITE (output_unit, "(T2,A,1X,ES20.10)") &
               "SKALA_GPW| Composite-direction total VXC contraction", total_contraction
         END IF
      END IF

      DO ispin = 1, nspin
         CALL pw_pool%give_back_pw(rho_test_g(ispin))
         CALL pw_pool%give_back_pw(rho_test_r(ispin))
         CALL pw_pool%give_back_pw(tau_test(ispin))
         CALL pw_pool%give_back_pw(tau_test_g(ispin))
         CALL pw_pool%give_back_pw(vxc_rho_g(ispin))
         CALL pw_pool%give_back_pw(vxc_tau_g(ispin))
      END DO
      DEALLOCATE (rho_test_g, rho_test_r, tau_test, tau_test_g, vxc_rho_g, vxc_tau_g)

   END SUBROUTINE diagnose_gapw_composite_direction

! **************************************************************************************************
!> \brief calculates the XC density: E_xc(r) - V_xc(r)*rho(r)  or  E_xc(r)/rho(r)
!> \param ks_env to get all the needed things
!> \param rho_struct density
!> \param xc_section ...
!> \param dispersion_env ...
!> \param xc_ener will contain the xc energy density E_xc(r) - V_xc(r)*rho(r)
!> \param xc_den will contain the xc energy density E_xc(r)/rho(r)
!> \param exc will contain the xc energy density E_xc(r)
!> \param vxc ...
!> \param vtau ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_xc_density(ks_env, rho_struct, xc_section, dispersion_env, &
                            xc_ener, xc_den, exc, vxc, vtau)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(qs_dispersion_type), OPTIONAL, POINTER        :: dispersion_env
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: xc_ener, xc_den
      TYPE(pw_r3d_rs_type), OPTIONAL                     :: exc
      TYPE(pw_r3d_rs_type), DIMENSION(:), OPTIONAL       :: vxc, vtau

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_xc_density'

      INTEGER                                            :: handle, ispin, mspin, myfun, nspins, vdw
      LOGICAL                                            :: rho_g_valid, tau_g_valid, tau_r_valid, &
                                                            uf_grid, vdW_nl
      REAL(KIND=dp)                                      :: edisp, excint, factor, rho_cutoff
      REAL(KIND=dp), DIMENSION(3, 3)                     :: vdum
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g, rho_struct_g, tau_g, tau_struct_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g, rho_nlcc_g_use, rho_nlcc_g_xc
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, vdw_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type)                               :: exc_r
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, rho_struct_r, tau_r, &
                                                            tau_struct_r, vxc_rho, vxc_tau
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc, rho_nlcc_use, rho_nlcc_xc, &
                                                            weights, weights_use, weights_xc

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, pw_env, auxbas_pw_pool, xc_pw_pool, vdw_pw_pool, cell, &
               rho_g, rho_struct_g, tau_g, tau_struct_g, rho_nlcc, rho_nlcc_g, &
               rho_nlcc_g_use, rho_nlcc_g_xc, rho_nlcc_use, rho_nlcc_xc, rho_r, &
               rho_struct_r, tau_r, tau_struct_r, vxc_rho, vxc_tau, weights, &
               weights_use, weights_xc)

      CALL get_ks_env(ks_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env, &
                      cell=cell, &
                      xcint_weights=weights, &
                      rho_nlcc=rho_nlcc, &
                      rho_nlcc_g=rho_nlcc_g)

      CALL qs_rho_get(rho_struct, &
                      tau_r_valid=tau_r_valid, &
                      tau_g_valid=tau_g_valid, &
                      rho_g_valid=rho_g_valid, &
                      rho_r=rho_struct_r, &
                      rho_g=rho_struct_g, &
                      tau_r=tau_struct_r, &
                      tau_g=tau_struct_g)
      nspins = dft_control%nspins
      mspin = SIZE(rho_struct_r)
      rho_r => rho_struct_r
      rho_g => rho_struct_g
      tau_r => tau_struct_r
      tau_g => tau_struct_g
      rho_nlcc_use => rho_nlcc
      rho_nlcc_g_use => rho_nlcc_g
      weights_use => weights

      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=myfun)
      CALL section_vals_val_get(xc_section, "VDW_POTENTIAL%POTENTIAL_TYPE", i_val=vdw)
      vdW_nl = (vdw == xc_vdw_fun_nonloc)
      IF (PRESENT(xc_ener)) THEN
         IF (tau_r_valid) THEN
            CALL cp_warn(__LOCATION__, "Tau contribution will not be correctly handled")
         END IF
      END IF
      IF (vdW_nl) THEN
         CALL cp_warn(__LOCATION__, "vdW functional contribution will be ignored")
      END IF

      CALL pw_env_get(pw_env, xc_pw_pool=xc_pw_pool, auxbas_pw_pool=auxbas_pw_pool)
      uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)

      IF (PRESENT(xc_ener)) THEN
         CALL pw_zero(xc_ener)
      END IF
      IF (PRESENT(xc_den)) THEN
         CALL pw_zero(xc_den)
      END IF
      IF (PRESENT(exc)) THEN
         CALL pw_zero(exc)
      END IF
      IF (PRESENT(vxc)) THEN
         DO ispin = 1, nspins
            CALL pw_zero(vxc(ispin))
         END DO
      END IF
      IF (PRESENT(vtau)) THEN
         DO ispin = 1, nspins
            CALL pw_zero(vtau(ispin))
         END DO
      END IF

      IF (myfun /= xc_none) THEN

         CPASSERT(ASSOCIATED(rho_struct))
         CPASSERT(dft_control%sic_method_id == sic_none)

         IF (uf_grid) THEN
            NULLIFY (rho_r, rho_g, tau_r, tau_g)
            IF (rho_g_valid) THEN
               CALL create_density_on_pool(xc_pw_pool, rho_struct_g, rho_r, rho_g)
            ELSE IF (ASSOCIATED(rho_struct_r)) THEN
               CALL create_density_on_pool_from_r(auxbas_pw_pool, xc_pw_pool, rho_struct_r, rho_r, rho_g)
            ELSE
               CPABORT("Fine Grid in qs_xc_density requires rho_r or rho_g")
            END IF
            IF (tau_r_valid) THEN
               IF (tau_g_valid) THEN
                  CALL create_density_on_pool(xc_pw_pool, tau_struct_g, tau_r, tau_g)
               ELSE IF (ASSOCIATED(tau_struct_r)) THEN
                  CALL create_density_on_pool_from_r(auxbas_pw_pool, xc_pw_pool, tau_struct_r, tau_r, tau_g)
               ELSE
                  CPABORT("Fine Grid in qs_xc_density requires tau_r or tau_g")
               END IF
            END IF
            IF (ASSOCIATED(weights)) THEN
               ALLOCATE (weights_xc)
               CALL xc_pw_pool%create_pw(weights_xc)
               CALL transfer_rspace_between_pools(auxbas_pw_pool, xc_pw_pool, weights, weights_xc)
               weights_use => weights_xc
            END IF
            IF (ASSOCIATED(rho_nlcc)) THEN
               CPASSERT(ASSOCIATED(rho_nlcc_g))
               ALLOCATE (rho_nlcc_g_xc, rho_nlcc_xc)
               CALL xc_pw_pool%create_pw(rho_nlcc_g_xc)
               CALL xc_pw_pool%create_pw(rho_nlcc_xc)
               CALL pw_transfer(rho_nlcc_g, rho_nlcc_g_xc)
               CALL pw_transfer(rho_nlcc_g_xc, rho_nlcc_xc)
               rho_nlcc_use => rho_nlcc_xc
               rho_nlcc_g_use => rho_nlcc_g_xc
            END IF
         END IF

         ! add the nlcc densities
         IF (ASSOCIATED(rho_nlcc_use)) THEN
            factor = 1.0_dp
            DO ispin = 1, mspin
               CALL pw_axpy(rho_nlcc_use, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g_use, rho_g(ispin), factor)
            END DO
         END IF
         NULLIFY (vxc_rho, vxc_tau)
         CALL xc_vxc_pw_create(vxc_rho=vxc_rho, vxc_tau=vxc_tau, rho_r=rho_r, &
                               rho_g=rho_g, tau=tau_r, exc=excint, &
                               xc_section=xc_section, &
                               weights=weights_use, pw_pool=xc_pw_pool, &
                               compute_virial=.FALSE., &
                               virial_xc=vdum, &
                               exc_r=exc_r)
         ! calclulate non-local vdW functional
         ! only if this XC_SECTION has it
         ! if yes, we use the dispersion_env from ks_env
         ! this is dangerous, as it assumes a special connection xc_section -> qs_env
         IF (vdW_nl) THEN
            CALL get_ks_env(ks_env=ks_env, para_env=para_env)
            ! no SIC functionals allowed
            CPASSERT(dft_control%sic_method_id == sic_none)
            !
            CALL pw_env_get(pw_env, vdw_pw_pool=vdw_pw_pool)
            CALL calculate_dispersion_nonloc(vxc_rho, rho_r, rho_g, edisp, dispersion_env, &
                                             .FALSE., vdw_pw_pool, xc_pw_pool, para_env)
         END IF

         ! remove the nlcc densities (keep stuff in original state)
         IF (ASSOCIATED(rho_nlcc_use)) THEN
            factor = -1.0_dp
            DO ispin = 1, mspin
               CALL pw_axpy(rho_nlcc_use, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g_use, rho_g(ispin), factor)
            END DO
         END IF
         !
         IF (PRESENT(xc_den)) THEN
            rho_cutoff = 1.E-14_dp
            IF (uf_grid) THEN
               BLOCK
                  TYPE(pw_r3d_rs_type) :: tmp_pw
                  CALL xc_pw_pool%create_pw(tmp_pw)
                  CALL pw_copy(exc_r, tmp_pw)
                  CALL calc_xc_density(tmp_pw, rho_r, rho_cutoff)
                  CALL transfer_rspace_between_pools(xc_pw_pool, auxbas_pw_pool, tmp_pw, xc_den)
                  CALL xc_pw_pool%give_back_pw(tmp_pw)
               END BLOCK
            ELSE
               CALL pw_copy(exc_r, xc_den)
               CALL calc_xc_density(xc_den, rho_r, rho_cutoff)
            END IF
         END IF
         IF (PRESENT(xc_ener)) THEN
            IF (uf_grid) THEN
               BLOCK
                  TYPE(pw_r3d_rs_type) :: tmp_pw
                  CALL xc_pw_pool%create_pw(tmp_pw)
                  CALL pw_copy(exc_r, tmp_pw)
                  DO ispin = 1, nspins
                     CALL pw_multiply(tmp_pw, vxc_rho(ispin), rho_r(ispin), alpha=-1.0_dp)
                  END DO
                  CALL transfer_rspace_between_pools(xc_pw_pool, auxbas_pw_pool, tmp_pw, xc_ener)
                  CALL xc_pw_pool%give_back_pw(tmp_pw)
               END BLOCK
            ELSE
               CALL pw_copy(exc_r, xc_ener)
               DO ispin = 1, nspins
                  CALL pw_multiply(xc_ener, vxc_rho(ispin), rho_r(ispin), alpha=-1.0_dp)
               END DO
            END IF
         END IF
         IF (PRESENT(exc)) THEN
            IF (uf_grid) THEN
               CALL transfer_rspace_between_pools(xc_pw_pool, auxbas_pw_pool, exc_r, exc)
            ELSE
               CALL pw_copy(exc_r, exc)
            END IF
         END IF
         IF (PRESENT(vxc)) THEN
            DO ispin = 1, nspins
               IF (uf_grid) THEN
                  CALL transfer_rspace_between_pools(xc_pw_pool, auxbas_pw_pool, vxc_rho(ispin), vxc(ispin))
               ELSE
                  CALL pw_copy(vxc_rho(ispin), vxc(ispin))
               END IF
            END DO
         END IF
         IF (PRESENT(vtau) .AND. ASSOCIATED(vxc_tau)) THEN
            DO ispin = 1, nspins
               IF (uf_grid) THEN
                  CALL transfer_rspace_between_pools(xc_pw_pool, auxbas_pw_pool, vxc_tau(ispin), vtau(ispin))
               ELSE
                  CALL pw_copy(vxc_tau(ispin), vtau(ispin))
               END IF
            END DO
         END IF
         ! remove arrays
         IF (ASSOCIATED(vxc_rho)) THEN
            DO ispin = 1, nspins
               CALL vxc_rho(ispin)%release()
            END DO
            DEALLOCATE (vxc_rho)
         END IF
         IF (ASSOCIATED(vxc_tau)) THEN
            DO ispin = 1, nspins
               CALL vxc_tau(ispin)%release()
            END DO
            DEALLOCATE (vxc_tau)
         END IF
         CALL exc_r%release()
         IF (uf_grid) THEN
            CALL give_back_density_on_pool(xc_pw_pool, rho_r, rho_g)
            IF (ASSOCIATED(tau_r)) CALL give_back_density_on_pool(xc_pw_pool, tau_r, tau_g)
            IF (ASSOCIATED(weights_xc)) THEN
               CALL xc_pw_pool%give_back_pw(weights_xc)
               DEALLOCATE (weights_xc)
            END IF
            IF (ASSOCIATED(rho_nlcc_xc)) THEN
               CALL xc_pw_pool%give_back_pw(rho_nlcc_xc)
               DEALLOCATE (rho_nlcc_xc)
            END IF
            IF (ASSOCIATED(rho_nlcc_g_xc)) THEN
               CALL xc_pw_pool%give_back_pw(rho_nlcc_g_xc)
               DEALLOCATE (rho_nlcc_g_xc)
            END IF
         END IF
         !
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_xc_density

! **************************************************************************************************
!> \brief transfers an r-space PW between two pools and writes into an existing target PW
!> \param source_pw_pool ...
!> \param target_pw_pool ...
!> \param source ...
!> \param TARGET ...
! **************************************************************************************************
   SUBROUTINE transfer_rspace_between_pools(source_pw_pool, target_pw_pool, source, TARGET)
      TYPE(pw_pool_type), POINTER                        :: source_pw_pool, target_pw_pool
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: source, TARGET

      TYPE(pw_c1d_gs_type)                               :: source_g, target_g

      CPASSERT(ASSOCIATED(source_pw_pool))
      CPASSERT(ASSOCIATED(target_pw_pool))

      IF (pw_grid_compare(source_pw_pool%pw_grid, target_pw_pool%pw_grid)) THEN
         CALL pw_copy(source, TARGET)
      ELSE
         CALL source_pw_pool%create_pw(source_g)
         CALL target_pw_pool%create_pw(target_g)
         CALL pw_transfer(source, source_g)
         CALL pw_transfer(source_g, target_g)
         CALL pw_transfer(target_g, TARGET)
         CALL target_pw_pool%give_back_pw(target_g)
         CALL source_pw_pool%give_back_pw(source_g)
      END IF

   END SUBROUTINE transfer_rspace_between_pools

! **************************************************************************************************
!> \brief transfers a g-space density to a given PW pool and creates its r-space representation
!> \param pw_pool ...
!> \param rho_g_in ...
!> \param rho_r_out ...
!> \param rho_g_out ...
! **************************************************************************************************
   SUBROUTINE create_density_on_pool(pw_pool, rho_g_in, rho_r_out, rho_g_out)
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g_in
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r_out
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g_out

      INTEGER                                            :: ispin, nspins

      CPASSERT(ASSOCIATED(pw_pool))
      CPASSERT(ASSOCIATED(rho_g_in))

      nspins = SIZE(rho_g_in)
      ALLOCATE (rho_r_out(nspins), rho_g_out(nspins))
      DO ispin = 1, nspins
         CALL pw_pool%create_pw(rho_g_out(ispin))
         CALL pw_pool%create_pw(rho_r_out(ispin))
         CALL pw_transfer(rho_g_in(ispin), rho_g_out(ispin))
         CALL pw_transfer(rho_g_out(ispin), rho_r_out(ispin))
      END DO

   END SUBROUTINE create_density_on_pool

! **************************************************************************************************
!> \brief transfers an r-space density to a given PW pool and creates its g-space representation
!> \param source_pw_pool ...
!> \param target_pw_pool ...
!> \param rho_r_in ...
!> \param rho_r_out ...
!> \param rho_g_out ...
! **************************************************************************************************
   SUBROUTINE create_density_on_pool_from_r(source_pw_pool, target_pw_pool, rho_r_in, rho_r_out, rho_g_out)
      TYPE(pw_pool_type), POINTER                        :: source_pw_pool, target_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r_in, rho_r_out
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g_out

      INTEGER                                            :: ispin, nspins
      TYPE(pw_c1d_gs_type)                               :: rho_g_in

      CPASSERT(ASSOCIATED(source_pw_pool))
      CPASSERT(ASSOCIATED(target_pw_pool))
      CPASSERT(ASSOCIATED(rho_r_in))

      nspins = SIZE(rho_r_in)
      ALLOCATE (rho_r_out(nspins), rho_g_out(nspins))
      DO ispin = 1, nspins
         CALL source_pw_pool%create_pw(rho_g_in)
         CALL target_pw_pool%create_pw(rho_g_out(ispin))
         CALL target_pw_pool%create_pw(rho_r_out(ispin))
         CALL pw_transfer(rho_r_in(ispin), rho_g_in)
         CALL pw_transfer(rho_g_in, rho_g_out(ispin))
         CALL pw_transfer(rho_g_out(ispin), rho_r_out(ispin))
         CALL source_pw_pool%give_back_pw(rho_g_in)
      END DO

   END SUBROUTINE create_density_on_pool_from_r

! **************************************************************************************************
!> \brief returns temporary density arrays to the given PW pool
!> \param pw_pool ...
!> \param rho_r ...
!> \param rho_g ...
! **************************************************************************************************
   SUBROUTINE give_back_density_on_pool(pw_pool, rho_r, rho_g)
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g

      INTEGER                                            :: ispin

      CPASSERT(ASSOCIATED(pw_pool))

      IF (ASSOCIATED(rho_r)) THEN
         DO ispin = 1, SIZE(rho_r)
            CALL pw_pool%give_back_pw(rho_r(ispin))
         END DO
         DEALLOCATE (rho_r)
      END IF
      IF (ASSOCIATED(rho_g)) THEN
         DO ispin = 1, SIZE(rho_g)
            CALL pw_pool%give_back_pw(rho_g(ispin))
         END DO
         DEALLOCATE (rho_g)
      END IF

   END SUBROUTINE give_back_density_on_pool

END MODULE qs_vxc
