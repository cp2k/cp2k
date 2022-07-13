!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate EXX in RPA
!> \par History
!>      07.2020 separated from mp2.F [F. Stein, code by Jan Wilhelm]
!> \author Jan Wilhelm, Frederick Stein
! **************************************************************************************************
MODULE rpa_hfx
   USE admm_methods,                    ONLY: admm_projection_derivative
   USE admm_types,                      ONLY: admm_env_create,&
                                              admm_env_release,&
                                              admm_type,&
                                              get_admm_env
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE hfx_admm_utils,                  ONLY: create_admm_xc_section
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_ri,                          ONLY: hfx_ri_update_forces,&
                                              hfx_ri_update_ks
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none
   USE input_section_types,             ONLY: section_vals_create,&
                                              section_vals_duplicate,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_set_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_hfx'

   PUBLIC :: calculate_exx, calc_ri_rpa_admm_xc_contributions, ri_rpa_pre_hfx, ri_rpa_post_hfx

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param unit_nr ...
!> \param do_gw ...
!> \param do_admm ...
!> \param calc_forces ...
!> \param E_ex_from_GW ...
!> \param E_admm_from_GW ...
!> \param t3 ...
! **************************************************************************************************
   SUBROUTINE calculate_exx(qs_env, unit_nr, do_gw, do_admm, calc_forces, E_ex_from_GW, E_admm_from_GW, t3)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: do_gw, do_admm, calc_forces
      REAL(KIND=dp), INTENT(IN)                          :: E_ex_from_GW, E_admm_from_GW(2), t3

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_exx', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, irep, ispin, mspin, n_rep_hf, &
                                                            nspins
      LOGICAL                                            :: calc_ints, hfx_treat_lsd_in_core, &
                                                            use_virial
      REAL(KIND=dp)                                      :: eh1, ehfx, t1, t2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux_fit, rho_ao, &
                                                            rho_ao_aux_fit, rho_ao_resp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_2d, rho_ao_2d
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      NULLIFY (hfx_sections, input, para_env, matrix_ks, matrix_ks_aux_fit, rho, rho_ao, virial, &
               dft_control, rho_aux_fit, rho_ao_aux_fit)

      CALL ri_rpa_pre_hfx(qs_env)

      CALL get_qs_env(qs_env=qs_env, &
                      input=input, &
                      para_env=para_env, &
                      energy=energy, &
                      rho=rho, &
                      matrix_ks=matrix_ks, &
                      virial=virial, &
                      dft_control=dft_control)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      IF (do_admm) THEN
         CALL get_admm_env(qs_env%admm_env, matrix_ks_aux_fit=matrix_ks_aux_fit, rho_aux_fit=rho_aux_fit)
         CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux_fit)

         IF (qs_env%admm_env%do_gapw) THEN
            CPABORT("RPA with ADMM EXX only implmented with GPW")
         END IF
      END IF

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                i_rep_section=1)

      ! put matrix_ks to zero
      DO i = 1, SIZE(matrix_ks)
         CALL dbcsr_set(matrix_ks(i)%matrix, 0.0_dp)
         IF (do_admm) THEN
            CALL dbcsr_set(matrix_ks_aux_fit(i)%matrix, 0.0_dp)
         END IF
      END DO

      ! take the exact exchange energy from GW or calculate it
      IF (do_gw) THEN

         IF (calc_forces) CPABORT("Not implemented")

         IF (qs_env%mp2_env%ri_g0w0%update_xc_energy) THEN
            CALL remove_exc_energy(energy)
            energy%total = energy%total + E_ex_from_GW
            energy%ex = E_ex_from_GW
            IF (do_admm) THEN
               energy%total = energy%total + E_admm_from_GW(1) + E_admm_from_GW(2)
               energy%exc = E_admm_from_GW(1)
               energy%exc_aux_fit = E_admm_from_GW(2)
            END IF
            t2 = m_walltime()

            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.6)') 'Total EXX Time=', t2 - t1 + t3
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'EXX energy  =   ', energy%ex
            IF (do_admm .AND. unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') &
               'EXX ADMM XC correction  =   ', E_admm_from_GW(1) + E_admm_from_GW(2)
         END IF

      ELSE

         CALL remove_exc_energy(energy)

         nspins = dft_control%nspins
         mspin = 1
         IF (hfx_treat_lsd_in_core) mspin = nspins

         calc_ints = .TRUE.
         IF (qs_env%mp2_env%ri_rpa%reuse_hfx) calc_ints = .FALSE.
         IF (calc_forces) calc_ints = .FALSE.

         ehfx = 0.0_dp
         DO irep = 1, n_rep_hf
            IF (do_admm) THEN
               matrix_ks_2d(1:nspins, 1:1) => matrix_ks_aux_fit(1:nspins)
               rho_ao_2d(1:nspins, 1:1) => rho_ao_aux_fit(1:nspins)
            ELSE
               matrix_ks_2d(1:nspins, 1:1) => matrix_ks(1:nspins)
               rho_ao_2d(1:nspins, 1:1) => rho_ao(1:nspins)
            END IF

            IF (qs_env%mp2_env%ri_rpa%x_data(irep, 1)%do_hfx_ri) THEN
               CALL hfx_ri_update_ks(qs_env, qs_env%mp2_env%ri_rpa%x_data(irep, 1)%ri_data, matrix_ks_2d, ehfx, &
                                     rho_ao=rho_ao_2d, geometry_did_change=calc_ints, nspins=nspins, &
                                     hf_fraction=qs_env%mp2_env%ri_rpa%x_data(irep, 1)%general_parameter%fraction)
            ELSE

               DO ispin = 1, mspin
                  CALL integrate_four_center(qs_env, qs_env%mp2_env%ri_rpa%x_data, matrix_ks_2d, eh1, &
                                             rho_ao_2d, hfx_sections, para_env, calc_ints, irep, .TRUE., ispin=ispin)
                  ehfx = ehfx + eh1
               END DO
            END IF
         END DO

         ! include the EXX contribution to the total energy
         energy%ex = ehfx
         energy%total = energy%total + energy%ex

         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
         IF (use_virial) THEN
            virial%pv_calculate = .TRUE.
            virial%pv_fock_4c = 0.0_dp
         END IF

         DO irep = 1, n_rep_hf
            IF (calc_forces) THEN
               !Note: no need to remove xc forces: they are not even calculated in the first place
               NULLIFY (rho_ao_resp)
               IF (do_admm) THEN
                  rho_ao_2d(1:nspins, 1:1) => rho_ao_aux_fit(1:nspins)
               ELSE
                  rho_ao_2d(1:nspins, 1:1) => rho_ao(1:nspins)
               END IF

               IF (qs_env%mp2_env%ri_rpa%x_data(irep, 1)%do_hfx_ri) THEN

                  CALL hfx_ri_update_forces(qs_env, qs_env%mp2_env%ri_rpa%x_data(irep, 1)%ri_data, nspins, &
                                            qs_env%mp2_env%ri_rpa%x_data(irep, 1)%general_parameter%fraction, &
                                            rho_ao=rho_ao_2d, rho_ao_resp=rho_ao_resp, use_virial=use_virial)

               ELSE

                  CALL derivatives_four_center(qs_env, rho_ao_2d, rho_ao_resp, hfx_sections, para_env, irep, &
                                               use_virial, external_x_data=qs_env%mp2_env%ri_rpa%x_data)

               END IF

            END IF
         END DO !irep

         IF (use_virial) THEN
            virial%pv_exx = virial%pv_exx - virial%pv_fock_4c
            virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
         END IF

         ! ADMM XC correction
         IF (do_admm) THEN

            CALL calc_ri_rpa_admm_xc_contributions(qs_env, matrix_ks, matrix_ks_aux_fit, energy%exc, &
                                                   energy%exc_aux_fit, calc_forces, use_virial)

            ! ADMM overlap forces
            IF (calc_forces) CALL admm_projection_derivative(qs_env, matrix_ks_aux_fit, rho_ao)

            energy%total = energy%total + energy%exc_aux_fit
            energy%total = energy%total + energy%exc

         END IF

         IF (use_virial) THEN
            virial%pv_calculate = .FALSE.
         END IF

         t2 = m_walltime()

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.6)') 'Total EXX Time=', t2 - t1 + t3
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'EXX energy  =   ', energy%ex
         IF (do_admm .AND. unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,A,T56,F25.14)') 'EXX ADMM XC correction  =   ', energy%exc + energy%exc_aux_fit
         END IF
      END IF

      CALL ri_rpa_post_hfx(qs_env)

      CALL timestop(handle)

   END SUBROUTINE calculate_exx

! **************************************************************************************************
!> \brief get the ADMM XC section from the ri_rpa type if available, create and store them otherwise
!> \param qs_env ...
!> \param xc_section_aux ...
!> \param xc_section_primary ...
! **************************************************************************************************
   SUBROUTINE get_ri_rpa_admm_xc_sections(qs_env, xc_section_aux, xc_section_primary)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: xc_section_aux, xc_section_primary

      INTEGER                                            :: natom
      TYPE(admm_type), POINTER                           :: qs_admm_env, tmp_admm_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: xc_fun, xc_fun_empty, xc_section, &
                                                            xc_section_empty

      NULLIFY (qs_admm_env, tmp_admm_env, para_env, xc_section, xc_section_empty, xc_fun_empty, &
               xc_fun, dft_control)

      IF (ASSOCIATED(qs_env%mp2_env%ri_rpa%xc_section_aux) .AND. &
          ASSOCIATED(qs_env%mp2_env%ri_rpa%xc_section_primary)) THEN
         xc_section_aux => qs_env%mp2_env%ri_rpa%xc_section_aux
         xc_section_primary => qs_env%mp2_env%ri_rpa%xc_section_primary
      ELSE

         CALL get_qs_env(qs_env, admm_env=qs_admm_env, natom=natom, para_env=para_env, dft_control=dft_control)
         CPASSERT(ASSOCIATED(qs_admm_env))

         !create XC section with XC_FUNCITONAL NONE (aka empty XC_FUNCTIONAL section)
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
         xc_fun => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
         CALL section_vals_duplicate(xc_section, xc_section_empty)
         CALL section_vals_create(xc_fun_empty, xc_fun%section)
         CALL section_vals_set_subs_vals(xc_section_empty, "XC_FUNCTIONAL", xc_fun_empty)

         CALL admm_env_create(tmp_admm_env, dft_control%admm_control, qs_admm_env%mos_aux_fit, &
                              para_env, natom, qs_admm_env%nao_aux_fit)

         CALL create_admm_xc_section(x_data=qs_env%mp2_env%ri_rpa%x_data, xc_section=xc_section_empty, &
                                     admm_env=tmp_admm_env)

         CALL section_vals_duplicate(tmp_admm_env%xc_section_aux, xc_section_aux)
         CALL section_vals_duplicate(tmp_admm_env%xc_section_primary, xc_section_primary)

         qs_env%mp2_env%ri_rpa%xc_section_aux => xc_section_aux
         qs_env%mp2_env%ri_rpa%xc_section_primary => xc_section_primary

         CALL section_vals_release(xc_section_empty)
         CALL section_vals_release(xc_fun_empty)
         CALL admm_env_release(tmp_admm_env)

      END IF

   END SUBROUTINE get_ri_rpa_admm_xc_sections

! **************************************************************************************************
!> \brief Calculate the RI_RPA%HF ADMM XC contributions to the KS matrices and the respective energies
!> \param qs_env ...
!> \param matrix_prim ...
!> \param matrix_aux ...
!> \param exc ...
!> \param exc_aux_fit ...
!> \param calc_forces ...
!> \param use_virial ...
! **************************************************************************************************
   SUBROUTINE calc_ri_rpa_admm_xc_contributions(qs_env, matrix_prim, matrix_aux, exc, exc_aux_fit, &
                                                calc_forces, use_virial)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_prim, matrix_aux
      REAL(dp), INTENT(INOUT)                            :: exc, exc_aux_fit
      LOGICAL, INTENT(IN)                                :: calc_forces, use_virial

      INTEGER                                            :: ispin, nspins
      REAL(dp), DIMENSION(3, 3)                          :: pv_loc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_ao_aux_fit
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_dummy, v_rspace, v_rspace_aux_fit
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: xc_section_aux, xc_section_primary
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (xc_section_aux, xc_section_primary, rho, rho_aux_fit, v_dummy, v_rspace, v_rspace_aux_fit, &
               auxbas_pw_pool, pw_env, rho_ao, rho_ao_aux_fit, dft_control, admm_env)

      CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env, rho=rho, admm_env=admm_env, virial=virial)
      CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit)

      nspins = dft_control%nspins
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux_fit)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL get_ri_rpa_admm_xc_sections(qs_env, xc_section_aux, xc_section_primary)

      IF (use_virial) virial%pv_xc = 0.0_dp
      CALL qs_vxc_create(qs_env%ks_env, rho_struct=rho_aux_fit, xc_section=xc_section_aux, &
                         vxc_rho=v_rspace_aux_fit, vxc_tau=v_dummy, exc=exc_aux_fit)
      IF (use_virial) THEN
         virial%pv_exc = virial%pv_exc - virial%pv_xc
         virial%pv_virial = virial%pv_virial - virial%pv_xc
      END IF

      IF (.NOT. dft_control%admm_control%aux_exch_func == do_admm_aux_exch_func_none) THEN
         IF (use_virial) pv_loc = virial%pv_virial
         DO ispin = 1, nspins
            CALL pw_scale(v_rspace_aux_fit(ispin)%pw, v_rspace_aux_fit(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_rspace_aux_fit(ispin), hmat=matrix_aux(ispin), &
                                    pmat=rho_ao_aux_fit(ispin), qs_env=qs_env, &
                                    basis_type="AUX_FIT", calculate_forces=calc_forces, &
                                    task_list_external=qs_env%admm_env%task_list_aux_fit)
         END DO
         IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
      END IF

      IF (ASSOCIATED(v_rspace_aux_fit)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_aux_fit(ispin)%pw)
         END DO
         DEALLOCATE (v_rspace_aux_fit)
      END IF
      IF (ASSOCIATED(v_dummy)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_dummy(ispin)%pw)
         END DO
         DEALLOCATE (v_dummy)
      END IF

      IF (use_virial) virial%pv_xc = 0.0_dp
      CALL qs_vxc_create(qs_env%ks_env, rho_struct=rho, xc_section=xc_section_primary, &
                         vxc_rho=v_rspace, vxc_tau=v_dummy, exc=exc)
      IF (use_virial) THEN
         virial%pv_exc = virial%pv_exc - virial%pv_xc
         virial%pv_virial = virial%pv_virial - virial%pv_xc
      END IF

      IF (.NOT. dft_control%admm_control%aux_exch_func == do_admm_aux_exch_func_none) THEN
         IF (use_virial) pv_loc = virial%pv_virial
         DO ispin = 1, nspins
            CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_rspace(ispin), hmat=matrix_prim(ispin), &
                                    pmat=rho_ao(ispin), qs_env=qs_env, &
                                    calculate_forces=calc_forces)
         END DO
         IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
      END IF

      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
         END DO
         DEALLOCATE (v_rspace)
      END IF
      IF (ASSOCIATED(v_dummy)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_dummy(ispin)%pw)
         END DO
         DEALLOCATE (v_dummy)
      END IF

   END SUBROUTINE calc_ri_rpa_admm_xc_contributions

! **************************************************************************************************
!> \brief Prepare the ri_rpa%x_data for integration. Simply change the HFX fraction in case the
!>        qs_env%x_data is reused
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE ri_rpa_pre_hfx(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: irep, n_rep_hf
      REAL(dp)                                           :: frac
      TYPE(section_vals_type), POINTER                   :: input, rpa_hfx_section

      IF (.NOT. qs_env%mp2_env%ri_rpa%reuse_hfx) RETURN

      CALL get_qs_env(qs_env, input=input)
      rpa_hfx_section => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
      CALL section_vals_get(rpa_hfx_section, n_repetition=n_rep_hf)

      DO irep = 1, n_rep_hf
         CALL section_vals_val_get(rpa_hfx_section, "FRACTION", r_val=frac, i_rep_section=irep)
         qs_env%mp2_env%ri_rpa%x_data(irep, :)%general_parameter%fraction = frac
      END DO

   END SUBROUTINE ri_rpa_pre_hfx

! **************************************************************************************************
!> \brief Revert back to the proper HFX fraction in case qs_env%x_data is reused
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE ri_rpa_post_hfx(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: irep, n_rep_hf
      REAL(dp)                                           :: frac
      TYPE(section_vals_type), POINTER                   :: input, qs_hfx_section

      IF (.NOT. qs_env%mp2_env%ri_rpa%reuse_hfx) RETURN

      CALL get_qs_env(qs_env, input=input)
      qs_hfx_section => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(qs_hfx_section, n_repetition=n_rep_hf)

      DO irep = 1, n_rep_hf
         CALL section_vals_val_get(qs_hfx_section, "FRACTION", r_val=frac, i_rep_section=irep)
         qs_env%mp2_env%ri_rpa%x_data(irep, :)%general_parameter%fraction = frac
      END DO

   END SUBROUTINE ri_rpa_post_hfx

! **************************************************************************************************
!> \brief ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE remove_exc_energy(energy)
      TYPE(qs_energy_type), INTENT(INOUT)                :: energy

      CHARACTER(len=*), PARAMETER :: routineN = 'remove_exc_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Remove the Exchange-correlation energy contributions from the total energy
      energy%total = energy%total - (energy%exc + energy%exc1 + energy%ex + &
                                     energy%exc_aux_fit + energy%exc1_aux_fit)

      energy%exc = 0.0_dp
      energy%exc1 = 0.0_dp
      energy%exc_aux_fit = 0.0_dp
      energy%exc1_aux_fit = 0.0_dp
      energy%ex = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE

END MODULE rpa_hfx

