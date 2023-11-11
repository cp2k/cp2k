!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate EXX in RPA and energy correction methods
!> \par History
!>      07.2020 separated from mp2.F [F. Stein, code by Jan Wilhelm]
!>      06.2022 EXX contribution to the forces [A. Bussy]
!>      03.2023 Generalized for energy correction methods
!> \author Jan Wilhelm, Frederick Stein, Augustin Bussy, Fabian Belleflamme
! **************************************************************************************************
MODULE hfx_exx
   USE admm_methods,                    ONLY: admm_projection_derivative
   USE admm_types,                      ONLY: admm_env_create,&
                                              admm_env_release,&
                                              admm_type,&
                                              get_admm_env
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_scale,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE hfx_admm_utils,                  ONLY: create_admm_xc_section,&
                                              tddft_hfx_matrix
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_ri,                          ONLY: hfx_ri_update_forces,&
                                              hfx_ri_update_ks
   USE hfx_types,                       ONLY: hfx_type
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
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_ks_reference,                 ONLY: ks_ref_potential
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE task_list_types,                 ONLY: task_list_type
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_exx'

   PUBLIC :: calculate_exx, add_exx_to_rhs, calc_exx_admm_xc_contributions, exx_pre_hfx, exx_post_hfx

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param unit_nr ...
!> \param hfx_sections ...
!> \param x_data ...
!> \param do_gw ...
!> \param do_admm ...
!> \param calc_forces ...
!> \param reuse_hfx ...
!> \param do_im_time ...
!> \param E_ex_from_GW ...
!> \param E_admm_from_GW ...
!> \param t3 ...
! **************************************************************************************************
   SUBROUTINE calculate_exx(qs_env, unit_nr, hfx_sections, x_data, &
                            do_gw, do_admm, calc_forces, reuse_hfx, do_im_time, &
                            E_ex_from_GW, E_admm_from_GW, t3)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(section_vals_type), POINTER                   :: hfx_sections
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      LOGICAL, INTENT(IN)                                :: do_gw, do_admm, calc_forces, reuse_hfx, &
                                                            do_im_time
      REAL(KIND=dp), INTENT(IN)                          :: E_ex_from_GW, E_admm_from_GW(2), t3

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calculate_exx'

      INTEGER                                            :: handle, i, irep, ispin, mspin, n_rep_hf, &
                                                            nspins
      LOGICAL                                            :: calc_ints, hfx_treat_lsd_in_core, &
                                                            use_virial
      REAL(KIND=dp)                                      :: eh1, ehfx, t1, t2, tf1, tf2
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux_fit, rho_ao, &
                                                            rho_ao_aux_fit, rho_ao_resp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_2d, rho_ao_2d
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      NULLIFY (input, para_env, matrix_ks, matrix_ks_aux_fit, rho, rho_ao, virial, &
               dft_control, rho_aux_fit, rho_ao_aux_fit)

      CALL exx_pre_hfx(hfx_sections, x_data, reuse_hfx)

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
         CALL get_admm_env(qs_env%admm_env, &
                           matrix_ks_aux_fit=matrix_ks_aux_fit, &
                           rho_aux_fit=rho_aux_fit)
         CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux_fit)

         IF (qs_env%admm_env%do_gapw) THEN
            CPABORT("ADMM EXX only implmented with GPW")
         END IF
      END IF

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
         IF (reuse_hfx) calc_ints = .FALSE.
         IF (calc_forces .AND. do_im_time) calc_ints = .FALSE.

         ehfx = 0.0_dp
         IF (do_admm) THEN
            matrix_ks_2d(1:nspins, 1:1) => matrix_ks_aux_fit(1:nspins)
            rho_ao_2d(1:nspins, 1:1) => rho_ao_aux_fit(1:nspins)
         ELSE
            matrix_ks_2d(1:nspins, 1:1) => matrix_ks(1:nspins)
            rho_ao_2d(1:nspins, 1:1) => rho_ao(1:nspins)
         END IF

         DO irep = 1, n_rep_hf

            IF (x_data(irep, 1)%do_hfx_ri) THEN
               CALL hfx_ri_update_ks(qs_env, x_data(irep, 1)%ri_data, matrix_ks_2d, ehfx, &
                                     rho_ao=rho_ao_2d, geometry_did_change=calc_ints, nspins=nspins, &
                                     hf_fraction=x_data(irep, 1)%general_parameter%fraction)
            ELSE

               DO ispin = 1, mspin

                  CALL integrate_four_center(qs_env, x_data, matrix_ks_2d, eh1, &
                                             rho_ao_2d, hfx_sections, para_env, &
                                             calc_ints, irep, .TRUE., ispin=ispin)
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

         IF (calc_forces) THEN
            tf1 = m_walltime()

            !Note: no need to remove xc forces: they are not even calculated in the first place
            NULLIFY (rho_ao_resp)

            DO irep = 1, n_rep_hf

               IF (x_data(irep, 1)%do_hfx_ri) THEN

                  CALL hfx_ri_update_forces(qs_env, x_data(irep, 1)%ri_data, nspins, &
                                            x_data(irep, 1)%general_parameter%fraction, &
                                            rho_ao=rho_ao_2d, rho_ao_resp=rho_ao_resp, &
                                            use_virial=use_virial)
               ELSE

                  CALL derivatives_four_center(qs_env, rho_ao_2d, rho_ao_resp, &
                                               hfx_sections, para_env, irep, &
                                               use_virial, external_x_data=x_data)

               END IF

            END DO !irep

            tf2 = m_walltime()
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.6)') 'Total EXX Force Time=', tf2 - tf1
         END IF

         IF (use_virial) THEN
            virial%pv_exx = virial%pv_exx - virial%pv_fock_4c
            virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
         END IF

         ! ADMM XC correction
         IF (do_admm) THEN

            CALL calc_exx_admm_xc_contributions(qs_env, matrix_ks, matrix_ks_aux_fit, x_data, &
                                                energy%exc, energy%exc_aux_fit, calc_forces, &
                                                use_virial)

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

      CALL exx_post_hfx(qs_env, x_data, reuse_hfx)

      CALL timestop(handle)

   END SUBROUTINE calculate_exx

! **************************************************************************************************
!> \brief Add the EXX contribution to the RHS of the Z-vector equation, namely the HF Hamiltonian
!> \param rhs ...
!> \param qs_env ...
!> \param ext_hfx_section ...
!> \param x_data ...
!> \param recalc_integrals ...
!> \param do_admm ...
!> \param do_ec ...
!> \param do_exx ...
!> \param reuse_hfx ...
! **************************************************************************************************
   SUBROUTINE add_exx_to_rhs(rhs, qs_env, ext_hfx_section, x_data, &
                             recalc_integrals, do_admm, do_ec, do_exx, reuse_hfx)

      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: rhs
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: ext_hfx_section
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      LOGICAL, INTENT(IN), OPTIONAL                      :: recalc_integrals, do_admm, do_ec, &
                                                            do_exx, reuse_hfx

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'add_exx_to_rhs'

      INTEGER                                            :: handle, ispin, nao, nao_aux, nspins, &
                                                            unit_nr
      LOGICAL                                            :: calc_ints, do_hfx, my_do_ec, my_do_exx, &
                                                            my_recalc_integrals
      REAL(dp)                                           :: dummy_real1, dummy_real2
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_work, matrix_ks, matrix_s_aux, &
                                                            rho_ao, rho_ao_aux, work_admm
      TYPE(dbcsr_type)                                   :: dbcsr_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: vh_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: vadmm_rspace, vtau_rspace, vxc_rspace
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit

      NULLIFY (dbcsr_work, matrix_ks, rho, hfx_section, pw_env, vadmm_rspace, vtau_rspace, vxc_rspace, &
               auxbas_pw_pool, dft_control, rho_ao, rho_aux_fit, rho_ao_aux, work_admm, &
               matrix_s_aux, admm_env, task_list_aux_fit)

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      CALL timeset(routineN, handle)

      my_recalc_integrals = .FALSE.
      IF (PRESENT(recalc_integrals)) my_recalc_integrals = recalc_integrals
      my_do_ec = .FALSE.
      IF (PRESENT(do_ec)) my_do_ec = do_ec
      my_do_exx = .TRUE.
      IF (PRESENT(do_exx)) my_do_exx = do_exx

      ! Strategy: we take the ks_matrix, remove the current xc contribution, and then add the RPA HF one
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, rho=rho, pw_env=pw_env, dft_control=dft_control)
      nspins = dft_control%nspins

      ! do_exx ; subtract XC and EX
      ! do_dcdft: subtract EX

      CALL dbcsr_allocate_matrix_set(dbcsr_work, nspins)
      DO ispin = 1, nspins
         ALLOCATE (dbcsr_work(ispin)%matrix)
         CALL dbcsr_copy(dbcsr_work(ispin)%matrix, matrix_ks(ispin)%matrix)
         CALL dbcsr_set(dbcsr_work(ispin)%matrix, 0.0_dp)
      END DO

      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL get_admm_env(admm_env, matrix_s_aux_fit=matrix_s_aux, task_list_aux_fit=task_list_aux_fit, &
                           rho_aux_fit=rho_aux_fit)
         CALL dbcsr_allocate_matrix_set(work_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (work_admm(ispin)%matrix)
            CALL dbcsr_create(work_admm(ispin)%matrix, template=matrix_s_aux(1)%matrix)
            CALL dbcsr_copy(work_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
            CALL dbcsr_set(work_admm(ispin)%matrix, 0.0_dp)
         END DO
         CALL dbcsr_create(dbcsr_tmp, template=matrix_ks(1)%matrix)
         CALL dbcsr_copy(dbcsr_tmp, matrix_ks(1)%matrix)

         nao = admm_env%nao_orb
         nao_aux = admm_env%nao_aux_fit
         CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux)
      END IF

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ! Remove the standard XC + HFX contribution
      CALL ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, dummy_real1, dummy_real2)

      ! Only remove standard XC and/or HFX
      ! (due to different requirements for RHS)
      IF (my_do_exx) THEN

         DO ispin = 1, nspins
            CALL dbcsr_copy(dbcsr_work(ispin)%matrix, matrix_ks(ispin)%matrix)
         END DO

         DO ispin = 1, nspins
            CALL pw_scale(vxc_rspace(ispin), -1.0_dp)
            CALL integrate_v_rspace(v_rspace=vxc_rspace(ispin), hmat=dbcsr_work(ispin), qs_env=qs_env, &
                                    calculate_forces=.FALSE.)
            IF (ASSOCIATED(vtau_rspace)) THEN
               CALL pw_scale(vtau_rspace(ispin), -1.0_dp)
               CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), hmat=dbcsr_work(ispin), qs_env=qs_env, &
                                       calculate_forces=.FALSE., compute_tau=.TRUE.)
            END IF
         END DO

      END IF ! do_exx

      ! Remove standard HFX
      IF (my_do_exx .OR. my_do_ec) THEN

         hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF")
         CALL section_vals_get(hfx_section, explicit=do_hfx)
         IF (do_hfx) THEN
            IF (dft_control%do_admm) THEN

               !note: factor -1.0 taken care of later, after HFX ADMM contribution is taken
               IF (.NOT. qs_env%admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
                  DO ispin = 1, nspins
                     CALL integrate_v_rspace(v_rspace=vadmm_rspace(ispin), hmat=work_admm(ispin), &
                                             qs_env=qs_env, calculate_forces=.FALSE., basis_type="AUX_FIT", &
                                             task_list_external=task_list_aux_fit)
                  END DO
               END IF

               CALL tddft_hfx_matrix(work_admm, rho_ao_aux, qs_env, .FALSE., my_recalc_integrals)

               DO ispin = 1, nspins
                  CALL copy_dbcsr_to_fm(work_admm(ispin)%matrix, admm_env%work_aux_aux)
                  CALL parallel_gemm('N', 'N', nao_aux, nao, nao_aux, 1.0_dp, admm_env%work_aux_aux, admm_env%A, &
                                     0.0_dp, admm_env%work_aux_orb)
                  CALL parallel_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_env%A, admm_env%work_aux_orb, &
                                     0.0_dp, admm_env%work_orb_orb)
                  CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbcsr_tmp, keep_sparsity=.TRUE.)
                  CALL dbcsr_add(dbcsr_work(ispin)%matrix, dbcsr_tmp, 1.0_dp, -1.0_dp)
               END DO
            ELSE
               DO ispin = 1, nspins
                  CALL dbcsr_scale(rho_ao(ispin)%matrix, -1.0_dp)
               END DO
               CALL tddft_hfx_matrix(dbcsr_work, rho_ao, qs_env, .FALSE., my_recalc_integrals)
               DO ispin = 1, nspins
                  CALL dbcsr_scale(rho_ao(ispin)%matrix, -1.0_dp)
               END DO
            END IF
         END IF !do_hfx

      END IF ! do_exx

      ! Clean
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL auxbas_pw_pool%give_back_pw(vh_rspace)
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%give_back_pw(vxc_rspace(ispin))
         IF (ASSOCIATED(vtau_rspace)) THEN
            CALL auxbas_pw_pool%give_back_pw(vtau_rspace(ispin))
         END IF
      END DO
      DEALLOCATE (vxc_rspace)
      IF (ASSOCIATED(vtau_rspace)) DEALLOCATE (vtau_rspace)

      IF (dft_control%do_admm) THEN
         DO ispin = 1, nspins
            IF (ASSOCIATED(vadmm_rspace)) THEN
               CALL auxbas_pw_pool%give_back_pw(vadmm_rspace(ispin))
            END IF
         END DO
         IF (ASSOCIATED(vadmm_rspace)) DEALLOCATE (vadmm_rspace)
      END IF

      !Add the HF contribution from RI_RPA/EC_ENV to the ks matrix
      CALL exx_pre_hfx(ext_hfx_section, x_data, reuse_hfx)

      calc_ints = .TRUE.
      IF (reuse_hfx) calc_ints = .FALSE.
      IF (do_admm) THEN
         DO ispin = 1, nspins
            CALL dbcsr_set(work_admm(ispin)%matrix, 0.0_dp)
         END DO

         CALL tddft_hfx_matrix(work_admm, rho_ao_aux, qs_env, .FALSE., calc_ints, ext_hfx_section, &
                               x_data)

         !ADMM XC correction
         CALL calc_exx_admm_xc_contributions(qs_env, dbcsr_work, work_admm, x_data, dummy_real1, &
                                             dummy_real2, .FALSE., .FALSE.)

         DO ispin = 1, nspins
            CALL copy_dbcsr_to_fm(work_admm(ispin)%matrix, admm_env%work_aux_aux)
            CALL parallel_gemm('N', 'N', nao_aux, nao, nao_aux, 1.0_dp, admm_env%work_aux_aux, admm_env%A, &
                               0.0_dp, admm_env%work_aux_orb)
            CALL parallel_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_env%A, admm_env%work_aux_orb, &
                               0.0_dp, admm_env%work_orb_orb)
            CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbcsr_tmp, keep_sparsity=.TRUE.)
            CALL dbcsr_add(dbcsr_work(ispin)%matrix, dbcsr_tmp, 1.0_dp, 1.0_dp)
         END DO

      ELSE
         CALL tddft_hfx_matrix(dbcsr_work, rho_ao, qs_env, .FALSE., calc_ints, ext_hfx_section, x_data)
      END IF
      CALL exx_post_hfx(qs_env, x_data, reuse_hfx)

      !Update the RHS
      DO ispin = 1, nspins
         CALL dbcsr_add(rhs(ispin)%matrix, dbcsr_work(ispin)%matrix, 1.0_dp, 1.0_dp)
      END DO

      CALL dbcsr_deallocate_matrix_set(dbcsr_work)
      IF (dft_control%do_admm) THEN
         CALL dbcsr_deallocate_matrix_set(work_admm)
         CALL dbcsr_release(dbcsr_tmp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE add_exx_to_rhs

! **************************************************************************************************
!> \brief get the ADMM XC section from the ri_rpa/ec_env type if available, create and store them otherwise
!> \param qs_env ...
!> \param x_data ...
!> \param xc_section_aux ...
!> \param xc_section_primary ...
! **************************************************************************************************
   SUBROUTINE get_exx_admm_xc_sections(qs_env, x_data, xc_section_aux, xc_section_primary)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(section_vals_type), POINTER                   :: xc_section_aux, xc_section_primary

      INTEGER                                            :: natom
      TYPE(admm_type), POINTER                           :: qs_admm_env, tmp_admm_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: xc_fun, xc_fun_empty, xc_section, &
                                                            xc_section_empty

      NULLIFY (qs_admm_env, tmp_admm_env, para_env, xc_section, xc_section_empty, xc_fun_empty, &
               xc_fun, dft_control)

      IF (ASSOCIATED(qs_env%mp2_env)) THEN
         IF (ASSOCIATED(qs_env%mp2_env%ri_rpa%xc_section_aux) .AND. &
             ASSOCIATED(qs_env%mp2_env%ri_rpa%xc_section_primary)) THEN
            xc_section_aux => qs_env%mp2_env%ri_rpa%xc_section_aux
            xc_section_primary => qs_env%mp2_env%ri_rpa%xc_section_primary
         END IF
      ELSEIF (qs_env%energy_correction) THEN
         IF (ASSOCIATED(qs_env%ec_env%xc_section_aux) .AND. &
             ASSOCIATED(qs_env%ec_env%xc_section_primary)) THEN
            xc_section_aux => qs_env%ec_env%xc_section_aux
            xc_section_primary => qs_env%ec_env%xc_section_primary
         END IF
      END IF

      IF (.NOT. ASSOCIATED(xc_section_aux) .OR. .NOT. ASSOCIATED(xc_section_primary)) THEN

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

         CALL create_admm_xc_section(x_data=x_data, xc_section=xc_section_empty, &
                                     admm_env=tmp_admm_env)

         CALL section_vals_duplicate(tmp_admm_env%xc_section_aux, xc_section_aux)
         CALL section_vals_duplicate(tmp_admm_env%xc_section_primary, xc_section_primary)

         IF (ASSOCIATED(qs_env%mp2_env)) THEN
            qs_env%mp2_env%ri_rpa%xc_section_aux => xc_section_aux
            qs_env%mp2_env%ri_rpa%xc_section_primary => xc_section_primary
         END IF
         IF (qs_env%energy_correction) THEN
            qs_env%ec_env%xc_section_aux => xc_section_aux
            qs_env%ec_env%xc_section_primary => xc_section_primary
         END IF

         CALL section_vals_release(xc_section_empty)
         CALL section_vals_release(xc_fun_empty)
         CALL admm_env_release(tmp_admm_env)

      END IF

   END SUBROUTINE get_exx_admm_xc_sections

! **************************************************************************************************
!> \brief Calculate the RI_RPA%HF / EC_ENV%HF ADMM XC contributions to the KS matrices and the respective energies
!> \param qs_env ...
!> \param matrix_prim ...
!> \param matrix_aux ...
!> \param x_data ...
!> \param exc ...
!> \param exc_aux_fit ...
!> \param calc_forces ...
!> \param use_virial ...
! **************************************************************************************************
   SUBROUTINE calc_exx_admm_xc_contributions(qs_env, matrix_prim, matrix_aux, x_data, exc, exc_aux_fit, &
                                             calc_forces, use_virial)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_prim, matrix_aux
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      REAL(dp), INTENT(INOUT)                            :: exc, exc_aux_fit
      LOGICAL, INTENT(IN)                                :: calc_forces, use_virial

      INTEGER                                            :: ispin, nspins
      REAL(dp), DIMENSION(3, 3)                          :: pv_loc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_ao_aux_fit
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), DIMENSION(:), POINTER               :: v_dummy, v_rspace, v_rspace_aux_fit
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

      CALL get_exx_admm_xc_sections(qs_env, x_data, xc_section_aux, xc_section_primary)

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
            CALL pw_scale(v_rspace_aux_fit(ispin), v_rspace_aux_fit(ispin)%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_rspace_aux_fit(ispin), hmat=matrix_aux(ispin), &
                                    pmat=rho_ao_aux_fit(ispin), qs_env=qs_env, &
                                    basis_type="AUX_FIT", calculate_forces=calc_forces, &
                                    task_list_external=qs_env%admm_env%task_list_aux_fit)
         END DO
         IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
      END IF

      IF (ASSOCIATED(v_rspace_aux_fit)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(v_rspace_aux_fit(ispin))
         END DO
         DEALLOCATE (v_rspace_aux_fit)
      END IF
      IF (ASSOCIATED(v_dummy)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(v_dummy(ispin))
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
            CALL pw_scale(v_rspace(ispin), v_rspace(ispin)%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_rspace(ispin), hmat=matrix_prim(ispin), &
                                    pmat=rho_ao(ispin), qs_env=qs_env, &
                                    calculate_forces=calc_forces)
         END DO
         IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
      END IF

      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(v_rspace(ispin))
         END DO
         DEALLOCATE (v_rspace)
      END IF
      IF (ASSOCIATED(v_dummy)) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(v_dummy(ispin))
         END DO
         DEALLOCATE (v_dummy)
      END IF

   END SUBROUTINE calc_exx_admm_xc_contributions

! **************************************************************************************************
!> \brief Prepare the external x_data for integration. Simply change the HFX fraction in case the
!>        qs_env%x_data is reused
!> \param ext_hfx_section ...
!> \param x_data ...
!> \param reuse_hfx ...
! **************************************************************************************************
   SUBROUTINE exx_pre_hfx(ext_hfx_section, x_data, reuse_hfx)
      TYPE(section_vals_type), POINTER                   :: ext_hfx_section
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      LOGICAL                                            :: reuse_hfx

      INTEGER                                            :: irep, n_rep_hf
      REAL(dp)                                           :: frac

      IF (.NOT. reuse_hfx) RETURN

      CALL section_vals_get(ext_hfx_section, n_repetition=n_rep_hf)

      DO irep = 1, n_rep_hf
         CALL section_vals_val_get(ext_hfx_section, "FRACTION", r_val=frac, i_rep_section=irep)
         x_data(irep, :)%general_parameter%fraction = frac
      END DO

   END SUBROUTINE exx_pre_hfx

! **************************************************************************************************
!> \brief Revert back to the proper HFX fraction in case qs_env%x_data is reused
!> \param qs_env ...
!> \param x_data ...
!> \param reuse_hfx ...
! **************************************************************************************************
   SUBROUTINE exx_post_hfx(qs_env, x_data, reuse_hfx)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      LOGICAL                                            :: reuse_hfx

      INTEGER                                            :: irep, n_rep_hf
      REAL(dp)                                           :: frac
      TYPE(section_vals_type), POINTER                   :: input, qs_hfx_section

      IF (.NOT. reuse_hfx) RETURN

      CALL get_qs_env(qs_env, input=input)
      qs_hfx_section => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(qs_hfx_section, n_repetition=n_rep_hf)

      DO irep = 1, n_rep_hf
         CALL section_vals_val_get(qs_hfx_section, "FRACTION", r_val=frac, i_rep_section=irep)
         x_data(irep, :)%general_parameter%fraction = frac
      END DO

   END SUBROUTINE exx_post_hfx

! **************************************************************************************************
!> \brief ...
!> \param energy ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE remove_exc_energy(energy)
      TYPE(qs_energy_type), INTENT(INOUT)                :: energy

      ! Remove the Exchange-correlation energy contributions from the total energy
      energy%total = energy%total - (energy%exc + energy%exc1 + energy%ex + &
                                     energy%exc_aux_fit + energy%exc1_aux_fit)

      energy%exc = 0.0_dp
      energy%exc1 = 0.0_dp
      energy%exc_aux_fit = 0.0_dp
      energy%exc1_aux_fit = 0.0_dp
      energy%ex = 0.0_dp

   END SUBROUTINE

END MODULE hfx_exx

