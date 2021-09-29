!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate EXX within GW
!> \par History
!>      07.2020 separated from mp2.F [F. Stein, code by Jan Wilhelm]
!> \author Jan Wilhelm, Frederick Stein
! **************************************************************************************************
MODULE rpa_gw_sigma
   USE admm_methods,                    ONLY: admm_mo_merge_ks_matrix
   USE admm_types,                      ONLY: admm_type
   USE admm_utils,                      ONLY: admm_correct_for_eigenvalues
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_get_diag, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_release_p, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_symmetric
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE input_constants,                 ONLY: do_admm_basis_projection,&
                                              do_admm_purify_none,&
                                              gw_print_exx,&
                                              gw_read_exx,&
                                              xc_none
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_sum
   USE mp2_types,                       ONLY: mp2_type
   USE physcon,                         ONLY: evolt
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_sigma'

   PUBLIC :: compute_vec_Sigma_x_minus_vxc_gw

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param mos_mp2 ...
!> \param energy_ex ...
!> \param t3 ...
!> \param unit_nr ...
!> \par History
!>      04.2015 created
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE compute_vec_Sigma_x_minus_vxc_gw(qs_env, mp2_env, mos_mp2, energy_ex, t3, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos_mp2
      REAL(KIND=dp), INTENT(OUT)                         :: energy_ex, t3
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_vec_Sigma_x_minus_vxc_gw', &
         routineP = moduleN//':'//routineN

      CHARACTER(4)                                       :: occ_virt
      CHARACTER(LEN=40)                                  :: line
      INTEGER :: dimen, gw_corr_lev_occ, gw_corr_lev_virt, handle, homo, i_img, ikp, irep, ispin, &
         iunit, myfun, myfun_aux, myfun_prim, n_level_gw, n_level_gw_ref, n_rep_hf, nkp, nmo, ns, &
         nspins, print_exx
      LOGICAL                                            :: charge_constrain_tmp, do_admm_rpa, &
                                                            do_hfx, do_kpoints_cubic_RPA, &
                                                            do_ri_Sigma_x, really_read_line
      REAL(KIND=dp) :: eh1, ehfx, energy_exc, energy_exc1, energy_exc1_aux_fit, &
         energy_exc_aux_fit, energy_total, exx_minus_vxc, hfx_fraction, t1, t2, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: vec_Sigma_x_minus_vxc_gw, &
                                                            vec_Sigma_x_minus_vxc_gw_im
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux_fit, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER :: matrix_ks_2d, matrix_ks_kp_im, &
         matrix_ks_kp_re, matrix_ks_transl, matrix_sigma_x_minus_vxc, matrix_sigma_x_minus_vxc_im, &
         rho_ao_2d
      TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_tmp_2, mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input, xc_section, &
                                                            xc_section_admm_aux, &
                                                            xc_section_admm_prim

      NULLIFY (admm_env, matrix_ks, matrix_ks_aux_fit, rho_ao, matrix_sigma_x_minus_vxc, input, &
               xc_section, xc_section_admm_aux, xc_section_admm_prim, hfx_sections, rho, &
               dft_control, para_env, ks_env, mo_coeff, matrix_sigma_x_minus_vxc_im)

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      do_admm_rpa = mp2_env%ri_rpa%do_admm
      do_ri_Sigma_x = mp2_env%ri_g0w0%do_ri_Sigma_x
      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      IF (do_kpoints_cubic_RPA) THEN
         CPASSERT(do_ri_Sigma_x)
      END IF

      ! In case we do admm for RPA/GW, get the aux. density matrix rho_aux_fit, otherwise
      ! just the normal density matrix rho
      IF (do_admm_rpa) THEN
         CALL get_qs_env(qs_env, &
                         admm_env=admm_env, &
                         matrix_ks=matrix_ks, &
                         matrix_ks_aux_fit=matrix_ks_aux_fit, &
                         rho_aux_fit=rho, &
                         input=input, &
                         dft_control=dft_control, &
                         para_env=para_env, &
                         ks_env=ks_env, &
                         energy=energy)

         ! RPA/GW with ADMM for EXX or the exchange self-energy only implemented for
         ! ADMM_PURIFICATION_METHOD  NONE
         ! METHOD                    BASIS_PROJECTION
         ! in the admm section
         CPASSERT(admm_env%purification_method == do_admm_purify_none)
         CPASSERT(dft_control%admm_control%method == do_admm_basis_projection)

      ELSE

         IF (do_kpoints_cubic_RPA) THEN

            CALL get_qs_env(qs_env, &
                            admm_env=admm_env, &
                            matrix_ks_kp=matrix_ks_transl, &
                            rho=rho, &
                            input=input, &
                            dft_control=dft_control, &
                            para_env=para_env, &
                            kpoints=kpoints, &
                            ks_env=ks_env, &
                            energy=energy)

            nkp = kpoints%nkp

         ELSE

            CALL get_qs_env(qs_env, &
                            admm_env=admm_env, &
                            matrix_ks=matrix_ks, &
                            rho=rho, &
                            input=input, &
                            dft_control=dft_control, &
                            para_env=para_env, &
                            ks_env=ks_env, &
                            energy=energy)

            nkp = 1

         END IF

      END IF

      nspins = dft_control%nspins

      IF (.NOT. do_kpoints_cubic_RPA) THEN
         CALL qs_rho_get(rho, rho_ao=rho_ao)
      END IF

      ! If ADMM we should make the ks matrix up-to-date
      IF (dft_control%do_admm) THEN
      DO ispin = 1, SIZE(matrix_ks)
         CALL admm_correct_for_eigenvalues(ispin, admm_env, matrix_ks(ispin)%matrix)
      END DO
      END IF

      IF (do_kpoints_cubic_RPA) THEN

         CALL allocate_matrix_ks_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)
         CALL transform_matrix_ks_to_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)

         DO ispin = 1, nspins
         DO i_img = 1, SIZE(matrix_ks_transl, 2)
            CALL dbcsr_set(matrix_ks_transl(ispin, i_img)%matrix, 0.0_dp)
         END DO
         END DO

      END IF

      ! initialize matrix_sigma_x_minus_vxc
      NULLIFY (matrix_sigma_x_minus_vxc)
      CALL dbcsr_allocate_matrix_set(matrix_sigma_x_minus_vxc, nspins, nkp)
      IF (do_kpoints_cubic_RPA) THEN
         NULLIFY (matrix_sigma_x_minus_vxc_im)
         CALL dbcsr_allocate_matrix_set(matrix_sigma_x_minus_vxc_im, nspins, nkp)
      END IF

      DO ispin = 1, nspins
         DO ikp = 1, nkp

            IF (do_kpoints_cubic_RPA) THEN

               ALLOCATE (matrix_sigma_x_minus_vxc(ispin, ikp)%matrix)
               CALL dbcsr_create(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, &
                                 template=matrix_ks_kp_re(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_symmetric)

               CALL dbcsr_copy(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, matrix_ks_kp_re(ispin, ikp)%matrix)
               CALL dbcsr_set(matrix_ks_kp_re(ispin, ikp)%matrix, 0.0_dp)

               ALLOCATE (matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix)
               CALL dbcsr_create(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix, &
                                 template=matrix_ks_kp_im(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)

               CALL dbcsr_copy(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix, matrix_ks_kp_im(ispin, ikp)%matrix)
               CALL dbcsr_set(matrix_ks_kp_im(ispin, ikp)%matrix, 0.0_dp)

            ELSE

               ALLOCATE (matrix_sigma_x_minus_vxc(ispin, ikp)%matrix)
               CALL dbcsr_create(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, &
                                 template=matrix_ks(1)%matrix)

               CALL dbcsr_copy(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, matrix_ks(ispin)%matrix)
               CALL dbcsr_set(matrix_ks(ispin)%matrix, 0.0_dp)

            END IF

         ENDDO
      ENDDO

      ! set DFT functional to none and hfx_fraction to zero
      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx)

      IF (do_hfx) THEN
         hfx_fraction = qs_env%x_data(1, 1)%general_parameter%fraction
         qs_env%x_data(:, :)%general_parameter%fraction = 0.0_dp
      END IF
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=xc_none)

      ! in ADMM, also set the XC functional for ADMM correction to none
      ! do not do this if we do ADMM for Sigma_x
      IF (dft_control%do_admm) THEN
         xc_section_admm_aux => section_vals_get_subs_vals(admm_env%xc_section_aux, &
                                                           "XC_FUNCTIONAL")
         CALL section_vals_val_get(xc_section_admm_aux, "_SECTION_PARAMETERS_", &
                                   i_val=myfun_aux)
         CALL section_vals_val_set(xc_section_admm_aux, "_SECTION_PARAMETERS_", &
                                   i_val=xc_none)

         ! the same for the primary basis
         xc_section_admm_prim => section_vals_get_subs_vals(admm_env%xc_section_primary, &
                                                            "XC_FUNCTIONAL")
         CALL section_vals_val_get(xc_section_admm_aux, "_SECTION_PARAMETERS_", &
                                   i_val=myfun_prim)
         CALL section_vals_val_set(xc_section_admm_prim, "_SECTION_PARAMETERS_", &
                                   i_val=xc_none)

         ! for ADMMQ/S, set the charge_constrain to false (otherwise wrong results)
         charge_constrain_tmp = .FALSE.
         IF (admm_env%charge_constrain) THEN
            admm_env%charge_constrain = .FALSE.
            charge_constrain_tmp = .TRUE.
         END IF

      END IF

      ! if we do ADMM for Sigma_x, we write the ADMM correction into matrix_ks_aux_fit
      ! and therefore we should set it to zero
      IF (do_admm_rpa) THEN
         DO ispin = 1, nspins
            CALL dbcsr_set(matrix_ks_aux_fit(ispin)%matrix, 0.0_dp)
         END DO
      END IF

      IF (.NOT. mp2_env%ri_g0w0%update_xc_energy) THEN
         energy_total = energy%total
         energy_exc = energy%exc
         energy_exc1 = energy%exc1
         energy_exc_aux_fit = energy%ex
         energy_exc1_aux_fit = energy%exc_aux_fit
         energy_ex = energy%exc1_aux_fit
      END IF

      ! Remove the Exchange-correlation energy contributions from the total energy
      energy%total = energy%total - (energy%exc + energy%exc1 + energy%ex + &
                                     energy%exc_aux_fit + energy%exc1_aux_fit)

      ! calculate KS-matrix without XC and without HF
      CALL qs_ks_build_kohn_sham_matrix(qs_env=qs_env, calculate_forces=.FALSE., &
                                        just_energy=.FALSE.)

      IF (.NOT. mp2_env%ri_g0w0%update_xc_energy) THEN
         energy%exc = energy_exc
         energy%exc1 = energy_exc1
         energy%exc_aux_fit = energy_ex
         energy%exc1_aux_fit = energy_exc_aux_fit
         energy%ex = energy_exc1_aux_fit
         energy%total = energy_total
      END IF

      ! set the DFT functional and HF fraction back
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)
      IF (do_hfx) THEN
         qs_env%x_data(:, :)%general_parameter%fraction = hfx_fraction
      END IF

      IF (dft_control%do_admm) THEN
         xc_section_admm_aux => section_vals_get_subs_vals(admm_env%xc_section_aux, &
                                                           "XC_FUNCTIONAL")
         xc_section_admm_prim => section_vals_get_subs_vals(admm_env%xc_section_primary, &
                                                            "XC_FUNCTIONAL")

         CALL section_vals_val_set(xc_section_admm_aux, "_SECTION_PARAMETERS_", &
                                   i_val=myfun_aux)
         CALL section_vals_val_set(xc_section_admm_prim, "_SECTION_PARAMETERS_", &
                                   i_val=myfun_prim)
         IF (charge_constrain_tmp) THEN
            admm_env%charge_constrain = .TRUE.
         END IF
      END IF

      IF (do_kpoints_cubic_RPA) THEN
         CALL transform_matrix_ks_to_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)
      END IF

      ! remove the single-particle part (kin. En + Hartree pot) and change the sign
      DO ispin = 1, nspins
         IF (do_kpoints_cubic_RPA) THEN
            DO ikp = 1, nkp
               CALL dbcsr_add(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, matrix_ks_kp_re(ispin, ikp)%matrix, -1.0_dp, 1.0_dp)
               CALL dbcsr_add(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix, matrix_ks_kp_im(ispin, ikp)%matrix, -1.0_dp, 1.0_dp)
            END DO
         ELSE
            CALL dbcsr_add(matrix_sigma_x_minus_vxc(ispin, 1)%matrix, matrix_ks(ispin)%matrix, -1.0_dp, 1.0_dp)
         END IF
      END DO

      IF (do_kpoints_cubic_RPA) THEN

         CALL transform_sigma_x_minus_vxc_to_MO_basis(kpoints, matrix_sigma_x_minus_vxc, &
                                                      matrix_sigma_x_minus_vxc_im, &
                                                      vec_Sigma_x_minus_vxc_gw, &
                                                      vec_Sigma_x_minus_vxc_gw_im, &
                                                      para_env, nmo, mp2_env)

      ELSE

         DO ispin = 1, nspins
            CALL dbcsr_set(matrix_ks(ispin)%matrix, 0.0_dp)
            IF (do_admm_rpa) THEN
               CALL dbcsr_set(matrix_ks_aux_fit(ispin)%matrix, 0.0_dp)
            END IF
         END DO

         hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")

         CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)

         ! in most cases, we calculate the exchange self-energy here. But if we do only RI for
         ! the exchange self-energy, we do not calculate exchange here
         ehfx = 0.0_dp
         IF (.NOT. do_ri_Sigma_x) THEN

            ! add here HFX (=Sigma_exchange) to matrix_sigma_x_minus_vxc
            DO irep = 1, n_rep_hf
               ns = SIZE(rho_ao)
               rho_ao_2d(1:ns, 1:1) => rho_ao(1:ns)
               ns = SIZE(matrix_ks)
               IF (do_admm_rpa) THEN
                  matrix_ks_2d(1:ns, 1:1) => matrix_ks_aux_fit(1:ns)
               ELSE
                  matrix_ks_2d(1:ns, 1:1) => matrix_ks(1:ns)
               END IF

               CALL integrate_four_center(qs_env, qs_env%mp2_env%ri_rpa%x_data, matrix_ks_2d, eh1, &
                                          rho_ao_2d, hfx_sections, &
                                          para_env, .TRUE., irep, .TRUE., &
                                          ispin=1)
               ehfx = ehfx + eh1
            END DO
         END IF
         energy_ex = ehfx

         ! transform Fock-Matrix (calculated in integrate_four_center, written in matrix_ks_aux_fit in case
         ! of ADMM) from ADMM basis to primary basis
         IF (do_admm_rpa) THEN
            CALL admm_mo_merge_ks_matrix(qs_env)
         END IF

         DO ispin = 1, nspins
            CALL dbcsr_add(matrix_sigma_x_minus_vxc(ispin, 1)%matrix, matrix_ks(ispin)%matrix, 1.0_dp, 1.0_dp)
         END DO

         CALL dbcsr_desymmetrize(matrix_ks(1)%matrix, mo_coeff_b)
         CALL dbcsr_set(mo_coeff_b, 0.0_dp)

         ! Transform matrix_sigma_x_minus_vxc to MO basis
         DO ispin = 1, nspins

            CALL get_mo_set(mo_set=mos_mp2(ispin)%mo_set, &
                            mo_coeff=mo_coeff, &
                            nmo=nmo, &
                            homo=homo, &
                            nao=dimen)

            IF (ispin == 1) THEN

               ALLOCATE (vec_Sigma_x_minus_vxc_gw(nmo, nspins, nkp))
               vec_Sigma_x_minus_vxc_gw = 0.0_dp
            END IF

            CALL dbcsr_set(mo_coeff_b, 0.0_dp)
            CALL copy_fm_to_dbcsr(mo_coeff, mo_coeff_b, keep_sparsity=.FALSE.)

            ! initialize matrix_tmp and matrix_tmp2
            IF (ispin == 1) THEN
               CALL dbcsr_create(matrix_tmp, template=mo_coeff_b)
               CALL dbcsr_copy(matrix_tmp, mo_coeff_b)
               CALL dbcsr_set(matrix_tmp, 0.0_dp)

               CALL dbcsr_create(matrix_tmp_2, template=mo_coeff_b)
               CALL dbcsr_copy(matrix_tmp_2, mo_coeff_b)
               CALL dbcsr_set(matrix_tmp_2, 0.0_dp)
            END IF

            gw_corr_lev_occ = mp2_env%ri_g0w0%corr_mos_occ
            gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_virt
            ! if corrected occ/virt levels exceed the number of occ/virt levels,
            ! correct all occ/virt level energies
            IF (gw_corr_lev_occ > homo) gw_corr_lev_occ = homo
            IF (gw_corr_lev_virt > dimen - homo) gw_corr_lev_virt = dimen - homo
            IF (ispin == 1) THEN
               mp2_env%ri_g0w0%corr_mos_occ = gw_corr_lev_occ
               mp2_env%ri_g0w0%corr_mos_virt = gw_corr_lev_virt
            ELSE IF (ispin == 2) THEN
               ! ensure that the total number of corrected MOs is the same for alpha and beta, important
               ! for parallelization
               IF (mp2_env%ri_g0w0%corr_mos_occ + mp2_env%ri_g0w0%corr_mos_virt /= &
                   gw_corr_lev_occ + gw_corr_lev_virt) THEN
                  gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_occ + mp2_env%ri_g0w0%corr_mos_virt - gw_corr_lev_occ
               END IF
               mp2_env%ri_g0w0%corr_mos_occ_beta = gw_corr_lev_occ
               mp2_env%ri_g0w0%corr_mos_virt_beta = gw_corr_lev_virt

            END IF

            CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sigma_x_minus_vxc(ispin, 1)%matrix, &
                                mo_coeff_b, 0.0_dp, matrix_tmp, first_column=homo + 1 - gw_corr_lev_occ, &
                                last_column=homo + gw_corr_lev_virt)

            CALL dbcsr_multiply('T', 'N', 1.0_dp, mo_coeff_b, &
                                matrix_tmp, 0.0_dp, matrix_tmp_2, first_row=homo + 1 - gw_corr_lev_occ, &
                                last_row=homo + gw_corr_lev_virt)

            CALL dbcsr_get_diag(matrix_tmp_2, vec_Sigma_x_minus_vxc_gw(:, ispin, 1))

            CALL dbcsr_set(matrix_tmp, 0.0_dp)
            CALL dbcsr_set(matrix_tmp_2, 0.0_dp)

         END DO

         CALL mp_sum(vec_Sigma_x_minus_vxc_gw, para_env%group)

      END IF

      CALL dbcsr_release(mo_coeff_b)
      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_tmp_2)
      IF (do_kpoints_cubic_RPA) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_ks_kp_re)
         CALL dbcsr_deallocate_matrix_set(matrix_ks_kp_im)
      END IF

      DO ispin = 1, nspins
         DO ikp = 1, nkp
            CALL dbcsr_release_p(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix)
            IF (do_kpoints_cubic_RPA) THEN
               CALL dbcsr_release_p(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix)
            END IF
         END DO
      END DO

      ALLOCATE (mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(nmo, nspins, nkp))

      print_exx = mp2_env%ri_g0w0%print_exx

      IF (print_exx == gw_print_exx) THEN

         IF (unit_nr > 0) THEN

            WRITE (unit_nr, '(T3,A)') 'Exchange energies'
            WRITE (unit_nr, '(T3,A)') '-----------------'
            WRITE (unit_nr, '(T3,A)') ''
            WRITE (unit_nr, '(T6,2A)') 'MO                    Sigma_x-vxc'
            DO n_level_gw = 1, gw_corr_lev_occ + gw_corr_lev_virt

               n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
               IF (n_level_gw <= gw_corr_lev_occ) THEN
                  occ_virt = 'occ'
               ELSE
                  occ_virt = 'vir'
               END IF

               exx_minus_vxc = REAL(vec_Sigma_x_minus_vxc_gw(n_level_gw_ref, 1, 1)*evolt, kind=dp)

               WRITE (unit_nr, '(T4,I4,3A,3F21.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')  ', exx_minus_vxc

            END DO
            WRITE (unit_nr, '(T3,A)') 'End of exchange energies'
            WRITE (unit_nr, '(T3,A)') '------------------------'
            WRITE (unit_nr, '(T3,A)') ''
         END IF

      END IF

      IF (print_exx == gw_read_exx) THEN

         CALL open_file(unit_number=iunit, file_name="exx.out")

         really_read_line = .FALSE.

         DO WHILE (.TRUE.)

            READ (iunit, '(A)') line

            IF (line == "  End of exchange energies              ") EXIT

            IF (really_read_line) THEN

               READ (line(1:7), *) n_level_gw_ref
               READ (line(17:40), *) tmp

               DO ikp = 1, SIZE(vec_Sigma_x_minus_vxc_gw, 3)
                  vec_Sigma_x_minus_vxc_gw(n_level_gw_ref, 1, ikp) = tmp/evolt
               END DO

            END IF

            IF (line == "     MO                    Sigma_x-vxc  ") really_read_line = .TRUE.

         END DO

         CALL close_file(iunit)

      END IF

      ! store vec_Sigma_x_minus_vxc_gw in the mp2_environment
      mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, :, :) = vec_Sigma_x_minus_vxc_gw(:, :, :)

      ! clean up
      DEALLOCATE (matrix_sigma_x_minus_vxc, vec_Sigma_x_minus_vxc_gw)
      IF (do_kpoints_cubic_RPA) THEN
         DEALLOCATE (matrix_sigma_x_minus_vxc_im)
      END IF

      t2 = m_walltime()

      t3 = t2 - t1

      CALL timestop(handle)

   END SUBROUTINE compute_vec_Sigma_x_minus_vxc_gw

! **************************************************************************************************
!> \brief ...
!> \param kpoints ...
!> \param matrix_sigma_x_minus_vxc ...
!> \param matrix_sigma_x_minus_vxc_im ...
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param vec_Sigma_x_minus_vxc_gw_im ...
!> \param para_env ...
!> \param nmo ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE transform_sigma_x_minus_vxc_to_MO_basis(kpoints, matrix_sigma_x_minus_vxc, &
                                                      matrix_sigma_x_minus_vxc_im, vec_Sigma_x_minus_vxc_gw, &
                                                      vec_Sigma_x_minus_vxc_gw_im, para_env, nmo, mp2_env)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_sigma_x_minus_vxc, &
                                                            matrix_sigma_x_minus_vxc_im
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: vec_Sigma_x_minus_vxc_gw, &
                                                            vec_Sigma_x_minus_vxc_gw_im
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: nmo
      TYPE(mp2_type), POINTER                            :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'transform_sigma_x_minus_vxc_to_MO_basis', &
         routineP = moduleN//':'//routineN

      INTEGER :: dimen, gw_corr_lev_occ, gw_corr_lev_virt, handle, homo, i_global, iiB, ikp, &
         ispin, j_global, jjB, ncol_local, nkp, nrow_local, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: imval, reval
      TYPE(cp_cfm_type), POINTER                         :: cfm_mos, cfm_sigma_x_minus_vxc, &
                                                            cfm_sigma_x_minus_vxc_mo_basis, cfm_tmp
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: fwork_im, fwork_re
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mo_set_type), POINTER                         :: mo_set, mo_set_im, mo_set_re

      CALL timeset(routineN, handle)

      mo_set => kpoints%kp_env(1)%kpoint_env%mos(1, 1)%mo_set
      CALL get_mo_set(mo_set, nmo=nmo)

      nspins = SIZE(matrix_sigma_x_minus_vxc, 1)
      CALL get_kpoint_info(kpoints, nkp=nkp, kp_range=kp_range)

      ! if this CPASSERT is wrong, please make sure that the kpoint group size PARALLEL_GROUP_SIZE
      ! in the kpoint environment &DFT &KPOINTS is -1
      CPASSERT(kp_range(1) == 1 .AND. kp_range(2) == nkp)

      ALLOCATE (vec_Sigma_x_minus_vxc_gw(nmo, nspins, nkp))
      vec_Sigma_x_minus_vxc_gw = 0.0_dp

      ALLOCATE (vec_Sigma_x_minus_vxc_gw_im(nmo, nspins, nkp))
      vec_Sigma_x_minus_vxc_gw_im = 0.0_dp

      NULLIFY (fwork_re, fwork_im, cfm_mos, cfm_sigma_x_minus_vxc, &
               cfm_sigma_x_minus_vxc_mo_basis, cfm_tmp)
      CALL cp_fm_get_info(mo_set%mo_coeff, matrix_struct=matrix_struct)
      CALL cp_fm_create(fwork_re, matrix_struct)
      CALL cp_fm_create(fwork_im, matrix_struct)
      CALL cp_cfm_create(cfm_mos, matrix_struct)
      CALL cp_cfm_create(cfm_sigma_x_minus_vxc, matrix_struct)
      CALL cp_cfm_create(cfm_sigma_x_minus_vxc_mo_basis, matrix_struct)
      CALL cp_cfm_create(cfm_tmp, matrix_struct)

      CALL cp_cfm_get_info(matrix=cfm_sigma_x_minus_vxc_mo_basis, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      ! Transform matrix_sigma_x_minus_vxc to MO basis
      DO ikp = 1, nkp

         kp => kpoints%kp_env(ikp)%kpoint_env

         DO ispin = 1, nspins

            ! v_xc_n to fm matrix
            CALL copy_dbcsr_to_fm(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, fwork_re)
            CALL copy_dbcsr_to_fm(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix, fwork_im)

            CALL cp_cfm_scale_and_add_fm(z_zero, cfm_sigma_x_minus_vxc, z_one, fwork_re)
            CALL cp_cfm_scale_and_add_fm(z_one, cfm_sigma_x_minus_vxc, gaussi, fwork_im)

            ! get real part (1) and imag. part (2) of the mo coeffs
            mo_set_re => kp%mos(1, ispin)%mo_set
            mo_set_im => kp%mos(2, ispin)%mo_set

            CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mos, z_one, mo_set_re%mo_coeff)
            CALL cp_cfm_scale_and_add_fm(z_one, cfm_mos, gaussi, mo_set_im%mo_coeff)

            ! tmp = V(k)*C(k)
            CALL cp_cfm_gemm('N', 'N', nmo, nmo, nmo, z_one, cfm_sigma_x_minus_vxc, &
                             cfm_mos, z_zero, cfm_tmp)

            ! V_n(k) = C^H(k)*tmp
            CALL cp_cfm_gemm('C', 'N', nmo, nmo, nmo, z_one, cfm_mos, cfm_tmp, &
                             z_zero, cfm_sigma_x_minus_vxc_mo_basis)

            DO jjB = 1, ncol_local

               j_global = col_indices(jjB)

               DO iiB = 1, nrow_local

                  i_global = row_indices(iiB)

                  IF (j_global == i_global .AND. i_global <= nmo) THEN

                     reval = REAL(cfm_sigma_x_minus_vxc_mo_basis%local_data(iiB, jjB), kind=dp)
                     imval = AIMAG(cfm_sigma_x_minus_vxc_mo_basis%local_data(iiB, jjB))

                     vec_Sigma_x_minus_vxc_gw(i_global, ispin, ikp) = reval
                     vec_Sigma_x_minus_vxc_gw_im(i_global, ispin, ikp) = imval

                  END IF

               END DO

            END DO

         END DO

      END DO

      CALL mp_sum(vec_Sigma_x_minus_vxc_gw, para_env%group)
      CALL mp_sum(vec_Sigma_x_minus_vxc_gw_im, para_env%group)

      ! also adjust in the case of kpoints too big gw_corr_lev_occ and gw_corr_lev_virt
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=kpoints%kp_env(1)%kpoint_env%mos(ispin, 1)%mo_set, &
                         homo=homo, nao=dimen)
         gw_corr_lev_occ = mp2_env%ri_g0w0%corr_mos_occ
         gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_virt
         ! if corrected occ/virt levels exceed the number of occ/virt levels,
         ! correct all occ/virt level energies
         IF (gw_corr_lev_occ > homo) gw_corr_lev_occ = homo
         IF (gw_corr_lev_virt > dimen - homo) gw_corr_lev_virt = dimen - homo
         IF (ispin == 1) THEN
            mp2_env%ri_g0w0%corr_mos_occ = gw_corr_lev_occ
            mp2_env%ri_g0w0%corr_mos_virt = gw_corr_lev_virt
         ELSE IF (ispin == 2) THEN
            ! ensure that the total number of corrected MOs is the same for alpha and beta, important
            ! for parallelization
            IF (mp2_env%ri_g0w0%corr_mos_occ + mp2_env%ri_g0w0%corr_mos_virt /= &
                gw_corr_lev_occ + gw_corr_lev_virt) THEN
               gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_occ + mp2_env%ri_g0w0%corr_mos_virt - gw_corr_lev_occ
            END IF
            mp2_env%ri_g0w0%corr_mos_occ_beta = gw_corr_lev_occ
            mp2_env%ri_g0w0%corr_mos_virt_beta = gw_corr_lev_virt
         END IF
      END DO

      CALL cp_fm_release(fwork_re)
      CALL cp_fm_release(fwork_im)
      CALL cp_cfm_release(cfm_mos)
      CALL cp_cfm_release(cfm_sigma_x_minus_vxc)
      CALL cp_cfm_release(cfm_sigma_x_minus_vxc_mo_basis)
      CALL cp_cfm_release(cfm_tmp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_ks_transl ...
!> \param matrix_ks_kp_re ...
!> \param matrix_ks_kp_im ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE transform_matrix_ks_to_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_transl, matrix_ks_kp_re, &
                                                            matrix_ks_kp_im
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_matrix_ks_to_kp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikp, ispin, nkp, nspin
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      NULLIFY (sab_nl)
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, sab_nl=sab_nl, cell_to_index=cell_to_index)

      CPASSERT(ASSOCIATED(sab_nl))

      nspin = SIZE(matrix_ks_transl, 1)

      DO ikp = 1, nkp
         DO ispin = 1, nspin

            CALL dbcsr_set(matrix_ks_kp_re(ispin, ikp)%matrix, 0.0_dp)
            CALL dbcsr_set(matrix_ks_kp_im(ispin, ikp)%matrix, 0.0_dp)
            CALL rskp_transform(rmatrix=matrix_ks_kp_re(ispin, ikp)%matrix, &
                                cmatrix=matrix_ks_kp_im(ispin, ikp)%matrix, &
                                rsmat=matrix_ks_transl, ispin=ispin, &
                                xkp=xkp(1:3, ikp), cell_to_index=cell_to_index, sab_nl=sab_nl)

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_ks_transl ...
!> \param matrix_ks_kp_re ...
!> \param matrix_ks_kp_im ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE allocate_matrix_ks_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_transl, matrix_ks_kp_re, &
                                                            matrix_ks_kp_im
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_matrix_ks_kp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikp, ispin, nkp, nspin
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      NULLIFY (sab_nl)
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, sab_nl=sab_nl, cell_to_index=cell_to_index)

      CPASSERT(ASSOCIATED(sab_nl))

      nspin = SIZE(matrix_ks_transl, 1)

      NULLIFY (matrix_ks_kp_re, matrix_ks_kp_im)
      CALL dbcsr_allocate_matrix_set(matrix_ks_kp_re, nspin, nkp)
      CALL dbcsr_allocate_matrix_set(matrix_ks_kp_im, nspin, nkp)

      DO ikp = 1, nkp
      DO ispin = 1, nspin

         ALLOCATE (matrix_ks_kp_re(ispin, ikp)%matrix)
         ALLOCATE (matrix_ks_kp_im(ispin, ikp)%matrix)

         CALL dbcsr_create(matrix_ks_kp_re(ispin, ikp)%matrix, &
                           template=matrix_ks_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_create(matrix_ks_kp_im(ispin, ikp)%matrix, &
                           template=matrix_ks_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)

         CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks_kp_re(ispin, ikp)%matrix, sab_nl)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks_kp_im(ispin, ikp)%matrix, sab_nl)

         CALL dbcsr_set(matrix_ks_kp_re(ispin, ikp)%matrix, 0.0_dp)
         CALL dbcsr_set(matrix_ks_kp_im(ispin, ikp)%matrix, 0.0_dp)

      END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

END MODULE rpa_gw_sigma

