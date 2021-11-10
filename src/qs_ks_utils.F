!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines that build the Kohn-Sham matrix (i.e calculate the coulomb
!>      and xc parts
!> \par History
!>      05.2002 moved from qs_scf (see there the history) [fawzi]
!>      JGH [30.08.02] multi-grid arrays independent from density and potential
!>      10.2002 introduced pools, uses updated rho as input,
!>              removed most temporary variables, renamed may vars,
!>              began conversion to LSD [fawzi]
!>      10.2004 moved calculate_w_matrix here [Joost VandeVondele]
!>              introduced energy derivative wrt MOs [Joost VandeVondele]
!> \author Fawzi Mohamed
! **************************************************************************************************

MODULE qs_ks_utils
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_ddapc,                        ONLY: cp_ddapc_apply_CD
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_deallocate_matrix, dbcsr_dot, dbcsr_get_info, dbcsr_init_p, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release_p, dbcsr_scale, dbcsr_scale_by_vector, &
        dbcsr_set, dbcsr_type
   USE input_constants,                 ONLY: &
        cdft_alpha_constraint, cdft_beta_constraint, cdft_charge_constraint, &
        cdft_magnetization_constraint, do_admm_aux_exch_func_none, do_ppl_grid, sic_ad, sic_eo, &
        sic_list_all, sic_list_unpaired, sic_mauri_spz, sic_mauri_us, sic_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kahan_sum,                       ONLY: accurate_dot_product,&
                                              accurate_sum
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE lri_environment_methods,         ONLY: v_int_ppl_update
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type,&
                                              lri_kind_type
   USE lri_forces,                      ONLY: calculate_lri_forces,&
                                              calculate_ri_forces
   USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix,&
                                              calculate_ri_ks_matrix
   USE message_passing,                 ONLY: mp_sum
   USE ps_implicit_types,               ONLY: MIXED_BC,&
                                              MIXED_PERIODIC_BC,&
                                              NEUMANN_BC,&
                                              PERIODIC_BC
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: PW_MODE_DISTRIBUTED
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_implicit,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_cdft_types,                   ONLY: cdft_control_type
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace,&
                                              integrate_v_rspace_diagonal,&
                                              integrate_v_rspace_one_center
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_qmmm_methods,              ONLY: qmmm_modify_hartree_pot
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE task_list_types,                 ONLY: task_list_type
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: xc_exc_calc,&
                                              xc_vxc_pw_create1
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_utils'

   PUBLIC :: low_spin_roks, sic_explicit_orbitals, calc_v_sic_rspace, print_densities, &
             print_detailed_energy, compute_matrix_vxc, sum_up_and_integrate, &
             calculate_zmp_potential, get_embed_potential_energy

CONTAINS

! **************************************************************************************************
!> \brief do ROKS calculations yielding low spin states
!> \param energy ...
!> \param qs_env ...
!> \param dft_control ...
!> \param just_energy ...
!> \param calculate_forces ...
!> \param auxbas_pw_pool ...
! **************************************************************************************************
   SUBROUTINE low_spin_roks(energy, qs_env, dft_control, just_energy, &
                            calculate_forces, auxbas_pw_pool)

      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL, INTENT(IN)                                :: just_energy, calculate_forces
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CHARACTER(*), PARAMETER                            :: routineN = 'low_spin_roks'

      INTEGER                                            :: handle, ispin, iterm, k, k_alpha, &
                                                            k_beta, n_rep, Nelectron, Nspin, Nterms
      INTEGER, DIMENSION(:), POINTER                     :: ivec
      INTEGER, DIMENSION(:, :, :), POINTER               :: occupations
      LOGICAL                                            :: compute_virial, in_range, &
                                                            uniform_occupation
      REAL(KIND=dp)                                      :: exc, total_rho
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc_tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: energy_scaling, rvec, scaling
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_p, mo_derivs, rho_ao
      TYPE(dbcsr_type), POINTER                          :: dbcsr_deriv, fm_deriv, fm_scaled, &
                                                            mo_coeff
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: work_v_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, tau, vxc, vxc_tau
      TYPE(pw_pool_type), POINTER                        :: xc_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, low_spin_roks_section, xc_section
      TYPE(virial_type), POINTER                         :: virial

      IF (.NOT. dft_control%low_spin_roks) RETURN
      NULLIFY (ks_env, rho_ao)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      mo_derivs=mo_derivs, &
                      mos=mo_array, &
                      rho=rho, &
                      pw_env=pw_env, &
                      input=input, &
                      cell=cell, &
                      virial=virial)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      compute_virial = virial%pv_calculate .AND. (.NOT. virial%pv_numer)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")

      ! some assumptions need to be checked
      ! we have two spins
      CPASSERT(SIZE(mo_array, 1) == 2)
      Nspin = 2
      ! we want uniform occupations
      CALL get_mo_set(mo_set=mo_array(1)%mo_set, uniform_occupation=uniform_occupation)
      CPASSERT(uniform_occupation)
      CALL get_mo_set(mo_set=mo_array(2)%mo_set, mo_coeff_b=mo_coeff, uniform_occupation=uniform_occupation)
      CPASSERT(uniform_occupation)

      NULLIFY (dbcsr_deriv)
      CALL dbcsr_init_p(dbcsr_deriv)
      CALL dbcsr_copy(dbcsr_deriv, mo_derivs(1)%matrix)
      CALL dbcsr_set(dbcsr_deriv, 0.0_dp)

      ! basic info
      CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff_b=mo_coeff)
      CALL dbcsr_get_info(mo_coeff, nfullcols_total=k_alpha)
      CALL get_mo_set(mo_set=mo_array(2)%mo_set, mo_coeff_b=mo_coeff)
      CALL dbcsr_get_info(mo_coeff, nfullcols_total=k_beta)

      ! read the input
      low_spin_roks_section => section_vals_get_subs_vals(input, "DFT%LOW_SPIN_ROKS")

      CALL section_vals_val_get(low_spin_roks_section, "ENERGY_SCALING", r_vals=rvec)
      Nterms = SIZE(rvec)
      ALLOCATE (energy_scaling(Nterms))
      energy_scaling = rvec !? just wondering, should this add up to 1, in which case we should cpp?

      CALL section_vals_val_get(low_spin_roks_section, "SPIN_CONFIGURATION", n_rep_val=n_rep)
      CPASSERT(n_rep == Nterms)
      CALL section_vals_val_get(low_spin_roks_section, "SPIN_CONFIGURATION", i_rep_val=1, i_vals=ivec)
      Nelectron = SIZE(ivec)
      CPASSERT(Nelectron == k_alpha - k_beta)
      ALLOCATE (occupations(2, Nelectron, Nterms))
      occupations = 0
      DO iterm = 1, Nterms
         CALL section_vals_val_get(low_spin_roks_section, "SPIN_CONFIGURATION", i_rep_val=iterm, i_vals=ivec)
         CPASSERT(Nelectron == SIZE(ivec))
         in_range = ALL(ivec >= 1) .AND. ALL(ivec <= 2)
         CPASSERT(in_range)
         DO k = 1, Nelectron
            occupations(ivec(k), k, iterm) = 1
         END DO
      END DO

      ! set up general data structures
      ! density matrices, kohn-sham matrices

      NULLIFY (matrix_p)
      CALL dbcsr_allocate_matrix_set(matrix_p, Nspin)
      DO ispin = 1, Nspin
         ALLOCATE (matrix_p(ispin)%matrix)
         CALL dbcsr_copy(matrix_p(ispin)%matrix, rho_ao(1)%matrix, &
                         name="density matrix low spin roks")
         CALL dbcsr_set(matrix_p(ispin)%matrix, 0.0_dp)
      END DO

      NULLIFY (matrix_h)
      CALL dbcsr_allocate_matrix_set(matrix_h, Nspin)
      DO ispin = 1, Nspin
         ALLOCATE (matrix_h(ispin)%matrix)
         CALL dbcsr_copy(matrix_h(ispin)%matrix, rho_ao(1)%matrix, &
                         name="KS matrix low spin roks")
         CALL dbcsr_set(matrix_h(ispin)%matrix, 0.0_dp)
      END DO

      ! grids in real and g space for rho and vxc
      ! tau functionals are not supported
      NULLIFY (tau, vxc_tau, vxc)
      CALL pw_env_get(pw_env, xc_pw_pool=xc_pw_pool)

      ALLOCATE (rho_r(Nspin))
      ALLOCATE (rho_g(Nspin))
      DO ispin = 1, Nspin
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_r(ispin)%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
      END DO
      CALL pw_pool_create_pw(auxbas_pw_pool, work_v_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ! get mo matrices needed to construct the density matrices
      ! we will base all on the alpha spin matrix, obviously possible in ROKS
      CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff_b=mo_coeff)
      NULLIFY (fm_scaled, fm_deriv)
      CALL dbcsr_init_p(fm_scaled)
      CALL dbcsr_init_p(fm_deriv)
      CALL dbcsr_copy(fm_scaled, mo_coeff)
      CALL dbcsr_copy(fm_deriv, mo_coeff)

      ALLOCATE (scaling(k_alpha))

      ! for each term, add it with the given scaling factor to the energy, and compute the required derivatives
      DO iterm = 1, Nterms

         DO ispin = 1, Nspin
            ! compute the proper density matrices with the required occupations
            CALL dbcsr_set(matrix_p(ispin)%matrix, 0.0_dp)
            scaling = 1.0_dp
            scaling(k_alpha - Nelectron + 1:k_alpha) = occupations(ispin, :, iterm)
            CALL dbcsr_copy(fm_scaled, mo_coeff)
            CALL dbcsr_scale_by_vector(fm_scaled, scaling, side='right')
            CALL dbcsr_multiply('n', 't', 1.0_dp, mo_coeff, fm_scaled, &
                                0.0_dp, matrix_p(ispin)%matrix, retain_sparsity=.TRUE.)
            ! compute the densities on the grid
            CALL calculate_rho_elec(matrix_p=matrix_p(ispin)%matrix, &
                                    rho=rho_r(ispin), rho_gspace=rho_g(ispin), total_rho=total_rho, &
                                    ks_env=ks_env)
         END DO

         ! compute the exchange energies / potential if needed
         IF (just_energy) THEN
            exc = xc_exc_calc(rho_r=rho_r, rho_g=rho_g, tau=tau, xc_section=xc_section, &
                              pw_pool=xc_pw_pool)
         ELSE
            CPASSERT(.NOT. compute_virial)
            CALL xc_vxc_pw_create1(vxc_rho=vxc, rho_r=rho_r, &
                                   rho_g=rho_g, tau=tau, vxc_tau=vxc_tau, exc=exc, xc_section=xc_section, &
                                   pw_pool=xc_pw_pool, compute_virial=.FALSE., virial_xc=virial_xc_tmp)
         END IF

         energy%exc = energy%exc + energy_scaling(iterm)*exc

         ! add the corresponding derivatives to the MO derivatives
         IF (.NOT. just_energy) THEN
            ! get the potential in matrix form
            DO ispin = 1, Nspin
               ! use a work_v_rspace
               work_v_rspace%pw%cr3d = (energy_scaling(iterm)*vxc(ispin)%pw%pw_grid%dvol)* &
                                       vxc(ispin)%pw%cr3d
               ! zero first ?!
               CALL dbcsr_set(matrix_h(ispin)%matrix, 0.0_dp)
               CALL integrate_v_rspace(v_rspace=work_v_rspace, pmat=matrix_p(ispin), hmat=matrix_h(ispin), &
                                       qs_env=qs_env, calculate_forces=calculate_forces)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc(ispin)%pw)
            END DO
            DEALLOCATE (vxc)

            ! add this to the mo_derivs, again based on the alpha mo_coeff
            DO ispin = 1, Nspin
               CALL dbcsr_multiply('n', 'n', 1.0_dp, matrix_h(ispin)%matrix, mo_coeff, &
                                   0.0_dp, dbcsr_deriv, last_column=k_alpha)

               scaling = 1.0_dp
               scaling(k_alpha - Nelectron + 1:k_alpha) = occupations(ispin, :, iterm)
               CALL dbcsr_scale_by_vector(dbcsr_deriv, scaling, side='right')
               CALL dbcsr_add(mo_derivs(1)%matrix, dbcsr_deriv, 1.0_dp, 1.0_dp)
            END DO

         END IF

      END DO

      ! release allocated memory
      DO ispin = 1, Nspin
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g(ispin)%pw)
      END DO
      DEALLOCATE (rho_r, rho_g)
      CALL dbcsr_deallocate_matrix_set(matrix_p)
      CALL dbcsr_deallocate_matrix_set(matrix_h)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, work_v_rspace%pw)

      CALL dbcsr_release_p(fm_deriv)
      CALL dbcsr_release_p(fm_scaled)

      DEALLOCATE (occupations)
      DEALLOCATE (energy_scaling)
      DEALLOCATE (scaling)

      CALL dbcsr_release_p(dbcsr_deriv)

      CALL timestop(handle)

   END SUBROUTINE low_spin_roks
! **************************************************************************************************
!> \brief do sic calculations on explicit orbitals
!> \param energy ...
!> \param qs_env ...
!> \param dft_control ...
!> \param poisson_env ...
!> \param just_energy ...
!> \param calculate_forces ...
!> \param auxbas_pw_pool ...
! **************************************************************************************************
   SUBROUTINE sic_explicit_orbitals(energy, qs_env, dft_control, poisson_env, just_energy, &
                                    calculate_forces, auxbas_pw_pool)

      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      LOGICAL, INTENT(IN)                                :: just_energy, calculate_forces
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CHARACTER(*), PARAMETER :: routineN = 'sic_explicit_orbitals'

      INTEGER                                            :: handle, i, Iorb, k_alpha, k_beta, Norb
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: sic_orbital_list
      LOGICAL                                            :: compute_virial, uniform_occupation
      REAL(KIND=dp)                                      :: ener, exc, total_rho
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc_tmp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mo_derivs_local
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: matrix_hv, matrix_v, mo_coeff
      TYPE(dbcsr_p_type)                                 :: orb_density_matrix_p, orb_h_p
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs, rho_ao, tmp_dbcsr
      TYPE(dbcsr_type), POINTER                          :: orb_density_matrix, orb_h
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: orb_rho_g, orb_rho_r, tmp_g, tmp_r, &
                                                            work_v_gspace, work_v_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, tau, vxc, vxc_tau
      TYPE(pw_pool_type), POINTER                        :: xc_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, xc_section
      TYPE(virial_type), POINTER                         :: virial

      IF (dft_control%sic_method_id .NE. sic_eo) RETURN

      CALL timeset(routineN, handle)

      NULLIFY (tau, vxc_tau, mo_derivs, ks_env, rho_ao)

      ! generate the lists of orbitals that need sic treatment
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      mo_derivs=mo_derivs, &
                      mos=mo_array, &
                      rho=rho, &
                      pw_env=pw_env, &
                      input=input, &
                      cell=cell, &
                      virial=virial)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      compute_virial = virial%pv_calculate .AND. (.NOT. virial%pv_numer)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")

      DO i = 1, SIZE(mo_array) !fm->dbcsr
         IF (mo_array(i)%mo_set%use_mo_coeff_b) THEN !fm->dbcsr
            CALL copy_dbcsr_to_fm(mo_array(i)%mo_set%mo_coeff_b, &
                                  mo_array(i)%mo_set%mo_coeff) !fm->dbcsr
         END IF !fm->dbcsr
      END DO !fm->dbcsr

      CALL pw_env_get(pw_env, xc_pw_pool=xc_pw_pool)

      ! we have two spins
      CPASSERT(SIZE(mo_array, 1) == 2)
      ! we want uniform occupations
      CALL get_mo_set(mo_set=mo_array(1)%mo_set, uniform_occupation=uniform_occupation)
      CPASSERT(uniform_occupation)
      CALL get_mo_set(mo_set=mo_array(2)%mo_set, mo_coeff=mo_coeff, uniform_occupation=uniform_occupation)
      CPASSERT(uniform_occupation)

      NULLIFY (tmp_dbcsr)
      CALL dbcsr_allocate_matrix_set(tmp_dbcsr, SIZE(mo_derivs, 1))
      DO i = 1, SIZE(mo_derivs, 1) !fm->dbcsr
         !
         NULLIFY (tmp_dbcsr(i)%matrix)
         CALL dbcsr_init_p(tmp_dbcsr(i)%matrix)
         CALL dbcsr_copy(tmp_dbcsr(i)%matrix, mo_derivs(i)%matrix)
         CALL dbcsr_set(tmp_dbcsr(i)%matrix, 0.0_dp)
      END DO !fm->dbcsr

      k_alpha = 0; k_beta = 0
      SELECT CASE (dft_control%sic_list_id)
      CASE (sic_list_all)

         CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, ncol_global=k_alpha)

         IF (SIZE(mo_array, 1) > 1) THEN
            CALL get_mo_set(mo_set=mo_array(2)%mo_set, mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, ncol_global=k_beta)
         END IF

         Norb = k_alpha + k_beta
         ALLOCATE (sic_orbital_list(3, Norb))

         iorb = 0
         DO i = 1, k_alpha
            iorb = iorb + 1
            sic_orbital_list(1, iorb) = 1
            sic_orbital_list(2, iorb) = i
            sic_orbital_list(3, iorb) = 1
         END DO
         DO i = 1, k_beta
            iorb = iorb + 1
            sic_orbital_list(1, iorb) = 2
            sic_orbital_list(2, iorb) = i
            IF (SIZE(mo_derivs, 1) == 1) THEN
               sic_orbital_list(3, iorb) = 1
            ELSE
               sic_orbital_list(3, iorb) = 2
            END IF
         END DO

      CASE (sic_list_unpaired)
         ! we have two spins
         CPASSERT(SIZE(mo_array, 1) == 2)
         ! we have them restricted
         CPASSERT(SIZE(mo_derivs, 1) == 1)
         CPASSERT(dft_control%restricted)

         CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, ncol_global=k_alpha)

         CALL get_mo_set(mo_set=mo_array(2)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, ncol_global=k_beta)

         Norb = k_alpha - k_beta
         ALLOCATE (sic_orbital_list(3, Norb))

         iorb = 0
         DO i = k_beta + 1, k_alpha
            iorb = iorb + 1
            sic_orbital_list(1, iorb) = 1
            sic_orbital_list(2, iorb) = i
            ! we are guaranteed to be restricted
            sic_orbital_list(3, iorb) = 1
         END DO

      CASE DEFAULT
         CPABORT("")
      END SELECT

      ! data needed for each of the orbs
      CALL pw_pool_create_pw(auxbas_pw_pool, orb_rho_r%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, tmp_r%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, orb_rho_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, tmp_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, work_v_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, work_v_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ALLOCATE (orb_density_matrix)
      CALL dbcsr_copy(orb_density_matrix, rho_ao(1)%matrix, &
                      name="orb_density_matrix")
      CALL dbcsr_set(orb_density_matrix, 0.0_dp)
      orb_density_matrix_p%matrix => orb_density_matrix

      ALLOCATE (orb_h)
      CALL dbcsr_copy(orb_h, rho_ao(1)%matrix, &
                      name="orb_density_matrix")
      CALL dbcsr_set(orb_h, 0.0_dp)
      orb_h_p%matrix => orb_h

      CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff=mo_coeff)

      CALL cp_fm_struct_create(fm_struct_tmp, ncol_global=1, &
                               template_fmstruct=mo_coeff%matrix_struct)
      CALL cp_fm_create(matrix_v, fm_struct_tmp, name="matrix_v")
      CALL cp_fm_create(matrix_hv, fm_struct_tmp, name="matrix_hv")
      CALL cp_fm_struct_release(fm_struct_tmp)

      ALLOCATE (mo_derivs_local(SIZE(mo_array, 1)))
      DO I = 1, SIZE(mo_array, 1)
         CALL get_mo_set(mo_set=mo_array(i)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_create(mo_derivs_local(I)%matrix, mo_coeff%matrix_struct)
      END DO

      ALLOCATE (rho_r(2))
      rho_r(1)%pw => orb_rho_r%pw
      rho_r(2)%pw => tmp_r%pw
      CALL pw_zero(tmp_r%pw)

      ALLOCATE (rho_g(2))
      rho_g(1)%pw => orb_rho_g%pw
      rho_g(2)%pw => tmp_g%pw
      CALL pw_zero(tmp_g%pw)

      NULLIFY (vxc)
      ! ALLOCATE(vxc(2),stat=stat)
      ! CPPostcondition(stat==0,cp_failure_level,routineP,failure)
      ! CALL pw_pool_create_pw(xc_pw_pool,vxc(1)%pw,&
      !         in_space=REALSPACE, use_data=REALDATA3D)
      ! CALL pw_pool_create_pw(xc_pw_pool,vxc(2)%pw,&
      !         in_space=REALSPACE, use_data=REALDATA3D)

      ! now apply to SIC correction to each selected orbital
      DO iorb = 1, Norb
         ! extract the proper orbital from the mo_coeff
         CALL get_mo_set(mo_set=mo_array(sic_orbital_list(1, iorb))%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_to_fm(mo_coeff, matrix_v, 1, sic_orbital_list(2, iorb), 1)

         ! construct the density matrix and the corresponding density
         CALL dbcsr_set(orb_density_matrix, 0.0_dp)
         CALL cp_dbcsr_plus_fm_fm_t(orb_density_matrix, matrix_v=matrix_v, ncol=1, &
                                    alpha=1.0_dp)

         CALL calculate_rho_elec(matrix_p=orb_density_matrix, &
                                 rho=orb_rho_r, rho_gspace=orb_rho_g, total_rho=total_rho, &
                                 ks_env=ks_env)

         ! write(*,*) 'Orbital ',sic_orbital_list(1,iorb),sic_orbital_list(2,iorb)
         ! write(*,*) 'Total orbital rho= ',total_rho

         ! compute the energy functional for this orbital and its derivative

         CALL pw_poisson_solve(poisson_env, orb_rho_g%pw, ener, work_v_gspace%pw)
         ! no PBC correction is done here, see "calc_v_sic_rspace" for SIC methods
         ! with PBC aware corrections
         energy%hartree = energy%hartree - dft_control%sic_scaling_a*ener
         IF (.NOT. just_energy) THEN
            CALL pw_transfer(work_v_gspace%pw, work_v_rspace%pw)
            CALL pw_scale(work_v_rspace%pw, -dft_control%sic_scaling_a*work_v_rspace%pw%pw_grid%dvol)
            CALL dbcsr_set(orb_h, 0.0_dp)
         END IF

         IF (just_energy) THEN
            exc = xc_exc_calc(rho_r=rho_r, rho_g=rho_g, tau=tau, xc_section=xc_section, &
                              pw_pool=xc_pw_pool)
         ELSE
            CPASSERT(.NOT. compute_virial)
            CALL xc_vxc_pw_create1(vxc_rho=vxc, rho_r=rho_r, &
                                   rho_g=rho_g, tau=tau, vxc_tau=vxc_tau, exc=exc, xc_section=xc_section, &
                                   pw_pool=xc_pw_pool, compute_virial=compute_virial, virial_xc=virial_xc_tmp)
            ! add to the existing work_v_rspace
            work_v_rspace%pw%cr3d = work_v_rspace%pw%cr3d - &
                                    dft_control%sic_scaling_b*vxc(1)%pw%pw_grid%dvol*vxc(1)%pw%cr3d
         END IF
         energy%exc = energy%exc - dft_control%sic_scaling_b*exc

         IF (.NOT. just_energy) THEN
            ! note, orb_h (which is being pointed to with orb_h_p) is zeroed above
            CALL integrate_v_rspace(v_rspace=work_v_rspace, pmat=orb_density_matrix_p, hmat=orb_h_p, &
                                    qs_env=qs_env, calculate_forces=calculate_forces)

            ! add this to the mo_derivs
            CALL cp_dbcsr_sm_fm_multiply(orb_h, matrix_v, matrix_hv, 1)
            ! silly trick, copy to an array of the right size and add to mo_derivs
            CALL cp_fm_set_all(mo_derivs_local(sic_orbital_list(3, iorb))%matrix, 0.0_dp)
            CALL cp_fm_to_fm(matrix_hv, mo_derivs_local(sic_orbital_list(3, iorb))%matrix, 1, 1, sic_orbital_list(2, iorb))
            CALL copy_fm_to_dbcsr(mo_derivs_local(sic_orbital_list(3, iorb))%matrix, &
                                  tmp_dbcsr(sic_orbital_list(3, iorb))%matrix)
            CALL dbcsr_add(mo_derivs(sic_orbital_list(3, iorb))%matrix, &
                           tmp_dbcsr(sic_orbital_list(3, iorb))%matrix, 1.0_dp, 1.0_dp)
            !
            ! need to deallocate vxc
            CALL pw_pool_give_back_pw(xc_pw_pool, vxc(1)%pw)
            CALL pw_pool_give_back_pw(xc_pw_pool, vxc(2)%pw)
            DEALLOCATE (vxc)

         END IF

      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, orb_rho_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, tmp_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, orb_rho_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, tmp_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, work_v_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, work_v_rspace%pw)

      CALL dbcsr_deallocate_matrix(orb_density_matrix)
      CALL dbcsr_deallocate_matrix(orb_h)
      CALL cp_fm_release(matrix_v)
      CALL cp_fm_release(matrix_hv)
      DO I = 1, SIZE(mo_derivs_local, 1)
         CALL cp_fm_release(mo_derivs_local(I)%matrix)
      END DO
      DEALLOCATE (mo_derivs_local)
      DEALLOCATE (rho_r)
      DEALLOCATE (rho_g)

      CALL dbcsr_deallocate_matrix_set(tmp_dbcsr) !fm->dbcsr

      CALL timestop(handle)

   END SUBROUTINE sic_explicit_orbitals

! **************************************************************************************************
!> \brief do sic calculations on the spin density
!> \param v_sic_rspace ...
!> \param energy ...
!> \param qs_env ...
!> \param dft_control ...
!> \param rho ...
!> \param poisson_env ...
!> \param just_energy ...
!> \param calculate_forces ...
!> \param auxbas_pw_pool ...
! **************************************************************************************************
   SUBROUTINE calc_v_sic_rspace(v_sic_rspace, energy, &
                                qs_env, dft_control, rho, poisson_env, just_energy, &
                                calculate_forces, auxbas_pw_pool)

      TYPE(pw_p_type), INTENT(INOUT)                     :: v_sic_rspace
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      LOGICAL, INTENT(IN)                                :: just_energy, calculate_forces
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      INTEGER                                            :: i, nelec, nelec_a, nelec_b, nforce
      REAL(kind=dp)                                      :: ener, full_scaling, scaling
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: store_forces
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(pw_p_type)                                    :: work_rho, work_v
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      NULLIFY (mo_array, rho_g)

      IF (dft_control%sic_method_id == sic_none) RETURN
      IF (dft_control%sic_method_id == sic_eo) RETURN

      IF (dft_control%qs_control%gapw) &
         CPABORT("sic and GAPW not yet compatible")

      ! OK, right now we like two spins to do sic, could be relaxed for AD
      CPASSERT(dft_control%nspins == 2)

      CALL pw_pool_create_pw(auxbas_pw_pool, work_rho%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, work_v%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL qs_rho_get(rho, rho_g=rho_g)

      ! Hartree sic corrections
      SELECT CASE (dft_control%sic_method_id)
      CASE (sic_mauri_us, sic_mauri_spz)
         CALL pw_copy(rho_g(1)%pw, work_rho%pw)
         CALL pw_axpy(rho_g(2)%pw, work_rho%pw, alpha=-1._dp)
         CALL pw_poisson_solve(poisson_env, work_rho%pw, ener, work_v%pw)
      CASE (sic_ad)
         ! find out how many elecs we have
         CALL get_qs_env(qs_env, mos=mo_array)
         CALL get_mo_set(mo_set=mo_array(1)%mo_set, nelectron=nelec_a)
         CALL get_mo_set(mo_set=mo_array(2)%mo_set, nelectron=nelec_b)
         nelec = nelec_a + nelec_b
         CALL pw_copy(rho_g(1)%pw, work_rho%pw)
         CALL pw_axpy(rho_g(2)%pw, work_rho%pw)
         scaling = 1.0_dp/REAL(nelec, KIND=dp)
         CALL pw_scale(work_rho%pw, scaling)
         CALL pw_poisson_solve(poisson_env, work_rho%pw, ener, work_v%pw)
      CASE DEFAULT
         CPABORT("Unknown sic method id")
      END SELECT

      ! Correct for  DDAP charges (if any)
      ! storing whatever force might be there from previous decoupling
      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, force=force)
         nforce = 0
         DO i = 1, SIZE(force)
            nforce = nforce + SIZE(force(i)%ch_pulay, 2)
         END DO
         ALLOCATE (store_forces(3, nforce))
         nforce = 0
         DO i = 1, SIZE(force)
            store_forces(1:3, nforce + 1:nforce + SIZE(force(i)%ch_pulay, 2)) = force(i)%ch_pulay(:, :)
            force(i)%ch_pulay(:, :) = 0.0_dp
            nforce = nforce + SIZE(force(i)%ch_pulay, 2)
         END DO
      END IF

      CALL cp_ddapc_apply_CD(qs_env, &
                             work_rho, &
                             ener, &
                             v_hartree_gspace=work_v, &
                             calculate_forces=calculate_forces, &
                             Itype_of_density="SPIN")

      SELECT CASE (dft_control%sic_method_id)
      CASE (sic_mauri_us, sic_mauri_spz)
         full_scaling = -dft_control%sic_scaling_a
      CASE (sic_ad)
         full_scaling = -dft_control%sic_scaling_a*nelec
      CASE DEFAULT
         CPABORT("Unknown sic method id")
      END SELECT
      energy%hartree = energy%hartree + full_scaling*ener

      ! add scaled forces, restoring the old
      IF (calculate_forces) THEN
         nforce = 0
         DO i = 1, SIZE(force)
            force(i)%ch_pulay(:, :) = force(i)%ch_pulay(:, :)*full_scaling + &
                                      store_forces(1:3, nforce + 1:nforce + SIZE(force(i)%ch_pulay, 2))
            nforce = nforce + SIZE(force(i)%ch_pulay, 2)
         END DO
      END IF

      IF (.NOT. just_energy) THEN
         CALL pw_pool_create_pw(auxbas_pw_pool, v_sic_rspace%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_transfer(work_v%pw, v_sic_rspace%pw)
         ! also take into account the scaling (in addition to the volume element)
         CALL pw_scale(v_sic_rspace%pw, &
                       dft_control%sic_scaling_a*v_sic_rspace%pw%pw_grid%dvol)
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, work_rho%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, work_v%pw)

   END SUBROUTINE calc_v_sic_rspace

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho ...
! **************************************************************************************************
   SUBROUTINE print_densities(qs_env, rho)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho

      INTEGER                                            :: img, ispin, n_electrons, output_unit
      REAL(dp)                                           :: tot1_h, tot1_s, tot_rho_r, trace, &
                                                            trace_tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r_arr
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input, scf_section

      NULLIFY (qs_charges, qs_kind_set, cell, input, logger, scf_section, matrix_s, &
               dft_control, tot_rho_r_arr, rho_ao)

      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, qs_charges=qs_charges, &
                      input=input, &
                      matrix_s_kp=matrix_s, &
                      dft_control=dft_control)

      CALL get_qs_kind_set(qs_kind_set, nelectron=n_electrons)

      scf_section => section_vals_get_subs_vals(input, "DFT%SCF")
      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%TOTAL_DENSITIES", &
                                         extension=".scfLog")

      CALL qs_rho_get(rho, tot_rho_r=tot_rho_r_arr, rho_ao_kp=rho_ao)
      n_electrons = n_electrons - dft_control%charge
      tot_rho_r = accurate_sum(tot_rho_r_arr)

      trace = 0
      IF (BTEST(cp_print_key_should_output(logger%iter_info, scf_section, "PRINT%TOTAL_DENSITIES"), cp_p_file)) THEN
         DO ispin = 1, dft_control%nspins
            DO img = 1, dft_control%nimages
               CALL dbcsr_dot(rho_ao(ispin, img)%matrix, matrix_s(1, img)%matrix, trace_tmp)
               trace = trace + trace_tmp
            END DO
         END DO
      END IF

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T3,A,T41,F20.10)") "Trace(PS):", trace
         WRITE (UNIT=output_unit, FMT="((T3,A,T41,2F20.10))") &
            "Electronic density on regular grids: ", &
            tot_rho_r, &
            tot_rho_r + &
            REAL(n_electrons, dp), &
            "Core density on regular grids:", &
            qs_charges%total_rho_core_rspace, &
            qs_charges%total_rho_core_rspace - REAL(n_electrons + dft_control%charge, dp)
      END IF
      IF (dft_control%qs_control%gapw) THEN
         tot1_h = qs_charges%total_rho1_hard(1)
         tot1_s = qs_charges%total_rho1_soft(1)
         DO ispin = 2, dft_control%nspins
            tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
            tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
         END DO
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="((T3,A,T41,2F20.10))") &
               "Hard and soft densities (Lebedev):", &
               tot1_h, tot1_s
            WRITE (UNIT=output_unit, FMT="(T3,A,T41,F20.10)") &
               "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ", &
               tot_rho_r + tot1_h - tot1_s, &
               "Total charge density (r-space):      ", &
               tot_rho_r + tot1_h - tot1_s &
               + qs_charges%total_rho_core_rspace, &
               "Total Rho_soft + Rho0_soft (g-space):", &
               qs_charges%total_rho_gspace
         END IF
         qs_charges%background = tot_rho_r + tot1_h - tot1_s + &
                                 qs_charges%total_rho_core_rspace
      ELSE IF (dft_control%qs_control%gapw_xc) THEN
         tot1_h = qs_charges%total_rho1_hard(1)
         tot1_s = qs_charges%total_rho1_soft(1)
         DO ispin = 2, dft_control%nspins
            tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
            tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
         END DO
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T41,2F20.10))") &
               "Hard and soft densities (Lebedev):", &
               tot1_h, tot1_s
            WRITE (UNIT=output_unit, FMT="(T3,A,T41,F20.10)") &
               "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ", &
               accurate_sum(tot_rho_r_arr) + tot1_h - tot1_s
         END IF
         qs_charges%background = tot_rho_r + &
                                 qs_charges%total_rho_core_rspace
      ELSE
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T41,F20.10)") &
               "Total charge density on r-space grids:     ", &
               tot_rho_r + &
               qs_charges%total_rho_core_rspace, &
               "Total charge density g-space grids:     ", &
               qs_charges%total_rho_gspace
         END IF
         qs_charges%background = tot_rho_r + &
                                 qs_charges%total_rho_core_rspace
      END IF
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="()")
      qs_charges%background = qs_charges%background/cell%deth

      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%TOTAL_DENSITIES")

   END SUBROUTINE print_densities

! **************************************************************************************************
!> \brief Print detailed energies
!>
!> \param qs_env ...
!> \param dft_control ...
!> \param input ...
!> \param energy ...
!> \param mulliken_order_p ...
!> \par History
!>    refactoring 04.03.2011 [MI]
!> \author
! **************************************************************************************************
   SUBROUTINE print_detailed_energy(qs_env, dft_control, input, energy, mulliken_order_p)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(qs_energy_type), POINTER                      :: energy
      REAL(KIND=dp), INTENT(IN)                          :: mulliken_order_p

      INTEGER                                            :: bc, n, output_unit, psolver
      REAL(KIND=dp)                                      :: ddapc_order_p, implicit_ps_ehartree, &
                                                            s2_order_p
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(pw_env_type), POINTER                         :: pw_env

      logger => cp_get_default_logger()

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      psolver = pw_env%poisson_env%parameters%solver

      output_unit = cp_print_key_unit_nr(logger, input, "DFT%SCF%PRINT%DETAILED_ENERGY", &
                                         extension=".scfLog")
      IF (output_unit > 0) THEN
         IF (dft_control%do_admm) THEN
            WRITE (UNIT=output_unit, FMT="((T3,A,T60,F20.10))") &
               "Wfn fit exchange-correlation energy:            ", energy%exc_aux_fit
            IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
               WRITE (UNIT=output_unit, FMT="((T3,A,T60,F20.10))") &
                  "Wfn fit soft/hard atomic rho1 Exc contribution: ", energy%exc1_aux_fit
            END IF
         END IF
         IF (dft_control%do_admm) THEN
            IF (psolver .EQ. pw_poisson_implicit) THEN
               implicit_ps_ehartree = pw_env%poisson_env%implicit_env%ehartree
               bc = pw_env%poisson_env%parameters%ps_implicit_params%boundary_condition
               SELECT CASE (bc)
               CASE (MIXED_PERIODIC_BC, MIXED_BC)
                  WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                     "Core Hamiltonian energy:                       ", energy%core, &
                     "Hartree energy:                                ", implicit_ps_ehartree, &
                     "Electric enthalpy:                             ", energy%hartree, &
                     "Exchange-correlation energy:                   ", energy%exc + energy%exc_aux_fit
               CASE (PERIODIC_BC, NEUMANN_BC)
                  WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                     "Core Hamiltonian energy:                       ", energy%core, &
                     "Hartree energy:                                ", energy%hartree, &
                     "Exchange-correlation energy:                   ", energy%exc + energy%exc_aux_fit
               END SELECT
            ELSE
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                  "Core Hamiltonian energy:                       ", energy%core, &
                  "Hartree energy:                                ", energy%hartree, &
                  "Exchange-correlation energy:                   ", energy%exc + energy%exc_aux_fit
            END IF
         ELSE
!ZMP to print some variables at each step
            IF (dft_control%apply_external_density) THEN
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                  "DOING ZMP CALCULATION FROM EXTERNAL DENSITY    "
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                  "Core Hamiltonian energy:                       ", energy%core, &
                  "Hartree energy:                                ", energy%hartree
            ELSE IF (dft_control%apply_external_vxc) THEN
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                  "DOING ZMP READING EXTERNAL VXC                 "
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                  "Core Hamiltonian energy:                       ", energy%core, &
                  "Hartree energy:                                ", energy%hartree
            ELSE
               IF (psolver .EQ. pw_poisson_implicit) THEN
                  implicit_ps_ehartree = pw_env%poisson_env%implicit_env%ehartree
                  bc = pw_env%poisson_env%parameters%ps_implicit_params%boundary_condition
                  SELECT CASE (bc)
                  CASE (MIXED_PERIODIC_BC, MIXED_BC)
                     WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                        "Core Hamiltonian energy:                       ", energy%core, &
                        "Hartree energy:                                ", implicit_ps_ehartree, &
                        "Electric enthalpy:                             ", energy%hartree, &
                        "Exchange-correlation energy:                   ", energy%exc
                  CASE (PERIODIC_BC, NEUMANN_BC)
                     WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                        "Core Hamiltonian energy:                       ", energy%core, &
                        "Hartree energy:                                ", energy%hartree, &
                        "Exchange-correlation energy:                   ", energy%exc
                  END SELECT
               ELSE
                  WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
                     "Core Hamiltonian energy:                       ", energy%core, &
                     "Hartree energy:                                ", energy%hartree, &
                     "Exchange-correlation energy:                   ", energy%exc
               END IF
            END IF
         END IF

         IF (dft_control%apply_external_density) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
               "Integral of the (density * v_xc):              ", energy%exc
         END IF

         IF (energy%e_hartree /= 0.0_dp) &
            WRITE (UNIT=output_unit, FMT="(T3,A,T61,F20.10)") &
            "Coulomb (electron-electron) energy:            ", energy%e_hartree
         IF (energy%dispersion /= 0.0_dp) &
            WRITE (UNIT=output_unit, FMT="(T3,A,T61,F20.10)") &
            "Dispersion energy:                             ", energy%dispersion
         IF (energy%gcp /= 0.0_dp) &
            WRITE (UNIT=output_unit, FMT="(T3,A,T61,F20.10)") &
            "gCP energy:                                    ", energy%gcp
         IF (dft_control%qs_control%gapw) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
               "GAPW| Exc from hard and soft atomic rho1:      ", energy%exc1 + energy%exc1_aux_fit, &
               "GAPW| local Eh = 1 center integrals:           ", energy%hartree_1c
         END IF
         IF (dft_control%qs_control%gapw_xc) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T61,F20.10))") &
               "GAPW| Exc from hard and soft atomic rho1:      ", energy%exc1 + energy%exc1_aux_fit
         END IF
         IF (dft_control%dft_plus_u) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T61,F20.10)") &
               "DFT+U energy:", energy%dft_plus_u
         END IF
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T61,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
            IF (qs_env%qmmm_env_qm%image_charge) THEN
               WRITE (UNIT=output_unit, FMT="(T3,A,T61,F20.10)") &
                  "QM/MM image charge energy:                ", energy%image_charge
            END IF
         END IF
         IF (dft_control%qs_control%mulliken_restraint) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T41,2F20.10)") &
               "Mulliken restraint (order_p,energy) : ", mulliken_order_p, energy%mulliken
         END IF
         IF (dft_control%qs_control%ddapc_restraint) THEN
            DO n = 1, SIZE(dft_control%qs_control%ddapc_restraint_control)
               ddapc_order_p = &
                  dft_control%qs_control%ddapc_restraint_control(n)%ddapc_restraint_control%ddapc_order_p
               WRITE (UNIT=output_unit, FMT="(T3,A,T41,2F20.10)") &
                  "DDAPC restraint (order_p,energy) : ", ddapc_order_p, energy%ddapc_restraint(n)
            END DO
         END IF
         IF (dft_control%qs_control%s2_restraint) THEN
            s2_order_p = dft_control%qs_control%s2_restraint_control%s2_order_p
            WRITE (UNIT=output_unit, FMT="(T3,A,T41,2F20.10)") &
               "S2 restraint (order_p,energy) : ", s2_order_p, energy%s2_restraint
         END IF

      END IF ! output_unit
      CALL cp_print_key_finished_output(output_unit, logger, input, &
                                        "DFT%SCF%PRINT%DETAILED_ENERGY")

   END SUBROUTINE print_detailed_energy

! **************************************************************************************************
!> \brief compute matrix_vxc, defined via the potential created by qs_vxc_create
!>        ignores things like tau functional, gapw, sic, ...
!>         so only OK for GGA & GPW right now
!> \param qs_env ...
!> \param v_rspace ...
!> \param matrix_vxc ...
!> \par History
!>    created 23.10.2012 [Joost VandeVondele]
!> \author
! **************************************************************************************************
   SUBROUTINE compute_matrix_vxc(qs_env, v_rspace, matrix_vxc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vxc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_matrix_vxc'

      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: gapw
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      ! create the matrix using matrix_ks as a template
      IF (ASSOCIATED(matrix_vxc)) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_vxc)
      END IF
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
      ALLOCATE (matrix_vxc(SIZE(matrix_ks)))
      DO ispin = 1, SIZE(matrix_ks)
         NULLIFY (matrix_vxc(ispin)%matrix)
         CALL dbcsr_init_p(matrix_vxc(ispin)%matrix)
         CALL dbcsr_copy(matrix_vxc(ispin)%matrix, matrix_ks(ispin)%matrix, &
                         name="Matrix VXC of spin "//cp_to_string(ispin))
         CALL dbcsr_set(matrix_vxc(ispin)%matrix, 0.0_dp)
      END DO

      ! and integrate
      CALL get_qs_env(qs_env, dft_control=dft_control)
      gapw = dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc
      DO ispin = 1, SIZE(matrix_ks)
         CALL integrate_v_rspace(v_rspace=v_rspace(ispin), &
                                 hmat=matrix_vxc(ispin), &
                                 qs_env=qs_env, &
                                 calculate_forces=.FALSE., &
                                 gapw=gapw)
         ! scale by the volume element... should really become part of integrate_v_rspace
         CALL dbcsr_scale(matrix_vxc(ispin)%matrix, v_rspace(ispin)%pw%pw_grid%dvol)
      END DO

      CALL timestop(handle)

   END SUBROUTINE compute_matrix_vxc

! **************************************************************************************************
!> \brief Sum up all potentials defined  on the grid and integrate
!>
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param my_rho ...
!> \param vppl_rspace ...
!> \param v_rspace_new ...
!> \param v_rspace_new_aux_fit ...
!> \param v_tau_rspace ...
!> \param v_tau_rspace_aux_fit ...
!> \param v_efield_rspace ...
!> \param v_sic_rspace ...
!> \param v_spin_ddapc_rest_r ...
!> \param v_sccs_rspace ...
!> \param v_rspace_embed ...
!> \param cdft_control ...
!> \param calculate_forces ...
!> \par History
!>      - refactoring 04.03.2011 [MI]
!>      - SCCS implementation (16.10.2013,MK)
!> \author
! **************************************************************************************************
   SUBROUTINE sum_up_and_integrate(qs_env, ks_matrix, rho, my_rho, &
                                   vppl_rspace, v_rspace_new, &
                                   v_rspace_new_aux_fit, v_tau_rspace, &
                                   v_tau_rspace_aux_fit, v_efield_rspace, &
                                   v_sic_rspace, v_spin_ddapc_rest_r, &
                                   v_sccs_rspace, v_rspace_embed, cdft_control, &
                                   calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: my_rho
      TYPE(pw_p_type), POINTER                           :: vppl_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace_new, v_rspace_new_aux_fit, &
                                                            v_tau_rspace, v_tau_rspace_aux_fit
      TYPE(pw_p_type)                                    :: v_efield_rspace, v_sic_rspace, &
                                                            v_spin_ddapc_rest_r, v_sccs_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace_embed
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      LOGICAL, INTENT(in)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sum_up_and_integrate'

      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER                                            :: handle, igroup, ikind, ispin, nkind, &
                                                            nspins
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_ppl, gapw, gapw_xc, lrigpw, rigpw
      REAL(dp)                                           :: dvol, sign
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ksmat, matrix_ks_aux_fit, &
                                                            matrix_ks_aux_fit_dft, rho_ao, &
                                                            rho_ao_aux, rho_ao_nokp, smat
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: v_rspace
      TYPE(pw_p_type), POINTER                           :: vee
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool, dft_control, pw_env, matrix_ks_aux_fit, &
               v_rspace%pw, rho_aux_fit, vee, rho_ao, rho_ao_kp, rho_ao_aux, &
               ksmat, matrix_ks_aux_fit_dft, lri_env, lri_density, atomic_kind_set, &
               rho_ao_nokp, ks_env, admm_env, task_list)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env, &
                      matrix_ks_aux_fit=matrix_ks_aux_fit, &
                      matrix_ks_aux_fit_dft=matrix_ks_aux_fit_dft, &
                      v_hartree_rspace=v_rspace%pw, &
                      rho_aux_fit=rho_aux_fit, &
                      vee=vee)

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
      CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux)
      CALL pw_env_get(pw_env, poisson_env=poisson_env, auxbas_pw_pool=auxbas_pw_pool)
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
      do_ppl = dft_control%qs_control%do_ppl_method == do_ppl_grid

      rigpw = dft_control%qs_control%rigpw
      lrigpw = dft_control%qs_control%lrigpw
      IF (lrigpw .OR. rigpw) THEN
         CALL get_qs_env(qs_env, &
                         lri_env=lri_env, &
                         lri_density=lri_density, &
                         atomic_kind_set=atomic_kind_set)
      END IF

      nspins = dft_control%nspins

      ! sum up potentials and integrate
      IF (ASSOCIATED(v_rspace_new)) THEN
         DO ispin = 1, nspins
            IF (gapw_xc) THEN
               ! SIC not implemented (or at least not tested)
               CPASSERT(dft_control%sic_method_id == sic_none)
               !Only the xc potential, because it has to be integrated with the soft basis
               v_rspace_new(ispin)%pw%cr3d = &
                  v_rspace_new(ispin)%pw%pw_grid%dvol* &
                  v_rspace_new(ispin)%pw%cr3d

               ! add the xc  part due to v_rspace soft
               rho_ao => rho_ao_kp(ispin, :)
               ksmat => ks_matrix(ispin, :)
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat_kp=rho_ao, hmat_kp=ksmat, &
                                       qs_env=qs_env, &
                                       calculate_forces=calculate_forces, &
                                       gapw=gapw_xc)

               ! Now the Hartree potential to be integrated with the full basis
               v_rspace_new(ispin)%pw%cr3d = v_rspace%pw%cr3d
            ELSE
               ! Add v_hartree + v_xc = v_rspace_new
               v_rspace_new(ispin)%pw%cr3d = &
                  v_rspace_new(ispin)%pw%pw_grid%dvol* &
                  v_rspace_new(ispin)%pw%cr3d + v_rspace%pw%cr3d
            END IF ! gapw_xc
            IF (dft_control%qs_control%ddapc_explicit_potential) THEN
               IF (dft_control%qs_control%ddapc_restraint_is_spin) THEN
                  IF (ispin == 1) THEN
                     v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d &
                                                   + v_spin_ddapc_rest_r%pw%cr3d
                  ELSE
                     v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d &
                                                   - v_spin_ddapc_rest_r%pw%cr3d
                  END IF
               ELSE
                  v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d &
                                                + v_spin_ddapc_rest_r%pw%cr3d
               END IF
            END IF
            ! CDFT constraint contribution
            IF (dft_control%qs_control%cdft) THEN
               DO igroup = 1, SIZE(cdft_control%group)
                  SELECT CASE (cdft_control%group(igroup)%constraint_type)
                  CASE (cdft_charge_constraint)
                     sign = 1.0_dp
                  CASE (cdft_magnetization_constraint)
                     IF (ispin == 1) THEN
                        sign = 1.0_dp
                     ELSE
                        sign = -1.0_dp
                     END IF
                  CASE (cdft_alpha_constraint)
                     sign = 1.0_dp
                     IF (ispin == 2) CYCLE
                  CASE (cdft_beta_constraint)
                     sign = 1.0_dp
                     IF (ispin == 1) CYCLE
                  CASE DEFAULT
                     CPABORT("Unknown constraint type.")
                  END SELECT
                  v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d &
                                                + sign*cdft_control%group(igroup)%weight%pw%cr3d* &
                                                cdft_control%strength(igroup)
               END DO
            END IF
            ! The efield contribution
            IF (dft_control%apply_efield_field) THEN
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d + &
                                             v_efield_rspace%pw%cr3d
            END IF
            ! functional derivative of the Hartree energy wrt the density in the presence of dielectric
            ! (vhartree + v_eps); v_eps is nonzero only if the dielectric constant is defind as a function
            ! of the charge density
            IF (poisson_env%parameters%solver .EQ. pw_poisson_implicit) THEN
               dvol = poisson_env%implicit_env%v_eps%pw_grid%dvol
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d + &
                                             poisson_env%implicit_env%v_eps%cr3d*dvol
            END IF
            ! Add SCCS contribution
            IF (dft_control%do_sccs) THEN
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d + &
                                             v_sccs_rspace%pw%cr3d
            END IF
            ! External electrostatic potential
            IF (dft_control%apply_external_potential) THEN
               CALL qmmm_modify_hartree_pot(v_hartree=v_rspace_new(ispin), &
                                            v_qmmm=vee, scale=-1.0_dp)
            END IF
            IF (do_ppl) THEN
               CPASSERT(.NOT. gapw)
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d + &
                                             vppl_rspace%pw%cr3d*vppl_rspace%pw%pw_grid%dvol
            END IF
            ! the electrostatic sic contribution
            SELECT CASE (dft_control%sic_method_id)
            CASE (sic_none)
               !
            CASE (sic_mauri_us, sic_mauri_spz)
               IF (ispin == 1) THEN
                  v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d - &
                                                v_sic_rspace%pw%cr3d
               ELSE
                  v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d + &
                                                v_sic_rspace%pw%cr3d
               END IF
            CASE (sic_ad)
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d - v_sic_rspace%pw%cr3d
            CASE (sic_eo)
               ! NOTHING TO BE DONE
            END SELECT
            ! DFT embedding
            IF (dft_control%apply_embed_pot) THEN
               v_rspace_new(ispin)%pw%cr3d = v_rspace_new(ispin)%pw%cr3d + &
                                             v_rspace_embed(ispin)%pw%pw_grid%dvol* &
                                             v_rspace_embed(ispin)%pw%cr3d
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_embed(ispin)%pw)
            END IF
            IF (lrigpw) THEN
               lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
               CALL get_qs_env(qs_env, nkind=nkind, para_env=para_env)
               DO ikind = 1, nkind
                  lri_v_int(ikind)%v_int = 0.0_dp
               END DO
               CALL integrate_v_rspace_one_center(v_rspace_new(ispin), qs_env, &
                                                  lri_v_int, calculate_forces, "LRI_AUX")
               DO ikind = 1, nkind
                  CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
               END DO
               IF (lri_env%exact_1c_terms) THEN
                  rho_ao => my_rho(ispin, :)
                  ksmat => ks_matrix(ispin, :)
                  CALL integrate_v_rspace_diagonal(v_rspace_new(ispin), ksmat(1)%matrix, &
                                                   rho_ao(1)%matrix, qs_env, &
                                                   calculate_forces, "ORB")
               END IF
               IF (lri_env%ppl_ri) THEN
                  CALL v_int_ppl_update(qs_env, lri_v_int, calculate_forces)
               END IF
            ELSEIF (rigpw) THEN
               lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
               CALL get_qs_env(qs_env, nkind=nkind, para_env=para_env)
               DO ikind = 1, nkind
                  lri_v_int(ikind)%v_int = 0.0_dp
               END DO
               CALL integrate_v_rspace_one_center(v_rspace_new(ispin), qs_env, &
                                                  lri_v_int, calculate_forces, "RI_HXC")
               DO ikind = 1, nkind
                  CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
               END DO
            ELSE
               rho_ao => my_rho(ispin, :)
               ksmat => ks_matrix(ispin, :)
               CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                       pmat_kp=rho_ao, hmat_kp=ksmat, &
                                       qs_env=qs_env, &
                                       calculate_forces=calculate_forces, &
                                       gapw=gapw)
            END IF
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_new(ispin)%pw)
         END DO ! ispin

         SELECT CASE (dft_control%sic_method_id)
         CASE (sic_none)
         CASE (sic_mauri_us, sic_mauri_spz, sic_ad)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_sic_rspace%pw)
         END SELECT
         DEALLOCATE (v_rspace_new)

      ELSE
         ! not implemented (or at least not tested)
         CPASSERT(dft_control%sic_method_id == sic_none)
         CPASSERT(.NOT. dft_control%qs_control%ddapc_restraint_is_spin)
         DO ispin = 1, nspins
            ! the efield contribution
            IF (dft_control%apply_efield_field) THEN
               v_rspace%pw%cr3d = v_rspace%pw%cr3d + v_efield_rspace%pw%cr3d
            END IF
            ! extra contribution attributed to the dependency of the dielectric constant to the charge density
            IF (poisson_env%parameters%solver .EQ. pw_poisson_implicit) THEN
               dvol = poisson_env%implicit_env%v_eps%pw_grid%dvol
               v_rspace%pw%cr3d = v_rspace%pw%cr3d + poisson_env%implicit_env%v_eps%cr3d*dvol
            END IF
            ! Add SCCS contribution
            IF (dft_control%do_sccs) THEN
               v_rspace%pw%cr3d = v_rspace%pw%cr3d + v_sccs_rspace%pw%cr3d
            END IF
            ! DFT embedding
            IF (dft_control%apply_embed_pot) THEN
               v_rspace%pw%cr3d = v_rspace%pw%cr3d + &
                                  v_rspace_embed(ispin)%pw%pw_grid%dvol* &
                                  v_rspace_embed(ispin)%pw%cr3d
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_embed(ispin)%pw)
            END IF
            IF (lrigpw) THEN
               lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
               CALL get_qs_env(qs_env, nkind=nkind, para_env=para_env)
               DO ikind = 1, nkind
                  lri_v_int(ikind)%v_int = 0.0_dp
               END DO
               CALL integrate_v_rspace_one_center(v_rspace, qs_env, &
                                                  lri_v_int, calculate_forces, "LRI_AUX")
               DO ikind = 1, nkind
                  CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
               END DO
               IF (lri_env%exact_1c_terms) THEN
                  rho_ao => my_rho(ispin, :)
                  ksmat => ks_matrix(ispin, :)
                  CALL integrate_v_rspace_diagonal(v_rspace, ksmat(1)%matrix, &
                                                   rho_ao(1)%matrix, qs_env, &
                                                   calculate_forces, "ORB")
               END IF
               IF (lri_env%ppl_ri) THEN
                  CALL v_int_ppl_update(qs_env, lri_v_int, calculate_forces)
               END IF
            ELSEIF (rigpw) THEN
               lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
               CALL get_qs_env(qs_env, nkind=nkind, para_env=para_env)
               DO ikind = 1, nkind
                  lri_v_int(ikind)%v_int = 0.0_dp
               END DO
               CALL integrate_v_rspace_one_center(v_rspace, qs_env, &
                                                  lri_v_int, calculate_forces, "RI_HXC")
               DO ikind = 1, nkind
                  CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
               END DO
            ELSE
               rho_ao => my_rho(ispin, :)
               ksmat => ks_matrix(ispin, :)
               CALL integrate_v_rspace(v_rspace=v_rspace, &
                                       pmat_kp=rho_ao, &
                                       hmat_kp=ksmat, &
                                       qs_env=qs_env, &
                                       calculate_forces=calculate_forces, &
                                       gapw=gapw)
            END IF
         END DO
      END IF ! ASSOCIATED(v_rspace_new)

      ! **** LRIGPW: KS matrix from integrated potential
      IF (lrigpw) THEN
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         DO ispin = 1, nspins
            ksmat => ks_matrix(ispin, :)
            CALL calculate_lri_ks_matrix(lri_env, lri_v_int, ksmat, atomic_kind_set, &
                                         cell_to_index=cell_to_index)
         END DO
         IF (calculate_forces) THEN
            CALL calculate_lri_forces(lri_env, lri_density, qs_env, rho_ao_kp, atomic_kind_set)
         END IF
      ELSEIF (rigpw) THEN
         CALL get_qs_env(qs_env, matrix_s=smat)
         DO ispin = 1, nspins
            CALL calculate_ri_ks_matrix(lri_env, lri_v_int, ks_matrix(ispin, 1)%matrix, &
                                        smat(1)%matrix, atomic_kind_set, ispin)
         END DO
         IF (calculate_forces) THEN
            rho_ao_nokp => rho_ao_kp(:, 1)
            CALL calculate_ri_forces(lri_env, lri_density, qs_env, rho_ao_nokp, atomic_kind_set)
         END IF
      END IF

      IF (ASSOCIATED(v_tau_rspace)) THEN
         IF (lrigpw .OR. rigpw) THEN
            CPABORT("LRIGPW/RIGPW not implemented for meta-GGAs")
         END IF
         DO ispin = 1, nspins
            v_tau_rspace(ispin)%pw%cr3d = &
               v_tau_rspace(ispin)%pw%pw_grid%dvol* &
               v_tau_rspace(ispin)%pw%cr3d

            rho_ao => rho_ao_kp(ispin, :)
            ksmat => ks_matrix(ispin, :)
            CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin), &
                                    pmat_kp=rho_ao, hmat_kp=ksmat, &
                                    qs_env=qs_env, &
                                    calculate_forces=calculate_forces, compute_tau=.TRUE., &
                                    gapw=gapw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
         END DO
         DEALLOCATE (v_tau_rspace)
      END IF

      ! Add contributions from ADMM if requested
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         IF (ASSOCIATED(v_rspace_new_aux_fit)) THEN
            DO ispin = 1, nspins
               ! Calculate the xc potential
               v_rspace_new_aux_fit(ispin)%pw%cr3d = &
                  v_rspace_new_aux_fit(ispin)%pw%pw_grid%dvol* &
                  v_rspace_new_aux_fit(ispin)%pw%cr3d

               ! set matrix_ks_aux_fit_dft = matrix_ks_aux_fit(k_HF)
               CALL dbcsr_copy(matrix_ks_aux_fit_dft(ispin)%matrix, matrix_ks_aux_fit(ispin)%matrix, &
                               name="DFT exch. part of matrix_ks_aux_fit")

               ! Add potential to ks_matrix aux_fit, skip integration if no DFT correction

               IF (admm_env%aux_exch_func .NE. do_admm_aux_exch_func_none) THEN

                  !GPW by default. IF GAPW, then take relevant task list and basis
                  CALL get_qs_env(qs_env, task_list_aux_fit=task_list)
                  basis_type = "AUX_FIT"
                  IF (admm_env%do_gapw) THEN
                     task_list => admm_env%admm_gapw_env%task_list
                     basis_type = "AUX_FIT_SOFT"
                  END IF

                  CALL integrate_v_rspace(v_rspace=v_rspace_new_aux_fit(ispin), &
                                          pmat=rho_ao_aux(ispin), &
                                          hmat=matrix_ks_aux_fit(ispin), &
                                          qs_env=qs_env, &
                                          calculate_forces=calculate_forces, &
                                          force_adm=.TRUE., &
                                          ispin=ispin, &
                                          gapw=.FALSE., & !even if actual GAPW calculation, want to use AUX_FIT_SOFT
                                          basis_type=basis_type, &
                                          task_list_external=task_list)
               END IF

               ! matrix_ks_aux_fit_dft(x_DFT)=matrix_ks_aux_fit_dft(old,k_HF)-matrix_ks_aux_fit(k_HF-x_DFT)
               CALL dbcsr_add(matrix_ks_aux_fit_dft(ispin)%matrix, &
                              matrix_ks_aux_fit(ispin)%matrix, 1.0_dp, -1.0_dp)

               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_new_aux_fit(ispin)%pw)
            END DO
            DEALLOCATE (v_rspace_new_aux_fit)
         END IF
         ! Clean up v_tau_rspace_aux_fit, which is actually not needed
         IF (ASSOCIATED(v_tau_rspace_aux_fit)) THEN
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace_aux_fit(ispin)%pw)
            END DO
            DEALLOCATE (v_tau_rspace_aux_fit)
         END IF
      END IF

      IF (dft_control%apply_embed_pot) DEALLOCATE (v_rspace_embed)

      CALL timestop(handle)

   END SUBROUTINE sum_up_and_integrate

!**************************************************************************
!> \brief Calculate the ZMP potential and energy as in Zhao, Morrison Parr
!> PRA 50i, 2138 (1994)
!> V_c^\lambda defined as int_rho-rho_0/r-r' or rho-rho_0 times a Lagrange
!> multiplier, plus Fermi-Amaldi potential that should give the V_xc in the
!> limit \lambda --> \infty
!>
!> \param qs_env ...
!> \param v_rspace_new ...
!> \param rho ...
!> \param exc ...
!> \author D. Varsano  [daniele.varsano@nano.cnr.it]
! **************************************************************************************************
   SUBROUTINE calculate_zmp_potential(qs_env, v_rspace_new, rho, exc)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace_new
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(KIND=dp)                                      :: exc

      CHARACTER(*), PARAMETER :: routineN = 'calculate_zmp_potential'

      INTEGER                                            :: handle, my_val, nelectron, nspins, stat
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      LOGICAL                                            :: do_zmp_read, fermi_amaldi
      REAL(KIND=dp)                                      :: factor, lambda, total_rho
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_ext_r
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_eff_gspace, v_xc_gspace, v_xc_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ext_g, rho_g, rho_r
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: ext_den_section, input

!, v_h_gspace, &

      CALL timeset(routineN, handle)
      NULLIFY (auxbas_pw_pool)
      NULLIFY (pw_env)
      NULLIFY (poisson_env)
      NULLIFY (v_rspace_new)
      NULLIFY (dft_control)
      NULLIFY (rho_r, rho_g, tot_rho_ext_r, rho_ext_g)
      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env, &
                      ks_env=ks_env, &
                      rho=rho, &
                      input=input, &
                      nelectron_spin=nelectron_spin, &
                      dft_control=dft_control)
      CALL pw_env_get(pw_env=pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g)
      nspins = 1
      ALLOCATE (v_rspace_new(nspins), stat=stat)
      CPASSERT(stat == 0)
      CALL pw_pool_create_pw(pool=auxbas_pw_pool, pw=v_rspace_new(1)%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(pool=auxbas_pw_pool, pw=v_xc_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL pw_zero(v_rspace_new(1)%pw)
      do_zmp_read = dft_control%apply_external_vxc
      IF (do_zmp_read) THEN
         CALL pw_copy(qs_env%external_vxc%pw, v_rspace_new(1)%pw)
         exc = 0.0_dp
         exc = accurate_dot_product(v_rspace_new(1)%pw%cr3d, rho_r(1)%pw%cr3d)* &
               v_rspace_new(1)%pw%pw_grid%dvol
      ELSE
         CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                                pw=rho_eff_gspace%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                                pw=v_xc_gspace%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_zero(rho_eff_gspace%pw)
         CALL pw_zero(v_xc_gspace%pw)
         CALL pw_zero(v_xc_rspace%pw)
         factor = pw_integrate_function(rho_g(1)%pw)
         CALL qs_rho_get(qs_env%rho_external, &
                         rho_g=rho_ext_g, &
                         tot_rho_r=tot_rho_ext_r)
         factor = tot_rho_ext_r(1)/factor

         CALL pw_axpy(rho_g(1)%pw, rho_eff_gspace%pw, alpha=factor)
         CALL pw_axpy(rho_ext_g(1)%pw, rho_eff_gspace%pw, alpha=-1.0_dp)
         total_rho = pw_integrate_function(rho_eff_gspace%pw, isign=1)
         ext_den_section => section_vals_get_subs_vals(input, "DFT%EXTERNAL_DENSITY")
         CALL section_vals_val_get(ext_den_section, "LAMBDA", r_val=lambda)
         CALL section_vals_val_get(ext_den_section, "ZMP_CONSTRAINT", i_val=my_val)
         CALL section_vals_val_get(ext_den_section, "FERMI_AMALDI", l_val=fermi_amaldi)

         CALL pw_scale(rho_eff_gspace%pw, a=lambda)
         nelectron = nelectron_spin(1)
         factor = -1.0_dp/nelectron
         CALL pw_axpy(rho_g(1)%pw, rho_eff_gspace%pw, alpha=factor)

         CALL pw_poisson_solve(poisson_env, rho_eff_gspace%pw, vhartree=v_xc_gspace%pw)
         CALL pw_transfer(v_xc_gspace%pw, v_rspace_new(1)%pw)
         CALL pw_copy(v_rspace_new(1)%pw, v_xc_rspace%pw)

         exc = 0.0_dp
         exc = accurate_sum(v_rspace_new(1)%pw%cr3d*rho_r(1)%pw%cr3d)* &
               v_rspace_new(1)%pw%pw_grid%dvol
         IF (v_rspace_new(1)%pw%pw_grid%para%mode == PW_MODE_DISTRIBUTED) THEN
            CALL mp_sum(exc, v_rspace_new(1)%pw%pw_grid%para%group)
         END IF

!     IF (v_rspace_new(1)%pw%pw_grid%para%my_pos==0) &
!           WRITE(*,FMT="(T3,A,T61,F20.10)") "ZMP| Integral of (electron density)*(v_xc):    " , exc
!Note that this is not the xc energy but \int(\rho*v_xc)
!Vxc---> v_rspace_new
!Exc---> energy%exc
         CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                   rho_eff_gspace%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                   v_xc_gspace%pw)
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc_rspace%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_zmp_potential

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho ...
!> \param v_rspace_embed ...
!> \param dft_control ...
!> \param embed_corr ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE get_embed_potential_energy(qs_env, rho, v_rspace_embed, dft_control, embed_corr, &
                                         just_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace_embed
      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(KIND=dp)                                      :: embed_corr
      LOGICAL                                            :: just_energy

      CHARACTER(*), PARAMETER :: routineN = 'get_embed_potential_energy'

      INTEGER                                            :: handle, ispin, ns1, ns2
      REAL(KIND=dp)                                      :: embed_corr_local
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool)
      NULLIFY (pw_env)
      NULLIFY (rho_r)
      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env, &
                      rho=rho)
      CALL pw_env_get(pw_env=pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL qs_rho_get(rho, rho_r=rho_r)
      ALLOCATE (v_rspace_embed(dft_control%nspins))

      embed_corr = 0.0_dp

      DO ispin = 1, dft_control%nspins
         CALL pw_pool_create_pw(pool=auxbas_pw_pool, pw=v_rspace_embed(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_rspace_embed(ispin)%pw)
         ns1 = SIZE(qs_env%embed_pot%pw%cr3d)
         ns2 = SIZE(v_rspace_embed(ispin)%pw%cr3d)
         IF (ns1 .NE. ns2) CPABORT("Subsystem grids must be identical")

         v_rspace_embed(ispin)%pw%cr3d(:, :, :) = qs_env%embed_pot%pw%cr3d
         embed_corr_local = 0.0_dp

         ! Spin embedding potential in open-shell case
         IF (dft_control%nspins .EQ. 2) THEN
            ns1 = SIZE(qs_env%spin_embed_pot%pw%cr3d)
            IF (ns1 .NE. ns2) CPABORT("Subsystem grids must be identical")
            IF (ispin .EQ. 1) v_rspace_embed(ispin)%pw%cr3d(:, :, :) = &
               v_rspace_embed(ispin)%pw%cr3d(:, :, :) + qs_env%spin_embed_pot%pw%cr3d
            IF (ispin .EQ. 2) v_rspace_embed(ispin)%pw%cr3d(:, :, :) = &
               v_rspace_embed(ispin)%pw%cr3d(:, :, :) - qs_env%spin_embed_pot%pw%cr3d
         END IF
         ! Integrate the density*potential
         embed_corr_local = accurate_dot_product(v_rspace_embed(ispin)%pw%cr3d, rho_r(ispin)%pw%cr3d)* &
                            v_rspace_embed(ispin)%pw%pw_grid%dvol

         IF (v_rspace_embed(ispin)%pw%pw_grid%para%mode == PW_MODE_DISTRIBUTED) THEN
            CALL mp_sum(embed_corr_local, v_rspace_embed(ispin)%pw%pw_grid%para%group)
         END IF

         embed_corr = embed_corr + embed_corr_local

      END DO

      ! If only energy requiested we delete the potential
      IF (just_energy) THEN
         DO ispin = 1, dft_control%nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_embed(ispin)%pw)
         END DO
         DEALLOCATE (v_rspace_embed)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_embed_potential_energy

END MODULE qs_ks_utils
