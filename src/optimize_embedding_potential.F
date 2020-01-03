!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE optimize_embedding_potential

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_element, cp_fm_get_info, cp_fm_release, &
        cp_fm_set_all, cp_fm_to_fm, cp_fm_to_fm_submat, cp_fm_type, cp_fm_write_unformatted
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_cube_to_pw,&
                                              cp_pw_to_cube,&
                                              cp_pw_to_simple_volumetric
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE embed_types,                     ONLY: opt_embed_pot_type
   USE force_env_types,                 ONLY: force_env_type
   USE input_constants,                 ONLY: &
        embed_diff, embed_fa, embed_grid_angstrom, embed_grid_bohr, embed_level_shift, embed_none, &
        embed_quasi_newton, embed_resp, embed_steep_desc
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_ivals,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE lri_environment_types,           ONLY: lri_kind_type
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_max,&
                                              mp_sum
   USE mixed_environment_utils,         ONLY: get_subsys_map_index
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: &
        pw_axpy, pw_copy, pw_derive, pw_dr2, pw_integral_ab, pw_integrate_function, pw_scale, &
        pw_transfer, pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_release,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_rho_resp_all,&
                                              calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_integrate_potential_single,   ONLY: integrate_v_rspace_one_center
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE xc,                              ONLY: smooth_cutoff
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_setall,&
                                              xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_embedding_potential'

   PUBLIC :: prepare_embed_opt, init_embed_pot, release_opt_embed, calculate_embed_pot_grad, &
             opt_embed_step, print_rho_diff, step_control, max_dens_diff, print_emb_opt_info, &
             conv_check_embed, make_subsys_embed_pot, print_embed_restart, find_aux_dimen, &
             read_embed_pot, understand_spin_states, given_embed_pot, print_rho_spin_diff, &
             print_pot_simple_grid, get_prev_density, get_max_subsys_diff, Coulomb_guess

CONTAINS

! **************************************************************************************************
!> \brief Find out whether we need to swap alpha- and beta- spind densities in the second subsystem
!> \brief It's only needed because by default alpha-spins go first in a subsystem.
!> \brief By swapping we impose the constraint:
!> \brief rho_1(alpha) + rho_2(alpha) = rho_total(alpha)
!> \brief rho_1(beta) + rho_2(beta) = rho_total(beta)
!> \param force_env ...
!> \param ref_subsys_number ...
!> \param change_spin ...
!> \param open_shell_embed ...
!> \param all_nspins ...
!> \return ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE understand_spin_states(force_env, ref_subsys_number, change_spin, open_shell_embed, all_nspins)
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER                                            :: ref_subsys_number
      LOGICAL                                            :: change_spin, open_shell_embed
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: all_nspins

      INTEGER                                            :: i_force_eval, nspins, sub_spin_1, &
                                                            sub_spin_2, total_spin
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      INTEGER, DIMENSION(2, 3)                           :: all_spins
      TYPE(dft_control_type), POINTER                    :: dft_control

      change_spin = .FALSE.
      open_shell_embed = .FALSE.
      ALLOCATE (all_nspins(ref_subsys_number))
      IF (ref_subsys_number .EQ. 3) THEN
         all_spins = 0
         DO i_force_eval = 1, ref_subsys_number
            CALL get_qs_env(qs_env=force_env%sub_force_env(i_force_eval)%force_env%qs_env, &
                            nelectron_spin=nelectron_spin, dft_control=dft_control)
            all_spins(:, i_force_eval) = nelectron_spin
            nspins = dft_control%nspins
            all_nspins(i_force_eval) = nspins
         ENDDO

         ! Find out whether we need a spin-dependend embedding potential
         IF (.NOT. ((all_nspins(1) .EQ. 1) .AND. (all_nspins(2) .EQ. 1) .AND. (all_nspins(3) .EQ. 1))) &
            open_shell_embed = .TRUE.

         ! If it's open shell, we need to check spin states
         IF (open_shell_embed) THEN

            IF (all_nspins(3) .EQ. 1) THEN
               total_spin = 0
            ELSE
               total_spin = all_spins(1, 3) - all_spins(2, 3)
            ENDIF
            IF (all_nspins(1) .EQ. 1) THEN
               sub_spin_1 = 0
            ELSE
               sub_spin_1 = all_spins(1, 1) - all_spins(2, 1)
            ENDIF
            IF (all_nspins(2) .EQ. 1) THEN
               sub_spin_2 = 0
            ELSE
               sub_spin_2 = all_spins(1, 2) - all_spins(2, 2)
            ENDIF
            IF ((sub_spin_1 + sub_spin_2) .EQ. total_spin) THEN
               change_spin = .FALSE.
            ELSE
               IF (ABS(sub_spin_1 - sub_spin_2) .EQ. total_spin) THEN
                  change_spin = .TRUE.
               ELSE
                  CPABORT("Spin states of subsystems are not compatible.")
               ENDIF
            ENDIF

         ENDIF ! not open_shell
      ELSE
         CPABORT("Reference subsystem must be the third FORCE_EVAL.")
      ENDIF

   END SUBROUTINE understand_spin_states

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param embed_pot ...
!> \param add_const_pot ...
!> \param Fermi_Amaldi ...
!> \param const_pot ...
!> \param open_shell_embed ...
!> \param spin_embed_pot ...
!> \param pot_diff ...
!> \param Coulomb_guess ...
!> \param grid_opt ...
! **************************************************************************************************
   SUBROUTINE init_embed_pot(qs_env, embed_pot, add_const_pot, Fermi_Amaldi, const_pot, open_shell_embed, &
                             spin_embed_pot, pot_diff, Coulomb_guess, grid_opt)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: embed_pot
      LOGICAL                                            :: add_const_pot, Fermi_Amaldi
      TYPE(pw_p_type), POINTER                           :: const_pot
      LOGICAL                                            :: open_shell_embed
      TYPE(pw_p_type), POINTER                           :: spin_embed_pot, pot_diff
      LOGICAL                                            :: Coulomb_guess, grid_opt

      INTEGER                                            :: nelectrons
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      REAL(KIND=dp)                                      :: factor
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: v_hartree_r_space
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      ! Extract  plane waves environment
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, &
                      nelectron_spin=nelectron_spin, &
                      v_hartree_rspace=v_hartree_r_space%pw)

      ! Prepare plane-waves pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! Create embedding potential and set to zero
      NULLIFY (embed_pot)
      ALLOCATE (embed_pot)
      NULLIFY (embed_pot%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, embed_pot%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_zero(embed_pot%pw)

      ! Spin embedding potential if asked
      IF (open_shell_embed) THEN
         NULLIFY (spin_embed_pot)
         ALLOCATE (spin_embed_pot)
         NULLIFY (spin_embed_pot%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, spin_embed_pot%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(spin_embed_pot%pw)
      ENDIF

      ! Coulomb guess/constant potential
      IF (Coulomb_guess) THEN
         NULLIFY (pot_diff)
         ALLOCATE (pot_diff)
         NULLIFY (pot_diff%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, pot_diff%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(pot_diff%pw)
      ENDIF

      ! Initialize constant part of the embedding potential
      IF (add_const_pot .AND. (.NOT. grid_opt)) THEN
         ! Now the constant potential is the Coulomb one
         NULLIFY (const_pot)
         ALLOCATE (const_pot)
         NULLIFY (const_pot%pw)

         CALL pw_pool_create_pw(auxbas_pw_pool, const_pot%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(const_pot%pw)
      ENDIF

      ! Add Fermi-Amaldi potential if requested
      IF (Fermi_Amaldi) THEN

         ! Extract  Hartree potential
         NULLIFY (v_hartree_r_space%pw)
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, &
                         v_hartree_rspace=v_hartree_r_space%pw)
         CALL pw_copy(v_hartree_r_space%pw, embed_pot%pw)

         ! Calculate the number of electrons
         nelectrons = nelectron_spin(1) + nelectron_spin(2)
         factor = (REAL(nelectrons) - 1.0_dp)/(REAL(nelectrons))

         ! Scale the Hartree potential to get Fermi-Amaldi
         CALL pw_scale(embed_pot%pw, a=factor)

         ! Copy Fermi-Amaldi to embedding potential for basis-based optimization
         IF (.NOT. grid_opt) CALL pw_copy(embed_pot%pw, embed_pot%pw)

      ENDIF

   END SUBROUTINE init_embed_pot

! **************************************************************************************************
!> \brief Creates and allocates objects for optimization of embedding potential
!> \param qs_env ...
!> \param opt_embed ...
!> \param opt_embed_section ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE prepare_embed_opt(qs_env, opt_embed, opt_embed_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(opt_embed_pot_type)                           :: opt_embed
      TYPE(section_vals_type), POINTER                   :: opt_embed_section

      INTEGER                                            :: diff_size, i_dens, size_prev_dens
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      !TYPE(pw_env_type), POINTER                         :: pw_env
      !TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      ! First, read the input

      CALL read_opt_embed_section(opt_embed, opt_embed_section)

      ! All these are needed for optimization in a finite Guassian basis
      IF (.NOT. opt_embed%grid_opt) THEN
         ! Create blacs environment
         CALL get_qs_env(qs_env=qs_env, &
                         para_env=para_env)
         NULLIFY (blacs_env)
         CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

         ! Reveal the dimention of the RI basis
         CALL find_aux_dimen(qs_env, opt_embed%dimen_aux)

         ! Prepare the object for integrals
         CALL make_lri_object(qs_env, opt_embed%lri)

         ! In case if spin embedding potential has to be optimized,
         ! the dimension of variational space is two times larger
         IF (opt_embed%open_shell_embed) THEN
            opt_embed%dimen_var_aux = 2*opt_embed%dimen_aux
         ELSE
            opt_embed%dimen_var_aux = opt_embed%dimen_aux
         ENDIF

         ! Allocate expansion coefficients and gradient
         NULLIFY (opt_embed%embed_pot_grad, opt_embed%embed_pot_coef, opt_embed%step, fm_struct)

         NULLIFY (opt_embed%prev_embed_pot_grad, opt_embed%prev_embed_pot_coef, opt_embed%prev_step)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=opt_embed%dimen_var_aux, ncol_global=1)
         CALL cp_fm_create(opt_embed%embed_pot_grad, fm_struct, name="pot_grad")
         CALL cp_fm_create(opt_embed%embed_pot_coef, fm_struct, name="pot_coef")
         CALL cp_fm_create(opt_embed%prev_embed_pot_grad, fm_struct, name="prev_pot_grad")
         CALL cp_fm_create(opt_embed%prev_embed_pot_coef, fm_struct, name="prev_pot_coef")
         CALL cp_fm_create(opt_embed%step, fm_struct, name="step")

         CALL cp_fm_create(opt_embed%prev_step, fm_struct, name="prev_step")

         CALL cp_fm_struct_release(fm_struct)
         CALL cp_fm_set_all(opt_embed%embed_pot_grad, 0.0_dp)
         CALL cp_fm_set_all(opt_embed%prev_embed_pot_grad, 0.0_dp)
         CALL cp_fm_set_all(opt_embed%embed_pot_coef, 0.0_dp)
         CALL cp_fm_set_all(opt_embed%prev_embed_pot_coef, 0.0_dp)
         CALL cp_fm_set_all(opt_embed%step, 0.0_dp)

         CALL cp_fm_set_all(opt_embed%prev_step, 0.0_dp)

         ! Allocate Hessian
         NULLIFY (opt_embed%embed_pot_hess, opt_embed%prev_embed_pot_hess, fm_struct)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=opt_embed%dimen_var_aux, ncol_global=opt_embed%dimen_var_aux)
         CALL cp_fm_create(opt_embed%embed_pot_hess, fm_struct, name="pot_Hess")
         CALL cp_fm_create(opt_embed%prev_embed_pot_hess, fm_struct, name="prev_pot_Hess")
         CALL cp_fm_struct_release(fm_struct)

         ! Special structure for the kinetic energy matrix
         NULLIFY (fm_struct, opt_embed%kinetic_mat)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=opt_embed%dimen_aux, ncol_global=opt_embed%dimen_aux)
         CALL cp_fm_create(opt_embed%kinetic_mat, fm_struct, name="kinetic_mat")
         CALL cp_fm_struct_release(fm_struct)
         CALL cp_fm_set_all(opt_embed%kinetic_mat, 0.0_dp)

         ! Hessian is set as a unit matrix
         CALL cp_fm_set_all(opt_embed%embed_pot_hess, 0.0_dp, -1.0_dp)
         CALL cp_fm_set_all(opt_embed%prev_embed_pot_hess, 0.0_dp, -1.0_dp)

         ! Release blacs environment
         CALL cp_blacs_env_release(blacs_env)

         ! ELSE ! Grid-based optimization
         !my_spins = 1
         !IF (opt_embed%open_shell_embed) my_spins = 2
         !CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
         !CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         !NULLIFY (opt_embed%prev_grid_grad, opt_embed%prev_grid_pot)
         !ALLOCATE (opt_embed%prev_grid_grad(my_spins))
         !ALLOCATE (opt_embed%prev_grid_pot(my_spins))
         !DO i_spin = 1, my_spins
         !   ALLOCATE (opt_embed%prev_grid_grad(i_spin))
         !   CALL pw_pool_create_pw(auxbas_pw_pool, opt_embed%prev_grid_grad(i_spin)%pw, &
         !                          use_data=REALDATA3D, &
         !                          in_space=REALSPACE)
         !   CALL pw_zero(opt_embed%prev_grid_grad(i_spin)%pw)
         !   ALLOCATE (opt_embed%prev_grid_pot(i_spin))
         !   CALL pw_pool_create_pw(auxbas_pw_pool, opt_embed%prev_grid_pot(i_spin)%pw, &
         !                          use_data=REALDATA3D, &
         !                          in_space=REALSPACE)
         !   CALL pw_zero(opt_embed%prev_grid_pot(i_spin)%pw)
         !ENDDO
      ENDIF

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      NULLIFY (opt_embed%prev_subsys_dens)
      size_prev_dens = SUM(opt_embed%all_nspins(1:(SIZE(opt_embed%all_nspins) - 1)))
      ALLOCATE (opt_embed%prev_subsys_dens(size_prev_dens))
      DO i_dens = 1, size_prev_dens
         CALL pw_pool_create_pw(auxbas_pw_pool, opt_embed%prev_subsys_dens(i_dens)%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(opt_embed%prev_subsys_dens(i_dens)%pw)
      ENDDO
      ALLOCATE (opt_embed%max_subsys_dens_diff(size_prev_dens))

      ! Array to store functional values
      ALLOCATE (opt_embed%w_func(opt_embed%n_iter))
      opt_embed%w_func = 0.0_dp

      ! Allocate max_diff and int_diff
      diff_size = 1
      IF (opt_embed%open_shell_embed) diff_size = 2
      ALLOCATE (opt_embed%max_diff(diff_size))
      ALLOCATE (opt_embed%int_diff(diff_size))
      ALLOCATE (opt_embed%int_diff_square(diff_size))

      ! FAB update
      IF (opt_embed%fab) THEN
         NULLIFY (opt_embed%prev_embed_pot)
         ALLOCATE (opt_embed%prev_embed_pot)
         NULLIFY (opt_embed%prev_embed_pot%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, opt_embed%prev_embed_pot%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(opt_embed%prev_embed_pot%pw)
         IF (opt_embed%open_shell_embed) THEN
            NULLIFY (opt_embed%prev_spin_embed_pot)
            ALLOCATE (opt_embed%prev_spin_embed_pot)
            NULLIFY (opt_embed%prev_spin_embed_pot%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, opt_embed%prev_spin_embed_pot%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(opt_embed%prev_spin_embed_pot%pw)
         ENDIF
      ENDIF

      ! Set allowed energy decrease parameter
      opt_embed%allowed_decrease = 0.0001_dp

      ! Regularization contribution is set to zero
      opt_embed%reg_term = 0.0_dp

      ! Step is accepted in the beginning
      opt_embed%accept_step = .TRUE.
      opt_embed%newton_step = .FALSE.
      opt_embed%last_accepted = 1

      ! Set maximum and minimum trust radii
      opt_embed%max_trad = opt_embed%trust_rad*7.900_dp
      opt_embed%min_trad = opt_embed%trust_rad*0.125*0.065_dp

   END SUBROUTINE prepare_embed_opt

! **************************************************************************************************
!> \brief ...
!> \param opt_embed ...
!> \param opt_embed_section ...
! **************************************************************************************************
   SUBROUTINE read_opt_embed_section(opt_embed, opt_embed_section)
      TYPE(opt_embed_pot_type)                           :: opt_embed
      TYPE(section_vals_type), POINTER                   :: opt_embed_section

      INTEGER                                            :: embed_guess, embed_optimizer

      ! Read keywords
      CALL section_vals_val_get(opt_embed_section, "REG_LAMBDA", &
                                r_val=opt_embed%lambda)

      CALL section_vals_val_get(opt_embed_section, "N_ITER", &
                                i_val=opt_embed%n_iter)

      CALL section_vals_val_get(opt_embed_section, "TRUST_RAD", &
                                r_val=opt_embed%trust_rad)

      CALL section_vals_val_get(opt_embed_section, "DENS_CONV_MAX", &
                                r_val=opt_embed%conv_max)

      CALL section_vals_val_get(opt_embed_section, "DENS_CONV_INT", &
                                r_val=opt_embed%conv_int)

      CALL section_vals_val_get(opt_embed_section, "SPIN_DENS_CONV_MAX", &
                                r_val=opt_embed%conv_max_spin)

      CALL section_vals_val_get(opt_embed_section, "SPIN_DENS_CONV_INT", &
                                r_val=opt_embed%conv_int_spin)

      CALL section_vals_val_get(opt_embed_section, "CHARGE_DISTR_WIDTH", &
                                r_val=opt_embed%eta)

      CALL section_vals_val_get(opt_embed_section, "READ_EMBED_POT", &
                                l_val=opt_embed%read_embed_pot)

      CALL section_vals_val_get(opt_embed_section, "READ_EMBED_POT_CUBE", &
                                l_val=opt_embed%read_embed_pot_cube)

      CALL section_vals_val_get(opt_embed_section, "GRID_OPT", &
                                l_val=opt_embed%grid_opt)

      CALL section_vals_val_get(opt_embed_section, "LEEUWEN-BAERENDS", &
                                l_val=opt_embed%leeuwen)

      CALL section_vals_val_get(opt_embed_section, "FAB", &
                                l_val=opt_embed%fab)

      CALL section_vals_val_get(opt_embed_section, "VW_CUTOFF", &
                                r_val=opt_embed%vw_cutoff)

      CALL section_vals_val_get(opt_embed_section, "VW_SMOOTH_CUT_RANGE", &
                                r_val=opt_embed%vw_smooth_cutoff_range)

      CALL section_vals_val_get(opt_embed_section, "OPTIMIZER", i_val=embed_optimizer)
      SELECT CASE (embed_optimizer)
      CASE (embed_steep_desc)
         opt_embed%steep_desc = .TRUE.
      CASE (embed_quasi_newton)
         opt_embed%steep_desc = .FALSE.
         opt_embed%level_shift = .FALSE.
      CASE (embed_level_shift)
         opt_embed%steep_desc = .FALSE.
         opt_embed%level_shift = .TRUE.
      CASE DEFAULT
         opt_embed%steep_desc = .TRUE.
      END SELECT

      CALL section_vals_val_get(opt_embed_section, "POT_GUESS", i_val=embed_guess)
      SELECT CASE (embed_guess)
      CASE (embed_none)
         opt_embed%add_const_pot = .FALSE.
         opt_embed%Fermi_Amaldi = .FALSE.
         opt_embed%Coulomb_guess = .FALSE.
         opt_embed%diff_guess = .FALSE.
      CASE (embed_diff)
         opt_embed%add_const_pot = .TRUE.
         opt_embed%Fermi_Amaldi = .FALSE.
         opt_embed%Coulomb_guess = .FALSE.
         opt_embed%diff_guess = .TRUE.
      CASE (embed_fa)
         opt_embed%add_const_pot = .TRUE.
         opt_embed%Fermi_Amaldi = .TRUE.
         opt_embed%Coulomb_guess = .FALSE.
         opt_embed%diff_guess = .FALSE.
      CASE (embed_resp)
         opt_embed%add_const_pot = .TRUE.
         opt_embed%Fermi_Amaldi = .TRUE.
         opt_embed%Coulomb_guess = .TRUE.
         opt_embed%diff_guess = .FALSE.
      CASE DEFAULT
         opt_embed%add_const_pot = .FALSE.
         opt_embed%Fermi_Amaldi = .FALSE.
         opt_embed%Coulomb_guess = .FALSE.
         opt_embed%diff_guess = .FALSE.
      END SELECT

   END SUBROUTINE read_opt_embed_section

! **************************************************************************************************
!> \brief Find the dimension of the auxiliary basis for the expansion of the embedding potential
!> \param qs_env ...
!> \param dimen_aux ...
! **************************************************************************************************
   SUBROUTINE find_aux_dimen(qs_env, dimen_aux)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: dimen_aux

      INTEGER                                            :: iatom, ikind, natom, nsgf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      ! First, reveal the dimention of the RI basis
      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      natom = SIZE(particle_set)
      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      dimen_aux = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type="RI_AUX")
         dimen_aux = dimen_aux + nsgf
      END DO

      DEALLOCATE (kind_of)

   END SUBROUTINE find_aux_dimen

! **************************************************************************************************
!> \brief Prepare the lri_kind_type object for integrals between density and aux. basis functions
!> \param qs_env ...
!> \param lri ...
! **************************************************************************************************
   SUBROUTINE make_lri_object(qs_env, lri)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri

      INTEGER                                            :: ikind, natom, nkind, nsgf
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (atomic_kind, lri)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set)
      nkind = SIZE(atomic_kind_set)

      ALLOCATE (lri(nkind))
      ! Here we need only v_int and acoef (the latter as dummies)
      DO ikind = 1, nkind
         NULLIFY (lri(ikind)%acoef)
         NULLIFY (lri(ikind)%v_int)
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom)
         CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type="RI_AUX")
         ALLOCATE (lri(ikind)%acoef(natom, nsgf))
         lri(ikind)%acoef = 0._dp
         ALLOCATE (lri(ikind)%v_int(natom, nsgf))
         lri(ikind)%v_int = 0._dp
      END DO

   END SUBROUTINE make_lri_object

! **************************************************************************************************
!> \brief Read the external embedding potential, not to be optimized
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE given_embed_pot(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      LOGICAL                                            :: open_shell_embed
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), POINTER                           :: embed_pot, spin_embed_pot
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool_subsys
      TYPE(section_vals_type), POINTER                   :: input, qs_section

      qs_env%given_embed_pot = .TRUE.
      NULLIFY (input, dft_control, embed_pot, spin_embed_pot, embed_pot, spin_embed_pot, &
               qs_section)
      CALL get_qs_env(qs_env=qs_env, &
                      input=input, &
                      dft_control=dft_control, &
                      pw_env=pw_env)
      qs_section => section_vals_get_subs_vals(input, "DFT%QS")
      open_shell_embed = .FALSE.
      IF (dft_control%nspins .EQ. 2) open_shell_embed = .TRUE.

      ! Prepare plane-waves pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool_subsys)

      ! Create embedding potential
      !CALL get_qs_env(qs_env=qs_env, &
      !                embed_pot=embed_pot)
      ALLOCATE (embed_pot)
      NULLIFY (embed_pot%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool_subsys, embed_pot%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      IF (open_shell_embed) THEN
         ! Create spin embedding potential
         !  CALL get_qs_env(qs_env=qs_env, &
         !               spin_embed_pot=spin_embed_pot)
         ALLOCATE (spin_embed_pot)
         NULLIFY (spin_embed_pot%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool_subsys, spin_embed_pot%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
      ENDIF
      ! Read the cubes
      CALL read_embed_pot_cube(embed_pot, spin_embed_pot, qs_section, open_shell_embed)

      IF (.NOT. open_shell_embed) THEN
         CALL set_qs_env(qs_env=qs_env, embed_pot=embed_pot)
      ELSE
         CALL set_qs_env(qs_env=qs_env, embed_pot=embed_pot, spin_embed_pot=spin_embed_pot)
      ENDIF

   END SUBROUTINE given_embed_pot

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param embed_pot ...
!> \param spin_embed_pot ...
!> \param section ...
!> \param opt_embed ...
! **************************************************************************************************
   SUBROUTINE read_embed_pot(qs_env, embed_pot, spin_embed_pot, section, opt_embed)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: embed_pot, spin_embed_pot
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(opt_embed_pot_type)                           :: opt_embed

      ! Read the potential as a vector in the auxiliary basis
      IF (opt_embed%read_embed_pot) &
         CALL read_embed_pot_vector(qs_env, embed_pot, spin_embed_pot, section, opt_embed%embed_pot_coef, &
                                    opt_embed%open_shell_embed)
      ! Read the potential as a cube (two cubes for open shell)
      IF (opt_embed%read_embed_pot_cube) &
         CALL read_embed_pot_cube(embed_pot, spin_embed_pot, section, opt_embed%open_shell_embed)

   END SUBROUTINE read_embed_pot

! **************************************************************************************************
!> \brief ...
!> \param embed_pot ...
!> \param spin_embed_pot ...
!> \param section ...
!> \param open_shell_embed ...
! **************************************************************************************************
   SUBROUTINE read_embed_pot_cube(embed_pot, spin_embed_pot, section, open_shell_embed)
      TYPE(pw_p_type), POINTER                           :: embed_pot, spin_embed_pot
      TYPE(section_vals_type), POINTER                   :: section
      LOGICAL                                            :: open_shell_embed

      CHARACTER(LEN=default_path_length)                 :: filename
      LOGICAL                                            :: exist
      REAL(KIND=dp)                                      :: scaling_factor

      exist = .FALSE.
      CALL section_vals_val_get(section, "EMBED_CUBE_FILE_NAME", c_val=filename)
      INQUIRE (FILE=filename, exist=exist)
      IF (.NOT. exist) &
         CPABORT("Embedding cube file not found. ")

      scaling_factor = 1.0_dp
      CALL cp_cube_to_pw(embed_pot%pw, filename, scaling_factor)

      ! Spin-dependent part of the potential
      IF (open_shell_embed) THEN
         exist = .FALSE.
         CALL section_vals_val_get(section, "EMBED_SPIN_CUBE_FILE_NAME", c_val=filename)
         INQUIRE (FILE=filename, exist=exist)
         IF (.NOT. exist) &
            CPABORT("Embedding spin cube file not found. ")

         scaling_factor = 1.0_dp
         CALL cp_cube_to_pw(spin_embed_pot%pw, filename, scaling_factor)
      ENDIF

   END SUBROUTINE read_embed_pot_cube

! **************************************************************************************************
!> \brief Read the embedding potential from the binary file
!> \param qs_env ...
!> \param embed_pot ...
!> \param spin_embed_pot ...
!> \param section ...
!> \param embed_pot_coef ...
!> \param open_shell_embed ...
! **************************************************************************************************
   SUBROUTINE read_embed_pot_vector(qs_env, embed_pot, spin_embed_pot, section, embed_pot_coef, &
                                    open_shell_embed)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: embed_pot, spin_embed_pot
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: embed_pot_coef
      LOGICAL                                            :: open_shell_embed

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: dimen_aux, dimen_restart_basis, &
                                                            dimen_var_aux, l_global, LLL, &
                                                            nrow_local, restart_unit
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coef, coef_read
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: my_embed_pot_coef
      TYPE(cp_para_env_type), POINTER                    :: para_env

      ! Get the vector dimension
      CALL find_aux_dimen(qs_env, dimen_aux)
      IF (open_shell_embed) THEN
         dimen_var_aux = dimen_aux*2
      ELSE
         dimen_var_aux = dimen_aux
      ENDIF

      ! We need a temporary vector of coefficients
      CALL get_qs_env(qs_env=qs_env, &
                      para_env=para_env)
      NULLIFY (blacs_env)
      NULLIFY (my_embed_pot_coef, fm_struct)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)
      CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen_var_aux, ncol_global=1)
      CALL cp_fm_create(my_embed_pot_coef, fm_struct, name="my_pot_coef")
      IF (.NOT. (PRESENT(embed_pot_coef))) THEN
         NULLIFY (embed_pot_coef)
         CALL cp_fm_create(my_embed_pot_coef, fm_struct, name="pot_coef")
      ENDIF

      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(my_embed_pot_coef, 0.0_dp)

      ! Read the coefficients vector
      restart_unit = -1

      ! Allocate the attay to read the coefficients
      ALLOCATE (coef(dimen_var_aux))
      coef = 0.0_dp

      IF (para_env%ionode) THEN

         ! Get the restart file name
         CALL embed_restart_file_name(filename, section)

         CALL open_file(file_name=filename, &
                        file_action="READ", &
                        file_form="UNFORMATTED", &
                        file_status="OLD", &
                        unit_number=restart_unit)

         READ (restart_unit) dimen_restart_basis
         ! Check the dimensions of the bases: the actual and the restart one
         IF (.NOT. (dimen_restart_basis == dimen_aux)) &
            CPABORT("Wrong dimension of the embedding basis in the restart file.")

         ALLOCATE (coef_read(dimen_var_aux))
         coef_read = 0.0_dp

         READ (restart_unit) coef_read
         coef(:) = coef_read(:)
         DEALLOCATE (coef_read)

         ! Close restart file
         CALL close_file(unit_number=restart_unit)

      ENDIF

      ! Broadcast the coefficients on all processes
      CALL mp_bcast(coef, para_env%source, para_env%group)

      ! Copy to fm_type structure
      ! Information about full matrix gradient
      CALL cp_fm_get_info(matrix=my_embed_pot_coef, &
                          nrow_local=nrow_local, &
                          row_indices=row_indices)

      DO LLL = 1, nrow_local
         l_global = row_indices(LLL)
         my_embed_pot_coef%local_data(LLL, 1) = coef(l_global)
      ENDDO

      DEALLOCATE (coef)

      ! Copy to the my_embed_pot_coef to embed_pot_coef
      CALL cp_fm_copy_general(my_embed_pot_coef, embed_pot_coef, para_env)

      ! Build the embedding potential on the grid
      CALL update_embed_pot(embed_pot_coef, dimen_aux, embed_pot, spin_embed_pot, &
                            qs_env, .FALSE., open_shell_embed)

      ! Release my_embed_pot_coef
      CALL cp_fm_release(my_embed_pot_coef)
      IF (.NOT. (PRESENT(embed_pot_coef))) CALL cp_fm_release(embed_pot_coef)

      ! Release blacs environment
      CALL cp_blacs_env_release(blacs_env)

   END SUBROUTINE read_embed_pot_vector

! **************************************************************************************************
!> \brief Find the embedding restart file name
!> \param filename ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE embed_restart_file_name(filename, section)
      CHARACTER(LEN=default_path_length), INTENT(OUT)    :: filename
      TYPE(section_vals_type), POINTER                   :: section

      LOGICAL                                            :: exist

      exist = .FALSE.
      CALL section_vals_val_get(section, "EMBED_RESTART_FILE_NAME", c_val=filename)
      INQUIRE (FILE=filename, exist=exist)
      IF (.NOT. exist) &
         CPABORT("Embedding restart file not found. ")

   END SUBROUTINE embed_restart_file_name

! **************************************************************************************************
!> \brief Deallocate stuff for optimizing embedding potential
!> \param opt_embed ...
! **************************************************************************************************
   SUBROUTINE release_opt_embed(opt_embed)
      TYPE(opt_embed_pot_type)                           :: opt_embed

      INTEGER                                            :: i_dens, i_spin, ikind

      IF (.NOT. opt_embed%grid_opt) THEN
         CALL cp_fm_release(opt_embed%embed_pot_grad)
         CALL cp_fm_release(opt_embed%embed_pot_coef)
         CALL cp_fm_release(opt_embed%step)

         CALL cp_fm_release(opt_embed%prev_step)

         CALL cp_fm_release(opt_embed%embed_pot_hess)
         CALL cp_fm_release(opt_embed%prev_embed_pot_grad)
         CALL cp_fm_release(opt_embed%prev_embed_pot_coef)
         CALL cp_fm_release(opt_embed%prev_embed_pot_hess)
         CALL cp_fm_release(opt_embed%kinetic_mat)
         DEALLOCATE (opt_embed%w_func)
         DEALLOCATE (opt_embed%max_diff)
         DEALLOCATE (opt_embed%int_diff)

         DO ikind = 1, SIZE(opt_embed%lri)
            DEALLOCATE (opt_embed%lri(ikind)%v_int)
            DEALLOCATE (opt_embed%lri(ikind)%acoef)
         ENDDO
         DEALLOCATE (opt_embed%lri)
      ENDIF

      DO i_dens = 1, SIZE(opt_embed%prev_subsys_dens)
         CALL pw_release(opt_embed%prev_subsys_dens(i_dens)%pw)
      ENDDO
      DEALLOCATE (opt_embed%prev_subsys_dens)
      DEALLOCATE (opt_embed%max_subsys_dens_diff)

      DEALLOCATE (opt_embed%all_nspins)

      IF (opt_embed%add_const_pot .AND. (.NOT. opt_embed%grid_opt)) THEN
         CALL pw_release(opt_embed%const_pot%pw)
         DEALLOCATE (opt_embed%const_pot)
      ENDIF

      IF (opt_embed%Coulomb_guess) THEN
         CALL pw_release(opt_embed%pot_diff%pw)
         DEALLOCATE (opt_embed%pot_diff)
      ENDIF

      IF (opt_embed%fab) THEN
         CALL pw_release(opt_embed%prev_embed_pot%pw)
         DEALLOCATE (opt_embed%prev_embed_pot)
         IF (opt_embed%open_shell_embed) THEN
            CALL pw_release(opt_embed%prev_spin_embed_pot%pw)
            DEALLOCATE (opt_embed%prev_spin_embed_pot)
         ENDIF
         DO i_spin = 1, SIZE(opt_embed%v_w)
            CALL pw_release(opt_embed%v_w(i_spin)%pw)
         ENDDO
         DEALLOCATE (opt_embed%v_w)
      ENDIF

   END SUBROUTINE release_opt_embed

! **************************************************************************************************
!> \brief Calculates subsystem Coulomb potential from the RESP charges of the total system
!> \param v_rspace ...
!> \param rhs ...
!> \param mapping_section ...
!> \param qs_env ...
!> \param nforce_eval ...
!> \param iforce_eval ...
!> \param eta ...
! **************************************************************************************************
   SUBROUTINE Coulomb_guess(v_rspace, rhs, mapping_section, qs_env, nforce_eval, iforce_eval, eta)
      TYPE(pw_p_type)                                    :: v_rspace
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rhs
      TYPE(section_vals_type), POINTER                   :: mapping_section
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: nforce_eval, iforce_eval
      REAL(KIND=dp)                                      :: eta

      INTEGER                                            :: iparticle, jparticle, natom
      INTEGER, DIMENSION(:), POINTER                     :: map_index
      REAL(KIND=dp)                                      :: dvol, normalize_factor
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rhs_subsys
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_resp, v_resp_gspace, v_resp_rspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_subsys_type), POINTER                      :: subsys

      ! Get available particles
      NULLIFY (subsys)
      CALL get_qs_env(qs_env=qs_env, subsys=subsys, pw_env=pw_env)
      CALL qs_subsys_get(subsys, particles=particles)
      natom = particles%n_els

      ALLOCATE (rhs_subsys(natom))

      NULLIFY (map_index)
      CALL get_subsys_map_index(mapping_section, natom, iforce_eval, nforce_eval, &
                                map_index, .TRUE.)

      ! Mapping particles from iforce_eval environment to the embed env
      DO iparticle = 1, natom
         jparticle = map_index(iparticle)
         rhs_subsys(iparticle) = rhs(jparticle)
      END DO

      ! Prepare plane waves
      NULLIFY (auxbas_pw_pool)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             v_resp_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             v_resp_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      CALL pw_pool_create_pw(auxbas_pw_pool, rho_resp%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ! Calculate charge density
      CALL pw_zero(rho_resp%pw)
      CALL calculate_rho_resp_all(rho_resp, rhs_subsys, natom, eta, qs_env)

      ! Calculate potential
      CALL pw_zero(v_resp_gspace%pw)
      CALL pw_poisson_solve(poisson_env, rho_resp%pw, &
                            vhartree=v_resp_gspace%pw)
      CALL pw_zero(v_resp_rspace%pw)
      CALL pw_transfer(v_resp_gspace%pw, v_resp_rspace%pw)
      dvol = v_resp_rspace%pw%pw_grid%dvol
      CALL pw_scale(v_resp_rspace%pw, dvol)
      normalize_factor = SQRT((eta/pi)**3)
      !normalize_factor = -2.0_dp
      CALL pw_scale(v_resp_rspace%pw, normalize_factor)

      ! Hard copy potential
      v_rspace%pw%cr3d(:, :, :) = v_resp_rspace%pw%cr3d(:, :, :)

      ! Release plane waves
      CALL pw_release(v_resp_gspace%pw)
      CALL pw_release(v_resp_rspace%pw)
      CALL pw_release(rho_resp%pw)

      ! Deallocate map_index array
      DEALLOCATE (map_index)
      ! Deallocate charges
      DEALLOCATE (rhs_subsys)

   END SUBROUTINE Coulomb_guess

! **************************************************************************************************
!> \brief Creates a subsystem embedding potential
!> \param qs_env ...
!> \param embed_pot ...
!> \param embed_pot_subsys ...
!> \param spin_embed_pot ...
!> \param spin_embed_pot_subsys ...
!> \param open_shell_embed ...
!> \param change_spin_sign ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE make_subsys_embed_pot(qs_env, embed_pot, embed_pot_subsys, &
                                    spin_embed_pot, spin_embed_pot_subsys, open_shell_embed, &
                                    change_spin_sign)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: embed_pot, embed_pot_subsys, &
                                                            spin_embed_pot, spin_embed_pot_subsys
      LOGICAL                                            :: open_shell_embed, change_spin_sign

      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool_subsys

      ! Extract  plane waves environment
      CALL get_qs_env(qs_env, pw_env=pw_env)

      ! Prepare plane-waves pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool_subsys)

      ! Create embedding potential and set to zero
      NULLIFY (embed_pot_subsys)
      ALLOCATE (embed_pot_subsys)
      NULLIFY (embed_pot_subsys%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool_subsys, embed_pot_subsys%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ! Hard copy the grid
      embed_pot_subsys%pw%cr3d(:, :, :) = embed_pot%pw%cr3d(:, :, :)

      IF (open_shell_embed) THEN
         NULLIFY (spin_embed_pot_subsys)
         ALLOCATE (spin_embed_pot_subsys)
         NULLIFY (spin_embed_pot_subsys%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool_subsys, spin_embed_pot_subsys%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         ! Hard copy the grid
         IF (change_spin_sign) THEN
            spin_embed_pot_subsys%pw%cr3d(:, :, :) = -spin_embed_pot%pw%cr3d(:, :, :)
         ELSE
            spin_embed_pot_subsys%pw%cr3d(:, :, :) = spin_embed_pot%pw%cr3d(:, :, :)
         ENDIF
      ENDIF

   END SUBROUTINE make_subsys_embed_pot

! **************************************************************************************************
!> \brief Calculates the derivative of the embedding potential wrt to the expansion coefficients
!> \param qs_env ...
!> \param diff_rho_r ...
!> \param diff_rho_spin ...
!> \param opt_embed ...
!> \author Vladimir Rybkin
! **************************************************************************************************

   SUBROUTINE calculate_embed_pot_grad(qs_env, diff_rho_r, diff_rho_spin, opt_embed)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: diff_rho_r, diff_rho_spin
      TYPE(opt_embed_pot_type)                           :: opt_embed

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_embed_pot_grad', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: embed_pot_coeff_spin, &
                                                            embed_pot_coeff_spinless, &
                                                            regular_term, spin_reg, spinless_reg
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      ! We destroy the previous gradient and Hessian:
      ! current data are now previous data
      CALL cp_fm_to_fm(opt_embed%embed_pot_grad, opt_embed%prev_embed_pot_grad)
      CALL cp_fm_to_fm(opt_embed%embed_pot_Hess, opt_embed%prev_embed_pot_Hess)

      NULLIFY (pw_env)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, para_env=para_env)

      ! Get plane waves pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! Calculate potential gradient coefficients
      CALL calculate_embed_pot_grad_inner(qs_env, opt_embed%dimen_aux, diff_rho_r, diff_rho_spin, &
                                          opt_embed%embed_pot_grad, &
                                          opt_embed%open_shell_embed, opt_embed%lri)

      ! Add regularization with kinetic matrix
      IF (opt_embed%i_iter .EQ. 1) THEN ! Else it is kept in memory
         CALL compute_kinetic_mat(qs_env, opt_embed%kinetic_mat)
      ENDIF

      NULLIFY (regular_term)
      CALL cp_fm_get_info(matrix=opt_embed%embed_pot_grad, &
                          matrix_struct=fm_struct)
      CALL cp_fm_create(regular_term, fm_struct, name="regular_term")
      CALL cp_fm_set_all(regular_term, 0.0_dp)

      ! In case of open shell embedding we need two terms of dimen_aux=dimen_var_aux/2 for
      ! the spinless and the spin parts
      IF (opt_embed%open_shell_embed) THEN
         ! Prepare auxiliary full matrices
         NULLIFY (embed_pot_coeff_spinless, embed_pot_coeff_spin, spinless_reg, spin_reg, fm_struct, blacs_env)

         !CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

         CALL cp_fm_get_info(matrix=opt_embed%embed_pot_coef, context=blacs_env)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=opt_embed%dimen_aux, ncol_global=1)
         CALL cp_fm_create(embed_pot_coeff_spinless, fm_struct, name="pot_coeff_spinless")
         CALL cp_fm_create(embed_pot_coeff_spin, fm_struct, name="pot_coeff_spin")
         CALL cp_fm_create(spinless_reg, fm_struct, name="spinless_reg")
         CALL cp_fm_create(spin_reg, fm_struct, name="spin_reg")
         CALL cp_fm_set_all(embed_pot_coeff_spinless, 0.0_dp)
         CALL cp_fm_set_all(embed_pot_coeff_spin, 0.0_dp)
         CALL cp_fm_set_all(spinless_reg, 0.0_dp)
         CALL cp_fm_set_all(spin_reg, 0.0_dp)
         CALL cp_fm_struct_release(fm_struct)

         ! Copy coefficients to the auxiliary structures
         CALL cp_fm_to_fm_submat(msource=opt_embed%embed_pot_coef, &
                                 mtarget=embed_pot_coeff_spinless, &
                                 nrow=opt_embed%dimen_aux, ncol=1, &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=1)
         CALL cp_fm_to_fm_submat(msource=opt_embed%embed_pot_coef, &
                                 mtarget=embed_pot_coeff_spin, &
                                 nrow=opt_embed%dimen_aux, ncol=1, &
                                 s_firstrow=opt_embed%dimen_aux + 1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=1)
         ! Multiply
         CALL cp_gemm(transa="N", transb="N", m=opt_embed%dimen_aux, n=1, &
                      k=opt_embed%dimen_aux, alpha=1.0_dp, &
                      matrix_a=opt_embed%kinetic_mat, matrix_b=embed_pot_coeff_spinless, &
                      beta=0.0_dp, matrix_c=spinless_reg)
         CALL cp_gemm(transa="N", transb="N", m=opt_embed%dimen_aux, n=1, &
                      k=opt_embed%dimen_aux, alpha=1.0_dp, &
                      matrix_a=opt_embed%kinetic_mat, matrix_b=embed_pot_coeff_spin, &
                      beta=0.0_dp, matrix_c=spin_reg)
         ! Copy from the auxiliary structures to the full regularization term
         CALL cp_fm_to_fm_submat(msource=spinless_reg, &
                                 mtarget=regular_term, &
                                 nrow=opt_embed%dimen_aux, ncol=1, &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=1)
         CALL cp_fm_to_fm_submat(msource=spin_reg, &
                                 mtarget=regular_term, &
                                 nrow=opt_embed%dimen_aux, ncol=1, &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=opt_embed%dimen_aux + 1, t_firstcol=1)
         ! Release internally used auxiliary structures
         CALL cp_fm_release(embed_pot_coeff_spinless)
         CALL cp_fm_release(embed_pot_coeff_spin)
         CALL cp_fm_release(spin_reg)
         CALL cp_fm_release(spinless_reg)

      ELSE ! Simply multiply
         CALL cp_gemm(transa="N", transb="N", m=opt_embed%dimen_var_aux, n=1, &
                      k=opt_embed%dimen_var_aux, alpha=1.0_dp, &
                      matrix_a=opt_embed%kinetic_mat, matrix_b=opt_embed%embed_pot_coef, &
                      beta=0.0_dp, matrix_c=regular_term)
      ENDIF

      ! Scale by the regularization parameter and add to the gradient
      CALL cp_fm_scale_and_add(1.0_dp, opt_embed%embed_pot_grad, 4.0_dp*opt_embed%lambda, regular_term)

      ! Calculate the regularization contribution to the energy functional
      CALL cp_fm_trace(opt_embed%embed_pot_coef, regular_term, opt_embed%reg_term)
      opt_embed%reg_term = 2.0_dp*opt_embed%lambda*opt_embed%reg_term

      ! Deallocate regular term
      CALL cp_fm_release(regular_term)

      CALL timestop(handle)

   END SUBROUTINE calculate_embed_pot_grad

! **************************************************************************************************
!> \brief Performs integration for the embedding potential gradient
!> \param qs_env ...
!> \param dimen_aux ...
!> \param rho_r ...
!> \param rho_spin ...
!> \param embed_pot_grad ...
!> \param open_shell_embed ...
!> \param lri ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE calculate_embed_pot_grad_inner(qs_env, dimen_aux, rho_r, rho_spin, embed_pot_grad, &
                                             open_shell_embed, lri)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: dimen_aux
      TYPE(pw_p_type), POINTER                           :: rho_r, rho_spin
      TYPE(cp_fm_type), POINTER                          :: embed_pot_grad
      LOGICAL                                            :: open_shell_embed
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_embed_pot_grad_inner', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ikind, l_global, LLL, &
                                                            nrow_local, nsgf, start_pos
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: pot_grad
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

! Needed to store integrals

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      cell=cell, &
                      atomic_kind_set=atomic_kind_set, &
                      para_env=para_env)

      ! Create wf_vector and gradient
      IF (open_shell_embed) THEN
         ALLOCATE (pot_grad(dimen_aux*2))
      ELSE
         ALLOCATE (pot_grad(dimen_aux))
      ENDIF

      ! Use lri subroutine
      DO ikind = 1, SIZE(lri)
         lri(ikind)%v_int = 0.0_dp
      END DO

      CALL integrate_v_rspace_one_center(rho_r, qs_env, lri, &
                                         .FALSE., "RI_AUX")
      DO ikind = 1, SIZE(lri)
         CALL mp_sum(lri(ikind)%v_int, para_env%group)
      END DO

      pot_grad = 0.0_dp
      start_pos = 1
      DO ikind = 1, SIZE(lri)
         DO iatom = 1, SIZE(lri(ikind)%v_int, DIM=1)
            nsgf = SIZE(lri(ikind)%v_int(iatom, :))
            pot_grad(start_pos:start_pos + nsgf - 1) = lri(ikind)%v_int(iatom, :)
            start_pos = start_pos + nsgf
         ENDDO
      ENDDO

      ! Open-shell embedding
      IF (open_shell_embed) THEN
         DO ikind = 1, SIZE(lri)
            lri(ikind)%v_int = 0.0_dp
         END DO

         CALL integrate_v_rspace_one_center(rho_spin, qs_env, lri, &
                                            .FALSE., "RI_AUX")
         DO ikind = 1, SIZE(lri)
            CALL mp_sum(lri(ikind)%v_int, para_env%group)
         END DO

         start_pos = dimen_aux + 1
         DO ikind = 1, SIZE(lri)
            DO iatom = 1, SIZE(lri(ikind)%v_int, DIM=1)
               nsgf = SIZE(lri(ikind)%v_int(iatom, :))
               pot_grad(start_pos:start_pos + nsgf - 1) = lri(ikind)%v_int(iatom, :)
               start_pos = start_pos + nsgf
            ENDDO
         ENDDO
      ENDIF

      ! Scale by the cell volume
      pot_grad = pot_grad*rho_r%pw%pw_grid%dvol

      ! Information about full matrix gradient
      CALL cp_fm_get_info(matrix=embed_pot_grad, &
                          nrow_local=nrow_local, &
                          row_indices=row_indices)

      ! Copy the gradient into the full matrix
      DO LLL = 1, nrow_local
         l_global = row_indices(LLL)
         embed_pot_grad%local_data(LLL, 1) = pot_grad(l_global)
      ENDDO

      DEALLOCATE (pot_grad)

      CALL timestop(handle)

   END SUBROUTINE calculate_embed_pot_grad_inner

! **************************************************************************************************
!> \brief Calculates kinetic energy matrix in auxiliary basis in the fm format
!> \param qs_env ...
!> \param kinetic_mat ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE compute_kinetic_mat(qs_env, kinetic_mat)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: kinetic_mat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_kinetic_mat', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_t
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      NULLIFY (ks_env, sab_orb, matrix_t)

      ! First, get the dbcsr structure from the overlap matrix
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, sab_orb=sab_orb)

      ! Calculate kinetic matrix
      CALL build_kinetic_matrix(ks_env, matrix_t=matrix_t, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="RI_AUX", &
                                sab_nl=sab_orb, calculate_forces=.FALSE.)

      ! Change to the fm format
      CALL copy_dbcsr_to_fm(matrix_t(1)%matrix, kinetic_mat)

      ! Release memory
      CALL dbcsr_deallocate_matrix_set(matrix_t)

      CALL timestop(handle)

   END SUBROUTINE compute_kinetic_mat

! **************************************************************************************************
!> \brief Regularizes the Wu-Yang potential on the grid
!> \param potential ...
!> \param pw_env ...
!> \param lambda ...
!> \param reg_term ...
! **************************************************************************************************
   SUBROUTINE grid_regularize(potential, pw_env, lambda, reg_term)
      TYPE(pw_p_type), POINTER                           :: potential
      TYPE(pw_env_type), POINTER                         :: pw_env
      REAL(KIND=dp)                                      :: lambda, reg_term

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, n, ub
      TYPE(pw_p_type), DIMENSION(3)                      :: dpot, dpot_g
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), POINTER                             :: dr2_pot, grid_reg, grid_reg_g, &
                                                            potential_g, square_norm_dpot

      !
      ! First, the contribution to the gradient
      !

      ! Get some of the grids ready
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      NULLIFY (potential_g)
      CALL pw_pool_create_pw(auxbas_pw_pool, potential_g, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      NULLIFY (dr2_pot)
      CALL pw_pool_create_pw(auxbas_pw_pool, dr2_pot, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      NULLIFY (grid_reg)
      CALL pw_pool_create_pw(auxbas_pw_pool, grid_reg, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      NULLIFY (grid_reg_g)
      CALL pw_pool_create_pw(auxbas_pw_pool, grid_reg_g, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_zero(grid_reg_g)

      ! Transfer potential to the reciprocal space
      CALL pw_transfer(potential%pw, potential_g)

      ! Calculate second derivatives: dx^2, dy^2, dz^2
      DO i = 1, 3
         CALL pw_dr2(potential_g, dr2_pot, i, i)
         CALL pw_axpy(dr2_pot, grid_reg_g, 1.0_dp)
      ENDDO
      ! Transfer potential to the real space
      CALL pw_transfer(grid_reg_g, grid_reg)

      ! Update the potential with a regularization term
      CALL pw_axpy(grid_reg, potential%pw, -4.0_dp*lambda)

      !
      ! Second, the contribution to the functional
      !
      DO i = 1, 3
         NULLIFY (dpot(i)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, dpot(i)%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         NULLIFY (dpot_g(i)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, dpot_g(i)%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
      ENDDO

      NULLIFY (square_norm_dpot)
      CALL pw_pool_create_pw(auxbas_pw_pool, square_norm_dpot, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      DO i = 1, 3
         n(:) = 0
         n(i) = 1
         CALL pw_copy(potential_g, dpot_g(i)%pw)
         CALL pw_derive(dpot_g(i)%pw, n(:))
         CALL pw_transfer(dpot_g(i)%pw, dpot(i)%pw)
      END DO

      lb(1:3) = square_norm_dpot%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = square_norm_dpot%pw_grid%bounds_local(2, 1:3)
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(dpot,lb,square_norm_dpot,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               square_norm_dpot%cr3d(i, j, k) = (dpot(1)%pw%cr3d(i, j, k)* &
                                                 dpot(1)%pw%cr3d(i, j, k) + &
                                                 dpot(2)%pw%cr3d(i, j, k)* &
                                                 dpot(2)%pw%cr3d(i, j, k) + &
                                                 dpot(3)%pw%cr3d(i, j, k)* &
                                                 dpot(3)%pw%cr3d(i, j, k))
            END DO
         END DO
      END DO
!$OMP    END PARALLEL DO

      reg_term = 2*lambda*pw_integrate_function(fun=square_norm_dpot)

      ! Release
      CALL pw_pool_give_back_pw(auxbas_pw_pool, potential_g)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, dr2_pot)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, grid_reg)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, grid_reg_g)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, square_norm_dpot)
      DO i = 1, 3
         CALL pw_pool_give_back_pw(auxbas_pw_pool, dpot(i)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, dpot_g(i)%pw)
      END DO

   END SUBROUTINE grid_regularize

! **************************************************************************************************
!> \brief Takes maximization step in embedding potential optimization
!> \param diff_rho_r ...
!> \param diff_rho_spin ...
!> \param opt_embed ...
!> \param embed_pot ...
!> \param spin_embed_pot ...
!> \param rho_r_ref ...
!> \param qs_env ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE opt_embed_step(diff_rho_r, diff_rho_spin, opt_embed, embed_pot, spin_embed_pot, rho_r_ref, qs_env)
      TYPE(pw_p_type), POINTER                           :: diff_rho_r, diff_rho_spin
      TYPE(opt_embed_pot_type)                           :: opt_embed
      TYPE(pw_p_type), POINTER                           :: embed_pot, spin_embed_pot
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r_ref
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'opt_embed_step', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: thresh = 0.000001_dp

      INTEGER                                            :: handle, l_global, LLL, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: diag_grad, diag_step, fm_U, fm_U_scale
      TYPE(pw_env_type), POINTER                         :: pw_env

      CALL timeset(routineN, handle)

      IF (opt_embed%grid_opt) THEN ! Grid based optimization

         opt_embed%step_len = opt_embed%trust_rad
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
         IF (opt_embed%leeuwen) THEN
            CALL Leeuwen_Baerends_potential_update(pw_env, embed_pot, spin_embed_pot, diff_rho_r, diff_rho_spin, &
                                                   rho_r_ref, opt_embed%open_shell_embed, opt_embed%trust_rad)
         ELSE
            IF (opt_embed%fab) THEN
               CALL FAB_update(qs_env, rho_r_ref, opt_embed%prev_embed_pot, opt_embed%prev_spin_embed_pot, &
                               embed_pot, spin_embed_pot, &
                               diff_rho_r, diff_rho_spin, opt_embed%v_w, opt_embed%i_iter, opt_embed%trust_rad, &
                               opt_embed%open_shell_embed, opt_embed%vw_cutoff, opt_embed%vw_smooth_cutoff_range)
            ELSE
               CALL grid_based_step(diff_rho_r, diff_rho_spin, pw_env, opt_embed, embed_pot, spin_embed_pot)
            ENDIF
         ENDIF

      ELSE ! Finite basis optimization
         ! If the previous step has been rejected, we go back to the previous expansion coefficients
         IF (.NOT. opt_embed%accept_step) &
            CALL cp_fm_scale_and_add(1.0_dp, opt_embed%embed_pot_coef, -1.0_dp, opt_embed%step)

         ! Do a simple steepest descent
         IF (opt_embed%steep_desc) THEN
            IF (opt_embed%i_iter .GT. 2) &
               opt_embed%trust_rad = Barzilai_Borwein(opt_embed%step, opt_embed%prev_step, &
                                                      opt_embed%embed_pot_grad, opt_embed%prev_embed_pot_grad)
            IF (ABS(opt_embed%trust_rad) .GT. opt_embed%max_trad) THEN
               IF (opt_embed%trust_rad .GT. 0.0_dp) THEN
                  opt_embed%trust_rad = opt_embed%max_trad
               ELSE
                  opt_embed%trust_rad = -opt_embed%max_trad
               ENDIF
            ENDIF

            CALL cp_fm_to_fm(opt_embed%step, opt_embed%prev_step)
            CALL cp_fm_scale_and_add(0.0_dp, opt_embed%prev_step, 1.0_dp, opt_embed%step)
            CALL cp_fm_set_all(opt_embed%step, 0.0_dp)
            CALL cp_fm_scale_and_add(1.0_dp, opt_embed%step, opt_embed%trust_rad, opt_embed%embed_pot_grad)
            opt_embed%step_len = opt_embed%trust_rad
         ELSE

            ! First, update the Hessian inverse if needed
            IF (opt_embed%i_iter > 1) THEN
               IF (opt_embed%accept_step) & ! We don't update Hessian if the step has been rejected
                  CALL symm_rank_one_update(opt_embed%embed_pot_grad, opt_embed%prev_embed_pot_grad, &
                                            opt_embed%step, opt_embed%prev_embed_pot_Hess, opt_embed%embed_pot_Hess)
            ENDIF

            ! Add regularization term to the Hessian
            !CALL cp_fm_scale_and_add(1.0_dp, opt_embed%embed_pot_Hess, 4.0_dp*opt_embed%lambda, &
            !                         opt_embed%kinetic_mat)

            ! Else use the first initial Hessian. Now it's just the unit matrix: embed_pot_hess
            ! Second, invert the Hessian
            ALLOCATE (eigenval(opt_embed%dimen_var_aux))
            eigenval = 0.0_dp
            NULLIFY (fm_U, fm_U_scale, diag_grad, diag_step)
            CALL cp_fm_get_info(matrix=opt_embed%embed_pot_hess, &
                                matrix_struct=fm_struct)
            CALL cp_fm_create(fm_U, fm_struct, name="fm_U")
            CALL cp_fm_create(fm_U_scale, fm_struct, name="fm_U")
            CALL cp_fm_set_all(fm_U, 0.0_dp)
            CALL cp_fm_set_all(fm_U_scale, 0.0_dp)
            CALL cp_fm_get_info(matrix=opt_embed%embed_pot_grad, &
                                matrix_struct=fm_struct)
            CALL cp_fm_create(diag_grad, fm_struct, name="diag_grad")
            CALL cp_fm_set_all(diag_grad, 0.0_dp)
            CALL cp_fm_create(diag_step, fm_struct, name="diag_step")
            CALL cp_fm_set_all(diag_step, 0.0_dp)

            ! Store the Hessian as it will be destroyed in diagonalization: use fm_U_scal for it
            CALL cp_fm_to_fm(opt_embed%embed_pot_hess, fm_U_scale)

            ! Diagonalize Hessian
            CALL choose_eigv_solver(opt_embed%embed_pot_hess, fm_U, eigenval)

            ! Copy the Hessian back
            CALL cp_fm_to_fm(fm_U_scale, opt_embed%embed_pot_hess)

            ! Find the step in diagonal representation, begin with gradient
            CALL cp_gemm(transa="T", transb="N", m=opt_embed%dimen_var_aux, n=1, &
                         k=opt_embed%dimen_var_aux, alpha=1.0_dp, &
                         matrix_a=fm_U, matrix_b=opt_embed%embed_pot_grad, beta=0.0_dp, &
                         matrix_c=diag_grad)

            CALL cp_fm_get_info(matrix=opt_embed%embed_pot_coef, &
                                nrow_local=nrow_local, &
                                row_indices=row_indices)

            DO LLL = 1, nrow_local
               l_global = row_indices(LLL)
               IF (ABS(eigenval(l_global)) .GE. thresh) THEN
                  diag_step%local_data(LLL, 1) = &
                     -diag_grad%local_data(LLL, 1)/(eigenval(l_global))
               ELSE
                  diag_step%local_data(LLL, 1) = 0.0_dp
               ENDIF
            ENDDO
            CALL cp_fm_trace(diag_step, diag_step, opt_embed%step_len)

            ! Transform step to a non-diagonal representation
            CALL cp_gemm(transa="N", transb="N", m=opt_embed%dimen_var_aux, n=1, &
                         k=opt_embed%dimen_var_aux, alpha=1.0_dp, &
                         matrix_a=fm_U, matrix_b=diag_step, beta=0.0_dp, &
                         matrix_c=opt_embed%step)

            ! Now use fm_U_scale for scaled eigenvectors
            CALL cp_fm_to_fm(fm_U, fm_U_scale)
            CALL cp_fm_column_scale(fm_U_scale, eigenval)

            CALL cp_fm_release(fm_U_scale)

            ! Scale the step to fit within the trust radius: it it's less already,
            ! then take the Newton step
            CALL cp_fm_trace(opt_embed%step, opt_embed%step, opt_embed%step_len)
            IF (opt_embed%step_len .GT. opt_embed%trust_rad) THEN

               IF (opt_embed%level_shift) THEN
                  ! Find a level shift parameter and apply it
                  CALL level_shift(opt_embed, diag_grad, eigenval, diag_step)
               ELSE ! Just scale
                  CALL cp_fm_trace(diag_step, diag_step, opt_embed%step_len)
                  CALL cp_fm_scale(opt_embed%trust_rad/opt_embed%step_len, diag_step)
               ENDIF
               CALL cp_fm_trace(diag_step, diag_step, opt_embed%step_len)
               ! Transform step to a non-diagonal representation
               CALL cp_gemm(transa="N", transb="N", m=opt_embed%dimen_var_aux, n=1, &
                            k=opt_embed%dimen_var_aux, alpha=1.0_dp, &
                            matrix_a=fm_U, matrix_b=diag_step, beta=0.0_dp, &
                            matrix_c=opt_embed%step)
               CALL cp_fm_trace(opt_embed%step, opt_embed%step, opt_embed%step_len)

               ! Recalculate step in diagonal representation
               opt_embed%newton_step = .FALSE.
            ELSE
               opt_embed%newton_step = .TRUE.
            ENDIF

            ! Release some memory
            DEALLOCATE (eigenval)
            ! Release more memory
            CALL cp_fm_release(diag_grad)
            CALL cp_fm_release(diag_step)
            CALL cp_fm_release(fm_U)

         ENDIF ! grad_descent

         ! Update the coefficients
         CALL cp_fm_scale_and_add(1.0_dp, opt_embed%embed_pot_coef, 1.0_dp, opt_embed%step)

         ! Update the embedding potential
         CALL update_embed_pot(opt_embed%embed_pot_coef, opt_embed%dimen_aux, embed_pot, &
                               spin_embed_pot, qs_env, opt_embed%add_const_pot, &
                               opt_embed%open_shell_embed, opt_embed%const_pot)
      ENDIF ! Grid-based optimization

      CALL timestop(handle)

   END SUBROUTINE opt_embed_step

!
! **************************************************************************************************
!> \brief ...
!> \param diff_rho_r ...
!> \param diff_rho_spin ...
!> \param pw_env ...
!> \param opt_embed ...
!> \param embed_pot ...
!> \param spin_embed_pot ...
! **************************************************************************************************
   SUBROUTINE grid_based_step(diff_rho_r, diff_rho_spin, pw_env, opt_embed, embed_pot, spin_embed_pot)
      TYPE(pw_p_type), POINTER                           :: diff_rho_r, diff_rho_spin
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(opt_embed_pot_type)                           :: opt_embed
      TYPE(pw_p_type), POINTER                           :: embed_pot, spin_embed_pot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grid_based_step', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: my_reg_term

      CALL timeset(routineN, handle)

      ! Take the step for spin-free part
      CALL pw_axpy(diff_rho_r%pw, embed_pot%pw, opt_embed%step_len)
      ! Regularize
      CALL grid_regularize(embed_pot, pw_env, opt_embed%lambda, my_reg_term)
      opt_embed%reg_term = opt_embed%reg_term + my_reg_term

      IF (opt_embed%open_shell_embed) THEN
         CALL pw_axpy(diff_rho_spin%pw, spin_embed_pot%pw, opt_embed%step_len)
         CALL grid_regularize(spin_embed_pot, pw_env, opt_embed%lambda, my_reg_term)
         opt_embed%reg_term = opt_embed%reg_term + my_reg_term
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE grid_based_step

! **************************************************************************************************
!> \brief ... Adds variable part of to the embedding potential
!> \param embed_pot_coef ...
!> \param dimen_aux ...
!> \param embed_pot ...
!> \param spin_embed_pot ...
!> \param qs_env ...
!> \param add_const_pot ...
!> \param open_shell_embed ...
!> \param const_pot ...
!> \author Vladimir Rybkin
! **************************************************************************************************

   SUBROUTINE update_embed_pot(embed_pot_coef, dimen_aux, embed_pot, spin_embed_pot, &
                               qs_env, add_const_pot, open_shell_embed, const_pot)
      TYPE(cp_fm_type), POINTER                          :: embed_pot_coef
      INTEGER                                            :: dimen_aux
      TYPE(pw_p_type), POINTER                           :: embed_pot, spin_embed_pot
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: add_const_pot, open_shell_embed
      TYPE(pw_p_type), OPTIONAL, POINTER                 :: const_pot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'update_embed_pot', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, l_global, LLL, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: wf_vector
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: embed_pot_coef_spin, &
                                                            embed_pot_coef_spinless, mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: psi_L, rho_g
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)
      ! Get MO coefficients: we need only the structure, therefore don't care about the spin
      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      cell=cell, &
                      atomic_kind_set=atomic_kind_set, &
                      pw_env=pw_env, mos=mos, para_env=para_env)
      CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff)

      ! Get plane waves pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! get some of the grids ready
      NULLIFY (rho_g%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      NULLIFY (psi_L%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, psi_L%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ! Create wf_vector and auxiliary wave functions
      ALLOCATE (wf_vector(dimen_aux))
      wf_vector = 0.0_dp

      ! Create auxiliary full matrices for open-shell case
      IF (open_shell_embed) THEN
         NULLIFY (blacs_env)
         CALL cp_fm_get_info(matrix=embed_pot_coef, context=blacs_env)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=dimen_aux, ncol_global=1)
         CALL cp_fm_create(embed_pot_coef_spinless, fm_struct, name="pot_coeff_spinless")
         CALL cp_fm_create(embed_pot_coef_spin, fm_struct, name="pot_coeff_spin")
         CALL cp_fm_set_all(embed_pot_coef_spinless, 0.0_dp)
         CALL cp_fm_set_all(embed_pot_coef_spin, 0.0_dp)
         CALL cp_fm_struct_release(fm_struct)

         ! Copy coefficients to the auxiliary structures
         CALL cp_fm_to_fm_submat(embed_pot_coef, &
                                 mtarget=embed_pot_coef_spinless, &
                                 nrow=dimen_aux, ncol=1, &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=1)
         CALL cp_fm_to_fm_submat(embed_pot_coef, &
                                 mtarget=embed_pot_coef_spin, &
                                 nrow=dimen_aux, ncol=1, &
                                 s_firstrow=dimen_aux + 1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=1)

         ! Spinless potential
         CALL cp_fm_get_info(matrix=embed_pot_coef_spinless, &
                             nrow_local=nrow_local, &
                             row_indices=row_indices)

         ! Copy fm_coeff to an array
         DO LLL = 1, nrow_local
            l_global = row_indices(LLL)
            wf_vector(l_global) = embed_pot_coef_spinless%local_data(LLL, 1)
         ENDDO
         CALL mp_sum(wf_vector, para_env%group)

         ! Calculate the variable part of the embedding potential
         CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env, &
                                     basis_type="RI_AUX", &
                                     external_vector=wf_vector)
         ! Update the full embedding potential
         IF (add_const_pot) THEN
            CALL pw_copy(const_pot%pw, embed_pot%pw)
         ELSE
            CALL pw_zero(embed_pot%pw)
         ENDIF

         CALL pw_axpy(psi_L%pw, embed_pot%pw)

         ! Spin-dependent potential
         wf_vector = 0.0_dp
         CALL cp_fm_get_info(matrix=embed_pot_coef_spin, &
                             nrow_local=nrow_local, &
                             row_indices=row_indices)

         ! Copy fm_coeff to an array
         DO LLL = 1, nrow_local
            l_global = row_indices(LLL)
            wf_vector(l_global) = embed_pot_coef_spin%local_data(LLL, 1)
         ENDDO
         CALL mp_sum(wf_vector, para_env%group)

         ! Calculate the variable part of the embedding potential
         CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env, &
                                     basis_type="RI_AUX", &
                                     external_vector=wf_vector)
         ! No constant potential for spin-dependent potential
         CALL pw_zero(spin_embed_pot%pw)
         CALL pw_axpy(psi_L%pw, spin_embed_pot%pw)

      ELSE ! Closed shell

         CALL cp_fm_get_info(matrix=embed_pot_coef, &
                             nrow_local=nrow_local, &
                             row_indices=row_indices)

         ! Copy fm_coeff to an array
         DO LLL = 1, nrow_local
            l_global = row_indices(LLL)
            wf_vector(l_global) = embed_pot_coef%local_data(LLL, 1)
         ENDDO
         CALL mp_sum(wf_vector, para_env%group)

         ! Calculate the variable part of the embedding potential
         CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env)

         CALL calculate_wavefunction(mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, pw_env, &
                                     basis_type="RI_AUX", &
                                     external_vector=wf_vector)

         ! Update the full embedding potential
         IF (add_const_pot) THEN
            CALL pw_copy(const_pot%pw, embed_pot%pw)
         ELSE
            CALL pw_zero(embed_pot%pw)
         ENDIF

         CALL pw_axpy(psi_L%pw, embed_pot%pw)
      ENDIF ! Open/closed shell

      ! Deallocate memory and release objects
      DEALLOCATE (wf_vector)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, psi_L%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)

      IF (open_shell_embed) THEN
         CALL cp_fm_release(embed_pot_coef_spin)
         CALL cp_fm_release(embed_pot_coef_spinless)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE update_embed_pot

! **************************************************************************************************
!> \brief BFGS update of the inverse Hessian in the full matrix format
!> \param grad ...
!> \param prev_grad ...
!> \param step ...
!> \param prev_inv_Hess ...
!> \param inv_Hess ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE inv_Hessian_update(grad, prev_grad, step, prev_inv_Hess, inv_Hess)
      TYPE(cp_fm_type), POINTER                          :: grad, prev_grad, step, prev_inv_Hess, &
                                                            inv_Hess

      INTEGER                                            :: mat_size
      REAL(KIND=dp)                                      :: factor1, s_dot_y, y_dot_B_inv_y
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mat, fm_struct_vec
      TYPE(cp_fm_type), POINTER                          :: B_inv_y, B_inv_y_s, s_s, s_y, s_y_B_inv, &
                                                            y

      ! Recover the dimension
      CALL cp_fm_get_info(matrix=inv_Hess, &
                          nrow_global=mat_size)

      CALL cp_fm_set_all(inv_Hess, 0.0_dp)
      CALL cp_fm_to_fm(prev_inv_Hess, inv_Hess)

      ! Get full matrix structures
      NULLIFY (fm_struct_mat, fm_struct_vec)

      CALL cp_fm_get_info(matrix=prev_inv_Hess, &
                          matrix_struct=fm_struct_mat)
      CALL cp_fm_get_info(matrix=grad, &
                          matrix_struct=fm_struct_vec)

      ! Allocate intermediates
      NULLIFY (B_inv_y, y, s_s, s_y, B_inv_y_s, s_y_B_inv)
      CALL cp_fm_create(B_inv_y, fm_struct_vec, name="B_inv_y")
      CALL cp_fm_create(y, fm_struct_vec, name="y")

      CALL cp_fm_create(s_s, fm_struct_mat, name="s_s")
      CALL cp_fm_create(s_y, fm_struct_mat, name="s_y")
      CALL cp_fm_create(B_inv_y_s, fm_struct_mat, name="B_inv_y_s")
      CALL cp_fm_create(s_y_B_inv, fm_struct_mat, name="s_y_B_inv")

      CALL cp_fm_set_all(B_inv_y, 0.0_dp)
      CALL cp_fm_set_all(s_s, 0.0_dp)
      CALL cp_fm_set_all(s_y, 0.0_dp)
      CALL cp_fm_set_all(B_inv_y_s, 0.0_dp)
      CALL cp_fm_set_all(s_y_B_inv, 0.0_dp)

      ! Calculate intermediates
      ! y the is gradient difference
      CALL cp_fm_get_info(matrix=grad)
      CALL cp_fm_to_fm(grad, y)
      CALL cp_fm_scale_and_add(1.0_dp, y, -1.0_dp, prev_grad)

      ! First term
      CALL cp_gemm(transa="N", transb="N", m=mat_size, n=1, &
                   k=mat_size, alpha=1.0_dp, &
                   matrix_a=prev_inv_Hess, matrix_b=y, beta=0.0_dp, &
                   matrix_c=B_inv_y)

      CALL cp_gemm(transa="N", transb="T", m=mat_size, n=mat_size, &
                   k=1, alpha=1.0_dp, &
                   matrix_a=step, matrix_b=step, beta=0.0_dp, &
                   matrix_c=s_s)

      CALL cp_gemm(transa="N", transb="T", m=mat_size, n=mat_size, &
                   k=1, alpha=1.0_dp, &
                   matrix_a=step, matrix_b=y, beta=0.0_dp, &
                   matrix_c=s_y)

      CALL cp_fm_trace(step, y, s_dot_y)

      CALL cp_fm_trace(y, y, s_dot_y)
      CALL cp_fm_trace(step, step, s_dot_y)

      CALL cp_fm_trace(y, B_inv_y, y_dot_B_inv_y)

      factor1 = (s_dot_y + y_dot_B_inv_y)/(s_dot_y)**2

      CALL cp_fm_scale_and_add(1.0_dp, inv_Hess, factor1, s_s)

      ! Second term
      CALL cp_gemm(transa="N", transb="T", m=mat_size, n=mat_size, &
                   k=1, alpha=1.0_dp, &
                   matrix_a=B_inv_y, matrix_b=step, beta=0.0_dp, &
                   matrix_c=B_inv_y_s)

      CALL cp_gemm(transa="N", transb="N", m=mat_size, n=mat_size, &
                   k=mat_size, alpha=1.0_dp, &
                   matrix_a=s_y, matrix_b=prev_inv_Hess, beta=0.0_dp, &
                   matrix_c=s_y_B_inv)

      CALL cp_fm_scale_and_add(1.0_dp, B_inv_y_s, 1.0_dp, s_y_B_inv)

      ! Assemble the new inverse Hessian
      CALL cp_fm_scale_and_add(1.0_dp, inv_Hess, -s_dot_y, B_inv_y_s)

      ! Deallocate intermediates
      CALL cp_fm_release(y)
      CALL cp_fm_release(B_inv_y)
      CALL cp_fm_release(s_s)
      CALL cp_fm_release(s_y)
      CALL cp_fm_release(B_inv_y_s)
      CALL cp_fm_release(s_y_B_inv)

   END SUBROUTINE inv_Hessian_update

! **************************************************************************************************
!> \brief ...
!> \param grad ...
!> \param prev_grad ...
!> \param step ...
!> \param prev_Hess ...
!> \param Hess ...
! **************************************************************************************************
   SUBROUTINE Hessian_update(grad, prev_grad, step, prev_Hess, Hess)
      TYPE(cp_fm_type), POINTER                          :: grad, prev_grad, step, prev_Hess, Hess

      INTEGER                                            :: mat_size
      REAL(KIND=dp)                                      :: s_b_s, y_t_s
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mat, fm_struct_vec, &
                                                            fm_struct_vec_t
      TYPE(cp_fm_type), POINTER                          :: B_s, B_s_s_B, s_t_B, y, y_y_t
      TYPE(cp_para_env_type), POINTER                    :: para_env

      ! Recover the dimension
      CALL cp_fm_get_info(matrix=Hess, &
                          nrow_global=mat_size, para_env=para_env)

      CALL cp_fm_set_all(Hess, 0.0_dp)
      CALL cp_fm_to_fm(prev_Hess, Hess)

      ! WARNING: our Hessian must be negative-definite, whereas BFGS makes it positive-definite!
      ! Therefore, we change sign in the beginning and in the end.
      CALL cp_fm_scale(-1.0_dp, Hess)

      ! Create blacs environment
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      ! Get full matrix structures
      NULLIFY (fm_struct_mat, fm_struct_vec, fm_struct_vec_t)

      CALL cp_fm_get_info(matrix=prev_Hess, &
                          matrix_struct=fm_struct_mat)
      CALL cp_fm_get_info(matrix=grad, &
                          matrix_struct=fm_struct_vec)
      CALL cp_fm_struct_create(fm_struct_vec_t, para_env=para_env, context=blacs_env, &
                               nrow_global=1, ncol_global=mat_size)

      ! Allocate intermediates
      NULLIFY (y_y_t, B_s, s_t_B, B_s_s_B, y)
      CALL cp_fm_create(B_s, fm_struct_vec, name="B_s")
      CALL cp_fm_create(s_t_B, fm_struct_vec_t, name="s_t_B")
      CALL cp_fm_create(y, fm_struct_vec, name="y")

      CALL cp_fm_create(y_y_t, fm_struct_mat, name="y_y_t")
      CALL cp_fm_create(B_s_s_B, fm_struct_mat, name="B_s_s_B")

      CALL cp_fm_set_all(y_y_t, 0.0_dp)
      CALL cp_fm_set_all(y, 0.0_dp)
      CALL cp_fm_set_all(B_s_s_B, 0.0_dp)
      CALL cp_fm_set_all(B_s, 0.0_dp)
      CALL cp_fm_set_all(s_t_B, 0.0_dp)

      ! Release the structure created only here
      CALL cp_fm_struct_release(fm_struct_vec_t)

      ! Calculate intermediates
      ! y the is gradient difference
      CALL cp_fm_to_fm(grad, y)
      CALL cp_fm_scale_and_add(1.0_dp, y, -1.0_dp, prev_grad)

      ! First term
      CALL cp_gemm(transa="N", transb="T", m=mat_size, n=mat_size, &
                   k=1, alpha=1.0_dp, &
                   matrix_a=y, matrix_b=y, beta=0.0_dp, &
                   matrix_c=y_y_t)

      CALL cp_fm_trace(y, step, y_t_s)

      CALL cp_fm_scale_and_add(1.0_dp, Hess, (1.0_dp/y_t_s), y_y_t)

      ! Second term
      CALL cp_gemm(transa="N", transb="N", m=mat_size, n=1, &
                   k=mat_size, alpha=1.0_dp, &
                   matrix_a=Hess, matrix_b=step, beta=0.0_dp, &
                   matrix_c=B_s)

      CALL cp_fm_trace(B_s, step, s_B_s)

      CALL cp_gemm(transa="T", transb="N", m=1, n=mat_size, &
                   k=mat_size, alpha=1.0_dp, &
                   matrix_a=step, matrix_b=Hess, beta=0.0_dp, &
                   matrix_c=s_t_B)

      CALL cp_gemm(transa="N", transb="N", m=mat_size, n=mat_size, &
                   k=1, alpha=1.0_dp, &
                   matrix_a=B_s, matrix_b=s_t_B, beta=0.0_dp, &
                   matrix_c=B_s_s_B)

      CALL cp_fm_scale_and_add(1.0_dp, Hess, -(1.0_dp/s_B_s), B_s_s_B)

      ! WARNING: our Hessian must be negative-definite, whereas BFGS makes it positive-definite!
      ! Therefore, we change sign in the beginning and in the end.
      CALL cp_fm_scale(-1.0_dp, Hess)

      ! Release blacs environment
      CALL cp_blacs_env_release(blacs_env)

      ! Deallocate intermediates
      CALL cp_fm_release(y_y_t)
      CALL cp_fm_release(B_s_s_B)
      CALL cp_fm_release(B_s)
      CALL cp_fm_release(s_t_B)
      CALL cp_fm_release(y)

   END SUBROUTINE Hessian_update

! **************************************************************************************************
!> \brief ...
!> \param grad ...
!> \param prev_grad ...
!> \param step ...
!> \param prev_Hess ...
!> \param Hess ...
! **************************************************************************************************
   SUBROUTINE symm_rank_one_update(grad, prev_grad, step, prev_Hess, Hess)
      TYPE(cp_fm_type), POINTER                          :: grad, prev_grad, step, prev_Hess, Hess

      INTEGER                                            :: mat_size
      REAL(KIND=dp)                                      :: factor
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mat, fm_struct_vec
      TYPE(cp_fm_type), POINTER                          :: B_x, y, y_B_x_y_B_x

      ! Recover the dimension
      CALL cp_fm_get_info(matrix=Hess, nrow_global=mat_size)

      CALL cp_fm_set_all(Hess, 0.0_dp)
      CALL cp_fm_to_fm(prev_Hess, Hess)

      ! Get full matrix structures
      NULLIFY (fm_struct_mat, fm_struct_vec)

      CALL cp_fm_get_info(matrix=prev_Hess, &
                          matrix_struct=fm_struct_mat)
      CALL cp_fm_get_info(matrix=grad, &
                          matrix_struct=fm_struct_vec)

      ! Allocate intermediates
      NULLIFY (y, B_x, y_B_x_y_B_x)
      CALL cp_fm_create(y, fm_struct_vec, name="y")
      CALL cp_fm_create(B_x, fm_struct_vec, name="B_x")
      CALL cp_fm_create(y_B_x_y_B_x, fm_struct_mat, name="y_B_x_y_B_x")

      CALL cp_fm_set_all(y, 0.0_dp)
      CALL cp_fm_set_all(B_x, 0.0_dp)
      CALL cp_fm_set_all(y_B_x_y_B_x, 0.0_dp)

      ! Calculate intermediates
      ! y the is gradient difference
      CALL cp_fm_to_fm(grad, y)
      CALL cp_fm_scale_and_add(1.0_dp, y, -1.0_dp, prev_grad)

      CALL cp_gemm(transa="N", transb="N", m=mat_size, n=1, &
                   k=mat_size, alpha=1.0_dp, &
                   matrix_a=Hess, matrix_b=step, beta=0.0_dp, &
                   matrix_c=B_x)

      CALL cp_fm_scale_and_add(1.0_dp, y, -1.0_dp, B_x)

      CALL cp_gemm(transa="N", transb="T", m=mat_size, n=mat_size, &
                   k=1, alpha=1.0_dp, &
                   matrix_a=y, matrix_b=y, beta=0.0_dp, &
                   matrix_c=y_B_x_y_B_x)

      ! Scaling factor
      CALL cp_fm_trace(y, step, factor)

      ! Assemble the Hessian
      CALL cp_fm_scale_and_add(1.0_dp, Hess, (1.0_dp/factor), y_B_x_y_B_x)

      ! Deallocate intermediates
      CALL cp_fm_release(y)
      CALL cp_fm_release(B_x)
      CALL cp_fm_release(y_B_x_y_B_x)

   END SUBROUTINE symm_rank_one_update

! **************************************************************************************************
!> \brief Controls the step, changes the trust radius if needed in maximization of the V_emb
!> \param opt_embed ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE step_control(opt_embed)
      TYPE(opt_embed_pot_type)                           :: opt_embed

      CHARACTER(LEN=*), PARAMETER :: routineN = 'step_control', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: actual_ener_change, ener_ratio, &
                                                            lin_term, pred_ener_change, quad_term
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: H_b

      CALL timeset(routineN, handle)

      NULLIFY (H_b, fm_struct)
      CALL cp_fm_get_info(matrix=opt_embed%embed_pot_grad, &
                          matrix_struct=fm_struct)
      CALL cp_fm_create(H_b, fm_struct, name="H_b")
      CALL cp_fm_set_all(H_b, 0.0_dp)

      ! Calculate the quadratic estimate for the energy
      ! Linear term
      CALL cp_fm_trace(opt_embed%step, opt_embed%embed_pot_grad, lin_term)

      ! Quadratic term
      CALL cp_gemm(transa="N", transb="N", m=opt_embed%dimen_aux, n=1, &
                   k=opt_embed%dimen_aux, alpha=1.0_dp, &
                   matrix_a=opt_embed%embed_pot_Hess, matrix_b=opt_embed%step, &
                   beta=0.0_dp, matrix_c=H_b)
      CALL cp_fm_trace(opt_embed%step, H_b, quad_term)

      pred_ener_change = lin_term + 0.5_dp*quad_term

      ! Reveal actual energy change
      actual_ener_change = opt_embed%w_func(opt_embed%i_iter) - &
                           opt_embed%w_func(opt_embed%last_accepted)

      ener_ratio = actual_ener_change/pred_ener_change

      CALL cp_fm_release(H_b)

      IF (actual_ener_change .GT. 0.0_dp) THEN ! If energy increases
         ! We accept step
         opt_embed%accept_step = .TRUE.
         ! If energy change is larger than the predicted one, increase trust radius twice
         ! Else (between 0 and 1) leave as it is, unless Newton step has been taken and if the step is less than max
         IF ((ener_ratio .GT. 1.0_dp) .AND. (.NOT. opt_embed%newton_step) .AND. &
             (opt_embed%trust_rad .LT. opt_embed%max_trad)) &
            opt_embed%trust_rad = 2.0_dp*opt_embed%trust_rad
      ELSE ! Energy decreases
         ! If the decrease is not too large we allow this step to be taken
         ! Otherwise, the step is rejected
         IF (ABS(actual_ener_change) .GE. opt_embed%allowed_decrease) THEN
            opt_embed%accept_step = .FALSE.
         ENDIF
         ! Trust radius is decreased 4 times unless it's smaller than the minimal allowed value
         IF (opt_embed%trust_rad .GE. opt_embed%min_trad) &
            opt_embed%trust_rad = 0.25_dp*opt_embed%trust_rad
      ENDIF

      IF (opt_embed%accept_step) opt_embed%last_accepted = opt_embed%i_iter

      CALL timestop(handle)

   END SUBROUTINE step_control

! **************************************************************************************************
!> \brief ...
!> \param opt_embed ...
!> \param diag_grad ...
!> \param eigenval ...
!> \param diag_step ...
! **************************************************************************************************
   SUBROUTINE level_shift(opt_embed, diag_grad, eigenval, diag_step)
      TYPE(opt_embed_pot_type)                           :: opt_embed
      TYPE(cp_fm_type), POINTER                          :: diag_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      TYPE(cp_fm_type), POINTER                          :: diag_step

      CHARACTER(LEN=*), PARAMETER :: routineN = 'level_shift', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: max_iter = 25
      REAL(KIND=dp), PARAMETER                           :: thresh = 0.00001_dp

      INTEGER                                            :: handle, i_iter, l_global, LLL, &
                                                            min_index, nrow_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: red_eigenval_map
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      LOGICAL                                            :: converged, do_shift
      REAL(KIND=dp) :: diag_grad_norm, grad_min, hess_min, shift, shift_max, shift_min, step_len, &
         step_minus_trad, step_minus_trad_first, step_minus_trad_max, step_minus_trad_min
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      ! Array properties
      CALL cp_fm_get_info(matrix=opt_embed%embed_pot_coef, &
                          nrow_local=nrow_local, &
                          row_indices=row_indices, &
                          para_env=para_env)

      min_index = MINLOC(ABS(eigenval), dim=1)
      hess_min = eigenval(min_index)
      CALL cp_fm_get_element(diag_grad, min_index, 1, grad_min)

      CALL cp_fm_trace(diag_grad, diag_grad, diag_grad_norm)

      IF (hess_min .LT. 0.0_dp) THEN
         !shift_min = -2.0_dp*(diag_grad_norm/opt_embed%trust_rad - min(hess_min, 0.0_dp))
         !shift_max = max(0.0_dp, -hess_min + 0.5_dp*grad_min/opt_embed%trust_rad)
         !shift_max = MIN(-hess_min+0.5_dp*grad_min/opt_embed%trust_rad, 0.0_dp)
         shift_max = hess_min + 0.1
         shift_min = diag_grad_norm/opt_embed%trust_rad
         shift_min = 10.0_dp
         !If (abs(shift_max) .LE. thresh) then
         !   shift_min = -20.0_dp*(diag_grad_norm/opt_embed%trust_rad - min(hess_min, 0.0_dp))
         !Else
         !   shift_min = 20.0_dp*shift_max
         !Endif

         ! The boundary values
         step_minus_trad_max = shifted_step(diag_grad, eigenval, shift_max, opt_embed%trust_rad)
         step_minus_trad_min = shifted_step(diag_grad, eigenval, shift_min, opt_embed%trust_rad)

         ! Find zero by bisection
         converged = .FALSE.
         do_shift = .FALSE.
         IF (ABS(step_minus_trad_max) .LE. thresh) THEN
            shift = shift_max
         ELSE
            IF (ABS(step_minus_trad_min) .LE. thresh) THEN
               shift = shift_min
            ELSE
               DO i_iter = 1, max_iter
                  shift = 0.5_dp*(shift_max + shift_min)
                  step_minus_trad = shifted_step(diag_grad, eigenval, shift, opt_embed%trust_rad)
                  IF (i_iter .EQ. 1) step_minus_trad_first = step_minus_trad
                  IF (step_minus_trad .GT. 0.0_dp) shift_max = shift
                  IF (step_minus_trad .LT. 0.0_dp) shift_min = shift
                  !IF (ABS(shift_max-shift_min) .LT. thresh) converged = .TRUE.
                  IF (ABS(step_minus_trad) .LT. thresh) converged = .TRUE.
                  IF (converged) EXIT
               ENDDO
               IF (ABS(step_minus_trad) .LT. ABS(step_minus_trad_first)) do_shift = .TRUE.
            ENDIF
         ENDIF
         ! Apply level-shifting
         IF (converged .OR. do_shift) THEN
            DO LLL = 1, nrow_local
               l_global = row_indices(LLL)
               IF (ABS(eigenval(l_global)) .GE. thresh) THEN
                  diag_step%local_data(LLL, 1) = &
                     -diag_grad%local_data(LLL, 1)/(eigenval(l_global) - shift)
               ELSE
                  diag_step%local_data(LLL, 1) = 0.0_dp
               ENDIF
            ENDDO
         ENDIF
         IF (.NOT. converged) THEN ! Scale if shift has not been found
            CALL cp_fm_trace(diag_step, diag_step, step_len)
            CALL cp_fm_scale(opt_embed%trust_rad/step_len, diag_step)
         ENDIF

         ! Special case
      ELSE ! Hess min .LT. 0.0_dp
         ! First, find all positive eigenvalues
         ALLOCATE (red_eigenval_map(opt_embed%dimen_var_aux))
         red_eigenval_map = 0
         DO LLL = 1, nrow_local
            l_global = row_indices(LLL)
            IF (eigenval(l_global) .GE. 0.0_dp) THEN
               red_eigenval_map(l_global) = 1
            ENDIF
         ENDDO
         CALL mp_sum(red_eigenval_map, para_env%group)

         ! Set shift as -hess_min and find step on the reduced space of negative-value
         ! eigenvectors
         shift = -hess_min
         DO LLL = 1, nrow_local
            l_global = row_indices(LLL)
            IF (red_eigenval_map(l_global) .EQ. 0) THEN
               IF (ABS(eigenval(l_global)) .GE. thresh) THEN
                  diag_step%local_data(LLL, 1) = &
                     -diag_grad%local_data(LLL, 1)/(eigenval(l_global) - shift)
               ELSE
                  diag_step%local_data(LLL, 1) = 0.0_dp
               ENDIF
            ELSE
               diag_step%local_data(LLL, 1) = 0.0_dp
            ENDIF
         ENDDO

         ! Find the step length of such a step
         CALL cp_fm_trace(diag_step, diag_step, step_len)

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE level_shift

! **************************************************************************************************
!> \brief ...
!> \param diag_grad ...
!> \param eigenval ...
!> \param shift ...
!> \param trust_rad ...
!> \return ...
! **************************************************************************************************
   FUNCTION shifted_step(diag_grad, eigenval, shift, trust_rad) RESULT(step_minus_trad)
      TYPE(cp_fm_type), POINTER                          :: diag_grad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      REAL(KIND=dp)                                      :: shift, trust_rad, step_minus_trad

      REAL(KIND=dp), PARAMETER                           :: thresh = 0.000001_dp

      INTEGER                                            :: l_global, LLL, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp)                                      :: step, step_1d
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL cp_fm_get_info(matrix=diag_grad, &
                          nrow_local=nrow_local, &
                          row_indices=row_indices, &
                          para_env=para_env)

      step = 0.0_dp
      DO LLL = 1, nrow_local
         l_global = row_indices(LLL)
         IF ((ABS(eigenval(l_global)) .GE. thresh) .AND. (ABS(diag_grad%local_data(LLL, 1)) .GE. thresh)) THEN
            step_1d = -diag_grad%local_data(LLL, 1)/(eigenval(l_global) + shift)
            step = step + step_1d**2
         ENDIF
      ENDDO

      CALL mp_sum(step, para_env%group)

      step_minus_trad = SQRT(step) - trust_rad

   END FUNCTION shifted_step

! **************************************************************************************************
!> \brief ...
!> \param step ...
!> \param prev_step ...
!> \param grad ...
!> \param prev_grad ...
!> \return ...
!> \retval length ...
! **************************************************************************************************
   FUNCTION Barzilai_Borwein(step, prev_step, grad, prev_grad) RESULT(length)
      TYPE(cp_fm_type), POINTER                          :: step, prev_step, grad, prev_grad
      REAL(KIND=dp)                                      :: length

      REAL(KIND=dp)                                      :: denominator, numerator
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: grad_diff, step_diff

      ! Get full matrix structures
      NULLIFY (fm_struct)

      CALL cp_fm_get_info(matrix=grad, &
                          matrix_struct=fm_struct)

      ! Allocate intermediates
      NULLIFY (grad_diff, step_diff)
      CALL cp_fm_create(grad_diff, fm_struct, name="grad_diff")
      CALL cp_fm_create(step_diff, fm_struct, name="step_diff")

      ! Calculate intermediates
      CALL cp_fm_to_fm(grad, grad_diff)
      CALL cp_fm_to_fm(step, step_diff)

      CALL cp_fm_scale_and_add(1.0_dp, grad_diff, -1.0_dp, prev_grad)
      CALL cp_fm_scale_and_add(1.0_dp, step_diff, -1.0_dp, prev_step)

      CALL cp_fm_trace(step_diff, grad_diff, numerator)
      CALL cp_fm_trace(grad_diff, grad_diff, denominator)

      ! Release intermediates
      CALL cp_fm_release(grad_diff)
      CALL cp_fm_release(step_diff)

      length = numerator/denominator

   END FUNCTION Barzilai_Borwein

! **************************************************************************************************
!> \brief ...
!> \param pw_env ...
!> \param embed_pot ...
!> \param spin_embed_pot ...
!> \param diff_rho_r ...
!> \param diff_rho_spin ...
!> \param rho_r_ref ...
!> \param open_shell_embed ...
!> \param step_len ...
! **************************************************************************************************
   SUBROUTINE Leeuwen_Baerends_potential_update(pw_env, embed_pot, spin_embed_pot, diff_rho_r, diff_rho_spin, &
                                                rho_r_ref, open_shell_embed, step_len)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), POINTER                           :: embed_pot, spin_embed_pot, diff_rho_r, &
                                                            diff_rho_spin
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r_ref
      LOGICAL                                            :: open_shell_embed
      REAL(KIND=dp)                                      :: step_len

      CHARACTER(LEN=*), PARAMETER :: routineN = 'Leeuwen_Baerends_potential_update', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, i_spin, j, k, nspins
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: my_rho, rho_cutoff
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: new_embed_pot, rho_n_1, temp_embed_pot
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      rho_cutoff = EPSILON(0.0_dp)

      ! Prepare plane-waves pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      NULLIFY (new_embed_pot)

      nspins = 1
      IF (open_shell_embed) nspins = 2
      NULLIFY (new_embed_pot)
      ALLOCATE (new_embed_pot(nspins))
      DO i_spin = 1, nspins
         NULLIFY (new_embed_pot(i_spin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, new_embed_pot(i_spin)%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(new_embed_pot(i_spin)%pw)
      ENDDO

      lb(1:3) = embed_pot%pw%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = embed_pot%pw%pw_grid%bounds_local(2, 1:3)

      IF (.NOT. open_shell_embed) THEN
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k, my_rho) &
!$OMP                SHARED(new_embed_pot, embed_pot, diff_rho_r, rho_r_ref, lb, ub, rho_cutoff, step_len)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  IF (rho_r_ref(1)%pw%cr3d(i, j, k) .GT. rho_cutoff) THEN
                     my_rho = rho_r_ref(1)%pw%cr3d(i, j, k)
                  ELSE
                     my_rho = rho_cutoff
                  ENDIF
                  new_embed_pot(1)%pw%cr3d(i, j, k) = step_len*embed_pot%pw%cr3d(i, j, k)* &
                                                      (diff_rho_r%pw%cr3d(i, j, k) + rho_r_ref(1)%pw%cr3d(i, j, k))/my_rho
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
         CALL pw_copy(new_embed_pot(1)%pw, embed_pot%pw)

      ELSE
         ! One has to work with spin components rather than with total and spin density
         NULLIFY (rho_n_1)
         ALLOCATE (rho_n_1(nspins))
         NULLIFY (temp_embed_pot)
         ALLOCATE (temp_embed_pot(nspins))
         DO i_spin = 1, nspins
            NULLIFY (rho_n_1(i_spin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_n_1(i_spin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(rho_n_1(i_spin)%pw)
            NULLIFY (temp_embed_pot(i_spin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, temp_embed_pot(i_spin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(temp_embed_pot(i_spin)%pw)
         ENDDO
         CALL pw_copy(diff_rho_r%pw, rho_n_1(1)%pw)
         CALL pw_copy(diff_rho_r%pw, rho_n_1(2)%pw)
         CALL pw_axpy(diff_rho_spin%pw, rho_n_1(1)%pw, 1.0_dp)
         CALL pw_axpy(diff_rho_spin%pw, rho_n_1(2)%pw, -1.0_dp)
         CALL pw_scale(rho_n_1(1)%pw, a=0.5_dp)
         CALL pw_scale(rho_n_1(2)%pw, a=0.5_dp)

         CALL pw_copy(embed_pot%pw, temp_embed_pot(1)%pw)
         CALL pw_copy(embed_pot%pw, temp_embed_pot(2)%pw)
         CALL pw_axpy(spin_embed_pot%pw, temp_embed_pot(1)%pw, 1.0_dp)
         CALL pw_axpy(spin_embed_pot%pw, temp_embed_pot(2)%pw, -1.0_dp)

         IF (SIZE(rho_r_ref) .EQ. 2) THEN
            CALL pw_axpy(rho_r_ref(1)%pw, rho_n_1(1)%pw, 1.0_dp)
            CALL pw_axpy(rho_r_ref(2)%pw, rho_n_1(2)%pw, 1.0_dp)

!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k, my_rho) &
!$OMP                SHARED(new_embed_pot, temp_embed_pot, rho_r_ref, rho_n_1, lb, ub, rho_cutoff, step_len)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (rho_r_ref(1)%pw%cr3d(i, j, k) .GT. rho_cutoff) THEN
                        my_rho = rho_r_ref(1)%pw%cr3d(i, j, k)
                     ELSE
                        my_rho = rho_cutoff
                     ENDIF
                     new_embed_pot(1)%pw%cr3d(i, j, k) = step_len*temp_embed_pot(1)%pw%cr3d(i, j, k)* &
                                                         (rho_n_1(1)%pw%cr3d(i, j, k))/my_rho
                     IF (rho_r_ref(2)%pw%cr3d(i, j, k) .GT. rho_cutoff) THEN
                        my_rho = rho_r_ref(2)%pw%cr3d(i, j, k)
                     ELSE
                        my_rho = rho_cutoff
                     ENDIF
                     new_embed_pot(2)%pw%cr3d(i, j, k) = step_len*temp_embed_pot(2)%pw%cr3d(i, j, k)* &
                                                         (rho_n_1(2)%pw%cr3d(i, j, k))/my_rho
                  END DO
               END DO
            END DO
!$OMP    END PARALLEL DO

         ELSE ! Reference system is closed-shell
            CALL pw_axpy(rho_r_ref(1)%pw, rho_n_1(1)%pw, 1.0_dp)
            ! The beta spin component is here equal to the difference: nothing to do

!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k, my_rho) &
!$OMP                SHARED(new_embed_pot, rho_r_ref, temp_embed_pot, rho_n_1, lb, ub, rho_cutoff, step_len)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (rho_r_ref(1)%pw%cr3d(i, j, k) .GT. rho_cutoff) THEN
                        my_rho = 0.5_dp*rho_r_ref(1)%pw%cr3d(i, j, k)
                     ELSE
                        my_rho = rho_cutoff
                     ENDIF
                     new_embed_pot(1)%pw%cr3d(i, j, k) = step_len*temp_embed_pot(1)%pw%cr3d(i, j, k)* &
                                                         (rho_n_1(1)%pw%cr3d(i, j, k))/my_rho
                     new_embed_pot(2)%pw%cr3d(i, j, k) = step_len*temp_embed_pot(2)%pw%cr3d(i, j, k)* &
                                                         (rho_n_1(2)%pw%cr3d(i, j, k))/my_rho
                  END DO
               END DO
            END DO
!$OMP    END PARALLEL DO
         ENDIF

         CALL pw_copy(new_embed_pot(1)%pw, embed_pot%pw)
         CALL pw_axpy(new_embed_pot(2)%pw, embed_pot%pw, 1.0_dp)
         CALL pw_scale(embed_pot%pw, a=0.5_dp)
         CALL pw_copy(new_embed_pot(1)%pw, spin_embed_pot%pw)
         CALL pw_axpy(new_embed_pot(2)%pw, spin_embed_pot%pw, -1.0_dp)
         CALL pw_scale(spin_embed_pot%pw, a=0.5_dp)

         DO i_spin = 1, nspins
            CALL pw_release(rho_n_1(i_spin)%pw)
            CALL pw_release(temp_embed_pot(i_spin)%pw)
         ENDDO
         DEALLOCATE (rho_n_1)
         DEALLOCATE (temp_embed_pot)
      ENDIF

      DO i_spin = 1, nspins
         CALL pw_release(new_embed_pot(i_spin)%pw)
      ENDDO

      DEALLOCATE (new_embed_pot)

      CALL timestop(handle)

   END SUBROUTINE Leeuwen_Baerends_potential_update

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho_r_ref ...
!> \param prev_embed_pot ...
!> \param prev_spin_embed_pot ...
!> \param embed_pot ...
!> \param spin_embed_pot ...
!> \param diff_rho_r ...
!> \param diff_rho_spin ...
!> \param v_w_ref ...
!> \param i_iter ...
!> \param step_len ...
!> \param open_shell_embed ...
!> \param vw_cutoff ...
!> \param vw_smooth_cutoff_range ...
! **************************************************************************************************
   SUBROUTINE FAB_update(qs_env, rho_r_ref, prev_embed_pot, prev_spin_embed_pot, embed_pot, spin_embed_pot, &
                         diff_rho_r, diff_rho_spin, v_w_ref, i_iter, step_len, open_shell_embed, &
                         vw_cutoff, vw_smooth_cutoff_range)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r_ref
      TYPE(pw_p_type), POINTER                           :: prev_embed_pot, prev_spin_embed_pot, &
                                                            embed_pot, spin_embed_pot, diff_rho_r, &
                                                            diff_rho_spin
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_w_ref
      INTEGER, INTENT(IN)                                :: i_iter
      REAL(KIND=dp)                                      :: step_len
      LOGICAL                                            :: open_shell_embed
      REAL(KIND=dp)                                      :: vw_cutoff, vw_smooth_cutoff_range

      CHARACTER(LEN=*), PARAMETER :: routineN = 'FAB_update', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_spin, nspins
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: curr_rho, new_embed_pot, temp_embed_pot, &
                                                            v_w
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      ! Update formula: v(n+1) = v(n-1) - v_w(ref) + v_w(n)

      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env)
      ! Get plane waves pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! We calculate von Weizsaecker potential for the reference density
      ! only at the first iteration
      IF (i_iter .LE. 1) THEN
         nspins = SIZE(rho_r_ref)
         NULLIFY (v_w_ref)
         ALLOCATE (v_w_ref(nspins))
         DO i_spin = 1, nspins
            NULLIFY (v_w_ref(i_spin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_w_ref(i_spin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
         ENDDO
         CALL Von_Weizsacker(rho_r_ref, v_w_ref, qs_env, vw_cutoff, vw_smooth_cutoff_range)
         ! For the first step previous are set to current
         CALL pw_copy(embed_pot%pw, prev_embed_pot%pw)
         CALL pw_axpy(diff_rho_r%pw, embed_pot%pw, 0.5_dp)
         IF (open_shell_embed) THEN
            CALL pw_copy(spin_embed_pot%pw, prev_spin_embed_pot%pw)
            CALL pw_axpy(diff_rho_r%pw, embed_pot%pw, 0.5_dp)
         ENDIF

      ELSE

         ! Reference can be closed shell, but total embedding - open shell:
         ! redefine nspins
         nspins = 1
         IF (open_shell_embed) nspins = 2
         NULLIFY (new_embed_pot)
         ALLOCATE (new_embed_pot(nspins))
         NULLIFY (v_w)
         ALLOCATE (v_w(nspins))
         NULLIFY (curr_rho)
         ALLOCATE (curr_rho(nspins))
         DO i_spin = 1, nspins
            NULLIFY (new_embed_pot(i_spin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, new_embed_pot(i_spin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(new_embed_pot(i_spin)%pw)

            NULLIFY (v_w(i_spin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, v_w(i_spin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(v_w(i_spin)%pw)

            NULLIFY (curr_rho(i_spin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, curr_rho(i_spin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(curr_rho(i_spin)%pw)
         ENDDO

         ! Now, deal with the current density

         IF (.NOT. open_shell_embed) THEN
            ! Reconstruct current density
            CALL pw_copy(diff_rho_r%pw, curr_rho(1)%pw)
            CALL pw_axpy(rho_r_ref(1)%pw, curr_rho(1)%pw, 1.0_dp)
            ! Compute von Weizsaecker potential
            CALL Von_Weizsacker(curr_rho, v_w, qs_env, vw_cutoff, vw_smooth_cutoff_range)
            ! Compute new embedding potential
            CALL pw_copy(prev_embed_pot%pw, new_embed_pot(1)%pw)
            CALL pw_axpy(v_w(1)%pw, new_embed_pot(1)%pw, step_len)
            CALL pw_axpy(v_w_ref(1)%pw, new_embed_pot(1)%pw, -step_len)
            ! Copy the potentials

            CALL pw_copy(embed_pot%pw, prev_embed_pot%pw)
            CALL pw_copy(new_embed_pot(1)%pw, embed_pot%pw)

         ELSE
            ! Reconstruct current density
            CALL pw_copy(diff_rho_r%pw, curr_rho(1)%pw)
            CALL pw_copy(diff_rho_r%pw, curr_rho(2)%pw)
            CALL pw_axpy(diff_rho_spin%pw, curr_rho(1)%pw, 1.0_dp)
            CALL pw_axpy(diff_rho_spin%pw, curr_rho(2)%pw, -1.0_dp)
            CALL pw_scale(curr_rho(1)%pw, a=0.5_dp)
            CALL pw_scale(curr_rho(2)%pw, a=0.5_dp)

            IF (SIZE(rho_r_ref) .EQ. 1) THEN ! If reference system is closed-shell
               CALL pw_axpy(rho_r_ref(1)%pw, curr_rho(1)%pw, 0.5_dp)
               CALL pw_axpy(rho_r_ref(1)%pw, curr_rho(2)%pw, 0.5_dp)
            ELSE ! If reference system is open-shell
               CALL pw_axpy(rho_r_ref(1)%pw, curr_rho(1)%pw, 1.0_dp)
               CALL pw_axpy(rho_r_ref(2)%pw, curr_rho(2)%pw, 1.0_dp)
            ENDIF

            ! Compute von Weizsaecker potential
            CALL Von_Weizsacker(curr_rho, v_w, qs_env, vw_cutoff, vw_smooth_cutoff_range)

            ! Reconstruct corrent spin components of the potential
            NULLIFY (temp_embed_pot)
            ALLOCATE (temp_embed_pot(nspins))
            DO i_spin = 1, nspins
               NULLIFY (temp_embed_pot(i_spin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, temp_embed_pot(i_spin)%pw, &
                                      use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_zero(temp_embed_pot(i_spin)%pw)
            ENDDO
            CALL pw_copy(embed_pot%pw, temp_embed_pot(1)%pw)
            CALL pw_copy(embed_pot%pw, temp_embed_pot(2)%pw)
            CALL pw_axpy(spin_embed_pot%pw, temp_embed_pot(1)%pw, 1.0_dp)
            CALL pw_axpy(spin_embed_pot%pw, temp_embed_pot(2)%pw, -1.0_dp)

            ! Compute new embedding potential
            IF (SIZE(v_w_ref) .EQ. 1) THEN ! Reference system is closed-shell
               CALL pw_copy(temp_embed_pot(1)%pw, new_embed_pot(1)%pw)
               CALL pw_axpy(v_w(1)%pw, new_embed_pot(1)%pw, 0.5_dp*step_len)
               CALL pw_axpy(v_w_ref(1)%pw, new_embed_pot(1)%pw, -0.5_dp*step_len)

               CALL pw_copy(temp_embed_pot(2)%pw, new_embed_pot(2)%pw)
               CALL pw_axpy(v_w(2)%pw, new_embed_pot(2)%pw, 0.5_dp)
               CALL pw_axpy(v_w_ref(1)%pw, new_embed_pot(2)%pw, -0.5_dp)

            ELSE ! Reference system is open-shell

               DO i_spin = 1, nspins
                  CALL pw_copy(temp_embed_pot(i_spin)%pw, new_embed_pot(i_spin)%pw)
                  CALL pw_axpy(v_w(1)%pw, new_embed_pot(i_spin)%pw, step_len)
                  CALL pw_axpy(v_w_ref(i_spin)%pw, new_embed_pot(i_spin)%pw, -step_len)
               ENDDO
            ENDIF

            ! Update embedding potentials
            CALL pw_copy(embed_pot%pw, prev_embed_pot%pw)
            CALL pw_copy(spin_embed_pot%pw, prev_spin_embed_pot%pw)

            CALL pw_copy(new_embed_pot(1)%pw, embed_pot%pw)
            CALL pw_axpy(new_embed_pot(2)%pw, embed_pot%pw, 1.0_dp)
            CALL pw_scale(embed_pot%pw, a=0.5_dp)
            CALL pw_copy(new_embed_pot(1)%pw, spin_embed_pot%pw)
            CALL pw_axpy(new_embed_pot(2)%pw, spin_embed_pot%pw, -1.0_dp)
            CALL pw_scale(spin_embed_pot%pw, a=0.5_dp)

            DO i_spin = 1, nspins
               CALL pw_release(temp_embed_pot(i_spin)%pw)
            ENDDO
            DEALLOCATE (temp_embed_pot)

         ENDIF

         DO i_spin = 1, nspins
            CALL pw_release(curr_rho(i_spin)%pw)
            CALL pw_release(new_embed_pot(i_spin)%pw)
            CALL pw_release(v_w(i_spin)%pw)
         ENDDO

         DEALLOCATE (new_embed_pot)
         DEALLOCATE (v_w)
         DEALLOCATE (curr_rho)

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE FAB_update

! **************************************************************************************************
!> \brief ...
!> \param rho_r ...
!> \param v_w ...
!> \param qs_env ...
!> \param vw_cutoff ...
!> \param vw_smooth_cutoff_range ...
! **************************************************************************************************
   SUBROUTINE Von_Weizsacker(rho_r, v_w, qs_env, vw_cutoff, vw_smooth_cutoff_range)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, v_w
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp)                                      :: vw_cutoff, vw_smooth_cutoff_range

      REAL(KIND=dp), PARAMETER                           :: one_4 = 0.25_dp, one_8 = 0.125_dp

      INTEGER                                            :: i, i_spin, j, k, nspins
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: density_smooth_cut_range, my_rho, &
                                                            rho_cutoff
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: rhoa, rhob
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, tau
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: input, xc_section
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho_set

      rho_cutoff = EPSILON(0.0_dp)

      nspins = SIZE(rho_r)

      NULLIFY (rho_set, xc_section)

      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env, &
                      input=input)

      ! Get plane waves pool
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! get some of the grids ready
      NULLIFY (rho_g)
      ALLOCATE (rho_g(nspins))
      DO i_spin = 1, nspins
         NULLIFY (rho_g(i_spin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_g(i_spin)%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_transfer(rho_r(i_spin)%pw, rho_g(i_spin)%pw)
      ENDDO

      xc_section => section_vals_get_subs_vals(input, "DFT%XC")

      CALL xc_rho_set_create(rho_set, &
                             rho_r(1)%pw%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "density_cutoff"), &
                             drho_cutoff=section_get_rval(xc_section, "gradient_cutoff"), &
                             tau_cutoff=section_get_rval(xc_section, "tau_cutoff"))

      CALL xc_rho_cflags_setall(needs, .FALSE.)

      IF (nspins .EQ. 2) THEN
         needs%rho_spin = .TRUE.
         needs%norm_drho_spin = .TRUE.
         needs%laplace_rho_spin = .TRUE.
      ELSE
         needs%rho = .TRUE.
         needs%norm_drho = .TRUE.
         needs%laplace_rho = .TRUE.
      ENDIF

      CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)

      CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", &
                                r_val=rho_cutoff)
      CALL section_vals_val_get(xc_section, "DENSITY_SMOOTH_CUTOFF_RANGE", &
                                r_val=density_smooth_cut_range)

      lb(1:3) = rho_r(1)%pw%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_r(1)%pw%pw_grid%bounds_local(2, 1:3)

      IF (nspins .EQ. 2) THEN
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k, my_rho) &
!$OMP                SHARED(v_w, rho_r, rho_set, lb, ub, rho_cutoff)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  IF (rho_r(1)%pw%cr3d(i, j, k) .GT. rho_cutoff) THEN
                     my_rho = rho_r(1)%pw%cr3d(i, j, k)
                  ELSE
                     my_rho = rho_cutoff
                  ENDIF
                  v_w(1)%pw%cr3d(i, j, k) = one_8*rho_set%norm_drhoa(i, j, k)**2/my_rho**2 - &
                                            one_4*rho_set%laplace_rhoa(i, j, k)/my_rho

                  IF (rho_r(2)%pw%cr3d(i, j, k) .GT. rho_cutoff) THEN
                     my_rho = rho_r(2)%pw%cr3d(i, j, k)
                  ELSE
                     my_rho = rho_cutoff
                  ENDIF
                  v_w(2)%pw%cr3d(i, j, k) = one_8*rho_set%norm_drhob(i, j, k)**2/my_rho**2 - &
                                            one_4*rho_set%laplace_rhob(i, j, k)/my_rho
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
      ELSE
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k, my_rho) &
!$OMP                SHARED(v_w, rho_r, rho_set, lb, ub, rho_cutoff)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  IF (rho_r(1)%pw%cr3d(i, j, k) .GT. rho_cutoff) THEN
                     my_rho = rho_r(1)%pw%cr3d(i, j, k)
                     v_w(1)%pw%cr3d(i, j, k) = one_8*rho_set%norm_drho(i, j, k)**2/my_rho**2 - &
                                               one_4*rho_set%laplace_rho(i, j, k)/my_rho
                  ELSE
                     v_w(1)%pw%cr3d(i, j, k) = 0.0_dp
                  ENDIF
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO

      ENDIF

      ! Smoothen the von Weizsaecker potential
      IF (nspins == 2) THEN
         density_smooth_cut_range = 0.5_dp*density_smooth_cut_range
         rho_cutoff = 0.5_dp*rho_cutoff
      ENDIF
      DO i_spin = 1, nspins
         CALL smooth_cutoff(pot=v_w(i_spin)%pw%cr3d, rho=rho_r(i_spin)%pw%cr3d, rhoa=rhoa, rhob=rhob, &
                            rho_cutoff=vw_cutoff, &
                            rho_smooth_cutoff_range=vw_smooth_cutoff_range)
      ENDDO

      CALL xc_rho_set_release(rho_set, pw_pool=auxbas_pw_pool)

      DO i_spin = 1, nspins
         CALL pw_release(rho_g(i_spin)%pw)
      ENDDO
      DEALLOCATE (rho_g)

   END SUBROUTINE Von_Weizsacker

! **************************************************************************************************
!> \brief ...
!> \param diff_rho_r ...
!> \return ...
! **************************************************************************************************
   FUNCTION max_dens_diff(diff_rho_r) RESULT(total_max_diff)
      TYPE(pw_p_type)                                    :: diff_rho_r
      REAL(KIND=dp)                                      :: total_max_diff

      INTEGER                                            :: size_x, size_y, size_z
      REAL(KIND=dp)                                      :: max_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: grid_3d

      !, i_x, i_y, i_z

      ! Get the sizes
      size_x = SIZE(diff_rho_r%pw%cr3d, 1)
      size_y = SIZE(diff_rho_r%pw%cr3d, 2)
      size_z = SIZE(diff_rho_r%pw%cr3d, 3)

      ! Allocate the density
      ALLOCATE (grid_3d(size_x, size_y, size_z))

      ! Copy density
      grid_3d(:, :, :) = diff_rho_r%pw%cr3d(:, :, :)

      ! Find the maximum absolute value
      max_diff = MAXVAL(ABS(grid_3d))
      total_max_diff = max_diff
      CALL mp_max(total_max_diff, diff_rho_r%pw%pw_grid%para%group)

      ! Deallocate the density
      DEALLOCATE (grid_3d)

   END FUNCTION max_dens_diff

! **************************************************************************************************
!> \brief Prints a cube for the (rho_A + rho_B - rho_ref) to be minimized in embedding
!> \param diff_rho_r ...
!> \param i_iter ...
!> \param qs_env ...
!> \param final_one ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE print_rho_diff(diff_rho_r, i_iter, qs_env, final_one)
      TYPE(pw_p_type), POINTER                           :: diff_rho_r
      INTEGER                                            :: i_iter
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: final_one

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube, title
      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      NULLIFY (subsys, input)

      CALL get_qs_env(qs_env=qs_env, &
                      subsys=subsys, &
                      input=input)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL qs_subsys_get(subsys, particles=particles)

      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_EMBED%EMBED_DENS_DIFF"), cp_p_file)) THEN
         my_pos_cube = "REWIND"
         IF (.NOT. final_one) THEN
            WRITE (filename, '(a5,I3.3,a1,I1.1)') "DIFF_", i_iter
         ELSE
            WRITE (filename, '(a5,I3.3,a1,I1.1)') "DIFF"
         ENDIF
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_EMBED%EMBED_DENS_DIFF", &
                                        extension=".cube", middle_name=TRIM(filename), file_position=my_pos_cube, &
                                        log_filename=.FALSE.)

         WRITE (title, *) "EMBEDDING DENSITY DIFFERENCE ", " optimization step ", i_iter
         CALL cp_pw_to_cube(diff_rho_r%pw, unit_nr, title, particles=particles, &
                            stride=section_get_ivals(dft_section, "QS%OPT_EMBED%EMBED_DENS_DIFF%STRIDE"))
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%QS%OPT_EMBED%EMBED_DENS_DIFF")
      ENDIF

   END SUBROUTINE print_rho_diff

! **************************************************************************************************
!> \brief Prints a cube for the (spin_rho_A + spin_rho_B - spin_rho_ref) to be minimized in embedding
!> \param spin_diff_rho_r ...
!> \param i_iter ...
!> \param qs_env ...
!> \param final_one ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE print_rho_spin_diff(spin_diff_rho_r, i_iter, qs_env, final_one)
      TYPE(pw_p_type), POINTER                           :: spin_diff_rho_r
      INTEGER                                            :: i_iter
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: final_one

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube, title
      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      NULLIFY (subsys, input)

      CALL get_qs_env(qs_env=qs_env, &
                      subsys=subsys, &
                      input=input)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL qs_subsys_get(subsys, particles=particles)

      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_EMBED%EMBED_DENS_DIFF"), cp_p_file)) THEN
         my_pos_cube = "REWIND"
         IF (.NOT. final_one) THEN
            WRITE (filename, '(a5,I3.3,a1,I1.1)') "SPIN_DIFF_", i_iter
         ELSE
            WRITE (filename, '(a9,I3.3,a1,I1.1)') "SPIN_DIFF"
         ENDIF
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_EMBED%EMBED_DENS_DIFF", &
                                        extension=".cube", middle_name=TRIM(filename), file_position=my_pos_cube, &
                                        log_filename=.FALSE.)

         WRITE (title, *) "EMBEDDING SPIN DENSITY DIFFERENCE ", " optimization step ", i_iter
         CALL cp_pw_to_cube(spin_diff_rho_r%pw, unit_nr, title, particles=particles, &
                            stride=section_get_ivals(dft_section, "QS%OPT_EMBED%EMBED_DENS_DIFF%STRIDE"))
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%QS%OPT_EMBED%EMBED_DENS_DIFF")
      ENDIF

   END SUBROUTINE print_rho_spin_diff
! **************************************************************************************************
!> \brief Print embedding potential as a cube and as a binary (for restarting)
!> \param qs_env ...
!> \param dimen_aux ...
!> \param embed_pot_coef ...
!> \param embed_pot ...
!> \param i_iter ...
!> \param embed_pot_spin ...
!> \param open_shell_embed ...
!> \param grid_opt ...
!> \param final_one ...
! **************************************************************************************************
   SUBROUTINE print_embed_restart(qs_env, dimen_aux, embed_pot_coef, embed_pot, i_iter, &
                                  embed_pot_spin, open_shell_embed, grid_opt, final_one)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: dimen_aux
      TYPE(cp_fm_type), POINTER                          :: embed_pot_coef
      TYPE(pw_p_type), POINTER                           :: embed_pot
      INTEGER                                            :: i_iter
      TYPE(pw_p_type), POINTER                           :: embed_pot_spin
      LOGICAL                                            :: open_shell_embed, grid_opt, final_one

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube, title
      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      NULLIFY (input)
      CALL get_qs_env(qs_env=qs_env, subsys=subsys, &
                      input=input)

      ! First we print an unformatted file
      IF (.NOT. grid_opt) THEN ! Only for finite basis optimization
         logger => cp_get_default_logger()
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                              "DFT%QS%OPT_EMBED%EMBED_POT_VECTOR"), cp_p_file)) THEN
            IF (.NOT. final_one) THEN
               WRITE (filename, '(a10,I3.3)') "embed_pot_", i_iter
            ELSE
               WRITE (filename, '(a10,I3.3)') "embed_pot"
            ENDIF
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_EMBED%EMBED_POT_VECTOR", extension=".wfn", &
                                           file_form="UNFORMATTED", middle_name=TRIM(filename), file_position="REWIND")
            IF (unit_nr > 0) THEN
               WRITE (unit_nr) dimen_aux
            ENDIF
            CALL cp_fm_write_unformatted(embed_pot_coef, unit_nr)
            IF (unit_nr > 0) THEN
               CALL close_file(unit_nr)
            ENDIF
         ENDIF
      ENDIF

      ! Second, cube files
      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL qs_subsys_get(subsys, particles=particles)

      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_EMBED%EMBED_POT_CUBE"), cp_p_file)) THEN
         my_pos_cube = "REWIND"
         IF (.NOT. final_one) THEN
            WRITE (filename, '(a10,I3.3)') "embed_pot_", i_iter
         ELSE
            WRITE (filename, '(a10,I3.3)') "embed_pot"
         ENDIF
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_EMBED%EMBED_POT_CUBE", &
                                        extension=".cube", middle_name=TRIM(filename), file_position=my_pos_cube, &
                                        log_filename=.FALSE.)

         WRITE (title, *) "EMBEDDING POTENTIAL at optimization step ", i_iter
         CALL cp_pw_to_cube(embed_pot%pw, unit_nr, title, particles=particles)
!, &
!                            stride=section_get_ivals(dft_section, "QS%OPT_EMBED%EMBED_POT_CUBE%STRIDE"))
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%QS%OPT_EMBED%EMBED_POT_CUBE")
         IF (open_shell_embed) THEN ! Print spin part of the embedding potential
            my_pos_cube = "REWIND"
            IF (.NOT. final_one) THEN
               WRITE (filename, '(a15,I3.3)') "spin_embed_pot_", i_iter
            ELSE
               WRITE (filename, '(a15,I3.3)') "spin_embed_pot"
            ENDIF
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_EMBED%EMBED_POT_CUBE", &
                                           extension=".cube", middle_name=TRIM(filename), file_position=my_pos_cube, &
                                           log_filename=.FALSE.)

            WRITE (title, *) "SPIN EMBEDDING POTENTIAL at optimization step ", i_iter
            CALL cp_pw_to_cube(embed_pot_spin%pw, unit_nr, title, particles=particles)
!,  &
!                               stride=section_get_ivals(dft_section, "QS%OPT_EMBED%EMBED_POT_CUBE%STRIDE"))
            CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                              "DFT%QS%OPT_EMBED%EMBED_POT_CUBE")
         ENDIF
      ENDIF

   END SUBROUTINE print_embed_restart

! **************************************************************************************************
!> \brief Prints a volumetric file: X Y Z value for interfacing with external programs.
!> \param qs_env ...
!> \param embed_pot ...
!> \param embed_pot_spin ...
!> \param i_iter ...
!> \param open_shell_embed ...
!> \param final_one ...
!> \param qs_env_cluster ...
! **************************************************************************************************
   SUBROUTINE print_pot_simple_grid(qs_env, embed_pot, embed_pot_spin, i_iter, open_shell_embed, &
                                    final_one, qs_env_cluster)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: embed_pot, embed_pot_spin
      INTEGER                                            :: i_iter
      LOGICAL                                            :: open_shell_embed, final_one
      TYPE(qs_environment_type), POINTER                 :: qs_env_cluster

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: my_units, unit_nr
      LOGICAL                                            :: angstrom, bohr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), POINTER                           :: pot_alpha, pot_beta
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      NULLIFY (input)
      CALL get_qs_env(qs_env=qs_env, input=input, pw_env=pw_env)

      ! Second, cube files
      dft_section => section_vals_get_subs_vals(input, "DFT")

      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_EMBED%WRITE_SIMPLE_GRID"), cp_p_file)) THEN

         ! Figure out the units
         angstrom = .FALSE.
         bohr = .TRUE.
         CALL section_vals_val_get(dft_section, "QS%OPT_EMBED%WRITE_SIMPLE_GRID%UNITS", i_val=my_units)
         SELECT CASE (my_units)
         CASE (embed_grid_bohr)
            bohr = .TRUE.
            angstrom = .FALSE.
         CASE (embed_grid_angstrom)
            bohr = .FALSE.
            angstrom = .TRUE.
         CASE DEFAULT
            bohr = .TRUE.
            angstrom = .FALSE.
         END SELECT

         ! Get alpha and beta potentials
         ! Prepare plane-waves pool
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

         ! Create embedding potential and set to zero
         NULLIFY (pot_alpha)
         ALLOCATE (pot_alpha)
         NULLIFY (pot_alpha%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, pot_alpha%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(pot_alpha%pw)

         CALL pw_copy(embed_pot%pw, pot_alpha%pw)

         IF (open_shell_embed) THEN
            NULLIFY (pot_beta)
            ALLOCATE (pot_beta)
            NULLIFY (pot_beta%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, pot_beta%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_copy(embed_pot%pw, pot_beta%pw)
            ! Add spin potential to the alpha, and subtract from the beta part
            CALL pw_axpy(embed_pot_spin%pw, pot_alpha%pw, 1.0_dp)
            CALL pw_axpy(embed_pot_spin%pw, pot_beta%pw, -1.0_dp)
         ENDIF

         IF (.NOT. final_one) THEN
            WRITE (filename, '(a10,I3.3)') "embed_pot_", i_iter
         ELSE
            WRITE (filename, '(a10,I3.3)') "embed_pot"
         ENDIF
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%QS%OPT_EMBED%WRITE_SIMPLE_GRID", extension=".dat", &
                                        middle_name=TRIM(filename), file_form="FORMATTED", file_position="REWIND")

         IF (open_shell_embed) THEN ! Print spin part of the embedding potential
            CALL cp_pw_to_simple_volumetric(pw=pot_alpha%pw, unit_nr=unit_nr, &
                                            stride=section_get_ivals(dft_section, "QS%OPT_EMBED%WRITE_SIMPLE_GRID%STRIDE"), &
                                            pw2=pot_beta%pw)
         ELSE
            CALL cp_pw_to_simple_volumetric(pot_alpha%pw, unit_nr, &
                                            stride=section_get_ivals(dft_section, "QS%OPT_EMBED%WRITE_SIMPLE_GRID%STRIDE"))
         ENDIF

         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%QS%OPT_EMBED%WRITE_SIMPLE_GRID")
         ! Release structures
         CALL pw_release(pot_alpha%pw)
         DEALLOCATE (pot_alpha)
         IF (open_shell_embed) THEN
            CALL pw_release(pot_beta%pw)
            DEALLOCATE (pot_beta)
         ENDIF

      ENDIF

      ! Fold the coordinates and write into separate file: needed to have the grid correspond to coordinates
      ! Needed for external software.
      CALL print_folded_coordinates(qs_env_cluster, input)

   END SUBROUTINE print_pot_simple_grid

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input ...
! **************************************************************************************************
   SUBROUTINE print_folded_coordinates(qs_env, input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:)        :: particles_el
      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: iat, n, unit_nr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: particles_r
      REAL(KIND=dp), DIMENSION(3)                        :: center, r_pbc, s
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(qs_subsys_type), POINTER                      :: subsys

      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%QS%OPT_EMBED%WRITE_SIMPLE_GRID/FOLD_COORD"), cp_p_file)) THEN
         CALL get_qs_env(qs_env=qs_env, cell=cell, subsys=subsys)
         CALL qs_subsys_get(subsys=subsys, particles=particles)

         ! Prepare the file
         WRITE (filename, '(a14)') "folded_cluster"
         unit_nr = cp_print_key_unit_nr(logger, input, &
                                        "DFT%QS%OPT_EMBED%WRITE_SIMPLE_GRID/FOLD_COORD", extension=".dat", &
                                        middle_name=TRIM(filename), file_form="FORMATTED", file_position="REWIND")
         IF (unit_nr > 0) THEN

            n = particles%n_els
            ALLOCATE (particles_el(n))
            ALLOCATE (particles_r(3, n))
            DO iat = 1, n
               CALL get_atomic_kind(particles%els(iat)%atomic_kind, element_symbol=particles_el(iat))
               particles_r(:, iat) = particles%els(iat)%r(:)
            END DO

            ! Fold the coordinates
            center(:) = cell%hmat(:, 1)/2.0_dp + cell%hmat(:, 2)/2.0_dp + cell%hmat(:, 3)/2.0_dp

            ! Print folded coordinates to file
            DO iat = 1, SIZE(particles_el)
               r_pbc(:) = particles_r(:, iat) - center
               s = MATMUL(cell%h_inv, r_pbc)
               s = s - ANINT(s)
               r_pbc = MATMUL(cell%hmat, s)
               r_pbc = r_pbc + center
               WRITE (unit_nr, '(a4,4f12.6)') particles_el(iat), r_pbc(:)
            END DO

            CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                              "DFT%QS%OPT_EMBED%WRITE_SIMPLE_GRID/FOLD_COORD")

            DEALLOCATE (particles_el)
            DEALLOCATE (particles_r)
         ENDIF

      ENDIF ! Should output

   END SUBROUTINE print_folded_coordinates

! **************************************************************************************************
!> \brief ...
!> \param output_unit ...
!> \param step_num ...
!> \param opt_embed ...
! **************************************************************************************************
   SUBROUTINE print_emb_opt_info(output_unit, step_num, opt_embed)
      INTEGER                                            :: output_unit, step_num
      TYPE(opt_embed_pot_type)                           :: opt_embed

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,8('-'),A,I5,1X,12('-'))") &
            "  Optimize embedding potential info at step = ", step_num
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Functional value         = ", opt_embed%w_func(step_num)
         IF (step_num .GT. 1) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Real energy change         = ", opt_embed%w_func(step_num) - &
               opt_embed%w_func(step_num - 1)

            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Step size                  = ", opt_embed%step_len

         END IF

         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Trust radius               = ", opt_embed%trust_rad

         WRITE (UNIT=output_unit, FMT="(T2,51('-'))")
      END IF

   END SUBROUTINE print_emb_opt_info

! **************************************************************************************************
!> \brief ...
!> \param opt_embed ...
!> \param force_env ...
!> \param subsys_num ...
! **************************************************************************************************
   SUBROUTINE get_prev_density(opt_embed, force_env, subsys_num)
      TYPE(opt_embed_pot_type)                           :: opt_embed
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER                                            :: subsys_num

      INTEGER                                            :: i_dens_start, i_spin, nspins
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (rho_r, rho)
      CALL get_qs_env(force_env%qs_env, rho=rho)
      CALL qs_rho_get(rho_struct=rho, rho_r=rho_r)

      nspins = opt_embed%all_nspins(subsys_num)

      i_dens_start = SUM(opt_embed%all_nspins(1:subsys_num)) - nspins + 1

      DO i_spin = 1, nspins
         opt_embed%prev_subsys_dens(i_dens_start + i_spin - 1)%pw%cr3d(:, :, :) = &
            rho_r(i_spin)%pw%cr3d(:, :, :)
      ENDDO

   END SUBROUTINE get_prev_density

! **************************************************************************************************
!> \brief ...
!> \param opt_embed ...
!> \param force_env ...
!> \param subsys_num ...
! **************************************************************************************************
   SUBROUTINE get_max_subsys_diff(opt_embed, force_env, subsys_num)
      TYPE(opt_embed_pot_type)                           :: opt_embed
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER                                            :: subsys_num

      INTEGER                                            :: i_dens_start, i_spin, nspins
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (rho_r, rho)
      CALL get_qs_env(force_env%qs_env, rho=rho)
      CALL qs_rho_get(rho_struct=rho, rho_r=rho_r)

      nspins = opt_embed%all_nspins(subsys_num)

      i_dens_start = SUM(opt_embed%all_nspins(1:subsys_num)) - nspins + 1

      DO i_spin = 1, nspins
         opt_embed%prev_subsys_dens(i_dens_start + i_spin - 1)%pw%cr3d(:, :, :) = &
            rho_r(i_spin)%pw%cr3d(:, :, :) - opt_embed%prev_subsys_dens(i_dens_start + i_spin - 1)%pw%cr3d(:, :, :)
         opt_embed%max_subsys_dens_diff(i_dens_start + i_spin - 1) = &
            max_dens_diff(opt_embed%prev_subsys_dens(i_dens_start + i_spin - 1))
      ENDDO

   END SUBROUTINE get_max_subsys_diff

! **************************************************************************************************
!> \brief ...
!> \param opt_embed ...
!> \param diff_rho_r ...
!> \param diff_rho_spin ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE conv_check_embed(opt_embed, diff_rho_r, diff_rho_spin, output_unit)
      TYPE(opt_embed_pot_type)                           :: opt_embed
      TYPE(pw_p_type), POINTER                           :: diff_rho_r, diff_rho_spin
      INTEGER                                            :: output_unit

      INTEGER                                            :: i_dens, i_dens_start, i_spin
      LOGICAL                                            :: conv_int_diff, conv_max_diff
      REAL(KIND=dp)                                      :: int_diff, int_diff_spin, &
                                                            int_diff_square, int_diff_square_spin, &
                                                            max_diff, max_diff_spin

      ! Calculate the convergence target values
      opt_embed%max_diff(1) = max_dens_diff(diff_rho_r)
      opt_embed%int_diff(1) = pw_integrate_function(fun=diff_rho_r%pw, oprt='ABS')
      opt_embed%int_diff_square(1) = pw_integral_ab(diff_rho_r%pw, diff_rho_r%pw)
      IF (opt_embed%open_shell_embed) THEN
         opt_embed%max_diff(2) = max_dens_diff(diff_rho_spin)
         opt_embed%int_diff(2) = pw_integrate_function(fun=diff_rho_spin%pw, oprt='ABS')
         opt_embed%int_diff_square(2) = pw_integral_ab(diff_rho_spin%pw, diff_rho_spin%pw)
      ENDIF

      ! Find out the convergence
      max_diff = opt_embed%max_diff(1)

      ! Maximum value criterium
      ! Open shell
      IF (opt_embed%open_shell_embed) THEN
         max_diff_spin = opt_embed%max_diff(2)
         IF ((max_diff .LE. opt_embed%conv_max) .AND. (max_diff_spin .LE. opt_embed%conv_max_spin)) THEN
            conv_max_diff = .TRUE.
         ELSE
            conv_max_diff = .FALSE.
         ENDIF
      ELSE
         ! Closed shell
         IF (max_diff .LE. opt_embed%conv_max) THEN
            conv_max_diff = .TRUE.
         ELSE
            conv_max_diff = .FALSE.
         ENDIF
      ENDIF

      ! Integrated value criterium
      int_diff = opt_embed%int_diff(1)
      ! Open shell
      IF (opt_embed%open_shell_embed) THEN
         int_diff_spin = opt_embed%int_diff(2)
         IF ((int_diff .LE. opt_embed%conv_int) .AND. (int_diff_spin .LE. opt_embed%conv_int_spin)) THEN
            conv_int_diff = .TRUE.
         ELSE
            conv_int_diff = .FALSE.
         ENDIF
      ELSE
         ! Closed shell
         IF (int_diff .LE. opt_embed%conv_int) THEN
            conv_int_diff = .TRUE.
         ELSE
            conv_int_diff = .FALSE.
         ENDIF
      ENDIF

      ! Integrated squared value criterium
      int_diff_square = opt_embed%int_diff_square(1)
      ! Open shell
      IF (opt_embed%open_shell_embed) int_diff_square_spin = opt_embed%int_diff_square(2)

      IF ((conv_max_diff) .AND. (conv_int_diff)) THEN
         opt_embed%converged = .TRUE.
      ELSE
         opt_embed%converged = .FALSE.
      ENDIF

      ! Print the information
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            " Convergence check :"

         ! Maximum value of density
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Maximum density difference                = ", max_diff
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Convergence limit for max. density diff.  = ", opt_embed%conv_max

         IF (opt_embed%open_shell_embed) THEN

            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Maximum spin density difference           = ", max_diff_spin
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Convergence limit for max. spin dens.diff.= ", opt_embed%conv_max_spin

         ENDIF

         IF (conv_max_diff) THEN
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in max. density diff.    =     ", &
               "             YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in max. density diff.    =     ", &
               "              NO"
         END IF

         ! Integrated abs. value of density
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Integrated density difference             = ", int_diff
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Conv. limit for integrated density diff.  = ", opt_embed%conv_int
         IF (opt_embed%open_shell_embed) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Integrated spin density difference        = ", int_diff_spin
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Conv. limit for integrated spin dens.diff.= ", opt_embed%conv_int_spin
         ENDIF

         IF (conv_int_diff) THEN
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in integrated density diff.    =     ", &
               "             YES"
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,2A)") &
               " Convergence in integrated density diff.    =     ", &
               "              NO"
         END IF

         ! Integrated squared value of density
         WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
            " Integrated squared density difference     = ", int_diff_square
         IF (opt_embed%open_shell_embed) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,F20.10)") &
               " Integrated squared spin density difference= ", int_diff_square_spin
         ENDIF

         ! Maximum subsystem density change
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            " Maximum density change in:"
         DO i_dens = 1, (SIZE(opt_embed%all_nspins) - 1)
            i_dens_start = SUM(opt_embed%all_nspins(1:i_dens)) - opt_embed%all_nspins(i_dens) + 1
            DO i_spin = 1, opt_embed%all_nspins(i_dens)
               WRITE (UNIT=output_unit, FMT="(T4,A10,I3,A6,I3,A1,F20.10)") &
                  " subsystem ", i_dens, ', spin', i_spin, ":", &
                  opt_embed%max_subsys_dens_diff(i_dens_start + i_spin - 1)
            ENDDO
         ENDDO

      ENDIF

      IF ((opt_embed%converged) .AND. (output_unit > 0)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") REPEAT("*", 79)
         WRITE (UNIT=output_unit, FMT="(T2,A,T25,A,T78,A)") &
            "***", "EMBEDDING POTENTIAL OPTIMIZATION COMPLETED", "***"
         WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", 79)
      END IF

   END SUBROUTINE conv_check_embed

END MODULE optimize_embedding_potential
