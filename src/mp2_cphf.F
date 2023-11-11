!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate CPHF like update and solve Z-vector equation
!>        for MP2 gradients (only GPW)
!> \par History
!>      11.2013 created [Mauro Del Ben]
! **************************************************************************************************
MODULE mp2_cphf
   USE admm_methods,                    ONLY: admm_projection_derivative
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_scale,&
                                              dbcsr_set
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_exx,                         ONLY: add_exx_to_rhs
   USE hfx_ri,                          ONLY: hfx_ri_update_forces
   USE hfx_types,                       ONLY: alloc_containers,&
                                              hfx_container_type,&
                                              hfx_init_container,&
                                              hfx_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              ot_precond_full_all,&
                                              z_solver_cg,&
                                              z_solver_pople,&
                                              z_solver_richardson,&
                                              z_solver_sd
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE kinds,                           ONLY: dp
   USE linear_systems,                  ONLY: solve_system
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: fourpi
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type,&
                                              ri_rpa_method_gpw
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_integral_ab,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_2nd_kernel_ao,                ONLY: apply_2nd_order_kernel
   USE qs_density_matrices,             ONLY: calculate_whz_matrix
   USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: deallocate_qs_force,&
                                              qs_force_type,&
                                              sum_qs_force,&
                                              zero_qs_force
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_reference,                 ONLY: ks_ref_potential
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_linres_types,                 ONLY: linres_control_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_p_env_methods,                ONLY: p_env_check_i_alloc,&
                                              p_env_create,&
                                              p_env_psi0_changed,&
                                              p_env_update_rho
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE task_list_types,                 ONLY: task_list_type
   USE virial_types,                    ONLY: virial_type,&
                                              zero_virial

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_cphf'
   LOGICAL, PARAMETER, PRIVATE :: debug_forces = .TRUE.

   PUBLIC :: solve_z_vector_eq, update_mp2_forces

CONTAINS

! **************************************************************************************************
!> \brief Solve Z-vector equations necessary for the calculation of the MP2
!>        gradients, in order to be consistent here the parameters for the
!>        calculation of the CPHF like updats have to be exactly equal to the
!>        SCF case
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param dft_control ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param Eigenval ...
!> \param unit_nr ...
!> \author Mauro Del Ben, Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE solve_z_vector_eq(qs_env, mp2_env, para_env, dft_control, &
                                mo_coeff, nmo, homo, Eigenval, unit_nr)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(dft_control_type), INTENT(IN), POINTER        :: dft_control
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, INTENT(IN)                                :: nmo
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_z_vector_eq'

      INTEGER :: bin, dimen, handle, handle2, i, i_global, i_thread, iiB, irep, ispin, j_global, &
         jjB, my_bin_size, n_rep_hf, n_threads, ncol_local, nrow_local, nspins, transf_type_in
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: virtual
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: alpha_beta, do_dynamic_load_balancing, &
                                                            do_exx, do_hfx, restore_p_screen
      REAL(KIND=dp)                                      :: focc
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: fm_back, fm_G_mu_nu
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: L_jb, mo_coeff_o, mo_coeff_v, P_ia, &
                                                            P_mo, W_Mo
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_p_mp2, &
                                                            matrix_p_mp2_admm, matrix_s, P_mu_nu, &
                                                            rho1_ao, rho_ao, rho_ao_aux_fit
      TYPE(hfx_container_type), DIMENSION(:), POINTER    :: integral_containers
      TYPE(hfx_container_type), POINTER                  :: maxval_container
      TYPE(hfx_type), POINTER                            :: actual_x_data
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_section, hfx_sections, input

      CALL timeset(routineN, handle)

      ! start collecting stuff
      dimen = nmo
      NULLIFY (input, matrix_s, blacs_env, rho, &
               matrix_p_mp2, matrix_p_mp2_admm, matrix_ks)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      input=input, &
                      matrix_s=matrix_s, &
                      matrix_ks=matrix_ks, &
                      matrix_p_mp2=matrix_p_mp2, &
                      matrix_p_mp2_admm=matrix_p_mp2_admm, &
                      blacs_env=blacs_env, &
                      rho=rho)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ! Get number of relevant spin states
      nspins = dft_control%nspins
      alpha_beta = (nspins == 2)

      CALL MOVE_ALLOC(mp2_env%ri_grad%P_mo, P_mo)
      CALL MOVE_ALLOC(mp2_env%ri_grad%W_mo, W_mo)
      CALL MOVE_ALLOC(mp2_env%ri_grad%L_jb, L_jb)

      ALLOCATE (virtual(nspins))
      virtual(:) = dimen - homo(:)

      NULLIFY (P_mu_nu)
      CALL dbcsr_allocate_matrix_set(P_mu_nu, nspins)
      DO ispin = 1, nspins
         ALLOCATE (P_mu_nu(ispin)%matrix)
         CALL dbcsr_copy(P_mu_nu(ispin)%matrix, rho_ao(1)%matrix, name="P_mu_nu")
         CALL dbcsr_set(P_mu_nu(ispin)%matrix, 0.0_dp)
      END DO

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen, ncol_global=dimen)
      CALL cp_fm_create(fm_G_mu_nu, fm_struct_tmp, name="G_mu_nu")
      CALL cp_fm_create(fm_back, fm_struct_tmp, name="fm_back")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_set_all(fm_G_mu_nu, 0.0_dp)
      CALL cp_fm_set_all(fm_back, 0.0_dp)

      ALLOCATE (mo_coeff_o(nspins), mo_coeff_v(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=dimen, ncol_global=homo(ispin))
         CALL cp_fm_create(mo_coeff_o(ispin), fm_struct_tmp, name="mo_coeff_o")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(mo_coeff_o(ispin), 0.0_dp)
         CALL cp_fm_to_fm_submat(msource=mo_coeff(ispin), mtarget=mo_coeff_o(ispin), &
                                 nrow=dimen, ncol=homo(ispin), &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=1)

         NULLIFY (fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=dimen, ncol_global=virtual(ispin))
         CALL cp_fm_create(mo_coeff_v(ispin), fm_struct_tmp, name="mo_coeff_v")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(mo_coeff_v(ispin), 0.0_dp)
         CALL cp_fm_to_fm_submat(msource=mo_coeff(ispin), mtarget=mo_coeff_v(ispin), &
                                 nrow=dimen, ncol=virtual(ispin), &
                                 s_firstrow=1, s_firstcol=homo(ispin) + 1, &
                                 t_firstrow=1, t_firstcol=1)
      END DO

      ! hfx section
      NULLIFY (hfx_sections)
      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx, n_repetition=n_rep_hf)
      IF (do_hfx) THEN
         ! here we check if we have to reallocate the HFX container
         IF (mp2_env%ri_grad%free_hfx_buffer .AND. (.NOT. qs_env%x_data(1, 1)%do_hfx_ri)) THEN
            CALL timeset(routineN//"_alloc_hfx", handle2)
            n_threads = 1
!$          n_threads = omp_get_max_threads()

            DO irep = 1, n_rep_hf
               DO i_thread = 0, n_threads - 1
                  actual_x_data => qs_env%x_data(irep, i_thread + 1)

                  do_dynamic_load_balancing = .TRUE.
                  IF (n_threads == 1 .OR. actual_x_data%memory_parameter%do_disk_storage) do_dynamic_load_balancing = .FALSE.

                  IF (do_dynamic_load_balancing) THEN
                     my_bin_size = SIZE(actual_x_data%distribution_energy)
                  ELSE
                     my_bin_size = 1
                  END IF

                  IF (.NOT. actual_x_data%memory_parameter%do_all_on_the_fly) THEN
                     CALL alloc_containers(actual_x_data%store_ints, my_bin_size)

                     DO bin = 1, my_bin_size
                        maxval_container => actual_x_data%store_ints%maxval_container(bin)
                        integral_containers => actual_x_data%store_ints%integral_containers(:, bin)
                        CALL hfx_init_container(maxval_container, actual_x_data%memory_parameter%actual_memory_usage, .FALSE.)
                        DO i = 1, 64
                           CALL hfx_init_container(integral_containers(i), &
                                                   actual_x_data%memory_parameter%actual_memory_usage, .FALSE.)
                        END DO
                     END DO
                  END IF
               END DO
            END DO
            CALL timestop(handle2)
         END IF

         ! set up parameters for P_screening
         restore_p_screen = qs_env%x_data(1, 1)%screening_parameter%do_initial_p_screening
         IF (qs_env%x_data(1, 1)%screening_parameter%do_initial_p_screening) THEN
            IF (mp2_env%ri_grad%free_hfx_buffer) THEN
               mp2_env%p_screen = .FALSE.
            ELSE
               mp2_env%p_screen = .TRUE.
            END IF
         END IF
      END IF

      ! Add exx part for RPA
      do_exx = .FALSE.
      IF (qs_env%mp2_env%method == ri_rpa_method_gpw) THEN
         hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
         CALL section_vals_get(hfx_section, explicit=do_exx)
      END IF
      IF (do_exx) THEN
         CALL add_exx_to_rhs(rhs=P_mu_nu, &
                             qs_env=qs_env, &
                             ext_hfx_section=hfx_section, &
                             x_data=mp2_env%ri_rpa%x_data, &
                             recalc_integrals=.FALSE., &
                             do_admm=mp2_env%ri_rpa%do_admm, &
                             do_exx=do_exx, &
                             reuse_hfx=mp2_env%ri_rpa%reuse_hfx)

         focc = 1.0_dp
         IF (nspins == 1) focc = 2.0_dp
         !focc = 0.0_dp
         DO ispin = 1, nspins
            CALL dbcsr_add(P_mu_nu(ispin)%matrix, matrix_ks(ispin)%matrix, 1.0_dp, -1.0_dp)
            CALL copy_dbcsr_to_fm(matrix=P_mu_nu(ispin)%matrix, fm=fm_G_mu_nu)
            CALL parallel_gemm("N", "N", dimen, homo(ispin), dimen, 1.0_dp, &
                               fm_G_mu_nu, mo_coeff_o(ispin), 0.0_dp, fm_back)
            CALL parallel_gemm("T", "N", homo(ispin), virtual(ispin), dimen, focc, &
                               fm_back, mo_coeff_v(ispin), 1.0_dp, L_jb(ispin))
            CALL parallel_gemm("T", "N", homo(ispin), homo(ispin), dimen, -focc, &
                               fm_back, mo_coeff_o(ispin), 1.0_dp, W_mo(ispin))
         END DO
      END IF

      ! Prepare arrays for linres code
      NULLIFY (linres_control)
      ALLOCATE (linres_control)
      linres_control%do_kernel = .TRUE.
      linres_control%lr_triplet = .FALSE.
      linres_control%linres_restart = .FALSE.
      linres_control%max_iter = mp2_env%ri_grad%cphf_max_num_iter
      linres_control%eps = mp2_env%ri_grad%cphf_eps_conv
      linres_control%eps_filter = mp2_env%mp2_gpw%eps_filter
      linres_control%restart_every = 50
      linres_control%preconditioner_type = ot_precond_full_all
      linres_control%energy_gap = 0.02_dp

      NULLIFY (p_env)
      ALLOCATE (p_env)
      CALL p_env_create(p_env, qs_env, p1_option=P_mu_nu, orthogonal_orbitals=.TRUE., linres_control=linres_control)
      CALL set_qs_env(qs_env, linres_control=linres_control)
      CALL p_env_psi0_changed(p_env, qs_env)
      p_env%new_preconditioner = .TRUE.
      CALL p_env_check_i_alloc(p_env, qs_env)
      mp2_env%ri_grad%p_env => p_env

      ! update Lagrangian with the CPHF like update, occ-occ block, first call (recompute hfx integrals if needed)
      transf_type_in = 1
      ! In alpha-beta case, L_bj_alpha has Coulomb and XC alpha-alpha part
      ! and (only) Coulomb alpha-beta part and vice versa.

      ! Complete in closed shell case, alpha-alpha (Coulomb and XC)
      ! part of L_bj(alpha) for open shell

      CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                            mo_coeff_o, &
                            mo_coeff_v, Eigenval, p_env, &
                            P_mo, fm_G_mu_nu, fm_back, transf_type_in, &
                            L_jb, &
                            recalc_hfx_integrals=(.NOT. do_exx .AND. mp2_env%ri_grad%free_hfx_buffer) &
                            .OR. (do_exx .AND. .NOT. mp2_env%ri_rpa%reuse_hfx))

      ! at this point Lagrangian is completed ready to solve the Z-vector equations
      ! P_ia will contain the solution of these equations
      ALLOCATE (P_ia(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=homo(ispin), ncol_global=virtual(ispin))
         CALL cp_fm_create(P_ia(ispin), fm_struct_tmp, name="P_ia")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(P_ia(ispin), 0.0_dp)
      END DO

      CALL solve_z_vector_eq_low(qs_env, mp2_env, homo, virtual, dimen, unit_nr, nspins, &
                                 mo_coeff_o, mo_coeff_v, Eigenval, p_env, &
                                 L_jb, fm_G_mu_nu, fm_back, P_ia)

      ! release fm stuff
      CALL cp_fm_release(fm_G_mu_nu)
      CALL cp_fm_release(L_jb)
      CALL cp_fm_release(mo_coeff_o)
      CALL cp_fm_release(mo_coeff_v)

      DO ispin = 1, nspins
         ! update the MP2-MO density matrix with the occ-virt block
         CALL cp_fm_to_fm_submat(msource=P_ia(ispin), mtarget=P_mo(ispin), &
                                 nrow=homo(ispin), ncol=virtual(ispin), &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=homo(ispin) + 1)
         ! transpose P_MO matrix (easy way to symmetrize)
         CALL cp_fm_set_all(fm_back, 0.0_dp)
         ! P_mo now is ready
         CALL cp_fm_upper_to_full(matrix=P_mo(ispin), work=fm_back)
      END DO
      CALL cp_fm_release(P_ia)

      ! do the final update to the energy weighted matrix W_MO
      DO ispin = 1, nspins
         CALL cp_fm_get_info(matrix=W_mo(ispin), &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            IF (j_global <= homo(ispin)) THEN
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  W_mo(ispin)%local_data(iiB, jjB) = W_mo(ispin)%local_data(iiB, jjB) &
                                                     - P_mo(ispin)%local_data(iiB, jjB)*Eigenval(j_global, ispin)
                  IF (i_global == j_global .AND. nspins == 1) W_mo(ispin)%local_data(iiB, jjB) = &
                     W_mo(ispin)%local_data(iiB, jjB) - 2.0_dp*Eigenval(j_global, ispin)
                  IF (i_global == j_global .AND. nspins == 2) W_mo(ispin)%local_data(iiB, jjB) = &
                     W_mo(ispin)%local_data(iiB, jjB) - Eigenval(j_global, ispin)
               END DO
            ELSE
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  IF (i_global <= homo(ispin)) THEN
                     ! virt-occ
                     W_mo(ispin)%local_data(iiB, jjB) = W_mo(ispin)%local_data(iiB, jjB) &
                                                        - P_mo(ispin)%local_data(iiB, jjB)*Eigenval(i_global, ispin)
                  ELSE
                     ! virt-virt
                     W_mo(ispin)%local_data(iiB, jjB) = W_mo(ispin)%local_data(iiB, jjB) &
                                                        - P_mo(ispin)%local_data(iiB, jjB)*Eigenval(j_global, ispin)
                  END IF
               END DO
            END IF
         END DO
      END DO

      ! create the MP2 energy weighted density matrix
      NULLIFY (p_env%w1)
      CALL dbcsr_allocate_matrix_set(p_env%w1, 1)
      ALLOCATE (p_env%w1(1)%matrix)
      CALL dbcsr_copy(p_env%w1(1)%matrix, matrix_s(1)%matrix, &
                      name="W MATRIX MP2")
      CALL dbcsr_set(p_env%w1(1)%matrix, 0.0_dp)

      ! backtnsform the collected parts of the energy-weighted density matrix into AO basis
      DO ispin = 1, nspins
         CALL parallel_gemm('N', 'N', dimen, dimen, dimen, 1.0_dp, &
                            mo_coeff(ispin), W_mo(ispin), 0.0_dp, fm_back)
         CALL cp_dbcsr_plus_fm_fm_t(p_env%w1(1)%matrix, fm_back, mo_coeff(ispin), dimen, 1.0_dp, .TRUE., 1)
      END DO
      CALL cp_fm_release(W_mo)

      CALL qs_rho_get(p_env%rho1, rho_ao=rho1_ao)

      DO ispin = 1, nspins
         CALL dbcsr_set(p_env%p1(ispin)%matrix, 0.0_dp)

         CALL parallel_gemm('N', 'N', dimen, dimen, dimen, 1.0_dp, &
                            mo_coeff(ispin), P_mo(ispin), 0.0_dp, fm_back)
         CALL cp_dbcsr_plus_fm_fm_t(p_env%p1(ispin)%matrix, fm_back, mo_coeff(ispin), dimen, 1.0_dp, .TRUE.)

         CALL dbcsr_copy(rho1_ao(ispin)%matrix, p_env%p1(ispin)%matrix)
      END DO
      CALL cp_fm_release(P_mo)
      CALL cp_fm_release(fm_back)

      CALL p_env_update_rho(p_env, qs_env)

      ! create mp2 DBCSR density
      CALL dbcsr_allocate_matrix_set(matrix_p_mp2, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_p_mp2(ispin)%matrix)
         CALL dbcsr_copy(matrix_p_mp2(ispin)%matrix, p_env%p1(ispin)%matrix, &
                         name="P MATRIX MP2")
      END DO

      IF (dft_control%do_admm) THEN
         CALL get_admm_env(qs_env%admm_env, rho_aux_fit=rho_aux_fit)
         CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux_fit)

         ! create mp2 DBCSR density in auxiliary basis
         CALL dbcsr_allocate_matrix_set(matrix_p_mp2_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (matrix_p_mp2_admm(ispin)%matrix)
            CALL dbcsr_copy(matrix_p_mp2_admm(ispin)%matrix, p_env%p1_admm(ispin)%matrix, &
                            name="P MATRIX MP2 ADMM")
         END DO
      END IF

      CALL set_ks_env(ks_env, matrix_p_mp2=matrix_p_mp2, matrix_p_mp2_admm=matrix_p_mp2_admm)

      ! We will need one more hfx calculation for HF gradient part
      mp2_env%not_last_hfx = .FALSE.
      mp2_env%p_screen = restore_p_screen

      CALL timestop(handle)

   END SUBROUTINE solve_z_vector_eq

! **************************************************************************************************
!> \brief Here we performe the CPHF like update using GPW,
!>        transf_type_in  defines the type of transformation for the matrix in input
!>        transf_type_in = 1 -> occ-occ back transformation
!>        transf_type_in = 2 -> virt-virt back transformation
!>        transf_type_in = 3 -> occ-virt back transformation including the
!>                              eigenvalues energy differences for the diagonal elements
!> \param qs_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param nspins ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param Eigenval ...
!> \param p_env ...
!> \param fm_mo ...
!> \param fm_ao ...
!> \param fm_back ...
!> \param transf_type_in ...
!> \param fm_mo_out ...
!> \param recalc_hfx_integrals ...
!> \author Mauro Del Ben, Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                               mo_coeff_o, mo_coeff_v, Eigenval, p_env, &
                               fm_mo, fm_ao, fm_back, transf_type_in, &
                               fm_mo_out, recalc_hfx_integrals)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      INTEGER, INTENT(IN)                                :: nspins, dimen
      INTEGER, DIMENSION(nspins), INTENT(IN)             :: virtual, homo
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), DIMENSION(dimen, nspins), &
         INTENT(IN)                                      :: Eigenval
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: fm_mo
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_ao, fm_back
      INTEGER, INTENT(IN)                                :: transf_type_in
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: fm_mo_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: recalc_hfx_integrals

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cphf_like_update'

      INTEGER                                            :: handle, i_global, iiB, ispin, j_global, &
                                                            jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho1_ao

      CALL timeset(routineN, handle)

      CALL qs_rho_get(p_env%rho1, rho_ao=rho1_ao)

      ! Determine the first-order density matrices in AO basis
      DO ispin = 1, nspins
         CALL dbcsr_set(p_env%p1(ispin)%matrix, 0.0_dp)

         ASSOCIATE (mat_in => fm_mo(ispin))

            ! perform back transformation
            SELECT CASE (transf_type_in)
            CASE (1)
               ! occ-occ block
               CALL parallel_gemm('N', 'N', dimen, homo(ispin), homo(ispin), 1.0_dp, &
                                  mo_coeff_o(ispin), mat_in, 0.0_dp, fm_back)
               CALL cp_dbcsr_plus_fm_fm_t(p_env%p1(ispin)%matrix, fm_back, mo_coeff_o(ispin), homo(ispin), 1.0_dp, .TRUE.)
               ! virt-virt block
               CALL parallel_gemm('N', 'N', dimen, virtual(ispin), virtual(ispin), 1.0_dp, &
                                  mo_coeff_v(ispin), mat_in, 0.0_dp, fm_back, &
                                  b_first_col=homo(ispin) + 1, &
                                  b_first_row=homo(ispin) + 1)
               CALL cp_dbcsr_plus_fm_fm_t(p_env%p1(ispin)%matrix, fm_back, mo_coeff_v(ispin), virtual(ispin), 1.0_dp, .TRUE.)

            CASE (3)
               ! virt-occ blocks
               CALL parallel_gemm('N', 'N', dimen, virtual(ispin), homo(ispin), 1.0_dp, &
                                  mo_coeff_o(ispin), mat_in, 0.0_dp, fm_back)
               CALL cp_dbcsr_plus_fm_fm_t(p_env%p1(ispin)%matrix, fm_back, mo_coeff_v(ispin), virtual(ispin), 1.0_dp, .TRUE.)
               ! and symmetrize (here again multiply instead of transposing)
               CALL parallel_gemm('N', 'T', dimen, homo(ispin), virtual(ispin), 1.0_dp, &
                                  mo_coeff_v(ispin), mat_in, 0.0_dp, fm_back)
               CALL cp_dbcsr_plus_fm_fm_t(p_env%p1(ispin)%matrix, fm_back, mo_coeff_o(ispin), homo(ispin), 1.0_dp, .TRUE.)

            CASE DEFAULT
               ! nothing
            END SELECT
         END ASSOCIATE

         CALL dbcsr_copy(rho1_ao(ispin)%matrix, p_env%p1(ispin)%matrix)
      END DO

      CALL p_env_update_rho(p_env, qs_env)

      IF (transf_type_in == 3) THEN
         DO ispin = 1, nspins

            ! scale for the orbital energy differences for the diagonal elements
            CALL cp_fm_get_info(matrix=fm_mo_out(ispin), &
                                nrow_local=nrow_local, &
                                ncol_local=ncol_local, &
                                row_indices=row_indices, &
                                col_indices=col_indices)
            DO jjB = 1, ncol_local
               j_global = col_indices(jjB)
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  fm_mo_out(ispin)%local_data(iiB, jjB) = fm_mo(ispin)%local_data(iiB, jjB)* &
                                                          (Eigenval(j_global + homo(ispin), ispin) - Eigenval(i_global, ispin))
               END DO
            END DO
         END DO
      END IF

      CALL apply_2nd_order_kernel(qs_env, p_env, recalc_hfx_integrals)

      DO ispin = 1, nspins
         ! copy back to fm
         CALL cp_fm_set_all(fm_ao, 0.0_dp)
         CALL copy_dbcsr_to_fm(matrix=p_env%kpp1(ispin)%matrix, fm=fm_ao)
         CALL cp_fm_set_all(fm_back, 0.0_dp)
         CALL cp_fm_upper_to_full(fm_ao, fm_back)

         ASSOCIATE (mat_out => fm_mo_out(ispin))

            ! transform to MO basis, here we always sum the result into the input matrix

            ! occ-virt block
            CALL parallel_gemm('T', 'N', homo(ispin), dimen, dimen, 1.0_dp, &
                               mo_coeff_o(ispin), fm_ao, 0.0_dp, fm_back)
            CALL parallel_gemm('N', 'N', homo(ispin), virtual(ispin), dimen, 1.0_dp, &
                               fm_back, mo_coeff_v(ispin), 1.0_dp, mat_out)
         END ASSOCIATE
      END DO

      CALL timestop(handle)

   END SUBROUTINE cphf_like_update

! **************************************************************************************************
!> \brief Low level subroutine for the iterative solution of a large
!>        system of linear equation
!> \param qs_env ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param unit_nr ...
!> \param nspins ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param Eigenval ...
!> \param p_env ...
!> \param L_jb ...
!> \param fm_G_mu_nu ...
!> \param fm_back ...
!> \param P_ia ...
!> \author Mauro Del Ben, Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE solve_z_vector_eq_low(qs_env, mp2_env, homo, virtual, dimen, unit_nr, nspins, &
                                    mo_coeff_o, mo_coeff_v, Eigenval, p_env, &
                                    L_jb, fm_G_mu_nu, fm_back, P_ia)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: nspins, unit_nr, dimen
      INTEGER, DIMENSION(nspins), INTENT(IN)             :: virtual, homo
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), DIMENSION(dimen, nspins), &
         INTENT(IN)                                      :: Eigenval
      TYPE(qs_p_env_type), INTENT(IN), POINTER           :: p_env
      TYPE(cp_fm_type), DIMENSION(dimen), INTENT(IN)     :: L_jb
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_G_mu_nu, fm_back
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: P_ia

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_z_vector_eq_low'

      INTEGER                                            :: handle, i_global, iiB, ispin, j_global, &
                                                            jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: precond

      CALL timeset(routineN, handle)

      ! Pople method
      ! change sign to L_jb
      DO ispin = 1, nspins
         L_jb(ispin)%local_data(:, :) = -L_jb(ispin)%local_data(:, :)
      END DO

      ! create fm structure
      ALLOCATE (precond(nspins))
      DO ispin = 1, nspins
         ! create preconditioner (for now only orbital energy differences)
         CALL cp_fm_create(precond(ispin), P_ia(ispin)%matrix_struct, name="precond")
         CALL cp_fm_set_all(precond(ispin), 1.0_dp)
         CALL cp_fm_get_info(matrix=precond(ispin), &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               precond(ispin)%local_data(iiB, jjB) = Eigenval(j_global + homo(ispin), ispin) - Eigenval(i_global, ispin)
            END DO
         END DO
      END DO

      DO ispin = 1, nspins
         precond(ispin)%local_data(:, :) = 1.0_dp/precond(ispin)%local_data(:, :)
      END DO

      SELECT CASE (mp2_env%ri_grad%z_solver_method)
      CASE (z_solver_pople)
         CALL solve_z_vector_pople(qs_env, mp2_env, homo, virtual, dimen, unit_nr, nspins, &
                                   mo_coeff_o, mo_coeff_v, Eigenval, p_env, &
                                   L_jb, fm_G_mu_nu, fm_back, P_ia, precond)
      CASE (z_solver_cg, z_solver_richardson, z_solver_sd)
         CALL solve_z_vector_cg(qs_env, mp2_env, homo, virtual, dimen, unit_nr, nspins, &
                                mo_coeff_o, mo_coeff_v, Eigenval, p_env, &
                                L_jb, fm_G_mu_nu, fm_back, P_ia, precond)
      CASE DEFAULT
         CPABORT("Unknown solver")
      END SELECT

      CALL cp_fm_release(precond)

      CALL timestop(handle)

   END SUBROUTINE solve_z_vector_eq_low

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param unit_nr ...
!> \param nspins ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param Eigenval ...
!> \param p_env ...
!> \param L_jb ...
!> \param fm_G_mu_nu ...
!> \param fm_back ...
!> \param P_ia ...
!> \param precond ...
! **************************************************************************************************
   SUBROUTINE solve_z_vector_pople(qs_env, mp2_env, homo, virtual, dimen, unit_nr, nspins, &
                                   mo_coeff_o, mo_coeff_v, Eigenval, p_env, &
                                   L_jb, fm_G_mu_nu, fm_back, P_ia, precond)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: nspins, unit_nr, dimen
      INTEGER, DIMENSION(nspins), INTENT(IN)             :: virtual, homo
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), DIMENSION(dimen, nspins), &
         INTENT(IN)                                      :: Eigenval
      TYPE(qs_p_env_type), INTENT(IN), POINTER           :: p_env
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: L_jb
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_G_mu_nu, fm_back
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: P_ia, precond

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_z_vector_pople'

      INTEGER                                            :: cycle_counter, handle, iiB, iiter, &
                                                            ispin, max_num_iter, transf_type_in
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: conv, eps_conv, scale_cphf, t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: proj_bi_xj, temp_vals, x_norm, xi_b
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: A_small, b_small, xi_Axi
      TYPE(cp_fm_struct_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: fm_struct_tmp
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: b_i, residual
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: Ax, xn

      CALL timeset(routineN, handle)

      eps_conv = mp2_env%ri_grad%cphf_eps_conv

      IF (SQRT(accurate_dot_product_spin(L_jb, L_jb)) >= eps_conv) THEN

         max_num_iter = mp2_env%ri_grad%cphf_max_num_iter
         scale_cphf = mp2_env%ri_grad%scale_step_size

         ! set the transformation type (equal for all methods all updates)
         transf_type_in = 3

         ! set convergence flag
         converged = .FALSE.

         ALLOCATE (fm_struct_tmp(nspins), b_i(nspins), residual(nspins))
         DO ispin = 1, nspins
            fm_struct_tmp(ispin)%struct => P_ia(ispin)%matrix_struct

            CALL cp_fm_create(b_i(ispin), fm_struct_tmp(ispin)%struct, name="b_i")
            CALL cp_fm_set_all(b_i(ispin), 0.0_dp)
            b_i(ispin)%local_data(:, :) = precond(ispin)%local_data(:, :)*L_jb(ispin)%local_data(:, :)

            ! create the residual vector (r), we check convergence on the norm of
            ! this vector r=(Ax-b)
            CALL cp_fm_create(residual(ispin), fm_struct_tmp(ispin)%struct, name="residual")
            CALL cp_fm_set_all(residual(ispin), 0.0_dp)
         END DO

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T3,A)') "MP2_CPHF| Iterative solution of Z-Vector equations (Pople's method)"
            WRITE (unit_nr, '(T3,A,T45,ES8.1)') 'MP2_CPHF| Convergence threshold:', eps_conv
            WRITE (unit_nr, '(T3,A,T45,I8)') 'MP2_CPHF| Maximum number of iterations: ', max_num_iter
            WRITE (unit_nr, '(T3,A,T45,ES8.1)') 'MP2_CPHF| Scaling of initial guess: ', scale_cphf
            WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
            WRITE (unit_nr, '(T4,A,T15,A,T33,A)') 'Step', 'Time', 'Convergence'
            WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
         END IF

         ALLOCATE (xn(nspins, max_num_iter))
         ALLOCATE (Ax(nspins, max_num_iter))
         ALLOCATE (x_norm(max_num_iter))
         ALLOCATE (xi_b(max_num_iter))
         ALLOCATE (xi_Axi(max_num_iter, 0:max_num_iter))
         x_norm = 0.0_dp
         xi_b = 0.0_dp
         xi_Axi = 0.0_dp

         cycle_counter = 0
         DO iiter = 1, max_num_iter
            cycle_counter = cycle_counter + 1

            t1 = m_walltime()

            ! create and update x_i (orthogonalization with previous vectors)
            DO ispin = 1, nspins
               CALL cp_fm_create(xn(ispin, iiter), fm_struct_tmp(ispin)%struct, name="xi")
               CALL cp_fm_set_all(xn(ispin, iiter), 0.0_dp)
            END DO

            ALLOCATE (proj_bi_xj(iiter - 1))
            proj_bi_xj = 0.0_dp
            ! first compute the projection of the actual b_i into all previous x_i
            ! already scaled with the norm of each x_i
            DO iiB = 1, iiter - 1
               proj_bi_xj(iiB) = proj_bi_xj(iiB) + accurate_dot_product_spin(b_i, xn(:, iiB))/x_norm(iiB)
            END DO

            ! update actual x_i
            DO ispin = 1, nspins
               xn(ispin, iiter)%local_data(:, :) = scale_cphf*b_i(ispin)%local_data(:, :)
               DO iiB = 1, iiter - 1
                  xn(ispin, iiter)%local_data(:, :) = xn(ispin, iiter)%local_data(:, :) - &
                                                      xn(ispin, iiB)%local_data(:, :)*proj_bi_xj(iiB)
               END DO
            END DO
            DEALLOCATE (proj_bi_xj)

            ! create Ax(iiter) that will store the matrix vector product for this cycle
            DO ispin = 1, nspins
               CALL cp_fm_create(Ax(ispin, iiter), fm_struct_tmp(ispin)%struct, name="Ai")
               CALL cp_fm_set_all(Ax(ispin, iiter), 0.0_dp)
            END DO

            CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                                  mo_coeff_o, &
                                  mo_coeff_v, Eigenval, p_env, &
                                  xn(:, iiter), fm_G_mu_nu, fm_back, transf_type_in, &
                                  Ax(:, iiter))

            ! in order to reduce the number of  parallel sums here we
            ! cluster all necessary scalar products into a single vector
            ! temp_vals contains:
            ! 1:iiter -> <Ax_i|x_j>
            ! iiter+1 -> <x_i|b>
            ! iiter+2 -> <x_i|x_i>

            ALLOCATE (temp_vals(iiter + 2))
            temp_vals = 0.0_dp
            ! <Ax_i|x_j>
            DO iiB = 1, iiter
               temp_vals(iiB) = temp_vals(iiB) + accurate_dot_product_spin(Ax(:, iiter), xn(:, iiB))
            END DO
            ! <x_i|b>
            temp_vals(iiter + 1) = temp_vals(iiter + 1) + accurate_dot_product_spin(xn(:, iiter), L_jb)
            ! norm
            temp_vals(iiter + 2) = temp_vals(iiter + 2) + accurate_dot_product_spin(xn(:, iiter), xn(:, iiter))
            ! update <Ax_i|x_j>,  <x_i|b> and norm <x_i|x_i>
            xi_Axi(iiter, 1:iiter) = temp_vals(1:iiter)
            xi_Axi(1:iiter, iiter) = temp_vals(1:iiter)
            xi_b(iiter) = temp_vals(iiter + 1)
            x_norm(iiter) = temp_vals(iiter + 2)
            DEALLOCATE (temp_vals)

            ! solve reduced system
            IF (ALLOCATED(A_small)) DEALLOCATE (A_small)
            IF (ALLOCATED(b_small)) DEALLOCATE (b_small)
            ALLOCATE (A_small(iiter, iiter))
            ALLOCATE (b_small(iiter, 1))
            A_small(1:iiter, 1:iiter) = xi_Axi(1:iiter, 1:iiter)
            b_small(1:iiter, 1) = xi_b(1:iiter)

            CALL solve_system(matrix=A_small, mysize=iiter, eigenvectors=b_small)

            ! check for convergence
            DO ispin = 1, nspins
               CALL cp_fm_set_all(residual(ispin), 0.0_dp)
               DO iiB = 1, iiter
                  residual(ispin)%local_data(:, :) = &
                     residual(ispin)%local_data(:, :) + &
                     b_small(iiB, 1)*Ax(ispin, iiB)%local_data(:, :)
               END DO

               residual(ispin)%local_data(:, :) = &
                  residual(ispin)%local_data(:, :) - &
                  L_jb(ispin)%local_data(:, :)
            END DO

            conv = SQRT(accurate_dot_product_spin(residual, residual))

            t2 = m_walltime()

            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T3,I5,T13,F6.1,11X,F14.8)') iiter, t2 - t1, conv
               CALL m_flush(unit_nr)
            END IF

            IF (conv <= eps_conv) THEN
               converged = .TRUE.
               EXIT
            END IF

            ! update b_i for the next round
            DO ispin = 1, nspins
               b_i(ispin)%local_data(:, :) = b_i(ispin)%local_data(:, :) &
                                             + precond(ispin)%local_data(:, :) &
                                             *Ax(ispin, iiter)%local_data(:, :)
            END DO

            scale_cphf = 1.0_dp

         END DO

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
            IF (converged) THEN
               WRITE (unit_nr, '(T3,A,I5,A)') 'Z-Vector equations converged in', cycle_counter, ' steps'
            ELSE
               WRITE (unit_nr, '(T3,A,I5,A)') 'Z-Vector equations NOT converged in', cycle_counter, ' steps'
            END IF
         END IF

         ! store solution into P_ia
         DO iiter = 1, cycle_counter
            DO ispin = 1, nspins
               P_ia(ispin)%local_data(:, :) = P_ia(ispin)%local_data(:, :) + &
                                              b_small(iiter, 1)*xn(ispin, iiter)%local_data(:, :)
            END DO
         END DO

         ! Release arrays
         DEALLOCATE (x_norm)
         DEALLOCATE (xi_b)
         DEALLOCATE (xi_Axi)

         CALL cp_fm_release(b_i)
         CALL cp_fm_release(residual)
         CALL cp_fm_release(Ax)
         CALL cp_fm_release(xn)

      ELSE
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
            WRITE (unit_nr, '(T3,A)') 'Residual smaller than EPS_CONV. Skip solution of Z-vector equation.'
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE solve_z_vector_pople

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param unit_nr ...
!> \param nspins ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param Eigenval ...
!> \param p_env ...
!> \param L_jb ...
!> \param fm_G_mu_nu ...
!> \param fm_back ...
!> \param P_ia ...
!> \param precond ...
! **************************************************************************************************
   SUBROUTINE solve_z_vector_cg(qs_env, mp2_env, homo, virtual, dimen, unit_nr, nspins, &
                                mo_coeff_o, mo_coeff_v, Eigenval, p_env, &
                                L_jb, fm_G_mu_nu, fm_back, P_ia, precond)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, INTENT(IN)                                :: nspins, unit_nr, dimen
      INTEGER, DIMENSION(nspins), INTENT(IN)             :: virtual, homo
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), DIMENSION(dimen, nspins), &
         INTENT(IN)                                      :: Eigenval
      TYPE(qs_p_env_type), INTENT(IN), POINTER           :: p_env
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: L_jb
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_G_mu_nu, fm_back
      TYPE(cp_fm_type), DIMENSION(nspins), INTENT(IN)    :: P_ia, precond

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_z_vector_cg'

      INTEGER :: cycles_passed, handle, iiter, ispin, max_num_iter, restart_counter, &
         restart_every, transf_type_in, z_solver_method
      LOGICAL                                            :: converged, do_restart
      REAL(KIND=dp) :: eps_conv, norm_residual, norm_residual_old, &
         residual_dot_diff_search_vec_new, residual_dot_diff_search_vec_old, &
         residual_dot_search_vec, residual_new_dot_diff_search_vec_old, scale_result, &
         scale_search, scale_step_size, search_vec_dot_A_search_vec, t1, t2
      TYPE(cp_fm_struct_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: fm_struct_tmp
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: A_dot_search_vector, diff_search_vector, &
                                                            residual, search_vector

      CALL timeset(routineN, handle)

      max_num_iter = mp2_env%ri_grad%cphf_max_num_iter
      eps_conv = mp2_env%ri_grad%cphf_eps_conv
      z_solver_method = mp2_env%ri_grad%z_solver_method
      restart_every = mp2_env%ri_grad%cphf_restart
      scale_step_size = mp2_env%ri_grad%scale_step_size
      transf_type_in = 3

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         SELECT CASE (z_solver_method)
         CASE (z_solver_cg)
            IF (mp2_env%ri_grad%polak_ribiere) THEN
               WRITE (unit_nr, '(T3,A)') 'MP2_CPHF| Iterative solution of Z-Vector equations (CG with Polak-Ribiere step)'
            ELSE
               WRITE (unit_nr, '(T3,A)') 'MP2_CPHF| Iterative solution of Z-Vector equations (CG with Fletcher-Reeves step)'
            END IF
         CASE (z_solver_richardson)
            WRITE (unit_nr, '(T3,A)') 'MP2_CPHF| Iterative solution of Z-Vector equations (Richardson method)'
         CASE (z_solver_sd)
            WRITE (unit_nr, '(T3,A)') 'MP2_CPHF| Iterative solution of Z-Vector equations (Steepest Descent method)'
         CASE DEFAULT
            CPABORT("Unknown solver")
         END SELECT
         WRITE (unit_nr, '(T3,A,T45,ES8.1)') 'MP2_CPHF| Convergence threshold:', eps_conv
         WRITE (unit_nr, '(T3,A,T45,I8)') 'MP2_CPHF| Maximum number of iterations: ', max_num_iter
         WRITE (unit_nr, '(T3,A,T45,I8)') 'MP2_CPHF| Number of steps for restart: ', restart_every
         WRITE (unit_nr, '(T3, A)') 'MP2_CPHF| Restart after no decrease'
         WRITE (unit_nr, '(T3,A,T45,ES8.1)') 'MP2_CPHF| Scaling factor of each step: ', scale_step_size
         WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
         WRITE (unit_nr, '(T4,A,T13,A,T28,A,T43,A)') 'Step', 'Restart', 'Time', 'Convergence'
         WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
      END IF

      ALLOCATE (fm_struct_tmp(nspins), residual(nspins), diff_search_vector(nspins), &
                search_vector(nspins), A_dot_search_vector(nspins))
      DO ispin = 1, nspins
         fm_struct_tmp(ispin)%struct => P_ia(ispin)%matrix_struct

         CALL cp_fm_create(residual(ispin), fm_struct_tmp(ispin)%struct, name="residual")
         CALL cp_fm_set_all(residual(ispin), 0.0_dp)

         CALL cp_fm_create(diff_search_vector(ispin), fm_struct_tmp(ispin)%struct, name="difference search vector")
         CALL cp_fm_set_all(diff_search_vector(ispin), 0.0_dp)

         CALL cp_fm_create(search_vector(ispin), fm_struct_tmp(ispin)%struct, name="search vector")
         CALL cp_fm_set_all(search_vector(ispin), 0.0_dp)

         CALL cp_fm_create(A_dot_search_vector(ispin), fm_struct_tmp(ispin)%struct, name="A times search vector")
         CALL cp_fm_set_all(A_dot_search_vector(ispin), 0.0_dp)
      END DO

      converged = .FALSE.
      cycles_passed = max_num_iter
      ! By that, we enforce the setup of the matrices
      do_restart = .TRUE.

      t1 = m_walltime()

      DO iiter = 1, max_num_iter

         ! During the first iteration, P_ia=0 such that the application of the 2nd order matrix is zero
         IF (do_restart) THEN
            ! We do not consider the first step to be a restart
            ! Do not recalculate residual if it is already enforced to save FLOPs
            IF (.NOT. mp2_env%ri_grad%recalc_residual .OR. (iiter == 1)) THEN
               IF (iiter > 1) THEN
                  CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                                        mo_coeff_o, &
                                        mo_coeff_v, Eigenval, p_env, &
                                        P_ia, fm_G_mu_nu, fm_back, transf_type_in, &
                                        residual)
               ELSE
                  do_restart = .FALSE.

                  DO ispin = 1, nspins
                     CALL cp_fm_set_all(residual(ispin), 0.0_dp)
                  END DO
               END IF

               DO ispin = 1, nspins
                  residual(ispin)%local_data(:, :) = L_jb(ispin)%local_data(:, :) &
                                                     - residual(ispin)%local_data(:, :)
               END DO
            END IF

            DO ispin = 1, nspins
               diff_search_vector(ispin)%local_data(:, :) = &
                  precond(ispin)%local_data(:, :)*residual(ispin)%local_data(:, :)
               search_vector(ispin)%local_data(:, :) = diff_search_vector(ispin)%local_data(:, :)
            END DO

            restart_counter = 1
         END IF

         norm_residual_old = SQRT(accurate_dot_product_spin(residual, residual))

         residual_dot_diff_search_vec_old = accurate_dot_product_spin(residual, diff_search_vector)

         CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                               mo_coeff_o, &
                               mo_coeff_v, Eigenval, p_env, &
                               search_vector, fm_G_mu_nu, fm_back, transf_type_in, &
                               A_dot_search_vector)

         IF (z_solver_method /= z_solver_richardson) THEN
            search_vec_dot_A_search_vec = accurate_dot_product_spin(search_vector, A_dot_search_vector)

            IF (z_solver_method == z_solver_cg) THEN
               scale_result = residual_dot_diff_search_vec_old/search_vec_dot_A_search_vec
            ELSE
               residual_dot_search_vec = accurate_dot_product_spin(residual, search_vector)
               scale_result = residual_dot_search_vec/search_vec_dot_A_search_vec
            END IF

            scale_result = scale_result*scale_step_size

         ELSE

            scale_result = scale_step_size

         END IF

         DO ispin = 1, nspins
            P_ia(ispin)%local_data(:, :) = P_ia(ispin)%local_data(:, :) &
                                           + scale_result*search_vector(ispin)%local_data(:, :)
         END DO

         IF (.NOT. mp2_env%ri_grad%recalc_residual) THEN

            DO ispin = 1, nspins
               residual(ispin)%local_data(:, :) = residual(ispin)%local_data(:, :) &
                                                  - scale_result*A_dot_search_vector(ispin)%local_data(:, :)
            END DO
         ELSE
            CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                                  mo_coeff_o, &
                                  mo_coeff_v, Eigenval, p_env, &
                                  P_ia, fm_G_mu_nu, fm_back, transf_type_in, &
                                  residual)

            DO ispin = 1, nspins
               residual(ispin)%local_data(:, :) = L_jb(ispin)%local_data(:, :) - residual(ispin)%local_data(:, :)
            END DO
         END IF

         norm_residual = SQRT(accurate_dot_product_spin(residual, residual))

         t2 = m_walltime()

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,I4,T16,L1,T26,F6.1,8X,F14.8)') iiter, do_restart, t2 - t1, norm_residual
            CALL m_flush(unit_nr)
         END IF

         IF (norm_residual <= eps_conv) THEN
            converged = .TRUE.
            cycles_passed = iiter
            EXIT
         END IF

         t1 = m_walltime()

         IF (z_solver_method == z_solver_richardson) THEN
            DO ispin = 1, nspins
               search_vector(ispin)%local_data(:, :) = &
                  scale_step_size*precond(ispin)%local_data(:, :)*residual(ispin)%local_data(:, :)
            END DO
         ELSE IF (z_solver_method == z_solver_sd) THEN
            DO ispin = 1, nspins
               search_vector(ispin)%local_data(:, :) = &
                  precond(ispin)%local_data(:, :)*residual(ispin)%local_data(:, :)
            END DO
         ELSE
            IF (mp2_env%ri_grad%polak_ribiere) &
               residual_new_dot_diff_search_vec_old = accurate_dot_product_spin(residual, diff_search_vector)

            DO ispin = 1, nspins
               diff_search_vector(ispin)%local_data(:, :) = &
                  precond(ispin)%local_data(:, :)*residual(ispin)%local_data(:, :)
            END DO

            residual_dot_diff_search_vec_new = accurate_dot_product_spin(residual, diff_search_vector)

            scale_search = residual_dot_diff_search_vec_new/residual_dot_diff_search_vec_old
            IF (mp2_env%ri_grad%polak_ribiere) scale_search = &
               scale_search - residual_new_dot_diff_search_vec_old/residual_dot_diff_search_vec_old

            DO ispin = 1, nspins
               search_vector(ispin)%local_data(:, :) = scale_search*search_vector(ispin)%local_data(:, :) &
                                                       + diff_search_vector(ispin)%local_data(:, :)
            END DO

            ! Make new to old
            residual_dot_diff_search_vec_old = residual_dot_diff_search_vec_new
         END IF

         ! Check whether the residual decrease or restart is enforced and ask for restart
         do_restart = (norm_residual >= norm_residual_old .OR. (MOD(restart_counter, restart_every) == 0))

         restart_counter = restart_counter + 1
         norm_residual_old = norm_residual

      END DO

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
         IF (converged) THEN
            WRITE (unit_nr, '(T3,A,I5,A)') 'Z-Vector equations converged in', cycles_passed, ' steps'
         ELSE
            WRITE (unit_nr, '(T3,A,I5,A)') 'Z-Vector equations NOT converged in', max_num_iter, ' steps'
         END IF
      END IF

      DEALLOCATE (fm_struct_tmp)
      CALL cp_fm_release(residual)
      CALL cp_fm_release(diff_search_vector)
      CALL cp_fm_release(search_vector)
      CALL cp_fm_release(A_dot_search_vector)

      CALL timestop(handle)

   END SUBROUTINE solve_z_vector_cg

! **************************************************************************************************
!> \brief ...
!> \param matrix1 ...
!> \param matrix2 ...
!> \return ...
! **************************************************************************************************
   FUNCTION accurate_dot_product_spin(matrix1, matrix2) RESULT(dotproduct)
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: matrix1, matrix2
      REAL(KIND=dp)                                      :: dotproduct

      INTEGER                                            :: ispin

      dotproduct = 0.0_dp
      DO ispin = 1, SIZE(matrix1)
         dotproduct = dotproduct + accurate_dot_product(matrix1(ispin)%local_data, matrix2(ispin)%local_data)
      END DO
      CALL matrix1(1)%matrix_struct%para_env%sum(dotproduct)

   END FUNCTION accurate_dot_product_spin

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE update_mp2_forces(qs_env)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_mp2_forces'

      INTEGER                                            :: alpha, beta, handle, idir, iounit, &
                                                            ispin, nocc, nspins
      INTEGER, DIMENSION(3)                              :: comp
      LOGICAL                                            :: do_exx, do_hfx, use_virial
      REAL(KIND=dp)                                      :: e_dummy, e_hartree, e_xc, ehartree, exc
      REAL(KIND=dp), DIMENSION(3)                        :: deb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_stress, pv_virial
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:), &
         TARGET                                          :: matrix_ks_aux
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_p_mp2, &
                                                            matrix_p_mp2_admm, matrix_s, rho1, &
                                                            rho_ao, rho_ao_aux, scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp, scrm_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: pot_g, pot_r, rho_tot_g, temp_pw_g, &
                                                            vh_rspace, vhxc_rspace
      TYPE(pw_type), ALLOCATABLE, DIMENSION(:)           :: dvg
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_g, rho_mp2_g, rho_mp2_r, &
                                                            rho_mp2_r_aux, rho_r, tau_mp2_r, &
                                                            vadmm_rspace, vtau_rspace, vxc_rspace
      TYPE(pw_type), POINTER                             :: rho_core
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux
      TYPE(section_vals_type), POINTER                   :: hfx_section, hfx_sections, input, &
                                                            xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (input, pw_env, matrix_s, rho, energy, force, virial, &
               matrix_p_mp2, matrix_p_mp2_admm, matrix_ks, rho_core)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env, &
                      input=input, &
                      mos=mos, &
                      para_env=para_env, &
                      matrix_s=matrix_s, &
                      matrix_ks=matrix_ks, &
                      matrix_p_mp2=matrix_p_mp2, &
                      matrix_p_mp2_admm=matrix_p_mp2_admm, &
                      rho=rho, &
                      cell=cell, &
                      force=force, &
                      virial=virial, &
                      sab_orb=sab_orb, &
                      energy=energy, &
                      rho_core=rho_core, &
                      x_data=x_data)

      logger => cp_get_default_logger()
      iounit = cp_print_key_unit_nr(logger, input, "DFT%XC%WF_CORRELATION%PRINT", &
                                    extension=".mp2Log")

      do_exx = .FALSE.
      IF (qs_env%mp2_env%method == ri_rpa_method_gpw) THEN
         hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
         CALL section_vals_get(hfx_section, explicit=do_exx)
      END IF

      p_env => qs_env%mp2_env%ri_grad%p_env

      CALL qs_rho_get(rho, rho_ao=rho_ao, rho_ao_kp=rho_ao_kp, rho_r=rho_r, rho_g=rho_g)
      nspins = SIZE(rho_ao)

      ! check if we have to calculate the virial
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (use_virial) virial%pv_calculate = .TRUE.

      CALL zero_qs_force(force)
      IF (use_virial) CALL zero_virial(virial, .FALSE.)

      DO ispin = 1, nspins
         CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, 1.0_dp)
      END DO

      IF (nspins == 2) THEN
         CALL dbcsr_add(rho_ao(1)%matrix, rho_ao(2)%matrix, 1.0_dp, 1.0_dp)
      END IF

      ! Kinetic energy matrix
      NULLIFY (scrm)
      IF (debug_forces) THEN
         deb(1:3) = force(1)%kinetic(1:3, 1)
         IF (use_virial) e_dummy = third_tr(virial%pv_virial)
      END IF
      CALL build_kinetic_matrix(ks_env, matrix_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=rho_ao(1)%matrix)
      IF (debug_forces) THEN
         deb(1:3) = force(1)%kinetic(1:3, 1) - deb(1:3)
         CALL para_env%sum(deb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dT       ", deb
         IF (use_virial) THEN
            e_dummy = third_tr(virial%pv_virial) - e_dummy
            CALL para_env%sum(e_dummy)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: T          ", e_dummy
         END IF
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (nspins == 2) THEN
         CALL dbcsr_add(rho_ao(1)%matrix, rho_ao(2)%matrix, 1.0_dp, -1.0_dp)
      END IF

      ! Add pseudo potential terms
      scrm_kp(1:nspins, 1:1) => matrix_ks(1:nspins)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, sac_ppl=sac_ppl, sap_ppnl=sap_ppnl)
      IF (ASSOCIATED(sac_ppl)) THEN
         IF (debug_forces) THEN
            deb(1:3) = force(1)%gth_ppl(1:3, 1)
            IF (use_virial) e_dummy = third_tr(virial%pv_virial)
         END IF
         CALL build_core_ppl(scrm_kp, rho_ao_kp, force, &
                             virial, .TRUE., use_virial, 1, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages=1, basis_type="ORB")
         IF (debug_forces) THEN
            deb(1:3) = force(1)%gth_ppl(1:3, 1) - deb(1:3)
            CALL para_env%sum(deb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dHppl    ", deb
            IF (use_virial) THEN
               e_dummy = third_tr(virial%pv_virial) - e_dummy
               CALL para_env%sum(e_dummy)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: Hppl       ", e_dummy
            END IF
         END IF
      END IF
      IF (ASSOCIATED(sap_ppnl)) THEN
         IF (debug_forces) THEN
            deb(1:3) = force(1)%gth_ppnl(1:3, 1)
            IF (use_virial) e_dummy = third_tr(virial%pv_virial)
         END IF
         CALL build_core_ppnl(scrm_kp, rho_ao_kp, force, &
                              virial, .TRUE., use_virial, 1, &
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                              nimages=1, basis_type="ORB")
         IF (debug_forces) THEN
            deb(1:3) = force(1)%gth_ppnl(1:3, 1) - deb(1:3)
            CALL para_env%sum(deb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dHppnl   ", deb
            IF (use_virial) THEN
               e_dummy = third_tr(virial%pv_virial) - e_dummy
               CALL para_env%sum(e_dummy)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: Hppnl      ", e_dummy
            END IF
         END IF
      END IF

      ! Get the different components of the KS potential
      NULLIFY (vxc_rspace, vtau_rspace, vadmm_rspace)
      IF (use_virial) THEN
         h_stress = 0.0_dp
         CALL ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, ehartree, exc, h_stress)
         ! Update virial
         virial%pv_ehartree = virial%pv_ehartree + h_stress/REAL(para_env%num_pe, dp)
         virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe, dp)
         IF (.NOT. do_exx) THEN
            virial%pv_exc = virial%pv_exc - virial%pv_xc
            virial%pv_virial = virial%pv_virial - virial%pv_xc
         ELSE
            virial%pv_xc = 0.0_dp
         END IF
         IF (debug_forces) THEN
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: xc       ", third_tr(h_stress)
            CALL para_env%sum(virial%pv_xc(1, 1))
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: Corr xc   ", third_tr(virial%pv_xc)
         END IF
      ELSE
         CALL ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, ehartree, exc)
      END IF

      ! Vhxc
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL auxbas_pw_pool%create_pw(vhxc_rspace, &
                                    use_data=REALDATA3D, in_space=REALSPACE)
      IF (use_virial) h_stress = virial%pv_virial
      IF (debug_forces) THEN
         deb(1:3) = force(1)%rho_elec(1:3, 1)
         IF (use_virial) e_dummy = third_tr(h_stress)
      END IF
      IF (do_exx) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(vh_rspace, vhxc_rspace)
            CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, -1.0_dp)
            CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                    hmat=matrix_ks(ispin), pmat=rho_ao(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE.)
            CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, 1.0_dp)
            CALL pw_axpy(vxc_rspace(ispin), vhxc_rspace)
            CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                    hmat=matrix_ks(ispin), pmat=matrix_p_mp2(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE.)
            IF (ASSOCIATED(vtau_rspace)) THEN
               CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                       hmat=matrix_ks(ispin), pmat=matrix_p_mp2(ispin), &
                                       qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE.)
            END IF
         END DO
      ELSE
         DO ispin = 1, nspins
            CALL pw_transfer(vh_rspace, vhxc_rspace)
            CALL pw_axpy(vxc_rspace(ispin), vhxc_rspace)
            CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                    hmat=matrix_ks(ispin), pmat=rho_ao(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE.)
            IF (ASSOCIATED(vtau_rspace)) THEN
               CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                       hmat=matrix_ks(ispin), pmat=rho_ao(ispin), &
                                       qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE.)
            END IF
         END DO
      END IF
      IF (debug_forces) THEN
         deb(1:3) = force(1)%rho_elec(1:3, 1) - deb(1:3)
         CALL para_env%sum(deb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dVhxc    ", deb
         IF (use_virial) THEN
            e_dummy = third_tr(virial%pv_virial) - e_dummy
            CALL para_env%sum(e_dummy)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: Vhxc      ", e_dummy
         END IF
      END IF
      IF (use_virial) THEN
         h_stress = virial%pv_virial - h_stress
         virial%pv_ehartree = virial%pv_ehartree + h_stress

         CALL qs_rho_get(p_env%rho1, rho_r=rho_mp2_r, tau_r=tau_mp2_r)
         e_xc = 0.0_dp
         DO ispin = 1, nspins
            ! The potentials have been scaled in ks_ref_potential, but for pw_integral_ab, we need the unscaled potentials
            CALL pw_scale(vxc_rspace(ispin), 1.0_dp/vxc_rspace(ispin)%pw_grid%dvol)
            e_xc = e_xc + pw_integral_ab(rho_mp2_r(ispin), vxc_rspace(ispin))
            IF (ASSOCIATED(vtau_rspace)) CALL pw_scale(vtau_rspace(ispin), 1.0_dp/vtau_rspace(ispin)%pw_grid%dvol)
            IF (ASSOCIATED(vtau_rspace)) e_xc = e_xc + pw_integral_ab(tau_mp2_r(ispin), vtau_rspace(ispin))
         END DO
         IF (debug_forces .AND. iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG VIRIAL:: vxc*d1   ", e_xc
         DO alpha = 1, 3
            virial%pv_exc(alpha, alpha) = virial%pv_exc(alpha, alpha) - e_xc/REAL(para_env%num_pe, dp)
            virial%pv_virial(alpha, alpha) = virial%pv_virial(alpha, alpha) - e_xc/REAL(para_env%num_pe, dp)
         END DO
      END IF
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%give_back_pw(vxc_rspace(ispin))
         IF (ASSOCIATED(vtau_rspace)) THEN
            CALL auxbas_pw_pool%give_back_pw(vtau_rspace(ispin))
         END IF
      END DO
      DEALLOCATE (vxc_rspace)
      CALL auxbas_pw_pool%give_back_pw(vhxc_rspace)
      IF (ASSOCIATED(vtau_rspace)) DEALLOCATE (vtau_rspace)

      DO ispin = 1, nspins
         CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, -1.0_dp)
      END DO

      ! pw stuff
      NULLIFY (poisson_env, auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ! get some of the grids ready
      CALL auxbas_pw_pool%create_pw(pot_r, &
                                    use_data=REALDATA3D, &
                                    in_space=REALSPACE)
      CALL auxbas_pw_pool%create_pw(pot_g, &
                                    use_data=COMPLEXDATA1D, &
                                    in_space=RECIPROCALSPACE)
      CALL auxbas_pw_pool%create_pw(rho_tot_g, &
                                    use_data=COMPLEXDATA1D, &
                                    in_space=RECIPROCALSPACE)

      CALL pw_zero(rho_tot_g)

      CALL qs_rho_get(p_env%rho1, rho_r=rho_mp2_r, rho_g=rho_mp2_g)
      DO ispin = 1, nspins
         CALL pw_axpy(rho_mp2_g(ispin), rho_tot_g)
      END DO

      IF (use_virial) THEN
         ALLOCATE (dvg(3))
         DO idir = 1, 3
            CALL auxbas_pw_pool%create_pw(dvg(idir), &
                                          use_data=COMPLEXDATA1D, &
                                          in_space=RECIPROCALSPACE)
         END DO
      END IF

      ! calculate the MP2 potential
      CALL pw_poisson_solve(poisson_env, rho_tot_g, vhartree=pot_g, dvhartree=dvg)
      CALL pw_transfer(pot_g, pot_r)
      CALL pw_scale(pot_r, pot_r%pw_grid%dvol)
      CALL pw_axpy(pot_r, vh_rspace)

      ! calculate core forces
      CALL integrate_v_core_rspace(vh_rspace, qs_env)
      IF (debug_forces) THEN
         deb(:) = force(1)%rho_core(:, 1)
         CALL para_env%sum(deb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: core      ", deb
         IF (use_virial) THEN
            e_dummy = third_tr(virial%pv_virial) - e_dummy
            CALL para_env%sum(e_dummy)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: core      ", e_dummy
         END IF
      END IF
      CALL auxbas_pw_pool%give_back_pw(vh_rspace)

      IF (use_virial) THEN
         ! update virial if necessary with the volume term
         ! first create pw auxiliary stuff
         CALL auxbas_pw_pool%create_pw(temp_pw_g, &
                                       use_data=COMPLEXDATA1D, &
                                       in_space=RECIPROCALSPACE)

         ! make a copy of the MP2 density in G space
         CALL pw_copy(rho_tot_g, temp_pw_g)

         ! calculate total SCF density and potential
         CALL pw_copy(rho_g(1), rho_tot_g)
         IF (nspins == 2) CALL pw_axpy(rho_g(2), rho_tot_g)
         CALL pw_axpy(rho_core, rho_tot_g)
         CALL pw_poisson_solve(poisson_env, rho_tot_g, vhartree=pot_g)

         ! finally update virial with the volume contribution
         e_hartree = pw_integral_ab(temp_pw_g, pot_g)
         IF (debug_forces .AND. iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: vh1*d0   ", e_hartree

         h_stress = 0.0_dp
         DO alpha = 1, 3
            comp = 0
            comp(alpha) = 1
            CALL pw_copy(pot_g, rho_tot_g)
            CALL pw_derive(rho_tot_g, comp)
            h_stress(alpha, alpha) = -e_hartree
            DO beta = alpha, 3
               h_stress(alpha, beta) = h_stress(alpha, beta) &
                                       - 2.0_dp*pw_integral_ab(rho_tot_g, dvg(beta))/fourpi
               h_stress(beta, alpha) = h_stress(alpha, beta)
            END DO
         END DO
         IF (debug_forces .AND. iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: Hartree  ", third_tr(h_stress)

         ! free stuff
         CALL auxbas_pw_pool%give_back_pw(temp_pw_g)
         DO idir = 1, 3
            CALL auxbas_pw_pool%give_back_pw(dvg(idir))
         END DO
         DEALLOCATE (dvg)

         virial%pv_ehartree = virial%pv_ehartree + h_stress/REAL(para_env%num_pe, dp)
         virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe, dp)

      END IF

      DO ispin = 1, nspins
         CALL dbcsr_set(p_env%kpp1(ispin)%matrix, 0.0_dp)
         IF (dft_control%do_admm) CALL dbcsr_set(p_env%kpp1_admm(ispin)%matrix, 0.0_dp)
      END DO

      CALL get_qs_env(qs_env=qs_env, linres_control=linres_control)

      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF

      IF (use_virial) THEN
         h_stress = 0.0_dp
         pv_virial = virial%pv_virial
      END IF
      IF (debug_forces) THEN
         deb = force(1)%rho_elec(1:3, 1)
         IF (use_virial) e_dummy = third_tr(pv_virial)
      END IF
      CALL apply_2nd_order_kernel(qs_env, p_env, .FALSE., .TRUE., use_virial, h_stress)
      IF (use_virial) THEN
         virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_virial)
         IF (debug_forces) THEN
            e_dummy = third_tr(virial%pv_virial - pv_virial)
            CALL para_env%sum(e_dummy)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: Kh       ", e_dummy
         END IF
         virial%pv_exc = virial%pv_exc + h_stress
         virial%pv_virial = virial%pv_virial + h_stress
         IF (debug_forces) THEN
            e_dummy = third_tr(h_stress)
            CALL para_env%sum(e_dummy)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: Kxc       ", e_dummy
         END IF
      END IF
      IF (debug_forces) THEN
         deb(:) = force(1)%rho_elec(:, 1) - deb
         CALL para_env%sum(deb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P0*Khxc    ", deb
         IF (use_virial) THEN
            e_dummy = third_tr(virial%pv_virial) - e_dummy
            CALL para_env%sum(e_dummy)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: Khxc      ", e_dummy
         END IF
      END IF

      ! hfx section
      NULLIFY (hfx_sections)
      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx)
      IF (do_hfx) THEN
         IF (do_exx) THEN
            IF (dft_control%do_admm) THEN
               CALL get_admm_env(qs_env%admm_env, rho_aux_fit=rho_aux)
               CALL qs_rho_get(rho_aux, rho_ao=rho_ao_aux, rho_ao_kp=rho_ao_kp)
               rho1 => p_env%p1_admm
            ELSE
               rho1 => p_env%p1
            END IF
         ELSE
            IF (dft_control%do_admm) THEN
               CALL get_admm_env(qs_env%admm_env, rho_aux_fit=rho_aux)
               CALL qs_rho_get(rho_aux, rho_ao=rho_ao_aux, rho_ao_kp=rho_ao_kp)
               DO ispin = 1, nspins
                  CALL dbcsr_add(rho_ao_aux(ispin)%matrix, p_env%p1_admm(ispin)%matrix, 1.0_dp, 1.0_dp)
               END DO
               rho1 => p_env%p1_admm
            ELSE
               DO ispin = 1, nspins
                  CALL dbcsr_add(rho_ao(ispin)%matrix, p_env%p1(ispin)%matrix, 1.0_dp, 1.0_dp)
               END DO
               rho1 => p_env%p1
            END IF
         END IF

         IF (x_data(1, 1)%do_hfx_ri) THEN

            CALL hfx_ri_update_forces(qs_env, x_data(1, 1)%ri_data, nspins, &
                                      x_data(1, 1)%general_parameter%fraction, &
                                      rho_ao=rho_ao_kp, rho_ao_resp=rho1, &
                                      use_virial=use_virial, resp_only=do_exx)

         ELSE
            CALL derivatives_four_center(qs_env, rho_ao_kp, rho1, hfx_sections, para_env, &
                                         1, use_virial, resp_only=do_exx)
         END IF

         IF (use_virial) THEN
            virial%pv_exx = virial%pv_exx - virial%pv_fock_4c
            virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
         END IF
         IF (debug_forces) THEN
            deb(1:3) = force(1)%fock_4c(1:3, 1) - deb(1:3)
            CALL para_env%sum(deb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*hfx  ", deb
            IF (use_virial) THEN
               e_dummy = third_tr(virial%pv_fock_4c)
               CALL para_env%sum(e_dummy)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: hfx    ", e_dummy
            END IF
         END IF

         IF (.NOT. do_exx) THEN
         IF (dft_control%do_admm) THEN
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            DO ispin = 1, nspins
               CALL dbcsr_add(rho_ao_aux(ispin)%matrix, p_env%p1_admm(ispin)%matrix, 1.0_dp, -1.0_dp)
            END DO
         ELSE
            DO ispin = 1, nspins
               CALL dbcsr_add(rho_ao(ispin)%matrix, p_env%p1(ispin)%matrix, 1.0_dp, -1.0_dp)
            END DO
         END IF
         END IF

         IF (dft_control%do_admm) THEN
            IF (debug_forces) THEN
               deb = force(1)%overlap_admm(1:3, 1)
               IF (use_virial) e_dummy = third_tr(virial%pv_virial)
            END IF
            ! The 2nd order kernel contains a factor of two in apply_xc_admm_ao which we don't need for the projection derivatives
            IF (nspins == 1) CALL dbcsr_scale(p_env%kpp1_admm(1)%matrix, 0.5_dp)
            CALL admm_projection_derivative(qs_env, p_env%kpp1_admm, rho_ao)
            IF (debug_forces) THEN
               deb(:) = force(1)%overlap_admm(:, 1) - deb
               CALL para_env%sum(deb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*KADMM*dS'", deb
               IF (use_virial) THEN
                  e_dummy = third_tr(virial%pv_virial) - e_dummy
                  CALL para_env%sum(e_dummy)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: KADMM*S'  ", e_dummy
               END IF
            END IF

            ALLOCATE (matrix_ks_aux(nspins))
            DO ispin = 1, nspins
               NULLIFY (matrix_ks_aux(ispin)%matrix)
               ALLOCATE (matrix_ks_aux(ispin)%matrix)
               CALL dbcsr_copy(matrix_ks_aux(ispin)%matrix, p_env%kpp1_admm(ispin)%matrix)
               CALL dbcsr_set(matrix_ks_aux(ispin)%matrix, 0.0_dp)
            END DO

            ! Calculate kernel
            CALL tddft_hfx_matrix(matrix_ks_aux, rho_ao_aux, qs_env, .FALSE.)

            IF (qs_env%admm_env%aux_exch_func /= do_admm_aux_exch_func_none) THEN
               CALL get_qs_env(qs_env, ks_env=ks_env)
               CALL get_admm_env(qs_env%admm_env, task_list_aux_fit=task_list_aux_fit)

               DO ispin = 1, nspins
                  CALL dbcsr_add(rho_ao_aux(ispin)%matrix, p_env%p1_admm(ispin)%matrix, 1.0_dp, 1.0_dp)
               END DO

               IF (use_virial) THEN
                  CALL qs_rho_get(p_env%rho1_admm, rho_r=rho_mp2_r_aux)
                  e_xc = 0.0_dp
                  DO ispin = 1, nspins
                     e_xc = e_xc + pw_integral_ab(rho_mp2_r_aux(ispin), vadmm_rspace(ispin))
                  END DO

                  e_xc = -e_xc/vadmm_rspace(1)%pw_grid%dvol/REAL(para_env%num_pe, dp)

                  ! Update the virial
                  DO alpha = 1, 3
                     virial%pv_exc(alpha, alpha) = virial%pv_exc(alpha, alpha) + e_xc
                     virial%pv_virial(alpha, alpha) = virial%pv_virial(alpha, alpha) + e_xc
                  END DO
                  IF (debug_forces) THEN
                     IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: P1*VADMM  ", e_xc
                  END IF
               END IF

               IF (use_virial) h_stress = virial%pv_virial
               IF (debug_forces) THEN
                  deb = force(1)%rho_elec(1:3, 1)
                  IF (use_virial) e_dummy = third_tr(virial%pv_virial)
               END IF
               DO ispin = 1, nspins
                  IF (do_exx) THEN
                     CALL integrate_v_rspace(v_rspace=vadmm_rspace(ispin), hmat=matrix_ks_aux(ispin), qs_env=qs_env, &
                                             calculate_forces=.TRUE., basis_type="AUX_FIT", &
                                             task_list_external=task_list_aux_fit, pmat=matrix_p_mp2_admm(ispin))
                  ELSE
                     CALL integrate_v_rspace(v_rspace=vadmm_rspace(ispin), hmat=matrix_ks_aux(ispin), qs_env=qs_env, &
                                             calculate_forces=.TRUE., basis_type="AUX_FIT", &
                                             task_list_external=task_list_aux_fit, pmat=rho_ao_aux(ispin))
                  END IF
                  CALL auxbas_pw_pool%give_back_pw(vadmm_rspace(ispin))
               END DO
               IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - h_stress)
               DEALLOCATE (vadmm_rspace)
               IF (debug_forces) THEN
                  deb(:) = force(1)%rho_elec(:, 1) - deb
                  CALL para_env%sum(deb)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*VADMM' ", deb
                  IF (use_virial) THEN
                     e_dummy = third_tr(virial%pv_virial) - e_dummy
                     CALL para_env%sum(e_dummy)
                     IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: VADMM'   ", e_dummy
                  END IF
               END IF

               DO ispin = 1, nspins
                  CALL dbcsr_add(rho_ao_aux(ispin)%matrix, p_env%p1_admm(ispin)%matrix, 1.0_dp, -1.0_dp)
               END DO

            END IF

            DO ispin = 1, nspins
               CALL dbcsr_add(rho_ao(ispin)%matrix, p_env%p1(ispin)%matrix, 1.0_dp, 1.0_dp)
            END DO

            IF (debug_forces) THEN
               deb = force(1)%overlap_admm(1:3, 1)
               IF (use_virial) e_dummy = third_tr(virial%pv_virial)
            END IF
            ! Add the second half of the projector deriatives contracting the first order density matrix with the fockian in the auxiliary basis
            IF (do_exx) THEN
               CALL admm_projection_derivative(qs_env, matrix_ks_aux, matrix_p_mp2)
            ELSE
               CALL admm_projection_derivative(qs_env, matrix_ks_aux, rho_ao)
            END IF
            IF (debug_forces) THEN
               deb(:) = force(1)%overlap_admm(:, 1) - deb
               CALL para_env%sum(deb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*VADMM*dS'", deb
               IF (use_virial) THEN
                  e_dummy = third_tr(virial%pv_virial) - e_dummy
                  CALL para_env%sum(e_dummy)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: VADMM*S'  ", e_dummy
               END IF
            END IF

            DO ispin = 1, nspins
               CALL dbcsr_add(rho_ao(ispin)%matrix, p_env%p1(ispin)%matrix, 1.0_dp, -1.0_dp)
            END DO

            DO ispin = 1, nspins
               CALL dbcsr_release(matrix_ks_aux(ispin)%matrix)
               DEALLOCATE (matrix_ks_aux(ispin)%matrix)
            END DO
            DEALLOCATE (matrix_ks_aux)
         END IF
      END IF

      CALL dbcsr_scale(p_env%w1(1)%matrix, -1.0_dp)

      ! Finish matrix_w_mp2 with occ-occ block
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), homo=nocc, nmo=alpha)
         CALL calculate_whz_matrix(mos(ispin)%mo_coeff, p_env%kpp1(ispin)%matrix, &
                                   p_env%w1(1)%matrix, 1.0_dp, nocc)
      END DO

      IF (debug_forces .AND. use_virial) e_dummy = third_tr(virial%pv_virial)

      NULLIFY (scrm)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=p_env%w1(1)%matrix)
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (debug_forces) THEN
         deb = force(1)%overlap(1:3, 1)
         CALL para_env%sum(deb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: W*dS     ", deb
         IF (use_virial) THEN
            e_dummy = third_tr(virial%pv_virial) - e_dummy
            CALL para_env%sum(e_dummy)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: S         ", e_dummy
         END IF
      END IF

      CALL auxbas_pw_pool%give_back_pw(pot_r)
      CALL auxbas_pw_pool%give_back_pw(pot_g)
      CALL auxbas_pw_pool%give_back_pw(rho_tot_g)

      ! Release linres stuff
      CALL p_env_release(p_env)
      DEALLOCATE (p_env)
      NULLIFY (qs_env%mp2_env%ri_grad%p_env)

      CALL sum_qs_force(force, qs_env%mp2_env%ri_grad%mp2_force)
      CALL deallocate_qs_force(qs_env%mp2_env%ri_grad%mp2_force)

      IF (use_virial) THEN
         virial%pv_mp2 = qs_env%mp2_env%ri_grad%mp2_virial
         virial%pv_virial = virial%pv_virial + qs_env%mp2_env%ri_grad%mp2_virial
         IF (debug_forces) THEN
            e_dummy = third_tr(virial%pv_mp2)
            CALL para_env%sum(e_dummy)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,F16.8)") "DEBUG VIRIAL:: MP2nonsep  ", e_dummy
         END IF
      END IF
      ! Rewind the change from the beginning
      IF (use_virial) virial%pv_calculate = .FALSE.

      ! Add the dispersion forces and virials
      CALL get_qs_env(qs_env, dispersion_env=dispersion_env)
      CALL calculate_dispersion_pairpot(qs_env, dispersion_env, e_dummy, .TRUE.)

      CALL cp_print_key_finished_output(iounit, logger, input, &
                                        "DFT%XC%WF_CORRELATION%PRINT")

      CALL timestop(handle)

   END SUBROUTINE update_mp2_forces

! **************************************************************************************************
!> \brief Calculates the third of the trace of a 3x3 matrix, for debugging purposes
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION third_tr(matrix)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: matrix
      REAL(KIND=dp)                                      :: third_tr

      third_tr = (matrix(1, 1) + matrix(2, 2) + matrix(3, 3))/3.0_dp

   END FUNCTION third_tr

END MODULE mp2_cphf
