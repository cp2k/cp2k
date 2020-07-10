!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
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
   USE admm_types,                      ONLY: admm_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_p_type,&
                                              dbcsr_scale,&
                                              dbcsr_set
   USE energy_corrections,              ONLY: calc_ks_response_potentials
   USE hfx_types,                       ONLY: alloc_containers,&
                                              hfx_container_type,&
                                              hfx_init_container,&
                                              hfx_type
   USE input_constants,                 ONLY: ot_precond_full_all
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE kinds,                           ONLY: dp
   USE linear_systems,                  ONLY: solve_system
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: fourpi
   USE message_passing,                 ONLY: mp_sum
   USE mp2_types,                       ONLY: mp2_type
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
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_2nd_kernel_ao,                ONLY: apply_hfx_ao,&
                                              apply_xc_admm_ao
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kpp1_env_methods,             ONLY: calc_kpp1
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_linres_types,                 ONLY: linres_control_create,&
                                              linres_control_release,&
                                              linres_control_type
   USE qs_p_env_methods,                ONLY: p_env_check_i_alloc,&
                                              p_env_create,&
                                              p_env_psi0_changed,&
                                              p_env_update_rho
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: xc_vxc_pw_create

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_cphf'

   PUBLIC :: solve_z_vector_eq

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
      TYPE(mp2_type), INTENT(IN), POINTER                :: mp2_env
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(dft_control_type), INTENT(IN), POINTER        :: dft_control
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: mo_coeff
      INTEGER, INTENT(IN)                                :: nmo
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'solve_z_vector_eq'

      INTEGER :: bin, dimen, handle, handle2, i, i_global, i_thread, iiB, irep, ispin, j_global, &
         jjB, my_bin_size, n_rep_hf, n_threads, ncol_local, nrow_local, nspins, transf_type_in, &
         transf_type_out
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: virtual
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: alpha_beta, do_dynamic_load_balancing, &
                                                            do_hfx, hfx_treat_lsd_in_core, &
                                                            restore_p_screen
      REAL(KIND=dp)                                      :: out_alpha
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: L_jb, mo_coeff_o, mo_coeff_v, P_ia, &
                                                            P_mo, W_mo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: fm_back, fm_G_mu_nu
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_p_mp2, matrix_s, &
                                                            matrix_w_mp2, P_mu_nu, rho_ao
      TYPE(hfx_container_type), DIMENSION(:), POINTER    :: integral_containers
      TYPE(hfx_container_type), POINTER                  :: maxval_container
      TYPE(hfx_type), POINTER                            :: actual_x_data
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL timeset(routineN, handle)

      ! start collecting stuff
      dimen = nmo
      NULLIFY (input, matrix_s, blacs_env, rho, matrix_w_mp2, &
               matrix_p_mp2, matrix_ks)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      input=input, &
                      matrix_s=matrix_s, &
                      matrix_ks=matrix_ks, &
                      matrix_p_mp2=matrix_p_mp2, &
                      matrix_w_mp2=matrix_w_mp2, &
                      blacs_env=blacs_env, &
                      rho=rho)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ! Get number of relevant spin states
      nspins = dft_control%nspins
      alpha_beta = (nspins == 2)

      ! mp2 matrices
      ALLOCATE (P_mo(nspins), W_mo(nspins), L_jb(nspins))

      ! Pack helper matrices
      DO ispin = 1, nspins
         P_mo(ispin)%matrix => mp2_env%ri_grad%P_mo(ispin)%matrix
         W_mo(ispin)%matrix => mp2_env%ri_grad%W_mo(ispin)%matrix
         L_jb(ispin)%matrix => mp2_env%ri_grad%L_jb(ispin)%matrix
      ENDDO
      ALLOCATE (virtual(nspins))
      virtual(:) = dimen - homo(:)

      ! hfx section
      NULLIFY (hfx_sections)
      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx, n_repetition=n_rep_hf)
      IF (do_hfx) THEN
         CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                   i_rep_section=1)
      END IF

      NULLIFY (P_mu_nu)
      CALL dbcsr_allocate_matrix_set(P_mu_nu, nspins)
      DO ispin = 1, nspins
         ALLOCATE (P_mu_nu(ispin)%matrix)
         CALL dbcsr_copy(P_mu_nu(ispin)%matrix, rho_ao(1)%matrix, name="P_mu_nu")
         CALL dbcsr_set(P_mu_nu(ispin)%matrix, 0.0_dp)
      END DO

      NULLIFY (fm_G_mu_nu, fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen, ncol_global=dimen)
      CALL cp_fm_create(fm_G_mu_nu, fm_struct_tmp, name="G_mu_nu")
      CALL cp_fm_create(fm_back, fm_struct_tmp, name="fm_back")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_set_all(fm_G_mu_nu, 0.0_dp)
      CALL cp_fm_set_all(fm_back, 0.0_dp)

      ALLOCATE (mo_coeff_o(nspins), mo_coeff_v(nspins))
      DO ispin = 1, nspins
         NULLIFY (mo_coeff_o(ispin)%matrix, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=dimen, ncol_global=homo(ispin))
         CALL cp_fm_create(mo_coeff_o(ispin)%matrix, fm_struct_tmp, name="mo_coeff_o")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(mo_coeff_o(ispin)%matrix, 0.0_dp)
         CALL cp_fm_to_fm_submat(msource=mo_coeff(ispin)%matrix, mtarget=mo_coeff_o(ispin)%matrix, &
                                 nrow=dimen, ncol=homo(ispin), &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=1)

         NULLIFY (mo_coeff_v(ispin)%matrix, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=dimen, ncol_global=virtual(ispin))
         CALL cp_fm_create(mo_coeff_v(ispin)%matrix, fm_struct_tmp, name="mo_coeff_v")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(mo_coeff_v(ispin)%matrix, 0.0_dp)
         CALL cp_fm_to_fm_submat(msource=mo_coeff(ispin)%matrix, mtarget=mo_coeff_v(ispin)%matrix, &
                                 nrow=dimen, ncol=virtual(ispin), &
                                 s_firstrow=1, s_firstcol=homo(ispin) + 1, &
                                 t_firstrow=1, t_firstcol=1)
      END DO

      IF (do_hfx) THEN
         ! here we check if we have to reallocate the HFX container
         IF (mp2_env%ri_mp2%free_hfx_buffer) THEN
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
            IF (mp2_env%ri_mp2%free_hfx_buffer) THEN
               mp2_env%p_screen = .FALSE.
            ELSE
               mp2_env%p_screen = .TRUE.
            ENDIF
         END IF
      END IF

      ! Prepare arrays for linres code
      NULLIFY (linres_control)
      CALL linres_control_create(linres_control)
      linres_control%do_kernel = .TRUE.
      linres_control%lr_triplet = .FALSE.
      linres_control%linres_restart = .FALSE.
      linres_control%max_iter = mp2_env%ri_grad%cphf_max_num_iter
      linres_control%eps = mp2_env%ri_grad%cphf_eps_conv
      linres_control%eps_filter = mp2_env%mp2_gpw%eps_filter
      linres_control%restart_every = 50
      linres_control%preconditioner_type = ot_precond_full_all
      linres_control%energy_gap = 0.02_dp

      CALL p_env_create(p_env, qs_env, p1_option=P_mu_nu, orthogonal_orbitals=.TRUE., linres_control=linres_control)
      CALL set_qs_env(qs_env, linres_control=linres_control)
      CALL p_env_psi0_changed(p_env, qs_env)
      p_env%new_preconditioner = .TRUE.
      CALL p_env_check_i_alloc(p_env, qs_env)

      ! update Lagrangian with the CPHF like update, occ-occ block, first call (recompute hfx integrals if needed)
      transf_type_in = 1
      transf_type_out = 1
      out_alpha = 0.5_dp
      ! In alpha-beta case, L_bj_alpha has Coulomb and XC alpha-alpha part
      ! and (only) Coulomb alpha-beta part and vice versa.

      ! Complete in closed shell case, alpha-alpha (Coulomb and XC)
      ! part of L_bj(alpha) for open shell

      CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                            mo_coeff, mo_coeff_o, &
                            mo_coeff_v, Eigenval, p_env, &
                            do_hfx, &
                            P_mo, fm_G_mu_nu, fm_back, transf_type_in, out_alpha, &
                            L_jb, transf_type_out, &
                            recalc_hfx_integrals=mp2_env%ri_mp2%free_hfx_buffer)

      ! update Lagrangian with the CPHF like update, virt-virt block
      transf_type_in = 2
      transf_type_out = 1
      out_alpha = 0.5_dp

      CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                            mo_coeff, mo_coeff_o, &
                            mo_coeff_v, Eigenval, p_env, &
                            do_hfx, &
                            P_mo, fm_G_mu_nu, fm_back, transf_type_in, out_alpha, &
                            L_jb, transf_type_out)

      ! at this point Lagrnagian is completed ready to solve the Z-vector equations
      ! P_ia will contain the solution of these equations
      ALLOCATE (P_ia(nspins))
      DO ispin = 1, nspins
         NULLIFY (P_ia(ispin)%matrix, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=homo(ispin), ncol_global=virtual(ispin))
         CALL cp_fm_create(P_ia(ispin)%matrix, fm_struct_tmp, name="P_ia")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_set_all(P_ia(ispin)%matrix, 0.0_dp)
      END DO

      CALL solve_z_vector_eq_low(qs_env, mp2_env, para_env, homo, virtual, dimen, unit_nr, nspins, &
                                 mo_coeff, mo_coeff_o, mo_coeff_v, Eigenval, p_env, blacs_env, &
                                 do_hfx, L_jb, fm_G_mu_nu, fm_back, P_ia)

      ! release Lagrangian
      DO ispin = 1, nspins
         CALL cp_fm_release(L_jb(ispin)%matrix)
      END DO
      DEALLOCATE (L_jb)

      DO ispin = 1, nspins
         ! update the MP2-MO density matrix with the occ-virt block
         CALL cp_fm_to_fm_submat(msource=P_ia(ispin)%matrix, mtarget=P_mo(ispin)%matrix, &
                                 nrow=homo(ispin), ncol=virtual(ispin), &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=homo(ispin) + 1)
         CALL cp_fm_release(P_ia(ispin)%matrix)
         ! transpose P_MO matrix (easy way to symmetrize)
         CALL cp_fm_set_all(fm_back, 0.0_dp)
         ! P_mo now is ready
         CALL cp_fm_upper_to_full(matrix=P_mo(ispin)%matrix, work=fm_back)
      END DO
      DEALLOCATE (P_ia)

      ! do the final update to MP2 energy weighted matrix W_MO
      DO ispin = 1, nspins
         CALL cp_fm_get_info(matrix=W_mo(ispin)%matrix, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            IF (j_global <= homo(ispin)) THEN
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  W_mo(ispin)%matrix%local_data(iiB, jjB) = W_mo(ispin)%matrix%local_data(iiB, jjB) &
                                                            - P_mo(ispin)%matrix%local_data(iiB, jjB)*Eigenval(j_global, ispin)
               END DO
            ELSE
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  IF (i_global <= homo(ispin)) THEN
                     ! virt-occ
                     W_mo(ispin)%matrix%local_data(iiB, jjB) = W_mo(ispin)%matrix%local_data(iiB, jjB) &
                                                               - P_mo(ispin)%matrix%local_data(iiB, jjB)*Eigenval(i_global, ispin)
                  ELSE
                     ! virt-virt
                     W_mo(ispin)%matrix%local_data(iiB, jjB) = W_mo(ispin)%matrix%local_data(iiB, jjB) &
                                                               - P_mo(ispin)%matrix%local_data(iiB, jjB)*Eigenval(j_global, ispin)
                  END IF
               END DO
            END IF
         END DO
      END DO

      ! complete the occ-occ block of W_mo with a CPHF like update
      transf_type_in = 4
      transf_type_out = 2
      out_alpha = -0.5_dp

      CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                            mo_coeff, mo_coeff_o, &
                            mo_coeff_v, Eigenval, p_env, &
                            do_hfx, &
                            P_mo, fm_G_mu_nu, fm_back, transf_type_in, out_alpha, &
                            W_mo, transf_type_out)

      ! release fm stuff
      CALL cp_fm_release(fm_G_mu_nu)
      DO ispin = 1, nspins
         CALL cp_fm_release(mo_coeff_o(ispin)%matrix)
         CALL cp_fm_release(mo_coeff_v(ispin)%matrix)
      END DO
      DEALLOCATE (mo_coeff_o, mo_coeff_v)

      ! backtransform into AO basis, since P_mo and W_mo
      ! are symmetric (in principle), no need to symmetrize
      ! first W_mo
      DO ispin = 1, nspins
         CALL cp_gemm('N', 'N', dimen, dimen, dimen, 1.0_dp, &
                      mo_coeff(ispin)%matrix, W_mo(ispin)%matrix, 0.0_dp, fm_back, &
                      a_first_col=1, &
                      a_first_row=1, &
                      b_first_col=1, &
                      b_first_row=1, &
                      c_first_col=1, &
                      c_first_row=1)
         CALL cp_gemm('N', 'T', dimen, dimen, dimen, 1.0_dp, &
                      fm_back, mo_coeff(ispin)%matrix, 0.0_dp, W_mo(ispin)%matrix, &
                      a_first_col=1, &
                      a_first_row=1, &
                      b_first_col=1, &
                      b_first_row=1, &
                      c_first_col=1, &
                      c_first_row=1)
         ! and P_mo
         CALL cp_gemm('N', 'N', dimen, dimen, dimen, 1.0_dp, &
                      mo_coeff(ispin)%matrix, P_mo(ispin)%matrix, 0.0_dp, fm_back, &
                      a_first_col=1, &
                      a_first_row=1, &
                      b_first_col=1, &
                      b_first_row=1, &
                      c_first_col=1, &
                      c_first_row=1)
         CALL cp_gemm('N', 'T', dimen, dimen, dimen, 1.0_dp, &
                      fm_back, mo_coeff(ispin)%matrix, 0.0_dp, P_mo(ispin)%matrix, &
                      a_first_col=1, &
                      a_first_row=1, &
                      b_first_col=1, &
                      b_first_row=1, &
                      c_first_col=1, &
                      c_first_row=1)
      END DO

      ! copy W_mo into dbcsr
      DO ispin = 1, nspins
         CALL copy_fm_to_dbcsr(W_mo(ispin)%matrix, matrix_w_mp2(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO

      ! create mp2 DBCSR density
      CALL dbcsr_allocate_matrix_set(matrix_p_mp2, nspins)
      DO ispin = 1, nspins
         ALLOCATE (matrix_p_mp2(ispin)%matrix)
         CALL dbcsr_copy(matrix_p_mp2(ispin)%matrix, rho_ao(ispin)%matrix, &
                         name="P MATRIX MP2")
         CALL dbcsr_set(matrix_p_mp2(ispin)%matrix, 0.0_dp)
         CALL copy_fm_to_dbcsr(P_mo(ispin)%matrix, matrix_p_mp2(ispin)%matrix, keep_sparsity=.TRUE.)
      END DO
      CALL set_ks_env(ks_env, matrix_p_mp2=matrix_p_mp2)

      CALL cp_fm_release(fm_back)

      ! release remaining fm stuff
      DO ispin = 1, nspins
         CALL cp_fm_release(W_mo(ispin)%matrix)
         CALL cp_fm_release(P_mo(ispin)%matrix)
      END DO
      DEALLOCATE (P_mo, W_mo)

      CALL update_mp2_forces(qs_env, para_env)

      ! We will need one more hfx calculation for HF gradient part
      mp2_env%not_last_hfx = .FALSE.
      mp2_env%p_screen = restore_p_screen

      ! Release linres stuff
      CALL p_env_release(p_env)
      CALL linres_control_release(linres_control)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Here we performe the CPHF like update using GPW,
!>        transf_type_in  defines the type of transformation for the matrix in input
!>        transf_type_in = 1 -> occ-occ back transformation
!>        transf_type_in = 2 -> virt-virt back transformation
!>        transf_type_in = 3 -> occ-virt back transformation including the
!>                              eigenvalues energy differences for the diagonal elements
!>        transf_type_in = 4 -> full range
!>        transf_type_out defines the type of transformation for the matrix in output
!>        transf_type_out = 1 -> occ-vit transformation
!>        transf_type_out = 2 -> occ-occ transformation
!> \param qs_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param nspins ...
!> \param mo_coeff ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param Eigenval ...
!> \param p_env ...
!> \param do_hfx ...
!> \param fm_mo ...
!> \param fm_ao ...
!> \param fm_back ...
!> \param transf_type_in ...
!> \param out_alpha ...
!> \param fm_mo_out ...
!> \param transf_type_out ...
!> \param recalc_hfx_integrals ...
!> \author Mauro Del Ben, Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                               mo_coeff, mo_coeff_o, mo_coeff_v, Eigenval, p_env, &
                               do_hfx, fm_mo, fm_ao, fm_back, transf_type_in, out_alpha, &
                               fm_mo_out, transf_type_out, recalc_hfx_integrals)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      INTEGER, INTENT(IN)                                :: nspins, dimen
      INTEGER, DIMENSION(nspins), INTENT(IN)             :: virtual, homo
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      TYPE(qs_p_env_type), INTENT(IN), POINTER           :: p_env
      LOGICAL, INTENT(IN)                                :: do_hfx
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mo
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_ao, fm_back
      INTEGER, INTENT(IN)                                :: transf_type_in
      REAL(KIND=dp), INTENT(IN)                          :: out_alpha
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: fm_mo_out
      INTEGER, INTENT(IN)                                :: transf_type_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: recalc_hfx_integrals

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cphf_like_update'

      INTEGER                                            :: handle, i_global, iiB, ispin, j_global, &
                                                            jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: lr_triplet, lrigpw, &
                                                            my_recalc_hfx_integrals
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: mat_in, mat_out
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:), &
         TARGET                                          :: mat_mu_nu, rho_work_ao, rho_work_ao_admm
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho1_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho1, rho1_xc
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      NULLIFY (ks_env)
      CALL timeset(routineN, handle)

      my_recalc_hfx_integrals = .FALSE.
      IF (PRESENT(recalc_hfx_integrals)) my_recalc_hfx_integrals = recalc_hfx_integrals

      CALL qs_rho_get(p_env%rho1, rho_ao=rho1_ao)

      ! Determine the first-order density matrices in AO basis
      DO ispin = 1, nspins

         mat_in => fm_mo(ispin)%matrix

         ! perform back transformation
         SELECT CASE (transf_type_in)
         CASE (1)
            ! occ-occ block
            CALL cp_gemm('N', 'N', dimen, homo(ispin), homo(ispin), 1.0_dp, &
                         mo_coeff_o(ispin)%matrix, mat_in, 0.0_dp, fm_back, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)
            CALL cp_gemm('N', 'T', dimen, dimen, homo(ispin), 1.0_dp, &
                         fm_back, mo_coeff_o(ispin)%matrix, 0.0_dp, fm_ao, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)

         CASE (2)
            ! virt-virt block
            CALL cp_gemm('N', 'N', dimen, virtual(ispin), virtual(ispin), 1.0_dp, &
                         mo_coeff_v(ispin)%matrix, mat_in, 0.0_dp, fm_back, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=homo(ispin) + 1, &
                         b_first_row=homo(ispin) + 1, &
                         c_first_col=1, &
                         c_first_row=1)
            CALL cp_gemm('N', 'T', dimen, dimen, virtual(ispin), 1.0_dp, &
                         fm_back, mo_coeff_v(ispin)%matrix, 0.0_dp, fm_ao, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)

         CASE (3)
            ! virt-occ blocks
            CALL cp_gemm('N', 'N', dimen, virtual(ispin), homo(ispin), 1.0_dp, &
                         mo_coeff_o(ispin)%matrix, mat_in, 0.0_dp, fm_back, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)
            CALL cp_gemm('N', 'T', dimen, dimen, virtual(ispin), 1.0_dp, &
                         fm_back, mo_coeff_v(ispin)%matrix, 0.0_dp, fm_ao, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)
            ! and symmetrize (here again multiply instead of transposing)
            CALL cp_gemm('N', 'T', dimen, homo(ispin), virtual(ispin), 1.0_dp, &
                         mo_coeff_v(ispin)%matrix, mat_in, 0.0_dp, fm_back, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)
            CALL cp_gemm('N', 'T', dimen, dimen, homo(ispin), 0.5_dp, &
                         fm_back, mo_coeff_o(ispin)%matrix, 0.5_dp, fm_ao, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)

         CASE (4)
            ! all-all block
            CALL cp_gemm('N', 'N', dimen, dimen, dimen, 1.0_dp, &
                         mo_coeff(ispin)%matrix, mat_in, 0.0_dp, fm_back, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)
            CALL cp_gemm('N', 'T', dimen, dimen, dimen, 1.0_dp, &
                         fm_back, mo_coeff(ispin)%matrix, 0.0_dp, fm_ao, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)

         CASE DEFAULT
            ! nothing
         END SELECT

         ! copy fm into DBCSR
         CALL dbcsr_set(p_env%p1(ispin)%matrix, 0.0_dp)
         CALL copy_fm_to_dbcsr(fm_ao, p_env%p1(ispin)%matrix, keep_sparsity=.TRUE.)

         CALL dbcsr_copy(rho1_ao(ispin)%matrix, p_env%p1(ispin)%matrix)
      END DO

      CALL get_qs_env(qs_env, ks_env=ks_env, x_data=x_data, dft_control=dft_control)

      CALL p_env_update_rho(p_env, qs_env)

      IF (transf_type_in == 3) THEN
         DO ispin = 1, nspins

            ! Add Fock contribution to Lagrangian
            mat_out => fm_mo_out(ispin)%matrix
            mat_in => fm_mo(ispin)%matrix
            ! scale for the orbital energy differences for the diagonal elements
            mat_out%local_data(:, :) = mat_in%local_data(:, :)
            CALL cp_fm_get_info(matrix=mat_out, &
                                nrow_local=nrow_local, &
                                ncol_local=ncol_local, &
                                row_indices=row_indices, &
                                col_indices=col_indices)
            DO jjB = 1, ncol_local
               j_global = col_indices(jjB)
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  mat_out%local_data(iiB, jjB) = mat_out%local_data(iiB, jjB)* &
                                                 (Eigenval(j_global + homo(ispin), ispin) - Eigenval(i_global, ispin))
               END DO
            END DO
         END DO
      END IF

      DO ispin = 1, nspins
         CALL dbcsr_set(p_env%kpp1(ispin)%matrix, 0.0_dp)
      END DO

      NULLIFY (input)

      rho1 => p_env%rho1
      rho1_xc => p_env%rho1_xc

      CALL get_qs_env(qs_env=qs_env, &
                      input=input, &
                      linres_control=linres_control)

      lrigpw = dft_control%qs_control%lrigpw

      lr_triplet = linres_control%lr_triplet

      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF

      CALL calc_kpp1(rho1_xc, rho1, xc_section, .FALSE., &
                     .FALSE., lrigpw, .TRUE., lr_triplet, &
                     qs_env, p_env)

      IF (do_hfx) THEN
         ! For HFX, we have to remap the different matrices
         ALLOCATE (rho_work_ao(nspins), mat_mu_nu(nspins))
         IF (dft_control%do_admm) ALLOCATE (rho_work_ao_admm(nspins))
         IF (nspins == 2) THEN
            ! In case of unrestricted calculations, we need a dummy dbcsr matrix for the hfx kernel
            NULLIFY (mat_mu_nu(2)%matrix)
            ALLOCATE (mat_mu_nu(2)%matrix)
            CALL dbcsr_create(mat_mu_nu(2)%matrix, template=p_env%kpp1(1)%matrix)
            rho_work_ao(2)%matrix => mat_mu_nu(2)%matrix
            IF (dft_control%do_admm) rho_work_ao_admm(2)%matrix => mat_mu_nu(2)%matrix
         END IF

         DO ispin = 1, nspins
            ! update with the exchange like contributions

            rho_work_ao(1)%matrix => p_env%p1(ispin)%matrix
            IF (dft_control%do_admm) rho_work_ao_admm(1)%matrix => p_env%p1_admm(ispin)%matrix
            mat_mu_nu(1)%matrix => p_env%kpp1(ispin)%matrix
            CALL apply_hfx_ao(qs_env, rho_work_ao, rho_work_ao_admm, mat_mu_nu, my_recalc_hfx_integrals)
         END DO

         IF (nspins == 2) THEN
            CALL dbcsr_deallocate_matrix(mat_mu_nu(2)%matrix)
         END IF
         DEALLOCATE (mat_mu_nu, rho_work_ao)
         IF (dft_control%do_admm) THEN
            DEALLOCATE (rho_work_ao_admm)
            CALL apply_xc_admm_ao(qs_env, p_env)
         END IF
      END IF

      DO ispin = 1, nspins
         CALL dbcsr_scale(p_env%kpp1(ispin)%matrix, 2.0_dp)

         ! copy back to fm
         CALL cp_fm_set_all(fm_ao, 0.0_dp)
         CALL copy_dbcsr_to_fm(matrix=p_env%kpp1(ispin)%matrix, fm=fm_ao)
         CALL cp_fm_set_all(fm_back, 0.0_dp)
         CALL cp_fm_upper_to_full(fm_ao, fm_back)

         mat_out => fm_mo_out(ispin)%matrix

         ! transform to MO basis, here we always sum the result into the input matrix

         SELECT CASE (transf_type_out)
            ! In alpha-beta case the density is contracted with the orbitals of different
            ! spin. Thus, there's branching at forward transform and no branching here,
            ! at back transform.
         CASE (1)
            ! occ-virt block
            CALL cp_gemm('T', 'N', homo(ispin), dimen, dimen, 1.0_dp, &
                         mo_coeff_o(ispin)%matrix, fm_ao, 0.0_dp, fm_back, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)
            CALL cp_gemm('N', 'N', homo(ispin), virtual(ispin), dimen, out_alpha, &
                         fm_back, mo_coeff_v(ispin)%matrix, 1.0_dp, mat_out, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)

         CASE (2)
            ! occ-occ block
            CALL cp_gemm('T', 'N', homo(ispin), dimen, dimen, 1.0_dp, &
                         mo_coeff_o(ispin)%matrix, fm_ao, 0.0_dp, fm_back, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)
            CALL cp_gemm('N', 'N', homo(ispin), homo(ispin), dimen, out_alpha, &
                         fm_back, mo_coeff_o(ispin)%matrix, 1.0_dp, mat_out, &
                         a_first_col=1, &
                         a_first_row=1, &
                         b_first_col=1, &
                         b_first_row=1, &
                         c_first_col=1, &
                         c_first_row=1)

         CASE DEFAULT
            ! nothing
         END SELECT

      END DO

      CALL timestop(handle)

   END SUBROUTINE cphf_like_update

! **************************************************************************************************
!> \brief Low level subroutine for the iterative solution of a large
!>        system of linear equation
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen ...
!> \param unit_nr ...
!> \param nspins ...
!> \param mo_coeff ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param Eigenval ...
!> \param p_env ...
!> \param blacs_env ...
!> \param do_hfx ...
!> \param L_jb ...
!> \param fm_G_mu_nu ...
!> \param fm_back ...
!> \param P_ia ...
!> \author Mauro Del Ben, Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE solve_z_vector_eq_low(qs_env, mp2_env, para_env, homo, virtual, dimen, unit_nr, nspins, &
                                    mo_coeff, mo_coeff_o, mo_coeff_v, Eigenval, p_env, blacs_env, &
                                    do_hfx, L_jb, fm_G_mu_nu, fm_back, P_ia)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(mp2_type), INTENT(IN), POINTER                :: mp2_env
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: nspins, unit_nr, dimen
      INTEGER, DIMENSION(nspins), INTENT(IN)             :: virtual, homo
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(nspins), INTENT(IN)  :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      TYPE(qs_p_env_type), INTENT(IN), POINTER           :: p_env
      TYPE(cp_blacs_env_type), INTENT(IN), POINTER       :: blacs_env
      LOGICAL, INTENT(IN)                                :: do_hfx
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: L_jb
      TYPE(cp_fm_type), INTENT(INOUT), POINTER           :: fm_G_mu_nu, fm_back
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: P_ia

      CHARACTER(LEN=*), PARAMETER :: routineN = 'solve_z_vector_eq_low'

      INTEGER :: cycle_counter, handle, i_global, iiB, iiter, ispin, j_global, jjB, max_num_iter, &
         ncol_local, nrow_local, transf_type_in, transf_type_out
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: eps_conv, out_alpha, t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: conv, proj_bi_xj, temp_vals, x_norm, xi_b
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: A_small, b_small, xi_Axi
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: b_i, precond, residual
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: Ax, xn
      TYPE(cp_fm_struct_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: fm_struct_tmp

      CALL timeset(routineN, handle)

      max_num_iter = mp2_env%ri_grad%cphf_max_num_iter
      eps_conv = mp2_env%ri_grad%cphf_eps_conv

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T3,A)') 'MP2_CPHF| Iterative solution of Z-Vector equations'
         WRITE (unit_nr, '(T3,A,T45,ES8.1)') 'MP2_CPHF| Convergence threshold:', eps_conv
         WRITE (unit_nr, '(T3,A,T45,I8)') 'MP2_CPHF| Maximum number of iterations: ', max_num_iter
         WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
         WRITE (unit_nr, '(T4,A,T15,A,T33,A)') 'Step', 'Time', 'Convergence'
         WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
      END IF

      ! set the transformation type (equal for all methods all updates)
      transf_type_in = 3
      transf_type_out = 1
      out_alpha = 1.0_dp

      ! set convergece flag
      converged = .FALSE.

      ! Pople method
      ! change sign to L_jb
      DO ispin = 1, nspins
         L_jb(ispin)%matrix%local_data(:, :) = -L_jb(ispin)%matrix%local_data(:, :)
      END DO

      ! allocate stuff
      ALLOCATE (xn(nspins, max_num_iter))
      ALLOCATE (Ax(nspins, max_num_iter))

      ! create fm structure
      ALLOCATE (fm_struct_tmp(nspins), precond(nspins), b_i(nspins), residual(ispin))
      DO ispin = 1, nspins
         NULLIFY (fm_struct_tmp(ispin)%struct)
         CALL cp_fm_struct_create(fm_struct_tmp(ispin)%struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=homo(ispin), ncol_global=virtual(ispin))

         ! create preconditioner (for now only orbital energy differences)
         NULLIFY (precond(ispin)%matrix)
         CALL cp_fm_create(precond(ispin)%matrix, fm_struct_tmp(ispin)%struct, name="precond")
         CALL cp_fm_set_all(precond(ispin)%matrix, 1.0_dp)
         CALL cp_fm_get_info(matrix=precond(ispin)%matrix, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               precond(ispin)%matrix%local_data(iiB, jjB) = precond(ispin)%matrix%local_data(iiB, jjB)/ &
                                                            (Eigenval(j_global + homo(ispin), ispin) - Eigenval(i_global, ispin))
            END DO
         END DO

         NULLIFY (b_i(ispin)%matrix)
         CALL cp_fm_create(b_i(ispin)%matrix, fm_struct_tmp(ispin)%struct, name="b_i")
         CALL cp_fm_set_all(b_i(ispin)%matrix, 0.0_dp)
         b_i(ispin)%matrix%local_data(:, :) = precond(ispin)%matrix%local_data(:, :)*L_jb(ispin)%matrix%local_data(:, :)

         ! create the residual vector (r), we check convergence on the norm of
         ! this vector r=(Ax-b)
         NULLIFY (residual(ispin)%matrix)
         CALL cp_fm_create(residual(ispin)%matrix, fm_struct_tmp(ispin)%struct, name="residual")
         CALL cp_fm_set_all(residual(ispin)%matrix, 0.0_dp)
      END DO

      ALLOCATE (conv(nspins))

      ! allocate array containing the various scalar products
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
            NULLIFY (xn(ispin, iiter)%matrix)
            CALL cp_fm_create(xn(ispin, iiter)%matrix, fm_struct_tmp(ispin)%struct, name="xi")
            CALL cp_fm_set_all(xn(ispin, iiter)%matrix, 0.0_dp)
         END DO

         ALLOCATE (proj_bi_xj(iiter - 1))
         proj_bi_xj = 0.0_dp
         ! first compute the projection of the actual b_i into all previous x_i
         ! already scaled with the norm of each x_i
         DO iiB = 1, iiter - 1
            DO ispin = 1, nspins
               proj_bi_xj(iiB) = proj_bi_xj(iiB) + accurate_dot_product(b_i(ispin)%matrix%local_data(:, :), &
                                                                        xn(ispin, iiB)%matrix%local_data(:, :))
            END DO
            proj_bi_xj(iiB) = proj_bi_xj(iiB)/x_norm(iiB)
         END DO

         CALL mp_sum(proj_bi_xj, para_env%group)

         ! update actual x_i
         DO ispin = 1, nspins
            xn(ispin, iiter)%matrix%local_data(:, :) = b_i(ispin)%matrix%local_data(:, :)
            DO iiB = 1, iiter - 1
               xn(ispin, iiter)%matrix%local_data(:, :) = xn(ispin, iiter)%matrix%local_data(:, :) - &
                                                          xn(ispin, iiB)%matrix%local_data(:, :)*proj_bi_xj(iiB)
            END DO
         END DO
         DEALLOCATE (proj_bi_xj)

         ! create Ax(iiter) that will store the matrix vector product for this cycle
         DO ispin = 1, nspins
            NULLIFY (Ax(ispin, iiter)%matrix)
            CALL cp_fm_create(Ax(ispin, iiter)%matrix, fm_struct_tmp(ispin)%struct, name="Ai")
            CALL cp_fm_set_all(Ax(ispin, iiter)%matrix, 0.0_dp)
         END DO

         CALL cphf_like_update(qs_env, homo, virtual, dimen, nspins, &
                               mo_coeff, mo_coeff_o, &
                               mo_coeff_v, Eigenval, p_env, &
                               do_hfx, &
                               xn(:, iiter), fm_G_mu_nu, fm_back, transf_type_in, out_alpha, &
                               Ax(:, iiter), transf_type_out)

         ! in order to reduce the number of calls to mp_sum here we
         ! cluster all necessary scalar products into a single vector
         ! temp_vals contains:
         ! 1:iiter -> <Ax_i|x_j>
         ! iiter+1 -> <x_i|b>
         ! iiter+2 -> <x_i|x_i>

         ALLOCATE (temp_vals(iiter + 2))
         temp_vals = 0.0_dp
         ! <Ax_i|x_j>
         DO ispin = 1, nspins
         DO iiB = 1, iiter
            temp_vals(iiB) = temp_vals(iiB) + &
                             accurate_dot_product(Ax(ispin, iiter)%matrix%local_data(:, :), &
                                                  xn(ispin, iiB)%matrix%local_data(:, :))
         END DO
         ! <x_i|b>
         temp_vals(iiter + 1) = temp_vals(iiter + 1) + accurate_dot_product(xn(ispin, iiter)%matrix%local_data(:, :), &
                                                                            L_jb(ispin)%matrix%local_data(:, :))
         ! norm
         temp_vals(iiter + 2) = temp_vals(iiter + 2) + accurate_dot_product(xn(ispin, iiter)%matrix%local_data(:, :), &
                                                                            xn(ispin, iiter)%matrix%local_data(:, :))
         END DO
         CALL mp_sum(temp_vals, para_env%group)
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
         conv = 0.0_dp
         DO ispin = 1, nspins
            CALL cp_fm_set_all(residual(ispin)%matrix, 0.0_dp)
            DO iiB = 1, iiter
               residual(ispin)%matrix%local_data(:, :) = &
                  residual(ispin)%matrix%local_data(:, :) + &
                  b_small(iiB, 1)*Ax(ispin, iiB)%matrix%local_data(:, :)
            END DO

            residual(ispin)%matrix%local_data(:, :) = &
               residual(ispin)%matrix%local_data(:, :) - &
               L_jb(ispin)%matrix%local_data(:, :)
            conv(ispin) = accurate_dot_product(residual(ispin)%matrix%local_data(:, :), &
                                               residual(ispin)%matrix%local_data(:, :))
         END DO

         CALL mp_sum(conv, para_env%group)
         conv(1) = SUM(SQRT(conv))

         t2 = m_walltime()

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,I5,T13,F6.1,11X,F14.8)') iiter, t2 - t1, conv(1)
         END IF

         IF (conv(1) <= eps_conv) THEN
            converged = .TRUE.
            EXIT
         END IF

         ! update b_i for the next round
         DO ispin = 1, nspins
            b_i(ispin)%matrix%local_data(:, :) = b_i(ispin)%matrix%local_data(:, :) &
                                                 + precond(ispin)%matrix%local_data(:, :) &
                                                 *Ax(ispin, iiter)%matrix%local_data(:, :)
         END DO

      END DO

      ! store solution into P_ia
      DO iiter = 1, cycle_counter
         DO ispin = 1, nspins
            P_ia(ispin)%matrix%local_data(:, :) = P_ia(ispin)%matrix%local_data(:, :) + &
                                                  b_small(iiter, 1)*xn(ispin, iiter)%matrix%local_data(:, :)
         END DO
      END DO

      ! Release arrays
      DEALLOCATE (x_norm)
      DEALLOCATE (xi_b)
      DEALLOCATE (xi_Axi)

      DO ispin = 1, nspins
         CALL cp_fm_release(precond(ispin)%matrix)
         CALL cp_fm_release(b_i(ispin)%matrix)
         CALL cp_fm_release(residual(ispin)%matrix)
         CALL cp_fm_struct_release(fm_struct_tmp(ispin)%struct)

         DO iiter = 1, cycle_counter
            CALL cp_fm_release(Ax(ispin, iiter)%matrix)
            CALL cp_fm_release(xn(ispin, iiter)%matrix)
         END DO
      END DO
      DEALLOCATE (precond, b_i, residual, fm_struct_tmp, Ax, xn)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T4,A)') REPEAT("-", 40)
         IF (converged) THEN
            WRITE (unit_nr, '(T3,A,I5,A)') 'Z-Vector equations converged in', cycle_counter, ' steps'
         ELSE
            WRITE (unit_nr, '(T3,A,I5,A)') 'Z-Vector equations NOT converged in', cycle_counter, ' steps'
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE solve_z_vector_eq_low

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE update_mp2_forces(qs_env, para_env)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_mp2_forces'

      INTEGER                                            :: alpha, beta, handle, handle2, i, ikind, &
                                                            ispin, nspins
      INTEGER, DIMENSION(3)                              :: comp
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: e_hartree, e_xc, pair_energy, tot_rho_r
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_dummy, h_stress
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_p_mp2, matrix_s, &
                                                            matrix_w_mp2, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: dvg(3), pot_g, pot_r, rho_tot_g, &
                                                            rho_tot_r, temp_pw_g
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:)         :: tau_mp2_g, tau_mp2_r
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_mp2_g, rho_mp2_r, rho_r, &
                                                            tau_r, v_rho, v_tau, v_xc
      TYPE(pw_p_type), POINTER                           :: rho_core
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, xc_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (input, pw_env, matrix_s, rho, energy, force, virial, matrix_w_mp2, &
               matrix_p_mp2, matrix_ks, rho_core)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env, &
                      input=input, &
                      matrix_s=matrix_s, &
                      matrix_ks=matrix_ks, &
                      matrix_p_mp2=matrix_p_mp2, &
                      matrix_w_mp2=matrix_w_mp2, &
                      rho=rho, &
                      energy=energy, &
                      force=force, &
                      virial=virial, &
                      rho_core=rho_core)

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      nspins = SIZE(rho_ao)

      ! check if we have to calculate the virial
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      ! pw stuff
      NULLIFY (poisson_env, auxbas_pw_pool)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ! get some of the grids ready
      NULLIFY (pot_r%pw, pot_g%pw, rho_tot_r%pw, rho_tot_g%pw, rho_mp2_g, rho_mp2_r)
      ALLOCATE (rho_mp2_r(nspins), rho_mp2_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (rho_mp2_r(ispin)%pw, rho_mp2_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_mp2_r(ispin)%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_mp2_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_zero(rho_mp2_r(ispin)%pw)
         CALL pw_zero(rho_mp2_g(ispin)%pw)
      END DO
      CALL pw_pool_create_pw(auxbas_pw_pool, pot_r%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, pot_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_r%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL pw_zero(rho_tot_r%pw)
      CALL pw_zero(rho_tot_g%pw)

      ! update the core-forces with the MP2-density contribution
      ! put MP2 density on the grid
      DO ispin = 1, nspins
         CALL calculate_rho_elec(matrix_p=matrix_p_mp2(ispin)%matrix, &
                                 rho=rho_mp2_r(ispin), &
                                 rho_gspace=rho_mp2_g(ispin), &
                                 total_rho=tot_rho_r, &
                                 ks_env=ks_env, &
                                 soft_valid=.FALSE.)

         CALL pw_axpy(rho_mp2_g(ispin)%pw, rho_tot_g%pw)
      END DO

      ! calculate the MP2 potential
      CALL pw_poisson_solve(poisson_env, rho_tot_g%pw, pair_energy, pot_g%pw)
      CALL pw_transfer(pot_g%pw, pot_r%pw)
      CALL pw_scale(pot_r%pw, pot_r%pw%pw_grid%dvol)

      ! calculate core forces
      CALL integrate_v_core_rspace(pot_r, qs_env)
      DO ikind = 1, SIZE(force)
         force(ikind)%mp2_sep = force(ikind)%rho_core
         force(ikind)%rho_core = 0.0_dp
      ENDDO

      ! Calculate XC part
      NULLIFY (v_xc)
      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO

      ! Get KS density

      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF
      IF (use_virial) h_stress = 0.0_dp
      CALL calc_ks_response_potentials(v_xc, xc_section, qs_env, rho_mp2_r, rho_mp2_g, compute_virial=use_virial, virial=h_stress)
      ! Update the virial
      IF (use_virial) THEN
         virial%pv_mp2 = virial%pv_mp2 + h_stress
         virial%pv_virial = virial%pv_virial + h_stress
      END IF
      !
      ! vtot = v_xc(ispin) + dv_hartree
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      DO ispin = 1, nspins
         ! Add v_hartree + v_xc = v_rspace
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(pot_r%pw, v_xc(ispin)%pw)
         ! integrate over potential <a|V|b>
         CALL integrate_v_rspace(v_rspace=v_xc(ispin), &
                                 hmat=matrix_ks(ispin), &
                                 pmat=rho_ao(ispin), &
                                 qs_env=qs_env, &
                                 calculate_forces=.TRUE.)

      END DO

      DO ikind = 1, SIZE(force)
         force(ikind)%mp2_sep = force(ikind)%mp2_sep + force(ikind)%rho_elec
         force(ikind)%rho_elec = 0.0_dp
      ENDDO

      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (v_xc)

      IF (use_virial) THEN
         ! update virial if necessary with the volume term
         ! first create pw auxiliary stuff
         CALL timeset(routineN//"_Virial", handle2)
         NULLIFY (temp_pw_g%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, temp_pw_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         DO i = 1, 3
            NULLIFY (dvg(i)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, dvg(i)%pw, &
                                   use_data=COMPLEXDATA1D, &
                                   in_space=RECIPROCALSPACE)
         END DO

         ! make a copy of the MP2 density in G space
         CALL pw_copy(rho_tot_g%pw, temp_pw_g%pw)
         ! calculate MP2-like-hartree potential derivatives
         DO i = 1, 3
            comp = 0
            comp(i) = 1
            CALL pw_copy(pot_g%pw, dvg(i)%pw)
            CALL pw_derive(dvg(i)%pw, comp)
         END DO

         CALL pw_zero(rho_tot_r%pw)
         CALL pw_zero(rho_tot_g%pw)

         NULLIFY (rho_r, rho_g)
         ALLOCATE (rho_r(nspins), rho_g(nspins))
         DO ispin = 1, nspins
            NULLIFY (rho_r(ispin)%pw, rho_g(ispin)%pw)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_r(ispin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_g(ispin)%pw, &
                                   use_data=COMPLEXDATA1D, &
                                   in_space=RECIPROCALSPACE)
         END DO

         ! calculate total SCF density and potential
         DO ispin = 1, nspins
            CALL calculate_rho_elec(matrix_p=rho_ao(ispin)%matrix, &
                                    rho=rho_r(ispin), &
                                    rho_gspace=rho_g(ispin), &
                                    total_rho=tot_rho_r, &
                                    ks_env=ks_env, &
                                    soft_valid=.FALSE.)

            CALL pw_axpy(rho_r(ispin)%pw, rho_tot_r%pw)
            CALL pw_axpy(rho_g(ispin)%pw, rho_tot_g%pw)
         END DO
         ! don't forget the core density
         CALL pw_axpy(rho_core%pw, rho_tot_g%pw)
         CALL pw_poisson_solve(poisson_env, rho_tot_g%pw, vhartree=pot_g%pw)

         ! finally update virial with the volume contribution
         e_hartree = pw_integral_ab(temp_pw_g%pw, pot_g%pw)
         h_stress = 0.0_dp
         DO alpha = 1, 3
            comp = 0
            comp(alpha) = 1
            CALL pw_copy(pot_g%pw, rho_tot_g%pw)
            CALL pw_derive(rho_tot_g%pw, comp)
            h_stress(alpha, alpha) = -e_hartree
            DO beta = alpha, 3
               h_stress(alpha, beta) = h_stress(alpha, beta) &
                                       - 2.0_dp*pw_integral_ab(rho_tot_g%pw, dvg(beta)%pw)/fourpi
               h_stress(beta, alpha) = h_stress(alpha, beta)
            END DO
         END DO

         ! free stuff
         CALL pw_pool_give_back_pw(auxbas_pw_pool, temp_pw_g%pw)
         DO i = 1, 3
            CALL pw_pool_give_back_pw(auxbas_pw_pool, dvg(i)%pw)
         END DO

         ! Add the contributions of the XC kernel
         NULLIFY (v_rho, v_tau)
         CALL qs_rho_get(rho, tau_r=tau_r)
         CALL xc_vxc_pw_create(v_rho, v_tau, e_xc, rho_r, rho_g, tau_r, xc_section, auxbas_pw_pool, .FALSE., h_dummy)

         ! For tau-dependent functionals, we need the perturbed kinetic energy densities
         IF (ASSOCIATED(v_tau)) THEN
            ALLOCATE (tau_mp2_r(nspins))
            ALLOCATE (tau_mp2_g(nspins))
            DO ispin = 1, nspins
               NULLIFY (tau_mp2_r(ispin)%pw)
               NULLIFY (tau_mp2_g(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, tau_mp2_r(ispin)%pw, &
                                      use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_pool_create_pw(auxbas_pw_pool, tau_mp2_g(ispin)%pw, &
                                      use_data=COMPLEXDATA1D, &
                                      in_space=RECIPROCALSPACE)

               CALL calculate_rho_elec(matrix_p=matrix_p_mp2(ispin)%matrix, &
                                       rho=tau_mp2_r(ispin), &
                                       rho_gspace=tau_mp2_g(ispin), &
                                       total_rho=tot_rho_r, &
                                       ks_env=ks_env, &
                                       soft_valid=.FALSE., &
                                       compute_tau=.TRUE.)

            END DO
         END IF

         e_xc = 0.0_dp
         DO ispin = 1, nspins
            e_xc = e_xc + pw_integral_ab(rho_mp2_r(ispin)%pw, v_rho(ispin)%pw)
            IF (ASSOCIATED(v_tau)) e_xc = e_xc + pw_integral_ab(tau_mp2_r(ispin)%pw, v_tau(ispin)%pw)
         END DO
WRITE(*,*) "exc", e_xc

         DO alpha = 1, 3
            h_stress(alpha, alpha) = h_stress(alpha, alpha) - e_xc
         END DO

         ! Update the virial
         virial%pv_mp2 = virial%pv_mp2 + h_stress/REAL(para_env%num_pe, dp)
         virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe, dp)

         ! release stuff
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r(ispin)%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g(ispin)%pw)
         END DO
         DEALLOCATE (rho_r, rho_g)

         IF (ASSOCIATED(v_rho)) THEN
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rho(ispin)%pw)
            END DO
            DEALLOCATE (v_rho)
         END IF
         IF (ASSOCIATED(v_tau)) THEN
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau(ispin)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, tau_mp2_r(ispin)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, tau_mp2_g(ispin)%pw)
            END DO
            DEALLOCATE (v_tau)
         END IF
         CALL timestop(handle2)
      END IF

      DO ispin = 1, nspins
         CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, 1.0_dp)
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_g%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_mp2_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_mp2_g(ispin)%pw)
      END DO
      DEALLOCATE (rho_mp2_r, rho_mp2_g)

      CALL timestop(handle)

   END SUBROUTINE update_mp2_forces

END MODULE mp2_cphf
