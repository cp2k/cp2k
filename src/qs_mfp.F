!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE qs_mfp
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_vnl_giao
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_scale,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_linres_methods,               ONLY: linres_solver
   USE qs_linres_types,                 ONLY: linres_control_type,&
                                              vcd_env_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_vcd_ao,                       ONLY: hr_mult_by_delta_1d
   USE qs_vcd_utils,                    ONLY: vcd_read_restart,&
                                              vcd_write_restart

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
!$ USE OMP_LIB, ONLY: omp_lock_kind, &
!$                    omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_destroy_lock
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: mfp_build_operator_gauge_independent
   PUBLIC :: mfp_build_operator_gauge_dependent
   PUBLIC :: mfp_response
   PUBLIC :: mfp_aat
   PUBLIC :: multiply_by_position

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mfp'

   REAL(dp), DIMENSION(3, 3, 3), PARAMETER  :: Levi_Civita = RESHAPE((/ &
                                                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                          0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, &
                                                        0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), &
                                                                     (/3, 3, 3/))
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param vcd_env ...
!> \param qs_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE mfp_aat(vcd_env, qs_env)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mfp_aat'
      INTEGER, PARAMETER                                 :: ispin = 1
      REAL(dp), DIMENSION(3), PARAMETER                  :: gauge_origin = 0._dp
      REAL(dp), PARAMETER                                :: f_spin = 2._dp

      INTEGER                                            :: alpha, delta, gamma, handle, ikind, nao, &
                                                            nmo
      LOGICAL                                            :: ghost
      REAL(dp)                                           :: aat_linmom, aat_moment, aat_moment_der, &
                                                            aat_overlap, aat_tmp, charge, lc_tmp, &
                                                            tmp_trace
      TYPE(cp_fm_type)                                   :: buf, buf2, matrix_dSdB_mo
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! Some setup
      nmo = vcd_env%dcdr_env%nmo(ispin)
      nao = vcd_env%dcdr_env%nao

      ASSOCIATE (mo_coeff => vcd_env%dcdr_env%mo_coeff(ispin), aat_atom => vcd_env%aat_atom_mfp)

         CALL get_qs_env(qs_env=qs_env, &
                         sab_all=sab_all, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set)

         ! Set up the matrices
         CALL cp_fm_create(buf, vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)
         CALL cp_fm_create(buf2, vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)
         CALL cp_fm_create(matrix_dSdB_mo, vcd_env%dcdr_env%momo_fm_struct(ispin)%struct)

         ! First prepare the coefficients dCB_prime = -dCB - 0.5 S1_ij
         !   The first negative because we get the negative coefficients out of the linres solver
         !   The second term due to the occ-occ block of the U matrix in C1 = U * C0
         DO alpha = 1, 3
            ! CALL cp_fm_scale_and_add(0._dp, vcd_env%dCB_prime(alpha), -1._dp, vcd_env%dCB(alpha))
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_dSdB(alpha)%matrix, mo_coeff, &
                                         buf, ncol=nmo, alpha=1._dp, beta=0._dp)
            CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                               1.0_dp, mo_coeff, buf, &
                               0.0_dp, matrix_dSdB_mo)

            CALL parallel_gemm("N", "N", nao, nmo, nmo, &
                               -0.5_dp, mo_coeff, matrix_dSdB_mo, &
                               0.0_dp, vcd_env%dCB_prime(alpha))

         END DO

         ! The AATs in MFP consist of four terms
         !
         ! i)     (i CB_alpha) * CR_beta * S0
         ! ii)    (i CB_alpha) * C0 * < mu | d/dR_beta^\lambda nu >
         ! iii)  -1/c * ε_{alpha gamma delta} * C0 * C0 * < mu | r_delta | d/dR_beta^\lambda nu > * R^mu_gamma
         ! iv)   -1/c * ε_{alpha gamma delta} * C0 * CR_beta * < mu | r_delta | nu > * R^mu_gamma

         ! iii)  -1/c * ε_{alpha gamma delta} * C0 * C0 * < mu | r_delta | d/dR_beta^\lambda nu > * R^mu_gamma
         !       +1/c * ε_{alpha gamma delta} * C0 * C0 * < mu | r_delta ∂_beta | nu > * R^mu_gamma * delta_nu^lambda
         !
         !
         ! vcd_env%moments_der_right(delta, beta)%matrix
         !    = +/- < mu | r_delta ∂_beta | nu > * (nu == lambda)
         DO alpha = 1, 3
            aat_moment_der = 0._dp

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE

                  CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                  vcd_env%moments_der_right(delta, vcd_env%dcdr_env%beta)%matrix)
                  CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                            qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                            gauge_origin=gauge_origin)

                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, ncol=nmo)
                  CALL cp_fm_trace(buf, mo_coeff, tmp_trace)

                  aat_moment_der = aat_moment_der - lc_tmp*tmp_trace
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_moment_der
         END DO

         ! And additionally we have the reference dependence
         !  + ε_{alpha gamma delta} * P0 * < mu | R^mu_gamma * r_delta ∂_beta | nu > * delta_nu^lambda
         !  - ε_{alpha gamma delta} * P0 * < mu | r_delta ∂_beta | nu >              * R^G_gamma * delta_nu^lambda
         !  - ε_{alpha gamma delta} * P0 * < mu | R^mu_gamma * ∂_beta | nu >         * R^G_delta * delta_nu^lambda
         !  + ε_{alpha gamma delta} * P0 * < mu | ∂_beta | nu >                      * R^G_gamma * R^G_delta * delta_nu^lambda
         DO alpha = 1, 3
            aat_tmp = 0._dp
            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE

                  ! - < mu | r_delta ∂_beta | nu > * R^G_gamma * delta_nu^lambda
                  ! I have
                  !  vcd_env%moments_der_right(delta, beta)%matrix = -< mu | r_delta ∂_beta | nu > * (nu == lambda)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%moments_der_right(delta, vcd_env%dcdr_env%beta)%matrix, &
                                               mo_coeff, buf, ncol=nmo)
                  CALL cp_fm_trace(buf, mo_coeff, tmp_trace)
                  aat_tmp = aat_tmp + lc_tmp*tmp_trace*vcd_env%magnetic_origin_atom(gamma)

                  ! - < mu | R^mu_gamma * ∂_beta | nu > * R^G_delta * delta_nu^lambda
                  ! I have
                  !  dipvel_ao(beta) = < a | ∂_beta | b >
                  CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%dipvel_ao(vcd_env%dcdr_env%beta)%matrix)
                  CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                            qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE.)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)

                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, &
                                               ncol=nmo, alpha=1._dp, beta=0._dp)
                  CALL cp_fm_trace(buf, mo_coeff, tmp_trace)
                  aat_tmp = aat_tmp - lc_tmp*tmp_trace*vcd_env%spatial_origin_atom(delta)

                  ! + < mu | ∂_beta | nu > * R^G_gamma * R^G_delta * delta_nu^lambda
                  ! I have
                  ! dipvel_ao(beta) = < a | ∂_beta | b >
                  CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%dipvel_ao(vcd_env%dcdr_env%beta)%matrix)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, &
                                               ncol=nmo, alpha=1._dp, beta=0._dp)
                  CALL cp_fm_trace(buf, mo_coeff, tmp_trace)
                  aat_tmp = aat_tmp + lc_tmp*tmp_trace &
                            *vcd_env%spatial_origin_atom(delta)*vcd_env%magnetic_origin_atom(gamma)

               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         ! ii)    (i CB_alpha) * C0 * < mu | d/dR_beta^\lambda nu >
         !         = - (i CB_alpha) * C0 * < mu | ∂_beta | nu > * delta_nu^lambda
         !
         ! dipvel_ao = + < a | ∂ | b >, so I can take that and multiply with the delta
         CALL dbcsr_desymmetrize(vcd_env%dipvel_ao(vcd_env%dcdr_env%beta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
         CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                  sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)

         DO alpha = 1, 3
            CALL cp_fm_set_all(buf, 0.0_dp)
            aat_linmom = 0.0_dp

            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                         mo_coeff, buf, ncol=nmo, alpha=1._dp, beta=0._dp)
            CALL cp_fm_trace(buf, vcd_env%dCB(alpha), aat_linmom)

            ! the minus, because dipvel_ao is the derivative wrt to the electronic coordinate,
            !   but we need the nuclear coordinate here.
            aat_linmom = -f_spin*aat_linmom

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_linmom
         END DO

         DO alpha = 1, 3
            CALL cp_fm_set_all(buf, 0.0_dp)
            aat_linmom = 0.0_dp

            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, ncol=nmo)
            CALL cp_fm_trace(buf, vcd_env%dCB_prime(alpha), aat_linmom)

            ! the minus, because dipvel_ao is the derivative wrt to the electronic coordinate,
            !   but we need the nuclear coordinate here.
            aat_linmom = -f_spin*aat_linmom

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_linmom
         END DO

         ! The reference dependent dSdB term is
         !   i/2c (- < mu | nu > * R^mu_gamma + < mu | nu > * R^nu_gamma) * R^G_delta
         DO alpha = 1, 3
            CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, 0._dp)

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE
                  CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, 0._dp)
                  CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, 0._dp)
                  CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix)
                  CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix)

                  ! < mu | nu > * R^mu_gamma
                  CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                                            qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                            gauge_origin=gauge_origin)

                  ! < mu | nu > * R^nu_gamma
                  CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, &
                                            qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                            gauge_origin=gauge_origin)

                  ! - < mu | nu > * R^mu_gamma + < mu | nu > * R^nu_gamma
                  CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                                 vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, -1._dp, +1._dp)

                  ! and accumulate the results
                  CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, &
                                 vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                                 1._dp, vcd_env%spatial_origin_atom(delta)*lc_tmp/2._dp)
               END DO
            END DO

            ! With the reference dependent overlap derivative in vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix,
            !  we now build the corresponding coefficients as buf2 = -0.5 * C0 * S1(ref)
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, mo_coeff, &
                                         buf, ncol=nmo, alpha=1._dp, beta=0._dp)
            CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                               1.0_dp, mo_coeff, buf, &
                               0.0_dp, matrix_dSdB_mo)

            CALL parallel_gemm("N", "N", nao, nmo, nmo, &
                               -0.5_dp, mo_coeff, matrix_dSdB_mo, &
                               0.0_dp, buf2)

            ! vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix still is < a | ∂ | b > * (nu==lambda)
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, ncol=nmo)
            CALL cp_fm_trace(buf, buf2, aat_linmom)

            aat_linmom = -f_spin*aat_linmom

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_linmom

         END DO

         ! iv)   -1/c * ε_{alpha gamma delta} * C0 * CR_beta * < mu | r_delta | nu > * R^mu_gamma
         DO alpha = 1, 3
            aat_moment = 0._dp

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE

                  CALL dbcsr_desymmetrize(vcd_env%dcdr_env%moments(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
                  CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                            qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                            gauge_origin=gauge_origin)

                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, ncol=nmo)
                  CALL cp_fm_trace(buf, vcd_env%dcdr_env%dCR_prime(ispin), tmp_trace)

                  aat_moment = aat_moment - lc_tmp*tmp_trace
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_moment
         END DO

         ! And additionally we have the reference dependence
         ! - ε_{alpha gamma delta} * C0 * CR_beta * < mu | R^mu_gamma r_delta | nu >
         ! + ε_{alpha gamma delta} * C0 * CR_beta * < mu | r_delta | nu > * R^G_gamma
         ! + ε_{alpha gamma delta} * C0 * CR_beta * < mu | R^mu_gamma | nu > * R^G_delta
         ! - ε_{alpha gamma delta} * C0 * CR_beta * < mu | nu > * R^G_gamma * R^G_delta
         DO alpha = 1, 3
            aat_moment = 0._dp

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE

                  ! < mu | r_delta | nu > * R^G_gamma
                  CALL dbcsr_desymmetrize(vcd_env%dcdr_env%moments(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, ncol=nmo)
                  CALL cp_fm_trace(buf, vcd_env%dcdr_env%dCR_prime(ispin), tmp_trace)

                  aat_moment = aat_moment + lc_tmp*tmp_trace*vcd_env%magnetic_origin_atom(gamma)

                  ! < mu | R^mu_gamma | nu > * R^G_delta
                  CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
                  CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                            qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE.)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, ncol=nmo)
                  CALL cp_fm_trace(buf, vcd_env%dcdr_env%dCR_prime(ispin), tmp_trace)

                  aat_moment = aat_moment + lc_tmp*tmp_trace*vcd_env%spatial_origin_atom(delta)

                  ! < mu | nu > * R^G_gamma * R^G_delta
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_s(1)%matrix, mo_coeff, buf, ncol=nmo)
                  CALL cp_fm_trace(buf, vcd_env%dcdr_env%dCR_prime(ispin), tmp_trace)

                  aat_moment = aat_moment + lc_tmp*tmp_trace &
                               *vcd_env%magnetic_origin_atom(gamma)*vcd_env%spatial_origin_atom(delta)
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_moment
         END DO

         ! i)       2 * (iCB_alpha) * CR_beta * S0
         !
         DO alpha = 1, 3
            CALL cp_fm_set_all(buf, 0.0_dp)
            aat_overlap = 0.0_dp

            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%dCR_prime(ispin), buf, ncol=nmo)
            CALL cp_fm_trace(buf, vcd_env%dCB(alpha), aat_overlap)

            aat_overlap = f_spin*aat_overlap

            ! dCB_prime * dCR_prime = + iP1
            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_overlap
         END DO

         DO alpha = 1, 3
            CALL cp_fm_set_all(buf, 0.0_dp)
            aat_overlap = 0.0_dp

            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%dCR_prime(ispin), buf, ncol=nmo)
            CALL cp_fm_trace(buf, vcd_env%dCB_prime(alpha), aat_overlap)

            aat_overlap = f_spin*aat_overlap

            ! dCB_prime * dCR_prime = + iP1
            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_overlap
         END DO

         ! The reference dependent dSdB term is
         !   i/2c (- < mu | nu > * R^mu_gamma + < mu | nu > * R^nu_gamma) * R^G_delta
         DO alpha = 1, 3
            CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, 0._dp)

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE
                  CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, 0._dp)
                  CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, 0._dp)
                  CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix)
                  CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix)

                  ! < mu | nu > * R^mu_gamma
                  CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                                            qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                            gauge_origin=gauge_origin)

                  ! < mu | nu > * R^nu_gamma
                  CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, &
                                            qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                            gauge_origin=gauge_origin)

                  ! !! A = alpha*A + beta*B or
                  ! - < mu | nu > * R^mu_gamma + < mu | nu > * R^nu_gamma
                  CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                                 vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, -1._dp, +1._dp)

                  ! and accumulate the results
                  CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                                 1._dp, -vcd_env%spatial_origin_atom(delta)*lc_tmp/2._dp)
               END DO
            END DO

            ! With the reference dependent overlap derivative in vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix,
            !  we now build the corresponding coefficients as buf2 = -0.5 * C0 * S1(ref)
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, mo_coeff, &
                                         buf, ncol=nmo, alpha=1._dp, beta=0._dp)
            CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                               1.0_dp, mo_coeff, buf, &
                               0.0_dp, matrix_dSdB_mo)

            CALL parallel_gemm("N", "N", nao, nmo, nmo, &
                               -0.5_dp, mo_coeff, matrix_dSdB_mo, &
                               0.0_dp, buf2)

            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%dCR_prime(ispin), buf, ncol=nmo)
            CALL cp_fm_trace(buf, buf2, aat_overlap)

            aat_overlap = f_spin*aat_overlap

            ! dCB_prime * dCR_prime = + iP1
            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_overlap
         END DO

         ! Nuclear contribution
         CALL get_atomic_kind(particle_set(vcd_env%dcdr_env%lambda)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost)
         IF (.NOT. ghost) THEN
            DO alpha = 1, 3
               aat_tmp = 0._dp
               DO gamma = 1, 3
                  IF (Levi_Civita(alpha, gamma, vcd_env%dcdr_env%beta) == 0._dp) CYCLE
                  aat_tmp = aat_tmp + charge &
                            *Levi_Civita(alpha, gamma, vcd_env%dcdr_env%beta) &
                            *(particle_set(vcd_env%dcdr_env%lambda)%r(gamma) - vcd_env%magnetic_origin_atom(gamma))

                  aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
                     = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp/4.
               END DO
            END DO
         END IF

      END ASSOCIATE

      CALL cp_fm_release(buf)
      CALL cp_fm_release(buf2)
      CALL cp_fm_release(matrix_dSdB_mo)

      CALL timestop(handle)
   END SUBROUTINE mfp_aat

! **************************************************************************************************
!> \brief ...
!> \param vcd_env ...
!> \param qs_env ...
!> \param alpha ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE mfp_build_operator_gauge_dependent(vcd_env, qs_env, alpha)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: alpha

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mfp_build_operator_gauge_dependent'
      INTEGER, PARAMETER                                 :: ispin = 1

      INTEGER                                            :: delta, gamma, handle, nao, nmo
      REAL(dp)                                           :: eps_ppnl, lc_tmp
      REAL(dp), DIMENSION(3)                             :: gauge_origin
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: buf
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      ! --------------------------------------------------------------- !
      ! The operator consists of four parts:
      !  1. R^mu x r * H^KS
      !  2. H^KS * R^nu x r
      !  3. [Vnl, (R^mu - Omag) x r]
      !  4. (r - Omag) x ∂
      ! and additionally the overlap derivative comes in as
      !  5. -dSdB * C0 * CHC
      ! --------------------------------------------------------------- !

      CALL timeset(routineN, handle)

      ! Some setup
      nmo = vcd_env%dcdr_env%nmo(ispin)
      nao = vcd_env%dcdr_env%nao
      ASSOCIATE (mo_coeff => vcd_env%dcdr_env%mo_coeff(ispin))

         CALL get_qs_env(qs_env=qs_env, &
                         sab_all=sab_all, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set, &
                         sap_ppnl=sap_ppnl, &
                         cell=cell, &
                         dft_control=dft_control, &     ! For the eps_ppnl
                         matrix_ks=matrix_ks)

         gauge_origin(:) = vcd_env%magnetic_origin_atom

         eps_ppnl = dft_control%qs_control%eps_ppnl

         CALL cp_fm_create(buf, vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)
         CALL cp_fm_set_all(vcd_env%op_dB(ispin), 0._dp)
         CALL cp_fm_set_all(buf, 0._dp)

         ! The first two terms in the operator are
         !    i/2c sum_{gamma delta} ε_{alpha gamma delta}
         !       (R^mu_gamma - O^mag_gamma) * < (r_delta - O^spatial) * H >
         !     - (R^nu_gamma - O^mag_gamma) * < H * (r_delta - O^spatial) >)

         ! they expand into
         ! i)     R^mu_gamma * matrix_rh(delta) - R^nu_gamma * matrix_hr(delta)
         ! ii)  - O^mag_gamma * (matrix_rh(delta) - matrix_hr(delta))
         ! iii) - (R^mu_gamma - R^nu_gamma) * O^spatial_delta * matrix_ks

         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! i)
               ! R^mu_gamma * matrix_rh(delta) - R^nu_gamma * matrix_hr(delta)
               ! R^mu_gamma * matrix_rh(delta)
               CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%matrix_rh(ispin, delta)%matrix)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                         gauge_origin=gauge_origin)

               ! - R^nu_gamma * matrix_hr(delta)
               CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, vcd_env%matrix_hr(ispin, delta)%matrix)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, 1._dp, -1._dp)

               ! and accumulate the results
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 1._dp, lc_tmp/2._dp)

               ! ii)
               ! - O^mag_gamma * (matrix_rh(delta) - matrix_hr(delta))
               CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%matrix_rh(ispin, delta)%matrix)
               CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, vcd_env%matrix_hr(ispin, delta)%matrix)
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, 1._dp, -1._dp)
               CALL dbcsr_scale(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, -vcd_env%magnetic_origin_atom(gamma))

               ! and accumulate the results
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 1._dp, lc_tmp/2._dp)

               ! iii)
               ! - (R^mu_gamma - R^nu_gamma) * O^spatial_delta * matrix_ks
               CALL dbcsr_desymmetrize(matrix_ks(ispin)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL dbcsr_desymmetrize(matrix_ks(ispin)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                         gauge_origin=gauge_origin)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                              1._dp, -1._dp)
               CALL dbcsr_scale(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, -vcd_env%spatial_origin_atom(delta))
               ! and accumulate the results
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 1._dp, lc_tmp/2._dp)

            END DO
         END DO

         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1._dp, beta=1.0_dp)

         ! The third term in the operator is
         !    i/2c sum_{gamma delta} ε_{alpha gamma delta}
         !       sum_\eta  < (R^eta_gamma - O^mag) * [Vnl^eta, r_delta] >
         !
         ! build_com_vnl_giao calculates
         !   matrix_rv(gamma, delta) =
         !       sum_\eta  < mu | R^eta_gamma * Vnl^eta * r_delta | nu >
         ! or the other way around: r_delta * Vnl^eta
         DO gamma = 1, 3
            DO delta = 1, 3
               CALL dbcsr_set(vcd_env%matrix_nosym_temp_33(gamma, delta)%matrix, 0._dp)
               CALL dbcsr_set(vcd_env%matrix_nosym_temp2_33(gamma, delta)%matrix, 0._dp)
            END DO
         END DO
         ! V * r_delta * R^eta_gamma
         CALL build_com_vnl_giao(qs_kind_set=qs_kind_set, sab_all=sab_all, sap_ppnl=sap_ppnl, &
                                 eps_ppnl=dft_control%qs_control%eps_ppnl, &
                                 particle_set=particle_set, matrix_rv=vcd_env%matrix_nosym_temp_33, &
                                 ref_point=[0._dp, 0._dp, 0._dp], cell=cell, direction_Or=.TRUE.)

         ! r_delta * V * R^eta_gamma
         CALL build_com_vnl_giao(qs_kind_set=qs_kind_set, sab_all=sab_all, sap_ppnl=sap_ppnl, &
                                 eps_ppnl=dft_control%qs_control%eps_ppnl, &
                                 particle_set=particle_set, matrix_rv=vcd_env%matrix_nosym_temp2_33, &
                                 ref_point=[0._dp, 0._dp, 0._dp], cell=cell, direction_Or=.FALSE.)

         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! + V * r_delta * R^eta_gamma
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, vcd_env%matrix_nosym_temp_33(gamma, delta)%matrix, &
                              1._dp, lc_tmp)

               ! - r_delta * V * R^eta_gamma
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, vcd_env%matrix_nosym_temp2_33(gamma, delta)%matrix, &
                              1._dp, -lc_tmp)
            END DO
         END DO

         ! Same factor 1/2 as above for Hr - rH
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1._dp/2._dp, beta=1.0_dp)

         ! And the O^mag dependent term:
         !    - sum_\eta  < O^mag * [Vnl^eta, r_delta] >
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! vcd_env%hcom(delta) = - < mu | [V, r_delta] | nu >
               CALL dbcsr_desymmetrize(vcd_env%hcom(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              1._dp, vcd_env%magnetic_origin_atom(gamma)*lc_tmp/(2._dp))
            END DO
         END DO
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1.0_dp, beta=1.0_dp)

         ! The fourth term in the operator is
         !  -i/2c sum_{gamma delta} ε_{alpha gamma delta}
         !      (r_gamma - O^mag_gamma) * ∂_delta

         ! It's also the P^1 term of the NVPT AAT with the reversed sign
         !  So in the first term alpha=+1
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! moments_der(gamma, delta) = - < a | r_gamma \partial_delta | b >
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, vcd_env%moments_der(gamma, delta)%matrix, &
                              1._dp, lc_tmp/(2._dp))

            END DO
         END DO
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1.0_dp, beta=1.0_dp)

         ! And the O^mag dependent term
         !   + i/2c sum_{gamma delta} ε_{alpha gamma delta}
         !      O^mag_gamma * ∂_delta
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! dipvel_ao(delta) = + < a | ∂_delta | b >
               CALL dbcsr_desymmetrize(vcd_env%dipvel_ao(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              1._dp, vcd_env%magnetic_origin_atom(gamma)*lc_tmp/(2._dp))
            END DO
         END DO
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1.0_dp, beta=1.0_dp)

         ! And the overlap derivative
         !   i/2c * sum_{gamma delta} ε_{alpha gamma delta}
         !     < r_delta > * (R^mu_gamma - R^nu_gamma)
         ! The reference independent part
         CALL dbcsr_set(vcd_env%matrix_dSdB(alpha)%matrix, 0._dp)

         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 0._dp)
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, 0._dp)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%moments(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%moments(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix)

               ! moments * R^mu_gamma
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                         gauge_origin=gauge_origin)

               ! moments * R^nu_gamma
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, 1._dp, -1._dp)

               ! and accumulate the results
               CALL dbcsr_add(vcd_env%matrix_dSdB(alpha)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 1._dp, lc_tmp/2._dp)
            END DO
         END DO

         ! the overlap derivative goes in as -S(1,B) * C0 * E0
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_dSdB(alpha)%matrix, mo_coeff, &
                                      buf, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL parallel_gemm('N', 'N', nao, nmo, nmo, &
                            -1.0_dp, buf, vcd_env%dcdr_env%chc(ispin), &   ! the negative sign is here
                            1.0_dp, vcd_env%op_dB(ispin))

         ! And the reference dependent part of the overlap
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, 0._dp)
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, 0._dp)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix)

               ! < mu | nu > * R^mu_gamma
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                         gauge_origin=gauge_origin)

               ! < mu | nu > * R^nu_gamma
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               ! !! A = alpha*A + beta*B or
               ! - < mu | nu > * R^mu_gamma + < mu | nu > * R^nu_gamma
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, -1._dp, +1._dp)

               ! and accumulate the results
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                              1._dp, vcd_env%spatial_origin_atom(delta)*lc_tmp/2._dp)
            END DO
         END DO

         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, mo_coeff, &
                                      buf, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL parallel_gemm('N', 'N', nao, nmo, nmo, &
                            -1.0_dp, buf, vcd_env%dcdr_env%chc(ispin), &   ! again, the negative sign is here
                            1.0_dp, vcd_env%op_dB(ispin))

         ! CALL cp_fm_set_all(vcd_env%op_dB(ispin), 0._dp)
         CALL cp_fm_release(buf)

      END ASSOCIATE

      ! We have built op_dB but plug -op_dB into the linres_solver
      ! CALL cp_fm_scale(-1.0_dp, vcd_env%op_dB(ispin))

      CALL timestop(handle)
   END SUBROUTINE mfp_build_operator_gauge_dependent

! **************************************************************************************************
!> \brief ...
!> \param vcd_env ...
!> \param qs_env ...
!> \param alpha ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE mfp_build_operator_gauge_independent(vcd_env, qs_env, alpha)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: alpha

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mfp_build_operator_gauge_independent'
      INTEGER, PARAMETER                                 :: ispin = 1

      INTEGER                                            :: delta, gamma, handle, nao, nmo
      REAL(dp)                                           :: eps_ppnl, lc_tmp
      REAL(dp), DIMENSION(3)                             :: gauge_origin
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: buf
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      ! --------------------------------------------------------------- !
      ! The operator consists of three parts:
      !  1. (R^mu - R^nu) x r * H_KS
      !  2. [Vnl, (R^mu - R^nu) x r]
      !  3. (r - R^nu) x ∂
      ! and additionally the overlap derivative comes in as
      !  4. -dSdB * C0 * CHC
      ! --------------------------------------------------------------- !

      CALL timeset(routineN, handle)

      ! Some setup
      nmo = vcd_env%dcdr_env%nmo(ispin)
      nao = vcd_env%dcdr_env%nao
      ASSOCIATE (mo_coeff => vcd_env%dcdr_env%mo_coeff(ispin))

         CALL get_qs_env(qs_env=qs_env, &
                         sab_all=sab_all, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set, &
                         sap_ppnl=sap_ppnl, &
                         cell=cell, &
                         dft_control=dft_control, &     ! For the eps_ppnl
                         matrix_ks=matrix_ks)

         gauge_origin(:) = 0._dp

         eps_ppnl = dft_control%qs_control%eps_ppnl

         CALL cp_fm_create(buf, vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)
         CALL cp_fm_set_all(vcd_env%op_dB(ispin), 0._dp)
         CALL cp_fm_set_all(buf, 0._dp)

         ! The first term in the operator is
         !    i/2c sum_{gamma delta} ε_{alpha gamma delta}
         !       (R^mu_gamma - R^nu_gamma) < (r_delta - Ospatial) H >

         ! In vcd_second_prepare we build
         !       vcd_env%matrix_hr(1, delta) = < mu | H * r_delta | nu >
         ! and   vcd_env%matrix_rh(1, delta) = < mu | r_delta * H | nu >

         ! So we need
         !   (R^mu_gamma - R^nu_gamma) * vcd_env%matrix_rh
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! The reference independent part
               ! vcd_env%matrix_rh * R^mu_gamma
               CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%matrix_rh(ispin, delta)%matrix)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                         gauge_origin=gauge_origin)

               ! - vcd_env%matrix_rh * R^nu_gamma
               CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, vcd_env%matrix_rh(ispin, delta)%matrix)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, 1._dp, -1._dp)

               ! and accumulate the results
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 1._dp, lc_tmp/2._dp)

               ! And the Ospatial dependent part
               ! - (R^mu_gamma - R^nu_gamma) * Ospatial * matrix_ks
               CALL dbcsr_desymmetrize(matrix_ks(ispin)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL dbcsr_desymmetrize(matrix_ks(ispin)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                         gauge_origin=gauge_origin)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                              1._dp, -1._dp)
               CALL dbcsr_scale(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, -vcd_env%spatial_origin_atom(delta))
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 1._dp, lc_tmp/2._dp)

            END DO
         END DO

         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1._dp, beta=1.0_dp)

         ! The second term in the operator is
         !    i/2c sum_{gamma delta} ε_{alpha gamma delta}
         !       sum_\eta  < (R^eta_gamma - R^nu_gamma) * [Vnl^eta, r_delta] >
         !
         ! build_com_vnl_giao calculates
         !   matrix_rv(gamma, delta) =
         !       sum_\eta  < mu | R^eta_gamma * Vnl^eta * r_delta | nu >
         ! or the other way around: r_delta * Vnl^eta
         DO gamma = 1, 3
            DO delta = 1, 3
               CALL dbcsr_set(vcd_env%matrix_nosym_temp_33(gamma, delta)%matrix, 0._dp)
               CALL dbcsr_set(vcd_env%matrix_nosym_temp2_33(gamma, delta)%matrix, 0._dp)
            END DO
         END DO
         ! V * r_delta * R^eta_gamma
         CALL build_com_vnl_giao(qs_kind_set=qs_kind_set, sab_all=sab_all, sap_ppnl=sap_ppnl, &
                                 eps_ppnl=dft_control%qs_control%eps_ppnl, &
                                 particle_set=particle_set, matrix_rv=vcd_env%matrix_nosym_temp_33, &
                                 ref_point=[0._dp, 0._dp, 0._dp], cell=cell, direction_Or=.TRUE.)

         ! r_delta * V * R^eta_gamma
         CALL build_com_vnl_giao(qs_kind_set=qs_kind_set, sab_all=sab_all, sap_ppnl=sap_ppnl, &
                                 eps_ppnl=dft_control%qs_control%eps_ppnl, &
                                 particle_set=particle_set, matrix_rv=vcd_env%matrix_nosym_temp2_33, &
                                 ref_point=[0._dp, 0._dp, 0._dp], cell=cell, direction_Or=.FALSE.)

         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! + V * r_delta * R^eta_gamma
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, &
                              vcd_env%matrix_nosym_temp_33(gamma, delta)%matrix, &
                              1._dp, lc_tmp)

               ! - r_delta * V * R^eta_gamma
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, &
                              vcd_env%matrix_nosym_temp2_33(gamma, delta)%matrix, &
                              1._dp, -lc_tmp)
            END DO
         END DO

         ! Same factor 1/2 as above for Hr - rH
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1._dp/2._dp, beta=1.0_dp)

         ! The third term in the operator is
         !  -i/2c sum_{gamma delta} ε_{alpha gamma delta}
         !      (r_gamma - R^nu_gamma) * ∂_delta

         ! It's also the P^1 term of the NVPT AAT with the reversed sign
         !  So in the first term alpha=+1
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! moments_der(gamma, delta) = - < a | r_gamma \partial_delta | b >
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, vcd_env%moments_der(gamma, delta)%matrix, &
                              1._dp, lc_tmp/(2._dp))

            END DO
         END DO
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1.0_dp, beta=1.0_dp)

         ! And the R^nu dependent term
         !   + i/2c sum_{gamma delta} ε_{alpha gamma delta}
         !      R^nu_gamma * ∂_delta
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE

               ! dipvel_ao(delta) = + < a | ∂_delta | b >
               CALL dbcsr_desymmetrize(vcd_env%dipvel_ao(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              1._dp, lc_tmp/(2._dp))
            END DO
         END DO
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(3)%matrix, mo_coeff, &
                                      vcd_env%op_dB(ispin), ncol=nmo, alpha=1.0_dp, beta=1.0_dp)

         ! And the overlap derivative
         !   i/2c * sum_{gamma delta} ε_{alpha gamma delta}
         !     < r_delta > * (R^mu_gamma - R^nu_gamma)
         ! The reference independent part
         CALL dbcsr_set(vcd_env%matrix_dSdB(alpha)%matrix, 0._dp)

         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 0._dp)
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, 0._dp)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%moments(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%moments(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix)

               ! moments * R^mu_gamma
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                         gauge_origin=gauge_origin)

               ! moments * R^nu_gamma
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               ! !! A = alpha*A + beta*B or
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, 1._dp, -1._dp)

               ! and accumulate the results
               CALL dbcsr_add(vcd_env%matrix_dSdB(alpha)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 1._dp, lc_tmp/2._dp)
            END DO
         END DO

         ! the overlap derivative goes in as -S(1,B) * C0 * E0
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_dSdB(alpha)%matrix, mo_coeff, &
                                      buf, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL parallel_gemm('N', 'N', nao, nmo, nmo, &
                            -1.0_dp, buf, vcd_env%dcdr_env%chc(ispin), &   ! the negative sign is here
                            1.0_dp, vcd_env%op_dB(ispin))

         ! And the reference dependent part of the overlap
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, 0._dp)
         DO gamma = 1, 3
            DO delta = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, delta)
               IF (lc_tmp == 0._dp) CYCLE
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, 0._dp)
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, 0._dp)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix)

               ! < mu | nu > * R^mu_gamma
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.FALSE., &
                                         gauge_origin=gauge_origin)

               ! < mu | nu > * R^nu_gamma
               CALL multiply_by_position(vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, &
                                         qs_kind_set, particle_set, "ORB", sab_all, gamma, basis_function_nu=.TRUE., &
                                         gauge_origin=gauge_origin)

               ! !! A = alpha*A + beta*B or
               ! - < mu | nu > * R^mu_gamma + < mu | nu > * R^nu_gamma
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                              vcd_env%dcdr_env%matrix_nosym_temp2(2)%matrix, -1._dp, +1._dp)

               ! and accumulate the results
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, vcd_env%dcdr_env%matrix_nosym_temp2(1)%matrix, &
                              1._dp, vcd_env%spatial_origin_atom(delta)*lc_tmp/2._dp)
            END DO
         END DO

         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp2(3)%matrix, mo_coeff, &
                                      buf, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL parallel_gemm('N', 'N', nao, nmo, nmo, &
                            -1.0_dp, buf, vcd_env%dcdr_env%chc(ispin), &   ! again, the negative sign is here
                            1.0_dp, vcd_env%op_dB(ispin))

         ! CALL cp_fm_set_all(vcd_env%op_dB(ispin), 0._dp)
         CALL cp_fm_release(buf)

      END ASSOCIATE

      ! We have built op_dB but plug -op_dB into the linres_solver
      ! CALL cp_fm_scale(-1.0_dp, vcd_env%op_dB(ispin))

      CALL timestop(handle)
   END SUBROUTINE mfp_build_operator_gauge_independent

! *****************************************************************************
!> \brief Get the dC/dB using the vcd_env%op_dB
!> \param vcd_env ...
!> \param p_env ...
!> \param qs_env ...
!> \param alpha ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE mfp_response(vcd_env, p_env, qs_env, alpha)

      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: alpha

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mfp_response'

      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: failure, should_stop
      TYPE(cp_fm_type), DIMENSION(1)                     :: h1_psi0, psi1
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(section_vals_type), POINTER                   :: lr_section, vcd_section

      CALL timeset(routineN, handle)
      failure = .FALSE.

      NULLIFY (linres_control, lr_section, logger)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      linres_control=linres_control, &
                      mos=mos)

      logger => cp_get_default_logger()
      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      vcd_section => section_vals_get_subs_vals(qs_env%input, &
                                                "PROPERTIES%LINRES%VCD")

      output_unit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".linresLog")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T10,A,/)") &
            "*** Self consistent optimization of the magnetic response wavefunction ***"
      END IF

      ! allocate the vectors
      ASSOCIATE (psi0_order => vcd_env%dcdr_env%mo_coeff)
         CALL cp_fm_create(psi1(1), vcd_env%dcdr_env%likemos_fm_struct(1)%struct)
         CALL cp_fm_create(h1_psi0(1), vcd_env%dcdr_env%likemos_fm_struct(1)%struct)

         ! Restart
         IF (linres_control%linres_restart) THEN
            CALL vcd_read_restart(qs_env, lr_section, psi1, vcd_env%dcdr_env%lambda, alpha, "dCdB")
         ELSE
            CALL cp_fm_set_all(psi1(1), 0.0_dp)
         END IF

         IF (output_unit > 0) THEN
            WRITE (output_unit, *) &
               "Response to the perturbation operator referring to the magnetic field in "//ACHAR(alpha + 119)
         END IF

         ! First response to get dCR
         ! (H0-E0) psi1 = (H1-E1) psi0
         ! psi1 = the perturbed wavefunction
         ! h1_psi0 = (H1-E1)
         ! psi0_order = the unperturbed wavefunction
         ! Second response to get dCB
         CALL cp_fm_set_all(vcd_env%dCB(alpha), 0.0_dp)
         CALL cp_fm_set_all(h1_psi0(1), 0.0_dp)
         CALL cp_fm_to_fm(vcd_env%op_dB(1), h1_psi0(1))

         linres_control%lr_triplet = .FALSE. ! we do singlet response
         linres_control%do_kernel = .FALSE. ! no coupled response since imaginary perturbation
         linres_control%converged = .FALSE.
         CALL linres_solver(p_env, qs_env, psi1, h1_psi0, psi0_order, &
                            output_unit, should_stop)
         CALL cp_fm_to_fm(psi1(1), vcd_env%dCB(alpha))

         ! Write the new result to the restart file
         IF (linres_control%linres_restart) THEN
            CALL vcd_write_restart(qs_env, lr_section, psi1, vcd_env%dcdr_env%lambda, alpha, "dCdB")
         END IF

         ! clean up
         CALL cp_fm_release(psi1(1))
         CALL cp_fm_release(h1_psi0(1))
      END ASSOCIATE

      CALL cp_print_key_finished_output(output_unit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)
   END SUBROUTINE mfp_response

! **************************************************************************************************
!> \brief Take matrix < mu | ^O^ | nu > and multiply the blocks with the positions of the basis functions.
!>        The matrix consists of blocks (iatom, jatom) corresponding to (mu, nu).
!>        With basis_function_nu = .TRUE.  we have to multiply each block by particle_set(jatom)%r(delta)
!>        With basis_function_nu = .FALSE. we have to multiply each block by particle_set(iatom)%r(delta)
!> \param matrix              The matrix we are operating on
!> \param qs_kind_set         Needed to set up the basis_sets
!> \param particle_set        Needed for the coordinates
!> \param basis_type          Needed to set up the basis_sets
!> \param sab_nl              The supplied neighborlist
!> \param direction           The index of the nuclear position we are multiplying by
!> \param basis_function_nu   Determines whether we multiply by R^nu or R^mu
!> \param gauge_origin        The gauge origin, if we do distributed gauge
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE multiply_by_position(matrix, qs_kind_set, particle_set, basis_type, sab_nl, &
                                   direction, basis_function_nu, gauge_origin)

      TYPE(dbcsr_type)                                   :: matrix
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER                                            :: direction
      LOGICAL                                            :: basis_function_nu
      REAL(dp), DIMENSION(3), OPTIONAL                   :: gauge_origin

      CHARACTER(len=*), PARAMETER :: routineN = 'multiply_by_position'

      INTEGER                                            :: handle, iatom, icol, ikind, inode, irow, &
                                                            jatom, jkind, last_jatom, mepos, &
                                                            nkind, nseta, nsetb, nthread
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: do_symmetric, found
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix_block, rpgfa, rpgfb, scon_a, &
                                                            scon_b, zeta, zetb
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CALL timeset(routineN, handle)

      ! check for symmetry
      CPASSERT(SIZE(sab_nl) > 0)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)

      ! prepare basis set and coordinates
      nkind = SIZE(qs_kind_set)
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, basis_type, qs_kind_set)

      nthread = 1
!$    nthread = omp_get_max_threads()
      ! Iterate of neighbor list
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,nl_iterator, do_symmetric) &
!$OMP SHARED (matrix,basis_set_list) &
!$OMP SHARED (basis_function_nu, direction, particle_set, gauge_origin) &
!$OMP PRIVATE (matrix_block,mepos,ikind,jkind,iatom,jatom,cell) &
!$OMP PRIVATE (basis_set_a,basis_set_b) &
!$OMP PRIVATE (first_sgfa, la_max, la_min, npgfa, nsgfa, nseta, rpgfa, set_radius_a) &
!$OMP PRIVATE (zeta, first_sgfb, lb_max, lb_min, npgfb, nsetb, rpgfb, set_radius_b, nsgfb) &
!$OMP PRIVATE (zetb, scon_a, scon_b, irow, icol, found, inode, last_jatom)
      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, inode=inode)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         scon_a => basis_set_a%scon
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         scon_b => basis_set_b%scon
         zetb => basis_set_b%zet

         IF (inode == 1) last_jatom = 0

         ! this guarantees minimum image convention
         ! anything else would not make sense
         IF (jatom == last_jatom) THEN
            CYCLE
         END IF

         last_jatom = jatom

         IF (do_symmetric) THEN
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
         ELSE
            irow = iatom
            icol = jatom
         END IF

         NULLIFY (matrix_block)
         CALL dbcsr_get_block_p(matrix, irow, icol, matrix_block, found)
         CPASSERT(found)

         IF (PRESENT(gauge_origin)) THEN
            IF (basis_function_nu) THEN
               !$OMP CRITICAL(blockadd)
               matrix_block(:, :) = matrix_block(:, :)*(particle_set(jatom)%r(direction) - gauge_origin(direction))
               !$OMP END CRITICAL(blockadd)
            ELSE
               !$OMP CRITICAL(blockadd)
               matrix_block(:, :) = matrix_block(:, :)*(particle_set(iatom)%r(direction) - gauge_origin(direction))
               !$OMP END CRITICAL(blockadd)
            END IF
         ELSE IF (.NOT. PRESENT(gauge_origin)) THEN
            IF (basis_function_nu) THEN
               !$OMP CRITICAL(blockadd)
               matrix_block(:, :) = matrix_block(:, :)*particle_set(jatom)%r(direction)
               !$OMP END CRITICAL(blockadd)
            ELSE
               !$OMP CRITICAL(blockadd)
               matrix_block(:, :) = matrix_block(:, :)*particle_set(iatom)%r(direction)
               !$OMP END CRITICAL(blockadd)
            END IF
         END IF
      END DO
!$OMP END PARALLEL
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Release work storage
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE multiply_by_position

END MODULE qs_mfp
