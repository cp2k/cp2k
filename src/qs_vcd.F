!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE qs_vcd
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_mom_nl
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_set
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
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
   USE qs_dcdr_ao,                      ONLY: hr_mult_by_delta_1d
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_linres_methods,               ONLY: linres_solver
   USE qs_linres_types,                 ONLY: linres_control_type,&
                                              vcd_env_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_moments,                      ONLY: dipole_velocity_deriv
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_vcd_ao,                       ONLY: build_dSdV_matrix,&
                                              build_dcom_rpnl,&
                                              build_drpnl_matrix,&
                                              build_matrix_hr_rh,&
                                              hr_mult_by_delta_3d
   USE qs_vcd_utils,                    ONLY: vcd_read_restart,&
                                              vcd_write_restart
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: prepare_per_atom_vcd
   PUBLIC :: vcd_build_op_dV
   PUBLIC :: vcd_response_dV
   PUBLIC :: apt_dV
   PUBLIC :: aat_dV

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vcd'

   REAL(dp), DIMENSION(3, 3, 3), PARAMETER  :: Levi_Civita = RESHAPE((/ &
                                                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                          0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, &
                                                        0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), &
                                                                     (/3, 3, 3/))
   INTEGER, DIMENSION(3, 3), PARAMETER :: multipole_2d_to_1d = RESHAPE([4, 5, 6, 5, 7, 8, 6, 8, 9], [3, 3])
CONTAINS

! **************************************************************************************************
!> \brief Compute I_{alpha beta}^lambda = d/dV^lambda_beta <m_alpha> = d/dV^lambda_beta < r x \dot{r} >
!>        The directions alpha, beta are stored in vcd_env%dcdr_env
!> \param vcd_env ...
!> \param qs_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE aat_dV(vcd_env, qs_env)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'aat_dV'
      INTEGER, PARAMETER                                 :: ispin = 1

      INTEGER                                            :: alpha, delta, gamma, handle, ikind, &
                                                            my_index, nao, nmo, nspins
      LOGICAL                                            :: ghost
      REAL(dp)                                           :: aat_prefactor, aat_tmp, charge, lc_tmp, &
                                                            tmp_trace
      REAL(dp), DIMENSION(3, 3)                          :: aat_tmp_33
      TYPE(cp_fm_type)                                   :: tmp_aomo
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      sap_ppnl=sap_ppnl, &
                      sab_orb=sab_orb, &
                      sab_all=sab_all, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      CALL cp_fm_create(tmp_aomo, vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)

      nspins = dft_control%nspins
      nmo = vcd_env%dcdr_env%nmo(ispin)
      nao = vcd_env%dcdr_env%nao
      ASSOCIATE (mo_coeff => vcd_env%dcdr_env%mo_coeff(ispin), aat_atom => vcd_env%aat_atom_nvpt)

         ! I_{alpha beta}^lambda = 1/2c \sum_j^occ ...
         aat_prefactor = 1.0_dp!/(c_light_au * 2._dp)
         IF (nspins .EQ. 1) aat_prefactor = aat_prefactor*2.0_dp

         ! The non-PP part of the AAT consists of four contributions:
         !  (A1):  + P^0 * ε_{alpha gamma delta} * < mu | r_beta r_gamma ∂_delta | nu > * (mu == lambda)
         !  (A2):  - P^0 * ε_{alpha gamma delta} * < mu | r_gamma r_beta ∂_delta | nu > * (nu == lambda)
         !  (B):   - P^0 * ε_{alpha gamma delta} * < mu | r_gamma | nu > * (delta == beta) * (nu == lambda)
         !  (C):   + iP^1 * ε_{alpha gamma delta} * < mu | r_gamma ∂_delta | nu >

         ! (A1) + P^0 * ε_{alpha gamma delta} * < mu | r_beta r_gamma ∂_delta | nu > * (mu == lambda)
         ! (A2) - P^0 * ε_{alpha gamma delta} * < mu | r_gamma r_beta ∂_delta | nu > * (nu == lambda)
         ! Conjecture : It doesn't matter that the beta and gamma are swapped around!
         !               We define o = | ∂_delta nu >
         !                 and then < a | r_beta r_gamma | o > = < a | r_gamma r_beta | o>
         ! (A) + P^0 * ε_{alpha gamma delta} * < mu | r_beta r_gamma ∂_delta | nu > * (mu == lambda - nu == lambda)
         ! We have built the matrices - < mu | r_beta r_gamma ∂_delta | nu > in vcd_env%moments_der
         ! moments_der(1:9; 1:3) = moments_der(x, y, z, xx, xy, xz, yy, yz, zz;
         !                                     x, y, z)

         aat_tmp_33 = 0._dp
         DO gamma = 1, 3
            my_index = multipole_2d_to_1d(vcd_env%dcdr_env%beta, gamma)
            DO delta = 1, 3
               ! moments_der(moment, delta) = - < a | moment \partial_\delta | b >
               ! matrix_nosym_temp = - < mu | r_beta r_gamma ∂_delta | nu > * (mu - nu)
               CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                               vcd_env%moments_der_right(my_index, delta)%matrix)
               CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                              vcd_env%moments_der_left(my_index, delta)%matrix, &
                              1._dp, -1._dp)

               CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
               CALL cp_fm_trace(mo_coeff, tmp_aomo, aat_tmp_33(gamma, delta))
            END DO
         END DO

         DO alpha = 1, 3
            aat_tmp = 0._dp

            ! There are two remaining combinations for gamma and delta.
            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE

                  ! moments_der(moment, delta) = - < a | moment \partial_\delta | b >
                  ! matrix_nosym_temp = - < mu | r_beta r_gamma ∂_delta | nu > * (mu - nu)
                  ! Because of the negative in moments_der, we need another negative sign here.
                  aat_tmp = aat_tmp + lc_tmp*aat_prefactor*aat_tmp_33(gamma, delta)
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (B):   - P^0 * ε_{alpha gamma delta} * < mu | r_gamma | nu > * (delta == beta) * (nu == lambda)
         !      =  - P^0 * ε_{alpha gamma beta} * < mu | r_gamma | nu > * (nu == lambda)

         DO alpha = 1, 3
            aat_tmp = 0._dp

            DO gamma = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, vcd_env%dcdr_env%beta)
               IF (lc_tmp == 0._dp) CYCLE

               ! matrix_nosym_temp = < mu | r_gamma | nu > * (nu == lambda)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%moments(gamma)%matrix, &
                                       vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                        sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)

               CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
               CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
               aat_tmp = aat_tmp - lc_tmp*aat_prefactor*tmp_trace
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (C):   + iP^1 * ε_{alpha gamma delta} * < mu | r_gamma ∂_delta | nu >
         DO alpha = 1, 3
            aat_tmp = 0._dp

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE

                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%moments_der(gamma, delta)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(tmp_aomo, vcd_env%dCV_prime(ispin), tmp_trace)

                  ! mo_coeff * dCV_prime = + iP1
                  ! moments_der(moment, delta) = - < a | moment \partial_\delta | b >
                  ! so we need the opposite sign.
                  aat_tmp = aat_tmp - 2._dp*aat_prefactor*tmp_trace*lc_tmp
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         ! The PP part consists of four contributions
         !  (D):  - P^0 * ε_{alpha gamma delta} * < mu | r_beta r_gamma [V, r_delta] | nu > * (mu == lambda)
         !  (E):  + P^0 * ε_{alpha gamma delta} * < mu | r_gamma [V, r_delta] r_beta | nu > * (nu == lambda)
         !  (F):  - P^0 * ε_{alpha gamma delta} * < mu | r_gamma [[V, r_beta], r_delta] | nu > * (eta == lambda)
         !  (G):  - iP^1 * ε_{alpha gamma delta} * < mu | r_gamma [V, r_delta] | nu >

         !  (D):  - P^0 * ε_{alpha gamma delta} * < mu | r_beta r_gamma [V, r_delta] | nu > * (mu == lambda)
         !    The negative of this is in vcd_env%matrix_r_rxvr
         DO alpha = 1, 3
            CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                            vcd_env%matrix_r_rxvr(alpha, vcd_env%dcdr_env%beta)%matrix)
            CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                     sab_all, direction_Or=.FALSE., lambda=vcd_env%dcdr_env%lambda)

            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
            CALL cp_fm_trace(mo_coeff, tmp_aomo, aat_tmp)
            aat_tmp = -aat_prefactor*aat_tmp

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (E):  + P^0 * ε_{alpha gamma delta} * < mu | r_gamma [V, r_delta] r_beta | nu > * (nu == lambda)
         !    This is in vcd_env%matrix_rxvr_r
         DO alpha = 1, 3
           CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%matrix_rxvr_r(alpha, vcd_env%dcdr_env%beta)%matrix)
            CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                     sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)

            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
            CALL cp_fm_trace(mo_coeff, tmp_aomo, aat_tmp)
            aat_tmp = aat_prefactor*aat_tmp

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (F):  - P^0 * ε_{alpha gamma delta} * < mu | r_gamma [[V, r_beta], r_delta] | nu > * (eta == lambda)
         !        + P^0 * ε_{alpha gamma delta} * < mu | [[V, r_beta], r_delta] | nu > * (eta == lambda) * R_gamma
         !    The negative is in vcd_env%matrix_r_doublecom
         DO alpha = 1, 3
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_r_doublecom(alpha, vcd_env%dcdr_env%beta)%matrix, &
                                         mo_coeff, tmp_aomo, ncol=nmo)
            CALL cp_fm_trace(mo_coeff, tmp_aomo, aat_tmp)
            aat_tmp = -aat_prefactor*aat_tmp

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (G):   - iP^1 * ε_{alpha gamma delta} * < mu | r_gamma [V, r_delta] | nu >
         DO alpha = 1, 3
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_rxrv(alpha)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
            CALL cp_fm_trace(tmp_aomo, vcd_env%dCV_prime(ispin), aat_tmp)

            !  I can take the positive, because build_com_mom_nl computes r x [r, V]
            aat_tmp = 2._dp*aat_prefactor*aat_tmp

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
                 + aat_tmp
         END DO

         ! All the reference dependent stuff
         ! (C) iP^1 * ε_{alpha gamma delta} * < mu | ∂_delta | nu > * (- R_gamma)
         DO alpha = 1, 3
            aat_tmp = 0._dp

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE
                  ! dipvel_ao = + < a | ∂ | b >
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dipvel_ao(delta)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(tmp_aomo, vcd_env%dCV_prime(ispin), tmp_trace)

                  ! The negative sign is due to (r - O^mag_gamma) and otherwise this is
                  !   exactly the APT dipvel(beta, delta) * (-O^mag_gamma)
                  aat_tmp = aat_tmp + 2._dp*aat_prefactor*tmp_trace*lc_tmp*(-vcd_env%magnetic_origin_atom(gamma))
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (G):  - iP^1 * ε_{alpha gamma delta} * < mu | [V, r_delta] | nu > * (- R_gamma)
         DO alpha = 1, 3
            aat_tmp = 0._dp
            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE
                  ! hcom = < a | [r, V] | b > = - < a | [V, r] | b >
                  ! mo_coeff * dCV_prime = + iP1
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%hcom(delta)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(tmp_aomo, vcd_env%dCV_prime(ispin), tmp_trace)

                  ! This is exactly APT hcom(beta, delta)
                  aat_tmp = aat_tmp + 2._dp*aat_prefactor*tmp_trace*lc_tmp*(-vcd_env%magnetic_origin_atom(gamma))
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  mag_vel, vel, mag
         ! Ai)   + ε_{alpha gamma delta} * R_beta R_gamma * < mu | ∂_delta | nu > * (mu - nu)
         ! Aii)  + ε_{alpha gamma delta} * (-R_beta) * < mu | r_gamma ∂_delta | nu > * (mu - nu)
         ! Aiii) + ε_{alpha gamma delta} * (-R_gamma) * < mu | r_beta ∂_delta | nu > * (mu - nu)
         DO alpha = 1, 3
            aat_tmp = 0._dp
            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE
                  ! iii) - R_gamma * < mu | r_beta ∂_delta | nu > * (mu - nu)
                  ! mag
                  ! matrix_difdip2(beta, alpha) = - < a | r_beta | ∂_alpha b >  * (mu - nu)
                  !   so I need matrix_difdip2(beta, delta)
                  ! Only this part correspond to the APT difdip(beta, alpha)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_difdip2(vcd_env%dcdr_env%beta, delta)%matrix, mo_coeff, &
                                               tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)

                  ! There is a negative sign here, because the routine dipole_velocity_deriv calculates
                  !   the derivatives with respect to nuclear positions and we need electronic positions
                  aat_tmp = aat_tmp - lc_tmp*aat_prefactor*tmp_trace*(-vcd_env%magnetic_origin_atom(gamma))

                  ! This part doesn't appear in the APT
                  ! ii)  - R_beta * < mu | r_gamma ∂_delta | nu > * (mu - nu)
                  ! vel
                  ! matrix_difdip2(beta, alpha) = - < a | r_beta | ∂_alpha b > * (mu - nu)
                  !   so I need matrix_difdip2(gamma, delta)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_difdip2(gamma, delta)%matrix, mo_coeff, &
                                               tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)

                  ! There is a negative sign here, because the routine dipole_velocity_deriv calculates
                  !   the derivatives with respect to nuclear positions and we need electronic positions
                  aat_tmp = aat_tmp - lc_tmp*aat_prefactor*tmp_trace*(-vcd_env%spatial_origin_atom(vcd_env%dcdr_env%beta))

                  ! i)   + R_beta R_gamma * < mu | ∂_delta | nu > * (mu - nu)
                  ! mag_vel
                  ! dipvel_ao = + < a | ∂ | b >
                  CALL dbcsr_desymmetrize(vcd_env%dipvel_ao(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
                  CALL dbcsr_desymmetrize(vcd_env%dipvel_ao(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.FALSE., lambda=vcd_env%dcdr_env%lambda)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)
                  CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                                 1._dp, -1._dp)

                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
                  aat_tmp = aat_tmp + lc_tmp*aat_prefactor*tmp_trace* &
                            (vcd_env%spatial_origin_atom(vcd_env%dcdr_env%beta)*vcd_env%magnetic_origin_atom(gamma))

               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         ! (B):  P^0 * ε_{alpha gamma beta} * < mu | nu > * (nu == lambda) * R_gamma
         DO alpha = 1, 3
            aat_tmp = 0._dp

            DO gamma = 1, 3
               lc_tmp = Levi_Civita(alpha, gamma, vcd_env%dcdr_env%beta)
               IF (lc_tmp == 0._dp) CYCLE
               CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(vcd_env%dcdr_env%beta)%matrix, 0.0_dp)
               CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s1(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
               CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", sab_all, &
                                        vcd_env%dcdr_env%lambda, direction_Or=.TRUE.)
               CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
               CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)

               ! This is in total positive because we are calculating
               !  -1/2c * P * < a | b > * (delta == beta) * (nu == lambda) * (-R_gamma)
               ! The whole term corresponds to difdip_s
               aat_tmp = aat_tmp + lc_tmp*aat_prefactor*tmp_trace*vcd_env%magnetic_origin_atom(gamma)
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (D):  - P^0 * ε_{alpha gamma delta} * < mu | r_gamma r_beta [V, r_delta] | nu > * (mu == lambda)
         !  mag, vel, mag_vel
         ! Di)   - ε_{alpha gamma delta} * (-R_gamma) * < mu | r_beta [V, r_delta] | nu > * (mu == lambda)
         ! Dii)  - ε_{alpha gamma delta} * (-R_beta) * < mu | r_gamma [V, r_delta] | nu > * (mu == lambda)
         ! Diii) - ε_{alpha gamma delta} * R_beta R_gamma * < mu | [V, r_delta] | nu > * (mu == lambda)

         DO alpha = 1, 3
            aat_tmp = 0._dp
            CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 0._dp)

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE
                  ! vcd_env%matrix_rrcom(alpha, beta) = r_beta * [V, r_alpha]

                  ! This corresponds to rcom
                  ! Di) mag
                  ! -(-R_gamma) * < mu | r_beta [V, r_delta] | nu > * (mu == lambda)
                  ! vcd_env%matrix_rrcom(alpha, beta) = r_beta * [V, r_alpha]
                  !  so I need vcd_env%matrix_rrcom(delta, beta)
                  !  The multiplication with delta was not done for all directions
                  CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                  vcd_env%matrix_rrcom(delta, vcd_env%dcdr_env%beta)%matrix)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.FALSE., lambda=vcd_env%dcdr_env%lambda)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
                  ! The sign is positive in total, because we have the negative coordinate and the whole term was negative
                  aat_tmp = aat_tmp + aat_prefactor*tmp_trace*lc_tmp*vcd_env%magnetic_origin_atom(gamma)

                  ! This doesn't appear in the APT formula
                  ! Dii) vel
                  ! -(-R_beta) * < mu | r_gamma [V, r_delta] | nu > * (mu == lambda)
                  ! vcd_env%matrix_rrcom(alpha, beta) = r_beta * [V, r_alpha]
                  !  so I need vcd_env%matrix_rrcom(delta, gamma)
                  !  The multiplication with delta was already done in SUBROUTINE apt_dV
                  CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%matrix_rrcom(delta, gamma)%matrix)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.FALSE., lambda=vcd_env%dcdr_env%lambda)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
                  aat_tmp = aat_tmp + aat_prefactor*tmp_trace*lc_tmp*vcd_env%spatial_origin_atom(vcd_env%dcdr_env%beta)

                  ! Diii) mag_vel
                  !  - R_beta R_gamma * < mu | [V, r_delta] | nu >
                  ! hcom(delta) = - [V, r_delta]
                  CALL dbcsr_desymmetrize(vcd_env%hcom(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.FALSE., lambda=vcd_env%dcdr_env%lambda)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
                  ! No need for a negative sign, because hcom already contains the negative sign.
                  aat_tmp = aat_tmp + &
                            aat_prefactor*tmp_trace*lc_tmp &
                            *(vcd_env%spatial_origin_atom(vcd_env%dcdr_env%beta)*vcd_env%magnetic_origin_atom(gamma))
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (E):  + P^0 * ε_{alpha gamma delta} * < mu | r_gamma [V, r_delta] r_beta | nu > * (nu == lambda)
         !  mag, vel, mag_vel
         ! Ei)   + ε_{alpha gamma delta} * (-R_gamma) * < mu | [V, r_delta] r_beta | nu > * (nu == lambda)
         ! Eii)  + ε_{alpha gamma delta} * (-R_beta) * < mu | r_gamma [V, r_delta] | nu > * (nu == lambda)
         ! Eiii) + ε_{alpha gamma delta} * R_beta R_gamma * < mu | [V, r_delta] | nu > * (nu == lambda)
         DO alpha = 1, 3
            aat_tmp = 0._dp
            CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, 0._dp)

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE
                  ! vcd_env%matrix_rrcom(alpha, beta) = r_beta * [V, r_alpha]
                  ! vcd_env%matrix_rcomr(alpha, beta) = [V, r_alpha] * r_beta

                  ! This corresponds to rcom
                  ! Ei) mag
                  ! (-R_gamma) * < mu | [V, r_delta] r_beta | nu > * (nu == lambda)
                  ! vcd_env%matrix_rcomr(alpha, beta) = [V, r_alpha] * r_beta
                  !  so I need vcd_env%matrix_rcomr(delta, beta)
                  CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                                  vcd_env%matrix_rcomr(delta, vcd_env%dcdr_env%beta)%matrix)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)

                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
                  aat_tmp = aat_tmp + aat_prefactor*tmp_trace*lc_tmp*(-vcd_env%magnetic_origin_atom(gamma))

                  ! This doesn't appear in the APT formula
                  ! E2) vel
                  ! (-R_beta) * < mu | r_gamma [V, r_delta] | nu > * (nu == lambda)
                  ! vcd_env%matrix_rrcom(alpha, beta) = r_beta * [V, r_alpha]
                  !  so I need vcd_env%matrix_rrcom(delta, gamma)
                  CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%matrix_rrcom(delta, gamma)%matrix)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)

                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
                  aat_tmp = aat_tmp + aat_prefactor*tmp_trace*lc_tmp*(-vcd_env%spatial_origin_atom(vcd_env%dcdr_env%beta))

                  ! E3) mag_vel
                  ! R_beta R_gamma * < mu | [V, r_delta] | nu > * (nu == lambda)
                  CALL dbcsr_desymmetrize(vcd_env%hcom(delta)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix)
                  CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                           sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)

                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
                  ! There has to be a minus here, because hcom = [r, V] = - [V, r]
                  aat_tmp = aat_tmp - &
                            aat_prefactor*tmp_trace*lc_tmp* &
                            (vcd_env%spatial_origin_atom(vcd_env%dcdr_env%beta)*vcd_env%magnetic_origin_atom(gamma))
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
         END DO

         !  (F):  - P^0 * ε_{alpha gamma delta} * < mu | [[V, r_beta], r_delta] | nu > * (eta == lambda) * (-R_gamma)
         ! This corresponds to APT dcom
         DO alpha = 1, 3
            aat_tmp = 0._dp

            DO gamma = 1, 3
               DO delta = 1, 3
                  lc_tmp = Levi_Civita(alpha, gamma, delta)
                  IF (lc_tmp == 0._dp) CYCLE
                  ! vcd_env%matrix_dcom(alpha, vcd_env%dcdr_env%beta) = - < mu | [ [V, r_beta], r_alpha ] | nu >
                  !  so I need matrix_dcom(delta, vcd_env%dcdr_env%beta)
                  CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_dcom(delta, vcd_env%dcdr_env%beta)%matrix, &
                                               mo_coeff, tmp_aomo, ncol=nmo)
                  CALL cp_fm_trace(mo_coeff, tmp_aomo, tmp_trace)
                  ! matrix_dcom has the negative sign and we include the negative sign of the coordinate
                  aat_tmp = aat_tmp + aat_prefactor*tmp_trace*lc_tmp*(-vcd_env%magnetic_origin_atom(gamma))
               END DO
            END DO

            aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
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
                     = aat_atom(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + aat_tmp
               END DO
            END DO
         END IF
      END ASSOCIATE

      CALL cp_fm_release(tmp_aomo)
      CALL timestop(handle)
   END SUBROUTINE aat_dV

! **************************************************************************************************
!> \brief Compute E_{alpha beta}^lambda = d/dV^lambda_beta <\mu_alpha> = d/dV^lambda_beta < \dot{r} >
!>        The directions alpha, beta are stored in vcd_env%dcdr_env
!> \param vcd_env ...
!> \param qs_env ...
!> \author Edward Ditler, Tomas Zimmermann
! **************************************************************************************************
   SUBROUTINE apt_dV(vcd_env, qs_env)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apt_dV'
      INTEGER, PARAMETER                                 :: ispin = 1
      REAL(dp), PARAMETER                                :: f_spin = 2._dp

      INTEGER                                            :: alpha, handle, ikind, nao, nmo
      LOGICAL                                            :: ghost
      REAL(dp)                                           :: charge
      REAL(KIND=dp)                                      :: apt_dcom, apt_difdip, apt_dipvel, &
                                                            apt_hcom, apt_rcom
      TYPE(cp_fm_type)                                   :: buf, matrix_dSdV_mo
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      sab_all=sab_all, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      nmo = vcd_env%dcdr_env%nmo(ispin)
      nao = vcd_env%dcdr_env%nao

      ASSOCIATE (apt_el => vcd_env%apt_el_nvpt, &
                 apt_nuc => vcd_env%apt_nuc_nvpt, &
                 apt_total => vcd_env%apt_total_nvpt, &
                 mo_coeff => vcd_env%dcdr_env%mo_coeff(ispin), &
                 deltaR => vcd_env%dcdr_env%deltaR)

         ! build the full matrices
         CALL cp_fm_create(buf, vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)
         CALL cp_fm_create(matrix_dSdV_mo, vcd_env%dcdr_env%momo_fm_struct(ispin)%struct)

         ! STEP 1: dCV contribution (dipvel + commutator)
         ! <mu|∂_alpha|nu> and <mu|[r_alpha, V]|nu> in AO basis
         ! We compute tr(c_1^* x ∂_munu x c_0) + tr(c_0 x ∂_munu x c_1)
         ! We compute tr(c_1^* x [,]_munu x c_0) + tr(c_0 x [,]_munu x c_1)
         CALL cp_fm_scale_and_add(0._dp, vcd_env%dCV_prime(ispin), -1._dp, vcd_env%dCV(ispin))

         ! Ref independent
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_dSdV(vcd_env%dcdr_env%beta)%matrix, mo_coeff, &
                                      buf, ncol=nmo)
         CALL parallel_gemm("T", "N", nmo, nmo, nao, &
                            1.0_dp, mo_coeff, buf, &
                            0.0_dp, matrix_dSdV_mo)

         CALL parallel_gemm("N", "N", nao, nmo, nmo, &
                            -0.5_dp, mo_coeff, matrix_dSdV_mo, &
                            1.0_dp, vcd_env%dCV_prime(ispin))

         ! + i∂ - i[Vnl, r]
         DO alpha = 1, 3
            CALL cp_fm_set_all(buf, 0.0_dp)
            apt_dipvel = 0.0_dp

            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dipvel_ao(alpha)%matrix, mo_coeff, buf, ncol=nmo)
            CALL cp_fm_trace(buf, vcd_env%dCV_prime(ispin), apt_dipvel)
            ! dipvel_ao = + < a | ∂ | b >
            ! mo_coeff * dCV_prime = + iP1
            apt_dipvel = 2._dp*apt_dipvel
            apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + f_spin*apt_dipvel
         END DO

         DO alpha = 1, 3
            CALL cp_fm_set_all(buf, 0.0_dp)
            apt_hcom = 0.0_dp
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%hcom(alpha)%matrix, mo_coeff, buf, ncol=nmo)
            CALL cp_fm_trace(buf, vcd_env%dCV_prime(ispin), apt_hcom)

            ! hcom = < a | [r, V] | b > = - < a | [V, r] | b >
            ! mo_coeff * dCV_prime = + iP1
            apt_hcom = +2._dp*apt_hcom

            apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + f_spin*apt_hcom
         END DO  !x/y/z

         ! STEP 2: basis function derivative contribution
      !! difdip_s
         CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(vcd_env%dcdr_env%beta)%matrix, 0.0_dp)
         CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s1(1)%matrix, &
                                 vcd_env%dcdr_env%matrix_nosym_temp(vcd_env%dcdr_env%beta)%matrix)
         CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(vcd_env%dcdr_env%beta)%matrix, qs_kind_set, "ORB", sab_all, &
                                  vcd_env%dcdr_env%lambda, direction_Or=.TRUE.)

         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(vcd_env%dcdr_env%beta)%matrix, mo_coeff, &
                                      buf, ncol=nmo, alpha=1._dp, beta=0._dp)
         CALL cp_fm_trace(mo_coeff, buf, apt_difdip)

         apt_difdip = -f_spin*apt_difdip
         apt_el(vcd_env%dcdr_env%beta, vcd_env%dcdr_env%beta, vcd_env%dcdr_env%lambda) &
            = apt_el(vcd_env%dcdr_env%beta, vcd_env%dcdr_env%beta, vcd_env%dcdr_env%lambda) + apt_difdip

      !! difdip(j, idir) = < a | r_j | ∂_idir b >
      !! matrix_difdip2(beta, alpha) = < a | r_beta | ∂_alpha b >
         DO alpha = 1, 3 ! x/y/z for differentiated AO
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_difdip2(vcd_env%dcdr_env%beta, alpha)%matrix, mo_coeff, &
                                         buf, ncol=nmo, alpha=1._dp, beta=0._dp)

            CALL cp_fm_trace(mo_coeff, buf, apt_difdip)
            ! There is a negative sign here, because the routine dipole_velocity_deriv calculates
            !   the derivatives with respect to nuclear positions and we need electronic positions
            apt_difdip = -f_spin*apt_difdip
            apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + apt_difdip

         END DO !alpha

         ! STEP 3: The terms r * [V, r]
         ! vcd_env%matrix_rrcom(alpha, beta) = r_beta * [V, r_alpha]
         ! vcd_env%matrix_rcomr(alpha, beta) = [V, r_alpha] * r_beta
         DO alpha = 1, 3 ! x/y/z for differentiated AO
            CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%matrix_rcomr(alpha, vcd_env%dcdr_env%beta)%matrix)
            CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, vcd_env%matrix_rrcom(alpha, vcd_env%dcdr_env%beta)%matrix)

            CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                                     sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)
            CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, qs_kind_set, "ORB", &
                                     sab_all, direction_Or=.FALSE., lambda=vcd_env%dcdr_env%lambda)

            CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                           1.0_dp, -1.0_dp)

            CALL cp_fm_set_all(buf, 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, buf, ncol=nmo)
            CALL cp_fm_trace(mo_coeff, buf, apt_rcom)

            apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + f_spin*apt_rcom
         END DO !alpha

         ! STEP 4: pseudopotential derivative contribution
         ! vcd_env%matrix_dcom(alpha, vcd_env%dcdr_env%beta) = - < mu | [ [V, r_beta], r_alpha ] | nu >
         DO alpha = 1, 3 !x/y/z for differentiated AO
            CALL cp_fm_set_all(buf, 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_dcom(alpha, vcd_env%dcdr_env%beta)%matrix, mo_coeff, buf, ncol=nmo)
            CALL cp_fm_trace(mo_coeff, buf, apt_dcom)
            apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + f_spin*apt_dcom
         END DO !alpha

         ! The reference point dependent terms:
      !! difdip_munu
         ! The additional term here is < a | db/dr(alpha)> * (delta_a - delta_b) * ref_point(beta)
         ! in qs_env%matrix_s1(2:4) there is < da/dR | b > = - < da/dr | b > = < a | db/dr >
         DO alpha = 1, 3
            CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix, 0._dp)
            CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp2(alpha)%matrix, 0._dp)
            CALL dbcsr_desymmetrize(vcd_env%dcdr_env%matrix_s(alpha + 1)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix)
            CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp2(alpha)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix)

            ! < a | db/dr(alpha) > * R^lambda_beta * delta^lambda_nu
            CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix, qs_kind_set, "ORB", sab_all, &
                                     vcd_env%dcdr_env%lambda, direction_Or=.TRUE.)
            ! < a | db/dr(alpha) > * R^lambda_beta * delta^lambda_mu
            CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp2(alpha)%matrix, qs_kind_set, "ORB", sab_all, &
                                     vcd_env%dcdr_env%lambda, direction_Or=.FALSE.)

            ! < a | db/dr > * R^lambda_beta * ( delta^lambda_mu - delta^lambda_nu )
            CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp2(alpha)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix, &
                           1._dp, -1._dp)

            CALL cp_fm_set_all(buf, 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp2(alpha)%matrix, mo_coeff, buf, ncol=nmo)
            CALL cp_fm_trace(mo_coeff, buf, apt_difdip)

            ! And the whole contribution is
            ! - < a | db/dr > * (mu - nu) * ref_point
            apt_difdip = -apt_difdip*vcd_env%spatial_origin_atom(vcd_env%dcdr_env%beta)

            apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + f_spin*apt_difdip
         END DO

         ! And the additional factor to rcom
         ! < mu | [V, r] | nu > * R^lambda_beta * delta^lambda_mu
         ! - < mu | [V, r] | nu > * R^lambda_beta * delta^lambda_nu
         !
         !                  vcd_env%hcom(alpha) = - < mu | [V, r_alpha] | nu >
         ! particle_set(lambda)%r(vcd_env%dcdr_env%beta) = R^lambda_beta

         DO alpha = 1, 3
            CALL dbcsr_set(vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix, 0._dp)
            CALL dbcsr_desymmetrize(vcd_env%hcom(alpha)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix)
            CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp2(alpha)%matrix, vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix)

            ! < mu | [V, r] | nu > * delta^lambda_nu
            CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix, qs_kind_set, "ORB", sab_all, &
                                     vcd_env%dcdr_env%lambda, direction_Or=.TRUE.)
            ! < mu | [V, r] | nu > * delta^lambda_mu
            CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp2(alpha)%matrix, qs_kind_set, "ORB", sab_all, &
                                     vcd_env%dcdr_env%lambda, direction_Or=.FALSE.)

            CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix, &
                           vcd_env%dcdr_env%matrix_nosym_temp2(alpha)%matrix, -1._dp, +1._dp)

            CALL cp_fm_set_all(buf, 0.0_dp)
            CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(alpha)%matrix, mo_coeff, buf, ncol=nmo)
            CALL cp_fm_trace(mo_coeff, buf, apt_rcom)
            apt_rcom = -vcd_env%spatial_origin_atom(vcd_env%dcdr_env%beta)*apt_rcom

            apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) &
               = apt_el(vcd_env%dcdr_env%beta, alpha, vcd_env%dcdr_env%lambda) + f_spin*apt_rcom
         END DO

         ! STEP 5: nuclear contribution
         ASSOCIATE (atomic_kind => particle_set(vcd_env%dcdr_env%lambda)%atomic_kind)
            CALL get_atomic_kind(atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost)
            IF (.NOT. ghost) THEN
               apt_nuc(vcd_env%dcdr_env%beta, vcd_env%dcdr_env%beta, vcd_env%dcdr_env%lambda) = &
                  apt_nuc(vcd_env%dcdr_env%beta, vcd_env%dcdr_env%beta, vcd_env%dcdr_env%lambda) + charge
            END IF
         END ASSOCIATE

         ! STEP 6: deallocations
         CALL cp_fm_release(buf)
         CALL cp_fm_release(matrix_dSdV_mo)

      END ASSOCIATE

      CALL timestop(handle)
   END SUBROUTINE apt_dV

! **************************************************************************************************
!> \brief Initialize the matrices for the NVPT calculation
!> \param vcd_env ...
!> \param qs_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE prepare_per_atom_vcd(vcd_env, qs_env)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_per_atom_vcd'

      INTEGER                                            :: handle, i, ispin, j
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
                      sab_orb=sab_orb, sab_all=sab_all, sap_ppnl=sap_ppnl, &
                      qs_kind_set=qs_kind_set, particle_set=particle_set, cell=cell)

      IF (vcd_env%distributed_origin) THEN
         vcd_env%magnetic_origin_atom(:) = particle_set(vcd_env%dcdr_env%lambda)%r(:) - vcd_env%magnetic_origin(:)
         vcd_env%spatial_origin_atom = particle_set(vcd_env%dcdr_env%lambda)%r(:) - vcd_env%spatial_origin(:)
      END IF

      ! Reset the matrices
      DO ispin = 1, dft_control%nspins
         DO j = 1, 3
            CALL dbcsr_set(vcd_env%matrix_dSdV(j)%matrix, 0._dp)
            CALL dbcsr_set(vcd_env%matrix_drpnl(j)%matrix, 0._dp)

            DO i = 1, 3
               CALL dbcsr_set(vcd_env%matrix_dcom(i, j)%matrix, 0.0_dp)
               CALL dbcsr_set(vcd_env%matrix_difdip2(i, j)%matrix, 0._dp)
            END DO
         END DO
         CALL cp_fm_set_all(vcd_env%op_dV(ispin), 0._dp)
         CALL dbcsr_set(vcd_env%matrix_hxc_dsdv(ispin)%matrix, 0._dp)
      END DO

      ! operator dV
      ! <mu|d/dV_beta [V, r_alpha]|nu>
      CALL build_dcom_rpnl(vcd_env%matrix_dcom, qs_kind_set, sab_orb, sap_ppnl, &
                           dft_control%qs_control%eps_ppnl, particle_set, vcd_env%dcdr_env%lambda)

      ! PP derivative
      CALL build_drpnl_matrix(vcd_env%matrix_drpnl, qs_kind_set, sab_all, sap_ppnl, &
                              dft_control%qs_control%eps_ppnl, particle_set, pseudoatom=vcd_env%dcdr_env%lambda)
      ! lin_mom
      DO i = 1, 3
         CALL dbcsr_set(vcd_env%dipvel_ao_delta(i)%matrix, 0._dp)
         CALL dbcsr_copy(vcd_env%dipvel_ao_delta(i)%matrix, vcd_env%dipvel_ao(i)%matrix)
      END DO

      CALL hr_mult_by_delta_3d(vcd_env%dipvel_ao_delta, qs_kind_set, "ORB", &
                               sab_all, vcd_env%dcdr_env%delta_basis_function, direction_Or=.TRUE.)

      ! dS/dV
      CALL build_dSdV_matrix(qs_env, vcd_env%matrix_dSdV, &
                             deltaR=vcd_env%dcdr_env%delta_basis_function, &
                             rcc=vcd_env%spatial_origin_atom)

      CALL dipole_velocity_deriv(qs_env, vcd_env%matrix_difdip2, 1, lambda=vcd_env%dcdr_env%lambda, &
                                 rc=[0._dp, 0._dp, 0._dp])
      ! AAT
      ! moments_throw: x, y, z, xx, xy, xz, yy, yz, zz
      ! moments_der:  (moment, xyz derivative)
      ! build_local_moments_der_matrix uses adbdr for calculating derivatives of the *primitive*
      !  on the right. So the resulting
      !  moments_der(moment, delta) = - < a | moment \partial_\delta | b >
      DO i = 1, 9 ! x, y, z, xx, xy, xz, yy, yz, zz
         DO j = 1, 3
            CALL dbcsr_set(vcd_env%moments_der_right(i, j)%matrix, 0.0_dp)
            CALL dbcsr_set(vcd_env%moments_der_left(i, j)%matrix, 0.0_dp)
         END DO
      END DO

      DO i = 1, 9
         DO j = 1, 3 ! derivatives
            CALL dbcsr_desymmetrize(vcd_env%moments_der(i, j)%matrix, vcd_env%moments_der_right(i, j)%matrix) ! A2
            CALL dbcsr_desymmetrize(vcd_env%moments_der(i, j)%matrix, vcd_env%moments_der_left(i, j)%matrix)  ! A1

            !  - < mu | r_beta r_gamma ∂_delta | nu > * (mu/nu == lambda)
            CALL hr_mult_by_delta_1d(vcd_env%moments_der_right(i, j)%matrix, qs_kind_set, "ORB", &
                                     sab_all, direction_Or=.TRUE., lambda=vcd_env%dcdr_env%lambda)
            CALL hr_mult_by_delta_1d(vcd_env%moments_der_left(i, j)%matrix, qs_kind_set, "ORB", &
                                     sab_all, direction_Or=.FALSE., lambda=vcd_env%dcdr_env%lambda)
         END DO
      END DO

      DO i = 1, 3
         DO j = 1, 3
            CALL dbcsr_set(vcd_env%matrix_r_doublecom(i, j)%matrix, 0._dp)
         END DO
      END DO

      CALL build_com_mom_nl(qs_kind_set, sab_all, sap_ppnl, dft_control%qs_control%eps_ppnl, &
                            particle_set, ref_point=[0._dp, 0._dp, 0._dp], cell=cell, &
                            matrix_r_doublecom=vcd_env%matrix_r_doublecom, &
                            pseudoatom=vcd_env%dcdr_env%lambda)

      CALL timestop(handle)

   END SUBROUTINE prepare_per_atom_vcd

! **************************************************************************************************
!> \brief What we are building here is the operator for the NVPT response:
!>     H0 * C1 - S0 * E0 * C1  = - op_dV
!>     linres_solver           = - [ H1 * C0 - S1 * C0 * E0 ]
!>   with
!>     H1 * C0 =   dH/dV * C0
!>               + i[∂]δ * C0
!>               - i S0 * C^(1,R)
!>               + i S0 * C0 * (C0 * S^(1,R) * C0)
!>               - S1 * C0 * E0
!>
!>     H1 * C0 = + i (Hr - rH) * C0                    [STEP 1]
!>               + i[∂]δ * C0                          [STEP 2]
!>               - i[V, r]δ * C0                       [STEP 3]
!>               - i S0 * C^(1,R)                      [STEP 4]
!>               - S1 * C0 * E0                        [STEP 5]
!> \param vcd_env ...
!> \param qs_env ...
!> \author Edward Ditler, Tomas Zimmermann
! **************************************************************************************************
   SUBROUTINE vcd_build_op_dV(vcd_env, qs_env)
      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vcd_build_op_dV'
      INTEGER, PARAMETER                                 :: ispin = 1

      INTEGER                                            :: handle, nao, nmo
      TYPE(cp_fm_type)                                   :: buf
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      sab_all=sab_all, &
                      qs_kind_set=qs_kind_set)

      nmo = vcd_env%dcdr_env%nmo(1)
      nao = vcd_env%dcdr_env%nao

      CALL build_matrix_hr_rh(vcd_env, qs_env, vcd_env%spatial_origin_atom)

      ! STEP 1: hr-rh
      CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, vcd_env%matrix_hr(ispin, vcd_env%dcdr_env%beta)%matrix)
      CALL dbcsr_copy(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, vcd_env%matrix_rh(ispin, vcd_env%dcdr_env%beta)%matrix)

      CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, qs_kind_set, "ORB", &
                               sab_all, vcd_env%dcdr_env%lambda, direction_or=.TRUE.)
      CALL hr_mult_by_delta_1d(vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, qs_kind_set, "ORB", &
                               sab_all, vcd_env%dcdr_env%lambda, direction_or=.FALSE.)
      CALL dbcsr_add(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, &
                     vcd_env%dcdr_env%matrix_nosym_temp(2)%matrix, &
                     1.0_dp, -1.0_dp)

      ASSOCIATE (mo_coeff => vcd_env%dcdr_env%mo_coeff(ispin))
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_nosym_temp(1)%matrix, mo_coeff, &
                                      vcd_env%op_dV(ispin), ncol=nmo, alpha=1.0_dp, beta=0.0_dp)

         ! STEP 2: electronic momentum operator contribution
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dipvel_ao_delta(vcd_env%dcdr_env%beta)%matrix, mo_coeff, &
                                      vcd_env%op_dV(ispin), &
                                      ncol=nmo, alpha=1.0_dp, beta=1.0_dp)

         ! STEP 3: +dV_ppnl/dV, but build_drpnl_matrix gives the negative of dV_ppnl
         ! The arguments (-1, 1) are swapped wrt to the hr-rh term, implying that
         ! direction_Or and direction_hr do what they should.
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_drpnl(vcd_env%dcdr_env%beta)%matrix, mo_coeff, &
                                      vcd_env%op_dV(ispin), &
                                      ncol=nmo, alpha=-1.0_dp, beta=1.0_dp)

         ! STEP 4: - S0 * C^(1,R)
         CALL cp_fm_create(buf, vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%dcdr_env%matrix_s1(1)%matrix, vcd_env%dcdr_env%dCR_prime(ispin), &
                                      vcd_env%op_dV(1), ncol=nmo, alpha=-1.0_dp, beta=1.0_dp)

         ! STEP 5: -S(1,V) * C0 * E0
         CALL cp_dbcsr_sm_fm_multiply(vcd_env%matrix_dSdV(vcd_env%dcdr_env%beta)%matrix, mo_coeff, &
                                      buf, nmo, alpha=1.0_dp, beta=0.0_dp)
         CALL parallel_gemm('N', 'N', nao, nmo, nmo, &
                            -1.0_dp, buf, vcd_env%dcdr_env%chc(ispin), &
                            1.0_dp, vcd_env%op_dV(ispin))

         CALL cp_fm_release(buf)
      END ASSOCIATE

      ! We have built op_dV but plug -op_dV into the linres_solver
      CALL cp_fm_scale(-1.0_dp, vcd_env%op_dV(1))

      ! Revert the matrices
      CALL build_matrix_hr_rh(vcd_env, qs_env, [0._dp, 0._dp, 0._dp])

      CALL timestop(handle)
   END SUBROUTINE vcd_build_op_dV

! *****************************************************************************
!> \brief Get the dC/dV using the vcd_env%op_dV
!> \param vcd_env ...
!> \param p_env ...
!> \param qs_env ...
!> \author Edward Ditler, Tomas Zimmermann
! **************************************************************************************************
   SUBROUTINE vcd_response_dV(vcd_env, p_env, qs_env)

      TYPE(vcd_env_type)                                 :: vcd_env
      TYPE(qs_p_env_type)                                :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vcd_response_dV'
      INTEGER, PARAMETER                                 :: ispin = 1

      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: failure, should_stop
      TYPE(cp_fm_type), DIMENSION(1)                     :: h1_psi0, psi1
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(section_vals_type), POINTER                   :: lr_section, vcd_section

      CALL timeset(routineN, handle)
      failure = .FALSE.

      NULLIFY (linres_control, lr_section, logger)

      CALL get_qs_env(qs_env=qs_env, &
                      linres_control=linres_control, &
                      mos=mos)

      logger => cp_get_default_logger()
      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      vcd_section => section_vals_get_subs_vals(qs_env%input, &
                                                "PROPERTIES%LINRES%vcd")

      output_unit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".linresLog")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T10,A,/)") &
            "*** Self consistent optimization of the response wavefunction ***"
      END IF

      ASSOCIATE (psi0_order => vcd_env%dcdr_env%mo_coeff)
         CALL cp_fm_create(psi1(ispin), vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)
         CALL cp_fm_create(h1_psi0(ispin), vcd_env%dcdr_env%likemos_fm_struct(ispin)%struct)

         ! Restart
         IF (linres_control%linres_restart) THEN
            CALL vcd_read_restart(qs_env, lr_section, psi1, vcd_env%dcdr_env%lambda, vcd_env%dcdr_env%beta, "dCdV")
         ELSE
            CALL cp_fm_set_all(psi1(ispin), 0.0_dp)
         END IF

         IF (output_unit > 0) THEN
            WRITE (output_unit, *) &
               "Response to the perturbation operator referring to the velocity of atom ", &
               vcd_env%dcdr_env%lambda, " in "//ACHAR(vcd_env%dcdr_env%beta + 119)
         END IF

         ! First response to get dCR
         ! (H0-E0) psi1 = (H1-E1) psi0
         ! psi1 = the perturbed wavefunction
         ! h1_psi0 = (H1-E1)
         ! psi0_order = the unperturbed wavefunction
         ! Second response to get dCV
         CALL cp_fm_set_all(vcd_env%dCV(ispin), 0.0_dp)
         CALL cp_fm_set_all(h1_psi0(ispin), 0.0_dp)
         CALL cp_fm_to_fm(vcd_env%op_dV(ispin), h1_psi0(ispin))

         linres_control%lr_triplet = .FALSE. ! we do singlet response
         linres_control%do_kernel = .FALSE. ! no coupled response since imaginary perturbation
         linres_control%converged = .FALSE.
         CALL linres_solver(p_env, qs_env, psi1, h1_psi0, psi0_order, &
                            output_unit, should_stop)
         CALL cp_fm_to_fm(psi1(ispin), vcd_env%dCV(ispin))

         ! Write the new result to the restart file
         IF (linres_control%linres_restart) THEN
            CALL vcd_write_restart(qs_env, lr_section, psi1, vcd_env%dcdr_env%lambda, vcd_env%dcdr_env%beta, "dCdV")
         END IF

      END ASSOCIATE

      ! clean up
      CALL cp_fm_release(psi1(ispin))
      CALL cp_fm_release(h1_psi0(ispin))

      CALL cp_print_key_finished_output(output_unit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      CALL timestop(handle)
   END SUBROUTINE vcd_response_dV

END MODULE qs_vcd
