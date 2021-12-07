!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for kpoint treatment in GW
!> \par History
!>      04.2019 created [Jan Wilhelm]
! **************************************************************************************************
MODULE rpa_gw_kpoints
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_invert,&
                                              cp_cfm_gemm,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_transpose
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_release
   USE input_constants,                 ONLY: kp_weights_W_auto,&
                                              kp_weights_W_tailored,&
                                              kp_weights_W_uniform
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: invmat
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE rpa_gw_im_time_util,             ONLY: compute_weight_re_im,&
                                              get_atom_index_from_basis_function_index
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_kpoints'

   PUBLIC :: compute_Wc_real_space_tau_GW, compute_Wc_kp_tau_GW, &
             compute_wkp_W

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_W_tau ...
!> \param cfm_mat_Q ...
!> \param fm_mat_L_re ...
!> \param fm_mat_L_im ...
!> \param dimen_RI ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param ikp ...
!> \param tj ...
!> \param tau_tj ...
!> \param weights_cos_tf_w_to_t ...
!> \param ikp_local ...
!> \param para_env ...
!> \param kpoints ...
!> \param qs_env ...
!> \param wkp_W ...
!> \param mat_SinvVSinv ...
!> \param do_W_and_not_V ...
! **************************************************************************************************
   SUBROUTINE compute_Wc_real_space_tau_GW(fm_mat_W_tau, cfm_mat_Q, fm_mat_L_re, fm_mat_L_im, &
                                           dimen_RI, num_integ_points, jquad, &
                                           ikp, tj, tau_tj, weights_cos_tf_w_to_t, ikp_local, &
                                           para_env, kpoints, qs_env, wkp_W, mat_SinvVSinv, do_W_and_not_V)

      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_W_tau
      TYPE(cp_cfm_type), INTENT(IN), POINTER             :: cfm_mat_Q
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_mat_L_re, fm_mat_L_im
      INTEGER, INTENT(IN)                                :: dimen_RI, num_integ_points, jquad, ikp
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tj
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weights_cos_tf_w_to_t
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ikp_local
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wkp_W
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_SinvVSinv
      LOGICAL, INTENT(IN)                                :: do_W_and_not_V

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Wc_real_space_tau_GW'

      INTEGER :: handle, handle2, i_global, iatom, iatom_old, iiB, iquad, irow, j_global, jatom, &
         jatom_old, jcol, jjB, jkp, ncol_local, nkp, nrow_local, num_cells
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      LOGICAL                                            :: do_V_and_not_W
      REAL(KIND=dp)                                      :: contribution, omega, tau, weight, &
                                                            weight_im, weight_re
      REAL(KIND=dp), DIMENSION(2)                        :: abs_cutoffs_chi_W
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_L, cfm_mat_work, cfm_mat_work_2
      TYPE(cp_fm_type), POINTER                          :: fm_dummy, fm_mat_work_global, &
                                                            fm_mat_work_local
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL timeset(routineN//"_1", handle2)

      NULLIFY (cfm_mat_work)
      CALL cp_cfm_create(cfm_mat_work, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      NULLIFY (cfm_mat_work_2)
      CALL cp_cfm_create(cfm_mat_work_2, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work_2, z_zero)

      NULLIFY (cfm_mat_L)
      CALL cp_cfm_create(cfm_mat_L, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_L, z_zero)

      ! Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      NULLIFY (fm_mat_work_global)
      CALL cp_fm_create(fm_mat_work_global, fm_mat_W_tau(1)%matrix%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_global, 0.0_dp)

      NULLIFY (fm_mat_work_local)
      CALL cp_fm_create(fm_mat_work_local, cfm_mat_Q%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_local, 0.0_dp)

      CALL timestop(handle2)

      IF (do_W_and_not_V) THEN

         CALL timeset(routineN//"_2", handle2)

         ! calculate [1+Q(iw')]^-1
         CALL cp_cfm_cholesky_invert(cfm_mat_Q)

         ! symmetrize the result
         CALL own_cfm_upper_to_full(cfm_mat_Q, cfm_mat_work)

         ! subtract exchange part by subtracing identity matrix from epsilon
         CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                              nrow_local=nrow_local, &
                              ncol_local=ncol_local, &
                              row_indices=row_indices, &
                              col_indices=col_indices)

         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
                  cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB) - z_one
               END IF
            END DO
         END DO

         CALL timestop(handle2)

         CALL timeset(routineN//"_3.1", handle2)

         ! work = epsilon(iw,k)*L^H(k)
         CALL cp_cfm_gemm('N', 'C', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                          z_zero, cfm_mat_work)

         ! W(iw,k) = L(k)*work
         CALL cp_cfm_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_work, &
                          z_zero, cfm_mat_work_2)

         CALL timestop(handle2)

      ELSE

         ! S^-1(k)V(k)S^-1(k) = L(k)*L(k)^H
         CALL cp_cfm_gemm('N', 'C', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_L, &
                          z_zero, cfm_mat_work_2)

      END IF

      CALL timeset(routineN//"_4", handle2)

      CALL get_kpoint_info(kpoints, xkp=xkp, wkp=wkp, nkp=nkp)
      index_to_cell => kpoints%index_to_cell
      num_cells = SIZE(index_to_cell, 2)

      abs_cutoffs_chi_W = qs_env%mp2_env%ri_rpa_im_time%abs_cutoffs_chi_W
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      ALLOCATE (atom_from_RI_index(dimen_RI))

      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_RI_index, dimen_RI, "RI_AUX")

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat)
      iatom_old = 0
      jatom_old = 0

      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      DO irow = 1, nrow_local
         DO jcol = 1, ncol_local

            iatom = atom_from_RI_index(row_indices(irow))
            jatom = atom_from_RI_index(col_indices(jcol))

            IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN

               CALL compute_weight_re_im(weight_re, weight_im, &
                                         num_cells, iatom, jatom, xkp(1:3, ikp), wkp_W(ikp), &
                                         cell, index_to_cell, hmat, particle_set, abs_cutoffs_chi_W)

               iatom_old = iatom
               jatom_old = jatom

            END IF

            contribution = weight_re*REAL(cfm_mat_work_2%local_data(irow, jcol)) + &
                           weight_im*AIMAG(cfm_mat_work_2%local_data(irow, jcol))

            fm_mat_work_local%local_data(irow, jcol) = fm_mat_work_local%local_data(irow, jcol) + contribution

         END DO
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_5", handle2)

      IF (do_W_and_not_V) THEN

         IF (SUM(ikp_local) > nkp) THEN

            CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)

            DO iquad = 1, num_integ_points

               omega = tj(jquad)
               tau = tau_tj(iquad)
               weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)

               IF (jquad == 1 .AND. ikp == 1) THEN
                  CALL cp_fm_set_all(matrix=fm_mat_W_tau(iquad)%matrix, alpha=0.0_dp)
               END IF

               CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_W_tau(iquad)%matrix, beta=weight, matrix_b=fm_mat_work_global)

            END DO

         ELSE

            DO jkp = 1, nkp

               IF (ANY(ikp_local(:) == jkp)) THEN
                  CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)
               ELSE
                  NULLIFY (fm_dummy)
                  CALL cp_fm_copy_general(fm_dummy, fm_mat_work_global, para_env)
               END IF

               DO iquad = 1, num_integ_points

                  omega = tj(jquad)
                  tau = tau_tj(iquad)
                  weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)

                  IF (jquad == 1 .AND. jkp == 1) THEN
                     CALL cp_fm_set_all(matrix=fm_mat_W_tau(iquad)%matrix, alpha=0.0_dp)
                  END IF

                  CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_W_tau(iquad)%matrix, beta=weight, &
                                           matrix_b=fm_mat_work_global)

               END DO

            END DO

         END IF

      END IF

      do_V_and_not_W = .NOT. do_W_and_not_V
      IF (do_V_and_not_W) THEN

         IF (SUM(ikp_local) > nkp) THEN
            CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)
            CALL fm_mat_work_global_to_mat_SinvVSinv(mat_SinvVSinv, fm_mat_work_global)
         ELSE
            DO jkp = 1, nkp
               IF (ANY(ikp_local(:) == jkp)) THEN
                  CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)
               ELSE
                  NULLIFY (fm_dummy)
                  CALL cp_fm_copy_general(fm_dummy, fm_mat_work_global, para_env)
               END IF
               CALL fm_mat_work_global_to_mat_SinvVSinv(mat_SinvVSinv, fm_mat_work_global)
            END DO
         END IF
      END IF

      CALL cp_cfm_release(cfm_mat_work)
      CALL cp_cfm_release(cfm_mat_work_2)
      CALL cp_cfm_release(cfm_mat_L)
      CALL cp_fm_release(fm_mat_work_global)
      CALL cp_fm_release(fm_mat_work_local)
      DEALLOCATE (atom_from_RI_index)

      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE compute_Wc_real_space_tau_GW

! **************************************************************************************************
!> \brief ...
!> \param mat_SinvVSinv ...
!> \param fm_mat_work_global ...
! **************************************************************************************************
   SUBROUTINE fm_mat_work_global_to_mat_SinvVSinv(mat_SinvVSinv, fm_mat_work_global)

      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_SinvVSinv
      TYPE(cp_fm_type), POINTER                          :: fm_mat_work_global

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_mat_work_global_to_mat_SinvVSinv'

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type)                                 :: mat_work

      CALL timeset(routineN, handle)

      NULLIFY (mat_work%matrix)
      ALLOCATE (mat_work%matrix)
      CALL dbcsr_create(mat_work%matrix, template=mat_SinvVSinv%matrix)

      CALL copy_fm_to_dbcsr(fm_mat_work_global, mat_work%matrix, keep_sparsity=.FALSE.)

      CALL dbcsr_add(mat_SinvVSinv%matrix, mat_work%matrix, 1.0_dp, 1.0_dp)

      CALL dbcsr_release(mat_work%matrix)
      DEALLOCATE (mat_work%matrix)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_W_kp_tau ...
!> \param cfm_mat_Q ...
!> \param fm_mat_L_re ...
!> \param fm_mat_L_im ...
!> \param dimen_RI ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param ikp ...
!> \param tj ...
!> \param tau_tj ...
!> \param weights_cos_tf_w_to_t ...
! **************************************************************************************************
   SUBROUTINE compute_Wc_kp_tau_GW(cfm_mat_W_kp_tau, cfm_mat_Q, fm_mat_L_re, fm_mat_L_im, &
                                   dimen_RI, num_integ_points, jquad, &
                                   ikp, tj, tau_tj, weights_cos_tf_w_to_t)

      TYPE(cp_cfm_p_type), DIMENSION(:, :), INTENT(IN)   :: cfm_mat_W_kp_tau
      TYPE(cp_cfm_type), INTENT(IN), POINTER             :: cfm_mat_Q
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_mat_L_re, fm_mat_L_im
      INTEGER, INTENT(IN)                                :: dimen_RI, num_integ_points, jquad, ikp
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tj
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weights_cos_tf_w_to_t

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Wc_kp_tau_GW'

      INTEGER                                            :: handle, handle2, i_global, iiB, iquad, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: omega, tau, weight
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_L, cfm_mat_work

      CALL timeset(routineN, handle)

      NULLIFY (cfm_mat_work)
      CALL cp_cfm_create(cfm_mat_work, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      NULLIFY (cfm_mat_L)
      CALL cp_cfm_create(cfm_mat_L, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_L, z_zero)

      CALL timeset(routineN//"_cholesky_inv", handle2)

      ! calculate [1+Q(iw')]^-1
      CALL cp_cfm_cholesky_invert(cfm_mat_Q)

      ! symmetrize the result
      CALL own_cfm_upper_to_full(cfm_mat_Q, cfm_mat_work)

      ! subtract exchange part by subtracing identity matrix from epsilon
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB) - z_one
            END IF
         END DO
      END DO

      CALL timestop(handle2)

      ! Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      ! work = epsilon(iw,k)*L^H(k)
      CALL cp_cfm_gemm('N', 'C', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                       z_zero, cfm_mat_work)

      ! W(iw,k) = L(k)*work
      CALL cp_cfm_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_work, &
                       z_zero, cfm_mat_Q)

      DO iquad = 1, num_integ_points
         omega = tj(jquad)
         tau = tau_tj(iquad)
         weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)
         CALL cp_cfm_scale_and_add(alpha=z_one, matrix_a=cfm_mat_W_kp_tau(ikp, iquad)%matrix, &
                                   beta=CMPLX(weight, KIND=dp), matrix_b=cfm_mat_Q)
      END DO

      CALL cp_cfm_release(cfm_mat_work)
      CALL cp_cfm_release(cfm_mat_L)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param wkp_W ...
!> \param wkp_V ...
!> \param kpoints ...
!> \param h_inv ...
!> \param periodic ...
! **************************************************************************************************
   SUBROUTINE compute_wkp_W(qs_env, wkp_W, wkp_V, kpoints, h_inv, periodic)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: wkp_W, wkp_V
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv
      INTEGER, DIMENSION(3)                              :: periodic

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_wkp_W'

      INTEGER                                            :: handle, i_x, ikp, info, j_y, k_z, &
                                                            kpoint_weights_W_method, n_x, n_y, &
                                                            n_z, nkp, nsuperfine, num_lin_eqs
      REAL(KIND=dp)                                      :: exp_kpoints, integral, k_sq, weight
      REAL(KIND=dp), DIMENSION(3)                        :: k_vec, x_vec
      REAL(KIND=dp), DIMENSION(:), POINTER               :: right_side, wkp, wkp_tmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix_lin_eqs, xkp

      CALL timeset(routineN, handle)

      kpoint_weights_W_method = qs_env%mp2_env%ri_rpa_im_time%kpoint_weights_W_method

      CALL get_kpoint_info(kpoints, xkp=xkp, wkp=wkp, nkp=nkp)

      ! we determine the kpoint weights of the Monkhors Pack mesh new
      ! such that the functions 1/k^2, 1/k and const are integrated exactly
      ! in the Brillouin zone
      ! this is done by minimizing sum_i |w_i|^2 where w_i are the weights of
      ! the i-th kpoint under the following constraints:
      ! 1) 1/k^2, 1/k and const are integrated exactly
      ! 2) the kpoint weights of kpoints with identical absolute value are
      !    the same, of e.g. (1/8,3/8,3/8) same weight as for (3/8,1/8,3/8)
      ! for 1d and 2d materials: we use ordinary Monkhorst-Pack weights, checked
      ! by SUM(periodic) == 3
      ALLOCATE (wkp_V(nkp), wkp_W(nkp))

      IF (kpoint_weights_W_method == kp_weights_W_uniform) THEN

         wkp_V(:) = wkp(:)
         wkp_W(:) = wkp(:)

      ELSE IF (kpoint_weights_W_method == kp_weights_W_tailored .OR. &
               kpoint_weights_W_method == kp_weights_W_auto) THEN

         IF (SUM(periodic) == 1 .OR. SUM(periodic) == 2) THEN

            wkp_V(:) = wkp(:)
            wkp_W(:) = wkp(:)

         ELSE IF (SUM(periodic) == 3) THEN

            IF (kpoint_weights_W_method == kp_weights_W_tailored) &
               exp_kpoints = qs_env%mp2_env%ri_rpa_im_time%exp_tailored_weights
            IF (kpoint_weights_W_method == kp_weights_W_auto) exp_kpoints = -2.0_dp

            ! first, compute the integral of f(k)=1/k^2 and 1/k on super fine grid
            nsuperfine = 500
            integral = 0.0_dp

            IF (periodic(1) == 1) THEN
               n_x = nsuperfine
            ELSE
               n_x = 1
            END IF
            IF (periodic(2) == 1) THEN
               n_y = nsuperfine
            ELSE
               n_y = 1
            END IF
            IF (periodic(3) == 1) THEN
               n_z = nsuperfine
            ELSE
               n_z = 1
            END IF

            ! actually, there is the factor *det_3x3(h_inv) missing to account for the
            ! integration volume but for wkp det_3x3(h_inv) is needed
            weight = 1.0_dp/(REAL(n_x, dp)*REAL(n_y, dp)*REAL(n_z, dp))
            DO i_x = 1, n_x
               DO j_y = 1, n_y
                  DO k_z = 1, n_z

                     IF (periodic(1) == 1) THEN
                        x_vec(1) = (REAL(i_x - nsuperfine/2, dp) - 0.5_dp)/REAL(nsuperfine, dp)
                     ELSE
                        x_vec(1) = 0.0_dp
                     END IF
                     IF (periodic(2) == 1) THEN
                        x_vec(2) = (REAL(j_y - nsuperfine/2, dp) - 0.5_dp)/REAL(nsuperfine, dp)
                     ELSE
                        x_vec(2) = 0.0_dp
                     END IF
                     IF (periodic(3) == 1) THEN
                        x_vec(3) = (REAL(k_z - nsuperfine/2, dp) - 0.5_dp)/REAL(nsuperfine, dp)
                     ELSE
                        x_vec(3) = 0.0_dp
                     END IF

                     k_vec = MATMUL(h_inv(1:3, 1:3), x_vec)
                     k_sq = k_vec(1)**2 + k_vec(2)**2 + k_vec(3)**2
                     integral = integral + weight*k_sq**(exp_kpoints*0.5_dp)

                  END DO
               END DO
            END DO

            num_lin_eqs = nkp + 2

            ALLOCATE (matrix_lin_eqs(num_lin_eqs, num_lin_eqs))
            matrix_lin_eqs(:, :) = 0.0_dp

            DO ikp = 1, nkp

               k_vec = MATMUL(h_inv(1:3, 1:3), xkp(1:3, ikp))
               k_sq = k_vec(1)**2 + k_vec(2)**2 + k_vec(3)**2

               matrix_lin_eqs(ikp, ikp) = 2.0_dp
               matrix_lin_eqs(ikp, nkp + 1) = 1.0_dp
               matrix_lin_eqs(nkp + 1, ikp) = 1.0_dp

               matrix_lin_eqs(ikp, nkp + 2) = k_sq**(exp_kpoints*0.5_dp)
               matrix_lin_eqs(nkp + 2, ikp) = k_sq**(exp_kpoints*0.5_dp)

            END DO

            CALL invmat(matrix_lin_eqs, info)
            ! check whether inversion was successful
            CPASSERT(info == 0)

            ALLOCATE (right_side(num_lin_eqs))
            right_side = 0.0_dp
            right_side(nkp + 1) = 1.0_dp
            ! divide integral by two because CP2K k-mesh already considers symmetry k <-> -k
            right_side(nkp + 2) = integral

            ALLOCATE (wkp_tmp(num_lin_eqs))

            wkp_tmp(1:num_lin_eqs) = MATMUL(matrix_lin_eqs, right_side)

            wkp_V(1:nkp) = wkp_tmp(1:nkp)

            DEALLOCATE (matrix_lin_eqs, right_side, wkp_tmp)

            wkp_W(:) = wkp_V(:)

         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
!> \param cfm_mat_work ...
! **************************************************************************************************
   SUBROUTINE own_cfm_upper_to_full(cfm_mat_Q, cfm_mat_work)

      TYPE(cp_cfm_type), INTENT(IN), POINTER             :: cfm_mat_Q, cfm_mat_work

      CHARACTER(LEN=*), PARAMETER :: routineN = 'own_cfm_upper_to_full'

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      ! get info of fm_mat_Q
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global < i_global) THEN
               cfm_mat_Q%local_data(iiB, jjB) = z_zero
            END IF
            IF (j_global == i_global) THEN
               cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB)/(2.0_dp, 0.0_dp)
            END IF
         END DO
      END DO

      CALL cp_cfm_transpose(cfm_mat_Q, 'C', cfm_mat_work)

      CALL cp_cfm_scale_and_add(z_one, cfm_mat_Q, z_one, cfm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE rpa_gw_kpoints
