!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines treating GW and RPA calculations with kpoints
!> \par History
!>      since 2018 continuous development [J. Wilhelm]
! **************************************************************************************************
MODULE rpa_gw_kpoints_util
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_cholesky_invert,&
                                              cp_cfm_column_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_transpose
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig,&
                                              cp_cfm_geeig_canon,&
                                              cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_filter, &
        dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_set, dbcsr_transposed, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE hfx_types,                       ONLY: hfx_release
   USE input_constants,                 ONLY: cholesky_off,&
                                              kp_weights_W_auto,&
                                              kp_weights_W_tailored,&
                                              kp_weights_W_uniform
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_band_structure,               ONLY: calculate_kpoints_for_bs
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE rpa_gw_im_time_util,             ONLY: compute_weight_re_im,&
                                              get_atom_index_from_basis_function_index
   USE rpa_im_time,                     ONLY: init_cell_index_rpa
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_kpoints_util'

   PUBLIC :: invert_eps_compute_W_and_Erpa_kp, cp_cfm_power, real_space_to_kpoint_transform_rpa, &
             get_mat_cell_T_from_mat_gamma, get_bandstruc_and_k_dependent_MOs, &
             compute_wkp_W, mat_kp_from_mat_gamma, cp_cfm_upper_to_full

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dimen_RI ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param nkp ...
!> \param count_ev_sc_GW ...
!> \param para_env ...
!> \param Erpa ...
!> \param tau_tj ...
!> \param tj ...
!> \param wj ...
!> \param weights_cos_tf_w_to_t ...
!> \param wkp_W ...
!> \param do_gw_im_time ...
!> \param do_ri_Sigma_x ...
!> \param do_kpoints_from_Gamma ...
!> \param cfm_mat_Q ...
!> \param ikp_local ...
!> \param mat_P_omega ...
!> \param mat_P_omega_kp ...
!> \param qs_env ...
!> \param eps_filter_im_time ...
!> \param unit_nr ...
!> \param kpoints ...
!> \param fm_mat_Minv_L_kpoints ...
!> \param fm_matrix_L_kpoints ...
!> \param fm_mat_W ...
!> \param fm_mat_RI_global_work ...
!> \param mat_MinvVMinv ...
!> \param fm_matrix_Minv ...
!> \param fm_matrix_Minv_Vtrunc_Minv ...
! **************************************************************************************************
   SUBROUTINE invert_eps_compute_W_and_Erpa_kp(dimen_RI, num_integ_points, jquad, nkp, count_ev_sc_GW, para_env, &
                                               Erpa, tau_tj, tj, wj, weights_cos_tf_w_to_t, wkp_W, do_gw_im_time, &
                                               do_ri_Sigma_x, do_kpoints_from_Gamma, &
                                               cfm_mat_Q, ikp_local, mat_P_omega, mat_P_omega_kp, &
                                               qs_env, eps_filter_im_time, unit_nr, kpoints, fm_mat_Minv_L_kpoints, &
                                               fm_matrix_L_kpoints, fm_mat_W, &
                                               fm_mat_RI_global_work, mat_MinvVMinv, fm_matrix_Minv, &
                                               fm_matrix_Minv_Vtrunc_Minv)

      INTEGER, INTENT(IN)                                :: dimen_RI, num_integ_points, jquad, nkp, &
                                                            count_ev_sc_GW
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(INOUT)                       :: Erpa
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tj, wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_w_to_t
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wkp_W
      LOGICAL, INTENT(IN)                                :: do_gw_im_time, do_ri_Sigma_x, &
                                                            do_kpoints_from_Gamma
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mat_Q
      INTEGER, INTENT(IN)                                :: ikp_local
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: mat_P_omega, mat_P_omega_kp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_mat_Minv_L_kpoints, &
                                                            fm_matrix_L_kpoints
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_W
      TYPE(cp_fm_type)                                   :: fm_mat_RI_global_work
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_MinvVMinv
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_matrix_Minv, &
                                                            fm_matrix_Minv_Vtrunc_Minv

      CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_eps_compute_W_and_Erpa_kp'

      INTEGER                                            :: handle, ikp
      LOGICAL                                            :: do_this_ikp
      REAL(KIND=dp)                                      :: t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: trace_Qomega

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      IF (do_kpoints_from_Gamma) THEN
         CALL get_mat_cell_T_from_mat_gamma(mat_P_omega(jquad, :), qs_env, kpoints, jquad, unit_nr)
      END IF

      CALL transform_P_from_real_space_to_kpoints(mat_P_omega, mat_P_omega_kp, &
                                                  kpoints, eps_filter_im_time, jquad)

      ALLOCATE (trace_Qomega(dimen_RI))

      IF (unit_nr > 0) WRITE (unit_nr, '(/T3,A,1X,I3)') &
         'GW_INFO| Computing chi and W frequency point', jquad

      DO ikp = 1, nkp

         ! parallization, we either have all kpoints on all processors or a single kpoint per group
         do_this_ikp = (ikp_local == -1) .OR. (ikp_local == 0 .AND. ikp == 1) .OR. (ikp_local == ikp)
         IF (.NOT. do_this_ikp) CYCLE

         ! 1. remove all spurious negative eigenvalues from P(iw,k), multiplication Q(iw,k) = K^H(k)P(iw,k)K(k)
         CALL compute_Q_kp_RPA(cfm_mat_Q, &
                               mat_P_omega_kp, &
                               fm_mat_Minv_L_kpoints(ikp, 1), &
                               fm_mat_Minv_L_kpoints(ikp, 2), &
                               fm_mat_RI_global_work, &
                               dimen_RI, ikp, nkp, ikp_local, para_env, &
                               qs_env%mp2_env%ri_rpa_im_time%make_chi_pos_definite)

         ! 2. Cholesky decomposition of Id + Q(iw,k)
         CALL cholesky_decomp_Q(cfm_mat_Q, para_env, trace_Qomega, dimen_RI)

         ! 3. Computing E_c^RPA = E_c^RPA + a_w/N_k*sum_k ln[det(1+Q(iw,k))-Tr(Q(iw,k))]
         CALL frequency_and_kpoint_integration(Erpa, cfm_mat_Q, para_env, trace_Qomega, &
                                               dimen_RI, wj(jquad), kpoints%wkp(ikp))

         IF (do_gw_im_time) THEN

            ! compute S^-1*V*S^-1 for exchange part of the self-energy in real space as W in real space
            IF (do_ri_Sigma_x .AND. jquad == 1 .AND. count_ev_sc_GW == 1 .AND. do_kpoints_from_Gamma) THEN

               CALL dbcsr_set(mat_MinvVMinv%matrix, 0.0_dp)
               CALL copy_fm_to_dbcsr(fm_matrix_Minv_Vtrunc_Minv(1, 1), mat_MinvVMinv%matrix, keep_sparsity=.FALSE.)

            END IF
            IF (do_kpoints_from_Gamma) THEN
               CALL compute_Wc_real_space_tau_GW(fm_mat_W, cfm_mat_Q, &
                                                 fm_matrix_L_kpoints(ikp, 1), &
                                                 fm_matrix_L_kpoints(ikp, 2), &
                                                 dimen_RI, num_integ_points, jquad, &
                                                 ikp, tj, tau_tj, weights_cos_tf_w_to_t, &
                                                 ikp_local, para_env, kpoints, qs_env, wkp_W)
            END IF

         END IF
      END DO

      ! after the transform of (eps(iw)-1)^-1 from iw to it is done, multiply with V^1/2 to obtain W(it)
      IF (do_gw_im_time .AND. do_kpoints_from_Gamma .AND. jquad == num_integ_points) THEN
         CALL Wc_to_Minv_Wc_Minv(fm_mat_W, fm_matrix_Minv, para_env, dimen_RI, num_integ_points)
         CALL deallocate_kp_matrices(fm_matrix_L_kpoints, fm_mat_Minv_L_kpoints)
      END IF

      DEALLOCATE (trace_Qomega)

      t2 = m_walltime()

      IF (unit_nr > 0) WRITE (unit_nr, '(T6,A,T56,F25.1)') 'Execution time (s):', t2 - t1

      CALL timestop(handle)

   END SUBROUTINE invert_eps_compute_W_and_Erpa_kp

! **************************************************************************************************
!> \brief ...
!> \param fm_matrix_L_kpoints ...
!> \param fm_mat_Minv_L_kpoints ...
! **************************************************************************************************
   SUBROUTINE deallocate_kp_matrices(fm_matrix_L_kpoints, fm_mat_Minv_L_kpoints)

      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_matrix_L_kpoints, &
                                                            fm_mat_Minv_L_kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_kp_matrices'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_release(fm_mat_Minv_L_kpoints)
      CALL cp_fm_release(fm_matrix_L_kpoints)

      CALL timestop(handle)

   END SUBROUTINE deallocate_kp_matrices

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param threshold ...
!> \param exponent ...
!> \param min_eigval ...
! **************************************************************************************************
   SUBROUTINE cp_cfm_power(matrix, threshold, exponent, min_eigval)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: matrix
      REAL(KIND=dp)                                      :: threshold, exponent
      REAL(KIND=dp), OPTIONAL                            :: min_eigval

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cp_cfm_power'
      COMPLEX(KIND=dp), PARAMETER :: czero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: eigenvalues_exponent
      INTEGER                                            :: handle, i, ncol_global, nrow_global
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_cfm_type)                                  :: cfm_work

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_work, matrix%matrix_struct)
      CALL cp_cfm_set_all(cfm_work, z_zero)

      ! Test that matrix is square
      CALL cp_cfm_get_info(matrix, nrow_global=nrow_global, ncol_global=ncol_global)
      CPASSERT(nrow_global == ncol_global)
      ALLOCATE (eigenvalues(nrow_global))
      eigenvalues(:) = 0.0_dp
      ALLOCATE (eigenvalues_exponent(nrow_global))
      eigenvalues_exponent(:) = czero

      ! Diagonalize matrix: get eigenvectors and eigenvalues
      CALL cp_cfm_heevd(matrix, cfm_work, eigenvalues)

      DO i = 1, nrow_global
         IF (eigenvalues(i) > threshold) THEN
            eigenvalues_exponent(i) = CMPLX((eigenvalues(i))**(0.5_dp*exponent), threshold, KIND=dp)
         ELSE
            IF (PRESENT(min_eigval)) THEN
               eigenvalues_exponent(i) = CMPLX(min_eigval, 0.0_dp, KIND=dp)
            ELSE
               eigenvalues_exponent(i) = czero
            END IF
         END IF
      END DO

      CALL cp_cfm_column_scale(cfm_work, eigenvalues_exponent)

      CALL parallel_gemm("N", "C", nrow_global, nrow_global, nrow_global, z_one, &
                         cfm_work, cfm_work, z_zero, matrix)

      DEALLOCATE (eigenvalues, eigenvalues_exponent)

      CALL cp_cfm_release(cfm_work)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_power

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
!> \param mat_P_omega_kp ...
!> \param fm_mat_L_re ...
!> \param fm_mat_L_im ...
!> \param fm_mat_RI_global_work ...
!> \param dimen_RI ...
!> \param ikp ...
!> \param nkp ...
!> \param ikp_local ...
!> \param para_env ...
!> \param make_chi_pos_definite ...
! **************************************************************************************************
   SUBROUTINE compute_Q_kp_RPA(cfm_mat_Q, mat_P_omega_kp, fm_mat_L_re, fm_mat_L_im, &
                               fm_mat_RI_global_work, dimen_RI, ikp, nkp, ikp_local, para_env, &
                               make_chi_pos_definite)

      TYPE(cp_cfm_type)                                  :: cfm_mat_Q
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: mat_P_omega_kp
      TYPE(cp_fm_type)                                   :: fm_mat_L_re, fm_mat_L_im, &
                                                            fm_mat_RI_global_work
      INTEGER, INTENT(IN)                                :: dimen_RI, ikp, nkp, ikp_local
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: make_chi_pos_definite

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Q_kp_RPA'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_mat_L, cfm_mat_work
      TYPE(cp_fm_type)                                   :: fm_mat_work

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_mat_work, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      CALL cp_cfm_create(cfm_mat_L, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_L, z_zero)

      CALL cp_fm_create(fm_mat_work, fm_mat_L_re%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

      ! 1. Convert the dbcsr matrix mat_P_omega_kp (that is chi(k,iw)) to a full matrix and
      !    distribute it to subgroups
      CALL mat_P_to_subgroup(mat_P_omega_kp, fm_mat_RI_global_work, &
                             fm_mat_work, cfm_mat_Q, ikp, nkp, ikp_local, para_env)

      ! 2. Remove all negative eigenvalues from chi(k,iw)
      IF (make_chi_pos_definite) THEN
         CALL cp_cfm_power(cfm_mat_Q, threshold=0.0_dp, exponent=1.0_dp)
      END IF

      ! 3. Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      ! 4. work = P(iw,k)*L(k)
      CALL parallel_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                         z_zero, cfm_mat_work)

      ! 5. Q(iw,k) = L^H(k)*work
      CALL parallel_gemm('C', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_work, &
                         z_zero, cfm_mat_Q)

      CALL cp_cfm_release(cfm_mat_work)
      CALL cp_cfm_release(cfm_mat_L)
      CALL cp_fm_release(fm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE compute_Q_kp_RPA

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega_kp ...
!> \param fm_mat_RI_global_work ...
!> \param fm_mat_work ...
!> \param cfm_mat_Q ...
!> \param ikp ...
!> \param nkp ...
!> \param ikp_local ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE mat_P_to_subgroup(mat_P_omega_kp, fm_mat_RI_global_work, &
                                fm_mat_work, cfm_mat_Q, ikp, nkp, ikp_local, para_env)

      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: mat_P_omega_kp
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_RI_global_work, fm_mat_work
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mat_Q
      INTEGER, INTENT(IN)                                :: ikp, nkp, ikp_local
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mat_P_to_subgroup'

      INTEGER                                            :: handle, jkp
      TYPE(cp_fm_type)                                   :: fm_dummy
      TYPE(dbcsr_type), POINTER                          :: mat_P_omega_im, mat_P_omega_re

      CALL timeset(routineN, handle)

      IF (ikp_local == -1) THEN

         mat_P_omega_re => mat_P_omega_kp(1, ikp)%matrix
         CALL cp_fm_set_all(fm_mat_work, 0.0_dp)
         CALL copy_dbcsr_to_fm(mat_P_omega_re, fm_mat_work)
         CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_Q, z_one, fm_mat_work)

         mat_P_omega_im => mat_P_omega_kp(2, ikp)%matrix
         CALL cp_fm_set_all(fm_mat_work, 0.0_dp)
         CALL copy_dbcsr_to_fm(mat_P_omega_im, fm_mat_work)
         CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_Q, gaussi, fm_mat_work)

      ELSE

         CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

         DO jkp = 1, nkp

            mat_P_omega_re => mat_P_omega_kp(1, jkp)%matrix

            CALL cp_fm_set_all(fm_mat_RI_global_work, 0.0_dp)
            CALL copy_dbcsr_to_fm(mat_P_omega_re, fm_mat_RI_global_work)

            CALL para_env%sync()

            IF (ikp_local == jkp) THEN
               CALL cp_fm_copy_general(fm_mat_RI_global_work, fm_mat_work, para_env)
            ELSE
               CALL cp_fm_copy_general(fm_mat_RI_global_work, fm_dummy, para_env)
            END IF

            CALL para_env%sync()

         END DO

         CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_Q, z_one, fm_mat_work)

         CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

         DO jkp = 1, nkp

            mat_P_omega_im => mat_P_omega_kp(2, jkp)%matrix

            CALL cp_fm_set_all(fm_mat_RI_global_work, 0.0_dp)
            CALL copy_dbcsr_to_fm(mat_P_omega_im, fm_mat_RI_global_work)

            CALL para_env%sync()

            IF (ikp_local == jkp) THEN
               CALL cp_fm_copy_general(fm_mat_RI_global_work, fm_mat_work, para_env)
            ELSE
               CALL cp_fm_copy_general(fm_mat_RI_global_work, fm_dummy, para_env)
            END IF

            CALL para_env%sync()

         END DO

         CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_Q, gaussi, fm_mat_work)

         CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

      END IF

      CALL para_env%sync()

      CALL timestop(handle)

   END SUBROUTINE mat_P_to_subgroup

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
!> \param para_env ...
!> \param trace_Qomega ...
!> \param dimen_RI ...
! **************************************************************************************************
   SUBROUTINE cholesky_decomp_Q(cfm_mat_Q, para_env, trace_Qomega, dimen_RI)

      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mat_Q
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: trace_Qomega
      INTEGER, INTENT(IN)                                :: dimen_RI

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cholesky_decomp_Q'

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_cfm_type)                                  :: cfm_mat_Q_tmp, cfm_mat_work

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_mat_work, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      CALL cp_cfm_create(cfm_mat_Q_tmp, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_Q_tmp, z_zero)

      ! get info of fm_mat_Q
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      ! calculate the trace of Q and add 1 on the diagonal
      trace_Qomega = 0.0_dp
!$OMP     PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                 SHARED(ncol_local,nrow_local,col_indices,row_indices,trace_Qomega,cfm_mat_Q,dimen_RI)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               trace_Qomega(i_global) = REAL(cfm_mat_Q%local_data(iiB, jjB))
               cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB) + z_one
            END IF
         END DO
      END DO
      CALL para_env%sum(trace_Qomega)

      CALL cp_cfm_to_cfm(cfm_mat_Q, cfm_mat_Q_tmp)

      CALL cp_cfm_cholesky_decompose(matrix=cfm_mat_Q, n=dimen_RI, info_out=info_chol)

      CPASSERT(info_chol == 0)

      CALL cp_cfm_release(cfm_mat_work)
      CALL cp_cfm_release(cfm_mat_Q_tmp)

      CALL timestop(handle)

   END SUBROUTINE cholesky_decomp_Q

! **************************************************************************************************
!> \brief ...
!> \param Erpa ...
!> \param cfm_mat_Q ...
!> \param para_env ...
!> \param trace_Qomega ...
!> \param dimen_RI ...
!> \param freq_weight ...
!> \param kp_weight ...
! **************************************************************************************************
   SUBROUTINE frequency_and_kpoint_integration(Erpa, cfm_mat_Q, para_env, trace_Qomega, &
                                               dimen_RI, freq_weight, kp_weight)

      REAL(KIND=dp), INTENT(INOUT)                       :: Erpa
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mat_Q
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: trace_Qomega
      INTEGER, INTENT(IN)                                :: dimen_RI
      REAL(KIND=dp), INTENT(IN)                          :: freq_weight, kp_weight

      CHARACTER(LEN=*), PARAMETER :: routineN = 'frequency_and_kpoint_integration'

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: FComega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Q_log

      CALL timeset(routineN, handle)

      ! get info of cholesky_decomposed(fm_mat_Q)
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      ALLOCATE (Q_log(dimen_RI))
      Q_log = 0.0_dp
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                SHARED(ncol_local,nrow_local,col_indices,row_indices,Q_log,cfm_mat_Q,dimen_RI)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               Q_log(i_global) = 2.0_dp*LOG(REAL(cfm_mat_Q%local_data(iiB, jjB)))
            END IF
         END DO
      END DO
      CALL para_env%sum(Q_log)

      FComega = 0.0_dp
      DO iiB = 1, dimen_RI
         IF (MODULO(iiB, para_env%num_pe) /= para_env%mepos) CYCLE
         ! FComega=FComega+(LOG(Q_log(iiB))-trace_Qomega(iiB))/2.0_dp
         FComega = FComega + (Q_log(iiB) - trace_Qomega(iiB))/2.0_dp
      END DO

      Erpa = Erpa + FComega*freq_weight*kp_weight

      DEALLOCATE (Q_log)

      CALL timestop(handle)

   END SUBROUTINE frequency_and_kpoint_integration

! **************************************************************************************************
!> \brief ...
!> \param tj_dummy ...
!> \param tau_tj_dummy ...
!> \param weights_cos_tf_w_to_t_dummy ...
! **************************************************************************************************
   SUBROUTINE get_dummys(tj_dummy, tau_tj_dummy, weights_cos_tf_w_to_t_dummy)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: tj_dummy, tau_tj_dummy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_w_to_t_dummy

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_dummys'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (weights_cos_tf_w_to_t_dummy(1, 1))
      ALLOCATE (tj_dummy(1))
      ALLOCATE (tau_tj_dummy(1))

      tj_dummy(1) = 0.0_dp
      tau_tj_dummy(1) = 0.0_dp
      weights_cos_tf_w_to_t_dummy(1, 1) = 1.0_dp

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param tj_dummy ...
!> \param tau_tj_dummy ...
!> \param weights_cos_tf_w_to_t_dummy ...
! **************************************************************************************************
   SUBROUTINE release_dummys(tj_dummy, tau_tj_dummy, weights_cos_tf_w_to_t_dummy)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: tj_dummy, tau_tj_dummy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_w_to_t_dummy

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'release_dummys'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      DEALLOCATE (weights_cos_tf_w_to_t_dummy)
      DEALLOCATE (tj_dummy)
      DEALLOCATE (tau_tj_dummy)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega ...
!> \param qs_env ...
!> \param kpoints ...
!> \param jquad ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE get_mat_cell_T_from_mat_gamma(mat_P_omega, qs_env, kpoints, jquad, unit_nr)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: mat_P_omega
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: jquad, unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_cell_T_from_mat_gamma'

      INTEGER                                            :: col, handle, i_cell, i_dim, j_cell, &
                                                            num_cells_P, num_integ_points, row
      INTEGER, DIMENSION(3)                              :: cell_grid_P, periodic
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_P
      LOGICAL :: i_cell_is_the_minimum_image_cell
      REAL(KIND=dp)                                      :: abs_rab_cell_i, abs_rab_cell_j
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector, cell_vector_j, rab_cell_i, &
                                                            rab_cell_j
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, &
                      particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat, periodic=periodic)

      DO i_dim = 1, 3
         ! we have at most 3 neigboring cells per dimension and at least one because
         ! the density response at Gamma is only divided to neighboring
         IF (periodic(i_dim) == 1) THEN
            cell_grid_P(i_dim) = MAX(MIN((kpoints%nkp_grid(i_dim)/2)*2 - 1, 1), 3)
         ELSE
            cell_grid_P(i_dim) = 1
         END IF
      END DO

      ! overwrite the cell indices in kpoints
      CALL init_cell_index_rpa(cell_grid_P, kpoints%cell_to_index, kpoints%index_to_cell, cell)

      index_to_cell_P => kpoints%index_to_cell

      num_cells_P = SIZE(index_to_cell_P, 2)

      num_integ_points = SIZE(mat_P_omega, 1)

      ! first, copy the Gamma-only result from mat_P_omega(1) into all other matrices and
      ! remove the blocks later which do not belong to the cell index
      DO i_cell = 2, num_cells_P
         CALL dbcsr_copy(mat_P_omega(i_cell)%matrix, &
                         mat_P_omega(1)%matrix)
      END DO

      IF (jquad == 1 .AND. unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A,T66,ES15.2)') 'GW_INFO| RI regularization parameter: ', &
            qs_env%mp2_env%ri_rpa_im_time%regularization_RI
         WRITE (unit_nr, '(T3,A,T66,ES15.2)') 'GW_INFO| eps_eigval_S: ', &
            qs_env%mp2_env%ri_rpa_im_time%eps_eigval_S
         IF (qs_env%mp2_env%ri_rpa_im_time%make_chi_pos_definite) THEN
            WRITE (unit_nr, '(T3,A,T81)') &
               'GW_INFO| Make chi(iw,k) positive definite?                                TRUE'
         ELSE
            WRITE (unit_nr, '(T3,A,T81)') &
               'GW_INFO| Make chi(iw,k) positive definite?                               FALSE'
         END IF

      END IF

      DO i_cell = 1, num_cells_P

         CALL dbcsr_iterator_start(iter, mat_P_omega(i_cell)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, data_block)

            cell_vector(1:3) = MATMUL(hmat, REAL(index_to_cell_P(1:3, i_cell), dp))
            rab_cell_i(1:3) = pbc(particle_set(row)%r(1:3), cell) - &
                              (pbc(particle_set(col)%r(1:3), cell) + cell_vector(1:3))
            abs_rab_cell_i = SQRT(rab_cell_i(1)**2 + rab_cell_i(2)**2 + rab_cell_i(3)**2)

            ! minimum image convention
            i_cell_is_the_minimum_image_cell = .TRUE.
            DO j_cell = 1, num_cells_P
               cell_vector_j(1:3) = MATMUL(hmat, REAL(index_to_cell_P(1:3, j_cell), dp))
               rab_cell_j(1:3) = pbc(particle_set(row)%r(1:3), cell) - &
                                 (pbc(particle_set(col)%r(1:3), cell) + cell_vector_j(1:3))
               abs_rab_cell_j = SQRT(rab_cell_j(1)**2 + rab_cell_j(2)**2 + rab_cell_j(3)**2)

               IF (abs_rab_cell_i > abs_rab_cell_j + 1.0E-6_dp) THEN
                  i_cell_is_the_minimum_image_cell = .FALSE.
               END IF
            END DO

            IF (.NOT. i_cell_is_the_minimum_image_cell) THEN
               data_block(:, :) = data_block(:, :)*0.0_dp
            END IF

         END DO
         CALL dbcsr_iterator_stop(iter)

      END DO

      CALL timestop(handle)

   END SUBROUTINE get_mat_cell_T_from_mat_gamma

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega ...
!> \param mat_P_omega_kp ...
!> \param kpoints ...
!> \param eps_filter_im_time ...
!> \param jquad ...
! **************************************************************************************************
   SUBROUTINE transform_P_from_real_space_to_kpoints(mat_P_omega, mat_P_omega_kp, &
                                                     kpoints, eps_filter_im_time, jquad)

      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: mat_P_omega, mat_P_omega_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(kind=dp), INTENT(IN)                          :: eps_filter_im_time
      INTEGER, INTENT(IN)                                :: jquad

      CHARACTER(LEN=*), PARAMETER :: routineN = 'transform_P_from_real_space_to_kpoints'

      INTEGER                                            :: handle, icell, nkp, num_integ_points

      CALL timeset(routineN, handle)

      num_integ_points = SIZE(mat_P_omega, 1)
      nkp = SIZE(mat_P_omega, 2)

      CALL real_space_to_kpoint_transform_rpa(mat_P_omega_kp(1, :), mat_P_omega_kp(2, :), mat_P_omega(jquad, :), &
                                              kpoints, eps_filter_im_time)

      DO icell = 1, SIZE(mat_P_omega, 2)
         CALL dbcsr_set(mat_P_omega(jquad, icell)%matrix, 0.0_dp)
         CALL dbcsr_filter(mat_P_omega(jquad, icell)%matrix, 1.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE transform_P_from_real_space_to_kpoints

! **************************************************************************************************
!> \brief ...
!> \param real_mat_kp ...
!> \param imag_mat_kp ...
!> \param mat_real_space ...
!> \param kpoints ...
!> \param eps_filter_im_time ...
!> \param real_mat_real_space ...
! **************************************************************************************************
   SUBROUTINE real_space_to_kpoint_transform_rpa(real_mat_kp, imag_mat_kp, mat_real_space, &
                                                 kpoints, eps_filter_im_time, real_mat_real_space)

      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: real_mat_kp, imag_mat_kp, mat_real_space
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time
      LOGICAL, INTENT(IN), OPTIONAL                      :: real_mat_real_space

      CHARACTER(LEN=*), PARAMETER :: routineN = 'real_space_to_kpoint_transform_rpa'

      INTEGER                                            :: handle, i_cell, ik, nkp, num_cells
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      LOGICAL                                            :: my_real_mat_real_space
      REAL(KIND=dp)                                      :: arg, coskl, sinkl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(dbcsr_type)                                   :: mat_work

      CALL timeset(routineN, handle)

      my_real_mat_real_space = .TRUE.
      IF (PRESENT(real_mat_real_space)) my_real_mat_real_space = real_mat_real_space

      CALL dbcsr_create(matrix=mat_work, &
                        template=real_mat_kp(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_reserve_all_blocks(mat_work)
      CALL dbcsr_set(mat_work, 0.0_dp)

      ! this kpoint environme t should be the kpoints for D(it) and X(it) created in init_cell_index_rpa
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp)

      NULLIFY (index_to_cell)
      index_to_cell => kpoints%index_to_cell

      num_cells = SIZE(index_to_cell, 2)

      CPASSERT(SIZE(mat_real_space) >= num_cells/2 + 1)

      DO ik = 1, nkp

         CALL dbcsr_reserve_all_blocks(real_mat_kp(ik)%matrix)
         CALL dbcsr_reserve_all_blocks(imag_mat_kp(ik)%matrix)

         CALL dbcsr_set(real_mat_kp(ik)%matrix, 0.0_dp)
         CALL dbcsr_set(imag_mat_kp(ik)%matrix, 0.0_dp)

         DO i_cell = 1, num_cells/2 + 1

            cell(:) = index_to_cell(:, i_cell)

            arg = REAL(cell(1), dp)*xkp(1, ik) + REAL(cell(2), dp)*xkp(2, ik) + REAL(cell(3), dp)*xkp(3, ik)
            coskl = COS(twopi*arg)
            sinkl = SIN(twopi*arg)

            IF (my_real_mat_real_space) THEN
               CALL dbcsr_add_local(real_mat_kp(ik)%matrix, mat_real_space(i_cell)%matrix, 1.0_dp, coskl)
               CALL dbcsr_add_local(imag_mat_kp(ik)%matrix, mat_real_space(i_cell)%matrix, 1.0_dp, sinkl)
            ELSE
               CALL dbcsr_add_local(real_mat_kp(ik)%matrix, mat_real_space(i_cell)%matrix, 1.0_dp, -sinkl)
               CALL dbcsr_add_local(imag_mat_kp(ik)%matrix, mat_real_space(i_cell)%matrix, 1.0_dp, coskl)
            END IF

            IF (.NOT. (cell(1) == 0 .AND. cell(2) == 0 .AND. cell(3) == 0)) THEN

               CALL dbcsr_transposed(mat_work, mat_real_space(i_cell)%matrix)

               IF (my_real_mat_real_space) THEN
                  CALL dbcsr_add_local(real_mat_kp(ik)%matrix, mat_work, 1.0_dp, coskl)
                  CALL dbcsr_add_local(imag_mat_kp(ik)%matrix, mat_work, 1.0_dp, -sinkl)
               ELSE
                  ! for an imaginary real-space matrix, we need to consider the imaginary unit
                  ! and we need to take into account that the transposed gives an extra "-" sign
                  ! because the transposed is actually Hermitian conjugate
                  CALL dbcsr_add_local(real_mat_kp(ik)%matrix, mat_work, 1.0_dp, -sinkl)
                  CALL dbcsr_add_local(imag_mat_kp(ik)%matrix, mat_work, 1.0_dp, -coskl)
               END IF

               CALL dbcsr_set(mat_work, 0.0_dp)

            END IF

         END DO

         CALL dbcsr_filter(real_mat_kp(ik)%matrix, eps_filter_im_time)
         CALL dbcsr_filter(imag_mat_kp(ik)%matrix, eps_filter_im_time)

      END DO

      CALL dbcsr_release(mat_work)

      CALL timestop(handle)

   END SUBROUTINE real_space_to_kpoint_transform_rpa

! **************************************************************************************************
!> \brief ...
!> \param mat_a ...
!> \param mat_b ...
!> \param alpha ...
!> \param beta ...
! **************************************************************************************************
   SUBROUTINE dbcsr_add_local(mat_a, mat_b, alpha, beta)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_a, mat_b
      REAL(kind=dp), INTENT(IN)                          :: alpha, beta

      INTEGER                                            :: col, row
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_to_compute, data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_iterator_start(iter, mat_b)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         NULLIFY (block_to_compute)
         CALL dbcsr_get_block_p(matrix=mat_a, &
                                row=row, col=col, block=block_to_compute, found=found)

         CPASSERT(found)

         block_to_compute(:, :) = alpha*block_to_compute(:, :) + beta*data_block(:, :)

      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE dbcsr_add_local

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
! **************************************************************************************************
   SUBROUTINE compute_Wc_real_space_tau_GW(fm_mat_W_tau, cfm_mat_Q, fm_mat_L_re, fm_mat_L_im, &
                                           dimen_RI, num_integ_points, jquad, &
                                           ikp, tj, tau_tj, weights_cos_tf_w_to_t, ikp_local, &
                                           para_env, kpoints, qs_env, wkp_W)

      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_W_tau
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mat_Q
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_L_re, fm_mat_L_im
      INTEGER, INTENT(IN)                                :: dimen_RI, num_integ_points, jquad, ikp
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tj
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weights_cos_tf_w_to_t
      INTEGER, INTENT(IN)                                :: ikp_local
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wkp_W

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Wc_real_space_tau_GW'

      INTEGER :: handle, handle2, i_global, iatom, iatom_old, iiB, iquad, irow, j_global, jatom, &
         jatom_old, jcol, jjB, jkp, ncol_local, nkp, nrow_local, num_cells
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      REAL(KIND=dp)                                      :: contribution, omega, tau, weight, &
                                                            weight_im, weight_re
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_type)                                  :: cfm_mat_L, cfm_mat_work, cfm_mat_work_2
      TYPE(cp_fm_type)                                   :: fm_dummy, fm_mat_work_global, &
                                                            fm_mat_work_local
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL timeset(routineN//"_1", handle2)

      CALL cp_cfm_create(cfm_mat_work, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      CALL cp_cfm_create(cfm_mat_work_2, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work_2, z_zero)

      CALL cp_cfm_create(cfm_mat_L, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_L, z_zero)

      ! Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      CALL cp_fm_create(fm_mat_work_global, fm_mat_W_tau(1)%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_global, 0.0_dp)

      CALL cp_fm_create(fm_mat_work_local, cfm_mat_Q%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_local, 0.0_dp)

      CALL timestop(handle2)

      CALL timeset(routineN//"_2", handle2)

      ! calculate [1+Q(iw')]^-1
      CALL cp_cfm_cholesky_invert(cfm_mat_Q)

      ! symmetrize the result
      CALL cp_cfm_upper_to_full(cfm_mat_Q)

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

      CALL timeset(routineN//"_3", handle2)

      ! work = epsilon(iw,k)*V^1/2(k)
      CALL parallel_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                         z_zero, cfm_mat_work)

      ! W(iw,k) = V^1/2(k)*work
      CALL parallel_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_work, &
                         z_zero, cfm_mat_work_2)

      CALL timestop(handle2)

      CALL timeset(routineN//"_4", handle2)

      CALL get_kpoint_info(kpoints, xkp=xkp, wkp=wkp, nkp=nkp)
      index_to_cell => kpoints%index_to_cell
      num_cells = SIZE(index_to_cell, 2)

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

               ! symmetrize=.FALSE. necessary since we already have a symmetrized index_to_cell
               CALL compute_weight_re_im(weight_re, weight_im, &
                                         num_cells, iatom, jatom, xkp(1:3, ikp), wkp_W(ikp), &
                                         cell, index_to_cell, hmat, particle_set)

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

      IF (ikp_local == -1) THEN

         CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)

         DO iquad = 1, num_integ_points

            omega = tj(jquad)
            tau = tau_tj(iquad)
            weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)

            IF (jquad == 1 .AND. ikp == 1) THEN
               CALL cp_fm_set_all(matrix=fm_mat_W_tau(iquad), alpha=0.0_dp)
            END IF

            CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_W_tau(iquad), beta=weight, matrix_b=fm_mat_work_global)

         END DO

      ELSE

         DO jkp = 1, nkp

            CALL para_env%sync()

            IF (ikp_local == jkp) THEN
               CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)
            ELSE
               CALL cp_fm_copy_general(fm_dummy, fm_mat_work_global, para_env)
            END IF

            CALL para_env%sync()

            DO iquad = 1, num_integ_points

               omega = tj(jquad)
               tau = tau_tj(iquad)
               weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)

               IF (jquad == 1 .AND. jkp == 1) THEN
                  CALL cp_fm_set_all(matrix=fm_mat_W_tau(iquad), alpha=0.0_dp)
               END IF

               CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_W_tau(iquad), beta=weight, &
                                        matrix_b=fm_mat_work_global)

            END DO

         END DO

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
!> \param fm_mat_W ...
!> \param fm_matrix_Minv ...
!> \param para_env ...
!> \param dimen_RI ...
!> \param num_integ_points ...
! **************************************************************************************************
   SUBROUTINE Wc_to_Minv_Wc_Minv(fm_mat_W, fm_matrix_Minv, para_env, dimen_RI, num_integ_points)
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_mat_W
      TYPE(cp_fm_type), DIMENSION(:, :)                  :: fm_matrix_Minv
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER                                            :: dimen_RI, num_integ_points

      CHARACTER(LEN=*), PARAMETER :: routineN = 'Wc_to_Minv_Wc_Minv'

      INTEGER                                            :: handle, jquad
      TYPE(cp_fm_type)                                   :: fm_work_Minv, fm_work_Minv_W

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_work_Minv, fm_mat_W(1)%matrix_struct)
      CALL cp_fm_copy_general(fm_matrix_Minv(1, 1), fm_work_Minv, para_env)

      CALL cp_fm_create(fm_work_Minv_W, fm_mat_W(1)%matrix_struct)

      DO jquad = 1, num_integ_points

         CALL parallel_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_work_Minv, fm_mat_W(jquad), &
                            0.0_dp, fm_work_Minv_W)
         CALL parallel_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_work_Minv_W, fm_work_Minv, &
                            0.0_dp, fm_mat_W(jquad))

      END DO

      CALL cp_fm_release(fm_work_Minv)

      CALL cp_fm_release(fm_work_Minv_W)

      CALL timestop(handle)

   END SUBROUTINE Wc_to_Minv_Wc_Minv

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

      ! for exchange part of self-energy, we use truncated Coulomb operator that should be fine
      ! with uniform weights (without k-point extrapolation)
      IF (ALLOCATED(qs_env%mp2_env%ri_rpa_im_time%wkp_V)) THEN
         wkp_V(:) = qs_env%mp2_env%ri_rpa_im_time%wkp_V(:)
      ELSE
         wkp_V(:) = wkp(:)
      END IF

      IF (kpoint_weights_W_method == kp_weights_W_uniform) THEN

         !  in the k-point weights wkp, there might be k-point extrapolation included
         wkp_W(:) = wkp(:)

      ELSE IF (kpoint_weights_W_method == kp_weights_W_tailored .OR. &
               kpoint_weights_W_method == kp_weights_W_auto) THEN

         IF (kpoint_weights_W_method == kp_weights_W_tailored) &
            exp_kpoints = qs_env%mp2_env%ri_rpa_im_time%exp_tailored_weights

         IF (kpoint_weights_W_method == kp_weights_W_auto) THEN
            IF (SUM(periodic) == 2) exp_kpoints = -1.0_dp
         END IF

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

         wkp_W(1:nkp) = wkp_tmp(1:nkp)

         DEALLOCATE (matrix_lin_eqs, right_side, wkp_tmp)

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
! **************************************************************************************************
   SUBROUTINE cp_cfm_upper_to_full(cfm_mat_Q)

      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mat_Q

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_cfm_upper_to_full'

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_cfm_type)                                  :: cfm_mat_work

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_mat_work, cfm_mat_Q%matrix_struct)

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

      CALL cp_cfm_release(cfm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_upper_to_full

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param Eigenval_kp ...
! **************************************************************************************************
   SUBROUTINE get_bandstruc_and_k_dependent_MOs(qs_env, Eigenval_kp)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Eigenval_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_bandstruc_and_k_dependent_MOs'

      INTEGER                                            :: handle, ikp, ispin, nmo, nspins
      INTEGER, DIMENSION(3)                              :: nkp_grid_G
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ev
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: kpgeneral
      TYPE(kpoint_type), POINTER                         :: kpoints_Sigma
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (qs_env%mp2_env%ri_rpa_im_time%kpoints_G, &
               qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma, &
               qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma_no_xc, &
               para_env)

      nkp_grid_G(1:3) = (/1, 1, 1/)

      CALL get_qs_env(qs_env=qs_env, para_env=para_env)

      CALL create_kp_and_calc_kp_orbitals(qs_env, qs_env%mp2_env%ri_rpa_im_time%kpoints_G, &
                                          "MONKHORST-PACK", para_env%num_pe, &
                                          mp_grid=nkp_grid_G(1:3))

      IF (qs_env%mp2_env%ri_g0w0%do_kpoints_Sigma) THEN

         ! set up k-points for GW band structure calculation, will be completed later
         CALL get_kpgeneral_for_Sigma_kpoints(qs_env, kpgeneral)

         CALL create_kp_and_calc_kp_orbitals(qs_env, qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma, &
                                             "GENERAL", para_env%num_pe, &
                                             kpgeneral=kpgeneral)

         CALL create_kp_and_calc_kp_orbitals(qs_env, qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma_no_xc, &
                                             "GENERAL", para_env%num_pe, &
                                             kpgeneral=kpgeneral, with_xc_terms=.FALSE.)

         kpoints_Sigma => qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma
         nmo = SIZE(Eigenval_kp, 1)
         nspins = SIZE(Eigenval_kp, 3)

         ALLOCATE (qs_env%mp2_env%ri_rpa_im_time%Eigenval_Gamma(nmo))
         qs_env%mp2_env%ri_rpa_im_time%Eigenval_Gamma(:) = Eigenval_kp(:, 1, 1)

         DEALLOCATE (Eigenval_kp)

         ALLOCATE (Eigenval_kp(nmo, kpoints_Sigma%nkp, nspins))

         DO ikp = 1, kpoints_Sigma%nkp

            DO ispin = 1, nspins

               ev => kpoints_Sigma%kp_env(ikp)%kpoint_env%mos(1, ispin)%eigenvalues

               Eigenval_kp(:, ikp, ispin) = ev(:)

            END DO

         END DO

         DEALLOCATE (kpgeneral)

      END IF

      CALL release_hfx_stuff(qs_env)

      CALL timestop(handle)

   END SUBROUTINE get_bandstruc_and_k_dependent_MOs

! **************************************************************************************************
!> \brief releases part of the given qs_env in order to save memory
!> \param qs_env the object to release
! **************************************************************************************************
   SUBROUTINE release_hfx_stuff(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      IF (ASSOCIATED(qs_env%x_data) .AND. .NOT. qs_env%mp2_env%ri_g0w0%do_ri_Sigma_x) THEN
         CALL hfx_release(qs_env%x_data)
      END IF

   END SUBROUTINE release_hfx_stuff

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpoints ...
!> \param scheme ...
!> \param group_size_ext ...
!> \param mp_grid ...
!> \param kpgeneral ...
!> \param with_xc_terms ...
! **************************************************************************************************
   SUBROUTINE create_kp_and_calc_kp_orbitals(qs_env, kpoints, scheme, &
                                             group_size_ext, mp_grid, kpgeneral, with_xc_terms)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*), INTENT(IN)                       :: scheme
      INTEGER                                            :: group_size_ext
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: mp_grid
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: kpgeneral
      LOGICAL, OPTIONAL                                  :: with_xc_terms

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_kp_and_calc_kp_orbitals'
      COMPLEX(KIND=dp), PARAMETER :: cone = CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
         czero = CMPLX(0.0_dp, 0.0_dp, KIND=dp), ione = CMPLX(0.0_dp, 1.0_dp, KIND=dp)

      INTEGER                                            :: handle, i_dim, i_re_im, ikp, ispin, nkp, &
                                                            nspins
      INTEGER, DIMENSION(3)                              :: cell_grid, periodic
      LOGICAL                                            :: my_with_xc_terms
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: cksmat, cmos, csmat, cwork
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(cp_fm_type), POINTER                          :: imos, rmos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_desymm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_ks_kp, mat_s_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      my_with_xc_terms = .TRUE.
      IF (PRESENT(with_xc_terms)) my_with_xc_terms = with_xc_terms

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      matrix_s=matrix_s, &
                      scf_env=scf_env, &
                      scf_control=scf_control, &
                      cell=cell)

      ! get kpoints
      CALL calculate_kpoints_for_bs(kpoints, scheme, kpgeneral=kpgeneral, mp_grid=mp_grid, &
                                    group_size_ext=group_size_ext)

      CALL kpoint_env_initialize(kpoints, para_env, blacs_env)

      ! calculate all MOs that are accessible in the given
      ! Gaussian AO basis, therefore nadd=1E10
      CALL kpoint_initialize_mos(kpoints, qs_env%mos, 2000000000)
      CALL kpoint_initialize_mo_set(kpoints)

      CALL get_cell(cell=cell, periodic=periodic)

      DO i_dim = 1, 3
         ! we have at most 3 neigboring cells per dimension and at least one because
         ! the density response at Gamma is only divided to neighboring
         IF (periodic(i_dim) == 1) THEN
            cell_grid(i_dim) = MAX(MIN((kpoints%nkp_grid(i_dim)/2)*2 - 1, 1), 3)
         ELSE
            cell_grid(i_dim) = 1
         END IF
      END DO
      CALL init_cell_index_rpa(cell_grid, kpoints%cell_to_index, kpoints%index_to_cell, cell)

      ! get S(k)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, scf_env=scf_env, scf_control=scf_control, dft_control=dft_control)

      NULLIFY (matrix_s_desymm)
      CALL dbcsr_allocate_matrix_set(matrix_s_desymm, 1)
      ALLOCATE (matrix_s_desymm(1)%matrix)
      CALL dbcsr_create(matrix=matrix_s_desymm(1)%matrix, template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(matrix_s(1)%matrix, matrix_s_desymm(1)%matrix)

      CALL mat_kp_from_mat_gamma(qs_env, mat_s_kp, matrix_s_desymm(1)%matrix, kpoints, 1)

      CALL get_kpoint_info(kpoints, nkp=nkp)

      matrix_struct => kpoints%kp_env(1)%kpoint_env%wmat(1, 1)%matrix_struct

      CALL cp_cfm_create(cksmat, matrix_struct)
      CALL cp_cfm_create(csmat, matrix_struct)
      CALL cp_cfm_create(cmos, matrix_struct)
      CALL cp_cfm_create(cwork, matrix_struct)
      CALL cp_fm_create(fm_work, matrix_struct)

      nspins = dft_control%nspins

      DO ispin = 1, nspins

         ! get H(k)
         IF (my_with_xc_terms) THEN
            CALL mat_kp_from_mat_gamma(qs_env, mat_ks_kp, qs_env%mp2_env%ri_g0w0%matrix_ks(ispin)%matrix, kpoints, ispin)
         ELSE
            CALL mat_kp_from_mat_gamma(qs_env, mat_ks_kp, qs_env%mp2_env%ri_g0w0%matrix_sigma_x_minus_vxc(ispin)%matrix, &
                                       kpoints, ispin)
         END IF

         DO ikp = 1, nkp

            CALL copy_dbcsr_to_fm(mat_ks_kp(ikp, 1)%matrix, kpoints%kp_env(ikp)%kpoint_env%wmat(1, ispin))
            CALL cp_cfm_scale_and_add_fm(czero, cksmat, cone, kpoints%kp_env(ikp)%kpoint_env%wmat(1, ispin))

            CALL copy_dbcsr_to_fm(mat_ks_kp(ikp, 2)%matrix, kpoints%kp_env(ikp)%kpoint_env%wmat(2, ispin))
            CALL cp_cfm_scale_and_add_fm(cone, cksmat, ione, kpoints%kp_env(ikp)%kpoint_env%wmat(2, ispin))

            CALL copy_dbcsr_to_fm(mat_s_kp(ikp, 1)%matrix, fm_work)
            CALL cp_cfm_scale_and_add_fm(czero, csmat, cone, fm_work)

            CALL copy_dbcsr_to_fm(mat_s_kp(ikp, 2)%matrix, fm_work)
            CALL cp_cfm_scale_and_add_fm(cone, csmat, ione, fm_work)

            kp => kpoints%kp_env(ikp)%kpoint_env

            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=rmos, eigenvalues=eigenvalues)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=imos)

            IF (scf_env%cholesky_method == cholesky_off .OR. &
                qs_env%mp2_env%ri_rpa_im_time%make_overlap_mat_ao_pos_definite) THEN
               CALL cp_cfm_geeig_canon(cksmat, csmat, cmos, eigenvalues, cwork, scf_control%eps_eigval)
            ELSE
               CALL cp_cfm_geeig(cksmat, csmat, cmos, eigenvalues, cwork)
            END IF

            CALL cp_cfm_to_fm(cmos, rmos, imos)

            kp%mos(2, ispin)%eigenvalues = eigenvalues

         END DO

      END DO

      DO ikp = 1, nkp
         DO i_re_im = 1, 2
            CALL dbcsr_deallocate_matrix(mat_ks_kp(ikp, i_re_im)%matrix)
         END DO
      END DO
      DEALLOCATE (mat_ks_kp)

      DO ikp = 1, nkp
         DO i_re_im = 1, 2
            CALL dbcsr_deallocate_matrix(mat_s_kp(ikp, i_re_im)%matrix)
         END DO
      END DO
      DEALLOCATE (mat_s_kp)

      CALL dbcsr_deallocate_matrix(matrix_s_desymm(1)%matrix)
      DEALLOCATE (matrix_s_desymm)

      CALL cp_cfm_release(cksmat)
      CALL cp_cfm_release(csmat)
      CALL cp_cfm_release(cwork)
      CALL cp_cfm_release(cmos)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE create_kp_and_calc_kp_orbitals

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mat_kp ...
!> \param mat_gamma ...
!> \param kpoints ...
!> \param ispin ...
!> \param real_mat_real_space ...
! **************************************************************************************************
   SUBROUTINE mat_kp_from_mat_gamma(qs_env, mat_kp, mat_gamma, kpoints, ispin, real_mat_real_space)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_kp
      TYPE(dbcsr_type)                                   :: mat_gamma
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ispin
      LOGICAL, INTENT(IN), OPTIONAL                      :: real_mat_real_space

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mat_kp_from_mat_gamma'

      INTEGER                                            :: handle, i_cell, i_re_im, ikp, nkp, &
                                                            num_cells
      INTEGER, DIMENSION(3)                              :: periodic
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_real_space

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, cell=cell)
      CALL get_cell(cell=cell, periodic=periodic)
      num_cells = 3**(periodic(1) + periodic(2) + periodic(3))

      NULLIFY (mat_real_space)
      CALL dbcsr_allocate_matrix_set(mat_real_space, num_cells)
      DO i_cell = 1, num_cells
         ALLOCATE (mat_real_space(i_cell)%matrix)
         CALL dbcsr_create(matrix=mat_real_space(i_cell)%matrix, &
                           template=mat_gamma)
         CALL dbcsr_reserve_all_blocks(mat_real_space(i_cell)%matrix)
         CALL dbcsr_set(mat_real_space(i_cell)%matrix, 0.0_dp)
      END DO

      CALL dbcsr_copy(mat_real_space(1)%matrix, mat_gamma)

      CALL get_mat_cell_T_from_mat_gamma(mat_real_space, qs_env, kpoints, 2, 0)

      NULLIFY (xkp, cell_to_index)
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, cell_to_index=cell_to_index)

      IF (ispin == 1) THEN
         NULLIFY (mat_kp)
         CALL dbcsr_allocate_matrix_set(mat_kp, nkp, 2)
         DO ikp = 1, nkp
            DO i_re_im = 1, 2
               ALLOCATE (mat_kp(ikp, i_re_im)%matrix)
               CALL dbcsr_create(matrix=mat_kp(ikp, i_re_im)%matrix, template=mat_gamma)
               CALL dbcsr_reserve_all_blocks(mat_kp(ikp, i_re_im)%matrix)
               CALL dbcsr_set(mat_kp(ikp, i_re_im)%matrix, 0.0_dp)
            END DO
         END DO
      END IF

      IF (PRESENT(real_mat_real_space)) THEN
         CALL real_space_to_kpoint_transform_rpa(mat_kp(:, 1), mat_kp(:, 2), mat_real_space, kpoints, 0.0_dp, &
                                                 real_mat_real_space)
      ELSE
         CALL real_space_to_kpoint_transform_rpa(mat_kp(:, 1), mat_kp(:, 2), mat_real_space, kpoints, 0.0_dp)
      END IF

      DO i_cell = 1, num_cells
         CALL dbcsr_deallocate_matrix(mat_real_space(i_cell)%matrix)
      END DO
      DEALLOCATE (mat_real_space)

      CALL timestop(handle)

   END SUBROUTINE mat_kp_from_mat_gamma

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpgeneral ...
! **************************************************************************************************
   SUBROUTINE get_kpgeneral_for_Sigma_kpoints(qs_env, kpgeneral)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: kpgeneral

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kpgeneral_for_Sigma_kpoints'

      INTEGER                                            :: handle, i_kp_in_kp_line, i_special_kp, &
                                                            i_x, ikk, j_y, k_z, n_kp_in_kp_line, &
                                                            n_special_kp
      INTEGER, DIMENSION(:), POINTER                     :: nkp_grid

      CALL timeset(routineN, handle)

      n_special_kp = qs_env%mp2_env%ri_g0w0%n_special_kp
      n_kp_in_kp_line = qs_env%mp2_env%ri_g0w0%n_kp_in_kp_line
      IF (n_special_kp > 0) THEN
         qs_env%mp2_env%ri_g0w0%nkp_self_energy_special_kp = n_kp_in_kp_line*(n_special_kp - 1) + 1
      ELSE
         qs_env%mp2_env%ri_g0w0%nkp_self_energy_special_kp = 0
      END IF

      qs_env%mp2_env%ri_g0w0%nkp_self_energy_monkh_pack = qs_env%mp2_env%ri_g0w0%kp_grid_Sigma(1)* &
                                                          qs_env%mp2_env%ri_g0w0%kp_grid_Sigma(2)* &
                                                          qs_env%mp2_env%ri_g0w0%kp_grid_Sigma(3)

      qs_env%mp2_env%ri_g0w0%nkp_self_energy = qs_env%mp2_env%ri_g0w0%nkp_self_energy_special_kp + &
                                               qs_env%mp2_env%ri_g0w0%nkp_self_energy_monkh_pack

      ALLOCATE (kpgeneral(3, qs_env%mp2_env%ri_g0w0%nkp_self_energy))

      IF (n_special_kp > 0) THEN

         kpgeneral(1:3, 1) = qs_env%mp2_env%ri_g0w0%xkp_special_kp(1:3, 1)

         ikk = 1

         DO i_special_kp = 2, n_special_kp
            DO i_kp_in_kp_line = 1, n_kp_in_kp_line

               ikk = ikk + 1
               kpgeneral(1:3, ikk) = qs_env%mp2_env%ri_g0w0%xkp_special_kp(1:3, i_special_kp - 1) + &
                                     REAL(i_kp_in_kp_line, KIND=dp)/REAL(n_kp_in_kp_line, KIND=dp)* &
                                     (qs_env%mp2_env%ri_g0w0%xkp_special_kp(1:3, i_special_kp) - &
                                      qs_env%mp2_env%ri_g0w0%xkp_special_kp(1:3, i_special_kp - 1))

            END DO
         END DO

      ELSE

         ikk = 0

      END IF

      nkp_grid => qs_env%mp2_env%ri_g0w0%kp_grid_Sigma

      DO i_x = 1, nkp_grid(1)
         DO j_y = 1, nkp_grid(2)
            DO k_z = 1, nkp_grid(3)
               ikk = ikk + 1
               kpgeneral(1, ikk) = REAL(2*i_x - nkp_grid(1) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(1), KIND=dp))
               kpgeneral(2, ikk) = REAL(2*j_y - nkp_grid(2) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(2), KIND=dp))
               kpgeneral(3, ikk) = REAL(2*k_z - nkp_grid(3) - 1, KIND=dp)/(2._dp*REAL(nkp_grid(3), KIND=dp))
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_kpgeneral_for_Sigma_kpoints

END MODULE rpa_gw_kpoints_util
