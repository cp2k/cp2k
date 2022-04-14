!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines treating RPA calculations with kpoints
!> \par History
!>      06.2019 Split from rpa_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE rpa_kpoints
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_column_scale,&
                                              cp_cfm_gemm,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_transpose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_filter, dbcsr_get_block_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_set, dbcsr_transposed, &
        dbcsr_type, dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              pi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_sum,&
                                              mp_sync
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE rpa_gw_kpoints,                  ONLY: compute_Wc_kp_tau_GW,&
                                              compute_Wc_real_space_tau_GW
   USE rpa_im_time,                     ONLY: init_cell_index_rpa
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_kpoints'

   PUBLIC :: RPA_postprocessing_kp, cp_cfm_robust_cholesky, real_space_to_kpoint_transform_rpa, get_mat_cell_T_from_mat_gamma

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dimen_RI ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param nkp ...
!> \param count_ev_sc_GW ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param Erpa ...
!> \param tau_tj ...
!> \param tj ...
!> \param wj ...
!> \param weights_cos_tf_w_to_t ...
!> \param wkp_W ...
!> \param wkp_V ...
!> \param do_gw_im_time ...
!> \param do_ri_Sigma_x ...
!> \param do_kpoints_from_Gamma ...
!> \param do_kpoints_cubic_RPA ...
!> \param cfm_mat_W_kp_tau ...
!> \param cfm_mat_Q ...
!> \param ikp_local ...
!> \param mat_P_omega ...
!> \param mat_P_omega_kp ...
!> \param qs_env ...
!> \param eps_filter_im_time ...
!> \param unit_nr ...
!> \param kpoints ...
!> \param fm_mat_L ...
!> \param fm_mat_W ...
!> \param fm_mat_RI_global_work ...
!> \param mat_SinvVSinv ...
! **************************************************************************************************
   SUBROUTINE RPA_postprocessing_kp(dimen_RI, num_integ_points, jquad, nkp, count_ev_sc_GW, para_env, para_env_RPA, &
                                    Erpa, tau_tj, tj, wj, weights_cos_tf_w_to_t, wkp_W, wkp_V, do_gw_im_time, &
                                    do_ri_Sigma_x, do_kpoints_from_Gamma, do_kpoints_cubic_RPA, &
                                    cfm_mat_W_kp_tau, cfm_mat_Q, ikp_local, mat_P_omega, mat_P_omega_kp, &
                                    qs_env, eps_filter_im_time, unit_nr, kpoints, fm_mat_L, fm_mat_W, &
                                    fm_mat_RI_global_work, mat_SinvVSinv, fm_mat_SiVtrSi)

      INTEGER, INTENT(IN)                                :: dimen_RI, num_integ_points, jquad, nkp, &
                                                            count_ev_sc_GW
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_RPA
      REAL(KIND=dp), INTENT(INOUT)                       :: Erpa
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tj, wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_w_to_t
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wkp_W, wkp_V
      LOGICAL, INTENT(IN)                                :: do_gw_im_time, do_ri_Sigma_x, &
                                                            do_kpoints_from_Gamma, &
                                                            do_kpoints_cubic_RPA
      TYPE(cp_cfm_p_type), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(INOUT)                  :: cfm_mat_W_kp_tau
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ikp_local
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: mat_P_omega, mat_P_omega_kp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L, fm_mat_SiVtrSi
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_W
      TYPE(cp_fm_type), POINTER                          :: fm_mat_RI_global_work
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_SinvVSinv

      CHARACTER(LEN=*), PARAMETER :: routineN = 'RPA_postprocessing_kp'

      INTEGER                                            :: handle, ikp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj_dummy, tj_dummy, trace_Qomega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_cos_tf_w_to_t_dummy

      CALL timeset(routineN, handle)

      IF (do_kpoints_cubic_RPA .AND. do_gw_im_time) THEN
         CALL allocate_Wc_kp_tau_GW(cfm_mat_W_kp_tau, cfm_mat_Q, num_integ_points, jquad, &
                                    nkp, ikp_local)
      END IF

      IF (do_kpoints_from_Gamma) THEN
         CALL get_mat_cell_T_from_mat_gamma(mat_P_omega(jquad,:), qs_env, kpoints, jquad, unit_nr)
      END IF

      CALL transform_P_from_real_space_to_kpoints(mat_P_omega, mat_P_omega_kp, &
                                                  kpoints, eps_filter_im_time, jquad)

      ALLOCATE (trace_Qomega(dimen_RI))

      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T69,I9,A,I2)') &
         'GW_INFO| Computing chi and W for frequency point: ', jquad, '/', num_integ_points

      DO ikp = 1, nkp

         ! parallization, we either have all kpoints on all processors or a single kpoint per group
         IF (ikp_local(ikp) .NE. ikp) CYCLE

         ! 1. remove all spurious negative eigenvalues from P(iw,k), multiplication Q(iw,k) = K^H(k)P(iw,k)K(k)
         CALL compute_Q_kp_RPA(cfm_mat_Q, &
                               mat_P_omega_kp, &
                               fm_mat_L(ikp, 1)%matrix, &
                               fm_mat_L(ikp, 2)%matrix, &
                               fm_mat_RI_global_work, &
                               dimen_RI, ikp, nkp, ikp_local, para_env, &
                               qs_env%mp2_env%ri_rpa_im_time%make_chi_pos_definite)

         ! 2. Cholesky decomposition of Id + Q(iw,k)
         CALL cholesky_decomp_Q(cfm_mat_Q, para_env_RPA, trace_Qomega, dimen_RI)

         ! 3. Computing E_c^RPA = E_c^RPA + a_w/N_k*sum_k ln[det(1+Q(iw,k))-Tr(Q(iw,k))]
         CALL frequency_and_kpoint_integration(Erpa, cfm_mat_Q, para_env_RPA, trace_Qomega, &
                                               dimen_RI, wj(jquad), kpoints%wkp(ikp))

         IF (do_gw_im_time) THEN

            ! compute S^-1*V*S^-1 for exchange part of the self-energy in real space as W in real space
            IF (do_ri_Sigma_x .AND. jquad == 1 .AND. count_ev_sc_GW == 1 &
                .AND. do_kpoints_from_Gamma) THEN

               CALL get_dummys(tj_dummy, tau_tj_dummy, weights_cos_tf_w_to_t_dummy)

               CALL compute_Wc_real_space_tau_GW(fm_mat_W, cfm_mat_Q, &
                                                 fm_mat_L(ikp, 1)%matrix, &
                                                 fm_mat_L(ikp, 2)%matrix, &
                                                 dimen_RI, 1, 1, &
                                                 ikp, tj_dummy, tau_tj_dummy, weights_cos_tf_w_to_t_dummy, &
                                                 ikp_local, para_env, kpoints, qs_env, wkp_V, &
                                                 mat_SinvVSinv, do_W_and_not_V=.FALSE., &
                                                 fm_mat_SiVtrSi_re=fm_mat_SiVtrSi(ikp, 1)%matrix, &
                                                 fm_mat_SiVtrSi_im=fm_mat_SiVtrSi(ikp, 2)%matrix)

               CALL release_dummys(tj_dummy, tau_tj_dummy, weights_cos_tf_w_to_t_dummy)

            END IF
            IF (do_kpoints_from_Gamma) THEN

               CALL compute_Wc_real_space_tau_GW(fm_mat_W, cfm_mat_Q, &
                                                 fm_mat_L(ikp, 1)%matrix, &
                                                 fm_mat_L(ikp, 2)%matrix, &
                                                 dimen_RI, num_integ_points, jquad, &
                                                 ikp, tj, tau_tj, weights_cos_tf_w_to_t, &
                                                 ikp_local, para_env, kpoints, qs_env, wkp_W, &
                                                 mat_SinvVSinv, do_W_and_not_V=.TRUE.)

            END IF

            IF (do_kpoints_cubic_RPA) THEN

               CALL compute_Wc_kp_tau_GW(cfm_mat_W_kp_tau, cfm_mat_Q, &
                                         fm_mat_L(ikp, 1)%matrix, &
                                         fm_mat_L(ikp, 2)%matrix, &
                                         dimen_RI, num_integ_points, jquad, &
                                         ikp, tj, tau_tj, &
                                         weights_cos_tf_w_to_t)

            END IF
         END IF
      END DO

      DEALLOCATE (trace_Qomega)

      CALL timestop(handle)

   END SUBROUTINE RPA_postprocessing_kp

! **************************************************************************************************
!> \brief Decomposition A = L^H*L of a Hermitian matrix A
!> \param matrix ...
!> \param work ...
!> \param threshold ...
!> \param exponent ...
!> \param min_eigval ...
!> \author Jan Wilhelm, 07.2021
! **************************************************************************************************
   SUBROUTINE cp_cfm_robust_cholesky(matrix, work, threshold, exponent, min_eigval)
      TYPE(cp_cfm_type), POINTER                         :: matrix, work
      REAL(KIND=dp)                                      :: threshold, exponent
      REAL(KIND=dp), OPTIONAL                            :: min_eigval

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_cfm_robust_cholesky'
      COMPLEX(KIND=dp), PARAMETER :: czero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: eigenvalues_exponent
      INTEGER                                            :: handle, i, ncol_global, nrow_global
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues

      CALL timeset(routineN, handle)

      ! Test that matrix is square
      CALL cp_cfm_get_info(matrix, nrow_global=nrow_global, ncol_global=ncol_global)
      CPASSERT(nrow_global == ncol_global)
      ALLOCATE (eigenvalues(nrow_global))
      eigenvalues(:) = 0.0_dp
      ALLOCATE (eigenvalues_exponent(nrow_global))
      eigenvalues_exponent(:) = czero

      ! Diagonalize matrix: get eigenvectors and eigenvalues
      CALL cp_cfm_heevd(matrix, work, eigenvalues)

      DO i = 1, nrow_global
         IF (eigenvalues(i) > threshold) THEN
            eigenvalues_exponent(i) = CMPLX((eigenvalues(i))**exponent, threshold, KIND=dp)
         ELSE
            IF (PRESENT(min_eigval)) THEN
               eigenvalues_exponent(i) = CMPLX(min_eigval, 0.0_dp, KIND=dp)
            ELSE
               eigenvalues_exponent(i) = czero
            END IF
         END IF
      END DO

      CALL cp_cfm_column_scale(work, eigenvalues_exponent)

      CALL cp_cfm_transpose(work, 'C', matrix)

      DEALLOCATE (eigenvalues, eigenvalues_exponent)

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_robust_cholesky

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

      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: mat_P_omega_kp
      TYPE(cp_fm_type), POINTER                          :: fm_mat_L_re, fm_mat_L_im, &
                                                            fm_mat_RI_global_work
      INTEGER, INTENT(IN)                                :: dimen_RI, ikp, nkp
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ikp_local
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: make_chi_pos_definite

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Q_kp_RPA'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_L, cfm_mat_work
      TYPE(cp_fm_type), POINTER                          :: fm_mat_work

      CALL timeset(routineN, handle)

      NULLIFY (cfm_mat_work)
      CALL cp_cfm_create(cfm_mat_work, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      NULLIFY (cfm_mat_L)
      CALL cp_cfm_create(cfm_mat_L, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_L, z_zero)

      NULLIFY (fm_mat_work)
      CALL cp_fm_create(fm_mat_work, fm_mat_L_re%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

      ! 1. Convert the dbcsr matrix mat_P_omega_kp (that is chi(k,iw)) to a full matrix and
      !    distribute it to subgroups
      CALL mat_P_to_subgroup(mat_P_omega_kp, fm_mat_RI_global_work, &
                             fm_mat_work, cfm_mat_Q, ikp, nkp, ikp_local, para_env)

      ! 2. Remove all negative eigenvalues from chi(k,iw)
      IF (make_chi_pos_definite) THEN
         CALL cp_cfm_robust_cholesky(cfm_mat_Q, cfm_mat_work, threshold=0.0_dp, exponent=0.5_dp)
         CALL cp_cfm_gemm("C", "N", dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_Q, &
                          z_zero, cfm_mat_work)
         CALL cp_cfm_to_cfm(cfm_mat_work, cfm_mat_Q)

      END IF

      ! 3. Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      ! 4. work = P(iw,k)*L(k)
      CALL cp_cfm_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                       z_zero, cfm_mat_work)

      ! 5. Q(iw,k) = L^H(k)*work
      CALL cp_cfm_gemm('C', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_work, &
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
      TYPE(cp_fm_type), POINTER                          :: fm_mat_RI_global_work, fm_mat_work
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      INTEGER, INTENT(IN)                                :: ikp, nkp
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ikp_local
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mat_P_to_subgroup'

      INTEGER                                            :: handle, jkp
      TYPE(cp_fm_type), POINTER                          :: fm_dummy
      TYPE(dbcsr_type), POINTER                          :: mat_P_omega_im, mat_P_omega_re

      CALL timeset(routineN, handle)

      IF (SUM(ikp_local) > nkp) THEN

         mat_P_omega_re => mat_P_omega_kp(1, ikp)%matrix
         CALL cp_fm_set_all(fm_mat_work, 0.0_dp)
         CALL copy_dbcsr_to_fm(mat_P_omega_re, fm_mat_work)
         CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_Q, z_one, fm_mat_work)

         mat_P_omega_im => mat_P_omega_kp(2, ikp)%matrix
         CALL cp_fm_set_all(fm_mat_work, 0.0_dp)
         CALL copy_dbcsr_to_fm(mat_P_omega_im, fm_mat_work)
         CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_Q, gaussi, fm_mat_work)

      ELSE

         DO jkp = 1, nkp

            mat_P_omega_re => mat_P_omega_kp(1, jkp)%matrix

            CALL cp_fm_set_all(fm_mat_RI_global_work, 0.0_dp)
            CALL copy_dbcsr_to_fm(mat_P_omega_re, fm_mat_RI_global_work)

            IF (ANY(ikp_local(:) == jkp)) THEN
               CALL cp_fm_set_all(fm_mat_work, 0.0_dp)
               CALL cp_fm_copy_general(fm_mat_RI_global_work, fm_mat_work, para_env)
               CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_Q, z_one, fm_mat_work)
            ELSE
               NULLIFY (fm_dummy)
               CALL cp_fm_copy_general(fm_mat_RI_global_work, fm_dummy, para_env)
            END IF

            mat_P_omega_im => mat_P_omega_kp(2, jkp)%matrix

            CALL cp_fm_set_all(fm_mat_RI_global_work, 0.0_dp)
            CALL copy_dbcsr_to_fm(mat_P_omega_im, fm_mat_RI_global_work)

            IF (ANY(ikp_local(:) == jkp)) THEN
               CALL cp_fm_set_all(fm_mat_work, 0.0_dp)
               CALL cp_fm_copy_general(fm_mat_RI_global_work, fm_mat_work, para_env)
               CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_Q, gaussi, fm_mat_work)
            ELSE
               NULLIFY (fm_dummy)
               CALL cp_fm_copy_general(fm_mat_RI_global_work, fm_dummy, para_env)
            END IF

         END DO

      END IF

      CALL mp_sync(para_env%group)

      CALL timestop(handle)

   END SUBROUTINE mat_P_to_subgroup

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
!> \param para_env_RPA ...
!> \param trace_Qomega ...
!> \param dimen_RI ...
! **************************************************************************************************
   SUBROUTINE cholesky_decomp_Q(cfm_mat_Q, para_env_RPA, trace_Qomega, dimen_RI)

      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_RPA
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: trace_Qomega
      INTEGER, INTENT(IN)                                :: dimen_RI

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cholesky_decomp_Q'

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q_tmp, cfm_mat_work

      CALL timeset(routineN, handle)

      NULLIFY (cfm_mat_work)
      CALL cp_cfm_create(cfm_mat_work, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      NULLIFY (cfm_mat_Q_tmp)
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
      CALL mp_sum(trace_Qomega, para_env_RPA%group)

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
!> \param para_env_RPA ...
!> \param trace_Qomega ...
!> \param dimen_RI ...
!> \param freq_weight ...
!> \param kp_weight ...
! **************************************************************************************************
   SUBROUTINE frequency_and_kpoint_integration(Erpa, cfm_mat_Q, para_env_RPA, trace_Qomega, &
                                               dimen_RI, freq_weight, kp_weight)

      REAL(KIND=dp), INTENT(INOUT)                       :: Erpa
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_RPA
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
      CALL mp_sum(Q_log, para_env_RPA%group)

      FComega = 0.0_dp
      DO iiB = 1, dimen_RI
         IF (MODULO(iiB, para_env_RPA%num_pe) /= para_env_RPA%mepos) CYCLE
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
!> \param cfm_mat_W_kp_tau ...
!> \param cfm_mat_Q ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param nkp ...
!> \param ikp_local ...
! **************************************************************************************************
   SUBROUTINE allocate_Wc_kp_tau_GW(cfm_mat_W_kp_tau, cfm_mat_Q, num_integ_points, jquad, &
                                    nkp, ikp_local)

      TYPE(cp_cfm_p_type), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: cfm_mat_W_kp_tau
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      INTEGER, INTENT(IN)                                :: num_integ_points, jquad, nkp
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ikp_local

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_Wc_kp_tau_GW'

      INTEGER                                            :: handle, iquad, jkp

      CALL timeset(routineN, handle)

      ! Fourier transform from w to t; allocate only one (at first time point and also at
      ! first kpoint in case we do not do
      IF (jquad == 1) THEN
         ALLOCATE (cfm_mat_W_kp_tau(nkp, num_integ_points))
         DO iquad = 1, num_integ_points
            DO jkp = 1, nkp
               NULLIFY (cfm_mat_W_kp_tau(jkp, iquad)%matrix)
               IF (.NOT. (ANY(ikp_local(:) == jkp))) CYCLE
               CALL cp_cfm_create(cfm_mat_W_kp_tau(jkp, iquad)%matrix, cfm_mat_Q%matrix_struct)
               CALL cp_cfm_set_all(matrix=cfm_mat_W_kp_tau(jkp, iquad)%matrix, alpha=z_zero)
            END DO
         END DO
      END IF

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
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)    :: mat_P_omega
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: jquad, unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_cell_T_from_mat_gamma'

      INTEGER                                            :: col, handle, i_cell, i_dim, num_cells_P, &
                                                            num_integ_points, row
      INTEGER, DIMENSION(3)                              :: cell_grid_P, periodic
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_P
      REAL(KIND=dp)                                      :: abs_rab_cell_i, first_cutoff, &
                                                            second_cutoff, weight
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector, rab_cell_i
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

      first_cutoff = qs_env%mp2_env%ri_rpa_im_time%abs_cutoffs_chi_W(1)
      second_cutoff = qs_env%mp2_env%ri_rpa_im_time%abs_cutoffs_chi_W(2)

      IF (jquad == 1 .AND. unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A,T72,F9.2)') 'GW_INFO| Relative cutoff 1: ', &
            qs_env%mp2_env%ri_rpa_im_time%rel_cutoffs_chi_W(1)
         WRITE (unit_nr, '(T3,A,T72,F9.2)') 'GW_INFO| Relative cutoff 2: ', &
            qs_env%mp2_env%ri_rpa_im_time%rel_cutoffs_chi_W(2)
         WRITE (unit_nr, '(T3,A,T72,F9.2)') 'GW_INFO| Cutoff radius 1 (Angstroem):', &
            cp_unit_from_cp2k(first_cutoff, "angstrom")
         WRITE (unit_nr, '(T3,A,T72,F9.2)') 'GW_INFO| Cutoff radius 2 (Angstroem):', &
            cp_unit_from_cp2k(second_cutoff, "angstrom")
         WRITE (unit_nr, '(T3,A,T72,F9.2)') &
            'GW_INFO| If Cutoff radius 1 = Cutoff radius 2 = 0.5: minimum image convention.'
         WRITE (unit_nr, '(T3,A,T66,ES15.2)') 'GW_INFO| RI regularization parameter: ', &
            qs_env%mp2_env%ri_rpa_im_time%regularization_RI
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

            ! smooth function that is only one if atoms are very close and becomes smoothly zero
            ! for larger distances
            IF (abs_rab_cell_i < first_cutoff) THEN
               weight = 1.0_dp
            ELSE IF (abs_rab_cell_i < second_cutoff) THEN
               weight = 0.5_dp*(COS(pi*(abs_rab_cell_i - first_cutoff)/ &
                                    (second_cutoff - first_cutoff)) + 1.0_dp)
            ELSE
               weight = 0.0_dp
            END IF

            data_block(:, :) = data_block(:, :)*weight

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
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm

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
!> \param index_to_cell ...
!> \param eps_filter_im_time ...
! **************************************************************************************************
   SUBROUTINE real_space_to_kpoint_transform_rpa(real_mat_kp, imag_mat_kp, mat_real_space, &
                                                 kpoints, eps_filter_im_time)

      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: real_mat_kp, imag_mat_kp, mat_real_space
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'real_space_to_kpoint_transform_rpa'

      INTEGER                                            :: handle, i_cell, ik, nkp, num_cells
      INTEGER, DIMENSION(3)                              :: cell
      REAL(KIND=dp)                                      :: arg, coskl, sinkl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(dbcsr_type)                                   :: mat_work

      CALL timeset(routineN, handle)

      CALL dbcsr_create(matrix=mat_work, &
                        template=real_mat_kp(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_reserve_all_blocks(mat_work)
      CALL dbcsr_set(mat_work, 0.0_dp)

      ! this kpoint environme t should be the kpoints for D(it) and X(it) created in init_cell_index_rpa
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp)

      NULLIFY(index_to_cell)
      index_to_cell => kpoints%index_to_cell

      num_cells = SIZE(index_to_cell, 2)

      CPASSERT(SIZE(mat_real_space) >= num_cells/2 + 1)

      DO ik = 1, nkp

         CALL dbcsr_set(real_mat_kp(ik)%matrix, 0.0_dp)
         CALL dbcsr_set(imag_mat_kp(ik)%matrix, 0.0_dp)

         CALL dbcsr_reserve_all_blocks(real_mat_kp(ik)%matrix)
         CALL dbcsr_reserve_all_blocks(imag_mat_kp(ik)%matrix)

         DO i_cell = 1, num_cells/2 + 1

            cell(:) = index_to_cell(:, i_cell)

            arg = REAL(cell(1), dp)*xkp(1, ik) + REAL(cell(2), dp)*xkp(2, ik) + REAL(cell(3), dp)*xkp(3, ik)
            coskl = COS(twopi*arg)
            sinkl = SIN(twopi*arg)

            CALL dbcsr_add_local(real_mat_kp(ik)%matrix, mat_real_space(i_cell)%matrix, 1.0_dp, coskl)
            CALL dbcsr_add_local(imag_mat_kp(ik)%matrix, mat_real_space(i_cell)%matrix, 1.0_dp, sinkl)

            IF (.NOT. (cell(1) == 0 .AND. cell(2) == 0 .AND. cell(3) == 0)) THEN

               CALL dbcsr_transposed(mat_work, mat_real_space(i_cell)%matrix)

               CALL dbcsr_add_local(real_mat_kp(ik)%matrix, mat_work, 1.0_dp, coskl)
               CALL dbcsr_add_local(imag_mat_kp(ik)%matrix, mat_work, 1.0_dp, -sinkl)

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

END MODULE rpa_kpoints
