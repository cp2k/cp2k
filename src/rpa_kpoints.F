!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines treating RPA calculations with kpoints
!> \par History
!>      06.2019 Split from rpa_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE rpa_kpoints
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_gemm,&
                                              cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_sum
   USE qs_environment_types,            ONLY: qs_environment_type
   USE rpa_gw_kpoints,                  ONLY: compute_Wc_kp_tau_GW,&
                                              compute_Wc_real_space_tau_GW
   USE rpa_im_time,                     ONLY: get_P_cell_T_from_P_gamma,&
                                              transform_P_from_real_space_to_kpoints
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_kpoints'

   PUBLIC :: RPA_postprocessing_kp

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
!> \param kpoints ...
!> \param fm_mat_L ...
!> \param fm_mat_W_tau ...
!> \param fm_mat_RI_global_work ...
!> \param mat_SinvVSinv ...
! **************************************************************************************************
   SUBROUTINE RPA_postprocessing_kp(dimen_RI, num_integ_points, jquad, nkp, count_ev_sc_GW, para_env, para_env_RPA, &
                                    Erpa, tau_tj, tj, wj, weights_cos_tf_w_to_t, wkp_W, do_gw_im_time, do_ri_Sigma_x, &
                                    do_kpoints_from_Gamma, do_kpoints_cubic_RPA, cfm_mat_W_kp_tau, cfm_mat_Q, ikp_local, &
                                    mat_P_omega, mat_P_omega_kp, qs_env, eps_filter_im_time, kpoints, fm_mat_L, fm_mat_W_tau, &
                                    fm_mat_RI_global_work, mat_SinvVSinv)

      INTEGER, INTENT(IN)                                :: dimen_RI, num_integ_points, jquad, nkp, &
                                                            count_ev_sc_GW
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      REAL(KIND=dp), INTENT(INOUT)                       :: Erpa
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tj, wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_w_to_t
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp_W
      LOGICAL, INTENT(IN)                                :: do_gw_im_time, do_ri_Sigma_x, &
                                                            do_kpoints_from_Gamma, &
                                                            do_kpoints_cubic_RPA
      TYPE(cp_cfm_p_type), DIMENSION(:, :), POINTER      :: cfm_mat_W_kp_tau
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ikp_local
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega, mat_P_omega_kp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_W_tau
      TYPE(cp_fm_type), POINTER                          :: fm_mat_RI_global_work
      TYPE(dbcsr_p_type)                                 :: mat_SinvVSinv

      CHARACTER(LEN=*), PARAMETER :: routineN = 'RPA_postprocessing_kp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj_dummy, tj_dummy, trace_Qomega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_cos_tf_w_to_t_dummy

      CALL timeset(routineN, handle)

      IF (do_kpoints_cubic_RPA .AND. do_gw_im_time) THEN
         CALL allocate_Wc_kp_tau_GW(cfm_mat_W_kp_tau, cfm_mat_Q, num_integ_points, jquad, &
                                    nkp, ikp_local)
      END IF

      IF (do_kpoints_from_Gamma) THEN
         CALL get_P_cell_T_from_P_gamma(mat_P_omega, qs_env, kpoints, jquad)
      END IF

      CALL transform_P_from_real_space_to_kpoints(mat_P_omega, mat_P_omega_kp, &
                                                  kpoints, eps_filter_im_time, jquad)

      ALLOCATE (trace_Qomega(dimen_RI))

      DO ikp = 1, nkp

         ! parallization
         IF (ikp_local(ikp) .NE. ikp) CYCLE

         ! 1. multiplication Q(iw,k) = K^H(k)P(iw,k)K(k)
         CALL compute_Q_kp_RPA(cfm_mat_Q, &
                               mat_P_omega_kp, &
                               fm_mat_L(ikp, 1)%matrix, &
                               fm_mat_L(ikp, 2)%matrix, &
                               fm_mat_RI_global_work, &
                               dimen_RI, ikp, nkp, ikp_local, para_env)

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

               CALL compute_Wc_real_space_tau_GW(fm_mat_W_tau, cfm_mat_Q, &
                                                 fm_mat_L(ikp, 1)%matrix, &
                                                 fm_mat_L(ikp, 2)%matrix, &
                                                 dimen_RI, 1, 1, &
                                                 ikp, tj_dummy, tau_tj_dummy, weights_cos_tf_w_to_t_dummy, &
                                                 ikp_local, para_env, kpoints, qs_env, wkp_W, &
                                                 mat_SinvVSinv, do_W_and_not_V=.FALSE.)

               CALL release_dummys(tj_dummy, tau_tj_dummy, weights_cos_tf_w_to_t_dummy)

            END IF
            IF (do_kpoints_from_Gamma) THEN

               CALL compute_Wc_real_space_tau_GW(fm_mat_W_tau, cfm_mat_Q, &
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
! **************************************************************************************************
   SUBROUTINE compute_Q_kp_RPA(cfm_mat_Q, mat_P_omega_kp, fm_mat_L_re, fm_mat_L_im, &
                               fm_mat_RI_global_work, dimen_RI, ikp, nkp, ikp_local, para_env)

      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega_kp
      TYPE(cp_fm_type), POINTER                          :: fm_mat_L_re, fm_mat_L_im, &
                                                            fm_mat_RI_global_work
      INTEGER, INTENT(IN)                                :: dimen_RI, ikp, nkp
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: ikp_local
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Q_kp_RPA', &
         routineP = moduleN//':'//routineN

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

      CALL mat_P_to_subgroup(mat_P_omega_kp, fm_mat_RI_global_work, &
                             fm_mat_work, cfm_mat_Q, ikp, nkp, ikp_local, para_env)

      ! 2. Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      ! 3. work = P(iw,k)*L(k)
      CALL cp_cfm_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                       z_zero, cfm_mat_work)

      ! 4. Q(iw,k) = L^H(k)*work
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

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega_kp
      TYPE(cp_fm_type), POINTER                          :: fm_mat_RI_global_work, fm_mat_work
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      INTEGER, INTENT(IN)                                :: ikp, nkp
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: ikp_local
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mat_P_to_subgroup', &
         routineP = moduleN//':'//routineN

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
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: trace_Qomega
      INTEGER, INTENT(IN)                                :: dimen_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cholesky_decomp_Q', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

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
               cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB)+z_one
            END IF
         END DO
      END DO
      CALL mp_sum(trace_Qomega, para_env_RPA%group)

      CALL cp_cfm_cholesky_decompose(matrix=cfm_mat_Q, n=dimen_RI, info_out=info_chol)
      CPASSERT(info_chol == 0)

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
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: trace_Qomega
      INTEGER, INTENT(IN)                                :: dimen_RI
      REAL(KIND=dp), INTENT(IN)                          :: freq_weight, kp_weight

      CHARACTER(LEN=*), PARAMETER :: routineN = 'frequency_and_kpoint_integration', &
         routineP = moduleN//':'//routineN

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
         FComega = FComega+(Q_log(iiB)-trace_Qomega(iiB))/2.0_dp
      END DO

      Erpa = Erpa+FComega*freq_weight*kp_weight

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_dummys', routineP = moduleN//':'//routineN

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_dummys', routineP = moduleN//':'//routineN

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

      TYPE(cp_cfm_p_type), DIMENSION(:, :), POINTER      :: cfm_mat_W_kp_tau
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      INTEGER, INTENT(IN)                                :: num_integ_points, jquad, nkp
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: ikp_local

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_Wc_kp_tau_GW', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iquad, jkp

      CALL timeset(routineN, handle)

      ! Fourier transform from w to t; allocate only one (at first time point and also at
      ! first kpoint in case we do not do
      IF (jquad == 1) THEN
         NULLIFY (cfm_mat_W_kp_tau)
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

END MODULE rpa_kpoints
