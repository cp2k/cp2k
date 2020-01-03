!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate weights for correlation methods
!> \par History
!>      05.2019 Refactored from rpa_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE mp2_weights
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_comm_split_direct,&
                                              mp_sum
   USE minimax_exp,                     ONLY: get_exp_minimax_coeff
   USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff
   USE mp2_types,                       ONLY: mp2_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_weights'

   PUBLIC :: get_minimax_weights, get_clenshaw_weights, test_least_square_ft

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param unit_nr ...
!> \param homo ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param do_im_time ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param do_print ...
!> \param tau_tj ...
!> \param tau_wj ...
!> \param qs_env ...
!> \param do_gw_im_time ...
!> \param do_kpoints_cubic_RPA ...
!> \param ext_scaling ...
!> \param a_scaling ...
!> \param e_fermi ...
!> \param tj ...
!> \param wj ...
!> \param mp2_env ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_cos_tf_w_to_t ...
!> \param weights_sin_tf_t_to_w ...
!> \param homo_beta ...
!> \param dimen_ia_beta ...
!> \param Eigenval_beta ...
! **************************************************************************************************
   SUBROUTINE get_minimax_weights(para_env, unit_nr, homo, Eigenval, num_integ_points, &
                                  do_im_time, do_ri_sos_laplace_mp2, do_print, tau_tj, tau_wj, qs_env, do_gw_im_time, &
                                  do_kpoints_cubic_RPA, ext_scaling, a_scaling, e_fermi, tj, wj, mp2_env, weights_cos_tf_t_to_w, &
                                  weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, &
                                  homo_beta, dimen_ia_beta, Eigenval_beta)

      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr, homo
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points
      LOGICAL, INTENT(IN)                                :: do_im_time, do_ri_sos_laplace_mp2, &
                                                            do_print
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tau_tj, tau_wj
      TYPE(qs_environment_type), OPTIONAL, POINTER       :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_gw_im_time, do_kpoints_cubic_RPA
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: ext_scaling
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: a_scaling, e_fermi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT), OPTIONAL                           :: tj, wj
      TYPE(mp2_type), OPTIONAL, POINTER                  :: mp2_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT), OPTIONAL                           :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w
      INTEGER, INTENT(IN), OPTIONAL                      :: homo_beta, dimen_ia_beta
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: Eigenval_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_minimax_weights', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ierr, jquad, &
                                                            num_points_per_magnitude
      LOGICAL                                            :: my_do_kpoints, my_open_shell
      REAL(KIND=dp)                                      :: E_Range, Emax, Emax_beta, Emin, &
                                                            Emin_beta, max_error_min, scaling
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x_tw

      CALL timeset(routineN, handle)

      ! Test for spin unrestricted
      my_open_shell = .FALSE.
      IF (PRESENT(homo_beta) .AND. PRESENT(dimen_ia_beta) .AND. PRESENT(Eigenval_beta)) THEN
         my_open_shell = .TRUE.
      END IF

      ! Test whether all necessary variables are available
      my_do_kpoints = .FALSE.
      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         IF (.NOT. (PRESENT(do_gw_im_time) .AND. PRESENT(do_kpoints_cubic_RPA) .AND. PRESENT(ext_scaling) .AND. PRESENT(a_scaling) &
                    .AND. PRESENT(e_fermi) .AND. PRESENT(tj) .AND. PRESENT(wj) .AND. PRESENT(weights_cos_tf_t_to_w) .AND. &
                    PRESENT(weights_cos_tf_w_to_t) .AND. PRESENT(weights_sin_tf_t_to_w) .AND. PRESENT(mp2_env))) THEN
            CPABORT("Need more parameters without SOS-MP2")
         END IF
         my_do_kpoints = do_kpoints_cubic_RPA
      END IF

      IF (my_do_kpoints) THEN

         CALL gap_and_max_eig_diff_kpoints(qs_env, para_env, Emin, Emax, e_fermi)

         E_Range = Emax/Emin

      ELSE

         Emin = Eigenval(homo + 1) - Eigenval(homo)
         Emax = MAXVAL(Eigenval) - MINVAL(Eigenval)
         IF (my_open_shell) THEN
            IF (homo_beta > 0) THEN
               Emin_beta = Eigenval_beta(homo_beta + 1) - Eigenval_beta(homo_beta)
               Emax_beta = MAXVAL(Eigenval_beta) - MINVAL(Eigenval_beta)
               Emin = MIN(Emin, Emin_beta)
               Emax = MAX(Emax, Emax_beta)
            END IF
         END IF
         E_Range = Emax/Emin

      END IF

      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp
         ierr = 0
         CALL get_rpa_minimax_coeff(num_integ_points, E_Range, x_tw, ierr)

         ALLOCATE (tj(num_integ_points))
         tj = 0.0_dp

         ALLOCATE (wj(num_integ_points))
         wj = 0.0_dp

         DO jquad = 1, num_integ_points
            tj(jquad) = x_tw(jquad)
            wj(jquad) = x_tw(jquad + num_integ_points)
         END DO

         DEALLOCATE (x_tw)

         IF (unit_nr > 0 .AND. do_print) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
               "INTEG_INFO| Range for the minimax approximation:", E_Range
            WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") "INTEG_INFO| Minimax parameters:", "Weights", "Abscissas"
            DO jquad = 1, num_integ_points
               WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") wj(jquad), tj(jquad)
            END DO
            CALL m_flush(unit_nr)
         END IF

         ! scale the minimax parameters
         tj(:) = tj(:)*Emin
         wj(:) = wj(:)*Emin
      ELSE
         ! When we perform SOS-MP2, we need an additional factor of 2 for the energies (compare with mp2_laplace.F)
         ! We do not need weights etc. for the cosine transform
         ! We do not scale Emax because it is not needed for SOS-MP2
         Emin = Emin*2.0_dp
      END IF

      ! set up the minimax time grid
      IF (do_im_time .OR. do_ri_sos_laplace_mp2) THEN

         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp

         CALL get_exp_minimax_coeff(num_integ_points, E_Range, x_tw)

         ! For RPA we include already a factor of two (see later steps)
         scaling = 2.0_dp
         IF (do_ri_sos_laplace_mp2) scaling = 1.0_dp

         ALLOCATE (tau_tj(0:num_integ_points))
         tau_tj = 0.0_dp

         ALLOCATE (tau_wj(num_integ_points))
         tau_wj = 0.0_dp

         DO jquad = 1, num_integ_points
            tau_tj(jquad) = x_tw(jquad)/scaling
            tau_wj(jquad) = x_tw(jquad + num_integ_points)/scaling
         END DO

         DEALLOCATE (x_tw)

         IF (unit_nr > 0 .AND. do_print) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
               "INTEG_INFO| Range for the minimax approximation:", E_Range
            ! For testing the gap
            WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
               "INTEG_INFO| Gap:", Emin
            WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") &
               "INTEG_INFO| Minimax parameters of the time grid:", "Weights", "Abscissas"
            DO jquad = 1, num_integ_points
               WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") tau_wj(jquad), tau_tj(jquad)
            END DO
            CALL m_flush(unit_nr)
         END IF

         ! scale grid from [1,R] to [Emin,Emax]
         tau_tj(:) = tau_tj(:)/Emin
         tau_wj(:) = tau_wj(:)/Emin

         IF (.NOT. do_ri_sos_laplace_mp2) THEN
            ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
            weights_cos_tf_t_to_w = 0.0_dp

            num_points_per_magnitude = mp2_env%ri_rpa_im_time%num_points_per_magnitude
            CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, tj, &
                                              Emin, Emax, max_error_min, num_points_per_magnitude)

            IF (do_gw_im_time) THEN

               ! get the weights for the cosine transform W^c(iw) -> W^c(it)
               ALLOCATE (weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
               weights_cos_tf_w_to_t = 0.0_dp

               CALL get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, tj, &
                                                 Emin, Emax, max_error_min, num_points_per_magnitude)

               ! get the weights for the sine transform Sigma^sin(it) -> Sigma^sin(iw) (PRB 94, 165109 (2016), Eq. 71)
               ALLOCATE (weights_sin_tf_t_to_w(num_integ_points, num_integ_points))
               weights_sin_tf_t_to_w = 0.0_dp

               CALL get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, tj, &
                                                 Emin, Emax, max_error_min, num_points_per_magnitude)

               IF (unit_nr > 0) THEN
                  WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
                     "INTEG_INFO| Maximum deviation of the imag. time fit:", max_error_min
               END IF
            END IF

         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE get_minimax_weights

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param unit_nr ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param fm_mat_S ...
!> \param my_do_gw ...
!> \param ext_scaling ...
!> \param a_scaling ...
!> \param tj ...
!> \param wj ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param dimen_ia_beta ...
!> \param Eigenval_beta ...
!> \param fm_mat_S_beta ...
! **************************************************************************************************
   SUBROUTINE get_clenshaw_weights(para_env, para_env_RPA, unit_nr, homo, virtual, Eigenval, num_integ_points, &
                                   num_integ_group, color_rpa_group, fm_mat_S, my_do_gw, &
                                   ext_scaling, a_scaling, tj, wj, &
                                   homo_beta, virtual_beta, dimen_ia_beta, Eigenval_beta, fm_mat_S_beta)

      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      INTEGER, INTENT(IN)                                :: unit_nr, homo, virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      LOGICAL, INTENT(IN)                                :: my_do_gw
      REAL(KIND=dp), INTENT(IN)                          :: ext_scaling
      REAL(KIND=dp), INTENT(OUT)                         :: a_scaling
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tj, wj
      INTEGER, INTENT(IN), OPTIONAL                      :: homo_beta, virtual_beta, dimen_ia_beta
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: Eigenval_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_clenshaw_weights', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, jquad
      LOGICAL                                            :: my_open_shell

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(homo_beta) .AND. PRESENT(virtual_beta) .AND. PRESENT(dimen_ia_beta) .AND. PRESENT(Eigenval_beta) .AND. &
          PRESENT(fm_mat_S_beta)) THEN
         my_open_shell = .TRUE.
      END IF

      ! Now, start to prepare the different grid
      ALLOCATE (tj(num_integ_points))
      tj = 0.0_dp

      ALLOCATE (wj(num_integ_points))
      wj = 0.0_dp

      DO jquad = 1, num_integ_points - 1
         tj(jquad) = jquad*pi/(2.0_dp*num_integ_points)
         wj(jquad) = pi/(num_integ_points*SIN(tj(jquad))**2)
      END DO
      tj(num_integ_points) = pi/2.0_dp
      wj(num_integ_points) = pi/(2.0_dp*num_integ_points*SIN(tj(num_integ_points))**2)

      a_scaling = 1.0_dp
      IF (my_open_shell) THEN
         CALL calc_scaling_factor(a_scaling, para_env, para_env_RPA, homo, virtual, Eigenval, &
                                  num_integ_points, num_integ_group, color_rpa_group, &
                                  tj, wj, fm_mat_S, &
                                  homo_beta, virtual_beta, dimen_ia_beta, Eigenval_beta, fm_mat_S_beta)
      ELSE
         CALL calc_scaling_factor(a_scaling, para_env, para_env_RPA, homo, virtual, Eigenval, &
                                  num_integ_points, num_integ_group, color_rpa_group, &
                                  tj, wj, fm_mat_S)
      END IF

      ! for G0W0, we may set the scaling factor by hand
      IF (my_do_gw .AND. ext_scaling > 0.0_dp) THEN
         a_scaling = ext_scaling
      END IF

      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.5)') 'INTEG_INFO| Scaling parameter:', a_scaling

      wj(:) = wj(:)*a_scaling

      CALL timestop(handle)

   END SUBROUTINE get_clenshaw_weights

! **************************************************************************************************
!> \brief ...
!> \param a_scaling_ext ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param tj_ext ...
!> \param wj_ext ...
!> \param fm_mat_S ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param dimen_ia_beta ...
!> \param Eigenval_beta ...
!> \param fm_mat_S_beta ...
! **************************************************************************************************
   SUBROUTINE calc_scaling_factor(a_scaling_ext, para_env, para_env_RPA, homo, virtual, Eigenval, &
                                  num_integ_points, num_integ_group, color_rpa_group, &
                                  tj_ext, wj_ext, fm_mat_S, &
                                  homo_beta, virtual_beta, dimen_ia_beta, Eigenval_beta, fm_mat_S_beta)
      REAL(KIND=dp), INTENT(INOUT)                       :: a_scaling_ext
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      INTEGER, INTENT(IN)                                :: homo, virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tj_ext, wj_ext
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      INTEGER, INTENT(IN), OPTIONAL                      :: homo_beta, virtual_beta, dimen_ia_beta
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: Eigenval_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_scaling_factor', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icycle, jquad, nrow_local, &
                                                            nrow_local_beta
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp) :: a_high, a_low, a_scaling, conv_param, eps, first_deriv, left_term, &
         right_term, right_term_ref, right_term_ref_beta, step
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cottj, D_ia, D_ia_beta, iaia_RI, &
                                                            iaia_RI_beta, M_ia, M_ia_beta
      TYPE(cp_para_env_type), POINTER                    :: para_env_row, para_env_row_beta

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(homo_beta) .AND. &
          PRESENT(virtual_beta) .AND. &
          PRESENT(dimen_ia_beta) .AND. &
          PRESENT(Eigenval_beta) .AND. &
          PRESENT(fm_mat_S_beta)) my_open_shell = .TRUE.

      eps = 1.0E-10_dp

      ALLOCATE (cottj(num_integ_points))

      ! calculate the cotangent of the abscissa tj
      DO jquad = 1, num_integ_points
         cottj(jquad) = 1.0_dp/TAN(tj_ext(jquad))
      END DO

      CALL calc_ia_ia_integrals(para_env_RPA, homo, virtual, nrow_local, right_term_ref, Eigenval, D_ia, iaia_RI, M_ia, fm_mat_S, &
                                para_env_row)

      ! In the open shell case do point 1-2-3 for the beta spin
      IF (my_open_shell) THEN
         CALL calc_ia_ia_integrals(para_env_RPA, homo_beta, virtual_beta, nrow_local_beta, right_term_ref_beta, Eigenval_beta, &
                                   D_ia_beta, iaia_RI_beta, M_ia_beta, fm_mat_S_beta, para_env_row_beta)

         right_term_ref = right_term_ref + right_term_ref_beta
      END IF

      ! bcast the result
      IF (para_env%mepos == 0) THEN
         CALL mp_bcast(right_term_ref, 0, para_env%group)
      ELSE
         right_term_ref = 0.0_dp
         CALL mp_bcast(right_term_ref, 0, para_env%group)
      END IF

      ! 5) start iteration for solving the non-linear equation by bisection
      ! find limit, here step=0.5 seems a good compromise
      conv_param = 100.0_dp*EPSILON(right_term_ref)
      step = 0.5_dp
      a_low = 0.0_dp
      a_high = step
      right_term = -right_term_ref
      DO icycle = 1, num_integ_points*2
         a_scaling = a_high

         CALL calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                M_ia, cottj, wj_ext, D_ia, D_ia_beta, M_ia_beta, &
                                nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                para_env, para_env_row, para_env_row_beta)
         left_term = left_term/4.0_dp/pi*a_scaling

         IF (ABS(left_term) > ABS(right_term) .OR. ABS(left_term + right_term) <= conv_param) EXIT
         a_low = a_high
         a_high = a_high + step

      END DO

      IF (ABS(left_term + right_term) >= conv_param) THEN
         IF (a_scaling >= 2*num_integ_points*step) THEN
            a_scaling = 1.0_dp
         ELSE

            DO icycle = 1, num_integ_points*2
               a_scaling = (a_low + a_high)/2.0_dp

               CALL calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                      M_ia, cottj, wj_ext, D_ia, D_ia_beta, M_ia_beta, &
                                      nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                      para_env, para_env_row, para_env_row_beta)
               left_term = left_term/4.0_dp/pi*a_scaling

               IF (ABS(left_term) > ABS(right_term)) THEN
                  a_high = a_scaling
               ELSE
                  a_low = a_scaling
               END IF

               IF (ABS(a_high - a_low) < 1.0e-5_dp) EXIT

            END DO

         END IF
      END IF

      a_scaling_ext = a_scaling
      CALL mp_bcast(a_scaling_ext, 0, para_env%group)

      DEALLOCATE (cottj)
      DEALLOCATE (iaia_RI)
      DEALLOCATE (D_ia)
      DEALLOCATE (M_ia)
      CALL cp_para_env_release(para_env_row)

      IF (my_open_shell) THEN
         DEALLOCATE (iaia_RI_beta)
         DEALLOCATE (D_ia_beta)
         DEALLOCATE (M_ia_beta)
         CALL cp_para_env_release(para_env_row_beta)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_scaling_factor

! **************************************************************************************************
!> \brief ...
!> \param para_env_RPA ...
!> \param homo ...
!> \param virtual ...
!> \param nrow_local ...
!> \param right_term_ref ...
!> \param Eigenval ...
!> \param D_ia ...
!> \param iaia_RI ...
!> \param M_ia ...
!> \param fm_mat_S ...
!> \param para_env_row ...
! **************************************************************************************************
   SUBROUTINE calc_ia_ia_integrals(para_env_RPA, homo, virtual, nrow_local, right_term_ref, Eigenval, &
                                   D_ia, iaia_RI, M_ia, fm_mat_S, para_env_row)

      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      INTEGER, INTENT(IN)                                :: homo, virtual
      INTEGER, INTENT(OUT)                               :: nrow_local
      REAL(KIND=dp), INTENT(OUT)                         :: right_term_ref
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: D_ia, iaia_RI, M_ia
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      TYPE(cp_para_env_type), POINTER                    :: para_env_row

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ia_ia_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: avirt, color_col, color_row, comm_col, &
                                                            comm_row, handle, i_global, iiB, iocc, &
                                                            jjB, ncol_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: eigen_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: iaia_RI_dp
      TYPE(cp_para_env_type), POINTER                    :: para_env_col

      CALL timeset(routineN, handle)

      ! calculate the (ia|ia) RI integrals
      ! ----------------------------------
      ! 1) get info fm_mat_S
      !XXX CALL cp_fm_to_fm(source=fm_mat_S,destination=fm_mat_G)
      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! allocate the local buffer of iaia_RI integrals (dp kind)
      ALLOCATE (iaia_RI_dp(nrow_local))
      iaia_RI_dp = 0.0_dp

      ! 2) perform the local multiplication SUM_K (ia|K)*(ia|K)
      DO jjB = 1, ncol_local
         DO iiB = 1, nrow_local
            iaia_RI_dp(iiB) = iaia_RI_dp(iiB) + fm_mat_S%local_data(iiB, jjB)*fm_mat_S%local_data(iiB, jjB)
         END DO
      END DO

      ! 3) sum the result with the processes of the RPA_group having the same rows
      !          _______K_______               _
      !         |   |   |   |   |             | |
      !     --> | 1 | 5 | 9 | 13|   SUM -->   | |
      !         |___|__ |___|___|             |_|
      !         |   |   |   |   |             | |
      !     --> | 2 | 6 | 10| 14|   SUM -->   | |
      !      ia |___|___|___|___|             |_|   (ia|ia)_RI
      !         |   |   |   |   |             | |
      !     --> | 3 | 7 | 11| 15|   SUM -->   | |
      !         |___|___|___|___|             |_|
      !         |   |   |   |   |             | |
      !     --> | 4 | 8 | 12| 16|   SUM -->   | |
      !         |___|___|___|___|             |_|
      !

      color_row = fm_mat_S%matrix_struct%context%mepos(1)
      CALL mp_comm_split_direct(para_env_RPA%group, comm_row, color_row)
      NULLIFY (para_env_row)
      CALL cp_para_env_create(para_env_row, comm_row)

      CALL mp_sum(iaia_RI_dp, para_env_row%group)

      ! convert the iaia_RI_dp into double-double precision
      ALLOCATE (iaia_RI(nrow_local))
      DO iiB = 1, nrow_local
         iaia_RI(iiB) = iaia_RI_dp(iiB)
      END DO
      DEALLOCATE (iaia_RI_dp)

      ! 4) calculate the right hand term, D_ia is the matrix containing the
      ! orbital energy differences, M_ia is the diagonal of the full RPA 'excitation'
      ! matrix
      ALLOCATE (D_ia(nrow_local))

      ALLOCATE (M_ia(nrow_local))

      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)

         iocc = MAX(1, i_global - 1)/virtual + 1
         avirt = i_global - (iocc - 1)*virtual
         eigen_diff = Eigenval(avirt + homo) - Eigenval(iocc)

         D_ia(iiB) = eigen_diff
      END DO

      DO iiB = 1, nrow_local
         M_ia(iiB) = D_ia(iiB)*D_ia(iiB) + 2.0_dp*D_ia(iiB)*iaia_RI(iiB)
      END DO

      right_term_ref = 0.0_dp
      DO iiB = 1, nrow_local
         right_term_ref = right_term_ref + (SQRT(M_ia(iiB)) - D_ia(iiB) - iaia_RI(iiB))
      END DO
      right_term_ref = right_term_ref/2.0_dp

      ! sum the result with the processes of the RPA_group having the same col
      color_col = fm_mat_S%matrix_struct%context%mepos(2)
      CALL mp_comm_split_direct(para_env_RPA%group, comm_col, color_col)
      NULLIFY (para_env_col)
      CALL cp_para_env_create(para_env_col, comm_col)

      ! allocate communication array for columns
      CALL mp_sum(right_term_ref, para_env_col%group)

      CALL cp_para_env_release(para_env_col)

      CALL timestop(handle)

   END SUBROUTINE calc_ia_ia_integrals

! **************************************************************************************************
!> \brief ...
!> \param a_scaling ...
!> \param left_term ...
!> \param first_deriv ...
!> \param num_integ_points ...
!> \param my_open_shell ...
!> \param M_ia ...
!> \param cottj ...
!> \param wj ...
!> \param D_ia ...
!> \param D_ia_beta ...
!> \param M_ia_beta ...
!> \param nrow_local ...
!> \param nrow_local_beta ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param para_env ...
!> \param para_env_row ...
!> \param para_env_row_beta ...
! **************************************************************************************************
   SUBROUTINE calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                M_ia, cottj, wj, D_ia, D_ia_beta, M_ia_beta, &
                                nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                para_env, para_env_row, para_env_row_beta)
      REAL(KIND=dp), INTENT(IN)                          :: a_scaling
      REAL(KIND=dp), INTENT(INOUT)                       :: left_term, first_deriv
      INTEGER, INTENT(IN)                                :: num_integ_points
      LOGICAL, INTENT(IN)                                :: my_open_shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: M_ia, cottj, wj, D_ia, D_ia_beta, &
                                                            M_ia_beta
      INTEGER, INTENT(IN)                                :: nrow_local, nrow_local_beta, &
                                                            num_integ_group, color_rpa_group
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_row, para_env_row_beta

      INTEGER                                            :: iiB, jquad
      REAL(KIND=dp)                                      :: first_deriv_beta, left_term_beta, omega

      left_term = 0.0_dp
      first_deriv = 0.0_dp
      left_term_beta = 0.0_dp
      first_deriv_beta = 0.0_dp
      DO jquad = 1, num_integ_points
         ! parallelize over integration points
         IF (MODULO(jquad, num_integ_group) /= color_rpa_group) CYCLE
         omega = a_scaling*cottj(jquad)

         DO iiB = 1, nrow_local
            ! parallelize over ia elements in the para_env_row group
            IF (MODULO(iiB, para_env_row%num_pe) /= para_env_row%mepos) CYCLE
            ! calculate left_term
            left_term = left_term + wj(jquad)* &
                        (LOG(1.0_dp + (M_ia(iiB) - D_ia(iiB)**2)/(omega**2 + D_ia(iiB)**2)) - &
                         (M_ia(iiB) - D_ia(iiB)**2)/(omega**2 + D_ia(iiB)**2))
            first_deriv = first_deriv + wj(jquad)*cottj(jquad)**2* &
                          ((-M_ia(iiB) + D_ia(iiB)**2)**2/((omega**2 + D_ia(iiB)**2)**2*(omega**2 + M_ia(iiB))))
         END DO

         IF (my_open_shell) THEN
            DO iiB = 1, nrow_local_beta
               ! parallelize over ia elements in the para_env_row group
               IF (MODULO(iiB, para_env_row_beta%num_pe) /= para_env_row_beta%mepos) CYCLE
               ! calculate left_term
               left_term_beta = left_term_beta + wj(jquad)* &
                                (LOG(1.0_dp + (M_ia_beta(iiB) - D_ia_beta(iiB)**2)/(omega**2 + D_ia_beta(iiB)**2)) - &
                                 (M_ia_beta(iiB) - D_ia_beta(iiB)**2)/(omega**2 + D_ia_beta(iiB)**2))
               first_deriv_beta = &
                  first_deriv_beta + wj(jquad)*cottj(jquad)**2* &
                  ((-M_ia_beta(iiB) + D_ia_beta(iiB)**2)**2/((omega**2 + D_ia_beta(iiB)**2)**2*(omega**2 + M_ia_beta(iiB))))
            END DO
         END IF

      END DO

      ! sum the contribution from all proc, starting form the row group
      CALL mp_sum(left_term, para_env%group)
      CALL mp_sum(first_deriv, para_env%group)

      IF (my_open_shell) THEN
         CALL mp_sum(left_term_beta, para_env%group)
         CALL mp_sum(first_deriv_beta, para_env%group)

         left_term = left_term + left_term_beta
         first_deriv = first_deriv + first_deriv_beta
      END IF

   END SUBROUTINE calculate_objfunc

! **************************************************************************************************
!> \brief Calculate integration weights for the tau grid (in dependency of the omega node)
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param weights_cos_tf_t_to_w ...
!> \param omega_tj ...
!> \param E_min ...
!> \param E_max ...
!> \param max_error ...
!> \param num_points_per_magnitude ...
! **************************************************************************************************
   SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude)

      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj
      REAL(KIND=dp), INTENT(IN)                          :: E_min, E_max
      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_t_to_w', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      CALL timeset(routineN, handle)

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (INT(LOG10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = MAX(num_x_nodes, num_integ_points)

      ALLOCATE (x_values(num_x_nodes))
      x_values = 0.0_dp
      ALLOCATE (y_values(num_x_nodes))
      y_values = 0.0_dp
      ALLOCATE (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0_dp
      ALLOCATE (tau_wj_work(num_integ_points))
      tau_wj_work = 0.0_dp
      ALLOCATE (work_array(2*num_integ_points))
      work_array = 0.0_dp
      ALLOCATE (sing_values(num_integ_points))
      sing_values = 0.0_dp
      ALLOCATE (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0_dp
      ALLOCATE (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0_dp
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      ALLOCATE (work(lwork))
      work = 0.0_dp
      ALLOCATE (iwork(8*num_integ_points))
      iwork = 0
      ALLOCATE (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0_dp
      ALLOCATE (vec_UTy(num_x_nodes))
      vec_UTy = 0.0_dp

      max_error = 0.0_dp

      ! loop over all omega frequency points
      DO jquad = 1, num_integ_points

         chi2_min_jquad = 100.0_dp

         ! set the x-values logarithmically in the interval [Emin,Emax]
         multiplicator = (E_max/E_min)**(1.0_dp/(REAL(num_x_nodes, KIND=dp) - 1.0_dp))
         DO iii = 1, num_x_nodes
            x_values(iii) = E_min*multiplicator**(iii - 1)
         END DO

         omega = omega_tj(jquad)

         ! y=2x/(x^2+omega_k^2)
         DO iii = 1, num_x_nodes
            y_values(iii) = 2.0_dp*x_values(iii)/((x_values(iii))**2 + omega**2)
         END DO

         ! calculate mat_A
         DO jjj = 1, num_integ_points
            DO iii = 1, num_x_nodes
               mat_A(iii, jjj) = COS(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
            END DO
         END DO

         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, info)

         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            END DO
         END DO

         ! 2) U^T y
         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0_dp, tau_wj_work, num_integ_points)

         weights_cos_tf_t_to_w(jquad, :) = tau_wj_work(:)

         CALL calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                      y_values, num_integ_points, num_x_nodes)

      END DO ! jquad

      DEALLOCATE (x_values, y_values, mat_A, tau_wj_work, work_array, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

      CALL timestop(handle)

   END SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w

! **************************************************************************************************
!> \brief Calculate integration weights for the tau grid (in dependency of the omega node)
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param weights_sin_tf_t_to_w ...
!> \param omega_tj ...
!> \param E_min ...
!> \param E_max ...
!> \param max_error ...
!> \param num_points_per_magnitude ...
! **************************************************************************************************
   SUBROUTINE get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude)

      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_sin_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj
      REAL(KIND=dp), INTENT(IN)                          :: E_min, E_max
      REAL(KIND=dp), INTENT(OUT)                         :: max_error
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_sin_tf_t_to_w', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      CALL timeset(routineN, handle)

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (INT(LOG10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = MAX(num_x_nodes, num_integ_points)

      ALLOCATE (x_values(num_x_nodes))
      x_values = 0.0_dp
      ALLOCATE (y_values(num_x_nodes))
      y_values = 0.0_dp
      ALLOCATE (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0_dp
      ALLOCATE (tau_wj_work(num_integ_points))
      tau_wj_work = 0.0_dp
      ALLOCATE (work_array(2*num_integ_points))
      work_array = 0.0_dp
      ALLOCATE (sing_values(num_integ_points))
      sing_values = 0.0_dp
      ALLOCATE (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0_dp
      ALLOCATE (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0_dp
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      ALLOCATE (work(lwork))
      work = 0.0_dp
      ALLOCATE (iwork(8*num_integ_points))
      iwork = 0
      ALLOCATE (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0_dp
      ALLOCATE (vec_UTy(num_x_nodes))
      vec_UTy = 0.0_dp

      max_error = 0.0_dp

      ! loop over all omega frequency points
      DO jquad = 1, num_integ_points

         chi2_min_jquad = 100.0_dp

         ! set the x-values logarithmically in the interval [Emin,Emax]
         multiplicator = (E_max/E_min)**(1.0_dp/(REAL(num_x_nodes, KIND=dp) - 1.0_dp))
         DO iii = 1, num_x_nodes
            x_values(iii) = E_min*multiplicator**(iii - 1)
         END DO

         omega = omega_tj(jquad)

         ! y=2x/(x^2+omega_k^2)
         DO iii = 1, num_x_nodes
!            y_values(iii) = 2.0_dp*x_values(iii)/((x_values(iii))**2+omega**2)
            y_values(iii) = 2.0_dp*omega/((x_values(iii))**2 + omega**2)
         END DO

         ! calculate mat_A
         DO jjj = 1, num_integ_points
            DO iii = 1, num_x_nodes
               mat_A(iii, jjj) = SIN(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
            END DO
         END DO

         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, info)

         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            END DO
         END DO

         ! 2) U^T y
         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0_dp, tau_wj_work, num_integ_points)

         weights_sin_tf_t_to_w(jquad, :) = tau_wj_work(:)

         CALL calc_max_error_fit_tau_grid_with_sine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                    y_values, num_integ_points, num_x_nodes)

      END DO ! jquad

      DEALLOCATE (x_values, y_values, mat_A, tau_wj_work, work_array, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

      CALL timestop(handle)

   END SUBROUTINE get_l_sq_wghts_sin_tf_t_to_w

! **************************************************************************************************
!> \brief ...
!> \param max_error ...
!> \param omega ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param x_values ...
!> \param y_values ...
!> \param num_integ_points ...
!> \param num_x_nodes ...
! **************************************************************************************************
   PURE SUBROUTINE calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                           y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp), INTENT(INOUT)                       :: max_error, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work, x_values, y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_tau_grid_with_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: kkk
      REAL(KIND=dp)                                      :: func_val, func_val_temp, max_error_tmp

      max_error_tmp = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_tau_grid_cosine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         IF (ABS(y_values(kkk) - func_val) > max_error_tmp) THEN
            max_error_tmp = ABS(y_values(kkk) - func_val)
            func_val_temp = func_val
         END IF

      END DO

      IF (max_error_tmp > max_error) THEN

         max_error = max_error_tmp

      END IF

   END SUBROUTINE calc_max_error_fit_tau_grid_with_cosine

! **************************************************************************************************
!> \brief Evaluate fit function when calculating tau grid for cosine transform
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param omega ...
! **************************************************************************************************
   PURE SUBROUTINE eval_fit_func_tau_grid_cosine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

      REAL(KIND=dp), INTENT(OUT)                         :: func_val
      REAL(KIND=dp), INTENT(IN)                          :: x_value
      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_fit_func_tau_grid_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val + tau_wj_work(iii)*COS(omega*tau_tj(iii))*EXP(-x_value*tau_tj(iii))

      END DO

   END SUBROUTINE eval_fit_func_tau_grid_cosine

! **************************************************************************************************
!> \brief Evaluate fit function when calculating tau grid for sine transform
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param omega ...
! **************************************************************************************************
   PURE SUBROUTINE eval_fit_func_tau_grid_sine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

      REAL(KIND=dp), INTENT(INOUT)                       :: func_val
      REAL(KIND=dp), INTENT(IN)                          :: x_value
      INTEGER, INTENT(in)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_fit_func_tau_grid_sine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val + tau_wj_work(iii)*SIN(omega*tau_tj(iii))*EXP(-x_value*tau_tj(iii))

      END DO

   END SUBROUTINE eval_fit_func_tau_grid_sine

! **************************************************************************************************
!> \brief ...
!> \param max_error ...
!> \param omega ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param x_values ...
!> \param y_values ...
!> \param num_integ_points ...
!> \param num_x_nodes ...
! **************************************************************************************************
   PURE SUBROUTINE calc_max_error_fit_tau_grid_with_sine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                         y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp), INTENT(INOUT)                       :: max_error, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work, x_values, y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_tau_grid_with_sine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: kkk
      REAL(KIND=dp)                                      :: func_val, func_val_temp, max_error_tmp

      max_error_tmp = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_tau_grid_sine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         IF (ABS(y_values(kkk) - func_val) > max_error_tmp) THEN
            max_error_tmp = ABS(y_values(kkk) - func_val)
            func_val_temp = func_val
         END IF

      END DO

      IF (max_error_tmp > max_error) THEN

         max_error = max_error_tmp

      END IF

   END SUBROUTINE calc_max_error_fit_tau_grid_with_sine

! **************************************************************************************************
!> \brief test the singular value decomposition for the computation of integration weights for the
!>         Fourier transform between time and frequency grid in cubic-scaling RPA
!> \param nR ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE test_least_square_ft(nR, iw)
      INTEGER, INTENT(IN)                                :: nR, iw

      INTEGER                                            :: ierr, iR, jquad, num_integ_points
      REAL(KIND=dp)                                      :: max_error, multiplicator, Rc, Rc_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tau_wj, tj, wj, x_tw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_cos_tf_t_to_w

      Rc_max = 1.0E+7

      multiplicator = Rc_max**(1.0_dp/(REAL(nR, KIND=dp) - 1.0_dp))

      DO num_integ_points = 1, 20

         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp
         ALLOCATE (tau_tj(num_integ_points))
         tau_tj = 0.0_dp
         ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
         weights_cos_tf_t_to_w = 0.0_dp
         ALLOCATE (tau_wj(num_integ_points))
         tau_wj = 0.0_dp
         ALLOCATE (tj(num_integ_points))
         tj = 0.0_dp
         ALLOCATE (wj(num_integ_points))
         wj = 0.0_dp

         DO iR = 0, nR - 1

            Rc = 2.0_dp*multiplicator**iR

            ierr = 0
            CALL get_rpa_minimax_coeff(num_integ_points, Rc, x_tw, ierr, print_warning=.FALSE.)

            DO jquad = 1, num_integ_points
               tj(jquad) = x_tw(jquad)
               wj(jquad) = x_tw(jquad + num_integ_points)
            END DO

            x_tw = 0.0_dp

            CALL get_exp_minimax_coeff(num_integ_points, Rc, x_tw)

            DO jquad = 1, num_integ_points
               tau_tj(jquad) = x_tw(jquad)/2.0_dp
               tau_wj(jquad) = x_tw(jquad + num_integ_points)/2.0_dp
            END DO

            CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, tj, &
                                              1.0_dp, Rc, max_error, 200)

            IF (iw > 0) THEN
               WRITE (iw, '(T2, I3, F12.1, ES12.3)') num_integ_points, Rc, max_error
            END IF

         END DO

         DEALLOCATE (x_tw, tau_tj, weights_cos_tf_t_to_w, tau_wj, wj, tj)

      END DO

   END SUBROUTINE test_least_square_ft

! **************************************************************************************************
!> \brief ...
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param weights_cos_tf_w_to_t ...
!> \param omega_tj ...
!> \param E_min ...
!> \param E_max ...
!> \param max_error ...
!> \param num_points_per_magnitude ...
! **************************************************************************************************
   SUBROUTINE get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude)

      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_w_to_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj
      REAL(KIND=dp), INTENT(IN)                          :: E_min, E_max
      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_w_to_t', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega, &
                                                            tau, x_value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: omega_wj_work, sing_values, vec_UTy, &
                                                            work, work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      CALL timeset(routineN, handle)

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (INT(LOG10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = MAX(num_x_nodes, num_integ_points)

      ALLOCATE (x_values(num_x_nodes))
      x_values = 0.0_dp
      ALLOCATE (y_values(num_x_nodes))
      y_values = 0.0_dp
      ALLOCATE (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0_dp
      ALLOCATE (omega_wj_work(num_integ_points))
      omega_wj_work = 0.0_dp
      ALLOCATE (work_array(2*num_integ_points))
      work_array = 0.0_dp
      ALLOCATE (sing_values(num_integ_points))
      sing_values = 0.0_dp
      ALLOCATE (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0_dp
      ALLOCATE (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0_dp
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points + 12*num_integ_points + 2*num_x_nodes
      ALLOCATE (work(lwork))
      work = 0.0_dp
      ALLOCATE (iwork(8*num_integ_points))
      iwork = 0
      ALLOCATE (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0_dp
      ALLOCATE (vec_UTy(num_x_nodes))
      vec_UTy = 0.0_dp

      ! set the x-values logarithmically in the interval [Emin,Emax]
      multiplicator = (E_max/E_min)**(1.0_dp/(REAL(num_x_nodes, KIND=dp) - 1.0_dp))
      DO iii = 1, num_x_nodes
         x_values(iii) = E_min*multiplicator**(iii - 1)
      END DO

      max_error = 0.0_dp

      ! loop over all tau time points
      DO jquad = 1, num_integ_points

         chi2_min_jquad = 100.0_dp

         tau = tau_tj(jquad)

         ! y=exp(-x*|tau_k|)
         DO iii = 1, num_x_nodes
            y_values(iii) = EXP(-x_values(iii)*tau)
         END DO

         ! calculate mat_A
         DO jjj = 1, num_integ_points
            DO iii = 1, num_x_nodes
               omega = omega_tj(jjj)
               x_value = x_values(iii)
               mat_A(iii, jjj) = COS(tau*omega)*2.0_dp*x_value/(x_value**2 + omega**2)
            END DO
         END DO

         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, info)

         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            END DO
         END DO

         ! 2) U^T y
         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0_dp, omega_wj_work, num_integ_points)

         weights_cos_tf_w_to_t(jquad, :) = omega_wj_work(:)

         CALL calc_max_error_fit_omega_grid_with_cosine(max_error, tau, omega_tj, omega_wj_work, x_values, &
                                                        y_values, num_integ_points, num_x_nodes)

      END DO ! jquad

      DEALLOCATE (x_values, y_values, mat_A, omega_wj_work, work_array, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

      CALL timestop(handle)

   END SUBROUTINE get_l_sq_wghts_cos_tf_w_to_t

! **************************************************************************************************
!> \brief ...
!> \param max_error ...
!> \param tau ...
!> \param omega_tj ...
!> \param omega_wj_work ...
!> \param x_values ...
!> \param y_values ...
!> \param num_integ_points ...
!> \param num_x_nodes ...
! **************************************************************************************************
   SUBROUTINE calc_max_error_fit_omega_grid_with_cosine(max_error, tau, omega_tj, omega_wj_work, x_values, &
                                                        y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      REAL(KIND=dp), INTENT(IN)                          :: tau
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj, omega_wj_work, x_values, &
                                                            y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_omega_grid_with_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, kkk
      REAL(KIND=dp)                                      :: func_val, func_val_temp, max_error_tmp

      CALL timeset(routineN, handle)

      max_error_tmp = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_omega_grid_cosine(func_val, x_values(kkk), num_integ_points, omega_tj, omega_wj_work, tau)

         IF (ABS(y_values(kkk) - func_val) > max_error_tmp) THEN
            max_error_tmp = ABS(y_values(kkk) - func_val)
            func_val_temp = func_val
         END IF

      END DO

      IF (max_error_tmp > max_error) THEN

         max_error = max_error_tmp

      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_max_error_fit_omega_grid_with_cosine

! **************************************************************************************************
!> \brief ...
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param omega_tj ...
!> \param omega_wj_work ...
!> \param tau ...
! **************************************************************************************************
   PURE SUBROUTINE eval_fit_func_omega_grid_cosine(func_val, x_value, num_integ_points, omega_tj, omega_wj_work, tau)
      REAL(KIND=dp), INTENT(OUT)                         :: func_val
      REAL(KIND=dp), INTENT(IN)                          :: x_value
      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj, omega_wj_work
      REAL(KIND=dp), INTENT(IN)                          :: tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_fit_func_omega_grid_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii
      REAL(KIND=dp)                                      :: omega

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         omega = omega_tj(iii)
         func_val = func_val + omega_wj_work(iii)*COS(tau*omega)*2.0_dp*x_value/(x_value**2 + omega**2)

      END DO

   END SUBROUTINE eval_fit_func_omega_grid_cosine

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env ...
!> \param gap ...
!> \param max_eig_diff ...
!> \param e_fermi ...
! **************************************************************************************************
   SUBROUTINE gap_and_max_eig_diff_kpoints(qs_env, para_env, gap, max_eig_diff, e_fermi)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(OUT)                         :: gap, max_eig_diff, e_fermi

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gap_and_max_eig_diff_kpoints', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, homo, ikpgr, ispin, kplocal, &
                                                            nmo, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      REAL(KIND=dp)                                      :: e_homo, e_homo_temp, e_lumo, e_lumo_temp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_homo_array, e_lumo_array, &
                                                            max_eig_diff_array
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      kpoints=kpoint)

      mo_set => kpoint%kp_env(1)%kpoint_env%mos(1, 1)%mo_set
      CALL get_mo_set(mo_set, nmo=nmo)

      CALL get_kpoint_info(kpoint, kp_range=kp_range)
      kplocal = kp_range(2) - kp_range(1) + 1

      gap = 1000.0_dp
      max_eig_diff = 0.0_dp
      e_homo = -1000.0_dp
      e_lumo = 1000.0_dp

      DO ikpgr = 1, kplocal
         kp => kpoint%kp_env(ikpgr)%kpoint_env
         nspin = SIZE(kp%mos, 2)
         DO ispin = 1, nspin
            mo_set => kp%mos(1, ispin)%mo_set
            CALL get_mo_set(mo_set, eigenvalues=eigenvalues, homo=homo)
            e_homo_temp = eigenvalues(homo)
            e_lumo_temp = eigenvalues(homo + 1)

            IF (e_homo_temp > e_homo) e_homo = e_homo_temp
            IF (e_lumo_temp < e_lumo) e_lumo = e_lumo_temp
            IF (eigenvalues(nmo) - eigenvalues(1) > max_eig_diff) max_eig_diff = eigenvalues(nmo) - eigenvalues(1)

         END DO
      END DO

      ALLOCATE (e_homo_array(0:para_env%num_pe - 1))
      e_homo_array = 0.0_dp
      e_homo_array(para_env%mepos) = e_homo

      ALLOCATE (e_lumo_array(0:para_env%num_pe - 1))
      e_lumo_array = 0.0_dp
      e_lumo_array(para_env%mepos) = e_lumo

      ALLOCATE (max_eig_diff_array(0:para_env%num_pe - 1))
      max_eig_diff_array = 0.0_dp
      max_eig_diff_array(para_env%mepos) = max_eig_diff

      CALL mp_sum(e_homo_array, para_env%group)

      CALL mp_sum(e_lumo_array, para_env%group)

      CALL mp_sum(max_eig_diff_array, para_env%group)

      gap = MINVAL(e_lumo_array) - MAXVAL(e_homo_array)

      e_fermi = (MAXVAL(e_homo_array) + MINVAL(e_lumo_array))/2.0_dp

      max_eig_diff = MAXVAL(max_eig_diff_array)

      DEALLOCATE (max_eig_diff_array, e_homo_array, e_lumo_array)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE mp2_weights
