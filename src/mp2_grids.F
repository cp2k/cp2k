!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate frequency and time grids (integration points and weights)
!>        for correlation methods
!> \par History
!>      05.2019 Refactored from rpa_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE mp2_grids
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE minimax_exp,                     ONLY: get_exp_minimax_coeff
   USE minimax_exp_gw,                  ONLY: get_exp_minimax_coeff_gw
   USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff,&
                                              get_rpa_minimax_coeff_larger_grid
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_grids'

   PUBLIC :: get_minimax_grid, get_clenshaw_grid, test_least_square_ft, get_l_sq_wghts_cos_tf_t_to_w, &
             get_l_sq_wghts_cos_tf_w_to_t, get_l_sq_wghts_sin_tf_t_to_w, get_matsubara_grid

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param unit_nr ...
!> \param homo ...
!> \param virtual ...
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
!> \param e_fermi ...
!> \param tj ...
!> \param wj ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_cos_tf_w_to_t ...
!> \param weights_sin_tf_t_to_w ...
!> \param regularization ...
! **************************************************************************************************
   SUBROUTINE get_minimax_grid(para_env, unit_nr, homo, virtual, Eigenval, num_integ_points, &
                               do_im_time, do_ri_sos_laplace_mp2, do_print, tau_tj, tau_wj, qs_env, do_gw_im_time, &
                               do_kpoints_cubic_RPA, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, regularization)

      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points
      LOGICAL, INTENT(IN)                                :: do_im_time, do_ri_sos_laplace_mp2, &
                                                            do_print
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tau_tj, tau_wj
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: do_gw_im_time, do_kpoints_cubic_RPA
      REAL(KIND=dp), INTENT(OUT)                         :: e_fermi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: tj, wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: regularization

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_minimax_grid'
      INTEGER, PARAMETER                                 :: num_points_per_magnitude = 200

      INTEGER                                            :: dimen, handle, ierr, ispin, jquad, &
                                                            nspins, orb1, orb2
      LOGICAL                                            :: my_do_kpoints, my_open_shell
      REAL(KIND=dp)                                      :: Emax, Emin, max_error_min, my_E_Range, &
                                                            my_regularization, scaling
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: x_tw
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ! Test for spin unrestricted
      nspins = SIZE(homo)
      my_open_shell = (nspins == 2)

      ! Test whether all necessary variables are available
      my_do_kpoints = .FALSE.
      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         my_do_kpoints = do_kpoints_cubic_RPA
      END IF

      my_regularization = 0.0_dp
      IF (PRESENT(regularization)) THEN
         my_regularization = regularization
      END IF

      IF (my_do_kpoints) THEN
         CALL gap_and_max_eig_diff_kpoints(qs_env, para_env, Emin, Emax, e_fermi)
         my_E_Range = Emax/Emin
      ELSE
         ! Consider partially occupied orbitals
         IF (qs_env%mp2_env%E_range <= 1.0_dp .OR. qs_env%mp2_env%E_gap <= 0.0_dp) THEN
            dimen = SIZE(Eigenval, 1)
            Emin = HUGE(dp)
            Emax = 0.0_dp
            DO ispin = 1, nspins
               IF (homo(ispin) > 0 .AND. virtual(ispin) > 0) THEN
                  ! First guess for minimal energy difference is highest fully occupied and the lowest fully unoccpied orbital
                  Emin = MIN(Emin, Eigenval(homo(ispin) + 1, 1, ispin) - Eigenval(dimen - virtual(ispin), 1, ispin))
                  ! Then, consider all pairs of partially occupied orbitals
                  ! Ignore excitations to the same orbital because of weight factors f_i-f_j
                  ! Keep only unique pairs because of factors (f_i-f_j)*(e_i-e_j)
                  DO orb1 = dimen - virtual(ispin) + 1, homo(ispin)
                     DO orb2 = orb1 + 1, homo(ispin)
                        Emin = MIN(Emin, Eigenval(orb2, 1, ispin) - Eigenval(orb1, 1, ispin))
                     END DO
                  END DO
                  Emax = MAX(Emax, MAXVAL(Eigenval(:, :, ispin)) - MINVAL(Eigenval(:, :, ispin)))
               END IF
            END DO
            my_E_Range = Emax/Emin
            qs_env%mp2_env%e_range = my_e_range
            qs_env%mp2_env%e_gap = Emin

            ! Store range and gap for successive RPA calculations for consistency (we ignore dependence on these parameters)
            CALL get_qs_env(qs_env, input=input)
            CALL section_vals_val_set(input, "DFT%XC%WF_CORRELATION%E_RANGE", r_val=my_e_range)
            CALL section_vals_val_set(input, "DFT%XC%WF_CORRELATION%E_GAP", r_val=emin)
         ELSE
            my_E_range = qs_env%mp2_env%E_range
            Emin = qs_env%mp2_env%E_gap
            Emax = Emin*my_E_range
         END IF
      END IF

      IF (num_integ_points > 20 .AND. my_E_Range < 100.0_dp) THEN
         IF (unit_nr > 0) &
            CALL cp_warn(__LOCATION__, &
                         "You requested a large minimax grid (> 20 points) for a small minimax range R (R < 100). "// &
                         "That may lead to numerical "// &
                         "instabilities when computing minimax grid weights. You can prevent small ranges by choosing "// &
                         "a larger basis set with higher angular momenta or alternatively using all-electron calculations.")
      END IF

      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp
         ierr = 0
         IF (num_integ_points .LE. 20) THEN
            CALL get_rpa_minimax_coeff(num_integ_points, my_E_Range, x_tw, ierr)
         ELSE
            CALL get_rpa_minimax_coeff_larger_grid(num_integ_points, my_E_Range, x_tw)
         END IF

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
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MINIMAX_INFO| Number of integration points:", num_integ_points
            WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
               "MINIMAX_INFO| Gap for the minimax approximation:", Emin
            WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
               "MINIMAX_INFO| Range for the minimax approximation:", my_E_Range
            WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") "MINIMAX_INFO| Minimax parameters:", "Weights", "Abscissas"
            DO jquad = 1, num_integ_points
               WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") wj(jquad), tj(jquad)
            END DO
            CALL m_flush(unit_nr)
         END IF

         ! scale the minimax parameters
         tj(:) = tj(:)*Emin
         wj(:) = wj(:)*Emin*0.5_dp
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

         IF (num_integ_points .LE. 20) THEN
            CALL get_exp_minimax_coeff(num_integ_points, my_E_Range, x_tw)
         ELSE
            CALL get_exp_minimax_coeff_gw(num_integ_points, my_E_Range, x_tw)
         END IF

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
               "MINIMAX_INFO| Range for the minimax approximation:", my_E_Range
            ! For testing the gap
            WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
               "MINIMAX_INFO| Gap:", Emin
            WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") &
               "MINIMAX_INFO| Minimax parameters of the time grid:", "Weights", "Abscissas"
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

            CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, tj, &
                                              Emin, Emax, max_error_min, num_points_per_magnitude, &
                                              my_regularization)

            ! get the weights for the cosine transform W^c(iw) -> W^c(it)
            ALLOCATE (weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
            weights_cos_tf_w_to_t = 0.0_dp

            CALL get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, tj, &
                                              Emin, Emax, max_error_min, num_points_per_magnitude, &
                                              my_regularization)

            IF (do_gw_im_time) THEN

               ! get the weights for the sine transform Sigma^sin(it) -> Sigma^sin(iw) (PRB 94, 165109 (2016), Eq. 71)
               ALLOCATE (weights_sin_tf_t_to_w(num_integ_points, num_integ_points))
               weights_sin_tf_t_to_w = 0.0_dp

               CALL get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, tj, &
                                                 Emin, Emax, max_error_min, num_points_per_magnitude, &
                                                 my_regularization)

               IF (unit_nr > 0) THEN
                  WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
                     "MINIMAX_INFO| Maximum deviation of the imag. time fit:", max_error_min
               END IF
            END IF

         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE get_minimax_grid

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param unit_nr ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param occupation_numbers ...
!> \param electronic_temperature ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param reference_grid_size ...
!> \param extrapolation_order ...
!> \param color_rpa_group ...
!> \param fm_mat_S ...
!> \param my_do_gw ...
!> \param ext_scaling ...
!> \param tj ...
!> \param wj ...
! **************************************************************************************************
   SUBROUTINE get_clenshaw_grid(para_env, para_env_RPA, unit_nr, homo, virtual, Eigenval, occupation_numbers, &
                                electronic_temperature, num_integ_points, &
                                num_integ_group, reference_grid_size, extrapolation_order, color_rpa_group, fm_mat_S, my_do_gw, &
                                ext_scaling, tj, wj)

      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_RPA
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Eigenval
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: occupation_numbers
      REAL(KIND=dp), INTENT(IN)                          :: electronic_temperature
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            reference_grid_size
      INTEGER, INTENT(INOUT)                             :: extrapolation_order
      INTEGER, INTENT(IN)                                :: color_rpa_group
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S
      LOGICAL, INTENT(IN)                                :: my_do_gw
      REAL(KIND=dp), INTENT(IN)                          :: ext_scaling
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tj, wj

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_clenshaw_grid'

      INTEGER                                            :: handle, jquad, nspins
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp)                                      :: a_scaling, tmp

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      my_open_shell = (nspins == 2)

      ! Now, start to prepare the different grid
      ALLOCATE (tj(num_integ_points))
      tj = 0.0_dp

      ALLOCATE (wj(num_integ_points))
      wj = 0.0_dp

      DO jquad = 1, num_integ_points - 1
         tmp = jquad*pi/(2.0_dp*num_integ_points)
         tj(jquad) = 1.0_dp/TAN(tmp)
         wj(jquad) = pi/(num_integ_points*SIN(tmp)**2)
      END DO
      tj(num_integ_points) = 0.0_dp
      wj(num_integ_points) = pi/(2.0_dp*num_integ_points)

      IF (my_do_gw .AND. ext_scaling > 0.0_dp) THEN
         a_scaling = ext_scaling
      ELSE
         CALL calc_scaling_factor(a_scaling, para_env, para_env_RPA, homo, virtual, Eigenval, occupation_numbers, &
                                  electronic_temperature, num_integ_points, reference_grid_size, extrapolation_order, &
                                  num_integ_group, color_rpa_group, tj, wj, fm_mat_S)
      END IF

      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.5)') 'INTEG_INFO| Scaling parameter:', a_scaling

      wj(:) = wj(:)*a_scaling
      tj(:) = tj*a_scaling
      IF (electronic_temperature > 0.0_dp) THEN
         wj(:) = wj/(2.0_dp*pi)
      END IF

      IF (unit_nr > 0) CALL m_flush(unit_nr)
      CALL para_env%sync()

      CALL timestop(handle)

   END SUBROUTINE get_clenshaw_grid

! **************************************************************************************************
!> \brief ...
!> \param a_scaling_ext ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param occupation_numbers ...
!> \param electronic_temperature ...
!> \param num_integ_points ...
!> \param reference_grid_size ...
!> \param extrapolation_order ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param tj_ext ...
!> \param wj_ext ...
!> \param fm_mat_S ...
! **************************************************************************************************
   SUBROUTINE calc_scaling_factor(a_scaling_ext, para_env, para_env_RPA, homo, virtual, Eigenval, occupation_numbers, &
                                  electronic_temperature, num_integ_points, reference_grid_size, extrapolation_order, &
                                  num_integ_group, color_rpa_group, &
                                  tj_ext, wj_ext, fm_mat_S)
      REAL(KIND=dp), INTENT(OUT)                         :: a_scaling_ext
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_RPA
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Eigenval
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: occupation_numbers
      REAL(KIND=dp), INTENT(IN)                          :: electronic_temperature
      INTEGER, INTENT(IN)                                :: num_integ_points, reference_grid_size
      INTEGER, INTENT(INOUT)                             :: extrapolation_order
      INTEGER, INTENT(IN)                                :: num_integ_group, color_rpa_group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tj_ext, wj_ext
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_scaling_factor'

      INTEGER                                            :: handle, icycle, nspins
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp)                                      :: a_high, a_low, a_scaling, conv_param, &
                                                            eps, left_term, right_term, &
                                                            right_term_ref, right_term_ref_beta, &
                                                            step
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: D_ia, D_ia_beta, iaia_RI, iaia_RI_beta, &
                                                            M_ia, M_ia_beta
      TYPE(mp_para_env_type), POINTER                    :: para_env_col, para_env_col_beta

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      my_open_shell = (nspins == 2)

      eps = 1.0E-10_dp

      CALL calc_ia_ia_integrals(para_env_RPA, virtual(1), right_term_ref, &
                                Eigenval(:, 1, 1), occupation_numbers(:, 1), reference_grid_size, &
                                extrapolation_order, electronic_temperature, &
                                D_ia, iaia_RI, M_ia, fm_mat_S(1), para_env_col)

      ! In the open shell case do point 1-2-3 for the beta spin
      IF (my_open_shell) THEN
         CALL calc_ia_ia_integrals(para_env_RPA, virtual(2), right_term_ref_beta, &
                                   Eigenval(:, 1, 2), occupation_numbers(:, 2), reference_grid_size, &
                                   extrapolation_order, electronic_temperature, &
                                   D_ia_beta, iaia_RI_beta, M_ia_beta, fm_mat_S(2), para_env_col_beta)

         right_term_ref = right_term_ref + right_term_ref_beta
      END IF

      ! bcast the result
      IF (para_env%mepos == 0) THEN
         CALL para_env%bcast(right_term_ref, 0)
      ELSE
         right_term_ref = 0.0_dp
         CALL para_env%bcast(right_term_ref, 0)
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

         CALL calculate_objfunc(a_scaling, left_term, num_integ_points, my_open_shell, &
                                M_ia, tj_ext, wj_ext, D_ia, D_ia_beta, M_ia_beta, &
                                num_integ_group, color_rpa_group, &
                                para_env, para_env_col, para_env_col_beta)
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

               CALL calculate_objfunc(a_scaling, left_term, num_integ_points, my_open_shell, &
                                      M_ia, tj_ext, wj_ext, D_ia, D_ia_beta, M_ia_beta, &
                                      num_integ_group, color_rpa_group, &
                                      para_env, para_env_col, para_env_col_beta)
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
      CALL para_env%bcast(a_scaling_ext, 0)

      DEALLOCATE (iaia_RI)
      DEALLOCATE (D_ia)
      DEALLOCATE (M_ia)
      CALL mp_para_env_release(para_env_col)

      IF (my_open_shell) THEN
         DEALLOCATE (iaia_RI_beta)
         DEALLOCATE (D_ia_beta)
         DEALLOCATE (M_ia_beta)
         CALL mp_para_env_release(para_env_col_beta)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_scaling_factor

! **************************************************************************************************
!> \brief ...
!> \param para_env_RPA ...
!> \param virtual ...
!> \param right_term_ref ...
!> \param Eigenval ...
!> \param occupation_numbers ...
!> \param reference_grid_size ...
!> \param extrapolation_order ...
!> \param electronic_temperature ...
!> \param D_ia ...
!> \param iaia_RI ...
!> \param M_ia ...
!> \param fm_mat_S ...
!> \param para_env_col ...
! **************************************************************************************************
   SUBROUTINE calc_ia_ia_integrals(para_env_RPA, virtual, right_term_ref, Eigenval, &
                                   occupation_numbers, reference_grid_size, extrapolation_order, electronic_temperature, &
                                   D_ia, iaia_RI, M_ia, fm_mat_S, para_env_col)

      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_RPA
      INTEGER, INTENT(IN)                                :: virtual
      REAL(KIND=dp), INTENT(OUT)                         :: right_term_ref
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval, occupation_numbers
      INTEGER, INTENT(IN)                                :: reference_grid_size
      INTEGER, INTENT(INOUT)                             :: extrapolation_order
      REAL(KIND=dp), INTENT(IN)                          :: electronic_temperature
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: D_ia, iaia_RI, M_ia
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S
      TYPE(mp_para_env_type), POINTER                    :: para_env_col

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ia_ia_integrals'

      INTEGER                                            :: avirt, color_col, color_row, dimen, &
                                                            handle, i_global, iiB, iocc, &
                                                            my_reference_grid_size, ncol_local, &
                                                            nrow_local, point
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: eigen_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tj_matsubara, wj_matsubara
      TYPE(mp_para_env_type), POINTER                    :: para_env_row

      CALL timeset(routineN, handle)

      ! calculate the (ia|ia) RI integrals
      ! ----------------------------------
      ! 1) get info fm_mat_S
      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      dimen = SIZE(Eigenval)

      ! allocate the local buffer of iaia_RI integrals (dp kind)
      ALLOCATE (iaia_RI(ncol_local))
      iaia_RI = 0.0_dp

      ! 2) perform the local multiplication SUM_K (ia|K)*(ia|K)
      IF (electronic_temperature > 0.0_dp) THEN
         DO iiB = 1, ncol_local
            i_global = col_indices(iiB)

            iocc = MAX(1, i_global - 1)/virtual + 1
            avirt = i_global - (iocc - 1)*virtual

            IF (avirt + dimen - virtual > iocc) THEN
               iaia_RI(iiB) = ABS(occupation_numbers(avirt + dimen - virtual) - occupation_numbers(iocc))* &
                              DOT_PRODUCT(fm_mat_S%local_data(:, iiB), fm_mat_S%local_data(:, iiB))
            ELSE
               iaia_RI(iiB) = 0.0_dp
            END IF
         END DO
      ELSE
         DO iiB = 1, ncol_local
            iaia_RI(iiB) = iaia_RI(iiB) + DOT_PRODUCT(fm_mat_S%local_data(:, iiB), fm_mat_S%local_data(:, iiB))
         END DO
      END IF

      ! 3) sum the result with the processes of the RPA_group having the same columns
      !          _______ia______               _
      !         |   |   |   |   |             | |
      !     --> | 1 | 5 | 9 | 13|   SUM -->   | |
      !         |___|__ |___|___|             |_|
      !         |   |   |   |   |             | |
      !     --> | 2 | 6 | 10| 14|   SUM -->   | |
      !       K |___|___|___|___|             |_|   (ia|ia)_RI
      !         |   |   |   |   |             | |
      !     --> | 3 | 7 | 11| 15|   SUM -->   | |
      !         |___|___|___|___|             |_|
      !         |   |   |   |   |             | |
      !     --> | 4 | 8 | 12| 16|   SUM -->   | |
      !         |___|___|___|___|             |_|
      !

      color_col = fm_mat_S%matrix_struct%context%mepos(2)
      ALLOCATE (para_env_col)
      CALL para_env_col%from_split(para_env_RPA, color_col)

      CALL para_env_col%sum(iaia_RI)

      ! 4) calculate the right hand term, D_ia is the matrix containing the
      ! orbital energy differences, M_ia is the diagonal of the full RPA 'excitation'
      ! matrix
      ALLOCATE (D_ia(ncol_local))

      ALLOCATE (M_ia(ncol_local))

      DO iiB = 1, ncol_local
         i_global = col_indices(iiB)

         iocc = MAX(1, i_global - 1)/virtual + 1
         avirt = i_global - (iocc - 1)*virtual
! Skip half of the partially occupied pairs to consider them only once
         IF (avirt + dimen - virtual > iocc) THEN
            eigen_diff = ABS(Eigenval(avirt + dimen - virtual) - Eigenval(iocc))
            D_ia(iiB) = eigen_diff
         ELSE
            D_ia(iiB) = 0.0_dp
         END IF
      END DO

      DO iiB = 1, ncol_local
         M_ia(iiB) = 2.0_dp*D_ia(iiB)*iaia_RI(iiB)
      END DO

      IF (electronic_temperature > 0.0_dp) THEN
         my_reference_grid_size = reference_grid_size
         IF (my_reference_grid_size < 1) THEN
            my_reference_grid_size = 2*CEILING(ABS(Eigenval(dimen) - Eigenval(1))/(2.0_dp*pi)) + 1
         END IF
         CALL get_matsubara_grid(-1, electronic_temperature, RESHAPE(Eigenval, [dimen, 1]), &
                                 my_reference_grid_size, extrapolation_order, &
                                 tj_matsubara, wj_matsubara)
         right_term_ref = wj_matsubara(1)*SUM(f_rpa(2.0_dp*iaia_RI/D_ia), MASK=(D_ia > 0.0_dp))
         DO point = 2, my_reference_grid_size
            right_term_ref = right_term_ref + wj_matsubara(point)*SUM(f_rpa(M_ia/(D_ia**2 + tj_matsubara(point)**2)))
         END DO
         right_term_ref = right_term_ref/2.0_dp
      ELSE
         right_term_ref = SUM(SQRT(M_ia + D_ia**2) - D_ia - iaia_RI)
         right_term_ref = right_term_ref/2.0_dp
      END IF

      ! sum the result with the processes of the RPA_group having the same row
      color_row = fm_mat_S%matrix_struct%context%mepos(1)
      ALLOCATE (para_env_row)
      CALL para_env_row%from_split(para_env_RPA, color_row)

      ! allocate communication array for rows
      CALL para_env_row%sum(right_term_ref)

      CALL mp_para_env_release(para_env_row)

      CALL timestop(handle)

   END SUBROUTINE calc_ia_ia_integrals

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL REAL(KIND=dp) FUNCTION f_rpa(x)
      REAL(KIND=dp), INTENT(IN)                          :: x

      IF (ABS(x) > 0.001_dp) THEN
         f_rpa = LOG(1.0_dp + x) - x
      ELSE
         f_rpa = x**2*(-0.5_dp + x/3.0_dp - x**2/4.0_dp)
      END IF
   END FUNCTION f_rpa

! **************************************************************************************************
!> \brief ...
!> \param a_scaling ...
!> \param left_term ...
!> \param num_integ_points ...
!> \param my_open_shell ...
!> \param M_ia ...
!> \param tj ...
!> \param wj ...
!> \param D_ia ...
!> \param D_ia_beta ...
!> \param M_ia_beta ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param para_env ...
!> \param para_env_col ...
!> \param para_env_col_beta ...
! **************************************************************************************************
   SUBROUTINE calculate_objfunc(a_scaling, left_term, num_integ_points, my_open_shell, &
                                M_ia, tj, wj, D_ia, D_ia_beta, M_ia_beta, &
                                num_integ_group, color_rpa_group, &
                                para_env, para_env_col, para_env_col_beta)
      REAL(KIND=dp), INTENT(IN)                          :: a_scaling
      REAL(KIND=dp), INTENT(INOUT)                       :: left_term
      INTEGER, INTENT(IN)                                :: num_integ_points
      LOGICAL, INTENT(IN)                                :: my_open_shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: M_ia, tj, wj, D_ia, D_ia_beta, M_ia_beta
      INTEGER, INTENT(IN)                                :: num_integ_group, color_rpa_group
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_col
      TYPE(mp_para_env_type), POINTER                    :: para_env_col_beta

      INTEGER                                            :: iiB, jquad
      REAL(KIND=dp)                                      :: left_term_beta, omega

      left_term = 0.0_dp
      left_term_beta = 0.0_dp
      DO jquad = 1, num_integ_points
         ! parallelize over integration points
         IF (MODULO(jquad, num_integ_group) /= color_rpa_group) CYCLE
         omega = a_scaling*tj(jquad)

         DO iiB = 1, SIZE(M_ia)
            ! parallelize over ia elements in the para_env_row group
            IF (MODULO(iiB, para_env_col%num_pe) /= para_env_col%mepos) CYCLE
            ! calculate left_term
! Branch required at finite temperatures
            IF (omega**2 + D_ia(iiB)**2 > 0.0_dp) THEN
               left_term = left_term + wj(jquad)* &
                           (LOG(1.0_dp + M_ia(iiB)/(omega**2 + D_ia(iiB)**2)) - &
                            M_ia(iiB)/(omega**2 + D_ia(iiB)**2))
            END IF
         END DO

         IF (my_open_shell) THEN
            DO iiB = 1, SIZE(M_ia_beta)
               ! parallelize over ia elements in the para_env_row group
               IF (MODULO(iiB, para_env_col_beta%num_pe) /= para_env_col_beta%mepos) CYCLE
               ! calculate left_term
               IF (omega**2 + D_ia_beta(iiB)**2 > 0.0_dp) THEN
                  left_term = left_term + wj(jquad)* &
                              (LOG(1.0_dp + M_ia_beta(iiB)/(omega**2 + D_ia_beta(iiB)**2)) - &
                               M_ia_beta(iiB)/(omega**2 + D_ia_beta(iiB)**2))
               END IF
            END DO
         END IF

      END DO

      ! sum the contribution from all proc, starting form the row group
      CALL para_env%sum(left_term)

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
!> \param regularization ...
! **************************************************************************************************
   SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude, &
                                           regularization)

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
      REAL(KIND=dp), INTENT(IN)                          :: regularization

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_t_to_w'

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            x_values, y_values
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
!               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)*sing_values(jjj) &
                                              /(regularization**2 + sing_values(jjj)**2)
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

      DEALLOCATE (x_values, y_values, mat_A, tau_wj_work, sing_values, mat_U, mat_SinvVSinvT, &
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
!> \param regularization ...
! **************************************************************************************************
   SUBROUTINE get_l_sq_wghts_sin_tf_t_to_w(num_integ_points, tau_tj, weights_sin_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude, regularization)

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
      REAL(KIND=dp), INTENT(IN)                          :: regularization

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_sin_tf_t_to_w'

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
!               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)*sing_values(jjj) &
                                              /(regularization**2 + sing_values(jjj)**2)
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

            CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, &
                                              weights_cos_tf_t_to_w, tj, &
                                              1.0_dp, Rc, max_error, 200, 0.0_dp)

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
!> \param regularization ...
! **************************************************************************************************
   SUBROUTINE get_l_sq_wghts_cos_tf_w_to_t(num_integ_points, tau_tj, weights_cos_tf_w_to_t, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude, regularization)

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
      REAL(KIND=dp), INTENT(IN)                          :: regularization

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_w_to_t'

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
!               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)*sing_values(jjj) &
                                              /(regularization**2 + sing_values(jjj)**2)
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_omega_grid_with_cosine'

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
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), INTENT(OUT)                         :: gap, max_eig_diff, e_fermi

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gap_and_max_eig_diff_kpoints'

      INTEGER                                            :: handle, homo, ikpgr, ispin, kplocal, &
                                                            nmo, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      REAL(KIND=dp)                                      :: e_homo, e_homo_temp, e_lumo, e_lumo_temp
      REAL(KIND=dp), DIMENSION(3)                        :: tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      kpoints=kpoint)

      mo_set => kpoint%kp_env(1)%kpoint_env%mos(1, 1)
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
            mo_set => kp%mos(1, ispin)
            CALL get_mo_set(mo_set, eigenvalues=eigenvalues, homo=homo)
            e_homo_temp = eigenvalues(homo)
            e_lumo_temp = eigenvalues(homo + 1)

            IF (e_homo_temp > e_homo) e_homo = e_homo_temp
            IF (e_lumo_temp < e_lumo) e_lumo = e_lumo_temp
            IF (eigenvalues(nmo) - eigenvalues(1) > max_eig_diff) max_eig_diff = eigenvalues(nmo) - eigenvalues(1)

         END DO
      END DO

      ! Collect all three numbers in an array
      ! Reverse sign of lumo to reduce number of MPI calls
      tmp(1) = e_homo
      tmp(2) = -e_lumo
      tmp(3) = max_eig_diff
      CALL para_env%max(tmp)

      gap = -tmp(2) - tmp(1)
      e_fermi = (tmp(1) - tmp(2))*0.5_dp
      max_eig_diff = tmp(3)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param unit_nr ...
!> \param electronic_temperature ...
!> \param eigenvalues ...
!> \param num_integ_points ...
!> \param extrapolation_order ...
!> \param tj ...
!> \param wj ...
! **************************************************************************************************
   SUBROUTINE get_matsubara_grid(unit_nr, electronic_temperature, eigenvalues, num_integ_points, extrapolation_order, tj, wj)
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), INTENT(IN)                          :: electronic_temperature
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eigenvalues
      INTEGER, INTENT(IN)                                :: num_integ_points
      INTEGER, INTENT(INOUT)                             :: extrapolation_order
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tj, wj

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_matsubara_grid'

      INTEGER                                            :: dimen, handle, j, &
                                                            recommended_number_of_points
      REAL(KIND=dp)                                      :: Emax, n_real, w_0

      CALL timeset(routineN, handle)

      w_0 = 2.0_dp*pi*electronic_temperature

      ALLOCATE (tj(num_integ_points), wj(num_integ_points))

      ! Matrix Q is an even function, so the first summand (omega=0) occurs only once, the others twice (+omega and -omega)
      wj(1) = electronic_temperature
      wj(2:) = 2.0_dp*electronic_temperature

      DO j = 1, num_integ_points
         tj(j) = REAL(j - 1, KIND=dp)*w_0
      END DO

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "MATSUBARA_INFO| Number of non-negative integration points:", num_integ_points
         WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
            "MATSUBARA_INFO| Size of frequency step:", w_0
      END IF

      ! Give hint on extrapolation
      dimen = SIZE(eigenvalues, 1)
      Emax = MAXVAL(eigenvalues(dimen, :) - eigenvalues(1, :))
      recommended_number_of_points = CEILING(Emax/w_0) + 1

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "MATSUBARA_INFO| Recommended number of points:", recommended_number_of_points
      END IF
      IF (recommended_number_of_points > num_integ_points) THEN
         CALL cp_warn(__LOCATION__, &
                      "The number of grid points is supposedly too low at the given "// &
                      "temperature! Check the convergence of the total energy thoroughly.")
      END IF

      IF (extrapolation_order <= 0) extrapolation_order = num_integ_points

      IF (extrapolation_order >= num_integ_points) THEN
         CALL cp_warn(__LOCATION__, &
                      "Not enough points for extrapolation (Hint: the point at the origin is not considered here). "// &
                      "CP2K is setting the extrapolation order to the number of integration points-1!")
         extrapolation_order = num_integ_points - 1
      END IF

      ! We perform a simple extrapolation of the total energy contribution
      ! We hard-code the rules to reduce numerical noise occurring at low temperatures and large number of required grid points
      IF (extrapolation_order >= 3) THEN
         IF (extrapolation_order > 3) THEN
            CALL cp_warn(__LOCATION__, &
                         "Unsupported extrapolation order. Extrapolation order was set to 3!")
            extrapolation_order = 3
         END IF
         n_real = REAL(num_integ_points - 1, dp)
         wj(num_integ_points - 1) = wj(num_integ_points - 1)*(n_real**2 - 2.0_dp*(n_real - 1.0_dp)**2 &
                                                        *(1.0_dp - 1.0_dp/n_real)**2)/(12.0_dp - 24.0_dp/n_real + 14.0_dp/n_real**2)
         wj(num_integ_points) = wj(num_integ_points)*n_real**2/(12.0_dp - 24.0_dp/n_real + 14.0_dp/n_real**2)
      ELSE IF (extrapolation_order == 2) THEN
         n_real = REAL(num_integ_points - 1, dp)
         wj(num_integ_points) = wj(num_integ_points)*(1.0_dp + n_real/(3.0_dp - 3.0_dp/n_real + n_real**(-2)))
      ELSE
         ! Nothing to do
      END IF

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "MATSUBARA_INFO| Extrapolation order:", extrapolation_order
         CALL m_flush(unit_nr)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_matsubara_grid

END MODULE mp2_grids
