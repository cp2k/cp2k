!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exchange and Correlation functional calculations
!> \par History
!>      (13-Feb-2001) JGH, based on earlier version of apsi
!>      02.2003 Many many changes [fawzi]
!>      03.2004 new xc interface [fawzi]
!>      04.2004 kinetic functionals [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE xc
   #:include 'xc.fypp'
   USE cp_array_utils, ONLY: cp_3d_r_p_type
   USE cp_linked_list_xc_deriv, ONLY: cp_sll_xc_deriv_next, &
                                      cp_sll_xc_deriv_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_unit_nr, &
                              cp_logger_type, &
                              cp_to_string
   USE input_section_types, ONLY: section_get_ival, &
                                  section_get_lval, &
                                  section_get_rval, &
                                  section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE kahan_sum, ONLY: accurate_dot_product, &
                        accurate_sum
   USE kinds, ONLY: default_path_length, &
                    dp
   USE message_passing, ONLY: mp_sum
   USE pw_grid_types, ONLY: PW_MODE_DISTRIBUTED, &
                            pw_grid_type
   USE pw_methods, ONLY: pw_axpy, &
                         pw_copy, &
                         pw_derive, &
                         pw_scale, &
                         pw_transfer, &
                         pw_zero, pw_laplace
   USE pw_pool_types, ONLY: pw_pool_create_pw, &
                            pw_pool_give_back_cr3d, &
                            pw_pool_give_back_pw, &
                            pw_pool_type
   USE pw_spline_utils, ONLY: &
      nn10_coeffs, nn10_deriv_coeffs, nn50_coeffs, nn50_deriv_coeffs, pw_nn_deriv_r, &
      pw_nn_smear_r, pw_spline2_deriv_g, pw_spline2_interpolate_values_g, pw_spline3_deriv_g, &
      pw_spline3_interpolate_values_g, pw_spline_scale_deriv, spline2_coeffs, &
      spline2_deriv_coeffs, spline3_coeffs, spline3_deriv_coeffs
   USE pw_types, ONLY: COMPLEXDATA1D, &
                       REALDATA3D, &
                       REALSPACE, &
                       RECIPROCALSPACE, &
                       pw_create, &
                       pw_p_type, &
                       pw_release, &
                       pw_type
   USE xc_derivative_desc, ONLY: MAX_DERIVATIVE_DESC_LENGTH, &
                                 MAX_LABEL_LENGTH
   USE xc_derivative_set_types, ONLY: xc_derivative_set_type, &
                                      xc_dset_create, &
                                      xc_dset_get_derivative, &
                                      xc_dset_release, &
                                      xc_dset_zero_all
   USE xc_derivative_types, ONLY: xc_derivative_get, &
                                  xc_derivative_type
   USE xc_derivatives, ONLY: xc_functionals_eval, &
                             xc_functionals_get_needs
   USE xc_input_constants, ONLY: &
      xc_debug_new_routine, xc_deriv_nn10_smooth, xc_deriv_nn50_smooth, xc_deriv_pw, &
      xc_deriv_spline2, xc_deriv_spline2_smooth, xc_deriv_spline3, xc_deriv_spline3_smooth, &
      xc_new_f_routine, xc_rho_nn10, xc_rho_nn50, xc_rho_no_smooth, xc_rho_spline2_smooth, &
      xc_rho_spline3_smooth, xc_test_lsd_f_routine
   USE xc_rho_cflags_types, ONLY: xc_rho_cflags_type
   USE xc_rho_set_types, ONLY: xc_rho_set_create, &
                               xc_rho_set_get, &
                               xc_rho_set_release, &
                               xc_rho_set_type, &
                               xc_rho_set_update
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: xc_vxc_pw_create1, xc_vxc_pw_create, &
             xc_exc_calc, xc_calc_2nd_deriv_analytical, xc_calc_2nd_deriv_numerical, &
             xc_calc_2nd_deriv, xc_prep_2nd_deriv, divide_by_norm_drho, smooth_cutoff, &
             xc_uses_kinetic_energy_density, xc_uses_norm_drho

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param xc_fun_section ...
!> \param lsd ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_uses_kinetic_energy_density(xc_fun_section, lsd) RESULT(res)
      TYPE(section_vals_type), POINTER, INTENT(IN) :: xc_fun_section
      LOGICAL, INTENT(IN) :: lsd
      LOGICAL :: res

      TYPE(xc_rho_cflags_type)                           :: needs

      needs = xc_functionals_get_needs(xc_fun_section, &
                                       lsd=lsd, &
                                       add_basic_components=.TRUE.)
      res = (needs%tau_spin .OR. needs%tau)

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param xc_fun_section ...
!> \param lsd ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_uses_norm_drho(xc_fun_section, lsd) RESULT(res)
      TYPE(section_vals_type), POINTER, INTENT(IN) :: xc_fun_section
      LOGICAL, INTENT(IN) :: lsd
      LOGICAL :: res

      TYPE(xc_rho_cflags_type)                           :: needs

      needs = xc_functionals_get_needs(xc_fun_section, &
                                       lsd=lsd, &
                                       add_basic_components=.TRUE.)
      res = (needs%norm_drho .OR. needs%norm_drho_spin)

   END FUNCTION

! **************************************************************************************************
!> \brief Exchange and Correlation functional calculations.
!>      depending on the selected functional_routine calls
!>      the correct routine
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xcthe functional
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param exc the xc energy
!> \param xc_section parameters selecting the xc and the method used to
!>        calculate it
!> \param pw_pool the pool for the grids
!> \param compute_virial should the virial be computed... if so virial_xc must be present
!> \param virial_xc for calculating the GGA part of the stress
!> \param exc_r ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE xc_vxc_pw_create1(vxc_rho, vxc_tau, rho_r, rho_g, tau, exc, xc_section, &
                                pw_pool, compute_virial, virial_xc, exc_r)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rho, vxc_tau, rho_r, rho_g, tau
      REAL(KIND=dp), INTENT(out)                         :: exc
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      LOGICAL, INTENT(IN)                                :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: virial_xc
      TYPE(pw_p_type), INTENT(INOUT), OPTIONAL           :: exc_r

      INTEGER                                            :: f_routine

      CPASSERT(ASSOCIATED(rho_r))
      CPASSERT(ASSOCIATED(xc_section))
      CPASSERT(.NOT. ASSOCIATED(vxc_rho))
      CPASSERT(.NOT. ASSOCIATED(vxc_tau))

      virial_xc = 0.0_dp

      CALL section_vals_val_get(xc_section, "FUNCTIONAL_ROUTINE", &
                                i_val=f_routine)
      SELECT CASE (f_routine)
      CASE (xc_new_f_routine)
         CALL xc_vxc_pw_create(vxc_rho=vxc_rho, vxc_tau=vxc_tau, tau=tau, exc_r=exc_r, &
                               rho_r=rho_r, rho_g=rho_g, exc=exc, xc_section=xc_section, &
                               pw_pool=pw_pool, compute_virial=compute_virial, virial_xc=virial_xc)
      CASE (xc_debug_new_routine)
         CPASSERT(.NOT. PRESENT(exc_r))
         CALL xc_vxc_pw_create_debug(vxc_rho=vxc_rho, vxc_tau=vxc_tau, tau=tau, &
                                     rho_r=rho_r, rho_g=rho_g, exc=exc, xc_section=xc_section, &
                                     pw_pool=pw_pool)
      CASE (xc_test_lsd_f_routine)
         CPASSERT(.NOT. PRESENT(exc_r))
         CALL xc_vxc_pw_create_test_lsd(vxc_rho=vxc_rho, vxc_tau=vxc_tau, &
                                        tau=tau, rho_r=rho_r, rho_g=rho_g, exc=exc, &
                                        xc_section=xc_section, pw_pool=pw_pool)
      CASE default
      END SELECT

   END SUBROUTINE xc_vxc_pw_create1

! **************************************************************************************************
!> \brief calculates vxc using lsd with rhoa=rhob=0.5*rho and compares
!>      with the lda result
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xc
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param exc the xc energy
!> \param xc_section which functional to calculate, and how
!> \param pw_pool the pool for the grids
!> \author Fawzi Mohamed
!> \note
!>      for debugging only: leaks, and non parallel
! **************************************************************************************************
   SUBROUTINE xc_vxc_pw_create_test_lsd(vxc_rho, vxc_tau, rho_r, rho_g, tau, &
                                        exc, xc_section, pw_pool)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rho, vxc_tau, rho_r, rho_g, tau
      REAL(KIND=dp), INTENT(out)                         :: exc
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      LOGICAL, PARAMETER                                 :: compute_virial = .FALSE.

      CHARACTER(len=default_path_length)                 :: filename
      CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), &
         POINTER                                         :: split_desc
      INTEGER                                            :: i, ii, ispin, j, k
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(kind=dp)                                      :: diff, exc2, maxdiff, tmp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: pot, pot2, pot3
      TYPE(cp_sll_xc_deriv_type), POINTER                :: deriv_iter
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho2_g, rho2_r, tau2, vxc_rho2, vxc_tau2
      TYPE(xc_derivative_set_type)                       :: dSet1, dSet2
      TYPE(xc_derivative_type), POINTER                  :: deriv, deriv2, deriv3
      TYPE(xc_rho_set_type), POINTER                     :: rho_set1, rho_set2

      NULLIFY (vxc_rho2, vxc_tau2, tau2, rho_set1, rho_set2, split_desc, pot, pot3, pot3, &
               deriv, deriv2, deriv3, rho2_g)

      bo = rho_r(1)%pw%pw_grid%bounds_local

      ALLOCATE (rho2_r(2))
      DO ispin = 1, 2
         NULLIFY (rho2_r(ispin)%pw)
         CALL pw_pool_create_pw(pw_pool, rho2_r(ispin)%pw, in_space=REALSPACE, &
                                use_data=REALDATA3D)
      END DO
      DO k = bo(1, 3), bo(2, 3)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 1), bo(2, 1)
               tmp = rho_r(1)%pw%cr3d(i, j, k)*0.5
               rho2_r(1)%pw%cr3d(i, j, k) = tmp
               rho2_r(2)%pw%cr3d(i, j, k) = tmp
            END DO
         END DO
      END DO

      IF (ASSOCIATED(tau)) THEN
         ALLOCATE (tau2(2))
         DO ispin = 1, 2
            NULLIFY (tau2(ispin)%pw)
            CALL pw_pool_create_pw(pw_pool, tau2(ispin)%pw, in_space=REALSPACE, &
                                   use_data=REALDATA3D)
         END DO

         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  tmp = tau(1)%pw%cr3d(i, j, k)*0.5
                  tau2(1)%pw%cr3d(i, j, k) = tmp
                  tau2(2)%pw%cr3d(i, j, k) = tmp
               END DO
            END DO
         END DO
      END IF

      PRINT *, "about to calculate xc (lda)"
      CALL xc_rho_set_and_dset_create(rho_r=rho_r, rho_g=rho_g, &
                                      tau=tau, xc_section=xc_section, &
                                      pw_pool=pw_pool, rho_set=rho_set1, &
                                      deriv_set=dSet1, deriv_order=1, &
                                      needs_basic_components=.FALSE.)
      CALL xc_vxc_pw_create(rho_r=rho_r, rho_g=rho_g, tau=tau, &
                            vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=exc, xc_section=xc_section, &
                            pw_pool=pw_pool, &
                            compute_virial=compute_virial, virial_xc=virial_xc)
      PRINT *, "did calculate xc (lda)"
      PRINT *, "about to calculate xc (lsd)"
      CALL xc_rho_set_and_dset_create(rho_set=rho_set2, deriv_set=dSet2, &
                                      rho_r=rho2_r, rho_g=rho2_g, tau=tau2, xc_section=xc_section, &
                                      pw_pool=pw_pool, deriv_order=1, &
                                      needs_basic_components=.FALSE.)
      CALL xc_vxc_pw_create(rho_r=rho2_r, rho_g=rho2_g, tau=tau2, &
                            vxc_rho=vxc_rho2, vxc_tau=vxc_tau2, exc=exc2, xc_section=xc_section, &
                            pw_pool=pw_pool, compute_virial=compute_virial, virial_xc=virial_xc)
      PRINT *, "did calculate xc (new)"
      PRINT *, "at (0,0,0) rho_r=", rho_r(1)%pw%cr3d(0, 0, 0), &
         "rho2_r(1)=", rho2_r(1)%pw%cr3d(0, 0, 0), &
         "rho2_r(2)=", rho2_r(2)%pw%cr3d(0, 0, 0), &
         "rho_r_sm=", rho_set1%rho(0, 0, 0), "rhoa2_r_sm=", rho_set2%rhoa(0, 0, 0), &
         "rhob2_r_sm=", rho_set2%rhob(0, 0, 0)
      OPEN (unit=120, file="rho.bindata", status="unknown", access='sequential', &
            form="unformatted", action="write")
      pot => rho_set1%rho
      WRITE (unit=120) pot(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(2, 3))
      CLOSE (unit=120)
      OPEN (unit=120, file="rhoa.bindata", status="unknown", access='sequential', &
            form="unformatted", action="write")
      pot => rho_set2%rhoa
      WRITE (unit=120) pot(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(2, 3))
      CLOSE (unit=120)
      OPEN (unit=120, file="rhob.bindata", status="unknown", access='sequential', &
            form="unformatted", action="write")
      pot => rho_set2%rhob
      WRITE (unit=120) pot(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(2, 3))
      CLOSE (unit=120)
      OPEN (unit=120, file="ndrho.bindata", status="unknown", access='sequential', &
            form="unformatted", action="write")
      pot => rho_set1%norm_drho
      WRITE (unit=120) pot(:, :, bo(2, 3))
      CLOSE (unit=120)
      OPEN (unit=120, file="ndrhoa.bindata", status="unknown", access='sequential', &
            form="unformatted", action="write")
      pot => rho_set2%norm_drhoa
      WRITE (unit=120) pot(:, :, bo(2, 3))
      CLOSE (unit=120)
      OPEN (unit=120, file="ndrhob.bindata", status="unknown", access='sequential', &
            form="unformatted", action="write")
      pot => rho_set2%norm_drhob
      WRITE (unit=120) pot(:, :, bo(2, 3))
      CLOSE (unit=120)
      IF (rho_set1%has%tau) THEN
         OPEN (unit=120, file="tau.bindata", status="unknown", access='sequential', &
               form="unformatted", action="write")
         pot => rho_set1%tau
         WRITE (unit=120) pot(:, :, bo(2, 3))
         CLOSE (unit=120)
      END IF
      IF (rho_set2%has%tau_spin) THEN
         OPEN (unit=120, file="tau_a.bindata", status="unknown", access='sequential', &
               form="unformatted", action="write")
         pot => rho_set2%tau_a
         WRITE (unit=120) pot(:, :, bo(2, 3))
         CLOSE (unit=120)
         OPEN (unit=120, file="tau_v.bindata", status="unknown", access='sequential', &
               form="unformatted", action="write")
         pot => rho_set2%tau_b
         WRITE (unit=120) pot(:, :, bo(2, 3))
         CLOSE (unit=120)
      END IF
      OPEN (unit=120, file="vxc.bindata", status="unknown", access='sequential', &
            form="unformatted", action="write")
      pot => vxc_rho(1)%pw%cr3d
      WRITE (unit=120) pot(:, :, bo(2, 3))
      CLOSE (unit=120)
      OPEN (unit=120, file="vxc2.bindata", status="unknown", access='sequential', &
            form="unformatted", action="write")
      pot => vxc_rho2(1)%pw%cr3d
      WRITE (unit=120) pot(:, :, bo(2, 3))
      CLOSE (unit=120)
      IF (ASSOCIATED(vxc_tau)) THEN
         OPEN (unit=120, file="vxc_tau.bindata", status="unknown", access='sequential', &
               form="unformatted", action="write")
         pot => vxc_tau(1)%pw%cr3d
         WRITE (unit=120) pot(:, :, bo(2, 3))
         CLOSE (unit=120)
      END IF
      IF (ASSOCIATED(vxc_tau2)) THEN
         OPEN (unit=120, file="vxc_tau2_a.bindata", status="unknown", access='sequential', &
               form="unformatted", action="write")
         pot => vxc_tau2(1)%pw%cr3d
         WRITE (unit=120) pot(:, :, bo(2, 3))
         CLOSE (unit=120)
         OPEN (unit=120, file="vxc_tau2_b.bindata", status="unknown", access='sequential', &
               form="unformatted", action="write")
         pot => vxc_tau2(2)%pw%cr3d
         WRITE (unit=120) pot(:, :, bo(2, 3))
         CLOSE (unit=120)
      END IF

      PRINT *, "calc diff on vxc"
      maxDiff = 0.0_dp
      DO ispin = 1, 1
         ii = 0
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  ii = ii + 1
                  diff = ABS(vxc_rho(ispin)%pw%cr3d(i, j, k) - &
                             vxc_rho2(ispin)%pw%cr3d(i, j, k))
                  IF (ii == 1) THEN
                     PRINT *, "vxc", ispin, "=", vxc_rho(ispin)%pw%cr3d(i, j, k), &
                        "vs", vxc_rho2(ispin)%pw%cr3d(i, j, k), "diff=", diff
                  END IF
                  IF (maxDiff < diff) THEN
                     maxDiff = diff
                     PRINT *, "diff=", diff, " at ", i, ",", j, ",", k, &
                        " spin=", ispin, "rho=", rho_set1%rho(i, j, k), &
                        " ndrho=", rho_set1%norm_drho(i, j, k)
                  END IF
               END DO
            END DO
         END DO
      END DO
      PRINT *, "diff exc=", ABS(exc - exc2), "diff vxc=", maxdiff
      ! CPASSERT(maxdiff<5.e-11)
      ! CPASSERT(ABS(exc-exc2)<1.e-14)

      IF (ASSOCIATED(vxc_tau)) THEN
         PRINT *, "calc diff on vxc_tau"
         maxDiff = 0.0_dp
         DO ispin = 1, 1
            ii = 0
            DO k = bo(1, 3), bo(2, 3)
               DO j = bo(1, 2), bo(2, 2)
                  DO i = bo(1, 1), bo(2, 1)
                     ii = ii + 1
                     diff = ABS(vxc_tau(ispin)%pw%cr3d(i, j, k) - &
                                vxc_tau2(ispin)%pw%cr3d(i, j, k))
                     IF (ii == 1) THEN
                        PRINT *, "vxc_tau", ispin, "=", vxc_tau(ispin)%pw%cr3d(i, j, k), &
                           "vs", vxc_tau2(ispin)%pw%cr3d(i, j, k), "diff=", diff
                     END IF
                     IF (maxDiff < diff) THEN
                        maxDiff = diff
                        PRINT *, "diff=", diff, " at ", i, ",", j, ",", k, &
                           " spin=", ispin, "rho=", rho_set1%rho(i, j, k), &
                           " ndrho=", rho_set1%norm_drho(i, j, k)
                     END IF
                  END DO
               END DO
            END DO
         END DO
         PRINT *, "diff exc=", ABS(exc - exc2), "diff vxc_tau=", maxdiff
      END IF
      deriv_iter => dSet1%derivs
      DO WHILE (cp_sll_xc_deriv_next(deriv_iter, el_att=deriv))
         CALL xc_derivative_get(deriv, &
                                split_desc=split_desc, deriv_data=pot)
         SELECT CASE (SIZE(split_desc))
         CASE (0)
            filename = "e_0.bindata"
            deriv2 => xc_dset_get_derivative(dSet2, "")
         CASE (1)
            filename = "e_"//TRIM(split_desc(1))//".bindata"
            IF (split_desc(1) == "rho") THEN
               deriv2 => xc_dset_get_derivative(dSet2, "(rhoa)")
            ELSEIF (split_desc(1) == "tau") THEN
               deriv2 => xc_dset_get_derivative(dSet2, "(tau_a)")
            ELSEIF (split_desc(1) == "norm_drho") THEN
               deriv2 => xc_dset_get_derivative(dSet2, "(norm_drhoa)")
               deriv3 => xc_dset_get_derivative(dSet2, "(norm_drho)")
               IF (ASSOCIATED(deriv3)) THEN
                  IF (ASSOCIATED(deriv2)) THEN
                     CALL xc_derivative_get(deriv2, &
                                            deriv_data=pot2)
                     CALL xc_derivative_get(deriv3, &
                                            deriv_data=pot3)
                     pot2 = pot2 + pot3
                  ELSE
                     deriv2 => deriv3
                  END IF
                  NULLIFY (deriv3, pot2, pot3)
               END IF
            ELSE
               CPABORT("")
            END IF
         CASE default
            CPABORT("")
         END SELECT
         CALL xc_derivative_get(deriv2, &
                                deriv_data=pot2)
         PRINT *, "checking ", filename
         maxDiff = 0.0_dp
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  diff = ABS(pot(i, j, k) - pot2(i, j, k))
                  IF (maxDiff < diff) THEN
                     maxDiff = diff
                     PRINT *, "ediff(", i, j, k, ")=", maxDiff, &
                        "rho=", rho_set1%rho(i, j, k), &
                        "ndrho=", rho_set1%norm_drho(i, j, k)
                  END IF
               END DO
            END DO
         END DO
         PRINT *, "maxdiff ", filename, "=", maxDiff
         OPEN (unit=120, file=TRIM(filename), status="unknown", &
               access='sequential', &
               form="unformatted")
         WRITE (unit=120) pot(:, :, bo(2, 3))
         CLOSE (unit=120)
      END DO
      deriv_iter => dSet2%derivs
      DO WHILE (cp_sll_xc_deriv_next(deriv_iter, el_att=deriv))
         CALL xc_derivative_get(deriv, &
                                split_desc=split_desc, deriv_data=pot)
         SELECT CASE (SIZE(split_desc))
         CASE (0)
            filename = "e_0-2.bindata"
         CASE (1)
            filename = "e_"//TRIM(split_desc(1))//"-2.bindata"
         CASE default
            CPABORT("")
         END SELECT
         OPEN (unit=120, file=TRIM(filename), status="unknown", &
               access='sequential', &
               form="unformatted")
         WRITE (unit=120) pot(:, :, bo(2, 3))
         CLOSE (unit=120)
      END DO
      CALL xc_rho_set_release(rho_set1)
      CALL xc_rho_set_release(rho_set2)
      CALL xc_dset_release(dSet2)
      CALL xc_dset_release(dSet1)
      DO ispin = 1, 2
         CALL pw_pool_give_back_pw(pw_pool, rho2_r(ispin)%pw)
         CALL pw_pool_give_back_pw(pw_pool, vxc_rho2(ispin)%pw)
         IF (ASSOCIATED(vxc_tau2)) THEN
            CALL pw_pool_give_back_pw(pw_pool, vxc_tau2(ispin)%pw)
         END IF
      END DO
      DEALLOCATE (vxc_rho2, rho2_r, rho2_g)
      IF (ASSOCIATED(vxc_tau2)) THEN
         DEALLOCATE (vxc_tau2)
      END IF

   END SUBROUTINE xc_vxc_pw_create_test_lsd

! **************************************************************************************************
!> \brief calculates vxc outputting the yz plane of rho, and of the various components
!>      of the the derivatives and of vxc
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xc
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param exc the xc energy
!> \param xc_section which functional should be used, and how to do it
!> \param pw_pool the pool for the grids
!> \author Fawzi Mohamed
!> \note
!>      for debugging only.
! **************************************************************************************************
   SUBROUTINE xc_vxc_pw_create_debug(vxc_rho, vxc_tau, rho_r, rho_g, tau, exc, &
                                     xc_section, pw_pool)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rho, vxc_tau, rho_r, rho_g, tau
      REAL(KIND=dp), INTENT(out)                         :: exc
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      LOGICAL, PARAMETER                                 :: compute_virial = .FALSE.

      CHARACTER(len=default_path_length)                 :: filename
      CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), &
         POINTER                                         :: split_desc
      INTEGER                                            :: i, ispin, j, k
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_xc
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: pot
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_sll_xc_deriv_type), POINTER                :: deriv_iter
      TYPE(xc_derivative_set_type)                       :: dSet1
      TYPE(xc_derivative_type), POINTER                  :: deriv
      TYPE(xc_rho_set_type), POINTER                     :: rho_set1

      NULLIFY (rho_set1, split_desc, pot, deriv)
      logger => cp_get_default_logger()

      bo = rho_r(1)%pw%pw_grid%bounds_local

      CALL xc_rho_set_and_dset_create(rho_r=rho_r, rho_g=rho_g, &
                                      tau=tau, xc_section=xc_section, &
                                      pw_pool=pw_pool, rho_set=rho_set1, &
                                      deriv_set=dSet1, deriv_order=1, &
                                      needs_basic_components=.FALSE.)

      ! outputs 0,:,: plane
      IF (bo(1, 1) <= 0 .AND. 0 <= bo(2, 1)) THEN
         IF (rho_set1%has%rho_spin) THEN
            OPEN (unit=120, file="rhoa.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%rhoa
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
            OPEN (unit=120, file="rhob.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%rhob
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
         END IF
         IF (rho_set1%has%norm_drho) THEN
            OPEN (unit=120, file="ndrho.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%norm_drho
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
         END IF
         IF (rho_set1%has%norm_drho_spin) THEN
            OPEN (unit=120, file="ndrhoa.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%norm_drhoa
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
            OPEN (unit=120, file="ndrhob.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%norm_drhob
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
         END IF
         IF (rho_set1%has%rho) THEN
            OPEN (unit=120, file="rho.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%rho
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
         END IF
         IF (rho_set1%has%tau) THEN
            OPEN (unit=120, file="tau.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%tau
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
         END IF
         IF (rho_set1%has%tau_spin) THEN
            OPEN (unit=120, file="tau_a.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%tau_a
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
            OPEN (unit=120, file="tau_b.bindata", status="unknown", access='sequential', &
                  form="unformatted", action="write")
            pot => rho_set1%tau_b
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
         END IF

         deriv_iter => dSet1%derivs
         DO WHILE (cp_sll_xc_deriv_next(deriv_iter, el_att=deriv))
            CALL xc_derivative_get(deriv, &
                                   split_desc=split_desc, deriv_data=pot)
            SELECT CASE (SIZE(split_desc))
            CASE (0)
               filename = "e_0.bindata"
            CASE (1)
               filename = "e_"//TRIM(split_desc(1))//".bindata"
            CASE default
               CPABORT("")
            END SELECT
            OPEN (unit=120, file=TRIM(filename), status="unknown", &
                  access='sequential', &
                  form="unformatted")
            WRITE (unit=120) pot(0, :, :)
            CLOSE (unit=120)
         END DO
      END IF

      CALL xc_vxc_pw_create(vxc_rho=vxc_rho, vxc_tau=vxc_tau, &
                            rho_r=rho_r, rho_g=rho_g, tau=tau, &
                            exc=exc, xc_section=xc_section, &
                            pw_pool=pw_pool, &
                            compute_virial=compute_virial, virial_xc=virial_xc)

      ! outputs 0,:,: plane
      IF (bo(1, 1) <= 0 .AND. 0 <= bo(2, 1)) THEN
         IF (ASSOCIATED(vxc_rho)) THEN
            DO ispin = 1, SIZE(vxc_rho)
               WRITE (filename, "('vxc-',i1,'.bindata')") ispin
               OPEN (unit=120, file=filename, status="unknown", access='sequential', &
                     form="unformatted", action="write")
               pot => vxc_rho(ispin)%pw%cr3d
               WRITE (unit=120) pot(0, :, :)
               CLOSE (unit=120)

               pot => vxc_rho(ispin)%pw%cr3d
               DO k = bo(1, 3), bo(2, 3)
                  DO j = bo(1, 2), bo(2, 2)
                     DO i = bo(1, 1), bo(2, 1)
                        IF (ABS(pot(i, j, k)) > 10.0_dp) THEN
                           WRITE (cp_logger_get_default_unit_nr(logger), &
                                  "(' vxc_',i1,'(',i6,',',i6,',',i6,')=',e11.4)", &
                                  advance="no") ispin, i, j, k, pot(i, j, k)
                           IF (rho_set1%has%rho_spin) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' rho=(',e11.4,',',e11.4,')')", advance="no") &
                                 rho_set1%rhoa(i, j, k), rho_set1%rhob(i, j, k)
                           ELSE IF (rho_set1%has%rho) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' rho=',e11.4)", advance="no") rho_set1%rho(i, j, k)
                           END IF
                           IF (rho_set1%has%norm_drho_spin) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' ndrho=(',e11.4,',',e11.4,')')", advance="no") &
                                 rho_set1%norm_drhoa(i, j, k), &
                                 rho_set1%norm_drhob(i, j, k)
                           ELSE IF (rho_set1%has%norm_drho) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' ndrho=',e11.4)", advance="no") rho_set1%norm_drho(i, j, k)
                           END IF
                           IF (rho_set1%has%tau_spin) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' tau=(',e11.4,',',e11.4,')')", advance="no") &
                                 rho_set1%tau_a(i, j, k), &
                                 rho_set1%tau_b(i, j, k)
                           ELSE IF (rho_set1%has%tau) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' tau=',e11.4)", advance="no") rho_set1%tau(i, j, k)
                           END IF

                           deriv_iter => dSet1%derivs
                           DO WHILE (cp_sll_xc_deriv_next(deriv_iter, el_att=deriv))
                              CALL xc_derivative_get(deriv, &
                                                     split_desc=split_desc, deriv_data=pot)
                              SELECT CASE (SIZE(split_desc))
                              CASE (0)
                                 filename = " e_0"
                              CASE (1)
                                 filename = " e_"//TRIM(split_desc(1))
                              CASE default
                                 CPABORT("")
                              END SELECT
                              WRITE (unit=cp_logger_get_default_unit_nr(logger), &
                                     fmt="(a,'=',e11.4)", advance="no") &
                                 TRIM(filename), pot(i, j, k)
                           END DO

                           WRITE (cp_logger_get_default_unit_nr(logger), &
                                  "()")
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF
         IF (ASSOCIATED(vxc_tau)) THEN
            DO ispin = 1, SIZE(vxc_tau)
               WRITE (filename, "('vxc_tau_',i1,'.bindata')") ispin
               OPEN (unit=120, file=filename, status="unknown", access='sequential', &
                     form="unformatted", action="write")
               pot => vxc_tau(ispin)%pw%cr3d
               WRITE (unit=120) pot(0, :, :)
               CLOSE (unit=120)

               pot => vxc_tau(ispin)%pw%cr3d
               DO k = bo(1, 3), bo(2, 3)
                  DO j = bo(1, 2), bo(2, 2)
                     DO i = bo(1, 1), bo(2, 1)
                        IF (ABS(pot(i, j, k)) > 10.0_dp) THEN
                           WRITE (cp_logger_get_default_unit_nr(logger), &
                                  "('vxc_tau_',i1,'(',i6,',',i6,',',i6,')=',e11.4)", &
                                  advance="no") ispin, i, j, k, pot(i, j, k)
                           IF (rho_set1%has%rho_spin) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' rho=(',e11.4,',',e11.4,')')", advance="no") &
                                 rho_set1%rhoa(i, j, k), rho_set1%rhob(i, j, k)
                           ELSE IF (rho_set1%has%rho) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' rho=',e11.4)", advance="no") rho_set1%rho(i, j, k)
                           END IF
                           IF (rho_set1%has%norm_drho_spin) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' ndrho=(',e11.4,',',e11.4,')')", advance="no") &
                                 rho_set1%norm_drhoa(i, j, k), &
                                 rho_set1%norm_drhob(i, j, k)
                           ELSE IF (rho_set1%has%norm_drho) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' ndrho=',e11.4)", advance="no") rho_set1%norm_drho(i, j, k)
                           END IF
                           IF (rho_set1%has%tau_spin) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' tau=(',e11.4,',',e11.4,')')", advance="no") &
                                 rho_set1%tau_a(i, j, k), &
                                 rho_set1%tau_b(i, j, k)
                           ELSE IF (rho_set1%has%tau) THEN
                              WRITE (cp_logger_get_default_unit_nr(logger), &
                                     "(' tau=',e11.4)", advance="no") rho_set1%tau(i, j, k)
                           END IF

                           deriv_iter => dSet1%derivs
                           DO WHILE (cp_sll_xc_deriv_next(deriv_iter, el_att=deriv))
                              CALL xc_derivative_get(deriv, &
                                                     split_desc=split_desc, deriv_data=pot)
                              SELECT CASE (SIZE(split_desc))
                              CASE (0)
                                 filename = " e_0"
                              CASE (1)
                                 filename = " e_"//TRIM(split_desc(1))
                              CASE default
                                 CPABORT("")
                              END SELECT
                              WRITE (unit=cp_logger_get_default_unit_nr(logger), &
                                     fmt="(a,'=',e11.4)", advance="no") &
                                 TRIM(filename), pot(i, j, k)
                           END DO

                           WRITE (cp_logger_get_default_unit_nr(logger), "()")
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END IF

      CALL xc_dset_release(dSet1)
      CALL xc_rho_set_release(rho_set1)
   END SUBROUTINE xc_vxc_pw_create_debug

! **************************************************************************************************
!> \brief creates a xc_rho_set and a derivative set containing the derivatives
!>      of the functionals with the given deriv_order.
!> \param rho_set will contain the rho set
!> \param deriv_set will contain the derivatives
!> \param deriv_order the order of the requested derivatives. If positive
!>        0:deriv_order are calculated, if negative only -deriv_order is
!>        guaranteed to be valid. Orders not requested might be present,
!>        but might contain garbage.
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (can be null, used only
!>        without smoothing of rho or deriv)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param xc_section the section describing the functional to use
!> \param pw_pool the pool for the grids
!> \param needs_basic_components if the basic components of the arguments
!>        should be kept in rho set (a basic component is for example drho
!>        when with lda a functional needs norm_drho)
!> \author fawzi
!> \note
!>      if any of the functionals is gradient corrected the full gradient is
!>      added to the rho set
! **************************************************************************************************
   SUBROUTINE xc_rho_set_and_dset_create(rho_set, deriv_set, deriv_order, &
                                         rho_r, rho_g, tau, xc_section, pw_pool, &
                                         needs_basic_components)

      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type)                       :: deriv_set
      INTEGER, INTENT(in)                                :: deriv_order
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, rho_g, tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      LOGICAL, INTENT(in)                                :: needs_basic_components

      CHARACTER(len=*), PARAMETER :: routineN = 'xc_rho_set_and_dset_create'

      INTEGER                                            :: handle, nspins
      LOGICAL                                            :: lsd
      TYPE(section_vals_type), POINTER                   :: xc_fun_sections

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(rho_set))
      CPASSERT(ASSOCIATED(pw_pool))

      nspins = SIZE(rho_r)
      lsd = (nspins /= 1)

      xc_fun_sections => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")

      CALL xc_dset_create(deriv_set, pw_pool)

      CALL xc_rho_set_create(rho_set, &
                             rho_r(1)%pw%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "density_cutoff"), &
                             drho_cutoff=section_get_rval(xc_section, "gradient_cutoff"), &
                             tau_cutoff=section_get_rval(xc_section, "tau_cutoff"))

      CALL xc_rho_set_update(rho_set, rho_r, rho_g, tau, &
                             xc_functionals_get_needs(xc_fun_sections, lsd, needs_basic_components), &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             pw_pool)

      CALL xc_functionals_eval(xc_fun_sections, &
                               lsd=lsd, &
                               rho_set=rho_set, &
                               deriv_set=deriv_set, &
                               deriv_order=deriv_order)

      CALL timestop(handle)

   END SUBROUTINE xc_rho_set_and_dset_create

! **************************************************************************************************
!> \brief smooths the cutoff on rho with a function smooth(rho) that is 0
!>      for rho<rho_cutoff and 1 for rho>rho_cutoff*rho_smooth_cutoff_range:
!>      E= integral e_0*smooth(rho) => dE/d...= de/d... * smooth,
!>      dE/drho = de/drho * smooth + e_0 * dsmooth/drho
!> \param pot the potential to smooth
!> \param rho , rhoa,rhob: the value of the density (used to apply the cutoff)
!> \param rhoa ...
!> \param rhob ...
!> \param rho_cutoff the value at whch the cutoff function must go to 0
!> \param rho_smooth_cutoff_range range of the smoothing
!> \param e_0 value of e_0, if given it is assumed that pot is the derivative
!>        wrt. to rho, and needs the dsmooth*e_0 contribution
!> \param e_0_scale_factor ...
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE smooth_cutoff(pot, rho, rhoa, rhob, rho_cutoff, &
                            rho_smooth_cutoff_range, e_0, e_0_scale_factor)
      REAL(kind=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: pot, rho, rhoa, rhob
      REAL(kind=dp), INTENT(in)                          :: rho_cutoff, rho_smooth_cutoff_range
      REAL(kind=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: e_0
      REAL(kind=dp), INTENT(in), OPTIONAL                :: e_0_scale_factor

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(kind=dp) :: my_e_0_scale_factor, my_rho, my_rho_n, my_rho_n2, rho_smooth_cutoff, &
                       rho_smooth_cutoff_2, rho_smooth_cutoff_range_2

      CPASSERT(ASSOCIATED(pot))
      bo(1, :) = LBOUND(pot)
      bo(2, :) = UBOUND(pot)
      my_e_0_scale_factor = 1.0_dp
      IF (PRESENT(e_0_scale_factor)) my_e_0_scale_factor = e_0_scale_factor
      rho_smooth_cutoff = rho_cutoff*rho_smooth_cutoff_range
      rho_smooth_cutoff_2 = (rho_cutoff + rho_smooth_cutoff)/2
      rho_smooth_cutoff_range_2 = rho_smooth_cutoff_2 - rho_cutoff

      IF (rho_smooth_cutoff_range > 0.0_dp) THEN
         IF (PRESENT(e_0)) THEN
            CPASSERT(ASSOCIATED(e_0))
            IF (ASSOCIATED(rho)) THEN
!$OMP              PARALLEL DO DEFAULT(NONE) SHARED(bo,e_0,pot,rho,&
!$OMP                          rho_cutoff,rho_smooth_cutoff,rho_smooth_cutoff_2,&
!$OMP                          rho_smooth_cutoff_range_2,my_e_0_scale_factor)&
!$OMP                    PRIVATE(k,j,i,my_rho,my_rho_n,my_rho_n2) COLLAPSE(3)
               DO k = bo(1, 3), bo(2, 3)
                  DO j = bo(1, 2), bo(2, 2)
                     DO i = bo(1, 1), bo(2, 1)
                        my_rho = rho(i, j, k)
                        IF (my_rho < rho_smooth_cutoff) THEN
                           IF (my_rho < rho_cutoff) THEN
                              pot(i, j, k) = 0.0_dp
                           ELSEIF (my_rho < rho_smooth_cutoff_2) THEN
                              my_rho_n = (my_rho - rho_cutoff)/rho_smooth_cutoff_range_2
                              my_rho_n2 = my_rho_n*my_rho_n
                              pot(i, j, k) = pot(i, j, k)* &
                                             my_rho_n2*(my_rho_n - 0.5_dp*my_rho_n2) + &
                                             my_e_0_scale_factor*e_0(i, j, k)* &
                                             my_rho_n2*(3.0_dp - 2.0_dp*my_rho_n) &
                                             /rho_smooth_cutoff_range_2
                           ELSE
                              my_rho_n = 2.0_dp - (my_rho - rho_cutoff)/rho_smooth_cutoff_range_2
                              my_rho_n2 = my_rho_n*my_rho_n
                              pot(i, j, k) = pot(i, j, k)* &
                                             (1.0_dp - my_rho_n2*(my_rho_n - 0.5_dp*my_rho_n2)) &
                                             + my_e_0_scale_factor*e_0(i, j, k)* &
                                             my_rho_n2*(3.0_dp - 2.0_dp*my_rho_n) &
                                             /rho_smooth_cutoff_range_2
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            ELSE
!$OMP              PARALLEL DO DEFAULT(NONE) SHARED(bo,pot,e_0,rhoa,rhob,&
!$OMP                          rho_cutoff,rho_smooth_cutoff,rho_smooth_cutoff_2,&
!$OMP                          rho_smooth_cutoff_range_2,my_e_0_scale_factor)&
!$OMP                    PRIVATE(k,j,i,my_rho,my_rho_n,my_rho_n2) COLLAPSE(3)
               DO k = bo(1, 3), bo(2, 3)
                  DO j = bo(1, 2), bo(2, 2)
                     DO i = bo(1, 1), bo(2, 1)
                        my_rho = rhoa(i, j, k) + rhob(i, j, k)
                        IF (my_rho < rho_smooth_cutoff) THEN
                           IF (my_rho < rho_cutoff) THEN
                              pot(i, j, k) = 0.0_dp
                           ELSEIF (my_rho < rho_smooth_cutoff_2) THEN
                              my_rho_n = (my_rho - rho_cutoff)/rho_smooth_cutoff_range_2
                              my_rho_n2 = my_rho_n*my_rho_n
                              pot(i, j, k) = pot(i, j, k)* &
                                             my_rho_n2*(my_rho_n - 0.5_dp*my_rho_n2) + &
                                             my_e_0_scale_factor*e_0(i, j, k)* &
                                             my_rho_n2*(3.0_dp - 2.0_dp*my_rho_n) &
                                             /rho_smooth_cutoff_range_2
                           ELSE
                              my_rho_n = 2.0_dp - (my_rho - rho_cutoff)/rho_smooth_cutoff_range_2
                              my_rho_n2 = my_rho_n*my_rho_n
                              pot(i, j, k) = pot(i, j, k)* &
                                             (1.0_dp - my_rho_n2*(my_rho_n - 0.5_dp*my_rho_n2)) &
                                             + my_e_0_scale_factor*e_0(i, j, k)* &
                                             my_rho_n2*(3.0_dp - 2.0_dp*my_rho_n) &
                                             /rho_smooth_cutoff_range_2
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END IF
         ELSE
            IF (ASSOCIATED(rho)) THEN
!$OMP              PARALLEL DO DEFAULT(NONE) SHARED(bo,pot,&
!$OMP                          rho_cutoff,rho_smooth_cutoff,rho_smooth_cutoff_2,&
!$OMP                          rho_smooth_cutoff_range_2,rho)&
!$OMP                    PRIVATE(k,j,i,my_rho,my_rho_n,my_rho_n2) COLLAPSE(3)
               DO k = bo(1, 3), bo(2, 3)
                  DO j = bo(1, 2), bo(2, 2)
                     DO i = bo(1, 1), bo(2, 1)
                        my_rho = rho(i, j, k)
                        IF (my_rho < rho_smooth_cutoff) THEN
                           IF (my_rho < rho_cutoff) THEN
                              pot(i, j, k) = 0.0_dp
                           ELSEIF (my_rho < rho_smooth_cutoff_2) THEN
                              my_rho_n = (my_rho - rho_cutoff)/rho_smooth_cutoff_range_2
                              my_rho_n2 = my_rho_n*my_rho_n
                              pot(i, j, k) = pot(i, j, k)* &
                                             my_rho_n2*(my_rho_n - 0.5_dp*my_rho_n2)
                           ELSE
                              my_rho_n = 2.0_dp - (my_rho - rho_cutoff)/rho_smooth_cutoff_range_2
                              my_rho_n2 = my_rho_n*my_rho_n
                              pot(i, j, k) = pot(i, j, k)* &
                                             (1.0_dp - my_rho_n2*(my_rho_n - 0.5_dp*my_rho_n2))
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            ELSE
!$OMP              PARALLEL DO DEFAULT(NONE) SHARED(bo,pot,&
!$OMP                          rho_cutoff,rho_smooth_cutoff,rho_smooth_cutoff_2,&
!$OMP                          rho_smooth_cutoff_range_2,rhoa,rhob)&
!$OMP                    PRIVATE(k,j,i,my_rho,my_rho_n,my_rho_n2) COLLAPSE(3)
               DO k = bo(1, 3), bo(2, 3)
                  DO j = bo(1, 2), bo(2, 2)
                     DO i = bo(1, 1), bo(2, 1)
                        my_rho = rhoa(i, j, k) + rhob(i, j, k)
                        IF (my_rho < rho_smooth_cutoff) THEN
                           IF (my_rho < rho_cutoff) THEN
                              pot(i, j, k) = 0.0_dp
                           ELSEIF (my_rho < rho_smooth_cutoff_2) THEN
                              my_rho_n = (my_rho - rho_cutoff)/rho_smooth_cutoff_range_2
                              my_rho_n2 = my_rho_n*my_rho_n
                              pot(i, j, k) = pot(i, j, k)* &
                                             my_rho_n2*(my_rho_n - 0.5_dp*my_rho_n2)
                           ELSE
                              my_rho_n = 2.0_dp - (my_rho - rho_cutoff)/rho_smooth_cutoff_range_2
                              my_rho_n2 = my_rho_n*my_rho_n
                              pot(i, j, k) = pot(i, j, k)* &
                                             (1.0_dp - my_rho_n2*(my_rho_n - 0.5_dp*my_rho_n2))
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END IF
         END IF
      END IF
   END SUBROUTINE smooth_cutoff

! **************************************************************************************************
!> \brief Exchange and Correlation functional calculations
!> \param vxc_rho will contain the v_xc part that depend on rho
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param vxc_tau will contain the kinetic (tau) part of v_xc
!>        (if one of the chosen xc functionals has it it is allocated and you
!>        are responsible for it)
!> \param exc the xc energy
!> \param rho_r the value of the density in the real space
!> \param rho_g value of the density in the g space (needs to be associated
!>        only for gradient corrections)
!> \param tau value of the kinetic density tau on the grid (can be null,
!>        used only with meta functionals)
!> \param xc_section which functional to calculate, and how to do it
!> \param pw_pool the pool for the grids
!> \param compute_virial ...
!> \param virial_xc ...
!> \param exc_r the value of the xc functional in the real space
!> \par History
!>      JGH (13-Jun-2002): adaptation to new functionals
!>      Fawzi (11.2002): drho_g(1:3)->drho_g
!>      Fawzi (1.2003). lsd version
!>      Fawzi (11.2003): version using the new xc interface
!>      Fawzi (03.2004): fft free for smoothed density and derivs, gga lsd
!>      Fawzi (04.2004): metafunctionals
!>      mguidon (12.2008) : laplace functionals
!> \author fawzi; based LDA version of JGH, based on earlier version of apsi
!> \note
!>      Beware: some really dirty pointer handling!
!>      energy should be kept consistent with xc_exc_calc
! **************************************************************************************************
   SUBROUTINE xc_vxc_pw_create(vxc_rho, vxc_tau, exc, rho_r, rho_g, tau, xc_section, &
                               pw_pool, compute_virial, virial_xc, exc_r)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rho, vxc_tau
      REAL(KIND=dp), INTENT(out)                         :: exc
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, rho_g, tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      LOGICAL                                            :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: virial_xc
      TYPE(pw_p_type), INTENT(INOUT), OPTIONAL           :: exc_r

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xc_vxc_pw_create'
      CHARACTER(len=30), DIMENSION(2), PARAMETER :: &
         norm_drho_spin_name = (/"(norm_drhoa)", "(norm_drhob)"/)
      INTEGER, DIMENSION(3, 3), PARAMETER ::        nd = RESHAPE((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/))

      CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH)          :: desc
      INTEGER                                            :: handle, i, idir, ispin, j, jdir, k, &
                                                            n_deriv, npoints, nspins, order, &
                                                            xc_deriv_method_id, xc_rho_smooth_id
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: dealloc_pw_to_deriv, has_laplace, &
                                                            has_tau, lsd, use_virial, zero_result
      REAL(KIND=dp)                                      :: density_smooth_cut_range, drho_cutoff, &
                                                            ndr, rho_cutoff
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: deriv_data, norm_drho, norm_drho_spin, &
                                                            rho, rhoa, rhob, tmp_cr3d
      TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER        :: drho, drho_spin, drhoa, drhob
      TYPE(cp_sll_xc_deriv_type), POINTER                :: pos
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_p_type), DIMENSION(2)                      :: vxc_to_deriv
      TYPE(pw_p_type), DIMENSION(3)                      :: pw_to_deriv, pw_to_deriv_rho
      TYPE(pw_type), POINTER                             :: tmp_g, v_drho_r, virial_pw, vxc_g
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_derivative_type), POINTER                  :: deriv_att
      TYPE(xc_rho_set_type), POINTER                     :: rho_set

      CALL timeset(routineN, handle)
      NULLIFY (tmp_g, v_drho_r, vxc_g, norm_drho_spin, norm_drho, drho_spin, drhoa, &
               drhob, pos, rho_set, virial_pw)
      DO idir = 1, 3
         NULLIFY (pw_to_deriv(idir)%pw, pw_to_deriv_rho(idir)%pw)
      END DO
      DO i = 1, 2
         NULLIFY (vxc_to_deriv(i)%pw)
      END DO

      pw_grid => rho_r(1)%pw%pw_grid

      CPASSERT(ASSOCIATED(xc_section))
      CPASSERT(ASSOCIATED(pw_pool))
      CPASSERT(.NOT. ASSOCIATED(vxc_rho))
      CPASSERT(.NOT. ASSOCIATED(vxc_tau))
      nspins = SIZE(rho_r)
      CPASSERT(ASSOCIATED(rho_r(SIZE(rho_r))%pw))
      lsd = (nspins /= 1)
      IF (lsd) THEN
         CPASSERT(nspins == 2)
      END IF

      use_virial = compute_virial
      virial_xc = 0.0_dp

      bo = rho_r(1)%pw%pw_grid%bounds_local
      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      ! calculate the potential derivatives
      CALL xc_rho_set_and_dset_create(rho_set=rho_set, deriv_set=deriv_set, &
                                      deriv_order=1, rho_r=rho_r, rho_g=rho_g, tau=tau, &
                                      xc_section=xc_section, &
                                      pw_pool=pw_pool, &
                                      needs_basic_components=.TRUE.)

      CALL section_vals_val_get(xc_section, "XC_GRID%XC_DERIV", &
                                i_val=xc_deriv_method_id)
      CALL section_vals_val_get(xc_section, "XC_GRID%XC_SMOOTH_RHO", &
                                i_val=xc_rho_smooth_id)
      CALL section_vals_val_get(xc_section, "DENSITY_SMOOTH_CUTOFF_RANGE", &
                                r_val=density_smooth_cut_range)

      CALL xc_rho_set_get(rho_set, rho_cutoff=rho_cutoff, &
                          drho_cutoff=drho_cutoff)

      has_tau = .FALSE.
      has_laplace = .FALSE.
      ! check for unknown derivatives
      n_deriv = 0
      pos => deriv_set%derivs
      DO WHILE (cp_sll_xc_deriv_next(pos, el_att=deriv_att))
         CALL xc_derivative_get(deriv_att, order=order, &
                                desc=desc)
         IF (order == 1) THEN
            IF (lsd) THEN
               SELECT CASE (desc)
               CASE ("(rho)", "(rhoa)", "(rhob)", "(norm_drho)", "(norm_drhoa)", &
                     "(norm_drhob)")
                  n_deriv = n_deriv + 1
               CASE ("(tau)", "(tau_a)", "(tau_b)")
                  has_tau = .TRUE.
                  n_deriv = n_deriv + 1
               CASE ("(laplace_rhoa)", "(laplace_rhob)")
                  has_laplace = .TRUE.
                  n_deriv = n_deriv + 1
               CASE default
                  !FM if you are looking at this error probably you are missing the
                  !FM cross term (drhoa_drhob), I never got round to implement it,
                  !FM in the functionals that I have implemented I use norm_drho
                  !FM instead, either do the same or implement it, it shouldn't be
                  !FM difficult (I am just lazy ;)
                  CPABORT("unknown functional derivative (LSD): '"//TRIM(desc)//"'")
               END SELECT
            ELSE
               SELECT CASE (desc)
               CASE ("(rho)", "(norm_drho)")
                  n_deriv = n_deriv + 1
               CASE ("(tau)")
                  has_tau = .TRUE.
                  n_deriv = n_deriv + 1
               CASE ("(laplace_rho)")
                  has_laplace = .TRUE.
                  n_deriv = n_deriv + 1
               CASE default
                  CPABORT("unknown functional derivative (LDA): '"//TRIM(desc)//"'")
               END SELECT
            END IF
         END IF
      END DO

      ALLOCATE (vxc_rho(nspins))
      DO ispin = 1, nspins
         NULLIFY (vxc_rho(ispin)%pw)
      END DO

      CALL xc_rho_set_get(rho_set, rho=rho, rhoa=rhoa, rhob=rhob, &
                          can_return_null=.TRUE.)

      ! recover the vxc arrays
      IF (lsd) THEN
         deriv_att => xc_dset_get_derivative(deriv_set, "(rhoa)")
         IF (ASSOCIATED(deriv_att)) THEN
            CALL pw_create(vxc_rho(1)%pw, pw_grid=pw_grid, &
                           cr3d_ptr=deriv_att%deriv_data, &
                           use_data=REALDATA3D, in_space=REALSPACE)
            NULLIFY (deriv_att%deriv_data)
         ELSE
            CALL pw_pool_create_pw(pw_pool, vxc_rho(1)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(vxc_rho(1)%pw)
         END IF

         deriv_att => xc_dset_get_derivative(deriv_set, "(rhob)")
         IF (ASSOCIATED(deriv_att)) THEN
            CALL pw_create(vxc_rho(2)%pw, pw_grid=pw_grid, &
                           cr3d_ptr=deriv_att%deriv_data, &
                           use_data=REALDATA3D, in_space=REALSPACE)
            NULLIFY (deriv_att%deriv_data)
         ELSE
            CALL pw_pool_create_pw(pw_pool, vxc_rho(2)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(vxc_rho(2)%pw)
         END IF

      ELSE
         deriv_att => xc_dset_get_derivative(deriv_set, "(rho)")
         IF (ASSOCIATED(deriv_att)) THEN
            CALL pw_create(vxc_rho(1)%pw, pw_grid=pw_grid, &
                           cr3d_ptr=deriv_att%deriv_data, &
                           use_data=REALDATA3D, in_space=REALSPACE)
            NULLIFY (deriv_att%deriv_data)
         ELSE
            CALL pw_pool_create_pw(pw_pool, vxc_rho(1)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(vxc_rho(1)%pw)
         END IF
      END IF

      deriv_att => xc_dset_get_derivative(deriv_set, "(rho)")
      IF (ASSOCIATED(deriv_att)) THEN
         IF (lsd) THEN ! otherwise already taken care in vxc recovery
            CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
!$OMP       PARALLEL DEFAULT(NONE) SHARED(vxc_rho,deriv_data)
            DO ispin = 1, 2
!$OMP          WORKSHARE
               vxc_rho(ispin)%pw%cr3d(:, :, :) = vxc_rho(ispin)%pw%cr3d(:, :, :) + deriv_data(:, :, :)
!$OMP          END WORKSHARE NOWAIT
            END DO
!$OMP       END PARALLEL
            CALL pw_pool_give_back_cr3d(pw_pool, deriv_att%deriv_data)
            NULLIFY (deriv_att%deriv_data)
         END IF
      END IF

      ! rhoa,rhob already taken care of in vxc recovery

      deriv_att => xc_dset_get_derivative(deriv_set, "(norm_drho)")
      IF (ASSOCIATED(deriv_att)) THEN
         CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)

         CALL xc_rho_set_get(rho_set, norm_drho=norm_drho, &
                             drho=drho, drhoa=drhoa, drhob=drhob, rho_cutoff=rho_cutoff, &
                             drho_cutoff=drho_cutoff, &
                             can_return_null=.TRUE.)

         IF (ASSOCIATED(norm_drho)) THEN
!$OMP       PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_data,norm_drho,drho_cutoff)
            deriv_data(:, :, :) = -deriv_data(:, :, :)/MAX(norm_drho(:, :, :), drho_cutoff)
!$OMP       END PARALLEL WORKSHARE
         ELSE IF (ASSOCIATED(drho)) THEN
!$OMP         PARALLEL DO DEFAULT(NONE) SHARED(bo,deriv_data,drho,drho_cutoff)&
!$OMP                     PRIVATE(k,j,i,ndr) COLLAPSE(3)
            DO k = bo(1, 3), bo(2, 3)
               DO j = bo(1, 2), bo(2, 2)
                  DO i = bo(1, 1), bo(2, 1)
                     ndr = SQRT(drho(1)%array(i, j, k)**2 + &
                                drho(2)%array(i, j, k)**2 + &
                                drho(3)%array(i, j, k)**2)
                     deriv_data(i, j, k) = -deriv_data(i, j, k)/MAX(ndr, drho_cutoff)
                  END DO
               END DO
            END DO
         ELSE
            CPASSERT(ASSOCIATED(drhoa))
            CPASSERT(ASSOCIATED(drhob))
!$OMP         PARALLEL DO DEFAULT(NONE) SHARED(bo,deriv_data,drhoa,drhob,drho_cutoff)&
!$OMP                     PRIVATE(k,j,i,ndr) COLLAPSE(3)
            DO k = bo(1, 3), bo(2, 3)
               DO j = bo(1, 2), bo(2, 2)
                  DO i = bo(1, 1), bo(2, 1)
                     ndr = SQRT((drhoa(1)%array(i, j, k) + drhob(1)%array(i, j, k))**2 + &
                                (drhoa(2)%array(i, j, k) + drhob(2)%array(i, j, k))**2 + &
                                (drhoa(3)%array(i, j, k) + drhob(3)%array(i, j, k))**2)
                     deriv_data(i, j, k) = -deriv_data(i, j, k)/MAX(ndr, drho_cutoff)
                  END DO
               END DO
            END DO
         END IF

         IF (ASSOCIATED(drho) .AND. ASSOCIATED(drho(1)%array)) THEN
            CPASSERT(ASSOCIATED(deriv_data))
            IF (use_virial) THEN
               CALL pw_pool_create_pw(pw_pool, virial_pw, &
                                      use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_zero(virial_pw)
               DO idir = 1, 3
!$OMP             PARALLEL WORKSHARE DEFAULT(NONE) SHARED(virial_pw,drho,deriv_data,idir)
                  virial_pw%cr3d(:, :, :) = drho(idir)%array(:, :, :)*deriv_data(:, :, :)
!$OMP             END PARALLEL WORKSHARE
                  DO jdir = 1, idir
                     virial_xc(idir, jdir) = pw_grid%dvol* &
                                             accurate_dot_product(virial_pw%cr3d(:, :, :), &
                                                                  drho(jdir)%array(:, :, :))
                     virial_xc(jdir, idir) = virial_xc(idir, jdir)
                  END DO
               END DO
               CALL pw_pool_give_back_pw(pw_pool, virial_pw)
            END IF ! use_virial
            DO idir = 1, 3
               CALL pw_create(pw_to_deriv_rho(idir)%pw, pw_grid=pw_grid, &
                              cr3d_ptr=drho(idir)%array, &
                              use_data=REALDATA3D, in_space=REALSPACE)
               CPASSERT(ASSOCIATED(drho(idir)%array))
!$OMP          PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_data,drho,idir)
               drho(idir)%array(:, :, :) = drho(idir)%array(:, :, :)*deriv_data(:, :, :)
!$OMP          END PARALLEL WORKSHARE
               NULLIFY (drho(idir)%array)
            END DO
         ELSE
            IF (use_virial) THEN
               CALL pw_pool_create_pw(pw_pool, virial_pw, &
                                      use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_zero(virial_pw)
               DO idir = 1, 3
!$OMP             PARALLEL WORKSHARE DEFAULT(NONE) SHARED(virial_pw,drhoa,drhob,deriv_data,idir)
                  virial_pw%cr3d(:, :, :) = (drhoa(idir)%array(:, :, :) + drhob(idir)%array(:, :, :))*deriv_data(:, :, :)
!$OMP             END PARALLEL WORKSHARE
                  DO jdir = 1, idir
                     virial_xc(idir, jdir) = pw_grid%dvol*accurate_sum(virial_pw%cr3d(:, :, :)* &
                                                                       (drhoa(jdir)%array(:, :, :) + drhob(jdir)%array(:, :, :)))
                     virial_xc(jdir, idir) = virial_xc(idir, jdir)
                  END DO
               END DO
               CALL pw_pool_give_back_pw(pw_pool, virial_pw)
            END IF ! use_virial
            DO idir = 1, 3
               CALL pw_pool_create_pw(pw_pool, pw_to_deriv_rho(idir)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
!$OMP            PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_data,&
!$OMP                           pw_to_deriv_rho,drhoa,drhob,idir)
               pw_to_deriv_rho(idir)%pw%cr3d(:, :, :) = (drhoa(idir)%array(:, :, :) + drhob(idir)%array(:, :, :))* &
                                                        deriv_data(:, :, :)
!$OMP          END PARALLEL WORKSHARE
            END DO
         END IF

         CALL pw_pool_give_back_cr3d(pw_pool, deriv_att%deriv_data)

      END IF

      DO ispin = 1, nspins

         IF (lsd) THEN
            IF (ispin == 1) THEN
               CALL xc_rho_set_get(rho_set, drhoa=drho_spin, &
                                   can_return_null=.TRUE.)
               CALL xc_rho_set_get(rho_set, norm_drhoa=norm_drho_spin, &
                                   can_return_null=.TRUE.)
            ELSE
               CALL xc_rho_set_get(rho_set, drhob=drho_spin, &
                                   can_return_null=.TRUE.)
               CALL xc_rho_set_get(rho_set, norm_drhob=norm_drho_spin, &
                                   can_return_null=.TRUE.)
            END IF

            deriv_att => xc_dset_get_derivative(deriv_set, norm_drho_spin_name(ispin))
            IF (ASSOCIATED(deriv_att)) THEN
               CPASSERT(lsd)
               CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
               CALL pw_create(v_drho_r, pw_grid, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              cr3d_ptr=deriv_data)

               IF (ASSOCIATED(norm_drho_spin)) THEN
!$OMP             PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_data,norm_drho_spin,drho_cutoff)
                  deriv_data(:, :, :) = -deriv_data(:, :, :)/ &
                                        MAX(norm_drho_spin(:, :, :), drho_cutoff)
!$OMP             END PARALLEL WORKSHARE
               ELSE
!$OMP               PARALLEL DO DEFAULT(NONE) SHARED(bo,deriv_data,drho_spin,drho_cutoff)&
!$OMP                           PRIVATE(k,j,i, ndr) COLLAPSE(3)
                  DO k = bo(1, 3), bo(2, 3)
                     DO j = bo(1, 2), bo(2, 2)
                        DO i = bo(1, 1), bo(2, 1)
                           ndr = SQRT(drho_spin(1)%array(i, j, k)**2 + &
                                      drho_spin(2)%array(i, j, k)**2 + &
                                      drho_spin(3)%array(i, j, k)**2)
                           deriv_data(i, j, k) = -deriv_data(i, j, k)/MAX(ndr, drho_cutoff)
                        END DO
                     END DO
                  END DO
               END IF

               IF (use_virial) THEN
                  CALL pw_pool_create_pw(pw_pool, virial_pw, &
                                         use_data=REALDATA3D, &
                                         in_space=REALSPACE)
                  CALL pw_zero(virial_pw)
                  DO idir = 1, 3
!$OMP                PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_data,drho_spin,virial_pw,idir)
                     virial_pw%cr3d(:, :, :) = drho_spin(idir)%array(:, :, :)*deriv_data(:, :, :)
!$OMP                END PARALLEL WORKSHARE
                     DO jdir = 1, idir
                        virial_xc(idir, jdir) = virial_xc(idir, jdir) + pw_grid%dvol* &
                                                accurate_dot_product(virial_pw%cr3d(:, :, :), &
                                                                     drho_spin(jdir)%array(:, :, :))
                        virial_xc(jdir, idir) = virial_xc(idir, jdir)
                     END DO
                  END DO
                  CALL pw_pool_give_back_pw(pw_pool, virial_pw)
               END IF ! use_virial

               vxc_to_deriv(ispin)%pw => v_drho_r
               NULLIFY (v_drho_r, deriv_att%deriv_data)

               DO idir = 1, 3
                  CPASSERT(ASSOCIATED(drho_spin(idir)%array))
                  CPASSERT(ASSOCIATED(vxc_to_deriv(ispin)%pw%cr3d))
!$OMP             PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_data,drho_spin,&
!$OMP                           ispin,idir,vxc_to_deriv,drho_cutoff)
                  drho_spin(idir)%array(:, :, :) = vxc_to_deriv(ispin)%pw%cr3d(:, :, :)* &
                                                   drho_spin(idir)%array(:, :, :)
!$OMP             END PARALLEL WORKSHARE

                  CALL pw_create(pw_to_deriv(idir)%pw, pw_grid=pw_grid, &
                                 cr3d_ptr=drho_spin(idir)%array, &
                                 use_data=REALDATA3D, in_space=REALSPACE)
                  NULLIFY (drho_spin(idir)%array)
               END DO

               dealloc_pw_to_deriv = .TRUE.
               CALL pw_pool_give_back_pw(pw_pool, vxc_to_deriv(ispin)%pw)
            END IF ! deriv_att

         END IF ! LSD

         IF (ASSOCIATED(pw_to_deriv_rho(1)%pw)) THEN
            IF (.NOT. ASSOCIATED(pw_to_deriv(1)%pw)) THEN
               DO idir = 1, 3
                  pw_to_deriv(idir)%pw => pw_to_deriv_rho(idir)%pw
               END DO
               dealloc_pw_to_deriv = (.NOT. lsd) .OR. (ispin == 2)
               IF (dealloc_pw_to_deriv) THEN
                  DO idir = 1, 3
                     NULLIFY (pw_to_deriv_rho(idir)%pw)
                  END DO
               END IF
            ELSE
               DO idir = 1, 3
                  CALL pw_axpy(pw_to_deriv_rho(idir)%pw, pw_to_deriv(idir)%pw)
                  IF (ispin == 2) THEN
                     CALL pw_pool_give_back_pw(pw_pool, pw_to_deriv_rho(idir)%pw)
                  END IF
               END DO
            END IF
         END IF

         IF (ASSOCIATED(pw_to_deriv(1)%pw)) THEN
            ! partial integration
            IF (xc_deriv_method_id /= xc_deriv_pw) THEN
               CALL pw_spline_scale_deriv(pw_to_deriv, &
                                          transpose=.TRUE.)
            END IF

            IF (xc_deriv_method_id == xc_deriv_pw .OR. &
                xc_deriv_method_id == xc_deriv_spline3 .OR. &
                xc_deriv_method_id == xc_deriv_spline2) THEN

               IF (.NOT. ASSOCIATED(vxc_g)) THEN
                  CALL pw_pool_create_pw(pw_pool, vxc_g, &
                                         use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                  zero_result = .TRUE.
               ELSE
                  zero_result = .FALSE.
               END IF

               DO idir = 1, 3
                  IF (zero_result .AND. idir == 1) THEN
                     tmp_g => vxc_g
                  ELSE
                     CALL pw_pool_create_pw(pw_pool, tmp_g, &
                                            use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                  END IF

                  CALL pw_transfer(pw_to_deriv(idir)%pw, tmp_g)

                  SELECT CASE (xc_deriv_method_id)
                  CASE (xc_deriv_pw)
                     CALL pw_derive(tmp_g, nd(:, idir))
                  CASE (xc_deriv_spline2)
                     CALL pw_spline2_interpolate_values_g(tmp_g)
                     CALL pw_spline2_deriv_g(tmp_g, idir=idir)
                  CASE (xc_deriv_spline3)
                     CALL pw_spline3_interpolate_values_g(tmp_g)
                     CALL pw_spline3_deriv_g(tmp_g, idir=idir)
                  CASE default
                     CPABORT("")
                  END SELECT

                  IF (zero_result .AND. idir == 1) THEN
                     NULLIFY (tmp_g)
                  ELSE
                     CALL pw_axpy(tmp_g, vxc_g)
                     CALL pw_pool_give_back_pw(pw_pool, tmp_g)
                  END IF
                  IF (dealloc_pw_to_deriv) THEN
                     CALL pw_pool_give_back_pw(pw_pool, pw_to_deriv(idir)%pw)
                  END IF
               END DO
               ! transfer vxc in real space
               CALL pw_pool_create_pw(pw_pool, v_drho_r, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_transfer(vxc_g, v_drho_r)
               CALL pw_axpy(v_drho_r, vxc_rho(ispin)%pw)
               CALL pw_pool_give_back_pw(pw_pool, v_drho_r)
               CALL pw_pool_give_back_pw(pw_pool, vxc_g)
            ELSE
               v_drho_r => vxc_rho(ispin)%pw
               DO idir = 1, 3
                  SELECT CASE (xc_deriv_method_id)
                  CASE (xc_deriv_spline2_smooth)
                     CALL pw_nn_deriv_r(pw_in=pw_to_deriv(idir)%pw, &
                                        pw_out=v_drho_r, coeffs=spline2_deriv_coeffs, &
                                        idir=idir)
                  CASE (xc_deriv_spline3_smooth)
                     CALL pw_nn_deriv_r(pw_in=pw_to_deriv(idir)%pw, &
                                        pw_out=v_drho_r, coeffs=spline3_deriv_coeffs, &
                                        idir=idir)
                  CASE (xc_deriv_nn10_smooth)
                     CALL pw_nn_deriv_r(pw_in=pw_to_deriv(idir)%pw, &
                                        pw_out=v_drho_r, coeffs=nn10_deriv_coeffs, &
                                        idir=idir)
                  CASE (xc_deriv_nn50_smooth)
                     CALL pw_nn_deriv_r(pw_in=pw_to_deriv(idir)%pw, &
                                        pw_out=v_drho_r, coeffs=nn50_deriv_coeffs, &
                                        idir=idir)
                  CASE default
                     CPABORT("")
                  END SELECT
                  IF (dealloc_pw_to_deriv) THEN
                     CALL pw_pool_give_back_pw(pw_pool, pw_to_deriv(idir)%pw)
                  END IF
               END DO
               NULLIFY (v_drho_r)
            END IF
         END IF

         ! ** Add laplace part to vxc_rho
         IF (has_laplace) THEN
            IF (lsd) THEN
               IF (ispin == 1) THEN
                  deriv_att => xc_dset_get_derivative(deriv_set, "(laplace_rhoa)")
                  CPASSERT(ASSOCIATED(deriv_att))
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
               ELSE
                  deriv_att => xc_dset_get_derivative(deriv_set, "(laplace_rhob)")
                  CPASSERT(ASSOCIATED(deriv_att))
                  CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
               END IF

            ELSE
               deriv_att => xc_dset_get_derivative(deriv_set, "(laplace_rho)")
               CPASSERT(ASSOCIATED(deriv_att))
               CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
            END IF

            IF (use_virial) CALL virial_laplace(rho_r(ispin)%pw, pw_pool, virial_xc, deriv_data)

            SELECT CASE (xc_deriv_method_id)
            CASE (xc_deriv_pw)
               CALL pw_pool_create_pw(pw_pool, pw_to_deriv(1)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_pool_create_pw(pw_pool, tmp_g, &
                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
               CALL pw_zero(tmp_g)
               CALL pw_pool_create_pw(pw_pool, v_drho_r, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_zero(v_drho_r)
!$OMP          PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_data,&
!$OMP                         pw_to_deriv,idir)
               pw_to_deriv(1)%pw%cr3d(:, :, :) = deriv_data(:, :, :)
!$OMP          END PARALLEL WORKSHARE
               CALL pw_transfer(pw_to_deriv(1)%pw, tmp_g)

               CALL pw_laplace(tmp_g)

               ! Add this to the potential
               CALL pw_transfer(tmp_g, v_drho_r)

               CALL pw_axpy(v_drho_r, vxc_rho(ispin)%pw)
               CALL pw_pool_give_back_pw(pw_pool, v_drho_r)
               CALL pw_pool_give_back_pw(pw_pool, pw_to_deriv(1)%pw)
               CALL pw_pool_give_back_pw(pw_pool, tmp_g)
            CASE default
               CPABORT("")
            END SELECT
         END IF

         IF (pw_grid%spherical) THEN
            ! filter vxc
            CALL pw_pool_create_pw(pw_pool, vxc_g, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            CALL pw_transfer(vxc_rho(ispin)%pw, vxc_g)
            CALL pw_transfer(vxc_g, vxc_rho(ispin)%pw)
            CALL pw_pool_give_back_pw(pw_pool, vxc_g)
         END IF
         CALL smooth_cutoff(pot=vxc_rho(ispin)%pw%cr3d, rho=rho, rhoa=rhoa, rhob=rhob, &
                            rho_cutoff=rho_cutoff*density_smooth_cut_range, &
                            rho_smooth_cutoff_range=density_smooth_cut_range)

         ! final smoothing if rho was smoothed
         IF (xc_rho_smooth_id /= xc_rho_no_smooth) THEN
            CALL pw_pool_create_pw(pw_pool, v_drho_r, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(v_drho_r)
            SELECT CASE (xc_rho_smooth_id)
            CASE (xc_rho_spline2_smooth)
               CALL pw_nn_smear_r(pw_in=vxc_rho(ispin)%pw, pw_out=v_drho_r, &
                                  coeffs=spline2_coeffs)
            CASE (xc_rho_spline3_smooth)
               CALL pw_nn_smear_r(pw_in=vxc_rho(ispin)%pw, pw_out=v_drho_r, &
                                  coeffs=spline3_coeffs)
            CASE (xc_rho_nn10)
               CALL pw_nn_smear_r(pw_in=vxc_rho(ispin)%pw, pw_out=v_drho_r, &
                                  coeffs=nn10_coeffs)
            CASE (xc_rho_nn50)
               CALL pw_nn_smear_r(pw_in=vxc_rho(ispin)%pw, pw_out=v_drho_r, &
                                  coeffs=nn50_coeffs)
            CASE default
               CPABORT("")
            END SELECT
            deriv_data => vxc_rho(ispin)%pw%cr3d
            vxc_rho(ispin)%pw%cr3d => v_drho_r%cr3d
            v_drho_r%cr3d => deriv_data
            CALL pw_pool_give_back_pw(pw_pool, v_drho_r)
         END IF
      END DO

      DO idir = 1, 3
         CPASSERT(.NOT. ASSOCIATED(pw_to_deriv(idir)%pw))
         CPASSERT(.NOT. ASSOCIATED(pw_to_deriv_rho(idir)%pw))
      END DO

      ! 0-deriv -> value of exc
      ! this has to be kept consistent with xc_exc_calc
      IF (n_deriv > 0) THEN
         deriv_att => xc_dset_get_derivative(deriv_set, "")
         CPASSERT(ASSOCIATED(deriv_att))
         CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)

         CALL pw_create(v_drho_r, pw_grid, &
                        use_data=REALDATA3D, in_space=REALSPACE, &
                        cr3d_ptr=deriv_data)
         NULLIFY (v_drho_r%cr3d)
         CALL pw_release(v_drho_r)

         CALL smooth_cutoff(pot=deriv_data, rho=rho, rhoa=rhoa, rhob=rhob, &
                            rho_cutoff=rho_cutoff, &
                            rho_smooth_cutoff_range=density_smooth_cut_range)

         exc = accurate_sum(deriv_data)*pw_grid%dvol
         IF (pw_grid%para%mode == PW_MODE_DISTRIBUTED) THEN
            CALL mp_sum(exc, pw_grid%para%group)
         END IF
         !
         ! return the xc functional value at the grid points
         !
         IF (PRESENT(exc_r)) THEN
            CPASSERT(.NOT. ASSOCIATED(exc_r%pw))
            deriv_att => xc_dset_get_derivative(deriv_set, "")
            CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
            CALL pw_create(exc_r%pw, pw_grid, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           cr3d_ptr=deriv_data)
            NULLIFY (deriv_att%deriv_data)
         END IF
      ELSE
         exc = 0.0_dp
      END IF

      CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)

      ! tau part
      IF (has_tau) THEN
         ALLOCATE (vxc_tau(nspins))
         DO ispin = 1, nspins
            NULLIFY (vxc_tau(ispin)%pw)
         END DO
         IF (lsd) THEN
            deriv_att => xc_dset_get_derivative(deriv_set, "(tau_a)")
            IF (ASSOCIATED(deriv_att)) THEN
               CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)

               CALL pw_create(vxc_tau(1)%pw, pw_grid, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              cr3d_ptr=deriv_data)
               NULLIFY (deriv_att%deriv_data)
            END IF
            deriv_att => xc_dset_get_derivative(deriv_set, "(tau_b)")
            IF (ASSOCIATED(deriv_att)) THEN
               CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)

               CALL pw_create(vxc_tau(2)%pw, pw_grid, &
                              use_data=REALDATA3D, in_space=REALSPACE, &
                              cr3d_ptr=deriv_data)
               NULLIFY (deriv_att%deriv_data)
            END IF
            deriv_att => xc_dset_get_derivative(deriv_set, "(tau)")
            IF (ASSOCIATED(deriv_att)) THEN
               CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
               IF (ASSOCIATED(vxc_tau(1)%pw)) THEN
                  DO ispin = 1, 2
                     CPASSERT(ASSOCIATED(vxc_tau(ispin)%pw))
                     tmp_cr3d => vxc_tau(ispin)%pw%cr3d
                     CALL daxpy(npoints, 1.0_dp, deriv_data, 1, tmp_cr3d, 1)
                  END DO
               ELSE
                  CALL pw_create(vxc_tau(1)%pw, pw_grid, &
                                 use_data=REALDATA3D, in_space=REALSPACE, &
                                 cr3d_ptr=deriv_data)
                  NULLIFY (deriv_att%deriv_data)
                  CALL pw_pool_create_pw(pw_pool, vxc_tau(2)%pw, &
                                         use_data=REALDATA3D, in_space=REALSPACE)
                  CALL pw_copy(vxc_tau(1)%pw, vxc_tau(2)%pw)
               END IF
            END IF
         ELSE
            deriv_att => xc_dset_get_derivative(deriv_set, "(tau)")
            CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
            CALL pw_create(vxc_tau(1)%pw, pw_grid, &
                           use_data=REALDATA3D, in_space=REALSPACE, &
                           cr3d_ptr=deriv_data)
            NULLIFY (deriv_att%deriv_data)
         END IF
         DO ispin = 1, nspins
            CPASSERT(ASSOCIATED(vxc_tau(ispin)%pw))
         END DO
      END IF
      CALL xc_dset_release(deriv_set)

      CALL timestop(handle)

   END SUBROUTINE xc_vxc_pw_create

! **************************************************************************************************
!> \brief calculates just the exchange and correlation energy
!>      (no vxc)
!> \param rho_r      realspace density on the grid
!> \param rho_g      g-space density on the grid
!> \param tau        kinetic energy density on the grid
!> \param xc_section XC parameters
!> \param pw_pool    pool of plain-wave grids
!> \return the XC energy
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
!> \note
!>      has to be kept consistent with xc_vxc_pw_create
! **************************************************************************************************
   FUNCTION xc_exc_calc(rho_r, rho_g, tau, xc_section, pw_pool) &
      RESULT(exc)
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, rho_g, tau
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(kind=dp)                                      :: exc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xc_exc_calc'

      INTEGER                                            :: handle
      REAL(dp)                                           :: density_smooth_cut_range, rho_cutoff
      REAL(dp), DIMENSION(:, :, :), POINTER              :: e_0
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_derivative_type), POINTER                  :: deriv
      TYPE(xc_rho_set_type), POINTER                     :: rho_set

      CALL timeset(routineN, handle)
      NULLIFY (rho_set, deriv, e_0)
      exc = 0.0_dp

      ! this has to be consistent with what is done in xc_vxc_pw_create
      CALL xc_rho_set_and_dset_create(rho_set=rho_set, &
                                      deriv_set=deriv_set, deriv_order=0, &
                                      rho_r=rho_r, rho_g=rho_g, tau=tau, xc_section=xc_section, &
                                      pw_pool=pw_pool, &
                                      needs_basic_components=.FALSE.)
      deriv => xc_dset_get_derivative(deriv_set, "")

      IF (ASSOCIATED(deriv)) THEN
         CALL xc_derivative_get(deriv, deriv_data=e_0)

         CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", &
                                   r_val=rho_cutoff)
         CALL section_vals_val_get(xc_section, "DENSITY_SMOOTH_CUTOFF_RANGE", &
                                   r_val=density_smooth_cut_range)
         CALL smooth_cutoff(pot=e_0, rho=rho_set%rho, &
                            rhoa=rho_set%rhoa, rhob=rho_set%rhob, &
                            rho_cutoff=rho_cutoff, &
                            rho_smooth_cutoff_range=density_smooth_cut_range)

         exc = accurate_sum(e_0)*rho_r(1)%pw%pw_grid%dvol
         IF (rho_r(1)%pw%pw_grid%para%mode == PW_MODE_DISTRIBUTED) THEN
            CALL mp_sum(exc, rho_r(1)%pw%pw_grid%para%group)
         END IF

         CALL xc_rho_set_release(rho_set, pw_pool=pw_pool)
         CALL xc_dset_release(deriv_set)
      END IF
      CALL timestop(handle)
   END FUNCTION xc_exc_calc

! **************************************************************************************************
!> \brief Caller routine to calculate the second order potential in the direction of rho1_r
!> \param v_xc XC potential, will be allocated, to be integrated with the KS density
!> \param v_xc_tau ...
!> \param deriv_set XC derivatives from xc_prep_2nd_deriv
!> \param rho_set XC rho set from KS rho from xc_prep_2nd_deriv
!> \param rho1_r first-order density in r space
!> \param rho1_g first-order density in g space
!> \param tau1_r ...
!> \param pw_pool pw pool to create new grids
!> \param xc_section XC section to calculate the derivatives from
!> \param gapw whether to carry out GAPW (not possible with numerical derivatives)
!> \param vxg GAPW potential
!> \param lsd_singlets ...
!> \param do_excitations ...
!> \param do_triplet ...
!> \param do_tddft ...
!> \param compute_virial ...
!> \param virial_xc virial terms will be collected here
! **************************************************************************************************
   SUBROUTINE xc_calc_2nd_deriv(v_xc, v_xc_tau, deriv_set, rho_set, rho1_r, rho1_g, tau1_r, &
                                pw_pool, xc_section, gapw, vxg, &
                                lsd_singlets, do_excitations, do_triplet, do_tddft, &
                                compute_virial, virial_xc)

      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc, v_xc_tau
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_r, rho1_g, tau1_r
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: pw_pool
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL, INTENT(IN)                                :: gapw
      REAL(KIND=dp), DIMENSION(:, :, :, :), OPTIONAL, &
         POINTER                                         :: vxg
      LOGICAL, INTENT(IN), OPTIONAL                      :: lsd_singlets, do_excitations, &
                                                            do_triplet, do_tddft, compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc

      CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_2nd_deriv', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd, my_compute_virial, &
                                                            my_do_excitations, my_do_tddft, &
                                                            my_do_triplet, my_lsd_singlets
      REAL(KIND=dp)                                      :: fac
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set

      CALL timeset(routineN, handle)

      my_compute_virial = .FALSE.
      IF (PRESENT(compute_virial)) my_compute_virial = compute_virial

      my_do_tddft = .FALSE.
      IF (PRESENT(do_tddft)) my_do_tddft = do_tddft

      my_do_excitations = .FALSE.
      IF (PRESENT(do_excitations)) my_do_excitations = do_excitations

      my_lsd_singlets = .FALSE.
      IF (PRESENT(lsd_singlets)) my_lsd_singlets = lsd_singlets

      my_do_triplet = .FALSE.
      IF (PRESENT(do_triplet)) my_do_triplet = do_triplet

      nspins = SIZE(rho1_r)
      lsd = (nspins == 2)
      IF (nspins == 1 .AND. my_do_excitations .AND. (my_lsd_singlets .OR. my_do_triplet)) THEN
         nspins = 2
         lsd = .TRUE.
      END IF

      NULLIFY (v_xc, v_xc_tau)
      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      IF (needs%tau .OR. needs%tau_spin) THEN
         IF (.NOT. ASSOCIATED(tau1_r)) &
            CPABORT("Tau-dependent functionals requires allocated kinetic energy density grid")
         ALLOCATE (v_xc_tau(nspins))
         DO ispin = 1, nspins
            NULLIFY (v_xc_tau(ispin)%pw)
            CALL pw_pool_create_pw(pw_pool, v_xc_tau(ispin)%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(v_xc_tau(ispin)%pw)
         END DO
      END IF

      IF (section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL") .AND. .NOT. my_do_tddft) THEN
         !------!
         ! rho1 !
         !------!
         bo = rho1_r(1)%pw%pw_grid%bounds_local
         NULLIFY (rho1_set)
         ! create the place where to store the argument for the functionals
         CALL xc_rho_set_create(rho1_set, bo, &
                                rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                                tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

         ! calculate the arguments needed by the functionals
         CALL xc_rho_set_update(rho1_set, rho1_r, rho1_g, tau1_r, needs, &
                                section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                pw_pool)

         fac = 0._dp
         IF (nspins == 1 .AND. my_do_excitations) THEN
            IF (my_lsd_singlets) fac = 1.0_dp
            IF (my_do_triplet) fac = -1.0_dp
         END IF

         CALL xc_calc_2nd_deriv_analytical(v_xc, v_xc_tau, deriv_set, rho_set, &
                                           rho1_set, pw_pool, xc_section, gapw, vxg=vxg, &
                                           tddfpt_fac=fac, compute_virial=compute_virial, virial_xc=virial_xc)

         CALL xc_rho_set_release(rho1_set)

      ELSE
         IF (gapw) CPABORT("Numerical 2nd derivatives not implemented with GAPW")

         CALL xc_calc_2nd_deriv_numerical(v_xc, v_xc_tau, rho_set, rho1_r, rho1_g, tau1_r, &
                                          pw_pool, xc_section, &
                                          my_do_excitations .AND. my_do_triplet, &
                                          compute_virial, virial_xc, deriv_set)
      END IF

      CALL timestop(handle)

   END SUBROUTINE xc_calc_2nd_deriv

! **************************************************************************************************
!> \brief calculates 2nd derivative numerically
!> \param v_xc potential to be calculated (has to be allocated already)
!> \param v_tau tau-part of the potential to be calculated (has to be allocated already)
!> \param rho_set KS density from xc_prep_2nd_deriv
!> \param rho1_r first-order density in r-space
!> \param rho1_g first-order density in g-space
!> \param tau1_r first-order kinetic-energy density in r-space
!> \param pw_pool pw pool for new grids
!> \param xc_section XC section to calculate the derivatives from
!> \param do_triplet ...
!> \param calc_virial whether to calculate virial terms
!> \param virial_xc collects stress tensor components (no metaGGAs!)
!> \param deriv_set deriv set from xc_prep_2nd_deriv (only for virials)
! **************************************************************************************************
   SUBROUTINE xc_calc_2nd_deriv_numerical(v_xc, v_tau, rho_set, rho1_r, rho1_g, tau1_r, &
                                          pw_pool, xc_section, &
                                          do_triplet, calc_virial, virial_xc, deriv_set)

      TYPE(pw_p_type), DIMENSION(:), INTENT(IN), POINTER :: v_xc, v_tau
      TYPE(xc_rho_set_type), INTENT(IN), POINTER         :: rho_set
      TYPE(pw_p_type), DIMENSION(:), INTENT(IN), POINTER :: rho1_r, rho1_g, tau1_r
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: pw_pool
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL, INTENT(IN)                                :: do_triplet
      LOGICAL, INTENT(IN), OPTIONAL                      :: calc_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc
      TYPE(xc_derivative_set_type), OPTIONAL             :: deriv_set

      CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_2nd_deriv_numerical', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, idir, ispin, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: gradient_f, lsd, my_calc_virial, tau_f, laplace_f
      REAL(KIND=dp)                                      :: exc, gradient_cut, h
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dr1dr, dra1dra, drb1drb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_dummy
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: norm_drho, norm_drho2, norm_drho2a, &
                                                            norm_drho2b, norm_drhoa, norm_drhob, &
                                                            rho, rho1, rho1a, rho1b, rhoa, rhob, &
                                                            tau_a, tau_b, tau, tau1, tau1a, tau1b, laplace, laplace1, &
                                                            laplacea, laplaceb, laplace1a, laplace1b, laplace2, laplace2a, laplace2b
      TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER        :: drho, drho1, drho1a, drho1b, drhoa, drhob
      TYPE(pw_p_type)                                    :: v_drho, v_drhoa, v_drhob
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, vxc_rho_1, &
                                                            vxc_rho_2, vxc_rho_3, vxc_rho_4, tau_r, &
                                                            vxc_tau_1, vxc_tau_2, vxc_tau_3, vxc_tau_4
      TYPE(pw_p_type), DIMENSION(:), ALLOCATABLE :: v_laplace
      TYPE(pw_type), POINTER                             :: virial_pw
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_derivative_set_type)                       :: deriv_set1, deriv_set2
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set, rho2_set

      CALL timeset(routineN, handle)

      my_calc_virial = .FALSE.
      IF (PRESENT(calc_virial) .AND. PRESENT(virial_xc)) my_calc_virial = calc_virial

      nspins = SIZE(v_xc)

      NULLIFY (tau, tau_r, tau_a, tau_b)

      h = section_get_rval(xc_section, "STEP_SIZE")

      IF (nspins == 2) THEN
         NULLIFY (vxc_rho_1, vxc_rho_2, rho_g, vxc_tau_1, vxc_tau_2)
         ALLOCATE (rho_r(2))
         DO ispin = 1, nspins
            NULLIFY (rho_r(ispin)%pw)
            CALL pw_pool_create_pw(pw_pool, rho_r(ispin)%pw, in_space=REALSPACE, use_data=REALDATA3D)
         END DO
         IF (ASSOCIATED(tau1_r) .AND. ASSOCIATED(v_tau)) THEN
            ALLOCATE (tau_r(2))
            DO ispin = 1, nspins
               NULLIFY (tau_r(ispin)%pw)
               CALL pw_pool_create_pw(pw_pool, tau_r(ispin)%pw, in_space=REALSPACE, use_data=REALDATA3D)
            END DO
         END IF
         CALL xc_rho_set_get(rho_set, can_return_null=.TRUE., rhoa=rhoa, rhob=rhob, tau_a=tau_a, tau_b=tau_b)
         CALL calc_resp_potential_numer_ab(rho_r, rho_g, rho1_r, rhoa, rhob, vxc_rho_1, &
                                           tau_r, tau1_r, tau_a, tau_b, vxc_tau_1, xc_section, pw_pool, h)
         CALL calc_resp_potential_numer_ab(rho_r, rho_g, rho1_r, rhoa, rhob, vxc_rho_2, &
                                           tau_r, tau1_r, tau_a, tau_b, vxc_tau_2, xc_section, pw_pool, -h)
         DO ispin = 1, nspins
!$OMP PARALLEL DEFAULT(NONE) SHARED(v_xc,vxc_rho_1,vxc_rho_2,h,v_tau,vxc_tau_1,vxc_tau_2,tau1_r,ispin)
!$OMP WORKSHARE
            v_xc(ispin)%pw%cr3d(:, :, :) = (vxc_rho_1(ispin)%pw%cr3d(:, :, :) - vxc_rho_2(ispin)%pw%cr3d(:, :, :))/h
!$OMP END WORKSHARE NOWAIT
            IF (ASSOCIATED(tau1_r) .AND. ASSOCIATED(v_tau) .AND. ASSOCIATED(vxc_tau_1) .AND. ASSOCIATED(vxc_tau_2)) THEN
!$OMP WORKSHARE
               v_tau(ispin)%pw%cr3d(:, :, :) = (vxc_tau_1(ispin)%pw%cr3d(:, :, :) - vxc_tau_2(ispin)%pw%cr3d(:, :, :))/h
!$OMP END WORKSHARE NOWAIT
            END IF
!$OMP END PARALLEL
         END DO
         DO ispin = 1, nspins
            CALL pw_release(vxc_rho_1(ispin)%pw)
            CALL pw_release(vxc_rho_2(ispin)%pw)
         END DO
         DEALLOCATE (vxc_rho_1, vxc_rho_2)
         IF (ASSOCIATED(vxc_tau_1)) THEN
            DO ispin = 1, nspins
               CALL pw_release(vxc_tau_1(ispin)%pw)
               CALL pw_release(vxc_tau_2(ispin)%pw)
            END DO
            DEALLOCATE (vxc_tau_1, vxc_tau_2)
         END IF
      ELSE IF (nspins == 1 .AND. do_triplet) THEN
         NULLIFY (vxc_rho_1, vxc_rho_2, vxc_rho_3, vxc_rho_4, vxc_tau_1, vxc_tau_2, vxc_tau_3, vxc_tau_4, rho_g)
         ALLOCATE (rho_r(2))
         DO ispin = 1, 2
            NULLIFY (rho_r(ispin)%pw)
            CALL pw_pool_create_pw(pw_pool, rho_r(ispin)%pw, in_space=REALSPACE, use_data=REALDATA3D)
         END DO
         IF (ASSOCIATED(tau1_r) .AND. ASSOCIATED(v_tau)) THEN
            ALLOCATE (tau_r(2))
            DO ispin = 1, nspins
               NULLIFY (tau_r(ispin)%pw)
               CALL pw_pool_create_pw(pw_pool, tau_r(ispin)%pw, in_space=REALSPACE, use_data=REALDATA3D)
            END DO
         END IF
         CALL xc_rho_set_get(rho_set, can_return_null=.TRUE., rhoa=rhoa, rhob=rhob, tau_a=tau_a, tau_b=tau_b)
         ! K(alpha,alpha)
!$OMP PARALLEL DEFAULT(NONE) SHARED(rho_r,rhoa,rhob,h,rho1_r,tau_r,tau_a,tau_b,tau1_r)
!$OMP WORKSHARE
         rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :) + 0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
         rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :)
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(tau1_r)) THEN
!$OMP WORKSHARE
            tau_r(1)%pw%cr3d(:, :, :) = tau_a(:, :, :) + 0.5_dp*h*tau1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
            tau_r(2)%pw%cr3d(:, :, :) = tau_b(:, :, :)
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP END PARALLEL
         CALL xc_vxc_pw_create(vxc_rho_1, vxc_tau_2, exc, rho_r, rho_g, tau_r, xc_section, &
                               pw_pool, .FALSE., virial_dummy)
!$OMP PARALLEL DEFAULT(NONE) SHARED(rho_r,rhoa,rhob,h,rho1_r,tau_r,tau_a,tau_b,tau1_r)
!$OMP WORKSHARE
         rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :) - 0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
         rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :)
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(tau1_r)) THEN
!$OMP WORKSHARE
            tau_r(1)%pw%cr3d(:, :, :) = tau_a(:, :, :) - 0.5_dp*h*tau1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
            tau_r(2)%pw%cr3d(:, :, :) = tau_b(:, :, :)
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP END PARALLEL
         CALL xc_vxc_pw_create(vxc_rho_2, vxc_tau_2, exc, rho_r, rho_g, tau_r, xc_section, &
                               pw_pool, .FALSE., virial_dummy)
!$OMP PARALLEL DEFAULT(NONE) SHARED(rho_r,rhoa,rhob,h,rho1_r,v_xc,vxc_rho_1,vxc_rho_2,&
!$OMP                               v_tau,vxc_tau_1,vxc_tau_2,tau_r,tau_a,tau_b,tau1_r)
!$OMP WORKSHARE
         v_xc(1)%pw%cr3d(:, :, :) = (vxc_rho_1(1)%pw%cr3d(:, :, :) - vxc_rho_2(1)%pw%cr3d(:, :, :))/h
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(vxc_tau_1) .AND. ASSOCIATED(vxc_tau_2) .AND. ASSOCIATED(v_tau)) THEN
!$OMP WORKSHARE
            v_tau(1)%pw%cr3d(:, :, :) = (vxc_tau_1(1)%pw%cr3d(:, :, :) - vxc_tau_2(1)%pw%cr3d(:, :, :))/h
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP WORKSHARE
         ! K(alpha,beta)
         rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
         rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :) + 0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(tau1_r)) THEN
!$OMP WORKSHARE
            tau_r(1)%pw%cr3d(:, :, :) = tau_a(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
            tau_r(2)%pw%cr3d(:, :, :) = tau_b(:, :, :) + 0.5_dp*h*tau1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP END PARALLEL
         CALL xc_vxc_pw_create(vxc_rho_3, vxc_tau_3, exc, rho_r, rho_g, tau_r, xc_section, &
                               pw_pool, .FALSE., virial_dummy)
!$OMP PARALLEL DEFAULT(NONE) SHARED(rho_r,rhoa,rhob,h,rho1_r,tau_r,tau_a,tau_b,tau1_r)
!$OMP WORKSHARE
         rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
         rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :) - 0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(tau1_r) .AND. ASSOCIATED(tau_a) .AND. ASSOCIATED(tau_b) .AND. ASSOCIATED(tau_r)) THEN
!$OMP WORKSHARE
            tau_r(1)%pw%cr3d(:, :, :) = tau_a(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
            tau_r(2)%pw%cr3d(:, :, :) = tau_b(:, :, :) - 0.5_dp*h*tau1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP END PARALLEL
         CALL xc_vxc_pw_create(vxc_rho_4, vxc_tau_4, exc, rho_r, rho_g, tau_r, xc_section, &
                               pw_pool, .FALSE., virial_dummy)
!$OMP PARALLEL DEFAULT(NONE) SHARED(v_xc,vxc_rho_3,vxc_rho_4,rhob,h,v_tau,vxc_tau_3,vxc_tau_4)
!$OMP WORKSHARE
         v_xc(1)%pw%cr3d(:, :, :) = v_xc(1)%pw%cr3d(:, :, :) - &
                                    (vxc_rho_3(1)%pw%cr3d(:, :, :) - vxc_rho_4(1)%pw%cr3d(:, :, :))/h
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(vxc_tau_3) .AND. ASSOCIATED(vxc_tau_4) .AND. ASSOCIATED(v_tau)) THEN
!$OMP WORKSHARE
            v_tau(1)%pw%cr3d(:, :, :) = v_tau(1)%pw%cr3d(:, :, :) - &
                                        (vxc_tau_3(1)%pw%cr3d(:, :, :) - vxc_tau_4(1)%pw%cr3d(:, :, :))/h
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP END PARALLEL
         DO ispin = 1, 2
            CALL pw_release(vxc_rho_1(ispin)%pw)
            CALL pw_release(vxc_rho_2(ispin)%pw)
            CALL pw_release(vxc_rho_3(ispin)%pw)
            CALL pw_release(vxc_rho_4(ispin)%pw)
         END DO
         DEALLOCATE (vxc_rho_1, vxc_rho_2, vxc_rho_3, vxc_rho_4)
         IF (ASSOCIATED(vxc_tau_1)) THEN
         DO ispin = 1, 2
            CALL pw_release(vxc_tau_1(ispin)%pw)
            CALL pw_release(vxc_tau_2(ispin)%pw)
            CALL pw_release(vxc_tau_3(ispin)%pw)
            CALL pw_release(vxc_tau_4(ispin)%pw)
         END DO
         DEALLOCATE (vxc_tau_1, vxc_tau_2, vxc_tau_3, vxc_tau_4)
         END IF
      ELSE
         NULLIFY (vxc_rho_1, vxc_rho_2, rho_r, rho_g, vxc_tau_1, vxc_tau_2, tau_r, tau)
         ALLOCATE (rho_r(1))
         NULLIFY (rho_r(1)%pw)
         CALL pw_pool_create_pw(pw_pool, rho_r(1)%pw, in_space=REALSPACE, use_data=REALDATA3D)
         IF (ASSOCIATED(tau1_r) .AND. ASSOCIATED(v_tau)) THEN
            ALLOCATE (tau_r(1))
            NULLIFY (tau_r(1)%pw)
            CALL pw_pool_create_pw(pw_pool, tau_r(1)%pw, in_space=REALSPACE, use_data=REALDATA3D)
         END IF
         CALL xc_rho_set_get(rho_set, can_return_null=.TRUE., rho=rho, tau=tau)
!$OMP PARALLEL DEFAULT(NONE) SHARED(rho_r,rho,h,rho1_r,tau1_r,tau,tau_r)
!$OMP WORKSHARE
         rho_r(1)%pw%cr3d(:, :, :) = rho(:, :, :) + 0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(tau1_r) .AND. ASSOCIATED(tau) .AND. ASSOCIATED(tau1_r)) THEN
!$OMP WORKSHARE
            tau_r(1)%pw%cr3d(:, :, :) = tau(:, :, :) + 0.5_dp*h*tau1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP END PARALLEL
         CALL xc_vxc_pw_create(vxc_rho_1, vxc_tau_1, exc, rho_r, rho_g, tau_r, xc_section, &
                               pw_pool, .FALSE., virial_dummy)
!$OMP PARALLEL DEFAULT(NONE) SHARED(rho_r,rho,h,rho1_r,tau_r,tau,tau1_r)
!$OMP WORKSHARE
         rho_r(1)%pw%cr3d(:, :, :) = rho(:, :, :) - 0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(tau1_r) .AND. ASSOCIATED(tau_r) .AND. ASSOCIATED(tau)) THEN
!$OMP WORKSHARE
            tau_r(1)%pw%cr3d(:, :, :) = tau(:, :, :) - 0.5_dp*h*tau1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP END PARALLEL
         CALL xc_vxc_pw_create(vxc_rho_2, vxc_tau_2, exc, rho_r, rho_g, tau_r, xc_section, &
                               pw_pool, .FALSE., virial_dummy)
!$OMP PARALLEL DEFAULT(NONE) SHARED(v_xc,vxc_rho_1,vxc_rho_2,h,v_tau,vxc_tau_1,vxc_tau_2)
!$OMP WORKSHARE
         v_xc(1)%pw%cr3d(:, :, :) = (vxc_rho_1(1)%pw%cr3d(:, :, :) - vxc_rho_2(1)%pw%cr3d(:, :, :))/h
!$OMP END WORKSHARE NOWAIT
         IF (ASSOCIATED(vxc_tau_1) .AND. ASSOCIATED(vxc_tau_2) .AND. ASSOCIATED(v_tau)) THEN
!$OMP WORKSHARE
            v_tau(1)%pw%cr3d(:, :, :) = (vxc_tau_1(1)%pw%cr3d(:, :, :) - vxc_tau_2(1)%pw%cr3d(:, :, :))/h
!$OMP END WORKSHARE NOWAIT
         END IF
!$OMP END PARALLEL
         CALL pw_release(vxc_rho_1(1)%pw)
         CALL pw_release(vxc_rho_2(1)%pw)
         DEALLOCATE (vxc_rho_1, vxc_rho_2)
         IF (ASSOCIATED(vxc_tau_1)) THEN
            CALL pw_release(vxc_tau_1(1)%pw)
            CALL pw_release(vxc_tau_2(1)%pw)
            DEALLOCATE (vxc_tau_1, vxc_tau_2)
         END IF
      END IF

      IF (my_calc_virial) THEN
         lsd = (nspins == 2)
         IF (nspins == 1 .AND. do_triplet) THEN
            lsd = .TRUE.
         END IF

         CALL check_for_gradients(deriv_set, (nspins == 2), gradient_f, tau_f, laplace_f)

         IF (laplace_f) THEN
            CPABORT("Stress tensors with Laplace-dependent functionals not supported!")
         END IF

         ! Calculate the virial terms
         ! Those arising from the first derivatives are treated like in xc_calc_2nd_deriv_analytical
         ! Those arising from the second derivatives are calculated numerically
         ! We assume that all metaGGA functionals require the gradient
         IF (gradient_f) THEN
            bo = rho_set%local_bounds

            ! Create the work grid for the virial terms
            CALL allocate_pw(virial_pw, pw_pool, bo)

            gradient_cut = section_get_rval(xc_section, "GRADIENT_CUTOFF")

            NULLIFY (rho1_set)
            ! create the container to store the argument of the functionals
            CALL xc_rho_set_create(rho1_set, bo, &
                                   rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                   drho_cutoff=gradient_cut, &
                                   tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

            xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
            needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

            ! calculate the arguments needed by the functionals
            CALL xc_rho_set_update(rho1_set, rho1_r, rho1_g, tau1_r, needs, &
                                   section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                   section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                   pw_pool)

            IF (lsd) THEN
               CALL xc_rho_set_get(rho_set, drhoa=drhoa, drhob=drhob, norm_drho=norm_drho, &
                                   norm_drhoa=norm_drhoa, norm_drhob=norm_drhob, tau_a=tau_a, tau_b=tau_b, &
                                   laplace_rhoa=laplacea, laplace_rhob=laplaceb, can_return_null=.TRUE.)
               CALL xc_rho_set_get(rho1_set, drhoa=drho1a, drhob=drho1b, laplace_rhoa=laplace1a, &
                                   laplace_rhob=laplace1b, can_return_null=.TRUE.)

               NULLIFY (drho, drho1)
               CALL calc_drho_from_ab(drho, drhoa, drhob)
               CALL calc_drho_from_ab(drho1, drho1a, drho1b)
            ELSE
               CALL xc_rho_set_get(rho_set, drho=drho, norm_drho=norm_drho, tau=tau, laplace_rho=laplace, can_return_null=.TRUE.)
               CALL xc_rho_set_get(rho1_set, drho=drho1, laplace_rho=laplace1, can_return_null=.TRUE.)
            END IF

            IF (laplace_f) THEN
               ALLOCATE (v_laplace(nspins))
               DO ispin = 1, nspins
                  NULLIFY (v_laplace(ispin)%pw)
                  CALL allocate_pw(v_laplace(ispin)%pw, pw_pool, bo)
               END DO

            END IF

            CALL prepare_dr1dr(dr1dr, drho, drho1)
            CALL divide_by_norm_drho(deriv_set, rho_set, lsd)

            IF (lsd) THEN
               CALL prepare_dr1dr(dra1dra, drhoa, drho1a)
               CALL prepare_dr1dr(drb1drb, drhob, drho1b)

               CALL allocate_pw(v_drho%pw, pw_pool, bo)
               CALL allocate_pw(v_drhoa%pw, pw_pool, bo)
               CALL allocate_pw(v_drhob%pw, pw_pool, bo)

               IF (ASSOCIATED(norm_drhoa)) CALL apply_drho(deriv_set, "(norm_drhoa)", virial_pw, drhoa, drho1a, virial_xc, &
                                                           norm_drhoa, gradient_cut, dra1dra, v_drhoa%pw%cr3d)
               IF (ASSOCIATED(norm_drhob)) CALL apply_drho(deriv_set, "(norm_drhob)", virial_pw, drhob, drho1b, virial_xc, &
                                                           norm_drhob, gradient_cut, drb1drb, v_drhob%pw%cr3d)
               IF (ASSOCIATED(norm_drho)) CALL apply_drho(deriv_set, "(norm_drho)", virial_pw, drho, drho1, virial_xc, &
                                                          norm_drho, gradient_cut, dr1dr, v_drho%pw%cr3d)

            ELSE

               ! Create the work grid for the potential of the gradient part
               CALL allocate_pw(v_drho%pw, pw_pool, bo)

               CALL apply_drho(deriv_set, "(norm_drho)", virial_pw, drho, drho1, virial_xc, &
                               norm_drho, gradient_cut, dr1dr, v_drho%pw%cr3d)
            END IF

            IF (lsd) THEN
               rho_r(1)%pw%cr3d = rhoa
               rho_r(2)%pw%cr3d = rhob
            ELSE
               rho_r(1)%pw%cr3d = rho
            END IF
            IF (ASSOCIATED(tau1_r)) THEN
            IF (lsd) THEN
               tau_r(1)%pw%cr3d = tau_a
               tau_r(2)%pw%cr3d = tau_b
            ELSE
               tau_r(1)%pw%cr3d = tau
            END IF
            END IF

            ! Create deriv sets with same densities but different gradients
            CALL xc_dset_create(deriv_set1, pw_pool)
            CALL xc_dset_create(deriv_set2, pw_pool)

            NULLIFY (rho2_set)
            ! create the place where to store the argument for the functionals
            CALL xc_rho_set_create(rho2_set, bo, &
                                   rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                   drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                                   tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

            ! calculate the arguments needed by the functionals
            CALL xc_rho_set_update(rho2_set, rho_r, rho_g, tau_r, needs, &
                                   section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                   section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                   pw_pool)

            IF (lsd) THEN
               CALL xc_rho_set_get(rho1_set, rhoa=rho1a, rhob=rho1b, tau_a=tau1a, tau_b=tau1b, &
                                   laplace_rhoa=laplace1a, laplace_rhob=laplace1b, can_return_null=.TRUE.)
               CALL xc_rho_set_get(rho2_set, norm_drhoa=norm_drho2a, norm_drhob=norm_drho2b, &
                                   norm_drho=norm_drho2, laplace_rhoa=laplace2a, laplace_rhob=laplace2b, can_return_null=.TRUE.)

               IF (ASSOCIATED(norm_drhoa)) THEN
                  CALL get_derivs_rho(norm_drho2a, norm_drhoa, h, xc_fun_section, lsd, rho2_set, deriv_set1, deriv_set2)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(rhoa)", bo, &
                                        norm_drhoa, gradient_cut, h, rho1a, v_drhoa%pw%cr3d)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(rhob)", bo, &
                                        norm_drhoa, gradient_cut, h, rho1b, v_drhoa%pw%cr3d)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(norm_drhoa)", bo, &
                                        norm_drhoa, gradient_cut, h, dra1dra, v_drhoa%pw%cr3d)
                  CALL update_deriv_drho_ab(deriv_set1, deriv_set2, "(norm_drhob)", bo, &
                                            norm_drhoa, gradient_cut, h, dra1dra, drb1drb, v_drhoa%pw%cr3d, v_drhob%pw%cr3d)
                  CALL update_deriv_drho_ab(deriv_set1, deriv_set2, "(norm_drho)", bo, &
                                            norm_drhoa, gradient_cut, h, dra1dra, dr1dr, v_drhoa%pw%cr3d, v_drho%pw%cr3d)
                  IF (ASSOCIATED(tau1a) .AND. ASSOCIATED(tau1b)) THEN
                     CALL update_deriv_rho(deriv_set1, deriv_set2, "(tau_a)", bo, &
                                           norm_drhoa, gradient_cut, h, tau1a, v_drhoa%pw%cr3d)
                     CALL update_deriv_rho(deriv_set1, deriv_set2, "(tau_b)", bo, &
                                           norm_drhoa, gradient_cut, h, tau1b, v_drhoa%pw%cr3d)
                  END IF
               END IF

               IF (ASSOCIATED(norm_drhob)) THEN
                  CALL get_derivs_rho(norm_drho2b, norm_drhob, h, xc_fun_section, lsd, rho2_set, deriv_set1, deriv_set2)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(rhoa)", bo, &
                                        norm_drhob, gradient_cut, h, rho1a, v_drhob%pw%cr3d)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(rhob)", bo, &
                                        norm_drhob, gradient_cut, h, rho1b, v_drhob%pw%cr3d)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(norm_drhob)", bo, &
                                        norm_drhob, gradient_cut, h, drb1drb, v_drhob%pw%cr3d)
                  CALL update_deriv_drho_ab(deriv_set1, deriv_set2, "(norm_drhoa)", bo, &
                                            norm_drhob, gradient_cut, h, drb1drb, dra1dra, v_drhob%pw%cr3d, v_drhoa%pw%cr3d)
                  CALL update_deriv_drho_ab(deriv_set1, deriv_set2, "(norm_drho)", bo, &
                                            norm_drhob, gradient_cut, h, drb1drb, dr1dr, v_drhob%pw%cr3d, v_drho%pw%cr3d)
                  IF (ASSOCIATED(tau1a) .AND. ASSOCIATED(tau1b)) THEN
                     CALL update_deriv_rho(deriv_set1, deriv_set2, "(tau_a)", bo, &
                                           norm_drhob, gradient_cut, h, tau1a, v_drhob%pw%cr3d)
                     CALL update_deriv_rho(deriv_set1, deriv_set2, "(tau_b)", bo, &
                                           norm_drhob, gradient_cut, h, tau1b, v_drhob%pw%cr3d)
                  END IF
               END IF

               IF (ASSOCIATED(norm_drho)) THEN
                  CALL get_derivs_rho(norm_drho2, norm_drho, h, xc_fun_section, lsd, rho2_set, deriv_set1, deriv_set2)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(rhoa)", bo, &
                                        norm_drho, gradient_cut, h, rho1a, v_drho%pw%cr3d)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(rhob)", bo, &
                                        norm_drho, gradient_cut, h, rho1b, v_drho%pw%cr3d)
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(norm_drho)", bo, &
                                        norm_drho, gradient_cut, h, dr1dr, v_drho%pw%cr3d)
                  CALL update_deriv_drho_ab(deriv_set1, deriv_set2, "(norm_drhoa)", bo, &
                                            norm_drho, gradient_cut, h, dr1dr, dra1dra, v_drho%pw%cr3d, v_drhoa%pw%cr3d)
                  CALL update_deriv_drho_ab(deriv_set1, deriv_set2, "(norm_drhob)", bo, &
                                            norm_drho, gradient_cut, h, dr1dr, drb1drb, v_drho%pw%cr3d, v_drhob%pw%cr3d)
                  IF (ASSOCIATED(tau1a) .AND. ASSOCIATED(tau1b)) THEN
                     CALL update_deriv_rho(deriv_set1, deriv_set2, "(tau_a)", bo, &
                                           norm_drho, gradient_cut, h, tau1a, v_drho%pw%cr3d)
                     CALL update_deriv_rho(deriv_set1, deriv_set2, "(tau_b)", bo, &
                                           norm_drho, gradient_cut, h, tau1b, v_drho%pw%cr3d)
                  END IF
               END IF

               CALL virial_drho_drho(virial_pw, drhoa, v_drhoa, virial_xc)
               CALL virial_drho_drho(virial_pw, drhob, v_drhob, virial_xc)
               CALL virial_drho_drho(virial_pw, drho, v_drho, virial_xc)

               CALL deallocate_pw(virial_pw, pw_pool)
               CALL deallocate_pw(v_drho%pw, pw_pool)
               CALL deallocate_pw(v_drhoa%pw, pw_pool)
               CALL deallocate_pw(v_drhob%pw, pw_pool)

               DO idir = 1, 3
                  DEALLOCATE (drho(idir)%array)
                  DEALLOCATE (drho1(idir)%array)
               END DO
               DEALLOCATE (drho, drho1)
               DEALLOCATE (dra1dra, drb1drb)

            ELSE
               CALL xc_rho_set_get(rho1_set, rho=rho1, tau=tau1, laplace_rho=laplace1, can_return_null=.TRUE.)
               CALL xc_rho_set_get(rho2_set, norm_drho=norm_drho2, laplace_rho=laplace2, can_return_null=.TRUE.)

               CALL get_derivs_rho(norm_drho2, norm_drho, h, xc_fun_section, lsd, rho2_set, deriv_set1, deriv_set2)

               ! Obtain the numerical 2nd derivatives w.r.t. to drho and collect the potential
               CALL update_deriv_rho(deriv_set1, deriv_set2, "(rho)", bo, &
                                     norm_drho, gradient_cut, h, rho1, v_drho%pw%cr3d)
               CALL update_deriv_rho(deriv_set1, deriv_set2, "(norm_drho)", bo, &
                                     norm_drho, gradient_cut, h, dr1dr, v_drho%pw%cr3d)

               IF (ASSOCIATED(tau1)) THEN
                  CALL update_deriv_rho(deriv_set1, deriv_set2, "(tau)", bo, &
                                        norm_drho, gradient_cut, h, tau1, v_drho%pw%cr3d)
               END IF

               ! Calculate the virial contribution from the potential
               CALL virial_drho_drho(virial_pw, drho, v_drho, virial_xc)

               CALL deallocate_pw(virial_pw, pw_pool)
               CALL deallocate_pw(v_drho%pw, pw_pool)

            END IF

            CALL xc_dset_release(deriv_set1)
            CALL xc_dset_release(deriv_set2)

            DEALLOCATE (dr1dr)

            CALL xc_rho_set_release(rho1_set)
            CALL xc_rho_set_release(rho2_set)
         END IF

      END IF

      DO ispin = 1, SIZE(rho_r)
         CALL pw_pool_give_back_pw(pw_pool, rho_r(ispin)%pw)
      END DO
      DEALLOCATE (rho_r)

      IF (ASSOCIATED(tau_r)) THEN
      DO ispin = 1, SIZE(tau_r)
         CALL pw_pool_give_back_pw(pw_pool, tau_r(ispin)%pw)
      END DO
      DEALLOCATE (tau_r)
      END IF

      CALL timestop(handle)

   END SUBROUTINE xc_calc_2nd_deriv_numerical

! **************************************************************************************************
!> \brief ...
!> \param rho_r ...
!> \param rho_g ...
!> \param rho1_r ...
!> \param rhoa ...
!> \param rhob ...
!> \param vxc_rho ...
!> \param tau_r ...
!> \param tau1_r ...
!> \param tau_a ...
!> \param tau_b ...
!> \param vxc_tau ...
!> \param xc_section ...
!> \param pw_pool ...
!> \param h ...
! **************************************************************************************************
   SUBROUTINE calc_resp_potential_numer_ab(rho_r, rho_g, rho1_r, rhoa, rhob, vxc_rho, &
                                           tau_r, tau1_r, tau_a, tau_b, vxc_tau, &
                                           xc_section, pw_pool, h)

      TYPE(pw_p_type), DIMENSION(:), POINTER, INTENT(IN) :: vxc_rho, vxc_tau
      TYPE(pw_p_type), DIMENSION(:), INTENT(IN), POINTER :: rho1_r, tau1_r
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: pw_pool
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_resp_potential_numer_ab', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: exc, h
      REAL(KIND=dp), DIMENSION(3, 3)                     :: virial_dummy
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rhoa, rhob, tau_a, tau_b
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r, tau_r

      CALL timeset(routineN, handle)

!$OMP PARALLEL DEFAULT(NONE) SHARED(rho_r,rhoa,rhob,h,rho1_r,tau_r,tau_a,tau_b,tau1_r)
!$OMP WORKSHARE
      rho_r(1)%pw%cr3d(:, :, :) = rhoa(:, :, :) + 0.5_dp*h*rho1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
      rho_r(2)%pw%cr3d(:, :, :) = rhob(:, :, :) + 0.5_dp*h*rho1_r(2)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
      IF (ASSOCIATED(tau1_r) .AND. ASSOCIATED(tau_r) .AND. ASSOCIATED(tau_a) .AND. ASSOCIATED(tau_b)) THEN
!$OMP WORKSHARE
         tau_r(1)%pw%cr3d(:, :, :) = tau_a(:, :, :) + 0.5_dp*h*tau1_r(1)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
         tau_r(2)%pw%cr3d(:, :, :) = tau_b(:, :, :) + 0.5_dp*h*tau1_r(2)%pw%cr3d(:, :, :)
!$OMP END WORKSHARE NOWAIT
      END IF
!$OMP END PARALLEL
      CALL xc_vxc_pw_create(vxc_rho, vxc_tau, exc, rho_r, rho_g, tau_r, xc_section, &
                            pw_pool, .FALSE., virial_dummy)

      CALL timestop(handle)

   END SUBROUTINE calc_resp_potential_numer_ab

! **************************************************************************************************
!> \brief calculates stress tensor and potential contributions from the first derivative
!> \param deriv_set ...
!> \param description ...
!> \param virial_pw ...
!> \param drho ...
!> \param drho1 ...
!> \param virial_xc ...
!> \param norm_drho ...
!> \param gradient_cut ...
!> \param dr1dr ...
!> \param v_drho ...
! **************************************************************************************************
   SUBROUTINE apply_drho(deriv_set, description, virial_pw, drho, drho1, virial_xc, norm_drho, gradient_cut, dr1dr, v_drho)

      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      CHARACTER(len=*), INTENT(in)                       :: description
      TYPE(pw_type), INTENT(IN)                          :: virial_pw
      TYPE(cp_3d_r_p_type), DIMENSION(3), INTENT(IN)     :: drho, drho1
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: norm_drho
      REAL(KIND=dp), INTENT(IN)                          :: gradient_cut
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: dr1dr
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v_drho

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_drho', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: deriv_data
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      CALL timeset(routineN, handle)

      deriv_att => xc_dset_get_derivative(deriv_set, TRIM(description))
      IF (ASSOCIATED(deriv_att)) THEN
         CALL xc_derivative_get(deriv_att, deriv_data=deriv_data)
         CALL virial_drho_drho1(virial_pw, drho, drho1, deriv_data, virial_xc)

!$OMP    PARALLEL WORKSHARE DEFAULT(NONE) SHARED(dr1dr,gradient_cut,norm_drho,v_drho,deriv_data)
         v_drho(:, :, :) = v_drho(:, :, :) + &
                           deriv_data(:, :, :)*dr1dr(:, :, :)/MAX(gradient_cut, norm_drho(:, :, :))**2
!$OMP    END PARALLEL WORKSHARE
      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_drho

! **************************************************************************************************
!> \brief adds potential contributions from derivatives of rho or diagonal terms of norm_drho
!> \param deriv_set1 ...
!> \param deriv_set2 ...
!> \param description ...
!> \param bo ...
!> \param norm_drho norm_drho of which derivative is calculated
!> \param gradient_cut ...
!> \param h ...
!> \param rho1 function to contract the derivative with (rho1 for rho, dr1dr for norm_drho)
!> \param v_drho ...
! **************************************************************************************************
   SUBROUTINE update_deriv_rho(deriv_set1, deriv_set2, description, bo, norm_drho, gradient_cut, h, rho1, v_drho)

      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set1, deriv_set2
      CHARACTER(len=*), INTENT(in)                       :: description
      INTEGER, DIMENSION(2, 3), INTENT(IN)               :: bo
      REAL(KIND=dp), DIMENSION(bo(1, 1):bo(2, 1), bo(1, &
                                                     2):bo(2, 2), bo(1, 3):bo(2, 3)), INTENT(IN)     :: norm_drho
      REAL(KIND=dp), INTENT(IN)                          :: gradient_cut, h
      REAL(KIND=dp), DIMENSION(bo(1, 1):bo(2, 1), bo(1, &
                                                     2):bo(2, 2), bo(1, 3):bo(2, 3)), INTENT(IN)     :: rho1
      REAL(KIND=dp), DIMENSION(bo(1, 1):bo(2, 1), bo(1, &
                                                     2):bo(2, 2), bo(1, 3):bo(2, 3)), INTENT(INOUT)  :: v_drho

      CHARACTER(len=*), PARAMETER :: routineN = 'update_deriv_rho', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k
      REAL(KIND=dp)                                      :: de
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: deriv_data1, deriv_data2
      TYPE(xc_derivative_type), POINTER                  :: deriv_att1, deriv_att2

      CALL timeset(routineN, handle)

      ! Obtain the numerical 2nd derivatives w.r.t. to drho and collect the potential
      deriv_att1 => xc_dset_get_derivative(deriv_set1, TRIM(description))
      deriv_att2 => xc_dset_get_derivative(deriv_set2, TRIM(description))
      IF (ASSOCIATED(deriv_att1) .AND. ASSOCIATED(deriv_att2)) THEN
         CALL xc_derivative_get(deriv_att1, deriv_data=deriv_data1)
         CALL xc_derivative_get(deriv_att2, deriv_data=deriv_data2)
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(bo,deriv_data1,deriv_data2,h,norm_drho,v_drho,rho1,gradient_cut) &
!$OMP PRIVATE(i,j,k,de) COLLAPSE(3)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  de = (deriv_data1(i, j, k) - deriv_data2(i, j, k))/(h*MAX(gradient_cut, norm_drho(i, j, k))**2)
                  v_drho(i, j, k) = v_drho(i, j, k) - de*rho1(i, j, k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE update_deriv_rho

! **************************************************************************************************
!> \brief adds potential contributions from derivatives of a component
!> \param deriv_set1 ...
!> \param deriv_set2 ...
!> \param description ...
!> \param bo ...
!> \param h ...
!> \param rho1 function to contract the derivative with (rho1 for rho, dr1dr for norm_drho)
!> \param v ...
! **************************************************************************************************
   SUBROUTINE update_deriv(deriv_set1, deriv_set2, description, bo, rho, h, rho1, v)

      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set1, deriv_set2
      CHARACTER(len=*), INTENT(in)                       :: description
      INTEGER, DIMENSION(2, 3), INTENT(IN)               :: bo
      REAL(KIND=dp), INTENT(IN)                          :: h
      REAL(KIND=dp), DIMENSION(bo(1, 1):bo(2, 1), bo(1, &
                                                     2):bo(2, 2), bo(1, 3):bo(2, 3)), INTENT(IN)     :: rho1, rho
      REAL(KIND=dp), DIMENSION(bo(1, 1):bo(2, 1), bo(1, &
                                                     2):bo(2, 2), bo(1, 3):bo(2, 3)), INTENT(INOUT)  :: v

      CHARACTER(len=*), PARAMETER :: routineN = 'update_deriv', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k
      REAL(KIND=dp)                                      :: de
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: deriv_data1, deriv_data2
      TYPE(xc_derivative_type), POINTER                  :: deriv_att1, deriv_att2

      CALL timeset(routineN, handle)

      ! Obtain the numerical 2nd derivatives w.r.t. to drho and collect the potential
      deriv_att1 => xc_dset_get_derivative(deriv_set1, TRIM(description))
      deriv_att2 => xc_dset_get_derivative(deriv_set2, TRIM(description))
      IF (ASSOCIATED(deriv_att1) .AND. ASSOCIATED(deriv_att2)) THEN
         CALL xc_derivative_get(deriv_att1, deriv_data=deriv_data1)
         CALL xc_derivative_get(deriv_att2, deriv_data=deriv_data2)
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(bo,deriv_data1,deriv_data2,h,v,rho1,rho) &
!$OMP PRIVATE(i,j,k,de) COLLAPSE(3)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  de = (deriv_data1(i, j, k) - deriv_data2(i, j, k))/h
                  v(i, j, k) = v(i, j, k) - de*rho1(i, j, k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE update_deriv

! **************************************************************************************************
!> \brief adds mixed derivatives of norm_drho
!> \param deriv_set1 ...
!> \param deriv_set2 ...
!> \param description ...
!> \param bo ...
!> \param norm_drhoa norm_drho of which derivatives is calculated
!> \param gradient_cut ...
!> \param h ...
!> \param dra1dra dr1dr corresponding to norm_drho
!> \param drb1drb ...
!> \param v_drhoa potential corresponding to norm_drho
!> \param v_drhob ...
! **************************************************************************************************
   SUBROUTINE update_deriv_drho_ab(deriv_set1, deriv_set2, description, bo, &
                                   norm_drhoa, gradient_cut, h, dra1dra, drb1drb, v_drhoa, v_drhob)

      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set1, deriv_set2
      CHARACTER(len=*), INTENT(in)                       :: description
      INTEGER, DIMENSION(2, 3), INTENT(IN)               :: bo
      REAL(KIND=dp), DIMENSION(bo(1, 1):bo(2, 1), bo(1, &
                                                     2):bo(2, 2), bo(1, 3):bo(2, 3)), INTENT(IN)     :: norm_drhoa
      REAL(KIND=dp), INTENT(IN)                          :: gradient_cut, h
      REAL(KIND=dp), DIMENSION(bo(1, 1):bo(2, 1), bo(1, &
                                                     2):bo(2, 2), bo(1, 3):bo(2, 3)), INTENT(IN)     :: dra1dra, drb1drb
      REAL(KIND=dp), DIMENSION(bo(1, 1):bo(2, 1), bo(1, &
                                                     2):bo(2, 2), bo(1, 3):bo(2, 3)), INTENT(INOUT)  :: v_drhoa, v_drhob

      CHARACTER(len=*), PARAMETER :: routineN = 'update_deriv_drho_ab', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, k
      REAL(KIND=dp)                                      :: de
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: deriv_data1, deriv_data2
      TYPE(xc_derivative_type), POINTER                  :: deriv_att1, deriv_att2

      CALL timeset(routineN, handle)

      deriv_att1 => xc_dset_get_derivative(deriv_set1, TRIM(description))
      deriv_att2 => xc_dset_get_derivative(deriv_set2, TRIM(description))
      IF (ASSOCIATED(deriv_att1) .AND. ASSOCIATED(deriv_att2)) THEN
         CALL xc_derivative_get(deriv_att1, deriv_data=deriv_data1)
         CALL xc_derivative_get(deriv_att1, deriv_data=deriv_data2)
!$OMP    PARALLEL DO PRIVATE(k,j,i,de)  DEFAULT(NONE) &
!$OMP    SHARED(bo,drb1drb,dra1dra,deriv_data1,deriv_data2,h,gradient_cut,norm_drhoa,v_drhoa,v_drhob) COLLAPSE(3)
         DO k = bo(1, 3), bo(2, 3)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 1), bo(2, 1)
                  ! We introduce a factor of two because we will average between both numerical derivatives
                  de = 0.5_dp*(deriv_data1(i, j, k) - deriv_data2(i, j, k))/(h*MAX(gradient_cut, norm_drhoa(i, j, k))**2)
                  v_drhoa(i, j, k) = v_drhoa(i, j, k) - de*drb1drb(i, j, k)
                  v_drhob(i, j, k) = v_drhob(i, j, k) - de*dra1dra(i, j, k)
               END DO
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE update_deriv_drho_ab

! **************************************************************************************************
!> \brief calculate derivative sets for helper points
!> \param norm_drho2 norm_drho of new points
!> \param norm_drho norm_drho of KS density
!> \param h ...
!> \param xc_fun_section ...
!> \param lsd ...
!> \param rho2_set rho_set for new points
!> \param deriv_set1 will contain derivatives for norm_drho*(1+0.5*h)
!> \param deriv_set2 will contain derivatives for norm_drho*(1-0.5*h)
! **************************************************************************************************
   SUBROUTINE get_derivs_rho(norm_drho2, norm_drho, h, xc_fun_section, lsd, rho2_set, deriv_set1, deriv_set2)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: norm_drho2
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: norm_drho
      REAL(KIND=dp), INTENT(IN)                          :: h
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_fun_section
      LOGICAL, INTENT(IN)                                :: lsd
      TYPE(xc_rho_set_type), INTENT(IN), POINTER         :: rho2_set
      TYPE(xc_derivative_set_type)                       :: deriv_set1, deriv_set2

      CHARACTER(len=*), PARAMETER :: routineN = 'get_derivs_rho', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Copy the densities, do one step into the direction of drho
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(norm_drho,norm_drho2,h)
      norm_drho2 = norm_drho*(1.0_dp + 0.5_dp*h)
!$OMP END PARALLEL WORKSHARE

      CALL xc_dset_zero_all(deriv_set1)

      ! Calculate the derivatives of the functional
      CALL xc_functionals_eval(xc_fun_section, &
                               lsd=lsd, &
                               rho_set=rho2_set, &
                               deriv_set=deriv_set1, &
                               deriv_order=1)

      ! Copy the densities, do one step into the direction of drho
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(norm_drho,norm_drho2,h)
      norm_drho2 = norm_drho*(1.0_dp - 0.5_dp*h)
!$OMP END PARALLEL WORKSHARE

      CALL xc_dset_zero_all(deriv_set2)

      ! Calculate the derivatives of the functional
      CALL xc_functionals_eval(xc_fun_section, &
                               lsd=lsd, &
                               rho_set=rho2_set, &
                               deriv_set=deriv_set2, &
                               deriv_order=1)

      ! Return to the original values
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(norm_drho,norm_drho2)
      norm_drho2 = norm_drho
!$OMP END PARALLEL WORKSHARE

      CALL divide_by_norm_drho(deriv_set1, rho2_set, lsd)
      CALL divide_by_norm_drho(deriv_set2, rho2_set, lsd)

      CALL timestop(handle)

   END SUBROUTINE get_derivs_rho

! **************************************************************************************************
!> \brief calculate derivative sets for helper points
!> \param norm_drho2 norm_drho of new points
!> \param norm_drho norm_drho of KS density
!> \param h ...
!> \param xc_fun_section ...
!> \param lsd ...
!> \param rho2_set rho_set for new points
!> \param deriv_set1 will contain derivatives for norm_drho*(1+0.5*h)
!> \param deriv_set2 will contain derivatives for norm_drho*(1-0.5*h)
! **************************************************************************************************
   SUBROUTINE get_derivs_laplace(laplace2, laplace, h, xc_fun_section, lsd, rho2_set, deriv_set1, deriv_set2)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: laplace2
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: laplace
      REAL(KIND=dp), INTENT(IN)                          :: h
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_fun_section
      LOGICAL, INTENT(IN)                                :: lsd
      TYPE(xc_rho_set_type), INTENT(IN), POINTER         :: rho2_set
      TYPE(xc_derivative_set_type)                       :: deriv_set1, deriv_set2

      CHARACTER(len=*), PARAMETER :: routineN = 'get_derivs_laplace', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Copy the densities, do one step into the direction of drho
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(laplace,laplace2,h)
      laplace2 = laplace + h*0.5_dp
!$OMP END PARALLEL WORKSHARE

      CALL xc_dset_zero_all(deriv_set1)

      ! Calculate the derivatives of the functional
      CALL xc_functionals_eval(xc_fun_section, &
                               lsd=lsd, &
                               rho_set=rho2_set, &
                               deriv_set=deriv_set1, &
                               deriv_order=1)

      ! Copy the densities, do one step into the direction of drho
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(laplace,laplace2,h)
      laplace2 = laplace - h*0.5_dp
!$OMP END PARALLEL WORKSHARE

      CALL xc_dset_zero_all(deriv_set2)

      ! Calculate the derivatives of the functional
      CALL xc_functionals_eval(xc_fun_section, &
                               lsd=lsd, &
                               rho_set=rho2_set, &
                               deriv_set=deriv_set2, &
                               deriv_order=1)

      ! Return to the original values
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(laplace,laplace2)
      laplace2 = laplace
!$OMP END PARALLEL WORKSHARE

      CALL divide_by_norm_drho(deriv_set1, rho2_set, lsd)
      CALL divide_by_norm_drho(deriv_set2, rho2_set, lsd)

      CALL timestop(handle)

   END SUBROUTINE get_derivs_laplace

! **************************************************************************************************
!> \brief Calculates the second derivative of E_xc at rho in the direction
!>      rho1  (if you see the second derivative as bilinear form)
!>      partial_rho|_(rho=rho) partial_rho|_(rho=rho) E_xc drho(rho1)drho
!>      The other direction is still undetermined, thus it returns
!>      a potential (partial integration is performed to reduce it to
!>      function of rho, removing the dependence from its partial derivs)
!>      Has to be called after the setup by xc_prep_2nd_deriv.
!> \param v_xc       exchange-correlation potential
!> \param v_xc_tau ...
!> \param deriv_set  derivatives of the exchange-correlation potential
!> \param rho_set    object containing the density at which the derivatives were calculated
!> \param rho1_set   object containing the density with which to fold
!> \param pw_pool    the pool for the grids
!> \param xc_section XC parameters
!> \param gapw       Gaussian and augmented plane waves calculation
!> \param vxg ...
!> \param tddfpt_fac factor that multiplies the crossterms (tddfpt triplets
!>        on a closed shell system it should be -1, defaults to 1)
!> \param compute_virial ...
!> \param virial_xc ...
!> \note
!>      The old version of this routine was smarter: it handled split_desc(1)
!>      and split_desc(2) separately, thus the code automatically handled all
!>      possible cross terms (you only had to check if it was diagonal to avoid
!>      double counting). I think that is the way to go if you want to add more
!>      terms (tau,rho in LSD,...). The problem with the old code was that it
!>      because of the old functional structure it sometime guessed wrongly
!>      which derivative was where. There were probably still bugs with gradient
!>      corrected functionals (never tested), and it didn't contain first
!>      derivatives with respect to drho (that contribute also to the second
!>      derivative wrt. rho).
!>      The code was a little complex because it really tried to handle any
!>      functional derivative in the most efficient way with the given contents of
!>      rho_set.
!>      Anyway I strongly encourage whoever wants to modify this code to give a
!>      look to the old version. [fawzi]
! **************************************************************************************************
   SUBROUTINE xc_calc_2nd_deriv_analytical(v_xc, v_xc_tau, deriv_set, rho_set, rho1_set, &
                                           pw_pool, xc_section, gapw, vxg, tddfpt_fac, &
                                           compute_virial, virial_xc)

      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc, v_xc_tau
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: rho_set, rho1_set
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: gapw
      REAL(kind=dp), DIMENSION(:, :, :, :), OPTIONAL, &
         POINTER                                         :: vxg
      REAL(kind=dp), INTENT(in), OPTIONAL                :: tddfpt_fac
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc

      CHARACTER(len=*), PARAMETER :: routineN = 'xc_calc_2nd_deriv_analytical'
      INTEGER, DIMENSION(3, 3), PARAMETER :: nd = RESHAPE((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/))

      INTEGER                                            :: handle, i, ia, idir, ir, ispin, j, jdir, &
                                                            k, nspins, xc_deriv_method_id
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: gradient_f, lsd, my_compute_virial, &
                                                            my_gapw, tau_f, laplace_f
      REAL(KIND=dp)                                      :: fac, gradient_cut, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dr1dr, dra1dra, drb1drb
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: deriv_data, e_drhoa, e_drhob, &
                                                            e_drho, norm_drho, norm_drhoa, &
                                                            norm_drhob, rho1, rho1a, rho1b, &
                                                            tau1, tau1a, tau1b, laplace1, laplace1a, laplace1b, &
                                                            rho, rhoa, rhob
      TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER        :: drho, drho1, drho1a, drho1b, drhoa, drhob
      TYPE(pw_p_type), DIMENSION(:), ALLOCATABLE         :: v_drhoa, v_drhob, v_drho, v_drho_r, v_laplace
      TYPE(pw_type), POINTER                             :: tmp_g, virial_pw
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      CALL timeset(routineN, handle)

      NULLIFY (e_drhoa, e_drhob, e_drho, tmp_g)

      my_gapw = .FALSE.
      IF (PRESENT(gapw)) my_gapw = gapw

      my_compute_virial = .FALSE.
      IF (PRESENT(compute_virial)) my_compute_virial = compute_virial

      CPASSERT(ASSOCIATED(rho_set))
      CPASSERT(ASSOCIATED(rho1_set))
      CPASSERT(ASSOCIATED(v_xc))
      CPASSERT(ASSOCIATED(xc_section))
      IF (my_gapw) THEN
         CPASSERT(PRESENT(vxg))
      END IF
      IF (my_compute_virial) THEN
         CPASSERT(PRESENT(virial_xc))
      END IF

      CALL section_vals_val_get(xc_section, "XC_GRID%XC_DERIV", &
                                i_val=xc_deriv_method_id)
      CALL xc_rho_set_get(rho_set, drho_cutoff=gradient_cut)
      nspins = SIZE(v_xc)
      lsd = ASSOCIATED(rho_set%rhoa)
      fac = 0.0_dp
      IF (PRESENT(tddfpt_fac)) fac = tddfpt_fac

      bo = rho_set%local_bounds

      CALL check_for_gradients(deriv_set, lsd, gradient_f, tau_f, laplace_f)

      IF (tau_f) THEN
         CPASSERT(ASSOCIATED(v_xc_tau))
      END IF

      IF (gradient_f) THEN
         ALLOCATE (v_drho_r(nspins), v_drho(nspins))
         DO ispin = 1, nspins
            NULLIFY (v_drho_r(ispin)%pw, v_drho(ispin)%pw)
            CALL allocate_pw(v_drho(ispin)%pw, pw_pool, bo)
            CALL allocate_pw(v_drho_r(ispin)%pw, pw_pool, bo)
         END DO

         IF ((xc_deriv_method_id == xc_deriv_pw .OR. &
              xc_deriv_method_id == xc_deriv_spline3 .OR. &
              xc_deriv_method_id == xc_deriv_spline2) .AND. .NOT. my_gapw) THEN
            IF (ASSOCIATED(pw_pool)) THEN
               CALL pw_pool_create_pw(pw_pool, tmp_g, &
                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            ELSE
               ! remember to refix for gapw
               CPABORT("XC_DERIV method is not implemented in GAPW")
            END IF
         END IF
      END IF

      IF (laplace_f .AND. my_gapw) &
         CPABORT("Laplace-dependent functional not implemented with GAPW!")

      IF (laplace_f .AND. my_compute_virial) THEN
         CPABORT("Laplace-dependent functionals not implemented with stress tensors!")
      END IF

      IF (laplace_f .AND. .NOT. ASSOCIATED(tmp_g)) THEN
         IF ((xc_deriv_method_id == xc_deriv_pw .OR. &
              xc_deriv_method_id == xc_deriv_spline3 .OR. &
              xc_deriv_method_id == xc_deriv_spline2) .AND. .NOT. my_gapw) THEN
            IF (ASSOCIATED(pw_pool)) THEN
               CALL pw_pool_create_pw(pw_pool, tmp_g, &
                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            ELSE
               ! remember to refix for gapw
               CPABORT("XC_DERIV method is not implemented in GAPW")
            END IF
         END IF
      END IF

      IF (my_compute_virial .AND. (gradient_f .OR. laplace_f)) CALL allocate_pw(virial_pw, pw_pool, bo)

      IF (lsd) THEN

         !-------------------!
         ! UNrestricted case !
         !-------------------!

         CALL xc_rho_set_get(rho1_set, rhoa=rho1a, rhob=rho1b)

         IF (gradient_f) THEN
            CALL xc_rho_set_get(rho_set, drhoa=drhoa, drhob=drhob, &
                                norm_drho=norm_drho, norm_drhoa=norm_drhoa, norm_drhob=norm_drhob)
            CALL xc_rho_set_get(rho1_set, drhoa=drho1a, drhob=drho1b)

            NULLIFY (drho, drho1)
            CALL calc_drho_from_ab(drho, drhoa, drhob)
            CALL calc_drho_from_ab(drho1, drho1a, drho1b)

            CALL prepare_dr1dr(dra1dra, drhoa, drho1a)
            IF (nspins /= 1) THEN
               CALL prepare_dr1dr(drb1drb, drhob, drho1b)
               CALL prepare_dr1dr(dr1dr, drho, drho1)
            ELSE
               CALL prepare_dr1dr_ab(dr1dr, drhoa, drhob, drho1a, drho1b, fac)
            END IF

            ALLOCATE (v_drhoa(nspins), v_drhob(nspins))
            DO ispin = 1, nspins
               NULLIFY (v_drhoa(ispin)%pw, v_drhob(ispin)%pw)
            END DO
            DO ispin = 1, nspins
               CALL allocate_pw(v_drhoa(ispin)%pw, pw_pool, bo)
               CALL allocate_pw(v_drhob(ispin)%pw, pw_pool, bo)
            END DO

         END IF

         IF (laplace_f) THEN
            CALL xc_rho_set_get(rho1_set, laplace_rhoa=laplace1a, laplace_rhob=laplace1b)

            ALLOCATE (v_laplace(nspins))
            DO ispin = 1, nspins
               NULLIFY (v_laplace(ispin)%pw)
               CALL allocate_pw(v_laplace(ispin)%pw, pw_pool, bo)
            END DO

            IF (my_compute_virial) CALL xc_rho_set_get(rho_set, rhoa=rhoa, rhob=rhob)
         END IF

         IF (tau_f) THEN
            CALL xc_rho_set_get(rho1_set, tau_a=tau1a, tau_b=tau1b)
         END IF

         IF (nspins /= 1) THEN

            $:add_2nd_derivative_terms(arguments_openshell)

         ELSE

            $:add_2nd_derivative_terms(arguments_triplet_outer, arguments_triplet_inner)

         END IF

         IF (gradient_f) THEN

            IF (my_compute_virial) THEN
               CALL virial_drho_drho(virial_pw, drhoa, v_drhoa(1), virial_xc)
               CALL virial_drho_drho(virial_pw, drhob, v_drhob(2), virial_xc)
               DO idir = 1, 3
!$OMP          PARALLEL WORKSHARE DEFAULT(NONE) SHARED(drho,idir,v_drho,virial_pw)
                  virial_pw%cr3d(:, :, :) = drho(idir)%array(:, :, :)*(v_drho(1)%pw%cr3d(:, :, :) + v_drho(2)%pw%cr3d(:, :, :))
!$OMP          END PARALLEL WORKSHARE
                  DO jdir = 1, idir
                     tmp = -0.5_dp*virial_pw%pw_grid%dvol*accurate_dot_product(virial_pw%cr3d(:, :, :), &
                                                                               drho(jdir)%array(:, :, :))
                     virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
                     virial_xc(idir, jdir) = virial_xc(jdir, idir)
                  END DO
               END DO
            END IF ! my_compute_virial

            IF (my_gapw) THEN
!$OMP            PARALLEL DO PRIVATE(ia,idir,ispin,ir) DEFAULT(NONE) &
!$OMP            SHARED(bo,nspins,vxg,drhoa,drhob,v_drhoa,v_drhob,v_drho, &
!$OMP                   e_drhoa,e_drhob,e_drho,drho1a,drho1b,fac,drho,drho1) COLLAPSE(3)
               DO ir = bo(1, 2), bo(2, 2)
                  DO ia = bo(1, 1), bo(2, 1)
                     DO idir = 1, 3
                        DO ispin = 1, nspins
                           vxg(idir, ia, ir, ispin) = &
                              v_drhoa(ispin)%pw%cr3d(ia, ir, 1)*drhoa(idir)%array(ia, ir, 1) + &
                              v_drhob(ispin)%pw%cr3d(ia, ir, 1)*drhob(idir)%array(ia, ir, 1) + &
                              v_drho(ispin)%pw%cr3d(ia, ir, 1)*drho(idir)%array(ia, ir, 1)
                        END DO
                        IF (ASSOCIATED(e_drhoa)) THEN
                           vxg(idir, ia, ir, 1) = vxg(idir, ia, ir, 1) - &
                                                  e_drhoa(ia, ir, 1)*drho1a(idir)%array(ia, ir, 1)
                        END IF
                        IF (nspins /= 1 .AND. ASSOCIATED(e_drhob)) THEN
                           vxg(idir, ia, ir, 2) = vxg(idir, ia, ir, 2) - &
                                                  e_drhob(ia, ir, 1)*drho1b(idir)%array(ia, ir, 1)
                        END IF
                        IF (ASSOCIATED(e_drho)) THEN
                           IF (nspins /= 1) THEN
                              vxg(idir, ia, ir, 1) = vxg(idir, ia, ir, 1) - &
                                                     e_drho(ia, ir, 1)*drho1(idir)%array(ia, ir, 1)
                              vxg(idir, ia, ir, 2) = vxg(idir, ia, ir, 2) - &
                                                     e_drho(ia, ir, 1)*drho1(idir)%array(ia, ir, 1)
                           ELSE
                              vxg(idir, ia, ir, 1) = vxg(idir, ia, ir, 1) - &
                                                     e_drho(ia, ir, 1)*(drho1a(idir)%array(ia, ir, 1) + &
                                                                        fac*drho1b(idir)%array(ia, ir, 1))
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            ELSE

               ! partial integration
               DO idir = 1, 3

                  DO ispin = 1, nspins
!$OMP                PARALLEL WORKSHARE DEFAULT(NONE) &
!$OMP                SHARED(v_drho_r,v_drhoa,v_drhob,v_drho,drhoa,drhob,drho,ispin,idir)
                     v_drho_r(ispin)%pw%cr3d(:, :, :) = &
                        v_drhoa(ispin)%pw%cr3d(:, :, :)*drhoa(idir)%array(:, :, :) + &
                        v_drhob(ispin)%pw%cr3d(:, :, :)*drhob(idir)%array(:, :, :) + &
                        v_drho(ispin)%pw%cr3d(:, :, :)*drho(idir)%array(:, :, :)
!$OMP                END PARALLEL WORKSHARE
                  END DO
                  IF (ASSOCIATED(e_drhoa)) THEN
!$OMP                PARALLEL WORKSHARE DEFAULT(NONE) &
!$OMP                SHARED(v_drho_r,e_drhoa,drho1a,idir)
                     v_drho_r(1)%pw%cr3d(:, :, :) = v_drho_r(1)%pw%cr3d(:, :, :) - &
                                                    e_drhoa(:, :, :)*drho1a(idir)%array(:, :, :)
!$OMP                END PARALLEL WORKSHARE
                  END IF
                  IF (nspins /= 1 .AND. ASSOCIATED(e_drhob)) THEN
!$OMP                PARALLEL WORKSHARE DEFAULT(NONE)&
!$OMP                SHARED(v_drho_r,e_drhob,drho1b,idir)
                     v_drho_r(2)%pw%cr3d(:, :, :) = v_drho_r(2)%pw%cr3d(:, :, :) - &
                                                    e_drhob(:, :, :)*drho1b(idir)%array(:, :, :)
!$OMP                END PARALLEL WORKSHARE
                  END IF
                  IF (ASSOCIATED(e_drho)) THEN
!$OMP                  PARALLEL DO PRIVATE(k,j,i) DEFAULT(NONE)&
!$OMP                  SHARED(bo,v_drho_r,e_drho,drho1a,drho1b,drho1,fac,idir,nspins) COLLAPSE(3)
                     DO k = bo(1, 3), bo(2, 3)
                        DO j = bo(1, 2), bo(2, 2)
                           DO i = bo(1, 1), bo(2, 1)
                              IF (nspins /= 1) THEN
                                 v_drho_r(1)%pw%cr3d(i, j, k) = v_drho_r(1)%pw%cr3d(i, j, k) - &
                                                                e_drho(i, j, k)*drho1(idir)%array(i, j, k)
                                 v_drho_r(2)%pw%cr3d(i, j, k) = v_drho_r(2)%pw%cr3d(i, j, k) - &
                                                                e_drho(i, j, k)*drho1(idir)%array(i, j, k)
                              ELSE
                                 v_drho_r(1)%pw%cr3d(i, j, k) = v_drho_r(1)%pw%cr3d(i, j, k) - &
                                                                e_drho(i, j, k)*(drho1a(idir)%array(i, j, k) + &
                                                                                 fac*drho1b(idir)%array(i, j, k))
                              END IF
                           END DO
                        END DO
                     END DO
                  END IF

                  DO ispin = 1, nspins

                     SELECT CASE (xc_deriv_method_id)

                     CASE (xc_deriv_pw)
                        CALL pw_transfer(v_drho_r(ispin)%pw, tmp_g)
                        CALL pw_derive(tmp_g, nd(:, idir))
                        CALL pw_transfer(tmp_g, v_drho_r(ispin)%pw)
                        CALL pw_axpy(v_drho_r(ispin)%pw, v_xc(ispin)%pw)

                     CASE (xc_deriv_spline2)
                        CALL pw_transfer(v_drho_r(ispin)%pw, tmp_g)
                        CALL pw_spline2_interpolate_values_g(tmp_g)
                        CALL pw_spline2_deriv_g(tmp_g, idir=idir)
                        CALL pw_transfer(tmp_g, v_drho_r(ispin)%pw)
                        CALL pw_axpy(v_drho_r(ispin)%pw, v_xc(ispin)%pw)

                     CASE (xc_deriv_spline3)
                        CALL pw_transfer(v_drho_r(ispin)%pw, tmp_g)
                        CALL pw_spline3_interpolate_values_g(tmp_g)
                        CALL pw_spline3_deriv_g(tmp_g, idir=idir)
                        CALL pw_transfer(tmp_g, v_drho_r(ispin)%pw)
                        CALL pw_axpy(v_drho_r(ispin)%pw, v_xc(ispin)%pw)

                     CASE (xc_deriv_spline2_smooth)
                        CALL pw_nn_deriv_r(pw_in=v_drho_r(ispin)%pw, &
                                           pw_out=v_xc(ispin)%pw, coeffs=spline2_deriv_coeffs, &
                                           idir=idir)
                     CASE (xc_deriv_spline3_smooth)
                        CALL pw_nn_deriv_r(pw_in=v_drho_r(ispin)%pw, &
                                           pw_out=v_xc(ispin)%pw, coeffs=spline3_deriv_coeffs, &
                                           idir=idir)
                     CASE (xc_deriv_nn10_smooth)
                        CALL pw_nn_deriv_r(pw_in=v_drho_r(ispin)%pw, &
                                           pw_out=v_xc(ispin)%pw, coeffs=nn10_deriv_coeffs, &
                                           idir=idir)
                     CASE (xc_deriv_nn50_smooth)
                        CALL pw_nn_deriv_r(pw_in=v_drho_r(ispin)%pw, &
                                           pw_out=v_xc(ispin)%pw, coeffs=nn50_deriv_coeffs, &
                                           idir=idir)
                     CASE default
                        CPABORT("XC_DERIV method is not implemented")
                     END SELECT
                  END DO ! ispin

               END DO ! idir

            END IF

            DO idir = 1, 3
               DEALLOCATE (drho(idir)%array)
               DEALLOCATE (drho1(idir)%array)
            END DO
            DEALLOCATE (drho, drho1)

            DO ispin = 1, nspins
               CALL deallocate_pw(v_drhoa(ispin)%pw, pw_pool)
               CALL deallocate_pw(v_drhob(ispin)%pw, pw_pool)
            END DO

            DEALLOCATE (v_drhoa, v_drhob)

         END IF ! gradient_f

      ELSE

         !-----------------!
         ! restricted case !
         !-----------------!

         CALL xc_rho_set_get(rho1_set, rho=rho1)

         IF (gradient_f) THEN
            CALL xc_rho_set_get(rho_set, drho=drho, norm_drho=norm_drho)
            CALL xc_rho_set_get(rho1_set, drho=drho1)
            CALL prepare_dr1dr(dr1dr, drho, drho1)
         END IF

         IF (laplace_f) THEN
            CALL xc_rho_set_get(rho1_set, laplace_rho=laplace1)

            ALLOCATE (v_laplace(nspins))
            DO ispin = 1, nspins
               NULLIFY (v_laplace(ispin)%pw)
               CALL allocate_pw(v_laplace(ispin)%pw, pw_pool, bo)
            END DO

            IF (my_compute_virial) CALL xc_rho_set_get(rho_set, rho=rho)
         END IF

         IF (tau_f) THEN
            CALL xc_rho_set_get(rho1_set, tau=tau1)
         END IF

         $:add_2nd_derivative_terms(arguments_closedshell)

         IF (gradient_f) THEN

            IF (my_compute_virial) THEN
               CALL virial_drho_drho(virial_pw, drho, v_drho(1), virial_xc)
            END IF ! my_compute_virial

            IF (my_gapw) THEN

               DO idir = 1, 3
!$OMP               PARALLEL DO PRIVATE(ia,ir) DEFAULT(NONE) &
!$OMP               SHARED(bo,vxg,drho,v_drho,e_drho,drho1,idir) COLLAPSE(2)
                  DO ia = bo(1, 1), bo(2, 1)
                     DO ir = bo(1, 2), bo(2, 2)
                        vxg(idir, ia, ir, 1) = drho(idir)%array(ia, ir, 1)*v_drho(1)%pw%cr3d(ia, ir, 1)
                        IF (ASSOCIATED(e_drho)) THEN
                           vxg(idir, ia, ir, 1) = vxg(idir, ia, ir, 1) - drho1(idir)%array(ia, ir, 1)*e_drho(ia, ir, 1)
                        END IF
                     END DO
                  END DO
               END DO

            ELSE
               ! partial integration
               DO idir = 1, 3
!$OMP             PARALLEL WORKSHARE DEFAULT(NONE)&
!$OMP             SHARED(v_drho_r,drho,v_drho,drho1,e_drho,idir)
                  v_drho_r(1)%pw%cr3d(:, :, :) = drho(idir)%array(:, :, :)*v_drho(1)%pw%cr3d(:, :, :) - &
                                                 drho1(idir)%array(:, :, :)*e_drho(:, :, :)
!$OMP             END PARALLEL WORKSHARE

                  SELECT CASE (xc_deriv_method_id)

                  CASE (xc_deriv_pw)
                     CALL pw_transfer(v_drho_r(1)%pw, tmp_g)
                     CALL pw_derive(tmp_g, nd(:, idir))
                     CALL pw_transfer(tmp_g, v_drho_r(1)%pw)
                     CALL pw_axpy(v_drho_r(1)%pw, v_xc(1)%pw)

                  CASE (xc_deriv_spline2)
                     CALL pw_transfer(v_drho_r(1)%pw, tmp_g)
                     CALL pw_spline2_interpolate_values_g(tmp_g)
                     CALL pw_spline2_deriv_g(tmp_g, idir=idir)
                     CALL pw_transfer(tmp_g, v_drho_r(1)%pw)
                     CALL pw_axpy(v_drho_r(1)%pw, v_xc(1)%pw)

                  CASE (xc_deriv_spline3)
                     CALL pw_transfer(v_drho_r(1)%pw, tmp_g)
                     CALL pw_spline3_interpolate_values_g(tmp_g)
                     CALL pw_spline3_deriv_g(tmp_g, idir=idir)
                     CALL pw_transfer(tmp_g, v_drho_r(1)%pw)
                     CALL pw_axpy(v_drho_r(1)%pw, v_xc(1)%pw)

                  CASE (xc_deriv_spline2_smooth)
                     CALL pw_nn_deriv_r(pw_in=v_drho_r(1)%pw, &
                                        pw_out=v_xc(1)%pw, coeffs=spline2_deriv_coeffs, &
                                        idir=idir)
                  CASE (xc_deriv_spline3_smooth)
                     CALL pw_nn_deriv_r(pw_in=v_drho_r(1)%pw, &
                                        pw_out=v_xc(1)%pw, coeffs=spline3_deriv_coeffs, &
                                        idir=idir)
                  CASE (xc_deriv_nn10_smooth)
                     CALL pw_nn_deriv_r(pw_in=v_drho_r(1)%pw, &
                                        pw_out=v_xc(1)%pw, coeffs=nn10_deriv_coeffs, &
                                        idir=idir)
                  CASE (xc_deriv_nn50_smooth)
                     CALL pw_nn_deriv_r(pw_in=v_drho_r(1)%pw, &
                                        pw_out=v_xc(1)%pw, coeffs=nn50_deriv_coeffs, &
                                        idir=idir)
                  CASE default
                     CPABORT("XC_DERIV method is not implemented")
                  END SELECT

               END DO
            END IF

         END IF

      END IF

      IF (laplace_f) THEN

         SELECT CASE (xc_deriv_method_id)

         CASE (xc_deriv_pw)
            DO ispin = 1, nspins
               CALL pw_transfer(v_laplace(ispin)%pw, tmp_g)
               CALL pw_laplace(tmp_g)
               CALL pw_transfer(tmp_g, v_laplace(ispin)%pw)
               CALL pw_axpy(v_laplace(ispin)%pw, v_xc(ispin)%pw)
            END DO
         CASE default
            CPABORT("XC_DERIV method not implemented with Laplace-dependent functionals!")
         END SELECT
      END IF

      IF (gradient_f) THEN

         DO ispin = 1, nspins
            CALL deallocate_pw(v_drho(ispin)%pw, pw_pool)
            CALL deallocate_pw(v_drho_r(ispin)%pw, pw_pool)
         END DO
         DEALLOCATE (v_drho, v_drho_r)

      END IF

      IF (laplace_f) THEN
      DO ispin = 1, nspins
         CALL deallocate_pw(v_laplace(ispin)%pw, pw_pool)
      END DO
      DEALLOCATE (v_laplace)
      END IF

      IF (ASSOCIATED(tmp_g) .AND. ASSOCIATED(pw_pool)) THEN
         CALL pw_pool_give_back_pw(pw_pool, tmp_g)
      END IF

      IF (my_compute_virial .AND. (gradient_f .OR. laplace_f)) THEN
         CALL deallocate_pw(virial_pw, pw_pool)
      END IF

      CALL timestop(handle)
   END SUBROUTINE xc_calc_2nd_deriv_analytical

! **************************************************************************************************
!> \brief allocates grids using pw_pool (if associated) or with bounds
!> \param pw ...
!> \param pw_pool ...
!> \param bo ...
! **************************************************************************************************
   SUBROUTINE allocate_pw(pw, pw_pool, bo)
      TYPE(pw_type), INTENT(INOUT), POINTER              :: pw
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: pw_pool
      INTEGER, DIMENSION(2, 3), INTENT(IN)               :: bo

      IF (ASSOCIATED(pw_pool)) THEN
         CALL pw_pool_create_pw(pw_pool, pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(pw)
      ELSE
         NULLIFY (pw)
         ALLOCATE (pw)
         ALLOCATE (pw%cr3d(bo(1, 1):bo(2, 1), bo(1, 2):bo(2, 2), bo(1, 3):bo(2, 3)))
         pw%cr3d = 0.0_dp
      END IF

   END SUBROUTINE allocate_pw

! **************************************************************************************************
!> \brief deallocates grid allocated with allocate_pw
!> \param pw ...
!> \param pw_pool ...
! **************************************************************************************************
   SUBROUTINE deallocate_pw(pw, pw_pool)
      TYPE(pw_type), INTENT(INOUT), POINTER              :: pw
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: pw_pool

      IF (ASSOCIATED(pw_pool)) THEN
         CALL pw_pool_give_back_pw(pw_pool, pw)
      ELSE
         DEALLOCATE (pw%cr3d)
         DEALLOCATE (pw)
      END IF

   END SUBROUTINE deallocate_pw

! **************************************************************************************************
!> \brief updates virial from first derivative w.r.t. norm_drho
!> \param virial_pw ...
!> \param drho ...
!> \param drho1 ...
!> \param deriv_data ...
!> \param virial_xc ...
! **************************************************************************************************
   SUBROUTINE virial_drho_drho1(virial_pw, drho, drho1, deriv_data, virial_xc)
      TYPE(pw_type), INTENT(IN)                          :: virial_pw
      TYPE(cp_3d_r_p_type), DIMENSION(3), INTENT(IN)     :: drho, drho1
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: deriv_data
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc

      INTEGER                                            :: idir, jdir
      REAL(KIND=dp)                                      :: tmp

      DO idir = 1, 3
!$OMP    PARALLEL WORKSHARE DEFAULT(NONE)&
!$OMP    SHARED(drho,idir,virial_pw,deriv_data)
         virial_pw%cr3d(:, :, :) = drho(idir)%array(:, :, :)*deriv_data(:, :, :)
!$OMP    END PARALLEL WORKSHARE
         DO jdir = 1, 3
            tmp = virial_pw%pw_grid%dvol*accurate_dot_product(virial_pw%cr3d(:, :, :), &
                                                              drho1(jdir)%array(:, :, :))
            virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
            virial_xc(idir, jdir) = virial_xc(idir, jdir) + tmp
         END DO
      END DO

   END SUBROUTINE virial_drho_drho1

! **************************************************************************************************
!> \brief Adds virial contribution from second order potential parts
!> \param virial_pw ...
!> \param drho ...
!> \param v_drho ...
!> \param virial_xc ...
! **************************************************************************************************
   SUBROUTINE virial_drho_drho(virial_pw, drho, v_drho, virial_xc)
      TYPE(pw_type), INTENT(IN)                          :: virial_pw
      TYPE(cp_3d_r_p_type), DIMENSION(3), INTENT(IN)     :: drho
      TYPE(pw_p_type), INTENT(IN)                        :: v_drho
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT)      :: virial_xc

      INTEGER                                            :: idir, jdir
      REAL(KIND=dp)                                      :: tmp

      DO idir = 1, 3
!$OMP PARALLEL WORKSHARE DEFAULT(NONE)&
!$OMP SHARED(drho,idir,v_drho,virial_pw)
         virial_pw%cr3d(:, :, :) = drho(idir)%array(:, :, :)*v_drho%pw%cr3d(:, :, :)
!$OMP END PARALLEL WORKSHARE
         DO jdir = 1, idir
            tmp = -virial_pw%pw_grid%dvol*accurate_dot_product(virial_pw%cr3d(:, :, :), &
                                                               drho(jdir)%array(:, :, :))
            virial_xc(jdir, idir) = virial_xc(jdir, idir) + tmp
            virial_xc(idir, jdir) = virial_xc(jdir, idir)
         END DO
      END DO

   END SUBROUTINE virial_drho_drho

! **************************************************************************************************
!> \brief ...
!> \param rho_r ...
!> \param pw_pool ...
!> \param virial_xc ...
!> \param deriv_data ...
! **************************************************************************************************
   SUBROUTINE virial_laplace(rho_r, pw_pool, virial_xc, deriv_data)
      TYPE(pw_type), TARGET :: rho_r
      TYPE(pw_pool_type), POINTER, INTENT(IN) :: pw_pool
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT) :: virial_xc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN) :: deriv_data

      CHARACTER(len=*), PARAMETER                        :: routineN = 'virial_laplace'

      INTEGER :: handle, idir, jdir
      TYPE(pw_type), POINTER :: virial_pw, tmp_g
      INTEGER, DIMENSION(3) :: my_deriv

      CALL timeset(routineN, handle)

      NULLIFY (virial_pw, tmp_g)

      CALL pw_pool_create_pw(pw_pool, virial_pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(pw_pool, tmp_g, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_zero(virial_pw)

      DO idir = 1, 3
         DO jdir = idir, 3
            CALL pw_transfer(rho_r, tmp_g)

            my_deriv = 0
            my_deriv(idir) = 1
            my_deriv(jdir) = my_deriv(jdir) + 1

            CALL pw_derive(tmp_g, my_deriv)
            CALL pw_transfer(tmp_g, virial_pw)
            virial_xc(idir, jdir) = virial_xc(idir, jdir) - 2.0_dp*virial_pw%pw_grid%dvol* &
                                    accurate_dot_product(virial_pw%cr3d(:, :, :), &
                                                         deriv_data(:, :, :))
            virial_xc(jdir, idir) = virial_xc(idir, jdir)
         END DO
      END DO
      CALL pw_pool_give_back_pw(pw_pool, virial_pw)
      CALL pw_pool_give_back_pw(pw_pool, tmp_g)

      CALL timestop(handle)

   END SUBROUTINE virial_laplace

! **************************************************************************************************
!> \brief Prepare objects for the calculation of the 2nd derivatives of the density functional.
!>      The calculation must then be performed with xc_calc_2nd_deriv.
!> \param deriv_set object containing the XC derivatives (out)
!> \param rho_set object that will contain the density at which the
!>        derivatives were calculated
!> \param rho_r the place where you evaluate the derivative
!> \param pw_pool the pool for the grids
!> \param xc_section which functional should be used and how to calculate it
!> \param tau_r kinetic energy density in real space
! **************************************************************************************************
   SUBROUTINE xc_prep_2nd_deriv(deriv_set, &
                                rho_set, rho_r, pw_pool, xc_section, tau_r)

      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(pw_p_type), DIMENSION(:), OPTIONAL, POINTER   :: tau_r

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xc_prep_2nd_deriv'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: lsd
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r_pw, tau

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(rho_set))
      CPASSERT(ASSOCIATED(xc_section))
      CPASSERT(ASSOCIATED(pw_pool))

      nspins = SIZE(rho_r)
      lsd = (nspins /= 1)

      ALLOCATE (rho_r_pw(nspins))
      DO ispin = 1, nspins
         rho_r_pw(ispin)%pw => rho_r(ispin)%pw
      END DO

      NULLIFY (rho_g, tau)
      IF (PRESENT(tau_r)) &
         tau => tau_r

      IF (section_get_lval(xc_section, "2ND_DERIV_ANALYTICAL")) THEN
         CALL xc_rho_set_and_dset_create(rho_set, deriv_set, 2, &
                                         rho_r_pw, rho_g, tau, xc_section, pw_pool, &
                                         needs_basic_components=.TRUE.)

         CALL divide_by_norm_drho(deriv_set, rho_set, lsd)
      ELSE
         CALL xc_rho_set_and_dset_create(rho_set, deriv_set, 1, &
                                         rho_r_pw, rho_g, tau, xc_section, pw_pool, &
                                         needs_basic_components=.TRUE.)
      END IF

      DEALLOCATE (rho_r_pw)

      CALL timestop(handle)

   END SUBROUTINE xc_prep_2nd_deriv

! **************************************************************************************************
!> \brief divides derivatives from deriv_set by norm_drho
!> \param deriv_set ...
!> \param rho_set ...
!> \param lsd ...
! **************************************************************************************************
   SUBROUTINE divide_by_norm_drho(deriv_set, rho_set, lsd)

      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      LOGICAL, INTENT(IN)                                :: lsd

      CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH)          :: desc
      CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), &
         POINTER                                         :: split_desc
      INTEGER                                            :: idesc, order
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: drho_cutoff
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: norm_drho, norm_drhoa, norm_drhob
      TYPE(cp_sll_xc_deriv_type), POINTER                :: pos
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

! check for unknown derivatives and divide by norm_drho where necessary

      bo = rho_set%local_bounds
      CALL xc_rho_set_get(rho_set, drho_cutoff=drho_cutoff, norm_drho=norm_drho, &
                          norm_drhoa=norm_drhoa, norm_drhob=norm_drhob, can_return_null=.TRUE.)

      pos => deriv_set%derivs
      DO WHILE (cp_sll_xc_deriv_next(pos, el_att=deriv_att))
         CALL xc_derivative_get(deriv_att, order=order, &
                                desc=desc, split_desc=split_desc)
         IF (order == 1 .OR. order == 2) THEN
            DO idesc = 1, SIZE(split_desc)
               SELECT CASE (split_desc(idesc))
               CASE ("norm_drho")
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_att,norm_drho,drho_cutoff)
                  deriv_att%deriv_data(:, :, :) = deriv_att%deriv_data(:, :, :)/ &
                                                  MAX(norm_drho(:, :, :), drho_cutoff)
!$OMP END PARALLEL WORKSHARE
               CASE ("norm_drhoa")
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_att,norm_drhoa,drho_cutoff)
                  deriv_att%deriv_data(:, :, :) = deriv_att%deriv_data(:, :, :)/ &
                                                  MAX(norm_drhoa(:, :, :), drho_cutoff)
!$OMP END PARALLEL WORKSHARE
               CASE ("norm_drhob")
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(deriv_att,norm_drhob,drho_cutoff)
                  deriv_att%deriv_data(:, :, :) = deriv_att%deriv_data(:, :, :)/ &
                                                  MAX(norm_drhob(:, :, :), drho_cutoff)
!$OMP END PARALLEL WORKSHARE
               CASE ("rho")
                  IF (lsd) &
                     CPABORT("rho not handled in lsd: '"//TRIM(desc)//"'")
               CASE ("tau")
                  IF (lsd) &
                     CPABORT("tau not handled in lsd: '"//TRIM(desc)//"'")
               CASE ("laplace_rho")
                  IF (lsd) &
                     CPABORT("laplace_rho not handled in lsd: '"//TRIM(desc)//"'")
               CASE ("rhoa", "rhob", "tau_a", "tau_b", "laplace_rhoa", "laplace_rhob")
               CASE default
                  CALL cp_abort(__LOCATION__, &
                                "unhandled derivative: '"// &
                                TRIM(split_desc(idesc))//"' in '"// &
                                TRIM(desc)//"'")
               END SELECT
            END DO
         END IF
      END DO

   END SUBROUTINE divide_by_norm_drho

! **************************************************************************************************
!> \brief allocates and calculates drho from given spin densities drhoa, drhob
!> \param drho ...
!> \param drhoa ...
!> \param drhob ...
! **************************************************************************************************
   SUBROUTINE calc_drho_from_ab(drho, drhoa, drhob)
      TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER        :: drho
      TYPE(cp_3d_r_p_type), DIMENSION(3), INTENT(IN)     :: drhoa, drhob

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_drho_from_ab'

      INTEGER                                            :: handle, idir

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(drho))

      ALLOCATE (drho(3))
      DO idir = 1, 3
         NULLIFY (drho(idir)%array)
         ALLOCATE (drho(idir)%array(LBOUND(drhoa(1)%array, 1):UBOUND(drhoa(1)%array, 1), &
                                    LBOUND(drhoa(1)%array, 2):UBOUND(drhoa(1)%array, 2), &
                                    LBOUND(drhoa(1)%array, 3):UBOUND(drhoa(1)%array, 3)))
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(drho,drhoa,drhob,idir)
         drho(idir)%array(:, :, :) = drhoa(idir)%array(:, :, :) + drhob(idir)%array(:, :, :)
!$OMP END PARALLEL WORKSHARE
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief allocates and calculates dot products of two density gradients
!> \param dr1dr ...
!> \param drho ...
!> \param drho1 ...
! **************************************************************************************************
   SUBROUTINE prepare_dr1dr(dr1dr, drho, drho1)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: dr1dr
      TYPE(cp_3d_r_p_type), DIMENSION(3), INTENT(IN)     :: drho, drho1

      CHARACTER(len=*), PARAMETER                        :: routineN = 'prepare_dr1dr'

      INTEGER                                            :: handle, idir

      CALL timeset(routineN, handle)

      ALLOCATE (dr1dr(LBOUND(drho(1)%array, 1):UBOUND(drho(1)%array, 1), &
                      LBOUND(drho(1)%array, 2):UBOUND(drho(1)%array, 2), &
                      LBOUND(drho(1)%array, 3):UBOUND(drho(1)%array, 3)))

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(dr1dr,drho,drho1)
      dr1dr(:, :, :) = drho(1)%array(:, :, :)*drho1(1)%array(:, :, :)
!$OMP END PARALLEL WORKSHARE
      DO idir = 2, 3
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(dr1dr,drho,drho1,idir)
         dr1dr(:, :, :) = dr1dr(:, :, :) + drho(idir)%array(:, :, :)*drho1(idir)%array(:, :, :)
!$OMP END PARALLEL WORKSHARE
      END DO

      CALL timestop(handle)

   END SUBROUTINE prepare_dr1dr

! **************************************************************************************************
!> \brief allocates and calculates dot product of two densities for triplets
!> \param dr1dr ...
!> \param drhoa ...
!> \param drhob ...
!> \param drho1a ...
!> \param drho1b ...
!> \param fac ...
! **************************************************************************************************
   SUBROUTINE prepare_dr1dr_ab(dr1dr, drhoa, drhob, drho1a, drho1b, fac)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: dr1dr
      TYPE(cp_3d_r_p_type), DIMENSION(3), INTENT(IN)     :: drhoa, drhob, drho1a, drho1b
      REAL(KIND=dp), INTENT(IN)                          :: fac

      CHARACTER(len=*), PARAMETER                        :: routineN = 'prepare_dr1dr_ab'

      INTEGER                                            :: handle, idir

      CALL timeset(routineN, handle)

      ALLOCATE (dr1dr(LBOUND(drhoa(1)%array, 1):UBOUND(drhoa(1)%array, 1), &
                      LBOUND(drhoa(1)%array, 2):UBOUND(drhoa(1)%array, 2), &
                      LBOUND(drhoa(1)%array, 3):UBOUND(drhoa(1)%array, 3)))

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(fac,dr1dr,drho1a,drho1b,drhoa,drhob)
      dr1dr(:, :, :) = drhoa(1)%array(:, :, :)*(drho1a(1)%array(:, :, :) + &
                                                fac*drho1b(1)%array(:, :, :)) + &
                       drhob(1)%array(:, :, :)*(fac*drho1a(1)%array(:, :, :) + &
                                                drho1b(1)%array(:, :, :))
!$OMP END PARALLEL WORKSHARE
      DO idir = 2, 3
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(fac,dr1dr,drho1a,drho1b,drhoa,drhob,idir)
         dr1dr(:, :, :) = dr1dr(:, :, :) + &
                          drhoa(idir)%array(:, :, :)*(drho1a(idir)%array(:, :, :) + &
                                                      fac*drho1b(idir)%array(:, :, :)) + &
                          drhob(idir)%array(:, :, :)*(fac*drho1a(idir)%array(:, :, :) + &
                                                      drho1b(idir)%array(:, :, :))
!$OMP END PARALLEL WORKSHARE
      END DO

      CALL timestop(handle)

   END SUBROUTINE prepare_dr1dr_ab

! **************************************************************************************************
!> \brief checks for gradients
!> \param deriv_set ...
!> \param lsd ...
!> \param gradient_f ...
!> \param tau_f ...
!> \param laplace_f ...
! **************************************************************************************************
   SUBROUTINE check_for_gradients(deriv_set, lsd, gradient_f, tau_f, laplace_f)
      TYPE(xc_derivative_set_type), INTENT(IN)           :: deriv_set
      LOGICAL, INTENT(IN)                                :: lsd
      LOGICAL, INTENT(OUT)                               :: gradient_f, tau_f, laplace_f

      CHARACTER(len=*), PARAMETER :: routineN = 'check_for_gradients'

      CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH)          :: desc
      CHARACTER(len=MAX_LABEL_LENGTH), DIMENSION(:), &
         POINTER                                         :: split_desc
      INTEGER                                            :: handle, iorder, order
      LOGICAL                                            :: unknown_deriv
      TYPE(cp_sll_xc_deriv_type), POINTER                :: pos
      TYPE(xc_derivative_type), POINTER                  :: deriv_att

      CALL timeset(routineN, handle)

      gradient_f = .FALSE.
      tau_f = .FALSE.
      laplace_f = .FALSE.
      ! check for unknown derivatives
      pos => deriv_set%derivs
      DO WHILE (cp_sll_xc_deriv_next(pos, el_att=deriv_att))
         unknown_deriv = .FALSE.
         CALL xc_derivative_get(deriv_att, order=order, &
                                desc=desc, split_desc=split_desc)
         IF (lsd) THEN
            DO iorder = 1, order
               SELECT CASE (split_desc(iorder))
               CASE ("rho", "rhoa", "rhob")
               CASE ("norm_drho", "norm_drhoa", "norm_drhob")
                  gradient_f = .TRUE.
               CASE ("tau_a", "tau_b", "tau")
                  tau_f = .TRUE.
               CASE ("laplace_rhoa", "laplace_rhob", "laplace_rho")
                  laplace_f = .TRUE.
               CASE default
                  unknown_deriv = .TRUE.
               END SELECT
            END DO
         ELSE
            DO iorder = 1, order
               SELECT CASE (split_desc(iorder))
               CASE ("rho")
               CASE ("tau")
                  tau_f = .TRUE.
               CASE ("norm_drho")
                  gradient_f = .TRUE.
               CASE ("laplace_rho")
                  laplace_f = .TRUE.
               CASE default
                  unknown_deriv = .TRUE.
               END SELECT
            END DO
         END IF
         IF (unknown_deriv) &
            CALL cp_abort(__LOCATION__, &
                          "unknown functional derivative (LSD="//TRIM(cp_to_string(lsd))// &
                          "): '"//TRIM(desc)//"'")
      END DO

      CALL timestop(handle)

   END SUBROUTINE check_for_gradients

END MODULE xc

