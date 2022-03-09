!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief contains utility functions for the xc package
!> \par History
!>      03.2022 created [F. Stein]
!> \author Frederick Stein
! **************************************************************************************************
MODULE xc_util
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_type
   USE pw_spline_utils,                 ONLY: nn10_coeffs,&
                                              nn50_coeffs,&
                                              pw_nn_smear_r,&
                                              spline2_coeffs,&
                                              spline3_coeffs
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_type
   USE xc_input_constants,              ONLY: xc_rho_nn10,&
                                              xc_rho_nn50,&
                                              xc_rho_no_smooth,&
                                              xc_rho_spline2_smooth,&
                                              xc_rho_spline3_smooth
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: xc_pw_smooth
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_util'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param pw_in ...
!> \param pw_out ...
!> \param xc_smooth_id ...
!> \param pw_pool ...
! **************************************************************************************************
   SUBROUTINE xc_pw_smooth(pw_in, pw_out, xc_smooth_id, pw_pool)
      TYPE(pw_type), POINTER                             :: pw_in, pw_out
      INTEGER, INTENT(IN)                                :: xc_smooth_id
      TYPE(pw_pool_type), INTENT(IN), OPTIONAL, POINTER  :: pw_pool

      IF (PRESENT(pw_pool)) THEN
         CALL pw_pool_create_pw(pw_pool, pw_out, &
                                use_data=REALDATA3D, in_space=REALSPACE)
      END IF

      SELECT CASE (xc_smooth_id)
      CASE (xc_rho_no_smooth)
         CALL pw_copy(pw_in, pw_out)
      CASE (xc_rho_spline2_smooth)
         CALL pw_zero(pw_out)
         CALL pw_nn_smear_r(pw_in=pw_in, &
                            pw_out=pw_out, &
                            coeffs=spline2_coeffs)
      CASE (xc_rho_spline3_smooth)
         CALL pw_zero(pw_out)
         CALL pw_nn_smear_r(pw_in=pw_in, &
                            pw_out=pw_out, &
                            coeffs=spline3_coeffs)
      CASE (xc_rho_nn10)
         CALL pw_zero(pw_out)
         CALL pw_nn_smear_r(pw_in=pw_in, &
                            pw_out=pw_out, &
                            coeffs=nn10_coeffs)
      CASE (xc_rho_nn50)
         CALL pw_zero(pw_out)
         CALL pw_nn_smear_r(pw_in=pw_in, &
                            pw_out=pw_out, &
                            coeffs=nn50_coeffs)
      CASE default
         CPABORT("Unsupported smoothing")
      END SELECT

   END SUBROUTINE xc_pw_smooth

END MODULE xc_util
