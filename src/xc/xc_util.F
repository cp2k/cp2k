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
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_laplace,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_spline_utils,                 ONLY: &
        nn10_coeffs, nn10_deriv_coeffs, nn50_coeffs, nn50_deriv_coeffs, pw_nn_deriv_r, &
        pw_nn_smear_r, pw_spline2_deriv_g, pw_spline2_interpolate_values_g, pw_spline3_deriv_g, &
        pw_spline3_interpolate_values_g, pw_spline_scale_deriv, spline2_coeffs, &
        spline2_deriv_coeffs, spline3_coeffs, spline3_deriv_coeffs
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE xc_input_constants,              ONLY: &
        xc_deriv_nn10_smooth, xc_deriv_nn50_smooth, xc_deriv_pw, xc_deriv_spline2, &
        xc_deriv_spline2_smooth, xc_deriv_spline3, xc_deriv_spline3_smooth, xc_rho_nn10, &
        xc_rho_nn50, xc_rho_no_smooth, xc_rho_spline2_smooth, xc_rho_spline3_smooth
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: xc_pw_smooth, xc_pw_laplace, xc_pw_divergence, xc_pw_derive, xc_requires_tmp_g
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_util'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param xc_deriv_id ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION xc_requires_tmp_g(xc_deriv_id) RESULT(requires)
      INTEGER, INTENT(IN)                                :: xc_deriv_id
      LOGICAL                                            :: requires

      requires = (xc_deriv_id == xc_deriv_spline2) .OR. &
                 (xc_deriv_id == xc_deriv_spline3) .OR. &
                 (xc_deriv_id == xc_deriv_pw)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param pw_in ...
!> \param pw_out ...
!> \param xc_smooth_id ...
! **************************************************************************************************
   SUBROUTINE xc_pw_smooth(pw_in, pw_out, xc_smooth_id)
      TYPE(pw_type), POINTER                             :: pw_in, pw_out
      INTEGER, INTENT(IN)                                :: xc_smooth_id

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xc_pw_smooth'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

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

      CALL timestop(handle)

   END SUBROUTINE xc_pw_smooth

! **************************************************************************************************
!> \brief Calculates the Laplacian of pw
!> \param pw on input: pw of which the Laplacian shall be calculated, on output if pw_out is absent: Laplacian of input
!> \param pw_pool ...
!> \param deriv_method_id ...
!> \param pw_out if present, save the Laplacian of pw here
!> \param tmp_g scratch grid in reciprocal space, used instead of the internal grid if given explicitly to save memory
! **************************************************************************************************
   SUBROUTINE xc_pw_laplace(pw, pw_pool, deriv_method_id, pw_out, tmp_g)
      TYPE(pw_type), INTENT(IN), POINTER                 :: pw
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: pw_pool
      INTEGER, INTENT(IN)                                :: deriv_method_id
      TYPE(pw_type), INTENT(IN), OPTIONAL, POINTER       :: pw_out, tmp_g

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xc_pw_laplace'

      INTEGER                                            :: handle
      LOGICAL                                            :: owns_tmp_g
      TYPE(pw_type), POINTER                             :: my_tmp_g

      CALL timeset(routineN, handle)

      SELECT CASE (deriv_method_id)
      CASE (xc_deriv_pw)

         NULLIFY (my_tmp_g)
         IF (PRESENT(tmp_g)) my_tmp_g => tmp_g

         owns_tmp_g = .FALSE.
         IF (.NOT. ASSOCIATED(my_tmp_g)) THEN
            CALL pw_pool_create_pw(pw_pool, my_tmp_g, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            owns_tmp_g = .TRUE.
         END IF
         CALL pw_zero(my_tmp_g)
         CALL pw_transfer(pw, my_tmp_g)

         CALL pw_laplace(my_tmp_g)

         IF (PRESENT(pw_out)) THEN
            CALL pw_transfer(my_tmp_g, pw_out)
         ELSE
            CALL pw_transfer(my_tmp_g, pw)
         END IF
         IF (owns_tmp_g) CALL pw_pool_give_back_pw(pw_pool, my_tmp_g)
      CASE default
         CPABORT("Unsupported derivative method")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE xc_pw_laplace

! **************************************************************************************************
!> \brief Calculates the divergence of pw_to_deriv
!> \param xc_deriv_method_id ...
!> \param pw_to_deriv ...
!> \param tmp_g ...
!> \param vxc_g ...
!> \param vxc_r ...
! **************************************************************************************************
   SUBROUTINE xc_pw_divergence(xc_deriv_method_id, pw_to_deriv, tmp_g, vxc_g, vxc_r)
      INTEGER, INTENT(IN)                                :: xc_deriv_method_id
      TYPE(pw_p_type), DIMENSION(3)                      :: pw_to_deriv
      TYPE(pw_type), POINTER                             :: tmp_g, vxc_g, vxc_r

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xc_pw_divergence'

      INTEGER                                            :: handle, idir

      CALL timeset(routineN, handle)

      ! partial integration
      IF (xc_deriv_method_id /= xc_deriv_pw) THEN
         CALL pw_spline_scale_deriv(pw_to_deriv, &
                                    transpose=.TRUE.)
      END IF

      IF (ASSOCIATED(vxc_g)) CALL pw_zero(vxc_g)

      DO idir = 1, 3
         CALL xc_pw_derive(pw_to_deriv(idir)%pw, tmp_g, vxc_r, idir, xc_deriv_method_id, copy_to_vxcr=.FALSE.)
         IF (ASSOCIATED(tmp_g) .AND. ASSOCIATED(vxc_g)) CALL pw_axpy(tmp_g, vxc_g)
      END DO

      IF (ASSOCIATED(vxc_g)) THEN
         CALL pw_transfer(vxc_g, pw_to_deriv(1)%pw)
         CALL pw_axpy(pw_to_deriv(1)%pw, vxc_r)
      END IF

      CALL timestop(handle)

   END SUBROUTINE xc_pw_divergence

! **************************************************************************************************
!> \brief Calculates the derivative of a function on a planewave grid in a given direction
!> \param pw function to derive
!> \param tmp_g temporary grid in reciprocal space, only required if derivative method is pw or spline
!> \param vxc_r if tmp_g is not required, add derivative here
!> \param idir direction of derivative
!> \param xc_deriv_method_id ...
!> \param copy_to_vxcr ...
!> \param pw_g ...
! **************************************************************************************************
   SUBROUTINE xc_pw_derive(pw, tmp_g, vxc_r, idir, xc_deriv_method_id, copy_to_vxcr, pw_g)
      TYPE(pw_type), INTENT(IN), POINTER                 :: pw
      TYPE(pw_type), POINTER                             :: tmp_g, vxc_r
      INTEGER, INTENT(IN)                                :: idir, xc_deriv_method_id
      LOGICAL, INTENT(IN), OPTIONAL                      :: copy_to_vxcr
      TYPE(pw_type), INTENT(IN), OPTIONAL, POINTER       :: pw_g

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xc_pw_derive'
      INTEGER, DIMENSION(3, 3), PARAMETER :: nd = RESHAPE((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3, 3/))

      INTEGER                                            :: handle
      LOGICAL                                            :: my_copy_to_vxcr

      CALL timeset(routineN, handle)

      my_copy_to_vxcr = .TRUE.
      IF (PRESENT(copy_to_vxcr)) my_copy_to_vxcr = copy_to_vxcr

      IF (xc_requires_tmp_g(xc_deriv_method_id)) THEN

         IF (PRESENT(pw_g)) THEN
            IF (ASSOCIATED(pw_g)) THEN
               CALL pw_copy(pw_g, tmp_g)
            ELSE
               CALL pw_transfer(pw, tmp_g)
            END IF
         ELSE
            CALL pw_transfer(pw, tmp_g)
         END IF

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
            CPABORT("Unsupported deriv method")
         END SELECT

         IF (my_copy_to_vxcr) CALL pw_transfer(tmp_g, vxc_r)
      ELSE
         SELECT CASE (xc_deriv_method_id)
         CASE (xc_deriv_spline2_smooth)
            CALL pw_nn_deriv_r(pw_in=pw, &
                               pw_out=vxc_r, coeffs=spline2_deriv_coeffs, &
                               idir=idir)
         CASE (xc_deriv_spline3_smooth)
            CALL pw_nn_deriv_r(pw_in=pw, &
                               pw_out=vxc_r, coeffs=spline3_deriv_coeffs, &
                               idir=idir)
         CASE (xc_deriv_nn10_smooth)
            CALL pw_nn_deriv_r(pw_in=pw, &
                               pw_out=vxc_r, coeffs=nn10_deriv_coeffs, &
                               idir=idir)
         CASE (xc_deriv_nn50_smooth)
            CALL pw_nn_deriv_r(pw_in=pw, &
                               pw_out=vxc_r, coeffs=nn50_deriv_coeffs, &
                               idir=idir)
         CASE default
            CPABORT("Unsupported derivative method")
         END SELECT
      END IF

      CALL timestop(handle)

   END SUBROUTINE xc_pw_derive

END MODULE xc_util
