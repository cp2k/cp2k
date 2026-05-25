!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for broadened DOS and PDOS output.
! **************************************************************************************************
MODULE qs_dos_utils
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dos_utils'

   INTEGER, PARAMETER, PUBLIC :: broadening_gaussian = 1, &
                                 broadening_lorentzian = 2, &
                                 broadening_pseudo_voigt = 3

   PUBLIC :: add_broadened_peak, broadening_cutoff, broadening_function, write_broadening_info

CONTAINS

! **************************************************************************************************
!> \brief Add a broadened spectral line to a DOS curve.
!> \param dos ...
!> \param occ_dos ...
!> \param emin ...
!> \param de ...
!> \param eig ...
!> \param occ ...
!> \param weight ...
!> \param broaden_type ...
!> \param broaden_width ...
!> \param voigt_mixing ...
! **************************************************************************************************
   SUBROUTINE add_broadened_peak(dos, occ_dos, emin, de, eig, occ, weight, broaden_type, &
                                 broaden_width, voigt_mixing)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: dos, occ_dos
      REAL(KIND=dp), INTENT(IN)                          :: emin, de, eig, occ, weight
      INTEGER, INTENT(IN)                                :: broaden_type
      REAL(KIND=dp), INTENT(IN)                          :: broaden_width, voigt_mixing

      INTEGER                                            :: i, iend, istart, nhist
      REAL(KIND=dp)                                      :: eval, line_shape

      nhist = SIZE(dos)
      istart = MAX(1, FLOOR((eig - broadening_cutoff(broaden_type, broaden_width) - emin)/de) + 1)
      iend = MIN(nhist, CEILING((eig + broadening_cutoff(broaden_type, broaden_width) - emin)/de) + 1)
      DO i = istart, iend
         eval = emin + (i - 1)*de
         line_shape = broadening_function(eval - eig, broaden_type, broaden_width, voigt_mixing)
         dos(i) = dos(i) + weight*line_shape
         occ_dos(i) = occ_dos(i) + weight*occ*line_shape
      END DO

   END SUBROUTINE add_broadened_peak

! **************************************************************************************************
!> \brief Broadening cutoff used for numerical accumulation.
!> \param broaden_type ...
!> \param broaden_width ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION broadening_cutoff(broaden_type, broaden_width) RESULT(cutoff)

      INTEGER, INTENT(IN)                                :: broaden_type
      REAL(KIND=dp), INTENT(IN)                          :: broaden_width
      REAL(KIND=dp)                                      :: cutoff

      IF (broaden_type == broadening_gaussian) THEN
         cutoff = 8.0_dp*broaden_width
      ELSE
         cutoff = 50.0_dp*broaden_width
      END IF

   END FUNCTION broadening_cutoff

! **************************************************************************************************
!> \brief Normalized broadening function. BROADEN_WIDTH is FWHM.
!> \param delta_e ...
!> \param broaden_type ...
!> \param broaden_width ...
!> \param voigt_mixing ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION broadening_function(delta_e, broaden_type, broaden_width, voigt_mixing) RESULT(value)

      REAL(KIND=dp), INTENT(IN)                          :: delta_e
      INTEGER, INTENT(IN)                                :: broaden_type
      REAL(KIND=dp), INTENT(IN)                          :: broaden_width, voigt_mixing
      REAL(KIND=dp)                                      :: value

      REAL(KIND=dp)                                      :: eta, gamma, gaussian_value, sigma

      IF (broaden_width <= 0.0_dp) THEN
         value = 0.0_dp
         RETURN
      END IF

      sigma = broaden_width/(2.0_dp*SQRT(2.0_dp*LOG(2.0_dp)))
      gamma = 0.5_dp*broaden_width
      gaussian_value = EXP(-0.5_dp*(delta_e/sigma)**2)/(sigma*SQRT(2.0_dp*pi))

      SELECT CASE (broaden_type)
      CASE (broadening_gaussian)
         value = gaussian_value
      CASE (broadening_lorentzian)
         value = gamma/(pi*(delta_e**2 + gamma**2))
      CASE (broadening_pseudo_voigt)
         eta = MIN(1.0_dp, MAX(0.0_dp, voigt_mixing))
         value = eta*gamma/(pi*(delta_e**2 + gamma**2)) + (1.0_dp - eta)*gaussian_value
      CASE DEFAULT
         value = gaussian_value
      END SELECT

   END FUNCTION broadening_function

! **************************************************************************************************
!> \brief Write broadening metadata.
!> \param iw ...
!> \param broaden_type ...
!> \param broaden_width ...
!> \param voigt_mixing ...
! **************************************************************************************************
   SUBROUTINE write_broadening_info(iw, broaden_type, broaden_width, voigt_mixing)

      INTEGER, INTENT(IN)                                :: iw, broaden_type
      REAL(KIND=dp), INTENT(IN)                          :: broaden_width, voigt_mixing

      REAL(KIND=dp)                                      :: broaden_width_ev

      broaden_width_ev = cp_unit_from_cp2k(value=broaden_width, unit_str="eV")

      SELECT CASE (broaden_type)
      CASE (broadening_gaussian)
         WRITE (UNIT=iw, FMT="(A,F10.8,A,F8.6,A)") &
            "# Gaussian broadening, FWHM = ", broaden_width, " a.u. = ", broaden_width_ev, " eV"
      CASE (broadening_lorentzian)
         WRITE (UNIT=iw, FMT="(A,F10.8,A,F8.6,A)") &
            "# Lorentzian broadening, FWHM = ", broaden_width, " a.u. = ", broaden_width_ev, " eV"
      CASE (broadening_pseudo_voigt)
         WRITE (UNIT=iw, FMT="(A,F10.8,A,F8.6,A,F8.4)") &
            "# Pseudo-Voigt broadening, FWHM = ", broaden_width, " a.u. = ", &
            broaden_width_ev, " eV, Lorentzian fraction = ", MIN(1.0_dp, MAX(0.0_dp, voigt_mixing))
      END SELECT

   END SUBROUTINE write_broadening_info

END MODULE qs_dos_utils
