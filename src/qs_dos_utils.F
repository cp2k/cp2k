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
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dos_utils'

   INTEGER, PARAMETER, PUBLIC :: broadening_gaussian = 1, &
                                 broadening_lorentzian = 2, &
                                 broadening_pseudo_voigt = 3
   INTEGER, PARAMETER, PUBLIC :: dos_energy_unit_hartree = 1, &
                                 dos_energy_unit_ev = 2
   INTEGER, PARAMETER, PUBLIC :: dos_energy_zero_auto = 1, &
                                 dos_energy_zero_absolute = 2, &
                                 dos_energy_zero_fermi = 3, &
                                 dos_energy_zero_hoco = 4

   PUBLIC :: add_broadened_peak, add_broadened_value, broadening_cutoff, broadening_function, &
             dos_density_scale, dos_energy_label, dos_energy_scale, dos_energy_zero_label, &
             dos_resolve_energy_zero, get_dos_pdos_flags, write_broadening_info

CONTAINS

! **************************************************************************************************
!> \brief Return the conversion factor from internal energy units to the selected DOS energy unit.
!> \param energy_unit ...
!> \return ...
! **************************************************************************************************
   FUNCTION dos_energy_scale(energy_unit) RESULT(scale)

      INTEGER, INTENT(IN)                                :: energy_unit
      REAL(KIND=dp)                                      :: scale

      SELECT CASE (energy_unit)
      CASE (dos_energy_unit_ev)
         scale = cp_unit_from_cp2k(value=1.0_dp, unit_str="eV")
      CASE DEFAULT
         scale = 1.0_dp
      END SELECT

   END FUNCTION dos_energy_scale

! **************************************************************************************************
!> \brief Return the DOS-density conversion factor for the selected energy unit.
!> \param energy_unit ...
!> \return ...
! **************************************************************************************************
   FUNCTION dos_density_scale(energy_unit) RESULT(scale)

      INTEGER, INTENT(IN)                                :: energy_unit
      REAL(KIND=dp)                                      :: scale

      scale = 1.0_dp/dos_energy_scale(energy_unit)

   END FUNCTION dos_density_scale

! **************************************************************************************************
!> \brief Return the energy-column label for DOS-like output.
!> \param energy_unit ...
!> \return ...
! **************************************************************************************************
   FUNCTION dos_energy_label(energy_unit) RESULT(label)

      INTEGER, INTENT(IN)                                :: energy_unit
      CHARACTER(LEN=16)                                  :: label

      SELECT CASE (energy_unit)
      CASE (dos_energy_unit_ev)
         label = "Energy[eV]"
      CASE DEFAULT
         label = "Energy[a.u.]"
      END SELECT

   END FUNCTION dos_energy_label

! **************************************************************************************************
!> \brief Return the label for the selected DOS energy zero.
!> \param energy_zero ...
!> \return ...
! **************************************************************************************************
   FUNCTION dos_energy_zero_label(energy_zero) RESULT(label)

      INTEGER, INTENT(IN)                                :: energy_zero
      CHARACTER(LEN=16)                                  :: label

      SELECT CASE (energy_zero)
      CASE (dos_energy_zero_absolute)
         label = "ABSOLUTE"
      CASE (dos_energy_zero_hoco)
         label = "HOCO"
      CASE (dos_energy_zero_auto)
         label = "AUTO"
      CASE DEFAULT
         label = "FERMI"
      END SELECT

   END FUNCTION dos_energy_zero_label

! **************************************************************************************************
!> \brief Resolve AUTO energy-zero selection for DOS-like output.
!> \param energy_zero ...
!> \param smearing_enabled ...
!> \param fractional_occupation ...
!> \return ...
! **************************************************************************************************
   FUNCTION dos_resolve_energy_zero(energy_zero, smearing_enabled, fractional_occupation) RESULT(resolved)

      INTEGER, INTENT(IN)                                :: energy_zero
      LOGICAL, INTENT(IN)                                :: smearing_enabled, fractional_occupation
      INTEGER                                            :: resolved

      IF (energy_zero == dos_energy_zero_auto) THEN
         IF (smearing_enabled .OR. fractional_occupation) THEN
            resolved = dos_energy_zero_fermi
         ELSE
            resolved = dos_energy_zero_hoco
         END IF
      ELSE
         resolved = energy_zero
      END IF

   END FUNCTION dos_resolve_energy_zero

! **************************************************************************************************
!> \brief Resolve projected-DOS requests from a DOS print section.
!> \param dos_section DOS print section
!> \param do_dos_output whether the DOS print key is active
!> \param do_projected_dos whether any projected DOS output is requested
!> \param do_pdos whether kind-resolved PDOS output is requested
!> \param do_pdos_raw whether raw kind-resolved PDOS output is requested
! **************************************************************************************************
   SUBROUTINE get_dos_pdos_flags(dos_section, do_dos_output, do_projected_dos, do_pdos, do_pdos_raw)

      TYPE(section_vals_type), POINTER                   :: dos_section
      LOGICAL, INTENT(IN)                                :: do_dos_output
      LOGICAL, INTENT(OUT)                               :: do_projected_dos, do_pdos, do_pdos_raw

      INTEGER                                            :: nrep
      LOGICAL                                            :: has_ldos, has_r_ldos
      TYPE(section_vals_type), POINTER                   :: ldos_section, pdos_section

      do_projected_dos = .FALSE.
      do_pdos = .FALSE.
      do_pdos_raw = .FALSE.
      has_ldos = .FALSE.
      has_r_ldos = .FALSE.

      IF (do_dos_output) THEN
         pdos_section => section_vals_get_subs_vals(dos_section, "PDOS")
         CALL section_vals_get(pdos_section, explicit=do_pdos)
         IF (do_pdos) CALL section_vals_val_get(pdos_section, "RAW", l_val=do_pdos_raw)
         ldos_section => section_vals_get_subs_vals(dos_section, "LDOS")
         CALL section_vals_get(ldos_section, n_repetition=nrep)
         has_ldos = (nrep > 0)
         ldos_section => section_vals_get_subs_vals(dos_section, "R_LDOS")
         CALL section_vals_get(ldos_section, n_repetition=nrep)
         has_r_ldos = (nrep > 0)
      END IF

      do_projected_dos = do_pdos .OR. has_ldos .OR. has_r_ldos

   END SUBROUTINE get_dos_pdos_flags

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
!> \brief Add a broadened spectral line with a scalar weight to a curve.
!> \param curve ...
!> \param emin ...
!> \param de ...
!> \param eig ...
!> \param weight ...
!> \param broaden_type ...
!> \param broaden_width ...
!> \param voigt_mixing ...
! **************************************************************************************************
   SUBROUTINE add_broadened_value(curve, emin, de, eig, weight, broaden_type, &
                                  broaden_width, voigt_mixing)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: curve
      REAL(KIND=dp), INTENT(IN)                          :: emin, de, eig, weight
      INTEGER, INTENT(IN)                                :: broaden_type
      REAL(KIND=dp), INTENT(IN)                          :: broaden_width, voigt_mixing

      INTEGER                                            :: i, iend, istart, nhist
      REAL(KIND=dp)                                      :: eval, line_shape

      nhist = SIZE(curve)
      istart = MAX(1, FLOOR((eig - broadening_cutoff(broaden_type, broaden_width) - emin)/de) + 1)
      iend = MIN(nhist, CEILING((eig + broadening_cutoff(broaden_type, broaden_width) - emin)/de) + 1)
      DO i = istart, iend
         eval = emin + (i - 1)*de
         line_shape = broadening_function(eval - eig, broaden_type, broaden_width, voigt_mixing)
         curve(i) = curve(i) + weight*line_shape
      END DO

   END SUBROUTINE add_broadened_value

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
