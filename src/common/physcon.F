!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition of physical constants:
!>
!>      a_bohr      : Bohr radius [m]
!>      a_fine      : Fine-structure constant
!>      a_mass      : Atomic mass unit [kg]; conversion factor [u] -> [kg]
!>      angstrom    : Conversion factor [Bohr] -> [Angstrom]
!>      bar         : Conversion factor [a.u.] -> [bar]
!>      bohr        : Conversion factor [Angstrom] -> [Bohr]
!>      boltzmann   : Boltzmann constant [J/K]
!>      c_light     : Speed of light in vacuum [m/s]
!>      c_light_au  : Speed of light in vacuum [a.u.]
!>      e_charge    : Elementary charge [C]
!>      e_mass      : Electron mass [kg]
!>      e_gfactor   : Electron g factor [ ]
!>      esu         : Conversion factors [a.u.] -> [esu]
!>      evolt       : Conversion factor [a.u.] -> [eV]
!>      femtoseconds: Conversion factor [a.u.] -> [fs]
!>      h_bar       : Planck constant [J*s]
!>      h_planck    : Planck constant [J*s]
!>      hertz       : Conversion factor [a.u.] -> [Hz]
!>      joule       : Conversion factor [a.u.] -> [J]
!>      kcalmol     : Conversion factor [a.u.] -> [kcal/mol]
!>      kelvin      : Conversion factor [a.u.] -> [K]
!>      kjmol       : Conversion factor [a.u.] -> [kJ/mol]
!>      massunit    : Conversion factor [u] -> [a.u.]
!>      mu_perm     : Magnetic constant or permeability of vacuum [N/A**2]
!>      n_avogadro  : Avogadro constant [1/mol]
!>      newton      : Conversion factor [a.u.] -> [N]
!>      pascal      : Conversion factor [a.u.] -> [Pa]
!>      permittivity: Electric constant or permittivity of vacuum [F/m]
!>      picoseconds : Conversion factor [a.u.] -> [ps]
!>      rydberg     : Rydberg constant [1/m]
!>      seconds     : Conversion factor [a.u.] -> [s]
!>      vibfac      : Conversion factor sqrt([hartree]/[u]/[Bohr**2]) -> [1/cm]
!>      wavenumbers : Conversion factor [a.u.] -> [1/cm]
!>      debye       : Conversion factor [a.u.] -> Debye
!> \note
!>      Fundamental physical constants (SI units)
!>      Literature: P. J. Mohr, D. B. Newell, B. N. Taylor, and E. Tiesinga,
!>                    "CODATA recommended values of the fundamental physical
!>                     constants: 2022", Rev. Mod. Phys. 97, 025002 (2025),
!>                    https://doi.org/10.1103/RevModPhys.97.025002
!>      This module is the authoritative source for CP2K.
!>      Independent post-processing tools maintain their own CODATA constants.
!>      Keep the Python conversions in python/src/cp2k/_units.py synchronized manually.
!> \par History
!>      - Adapted for use in CP2K (JGH)
!>      - Updated to CODATA 1998 and cleaned (05.09.2003,MK)
!>      - Updated to CODATA 2006. (26.03.2008,AK)
!>      - Updated to CODATA 2010. (10.12.2012,MK)
!>      - Turned constants into Fortran parameters (2014, Ole Schuett)
!>      - Remove all but CODATA 2006 (2015, Ole Schuett)
!>      - Updated to CODATA 2022 and the revised SI (2026)
!> \author Matthias Krack (MK)
! **************************************************************************************************
MODULE physcon

   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: a_bohr, a_fine, a_mass, angstrom, atm, bar, bohr, boltzmann, c_light_au, &
             debye, e_charge, e_gfactor, e_mass, evolt, femtoseconds, h_bar, h_planck, &
             hertz, joule, kcalmol, kelvin, kjmol, massunit, mu_perm, n_avogadro, newton, &
             p_mass, pascal, picoseconds, seconds, vibfac, wavenumbers, c_light

   PUBLIC :: codata_year, write_physcon

   INTEGER, PARAMETER :: codata_year = 2022

   ! Exact SI defining constants, including the 2019 redefinitions.
   ! Speed of light in vacuum [m/s]
   REAL(KIND=dp), PARAMETER :: c_light = 299792458.0_dp

   ! Planck constant [J*s]
   REAL(KIND=dp), PARAMETER :: h_planck = 6.62607015E-34_dp

   ! Elementary charge [C]
   REAL(KIND=dp), PARAMETER :: e_charge = 1.602176634E-19_dp

   ! Avogadro constant [1/mol]
   REAL(KIND=dp), PARAMETER :: n_avogadro = 6.02214076E+23_dp

   ! Boltzmann constant [J/K]
   REAL(KIND=dp), PARAMETER :: boltzmann = 1.380649E-23_dp

   ! Recommended CODATA table entries. Some quantities are correlated or derived in the
   ! adjustment; keep their published rounded values here, separately from CP2K conversions.

   ! Electron mass [kg]
   REAL(KIND=dp), PARAMETER :: e_mass = 9.1093837139E-31_dp

   ! Proton mass [kg]
   REAL(KIND=dp), PARAMETER :: p_mass = 1.67262192595E-27_dp

   ! Electron g factor [ ]
   REAL(KIND=dp), PARAMETER :: e_gfactor = -2.00231930436092_dp

   ! Fine-structure constant
   REAL(KIND=dp), PARAMETER :: a_fine = 7.2973525643E-3_dp

   ! Rydberg constant [1/m]
   REAL(KIND=dp), PARAMETER :: rydberg = 10973731.568157_dp

   ! Atomic mass unit [kg]; conversion factor [u] -> [kg]
   REAL(KIND=dp), PARAMETER :: a_mass = 1.66053906892E-27_dp

   ! Bohr radius [m]
   REAL(KIND=dp), PARAMETER :: a_bohr = 5.29177210544E-11_dp

   ! Derived constants and unit conversions; all numerical inputs are defined above.

   ! Speed of light in atomic units
   REAL(KIND=dp), PARAMETER :: c_light_au = 1.0_dp/a_fine

   ! Reduced Planck constant [J*s]
   REAL(KIND=dp), PARAMETER :: h_bar = h_planck/(2.0_dp*pi)

   ! Magnetic constant or permeability of vacuum [N/A**2]; no longer exact in the revised SI.
   REAL(KIND=dp), PARAMETER :: mu_perm = 2.0_dp*a_fine*h_planck/(e_charge**2*c_light)

   ! Electric constant or permittivity of vacuum [F/m]
   REAL(KIND=dp), PARAMETER :: permittivity = 1.0_dp/(mu_perm*c_light**2)

   ! [u] -> [a.u.]
   REAL(KIND=dp), PARAMETER :: massunit = a_mass/e_mass

   ! [Bohr] -> [Angstrom]
   REAL(KIND=dp), PARAMETER :: angstrom = 1.0E+10_dp*a_bohr

   ! [Angstrom] -> [Bohr]
   REAL(KIND=dp), PARAMETER :: bohr = 1.0_dp/angstrom

   ! [a.u.] -> [s]
   REAL(KIND=dp), PARAMETER :: seconds = 1.0_dp/(4.0_dp*pi*rydberg*c_light)

   ! [a.u.] -> [fs]
   REAL(KIND=dp), PARAMETER :: femtoseconds = 1.0E+15_dp*seconds

   ! [a.u.] -> [ps]
   REAL(KIND=dp), PARAMETER :: picoseconds = 1.0E+12_dp*seconds

   ! [a.u.] -> [J]
   REAL(KIND=dp), PARAMETER :: joule = 2.0_dp*rydberg*h_planck*c_light

   ! [a.u.] -> [N]
   REAL(KIND=dp), PARAMETER :: newton = joule/a_bohr

   ! [a.u.] -> [K]
   REAL(KIND=dp), PARAMETER :: kelvin = joule/boltzmann

   ! [a.u.] -> [kJ/mol]
   REAL(KIND=dp), PARAMETER :: kjmol = 0.001_dp*joule*n_avogadro

   ! [a.u.] -> [kcal/mol]
   REAL(KIND=dp), PARAMETER :: kcalmol = kjmol/4.184_dp

   ! [a.u.] -> [Pa]
   REAL(KIND=dp), PARAMETER :: pascal = joule/a_bohr**3

   ! [a.u.] -> [bar]
   REAL(KIND=dp), PARAMETER :: bar = pascal/1.0E+5_dp

   ! [a.u.] -> [atm]
   REAL(KIND=dp), PARAMETER :: atm = pascal/1.013250E+5_dp

   ! [a.u.] -> [eV]
   REAL(KIND=dp), PARAMETER :: evolt = joule/e_charge

   ! [a.u.] -> [Hz]
   REAL(KIND=dp), PARAMETER :: hertz = joule/h_planck

   ! sqrt([hartree]/[u]/[Bohr**2]) -> [1/cm] (wave numbers).
   ! Use a_mass explicitly: n_avogadro*a_mass is no longer exactly 1 g/mol in the revised SI.
   REAL(KIND=dp), PARAMETER :: vibfac = SQRT(joule/a_mass)/(200.0_dp*pi*a_bohr*c_light)

   ! [a.u.] -> [1/cm] (wave numbers)
   REAL(KIND=dp), PARAMETER :: wavenumbers = 0.02_dp*rydberg

   ! [a.u.] -> [esu] (electrostatic units)
   REAL(KIND=dp), PARAMETER :: esu_1 = 1.0E+21_dp*a_bohr*c_light*e_charge

   ! [a.u.] -> [debye] (electrostatic units)
   REAL(KIND=dp), PARAMETER :: debye = esu_1

CONTAINS

! **************************************************************************************************
!> \brief  Write all basic physical constants used by CP2K to a logical
!>           output unit.
!> \param output_unit ...
!> \date    14.11.2000
!> \par History
!>       - Updated to CODATA 1998 and cleaned (05.09.2003,MK)
!>       - Updated to CODATA 2006. (26.03.2008,AK)
!> \author  JGH
!> \version 1.1
! **************************************************************************************************
   SUBROUTINE write_physcon(output_unit)

      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=4)                                   :: year

      WRITE (year, '(I4)') codata_year

      WRITE (UNIT=output_unit, FMT="(T2,/,T2,A,/,/,(T2,A))") &
         "*** Fundamental physical constants (SI units) ***", &
         "*** Literature: P. J. Mohr, D. B. Newell, B. N. Taylor, and E. Tiesinga,", &
         "***             CODATA recommended values of the fundamental physical", &
         "***             constants: "//year, &
         "***             https://physics.nist.gov/constants"

      WRITE (UNIT=output_unit, FMT="(/,T2,A,T61,ES20.14)") &
         "Speed of light in vacuum [m/s]", c_light
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Magnetic constant or permeability of vacuum [N/A**2]", mu_perm
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Electric constant or permittivity of vacuum [F/m]", permittivity
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Planck constant (h) [J*s]", h_planck
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Planck constant (h-bar) [J*s]", h_bar
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Elementary charge [C]", e_charge
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Electron mass [kg]", e_mass
      WRITE (UNIT=output_unit, FMT="(T2,A,T60,ES21.14)") &
         "Electron g factor [ ]", e_gfactor
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Proton mass [kg]", p_mass
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Fine-structure constant", a_fine
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Rydberg constant [1/m]", rydberg
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Avogadro constant [1/mol]", n_avogadro
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Boltzmann constant [J/K]", boltzmann
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Atomic mass unit [kg]", a_mass
      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "Bohr radius [m]", a_bohr

      ! Conversion factors

      WRITE (UNIT=output_unit, FMT="(/,T2,A,/)") &
         "*** Conversion factors ***"

      WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES20.14)") &
         "[u] -> [a.u.]", massunit, &
         "[Angstrom] -> [Bohr] = [a.u.]", bohr, &
         "[a.u.] = [Bohr] -> [Angstrom]", angstrom, &
         "[a.u.] -> [s]", seconds, &
         "[a.u.] -> [fs]", femtoseconds, &
         "[a.u.] -> [J]", joule, &
         "[a.u.] -> [N]", newton, &
         "[a.u.] -> [K]", kelvin, &
         "[a.u.] -> [kJ/mol]", kjmol, &
         "[a.u.] -> [kcal/mol]", kcalmol, &
         "[a.u.] -> [Pa]", pascal, &
         "[a.u.] -> [bar]", bar, &
         "[a.u.] -> [atm]", atm, &
         "[a.u.] -> [eV]", evolt, &
         "[a.u.] -> [Hz]", hertz, &
         "[a.u.] -> [1/cm] (wave numbers)", wavenumbers, &
         "sqrt([hartree]/[u]/[Bohr**2]) -> [1/cm]", vibfac
      WRITE (UNIT=output_unit, FMT="(T2,A)") ""

   END SUBROUTINE write_physcon

END MODULE physcon
