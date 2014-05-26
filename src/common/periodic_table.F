!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Periodic Table related data definitions
!> \par Literature
!>      M. E. Wieser and M. Berglund,
!>      Atomic weights of the elements 2007 (IUPAC Technical Report),
!>      Pure Appl. Chem. 81, 2131-2156 (2009)
!> \par History
!>      27.05.2010 Update of the atomic weights
!>      20.10.2009 Use vdW Radii from the Cambridge Structural Database
!> \author JGH, MK
! *****************************************************************************
MODULE periodic_table

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_periodic_table, ptable, nelem,&
            get_ptable_info

  TYPE atom
    CHARACTER(LEN=2)  :: symbol
    CHARACTER(LEN=14) :: name
    INTEGER           :: number
    REAL(KIND=dp)     :: amass ! Standard atomic weight, i.e. the relative atomic mass
                               ! of an element derived by averaging the relative
                               ! atomic masses of the (natural) isotopes of that
                               ! element or the relative (possibly estimated)
                               ! longest half-life period in formula units.
    REAL(KIND=dp)     :: covalent_radius ! in Angstroms
    REAL(KIND=dp)     :: vdw_radius ! in Angstroms
    INTEGER           :: e_conv(0:3)
    REAL(KIND=dp)     :: heat_of_formation ! in kcal/mol
    REAL(KIND=dp)     :: eht_param (0:3) ! in eV
    REAL(KIND=dp)     :: gyrom_ratio ! in MHz/Tesla
    INTEGER           :: gyrom_ratio_isotope ! isotope number corresponding to gyrom_ratio
  END TYPE atom

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'periodic_table'
  REAL(KIND=dp), PARAMETER             :: z = 0.0_dp
  INTEGER, PARAMETER                   :: nelem = 118

  TYPE(atom) :: ptable(0:nelem)

CONTAINS


! *****************************************************************************
!> \brief Pass information about the kind given the element symbol
!> \par History
!>      none
!> \author MI
! *****************************************************************************
  SUBROUTINE get_ptable_info(symbol,number,amass,ielement,covalent_radius,&
                             vdw_radius,found)

    CHARACTER(LEN=2), INTENT(IN)             :: symbol
    INTEGER, INTENT(OUT), OPTIONAL           :: number
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: amass
    INTEGER, INTENT(OUT), OPTIONAL           :: ielement
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: covalent_radius, vdw_radius
    LOGICAL, OPTIONAL                        :: found

    CHARACTER(len=*), PARAMETER :: routineN = 'get_ptable_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: symb_ielem, symbol_in
    INTEGER                                  :: ielem, ielem_found

    ielem_found = -1
    symbol_in = symbol
    CALL uppercase(symbol_in)
    DO ielem=1,nelem
      symb_ielem = ptable(ielem)%symbol
      CALL uppercase(symb_ielem)
      IF (symbol_in == symb_ielem) THEN
        ielem_found = ielem
        EXIT
      END IF
    END DO

    ! If requested give back the value and avoid failing
    IF (PRESENT(found)) THEN
      found = (ielem_found > 0)
    ELSE
      IF (ielem_found < 0) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Unknown element symbol <"//TRIM(symbol)//"> found.")
      END IF
    END IF

    ! If the element has been identified provide the requested information
    IF (ielem_found > 0) THEN
      IF(PRESENT(ielement)) ielement = ielem_found
      IF(PRESENT(number)) number = ptable(ielem_found)%number
      IF(PRESENT(amass)) amass = ptable(ielem_found)%amass
      IF(PRESENT(covalent_radius)) covalent_radius = ptable(ielem_found)%covalent_radius
      IF(PRESENT(vdw_radius)) vdw_radius = ptable(ielem_found)%vdw_radius
    END IF

  END SUBROUTINE get_ptable_info

! *****************************************************************************
!> \brief Initialization of Periodic Table related data
!> \par History
!>      27.05.2010 Update of the atomic weights
!> \author JGH, MK
! *****************************************************************************
  SUBROUTINE init_periodic_table()

    ! Dummy
    ptable(0) % symbol = 'X '
    ptable(0) % name = 'Dummy'
    ptable(0) % number = 0
    ptable(0) % amass = z
    ptable(0) % covalent_radius = z
    ptable(0) % vdw_radius = z
    ptable(0) % e_conv(0:3) = (/ 0, 0, 0, 0 /)
    ptable(0) % eht_param(0:3) = (/ z, z, z, z/)

    ! Hydrogen
    ptable(1) % symbol = 'H '
    ptable(1) % name = 'Hydrogen'
    ptable(1) % number = 1
    ptable(1) % amass = 1.00794_dp
    ptable(1) % covalent_radius = 0.32_dp
    ptable(1) % vdw_radius = 1.09_dp
    ptable(1) % e_conv(0:3) = (/ 1, 0, 0, 0 /)
    ptable(1) % eht_param(0:3) = (/ -13.60_dp, z, z, z/)

    ! Helium
    ptable(2) % symbol = 'He'
    ptable(2) % name = 'Helium'
    ptable(2) % number = 2
    ptable(2) % amass = 4.002602_dp
    ptable(2) % covalent_radius = 0.9300_dp
    ptable(2) % vdw_radius = 1.40_dp
    ptable(2) % e_conv(0:3) = (/ 2, 0, 0, 0 /)
    ptable(2) % eht_param(0:3) = (/ -23.40_dp, z, z, z/)

    ! Lithium
    ptable(3) % symbol = 'Li'
    ptable(3) % name = 'Lithium'
    ptable(3) % number = 3
    ptable(3) % amass = 6.941_dp
    ptable(3) % covalent_radius = 1.2300_dp
    ptable(3) % vdw_radius = 1.82_dp
    ptable(3) % e_conv(0:3) = (/ 3, 0, 0, 0 /)
    ptable(3) % eht_param(0:3) = (/ -5.40_dp, -3.50_dp, z, z/)

    ! Beryllium
    ptable(4) % symbol = 'Be'
    ptable(4) % name = 'Beryllium'
    ptable(4) % number = 4
    ptable(4) % amass = 9.012182_dp
    ptable(4) % covalent_radius = 0.9000_dp
    ptable(4) % vdw_radius = 2.00_dp
    ptable(4) % e_conv(0:3) = (/ 4, 0, 0, 0 /)
    ptable(4) % eht_param(0:3) = (/ -10.00_dp, -6.00_dp, z, z/)

    ! Boron
    ptable(5) % symbol = 'B '
    ptable(5) % name = 'Boron'
    ptable(5) % number = 5
    ptable(5) % amass = 10.811_dp
    ptable(5) % covalent_radius = 0.8200_dp
    ptable(5) % vdw_radius = 2.00_dp
    ptable(5) % e_conv(0:3) = (/ 4, 1, 0, 0 /)
    ptable(5) % eht_param(0:3) = (/ -15.20_dp, -8.50_dp, z, z/)

    ! Carbon
    ptable(6) % symbol = 'C '
    ptable(6) % name = 'Carbon'
    ptable(6) % number = 6
    ptable(6) % amass = 12.0107_dp
    ptable(6) % covalent_radius = 0.7700_dp
    ptable(6) % vdw_radius = 1.7_dp
    ptable(6) % e_conv(0:3) = (/ 4, 2, 0, 0 /)
    ptable(6) % eht_param(0:3) = (/ -21.40_dp, -11.40_dp, z, z/)

    ! Nitrogen
    ptable(7) % symbol = 'N '
    ptable(7) % name = 'Nitrogen'
    ptable(7) % number = 7
    ptable(7) % amass = 14.0067_dp
    ptable(7) % covalent_radius = 0.7500_dp
    ptable(7) % vdw_radius = 1.55_dp
    ptable(7) % e_conv(0:3) = (/ 4, 3, 0, 0 /)
    ptable(7) % eht_param(0:3) = (/ -26.00_dp, -13.40_dp, z, z/)

    ! Oxygen
    ptable(8) % symbol = 'O '
    ptable(8) % name = 'Oxygen'
    ptable(8) % number = 8
    ptable(8) % amass = 15.9994_dp
    ptable(8) % covalent_radius = 0.7300_dp
    ptable(8) % vdw_radius = 1.52_dp
    ptable(8) % e_conv(0:3) = (/ 4, 4, 0, 0 /)
    ptable(8) % eht_param(0:3) = (/ -32.30_dp, -14.80_dp, z, z/)

    ! Fluorine
    ptable(9) % symbol = 'F '
    ptable(9) % name = 'Fluorine'
    ptable(9) % number = 9
    ptable(9) % amass = 18.9984032_dp
    ptable(9) % covalent_radius = 0.7200_dp
    ptable(9) % vdw_radius = 1.47_dp
    ptable(9) % e_conv(0:3) = (/ 4, 5, 0, 0 /)
    ptable(9) % eht_param(0:3) = (/ -40.00_dp, -18.10_dp, z, z/)

    ! Neon
    ptable(10) % symbol = 'Ne'
    ptable(10) % name = 'Neon'
    ptable(10) % number = 10
    ptable(10) % amass = 20.1797_dp
    ptable(10) % covalent_radius = 0.7100_dp
    ptable(10) % vdw_radius = 1.54_dp
    ptable(10) % e_conv(0:3) = (/ 4, 6, 0, 0 /)
    ptable(10) % eht_param(0:3) = (/ -43.20_dp, -20.00_dp, z, z/)

    ! Sodium
    ptable(11) % symbol = 'Na'
    ptable(11) % name = 'Sodium'
    ptable(11) % number = 11
    ptable(11) % amass = 22.98976928_dp
    ptable(11) % covalent_radius = 1.5400_dp
    ptable(11) % vdw_radius = 2.27_dp
    ptable(11) % e_conv(0:3) = (/ 5, 6, 0, 0 /)
    ptable(11) % eht_param(0:3) = (/ -5.10_dp, -3.00_dp, z, z/)

    ! Magnesium
    ptable(12) % symbol = 'Mg'
    ptable(12) % name = 'Magnesium'
    ptable(12) % number = 12
    ptable(12) % amass = 24.305_dp
    ptable(12) % covalent_radius = 1.3600_dp
    ptable(12) % vdw_radius = 1.73_dp
    ptable(12) % e_conv(0:3) = (/ 6, 6, 0, 0 /)
    ptable(12) % eht_param(0:3) = (/ -9.00_dp, -4.50_dp, z, z/)

    ! Aluminium
    ptable(13) % symbol = 'Al'
    ptable(13) % name = 'Aluminium'
    ptable(13) % number = 13
    ptable(13) % amass = 26.9815386_dp
    ptable(13) % covalent_radius = 1.1800_dp
    ptable(13) % vdw_radius = 2.00_dp
    ptable(13) % e_conv(0:3) = (/ 6, 7, 0, 0 /)
    ptable(13) % eht_param(0:3) = (/ -12.30_dp, -6.50_dp, z, z/)

    ! Silicon
    ptable(14) % symbol = 'Si'
    ptable(14) % name = 'Silicon'
    ptable(14) % number = 14
    ptable(14) % amass = 28.0855_dp
    ptable(14) % covalent_radius = 1.1100_dp
    ptable(14) % vdw_radius = 2.10_dp
    ptable(14) % e_conv(0:3) = (/ 6, 8, 0, 0 /)
    ptable(14) % eht_param(0:3) = (/ -17.30_dp, -9.20_dp, z, z/)

    ! Phosphorus
    ptable(15) % symbol = 'P '
    ptable(15) % name = 'Phosphorus'
    ptable(15) % number = 15
    ptable(15) % amass = 30.973762_dp
    ptable(15) % covalent_radius = 1.0600_dp
    ptable(15) % vdw_radius = 1.80_dp
    ptable(15) % e_conv(0:3) = (/ 6, 9, 0, 0 /)
    ptable(15) % eht_param(0:3) = (/ -18.60_dp, -14.00_dp, z, z/)

    ! Sulfur
    ptable(16) % symbol = 'S '
    ptable(16) % name = 'Sulfur'
    ptable(16) % number = 16
    ptable(16) % amass = 32.065_dp
    ptable(16) % covalent_radius = 1.0200_dp
    ptable(16) % vdw_radius = 1.80_dp
    ptable(16) % e_conv(0:3) = (/ 6, 10, 0, 0 /)
    ptable(16) % eht_param(0:3) = (/ -20.00_dp, -11.00_dp, z, z/)

    ! Chlorine
    ptable(17) % symbol = 'Cl'
    ptable(17) % name = 'Chlorine'
    ptable(17) % number = 17
    ptable(17) % amass = 35.453_dp
    ptable(17) % covalent_radius = 0.9900_dp
    ptable(17) % vdw_radius = 1.75_dp
    ptable(17) % e_conv(0:3) = (/ 6, 11, 0, 0 /)
    ptable(17) % eht_param(0:3) = (/ -26.30_dp, -14.20_dp, z, z/)

    ! Argon
    ptable(18) % symbol = 'Ar'
    ptable(18) % name = 'Argon'
    ptable(18) % number = 18
    ptable(18) % amass = 39.948_dp
    ptable(18) % covalent_radius = 0.9800_dp
    ptable(18) % vdw_radius = 1.88_dp
    ptable(18) % e_conv(0:3) = (/ 6, 12, 0, 0 /)
    ptable(18) % eht_param(0:3) = (/ z, z, z, z/)

    ! Potassium
    ptable(19) % symbol = 'K '
    ptable(19) % name = 'Potassium'
    ptable(19) % number = 19
    ptable(19) % amass = 39.0983_dp
    ptable(19) % covalent_radius = 2.0300_dp
    ptable(19) % vdw_radius = 2.75_dp
    ptable(19) % e_conv(0:3) = (/ 7, 12, 0, 0 /)
    ptable(19) % eht_param(0:3) = (/ -4.34_dp, -2.73_dp, z, z/)

    ! Calcium
    ptable(20) % symbol = 'Ca'
    ptable(20) % name = 'Calcium'
    ptable(20) % number = 20
    ptable(20) % amass = 40.078_dp
    ptable(20) % covalent_radius = 1.7400_dp
    ptable(20) % vdw_radius = 2.00_dp
    ptable(20) % e_conv(0:3) = (/ 8, 12, 0, 0 /)
    ptable(20) % eht_param(0:3) = (/ -7.00_dp, -4.00_dp, z, z/)

    ! Scandium
    ptable(21) % symbol = 'Sc'
    ptable(21) % name = 'Scandium'
    ptable(21) % number = 21
    ptable(21) % amass = 44.955912_dp
    ptable(21) % covalent_radius = 1.4400_dp
    ptable(21) % vdw_radius = 2.00_dp
    ptable(21) % e_conv(0:3) = (/ 8, 12, 1, 0 /)
    ptable(21) % eht_param(0:3) = (/ -8.87_dp, -2.75_dp, -8.51_dp, z/)

    ! Titanium
    ptable(22) % symbol = 'Ti'
    ptable(22) % name = 'Titanium'
    ptable(22) % number = 22
    ptable(22) % amass = 47.867_dp
    ptable(22) % covalent_radius = 1.3200_dp
    ptable(22) % vdw_radius = 2.00_dp
    ptable(22) % e_conv(0:3) = (/ 8, 12, 2, 0 /)
    ptable(22) % eht_param(0:3) = (/ -8.97_dp, -5.44_dp, -10.81_dp, z/)

    ! Vanadium
    ptable(23) % symbol = 'V '
    ptable(23) % name = 'Vanadium'
    ptable(23) % number = 23
    ptable(23) % amass = 50.9415_dp
    ptable(23) % covalent_radius = 1.2200_dp
    ptable(23) % vdw_radius = 2.00_dp
    ptable(23) % e_conv(0:3) = (/ 8, 12, 3, 0 /)
    ptable(23) % eht_param(0:3) = (/ -8.81_dp, -5.52_dp, -11.00_dp, z/)

    ! Chromium
    ptable(24) % symbol = 'Cr'
    ptable(24) % name = 'Chromium'
    ptable(24) % number = 24
    ptable(24) % amass = 51.9961_dp
    ptable(24) % covalent_radius = 1.1800_dp
    ptable(24) % vdw_radius = 2.00_dp
    ptable(24) % e_conv(0:3) = (/ 7, 12, 5, 0 /)
    ptable(24) % eht_param(0:3) = (/ -8.66_dp, -5.24_dp, -11.22_dp, z/)

    ! Manganese
    ptable(25) % symbol = 'Mn'
    ptable(25) % name = 'Manganese'
    ptable(25) % number = 25
    ptable(25) % amass = 54.938045_dp
    ptable(25) % covalent_radius = 1.1700_dp
    ptable(25) % vdw_radius = 2.00_dp
    ptable(25) % e_conv(0:3) = (/ 8, 12, 5, 0 /)
    ptable(25) % eht_param(0:3) = (/ -9.75_dp, -5.89_dp, -11.67_dp, z/)

    ! Iron
    ptable(26) % symbol = 'Fe'
    ptable(26) % name = 'Iron'
    ptable(26) % number = 26
    ptable(26) % amass = 55.845_dp
    ptable(26) % covalent_radius = 1.1700_dp
    ptable(26) % vdw_radius = 2.00_dp
    ptable(26) % e_conv(0:3) = (/ 8, 12, 6, 0 /)
    ptable(26) % eht_param(0:3) = (/ -9.10_dp, -5.32_dp, -12.60_dp, z/)

    ! Cobalt
    ptable(27) % symbol = 'Co'
    ptable(27) % name = 'Cobalt'
    ptable(27) % number = 27
    ptable(27) % amass = 58.933195_dp
    ptable(27) % covalent_radius = 1.1600_dp
    ptable(27) % vdw_radius = 2.00_dp
    ptable(27) % e_conv(0:3) = (/ 8, 12, 7, 0 /)
    ptable(27) % eht_param(0:3) = (/ -9.21_dp, -5.29_dp, -13.18_dp, z/)

    ! Nickel
    ptable(28) % symbol = 'Ni'
    ptable(28) % name = 'Nickel'
    ptable(28) % number = 28
    ptable(28) % amass = 58.6934_dp
    ptable(28) % covalent_radius = 1.1500_dp
    ptable(28) % vdw_radius = 1.63_dp
    ptable(28) % e_conv(0:3) = (/ 8, 12, 8, 0 /)
    ptable(28) % eht_param(0:3) = (/ -9.17_dp, -5.15_dp, -13.49_dp, z/)

    ! Copper
    ptable(29) % symbol = 'Cu'
    ptable(29) % name = 'Copper'
    ptable(29) % number = 29
    ptable(29) % amass = 63.546_dp
    ptable(29) % covalent_radius = 1.1700_dp
    ptable(29) % vdw_radius = 1.40_dp
    ptable(29) % e_conv(0:3) = (/ 7, 12, 10, 0 /)
    ptable(29) % eht_param(0:3) = (/ -11.40_dp, -6.06_dp, -14.00_dp, z/)

    ! Zinc
    ptable(30) % symbol = 'Zn'
    ptable(30) % name = 'Zinc'
    ptable(30) % number = 30
    ptable(30) % amass = 65.38_dp
    ptable(30) % covalent_radius = 1.2500_dp
    ptable(30) % vdw_radius = 1.39_dp
    ptable(30) % e_conv(0:3) = (/ 8, 12, 10, 0 /)
    ptable(30) % eht_param(0:3) = (/ -12.41_dp, -6.53_dp, z, z/)

    ! Gallium
    ptable(31) % symbol = 'Ga'
    ptable(31) % name = 'Gallium'
    ptable(31) % number = 31
    ptable(31) % amass = 69.723_dp
    ptable(31) % covalent_radius = 1.2600_dp
    ptable(31) % vdw_radius = 1.87_dp
    ptable(31) % e_conv(0:3) = (/ 8, 13, 10, 0 /)
    ptable(31) % eht_param(0:3) = (/ -14.58_dp, -6.75_dp, z, z/)

    ! Germanium
    ptable(32) % symbol = 'Ge'
    ptable(32) % name = 'Germanium'
    ptable(32) % number = 32
    ptable(32) % amass = 72.64_dp
    ptable(32) % covalent_radius = 1.2200_dp
    ptable(32) % vdw_radius = 2.00_dp
    ptable(32) % e_conv(0:3) = (/ 8, 14, 10, 0 /)
    ptable(32) % eht_param(0:3) = (/ -16.00_dp, -9.00_dp, z, z/)

    ! Arsenic
    ptable(33) % symbol = 'As'
    ptable(33) % name = 'Arsenic'
    ptable(33) % number = 33
    ptable(33) % amass = 74.9216_dp
    ptable(33) % covalent_radius = 1.2000_dp
    ptable(33) % vdw_radius = 1.85_dp
    ptable(33) % e_conv(0:3) = (/ 8, 15, 10, 0 /)
    ptable(33) % eht_param(0:3) = (/ -16.22_dp, -12.16_dp, z, z/)

    ! Selenium
    ptable(34) % symbol = 'Se'
    ptable(34) % name = 'Selenium'
    ptable(34) % number = 34
    ptable(34) % amass = 78.96_dp
    ptable(34) % covalent_radius = 1.1600_dp
    ptable(34) % vdw_radius = 1.90_dp
    ptable(34) % e_conv(0:3) = (/ 8, 16, 10, 0 /)
    ptable(34) % eht_param(0:3) = (/ -20.50_dp, -14.40_dp, z, z/)

    ! Bromine
    ptable(35) % symbol = 'Br'
    ptable(35) % name = 'Bromine'
    ptable(35) % number = 35
    ptable(35) % amass = 79.904_dp
    ptable(35) % covalent_radius = 1.1400_dp
    ptable(35) % vdw_radius = 1.85_dp
    ptable(35) % e_conv(0:3) = (/ 8, 17, 10, 0 /)
    ptable(35) % eht_param(0:3) = (/ -22.07_dp, -13.10_dp, z, z/)

    ! Krypton
    ptable(36) % symbol = 'Kr'
    ptable(36) % name = 'Krypton'
    ptable(36) % number = 36
    ptable(36) % amass = 83.798_dp
    ptable(36) % covalent_radius = 1.1200_dp
    ptable(36) % vdw_radius = 2.02_dp
    ptable(36) % e_conv(0:3) = (/ 8, 18, 10, 0 /)
    ptable(36) % eht_param(0:3) = (/ z, z, z, z/)

    ! Rubidium
    ptable(37) % symbol = 'Rb'
    ptable(37) % name = 'Rubidium'
    ptable(37) % number = 37
    ptable(37) % amass = 85.4678_dp
    ptable(37) % covalent_radius = 2.1600_dp
    ptable(37) % vdw_radius = 2.00_dp
    ptable(37) % e_conv(0:3) = (/ 9, 18, 10, 0 /)
    ptable(37) % eht_param(0:3) = (/ -4.18_dp, -2.60_dp, z, z/)

    ! Strontium
    ptable(38) % symbol = 'Sr'
    ptable(38) % name = 'Strontium'
    ptable(38) % number = 38
    ptable(38) % amass = 87.62_dp
    ptable(38) % covalent_radius = 1.9100_dp
    ptable(38) % vdw_radius = 2.00_dp
    ptable(38) % e_conv(0:3) = (/ 10, 18, 10, 0 /)
    ptable(38) % eht_param(0:3) = (/ -6.62_dp, -3.92_dp, z, z/)

    ! Yttrium
    ptable(39) % symbol = 'Y '
    ptable(39) % name = 'Yttrium'
    ptable(39) % number = 39
    ptable(39) % amass = 88.90585_dp
    ptable(39) % covalent_radius = 1.6200_dp
    ptable(39) % vdw_radius = 2.00_dp
    ptable(39) % e_conv(0:3) = (/ 10, 18, 11, 0 /)
    ptable(39) % eht_param(0:3) = (/ z, z, z, z/)

    ! Zirconium
    ptable(40) % symbol = 'Zr'
    ptable(40) % name = 'Zirconium'
    ptable(40) % number = 40
    ptable(40) % amass = 91.224_dp
    ptable(40) % covalent_radius = 1.4500_dp
    ptable(40) % vdw_radius = 2.00_dp
    ptable(40) % e_conv(0:3) = (/ 10, 18, 12, 0 /)
    ptable(40) % eht_param(0:3) = (/ -8.00_dp, -5.40_dp, -10.20_dp, z/)

    ! Niobium
    ptable(41) % symbol = 'Nb'
    ptable(41) % name = 'Niobium'
    ptable(41) % number = 41
    ptable(41) % amass = 92.90638_dp
    ptable(41) % covalent_radius = 1.3400_dp
    ptable(41) % vdw_radius = 2.00_dp
    ptable(41) % e_conv(0:3) = (/ 9, 18, 14, 0 /)
    ptable(41) % eht_param(0:3) = (/ -10.10_dp, -6.86_dp, -12.10_dp, z/)

    ! Molybdenum
    ptable(42) % symbol = 'Mo'
    ptable(42) % name = 'Molybdenum'
    ptable(42) % number = 42
    ptable(42) % amass = 95.96_dp
    ptable(42) % covalent_radius = 1.3000_dp
    ptable(42) % vdw_radius = 2.00_dp
    ptable(42) % e_conv(0:3) = (/ 9, 18, 15, 0 /)
    ptable(42) % eht_param(0:3) = (/ -8.34_dp, -5.25_dp, -10.50_dp, z/)

    ! Technetium
    ptable(43) % symbol = 'Tc'
    ptable(43) % name = 'Technetium'
    ptable(43) % number = 43
    ptable(43) % amass = 97.9072_dp ! (98)
    ptable(43) % covalent_radius = 1.2700_dp
    ptable(43) % vdw_radius = 2.00_dp
    ptable(43) % e_conv(0:3) = (/10, 18, 15, 0 /)
    ptable(43) % eht_param(0:3) = (/ -10.07_dp, -5.40_dp, -12.82_dp, z/)

    ! Ruthenium
    ptable(44) % symbol = 'Ru'
    ptable(44) % name = 'Ruthenium'
    ptable(44) % number = 44
    ptable(44) % amass = 101.07_dp
    ptable(44) % covalent_radius = 1.2500_dp
    ptable(44) % vdw_radius = 2.00_dp
    ptable(44) % e_conv(0:3) = (/ 9, 18, 17, 0 /)
    ptable(44) % eht_param(0:3) = (/ -10.40_dp, -6.87_dp, -14.90_dp, z/)

    ! Rhodium
    ptable(45) % symbol = 'Rh'
    ptable(45) % name = 'Rhodium'
    ptable(45) % number = 45
    ptable(45) % amass = 102.9055_dp
    ptable(45) % covalent_radius = 1.2500_dp
    ptable(45) % vdw_radius = 2.00_dp
    ptable(45) % e_conv(0:3) = (/ 9, 18, 18, 0 /)
    ptable(45) % eht_param(0:3) = (/ -8.09_dp, -4.57_dp, -12.50_dp, z/)

    ! Palladium
    ptable(46) % symbol = 'Pd'
    ptable(46) % name = 'Palladium'
    ptable(46) % number = 46
    ptable(46) % amass = 106.42_dp
    ptable(46) % covalent_radius = 1.2800_dp
    ptable(46) % vdw_radius = 1.63_dp
    ptable(46) % e_conv(0:3) = (/ 8, 18, 20, 0 /)
    ptable(46) % eht_param(0:3) = (/ -7.32_dp, -3.75_dp, -12.02_dp, z/)

    ! Silver
    ptable(47) % symbol = 'Ag'
    ptable(47) % name = 'Silver'
    ptable(47) % number = 47
    ptable(47) % amass = 107.8682_dp
    ptable(47) % covalent_radius = 1.3400_dp
    ptable(47) % vdw_radius = 1.72_dp
    ptable(47) % e_conv(0:3) = (/ 9, 18, 20, 0 /)
    ptable(47) % eht_param(0:3) = (/ z, z, z, z/)

    ! Cadmium
    ptable(48) % symbol = 'Cd'
    ptable(48) % name = 'Cadmium'
    ptable(48) % number = 48
    ptable(48) % amass = 112.411_dp
    ptable(48) % covalent_radius = 1.4800_dp
    ptable(48) % vdw_radius = 1.58_dp
    ptable(48) % e_conv(0:3) = (/ 10, 18, 20, 0 /)
    ptable(48) % eht_param(0:3) = (/ z, z, z, z/)

    ! Indium
    ptable(49) % symbol = 'In'
    ptable(49) % name = 'Indium'
    ptable(49) % number = 49
    ptable(49) % amass = 114.818_dp
    ptable(49) % covalent_radius = 1.4400_dp
    ptable(49) % vdw_radius = 1.93_dp
    ptable(49) % e_conv(0:3) = (/ 10, 19, 20, 0 /)
    ptable(49) % eht_param(0:3) = (/ -12.60_dp, -6.19_dp, z, z/)

    ! Tin
    ptable(50) % symbol = 'Sn'
    ptable(50) % name = 'Tin'
    ptable(50) % number = 50
    ptable(50) % amass = 118.71_dp
    ptable(50) % covalent_radius = 1.4100_dp
    ptable(50) % vdw_radius = 2.17_dp
    ptable(50) % e_conv(0:3) = (/ 10, 20, 20, 0 /)
    ptable(50) % eht_param(0:3) = (/ -16.16_dp, -8.32_dp, z, z/)

    ! Antimony
    ptable(51) % symbol = 'Sb'
    ptable(51) % name = 'Antimony'
    ptable(51) % number = 51
    ptable(51) % amass = 121.76_dp
    ptable(51) % covalent_radius = 1.4000_dp
    ptable(51) % vdw_radius = 2.00_dp
    ptable(51) % e_conv(0:3) = (/ 10, 21, 20, 0 /)
    ptable(51) % eht_param(0:3) = (/ -18.80_dp, -11.70_dp, z, z/)

    ! Tellurium
    ptable(52) % symbol = 'Te'
    ptable(52) % name = 'Tellurium'
    ptable(52) % number = 52
    ptable(52) % amass = 127.6_dp
    ptable(52) % covalent_radius = 1.3600_dp
    ptable(52) % vdw_radius = 2.06_dp
    ptable(52) % e_conv(0:3) = (/ 10, 22, 20, 0 /)
    ptable(52) % eht_param(0:3) = (/ -20.80_dp, -13.20_dp, z, z/)

    ! Iodine
    ptable(53) % symbol = 'I '
    ptable(53) % name = 'Iodine'
    ptable(53) % number = 53
    ptable(53) % amass = 126.90447_dp
    ptable(53) % covalent_radius = 1.3300_dp
    ptable(53) % vdw_radius = 1.98_dp
    ptable(53) % e_conv(0:3) = (/ 10, 23, 20, 0 /)
    ptable(53) % eht_param(0:3) = (/ -18.00_dp, -12.70_dp, z, z/)

    ! Xenon
    ptable(54) % symbol = 'Xe'
    ptable(54) % name = 'Xenon'
    ptable(54) % number = 54
    ptable(54) % amass = 131.293_dp
    ptable(54) % covalent_radius = 1.3100_dp
    ptable(54) % vdw_radius = 2.16_dp
    ptable(54) % e_conv(0:3) = (/ 10, 24, 20, 0 /)
    ptable(54) % eht_param(0:3) = (/ z, z, z, z/)

    ! Cesium
    ptable(55) % symbol = 'Cs'
    ptable(55) % name = 'Cesium'
    ptable(55) % number = 55
    ptable(55) % amass = 132.9054519_dp
    ptable(55) % covalent_radius = 2.3500_dp
    ptable(55) % vdw_radius = 2.00_dp
    ptable(55) % e_conv(0:3) = (/ 11, 24, 20, 0 /)
    ptable(55) % eht_param(0:3) = (/ -3.88_dp, -2.49_dp, z, z/)

    ! Barium
    ptable(56) % symbol = 'Ba'
    ptable(56) % name = 'Barium'
    ptable(56) % number = 56
    ptable(56) % amass = 137.327_dp
    ptable(56) % covalent_radius = 1.9800_dp
    ptable(56) % vdw_radius = 2.00_dp
    ptable(56) % e_conv(0:3) = (/ 12, 24, 20, 0 /)
    ptable(56) % eht_param(0:3) = (/ z, z, z, z/)

    ! Lanthanum
    ptable(57) % symbol = 'La'
    ptable(57) % name = 'Lanthanum'
    ptable(57) % number = 57
    ptable(57) % amass = 138.90547_dp
    ptable(57) % covalent_radius = 1.6900_dp
    ptable(57) % vdw_radius = 2.00_dp
    ptable(57) % e_conv(0:3) = (/ 12, 24, 21, 0 /)
    ptable(57) % eht_param(0:3) = (/ -7.67_dp, -5.01_dp, -8.21_dp, z/)

    ! Cerium
    ptable(58) % symbol = 'Ce'
    ptable(58) % name = 'Cerium'
    ptable(58) % number = 58
    ptable(58) % amass = 140.116_dp
    ptable(58) % covalent_radius = 1.6500_dp
    ptable(58) % vdw_radius = 2.00_dp
    ptable(58) % e_conv(0:3) = (/ 12, 24, 21, 1 /)
    ptable(58) % eht_param(0:3) = (/ z, z, z, z/)

    ! Praseodymium
    ptable(59) % symbol = 'Pr'
    ptable(59) % name = 'Praseodymium'
    ptable(59) % number = 59
    ptable(59) % amass = 140.90765_dp
    ptable(59) % covalent_radius = 1.6500_dp
    ptable(59) % vdw_radius = 2.00_dp
    ptable(59) % e_conv(0:3) = (/ 12, 24, 20, 3 /)
    ptable(59) % eht_param(0:3) = (/ z, z, z, z/)

    ! Neodymium
    ptable(60) % symbol = 'Nd'
    ptable(60) % name = 'Neodymium'
    ptable(60) % number = 60
    ptable(60) % amass = 144.242_dp
    ptable(60) % covalent_radius = 1.6400_dp
    ptable(60) % vdw_radius = 2.00_dp
    ptable(60) % e_conv(0:3) = (/ 12, 24, 20, 4 /)
    ptable(60) % eht_param(0:3) = (/ z, z, z, z/)

    ! Promethium
    ptable(61) % symbol = 'Pm'
    ptable(61) % name = 'Promethium'
    ptable(61) % number = 61
    ptable(61) % amass = 144.9127_dp ! (145)
    ptable(61) % covalent_radius = 1.6300_dp
    ptable(61) % vdw_radius = 2.00_dp
    ptable(61) % e_conv(0:3) = (/ 12, 24, 20, 5 /)
    ptable(61) % eht_param(0:3) = (/ z, z, z, z/)

    ! Samarium
    ptable(62) % symbol = 'Sm'
    ptable(62) % name = 'Samarium'
    ptable(62) % number = 62
    ptable(62) % amass = 150.36_dp
    ptable(62) % covalent_radius = 1.6200_dp
    ptable(62) % vdw_radius = 2.00_dp
    ptable(62) % e_conv(0:3) = (/ 12, 24, 20, 6 /)
    ptable(62) % eht_param(0:3) = (/ -4.86_dp, -4.86_dp, -6.06_dp, -11.28_dp/)

    ! Europium
    ptable(63) % symbol = 'Eu'
    ptable(63) % name = 'Europium'
    ptable(63) % number = 63
    ptable(63) % amass = 151.964_dp
    ptable(63) % covalent_radius = 1.8500_dp
    ptable(63) % vdw_radius = 2.00_dp
    ptable(63) % e_conv(0:3) = (/ 12, 24, 20, 7 /)
    ptable(63) % eht_param(0:3) = (/ z, z, z, z/)

    ! Gadolinium
    ptable(64) % symbol = 'Gd'
    ptable(64) % name = 'Gadolinium'
    ptable(64) % number = 64
    ptable(64) % amass = 157.25_dp
    ptable(64) % covalent_radius = 1.6100_dp
    ptable(64) % vdw_radius = 2.00_dp
    ptable(64) % e_conv(0:3) = (/ 12, 24, 21, 7 /)
    ptable(64) % eht_param(0:3) = (/ z, z, z, z/)

    ! Terbium
    ptable(65) % symbol = 'Tb'
    ptable(65) % name = 'Terbium'
    ptable(65) % number = 65
    ptable(65) % amass = 158.92535_dp
    ptable(65) % covalent_radius = 1.5900_dp
    ptable(65) % vdw_radius = 2.00_dp
    ptable(65) % e_conv(0:3) = (/ 12, 24, 20, 9 /)
    ptable(65) % eht_param(0:3) = (/ z, z, z, z/)

    ! Dysprosium
    ptable(66) % symbol = 'Dy'
    ptable(66) % name = 'Dysprosium'
    ptable(66) % number = 66
    ptable(66) % amass = 162.5_dp
    ptable(66) % covalent_radius = 1.5900_dp
    ptable(66) % vdw_radius = 2.00_dp
    ptable(66) % e_conv(0:3) = (/ 12, 24, 20, 10 /)
    ptable(66) % eht_param(0:3) = (/ z, z, z, z/)

    ! Holmium
    ptable(67) % symbol = 'Ho'
    ptable(67) % name = 'Holmium'
    ptable(67) % number = 67
    ptable(67) % amass = 164.93032_dp
    ptable(67) % covalent_radius = 1.5800_dp
    ptable(67) % vdw_radius = 2.00_dp
    ptable(67) % e_conv(0:3) = (/ 12, 24, 20, 11 /)
    ptable(67) % eht_param(0:3) = (/ z, z, z, z/)

    ! Erbium
    ptable(68) % symbol = 'Er'
    ptable(68) % name = 'Erbium'
    ptable(68) % number = 68
    ptable(68) % amass = 167.259_dp
    ptable(68) % covalent_radius = 1.5700_dp
    ptable(68) % vdw_radius = 2.00_dp
    ptable(68) % e_conv(0:3) = (/ 12, 24, 20, 12 /)
    ptable(68) % eht_param(0:3) = (/ z, z, z, z/)

    ! Thulium
    ptable(69) % symbol = 'Tm'
    ptable(69) % name = 'Thulium'
    ptable(69) % number = 69
    ptable(69) % amass = 168.93421_dp
    ptable(69) % covalent_radius = 1.5600_dp
    ptable(69) % vdw_radius = 2.00_dp
    ptable(69) % e_conv(0:3) = (/ 12, 24, 20, 13 /)
    ptable(69) % eht_param(0:3) = (/ z, z, z, z/)

    ! Ytterbium
    ptable(70) % symbol = 'Yb'
    ptable(70) % name = 'Ytterbium'
    ptable(70) % number = 70
    ptable(70) % amass = 173.054_dp
    ptable(70) % covalent_radius = 1.5600_dp
    ptable(70) % vdw_radius = 2.00_dp
    ptable(70) % e_conv(0:3) = (/ 12, 24, 20, 14 /)
    ptable(70) % eht_param(0:3) = (/ -5.35_dp, -5.35_dp, -5.21_dp, &
         -13.86_dp/)

    ! Lutetium
    ptable(71) % symbol = 'Lu'
    ptable(71) % name = 'Lutetium'
    ptable(71) % number = 71
    ptable(71) % amass = 174.9668_dp
    ptable(71) % covalent_radius = 1.5600_dp
    ptable(71) % vdw_radius = 2.00_dp
    ptable(71) % e_conv(0:3) = (/ 12, 24, 21, 14 /)
    ptable(71) % eht_param(0:3) = (/ -6.05_dp, -6.05_dp, -5.12_dp, -22.40_dp/)

    ! Hafnium
    ptable(72) % symbol = 'Hf'
    ptable(72) % name = 'Hafnium'
    ptable(72) % number = 72
    ptable(72) % amass = 178.49_dp
    ptable(72) % covalent_radius = 1.4400_dp
    ptable(72) % vdw_radius = 2.00_dp
    ptable(72) % e_conv(0:3) = (/ 12, 24, 22, 14 /)
    ptable(72) % eht_param(0:3) = (/ z, z, z, z/)

    ! Tantalum
    ptable(73) % symbol = 'Ta'
    ptable(73) % name = 'Tantalum'
    ptable(73) % number = 73
    ptable(73) % amass = 180.94788_dp
    ptable(73) % covalent_radius = 1.3400_dp
    ptable(73) % vdw_radius = 2.00_dp
    ptable(73) % e_conv(0:3) = (/ 12, 24, 23, 14 /)
    ptable(73) % eht_param(0:3) = (/ -10.10_dp, -6.86_dp, -12.10_dp, z/)

    ! Tungsten
    ptable(74) % symbol = 'W '
    ptable(74) % name = 'Tungsten'
    ptable(74) % number = 74
    ptable(74) % amass = 183.84_dp
    ptable(74) % covalent_radius = 1.3000_dp
    ptable(74) % vdw_radius = 2.00_dp
    ptable(74) % e_conv(0:3) = (/ 12, 24, 24, 14 /)
    ptable(74) % eht_param(0:3) = (/ -8.26_dp, -5.17_dp, -10.37_dp, z/)

    ! Rhenium
    ptable(75) % symbol = 'Re'
    ptable(75) % name = 'Rhenium'
    ptable(75) % number = 75
    ptable(75) % amass = 186.207_dp
    ptable(75) % covalent_radius = 1.2800_dp
    ptable(75) % vdw_radius = 2.00_dp
    ptable(75) % e_conv(0:3) = (/ 12, 24, 25, 14 /)
    ptable(75) % eht_param(0:3) = (/ -9.36_dp, -5.96_dp, -12.66_dp, z/)

    ! Osmium
    ptable(76) % symbol = 'Os'
    ptable(76) % name = 'Osmium'
    ptable(76) % number = 76
    ptable(76) % amass = 190.23_dp
    ptable(76) % covalent_radius = 1.2600_dp
    ptable(76) % vdw_radius = 2.00_dp
    ptable(76) % e_conv(0:3) = (/ 12, 24, 26, 14 /)
    ptable(76) % eht_param(0:3) = (/ -8.17_dp, -4.81_dp, -11.84_dp, z/)

    ! Iridium
    ptable(77) % symbol = 'Ir'
    ptable(77) % name = 'Iridium'
    ptable(77) % number = 77
    ptable(77) % amass = 192.217_dp
    ptable(77) % covalent_radius = 1.2700_dp
    ptable(77) % vdw_radius = 2.00_dp
    ptable(77) % e_conv(0:3) = (/ 12, 24, 27, 14 /)
    ptable(77) % eht_param(0:3) = (/ -11.36_dp, -4.50_dp, -12.17_dp, z/)

    ! Platinum
    ptable(78) % symbol = 'Pt'
    ptable(78) % name = 'Platinum'
    ptable(78) % number = 78
    ptable(78) % amass = 195.084_dp
    ptable(78) % covalent_radius = 1.3000_dp
    ptable(78) % vdw_radius = 1.72_dp
    ptable(78) % e_conv(0:3) = (/ 11, 24, 29, 14 /)
    ptable(78) % eht_param(0:3) = (/ -9.077_dp, -5.475_dp, -12.59_dp, z/)

    ! Gold
    ptable(79) % symbol = 'Au'
    ptable(79) % name = 'Gold'
    ptable(79) % number = 79
    ptable(79) % amass = 196.966569_dp
    ptable(79) % covalent_radius = 1.3400_dp
    ptable(79) % vdw_radius = 1.66_dp
    ptable(79) % e_conv(0:3) = (/ 11, 24, 30, 14 /)
    ptable(79) % eht_param(0:3) = (/ -10.92_dp, -5.55_dp, -15.076_dp, z/)

    ! Mercury
    ptable(80) % symbol = 'Hg'
    ptable(80) % name = 'Mercury'
    ptable(80) % number = 80
    ptable(80) % amass = 200.59_dp
    ptable(80) % covalent_radius = 1.4900_dp
    ptable(80) % vdw_radius = 1.55_dp
    ptable(80) % e_conv(0:3) = (/ 12, 24, 30, 14 /)
    ptable(80) % eht_param(0:3) = (/ -13.68_dp, -8.47_dp, -17.50_dp, z/)

    ! Thallium
    ptable(81) % symbol = 'Tl'
    ptable(81) % name = 'Thallium'
    ptable(81) % number = 81
    ptable(81) % amass = 204.3833_dp
    ptable(81) % covalent_radius = 1.4800_dp
    ptable(81) % vdw_radius = 1.96_dp
    ptable(81) % e_conv(0:3) = (/ 12, 25, 30, 14 /)
    ptable(81) % eht_param(0:3) = (/ -11.60_dp, -5.80_dp, z, z/)

    ! Lead
    ptable(82) % symbol = 'Pb'
    ptable(82) % name = 'Lead'
    ptable(82) % number = 82
    ptable(82) % amass = 207.2_dp
    ptable(82) % covalent_radius = 1.4700_dp
    ptable(82) % vdw_radius = 2.02_dp
    ptable(82) % e_conv(0:3) = (/ 12, 26, 30, 14 /)
    ptable(82) % eht_param(0:3) = (/ -15.70_dp, -8.00_dp, z, z/)

    ! Bismuth
    ptable(83) % symbol = 'Bi'
    ptable(83) % name = 'Bismuth'
    ptable(83) % number = 83
    ptable(83) % amass = 208.9804_dp
    ptable(83) % covalent_radius = 1.4600_dp
    ptable(83) % vdw_radius = 2.00_dp
    ptable(83) % e_conv(0:3) = (/ 12, 27, 30, 14 /)
    ptable(83) % eht_param(0:3) = (/ -15.19_dp, -7.79_dp, z, z/)

    ! Polonium
    ptable(84) % symbol = 'Po'
    ptable(84) % name = 'Polonium'
    ptable(84) % number = 84
    ptable(84) % amass = 208.9824_dp ! (209)
    ptable(84) % covalent_radius = 1.4600_dp
    ptable(84) % vdw_radius = 2.00_dp
    ptable(84) % e_conv(0:3) = (/ 12, 28, 30, 14 /)
    ptable(84) % eht_param(0:3) = (/ z, z, z, z/)

    ! Astatine
    ptable(85) % symbol = 'At'
    ptable(85) % name = 'Astatine'
    ptable(85) % number = 85
    ptable(85) % amass = 209.9871_dp ! (210)
    ptable(85) % covalent_radius = 1.4500_dp
    ptable(85) % vdw_radius = 2.00_dp
    ptable(85) % e_conv(0:3) = (/ 12, 29, 30, 14 /)
    ptable(85) % eht_param(0:3) = (/ z, z, z, z/)

    ! Radon
    ptable(86) % symbol = 'Rn'
    ptable(86) % name = 'Radon'
    ptable(86) % number = 86
    ptable(86) % amass = 222.0176_dp ! (222)
    ptable(86) % covalent_radius = 1.50_dp
    ptable(86) % vdw_radius = 2.00_dp
    ptable(86) % e_conv(0:3) = (/ 12, 30, 30, 14 /)
    ptable(86) % eht_param(0:3) = (/ z, z, z, z/)

    ! Francium
    ptable(87) % symbol = 'Fr'
    ptable(87) % name = 'Francium'
    ptable(87) % number = 87
    ptable(87) % amass = 223.0197_dp ! (223)
    ptable(87) % covalent_radius = 1.50_dp
    ptable(87) % vdw_radius = 2.00_dp
    ptable(87) % e_conv(0:3) = (/ 13, 30, 30, 14 /)
    ptable(87) % eht_param(0:3) = (/ z, z, z, z/)

    ! Radium
    ptable(88) % symbol = 'Ra'
    ptable(88) % name = 'Radium'
    ptable(88) % number = 88
    ptable(88) % amass = 226.0254_dp ! (226)
    ptable(88) % covalent_radius = 1.90_dp
    ptable(88) % vdw_radius = 2.00_dp
    ptable(88) % e_conv(0:3) = (/ 14, 30, 30, 14 /)
    ptable(88) % eht_param(0:3) = (/ z, z, z, z/)

    ! Actinium
    ptable(89) % symbol = 'Ac'
    ptable(89) % name = 'Actinium'
    ptable(89) % number = 89
    ptable(89) % amass = 227.0278_dp ! (227)
    ptable(89) % covalent_radius = 1.88_dp
    ptable(89) % vdw_radius = 2.00_dp
    ptable(89) % e_conv(0:3) = (/ 14, 30, 31, 14 /)
    ptable(89) % eht_param(0:3) = (/ z, z, z, z/)

    ! Thorium
    ptable(90) % symbol = 'Th'
    ptable(90) % name = 'Thorium'
    ptable(90) % number = 90
    ptable(90) % amass = 232.03806_dp
    ptable(90) % covalent_radius = 1.6500_dp
    ptable(90) % vdw_radius = 2.00_dp
    ptable(90) % e_conv(0:3) = (/ 14, 30, 32, 14 /)
    ptable(90) % eht_param(0:3) = (/ -5.39_dp, -5.39_dp, -10.11_dp, -9.64_dp/)

    ! Proctactinium
    ptable(91) % symbol = 'Pa'
    ptable(91) % name = 'Proctactinium'
    ptable(91) % number = 91
    ptable(91) % amass = 231.03588_dp
    ptable(91) % covalent_radius = 1.61_dp
    ptable(91) % vdw_radius = 2.00_dp
    ptable(91) % e_conv(0:3) = (/ 14, 30, 31, 16 /)
    ptable(91) % eht_param(0:3) = (/ z, z, z, z/)

    ! Uranium
    ptable(92) % symbol = 'U '
    ptable(92) % name = 'Uranium'
    ptable(92) % number = 92
    ptable(92) % amass = 238.02891_dp
    ptable(92) % covalent_radius = 1.4200_dp
    ptable(92) % vdw_radius = 1.86_dp
    ptable(92) % e_conv(0:3) = (/ 14, 30, 31, 17 /)
    ptable(92) % eht_param(0:3) = (/ -5.50_dp, -5.50_dp, -9.19_dp, -10.62_dp/)

    ! Neptunium
    ptable(93) % symbol = 'Np'
    ptable(93) % name = 'Neptunium'
    ptable(93) % number = 93
    ptable(93) % amass = 237.0482_dp ! (237)
    ptable(93) % covalent_radius = 1.55_dp
    ptable(93) % vdw_radius = 2.00_dp
    ptable(93) % e_conv(0:3) = (/ 14, 30, 31, 18 /)
    ptable(93) % eht_param(0:3) = (/ z, z, z, z/)

    ! Plutonium
    ptable(94) % symbol = 'Pu'
    ptable(94) % name = 'Plutonium'
    ptable(94) % number = 94
    ptable(94) % amass = 244.0642_dp ! (244)
    ptable(94) % covalent_radius = 1.53_dp
    ptable(94) % vdw_radius = 2.00_dp
    ptable(94) % e_conv(0:3) = (/ 14, 30, 30, 20 /)
    ptable(94) % eht_param(0:3) = (/ z, z, z, z/)

    ! Americum
    ptable(95) % symbol = 'Am'
    ptable(95) % name = 'Americum'
    ptable(95) % number = 95
    ptable(95) % amass = 243.0614_dp ! (243)
    ptable(95) % covalent_radius = 1.51_dp
    ptable(95) % vdw_radius = 2.00_dp
    ptable(95) % e_conv(0:3) = (/ 14, 30, 30, 21 /)
    ptable(95) % eht_param(0:3) = (/ z, z, z, z/)

    ! Curium
    ptable(96) % symbol = 'Cm'
    ptable(96) % name = 'Curium'
    ptable(96) % number = 96
    ptable(96) % amass = 247.0704_dp ! (247)
    ptable(96) % covalent_radius = 0.99_dp
    ptable(96) % vdw_radius = 2.00_dp
    ptable(96) % e_conv(0:3) = (/ 14, 30, 31, 21 /)
    ptable(96) % eht_param(0:3) = (/ z, z, z, z/)

    ! Berkelium
    ptable(97) % symbol = 'Bk'
    ptable(97) % name = 'Berkelium'
    ptable(97) % number = 97
    ptable(97) % amass = 247.0703_dp ! (247)
    ptable(97) % covalent_radius = 1.54_dp
    ptable(97) % vdw_radius = 2.00_dp
    ptable(97) % e_conv(0:3) = (/ 14, 30, 30, 23 /)
    ptable(97) % eht_param(0:3) = (/ z, z, z, z/)

    ! Californium
    ptable(98) % symbol = 'Cf'
    ptable(98) % name = 'Californium'
    ptable(98) % number = 98
    ptable(98) % amass = 251.0796_dp ! (251)
    ptable(98) % covalent_radius = 1.83_dp
    ptable(98) % vdw_radius = 2.00_dp
    ptable(98) % e_conv(0:3) = (/ 14, 30, 30, 24 /)
    ptable(98) % eht_param(0:3) = (/ z, z, z, z/)

    ! Einsteinium
    ptable(99) % symbol = 'Es'
    ptable(99) % name = 'Einsteinium'
    ptable(99) % number = 99
    ptable(99) % amass = 252.083_dp ! (252)
    ptable(99) % covalent_radius = 1.50_dp
    ptable(99) % vdw_radius = 2.00_dp
    ptable(99) % e_conv(0:3) = (/ 14, 30, 30, 25 /)
    ptable(99) % eht_param(0:3) = (/ z, z, z, z/)

    ! Fermium
    ptable(100) % symbol = 'Fm'
    ptable(100) % name = 'Fermium'
    ptable(100) % number = 100
    ptable(100) % amass = 257.0951_dp ! (257)
    ptable(100) % covalent_radius = 1.50_dp
    ptable(100) % vdw_radius = 2.00_dp
    ptable(100) % e_conv(0:3) = (/ 14, 30, 30, 26 /)
    ptable(100) % eht_param(0:3) = (/ z, z, z, z/)

    ! Mendelevium
    ptable(101) % symbol = 'Md'
    ptable(101) % name = 'Mendelevium'
    ptable(101) % number = 101
    ptable(101) % amass = 258.0984_dp ! (258)
    ptable(101) % covalent_radius = 1.50_dp
    ptable(101) % vdw_radius = 2.00_dp
    ptable(101) % e_conv(0:3) = (/ 14, 30, 30, 27 /)
    ptable(101) % eht_param(0:3) = (/ z, z, z, z/)

    ! Nobelium
    ptable(102) % symbol = 'No'
    ptable(102) % name = 'Nobelium'
    ptable(102) % number = 102
    ptable(102) % amass = 259.101_dp ! (259)
    ptable(102) % covalent_radius = 1.50_dp
    ptable(102) % vdw_radius = 2.00_dp
    ptable(102) % e_conv(0:3) = (/ 14, 30, 30, 28 /)
    ptable(102) % eht_param(0:3) = (/ z, z, z, z/)

    ! Lawrencium
    ptable(103) % symbol = 'Lr'
    ptable(103) % name = 'Lawrencium'
    ptable(103) % number = 103
    ptable(103) % amass = 262.1096_dp ! (262)
    ptable(103) % covalent_radius = 1.50_dp
    ptable(103) % vdw_radius = 2.00_dp
    ptable(103) % e_conv(0:3) = (/ 14, 31, 30, 28 /)
    ptable(103) % eht_param(0:3) = (/ z, z, z, z/)

    ! Rutherfordium
    ptable(104) % symbol = 'Rf'
    ptable(104) % name = 'Rutherfordium'
    ptable(104) % number = 104
    ptable(104) % amass = 267.122_dp ! (267)
    ptable(104) % covalent_radius = 1.50_dp
    ptable(104) % vdw_radius = 2.00_dp
    ptable(104) % e_conv(0:3) = (/ 14, 30, 32, 28 /)
    ptable(104) % eht_param(0:3) = (/ z, z, z, z/)

    ! Dubnium
    ptable(105) % symbol = 'Db'
    ptable(105) % name = 'Dubnium'
    ptable(105) % number = 105
    ptable(105) % amass = 268.125_dp ! (268)
    ptable(105) % covalent_radius = 1.50_dp
    ptable(105) % vdw_radius = 2.00_dp
    ptable(105) % e_conv(0:3) = (/ 14, 30, 33, 28 /)
    ptable(105) % eht_param(0:3) = (/ z, z, z, z/)

    ! Seaborgium
    ptable(106) % symbol = 'Sg'
    ptable(106) % name = 'Seaborgium'
    ptable(106) % number = 106
    ptable(106) % amass = 271.133_dp ! (271)
    ptable(106) % covalent_radius = 1.50_dp
    ptable(106) % vdw_radius = 2.00_dp
    ptable(106) % e_conv(0:3) = (/ 14, 30, 34, 28 /)
    ptable(106) % eht_param(0:3) = (/ z, z, z, z/)

    ! Bohrium
    ptable(107) % symbol = 'Bh'
    ptable(107) % name = 'Bohrium'
    ptable(107) % number = 107
    ptable(107) % amass = 272.0_dp ! (272)
    ptable(107) % covalent_radius = 1.50_dp
    ptable(107) % vdw_radius = 2.00_dp
    ptable(107) % e_conv(0:3) = (/ 14, 30, 35, 28 /)
    ptable(107) % eht_param(0:3) = (/ z, z, z, z/)

    ! Hassium
    ptable(108) % symbol = 'Hs'
    ptable(108) % name = 'Hassium'
    ptable(108) % number = 108
    ptable(108) % amass = 277.15_dp ! (277)
    ptable(108) % covalent_radius = 1.50_dp
    ptable(108) % vdw_radius = 2.00_dp
    ptable(108) % e_conv(0:3) = (/ 14, 30, 36, 28 /)
    ptable(108) % eht_param(0:3) = (/ z, z, z, z/)

    ! Meitnerium
    ptable(109) % symbol = 'Mt'
    ptable(109) % name = 'Meitnerium'
    ptable(109) % number = 109
    ptable(109) % amass = 276.151_dp ! (276)
    ptable(109) % covalent_radius = 1.50_dp
    ptable(109) % vdw_radius = 2.00_dp
    ptable(109) % e_conv(0:3) = (/ 14, 30, 37, 28 /)
    ptable(109) % eht_param(0:3) = (/ z, z, z, z/)

    ! Darmstadtium
    ptable(110) % symbol = 'Ds'
    ptable(110) % name = 'Darmstadtium'
    ptable(110) % number = 110
    ptable(110) % amass = 281.162_dp ! (281)
    ptable(110) % covalent_radius = 1.50_dp
    ptable(110) % vdw_radius = 2.00_dp
    ptable(110) % e_conv(0:3) = (/ 13, 30, 39, 28 /)
    ptable(110) % eht_param(0:3) = (/ z, z, z, z/)

    ! Roentgenium
    ptable(111) % symbol = 'Rg'
    ptable(111) % name = 'Roentgenium'
    ptable(111) % number = 111
    ptable(111) % amass = 280.164_dp ! (280)
    ptable(111) % covalent_radius = 1.50_dp
    ptable(111) % vdw_radius = 2.00_dp
    ptable(111) % e_conv(0:3) = (/ 13, 30, 40, 28 /)
    ptable(111) % eht_param(0:3) = (/ z, z, z, z/)

    ! Copernicium
    ptable(112) % symbol = 'Cn'
    ptable(112) % name = 'Copernicium'
    ptable(112) % number = 112
    ptable(112) % amass = 285.174_dp ! (285)
    ptable(112) % covalent_radius = 1.50_dp
    ptable(112) % vdw_radius = 2.00_dp
    ptable(112) % e_conv(0:3) = (/ 14, 30, 40, 28 /)
    ptable(112) % eht_param(0:3) = (/ z, z, z, z/)

    ! Ununtrium
    ptable(113) % symbol = 'Uut'
    ptable(113) % name = 'Ununtrium'
    ptable(113) % number = 113
    ptable(113) % amass = 284.178_dp ! (284)
    ptable(113) % covalent_radius = 1.50_dp
    ptable(113) % vdw_radius = 2.00_dp
    ptable(113) % e_conv(0:3) = (/ 14, 31, 40, 28 /)
    ptable(113) % eht_param(0:3) = (/ z, z, z, z/)

    ! Ununquadium
    ptable(114) % symbol = 'Uuq'
    ptable(114) % name = 'Ununquadium'
    ptable(114) % number = 114
    ptable(114) % amass = 289.189_dp ! (289)
    ptable(114) % covalent_radius = 1.50_dp
    ptable(114) % vdw_radius = 2.00_dp
    ptable(114) % e_conv(0:3) = (/ 14, 32, 40, 28 /)
    ptable(114) % eht_param(0:3) = (/ z, z, z, z/)

    ! Ununpentium
    ptable(115) % symbol = 'Uup'
    ptable(115) % name = 'Ununpentium'
    ptable(115) % number = 115
    ptable(115) % amass = 288.192_dp ! (288)
    ptable(115) % covalent_radius = 1.50_dp
    ptable(115) % vdw_radius = 2.00_dp
    ptable(115) % e_conv(0:3) = (/ 14, 33, 40, 28 /)
    ptable(115) % eht_param(0:3) = (/ z, z, z, z/)

    ! Ununhexium
    ptable(116) % symbol = 'Uuh'
    ptable(116) % name = 'Ununhexium'
    ptable(116) % number = 116
    ptable(116) % amass = 293.0_dp ! (293)
    ptable(116) % covalent_radius = 1.50_dp
    ptable(116) % vdw_radius = 2.00_dp
    ptable(116) % e_conv(0:3) = (/ 14, 34, 40, 28 /)
    ptable(116) % eht_param(0:3) = (/ z, z, z, z/)

    ! Ununseptium
    ptable(117) % symbol = 'Uus'
    ptable(117) % name = 'Ununseptium'
    ptable(117) % number = 117
    ptable(117) % amass = 294.0_dp ! (294)
    ptable(117) % covalent_radius = 1.50_dp
    ptable(117) % vdw_radius = 2.00_dp
    ptable(117) % e_conv(0:3) = (/ 14, 35, 40, 28 /)
    ptable(117) % eht_param(0:3) = (/ z, z, z, z/)

    ! Ununoctium
    ptable(118) % symbol = 'Uuo'
    ptable(118) % name = 'Ununoctium'
    ptable(118) % number = 118
    ptable(118) % amass = 294.0_dp ! (294)
    ptable(118) % covalent_radius = 1.50_dp
    ptable(118) % vdw_radius = 2.00_dp
    ptable(118) % e_conv(0:3) = (/ 14, 36, 40, 28 /)
    ptable(118) % eht_param(0:3) = (/ z, z, z, z/)

#if defined(__USE_LEGACY_WEIGHTS)
    ptable(1) % amass = 1.00797_dp
    ptable(2) % amass = 4.00260_dp
    ptable(3) % amass = 6.93900_dp
    ptable(4) % amass = 9.01220_dp
    ptable(5) % amass = 10.81100_dp
    ptable(6) % amass = 12.01115_dp
    ptable(7) % amass = 14.00670_dp
    ptable(8) % amass = 15.99940_dp
    ptable(9) % amass = 18.99840_dp
    ptable(10) % amass = 20.18300_dp
    ptable(11) % amass = 22.98980_dp
    ptable(12) % amass = 24.31200_dp
    ptable(13) % amass = 26.98153_dp
    ptable(14) % amass = 28.08600_dp
    ptable(15) % amass = 30.97380_dp
    ptable(16) % amass = 32.06400_dp
    ptable(17) % amass = 35.45300_dp
    ptable(18) % amass = 39.94800_dp
    ptable(19) % amass = 39.10200_dp
    ptable(20) % amass = 40.08000_dp
    ptable(21) % amass = 44.95600_dp
    ptable(22) % amass = 47.90000_dp
    ptable(23) % amass = 50.94200_dp
    ptable(24) % amass = 51.99600_dp
    ptable(25) % amass = 54.93800_dp
    ptable(26) % amass = 55.84700_dp
    ptable(27) % amass = 58.93300_dp
    ptable(28) % amass = 58.71000_dp
    ptable(29) % amass = 63.54000_dp
    ptable(30) % amass = 65.37000_dp
    ptable(31) % amass = 69.72000_dp
    ptable(32) % amass = 72.59000_dp
    ptable(33) % amass = 74.92200_dp
    ptable(34) % amass = 78.96000_dp
    ptable(35) % amass = 79.90900_dp
    ptable(36) % amass = 83.80000_dp
    ptable(37) % amass = 85.47000_dp
    ptable(38) % amass = 87.62000_dp
    ptable(39) % amass = 88.90500_dp
    ptable(40) % amass = 91.22000_dp
    ptable(41) % amass = 92.90600_dp
    ptable(42) % amass = 95.94000_dp
    ptable(43) % amass = 98.90600_dp
    ptable(44) % amass = 101.07000_dp
    ptable(45) % amass = 102.90500_dp
    ptable(46) % amass = 106.40000_dp
    ptable(47) % amass = 107.87000_dp
    ptable(48) % amass = 112.40000_dp
    ptable(49) % amass = 114.82000_dp
    ptable(50) % amass = 118.69000_dp
    ptable(51) % amass = 121.75000_dp
    ptable(52) % amass = 127.60000_dp
    ptable(53) % amass = 126.90440_dp
    ptable(54) % amass = 131.30000_dp
    ptable(55) % amass = 132.90500_dp
    ptable(56) % amass = 137.34000_dp
    ptable(57) % amass = 138.91000_dp
    ptable(58) % amass = 140.12000_dp
    ptable(59) % amass = 140.90700_dp
    ptable(60) % amass = 144.24000_dp
    ptable(61) % amass = 144.91300_dp
    ptable(62) % amass = 150.35000_dp
    ptable(63) % amass = 151.96000_dp
    ptable(64) % amass = 157.25000_dp
    ptable(65) % amass = 158.92400_dp
    ptable(66) % amass = 162.50000_dp
    ptable(67) % amass = 164.93000_dp
    ptable(68) % amass = 167.26000_dp
    ptable(69) % amass = 168.93400_dp
    ptable(70) % amass = 173.04000_dp
    ptable(71) % amass = 174.97000_dp
    ptable(72) % amass = 178.49000_dp
    ptable(73) % amass = 180.94800_dp
    ptable(74) % amass = 183.85000_dp
    ptable(75) % amass = 186.20000_dp
    ptable(76) % amass = 190.20000_dp
    ptable(77) % amass = 192.20000_dp
    ptable(78) % amass = 195.09000_dp
    ptable(79) % amass = 196.96700_dp
    ptable(80) % amass = 200.59000_dp
    ptable(81) % amass = 204.37000_dp
    ptable(82) % amass = 207.19000_dp
    ptable(83) % amass = 208.98000_dp
    ptable(84) % amass = 209.98290_dp
    ptable(85) % amass = 209.98700_dp
    ptable(86) % amass = 222.01750_dp
    ptable(87) % amass = 223.01980_dp
    ptable(88) % amass = 226.02540_dp
    ptable(89) % amass = 227.02780_dp
    ptable(90) % amass = 232.03810_dp
    ptable(91) % amass = 231.03590_dp
    ptable(92) % amass = 238.05080_dp
    ptable(93) % amass = 237.04820_dp
    ptable(94) % amass = 244.0640_dp
    ptable(95) % amass = 243.0614_dp
    ptable(96) % amass = 247.0700_dp
    ptable(97) % amass = 251.0800_dp
    ptable(98) % amass = 252.0820_dp
    ptable(99) % amass = 252.0829_dp
    ptable(100) % amass = 257.0950_dp
    ptable(101) % amass = 256.0000_dp
    ptable(102) % amass = 254.0000_dp
    ptable(103) % amass = 257.0000_dp
    ptable(104) % amass = 261.0000_dp
    ptable(105) % amass = 262.0000_dp
    ptable(106:nelem) % amass = 263.0000_dp
#endif

    ! Initialize heat of formation
    CALL init_eheat

    ! Initialize gyromagnetic ratio
    CALL init_gratio

  END SUBROUTINE init_periodic_table

! *****************************************************************************
  SUBROUTINE init_eheat
    ! All values in kcal/mol
    ! Dummy
    ptable(0) % heat_of_formation = z
    ! Hydrogen
    ptable(1) % heat_of_formation = 52.102_dp
    ! Helium
    ptable(2) % heat_of_formation = z
    ! Lithium
    ptable(3) % heat_of_formation = 38.410_dp
    ! Beryllium
    ptable(4) % heat_of_formation = 76.960_dp
    ! Boron
    ptable(5) % heat_of_formation = 135.700_dp
    ! Carbon
    ptable(6) % heat_of_formation = 170.890_dp
    ! Nitrogen
    ptable(7) % heat_of_formation = 113.000_dp
    ! Oxygen
    ptable(8) % heat_of_formation = 59.559_dp
    ! Fluorine
    ptable(9) % heat_of_formation = 18.890_dp
    ! Neon
    ptable(10) % heat_of_formation = z
    ! Sodium
    ptable(11) % heat_of_formation = 25.850_dp
    ! Magnesium
    ptable(12) % heat_of_formation = 35.000_dp
    ! Aluminium
    ptable(13) % heat_of_formation = 79.490_dp
    ! Silicon
    ptable(14) % heat_of_formation = 108.390_dp
    ! Phosphorus
    ptable(15) % heat_of_formation = 75.570_dp
    ! Sulfur
    ptable(16) % heat_of_formation = 66.400_dp
    ! Chlorine
    ptable(17) % heat_of_formation = 28.990_dp
    ! Argon
    ptable(18) % heat_of_formation = z
    ! Potassium
    ptable(19) % heat_of_formation = 21.420_dp
    ! Calcium
    ptable(20) % heat_of_formation = 42.600_dp
    ! Scandium
    ptable(21) % heat_of_formation = 90.300_dp
    ! Titanium
    ptable(22) % heat_of_formation = 112.300_dp
    ! Vanadium
    ptable(23) % heat_of_formation = 122.900_dp
    ! Chromium
    ptable(24) % heat_of_formation = 95.000_dp
    ! Manganese
    ptable(25) % heat_of_formation = 67.700_dp
    ! Iron
    ptable(26) % heat_of_formation = 99.300_dp
    ! Cobalt
    ptable(27) % heat_of_formation = 102.400_dp
    ! Nickel
    ptable(28) % heat_of_formation = 102.800_dp
    ! Copper
    ptable(29) % heat_of_formation = 80.700_dp
    ! Zinc
    ptable(30) % heat_of_formation = 31.170_dp
    ! Gallium
    ptable(31) % heat_of_formation = 65.400_dp
    ! Germanium
    ptable(32) % heat_of_formation = 89.500_dp
    ! Arsenic
    ptable(33) % heat_of_formation = 72.300_dp
    ! Selenium
    ptable(34) % heat_of_formation = 54.300_dp
    ! Bromine
    ptable(35) % heat_of_formation = 26.740_dp
    ! Krypton
    ptable(36) % heat_of_formation = z
    ! Rubidium
    ptable(37) % heat_of_formation = 19.600_dp
    ! Strontium
    ptable(38) % heat_of_formation = 39.100_dp
    ! Yttrium
    ptable(39) % heat_of_formation = 101.500_dp
    ! Zirconium
    ptable(40) % heat_of_formation = 145.500_dp
    ! Niobium
    ptable(41) % heat_of_formation = 172.400_dp
    ! Molybdenum
    ptable(42) % heat_of_formation = 157.300_dp
    ! Technetium
    ptable(43) % heat_of_formation = z
    ! Ruthenium
    ptable(44) % heat_of_formation = 155.500_dp
    ! Rhodium
    ptable(45) % heat_of_formation = 133.000_dp
    ! Palladium
    ptable(46) % heat_of_formation = 90.000_dp
    ! Silver
    ptable(47) % heat_of_formation = 68.100_dp
    ! Cadmium
    ptable(48) % heat_of_formation = 26.720_dp
    ! Indium
    ptable(49) % heat_of_formation = 58.000_dp
    ! Tin
    ptable(50) % heat_of_formation = 72.200_dp
    ! Antimony
    ptable(51) % heat_of_formation = 63.200_dp
    ! Tellurium
    ptable(52) % heat_of_formation = 47.000_dp
    ! Iodine
    ptable(53) % heat_of_formation = 25.517_dp
    ! Xenon
    ptable(54) % heat_of_formation = z
    ! Cesium
    ptable(55) % heat_of_formation = 18.700_dp
    ! Barium
    ptable(56) % heat_of_formation = 42.500_dp
    ! Lanthanum
    ptable(57) % heat_of_formation = z
    ! Cerium
    ptable(58) % heat_of_formation = 101.300_dp
    ! Praseodymium
    ptable(59) % heat_of_formation = z
    ! Neodymium
    ptable(60) % heat_of_formation = z
    ! Promethium
    ptable(61) % heat_of_formation = z
    ! Samarium
    ptable(62) % heat_of_formation = 49.400_dp
    ! Europium
    ptable(63) % heat_of_formation = z
    ! Gadolinium
    ptable(64) % heat_of_formation = z
    ! Terbium
    ptable(65) % heat_of_formation = z
    ! Dysprosium
    ptable(66) % heat_of_formation = z
    ! Holmium
    ptable(67) % heat_of_formation = z
    ! Erbium
    ptable(68) % heat_of_formation = 75.800_dp
    ! Thulium
    ptable(69) % heat_of_formation = z
    ! Ytterbium
    ptable(70) % heat_of_formation = 36.350_dp
    ! Lutetium
    ptable(71) % heat_of_formation = z
    ! Hafnium
    ptable(72) % heat_of_formation = 148.000_dp
    ! Tantalum
    ptable(73) % heat_of_formation = 186.900_dp
    ! Tungsten
    ptable(74) % heat_of_formation = 203.100_dp
    ! Rhenium
    ptable(75) % heat_of_formation = 185.000_dp
    ! Osmium
    ptable(76) % heat_of_formation = 188.000_dp
    ! Iridium
    ptable(77) % heat_of_formation = 160.000_dp
    ! Platinum
    ptable(78) % heat_of_formation = 135.200_dp
    ! Gold
    ptable(79) % heat_of_formation = 88.000_dp
    ! Mercury
    ptable(80) % heat_of_formation = 14.690_dp
    ! Thallium
    ptable(81) % heat_of_formation = 43.550_dp
    ! Lead
    ptable(82) % heat_of_formation = 46.620_dp
    ! Bismuth
    ptable(83) % heat_of_formation = 50.100_dp
    ! Polonium
    ptable(84) % heat_of_formation = z
    ! Astatine
    ptable(85) % heat_of_formation = z
    ! Radon
    ! from Radon no parametrisation in dynamo
    ptable(86:nelem) % heat_of_formation = z

  END SUBROUTINE init_eheat

! *****************************************************************************
  SUBROUTINE init_gratio
    ! D. M. Granty and R. K. Harris, Encyclopedia of Nuclear
    ! Magnetic Resonance vol. 5 (Wiley, Chichester UK, 1996)
    !
    ! gyrom_ratio values in MHz/Tesla
    !
    ! Dummy
    ptable(0) % gyrom_ratio = 0.0_dp
    ptable(0) % gyrom_ratio_isotope = 0
    ! Hydrogen
    ptable(1) % gyrom_ratio = 42.5774690577_dp
    ptable(1) % gyrom_ratio_isotope = 1
    ! Helium
    ptable(2) % gyrom_ratio = -32.4360299810_dp
    ptable(2) % gyrom_ratio_isotope = 3
    ! Lithium
    ptable(3) % gyrom_ratio = 16.5484555869_dp
    ptable(3) % gyrom_ratio_isotope = 7
    ! Beryllium
    ptable(4) % gyrom_ratio = -5.9836942827_dp
    ptable(4) % gyrom_ratio_isotope = 9
    ! Boron
    ptable(5) % gyrom_ratio = 13.6629814024_dp
    ptable(5) % gyrom_ratio_isotope = 11
    ! Carbon
    ptable(6) % gyrom_ratio = 10.7083965713_dp
    ptable(6) % gyrom_ratio_isotope = 13
    ! Nitrogen
    ptable(7) % gyrom_ratio = 3.0777051853_dp
    ptable(7) % gyrom_ratio_isotope = 14
    ! Oxygen
    ptable(8) % gyrom_ratio = -5.7742686593_dp
    ptable(8) % gyrom_ratio_isotope = 17
    ! Fluorine
    ptable(9) % gyrom_ratio = 40.0775701637_dp
    ptable(9) % gyrom_ratio_isotope = 19
    ! Neon
    ptable(10) % gyrom_ratio = -3.3630712715_dp
    ptable(10) % gyrom_ratio_isotope = 21
    ! Sodium
    ptable(11) % gyrom_ratio = 11.2695216738_dp
    ptable(11) % gyrom_ratio_isotope = 23
    ! Magnesium
    ptable(12) % gyrom_ratio = -2.6083426159_dp
    ptable(12) % gyrom_ratio_isotope = 25
    ! Aluminium
    ptable(13) % gyrom_ratio = 11.1030809358_dp
    ptable(13) % gyrom_ratio_isotope = 27
    ! Silicon
    ptable(14) % gyrom_ratio = -8.4654514231_dp
    ptable(14) % gyrom_ratio_isotope = 29
    ! Phosphorus
    ptable(15) % gyrom_ratio = 17.2514409015_dp
    ptable(15) % gyrom_ratio_isotope = 31
    ! Sulfur
    ptable(16) % gyrom_ratio = 3.2717242919_dp
    ptable(16) % gyrom_ratio_isotope = 33
    ! Chlorine
    ptable(17) % gyrom_ratio = 4.1765408335_dp
    ptable(17) % gyrom_ratio_isotope = 35
    ! Argon
    ptable(18) % gyrom_ratio = 0.0_dp
    ptable(18) % gyrom_ratio_isotope = 0
    ! Potassium
    ptable(19) % gyrom_ratio = 1.9895335549_dp
    ptable(19) % gyrom_ratio_isotope = 39
    ! Calcium
    ptable(20) % gyrom_ratio = -2.8696734409_dp
    ptable(20) % gyrom_ratio_isotope = 43
    ! Scandium
    ptable(21) % gyrom_ratio = 10.3590726388_dp
    ptable(21) % gyrom_ratio_isotope = 45
    ! Titanium
    ptable(22) % gyrom_ratio = -2.4040354154_dp
    ptable(22) % gyrom_ratio_isotope = 47
    ! Vanadium
    ptable(23) % gyrom_ratio = 11.2132801367_dp
    ptable(23) % gyrom_ratio_isotope = 51
    ! Chromium
    ptable(24) % gyrom_ratio = -2.4115156977_dp
    ptable(24) % gyrom_ratio_isotope = 53
    ! Manganese
    ptable(25) % gyrom_ratio = 10.5762511769_dp
    ptable(25) % gyrom_ratio_isotope = 55
    ! Iron
    ptable(26) % gyrom_ratio = 1.3815642187_dp
    ptable(26) % gyrom_ratio_isotope = 57
    ! Cobalt
    ptable(27) % gyrom_ratio = 10.0776909966_dp
    ptable(27) % gyrom_ratio_isotope = 59
    ! Nickel
    ptable(28) % gyrom_ratio = -3.8114425772_dp
    ptable(28) % gyrom_ratio_isotope = 61
    ! Copper
    ptable(29) % gyrom_ratio = 11.3187637358_dp
    ptable(29) % gyrom_ratio_isotope = 63
    ! Zinc
    ptable(30) % gyrom_ratio = 2.6685318322_dp
    ptable(30) % gyrom_ratio_isotope = 67
    ! Gallium
    ptable(31) % gyrom_ratio = 10.2477560110_dp
    ptable(31) % gyrom_ratio_isotope = 69
    ! Germanium
    ptable(32) % gyrom_ratio = -1.4897384913_dp
    ptable(32) % gyrom_ratio_isotope = 73
    ! Arsenic
    ptable(33) % gyrom_ratio = 7.3150206071_dp
    ptable(33) % gyrom_ratio_isotope = 75
    ! Selenium
    ptable(34) % gyrom_ratio = 8.1573046941_dp
    ptable(34) % gyrom_ratio_isotope = 77
    ! Bromine
    ptable(35) % gyrom_ratio = 10.7041503174_dp
    ptable(35) % gyrom_ratio_isotope = 79
    ! Krypton
    ptable(36) % gyrom_ratio = -1.6442297171_dp
    ptable(36) % gyrom_ratio_isotope = 83
    ! Rubidium
    ptable(37) % gyrom_ratio = 4.1264181673_dp
    ptable(37) % gyrom_ratio_isotope = 85
    ! Strontium
    ptable(38) % gyrom_ratio = -1.8524642249_dp
    ptable(38) % gyrom_ratio_isotope = 87
    ! Yttrium
    ptable(39) % gyrom_ratio = -2.0949232525_dp
    ptable(39) % gyrom_ratio_isotope = 89
    ! Zirconium
    ptable(40) % gyrom_ratio = -3.9747832953_dp
    ptable(40) % gyrom_ratio_isotope = 91
    ! Niobium
    ptable(41) % gyrom_ratio = 10.4523417326_dp
    ptable(41) % gyrom_ratio_isotope = 93
    ! Molybdenum
    ptable(42) % gyrom_ratio = -2.7868030535_dp
    ptable(42) % gyrom_ratio_isotope = 95
    ! Technetium
    ptable(43) % gyrom_ratio = 9.6225078593_dp
    ptable(43) % gyrom_ratio_isotope = 99
    ! Ruthenium
    ptable(44) % gyrom_ratio = -2.1915635664_dp
    ptable(44) % gyrom_ratio_isotope = 101
    ! Rhodium
    ptable(45) % gyrom_ratio = -1.3477240581_dp
    ptable(45) % gyrom_ratio_isotope = 103
    ! Palladium
    ptable(46) % gyrom_ratio = -1.9576058000_dp
    ptable(46) % gyrom_ratio_isotope = 105
    ! Silver
    ptable(47) % gyrom_ratio = -1.7330669824_dp
    ptable(47) % gyrom_ratio_isotope = 107
    ! Cadmium
    ptable(48) % gyrom_ratio = -9.0691469715_dp
    ptable(48) % gyrom_ratio_isotope = 111
    ! Indium
    ptable(49) % gyrom_ratio = 9.3856853040_dp
    ptable(49) % gyrom_ratio_isotope = 115
    ! Tin
    ptable(50) % gyrom_ratio = -15.9659464261_dp
    ptable(50) % gyrom_ratio_isotope = 119
    ! Antimony
    ptable(51) % gyrom_ratio = 10.2551487581_dp
    ptable(51) % gyrom_ratio_isotope = 121
    ! Tellurium
    ptable(52) % gyrom_ratio = -13.5454231953_dp
    ptable(52) % gyrom_ratio_isotope = 125
    ! Iodine
    ptable(53) % gyrom_ratio = 8.5777718410_dp
    ptable(53) % gyrom_ratio_isotope = 127
    ! Xenon
    ptable(54) % gyrom_ratio = -11.8603902888_dp
    ptable(54) % gyrom_ratio_isotope = 129
    ! Cesium
    ptable(55) % gyrom_ratio = 5.6233482338_dp
    ptable(55) % gyrom_ratio_isotope = 133
    ! Barium
    ptable(56) % gyrom_ratio = 4.7634278693_dp
    ptable(56) % gyrom_ratio_isotope = 137
    ! Lantanum
    ptable(57) % gyrom_ratio = 6.0611483090_dp
    ptable(57) % gyrom_ratio_isotope = 139
    ! Cerium
    ptable(58) % gyrom_ratio = 0.0_dp
    ptable(58) % gyrom_ratio_isotope = 0
    ! Praseodymium
    ptable(59) % gyrom_ratio = 13.0359039238_dp
    ptable(59) % gyrom_ratio_isotope = 141
    ! Neodymium
    ptable(60) % gyrom_ratio = -2.3188875208_dp
    ptable(60) % gyrom_ratio_isotope = 143
    ! Promethium
    ptable(61) % gyrom_ratio = 5.7502680939_dp
    ptable(61) % gyrom_ratio_isotope = 147
    ! Samarium
    ptable(62) % gyrom_ratio = -1.7745776155_dp
    ptable(62) % gyrom_ratio_isotope = 147
    ! Europium
    ptable(63) % gyrom_ratio = 4.6742215237_dp
    ptable(63) % gyrom_ratio_isotope = 153
    ! Gadolinium
    ptable(64) % gyrom_ratio = -1.7139395822_dp
    ptable(64) % gyrom_ratio_isotope = 157
    ! Terbium
    ptable(65) % gyrom_ratio = 10.2352543902_dp
    ptable(65) % gyrom_ratio_isotope = 159
    ! Dysprosium
    ptable(66) % gyrom_ratio = 2.0515072165_dp
    ptable(66) % gyrom_ratio_isotope = 163
    ! Holmium
    ptable(67) % gyrom_ratio = 9.0877472505_dp
    ptable(67) % gyrom_ratio_isotope = 165
    ! Erbium
    ptable(68) % gyrom_ratio = -1.2279917944_dp
    ptable(68) % gyrom_ratio_isotope = 167
    ! Thulium
    ptable(69) % gyrom_ratio = -3.5300566378_dp
    ptable(69) % gyrom_ratio_isotope = 169
    ! Ytterbium
    ptable(70) % gyrom_ratio = -2.0729931338_dp
    ptable(70) % gyrom_ratio_isotope = 173
    ! Lutetium
    ptable(71) % gyrom_ratio = 4.8625018213_dp
    ptable(71) % gyrom_ratio_isotope = 175
    ! Hafnium
    ptable(72) % gyrom_ratio = 1.7284226820_dp
    ptable(72) % gyrom_ratio_isotope = 177
    ! Tantalum
    ptable(73) % gyrom_ratio = 5.1626680440_dp
    ptable(73) % gyrom_ratio_isotope = 181
    ! Tungsten
    ptable(74) % gyrom_ratio = 1.7956502074_dp
    ptable(74) % gyrom_ratio_isotope = 183
    ! Rhenium
    ptable(75) % gyrom_ratio = 9.8169951998_dp
    ptable(75) % gyrom_ratio_isotope = 187
    ! Osmium
    ptable(76) % gyrom_ratio = 3.3536015524_dp
    ptable(76) % gyrom_ratio_isotope = 189
    ! Iridium
    ptable(77) % gyrom_ratio = 0.8319028875_dp
    ptable(77) % gyrom_ratio_isotope = 193
    ! Platinum
    ptable(78) % gyrom_ratio = 9.2922613524_dp
    ptable(78) % gyrom_ratio_isotope = 195
    ! Gold
    ptable(79) % gyrom_ratio = 0.7528983738_dp
    ptable(79) % gyrom_ratio_isotope = 197
    ! Mercury
    ptable(80) % gyrom_ratio = 7.7123168633_dp
    ptable(80) % gyrom_ratio_isotope = 199
    ! Thallium
    ptable(81) % gyrom_ratio = 24.9748814221_dp
    ptable(81) % gyrom_ratio_isotope = 205
    ! Lead
    ptable(82) % gyrom_ratio = 8.8815779373_dp
    ptable(82) % gyrom_ratio_isotope = 207
    ! Bismuth
    ptable(83) % gyrom_ratio = 6.9630287603_dp
    ptable(83) % gyrom_ratio_isotope = 209
    ! Polonium
    ptable(84) % gyrom_ratio = 11.7774657888_dp
    ptable(84) % gyrom_ratio_isotope = 209
    ! Astatine
    ptable(85) % gyrom_ratio = 0.0_dp
    ptable(85) % gyrom_ratio_isotope = 0
    ! Radon
    ptable(86) % gyrom_ratio = 0.0_dp
    ptable(86) % gyrom_ratio_isotope = 0
    ! Francium
    ptable(87) % gyrom_ratio = 0.0_dp
    ptable(87) % gyrom_ratio_isotope = 0
    ! Radium
    ptable(88) % gyrom_ratio = 0.0_dp
    ptable(88) % gyrom_ratio_isotope = 0
    ! Actinium
    ptable(89) % gyrom_ratio = 5.5704230082_dp
    ptable(89) % gyrom_ratio_isotope = 227
    ! Thorium
    ptable(90) % gyrom_ratio = 0.6366197724_dp
    ptable(90) % gyrom_ratio_isotope = 229
    ! Proctactinium
    ptable(91) % gyrom_ratio = 5.1088736732_dp
    ptable(91) % gyrom_ratio_isotope = 231
    ! Uranium
    ptable(92) % gyrom_ratio = -0.8276057041_dp
    ptable(92) % gyrom_ratio_isotope = 235
    ! Neptunium
    ptable(93) % gyrom_ratio = 4.9338032358_dp
    ptable(93) % gyrom_ratio_isotope = 237
    ! Plutonium
    ptable(94) % gyrom_ratio = 1.5469860469_dp
    ptable(94) % gyrom_ratio_isotope = 239
    ! Americium
    ptable(95) % gyrom_ratio = 2.4509861236_dp
    ptable(95) % gyrom_ratio_isotope = 243
    ! Curium
    ptable(96) % gyrom_ratio = 0.3183098862_dp
    ptable(96) % gyrom_ratio_isotope = 247
    ! No data available
    ptable(97:nelem) % gyrom_ratio = 0.0_dp
    ptable(97:nelem) % gyrom_ratio_isotope = 0

  END SUBROUTINE init_gratio

END MODULE periodic_table
