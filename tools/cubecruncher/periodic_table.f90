!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/periodic_table [1.0] *
!!
!!   NAME
!!     periodic_table
!!
!!   FUNCTION
!!     Periodic Table related data definitions
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE periodic_table
  USE UTILS, ONLY: toupper
  IMPLICIT NONE
  INTEGER, PARAMETER :: dbl=KIND(0.0D0)

  PRIVATE
  PUBLIC :: init_periodic_table, ptable, atom, nelem, label2Z

  TYPE atom
     CHARACTER ( LEN = 2 ) :: symbol
     CHARACTER ( LEN = 14 ) :: name
     INTEGER :: number
     REAL ( dbl ) :: amass ! average mass in formula units
     REAL ( dbl ) :: mass ! mass of most abundant isotope in formula units
     REAL ( dbl ) :: covalent_radius ! in Angstroms
     REAL ( dbl ) :: vdw_radius ! in Angstroms
     INTEGER :: e_conv ( 0:3 )
     REAL ( dbl ) :: eht_param ( 0:3 ) ! in eV
  END TYPE atom

  INTEGER, PARAMETER :: nelem = 106
  TYPE ( atom ) :: ptable ( 0:nelem )

!!*****
!******************************************************************************

CONTAINS

!-------------------------------------------------------------------------------
! gives back the Z value for a given atom label
!-------------------------------------------------------------------------------
   FUNCTION label2Z(label)
     CHARACTER(LEN=2) :: label
     INTEGER          :: label2Z

     CHARACTER(LEN=2) :: mylabel
     INTEGER          :: I
     mylabel=toupper(label)
     label2Z=0
     DO I=1,nelem
        IF (mylabel.EQ.toupper(ptable(I)%symbol)) THEN
           label2Z=I
           EXIT
        ENDIF
     ENDDO
   END FUNCTION label2Z

!******************************************************************************
!!****** periodic_table/init_periodic_table [1.0] *
!!
!!   NAME
!!     init_periodic_table
!!
!!   FUNCTION
!!     Initialization of Periodic Table related data
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE init_periodic_table()

    REAL(dbl), PARAMETER                     :: z = 0.0_dbl

!------------------------------------------------------------------------------
! Dummy

  ptable(0) % symbol = 'X '
  ptable(0) % name = 'Dummy'
  ptable(0) % number = 0
  ptable(0) % amass = z
  ptable(0) % mass = z
  ptable(0) % covalent_radius = z
  ptable(0) % vdw_radius = z
  ptable(0) % e_conv(0:3) = (/ 0, 0, 0, 0 /)
  ptable(0) % eht_param(0:3) = (/ z, z, z, z/)

! Hydrogen
  ptable(1) % symbol = 'H '
  ptable(1) % name = 'Hydrogen'
  ptable(1) % number = 1
  ptable(1) % amass = 1.00797_dbl
  ptable(1) % mass = 1.007825_dbl
  ptable(1) % covalent_radius = 0.32_dbl
  ptable(1) % vdw_radius = 1.2_dbl
  ptable(1) % e_conv(0:3) = (/ 1, 0, 0, 0 /)
  ptable(1) % eht_param(0:3) = (/ -13.60_dbl, z, z, z/)

! Helium
  ptable(2) % symbol = 'He'
  ptable(2) % name = 'Helium'
  ptable(2) % number = 2
  ptable(2) % amass = 4.00260_dbl
  ptable(2) % mass = 4.00260_dbl
  ptable(2) % covalent_radius = 0.9300_dbl
  ptable(2) % vdw_radius = z
  ptable(2) % e_conv(0:3) = (/ 2, 0, 0, 0 /)
  ptable(2) % eht_param(0:3) = (/ -23.40_dbl, z, z, z/)

! Lithium
  ptable(3) % symbol = 'Li'
  ptable(3) % name = 'Lithium'
  ptable(3) % number = 3
  ptable(3) % amass = 6.93900_dbl
  ptable(3) % mass = 7.01600_dbl
  ptable(3) % covalent_radius = 1.2300_dbl
  ptable(3) % vdw_radius = z
  ptable(3) % e_conv(0:3) = (/ 3, 0, 0, 0 /)
  ptable(3) % eht_param(0:3) = (/ -5.40_dbl, -3.50_dbl, z, z/)

! Beryllium
  ptable(4) % symbol = 'Be'
  ptable(4) % name = 'Beryllium'
  ptable(4) % number = 4
  ptable(4) % amass = 9.01220_dbl
  ptable(4) % mass = 9.01218_dbl
  ptable(4) % covalent_radius = 0.9000_dbl
  ptable(4) % vdw_radius = z
  ptable(4) % e_conv(0:3) = (/ 4, 0, 0, 0 /)
  ptable(4) % eht_param(0:3) = (/ -10.00_dbl, -6.00_dbl, z, z/)

! Boron
  ptable(5) % symbol = 'B '
  ptable(5) % name = 'Boron'
  ptable(5) % number = 5
  ptable(5) % amass = 10.81100_dbl
  ptable(5) % mass = 11.00931_dbl
  ptable(5) % covalent_radius = 0.8200_dbl
  ptable(5) % vdw_radius = z
  ptable(5) % e_conv(0:3) = (/ 4, 1, 0, 0 /)
  ptable(5) % eht_param(0:3) = (/ -15.20_dbl, -8.50_dbl, z, z/)

! Carbon
  ptable(6) % symbol = 'C '
  ptable(6) % name = 'Carbon'
  ptable(6) % number = 6
  ptable(6) % amass = 12.01115_dbl
  ptable(6) % mass = 12.0000_dbl
  ptable(6) % covalent_radius = 0.7700_dbl
  ptable(6) % vdw_radius = 1.7_dbl           ! obtained from www.webelement.com
  ptable(6) % e_conv(0:3) = (/ 4, 2, 0, 0 /)
  ptable(6) % eht_param(0:3) = (/ -21.40_dbl, -11.40_dbl, z, z/)

! Nitrogen
  ptable(7) % symbol = 'N '
  ptable(7) % name = 'Nitrogen'
  ptable(7) % number = 7
  ptable(7) % amass = 14.00670_dbl
  ptable(7) % mass = 14.00307_dbl
  ptable(7) % covalent_radius = 0.7500_dbl
  ptable(7) % vdw_radius = 1.5_dbl
  ptable(7) % e_conv(0:3) = (/ 4, 3, 0, 0 /)
  ptable(7) % eht_param(0:3) = (/ -26.00_dbl, -13.40_dbl, z, z/)

! Oxygen
  ptable(8) % symbol = 'O '
  ptable(8) % name = 'Oxygen'
  ptable(8) % number = 8
  ptable(8) % amass = 15.99940_dbl
  ptable(8) % mass = 15.99491_dbl
  ptable(8) % covalent_radius = 0.7300_dbl
  ptable(8) % vdw_radius = 1.40_dbl
  ptable(8) % e_conv(0:3) = (/ 4, 4, 0, 0 /)
  ptable(8) % eht_param(0:3) = (/ -32.30_dbl, -14.80_dbl, z, z/)

! Fluorine
  ptable(9) % symbol = 'F '
  ptable(9) % name = 'Fluorine'
  ptable(9) % number = 9
  ptable(9) % amass = 18.99840_dbl
  ptable(9) % mass = 18.99840_dbl
  ptable(9) % covalent_radius = 0.7200_dbl
  ptable(9) % vdw_radius = 1.35_dbl
  ptable(9) % e_conv(0:3) = (/ 4, 5, 0, 0 /)
  ptable(9) % eht_param(0:3) = (/ -40.00_dbl, -18.10_dbl, z, z/)

! Neon
  ptable(10) % symbol = 'Ne'
  ptable(10) % name = 'Neon'
  ptable(10) % number = 10
  ptable(10) % amass = 20.18300_dbl
  ptable(10) % mass = 19.99244_dbl
  ptable(10) % covalent_radius = 0.7100_dbl
  ptable(10) % vdw_radius = z
  ptable(10) % e_conv(0:3) = (/ 4, 6, 0, 0 /)
  ptable(10) % eht_param(0:3) = (/ -43.20_dbl, -20.00_dbl, z, z/)

! Sodium
  ptable(11) % symbol = 'Na'
  ptable(11) % name = 'Sodium'
  ptable(11) % number = 11
  ptable(11) % amass = 22.98980_dbl
  ptable(11) % mass = 22.9898_dbl
  ptable(11) % covalent_radius = 1.5400_dbl
  ptable(11) % vdw_radius = z
  ptable(11) % e_conv(0:3) = (/ 5, 6, 0, 0 /)
  ptable(11) % eht_param(0:3) = (/ -5.10_dbl, -3.00_dbl, z, z/)

! Magnesium
  ptable(12) % symbol = 'Mg'
  ptable(12) % name = 'Magnesium'
  ptable(12) % number = 12
  ptable(12) % amass = 24.31200_dbl
  ptable(12) % mass = 23.98504_dbl
  ptable(12) % covalent_radius = 1.3600_dbl
  ptable(12) % vdw_radius = z
  ptable(12) % e_conv(0:3) = (/ 6, 6, 0, 0 /)
  ptable(12) % eht_param(0:3) = (/ -9.00_dbl, -4.50_dbl, z, z/)

! Aluminium
  ptable(13) % symbol = 'Al'
  ptable(13) % name = 'Aluminium'
  ptable(13) % number = 13
  ptable(13) % amass = 26.98153_dbl
  ptable(13) % mass = 26.98153_dbl
  ptable(13) % covalent_radius = 1.1800_dbl
  ptable(13) % vdw_radius = z
  ptable(13) % e_conv(0:3) = (/ 6, 7, 0, 0 /)
  ptable(13) % eht_param(0:3) = (/ -12.30_dbl, -6.50_dbl, z, z/)

! Silicon
  ptable(14) % symbol = 'Si'
  ptable(14) % name = 'Silicon'
  ptable(14) % number = 14
  ptable(14) % amass = 28.08600_dbl
  ptable(14) % mass = 27.97693_dbl
  ptable(14) % covalent_radius = 1.1100_dbl
  ptable(14) % vdw_radius = z
  ptable(14) % e_conv(0:3) = (/ 6, 8, 0, 0 /)
  ptable(14) % eht_param(0:3) = (/ -17.30_dbl, -9.20_dbl, z, z/)

! Phosphorus
  ptable(15) % symbol = 'P '
  ptable(15) % name = 'Phosphorus'
  ptable(15) % number = 15
  ptable(15) % amass = 30.97380_dbl
  ptable(15) % mass = 30.97376_dbl
  ptable(15) % covalent_radius = 1.0600_dbl
  ptable(15) % vdw_radius = 1.9_dbl
  ptable(15) % e_conv(0:3) = (/ 6, 9, 0, 0 /)
  ptable(15) % eht_param(0:3) = (/ -18.60_dbl, -14.00_dbl, z, z/)

! Sulfur
  ptable(16) % symbol = 'S '
  ptable(16) % name = 'Sulfur'
  ptable(16) % number = 16
  ptable(16) % amass = 32.06400_dbl
  ptable(16) % mass = 31.97207_dbl
  ptable(16) % covalent_radius = 1.0200_dbl
  ptable(16) % vdw_radius = 1.85_dbl
  ptable(16) % e_conv(0:3) = (/ 6, 10, 0, 0 /)
  ptable(16) % eht_param(0:3) = (/ -20.00_dbl, -11.00_dbl, z, z/)

! Chlorine
  ptable(17) % symbol = 'Cl'
  ptable(17) % name = 'Chlorine'
  ptable(17) % number = 17
  ptable(17) % amass = 35.45300_dbl
  ptable(17) % mass = 34.96885_dbl
  ptable(17) % covalent_radius = 0.9900_dbl
  ptable(17) % vdw_radius = 1.80_dbl
  ptable(17) % e_conv(0:3) = (/ 6, 11, 0, 0 /)
  ptable(17) % eht_param(0:3) = (/ -26.30_dbl, -14.20_dbl, z, z/)

! Argon
  ptable(18) % symbol = 'Ar'
  ptable(18) % name = 'Argon'
  ptable(18) % number = 18
  ptable(18) % amass = 39.94800_dbl
  ptable(18) % mass = 39.94800_dbl
  ptable(18) % covalent_radius = 0.9800_dbl
  ptable(18) % vdw_radius = 3.83_dbl
  ptable(18) % e_conv(0:3) = (/ 6, 12, 0, 0 /)
  ptable(18) % eht_param(0:3) = (/ z, z, z, z/)

! Potassium
  ptable(19) % symbol = 'K '
  ptable(19) % name = 'Potassium'
  ptable(19) % number = 19
  ptable(19) % amass = 39.10200_dbl
  ptable(19) % mass = 38.96371_dbl
  ptable(19) % covalent_radius = 2.0300_dbl
  ptable(19) % vdw_radius = z
  ptable(19) % e_conv(0:3) = (/ 7, 12, 0, 0 /)
  ptable(19) % eht_param(0:3) = (/ -4.34_dbl, -2.73_dbl, z, z/)

! Calcium
  ptable(20) % symbol = 'Ca'
  ptable(20) % name = 'Calcium'
  ptable(20) % number = 20
  ptable(20) % amass = 40.08000_dbl
  ptable(20) % mass = 39.96259_dbl
  ptable(20) % covalent_radius = 1.7400_dbl
  ptable(20) % vdw_radius = z
  ptable(20) % e_conv(0:3) = (/ 8, 12, 0, 0 /)
  ptable(20) % eht_param(0:3) = (/ -7.00_dbl, -4.00_dbl, z, z/)

! Scandium
  ptable(21) % symbol = 'Sc'
  ptable(21) % name = 'Scandium'
  ptable(21) % number = 21
  ptable(21) % amass = 44.95600_dbl
  ptable(21) % mass = 44.95592_dbl
  ptable(21) % covalent_radius = 1.4400_dbl
  ptable(21) % vdw_radius = z
  ptable(21) % e_conv(0:3) = (/ 8, 12, 1, 0 /)
  ptable(21) % eht_param(0:3) = (/ -8.87_dbl, -2.75_dbl, -8.51_dbl, z/)

! Titanium
  ptable(22) % symbol = 'Ti'
  ptable(22) % name = 'Titanium'
  ptable(22) % number = 22
  ptable(22) % amass = 47.90000_dbl
  ptable(22) % mass = 48.00000_dbl
  ptable(22) % covalent_radius = 1.3200_dbl
  ptable(22) % vdw_radius = z
  ptable(22) % e_conv(0:3) = (/ 8, 12, 2, 0 /)
  ptable(22) % eht_param(0:3) = (/ -8.97_dbl, -5.44_dbl, -10.81_dbl, z/)

! Vanadium
  ptable(23) % symbol = 'V '
  ptable(23) % name = 'Vanadium'
  ptable(23) % number = 23
  ptable(23) % amass = 50.94200_dbl
  ptable(23) % mass = 50.94400_dbl
  ptable(23) % covalent_radius = 1.2200_dbl
  ptable(23) % vdw_radius = z
  ptable(23) % e_conv(0:3) = (/ 8, 12, 3, 0 /)
  ptable(23) % eht_param(0:3) = (/ -8.81_dbl, -5.52_dbl, -11.00_dbl, z/)

! Chromium
  ptable(24) % symbol = 'Cr'
  ptable(24) % name = 'Chromium'
  ptable(24) % number = 24
  ptable(24) % amass = 51.99600_dbl
  ptable(24) % mass = 51.94050_dbl
  ptable(24) % covalent_radius = 1.1800_dbl
  ptable(24) % vdw_radius = z
  ptable(24) % e_conv(0:3) = (/ 7, 12, 5, 0 /)
  ptable(24) % eht_param(0:3) = (/ -8.66_dbl, -5.24_dbl, -11.22_dbl, z/)

! Manganese
  ptable(25) % symbol = 'Mn'
  ptable(25) % name = 'Manganese'
  ptable(25) % number = 25
  ptable(25) % amass = 54.93800_dbl
  ptable(25) % mass = 54.93810_dbl
  ptable(25) % covalent_radius = 1.1700_dbl
  ptable(25) % vdw_radius = z
  ptable(25) % e_conv(0:3) = (/ 8, 12, 5, 0 /)
  ptable(25) % eht_param(0:3) = (/ -9.75_dbl, -5.89_dbl, -11.67_dbl, z/)

! Iron
  ptable(26) % symbol = 'Fe'
  ptable(26) % name = 'Iron'
  ptable(26) % number = 26
  ptable(26) % amass = 55.84700_dbl
  ptable(26) % mass = 55.93490_dbl
  ptable(26) % covalent_radius = 1.1700_dbl
  ptable(26) % vdw_radius = z
  ptable(26) % e_conv(0:3) = (/ 8, 12, 6, 0 /)
  ptable(26) % eht_param(0:3) = (/ -9.10_dbl, -5.32_dbl, -12.60_dbl, z/)

! Cobalt
  ptable(27) % symbol = 'Co'
  ptable(27) % name = 'Cobalt'
  ptable(27) % number = 27
  ptable(27) % amass = 58.93300_dbl
  ptable(27) % mass = 58.93320_dbl
  ptable(27) % covalent_radius = 1.1600_dbl
  ptable(27) % vdw_radius = z
  ptable(27) % e_conv(0:3) = (/ 8, 12, 7, 0 /)
  ptable(27) % eht_param(0:3) = (/ -9.21_dbl, -5.29_dbl, -13.18_dbl, z/)

! Nickel
  ptable(28) % symbol = 'Ni'
  ptable(28) % name = 'Nickel'
  ptable(28) % number = 28
  ptable(28) % amass = 58.71000_dbl
  ptable(28) % mass = 57.93530_dbl
  ptable(28) % covalent_radius = 1.1500_dbl
  ptable(28) % vdw_radius = z
  ptable(28) % e_conv(0:3) = (/ 8, 12, 8, 0 /)
  ptable(28) % eht_param(0:3) = (/ -9.17_dbl, -5.15_dbl, -13.49_dbl, z/)

! Copper
  ptable(29) % symbol = 'Cu'
  ptable(29) % name = 'Copper'
  ptable(29) % number = 29
  ptable(29) % amass = 63.54000_dbl
  ptable(29) % mass = 63.20000_dbl
  ptable(29) % covalent_radius = 1.1700_dbl
  ptable(29) % vdw_radius = z
  ptable(29) % e_conv(0:3) = (/ 7, 12, 10, 0 /)
  ptable(29) % eht_param(0:3) = (/ -11.40_dbl, -6.06_dbl, -14.00_dbl, z/)

! Zinc
  ptable(30) % symbol = 'Zn'
  ptable(30) % name = 'Zinc'
  ptable(30) % number = 30
  ptable(30) % amass = 65.37000_dbl
  ptable(30) % mass = 63.92910_dbl
  ptable(30) % covalent_radius = 1.2500_dbl
  ptable(30) % vdw_radius = z
  ptable(30) % e_conv(0:3) = (/ 8, 12, 10, 0 /)
  ptable(30) % eht_param(0:3) = (/ -12.41_dbl, -6.53_dbl, z, z/)

! Gallium
  ptable(31) % symbol = 'Ga'
  ptable(31) % name = 'Gallium'
  ptable(31) % number = 31
  ptable(31) % amass = 69.72000_dbl
  ptable(31) % mass = 68.92570_dbl
  ptable(31) % covalent_radius = 1.2600_dbl
  ptable(31) % vdw_radius = z
  ptable(31) % e_conv(0:3) = (/ 8, 13, 10, 0 /)
  ptable(31) % eht_param(0:3) = (/ -14.58_dbl, -6.75_dbl, z, z/)

! Germanium
  ptable(32) % symbol = 'Ge'
  ptable(32) % name = 'Germanium'
  ptable(32) % number = 32
  ptable(32) % amass = 72.59000_dbl
  ptable(32) % mass = 73.92140_dbl
  ptable(32) % covalent_radius = 1.2200_dbl
  ptable(32) % vdw_radius = z
  ptable(32) % e_conv(0:3) = (/ 8, 14, 10, 0 /)
  ptable(32) % eht_param(0:3) = (/ -16.00_dbl, -9.00_dbl, z, z/)

! Arsenic
  ptable(33) % symbol = 'As'
  ptable(33) % name = 'Arsenic'
  ptable(33) % number = 33
  ptable(33) % amass = 74.92200_dbl
  ptable(33) % mass = 74.92160_dbl
  ptable(33) % covalent_radius = 1.2000_dbl
  ptable(33) % vdw_radius = 2.0_dbl
  ptable(33) % e_conv(0:3) = (/ 8, 15, 10, 0 /)
  ptable(33) % eht_param(0:3) = (/ -16.22_dbl, -12.16_dbl, z, z/)

! Selenium
  ptable(34) % symbol = 'Se'
  ptable(34) % name = 'Selenium'
  ptable(34) % number = 34
  ptable(34) % amass = 78.96000_dbl
  ptable(34) % mass = 79.91650_dbl
  ptable(34) % covalent_radius = 1.1600_dbl
  ptable(34) % vdw_radius = 2.00_dbl
  ptable(34) % e_conv(0:3) = (/ 8, 16, 10, 0 /)
  ptable(34) % eht_param(0:3) = (/ -20.50_dbl, -14.40_dbl, z, z/)

! Bromine
  ptable(35) % symbol = 'Br'
  ptable(35) % name = 'Bromine'
  ptable(35) % number = 35
  ptable(35) % amass = 79.90900_dbl
  ptable(35) % mass = 78.91830_dbl
  ptable(35) % covalent_radius = 1.1400_dbl
  ptable(35) % vdw_radius = 1.95_dbl
  ptable(35) % e_conv(0:3) = (/ 8, 17, 10, 0 /)
  ptable(35) % eht_param(0:3) = (/ -22.07_dbl, -13.10_dbl, z, z/)

! Krypton
  ptable(36) % symbol = 'Kr'
  ptable(36) % name = 'Krypton'
  ptable(36) % number = 36
  ptable(36) % amass = 83.80000_dbl
  ptable(36) % mass = 84.00000_dbl
  ptable(36) % covalent_radius = 1.1200_dbl
  ptable(36) % vdw_radius = z
  ptable(36) % e_conv(0:3) = (/ 8, 18, 10, 0 /)
  ptable(36) % eht_param(0:3) = (/ z, z, z, z/)

! Rubidium
  ptable(37) % symbol = 'Rb'
  ptable(37) % name = 'Rubidium'
  ptable(37) % number = 37
  ptable(37) % amass = 85.47000_dbl
  ptable(37) % mass = 84.91170_dbl
  ptable(37) % covalent_radius = 2.1600_dbl
  ptable(37) % vdw_radius = z
  ptable(37) % e_conv(0:3) = (/ 9, 18, 10, 0 /)
  ptable(37) % eht_param(0:3) = (/ -4.18_dbl, -2.60_dbl, z, z/)

! Strontium
  ptable(38) % symbol = 'Sr'
  ptable(38) % name = 'Strontium'
  ptable(38) % number = 38
  ptable(38) % amass = 87.62000_dbl
  ptable(38) % mass = 87.90560_dbl
  ptable(38) % covalent_radius = 1.9100_dbl
  ptable(38) % vdw_radius = z
  ptable(38) % e_conv(0:3) = (/ 10, 18, 10, 0 /)
  ptable(38) % eht_param(0:3) = (/ -6.62_dbl, -3.92_dbl, z, z/)

! Yttrium
  ptable(39) % symbol = 'Y '
  ptable(39) % name = 'Yttrium'
  ptable(39) % number = 39
  ptable(39) % amass = 88.90500_dbl
  ptable(39) % mass = 88.90590_dbl
  ptable(39) % covalent_radius = 1.6200_dbl
  ptable(39) % vdw_radius = z
  ptable(39) % e_conv(0:3) = (/ 10, 18, 11, 0 /)
  ptable(39) % eht_param(0:3) = (/ z, z, z, z/)

! Zirconium
  ptable(40) % symbol = 'Zr'
  ptable(40) % name = 'Zirconium'
  ptable(40) % number = 40
  ptable(40) % amass = 91.22000_dbl
  ptable(40) % mass = 89.90430_dbl
  ptable(40) % covalent_radius = 1.4500_dbl
  ptable(40) % vdw_radius = z
  ptable(40) % e_conv(0:3) = (/ 10, 18, 12, 0 /)
  ptable(40) % eht_param(0:3) = (/ -8.00_dbl, -5.40_dbl, -10.20_dbl, z/)

! Niobium
  ptable(41) % symbol = 'Nb'
  ptable(41) % name = 'Niobium'
  ptable(41) % number = 41
  ptable(41) % amass = 92.90600_dbl
  ptable(41) % mass = 92.90600_dbl
  ptable(41) % covalent_radius = 1.3400_dbl
  ptable(41) % vdw_radius = z
  ptable(41) % e_conv(0:3) = (/ 9, 18, 14, 0 /)
  ptable(41) % eht_param(0:3) = (/ -10.10_dbl, -6.86_dbl, -12.10_dbl, z/)

! Molybdenum
  ptable(42) % symbol = 'Mo'
  ptable(42) % name = 'Molybdenum'
  ptable(42) % number = 42
  ptable(42) % amass = 95.94000_dbl
  ptable(42) % mass = 97.90550_dbl
  ptable(42) % covalent_radius = 1.3000_dbl
  ptable(42) % vdw_radius = z
  ptable(42) % e_conv(0:3) = (/ 9, 18, 15, 0 /)
  ptable(42) % eht_param(0:3) = (/ -8.34_dbl, -5.25_dbl, -10.50_dbl, z/)

! Technetium
  ptable(43) % symbol = 'Tc'
  ptable(43) % name = 'Technetium'
  ptable(43) % number = 43
  ptable(43) % amass = 98.90600_dbl
  ptable(43) % mass = 98.90600_dbl
  ptable(43) % covalent_radius = 1.2700_dbl
  ptable(43) % vdw_radius = z
  ptable(43) % e_conv(0:3) = (/ 9, 18, 16, 0 /)
  ptable(43) % eht_param(0:3) = (/ -10.07_dbl, -5.40_dbl, -12.82_dbl, z/)

! Ruthenium
  ptable(44) % symbol = 'Ru'
  ptable(44) % name = 'Ruthenium'
  ptable(44) % number = 44
  ptable(44) % amass = 101.07000_dbl
  ptable(44) % mass = 101.90370_dbl
  ptable(44) % covalent_radius = 1.2500_dbl
  ptable(44) % vdw_radius = z
  ptable(44) % e_conv(0:3) = (/ 9, 18, 17, 0 /)
  ptable(44) % eht_param(0:3) = (/ -10.40_dbl, -6.87_dbl, -14.90_dbl, z/)

! Rhodium
  ptable(45) % symbol = 'Rh'
  ptable(45) % name = 'Rhodium'
  ptable(45) % number = 45
  ptable(45) % amass = 102.90500_dbl
  ptable(45) % mass = 102.90480_dbl
  ptable(45) % covalent_radius = 1.2500_dbl
  ptable(45) % vdw_radius = z
  ptable(45) % e_conv(0:3) = (/ 9, 18, 18, 0 /)
  ptable(45) % eht_param(0:3) = (/ -8.09_dbl, -4.57_dbl, -12.50_dbl, z/)

! Palladium
  ptable(46) % symbol = 'Pd'
  ptable(46) % name = 'Palladium'
  ptable(46) % number = 46
  ptable(46) % amass = 106.40000_dbl
  ptable(46) % mass = 105.90320_dbl
  ptable(46) % covalent_radius = 1.2800_dbl
  ptable(46) % vdw_radius = z
  ptable(46) % e_conv(0:3) = (/ 8, 18, 20, 0 /)
  ptable(46) % eht_param(0:3) = (/ -7.32_dbl, -3.75_dbl, -12.02_dbl, z/)

! Silver
  ptable(47) % symbol = 'Ag'
  ptable(47) % name = 'Silver'
  ptable(47) % number = 47
  ptable(47) % amass = 107.87000_dbl
  ptable(47) % mass = 106.90509_dbl
  ptable(47) % covalent_radius = 1.3400_dbl
  ptable(47) % vdw_radius = z
  ptable(47) % e_conv(0:3) = (/ 9, 18, 20, 0 /)
  ptable(47) % eht_param(0:3) = (/ z, z, z, z/)

! Cadmium
  ptable(48) % symbol = 'Cd'
  ptable(48) % name = 'Cadmium'
  ptable(48) % number = 48
  ptable(48) % amass = 112.40000_dbl
  ptable(48) % mass = 113.90360_dbl
  ptable(48) % covalent_radius = 1.4800_dbl
  ptable(48) % vdw_radius = z
  ptable(48) % e_conv(0:3) = (/ 10, 18, 20, 0 /)
  ptable(48) % eht_param(0:3) = (/ z, z, z, z/)

! Indium
  ptable(49) % symbol = 'In'
  ptable(49) % name = 'Indium'
  ptable(49) % number = 49
  ptable(49) % amass = 114.82000_dbl
  ptable(49) % mass = 114.90410_dbl
  ptable(49) % covalent_radius = 1.4400_dbl
  ptable(49) % vdw_radius = z
  ptable(49) % e_conv(0:3) = (/ 10, 19, 20, 0 /)
  ptable(49) % eht_param(0:3) = (/ -12.60_dbl, -6.19_dbl, z, z/)

! Tin
  ptable(50) % symbol = 'Sn'
  ptable(50) % name = 'Tin'
  ptable(50) % number = 50
  ptable(50) % amass = 118.69000_dbl
  ptable(50) % mass = 120.00000_dbl
  ptable(50) % covalent_radius = 1.4100_dbl
  ptable(50) % vdw_radius = z
  ptable(50) % e_conv(0:3) = (/ 10, 20, 20, 0 /)
  ptable(50) % eht_param(0:3) = (/ -16.16_dbl, -8.32_dbl, z, z/)

! Antimony
  ptable(51) % symbol = 'Sb'
  ptable(51) % name = 'Antimony'
  ptable(51) % number = 51
  ptable(51) % amass = 121.75000_dbl
  ptable(51) % mass = 120.90380_dbl
  ptable(51) % covalent_radius = 1.4000_dbl
  ptable(51) % vdw_radius = 2.2_dbl
  ptable(51) % e_conv(0:3) = (/ 10, 21, 20, 0 /)
  ptable(51) % eht_param(0:3) = (/ -18.80_dbl, -11.70_dbl, z, z/)

! Tellurium
  ptable(52) % symbol = 'Te'
  ptable(52) % name = 'Tellurium'
  ptable(52) % number = 52
  ptable(52) % amass = 127.60000_dbl
  ptable(52) % mass = 129.90670_dbl
  ptable(52) % covalent_radius = 1.3600_dbl
  ptable(52) % vdw_radius = 2.20_dbl
  ptable(52) % e_conv(0:3) = (/ 10, 22, 20, 0 /)
  ptable(52) % eht_param(0:3) = (/ -20.80_dbl, -13.20_dbl, z, z/)

! Iodine
  ptable(53) % symbol = 'I '
  ptable(53) % name = 'Iodine'
  ptable(53) % number = 53
  ptable(53) % amass = 126.90440_dbl
  ptable(53) % mass = 126.90440_dbl
  ptable(53) % covalent_radius = 1.3300_dbl
  ptable(53) % vdw_radius = 2.15_dbl
  ptable(53) % e_conv(0:3) = (/ 10, 23, 20, 0 /)
  ptable(53) % eht_param(0:3) = (/ -18.00_dbl, -12.70_dbl, z, z/)

! Xenon
  ptable(54) % symbol = 'Xe'
  ptable(54) % name = 'Xenon'
  ptable(54) % number = 54
  ptable(54) % amass = 131.30000_dbl
  ptable(54) % mass = 131.90420_dbl
  ptable(54) % covalent_radius = 1.3100_dbl
  ptable(54) % vdw_radius = z
  ptable(54) % e_conv(0:3) = (/ 10, 24, 20, 0 /)
  ptable(54) % eht_param(0:3) = (/ z, z, z, z/)

! Cesium
  ptable(55) % symbol = 'Cs'
  ptable(55) % name = 'Cesium'
  ptable(55) % number = 55
  ptable(55) % amass = 132.90500_dbl
  ptable(55) % mass = 132.90510_dbl
  ptable(55) % covalent_radius = 2.3500_dbl
  ptable(55) % vdw_radius = z
  ptable(55) % e_conv(0:3) = (/ 11, 24, 20, 0 /)
  ptable(55) % eht_param(0:3) = (/ -3.88_dbl, -2.49_dbl, z, z/)

! Barium
  ptable(56) % symbol = 'Ba'
  ptable(56) % name = 'Barium'
  ptable(56) % number = 56
  ptable(56) % amass = 137.34000_dbl
  ptable(56) % mass = 137.90500_dbl
  ptable(56) % covalent_radius = 1.9800_dbl
  ptable(56) % vdw_radius = z
  ptable(56) % e_conv(0:3) = (/ 12, 24, 20, 0 /)
  ptable(56) % eht_param(0:3) = (/ z, z, z, z/)

! Lantanum
  ptable(57) % symbol = 'La'
  ptable(57) % name = 'Lantanum'
  ptable(57) % number = 57
  ptable(57) % amass = 138.91000_dbl
  ptable(57) % mass = 138.90610_dbl
  ptable(57) % covalent_radius = 1.6900_dbl
  ptable(57) % vdw_radius = z
  ptable(57) % e_conv(0:3) = (/ 12, 24, 21, 0 /)
  ptable(57) % eht_param(0:3) = (/ -7.67_dbl, -5.01_dbl, -8.21_dbl, z/)

! Cerium
  ptable(58) % symbol = 'Ce'
  ptable(58) % name = 'Cerium'
  ptable(58) % number = 58
  ptable(58) % amass = 140.12000_dbl
  ptable(58) % mass = 139.90530_dbl
  ptable(58) % covalent_radius = 1.6500_dbl
  ptable(58) % vdw_radius = z
  ptable(58) % e_conv(0:3) = (/ 12, 24, 20, 2 /)
  ptable(58) % eht_param(0:3) = (/ z, z, z, z/)

! Praseodymium
  ptable(59) % symbol = 'Pr'
  ptable(59) % name = 'Praseodymium'
  ptable(59) % number = 59
  ptable(59) % amass = 140.90700_dbl
  ptable(59) % mass = 140.90740_dbl
  ptable(59) % covalent_radius = 1.6500_dbl
  ptable(59) % vdw_radius = z
  ptable(59) % e_conv(0:3) = (/ 12, 24, 20, 3 /)
  ptable(59) % eht_param(0:3) = (/ z, z, z, z/)

! Neodymium
  ptable(60) % symbol = 'Nd'
  ptable(60) % name = 'Neodymium'
  ptable(60) % number = 60
  ptable(60) % amass = 144.24000_dbl
  ptable(60) % mass = 141.90750_dbl
  ptable(60) % covalent_radius = 1.6400_dbl
  ptable(60) % vdw_radius = z
  ptable(60) % e_conv(0:3) = (/ 12, 24, 20, 4 /)
  ptable(60) % eht_param(0:3) = (/ z, z, z, z/)

! Promethium
  ptable(61) % symbol = 'Pm'
  ptable(61) % name = 'Promethium'
  ptable(61) % number = 61
  ptable(61) % amass = 144.91300_dbl
  ptable(61) % mass = 144.91300_dbl
  ptable(61) % covalent_radius = 1.6300_dbl
  ptable(61) % vdw_radius = z
  ptable(61) % e_conv(0:3) = (/ 12, 24, 20, 5 /)
  ptable(61) % eht_param(0:3) = (/ z, z, z, z/)

! Samarium
  ptable(62) % symbol = 'Sm'
  ptable(62) % name = 'Samarium'
  ptable(62) % number = 62
  ptable(62) % amass = 150.35000_dbl
  ptable(62) % mass = 151.91950_dbl
  ptable(62) % covalent_radius = 1.6200_dbl
  ptable(62) % vdw_radius = z
  ptable(62) % e_conv(0:3) = (/ 12, 24, 20, 6 /)
  ptable(62) % eht_param(0:3) = (/ -4.86_dbl, -4.86_dbl, -6.06_dbl, &
          -11.28_dbl/)

! Europium
  ptable(63) % symbol = 'Eu'
  ptable(63) % name = 'Europium'
  ptable(63) % number = 63
  ptable(63) % amass = 151.96000_dbl
  ptable(63) % mass = 152.92090_dbl
  ptable(63) % covalent_radius = 1.8500_dbl
  ptable(63) % vdw_radius = z
  ptable(63) % e_conv(0:3) = (/ 12, 24, 20, 7 /)
  ptable(63) % eht_param(0:3) = (/ z, z, z, z/)

! Gadolinium
  ptable(64) % symbol = 'Gd'
  ptable(64) % name = 'Gadolinium'
  ptable(64) % number = 64
  ptable(64) % amass = 157.25000_dbl
  ptable(64) % mass = 157.92410_dbl
  ptable(64) % covalent_radius = 1.6100_dbl
  ptable(64) % vdw_radius = z
  ptable(64) % e_conv(0:3) = (/ 12, 24, 21, 7 /)
  ptable(64) % eht_param(0:3) = (/ z, z, z, z/)

! Terbium
  ptable(65) % symbol = 'Tb'
  ptable(65) % name = 'Terbium'
  ptable(65) % number = 65
  ptable(65) % amass = 158.92400_dbl
  ptable(65) % mass = 158.92500_dbl
  ptable(65) % covalent_radius = 1.5900_dbl
  ptable(65) % vdw_radius = z
  ptable(65) % e_conv(0:3) = (/ 12, 24, 20, 9 /)
  ptable(65) % eht_param(0:3) = (/ z, z, z, z/)

! Dysprosium
  ptable(66) % symbol = 'Dy'
  ptable(66) % name = 'Dysprosium'
  ptable(66) % number = 66
  ptable(66) % amass = 162.50000_dbl
  ptable(66) % mass = 163.92880_dbl
  ptable(66) % covalent_radius = 1.5900_dbl
  ptable(66) % vdw_radius = z
  ptable(66) % e_conv(0:3) = (/ 12, 24, 20, 10 /)
  ptable(66) % eht_param(0:3) = (/ z, z, z, z/)

! Holmium
  ptable(67) % symbol = 'Ho'
  ptable(67) % name = 'Holmium'
  ptable(67) % number = 67
  ptable(67) % amass = 164.93000_dbl
  ptable(67) % mass = 164.93000_dbl
  ptable(67) % covalent_radius = 1.5800_dbl
  ptable(67) % vdw_radius = z
  ptable(67) % e_conv(0:3) = (/ 12, 24, 20, 11 /)
  ptable(67) % eht_param(0:3) = (/ z, z, z, z/)

! Erbium
  ptable(68) % symbol = 'Er'
  ptable(68) % name = 'Erbium'
  ptable(68) % number = 68
  ptable(68) % amass = 167.26000_dbl
  ptable(68) % mass = 165.93040_dbl
  ptable(68) % covalent_radius = 1.5700_dbl
  ptable(68) % vdw_radius = z
  ptable(68) % e_conv(0:3) = (/ 12, 24, 20, 12 /)
  ptable(68) % eht_param(0:3) = (/ z, z, z, z/)

! Thulium
  ptable(69) % symbol = 'Tm'
  ptable(69) % name = 'Thulium'
  ptable(69) % number = 69
  ptable(69) % amass = 168.93400_dbl
  ptable(69) % mass = 168.93440_dbl
  ptable(69) % covalent_radius = 1.5600_dbl
  ptable(69) % vdw_radius = z
  ptable(69) % e_conv(0:3) = (/ 12, 24, 20, 13 /)
  ptable(69) % eht_param(0:3) = (/ z, z, z, z/)

! Ytterbium
  ptable(70) % symbol = 'Yb'
  ptable(70) % name = 'Ytterbium'
  ptable(70) % number = 70
  ptable(70) % amass = 173.04000_dbl
  ptable(70) % mass = 173.93900_dbl
  ptable(70) % covalent_radius = 1.5600_dbl
  ptable(70) % vdw_radius = z
  ptable(70) % e_conv(0:3) = (/ 12, 24, 20, 14 /)
  ptable(70) % eht_param(0:3) = (/ -5.35_dbl, -5.35_dbl, -5.21_dbl, &
          -13.86_dbl/)

! Lutetium
  ptable(71) % symbol = 'Lu'
  ptable(71) % name = 'Lutetium'
  ptable(71) % number = 71
  ptable(71) % amass = 174.97000_dbl
  ptable(71) % mass = 174.94090_dbl
  ptable(71) % covalent_radius = 1.5600_dbl
  ptable(71) % vdw_radius = z
  ptable(71) % e_conv(0:3) = (/ 12, 24, 21, 14 /)
  ptable(71) % eht_param(0:3) = (/ -6.05_dbl, -6.05_dbl, -5.12_dbl, &
          -22.40_dbl/)

! Hafnium
  ptable(72) % symbol = 'Hf'
  ptable(72) % name = 'Hafnium'
  ptable(72) % number = 72
  ptable(72) % amass = 178.49000_dbl
  ptable(72) % mass = 179.94680_dbl
  ptable(72) % covalent_radius = 1.4400_dbl
  ptable(72) % vdw_radius = z
  ptable(72) % e_conv(0:3) = (/ 12, 24, 22, 14 /)
  ptable(72) % eht_param(0:3) = (/ z, z, z, z/)

! Tantalum
  ptable(73) % symbol = 'Ta'
  ptable(73) % name = 'Tantalum'
  ptable(73) % number = 73
  ptable(73) % amass = 180.94800_dbl
  ptable(73) % mass = 180.94800_dbl
  ptable(73) % covalent_radius = 1.3400_dbl
  ptable(73) % vdw_radius = z
  ptable(73) % e_conv(0:3) = (/ 12, 24, 23, 14 /)
  ptable(73) % eht_param(0:3) = (/ -10.10_dbl, -6.86_dbl, -12.10_dbl, z/)

! Tungsten
  ptable(74) % symbol = 'W '
  ptable(74) % name = 'Tungsten'
  ptable(74) % number = 74
  ptable(74) % amass = 183.85000_dbl
  ptable(74) % mass = 183.95100_dbl
  ptable(74) % covalent_radius = 1.3000_dbl
  ptable(74) % vdw_radius = z
  ptable(74) % e_conv(0:3) = (/ 12, 24, 24, 14 /)
  ptable(74) % eht_param(0:3) = (/ -8.26_dbl, -5.17_dbl, -10.37_dbl, z/)

! Rhenium
  ptable(75) % symbol = 'Re'
  ptable(75) % name = 'Rhenium'
  ptable(75) % number = 75
  ptable(75) % amass = 186.20000_dbl
  ptable(75) % mass = 186.95600_dbl
  ptable(75) % covalent_radius = 1.2800_dbl
  ptable(75) % vdw_radius = z
  ptable(75) % e_conv(0:3) = (/ 12, 24, 25, 14 /)
  ptable(75) % eht_param(0:3) = (/ -9.36_dbl, -5.96_dbl, -12.66_dbl, z/)

! Osmium
  ptable(76) % symbol = 'Os'
  ptable(76) % name = 'Osmium'
  ptable(76) % number = 76
  ptable(76) % amass = 190.20000_dbl
  ptable(76) % mass = 192.00000_dbl
  ptable(76) % covalent_radius = 1.2600_dbl
  ptable(76) % vdw_radius = z
  ptable(76) % e_conv(0:3) = (/ 12, 24, 26, 14 /)
  ptable(76) % eht_param(0:3) = (/ -8.17_dbl, -4.81_dbl, -11.84_dbl, z/)

! Iridium
  ptable(77) % symbol = 'Ir'
  ptable(77) % name = 'Iridium'
  ptable(77) % number = 77
  ptable(77) % amass = 192.20000_dbl
  ptable(77) % mass = 192.96330_dbl
  ptable(77) % covalent_radius = 1.2700_dbl
  ptable(77) % vdw_radius = z
  ptable(77) % e_conv(0:3) = (/ 12, 24, 27, 14 /)
  ptable(77) % eht_param(0:3) = (/ -11.36_dbl, -4.50_dbl, -12.17_dbl, z/)

! Platinum
  ptable(78) % symbol = 'Pt'
  ptable(78) % name = 'Platinum'
  ptable(78) % number = 78
  ptable(78) % amass = 195.09000_dbl
  ptable(78) % mass = 194.96480_dbl
  ptable(78) % covalent_radius = 1.3000_dbl
  ptable(78) % vdw_radius = z
  ptable(78) % e_conv(0:3) = (/ 11, 24, 29, 14 /)
  ptable(78) % eht_param(0:3) = (/ -9.077_dbl, -5.475_dbl, -12.59_dbl, &
          z/)

! Gold
  ptable(79) % symbol = 'Au'
  ptable(79) % name = 'Gold'
  ptable(79) % number = 79
  ptable(79) % amass = 196.96700_dbl
  ptable(79) % mass = 196.96660_dbl
  ptable(79) % covalent_radius = 1.3400_dbl
  ptable(79) % vdw_radius = z
  ptable(79) % e_conv(0:3) = (/ 11, 24, 30, 14 /)
  ptable(79) % eht_param(0:3) = (/ -10.92_dbl, -5.55_dbl, -15.076_dbl, &
          z/)

! Mercury
  ptable(80) % symbol = 'Hg'
  ptable(80) % name = 'Mercury'
  ptable(80) % number = 80
  ptable(80) % amass = 200.59000_dbl
  ptable(80) % mass = 201.97060_dbl
  ptable(80) % covalent_radius = 1.4900_dbl
  ptable(80) % vdw_radius = z
  ptable(80) % e_conv(0:3) = (/ 12, 24, 30, 14 /)
  ptable(80) % eht_param(0:3) = (/ -13.68_dbl, -8.47_dbl, -17.50_dbl, z/)

! Thallium
  ptable(81) % symbol = 'Tl'
  ptable(81) % name = 'Thallium'
  ptable(81) % number = 81
  ptable(81) % amass = 204.37000_dbl
  ptable(81) % mass = 204.97450_dbl
  ptable(81) % covalent_radius = 1.4800_dbl
  ptable(81) % vdw_radius = z
  ptable(81) % e_conv(0:3) = (/ 12, 25, 30, 14 /)
  ptable(81) % eht_param(0:3) = (/ -11.60_dbl, -5.80_dbl, z, z/)

! Lead
  ptable(82) % symbol = 'Pb'
  ptable(82) % name = 'Lead'
  ptable(82) % number = 82
  ptable(82) % amass = 207.19000_dbl
  ptable(82) % mass = 207.97660_dbl
  ptable(82) % covalent_radius = 1.4700_dbl
  ptable(82) % vdw_radius = z
  ptable(82) % e_conv(0:3) = (/ 12, 26, 30, 14 /)
  ptable(82) % eht_param(0:3) = (/ -15.70_dbl, -8.00_dbl, z, z/)

! Bismuth
  ptable(83) % symbol = 'Bi'
  ptable(83) % name = 'Bismuth'
  ptable(83) % number = 83
  ptable(83) % amass = 208.98000_dbl
  ptable(83) % mass = 208.98040_dbl
  ptable(83) % covalent_radius = 1.4600_dbl
  ptable(83) % vdw_radius = z
  ptable(83) % e_conv(0:3) = (/ 12, 27, 30, 14 /)
  ptable(83) % eht_param(0:3) = (/ -15.19_dbl, -7.79_dbl, z, z/)

! Polonium
  ptable(84) % symbol = 'Po'
  ptable(84) % name = 'Polonium'
  ptable(84) % number = 84
  ptable(84) % amass = 209.98290_dbl
  ptable(84) % mass = 209.98290_dbl
  ptable(84) % covalent_radius = 1.4600_dbl
  ptable(84) % vdw_radius = z
  ptable(84) % e_conv(0:3) = (/ 12, 28, 30, 14 /)
  ptable(84) % eht_param(0:3) = (/ z, z, z, z/)

! Astatine
  ptable(85) % symbol = 'At'
  ptable(85) % name = 'Astatine'
  ptable(85) % number = 85
  ptable(85) % amass = 209.98700_dbl
  ptable(85) % mass = 209.98700_dbl
  ptable(85) % covalent_radius = 1.4500_dbl
  ptable(85) % vdw_radius = z
  ptable(85) % e_conv(0:3) = (/ 12, 29, 30, 14 /)
  ptable(85) % eht_param(0:3) = (/ z, z, z, z/)

! Radon
  ptable(86) % symbol = 'Rn'
  ptable(86) % name = 'Radon'
  ptable(86) % number = 86
  ptable(86) % amass = 222.01750_dbl
  ptable(86) % mass = 222.01750_dbl
  ptable(86) % covalent_radius = z
  ptable(86) % vdw_radius = z
  ptable(86) % e_conv(0:3) = (/ 12, 30, 30, 14 /)
  ptable(86) % eht_param(0:3) = (/ z, z, z, z/)

! Francium
  ptable(87) % symbol = 'Fr'
  ptable(87) % name = 'Francium'
  ptable(87) % number = 87
  ptable(87) % amass = 223.01980_dbl
  ptable(87) % mass = 223.01980_dbl
  ptable(87) % covalent_radius = z
  ptable(87) % vdw_radius = z
  ptable(87) % e_conv(0:3) = (/ 13, 30, 30, 14 /)
  ptable(87) % eht_param(0:3) = (/ z, z, z, z/)

! Radium
  ptable(88) % symbol = 'Ra'
  ptable(88) % name = 'Radium'
  ptable(88) % number = 88
  ptable(88) % amass = 226.02540_dbl
  ptable(88) % mass = 226.02540_dbl
  ptable(88) % covalent_radius = z
  ptable(88) % vdw_radius = z
  ptable(88) % e_conv(0:3) = (/ 14, 30, 30, 14 /)
  ptable(88) % eht_param(0:3) = (/ z, z, z, z/)

! Actinium
  ptable(89) % symbol = 'Ac'
  ptable(89) % name = 'Actinium'
  ptable(89) % number = 89
  ptable(89) % amass = 227.02780_dbl
  ptable(89) % mass = 227.02780_dbl
  ptable(89) % covalent_radius = z
  ptable(89) % vdw_radius = z
  ptable(89) % e_conv(0:3) = (/ 14, 30, 31, 14 /)
  ptable(89) % eht_param(0:3) = (/ z, z, z, z/)

! Thorium
  ptable(90) % symbol = 'Th'
  ptable(90) % name = 'Thorium'
  ptable(90) % number = 90
  ptable(90) % amass = 232.03810_dbl
  ptable(90) % mass = 232.03810_dbl
  ptable(90) % covalent_radius = 1.6500_dbl
  ptable(90) % vdw_radius = z
  ptable(90) % e_conv(0:3) = (/ 14, 30, 32, 14 /)
  ptable(90) % eht_param(0:3) = (/ -5.39_dbl, -5.39_dbl, -10.11_dbl, &
          -9.64_dbl/)

! Proctactinium
  ptable(91) % symbol = 'Pa'
  ptable(91) % name = 'Proctactinium'
  ptable(91) % number = 91
  ptable(91) % amass = 231.03590_dbl
  ptable(91) % mass = 231.03590_dbl
  ptable(91) % covalent_radius = z
  ptable(91) % vdw_radius = z
  ptable(91) % e_conv(0:3) = (/ 14, 30, 31, 16 /)
  ptable(91) % eht_param(0:3) = (/ z, z, z, z/)

! Uranium
  ptable(92) % symbol = 'U '
  ptable(92) % name = 'Uranium'
  ptable(92) % number = 92
  ptable(92) % amass = 238.05080_dbl
  ptable(92) % mass = 238.05080_dbl
  ptable(92) % covalent_radius = 1.4200_dbl
  ptable(92) % vdw_radius = z
  ptable(92) % e_conv(0:3) = (/ 14, 30, 31, 17 /)
  ptable(92) % eht_param(0:3) = (/ -5.50_dbl, -5.50_dbl, -9.19_dbl, &
          -10.62_dbl/)

! Neptunium
  ptable(93) % symbol = 'Np'
  ptable(93) % name = 'Neptunium'
  ptable(93) % number = 93
  ptable(93) % amass = 237.04820_dbl
  ptable(93) % mass = 237.04820_dbl
  ptable(93) % covalent_radius = z
  ptable(93) % vdw_radius = z
  ptable(93) % e_conv(0:3) = (/ 14, 30, 31, 18 /)
  ptable(93) % eht_param(0:3) = (/ z, z, z, z/)

! Plutonium
  ptable(94) % symbol = 'Pu'
  ptable(94) % name = 'Plutonium'
  ptable(94) % number = 94
  ptable(94) % amass = 244.0640_dbl
  ptable(94) % mass = 244.0640_dbl
  ptable(94) % covalent_radius = z
  ptable(94) % vdw_radius = z
  ptable(94) % e_conv(0:3) = (/ 14, 30, 30, 20 /)
  ptable(94) % eht_param(0:3) = (/ z, z, z, z/)

! Americum
  ptable(95) % symbol = 'Am'
  ptable(95) % name = 'Americum'
  ptable(95) % number = 95
  ptable(95) % amass = 243.0614_dbl
  ptable(95) % mass = 243.0614_dbl
  ptable(95) % covalent_radius = z
  ptable(95) % vdw_radius = z
  ptable(95) % e_conv(0:3) = (/ 14, 30, 30, 21 /)
  ptable(95) % eht_param(0:3) = (/ z, z, z, z/)

! Curium
  ptable(96) % symbol = 'Cm'
  ptable(96) % name = 'Curium'
  ptable(96) % number = 96
  ptable(96) % amass = 247.0700_dbl
  ptable(96) % mass = 247.0700_dbl
  ptable(96) % covalent_radius = z
  ptable(96) % vdw_radius = z
  ptable(96) % e_conv(0:3) = (/ 14, 30, 31, 21 /)
  ptable(96) % eht_param(0:3) = (/ z, z, z, z/)

! Berkelium
  ptable(97) % symbol = 'Bk'
  ptable(97) % name = 'Berkelium'
  ptable(97) % number = 97
  ptable(97) % amass = 251.0800_dbl
  ptable(97) % mass = 251.0800_dbl
  ptable(97) % covalent_radius = z
  ptable(97) % vdw_radius = z
  ptable(97) % e_conv(0:3) = (/ 14, 30, 31, 22 /)
  ptable(97) % eht_param(0:3) = (/ z, z, z, z/)

! Californium
  ptable(98) % symbol = 'Cf'
  ptable(98) % name = 'Californium'
  ptable(98) % number = 98
  ptable(98) % amass = 252.0820_dbl
  ptable(98) % mass = 252.0820_dbl
  ptable(98) % covalent_radius = z
  ptable(98) % vdw_radius = z
  ptable(98) % e_conv(0:3) = (/ 14, 30, 30, 24 /)
  ptable(98) % eht_param(0:3) = (/ z, z, z, z/)

! Einsteinium
  ptable(99) % symbol = 'Es'
  ptable(99) % name = 'Einsteinium'
  ptable(99) % number = 99
  ptable(99) % amass = 252.0829_dbl
  ptable(99) % mass = 252.0829_dbl
  ptable(99) % covalent_radius = z
  ptable(99) % vdw_radius = z
  ptable(99) % e_conv(0:3) = (/ 14, 30, 30, 25 /)
  ptable(99) % eht_param(0:3) = (/ z, z, z, z/)

! Fermium
  ptable(100) % symbol = 'Fm'
  ptable(100) % name = 'Fermium'
  ptable(100) % number = 100
  ptable(100) % amass = 257.0950_dbl
  ptable(100) % mass = 257.0950_dbl
  ptable(100) % covalent_radius = z
  ptable(100) % vdw_radius = z
  ptable(100) % e_conv(0:3) = (/ 14, 30, 30, 26 /)
  ptable(100) % eht_param(0:3) = (/ z, z, z, z/)

! Mendelevium
  ptable(101) % symbol = 'Md'
  ptable(101) % name = 'Mendelevium'
  ptable(101) % number = 101
  ptable(101) % amass = 256.0000_dbl
  ptable(101) % mass = 256.0000_dbl
  ptable(101) % covalent_radius = z
  ptable(101) % vdw_radius = z
  ptable(101) % e_conv(0:3) = (/ 14, 30, 30, 27 /)
  ptable(101) % eht_param(0:3) = (/ z, z, z, z/)

! Nobelium
  ptable(102) % symbol = 'No'
  ptable(102) % name = 'Nobelium'
  ptable(102) % number = 102
  ptable(102) % amass = 254.0000_dbl
  ptable(102) % mass = 254.0000_dbl
  ptable(102) % covalent_radius = z
  ptable(102) % vdw_radius = z
  ptable(102) % e_conv(0:3) = (/ 14, 30, 30, 28 /)
  ptable(102) % eht_param(0:3) = (/ z, z, z, z/)

! Lawrencium
  ptable(103) % symbol = 'Lr'
  ptable(103) % name = 'Lawrencium'
  ptable(103) % number = 103
  ptable(103) % amass = 257.0000_dbl
  ptable(103) % mass = 257.0000_dbl
  ptable(103) % covalent_radius = z
  ptable(103) % vdw_radius = z
  ptable(103) % e_conv(0:3) = (/ 14, 30, 31, 28 /)
  ptable(103) % eht_param(0:3) = (/ z, z, z, z/)

! Unnilquadium
  ptable(104) % symbol = 'Uq'
  ptable(104) % name = 'Unnilquadium'
  ptable(104) % number = 104
  ptable(104) % amass = 261.0000_dbl
  ptable(104) % mass = 261.0000_dbl
  ptable(104) % covalent_radius = z
  ptable(104) % vdw_radius = z
  ptable(104) % e_conv(0:3) = (/ 14, 30, 32, 28 /)
  ptable(104) % eht_param(0:3) = (/ z, z, z, z/)

! Unnilpentium
  ptable(105) % symbol = 'Up'
  ptable(105) % name = 'Unnilpentium'
  ptable(105) % number = 105
  ptable(105) % amass = 262.0000_dbl
  ptable(105) % mass = 262.0000_dbl
  ptable(105) % covalent_radius = z
  ptable(105) % vdw_radius = z
  ptable(105) % e_conv(0:3) = (/ 14, 30, 33, 28 /)
  ptable(105) % eht_param(0:3) = (/ z, z, z, z/)

! Unnilhexium
  ptable(106) % symbol = 'Uh'
  ptable(106) % name = 'Unnilhexium'
  ptable(106) % number = 106
  ptable(106) % amass = 263.0000_dbl
  ptable(106) % mass = 263.0000_dbl
  ptable(106) % covalent_radius = z
  ptable(106) % vdw_radius = z
  ptable(106) % e_conv(0:3) = (/ 14, 30, 34, 28 /)
  ptable(106) % eht_param(0:3) = (/ z, z, z, z/)

END SUBROUTINE init_periodic_table

!!*****
!******************************************************************************

END MODULE periodic_table

!******************************************************************************
