!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Default parameter sets for SCPTB method
!>      key, atomname            key, atomic name
!>      z, zeff                  atomic number, core charge
!>      lmaxorb, lmaxscp         max l basis, max l polarization
!>      norb                     number of orbitals (0,1,2,3)
!>      zeta                     exponents of orbitals (0,1,2,3)
!>      hcore                    two body terms
!>      occupation               orbital occupation in free atom
!>      energy                   total energy of free atom
!>      crep                     core repulsion terms (1=a,2=c,3=d)
!>      pol                      polarization constants (1,2,3)
!>      ag                       Gaussian exponents
!> \par History
!> \author JGH (12.2012)
! *****************************************************************************
MODULE scptb_parameters
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE periodic_table,                  ONLY: ptable
  USE scptb_types,                     ONLY: scptb_parameter_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scptb_parameters'

  INTEGER, PARAMETER, PRIVATE :: nelem = 106

  PUBLIC :: scptb_default_parameter

CONTAINS

! *****************************************************************************
!> \brief Default parameter sets for SCPTB method
!> \param param ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE scptb_default_parameter ( param, error)
    TYPE(scptb_parameter_type), &
      INTENT(INOUT)                          :: param
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scptb_default_parameter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, l
    LOGICAL                                  :: failure

    param%parameterization = "UZH.00.00"
    SELECT CASE (TRIM(param%key))
    CASE DEFAULT
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    CASE ("XX")
       param%defined          = .TRUE.
    CASE ("H")
       !Hydrogen
       param%defined               = .TRUE.
       param%key                   = "H"
       param%z                     = 1
       param%zeff                  = 1
       param%lmaxorb               = 0
       param%lmaxscp               = 1
       param%norb(0:3)             = (/1,0,0,0/)
       param%nqm(1,0)              = 1
       param%zeta(1,0)             = 1.2_dp
       param%hcore(1,0)            = 1.0_dp
       param%occupation(1,0)       = 1
       param%crep(1:3)             = (/1.0_dp,1.0_dp,1.0_dp/)
       param%pol(1:3)              = (/1.0_dp,1.0_dp,1.0_dp/)
       param%ag                    = 4.0_dp
       param%energy                = -0.5_dp
!   CASE (2)
!      !Helium
!   CASE (3)
!      !Lithium
!   CASE (4)
!      !Beryllium
!   CASE (5)
!      !Boron
!   CASE (6)
!      !Carbon
!   CASE (7)
!      !Nitrogen
    CASE ("O")
       !Oxygen
       param%defined               = .TRUE.
       param%key                   = "O"
       param%z                     = 8
       param%zeff                  = 6
       param%lmaxorb               = 1
       param%lmaxscp               = 2
       param%norb(0:3)             = (/1,1,0,0/)
       param%nqm(1,0:1)            = (/2,2/)
       param%zeta(1,0:1)           = (/1.8_dp,1.8_dp/)
       param%hcore(1,0:1)          = (/1.0_dp,1.0_dp/)
       param%occupation(1,0:1)     = (/2,4/)
       param%crep(1:3)             = (/1.0_dp,1.0_dp,1.0_dp/)
       param%pol(1:3)              = (/1.0_dp,1.0_dp,1.0_dp/)
       param%ag                    = 1.0_dp
       param%energy                = -0.5_dp
    CASE ("F")
       !Fluorine
       param%defined               = .TRUE.
       param%key                   = "F"
       param%z                     = 9
       param%zeff                  = 7
       param%lmaxorb               = 1
       param%lmaxscp               = 2
       param%norb(0:3)             = (/1,1,0,0/)
       param%nqm(1,0:1)            = (/2,2/)
       param%zeta(1,0:1)           = (/1.8_dp,1.8_dp/)
       param%hcore(1,0:1)          = (/1.0_dp,1.0_dp/)
       param%occupation(1,0:1)     = (/2,5/)
       param%crep(1:3)             = (/1.0_dp,1.0_dp,1.0_dp/)
       param%pol(1:3)              = (/1.0_dp,1.0_dp,1.0_dp/)
       param%ag                    = 4.0_dp
       param%energy                = -0.5_dp
!   CASE (10)
!      !Neon
!   CASE (11)
!      !Sodium
!   CASE (12)
!      !Magnesium
!   CASE (13)
!      !Aluminium
!   CASE (14)
!      !Silicon
!   CASE (15)
!      !Phosphorus
!   CASE (16)
!      !Sulphur
!   CASE (17)
!      !Chlorine
!   CASE (18)
!      !Argon
!   CASE (19)
!      !Potassium
!   CASE (20)
!      !Calcium
!   CASE (21)
!   CASE (22)
!   CASE (23)
!   CASE (24)
!   CASE (25)
!   CASE (26)
!   CASE (27)
!   CASE (28)
!   CASE (29)
!   CASE (30)
!      !Zinc
!   CASE (31)
!      !Gallium
!   CASE (32)
!      !Germanium
!   CASE (33)
!      !Arsenic
!   CASE (34)
!      !Selenium
!   CASE (35)
!      !Bromine
!   CASE (36)
!   CASE (37)
!      !Rubidium
!   CASE (38)
!      !Strontium
!   CASE (39)
!   CASE (40)
!   CASE (41)
!   CASE (42)
!   CASE (43)
!   CASE (44)
!   CASE (45)
!   CASE (46)
!   CASE (47)
!   CASE (48)
!      !Cadmium
!   CASE (49)
!      !Indium
!   CASE (50)
!      !Tin
!   CASE (51)
!      !Antimony
!   CASE (52)
!      !Tellurium
!   CASE (53)
!      !Iodine
!   CASE (54)
!   CASE (55)
!      !Cesium
!   CASE (56)
!      !Barium
!   CASE (57)
!   CASE (58)
!   CASE (59)
!   CASE (60)
!   CASE (61)
!   CASE (62)
!   CASE (63)
!   CASE (64)
!   CASE (65)
!   CASE (66)
!   CASE (67)
!   CASE (68)
!   CASE (69)
!   CASE (70)
!   CASE (71)
!   CASE (72)
!   CASE (73)
!   CASE (74)
!   CASE (75)
!   CASE (76)
!   CASE (77)
!   CASE (78)
!   CASE (79)
!   CASE (80)
!      !Mercury
!   CASE (81)
!      !Thallium
!   CASE (82)
!      !Lead
!   CASE (83)
!      !Bismuth
!   CASE (84)
!   CASE (85)
!   CASE (86)
!   CASE (87)
!   CASE (88)
!   CASE (89)
!   CASE (90)
!   CASE (91)
!   CASE (92)
!   CASE (93)
!   CASE (94)
!   CASE (95)
!   CASE (96)
!   CASE (97)
!   CASE (98)
!   CASE (99:nelem)
    END SELECT

    IF(param%defined) THEN
      param%natorb = 0
      DO l=0,3
        param%natorb = param%natorb + ((l+1)*(l+2))/2 * param%norb(l)
      END DO
      ! check atomic name
      i = param%z
      param%atomname = ptable(i)%symbol
    END IF

  END SUBROUTINE scptb_default_parameter

END MODULE scptb_parameters

