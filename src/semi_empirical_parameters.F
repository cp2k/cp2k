!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/semi_empirical_parameters [1.0] *
!!
!!   NAME
!!     semi_empirical_parameters
!!
!!   FUNCTION
!!     Default parameter sets for semi empirical models: sep%...
!!     ad, am, aq,    a.u.      parameters for the multipole approximation to
!!     dd, qq         a.u.      the two-centre, two-electron integrals.
!!     alp            A**(-1)   the exponents used to calculate the core-core
!!                              repulsion terms.
!!     beta           eV        parameters for the two-centre, one-electron core
!!                              resonance integrals.
!!     z                        the CORE charge on each element.
!!     eheat          kcal/mol  the experimental heat of formation for each element.
!!     eisol          eV        the calculated electronic energies for each element.
!!     fn1            A*eV      AM1/PM3 specific parameters
!!     fn2            A**(-2)   for the core-core repulsion
!!     fn3            A         interaction.                       
!!     pre            eV        PDDG function                      
!!     d              A         PDDG function                      
!!     bfn1-3         special core-core interactions
!!                        
!!     gdd, gpd, gsd, eV        parameters for the Coulomb and exchange one-centre,
!!     gp2, gpp, gsp, eV        two-electron integrals.
!!     gss, hsp       eV
!!     uss, upp, udd  eV        electron kinetic energy integrals.
!!     sto_exponents  a.u.      the Slater exponents of the basis functions.
!!
!!     all dimensions will be converted to a.u. at the end of the routines
!!
!!   AUTHOR
!!     JGH (14.08.2004)
!!
!!   MODIFICATION HISTORY
!!     BK   (12.10.2004)
!!     JGH  (36.03.2006) added PDDG function
!!
!!   SOURCE
!******************************************************************************

MODULE semi_empirical_parameters
  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_mndo,&
                                             do_method_pchg,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_rm1
  USE kinds,                           ONLY: dp
  USE periodic_table,                  ONLY: nelem,&
                                             ptable
  USE physcon,                         ONLY: bohr,&
                                             evolt,&
                                             kcalmol
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_parameters'

! *** Public subroutines ***

  PUBLIC :: mndo_default_parameter, &
            am1_default_parameter, &
            rm1_default_parameter, &
            pdg_default_parameter, &
            pm3_default_parameter,&
            pcharge_default_parameter,&
            slater_exponents,srules

!              H           Initial "s" Orbital Occupancies                     He
!              Li Be                                            B  C  N  O  F  Ne
!              Na Mg                                            Al Si P  S  Cl Ar
!              K  Ca Sc            Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y             Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu      Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U    Np Pu Am Cm Bk Cf            Cb ++ +  -- -  Tv
!                                      "s" shell
  INTEGER, DIMENSION(107), PRIVATE :: Nos =(/&
       1,                                                                2, &!    2
       1, 2,                                              2, 2, 2, 2, 2, 0, &!   10
       1, 2,                                              2, 2, 2, 2, 2, 0, &!   18
       1, 2, 2,              2, 2, 1, 2, 2, 2, 2, 1, 2,   2, 2, 2, 2, 2, 0, &!   36
       1, 2, 2,              2, 1, 1, 2, 1, 1, 0, 1, 2,   2, 2, 2, 2, 2, 0, &!   54
       1, 2, 2, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2,                   &
                             2, 2, 1, 2, 2, 2, 1, 1, 2,   2, 2, 2, 2, 2, 0, &!   86
       1, 1, 2, 4, 2, 2,     2, 2, 2, 2, 2, 2, 0, 3,-3,   1, 2, 1,-2,-1, 0 /)
!
!
!              H           Initial "p" Orbital Occupancies                   He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf (The rest are reserved for MOPAC)
!                                      "p" shell
  INTEGER, DIMENSION(107), PRIVATE :: Nop =(/&
          0 ,                                                           0, &!    2
          0, 0,                                          1, 2, 3, 4, 5, 6, &!   10
          0, 0,                                          1, 2, 3, 4, 5, 6, &!   18
          0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   36
          0, 0, 0,          0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   54
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,               &
                            0, 0, 0, 0, 0, 0, 0, 0, 0,   1, 2, 3, 4, 5, 6, &!   86
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
!
!              H           Initial "d" Orbital Occupancies                   He
!              Li Be                                          B  C  N  O  F  Ne
!              Na Mg                                          Al Si P  S  Cl Ar
!              K  Ca Sc          Ti V  Cr Mn Fe Co Ni Cu Zn   Ga Ge As Se Br Kr
!              Rb Sr Y           Zr Nb Mo Tc Ru Rh Pd Ag Cd   In Sn Sb Te I  Xe
!              Cs Ba La Ce-Lu    Hf Ta W  Re Os Ir Pt Au Hg   Tl Pb Bi Po At Rn
!              Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf (The rest are reserved for MOPAC)
!                                      "d" shell
  INTEGER, DIMENSION(107), PRIVATE :: Nod =(/&
          0,                                                           0, &!    2
          0, 0,                                         0, 0, 0, 0, 0, 0, &!   10
          0, 0,                                         0, 0, 0, 0, 0, 0, &!   18
          0, 0, 1,          2, 3, 5, 5, 6, 7, 8, 10, 0, 0, 0, 0, 0, 0, 0, &!   36
          0, 0, 1,          2, 4, 5, 5, 7, 8,10, 10, 0, 0, 0, 0, 0, 0, 0, &!   54
          0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 1,             &
                            2, 3, 5, 5, 6, 7, 9, 10, 0, 0, 0, 0, 0, 0, 0, &!   86
          0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0 /)


!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE pcharge_default_parameter(sep, z)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z

    CHARACTER(len=*), PARAMETER :: routineN = 'pcharge_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ  = do_method_pchg
    sep%name = "CLASSICAL"
    sep%eheat = 0.0_dp
    sep%z = z
    CALL valence_electrons ( z, sep%zeff, sep%natorb )
    sep%ngauss = 0
    sep%defined = .TRUE.
    sep%sto_exponents(0:3) = 0._dp !eV
    sep%beta(0:3) = 0._dp !eV
    sep%uss = 0._dp; sep%upp = 0._dp; sep%udd = 0._dp !eV
    sep%alp = 0._dp
    sep%eisol = 0._dp
    sep%gss = 0._dp; sep%gsp = 0._dp; sep%gpp = 0._dp; sep%gp2 = 0._dp
    sep%gsd = 0._dp; sep%gpd = 0._dp; sep%gdd = 0._dp
    sep%hsp = 0._dp
    sep%dd = 0._dp; sep%qq = 0._dp; sep%am = 0._dp; sep%ad = 0._dp; sep%aq = 0._dp
    sep%fn1(1:4) = 0._dp; sep%fn2(1:4) = 0._dp; sep%fn3(1:4) = 0._dp
    sep%bfn1 = 0._dp; sep%bfn2 = 0._dp; sep%bfn3 = 0._dp

  END SUBROUTINE pcharge_default_parameter

! *****************************************************************************

  SUBROUTINE mndo_default_parameter ( sep, z )

    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z

    CHARACTER(len=*), PARAMETER :: routineN = 'mndo_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ = do_method_mndo
    sep%name = ptable(z)%name
    sep%eheat = ptable(z)%heat_of_formation / kcalmol
    sep%z = z
    CALL valence_electrons ( z, sep%zeff, sep%natorb )
    sep%ngauss = 6
    sep%defined = .FALSE.
    sep%sto_exponents(0:3) = 0._dp !eV
    sep%beta(0:3) = 0._dp !eV
    sep%uss = 0._dp; sep%upp = 0._dp; sep%udd = 0._dp !eV
    sep%alp = 0._dp
    sep%eisol = 0._dp
    sep%gss = 0._dp; sep%gsp = 0._dp; sep%gpp = 0._dp; sep%gp2 = 0._dp
    sep%gsd = 0._dp; sep%gpd = 0._dp; sep%gdd = 0._dp
    sep%hsp = 0._dp
    sep%dd = 0._dp; sep%qq = 0._dp; sep%am = 0._dp; sep%ad = 0._dp; sep%aq = 0._dp
    sep%fn1(1:4) = 0._dp; sep%fn2(1:4) = 0._dp; sep%fn3(1:4) = 0._dp
    sep%bfn1 = 0._dp; sep%bfn2 = 0._dp; sep%bfn3 = 0._dp

    SELECT CASE (z)
      CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "invalid atomic number")
      CASE (0)
        sep%defined = .TRUE.
      CASE (1)
        !Hydrogen
        sep%defined          = .TRUE.
        sep%alp              =   2.5441341_dp
        sep%beta(0)          =  -6.9890640_dp
        sep%sto_exponents(0) =   1.3319670_dp
        sep%uss              = -11.9062760_dp
        sep%gss              =  12.8480000_dp
      CASE (2)
      CASE (3)
        !Lithium
        sep%defined          = .TRUE.
        sep%alp              =   1.2501400_dp
        sep%beta(0)          =  -1.3500400_dp
        sep%beta(1)          =  -1.3500400_dp
        sep%sto_exponents(0) =   0.7023800_dp 
        sep%sto_exponents(1) =   0.7023800_dp
        sep%uss              =  -5.1280000_dp 
        sep%upp              =  -2.7212000_dp
        sep%gss              =   7.3000000_dp 
        sep%gsp              =   5.4200000_dp 
        sep%gpp              =   5.0000000_dp
        sep%gp2              =   4.5200000_dp
        sep%hsp              =   0.8300000_dp
      CASE (4)
      !Beryllium
        sep%defined          = .TRUE.
        sep%alp              =   1.6694340_dp
        sep%beta(0) = -4.0170960_dp
        sep%beta(1) = -4.0170960_dp
        sep%sto_exponents(0) = 1.0042100_dp 
        sep%sto_exponents(1) = 1.0042100_dp
        sep%uss = -16.6023780_dp 
        sep%upp = -10.7037710_dp 
        sep%gss = 9.0000000_dp 
        sep%gsp = 7.4300000_dp 
        sep%gpp = 6.9700000_dp 
        sep%gp2 = 6.2200000_dp
        sep%hsp              =   1.2800000_dp
      CASE (5)
        !Boron
        sep%defined = .TRUE.
        sep%alp   =   2.1349930_dp
        sep%beta(0) =  -8.2520540_dp
        sep%beta(1) =  -8.2520540_dp
        sep%sto_exponents(0)    =   1.5068010_dp
        sep%sto_exponents(1)    =   1.5068010_dp
        sep%uss   = -34.5471300_dp
        sep%upp   = -23.1216900_dp
        sep%gss   =  10.5900000_dp
        sep%gsp   =   9.5600000_dp
        sep%gpp   =   8.8600000_dp
        sep%gp2   =   7.8600000_dp
        sep%hsp   =   1.8100000_dp
      CASE (6)
        !Carbon
        sep%defined = .TRUE.
        sep%alp   =    2.5463800_dp
        sep%beta(0) =  -18.9850440_dp
        sep%beta(1) =   -7.9341220_dp
        sep%sto_exponents(0)    =    1.7875370_dp
        sep%sto_exponents(1)    =    1.7875370_dp
        sep%uss   =  -52.2797450_dp
        sep%upp   =  -39.2055580_dp
        sep%gss   =   12.2300000_dp
        sep%gsp   =   11.4700000_dp
        sep%gpp   =   11.0800000_dp
        sep%gp2   =    9.8400000_dp
        sep%hsp   =    2.4300000_dp
      CASE (7)
      !Nitrogen
        sep%defined = .TRUE.
        sep%alp   =    2.8613420_dp
        sep%beta(0) =  -20.4957580_dp
        sep%beta(1) =  -20.4957580_dp
        sep%sto_exponents(0)    =    2.2556140_dp
        sep%sto_exponents(1)    =    2.2556140_dp
        sep%uss   =  -71.9321220_dp
        sep%upp   =  -57.1723190_dp
        sep%gss   =   13.5900000_dp
        sep%gsp   =   12.6600000_dp
        sep%gpp   =   12.9800000_dp
        sep%gp2   =   11.5900000_dp
        sep%hsp   =    3.1400000_dp
      CASE (8)
      !Oxygen
        sep%defined = .TRUE.
        sep%alp   =    3.1606040_dp
        sep%beta(0) =  -32.6880820_dp
        sep%beta(1) =  -32.6880820_dp
        sep%sto_exponents(0)    =    2.6999050_dp
        sep%sto_exponents(1)    =    2.6999050_dp
        sep%uss   =  -99.6443090_dp
        sep%upp   =  -77.7974720_dp
        sep%gss   =   15.4200000_dp
        sep%gsp   =   14.4800000_dp
        sep%gpp   =   14.5200000_dp
        sep%gp2   =   12.9800000_dp
        sep%hsp   =    3.9400000_dp
      CASE (9)
      !Fluorine
        sep%defined = .TRUE.
        sep%alp   =    3.4196606_dp
        sep%beta(0) =  -48.2904660_dp
        sep%beta(1) =  -36.5085400_dp
        sep%sto_exponents(0)    =    2.8484870_dp
        sep%sto_exponents(1)    =    2.8484870_dp
        sep%uss   = -131.0715480_dp
        sep%upp   = -105.7821370_dp
        sep%gss   =   16.9200000_dp
        sep%gsp   =   17.2500000_dp
        sep%gpp   =   16.7100000_dp
        sep%gp2   =   14.9100000_dp
        sep%hsp   =    4.8300000_dp
      CASE (10:12)
      CASE (13)
       !Aluminium
        sep%defined = .TRUE.
        sep%alp   =   1.8688394_dp
        sep%beta(0) =  -2.6702840_dp
        sep%beta(1) =  -2.6702840_dp
        sep%sto_exponents(0) =   1.4441610_dp
        sep%sto_exponents(1) =   1.4441610_dp
        sep%sto_exponents(2) =   1.0000000_dp
        sep%uss   = -23.8070970_dp
        sep%upp   = -17.5198780_dp
        sep%gss   =   8.0900000_dp
        sep%gsp   =   6.6300000_dp
        sep%gpp   =   5.9800000_dp
        sep%gp2   =   5.4000000_dp
        sep%hsp   =   0.7000000_dp
      CASE (14)
      !Silicon
        sep%defined = .TRUE.
        sep%alp   =   2.2053160_dp
        sep%beta(0) =  -9.0868040_dp
        sep%beta(1) =  -1.0758270_dp
        sep%sto_exponents(0)    =   1.3159860_dp
        sep%sto_exponents(1)    =   1.7099430_dp
        sep%sto_exponents(2)    =   1.0000000_dp
        sep%uss   = -37.0375330_dp
        sep%upp   = -27.7696780_dp
        sep%gss   =   9.8200000_dp
        sep%gsp   =   8.3600000_dp
        sep%gpp   =   7.3100000_dp
        sep%gp2   =   6.5400000_dp
        sep%hsp   =   1.3200000_dp
      CASE (15)
      !Phosphorus
        sep%defined = .TRUE.
        sep%alp   =    2.4152800_dp
        sep%beta(0) =   -6.7916000_dp
        sep%beta(1) =   -6.7916000_dp
        sep%sto_exponents(0)    =    2.1087200_dp
        sep%sto_exponents(1)    =    1.7858100_dp
        sep%sto_exponents(2)    =    1.0000000_dp
        sep%uss   =  -56.1433600_dp
        sep%upp   =  -42.8510800_dp
        sep%gss   =   11.5600000_dp
        sep%gsp   =   10.0800000_dp
        sep%gpp   =    8.6400000_dp
        sep%gp2   =    7.6800000_dp
        sep%hsp   =    1.9200000_dp
      CASE (16)
      !Sulphur
        sep%defined = .TRUE.
        sep%alp   =    2.4780260_dp
        sep%beta(0) =  -10.7616700_dp
        sep%beta(1) =  -10.1084330_dp
        sep%sto_exponents(0)    =    2.3129620_dp
        sep%sto_exponents(1)    =    2.0091460_dp
        sep%sto_exponents(2)    =    1.0000000_dp
        sep%uss   =  -72.2422810_dp
        sep%upp   =  -56.9732070_dp
        sep%gss   =   12.8800000_dp
        sep%gsp   =   11.2600000_dp
        sep%gpp   =    9.9000000_dp
        sep%gp2   =    8.8300000_dp
        sep%hsp   =    2.2600000_dp
      CASE (17)
      !Chlorine
        sep%defined = .TRUE.
        sep%alp   =    2.5422010_dp
        sep%beta(0) =  -14.2623200_dp
        sep%beta(1) =  -14.2623200_dp
        sep%sto_exponents(0)    =    3.7846450_dp
        sep%sto_exponents(1)    =    2.0362630_dp
        sep%sto_exponents(2)    =    1.0000000_dp
        sep%uss   = -100.2271660_dp
        sep%upp   =  -77.3786670_dp
        sep%gss   =   15.0300000_dp
        sep%gsp   =   13.1600000_dp
        sep%gpp   =   11.3000000_dp
        sep%gp2   =    9.9700000_dp
        sep%hsp   =    2.4200000_dp
      CASE (18:23)
      CASE (24)
        !Chromium.
        sep%defined = .TRUE.
        sep%alp   =    3.0683070_dp
        sep%beta(0) =   -0.1000000_dp
        sep%beta(1) =   -0.1000000_dp
        sep%beta(2) =   -8.7766360_dp
        sep%sto_exponents(0)    =    1.5000000_dp
        sep%sto_exponents(1)    =    1.5000000_dp
        sep%sto_exponents(2)    =    2.8845490_dp
        sep%uss   =  -17.5170270_dp
        sep%upp   =  -12.5337290_dp
        sep%udd   =  -44.1249280_dp
        sep%gss   =    6.0000000_dp
        sep%gsp   =    4.1500000_dp
        sep%gpp   =    5.0000000_dp
        sep%gp2   =    3.5000000_dp
        sep%hsp   =    1.0000000_dp
        sep%gsd   =    2.8746410_dp
        sep%gpd   =    3.0000000_dp
        sep%gdd   =    8.8949670_dp
      CASE (25:29)
      CASE (30)
        !Zinc
        sep%defined = .TRUE.
        sep%alp   =    1.5064570_dp
        sep%beta(0) =   -1.0000000_dp
        sep%beta(1) =   -2.0000000_dp
        sep%sto_exponents(0)    =    2.0473590_dp
        sep%sto_exponents(1)    =    1.4609460_dp
        sep%sto_exponents(2)    =    1.0000000_dp
        sep%uss   =  -20.8397160_dp
        sep%upp   =  -19.6252240_dp
        sep%gss   =   11.8000000_dp
        sep%gsp   =   11.1820180_dp
        sep%gpp   =   13.3000000_dp
        sep%gp2   =   12.9305200_dp
        sep%hsp   =    0.4846060_dp
      CASE (31)
      CASE (32)
        !Germanium
        sep%defined = .TRUE.
        sep%alp   =   1.9784980_dp
        sep%beta(0) =  -4.5164790_dp
        sep%beta(1) =  -1.7555170_dp
        sep%sto_exponents(0)    =   1.2931800_dp
        sep%sto_exponents(1)    =   2.0205640_dp
        sep%uss   = -33.9493670_dp
        sep%upp   = -27.4251050_dp
        sep%gss   =   9.8000000_dp
        sep%gsp   =   8.3000000_dp
        sep%gpp   =   7.3000000_dp
        sep%gp2   =   6.5000000_dp
        sep%hsp   =   1.3000000_dp
      CASE (33:34)
      CASE (35)
      !Bromine
        sep%defined = .TRUE.
        sep%alp   =    2.44570510_dp
        sep%beta(0) =   -8.91710700_dp
        sep%beta(1) =   -9.94374000_dp
        sep%sto_exponents(0)    =    3.85430190_dp
        sep%sto_exponents(1)    =    2.19920910_dp
        sep%sto_exponents(2)    =    1.00000000_dp
        sep%uss   =  -99.98644050_dp
        sep%upp   =  -75.67130750_dp
        sep%gss   =   15.03643948_dp
        sep%gsp   =   13.03468242_dp
        sep%gpp   =   11.27632539_dp
        sep%gp2   =    9.85442552_dp
        sep%hsp   =    2.45586832_dp
      CASE (36:49)
      CASE (50)
        !Tin
        sep%defined = .TRUE.
        sep%alp   =   1.8008140_dp
        sep%beta(0) =  -3.2351470_dp
        sep%beta(1) =  -4.2904160_dp
        sep%sto_exponents(0)    =   2.0803800_dp
        sep%sto_exponents(1)    =   1.9371060_dp
        sep%uss   = -40.8518020_dp
        sep%upp   = -28.5602490_dp
        sep%gss   =   9.8000000_dp
        sep%gsp   =   8.3000000_dp
        sep%gpp   =   7.3000000_dp
        sep%gp2   =   6.5000000_dp
        sep%hsp   =   1.3000000_dp
      CASE (51:52)
      CASE (53)
        !Iodine
        sep%defined = .TRUE.
        sep%uss   = -100.00305380_dp
        sep%upp   =  -74.61146920_dp
        sep%beta(0) =   -7.41445100_dp
        sep%beta(1) =   -6.19678100_dp
        sep%sto_exponents(0)    =    2.27296100_dp
        sep%sto_exponents(1)    =    2.16949800_dp
        sep%sto_exponents(2)    =    1.00000000_dp
        sep%alp   =    2.20732000_dp
        sep%gss   =   15.04044855_dp
        sep%gsp   =   13.05655798_dp
        sep%gpp   =   11.14778369_dp
        sep%gp2   =    9.91409071_dp
        sep%hsp   =    2.45638202_dp
      CASE (54:79)
      CASE (80)
        !Mercury
        sep%defined = .TRUE.
        sep%alp   =   1.3356410_dp
        sep%beta(0) =  -0.4045250_dp
        sep%beta(1) =  -6.2066830_dp
        sep%sto_exponents(0)    =   2.2181840_dp
        sep%sto_exponents(1)    =   2.0650380_dp
        sep%uss   = -19.8095740_dp
        sep%upp   = -13.1025300_dp
        sep%gss   =  10.8000000_dp
        sep%gsp   =   9.3000000_dp
        sep%gpp   =  14.3000000_dp
        sep%gp2   =  13.5000000_dp
        sep%hsp   =   1.3000000_dp
      CASE (81)
      CASE (82)
        !Lead
        sep%defined = .TRUE.
        sep%alp   =    1.7283330_dp
        sep%beta(0) =   -8.0423870_dp
        sep%beta(1) =   -3.0000000_dp
        sep%sto_exponents(0)    =    2.4982860_dp
        sep%sto_exponents(1)    =    2.0820710_dp
        sep%uss   =  -47.3196920_dp
        sep%upp   =  -28.8475600_dp
        sep%gss   =    9.8000000_dp
        sep%gsp   =    8.3000000_dp
        sep%gpp   =    7.3000000_dp
        sep%gp2   =    6.5000000_dp
        sep%hsp   =    1.3000000_dp
      CASE (83:106)
    END SELECT
    CALL calpar(z,sep)

    sep%beta = sep%beta/evolt
    sep%uss = sep%uss/evolt; sep%upp = sep%upp/evolt; sep%udd = sep%udd/evolt
    sep%alp = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss = sep%gss/evolt; sep%gsp = sep%gsp/evolt; sep%gpp = sep%gpp/evolt; sep%gp2 = sep%gp2/evolt 
    sep%gsd = sep%gsd/evolt; sep%gpd = sep%gpd/evolt; sep%gdd = sep%gdd/evolt
    sep%hsp = sep%hsp/evolt
    sep%fn1 = sep%fn1*bohr/evolt; sep%fn2 = sep%fn2/bohr/bohr; sep%fn3 = sep%fn3*bohr

  END SUBROUTINE mndo_default_parameter

! *****************************************************************************

  SUBROUTINE am1_default_parameter ( sep, z )
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z

    CHARACTER(len=*), PARAMETER :: routineN = 'am1_default_parameter', &
      routineP = moduleN//':'//routineN

! -----------------------------------------------------------------------------

    sep%typ = do_method_am1
    sep%name = ptable(z)%name
    sep%eheat = ptable(z)%heat_of_formation / kcalmol
    sep%z = z
    CALL valence_electrons ( z, sep%zeff, sep%natorb )
    sep%ngauss = 6
    sep%defined = .FALSE.
    sep%sto_exponents(0:3) = 0._dp
    sep%beta(0:3) = 0._dp
    sep%uss = 0._dp; sep%upp = 0._dp; sep%udd = 0._dp
    sep%alp = 0._dp
    sep%eisol = 0._dp
    sep%gss = 0._dp; sep%gsp = 0._dp; sep%gpp = 0._dp; sep%gp2 = 0._dp
    sep%gsd = 0._dp; sep%gpd = 0._dp; sep%gdd = 0._dp
    sep%hsp = 0._dp
    sep%dd = 0._dp; sep%qq = 0._dp; sep%am = 0._dp; sep%ad = 0._dp; sep%aq = 0._dp
    sep%fn1 = 0._dp; sep%fn2 = 0._dp; sep%fn3 = 0._dp
    sep%bfn1 = 0._dp; sep%bfn2 = 0._dp; sep%bfn3 = 0._dp

    SELECT CASE (z)
      CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "invalid atomic number")
      CASE (0)
        sep%defined = .TRUE.
      CASE (1)
        !Hydrogen
        sep%defined          = .TRUE.
        sep%alp   =   2.8823240_dp
        sep%beta(0) =  -6.1737870_dp
        sep%sto_exponents(0)    =   1.1880780_dp
        sep%uss   = -11.3964270_dp
        sep%gss   =  12.8480000_dp
        sep%fn1(1) =  0.1227960_dp; sep%fn2(1) = 5.0000000_dp; sep%fn3(1) = 1.2000000_dp
        sep%fn1(2) =  0.0050900_dp; sep%fn2(2) = 5.0000000_dp; sep%fn3(2) = 1.8000000_dp
        sep%fn1(3) = -0.0183360_dp; sep%fn2(3) = 2.0000000_dp; sep%fn3(3) = 2.1000000_dp
      CASE (2)
      CASE (3)
        !Lithium
        sep%defined          = .TRUE.
        sep%alp   =   1.2501400_dp
        sep%beta(0) =  -1.3500400_dp
        sep%beta(1) =  -1.3500400_dp
        sep%sto_exponents(0)    =   0.7023800_dp
        sep%sto_exponents(1)   =   0.7023800_dp
        sep%uss   =  -5.1280000_dp
        sep%upp   =  -2.7212000_dp
        sep%gss   =   7.3000000_dp
        sep%gsp  =   5.4200000_dp
        sep%gpp  =   5.0000000_dp
        sep%gp2  =   4.5200000_dp
        sep%hsp  =   0.8300000_dp
      CASE (4)
      !Beryllium
        sep%defined          = .TRUE.
        sep%alp   =   1.6694340_dp
        sep%beta(0) =  -4.0170960_dp
        sep%beta(1)=  -4.0170960_dp
        sep%sto_exponents(0)    =   1.0042100_dp
        sep%sto_exponents(1)   =   1.0042100_dp
        sep%uss   = -16.6023780_dp
        sep%upp   = -10.7037710_dp
        sep%gss   =   9.0000000_dp
        sep%gsp  =   7.4300000_dp
        sep%gpp  =   6.9700000_dp
        sep%gp2  =   6.2200000_dp
        sep%hsp  =   1.2800000_dp
      CASE (5)
        !Boron
        sep%defined = .TRUE.
        sep%alp   =   2.4469090_dp
        sep%beta(0) =  -9.5991140_dp
        sep%beta(1)=  -6.2737570_dp
        sep%sto_exponents(0)    =   1.6117090_dp
        sep%sto_exponents(1)   =   1.5553850_dp
        sep%uss   = -34.4928700_dp
        sep%upp   = -22.6315250_dp
        sep%gss   =  10.5900000_dp
        sep%gsp  =   9.5600000_dp
        sep%gpp  =   8.8600000_dp
        sep%gp2  =   7.8600000_dp
        sep%hsp  =   1.8100000_dp
        sep%bfn1(1,1) =  0.1826130_dp; sep%bfn2(1,1) = 6.0000000_dp; sep%bfn3(1,1) = 0.7275920_dp
        sep%bfn1(2,1) =  0.1185870_dp; sep%bfn2(2,1) = 6.0000000_dp; sep%bfn3(2,1) = 1.4666390_dp
        sep%bfn1(3,1) = -0.0732800_dp; sep%bfn2(3,1) = 5.0000000_dp; sep%bfn3(3,1) = 1.5709750_dp
        sep%bfn1(1,2) =  0.4122530_dp; sep%bfn2(1,2) =10.0000000_dp; sep%bfn3(1,2) = 0.8325860_dp
        sep%bfn1(2,2) = -0.1499170_dp; sep%bfn2(2,2) = 6.0000000_dp; sep%bfn3(2,2) = 1.1862200_dp
        sep%bfn1(3,2) =  0.0000000_dp; sep%bfn2(3,2) = 0.0000000_dp; sep%bfn3(3,2) = 0.0000000_dp
        sep%bfn1(1,3) =  0.2617510_dp; sep%bfn2(1,3) = 8.0000000_dp; sep%bfn3(1,3) = 1.0639950_dp
        sep%bfn1(2,3) =  0.0502750_dp; sep%bfn2(2,3) = 5.0000000_dp; sep%bfn3(2,3) = 1.9364920_dp
        sep%bfn1(3,3) =  0.0000000_dp; sep%bfn2(3,3) = 0.0000000_dp; sep%bfn3(3,3) = 0.0000000_dp
        sep%bfn1(1,4) =  0.3592440_dp; sep%bfn2(1,4) = 9.0000000_dp; sep%bfn3(1,4) = 0.8193510_dp
        sep%bfn1(2,4) =  0.0747290_dp; sep%bfn2(2,4) = 9.0000000_dp; sep%bfn3(2,4) = 1.5744140_dp
        sep%bfn1(3,4) =  0.0000000_dp; sep%bfn2(3,4) = 0.0000000_dp; sep%bfn3(3,4) = 0.0000000_dp
      CASE (6)
        !Carbon
        sep%defined = .TRUE.
        sep%alp   =    2.6482740_dp
        sep%beta(0) =  -15.7157830_dp
        sep%beta(1) =   -7.7192830_dp
        sep%sto_exponents(0)     =    1.8086650_dp
        sep%sto_exponents(1)     =    1.6851160_dp
        sep%uss   =  -52.0286580_dp
        sep%upp   =  -39.6142390_dp
        sep%gss   =   12.2300000_dp
        sep%gsp   =   11.4700000_dp
        sep%gpp   =   11.0800000_dp
        sep%gp2   =    9.8400000_dp
        sep%hsp   =    2.4300000_dp
        sep%fn1(1) =  0.0113550_dp; sep%fn2(1) =  5.0000000_dp; sep%fn3(1) =  1.6000000_dp
        sep%fn1(2) =  0.0459240_dp; sep%fn2(2) =  5.0000000_dp; sep%fn3(2) =  1.8500000_dp
        sep%fn1(3) = -0.0200610_dp; sep%fn2(3) =  5.0000000_dp; sep%fn3(3) =  2.0500000_dp
        sep%fn1(4) = -0.0012600_dp; sep%fn2(4) =  5.0000000_dp; sep%fn3(4) =  2.6500000_dp
      CASE (7)
      !Nitrogen
        sep%defined = .TRUE.
        sep%alp   =    2.9472860_dp
        sep%beta(0) =  -20.2991100_dp
        sep%beta(1) =  -18.2386660_dp
        sep%sto_exponents(0)   =    2.3154100_dp
        sep%sto_exponents(1)    =    2.1579400_dp
        sep%uss   =  -71.8600000_dp
        sep%upp   =  -57.1675810_dp
        sep%gss   =   13.5900000_dp
        sep%gsp   =   12.6600000_dp
        sep%gpp   =   12.9800000_dp
        sep%gp2   =   11.5900000_dp
        sep%hsp   =    3.1400000_dp
        sep%fn1(1) =  0.0252510_dp; sep%fn2(1) =  5.0000000_dp; sep%fn3(1) =  1.5000000_dp
        sep%fn1(2) =  0.0289530_dp; sep%fn2(2) =  5.0000000_dp; sep%fn3(2) =  2.1000000_dp
        sep%fn1(3) = -0.0058060_dp; sep%fn2(3) =  2.0000000_dp; sep%fn3(3) =  2.4000000_dp
      CASE (8)
      !Oxygen
        sep%defined = .TRUE.
        sep%alp   =    4.4553710_dp
        sep%beta(0) =  -29.2727730_dp
        sep%beta(1) =  -29.2727730_dp
        sep%sto_exponents(0)   =    3.1080320_dp
        sep%sto_exponents(1)    =    2.5240390_dp
        sep%uss   =  -97.8300000_dp
        sep%upp   =  -78.2623800_dp
        sep%gss   =   15.4200000_dp
        sep%gsp   =   14.4800000_dp
        sep%gpp   =   14.5200000_dp
        sep%gp2   =   12.9800000_dp
        sep%hsp   =    3.9400000_dp
        sep%fn1(1) = 0.2809620_dp; sep%fn2(1) = 5.0000000_dp; sep%fn3(1) = 0.8479180_dp
        sep%fn1(2) = 0.0814300_dp; sep%fn2(2) = 7.0000000_dp; sep%fn3(2) = 1.4450710_dp
      CASE (9)
      !Fluorine
        sep%defined = .TRUE.
        sep%alp   =    5.5178000_dp
        sep%beta(0) =  -69.5902770_dp
        sep%beta(1) =  -27.9223600_dp
        sep%sto_exponents(0)   =    3.7700820_dp
        sep%sto_exponents(1)    =    2.4946700_dp
        sep%uss   = -136.1055790_dp
        sep%upp   = -104.8898850_dp
        sep%gss   =   16.9200000_dp
        sep%gsp   =   17.2500000_dp
        sep%gpp   =   16.7100000_dp
        sep%gp2   =   14.9100000_dp
        sep%hsp   =    4.8300000_dp
        sep%fn1(1) = 0.2420790_dp; sep%fn2(1) = 4.8000000_dp; sep%fn3(1) = 0.9300000_dp
        sep%fn1(2) = 0.0036070_dp; sep%fn2(2) = 4.6000000_dp; sep%fn3(2) = 1.6600000_dp
      CASE (10:11)
      CASE (12)
        !Magnesium
        sep%defined = .TRUE.
        sep%uss   =  -14.96959313_dp
        sep%upp   =  -11.56229248_dp
        sep%beta(0) =   -1.25974355_dp
        sep%beta(1)=   -0.77836604_dp
        sep%sto_exponents(0)    =    1.22339270_dp
        sep%sto_exponents(1)   =    1.02030798_dp
        sep%alp   =    1.67049799_dp
        sep%gss   =    7.50132277_dp
        sep%gsp  =    6.34591536_dp
        sep%gpp  =    4.77534467_dp
        sep%gp2  =    4.34017279_dp
        sep%hsp  =    0.48930466_dp
        sep%fn1(1) =  2.55017735_dp; sep%fn2(1) = 4.29397225_dp; sep%fn3(1) = 0.79989601_dp
        sep%fn1(2) = -0.00565806_dp; sep%fn2(2) = 2.96053910_dp; sep%fn3(2) = 1.47499983_dp
        sep%fn1(3) = -0.00610286_dp; sep%fn2(3) = 2.61416919_dp; sep%fn3(3) = 2.42604040_dp
      CASE (13)
       !Aluminium
        sep%defined = .TRUE.
        sep%alp   =   1.9765860_dp
        sep%beta(0) =  -3.8668220_dp
        sep%beta(1)=  -2.3171460_dp
        sep%sto_exponents(0)   =   1.5165930_dp
        sep%sto_exponents(1)   =   1.3063470_dp
        sep%sto_exponents(2)   =   1.0000000_dp
        sep%uss   = -24.3535850_dp
        sep%upp   = -18.3636450_dp
        sep%gss   =   8.0900000_dp
        sep%gsp  =   6.6300000_dp
        sep%gpp  =   5.9800000_dp
        sep%gp2  =   5.4000000_dp
        sep%hsp  =   0.7000000_dp
        sep%fn1(1) = 0.0900000_dp; sep%fn2(1) = 12.3924430_dp; sep%fn3(1) = 2.0503940_dp
      CASE (14)
      !Silicon
        sep%defined = .TRUE.
        sep%alp   =   2.257816_dp
        sep%beta(0) =  -3.784852_dp
        sep%beta(1)=  -1.968123_dp
        sep%sto_exponents(0)   =   1.830697_dp
        sep%sto_exponents(1)   =   1.2849530_dp
        sep%sto_exponents(2)   =   1.0000000_dp
        sep%uss   = -33.9536220_dp
        sep%upp   = -28.9347490_dp
        sep%gss   =   9.8200000_dp
        sep%gsp  =   8.3600000_dp
        sep%gpp  =   7.3100000_dp
        sep%gp2  =   6.5400000_dp
        sep%hsp  =   1.3200000_dp
        sep%fn1(1) =  0.25000000_dp; sep%fn2(1) = 9.00000000_dp; sep%fn3(1) = 0.91145300_dp
        sep%fn1(2) =  0.06151300_dp; sep%fn2(2) = 5.00000000_dp; sep%fn3(2) = 1.99556900_dp
        sep%fn1(3) =  0.00207890_dp; sep%fn2(3) = 5.00000000_dp; sep%fn3(3) = 2.99061000_dp
      CASE (15)
      !Phosphorus
        sep%defined = .TRUE.
        sep%alp = 2.4553220_dp
        sep%beta(0) =   -6.3537640_dp
        sep%beta(1)=   -6.5907090_dp
        sep%sto_exponents(0)   =    1.9812800_dp
        sep%sto_exponents(1)   =    1.8751500_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%uss   =  -42.0298630_dp
        sep%upp   =  -34.0307090_dp
        sep%gss   =   11.5600050_dp
        sep%gsp  =   5.2374490_dp
        sep%gpp  =    7.8775890_dp
        sep%gp2  =    7.3076480_dp
        sep%hsp  =    0.7792380_dp
        sep%fn1(1) = -0.03182700_dp; sep%fn2(1) = 6.00000000_dp; sep%fn3(1) = 1.47432300_dp
        sep%fn1(2) =  0.01847000_dp; sep%fn2(2) = 7.00000000_dp; sep%fn3(2) = 0.03329000_dp
        sep%fn1(3) =  0.03329000_dp; sep%fn2(3) = 9.00000000_dp; sep%fn3(3) = 3.00657600_dp
      CASE (16)
        !Sulphur
        sep%defined = .TRUE.
        sep%alp   =    2.4616480_dp
        sep%beta(0) =  -3.9205660_dp
        sep%beta(1)=  -7.9052780_dp
        sep%sto_exponents(0)   =    2.3665150_dp
        sep%sto_exponents(1)   =    1.6672630_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%uss   =  -56.6940560_dp
        sep%upp   =  -48.7170490_dp
        sep%gss   =   11.7863290_dp
        sep%gsp  =    8.6631270_dp
        sep%gpp  =   10.0393080_dp
        sep%gp2  =    7.7816880_dp
        sep%hsp  =    2.5321370_dp
        sep%fn1(1) = -0.50919500_dp; sep%fn2(1) = 4.59369100_dp; sep%fn3(1) = 0.77066500_dp
        sep%fn1(2) = -0.01186300_dp; sep%fn2(2) = 5.86573100_dp; sep%fn3(2) = 1.50331300_dp
        sep%fn1(3) =  0.01233400_dp; sep%fn2(3) =13.55733600_dp; sep%fn3(3) = 2.00917300_dp
      CASE (17)
        !Chlorine
        sep%defined = .TRUE.
        sep%alp   =    2.9193680_dp
        sep%beta(0) =  -24.5946700_dp
        sep%beta(1)=  -14.6372160_dp
        sep%sto_exponents(0)    =    3.6313760_dp
        sep%sto_exponents(1)   =    2.0767990_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%uss   = -111.6139480_dp
        sep%upp   =  -76.6401070_dp
        sep%gss   =   15.0300000_dp
        sep%gsp  =   13.1600000_dp
        sep%gpp  =   11.3000000_dp
        sep%gp2  =    9.9700000_dp
        sep%hsp  =    2.4200000_dp
        sep%fn1(1) = 0.0942430_dp; sep%fn2(1) = 4.0000000_dp; sep%fn3(1) = 1.3000000_dp
        sep%fn1(2) = 0.0271680_dp; sep%fn2(2) = 4.0000000_dp; sep%fn3(2) = 2.1000000_dp
      CASE (18:29)
      CASE (30)
        !Zinc
        sep%defined = .TRUE.
        sep%uss    =    -21.0400080_dp
        sep%upp    =    -17.6555740_dp
        sep%beta(0)  =     -1.9974290_dp
        sep%beta(1) =     -4.7581190_dp
        sep%sto_exponents(0)     =      1.9542990_dp
        sep%sto_exponents(1)    =      1.3723650_dp
        sep%sto_exponents(2)    =      1.0000000_dp
        sep%alp    =      1.4845630_dp
        sep%gss    =     11.8000000_dp
        sep%gsp   =     11.1820180_dp
        sep%gpp   =     13.3000000_dp
        sep%gp2   =     12.9305200_dp
        sep%hsp   =      0.4846060_dp
      CASE (31)
      CASE (32)
        !Germanium
        sep%defined = .TRUE.
        sep%uss    =    -34.1838890_dp
        sep%upp    =    -28.6408110_dp
        sep%beta(0)  =     -4.3566070_dp
        sep%beta(1) =     -0.9910910_dp
        sep%sto_exponents(0)     =      1.2196310_dp
        sep%sto_exponents(1)    =      1.9827940_dp
        sep%sto_exponents(2)    =      1.0000000_dp
        sep%alp    =      2.1364050_dp
        sep%gss    =     10.1686050_dp
        sep%gsp   =      8.1444730_dp
        sep%gpp   =      6.6719020_dp
        sep%gp2   =      6.2697060_dp
        sep%hsp   =      0.9370930_dp
      CASE (33:34)
!see calpar.F90 for additional constraint params
!                    DATA FOR ELEMENT 33      ARSENIC
!     data ussam1(33)/  - 41.6817510e0_dp/
!     data uppam1(33)/  - 33.4506152e0_dp/
!     data betasa(33)/  - 5.6481504e0_dp/
!     data betapa(33)/  - 4.9979109e0_dp/
!     data zsam1(33)/ 2.2576897e0_dp/
!     data zpam1(33)/ 1.7249710e0_dp/
!     data alpam1(33)/ 2.2405380e0_dp/
!     data gssam1(33)/ 11.0962258e0_dp/
!     data gspam1(33)/ 4.9259328e0_dp/
!     data gppam1(33)/ 7.8781648e0_dp/
!     data gp2am1(33)/ 7.5961088e0_dp/
!     data hspam1(33)/ 0.6246173e0_dp/
!     data guesa1(33,1)/  - 0.0073614e0_dp/
!     data guesa2(33,1)/ 4.9433993e0_dp/
!     data guesa3(33,1)/ 1.4544264e0_dp/
!     data guesa1(33,2)/ 0.0437629e0_dp/
!     data guesa2(33,2)/ 3.1944613e0_dp/
!     data guesa3(33,2)/ 2.0144939e0_dp/
!                    DATA FOR ELEMENT 34        SELENIUM
!     data ussam1(34)/  - 41.9984056e0_dp/
!     data uppam1(34)/  - 32.8575485e0_dp/
!     data betasa(34)/  - 3.1470826e0_dp/
!     data betapa(34)/  - 6.1468406e0_dp/
!     data zsam1(34)/ 2.6841570e0_dp/
!     data zpam1(34)/ 2.0506164e0_dp/
!     data alpam1(34)/ 2.6375694e0_dp/
!     data gssam1(34)/ 6.7908891e0_dp/
!     data gspam1(34)/ 6.4812786e0_dp/
!     data gppam1(34)/ 6.4769273e0_dp/
!     data gp2am1(34)/ 5.2796993e0_dp/
!     data hspam1(34)/ 4.4548356e0_dp/
!     data guesa1(34,1)/ 0.1116681e0_dp/
!     data guesa2(34,1)/ 6.5086644e0_dp/
!     data guesa3(34,1)/ 1.4981077e0_dp/
!     data guesa1(34,2)/ 0.0396143e0_dp/
!     data guesa2(34,2)/ 6.5241228e0_dp/
!     data guesa3(34,2)/ 2.0751916e0_dp/
!
      CASE (35)
        !Bromine
        sep%defined = .TRUE.
        sep%alp   =    2.5765460_dp
        sep%beta(0) =  -19.3998800_dp
        sep%beta(1)=   -8.9571950_dp
        sep%sto_exponents(0)    =    3.0641330_dp
        sep%sto_exponents(1)   =    2.0383330_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%uss   = -104.6560630_dp
        sep%upp   =  -74.9300520_dp
        sep%gss   =   15.0364395_dp
        sep%gsp  =   13.0346824_dp
        sep%gpp  =   11.2763254_dp
        sep%gp2  =    9.8544255_dp
        sep%hsp  =    2.4558683_dp
        sep%fn1(1) = 0.0666850_dp; sep%fn2(1) = 4.0000000_dp; sep%fn3(1) = 1.5000000_dp
        sep%fn1(2) = 0.0255680_dp; sep%fn2(2) = 4.0000000_dp; sep%fn3(2) = 2.3000000_dp
      CASE (36:49)
      CASE (50)
        !Tin
        sep%defined = .TRUE.
        sep%alp   =    1.8369360_dp
        sep%beta(0) =   -3.2350000_dp
        sep%beta(1)=   -2.5778900_dp
        sep%sto_exponents(0)    =    2.5993760_dp
        sep%sto_exponents(1)   =    1.6959620_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%uss   =  -35.4967410_dp
        sep%upp   =  -28.0976360_dp
        sep%gss   =    9.8000000_dp
        sep%gsp  =    8.3000000_dp
        sep%gpp  =    7.3000000_dp
        sep%gp2  =    6.5000000_dp
        sep%hsp  =    1.3000000_dp
      CASE (51:52)
!                    DATA FOR ELEMENT 51        ANTIMONY
!     data ussam1(51)/  - 44.4381620e0_dp/
!     data uppam1(51)/  - 32.3895140e0_dp/
!     data betasa(51)/  - 7.3823300e0_dp/
!     data betapa(51)/  - 3.6331190e0_dp/
!     data zsam1(51)/ 2.2548230e0_dp/
!     data zpam1(51)/ 2.2185920e0_dp/
!     data alpam1(51)/ 2.2763310e0_dp/
!     data gssam1(51)/ 11.4302510e0_dp/
!     data gspam1(51)/ 5.7879220e0_dp/
!     data gppam1(51)/ 6.4240940e0_dp/
!     data gp2am1(51)/ 6.8491810e0_dp/
!     data hspam1(51)/ 0.5883400e0_dp/
!     data guesa1(51,1)/  - 0.5964470e0_dp/
!     data guesa2(51,1)/ 6.0279500e0_dp/
!     data guesa3(51,1)/ 1.7103670e0_dp/
!     data guesa1(51,2)/ 0.8955130e0_dp/
!     data guesa2(51,2)/ 3.0281090e0_dp/
!     data guesa3(51,2)/ 1.5383180e0_dp/
!                    DATA FOR ELEMENT 52        TELLURIUM
!     data ussam1(52)/  - 39.2454230e0_dp/
!     data uppam1(52)/  - 30.8515845e0_dp/
!     data betasa(52)/  - 8.3897294e0_dp/
!     data betapa(52)/  - 5.1065429e0_dp/
!     data zsam1(52)/ 2.1321165e0_dp/
!     data zpam1(52)/ 1.9712680e0_dp/
!     data alpam1(52)/ 6.0171167e0_dp/
!     data gssam1(52)/ 4.9925231e0_dp/
!     data gspam1(52)/ 4.9721484e0_dp/
!     data gppam1(52)/ 7.2097852e0_dp/
!     data gp2am1(52)/ 5.6211521e0_dp/
!     data hspam1(52)/ 4.0071821e0_dp/
!     data guesa1(52,1)/ 0.4873378e0_dp/
!     data guesa2(52,1)/ 6.0519413e0_dp/
!     data guesa3(52,1)/ 1.3079857e0_dp/
!     data guesa1(52,2)/ 0.1520464e0_dp/
!     data guesa2(52,2)/ 3.8304067e0_dp/
!     data guesa3(52,2)/ 2.0899707e0_dp/

      CASE (53)
        !Iodine
        sep%defined = .TRUE.
        sep%uss   = -103.5896630_dp
        sep%upp   =  -74.4299970_dp
        sep%beta(0) =   -8.4433270_dp
        sep%beta(1)=   -6.3234050_dp
        sep%sto_exponents(0)    =    2.1028580_dp
        sep%sto_exponents(1)   =    2.1611530_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    2.2994240_dp
        sep%gss   =   15.0404486_dp
        sep%gsp  =   13.0565580_dp
        sep%gpp  =   11.1477837_dp
        sep%gp2  =    9.9140907_dp
        sep%hsp  =    2.4563820_dp
        sep%fn1(1) = 0.0043610_dp; sep%fn2(1) = 2.3000000_dp; sep%fn3(1) = 1.8000000_dp
        sep%fn1(2) = 0.0157060_dp; sep%fn2(2) = 3.0000000_dp; sep%fn3(2) = 2.2400000_dp
      CASE (80)
        !Mercury
        sep%defined = .TRUE.
        sep%uss   = -19.9415780_dp
        sep%upp   =  -11.1108700_dp
        sep%beta(0) =   -0.9086570_dp
        sep%beta(1)=   -4.9093840_dp
        sep%sto_exponents(0)    =    2.0364130_dp
        sep%sto_exponents(1)   =    1.9557660_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    1.4847340_dp
        sep%gss   =   10.8000000_dp
        sep%gsp  =    9.3000000_dp
        sep%gpp  =   14.3000000_dp
        sep%gp2  =   13.5000000_dp
        sep%hsp  =    1.3000000_dp
      CASE (81)
      CASE (82)
        !Lead
        sep%defined = .TRUE.
        sep%uss   = -47.3196920_dp
        sep%upp   =  -28.8475600_dp
        sep%beta(0) =   -8.0423870_dp
        sep%beta(1)=   -3.0000000_dp
        sep%sto_exponents(0)    =    2.4982860_dp
        sep%sto_exponents(1)   =    1.7283330_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    1.7283330_dp
        sep%gss   =    9.8000000_dp
        sep%gsp  =    8.3000000_dp
        sep%gpp  =    7.3000000_dp
        sep%gp2  =    6.5000000_dp
        sep%hsp  =    1.3000000_dp
      CASE (83:106)
    END SELECT
    CALL calpar(z,sep)

    sep%beta = sep%beta/evolt
    sep%uss = sep%uss/evolt; sep%upp = sep%upp/evolt; sep%udd = sep%udd/evolt
    sep%alp = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss = sep%gss/evolt; sep%gsp = sep%gsp/evolt; sep%gpp = sep%gpp/evolt; sep%gp2 = sep%gp2/evolt 
    sep%gsd = sep%gsd/evolt; sep%gpd = sep%gpd/evolt; sep%gdd = sep%gdd/evolt
    sep%hsp = sep%hsp/evolt
    sep%fn1 = sep%fn1*bohr/evolt; sep%fn2 = sep%fn2/bohr/bohr; sep%fn3 = sep%fn3*bohr
    sep%bfn1 = sep%bfn1*bohr/evolt; sep%bfn2 = sep%bfn2/bohr/bohr; sep%bfn3 = sep%bfn3*bohr

  END SUBROUTINE am1_default_parameter

! *****************************************************************************

  SUBROUTINE rm1_default_parameter ( sep, z )
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z

    CHARACTER(len=*), PARAMETER :: routineN = 'rm1_default_parameter', &
      routineP = moduleN//':'//routineN

! -----------------------------------------------------------------------------

    sep%typ = do_method_rm1
    sep%name = ptable(z)%name
    sep%eheat = ptable(z)%heat_of_formation / kcalmol
    sep%z = z
    CALL valence_electrons ( z, sep%zeff, sep%natorb )
    sep%ngauss = 6
    sep%defined = .FALSE.
    sep%sto_exponents(0:3) = 0._dp
    sep%beta(0:3) = 0._dp
    sep%uss = 0._dp; sep%upp = 0._dp; sep%udd = 0._dp
    sep%alp = 0._dp
    sep%eisol = 0._dp
    sep%gss = 0._dp; sep%gsp = 0._dp; sep%gpp = 0._dp; sep%gp2 = 0._dp
    sep%gsd = 0._dp; sep%gpd = 0._dp; sep%gdd = 0._dp
    sep%hsp = 0._dp
    sep%dd = 0._dp; sep%qq = 0._dp; sep%am = 0._dp; sep%ad = 0._dp; sep%aq = 0._dp
    sep%fn1 = 0._dp; sep%fn2 = 0._dp; sep%fn3 = 0._dp
    sep%bfn1 = 0._dp; sep%bfn2 = 0._dp; sep%bfn3 = 0._dp

    SELECT CASE (z)
      CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "invalid atomic number")
      CASE (0)
        sep%defined = .TRUE.
      CASE (1)
         !Hydrogen
         sep%defined = .TRUE.
         sep%alp     =   3.06835947_dp
         sep%beta(0) =  -5.76544469_dp
         sep%sto_exponents(0)    =   1.08267366_dp
         sep%uss    = -11.96067697_dp
         sep%gss    =  13.98321296_dp
         sep%fn1(1) =   0.10288875_dp; sep%fn2(1) = 5.90172268_dp; sep%fn3(1) = 1.17501185_dp
         sep%fn1(2) =   0.06457449_dp; sep%fn2(2) = 6.41785671_dp; sep%fn3(2) = 1.93844484_dp
         sep%fn1(3) =  -0.03567387_dp; sep%fn2(3) = 2.80473127_dp; sep%fn3(3) = 1.63655241_dp
      CASE (2:5)
      CASE (6)
         !Carbon
         sep%defined = .TRUE.
         sep%alp     =    2.79282078_dp
         sep%beta(0) =  -15.45932428_dp
         sep%beta(1) =   -8.23608638_dp
         sep%sto_exponents(0)     =    1.85018803_dp
         sep%sto_exponents(1)     =    1.76830093_dp
         sep%uss   =  -51.72556032_dp
         sep%upp   =  -39.40728943_dp
         sep%gss   =   13.05312440_dp
         sep%gsp   =   11.33479389_dp
         sep%gpp   =   10.95113739_dp
         sep%gp2   =    9.72395099_dp
         sep%hsp   =    1.55215133_dp
         sep%fn1(1) =  0.07462271_dp; sep%fn2(1) =  5.73921605_dp; sep%fn3(1) =  1.04396983_dp
         sep%fn1(2) =  0.01177053_dp; sep%fn2(2) =  6.92401726_dp; sep%fn3(2) =  1.66159571_dp
         sep%fn1(3) =  0.03720662_dp; sep%fn2(3) =  6.26158944_dp; sep%fn3(3) =  1.63158721_dp
         sep%fn1(4) = -0.00270657_dp; sep%fn2(4) =  9.00003735_dp; sep%fn3(4) =  2.79557901_dp
      CASE (7)
         !Nitrogen
         sep%defined = .TRUE.
         sep%alp     =    2.96422542_dp
         sep%beta(0) =  -20.87124548_dp
         sep%beta(1) =  -16.67171853_dp
         sep%sto_exponents(0)   =    2.37447159_dp
         sep%sto_exponents(1)   =    1.97812569_dp
         sep%uss   =  -70.85123715_dp
         sep%upp   =  -57.97730920_dp
         sep%gss   =   13.08736234_dp
         sep%gsp   =   13.21226834_dp
         sep%gpp   =   13.69924324_dp
         sep%gp2   =   11.94103953_dp
         sep%hsp   =    5.00000846_dp
         sep%fn1(1)=  0.06073380_dp; sep%fn2(1)=  4.58892946_dp; sep%fn3(1)=  1.37873881_dp
         sep%fn1(2)=  0.02438558_dp; sep%fn2(2)=  4.62730519_dp; sep%fn3(2)=  2.08370698_dp
         sep%fn1(3)= -0.02283430_dp; sep%fn2(3)=  2.05274659_dp; sep%fn3(3)=  1.86763816_dp
      CASE (8)
         !Oxygen
         sep%defined = .TRUE.
         sep%alp     =    4.17196717_dp
         sep%beta(0) =  -29.85101212_dp
         sep%beta(1) =  -29.15101314_dp
         sep%sto_exponents(0)   =    3.17936914_dp
         sep%sto_exponents(1)   =    2.55361907_dp
         sep%uss   =  -96.94948069_dp
         sep%upp   =  -77.89092978_dp
         sep%gss   =   14.00242788_dp
         sep%gsp   =   14.95625043_dp
         sep%gpp   =   14.14515138_dp
         sep%gp2   =   12.70325497_dp
         sep%hsp   =    3.93217161_dp
         sep%fn1(1) = 0.23093552_dp; sep%fn2(1) = 5.21828736_dp; sep%fn3(1) = 0.90363555_dp
         sep%fn1(2) = 0.05859873_dp; sep%fn2(2) = 7.42932932_dp; sep%fn3(2) = 1.51754610_dp
      CASE (9)
         !Fluorine
         sep%defined = .TRUE.
         sep%alp     =    6.00000062_dp
         sep%beta(0) =  -70.00000512_dp
         sep%beta(1) =  -32.67982711_dp
         sep%sto_exponents(0)   =    4.40337913_dp
         sep%sto_exponents(1)   =    2.64841556_dp
         sep%uss   = -134.18369591_dp
         sep%upp   = -107.84660920_dp
         sep%gss   =   16.72091319_dp
         sep%gsp   =   16.76142629_dp
         sep%gpp   =   15.22581028_dp
         sep%gp2   =   14.86578679_dp
         sep%hsp   =    1.99766171_dp
         sep%fn1(1) = 0.40302025_dp; sep%fn2(1) = 7.20441959_dp; sep%fn3(1) = 0.81653013_dp
         sep%fn1(2) = 0.07085831_dp; sep%fn2(2) = 9.00001562_dp; sep%fn3(2) = 1.43802381_dp
      CASE (10:14)
      CASE (15)
         !Phosphorus
         sep%defined = .TRUE.
         sep%alp     =    1.90993294_dp
         sep%beta(0) =   -6.13514969_dp
         sep%beta(1) =   -5.94442127_dp
         sep%sto_exponents(0)   =    2.12240118_dp
         sep%sto_exponents(1)   =    1.74327954_dp 
         sep%sto_exponents(2)   =    1.00000000_dp
         sep%uss   =  -41.81533184_dp
         sep%upp   =  -34.38342529_dp 
         sep%gss   =   11.08059265_dp
         sep%gsp   =    5.68339201_dp
         sep%gpp   =    7.60417563_dp
         sep%gp2   =    7.40265182_dp
         sep%hsp   =    1.16181792_dp
         sep%fn1(1) = -0.41063467_dp; sep%fn2(1) = 6.08752832_dp; sep%fn3(1) = 1.31650261_dp
         sep%fn1(2) = -0.16299288_dp; sep%fn2(2) = 7.09472602_dp; sep%fn3(2) = 1.90721319_dp
         sep%fn1(3) = -0.04887125_dp; sep%fn2(3) = 8.99979308_dp; sep%fn3(3) = 2.65857780_dp
      CASE (16)
         !Sulphur
         sep%defined = .TRUE.
         sep%alp     =    2.44015636_dp
         sep%beta(0) =  -1.95910719_dp
         sep%beta(1) =  -8.77430652_dp
         sep%sto_exponents(0)   =    2.13344308_dp
         sep%sto_exponents(1)   =    1.87460650_dp
         sep%sto_exponents(2)   =    1.00000000_dp
         sep%uss   =  -55.16775121_dp
         sep%upp   =  -46.52930422_dp
         sep%gss   =   12.48828408_dp
         sep%gsp   =    8.56910574_dp
         sep%gpp   =    8.52301167_dp
         sep%gp2   =    7.66863296_dp
         sep%hsp   =    3.88978932_dp
         sep%fn1(1) = -0.74601055_dp; sep%fn2(1) = 4.81038002_dp; sep%fn3(1) = 0.59380129_dp
         sep%fn1(2) = -0.06519286_dp; sep%fn2(2) = 7.20760864_dp; sep%fn3(2) = 1.29492008_dp
         sep%fn1(3) = -0.00655977_dp; sep%fn2(3) = 9.00000180_dp; sep%fn3(3) = 1.80060151_dp
      CASE (17)
         !Chlorine
         sep%defined = .TRUE.
         sep%alp     =    3.69358828_dp
         sep%beta(0) =  -19.92430432_dp
         sep%beta(1) =  -11.52935197_dp
         sep%sto_exponents(0)    =   3.86491071_dp
         sep%sto_exponents(1)   =    1.89593144_dp
         sep%sto_exponents(2)   =    1.00000000_dp
         sep%uss   = -118.47306918_dp
         sep%upp   =  -76.35330340_dp
         sep%gss   =   15.36023105_dp
         sep%gsp   =   13.30671171_dp
         sep%gpp   =   12.56502640_dp
         sep%gp2   =    9.66397083_dp
         sep%hsp   =    1.76489897_dp
         sep%fn1(1) = 0.12947108_dp; sep%fn2(1) = 2.97724424_dp; sep%fn3(1) = 1.46749784_dp
         sep%fn1(2) = 0.00288899_dp; sep%fn2(2) = 7.09827589_dp; sep%fn3(2) = 2.50002723_dp
      CASE (18:34)
      CASE (35)
         !Bromine
         sep%defined = .TRUE.
         sep%alp     =    2.86710532_dp
         sep%beta(0) =  -1.34139841_dp
         sep%beta(1) =  -8.20225991_dp
         sep%sto_exponents(0)   =    5.73157215_dp
         sep%sto_exponents(1)   =    2.03147582_dp
         sep%sto_exponents(2)   =    1.00000000_dp
         sep%uss   = -113.48398183_dp
         sep%upp   =  -76.18720023_dp
         sep%gss   =   17.11563074_dp
         sep%gsp   =   15.62419251_dp
         sep%gpp   =   10.73546293_dp
         sep%gp2   =    8.86056199_dp
         sep%hsp   =    2.23512762_dp
         sep%fn1(1) =  0.98689937_dp; sep%fn2(1) =  4.28484191_dp; sep%fn3(1) = 2.00019696_dp
         sep%fn1(2) = -0.92731247_dp; sep%fn2(2) =  4.54005910_dp; sep%fn3(2) = 2.01617695_dp
      CASE (36:49)
      CASE (50)
      CASE (51:52)
      CASE (53)
         !Iodine
         sep%defined = .TRUE.
         sep%alp     =    2.14157092_dp
         sep%beta(0) =   -4.19316149_dp
         sep%beta(1) =   -4.40038412_dp
         sep%sto_exponents(0)   =    2.53003753_dp
         sep%sto_exponents(1)   =    2.31738678_dp
         sep%sto_exponents(2)   =    1.0000000_dp
         sep%uss   =  -74.89997837_dp
         sep%upp   =  -51.41023805_dp
         sep%gss   =   19.99974131_dp
         sep%gsp   =    7.68957672_dp
         sep%gpp   =    7.30488343_dp
         sep%gp2   =    6.85424614_dp
         sep%hsp   =    1.41602940_dp
         sep%fn1(1) = -0.08147724_dp; sep%fn2(1) = 1.56065072_dp; sep%fn3(1) = 2.00002063_dp
         sep%fn1(2) =  0.05914991_dp; sep%fn2(2) = 5.76111270_dp; sep%fn3(2) = 2.20488800_dp
      CASE (80:106)
    END SELECT
    CALL calpar(z,sep)

    sep%beta = sep%beta/evolt
    sep%uss = sep%uss/evolt; sep%upp = sep%upp/evolt; sep%udd = sep%udd/evolt
    sep%alp = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss = sep%gss/evolt; sep%gsp = sep%gsp/evolt; sep%gpp = sep%gpp/evolt; sep%gp2 = sep%gp2/evolt 
    sep%gsd = sep%gsd/evolt; sep%gpd = sep%gpd/evolt; sep%gdd = sep%gdd/evolt
    sep%hsp = sep%hsp/evolt
    sep%fn1 = sep%fn1*bohr/evolt; sep%fn2 = sep%fn2/bohr/bohr; sep%fn3 = sep%fn3*bohr
    sep%bfn1 = sep%bfn1*bohr/evolt; sep%bfn2 = sep%bfn2/bohr/bohr; sep%bfn3 = sep%bfn3*bohr

  END SUBROUTINE rm1_default_parameter

! *****************************************************************************

  SUBROUTINE pm3_default_parameter ( sep, z )
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z

    CHARACTER(len=*), PARAMETER :: routineN = 'pm3_default_parameter', &
      routineP = moduleN//':'//routineN

! -----------------------------------------------------------------------------

    sep%typ = do_method_pm3
    sep%name = ptable(z)%name
    sep%eheat = ptable(z)%heat_of_formation / kcalmol
    sep%z = z
    CALL valence_electrons ( z, sep%zeff, sep%natorb )
    sep%ngauss = 6
    sep%defined = .FALSE.
    sep%sto_exponents(0:3) = 0._dp
    sep%beta(0:3) = 0._dp
    sep%uss = 0._dp; sep%upp = 0._dp; sep%udd = 0._dp
    sep%alp = 0._dp
    sep%eisol = 0._dp
    sep%gss = 0._dp; sep%gsp = 0._dp; sep%gpp = 0._dp; sep%gp2 = 0._dp
    sep%gsd = 0._dp; sep%gpd = 0._dp; sep%gdd = 0._dp
    sep%hsp = 0._dp
    sep%dd = 0._dp; sep%qq = 0._dp; sep%am = 0._dp; sep%ad = 0._dp; sep%aq = 0._dp
    sep%fn1 = 0._dp; sep%fn2 = 0._dp; sep%fn3 = 0._dp
    sep%bfn1 = 0._dp; sep%bfn2 = 0._dp; sep%bfn3 = 0._dp

    SELECT CASE (z)
      CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "invalid atomic number")
      CASE (0)
        sep%defined = .TRUE.
      CASE (1)
        !Hydrogen
        sep%defined          = .TRUE.
        sep%uss  =  -13.0733210_dp
        sep%beta(0) =   -5.6265120_dp
        sep%sto_exponents(0)    =    0.9678070_dp
        sep%alp   =    3.3563860_dp
        sep%gss   =   14.7942080_dp
        sep%fn1(1) =  1.1287500_dp; sep%fn2(1) = 5.0962820_dp; sep%fn3(1) = 1.5374650_dp
        sep%fn1(2) = -1.0603290_dp; sep%fn2(2) = 6.0037880_dp; sep%fn3(2) = 1.5701890_dp
      CASE (2)
      CASE (3)
        !Lithium
        sep%defined          = .TRUE.
        sep%uss   =   -5.3000000_dp
        sep%upp   =   -3.4000000_dp
        sep%beta(0) =   -0.5500000_dp
        sep%beta(1)=   -1.5000000_dp
        sep%sto_exponents(0)    =    0.6500000_dp
        sep%sto_exponents(1)   =    0.7500000_dp
        sep%alp   =    1.2550000_dp
        sep%gss   =    9.0128510_dp
        sep%gsp  =    6.5761990_dp
        sep%gpp  =    6.0571820_dp
        sep%gp2  =    9.0052190_dp
        sep%hsp  =    0.5446790_dp
        sep%fn1(1) = -0.4500000_dp; sep%fn2(1) = 5.0000000_dp; sep%fn3(1) = 1.0000000_dp
        sep%fn1(2) =  0.8000000_dp; sep%fn2(2) = 6.5000000_dp; sep%fn3(2) = 1.0000000_dp
      CASE (4)
        !Beryllium
        sep%defined          = .TRUE.
        sep%uss   =  -17.2647520_dp
        sep%upp   =  -11.3042430_dp
        sep%beta(0) =   -3.9620530_dp
        sep%beta(1)=   -2.7806840_dp
        sep%sto_exponents(0)    =    0.8774390_dp
        sep%sto_exponents(1)   =    1.5087550_dp
        sep%alp   =    1.5935360_dp
        sep%gss   =    9.0128510_dp
        sep%gsp  =    6.5761990_dp
        sep%gpp  =    6.0571820_dp
        sep%gp2  =    9.0052190_dp
        sep%hsp  =    0.5446790_dp
        sep%fn1(1) =  1.6315720_dp; sep%fn2(1) = 2.6729620_dp; sep%fn3(1) = 1.7916860_dp
        sep%fn1(2) = -2.1109590_dp; sep%fn2(2) = 1.9685940_dp; sep%fn3(2) = 1.7558710_dp
      CASE (5)
        !Boron 
        sep%defined = .TRUE.
        sep%uss   =  -28.3539160_dp
        sep%upp   =  -23.9738780_dp
        sep%beta(0) =   -5.5215820_dp
        sep%beta(1)=   -5.9611770_dp
        sep%sto_exponents(0)    =    0.9885410_dp
        sep%sto_exponents(1)   =    1.9094410_dp
        sep%alp   =    2.4981070_dp
        sep%gss   =    8.8104530_dp
        sep%gsp  =   10.2992280_dp
        sep%gpp  =    4.9891680_dp
        sep%gp2  =    7.2899770_dp
        sep%hsp  =    0.8370320_dp
        sep%fn1(1) = 3.0046640_dp; sep%fn2(1) = 6.0288690_dp; sep%fn3(1) = 0.3472230_dp
        sep%fn1(2) = 0.0201950_dp; sep%fn2(2) = 6.0106580_dp; sep%fn3(2) = 1.2140670_dp
      CASE (6)
        !Carbon
        sep%defined = .TRUE.
        sep%uss   =  -47.2703200_dp
        sep%upp   =  -36.2669180_dp
        sep%beta(0) =  -11.9100150_dp
        sep%beta(1)=   -9.8027550_dp
        sep%sto_exponents(0)    =    1.5650850_dp
        sep%sto_exponents(1)   =    1.8423450_dp
        sep%alp   =    2.7078070_dp
        sep%gss   =   11.2007080_dp
        sep%gsp  =   10.2650270_dp
        sep%gpp  =   10.7962920_dp
        sep%gp2  =    9.0425660_dp
        sep%hsp  =    2.2909800_dp
        sep%fn1(1) = 0.0501070_dp; sep%fn2(1) = 6.0031650_dp; sep%fn3(1) = 1.6422140_dp
        sep%fn1(2) = 0.0507330_dp; sep%fn2(2) = 6.0029790_dp; sep%fn3(2) = 0.8924880_dp
      CASE (7)
      !Nitrogen
        sep%defined = .TRUE.
        sep%uss   =  -49.3356720_dp
        sep%upp   =  -47.5097360_dp
        sep%beta(0) =  -14.0625210_dp
        sep%beta(1)=  -20.0438480_dp
        sep%sto_exponents(0)    =    2.0280940_dp
        sep%sto_exponents(1)   =    2.3137280_dp
        sep%alp   =    2.8305450_dp
        sep%gss   =   11.9047870_dp
        sep%gsp  =    7.3485650_dp
        sep%gpp  =   11.7546720_dp
        sep%gp2  =   10.8072770_dp
        sep%hsp  =    1.1367130_dp
        sep%fn1(1) =  1.5016740_dp; sep%fn2(1) = 5.9011480_dp; sep%fn3(1) = 1.7107400_dp
        sep%fn1(2) = -1.5057720_dp; sep%fn2(2) = 6.0046580_dp; sep%fn3(2) = 1.7161490_dp
      CASE (8)
      !Oxygen
        sep%defined = .TRUE.
        sep%uss   =  -86.9930020_dp
        sep%upp   =  -71.8795800_dp
        sep%beta(0) =  -45.2026510_dp
        sep%beta(1)=  -24.7525150_dp
        sep%sto_exponents(0)    =    3.7965440_dp
        sep%sto_exponents(1)   =    2.3894020_dp
        sep%alp   =    3.2171020_dp
        sep%gss   =   15.7557600_dp
        sep%gsp  =   10.6211600_dp
        sep%gpp  =   13.6540160_dp
        sep%gp2  =   12.4060950_dp
        sep%hsp  =    0.5938830_dp
        sep%fn1(1) = -1.1311280_dp; sep%fn2(1) = 6.0024770_dp; sep%fn3(1) = 1.6073110_dp
        sep%fn1(2) =  1.1378910_dp; sep%fn2(2) = 5.9505120_dp; sep%fn3(2) = 1.5983950_dp
      CASE (9)
        !Fluorine
        sep%defined = .TRUE.
        sep%uss   = -110.4353030_dp
        sep%upp   = -105.6850470_dp
        sep%beta(0) =  -48.4059390_dp
        sep%beta(1)=  -27.7446600_dp
        sep%sto_exponents(0)    =    4.7085550_dp
        sep%sto_exponents(1)   =    2.4911780_dp
        sep%alp   =    3.3589210_dp
        sep%gss   =   10.4966670_dp
        sep%gsp  =   16.0736890_dp
        sep%gpp  =   14.8172560_dp
        sep%gp2  =   14.4183930_dp
        sep%hsp  =    0.7277630_dp
        sep%fn1(1) = -0.0121660_dp; sep%fn2(1) = 6.0235740_dp; sep%fn3(1) = 1.8568590_dp
        sep%fn1(2) = -0.0028520_dp; sep%fn2(2) = 6.0037170_dp; sep%fn3(2) = 2.6361580_dp
      CASE (10)
      CASE (11)
       !Sodium
        sep%defined = .TRUE.
        sep%defined = .TRUE.
        sep%uss   =  -4.7683450_dp
        sep%upp   =  -4.7703090_dp
        sep%beta(0) =   -0.5354170_dp
        sep%beta(1)=    0.0548760_dp
        sep%sto_exponents(0)    =    0.6767940_dp
        sep%sto_exponents(1)   =    0.9876530_dp
        sep%alp   =    1.6810000_dp
        sep%gss   =   5.1531420_dp
        sep%gsp  =    3.1778570_dp
        sep%gpp  =    3.7335330_dp
        sep%gp2  =   11.5386740_dp
        sep%hsp  =    0.5455510_dp
        sep%fn1(1) =  0.4175960_dp; sep%fn2(1) = 8.5969780_dp; sep%fn3(1) = 1.6731780_dp
        sep%fn1(2) =  0.4432600_dp; sep%fn2(2) = 3.3594960_dp; sep%fn3(2) = 2.2315110_dp
      CASE (12)
       !Magnesium
        sep%defined = .TRUE.
        sep%uss   =  -14.6236880_dp
        sep%upp   =  -14.1734600_dp
        sep%beta(0) =   -2.0716910_dp
        sep%beta(1)=   -0.5695810_dp
        sep%sto_exponents(0)    =    0.6985520_dp
        sep%sto_exponents(1)   =    1.4834530_dp
        sep%alp   =    1.3291470_dp
        sep%gss   =    6.6943000_dp
        sep%gsp  =    6.7939950_dp
        sep%gpp  =    6.9104460_dp
        sep%gp2  =    7.0908230_dp
        sep%hsp  =    0.5433000_dp
        sep%fn1(1) =  2.1170500_dp; sep%fn2(1) = 6.0094770_dp; sep%fn3(1) = 2.0844060_dp
        sep%fn1(2) = -2.5477670_dp; sep%fn2(2) = 4.3953700_dp; sep%fn3(2) = 2.0636740_dp
      CASE (13)
        !Aluminium
        sep%defined = .TRUE.
        sep%uss   =  -24.8454040_dp
        sep%upp   =  -22.2641590_dp
        sep%beta(0) =   -0.5943010_dp
        sep%beta(1)=   -0.9565500_dp
        sep%sto_exponents(0)    =    1.7028880_dp
        sep%sto_exponents(1)   =    1.0736290_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    1.5217030_dp
        sep%gss   =    5.7767370_dp
        sep%gsp  =   11.6598560_dp
        sep%gpp  =    6.3477900_dp
        sep%gp2  =    6.1210770_dp
        sep%hsp  =    4.0062450_dp
        sep%fn1(1) = -0.4730900_dp; sep%fn2(1) = 1.9158250_dp; sep%fn3(1) = 1.4517280_dp
        sep%fn1(2) = -0.1540510_dp; sep%fn2(2) = 6.0050860_dp; sep%fn3(2) = 2.5199970_dp
      CASE (14)
        !Silicon
        sep%defined = .TRUE.
        sep%uss   =  -26.7634830_dp
        sep%upp   =  -22.8136350_dp
        sep%beta(0) =   -2.8621450_dp
        sep%beta(1)=   -3.9331480_dp
        sep%sto_exponents(0)    =    1.6350750_dp
        sep%sto_exponents(1)   =    1.3130880_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    2.1358090_dp
        sep%gss   =    5.0471960_dp
        sep%gsp  =    5.9490570_dp
        sep%gpp  =    6.7593670_dp
        sep%gp2  =    5.1612970_dp
        sep%hsp  =    0.9198320_dp
        sep%fn1(1) = -0.3906000_dp; sep%fn2(1) = 6.0000540_dp; sep%fn3(1) = 0.6322620_dp
        sep%fn1(2) =  0.0572590_dp; sep%fn2(2) = 6.0071830_dp; sep%fn3(2) = 2.0199870_dp
      CASE (15)
        !Phosphorus
        sep%defined = .TRUE.
        sep%uss   =  -40.4130960_dp
        sep%upp   =  -29.5930520_dp
        sep%beta(0) =  -12.6158790_dp
        sep%beta(1)=   -4.1600400_dp
        sep%sto_exponents(0)    =    2.0175630_dp
        sep%sto_exponents(1)   =    1.5047320_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    1.9405340_dp
        sep%gss   =    7.8016150_dp
        sep%gsp  =    5.1869490_dp
        sep%gpp  =    6.6184780_dp
        sep%gp2  =    6.0620020_dp
        sep%hsp  =    1.5428090_dp
        sep%fn1(1) = -0.6114210_dp; sep%fn2(1) = 1.9972720_dp; sep%fn3(1) = 0.7946240_dp
        sep%fn1(2) = -0.0939350_dp; sep%fn2(2) = 1.9983600_dp; sep%fn3(2) = 1.9106770_dp
      CASE (16)
      !Sulphur
        sep%defined = .TRUE.
        sep%uss   =  -49.8953710_dp
        sep%upp   =  -44.3925830_dp
        sep%beta(0) =   -8.8274650_dp
        sep%beta(1)=   -8.0914150_dp
        sep%sto_exponents(0)    =    1.8911850_dp
        sep%sto_exponents(1)   =    1.6589720_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    2.2697060_dp
        sep%gss   =    8.9646670_dp
        sep%gsp  =    6.7859360_dp
        sep%gpp  =    9.9681640_dp
        sep%gp2  =    7.9702470_dp
        sep%hsp  =    4.0418360_dp
        sep%fn1(1) = -0.3991910_dp; sep%fn2(1) = 6.0006690_dp; sep%fn3(1) = 0.9621230_dp
        sep%fn1(2) = -0.0548990_dp; sep%fn2(2) = 6.0018450_dp; sep%fn3(2) = 1.5799440_dp
      CASE (17)
      !Chlorine
        sep%defined = .TRUE.
        sep%uss   = -100.6267470_dp
        sep%upp   =  -53.6143960_dp
        sep%beta(0) =  -27.5285600_dp
        sep%beta(1)=  -11.5939220_dp
        sep%sto_exponents(0)    =    2.2462100_dp
        sep%sto_exponents(1)   =    2.1510100_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    2.5172960_dp
        sep%gss   =   16.0136010_dp
        sep%gsp  =    8.0481150_dp
        sep%gpp  =    7.5222150_dp
        sep%gp2  =    7.5041540_dp
        sep%hsp  =    3.4811530_dp
        sep%fn1(1) = -0.1715910_dp; sep%fn2(1) = 6.0008020_dp; sep%fn3(1) = 1.0875020_dp
        sep%fn1(2) = -0.0134580_dp; sep%fn2(2) = 1.9666180_dp; sep%fn3(2) = 2.2928910_dp
      CASE (18)
      CASE (19)
        !Potassium
        sep%defined = .TRUE.
        sep%uss   = -4.1079950_dp
        sep%upp   =  -1.9126360_dp
        sep%beta(0) =  -0.7235020_dp
        sep%beta(1)=    0.0109380_dp
        sep%sto_exponents(0)   =    0.9158730_dp
        sep%sto_exponents(1)   =    1.3430990_dp
        sep%alp   =    4.4452360_dp
        sep%gss   =    4.0502660_dp
        sep%gsp  =    2.9879670_dp
        sep%gpp  =    2.9983050_dp
        sep%gp2  =    4.9936140_dp
        sep%hsp  =    0.5653070_dp
        sep%fn1(1) = 6.0409170_dp; sep%fn2(1) = 10.0465540_dp; sep%fn3(1) = 1.6005120_dp
        sep%fn1(2) = 0.7242240_dp; sep%fn2(2) =  3.0082750_dp; sep%fn3(2) = 2.4643520_dp
      CASE (20)
        !Calcium   
        sep%defined = .TRUE.
        sep%uss   =-11.6466770_dp
        sep%upp   = -11.7003040_dp
        sep%beta(0) =  -1.3663020_dp
        sep%beta(1)=   -0.6800840_dp
        sep%sto_exponents(0)   =    0.8224530_dp
        sep%sto_exponents(1)   =    1.2791170_dp
        sep%alp   =    3.2251430_dp
        sep%gss   =    5.5661490_dp
        sep%gsp  =   15.8559960_dp
        sep%gpp  =    6.7057840_dp
        sep%gp2  =    4.8642700_dp
        sep%hsp  =    1.0053200_dp
        sep%fn1(1) = 0.7940250_dp; sep%fn2(1) =  3.8548560_dp; sep%fn3(1) = 1.5284410_dp
        sep%fn1(2) = 0.3400220_dp; sep%fn2(2) =  5.7455560_dp; sep%fn3(2) = 2.3571110_dp
      CASE (21:23)
!                    DATA FOR ELEMENT 22 Titanium (Russell-Saunders only)
!  data alppm3(22)/3.0683070e0_dp/
!  data betapp(22)/-0.1000000e0_dp/
!  data betasp(22)/-0.1000000e0_dp/
!  data gp2pm3(22)/3.5000000e0_dp/
!  data gpppm3(22)/5.0000000e0_dp/
!  data gsppm3(22)/4.1500000e0_dp/
!  data gsspm3(22)/6.0000000e0_dp/
!  data hsppm3(22)/1.0000000e0_dp/
!  data uddpm3(22)/-30.0000000e0_dp/
!  data upppm3(22)/10.0000000e0_dp/
!  data usspm3(22)/10.0000000e0_dp/
!  data zdpm3(22)/2.8845490e0_dp/
!  data zppm3(22)/1.5000000e0_dp/
!  data zspm3(22)/1.5000000e0_dp/
      CASE (24)
        !Chromuim
        sep%defined = .TRUE.
        sep%uss   =-17.5170270_dp
        sep%upp   = -12.5337290_dp
        sep%udd   = -44.1249280_dp
        sep%beta(0) =  -0.1000000_dp
        sep%beta(1)=   -0.1000000_dp
        sep%beta(2)=   -8.7766360_dp
        sep%sto_exponents(0)   =    1.5000000_dp
        sep%sto_exponents(1)   =    1.5000000_dp
        sep%sto_exponents(1)   =    2.8845490_dp
        sep%alp   =    3.0683070_dp
        sep%gss   =    6.0000000_dp
        sep%gsp  =    4.1500000_dp
        sep%gpp  =    5.0000000_dp
        sep%gp2  =    3.5000000_dp
        sep%hsp  =    1.0000000_dp
      CASE (25:29)
      CASE (30)
        !Zinc
        sep%defined = .TRUE.
        sep%uss   =  -18.5321980_dp
        sep%upp   =  -11.0474090_dp
        sep%beta(0) =   -0.7155780_dp
        sep%beta(1)=   -6.3518640_dp
        sep%sto_exponents(0)    =    1.8199890_dp
        sep%sto_exponents(1)   =    1.5069220_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    1.3501260_dp
        sep%gss   =    9.6771960_dp
        sep%gsp  =    7.7362040_dp
        sep%gpp  =    4.9801740_dp
        sep%gp2  =    4.6696560_dp
        sep%hsp  =    0.6004130_dp
        sep%fn1(1) = -0.1112340_dp; sep%fn2(1) = 6.0014780_dp; sep%fn3(1) = 1.5160320_dp
        sep%fn1(2) = -0.1323700_dp; sep%fn2(2) = 1.9958390_dp; sep%fn3(2) = 2.5196420_dp
      CASE (31)
        !Gallium
        sep%defined = .TRUE.
        sep%uss   =  -29.8555930_dp
        sep%upp   =  -21.8753710_dp
        sep%beta(0) =   -4.9456180_dp
        sep%beta(1)=   -0.4070530_dp
        sep%sto_exponents(0)    =    1.8470400_dp
        sep%sto_exponents(1)   =    0.8394110_dp
        sep%alp   =    1.6051150_dp
        sep%gss   =    8.4585540_dp
        sep%gsp  =    8.9256190_dp
        sep%gpp  =    5.0868550_dp
        sep%gp2  =    4.9830450_dp
        sep%hsp  =    2.0512600_dp
        sep%fn1(1) = -0.5601790_dp; sep%fn2(1) = 5.6232730_dp; sep%fn3(1) = 1.5317800_dp
        sep%fn1(2) = -0.2727310_dp; sep%fn2(2) = 1.9918430_dp; sep%fn3(2) = 2.1838640_dp
      CASE (32)
        !Germanium
        sep%defined = .TRUE.
        sep%uss   =  -35.4671955_dp
        sep%upp   =  -31.5863583_dp
        sep%beta(0) =   -5.3250024_dp
        sep%beta(1)=   -2.2501567_dp
        sep%sto_exponents(0)    =    2.2373526_dp
        sep%sto_exponents(1)   =    1.5924319_dp
        sep%alp   =    1.9723370_dp
        sep%gss   =    5.3769635_dp
        sep%gsp  =   10.2095293_dp
        sep%gpp  =    7.6718647_dp
        sep%gp2  =    6.9242663_dp
        sep%hsp  =    1.3370204_dp
        sep%fn1(1) =  0.9631726_dp; sep%fn2(1) = 6.0120134_dp; sep%fn3(1) = 2.1633655_dp
        sep%fn1(2) = -0.9593891_dp; sep%fn2(2) = 5.7491802_dp; sep%fn3(2) = 2.1693724_dp
      CASE (33)
        !Arsenic
        sep%defined = .TRUE.
        sep%uss   =  -38.5074240_dp
        sep%upp   =  -35.1524150_dp
        sep%beta(0) =   -8.2321650_dp
        sep%beta(1)=   -5.0173860_dp
        sep%sto_exponents(0)    =    2.6361770_dp
        sep%sto_exponents(1)   =    1.7038890_dp
        sep%alp   =    1.7944770_dp
        sep%gss   =    8.7890010_dp
        sep%gsp  =    5.3979830_dp
        sep%gpp  =    8.2872500_dp
        sep%gp2  =    8.2103460_dp
        sep%hsp  =    1.9510340_dp
        sep%fn1(1) = -0.4600950_dp; sep%fn2(1) = 1.9831150_dp; sep%fn3(1) = 1.0867930_dp
        sep%fn1(2) = -0.0889960_dp; sep%fn2(2) = 1.9929440_dp; sep%fn3(2) = 2.1400580_dp
      CASE (34)
        !Selenium
        sep%defined = .TRUE.
        sep%uss   =  -55.3781350_dp
        sep%upp   =  -49.8230760_dp
        sep%beta(0) =   -6.1578220_dp
        sep%beta(1)=   -5.4930390_dp
        sep%sto_exponents(0)    =    2.8280510_dp
        sep%sto_exponents(1)   =    1.7325360_dp
        sep%alp   =    3.0439570_dp
        sep%gss   =    7.4325910_dp
        sep%gsp  =   10.0604610_dp
        sep%gpp  =    9.5683260_dp
        sep%gp2  =    7.7242890_dp
        sep%hsp  =    4.0165580_dp
        sep%fn1(1) = 0.0478730_dp; sep%fn2(1) = 6.0074000_dp; sep%fn3(1) = 2.0817170_dp
        sep%fn1(2) = 0.1147200_dp; sep%fn2(2) = 6.0086720_dp; sep%fn3(2) = 1.5164230_dp
      CASE (35)
        !Bromine
        sep%defined = .TRUE.
        sep%uss   = -116.6193110_dp
        sep%upp   =  -74.2271290_dp
        sep%beta(0) =  -31.1713420_dp
        sep%beta(1)=   -6.8140130_dp
        sep%sto_exponents(0)    =    5.3484570_dp
        sep%sto_exponents(1)   =    2.1275900_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    2.5118420_dp
        sep%gss   =   15.9434250_dp
        sep%gsp  =   16.0616800_dp
        sep%gpp  =    8.2827630_dp
        sep%gp2  =    7.8168490_dp
        sep%hsp  =    0.5788690_dp
        sep%fn1(1) =  0.9604580_dp; sep%fn2(1) = 5.9765080_dp; sep%fn3(1) = 2.3216540_dp
        sep%fn1(2) = -0.9549160_dp; sep%fn2(2) = 5.9447030_dp; sep%fn3(2) = 2.3281420_dp
      CASE (36)
      CASE (37)
!     DATA USSPM3( 37)/      -4.4841970e0_dp/
!     DATA UPPPM3( 37)/      -3.0722520e0_dp/
!     DATA BETASP( 37)/      -0.7442240e0_dp/
!     DATA BETAPP( 37)/      -0.5677060e0_dp/
!     DATA ZSPM3 ( 37)/       0.6801340e0_dp/
!     DATA ZPPM3 ( 37)/       0.9658140e0_dp/
!     DATA ALPPM3( 37)/       5.7668440e0_dp/
!     DATA EISOLP( 37)/      -4.4841970e0_dp/
!     DATA GSSPM3( 37)/       4.9848710e0_dp/
!     DATA GSPPM3( 37)/       4.9983590e0_dp/
!     DATA GPPPM3( 37)/       6.9857170e0_dp/
!     DATA GP2PM3( 37)/       5.5079160e0_dp/
!     DATA HSPPM3( 37)/       0.6211430e0_dp/
!     DATA DDPM3 ( 37)/       3.2610134e0_dp/
!     DATA QQPM3 ( 37)/       2.6599806e0_dp/
!     DATA AMPM3 ( 37)/       0.1832000e0_dp/
!     DATA ADPM3 ( 37)/       0.1544789e0_dp/
!     DATA AQPM3 ( 37)/       0.3073402e0_dp/
!     DATA GuesP1( 37,1)/       4.7724480e0_dp/
!     DATA GuesP2( 37,1)/       5.8039230e0_dp/
!     DATA GuesP3( 37,1)/       1.4442230e0_dp/
!     DATA GuesP1( 37,2)/       0.5062010e0_dp/
!     DATA GuesP2( 37,2)/       2.8303820e0_dp/
!     DATA GuesP3( 37,2)/       2.6433260e0_dp/
!     Data (GuesP1(37,J),J=3,10)/8*0.0e0_dp/
!     Data (GuesP2(37,J),J=3,10)/8*0.0e0_dp/
!     Data (GuesP3(37,J),J=3,10)/8*0.0e0_dp/
      CASE (38)
!     DATA USSPM3( 38)/     -11.4870450e0_dp/
!     DATA UPPPM3( 38)/      -7.7152160e0_dp/
!     DATA BETASP( 38)/      -0.7533420e0_dp/
!     DATA BETAPP( 38)/      -0.5760720e0_dp/
!     DATA ZSPM3 ( 38)/       1.5592340e0_dp/
!     DATA ZPPM3 ( 38)/       1.2129820e0_dp/
!     DATA ALPPM3( 38)/       4.5905860e0_dp/
!     DATA EISOLP( 38)/     -16.9990050e0_dp/
!     DATA GSSPM3( 38)/       5.9750850e0_dp/
!     DATA GSPPM3( 38)/       4.8869600e0_dp/
!     DATA GPPPM3( 38)/       7.3407630e0_dp/
!     DATA GP2PM3( 38)/      15.8140850e0_dp/
!     DATA HSPPM3( 38)/       4.0512180e0_dp/
!     DATA DDPM3 ( 38)/       2.1011077e0_dp/
!     DATA QQPM3 ( 38)/       2.1179593e0_dp/
!     DATA AMPM3 ( 38)/       0.2195915e0_dp/
!     DATA ADPM3 ( 38)/       0.5137005e0_dp/
!     DATA AQPM3 ( 38)/       0.1752812e0_dp/
!     DATA GuesP1( 38,1)/       0.4398870e0_dp/
!     DATA GuesP2( 38,1)/       5.7510950e0_dp/
!     DATA GuesP3( 38,1)/       1.6731940e0_dp/
!     DATA GuesP1( 38,2)/       0.3957700e0_dp/
!     DATA GuesP2( 38,2)/       4.4993190e0_dp/
!     DATA GuesP3( 38,2)/       2.3281550e0_dp/
!     Data (GuesP1(38,J),J=3,10)/8*0.0e0_dp/
!     Data (GuesP2(38,J),J=3,10)/8*0.0e0_dp/
!     Data (GuesP3(38,J),J=3,10)/8*0.0e0_dp/
      CASE (39:47)
      CASE (48)
        !Cadmium
        sep%defined = .TRUE.
        sep%uss   =  -15.8285840_dp
        sep%upp   =    8.7497950_dp
        sep%beta(0) =   -8.5819440_dp
        sep%beta(1)=   -0.6010340_dp
        sep%sto_exponents(0)    =    1.6793510_dp
        sep%sto_exponents(1)   =    2.0664120_dp
        sep%alp   =    1.5253820_dp
        sep%gss   =    9.2069600_dp
        sep%gsp  =    8.2315390_dp
        sep%gpp  =    4.9481040_dp
        sep%gp2  =    4.6696560_dp
        sep%hsp  =    1.6562340_dp
      CASE (49)
        !Indium
        sep%defined = .TRUE.
        sep%uss   =  -26.1762050_dp
        sep%upp   =  -20.0058220_dp
        sep%beta(0) =   -2.9933190_dp
        sep%beta(1)=   -1.8289080_dp
        sep%sto_exponents(0)    =    2.0161160_dp
        sep%sto_exponents(1)   =    1.4453500_dp
        sep%alp   =    1.4183850_dp
        sep%gss   =    6.5549000_dp
        sep%gsp  =    8.2298730_dp
        sep%gpp  =    6.2992690_dp
        sep%gp2  =    4.9842110_dp
        sep%hsp  =    2.6314610_dp
        sep%fn1(1) = -0.3431380_dp; sep%fn2(1) = 1.9940340_dp; sep%fn3(1) = 1.6255160_dp
        sep%fn1(2) = -0.1095320_dp; sep%fn2(2) = 5.6832170_dp; sep%fn3(2) = 2.8670090_dp
      CASE (50)
        !Tin
        sep%defined = .TRUE.
        sep%uss   =  -34.5501920_dp
        sep%upp   =  -25.8944190_dp
        sep%beta(0) =   -2.7858020_dp
        sep%beta(1)=   -2.0059990_dp
        sep%sto_exponents(0)    =    2.3733280_dp
        sep%sto_exponents(1)   =    1.6382330_dp
        sep%alp   =    1.6996500_dp
        sep%gss   =   10.1900330_dp
        sep%gsp  =    7.2353270_dp
        sep%gpp  =    5.6738100_dp
        sep%gp2  =    5.1822140_dp
        sep%hsp  =    1.0331570_dp
        sep%fn1(1) = -0.1503530_dp; sep%fn2(1) = 6.0056940_dp; sep%fn3(1) = 1.7046420_dp
        sep%fn1(2) = -0.0444170_dp; sep%fn2(2) = 2.2573810_dp; sep%fn3(2) = 2.4698690_dp
      CASE (51)
        !Antimony
        sep%defined = .TRUE.
        sep%uss   =  -56.4321960_dp
        sep%upp   =  -29.4349540_dp
        sep%beta(0) =  -14.7942170_dp
        sep%beta(1)=   -2.8179480_dp
        sep%sto_exponents(0)    =    2.3430390_dp
        sep%sto_exponents(1)   =    1.8999920_dp
        sep%alp   =    2.0343010_dp
        sep%gss   =    9.2382770_dp
        sep%gsp  =    5.2776800_dp
        sep%gpp  =    6.3500000_dp
        sep%gp2  =    6.2500000_dp
        sep%hsp  =    2.4244640_dp
        sep%fn1(1) =  3.0020280_dp; sep%fn2(1) = 6.0053420_dp; sep%fn3(1) = 0.8530600_dp
        sep%fn1(2) = -0.0188920_dp; sep%fn2(2) = 6.0114780_dp; sep%fn3(2) = 2.7933110_dp
      CASE (52)
        !Tellurium
        sep%defined = .TRUE.
        sep%uss   =  -44.9380360_dp
        sep%upp   =  -46.3140990_dp
        sep%beta(0) =   -2.6651460_dp
        sep%beta(1)=   -3.8954300_dp
        sep%sto_exponents(0)    =    4.1654920_dp
        sep%sto_exponents(1)   =    1.6475550_dp
        sep%alp   =    2.4850190_dp
        sep%gss   =   10.2550730_dp
        sep%gsp  =    8.1691450_dp
        sep%gpp  =    7.7775920_dp
        sep%gp2  =    7.7551210_dp
        sep%hsp  =    3.7724620_dp
        sep%fn1(1) =  0.0333910_dp; sep%fn2(1) = 5.9563790_dp; sep%fn3(1) = 2.2775750_dp 
        sep%fn1(2) = -1.9218670_dp; sep%fn2(2) = 4.9732190_dp; sep%fn3(2) = 0.5242430_dp
      CASE (53)
        !Iodine
        sep%defined = .TRUE.
        sep%uss   =  -96.4540370_dp
        sep%upp   =  -61.0915820_dp
        sep%beta(0) =  -14.4942340_dp
        sep%beta(1)=   -5.8947030_dp
        sep%sto_exponents(0)    =    7.0010130_dp
        sep%sto_exponents(1)   =    2.4543540_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    1.9901850_dp
        sep%gss   =   13.6319430_dp
        sep%gsp  =   14.9904060_dp
        sep%gpp  =    7.2883300_dp
        sep%gp2  =    5.9664070_dp
        sep%hsp  =    2.6300350_dp
        sep%fn1(1) = -0.1314810_dp; sep%fn2(1) = 5.2064170_dp; sep%fn3(1) = 1.7488240_dp
        sep%fn1(2) = -0.0368970_dp; sep%fn2(2) = 6.0101170_dp; sep%fn3(2) = 2.7103730_dp
      CASE (54)
      CASE (55)
!     DATA USSPM3( 55)/      -4.3207390e0_dp/
!     DATA UPPPM3( 55)/      -3.3475800e0_dp/
!     DATA BETASP( 55)/      -0.9775130e0_dp/
!     DATA BETAPP( 55)/      -0.5484990e0_dp/
!     DATA ZSPM3 ( 55)/       0.9188020e0_dp/
!     DATA ZPPM3 ( 55)/       1.5522530e0_dp/
!     DATA ALPPM3( 55)/       2.5439870e0_dp/
!     DATA EISOLP( 55)/      -4.3207390e0_dp/
!     DATA GSSPM3( 55)/       4.9615320e0_dp/
!     DATA GSPPM3( 55)/       4.9943210e0_dp/
!     DATA GPPPM3( 55)/       5.0115170e0_dp/
!     DATA GP2PM3( 55)/      16.0009950e0_dp/
!     DATA HSPPM3( 55)/       0.5830050e0_dp/
!     DATA DDPM3 ( 55)/       1.9526201e0_dp/
!     DATA QQPM3 ( 55)/       1.9433821e0_dp/
!     DATA AMPM3 ( 55)/       0.1823422e0_dp/
!     DATA ADPM3 ( 55)/       0.2009234e0_dp/
!     DATA AQPM3 ( 55)/       0.1864057e0_dp/
!     DATA GuesP1( 55,1)/       3.3960620e0_dp/
!     DATA GuesP2( 55,1)/       6.0007540e0_dp/
!     DATA GuesP3( 55,1)/       1.6762910e0_dp/
!     DATA GuesP1( 55,2)/       0.6677680e0_dp/
!     DATA GuesP2( 55,2)/       3.0084010e0_dp/
!     DATA GuesP3( 55,2)/       2.7188720e0_dp/
!     Data (GuesP1(55,J),J=3,10)/8*0.0e0_dp/
!     Data (GuesP2(55,J),J=3,10)/8*0.0e0_dp/
!     Data (GuesP3(55,J),J=3,10)/8*0.0e0_dp/
      CASE (56)
!     DATA USSPM3( 56)/     -13.5286010e0_dp/
!     DATA UPPPM3( 56)/     -10.1346170e0_dp/
!     DATA BETASP( 56)/      -3.9502870e0_dp/
!     DATA BETAPP( 56)/      -0.5749730e0_dp/
!     DATA ZSPM3 ( 56)/       1.7594520e0_dp/
!     DATA ZPPM3 ( 56)/       1.4916880e0_dp/
!     DATA ALPPM3( 56)/       1.6584300e0_dp/
!     DATA EISOLP( 56)/     -16.9427990e0_dp/
!     DATA GSSPM3( 56)/      10.1144030e0_dp/
!     DATA GSPPM3( 56)/       4.8886780e0_dp/
!     DATA GPPPM3( 56)/       8.8650930e0_dp/
!     DATA GP2PM3( 56)/       8.7503680e0_dp/
!     DATA HSPPM3( 56)/       1.4151370e0_dp/
!     DATA DDPM3 ( 56)/       2.2086835e0_dp/
!     DATA QQPM3 ( 56)/       2.0222866e0_dp/
!     DATA AMPM3 ( 56)/       0.3717164e0_dp/
!     DATA ADPM3 ( 56)/       0.2800750e0_dp/
!     DATA AQPM3 ( 56)/       0.1811672e0_dp/
!     DATA GuesP1( 56,1)/       0.3501710e0_dp/
!     DATA GuesP2( 56,1)/       5.8729050e0_dp/
!     DATA GuesP3( 56,1)/       1.9750730e0_dp/
!     DATA GuesP1( 56,2)/       0.3339000e0_dp/
!     DATA GuesP2( 56,2)/       4.0328190e0_dp/
!     DATA GuesP3( 56,2)/       2.6276060e0_dp/
!     Data (GuesP1(56,J),J=3,10)/8*0.0e0_dp/
!     Data (GuesP2(56,J),J=3,10)/8*0.0e0_dp/
!     Data (GuesP3(56,J),J=3,10)/8*0.0e0_dp/
      CASE (57:79)
      CASE (80)
        !Mercury
        sep%defined = .TRUE.
        sep%uss   =  -17.7622290_dp
        sep%upp   =  -18.3307510_dp
        sep%beta(0) =   -3.1013650_dp
        sep%beta(1)=   -3.4640310_dp
        sep%sto_exponents(0)    =    1.4768850_dp
        sep%sto_exponents(1)   =    2.4799510_dp
        sep%alp   =    1.5293770_dp
        sep%gss   =    6.6247200_dp
        sep%gsp  =   10.6392970_dp
        sep%gpp  =   14.7092830_dp
        sep%gp2  =   16.0007400_dp
        sep%hsp  =    2.0363110_dp
        sep%fn1(1) =  1.0827200_dp; sep%fn2(1) = 6.4965980_dp; sep%fn3(1) = 1.1951460_dp
        sep%fn1(2) = -0.0965530_dp; sep%fn2(2) = 3.9262810_dp; sep%fn3(2) = 2.6271600_dp
      CASE (81)
        !Thallium
        sep%defined = .TRUE.
        sep%uss   =  -30.0531700_dp
        sep%upp   =  -26.9206370_dp
        sep%beta(0) =   -1.0844950_dp
        sep%beta(1)=   -7.9467990_dp
        sep%sto_exponents(0)    =    6.8679210_dp
        sep%sto_exponents(1)   =    1.9694450_dp
        sep%alp   =    1.3409510_dp
        sep%gss   =   10.4604120_dp
        sep%gsp  =   11.2238830_dp
        sep%gpp  =    4.9927850_dp
        sep%gp2  =    8.9627270_dp
        sep%hsp  =    2.5304060_dp
        sep%fn1(1) = -1.3613990_dp; sep%fn2(1) = 3.5572260_dp; sep%fn3(1) = 1.0928020_dp
        sep%fn1(2) = -0.0454010_dp; sep%fn2(2) = 2.3069950_dp; sep%fn3(2) = 2.9650290_dp
      CASE (82)
        !Lead
        sep%defined = .TRUE.
        sep%uss   =  -30.3227560_dp
        sep%upp   =  -24.4258340_dp
        sep%beta(0) =   -6.1260240_dp
        sep%beta(1)=   -1.3954300_dp
        sep%sto_exponents(0)    =    3.1412890_dp
        sep%sto_exponents(1)   =    1.8924180_dp
        sep%alp   =    1.6200450_dp
        sep%gss   =    7.0119920_dp
        sep%gsp  =    6.7937820_dp
        sep%gpp  =    5.1837800_dp
        sep%gp2  =    5.0456510_dp
        sep%hsp  =    1.5663020_dp
        sep%fn1(1) = -0.1225760_dp; sep%fn2(1) = 6.0030620_dp; sep%fn3(1) = 1.9015970_dp
        sep%fn1(2) = -0.0566480_dp; sep%fn2(2) = 4.7437050_dp; sep%fn3(2) = 2.8618790_dp
      CASE (83)
        !Bismuth
        sep%defined = .TRUE.
        sep%uss   =  -33.4959380_dp
        sep%upp   =  -35.5210260_dp
        sep%beta(0) =   -5.6072830_dp
        sep%beta(1) =   -5.8001520_dp
        sep%sto_exponents(0)    =    4.9164510_dp
        sep%sto_exponents(1)   =    1.9349350_dp
        sep%alp   =    1.8574310_dp
        sep%gss   =    4.9894800_dp
        sep%gsp  =    6.1033080_dp
        sep%gpp  =    8.6960070_dp
        sep%gp2  =    8.3354470_dp
        sep%hsp  =    0.5991220_dp
        sep%fn1(1) = 2.5816930_dp; sep%fn2(1) = 5.0940220_dp; sep%fn3(1) = 0.4997870_dp
        sep%fn1(2) = 0.0603200_dp; sep%fn2(2) = 6.0015380_dp; sep%fn3(2) = 2.4279700_dp
      CASE (84:89)
      CASE (90)
!     DATA USSPM3( 90)/     -40.5682920e0_dp/
!     DATA UPPPM3( 90)/     -28.0891870e0_dp/
!     DATA BETASP( 90)/      -4.2562180e0_dp/
!     DATA BETAPP( 90)/      -4.2562180e0_dp/
!     DATA ZSPM3 ( 90)/       1.4353060e0_dp/
!     DATA ZPPM3 ( 90)/       1.4353060e0_dp/
!     DATA ZDPM3 ( 90)/       1.0000000e0_dp/
!     DATA ALPPM3( 90)/       2.1961078e0_dp/
!     DATA EISOLP( 90)/       0.0000000e0_dp/
!     DATA GSSPM3( 90)/       9.8200000e0_dp/
!     DATA GSPPM3( 90)/       8.3600000e0_dp/
!     DATA GPPPM3( 90)/       7.3100000e0_dp/
!     DATA GP2PM3( 90)/       6.5400000e0_dp/
!     DATA HSPPM3( 90)/       1.3200000e0_dp/
!     DATA DDPM3 ( 90)/       0.2011245e0_dp/
!     DATA QQPM3 ( 90)/       0.2203208e0_dp/
!     DATA AMPM3 ( 90)/       0.3608967e0_dp/
!     DATA ADPM3 ( 90)/       1.1124287e0_dp/
!     DATA AQPM3 ( 90)/       1.3153602e0_dp/
!     Data (GuesP1(90,J),J=1,10)/10*0.0e0_dp/
!     Data (GuesP2(90,J),J=1,10)/10*0.0e0_dp/
!     Data (GuesP3(90,J),J=1,10)/10*0.0e0_dp/
      CASE (91)
!     DATA USSPM3( 91)/     -75.2391520e0_dp/
!     DATA UPPPM3( 91)/     -57.8320130e0_dp/
!     DATA BETASP( 91)/     -11.1422310e0_dp/
!     DATA BETAPP( 91)/     -11.1422310e0_dp/
!     DATA ZSPM3 ( 91)/       2.6135910e0_dp/
!     DATA ZPPM3 ( 91)/       2.0343930e0_dp/
!     DATA ZDPM3 ( 91)/       1.0000000e0_dp/
!     DATA ALPPM3( 91)/       2.4916445e0_dp/
!     DATA EISOLP( 91)/       0.0000000e0_dp/
!     DATA HSPPM3( 91)/       0.0000001e0_dp/
!     DATA DDPM3 ( 91)/       0.1232470e0_dp/
!     DATA QQPM3 ( 91)/       0.1554408e0_dp/
!     DATA AMPM3 ( 91)/       0.0000000e0_dp/
!     DATA ADPM3 ( 91)/       0.0062312e0_dp/
!     DATA AQPM3 ( 91)/       1.2399130e0_dp/
!     Data (GuesP1(91,J),J=1,10)/10*0.0e0_dp/
!     Data (GuesP2(91,J),J=1,10)/10*0.0e0_dp/
!     Data (GuesP3(91,J),J=1,10)/10*0.0e0_dp/
      CASE (92:106)
    END SELECT
    CALL calpar(z,sep)

    sep%beta = sep%beta/evolt
    sep%uss = sep%uss/evolt; sep%upp = sep%upp/evolt; sep%udd = sep%udd/evolt
    sep%alp = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss = sep%gss/evolt; sep%gsp = sep%gsp/evolt; sep%gpp = sep%gpp/evolt; sep%gp2 = sep%gp2/evolt 
    sep%gsd = sep%gsd/evolt; sep%gpd = sep%gpd/evolt; sep%gdd = sep%gdd/evolt
    sep%hsp = sep%hsp/evolt
    sep%fn1 = sep%fn1*bohr/evolt; sep%fn2 = sep%fn2/bohr/bohr; sep%fn3 = sep%fn3*bohr

  END SUBROUTINE pm3_default_parameter

! *****************************************************************************

  SUBROUTINE pdg_default_parameter ( sep, z )
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z

    CHARACTER(len=*), PARAMETER :: routineN = 'pdg_default_parameter', &
      routineP = moduleN//':'//routineN

! -----------------------------------------------------------------------------

    sep%typ = do_method_pdg
    sep%name = ptable(z)%name
    sep%eheat = ptable(z)%heat_of_formation / kcalmol
    sep%z = z
    CALL valence_electrons ( z, sep%zeff, sep%natorb )
    sep%ngauss = 6
    sep%defined = .FALSE.
    sep%sto_exponents(0:3) = 0._dp
    sep%beta(0:3) = 0._dp
    sep%uss = 0._dp; sep%upp = 0._dp; sep%udd = 0._dp
    sep%alp = 0._dp
    sep%eisol = 0._dp
    sep%gss = 0._dp; sep%gsp = 0._dp; sep%gpp = 0._dp; sep%gp2 = 0._dp
    sep%gsd = 0._dp; sep%gpd = 0._dp; sep%gdd = 0._dp
    sep%hsp = 0._dp
    sep%dd = 0._dp; sep%qq = 0._dp; sep%am = 0._dp; sep%ad = 0._dp; sep%aq = 0._dp
    sep%fn1 = 0._dp; sep%fn2 = 0._dp; sep%fn3 = 0._dp
    sep%bfn1 = 0._dp; sep%bfn2 = 0._dp; sep%bfn3 = 0._dp

    SELECT CASE (z)
      CASE DEFAULT
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "invalid atomic number")
      CASE (0)
        sep%defined = .TRUE.
      CASE (1)
        !Hydrogen
        sep%defined          = .TRUE.
        sep%uss  =  -12.893272_dp
        sep%beta(0) =   -6.152654
        sep%sto_exponents(0)    =    0.972786_dp
        sep%alp   =    3.381686_dp
        sep%eisol =  -13.120566_dp
        sep%gss   =   14.7942080_dp
        sep%am    =    0.5_dp/0.919616_dp
        sep%ad    =    0.5_dp/0.919616_dp
        sep%aq   =    0.5_dp/0.919616_dp
        sep%fn1(1) =  1.122244_dp 
        sep%fn2(1) = 4.707790_dp
        sep%fn3(1) = 1.547099_dp
        sep%fn1(2) = -1.069737_dp
        sep%fn2(2) = 5.857995_dp
        sep%fn3(2) = 1.567893_dp
        sep%pre(1) = 0.057193_dp
        sep%pre(2) = -0.034823_dp
        sep%d(1) = 0.663395_dp
        sep%d(2) = 1.081901_dp
      CASE (2:5)
      CASE (6)
        !Carbon
        sep%defined = .TRUE.
        sep%uss   =  -48.241241_dp
        sep%upp   =  -36.461256_dp
        sep%beta(0) =  -11.952818_dp
        sep%beta(1)=   -9.922411_dp
        sep%sto_exponents(0)    =    1.567864_dp
        sep%sto_exponents(1)   =    1.846659_dp
        sep%alp   =    2.725772_dp
        sep%eisol = -113.428242_dp
        sep%gss   =   11.2007080_dp
        sep%gsp  =   10.2650270_dp
        sep%gpp  =   10.7962920_dp
        sep%gp2  =    9.0425660_dp
        sep%hsp  =    2.2909800_dp
        sep%dd   =    0.831413_dp
        sep%qq    =   0.663222_dp
        sep%am    =   0.5_dp/1.214657_dp
        sep%ad    =   0.5_dp/0.848467_dp
        sep%aq   =    0.5_dp/0.652785_dp
        sep%fn1(1) = 0.048906_dp
        sep%fn2(1) = 5.765340_dp
        sep%fn3(1) = 1.682232_dp
        sep%fn1(2) = 0.047697_dp
        sep%fn2(2) = 5.973721_dp
        sep%fn3(2) = 0.894406_dp
        sep%pre(1) = -0.000743_dp
        sep%pre(2) = 0.000985_dp
        sep%d(1) = 0.836915_dp
        sep%d(2) = 1.585236_dp
      CASE (7)
      !Nitrogen
        sep%defined = .TRUE.
        sep%uss   =  -49.454546_dp
        sep%upp   =  -47.757406_dp
        sep%beta(0) =  -14.117230_dp
        sep%beta(1)=  -19.938509_dp
        sep%sto_exponents(0)    =    2.035807_dp
        sep%sto_exponents(1)   =    2.324327_dp
        sep%alp   =    2.849124_dp
        sep%eisol = -158.416205_dp
        sep%gss   =   11.9047870_dp
        sep%gsp  =    7.3485650_dp
        sep%gpp  =   11.7546720_dp
        sep%gp2  =   10.8072770_dp
        sep%hsp  =    1.1367130_dp
        sep%dd   =    0.654855_dp
        sep%qq    =    0.526924_dp
        sep%am    =    0.5_dp/1.142818_dp
        sep%ad    =    0.5_dp/0.991235_dp
        sep%aq   =    0.5_dp/0.676704_dp
        sep%fn1(1) =  1.513320_dp
        sep%fn2(1) = 5.904394_dp
        sep%fn3(1) = 1.728376_dp
        sep%fn1(2) = -1.511892_dp
        sep%fn2(2) = 6.030014_dp
        sep%fn3(2) = 1.734108_dp
        sep%pre(1) = -0.003160_dp
        sep%pre(2) = 0.012501_dp
        sep%d(1) = 1.004172_dp
        sep%d(2) = 1.516336_dp
      CASE (8)
      !Oxygen
        sep%defined = .TRUE.
        sep%uss   =  -87.412505_dp
        sep%upp   =  -72.183070_dp
        sep%beta(0) =  -44.874553_dp
        sep%beta(1)=  -24.601939_dp
        sep%sto_exponents(0)    =    3.814565_dp
        sep%sto_exponents(1)   =    2.318011_dp
        sep%alp   =    3.225309_dp
        sep%eisol = -292.188766_dp
        sep%gss   =   15.7557600_dp
        sep%gsp  =   10.6211600_dp
        sep%gpp  =   13.6540160_dp
        sep%gp2  =   12.4060950_dp
        sep%hsp  =    0.5938830_dp
        sep%dd   =    0.403741_dp
        sep%qq    =    0.528360_dp
        sep%am    =    0.5_dp/0.863494_dp
        sep%ad    =    0.5_dp/0.936266_dp
        sep%aq   =    0.5_dp/0.624291_dp
        sep%fn1(1) = -1.138455_dp
        sep%fn2(1) = 6.000043_dp
        sep%fn3(1) = 1.622362_dp
        sep%fn1(2) =  1.146007_dp
        sep%fn2(2) = 5.963494_dp
        sep%fn3(2) = 1.614788_dp
        sep%pre(1) = -0.001000_dp
        sep%pre(2) = -0.0015220_dp
        sep%d(1) = 1.360685_dp
        sep%d(2) = 1.366407_dp
      CASE (9)
        !Fluorine
        sep%defined = .TRUE.
        sep%uss   = -110.4004320_dp
        sep%upp   = -106.3952640_dp
        sep%beta(0) =  -50.9373010_dp
        sep%beta(1)=  -31.6369760_dp
        sep%sto_exponents(0)    =    5.5380330_dp
        sep%sto_exponents(1)   =    2.5380660_dp
        sep%alp   =    3.2005710_dp
        sep%eisol = -442.4571330_dp
        sep%gss   =   10.4966670_dp
        sep%gsp  =   16.0736890_dp
        sep%gpp  =   14.8172560_dp
        sep%gp2  =   14.4183930_dp
        sep%hsp  =    0.7277630_dp
        sep%dd   =    0.246601_dp
        sep%qq    =    0.482551_dp
        sep%am    =    0.5_dp/1.296126_dp
        sep%ad    =    0.5_dp/0.6346333_dp
        sep%aq   =    0.5_dp/0.805802_dp
        sep%fn1(1) = -0.008079_dp
        sep%fn2(1) = 5.9389690_dp
        sep%fn3(1) = 1.863949_dp
        sep%fn1(2) = -0.0026590_dp
        sep%fn2(2) = 5.9251050_dp
        sep%fn3(2) = 2.3888640_dp
        sep%pre(1) = -0.012866_dp
        sep%pre(2) = -0.007315_dp
        sep%d(1) = 1.305681_dp
        sep%d(2) = 1.842572_dp
      CASE (10:13)
      CASE (14)
      !Silicon
        sep%defined = .TRUE.
        sep%uss   =  -26.332522_dp
        sep%upp   =  -22.602540_dp
        sep%beta(0) =  -3.376445_dp
        sep%beta(1)=  -3.6209690_dp
        sep%sto_exponents(0)    =    1.586389_dp
        sep%sto_exponents(1)   =    1.485958_dp
        sep%alp   =    2.215157_dp
        sep%eisol = -66.839000_dp
        sep%gss   =    5.0471960_dp
        sep%gsp  =    5.9490570_dp
        sep%gpp  =    6.7593670_dp
        sep%gp2  =    5.1612970_dp
        sep%hsp  =    0.9198320_dp
        sep%dd   =    1.310515_dp
        sep%qq    =    1.126089_dp
        sep%am    =    0.5_dp/2.695556_dp
        sep%ad    =    0.5_dp/1.630757_dp
        sep%aq   =    0.5_dp/0.949200_dp
        sep%fn1(1) = -0.071314_dp
        sep%fn2(1) = 6.000000_dp
        sep%fn3(1) = 0.237995_dp
        sep%fn1(2) =  0.089451_dp
        sep%fn2(2) = 6.000000_dp
        sep%fn3(2) = 1.897728_dp
        sep%pre(1) = -0.091928_dp
        sep%pre(2) = -0.040753_dp
        sep%d(1) = 1.163190_dp
        sep%d(2) = 2.190526_dp
      CASE (15)
      !Phospherus
        sep%defined = .TRUE.
        sep%uss   =  -37.882113_dp
        sep%upp   =  -30.312979_dp
        sep%beta(0) =  -12.676297_dp
        sep%beta(1)=  -7.093318_dp
        sep%sto_exponents(0)    =    2.395882_dp
        sep%sto_exponents(1)   =    1.742213_dp
        sep%alp   =    2.005294_dp
        sep%eisol = -117.212854_dp
        sep%gss   =    7.8016150_dp
        sep%gsp  =    5.1869490_dp
        sep%gpp  =    6.6184780_dp
        sep%gp2  =    6.0620020_dp
        sep%hsp  =    1.5428090_dp
        sep%dd   =    0.893978_dp
        sep%qq    =    0.960457_dp
        sep%am    =    0.5_dp/1.743870_dp
        sep%ad    =    0.5_dp/1.050851_dp
        sep%aq   =    0.5_dp/1.208907_dp
        sep%fn1(1) = -0.398055_dp
        sep%fn2(1) = 1.997272_dp
        sep%fn3(1) = 0.950073_dp
        sep%fn1(2) = -0.079653_dp
        sep%fn2(2) = 1.998360_dp
        sep%fn3(2) = 2.336959_dp
        sep%pre(1) = 0.462741_dp
        sep%pre(2) = -0.020444_dp
        sep%d(1) = 0.714296_dp
        sep%d(2) = 2.041209_dp
      CASE (16)
      !Sulfur
        sep%defined = .TRUE.
        sep%uss   =  -43.906366_dp
        sep%upp   =  -43.461348_dp
        sep%beta(0) =  -2.953912_dp
        sep%beta(1)=  -8.507779_dp
        sep%sto_exponents(0)    =    1.012002_dp
        sep%sto_exponents(1)   =    1.876999_dp
        sep%alp   =    2.539751_dp
        sep%eisol = -166.336554_dp
        sep%gss   =    8.9646670_dp
        sep%gsp  =    6.7859360_dp
        sep%gpp  =    9.9681640_dp
        sep%gp2  =    7.9702470_dp
        sep%hsp  =    4.0418360_dp
        sep%dd   =    1.006989_dp
        sep%qq    =    0.891487_dp
        sep%am    =    0.5_dp/1.517625_dp
        sep%ad    =    0.5_dp/0.711672_dp
        sep%aq   =    0.5_dp/0.754336_dp
        sep%fn1(1) = -0.330692_dp
        sep%fn2(1) = 6.000000_dp
        sep%fn3(1) = 0.823837_dp
        sep%fn1(2) =  0.024171_dp
        sep%fn2(2) = 6.000000_dp
        sep%fn3(2) = 2.017756_dp
        sep%pre(1) = 0.120434_dp
        sep%pre(2) = -0.002663_dp
        sep%d(1) = 0.672870_dp
        sep%d(2) = 2.032340_dp
      CASE (17)
      !Chlorine
        sep%defined = .TRUE.
        sep%uss   = -95.094434_dp
        sep%upp   =  -53.921651_dp
        sep%beta(0) =  -26.913129_dp
        sep%beta(1)=  -14.991178_dp
        sep%sto_exponents(0)    =    2.5482680_dp
        sep%sto_exponents(1)   =    2.2846240_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    2.4976170_dp
        sep%eisol = -305.7152010_dp
        sep%gss   =   16.0136010_dp
        sep%gsp  =    8.0481150_dp
        sep%gpp  =    7.5222150_dp
        sep%gp2  =    7.5041540_dp
        sep%hsp  =    3.4811530_dp
        sep%dd   =    0.827561_dp
        sep%qq    =    0.7324270_dp
        sep%am    =    0.5_dp/0.8495900_dp
        sep%ad    =    0.5_dp/0.696164_dp
        sep%aq   =    0.5_dp/2.299104_dp
        sep%fn1(1) = -0.112222_dp
        sep%fn2(1) = 5.9637190_dp
        sep%fn3(1) = 1.0277190_dp
        sep%fn1(2) = -0.0130610_dp
        sep%fn2(2) = 1.9995560_dp
        sep%fn3(2) = 2.2863770_dp
        sep%pre(1) = -0.016552_dp
        sep%pre(2) = -0.016646_dp
        sep%d(1) = 1.727690_dp
        sep%d(2) = 1.784655_dp
      CASE (18:34)
      CASE (35)
        !Bromine
        sep%defined = .TRUE.
        sep%uss   = -115.8419630_dp
        sep%upp   =  -74.2051460_dp
        sep%beta(0) =  -21.5380440_dp
        sep%beta(1)=   -8.5247640_dp
        sep%sto_exponents(0)    =    4.3450790_dp
        sep%sto_exponents(1)   =    2.1909610_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    2.4246730_dp
        sep%eisol = -351.0138870_dp
        sep%gss   =   15.9434250_dp
        sep%gsp  =   16.0616800_dp
        sep%gpp  =    8.2827630_dp
        sep%gp2  =    7.8168490_dp
        sep%hsp  =    0.5788690_dp
        sep%dd   =    0.473860_dp
        sep%qq    =    0.968214_dp
        sep%am    =    0.5_dp/0.853330_dp
        sep%ad    =    0.5_dp/1.046430_dp
        sep%aq   =    0.5_dp/1.280643_dp
        sep%fn1(1) =  0.9613620_dp
        sep%fn2(1) = 6.0136000_dp
        sep%fn3(1) = 2.3404450_dp
        sep%fn1(2) = -0.9548834_dp
        sep%fn2(2) = 5.9763290_dp
        sep%fn3(2) = 2.3487450_dp
        sep%pre(1) = -0.013772_dp
        sep%pre(2) =  0.008849_dp
        sep%d(1) = 1.852030_dp
        sep%d(2) = 2.338958_dp
      CASE (36:52)
      CASE (53)
        !Iodine
        sep%defined = .TRUE.
        sep%uss   =  -96.6641740_dp
        sep%upp   =  -61.1671370_dp
        sep%beta(0) =  -16.5926210_dp
        sep%beta(1)=   -6.5998160_dp
        sep%sto_exponents(0)    =    5.0628010_dp
        sep%sto_exponents(1)   =    2.4177570_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%alp   =    1.9781700_dp
        sep%eisol = -291.5378690_dp
        sep%gss   =   13.6319430_dp
        sep%gsp  =   14.9904060_dp
        sep%gpp  =    7.2883300_dp
        sep%gp2  =    5.9664070_dp
        sep%hsp  =    2.6300350_dp
        sep%dd   =    0.407261_dp
        sep%qq    =    1.062574_dp
        sep%am    =    0.5_dp/0.998024_dp
        sep%ad    =    0.5_dp/0.532290_dp
        sep%aq   =    0.5_dp/0.979783_dp
        sep%fn1(1) = -0.136003_dp
        sep%fn2(1) = 3.8529120_dp
        sep%fn3(1) = 1.6974550_dp
        sep%fn1(2) = -0.0372870_dp
        sep%fn2(2) = 5.2292640_dp
        sep%fn3(2) = 2.7686690_dp
        sep%pre(1) =  0.012901_dp
        sep%pre(2) = -0.012825_dp
        sep%d(1) = 1.994299_dp
        sep%d(2) = 2.263417_dp
      CASE (54:106)
    END SELECT

    sep%beta = sep%beta/evolt
    sep%uss = sep%uss/evolt; sep%upp = sep%upp/evolt; sep%udd = sep%udd/evolt
    sep%alp = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss = sep%gss/evolt; sep%gsp = sep%gsp/evolt; sep%gpp = sep%gpp/evolt; sep%gp2 = sep%gp2/evolt 
    sep%gsd = sep%gsd/evolt; sep%gpd = sep%gpd/evolt; sep%gdd = sep%gdd/evolt
    sep%hsp = sep%hsp/evolt
    sep%fn1 = sep%fn1*bohr/evolt; sep%fn2 = sep%fn2/bohr/bohr; sep%fn3 = sep%fn3*bohr
    sep%pre = sep%pre/evolt; sep%d = sep%d/bohr

  END SUBROUTINE pdg_default_parameter
! *****************************************************************************
  
  SUBROUTINE valence_electrons ( z, zeff, natorb )
    INTEGER, INTENT(IN)                      :: z
    REAL(KIND=dp)                            :: zeff
    INTEGER, INTENT(OUT)                     :: natorb

    INTEGER, DIMENSION(0:nelem), PARAMETER :: norb = (/ -1, 1,  1,  4,  4,  4,&
      4,  4,  4,  4,  4,  0,  4,  4,  4,  4,  4,  4,  4,0,  4,  9,  9,  9,  9,&
      9,  9,  9,  9,  9,  4,  4,  4,  4,  4,  4,  4,4,  4,  9,  9,  9,  9,  9,&
      9,  9,  9,  9,  4,  4,  4,  4,  4,  4,  4,2,  2,  8,  8,  8,  8,  8,  8,&
      8,  8,  8,  8,  8,  8,  8,  8,  9,  9,9,  9,  9,  9,  9,  9,  9,  4,  4,&
      4,  4,  4,  4, -1, -1, -1, -1,  9,-1,  9, -1, -1, -1, -1, -1, -1, -1, -1&
      , -1, -1, -1, -1, -1, -1         /), zval = (/ -1, 1,  0,  1,  2,  3,  4&
      ,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,  0,1,  2,  3,  4,  5,  6, &
      7,  8,  9, 10, 11,  2,  3,  4,  5,  6,  7,  0,1,  2,  3,  4,  5,  6,  7,&
      8,  9, 10, 11,  2,  3,  4,  5,  6,  7,  0,1,  2,  3,  4,  5,  6,  7,  8,&
      9, 10, 11, 12, 13, 14, 15, 16,  3,  4,5,  6,  7,  8,  9, 10, 11,  2,  3,&
      4,  5,  6,  7, -1, -1, -1, -1,  4,-1,  6, -1, -1, -1, -1, -1, -1, -1, -1&
      , -1, -1, -1, -1, -1, -1         /)

    zeff = REAL(zval(z),dp)
    natorb = norb(z)

  END SUBROUTINE valence_electrons 

! *****************************************************************************

  SUBROUTINE slater_exponents (z,norb,sto,itype) 
    INTEGER, INTENT(IN)                      :: z
    INTEGER, DIMENSION(0:3), INTENT(IN)      :: norb
    REAL(dp), DIMENSION(:, 0:)               :: sto
    INTEGER, INTENT(IN)                      :: itype

    REAL(dp), DIMENSION(18), PARAMETER :: crtab1s = (/ 1.0000, 1.6875, 2.6906,&
      3.6848, 4.6795, 5.6727, 6.6651,7.6579, 8.6501, 9.6421,10.6259,11.6089, &
      12.5910,13.5745,14.5578,15.5409,16.5239,17.5075 /), crtab2p = (/ 0.0000,&
      0.0000, 0.0000, 0.0000, 1.2107, 1.5679, 1.9170,2.2266, 2.5500, 2.8792, &
      3.4009, 3.9129, 4.4817, 4.9725,5.4806, 5.9885, 6.4966, 7.0041 /), &
      crtab2s = (/ 0.0000, 0.0000, 0.6396, 0.9560, 1.2881, 1.6083, 1.9237,&
      2.2458, 2.5638, 2.8792,3.2857, 3.6960, 4.1068, 5.5100,4.9125, 5.3144, &
      5.7152, 6.1152 /), crtab3p = (/ 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
      0.0000, 0.0000,0.0000, 0.0000, 0.0000,0.0000, 0.0000, 1.3552, 1.4284,&
      1.6288, 1.8273, 2.0387, 2.2547 /) , crtab3s = (/ 0.0000, 0.0000, 0.0000,&
      0.0000, 0.0000, 0.0000, 0.0000,0.0000, 0.0000, 0.0000,0.8358, 1.1025, &
      1.3724, 1.6344,1.8806, 2.1223, 2.3561, 2.5856 /)

    INTEGER                                  :: i, j, l, nj
    INTEGER, DIMENSION(1:4)                  :: no
    INTEGER, DIMENSION(1:4, 1:7)             :: ne

!
!   Slater exponents according to
!   itype = 0 : Clementi-Raimondi if available, Slater rules for others
!   itype = 1 : Slater rules
!   itype = 2 : Clementi-Raimondi
!
!   Clementi-Raimondi Table

    sto = 0._dp
    ne = 0
    DO l=1,4
      nj = 2*(l-1)+1
      no(l) = l-1
      DO i=l,7
        ne(l,i) = ptable(z)%e_conv(l-1) - 2*nj*(i-l)
        ne(l,i) = MAX(ne(l,i),0)
        ne(l,i) = MIN(ne(l,i),2*nj)
        IF(ne(l,i) > 0) no(l)=no(l)+1
      END DO
    END DO
    DO l=0,3
      DO i=1,norb(l)
        j=no(l+1)-norb(l)+1
        IF ( itype == 0 .OR. itype == 1 ) THEN
          sto(i,l)=srules(z,ne,j,l)
        END IF
        IF ( itype == 0 .OR. itype == 2 ) THEN
          IF ( z <= 18 .AND. z > 0 ) THEN
            IF ( l == 0 ) THEN
              IF ( j == 1 ) sto(i,l)=crtab1s(z)
              IF ( j == 2 ) sto(i,l)=crtab2s(z)
              IF ( j == 3 ) sto(i,l)=crtab3s(z)
            ELSE IF ( l == 1 ) THEN
              IF ( j == 2 ) sto(i,l)=crtab2p(z)
              IF ( j == 3 ) sto(i,l)=crtab3p(z)
            END IF
          END IF
        END IF
      END DO
    END DO
    
  END SUBROUTINE slater_exponents  

! *****************************************************************************

  FUNCTION srules(z,ne,n,l)
    ! Slater rules
    INTEGER                                  :: z
    INTEGER, DIMENSION(4, 7)                 :: ne
    INTEGER                                  :: n, l
    REAL(dp)                                 :: srules

    REAL(dp), DIMENSION(7), PARAMETER :: &
      xns = (/ 1.0_dp,2.0_dp,3.0_dp,3.7_dp,4.0_dp,4.2_dp,4.4_dp /)

    INTEGER                                  :: i, l1, l2, m, m1, m2
    REAL(dp)                                 :: s

    s = 0.0_dp
    ! The complete shell
    l1=l+1
    IF(l1 == 1) l2=2
    IF(l1 == 2) l2=1
    IF(l1 == 3) l2=4
    IF(l1 == 4) l2=3
    ! Rule a) no contribution from shells further out
    ! Rule b) 0.35 (1s 0.3) from each other electron in the same shell
    IF(n == 1) THEN
       m=ne(1,1)
       s=s+0.3_dp*REAL(m-1,dp)
    ELSE
       m=ne(l1,n)+ne(l2,n)
       s=s+0.35_dp*REAL(m-1,dp)
    END IF
    ! Rule c) if (s,p) shell 0.85 from each electron with n-1, and 1.0
    !      from all electrons further in
    IF(l1+l2 == 3) THEN
       IF(n > 1) THEN
          m1=ne(1,n-1)+ne(2,n-1)+ne(3,n-1)+ne(4,n-1)
          m2=0
          DO i=1,n-2
             m2=m2+ne(1,i)+ne(2,i)+ne(3,i)+ne(4,I)
          END DO
          s=s+0.85_dp*REAL(m1,dp)+1._dp*REAL(m2,dp)
       END IF
    ELSE
       ! Rule d) if (d,f) shell 1.0 from each electron inside
       m=0
       DO i=1,n-1
          m=m+ne(1,i)+ne(2,i)+ne(3,i)+ne(4,i)
       END DO
       s=s+1._dp*REAL(m,dp)
    END IF
    ! Slater exponent is (Z-S)/NS
    srules = (REAL(z,dp) - s)/xns(n)

  END FUNCTION srules

! *****************************************************************************

  SUBROUTINE calpar(z,sep)
    INTEGER                                  :: z
    TYPE(semi_empirical_type), POINTER       :: sep

    CHARACTER(len=*), PARAMETER :: routineN = 'calpar', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: j
    INTEGER, DIMENSION(107), PARAMETER :: nspqn = (/(1,j=1,2), (2,j=1,8), (3,j&
      =1,8), (4,j=1,18), (5,j=1,18), (6,j=1,32), (0,j=1,21)/)

    INTEGER                                  :: iod, iop, ios, jmax, k, l
    REAL(KIND=dp) :: ad, am, aq, d1, d2, d3, dd, df, eisol, gdd1, gp2, gp2c, &
      gpp, gppc, gqq, gsp, gspc, gss, gssc, hpp, hpp1, hpp2, hsp, hsp1, hsp2, &
      hspc, p, p4, q1, q2, q3, qf, qn, qq, udd, upp, uss, zp, zs

    IF (.NOT.sep%defined) RETURN
    uss = sep%uss
    upp = sep%upp
    udd = sep%udd
    gss = sep%gss
    gpp = sep%gpp
    gsp = sep%gsp
    gp2 = sep%gp2
    hsp = sep%hsp
    zs  = sep%sto_exponents(0)
    zp  = sep%sto_exponents(1)
    ios = Nos(z)
    iop = Nop(z)
    iod = Nod(z)

    p = 2.e0_dp 
    p4 = p**4
    !  GSSC is the number of two-electron terms of type <SS|SS>
    gssc  = REAL(MAX(ios-1,0),KINd=dp)
    k = iop
    !  GSPC is the number of two-electron terms of type <SS|PP>
    gspc = REAL(ios*k,KIND=dp) 
    l = MIN(k,6-k) 
    !  GP2C is the number of two-electron terms of type <PP|PP>
    !       plus 0.5 of the number of HPP integrals.
    !  (HPP is not used; instead it is replaced by 0.5(GPP-GP2))
    gp2c = REAL((k*(k - 1))/2 + 0.5e0_dp*(l*(l - 1))/2,KIND=dp)
    !  GPPC is minus 0.5 times the number of HPP integrals.
    gppc = REAL(-0.5e0_dp*(l*(l - 1))/2,KIND=dp)
    !  HSPC is the number of two-electron terms of type <SP|SP>.
    !       (S and P must have the same spin.  In all cases, if
    !  P is non-zero, there are two S electrons)
    hspc = REAL(-k ,KIND=dp)
    !   CONSTRAINTS ON THE POSSIBLE VALUES OF PARAMETERS
    zp = MAX(0.3e0_dp,zp) 
    !  PUT IN ANY CONSTRAINTS AT THIS POINT
    hpp = 0.5e0_dp*(gpp-gp2) 
    hpp = MAX(0.1e0_dp,hpp) 
    hsp = MAX(1.e-7_dp,hsp) 
    eisol = uss*ios + upp*iop + udd*iod + gss*gssc + gpp*gppc + gsp*gspc + gp2*gp2c + hsp*hspc 
    qn = REAL(nspqn(z),KIND=dp)
    dd = (2.e0_dp*qn + 1)*(4.e0_dp*zs*zp)**(qn + 0.5e0_dp)/(zs+zp)**(2.e0_dp*qn + 2)/SQRT(3.e0_dp) 
    qq = SQRT((4.e0_dp*qn*qn + 6.e0_dp*qn + 2.e0_dp)/20.e0_dp)/zp
    !     CALCULATE ADDITIVE TERMS, IN ATOMIC UNITS.
    jmax = 5 
    gdd1 = (hsp/(evolt*dd**2))**(1.e0_dp/3.e0_dp) 
    d1 = gdd1 
    d2 = gdd1 + 0.04e0_dp 
    DO j = 1, jmax 
       df = d2 - d1 
       hsp1 = 0.5e0_dp*d1 - 0.5e0_dp/SQRT(4.e0_dp*dd**2+1.e0_dp/d1**2) 
       hsp2 = 0.5e0_dp*d2 - 0.5e0_dp/SQRT(4.e0_dp*dd**2+1.e0_dp/d2**2) 
       IF (ABS(hsp2 - hsp1) < 1.e-25_dp) EXIT  
       d3 = d1 + df*(hsp/evolt-hsp1)/(hsp2 - hsp1) 
       d1 = d2 
       d2 = d3 
    END DO
    gqq = (p4*hpp/(evolt*48.e0_dp*qq**4))**0.2e0_dp 
    q1 = gqq 
    q2 = gqq + 0.04e0_dp 
    DO j = 1, jmax 
       qf = q2 - q1 
       hpp1 = 0.25e0_dp*q1 - 0.5e0_dp/SQRT(4.e0_dp*qq**2+1.e0_dp/q1**2) + 0.25e0_dp/SQRT(8.e0_dp*qq**2+1.e0_dp/q1**2) 
       hpp2 = 0.25e0_dp*q2 - 0.5e0_dp/SQRT(4.e0_dp*qq**2+1.e0_dp/q2**2) + 0.25e0_dp/SQRT(8.e0_dp*qq**2+1.e0_dp/q2**2) 
       IF (ABS(hpp2 - hpp1) < 1.e-25_dp) EXIT  
       q3 = q1 + qf*(hpp/evolt - hpp1)/(hpp2 - hpp1) 
       q1 = q2 
       q2 = q3 
    END DO
    am = gss/evolt 
    ad = d2 
    aq = q2 
    IF (z==1) THEN
       ad = am
       aq = am
       dd = 0.0_dp
       qq = 0.0_dp
    END IF
    sep%eisol = eisol
    sep%dd = dd
    sep%qq = qq
    sep%am = am
    sep%ad = ad
    sep%aq = aq
  END SUBROUTINE calpar
  
END MODULE semi_empirical_parameters

! *****************************************************************************
