!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 CP2K developers group                                  !
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
!!     fn1            none      AM1/PM3 specific parameters
!!     fn2            A**(-2)   for the core-core repulsion
!!     fn3            A         interaction.                       
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
!!     BK  (12.10.2004)
!!
!!   SOURCE
!******************************************************************************

MODULE semi_empirical_parameters

  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE periodic_table,                  ONLY: ptable, nelem
  USE physcon,                         ONLY: evolt, kcalmol, bohr
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "semi_empirical_parameters"

! *** Public subroutines ***

  PUBLIC :: eht_default_parameter, &
            mndo_default_parameter, &
            am1_default_parameter, &
            pm3_default_parameter

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE eht_default_parameter ( sep, z )

    TYPE(semi_empirical_type),POINTER    :: sep
    INTEGER                              :: z

    CHARACTER(LEN=*), PARAMETER :: routine_name = "eht_default_parameter"

    INTEGER, DIMENSION(0:3)              :: norb
    REAL(dp), DIMENSION(1,0:3)           :: alpha

! -----------------------------------------------------------------------------

    sep%typ = "EHT"
    sep%name = ptable(z)%name
    sep%z = z
    CALL valence_electrons ( z, sep%zeff, sep%natorb )
    norb=0
    IF(sep%natorb >= 1) norb(0)=1
    IF(sep%natorb >= 4) norb(1)=1
    IF(sep%natorb >= 9) norb(2)=1
    IF(sep%natorb >=16) norb(3)=1
    CALL slater_exponents ( z, norb, alpha, 1 ) 
    sep%sto_exponents(0:3) = alpha(1,0:3)
    sep%beta(0:3) = ptable(z)%eht_param(0:3)/evolt
    IF ( ABS ( SUM (sep%beta(0:3)) ) > 0._dp ) sep%defined = .TRUE.
    sep%ngauss = 6
    sep%eheat = ptable(z)%heat_of_formation / kcalmol

  END SUBROUTINE eht_default_parameter

! *****************************************************************************

  SUBROUTINE mndo_default_parameter ( sep, z )

    TYPE(semi_empirical_type),POINTER    :: sep
    INTEGER                              :: z

    CHARACTER(LEN=*), PARAMETER :: routine_name = "mndo_default_parameter"

! -----------------------------------------------------------------------------

    sep%typ = "MNDO"
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

    SELECT CASE (z)
      CASE DEFAULT
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "invalid atomic number")
      CASE (0)
        sep%defined = .TRUE.
      CASE (1)
        !Hydrogen
        sep%defined          = .TRUE.
        sep%alp              =   2.5441341_dp
        sep%eisol            = -11.9062760_dp
        sep%beta(0)          =  -6.9890640_dp
        sep%sto_exponents(0) =   1.3319670_dp
        sep%am               =   0.4721793_dp
        sep%ad               =   0.4721793_dp 
        sep%aq               =   0.4721793_dp
        sep%uss              = -11.9062760_dp
        sep%gss              =  12.8480000_dp
      CASE (2)
      CASE (3)
        !Lithium
        sep%defined          = .TRUE.
        sep%alp              =   1.2501400_dp
        sep%eisol            =  -5.1280000_dp
        sep%beta(0)          =  -1.3500400_dp
        sep%beta(1)          =  -1.3500400_dp
        sep%sto_exponents(0) =   0.7023800_dp 
        sep%sto_exponents(1) =   0.7023800_dp
        sep%dd               =   2.0549783_dp 
        sep%qq               =   1.7437069_dp
        sep%am               =   0.2682837_dp 
        sep%ad               =   0.2269793_dp 
        sep%aq               =   0.2614581_dp
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
        sep%eisol            = -24.2047560_dp
        sep%beta(0) = -4.0170960_dp
        sep%beta(1) = -4.0170960_dp
        sep%sto_exponents(0) = 1.0042100_dp 
        sep%sto_exponents(1) = 1.0042100_dp
        sep%dd = 1.4373245_dp 
        sep%qq = 1.2196103_dp 
        sep%am = 0.3307607_dp 
        sep%ad = 0.3356142_dp 
        sep%aq = 0.3846373_dp
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
        sep%eisol = -64.3159500_dp
        sep%beta(0) =  -8.2520540_dp
        sep%beta(1) =  -8.2520540_dp
        sep%sto_exponents(0)    =   1.5068010_dp
        sep%sto_exponents(1)    =   1.5068010_dp
        sep%dd    =   0.9579073_dp
        sep%qq    =   0.8128113_dp
        sep%am    =   0.3891951_dp
        sep%ad    =   0.4904730_dp
        sep%aq    =   0.5556979_dp
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
        sep%eisol = -120.5006060_dp
        sep%beta(0) =  -18.9850440_dp
        sep%beta(1) =   -7.9341220_dp
        sep%sto_exponents(0)    =    1.7875370_dp
        sep%sto_exponents(1)    =    1.7875370_dp
        sep%dd    =    0.8074662_dp
        sep%qq    =    0.6851578_dp
        sep%am    =    0.4494671_dp
        sep%ad    =    0.6149474_dp
        sep%aq    =    0.6685897_dp
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
        sep%eisol = -202.5812010_dp
        sep%beta(0) =  -20.4957580_dp
        sep%beta(1) =  -20.4957580_dp
        sep%sto_exponents(0)    =    2.2556140_dp
        sep%sto_exponents(1)    =    2.2556140_dp
        sep%dd    =    0.6399037_dp
        sep%qq    =    0.5429763_dp
        sep%am    =    0.4994487_dp
        sep%ad    =    0.7843643_dp
        sep%aq    =    0.8144720_dp
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
        sep%eisol = -317.8685060_dp
        sep%beta(0) =  -32.6880820_dp
        sep%beta(1) =  -32.6880820_dp
        sep%sto_exponents(0)    =    2.6999050_dp
        sep%sto_exponents(1)    =    2.6999050_dp
        sep%dd    =    0.5346024_dp
        sep%qq    =    0.4536252_dp
        sep%am    =    0.5667034_dp
        sep%ad    =    0.9592562_dp
        sep%aq    =    0.9495934_dp
        sep%uss   =  -99.6443090_dp
        sep%upp   =  -77.7974720_dp
        sep%gss   =   15.4200000_dp
        sep%gsp   =   14.4800000_dp
        sep%gpp   =   14.5200000_dp
        sep%gp2   =   12.9800000_dp
        sep%hsp   =    3.9400000_dp
      CASE (9)
      CASE (10)
      !Fluorine
        sep%defined = .TRUE.
        sep%alp   =    3.4196606_dp
        sep%eisol = -476.6837810_dp
        sep%beta(0) =  -48.2904660_dp
        sep%beta(1) =  -36.5085400_dp
        sep%sto_exponents(0)    =    2.8484870_dp
        sep%sto_exponents(1)    =    2.8484870_dp
        sep%dd    =    0.5067166_dp
        sep%qq    =    0.4299633_dp
        sep%am    =    0.6218302_dp
        sep%ad    =    1.0850301_dp
        sep%aq    =    1.0343643_dp
        sep%uss   = -131.0715480_dp
        sep%upp   = -105.7821370_dp
        sep%gss   =   16.9200000_dp
        sep%gsp   =   17.2500000_dp
        sep%gpp   =   16.7100000_dp
        sep%gp2   =   14.9100000_dp
        sep%hsp   =    4.8300000_dp
      CASE (11)
      !Sodium
        sep%defined = .TRUE.
        sep%alp   =  1.66_dp
        sep%eisol =  0.00_dp
        sep%am    =  0.50_dp
      CASE (12)
      CASE (13)
       !Aluminium
        sep%defined = .TRUE.
        sep%alp   =   1.8688394_dp
        sep%eisol = -44.4840720_dp
        sep%beta(0) =  -2.6702840_dp
        sep%beta(1) =  -2.6702840_dp
        sep%sto_exponents(0) =   1.4441610_dp
        sep%sto_exponents(1) =   1.4441610_dp
        sep%sto_exponents(2) =   1.0000000_dp
        sep%dd    =   1.3992387_dp
        sep%qq    =   1.1586797_dp
        sep%am    =   0.2973172_dp
        sep%ad    =   0.2635574_dp
        sep%aq    =   0.3673560_dp
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
        sep%eisol = -82.8394220_dp
        sep%beta(0) =  -9.0868040_dp
        sep%beta(1) =  -1.0758270_dp
        sep%sto_exponents(0)    =   1.3159860_dp
        sep%sto_exponents(1)    =   1.7099430_dp
        sep%sto_exponents(2)    =   1.0000000_dp
        sep%dd    =   1.2580349_dp
        sep%qq    =   0.9785824_dp
        sep%am    =   0.3608967_dp
        sep%ad    =   0.3664244_dp
        sep%aq    =   0.4506740_dp
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
        sep%eisol = -152.9599600_dp
        sep%beta(0) =   -6.7916000_dp
        sep%beta(1) =   -6.7916000_dp
        sep%sto_exponents(0)    =    2.1087200_dp
        sep%sto_exponents(1)    =    1.7858100_dp
        sep%sto_exponents(2)    =    1.0000000_dp
        sep%dd    =    1.0129699_dp
        sep%qq    =    0.9370090_dp
        sep%am    =    0.4248438_dp
        sep%ad    =    0.4882420_dp
        sep%aq    =    0.4979406_dp
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
        sep%eisol = -226.0123900_dp
        sep%beta(0) =  -10.7616700_dp
        sep%beta(1) =  -10.1084330_dp
        sep%sto_exponents(0)    =    2.3129620_dp
        sep%sto_exponents(1)    =    2.0091460_dp
        sep%sto_exponents(2)    =    1.0000000_dp
        sep%dd    =    0.9189935_dp
        sep%qq    =    0.8328514_dp
        sep%am    =    0.4733554_dp
        sep%ad    =    0.5544502_dp
        sep%aq    =    0.5585244_dp
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
        sep%eisol = -353.1176670_dp
        sep%beta(0) =  -14.2623200_dp
        sep%beta(1) =  -14.2623200_dp
        sep%sto_exponents(0)    =    3.7846450_dp
        sep%sto_exponents(1)    =    2.0362630_dp
        sep%sto_exponents(2)    =    1.0000000_dp
        sep%dd    =    0.4986870_dp
        sep%qq    =    0.8217603_dp
        sep%am    =    0.5523705_dp
        sep%ad    =    0.8061220_dp
        sep%aq    =    0.6053435_dp
        sep%uss   = -100.2271660_dp
        sep%upp   =  -77.3786670_dp
        sep%gss   =   15.0300000_dp
        sep%gsp   =   13.1600000_dp
        sep%gpp   =   11.3000000_dp
        sep%gp2   =    9.9700000_dp
        sep%hsp   =    2.4200000_dp
      CASE (18)
      CASE (19)
        !Potassium
        sep%defined = .TRUE.
        sep%alp   =  1.16_dp
        sep%eisol =  0.00_dp
        sep%am    =  0.50_dp
      CASE (20:23)
      CASE (24)
        !Chromium.
        sep%defined = .TRUE.
        sep%alp   =    3.0683070_dp
        sep%eisol = -134.8187920_dp
        sep%beta(0) =   -0.1000000_dp
        sep%beta(1) =   -0.1000000_dp
        sep%beta(2) =   -8.7766360_dp
        sep%sto_exponents(0)    =    1.5000000_dp
        sep%sto_exponents(1)    =    1.5000000_dp
        sep%sto_exponents(2)    =    2.8845490_dp
        sep%dd    =    1.7320508_dp
        sep%qq    =    1.4142136_dp
        sep%am    =    0.2205072_dp
        sep%ad    =    0.2711332_dp
        sep%aq    =    0.4464656_dp
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
        sep%eisol =  -29.8794320_dp
        sep%beta(0) =   -1.0000000_dp
        sep%beta(1) =   -2.0000000_dp
        sep%sto_exponents(0)    =    2.0473590_dp
        sep%sto_exponents(1)    =    1.4609460_dp
        sep%sto_exponents(2)    =    1.0000000_dp
        sep%dd    =    1.3037826_dp
        sep%qq    =    1.4520183_dp
        sep%am    =    0.4336641_dp
        sep%ad    =    0.2375912_dp
        sep%aq    =    0.2738858_dp
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
        sep%eisol = -76.2489440_dp
        sep%beta(0) =  -4.5164790_dp
        sep%beta(1) =  -1.7555170_dp
        sep%sto_exponents(0)    =   1.2931800_dp
        sep%sto_exponents(1)    =   2.0205640_dp
        sep%dd    =   1.2556091_dp
        sep%qq    =   1.0498655_dp
        sep%am    =   0.3601617_dp
        sep%ad    =   0.3643722_dp
        sep%aq    =   0.4347337_dp
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
        sep%eisol = -346.68125000_dp
        sep%beta(0) =   -8.91710700_dp
        sep%beta(1) =   -9.94374000_dp
        sep%sto_exponents(0)    =    3.85430190_dp
        sep%sto_exponents(1)    =    2.19920910_dp
        sep%sto_exponents(2)    =    1.00000000_dp
        sep%dd    =    0.60510740_dp
        sep%qq    =    0.96458730_dp
        sep%am    =    0.55260680_dp
        sep%ad    =    0.72583300_dp
        sep%aq    =    0.55745890_dp
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
        sep%eisol = -92.3241020_dp
        sep%beta(0) =  -3.2351470_dp
        sep%beta(1) =  -4.2904160_dp
        sep%sto_exponents(0)    =   2.0803800_dp
        sep%sto_exponents(1)    =   1.9371060_dp
        sep%dd    =   1.5697766_dp
        sep%qq    =   1.3262292_dp
        sep%am    =   0.3601617_dp
        sep%ad    =   0.3219998_dp
        sep%aq    =   0.3713827_dp
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
        sep%eisol = -340.59836000_dp
        sep%gss   =   15.04044855_dp
        sep%gsp   =   13.05655798_dp
        sep%gpp   =   11.14778369_dp
        sep%gp2   =    9.91409071_dp
        sep%hsp   =    2.45638202_dp
        sep%dd    =    1.42532330_dp
        sep%qq    =    1.18417070_dp
        sep%am    =    0.55275410_dp
        sep%ad    =    0.45934510_dp
        sep%aq    =    0.45853760_dp
      CASE (54:79)
      CASE (80)
        !Mercury
        sep%defined = .TRUE.
        sep%alp   =   1.3356410_dp
        sep%eisol = -28.8191480_dp
        sep%beta(0) =  -0.4045250_dp
        sep%beta(1) =  -6.2066830_dp
        sep%sto_exponents(0)    =   2.2181840_dp
        sep%sto_exponents(1)    =   2.0650380_dp
        sep%dd    =   1.7378048_dp
        sep%qq    =   1.4608064_dp
        sep%am    =   0.3969129_dp
        sep%ad    =   0.3047694_dp
        sep%aq    =   0.3483102_dp
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
        sep%eisol = -105.8345040_dp
        sep%beta(0) =   -8.0423870_dp
        sep%beta(1) =   -3.0000000_dp
        sep%sto_exponents(0)    =    2.4982860_dp
        sep%sto_exponents(1)    =    2.0820710_dp
        sep%dd    =    1.5526624_dp
        sep%qq    =    1.4488558_dp
        sep%am    =    0.3601617_dp
        sep%ad    =    0.3239309_dp
        sep%aq    =    0.3502057_dp
        sep%uss   =  -47.3196920_dp
        sep%upp   =  -28.8475600_dp
        sep%gss   =    9.8000000_dp
        sep%gsp   =    8.3000000_dp
        sep%gpp   =    7.3000000_dp
        sep%gp2   =    6.5000000_dp
        sep%hsp   =    1.3000000_dp
      CASE (83:106)
    END SELECT

    sep%beta = sep%beta/evolt
    sep%uss = sep%uss/evolt; sep%upp = sep%upp/evolt; sep%udd = sep%udd/evolt
    sep%alp = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss = sep%gss/evolt; sep%gsp = sep%gsp/evolt; sep%gpp = sep%gpp/evolt; sep%gp2 = sep%gp2/evolt 
    sep%gsd = sep%gsd/evolt; sep%gpd = sep%gpd/evolt; sep%gdd = sep%gdd/evolt
    sep%hsp = sep%hsp/evolt
    sep%fn1 = sep%fn1; sep%fn2 = sep%fn2/bohr/bohr; sep%fn3 = sep%fn3*bohr

  END SUBROUTINE mndo_default_parameter

! *****************************************************************************

  SUBROUTINE am1_default_parameter ( sep, z )
    TYPE(semi_empirical_type),POINTER    :: sep
    INTEGER                              :: z

    CHARACTER(LEN=*), PARAMETER :: routine_name = "am1_default_parameter"

! -----------------------------------------------------------------------------

    sep%typ = "AM1"
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


    SELECT CASE (z)
      CASE DEFAULT
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "invalid atomic number")
      CASE (0)
        sep%defined = .TRUE.
      CASE (1)
        !Hydrogen
        sep%defined          = .TRUE.
        sep%alp   =   2.8823240_dp
        sep%eisol = -11.3964270_dp
        sep%beta(0) =  -6.1737870_dp
        sep%sto_exponents(0)    =   1.1880780_dp
        sep%am    =   0.4721793_dp
        sep%ad    =   0.4721793_dp
        sep%aq   =   0.4721793_dp
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
        sep%eisol =  -5.1280000_dp
        sep%beta(0) =  -1.3500400_dp
        sep%beta(1) =  -1.3500400_dp
        sep%sto_exponents(0)    =   0.7023800_dp
        sep%sto_exponents(1)   =   0.7023800_dp
        sep%dd   =   2.0549783_dp
        sep%qq    =   1.7437069_dp
        sep%am    =   0.2682837_dp
        sep%ad    =   0.2269793_dp
        sep%aq   =   0.2614581_dp
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
        sep%eisol = -24.2047560_dp
        sep%beta(0) =  -4.0170960_dp
        sep%beta(1)=  -4.0170960_dp
        sep%sto_exponents(0)    =   1.0042100_dp
        sep%sto_exponents(1)   =   1.0042100_dp
        sep%dd   =   1.4373245_dp
        sep%qq    =   1.2196103_dp
        sep%am    =   0.3307607_dp
        sep%ad    =   0.3356142_dp
        sep%aq   =   0.3846373_dp
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
        sep%eisol = -63.7172650_dp
        sep%beta(0) =  -9.5991140_dp
        sep%beta(1)=  -6.2737570_dp
        sep%sto_exponents(0)    =   1.6117090_dp
        sep%sto_exponents(1)   =   1.5553850_dp
        sep%dd   =   0.9107622_dp
        sep%qq    =   0.7874223_dp
        sep%am    =   0.3891951_dp
        sep%ad    =   0.5045152_dp
        sep%aq   =   0.5678856_dp
        sep%uss   = -34.4928700_dp
        sep%upp   = -22.6315250_dp
        sep%gss   =  10.5900000_dp
        sep%gsp  =   9.5600000_dp
        sep%gpp  =   8.8600000_dp
        sep%gp2  =   7.8600000_dp
        sep%hsp  =   1.8100000_dp
      CASE (6)
        !Carbon
        sep%defined = .TRUE.
        sep%alp   =    2.6482740_dp
        sep%eisol = -120.8157940_dp
        sep%beta(0) =  -15.7157830_dp
        sep%beta(1) =   -7.7192830_dp
        sep%sto_exponents(0)     =    1.8086650_dp
        sep%sto_exponents(1)     =    1.6851160_dp
        sep%dd    =    0.8236736_dp
        sep%qq    =    0.7268015_dp
        sep%am    =    0.4494671_dp
        sep%ad    =    0.6082946_dp
        sep%aq    =    0.6423492_dp
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
        sep%eisol= -202.4077430_dp
        sep%beta(0) =  -20.2991100_dp
        sep%beta(1) =  -18.2386660_dp
        sep%sto_exponents(0)   =    2.3154100_dp
        sep%sto_exponents(1)    =    2.1579400_dp
        sep%dd   =    0.6433247_dp
        sep%qq    =    0.5675528_dp
        sep%am    =    0.4994487_dp
        sep%ad    =    0.7820840_dp
        sep%aq    =    0.7883498_dp
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
        sep%eisol= -316.0995200_dp
        sep%beta(0) =  -29.2727730_dp
        sep%beta(1) =  -29.2727730_dp
        sep%sto_exponents(0)   =    3.1080320_dp
        sep%sto_exponents(1)    =    2.5240390_dp
        sep%dd    =    0.4988896_dp
        sep%qq    =    0.4852322_dp
        sep%am    =    0.5667034_dp
        sep%ad    =    0.9961066_dp
        sep%aq    =    0.9065223_dp
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
      CASE (10)
      !Fluorine
        sep%defined = .TRUE.
        sep%alp   =    5.5178000_dp
        sep%eisol= -482.2905830_dp
        sep%beta(0) =  -69.5902770_dp
        sep%beta(1) =  -27.9223600_dp
        sep%sto_exponents(0)   =    3.7700820_dp
        sep%sto_exponents(1)    =    2.4946700_dp
        sep%dd    =    0.4145203_dp
        sep%qq    =    0.4909446_dp
        sep%am    =    0.6218302_dp
        sep%ad    =    1.2088792_dp
        sep%aq    =    0.9449355_dp
        sep%uss   = -136.1055790_dp
        sep%upp   = -104.8898850_dp
        sep%gss   =   16.9200000_dp
        sep%gsp   =   17.2500000_dp
        sep%gpp   =   16.7100000_dp
        sep%gp2   =   14.9100000_dp
        sep%hsp   =    4.8300000_dp
        sep%fn1(1) = 0.2420790_dp; sep%fn2(1) = 4.8000000_dp; sep%fn3(1) = 0.9300000_dp
        sep%fn1(2) = 0.0036070_dp; sep%fn2(2) = 4.6000000_dp; sep%fn3(2) = 1.6600000_dp
      CASE (11)
        !Sodium
        sep%defined = .TRUE.
        sep%alp   =  1.32_dp
        sep%eisol=  0.00_dp
        sep%am    =  0.50_dp
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
        sep%eisol =  -22.43786349_dp 
        sep%gss   =    7.50132277_dp
        sep%gsp  =    6.34591536_dp
        sep%gpp  =    4.77534467_dp
        sep%gp2  =    4.34017279_dp
        sep%hsp  =    0.48930466_dp
        sep%dd   =    1.75012115_dp 
        sep%qq    =    1.64001467_dp 
        sep%am    =    0.27568257_dp 
        sep%ad    =    0.19957236_dp 
        sep%aq   =    0.26440152_dp 
        sep%fn1(1) =  2.55017735_dp; sep%fn2(1) = 4.29397225_dp; sep%fn3(1) = 0.79989601_dp
        sep%fn1(2) = -0.00565806_dp; sep%fn2(2) = 2.96053910_dp; sep%fn3(2) = 1.47499983_dp
        sep%fn1(3) = -0.00610286_dp; sep%fn2(3) = 2.61416919_dp; sep%fn3(3) = 2.42604040_dp
      CASE (13)
       !Aluminium
        sep%defined = .TRUE.
        sep%alp   =   1.9765860_dp
        sep%eisol = -46.4208150_dp
        sep%beta(0) =  -3.8668220_dp
        sep%beta(1)=  -2.3171460_dp
        sep%sto_exponents(0)   =   1.5165930_dp
        sep%sto_exponents(1)   =   1.3063470_dp
        sep%sto_exponents(2)   =   1.0000000_dp
        sep%dd   =   1.4040443_dp
        sep%qq    =   1.2809154_dp
        sep%am    =   0.2973172_dp
        sep%ad    =   0.2630229_dp
        sep%aq   =   0.3427832_dp
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
        sep%eisol = -79.0017420_dp
        sep%beta(0) =  -3.784852_dp
        sep%beta(1)=  -1.968123_dp
        sep%sto_exponents(0)   =   1.830697_dp
        sep%sto_exponents(1)   =   1.2849530_dp
        sep%sto_exponents(2)   =   1.0000000_dp
        sep%dd   =   1.1631107_dp
        sep%qq    =   1.3022422_dp
        sep%am    =   0.3608967_dp
        sep%ad    =   0.3829813_dp
        sep%aq   =   0.3999442_dp
        sep%uss   = -33.9536220_dp
        sep%upp   = -28.9347490_dp
        sep%gss   =   9.8200000_dp
        sep%gsp  =   8.3600000_dp
        sep%gpp  =   7.3100000_dp
        sep%gp2  =   6.5400000_dp
        sep%hsp  =   1.3200000_dp
      CASE (15)
      !Phosphorus
        sep%defined = .TRUE.
        sep%alp   =    2.4553220_dp
        sep%eisol = -124.4368355_dp
        sep%beta(0) =   -6.3537640_dp
        sep%beta(1)=   -6.5907090_dp
        sep%sto_exponents(0)   =    1.9812800_dp
        sep%sto_exponents(1)   =    1.8751500_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%dd   =    1.0452022_dp
        sep%qq    =    0.8923660_dp
        sep%am    =    0.4248440_dp
        sep%ad    =    0.3275319_dp
        sep%aq   =    0.4386854_dp
        sep%uss   =  -42.0298630_dp
        sep%upp   =  -34.0307090_dp
        sep%gss   =   11.5600050_dp
        sep%gsp  =   5.2374490_dp
        sep%gpp  =    7.8775890_dp
        sep%gp2  =    7.3076480_dp
        sep%hsp  =    0.7792380_dp
      CASE (16)
        !Sulphur
        sep%defined = .TRUE.
        sep%alp   =    2.4616480_dp
        sep%eisol = -191.7321930_dp
        sep%beta(0) =  -3.9205660_dp
        sep%beta(1)=  -7.9052780_dp
        sep%sto_exponents(0)   =    2.3665150_dp
        sep%sto_exponents(1)   =    1.6672630_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%dd   =    0.9004265_dp
        sep%qq    =    1.0036329_dp
        sep%am    =    0.4331617_dp
        sep%ad    =    0.5907115_dp
        sep%aq   =    0.6454943_dp
        sep%uss   =  -56.6940560_dp
        sep%upp   =  -48.7170490_dp
        sep%gss   =   11.7863290_dp
        sep%gsp  =    8.6631270_dp
        sep%gpp  =   10.0393080_dp
        sep%gp2  =    7.7816880_dp
        sep%hsp  =    2.5321370_dp
      CASE (17)
        !Chlorine
        sep%defined = .TRUE.
        sep%alp   =    2.9193680_dp
        sep%eisol = -372.1984310_dp
        sep%beta(0) =  -24.5946700_dp
        sep%beta(1)=  -14.6372160_dp
        sep%sto_exponents(0)    =    3.6313760_dp
        sep%sto_exponents(1)   =    2.0767990_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%dd   =    0.5406286_dp
        sep%qq    =    0.8057208_dp
        sep%am    =    0.5523705_dp
        sep%ad    =    0.7693200_dp
        sep%aq   =    0.6133369_dp
        sep%uss   = -111.6139480_dp
        sep%upp   =  -76.6401070_dp
        sep%gss   =   15.0300000_dp
        sep%gsp  =   13.1600000_dp
        sep%gpp  =   11.3000000_dp
        sep%gp2  =    9.9700000_dp
        sep%hsp  =    2.4200000_dp
        sep%fn1(1) = 0.0942430_dp; sep%fn2(1) = 4.0000000_dp; sep%fn3(1) = 1.3000000_dp
        sep%fn1(2) = 0.0271680_dp; sep%fn2(2) = 4.0000000_dp; sep%fn3(2) = 2.1000000_dp
      CASE (18)
      CASE (19)
        !Potassium
        sep%defined = .TRUE.
        sep%alp   =  1.16_dp
        sep%eisol =  0.00_dp
        sep%am    =  0.50_dp
      CASE (20:29)
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
        sep%eisol  =    -30.2800160_dp
        sep%gss    =     11.8000000_dp
        sep%gsp   =     11.1820180_dp
        sep%gpp   =     13.3000000_dp
        sep%gp2   =     12.9305200_dp
        sep%hsp   =      0.4846060_dp
        sep%dd    =      1.3581113_dp
        sep%qq     =      1.5457406_dp
        sep%am     =      0.4336641_dp
        sep%ad     =      0.2317423_dp
        sep%aq    =      0.2621165_dp
      CASE (31:34)
      CASE (35)
        !Bromine
        sep%defined = .TRUE.
        sep%alp   =    2.5765460_dp
        sep%eisol = -352.3142087_dp
        sep%beta(0) =  -19.3998800_dp
        sep%beta(1)=   -8.9571950_dp
        sep%sto_exponents(0)    =    3.0641330_dp
        sep%sto_exponents(1)   =    2.0383330_dp
        sep%sto_exponents(2)   =    1.0000000_dp
        sep%dd   =    0.8458104_dp
        sep%qq    =    1.0407133_dp
        sep%am    =    0.5526071_dp
        sep%ad    =    0.6024598_dp
        sep%aq   =    0.5307555_dp
        sep%uss   = -104.6560630_dp
        sep%upp   =  -74.9300520_dp
        sep%gss   =   15.0364395_dp
        sep%gsp  =   13.0346824_dp
        sep%gpp  =   11.2763254_dp
        sep%gp2  =    9.8544255_dp
        sep%hsp  =    2.4558683_dp
        sep%fn1(1) = 0.0666850_dp; sep%fn2(1) = 4.0000000_dp; sep%fn3(1) = 1.5000000_dp
        sep%fn1(2) = 0.0255680_dp; sep%fn2(2) = 4.0000000_dp; sep%fn3(2) = 2.3000000_dp
      CASE (36:52)
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
        sep%eisol = -346.8642857_dp
        sep%gss   =   15.0404486_dp
        sep%gsp  =   13.0565580_dp
        sep%gpp  =   11.1477837_dp
        sep%gp2  =    9.9140907_dp
        sep%hsp  =    2.4563820_dp
        sep%dd   =    1.4878778_dp
        sep%qq    =    1.1887388_dp
        sep%am    =    0.5527544_dp
        sep%ad    =    0.4497523_dp
        sep%aq   =    0.4631775_dp
        sep%fn1(1) = 0.0043610_dp; sep%fn2(1) = 2.3000000_dp; sep%fn3(1) = 1.8000000_dp
        sep%fn1(2) = 0.0157060_dp; sep%fn2(2) = 3.0000000_dp; sep%fn3(2) = 2.2400000_dp
      CASE (54:106)
    END SELECT

    sep%beta = sep%beta/evolt
    sep%uss = sep%uss/evolt; sep%upp = sep%upp/evolt; sep%udd = sep%udd/evolt
    sep%alp = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss = sep%gss/evolt; sep%gsp = sep%gsp/evolt; sep%gpp = sep%gpp/evolt; sep%gp2 = sep%gp2/evolt 
    sep%gsd = sep%gsd/evolt; sep%gpd = sep%gpd/evolt; sep%gdd = sep%gdd/evolt
    sep%hsp = sep%hsp/evolt
    sep%fn1 = sep%fn1; sep%fn2 = sep%fn2/bohr/bohr; sep%fn3 = sep%fn3*bohr

  END SUBROUTINE am1_default_parameter

! *****************************************************************************

  SUBROUTINE pm3_default_parameter ( sep, z )
    TYPE(semi_empirical_type),POINTER    :: sep
    INTEGER                              :: z

    CHARACTER(LEN=*), PARAMETER :: routine_name = "pm3_default_parameter"

! -----------------------------------------------------------------------------

    sep%typ = "PM3"
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


    SELECT CASE (z)
      CASE DEFAULT
        CALL stop_program(routine_name,module_name,__LINE__,&
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
        sep%eisol =  -13.0733210_dp
        sep%gss   =   14.7942080_dp
        sep%am    =    0.5437048_dp
        sep%ad    =    0.5437048_dp
        sep%aq   =    0.5437048_dp
        sep%fn1(1) =  1.1287500_dp; sep%fn2(1) = 5.0962820_dp; sep%fn3(1) = 1.5374650_dp
        sep%fn1(2) = -1.0603290_dp; sep%fn2(2) = 6.0037880_dp; sep%fn3(2) = 1.5701890_dp
      CASE (2:3)
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
        sep%eisol =  -25.5166530_dp
        sep%gss   =    9.0128510_dp
        sep%gsp  =    6.5761990_dp
        sep%gpp  =    6.0571820_dp
        sep%gp2  =    9.0052190_dp
        sep%hsp  =    0.5446790_dp
        sep%dd   =    1.0090531_dp
        sep%qq    =    0.8117586_dp
        sep%am    =    0.3312330_dp
        sep%ad    =    0.2908996_dp
        sep%aq   =    0.3530008_dp
        sep%fn1(1) =  1.6315720_dp; sep%fn2(1) = 2.6729620_dp; sep%fn3(1) = 1.7916860_dp
        sep%fn1(2) = -2.1109590_dp; sep%fn2(2) = 1.9685940_dp; sep%fn3(2) = 1.7558710_dp
      CASE (5)
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
        sep%eisol = -111.2299170_dp
        sep%gss   =   11.2007080_dp
        sep%gsp  =   10.2650270_dp
        sep%gpp  =   10.7962920_dp
        sep%gp2  =    9.0425660_dp
        sep%hsp  =    2.2909800_dp
        sep%dd   =    0.8332396_dp
        sep%qq    =    0.6647750_dp
        sep%am    =    0.4116394_dp
        sep%ad    =    0.5885862_dp
        sep%aq   =    0.7647667_dp
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
        sep%eisol = -157.6137755_dp
        sep%gss   =   11.9047870_dp
        sep%gsp  =    7.3485650_dp
        sep%gpp  =   11.7546720_dp
        sep%gp2  =   10.8072770_dp
        sep%hsp  =    1.1367130_dp
        sep%dd   =    0.6577006_dp
        sep%qq    =    0.5293383_dp
        sep%am    =    0.4375151_dp
        sep%ad    =    0.5030995_dp
        sep%aq   =    0.7364933_dp
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
        sep%eisol = -289.3422065_dp
        sep%gss   =   15.7557600_dp
        sep%gsp  =   10.6211600_dp
        sep%gpp  =   13.6540160_dp
        sep%gp2  =   12.4060950_dp
        sep%hsp  =    0.5938830_dp
        sep%dd   =    0.4086173_dp
        sep%qq    =    0.5125738_dp
        sep%am    =    0.5790430_dp
        sep%ad    =    0.5299630_dp
        sep%aq   =    0.8179630_dp
        sep%fn1(1) = -1.1311280_dp; sep%fn2(1) = 6.0024770_dp; sep%fn3(1) = 1.6073110_dp
        sep%fn1(2) =  1.1378910_dp; sep%fn2(2) = 5.9505120_dp; sep%fn3(2) = 1.5983950_dp
      CASE (9)
      CASE (10)
        !Fluorine
        sep%defined = .TRUE.
        sep%uss   = -110.4353030_dp
        sep%upp   = -105.6850470_dp
        sep%beta(0) =  -48.4059390_dp
        sep%beta(1)=  -27.7446600_dp
        sep%sto_exponents(0)    =    4.7085550_dp
        sep%sto_exponents(1)   =    2.4911780_dp
        sep%alp   =    3.3589210_dp
        sep%eisol = -437.5171690_dp
        sep%gss   =   10.4966670_dp
        sep%gsp  =   16.0736890_dp
        sep%gpp  =   14.8172560_dp
        sep%gp2  =   14.4183930_dp
        sep%hsp  =    0.7277630_dp
        sep%dd   =    0.3125302_dp
        sep%qq    =    0.4916328_dp
        sep%am    =    0.3857650_dp
        sep%ad    =    0.6768503_dp
        sep%aq   =    0.6120047_dp
        sep%fn1(1) = -0.0121660_dp; sep%fn2(1) = 6.0235740_dp; sep%fn3(1) = 1.8568590_dp
        sep%fn1(2) = -0.0028520_dp; sep%fn2(2) = 6.0037170_dp; sep%fn3(2) = 2.6361580_dp
      CASE (11)
       !Sodium
        sep%defined = .TRUE.
        sep%alp   = 1.681_dp
        sep%eisol = 0.000_dp
        sep%am    = 0.500_dp
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
        sep%eisol =  -22.5530760_dp
        sep%gss   =    6.6943000_dp
        sep%gsp  =    6.7939950_dp
        sep%gpp  =    6.9104460_dp
        sep%gp2  =    7.0908230_dp
        sep%hsp  =    0.5433000_dp
        sep%dd   =    1.1403950_dp
        sep%qq    =    1.1279899_dp
        sep%am    =    0.2460235_dp
        sep%ad    =    0.2695751_dp
        sep%aq   =    0.2767522_dp
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
        sep%eisol =  -46.8647630_dp
        sep%gss   =    5.7767370_dp
        sep%gsp  =   11.6598560_dp
        sep%gpp  =    6.3477900_dp
        sep%gp2  =    6.1210770_dp
        sep%hsp  =    4.0062450_dp
        sep%dd   =    1.2102799_dp
        sep%qq    =    1.5585645_dp
        sep%am    =    0.2123020_dp
        sep%ad    =    0.6418584_dp
        sep%aq   =    0.2262838_dp
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
        sep%eisol =  -67.7882140_dp
        sep%gss   =    5.0471960_dp
        sep%gsp  =    5.9490570_dp
        sep%gpp  =    6.7593670_dp
        sep%gp2  =    5.1612970_dp
        sep%hsp  =    0.9198320_dp
        sep%dd   =    1.3144550_dp
        sep%qq    =    1.2743396_dp
        sep%am    =    0.1854905_dp
        sep%ad    =    0.3060715_dp
        sep%aq   =    0.4877432_dp
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
        sep%eisol = -117.9591740_dp
        sep%gss   =    7.8016150_dp
        sep%gsp  =    5.1869490_dp
        sep%gpp  =    6.6184780_dp
        sep%gp2  =    6.0620020_dp
        sep%hsp  =    1.5428090_dp
        sep%dd   =    1.0644947_dp
        sep%qq    =    1.1120386_dp
        sep%am    =    0.2867187_dp
        sep%ad    =    0.4309446_dp
        sep%aq   =    0.3732517_dp
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
        sep%eisol = -183.4537395_dp
        sep%gss   =    8.9646670_dp
        sep%gsp  =    6.7859360_dp
        sep%gpp  =    9.9681640_dp
        sep%gp2  =    7.9702470_dp
        sep%hsp  =    4.0418360_dp
        sep%dd   =    1.1214313_dp
        sep%qq    =    1.0086488_dp
        sep%am    =    0.3294622_dp
        sep%ad    =    0.6679118_dp
        sep%aq   =    0.6137472_dp
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
        sep%eisol = -315.1949480_dp
        sep%gss   =   16.0136010_dp
        sep%gsp  =    8.0481150_dp
        sep%gpp  =    7.5222150_dp
        sep%gp2  =    7.5041540_dp
        sep%hsp  =    3.4811530_dp
        sep%dd   =    0.9175856_dp
        sep%qq    =    0.7779230_dp
        sep%am    =    0.5885190_dp
        sep%ad    =    0.6814522_dp
        sep%aq   =    0.3643694_dp
        sep%fn1(1) = -0.1715910_dp; sep%fn2(1) = 6.0008020_dp; sep%fn3(1) = 1.0875020_dp
        sep%fn1(2) = -0.0134580_dp; sep%fn2(2) = 1.9666180_dp; sep%fn3(2) = 2.2928910_dp
      CASE (18)
      CASE (19)
        !Potassium
        sep%defined = .TRUE.
        sep%alp   = 1.400_dp
        sep%eisol = 0.000_dp
        sep%am    = 0.500_dp
      CASE (20:29)
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
        sep%eisol =  -27.3872000_dp
        sep%gss   =    9.6771960_dp
        sep%gsp  =    7.7362040_dp
        sep%gpp  =    4.9801740_dp
        sep%gp2  =    4.6696560_dp
        sep%hsp  =    0.6004130_dp
        sep%dd   =    1.5005758_dp
        sep%qq    =    1.4077174_dp
        sep%am    =    0.3556485_dp
        sep%ad    =    0.2375689_dp
        sep%aq   =    0.2661069_dp
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
        sep%eisol =  -57.3280250_dp
        sep%gss   =    8.4585540_dp
        sep%gsp  =    8.9256190_dp
        sep%gpp  =    5.0868550_dp
        sep%gp2  =    4.9830450_dp
        sep%hsp  =    2.0512600_dp
        sep%dd   =    0.9776692_dp
        sep%qq    =    2.5271534_dp
        sep%am    =    0.3108620_dp
        sep%ad    =    0.5129360_dp
        sep%aq   =    0.1546208_dp
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
        sep%eisol =  -84.0156006_dp
        sep%gss   =    5.3769635_dp
        sep%gsp  =   10.2095293_dp
        sep%gpp  =    7.6718647_dp
        sep%gp2  =    6.9242663_dp
        sep%hsp  =    1.3370204_dp
        sep%dd   =    1.1920304_dp
        sep%qq    =    1.3321263_dp
        sep%am    =    0.1976098_dp
        sep%ad    =    0.3798182_dp
        sep%aq   =    0.3620669_dp
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
        sep%eisol = -122.6326140_dp
        sep%gss   =    8.7890010_dp
        sep%gsp  =    5.3979830_dp
        sep%gpp  =    8.2872500_dp
        sep%gp2  =    8.2103460_dp
        sep%hsp  =    1.9510340_dp
        sep%dd   =    0.9679655_dp
        sep%qq    =    1.2449874_dp
        sep%am    =    0.3230063_dp
        sep%ad    =    0.5042239_dp
        sep%aq   =    0.2574219_dp
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
        sep%eisol = -192.7748115_dp
        sep%gss   =    7.4325910_dp
        sep%gsp  =   10.0604610_dp
        sep%gpp  =    9.5683260_dp
        sep%gp2  =    7.7242890_dp
        sep%hsp  =    4.0165580_dp
        sep%dd   =    0.8719813_dp
        sep%qq    =    1.2244019_dp
        sep%am    =    0.2731566_dp
        sep%ad    =    0.7509697_dp
        sep%aq   =    0.5283737_dp
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
        sep%eisol = -352.5398970_dp
        sep%gss   =   15.9434250_dp
        sep%gsp  =   16.0616800_dp
        sep%gpp  =    8.2827630_dp
        sep%gp2  =    7.8168490_dp
        sep%hsp  =    0.5788690_dp
        sep%dd   =    0.2759025_dp
        sep%qq    =    0.9970532_dp
        sep%am    =    0.5859399_dp
        sep%ad    =    0.6755383_dp
        sep%aq   =    0.3823719_dp
        sep%fn1(1) =  0.9604580_dp; sep%fn2(1) = 5.9765080_dp; sep%fn3(1) = 2.3216540_dp
        sep%fn1(2) = -0.9549160_dp; sep%fn2(2) = 5.9447030_dp; sep%fn3(2) = 2.3281420_dp
      CASE (36:47)
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
        sep%eisol =  -22.4502080_dp
        sep%gss   =    9.2069600_dp
        sep%gsp  =    8.2315390_dp
        sep%gpp  =    4.9481040_dp
        sep%gp2  =    4.6696560_dp
        sep%hsp  =    1.6562340_dp
        sep%dd   =    1.5982681_dp
        sep%qq    =    1.2432402_dp
        sep%am    =    0.3383668_dp
        sep%ad    =    0.3570290_dp
        sep%aq   =    0.2820582_dp
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
        sep%eisol =  -51.9750470_dp
        sep%gss   =    6.5549000_dp
        sep%gsp  =    8.2298730_dp
        sep%gpp  =    6.2992690_dp
        sep%gp2  =    4.9842110_dp
        sep%hsp  =    2.6314610_dp
        sep%dd   =    1.5766241_dp
        sep%qq    =    1.7774563_dp
        sep%am    =    0.2409004_dp
        sep%ad    =    0.4532655_dp
        sep%aq   =    0.3689812_dp
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
        sep%eisol =  -78.8877790_dp
        sep%gss   =   10.1900330_dp
        sep%gsp  =    7.2353270_dp
        sep%gpp  =    5.6738100_dp
        sep%gp2  =    5.1822140_dp
        sep%hsp  =    1.0331570_dp
        sep%dd   =    1.3120038_dp
        sep%qq    =    1.5681814_dp
        sep%am    =    0.3744959_dp
        sep%ad    =    0.3218163_dp
        sep%aq   =    0.2832529_dp
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
        sep%eisol = -148.9382890_dp
        sep%gss   =    9.2382770_dp
        sep%gsp  =    5.2776800_dp
        sep%gpp  =    6.3500000_dp
        sep%gp2  =    6.2500000_dp
        sep%hsp  =    2.4244640_dp
        sep%dd   =    1.4091903_dp
        sep%qq    =    1.3521354_dp
        sep%am    =    0.3395177_dp
        sep%ad    =    0.4589010_dp
        sep%aq   =    0.2423472_dp
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
        sep%eisol = -168.0945925_dp
        sep%gss   =   10.2550730_dp
        sep%gsp  =    8.1691450_dp
        sep%gpp  =    7.7775920_dp
        sep%gp2  =    7.7551210_dp
        sep%hsp  =    3.7724620_dp
        sep%dd   =    0.3484177_dp
        sep%qq    =    1.5593085_dp
        sep%am    =    0.3768862_dp
        sep%ad    =    1.1960743_dp
        sep%aq   =    0.2184786_dp
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
        sep%eisol = -288.3160860_dp
        sep%gss   =   13.6319430_dp
        sep%gsp  =   14.9904060_dp
        sep%gpp  =    7.2883300_dp
        sep%gp2  =    5.9664070_dp
        sep%hsp  =    2.6300350_dp
        sep%dd   =    0.1581469_dp
        sep%qq    =    1.0467302_dp
        sep%am    =    0.5009902_dp
        sep%ad    =    1.6699104_dp
        sep%aq   =    0.5153082_dp
        sep%fn1(1) = -0.1314810_dp; sep%fn2(1) = 5.2064170_dp; sep%fn3(1) = 1.7488240_dp
        sep%fn1(2) = -0.0368970_dp; sep%fn2(2) = 6.0101170_dp; sep%fn3(2) = 2.7103730_dp
      CASE (54:79)
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
        sep%eisol =  -28.8997380_dp
        sep%gss   =    6.6247200_dp
        sep%gsp  =   10.6392970_dp
        sep%gpp  =   14.7092830_dp
        sep%gp2  =   16.0007400_dp
        sep%hsp  =    2.0363110_dp
        sep%dd   =    1.2317811_dp
        sep%qq    =    1.2164033_dp
        sep%am    =    0.2434664_dp
        sep%ad    =    0.4515472_dp
        sep%aq   =    0.2618394_dp
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
        sep%eisol =  -56.6492050_dp
        sep%gss   =   10.4604120_dp
        sep%gsp  =   11.2238830_dp
        sep%gpp  =    4.9927850_dp
        sep%gp2  =    8.9627270_dp
        sep%hsp  =    2.5304060_dp
        sep%dd   =    0.0781362_dp
        sep%qq    =    1.5317110_dp
        sep%am    =    0.3844326_dp
        sep%ad    =    2.5741815_dp
        sep%aq   =    0.2213264_dp
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
        sep%eisol =  -73.4660775_dp
        sep%gss   =    7.0119920_dp
        sep%gsp  =    6.7937820_dp
        sep%gpp  =    5.1837800_dp
        sep%gp2  =    5.0456510_dp
        sep%hsp  =    1.5663020_dp
        sep%dd   =    0.9866290_dp
        sep%qq    =    1.5940562_dp
        sep%am    =    0.2576991_dp
        sep%ad    =    0.4527678_dp
        sep%aq   =    0.2150175_dp
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
        sep%eisol = -109.2774910_dp
        sep%gss   =    4.9894800_dp
        sep%gsp  =    6.1033080_dp
        sep%gpp  =    8.6960070_dp
        sep%gp2  =    8.3354470_dp
        sep%hsp  =    0.5991220_dp
        sep%dd   =    0.2798609_dp
        sep%qq    =    1.5590294_dp
        sep%am    =    0.1833693_dp
        sep%ad    =    0.6776013_dp
        sep%aq   =    0.2586520_dp
        sep%fn1(1) = 2.5816930_dp; sep%fn2(1) = 5.0940220_dp; sep%fn3(1) = 0.4997870_dp
        sep%fn1(2) = 0.0603200_dp; sep%fn2(2) = 6.0015380_dp; sep%fn3(2) = 2.4279700_dp
      CASE (84:106)
    END SELECT

    sep%beta = sep%beta/evolt
    sep%uss = sep%uss/evolt; sep%upp = sep%upp/evolt; sep%udd = sep%udd/evolt
    sep%alp = sep%alp/bohr
    sep%eisol = sep%eisol/evolt
    sep%gss = sep%gss/evolt; sep%gsp = sep%gsp/evolt; sep%gpp = sep%gpp/evolt; sep%gp2 = sep%gp2/evolt 
    sep%gsd = sep%gsd/evolt; sep%gpd = sep%gpd/evolt; sep%gdd = sep%gdd/evolt
    sep%hsp = sep%hsp/evolt
    sep%fn1 = sep%fn1; sep%fn2 = sep%fn2/bohr/bohr; sep%fn3 = sep%fn3*bohr

  END SUBROUTINE pm3_default_parameter

! *****************************************************************************

  
  SUBROUTINE valence_electrons ( z, zeff, natorb )
    INTEGER, INTENT ( IN )              :: z
    INTEGER, INTENT ( OUT )             :: zeff
    INTEGER, INTENT ( OUT )             :: natorb

    INTEGER, DIMENSION(0:nelem), PARAMETER  :: zval = (/ -1, &
       1,  0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,  0,&
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,  2,  3,  4,  5,  6,  7,  0,&
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,  2,  3,  4,  5,  6,  7,  0,&
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,  3,  4,&
       5,  6,  7,  8,  9, 10, 11,  2,  3,  4,  5,  6,  7, -1, -1, -1, -1,  4,&
      -1,  6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1         &
       /)
    INTEGER, DIMENSION(0:nelem), PARAMETER  :: norb = (/ -1, &
       1,  1,  4,  4,  4,  4,  4,  4,  4,  4,  0,  4,  4,  4,  4,  4,  4,  4,&
       0,  4,  9,  9,  9,  9,  9,  9,  9,  9,  9,  4,  4,  4,  4,  4,  4,  4,&
       4,  4,  9,  9,  9,  9,  9,  9,  9,  9,  9,  4,  4,  4,  4,  4,  4,  4,&
       2,  2,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  9,  9,&
       9,  9,  9,  9,  9,  9,  9,  4,  4,  4,  4,  4,  4, -1, -1, -1, -1,  9,&
      -1,  9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1         &
       /)

    zeff = zval(z)
    natorb = norb(z)

  END SUBROUTINE valence_electrons 

! *****************************************************************************

  SUBROUTINE slater_exponents (z,norb,sto,itype) 
    INTEGER, INTENT (IN)                 :: z
    INTEGER, DIMENSION(0:3), INTENT (IN) :: norb
    REAL(dp), DIMENSION(:,0:)            :: sto
    INTEGER, INTENT (IN)                 :: itype
!
!   Slater exponents according to
!   itype = 0 : Clementi-Raimondi if available, Slater rules for others
!   itype = 1 : Slater rules
!   itype = 2 : Clementi-Raimondi
!
!   Clementi-Raimondi Table
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab1s = (/ &
     1.0000, 1.6875, 2.6906, 3.6848, 4.6795, 5.6727, 6.6651, 7.6579, 8.6501, 9.6421,&
    10.6259,11.6089,12.5910,13.5745,14.5578,15.5409,16.5239,17.5075 /)
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab2s = (/ &
     0.0000, 0.0000, 0.6396, 0.9560, 1.2881, 1.6083, 1.9237, 2.2458, 2.5638, 2.8792,&
     3.2857, 3.6960, 4.1068, 5.5100, 4.9125, 5.3144, 5.7152, 6.1152 /)
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab3s = (/ &
     0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,&
     0.8358, 1.1025, 1.3724, 1.6344, 1.8806, 2.1223, 2.3561, 2.5856 /)
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab2p = (/ &
     0.0000, 0.0000, 0.0000, 0.0000, 1.2107, 1.5679, 1.9170, 2.2266, 2.5500, 2.8792,&
     3.4009, 3.9129, 4.4817, 4.9725, 5.4806, 5.9885, 6.4966, 7.0041 /)
    REAL(dp),DIMENSION(18),PARAMETER   :: crtab3p = (/ &
     0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,&
     0.0000, 0.0000, 1.3552, 1.4284, 1.6288, 1.8273, 2.0387, 2.2547 /)

    INTEGER, DIMENSION(1:4,1:7)   :: ne
    INTEGER, DIMENSION(1:4)       :: no
    INTEGER                       :: i,l,nj,j

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

  FUNCTION srules(z,ne,n,l)
    ! Slater rules
    REAL(dp)                             :: srules
    INTEGER                              :: z
    INTEGER, DIMENSION(4,7)              :: ne
    INTEGER                              :: n,l

    INTEGER                              :: l1,l2,m,m1,m2,i
    REAL(dp)                             :: s
    REAL(dp), DIMENSION(7), PARAMETER    :: xns = (/ &
              1.0_dp,2.0_dp,3.0_dp,3.7_dp,4.0_dp,4.2_dp,4.4_dp /)

    ! Calculate the shielding
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
  
END MODULE semi_empirical_parameters
