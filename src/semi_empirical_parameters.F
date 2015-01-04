!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Default parameter sets for semi empirical models: sep%...
!>      ass, asp, app, a.u.      parameters for the SCP-NDDO calculations
!>      ad, am, aq,    a.u.      parameters for the multipole approximation to
!>      dd, qq         a.u.      the two-centre, two-electron integrals.
!>      alp            A**(-1)   the exponents used to calculate the core-core
!>                               repulsion terms.
!>      beta           eV        parameters for the two-centre, one-electron core
!>                               resonance integrals.
!>      z                        the CORE charge on each element.
!>      eheat          kcal/mol  the experimental heat of formation for each element.
!>      eisol          eV        the calculated electronic energies for each element.
!>      fn1            A*eV      AM1/PM3 specific parameters
!>      fn2            A**(-2)   for the core-core repulsion
!>      fn3            A         interaction.
!>      pre            eV        PDDG function
!>      d              A         PDDG function
!>      bfn1-3         special core-core interactions
!>
!>      gdd, gpd, gsd, eV        parameters for the Coulomb and exchange one-centre,
!>      gp2, gpp, gsp, eV        two-electron integrals.
!>      gss, hsp       eV
!>      uss, upp, udd  eV        electron kinetic energy integrals.
!>      sto_exponents  a.u.      the Slater exponents of the basis functions.
!>
!>      all dimensions will be converted to a.u. at the end of the routines
!> \par History
!>      BK   (12.10.2004)
!>      JGH  (36.03.2006) added PDDG function
!>      Samuele Giani - Teo (28.01.2008) : Parameters for PM6
!>      Teodoro Laino [tlaino] - 03.2008 - University of Zurich: preparing for d-orbitals
!> \author JGH (14.08.2004)
! *****************************************************************************
MODULE semi_empirical_parameters
  
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_mndod,&
                                             do_method_pchg,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_pnnl,&
                                             do_method_rm1
  USE kinds,                           ONLY: dp
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: angstrom,&
                                             evolt,&
                                             kcalmol
  USE semi_empirical_par_utils,        ONLY: calpar,&
                                             convert_param_to_cp2k,&
                                             valence_electrons
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE termination,                     ONLY: stop_program
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_parameters'

  INTEGER, PARAMETER, PRIVATE :: nelem = 106

  REAL(KIND=dp), DIMENSION(0:nelem,0:nelem) :: xab_pm6, aab_pm6

  PUBLIC :: mndo_default_parameter,&
            am1_default_parameter, &
            rm1_default_parameter,&
            pdg_default_parameter,&
            pm3_default_parameter,&
            pm6_default_parameter,&
            pcharge_default_parameter,&
            pnnl_default_parameter

CONTAINS

! *****************************************************************************
!> \brief Default parameter sets for semi empirical models: POINT_CHARGE
!> \param sep ...
!> \param z ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pcharge_default_parameter(sep, z, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pcharge_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ     = do_method_pchg
    sep%name    = "CLASSICAL"
    sep%eheat   = 0.0_dp
    sep%z       = z
    sep%defined = .TRUE.
    CALL valence_electrons (sep, extended_basis_set=.FALSE., error=error)

  END SUBROUTINE pcharge_default_parameter

! *****************************************************************************
!> \brief Default parameter sets for semi empirical models: MNDO
!> \param sep ...
!> \param z ...
!> \param itype ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mndo_default_parameter ( sep, z, itype, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER, INTENT(IN)                      :: z, itype
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mndo_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ    = itype
    sep%name   = ptable(z)%name
    sep%eheat  = ptable(z)%heat_of_formation / kcalmol
    sep%z      = z

    SELECT CASE (z)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            "invalid atomic number")
    CASE (0)
       sep%defined          = .TRUE.
    CASE (1)
       !Hydrogen
       sep%defined          = .TRUE.
       sep%alp              =   2.5441341_dp
       sep%beta(0)          =  -6.9890640_dp
       sep%sto_exponents(0) =   1.3319670_dp
       sep%uss              = -11.9062760_dp
       sep%gss              =  12.8480000_dp
       sep%ass              =  1.0e-16_dp
       sep%asp              =  1.0e-16_dp
       sep%app              =  1.0e-16_dp
    CASE (2)
    CASE (3)
       !Lithium
       sep%defined          = .TRUE.
       sep%alp              =   1.2083987_dp
       sep%beta(0)          =  -0.1904126_dp
       sep%beta(1)          =  -0.1904126_dp
       sep%sto_exponents(0) =   0.4296141_dp
       sep%sto_exponents(1) =   0.7554884_dp
       sep%uss              =  -4.8578570_dp
       sep%upp              =  -2.0266084_dp
       sep%gss              =   7.5947069_dp
       sep%gsp              =   6.7259856_dp
       sep%gpp              =   8.6596829_dp
       sep%gp2              =   3.8714751_dp
       sep%hsp              =   5.0003381_dp
    CASE (4)
       !Beryllium
       sep%defined          = .TRUE.
       sep%uss              = -16.6023780_dp
       sep%upp              = -10.7037710_dp
       sep%beta(0)          =  -4.0170960_dp
       sep%beta(1)          =  -4.0170960_dp
       sep%sto_exponents(0) =   1.0042100_dp
       sep%sto_exponents(1) =   1.0042100_dp
       sep%alp              =   1.6694340_dp
       sep%gss              =   9.0000000_dp
       sep%gsp              =   7.4300000_dp
       sep%gpp              =   6.9700000_dp
       sep%gp2              =   6.2200000_dp
       sep%hsp              =   1.2800000_dp
    CASE (5)
       !Boron
       sep%defined          = .TRUE.
       sep%uss              = -34.5471300_dp
       sep%upp              = -23.1216900_dp
       sep%beta(0)          =  -8.2520540_dp
       sep%beta(1)          =  -8.2520540_dp
       sep%sto_exponents(0) =   1.5068010_dp
       sep%sto_exponents(1) =   1.5068010_dp
       sep%alp              =   2.1349930_dp
       sep%gss              =  10.5900000_dp
       sep%gsp              =   9.5600000_dp
       sep%gpp              =   8.8600000_dp
       sep%gp2              =   7.8600000_dp
       sep%hsp              =   1.8100000_dp
    CASE (6)
       !Carbon
       sep%defined          = .TRUE.
       sep%uss              = -52.2797450_dp
       sep%upp              = -39.2055580_dp
       sep%beta(0)          = -18.9850440_dp
       sep%beta(1)          =  -7.9341220_dp
       sep%sto_exponents(0) =   1.7875370_dp
       sep%sto_exponents(1) =   1.7875370_dp
       sep%alp              =   2.5463800_dp
       sep%gss              =  12.2300000_dp
       sep%gsp              =  11.4700000_dp
       sep%gpp              =  11.0800000_dp
       sep%gp2              =   9.8400000_dp
       sep%hsp              =   2.4300000_dp
    CASE (7)
       !Nitrogen
       sep%defined          = .TRUE.
       sep%uss              = -71.9321220_dp
       sep%upp              = -57.1723190_dp
       sep%beta(0)          = -20.4957580_dp
       sep%beta(1)          = -20.4957580_dp
       sep%sto_exponents(0) =   2.2556140_dp
       sep%sto_exponents(1) =   2.2556140_dp
       sep%alp              =   2.8613420_dp
       sep%gss              =  13.5900000_dp
       sep%gsp              =  12.6600000_dp
       sep%gpp              =  12.9800000_dp
       sep%gp2              =  11.5900000_dp
       sep%hsp              =   3.1400000_dp
    CASE (8)
       !Oxygen
       sep%defined          = .TRUE.
       sep%uss              = -99.6443090_dp
       sep%upp              = -77.7974720_dp
       sep%beta(0)          = -32.6880820_dp
       sep%beta(1)          = -32.6880820_dp
       sep%sto_exponents(0) =   2.6999050_dp
       sep%sto_exponents(1) =   2.6999050_dp
       sep%alp              =   3.1606040_dp
       sep%gss              =  15.4200000_dp
       sep%gsp              =  14.4800000_dp
       sep%gpp              =  14.5200000_dp
       sep%gp2              =  12.9800000_dp
       sep%hsp              =   3.9400000_dp
       sep%ass              =  1.0e-16_dp
       sep%asp              =  6.0_dp
       sep%app              =  1.0e-16_dp
    CASE (9)
       !Fluorine
       sep%defined          = .TRUE.
       sep%uss              =-131.0715480_dp
       sep%upp              =-105.7821370_dp
       sep%beta(0)          = -48.2904660_dp
       sep%beta(1)          = -36.5085400_dp
       sep%sto_exponents(0) =   2.8484870_dp
       sep%sto_exponents(1) =   2.8484870_dp
       sep%alp              =   3.4196606_dp
       sep%gss              =  16.9200000_dp
       sep%gsp              =  17.2500000_dp
       sep%gpp              =  16.7100000_dp
       sep%gp2              =  14.9100000_dp
       sep%hsp              =   4.8300000_dp
    CASE (10)
    CASE (11)
       !Sodium
       sep%defined          = .TRUE.
       sep%uss              = -5.1235942_dp
       sep%upp              = -3.0124713_dp
       sep%beta(0)          = -1.4916657_dp
       sep%beta(1)          = -0.2208233_dp
       sep%sto_exponents(0) =  0.8213124_dp
       sep%sto_exponents(1) =  1.0303270_dp
       sep%alp              =  5.9940638_dp
       sep%gss              =  6.9934520_dp
       sep%gsp              =  5.4380530_dp
       sep%gpp              =  6.9285408_dp
       sep%gp2              =  2.4299952_dp
       sep%hsp              =  3.1374231_dp
    CASE (12)
       !Magnesium
       sep%defined          = .TRUE.
       sep%uss              = -14.8057921_dp
       sep%upp              = -12.7545494_dp
       sep%beta(0)          =  -0.0999626_dp
       sep%beta(1)          =  -7.6952692_dp
       sep%sto_exponents(0) =   0.9394811_dp
       sep%sto_exponents(1) =   1.3103428_dp
       sep%alp              =   0.5354201_dp
       sep%gss              =   6.9632774_dp
       sep%gsp              =   8.2410934_dp
       sep%gpp              =  10.0000349_dp
       sep%gp2              =   9.0593304_dp
       sep%hsp              =   0.7165063_dp
    CASE (13)
       !Aluminium
       sep%defined          = .TRUE.
       sep%uss              = -23.8070970_dp
       sep%upp              = -17.5198780_dp
       sep%beta(0)          =  -2.6702840_dp
       sep%beta(1)          =  -2.6702840_dp
       sep%sto_exponents(0) =   1.4441610_dp
       sep%sto_exponents(1) =   1.4441610_dp
       sep%alp              =   1.8688394_dp
       sep%gss              =   8.0900000_dp
       sep%gsp              =   6.6300000_dp
       sep%gpp              =   5.9800000_dp
       sep%gp2              =   5.4000000_dp
       sep%hsp              =   0.7000000_dp
    CASE (14)
       !Silicon
       sep%defined          = .TRUE.
       sep%uss              = -37.0375330_dp
       sep%upp              = -27.7696780_dp
       sep%beta(0)          =  -9.0868040_dp
       sep%beta(1)          =  -1.0758270_dp
       sep%sto_exponents(0) =   1.3159860_dp
       sep%sto_exponents(1) =   1.7099430_dp
       sep%alp              =   2.2053160_dp
       sep%gss              =   9.8200000_dp
       sep%gsp              =   8.3600000_dp
       sep%gpp              =   7.3100000_dp
       sep%gp2              =   6.5400000_dp
       sep%hsp              =   1.3200000_dp
    CASE (15)
       !Phosphorus
       sep%defined          = .TRUE.
       sep%uss              = -56.1433600_dp
       sep%upp              = -42.8510800_dp
       sep%beta(0)          =  -6.7916000_dp
       sep%beta(1)          =  -6.7916000_dp
       sep%sto_exponents(0) =   2.1087200_dp
       sep%sto_exponents(1) =   1.7858100_dp
       sep%alp              =   2.4152800_dp
       sep%gss              =  11.5600000_dp
       sep%gsp              =  10.0800000_dp
       sep%gpp              =   8.6400000_dp
       sep%gp2              =   7.6800000_dp
       sep%hsp              =   1.9200000_dp
    CASE (16)
       !Sulphur
       sep%defined          = .TRUE.
       sep%uss              = -72.2422810_dp
       sep%upp              = -56.9732070_dp
       sep%beta(0)          = -10.7616700_dp
       sep%beta(1)          = -10.1084330_dp
       sep%sto_exponents(0) =   2.3129620_dp
       sep%sto_exponents(1) =   2.0091460_dp
       sep%alp              =   2.4780260_dp
       sep%gss              =  12.8800000_dp
       sep%gsp              =  11.2600000_dp
       sep%gpp              =   9.9000000_dp
       sep%gp2              =   8.8300000_dp
       sep%hsp              =   2.2600000_dp
    CASE (17)
       !Chlorine
       sep%defined          = .TRUE.
       sep%uss              =-100.2271660_dp
       sep%upp              = -77.3786670_dp
       sep%beta(0)          = -14.2623200_dp
       sep%beta(1)          = -14.2623200_dp
       sep%sto_exponents(0) =   3.7846450_dp
       sep%sto_exponents(1) =   2.0362630_dp
       sep%alp              =   2.5422010_dp
       sep%gss              =  15.0300000_dp
       sep%gsp              =  13.1600000_dp
       sep%gpp              =  11.3000000_dp
       sep%gp2              =   9.9700000_dp
       sep%hsp              =   2.4200000_dp
    CASE (18)
    CASE (19)
       !Potassium
       sep%defined          = .TRUE.
       sep%uss              =  -3.6401731_dp
       sep%upp              =  -2.0249253_dp
       sep%beta(0)          =  -0.1361851_dp
       sep%beta(1)          =  -2.8142350_dp
       sep%sto_exponents(0) =   0.7276039_dp
       sep%sto_exponents(1) =   0.9871174_dp
       sep%alp              =   0.5616422_dp
       sep%gss              =   3.7939792_dp
       sep%gsp              =   6.4170233_dp
       sep%gpp              =   5.0972823_dp
       sep%gp2              =   2.1945567_dp
       sep%hsp              =   1.5788130_dp
    CASE (20)
       !Calcium
       sep%defined          = .TRUE.
       sep%uss              = -12.3919094_dp
       sep%upp              =  -9.9348289_dp
       sep%beta(0)          =  -8.6404687_dp
       sep%beta(1)          =  -9.9515712_dp
       sep%sto_exponents(0) =   1.0034161_dp
       sep%sto_exponents(1) =   1.3102564_dp
       sep%alp              =   0.4999997_dp
       sep%gss              =   6.5321649_dp
       sep%gsp              =   6.5424442_dp
       sep%gpp              =   6.4627059_dp
       sep%gp2              =   6.3842472_dp
       sep%hsp              =   0.5789676_dp
    CASE (21:29)
    CASE (30)
       !Zinc
       sep%defined          = .TRUE.
       sep%uss              = -20.8397160_dp
       sep%upp              = -19.6252240_dp
       sep%beta(0)          =  -1.0000000_dp
       sep%beta(1)          =  -2.0000000_dp
       sep%sto_exponents(0) =   2.0473590_dp
       sep%sto_exponents(1) =   1.4609460_dp
       sep%alp              =   1.5064570_dp
       sep%gss              =  11.8000000_dp
       sep%gsp              =  11.1820180_dp
       sep%gpp              =  13.3000000_dp
       sep%gp2              =  12.9305200_dp
       sep%hsp              =   0.4846060_dp
    CASE (31)
       !Gallium
       sep%defined          = .TRUE.
       sep%uss              = -28.3044924_dp
       sep%upp              = -27.2063910_dp
       sep%beta(0)          =  -3.9987435_dp
       sep%beta(1)          =  -4.3146711_dp
       sep%sto_exponents(0) =   0.6986316_dp
       sep%sto_exponents(1) =   1.8386933_dp
       sep%alp              =   2.7577991_dp
       sep%gss              =   7.5468114_dp
       sep%gsp              =  10.4697612_dp
       sep%gpp              =   8.4599454_dp
       sep%gp2              =  10.4251148_dp
       sep%hsp              =   1.0628013_dp
    CASE (32)
       !Germanium
       sep%defined          = .TRUE.
       sep%uss              = -33.9493670_dp
       sep%upp              = -27.4251050_dp
       sep%beta(0)          =  -4.5164790_dp
       sep%beta(1)          =  -1.7555170_dp
       sep%sto_exponents(0) =   1.2931800_dp
       sep%sto_exponents(1) =   2.0205640_dp
       sep%alp              =   1.9784980_dp
       sep%gss              =   9.8000000_dp
       sep%gsp              =   8.3000000_dp
       sep%gpp              =   7.3000000_dp
       sep%gp2              =   6.5000000_dp
       sep%hsp              =   1.3000000_dp
    CASE (33)
       !Arsenic
       sep%defined          = .TRUE.
       sep%uss              = -38.6240790_dp
       sep%upp              = -33.4995395_dp
       sep%beta(0)          =  -3.9998231_dp
       sep%beta(1)          =  -4.9056176_dp
       sep%sto_exponents(0) =   2.5614338_dp
       sep%sto_exponents(1) =   1.6117315_dp
       sep%alp              =   1.9381219_dp
       sep%gss              =   6.7464011_dp
       sep%gsp              =   5.6174985_dp
       sep%gpp              =   6.9333925_dp
       sep%gp2              =   6.3054798_dp
       sep%hsp              =   0.5994346_dp
    CASE (34)
       !Selenium
       sep%defined          = .TRUE.
       sep%uss              = -49.8117347_dp
       sep%upp              = -38.0475911_dp
       sep%beta(0)          = -12.4685178_dp
       sep%beta(1)          =  -5.1744376_dp
       sep%sto_exponents(0) =   0.7242956_dp
       sep%sto_exponents(1) =   1.9267288_dp
       sep%alp              =   2.3513155_dp
       sep%gss              =  10.3549483_dp
       sep%gsp              =   5.2801360_dp
       sep%gpp              =   7.3611317_dp
       sep%gp2              =   6.1897284_dp
       sep%hsp              =   0.5996560_dp
    CASE (35)
       !Bromine
       sep%defined          = .TRUE.
       sep%uss              = -99.98644050_dp
       sep%upp              = -75.67130750_dp
       sep%beta(0)          =  -8.91710700_dp
       sep%beta(1)          =  -9.94374000_dp
       sep%sto_exponents(0) =   3.85430190_dp
       sep%sto_exponents(1) =   2.19920910_dp
       sep%alp              =   2.44570510_dp
       sep%gss              =  15.03643948_dp
       sep%gsp              =  13.03468242_dp
       sep%gpp              =  11.27632539_dp
       sep%gp2              =   9.85442552_dp
       sep%hsp              =   2.45586832_dp
    CASE (36)
    CASE (37)
       !Rubidium
       sep%defined          = .TRUE.
       sep%uss              =  -4.3098071_dp
       sep%upp              =  -2.7381921_dp
       sep%beta(0)          =  -2.2062173_dp
       sep%beta(1)          =  -6.2176392_dp
       sep%sto_exponents(0) =   4.0001632_dp
       sep%sto_exponents(1) =   0.9187408_dp
       sep%alp              =   0.9976197_dp
       sep%gss              =  10.7409462_dp
       sep%gsp              =  11.4853623_dp
       sep%gpp              =   8.9878480_dp
       sep%gp2              =   7.7258368_dp
       sep%hsp              =   0.1999294_dp
    CASE (38)
       !Strontium
       sep%defined          = .TRUE.
       sep%uss              =  -10.8451287_dp
       sep%upp              =   -8.3129821_dp
       sep%beta(0)          =   -9.9683427_dp
       sep%beta(1)          =   -9.9946390_dp
       sep%sto_exponents(0) =    1.3729266_dp
       sep%sto_exponents(1) =    1.1118128_dp
       sep%alp              =    0.5082703_dp
       sep%gss              =    4.9305520_dp
       sep%gsp              =    4.4249843_dp
       sep%gpp              =    4.0178045_dp
       sep%gp2              =    4.0335640_dp
       sep%hsp              =    0.6272993_dp
    CASE (39:47)
    CASE(48)
       !Cadmium
       sep%defined          = .TRUE.
       sep%uss              =  -26.1908325_dp
       sep%upp              =  -22.4523736_dp
       sep%beta(0)          =  -11.9610608_dp
       sep%beta(1)          =   -3.9999848_dp
       sep%sto_exponents(0) =    1.4192491_dp
       sep%sto_exponents(1) =    1.0480637_dp
       sep%alp              =    1.1507745_dp
       sep%gss              =   17.2196544_dp
       sep%gsp              =   17.9900180_dp
       sep%gpp              =   19.2007795_dp
       sep%gp2              =   19.0000212_dp
       sep%hsp              =    0.7384511_dp
    CASE (49)
       !Indium
       sep%defined          = .TRUE.
       sep%uss              =  -28.0980892_dp
       sep%upp              =  -19.2780588_dp
       sep%beta(0)          =   -6.3107479_dp
       sep%beta(1)          =   -2.7025837_dp
       sep%sto_exponents(0) =    1.7625740_dp
       sep%sto_exponents(1) =    1.8648962_dp
       sep%alp              =    2.3438756_dp
       sep%gss              =    9.4928794_dp
       sep%gsp              =    7.0094241_dp
       sep%gpp              =    9.6640986_dp
       sep%gp2              =    7.0100315_dp
       sep%hsp              =    0.5995894_dp
    CASE (50)
       !Tin
       sep%defined          = .TRUE.
       sep%uss              = -40.8518020_dp
       sep%upp              = -28.5602490_dp
       sep%beta(0)          =  -3.2351470_dp
       sep%beta(1)          =  -4.2904160_dp
       sep%sto_exponents(0) =   2.0803800_dp
       sep%sto_exponents(1) =   1.9371060_dp
       sep%alp              =   1.8008140_dp
       sep%gss              =   9.8000000_dp
       sep%gsp              =   8.3000000_dp
       sep%gpp              =   7.3000000_dp
       sep%gp2              =   6.5000000_dp
       sep%hsp              =   1.3000000_dp
    CASE (51)
       !Antimony
       sep%defined          = .TRUE.
       sep%uss              = -42.0643435_dp
       sep%upp              = -35.0626031_dp
       sep%beta(0)          =  -0.9999715_dp
       sep%beta(1)          =  -4.0920176_dp
       sep%sto_exponents(0) =   3.6458835_dp
       sep%sto_exponents(1) =   1.9733156_dp
       sep%alp              =   1.9763403_dp
       sep%gss              =  10.6739308_dp
       sep%gsp              =   7.0477648_dp
       sep%gpp              =   6.7446162_dp
       sep%gp2              =   6.3408531_dp
       sep%hsp              =   0.5997512_dp
    CASE (52)
       !Tellurium
       sep%defined          = .TRUE.
       sep%uss              = -84.2274722_dp
       sep%upp              = -46.5332871_dp
       sep%beta(0)          =  -8.5622652_dp
       sep%beta(1)          =  -2.6942963_dp
       sep%sto_exponents(0) =   2.7461609_dp
       sep%sto_exponents(1) =   1.6160376_dp
       sep%alp              =   2.2924145_dp
       sep%gss              =   5.1367706_dp
       sep%gsp              =  11.0720752_dp
       sep%gpp              =   5.8447934_dp
       sep%gp2              =   5.0720495_dp
       sep%hsp              =   0.5997994_dp
    CASE (53)
       !Iodine
       sep%defined          = .TRUE.
       sep%uss              =-100.00305380_dp
       sep%upp              = -74.61146920_dp
       sep%beta(0)          =  -7.41445100_dp
       sep%beta(1)          =  -6.19678100_dp
       sep%sto_exponents(0) =   2.27296100_dp
       sep%sto_exponents(1) =   2.16949800_dp
       sep%alp              =   2.20732000_dp
       sep%gss              =  15.04044855_dp
       sep%gsp              =  13.05655798_dp
       sep%gpp              =  11.14778369_dp
       sep%gp2              =   9.91409071_dp
       sep%hsp              =   2.45638202_dp
    CASE (54)
    CASE (55)
       !Cesium
       sep%defined          = .TRUE.
       sep%uss              = -3.2184078_dp
       sep%upp              = -1.7699113_dp
       sep%beta(0)          = -1.6043600_dp
       sep%beta(1)          = -4.2698040_dp
       sep%sto_exponents(0) =  6.0004170_dp
       sep%sto_exponents(1) =  0.8986916_dp
       sep%alp              =  0.4981646_dp
       sep%gss              =  7.6447851_dp
       sep%gsp              =  3.0454989_dp
       sep%gpp              = 10.0000745_dp
       sep%gp2              =  6.1761092_dp
       sep%hsp              =  0.4647853_dp
    CASE (56)
       !Barium
       sep%defined          = .TRUE.
       sep%uss              = -10.1125345_dp
       sep%upp              =  -8.2347224_dp
       sep%beta(0)          =  -9.9994459_dp
       sep%beta(1)          =  -9.6197255_dp
       sep%sto_exponents(0) =   1.9765973_dp
       sep%sto_exponents(1) =   1.3157348_dp
       sep%alp              =   0.8594840_dp
       sep%gss              =   4.8486178_dp
       sep%gsp              =   4.5659982_dp
       sep%gpp              =   5.0937708_dp
       sep%gp2              =   5.2125824_dp
       sep%hsp              =   0.5237082_dp
    CASE (57:79)
    CASE (80)
       !Mercury
       sep%defined          = .TRUE.
       sep%uss              = -19.8095740_dp
       sep%upp              = -13.1025300_dp
       sep%beta(0)          =  -0.4045250_dp
       sep%beta(1)          =  -6.2066830_dp
       sep%sto_exponents(0) =   2.2181840_dp
       sep%sto_exponents(1) =   2.0650380_dp
       sep%alp              =   1.3356410_dp
       sep%gss              =  10.8000000_dp
       sep%gsp              =   9.3000000_dp
       sep%gpp              =  14.3000000_dp
       sep%gp2              =  13.5000000_dp
       sep%hsp              =   1.3000000_dp
    CASE (81)
       !Thallium
       sep%defined          = .TRUE.
       sep%uss              = -29.7009655_dp
       sep%upp              = -29.5937539_dp
       sep%beta(0)          =  -4.9667442_dp
       sep%beta(1)          =  -7.7616060_dp
       sep%sto_exponents(0) =   4.0000447_dp
       sep%sto_exponents(1) =   1.8076332_dp
       sep%alp              =   1.3116968_dp
       sep%gss              =   8.8675337_dp
       sep%gsp              =  12.1148290_dp
       sep%gpp              =  10.6532769_dp
       sep%gp2              =  13.5333191_dp
       sep%hsp              =   0.5997565_dp
    CASE (82)
       !Lead
       sep%defined          = .TRUE.
       sep%uss              = -47.3196920_dp
       sep%upp              = -28.8475600_dp
       sep%beta(0)          =  -8.0423870_dp
       sep%beta(1)          =  -3.0000000_dp
       sep%sto_exponents(0) =   2.4982860_dp
       sep%sto_exponents(1) =   2.0820710_dp
       sep%alp              =   1.7283330_dp
       sep%gss              =   9.8000000_dp
       sep%gsp              =   8.3000000_dp
       sep%gpp              =   7.3000000_dp
       sep%gp2              =   6.5000000_dp
       sep%hsp              =   1.3000000_dp
    CASE (83)
       !Bismuth
       sep%defined          = .TRUE.
       sep%uss              = -53.5827147_dp
       sep%upp              = -39.4572213_dp
       sep%beta(0)          =  -9.0000249_dp
       sep%beta(1)          =  -1.9830269_dp
       sep%sto_exponents(0) =   2.6772255_dp
       sep%sto_exponents(1) =   0.6936864_dp
       sep%alp              =   5.7660628_dp
       sep%gss              =   8.3702778_dp
       sep%gsp              =   7.7974668_dp
       sep%gpp              =   9.8303165_dp
       sep%gp2              =   8.9291355_dp
       sep%hsp              =   0.5999908_dp
    CASE (84:nelem)
    END SELECT
    ! In case overwrite with MNDO-d parameters
    IF (itype==do_method_mndod) CALL  mndod_default_parameter(sep, z, error)
    ! Finalize parameters for the element
    CALL valence_electrons(sep, extended_basis_set=(itype==do_method_mndod),error=error)
    CALL calpar(z,sep,error)
    CALL convert_param_to_cp2k(sep,error)

  END SUBROUTINE mndo_default_parameter


! *****************************************************************************
!> \brief Default parameter sets for semi empirical models developed at PNNL
!> \param sep ...
!> \param z ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pnnl_default_parameter ( sep, z, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER, INTENT(IN)                      :: z
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pnnl_default_parameter', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: Ang = 0.529167_dp, &
                                                eV = 27.21_dp

    sep%typ    = do_method_pnnl
    sep%name   = ptable(z)%name
    sep%eheat  = ptable(z)%heat_of_formation / kcalmol
    sep%z      = z

    SELECT CASE (z)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            "invalid atomic number")
    CASE (0)
       sep%defined          = .TRUE.
    CASE (1)
       !Hydrogen
       sep%defined          = .TRUE.
       sep%alp              =   3.323150026265647_dp*Ang/angstrom
       sep%beta(0)          =  -6.070138572957370_dp*evolt/eV
       sep%beta(1)          =   0.0_dp              *evolt/eV
       sep%sto_exponents(0) =   1.220643961907751_dp
       sep%sto_exponents(1) =   1.351277996861119_dp
       sep%uss              =  -11.185667_dp        *evolt/eV
       sep%upp              =   1000000.0_dp        *evolt/eV
       sep%gss              =   14.794208_dp        *evolt/eV
       sep%gsp              =   0.0_dp              *evolt/eV
       sep%gpp              =   0.0_dp              *evolt/eV
       sep%gp2              =   0.0_dp              *evolt/eV
       sep%hsp              =   0.0_dp              *evolt/eV
       sep%ass              =   1.0e-16_dp
       sep%asp              =   3.355792629670590_dp
       sep%app              =   0.03941680815458009_dp
       sep%de               =   20.06411392885370_dp      /eV
       sep%nr               =   1
       sep%acoul            =   16.27368453086666_dp      /eV
       ! This parameterization explicitly supports
       ! p-orbitals for H. So by default we turn them
       ! on for this element.
       sep%p_orbitals_on_h  = .TRUE.
    CASE (8)
       !Oxygen
       sep%defined          = .TRUE.
       sep%alp              =   2.980042911001686_dp*Ang/angstrom
       sep%beta(0)          = -40.62654394990248_dp *evolt/eV
       sep%beta(1)          = -28.39467965792422_dp *evolt/eV
       sep%sto_exponents(0) =   3.814650982832301_dp
       sep%sto_exponents(1) =   2.446497244623426_dp
       sep%uss              = -79.081852_dp         *evolt/eV
       sep%upp              = -72.59032_dp          *evolt/eV
       sep%gss              =  15.75576_dp          *evolt/eV
       sep%gsp              =  10.62116_dp          *evolt/eV
       sep%gpp              =  13.654016_dp         *evolt/eV
       sep%gp2              =  12.406095_dp         *evolt/eV
       sep%hsp              =   0.593883_dp         *evolt/eV
       sep%ass              =  1.0e-16_dp
       sep%asp              =  3.583554483196465_dp
       sep%app              =  0.0008605218877175237_dp
       sep%de               =  7.159002319741392_dp       /eV
       sep%nr               =  2
       sep%acoul            =  16.22361790703965_dp       /eV
    END SELECT
    ! Finalize parameters for the element
    CALL valence_electrons(sep, extended_basis_set=.FALSE., error=error)
    CALL calpar(z,sep,error)
    CALL convert_param_to_cp2k(sep,error)

  END SUBROUTINE pnnl_default_parameter

! *****************************************************************************
!> \brief Additional default parameter sets for semi empirical models: MNDO-d
!> \param sep ...
!> \param z ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE mndod_default_parameter ( sep, z, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mndod_default_parameter', &
      routineP = moduleN//':'//routineN

    SELECT CASE (z)
    CASE (11)
       ! Sodium
       sep%defined          = .TRUE.
       sep%uss              =  -5.20100000_dp
       sep%upp              =  -2.71257317_dp
       sep%beta(0)          =  -1.08738166_dp
       sep%beta(1)          =  -0.48623935_dp
       sep%sto_exponents(0) =   0.98750834_dp
       sep%sto_exponents(1) =   0.89334983_dp
       sep%zn(0)            =   0.65411258_dp
       sep%zn(1)            =   0.56440874_dp
       sep%alp              =   1.17010202_dp
       sep%gss              =   4.59444476_dp
       sep%gsp              =   4.14757400_dp
       sep%gpp              =   4.29919761_dp
       sep%gp2              =   3.79695732_dp
       sep%hsp              =   0.53440874_dp
    CASE (12)
       ! Magnesium
       sep%defined          = .TRUE.
       sep%uss              = -15.09700000_dp
       sep%upp              = -10.65000000_dp
       sep%beta(0)          =  -1.89588355_dp
       sep%beta(1)          =  -2.14108943_dp
       sep%sto_exponents(0) =   1.44890446_dp
       sep%sto_exponents(1) =   0.95293002_dp
       sep%zn(0)            =   1.05000000_dp
       sep%zn(1)            =   0.92527190_dp
       sep%alp              =   1.62146984_dp
       sep%gss              =   7.37513258_dp
       sep%gsp              =   6.88890741_dp
       sep%gpp              =   7.04795383_dp
       sep%gp2              =   6.22459871_dp
       sep%hsp              =   0.72673390_dp
    CASE (13)
       !Aluminium
       sep%defined          = .TRUE.
       sep%uss              = -24.01792910_dp
       sep%upp              = -20.79597967_dp
       sep%udd              =  -5.22082737_dp
       sep%beta(0)          =  -7.10185851_dp
       sep%beta(1)          =  -2.31809618_dp
       sep%beta(2)          =  -3.35638545_dp
       sep%sto_exponents(0) =   1.79402273_dp
       sep%sto_exponents(1) =   1.37130919_dp
       sep%sto_exponents(2) =   0.80591133_dp
       sep%zn(0)            =   1.22249269_dp
       sep%zn(1)            =   1.09291990_dp
       sep%zn(2)            =   0.80038285_dp
       sep%alp              =   1.44301676_dp
       sep%gss              =   8.58671016_dp
       sep%gsp              =   7.66469306_dp
       sep%gpp              =   8.32495723_dp
       sep%gp2              =   7.35242020_dp
       sep%hsp              =   0.54401293_dp
    CASE (14)
       !Silicon
       sep%defined          = .TRUE.
       sep%uss              = -36.05153000_dp
       sep%upp              = -27.53569100_dp
       sep%udd              = -14.67743900_dp
       sep%beta(0)          =  -8.21073420_dp
       sep%beta(1)          =  -4.88462030_dp
       sep%beta(2)          =  -2.60801150_dp
       sep%sto_exponents(0) =   1.91565460_dp
       sep%sto_exponents(1) =   1.68161130_dp
       sep%sto_exponents(2) =   0.96677166_dp
       sep%zn(0)            =   1.52929180_dp
       sep%zn(1)            =   0.97628075_dp
       sep%zn(2)            =   0.93816441_dp
       sep%alp              =   1.66006930_dp
       sep%gss              =  10.74164700_dp
       sep%gsp              =   7.56066640_dp
       sep%gpp              =   7.43649690_dp
       sep%gp2              =   6.56775150_dp
       sep%hsp              =   0.87753880_dp
    CASE (15)
       !Phosphorus
       sep%defined          = .TRUE.
       sep%uss              = -47.05552900_dp
       sep%upp              = -38.06705900_dp
       sep%udd              = -23.69159700_dp
       sep%beta(0)          =  -8.90210430_dp
       sep%beta(1)          =  -9.38611080_dp
       sep%beta(2)          =  -2.09170080_dp
       sep%sto_exponents(0) =   2.26646290_dp
       sep%sto_exponents(1) =   1.94001490_dp
       sep%sto_exponents(2) =   1.10010900_dp
       sep%zn(0)            =   1.63437610_dp
       sep%zn(1)            =   1.08291170_dp
       sep%zn(2)            =   1.00651470_dp
       sep%alp              =   1.85255120_dp
       sep%gss              =   11.4797530_dp
       sep%gsp              =   8.55756910_dp
       sep%gpp              =   8.24872280_dp
       sep%gp2              =   7.28509170_dp
       sep%hsp              =   2.10780440_dp
    CASE (16)
       !Sulphur
       sep%defined          = .TRUE.
       sep%uss              = -56.88912800_dp
       sep%upp              = -47.27474500_dp
       sep%udd              = -25.09511800_dp
       sep%beta(0)          = -10.99954500_dp
       sep%beta(1)          = -12.21543700_dp
       sep%beta(2)          =  -1.88066950_dp
       sep%sto_exponents(0) =   2.22585050_dp
       sep%sto_exponents(1) =   2.09970560_dp
       sep%sto_exponents(2) =   1.23147250_dp
       sep%zn(0)            =   1.73639140_dp
       sep%zn(1)            =   1.12118170_dp
       sep%zn(2)            =   1.05084670_dp
       sep%alp              =   2.02305950_dp
       sep%gss              =  12.19630200_dp
       sep%gsp              =   8.85390090_dp
       sep%gpp              =   8.54023240_dp
       sep%gp2              =   7.54254650_dp
       sep%hsp              =   2.64635230_dp
    CASE (17)
       !Chlorine
       sep%defined          = .TRUE.
       sep%uss              = -69.62297275_dp
       sep%upp              = -59.10072899_dp
       sep%udd              = -36.67457320_dp
       sep%beta(0)          =  -6.03729165_dp
       sep%beta(1)          = -19.18338497_dp
       sep%beta(2)          =  -1.87778198_dp
       sep%sto_exponents(0) =   2.56161065_dp
       sep%sto_exponents(1) =   2.38933800_dp
       sep%sto_exponents(2) =   1.25139777_dp
       sep%zn(0)            =   1.88087547_dp
       sep%zn(1)            =   1.18104227_dp
       sep%zn(2)            =   1.14061555_dp
       sep%alp              =   2.18030019_dp
       sep%gss              =  13.21114854_dp
       sep%gsp              =   9.41949513_dp
       sep%gpp              =   8.99620033_dp
       sep%gp2              =   7.94524745_dp
       sep%hsp              =   3.08149862_dp
    CASE (30)
       !Zinc
       sep%defined          = .TRUE.
       sep%uss              = -18.02300000_dp
       sep%upp              = -12.24216585_dp
       sep%beta(0)          =  -5.01726076_dp
       sep%beta(1)          =  -0.71205972_dp
       sep%sto_exponents(0) =   1.73150352_dp
       sep%sto_exponents(1) =   1.39358305_dp
       sep%zn(0)            =   1.56600000_dp
       sep%zn(1)            =   0.862839817_dp
       sep%alp              =   1.51763697_dp
       sep%gss              =   8.56072836_dp
       sep%gsp              =   7.49003598_dp
       sep%gpp              =   5.13964830_dp
       sep%gp2              =   4.50540309_dp
       sep%hsp              =   0.53294610_dp
    CASE (35)
       !Bromine
       sep%defined          = .TRUE.
       sep%uss              = -65.40277790_dp
       sep%upp              = -54.55375352_dp
       sep%udd              = -13.72809929_dp
       sep%beta(0)          =  -8.31497607_dp
       sep%beta(1)          = -10.50704145_dp
       sep%beta(2)          =  -0.96259930_dp
       sep%sto_exponents(0) =   2.59054101_dp
       sep%sto_exponents(1) =   2.33085649_dp
       sep%sto_exponents(2) =   1.35736115_dp
       sep%zn(0)            =   2.23581544_dp
       sep%zn(1)            =   1.43292654_dp
       sep%zn(2)            =   1.24257826_dp
       sep%alp              =   2.09105000_dp
       sep%gss              =  12.22235546_dp
       sep%gsp              =   8.26372010_dp
       sep%gpp              =   8.53546437_dp
       sep%gp2              =   7.48216712_dp
       sep%hsp              =   2.74952230_dp
    CASE (48)
       ! Cadmium
       sep%defined          = .TRUE.
       sep%uss              = -16.96970000_dp
       sep%upp              = -12.40096476_dp
       sep%beta(0)          =  -2.77154379_dp
       sep%beta(1)          =  -1.80565019_dp
       sep%sto_exponents(0) =   1.74880559_dp
       sep%sto_exponents(1) =   1.56321473_dp
       sep%zn(0)            =   1.76314840_dp
       sep%zn(1)            =   1.52551900_dp
       sep%alp              =   1.42461329_dp
       sep%gss              =   7.90443438_dp
       sep%gsp              =   7.51570687_dp
       sep%gpp              =   7.47999993_dp
       sep%gp2              =   6.51866416_dp
       sep%hsp              =   0.63674441_dp
    CASE (53)
       !Iodine
       sep%defined          = .TRUE.
       sep%uss              = -62.76535256_dp
       sep%upp              = -50.29211568_dp
       sep%udd              = -12.24830501_dp
       sep%beta(0)          = -10.69948666_dp
       sep%beta(1)          =  -4.94117776_dp
       sep%beta(2)          =  -2.35046098_dp
       sep%sto_exponents(0) =   2.75654324_dp
       sep%sto_exponents(1) =   2.25307954_dp
       sep%sto_exponents(2) =   1.50233509_dp
       sep%zn(0)            =   2.67241100_dp
       sep%zn(1)            =   1.57229871_dp
       sep%zn(2)            =   1.25884802_dp
       sep%alp              =   1.90617441_dp
       sep%gss              =  11.98078196_dp
       sep%gsp              =   7.85590192_dp
       sep%gpp              =   7.70937227_dp
       sep%gp2              =   6.71855729_dp
       sep%hsp              =   2.07147462_dp
    CASE (80)
       !Mercury
       sep%defined          = .TRUE.
       sep%uss              = -18.81564903_dp
       sep%upp              = -13.39711352_dp
       sep%beta(0)          =  -2.21872239_dp
       sep%beta(1)          =  -2.90978573_dp
       sep%sto_exponents(0) =   2.33310757_dp
       sep%sto_exponents(1) =   1.70831069_dp
       sep%zn(0)            =   2.18600011_dp
       sep%zn(1)            =   1.70500462_dp
       sep%alp              =   1.38224172_dp
       sep%gss              =   8.31564948_dp
       sep%gsp              =   8.21217300_dp
       sep%gpp              =   7.11525878_dp
       sep%gp2              =   6.17124983_dp
       sep%hsp              =   0.83594100_dp
    END SELECT

  END SUBROUTINE mndod_default_parameter

! *****************************************************************************
!> \brief Default parameter sets for semi empirical models: AM1
!> \param sep ...
!> \param z ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE am1_default_parameter ( sep, z, error )
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'am1_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ    = do_method_am1
    sep%name   = ptable(z)%name
    sep%eheat  = ptable(z)%heat_of_formation / kcalmol
    sep%z      = z

    SELECT CASE (z)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            "invalid atomic number")
    CASE (0)
       sep%defined = .TRUE.
    CASE (1)
       !Hydrogen
       sep%defined          = .TRUE.
       sep%alp              =   2.8823240_dp
       sep%beta(0)          =  -6.1737870_dp
       sep%sto_exponents(0) =   1.1880780_dp
       sep%uss    = -11.3964270_dp
       sep%gss    =  12.8480000_dp
       sep%fn1(1) =  0.1227960_dp; sep%fn2(1) = 5.0000000_dp; sep%fn3(1) = 1.2000000_dp
       sep%fn1(2) =  0.0050900_dp; sep%fn2(2) = 5.0000000_dp; sep%fn3(2) = 1.8000000_dp
       sep%fn1(3) = -0.0183360_dp; sep%fn2(3) = 2.0000000_dp; sep%fn3(3) = 2.1000000_dp
    CASE (2)
    CASE (3)
       !Lithium
       sep%defined          = .TRUE.
       sep%alp              =   1.2501400_dp
       sep%beta(0)          =  -1.3500400_dp
       sep%beta(1)          =  -1.3500400_dp
       sep%sto_exponents(0) =   0.7023800_dp
       sep%sto_exponents(1) =   0.7023800_dp
       sep%uss    =  -5.1280000_dp
       sep%upp    =  -2.7212000_dp
       sep%gss    =   7.3000000_dp
       sep%gsp    =   5.4200000_dp
       sep%gpp    =   5.0000000_dp
       sep%gp2    =   4.5200000_dp
       sep%hsp    =   0.8300000_dp
    CASE (4)
       !Beryllium
       sep%defined          = .TRUE.
       sep%alp              =   1.6694340_dp
       sep%beta(0)          =  -4.0170960_dp
       sep%beta(1)          =  -4.0170960_dp
       sep%sto_exponents(0) =   1.0042100_dp
       sep%sto_exponents(1) =   1.0042100_dp
       sep%uss    = -16.6023780_dp
       sep%upp    = -10.7037710_dp
       sep%gss    =   9.0000000_dp
       sep%gsp    =   7.4300000_dp
       sep%gpp    =   6.9700000_dp
       sep%gp2    =   6.2200000_dp
       sep%hsp    =   1.2800000_dp
    CASE (5)
       !Boron
       sep%defined          = .TRUE.
       sep%alp              =   2.4469090_dp
       sep%beta(0)          =  -9.5991140_dp
       sep%beta(1)          =  -6.2737570_dp
       sep%sto_exponents(0) =   1.6117090_dp
       sep%sto_exponents(1) =   1.5553850_dp
       sep%uss   = -34.4928700_dp
       sep%upp   = -22.6315250_dp
       sep%gss   =  10.5900000_dp
       sep%gsp   =   9.5600000_dp
       sep%gpp   =   8.8600000_dp
       sep%gp2   =   7.8600000_dp
       sep%hsp   =   1.8100000_dp
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
       sep%defined          = .TRUE.
       sep%alp              =    2.6482740_dp
       sep%beta(0)          =  -15.7157830_dp
       sep%beta(1)          =   -7.7192830_dp
       sep%sto_exponents(0) =    1.8086650_dp
       sep%sto_exponents(1) =    1.6851160_dp
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
       sep%defined          = .TRUE.
       sep%alp              =    2.9472860_dp
       sep%beta(0)          =  -20.2991100_dp
       sep%beta(1)          =  -18.2386660_dp
       sep%sto_exponents(0) =    2.3154100_dp
       sep%sto_exponents(1) =    2.1579400_dp
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
       sep%defined          = .TRUE.
       sep%alp              =    4.4553710_dp
       sep%beta(0)          =  -29.2727730_dp
       sep%beta(1)          =  -29.2727730_dp
       sep%sto_exponents(0) =    3.1080320_dp
       sep%sto_exponents(1) =    2.5240390_dp
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
       sep%defined          = .TRUE.
       sep%alp              =    5.5178000_dp
       sep%beta(0)          =  -69.5902770_dp
       sep%beta(1)          =  -27.9223600_dp
       sep%sto_exponents(0) =    3.7700820_dp
       sep%sto_exponents(1) =    2.4946700_dp
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
       sep%defined          = .TRUE.
       sep%uss              =  -14.96959313_dp
       sep%upp              =  -11.56229248_dp
       sep%beta(0)          =   -1.25974355_dp
       sep%beta(1)          =   -0.77836604_dp
       sep%sto_exponents(0) =    1.22339270_dp
       sep%sto_exponents(1) =    1.02030798_dp
       sep%alp   =    1.67049799_dp
       sep%gss   =    7.50132277_dp
       sep%gsp   =    6.34591536_dp
       sep%gpp   =    4.77534467_dp
       sep%gp2   =    4.34017279_dp
       sep%hsp   =    0.48930466_dp
       sep%fn1(1) =  2.55017735_dp; sep%fn2(1) = 4.29397225_dp; sep%fn3(1) = 0.79989601_dp
       sep%fn1(2) = -0.00565806_dp; sep%fn2(2) = 2.96053910_dp; sep%fn3(2) = 1.47499983_dp
       sep%fn1(3) = -0.00610286_dp; sep%fn2(3) = 2.61416919_dp; sep%fn3(3) = 2.42604040_dp
    CASE (13)
       !Aluminium
       sep%defined          = .TRUE.
       sep%alp              =   1.9765860_dp
       sep%beta(0)          =  -3.8668220_dp
       sep%beta(1)          =  -2.3171460_dp
       sep%sto_exponents(0) =   1.5165930_dp
       sep%sto_exponents(1) =   1.3063470_dp
       sep%sto_exponents(2) =   1.0000000_dp
       sep%uss   = -24.3535850_dp
       sep%upp   = -18.3636450_dp
       sep%gss   =   8.0900000_dp
       sep%gsp   =   6.6300000_dp
       sep%gpp   =   5.9800000_dp
       sep%gp2   =   5.4000000_dp
       sep%hsp   =   0.7000000_dp
       sep%fn1(1) = 0.0900000_dp; sep%fn2(1) = 12.3924430_dp; sep%fn3(1) = 2.0503940_dp
    CASE (14)
       !Silicon
       sep%defined          = .TRUE.
       sep%alp              =   2.257816_dp
       sep%beta(0)          =  -3.784852_dp
       sep%beta(1)          =  -1.968123_dp
       sep%sto_exponents(0) =   1.830697_dp
       sep%sto_exponents(1) =   1.2849530_dp
       sep%sto_exponents(2) =   1.0000000_dp
       sep%uss   = -33.9536220_dp
       sep%upp   = -28.9347490_dp
       sep%gss   =   9.8200000_dp
       sep%gsp   =   8.3600000_dp
       sep%gpp   =   7.3100000_dp
       sep%gp2   =   6.5400000_dp
       sep%hsp   =   1.3200000_dp
       sep%fn1(1) =  0.25000000_dp; sep%fn2(1) = 9.00000000_dp; sep%fn3(1) = 0.91145300_dp
       sep%fn1(2) =  0.06151300_dp; sep%fn2(2) = 5.00000000_dp; sep%fn3(2) = 1.99556900_dp
       sep%fn1(3) =  0.00207890_dp; sep%fn2(3) = 5.00000000_dp; sep%fn3(3) = 2.99061000_dp
    CASE (15)
       !Phosphorus
       sep%defined          = .TRUE.
       sep%alp              =   2.4553220_dp
       sep%beta(0)          =  -6.3537640_dp
       sep%beta(1)          =  -6.5907090_dp
       sep%sto_exponents(0) =   1.9812800_dp
       sep%sto_exponents(1) =   1.8751500_dp
       sep%sto_exponents(2) =   1.0000000_dp
       sep%uss   =  -42.0298630_dp
       sep%upp   =  -34.0307090_dp
       sep%gss   =   11.5600050_dp
       sep%gsp   =    5.2374490_dp
       sep%gpp   =    7.8775890_dp
       sep%gp2   =    7.3076480_dp
       sep%hsp   =    0.7792380_dp
       sep%fn1(1) = -0.03182700_dp; sep%fn2(1) = 6.00000000_dp; sep%fn3(1) = 1.47432300_dp
       sep%fn1(2) =  0.01847000_dp; sep%fn2(2) = 7.00000000_dp; sep%fn3(2) = 1.77935400_dp
       sep%fn1(3) =  0.03329000_dp; sep%fn2(3) = 9.00000000_dp; sep%fn3(3) = 3.00657600_dp
    CASE (16)
       !Sulphur
       sep%defined          = .TRUE.
       sep%alp              =   2.4616480_dp
       sep%beta(0)          =  -3.9205660_dp
       sep%beta(1)          =  -7.9052780_dp
       sep%sto_exponents(0) =   2.3665150_dp
       sep%sto_exponents(1) =   1.6672630_dp
       sep%sto_exponents(2) =   1.0000000_dp
       sep%uss   =  -56.6940560_dp
       sep%upp   =  -48.7170490_dp
       sep%gss   =   11.7863290_dp
       sep%gsp   =    8.6631270_dp
       sep%gpp   =   10.0393080_dp
       sep%gp2   =    7.7816880_dp
       sep%hsp   =    2.5321370_dp
       sep%fn1(1) = -0.50919500_dp; sep%fn2(1) = 4.59369100_dp; sep%fn3(1) = 0.77066500_dp
       sep%fn1(2) = -0.01186300_dp; sep%fn2(2) = 5.86573100_dp; sep%fn3(2) = 1.50331300_dp
       sep%fn1(3) =  0.01233400_dp; sep%fn2(3) =13.55733600_dp; sep%fn3(3) = 2.00917300_dp
    CASE (17)
       !Chlorine
       sep%defined          = .TRUE.
       sep%alp              =    2.9193680_dp
       sep%beta(0)          =  -24.5946700_dp
       sep%beta(1)          =  -14.6372160_dp
       sep%sto_exponents(0) =    3.6313760_dp
       sep%sto_exponents(1) =    2.0767990_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%uss   = -111.6139480_dp
       sep%upp   =  -76.6401070_dp
       sep%gss   =   15.0300000_dp
       sep%gsp   =   13.1600000_dp
       sep%gpp   =   11.3000000_dp
       sep%gp2   =    9.9700000_dp
       sep%hsp   =    2.4200000_dp
       sep%fn1(1) = 0.0942430_dp; sep%fn2(1) = 4.0000000_dp; sep%fn3(1) = 1.3000000_dp
       sep%fn1(2) = 0.0271680_dp; sep%fn2(2) = 4.0000000_dp; sep%fn3(2) = 2.1000000_dp
    CASE (18:29)
    CASE (30)
       !Zinc
       sep%defined          = .TRUE.
       sep%uss              =  -21.0400080_dp
       sep%upp              =  -17.6555740_dp
       sep%beta(0)          =   -1.9974290_dp
       sep%beta(1)          =   -4.7581190_dp
       sep%sto_exponents(0) =    1.9542990_dp
       sep%sto_exponents(1) =    1.3723650_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp    =      1.4845630_dp
       sep%gss    =     11.8000000_dp
       sep%gsp    =     11.1820180_dp
       sep%gpp    =     13.3000000_dp
       sep%gp2    =     12.9305200_dp
       sep%hsp    =      0.4846060_dp
    CASE (31)
    CASE (32)
       !Germanium
       sep%defined          = .TRUE.
       sep%uss              = -34.1838890_dp
       sep%upp              = -28.6408110_dp
       sep%beta(0)          =  -4.3566070_dp
       sep%beta(1)          =  -0.9910910_dp
       sep%sto_exponents(0) =   1.2196310_dp
       sep%sto_exponents(1) =   1.9827940_dp
       sep%sto_exponents(2) =   1.0000000_dp
       sep%alp    =      2.1364050_dp
       sep%gss    =     10.1686050_dp
       sep%gsp    =      8.1444730_dp
       sep%gpp    =      6.6719020_dp
       sep%gp2    =      6.2697060_dp
       sep%hsp    =      0.9370930_dp
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
       sep%defined          = .TRUE.
       sep%alp              =    2.5765460_dp
       sep%beta(0)          =  -19.3998800_dp
       sep%beta(1)          =   -8.9571950_dp
       sep%sto_exponents(0) =    3.0641330_dp
       sep%sto_exponents(1) =    2.0383330_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%uss   = -104.6560630_dp
       sep%upp   =  -74.9300520_dp
       sep%gss   =   15.0364395_dp
       sep%gsp   =   13.0346824_dp
       sep%gpp   =   11.2763254_dp
       sep%gp2   =    9.8544255_dp
       sep%hsp   =    2.4558683_dp
       sep%fn1(1) = 0.0666850_dp; sep%fn2(1) = 4.0000000_dp; sep%fn3(1) = 1.5000000_dp
       sep%fn1(2) = 0.0255680_dp; sep%fn2(2) = 4.0000000_dp; sep%fn3(2) = 2.3000000_dp
    CASE (36:49)
    CASE (50)
       !Tin
       sep%defined          = .TRUE.
       sep%alp              =    1.8369360_dp
       sep%beta(0)          =   -3.2350000_dp
       sep%beta(1)          =   -2.5778900_dp
       sep%sto_exponents(0) =    2.5993760_dp
       sep%sto_exponents(1) =    1.6959620_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%uss   =  -35.4967410_dp
       sep%upp   =  -28.0976360_dp
       sep%gss   =    9.8000000_dp
       sep%gsp   =    8.3000000_dp
       sep%gpp   =    7.3000000_dp
       sep%gp2   =    6.5000000_dp
       sep%hsp   =    1.3000000_dp
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
       sep%defined          = .TRUE.
       sep%uss              = -103.5896630_dp
       sep%upp              =  -74.4299970_dp
       sep%beta(0)          =   -8.4433270_dp
       sep%beta(1)          =   -6.3234050_dp
       sep%sto_exponents(0) =    2.1028580_dp
       sep%sto_exponents(1) =    2.1611530_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    2.2994240_dp
       sep%gss   =   15.0404486_dp
       sep%gsp   =   13.0565580_dp
       sep%gpp   =   11.1477837_dp
       sep%gp2   =    9.9140907_dp
       sep%hsp   =    2.4563820_dp
       sep%fn1(1) = 0.0043610_dp; sep%fn2(1) = 2.3000000_dp; sep%fn3(1) = 1.8000000_dp
       sep%fn1(2) = 0.0157060_dp; sep%fn2(2) = 3.0000000_dp; sep%fn3(2) = 2.2400000_dp
    CASE (80)
       !Mercury
       sep%defined          = .TRUE.
       sep%uss              =  -19.9415780_dp
       sep%upp              =  -11.1108700_dp
       sep%beta(0)          =   -0.9086570_dp
       sep%beta(1)          =   -4.9093840_dp
       sep%sto_exponents(0) =    2.0364130_dp
       sep%sto_exponents(1) =    1.9557660_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    1.4847340_dp
       sep%gss   =   10.8000000_dp
       sep%gsp   =    9.3000000_dp
       sep%gpp   =   14.3000000_dp
       sep%gp2   =   13.5000000_dp
       sep%hsp   =    1.3000000_dp
    CASE (81)
    CASE (82)
       !Lead
       sep%defined          = .TRUE.
       sep%uss              =  -47.3196920_dp
       sep%upp              =  -28.8475600_dp
       sep%beta(0)          =   -8.0423870_dp
       sep%beta(1)          =   -3.0000000_dp
       sep%sto_exponents(0) =    2.4982860_dp
       sep%sto_exponents(1) =    1.7283330_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    1.7283330_dp
       sep%gss   =    9.8000000_dp
       sep%gsp   =    8.3000000_dp
       sep%gpp   =    7.3000000_dp
       sep%gp2   =    6.5000000_dp
       sep%hsp   =    1.3000000_dp
    CASE (83:nelem)
    END SELECT
    ! Finalize parameters for the element
    CALL valence_electrons(sep, extended_basis_set=.FALSE., error=error)
    CALL calpar(z,sep,error)
    CALL convert_param_to_cp2k(sep,error)

  END SUBROUTINE am1_default_parameter

! *****************************************************************************
!> \brief Default parameter sets for semi empirical models: RM1
!> \param sep ...
!> \param z ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rm1_default_parameter ( sep, z, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rm1_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ    = do_method_rm1
    sep%name   = ptable(z)%name
    sep%eheat  = ptable(z)%heat_of_formation / kcalmol
    sep%z      = z

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
      CASE (80:nelem)
    END SELECT
    ! Finalize parameters for the element
    CALL valence_electrons(sep, extended_basis_set=.FALSE., error=error)
    CALL calpar(z,sep,error)
    CALL convert_param_to_cp2k(sep,error)

  END SUBROUTINE rm1_default_parameter

! *****************************************************************************
!> \brief Default parameter sets for semi empirical models: PM3
!> \param sep ...
!> \param z ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pm3_default_parameter ( sep, z, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pm3_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ    = do_method_pm3
    sep%name   = ptable(z)%name
    sep%eheat  = ptable(z)%heat_of_formation / kcalmol
    sep%z      = z

    SELECT CASE (z)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            "invalid atomic number")
    CASE (0)
       sep%defined = .TRUE.
    CASE (1)
       !Hydrogen
       sep%defined = .TRUE.
       sep%uss     =  -13.0733210_dp
       sep%beta(0) =   -5.6265120_dp
       sep%sto_exponents(0) =    0.9678070_dp
       sep%alp   =    3.3563860_dp
       sep%gss   =   14.7942080_dp
       sep%fn1(1) =  1.1287500_dp; sep%fn2(1) = 5.0962820_dp; sep%fn3(1) = 1.5374650_dp
       sep%fn1(2) = -1.0603290_dp; sep%fn2(2) = 6.0037880_dp; sep%fn3(2) = 1.5701890_dp
       sep%ass=1.0e-10_dp
       sep%asp=3.13469_dp
       sep%app=0.604262_dp
    CASE (2)
    CASE (3)
       !Lithium
       sep%defined = .TRUE.
       sep%uss     =   -5.3000000_dp
       sep%upp     =   -3.4000000_dp
       sep%beta(0) =   -0.5500000_dp
       sep%beta(1) =   -1.5000000_dp
       sep%sto_exponents(0) =    0.6500000_dp
       sep%sto_exponents(1) =    0.7500000_dp
       sep%alp   =    1.2550000_dp
       sep%gss   =    4.5000000_dp
       sep%gsp   =    3.0000000_dp
       sep%gpp   =    5.2500000_dp
       sep%gp2   =    4.5000000_dp
       sep%hsp   =    0.1500000_dp
       sep%fn1(1) = -0.4500000_dp; sep%fn2(1) = 5.0000000_dp; sep%fn3(1) = 1.0000000_dp
       sep%fn1(2) =  0.8000000_dp; sep%fn2(2) = 6.5000000_dp; sep%fn3(2) = 1.0000000_dp
    CASE (4)
       !Beryllium
       sep%defined = .TRUE.
       sep%uss     =  -17.2647520_dp
       sep%upp     =  -11.3042430_dp
       sep%beta(0) =   -3.9620530_dp
       sep%beta(1) =   -2.7806840_dp
       sep%sto_exponents(0) =    0.8774390_dp
       sep%sto_exponents(1) =    1.5087550_dp
       sep%alp   =    1.5935360_dp
       sep%gss   =    9.0128510_dp
       sep%gsp   =    6.5761990_dp
       sep%gpp   =    6.0571820_dp
       sep%gp2   =    9.0052190_dp
       sep%hsp   =    0.5446790_dp
       sep%fn1(1) =  1.6315720_dp; sep%fn2(1) = 2.6729620_dp; sep%fn3(1) = 1.7916860_dp
       sep%fn1(2) = -2.1109590_dp; sep%fn2(2) = 1.9685940_dp; sep%fn3(2) = 1.7558710_dp
    CASE (5)
       !Boron
       sep%defined = .TRUE.
       sep%uss     =  -28.3539160_dp
       sep%upp     =  -23.9738780_dp
       sep%beta(0) =   -5.5215820_dp
       sep%beta(1) =   -5.9611770_dp
       sep%sto_exponents(0) =    0.9885410_dp
       sep%sto_exponents(1) =    1.9094410_dp
       sep%alp   =    2.4981070_dp
       sep%gss   =    8.8104530_dp
       sep%gsp   =   10.2992280_dp
       sep%gpp   =    4.9891680_dp
       sep%gp2   =    7.2899770_dp
       sep%hsp   =    0.8370320_dp
       sep%fn1(1) = 3.0046640_dp; sep%fn2(1) = 6.0288690_dp; sep%fn3(1) = 0.3472230_dp
       sep%fn1(2) = 0.0201950_dp; sep%fn2(2) = 6.0106580_dp; sep%fn3(2) = 1.2140670_dp
    CASE (6)
       !Carbon
       sep%defined = .TRUE.
       sep%uss     =  -47.2703200_dp
       sep%upp     =  -36.2669180_dp
       sep%beta(0) =  -11.9100150_dp
       sep%beta(1) =   -9.8027550_dp
       sep%sto_exponents(0) =    1.5650850_dp
       sep%sto_exponents(1) =    1.8423450_dp
       sep%alp   =    2.7078070_dp
       sep%gss   =   11.2007080_dp
       sep%gsp   =   10.2650270_dp
       sep%gpp   =   10.7962920_dp
       sep%gp2   =    9.0425660_dp
       sep%hsp   =    2.2909800_dp
       sep%fn1(1) = 0.0501070_dp; sep%fn2(1) = 6.0031650_dp; sep%fn3(1) = 1.6422140_dp
       sep%fn1(2) = 0.0507330_dp; sep%fn2(2) = 6.0029790_dp; sep%fn3(2) = 0.8924880_dp
    CASE (7)
       !Nitrogen
       sep%defined = .TRUE.
       sep%uss     =  -49.3356720_dp
       sep%upp     =  -47.5097360_dp
       sep%beta(0) =  -14.0625210_dp
       sep%beta(1) =  -20.0438480_dp
       sep%sto_exponents(0) =    2.0280940_dp
       sep%sto_exponents(1) =    2.3137280_dp
       sep%alp   =    2.8305450_dp
       sep%gss   =   11.9047870_dp
       sep%gsp   =    7.3485650_dp
       sep%gpp   =   11.7546720_dp
       sep%gp2   =   10.8072770_dp
       sep%hsp   =    1.1367130_dp
       sep%fn1(1) =  1.5016740_dp; sep%fn2(1) = 5.9011480_dp; sep%fn3(1) = 1.7107400_dp
       sep%fn1(2) = -1.5057720_dp; sep%fn2(2) = 6.0046580_dp; sep%fn3(2) = 1.7161490_dp
    CASE (8)
       !Oxygen
       sep%defined = .TRUE.
       sep%uss     =  -86.9930020_dp
       sep%upp     =  -71.8795800_dp
       sep%beta(0) =  -45.2026510_dp
       sep%beta(1) =  -24.7525150_dp
       sep%sto_exponents(0) =    3.7965440_dp
       sep%sto_exponents(1) =    2.3894020_dp
       sep%alp   =    3.2171020_dp
       sep%gss   =   15.7557600_dp
       sep%gsp   =   10.6211600_dp
       sep%gpp   =   13.6540160_dp
       sep%gp2   =   12.4060950_dp
       sep%hsp   =    0.5938830_dp
       sep%fn1(1) = -1.1311280_dp; sep%fn2(1) = 6.0024770_dp; sep%fn3(1) = 1.6073110_dp
       sep%fn1(2) =  1.1378910_dp; sep%fn2(2) = 5.9505120_dp; sep%fn3(2) = 1.5983950_dp
       sep%ass=1.0e-10_dp
       sep%asp=6.906169_dp
       sep%app=0.911779_dp
    CASE (9)
       !Fluorine
       sep%defined = .TRUE.
       sep%uss     = -110.4353030_dp
       sep%upp     = -105.6850470_dp
       sep%beta(0) =  -48.4059390_dp
       sep%beta(1) =  -27.7446600_dp
       sep%sto_exponents(0) =    4.7085550_dp
       sep%sto_exponents(1) =    2.4911780_dp
       sep%alp   =    3.3589210_dp
       sep%gss   =   10.4966670_dp
       sep%gsp   =   16.0736890_dp
       sep%gpp   =   14.8172560_dp
       sep%gp2   =   14.4183930_dp
       sep%hsp   =    0.7277630_dp
       sep%fn1(1) = -0.0121660_dp; sep%fn2(1) = 6.0235740_dp; sep%fn3(1) = 1.8568590_dp
       sep%fn1(2) = -0.0028520_dp; sep%fn2(2) = 6.0037170_dp; sep%fn3(2) = 2.6361580_dp
    CASE (10)
    CASE (11)
       !Sodium
       sep%defined = .TRUE.
       sep%uss     =  -4.7683450_dp
       sep%upp     =  -4.7703090_dp
       sep%beta(0) =  -0.5354170_dp
       sep%beta(1) =   0.0548760_dp
       sep%sto_exponents(0) =    0.6767940_dp
       sep%sto_exponents(1) =    0.9876530_dp
       sep%alp   =    1.6810000_dp
       sep%gss   =    5.1531420_dp
       sep%gsp   =    3.1778570_dp
       sep%gpp   =    3.7335330_dp
       sep%gp2   =   11.5386740_dp
       sep%hsp   =    0.5455510_dp
       sep%fn1(1) =  0.4175960_dp; sep%fn2(1) = 8.5969780_dp; sep%fn3(1) = 1.6731780_dp
       sep%fn1(2) =  0.4432600_dp; sep%fn2(2) = 3.3594960_dp; sep%fn3(2) = 2.2315110_dp
    CASE (12)
       !Magnesium
       sep%defined = .TRUE.
       sep%uss     =  -14.6236880_dp
       sep%upp     =  -14.1734600_dp
       sep%beta(0) =   -2.0716910_dp
       sep%beta(1) =   -0.5695810_dp
       sep%sto_exponents(0) =    0.6985520_dp
       sep%sto_exponents(1) =    1.4834530_dp
       sep%alp   =    1.3291470_dp
       sep%gss   =    6.6943000_dp
       sep%gsp   =    6.7939950_dp
       sep%gpp   =    6.9104460_dp
       sep%gp2   =    7.0908230_dp
       sep%hsp   =    0.5433000_dp
       sep%fn1(1) =  2.1170500_dp; sep%fn2(1) = 6.0094770_dp; sep%fn3(1) = 2.0844060_dp
       sep%fn1(2) = -2.5477670_dp; sep%fn2(2) = 4.3953700_dp; sep%fn3(2) = 2.0636740_dp
    CASE (13)
       !Aluminium
       sep%defined = .TRUE.
       sep%uss     =  -24.8454040_dp
       sep%upp     =  -22.2641590_dp
       sep%beta(0) =   -0.5943010_dp
       sep%beta(1) =   -0.9565500_dp
       sep%sto_exponents(0) =    1.7028880_dp
       sep%sto_exponents(1) =    1.0736290_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    1.5217030_dp
       sep%gss   =    5.7767370_dp
       sep%gsp   =   11.6598560_dp
       sep%gpp   =    6.3477900_dp
       sep%gp2   =    6.1210770_dp
       sep%hsp   =    4.0062450_dp
       sep%fn1(1) = -0.4730900_dp; sep%fn2(1) = 1.9158250_dp; sep%fn3(1) = 1.4517280_dp
       sep%fn1(2) = -0.1540510_dp; sep%fn2(2) = 6.0050860_dp; sep%fn3(2) = 2.5199970_dp
    CASE (14)
       !Silicon
       sep%defined = .TRUE.
       sep%uss     =  -26.7634830_dp
       sep%upp     =  -22.8136350_dp
       sep%beta(0) =   -2.8621450_dp
       sep%beta(1) =   -3.9331480_dp
       sep%sto_exponents(0) =    1.6350750_dp
       sep%sto_exponents(1) =    1.3130880_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    2.1358090_dp
       sep%gss   =    5.0471960_dp
       sep%gsp   =    5.9490570_dp
       sep%gpp   =    6.7593670_dp
       sep%gp2   =    5.1612970_dp
       sep%hsp   =    0.9198320_dp
       sep%fn1(1) = -0.3906000_dp; sep%fn2(1) = 6.0000540_dp; sep%fn3(1) = 0.6322620_dp
       sep%fn1(2) =  0.0572590_dp; sep%fn2(2) = 6.0071830_dp; sep%fn3(2) = 2.0199870_dp
    CASE (15)
       !Phosphorus
       sep%defined = .TRUE.
       sep%uss     =  -40.4130960_dp
       sep%upp     =  -29.5930520_dp
       sep%beta(0) =  -12.6158790_dp
       sep%beta(1) =   -4.1600400_dp
       sep%sto_exponents(0) =    2.0175630_dp
       sep%sto_exponents(1) =    1.5047320_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    1.9405340_dp
       sep%gss   =    7.8016150_dp
       sep%gsp   =    5.1869490_dp
       sep%gpp   =    6.6184780_dp
       sep%gp2   =    6.0620020_dp
       sep%hsp   =    1.5428090_dp
       sep%fn1(1) = -0.6114210_dp; sep%fn2(1) = 1.9972720_dp; sep%fn3(1) = 0.7946240_dp
       sep%fn1(2) = -0.0939350_dp; sep%fn2(2) = 1.9983600_dp; sep%fn3(2) = 1.9106770_dp
    CASE (16)
       !Sulphur
       sep%defined = .TRUE.
       sep%uss     =  -49.8953710_dp
       sep%upp     =  -44.3925830_dp
       sep%beta(0) =   -8.8274650_dp
       sep%beta(1) =   -8.0914150_dp
       sep%sto_exponents(0) =    1.8911850_dp
       sep%sto_exponents(1) =    1.6589720_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    2.2697060_dp
       sep%gss   =    8.9646670_dp
       sep%gsp   =    6.7859360_dp
       sep%gpp   =    9.9681640_dp
       sep%gp2   =    7.9702470_dp
       sep%hsp   =    4.0418360_dp
       sep%fn1(1) = -0.3991910_dp; sep%fn2(1) = 6.0006690_dp; sep%fn3(1) = 0.9621230_dp
       sep%fn1(2) = -0.0548990_dp; sep%fn2(2) = 6.0018450_dp; sep%fn3(2) = 1.5799440_dp
    CASE (17)
       !Chlorine
       sep%defined = .TRUE.
       sep%uss     = -100.6267470_dp
       sep%upp     =  -53.6143960_dp
       sep%beta(0) =  -27.5285600_dp
       sep%beta(1) =  -11.5939220_dp
       sep%sto_exponents(0) =    2.2462100_dp
       sep%sto_exponents(1) =    2.1510100_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    2.5172960_dp
       sep%gss   =   16.0136010_dp
       sep%gsp   =    8.0481150_dp
       sep%gpp   =    7.5222150_dp
       sep%gp2   =    7.5041540_dp
       sep%hsp   =    3.4811530_dp
       sep%fn1(1) = -0.1715910_dp; sep%fn2(1) = 6.0008020_dp; sep%fn3(1) = 1.0875020_dp
       sep%fn1(2) = -0.0134580_dp; sep%fn2(2) = 1.9666180_dp; sep%fn3(2) = 2.2928910_dp
    CASE (18)
    CASE (19)
       !Potassium
       sep%defined = .TRUE.
       sep%uss     = -4.1079950_dp
       sep%upp     = -1.9126360_dp
       sep%beta(0) =  -0.7235020_dp
       sep%beta(1) =   0.0109380_dp
       sep%sto_exponents(0) =    0.9158730_dp
       sep%sto_exponents(1) =    1.3430990_dp
       sep%alp   =    4.4452360_dp
       sep%gss   =    4.0502660_dp
       sep%gsp   =    2.9879670_dp
       sep%gpp   =    2.9983050_dp
       sep%gp2   =    4.9936140_dp
       sep%hsp   =    0.5653070_dp
       sep%fn1(1) = 6.0409170_dp; sep%fn2(1) = 10.0465540_dp; sep%fn3(1) = 1.6005120_dp
       sep%fn1(2) = 0.7242240_dp; sep%fn2(2) =  3.0082750_dp; sep%fn3(2) = 2.4643520_dp
    CASE (20)
       !Calcium
       sep%defined = .TRUE.
       sep%uss     = -11.6466770_dp
       sep%upp     = -11.7003040_dp
       sep%beta(0) =  -1.3663020_dp
       sep%beta(1) =  -0.6800840_dp
       sep%sto_exponents(0) =    0.8224530_dp
       sep%sto_exponents(1) =    1.2791170_dp
       sep%alp   =    3.2251430_dp
       sep%gss   =    5.5661490_dp
       sep%gsp   =   15.8559960_dp
       sep%gpp   =    6.7057840_dp
       sep%gp2   =    4.8642700_dp
       sep%hsp   =    1.0053200_dp
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
       sep%uss     = -17.5170270_dp
       sep%upp     = -12.5337290_dp
       sep%udd     = -44.1249280_dp
       sep%beta(0) =  -0.1000000_dp
       sep%beta(1) =  -0.1000000_dp
       sep%beta(2) =  -8.7766360_dp
       sep%sto_exponents(0) =    1.5000000_dp
       sep%sto_exponents(1) =    1.5000000_dp
       sep%sto_exponents(1) =    2.8845490_dp
       sep%alp   =    3.0683070_dp
       sep%gss   =    6.0000000_dp
       sep%gsp   =    4.1500000_dp
       sep%gpp   =    5.0000000_dp
       sep%gp2   =    3.5000000_dp
       sep%hsp   =    1.0000000_dp
    CASE (25:29)
    CASE (30)
       !Zinc
       sep%defined = .TRUE.
       sep%uss     =  -18.5321980_dp
       sep%upp     =  -11.0474090_dp
       sep%beta(0) =   -0.7155780_dp
       sep%beta(1) =   -6.3518640_dp
       sep%sto_exponents(0) =    1.8199890_dp
       sep%sto_exponents(1) =    1.5069220_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    1.3501260_dp
       sep%gss   =    9.6771960_dp
       sep%gsp   =    7.7362040_dp
       sep%gpp   =    4.9801740_dp
       sep%gp2   =    4.6696560_dp
       sep%hsp   =    0.6004130_dp
       sep%fn1(1) = -0.1112340_dp; sep%fn2(1) = 6.0014780_dp; sep%fn3(1) = 1.5160320_dp
       sep%fn1(2) = -0.1323700_dp; sep%fn2(2) = 1.9958390_dp; sep%fn3(2) = 2.5196420_dp
    CASE (31)
       !Gallium
       sep%defined = .TRUE.
       sep%uss     =  -29.8555930_dp
       sep%upp     =  -21.8753710_dp
       sep%beta(0) =   -4.9456180_dp
       sep%beta(1) =   -0.4070530_dp
       sep%sto_exponents(0) =    1.8470400_dp
       sep%sto_exponents(1) =    0.8394110_dp
       sep%alp   =    1.6051150_dp
       sep%gss   =    8.4585540_dp
       sep%gsp   =    8.9256190_dp
       sep%gpp   =    5.0868550_dp
       sep%gp2   =    4.9830450_dp
       sep%hsp   =    2.0512600_dp
       sep%fn1(1) = -0.5601790_dp; sep%fn2(1) = 5.6232730_dp; sep%fn3(1) = 1.5317800_dp
       sep%fn1(2) = -0.2727310_dp; sep%fn2(2) = 1.9918430_dp; sep%fn3(2) = 2.1838640_dp
    CASE (32)
       !Germanium
       sep%defined = .TRUE.
       sep%uss     =  -35.4671955_dp
       sep%upp     =  -31.5863583_dp
       sep%beta(0) =   -5.3250024_dp
       sep%beta(1) =   -2.2501567_dp
       sep%sto_exponents(0) =    2.2373526_dp
       sep%sto_exponents(1) =    1.5924319_dp
       sep%alp   =    1.9723370_dp
       sep%gss   =    5.3769635_dp
       sep%gsp   =   10.2095293_dp
       sep%gpp   =    7.6718647_dp
       sep%gp2   =    6.9242663_dp
       sep%hsp   =    1.3370204_dp
       sep%fn1(1) =  0.9631726_dp; sep%fn2(1) = 6.0120134_dp; sep%fn3(1) = 2.1633655_dp
       sep%fn1(2) = -0.9593891_dp; sep%fn2(2) = 5.7491802_dp; sep%fn3(2) = 2.1693724_dp
    CASE (33)
       !Arsenic
       sep%defined = .TRUE.
       sep%uss     =  -38.5074240_dp
       sep%upp     =  -35.1524150_dp
       sep%beta(0) =   -8.2321650_dp
       sep%beta(1) =   -5.0173860_dp
       sep%sto_exponents(0) =    2.6361770_dp
       sep%sto_exponents(1) =    1.7038890_dp
       sep%alp   =    1.7944770_dp
       sep%gss   =    8.7890010_dp
       sep%gsp   =    5.3979830_dp
       sep%gpp   =    8.2872500_dp
       sep%gp2   =    8.2103460_dp
       sep%hsp   =    1.9510340_dp
       sep%fn1(1) = -0.4600950_dp; sep%fn2(1) = 1.9831150_dp; sep%fn3(1) = 1.0867930_dp
       sep%fn1(2) = -0.0889960_dp; sep%fn2(2) = 1.9929440_dp; sep%fn3(2) = 2.1400580_dp
    CASE (34)
       !Selenium
       sep%defined = .TRUE.
       sep%uss     =  -55.3781350_dp
       sep%upp     =  -49.8230760_dp
       sep%beta(0) =   -6.1578220_dp
       sep%beta(1) =   -5.4930390_dp
       sep%sto_exponents(0) =    2.8280510_dp
       sep%sto_exponents(1) =    1.7325360_dp
       sep%alp   =    3.0439570_dp
       sep%gss   =    7.4325910_dp
       sep%gsp   =   10.0604610_dp
       sep%gpp   =    9.5683260_dp
       sep%gp2   =    7.7242890_dp
       sep%hsp   =    4.0165580_dp
       sep%fn1(1) = 0.0478730_dp; sep%fn2(1) = 6.0074000_dp; sep%fn3(1) = 2.0817170_dp
       sep%fn1(2) = 0.1147200_dp; sep%fn2(2) = 6.0086720_dp; sep%fn3(2) = 1.5164230_dp
    CASE (35)
       !Bromine
       sep%defined = .TRUE.
       sep%uss     = -116.6193110_dp
       sep%upp     =  -74.2271290_dp
       sep%beta(0) =  -31.1713420_dp
       sep%beta(1) =   -6.8140130_dp
       sep%sto_exponents(0) =    5.3484570_dp
       sep%sto_exponents(1) =    2.1275900_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    2.5118420_dp
       sep%gss   =   15.9434250_dp
       sep%gsp   =   16.0616800_dp
       sep%gpp   =    8.2827630_dp
       sep%gp2   =    7.8168490_dp
       sep%hsp   =    0.5788690_dp
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
       sep%uss     =  -15.8285840_dp
       sep%upp     =    8.7497950_dp
       sep%beta(0) =   -8.5819440_dp
       sep%beta(1) =   -0.6010340_dp
       sep%sto_exponents(0) =    1.6793510_dp
       sep%sto_exponents(1) =    2.0664120_dp
       sep%alp  =    1.5253820_dp
       sep%gss  =    9.2069600_dp
       sep%gsp  =    8.2315390_dp
       sep%gpp  =    4.9481040_dp
       sep%gp2  =    4.6696560_dp
       sep%hsp  =    1.6562340_dp
    CASE (49)
       !Indium
       sep%defined = .TRUE.
       sep%uss     =  -26.1762050_dp
       sep%upp     =  -20.0058220_dp
       sep%beta(0) =   -2.9933190_dp
       sep%beta(1) =   -1.8289080_dp
       sep%sto_exponents(0) =    2.0161160_dp
       sep%sto_exponents(1) =    1.4453500_dp
       sep%alp   =    1.4183850_dp
       sep%gss   =    6.5549000_dp
       sep%gsp   =    8.2298730_dp
       sep%gpp   =    6.2992690_dp
       sep%gp2   =    4.9842110_dp
       sep%hsp   =    2.6314610_dp
       sep%fn1(1) = -0.3431380_dp; sep%fn2(1) = 1.9940340_dp; sep%fn3(1) = 1.6255160_dp
       sep%fn1(2) = -0.1095320_dp; sep%fn2(2) = 5.6832170_dp; sep%fn3(2) = 2.8670090_dp
    CASE (50)
       !Tin
       sep%defined = .TRUE.
       sep%uss     =  -34.5501920_dp
       sep%upp     =  -25.8944190_dp
       sep%beta(0) =   -2.7858020_dp
       sep%beta(1) =   -2.0059990_dp
       sep%sto_exponents(0) =    2.3733280_dp
       sep%sto_exponents(1) =    1.6382330_dp
       sep%alp   =    1.6996500_dp
       sep%gss   =   10.1900330_dp
       sep%gsp   =    7.2353270_dp
       sep%gpp   =    5.6738100_dp
       sep%gp2   =    5.1822140_dp
       sep%hsp   =    1.0331570_dp
       sep%fn1(1) = -0.1503530_dp; sep%fn2(1) = 6.0056940_dp; sep%fn3(1) = 1.7046420_dp
       sep%fn1(2) = -0.0444170_dp; sep%fn2(2) = 2.2573810_dp; sep%fn3(2) = 2.4698690_dp
    CASE (51)
       !Antimony
       sep%defined = .TRUE.
       sep%uss     =  -56.4321960_dp
       sep%upp     =  -29.4349540_dp
       sep%beta(0) =  -14.7942170_dp
       sep%beta(1) =   -2.8179480_dp
       sep%sto_exponents(0) =    2.3430390_dp
       sep%sto_exponents(1) =    1.8999920_dp
       sep%alp   =    2.0343010_dp
       sep%gss   =    9.2382770_dp
       sep%gsp   =    5.2776800_dp
       sep%gpp   =    6.3500000_dp
       sep%gp2   =    6.2500000_dp
       sep%hsp   =    2.4244640_dp
       sep%fn1(1) =  3.0020280_dp; sep%fn2(1) = 6.0053420_dp; sep%fn3(1) = 0.8530600_dp
       sep%fn1(2) = -0.0188920_dp; sep%fn2(2) = 6.0114780_dp; sep%fn3(2) = 2.7933110_dp
    CASE (52)
       !Tellurium
       sep%defined = .TRUE.
       sep%uss     =  -44.9380360_dp
       sep%upp     =  -46.3140990_dp
       sep%beta(0) =   -2.6651460_dp
       sep%beta(1) =   -3.8954300_dp
       sep%sto_exponents(0) =    4.1654920_dp
       sep%sto_exponents(1) =    1.6475550_dp
       sep%alp   =    2.4850190_dp
       sep%gss   =   10.2550730_dp
       sep%gsp   =    8.1691450_dp
       sep%gpp   =    7.7775920_dp
       sep%gp2   =    7.7551210_dp
       sep%hsp   =    3.7724620_dp
       sep%fn1(1) =  0.0333910_dp; sep%fn2(1) = 5.9563790_dp; sep%fn3(1) = 2.2775750_dp
       sep%fn1(2) = -1.9218670_dp; sep%fn2(2) = 4.9732190_dp; sep%fn3(2) = 0.5242430_dp
    CASE (53)
       !Iodine
       sep%defined = .TRUE.
       sep%uss     =  -96.4540370_dp
       sep%upp     =  -61.0915820_dp
       sep%beta(0) =  -14.4942340_dp
       sep%beta(1) =   -5.8947030_dp
       sep%sto_exponents(0) =    7.0010130_dp
       sep%sto_exponents(1) =    2.4543540_dp
       sep%sto_exponents(2) =    1.0000000_dp
       sep%alp   =    1.9901850_dp
       sep%gss   =   13.6319430_dp
       sep%gsp   =   14.9904060_dp
       sep%gpp   =    7.2883300_dp
       sep%gp2   =    5.9664070_dp
       sep%hsp   =    2.6300350_dp
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
       sep%uss     =  -17.7622290_dp
       sep%upp     =  -18.3307510_dp
       sep%beta(0) =   -3.1013650_dp
       sep%beta(1) =   -3.4640310_dp
       sep%sto_exponents(0) =    1.4768850_dp
       sep%sto_exponents(1) =    2.4799510_dp
       sep%alp   =    1.5293770_dp
       sep%gss   =    6.6247200_dp
       sep%gsp   =   10.6392970_dp
       sep%gpp   =   14.7092830_dp
       sep%gp2   =   16.0007400_dp
       sep%hsp   =    2.0363110_dp
       sep%fn1(1) =  1.0827200_dp; sep%fn2(1) = 6.4965980_dp; sep%fn3(1) = 1.1951460_dp
       sep%fn1(2) = -0.0965530_dp; sep%fn2(2) = 3.9262810_dp; sep%fn3(2) = 2.6271600_dp
    CASE (81)
       !Thallium
       sep%defined = .TRUE.
       sep%uss     =  -30.0531700_dp
       sep%upp     =  -26.9206370_dp
       sep%beta(0) =   -1.0844950_dp
       sep%beta(1) =   -7.9467990_dp
       sep%sto_exponents(0) =    6.8679210_dp
       sep%sto_exponents(1) =    1.9694450_dp
       sep%alp   =    1.3409510_dp
       sep%gss   =   10.4604120_dp
       sep%gsp   =   11.2238830_dp
       sep%gpp   =    4.9927850_dp
       sep%gp2   =    8.9627270_dp
       sep%hsp   =    2.5304060_dp
       sep%fn1(1) = -1.3613990_dp; sep%fn2(1) = 3.5572260_dp; sep%fn3(1) = 1.0928020_dp
       sep%fn1(2) = -0.0454010_dp; sep%fn2(2) = 2.3069950_dp; sep%fn3(2) = 2.9650290_dp
    CASE (82)
       !Lead
       sep%defined = .TRUE.
       sep%uss     =  -30.3227560_dp
       sep%upp     =  -24.4258340_dp
       sep%beta(0) =   -6.1260240_dp
       sep%beta(1) =   -1.3954300_dp
       sep%sto_exponents(0) =    3.1412890_dp
       sep%sto_exponents(1) =    1.8924180_dp
       sep%alp   =    1.6200450_dp
       sep%gss   =    7.0119920_dp
       sep%gsp   =    6.7937820_dp
       sep%gpp   =    5.1837800_dp
       sep%gp2   =    5.0456510_dp
       sep%hsp   =    1.5663020_dp
       sep%fn1(1) = -0.1225760_dp; sep%fn2(1) = 6.0030620_dp; sep%fn3(1) = 1.9015970_dp
       sep%fn1(2) = -0.0566480_dp; sep%fn2(2) = 4.7437050_dp; sep%fn3(2) = 2.8618790_dp
    CASE (83)
       !Bismuth
       sep%defined = .TRUE.
       sep%uss     =  -33.4959380_dp
       sep%upp     =  -35.5210260_dp
       sep%beta(0) =   -5.6072830_dp
       sep%beta(1) =   -5.8001520_dp
       sep%sto_exponents(0) =    4.9164510_dp
       sep%sto_exponents(1) =    1.9349350_dp
       sep%alp   =    1.8574310_dp
       sep%gss   =    4.9894800_dp
       sep%gsp   =    6.1033080_dp
       sep%gpp   =    8.6960070_dp
       sep%gp2   =    8.3354470_dp
       sep%hsp   =    0.5991220_dp
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
    CASE (92:nelem)
    END SELECT
    ! Finalize parameters for the element
    CALL valence_electrons(sep, extended_basis_set=.FALSE., error=error)
    CALL calpar(z,sep,error)
    CALL convert_param_to_cp2k(sep,error)

  END SUBROUTINE pm3_default_parameter

! *****************************************************************************
!> \brief Default parameter sets for semi empirical models: PM6
!> \param sep ...
!> \param z ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pm6_default_parameter ( sep, z, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pm6_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ    = do_method_pm6
    sep%name   = ptable(z)%name
    sep%eheat  = ptable(z)%heat_of_formation / kcalmol
    sep%z      = z

    SELECT CASE (z)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            "invalid atomic number")
    CASE (0)
       sep%defined = .TRUE.
    CASE(1)
       ! Hydrogen
       sep%defined = .TRUE.
       sep%uss         =  -11.246958_dp
       sep%beta(0)     =   -8.352984_dp
       sep%sto_exponents(0) =    1.268641_dp
       sep%gss         =   14.448686_dp
       sep%a           =    0.024184_dp
       sep%b           =    3.055953_dp
       sep%c           =    1.786011_dp
    CASE(2)
       ! Helium
       sep%defined = .TRUE.
       sep%uss         =  -31.770969_dp
       sep%upp         =   -5.856352_dp
       sep%beta(0)     =  -58.903774_dp
       sep%beta(1)     =  -37.039974_dp
       sep%sto_exponents(0) =    3.313204_dp
       sep%sto_exponents(1) =    3.657133_dp
       sep%gss         =    9.445299_dp
       sep%gsp         =   11.201419_dp
       sep%gpp         =    9.214548_dp
       sep%gp2         =   13.046115_dp
       sep%hsp         =    0.299954_dp
    CASE(3)
       ! Lithium
       sep%defined = .TRUE.
       sep%uss         =   -4.709912_dp
       sep%upp         =   -2.722581_dp
       sep%beta(0)     =   -2.283946_dp
       sep%beta(1)     =   -7.535573_dp
       sep%sto_exponents(0) =    0.981041_dp
       sep%sto_exponents(1) =    2.953445_dp
       sep%gss         =   11.035907_dp
       sep%gsp         =   19.998647_dp
       sep%gpp         =   11.543650_dp
       sep%gp2         =    9.059036_dp
       sep%hsp         =    1.641886_dp
    CASE(4)
       ! Beryllium
       sep%defined = .TRUE.
       sep%uss         =  -16.360315_dp
       sep%upp         =  -16.339216_dp
       sep%beta(0)     =   -3.199549_dp
       sep%beta(1)     =   -4.451920_dp
       sep%sto_exponents(0) =    1.212539_dp
       sep%sto_exponents(1) =    1.276487_dp
       sep%gss         =    7.552804_dp
       sep%gsp         =   10.203146_dp
       sep%gpp         =   12.862153_dp
       sep%gp2         =   13.602858_dp
       sep%hsp         =    1.501452_dp
       sep%a           =    0.164180_dp
       sep%b           =    1.704828_dp
       sep%c           =    1.785591_dp
    CASE(5)
       ! Boron
       sep%defined = .TRUE.
       sep%uss         =  -25.967679_dp
       sep%upp         =  -19.115864_dp
       sep%beta(0)     =   -4.959706_dp
       sep%beta(1)     =   -4.656753_dp
       sep%sto_exponents(0) =    1.634174_dp
       sep%sto_exponents(1) =    1.479195_dp
       sep%gss         =    8.179341_dp
       sep%gsp         =    7.294021_dp
       sep%gpp         =    7.829395_dp
       sep%gp2         =    6.401072_dp
       sep%hsp         =    1.252845_dp
    CASE(6)
       ! Carbon
       sep%defined = .TRUE.
       sep%uss         =  -51.089653_dp
       sep%upp         =  -39.937920_dp
       sep%beta(0)     =  -15.385236_dp
       sep%beta(1)     =   -7.471929_dp
       sep%sto_exponents(0) =    2.047558_dp
       sep%sto_exponents(1) =    1.702841_dp
       sep%gss         =   13.335519_dp
       sep%gsp         =   11.528134_dp
       sep%gpp         =   10.778326_dp
       sep%gp2         =    9.486212_dp
       sep%hsp         =    0.717322_dp
       sep%a           =    0.046302_dp
       sep%b           =    2.100206_dp
       sep%c           =    1.333959_dp
    CASE(7)
       ! Nitrogen
       sep%defined = .TRUE.
       sep%uss         =  -57.784823_dp
       sep%upp         =  -49.893036_dp
       sep%beta(0)     =  -17.979377_dp
       sep%beta(1)     =  -15.055017_dp
       sep%sto_exponents(0) =    2.380406_dp
       sep%sto_exponents(1) =    1.999246_dp
       sep%gss         =   12.357026_dp
       sep%gsp         =    9.636190_dp
       sep%gpp         =   12.570756_dp
       sep%gp2         =   10.576425_dp
       sep%hsp         =    2.871545_dp
       sep%a           =   -0.001436_dp
       sep%b           =    0.495196_dp
       sep%c           =    1.704857_dp
    CASE(8)
       ! Oxygen
       sep%defined = .TRUE.
       sep%uss         =  -91.678761_dp
       sep%upp         =  -70.460949_dp
       sep%beta(0)     =  -65.635137_dp
       sep%beta(1)     =  -21.622604_dp
       sep%sto_exponents(0) =    5.421751_dp
       sep%sto_exponents(1) =    2.270960_dp
       sep%gss         =   11.304042_dp
       sep%gsp         =   15.807424_dp
       sep%gpp         =   13.618205_dp
       sep%gp2         =   10.332765_dp
       sep%hsp         =    5.010801_dp
       sep%a           =   -0.017771_dp
       sep%b           =    3.058310_dp
       sep%c           =    1.896435_dp
    CASE(9)
       ! Fluorine
       sep%defined = .TRUE.
       sep%uss         = -140.225626_dp
       sep%upp         =  -98.778044_dp
       sep%beta(0)     =  -69.922593_dp
       sep%beta(1)     =  -30.448165_dp
       sep%sto_exponents(0) =    6.043849_dp
       sep%sto_exponents(1) =    2.906722_dp
       sep%gss         =   12.446818_dp
       sep%gsp         =   18.496082_dp
       sep%gpp         =    8.417366_dp
       sep%gp2         =   12.179816_dp
       sep%hsp         =    2.604382_dp
       sep%a           =   -0.010792_dp
       sep%b           =    6.004648_dp
       sep%c           =    1.847724_dp
    CASE(10)
       ! Neon
       sep%defined = .TRUE.
       sep%uss         =   -2.978729_dp
       sep%upp         =  -85.441118_dp
       sep%beta(0)     =  -69.793475_dp
       sep%beta(1)     =  -33.261962_dp
       sep%sto_exponents(0) =    6.000148_dp
       sep%sto_exponents(1) =    3.834528_dp
       sep%gss         =   19.999574_dp
       sep%gsp         =   16.896951_dp
       sep%gpp         =    8.963560_dp
       sep%gp2         =   16.027799_dp
       sep%hsp         =    1.779280_dp
    CASE(11)
       ! Sodium
       sep%defined = .TRUE.
       sep%uss         =   -4.537153_dp
       sep%upp         =   -2.433015_dp
       sep%beta(0)     =    0.244853_dp
       sep%beta(1)     =    0.491998_dp
       sep%sto_exponents(0) =    0.686327_dp
       sep%sto_exponents(1) =    0.950068_dp
       sep%gss         =    4.059972_dp
       sep%gsp         =    7.061183_dp
       sep%gpp         =    9.283540_dp
       sep%gp2         =   17.034978_dp
       sep%hsp         =    0.640715_dp
       sep%a           =   -1.026036_dp
       sep%b           =    2.014506_dp
       sep%c           =    1.271202_dp
    CASE(12)
       ! Magnesium
       sep%defined = .TRUE.
       sep%uss         =  -14.574226_dp
       sep%upp         =   -7.583850_dp
       sep%beta(0)     =   -9.604932_dp
       sep%beta(1)     =    3.416908_dp
       sep%sto_exponents(0) =    1.310830_dp
       sep%sto_exponents(1) =    1.388897_dp
       sep%gss         =    7.115328_dp
       sep%gsp         =    3.253024_dp
       sep%gpp         =    4.737311_dp
       sep%gp2         =    8.428485_dp
       sep%hsp         =    0.877379_dp
    CASE(13)
       ! Aluminium
       sep%defined = .TRUE.
       sep%uss         =  -24.546778_dp
       sep%upp         =  -20.104434_dp
       sep%udd         =    8.004394_dp
       sep%beta(0)     =  -18.375229_dp
       sep%beta(1)     =   -9.382700_dp
       sep%beta(2)     =  -20.840474_dp
       sep%sto_exponents(0) =    2.364264_dp
       sep%sto_exponents(1) =    1.749102_dp
       sep%sto_exponents(2) =    1.269384_dp
       sep%zn(0)       =    4.742341_dp
       sep%zn(1)       =    4.669626_dp
       sep%zn(2)       =    7.131138_dp
       sep%gss         =    6.652155_dp
       sep%gsp         =    7.459435_dp
       sep%gpp         =    7.668857_dp
       sep%gp2         =    6.673299_dp
       sep%hsp         =    0.435060_dp
       sep%a           =    1.002222_dp
       sep%b           =    1.517400_dp
       sep%c           =    0.659101_dp
    CASE(14)
       ! Silicon
       sep%defined = .TRUE.
       sep%uss         =  -27.358058_dp
       sep%upp         =  -20.490578_dp
       sep%udd         =  -22.751900_dp
       sep%beta(0)     =   -8.686909_dp
       sep%beta(1)     =   -1.856482_dp
       sep%beta(2)     =   -6.360627_dp
       sep%sto_exponents(0) =    1.752741_dp
       sep%sto_exponents(1) =    1.198413_dp
       sep%sto_exponents(2) =    2.128593_dp
       sep%zn(0)       =    8.388111_dp
       sep%zn(1)       =    1.843048_dp
       sep%zn(2)       =    0.708600_dp
       sep%gss         =    5.194805_dp
       sep%gsp         =    5.090534_dp
       sep%gpp         =    5.185150_dp
       sep%gp2         =    4.769775_dp
       sep%hsp         =    1.425012_dp
       sep%a           =    0.208571_dp
       sep%b           =    6.000483_dp
       sep%c           =    1.185245_dp
    CASE(15)
       ! Phosphorus
       sep%defined = .TRUE.
       sep%uss         =  -48.729905_dp
       sep%upp         =  -40.354689_dp
       sep%udd         =   -7.349246_dp
       sep%beta(0)     =  -14.583780_dp
       sep%beta(1)     =  -11.744725_dp
       sep%beta(2)     =  -20.099893_dp
       sep%sto_exponents(0) =    2.158033_dp
       sep%sto_exponents(1) =    1.805343_dp
       sep%sto_exponents(2) =    1.230358_dp
       sep%zn(0)       =    6.042706_dp
       sep%zn(1)       =    2.376473_dp
       sep%zn(2)       =    7.147750_dp
       sep%gss         =    8.758856_dp
       sep%gsp         =    8.483679_dp
       sep%gpp         =    8.662754_dp
       sep%gp2         =    7.734264_dp
       sep%hsp         =    0.871681_dp
       sep%a           =   -0.034320_dp
       sep%b           =    6.001394_dp
       sep%c           =    2.296737_dp
    CASE(16)
       ! Sulfur
       sep%defined = .TRUE.
       sep%uss         =  -47.530706_dp
       sep%upp         =  -39.191045_dp
       sep%udd         =  -46.306944_dp
       sep%beta(0)     =  -13.827440_dp
       sep%beta(1)     =   -7.664613_dp
       sep%beta(2)     =   -9.986172_dp
       sep%sto_exponents(0) =    2.192844_dp
       sep%sto_exponents(1) =    1.841078_dp
       sep%sto_exponents(2) =    3.109401_dp
       sep%zn(0)       =    0.479722_dp
       sep%zn(1)       =    1.015507_dp
       sep%zn(2)       =    4.317470_dp
       sep%gss         =    9.170350_dp
       sep%gsp         =    5.944296_dp
       sep%gpp         =    8.165473_dp
       sep%gp2         =    7.301878_dp
       sep%hsp         =    5.005404_dp
       sep%a           =   -0.036928_dp
       sep%b           =    1.795067_dp
       sep%c           =    2.082618_dp
    CASE(17)
       ! Chlorine
       sep%defined = .TRUE.
       sep%uss         =  -61.389930_dp
       sep%upp         =  -54.482801_dp
       sep%udd         =  -38.258155_dp
       sep%beta(0)     =   -2.367988_dp
       sep%beta(1)     =  -13.802139_dp
       sep%beta(2)     =   -4.037751_dp
       sep%sto_exponents(0) =    2.637050_dp
       sep%sto_exponents(1) =    2.118146_dp
       sep%sto_exponents(2) =    1.324033_dp
       sep%zn(0)       =    0.956297_dp
       sep%zn(1)       =    2.464067_dp
       sep%zn(2)       =    6.410325_dp
       sep%gss         =   11.142654_dp
       sep%gsp         =    7.487881_dp
       sep%gpp         =    9.551886_dp
       sep%gp2         =    8.128436_dp
       sep%hsp         =    5.004267_dp
       sep%a           =   -0.013213_dp
       sep%b           =    3.687022_dp
       sep%c           =    2.544635_dp
    CASE(18)
       ! Argon
       sep%defined = .TRUE.
       sep%uss         =   -7.797931_dp
       sep%upp         =  -83.211487_dp
       sep%beta(0)     =   -8.839842_dp
       sep%beta(1)     =  -28.427303_dp
       sep%sto_exponents(0) =    6.000272_dp
       sep%sto_exponents(1) =    5.949170_dp
       sep%gss         =   17.858776_dp
       sep%gsp         =    4.168451_dp
       sep%gpp         =   11.852500_dp
       sep%gp2         =   15.669543_dp
       sep%hsp         =    4.574549_dp
    CASE(19)
       ! Potassium
       sep%defined = .TRUE.
       sep%uss         =   -3.801108_dp
       sep%upp         =   -3.339656_dp
       sep%beta(0)     =   -8.755195_dp
       sep%beta(1)     =   -1.788061_dp
       sep%sto_exponents(0) =    6.000478_dp
       sep%sto_exponents(1) =    1.127503_dp
       sep%gss         =    3.369251_dp
       sep%gsp         =    6.129351_dp
       sep%gpp         =    0.999505_dp
       sep%gp2         =   18.999148_dp
       sep%hsp         =    0.300325_dp
       sep%a           =    0.157519_dp
       sep%b           =    6.000566_dp
       sep%c           =    2.047539_dp
    CASE(20)
       ! Calcium
       sep%defined = .TRUE.
       sep%uss         =  -10.770058_dp
       sep%upp         =   -9.754177_dp
       sep%beta(0)     =   -4.343881_dp
       sep%beta(1)     =   -1.296612_dp
       sep%sto_exponents(0) =    1.528258_dp
       sep%sto_exponents(1) =    2.060094_dp
       sep%gss         =    5.725773_dp
       sep%gsp         =    4.781065_dp
       sep%gpp         =    7.172103_dp
       sep%gp2         =    7.431876_dp
       sep%hsp         =    1.240572_dp
       sep%a           =   -0.025275_dp
       sep%b           =    0.500017_dp
       sep%c           =    2.329051_dp
    CASE(21)
       ! Scandium
       sep%defined = .TRUE.
       sep%uss         =  -15.544461_dp
       sep%upp         =  -18.646295_dp
       sep%udd         =  -16.069444_dp
       sep%beta(0)     =   -8.620944_dp
       sep%beta(1)     =    3.075948_dp
       sep%beta(2)     =   -9.768661_dp
       sep%sto_exponents(0) =    1.402469_dp
       sep%sto_exponents(1) =    1.345196_dp
       sep%sto_exponents(2) =    1.859012_dp
       sep%zn(0)       =    0.848418_dp
       sep%zn(1)       =    2.451729_dp
       sep%zn(2)       =    0.789372_dp
       sep%gss         =    4.638215_dp
       sep%gsp         =    5.739164_dp
       sep%gpp         =   14.604872_dp
       sep%gp2         =   12.802595_dp
       sep%hsp         =    0.193835_dp
       sep%rho         =    3.173734_dp
       sep%f0sd        =    4.798313_dp
       sep%g2sd        =    5.380136_dp
    CASE(22)
       ! Titanium
       sep%defined = .TRUE.
       sep%uss         =  -25.507973_dp
       sep%upp         =  -17.260909_dp
       sep%udd         =  -23.809486_dp
       sep%beta(0)     =    3.389142_dp
       sep%beta(1)     =   -3.355350_dp
       sep%beta(2)     =   -1.842829_dp
       sep%sto_exponents(0) =    5.324777_dp
       sep%sto_exponents(1) =    1.164068_dp
       sep%sto_exponents(2) =    1.418280_dp
       sep%zn(0)       =    1.045904_dp
       sep%zn(1)       =    1.076844_dp
       sep%zn(2)       =    0.717945_dp
       sep%gss         =    5.717851_dp
       sep%gsp         =    5.800015_dp
       sep%gpp         =    6.414726_dp
       sep%gp2         =    5.623133_dp
       sep%hsp         =    1.403732_dp
       sep%f0sd        =    6.560562_dp
       sep%g2sd        =    3.396235_dp
    CASE(23)
       ! Vanadium
       sep%defined = .TRUE.
       sep%uss         =  -32.162276_dp
       sep%upp         =  -21.572501_dp
       sep%udd         =  -34.506245_dp
       sep%beta(0)     =   -1.211330_dp
       sep%beta(1)     =    0.740746_dp
       sep%beta(2)     =    3.153669_dp
       sep%sto_exponents(0) =    1.974330_dp
       sep%sto_exponents(1) =    1.063106_dp
       sep%sto_exponents(2) =    1.394806_dp
       sep%zn(0)       =    1.094426_dp
       sep%zn(1)       =    0.755378_dp
       sep%zn(2)       =    1.099367_dp
       sep%gss         =    5.983116_dp
       sep%gsp         =    4.736769_dp
       sep%gpp         =    4.499763_dp
       sep%gp2         =    3.944481_dp
       sep%hsp         =    0.901105_dp
       sep%f0sd        =    6.810021_dp
       sep%g2sd        =    1.831407_dp
    CASE(24)
       ! Chromium
       sep%defined = .TRUE.
       sep%uss         =  -34.864339_dp
       sep%upp         =  -26.978615_dp
       sep%udd         =  -54.431036_dp
       sep%beta(0)     =   -5.122615_dp
       sep%beta(1)     =    3.926711_dp
       sep%beta(2)     =   -4.230550_dp
       sep%sto_exponents(0) =    3.283460_dp
       sep%sto_exponents(1) =    1.029394_dp
       sep%sto_exponents(2) =    1.623119_dp
       sep%zn(0)       =    1.619853_dp
       sep%zn(1)       =    0.848266_dp
       sep%zn(2)       =    1.405015_dp
       sep%gss         =    8.855572_dp
       sep%gsp         =    5.588631_dp
       sep%gpp         =    5.053094_dp
       sep%gp2         =    4.429530_dp
       sep%hsp         =    0.648039_dp
       sep%f0sd        =    6.150136_dp
       sep%g2sd        =    2.000300_dp
    CASE(25)
       ! Manganese
       sep%defined = .TRUE.
       sep%uss         =  -51.460000_dp
       sep%upp         =  -37.543990_dp
       sep%udd         =  -47.655370_dp
       sep%beta(0)     =   -4.185290_dp
       sep%beta(1)     =   -3.479630_dp
       sep%beta(2)     =  -13.473190_dp
       sep%sto_exponents(0) =    2.131680_dp
       sep%sto_exponents(1) =    1.525880_dp
       sep%sto_exponents(2) =    2.607800_dp
       sep%zn(0)       =    1.132450_dp
       sep%zn(1)       =    1.390740_dp
       sep%zn(2)       =    0.962550_dp
       sep%gss         =    6.190990_dp
       sep%gsp         =    6.757427_dp
       sep%gpp         =    8.284594_dp
       sep%gp2         =    7.262255_dp
       sep%hsp         =    1.520518_dp
       sep%f0sd        =    7.690920_dp
       sep%g2sd        =    1.105330_dp
    CASE(26)
       ! Iron
       sep%defined = .TRUE.
       sep%uss         =  -70.515047_dp
       sep%upp         =  -62.963069_dp
       sep%udd         = -103.631790_dp
       sep%beta(0)     =    8.027621_dp
       sep%beta(1)     =   -1.125760_dp
       sep%beta(2)     =   -3.507531_dp
       sep%sto_exponents(0) =    1.479150_dp
       sep%sto_exponents(1) =    6.002246_dp
       sep%sto_exponents(2) =    1.080747_dp
       sep%zn(0)       =    1.459152_dp
       sep%zn(1)       =    1.392614_dp
       sep%zn(2)       =    2.161909_dp
       sep%gss         =    7.977036_dp
       sep%gsp         =    7.786867_dp
       sep%gpp         =    8.295758_dp
       sep%gp2         =    7.272041_dp
       sep%hsp         =    1.880189_dp
       sep%rho         =    1.272092_dp
       sep%f0sd        =    9.300165_dp
       sep%g2sd        =    1.601345_dp
    CASE(27)
       ! Cobalt
       sep%defined = .TRUE.
       sep%uss         =  -21.039413_dp
       sep%upp         =   10.000000_dp
       sep%udd         =  -28.068971_dp
       sep%beta(0)     =   -8.992062_dp
       sep%beta(1)     =   -0.100000_dp
       sep%beta(2)     =   -2.481509_dp
       sep%sto_exponents(0) =    1.166613_dp
       sep%sto_exponents(1) =    3.000000_dp
       sep%sto_exponents(2) =    1.860218_dp
       sep%zn(0)       =    0.519518_dp
       sep%zn(1)       =    1.000000_dp
       sep%zn(2)       =    0.352115_dp
       sep%gss         =    2.840152_dp
       sep%gsp         =    3.425933_dp
       sep%gpp         =    5.956968_dp
       sep%gp2         =    5.221864_dp
       sep%hsp         =    0.390087_dp
       sep%f0sd        =    1.446283_dp
       sep%g2sd        =    1.680225_dp
    CASE(28)
       ! Nickel
       sep%defined = .TRUE.
       sep%uss         =  -47.620247_dp
       sep%upp         =  -32.878408_dp
       sep%udd         =  -93.026395_dp
       sep%beta(0)     =   -9.151521_dp
       sep%beta(1)     =   -8.086696_dp
       sep%beta(2)     =   -8.655910_dp
       sep%sto_exponents(0) =    1.591828_dp
       sep%sto_exponents(1) =    2.304739_dp
       sep%sto_exponents(2) =    2.514761_dp
       sep%zn(0)       =    0.746470_dp
       sep%zn(1)       =    0.753327_dp
       sep%zn(2)       =    1.461345_dp
       sep%gss         =    4.080876_dp
       sep%gsp         =    4.099452_dp
       sep%gpp         =    4.487545_dp
       sep%gp2         =    3.933771_dp
       sep%hsp         =    0.993498_dp
       sep%rho         =    1.586979_dp
       sep%f0sd        =    4.651664_dp
       sep%g2sd        =    1.880502_dp
    CASE(29)
       ! Copper
       sep%defined = .TRUE.
       sep%uss         =  -97.002205_dp
       sep%upp         =   -1.000000_dp
       sep%udd         = -110.442592_dp
       sep%beta(0)     =   -9.369508_dp
       sep%beta(1)     =   -0.100000_dp
       sep%beta(2)     =  -16.982092_dp
       sep%sto_exponents(0) =    1.669096_dp
       sep%sto_exponents(1) =    3.000000_dp
       sep%sto_exponents(2) =    2.734990_dp
       sep%zn(0)       =    1.899598_dp
       sep%zn(1)       =    3.000000_dp
       sep%zn(2)       =    1.484317_dp
       sep%gss         =   10.384910_dp
       sep%gsp         =   12.145361_dp
       sep%gpp         =   17.870905_dp
       sep%gp2         =   15.665592_dp
       sep%hsp         =    2.037394_dp
       sep%f0sd        =    9.848807_dp
       sep%g2sd        =    9.847577_dp
    CASE(30)
       ! Zinc
       sep%defined = .TRUE.
       sep%uss         =  -18.040862_dp
       sep%upp         =   -7.834895_dp
       sep%beta(0)     =  -13.276583_dp
       sep%beta(1)     =    1.479642_dp
       sep%sto_exponents(0) =    1.512875_dp
       sep%sto_exponents(1) =    1.789482_dp
       sep%gss         =    8.707424_dp
       sep%gsp         =    3.436116_dp
       sep%gpp         =   20.000041_dp
       sep%gp2         =    6.782785_dp
       sep%hsp         =    0.662036_dp
    CASE(31)
       ! Gallium
       sep%defined = .TRUE.
       sep%uss         =  -30.600226_dp
       sep%upp         =  -21.032425_dp
       sep%beta(0)     =  -10.808320_dp
       sep%beta(1)     =   -4.185500_dp
       sep%sto_exponents(0) =    2.339067_dp
       sep%sto_exponents(1) =    1.729592_dp
       sep%gss         =   10.354885_dp
       sep%gsp         =    7.993674_dp
       sep%gpp         =    6.090184_dp
       sep%gp2         =    6.299226_dp
       sep%hsp         =    1.295974_dp
    CASE(32)
       ! Germanium
       sep%defined = .TRUE.
       sep%uss         =  -32.747338_dp
       sep%upp         =  -24.709016_dp
       sep%beta(0)     =  -14.854297_dp
       sep%beta(1)     =   -2.591260_dp
       sep%sto_exponents(0) =    2.546073_dp
       sep%sto_exponents(1) =    1.709130_dp
       sep%gss         =    7.518301_dp
       sep%gsp         =    6.594443_dp
       sep%gpp         =    6.066801_dp
       sep%gp2         =    5.305947_dp
       sep%hsp         =    0.290742_dp
    CASE(33)
       ! Arsenic
       sep%defined = .TRUE.
       sep%uss         =  -37.956965_dp
       sep%upp         =  -38.453701_dp
       sep%udd         =  -30.282658_dp
       sep%beta(0)     =  -11.963725_dp
       sep%beta(1)     =   -7.340073_dp
       sep%beta(2)     =    3.753005_dp
       sep%sto_exponents(0) =    2.926171_dp
       sep%sto_exponents(1) =    1.765191_dp
       sep%sto_exponents(2) =    1.392142_dp
       sep%zn(0)       =    2.006543_dp
       sep%zn(1)       =    3.316832_dp
       sep%zn(2)       =    4.653440_dp
       sep%gss         =    6.665030_dp
       sep%gsp         =    6.213867_dp
       sep%gpp         =    9.310836_dp
       sep%gp2         =    8.712542_dp
       sep%hsp         =    0.280662_dp
    CASE(34)
       ! Selenium
       sep%defined = .TRUE.
       sep%uss         =  -32.671088_dp
       sep%upp         =  -32.522220_dp
       sep%beta(0)     =    2.636001_dp
       sep%beta(1)     =   -9.557700_dp
       sep%sto_exponents(0) =    2.512366_dp
       sep%sto_exponents(1) =    2.007576_dp
       sep%gss         =    5.522356_dp
       sep%gsp         =    2.907562_dp
       sep%gpp         =    8.042391_dp
       sep%gp2         =    6.735106_dp
       sep%hsp         =    3.095789_dp
    CASE(35)
       ! Bromine
       sep%defined = .TRUE.
       sep%uss         =  -45.834364_dp
       sep%upp         =  -50.293675_dp
       sep%udd         =    7.086738_dp
       sep%beta(0)     =  -32.131665_dp
       sep%beta(1)     =   -9.514484_dp
       sep%beta(2)     =   -9.839124_dp
       sep%sto_exponents(0) =    4.670684_dp
       sep%sto_exponents(1) =    2.035626_dp
       sep%sto_exponents(2) =    1.521031_dp
       sep%zn(0)       =    3.094777_dp
       sep%zn(1)       =    3.065764_dp
       sep%zn(2)       =    2.820003_dp
       sep%gss         =    7.616791_dp
       sep%gsp         =    5.010425_dp
       sep%gpp         =    9.649216_dp
       sep%gp2         =    8.343792_dp
       sep%hsp         =    4.996553_dp
       sep%a           =   -0.004996_dp
       sep%b           =    6.001292_dp
       sep%c           =    2.895153_dp
    CASE(36)
       ! Krypton
       sep%defined = .TRUE.
       sep%uss         =    8.535384_dp
       sep%upp         =  -80.484321_dp
       sep%beta(0)     =   -2.727088_dp
       sep%beta(1)     =  -16.142951_dp
       sep%sto_exponents(0) =    1.312248_dp
       sep%sto_exponents(1) =    4.491371_dp
       sep%gss         =   19.999857_dp
       sep%gsp         =    1.175304_dp
       sep%gpp         =    9.174784_dp
       sep%gp2         =   14.926948_dp
       sep%hsp         =    0.299867_dp
    CASE(37)
       ! Rubidium
       sep%defined = .TRUE.
       sep%uss         =   -3.636505_dp
       sep%upp         =   -2.500671_dp
       sep%beta(0)     =    9.998744_dp
       sep%beta(1)     =    1.343004_dp
       sep%sto_exponents(0) =    5.510145_dp
       sep%sto_exponents(1) =    1.335170_dp
       sep%gss         =    6.680824_dp
       sep%gsp         =   20.001098_dp
       sep%gpp         =    5.068874_dp
       sep%gp2         =    2.747860_dp
       sep%hsp         =    3.602834_dp
    CASE(38)
       ! Strontium
       sep%defined = .TRUE.
       sep%uss         =  -10.427671_dp
       sep%upp         =   -9.943751_dp
       sep%beta(0)     =   -6.253108_dp
       sep%beta(1)     =   -9.844498_dp
       sep%sto_exponents(0) =    2.197303_dp
       sep%sto_exponents(1) =    1.730137_dp
       sep%gss         =    4.603664_dp
       sep%gsp         =    5.716069_dp
       sep%gpp         =    7.334620_dp
       sep%gp2         =    7.443088_dp
       sep%hsp         =    0.831527_dp
       sep%a           =   -0.012948_dp
       sep%b           =    6.000126_dp
       sep%c           =    3.011964_dp
    CASE(39)
       ! Yttrium
       sep%defined = .TRUE.
       sep%uss         =  -14.247809_dp
       sep%upp         =  -14.817140_dp
       sep%udd         =  -16.394302_dp
       sep%beta(0)     =    0.343336_dp
       sep%beta(1)     =   -3.180807_dp
       sep%beta(2)     =   -4.508957_dp
       sep%sto_exponents(0) =    0.593368_dp
       sep%sto_exponents(1) =    1.490422_dp
       sep%sto_exponents(2) =    1.650893_dp
       sep%zn(0)       =    0.902611_dp
       sep%zn(1)       =    1.484400_dp
       sep%zn(2)       =    1.384238_dp
       sep%gss         =    4.046733_dp
       sep%gsp         =    4.726277_dp
       sep%gpp         =    7.278752_dp
       sep%gp2         =    6.343281_dp
       sep%hsp         =    0.679228_dp
       sep%rho         =    2.773703_dp
       sep%f0sd        =    4.972716_dp
       sep%g2sd        =    5.016364_dp
    CASE(40)
       ! Zirconium
       sep%defined = .TRUE.
       sep%uss         =  -20.008884_dp
       sep%upp         =  -14.559692_dp
       sep%udd         =  -21.302657_dp
       sep%beta(0)     =    9.551952_dp
       sep%beta(1)     =   -4.551915_dp
       sep%beta(2)     =   -3.213274_dp
       sep%sto_exponents(0) =    1.692590_dp
       sep%sto_exponents(1) =    1.694916_dp
       sep%sto_exponents(2) =    1.567392_dp
       sep%zn(0)       =    1.189109_dp
       sep%zn(1)       =    0.809092_dp
       sep%zn(2)       =    1.190249_dp
       sep%gss         =    5.331208_dp
       sep%gsp         =    4.150579_dp
       sep%gpp         =    3.967381_dp
       sep%gp2         =    3.457490_dp
       sep%hsp         =    0.743676_dp
       sep%f0sd        =    5.010704_dp
       sep%g2sd        =    2.943652_dp
    CASE(41)
       ! Niobium
       sep%defined = .TRUE.
       sep%uss         =  -31.269298_dp
       sep%upp         =  -20.151277_dp
       sep%udd         =  -35.893116_dp
       sep%beta(0)     =  -12.045244_dp
       sep%beta(1)     =    1.465762_dp
       sep%beta(2)     =   -5.920160_dp
       sep%sto_exponents(0) =    2.355562_dp
       sep%sto_exponents(1) =    1.386907_dp
       sep%sto_exponents(2) =    1.977324_dp
       sep%zn(0)       =    1.490754_dp
       sep%zn(1)       =    0.892760_dp
       sep%zn(2)       =    1.443837_dp
       sep%gss         =    6.683592_dp
       sep%gsp         =    4.685339_dp
       sep%gpp         =    4.377647_dp
       sep%gp2         =    3.815028_dp
       sep%hsp         =    0.650679_dp
       sep%f0sd        =    6.550674_dp
       sep%g2sd        =    1.065577_dp
    CASE(42)
       ! Molybdenum
       sep%defined = .TRUE.
       sep%uss         =  -53.467728_dp
       sep%upp         =  -35.291951_dp
       sep%udd         =  -55.836977_dp
       sep%beta(0)     =   -0.189344_dp
       sep%beta(1)     =    7.017762_dp
       sep%beta(2)     =  -10.941126_dp
       sep%sto_exponents(0) =    1.060429_dp
       sep%sto_exponents(1) =    1.350412_dp
       sep%sto_exponents(2) =    1.827152_dp
       sep%zn(0)       =    1.912995_dp
       sep%zn(1)       =    1.355055_dp
       sep%zn(2)       =    1.876231_dp
       sep%gss         =    8.576652_dp
       sep%gsp         =    6.888293_dp
       sep%gpp         =    6.644509_dp
       sep%gp2         =    5.790552_dp
       sep%hsp         =    1.317368_dp
       sep%f0sd        =   10.000608_dp
       sep%g2sd        =    1.216752_dp
    CASE(43)
       ! Technetium
       sep%defined = .TRUE.
       sep%uss         =  -41.850292_dp
       sep%upp         =  -34.910293_dp
       sep%udd         =  -45.530412_dp
       sep%beta(0)     =   -2.791024_dp
       sep%beta(1)     =   -8.086697_dp
       sep%beta(2)     =   -5.724335_dp
       sep%sto_exponents(0) =    1.956245_dp
       sep%sto_exponents(1) =    6.006299_dp
       sep%sto_exponents(2) =    1.767360_dp
       sep%zn(0)       =    1.411033_dp
       sep%zn(1)       =    1.141313_dp
       sep%zn(2)       =    1.159312_dp
       sep%gss         =    6.326174_dp
       sep%gsp         =    5.587138_dp
       sep%gpp         =    5.596426_dp
       sep%gp2         =    4.877169_dp
       sep%hsp         =    1.258989_dp
       sep%f0sd        =    5.434886_dp
       sep%g2sd        =    1.106875_dp
    CASE(44)
       ! Ruthenium
       sep%defined = .TRUE.
       sep%uss         =  -44.901521_dp
       sep%upp         =  -41.424409_dp
       sep%udd         =  -37.934514_dp
       sep%beta(0)     =  -12.859508_dp
       sep%beta(1)     =   -8.475518_dp
       sep%beta(2)     =   -3.830797_dp
       sep%sto_exponents(0) =    1.459195_dp
       sep%sto_exponents(1) =    5.537201_dp
       sep%sto_exponents(2) =    2.093164_dp
       sep%zn(0)       =    0.984449_dp
       sep%zn(1)       =    4.586613_dp
       sep%zn(2)       =    0.765332_dp
       sep%gss         =    4.413643_dp
       sep%gsp         =    5.356996_dp
       sep%gpp         =   22.490448_dp
       sep%gp2         =   19.599957_dp
       sep%hsp         =    0.008058_dp
       sep%f0sd        =    5.917404_dp
       sep%g2sd        =    5.859738_dp
    CASE(45)
       ! Rhodium
       sep%defined = .TRUE.
       sep%uss         =  -20.513756_dp
       sep%upp         =  -40.045431_dp
       sep%udd         =  -35.818492_dp
       sep%beta(0)     =   -8.222141_dp
       sep%beta(1)     =  -15.556691_dp
       sep%beta(2)     =  -13.396182_dp
       sep%sto_exponents(0) =    1.324919_dp
       sep%sto_exponents(1) =    4.306111_dp
       sep%sto_exponents(2) =    2.901406_dp
       sep%zn(0)       =    0.809923_dp
       sep%zn(1)       =    6.898259_dp
       sep%zn(2)       =    0.643134_dp
       sep%gss         =    3.631179_dp
       sep%gsp         =    4.407820_dp
       sep%gpp         =   33.825599_dp
       sep%gp2         =   29.478305_dp
       sep%hsp         =    0.000092_dp
       sep%f0sd        =    1.775497_dp
       sep%g2sd        =    1.851571_dp
    CASE(46)
       ! Palladium
       sep%defined = .TRUE.
       sep%uss         =  -76.140196_dp
       sep%upp         =  -21.073362_dp
       sep%udd         =  -85.325301_dp
       sep%beta(0)     =   -8.038245_dp
       sep%beta(1)     =    0.740037_dp
       sep%beta(2)     =   -2.394498_dp
       sep%sto_exponents(0) =    1.658503_dp
       sep%sto_exponents(1) =    1.156718_dp
       sep%sto_exponents(2) =    2.219861_dp
       sep%zn(0)       =    1.794085_dp
       sep%zn(1)       =    6.158778_dp
       sep%zn(2)       =    1.630913_dp
       sep%gss         =    8.043535_dp
       sep%gsp         =    9.755042_dp
       sep%gpp         =   30.199556_dp
       sep%gp2         =   26.318284_dp
       sep%hsp         =    0.086121_dp
       sep%f0sd        =    8.004447_dp
       sep%g2sd        =    2.613148_dp
    CASE(47)
       ! Silver
       sep%defined = .TRUE.
       sep%uss         =  -25.484137_dp
       sep%upp         =  -36.116023_dp
       sep%udd         =  -35.668272_dp
       sep%beta(0)     =   -6.129623_dp
       sep%beta(1)     =    1.004115_dp
       sep%beta(2)     =  -69.238347_dp
       sep%sto_exponents(0) =    1.994004_dp
       sep%sto_exponents(1) =    0.681817_dp
       sep%sto_exponents(2) =    6.007328_dp
       sep%zn(0)       =    0.695514_dp
       sep%zn(1)       =    4.729949_dp
       sep%zn(2)       =    0.506522_dp
       sep%gss         =    3.118242_dp
       sep%gsp         =    3.785152_dp
       sep%gpp         =   23.193295_dp
       sep%gp2         =   20.212474_dp
       sep%hsp         =    0.000432_dp
       sep%f0sd        =    1.938327_dp
       sep%g2sd        =    1.071901_dp
    CASE(48)
       ! Cadmium
       sep%defined = .TRUE.
       sep%uss         =  -14.645792_dp
       sep%upp         =   -9.318664_dp
       sep%beta(0)     =  -11.613183_dp
       sep%beta(1)     =    1.663168_dp
       sep%sto_exponents(0) =    1.384108_dp
       sep%sto_exponents(1) =    1.957413_dp
       sep%gss         =    6.677284_dp
       sep%gsp         =    5.953373_dp
       sep%gpp         =   18.729843_dp
       sep%gp2         =    9.917452_dp
       sep%hsp         =    0.825192_dp
    CASE(49)
       ! Indium
       sep%defined = .TRUE.
       sep%uss         =  -28.339246_dp
       sep%upp         =  -23.373875_dp
       sep%beta(0)     =   -1.982376_dp
       sep%beta(1)     =   -3.330294_dp
       sep%sto_exponents(0) =    2.023087_dp
       sep%sto_exponents(1) =    2.106618_dp
       sep%gss         =    9.906091_dp
       sep%gsp         =   10.520060_dp
       sep%gpp         =    4.826006_dp
       sep%gp2         =    7.906563_dp
       sep%hsp         =    3.500299_dp
    CASE(50)
       ! Tin
       sep%defined = .TRUE.
       sep%uss         =  -29.888217_dp
       sep%upp         =  -22.156954_dp
       sep%beta(0)     =   -8.621087_dp
       sep%beta(1)     =   -4.989752_dp
       sep%sto_exponents(0) =    2.383941_dp
       sep%sto_exponents(1) =    2.057908_dp
       sep%gss         =    8.269655_dp
       sep%gsp         =    5.013349_dp
       sep%gpp         =    6.584874_dp
       sep%gp2         =    5.855159_dp
       sep%hsp         =    0.531212_dp
       sep%a           =   -1.004587_dp
       sep%b           =    4.706252_dp
       sep%c           =    1.180218_dp
    CASE(51)
       ! Antimony
       sep%defined = .TRUE.
       sep%uss         =  -41.688879_dp
       sep%upp         =  -39.541180_dp
       sep%udd         =   -6.581663_dp
       sep%beta(0)     =   -7.472322_dp
       sep%beta(1)     =   -5.940750_dp
       sep%beta(2)     =   -3.979108_dp
       sep%sto_exponents(0) =    2.391178_dp
       sep%sto_exponents(1) =    1.773006_dp
       sep%sto_exponents(2) =    2.465590_dp
       sep%zn(0)       =    5.993591_dp
       sep%zn(1)       =    6.145086_dp
       sep%zn(2)       =    5.704031_dp
       sep%gss         =   10.588832_dp
       sep%gsp         =    7.310023_dp
       sep%gpp         =    9.281609_dp
       sep%gp2         =    8.954081_dp
       sep%hsp         =    0.779112_dp
    CASE(52)
       ! Tellurium
       sep%defined = .TRUE.
       sep%uss         = -114.733316_dp
       sep%upp         =  -50.096389_dp
       sep%beta(0)     =  -70.001062_dp
       sep%beta(1)     =   -6.151642_dp
       sep%sto_exponents(0) =    2.769862_dp
       sep%sto_exponents(1) =    1.731319_dp
       sep%gss         =    7.030626_dp
       sep%gsp         =   12.601389_dp
       sep%gpp         =    7.883479_dp
       sep%gp2         =    6.973163_dp
       sep%hsp         =    5.000826_dp
    CASE(53)
       ! Iodine
       sep%defined = .TRUE.
       sep%uss         =  -59.973232_dp
       sep%upp         =  -56.459835_dp
       sep%udd         =  -28.822603_dp
       sep%beta(0)     =  -30.522481_dp
       sep%beta(1)     =   -5.942120_dp
       sep%beta(2)     =   -7.676107_dp
       sep%sto_exponents(0) =    4.498653_dp
       sep%sto_exponents(1) =    1.917072_dp
       sep%sto_exponents(2) =    1.875175_dp
       sep%zn(0)       =    9.135244_dp
       sep%zn(1)       =    6.888191_dp
       sep%zn(2)       =    3.791523_dp
       sep%gss         =    7.234759_dp
       sep%gsp         =    9.154406_dp
       sep%gpp         =    9.877466_dp
       sep%gp2         =    8.035916_dp
       sep%hsp         =    5.004215_dp
       sep%a           =   -0.035519_dp
       sep%b           =    1.744389_dp
       sep%c           =    1.223844_dp
    CASE(54)
       ! Xenon
       sep%defined = .TRUE.
       sep%uss         =  -18.270227_dp
       sep%upp         = -167.163063_dp
       sep%beta(0)     =   -3.980622_dp
       sep%beta(1)     =  -38.822792_dp
       sep%sto_exponents(0) =    2.759787_dp
       sep%sto_exponents(1) =    1.977446_dp
       sep%gss         =   20.000252_dp
       sep%gsp         =    4.175902_dp
       sep%gpp         =    2.305787_dp
       sep%gp2         =    4.063220_dp
       sep%hsp         =    4.418843_dp
    CASE(55)
       ! Cesium
       sep%defined = .TRUE.
       sep%uss         =   -3.748609_dp
       sep%upp         =   -2.348109_dp
       sep%beta(0)     =    2.287838_dp
       sep%beta(1)     =   -5.908071_dp
       sep%sto_exponents(0) =    5.956008_dp
       sep%sto_exponents(1) =    1.619485_dp
       sep%gss         =    6.464751_dp
       sep%gsp         =    4.004501_dp
       sep%gpp         =   13.775390_dp
       sep%gp2         =   12.912537_dp
       sep%hsp         =    1.026928_dp
    CASE(56)
       ! Barium
       sep%defined = .TRUE.
       sep%uss         =   -9.306985_dp
       sep%upp         =   -8.826713_dp
       sep%beta(0)     =   10.003125_dp
       sep%beta(1)     =   -6.335160_dp
       sep%sto_exponents(0) =    1.395379_dp
       sep%sto_exponents(1) =    1.430139_dp
       sep%gss         =    3.600823_dp
       sep%gsp         =    4.740579_dp
       sep%gpp         =    3.345166_dp
       sep%gp2         =    3.142783_dp
       sep%hsp         =    0.929429_dp
    CASE(57)
       ! Lanthanum
       sep%defined = .TRUE.
       sep%uss         =  -19.641953_dp
       sep%upp         =  -22.059431_dp
       sep%udd         =  -22.638986_dp
       sep%beta(0)     =    0.796727_dp
       sep%beta(1)     =  -10.856056_dp
       sep%beta(2)     =   -0.484922_dp
       sep%sto_exponents(0) =    2.673780_dp
       sep%sto_exponents(1) =    1.248192_dp
       sep%sto_exponents(2) =    1.688562_dp
       sep%zn(0)       =    1.617784_dp
       sep%zn(1)       =    4.331620_dp
       sep%zn(2)       =    2.285738_dp
       sep%gss         =    6.154440_dp
       sep%gsp         =    7.322704_dp
       sep%gpp         =   18.077465_dp
       sep%gp2         =   15.679057_dp
       sep%hsp         =    0.138601_dp
       sep%rho         =    2.511701_dp
       sep%f0sd        =    8.856858_dp
       sep%g2sd        =    7.925585_dp
    CASE(58:70)
    CASE(71)
       ! Lutetium
       sep%defined = .TRUE.
       sep%uss         =  -15.954994_dp
       sep%upp         =  -11.606213_dp
       sep%udd         =  -13.050056_dp
       sep%beta(0)     =   -5.590778_dp
       sep%beta(1)     =   -0.937679_dp
       sep%beta(2)     =   -7.737752_dp
       sep%sto_exponents(0) =    5.471741_dp
       sep%sto_exponents(1) =    1.712296_dp
       sep%sto_exponents(2) =    2.225892_dp
       sep%zn(0)       =    1.632335_dp
       sep%zn(1)       =    4.033128_dp
       sep%zn(2)       =    0.921999_dp
       sep%gss         =    6.209796_dp
       sep%gsp         =    7.379102_dp
       sep%gpp         =   16.831746_dp
       sep%gp2         =   14.598613_dp
       sep%hsp         =    0.209008_dp
       sep%rho         =    2.743262_dp
       sep%f0sd        =    3.924927_dp
       sep%g2sd        =    1.000946_dp
    CASE(72)
       ! Hafnium
       sep%defined = .TRUE.
       sep%uss         =  -22.375140_dp
       sep%upp         =  -13.081670_dp
       sep%udd         =  -20.637741_dp
       sep%beta(0)     =   -5.366351_dp
       sep%beta(1)     =  -21.550119_dp
       sep%beta(2)     =   -3.884443_dp
       sep%sto_exponents(0) =    3.085344_dp
       sep%sto_exponents(1) =    1.575819_dp
       sep%sto_exponents(2) =    1.840840_dp
       sep%zn(0)       =    0.946927_dp
       sep%zn(1)       =    3.538911_dp
       sep%zn(2)       =    0.940283_dp
       sep%gss         =    3.602338_dp
       sep%gsp         =    4.293729_dp
       sep%gpp         =   14.769194_dp
       sep%gp2         =   12.809708_dp
       sep%hsp         =    0.011028_dp
       sep%f0sd        =    4.842900_dp
       sep%g2sd        =    4.386101_dp
    CASE(73)
       ! Tantalum
       sep%defined = .TRUE.
       sep%uss         =  -39.009984_dp
       sep%upp         =    1.163975_dp
       sep%udd         =  -43.266315_dp
       sep%beta(0)     =  -17.199605_dp
       sep%beta(1)     =   -5.818839_dp
       sep%beta(2)     =   -9.816794_dp
       sep%sto_exponents(0) =    4.578087_dp
       sep%sto_exponents(1) =    4.841244_dp
       sep%sto_exponents(2) =    1.838249_dp
       sep%zn(0)       =    1.741367_dp
       sep%zn(1)       =    3.430157_dp
       sep%zn(2)       =    2.311198_dp
       sep%gss         =    6.624580_dp
       sep%gsp         =    7.805321_dp
       sep%gpp         =   14.315323_dp
       sep%gp2         =   12.416054_dp
       sep%hsp         =    0.577263_dp
       sep%f0sd        =    8.544427_dp
       sep%g2sd        =    2.074254_dp
    CASE(74)
       ! Tungsten
       sep%defined = .TRUE.
       sep%uss         =  -44.524950_dp
       sep%upp         =  -40.011500_dp
       sep%udd         =  -46.490410_dp
       sep%beta(0)     =  -16.946460_dp
       sep%beta(1)     =    5.623170_dp
       sep%beta(2)     =   -2.947340_dp
       sep%sto_exponents(0) =    2.664560_dp
       sep%sto_exponents(1) =    1.624010_dp
       sep%sto_exponents(2) =    1.794400_dp
       sep%zn(0)       =    1.498860_dp
       sep%zn(1)       =    1.965900_dp
       sep%zn(2)       =    1.876450_dp
       sep%gss         =    5.702025_dp
       sep%gsp         =    6.323145_dp
       sep%gpp         =    8.204433_dp
       sep%gp2         =    7.115919_dp
       sep%hsp         =    1.319912_dp
       sep%f0sd        =    7.788180_dp
       sep%g2sd        =    1.684940_dp
    CASE(75)
       ! Rhenium
       sep%defined = .TRUE.
       sep%uss         =  -41.291342_dp
       sep%upp         =  -35.089592_dp
       sep%udd         =  -44.178985_dp
       sep%beta(0)     =    3.830075_dp
       sep%beta(1)     =   -1.638530_dp
       sep%beta(2)     =   -1.414411_dp
       sep%sto_exponents(0) =    2.411839_dp
       sep%sto_exponents(1) =    1.815351_dp
       sep%sto_exponents(2) =    2.522766_dp
       sep%zn(0)       =    1.680823_dp
       sep%zn(1)       =    1.331218_dp
       sep%zn(2)       =    1.490623_dp
       sep%gss         =    6.394256_dp
       sep%gsp         =    5.555571_dp
       sep%gpp         =    5.555669_dp
       sep%gp2         =    4.818577_dp
       sep%hsp         =    1.220913_dp
       sep%f0sd        =    5.442812_dp
       sep%g2sd        =    2.376279_dp
    CASE(76)
       ! Osmium
       sep%defined = .TRUE.
       sep%uss         =  -26.434080_dp
       sep%upp         =  -48.739500_dp
       sep%udd         =  -55.837880_dp
       sep%beta(0)     =  -12.508730_dp
       sep%beta(1)     =    0.846880_dp
       sep%beta(2)     =    5.164360_dp
       sep%sto_exponents(0) =    3.031000_dp
       sep%sto_exponents(1) =    1.593960_dp
       sep%sto_exponents(2) =    1.775570_dp
       sep%zn(0)       =    1.844700_dp
       sep%zn(1)       =    1.564220_dp
       sep%zn(2)       =    1.770010_dp
       sep%gss         =    7.017683_dp
       sep%gsp         =    6.384200_dp
       sep%gpp         =    6.528073_dp
       sep%gp2         =    5.661968_dp
       sep%hsp         =    1.508926_dp
       sep%f0sd        =    2.021170_dp
       sep%g2sd        =    1.392130_dp
    CASE(77)
       ! Iridium
       sep%defined = .TRUE.
       sep%uss         =  -29.703974_dp
       sep%upp         =  -38.210924_dp
       sep%udd         =  -32.538202_dp
       sep%beta(0)     =  -10.943427_dp
       sep%beta(1)     =    2.908880_dp
       sep%beta(2)     =   -3.791731_dp
       sep%sto_exponents(0) =    1.500907_dp
       sep%sto_exponents(1) =    4.106373_dp
       sep%sto_exponents(2) =    2.676047_dp
       sep%zn(0)       =    0.927246_dp
       sep%zn(1)       =    3.191892_dp
       sep%zn(2)       =    0.662007_dp
       sep%gss         =    3.527467_dp
       sep%gsp         =    4.203820_dp
       sep%gpp         =   13.320955_dp
       sep%gp2         =   11.553612_dp
       sep%hsp         =    0.018501_dp
       sep%f0sd        =    2.627170_dp
       sep%g2sd        =    2.996029_dp
    CASE(78)
       ! Platinum
       sep%defined = .TRUE.
       sep%uss         =  -73.516173_dp
       sep%upp         =  -68.320056_dp
       sep%udd         =  -76.598873_dp
       sep%beta(0)     =    1.151418_dp
       sep%beta(1)     =    3.298694_dp
       sep%beta(2)     =  -18.044737_dp
       sep%sto_exponents(0) =    2.301264_dp
       sep%sto_exponents(1) =    1.662404_dp
       sep%sto_exponents(2) =    3.168852_dp
       sep%zn(0)       =    2.270699_dp
       sep%zn(1)       =    1.949896_dp
       sep%zn(2)       =    1.713856_dp
       sep%gss         =    8.638286_dp
       sep%gsp         =    7.922254_dp
       sep%gpp         =    8.137643_dp
       sep%gp2         =    7.057990_dp
       sep%hsp         =    1.892617_dp
       sep%f0sd        =    7.098591_dp
       sep%g2sd        =    4.484183_dp
    CASE(79)
       ! Gold
       sep%defined = .TRUE.
       sep%uss         =  -95.041846_dp
       sep%upp         =  -63.890158_dp
       sep%udd         =  -88.066087_dp
       sep%beta(0)     =   -7.479625_dp
       sep%beta(1)     =    3.664356_dp
       sep%beta(2)     =  -61.715468_dp
       sep%sto_exponents(0) =    1.814169_dp
       sep%sto_exponents(1) =    1.618657_dp
       sep%sto_exponents(2) =    5.053167_dp
       sep%zn(0)       =    2.444680_dp
       sep%zn(1)       =    7.014990_dp
       sep%zn(2)       =    1.777089_dp
       sep%gss         =    9.300152_dp
       sep%gsp         =   11.073443_dp
       sep%gpp         =   29.276168_dp
       sep%gp2         =   25.391984_dp
       sep%hsp         =    0.144384_dp
       sep%f0sd        =    8.827257_dp
       sep%g2sd        =    4.915625_dp
    CASE(80)
       ! Mercury
       sep%defined = .TRUE.
       sep%uss         =  -17.608732_dp
       sep%upp         =  -18.369417_dp
       sep%beta(0)     =   -3.045239_dp
       sep%beta(1)     =   -5.693556_dp
       sep%sto_exponents(0) =    2.104896_dp
       sep%sto_exponents(1) =    1.516293_dp
       sep%gss         =    6.372822_dp
       sep%gsp         =   10.143176_dp
       sep%gpp         =   10.397393_dp
       sep%gp2         =   14.794056_dp
       sep%hsp         =    0.926128_dp
    CASE(81)
       ! Thallium
       sep%defined = .TRUE.
       sep%uss         =  -29.518621_dp
       sep%upp         =  -29.826907_dp
       sep%beta(0)     =   -7.230170_dp
       sep%beta(1)     =   -7.575544_dp
       sep%sto_exponents(0) =    3.335883_dp
       sep%sto_exponents(1) =    1.766141_dp
       sep%gss         =    5.015118_dp
       sep%gsp         =   13.932049_dp
       sep%gpp         =   10.495551_dp
       sep%gp2         =   10.526198_dp
       sep%hsp         =    0.293760_dp
    CASE(82)
       ! Lead
       sep%defined = .TRUE.
       sep%uss         =  -35.038145_dp
       sep%upp         =  -25.413401_dp
       sep%beta(0)     =   -8.323792_dp
       sep%beta(1)     =   -2.237891_dp
       sep%sto_exponents(0) =    2.368901_dp
       sep%sto_exponents(1) =    1.685246_dp
       sep%gss         =    5.254128_dp
       sep%gsp         =    7.061016_dp
       sep%gpp         =    6.818551_dp
       sep%gp2         =    5.603019_dp
       sep%hsp         =    1.018819_dp
       sep%a           =   -0.239463_dp
       sep%b           =    5.444338_dp
       sep%c           =    1.613682_dp
    CASE(83)
       ! Bismuth
       sep%defined = .TRUE.
       sep%uss         =  -42.409177_dp
       sep%upp         =  -36.393746_dp
       sep%beta(0)     =  -34.951578_dp
       sep%beta(1)     =   -7.359060_dp
       sep%sto_exponents(0) =    3.702377_dp
       sep%sto_exponents(1) =    1.872327_dp
       sep%gss         =    5.851803_dp
       sep%gsp         =    6.790583_dp
       sep%gpp         =    8.389442_dp
       sep%gp2         =    7.724219_dp
       sep%hsp         =    0.295606_dp
    CASE(84:nelem)
    END SELECT
    ! Finalize parameters for the element
    CALL valence_electrons(sep, extended_basis_set=.TRUE., error=error)
    CALL calpar(z,sep,error)
    CALL convert_param_to_cp2k(sep, error)
    CALL init_pm6_pair_params(error)

    sep%aab(0:nelem) = aab_pm6(z,0:nelem)
    sep%xab(0:nelem) = xab_pm6(z,0:nelem)
    !    ! We seem to have also C-H as a special case
    !    IF (z == 1) sep%aab(6)= sep%aab(6)/bohr
    !    IF (z == 1) sep%aab(7)= sep%aab(7)/bohr
    !    IF (z == 1) sep%aab(8)= sep%aab(8)/bohr
    !    IF (z == 6) sep%aab(1)= sep%aab(1)/bohr
    !    IF (z == 7) sep%aab(1)= sep%aab(1)/bohr
    !    IF (z == 8) sep%aab(1)= sep%aab(1)/bohr

  END SUBROUTINE pm6_default_parameter

! *****************************************************************************
!> \brief Default parameter sets for semi empirical models: PDDG
!> \param sep ...
!> \param z ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pdg_default_parameter ( sep, z, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER                                  :: z
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pdg_default_parameter', &
      routineP = moduleN//':'//routineN

    sep%typ    = do_method_pdg
    sep%name   = ptable(z)%name
    sep%eheat  = ptable(z)%heat_of_formation / kcalmol
    sep%z      = z

    SELECT CASE (z)
    CASE DEFAULT
       CALL stop_program(routineN,moduleN,__LINE__,&
            "invalid atomic number")
    CASE (0)
       sep%defined = .TRUE.
    CASE (1)
       !Hydrogen
       sep%defined          = .TRUE.
       sep%uss              =  -12.893272_dp
       sep%beta(0)          =   -6.152654_dp
       sep%sto_exponents(0) =    0.972786_dp
       sep%alp              =    3.381686_dp
       sep%eisol            =  -13.120566_dp
       sep%gss              =   14.794208_dp
       sep%am               =    0.5_dp/0.919616_dp
       sep%ad               =    0.5_dp/0.919616_dp
       sep%aq               =    0.5_dp/0.919616_dp
       sep%fn1(1)           =    1.122244_dp
       sep%fn2(1)           =    4.707790_dp
       sep%fn3(1)           =    1.547099_dp
       sep%fn1(2)           =   -1.069737_dp
       sep%fn2(2)           =    5.857995_dp
       sep%fn3(2)           =    1.567893_dp
       sep%pre(1)           =    0.057193_dp
       sep%pre(2)           =   -0.034823_dp
       sep%d(1)             =    0.663395_dp
       sep%d(2)             =    1.081901_dp
    CASE (2:5)
    CASE (6)
       !Carbon
       sep%defined          = .TRUE.
       sep%uss              =  -48.241241_dp
       sep%upp              =  -36.461256_dp
       sep%beta(0)          =  -11.952818_dp
       sep%beta(1)          =   -9.922411_dp
       sep%sto_exponents(0) =    1.567864_dp
       sep%sto_exponents(1) =    1.846659_dp
       sep%alp              =    2.725772_dp
       sep%eisol            = -113.428242_dp
       sep%gss              =   11.200708_dp
       sep%gsp              =   10.265027_dp
       sep%gpp              =   10.796292_dp
       sep%gp2              =    9.042566_dp
       sep%hsp              =    2.290980_dp
       sep%dd               =    0.831413_dp
       sep%qq               =    0.663222_dp
       sep%am               =    0.5_dp/1.214657_dp
       sep%ad               =    0.5_dp/0.848467_dp
       sep%aq               =    0.5_dp/0.652785_dp
       sep%fn1(1)           =    0.048906_dp
       sep%fn2(1)           =    5.765340_dp
       sep%fn3(1)           =    1.682232_dp
       sep%fn1(2)           =    0.047697_dp
       sep%fn2(2)           =    5.973721_dp
       sep%fn3(2)           =    0.894406_dp
       sep%pre(1)           =   -0.000743_dp
       sep%pre(2)           =    0.000985_dp
       sep%d(1)             =    0.836915_dp
       sep%d(2)             =    1.585236_dp
    CASE (7)
       !Nitrogen
       sep%defined          = .TRUE.
       sep%uss              =  -49.454546_dp
       sep%upp              =  -47.757406_dp
       sep%beta(0)          =  -14.117230_dp
       sep%beta(1)          =  -19.938509_dp
       sep%sto_exponents(0) =    2.035807_dp
       sep%sto_exponents(1) =    2.324327_dp
       sep%alp              =    2.849124_dp
       sep%eisol            = -158.416205_dp
       sep%gss              =   11.904787_dp
       sep%gsp              =    7.348565_dp
       sep%gpp              =   11.754672_dp
       sep%gp2              =   10.807277_dp
       sep%hsp              =    1.136713_dp
       sep%dd               =    0.654855_dp
       sep%qq               =    0.526924_dp
       sep%am               =    0.5_dp/1.142818_dp
       sep%ad               =    0.5_dp/0.991235_dp
       sep%aq               =    0.5_dp/0.676704_dp
       sep%fn1(1)           =    1.513320_dp
       sep%fn2(1)           =    5.904394_dp
       sep%fn3(1)           =    1.728376_dp
       sep%fn1(2)           =   -1.511892_dp
       sep%fn2(2)           =    6.030014_dp
       sep%fn3(2)           =    1.734108_dp
       sep%pre(1)           =   -0.003160_dp
       sep%pre(2)           =    0.012501_dp
       sep%d(1)             =    1.004172_dp
       sep%d(2)             =    1.516336_dp
    CASE (8)
       !Oxygen
       sep%defined          = .TRUE.
       sep%uss              =  -87.412505_dp
       sep%upp              =  -72.183070_dp
       sep%beta(0)          =  -44.874553_dp
       sep%beta(1)          =  -24.601939_dp
       sep%sto_exponents(0) =    3.814565_dp
       sep%sto_exponents(1) =    2.318011_dp
       sep%alp              =    3.225309_dp
       sep%eisol            = -292.188766_dp
       sep%gss              =   15.755760_dp
       sep%gsp              =   10.621160_dp
       sep%gpp              =   13.654016_dp
       sep%gp2              =   12.406095_dp
       sep%hsp              =    0.593883_dp
       sep%dd               =    0.403741_dp
       sep%qq               =    0.528360_dp
       sep%am               =    0.5_dp/0.863494_dp
       sep%ad               =    0.5_dp/0.936266_dp
       sep%aq               =    0.5_dp/0.624291_dp
       sep%fn1(1)           =   -1.138455_dp
       sep%fn2(1)           =    6.000043_dp
       sep%fn3(1)           =    1.622362_dp
       sep%fn1(2)           =    1.146007_dp
       sep%fn2(2)           =    5.963494_dp
       sep%fn3(2)           =    1.614788_dp
       sep%pre(1)           =   -0.001000_dp
       sep%pre(2)           =   -0.001522_dp
       sep%d(1)             =    1.360685_dp
       sep%d(2)             =    1.366407_dp
    CASE (9)
       !Fluorine
       sep%defined          = .TRUE.
       sep%uss              = -110.400432_dp
       sep%upp              = -106.395264_dp
       sep%beta(0)          =  -50.937301_dp
       sep%beta(1)          =  -31.636976_dp
       sep%sto_exponents(0) =    5.538033_dp
       sep%sto_exponents(1) =    2.538066_dp
       sep%alp              =    3.200571_dp
       sep%eisol            = -442.457133_dp
       sep%gss              =   10.496667_dp
       sep%gsp              =   16.073689_dp
       sep%gpp              =   14.817256_dp
       sep%gp2              =   14.418393_dp
       sep%hsp              =    0.727763_dp
       sep%dd               =    0.246601_dp
       sep%qq               =    0.482551_dp
       sep%am               =    0.5_dp/1.296126_dp
       sep%ad               =    0.5_dp/0.6346333_dp
       sep%aq               =    0.5_dp/0.805802_dp
       sep%fn1(1)           =   -0.008079_dp
       sep%fn2(1)           =    5.938969_dp
       sep%fn3(1)           =    1.863949_dp
       sep%fn1(2)           =   -0.002659_dp
       sep%fn2(2)           =    5.925105_dp
       sep%fn3(2)           =    2.388864_dp
       sep%pre(1)           =   -0.012866_dp
       sep%pre(2)           =   -0.007315_dp
       sep%d(1)             =    1.305681_dp
       sep%d(2)             =    1.842572_dp
    CASE (10:13)
    CASE (14)
       !Silicon
       sep%defined          = .TRUE.
       sep%uss              =  -26.332522_dp
       sep%upp              =  -22.602540_dp
       sep%beta(0)          =   -3.376445_dp
       sep%beta(1)          =   -3.620969_dp
       sep%sto_exponents(0) =    1.586389_dp
       sep%sto_exponents(1) =    1.485958_dp
       sep%alp              =    2.215157_dp
       sep%eisol            =  -66.839000_dp
       sep%gss              =    5.047196_dp
       sep%gsp              =    5.949057_dp
       sep%gpp              =    6.759367_dp
       sep%gp2              =    5.161297_dp
       sep%hsp              =    0.919832_dp
       sep%dd               =    1.310515_dp
       sep%qq               =    1.126089_dp
       sep%am               =    0.5_dp/2.695556_dp
       sep%ad               =    0.5_dp/1.630757_dp
       sep%aq               =    0.5_dp/0.949200_dp
       sep%fn1(1)           =   -0.071314_dp
       sep%fn2(1)           =    6.000000_dp
       sep%fn3(1)           =    0.237995_dp
       sep%fn1(2)           =    0.089451_dp
       sep%fn2(2)           =    6.000000_dp
       sep%fn3(2)           =    1.897728_dp
       sep%pre(1)           =   -0.091928_dp
       sep%pre(2)           =   -0.040753_dp
       sep%d(1)             =    1.163190_dp
       sep%d(2)             =    2.190526_dp
    CASE (15)
       !Phospherus
       sep%defined          = .TRUE.
       sep%uss              =  -37.882113_dp
       sep%upp              =  -30.312979_dp
       sep%beta(0)          =  -12.676297_dp
       sep%beta(1)          =   -7.093318_dp
       sep%sto_exponents(0) =    2.395882_dp
       sep%sto_exponents(1) =    1.742213_dp
       sep%alp              =    2.005294_dp
       sep%eisol            = -117.212854_dp
       sep%gss              =    7.801615_dp
       sep%gsp              =    5.186949_dp
       sep%gpp              =    6.618478_dp
       sep%gp2              =    6.062002_dp
       sep%hsp              =    1.542809_dp
       sep%dd               =    0.893978_dp
       sep%qq               =    0.960457_dp
       sep%am               =    0.5_dp/1.743870_dp
       sep%ad               =    0.5_dp/1.050851_dp
       sep%aq               =    0.5_dp/1.208907_dp
       sep%fn1(1)           =   -0.398055_dp
       sep%fn2(1)           =    1.997272_dp
       sep%fn3(1)           =    0.950073_dp
       sep%fn1(2)           =   -0.079653_dp
       sep%fn2(2)           =    1.998360_dp
       sep%fn3(2)           =    2.336959_dp
       sep%pre(1)           =    0.462741_dp
       sep%pre(2)           =   -0.020444_dp
       sep%d(1)             =    0.714296_dp
       sep%d(2)             =    2.041209_dp
    CASE (16)
       !Sulfur
       sep%defined          = .TRUE.
       sep%uss              =  -43.906366_dp
       sep%upp              =  -43.461348_dp
       sep%beta(0)          =   -2.953912_dp
       sep%beta(1)          =   -8.507779_dp
       sep%sto_exponents(0) =    1.012002_dp
       sep%sto_exponents(1) =    1.876999_dp
       sep%alp              =    2.539751_dp
       sep%eisol            = -166.336554_dp
       sep%gss              =    8.964667_dp
       sep%gsp              =    6.785936_dp
       sep%gpp              =    9.968164_dp
       sep%gp2              =    7.970247_dp
       sep%hsp              =    4.041836_dp
       sep%dd               =    1.006989_dp
       sep%qq               =    0.891487_dp
       sep%am               =    0.5_dp/1.517625_dp
       sep%ad               =    0.5_dp/0.711672_dp
       sep%aq               =    0.5_dp/0.754336_dp
       sep%fn1(1)           =   -0.330692_dp
       sep%fn2(1)           =    6.000000_dp
       sep%fn3(1)           =    0.823837_dp
       sep%fn1(2)           =    0.024171_dp
       sep%fn2(2)           =    6.000000_dp
       sep%fn3(2)           =    2.017756_dp
       sep%pre(1)           =    0.120434_dp
       sep%pre(2)           =   -0.002663_dp
       sep%d(1)             =    0.672870_dp
       sep%d(2)             =    2.032340_dp
    CASE (17)
       !Chlorine
       sep%defined          = .TRUE.
       sep%uss              =  -95.094434_dp
       sep%upp              =  -53.921651_dp
       sep%beta(0)          =  -26.913129_dp
       sep%beta(1)          =  -14.991178_dp
       sep%sto_exponents(0) =    2.548268_dp
       sep%sto_exponents(1) =    2.284624_dp
       sep%sto_exponents(2) =    1.000000_dp
       sep%alp              =    2.497617_dp
       sep%eisol            = -305.715201_dp
       sep%gss              =   16.013601_dp
       sep%gsp              =    8.048115_dp
       sep%gpp              =    7.522215_dp
       sep%gp2              =    7.504154_dp
       sep%hsp              =    3.481153_dp
       sep%dd               =    0.827561_dp
       sep%qq               =    0.732427_dp
       sep%am               =    0.5_dp/0.8495900_dp
       sep%ad               =    0.5_dp/0.696164_dp
       sep%aq               =    0.5_dp/2.299104_dp
       sep%fn1(1)           =   -0.112222_dp
       sep%fn2(1)           =    5.963719_dp
       sep%fn3(1)           =    1.027719_dp
       sep%fn1(2)           =   -0.013061_dp
       sep%fn2(2)           =    1.999556_dp
       sep%fn3(2)           =    2.286377_dp
       sep%pre(1)           =   -0.016552_dp
       sep%pre(2)           =   -0.016646_dp
       sep%d(1)             =    1.727690_dp
       sep%d(2)             =    1.784655_dp
    CASE (18:34)
    CASE (35)
       !Bromine
       sep%defined          = .TRUE.
       sep%uss              = -115.841963_dp
       sep%upp              =  -74.205146_dp
       sep%beta(0)          =  -21.538044_dp
       sep%beta(1)          =   -8.524764_dp
       sep%sto_exponents(0) =    4.345079_dp
       sep%sto_exponents(1) =    2.190961_dp
       sep%sto_exponents(2) =    1.000000_dp
       sep%alp              =    2.424673_dp
       sep%eisol            = -351.013887_dp
       sep%gss              =   15.943425_dp
       sep%gsp              =   16.061680_dp
       sep%gpp              =    8.282763_dp
       sep%gp2              =    7.816849_dp
       sep%hsp              =    0.578869_dp
       sep%dd               =    0.473860_dp
       sep%qq               =    0.968214_dp
       sep%am               =    0.5_dp/0.853330_dp
       sep%ad               =    0.5_dp/1.046430_dp
       sep%aq               =    0.5_dp/1.280643_dp
       sep%fn1(1)           =    0.9613620_dp
       sep%fn2(1)           =    6.0136000_dp
       sep%fn3(1)           =    2.3404450_dp
       sep%fn1(2)           =   -0.9548834_dp
       sep%fn2(2)           =    5.9763290_dp
       sep%fn3(2)           =    2.3487450_dp
       sep%pre(1)           =   -0.013772_dp
       sep%pre(2)           =    0.008849_dp
       sep%d(1)             =    1.852030_dp
       sep%d(2)             =    2.338958_dp
    CASE (36:52)
    CASE (53)
       !Iodine
       sep%defined          = .TRUE.
       sep%uss              =  -96.664174_dp
       sep%upp              =  -61.167137_dp
       sep%beta(0)          =  -16.592621_dp
       sep%beta(1)          =   -6.599816_dp
       sep%sto_exponents(0) =    5.062801_dp
       sep%sto_exponents(1) =    2.417757_dp
       sep%sto_exponents(2) =    1.000000_dp
       sep%alp              =    1.978170_dp
       sep%eisol            = -291.537869_dp
       sep%gss              =   13.631943_dp
       sep%gsp              =   14.990406_dp
       sep%gpp              =    7.288330_dp
       sep%gp2              =    5.966407_dp
       sep%hsp              =    2.630035_dp
       sep%dd               =    0.407261_dp
       sep%qq               =    1.062574_dp
       sep%am               =    0.5_dp/0.998024_dp
       sep%ad               =    0.5_dp/0.532290_dp
       sep%aq               =    0.5_dp/0.979783_dp
       sep%fn1(1)           =   -0.136003_dp
       sep%fn2(1)           =    3.852912_dp
       sep%fn3(1)           =    1.697455_dp
       sep%fn1(2)           =   -0.037287_dp
       sep%fn2(2)           =    5.229264_dp
       sep%fn3(2)           =    2.768669_dp
       sep%pre(1)           =    0.012901_dp
       sep%pre(2)           =   -0.012825_dp
       sep%d(1)             =    1.994299_dp
       sep%d(2)             =    2.263417_dp
    CASE (54:nelem)
    END SELECT
    ! Finalize parameters for the element
    CALL valence_electrons(sep, extended_basis_set=.FALSE., error=error)
    CALL calpar(z,sep, error)
    CALL convert_param_to_cp2k(sep, error)

  END SUBROUTINE pdg_default_parameter

! *****************************************************************************
!> \brief Default pair parameter for: PM6
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_pm6_pair_params(error)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_pm6_pair_params', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    INTEGER, SAVE                            :: init = 0
    REAL(dp)                                 :: aab, xab

    IF (init /= 0) RETURN

    xab_pm6 = 0._dp
    aab_pm6 = 0._dp

    xab_pm6(  1,  1)=        2.2435870171_dp
    xab_pm6(  2,  1)=        2.3711988926_dp
    xab_pm6(  2,  2)=        3.4509000778_dp
    xab_pm6(  3,  1)=        2.1919848919_dp
    xab_pm6(  3,  2)=        9.2736759186_dp
    xab_pm6(  3,  3)=       16.1163845062_dp
    xab_pm6(  4,  1)=        2.5628309250_dp
    xab_pm6(  4,  2)=       12.5448780060_dp
    xab_pm6(  4,  3)=        3.2871649265_dp
    xab_pm6(  4,  4)=        0.2386330068_dp
    xab_pm6(  5,  1)=        1.3213939667_dp
    xab_pm6(  5,  2)=        1.9741699696_dp
    xab_pm6(  5,  3)=        7.8860177994_dp
    xab_pm6(  5,  4)=        1.1517920494_dp
    xab_pm6(  5,  5)=        3.5936191082_dp
    xab_pm6(  6,  1)=        0.2165060000_dp
    xab_pm6(  6,  2)=        3.2139708996_dp
    xab_pm6(  6,  3)=       16.1800022125_dp
    xab_pm6(  6,  4)=       25.0358791351_dp
    xab_pm6(  6,  5)=        1.8748589754_dp
    xab_pm6(  6,  6)=        0.8135100007_dp
    xab_pm6(  7,  1)=        0.1755059958_dp
    xab_pm6(  7,  2)=        1.0778609514_dp
    xab_pm6(  7,  3)=        2.8234028816_dp
    xab_pm6(  7,  4)=        1.7406049967_dp
    xab_pm6(  7,  5)=        0.9528819919_dp
    xab_pm6(  7,  6)=        0.8599489927_dp
    xab_pm6(  7,  7)=        0.6753129959_dp
    xab_pm6(  8,  1)=        0.1922950000_dp
    xab_pm6(  8,  2)=        6.6845250130_dp
    xab_pm6(  8,  3)=        1.9685980082_dp
    xab_pm6(  8,  4)=        3.2181549072_dp
    xab_pm6(  8,  5)=        1.2698010206_dp
    xab_pm6(  8,  6)=        0.9902110100_dp
    xab_pm6(  8,  7)=        0.7647560239_dp
    xab_pm6(  8,  8)=        0.5351120234_dp
    xab_pm6(  9,  1)=        0.8158019781_dp
    xab_pm6(  9,  2)=        0.7451069951_dp
    xab_pm6(  9,  3)=        1.9759850502_dp
    xab_pm6(  9,  4)=        3.8829929829_dp
    xab_pm6(  9,  5)=        0.8627610207_dp
    xab_pm6(  9,  6)=        0.7329679728_dp
    xab_pm6(  9,  7)=        0.6358540058_dp
    xab_pm6(  9,  8)=        0.6742510200_dp
    xab_pm6(  9,  9)=        0.6813430190_dp
    xab_pm6( 10,  1)=        5.5350208282_dp
    xab_pm6( 10,  2)=        1.9609240294_dp
    xab_pm6( 10,  3)=        0.7049580216_dp
    xab_pm6( 10,  4)=        0.3926280141_dp
    xab_pm6( 10,  5)=        2.7641398907_dp
    xab_pm6( 10,  6)=        5.4687800407_dp
    xab_pm6( 10,  7)=       29.9996089935_dp
    xab_pm6( 10,  8)=        0.7638990283_dp
    xab_pm6( 10,  9)=        2.7067539692_dp
    xab_pm6( 10, 10)=        2.7948300838_dp
    xab_pm6( 11,  1)=        0.2078309953_dp
    xab_pm6( 11,  2)=        4.2825169563_dp
    xab_pm6( 11,  3)=        0.8814820051_dp
    xab_pm6( 11,  4)=        3.1216199398_dp
    xab_pm6( 11,  5)=        3.1886079311_dp
    xab_pm6( 11,  6)=        4.5204291344_dp
    xab_pm6( 11,  7)=        8.5863866806_dp
    xab_pm6( 11,  8)=        3.2700788975_dp
    xab_pm6( 11,  9)=        7.0473508835_dp
    xab_pm6( 11, 10)=        1.3430370092_dp
    xab_pm6( 11, 11)=        0.2871370018_dp
    xab_pm6( 12,  1)=        7.7582368851_dp
    xab_pm6( 12,  2)=        3.7258501053_dp
    xab_pm6( 12,  3)=        2.4902501106_dp
    xab_pm6( 12,  4)=        2.0663919449_dp
    xab_pm6( 12,  5)=        6.1467008591_dp
    xab_pm6( 12,  6)=       10.5176897049_dp
    xab_pm6( 12,  7)=        1.2080750465_dp
    xab_pm6( 12,  8)=        1.5357340574_dp
    xab_pm6( 12,  9)=        5.8590230942_dp
    xab_pm6( 12, 10)=        1.2148590088_dp
    xab_pm6( 12, 11)=        8.6756191254_dp
    xab_pm6( 12, 12)=        0.4656449854_dp
    xab_pm6( 13,  1)=        2.9583790302_dp
    xab_pm6( 13,  2)=        2.7014000416_dp
    xab_pm6( 13,  3)=        1.1068190336_dp
    xab_pm6( 13,  4)=        5.0372138023_dp
    xab_pm6( 13,  5)=        2.7414789200_dp
    xab_pm6( 13,  6)=        2.9280560017_dp
    xab_pm6( 13,  7)=        1.3452019691_dp
    xab_pm6( 13,  8)=        2.1313960552_dp
    xab_pm6( 13,  9)=        1.9756350517_dp
    xab_pm6( 13, 10)=        1.7092000246_dp
    xab_pm6( 13, 11)=        2.0718469620_dp
    xab_pm6( 13, 12)=       13.4724426270_dp
    xab_pm6( 13, 13)=        2.1391999722_dp
    xab_pm6( 14,  1)=        0.9241960049_dp
    xab_pm6( 14,  2)=        1.8535829782_dp
    xab_pm6( 14,  3)=        3.0907909870_dp
    xab_pm6( 14,  4)=        0.6234329939_dp
    xab_pm6( 14,  5)=        1.0282870531_dp
    xab_pm6( 14,  6)=        0.7857450247_dp
    xab_pm6( 14,  7)=        0.5929719806_dp
    xab_pm6( 14,  8)=        0.7510949969_dp
    xab_pm6( 14,  9)=        0.5435159802_dp
    xab_pm6( 14, 10)=       14.3786764145_dp
    xab_pm6( 14, 11)=        9.2376441956_dp
    xab_pm6( 14, 12)=       29.9945201874_dp
    xab_pm6( 14, 13)=        2.0000000000_dp
    xab_pm6( 14, 14)=        0.2734769881_dp
    xab_pm6( 15,  1)=        1.2349859476_dp
    xab_pm6( 15,  2)=        1.4902180433_dp
    xab_pm6( 15,  3)=        1.1229499578_dp
    xab_pm6( 15,  4)=        1.6848310232_dp
    xab_pm6( 15,  5)=        1.4508860111_dp
    xab_pm6( 15,  6)=        0.9795119762_dp
    xab_pm6( 15,  7)=        0.9721540213_dp
    xab_pm6( 15,  8)=        0.8787050247_dp
    xab_pm6( 15,  9)=        0.5145750046_dp
    xab_pm6( 15, 10)=        0.7749540210_dp
    xab_pm6( 15, 11)=        2.8370950222_dp
    xab_pm6( 15, 12)=        1.1778810024_dp
    xab_pm6( 15, 13)=        5.0508160591_dp
    xab_pm6( 15, 14)=       13.2391214371_dp
    xab_pm6( 15, 15)=        0.9025009871_dp
    xab_pm6( 16,  1)=        0.8497120142_dp
    xab_pm6( 16,  2)=        0.4376179874_dp
    xab_pm6( 16,  3)=        2.6425020695_dp
    xab_pm6( 16,  4)=        3.7915649414_dp
    xab_pm6( 16,  5)=        1.1253939867_dp
    xab_pm6( 16,  6)=        0.6668490171_dp
    xab_pm6( 16,  7)=        0.7387099862_dp
    xab_pm6( 16,  8)=        0.7472149730_dp
    xab_pm6( 16,  9)=        0.3752509952_dp
    xab_pm6( 16, 10)=        3.2961599827_dp
    xab_pm6( 16, 11)=        0.8524339795_dp
    xab_pm6( 16, 12)=        0.5007479787_dp
    xab_pm6( 16, 13)=        2.3473839760_dp
    xab_pm6( 16, 14)=        0.8766580224_dp
    xab_pm6( 16, 15)=        0.5622659922_dp
    xab_pm6( 16, 16)=        0.4738560021_dp
    xab_pm6( 17,  1)=        0.7548310161_dp
    xab_pm6( 17,  2)=        0.2729640007_dp
    xab_pm6( 17,  3)=        4.2277941704_dp
    xab_pm6( 17,  4)=        2.5072751045_dp
    xab_pm6( 17,  5)=        0.8221290112_dp
    xab_pm6( 17,  6)=        0.5157870054_dp
    xab_pm6( 17,  7)=        0.5207449794_dp
    xab_pm6( 17,  8)=        0.5855100155_dp
    xab_pm6( 17,  9)=        0.4111239910_dp
    xab_pm6( 17, 10)=        0.1251329929_dp
    xab_pm6( 17, 11)=        1.3578939438_dp
    xab_pm6( 17, 12)=        2.4308559895_dp
    xab_pm6( 17, 13)=        2.1534509659_dp
    xab_pm6( 17, 14)=        0.5130000114_dp
    xab_pm6( 17, 15)=        0.3523609936_dp
    xab_pm6( 17, 16)=        0.3569709957_dp
    xab_pm6( 17, 17)=        0.3329190016_dp
    xab_pm6( 18,  1)=        3.9334449768_dp
    xab_pm6( 18,  2)=        1.1772110462_dp
    xab_pm6( 18,  3)=        3.3629100323_dp
    xab_pm6( 18,  4)=        2.7554919720_dp
    xab_pm6( 18,  5)=        1.9315860271_dp
    xab_pm6( 18,  6)=        0.1223089993_dp
    xab_pm6( 18,  7)=        0.5625810027_dp
    xab_pm6( 18,  8)=        0.3557949960_dp
    xab_pm6( 18,  9)=        9.2697153091_dp
    xab_pm6( 18, 10)=        1.3046970367_dp
    xab_pm6( 18, 11)=        3.3981380463_dp
    xab_pm6( 18, 12)=        1.9706380367_dp
    xab_pm6( 18, 13)=        1.8520090580_dp
    xab_pm6( 18, 14)=        1.0678900480_dp
    xab_pm6( 18, 15)=        0.1710139960_dp
    xab_pm6( 18, 16)=        0.6537690163_dp
    xab_pm6( 18, 17)=        2.2560939789_dp
    xab_pm6( 18, 18)=        0.9726989865_dp
    xab_pm6( 19,  1)=        0.3693400025_dp
    xab_pm6( 19,  2)=        2.8950450420_dp
    xab_pm6( 19,  3)=        4.3745670319_dp
    xab_pm6( 19,  4)=        6.7322211266_dp
    xab_pm6( 19,  5)=        8.9005413055_dp
    xab_pm6( 19,  6)=       10.3179874420_dp
    xab_pm6( 19,  7)=        7.9777069092_dp
    xab_pm6( 19,  8)=        1.0122749805_dp
    xab_pm6( 19,  9)=        6.5929708481_dp
    xab_pm6( 19, 10)=        0.2339950055_dp
    xab_pm6( 19, 11)=        5.5630269051_dp
    xab_pm6( 19, 12)=        3.2905020714_dp
    xab_pm6( 19, 13)=       29.9447078705_dp
    xab_pm6( 19, 14)=        8.2791996002_dp
    xab_pm6( 19, 15)=        4.4753842354_dp
    xab_pm6( 19, 16)=       29.5289516449_dp
    xab_pm6( 19, 17)=        1.2314809561_dp
    xab_pm6( 19, 18)=        9.7105083466_dp
    xab_pm6( 19, 19)=        5.9343290329_dp
    xab_pm6( 20,  1)=        7.7286062241_dp
    xab_pm6( 20,  2)=        2.9138519764_dp
    xab_pm6( 20,  5)=        1.7000099421_dp
    xab_pm6( 20,  6)=        0.1484500021_dp
    xab_pm6( 20,  7)=        2.9880740643_dp
    xab_pm6( 20,  8)=       17.0289459229_dp
    xab_pm6( 20,  9)=        3.4828209877_dp
    xab_pm6( 20, 10)=        0.3325859904_dp
    xab_pm6( 20, 11)=        9.6575088501_dp
    xab_pm6( 20, 12)=        8.5998001099_dp
    xab_pm6( 20, 13)=        4.1885547638_dp
    xab_pm6( 20, 14)=        0.3362329900_dp
    xab_pm6( 20, 15)=        0.4108400047_dp
    xab_pm6( 20, 16)=        0.3257389963_dp
    xab_pm6( 20, 17)=        5.9561438560_dp
    xab_pm6( 20, 18)=        0.2910720110_dp
    xab_pm6( 20, 19)=        1.2403199673_dp
    xab_pm6( 20, 20)=       30.0035915375_dp
    xab_pm6( 21,  1)=        0.3511990011_dp
    xab_pm6( 21,  6)=        8.6080522537_dp
    xab_pm6( 21,  7)=        3.2318809032_dp
    xab_pm6( 21,  8)=        3.0586719513_dp
    xab_pm6( 21,  9)=        7.2523469925_dp
    xab_pm6( 21, 13)=        0.5006200075_dp
    xab_pm6( 21, 14)=        3.2190699577_dp
    xab_pm6( 21, 15)=        0.6267489791_dp
    xab_pm6( 21, 16)=        0.2118500024_dp
    xab_pm6( 21, 17)=        2.9961290359_dp
    xab_pm6( 21, 21)=        2.5981659889_dp
    xab_pm6( 22,  1)=        0.1437219977_dp
    xab_pm6( 22,  5)=        0.6493600011_dp
    xab_pm6( 22,  6)=        0.4167059958_dp
    xab_pm6( 22,  7)=        0.5454609990_dp
    xab_pm6( 22,  8)=        0.7994859815_dp
    xab_pm6( 22,  9)=        1.0857419968_dp
    xab_pm6( 22, 12)=        4.3302397728_dp
    xab_pm6( 22, 13)=        2.0918409824_dp
    xab_pm6( 22, 14)=        6.7738151550_dp
    xab_pm6( 22, 15)=        4.1504998207_dp
    xab_pm6( 22, 16)=        0.9437839985_dp
    xab_pm6( 22, 17)=        0.3332970142_dp
    xab_pm6( 22, 20)=        4.1091408730_dp
    xab_pm6( 22, 22)=        2.0000000000_dp
    xab_pm6( 23,  1)=        0.1052040011_dp
    xab_pm6( 23,  6)=        1.9387600422_dp
    xab_pm6( 23,  7)=        0.2767249942_dp
    xab_pm6( 23,  8)=        0.4153119922_dp
    xab_pm6( 23,  9)=        0.3428150117_dp
    xab_pm6( 23, 11)=        8.2760200500_dp
    xab_pm6( 23, 15)=        6.2506241798_dp
    xab_pm6( 23, 16)=        2.0350389481_dp
    xab_pm6( 23, 17)=        0.2436569929_dp
    xab_pm6( 23, 19)=        2.0265901089_dp
    xab_pm6( 23, 23)=       10.7798919678_dp
    xab_pm6( 24,  1)=        0.0444689989_dp
    xab_pm6( 24,  6)=        6.1101870537_dp
    xab_pm6( 24,  7)=        1.9203239679_dp
    xab_pm6( 24,  8)=        1.0555109978_dp
    xab_pm6( 24,  9)=        0.7376070023_dp
    xab_pm6( 24, 11)=        8.3642740250_dp
    xab_pm6( 24, 14)=        1.0291099548_dp
    xab_pm6( 24, 15)=        0.6001769900_dp
    xab_pm6( 24, 16)=        0.5503339767_dp
    xab_pm6( 24, 17)=        0.3690730035_dp
    xab_pm6( 24, 19)=        2.0000000000_dp
    xab_pm6( 24, 24)=       10.3186073303_dp
    xab_pm6( 25,  1)=        1.2692099810_dp
    xab_pm6( 25,  6)=        2.5831100941_dp
    xab_pm6( 25,  7)=        1.9567500353_dp
    xab_pm6( 25,  8)=        1.2856199741_dp
    xab_pm6( 25,  9)=        1.1130700111_dp
    xab_pm6( 25, 13)=        1.0407899618_dp
    xab_pm6( 25, 14)=        0.9505800009_dp
    xab_pm6( 25, 15)=        1.1303199530_dp
    xab_pm6( 25, 16)=        1.6126500368_dp
    xab_pm6( 25, 17)=        0.2018499970_dp
    xab_pm6( 25, 20)=        0.6201800108_dp
    xab_pm6( 25, 25)=        2.4600400925_dp
    xab_pm6( 26,  1)=        0.0251950007_dp
    xab_pm6( 26,  6)=        0.3668349981_dp
    xab_pm6( 26,  7)=        0.1553419977_dp
    xab_pm6( 26,  8)=        0.1364219934_dp
    xab_pm6( 26,  9)=        3.6573500633_dp
    xab_pm6( 26, 15)=        0.4312910140_dp
    xab_pm6( 26, 16)=        0.0334779993_dp
    xab_pm6( 26, 17)=        0.0194729995_dp
    xab_pm6( 26, 19)=        6.0000000000_dp
    xab_pm6( 26, 26)=        1.8468899727_dp
    xab_pm6( 27,  1)=        2.4724650383_dp
    xab_pm6( 27,  6)=        2.1239299774_dp
    xab_pm6( 27,  7)=        2.6538360119_dp
    xab_pm6( 27,  8)=        5.2520217896_dp
    xab_pm6( 27,  9)=        4.5850300789_dp
    xab_pm6( 27, 14)=        1.0902400017_dp
    xab_pm6( 27, 15)=        0.1059359983_dp
    xab_pm6( 27, 16)=        0.4367069900_dp
    xab_pm6( 27, 17)=        1.0334140062_dp
    xab_pm6( 27, 27)=        3.9196178913_dp
    xab_pm6( 28,  1)=        1.7631239891_dp
    xab_pm6( 28,  6)=        7.1333241463_dp
    xab_pm6( 28,  7)=        4.2867999077_dp
    xab_pm6( 28,  8)=        0.8859419823_dp
    xab_pm6( 28,  9)=        1.0882079601_dp
    xab_pm6( 28, 14)=        0.9386460185_dp
    xab_pm6( 28, 15)=        5.0597271919_dp
    xab_pm6( 28, 16)=        0.2748520076_dp
    xab_pm6( 28, 17)=        0.2023130059_dp
    xab_pm6( 28, 28)=        0.0354739986_dp
    xab_pm6( 29,  1)=        0.6035910249_dp
    xab_pm6( 29,  6)=        7.0677938461_dp
    xab_pm6( 29,  7)=        3.2286670208_dp
    xab_pm6( 29,  8)=        2.0000000000_dp
    xab_pm6( 29,  9)=        1.2821079493_dp
    xab_pm6( 29, 15)=        0.0201260000_dp
    xab_pm6( 29, 16)=        0.0052479999_dp
    xab_pm6( 29, 17)=        0.1390649974_dp
    xab_pm6( 29, 29)=        5.1843757629_dp
    xab_pm6( 30,  1)=        3.1091930866_dp
    xab_pm6( 30,  6)=        0.9914649725_dp
    xab_pm6( 30,  7)=        0.9524760246_dp
    xab_pm6( 30,  8)=        2.2653129101_dp
    xab_pm6( 30,  9)=        1.2255450487_dp
    xab_pm6( 30, 14)=        3.7839050293_dp
    xab_pm6( 30, 15)=        0.5815299749_dp
    xab_pm6( 30, 16)=        0.6480000019_dp
    xab_pm6( 30, 17)=        0.7213510275_dp
    xab_pm6( 30, 20)=        1.2402900457_dp
    xab_pm6( 30, 30)=        0.4650000036_dp
    xab_pm6( 31,  1)=        1.3866519928_dp
    xab_pm6( 31,  6)=        1.9629900455_dp
    xab_pm6( 31,  7)=        1.1883380413_dp
    xab_pm6( 31,  8)=        1.5236439705_dp
    xab_pm6( 31,  9)=        1.4169420004_dp
    xab_pm6( 31, 14)=        1.0022900105_dp
    xab_pm6( 31, 15)=        0.5000000000_dp
    xab_pm6( 31, 16)=        2.4562840462_dp
    xab_pm6( 31, 17)=        1.1866610050_dp
    xab_pm6( 31, 31)=        1.1983939409_dp
    xab_pm6( 32,  1)=        1.7332259417_dp
    xab_pm6( 32,  6)=        1.2975100279_dp
    xab_pm6( 32,  7)=        0.6375060081_dp
    xab_pm6( 32,  8)=        0.8269640207_dp
    xab_pm6( 32,  9)=        0.6519770026_dp
    xab_pm6( 32, 14)=        0.1786800027_dp
    xab_pm6( 32, 15)=        5.6163492203_dp
    xab_pm6( 32, 16)=        1.1609569788_dp
    xab_pm6( 32, 17)=        0.5452389717_dp
    xab_pm6( 32, 25)=        2.2551510334_dp
    xab_pm6( 32, 27)=        2.1518499851_dp
    xab_pm6( 32, 32)=        3.0230000019_dp
    xab_pm6( 33,  1)=        1.0905890465_dp
    xab_pm6( 33,  6)=        0.5790979862_dp
    xab_pm6( 33,  7)=        0.2733370066_dp
    xab_pm6( 33,  8)=        0.7016140223_dp
    xab_pm6( 33,  9)=        0.4026280046_dp
    xab_pm6( 33, 13)=        1.0035799742_dp
    xab_pm6( 33, 14)=        1.4307060242_dp
    xab_pm6( 33, 16)=        1.0337840319_dp
    xab_pm6( 33, 17)=        0.4544329941_dp
    xab_pm6( 33, 22)=        1.5813169479_dp
    xab_pm6( 33, 27)=        1.6752400398_dp
    xab_pm6( 33, 30)=        3.1565709114_dp
    xab_pm6( 33, 31)=        1.6862980127_dp
    xab_pm6( 33, 33)=        0.7373070121_dp
    xab_pm6( 34,  1)=        0.8479980230_dp
    xab_pm6( 34,  6)=        1.1147869825_dp
    xab_pm6( 34,  7)=        0.4828400016_dp
    xab_pm6( 34,  8)=        0.9551609755_dp
    xab_pm6( 34,  9)=        0.4448060095_dp
    xab_pm6( 34, 14)=        0.5182269812_dp
    xab_pm6( 34, 15)=        0.2920520008_dp
    xab_pm6( 34, 16)=        0.3917210102_dp
    xab_pm6( 34, 17)=        0.9810670018_dp
    xab_pm6( 34, 25)=        2.1807200909_dp
    xab_pm6( 34, 27)=        2.2024099827_dp
    xab_pm6( 34, 30)=        0.5115939975_dp
    xab_pm6( 34, 32)=        5.8720507622_dp
    xab_pm6( 34, 33)=        0.7112609744_dp
    xab_pm6( 34, 34)=        0.8218230009_dp
    xab_pm6( 35,  1)=        0.8503779769_dp
    xab_pm6( 35,  2)=        1.0620429516_dp
    xab_pm6( 35,  3)=        1.8588659763_dp
    xab_pm6( 35,  4)=        1.9409329891_dp
    xab_pm6( 35,  5)=        1.2264200449_dp
    xab_pm6( 35,  6)=        0.5706859827_dp
    xab_pm6( 35,  7)=       30.0001335144_dp
    xab_pm6( 35,  8)=        0.7065839767_dp
    xab_pm6( 35,  9)=        0.2935000062_dp
    xab_pm6( 35, 10)=        1.0061589479_dp
    xab_pm6( 35, 11)=        1.7529369593_dp
    xab_pm6( 35, 12)=        2.9162800312_dp
    xab_pm6( 35, 13)=        2.3571300507_dp
    xab_pm6( 35, 14)=        0.5895109773_dp
    xab_pm6( 35, 15)=        0.4565210044_dp
    xab_pm6( 35, 16)=        0.2866880000_dp
    xab_pm6( 35, 17)=        0.3892379999_dp
    xab_pm6( 35, 18)=        3.2626678944_dp
    xab_pm6( 35, 19)=        3.3227949142_dp
    xab_pm6( 35, 20)=        4.0529098511_dp
    xab_pm6( 35, 21)=        2.0982511044_dp
    xab_pm6( 35, 22)=        0.8834339976_dp
    xab_pm6( 35, 23)=        0.6126980186_dp
    xab_pm6( 35, 24)=        0.2178529948_dp
    xab_pm6( 35, 25)=        1.1835800409_dp
    xab_pm6( 35, 26)=        6.0619211197_dp
    xab_pm6( 35, 27)=        0.4251480103_dp
    xab_pm6( 35, 28)=        0.6321449876_dp
    xab_pm6( 35, 29)=        0.7685170174_dp
    xab_pm6( 35, 30)=        0.7470269799_dp
    xab_pm6( 35, 31)=        1.2610360384_dp
    xab_pm6( 35, 32)=        0.6277369857_dp
    xab_pm6( 35, 33)=        0.5141530037_dp
    xab_pm6( 35, 34)=        0.3193419874_dp
    xab_pm6( 35, 35)=        0.6153079867_dp
    xab_pm6( 36,  1)=        5.1258969307_dp
    xab_pm6( 36,  2)=        0.6277009845_dp
    xab_pm6( 36,  3)=        8.7586965561_dp
    xab_pm6( 36,  4)=       10.2377958298_dp
    xab_pm6( 36,  5)=        2.9467809200_dp
    xab_pm6( 36,  6)=        0.6526229978_dp
    xab_pm6( 36,  7)=        0.1996060014_dp
    xab_pm6( 36,  8)=        0.0067329998_dp
    xab_pm6( 36,  9)=        4.1344070435_dp
    xab_pm6( 36, 10)=        1.4337220192_dp
    xab_pm6( 36, 11)=        8.3544483185_dp
    xab_pm6( 36, 12)=        0.8884360194_dp
    xab_pm6( 36, 13)=        5.0917158127_dp
    xab_pm6( 36, 14)=        0.5542500019_dp
    xab_pm6( 36, 17)=        0.5202170014_dp
    xab_pm6( 36, 18)=        0.5548740029_dp
    xab_pm6( 36, 19)=        8.6097822189_dp
    xab_pm6( 36, 20)=        0.8788909912_dp
    xab_pm6( 36, 35)=        0.3080979884_dp
    xab_pm6( 36, 36)=        0.0520990007_dp
    xab_pm6( 37,  1)=       29.8616313934_dp
    xab_pm6( 37,  2)=        1.8625849485_dp
    xab_pm6( 37,  5)=        9.0404930115_dp
    xab_pm6( 37,  6)=       29.9740314484_dp
    xab_pm6( 37,  7)=        0.0246360004_dp
    xab_pm6( 37,  8)=        1.1253499985_dp
    xab_pm6( 37,  9)=       28.8152770996_dp
    xab_pm6( 37, 10)=        7.7365632057_dp
    xab_pm6( 37, 13)=        2.9924569130_dp
    xab_pm6( 37, 16)=        0.9644110203_dp
    xab_pm6( 37, 17)=       10.3844861984_dp
    xab_pm6( 37, 18)=       18.4333286285_dp
    xab_pm6( 37, 35)=        5.1762142181_dp
    xab_pm6( 37, 36)=       15.3075027466_dp
    xab_pm6( 37, 37)=       20.1476097107_dp
    xab_pm6( 38,  1)=       12.9733161926_dp
    xab_pm6( 38,  6)=        6.6546568871_dp
    xab_pm6( 38,  7)=        6.8538661003_dp
    xab_pm6( 38,  8)=        3.5613958836_dp
    xab_pm6( 38,  9)=       10.9717054367_dp
    xab_pm6( 38, 14)=        2.7647500038_dp
    xab_pm6( 38, 15)=        2.5520999432_dp
    xab_pm6( 38, 16)=        3.1296029091_dp
    xab_pm6( 38, 17)=        3.7839550972_dp
    xab_pm6( 38, 22)=        2.8172500134_dp
    xab_pm6( 38, 35)=        2.7665669918_dp
    xab_pm6( 38, 38)=        5.3721199036_dp
    xab_pm6( 39,  1)=        0.6123989820_dp
    xab_pm6( 39,  6)=        0.5043060184_dp
    xab_pm6( 39,  7)=        1.6279029846_dp
    xab_pm6( 39,  8)=        1.7429219484_dp
    xab_pm6( 39,  9)=        4.4338088036_dp
    xab_pm6( 39, 13)=        0.5006700158_dp
    xab_pm6( 39, 14)=        3.2190299034_dp
    xab_pm6( 39, 15)=        0.5416600108_dp
    xab_pm6( 39, 16)=        0.3182219863_dp
    xab_pm6( 39, 17)=        1.1549589634_dp
    xab_pm6( 39, 35)=        1.0543160439_dp
    xab_pm6( 39, 39)=        1.6917250156_dp
    xab_pm6( 40,  1)=        0.5937319994_dp
    xab_pm6( 40,  6)=        1.9991819859_dp
    xab_pm6( 40,  7)=        0.9950450063_dp
    xab_pm6( 40,  8)=        1.0575250387_dp
    xab_pm6( 40,  9)=        0.8611419797_dp
    xab_pm6( 40, 13)=        0.8740599751_dp
    xab_pm6( 40, 14)=        1.7233430147_dp
    xab_pm6( 40, 15)=        0.7483760118_dp
    xab_pm6( 40, 16)=        2.4293239117_dp
    xab_pm6( 40, 17)=        0.4430989921_dp
    xab_pm6( 40, 35)=        0.8589090109_dp
    xab_pm6( 40, 40)=        3.0777730942_dp
    xab_pm6( 41,  1)=        3.6037790775_dp
    xab_pm6( 41,  6)=        4.5754809380_dp
    xab_pm6( 41,  7)=        1.2135870457_dp
    xab_pm6( 41,  8)=        1.1847189665_dp
    xab_pm6( 41,  9)=        3.6636819839_dp
    xab_pm6( 41, 11)=        8.2760200500_dp
    xab_pm6( 41, 15)=        6.2015070915_dp
    xab_pm6( 41, 16)=        2.4600200653_dp
    xab_pm6( 41, 17)=        1.8915569782_dp
    xab_pm6( 41, 19)=        2.0265901089_dp
    xab_pm6( 41, 35)=        1.9212690592_dp
    xab_pm6( 41, 41)=        2.1223878860_dp
    xab_pm6( 42,  1)=        0.9346860051_dp
    xab_pm6( 42,  6)=        1.1907420158_dp
    xab_pm6( 42,  7)=        0.6082680225_dp
    xab_pm6( 42,  8)=        0.5112670064_dp
    xab_pm6( 42,  9)=        0.6104289889_dp
    xab_pm6( 42, 11)=        8.2865495682_dp
    xab_pm6( 42, 15)=        1.5228459835_dp
    xab_pm6( 42, 16)=        0.8304280043_dp
    xab_pm6( 42, 17)=        0.4743250012_dp
    xab_pm6( 42, 19)=        2.1423900127_dp
    xab_pm6( 42, 24)=        1.7419430017_dp
    xab_pm6( 42, 35)=        0.2259179950_dp
    xab_pm6( 42, 42)=        0.6264619827_dp
    xab_pm6( 43,  1)=        6.3103342056_dp
    xab_pm6( 43,  6)=        3.9724390507_dp
    xab_pm6( 43,  7)=        0.7271299958_dp
    xab_pm6( 43,  8)=        1.0246160030_dp
    xab_pm6( 43,  9)=        5.8117837906_dp
    xab_pm6( 43, 16)=        1.4965020418_dp
    xab_pm6( 43, 17)=        1.6515829563_dp
    xab_pm6( 43, 32)=        2.1520600319_dp
    xab_pm6( 43, 34)=        2.2026200294_dp
    xab_pm6( 43, 35)=        3.8201301098_dp
    xab_pm6( 44,  1)=        7.1379761696_dp
    xab_pm6( 44,  6)=        1.1349359751_dp
    xab_pm6( 44,  7)=        2.3340940475_dp
    xab_pm6( 44,  8)=        2.9762790203_dp
    xab_pm6( 44,  9)=        6.9471278191_dp
    xab_pm6( 44, 14)=        0.8494300246_dp
    xab_pm6( 44, 15)=        0.0569740012_dp
    xab_pm6( 44, 16)=        1.0066829920_dp
    xab_pm6( 44, 17)=        0.1265860051_dp
    xab_pm6( 44, 32)=        2.1515600681_dp
    xab_pm6( 44, 34)=        2.2021200657_dp
    xab_pm6( 44, 35)=        0.6598809958_dp
    xab_pm6( 44, 44)=        0.0978050008_dp
    xab_pm6( 45,  1)=        2.3061070442_dp
    xab_pm6( 45,  6)=        3.4880790710_dp
    xab_pm6( 45,  7)=        4.0009469986_dp
    xab_pm6( 45,  8)=       10.2986755371_dp
    xab_pm6( 45,  9)=        9.0653839111_dp
    xab_pm6( 45, 14)=        0.8500099778_dp
    xab_pm6( 45, 15)=        1.0381410122_dp
    xab_pm6( 45, 16)=        4.8164100647_dp
    xab_pm6( 45, 17)=        3.5868649483_dp
    xab_pm6( 45, 32)=        2.1521399021_dp
    xab_pm6( 45, 34)=        2.2026998997_dp
    xab_pm6( 45, 35)=        1.5101490021_dp
    xab_pm6( 45, 45)=        2.0701138973_dp
    xab_pm6( 46,  1)=        0.4432690144_dp
    xab_pm6( 46,  6)=        9.8537149429_dp
    xab_pm6( 46,  7)=        0.2497030050_dp
    xab_pm6( 46,  8)=        0.2164029926_dp
    xab_pm6( 46,  9)=        6.9453120232_dp
    xab_pm6( 46, 13)=        1.0572899580_dp
    xab_pm6( 46, 14)=        2.2251040936_dp
    xab_pm6( 46, 15)=        0.0450170003_dp
    xab_pm6( 46, 16)=        0.2552289963_dp
    xab_pm6( 46, 17)=        2.9698910713_dp
    xab_pm6( 46, 35)=        4.6380510330_dp
    xab_pm6( 46, 46)=        0.0519560017_dp
    xab_pm6( 47,  1)=        1.9951679707_dp
    xab_pm6( 47,  6)=       11.3354558945_dp
    xab_pm6( 47,  7)=       19.8037109375_dp
    xab_pm6( 47,  8)=        0.1656610072_dp
    xab_pm6( 47,  9)=       12.6958837509_dp
    xab_pm6( 47, 13)=        0.8965139985_dp
    xab_pm6( 47, 15)=        0.0499319993_dp
    xab_pm6( 47, 16)=       11.1880216599_dp
    xab_pm6( 47, 17)=       23.7654590607_dp
    xab_pm6( 47, 35)=        1.7143690586_dp
    xab_pm6( 47, 47)=        0.5577419996_dp
    xab_pm6( 48,  1)=       11.9142007828_dp
    xab_pm6( 48,  6)=        0.6034410000_dp
    xab_pm6( 48,  7)=        0.1806630045_dp
    xab_pm6( 48,  8)=        0.9261459708_dp
    xab_pm6( 48,  9)=        1.3536649942_dp
    xab_pm6( 48, 14)=        2.2533459663_dp
    xab_pm6( 48, 16)=        0.3613890111_dp
    xab_pm6( 48, 17)=        0.1404239982_dp
    xab_pm6( 48, 35)=        0.2722670138_dp
    xab_pm6( 48, 48)=       18.6179981232_dp
    xab_pm6( 49,  1)=       14.9752931595_dp
    xab_pm6( 49,  6)=        2.1873850822_dp
    xab_pm6( 49,  7)=        3.3699929714_dp
    xab_pm6( 49,  8)=        4.1285829544_dp
    xab_pm6( 49,  9)=        3.7010159492_dp
    xab_pm6( 49, 16)=        6.3411049843_dp
    xab_pm6( 49, 17)=        2.3885519505_dp
    xab_pm6( 49, 31)=        2.4219870567_dp
    xab_pm6( 49, 33)=        6.2083501816_dp
    xab_pm6( 49, 34)=        2.3193230629_dp
    xab_pm6( 49, 35)=        3.7285981178_dp
    xab_pm6( 49, 49)=        8.0634908676_dp
    xab_pm6( 50,  1)=        6.5351619720_dp
    xab_pm6( 50,  6)=        3.3743550777_dp
    xab_pm6( 50,  7)=        1.3918999434_dp
    xab_pm6( 50,  8)=        4.3740167618_dp
    xab_pm6( 50,  9)=       18.5986633301_dp
    xab_pm6( 50, 16)=        2.3148701191_dp
    xab_pm6( 50, 17)=        0.8077819943_dp
    xab_pm6( 50, 32)=       12.3434114456_dp
    xab_pm6( 50, 34)=        3.0618851185_dp
    xab_pm6( 50, 35)=        0.6687980294_dp
    xab_pm6( 50, 50)=        0.2870000005_dp
    xab_pm6( 51,  1)=        0.7953429818_dp
    xab_pm6( 51,  6)=        0.5792120099_dp
    xab_pm6( 51,  7)=        0.0820650011_dp
    xab_pm6( 51,  8)=        0.6342340112_dp
    xab_pm6( 51,  9)=        0.6502770185_dp
    xab_pm6( 51, 13)=        1.6166900396_dp
    xab_pm6( 51, 14)=        8.7137489319_dp
    xab_pm6( 51, 16)=        0.3969689906_dp
    xab_pm6( 51, 17)=        0.1564749926_dp
    xab_pm6( 51, 25)=        2.2367100716_dp
    xab_pm6( 51, 27)=        2.2760500908_dp
    xab_pm6( 51, 35)=        0.1980440021_dp
    xab_pm6( 51, 43)=        2.2762598991_dp
    xab_pm6( 51, 44)=        2.2757599354_dp
    xab_pm6( 51, 45)=        2.2763400078_dp
    xab_pm6( 51, 49)=        6.6608009338_dp
    xab_pm6( 51, 51)=        0.7248849869_dp
    xab_pm6( 52,  1)=        1.8076790571_dp
    xab_pm6( 52,  6)=        0.9704939723_dp
    xab_pm6( 52,  7)=        0.3585929871_dp
    xab_pm6( 52,  8)=        0.3829259872_dp
    xab_pm6( 52,  9)=        0.2008219957_dp
    xab_pm6( 52, 13)=        2.1068120003_dp
    xab_pm6( 52, 15)=        1.1092890501_dp
    xab_pm6( 52, 16)=        0.9439250231_dp
    xab_pm6( 52, 17)=        0.2854779959_dp
    xab_pm6( 52, 30)=        1.7560700178_dp
    xab_pm6( 52, 32)=        7.0190491676_dp
    xab_pm6( 52, 33)=        0.6857740283_dp
    xab_pm6( 52, 34)=        1.1878260374_dp
    xab_pm6( 52, 35)=        0.3942019939_dp
    xab_pm6( 52, 48)=        1.0859190226_dp
    xab_pm6( 52, 49)=        2.0395820141_dp
    xab_pm6( 52, 50)=        2.9519760609_dp
    xab_pm6( 52, 52)=        0.6424859762_dp
    xab_pm6( 53,  1)=        0.9818980098_dp
    xab_pm6( 53,  2)=        1.6307209730_dp
    xab_pm6( 53,  3)=        4.1685991287_dp
    xab_pm6( 53,  4)=        2.3518979549_dp
    xab_pm6( 53,  5)=        3.1613850594_dp
    xab_pm6( 53,  6)=        0.8101559877_dp
    xab_pm6( 53,  7)=        0.2649030089_dp
    xab_pm6( 53,  8)=        0.8662040234_dp
    xab_pm6( 53,  9)=        0.3924250007_dp
    xab_pm6( 53, 10)=        1.5035680532_dp
    xab_pm6( 53, 11)=        1.9861119986_dp
    xab_pm6( 53, 12)=        3.2769138813_dp
    xab_pm6( 53, 13)=        2.9290800095_dp
    xab_pm6( 53, 14)=        0.7002990246_dp
    xab_pm6( 53, 15)=        3.0472071171_dp
    xab_pm6( 53, 16)=        0.7099289894_dp
    xab_pm6( 53, 17)=        0.3104740083_dp
    xab_pm6( 53, 18)=        0.3053669930_dp
    xab_pm6( 53, 19)=        4.8243532181_dp
    xab_pm6( 53, 20)=        7.6899209023_dp
    xab_pm6( 53, 21)=        3.1142818928_dp
    xab_pm6( 53, 22)=        2.4267470837_dp
    xab_pm6( 53, 23)=        6.1981120110_dp
    xab_pm6( 53, 24)=        2.5985898972_dp
    xab_pm6( 53, 25)=        1.1934100389_dp
    xab_pm6( 53, 26)=        0.5326219797_dp
    xab_pm6( 53, 27)=        1.1052390337_dp
    xab_pm6( 53, 28)=        0.0174589995_dp
    xab_pm6( 53, 29)=        0.0067810002_dp
    xab_pm6( 53, 30)=        0.9766070247_dp
    xab_pm6( 53, 31)=        1.2521680593_dp
    xab_pm6( 53, 32)=        1.3232669830_dp
    xab_pm6( 53, 33)=        0.3108240068_dp
    xab_pm6( 53, 35)=        0.4830540121_dp
    xab_pm6( 53, 36)=        0.2011359930_dp
    xab_pm6( 53, 37)=        4.0924458504_dp
    xab_pm6( 53, 38)=        2.1039409637_dp
    xab_pm6( 53, 39)=        1.0214020014_dp
    xab_pm6( 53, 40)=        4.5139431953_dp
    xab_pm6( 53, 41)=        2.3992979527_dp
    xab_pm6( 53, 42)=        0.1246950030_dp
    xab_pm6( 53, 43)=        0.1105939969_dp
    xab_pm6( 53, 44)=        7.8377099037_dp
    xab_pm6( 53, 45)=        2.3579440117_dp
    xab_pm6( 53, 46)=        3.5220170021_dp
    xab_pm6( 53, 47)=        0.0489040017_dp
    xab_pm6( 53, 48)=        0.3967840075_dp
    xab_pm6( 53, 49)=        5.9478211403_dp
    xab_pm6( 53, 50)=        1.7831629515_dp
    xab_pm6( 53, 51)=        0.3181900084_dp
    xab_pm6( 53, 52)=        1.1011159420_dp
    xab_pm6( 53, 53)=        0.5105419755_dp
    xab_pm6( 54,  1)=        0.7010160089_dp
    xab_pm6( 54,  2)=        2.5994710922_dp
    xab_pm6( 54,  3)=        4.5820808411_dp
    xab_pm6( 54,  4)=        0.6605250239_dp
    xab_pm6( 54,  5)=        1.1006120443_dp
    xab_pm6( 54,  6)=        0.8267269731_dp
    xab_pm6( 54,  7)=        0.9256240129_dp
    xab_pm6( 54,  8)=        0.0353560001_dp
    xab_pm6( 54,  9)=        0.0650110021_dp
    xab_pm6( 54, 10)=        0.2938619852_dp
    xab_pm6( 54, 11)=        8.3682041168_dp
    xab_pm6( 54, 12)=        9.7235717773_dp
    xab_pm6( 54, 13)=        7.4044651985_dp
    xab_pm6( 54, 14)=       16.0919990540_dp
    xab_pm6( 54, 17)=        0.4637579918_dp
    xab_pm6( 54, 18)=        0.0492659993_dp
    xab_pm6( 54, 19)=        1.2248890400_dp
    xab_pm6( 54, 20)=        1.7171210051_dp
    xab_pm6( 54, 35)=        0.4751160145_dp
    xab_pm6( 54, 36)=        0.0497930013_dp
    xab_pm6( 54, 37)=        0.9749649763_dp
    xab_pm6( 54, 53)=        0.1120899990_dp
    xab_pm6( 54, 54)=        0.3444739878_dp
    xab_pm6( 55,  1)=        0.0969009995_dp
    xab_pm6( 55,  5)=       10.3926095963_dp
    xab_pm6( 55,  6)=       24.5146236420_dp
    xab_pm6( 55,  7)=       29.7110767365_dp
    xab_pm6( 55,  8)=        8.1768426895_dp
    xab_pm6( 55,  9)=       22.2334156036_dp
    xab_pm6( 55, 15)=        0.5065119863_dp
    xab_pm6( 55, 16)=        0.0917429999_dp
    xab_pm6( 55, 17)=        4.5319647789_dp
    xab_pm6( 55, 35)=        1.6584270000_dp
    xab_pm6( 55, 53)=        1.0721780062_dp
    xab_pm6( 55, 55)=       25.3200550079_dp
    xab_pm6( 56,  1)=        2.0400040150_dp
    xab_pm6( 56,  6)=        0.1197929978_dp
    xab_pm6( 56,  7)=        0.2079340070_dp
    xab_pm6( 56,  8)=        0.3489449918_dp
    xab_pm6( 56,  9)=        5.5752549171_dp
    xab_pm6( 56, 13)=        9.5390987396_dp
    xab_pm6( 56, 14)=        1.2126599550_dp
    xab_pm6( 56, 16)=        0.2153860033_dp
    xab_pm6( 56, 17)=        0.1601770073_dp
    xab_pm6( 56, 22)=        9.4935302734_dp
    xab_pm6( 56, 35)=        0.8287940025_dp
    xab_pm6( 56, 53)=        0.8355969787_dp
    xab_pm6( 56, 56)=        0.3561860025_dp
    xab_pm6( 57,  1)=        0.6235010028_dp
    xab_pm6( 57,  6)=        0.1086490005_dp
    xab_pm6( 57,  7)=        0.1047779992_dp
    xab_pm6( 57,  8)=        0.5579569936_dp
    xab_pm6( 57,  9)=        2.4019029140_dp
    xab_pm6( 57, 13)=        0.5005400181_dp
    xab_pm6( 57, 14)=        3.2190299034_dp
    xab_pm6( 57, 15)=        0.5416600108_dp
    xab_pm6( 57, 16)=        2.6824119091_dp
    xab_pm6( 57, 17)=        0.2302030027_dp
    xab_pm6( 57, 35)=        0.2385820001_dp
    xab_pm6( 57, 53)=        0.2268829942_dp
    xab_pm6( 57, 57)=        5.1751618385_dp
    xab_pm6( 64,  1)=        0.1358100027_dp
    xab_pm6( 64,  6)=        0.0530400015_dp
    xab_pm6( 64,  7)=        0.2050500065_dp
    xab_pm6( 64,  8)=        0.1757999957_dp
    xab_pm6( 64,  9)=        0.3346300125_dp
    xab_pm6( 64, 13)=        0.5005400181_dp
    xab_pm6( 64, 14)=        3.2190299034_dp
    xab_pm6( 64, 15)=        0.5416600108_dp
    xab_pm6( 64, 16)=        2.6554000378_dp
    xab_pm6( 64, 17)=        0.0899700001_dp
    xab_pm6( 64, 35)=        0.2407400012_dp
    xab_pm6( 64, 53)=        0.2782399952_dp
    xab_pm6( 64, 64)=        2.6703999043_dp
    xab_pm6( 71,  1)=        0.7879199982_dp
    xab_pm6( 71,  6)=        4.4538249969_dp
    xab_pm6( 71,  7)=        2.8608279228_dp
    xab_pm6( 71,  8)=        2.9170761108_dp
    xab_pm6( 71, 15)=        0.5000000000_dp
    xab_pm6( 71, 17)=       12.7570991516_dp
    xab_pm6( 71, 35)=        8.6482744217_dp
    xab_pm6( 71, 53)=       10.0823154449_dp
    xab_pm6( 72,  1)=        3.4273118973_dp
    xab_pm6( 72,  5)=        0.6592699885_dp
    xab_pm6( 72,  6)=        0.3785789907_dp
    xab_pm6( 72,  7)=        0.6557949781_dp
    xab_pm6( 72,  8)=        1.0347180367_dp
    xab_pm6( 72,  9)=        1.6793349981_dp
    xab_pm6( 72, 12)=        4.3302497864_dp
    xab_pm6( 72, 13)=        0.6225200295_dp
    xab_pm6( 72, 14)=        3.3822999001_dp
    xab_pm6( 72, 15)=        0.5055299997_dp
    xab_pm6( 72, 16)=        1.6667599678_dp
    xab_pm6( 72, 17)=        0.7064210176_dp
    xab_pm6( 72, 20)=        4.3195099831_dp
    xab_pm6( 72, 33)=        1.2808200121_dp
    xab_pm6( 72, 35)=        0.6924560070_dp
    xab_pm6( 72, 53)=        0.8209480047_dp
    xab_pm6( 72, 56)=        9.0225200653_dp
    xab_pm6( 72, 72)=        1.0589109659_dp
    xab_pm6( 73,  1)=        2.8276689053_dp
    xab_pm6( 73,  6)=        0.8474389911_dp
    xab_pm6( 73,  7)=        1.0154609680_dp
    xab_pm6( 73,  8)=        1.7510830164_dp
    xab_pm6( 73,  9)=        3.1465198994_dp
    xab_pm6( 73, 11)=        8.2761297226_dp
    xab_pm6( 73, 15)=        6.2618799210_dp
    xab_pm6( 73, 16)=        2.9759800434_dp
    xab_pm6( 73, 17)=        0.5164129734_dp
    xab_pm6( 73, 19)=        2.0267000198_dp
    xab_pm6( 73, 35)=        0.7914450169_dp
    xab_pm6( 73, 53)=        6.5515508652_dp
    xab_pm6( 73, 73)=       10.9870529175_dp
    xab_pm6( 74,  1)=        1.8322700262_dp
    xab_pm6( 74,  6)=        1.1607700586_dp
    xab_pm6( 74,  7)=        0.4783500135_dp
    xab_pm6( 74,  8)=        0.3490099907_dp
    xab_pm6( 74,  9)=        0.2138900012_dp
    xab_pm6( 74, 11)=        8.2760400772_dp
    xab_pm6( 74, 15)=        5.9538598061_dp
    xab_pm6( 74, 16)=        0.4886299968_dp
    xab_pm6( 74, 17)=        0.2779999971_dp
    xab_pm6( 74, 19)=        2.0266098976_dp
    xab_pm6( 74, 35)=        0.3723900020_dp
    xab_pm6( 74, 53)=        1.0773700476_dp
    xab_pm6( 74, 74)=        7.4713897705_dp
    xab_pm6( 75,  1)=        0.3458940089_dp
    xab_pm6( 75,  6)=        0.6906870008_dp
    xab_pm6( 75,  7)=        0.4452129900_dp
    xab_pm6( 75,  8)=        0.6359599829_dp
    xab_pm6( 75,  9)=        0.5359659791_dp
    xab_pm6( 75, 14)=        0.8494499922_dp
    xab_pm6( 75, 15)=        0.9669420123_dp
    xab_pm6( 75, 16)=        0.0688740015_dp
    xab_pm6( 75, 17)=        0.1463190019_dp
    xab_pm6( 75, 32)=        2.1515800953_dp
    xab_pm6( 75, 34)=        2.2021400928_dp
    xab_pm6( 75, 35)=        0.2875280082_dp
    xab_pm6( 75, 51)=        2.2757799625_dp
    xab_pm6( 75, 53)=        3.5592861176_dp
    xab_pm6( 75, 75)=        4.4888520241_dp
    xab_pm6( 76,  1)=        4.3938698769_dp
    xab_pm6( 76,  6)=        0.4984099865_dp
    xab_pm6( 76,  7)=        0.0808700025_dp
    xab_pm6( 76,  8)=        0.1843000054_dp
    xab_pm6( 76,  9)=        0.1400499940_dp
    xab_pm6( 76, 11)=        8.2757501602_dp
    xab_pm6( 76, 15)=        6.0583000183_dp
    xab_pm6( 76, 16)=        4.1860499382_dp
    xab_pm6( 76, 17)=        0.3279199898_dp
    xab_pm6( 76, 19)=        2.0263199806_dp
    xab_pm6( 76, 35)=        0.3824299872_dp
    xab_pm6( 76, 53)=        2.1991899014_dp
    xab_pm6( 76, 76)=        0.8304399848_dp
    xab_pm6( 77,  1)=        0.0580470003_dp
    xab_pm6( 77,  6)=        0.1150470003_dp
    xab_pm6( 77,  7)=        8.5186395645_dp
    xab_pm6( 77,  8)=        9.6254024506_dp
    xab_pm6( 77,  9)=        1.4996390343_dp
    xab_pm6( 77, 11)=        8.2758302689_dp
    xab_pm6( 77, 15)=        6.2846698761_dp
    xab_pm6( 77, 16)=        4.1357321739_dp
    xab_pm6( 77, 17)=        0.2589159906_dp
    xab_pm6( 77, 19)=        2.0264000893_dp
    xab_pm6( 77, 35)=        0.1718789935_dp
    xab_pm6( 77, 53)=        1.4971480370_dp
    xab_pm6( 77, 77)=       11.1751928329_dp
    xab_pm6( 78,  1)=        8.9240150452_dp
    xab_pm6( 78,  6)=        3.4934029579_dp
    xab_pm6( 78,  7)=        0.5407299995_dp
    xab_pm6( 78,  8)=        0.4877560139_dp
    xab_pm6( 78,  9)=        5.6170139313_dp
    xab_pm6( 78, 13)=        1.0569299459_dp
    xab_pm6( 78, 14)=        0.0999900028_dp
    xab_pm6( 78, 15)=        0.2337120026_dp
    xab_pm6( 78, 16)=        2.2242629528_dp
    xab_pm6( 78, 17)=        0.3410010040_dp
    xab_pm6( 78, 35)=        0.5203610063_dp
    xab_pm6( 78, 53)=        4.6012477875_dp
    xab_pm6( 78, 78)=        9.0102519989_dp
    xab_pm6( 79,  1)=        2.6052830219_dp
    xab_pm6( 79,  6)=       21.4856338501_dp
    xab_pm6( 79,  7)=        0.2220589966_dp
    xab_pm6( 79,  8)=        0.0771920010_dp
    xab_pm6( 79,  9)=        9.5943841934_dp
    xab_pm6( 79, 13)=        1.0571399927_dp
    xab_pm6( 79, 15)=        0.0670010000_dp
    xab_pm6( 79, 16)=       21.6191444397_dp
    xab_pm6( 79, 17)=        2.2577021122_dp
    xab_pm6( 79, 35)=        0.0042369999_dp
    xab_pm6( 79, 53)=        0.0088160001_dp
    xab_pm6( 79, 79)=        0.0130909998_dp
    xab_pm6( 80,  1)=        0.7993990183_dp
    xab_pm6( 80,  6)=        0.1471280009_dp
    xab_pm6( 80,  7)=        0.0502399988_dp
    xab_pm6( 80,  8)=        0.2407200038_dp
    xab_pm6( 80,  9)=        0.1138269976_dp
    xab_pm6( 80, 14)=        3.6807401180_dp
    xab_pm6( 80, 15)=        0.2149509937_dp
    xab_pm6( 80, 16)=        0.3473829925_dp
    xab_pm6( 80, 17)=        0.0536600016_dp
    xab_pm6( 80, 22)=        2.9572000504_dp
    xab_pm6( 80, 35)=        0.1723629981_dp
    xab_pm6( 80, 52)=        0.2127320021_dp
    xab_pm6( 80, 53)=        0.3420580029_dp
    xab_pm6( 80, 80)=        0.4232760072_dp
    xab_pm6( 81,  1)=        0.1382050067_dp
    xab_pm6( 81,  5)=       10.5043382645_dp
    xab_pm6( 81,  6)=        0.5828949809_dp
    xab_pm6( 81,  7)=        0.1588120013_dp
    xab_pm6( 81,  8)=        0.6369060278_dp
    xab_pm6( 81,  9)=        0.2261659950_dp
    xab_pm6( 81, 16)=        0.3034259975_dp
    xab_pm6( 81, 17)=        0.1620370001_dp
    xab_pm6( 81, 35)=        0.2968359888_dp
    xab_pm6( 81, 53)=        0.4300329983_dp
    xab_pm6( 81, 81)=        9.5351266861_dp
    xab_pm6( 82,  1)=        0.8400959969_dp
    xab_pm6( 82,  3)=        1.2850639820_dp
    xab_pm6( 82,  5)=        1.1381570101_dp
    xab_pm6( 82,  6)=        0.4046559930_dp
    xab_pm6( 82,  7)=        0.3357869983_dp
    xab_pm6( 82,  8)=        0.7825059891_dp
    xab_pm6( 82,  9)=        8.3685617447_dp
    xab_pm6( 82, 15)=        5.0331997871_dp
    xab_pm6( 82, 16)=        0.1751500070_dp
    xab_pm6( 82, 17)=        0.1648139954_dp
    xab_pm6( 82, 23)=        1.0000000000_dp
    xab_pm6( 82, 24)=        1.0291099548_dp
    xab_pm6( 82, 30)=        1.0000000000_dp
    xab_pm6( 82, 34)=        0.1111949980_dp
    xab_pm6( 82, 35)=        0.1482290030_dp
    xab_pm6( 82, 41)=        1.0000000000_dp
    xab_pm6( 82, 42)=        5.0000000000_dp
    xab_pm6( 82, 52)=        0.8090419769_dp
    xab_pm6( 82, 53)=        0.2674260139_dp
    xab_pm6( 82, 82)=        2.3623430729_dp
    xab_pm6( 83,  1)=        1.3974620104_dp
    xab_pm6( 83,  3)=        0.6953200102_dp
    xab_pm6( 83,  6)=        0.5761790276_dp
    xab_pm6( 83,  7)=        0.1527380049_dp
    xab_pm6( 83,  8)=        0.3330419958_dp
    xab_pm6( 83,  9)=        1.0353239775_dp
    xab_pm6( 83, 16)=        0.6209970117_dp
    xab_pm6( 83, 17)=        0.3268710077_dp
    xab_pm6( 83, 34)=        0.6512079835_dp
    xab_pm6( 83, 35)=        0.3811700046_dp
    xab_pm6( 83, 53)=        0.8623769879_dp
    xab_pm6( 83, 83)=        1.1682139635_dp

    aab_pm6(  1,  1)=        3.5409419537_dp
    aab_pm6(  2,  1)=        2.9898810387_dp
    aab_pm6(  2,  2)=        3.7835590839_dp
    aab_pm6(  3,  1)=        2.1362650394_dp
    aab_pm6(  3,  2)=        3.1124029160_dp
    aab_pm6(  3,  3)=        4.7146739960_dp
    aab_pm6(  4,  1)=        2.4754180908_dp
    aab_pm6(  4,  2)=        3.3067018986_dp
    aab_pm6(  4,  3)=        2.2367279530_dp
    aab_pm6(  4,  4)=        1.4999070168_dp
    aab_pm6(  5,  1)=        2.6152310371_dp
    aab_pm6(  5,  2)=        3.1631400585_dp
    aab_pm6(  5,  3)=        3.7593970299_dp
    aab_pm6(  5,  4)=        1.8889980316_dp
    aab_pm6(  5,  5)=        3.3186240196_dp
    aab_pm6(  6,  1)=        1.0278060000_dp
    aab_pm6(  6,  2)=        3.0427050591_dp
    aab_pm6(  6,  3)=        3.2418739796_dp
    aab_pm6(  6,  4)=        4.2128820419_dp
    aab_pm6(  6,  5)=        2.9190070629_dp
    aab_pm6(  6,  6)=        2.6137130260_dp
    aab_pm6(  7,  1)=        0.9694060087_dp
    aab_pm6(  7,  2)=        2.8143389225_dp
    aab_pm6(  7,  3)=        2.6406230927_dp
    aab_pm6(  7,  4)=        2.5808949471_dp
    aab_pm6(  7,  5)=        2.4770040512_dp
    aab_pm6(  7,  6)=        2.6861081123_dp
    aab_pm6(  7,  7)=        2.5745019913_dp
    aab_pm6(  8,  1)=        1.2609419823_dp
    aab_pm6(  8,  2)=        3.6537749767_dp
    aab_pm6(  8,  3)=        2.5844419003_dp
    aab_pm6(  8,  4)=        3.0518670082_dp
    aab_pm6(  8,  5)=        2.6953508854_dp
    aab_pm6(  8,  6)=        2.8896069527_dp
    aab_pm6(  8,  7)=        2.7842919827_dp
    aab_pm6(  8,  8)=        2.6239979267_dp
    aab_pm6(  9,  1)=        3.1367399693_dp
    aab_pm6(  9,  2)=        2.8565430641_dp
    aab_pm6(  9,  3)=        3.0439009666_dp
    aab_pm6(  9,  4)=        3.7269229889_dp
    aab_pm6(  9,  5)=        2.8238370419_dp
    aab_pm6(  9,  6)=        3.0276000500_dp
    aab_pm6(  9,  7)=        2.8566460609_dp
    aab_pm6(  9,  8)=        3.0154440403_dp
    aab_pm6(  9,  9)=        3.1757590771_dp
    aab_pm6( 10,  1)=        5.9996800423_dp
    aab_pm6( 10,  2)=        3.6777579784_dp
    aab_pm6( 10,  3)=        2.1936659813_dp
    aab_pm6( 10,  4)=        1.3165880442_dp
    aab_pm6( 10,  5)=        2.7561900616_dp
    aab_pm6( 10,  6)=        3.4411880970_dp
    aab_pm6( 10,  7)=        4.4263701439_dp
    aab_pm6( 10,  8)=        2.8895869255_dp
    aab_pm6( 10,  9)=        3.6756110191_dp
    aab_pm6( 10, 10)=        3.9745669365_dp
    aab_pm6( 11,  1)=        0.5003259778_dp
    aab_pm6( 11,  2)=        1.7030290365_dp
    aab_pm6( 11,  3)=        1.2672990561_dp
    aab_pm6( 11,  4)=        1.2554800510_dp
    aab_pm6( 11,  5)=        1.5699609518_dp
    aab_pm6( 11,  6)=        2.1960499287_dp
    aab_pm6( 11,  7)=        2.4943840504_dp
    aab_pm6( 11,  8)=        1.9814490080_dp
    aab_pm6( 11,  9)=        2.6195509434_dp
    aab_pm6( 11, 10)=        1.7742359638_dp
    aab_pm6( 11, 11)=        0.4464350045_dp
    aab_pm6( 12,  1)=        2.6515939236_dp
    aab_pm6( 12,  2)=        2.2106029987_dp
    aab_pm6( 12,  3)=        1.1843800545_dp
    aab_pm6( 12,  4)=        1.5575909615_dp
    aab_pm6( 12,  5)=        2.5274410248_dp
    aab_pm6( 12,  6)=        3.0409460068_dp
    aab_pm6( 12,  7)=        2.0791249275_dp
    aab_pm6( 12,  8)=        2.2515199184_dp
    aab_pm6( 12,  9)=        3.3622078896_dp
    aab_pm6( 12, 10)=        2.0316760540_dp
    aab_pm6( 12, 11)=        1.5067729950_dp
    aab_pm6( 12, 12)=        1.0935729742_dp
    aab_pm6( 13,  1)=        2.0259959698_dp
    aab_pm6( 13,  2)=        2.2558300495_dp
    aab_pm6( 13,  3)=        1.5815930367_dp
    aab_pm6( 13,  4)=        1.9382369518_dp
    aab_pm6( 13,  5)=        2.0595688820_dp
    aab_pm6( 13,  6)=        2.2674400806_dp
    aab_pm6( 13,  7)=        2.0097539425_dp
    aab_pm6( 13,  8)=        2.4986600876_dp
    aab_pm6( 13,  9)=        3.0842580795_dp
    aab_pm6( 13, 10)=        2.4478690624_dp
    aab_pm6( 13, 11)=        1.2028709650_dp
    aab_pm6( 13, 12)=        1.9725300074_dp
    aab_pm6( 13, 13)=        1.3877140284_dp
    aab_pm6( 14,  1)=        1.8969500065_dp
    aab_pm6( 14,  2)=        2.0404980183_dp
    aab_pm6( 14,  3)=        1.7896089554_dp
    aab_pm6( 14,  4)=        1.2631319761_dp
    aab_pm6( 14,  5)=        1.9826530218_dp
    aab_pm6( 14,  6)=        1.9844980240_dp
    aab_pm6( 14,  7)=        1.8189879656_dp
    aab_pm6( 14,  8)=        1.9235999584_dp
    aab_pm6( 14,  9)=        2.1310279369_dp
    aab_pm6( 14, 10)=        2.8677840233_dp
    aab_pm6( 14, 11)=        2.0076150894_dp
    aab_pm6( 14, 12)=        3.1397490501_dp
    aab_pm6( 14, 13)=        1.8999999762_dp
    aab_pm6( 14, 14)=        1.3289999962_dp
    aab_pm6( 15,  1)=        1.9265370369_dp
    aab_pm6( 15,  2)=        2.0931580067_dp
    aab_pm6( 15,  3)=        1.3945440054_dp
    aab_pm6( 15,  4)=        1.8000700474_dp
    aab_pm6( 15,  5)=        1.9231679440_dp
    aab_pm6( 15,  6)=        1.9946529865_dp
    aab_pm6( 15,  7)=        2.1470420361_dp
    aab_pm6( 15,  8)=        2.2207679749_dp
    aab_pm6( 15,  9)=        2.2343559265_dp
    aab_pm6( 15, 10)=        2.2190361023_dp
    aab_pm6( 15, 11)=        1.5003199577_dp
    aab_pm6( 15, 12)=        1.3837729692_dp
    aab_pm6( 15, 13)=        1.9807269573_dp
    aab_pm6( 15, 14)=        3.3134660721_dp
    aab_pm6( 15, 15)=        1.5057920218_dp
    aab_pm6( 16,  1)=        2.2159750462_dp
    aab_pm6( 16,  2)=        1.9591490030_dp
    aab_pm6( 16,  3)=        2.2942750454_dp
    aab_pm6( 16,  4)=        2.7817358971_dp
    aab_pm6( 16,  5)=        2.4036960602_dp
    aab_pm6( 16,  6)=        2.2103049755_dp
    aab_pm6( 16,  7)=        2.2899899483_dp
    aab_pm6( 16,  8)=        2.3832890987_dp
    aab_pm6( 16,  9)=        2.1871860027_dp
    aab_pm6( 16, 10)=        2.7870581150_dp
    aab_pm6( 16, 11)=        1.4008500576_dp
    aab_pm6( 16, 12)=        1.5001629591_dp
    aab_pm6( 16, 13)=        1.9767049551_dp
    aab_pm6( 16, 14)=        1.8859159946_dp
    aab_pm6( 16, 15)=        1.5953249931_dp
    aab_pm6( 16, 16)=        1.7945560217_dp
    aab_pm6( 17,  1)=        2.4028859138_dp
    aab_pm6( 17,  2)=        1.6716769934_dp
    aab_pm6( 17,  3)=        2.7830009460_dp
    aab_pm6( 17,  4)=        2.8226759434_dp
    aab_pm6( 17,  5)=        2.2593228817_dp
    aab_pm6( 17,  6)=        2.1621971130_dp
    aab_pm6( 17,  7)=        2.1721339226_dp
    aab_pm6( 17,  8)=        2.3232359886_dp
    aab_pm6( 17,  9)=        2.3132700920_dp
    aab_pm6( 17, 10)=        1.7031509876_dp
    aab_pm6( 17, 11)=        1.8164290190_dp
    aab_pm6( 17, 12)=        2.3918058872_dp
    aab_pm6( 17, 13)=        2.1259388924_dp
    aab_pm6( 17, 14)=        1.6849780083_dp
    aab_pm6( 17, 15)=        1.4683059454_dp
    aab_pm6( 17, 16)=        1.7154350281_dp
    aab_pm6( 17, 17)=        1.8232389688_dp
    aab_pm6( 18,  1)=        4.0561671257_dp
    aab_pm6( 18,  2)=        2.7165620327_dp
    aab_pm6( 18,  3)=        3.1228950024_dp
    aab_pm6( 18,  4)=        3.0440070629_dp
    aab_pm6( 18,  5)=        2.4154710770_dp
    aab_pm6( 18,  6)=        1.4713089466_dp
    aab_pm6( 18,  7)=        2.3268051147_dp
    aab_pm6( 18,  8)=        2.2406730652_dp
    aab_pm6( 18,  9)=        3.9206581116_dp
    aab_pm6( 18, 10)=        2.9637470245_dp
    aab_pm6( 18, 11)=        2.1676769257_dp
    aab_pm6( 18, 12)=        2.0926640034_dp
    aab_pm6( 18, 13)=        2.6451649666_dp
    aab_pm6( 18, 14)=        1.7803499699_dp
    aab_pm6( 18, 15)=        4.3725161552_dp
    aab_pm6( 18, 16)=        2.0493979454_dp
    aab_pm6( 18, 17)=        2.5544490814_dp
    aab_pm6( 18, 18)=        2.3064320087_dp
    aab_pm6( 19,  1)=        0.6481729746_dp
    aab_pm6( 19,  2)=        1.4185010195_dp
    aab_pm6( 19,  3)=        1.0364869833_dp
    aab_pm6( 19,  4)=        1.9318879843_dp
    aab_pm6( 19,  5)=        2.0317680836_dp
    aab_pm6( 19,  6)=        2.2417569160_dp
    aab_pm6( 19,  7)=        2.3258590698_dp
    aab_pm6( 19,  8)=        1.5085710287_dp
    aab_pm6( 19,  9)=        3.1828169823_dp
    aab_pm6( 19, 10)=        1.1380209923_dp
    aab_pm6( 19, 11)=        0.8843070269_dp
    aab_pm6( 19, 12)=        0.8848099709_dp
    aab_pm6( 19, 13)=        1.9760760069_dp
    aab_pm6( 19, 14)=        1.6759300232_dp
    aab_pm6( 19, 15)=        1.4437379837_dp
    aab_pm6( 19, 16)=        2.5121560097_dp
    aab_pm6( 19, 17)=        1.6221630573_dp
    aab_pm6( 19, 18)=        2.3028030396_dp
    aab_pm6( 19, 19)=        1.4355139732_dp
    aab_pm6( 20,  1)=        2.1418590546_dp
    aab_pm6( 20,  2)=        1.7198469639_dp
    aab_pm6( 20,  5)=        1.7000099421_dp
    aab_pm6( 20,  6)=        1.0353050232_dp
    aab_pm6( 20,  7)=        2.3866000175_dp
    aab_pm6( 20,  8)=        3.2638969421_dp
    aab_pm6( 20,  9)=        2.6450529099_dp
    aab_pm6( 20, 10)=        0.9545300007_dp
    aab_pm6( 20, 11)=        3.1071040630_dp
    aab_pm6( 20, 12)=        2.2997999191_dp
    aab_pm6( 20, 13)=        1.6125650406_dp
    aab_pm6( 20, 14)=        1.2187880278_dp
    aab_pm6( 20, 15)=        1.0241420269_dp
    aab_pm6( 20, 16)=        0.9581710100_dp
    aab_pm6( 20, 17)=        2.3833909035_dp
    aab_pm6( 20, 18)=        1.0348809958_dp
    aab_pm6( 20, 19)=        1.1191999912_dp
    aab_pm6( 20, 20)=        1.8896739483_dp
    aab_pm6( 21,  1)=        1.1794849634_dp
    aab_pm6( 21,  6)=        2.6304900646_dp
    aab_pm6( 21,  7)=        2.2700040340_dp
    aab_pm6( 21,  8)=        2.2565159798_dp
    aab_pm6( 21,  9)=        3.1079850197_dp
    aab_pm6( 21, 13)=        1.0035500526_dp
    aab_pm6( 21, 14)=        2.0168700218_dp
    aab_pm6( 21, 15)=        0.8681650162_dp
    aab_pm6( 21, 16)=        0.4229390025_dp
    aab_pm6( 21, 17)=        2.1414740086_dp
    aab_pm6( 21, 21)=        1.1328380108_dp
    aab_pm6( 22,  1)=        0.8326690197_dp
    aab_pm6( 22,  5)=        1.6287100315_dp
    aab_pm6( 22,  6)=        1.5979729891_dp
    aab_pm6( 22,  7)=        1.6786860228_dp
    aab_pm6( 22,  8)=        1.7891180515_dp
    aab_pm6( 22,  9)=        2.3070869446_dp
    aab_pm6( 22, 12)=        1.9113399982_dp
    aab_pm6( 22, 13)=        1.3694859743_dp
    aab_pm6( 22, 14)=        2.8560380936_dp
    aab_pm6( 22, 15)=        2.1519289017_dp
    aab_pm6( 22, 16)=        1.8464390039_dp
    aab_pm6( 22, 17)=        1.4610340595_dp
    aab_pm6( 22, 20)=        2.0000000000_dp
    aab_pm6( 22, 22)=        2.6485970020_dp
    aab_pm6( 23,  1)=        1.2801330090_dp
    aab_pm6( 23,  6)=        2.7898550034_dp
    aab_pm6( 23,  7)=        1.6075400114_dp
    aab_pm6( 23,  8)=        1.6239730120_dp
    aab_pm6( 23,  9)=        1.8251600266_dp
    aab_pm6( 23, 11)=        2.5510098934_dp
    aab_pm6( 23, 15)=        2.5491540432_dp
    aab_pm6( 23, 16)=        2.7041239738_dp
    aab_pm6( 23, 17)=        1.6885290146_dp
    aab_pm6( 23, 19)=        4.5213599205_dp
    aab_pm6( 23, 23)=        4.8323907852_dp
    aab_pm6( 24,  1)=        0.8826609850_dp
    aab_pm6( 24,  6)=        3.6567540169_dp
    aab_pm6( 24,  7)=        3.0291860104_dp
    aab_pm6( 24,  8)=        2.5000000000_dp
    aab_pm6( 24,  9)=        2.7165210247_dp
    aab_pm6( 24, 11)=        2.2950561047_dp
    aab_pm6( 24, 14)=        1.8607599735_dp
    aab_pm6( 24, 15)=        1.6953829527_dp
    aab_pm6( 24, 16)=        2.2609779835_dp
    aab_pm6( 24, 17)=        2.1526179314_dp
    aab_pm6( 24, 19)=        2.0000000000_dp
    aab_pm6( 24, 24)=        4.6554188728_dp
    aab_pm6( 25,  1)=        2.3099400997_dp
    aab_pm6( 25,  6)=        3.0007500648_dp
    aab_pm6( 25,  7)=        2.9214699268_dp
    aab_pm6( 25,  8)=        2.5775399208_dp
    aab_pm6( 25,  9)=        2.7919499874_dp
    aab_pm6( 25, 13)=        1.7683600187_dp
    aab_pm6( 25, 14)=        1.9379589558_dp
    aab_pm6( 25, 15)=        1.9470200539_dp
    aab_pm6( 25, 16)=        2.4825100899_dp
    aab_pm6( 25, 17)=        1.6570099592_dp
    aab_pm6( 25, 20)=        1.4914400578_dp
    aab_pm6( 25, 25)=        2.6654200554_dp
    aab_pm6( 26,  1)=        0.8544880152_dp
    aab_pm6( 26,  6)=        3.9913430214_dp
    aab_pm6( 26,  7)=        2.5004858971_dp
    aab_pm6( 26,  8)=        1.7263129950_dp
    aab_pm6( 26,  9)=        4.2947068214_dp
    aab_pm6( 26, 15)=        2.5675339699_dp
    aab_pm6( 26, 16)=        0.9889910221_dp
    aab_pm6( 26, 17)=        1.2297929525_dp
    aab_pm6( 26, 19)=        2.0000000000_dp
    aab_pm6( 26, 26)=        2.7207849026_dp
    aab_pm6( 27,  1)=        2.9665179253_dp
    aab_pm6( 27,  6)=        3.7162330151_dp
    aab_pm6( 27,  7)=        3.6186380386_dp
    aab_pm6( 27,  8)=        3.7269110680_dp
    aab_pm6( 27,  9)=        3.9563469887_dp
    aab_pm6( 27, 14)=        2.4698050022_dp
    aab_pm6( 27, 15)=        1.1525050402_dp
    aab_pm6( 27, 16)=        2.4292550087_dp
    aab_pm6( 27, 17)=        3.2174971104_dp
    aab_pm6( 27, 27)=        3.2881660461_dp
    aab_pm6( 28,  1)=        2.6352798939_dp
    aab_pm6( 28,  6)=        4.2855129242_dp
    aab_pm6( 28,  7)=        3.8452150822_dp
    aab_pm6( 28,  8)=        2.9372320175_dp
    aab_pm6( 28,  9)=        3.4402410984_dp
    aab_pm6( 28, 14)=        2.0688810349_dp
    aab_pm6( 28, 15)=        3.2602829933_dp
    aab_pm6( 28, 16)=        2.0027520657_dp
    aab_pm6( 28, 17)=        2.2005119324_dp
    aab_pm6( 28, 28)=        1.0979599953_dp
    aab_pm6( 29,  1)=        2.3353590965_dp
    aab_pm6( 29,  6)=        4.6387729645_dp
    aab_pm6( 29,  7)=        4.2143368721_dp
    aab_pm6( 29,  8)=        3.9599509239_dp
    aab_pm6( 29,  9)=        4.4788317680_dp
    aab_pm6( 29, 15)=        0.2106399983_dp
    aab_pm6( 29, 16)=        0.2731119990_dp
    aab_pm6( 29, 17)=        2.7765309811_dp
    aab_pm6( 29, 29)=        3.6168460846_dp
    aab_pm6( 30,  1)=        1.9878909588_dp
    aab_pm6( 30,  6)=        1.8023270369_dp
    aab_pm6( 30,  7)=        1.8445789814_dp
    aab_pm6( 30,  8)=        2.3350539207_dp
    aab_pm6( 30,  9)=        2.4100210667_dp
    aab_pm6( 30, 14)=        1.8320579529_dp
    aab_pm6( 30, 15)=        1.2204799652_dp
    aab_pm6( 30, 16)=        1.4550000429_dp
    aab_pm6( 30, 17)=        1.6251759529_dp
    aab_pm6( 30, 20)=        1.1191799641_dp
    aab_pm6( 30, 30)=        0.9290000200_dp
    aab_pm6( 31,  1)=        1.8473500013_dp
    aab_pm6( 31,  6)=        2.3254098892_dp
    aab_pm6( 31,  7)=        2.1218199730_dp
    aab_pm6( 31,  8)=        2.3483469486_dp
    aab_pm6( 31,  9)=        2.6798689365_dp
    aab_pm6( 31, 14)=        1.9137799740_dp
    aab_pm6( 31, 15)=        2.9796500206_dp
    aab_pm6( 31, 16)=        2.2321081161_dp
    aab_pm6( 31, 17)=        2.0247099400_dp
    aab_pm6( 31, 31)=        1.3346430063_dp
    aab_pm6( 32,  1)=        2.2067930698_dp
    aab_pm6( 32,  6)=        2.2574689388_dp
    aab_pm6( 32,  7)=        1.9882260561_dp
    aab_pm6( 32,  8)=        2.1394131184_dp
    aab_pm6( 32,  9)=        2.3847770691_dp
    aab_pm6( 32, 14)=        0.2997210026_dp
    aab_pm6( 32, 15)=        2.4692909718_dp
    aab_pm6( 32, 16)=        2.0245881081_dp
    aab_pm6( 32, 17)=        1.7712279558_dp
    aab_pm6( 32, 25)=        2.3828339577_dp
    aab_pm6( 32, 27)=        2.8526101112_dp
    aab_pm6( 32, 32)=        2.0190000534_dp
    aab_pm6( 33,  1)=        1.9935270548_dp
    aab_pm6( 33,  6)=        1.8550690413_dp
    aab_pm6( 33,  7)=        1.4965430498_dp
    aab_pm6( 33,  8)=        2.0039501190_dp
    aab_pm6( 33,  9)=        2.0125830173_dp
    aab_pm6( 33, 13)=        1.1527860165_dp
    aab_pm6( 33, 14)=        1.9155999422_dp
    aab_pm6( 33, 16)=        1.9543679953_dp
    aab_pm6( 33, 17)=        1.6910699606_dp
    aab_pm6( 33, 22)=        1.9329110384_dp
    aab_pm6( 33, 27)=        3.3681399822_dp
    aab_pm6( 33, 30)=        1.4591300488_dp
    aab_pm6( 33, 31)=        1.7309770584_dp
    aab_pm6( 33, 33)=        1.5882639885_dp
    aab_pm6( 34,  1)=        2.0350680351_dp
    aab_pm6( 34,  6)=        2.3871181011_dp
    aab_pm6( 34,  7)=        1.9377640486_dp
    aab_pm6( 34,  8)=        2.4842629433_dp
    aab_pm6( 34,  9)=        2.3021800518_dp
    aab_pm6( 34, 14)=        1.5298169851_dp
    aab_pm6( 34, 15)=        1.0481829643_dp
    aab_pm6( 34, 16)=        1.4796060324_dp
    aab_pm6( 34, 17)=        2.1288609505_dp
    aab_pm6( 34, 25)=        2.6480379105_dp
    aab_pm6( 34, 27)=        2.5234498978_dp
    aab_pm6( 34, 30)=        1.1862419844_dp
    aab_pm6( 34, 32)=        2.6690568924_dp
    aab_pm6( 34, 33)=        1.6652799845_dp
    aab_pm6( 34, 34)=        1.7958940268_dp
    aab_pm6( 35,  1)=        2.1928029060_dp
    aab_pm6( 35,  2)=        2.1282749176_dp
    aab_pm6( 35,  3)=        2.0744409561_dp
    aab_pm6( 35,  4)=        2.3671460152_dp
    aab_pm6( 35,  5)=        2.3078899384_dp
    aab_pm6( 35,  6)=        2.0150859356_dp
    aab_pm6( 35,  7)=        4.2249011993_dp
    aab_pm6( 35,  8)=        2.2830460072_dp
    aab_pm6( 35,  9)=        2.0317649841_dp
    aab_pm6( 35, 10)=        2.4641718864_dp
    aab_pm6( 35, 11)=        1.6222180128_dp
    aab_pm6( 35, 12)=        2.1956970692_dp
    aab_pm6( 35, 13)=        1.8941409588_dp
    aab_pm6( 35, 14)=        1.5708249807_dp
    aab_pm6( 35, 15)=        1.4021389484_dp
    aab_pm6( 35, 16)=        1.5098739862_dp
    aab_pm6( 35, 17)=        1.7103309631_dp
    aab_pm6( 35, 18)=        2.4508008957_dp
    aab_pm6( 35, 19)=        1.6160930395_dp
    aab_pm6( 35, 20)=        2.0784049034_dp
    aab_pm6( 35, 21)=        1.7934859991_dp
    aab_pm6( 35, 22)=        1.6748470068_dp
    aab_pm6( 35, 23)=        1.9029040337_dp
    aab_pm6( 35, 24)=        1.5660279989_dp
    aab_pm6( 35, 25)=        2.2838199139_dp
    aab_pm6( 35, 26)=        3.6417820454_dp
    aab_pm6( 35, 27)=        2.6326880455_dp
    aab_pm6( 35, 28)=        2.7721359730_dp
    aab_pm6( 35, 29)=        5.8264069557_dp
    aab_pm6( 35, 30)=        1.4161200523_dp
    aab_pm6( 35, 31)=        1.8191050291_dp
    aab_pm6( 35, 32)=        1.6023659706_dp
    aab_pm6( 35, 33)=        1.5201699734_dp
    aab_pm6( 35, 34)=        1.4837130308_dp
    aab_pm6( 35, 35)=        1.7581460476_dp
    aab_pm6( 36,  1)=        3.7704529762_dp
    aab_pm6( 36,  2)=        1.9969429970_dp
    aab_pm6( 36,  3)=        3.3145620823_dp
    aab_pm6( 36,  4)=        3.2530479431_dp
    aab_pm6( 36,  5)=        2.3631689548_dp
    aab_pm6( 36,  6)=        2.0767381191_dp
    aab_pm6( 36,  7)=        1.6440520287_dp
    aab_pm6( 36,  8)=        0.2922999859_dp
    aab_pm6( 36,  9)=        3.4523210526_dp
    aab_pm6( 36, 10)=        2.8136789799_dp
    aab_pm6( 36, 11)=        2.4805979729_dp
    aab_pm6( 36, 12)=        1.3914870024_dp
    aab_pm6( 36, 13)=        2.4671308994_dp
    aab_pm6( 36, 14)=        1.7640999556_dp
    aab_pm6( 36, 17)=        1.8849740028_dp
    aab_pm6( 36, 18)=        1.9951250553_dp
    aab_pm6( 36, 19)=        2.1824870110_dp
    aab_pm6( 36, 20)=        1.3051970005_dp
    aab_pm6( 36, 35)=        1.5290060043_dp
    aab_pm6( 36, 36)=        1.1353189945_dp
    aab_pm6( 37,  1)=        2.4435560703_dp
    aab_pm6( 37,  2)=        1.2707409859_dp
    aab_pm6( 37,  5)=        5.5322389603_dp
    aab_pm6( 37,  6)=        2.7658300400_dp
    aab_pm6( 37,  7)=        0.7610470057_dp
    aab_pm6( 37,  8)=        1.3349080086_dp
    aab_pm6( 37,  9)=        3.6381220818_dp
    aab_pm6( 37, 10)=        2.2675909996_dp
    aab_pm6( 37, 13)=        0.7987740040_dp
    aab_pm6( 37, 16)=        1.3031840324_dp
    aab_pm6( 37, 17)=        2.2744109631_dp
    aab_pm6( 37, 18)=        2.5109770298_dp
    aab_pm6( 37, 35)=        1.7977659702_dp
    aab_pm6( 37, 36)=        2.2687530518_dp
    aab_pm6( 37, 37)=        1.1808179617_dp
    aab_pm6( 38,  1)=        2.1059141159_dp
    aab_pm6( 38,  6)=        1.9866880178_dp
    aab_pm6( 38,  7)=        2.1836290359_dp
    aab_pm6( 38,  8)=        2.1383988857_dp
    aab_pm6( 38,  9)=        3.0506660938_dp
    aab_pm6( 38, 14)=        2.9697799683_dp
    aab_pm6( 38, 15)=        2.7891499996_dp
    aab_pm6( 38, 16)=        1.5981060266_dp
    aab_pm6( 38, 17)=        1.8541899920_dp
    aab_pm6( 38, 22)=        2.8800299168_dp
    aab_pm6( 38, 35)=        1.5243159533_dp
    aab_pm6( 38, 38)=        1.0000400543_dp
    aab_pm6( 39,  1)=        1.1890530586_dp
    aab_pm6( 39,  6)=        1.3360940218_dp
    aab_pm6( 39,  7)=        1.7787959576_dp
    aab_pm6( 39,  8)=        1.8510299921_dp
    aab_pm6( 39,  9)=        2.6480460167_dp
    aab_pm6( 39, 13)=        1.0034999847_dp
    aab_pm6( 39, 14)=        2.0168199539_dp
    aab_pm6( 39, 15)=        0.9544500113_dp
    aab_pm6( 39, 16)=        0.9716879725_dp
    aab_pm6( 39, 17)=        1.6301519871_dp
    aab_pm6( 39, 35)=        1.4012080431_dp
    aab_pm6( 39, 39)=        1.0126810074_dp
    aab_pm6( 40,  1)=        1.3797030449_dp
    aab_pm6( 40,  6)=        2.0294270515_dp
    aab_pm6( 40,  7)=        1.7070829868_dp
    aab_pm6( 40,  8)=        1.7095700502_dp
    aab_pm6( 40,  9)=        1.9009250402_dp
    aab_pm6( 40, 13)=        1.2706199884_dp
    aab_pm6( 40, 14)=        1.7508330345_dp
    aab_pm6( 40, 15)=        1.0918580294_dp
    aab_pm6( 40, 16)=        2.1297609806_dp
    aab_pm6( 40, 17)=        1.3288350105_dp
    aab_pm6( 40, 35)=        1.4468679428_dp
    aab_pm6( 40, 40)=        3.8659679890_dp
    aab_pm6( 41,  1)=        2.5059120655_dp
    aab_pm6( 41,  6)=        2.6210119724_dp
    aab_pm6( 41,  7)=        2.0238630772_dp
    aab_pm6( 41,  8)=        2.0494890213_dp
    aab_pm6( 41,  9)=        3.0031569004_dp
    aab_pm6( 41, 11)=        2.5510098934_dp
    aab_pm6( 41, 15)=        2.2216079235_dp
    aab_pm6( 41, 16)=        2.2494819164_dp
    aab_pm6( 41, 17)=        2.2152750492_dp
    aab_pm6( 41, 19)=        4.5213599205_dp
    aab_pm6( 41, 35)=        2.0066781044_dp
    aab_pm6( 41, 41)=        1.7279410362_dp
    aab_pm6( 42,  1)=        2.0357480049_dp
    aab_pm6( 42,  6)=        2.1986720562_dp
    aab_pm6( 42,  7)=        1.8694750071_dp
    aab_pm6( 42,  8)=        1.7554240227_dp
    aab_pm6( 42,  9)=        2.2025930882_dp
    aab_pm6( 42, 11)=        2.4407699108_dp
    aab_pm6( 42, 15)=        1.8504409790_dp
    aab_pm6( 42, 16)=        1.9396580458_dp
    aab_pm6( 42, 17)=        1.7833620310_dp
    aab_pm6( 42, 19)=        3.9394199848_dp
    aab_pm6( 42, 24)=        2.6746160984_dp
    aab_pm6( 42, 35)=        1.2833340168_dp
    aab_pm6( 42, 42)=        2.0342540741_dp
    aab_pm6( 43,  1)=        2.8303449154_dp
    aab_pm6( 43,  6)=        3.1983261108_dp
    aab_pm6( 43,  7)=        2.3154170513_dp
    aab_pm6( 43,  8)=        2.4051899910_dp
    aab_pm6( 43,  9)=        3.6048150063_dp
    aab_pm6( 43, 16)=        2.4634010792_dp
    aab_pm6( 43, 17)=        2.5720429420_dp
    aab_pm6( 43, 32)=        2.8528199196_dp
    aab_pm6( 43, 34)=        2.5236599445_dp
    aab_pm6( 43, 35)=        2.8282639980_dp
    aab_pm6( 44,  1)=        2.8928990364_dp
    aab_pm6( 44,  6)=        2.7848329544_dp
    aab_pm6( 44,  7)=        3.0555040836_dp
    aab_pm6( 44,  8)=        3.1349399090_dp
    aab_pm6( 44,  9)=        3.8787109852_dp
    aab_pm6( 44, 14)=        2.7759099007_dp
    aab_pm6( 44, 15)=        0.2989160120_dp
    aab_pm6( 44, 16)=        2.5080759525_dp
    aab_pm6( 44, 17)=        1.7598830462_dp
    aab_pm6( 44, 32)=        2.8523199558_dp
    aab_pm6( 44, 34)=        2.5231599808_dp
    aab_pm6( 44, 35)=        2.5847349167_dp
    aab_pm6( 44, 44)=        0.5720559955_dp
    aab_pm6( 45,  1)=        3.1041650772_dp
    aab_pm6( 45,  6)=        3.4159910679_dp
    aab_pm6( 45,  7)=        3.5854620934_dp
    aab_pm6( 45,  8)=        3.9278299809_dp
    aab_pm6( 45,  9)=        4.0516538620_dp
    aab_pm6( 45, 14)=        2.7764899731_dp
    aab_pm6( 45, 15)=        2.3346068859_dp
    aab_pm6( 45, 16)=        3.1540060043_dp
    aab_pm6( 45, 17)=        3.3001298904_dp
    aab_pm6( 45, 32)=        2.8529000282_dp
    aab_pm6( 45, 34)=        2.5237400532_dp
    aab_pm6( 45, 35)=        2.9280819893_dp
    aab_pm6( 45, 45)=        2.4973280430_dp
    aab_pm6( 46,  1)=        2.1837608814_dp
    aab_pm6( 46,  6)=        4.7771921158_dp
    aab_pm6( 46,  7)=        2.3280460835_dp
    aab_pm6( 46,  8)=        2.1548669338_dp
    aab_pm6( 46,  9)=        4.2373118401_dp
    aab_pm6( 46, 13)=        1.5727200508_dp
    aab_pm6( 46, 14)=        2.9481999874_dp
    aab_pm6( 46, 15)=        0.8036299944_dp
    aab_pm6( 46, 16)=        2.1778008938_dp
    aab_pm6( 46, 17)=        3.8712430000_dp
    aab_pm6( 46, 35)=        5.9948787689_dp
    aab_pm6( 46, 46)=        1.0643750429_dp
    aab_pm6( 47,  1)=        2.8959360123_dp
    aab_pm6( 47,  6)=        4.4043359756_dp
    aab_pm6( 47,  7)=        4.6598711014_dp
    aab_pm6( 47,  8)=        1.8938740492_dp
    aab_pm6( 47,  9)=        4.6284232140_dp
    aab_pm6( 47, 13)=        1.9287999868_dp
    aab_pm6( 47, 15)=        6.0000061989_dp
    aab_pm6( 47, 16)=        3.6531209946_dp
    aab_pm6( 47, 17)=        4.4411759377_dp
    aab_pm6( 47, 35)=        3.6774909496_dp
    aab_pm6( 47, 47)=        2.1276450157_dp
    aab_pm6( 48,  1)=        2.6287479401_dp
    aab_pm6( 48,  6)=        1.4256780148_dp
    aab_pm6( 48,  7)=        0.9704229832_dp
    aab_pm6( 48,  8)=        1.6966730356_dp
    aab_pm6( 48,  9)=        2.3121349812_dp
    aab_pm6( 48, 14)=        1.3712249994_dp
    aab_pm6( 48, 16)=        1.1822019815_dp
    aab_pm6( 48, 17)=        0.9435470104_dp
    aab_pm6( 48, 35)=        1.0014510155_dp
    aab_pm6( 48, 48)=        1.5640439987_dp
    aab_pm6( 49,  1)=        3.0641438961_dp
    aab_pm6( 49,  6)=        2.1892719269_dp
    aab_pm6( 49,  7)=        2.4698679447_dp
    aab_pm6( 49,  8)=        2.6620950699_dp
    aab_pm6( 49,  9)=        2.9487969875_dp
    aab_pm6( 49, 16)=        2.5421309471_dp
    aab_pm6( 49, 17)=        2.2334051132_dp
    aab_pm6( 49, 31)=        1.6288700104_dp
    aab_pm6( 49, 33)=        2.2995519638_dp
    aab_pm6( 49, 34)=        1.9065719843_dp
    aab_pm6( 49, 35)=        2.2579569817_dp
    aab_pm6( 49, 49)=        2.0732409954_dp
    aab_pm6( 50,  1)=        2.6489100456_dp
    aab_pm6( 50,  6)=        2.4405379295_dp
    aab_pm6( 50,  7)=        2.0855889320_dp
    aab_pm6( 50,  8)=        2.7272601128_dp
    aab_pm6( 50,  9)=        3.7242860794_dp
    aab_pm6( 50, 16)=        2.1315419674_dp
    aab_pm6( 50, 17)=        1.7715220451_dp
    aab_pm6( 50, 32)=        2.5246329308_dp
    aab_pm6( 50, 34)=        2.1273770332_dp
    aab_pm6( 50, 35)=        1.5350890160_dp
    aab_pm6( 50, 50)=        0.9210000038_dp
    aab_pm6( 51,  1)=        1.5712720156_dp
    aab_pm6( 51,  6)=        1.6962059736_dp
    aab_pm6( 51,  7)=        0.6761149764_dp
    aab_pm6( 51,  8)=        1.8463840485_dp
    aab_pm6( 51,  9)=        2.1829218864_dp
    aab_pm6( 51, 13)=        1.4226410389_dp
    aab_pm6( 51, 14)=        2.6865899563_dp
    aab_pm6( 51, 16)=        1.4188369513_dp
    aab_pm6( 51, 17)=        1.1172870398_dp
    aab_pm6( 51, 25)=        2.4003200531_dp
    aab_pm6( 51, 27)=        2.2046298981_dp
    aab_pm6( 51, 35)=        1.0639159679_dp
    aab_pm6( 51, 43)=        2.2048499584_dp
    aab_pm6( 51, 44)=        2.2043499947_dp
    aab_pm6( 51, 45)=        2.2049300671_dp
    aab_pm6( 51, 49)=        2.1419329643_dp
    aab_pm6( 51, 51)=        1.3485349417_dp
    aab_pm6( 52,  1)=        2.0391299725_dp
    aab_pm6( 52,  6)=        1.9928159714_dp
    aab_pm6( 52,  7)=        1.7222690582_dp
    aab_pm6( 52,  8)=        1.8530639410_dp
    aab_pm6( 52,  9)=        1.9985760450_dp
    aab_pm6( 52, 13)=        1.3875410557_dp
    aab_pm6( 52, 15)=        1.4537179470_dp
    aab_pm6( 52, 16)=        1.8301700354_dp
    aab_pm6( 52, 17)=        1.3002599478_dp
    aab_pm6( 52, 30)=        1.2189290524_dp
    aab_pm6( 52, 32)=        2.3423719406_dp
    aab_pm6( 52, 33)=        1.1892529726_dp
    aab_pm6( 52, 34)=        1.5660079718_dp
    aab_pm6( 52, 35)=        1.2509399652_dp
    aab_pm6( 52, 48)=        1.3072619438_dp
    aab_pm6( 52, 49)=        1.5409879684_dp
    aab_pm6( 52, 50)=        1.7639410496_dp
    aab_pm6( 52, 52)=        1.1649780273_dp
    aab_pm6( 53,  1)=        2.1399130821_dp
    aab_pm6( 53,  2)=        2.1729838848_dp
    aab_pm6( 53,  3)=        2.1212511063_dp
    aab_pm6( 53,  4)=        2.2880229950_dp
    aab_pm6( 53,  5)=        2.6676049232_dp
    aab_pm6( 53,  6)=        2.0687100887_dp
    aab_pm6( 53,  7)=        1.6775180101_dp
    aab_pm6( 53,  8)=        2.2889189720_dp
    aab_pm6( 53,  9)=        2.2035799026_dp
    aab_pm6( 53, 10)=        2.4144148827_dp
    aab_pm6( 53, 11)=        1.4030900002_dp
    aab_pm6( 53, 12)=        2.0451369286_dp
    aab_pm6( 53, 13)=        1.8160680532_dp
    aab_pm6( 53, 14)=        1.5595790148_dp
    aab_pm6( 53, 15)=        2.1315929890_dp
    aab_pm6( 53, 16)=        1.8551100492_dp
    aab_pm6( 53, 17)=        1.5741610527_dp
    aab_pm6( 53, 18)=        1.5765869617_dp
    aab_pm6( 53, 19)=        1.5397139788_dp
    aab_pm6( 53, 20)=        2.1964900494_dp
    aab_pm6( 53, 21)=        1.8148839474_dp
    aab_pm6( 53, 22)=        1.9334690571_dp
    aab_pm6( 53, 23)=        2.6835200787_dp
    aab_pm6( 53, 24)=        2.6342239380_dp
    aab_pm6( 53, 25)=        2.2665998936_dp
    aab_pm6( 53, 26)=        1.9128290415_dp
    aab_pm6( 53, 27)=        3.2352039814_dp
    aab_pm6( 53, 28)=        1.0853430033_dp
    aab_pm6( 53, 29)=        0.8343049884_dp
    aab_pm6( 53, 30)=        1.3947620392_dp
    aab_pm6( 53, 31)=        1.6717289686_dp
    aab_pm6( 53, 32)=        1.8174250126_dp
    aab_pm6( 53, 33)=        1.2452620268_dp
    aab_pm6( 53, 35)=        1.5793759823_dp
    aab_pm6( 53, 36)=        1.2385740280_dp
    aab_pm6( 53, 37)=        1.4326750040_dp
    aab_pm6( 53, 38)=        1.2620420456_dp
    aab_pm6( 53, 39)=        1.2791099548_dp
    aab_pm6( 53, 40)=        1.9951820374_dp
    aab_pm6( 53, 41)=        1.9672509432_dp
    aab_pm6( 53, 42)=        0.9484609962_dp
    aab_pm6( 53, 43)=        1.2923120260_dp
    aab_pm6( 53, 44)=        3.9532029629_dp
    aab_pm6( 53, 45)=        3.7081699371_dp
    aab_pm6( 53, 46)=        5.1445441246_dp
    aab_pm6( 53, 47)=        2.5931611061_dp
    aab_pm6( 53, 48)=        0.9962379932_dp
    aab_pm6( 53, 49)=        2.3517580032_dp
    aab_pm6( 53, 50)=        1.8556330204_dp
    aab_pm6( 53, 51)=        1.1553150415_dp
    aab_pm6( 53, 52)=        1.4939509630_dp
    aab_pm6( 53, 53)=        1.5199249983_dp
    aab_pm6( 54,  1)=        1.3568609953_dp
    aab_pm6( 54,  2)=        2.4978320599_dp
    aab_pm6( 54,  3)=        2.4668951035_dp
    aab_pm6( 54,  4)=        6.0000028610_dp
    aab_pm6( 54,  5)=        5.0519571304_dp
    aab_pm6( 54,  6)=        1.7044399977_dp
    aab_pm6( 54,  7)=        1.9329520464_dp
    aab_pm6( 54,  8)=        0.8392329812_dp
    aab_pm6( 54,  9)=        1.1288119555_dp
    aab_pm6( 54, 10)=        1.3302019835_dp
    aab_pm6( 54, 11)=        2.1030030251_dp
    aab_pm6( 54, 12)=        2.6984140873_dp
    aab_pm6( 54, 13)=        2.4120390415_dp
    aab_pm6( 54, 14)=        3.0870599747_dp
    aab_pm6( 54, 17)=        1.5463960171_dp
    aab_pm6( 54, 18)=        0.5915200114_dp
    aab_pm6( 54, 19)=        1.1712499857_dp
    aab_pm6( 54, 20)=        1.5106530190_dp
    aab_pm6( 54, 35)=        1.4396179914_dp
    aab_pm6( 54, 36)=        0.5515609980_dp
    aab_pm6( 54, 37)=        1.0878230333_dp
    aab_pm6( 54, 53)=        0.7991549969_dp
    aab_pm6( 54, 54)=        1.2447619438_dp
    aab_pm6( 55,  1)=        0.2648819983_dp
    aab_pm6( 55,  5)=        1.4871100187_dp
    aab_pm6( 55,  6)=        2.1471040249_dp
    aab_pm6( 55,  7)=        2.4465320110_dp
    aab_pm6( 55,  8)=        2.0851390362_dp
    aab_pm6( 55,  9)=        2.8341000080_dp
    aab_pm6( 55, 15)=        2.9249529839_dp
    aab_pm6( 55, 16)=        0.2894119918_dp
    aab_pm6( 55, 17)=        1.6736630201_dp
    aab_pm6( 55, 35)=        1.1671890020_dp
    aab_pm6( 55, 53)=        0.9195619822_dp
    aab_pm6( 55, 55)=        1.1708430052_dp
    aab_pm6( 56,  1)=        6.0001349449_dp
    aab_pm6( 56,  6)=        0.7706260085_dp
    aab_pm6( 56,  7)=        1.1482330561_dp
    aab_pm6( 56,  8)=        1.2830179930_dp
    aab_pm6( 56,  9)=        3.0006179810_dp
    aab_pm6( 56, 13)=        2.1059238911_dp
    aab_pm6( 56, 14)=        1.2404199839_dp
    aab_pm6( 56, 16)=        0.7051879764_dp
    aab_pm6( 56, 17)=        1.0710439682_dp
    aab_pm6( 56, 22)=        2.1760399342_dp
    aab_pm6( 56, 35)=        1.1903460026_dp
    aab_pm6( 56, 53)=        0.9825279713_dp
    aab_pm6( 56, 56)=        0.3392690122_dp
    aab_pm6( 57,  1)=        0.8336669803_dp
    aab_pm6( 57,  6)=        0.6048690081_dp
    aab_pm6( 57,  7)=        0.7588809729_dp
    aab_pm6( 57,  8)=        1.3183330297_dp
    aab_pm6( 57,  9)=        2.3793349266_dp
    aab_pm6( 57, 13)=        1.0035099983_dp
    aab_pm6( 57, 14)=        2.0168199539_dp
    aab_pm6( 57, 15)=        0.9544500113_dp
    aab_pm6( 57, 16)=        1.8341289759_dp
    aab_pm6( 57, 17)=        0.9937530160_dp
    aab_pm6( 57, 35)=        0.7581840158_dp
    aab_pm6( 57, 53)=        0.5926659703_dp
    aab_pm6( 57, 57)=        4.2480669022_dp
    aab_pm6( 64,  1)=        0.3908700049_dp
    aab_pm6( 64,  6)=        0.4468699992_dp
    aab_pm6( 64,  7)=        1.1594099998_dp
    aab_pm6( 64,  8)=        0.8620399833_dp
    aab_pm6( 64,  9)=        1.4979799986_dp
    aab_pm6( 64, 13)=        1.0035099983_dp
    aab_pm6( 64, 14)=        2.0168199539_dp
    aab_pm6( 64, 15)=        0.9544500113_dp
    aab_pm6( 64, 16)=        2.0039300919_dp
    aab_pm6( 64, 17)=        0.8068100214_dp
    aab_pm6( 64, 35)=        0.7158100009_dp
    aab_pm6( 64, 53)=        0.5853599906_dp
    aab_pm6( 64, 64)=        3.3481800556_dp
    aab_pm6( 71,  1)=        1.4157899618_dp
    aab_pm6( 71,  6)=        2.3128130436_dp
    aab_pm6( 71,  7)=        2.1413021088_dp
    aab_pm6( 71,  8)=        2.1924860477_dp
    aab_pm6( 71, 15)=        5.6188201904_dp
    aab_pm6( 71, 17)=        2.7536358833_dp
    aab_pm6( 71, 35)=        2.3226180077_dp
    aab_pm6( 71, 53)=        2.2483479977_dp
    aab_pm6( 72,  1)=        1.4237879515_dp
    aab_pm6( 72,  5)=        1.6334999800_dp
    aab_pm6( 72,  6)=        1.0021940470_dp
    aab_pm6( 72,  7)=        1.3324099779_dp
    aab_pm6( 72,  8)=        1.6332889795_dp
    aab_pm6( 72,  9)=        2.2908029556_dp
    aab_pm6( 72, 12)=        1.9113500118_dp
    aab_pm6( 72, 13)=        0.9491500258_dp
    aab_pm6( 72, 14)=        2.1893000603_dp
    aab_pm6( 72, 15)=        1.2312200069_dp
    aab_pm6( 72, 16)=        2.3271100521_dp
    aab_pm6( 72, 17)=        1.2971169949_dp
    aab_pm6( 72, 20)=        2.0545001030_dp
    aab_pm6( 72, 33)=        1.7994999886_dp
    aab_pm6( 72, 35)=        1.0907590389_dp
    aab_pm6( 72, 53)=        1.0140960217_dp
    aab_pm6( 72, 56)=        2.2648301125_dp
    aab_pm6( 72, 72)=        0.5441439748_dp
    aab_pm6( 73,  1)=        2.2880139351_dp
    aab_pm6( 73,  6)=        1.8389489651_dp
    aab_pm6( 73,  7)=        2.0536789894_dp
    aab_pm6( 73,  8)=        2.4126288891_dp
    aab_pm6( 73,  9)=        3.1073899269_dp
    aab_pm6( 73, 11)=        2.5511200428_dp
    aab_pm6( 73, 15)=        2.5137999058_dp
    aab_pm6( 73, 16)=        2.2467229366_dp
    aab_pm6( 73, 17)=        1.6088049412_dp
    aab_pm6( 73, 19)=        4.5214700699_dp
    aab_pm6( 73, 35)=        1.6403759718_dp
    aab_pm6( 73, 53)=        2.4010529518_dp
    aab_pm6( 73, 73)=        2.0828630924_dp
    aab_pm6( 74,  1)=        2.1308801174_dp
    aab_pm6( 74,  6)=        2.0974800587_dp
    aab_pm6( 74,  7)=        1.5960400105_dp
    aab_pm6( 74,  8)=        1.3590199947_dp
    aab_pm6( 74,  9)=        1.4460500479_dp
    aab_pm6( 74, 11)=        2.5510299206_dp
    aab_pm6( 74, 15)=        2.3380599022_dp
    aab_pm6( 74, 16)=        1.5425699949_dp
    aab_pm6( 74, 17)=        1.3106900454_dp
    aab_pm6( 74, 19)=        4.5213799477_dp
    aab_pm6( 74, 35)=        1.2932599783_dp
    aab_pm6( 74, 53)=        1.5735700130_dp
    aab_pm6( 74, 74)=        2.9408700466_dp
    aab_pm6( 75,  1)=        1.6345000267_dp
    aab_pm6( 75,  6)=        2.3062849045_dp
    aab_pm6( 75,  7)=        1.9183319807_dp
    aab_pm6( 75,  8)=        1.9677469730_dp
    aab_pm6( 75,  9)=        2.1542189121_dp
    aab_pm6( 75, 14)=        2.7759299278_dp
    aab_pm6( 75, 15)=        1.8041679859_dp
    aab_pm6( 75, 16)=        1.0839190483_dp
    aab_pm6( 75, 17)=        1.4338749647_dp
    aab_pm6( 75, 32)=        2.8523399830_dp
    aab_pm6( 75, 34)=        2.5231699944_dp
    aab_pm6( 75, 35)=        1.6030600071_dp
    aab_pm6( 75, 51)=        2.2043600082_dp
    aab_pm6( 75, 53)=        2.6101191044_dp
    aab_pm6( 75, 75)=        6.0002579689_dp
    aab_pm6( 76,  1)=        3.4041800499_dp
    aab_pm6( 76,  6)=        2.3364999294_dp
    aab_pm6( 76,  7)=        1.1430900097_dp
    aab_pm6( 76,  8)=        1.3503600359_dp
    aab_pm6( 76,  9)=        1.5076199770_dp
    aab_pm6( 76, 11)=        2.5507400036_dp
    aab_pm6( 76, 15)=        2.8360900879_dp
    aab_pm6( 76, 16)=        2.8094999790_dp
    aab_pm6( 76, 17)=        1.8330700397_dp
    aab_pm6( 76, 19)=        4.5210900307_dp
    aab_pm6( 76, 35)=        1.7668800354_dp
    aab_pm6( 76, 53)=        2.2037599087_dp
    aab_pm6( 76, 76)=        2.0216300488_dp
    aab_pm6( 77,  1)=        1.0339000225_dp
    aab_pm6( 77,  6)=        1.6902949810_dp
    aab_pm6( 77,  7)=        3.9345080853_dp
    aab_pm6( 77,  8)=        3.7482719421_dp
    aab_pm6( 77,  9)=        2.9827990532_dp
    aab_pm6( 77, 11)=        2.5508201122_dp
    aab_pm6( 77, 15)=        2.7140600681_dp
    aab_pm6( 77, 16)=        3.2048339844_dp
    aab_pm6( 77, 17)=        2.0097699165_dp
    aab_pm6( 77, 19)=        4.5211701393_dp
    aab_pm6( 77, 35)=        2.0381419659_dp
    aab_pm6( 77, 53)=        3.4109139442_dp
    aab_pm6( 77, 77)=        5.7716631889_dp
    aab_pm6( 78,  1)=        4.0011978149_dp
    aab_pm6( 78,  6)=        3.3067219257_dp
    aab_pm6( 78,  7)=        2.3079230785_dp
    aab_pm6( 78,  8)=        2.1105630398_dp
    aab_pm6( 78,  9)=        3.7144410610_dp
    aab_pm6( 78, 13)=        1.5723600388_dp
    aab_pm6( 78, 14)=        0.9999899864_dp
    aab_pm6( 78, 15)=        1.4032390118_dp
    aab_pm6( 78, 16)=        2.7915000916_dp
    aab_pm6( 78, 17)=        2.1085259914_dp
    aab_pm6( 78, 35)=        2.1853070259_dp
    aab_pm6( 78, 53)=        3.0773379803_dp
    aab_pm6( 78, 78)=        3.4042758942_dp
    aab_pm6( 79,  1)=        3.3690409660_dp
    aab_pm6( 79,  6)=        4.5800161362_dp
    aab_pm6( 79,  7)=        2.1380949020_dp
    aab_pm6( 79,  8)=        1.5487630367_dp
    aab_pm6( 79,  9)=        4.4531450272_dp
    aab_pm6( 79, 13)=        1.5725699663_dp
    aab_pm6( 79, 15)=        1.6187130213_dp
    aab_pm6( 79, 16)=        4.3062381744_dp
    aab_pm6( 79, 17)=        3.5394139290_dp
    aab_pm6( 79, 35)=        0.5819110274_dp
    aab_pm6( 79, 53)=        0.5779160261_dp
    aab_pm6( 79, 79)=        0.9031620026_dp
    aab_pm6( 80,  1)=        1.1365870237_dp
    aab_pm6( 80,  6)=        0.7958160043_dp
    aab_pm6( 80,  7)=        0.3321520090_dp
    aab_pm6( 80,  8)=        1.0521450043_dp
    aab_pm6( 80,  9)=        1.2405719757_dp
    aab_pm6( 80, 14)=        2.7708599567_dp
    aab_pm6( 80, 15)=        0.6086040139_dp
    aab_pm6( 80, 16)=        1.0416820049_dp
    aab_pm6( 80, 17)=        0.4307309985_dp
    aab_pm6( 80, 22)=        3.4146299362_dp
    aab_pm6( 80, 35)=        0.6387169957_dp
    aab_pm6( 80, 52)=        0.2915000021_dp
    aab_pm6( 80, 53)=        0.7581620216_dp
    aab_pm6( 80, 80)=        0.4744130075_dp
    aab_pm6( 81,  1)=        0.6736580133_dp
    aab_pm6( 81,  5)=        1.5283470154_dp
    aab_pm6( 81,  6)=        1.3903449774_dp
    aab_pm6( 81,  7)=        0.9823349714_dp
    aab_pm6( 81,  8)=        1.5500680208_dp
    aab_pm6( 81,  9)=        1.4695160389_dp
    aab_pm6( 81, 16)=        0.9948509932_dp
    aab_pm6( 81, 17)=        0.8461930156_dp
    aab_pm6( 81, 35)=        0.8744189739_dp
    aab_pm6( 81, 53)=        0.9020119905_dp
    aab_pm6( 81, 81)=        1.1916840076_dp
    aab_pm6( 82,  1)=        1.5226759911_dp
    aab_pm6( 82,  3)=        1.0018099546_dp
    aab_pm6( 82,  5)=        0.9111970067_dp
    aab_pm6( 82,  6)=        1.5255930424_dp
    aab_pm6( 82,  7)=        1.3173940182_dp
    aab_pm6( 82,  8)=        1.7632100582_dp
    aab_pm6( 82,  9)=        3.2889020443_dp
    aab_pm6( 82, 15)=        4.5167999268_dp
    aab_pm6( 82, 16)=        1.0275189877_dp
    aab_pm6( 82, 17)=        1.0941230059_dp
    aab_pm6( 82, 23)=        1.5000000000_dp
    aab_pm6( 82, 24)=        1.8607599735_dp
    aab_pm6( 82, 30)=        1.5000000000_dp
    aab_pm6( 82, 34)=        2.0000000000_dp
    aab_pm6( 82, 35)=        0.8655499816_dp
    aab_pm6( 82, 41)=        1.5000000000_dp
    aab_pm6( 82, 42)=        2.0000000000_dp
    aab_pm6( 82, 52)=        1.0025589466_dp
    aab_pm6( 82, 53)=        0.9834740162_dp
    aab_pm6( 82, 82)=        1.8817640543_dp
    aab_pm6( 83,  1)=        1.6799050570_dp
    aab_pm6( 83,  3)=        0.3401400149_dp
    aab_pm6( 83,  6)=        1.5340249538_dp
    aab_pm6( 83,  7)=        1.1438759565_dp
    aab_pm6( 83,  8)=        1.5532970428_dp
    aab_pm6( 83,  9)=        2.3554000854_dp
    aab_pm6( 83, 16)=        1.4668790102_dp
    aab_pm6( 83, 17)=        1.2729749680_dp
    aab_pm6( 83, 34)=        1.3447459936_dp
    aab_pm6( 83, 35)=        1.1462329626_dp
    aab_pm6( 83, 53)=        1.3021709919_dp
    aab_pm6( 83, 83)=        1.0740640163_dp

    DO i=0,nelem
       DO j=i+1,nelem
          xab = xab_pm6(i,j)+xab_pm6(j,i)
          xab_pm6(i,j) = xab
          xab_pm6(j,i) = xab
          aab = aab_pm6(i,j)+aab_pm6(j,i)
          aab_pm6(i,j) = aab
          aab_pm6(j,i) = aab
       END DO
    END DO

    init = 1

  END SUBROUTINE init_pm6_pair_params

END MODULE semi_empirical_parameters

