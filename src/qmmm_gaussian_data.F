!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_gaussian_data [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_data
!!
!!   FUNCTION
!!     GEEP coefficients for fixed radii. THIS IS THE GEEP LIB.
!!     They can be used together with the "omotetia" rule for erf function
!!     to obtain the expansion for whatever radius one desires 
!!     Radius value are in Angstrom
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     7.2005 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_gaussian_data

  USE kinds,                           ONLY: dp  

  IMPLICIT NONE
  ! G5
  REAL(KIND=dp), PARAMETER, PUBLIC :: g5_rc = 0.9_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g5_a1 = 0.0728205_dp, g5_b1 = 4.98419_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g5_a2 = 0.0797713_dp, g5_b2 = 1.78985_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g5_a3 = 0.0885763_dp, g5_b3 = 3.68614_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g5_a4 = 0.1104310_dp, g5_b4 = 2.78214_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g5_a5 = 0.1211280_dp, g5_b5 = 2.15508_dp
  ! G6
  REAL(KIND=dp), PARAMETER, PUBLIC :: g6_rc = 0.9_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g6_a1 = 0.0863850_dp, g6_b1 = 3.20872_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g6_a2 = 0.0704031_dp, g6_b2 = 4.12508_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g6_a3 = 0.0604464_dp, g6_b3 = 5.40108_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g6_a4 = 0.1029610_dp, g6_b4 = 2.04751_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g6_a5 = 0.1020950_dp, g6_b5 = 2.52820_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g6_a6 = 0.0621064_dp, g6_b6 = 1.76831_dp
  ! G7
  REAL(KIND=dp), PARAMETER, PUBLIC :: g7_rc = 0.85_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g7_a1 = 0.0969005_dp, g7_b1 = 1.88076_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g7_a2 = 0.0556821_dp, g7_b2 = 1.65955_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g7_a3 = 0.0607934_dp, g7_b3 = 4.42485_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g7_a4 = 0.0888440_dp, g7_b4 = 2.79941_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g7_a5 = 0.0535919_dp, g7_b5 = 5.68373_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g7_a6 = 0.0736222_dp, g7_b6 = 3.50534_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g7_a7 = 0.1011310_dp, g7_b7 = 2.26261_dp
  ! G8
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_rc = 0.80_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_a1 = 0.0557370_dp, g8_b1 = 1.55883_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_a2 = 0.0450800_dp, g8_b2 = 6.02024_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_a3 = 0.0496262_dp, g8_b3 = 4.82212_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_a4 = 0.0757863_dp, g8_b4 = 3.16120_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_a5 = 0.0608617_dp, g8_b5 = 3.90652_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_a6 = 0.1040810_dp, g8_b6 = 2.09075_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_a7 = 0.0983610_dp, g8_b7 = 1.75423_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g8_a8 = 0.0920273_dp, g8_b8 = 2.55933_dp
  ! G9
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_rc = 0.70_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a1 = 0.1152140_dp, g9_b1 = 1.54224_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a2 = 0.0280838_dp, g9_b2 = 5.91022_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a3 = 0.0718097_dp, g9_b3 = 3.54224_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a4 = 0.0585113_dp, g9_b4 = 4.44860_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a5 = 0.1068470_dp, g9_b5 = 2.27491_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a6 = 0.0883974_dp, g9_b6 = 2.83230_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a7 = 0.1211100_dp, g9_b7 = 1.84735_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a8 = 0.0273353_dp, g9_b8 = 5.52297_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g9_a9 = 0.0657478_dp, g9_b9 = 1.36539_dp
  ! G10
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_rc  = 0.60_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a1  = 0.1005460_dp, g10_b1  = 2.49407_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a2  = 0.1246370_dp, g10_b2  = 1.99718_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a3  = 0.0492434_dp, g10_b3  = 4.70767_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a4  = 0.0785904_dp, g10_b4  = 3.10735_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a5  = 0.0242791_dp, g10_b5  = 5.81268_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a6  = 0.0229056_dp, g10_b6  = 5.92886_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a7  = 0.0606229_dp, g10_b7  = 3.83892_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a8  = 0.1442920_dp, g10_b8  = 1.61113_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a9  = 0.0814189_dp, g10_b9  = 1.17276_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g10_a10 = 0.1401800_dp, g10_b10 = 1.33402_dp
  ! G11
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_rc  = 0.50_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a1  = 0.1041690_dp, g11_b1  = 0.979627_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a2  = 0.0271666_dp, g11_b2  = 5.799740_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a3  = 0.0200111_dp, g11_b3  = 5.799660_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a4  = 0.0483906_dp, g11_b4  = 4.203330_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a5  = 0.1777430_dp, g11_b5  = 1.370190_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a6  = 0.0731972_dp, g11_b6  = 3.415900_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a7  = 0.0950076_dp, g11_b7  = 2.724470_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a8  = 0.0295802_dp, g11_b8  = 4.825220_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a9  = 0.1212380_dp, g11_b9  = 2.162170_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a10 = 0.1513510_dp, g11_b10 = 1.714860_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g11_a11 = 0.1760780_dp, g11_b11 = 1.123400_dp 
  ! G12
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_rc  = 0.40_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a1  = 0.1420750_dp, g12_b1  = 0.786470_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a2  = 0.0880114_dp, g12_b2  = 2.942550_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a3  = 0.0645509_dp, g12_b3  = 3.669640_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a4  = 0.1505330_dp, g12_b4  = 1.826690_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a5  = 0.0183100_dp, g12_b5  = 5.965520_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a6  = 0.0182645_dp, g12_b6  = 5.965520_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a7  = 0.1160340_dp, g12_b7  = 2.326450_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a8  = 0.0449121_dp, g12_b8  = 4.458350_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a9  = 0.0264801_dp, g12_b9  = 5.182710_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a10 = 0.2338230_dp, g12_b10 = 0.912631_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a11 = 0.2295300_dp, g12_b11 = 1.128850_dp 
  REAL(KIND=dp), PARAMETER, PUBLIC :: g12_a12 = 0.1911980_dp, g12_b12 = 1.431550_dp 
  ! G13
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_rc  = 0.30_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a1  = 0.2001960_dp, g13_b1  = 0.591285_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a2  = 0.0884682_dp, g13_b2  = 2.958380_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a3  = 0.1999450_dp, g13_b3  = 1.420640_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a4  = 0.3235490_dp, g13_b4  = 0.691636_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a5  = 0.0636926_dp, g13_b5  = 3.693610_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a6  = 0.0409795_dp, g13_b6  = 4.450810_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a7  = 0.1172330_dp, g13_b7  = 2.332480_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a8  = 0.0278787_dp, g13_b8  = 5.119550_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a9  = 0.0187190_dp, g13_b9  = 5.952110_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a10 = 0.0187159_dp, g13_b10 = 5.952110_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a11 = 0.3118340_dp, g13_b11 = 0.863435_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a12 = 0.2562320_dp, g13_b12 = 1.104410_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g13_a13 = 0.1535690_dp, g13_b13 = 1.824620_dp
  ! G14
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_rc  = 0.20_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a1  = 0.3778350_dp, g14_b1  = 0.786779_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a2  = 0.0944711_dp, g14_b2  = 2.781090_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a3  = 0.0204159_dp, g14_b3  = 5.876050_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a4  = 0.0204239_dp, g14_b4  = 5.876050_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a5  = 0.2884320_dp, g14_b5  = 1.023350_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a6  = 0.2186240_dp, g14_b6  = 1.326440_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a7  = 0.0700271_dp, g14_b7  = 3.483220_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a8  = 0.1656890_dp, g14_b8  = 1.710300_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a9  = 0.0493301_dp, g14_b9  = 4.268870_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a10 = 0.1254840_dp, g14_b10 = 2.190730_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a11 = 0.0303071_dp, g14_b11 = 4.999930_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a12 = 0.5269670_dp, g14_b12 = 0.475269_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a13 = 0.3496170_dp, g14_b13 = 0.397232_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g14_a14 = 0.4776640_dp, g14_b14 = 0.605916_dp
  ! G15
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_rc  = 0.10_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a1  = 1.3587200_dp, g15_b1  = 0.309966_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a2  = 0.0273669_dp, g15_b2  = 5.820070_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a3  = 0.2845940_dp, g15_b3  = 1.061340_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a4  = 0.1599350_dp, g15_b4  = 1.784880_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a5  = 0.8510690_dp, g15_b5  = 0.441778_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a6  = 0.5579850_dp, g15_b6  = 0.605330_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a7  = 0.0256516_dp, g15_b7  = 4.901890_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a8  = 0.3903200_dp, g15_b8  = 0.807711_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a9  = 0.0164043_dp, g15_b9  = 5.819900_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a10 = 0.0368220_dp, g15_b10 = 4.412630_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a11 = 0.1211020_dp, g15_b11 = 2.290750_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a12 = 0.2121480_dp, g15_b12 = 1.381320_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a13 = 0.0662304_dp, g15_b13 = 3.657730_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a14 = 0.0912135_dp, g15_b14 = 2.915640_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a15 = 1.5942000_dp, g15_b15 = 0.216197_dp

END MODULE qmmm_gaussian_data
