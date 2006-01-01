!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


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
  INTEGER, PARAMETER, PUBLIC       :: min_geep_lib_gauss = 2
  INTEGER, PARAMETER, PUBLIC       :: max_geep_lib_gauss = 18
  ! G2
  REAL(KIND=dp), PARAMETER, PUBLIC :: g2_rc = 1.30_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g2_a1 = 0.121492_dp, g2_b1 = 3.85695_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g2_a2 = 0.104509_dp, g2_b2 = 2.72971_dp
  ! G3
  REAL(KIND=dp), PARAMETER, PUBLIC :: g3_rc = 1.17_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g3_a1 = 0.1142400_dp, g3_b1 = 3.08635_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g3_a2 = 0.0867876_dp, g3_b2 = 2.38345_dp  
  REAL(KIND=dp), PARAMETER, PUBLIC :: g3_a3 = 0.0971514_dp, g3_b3 = 4.35065_dp  
  ! G4
  REAL(KIND=dp), PARAMETER, PUBLIC :: g4_rc = 1.05_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g4_a1 = 0.0776927_dp, g4_b1 = 2.10429_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g4_a2 = 0.1123350_dp, g4_b2 = 2.59523_dp  
  REAL(KIND=dp), PARAMETER, PUBLIC :: g4_a3 = 0.0812506_dp, g4_b3 = 4.73860_dp  
  REAL(KIND=dp), PARAMETER, PUBLIC :: g4_a4 = 0.0989222_dp, g4_b4 = 3.44579_dp  
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
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a1  = 1.3606800_dp, g15_b1  = 0.310940_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a2  = 0.0266233_dp, g15_b2  = 5.801190_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a3  = 0.3896280_dp, g15_b3  = 0.813967_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a4  = 0.0155094_dp, g15_b4  = 5.800280_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a5  = 0.5574040_dp, g15_b5  = 0.609044_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a6  = 0.8507370_dp, g15_b6  = 0.443840_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a7  = 0.1590220_dp, g15_b7  = 1.808170_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a8  = 0.2837530_dp, g15_b8  = 1.071380_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a9  = 0.0645391_dp, g15_b9  = 3.723390_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a10 = 0.0241980_dp, g15_b10 = 4.946500_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a11 = 0.0344342_dp, g15_b11 = 4.462030_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a12 = 0.2112170_dp, g15_b12 = 1.396840_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a13 = 0.1202950_dp, g15_b13 = 2.325120_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a14 = 0.0904425_dp, g15_b14 = 2.965550_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g15_a15 = 1.6040000_dp, g15_b15 = 0.216434_dp
  ! G16
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_rc  = 0.10_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a1  = 1.5843100_dp, g16_b1  = 0.215959_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a2  = 0.0253725_dp, g16_b2  = 5.737070_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a3  = 0.1191710_dp, g16_b3  = 2.233620_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a4  = 0.2848860_dp, g16_b4  = 1.050490_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a5  = 0.0134257_dp, g16_b5  = 5.720440_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a6  = 0.0269327_dp, g16_b6  = 4.542090_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a7  = 0.0888076_dp, g16_b7  = 2.807370_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a8  = 0.3907140_dp, g16_b8  = 0.801210_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a9  = 0.5584020_dp, g16_b9  = 0.601540_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a10 = 0.0674428_dp, g16_b10 = 3.477020_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a11 = 0.0167448_dp, g16_b11 = 5.320770_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a12 = 0.2120830_dp, g16_b12 = 1.363230_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a13 = 0.0363796_dp, g16_b13 = 4.204950_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a14 = 0.8512990_dp, g16_b14 = 0.439692_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a15 = 0.1591170_dp, g16_b15 = 1.753690_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g16_a16 = 1.3566500_dp, g16_b16 = 0.308985_dp
  ! G17
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_rc  = 0.10_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a1  = 0.0286263_dp, g17_b1  = 5.897570_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a2  = 1.3621700_dp, g17_b2  = 0.210930_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a3  = 0.8582340_dp, g17_b3  = 0.396804_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a4  = 0.0918589_dp, g17_b4  = 2.514930_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a5  = 0.2314390_dp, g17_b5  = 1.077620_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a6  = 1.3030900_dp, g17_b6  = 0.288252_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a7  = 0.1437920_dp, g17_b7  = 1.658970_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a8  = 0.3010060_dp, g17_b8  = 0.860535_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a9  = 0.0274321_dp, g17_b9  = 5.099500_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a10 = 0.1814220_dp, g17_b10 = 1.340220_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a11 = 0.4046220_dp, g17_b11 = 0.679621_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a12 = 0.0735338_dp, g17_b12 = 3.080580_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a13 = 0.0577419_dp, g17_b13 = 3.752150_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a14 = 0.0388557_dp, g17_b14 = 4.481620_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a15 = 0.1147140_dp, g17_b15 = 2.045960_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a16 = 0.5719720_dp, g17_b16 = 0.526925_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g17_a17 = 0.0222331_dp, g17_b17 = 6.674390_dp
  ! G18
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_rc  = 0.10_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a1  = 0.0326473_dp, g18_b1  = 5.962800_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a2  = 0.0181199_dp, g18_b2  = 6.674950_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a3  = 0.8582970_dp, g18_b3  = 0.396103_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a4  = 0.0880221_dp, g18_b4  = 2.468750_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a5  = 0.1812460_dp, g18_b5  = 1.332910_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a6  = 0.2314380_dp, g18_b6  = 1.072960_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a7  = 0.0449632_dp, g18_b7  = 3.435210_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a8  = 1.3583000_dp, g18_b8  = 0.210847_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a9  = 0.1132310_dp, g18_b9  = 2.024200_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a10 = 0.0251085_dp, g18_b10 = 4.226690_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a11 = 0.4047630_dp, g18_b11 = 0.677704_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a12 = 0.3010870_dp, g18_b12 = 0.857535_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a13 = 1.3019600_dp, g18_b13 = 0.287909_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a14 = 0.1432190_dp, g18_b14 = 1.646940_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a15 = 0.0634725_dp, g18_b15 = 2.959960_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a16 = 0.5721690_dp, g18_b16 = 0.525730_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a17 = 0.0309731_dp, g18_b17 = 3.918400_dp
  REAL(KIND=dp), PARAMETER, PUBLIC :: g18_a18 = 0.0429091_dp, g18_b18 = 4.928750_dp

END MODULE qmmm_gaussian_data


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
