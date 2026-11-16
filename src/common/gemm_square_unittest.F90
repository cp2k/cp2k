!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM gemm_square_unittest
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: gemm_square
#include "../base/base_uses.f90"

   IMPLICIT NONE

   COMPLEX(kind=dp), DIMENSION(3, 3) :: A_in, B_in, C_in, res_c, res_c_ref
   REAL(kind=dp), DIMENSION(3, 3) :: X_in, Y_in, Z_in, res_r, res_r_ref
   REAL(kind=dp) :: tolerance = 1.0e-6_dp

   ! Prepare inputs
   A_in(1, 1) = CMPLX(0.8815928086074307_dp, 0.6726432190297216_dp, kind=dp)
   A_in(1, 2) = CMPLX(0.7660079579530265_dp, 0.6663301208376479_dp, kind=dp)
   A_in(1, 3) = CMPLX(0.8910730680466552_dp, 0.6447684662974965_dp, kind=dp)
   A_in(2, 1) = CMPLX(0.270178070784315_dp, 0.9380895276020503_dp, kind=dp)
   A_in(2, 2) = CMPLX(0.4365740872106577_dp, 0.5843460996868933_dp, kind=dp)
   A_in(2, 3) = CMPLX(0.07466461985206008_dp, 0.6899750234684598_dp, kind=dp)
   A_in(3, 1) = CMPLX(0.840974290337725_dp, 0.8395064317543346_dp, kind=dp)
   A_in(3, 2) = CMPLX(0.5872667752635958_dp, 0.6233467352665024_dp, kind=dp)
   A_in(3, 3) = CMPLX(0.5024930933188588_dp, 0.7727803824712417_dp, kind=dp)

   B_in(1, 1) = CMPLX(0.08269815253296364_dp, 0.34184561260312574_dp, kind=dp)
   B_in(1, 2) = CMPLX(0.9876346392802493_dp, 0.26436123295003866_dp, kind=dp)
   B_in(1, 3) = CMPLX(0.780810836185207_dp, 0.376036133357872_dp, kind=dp)
   B_in(2, 1) = CMPLX(0.4787818411690774_dp, 0.7596241356044092_dp, kind=dp)
   B_in(2, 2) = CMPLX(0.4298758196722595_dp, 0.4813479548810141_dp, kind=dp)
   B_in(2, 3) = CMPLX(0.2086685419449945_dp, 0.3860478932514133_dp, kind=dp)
   B_in(3, 1) = CMPLX(0.34008386216308817_dp, 0.8353095227337101_dp, kind=dp)
   B_in(3, 2) = CMPLX(0.7379600798045334_dp, 0.7634442366598211_dp, kind=dp)
   B_in(3, 3) = CMPLX(0.9840849653895581_dp, 0.9273454280026875_dp, kind=dp)

   C_in(1, 1) = CMPLX(0.731192921191078_dp, 0.9732725403607281_dp, kind=dp)
   C_in(1, 2) = CMPLX(0.07386957805916261_dp, 0.14228952898305391_dp, kind=dp)
   C_in(1, 3) = CMPLX(0.12229506342104235_dp, 0.6298697123856768_dp, kind=dp)
   C_in(2, 1) = CMPLX(0.007352653494114958_dp, 0.29359318766569575_dp, kind=dp)
   C_in(2, 2) = CMPLX(0.29087841717040863_dp, 0.48194825561460775_dp, kind=dp)
   C_in(2, 3) = CMPLX(0.22558916232632764_dp, 0.9229223568661166_dp, kind=dp)
   C_in(3, 1) = CMPLX(0.5728946948517463_dp, 0.9149335302204014_dp, kind=dp)
   C_in(3, 2) = CMPLX(0.20475976494474424_dp, 0.6082208447082643_dp, kind=dp)
   C_in(3, 3) = CMPLX(0.9060121198373113_dp, 0.008565705864987172_dp, kind=dp)

   X_in(1, 1) = 0.42929014430726375_dp
   X_in(1, 2) = 0.21820709659663573_dp
   X_in(1, 3) = 0.5394292090282415_dp
   X_in(2, 1) = 0.7828031363115503_dp
   X_in(2, 2) = 0.1422677264194132_dp
   X_in(2, 3) = 0.25344520034350637_dp
   X_in(3, 1) = 0.5044049742159297_dp
   X_in(3, 2) = 0.6969177100349894_dp
   X_in(3, 3) = 0.6999162742203425_dp

   Y_in(1, 1) = 0.5331333823513378_dp
   Y_in(1, 2) = 0.8001773249628732_dp
   Y_in(1, 3) = 0.2850504760853374_dp
   Y_in(2, 1) = 0.23062673571851455_dp
   Y_in(2, 2) = 0.5013417881822918_dp
   Y_in(2, 3) = 0.07530315834987644_dp
   Y_in(3, 1) = 0.2267846125008932_dp
   Y_in(3, 2) = 0.19831160340777076_dp
   Y_in(3, 3) = 0.3050258528838238_dp

   Z_in(1, 1) = 0.5400800562659297_dp
   Z_in(1, 2) = 0.506259700373107_dp
   Z_in(1, 3) = 0.24342576996957088_dp
   Z_in(2, 1) = 0.3517364012861689_dp
   Z_in(2, 2) = 0.04901381134580918_dp
   Z_in(2, 3) = 0.31263102401008236_dp
   Z_in(3, 1) = 0.20684120795408456_dp
   Z_in(3, 2) = 0.8051322416754273_dp
   Z_in(3, 3) = 0.5860282518273413_dp

   ! Test X * Y

   CALL gemm_square(X_in, 'N', Y_in, 'N', res_r)

   res_r_ref(1, 1) = 0.4015275411844552_dp
   res_r_ref(1, 2) = 0.5598796466739117_dp
   res_r_ref(1, 3) = 0.3033408981158978_dp
   res_r_ref(2, 1) = 0.5076266966693295_dp
   res_r_ref(2, 2) = 0.7479672000061859_dp
   res_r_ref(2, 3) = 0.31115895421143025_dp
   res_r_ref(3, 1) = 0.588373227540499_dp
   res_r_ref(3, 2) = 0.8918089125227482_dp
   res_r_ref(3, 3) = 0.4097535412069894_dp

   CALL check_ref_r(res_r, res_r_ref, tolerance)

   ! Test A * B

   CALL gemm_square(A_in, 'N', B_in, 'N', res_c)

   res_c_ref(1, 1) = CMPLX(-0.5319854477298264_dp, 2.221497049670068_dp, kind=dp)
   res_c_ref(1, 2) = CMPLX(0.8667540454951561_dp, 2.708638262772094_dp, kind=dp)
   res_c_ref(1, 3) = CMPLX(0.6169940071822635_dp, 2.7523152479106017_dp, kind=dp)
   res_c_ref(2, 1) = CMPLX(-1.0841486927833452_dp, 1.0783614128260843_dp, kind=dp)
   res_c_ref(2, 2) = CMPLX(-0.5464163853751492_dp, 2.025430919812893_dp, kind=dp)
   res_c_ref(2, 3) = CMPLX(-0.8426527494261556_dp, 1.872774281691093_dp, kind=dp)
   res_c_ref(3, 1) = CMPLX(-0.884392147923063_dp, 1.78400551956788_dp, kind=dp)
   res_c_ref(3, 2) = CMPLX(0.3418926087394045_dp, 2.5559945109530244_dp, kind=dp)
   res_c_ref(3, 3) = CMPLX(0.000721037947416292_dp, 2.5549846237243488_dp, kind=dp)

   CALL check_ref_c(res_c, res_c_ref, tolerance)

   ! Test X * Y * Z

   CALL gemm_square(X_in, 'N', Y_in, 'N', Z_in, 'N', res_r)

   res_r_ref(1, 1) = 0.4765304668978436_dp
   res_r_ref(1, 2) = 0.47494858536191625_dp
   res_r_ref(1, 3) = 0.45054423436947794_dp
   res_r_ref(2, 1) = 0.6016068400643494_dp
   res_r_ref(2, 2) = 0.5441757689127916_dp
   res_r_ref(2, 3) = 0.5397551091346775_dp
   res_r_ref(3, 1) = 0.7162042207878401_dp
   res_r_ref(3, 2) = 0.6714863948435401_dp
   res_r_ref(3, 3) = 0.6621594909204267_dp

   CALL check_ref_r(res_r, res_r_ref, tolerance)

   ! Test A * B * C

   CALL gemm_square(A_in, 'N', B_in, 'N', C_in, 'N', res_c)

   res_c_ref(1, 1) = CMPLX(-5.504683782712595_dp, 3.5222601599484484_dp, kind=dp)
   res_c_ref(1, 2) = CMPLX(-2.9563767075011937_dp, 2.232852141215477_dp, kind=dp)
   res_c_ref(1, 3) = CMPLX(-3.233216866606864_dp, 3.8464986862655572_dp, kind=dp)
   res_c_ref(2, 1) = CMPLX(-4.63714700646176_dp, -0.11028256160215588_dp, kind=dp)
   res_c_ref(2, 2) = CMPLX(-2.680220510420048_dp, 0.12215466665130881_dp, kind=dp)
   res_c_ref(2, 3) = CMPLX(-3.583889556370547_dp, 1.091159499729193_dp, kind=dp)
   res_c_ref(3, 1) = CMPLX(-5.468121643036268_dp, 2.0272651357548415_dp, kind=dp)
   res_c_ref(3, 2) = CMPLX(-3.0054301614279586_dp, 1.4377987781538424_dp, kind=dp)
   res_c_ref(3, 3) = CMPLX(-3.534937025955706_dp, 2.8681214444912606_dp, kind=dp)

   CALL check_ref_c(res_c, res_c_ref, tolerance)

   ! Test X^T * Y * Z

   CALL gemm_square(X_in, 'T', Y_in, 'N', Z_in, 'N', res_r)

   res_r_ref(1, 1) = 0.6462671475738445_dp
   res_r_ref(1, 2) = 0.5760105624808499_dp
   res_r_ref(1, 3) = 0.5852827099256332_dp
   res_r_ref(2, 1) = 0.36007554144181786_dp
   res_r_ref(2, 2) = 0.40420627668516457_dp
   res_r_ref(2, 3) = 0.36217776140102986_dp
   res_r_ref(3, 1) = 0.5978645617240842_dp
   res_r_ref(3, 2) = 0.600788264751924_dp
   res_r_ref(3, 3) = 0.5673424971463421_dp

   CALL check_ref_r(res_r, res_r_ref, tolerance)

   ! Test A^H * B * C

   CALL gemm_square(A_in, 'C', B_in, 'N', C_in, 'N', res_c)

   res_c_ref(1, 1) = CMPLX(3.375089298965469_dp, 5.744913993063936_dp, kind=dp)
   res_c_ref(1, 2) = CMPLX(2.0725172551868294_dp, 3.258926327791143_dp, kind=dp)
   res_c_ref(1, 3) = CMPLX(3.965529787950442_dp, 3.621340775428089_dp, kind=dp)
   res_c_ref(2, 1) = CMPLX(2.4231309591599897_dp, 4.665551869666368_dp, kind=dp)
   res_c_ref(2, 2) = CMPLX(1.5937647760286848_dp, 2.6021783330446246_dp, kind=dp)
   res_c_ref(2, 3) = CMPLX(2.9609793918714686_dp, 2.92153954960111_dp, kind=dp)
   res_c_ref(3, 1) = CMPLX(3.278689562249669_dp, 4.308656958132163_dp, kind=dp)
   res_c_ref(3, 2) = CMPLX(2.05357432643753_dp, 2.5060755291807237_dp, kind=dp)
   res_c_ref(3, 3) = CMPLX(3.646272530313196_dp, 2.5667051324585874_dp, kind=dp)

   CALL check_ref_c(res_c, res_c_ref, tolerance)

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param mat ...
!> \param ref ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check_ref_r(mat, ref, tolerance)
      REAL(kind=dp), DIMENSION(3, 3)                     :: mat, ref
      REAL(kind=dp)                                      :: tolerance

      INTEGER                                            :: i, j

      DO i = 1, 3
         DO j = 1, 3
            CPASSERT(ABS(mat(i, j) - ref(i, j)) <= tolerance)
         END DO
      END DO
   END SUBROUTINE check_ref_r
! **************************************************************************************************
!> \brief ...
!> \param mat ...
!> \param ref ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check_ref_c(mat, ref, tolerance)
      COMPLEX(kind=dp), DIMENSION(3, 3)                  :: mat, ref
      REAL(kind=dp)                                      :: tolerance

      INTEGER                                            :: i, j

      DO i = 1, 3
         DO j = 1, 3
            CPASSERT(ABS(mat(i, j) - ref(i, j)) <= tolerance)
         END DO
      END DO
   END SUBROUTINE check_ref_c
END PROGRAM gemm_square_unittest
