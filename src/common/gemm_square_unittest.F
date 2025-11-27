!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
   A_in(1, 1) = CMPLX(0.8815928086074307, 0.6726432190297216, kind=dp)
   A_in(1, 2) = CMPLX(0.7660079579530265, 0.6663301208376479, kind=dp)
   A_in(1, 3) = CMPLX(0.8910730680466552, 0.6447684662974965, kind=dp)
   A_in(2, 1) = CMPLX(0.270178070784315, 0.9380895276020503, kind=dp)
   A_in(2, 2) = CMPLX(0.4365740872106577, 0.5843460996868933, kind=dp)
   A_in(2, 3) = CMPLX(0.07466461985206008, 0.6899750234684598, kind=dp)
   A_in(3, 1) = CMPLX(0.840974290337725, 0.8395064317543346, kind=dp)
   A_in(3, 2) = CMPLX(0.5872667752635958, 0.6233467352665024, kind=dp)
   A_in(3, 3) = CMPLX(0.5024930933188588, 0.7727803824712417, kind=dp)

   B_in(1, 1) = CMPLX(0.08269815253296364, 0.34184561260312574, kind=dp)
   B_in(1, 2) = CMPLX(0.9876346392802493, 0.26436123295003866, kind=dp)
   B_in(1, 3) = CMPLX(0.780810836185207, 0.376036133357872, kind=dp)
   B_in(2, 1) = CMPLX(0.4787818411690774, 0.7596241356044092, kind=dp)
   B_in(2, 2) = CMPLX(0.4298758196722595, 0.4813479548810141, kind=dp)
   B_in(2, 3) = CMPLX(0.2086685419449945, 0.3860478932514133, kind=dp)
   B_in(3, 1) = CMPLX(0.34008386216308817, 0.8353095227337101, kind=dp)
   B_in(3, 2) = CMPLX(0.7379600798045334, 0.7634442366598211, kind=dp)
   B_in(3, 3) = CMPLX(0.9840849653895581, 0.9273454280026875, kind=dp)

   C_in(1, 1) = CMPLX(0.731192921191078, 0.9732725403607281, kind=dp)
   C_in(1, 2) = CMPLX(0.07386957805916261, 0.14228952898305391, kind=dp)
   C_in(1, 3) = CMPLX(0.12229506342104235, 0.6298697123856768, kind=dp)
   C_in(2, 1) = CMPLX(0.007352653494114958, 0.29359318766569575, kind=dp)
   C_in(2, 2) = CMPLX(0.29087841717040863, 0.48194825561460775, kind=dp)
   C_in(2, 3) = CMPLX(0.22558916232632764, 0.9229223568661166, kind=dp)
   C_in(3, 1) = CMPLX(0.5728946948517463, 0.9149335302204014, kind=dp)
   C_in(3, 2) = CMPLX(0.20475976494474424, 0.6082208447082643, kind=dp)
   C_in(3, 3) = CMPLX(0.9060121198373113, 0.008565705864987172, kind=dp)

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

   res_c_ref(1, 1) = CMPLX(-0.5319854477298264, 2.221497049670068, kind=dp)
   res_c_ref(1, 2) = CMPLX(0.8667540454951561, 2.708638262772094, kind=dp)
   res_c_ref(1, 3) = CMPLX(0.6169940071822635, 2.7523152479106017, kind=dp)
   res_c_ref(2, 1) = CMPLX(-1.0841486927833452, 1.0783614128260843, kind=dp)
   res_c_ref(2, 2) = CMPLX(-0.5464163853751492, 2.025430919812893, kind=dp)
   res_c_ref(2, 3) = CMPLX(-0.8426527494261556, 1.872774281691093, kind=dp)
   res_c_ref(3, 1) = CMPLX(-0.884392147923063, 1.78400551956788, kind=dp)
   res_c_ref(3, 2) = CMPLX(0.3418926087394045, 2.5559945109530244, kind=dp)
   res_c_ref(3, 3) = CMPLX(0.000721037947416292, 2.5549846237243488, kind=dp)

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

   res_c_ref(1, 1) = CMPLX(-5.504683782712595, 3.5222601599484484, kind=dp)
   res_c_ref(1, 2) = CMPLX(-2.9563767075011937, 2.232852141215477, kind=dp)
   res_c_ref(1, 3) = CMPLX(-3.233216866606864, 3.8464986862655572, kind=dp)
   res_c_ref(2, 1) = CMPLX(-4.63714700646176, -0.11028256160215588, kind=dp)
   res_c_ref(2, 2) = CMPLX(-2.680220510420048, 0.12215466665130881, kind=dp)
   res_c_ref(2, 3) = CMPLX(-3.583889556370547, 1.091159499729193, kind=dp)
   res_c_ref(3, 1) = CMPLX(-5.468121643036268, 2.0272651357548415, kind=dp)
   res_c_ref(3, 2) = CMPLX(-3.0054301614279586, 1.4377987781538424, kind=dp)
   res_c_ref(3, 3) = CMPLX(-3.534937025955706, 2.8681214444912606, kind=dp)

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

   res_c_ref(1, 1) = CMPLX(3.375089298965469, 5.744913993063936, kind=dp)
   res_c_ref(1, 2) = CMPLX(2.0725172551868294, 3.258926327791143, kind=dp)
   res_c_ref(1, 3) = CMPLX(3.965529787950442, 3.621340775428089, kind=dp)
   res_c_ref(2, 1) = CMPLX(2.4231309591599897, 4.665551869666368, kind=dp)
   res_c_ref(2, 2) = CMPLX(1.5937647760286848, 2.6021783330446246, kind=dp)
   res_c_ref(2, 3) = CMPLX(2.9609793918714686, 2.92153954960111, kind=dp)
   res_c_ref(3, 1) = CMPLX(3.278689562249669, 4.308656958132163, kind=dp)
   res_c_ref(3, 2) = CMPLX(2.05357432643753, 2.5060755291807237, kind=dp)
   res_c_ref(3, 3) = CMPLX(3.646272530313196, 2.5667051324585874, kind=dp)

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
