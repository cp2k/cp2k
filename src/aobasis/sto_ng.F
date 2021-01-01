!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      06.2004 created
!> \author jgh
! **************************************************************************************************
MODULE sto_ng

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sto_ng'

   PUBLIC :: get_sto_ng

! **************************************************************************************************
!> \brief return STO-NG parameters;
!>        INPUT: zeta (Slater exponent)
!>               n    (Expansion length)
!>               nq   (principle quantum number)
!>               lq   (angular quantum number)
!>        OUTPUT: alpha (n Gaussian exponents)
!>                coef  (n corresponding expansion coefficients)
!>
!>        Error of expansion: Int ( |Phi - Chi|^2 )
!>
!> \par History
!>      06.2004 created using parameters from : RF Stewart JCP 52 p431 1970
!>      01.2008 added newly generated parameters for 6s - 7i
!>      08.2018 added newly generated parameters for 8s - 9l
!> \author jgh
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param zeta ...
!> \param n ...
!> \param nq ...
!> \param lq ...
!> \param alpha ...
!> \param coef ...
! **************************************************************************************************
   SUBROUTINE get_sto_ng(zeta, n, nq, lq, alpha, coef)

      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: n, nq, lq
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: alpha, coef

      INTEGER                                            :: m, na, nc

      na = SIZE(alpha)
      nc = SIZE(coef)

      CPASSERT(n >= 1 .AND. n <= 6)
      CPASSERT(na >= n)
      CPASSERT(nc >= n)
      CPASSERT(zeta > 0.0_dp)
      CPASSERT(nq >= 1)
      CPASSERT(lq >= 0)

      m = (nq*(nq - 1))/2 + lq + 1
      SELECT CASE (m)
      CASE (1) ! 1s
         SELECT CASE (n)
         CASE (1) ! 4.3191E-2
            alpha(1) = 2.709496091e-1_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 3.1603E-3
            alpha(1) = 8.518186635e-1_dp; coef(1) = 4.301284983e-1_dp
            alpha(2) = 1.516232927e-1_dp; coef(2) = 6.789135305e-1_dp
         CASE (3) ! 3.3053E-4
            alpha(1) = 2.227660584e+0_dp; coef(1) = 1.543289673e-1_dp
            alpha(2) = 4.057711562e-1_dp; coef(2) = 5.353281423e-1_dp
            alpha(3) = 1.098175104e-1_dp; coef(3) = 4.446345422e-1_dp
         CASE (4) ! 4.3763E-5
            alpha(1) = 5.216844534e+0_dp; coef(1) = 5.675242080e-2_dp
            alpha(2) = 9.546182760e-1_dp; coef(2) = 2.601413550e-1_dp
            alpha(3) = 2.652034102e-1_dp; coef(3) = 5.328461143e-1_dp
            alpha(4) = 8.801862774e-2_dp; coef(4) = 2.916254405e-1_dp
         CASE (5) ! 6.8840E-6
            alpha(1) = 1.130563696e+1_dp; coef(1) = 2.214055312e-2_dp
            alpha(2) = 2.071728178e+0_dp; coef(2) = 1.135411520e-1_dp
            alpha(3) = 5.786484833e-1_dp; coef(3) = 3.318161484e-1_dp
            alpha(4) = 1.975724573e-1_dp; coef(4) = 4.825700713e-1_dp
            alpha(5) = 7.445271746e-2_dp; coef(5) = 1.935721966e-1_dp
         CASE (6) ! 1.2372E-6
            alpha(1) = 2.310303149e+1_dp; coef(1) = 9.163596280e-3_dp
            alpha(2) = 4.235915534e+0_dp; coef(2) = 4.936149294e-2_dp
            alpha(3) = 1.185056519e+0_dp; coef(3) = 1.685383049e-1_dp
            alpha(4) = 4.070988982e-1_dp; coef(4) = 3.705627997e-1_dp
            alpha(5) = 1.580884151e-1_dp; coef(5) = 4.164915298e-1_dp
            alpha(6) = 6.510953954e-2_dp; coef(6) = 1.303340841e-1_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (2) ! 2s
         SELECT CASE (n)
         CASE (1) ! 6.2923E-3
            alpha(1) = 1.012151084e-1_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 2.6728E-3
            alpha(1) = 1.292278611e-1_dp; coef(1) = 7.470867124e-1_dp
            alpha(2) = 4.908584205e-2_dp; coef(2) = 2.855980556e-1_dp
         CASE (3) ! 6.8661E-5
            alpha(1) = 2.581578398e+0_dp; coef(1) = -5.994474934e-2_dp
            alpha(2) = 1.567622104e-1_dp; coef(2) = 5.960385398e-1_dp
            alpha(3) = 6.018332272e-2_dp; coef(3) = 4.581786291e-1_dp
         CASE (4) ! 2.6974E-5
            alpha(1) = 1.161525551e+1_dp; coef(1) = -1.198411747e-2_dp
            alpha(2) = 2.000243111e+0_dp; coef(2) = -5.472052539e-2_dp
            alpha(3) = 1.607280687e-1_dp; coef(3) = 5.805587176e-1_dp
            alpha(4) = 6.125744532e-2_dp; coef(4) = 4.770079976e-1_dp
         CASE (5) ! 2.7009E-6
            alpha(1) = 8.984956862e+0_dp; coef(1) = -1.596349096e-2_dp
            alpha(2) = 1.673710636e+0_dp; coef(2) = -5.685884883e-2_dp
            alpha(3) = 1.944726668e-1_dp; coef(3) = 3.698265599e-1_dp
            alpha(4) = 8.806345634e-2_dp; coef(4) = 5.480512593e-1_dp
            alpha(5) = 4.249068522e-2_dp; coef(5) = 1.472634893e-1_dp
         CASE (6) ! 3.2726E-7
            alpha(1) = 2.768496241e+1_dp; coef(1) = -4.151277819e-3_dp
            alpha(2) = 5.077140627e+0_dp; coef(2) = -2.067024148e-2_dp
            alpha(3) = 1.426786050e+0_dp; coef(3) = -5.150303337e-2_dp
            alpha(4) = 2.040335729e-1_dp; coef(4) = 3.346271174e-1_dp
            alpha(5) = 9.260298399e-2_dp; coef(5) = 5.621061301e-1_dp
            alpha(6) = 4.416183978e-2_dp; coef(6) = 1.712994697e-1_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (3) ! 2p
         SELECT CASE (n)
         CASE (1) ! 4.8232E-2
            alpha(1) = 1.759666885e-1_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 3.0947E-3
            alpha(1) = 4.323908358e-1_dp; coef(1) = 4.522627513e-1_dp
            alpha(2) = 1.069439065e-1_dp; coef(2) = 6.713122642e-1_dp
         CASE (3) ! 2.6858E-4
            alpha(1) = 9.192379002e-1_dp; coef(1) = 1.623948553e-1_dp
            alpha(2) = 2.359194503e-1_dp; coef(2) = 5.661708862e-1_dp
            alpha(3) = 8.009805746e-2_dp; coef(3) = 4.223071752e-1_dp
         CASE (4) ! 2.9037E-5
            alpha(1) = 1.798260992e+0_dp; coef(1) = 5.713170255e-2_dp
            alpha(2) = 4.662622228e-1_dp; coef(2) = 2.857455515e-1_dp
            alpha(3) = 1.643718620e-1_dp; coef(3) = 5.517873105e-1_dp
            alpha(4) = 6.543927065e-2_dp; coef(4) = 2.632314924e-1_dp
         CASE (5) ! 3.7171E-6
            alpha(1) = 3.320386533e+0_dp; coef(1) = 2.079051117e-2_dp
            alpha(2) = 8.643257633e-1_dp; coef(2) = 1.235472099e-1_dp
            alpha(3) = 3.079819284e-1_dp; coef(3) = 3.667738886e-1_dp
            alpha(4) = 1.273309895e-1_dp; coef(4) = 4.834930290e-1_dp
            alpha(5) = 5.606243164e-2_dp; coef(5) = 1.653444074e-1_dp
         CASE (6) ! 5.4444E-7
            alpha(1) = 5.868285913e+0_dp; coef(1) = 7.924233646e-3_dp
            alpha(2) = 1.530329631e+0_dp; coef(2) = 5.144104825e-2_dp
            alpha(3) = 5.475665231e-1_dp; coef(3) = 1.898400060e-1_dp
            alpha(4) = 2.288932733e-1_dp; coef(4) = 4.049863191e-1_dp
            alpha(5) = 1.046655969e-1_dp; coef(5) = 4.012362861e-1_dp
            alpha(6) = 4.948220127e-2_dp; coef(6) = 1.051855189e-1_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (4) ! 3s
         SELECT CASE (n)
         CASE (1) ! 1.6764E-2
            alpha(1) = 5.296881757e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 7.6424E-4
            alpha(1) = 6.694095822e-1_dp; coef(1) = -1.529645716e-1_dp
            alpha(2) = 5.837135094e-2_dp; coef(2) = 1.051370110e+0_dp
         CASE (3) ! 8.0718E-5
            alpha(1) = 5.641487709e-1_dp; coef(1) = -1.782577972e-1_dp
            alpha(2) = 6.924421391e-2_dp; coef(2) = 8.612761663e-1_dp
            alpha(3) = 3.269529097e-2_dp; coef(3) = 2.261841969e-1_dp
         CASE (4) ! 1.7254E-6
            alpha(1) = 1.513265591e+0_dp; coef(1) = -3.295496352e-2_dp
            alpha(2) = 4.262497508e-1_dp; coef(2) = -1.724516959e-1_dp
            alpha(3) = 7.643320863e-2_dp; coef(3) = 7.518511194e-1_dp
            alpha(4) = 3.760545063e-2_dp; coef(4) = 3.589627310e-1_dp
         CASE (5) ! 7.9816E-7
            alpha(1) = 4.275877914e+0_dp; coef(1) = -3.920358850e-3_dp
            alpha(2) = 1.132409433e+0_dp; coef(2) = -4.168430506e-2_dp
            alpha(3) = 4.016256968e-1_dp; coef(3) = -1.637440990e-1_dp
            alpha(4) = 7.732370620e-2_dp; coef(4) = 7.419373723e-1_dp
            alpha(5) = 3.800708627e-2_dp; coef(5) = 3.724364929e-1_dp
         CASE (6) ! 4.0662E-8
            alpha(1) = 3.273031938e+0_dp; coef(1) = -6.775596947e-3_dp
            alpha(2) = 9.200611311e-1_dp; coef(2) = -5.639325779e-2_dp
            alpha(3) = 3.593349765e-1_dp; coef(3) = -1.587856086e-1_dp
            alpha(4) = 8.636686991e-2_dp; coef(4) = 5.534527651e-1_dp
            alpha(5) = 4.797373812e-2_dp; coef(5) = 5.015351020e-1_dp
            alpha(6) = 2.724741144e-2_dp; coef(6) = 7.223633674e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (5) ! 3p
         SELECT CASE (n)
         CASE (1) ! 1.2745E-2
            alpha(1) = 9.113614253e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 2.3771E-4
            alpha(1) = 1.458620964e-1_dp; coef(1) = 5.349653144e-1_dp
            alpha(2) = 5.664210742e-2_dp; coef(2) = 5.299607212e-1_dp
         CASE (3) ! 1.3487E-4
            alpha(1) = 2.692880368e+0_dp; coef(1) = -1.061945788e-2_dp
            alpha(2) = 1.489359592e-1_dp; coef(2) = 5.218564264e-1_dp
            alpha(3) = 5.739585040e-2_dp; coef(3) = 5.450015143e-1_dp
         CASE (4) ! 2.9785E-6
            alpha(1) = 1.853180239e+0_dp; coef(1) = -1.434249391e-2_dp
            alpha(2) = 1.915075719e-1_dp; coef(2) = 2.755177589e-1_dp
            alpha(3) = 8.655487938e-2_dp; coef(3) = 5.846750879e-1_dp
            alpha(4) = 4.184253862e-2_dp; coef(4) = 2.144986514e-1_dp
         CASE (5) ! 1.3387E-6
            alpha(1) = 6.466803859e+0_dp; coef(1) = -2.329023747e-3_dp
            alpha(2) = 1.555914802e+0_dp; coef(2) = -1.357395221e-2_dp
            alpha(3) = 1.955925255e-1_dp; coef(3) = 2.632185383e-1_dp
            alpha(4) = 8.809647701e-2_dp; coef(4) = 5.880427024e-1_dp
            alpha(5) = 4.234835707e-2_dp; coef(5) = 2.242794445e-1_dp
         CASE (6) ! 7.9285E-8
            alpha(1) = 5.077973607e+0_dp; coef(1) = -3.329929840e-3_dp
            alpha(2) = 1.340786940e+0_dp; coef(2) = -1.419488340e-2_dp
            alpha(3) = 2.248434849e-1_dp; coef(3) = 1.639395770e-1_dp
            alpha(4) = 1.131741848e-1_dp; coef(4) = 4.485358256e-1_dp
            alpha(5) = 6.076408893e-2_dp; coef(5) = 3.908813050e-1_dp
            alpha(6) = 3.315424265e-2_dp; coef(6) = 7.411456232e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (6) ! 3d
         SELECT CASE (n)
         CASE (1) ! 5.0728E-2
            alpha(1) = 1.302270363e-1_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 2.9834E-3
            alpha(1) = 2.777427345e-1_dp; coef(1) = 4.666137923e-1_dp
            alpha(2) = 8.336507714e-2_dp; coef(2) = 6.644706516e-1_dp
         CASE (3) ! 2.2687E-4
            alpha(1) = 5.229112225e-1_dp; coef(1) = 1.696596060e-1_dp
            alpha(2) = 1.639595876e-1_dp; coef(2) = 5.847984817e-1_dp
            alpha(3) = 6.386630021e-2_dp; coef(3) = 4.056779523e-1_dp
         CASE (4) ! 2.1116E-5
            alpha(1) = 9.185846715e-1_dp; coef(1) = 5.799057705e-2_dp
            alpha(2) = 2.920461109e-1_dp; coef(2) = 3.045581349e-1_dp
            alpha(3) = 1.187568890e-1_dp; coef(3) = 5.601358038e-1_dp
            alpha(4) = 5.286755896e-2_dp; coef(4) = 2.432423313e-1_dp
         CASE (5) ! 2.3115E-6
            alpha(1) = 1.539033958e+0_dp; coef(1) = 2.020869128e-2_dp
            alpha(2) = 4.922090297e-1_dp; coef(2) = 1.321157923e-1_dp
            alpha(3) = 2.029756928e-1_dp; coef(3) = 3.911240346e-1_dp
            alpha(4) = 9.424112917e-2_dp; coef(4) = 4.779609701e-1_dp
            alpha(5) = 4.569058269e-2_dp; coef(5) = 1.463662294e-1_dp
         CASE (6) ! 2.8899E-7
            alpha(1) = 2.488296923e+0_dp; coef(1) = 7.283828112e-3_dp
            alpha(2) = 7.981487853e-1_dp; coef(2) = 5.386799363e-2_dp
            alpha(3) = 3.311327490e-1_dp; coef(3) = 2.072139149e-1_dp
            alpha(4) = 1.559114463e-1_dp; coef(4) = 4.266269092e-1_dp
            alpha(5) = 7.877734732e-2_dp; coef(5) = 3.843100204e-1_dp
            alpha(6) = 4.058484363e-2_dp; coef(6) = 8.902827546e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (7) ! 4s
         SELECT CASE (n)
         CASE (1) ! 4.1760E-2
            alpha(1) = 3.264600274e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 1.7480E-4
            alpha(1) = 2.441785453e-1_dp; coef(1) = -3.046656896e-1_dp
            alpha(2) = 4.051097664e-2_dp; coef(2) = 1.146877294e+0_dp
         CASE (3) ! 1.6867E-5
            alpha(1) = 2.267938753e-1_dp; coef(1) = -3.349048323e-1_dp
            alpha(2) = 4.448178019e-2_dp; coef(2) = 1.056744667e+0_dp
            alpha(3) = 2.195294664e-2_dp; coef(3) = 1.256661680e-1_dp
         CASE (4) ! 8.4105E-7
            alpha(1) = 3.242212833e-1_dp; coef(1) = -1.120682822e-1_dp
            alpha(2) = 1.663217177e-1_dp; coef(2) = -2.845426863e-1_dp
            alpha(3) = 5.081097451e-2_dp; coef(3) = 8.909873788e-1_dp
            alpha(4) = 2.829066600e-2_dp; coef(4) = 3.517811205e-1_dp
         CASE (5) ! 6.0273E-8
            alpha(1) = 2.980263783e+0_dp; coef(1) = 1.513948997e-3_dp
            alpha(2) = 3.792228833e-1_dp; coef(2) = -7.316801548e-2_dp
            alpha(3) = 1.789717224e-1_dp; coef(3) = -3.143703799e-1_dp
            alpha(4) = 5.002110360e-2_dp; coef(4) = 9.032615169e-1_dp
            alpha(5) = 2.789361681e-2_dp; coef(5) = 3.294210848e-1_dp
         CASE (6) ! 5.9700E-9
            alpha(1) = 3.232838646e+0_dp; coef(1) = 1.374817488e-3_dp
            alpha(2) = 3.605788802e-1_dp; coef(2) = -8.666390043e-2_dp
            alpha(3) = 1.717905487e-1_dp; coef(3) = -3.130627309e-1_dp
            alpha(4) = 5.277666487e-2_dp; coef(4) = 7.812787397e-1_dp
            alpha(5) = 3.163400284e-2_dp; coef(5) = 4.389247988e-1_dp
            alpha(6) = 1.874093091e-2_dp; coef(6) = 2.487178756e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (8) ! 4p
         SELECT CASE (n)
         CASE (1) ! 3.9516E-3
            alpha(1) = 5.578350235e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 2.1579E-3
            alpha(1) = 6.190052680e-2_dp; coef(1) = 8.743116767e-1_dp
            alpha(2) = 2.648418407e-2_dp; coef(2) = 1.513640107e-1_dp
         CASE (3) ! 1.2621E-5
            alpha(1) = 4.859692220e-1_dp; coef(1) = -6.147823411e-2_dp
            alpha(2) = 7.430216918e-2_dp; coef(2) = 6.604172234e-1_dp
            alpha(3) = 3.653340923e-2_dp; coef(3) = 3.932639495e-1_dp
         CASE (4) ! 5.4740E-6
            alpha(1) = 1.492607880e+0_dp; coef(1) = -6.035216774e-3_dp
            alpha(2) = 4.327619272e-1_dp; coef(2) = -6.013310874e-2_dp
            alpha(3) = 7.553156064e-2_dp; coef(3) = 6.451518200e-1_dp
            alpha(4) = 3.706272183e-2_dp; coef(4) = 4.117923820e-1_dp
         CASE (5) ! 1.3595E-7
            alpha(1) = 1.091977298e+0_dp; coef(1) = -1.143929558e-2_dp
            alpha(2) = 3.719985051e-1_dp; coef(2) = -6.322651538e-2_dp
            alpha(3) = 8.590019352e-2_dp; coef(3) = 4.398907721e-1_dp
            alpha(4) = 4.786503860e-2_dp; coef(4) = 5.245859166e-1_dp
            alpha(5) = 2.730479990e-2_dp; coef(5) = 1.017072253e-1_dp
         CASE (6) ! 1.3897E-8
            alpha(1) = 2.389722618e+0_dp; coef(1) = -1.665913575e-3_dp
            alpha(2) = 7.960947826e-1_dp; coef(2) = -1.657464971e-2_dp
            alpha(3) = 3.415541380e-1_dp; coef(3) = -5.958513378e-2_dp
            alpha(4) = 8.847434525e-2_dp; coef(4) = 4.053115554e-1_dp
            alpha(5) = 4.958248334e-2_dp; coef(5) = 5.433958189e-1_dp
            alpha(6) = 2.816929784e-2_dp; coef(6) = 1.204970491e-1_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (9) ! 4d
         SELECT CASE (n)
         CASE (1) ! 2.0360E-2
            alpha(1) = 7.941656339e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 3.7488E-4
            alpha(1) = 1.330958892e-1_dp; coef(1) = 4.932764167e-1_dp
            alpha(2) = 5.272119659e-2_dp; coef(2) = 5.918727866e-1_dp
         CASE (3) ! 1.2729E-5
            alpha(1) = 1.777717219e-1_dp; coef(1) = 2.308552718e-1_dp
            alpha(2) = 8.040647350e-2_dp; coef(2) = 6.042409177e-1_dp
            alpha(3) = 3.949855551e-2_dp; coef(3) = 2.595768926e-1_dp
         CASE (4) ! 5.7655E-6
            alpha(1) = 1.995825422e+0_dp; coef(1) = -2.816702620e-3_dp
            alpha(2) = 1.823461280e-1_dp; coef(2) = 2.177095871e-1_dp
            alpha(3) = 8.197240896e-2_dp; coef(3) = 6.058047348e-1_dp
            alpha(4) = 4.000634951e-2_dp; coef(4) = 2.717811257e-1_dp
         CASE (5) ! 1.5122E-7
            alpha(1) = 1.522122079e+0_dp; coef(1) = -3.673711876e-3_dp
            alpha(2) = 2.173041823e-1_dp; coef(2) = 1.167122499e-1_dp
            alpha(3) = 1.084876577e-1_dp; coef(3) = 4.216476416e-1_dp
            alpha(4) = 5.836797641e-2_dp; coef(4) = 4.547673415e-1_dp
            alpha(5) = 3.206682246e-2_dp; coef(5) = 1.037803318e-1_dp
         CASE (6) ! 8.0888E-8
            alpha(1) = 4.634239420e+0_dp; coef(1) = -4.749842876e-4_dp
            alpha(2) = 1.341648295e+0_dp; coef(2) = -3.566777891e-3_dp
            alpha(3) = 2.209593028e-1_dp; coef(3) = 1.108670481e-1_dp
            alpha(4) = 1.101467943e-1_dp; coef(4) = 4.159646930e-1_dp
            alpha(5) = 5.904190370e-2_dp; coef(5) = 4.621672517e-1_dp
            alpha(6) = 3.232628887e-2_dp; coef(6) = 1.081250196e-1_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (10) ! 4f
         SELECT CASE (n)
         CASE (1) ! 5.2215E-2
            alpha(1) = 1.033434062e-1_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 2.8877E-3
            alpha(1) = 2.006693538e-1_dp; coef(1) = 4.769346276e-1_dp
            alpha(2) = 6.865384900e-2_dp; coef(2) = 6.587383976e-1_dp
         CASE (3) ! 1.9904E-4
            alpha(1) = 3.483826963e-1_dp; coef(1) = 1.737856685e-1_dp
            alpha(2) = 1.249380537e-1_dp; coef(2) = 5.973380628e-1_dp
            alpha(3) = 5.349995725e-2_dp; coef(3) = 3.929395614e-1_dp
         CASE (4) ! 1.6508E-5
            alpha(1) = 5.691670217e-1_dp; coef(1) = 5.902730589e-2_dp
            alpha(2) = 2.074585819e-1_dp; coef(2) = 3.191828952e-1_dp
            alpha(3) = 9.298346885e-2_dp; coef(3) = 5.639423893e-1_dp
            alpha(4) = 4.473508853e-2_dp; coef(4) = 2.284796537e-1_dp
         CASE (5) ! 1.5974E-6
            alpha(1) = 8.925960415e-1_dp; coef(1) = 1.999839052e-2_dp
            alpha(2) = 3.277589120e-1_dp; coef(2) = 1.395427440e-1_dp
            alpha(3) = 1.492869962e-1_dp; coef(3) = 4.091508237e-1_dp
            alpha(4) = 7.506099109e-2_dp; coef(4) = 4.708252119e-1_dp
            alpha(5) = 3.892475795e-2_dp; coef(5) = 1.328082566e-1_dp
         CASE (6) ! 1.7590E-7
            alpha(1) = 1.357718039e+0_dp; coef(1) = 6.930234381e-3_dp
            alpha(2) = 5.004907278e-1_dp; coef(2) = 5.634263745e-2_dp
            alpha(3) = 2.296565064e-1_dp; coef(3) = 2.217065797e-1_dp
            alpha(4) = 1.173146814e-1_dp; coef(4) = 4.411388883e-1_dp
            alpha(5) = 6.350097171e-2_dp; coef(5) = 3.688112625e-1_dp
            alpha(6) = 3.474556673e-2_dp; coef(6) = 7.787514504e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (11) ! 5s
         SELECT CASE (n)
         CASE (1) ! 7.1317E-2
            alpha(1) = 2.216912938e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 1.4289E-4
            alpha(1) = 1.213425654e-1_dp; coef(1) = -5.114756049e-1_dp
            alpha(2) = 3.133152144e-2_dp; coef(2) = 1.307377277e+0_dp
         CASE (3) ! 4.1124E-5
            alpha(1) = 1.080198458e-1_dp; coef(1) = -6.617401158e-1_dp
            alpha(2) = 4.408119328e-2_dp; coef(2) = 7.467595004e-1_dp
            alpha(3) = 2.610811810e-2_dp; coef(3) = 7.146490945e-1_dp
         CASE (4) ! 5.4159E-7
            alpha(1) = 8.602284252e-1_dp; coef(1) = 1.103657561e-2_dp
            alpha(2) = 1.189050200e-1_dp; coef(2) = -5.606519023e-1_dp
            alpha(3) = 3.446076176e-2_dp; coef(3) = 1.179429987e+0_dp
            alpha(4) = 1.974798796e-2_dp; coef(4) = 1.734974376e-1_dp
         CASE (5) ! 7.0816E-8
            alpha(1) = 7.403763257e-1_dp; coef(1) = 1.375523371e-2_dp
            alpha(2) = 1.367990863e-1_dp; coef(2) = -3.097344179e-1_dp
            alpha(3) = 9.135301779e-2_dp; coef(3) = -3.199192259e-1_dp
            alpha(4) = 3.726907315e-2_dp; coef(4) = 1.084547038e+0_dp
            alpha(5) = 2.241490836e-2_dp; coef(5) = 3.345288361e-1_dp
         CASE (6) ! 7.9988E-9
            alpha(1) = 1.410128298e+0_dp; coef(1) = 2.695439582e-3_dp
            alpha(2) = 5.077878915e-1_dp; coef(2) = 1.850157487e-2_dp
            alpha(3) = 1.847926858e-1_dp; coef(3) = -9.588628125e-2_dp
            alpha(4) = 1.061070594e-1_dp; coef(4) = -5.200673560e-1_dp
            alpha(5) = 3.669584901e-2_dp; coef(5) = 1.087619490e+0_dp
            alpha(6) = 2.213558430e-2_dp; coef(6) = 3.103964343e-1_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (12) ! 5p
         SELECT CASE (n)
         CASE (1) ! 8.0406E-3
            alpha(1) = 3.769845216e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 6.8912E-4
            alpha(1) = 2.691294191e-1_dp; coef(1) = -1.034227010e-1_dp
            alpha(2) = 3.980805011e-2_dp; coef(2) = 1.033376378e+0_dp
         CASE (3) ! 6.5766E-6
            alpha(1) = 2.127482317e-1_dp; coef(1) = -1.389529695e-1_dp
            alpha(2) = 4.729648620e-2_dp; coef(2) = 8.076691064e-1_dp
            alpha(3) = 2.604865324e-2_dp; coef(3) = 2.726029342e-1_dp
         CASE (4) ! 4.4493E-7
            alpha(1) = 3.962838833e-1_dp; coef(1) = -1.801459207e-2_dp
            alpha(2) = 1.838858552e-1_dp; coef(2) = -1.360777372e-1_dp
            alpha(3) = 4.943555157e-2_dp; coef(3) = 7.533973719e-1_dp
            alpha(4) = 2.750222273e-2_dp; coef(4) = 3.409304859e-1_dp
         CASE (5) ! 1.3483E-8
            alpha(1) = 3.422168934e-1_dp; coef(1) = -3.113958289e-2_dp
            alpha(2) = 1.665099900e-1_dp; coef(2) = -1.374007017e-1_dp
            alpha(3) = 5.443732013e-2_dp; coef(3) = 5.573881018e-1_dp
            alpha(4) = 3.367775277e-2_dp; coef(4) = 4.855428100e-1_dp
            alpha(5) = 2.091949042e-2_dp; coef(5) = 6.605423564e-2_dp
         CASE (6) ! 3.9301E-9
            alpha(1) = 3.778623374e+0_dp; coef(1) = 1.163246387e-4_dp
            alpha(2) = 3.499121109e-1_dp; coef(2) = -2.920771322e-2_dp
            alpha(3) = 1.683175469e-1_dp; coef(3) = -1.381051233e-1_dp
            alpha(4) = 5.404070736e-2_dp; coef(4) = 5.706134877e-1_dp
            alpha(5) = 3.328911801e-2_dp; coef(5) = 4.768808140e-1_dp
            alpha(6) = 2.063815019e-2_dp; coef(6) = 6.021665516e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (13) ! 5d
         SELECT CASE (n)
         CASE (1) ! 6.7308E-3
            alpha(1) = 5.352200793e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 3.1257E-4
            alpha(1) = 6.906014388e-2_dp; coef(1) = 6.539405185e-1_dp
            alpha(2) = 3.399457777e-2_dp; coef(2) = 3.948945302e-1_dp
         CASE (3) ! 2.4664E-5
            alpha(1) = 4.913352950e-1_dp; coef(1) = -2.010175008e-2_dp
            alpha(2) = 7.329090601e-2_dp; coef(2) = 5.899370608e-1_dp
            alpha(3) = 3.594209290e-2_dp; coef(3) = 4.658445960e-1_dp
         CASE (4) ! 1.2737E-6
            alpha(1) = 4.230617826e-1_dp; coef(1) = -2.421626009e-2_dp
            alpha(2) = 8.293863702e-2_dp; coef(2) = 3.937644956e-1_dp
            alpha(3) = 4.590326388e-2_dp; coef(3) = 5.489520286e-1_dp
            alpha(4) = 2.628744797e-2_dp; coef(4) = 1.190436963e-1_dp
         CASE (5) ! 8.6929E-8
            alpha(1) = 9.702946470e-1_dp; coef(1) = -3.231527611e-3_dp
            alpha(2) = 3.603270196e-1_dp; coef(2) = -2.434931372e-2_dp
            alpha(3) = 8.668717752e-2_dp; coef(3) = 3.440817054e-1_dp
            alpha(4) = 4.833708379e-2_dp; coef(4) = 5.693674376e-1_dp
            alpha(5) = 2.751899341e-2_dp; coef(5) = 1.511340183e-1_dp
         CASE (6) ! 1.4086E-8
            alpha(1) = 8.820520428e-1_dp; coef(1) = -4.097377019e-3_dp
            alpha(2) = 3.410838409e-1_dp; coef(2) = -2.508271857e-2_dp
            alpha(3) = 9.204308840e-2_dp; coef(3) = 2.648458555e-1_dp
            alpha(4) = 5.472831774e-2_dp; coef(4) = 5.097437054e-1_dp
            alpha(5) = 3.391202830e-2_dp; coef(5) = 2.654483467e-1_dp
            alpha(6) = 2.108227374e-2_dp; coef(6) = 2.623132212e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (14) ! 5f
         SELECT CASE (n)
         CASE (1) ! 2.6196E-2
            alpha(1) = 6.952785407e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 6.5045E-4
            alpha(1) = 1.156094555e-1_dp; coef(1) = 4.856637346e-1_dp
            alpha(2) = 4.778940916e-2_dp; coef(2) = 6.125980914e-1_dp
         CASE (3) ! 1.5513E-5
            alpha(1) = 1.649233885e-1_dp; coef(1) = 1.909729355e-1_dp
            alpha(2) = 7.487066646e-2_dp; coef(2) = 6.146060459e-1_dp
            alpha(3) = 3.735787219e-2_dp; coef(3) = 3.059611271e-1_dp
         CASE (4) ! 7.5442E-7
            alpha(1) = 2.017831152e-1_dp; coef(1) = 9.174268830e-2_dp
            alpha(2) = 1.001952178e-1_dp; coef(2) = 4.023496947e-1_dp
            alpha(3) = 5.447006630e-2_dp; coef(3) = 4.937432100e-1_dp
            alpha(4) = 3.037569283e-2_dp; coef(4) = 1.254001522e-1_dp
         CASE (5) ! 3.0371E-7
            alpha(1) = 1.670735676e+0_dp; coef(1) = -7.301193568e-4_dp
            alpha(2) = 2.072477219e-1_dp; coef(2) = 8.414991343e-2_dp
            alpha(3) = 1.024709357e-1_dp; coef(3) = 3.923683153e-1_dp
            alpha(4) = 5.537913898e-2_dp; coef(4) = 5.040033146e-1_dp
            alpha(5) = 3.072866652e-2_dp; coef(5) = 1.328979300e-1_dp
         CASE (6) ! 8.4380E-9
            alpha(1) = 1.334096840e+0_dp; coef(1) = -9.486751531e-4_dp
            alpha(2) = 2.372312347e-1_dp; coef(2) = 4.624275998e-2_dp
            alpha(3) = 1.269485744e-1_dp; coef(3) = 2.373699784e-1_dp
            alpha(4) = 7.290318381e-2_dp; coef(4) = 4.589112231e-1_dp
            alpha(5) = 4.351355997e-2_dp; coef(5) = 3.205010548e-1_dp
            alpha(6) = 2.598071843e-2_dp; coef(6) = 5.077063693e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (15) ! 5g
         SELECT CASE (n)
         CASE (1) ! 5.3199E-2
            alpha(1) = 8.565417784e-2_dp; coef(1) = 1.000000000e+0_dp
         CASE (2) ! 2.8109E-3
            alpha(1) = 1.554531559e-1_dp; coef(1) = 4.848298074e-1_dp
            alpha(2) = 5.854079811e-2_dp; coef(2) = 6.539381621e-1_dp
         CASE (3) ! 1.7959E-4
            alpha(1) = 2.545432122e-1_dp; coef(1) = 1.780980905e-1_dp
            alpha(2) = 1.006544376e-1_dp; coef(2) = 6.063757846e-1_dp
            alpha(3) = 4.624463922e-2_dp; coef(3) = 3.828552923e-1_dp
         CASE (4) ! 1.3592E-5
            alpha(1) = 3.945205573e-1_dp; coef(1) = 6.010484250e-2_dp
            alpha(2) = 1.588100623e-1_dp; coef(2) = 3.309738329e-1_dp
            alpha(3) = 7.646521729e-2_dp; coef(3) = 5.655207585e-1_dp
            alpha(4) = 3.898703611e-2_dp; coef(4) = 2.171122608e-1_dp
         CASE (5) ! 1.1901E-6
            alpha(1) = 5.895429375e-1_dp; coef(1) = 1.998085812e-2_dp
            alpha(2) = 2.393343780e-1_dp; coef(2) = 1.460384050e-1_dp
            alpha(3) = 1.172646904e-1_dp; coef(3) = 4.230565459e-1_dp
            alpha(4) = 6.254074479e-2_dp; coef(4) = 4.635699665e-1_dp
            alpha(5) = 3.411243214e-2_dp; coef(5) = 1.226411691e-1_dp
         CASE (6) ! 1.1806E-7
            alpha(1) = 8.574668996e-1_dp; coef(1) = 6.729778096e-3_dp
            alpha(2) = 3.497184772e-1_dp; coef(2) = 5.874145170e-2_dp
            alpha(3) = 1.727917060e-1_dp; coef(3) = 2.339955227e-1_dp
            alpha(4) = 9.373643151e-2_dp; coef(4) = 4.512983737e-1_dp
            alpha(5) = 5.340032759e-2_dp; coef(5) = 3.552053926e-1_dp
            alpha(6) = 3.057364464e-2_dp; coef(6) = 6.974153145e-2_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (16) ! 6s
         SELECT CASE (n)
         CASE (1) ! 0.10174284692137370
            alpha(1) = 0.016053956888690133_dp; coef(1) = 1.000000000000000_dp
         CASE (2) ! 4.7289518343962555E-004
            alpha(1) = 0.026188993287017705_dp; coef(1) = 1.6061745480037777_dp
            alpha(2) = 0.06916086082602657_dp; coef(2) = -0.8493191121412285_dp
         CASE (3) ! 2.2051964796210299E-005
            alpha(1) = 0.025012605890022674_dp; coef(1) = 1.4570914918085252_dp
            alpha(2) = 0.08530028544682955_dp; coef(2) = -0.7792202598134607_dp
            alpha(3) = 0.23664111527549864_dp; coef(3) = 0.08944830584494995_dp
         CASE (4) ! 7.8703519701406321E-008
            alpha(1) = 0.018111879224996335_dp; coef(1) = 0.2742426110211285_dp
            alpha(2) = 0.028489853269487397_dp; coef(2) = 1.3154513257751648_dp
            alpha(3) = 0.07575858742450885_dp; coef(3) = -0.8700267823821266_dp
            alpha(4) = 0.3135211502054625_dp; coef(4) = 0.04933652686369165_dp
         CASE (5) ! 1.4679699233965948E-008
            alpha(1) = 0.017704122931955696_dp; coef(1) = 0.23979549908360312_dp
            alpha(2) = 0.02816438780280081_dp; coef(2) = 1.339995854967376_dp
            alpha(3) = 0.07635157423691516_dp; coef(3) = -0.8631396815113626_dp
            alpha(4) = 0.29378907244008734_dp; coef(4) = 0.05070584407672884_dp
            alpha(5) = 0.7609367674925651_dp; coef(5) = 0.0017452546896168622_dp
         CASE (6) ! 1.1007299516307967E-009
            alpha(1) = 0.018375086251195376_dp; coef(1) = 0.3019609612694464_dp
            alpha(2) = 0.02883213591739931_dp; coef(2) = 1.3047682665168656_dp
            alpha(3) = 0.0735171952084867_dp; coef(3) = -0.8423706115013339_dp
            alpha(4) = 0.13718508837897095_dp; coef(4) = -0.08693662797569257_dp
            alpha(5) = 0.21928060661775448_dp; coef(5) = 0.08040482815742164_dp
            alpha(6) = 0.47255402578342814_dp; coef(6) = 0.011258177324556012_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (17) ! 6p
         SELECT CASE (n)
         CASE (1) ! 1.8826022835702449E-002
            alpha(1) = 0.02720326117498455_dp; coef(1) = 1.0000000000000000_dp
         CASE (2) ! 2.3458059254877206E-004
            alpha(1) = 0.03068554921378584_dp; coef(1) = 1.0932524470142648_dp
            alpha(2) = 0.13327144868476606_dp; coef(2) = -0.20022333756955418_dp
         CASE (3) ! 2.0181514870376560E-006
            alpha(1) = 0.021306449714000038_dp; coef(1) = 0.2589896573574047_dp
            alpha(2) = 0.03558737683798757_dp; coef(2) = 0.8948520093263641_dp
            alpha(3) = 0.1159879775604416_dp; coef(3) = -0.2541710812641519_dp
         CASE (4) ! 5.8503101629803211E-007
            alpha(1) = 0.02094853126286551_dp; coef(1) = 0.23874968368587424_dp
            alpha(2) = 0.03524009985593678_dp; coef(2) = 0.9113069782415131_dp
            alpha(3) = 0.1172381116493136_dp; coef(3) = -0.25095169849099785_dp
            alpha(4) = 1.074757975437711_dp; coef(4) = 0.0014270633206840072_dp
         CASE (5) ! 5.5095088491441402E-008
            alpha(1) = 0.02191796448535917_dp; coef(1) = 0.29958347534070895_dp
            alpha(2) = 0.036406617593622265_dp; coef(2) = 0.8669095451744081_dp
            alpha(3) = 0.10823682952141411_dp; coef(3) = -0.24403207325952317_dp
            alpha(4) = 0.1912844969312091_dp; coef(4) = -0.025610699331187835_dp
            alpha(5) = 0.605311988824089_dp; coef(5) = 0.003533893449371835_dp
         CASE (6) ! 3.1217388674065205E-009
            alpha(1) = 0.01836987035282947_dp; coef(1) = 0.0881865874022419_dp
            alpha(2) = 0.028283616813330213_dp; coef(2) = 0.58260237820929_dp
            alpha(3) = 0.0424808315746305_dp; coef(3) = 0.5465800645782257_dp
            alpha(4) = 0.08775028399997764_dp; coef(4) = -0.20802232268176662_dp
            alpha(5) = 0.14284922158351138_dp; coef(5) = -0.1117378797434228_dp
            alpha(6) = 0.6631101746841396_dp; coef(6) = 0.0028390610573703423_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (18) ! 6d
         SELECT CASE (n)
         CASE (1) ! 2.8755369986284141E-003
            alpha(1) = 0.038539962557578994_dp; coef(1) = 1.0000000000000000_dp
         CASE (2) ! 1.7076604937522610E-003
            alpha(1) = 0.018343252403866257_dp; coef(1) = 0.1002697090882244_dp
            alpha(2) = 0.04095508659742878_dp; coef(2) = 0.921707450678634_dp
         CASE (3) ! 3.0182404084921899E-006
            alpha(1) = 0.026815363906634464_dp; coef(1) = 0.36430995341001843_dp
            alpha(2) = 0.048285664795074575_dp; coef(2) = 0.6902756637148024_dp
            alpha(3) = 0.20398527022566643_dp; coef(3) = -0.061643555476797514_dp
         CASE (4) ! 1.8149850309345794E-006
            alpha(1) = 0.014812592588885804_dp; coef(1) = 0.011887725165713237_dp
            alpha(2) = 0.028604392422323988_dp; coef(2) = 0.4126202273551119_dp
            alpha(3) = 0.049571759954916264_dp; coef(3) = 0.6324841000859261_dp
            alpha(4) = 0.20042444873055187_dp; coef(4) = -0.06355840930231531_dp
         CASE (5) ! 3.5960594502171261E-009
            alpha(1) = 0.02104884011402785_dp; coef(1) = 0.09247565169723151_dp
            alpha(2) = 0.03387333992146401_dp; coef(2) = 0.5183469446932626_dp
            alpha(3) = 0.05464871753445727_dp; coef(3) = 0.4583625355467117_dp
            alpha(4) = 0.16559124243612725_dp; coef(4) = -0.06557358037625106_dp
            alpha(5) = 0.3388543915705088_dp; coef(5) = -0.01073326765431463_dp
         CASE (6) ! 2.2458352955112559E-009
            alpha(1) = 0.012538741501253952_dp; coef(1) = 0.0013891559057176176_dp
            alpha(2) = 0.022411877391156923_dp; coef(2) = 0.12323543239586941_dp
            alpha(3) = 0.03502206086139318_dp; coef(3) = 0.516335137452715_dp
            alpha(4) = 0.05542414129417389_dp; coef(4) = 0.42911284700977503_dp
            alpha(5) = 0.16410056143462523_dp; coef(5) = -0.065949692922377_dp
            alpha(6) = 0.33496788346358164_dp; coef(6) = -0.01123082945384288_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (19) ! 6f
         SELECT CASE (n)
         CASE (1) ! 1.1589084326864452E-002
            alpha(1) = 0.04999708511585302_dp; coef(1) = 1._dp
         CASE (2) ! 1.0223625463345876E-004
            alpha(1) = 0.03443190219402774_dp; coef(1) = 0.5095525299458369_dp
            alpha(2) = 0.06993938754199006_dp; coef(2) = 0.556840322894986_dp
         CASE (3) ! 4.5989236926224919E-005
            alpha(1) = 0.022913376030207443_dp; coef(1) = 0.0844980697596352_dp
            alpha(2) = 0.040621057779302076_dp; coef(2) = 0.558907102542338_dp
            alpha(3) = 0.07551452585114922_dp; coef(3) = 0.4268730468790412_dp
         CASE (4) ! 4.5763401113063651E-007
            alpha(1) = 0.026723978066792376_dp; coef(1) = 0.18313299355518015_dp
            alpha(2) = 0.04675842108187079_dp; coef(2) = 0.5931630352672648_dp
            alpha(3) = 0.08424111895397314_dp; coef(3) = 0.29815135014219757_dp
            alpha(4) = 0.40182223710248033_dp; coef(4) = -0.009003129808147884_dp
         CASE (5) ! 1.5904061800853242E-007
            alpha(1) = 0.019820990273457767_dp; coef(1) = 0.026554164306668928_dp
            alpha(2) = 0.031943455135755784_dp; coef(2) = 0.2857152467032053_dp
            alpha(3) = 0.05180437300600356_dp; coef(3) = 0.5249457199886216_dp
            alpha(4) = 0.08839526561404966_dp; coef(4) = 0.23802584890828332_dp
            alpha(5) = 0.3861344493937384_dp; coef(5) = -0.00959390690588398_dp
         CASE (6) ! 2.7517637057883348E-009
            alpha(1) = 0.021979238503298942_dp; coef(1) = 0.05529549645466706_dp
            alpha(2) = 0.0353635379450325_dp; coef(2) = 0.36260443662220726_dp
            alpha(3) = 0.056802330945072296_dp; coef(3) = 0.47981873895702293_dp
            alpha(4) = 0.09436112630657385_dp; coef(4) = 0.17913120326169676_dp
            alpha(5) = 0.3258167395765562_dp; coef(5) = -0.010025342773471209_dp
            alpha(6) = 0.7758137452357111_dp; coef(6) = -0.0013665302206626205_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (20) ! 6g
         SELECT CASE (n)
         CASE (1) ! 3.0585704257612756E-002
            alpha(1) = 0.061536068942275354_dp; coef(1) = 1._dp
         CASE (2) ! 8.8136922735637668E-004
            alpha(1) = 0.04332629465536664_dp; coef(1) = 0.6209033605306379_dp
            alpha(2) = 0.10013340560986961_dp; coef(2) = 0.4864065117811783_dp
         CASE (3) ! 2.7077377396933500E-005
            alpha(1) = 0.03454561337035713_dp; coef(1) = 0.3227633523004909_dp
            alpha(2) = 0.06729191664458245_dp; coef(2) = 0.6193784317147506_dp
            alpha(3) = 0.14454284915249543_dp; coef(3) = 0.1810539000206162_dp
         CASE (4) ! 7.7121157393555961E-007
            alpha(1) = 0.028983896307773004_dp; coef(1) = 0.15152805717139864_dp
            alpha(2) = 0.05135519389407211_dp; coef(2) = 0.5309701439600488_dp
            alpha(3) = 0.09369371519055944_dp; coef(3) = 0.37457949569535787_dp
            alpha(4) = 0.18939419012133174_dp; coef(4) = 0.06851383804287256_dp
         CASE (5) ! 4.5646906920993047E-008
            alpha(1) = 0.024699673908324754_dp; coef(1) = 0.06048210643335771_dp
            alpha(2) = 0.04083036114944879_dp; coef(2) = 0.3536959373883507_dp
            alpha(3) = 0.06763537320097322_dp; coef(3) = 0.46434612683985144_dp
            alpha(4) = 0.1170091091566765_dp; coef(4) = 0.2130170165201656_dp
            alpha(5) = 0.22066124899602302_dp; coef(5) = 0.0343912188428939_dp
         CASE (6) ! 1.8087060515981079E-008
            alpha(1) = 0.024944035346224053_dp; coef(1) = 0.06443914638193514_dp
            alpha(2) = 0.04138718127148023_dp; coef(2) = 0.36507904601268215_dp
            alpha(3) = 0.06886079537637349_dp; coef(3) = 0.4624056701486661_dp
            alpha(4) = 0.11972435105059721_dp; coef(4) = 0.20306165490832248_dp
            alpha(5) = 0.22639206282070937_dp; coef(5) = 0.03099725943908058_dp
            alpha(6) = 1.4892687021530484_dp; coef(6) = -0.00018343848780163463_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (21) ! 6h
         SELECT CASE (n)
         CASE (1) ! 5.3899613359310816E-002
            alpha(1) = 0.07313314184000295_dp; coef(1) = 1._dp
         CASE (2) ! 2.7493626733567922E-003
            alpha(1) = 0.05113303517913134_dp; coef(1) = 0.6498709211538302_dp
            alpha(2) = 0.12607272311510218_dp; coef(2) = 0.4911257514330729_dp
         CASE (3) ! 1.6537265492422293E-004
            alpha(1) = 0.040851454032819305_dp; coef(1) = 0.3746467681927208_dp
            alpha(2) = 0.08415648437016819_dp; coef(2) = 0.6132065860909666_dp
            alpha(3) = 0.1974425255793456_dp; coef(3) = 0.18179709985670622_dp
         CASE (4) ! 1.1617584376732992E-005
            alpha(1) = 0.034679493065603105_dp; coef(1) = 0.20806488992498948_dp
            alpha(2) = 0.06498507795614539_dp; coef(2) = 0.5659152898854298_dp
            alpha(3) = 0.1276799968878006_dp; coef(3) = 0.3407357791588573_dp
            alpha(4) = 0.293954336816017_dp; coef(4) = 0.061163634708810435_dp
         CASE (5) ! 9.3645049425994387E-007
            alpha(1) = 0.030487493495559674_dp; coef(1) = 0.1147257294844731_dp
            alpha(2) = 0.05371457065868773_dp; coef(2) = 0.45668261354533374_dp
            alpha(3) = 0.09619925104884211_dp; coef(3) = 0.434116528713387_dp
            alpha(4) = 0.18557752006634115_dp; coef(4) = 0.15177021621261605_dp
            alpha(5) = 0.42330005591830827_dp; coef(5) = 0.02006866296723323_dp
         CASE (6) ! 8.5105151370612475E-008
            alpha(1) = 0.027419761875036298_dp; coef(1) = 0.06355638710922437_dp
            alpha(2) = 0.04620886209790054_dp; coef(2) = 0.34336050905161936_dp
            alpha(3) = 0.07793999480997445_dp; coef(3) = 0.45863493125119364_dp
            alpha(4) = 0.13712767813357474_dp; coef(4) = 0.24456067327385683_dp
            alpha(5) = 0.2621356941062822_dp; coef(5) = 0.061018224093793164_dp
            alpha(6) = 0.5948121608710077_dp; coef(6) = 0.006618011961621055_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (22) ! 7s
         SELECT CASE (n)
         CASE (1) ! 0.13154539754786043
            alpha(1) = 0.012170316100569041_dp; coef(1) = 1._dp
         CASE (2) ! 1.2167263197611078E-003
            alpha(1) = 0.02380792064880498_dp; coef(1) = 2.4733498823505897_dp
            alpha(2) = 0.0413963085083486_dp; coef(2) = -1.7500626134752972_dp
         CASE (3) ! 9.2369221911248189E-006
            alpha(1) = 0.020979847986630586_dp; coef(1) = 1.675236183601737_dp
            alpha(2) = 0.061178936229820555_dp; coef(2) = -1.147227387784115_dp
            alpha(3) = 0.1401881324437541_dp; coef(3) = 0.2121061484210291_dp
         CASE (4) ! 2.7376448663574138E-008
            alpha(1) = 0.015915368005026652_dp; coef(1) = 0.3275001969711621_dp
            alpha(2) = 0.023863997843078046_dp; coef(2) = 1.5476507021661345_dp
            alpha(3) = 0.054562166429946386_dp; coef(3) = -1.2697132490611722_dp
            alpha(4) = 0.16342646591350343_dp; coef(4) = 0.13574231777373438_dp
         CASE (5) ! 1.0587942078643664E-008
            alpha(1) = 0.01606210559622015_dp; coef(1) = 0.3489870695674181_dp
            alpha(2) = 0.024001717565201447_dp; coef(2) = 1.5335261688987925_dp
            alpha(3) = 0.05439364645164688_dp; coef(3) = -1.275299080068857_dp
            alpha(4) = 0.16456835599844608_dp; coef(4) = 0.13406783441420897_dp
            alpha(5) = 1.8632404725730192_dp; coef(5) = -0.0001818568135349486_dp
         CASE (6) ! 4.4852654923488444E-010
            alpha(1) = 0.01565231604983483_dp; coef(1) = 0.29154442199853614_dp
            alpha(2) = 0.023623292255693457_dp; coef(2) = 1.5691344169994597_dp
            alpha(3) = 0.05497538487325693_dp; coef(3) = -1.2618704056879249_dp
            alpha(4) = 0.15604786033627724_dp; coef(4) = 0.13683678315315853_dp
            alpha(5) = 0.30526699650551353_dp; coef(5) = 0.006442390493458186_dp
            alpha(6) = 0.7613613898332257_dp; coef(6) = -0.0009736214432373223_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (23) ! 7p
         SELECT CASE (n)
         CASE (1) ! 3.3186835107151591E-002
            alpha(1) = 0.02056749089204262_dp; coef(1) = 1.0000000000000000_dp
         CASE (2) ! 1.7064265376842691E-004
            alpha(1) = 0.02509434622887706_dp; coef(1) = 1.1953936450048124_dp
            alpha(2) = 0.07932687651759715_dp; coef(2) = -0.33554981426999153_dp
         CASE (3) ! 2.6561728884821889E-005
            alpha(1) = 0.02216294295259616_dp; coef(1) = 0.8184285461616602_dp
            alpha(2) = 0.042045475539838104_dp; coef(2) = 0.8328277926936017_dp
            alpha(3) = 0.0622825034958455_dp; coef(3) = -0.7865628896072743_dp
         CASE (4) ! 2.4336365811450378E-008
            alpha(1) = 0.01798640692496774_dp; coef(1) = 0.2517796607641491_dp
            alpha(2) = 0.02838134507613473_dp; coef(2) = 1.0084655616939668_dp
            alpha(3) = 0.07596133909152862_dp; coef(3) = -0.40416345513837804_dp
            alpha(4) = 0.30794728288131357_dp; coef(4) = 0.011606225585276365_dp
         CASE (5) ! 1.6583495687783056E-008
            alpha(1) = 0.018253221396022562_dp; coef(1) = 0.27532738928860534_dp
            alpha(2) = 0.028708043693626412_dp; coef(2) = 0.9951443329321299_dp
            alpha(3) = 0.07061876309146031_dp; coef(3) = -0.2663034712401708_dp
            alpha(4) = 0.08392634608758946_dp; coef(4) = -0.1489548148570699_dp
            alpha(5) = 0.298415801799801_dp; coef(5) = 0.012451440490101195_dp
         CASE (6) ! 7.4004156136453503E-009
            alpha(1) = 0.182855762756E-01_dp; coef(1) = 0.278085459258E+00_dp
            alpha(2) = 0.287346322524E-01_dp; coef(2) = 0.992180661793E+00_dp
            alpha(3) = 0.738243712945E-01_dp; coef(3) = -0.397256781529E+00_dp
            alpha(4) = 0.127917472538E+00_dp; coef(4) = -0.238118471506E-01_dp
            alpha(5) = 0.240852415504E+00_dp; coef(5) = 0.170287124112E-01_dp
            alpha(6) = 0.507988199812E+00_dp; coef(6) = 0.146798451062E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (24) ! 7d
         SELECT CASE (n)
         CASE (1) ! 5.0327602037498753E-003
            alpha(1) = 0.02909039015051073_dp; coef(1) = 1.0000000000000002_dp
         CASE (2) ! 6.3288899527891296E-004
            alpha(1) = 0.03018482426518148_dp; coef(1) = 1.025499895914138_dp
            alpha(2) = 0.15242731708213653_dp; coef(2) = -0.0801196219224754_dp
         CASE (3) ! 5.3862101490054215E-007
            alpha(1) = 0.021405281219812042_dp; coef(1) = 0.2999132044312777_dp
            alpha(2) = 0.03568228986173447_dp; coef(2) = 0.7809338560049855_dp
            alpha(3) = 0.11574152434647746_dp; coef(3) = -0.1242352213758718_dp
         CASE (4) ! 3.8222040621427311E-007
            alpha(1) = 0.021660976213634955_dp; coef(1) = 0.3177854736353611_dp
            alpha(2) = 0.03606163962267541_dp; coef(2) = 0.7675316218690252_dp
            alpha(3) = 0.10985898319545016_dp; coef(3) = -0.11712958841599931_dp
            alpha(4) = 0.15908195611925_dp; coef(4) = -0.011747385912275846_dp
         CASE (5) ! 2.5139247838978918E-007
            alpha(1) = 0.216681757691E-01_dp; coef(1) = 0.318168815768E+00_dp
            alpha(2) = 0.360617800244E-01_dp; coef(2) = 0.766844168160E+00_dp
            alpha(3) = 0.111939820291E+00_dp; coef(3) = -0.125488126534E+00_dp
            alpha(4) = 0.276558830755E+00_dp; coef(4) = -0.448229683916E-02_dp
            alpha(5) = 0.502378300635E+00_dp; coef(5) = 0.158547872047E-02_dp
         CASE (6) ! 3.6283775983747546E-010
            alpha(1) = 0.017457157083607015_dp; coef(1) = 0.07250557935571439_dp
            alpha(2) = 0.026593769405605233_dp; coef(2) = 0.5007763109536889_dp
            alpha(3) = 0.04008932402169768_dp; coef(3) = 0.5339957399786912_dp
            alpha(4) = 0.09560872001992007_dp; coef(4) = -0.11880032574322265_dp
            alpha(5) = 0.15425745053701811_dp; coef(5) = -0.03241452360588799_dp
            alpha(6) = 0.5964979973145561_dp; coef(6) = 0.0007012418324482315_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (25) ! 7f
         SELECT CASE (n)
         CASE (1) ! 4.4122414620793027E-003
            alpha(1) = 0.037696619888011394_dp; coef(1) = 1._dp
         CASE (2) ! 3.8052356073725591E-004
            alpha(1) = 0.02399618488506292_dp; coef(1) = 0.2938876408614713_dp
            alpha(2) = 0.0441727483070652_dp; coef(2) = 0.7461850531456674_dp
         CASE (3) ! 8.2161018017945509E-006
            alpha(1) = 0.02649072884749045_dp; coef(1) = 0.41965629582857367_dp
            alpha(2) = 0.04781482493797306_dp; coef(2) = 0.6320266162094202_dp
            alpha(3) = 0.20807007812612838_dp; coef(3) = -0.025856724476949224_dp
         CASE (4) ! 4.3613902489525458E-007
            alpha(1) = 0.019516667217357033_dp; coef(1) = 0.07760665191295739_dp
            alpha(2) = 0.031536542788218214_dp; coef(2) = 0.5138599419165162_dp
            alpha(3) = 0.0518809354792556_dp; coef(3) = 0.4658306981238731_dp
            alpha(4) = 0.18978036627517764_dp; coef(4) = -0.03032718718680631_dp
         CASE (5) ! 1.3995496761509685E-008
            alpha(1) = 0.02085272354115486_dp; coef(1) = 0.11822768394863675_dp
            alpha(2) = 0.033574888134901455_dp; coef(2) = 0.5490677063820284_dp
            alpha(3) = 0.05432263605261335_dp; coef(3) = 0.3941968220546063_dp
            alpha(4) = 0.1670441095948162_dp; coef(4) = -0.031436205691901344_dp
            alpha(5) = 0.3442595341279306_dp; coef(5) = -0.0033909812829636394_dp
         CASE (6) ! 8.1826723175026927E-010
            alpha(1) = 0.01641594422436568_dp; coef(1) = 0.01709822656839514_dp
            alpha(2) = 0.024762426761735977_dp; coef(2) = 0.22140165119709196_dp
            alpha(3) = 0.03711681387997685_dp; coef(3) = 0.5157825082425954_dp
            alpha(4) = 0.057006537132452996_dp; coef(4) = 0.30947877900909615_dp
            alpha(5) = 0.15940261392199198_dp; coef(5) = -0.032377223278152595_dp
            alpha(6) = 0.3168807735151213_dp; coef(6) = -0.004680397785899489_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (26) ! 7g
         SELECT CASE (n)
         CASE (1) ! 1.6232140186735133E-002
            alpha(1) = 0.04636088534182354_dp; coef(1) = 1.0000000000000002_dp
         CASE (2) ! 1.8148563973063858E-004
            alpha(1) = 0.033089231025105_dp; coef(1) = 0.5580805992890315_dp
            alpha(2) = 0.06689561975474233_dp; coef(2) = 0.5213501458304487_dp
         CASE (3) ! 7.1118265680958359E-006
            alpha(1) = 0.02536613332751714_dp; coef(1) = 0.1958347745511226_dp
            alpha(2) = 0.04392292490201606_dp; coef(2) = 0.6069779654592313_dp
            alpha(3) = 0.07935572931927462_dp; coef(3) = 0.28358299617643856_dp
         CASE (4) ! 1.3483475262976441E-006
            alpha(1) = 0.026003867919966516_dp; coef(1) = 0.2209616450163026_dp
            alpha(2) = 0.045349648935287534_dp; coef(2) = 0.6131566254139853_dp
            alpha(3) = 0.08207336886154683_dp; coef(3) = 0.25321528599323867_dp
            alpha(4) = 0.43630841028097767_dp; coef(4) = -0.0028051508587602996_dp
         CASE (5) ! 2.3403034621338747E-008
            alpha(1) = 0.021295076240560495_dp; coef(1) = 0.06733674728927509_dp
            alpha(2) = 0.03407297844465803_dp; coef(2) = 0.39912419436935137_dp
            alpha(3) = 0.05463343066262556_dp; coef(3) = 0.47165433606934565_dp
            alpha(4) = 0.09142718539759649_dp; coef(4) = 0.15054429419419305_dp
            alpha(5) = 0.3631735814032239_dp; coef(5) = -0.003735801628196915_dp
         CASE (6) ! 8.6698974666177264E-009
            alpha(1) = 0.02154419398101922_dp; coef(1) = 0.07294053032928545_dp
            alpha(2) = 0.03458506500828156_dp; coef(2) = 0.41255126079553106_dp
            alpha(3) = 0.055591221080949486_dp; coef(3) = 0.4632159545572935_dp
            alpha(4) = 0.0929128679746364_dp; coef(4) = 0.14023937661406255_dp
            alpha(5) = 0.3353651542417124_dp; coef(5) = -0.00381264707490568_dp
            alpha(6) = 0.8332537115987487_dp; coef(6) = -0.0002777669788215993_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (27) ! 7h
         SELECT CASE (n)
         CASE (1) ! 3.3956902793941257E-002
            alpha(1) = 0.055067144574685124_dp; coef(1) = 1.0000000000000002_dp
         CASE (2) ! 1.0620083024626137E-003
            alpha(1) = 0.039504306298949056_dp; coef(1) = 0.6244224977499361_dp
            alpha(2) = 0.08748752023471955_dp; coef(2) = 0.48947010311740946_dp
         CASE (3) ! 3.6773022004554790E-005
            alpha(1) = 0.031898703524951964_dp; coef(1) = 0.329315310335007_dp
            alpha(2) = 0.060322239914691864_dp; coef(2) = 0.6230689660925297_dp
            alpha(3) = 0.12542722211419774_dp; coef(3) = 0.1792628471936351_dp
         CASE (4) ! 1.3161279064366482E-006
            alpha(1) = 0.027133297390849163_dp; coef(1) = 0.16105237279351875_dp
            alpha(2) = 0.04709911293196238_dp; coef(2) = 0.5444503553138231_dp
            alpha(3) = 0.08417053216258696_dp; coef(3) = 0.36701928787783233_dp
            alpha(4) = 0.16769024925679932_dp; coef(4) = 0.06200372233567246_dp
         CASE (5) ! 4.2949577450102083E-008
            alpha(1) = 0.0236897013320386_dp; coef(1) = 0.07357940603188007_dp
            alpha(2) = 0.0387782492773279_dp; coef(2) = 0.39203903619828495_dp
            alpha(3) = 0.0637197820910416_dp; coef(3) = 0.46247919691910666_dp
            alpha(4) = 0.10978567271801772_dp; coef(4) = 0.1832898128285397_dp
            alpha(5) = 0.2085364358880184_dp; coef(5) = 0.023617413014962847_dp
         CASE (6) ! 2.7450561823627595E-009
            alpha(1) = 0.020820906789523693_dp; coef(1) = 0.029194538236458435_dp
            alpha(2) = 0.032494198134763036_dp; coef(2) = 0.23574088940165333_dp
            alpha(3) = 0.05023584266268183_dp; coef(3) = 0.44360291827793874_dp
            alpha(4) = 0.0794251556241889_dp; coef(4) = 0.3155567845762703_dp
            alpha(5) = 0.13148712880820534_dp; coef(5) = 0.09872391034903429_dp
            alpha(6) = 0.2359872300875659_dp; coef(6) = 0.012260115931631366_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (28) ! 7i
         SELECT CASE (n)
         CASE (1) ! 5.4422862380136161E-002
            alpha(1) = 0.06380473723420682_dp; coef(1) = 0.9999999999999998_dp
         CASE (2) ! 2.6995295321672330E-003
            alpha(1) = 0.045457834184705205_dp; coef(1) = 0.6463784604860205_dp
            alpha(2) = 0.10559967649617007_dp; coef(2) = 0.49629977931368163_dp
         CASE (3) ! 1.5458937276235929E-004
            alpha(1) = 0.03666777792117483_dp; coef(1) = 0.3678119967003228_dp
            alpha(2) = 0.07224640498805254_dp; coef(2) = 0.618554233547922_dp
            alpha(3) = 0.15965614371413928_dp; coef(3) = 0.18501800416077063_dp
         CASE (4) ! 1.0208819409163183E-005
            alpha(1) = 0.031314605885208034_dp; coef(1) = 0.20067272567322478_dp
            alpha(2) = 0.05654362551330421_dp; coef(2) = 0.5656501384781_dp
            alpha(3) = 0.1062477818701912_dp; coef(3) = 0.3489852946890129_dp
            alpha(4) = 0.230241610829344_dp; coef(4) = 0.06217858499289526_dp
         CASE (5) ! 7.6749801847242338E-007
            alpha(1) = 0.02764326937219174_dp; coef(1) = 0.1083794411727515_dp
            alpha(2) = 0.047149043197080166_dp; coef(2) = 0.4503068342384467_dp
            alpha(3) = 0.08137277011103897_dp; coef(3) = 0.4431260931378365_dp
            alpha(4) = 0.1500423393404324_dp; coef(4) = 0.15687089041229385_dp
            alpha(5) = 0.32194304635741355_dp; coef(5) = 0.020215660073712404_dp
         CASE (6) ! 6.4734517968290106E-008
            alpha(1) = 0.024937530903194137_dp; coef(1) = 0.058696755716280226_dp
            alpha(2) = 0.04081596003413946_dp; coef(2) = 0.3330255270864295_dp
            alpha(3) = 0.06665805855772095_dp; coef(3) = 0.4640568401312629_dp
            alpha(4) = 0.11296341467059927_dp; coef(4) = 0.253753255644_dp
            alpha(5) = 0.2062966749989214_dp; coef(5) = 0.06316006418978992_dp
            alpha(6) = 0.4401320366389093_dp; coef(6) = 0.0065607991581329616_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (29) ! 8s
         SELECT CASE (n)
         CASE (1) !Final Error     0.16012012792726149              Nfun:          26
            alpha(1) = 0.954817062770E-02_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) !Final Error      2.5381935191433058E-003         Nfun:        5000
            alpha(1) = 0.233437156506E-01_dp; coef(1) = 0.119129429454E+02_dp
            alpha(2) = 0.258428856774E-01_dp; coef(2) = -0.112184681311E+02_dp
         CASE (3) !Final Error      3.7901704366127831E-006         Nfun:          84
            alpha(1) = 0.180016350547E-01_dp; coef(1) = 0.192251210431E+01_dp
            alpha(2) = 0.483977532005E-01_dp; coef(2) = -0.174735592554E+01_dp
            alpha(3) = 0.862089558510E-01_dp; coef(3) = 0.542216684260E+00_dp
         CASE (4) !Final Error      8.7725764097434933E-007         Nfun:         124
            alpha(1) = 0.112292312447E-01_dp; coef(1) = 0.800230664629E-01_dp
            alpha(2) = 0.189551778937E-01_dp; coef(2) = 0.198501632216E+01_dp
            alpha(3) = 0.447185315682E-01_dp; coef(3) = -0.172756054240E+01_dp
            alpha(4) = 0.936481128491E-01_dp; coef(4) = 0.380520807030E+00_dp
         CASE (5) !Final Error      6.5327121490099671E-010         Nfun:         156
            alpha(1) = 0.140923696754E-01_dp; coef(1) = 0.372819666670E+00_dp
            alpha(2) = 0.204188098143E-01_dp; coef(2) = 0.184957327860E+01_dp
            alpha(3) = 0.421384023783E-01_dp; coef(3) = -0.179644496345E+01_dp
            alpha(4) = 0.103255591497E+00_dp; coef(4) = 0.295269129687E+00_dp
            alpha(5) = 0.407317552947E+00_dp; coef(5) = -0.363956012755E-02_dp
         CASE (6) !Final Error      3.5843905621391059E-010         Nfun:         202
            alpha(1) = 0.774373712240E-02_dp; coef(1) = 0.177912928441E-02_dp
            alpha(2) = 0.146710572275E-01_dp; coef(2) = 0.473873521097E+00_dp
            alpha(3) = 0.207400451684E-01_dp; coef(3) = 0.176348848069E+01_dp
            alpha(4) = 0.420012824724E-01_dp; coef(4) = -0.181118198359E+01_dp
            alpha(5) = 0.103495646254E+00_dp; coef(5) = 0.293340582373E+00_dp
            alpha(6) = 0.403606855142E+00_dp; coef(6) = -0.371676924018E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (30) ! 8p
         SELECT CASE (n)
         CASE (1) ! Final Error      4.9424644896695069E-002         Nfun:          32
            alpha(1) = 0.161028140061E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      2.8975996622193634E-004         Nfun:          56
            alpha(1) = 0.214927884597E-01_dp; coef(1) = 0.137340335964E+01_dp
            alpha(2) = 0.517288456046E-01_dp; coef(2) = -0.543362689377E+00_dp
         CASE (3) ! Final Error      2.4112893364236498E-005         Nfun:        5000
            alpha(1) = 0.207032799112E-01_dp; coef(1) = 0.125897435808E+01_dp
            alpha(2) = 0.775417675460E-01_dp; coef(2) = -0.175313428477E+01_dp
            alpha(3) = 0.853618175731E-01_dp; coef(3) = 0.133222160068E+01_dp
         CASE (4) ! Final Error      5.1550181900950065E-009         Nfun:         115
            alpha(1) = 0.157817873157E-01_dp; coef(1) = 0.275272236164E+00_dp
            alpha(2) = 0.237507657245E-01_dp; coef(2) = 0.113175870878E+01_dp
            alpha(3) = 0.546913738102E-01_dp; coef(3) = -0.604838325903E+00_dp
            alpha(4) = 0.162719740193E+00_dp; coef(4) = 0.379559818155E-01_dp
         CASE (5) ! Final Error      3.3981357727697059E-009         Nfun:         208
            alpha(1) = 0.831136251369E-02_dp; coef(1) = 0.120783003696E-02_dp
            alpha(2) = 0.162032540236E-01_dp; coef(2) = 0.311803068732E+00_dp
            alpha(3) = 0.239700730039E-01_dp; coef(3) = 0.109809634731E+01_dp
            alpha(4) = 0.545591646588E-01_dp; coef(4) = -0.608585838680E+00_dp
            alpha(5) = 0.163206247987E+00_dp; coef(5) = 0.376432425906E-01_dp
         CASE (6) ! Final Error      1.7629542270469756E-009         Nfun:         181
            alpha(1) = 0.762583234870E-02_dp; coef(1) = 0.608041238574E-03_dp
            alpha(2) = 0.159886043406E-01_dp; coef(2) = 0.290715964586E+00_dp
            alpha(3) = 0.238233242498E-01_dp; coef(3) = 0.111572208650E+01_dp
            alpha(4) = 0.547475656337E-01_dp; coef(4) = -0.605573011155E+00_dp
            alpha(5) = 0.159838909861E+00_dp; coef(5) = 0.381832212541E-01_dp
            alpha(6) = 0.303650029843E+00_dp; coef(6) = 0.521829479449E-03_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (31) ! 8d
         SELECT CASE (n)
         CASE (1) ! Final Error      1.1007451787931499E-002         Nfun:          30
            alpha(1) = 0.227452786936E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      2.6261780162051274E-004         Nfun:          59
            alpha(1) = 0.246415518724E-01_dp; coef(1) = 0.106888910159E+01_dp
            alpha(2) = 0.882423565832E-01_dp; coef(2) = -0.150832991506E+00_dp
         CASE (3) ! Final Error      3.1438629715108135E-006         Nfun:          79
            alpha(1) = 0.187908325877E-01_dp; coef(1) = 0.345753063788E+00_dp
            alpha(2) = 0.294062808451E-01_dp; coef(2) = 0.805900635094E+00_dp
            alpha(3) = 0.729972483983E-01_dp; coef(3) = -0.226419150580E+00_dp
         CASE (4) ! Final Error      4.6681045251517617E-008         Nfun:         155
            alpha(1) = 0.180591659594E-01_dp; coef(1) = 0.276777914435E+00_dp
            alpha(2) = 0.284516297282E-01_dp; coef(2) = 0.859179877241E+00_dp
            alpha(3) = 0.757835234545E-01_dp; coef(3) = -0.212512286908E+00_dp
            alpha(4) = 0.316179197045E+00_dp; coef(4) = 0.290411896524E-02_dp
         CASE (5) ! Final Error      2.8592438860997049E-008         Nfun:         157
            alpha(1) = 0.788957125965E-02_dp; coef(1) = 0.922033355948E-03_dp
            alpha(2) = 0.184048411402E-01_dp; coef(2) = 0.297899483367E+00_dp
            alpha(3) = 0.286569642469E-01_dp; coef(3) = 0.838829716148E+00_dp
            alpha(4) = 0.755539142959E-01_dp; coef(4) = -0.214060069719E+00_dp
            alpha(5) = 0.321295771453E+00_dp; coef(5) = 0.281601399193E-02_dp
         CASE (6) ! Final Error      2.5824181237865673E-008         Nfun:         360
            alpha(1) = 0.799566301376E-02_dp; coef(1) = 0.989916078242E-03_dp
            alpha(2) = 0.184314393623E-01_dp; coef(2) = 0.299806753776E+00_dp
            alpha(3) = 0.286793174182E-01_dp; coef(3) = 0.837161648054E+00_dp
            alpha(4) = 0.754861483294E-01_dp; coef(4) = -0.214279054211E+00_dp
            alpha(5) = 0.347665717476E+00_dp; coef(5) = 0.316665708546E-02_dp
            alpha(6) = 0.519145041687E+00_dp; coef(6) = -0.456310557032E-03_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (32) ! 8f
         SELECT CASE (n)
         CASE (1) ! Final Error      2.2576226875590244E-003         Nfun:          34
            alpha(1) = 0.294475259826E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      1.5251988552043017E-003         Nfun:          57
            alpha(1) = 0.296782812720E-01_dp; coef(1) = 0.100475463592E+01_dp
            alpha(2) = 0.190795493774E+00_dp; coef(2) = -0.287320744182E-01_dp
         CASE (3) ! Final Error      1.1028715960481605E-006         Nfun:          81
            alpha(1) = 0.213883927368E-01_dp; coef(1) = 0.349292158371E+00_dp
            alpha(2) = 0.356767971245E-01_dp; coef(2) = 0.706980758418E+00_dp
            alpha(3) = 0.115522086790E+00_dp; coef(3) = -0.624625138479E-01_dp
         CASE (4) ! Final Error      2.5905283851379579E-007         Nfun:         114
            alpha(1) = 0.135735156650E-01_dp; coef(1) = 0.160565058328E-01_dp
            alpha(2) = 0.230188101338E-01_dp; coef(2) = 0.412074616105E+00_dp
            alpha(3) = 0.367683495092E-01_dp; coef(3) = 0.631581445296E+00_dp
            alpha(4) = 0.113172928529E+00_dp; coef(4) = -0.654897582475E-01_dp
         CASE (5) ! Final Error      9.1994127870975717E-009         Nfun:         138
            alpha(1) = 0.174175453826E-01_dp; coef(1) = 0.926732388033E-01_dp
            alpha(2) = 0.266164388803E-01_dp; coef(2) = 0.530370002863E+00_dp
            alpha(3) = 0.402882842937E-01_dp; coef(3) = 0.452897093313E+00_dp
            alpha(4) = 0.909981814695E-01_dp; coef(4) = -0.598885399214E-01_dp
            alpha(5) = 0.140524344225E+00_dp; coef(5) = -0.220199882869E-01_dp
         CASE (6) ! Final Error      9.2777363569496174E-010         Nfun:         201
            alpha(1) = 0.171461611966E-01_dp; coef(1) = 0.818763681195E-01_dp
            alpha(2) = 0.261005461673E-01_dp; coef(2) = 0.507613187845E+00_dp
            alpha(3) = 0.395039423214E-01_dp; coef(3) = 0.481514356460E+00_dp
            alpha(4) = 0.991873066686E-01_dp; coef(4) = -0.676983741468E-01_dp
            alpha(5) = 0.169263793429E+00_dp; coef(5) = -0.951176902180E-02_dp
            alpha(6) = 0.446424281621E+00_dp; coef(6) = 0.287939510863E-03_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (33) ! 8g
         SELECT CASE (n)
         CASE (1) ! Final Error      7.6546829464960719E-003         Nfun:          32
            alpha(1) = 0.361922956953E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      9.1507803362755880E-005         Nfun:          56
            alpha(1) = 0.253468957727E-01_dp; coef(1) = 0.432183514630E+00_dp
            alpha(2) = 0.457296941015E-01_dp; coef(2) = 0.622891258660E+00_dp
         CASE (3) ! Final Error      2.5033245749872535E-005         Nfun:          83
            alpha(1) = 0.258871489345E-01_dp; coef(1) = 0.465225762875E+00_dp
            alpha(2) = 0.467908293387E-01_dp; coef(2) = 0.592676637235E+00_dp
            alpha(3) = 0.235898727433E+00_dp; coef(3) = -0.936999722044E-02_dp
         CASE (4) ! Final Error      1.2374131963532875E-007         Nfun:          99
            alpha(1) = 0.202912068011E-01_dp; coef(1) = 0.134834556651E+00_dp
            alpha(2) = 0.325870891557E-01_dp; coef(2) = 0.568204416951E+00_dp
            alpha(3) = 0.529799382791E-01_dp; coef(3) = 0.361768272513E+00_dp
            alpha(4) = 0.182951379245E+00_dp; coef(4) = -0.143138058545E-01_dp
         CASE (5) ! Final Error      6.1520559757610727E-008         Nfun:         137
            alpha(1) = 0.204904098441E-01_dp; coef(1) = 0.143247104144E+00_dp
            alpha(2) = 0.329662129320E-01_dp; coef(2) = 0.574577306053E+00_dp
            alpha(3) = 0.535480229206E-01_dp; coef(3) = 0.347771773873E+00_dp
            alpha(4) = 0.173905195437E+00_dp; coef(4) = -0.146105870157E-01_dp
            alpha(5) = 0.411974802698E+00_dp; coef(5) = -0.630478427757E-03_dp
         CASE (6) ! Final Error      2.8878899271944647E-010         Nfun:         224
            alpha(1) = 0.169892734782E-01_dp; coef(1) = 0.343634181341E-01_dp
            alpha(2) = 0.255306653359E-01_dp; coef(2) = 0.297490048232E+00_dp
            alpha(3) = 0.380226605577E-01_dp; coef(3) = 0.504330081770E+00_dp
            alpha(4) = 0.578237162858E-01_dp; coef(4) = 0.232028681855E+00_dp
            alpha(5) = 0.156319297580E+00_dp; coef(5) = -0.158757742665E-01_dp
            alpha(6) = 0.303215876261E+00_dp; coef(6) = -0.190549114536E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (34) ! 8h
         SELECT CASE (n)
         CASE (1) ! Final Error      2.0263452932776316E-002         Nfun:          30
            alpha(1) = 0.429683678636E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      3.1675309856771250E-004         Nfun:          55
            alpha(1) = 0.313343436418E-01_dp; coef(1) = 0.581445295878E+00_dp
            alpha(2) = 0.624854594695E-01_dp; coef(2) = 0.507843999042E+00_dp
         CASE (3) ! Final Error      3.9852674413864975E-006         Nfun:          84
            alpha(1) = 0.250505193613E-01_dp; coef(1) = 0.249866933077E+00_dp
            alpha(2) = 0.433715388861E-01_dp; coef(2) = 0.624281652229E+00_dp
            alpha(3) = 0.785922725627E-01_dp; coef(3) = 0.224723740642E+00_dp
         CASE (4) ! Final Error      9.6134466009978325E-007         Nfun:         102
            alpha(1) = 0.197053704551E-01_dp; coef(1) = 0.587329016952E-01_dp
            alpha(2) = 0.311703214507E-01_dp; coef(2) = 0.391556705376E+00_dp
            alpha(3) = 0.498252109293E-01_dp; coef(3) = 0.495622127752E+00_dp
            alpha(4) = 0.847225051968E-01_dp; coef(4) = 0.153759464550E+00_dp
         CASE (5) ! Final Error      3.1387072496613655E-008         Nfun:         171
            alpha(1) = 0.209637027387E-01_dp; coef(1) = 0.896162511950E-01_dp
            alpha(2) = 0.335031104007E-01_dp; coef(2) = 0.451988472126E+00_dp
            alpha(3) = 0.537855393028E-01_dp; coef(3) = 0.445684549220E+00_dp
            alpha(4) = 0.904685841856E-01_dp; coef(4) = 0.112975621737E+00_dp
            alpha(5) = 0.367009548127E+00_dp; coef(5) = -0.128702498903E-02_dp
         CASE (6) ! Final Error      2.3379744806817371E-009         Nfun:         178
            alpha(1) = 0.176131589633E-01_dp; coef(1) = 0.228801729726E-01_dp
            alpha(2) = 0.265241759564E-01_dp; coef(2) = 0.219935173899E+00_dp
            alpha(3) = 0.395682502275E-01_dp; coef(3) = 0.457742042442E+00_dp
            alpha(4) = 0.603163601880E-01_dp; coef(4) = 0.322355672010E+00_dp
            alpha(5) = 0.962509769419E-01_dp; coef(5) = 0.776110865189E-01_dp
            alpha(6) = 0.340456560556E+00_dp; coef(6) = -0.149039074155E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (35) ! 8i
         SELECT CASE (n)
         CASE (1) ! Final Error      3.6610861464450251E-002         Nfun:          28
            alpha(1) = 0.497681342869E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      1.2037998634699676E-003         Nfun:          54
            alpha(1) = 0.362582580722E-01_dp; coef(1) = 0.625797636897E+00_dp
            alpha(2) = 0.772789743491E-01_dp; coef(2) = 0.493087055190E+00_dp
         CASE (3) ! Final Error      4.4138907268287397E-005         Nfun:          84
            alpha(1) = 0.295565302253E-01_dp; coef(1) = 0.331648536218E+00_dp
            alpha(2) = 0.543650167678E-01_dp; coef(2) = 0.626252487836E+00_dp
            alpha(3) = 0.109385279185E+00_dp; coef(3) = 0.180132637952E+00_dp
         CASE (4) ! Final Error      1.7565340795222539E-006         Nfun:         109
            alpha(1) = 0.253506952420E-01_dp; coef(1) = 0.164411803656E+00_dp
            alpha(2) = 0.430745189784E-01_dp; coef(2) = 0.550397809114E+00_dp
            alpha(3) = 0.752046515091E-01_dp; coef(3) = 0.366301098603E+00_dp
            alpha(4) = 0.146328151358E+00_dp; coef(4) = 0.601407102534E-01_dp
         CASE (5) ! Final Error      7.0760998394092667E-008         Nfun:         135
            alpha(1) = 0.223631558595E-01_dp; coef(1) = 0.782316855563E-01_dp
            alpha(2) = 0.360324076062E-01_dp; coef(2) = 0.405608377951E+00_dp
            alpha(3) = 0.582625147882E-01_dp; coef(3) = 0.463320173139E+00_dp
            alpha(4) = 0.988379770544E-01_dp; coef(4) = 0.174224918486E+00_dp
            alpha(5) = 0.186254974327E+00_dp; coef(5) = 0.204375851962E-01_dp
         CASE (6) ! Final Error      2.5812996184981785E-009         Nfun:         154
            alpha(1) = 0.200353504400E-01_dp; coef(1) = 0.353572119094E-01_dp
            alpha(2) = 0.310058228263E-01_dp; coef(2) = 0.264986349092E+00_dp
            alpha(3) = 0.475789622751E-01_dp; coef(3) = 0.459604400295E+00_dp
            alpha(4) = 0.747945594816E-01_dp; coef(4) = 0.294520734067E+00_dp
            alpha(5) = 0.123608348534E+00_dp; coef(5) = 0.795678135397E-01_dp
            alpha(6) = 0.223837832433E+00_dp; coef(6) = 0.788911126553E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (36) ! 8k
         SELECT CASE (n)
         CASE (1) ! Final Error      5.4828647243574391E-002         Nfun:          31
            alpha(1) = 0.565862527158E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      2.6585281818616302E-003         Nfun:          58
            alpha(1) = 0.409621413283E-01_dp; coef(1) = 0.643341930303E+00_dp
            alpha(2) = 0.905906355164E-01_dp; coef(2) = 0.500650234401E+00_dp
         CASE (3) ! Final Error      1.4615614325030890E-004         Nfun:          79
            alpha(1) = 0.333174897003E-01_dp; coef(1) = 0.362014818487E+00_dp
            alpha(2) = 0.632575769483E-01_dp; coef(2) = 0.622856136324E+00_dp
            alpha(3) = 0.133092573863E+00_dp; coef(3) = 0.187857098375E+00_dp
         CASE (4) ! Final Error      9.1612612234115431E-006         Nfun:         101
            alpha(1) = 0.286033556743E-01_dp; coef(1) = 0.194503877864E+00_dp
            alpha(2) = 0.500739500873E-01_dp; coef(2) = 0.565008796654E+00_dp
            alpha(3) = 0.906868771708E-01_dp; coef(3) = 0.356072167996E+00_dp
            alpha(4) = 0.187014077875E+00_dp; coef(4) = 0.631401639920E-01_dp
         CASE (5) ! Final Error      6.4892919304782026E-007         Nfun:         129
            alpha(1) = 0.253426041474E-01_dp; coef(1) = 0.103171589032E+00_dp
            alpha(2) = 0.420679943919E-01_dp; coef(2) = 0.444461963838E+00_dp
            alpha(3) = 0.704104195426E-01_dp; coef(3) = 0.450607242799E+00_dp
            alpha(4) = 0.125087213598E+00_dp; coef(4) = 0.161443127250E+00_dp
            alpha(5) = 0.255297703738E+00_dp; coef(5) = 0.203951623802E-01_dp
         CASE (6) ! Final Error      5.1319749783118596E-008         Nfun:         167
            alpha(1) = 0.229249907854E-01_dp; coef(1) = 0.547921797808E-01_dp
            alpha(2) = 0.366159972070E-01_dp; coef(2) = 0.324017832257E+00_dp
            alpha(3) = 0.582200416027E-01_dp; coef(3) = 0.468159342618E+00_dp
            alpha(4) = 0.956611769853E-01_dp; coef(4) = 0.261784418007E+00_dp
            alpha(5) = 0.168269053498E+00_dp; coef(5) = 0.651261040640E-01_dp
            alpha(6) = 0.341392440896E+00_dp; coef(6) = 0.652958748406E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (37) ! 9s
         SELECT CASE (n)
         CASE (1) ! Final Error     0.18725174229991048              Nfun:          27
            alpha(1) = 0.769349917851E-02_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      5.4066663410257121E-003         Nfun:        5000
            alpha(1) = 0.189568723452E-01_dp; coef(1) = 0.112741628118E+02_dp
            alpha(2) = 0.211773359716E-01_dp; coef(2) = -0.106026364073E+02_dp
         CASE (3) ! Final Error      5.7494822049974914E-006         Nfun:        5000
            alpha(1) = 0.156765999362E-01_dp; coef(1) = 0.217812259686E+01_dp
            alpha(2) = 0.452235262760E-01_dp; coef(2) = -0.775211723209E+01_dp
            alpha(3) = 0.501515376761E-01_dp; coef(3) = 0.627207417672E+01_dp
         CASE (4) ! Final Error      5.3790535541331508E-006         Nfun:         138
            alpha(1) = 0.456007702938E-02_dp; coef(1) = 0.273356915234E-02_dp
            alpha(2) = 0.157626429250E-01_dp; coef(2) = 0.219752378246E+01_dp
            alpha(3) = 0.439949979621E-01_dp; coef(3) = -0.553333931604E+01_dp
            alpha(4) = 0.513218824088E-01_dp; coef(4) = 0.403137965730E+01_dp
         CASE (5) ! Final Error      1.1499532903691545E-007         Nfun:        5000
            alpha(1) = 0.145726040912E-01_dp; coef(1) = 0.144719347822E+01_dp
            alpha(2) = 0.254992984481E-01_dp; coef(2) = 0.798252833958E+01_dp
            alpha(3) = 0.281130573683E-01_dp; coef(3) = -0.912426799710E+01_dp
            alpha(4) = 0.864957612715E-01_dp; coef(4) = 0.480287729936E+00_dp
            alpha(5) = 0.139046050110E+00_dp; coef(5) = -0.883281830654E-01_dp
         CASE (6) ! Final Error      5.3905318209146458E-008         Nfun:        5000
            alpha(1) = 0.145726040912E-01_dp; coef(1) = 0.144495498475E+01_dp
            alpha(2) = 0.254992984481E-01_dp; coef(2) = 0.801876742489E+01_dp
            alpha(3) = 0.281130573683E-01_dp; coef(3) = -0.916277091786E+01_dp
            alpha(4) = 0.854957612715E-01_dp; coef(4) = 0.483824740433E+00_dp
            alpha(5) = 0.133046050110E+00_dp; coef(5) = -0.854674447612E-01_dp
            alpha(6) = 0.338046050110E+00_dp; coef(6) = -0.199498664467E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (38) ! 9p
         SELECT CASE (n)
         CASE (1) ! Final Error      6.6567424680416876E-002         Nfun:          27
            alpha(1) = 0.129537478531E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      5.7874505043797875E-004         Nfun:          59
            alpha(1) = 0.192770197819E-01_dp; coef(1) = 0.175728933087E+01_dp
            alpha(2) = 0.352243700375E-01_dp; coef(2) = -0.954147437837E+00_dp
         CASE (3) ! Final Error      1.0318638513195921E-005         Nfun:          81
            alpha(1) = 0.178673196113E-01_dp; coef(1) = 0.139343045871E+01_dp
            alpha(2) = 0.527645536481E-01_dp; coef(2) = -0.991286262027E+00_dp
            alpha(3) = 0.738583901124E-01_dp; coef(3) = 0.413150119926E+00_dp
         CASE (4) ! Final Error      4.8420569775942113E-008         Nfun:         110
            alpha(1) = 0.136234961681E-01_dp; coef(1) = 0.242173247853E+00_dp
            alpha(2) = 0.200005112607E-01_dp; coef(2) = 0.132452980029E+01_dp
            alpha(3) = 0.427996048691E-01_dp; coef(3) = -0.849105094814E+00_dp
            alpha(4) = 0.996245179618E-01_dp; coef(4) = 0.990043491665E-01_dp
         CASE (5) ! Final Error      4.2303451852632179E-009         Nfun:         155
            alpha(1) = 0.140237921640E-01_dp; coef(1) = 0.300921484527E+00_dp
            alpha(2) = 0.203711629454E-01_dp; coef(2) = 0.129026895819E+01_dp
            alpha(3) = 0.421355382404E-01_dp; coef(3) = -0.864335601306E+00_dp
            alpha(4) = 0.103998185340E+00_dp; coef(4) = 0.905676566706E-01_dp
            alpha(5) = 0.296838493801E+00_dp; coef(5) = -0.921081091161E-03_dp
         CASE (6) ! Final Error      7.8384676527321062E-010         Nfun:         209
            alpha(1) = 0.137404800610E-01_dp; coef(1) = 0.256851266967E+00_dp
            alpha(2) = 0.200717709084E-01_dp; coef(2) = 0.131145904095E+01_dp
            alpha(3) = 0.430180557283E-01_dp; coef(3) = -0.867675787268E+00_dp
            alpha(4) = 0.725488993672E-01_dp; coef(4) = 0.503842041627E-01_dp
            alpha(5) = 0.110730465365E+00_dp; coef(5) = 0.665510025893E-01_dp
            alpha(6) = 0.312803980774E+00_dp; coef(6) = -0.107397924374E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (39) ! 9d
         SELECT CASE (n)
         CASE (1) ! Final Error      1.9446423107279065E-002         Nfun:          29
            alpha(1) = 0.182773605707E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      1.7827653387159081E-004         Nfun:          60
            alpha(1) = 0.208881247198E-01_dp; coef(1) = 0.114298696342E+01_dp
            alpha(2) = 0.580996673168E-01_dp; coef(2) = -0.250564466849E+00_dp
         CASE (3) ! Final Error      1.8158482513186414E-005         Nfun:        5000
            alpha(1) = 0.187491333354E-01_dp; coef(1) = 0.799553035255E+00_dp
            alpha(2) = 0.366808493267E-01_dp; coef(2) = 0.203831143189E+01_dp
            alpha(3) = 0.412043605803E-01_dp; coef(3) = -0.194036514672E+01_dp
         CASE (4) ! Final Error      8.8698213218663113E-009         Nfun:         112
            alpha(1) = 0.156726689461E-01_dp; coef(1) = 0.268520959785E+00_dp
            alpha(2) = 0.236479821892E-01_dp; coef(2) = 0.947834833544E+00_dp
            alpha(3) = 0.548782222741E-01_dp; coef(3) = -0.328498127965E+00_dp
            alpha(4) = 0.160594602797E+00_dp; coef(4) = 0.123971623459E-01_dp
         CASE (5) ! Final Error      1.5829011434931317E-009         Nfun:         131
            alpha(1) = 0.888001143555E-02_dp; coef(1) = 0.206436961917E-02_dp
            alpha(2) = 0.161772821811E-01_dp; coef(2) = 0.308460403880E+00_dp
            alpha(3) = 0.239306136185E-01_dp; coef(3) = 0.910003374354E+00_dp
            alpha(4) = 0.546247989358E-01_dp; coef(4) = -0.332279687588E+00_dp
            alpha(5) = 0.162395488787E+00_dp; coef(5) = 0.120409054874E-01_dp
         CASE (6) ! Final Error      1.5047603163509393E-009         Nfun:         197
            alpha(1) = 0.870001973021E-02_dp; coef(1) = 0.179237366982E-02_dp
            alpha(2) = 0.161294730998E-01_dp; coef(2) = 0.304341866264E+00_dp
            alpha(3) = 0.238974541743E-01_dp; coef(3) = 0.913709005226E+00_dp
            alpha(4) = 0.546831680853E-01_dp; coef(4) = -0.331785002293E+00_dp
            alpha(5) = 0.159962733647E+00_dp; coef(5) = 0.120202087140E-01_dp
            alpha(6) = 0.251779315045E+00_dp; coef(6) = 0.214722618963E-03_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (40) ! 9f
         SELECT CASE (n)
         CASE (1) ! Final Error      3.5866854898298772E-003         Nfun:          30
            alpha(1) = 0.236451214456E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      5.7884916918893303E-004         Nfun:          56
            alpha(1) = 0.243065324818E-01_dp; coef(1) = 0.102101453863E+01_dp
            alpha(2) = 0.101415832337E+00_dp; coef(2) = -0.663738453708E-01_dp
         CASE (3) ! Final Error      5.0962693998179986E-007         Nfun:          82
            alpha(1) = 0.181242665872E-01_dp; coef(1) = 0.318224500768E+00_dp
            alpha(2) = 0.285618680076E-01_dp; coef(2) = 0.765670747468E+00_dp
            alpha(3) = 0.749787638189E-01_dp; coef(3) = -0.118327893498E+00_dp
         CASE (4) ! Final Error      3.6966872118959770E-007         Nfun:        5000
            alpha(1) = 0.179393655766E-01_dp; coef(1) = 0.299539260646E+00_dp
            alpha(2) = 0.282203526140E-01_dp; coef(2) = 0.775036811710E+00_dp
            alpha(3) = 0.627633883519E-01_dp; coef(3) = 0.110203506108E+00_dp
            alpha(4) = 0.701650018606E-01_dp; coef(4) = -0.219105186325E+00_dp
         CASE (5) ! Final Error      2.8785715144863389E-008         Nfun:         208
            alpha(1) = 0.112002739189E-01_dp; coef(1) = 0.783936105666E-02_dp
            alpha(2) = 0.189972972767E-01_dp; coef(2) = 0.362667355976E+00_dp
            alpha(3) = 0.290110276248E-01_dp; coef(3) = 0.714799764594E+00_dp
            alpha(4) = 0.750038510116E-01_dp; coef(4) = -0.119828548388E+00_dp
            alpha(5) = 0.364343520874E+00_dp; coef(5) = 0.632058496643E-03_dp
         CASE (6) ! Final Error      2.7361617638277380E-008         Nfun:        5000
            alpha(1) = 0.747151007821E-02_dp; coef(1) = -0.508489048126E-02_dp
            alpha(2) = 0.834341100476E-02_dp; coef(2) = 0.919759286561E-02_dp
            alpha(3) = 0.188388804576E-01_dp; coef(3) = 0.359513481163E+00_dp
            alpha(4) = 0.289596243605E-01_dp; coef(4) = 0.721522410714E+00_dp
            alpha(5) = 0.750304334833E-01_dp; coef(5) = -0.119686679826E+00_dp
            alpha(6) = 0.363852001112E+00_dp; coef(6) = 0.634389826897E-03_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (41) ! 9g
         SELECT CASE (n)
         CASE (1) ! Final Error      3.2311360041026305E-003         Nfun:          28
            alpha(1) = 0.290448768988E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      4.0309036684949007E-004         Nfun:          58
            alpha(1) = 0.184555219268E-01_dp; coef(1) = 0.224411382739E+00_dp
            alpha(2) = 0.323633375144E-01_dp; coef(2) = 0.810028895258E+00_dp
         CASE (3) ! Final Error      3.7960899590405717E-006         Nfun:          85
            alpha(1) = 0.211235010765E-01_dp; coef(1) = 0.390218406191E+00_dp
            alpha(2) = 0.353165989409E-01_dp; coef(2) = 0.659927541501E+00_dp
            alpha(3) = 0.118631219816E+00_dp; coef(3) = -0.300201339669E-01_dp
         CASE (4) ! Final Error      1.0936530059346694E-007         Nfun:         115
            alpha(1) = 0.159350384599E-01_dp; coef(1) = 0.627750934903E-01_dp
            alpha(2) = 0.244502567718E-01_dp; coef(2) = 0.492305494010E+00_dp
            alpha(3) = 0.377511852441E-01_dp; coef(3) = 0.501102410850E+00_dp
            alpha(4) = 0.110460493417E+00_dp; coef(4) = -0.352406429893E-01_dp
         CASE (5) ! Final Error      1.0707199282755653E-007         Nfun:        5000
            alpha(1) = 0.128314035490E-01_dp; coef(1) = -0.246201451659E-01_dp
            alpha(2) = 0.141795688822E-01_dp; coef(2) = 0.707044826121E-01_dp
            alpha(3) = 0.241033889047E-01_dp; coef(3) = 0.498271439917E+00_dp
            alpha(4) = 0.376170364805E-01_dp; coef(4) = 0.511689743479E+00_dp
            alpha(5) = 0.110592915547E+00_dp; coef(5) = -0.351198855531E-01_dp
         CASE (6) ! Final Error      4.3275267813669416E-009         Nfun:         174
            alpha(1) = 0.169601718642E-01_dp; coef(1) = 0.992638967613E-01_dp
            alpha(2) = 0.258271395411E-01_dp; coef(2) = 0.532122570919E+00_dp
            alpha(3) = 0.392098811035E-01_dp; coef(3) = 0.429410181833E+00_dp
            alpha(4) = 0.100427254790E+00_dp; coef(4) = -0.364948541963E-01_dp
            alpha(5) = 0.170035053051E+00_dp; coef(5) = -0.351897199401E-02_dp
            alpha(6) = 0.793628524125E+00_dp; coef(6) = 0.366533419657E-04_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (42) ! 9h
         SELECT CASE (n)
         CASE (1) ! Final Error      1.1167176664406142E-002         Nfun:          28
            alpha(1) = 0.344685478967E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      6.8510386961317948E-005         Nfun:          58
            alpha(1) = 0.251457484930E-01_dp; coef(1) = 0.503639183505E+00_dp
            alpha(2) = 0.454325437047E-01_dp; coef(2) = 0.563469443918E+00_dp
         CASE (3) ! Final Error      1.8408975835715324E-005         Nfun:          77
            alpha(1) = 0.178507120540E-01_dp; coef(1) = 0.843155675708E-01_dp
            alpha(2) = 0.290680331860E-01_dp; coef(2) = 0.562740517373E+00_dp
            alpha(3) = 0.486416415762E-01_dp; coef(3) = 0.423883783920E+00_dp
         CASE (4) ! Final Error      2.1966773844361853E-007         Nfun:         114
            alpha(1) = 0.200910816379E-01_dp; coef(1) = 0.170082883593E+00_dp
            alpha(2) = 0.322821386101E-01_dp; coef(2) = 0.596685776818E+00_dp
            alpha(3) = 0.526413846701E-01_dp; coef(3) = 0.307378666817E+00_dp
            alpha(4) = 0.184393487575E+00_dp; coef(4) = -0.620858825903E-02_dp
         CASE (5) ! Final Error      1.9557704344919102E-008         Nfun:         137
            alpha(1) = 0.159683939401E-01_dp; coef(1) = 0.297824322771E-01_dp
            alpha(2) = 0.239574844043E-01_dp; coef(2) = 0.293299673556E+00_dp
            alpha(3) = 0.358631388520E-01_dp; coef(3) = 0.523064519812E+00_dp
            alpha(4) = 0.554715180474E-01_dp; coef(4) = 0.229020307265E+00_dp
            alpha(5) = 0.174061705184E+00_dp; coef(5) = -0.705759779206E-02_dp
         CASE (6) ! Final Error      6.9780936584606934E-010         Nfun:         185
            alpha(1) = 0.168559281324E-01_dp; coef(1) = 0.469638749404E-01_dp
            alpha(2) = 0.253308929084E-01_dp; coef(2) = 0.346139318258E+00_dp
            alpha(3) = 0.377673118695E-01_dp; coef(3) = 0.495506766374E+00_dp
            alpha(4) = 0.575806728406E-01_dp; coef(4) = 0.187620357184E+00_dp
            alpha(5) = 0.157214165657E+00_dp; coef(5) = -0.756380602750E-02_dp
            alpha(6) = 0.306486717171E+00_dp; coef(6) = -0.623803960659E-03_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (43) ! 9i
         SELECT CASE (n)
         CASE (1) ! Final Error      2.3689235998586300E-002         Nfun:          31
            alpha(1) = 0.399105521873E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      4.5496038034653630E-004         Nfun:          55
            alpha(1) = 0.295580534767E-01_dp; coef(1) = 0.593981286347E+00_dp
            alpha(2) = 0.578889901228E-01_dp; coef(2) = 0.502956913677E+00_dp
         CASE (3) ! Final Error      7.5259707621455618E-006         Nfun:          81
            alpha(1) = 0.240893745271E-01_dp; coef(1) = 0.275981205704E+00_dp
            alpha(2) = 0.413540745940E-01_dp; coef(2) = 0.630778447691E+00_dp
            alpha(3) = 0.748077719001E-01_dp; coef(3) = 0.201877111693E+00_dp
         CASE (4) ! Final Error      1.8243439536824724E-007         Nfun:         102
            alpha(1) = 0.201737742424E-01_dp; coef(1) = 0.100332915029E+00_dp
            alpha(2) = 0.319785024939E-01_dp; coef(2) = 0.475110362981E+00_dp
            alpha(3) = 0.510717423884E-01_dp; coef(3) = 0.435843680539E+00_dp
            alpha(4) = 0.862948975252E-01_dp; coef(4) = 0.986301945129E-01_dp
         CASE (5) ! Final Error      8.8678959953014669E-008         Nfun:         132
            alpha(1) = 0.203495702186E-01_dp; coef(1) = 0.106257112776E+00_dp
            alpha(2) = 0.323489564862E-01_dp; coef(2) = 0.485290988284E+00_dp
            alpha(3) = 0.518116731159E-01_dp; coef(3) = 0.426725886325E+00_dp
            alpha(4) = 0.876363333272E-01_dp; coef(4) = 0.917412903436E-01_dp
            alpha(5) = 0.429175342128E+00_dp; coef(5) = -0.343597906434E-03_dp
         CASE (6) ! Final Error      9.8620400734716895E-010         Nfun:         163
            alpha(1) = 0.176126553735E-01_dp; coef(1) = 0.356142042204E-01_dp
            alpha(2) = 0.265571556005E-01_dp; coef(2) = 0.280406572602E+00_dp
            alpha(3) = 0.396805511195E-01_dp; coef(3) = 0.475897619059E+00_dp
            alpha(4) = 0.605783936792E-01_dp; coef(4) = 0.268114480926E+00_dp
            alpha(5) = 0.966879484399E-01_dp; coef(5) = 0.502715249385E-01_dp
            alpha(6) = 0.333829877064E+00_dp; coef(6) = -0.541522750889E-03_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (44) ! 9k
         SELECT CASE (n)
         CASE (1) ! Final Error      3.8748038276702701E-002         Nfun:          30
            alpha(1) = 0.453669165294E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      1.3169626327722561E-003         Nfun:          59
            alpha(1) = 0.334895210540E-01_dp; coef(1) = 0.626118257597E+00_dp
            alpha(2) = 0.689848257957E-01_dp; coef(2) = 0.496677586083E+00_dp
         CASE (3) ! Final Error      4.9719370870615620E-005         Nfun:          75
            alpha(1) = 0.275115841886E-01_dp; coef(1) = 0.332044884932E+00_dp
            alpha(2) = 0.493456345687E-01_dp; coef(2) = 0.629036646023E+00_dp
            alpha(3) = 0.962703662345E-01_dp; coef(3) = 0.181940947342E+00_dp
         CASE (4) ! Final Error      2.0711119399496880E-006         Nfun:         111
            alpha(1) = 0.237393163999E-01_dp; coef(1) = 0.165166797322E+00_dp
            alpha(2) = 0.395303853809E-01_dp; coef(2) = 0.553244937874E+00_dp
            alpha(3) = 0.674713019152E-01_dp; coef(3) = 0.368187388808E+00_dp
            alpha(4) = 0.127917815446E+00_dp; coef(4) = 0.598886690687E-01_dp
         CASE (5) ! Final Error      9.1917908173755336E-008         Nfun:         139
            alpha(1) = 0.210689085098E-01_dp; coef(1) = 0.796318173638E-01_dp
            alpha(2) = 0.333936927595E-01_dp; coef(2) = 0.410598293202E+00_dp
            alpha(3) = 0.530601295988E-01_dp; coef(3) = 0.465704656520E+00_dp
            alpha(4) = 0.883303285687E-01_dp; coef(4) = 0.171826635928E+00_dp
            alpha(5) = 0.163567598229E+00_dp; coef(5) = 0.193597808880E-01_dp
         CASE (6) ! Final Error      4.0797716316376409E-009         Nfun:         176
            alpha(1) = 0.190356002937E-01_dp; coef(1) = 0.375463358934E-01_dp
            alpha(2) = 0.290884695505E-01_dp; coef(2) = 0.275392194343E+00_dp
            alpha(3) = 0.440695090027E-01_dp; coef(3) = 0.466380998480E+00_dp
            alpha(4) = 0.683849761340E-01_dp; coef(4) = 0.288141727085E+00_dp
            alpha(5) = 0.111687327365E+00_dp; coef(5) = 0.732462424004E-01_dp
            alpha(6) = 0.201494488599E+00_dp; coef(6) = 0.653481769685E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT
      CASE (45) ! 9l
         SELECT CASE (n)
         CASE (1) ! Final Error      5.5152499658008169E-002         Nfun:          28
            alpha(1) = 0.508347257404E-01_dp; coef(1) = 0.100000000000E+01_dp
         CASE (2) ! Final Error      2.6242963416316378E-003         Nfun:          54
            alpha(1) = 0.373073954899E-01_dp; coef(1) = 0.640672393832E+00_dp
            alpha(2) = 0.791555644687E-01_dp; coef(2) = 0.504374727065E+00_dp
         CASE (3) ! Final Error      1.3939440085208865E-004         Nfun:          81
            alpha(1) = 0.305676686410E-01_dp; coef(1) = 0.357021646712E+00_dp
            alpha(2) = 0.562394587234E-01_dp; coef(2) = 0.626392745228E+00_dp
            alpha(3) = 0.113548661618E+00_dp; coef(3) = 0.190385558084E+00_dp
         CASE (4) ! Final Error      8.3561457435710196E-006         Nfun:         106
            alpha(1) = 0.263657375595E-01_dp; coef(1) = 0.189266601464E+00_dp
            alpha(2) = 0.449557719142E-01_dp; coef(2) = 0.564151044943E+00_dp
            alpha(3) = 0.789244640149E-01_dp; coef(3) = 0.362243128056E+00_dp
            alpha(4) = 0.156137313529E+00_dp; coef(4) = 0.640470311719E-01_dp
         CASE (5) ! Final Error      5.6219615807151513E-007         Nfun:         136
            alpha(1) = 0.234378071768E-01_dp; coef(1) = 0.988329543303E-01_dp
            alpha(2) = 0.380161165399E-01_dp; coef(2) = 0.439153630304E+00_dp
            alpha(3) = 0.619990036683E-01_dp; coef(3) = 0.456902356995E+00_dp
            alpha(4) = 0.106750312347E+00_dp; coef(4) = 0.165539804648E+00_dp
            alpha(5) = 0.208942371572E+00_dp; coef(5) = 0.205854955241E-01_dp
         CASE (6) ! Final Error      4.2032549218973259E-008         Nfun:         179
            alpha(1) = 0.212582275797E-01_dp; coef(1) = 0.516510174658E-01_dp
            alpha(2) = 0.332541608621E-01_dp; coef(2) = 0.316300050287E+00_dp
            alpha(3) = 0.516918539766E-01_dp; coef(3) = 0.471298406512E+00_dp
            alpha(4) = 0.827531017668E-01_dp; coef(4) = 0.268695845454E+00_dp
            alpha(5) = 0.141048916908E+00_dp; coef(5) = 0.668686375578E-01_dp
            alpha(6) = 0.274360292511E+00_dp; coef(6) = 0.651128382051E-02_dp
         CASE DEFAULT
            CPABORT("Illegal expansion")
         END SELECT

      CASE DEFAULT
         CPABORT("Illegal STO")
      END SELECT
      ! scaling: a = zeta^2 * a
      alpha(1:n) = zeta**2*alpha(1:n)

   END SUBROUTINE get_sto_ng

END MODULE sto_ng

