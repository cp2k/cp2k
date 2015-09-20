!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      06.2004 created
!> \author jgh
! *****************************************************************************
MODULE sto_ng
  
  USE kinds,                           ONLY: dp
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sto_ng'

  PUBLIC :: get_sto_ng

! *****************************************************************************
!> \brief return STO-NG parameters;
!>        INPUT: zeta (Slater exponent)
!>               n    (Expansion length)
!>               nq   (principle quantum number)
!>               lq   (angular quantum number)
!>        OUTPUT: alpha (n Gaussian exponents)
!>                coef  (n corresponding expansion coefficients)
!> \par History
!>      06.2004 created using parameters from : RF Stewart JCP 52 p431 1970
!>      01.2008 added newly generated parameters for 6s - 7i
!> \author jgh
CONTAINS

! *****************************************************************************
!> \brief ...
!> \param zeta ...
!> \param n ...
!> \param nq ...
!> \param lq ...
!> \param alpha ...
!> \param coef ...
! *****************************************************************************
SUBROUTINE get_sto_ng ( zeta, n, nq, lq, alpha, coef)

    REAL(KIND=dp), INTENT(IN)                :: zeta
    INTEGER, INTENT(IN)                      :: n, nq, lq
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: alpha, coef

    CHARACTER(len=*), PARAMETER :: routineN = 'get_sto_ng', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: m, na, nc
    LOGICAL                                  :: failure

  failure=.FALSE.

  na = SIZE ( alpha )
  nc = SIZE ( coef )

  CPASSERT(n >= 1 .AND. n <= 6)
  CPASSERT(na >= n)
  CPASSERT(nc >= n)
  CPASSERT(zeta > 0.0_dp)
  CPASSERT(nq >= 1)
  CPASSERT(lq >= 0)

  m = (nq*(nq-1))/2 + lq + 1
  SELECT CASE ( m )
    CASE (  1 ) ! 1s
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 2.709496091e-1_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 8.518186635e-1_dp; coef ( 1 ) = 4.301284983e-1_dp
           alpha( 2 ) = 1.516232927e-1_dp; coef ( 2 ) = 6.789135305e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 2.227660584e+0_dp; coef ( 1 ) = 1.543289673e-1_dp
           alpha( 2 ) = 4.057711562e-1_dp; coef ( 2 ) = 5.353281423e-1_dp
           alpha( 3 ) = 1.098175104e-1_dp; coef ( 3 ) = 4.446345422e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 5.216844534e+0_dp; coef ( 1 ) = 5.675242080e-2_dp
           alpha( 2 ) = 9.546182760e-1_dp; coef ( 2 ) = 2.601413550e-1_dp
           alpha( 3 ) = 2.652034102e-1_dp; coef ( 3 ) = 5.328461143e-1_dp
           alpha( 4 ) = 8.801862774e-2_dp; coef ( 4 ) = 2.916254405e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 1.130563696e+1_dp; coef ( 1 ) = 2.214055312e-2_dp
           alpha( 2 ) = 2.071728178e+0_dp; coef ( 2 ) = 1.135411520e-1_dp
           alpha( 3 ) = 5.786484833e-1_dp; coef ( 3 ) = 3.318161484e-1_dp
           alpha( 4 ) = 1.975724573e-1_dp; coef ( 4 ) = 4.825700713e-1_dp
           alpha( 5 ) = 7.445271746e-2_dp; coef ( 5 ) = 1.935721966e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 2.310303149e+1_dp; coef ( 1 ) = 9.163596280e-3_dp
           alpha( 2 ) = 4.235915534e+0_dp; coef ( 2 ) = 4.936149294e-2_dp
           alpha( 3 ) = 1.185056519e+0_dp; coef ( 3 ) = 1.685383049e-1_dp
           alpha( 4 ) = 4.070988982e-1_dp; coef ( 4 ) = 3.705627997e-1_dp
           alpha( 5 ) = 1.580884151e-1_dp; coef ( 5 ) = 4.164915298e-1_dp
           alpha( 6 ) = 6.510953954e-2_dp; coef ( 6 ) = 1.303340841e-1_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE (  2 ) ! 2s
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 1.012151084e-1_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 1.292278611e-1_dp; coef ( 1 ) = 7.470867124e-1_dp
           alpha( 2 ) = 4.908584205e-2_dp; coef ( 2 ) = 2.855980556e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 2.581578398e+0_dp; coef ( 1 ) =-5.994474934e-2_dp
           alpha( 2 ) = 1.567622104e-1_dp; coef ( 2 ) = 5.960385398e-1_dp
           alpha( 3 ) = 6.018332272e-2_dp; coef ( 3 ) = 4.581786291e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 1.161525551e+1_dp; coef ( 1 ) =-1.198411747e-2_dp
           alpha( 2 ) = 2.000243111e+0_dp; coef ( 2 ) =-5.472052539e-2_dp
           alpha( 3 ) = 1.607280687e-1_dp; coef ( 3 ) = 5.805587176e-1_dp
           alpha( 4 ) = 6.125744532e-2_dp; coef ( 4 ) = 4.770079976e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 8.984956862e+0_dp; coef ( 1 ) =-1.596349096e-2_dp
           alpha( 2 ) = 1.673710636e+0_dp; coef ( 2 ) =-5.685884883e-2_dp
           alpha( 3 ) = 1.944726668e-1_dp; coef ( 3 ) = 3.698265599e-1_dp
           alpha( 4 ) = 8.806345634e-2_dp; coef ( 4 ) = 5.480512593e-1_dp
           alpha( 5 ) = 4.249068522e-2_dp; coef ( 5 ) = 1.472634893e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 2.768496241e+1_dp; coef ( 1 ) =-4.151277819e-3_dp
           alpha( 2 ) = 5.077140627e+0_dp; coef ( 2 ) =-2.067024148e-2_dp
           alpha( 3 ) = 1.426786050e+0_dp; coef ( 3 ) =-5.150303337e-2_dp
           alpha( 4 ) = 2.040335729e-1_dp; coef ( 4 ) = 3.346271174e-1_dp
           alpha( 5 ) = 9.260298399e-2_dp; coef ( 5 ) = 5.621061301e-1_dp
           alpha( 6 ) = 4.416183978e-2_dp; coef ( 6 ) = 1.712994697e-1_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE (  3 ) ! 2p
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 1.759666885e-1_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 4.323908358e-1_dp; coef ( 1 ) = 4.522627513e-1_dp
           alpha( 2 ) = 1.069439065e-1_dp; coef ( 2 ) = 6.713122642e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 9.192379002e-1_dp; coef ( 1 ) = 1.623948553e-1_dp
           alpha( 2 ) = 2.359194503e-1_dp; coef ( 2 ) = 5.661708862e-1_dp
           alpha( 3 ) = 8.009805746e-2_dp; coef ( 3 ) = 4.223071752e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 1.798260992e+0_dp; coef ( 1 ) = 5.713170255e-2_dp
           alpha( 2 ) = 4.662622228e-1_dp; coef ( 2 ) = 2.857455515e-1_dp
           alpha( 3 ) = 1.643718620e-1_dp; coef ( 3 ) = 5.517873105e-1_dp
           alpha( 4 ) = 6.543927065e-2_dp; coef ( 4 ) = 2.632314924e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 3.320386533e+0_dp; coef ( 1 ) = 2.079051117e-2_dp
           alpha( 2 ) = 8.643257633e-1_dp; coef ( 2 ) = 1.235472099e-1_dp
           alpha( 3 ) = 3.079819284e-1_dp; coef ( 3 ) = 3.667738886e-1_dp
           alpha( 4 ) = 1.273309895e-1_dp; coef ( 4 ) = 4.834930290e-1_dp
           alpha( 5 ) = 5.606243164e-2_dp; coef ( 5 ) = 1.653444074e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 5.868285913e+0_dp; coef ( 1 ) = 7.924233646e-3_dp
           alpha( 2 ) = 1.530329631e+0_dp; coef ( 2 ) = 5.144104825e-2_dp
           alpha( 3 ) = 5.475665231e-1_dp; coef ( 3 ) = 1.898400060e-1_dp
           alpha( 4 ) = 2.288932733e-1_dp; coef ( 4 ) = 4.049863191e-1_dp
           alpha( 5 ) = 1.046655969e-1_dp; coef ( 5 ) = 4.012362861e-1_dp
           alpha( 6 ) = 4.948220127e-2_dp; coef ( 6 ) = 1.051855189e-1_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE (  4 ) ! 3s
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 5.296881757e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 6.694095822e-1_dp; coef ( 1 ) =-1.529645716e-1_dp
           alpha( 2 ) = 5.837135094e-2_dp; coef ( 2 ) = 1.051370110e+0_dp
         CASE ( 3 )
           alpha( 1 ) = 5.641487709e-1_dp; coef ( 1 ) =-1.782577972e-1_dp
           alpha( 2 ) = 6.924421391e-2_dp; coef ( 2 ) = 8.612761663e-1_dp
           alpha( 3 ) = 3.269529097e-2_dp; coef ( 3 ) = 2.261841969e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 1.513265591e+0_dp; coef ( 1 ) =-3.295496352e-2_dp
           alpha( 2 ) = 4.262497508e-1_dp; coef ( 2 ) =-1.724516959e-1_dp
           alpha( 3 ) = 7.643320863e-2_dp; coef ( 3 ) = 7.518511194e-1_dp
           alpha( 4 ) = 3.760545063e-2_dp; coef ( 4 ) = 3.589627310e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 4.275877914e+0_dp; coef ( 1 ) =-3.920358850e-3_dp
           alpha( 2 ) = 1.132409433e+0_dp; coef ( 2 ) =-4.168430506e-2_dp
           alpha( 3 ) = 4.016256968e-1_dp; coef ( 3 ) =-1.637440990e-1_dp
           alpha( 4 ) = 7.732370620e-2_dp; coef ( 4 ) = 7.419373723e-1_dp
           alpha( 5 ) = 3.800708627e-2_dp; coef ( 5 ) = 3.724364929e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 3.273031938e+0_dp; coef ( 1 ) =-6.775596947e-3_dp
           alpha( 2 ) = 9.200611311e-1_dp; coef ( 2 ) =-5.639325779e-2_dp
           alpha( 3 ) = 3.593349765e-1_dp; coef ( 3 ) =-1.587856086e-1_dp
           alpha( 4 ) = 8.636686991e-2_dp; coef ( 4 ) = 5.534527651e-1_dp
           alpha( 5 ) = 4.797373812e-2_dp; coef ( 5 ) = 5.015351020e-1_dp
           alpha( 6 ) = 2.724741144e-2_dp; coef ( 6 ) = 7.223633674e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE (  5 ) ! 3p
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 9.113614253e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 1.458620964e-1_dp; coef ( 1 ) = 5.349653144e-1_dp
           alpha( 2 ) = 5.664210742e-2_dp; coef ( 2 ) = 5.299607212e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 2.692880368e+0_dp; coef ( 1 ) =-1.061945788e-2_dp
           alpha( 2 ) = 1.489359592e-1_dp; coef ( 2 ) = 5.218564264e-1_dp
           alpha( 3 ) = 5.739585040e-2_dp; coef ( 3 ) = 5.450015143e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 1.853180239e+0_dp; coef ( 1 ) =-1.434249391e-2_dp
           alpha( 2 ) = 1.915075719e-1_dp; coef ( 2 ) = 2.755177589e-1_dp
           alpha( 3 ) = 8.655487938e-2_dp; coef ( 3 ) = 5.846750879e-1_dp
           alpha( 4 ) = 4.184253862e-2_dp; coef ( 4 ) = 2.144986514e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 6.466803859e+0_dp; coef ( 1 ) =-2.329023747e-3_dp
           alpha( 2 ) = 1.555914802e+0_dp; coef ( 2 ) =-1.357395221e-2_dp
           alpha( 3 ) = 1.955925255e-1_dp; coef ( 3 ) = 2.632185383e-1_dp
           alpha( 4 ) = 8.809647701e-2_dp; coef ( 4 ) = 5.880427024e-1_dp
           alpha( 5 ) = 4.234835707e-2_dp; coef ( 5 ) = 2.242794445e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 5.077973607e+0_dp; coef ( 1 ) =-3.329929840e-3_dp
           alpha( 2 ) = 1.340786940e+0_dp; coef ( 2 ) =-1.419488340e-2_dp
           alpha( 3 ) = 2.248434849e-1_dp; coef ( 3 ) = 1.639395770e-1_dp
           alpha( 4 ) = 1.131741848e-1_dp; coef ( 4 ) = 4.485358256e-1_dp
           alpha( 5 ) = 6.076408893e-2_dp; coef ( 5 ) = 3.908813050e-1_dp
           alpha( 6 ) = 3.315424265e-2_dp; coef ( 6 ) = 7.411456232e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE (  6 ) ! 3d
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 1.302270363e-1_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 2.777427345e-1_dp; coef ( 1 ) = 4.666137923e-1_dp
           alpha( 2 ) = 8.336507714e-2_dp; coef ( 2 ) = 6.644706516e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 5.229112225e-1_dp; coef ( 1 ) = 1.696596060e-1_dp
           alpha( 2 ) = 1.639595876e-1_dp; coef ( 2 ) = 5.847984817e-1_dp
           alpha( 3 ) = 6.386630021e-2_dp; coef ( 3 ) = 4.056779523e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 9.185846715e-1_dp; coef ( 1 ) = 5.799057705e-2_dp
           alpha( 2 ) = 2.920461109e-1_dp; coef ( 2 ) = 3.045581349e-1_dp
           alpha( 3 ) = 1.187568890e-1_dp; coef ( 3 ) = 5.601358038e-1_dp
           alpha( 4 ) = 5.286755896e-2_dp; coef ( 4 ) = 2.432423313e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 1.539033958e+0_dp; coef ( 1 ) = 2.020869128e-2_dp
           alpha( 2 ) = 4.922090297e-1_dp; coef ( 2 ) = 1.321157923e-1_dp
           alpha( 3 ) = 2.029756928e-1_dp; coef ( 3 ) = 3.911240346e-1_dp
           alpha( 4 ) = 9.424112917e-2_dp; coef ( 4 ) = 4.779609701e-1_dp
           alpha( 5 ) = 4.569058269e-2_dp; coef ( 5 ) = 1.463662294e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 2.488296923e+0_dp; coef ( 1 ) = 7.283828112e-3_dp
           alpha( 2 ) = 7.981487853e-1_dp; coef ( 2 ) = 5.386799363e-2_dp
           alpha( 3 ) = 3.311327490e-1_dp; coef ( 3 ) = 2.072139149e-1_dp
           alpha( 4 ) = 1.559114463e-1_dp; coef ( 4 ) = 4.266269092e-1_dp
           alpha( 5 ) = 7.877734732e-2_dp; coef ( 5 ) = 3.843100204e-1_dp
           alpha( 6 ) = 4.058484363e-2_dp; coef ( 6 ) = 8.902827546e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE (  7 ) ! 4s
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 3.264600274e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 2.441785453e-1_dp; coef ( 1 ) =-3.046656896e-1_dp
           alpha( 2 ) = 4.051097664e-2_dp; coef ( 2 ) = 1.146877294e+0_dp
         CASE ( 3 )
           alpha( 1 ) = 2.267938753e-1_dp; coef ( 1 ) =-3.349048323e-1_dp
           alpha( 2 ) = 4.448178019e-2_dp; coef ( 2 ) = 1.056744667e+0_dp
           alpha( 3 ) = 2.195294664e-2_dp; coef ( 3 ) = 1.256661680e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 3.242212833e-1_dp; coef ( 1 ) =-1.120682822e-1_dp
           alpha( 2 ) = 1.663217177e-1_dp; coef ( 2 ) =-2.845426863e-1_dp
           alpha( 3 ) = 5.081097451e-2_dp; coef ( 3 ) = 8.909873788e-1_dp
           alpha( 4 ) = 2.829066600e-2_dp; coef ( 4 ) = 3.517811205e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 2.980263783e+0_dp; coef ( 1 ) = 1.513948997e-3_dp
           alpha( 2 ) = 3.792228833e-1_dp; coef ( 2 ) =-7.316801548e-2_dp
           alpha( 3 ) = 1.789717224e-1_dp; coef ( 3 ) =-3.143703799e-1_dp
           alpha( 4 ) = 5.002110360e-2_dp; coef ( 4 ) = 9.032615169e-1_dp
           alpha( 5 ) = 2.789361681e-2_dp; coef ( 5 ) = 3.294210848e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 3.232838646e+0_dp; coef ( 1 ) = 1.374817488e-3_dp
           alpha( 2 ) = 3.605788802e-1_dp; coef ( 2 ) =-8.666390043e-2_dp
           alpha( 3 ) = 1.717905487e-1_dp; coef ( 3 ) =-3.130627309e-1_dp
           alpha( 4 ) = 5.277666487e-2_dp; coef ( 4 ) = 7.812787397e-1_dp
           alpha( 5 ) = 3.163400284e-2_dp; coef ( 5 ) = 4.389247988e-1_dp
           alpha( 6 ) = 1.874093091e-2_dp; coef ( 6 ) = 2.487178756e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE (  8 ) ! 4p
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 5.578350235e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 6.190052680e-2_dp; coef ( 1 ) = 8.743116767e-1_dp
           alpha( 2 ) = 2.648418407e-2_dp; coef ( 2 ) = 1.513640107e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 4.859692220e-1_dp; coef ( 1 ) =-6.147823411e-2_dp
           alpha( 2 ) = 7.430216918e-2_dp; coef ( 2 ) = 6.604172234e-1_dp
           alpha( 3 ) = 3.653340923e-2_dp; coef ( 3 ) = 3.932639495e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 1.492607880e+0_dp; coef ( 1 ) =-6.035216774e-3_dp
           alpha( 2 ) = 4.327619272e-1_dp; coef ( 2 ) =-6.013310874e-2_dp
           alpha( 3 ) = 7.553156064e-2_dp; coef ( 3 ) = 6.451518200e-1_dp
           alpha( 4 ) = 3.706272183e-2_dp; coef ( 4 ) = 4.117923820e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 1.091977298e+0_dp; coef ( 1 ) =-1.143929558e-2_dp
           alpha( 2 ) = 3.719985051e-1_dp; coef ( 2 ) =-6.322651538e-2_dp
           alpha( 3 ) = 8.590019352e-2_dp; coef ( 3 ) = 4.398907721e-1_dp
           alpha( 4 ) = 4.786503860e-2_dp; coef ( 4 ) = 5.245859166e-1_dp
           alpha( 5 ) = 2.730479990e-2_dp; coef ( 5 ) = 1.017072253e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 2.389722618e+0_dp; coef ( 1 ) =-1.665913575e-3_dp
           alpha( 2 ) = 7.960947826e-1_dp; coef ( 2 ) =-1.657464971e-2_dp
           alpha( 3 ) = 3.415541380e-1_dp; coef ( 3 ) =-5.958513378e-2_dp
           alpha( 4 ) = 8.847434525e-2_dp; coef ( 4 ) = 4.053115554e-1_dp
           alpha( 5 ) = 4.958248334e-2_dp; coef ( 5 ) = 5.433958189e-1_dp
           alpha( 6 ) = 2.816929784e-2_dp; coef ( 6 ) = 1.204970491e-1_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE (  9 ) ! 4d
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 7.941656339e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 1.330958892e-1_dp; coef ( 1 ) = 4.932764167e-1_dp
           alpha( 2 ) = 5.272119659e-2_dp; coef ( 2 ) = 5.918727866e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 1.777717219e-1_dp; coef ( 1 ) = 2.308552718e-1_dp
           alpha( 2 ) = 8.040647350e-2_dp; coef ( 2 ) = 6.042409177e-1_dp
           alpha( 3 ) = 3.949855551e-2_dp; coef ( 3 ) = 2.595768926e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 1.995825422e+0_dp; coef ( 1 ) =-2.816702620e-3_dp
           alpha( 2 ) = 1.823461280e-1_dp; coef ( 2 ) = 2.177095871e-1_dp
           alpha( 3 ) = 8.197240896e-2_dp; coef ( 3 ) = 6.058047348e-1_dp
           alpha( 4 ) = 4.000634951e-2_dp; coef ( 4 ) = 2.717811257e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 1.522122079e+0_dp; coef ( 1 ) =-3.673711876e-3_dp
           alpha( 2 ) = 2.173041823e-1_dp; coef ( 2 ) = 1.167122499e-1_dp
           alpha( 3 ) = 1.084876577e-1_dp; coef ( 3 ) = 4.216476416e-1_dp
           alpha( 4 ) = 5.836797641e-2_dp; coef ( 4 ) = 4.547673415e-1_dp
           alpha( 5 ) = 3.206682246e-2_dp; coef ( 5 ) = 1.037803318e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 4.634239420e+0_dp; coef ( 1 ) =-4.749842876e-4_dp
           alpha( 2 ) = 1.341648295e+0_dp; coef ( 2 ) =-3.566777891e-3_dp
           alpha( 3 ) = 2.209593028e-1_dp; coef ( 3 ) = 1.108670481e-1_dp
           alpha( 4 ) = 1.101467943e-1_dp; coef ( 4 ) = 4.159646930e-1_dp
           alpha( 5 ) = 5.904190370e-2_dp; coef ( 5 ) = 4.621672517e-1_dp
           alpha( 6 ) = 3.232628887e-2_dp; coef ( 6 ) = 1.081250196e-1_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 10 ) ! 4f
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 1.033434062e-1_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 2.006693538e-1_dp; coef ( 1 ) = 4.769346276e-1_dp
           alpha( 2 ) = 6.865384900e-2_dp; coef ( 2 ) = 6.587383976e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 3.483826963e-1_dp; coef ( 1 ) = 1.737856685e-1_dp
           alpha( 2 ) = 1.249380537e-1_dp; coef ( 2 ) = 5.973380628e-1_dp
           alpha( 3 ) = 5.349995725e-2_dp; coef ( 3 ) = 3.929395614e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 5.691670217e-1_dp; coef ( 1 ) = 5.902730589e-2_dp
           alpha( 2 ) = 2.074585819e-1_dp; coef ( 2 ) = 3.191828952e-1_dp
           alpha( 3 ) = 9.298346885e-2_dp; coef ( 3 ) = 5.639423893e-1_dp
           alpha( 4 ) = 4.473508853e-2_dp; coef ( 4 ) = 2.284796537e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 8.925960415e-1_dp; coef ( 1 ) = 1.999839052e-2_dp
           alpha( 2 ) = 3.277589120e-1_dp; coef ( 2 ) = 1.395427440e-1_dp
           alpha( 3 ) = 1.492869962e-1_dp; coef ( 3 ) = 4.091508237e-1_dp
           alpha( 4 ) = 7.506099109e-2_dp; coef ( 4 ) = 4.708252119e-1_dp
           alpha( 5 ) = 3.892475795e-2_dp; coef ( 5 ) = 1.328082566e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 1.357718039e+0_dp; coef ( 1 ) = 6.930234381e-3_dp
           alpha( 2 ) = 5.004907278e-1_dp; coef ( 2 ) = 5.634263745e-2_dp
           alpha( 3 ) = 2.296565064e-1_dp; coef ( 3 ) = 2.217065797e-1_dp
           alpha( 4 ) = 1.173146814e-1_dp; coef ( 4 ) = 4.411388883e-1_dp
           alpha( 5 ) = 6.350097171e-2_dp; coef ( 5 ) = 3.688112625e-1_dp
           alpha( 6 ) = 3.474556673e-2_dp; coef ( 6 ) = 7.787514504e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 11 ) ! 5s
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 2.216912938e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 1.213425654e-1_dp; coef ( 1 ) =-5.114756049e-1_dp
           alpha( 2 ) = 3.133152144e-2_dp; coef ( 2 ) = 1.307377277e+0_dp
         CASE ( 3 )
           alpha( 1 ) = 1.080198458e-1_dp; coef ( 1 ) =-6.617401158e-1_dp
           alpha( 2 ) = 4.408119328e-2_dp; coef ( 2 ) = 7.467595004e-1_dp
           alpha( 3 ) = 2.610811810e-2_dp; coef ( 3 ) = 7.146490945e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 8.602284252e-1_dp; coef ( 1 ) = 1.103657561e-2_dp
           alpha( 2 ) = 1.189050200e-1_dp; coef ( 2 ) =-5.606519023e-1_dp
           alpha( 3 ) = 3.446076176e-2_dp; coef ( 3 ) = 1.179429987e+0_dp
           alpha( 4 ) = 1.974798796e-2_dp; coef ( 4 ) = 1.734974376e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 7.403763257e-1_dp; coef ( 1 ) = 1.375523371e-2_dp
           alpha( 2 ) = 1.367990863e-1_dp; coef ( 2 ) =-3.097344179e-1_dp
           alpha( 3 ) = 9.135301779e-2_dp; coef ( 3 ) =-3.199192259e-1_dp
           alpha( 4 ) = 3.726907315e-2_dp; coef ( 4 ) = 1.084547038e+0_dp
           alpha( 5 ) = 2.241490836e-2_dp; coef ( 5 ) = 3.345288361e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 1.410128298e+0_dp; coef ( 1 ) = 2.695439582e-3_dp
           alpha( 2 ) = 5.077878915e-1_dp; coef ( 2 ) = 1.850157487e-2_dp
           alpha( 3 ) = 1.847926858e-1_dp; coef ( 3 ) =-9.588628125e-2_dp
           alpha( 4 ) = 1.061070594e-1_dp; coef ( 4 ) =-5.200673560e-1_dp
           alpha( 5 ) = 3.669584901e-2_dp; coef ( 5 ) = 1.087619490e+0_dp
           alpha( 6 ) = 2.213558430e-2_dp; coef ( 6 ) = 3.103964343e-1_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 12 ) ! 5p
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 3.769845216e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 2.691294191e-1_dp; coef ( 1 ) =-1.034227010e-1_dp
           alpha( 2 ) = 3.980805011e-2_dp; coef ( 2 ) = 1.033376378e+0_dp
         CASE ( 3 )
           alpha( 1 ) = 2.127482317e-1_dp; coef ( 1 ) =-1.389529695e-1_dp
           alpha( 2 ) = 4.729648620e-2_dp; coef ( 2 ) = 8.076691064e-1_dp
           alpha( 3 ) = 2.604865324e-2_dp; coef ( 3 ) = 2.726029342e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 3.962838833e-1_dp; coef ( 1 ) =-1.801459207e-2_dp
           alpha( 2 ) = 1.838858552e-1_dp; coef ( 2 ) =-1.360777372e-1_dp
           alpha( 3 ) = 4.943555157e-2_dp; coef ( 3 ) = 7.533973719e-1_dp
           alpha( 4 ) = 2.750222273e-2_dp; coef ( 4 ) = 3.409304859e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 3.422168934e-1_dp; coef ( 1 ) =-3.113958289e-2_dp
           alpha( 2 ) = 1.665099900e-1_dp; coef ( 2 ) =-1.374007017e-1_dp
           alpha( 3 ) = 5.443732013e-2_dp; coef ( 3 ) = 5.573881018e-1_dp
           alpha( 4 ) = 3.367775277e-2_dp; coef ( 4 ) = 4.855428100e-1_dp
           alpha( 5 ) = 2.091949042e-2_dp; coef ( 5 ) = 6.605423564e-2_dp
         CASE ( 6 )
           alpha( 1 ) = 3.778623374e+0_dp; coef ( 1 ) = 1.163246387e-4_dp
           alpha( 2 ) = 3.499121109e-1_dp; coef ( 2 ) =-2.920771322e-2_dp
           alpha( 3 ) = 1.683175469e-1_dp; coef ( 3 ) =-1.381051233e-1_dp
           alpha( 4 ) = 5.404070736e-2_dp; coef ( 4 ) = 5.706134877e-1_dp
           alpha( 5 ) = 3.328911801e-2_dp; coef ( 5 ) = 4.768808140e-1_dp
           alpha( 6 ) = 2.063815019e-2_dp; coef ( 6 ) = 6.021665516e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 13 ) ! 5d
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 5.352200793e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 6.906014388e-2_dp; coef ( 1 ) = 6.539405185e-1_dp
           alpha( 2 ) = 3.399457777e-2_dp; coef ( 2 ) = 3.948945302e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 4.913352950e-1_dp; coef ( 1 ) =-2.010175008e-2_dp
           alpha( 2 ) = 7.329090601e-2_dp; coef ( 2 ) = 5.899370608e-1_dp
           alpha( 3 ) = 3.594209290e-2_dp; coef ( 3 ) = 4.658445960e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 4.230617826e-1_dp; coef ( 1 ) =-2.421626009e-2_dp
           alpha( 2 ) = 8.293863702e-2_dp; coef ( 2 ) = 3.937644956e-1_dp
           alpha( 3 ) = 4.590326388e-2_dp; coef ( 3 ) = 5.489520286e-1_dp
           alpha( 4 ) = 2.628744797e-2_dp; coef ( 4 ) = 1.190436963e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 9.702946470e-1_dp; coef ( 1 ) =-3.231527611e-3_dp
           alpha( 2 ) = 3.603270196e-1_dp; coef ( 2 ) =-2.434931372e-2_dp
           alpha( 3 ) = 8.668717752e-2_dp; coef ( 3 ) = 3.440817054e-1_dp
           alpha( 4 ) = 4.833708379e-2_dp; coef ( 4 ) = 5.693674376e-1_dp
           alpha( 5 ) = 2.751899341e-2_dp; coef ( 5 ) = 1.511340183e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 8.820520428e-1_dp; coef ( 1 ) =-4.097377019e-3_dp
           alpha( 2 ) = 3.410838409e-1_dp; coef ( 2 ) =-2.508271857e-2_dp
           alpha( 3 ) = 9.204308840e-2_dp; coef ( 3 ) = 2.648458555e-1_dp
           alpha( 4 ) = 5.472831774e-2_dp; coef ( 4 ) = 5.097437054e-1_dp
           alpha( 5 ) = 3.391202830e-2_dp; coef ( 5 ) = 2.654483467e-1_dp
           alpha( 6 ) = 2.108227374e-2_dp; coef ( 6 ) = 2.623132212e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 14 ) ! 5f
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 6.952785407e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 1.156094555e-1_dp; coef ( 1 ) = 4.856637346e-1_dp
           alpha( 2 ) = 4.778940916e-2_dp; coef ( 2 ) = 6.125980914e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 1.649233885e-1_dp; coef ( 1 ) = 1.909729355e-1_dp
           alpha( 2 ) = 7.487066646e-2_dp; coef ( 2 ) = 6.146060459e-1_dp
           alpha( 3 ) = 3.735787219e-2_dp; coef ( 3 ) = 3.059611271e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 2.017831152e-1_dp; coef ( 1 ) = 9.174268830e-2_dp
           alpha( 2 ) = 1.001952178e-1_dp; coef ( 2 ) = 4.023496947e-1_dp
           alpha( 3 ) = 5.447006630e-2_dp; coef ( 3 ) = 4.937432100e-1_dp
           alpha( 4 ) = 3.037569283e-2_dp; coef ( 4 ) = 1.254001522e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 1.670735676e+0_dp; coef ( 1 ) =-7.301193568e-4_dp
           alpha( 2 ) = 2.072477219e-1_dp; coef ( 2 ) = 8.414991343e-2_dp
           alpha( 3 ) = 1.024709357e-1_dp; coef ( 3 ) = 3.923683153e-1_dp
           alpha( 4 ) = 5.537913898e-2_dp; coef ( 4 ) = 5.040033146e-1_dp
           alpha( 5 ) = 3.072866652e-2_dp; coef ( 5 ) = 1.328979300e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 1.334096840e+0_dp; coef ( 1 ) =-9.486751531e-4_dp
           alpha( 2 ) = 2.372312347e-1_dp; coef ( 2 ) = 4.624275998e-2_dp
           alpha( 3 ) = 1.269485744e-1_dp; coef ( 3 ) = 2.373699784e-1_dp
           alpha( 4 ) = 7.290318381e-2_dp; coef ( 4 ) = 4.589112231e-1_dp
           alpha( 5 ) = 4.351355997e-2_dp; coef ( 5 ) = 3.205010548e-1_dp
           alpha( 6 ) = 2.598071843e-2_dp; coef ( 6 ) = 5.077063693e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 15 ) ! 5g
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 8.565417784e-2_dp; coef ( 1 ) = 1.000000000e+0_dp
         CASE ( 2 )
           alpha( 1 ) = 1.554531559e-1_dp; coef ( 1 ) = 4.848298074e-1_dp
           alpha( 2 ) = 5.854079811e-2_dp; coef ( 2 ) = 6.539381621e-1_dp
         CASE ( 3 )
           alpha( 1 ) = 2.545432122e-1_dp; coef ( 1 ) = 1.780980905e-1_dp
           alpha( 2 ) = 1.006544376e-1_dp; coef ( 2 ) = 6.063757846e-1_dp
           alpha( 3 ) = 4.624463922e-2_dp; coef ( 3 ) = 3.828552923e-1_dp
         CASE ( 4 )
           alpha( 1 ) = 3.945205573e-1_dp; coef ( 1 ) = 6.010484250e-2_dp
           alpha( 2 ) = 1.588100623e-1_dp; coef ( 2 ) = 3.309738329e-1_dp
           alpha( 3 ) = 7.646521729e-2_dp; coef ( 3 ) = 5.655207585e-1_dp
           alpha( 4 ) = 3.898703611e-2_dp; coef ( 4 ) = 2.171122608e-1_dp
         CASE ( 5 )
           alpha( 1 ) = 5.895429375e-1_dp; coef ( 1 ) = 1.998085812e-2_dp
           alpha( 2 ) = 2.393343780e-1_dp; coef ( 2 ) = 1.460384050e-1_dp
           alpha( 3 ) = 1.172646904e-1_dp; coef ( 3 ) = 4.230565459e-1_dp
           alpha( 4 ) = 6.254074479e-2_dp; coef ( 4 ) = 4.635699665e-1_dp
           alpha( 5 ) = 3.411243214e-2_dp; coef ( 5 ) = 1.226411691e-1_dp
         CASE ( 6 )
           alpha( 1 ) = 8.574668996e-1_dp; coef ( 1 ) = 6.729778096e-3_dp
           alpha( 2 ) = 3.497184772e-1_dp; coef ( 2 ) = 5.874145170e-2_dp
           alpha( 3 ) = 1.727917060e-1_dp; coef ( 3 ) = 2.339955227e-1_dp
           alpha( 4 ) = 9.373643151e-2_dp; coef ( 4 ) = 4.512983737e-1_dp
           alpha( 5 ) = 5.340032759e-2_dp; coef ( 5 ) = 3.552053926e-1_dp
           alpha( 6 ) = 3.057364464e-2_dp; coef ( 6 ) = 6.974153145e-2_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 16 ) ! 6s
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.016053956888690133_dp; coef ( 1 ) = 1.000000000000000_dp
         CASE ( 2 )
           alpha( 1 ) = 0.026188993287017705_dp; coef ( 1 ) = 1.6061745480037777_dp
           alpha( 2 ) = 0.06916086082602657_dp;  coef ( 2 ) = -0.8493191121412285_dp
         CASE ( 3 )
           alpha( 1 ) = 0.025012605890022674_dp; coef ( 1 ) = 1.4570914918085252_dp
           alpha( 2 ) = 0.08530028544682955_dp;  coef ( 2 ) = -0.7792202598134607_dp
           alpha( 3 ) = 0.23664111527549864_dp;  coef ( 3 ) = 0.08944830584494995_dp
         CASE ( 4 )
           alpha( 1 ) = 0.018111879224996335_dp; coef ( 1 ) = 0.2742426110211285_dp
           alpha( 2 ) = 0.028489853269487397_dp; coef ( 2 ) = 1.3154513257751648_dp
           alpha( 3 ) = 0.07575858742450885_dp;  coef ( 3 ) = -0.8700267823821266_dp
           alpha( 4 ) = 0.3135211502054625_dp;   coef ( 4 ) = 0.04933652686369165_dp
         CASE ( 5 )
           alpha( 1 ) = 0.017704122931955696_dp; coef ( 1 ) = 0.23979549908360312_dp
           alpha( 2 ) = 0.02816438780280081_dp;  coef ( 2 ) = 1.339995854967376_dp
           alpha( 3 ) = 0.07635157423691516_dp;  coef ( 3 ) = -0.8631396815113626_dp
           alpha( 4 ) = 0.29378907244008734_dp;  coef ( 4 ) = 0.05070584407672884_dp
           alpha( 5 ) = 0.7609367674925651_dp;   coef ( 5 ) = 0.0017452546896168622_dp
         CASE ( 6 )
           alpha( 1 ) = 0.018375086251195376_dp; coef ( 1 ) = 0.3019609612694464_dp
           alpha( 2 ) = 0.02883213591739931_dp;  coef ( 2 ) = 1.3047682665168656_dp
           alpha( 3 ) = 0.0735171952084867_dp;   coef ( 3 ) = -0.8423706115013339_dp
           alpha( 4 ) = 0.13718508837897095_dp;  coef ( 4 ) = -0.08693662797569257_dp
           alpha( 5 ) = 0.21928060661775448_dp;  coef ( 5 ) = 0.08040482815742164_dp
           alpha( 6 ) = 0.47255402578342814_dp;  coef ( 6 ) = 0.011258177324556012_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 17 ) ! 6p
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.02720326117498455_dp;  coef ( 1 ) = 1.0000000000000000_dp
         CASE ( 2 )
           alpha( 1 ) = 0.03068554921378584_dp;  coef ( 1 ) = 1.0932524470142648_dp
           alpha( 2 ) = 0.13327144868476606_dp;  coef ( 2 ) = -0.20022333756955418_dp
         CASE ( 3 )
           alpha( 1 ) = 0.021306449714000038_dp; coef ( 1 ) = 0.2589896573574047_dp
           alpha( 2 ) = 0.03558737683798757_dp;  coef ( 2 ) = 0.8948520093263641_dp
           alpha( 3 ) = 0.1159879775604416_dp;   coef ( 3 ) = -0.2541710812641519_dp
         CASE ( 4 )
           alpha( 1 ) = 0.02094853126286551_dp;  coef ( 1 ) = 0.23874968368587424_dp
           alpha( 2 ) = 0.03524009985593678_dp;  coef ( 2 ) = 0.9113069782415131_dp
           alpha( 3 ) = 0.1172381116493136_dp;   coef ( 3 ) = -0.25095169849099785_dp
           alpha( 4 ) = 1.074757975437711_dp;    coef ( 4 ) = 0.0014270633206840072_dp
         CASE ( 5 )
           alpha( 1 ) = 0.02191796448535917_dp;  coef ( 1 ) = 0.29958347534070895_dp
           alpha( 2 ) = 0.036406617593622265_dp; coef ( 2 ) = 0.8669095451744081_dp
           alpha( 3 ) = 0.10823682952141411_dp;  coef ( 3 ) = -0.24403207325952317_dp
           alpha( 4 ) = 0.1912844969312091_dp;   coef ( 4 ) = -0.025610699331187835_dp
           alpha( 5 ) = 0.605311988824089_dp;    coef ( 5 ) = 0.003533893449371835_dp
         CASE ( 6 )
           alpha( 1 ) = 0.01836987035282947_dp;  coef ( 1 ) = 0.0881865874022419_dp
           alpha( 2 ) = 0.028283616813330213_dp; coef ( 2 ) = 0.58260237820929_dp
           alpha( 3 ) = 0.0424808315746305_dp;   coef ( 3 ) = 0.5465800645782257_dp
           alpha( 4 ) = 0.08775028399997764_dp;  coef ( 4 ) = -0.20802232268176662_dp
           alpha( 5 ) = 0.14284922158351138_dp;  coef ( 5 ) = -0.1117378797434228_dp
           alpha( 6 ) = 0.6631101746841396_dp;   coef ( 6 ) = 0.0028390610573703423_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 18 ) ! 6d
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.038539962557578994_dp; coef ( 1 ) = 1.0000000000000000_dp
         CASE ( 2 )
           alpha( 1 ) = 0.018343252403866257_dp; coef ( 1 ) = 0.1002697090882244_dp
           alpha( 2 ) = 0.04095508659742878_dp;  coef ( 2 ) = 0.921707450678634_dp
         CASE ( 3 )
           alpha( 1 ) = 0.026815363906634464_dp; coef ( 1 ) = 0.36430995341001843_dp
           alpha( 2 ) = 0.048285664795074575_dp; coef ( 2 ) = 0.6902756637148024_dp
           alpha( 3 ) = 0.20398527022566643_dp;  coef ( 3 ) = -0.061643555476797514_dp
         CASE ( 4 )
           alpha( 1 ) = 0.014812592588885804_dp; coef ( 1 ) = 0.011887725165713237_dp
           alpha( 2 ) = 0.028604392422323988_dp; coef ( 2 ) = 0.4126202273551119_dp
           alpha( 3 ) = 0.049571759954916264_dp; coef ( 3 ) = 0.6324841000859261_dp
           alpha( 4 ) = 0.20042444873055187_dp;  coef ( 4 ) = -0.06355840930231531_dp
         CASE ( 5 )
           alpha( 1 ) = 0.02104884011402785_dp;  coef ( 1 ) = 0.09247565169723151_dp
           alpha( 2 ) = 0.03387333992146401_dp;  coef ( 2 ) = 0.5183469446932626_dp
           alpha( 3 ) = 0.05464871753445727_dp;  coef ( 3 ) = 0.4583625355467117_dp
           alpha( 4 ) = 0.16559124243612725_dp;  coef ( 4 ) = -0.06557358037625106_dp
           alpha( 5 ) = 0.3388543915705088_dp;   coef ( 5 ) = -0.01073326765431463_dp
         CASE ( 6 )
           alpha( 1 ) = 0.012538741501253952_dp; coef ( 1 ) = 0.0013891559057176176_dp
           alpha( 2 ) = 0.022411877391156923_dp; coef ( 2 ) = 0.12323543239586941_dp
           alpha( 3 ) = 0.03502206086139318_dp;  coef ( 3 ) = 0.516335137452715_dp
           alpha( 4 ) = 0.05542414129417389_dp;  coef ( 4 ) = 0.42911284700977503_dp
           alpha( 5 ) = 0.16410056143462523_dp;  coef ( 5 ) = -0.065949692922377_dp
           alpha( 6 ) = 0.33496788346358164_dp;  coef ( 6 ) = -0.01123082945384288_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 19 ) ! 6f
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.04999708511585302_dp;  coef ( 1 ) = 1._dp
         CASE ( 2 )
           alpha( 1 ) = 0.03443190219402774_dp;  coef ( 1 ) = 0.5095525299458369_dp
           alpha( 2 ) = 0.06993938754199006_dp;  coef ( 2 ) = 0.556840322894986_dp
         CASE ( 3 )
           alpha( 1 ) = 0.022913376030207443_dp; coef ( 1 ) = 0.0844980697596352_dp
           alpha( 2 ) = 0.040621057779302076_dp; coef ( 2 ) = 0.558907102542338_dp
           alpha( 3 ) = 0.07551452585114922_dp;  coef ( 3 ) = 0.4268730468790412_dp
         CASE ( 4 )
           alpha( 1 ) = 0.026723978066792376_dp; coef ( 1 ) = 0.18313299355518015_dp
           alpha( 2 ) = 0.04675842108187079_dp;  coef ( 2 ) = 0.5931630352672648_dp
           alpha( 3 ) = 0.08424111895397314_dp;  coef ( 3 ) = 0.29815135014219757_dp
           alpha( 4 ) = 0.40182223710248033_dp;  coef ( 4 ) = -0.009003129808147884_dp
         CASE ( 5 )
           alpha( 1 ) = 0.019820990273457767_dp; coef ( 1 ) = 0.026554164306668928_dp
           alpha( 2 ) = 0.031943455135755784_dp; coef ( 2 ) = 0.2857152467032053_dp
           alpha( 3 ) = 0.05180437300600356_dp;  coef ( 3 ) = 0.5249457199886216_dp
           alpha( 4 ) = 0.08839526561404966_dp;  coef ( 4 ) = 0.23802584890828332_dp
           alpha( 5 ) = 0.3861344493937384_dp;   coef ( 5 ) = -0.00959390690588398_dp
         CASE ( 6 )
           alpha( 1 ) = 0.021979238503298942_dp; coef ( 1 ) = 0.05529549645466706_dp
           alpha( 2 ) = 0.0353635379450325_dp;   coef ( 2 ) = 0.36260443662220726_dp
           alpha( 3 ) = 0.056802330945072296_dp; coef ( 3 ) = 0.47981873895702293_dp
           alpha( 4 ) = 0.09436112630657385_dp;  coef ( 4 ) = 0.17913120326169676_dp
           alpha( 5 ) = 0.3258167395765562_dp;   coef ( 5 ) = -0.010025342773471209_dp
           alpha( 6 ) = 0.7758137452357111_dp;   coef ( 6 ) = -0.0013665302206626205_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 20 ) ! 6g
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.061536068942275354_dp; coef ( 1 ) = 1._dp
         CASE ( 2 )
           alpha( 1 ) = 0.04332629465536664_dp;  coef ( 1 ) = 0.6209033605306379_dp
           alpha( 2 ) = 0.10013340560986961_dp;  coef ( 2 ) = 0.4864065117811783_dp
         CASE ( 3 )
           alpha( 1 ) = 0.03454561337035713_dp;  coef ( 1 ) = 0.3227633523004909_dp
           alpha( 2 ) = 0.06729191664458245_dp;  coef ( 2 ) = 0.6193784317147506_dp
           alpha( 3 ) = 0.14454284915249543_dp;  coef ( 3 ) = 0.1810539000206162_dp
         CASE ( 4 )
           alpha( 1 ) = 0.028983896307773004_dp; coef ( 1 ) = 0.15152805717139864_dp
           alpha( 2 ) = 0.05135519389407211_dp;  coef ( 2 ) = 0.5309701439600488_dp
           alpha( 3 ) = 0.09369371519055944_dp;  coef ( 3 ) = 0.37457949569535787_dp
           alpha( 4 ) = 0.18939419012133174_dp;  coef ( 4 ) = 0.06851383804287256_dp
         CASE ( 5 )
           alpha( 1 ) = 0.024699673908324754_dp; coef ( 1 ) = 0.06048210643335771_dp
           alpha( 2 ) = 0.04083036114944879_dp;  coef ( 2 ) = 0.3536959373883507_dp
           alpha( 3 ) = 0.06763537320097322_dp;  coef ( 3 ) = 0.46434612683985144_dp
           alpha( 4 ) = 0.1170091091566765_dp;   coef ( 4 ) = 0.2130170165201656_dp
           alpha( 5 ) = 0.22066124899602302_dp;  coef ( 5 ) = 0.0343912188428939_dp
         CASE ( 6 )
           alpha( 1 ) = 0.024944035346224053_dp; coef ( 1 ) = 0.06443914638193514_dp
           alpha( 2 ) = 0.04138718127148023_dp;  coef ( 2 ) = 0.36507904601268215_dp
           alpha( 3 ) = 0.06886079537637349_dp;  coef ( 3 ) = 0.4624056701486661_dp
           alpha( 4 ) = 0.11972435105059721_dp;  coef ( 4 ) = 0.20306165490832248_dp
           alpha( 5 ) = 0.22639206282070937_dp;  coef ( 5 ) = 0.03099725943908058_dp
           alpha( 6 ) = 1.4892687021530484_dp;   coef ( 6 ) = -0.00018343848780163463_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 21 ) ! 6h
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.07313314184000295_dp;  coef ( 1 ) = 1._dp
         CASE ( 2 )
           alpha( 1 ) = 0.05113303517913134_dp;  coef ( 1 ) = 0.6498709211538302_dp
           alpha( 2 ) = 0.12607272311510218_dp;  coef ( 2 ) = 0.4911257514330729_dp
         CASE ( 3 )
           alpha( 1 ) = 0.040851454032819305_dp; coef ( 1 ) = 0.3746467681927208_dp
           alpha( 2 ) = 0.08415648437016819_dp;  coef ( 2 ) = 0.6132065860909666_dp
           alpha( 3 ) = 0.1974425255793456_dp;   coef ( 3 ) = 0.18179709985670622_dp
         CASE ( 4 )
           alpha( 1 ) = 0.034679493065603105_dp; coef ( 1 ) = 0.20806488992498948_dp
           alpha( 2 ) = 0.06498507795614539_dp;  coef ( 2 ) = 0.5659152898854298_dp
           alpha( 3 ) = 0.1276799968878006_dp;   coef ( 3 ) = 0.3407357791588573_dp
           alpha( 4 ) = 0.293954336816017_dp;    coef ( 4 ) = 0.061163634708810435_dp
         CASE ( 5 )
           alpha( 1 ) = 0.030487493495559674_dp; coef ( 1 ) = 0.1147257294844731_dp
           alpha( 2 ) = 0.05371457065868773_dp;  coef ( 2 ) = 0.45668261354533374_dp
           alpha( 3 ) = 0.09619925104884211_dp;  coef ( 3 ) = 0.434116528713387_dp
           alpha( 4 ) = 0.18557752006634115_dp;  coef ( 4 ) = 0.15177021621261605_dp
           alpha( 5 ) = 0.42330005591830827_dp;  coef ( 5 ) = 0.02006866296723323_dp
         CASE ( 6 )
           alpha( 1 ) = 0.027419761875036298_dp; coef ( 1 ) = 0.06355638710922437_dp
           alpha( 2 ) = 0.04620886209790054_dp;  coef ( 2 ) = 0.34336050905161936_dp
           alpha( 3 ) = 0.07793999480997445_dp;  coef ( 3 ) = 0.45863493125119364_dp
           alpha( 4 ) = 0.13712767813357474_dp;  coef ( 4 ) = 0.24456067327385683_dp
           alpha( 5 ) = 0.2621356941062822_dp;   coef ( 5 ) = 0.061018224093793164_dp
           alpha( 6 ) = 0.5948121608710077_dp;   coef ( 6 ) = 0.006618011961621055_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 22 ) ! 7s
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.012170316100569041_dp; coef ( 1 ) = 1._dp

         CASE ( 2 )
           alpha( 1 ) = 0.02380792064880498_dp;  coef ( 1 ) = 2.4733498823505897_dp
           alpha( 2 ) = 0.0413963085083486_dp;   coef ( 2 ) = -1.7500626134752972_dp
         CASE ( 3 )
           alpha( 1 ) = 0.020979847986630586_dp; coef ( 1 ) = 1.675236183601737_dp
           alpha( 2 ) = 0.061178936229820555_dp; coef ( 2 ) = -1.147227387784115_dp
           alpha( 3 ) = 0.1401881324437541_dp;   coef ( 3 ) = 0.2121061484210291_dp
         CASE ( 4 )
           alpha( 1 ) = 0.015915368005026652_dp; coef ( 1 ) = 0.3275001969711621_dp
           alpha( 2 ) = 0.023863997843078046_dp; coef ( 2 ) = 1.5476507021661345_dp
           alpha( 3 ) = 0.054562166429946386_dp; coef ( 3 ) = -1.2697132490611722_dp
           alpha( 4 ) = 0.16342646591350343_dp;  coef ( 4 ) = 0.13574231777373438_dp
         CASE ( 5 )
           alpha( 1 ) = 0.01606210559622015_dp;  coef ( 1 ) = 0.3489870695674181_dp
           alpha( 2 ) = 0.024001717565201447_dp; coef ( 2 ) = 1.5335261688987925_dp
           alpha( 3 ) = 0.05439364645164688_dp;  coef ( 3 ) = -1.275299080068857_dp
           alpha( 4 ) = 0.16456835599844608_dp;  coef ( 4 ) = 0.13406783441420897_dp
           alpha( 5 ) = 1.8632404725730192_dp;   coef ( 5 ) = -0.0001818568135349486_dp
         CASE ( 6 )
           alpha( 1 ) = 0.01565231604983483_dp;  coef ( 1 ) = 0.29154442199853614_dp
           alpha( 2 ) = 0.023623292255693457_dp; coef ( 2 ) = 1.5691344169994597_dp
           alpha( 3 ) = 0.05497538487325693_dp;  coef ( 3 ) = -1.2618704056879249_dp
           alpha( 4 ) = 0.15604786033627724_dp;  coef ( 4 ) = 0.13683678315315853_dp
           alpha( 5 ) = 0.30526699650551353_dp;  coef ( 5 ) = 0.006442390493458186_dp
           alpha( 6 ) = 0.7613613898332257_dp;   coef ( 6 ) = -0.0009736214432373223_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 23 ) ! 7p
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.02056749089204262_dp;  coef ( 1 ) = 1.0000000000000000_dp
         CASE ( 2 )
           alpha( 1 ) = 0.02509434622887706_dp;  coef ( 1 ) = 1.1953936450048124_dp
           alpha( 2 ) = 0.07932687651759715_dp;  coef ( 2 ) = -0.33554981426999153_dp
         CASE ( 3 )
           alpha( 1 ) = 0.02216294295259616_dp;  coef ( 1 ) = 0.8184285461616602_dp
           alpha( 2 ) = 0.042045475539838104_dp; coef ( 2 ) = 0.8328277926936017_dp
           alpha( 3 ) = 0.0622825034958455_dp;   coef ( 3 ) = -0.7865628896072743_dp
         CASE ( 4 )
           alpha( 1 ) = 0.01798640692496774_dp;  coef ( 1 ) = 0.2517796607641491_dp
           alpha( 2 ) = 0.02838134507613473_dp;  coef ( 2 ) = 1.0084655616939668_dp
           alpha( 3 ) = 0.07596133909152862_dp;  coef ( 3 ) = -0.40416345513837804_dp
           alpha( 4 ) = 0.30794728288131357_dp;  coef ( 4 ) = 0.011606225585276365_dp
         CASE ( 5 )
           alpha( 1 ) = 0.018253221396022562_dp; coef ( 1 ) = 0.27532738928860534_dp
           alpha( 2 ) = 0.028708043693626412_dp; coef ( 2 ) = 0.9951443329321299_dp
           alpha( 3 ) = 0.07061876309146031_dp;  coef ( 3 ) = -0.2663034712401708_dp
           alpha( 4 ) = 0.08392634608758946_dp;  coef ( 4 ) = -0.1489548148570699_dp
           alpha( 5 ) = 0.298415801799801_dp;    coef ( 5 ) = 0.012451440490101195_dp
         CASE ( 6 )
           ! Warning: optimization not complete
           alpha( 1 ) = 0.018266628106424383_dp; coef ( 1 ) = 0.2762223367567883_dp
           alpha( 2 ) = 0.02870452146813995_dp;  coef ( 2 ) = 0.992858662003434_dp
           alpha( 3 ) = 0.07430061134905402_dp;  coef ( 3 ) = -0.4029967818464786_dp
           alpha( 4 ) = 0.18615871115244434_dp;  coef ( 4 ) = -41.116221124256434_dp
           alpha( 5 ) = 0.1861966605860832_dp;   coef ( 5 ) = 41.115741037337784_dp
           alpha( 6 ) = 0.476774875907382_dp;    coef ( 6 ) = 0.0020913658772264638_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 24 ) ! 7d
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.02909039015051073_dp;  coef ( 1 ) = 1.0000000000000002_dp
         CASE ( 2 )
           alpha( 1 ) = 0.03018482426518148_dp;  coef ( 1 ) = 1.025499895914138_dp
           alpha( 2 ) = 0.15242731708213653_dp;  coef ( 2 ) = -0.0801196219224754_dp
         CASE ( 3 )
           alpha( 1 ) = 0.021405281219812042_dp; coef ( 1 ) = 0.2999132044312777_dp
           alpha( 2 ) = 0.03568228986173447_dp;  coef ( 2 ) = 0.7809338560049855_dp
           alpha( 3 ) = 0.11574152434647746_dp;  coef ( 3 ) = -0.1242352213758718_dp
         CASE ( 4 )
           alpha( 1 ) = 0.021660976213634955_dp; coef ( 1 ) = 0.3177854736353611_dp
           alpha( 2 ) = 0.03606163962267541_dp;  coef ( 2 ) = 0.7675316218690252_dp
           alpha( 3 ) = 0.10985898319545016_dp;  coef ( 3 ) = -0.11712958841599931_dp
           alpha( 4 ) = 0.15908195611925_dp;     coef ( 4 ) = -0.011747385912275846_dp
         CASE ( 5 )
           ! Warning: optimization not complete
           alpha( 1 ) = 0.018501533387022245_dp; coef ( 1 ) = 0.11676114014037238_dp
           alpha( 2 ) = 0.02923074108604399_dp;  coef ( 2 ) = 0.6656817970408199_dp
           alpha( 3 ) = 0.05423405317769654_dp;  coef ( 3 ) = 600.6430437919528_dp
           alpha( 4 ) = 0.054241824018008936_dp; coef ( 4 ) = -600.3759617830652_dp
           alpha( 5 ) = 0.12494194517082237_dp;  coef ( 5 ) = -0.09295829453419228_dp
         CASE ( 6 )
           alpha( 1 ) = 0.017457157083607015_dp; coef ( 1 ) = 0.07250557935571439_dp
           alpha( 2 ) = 0.026593769405605233_dp; coef ( 2 ) = 0.5007763109536889_dp
           alpha( 3 ) = 0.04008932402169768_dp;  coef ( 3 ) = 0.5339957399786912_dp
           alpha( 4 ) = 0.09560872001992007_dp;  coef ( 4 ) = -0.11880032574322265_dp
           alpha( 5 ) = 0.15425745053701811_dp;  coef ( 5 ) = -0.03241452360588799_dp
           alpha( 6 ) = 0.5964979973145561_dp;   coef ( 6 ) = 0.0007012418324482315_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 25 ) ! 7f
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.037696619888011394_dp; coef ( 1 ) = 1._dp
         CASE ( 2 )
           alpha( 1 ) = 0.02399618488506292_dp;  coef ( 1 ) = 0.2938876408614713_dp
           alpha( 2 ) = 0.0441727483070652_dp;   coef ( 2 ) = 0.7461850531456674_dp
         CASE ( 3 )
           alpha( 1 ) = 0.02649072884749045_dp;  coef ( 1 ) = 0.41965629582857367_dp
           alpha( 2 ) = 0.04781482493797306_dp;  coef ( 2 ) = 0.6320266162094202_dp
           alpha( 3 ) = 0.20807007812612838_dp;  coef ( 3 ) = -0.025856724476949224_dp
         CASE ( 4 )
           alpha( 1 ) = 0.019516667217357033_dp; coef ( 1 ) = 0.07760665191295739_dp
           alpha( 2 ) = 0.031536542788218214_dp; coef ( 2 ) = 0.5138599419165162_dp
           alpha( 3 ) = 0.0518809354792556_dp;   coef ( 3 ) = 0.4658306981238731_dp
           alpha( 4 ) = 0.18978036627517764_dp;  coef ( 4 ) = -0.03032718718680631_dp
         CASE ( 5 )
           alpha( 1 ) = 0.02085272354115486_dp;  coef ( 1 ) = 0.11822768394863675_dp
           alpha( 2 ) = 0.033574888134901455_dp; coef ( 2 ) = 0.5490677063820284_dp
           alpha( 3 ) = 0.05432263605261335_dp;  coef ( 3 ) = 0.3941968220546063_dp
           alpha( 4 ) = 0.1670441095948162_dp;   coef ( 4 ) = -0.031436205691901344_dp
           alpha( 5 ) = 0.3442595341279306_dp;   coef ( 5 ) = -0.0033909812829636394_dp
         CASE ( 6 )
           alpha( 1 ) = 0.01641594422436568_dp;  coef ( 1 ) = 0.01709822656839514_dp
           alpha( 2 ) = 0.024762426761735977_dp; coef ( 2 ) = 0.22140165119709196_dp
           alpha( 3 ) = 0.03711681387997685_dp;  coef ( 3 ) = 0.5157825082425954_dp
           alpha( 4 ) = 0.057006537132452996_dp; coef ( 4 ) = 0.30947877900909615_dp
           alpha( 5 ) = 0.15940261392199198_dp;  coef ( 5 ) = -0.032377223278152595_dp
           alpha( 6 ) = 0.3168807735151213_dp;   coef ( 6 ) = -0.004680397785899489_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 26 ) ! 7g
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.04636088534182354_dp;  coef ( 1 ) = 1.0000000000000002_dp
         CASE ( 2 )
           alpha( 1 ) = 0.033089231025105_dp;    coef ( 1 ) = 0.5580805992890315_dp
           alpha( 2 ) = 0.06689561975474233_dp;  coef ( 2 ) = 0.5213501458304487_dp
         CASE ( 3 )
           alpha( 1 ) = 0.02536613332751714_dp;  coef ( 1 ) = 0.1958347745511226_dp
           alpha( 2 ) = 0.04392292490201606_dp;  coef ( 2 ) = 0.6069779654592313_dp
           alpha( 3 ) = 0.07935572931927462_dp;  coef ( 3 ) = 0.28358299617643856_dp
         CASE ( 4 )
           alpha( 1 ) = 0.026003867919966516_dp; coef ( 1 ) = 0.2209616450163026_dp
           alpha( 2 ) = 0.045349648935287534_dp; coef ( 2 ) = 0.6131566254139853_dp
           alpha( 3 ) = 0.08207336886154683_dp;  coef ( 3 ) = 0.25321528599323867_dp
           alpha( 4 ) = 0.43630841028097767_dp;  coef ( 4 ) = -0.0028051508587602996_dp
         CASE ( 5 )
           alpha( 1 ) = 0.021295076240560495_dp; coef ( 1 ) = 0.06733674728927509_dp
           alpha( 2 ) = 0.03407297844465803_dp;  coef ( 2 ) = 0.39912419436935137_dp
           alpha( 3 ) = 0.05463343066262556_dp;  coef ( 3 ) = 0.47165433606934565_dp
           alpha( 4 ) = 0.09142718539759649_dp;  coef ( 4 ) = 0.15054429419419305_dp
           alpha( 5 ) = 0.3631735814032239_dp;   coef ( 5 ) = -0.003735801628196915_dp
         CASE ( 6 )
           alpha( 1 ) = 0.02154419398101922_dp;  coef ( 1 ) = 0.07294053032928545_dp
           alpha( 2 ) = 0.03458506500828156_dp;  coef ( 2 ) = 0.41255126079553106_dp
           alpha( 3 ) = 0.055591221080949486_dp; coef ( 3 ) = 0.4632159545572935_dp
           alpha( 4 ) = 0.0929128679746364_dp;   coef ( 4 ) = 0.14023937661406255_dp
           alpha( 5 ) = 0.3353651542417124_dp;   coef ( 5 ) = -0.00381264707490568_dp
           alpha( 6 ) = 0.8332537115987487_dp;   coef ( 6 ) = -0.0002777669788215993_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 27 ) ! 7h
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.055067144574685124_dp; coef ( 1 ) = 1.0000000000000002_dp
         CASE ( 2 )
           alpha( 1 ) = 0.039504306298949056_dp; coef ( 1 ) = 0.6244224977499361_dp
           alpha( 2 ) = 0.08748752023471955_dp;  coef ( 2 ) = 0.48947010311740946_dp
         CASE ( 3 )
           alpha( 1 ) = 0.031898703524951964_dp; coef ( 1 ) = 0.329315310335007_dp
           alpha( 2 ) = 0.060322239914691864_dp; coef ( 2 ) = 0.6230689660925297_dp
           alpha( 3 ) = 0.12542722211419774_dp;  coef ( 3 ) = 0.1792628471936351_dp
         CASE ( 4 )
           alpha( 1 ) = 0.027133297390849163_dp; coef ( 1 ) = 0.16105237279351875_dp
           alpha( 2 ) = 0.04709911293196238_dp;  coef ( 2 ) = 0.5444503553138231_dp
           alpha( 3 ) = 0.08417053216258696_dp;  coef ( 3 ) = 0.36701928787783233_dp
           alpha( 4 ) = 0.16769024925679932_dp;  coef ( 4 ) = 0.06200372233567246_dp
         CASE ( 5 )
           alpha( 1 ) = 0.0236897013320386_dp;   coef ( 1 ) = 0.07357940603188007_dp
           alpha( 2 ) = 0.0387782492773279_dp;   coef ( 2 ) = 0.39203903619828495_dp
           alpha( 3 ) = 0.0637197820910416_dp;   coef ( 3 ) = 0.46247919691910666_dp
           alpha( 4 ) = 0.10978567271801772_dp;  coef ( 4 ) = 0.1832898128285397_dp
           alpha( 5 ) = 0.2085364358880184_dp;   coef ( 5 ) = 0.023617413014962847_dp
         CASE ( 6 )
           alpha( 1 ) = 0.020820906789523693_dp; coef ( 1 ) = 0.029194538236458435_dp
           alpha( 2 ) = 0.032494198134763036_dp; coef ( 2 ) = 0.23574088940165333_dp
           alpha( 3 ) = 0.05023584266268183_dp;  coef ( 3 ) = 0.44360291827793874_dp
           alpha( 4 ) = 0.0794251556241889_dp;   coef ( 4 ) = 0.3155567845762703_dp
           alpha( 5 ) = 0.13148712880820534_dp;  coef ( 5 ) = 0.09872391034903429_dp
           alpha( 6 ) = 0.2359872300875659_dp;   coef ( 6 ) = 0.012260115931631366_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE ( 28 ) ! 7i
       SELECT CASE ( n )
         CASE ( 1 )
           alpha( 1 ) = 0.06380473723420682_dp;  coef ( 1 ) = 0.9999999999999998_dp
         CASE ( 2 )
           alpha( 1 ) = 0.045457834184705205_dp; coef ( 1 ) = 0.6463784604860205_dp
           alpha( 2 ) = 0.10559967649617007_dp;  coef ( 2 ) = 0.49629977931368163_dp
         CASE ( 3 )
           alpha( 1 ) = 0.03666777792117483_dp;  coef ( 1 ) = 0.3678119967003228_dp
           alpha( 2 ) = 0.07224640498805254_dp;  coef ( 2 ) = 0.618554233547922_dp
           alpha( 3 ) = 0.15965614371413928_dp;  coef ( 3 ) = 0.18501800416077063_dp
         CASE ( 4 )
           alpha( 1 ) = 0.031314605885208034_dp; coef ( 1 ) = 0.20067272567322478_dp
           alpha( 2 ) = 0.05654362551330421_dp;  coef ( 2 ) = 0.5656501384781_dp
           alpha( 3 ) = 0.1062477818701912_dp;   coef ( 3 ) = 0.3489852946890129_dp
           alpha( 4 ) = 0.230241610829344_dp;    coef ( 4 ) = 0.06217858499289526_dp
         CASE ( 5 )
           alpha( 1 ) = 0.02764326937219174_dp;  coef ( 1 ) = 0.1083794411727515_dp
           alpha( 2 ) = 0.047149043197080166_dp; coef ( 2 ) = 0.4503068342384467_dp
           alpha( 3 ) = 0.08137277011103897_dp;  coef ( 3 ) = 0.4431260931378365_dp
           alpha( 4 ) = 0.1500423393404324_dp;   coef ( 4 ) = 0.15687089041229385_dp
           alpha( 5 ) = 0.32194304635741355_dp;  coef ( 5 ) = 0.020215660073712404_dp
         CASE ( 6 )
           alpha( 1 ) = 0.024937530903194137_dp; coef ( 1 ) = 0.058696755716280226_dp
           alpha( 2 ) = 0.04081596003413946_dp;  coef ( 2 ) = 0.3330255270864295_dp
           alpha( 3 ) = 0.06665805855772095_dp;  coef ( 3 ) = 0.4640568401312629_dp
           alpha( 4 ) = 0.11296341467059927_dp;  coef ( 4 ) = 0.253753255644_dp
           alpha( 5 ) = 0.2062966749989214_dp;   coef ( 5 ) = 0.06316006418978992_dp
           alpha( 6 ) = 0.4401320366389093_dp;   coef ( 6 ) = 0.0065607991581329616_dp
         CASE DEFAULT
           CPABORT("Illegal expansion")
       END SELECT
    CASE DEFAULT
      CPABORT("Illegal STO")
  END SELECT
  ! scaling: a = zeta^2 * a
  alpha (1:n) = zeta**2 * alpha (1:n)

END SUBROUTINE get_sto_ng

END MODULE sto_ng

