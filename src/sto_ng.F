!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/sto_ng *
!!
!!   NAME
!!     sto_ng
!!
!!   AUTHOR
!!     jgh
!!
!!   MODIFICATION HISTORY
!!     06.2004 created 
!!
!***
!****************************************************************************
MODULE sto_ng
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dbl

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='sto_ng'

  PUBLIC :: get_sto_ng

!****************************************************************************

!!****s* sto_ng/get_sto_ng [1.0] *
!!
!!   NAME
!!     get_sto_ng
!!
!!   FUNCTION
!!     return STO-NG parameters;
!!       INPUT: zeta (Slater exponent)
!!              n    (Expansion length)
!!              nq   (principle quantum number)
!!              lq   (angular quantum number)
!!       OUTPUT: alpha (n Gaussian exponents)
!!               coef  (n corresponding expansion coefficients)
!!
!!   AUTHOR
!!     jgh 
!!
!!   MODIFICATION HISTORY
!!     06.2004 created
!!
!***
!****************************************************************************

CONTAINS

!***************************************************************************
SUBROUTINE get_sto_ng ( zeta, n, nq, lq, alpha, coef ) 

  REAL ( dbl ), INTENT ( IN )                 :: zeta
  INTEGER, INTENT ( IN )                      :: n, nq, lq
  REAL ( dbl ), DIMENSION (:), INTENT ( OUT ) :: alpha, coef

  CHARACTER(len=*), PARAMETER :: routineN = 'get_sto_ng', &
      routineP = moduleN//':'//routineN

  INTEGER                                     :: na, nc, m
  LOGICAL                                     :: failure

  TYPE(cp_error_type)                         :: error

!---------------------------------------------------------------------------

  failure=.FALSE.

  na = SIZE ( alpha )
  nc = SIZE ( coef )

  CPPrecondition((n < 1 .OR. n > 6),cp_failure_level,routineP,error,failure)
  CPPrecondition((na < n),cp_failure_level,routineP,error,failure)
  CPPrecondition((nc < n),cp_failure_level,routineP,error,failure)
  CPPrecondition((zeta < 0._dbl),cp_failure_level,routineP,error,failure)
  CPPrecondition((nq < 1),cp_failure_level,routineP,error,failure)
  CPPrecondition((lq < 0),cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
    m = (nq*(nq-1))/2 + lq + 1
    SELECT CASE ( m )
      CASE (  1 ) ! 1s
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 2.709496091e-1_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 8.518186635e-1_dbl; coef ( 1 ) = 4.301284983e-1_dbl
             alpha( 2 ) = 1.516232927e-1_dbl; coef ( 2 ) = 6.789135305e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 2.227660584e+0_dbl; coef ( 1 ) = 1.543289673e-1_dbl
             alpha( 2 ) = 4.057711562e-1_dbl; coef ( 2 ) = 5.353281423e-1_dbl
             alpha( 3 ) = 1.098175104e-1_dbl; coef ( 3 ) = 4.446345422e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 5.216844534e+0_dbl; coef ( 1 ) = 5.675242080e-2_dbl
             alpha( 2 ) = 9.546182760e-1_dbl; coef ( 2 ) = 2.601413550e-1_dbl
             alpha( 3 ) = 2.652034102e-1_dbl; coef ( 3 ) = 5.328461143e-1_dbl
             alpha( 4 ) = 8.801862774e-2_dbl; coef ( 4 ) = 2.916254405e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 1.130563696e+1_dbl; coef ( 1 ) = 2.214055312e-2_dbl
             alpha( 2 ) = 2.071728178e+0_dbl; coef ( 2 ) = 1.135411520e-1_dbl
             alpha( 3 ) = 5.786484833e-1_dbl; coef ( 3 ) = 3.318161484e-1_dbl
             alpha( 4 ) = 1.975724573e-1_dbl; coef ( 4 ) = 4.825700713e-1_dbl
             alpha( 5 ) = 7.445271746e-2_dbl; coef ( 5 ) = 1.935721966e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 2.310303149e+1_dbl; coef ( 1 ) = 9.163596280e-3_dbl
             alpha( 2 ) = 4.235915534e+0_dbl; coef ( 2 ) = 4.936149294e-2_dbl
             alpha( 3 ) = 1.185056519e+0_dbl; coef ( 3 ) = 1.685383049e-1_dbl
             alpha( 4 ) = 4.070988982e-1_dbl; coef ( 4 ) = 3.705627997e-1_dbl
             alpha( 5 ) = 1.580884151e-1_dbl; coef ( 5 ) = 4.164915298e-1_dbl
             alpha( 6 ) = 6.510953954e-2_dbl; coef ( 6 ) = 1.303340841e-1_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE (  2 ) ! 2s
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 1.012151084e-1_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 1.292278611e-1_dbl; coef ( 1 ) = 7.470867124e-1_dbl
             alpha( 2 ) = 4.908584205e-2_dbl; coef ( 2 ) = 2.855980556e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 2.581578398e+0_dbl; coef ( 1 ) =-5.994474934e-2_dbl
             alpha( 2 ) = 1.567622104e-1_dbl; coef ( 2 ) = 5.960385398e-1_dbl
             alpha( 3 ) = 6.018332272e-2_dbl; coef ( 3 ) = 4.581786291e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 1.161525551e+1_dbl; coef ( 1 ) =-1.198411747e-2_dbl
             alpha( 2 ) = 2.000243111e+0_dbl; coef ( 2 ) =-5.472052539e-2_dbl
             alpha( 3 ) = 1.607280687e-1_dbl; coef ( 3 ) = 5.805587176e-1_dbl
             alpha( 4 ) = 6.125744532e-2_dbl; coef ( 4 ) = 4.770079976e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 8.984956862e+0_dbl; coef ( 1 ) =-1.596349096e-2_dbl
             alpha( 2 ) = 1.673710636e+0_dbl; coef ( 2 ) =-5.685884883e-2_dbl
             alpha( 3 ) = 1.944726668e-1_dbl; coef ( 3 ) = 3.698265599e-1_dbl
             alpha( 4 ) = 8.806345634e-2_dbl; coef ( 4 ) = 5.480512593e-1_dbl
             alpha( 5 ) = 4.249068522e-2_dbl; coef ( 5 ) = 1.472634893e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 2.768496241e+1_dbl; coef ( 1 ) =-4.151277819e-3_dbl
             alpha( 2 ) = 5.077140627e+0_dbl; coef ( 2 ) =-2.067024148e-2_dbl
             alpha( 3 ) = 1.426786050e+0_dbl; coef ( 3 ) =-5.150303337e-2_dbl
             alpha( 4 ) = 2.040335729e-1_dbl; coef ( 4 ) = 3.346271174e-1_dbl
             alpha( 5 ) = 9.260298399e-2_dbl; coef ( 5 ) = 5.621061301e-1_dbl
             alpha( 6 ) = 4.416183978e-2_dbl; coef ( 6 ) = 1.712994697e-1_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE (  3 ) ! 2p
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 1.759666885e-1_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 4.323908358e-1_dbl; coef ( 1 ) = 4.522627513e-1_dbl
             alpha( 2 ) = 1.069439065e-1_dbl; coef ( 2 ) = 6.713122642e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 9.192379002e-1_dbl; coef ( 1 ) = 1.623948553e-1_dbl
             alpha( 2 ) = 2.359194503e-1_dbl; coef ( 2 ) = 5.661708862e-1_dbl
             alpha( 3 ) = 8.009805746e-2_dbl; coef ( 3 ) = 4.223071752e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 1.798260992e+0_dbl; coef ( 1 ) = 5.713170255e-2_dbl
             alpha( 2 ) = 4.662622228e-1_dbl; coef ( 2 ) = 2.857455515e-1_dbl
             alpha( 3 ) = 1.643718620e-1_dbl; coef ( 3 ) = 5.517873105e-1_dbl
             alpha( 4 ) = 6.543927065e-2_dbl; coef ( 4 ) = 2.632314924e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 3.320386533e+0_dbl; coef ( 1 ) = 2.079051117e-2_dbl
             alpha( 2 ) = 8.643257633e-1_dbl; coef ( 2 ) = 1.235472099e-1_dbl
             alpha( 3 ) = 3.079819284e-1_dbl; coef ( 3 ) = 3.667738886e-1_dbl
             alpha( 4 ) = 1.273309895e-1_dbl; coef ( 4 ) = 4.834930290e-1_dbl
             alpha( 5 ) = 5.606243164e-2_dbl; coef ( 5 ) = 1.653444074e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 5.868285913e+0_dbl; coef ( 1 ) = 7.924233646e-3_dbl
             alpha( 2 ) = 1.530329631e+0_dbl; coef ( 2 ) = 5.144104825e-2_dbl
             alpha( 3 ) = 5.475665231e-1_dbl; coef ( 3 ) = 1.898400060e-1_dbl
             alpha( 4 ) = 2.288932733e-1_dbl; coef ( 4 ) = 4.049863191e-1_dbl
             alpha( 5 ) = 1.046655969e-1_dbl; coef ( 5 ) = 4.012362861e-1_dbl
             alpha( 6 ) = 4.948220127e-2_dbl; coef ( 6 ) = 1.051855189e-1_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE (  4 ) ! 3s
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 5.296881757e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 6.694095822e-1_dbl; coef ( 1 ) =-1.529645716e-1_dbl
             alpha( 2 ) = 5.837135094e-2_dbl; coef ( 2 ) = 1.051370110e+0_dbl
           CASE ( 3 )
             alpha( 1 ) = 5.641487709e-1_dbl; coef ( 1 ) =-1.782577972e-1_dbl
             alpha( 2 ) = 6.924421391e-2_dbl; coef ( 2 ) = 8.612761663e-1_dbl
             alpha( 3 ) = 3.269529097e-2_dbl; coef ( 3 ) = 2.261841969e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 1.513265591e+0_dbl; coef ( 1 ) =-3.295496352e-2_dbl
             alpha( 2 ) = 4.262497508e-1_dbl; coef ( 2 ) =-1.724516959e-1_dbl
             alpha( 3 ) = 7.643320863e-2_dbl; coef ( 3 ) = 7.518511194e-1_dbl
             alpha( 4 ) = 3.760545063e-2_dbl; coef ( 4 ) = 3.589627310e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 4.275877914e+0_dbl; coef ( 1 ) =-3.920358850e-3_dbl
             alpha( 2 ) = 1.132409433e+0_dbl; coef ( 2 ) =-4.168430506e-2_dbl
             alpha( 3 ) = 4.016256968e-1_dbl; coef ( 3 ) =-1.637440990e-1_dbl
             alpha( 4 ) = 7.732370620e-2_dbl; coef ( 4 ) = 7.419373723e-1_dbl
             alpha( 5 ) = 3.800708627e-2_dbl; coef ( 5 ) = 3.724364929e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 3.273031938e+0_dbl; coef ( 1 ) =-6.775596947e-3_dbl
             alpha( 2 ) = 9.200611311e-1_dbl; coef ( 2 ) =-5.639325779e-2_dbl
             alpha( 3 ) = 3.593349765e-1_dbl; coef ( 3 ) =-1.587856086e-1_dbl
             alpha( 4 ) = 8.636686991e-2_dbl; coef ( 4 ) = 5.534527651e-1_dbl
             alpha( 5 ) = 4.797373812e-2_dbl; coef ( 5 ) = 5.015351020e-1_dbl
             alpha( 6 ) = 2.724741144e-2_dbl; coef ( 6 ) = 7.223633674e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE (  5 ) ! 3p
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 9.113614253e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 1.458620964e-1_dbl; coef ( 1 ) = 5.349653144e-1_dbl
             alpha( 2 ) = 5.664210742e-2_dbl; coef ( 2 ) = 5.299607212e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 2.692880368e+0_dbl; coef ( 1 ) =-1.061945788e-2_dbl
             alpha( 2 ) = 1.489359592e-1_dbl; coef ( 2 ) = 5.218564264e-1_dbl
             alpha( 3 ) = 5.739585040e-2_dbl; coef ( 3 ) = 5.450015143e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 1.853180239e+0_dbl; coef ( 1 ) =-1.434249391e-2_dbl
             alpha( 2 ) = 1.915075719e-1_dbl; coef ( 2 ) = 2.755177589e-1_dbl
             alpha( 3 ) = 8.655487938e-2_dbl; coef ( 3 ) = 5.846750879e-1_dbl
             alpha( 4 ) = 4.184253862e-2_dbl; coef ( 4 ) = 2.144986514e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 6.466803859e+0_dbl; coef ( 1 ) =-2.329023747e-3_dbl
             alpha( 2 ) = 1.555914802e+0_dbl; coef ( 2 ) =-1.357395221e-2_dbl
             alpha( 3 ) = 1.955925255e-1_dbl; coef ( 3 ) = 2.632185383e-1_dbl
             alpha( 4 ) = 8.809647701e-2_dbl; coef ( 4 ) = 5.880427024e-1_dbl
             alpha( 5 ) = 4.234835707e-2_dbl; coef ( 5 ) = 2.242794445e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 5.077973607e+0_dbl; coef ( 1 ) =-3.329929840e-3_dbl
             alpha( 2 ) = 1.340786940e+0_dbl; coef ( 2 ) =-1.419488340e-2_dbl
             alpha( 3 ) = 2.248434849e-1_dbl; coef ( 3 ) = 1.639395770e-1_dbl
             alpha( 4 ) = 1.131741848e-1_dbl; coef ( 4 ) = 4.485358256e-1_dbl
             alpha( 5 ) = 6.076408893e-2_dbl; coef ( 5 ) = 3.908813050e-1_dbl
             alpha( 6 ) = 3.315424265e-2_dbl; coef ( 6 ) = 7.411456232e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE (  6 ) ! 3d
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 1.302270363e-1_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 2.777427345e-1_dbl; coef ( 1 ) = 4.666137923e-1_dbl
             alpha( 2 ) = 8.336507714e-2_dbl; coef ( 2 ) = 6.644706516e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 5.229112225e-1_dbl; coef ( 1 ) = 1.696596060e-1_dbl
             alpha( 2 ) = 1.639595876e-1_dbl; coef ( 2 ) = 5.847984817e-1_dbl
             alpha( 3 ) = 6.386630021e-2_dbl; coef ( 3 ) = 4.056779523e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 9.185846715e-1_dbl; coef ( 1 ) = 5.799057705e-2_dbl
             alpha( 2 ) = 2.920461109e-1_dbl; coef ( 2 ) = 3.045581349e-1_dbl
             alpha( 3 ) = 1.187568890e-1_dbl; coef ( 3 ) = 5.601358038e-1_dbl
             alpha( 4 ) = 5.286755896e-2_dbl; coef ( 4 ) = 2.432423313e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 1.539033958e+0_dbl; coef ( 1 ) = 2.020869128e-2_dbl
             alpha( 2 ) = 4.922090297e-1_dbl; coef ( 2 ) = 1.321157923e-1_dbl
             alpha( 3 ) = 2.029756928e-1_dbl; coef ( 3 ) = 3.911240346e-1_dbl
             alpha( 4 ) = 9.424112917e-2_dbl; coef ( 4 ) = 4.779609701e-1_dbl
             alpha( 5 ) = 4.569058269e-2_dbl; coef ( 5 ) = 1.463662294e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 2.488296923e+0_dbl; coef ( 1 ) = 7.283828112e-3_dbl
             alpha( 2 ) = 7.981487853e-1_dbl; coef ( 2 ) = 5.386799363e-2_dbl
             alpha( 3 ) = 3.311327490e-1_dbl; coef ( 3 ) = 2.072139149e-1_dbl
             alpha( 4 ) = 1.559114463e-1_dbl; coef ( 4 ) = 4.266269092e-1_dbl
             alpha( 5 ) = 7.877734732e-2_dbl; coef ( 5 ) = 3.8431000204-1_dbl
             alpha( 6 ) = 4.058484363e-2_dbl; coef ( 6 ) = 8.902827546e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE (  7 ) ! 4s
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 3.264600274e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 2.441785453e-1_dbl; coef ( 1 ) =-3.046656896e-1_dbl
             alpha( 2 ) = 4.051097664e-2_dbl; coef ( 2 ) = 1.146877294e+0_dbl
           CASE ( 3 )
             alpha( 1 ) = 2.267938753e-1_dbl; coef ( 1 ) =-3.349048323e-1_dbl
             alpha( 2 ) = 4.448178019e-2_dbl; coef ( 2 ) = 1.056744667e+0_dbl
             alpha( 3 ) = 2.195294664e-2_dbl; coef ( 3 ) = 1.256661680e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 3.242212833e-1_dbl; coef ( 1 ) =-1.120682822e-1_dbl
             alpha( 2 ) = 1.663217177e-1_dbl; coef ( 2 ) =-2.845426863e-1_dbl
             alpha( 3 ) = 5.081097451e-2_dbl; coef ( 3 ) = 8.909873788e-1_dbl
             alpha( 4 ) = 2.829066600e-2_dbl; coef ( 4 ) = 3.517811205e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 2.980263783e+0_dbl; coef ( 1 ) = 1.513948997e-3_dbl
             alpha( 2 ) = 3.792228833e-1_dbl; coef ( 2 ) =-7.316801548e-2_dbl
             alpha( 3 ) = 1.789717224e-1_dbl; coef ( 3 ) =-3.143703799e-1_dbl
             alpha( 4 ) = 5.002110360e-2_dbl; coef ( 4 ) = 9.032615169e-1_dbl
             alpha( 5 ) = 2.789361681e-2_dbl; coef ( 5 ) = 3.294210848e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 3.232838646e+0_dbl; coef ( 1 ) = 1.374817488e-3_dbl
             alpha( 2 ) = 3.605788802e-1_dbl; coef ( 2 ) =-8.666390043e-2_dbl
             alpha( 3 ) = 1.717905487e-1_dbl; coef ( 3 ) =-3.130627309e-1_dbl
             alpha( 4 ) = 5.277666487e-2_dbl; coef ( 4 ) = 7.812787397e-1_dbl
             alpha( 5 ) = 3.163400284e-2_dbl; coef ( 5 ) = 4.389247988e-1_dbl
             alpha( 6 ) = 1.874093091e-2_dbl; coef ( 6 ) = 2.487178756e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE (  8 ) ! 4p
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 5.578350235e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 6.190052680e-2_dbl; coef ( 1 ) = 8.743116767e-1_dbl
             alpha( 2 ) = 2.648418407e-2_dbl; coef ( 2 ) = 1.513640107e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 4.859692220e-1_dbl; coef ( 1 ) =-6.147823411e-2_dbl
             alpha( 2 ) = 7.430216918e-2_dbl; coef ( 2 ) = 6.604172234e-1_dbl
             alpha( 3 ) = 3.653340923e-2_dbl; coef ( 3 ) = 3.932639495e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 1.492607880e+0_dbl; coef ( 1 ) =-6.035216774e-3_dbl
             alpha( 2 ) = 4.327619272e-1_dbl; coef ( 2 ) =-6.013310874e-2_dbl
             alpha( 3 ) = 7.553156064e-2_dbl; coef ( 3 ) = 6.451518200e-1_dbl
             alpha( 4 ) = 3.706272183e-2_dbl; coef ( 4 ) = 4.117923820e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 1.091977298e+0_dbl; coef ( 1 ) =-1.143929558e-2_dbl
             alpha( 2 ) = 3.719985051e-1_dbl; coef ( 2 ) =-6.322651538e-2_dbl
             alpha( 3 ) = 8.590019352e-2_dbl; coef ( 3 ) = 4.398907721e-1_dbl
             alpha( 4 ) = 4.786503860e-2_dbl; coef ( 4 ) = 5.245859166e-1_dbl
             alpha( 5 ) = 2.730479990e-2_dbl; coef ( 5 ) = 1.017072253e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 2.389722618e+0_dbl; coef ( 1 ) =-1.665913575e-3_dbl
             alpha( 2 ) = 7.960947826e-1_dbl; coef ( 2 ) =-1.657464971e-2_dbl
             alpha( 3 ) = 3.415541380e-1_dbl; coef ( 3 ) =-5.958513378e-2_dbl
             alpha( 4 ) = 8.847434525e-2_dbl; coef ( 4 ) = 4.053115554e-1_dbl
             alpha( 5 ) = 4.958248334e-2_dbl; coef ( 5 ) = 5.433958189e-1_dbl
             alpha( 6 ) = 2.816929784e-2_dbl; coef ( 6 ) = 1.204970491e-1_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE (  9 ) ! 4d
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 7.941656339e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 1.330958892e-1_dbl; coef ( 1 ) = 4.932764167e-1_dbl
             alpha( 2 ) = 5.272119659e-2_dbl; coef ( 2 ) = 5.9187277866-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 1.777717219e-1_dbl; coef ( 1 ) = 2.308552718e-1_dbl
             alpha( 2 ) = 8.040647350e-2_dbl; coef ( 2 ) = 6.042409177e-1_dbl
             alpha( 3 ) = 3.949855551e-2_dbl; coef ( 3 ) = 2.595768926e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 1.995825422e+0_dbl; coef ( 1 ) =-2.816702620e-3_dbl
             alpha( 2 ) = 1.823461280e-1_dbl; coef ( 2 ) = 2.177095871e-1_dbl
             alpha( 3 ) = 8.197240896e-2_dbl; coef ( 3 ) = 6.058047348e-1_dbl
             alpha( 4 ) = 4.000634951e-2_dbl; coef ( 4 ) = 2.717811257e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 1.522122079e+0_dbl; coef ( 1 ) =-3.673711876e-3_dbl
             alpha( 2 ) = 2.173041823e-1_dbl; coef ( 2 ) = 1.167122499e-1_dbl
             alpha( 3 ) = 1.084876577e-1_dbl; coef ( 3 ) = 4.216476416e-1_dbl
             alpha( 4 ) = 5.836797641e-2_dbl; coef ( 4 ) = 4.547673415e-1_dbl
             alpha( 5 ) = 3.206682246e-2_dbl; coef ( 5 ) = 1.037803318e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 4.634239420e+0_dbl; coef ( 1 ) =-4.749842876e-4_dbl
             alpha( 2 ) = 1.341648295e+0_dbl; coef ( 2 ) =-3.566777891e-3_dbl
             alpha( 3 ) = 2.209593028e-1_dbl; coef ( 3 ) = 1.108670481e-1_dbl
             alpha( 4 ) = 1.101467943e-1_dbl; coef ( 4 ) = 4.159646930e-1_dbl
             alpha( 5 ) = 5.904190370e-2_dbl; coef ( 5 ) = 4.621672517e-1_dbl
             alpha( 6 ) = 3.232628887e-2_dbl; coef ( 6 ) = 1.081250196e-1_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE ( 10 ) ! 4f
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 1.033434062e-1_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 2.006693538e-1_dbl; coef ( 1 ) = 4.769346276e-1_dbl
             alpha( 2 ) = 6.865384900e-2_dbl; coef ( 2 ) = 6.587383976e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 3.483826963e-1_dbl; coef ( 1 ) = 1.737856685e-1_dbl
             alpha( 2 ) = 1.249380537e-1_dbl; coef ( 2 ) = 5.973380628e-1_dbl
             alpha( 3 ) = 5.349995725e-2_dbl; coef ( 3 ) = 3.929395614e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 5.691670217e-1_dbl; coef ( 1 ) = 5.902730589e-2_dbl
             alpha( 2 ) = 2.074585819e-1_dbl; coef ( 2 ) = 3.191828952e-1_dbl
             alpha( 3 ) = 9.298346885e-2_dbl; coef ( 3 ) = 5.639423893e-1_dbl
             alpha( 4 ) = 4.473508853e-2_dbl; coef ( 4 ) = 2.284796537e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 8.925960415e-1_dbl; coef ( 1 ) = 1.999839052e-2_dbl
             alpha( 2 ) = 3.277589120e-1_dbl; coef ( 2 ) = 1.395427440e-1_dbl
             alpha( 3 ) = 1.492869962e-1_dbl; coef ( 3 ) = 4.091508237e-1_dbl
             alpha( 4 ) = 7.506099109e-2_dbl; coef ( 4 ) = 4.708252119e-1_dbl
             alpha( 5 ) = 3.892475795e-2_dbl; coef ( 5 ) = 1.328082566e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 1.357718039e+0_dbl; coef ( 1 ) = 6.930234381e-3_dbl
             alpha( 2 ) = 5.004907278e-1_dbl; coef ( 2 ) = 5.634263745e-2_dbl
             alpha( 3 ) = 2.296565064e-1_dbl; coef ( 3 ) = 2.217065797e-1_dbl
             alpha( 4 ) = 1.173146814e-1_dbl; coef ( 4 ) = 4.411388883e-1_dbl
             alpha( 5 ) = 6.350097171e-2_dbl; coef ( 5 ) = 3.688112625e-1_dbl
             alpha( 6 ) = 3.474556673e-2_dbl; coef ( 6 ) = 7.787514504e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE ( 11 ) ! 5s
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 2.216912938e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 1.213425654e-1_dbl; coef ( 1 ) =-5.114756049e-1_dbl
             alpha( 2 ) = 3.133152144e-2_dbl; coef ( 2 ) = 1.307377277e+0_dbl
           CASE ( 3 )
             alpha( 1 ) = 1.080198458e-1_dbl; coef ( 1 ) =-6.617401158e-1_dbl
             alpha( 2 ) = 4.408119328e-2_dbl; coef ( 2 ) = 7.467595004e-1_dbl
             alpha( 3 ) = 2.610811810e-2_dbl; coef ( 3 ) = 7.146490945e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 8.602284252e-1_dbl; coef ( 1 ) = 1.103657561e-2_dbl
             alpha( 2 ) = 1.189050200e-1_dbl; coef ( 2 ) =-5.606519023e-1_dbl
             alpha( 3 ) = 3.446076176e-2_dbl; coef ( 3 ) = 1.179429987e+0_dbl
             alpha( 4 ) = 1.974798796e-2_dbl; coef ( 4 ) = 1.734974376e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 7.403763257e-1_dbl; coef ( 1 ) = 1.375523371e-2_dbl
             alpha( 2 ) = 1.367990863e-1_dbl; coef ( 2 ) =-3.097344179e-1_dbl
             alpha( 3 ) = 9.135301779e-2_dbl; coef ( 3 ) =-3.199192259e-1_dbl
             alpha( 4 ) = 3.726907315e-2_dbl; coef ( 4 ) = 1.084547038e+0_dbl
             alpha( 5 ) = 2.241490836e-2_dbl; coef ( 5 ) = 3.345288361e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 1.410128298e+0_dbl; coef ( 1 ) = 2.695439582e-3_dbl
             alpha( 2 ) = 5.077878915e-1_dbl; coef ( 2 ) = 1.850157487e-2_dbl
             alpha( 3 ) = 1.847926858e-1_dbl; coef ( 3 ) =-9.588628125e-2_dbl
             alpha( 4 ) = 1.061070594e-1_dbl; coef ( 4 ) =-5.200673560e-1_dbl
             alpha( 5 ) = 3.669584901e-2_dbl; coef ( 5 ) = 1.087619490e+0_dbl
             alpha( 6 ) = 2.213558430e-2_dbl; coef ( 6 ) = 3.103964343e-1_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE ( 12 ) ! 5p
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 3.769845216e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 2.691294191e-1_dbl; coef ( 1 ) =-1.034227010e-1_dbl
             alpha( 2 ) = 3.980805011e-2_dbl; coef ( 2 ) = 1.033376378e+0_dbl
           CASE ( 3 )
             alpha( 1 ) = 2.127482317e-1_dbl; coef ( 1 ) =-1.389529695e-1_dbl
             alpha( 2 ) = 4.729648620e-2_dbl; coef ( 2 ) = 8.076691064e-1_dbl
             alpha( 3 ) = 2.604865324e-2_dbl; coef ( 3 ) = 2.726029342e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 3.962838833e-1_dbl; coef ( 1 ) =-1.801459207e-2_dbl
             alpha( 2 ) = 1.838858552e-1_dbl; coef ( 2 ) =-1.360777372e-1_dbl
             alpha( 3 ) = 4.943555157e-2_dbl; coef ( 3 ) = 7.533973719e-1_dbl
             alpha( 4 ) = 2.750222273e-2_dbl; coef ( 4 ) = 3.409304859e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 3.422168934e-1_dbl; coef ( 1 ) =-3.113958289e-2_dbl
             alpha( 2 ) = 1.665099900e-1_dbl; coef ( 2 ) =-1.374007017e-1_dbl
             alpha( 3 ) = 5.443732013e-2_dbl; coef ( 3 ) = 5.573881018e-1_dbl
             alpha( 4 ) = 3.367775277e-2_dbl; coef ( 4 ) = 4.855428100e-1_dbl
             alpha( 5 ) = 2.091949042e-2_dbl; coef ( 5 ) = 6.605423564e-2_dbl
           CASE ( 6 )
             alpha( 1 ) = 3.778623374e+0_dbl; coef ( 1 ) = 1.163246387e-4_dbl
             alpha( 2 ) = 3.499121109e-1_dbl; coef ( 2 ) =-2.920771322e-2_dbl
             alpha( 3 ) = 1.683175469e-1_dbl; coef ( 3 ) =-1.381051233e-1_dbl
             alpha( 4 ) = 5.404070736e-2_dbl; coef ( 4 ) = 5.706134877e-1_dbl
             alpha( 5 ) = 3.328911801e-2_dbl; coef ( 5 ) = 4.768808140e-1_dbl
             alpha( 6 ) = 2.063815019e-2_dbl; coef ( 6 ) = 6.021665516e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE ( 13 ) ! 5d
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 5.352200793e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 6.906014388e-2_dbl; coef ( 1 ) = 6.539405185e-1_dbl
             alpha( 2 ) = 3.399457777e-2_dbl; coef ( 2 ) = 3.948945302e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 4.913352950e-1_dbl; coef ( 1 ) =-2.010175008e-2_dbl
             alpha( 2 ) = 7.329090601e-2_dbl; coef ( 2 ) = 5.899370608e-1_dbl
             alpha( 3 ) = 3.594209290e-2_dbl; coef ( 3 ) = 4.658445960e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 4.230617826e-1_dbl; coef ( 1 ) =-2.421626009e-2_dbl
             alpha( 2 ) = 8.293863702e-2_dbl; coef ( 2 ) = 3.937644956e-1_dbl
             alpha( 3 ) = 4.590326388e-2_dbl; coef ( 3 ) = 5.489520286e-1_dbl
             alpha( 4 ) = 2.628744797e-2_dbl; coef ( 4 ) = 1.190436963e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 9.702946470e-1_dbl; coef ( 1 ) =-3.231527611e-3_dbl
             alpha( 2 ) = 3.603270196e-1_dbl; coef ( 2 ) =-2.434931372e-2_dbl
             alpha( 3 ) = 8.668717752e-2_dbl; coef ( 3 ) = 3.440817054e-1_dbl
             alpha( 4 ) = 4.833708379e-2_dbl; coef ( 4 ) = 5.693674376e-1_dbl
             alpha( 5 ) = 2.751899341e-2_dbl; coef ( 5 ) = 1.511340183e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 8.820520428e-1_dbl; coef ( 1 ) =-4.097377019e-3_dbl
             alpha( 2 ) = 3.410838409e-1_dbl; coef ( 2 ) =-2.508271857e-2_dbl
             alpha( 3 ) = 9.204308840e-2_dbl; coef ( 3 ) = 2.648458555e-1_dbl
             alpha( 4 ) = 5.472831774e-2_dbl; coef ( 4 ) = 5.097437054e-1_dbl
             alpha( 5 ) = 3.391202830e-2_dbl; coef ( 5 ) = 2.654483467e-1_dbl
             alpha( 6 ) = 2.108227374e-2_dbl; coef ( 6 ) = 2.623132212e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE ( 14 ) ! 5f
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 6.952785407e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 1.156094555e-1_dbl; coef ( 1 ) = 4.856637346e-1_dbl
             alpha( 2 ) = 4.778940916e-2_dbl; coef ( 2 ) = 6.125980914e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 1.649233885e-1_dbl; coef ( 1 ) = 1.909729355e-1_dbl
             alpha( 2 ) = 7.487066646e-2_dbl; coef ( 2 ) = 6.146060459e-1_dbl
             alpha( 3 ) = 3.735787219e-2_dbl; coef ( 3 ) = 3.059611271e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 2.017831152e-1_dbl; coef ( 1 ) = 9.174268830e-2_dbl
             alpha( 2 ) = 1.001952178e-1_dbl; coef ( 2 ) = 4.023496947e-1_dbl
             alpha( 3 ) = 5.447006630e-2_dbl; coef ( 3 ) = 4.937432100e-1_dbl
             alpha( 4 ) = 3.037569283e-2_dbl; coef ( 4 ) = 1.254001522e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 1.670735676e+0_dbl; coef ( 1 ) =-7.301193568e-4_dbl
             alpha( 2 ) = 2.072477219e-1_dbl; coef ( 2 ) = 8.414991343e-2_dbl
             alpha( 3 ) = 1.024709357e-1_dbl; coef ( 3 ) = 3.923683153e-1_dbl
             alpha( 4 ) = 5.537913898e-2_dbl; coef ( 4 ) = 5.040033146e-1_dbl
             alpha( 5 ) = 3.072866652e-2_dbl; coef ( 5 ) = 1.328979300e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 1.334096840e+0_dbl; coef ( 1 ) =-9.486751531e-4_dbl
             alpha( 2 ) = 2.372312347e-1_dbl; coef ( 2 ) = 4.624275998e-2_dbl
             alpha( 3 ) = 1.269485744e-1_dbl; coef ( 3 ) = 2.373699784e-1_dbl
             alpha( 4 ) = 7.290318381e-2_dbl; coef ( 4 ) = 4.589112231e-1_dbl
             alpha( 5 ) = 4.351355997e-2_dbl; coef ( 5 ) = 3.205010548e-1_dbl
             alpha( 6 ) = 2.598071843e-2_dbl; coef ( 6 ) = 5.077063693e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE ( 15 ) ! 5g
         SELECT CASE ( n )
           CASE ( 1 )
             alpha( 1 ) = 8.565417784e-2_dbl; coef ( 1 ) = 1.000000000e+0_dbl
           CASE ( 2 )
             alpha( 1 ) = 1.554531559e-1_dbl; coef ( 1 ) = 4.848298074e-1_dbl
             alpha( 2 ) = 5.854079811e-2_dbl; coef ( 2 ) = 6.539381621e-1_dbl
           CASE ( 3 )
             alpha( 1 ) = 2.545432122e-1_dbl; coef ( 1 ) = 1.780980905e-1_dbl
             alpha( 2 ) = 1.006544376e-1_dbl; coef ( 2 ) = 6.063757846e-1_dbl
             alpha( 3 ) = 4.624463922e-2_dbl; coef ( 3 ) = 3.828552923e-1_dbl
           CASE ( 4 )
             alpha( 1 ) = 3.945205573e-1_dbl; coef ( 1 ) = 6.010484250e-2_dbl
             alpha( 2 ) = 1.588100623e-1_dbl; coef ( 2 ) = 3.309738329e-1_dbl
             alpha( 3 ) = 7.646521729e-2_dbl; coef ( 3 ) = 5.655207585e-1_dbl
             alpha( 4 ) = 3.898703611e-2_dbl; coef ( 4 ) = 2.171122608e-1_dbl
           CASE ( 5 )
             alpha( 1 ) = 5.895429375e-1_dbl; coef ( 1 ) = 1.998085812e-2_dbl
             alpha( 2 ) = 2.393343780e-1_dbl; coef ( 2 ) = 1.460384050e-1_dbl
             alpha( 3 ) = 1.172646904e-1_dbl; coef ( 3 ) = 4.230565459e-1_dbl
             alpha( 4 ) = 6.254074479e-2_dbl; coef ( 4 ) = 4.635699665e-1_dbl
             alpha( 5 ) = 3.411243214e-2_dbl; coef ( 5 ) = 1.226411691e-1_dbl
           CASE ( 6 )
             alpha( 1 ) = 8.574668996e-1_dbl; coef ( 1 ) = 6.729778096e-3_dbl
             alpha( 2 ) = 3.497184772e-1_dbl; coef ( 2 ) = 5.874145170e-2_dbl
             alpha( 3 ) = 1.727917060e-1_dbl; coef ( 3 ) = 2.339955227e-1_dbl
             alpha( 4 ) = 9.373643151e-2_dbl; coef ( 4 ) = 4.512983737e-1_dbl
             alpha( 5 ) = 5.340032759e-2_dbl; coef ( 5 ) = 3.552053926e-1_dbl
             alpha( 6 ) = 3.057364464e-2_dbl; coef ( 6 ) = 6.974153145e-2_dbl
           CASE DEFAULT
             CPErrorMessage(-100,routineP,"Illegal expansion",error)
         END SELECT
      CASE DEFAULT
        CPErrorMessage(-101,routineP,"Illegal STO",error)
    END SELECT
    ! scaling: a = zeta^2 * a
    alpha (1:n) = zeta**2 * alpha (1:n)
  END IF

END SUBROUTINE get_sto_ng 

!***************************************************************************

END MODULE sto_ng

!***************************************************************************
