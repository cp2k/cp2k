!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/parallel_rng_types [1.0] *
!!
!!   NAME
!!     parallel_rng_types
!!
!!   FUNCTION
!!     Parallel (pseudo)random number generator (RNG) for multiple streams
!!     and substreams of random numbers.
!!
!!     In detail, this RNG provides 2**64 random number streams each with a
!!     length of 2**127 resulting in a length of 2**191 for the total RNG.
!!     Moreover, each stream is divided in 2*51 substream of length 2**76.
!!     The stream lengths refer to the default precision of 32 bit random
!!     number, but also an extended precision of 53 bit per random number
!!     can be requested. In this case, two 32 bit random numbers are used
!!     to generate a 53 bit random number and therefore the stream length
!!     is halved when extended precision are requested.
!!
!!   LITERATURE
!!     P. L'Ecuyer, R. Simard, E. J. Chen, and W. D. Kelton,
!!     "An object-oriented random-number package with many long streams
!!      and substreams", Operations Research 50(6), 1073-1075 (2002)
!!
!!   AUTHOR
!!     C++ code converted to Fortran 90/95 (18.05.2005, Matthias Krack)
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE parallel_rng_types

  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE

  ! Global parameters in this module

  CHARACTER(LEN=*), PARAMETER :: module_name = "parallel_rng_types"

  REAL(KIND=dp), PARAMETER :: norm  = 2.328306549295727688e-10_dp,&
                              m1    = 4294967087.0_dp,&
                              m2    = 4294944443.0_dp,&
                              a12   = 1403580.0_dp,&
                              a13n  = 810728.0_dp,&
                              a21   = 527612.0_dp,&
                              a23n  = 1370589.0_dp,&
                              two17 = 131072.0_dp,&           ! 2**17
                              two53 = 9007199254740992.0_dp,& ! 2**53
                              fact  = 5.9604644775390625e-8   ! 1/2**24

  ! Data type definitions

  ! Information on a stream. The arrays bg, cg, and ig contain the current
  ! state of the stream, the starting state of the current substream, and the
  ! starting state of the stream. This stream generates antithetic variates
  ! if antithetic = .TRUE.. It also generates numbers with extended precision
  ! (53 bits, if machine follows IEEE 754 standard), if extended_precision =
  ! .TRUE., otherwise, numbers with 32 bits precision are generated.

  TYPE rng_stream_type
    PRIVATE
    CHARACTER(LEN=40)             :: name
    REAL(KIND=dp), DIMENSION(3,2) :: bg,cg,ig
    LOGICAL                       :: antithetic,extended_precision
  END TYPE rng_stream_type

  TYPE rng_stream_p_type
    PRIVATE
    TYPE(rng_stream_type), POINTER :: rng_stream
  END TYPE rng_stream_p_type

  ! The following are the transition matrices of the two MRG components
  ! (in matrix form), raised to the powers -1, 1, 2**76, and 2**127, resp.

  REAL(KIND=dp), DIMENSION(3,3) :: a1p0,a1p76,a1p127,&
                                   a2p0,a2p76,a2p127,&
                                   inv_a1,inv_a2

  ! Interfaces

  INTERFACE next_random_number
    MODULE PROCEDURE next_integer_random_number,&
                     next_real_random_number
  END INTERFACE

  ! Public data types

  PUBLIC :: rng_stream_type

  ! Public subroutines

  PUBLIC :: advance_rng_state,&
            create_rng_stream,&
            delete_rng_stream,&
            dump_rng_state,&
            get_rng_stream,&
            init_rng,&
            next_random_number,&
            next_rng_seed,&
            reset_rng_stream,&
            reset_rng_substream,&
            reset_to_next_rng_substream,&
            set_rng_stream,&
            write_rng_matrices,&
            write_rng_stream

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE advance_rng_state(rng_stream,e,c)

    ! Advance the state by n steps, i.e. jump n steps forward, if n > 0, or
    ! backward if n < 0.
    ! IF e > 0, let n = 2**e + c
    ! IF e < 0, let n = -2**(-e) + c
    ! IF e = 0, let n = c

    ! NOTE: The use of this method is discouraged.

    TYPE(rng_stream_type), POINTER :: rng_stream
    INTEGER, INTENT(IN)            :: e,c

    ! Local variables

    REAL(KIND=dp), DIMENSION(3,3) :: u1,u2,v1,v2,w1,w2
    REAL(KIND=dp), DIMENSION(3,2) :: x

    ! -------------------------------------------------------------------------

    u1 = 0.0_dp
    u2 = 0.0_dp
    v1 = 0.0_dp
    v2 = 0.0_dp
    w1 = 0.0_dp
    w2 = 0.0_dp

    IF (e > 0) THEN
      CALL mat_two_pow_mod_m(a1p0,u1,m1,e)
      CALL mat_two_pow_mod_m(a2p0,u2,m2,e)
    ELSE IF (e < 0) THEN
      CALL mat_two_pow_mod_m(inv_a1,u1,m1,-e)
      CALL mat_two_pow_mod_m(inv_a2,u2,m2,-e)
    END IF

    IF (c >= 0) THEN
      CALL mat_pow_mod_m(a1p0,v1,m1,c)
      CALL mat_pow_mod_m(a2p0,v2,m2,c)
    ELSE
      CALL mat_pow_mod_m(inv_a1,v1,m1,-c)
      CALL mat_pow_mod_m(inv_a2,v2,m2,-c)
    END IF

    IF (e == 0) THEN
      w1 = v1
      w2 = v2
    ELSE
      CALL mat_mat_mod_m(u1,v1,w1,m1)
      CALL mat_mat_mod_m(u2,v2,w2,m2)
    END IF

    x = 0.0_dp

    CALL mat_vec_mod_m(w1,rng_stream%cg(:,1),x(:,1),m1)
    CALL mat_vec_mod_m(w2,rng_stream%cg(:,2),x(:,2),m2)

    rng_stream%cg = x

  END SUBROUTINE advance_rng_state

  ! ***************************************************************************

  SUBROUTINE check_seed(seed)

    ! Check that the seeds are legitimate values.

    INTEGER, DIMENSION(3,2), INTENT(IN) :: seed

    ! Local variables

    INTEGER :: i

    ! -------------------------------------------------------------------------

    DO i=1,3
      IF (seed(i,1) >= m1) THEN
        WRITE (*,*) "ERROR: seed(",i,",1) >=",m1
        STOP
      END IF
      IF (seed(i,2) >= m2) THEN
        WRITE (*,*) "ERROR: seed(",i,",2) >=",m2
        STOP
      END IF
    END DO

    IF ((seed(1,1) == 0).AND.&
        (seed(2,1) == 0).AND.&
        (seed(3,1) == 0)) THEN
      WRITE (*,*) "ERROR: First seed = 0"
      STOP
    END IF

    IF ((seed(1,2) == 0).AND.&
        (seed(2,2) == 0).AND.&
        (seed(3,2) == 0)) THEN
      WRITE (*,*) "ERROR: Second seed = 0"
      STOP
    END IF

  END SUBROUTINE check_seed

  ! ***************************************************************************

  SUBROUTINE create_rng_stream(rng_stream,name,last_rng_stream,seed,&
                               antithetic,extended_precision)

    ! Create a new RNG stream.
    ! last_rng_stream is used as a reference stream, if it is specified.

    ! Usage hint:
    ! CALL create_rng_stream(rng_stream,name) to generate the first stream.
    ! Then CALL create_rng_stream(next_rng_stream,name,name,rng_stream) to
    ! create all the following RNG stream.

    TYPE(rng_stream_type), POINTER                :: rng_stream
    CHARACTER(LEN=*), INTENT(IN)                  :: name
    TYPE(rng_stream_type), OPTIONAL, POINTER      :: last_rng_stream
    INTEGER, DIMENSION(3,2), OPTIONAL, INTENT(IN) :: seed
    LOGICAL, OPTIONAL, INTENT(IN)                 :: antithetic,&
                                                     extended_precision

    ! Local variables

    INTEGER :: istat

    ! -------------------------------------------------------------------------

    ALLOCATE (rng_stream,STAT=istat)
    IF (istat /= 0) STOP "Error allocating rng_stream"

    rng_stream%name = name

    IF (PRESENT(seed)) THEN
      CALL check_seed(seed)
      rng_stream%ig = seed
    ELSE IF (PRESENT(last_rng_stream)) THEN
      rng_stream%ig = next_rng_seed(last_rng_stream%ig)
    ELSE
      rng_stream%ig = next_rng_seed()
    END IF

    rng_stream%cg = rng_stream%ig
    rng_stream%bg = rng_stream%ig

    IF (PRESENT(antithetic)) THEN
      rng_stream%antithetic = antithetic
    ELSE IF (PRESENT(last_rng_stream)) THEN
      rng_stream%antithetic = last_rng_stream%antithetic
    ELSE
      rng_stream%antithetic = .FALSE.
    END IF

    IF (PRESENT(extended_precision)) THEN
      rng_stream%extended_precision = extended_precision
    ELSE IF (PRESENT(last_rng_stream)) THEN
      rng_stream%extended_precision = last_rng_stream%extended_precision
    ELSE
      rng_stream%extended_precision = .FALSE.
    END IF

  END SUBROUTINE create_rng_stream

  ! ***************************************************************************

  SUBROUTINE delete_rng_stream(rng_stream)

    ! Delete a random number stream.

    TYPE(rng_stream_type), POINTER :: rng_stream

    ! Local variables

    INTEGER :: istat

    ! -------------------------------------------------------------------------

    DEALLOCATE (rng_stream,STAT=istat)
    IF (istat /= 0) STOP "Error deallocating rng_stream"

  END SUBROUTINE delete_rng_stream

  ! ***************************************************************************

  SUBROUTINE dump_rng_state(rng_stream,output_unit)

    ! Dump the state of a RNG stream to a logical output unit.

    TYPE(rng_stream_type), POINTER :: rng_stream
    INTEGER, INTENT(IN)            :: output_unit

    ! -------------------------------------------------------------------------

    WRITE (UNIT=output_unit,FMT=*)&
      TRIM(rng_stream%name),&
      rng_stream%ig,&
      rng_stream%bg,&
      rng_stream%cg,&
      rng_stream%antithetic,&
      rng_stream%extended_precision

  END SUBROUTINE dump_rng_state

  ! ***************************************************************************

  SUBROUTINE get_rng_stream(rng_stream,name,bg,cg,ig,antithetic,&
                            extended_precision)

    ! Get the components of a RNG stream.

    TYPE(rng_stream_type), POINTER                 :: rng_stream
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT)        :: name
    INTEGER, DIMENSION(3,2), OPTIONAL, INTENT(OUT) :: bg,cg,ig
    LOGICAL, OPTIONAL, INTENT(OUT)                 :: antithetic,&
                                                      extended_precision

    ! -------------------------------------------------------------------------

    IF (PRESENT(name)) name = rng_stream%name
    IF (PRESENT(bg)) bg = rng_stream%bg
    IF (PRESENT(cg)) cg = rng_stream%cg
    IF (PRESENT(ig)) ig = rng_stream%ig
    IF (PRESENT(antithetic)) antithetic = rng_stream%antithetic
    IF (PRESENT(extended_precision)) THEN
      extended_precision = rng_stream%extended_precision
    END IF

  END SUBROUTINE get_rng_stream

  ! ***************************************************************************

  SUBROUTINE init_rng()

    ! Initialize the transformation matrices for the RNG.
    ! This matrices could be stored also as parameters, however, as long as
    ! RESHAPE is not processed properly by forpar.x, we have to initialize
    ! them during the program startup.

    ! -------------------------------------------------------------------------

    ! The following are the transition matrices of the two MRG components
    ! (in matrix form), raised to the powers 1, 2**76, 2**127, and -1

    ! Transition matrix for the first component raised to the power 2**0

    a1p0(1,1) =        0.0_dp
    a1p0(2,1) =        0.0_dp
    a1p0(3,1) =  -810728.0_dp
    a1p0(1,2) =        1.0_dp
    a1p0(2,2) =        0.0_dp
    a1p0(3,2) =  1403580.0_dp
    a1p0(1,3) =        0.0_dp
    a1p0(2,3) =        1.0_dp
    a1p0(3,3) =        0.0_dp

    ! Transition matrix for the second component raised to the power 2**0

    a2p0(1,1) =        0.0_dp
    a2p0(2,1) =        0.0_dp
    a2p0(3,1) = -1370589.0_dp
    a2p0(1,2) =        1.0_dp
    a2p0(2,2) =        0.0_dp
    a2p0(3,2) =        0.0_dp
    a2p0(1,3) =        0.0_dp
    a2p0(2,3) =        1.0_dp
    a2p0(3,3) =   527612.0_dp

    ! Transition matrix for the first component raised to the power 2**76

    a1p76(1,1) =   82758667.0_dp
    a1p76(2,1) = 3672831523.0_dp
    a1p76(3,1) = 3672091415.0_dp
    a1p76(1,2) = 1871391091.0_dp
    a1p76(2,2) =   69195019.0_dp
    a1p76(3,2) = 3528743235.0_dp
    a1p76(1,3) = 4127413238.0_dp
    a1p76(2,3) = 1871391091.0_dp
    a1p76(3,3) =   69195019.0_dp

    ! Transition matrix for the second component raised to the power 2**76

    a2p76(1,1) = 1511326704.0_dp
    a2p76(2,1) = 4292754251.0_dp
    a2p76(3,1) = 3859662829.0_dp
    a2p76(1,2) = 3759209742.0_dp
    a2p76(2,2) = 1511326704.0_dp
    a2p76(3,2) = 4292754251.0_dp
    a2p76(1,3) = 1610795712.0_dp
    a2p76(2,3) = 3889917532.0_dp
    a2p76(3,3) = 3708466080.0_dp

    ! Transition matrix for the first component raised to the power 2**127

    a1p127(1,1) = 2427906178.0_dp
    a1p127(2,1) =  226153695.0_dp
    a1p127(3,1) = 1988835001.0_dp
    a1p127(1,2) = 3580155704.0_dp
    a1p127(2,2) = 1230515664.0_dp
    a1p127(3,2) =  986791581.0_dp
    a1p127(1,3) =  949770784.0_dp
    a1p127(2,3) = 3580155704.0_dp
    a1p127(3,3) = 1230515664.0_dp

    ! Transition matrix for the second component raised to the power 2**127

    a2p127(1,1) = 1464411153.0_dp
    a2p127(2,1) =   32183930.0_dp
    a2p127(3,1) = 2824425944.0_dp
    a2p127(1,2) =  277697599.0_dp
    a2p127(2,2) = 1464411153.0_dp
    a2p127(3,2) =   32183930.0_dp
    a2p127(1,3) = 1610723613.0_dp
    a2p127(2,3) = 1022607788.0_dp
    a2p127(3,3) = 2093834863.0_dp

    ! Inverse of a1p0

    inv_a1(1,1) =  184888585.0_dp
    inv_a1(2,1) =          1.0_dp
    inv_a1(3,1) =          0.0_dp
    inv_a1(1,2) =          0.0_dp
    inv_a1(2,2) =          0.0_dp
    inv_a1(3,2) =          1.0_dp
    inv_a1(1,3) = 1945170933.0_dp
    inv_a1(2,3) =          0.0_dp
    inv_a1(3,3) =          0.0_dp

    ! Inverse of a2p0

    inv_a2(1,1) =          0.0_dp
    inv_a2(2,1) =          1.0_dp
    inv_a2(3,1) =          0.0_dp
    inv_a2(1,2) =  360363334.0_dp
    inv_a2(2,2) =          0.0_dp
    inv_a2(3,2) =          1.0_dp
    inv_a2(1,3) = 4225571728.0_dp
    inv_a2(2,3) =          0.0_dp
    inv_a2(3,3) =          0.0_dp

  END SUBROUTINE init_rng

  ! ***************************************************************************

  SUBROUTINE mat_mat_mod_m(a,b,c,m)

    ! Returns c = MODULO(a*b,m).

    REAL(KIND=dp), DIMENSION(3,3), INTENT(IN)  :: a,b
    REAL(KIND=dp), DIMENSION(3,3), INTENT(OUT) :: c
    REAL(KIND=dp), INTENT(IN)                  :: m

    ! Local variables

    INTEGER :: i

    ! -------------------------------------------------------------------------

    DO i=1,3
      CALL mat_vec_mod_m(a,b(:,i),c(:,i),m)
    END DO

  END SUBROUTINE mat_mat_mod_m

  ! ***************************************************************************

  SUBROUTINE mat_pow_mod_m(a,b,m,n)

    ! Compute matrix b = MODULO(a**n,m)

    REAL(KIND=dp), DIMENSION(3,3), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(3,3), INTENT(OUT) :: b
    REAL(KIND=dp), INTENT(IN)                  :: m
    INTEGER, INTENT(IN)                        :: n

    ! Local variables

    REAL(KIND=dp), DIMENSION(3,3) :: u,v,w
    INTEGER                       :: i

    ! -------------------------------------------------------------------------

    ! Initialize: u = v = a; b = I

    w = a

    b(1,1) = 1.0_dp
    b(2,1) = 0.0_dp
    b(3,1) = 0.0_dp
    b(1,2) = 0.0_dp
    b(2,2) = 1.0_dp
    b(3,2) = 0.0_dp
    b(1,3) = 0.0_dp
    b(2,3) = 0.0_dp
    b(3,3) = 1.0_dp

    ! Compute b = MODULO(a**n,m) using the binary decomposition of n

    i = n

    DO
      IF (MODULO(i,2) /= 0) THEN
        u = w
        v = b
        CALL mat_mat_mod_m(u,v,b,m)
      END IF
      i = i/2
      IF (i == 0) EXIT
      u = w
      v = w
      CALL mat_mat_mod_m(u,v,w,m)
    END DO

  END SUBROUTINE mat_pow_mod_m

  ! ***************************************************************************

  SUBROUTINE mat_two_pow_mod_m(a,b,m,e)

    ! Compute matrix b = MODULO(a**(2**e),m)

    REAL(KIND=dp), DIMENSION(3,3), INTENT(IN)  :: a
    REAL(KIND=dp), DIMENSION(3,3), INTENT(OUT) :: b
    REAL(KIND=dp), INTENT(IN)                  :: m
    INTEGER, INTENT(IN)                        :: e

    ! Local variables

    REAL(KIND=dp), DIMENSION(3,3) :: u,v
    INTEGER                       :: i

    ! -------------------------------------------------------------------------

    b = a

    ! Compute b = MODULO(a**(2**e),m)

    DO i=1,e
      u = b
      v = b
      CALL mat_mat_mod_m(u,v,b,m)
    END DO

  END SUBROUTINE mat_two_pow_mod_m

  ! ***************************************************************************

  SUBROUTINE mat_vec_mod_m(a,s,v,m)

    ! Returns v = MODULO(a*s,m). Assumes that -m < s(i) < m.

    REAL(KIND=dp), DIMENSION(3,3), INTENT(IN) :: a
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)   :: s
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT)  :: v
    REAL(KIND=dp), INTENT(IN)                 :: m

    ! Local variables

    REAL(KIND=dp) :: a1,a2,c
    INTEGER       :: i,j

    ! -------------------------------------------------------------------------

    v = 0.0_dp

    DO i=1,3
      DO j=1,3
        a2 = a(i,j)
        c = v(i)
        v(i) = a2*s(j) + c
        IF ((v(i) >= two53).OR.(v(i) <= -two53)) THEN
          a1 = INT(a2/two17)
          a2 = a2 - a1*two17
          v(i) = a1*s(j)
          a1 = INT(v(i)/m)
          v(i) = v(i) - a1*m
          v(i) = v(i)*two17 + a2*s(j) + c
        END IF
        a1 = INT(v(i)/m)
        v(i) = v(i) - a1*m
        IF (v(i) < 0.0_dp) v(i) = v(i) + m
      END DO
    END DO

  END SUBROUTINE mat_vec_mod_m

  ! ***************************************************************************

  FUNCTION next_integer_random_number(rng_stream,low,high) RESULT(u)

    ! Get the next integer random number between low and high from the stream
    ! rng_stream.

    TYPE(rng_stream_type), POINTER :: rng_stream
    INTEGER, INTENT(IN)            :: low,high

    INTEGER :: u

    ! -------------------------------------------------------------------------

    u = low + INT(next_real_random_number(rng_stream)*REAL(high - low + 1,dp))

  END FUNCTION next_integer_random_number

  ! ***************************************************************************

  FUNCTION next_real_random_number(rng_stream) RESULT(u)

    ! Get the next real random number between 0 and 1 from the stream
    ! rng_stream.

    TYPE(rng_stream_type), POINTER :: rng_stream

    REAL(KIND=dp) :: u

    ! -------------------------------------------------------------------------

    IF (rng_stream%extended_precision) THEN
      u = rn53(rng_stream)
    ELSE
      u = rn32(rng_stream)
    END IF

  END FUNCTION next_real_random_number

  ! ***************************************************************************

  FUNCTION next_rng_seed(seed) RESULT(next_seed)

    ! Get the seed for the next RNG stream w.r.t. a given seed.
    ! If the optional argument seed is missing, then the default seed is
    ! returned.

    REAL(KIND=dp), DIMENSION(3,2), OPTIONAL, INTENT(IN) :: seed

    REAL(KIND=dp), DIMENSION(3,2) :: next_seed

    ! -------------------------------------------------------------------------

    IF (PRESENT(seed)) THEN
      CALL mat_vec_mod_m(a1p127,seed(:,1),next_seed(:,1),m1)
      CALL mat_vec_mod_m(a2p127,seed(:,2),next_seed(:,2),m2)
    ELSE
      next_seed = 12345.0_dp ! default seed
    END IF

  END FUNCTION next_rng_seed

  ! ***************************************************************************

  SUBROUTINE reset_rng_stream(rng_stream)

    ! Reset a random number stream to its initial state.

    TYPE(rng_stream_type), POINTER :: rng_stream

    ! -------------------------------------------------------------------------

    rng_stream%cg = rng_stream%ig
    rng_stream%bg = rng_stream%ig

  END SUBROUTINE reset_rng_stream

  ! ***************************************************************************

  SUBROUTINE reset_rng_substream(rng_stream)

    ! Reset a random number stream to the beginning of its current substream.

    TYPE(rng_stream_type), POINTER :: rng_stream

    ! -------------------------------------------------------------------------

    rng_stream%cg = rng_stream%bg

  END SUBROUTINE reset_rng_substream

  ! ***************************************************************************

  SUBROUTINE reset_to_next_rng_substream(rng_stream)

    ! Reset a random number stream to the beginning of its next substream.

    TYPE(rng_stream_type), POINTER :: rng_stream

    ! Local variables

    REAL(KIND=dp), DIMENSION(3,2) :: u

    ! -------------------------------------------------------------------------

    u = 0.0_dp

    CALL mat_vec_mod_m(a1p76,rng_stream%bg(:,1),u(:,1),m1)
    CALL mat_vec_mod_m(a2p76,rng_stream%bg(:,2),u(:,2),m2)

    rng_stream%bg = u
    rng_stream%cg = u

  END SUBROUTINE reset_to_next_rng_substream

  ! ***************************************************************************

  FUNCTION rn32(rng_stream) RESULT(u)

    ! Generate the next random number with standard precision (32 bits).

    TYPE(rng_stream_type), POINTER :: rng_stream

    REAL(KIND=dp) :: u

    ! Local variables

    REAL(KIND=dp) :: p1,p2
    INTEGER       :: k

    ! -------------------------------------------------------------------------

    ! Component 1

    p1 = a12*rng_stream%cg(2,1) - a13n*rng_stream%cg(1,1)
    k = INT(p1/m1)
    p1 = p1 - k*m1
    IF (p1 < 0.0_dp) p1 = p1 + m1
    rng_stream%cg(1,1) = rng_stream%cg(2,1)
    rng_stream%cg(2,1) = rng_stream%cg(3,1)
    rng_stream%cg(3,1) = p1

    ! Component 2

    p2 = a21*rng_stream%cg(3,2) - a23n*rng_stream%cg(1,2)
    k = INT(p2/m2)
    p2 = p2 - k*m2
    IF (p2 < 0.0_dp) p2 = p2 + m2
    rng_stream%cg(1,2) = rng_stream%cg(2,2)
    rng_stream%cg(2,2) = rng_stream%cg(3,2)
    rng_stream%cg(3,2) = p2

    ! Combination

    IF (p1 > p2) THEN
      u = (p1 - p2)*norm
    ELSE
      u = (p1 - p2 + m1)*norm
    END IF

    IF (rng_stream%antithetic) u = 1.0_dp - u

  END FUNCTION rn32

  ! ***************************************************************************

  FUNCTION rn53(rng_stream) RESULT(u)

    ! Generate the next random number with extended precision (53 bits).

    TYPE(rng_stream_type), POINTER :: rng_stream

    REAL(KIND=dp) :: u

    ! -------------------------------------------------------------------------

    u = rn32(rng_stream)

    ! Note: rn32 returns 1 - u in the antithetic case

    IF (rng_stream%antithetic) THEN
      u = u + (rn32(rng_stream) - 1.0_dp)*fact
      IF (u < 0.0_dp) u = u + 1.0_dp
    ELSE
      u = u + rn32(rng_stream)*fact
      IF (u >= 1.0_dp) u = u - 1.0_dp
    END IF

  END FUNCTION rn53

  ! ***************************************************************************

  SUBROUTINE set_rng_stream(rng_stream,name,bg,cg,ig,seed,antithetic,&
                            extended_precision)

    ! Set the components of a RNG stream.
    ! NOTE: The manipulation of an active RNG stream is discouraged.

    TYPE(rng_stream_type), POINTER                :: rng_stream
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)        :: name
    INTEGER, DIMENSION(3,2), OPTIONAL, INTENT(IN) :: bg,cg,ig,seed
    LOGICAL, OPTIONAL, INTENT(IN)                 :: antithetic,&
                                                     extended_precision

    ! -------------------------------------------------------------------------

    IF (PRESENT(name)) rng_stream%name = name
    IF (PRESENT(bg)) rng_stream%bg = bg
    IF (PRESENT(cg)) rng_stream%cg = cg
    IF (PRESENT(ig)) rng_stream%ig = ig
    IF (PRESENT(seed)) THEN
      ! Sets the initial seed of the stream to seed
      ! NOTE: The use of this method is discouraged
      CALL check_seed(seed)
      rng_stream%ig = seed
      rng_stream%cg = seed
      rng_stream%bg = seed
    END IF
    IF (PRESENT(antithetic)) rng_stream%antithetic = antithetic
    IF (PRESENT(extended_precision)) THEN
      rng_stream%extended_precision = extended_precision
    END IF

  END SUBROUTINE set_rng_stream

  ! ***************************************************************************

  SUBROUTINE write_rng_matrices(output_unit)

    ! Write the transformation matrices of the two MRG components (raised to
    ! the powers -1, 1, 2**76, and 2**127) the a logical output unit.

    INTEGER, INTENT(IN) :: output_unit

    ! Local variables

    CHARACTER(LEN=40) :: fmtstr
    INTEGER           :: i,j

    ! -------------------------------------------------------------------------

    ! Print the transformation matrices for both components

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
      "TRANSFORMATION MATRICES FOR THE PARALLEL (PSEUDO)RANDOM NUMBER "//&
      "GENERATOR"

    fmtstr = "(/,T4,A,/,/,(2X,3F14.1))"

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A1",((a1p0(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A2",((a2p0(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A1**(2**76)",((a1p76(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A2**(2**76)",((a2p76(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A1**(2**127)",((a1p127(i,j),j=1,3),i=1,3)

    WRITE (UNIT=output_unit,FMT=fmtstr)&
      "A2**(2**127)",((a2p127(i,j),j=1,3),i=1,3)

  END SUBROUTINE write_rng_matrices

  ! ***************************************************************************

  SUBROUTINE write_rng_stream(rng_stream,output_unit)

    TYPE(rng_stream_type), POINTER :: rng_stream
    INTEGER, INTENT(IN)            :: output_unit

    ! -------------------------------------------------------------------------

    WRITE (UNIT=output_unit,FMT="(/,T2,A,/)")&
      "Random number stream <"//TRIM(rng_stream%name)//">:"

    IF (rng_stream%antithetic) THEN
      WRITE (UNIT=output_unit,FMT="(T4,A)") "Antithetic:    yes"
    ELSE
      WRITE (UNIT=output_unit,FMT="(T4,A)") "Antithetic:     no"
    END IF

    IF (rng_stream%extended_precision) THEN
      WRITE (UNIT=output_unit,FMT="(T4,A)") "Precision:  53 Bit"
    ELSE
      WRITE (UNIT=output_unit,FMT="(T4,A)") "Precision:  32 Bit"
    END IF

    WRITE (UNIT=output_unit,FMT="(/,T4,A,/,/,(T4,A,3F20.1))")&
      "Initial state of the stream:",&
      "Component 1:",rng_stream%ig(:,1),&
      "Component 2:",rng_stream%ig(:,2)

    WRITE (UNIT=output_unit,FMT="(/,T4,A,/,/,(T4,A,3F20.1))")&
      "Initial state of the current substream:",&
      "Component 1:",rng_stream%bg(:,1),&
      "Component 2:",rng_stream%bg(:,2)

    WRITE (UNIT=output_unit,FMT="(/,T4,A,/,/,(T4,A,3F20.1))")&
      "Current state of the stream:",&
      "Component 1:",rng_stream%cg(:,1),&
      "Component 2:",rng_stream%cg(:,2)

  END SUBROUTINE write_rng_stream

  ! ***************************************************************************

END MODULE parallel_rng_types
