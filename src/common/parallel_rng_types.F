!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Parallel (pseudo)random number generator (RNG) for multiple streams
!>      and substreams of random numbers.
!>
!>      In detail, this RNG provides 2**64 random number streams each with a
!>      length of 2**127 resulting in a length of 2**191 for the total RNG.
!>      Moreover, each stream is divided in 2**51 substream of length 2**76.
!>      The stream lengths refer to the default precision of 32 bit random
!>      number, but also an extended precision of 53 bit per random number
!>      can be requested. In this case, two 32 bit random numbers are used
!>      to generate a 53 bit random number and therefore the stream length
!>      is halved when extended precision are requested.
!>
!>      Usage hint:
!>
!>      type(rng_stream_type) :: rng_stream
!>      rng_stream = rng_stream_type(name, ..., error=error)
!>
!>      to generate the first stream. Optionally, you may define a different
!>      seed or create a stream of extended precision (53 bits). Then
!>
!>      type(rng_stream_type) :: next_rng_stream
!>      next_rng_stream = rng_stream_type(name, last_rng_stream=rng_stream)
!>
!>      to create all the following RNG streams w.r.t. the previous stream.
!>      The command line
!>
!>      x = rng_stream%next(error=error)
!>
!>      will provide the next real random number x between 0 and 1 and
!>
!>      ix = rng_stream%next(low, high, error=error)
!>
!>      the next integer random number ix between low and high from stream
!>      rng_stream. The default distribution type is a uniform distribution
!>      [0,1], but also other distribution types are available (see below).
!>
!> \par Literature
!>      P. L'Ecuyer, R. Simard, E. J. Chen, and W. D. Kelton,
!>      "An object-oriented random-number package with many long streams
!>       and substreams", Operations Research 50(6), 1073-1075 (2002)
!> \author C++ code converted to Fortran 90/95 (18.05.2005, Matthias Krack)
! **************************************************************************************************
MODULE parallel_rng_types

   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE string_utilities,                ONLY: compress
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters in this module

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: rng_record_format = "(A40,I2,3L2,ES25.16,18F20.1)"
   INTEGER, PARAMETER                   :: rng_record_length = 433
   INTEGER, PARAMETER                   :: rng_name_length = 40

   ! Distribution types:

   ! GAUSSIAN: Gaussian distribution with zero mean and unit variance
   ! UNIFORM:  Uniform distribution [0,1] with 1/2 mean (default)

   INTEGER, PARAMETER       :: GAUSSIAN = 1, &
                               UNIFORM = 2

   REAL(KIND=dp), PARAMETER :: norm = 2.328306549295727688e-10_dp, &
                               m1 = 4294967087.0_dp, &
                               m2 = 4294944443.0_dp, &
                               a12 = 1403580.0_dp, &
                               a13n = 810728.0_dp, &
                               a21 = 527612.0_dp, &
                               a23n = 1370589.0_dp, &
                               two17 = 131072.0_dp, & ! 2**17
                               two53 = 9007199254740992.0_dp, & ! 2**53
                               fact = 5.9604644775390625e-8_dp ! 1/2**24

   !&<
   ! The following are the transition matrices of the two MRG components
   ! (in matrix form), raised to the powers 1, 2**76, 2**127, and -1

   ! Transition matrix for the first component raised to the power 2**0
   REAL(KIND=dp), DIMENSION(3, 3), PARAMETER :: a1p0 = RESHAPE([ &
         0.0_dp, 0.0_dp, -810728.0_dp, &
         1.0_dp, 0.0_dp, 1403580.0_dp, &
         0.0_dp, 1.0_dp,       0.0_dp &
         ], [3,3])

   ! Transition matrix for the second component raised to the power 2**0
   REAL(KIND=dp), DIMENSION(3, 3), PARAMETER :: a2p0 = RESHAPE([ &
         0.0_dp, 0.0_dp, -1370589.0_dp, &
         1.0_dp, 0.0_dp,        0.0_dp, &
         0.0_dp, 1.0_dp,   527612.0_dp &
         ], [3,3])

   ! Transition matrix for the first component raised to the power 2**76
   REAL(KIND=dp), DIMENSION(3, 3), PARAMETER :: a1p76 = RESHAPE([ &
           82758667.0_dp, 3672831523.0_dp, 3672091415.0_dp, &
         1871391091.0_dp,   69195019.0_dp, 3528743235.0_dp, &
         4127413238.0_dp, 1871391091.0_dp,   69195019.0_dp &
         ], [3,3])

   ! Transition matrix for the second component raised to the power 2**76
   REAL(KIND=dp), DIMENSION(3, 3), PARAMETER :: a2p76 = RESHAPE([ &
         1511326704.0_dp, 4292754251.0_dp, 3859662829.0_dp, &
         3759209742.0_dp, 1511326704.0_dp, 4292754251.0_dp, &
         1610795712.0_dp, 3889917532.0_dp, 3708466080.0_dp &
         ], [3,3])

   ! Transition matrix for the first component raised to the power 2**127
   REAL(KIND=dp), DIMENSION(3, 3), PARAMETER :: a1p127 = RESHAPE([ &
         2427906178.0_dp,  226153695.0_dp, 1988835001.0_dp, &
         3580155704.0_dp, 1230515664.0_dp,  986791581.0_dp, &
          949770784.0_dp, 3580155704.0_dp, 1230515664.0_dp &
         ], [3,3])

   ! Transition matrix for the second component raised to the power 2**127
   REAL(KIND=dp), DIMENSION(3, 3), PARAMETER :: a2p127 = RESHAPE([ &
         1464411153.0_dp,   32183930.0_dp, 2824425944.0_dp, &
          277697599.0_dp, 1464411153.0_dp,   32183930.0_dp, &
         1610723613.0_dp, 1022607788.0_dp, 2093834863.0_dp &
         ], [3,3])

   ! Inverse of a1p0
   REAL(KIND=dp), DIMENSION(3, 3), PARAMETER :: inv_a1 = RESHAPE([ &
          184888585.0_dp, 1.0_dp, 0.0_dp, &
                  0.0_dp, 0.0_dp, 1.0_dp, &
         1945170933.0_dp, 0.0_dp, 0.0_dp &
         ], [3,3])

   ! Inverse of a2p0
   REAL(KIND=dp), DIMENSION(3, 3), PARAMETER :: inv_a2 = RESHAPE([ &
                  0.0_dp, 1.0_dp, 0.0_dp, &
          360363334.0_dp, 0.0_dp, 1.0_dp, &
         4225571728.0_dp, 0.0_dp, 0.0_dp &
         ], [3,3])
   !&>

   ! Data type definitions

   ! Information on a stream. The arrays bg, cg, and ig contain the current
   ! state of the stream, the starting state of the current substream, and the
   ! starting state of the stream. This stream generates antithetic variates
   ! if antithetic = .TRUE.. It also generates numbers with extended precision
   ! (53 bits, if machine follows IEEE 754 standard), if extended_precision =
   ! .TRUE., otherwise, numbers with 32 bits precision are generated.

   TYPE rng_stream_type
      PRIVATE
      ! the name could be an allocatable, but gfortran (even with 9.1) does not properly implement
      ! automatic deallocation of it and a `final`routine which would do it triggers an ICE in 7.4.1
      CHARACTER(LEN=rng_name_length) :: name
      INTEGER                        :: distribution_type = UNIFORM
      ! ig: initial state, cg: current state, bg: initial state of the substream
      REAL(KIND=dp), DIMENSION(3, 2) :: bg, cg, ig
      LOGICAL                        :: antithetic = .FALSE., extended_precision = .FALSE.
      ! only used for distribution type GAUSSIAN
      REAL(KIND=dp)                  :: buffer = 0.0_dp
      LOGICAL                        :: buffer_filled = .FALSE.

   CONTAINS
      PROCEDURE, PASS(self) :: fill_1
      PROCEDURE, PASS(self) :: fill_2
      PROCEDURE, PASS(self) :: fill_3
      GENERIC, PUBLIC :: fill => fill_1, fill_2, fill_3

      PROCEDURE, PASS(self) :: next_int
      PROCEDURE, PASS(self) :: next_real
      GENERIC, PUBLIC :: next => next_int, next_real

      PROCEDURE, PASS(self), PUBLIC :: dump
      PROCEDURE, PASS(self), PUBLIC :: write
      PROCEDURE, PASS(self), PUBLIC :: advance
      PROCEDURE, PASS(self), PUBLIC :: set
      PROCEDURE, PASS(self), PUBLIC :: get
      PROCEDURE, PASS(self), PUBLIC :: reset
      PROCEDURE, PASS(self), PUBLIC :: reset_to_substream
      PROCEDURE, PASS(self), PUBLIC :: reset_to_next_substream
   END TYPE rng_stream_type

   INTERFACE rng_stream_type
      MODULE PROCEDURE :: rng_stream_constructor
   END INTERFACE

   TYPE rng_stream_p_type
      TYPE(rng_stream_type), POINTER :: stream
   END TYPE rng_stream_p_type

   ! Public parameters

   PUBLIC :: rng_record_length, &
             rng_name_length, &
             GAUSSIAN, &
             UNIFORM

   ! Public data types

   PUBLIC :: rng_stream_p_type, &
             rng_stream_type

   ! Public subroutines

   PUBLIC :: check_rng, &
             next_rng_seed, &
             write_rng_matrices, &
             rng_stream_type_from_record

CONTAINS

! **************************************************************************************************
!> \brief Advance the state by n steps, i.e. jump n steps forward, if n > 0, or backward if n < 0.
!> \param self ...
!> \param e IF e > 0, let n = 2**e + c, IF e < 0, let n = -2**(-e) + c, IF e = 0, let n = c
!> \param c ...
!> \note The use of this method is discouraged
! **************************************************************************************************
   SUBROUTINE advance(self, e, c)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self
      INTEGER, INTENT(IN)                                :: e, c

      REAL(KIND=dp), DIMENSION(3, 2)                     :: x
      REAL(KIND=dp), DIMENSION(3, 3)                     :: u1, u2, v1, v2, w1, w2

      u1 = 0.0_dp
      u2 = 0.0_dp
      v1 = 0.0_dp
      v2 = 0.0_dp
      w1 = 0.0_dp
      w2 = 0.0_dp

      IF (e > 0) THEN
         CALL mat_two_pow_mod_m(a1p0, u1, m1, e)
         CALL mat_two_pow_mod_m(a2p0, u2, m2, e)
      ELSE IF (e < 0) THEN
         CALL mat_two_pow_mod_m(inv_a1, u1, m1, -e)
         CALL mat_two_pow_mod_m(inv_a2, u2, m2, -e)
      END IF

      IF (c >= 0) THEN
         CALL mat_pow_mod_m(a1p0, v1, m1, c)
         CALL mat_pow_mod_m(a2p0, v2, m2, c)
      ELSE
         CALL mat_pow_mod_m(inv_a1, v1, m1, -c)
         CALL mat_pow_mod_m(inv_a2, v2, m2, -c)
      END IF

      IF (e == 0) THEN
         w1 = v1
         w2 = v2
      ELSE
         CALL mat_mat_mod_m(u1, v1, w1, m1)
         CALL mat_mat_mod_m(u2, v2, w2, m2)
      END IF

      x = 0.0_dp

      CALL mat_vec_mod_m(w1, self%cg(:, 1), x(:, 1), m1)
      CALL mat_vec_mod_m(w2, self%cg(:, 2), x(:, 2), m2)

      self%cg = x
   END SUBROUTINE advance

! **************************************************************************************************
!> \brief ...
!> \param output_unit ...
!> \param ionode ...
! **************************************************************************************************
   SUBROUTINE check_rng(output_unit, ionode)

      ! Check the parallel (pseudo)random number generator (RNG).

      INTEGER, INTENT(IN)                                :: output_unit
      LOGICAL, INTENT(IN)                                :: ionode

      INTEGER                                            :: i, sumi
      REAL(KIND=dp)                                      :: sum, sum3
      REAL(KIND=dp), DIMENSION(3, 2)                     :: germe
      TYPE(rng_stream_type)                              :: cantor, g1, g2, g3, galois, laplace, &
                                                            poisson

      ! -------------------------------------------------------------------------
      ! Test 1

      ! Create RNG test streams

      g1 = rng_stream_type("g1")
      g2 = rng_stream_type("g2", g1)
      g3 = rng_stream_type("g3", g2)

      IF (ionode) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            "RESULTS OF THE PSEUDO(RANDOM) NUMBER GENERATOR TEST RUNS", &
            "Initial states of the (pseudo)random number streams (test 1):"
         CALL g1%write(output_unit)
         CALL g2%write(output_unit)
         CALL g3%write(output_unit)
      END IF

      sum = g2%next() + g3%next()

      CALL g1%advance(5, 3)
      sum = sum + g1%next()

      CALL g1%reset()
      DO i = 1, 35
         CALL g1%advance(0, 1)
      END DO
      sum = sum + g1%next()

      CALL g1%reset()

      sumi = 0
      DO i = 1, 35
         sumi = sumi + g1%next(1, 10)
      END DO
      sum = sum + sumi/100.0_dp

      sum3 = 0.0_dp
      DO i = 1, 100
         sum3 = sum3 + g3%next()
      END DO
      sum = sum + sum3/10.0_dp

      CALL g3%reset()
      DO i = 1, 5
         sum = sum + g3%next()
      END DO

      CALL g3%reset()
      DO i = 1, 4
         CALL g3%reset_to_next_substream()
      END DO
      DO i = 1, 5
         sum = sum + g3%next()
      END DO

      CALL g3%reset_to_substream()
      DO i = 1, 5
         sum = sum + g3%next()
      END DO

      CALL g2%reset_to_next_substream()
      sum3 = 0.0_dp
      DO i = 1, 100000
         sum3 = sum3 + g2%next()
      END DO
      sum = sum + sum3/10000.0_dp

      CALL g3%set(antithetic=.TRUE.)
      sum3 = 0.0_dp
      DO i = 1, 100000
         sum3 = sum3 + g3%next()
      END DO
      sum = sum + sum3/10000.0_dp

      IF (ionode) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            "Final states of the (pseudo)random number streams (test 1):"
         CALL g1%write(output_unit)
         CALL g2%write(output_unit)
         CALL g3%write(output_unit)
         WRITE (UNIT=output_unit, FMT="(/,(T2,A))") &
            "This test routine should print for test 1 the number 25.342059"
         WRITE (UNIT=output_unit, FMT="(T2,A,F10.6)") &
            "The actual result of test 1 is                      ", sum
      END IF

      ! -------------------------------------------------------------------------
      ! Test 2

      germe(:, :) = 1

      poisson = rng_stream_type("Poisson", seed=germe)
      laplace = rng_stream_type("Laplace", poisson)
      galois = rng_stream_type("Galois", laplace)
      cantor = rng_stream_type("Cantor", galois)

      IF (ionode) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            "Initial states of the (pseudo)random number streams (test 2):"
         CALL poisson%write(output_unit)
         CALL laplace%write(output_unit)
         CALL galois%write(output_unit)
         CALL cantor%write(output_unit)
      END IF

      sum = sum + poisson%next() + laplace%next() + galois%next() + cantor%next()

      CALL galois%advance(-127, 0)
      sum = sum + galois%next()

      CALL galois%reset_to_next_substream()
      CALL galois%set(extended_precision=.TRUE.)
      sum3 = 0.0_dp
      DO i = 1, 100000
         sum3 = sum3 + galois%next()
      END DO
      sum = sum + sum3/10000.0_dp

      CALL galois%set(antithetic=.TRUE.)
      sum3 = 0.0_dp
      DO i = 1, 100000
         sum3 = sum3 + galois%next()
      END DO
      sum = sum + sum3/10000.0_dp
      CALL galois%set(antithetic=.FALSE.)

      CALL galois%set(extended_precision=.FALSE.)
      sum = sum + poisson%next() + laplace%next() + galois%next() + cantor%next()

      IF (ionode) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            "Final states of the (pseudo)random number streams (test 2):"
         CALL poisson%write(output_unit)
         CALL laplace%write(output_unit)
         CALL galois%write(output_unit)
         CALL cantor%write(output_unit)
         WRITE (UNIT=output_unit, FMT="(/,(T2,A))") &
            "This test routine should print for test 2 the number 39.697547"
         WRITE (UNIT=output_unit, FMT="(T2,A,F10.6)") &
            "The actual result of test 2 is                      ", sum
      END IF

   END SUBROUTINE check_rng

! **************************************************************************************************
!> \brief Check that the seeds are legitimate values.
!> \param seed ...
! **************************************************************************************************
   SUBROUTINE check_seed(seed)
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(IN)         :: seed

      CHARACTER(LEN=*), PARAMETER :: fmtstr = "(A,I1,A,ES23.14,A,ES23.14)"

      CHARACTER(LEN=default_string_length)               :: message
      INTEGER                                            :: i

      DO i = 1, 3

         ! Check condition: 0 <= seed(:,1) < m1

         IF (seed(i, 1) < 0.0_dp) THEN
            WRITE (UNIT=message, FMT=fmtstr) &
               "seed(", i, ",1) = ", seed(i, 1), " < ", 0.0_dp
            CALL compress(message)
            CPABORT(message)
         END IF
         IF (seed(i, 1) >= m1) THEN
            WRITE (UNIT=message, FMT=fmtstr) &
               "seed(", i, ",1) = ", seed(i, 1), " >= ", m1
            CALL compress(message)
            CPABORT(message)
         END IF

         ! Check condition: 0 <= seed(:,2) < m2

         IF (seed(i, 2) < 0.0_dp) THEN
            WRITE (UNIT=message, FMT=fmtstr) &
               "seed(", i, ",2) = ", seed(i, 2), " < ", 0.0_dp
            CALL compress(message)
            CPABORT(message)
         END IF
         IF (seed(i, 2) >= m2) THEN
            WRITE (UNIT=message, FMT=fmtstr) &
               "seed(", i, ",2) = ", seed(i, 2), " >= ", m2
            CALL compress(message)
            CPABORT(message)
         END IF

      END DO

      ! Check condition: first or second seed is 0

      IF (ALL(seed(:, 1) < 1.0_dp)) THEN
         CPABORT("First seed = 0")
      END IF

      IF (ALL(seed(:, 2) < 1.0_dp)) THEN
         CPABORT("Second seed = 0")
      END IF

   END SUBROUTINE check_seed

! **************************************************************************************************
!> \brief Create a new RNG stream.
!> \param name ...
!> \param last_rng_stream ...
!> \param distribution_type ...
!> \param seed ...
!> \param antithetic ...
!> \param extended_precision ...
!> \return ...
! **************************************************************************************************
   FUNCTION rng_stream_constructor(name, last_rng_stream, distribution_type, seed, antithetic, extended_precision) &
      RESULT(rng_stream)

      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(rng_stream_type), INTENT(IN), OPTIONAL        :: last_rng_stream
      INTEGER, INTENT(IN), OPTIONAL                      :: distribution_type
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(IN), &
         OPTIONAL                                        :: seed
      LOGICAL, INTENT(IN), OPTIONAL                      :: antithetic, extended_precision
      TYPE(rng_stream_type)                              :: rng_stream

      IF (LEN_TRIM(name) .GT. rng_name_length) &
         CPABORT("given random number generator name is too long")

      rng_stream%name = TRIM(name)

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

      IF (PRESENT(distribution_type)) THEN
         SELECT CASE (distribution_type)
         CASE (GAUSSIAN)
            rng_stream%distribution_type = GAUSSIAN
         CASE (UNIFORM)
            rng_stream%distribution_type = UNIFORM
         CASE DEFAULT
            CPABORT("Invalid distribution type specified")
         END SELECT
      ELSE IF (PRESENT(last_rng_stream)) THEN
         rng_stream%distribution_type = last_rng_stream%distribution_type
      END IF

      IF (PRESENT(antithetic)) THEN
         rng_stream%antithetic = antithetic
      ELSE IF (PRESENT(last_rng_stream)) THEN
         rng_stream%antithetic = last_rng_stream%antithetic
      END IF

      IF (PRESENT(extended_precision)) THEN
         rng_stream%extended_precision = extended_precision
      ELSE IF (PRESENT(last_rng_stream)) THEN
         rng_stream%extended_precision = last_rng_stream%extended_precision
      END IF
   END FUNCTION rng_stream_constructor

! **************************************************************************************************
!> \brief Create a RNG stream from a record given as an internal file (string).
!> \param rng_record ...
!> \return ...
! **************************************************************************************************
   FUNCTION rng_stream_type_from_record(rng_record) RESULT(rng_stream)
      CHARACTER(LEN=rng_record_length), INTENT(IN)       :: rng_record
      TYPE(rng_stream_type)                              :: rng_stream

      READ (UNIT=rng_record, FMT=rng_record_format) &
         rng_stream%name, &
         rng_stream%distribution_type, &
         rng_stream%antithetic, &
         rng_stream%extended_precision, &
         rng_stream%buffer_filled, &
         rng_stream%buffer, &
         rng_stream%cg, &
         rng_stream%bg, &
         rng_stream%ig
   END FUNCTION rng_stream_type_from_record

! **************************************************************************************************
!> \brief Dump a RNG stream as a record given as an internal file (string).
!> \param self ...
!> \param rng_record ...
! **************************************************************************************************
   SUBROUTINE dump(self, rng_record)
      CLASS(rng_stream_type), INTENT(IN)                 :: self
      CHARACTER(LEN=rng_record_length), INTENT(OUT)      :: rng_record

      rng_record = " "
      WRITE (UNIT=rng_record, FMT=rng_record_format) &
         self%name, &
         self%distribution_type, &
         self%antithetic, &
         self%extended_precision, &
         self%buffer_filled, &
         self%buffer, &
         self%cg, &
         self%bg, &
         self%ig
   END SUBROUTINE dump

! **************************************************************************************************
!> \brief Get the components of a RNG stream.
!> \param self ...
!> \param name ...
!> \param distribution_type ...
!> \param bg ...
!> \param cg ...
!> \param ig ...
!> \param antithetic ...
!> \param extended_precision ...
!> \param buffer ...
!> \param buffer_filled ...
!> \par History
!>      2009-11-04 changed bg, cg and ig type from INTEGER, DIMENSION(3, 2)
!>      to REAL(KIND=dp), DIMENSION(3, 2) [lwalewski]
!>      2009-11-09 getting the buffer and buffer_filled components
!>      added [lwalewski]
! **************************************************************************************************
   SUBROUTINE get(self, name, distribution_type, bg, cg, ig, &
                  antithetic, extended_precision, &
                  buffer, buffer_filled)

      CLASS(rng_stream_type), INTENT(IN)                 :: self
      CHARACTER(LEN=rng_name_length), INTENT(OUT), OPTIONAL :: name
      INTEGER, INTENT(OUT), OPTIONAL                     :: distribution_type
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(OUT), &
         OPTIONAL                                        :: bg, cg, ig
      LOGICAL, INTENT(OUT), OPTIONAL                     :: antithetic, extended_precision
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: buffer
      LOGICAL, INTENT(OUT), OPTIONAL                     :: buffer_filled

      IF (PRESENT(name)) name = self%name
      IF (PRESENT(distribution_type)) &
         distribution_type = self%distribution_type
      IF (PRESENT(bg)) bg = self%bg
      IF (PRESENT(cg)) cg = self%cg
      IF (PRESENT(ig)) ig = self%ig
      IF (PRESENT(antithetic)) antithetic = self%antithetic
      IF (PRESENT(extended_precision)) &
         extended_precision = self%extended_precision
      IF (PRESENT(buffer)) buffer = self%buffer
      IF (PRESENT(buffer_filled)) buffer_filled = self%buffer_filled
   END SUBROUTINE get

! **************************************************************************************************
!> \brief Returns c = MODULO(a*b,m)
!> \param a ...
!> \param b ...
!> \param c ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE mat_mat_mod_m(a, b, c, m)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: a, b
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: c
      REAL(KIND=dp), INTENT(IN)                          :: m

      INTEGER                                            :: i

      DO i = 1, 3
         CALL mat_vec_mod_m(a, b(:, i), c(:, i), m)
      END DO

   END SUBROUTINE mat_mat_mod_m

! **************************************************************************************************
!> \brief Compute matrix b = MODULO(a**n,m)
!> \param a ...
!> \param b ...
!> \param m ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE mat_pow_mod_m(a, b, m, n)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: a
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: b
      REAL(KIND=dp), INTENT(IN)                          :: m
      INTEGER, INTENT(IN)                                :: n

      INTEGER                                            :: i
      REAL(KIND=dp), DIMENSION(3, 3)                     :: u, v, w

      ! Initialize: u = v = a; b = I

      w = a

      b(1, 1) = 1.0_dp
      b(2, 1) = 0.0_dp
      b(3, 1) = 0.0_dp
      b(1, 2) = 0.0_dp
      b(2, 2) = 1.0_dp
      b(3, 2) = 0.0_dp
      b(1, 3) = 0.0_dp
      b(2, 3) = 0.0_dp
      b(3, 3) = 1.0_dp

      ! Compute b = MODULO(a**n,m) using the binary decomposition of n

      i = n

      DO
         IF (MODULO(i, 2) /= 0) THEN
            u = w
            v = b
            CALL mat_mat_mod_m(u, v, b, m)
         END IF
         i = i/2
         IF (i == 0) EXIT
         u = w
         v = w
         CALL mat_mat_mod_m(u, v, w, m)
      END DO
   END SUBROUTINE mat_pow_mod_m

! **************************************************************************************************
!> \brief Compute matrix b = MODULO(a**(2**e),m)
!> \param a ...
!> \param b ...
!> \param m ...
!> \param e ...
! **************************************************************************************************
   SUBROUTINE mat_two_pow_mod_m(a, b, m, e)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: a
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: b
      REAL(KIND=dp), INTENT(IN)                          :: m
      INTEGER, INTENT(IN)                                :: e

      INTEGER                                            :: i
      REAL(KIND=dp), DIMENSION(3, 3)                     :: u, v

      b = a

      DO i = 1, e
         u = b
         v = b
         CALL mat_mat_mod_m(u, v, b, m)
      END DO

   END SUBROUTINE mat_two_pow_mod_m

! **************************************************************************************************
!> \brief Returns v = MODULO(a*s,m). Assumes that -m < s(i) < m.
!> \param a ...
!> \param s ...
!> \param v ...
!> \param m ...
! **************************************************************************************************
   SUBROUTINE mat_vec_mod_m(a, s, v, m)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: a
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: s
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: v
      REAL(KIND=dp), INTENT(IN)                          :: m

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: a1, a2, c

      v = 0.0_dp

      DO i = 1, 3
         DO j = 1, 3
            a2 = a(i, j)
            c = v(i)
            v(i) = a2*s(j) + c
            IF ((v(i) >= two53) .OR. (v(i) <= -two53)) THEN
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

! **************************************************************************************************
!> \brief Get the next integer random number between low and high from the stream
!> \param self ...
!> \param low ...
!> \param high ...
!> \return ...
! **************************************************************************************************
   FUNCTION next_int(self, low, high) RESULT(u)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self
      INTEGER, INTENT(IN)                                :: low, high
      INTEGER                                            :: u

      REAL(KIND=dp)                                      :: r

      CPASSERT(self%distribution_type == UNIFORM)

      r = self%next_real()
      u = low + INT(r*REAL(high - low + 1, dp))
   END FUNCTION next_int

! **************************************************************************************************
!> \brief Get the next real random number from the stream rng_stream.
!> \param self ...
!> \param variance variance of the Gaussian distribution (defaults to 1)
!> \return ...
! **************************************************************************************************
   FUNCTION next_real(self, variance) RESULT(u)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: variance
      REAL(KIND=dp)                                      :: u

      REAL(KIND=dp)                                      :: f, r, u1, u2, var

      SELECT CASE (self%distribution_type)
      CASE (GAUSSIAN)
         var = 1.0_dp
         IF (PRESENT(variance)) var = variance
         ! take the random number from the buffer, if the buffer is filled
         IF (self%buffer_filled) THEN
            u = SQRT(var)*self%buffer
            self%buffer_filled = .FALSE.
         ELSE
            DO
               IF (self%extended_precision) THEN
                  u1 = 2.0_dp*rn53(self) - 1.0_dp
                  u2 = 2.0_dp*rn53(self) - 1.0_dp
               ELSE
                  u1 = 2.0_dp*rn32(self) - 1.0_dp
                  u2 = 2.0_dp*rn32(self) - 1.0_dp
               END IF
               r = u1*u1 + u2*u2
               IF ((r > 0.0_dp) .AND. (r < 1.0_dp)) EXIT
            END DO
            ! Box-Muller transformation
            f = SQRT(-2.0_dp*LOG(r)/r)
            u = SQRT(var)*f*u1
            ! save the second random number for the next call
            self%buffer = f*u2
            self%buffer_filled = .TRUE.
         END IF
      CASE (UNIFORM)
         IF (self%extended_precision) THEN
            u = rn53(self)
         ELSE
            u = rn32(self)
         END IF
      END SELECT
   END FUNCTION next_real

! **************************************************************************************************
!> \brief Get the seed for the next RNG stream w.r.t. a given seed.
!> \param seed If the optional argument seed is missing, then the default seed is returned.
!> \return ...
! **************************************************************************************************
   FUNCTION next_rng_seed(seed) RESULT(next_seed)
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(IN), &
         OPTIONAL                                        :: seed
      REAL(KIND=dp), DIMENSION(3, 2)                     :: next_seed

      IF (PRESENT(seed)) THEN
         CALL check_seed(seed)
         CALL mat_vec_mod_m(a1p127, seed(:, 1), next_seed(:, 1), m1)
         CALL mat_vec_mod_m(a2p127, seed(:, 2), next_seed(:, 2), m2)
      ELSE
         next_seed = 12345.0_dp ! default seed
      END IF

   END FUNCTION next_rng_seed

! **************************************************************************************************
!> \brief Fill entity array with random numbers from the RNG stream rng_stream
!> \param self ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE fill_1(self, array)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: array

      INTEGER                                            :: i

      DO i = 1, SIZE(array)
         array(i) = self%next()
      END DO
   END SUBROUTINE fill_1

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE fill_2(self, array)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)         :: array

      INTEGER                                            :: i, j

      DO j = 1, SIZE(array, 2)
         DO i = 1, SIZE(array, 1)
            array(i, j) = self%next()
         END DO
      END DO
   END SUBROUTINE fill_2

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE fill_3(self, array)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)       :: array

      INTEGER                                            :: i, j, k

      DO k = 1, SIZE(array, 3)
         DO j = 1, SIZE(array, 2)
            DO i = 1, SIZE(array, 1)
               array(i, j, k) = self%next()
            END DO
         END DO
      END DO
   END SUBROUTINE fill_3

! **************************************************************************************************
!> \brief Reset a random number stream to its initial state.
!> \param self ...
! **************************************************************************************************
   SUBROUTINE reset(self)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self

      self%cg = self%ig
      self%bg = self%ig
   END SUBROUTINE reset

! **************************************************************************************************
!> \brief Reset a random number stream to the beginning of its current substream.
!> \param self ...
! **************************************************************************************************
   SUBROUTINE reset_to_substream(self)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self

      self%cg = self%bg
   END SUBROUTINE reset_to_substream

! **************************************************************************************************
!> \brief Reset a random number stream to the beginning of its next substream.
!> \param self ...
! **************************************************************************************************
   SUBROUTINE reset_to_next_substream(self)
      CLASS(rng_stream_type), INTENT(INOUT)              :: self

      REAL(KIND=dp), DIMENSION(3, 2)                     :: u

      u = 0.0_dp

      CALL mat_vec_mod_m(a1p76, self%bg(:, 1), u(:, 1), m1)
      CALL mat_vec_mod_m(a2p76, self%bg(:, 2), u(:, 2), m2)

      self%bg = u
      self%cg = u
   END SUBROUTINE reset_to_next_substream

! **************************************************************************************************
!> \brief Generate the next random number with standard precision (32 bits)
!> \param rng_stream ...
!> \return ...
! **************************************************************************************************
   FUNCTION rn32(rng_stream) RESULT(u)
      TYPE(rng_stream_type)                              :: rng_stream
      REAL(KIND=dp)                                      :: u

      INTEGER                                            :: k
      REAL(KIND=dp)                                      :: p1, p2

      ! Component 1

      p1 = a12*rng_stream%cg(2, 1) - a13n*rng_stream%cg(1, 1)
      k = INT(p1/m1)
      p1 = p1 - k*m1
      IF (p1 < 0.0_dp) p1 = p1 + m1
      rng_stream%cg(1, 1) = rng_stream%cg(2, 1)
      rng_stream%cg(2, 1) = rng_stream%cg(3, 1)
      rng_stream%cg(3, 1) = p1

      ! Component 2

      p2 = a21*rng_stream%cg(3, 2) - a23n*rng_stream%cg(1, 2)
      k = INT(p2/m2)
      p2 = p2 - k*m2
      IF (p2 < 0.0_dp) p2 = p2 + m2
      rng_stream%cg(1, 2) = rng_stream%cg(2, 2)
      rng_stream%cg(2, 2) = rng_stream%cg(3, 2)
      rng_stream%cg(3, 2) = p2

      ! Combination

      IF (p1 > p2) THEN
         u = (p1 - p2)*norm
      ELSE
         u = (p1 - p2 + m1)*norm
      END IF

      IF (rng_stream%antithetic) u = 1.0_dp - u

   END FUNCTION rn32

! **************************************************************************************************
!> \brief Generate the next random number with extended precision (53 bits)
!> \param rng_stream ...
!> \return ...
! **************************************************************************************************
   FUNCTION rn53(rng_stream) RESULT(u)
      TYPE(rng_stream_type)                              :: rng_stream
      REAL(KIND=dp)                                      :: u

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

! **************************************************************************************************
!> \brief Set the components of a RNG stream.
!> \param self ...
!> \param name ...
!> \param distribution_type ...
!> \param bg ...
!> \param cg ...
!> \param ig ...
!> \param seed ...
!> \param antithetic ...
!> \param extended_precision ...
!> \param buffer ...
!> \param buffer_filled ...
!> \par History
!>      2009-11-09 setting the buffer and buffer_filled components
!>      added [lwalewski]
! **************************************************************************************************
   SUBROUTINE set(self, name, distribution_type, bg, cg, ig, &
                  seed, antithetic, extended_precision, &
                  buffer, buffer_filled)

      ! NOTE: The manipulation of an active RNG stream is discouraged.

      CLASS(rng_stream_type), INTENT(INOUT)              :: self
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: name
      INTEGER, INTENT(IN), OPTIONAL                      :: distribution_type
      REAL(KIND=dp), DIMENSION(3, 2), INTENT(IN), &
         OPTIONAL                                        :: bg, cg, ig, seed
      LOGICAL, INTENT(IN), OPTIONAL                      :: antithetic, extended_precision
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: buffer
      LOGICAL, INTENT(IN), OPTIONAL                      :: buffer_filled

      IF (PRESENT(name)) self%name = name
      IF (PRESENT(distribution_type)) THEN
         self%distribution_type = distribution_type
      END IF
      IF (PRESENT(bg)) self%bg = bg
      IF (PRESENT(cg)) self%cg = cg
      IF (PRESENT(ig)) self%ig = ig
      IF (PRESENT(seed)) THEN
         ! Sets the initial seed of the stream to seed
         ! NOTE: The use of this method is discouraged
         CALL check_seed(seed)
         self%ig = seed
         self%cg = seed
         self%bg = seed
      END IF
      IF (PRESENT(antithetic)) self%antithetic = antithetic
      IF (PRESENT(extended_precision)) THEN
         self%extended_precision = extended_precision
      END IF
      IF (PRESENT(buffer)) self%buffer = buffer
      IF (PRESENT(buffer_filled)) self%buffer_filled = buffer_filled
   END SUBROUTINE set

! **************************************************************************************************
!> \brief Write the transformation matrices of the two MRG components (raised to the specified output)
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE write_rng_matrices(output_unit)
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=40)                                  :: fmtstr
      INTEGER                                            :: i, j

      ! Print the transformation matrices for both components

      WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
         "TRANSFORMATION MATRICES FOR THE PARALLEL (PSEUDO)RANDOM NUMBER "// &
         "GENERATOR"

      fmtstr = "(/,T4,A,/,/,(2X,3F14.1))"

      WRITE (UNIT=output_unit, FMT=fmtstr) &
         "A1", ((a1p0(i, j), j=1, 3), i=1, 3)

      WRITE (UNIT=output_unit, FMT=fmtstr) &
         "A2", ((a2p0(i, j), j=1, 3), i=1, 3)

      WRITE (UNIT=output_unit, FMT=fmtstr) &
         "A1**(2**76)", ((a1p76(i, j), j=1, 3), i=1, 3)

      WRITE (UNIT=output_unit, FMT=fmtstr) &
         "A2**(2**76)", ((a2p76(i, j), j=1, 3), i=1, 3)

      WRITE (UNIT=output_unit, FMT=fmtstr) &
         "A1**(2**127)", ((a1p127(i, j), j=1, 3), i=1, 3)

      WRITE (UNIT=output_unit, FMT=fmtstr) &
         "A2**(2**127)", ((a2p127(i, j), j=1, 3), i=1, 3)

   END SUBROUTINE write_rng_matrices

! **************************************************************************************************
!> \brief ...
!> \param self ...
!> \param output_unit ...
!> \param write_all if .TRUE., then print all stream informations (the default is .FALSE.).
! **************************************************************************************************
   SUBROUTINE write (self, output_unit, write_all)
      CLASS(rng_stream_type), INTENT(IN)                 :: self
      INTEGER, INTENT(IN)                                :: output_unit
      LOGICAL, INTENT(IN), OPTIONAL                      :: write_all

      LOGICAL                                            :: my_write_all

      my_write_all = .FALSE.

      IF (PRESENT(write_all)) &
         my_write_all = write_all

      WRITE (UNIT=output_unit, FMT="(/,T2,A,/)") &
         "Random number stream <"//TRIM(self%name)//">:"

      SELECT CASE (self%distribution_type)
      CASE (GAUSSIAN)
         WRITE (UNIT=output_unit, FMT="(T4,A)") &
            "Distribution type: "// &
            "Normal Gaussian distribution with zero mean"
      CASE (UNIFORM)
         WRITE (UNIT=output_unit, FMT="(T4,A)") &
            "Distribution type: "// &
            "Uniform distribution [0,1] with 1/2 mean"
      END SELECT

      IF (self%antithetic) THEN
         WRITE (UNIT=output_unit, FMT="(T4,A)") "Antithetic:        yes"
      ELSE
         WRITE (UNIT=output_unit, FMT="(T4,A)") "Antithetic:        no"
      END IF

      IF (self%extended_precision) THEN
         WRITE (UNIT=output_unit, FMT="(T4,A)") "Precision:         53 Bit"
      ELSE
         WRITE (UNIT=output_unit, FMT="(T4,A)") "Precision:         32 Bit"
      END IF

      IF (my_write_all) THEN

         WRITE (UNIT=output_unit, FMT="(/,T4,A,/,/,(T4,A,3F20.1))") &
            "Initial state of the stream:", &
            "Component 1:", self%ig(:, 1), &
            "Component 2:", self%ig(:, 2)

         WRITE (UNIT=output_unit, FMT="(/,T4,A,/,/,(T4,A,3F20.1))") &
            "Initial state of the current substream:", &
            "Component 1:", self%bg(:, 1), &
            "Component 2:", self%bg(:, 2)

      END IF

      WRITE (UNIT=output_unit, FMT="(/,T4,A,/,/,(T4,A,3F20.1))") &
         "Current state of the stream:", &
         "Component 1:", self%cg(:, 1), &
         "Component 2:", self%cg(:, 2)
   END SUBROUTINE write

END MODULE parallel_rng_types
