!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Functions which are common to different hash tables
!>         Derived from qs_fb_hash_table_types and qs_fb_hash_table_types (Mark Tucker, Jun 2016)
! **************************************************************************************************
MODULE qs_hash_table_functions

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

! public methods
   PUBLIC :: hash_table_matching_prime

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_hash_table_functions'

CONTAINS

! **************************************************************************************************
!> \brief Find a prime number equal or larger than ii
!> \param ii   : input integer
!> \retval res : the prime number
! **************************************************************************************************
   PURE FUNCTION hash_table_matching_prime(ii) RESULT(res)
      INTEGER, INTENT(IN)                                :: ii
      INTEGER                                            :: res

      ! even numbers are not prime, so no point testing them, so increment by 2 each time starting
      ! from an odd number greater or equal to ii (as noted in \brief)
      res = ii+1-MOD(ii, 2)

      DO WHILE (.NOT. is_positive_prime(res))
         res = res+2
      END DO
   END FUNCTION hash_table_matching_prime

! **************************************************************************************************
!> \brief Check if a number is a positive prime
!> \param num  : number to check
!> \retval res : returns TRUE if num is a positive prime, FALSE otherwise
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! **************************************************************************************************
   PURE FUNCTION is_positive_prime(num) RESULT(res)
      INTEGER, INTENT(IN)                                :: num
      LOGICAL                                            :: res

      INTEGER                                            :: ii

      IF (num .LE. 3) THEN
         res = .FALSE.
         RETURN
      END IF
      IF (MOD(num, 2) == 0 .OR. MOD(num, 3) == 0) THEN
         res = .FALSE.
         RETURN
      END IF

      ! all primes > 3 are of the form 6*kk +/- 1, kk=1,2,3...
      ! (although not all 6*kk +/- 1 is a prime);
      ! and we only have to check factors less than and equal to SQRT(num)
      ii = 5
      DO WHILE (ii*ii .LE. num)
         IF (MOD(num, ii) == 0 .OR. MOD(num, ii+2) == 0) THEN
            res = .FALSE.
            RETURN
         END IF
         ii = ii+6
      END DO
      res = .TRUE.
   END FUNCTION is_positive_prime

END MODULE qs_hash_table_functions

