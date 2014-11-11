!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Tools usually found in a standard library.
!> \author  Urban Borstnik
!> \date    2011-11-02
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-11-02
! *****************************************************************************
MODULE dbcsr_toollib

  USE dbcsr_array_d_sort,              ONLY: dbcsr_1d_d_sort
  USE dbcsr_array_i4_sort,             ONLY: dbcsr_1d_i4_sort
  USE dbcsr_array_i8_sort,             ONLY: dbcsr_1d_i8_sort
  USE dbcsr_array_s_sort,              ONLY: dbcsr_1d_s_sort
  USE kinds,                           ONLY: int_4,&
                                             int_8

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_toollib'

  ! Block sizes and arrays
  PUBLIC :: dbcsr_unpack_i8_2i4, make_coordinate_tuple
  PUBLIC :: swap
  ! string routines
  PUBLIC :: uppercase
  ! math routines
  PUBLIC :: gcd, lcm
  ! utility routines
  PUBLIC :: sort, joaat_hash

  INTERFACE swap
     MODULE PROCEDURE iswap, bswap
  END INTERFACE

  INTERFACE sort
     MODULE PROCEDURE dbcsr_1d_i4_sort, dbcsr_1d_i8_sort
     MODULE PROCEDURE dbcsr_1d_s_sort, dbcsr_1d_d_sort
  END INTERFACE


CONTAINS

! *****************************************************************************
!> \brief ...
!> \param most ...
!> \param least ...
!> \retval tuple ...
! *****************************************************************************
  ELEMENTAL FUNCTION make_coordinate_tuple (most, least) RESULT (tuple)
    INTEGER, INTENT(IN)                      :: most, least
    INTEGER(KIND=int_8)                      :: tuple

!tuple = IOR (ISHFT (most, 32), least)

    tuple = most
    tuple = IOR (ISHFT(tuple, 32), INT(least, int_8))
  END FUNCTION make_coordinate_tuple


! *****************************************************************************
!> \brief Swaps two integers
!> \param a Integers to swap
!> \param b Integers to swap
! *****************************************************************************
  ELEMENTAL SUBROUTINE iswap(a, b)
    INTEGER, INTENT(INOUT)                   :: a, b

    INTEGER                                  :: tmp

    tmp = a
    a = b
    b = tmp
  END SUBROUTINE iswap


! *****************************************************************************
!> \brief Swaps two logicals
!> \param a Logicals to swap
!> \param b Logicals to swap
! *****************************************************************************
  ELEMENTAL SUBROUTINE bswap(a, b)
    LOGICAL, INTENT(INOUT)                   :: a, b

    LOGICAL                                  :: tmp

    tmp = a
    a = b
    b = tmp
  END SUBROUTINE bswap


! *****************************************************************************
!> \brief Splits an array of int8 values into two int4 arrays.
!> \param[in] merged          array of merged values
!> \param[out] array_upper    array to fill with the upper bytes of the merged
!>                            values
!> \param[out] array_lower    array to fill with the lower bytes of the merged
!>                            values
! *****************************************************************************
  SUBROUTINE dbcsr_unpack_i8_2i4 (merged, array_upper, array_lower)
    INTEGER(KIND=int_8), DIMENSION(:), &
      INTENT(IN)                             :: merged
    INTEGER(KIND=int_4), DIMENSION(:), &
      INTENT(OUT)                            :: array_upper, array_lower

    INTEGER(KIND=int_8), PARAMETER           :: lmask8 = 4294967295_int_8

    INTEGER                                  :: i

!
!   ---------------------------------------------------------------------------
! Lmask is used to filter in the lower 4 bytes and so its lower 32 bits are
! set to 1: lmask8 = 2^32-1.
! Umask is used to filter in the higher 4 bytes and so its higher 32 bits
! are set to 1: umask8 = 2^32-1 << 32
!lmask8 = 4294967295 ! 2^32-1
!umask8 = 18446744069414584320 ! (2^32-1) * 2^32 = (2^64-1)-(2^32-1)

    DO i = 1, SIZE (merged)
       array_upper(i) = INT (ISHFT (merged(i), -32), KIND=int_4)
       array_lower(i) = INT (IAND (merged(i), lmask8), KIND=int_4)
    ENDDO
  END SUBROUTINE dbcsr_unpack_i8_2i4


! *****************************************************************************
!> \brief ...
!> \param string ...
! *****************************************************************************
  SUBROUTINE uppercase(string)
    CHARACTER(LEN=*), INTENT(INOUT)          :: string

    INTEGER                                  :: i, iascii

    DO i=1,LEN_TRIM(string)
       iascii = ICHAR(string(i:i))
       IF ((iascii >= 97).AND.(iascii <= 122)) THEN
          string(i:i) = CHAR(iascii - 32)
       END IF
    END DO

  END SUBROUTINE uppercase


! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
! *****************************************************************************
  ELEMENTAL FUNCTION gcd(a,b)
    INTEGER, INTENT(IN)                      :: a, b
    INTEGER                                  :: gcd

    INTEGER                                  :: aa, ab, l, rem, s

    aa=ABS(a)
    ab=ABS(b)
    IF (aa<ab) THEN
       s=aa
       l=ab
    ELSE
       s=ab
       l=aa
    ENDIF
    IF (s.NE.0) THEN
       DO
          rem=MOD(l,s)
          IF (rem==0) EXIT
          l=s
          s=rem
       ENDDO
       GCD=s
    ELSE
       GCD=l
    ENDIF
  END FUNCTION gcd

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
! *****************************************************************************
  ELEMENTAL FUNCTION lcm(a,b)
    INTEGER, INTENT(IN)                      :: a, b
    INTEGER                                  :: lcm

    INTEGER                                  :: tmp

    tmp = gcd(a,b)
    IF (tmp==0) THEN
       lcm = 0
    ELSE
       ! could still overflow if the true lcm is larger than maxint
       lcm = ABS((a/tmp)*b)
    END IF
  END FUNCTION lcm

! *****************************************************************************
!> \brief generates the hash of a string and the index in the table
!> \param key a string of any length
!> \retval hash_index ...
!> \par History
!>       09.2006 created [Joost VandeVondele]
!>       2010-11-23 copied from CP2K and converted to take integers
!> \note
!>       http://en.wikipedia.org/wiki/Hash_table
!>       http://www.burtleburtle.net/bob/hash/doobs.html
!>       However, since fortran doesn't have an unsigned 4 byte int
!>       we compute it using an integer with the appropriate range
!>       we return already the index in the table as a final result
! *****************************************************************************
  FUNCTION joaat_hash(key) RESULT(hash_index)
    INTEGER, DIMENSION(:), INTENT(IN)        :: key
    INTEGER                                  :: hash_index

    INTEGER(KIND=int_8), PARAMETER           :: b32 = 2_int_8**32-1_int_8

    INTEGER                                  :: i, j
    INTEGER(KIND=int_8)                      :: byte, hash

    hash=0_int_8
    DO i=1,SIZE(key)
       DO j = 0, 3
          byte = IAND (ISHFT (key(i), -j*8), 255)
          hash=IAND(hash+byte                           ,b32)
          hash=IAND(     hash+IAND(ISHFT(hash,10),b32)  ,b32)
          hash=IAND(IEOR(hash,IAND(ISHFT(hash,-6),b32)) ,b32)
       ENDDO
    ENDDO
    hash=IAND(     hash+IAND(ISHFT(hash,  3),b32)  ,b32)
    hash=IAND(IEOR(hash,IAND(ISHFT(hash,-11),b32)) ,b32)
    hash=IAND(     hash+IAND(ISHFT(hash, 15),b32)  ,b32)
    ! In fortran 4-byte-integers have only 31 bits because they are signed
    ! In fortran the rightmost (least significant) bit is in position 0
    hash_index = IBCLR(hash, 31)
  END FUNCTION joaat_hash

END MODULE dbcsr_toollib
