!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE routine_map
   USE kinds, ONLY: default_string_length, int_4, int_8
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   #:include 'hash_map.fypp'

   #:call hash_map(prefix='routine_map', &
      key_type='CHARACTER(LEN=default_string_length)', &
      value_type='INTEGER(kind=int_4)')
   #:endcall hash_map

! **************************************************************************************************
! This is joaat_hash from string_table.F
!
!> \brief generates the hash of a given string
!> \param key a string of any length
!> \return ...
!> \par History
!>       09.2006 created [Joost VandeVondele]
!>       12.2012 copied and adopted [ole]
!> \note
!>       http://en.wikipedia.org/wiki/Hash_table
!>       http://www.burtleburtle.net/bob/hash/doobs.html
! **************************************************************************************************
   PURE FUNCTION routine_map_hash_function(key) RESULT(hash)
      CHARACTER(LEN=*), INTENT(IN)                       :: key
      INTEGER(KIND=int_8)                                :: hash

      INTEGER(KIND=int_8), PARAMETER                     :: b32 = 2_int_8**32 - 1_int_8

      INTEGER                                            :: i

      hash = 0_int_8
      DO i = 1, LEN(key)
         hash = IAND(hash + ICHAR(key(i:i)), b32)
         hash = IAND(hash + IAND(ISHFT(hash, 10), b32), b32)
         hash = IAND(IEOR(hash, IAND(ISHFT(hash, -6), b32)), b32)
      END DO
      hash = IAND(hash + IAND(ISHFT(hash, 3), b32), b32)
      hash = IAND(IEOR(hash, IAND(ISHFT(hash, -11), b32)), b32)
      hash = IAND(hash + IAND(ISHFT(hash, 15), b32), b32)
   END FUNCTION routine_map_hash_function

! **************************************************************************************************
!> \brief ...
!> \param key ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION routine_map_keys_equal(key1, key2) RESULT(res)
      CHARACTER(LEN=*), INTENT(IN)                      :: key1, key2
      LOGICAL                                            :: res

      res = (key1 == key2)
   END FUNCTION routine_map_keys_equal

END MODULE routine_map
