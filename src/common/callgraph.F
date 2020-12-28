!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE callgraph
   USE kinds, ONLY: int_4, int_8
   USE timings_base_type, ONLY: call_stat_type
#include "base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

#:include 'hash_map.fypp'

#:call hash_map(prefix='callgraph', &
                key_type='INTEGER(kind=int_4), DIMENSION(2)', &
                value_type='TYPE(call_stat_type), POINTER')
#:endcall hash_map

! **************************************************************************************************
!> \brief ...
!> \param key ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION callgraph_hash_function(key) RESULT(hash)
      INTEGER(kind=int_4), DIMENSION(2), INTENT(in)      :: key
      INTEGER(KIND=int_8)                                :: hash

      INTEGER(kind=int_8)                                :: k1, k2
      k1 = key(1) ! cast to int_8
      k2 = key(2)
      hash = IOR(k1, ISHFT(k2, 32))
   END FUNCTION callgraph_hash_function

! **************************************************************************************************
!> \brief ...
!> \param key ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION callgraph_keys_equal(key1, key2) RESULT(res)
      INTEGER(kind=int_4), DIMENSION(2), INTENT(in)      :: key1, key2
      LOGICAL                                            :: res

      res = ALL(key1 == key2)
   END FUNCTION callgraph_keys_equal

END MODULE callgraph
