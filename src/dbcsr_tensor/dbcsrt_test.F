!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief General methods for testing DBCSR tensors.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsrt_test
   USE dbcsr_api,                       ONLY: dbcsr_type_complex_4,&
                                              dbcsr_type_complex_8,&
                                              dbcsr_type_real_4,&
                                              dbcsr_type_real_8
   USE dbcsrt,                          ONLY: &
        dbcsrt_copy, dbcsrt_get_block, dbcsrt_iterator, dbcsrt_iterator_blocks_left, &
        dbcsrt_iterator_next_block, dbcsrt_iterator_start, dbcsrt_iterator_stop, dbcsrt_ndims, &
        dbcsrt_reserve_blocks
   USE dbcsrt_block,                    ONLY: block_nd
   USE dbcsrt_types,                    ONLY: dbcsrt_create,&
                                              dbcsrt_destroy,&
                                              dbcsrt_type
   USE kinds,                           ONLY: real_4,&
                                              real_8
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsrt_test'

   PUBLIC :: &
      dbcsrt_equal

CONTAINS

! **************************************************************************************************
!> \brief check if two (arbitrarily mapped and distributed) tensors are equal.
!> \param tensor1 ...
!> \param tensor2 ...
!> \retval dbcsrt_equal ...
! **************************************************************************************************
   FUNCTION dbcsrt_equal(tensor1, tensor2)
      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor1, tensor2
      LOGICAL                                            :: dbcsrt_equal

      INTEGER                                            :: blk
      TYPE(dbcsrt_type)                                  :: tensor2_tmp
      TYPE(dbcsrt_iterator)                              :: iter
      TYPE(block_nd)                                     :: blk_data1, blk_data2
      INTEGER, DIMENSION(dbcsrt_ndims(tensor1))          :: blk_size, ind_nd

      ! create a copy of tensor2 that has exact same data format as tensor1
      CALL dbcsrt_create(tensor1, tensor2_tmp)
      CALL dbcsrt_reserve_blocks(tensor1, tensor2_tmp)
      CALL dbcsrt_copy(tensor2, tensor2_tmp)

      dbcsrt_equal = .TRUE.
      MARK_USED(tensor2)

      CALL dbcsrt_iterator_start(iter, tensor1)

      DO WHILE (dbcsrt_iterator_blocks_left(iter))
         CALL dbcsrt_iterator_next_block(iter, ind_nd, blk, blk_size=blk_size)
         CALL dbcsrt_get_block(tensor1, ind_nd, blk_data1)
         CALL dbcsrt_get_block(tensor2_tmp, ind_nd, blk_data2)
         IF (.NOT. blocks_equal(blk_data1, blk_data2)) THEN
            dbcsrt_equal = .FALSE.
         ENDIF
      ENDDO

      CALL dbcsrt_iterator_stop(iter)
      CALL dbcsrt_destroy(tensor2_tmp)
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param block1 ...
!> \param block2 ...
!> \retval blocks_equal ...
! **************************************************************************************************
   PURE FUNCTION blocks_equal(block1, block2)
      TYPE(block_nd), INTENT(IN)                         :: block1, block2
      LOGICAL                                            :: blocks_equal

      SELECT CASE (block1%data_type)
      CASE (dbcsr_type_real_4)
         blocks_equal = MAXVAL(ABS(block1%r_sp%blk-block2%r_sp%blk)) .LT. 1.0E-12_real_4
      CASE (dbcsr_type_real_8)
         blocks_equal = MAXVAL(ABS(block1%r_dp%blk-block2%r_dp%blk)) .LT. 1.0E-12_real_8
      CASE (dbcsr_type_complex_4)
         blocks_equal = MAXVAL(ABS(block1%c_sp%blk-block2%c_sp%blk)) .LT. 1.0E-12_real_4
      CASE (dbcsr_type_complex_8)
         blocks_equal = MAXVAL(ABS(block1%c_dp%blk-block2%c_dp%blk)) .LT. 1.0E-12_real_8
      END SELECT

   END FUNCTION

END MODULE dbcsrt_test
