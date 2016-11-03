!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   DBCSR tensor unit test.
!> create the same 3-rank tensor with two different distributions and matrix representations,
!> and check for equality. This is an early tester, more complete testing will come soon.
!> \author  Patrick Seewald
! **************************************************************************************************
PROGRAM dbcsr_tensor_unittest
   USE dbcsr_api,                       ONLY: dbcsr_finalize_lib,&
                                              dbcsr_init_lib,&
                                              dbcsr_set_default_config,&
                                              dbcsr_type_real_8
   USE dbcsr_tensor,                    ONLY: &
        dbcsr_t_get_block, dbcsr_t_get_stored_coordinates, dbcsr_t_iterator, &
        dbcsr_t_iterator_blocks_left, dbcsr_t_iterator_next_block, dbcsr_t_iterator_start, &
        dbcsr_t_iterator_stop, dbcsr_t_put_block, dbcsr_t_reserve_blocks
   USE dbcsr_tensor_test,               ONLY: dbcsr_t_equal
   USE dbcsr_tensor_types,              ONLY: dbcsr_t_create,&
                                              dbcsr_t_destroy,&
                                              dbcsr_t_distribution_destroy,&
                                              dbcsr_t_distribution_new,&
                                              dbcsr_t_distribution_type,&
                                              dbcsr_t_type,&
                                              dims_tensor
   USE kinds,                           ONLY: real_8
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_cart_create,&
                                              mp_comm_free,&
                                              mp_environ,&
                                              mp_world_finalize,&
                                              mp_world_init
#include "../base/base_uses.f90"

   TYPE(dbcsr_t_distribution_type)             :: dist1, dist2
   TYPE(dbcsr_t_type)                          :: tensor1, tensor2
   INTEGER                                    :: mp_comm, group, numnodes, mynode, &
                                                 io_unit

   INTEGER, DIMENSION(:), ALLOCATABLE         :: dist1_1, dist1_2, dist1_3, dist2_1, dist2_2, dist2_3
   INTEGER, DIMENSION(3)                      :: pdims, myploc
   INTEGER                                    :: ndims, nblks_alloc, nblks_1, nblks_2, nblks_3, blk
   INTEGER, DIMENSION(:), ALLOCATABLE         :: size_1, size_2, size_3
   INTEGER, DIMENSION(:), ALLOCATABLE         :: blk_ind_1, blk_ind_2, blk_ind_3
   LOGICAL                                    :: eql

   ndims = 3

   CALL mp_world_init(mp_comm)
   CALL mp_environ(numnodes, mynode, mp_comm)

   pdims(:) = 0
   CALL mp_cart_create(mp_comm, ndims, pdims, myploc, group)

   ! initialize libdbcsr
   CALL dbcsr_init_lib()
   CALL dbcsr_set_default_config()

   ! set standard output parameters
   io_unit = 0
   IF (mynode .EQ. 0) io_unit = default_output_unit

   nblks_1 = 4
   nblks_2 = 6
   nblks_3 = 3
   nblks_alloc = 6

   ALLOCATE (size_1(nblks_1), size_2(nblks_2), size_3(nblks_3))
   ALLOCATE (blk_ind_1(nblks_alloc), blk_ind_2(nblks_alloc), blk_ind_3(nblks_alloc))
   blk_ind_1(:) = [1, 1, 1, 2, 2, 2]
   blk_ind_2(:) = [2, 2, 4, 1, 1, 2]
   blk_ind_3(:) = [1, 3, 3, 2, 3, 2]

   size_1(:) = [3, 1, 5, 2]
   size_2(:) = [1, 2, 5, 3, 2, 4]
   size_3(:) = [4, 2, 10]

   CALL random_dist(dist1_1, nblks_1, 1, pdims(1))
   CALL random_dist(dist1_2, nblks_2, 2, pdims(2))
   CALL random_dist(dist1_3, nblks_3, 3, pdims(3))
   CALL random_dist(dist2_1, nblks_1, 3, pdims(1))
   CALL random_dist(dist2_2, nblks_2, 2, pdims(2))
   CALL random_dist(dist2_3, nblks_3, 1, pdims(3))

   CALL dbcsr_t_distribution_new(dist1, group, [1], [2, 3], dist1_1, dist1_2, dist1_3)
   CALL dbcsr_t_distribution_new(dist2, group, [1, 2], [3], dist2_1, dist2_2, dist2_3)
   CALL dbcsr_t_create(tensor1, "(1|23)", dist1, [1], [2, 3], &
                       dbcsr_type_real_8, size_1, size_2, size_3)
   CALL dbcsr_t_create(tensor2, "(12|3)", dist2, [1, 2], [3], &
                       dbcsr_type_real_8, size_1, size_2, size_3)

   CALL allocate_and_fill_test_tensor(tensor1, blk_ind_1, blk_ind_2, blk_ind_3, mynode)
   CALL allocate_and_fill_test_tensor(tensor2, blk_ind_1, blk_ind_2, blk_ind_3, mynode)

   eql = dbcsr_t_equal(tensor1, tensor2)
   IF (.NOT. eql) CPABORT("Tensors are not the same ;-(")

   DEALLOCATE (dist1_1, dist1_2, dist1_3, dist2_1, dist2_2, dist2_3, size_1, size_2, size_3, blk_ind_1, blk_ind_2, blk_ind_3)

   CALL dbcsr_t_destroy(tensor1)
   CALL dbcsr_t_distribution_destroy(dist1)
   CALL dbcsr_t_destroy(tensor2)
   CALL dbcsr_t_distribution_destroy(dist2)

   ! finalize libdbcsr
   CALL dbcsr_finalize_lib(mp_comm, io_unit)

   ! free comm
   CALL mp_comm_free(group)

   !
   ! finalize mpi
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param blk_ind_1 ...
!> \param blk_ind_2 ...
!> \param blk_ind_3 ...
!> \param mynode ...
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_test_tensor(tensor, blk_ind_1, blk_ind_2, blk_ind_3, mynode)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, DIMENSION(:)                              :: blk_ind_1, blk_ind_2, blk_ind_3
      INTEGER, INTENT(IN)                                :: mynode

      INTEGER                                            :: i, ib, my_nblks_alloc, nblks_alloc, proc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_blk_ind_1, my_blk_ind_2, my_blk_ind_3
      INTEGER, DIMENSION(3)                              :: blk_index, blk_offset, blk_size, ind_nd, &
                                                            tensor_dims
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :, :) :: blk_values
      TYPE(dbcsr_t_iterator)                             :: iterator

      nblks_alloc = SIZE(blk_ind_1)

      my_nblks_alloc = 0
      DO ib = 1, nblks_alloc
         ind_nd(:) = [blk_ind_1(ib), blk_ind_2(ib), blk_ind_3(ib)]
         CALL dbcsr_t_get_stored_coordinates(tensor, ind_nd, proc)
         IF (proc == mynode) THEN
            my_nblks_alloc = my_nblks_alloc+1
         ENDIF
      ENDDO

      ALLOCATE (my_blk_ind_1(my_nblks_alloc))
      ALLOCATE (my_blk_ind_2(my_nblks_alloc))
      ALLOCATE (my_blk_ind_3(my_nblks_alloc))

      i = 0
      DO ib = 1, nblks_alloc
         ind_nd(:) = [blk_ind_1(ib), blk_ind_2(ib), blk_ind_3(ib)]
         CALL dbcsr_t_get_stored_coordinates(tensor, ind_nd, proc)
         IF (proc == mynode) THEN
            i = i+1
            my_blk_ind_1(i) = blk_ind_1(ib)
            my_blk_ind_2(i) = blk_ind_2(ib)
            my_blk_ind_3(i) = blk_ind_3(ib)
         ENDIF
      ENDDO

      CALL dbcsr_t_reserve_blocks(tensor, my_blk_ind_1, my_blk_ind_2, my_blk_ind_3)

      CALL dbcsr_t_iterator_start(iterator, tensor)
      DO WHILE (dbcsr_t_iterator_blocks_left(iterator))
         CALL dbcsr_t_iterator_next_block(iterator, blk_index, blk, blk_size=blk_size, blk_offset=blk_offset)
         ALLOCATE (blk_values(blk_size(1), blk_size(2), blk_size(3)))
         CALL dims_tensor(tensor, tensor_dims)
         CALL enumerate_block_elements(blk_values, blk_offset, tensor_dims)
         CALL dbcsr_t_put_block(tensor, blk_index, blk_size, blk_values)
         DEALLOCATE (blk_values)
      ENDDO
      CALL dbcsr_t_iterator_stop(iterator)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param blk_size ...
!> \param blk_index ...
!> \param blk_values ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE write_block(blk_size, blk_index, blk_values, unit_nr)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: blk_size, blk_index
      REAL(KIND=real_8), DIMENSION(blk_size(1), blk_size&
         (2), blk_size(3)), INTENT(IN)                   :: blk_values
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: i1, i2, i3

      WRITE (unit_nr, *) "writing block", blk_index

      DO i3 = 1, blk_size(3)
      DO i2 = 1, blk_size(2)
      DO i1 = 1, blk_size(1)
         WRITE (unit_nr, *) i1, i2, i3, blk_values(i1, i2, i3)
      ENDDO
      ENDDO
      ENDDO
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE write_tensor_blocks(tensor, unit_nr)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: blk
      INTEGER, DIMENSION(3)                              :: blk_index, blk_size
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :, :) :: blk_values
      TYPE(dbcsr_t_iterator)                             :: iterator

      CALL dbcsr_t_iterator_start(iterator, tensor)
      DO WHILE (dbcsr_t_iterator_blocks_left(iterator))
         CALL dbcsr_t_iterator_next_block(iterator, blk_index, blk, blk_size=blk_size)
         ALLOCATE (blk_values(blk_size(1), blk_size(2), blk_size(3)))
         CALL dbcsr_t_get_block(tensor, blk_index, blk_size, blk_values)
         CALL write_block(blk_size, blk_index, blk_values, unit_nr)
         DEALLOCATE (blk_values)
      ENDDO
      CALL dbcsr_t_iterator_stop(iterator)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param dist_array ...
!> \param dist_size ...
!> \param seed ...
!> \param nbins ...
! **************************************************************************************************
   SUBROUTINE random_dist(dist_array, dist_size, seed, nbins)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out)    :: dist_array
      INTEGER, INTENT(in)                                :: dist_size, seed, nbins

      INTEGER                                            :: i

      ALLOCATE (dist_array(dist_size))
      DO i = 1, dist_size
         dist_array(i) = MODULO(nbins-(i-1+seed), nbins)
      END DO

   END SUBROUTINE random_dist

! **************************************************************************************************
!> \brief ...
!> \param blk_values ...
!> \param blk_offset ...
!> \param tensor_size ...
! **************************************************************************************************
   SUBROUTINE enumerate_block_elements(blk_values, blk_offset, tensor_size)
      REAL(KIND=real_8), DIMENSION(:, :, :), INTENT(OUT) :: blk_values
      INTEGER, DIMENSION(3), INTENT(IN)                  :: blk_offset, tensor_size

      INTEGER                                            :: i, num_offset, num_val

      num_offset = ((blk_offset(1)-1)* &
                    tensor_size(1)+blk_offset(2)-1)* &
                   tensor_size(2)+blk_offset(3)-1+1

      num_val = PRODUCT(SHAPE(blk_values))
      blk_values(:, :, :) = RESHAPE((/(i, i=num_offset, num_offset+num_val, 1)/), SHAPE(blk_values))
   END SUBROUTINE
END PROGRAM
