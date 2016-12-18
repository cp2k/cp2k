!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR tensor unit test.
!> \author Patrick Seewald
! **************************************************************************************************
PROGRAM dbcsr_tensor_unittest
   USE dbcsr_api,                       ONLY: dbcsr_finalize_lib,&
                                              dbcsr_init_lib
   USE dbcsr_tensor_test,               ONLY: test_tensor_formats
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_cart_create,&
                                              mp_comm_free,&
                                              mp_environ,&
                                              mp_world_finalize,&
                                              mp_world_init
#include "../base/base_uses.f90"

   INTEGER                            :: mp_comm, group, numnodes, mynode, io_unit
   INTEGER, DIMENSION(:), ALLOCATABLE :: pdims, myploc
   INTEGER                            :: ndims, nblks_alloc, nblks_1, nblks_2, nblks_3, nblks_4
   INTEGER, DIMENSION(:), ALLOCATABLE :: size_1, size_2, size_3, size_4
   INTEGER, DIMENSION(:), ALLOCATABLE :: blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4

   LOGICAL, PARAMETER                 :: verbose = .FALSE.

   CALL mp_world_init(mp_comm)
   CALL mp_environ(numnodes, mynode, mp_comm)

   ! set standard output parameters
   io_unit = 0
   IF (mynode .EQ. 0) io_unit = default_output_unit

   ! initialize libdbcsr
   CALL dbcsr_init_lib()

   ! Process grid
   ALLOCATE (pdims(2), myploc(2))
   pdims(:) = 0
   CALL mp_cart_create(mp_comm, 2, pdims, myploc, group)
   DEALLOCATE (pdims, myploc)

!--------------------------------------------------------------------------------------------------!
! Test 1: Testing matrix representations of tensor rank 2                                                  !
!--------------------------------------------------------------------------------------------------!
   ndims = 2

   ! Number of blocks in each dimension
   nblks_1 = 14
   nblks_2 = 21

   ! Block sizes in each dimension
   ALLOCATE (size_1(nblks_1), size_2(nblks_2))

   size_1(:) = [3, 5, 1, 23, 2, 3, 1, 6, 3, 8, 2, 3, 5, 1]
   size_2(:) = [4, 2, 5, 3, 1, 5, 13, 5, 2, 4, 5, 6, 7, 2, 3, 1, 2, 6, 9, 12, 21]

   ! Number of non-zero blocks
   nblks_alloc = 12
   ALLOCATE (blk_ind_1(nblks_alloc), blk_ind_2(nblks_alloc))

   ! Indices of non-zero blocks (s.t. index of ith block is [blk_ind_1(i), blk_ind_2(i), ...])
   blk_ind_1(:) = [1, 1,  1,  2, 4,  4,  7,  10, 10, 10, 10, 13] !&
   blk_ind_2(:) = [1, 3, 11, 15, 4, 17, 21,   6,  9, 13, 19,  7] !&

   ! Test tensor formats
   CALL test_tensor_formats(ndims, group, default_output_unit, verbose, &
                            blk_size_1=size_1, blk_size_2=size_2, &
                            blk_ind_1=blk_ind_1, blk_ind_2=blk_ind_2)

   DEALLOCATE (size_1, size_2)
   DEALLOCATE (blk_ind_1, blk_ind_2)

!--------------------------------------------------------------------------------------------------!
! Test 2: Testing matrix representations of tensor rank 3                                                  !
!--------------------------------------------------------------------------------------------------!
   ndims = 3

   ! Number of blocks in each dimension
   nblks_1 = 4
   nblks_2 = 6
   nblks_3 = 3

   ! Block sizes in each dimension
   ALLOCATE (size_1(nblks_1), size_2(nblks_2), size_3(nblks_3))

   size_1(:) = [3, 1, 5, 2]
   size_2(:) = [1, 2, 5, 3, 2, 4]
   size_3(:) = [4, 2, 10]

   ! Number of non-zero blocks
   nblks_alloc = 6
   ALLOCATE (blk_ind_1(nblks_alloc), blk_ind_2(nblks_alloc), blk_ind_3(nblks_alloc))

   ! Indices of non-zero blocks (s.t. index of ith block is [blk_ind_1(i), blk_ind_2(i), ...])
   blk_ind_1(:) = [1, 1, 1, 2, 2, 2] !&
   blk_ind_2(:) = [2, 2, 4, 1, 1, 2] !&
   blk_ind_3(:) = [1, 3, 3, 2, 3, 2] !&

   ! Test tensor formats
   CALL test_tensor_formats(ndims, group, default_output_unit, verbose, &
                            blk_size_1=size_1, blk_size_2=size_2, blk_size_3=size_3, &
                            blk_ind_1=blk_ind_1, blk_ind_2=blk_ind_2, blk_ind_3=blk_ind_3)

   DEALLOCATE (size_1, size_2, size_3)
   DEALLOCATE (blk_ind_1, blk_ind_2, blk_ind_3)

!--------------------------------------------------------------------------------------------------!
! Test 3: Testing matrix representations of tensor rank 4                                                  !
!--------------------------------------------------------------------------------------------------!
   ndims = 4

   ! Number of blocks in each dimension
   nblks_1 = 2
   nblks_2 = 13
   nblks_3 = 7
   nblks_4 = 3

   ! Block sizes in each dimension
   ALLOCATE (size_1(nblks_1), size_2(nblks_2), size_3(nblks_3), size_4(nblks_4))

   size_1(:) = [5, 9]
   size_2(:) = [6, 2, 5, 12, 3, 1, 7, 2, 5, 17, 9, 3, 4]
   size_3(:) = [2, 7, 3, 8, 5, 15, 1]
   size_4(:) = [12, 5, 3]

   ! Number of non-zero blocks
   nblks_alloc = 19
   ALLOCATE (blk_ind_1(nblks_alloc), blk_ind_2(nblks_alloc), blk_ind_3(nblks_alloc), blk_ind_4(nblks_alloc))

   ! Indices of non-zero blocks (s.t. index of ith block is [blk_ind_1(i), blk_ind_2(i), ...])
   blk_ind_1(:) = [1, 1, 1, 1, 1, 1,  1,  1,  1,  1,  1, 2, 2, 2, 2, 2, 2, 2,  2] !&
   blk_ind_2(:) = [2, 2, 3, 4, 7, 7, 10, 11, 11, 12, 12, 1, 1, 3, 5, 6, 6, 9, 12] !&
   blk_ind_3(:) = [1, 4, 6, 3, 1, 4,  2,  5,  7,  3,  3, 1, 4, 7, 6, 4, 5, 2,  3] !&
   blk_ind_4(:) = [3, 2, 3, 1, 1, 2,  1,  3,  2,  2,  3, 1, 3, 2, 1, 1, 3, 2,  2] !&

   ! Test tensor formats
   CALL test_tensor_formats(ndims, group, default_output_unit, verbose, &
                            blk_size_1=size_1, blk_size_2=size_2, blk_size_3=size_3, blk_size_4=size_4, &
                            blk_ind_1=blk_ind_1, blk_ind_2=blk_ind_2, blk_ind_3=blk_ind_3, blk_ind_4=blk_ind_4)

   DEALLOCATE (size_1, size_2, size_3, size_4)
   DEALLOCATE (blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4)

!--------------------------------------------------------------------------------------------------!
! End tests                                                                                        !
!--------------------------------------------------------------------------------------------------!

   CALL mp_comm_free(group)

   ! finalize libdbcsr
   CALL dbcsr_finalize_lib(mp_comm, io_unit)

   !
   ! finalize mpi
   CALL mp_world_finalize()

END PROGRAM
