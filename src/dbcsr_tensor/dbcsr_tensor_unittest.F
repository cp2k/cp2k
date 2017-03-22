!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR tensor unit test.
!> \author Patrick Seewald
! **************************************************************************************************
PROGRAM dbcsr_tensor_unittest
   USE dbcsr_api,                       ONLY: dbcsr_finalize_lib,&
                                              dbcsr_init_lib,&
                                              dbcsr_type_real_8
   USE dbcsr_tensor_test,               ONLY: dbcsr_t_contract_test,&
                                              dbcsr_t_random_dist,&
                                              dbcsr_t_setup_test_tensor,&
                                              dbcsr_t_test_formats
   USE dbcsr_tensor_types,              ONLY: dbcsr_t_create,&
                                              dbcsr_t_destroy,&
                                              dbcsr_t_distribution_destroy,&
                                              dbcsr_t_distribution_new,&
                                              dbcsr_t_distribution_type,&
                                              dbcsr_t_nd_mp_comm,&
                                              dbcsr_t_type
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_cart_create,&
                                              mp_comm_free,&
                                              mp_environ,&
                                              mp_world_finalize,&
                                              mp_world_init
#include "../base/base_uses.f90"

   IMPLICIT NONE

   INTEGER                            :: mp_comm, group, numnodes, mynode, io_unit
   INTEGER, DIMENSION(:), ALLOCATABLE :: pdims, myploc
   INTEGER                            :: ndims, nblks_alloc, nblks_1, nblks_2, nblks_3, nblks_4, &
                                         mp_comm_1, mp_comm_2, mp_comm_3, nblks_alloc_1, nblks_alloc_2
   INTEGER, DIMENSION(:), ALLOCATABLE :: size_1, size_2, size_3, size_4, dist1_1, dist1_2, dist1_3, &
                                         dist2_1, dist2_2, dist3_1, dist3_2, dist3_3
   INTEGER, DIMENSION(:), ALLOCATABLE :: blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_3_1, blk_ind_3_2, blk_ind_4, &
                                         myploc_1, myploc_2, myploc_3, pdims_1, pdims_2, pdims_3
   INTEGER, DIMENSION(:), ALLOCATABLE :: map11, map31, map12, map32, map21, map22

   LOGICAL, PARAMETER                 :: verbose = .FALSE.
   LOGICAL, DIMENSION(:), ALLOCATABLE :: periods_1, periods_2, periods_3
   TYPE(dbcsr_t_distribution_type)         :: dist1, dist2, dist3
   TYPE(dbcsr_t_type)            :: tensor_A, tensor_B, tensor_C

   LOGICAL, PARAMETER                 :: test_format = .TRUE.
   LOGICAL, PARAMETER                 :: test_contraction = .TRUE.

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

   IF (test_format) THEN
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
      CALL dbcsr_t_test_formats(ndims, group, default_output_unit, verbose, &
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
      CALL dbcsr_t_test_formats(ndims, group, default_output_unit, verbose, &
                                blk_size_1=size_1, blk_size_2=size_2, blk_size_3=size_3, &
                                blk_ind_1=blk_ind_1, blk_ind_2=blk_ind_2, blk_ind_3=blk_ind_3)

      DEALLOCATE (size_1, size_2, size_3)
      DEALLOCATE (blk_ind_1, blk_ind_2, blk_ind_3)

!--------------------------------------------------------------------------------------------------!
! Test 3: Testing matrix representations of tensor rank 4                                          !
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
      CALL dbcsr_t_test_formats(ndims, group, default_output_unit, verbose, &
                                blk_size_1=size_1, blk_size_2=size_2, blk_size_3=size_3, blk_size_4=size_4, &
                                blk_ind_1=blk_ind_1, blk_ind_2=blk_ind_2, blk_ind_3=blk_ind_3, blk_ind_4=blk_ind_4)

      DEALLOCATE (size_1, size_2, size_3, size_4)
      DEALLOCATE (blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4)

   ENDIF
   IF (test_contraction) THEN

!--------------------------------------------------------------------------------------------------!
! Preparations for tensor contraction tests                                                        !
!--------------------------------------------------------------------------------------------------!

      nblks_1 = 4
      nblks_2 = 11
      nblks_3 = 9
      nblks_4 = 5

      ! Block sizes in each dimension
      ALLOCATE (size_1(nblks_1), size_2(nblks_2), size_3(nblks_3), size_4(nblks_4))

      size_1(:) = [3, 9, 12, 1]
      size_2(:) = [4, 2, 3, 1, 9, 2, 32, 10, 5, 8, 7]
      size_3(:) = [7, 3, 8, 7, 9, 5, 10, 23, 2]
      size_4(:) = [8, 1, 4, 13, 6]

      nblks_alloc_1 = 32
      ALLOCATE (blk_ind_1(nblks_alloc_1), blk_ind_2(nblks_alloc_1), blk_ind_3_1(nblks_alloc_1))

      blk_ind_1(:) = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, & !&
                      1, 2, 2, 2, 2, 2, 2, 2, 3, 3, & !&
                      3, 3, 3, 3, 3, 3, 3, 4, 4, 4, & !&
                      4, 4] !&

      blk_ind_2(:) = [ 3, 5, 5, 5, 6,  6,  7,  8, 10, 11, & !&
                      11, 1, 1, 4, 7,  7,  9, 10 , 2,  2, & !&
                       5, 6, 8, 8, 9, 11, 11,  2 , 4,  5, & !&
                       5, 8] !&

      blk_ind_3_1(:) = [7, 3, 5, 9, 6, 8, 2, 8, 3, 2, & !&
                        3, 1, 4, 6, 2, 7, 5, 8, 3, 7, & !&
                        1, 4, 3, 7, 8, 5, 8, 9, 6, 1, & !&
                        2, 7] !&

      nblks_alloc_2 = 12
      ALLOCATE (blk_ind_3_2(nblks_alloc_2), blk_ind_4(nblks_alloc_2))

      blk_ind_3_2(:) = [1, 1, 2, 2, 2, 4, 4, 5, 5, 6, & !&
                        8, 8] !&
      blk_ind_4(:) = [2, 3, 2, 4, 5, 3, 5, 1, 3, 3, & !&
                      1, 4] !&

!--------------------------------------------------------------------------------------------------!
! Test 4: Testing tensor contraction (12|3)x(3|4)=(12|4)                                           !
!--------------------------------------------------------------------------------------------------!

      ALLOCATE (map11(2), map12(1), map21(1), map22(1), map31(2), map32(1))
      map11(:) = [1, 2]
      map12(:) = [3]
      map21(:) = [1]
      map22(:) = [2]
      map31(:) = [1, 2]
      map32(:) = [3]

      ALLOCATE (myploc_1(3), pdims_1(3), periods_1(3))
      ALLOCATE (myploc_2(2), pdims_2(2), periods_2(2))
      ALLOCATE (myploc_3(3), pdims_3(3), periods_3(3))

      mp_comm_1 = dbcsr_t_nd_mp_comm(group, map11, map12)
      CALL mp_environ(mp_comm_1, 3, pdims_1, myploc_1, periods_1)

      mp_comm_2 = dbcsr_t_nd_mp_comm(group, map21, map22)
      CALL mp_environ(mp_comm_2, 2, pdims_2, myploc_2, periods_2)

      mp_comm_3 = dbcsr_t_nd_mp_comm(group, map31, map32)
      CALL mp_environ(mp_comm_3, 3, pdims_3, myploc_3, periods_3)

      CALL dbcsr_t_random_dist(dist1_1, nblks_1, pdims_1(1), mp_comm_1)
      CALL dbcsr_t_random_dist(dist1_2, nblks_2, pdims_1(2), mp_comm_1)
      CALL dbcsr_t_random_dist(dist1_3, nblks_3, pdims_1(3), mp_comm_1)

      CALL dbcsr_t_random_dist(dist2_1, nblks_3, pdims_2(1), mp_comm_2)
      CALL dbcsr_t_random_dist(dist2_2, nblks_4, pdims_2(2), mp_comm_2)

      CALL dbcsr_t_random_dist(dist3_1, nblks_1, pdims_3(1), mp_comm_3)
      CALL dbcsr_t_random_dist(dist3_2, nblks_2, pdims_3(2), mp_comm_3)
      CALL dbcsr_t_random_dist(dist3_3, nblks_4, pdims_3(3), mp_comm_3)

      CALL dbcsr_t_distribution_new(dist1, mp_comm_1, group, map11, map12, dist1_1, dist1_2, dist1_3)
      CALL dbcsr_t_distribution_new(dist2, mp_comm_2, group, map21, map22, dist2_1, dist2_2)
      CALL dbcsr_t_distribution_new(dist3, mp_comm_3, group, map31, map32, dist3_1, dist3_2, dist3_3)

      CALL dbcsr_t_create(tensor_A, "(12|3)", dist1, map11, map12, dbcsr_type_real_8, &
                          size_1, size_2, size_3)
      CALL dbcsr_t_create(tensor_B, "(3|4)", dist2, map21, map22, dbcsr_type_real_8, &
                          size_3, size_4)
      CALL dbcsr_t_create(tensor_C, "(12|4)", dist3, map31, map32, dbcsr_type_real_8, &
                          size_1, size_2, size_4)

      CALL dbcsr_t_setup_test_tensor(tensor_A, mp_comm_1, blk_ind_1, blk_ind_2, blk_ind_3_1)
      CALL dbcsr_t_setup_test_tensor(tensor_B, mp_comm_2, blk_ind_3_2, blk_ind_4)

      CALL dbcsr_t_contract_test(tensor_A, tensor_B, tensor_C, &
                                 [3], [2, 1], &
                                 [1], [2], &
                                 [2, 1], [3], &
                                 default_output_unit)

      DEALLOCATE (myploc_1, myploc_2, myploc_3)
      DEALLOCATE (pdims_1, pdims_2, pdims_3)
      DEALLOCATE (periods_1, periods_2, periods_3)
      DEALLOCATE (dist1_1, dist1_2, dist1_3)
      DEALLOCATE (dist2_1, dist2_2)
      DEALLOCATE (dist3_1, dist3_2, dist3_3)
      DEALLOCATE (map11, map12, map21, map22, map31, map32)

      CALL dbcsr_t_destroy(tensor_A)
      CALL dbcsr_t_destroy(tensor_B)
      CALL dbcsr_t_destroy(tensor_C)
      CALL dbcsr_t_distribution_destroy(dist1)
      CALL dbcsr_t_distribution_destroy(dist2)
      CALL dbcsr_t_distribution_destroy(dist3)
      CALL mp_comm_free(mp_comm_1)
      CALL mp_comm_free(mp_comm_2)
      CALL mp_comm_free(mp_comm_3)

!--------------------------------------------------------------------------------------------------!
! Test 5: Testing tensor contraction (2|31)x(4|3)=(24|1)                                           !
!--------------------------------------------------------------------------------------------------!

      ALLOCATE (map11(1), map12(2), map21(1), map22(1), map31(2), map32(1))
      map11(:) = [2]
      map12(:) = [3, 1]
      map21(:) = [2]
      map22(:) = [1]
      map31(:) = [2, 3]
      map32(:) = [1]

      ALLOCATE (myploc_1(3), pdims_1(3), periods_1(3))
      ALLOCATE (myploc_2(2), pdims_2(2), periods_2(2))
      ALLOCATE (myploc_3(3), pdims_3(3), periods_3(3))

      mp_comm_1 = dbcsr_t_nd_mp_comm(group, map11, map12)
      CALL mp_environ(mp_comm_1, 3, pdims_1, myploc_1, periods_1)

      mp_comm_2 = dbcsr_t_nd_mp_comm(group, map21, map22)
      CALL mp_environ(mp_comm_2, 2, pdims_2, myploc_2, periods_2)

      mp_comm_3 = dbcsr_t_nd_mp_comm(group, map31, map32)
      CALL mp_environ(mp_comm_3, 3, pdims_3, myploc_3, periods_3)

      CALL dbcsr_t_random_dist(dist1_1, nblks_1, pdims_1(1), mp_comm_1)
      CALL dbcsr_t_random_dist(dist1_2, nblks_2, pdims_1(2), mp_comm_1)
      CALL dbcsr_t_random_dist(dist1_3, nblks_3, pdims_1(3), mp_comm_1)

      CALL dbcsr_t_random_dist(dist2_1, nblks_3, pdims_2(1), mp_comm_2)
      CALL dbcsr_t_random_dist(dist2_2, nblks_4, pdims_2(2), mp_comm_2)

      CALL dbcsr_t_random_dist(dist3_1, nblks_1, pdims_3(1), mp_comm_3)
      CALL dbcsr_t_random_dist(dist3_2, nblks_2, pdims_3(2), mp_comm_3)
      CALL dbcsr_t_random_dist(dist3_3, nblks_4, pdims_3(3), mp_comm_3)

      CALL dbcsr_t_distribution_new(dist1, mp_comm_1, group, map11, map12, dist1_1, dist1_2, dist1_3)
      CALL dbcsr_t_distribution_new(dist2, mp_comm_2, group, map21, map22, dist2_1, dist2_2)
      CALL dbcsr_t_distribution_new(dist3, mp_comm_3, group, map31, map32, dist3_1, dist3_2, dist3_3)

      CALL dbcsr_t_create(tensor_A, "(2|31)", dist1, map11, map12, dbcsr_type_real_8, &
                          size_1, size_2, size_3)
      CALL dbcsr_t_create(tensor_B, "(4|3)", dist2, map21, map22, dbcsr_type_real_8, &
                          size_3, size_4)
      CALL dbcsr_t_create(tensor_C, "(24|1)", dist3, map31, map32, dbcsr_type_real_8, &
                          size_1, size_2, size_4)

      CALL dbcsr_t_setup_test_tensor(tensor_A, mp_comm_1, blk_ind_1, blk_ind_2, blk_ind_3_1)
      CALL dbcsr_t_setup_test_tensor(tensor_B, mp_comm_2, blk_ind_3_2, blk_ind_4)

      CALL dbcsr_t_contract_test(tensor_A, tensor_B, tensor_C, &
                                 [3], [1, 2], &
                                 [1], [2], &
                                 [1, 2], [3], &
                                 default_output_unit)

      DEALLOCATE (myploc_1, myploc_2, myploc_3)
      DEALLOCATE (pdims_1, pdims_2, pdims_3)
      DEALLOCATE (periods_1, periods_2, periods_3)
      DEALLOCATE (dist1_1, dist1_2, dist1_3)
      DEALLOCATE (dist2_1, dist2_2)
      DEALLOCATE (dist3_1, dist3_2, dist3_3)
      DEALLOCATE (map11, map12, map21, map22, map31, map32)

      CALL dbcsr_t_destroy(tensor_A)
      CALL dbcsr_t_destroy(tensor_B)
      CALL dbcsr_t_destroy(tensor_C)
      CALL dbcsr_t_distribution_destroy(dist1)
      CALL dbcsr_t_distribution_destroy(dist2)
      CALL dbcsr_t_distribution_destroy(dist3)
      CALL mp_comm_free(mp_comm_1)
      CALL mp_comm_free(mp_comm_2)
      CALL mp_comm_free(mp_comm_3)

!--------------------------------------------------------------------------------------------------!
! Test 6: Testing tensor contraction (4|3)x(1|32)=(24|1)                                           !
!--------------------------------------------------------------------------------------------------!

      ALLOCATE (map11(1), map12(2), map21(1), map22(1), map31(2), map32(1))
      map11(:) = [1]
      map12(:) = [3, 2]
      map21(:) = [2]
      map22(:) = [1]
      map31(:) = [2, 3]
      map32(:) = [1]

      ALLOCATE (myploc_1(3), pdims_1(3), periods_1(3))
      ALLOCATE (myploc_2(2), pdims_2(2), periods_2(2))
      ALLOCATE (myploc_3(3), pdims_3(3), periods_3(3))

      mp_comm_1 = dbcsr_t_nd_mp_comm(group, map11, map12)
      CALL mp_environ(mp_comm_1, 3, pdims_1, myploc_1, periods_1)

      mp_comm_2 = dbcsr_t_nd_mp_comm(group, map21, map22)
      CALL mp_environ(mp_comm_2, 2, pdims_2, myploc_2, periods_2)

      mp_comm_3 = dbcsr_t_nd_mp_comm(group, map31, map32)
      CALL mp_environ(mp_comm_3, 3, pdims_3, myploc_3, periods_3)

      CALL dbcsr_t_random_dist(dist1_1, nblks_1, pdims_1(1), mp_comm_1)
      CALL dbcsr_t_random_dist(dist1_2, nblks_2, pdims_1(2), mp_comm_1)
      CALL dbcsr_t_random_dist(dist1_3, nblks_3, pdims_1(3), mp_comm_1)

      CALL dbcsr_t_random_dist(dist2_1, nblks_3, pdims_2(1), mp_comm_2)
      CALL dbcsr_t_random_dist(dist2_2, nblks_4, pdims_2(2), mp_comm_2)

      CALL dbcsr_t_random_dist(dist3_1, nblks_1, pdims_3(1), mp_comm_3)
      CALL dbcsr_t_random_dist(dist3_2, nblks_2, pdims_3(2), mp_comm_3)
      CALL dbcsr_t_random_dist(dist3_3, nblks_4, pdims_3(3), mp_comm_3)

      CALL dbcsr_t_distribution_new(dist1, mp_comm_1, group, map11, map12, dist1_1, dist1_2, dist1_3)
      CALL dbcsr_t_distribution_new(dist2, mp_comm_2, group, map21, map22, dist2_1, dist2_2)
      CALL dbcsr_t_distribution_new(dist3, mp_comm_3, group, map31, map32, dist3_1, dist3_2, dist3_3)

      CALL dbcsr_t_create(tensor_A, "(1|32)", dist1, map11, map12, dbcsr_type_real_8, &
                          size_1, size_2, size_3)
      CALL dbcsr_t_create(tensor_B, "(4|3)", dist2, map21, map22, dbcsr_type_real_8, &
                          size_3, size_4)
      CALL dbcsr_t_create(tensor_C, "(24|1)", dist3, map31, map32, dbcsr_type_real_8, &
                          size_1, size_2, size_4)

      CALL dbcsr_t_setup_test_tensor(tensor_A, mp_comm_1, blk_ind_1, blk_ind_2, blk_ind_3_1)
      CALL dbcsr_t_setup_test_tensor(tensor_B, mp_comm_2, blk_ind_3_2, blk_ind_4)

      CALL dbcsr_t_contract_test(tensor_B, tensor_A, tensor_C, &
                                 [1], [2], &
                                 [3], [1, 2], &
                                 [3], [1, 2], &
                                 default_output_unit)

      DEALLOCATE (myploc_1, myploc_2, myploc_3)
      DEALLOCATE (pdims_1, pdims_2, pdims_3)
      DEALLOCATE (periods_1, periods_2, periods_3)
      DEALLOCATE (dist1_1, dist1_2, dist1_3)
      DEALLOCATE (dist2_1, dist2_2)
      DEALLOCATE (dist3_1, dist3_2, dist3_3)
      DEALLOCATE (map11, map12, map21, map22, map31, map32)

      CALL dbcsr_t_destroy(tensor_A)
      CALL dbcsr_t_destroy(tensor_B)
      CALL dbcsr_t_destroy(tensor_C)
      CALL dbcsr_t_distribution_destroy(dist1)
      CALL dbcsr_t_distribution_destroy(dist2)
      CALL dbcsr_t_distribution_destroy(dist3)
      CALL mp_comm_free(mp_comm_1)
      CALL mp_comm_free(mp_comm_2)
      CALL mp_comm_free(mp_comm_3)

!--------------------------------------------------------------------------------------------------!
! Cleanup for tensor contraction tests                                                             !
!--------------------------------------------------------------------------------------------------!

      DEALLOCATE (size_1, size_2, size_3, size_4)
      DEALLOCATE (blk_ind_1, blk_ind_2, blk_ind_3_1)
      DEALLOCATE (blk_ind_3_2, blk_ind_4)

   ENDIF

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
