!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Generated from dbcsr_tensor_test.fypp using Fypp.                                                !
! **DO NOT** modify this file, edit dbcsr_tensor_test.fypp instead.                                !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief General methods for testing DBCSR tensors.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor_test
   USE allocate_wrap,                   ONLY: allocate_any
   USE dbcsr_api,                       ONLY: dbcsr_type_complex_4,&
                                              dbcsr_type_complex_8,&
                                              dbcsr_type_real_4,&
                                              dbcsr_type_real_8
   USE dbcsr_tensor,                    ONLY: &
        dbcsr_t_copy, dbcsr_t_get_block, dbcsr_t_get_stored_coordinates, &
        dbcsr_t_iterator_blocks_left, dbcsr_t_iterator_next_block, dbcsr_t_iterator_start, &
        dbcsr_t_iterator_stop, dbcsr_t_iterator_type, dbcsr_t_ndims, dbcsr_t_put_block, &
        dbcsr_t_reserve_blocks
   USE dbcsr_tensor_block,              ONLY: block_nd
   USE dbcsr_tensor_index,              ONLY: combine_index
   USE dbcsr_tensor_types,              ONLY: dbcsr_t_create,&
                                              dbcsr_t_destroy,&
                                              dbcsr_t_distribution_destroy,&
                                              dbcsr_t_distribution_new,&
                                              dbcsr_t_distribution_type,&
                                              dbcsr_t_nd_mp_comm,&
                                              dbcsr_t_type,&
                                              dims_tensor
   USE kinds,                           ONLY: default_string_length,&
                                              real_4,&
                                              real_8
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_comm_free,&
                                              mp_environ
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_test'

   PUBLIC :: &
      test_tensor_formats

CONTAINS

! **************************************************************************************************
!> \brief check if two (arbitrarily mapped and distributed) tensors are equal.
!> \param tensor1 ...
!> \param tensor2 ...
!> \retval dbcsr_t_equal ...
! **************************************************************************************************
   FUNCTION dbcsr_t_equal(tensor1, tensor2)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor1, tensor2
      LOGICAL                                            :: dbcsr_t_equal

      INTEGER                                            :: blk
      TYPE(dbcsr_t_type)                                 :: tensor2_tmp
      TYPE(dbcsr_t_iterator_type)                        :: iter
      TYPE(block_nd)                                     :: blk_data1, blk_data2
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor1))         :: blk_size, ind_nd

      ! create a copy of tensor2 that has exact same data format as tensor1
      CALL dbcsr_t_create(tensor1, tensor2_tmp)
      CALL dbcsr_t_reserve_blocks(tensor1, tensor2_tmp)
      CALL dbcsr_t_copy(tensor2, tensor2_tmp)

      dbcsr_t_equal = .TRUE.

      CALL dbcsr_t_iterator_start(iter, tensor1)

      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind_nd, blk, blk_size=blk_size)
         CALL dbcsr_t_get_block(tensor1, ind_nd, blk_data1)
         CALL dbcsr_t_get_block(tensor2_tmp, ind_nd, blk_data2)

         IF (.NOT. blocks_equal(blk_data1, blk_data2)) THEN
            dbcsr_t_equal = .FALSE.
         ENDIF
      ENDDO

      CALL dbcsr_t_iterator_stop(iter)
      CALL dbcsr_t_destroy(tensor2_tmp)
   END FUNCTION

! **************************************************************************************************
!> \brief check if two blocks are equal
!> \param block1 ...
!> \param block2 ...
!> \retval blocks_equal ...
! **************************************************************************************************
   PURE FUNCTION blocks_equal(block1, block2)
      TYPE(block_nd), INTENT(IN)                         :: block1, block2
      LOGICAL                                            :: blocks_equal

      SELECT CASE (block1%data_type)
      CASE (dbcsr_type_real_8)
         blocks_equal = MAXVAL(ABS(block1%r_dp%blk-block2%r_dp%blk)) .LT. 1.0E-12_real_8
      CASE (dbcsr_type_real_4)
         blocks_equal = MAXVAL(ABS(block1%r_sp%blk-block2%r_sp%blk)) .LT. 1.0E-12_real_4
      CASE (dbcsr_type_complex_8)
         blocks_equal = MAXVAL(ABS(block1%c_dp%blk-block2%c_dp%blk)) .LT. 1.0E-12_real_8
      CASE (dbcsr_type_complex_4)
         blocks_equal = MAXVAL(ABS(block1%c_sp%blk-block2%c_sp%blk)) .LT. 1.0E-12_real_4
      END SELECT

   END FUNCTION

! **************************************************************************************************
!> \brief Compute factorial
!> \param n ...
!> \retval factorial ...
! **************************************************************************************************
   PURE FUNCTION factorial(n)
      INTEGER, INTENT(IN)                                :: n
      INTEGER                                            :: factorial

      INTEGER                                            :: k

      factorial = PRODUCT((/(k, k=1, n)/))
   END FUNCTION

! **************************************************************************************************
!> \brief Compute all permutations p of (1, 2, ..., n)
!> \param n ...
!> \param p ...
! **************************************************************************************************
   SUBROUTINE permute(n, p)
      INTEGER, INTENT(IN)                              :: n
      INTEGER                                          :: i, c
      INTEGER, DIMENSION(n)                            :: pp
      INTEGER, DIMENSION(n, factorial(n)), INTENT(OUT) :: p

      pp = [(i, i=1, n)]
      c = 1
      CALL perm(1)
   CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param i ...
! **************************************************************************************************
      RECURSIVE SUBROUTINE perm(i)
      INTEGER, INTENT(IN)                                :: i

      INTEGER                                            :: j, t

         IF (i == n) THEN
            p(:, c) = pp(:)
            c = c+1
         ELSE
            DO j = i, n
               t = pp(i)
               pp(i) = pp(j)
               pp(j) = t
               CALL perm(i+1)
               t = pp(i)
               pp(i) = pp(j)
               pp(j) = t
            END DO
         END IF
      END SUBROUTINE
   END SUBROUTINE

! **************************************************************************************************
!> \brief Test equivalence of all tensor formats, using a random distribution.
!> \param ndims tensor rank
!> \param mp_comm_2d communicator
!> \param unit_nr output unit, needs to be a valid unit number on all mpi ranks
!> \param verbose if .TRUE., print all tensor blocks
!> \param blk_size_1 block sizes along dimension 1
!> \param blk_size_2 ...
!> \param blk_size_3 ...
!> \param blk_size_4 ...
!> \param blk_ind_1 index along dimension 1 of non-zero blocks
!> \param blk_ind_2 ...
!> \param blk_ind_3 ...
!> \param blk_ind_4 ...
! **************************************************************************************************
   SUBROUTINE test_tensor_formats(ndims, mp_comm_2d, unit_nr, verbose, &
                                  blk_size_1, blk_size_2, blk_size_3, blk_size_4, &
                                  blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4)
      INTEGER, INTENT(IN)                                :: ndims, mp_comm_2d, unit_nr
      LOGICAL, INTENT(IN)                                :: verbose
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: blk_size_1, blk_size_2, blk_size_3, &
                                                            blk_size_4, blk_ind_1, blk_ind_2, &
                                                            blk_ind_3, blk_ind_4

      CHARACTER(LEN=default_string_length)               :: tensor_name
      INTEGER                                            :: iblk, icount, idist, imap, io_unit, &
                                                            iperm, isep, mp_comm_ref, &
                                                            mp_comm_test, mynode, nblks, numnodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1_1, dist1_2, dist1_3, dist1_4, &
                                                            dist2_1, dist2_2, dist2_3, dist2_4, &
                                                            map1, map1_ref, map2, map2_ref
      INTEGER, DIMENSION(ndims, factorial(ndims))        :: perm
      INTEGER, DIMENSION(ndims)                          :: myploc, pdims_ref, pdims_test
      LOGICAL                                            :: eql
      LOGICAL, DIMENSION(ndims)                          :: periods
      TYPE(dbcsr_t_distribution_type)                    :: dist1, dist2
      TYPE(dbcsr_t_type)                                 :: tensor1, tensor2

      CALL mp_environ(numnodes, mynode, mp_comm_2d)

      io_unit = 0
      IF (mynode .EQ. 0) io_unit = unit_nr

      CALL permute(ndims, perm)
      ALLOCATE (map1_ref(ndims/2), map2_ref(ndims-ndims/2))
      map1_ref(:) = perm(1:ndims/2, 1)
      map2_ref(:) = perm(ndims/2+1:ndims, 1)

      mp_comm_ref = dbcsr_t_nd_mp_comm(mp_comm_2d, map1_ref, map2_ref)
      CALL mp_environ(numnodes, mynode, mp_comm_ref)
      CALL mp_environ(mp_comm_ref, ndims, pdims_ref, myploc, periods)

      IF (io_unit > 0) THEN
         WRITE (io_unit, *)
         WRITE (io_unit, '(A)') REPEAT("-", 80)
         WRITE (io_unit, '(A,1X,I1)') "Testing matrix representations of tensor rank", ndims
         WRITE (io_unit, '(A)') REPEAT("-", 80)
         WRITE (io_unit, '(A)') "Block sizes:"

         IF (ndims >= 1) THEN
            WRITE (io_unit, '(T4,A,1X,I1,A,1X)', advance='no') 'Dim', 1, ':'
            DO iblk = 1, SIZE(blk_size_1)
               WRITE (io_unit, '(I2,1X)', advance='no') blk_size_1(iblk)
            ENDDO
            WRITE (io_unit, *)
         ENDIF
         IF (ndims >= 2) THEN
            WRITE (io_unit, '(T4,A,1X,I1,A,1X)', advance='no') 'Dim', 2, ':'
            DO iblk = 1, SIZE(blk_size_2)
               WRITE (io_unit, '(I2,1X)', advance='no') blk_size_2(iblk)
            ENDDO
            WRITE (io_unit, *)
         ENDIF
         IF (ndims >= 3) THEN
            WRITE (io_unit, '(T4,A,1X,I1,A,1X)', advance='no') 'Dim', 3, ':'
            DO iblk = 1, SIZE(blk_size_3)
               WRITE (io_unit, '(I2,1X)', advance='no') blk_size_3(iblk)
            ENDDO
            WRITE (io_unit, *)
         ENDIF
         IF (ndims >= 4) THEN
            WRITE (io_unit, '(T4,A,1X,I1,A,1X)', advance='no') 'Dim', 4, ':'
            DO iblk = 1, SIZE(blk_size_4)
               WRITE (io_unit, '(I2,1X)', advance='no') blk_size_4(iblk)
            ENDDO
            WRITE (io_unit, *)
         ENDIF

         WRITE (io_unit, '(A)') "Non-zero blocks:"
         DO iblk = 1, SIZE(blk_ind_1)
            IF (ndims == 2) THEN
               WRITE (io_unit, '(T4,A, I3, A, 2I3, 1X, A)') &
                  'Block', iblk, ': (', blk_ind_1(iblk), blk_ind_2(iblk), ')'
            ENDIF
            IF (ndims == 3) THEN
               WRITE (io_unit, '(T4,A, I3, A, 3I3, 1X, A)') &
                  'Block', iblk, ': (', blk_ind_1(iblk), blk_ind_2(iblk), blk_ind_3(iblk), ')'
            ENDIF
            IF (ndims == 4) THEN
               WRITE (io_unit, '(T4,A, I3, A, 4I3, 1X, A)') &
                  'Block', iblk, ': (', blk_ind_1(iblk), blk_ind_2(iblk), blk_ind_3(iblk), blk_ind_4(iblk), ')'
            ENDIF
         ENDDO

         WRITE (io_unit, *)
         WRITE (io_unit, '(A,1X)', advance='no') "Reference map:"
         WRITE (io_unit, '(A1,1X)', advance='no') "("
         DO imap = 1, SIZE(map1_ref)
            WRITE (io_unit, '(I1,1X)', advance='no') map1_ref(imap)
         ENDDO
         WRITE (io_unit, '(A1,1X)', advance='no') "|"
         DO imap = 1, SIZE(map2_ref)
            WRITE (io_unit, '(I1,1X)', advance='no') map2_ref(imap)
         ENDDO
         WRITE (io_unit, '(A1)') ")"

      ENDIF

      icount = 0
      DO iperm = 1, factorial(ndims)
         DO isep = 1, ndims-1
            icount = icount+1

            ALLOCATE (map1(isep), map2(ndims-isep))
            map1(:) = perm(1:isep, iperm)
            map2(:) = perm(isep+1:ndims, iperm)

            mp_comm_test = dbcsr_t_nd_mp_comm(mp_comm_2d, map1, map2)
            CALL mp_environ(numnodes, mynode, mp_comm_test)
            CALL mp_environ(mp_comm_test, ndims, pdims_test, myploc, periods)

            IF (1 <= ndims) THEN
               nblks = SIZE(blk_size_1)
               CALL random_dist(dist1_1, nblks, pdims_test(1), mp_comm_test)
               CALL random_dist(dist2_1, nblks, pdims_ref(1), mp_comm_ref)
            ENDIF
            IF (2 <= ndims) THEN
               nblks = SIZE(blk_size_2)
               CALL random_dist(dist1_2, nblks, pdims_test(2), mp_comm_test)
               CALL random_dist(dist2_2, nblks, pdims_ref(2), mp_comm_ref)
            ENDIF
            IF (3 <= ndims) THEN
               nblks = SIZE(blk_size_3)
               CALL random_dist(dist1_3, nblks, pdims_test(3), mp_comm_test)
               CALL random_dist(dist2_3, nblks, pdims_ref(3), mp_comm_ref)
            ENDIF
            IF (4 <= ndims) THEN
               nblks = SIZE(blk_size_4)
               CALL random_dist(dist1_4, nblks, pdims_test(4), mp_comm_test)
               CALL random_dist(dist2_4, nblks, pdims_ref(4), mp_comm_ref)
            ENDIF

            IF (io_unit > 0) THEN
               WRITE (io_unit, *)
               WRITE (tensor_name, '(A,1X,I3,1X)') "Test", icount
               WRITE (io_unit, '(A,A,1X)', advance='no') TRIM(tensor_name), ':'
               WRITE (io_unit, '(A1,1X)', advance='no') "("
               DO imap = 1, SIZE(map1)
                  WRITE (io_unit, '(I1,1X)', advance='no') map1(imap)
               ENDDO
               WRITE (io_unit, '(A1,1X)', advance='no') "|"
               DO imap = 1, SIZE(map2)
                  WRITE (io_unit, '(I1,1X)', advance='no') map2(imap)
               ENDDO
               WRITE (io_unit, '(A1)') ")"

               WRITE (io_unit, '(T4,A)') "Reference distribution:"
               IF (1 <= ndims) THEN
                  WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec 1:"
                  DO idist = 1, SIZE(dist2_1)
                     WRITE (io_unit, '(I2,1X)', advance='no') dist2_1(idist)
                  ENDDO
                  WRITE (io_unit, *)
               ENDIF
               IF (2 <= ndims) THEN
                  WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec 2:"
                  DO idist = 1, SIZE(dist2_2)
                     WRITE (io_unit, '(I2,1X)', advance='no') dist2_2(idist)
                  ENDDO
                  WRITE (io_unit, *)
               ENDIF
               IF (3 <= ndims) THEN
                  WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec 3:"
                  DO idist = 1, SIZE(dist2_3)
                     WRITE (io_unit, '(I2,1X)', advance='no') dist2_3(idist)
                  ENDDO
                  WRITE (io_unit, *)
               ENDIF
               IF (4 <= ndims) THEN
                  WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec 4:"
                  DO idist = 1, SIZE(dist2_4)
                     WRITE (io_unit, '(I2,1X)', advance='no') dist2_4(idist)
                  ENDDO
                  WRITE (io_unit, *)
               ENDIF

               WRITE (io_unit, '(T4,A)') "Test distribution:"
               IF (1 <= ndims) THEN
                  WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec 1:"
                  DO idist = 1, SIZE(dist2_1)
                     WRITE (io_unit, '(I2,1X)', advance='no') dist1_1(idist)
                  ENDDO
                  WRITE (io_unit, *)
               ENDIF
               IF (2 <= ndims) THEN
                  WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec 2:"
                  DO idist = 1, SIZE(dist2_2)
                     WRITE (io_unit, '(I2,1X)', advance='no') dist1_2(idist)
                  ENDDO
                  WRITE (io_unit, *)
               ENDIF
               IF (3 <= ndims) THEN
                  WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec 3:"
                  DO idist = 1, SIZE(dist2_3)
                     WRITE (io_unit, '(I2,1X)', advance='no') dist1_3(idist)
                  ENDDO
                  WRITE (io_unit, *)
               ENDIF
               IF (4 <= ndims) THEN
                  WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec 4:"
                  DO idist = 1, SIZE(dist2_4)
                     WRITE (io_unit, '(I2,1X)', advance='no') dist1_4(idist)
                  ENDDO
                  WRITE (io_unit, *)
               ENDIF
            ENDIF

            IF (ndims == 2) THEN
               CALL dbcsr_t_distribution_new(dist2, mp_comm_ref, mp_comm_2d, map1_ref, map2_ref, dist2_1, dist2_2)
               CALL dbcsr_t_create(tensor2, "Ref", dist2, map1_ref, map2_ref, &
                                   dbcsr_type_real_8, blk_size_1, blk_size_2)
               CALL allocate_and_fill_test_tensor(tensor2, mp_comm_ref, blk_ind_1, blk_ind_2)
            ENDIF
            IF (ndims == 3) THEN
               CALL dbcsr_t_distribution_new(dist2, mp_comm_ref, mp_comm_2d, map1_ref, map2_ref, dist2_1, dist2_2, dist2_3)
               CALL dbcsr_t_create(tensor2, "Ref", dist2, map1_ref, map2_ref, &
                                   dbcsr_type_real_8, blk_size_1, blk_size_2, blk_size_3)
               CALL allocate_and_fill_test_tensor(tensor2, mp_comm_ref, blk_ind_1, blk_ind_2, blk_ind_3)
            ENDIF
            IF (ndims == 4) THEN
               CALL dbcsr_t_distribution_new(dist2, mp_comm_ref, mp_comm_2d, map1_ref, map2_ref, dist2_1, dist2_2, dist2_3, dist2_4)
               CALL dbcsr_t_create(tensor2, "Ref", dist2, map1_ref, map2_ref, &
                                   dbcsr_type_real_8, blk_size_1, blk_size_2, blk_size_3, blk_size_4)
               CALL allocate_and_fill_test_tensor(tensor2, mp_comm_ref, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4)
            ENDIF

            IF (verbose) CALL write_tensor_blocks(tensor2, io_unit, unit_nr)

            IF (ndims == 2) THEN
               CALL dbcsr_t_distribution_new(dist1, mp_comm_test, mp_comm_2d, map1, map2, dist1_1, dist1_2)
               CALL dbcsr_t_create(tensor1, tensor_name, dist1, map1, map2, &
                                   dbcsr_type_real_8, blk_size_1, blk_size_2)
               CALL allocate_and_fill_test_tensor(tensor1, mp_comm_test, blk_ind_1, blk_ind_2)
            ENDIF
            IF (ndims == 3) THEN
               CALL dbcsr_t_distribution_new(dist1, mp_comm_test, mp_comm_2d, map1, map2, dist1_1, dist1_2, dist1_3)
               CALL dbcsr_t_create(tensor1, tensor_name, dist1, map1, map2, &
                                   dbcsr_type_real_8, blk_size_1, blk_size_2, blk_size_3)
               CALL allocate_and_fill_test_tensor(tensor1, mp_comm_test, blk_ind_1, blk_ind_2, blk_ind_3)
            ENDIF
            IF (ndims == 4) THEN
               CALL dbcsr_t_distribution_new(dist1, mp_comm_test, mp_comm_2d, map1, map2, dist1_1, dist1_2, dist1_3, dist1_4)
               CALL dbcsr_t_create(tensor1, tensor_name, dist1, map1, map2, &
                                   dbcsr_type_real_8, blk_size_1, blk_size_2, blk_size_3, blk_size_4)
               CALL allocate_and_fill_test_tensor(tensor1, mp_comm_test, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4)
            ENDIF

            IF (verbose) CALL write_tensor_blocks(tensor1, io_unit, unit_nr)

            eql = dbcsr_t_equal(tensor1, tensor2)

            IF (.NOT. eql) THEN
               IF (io_unit > 0) WRITE (io_unit, '(A,1X,A)') TRIM(tensor_name), 'failed!'
               CPABORT('')
            ELSE
               IF (io_unit > 0) WRITE (io_unit, '(A,1X,A)') TRIM(tensor_name), 'passed!'
            ENDIF
            DEALLOCATE (map1, map2)

            CALL dbcsr_t_destroy(tensor1)
            CALL dbcsr_t_distribution_destroy(dist1)

            CALL mp_comm_free(mp_comm_test)
            CALL dbcsr_t_destroy(tensor2)
            CALL dbcsr_t_distribution_destroy(dist2)

         ENDDO
      ENDDO
      CALL mp_comm_free(mp_comm_ref)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create test distribution
!> \param dist_array ...
!> \param dist_size ...
!> \param nbins ...
!> \param mp_comm ...
! **************************************************************************************************
   SUBROUTINE random_dist(dist_array, dist_size, nbins, mp_comm)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out)    :: dist_array
      INTEGER, INTENT(IN)                                :: dist_size, nbins, mp_comm

      INTEGER                                            :: mynode, numnodes
      REAL, DIMENSION(dist_size)                         :: rn

      CALL mp_environ(numnodes, mynode, mp_comm)

      IF (mynode .EQ. 0) THEN
         CALL RANDOM_NUMBER(rn)
      ENDIF
      CALL mp_bcast(rn, 0, mp_comm)

      ALLOCATE (dist_array(dist_size))

      dist_array(:) = FLOOR(rn*nbins)

   END SUBROUTINE random_dist

! **************************************************************************************************
!> \brief Allocate and fill test tensor - entries are enumerated by their index s.t. they only depend
!>        on global properties of the tensor but not on distribution, matrix representation, etc.
!> \param tensor ...
!> \param mp_comm communicator
!> \param blk_ind_1 index along dimension 1 of non-zero blocks
!> \param blk_ind_2 ...
!> \param blk_ind_3 ...
!> \param blk_ind_4 ...
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_test_tensor(tensor, mp_comm, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: blk_ind_1, blk_ind_2, blk_ind_3, &
                                                            blk_ind_4

      INTEGER                                            :: blk, i, ib, my_nblks_alloc, mynode, &
                                                            nblks_alloc, numnodes, proc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_blk_ind_1, my_blk_ind_2, &
                                                            my_blk_ind_3, my_blk_ind_4
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ind_nd
      TYPE(dbcsr_t_iterator_type)                        :: iterator
      REAL(KIND=real_8), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: blk_values_4
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :, :) :: blk_values_3
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :)    :: blk_values_2
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor))          :: blk_index, blk_offset, blk_size, &
                                                            tensor_dims

      nblks_alloc = SIZE(blk_ind_1)
      CALL mp_environ(numnodes, mynode, mp_comm)

      ALLOCATE (ind_nd(nblks_alloc, dbcsr_t_ndims(tensor)))
      my_nblks_alloc = 0
      DO ib = 1, nblks_alloc
         IF (dbcsr_t_ndims(tensor) == 2) THEN
            ind_nd(ib, :) = [blk_ind_1(ib), blk_ind_2(ib)]
         ENDIF
         IF (dbcsr_t_ndims(tensor) == 3) THEN
            ind_nd(ib, :) = [blk_ind_1(ib), blk_ind_2(ib), blk_ind_3(ib)]
         ENDIF
         IF (dbcsr_t_ndims(tensor) == 4) THEN
            ind_nd(ib, :) = [blk_ind_1(ib), blk_ind_2(ib), blk_ind_3(ib), blk_ind_4(ib)]
         ENDIF
         CALL dbcsr_t_get_stored_coordinates(tensor, ind_nd(ib, :), proc)
         IF (proc == mynode) THEN
            my_nblks_alloc = my_nblks_alloc+1
         ENDIF
      ENDDO

      IF (dbcsr_t_ndims(tensor) >= 1) THEN
         ALLOCATE (my_blk_ind_1(my_nblks_alloc))
      ENDIF
      IF (dbcsr_t_ndims(tensor) >= 2) THEN
         ALLOCATE (my_blk_ind_2(my_nblks_alloc))
      ENDIF
      IF (dbcsr_t_ndims(tensor) >= 3) THEN
         ALLOCATE (my_blk_ind_3(my_nblks_alloc))
      ENDIF
      IF (dbcsr_t_ndims(tensor) >= 4) THEN
         ALLOCATE (my_blk_ind_4(my_nblks_alloc))
      ENDIF

      i = 0
      DO ib = 1, nblks_alloc
         CALL dbcsr_t_get_stored_coordinates(tensor, ind_nd(ib, :), proc)
         IF (proc == mynode) THEN
            i = i+1
            IF (dbcsr_t_ndims(tensor) >= 1) THEN
               my_blk_ind_1(i) = blk_ind_1(ib)
            ENDIF
            IF (dbcsr_t_ndims(tensor) >= 2) THEN
               my_blk_ind_2(i) = blk_ind_2(ib)
            ENDIF
            IF (dbcsr_t_ndims(tensor) >= 3) THEN
               my_blk_ind_3(i) = blk_ind_3(ib)
            ENDIF
            IF (dbcsr_t_ndims(tensor) >= 4) THEN
               my_blk_ind_4(i) = blk_ind_4(ib)
            ENDIF
         ENDIF
      ENDDO

      IF (dbcsr_t_ndims(tensor) == 2) THEN
         CALL dbcsr_t_reserve_blocks(tensor, my_blk_ind_1, my_blk_ind_2)
      ENDIF
      IF (dbcsr_t_ndims(tensor) == 3) THEN
         CALL dbcsr_t_reserve_blocks(tensor, my_blk_ind_1, my_blk_ind_2, my_blk_ind_3)
      ENDIF
      IF (dbcsr_t_ndims(tensor) == 4) THEN
         CALL dbcsr_t_reserve_blocks(tensor, my_blk_ind_1, my_blk_ind_2, my_blk_ind_3, my_blk_ind_4)
      ENDIF

      CALL dbcsr_t_iterator_start(iterator, tensor)
      DO WHILE (dbcsr_t_iterator_blocks_left(iterator))
         CALL dbcsr_t_iterator_next_block(iterator, blk_index, blk, blk_size=blk_size, blk_offset=blk_offset)

         IF (dbcsr_t_ndims(tensor) == 2) THEN
            CALL allocate_any(blk_size, blk_values_2)
            CALL dims_tensor(tensor, tensor_dims)
            CALL enumerate_block_elements(blk_size, blk_offset, tensor_dims, blk_2=blk_values_2)
            CALL dbcsr_t_put_block(tensor, blk_index, blk_size, blk_values_2)
            DEALLOCATE (blk_values_2)
         ENDIF
         IF (dbcsr_t_ndims(tensor) == 3) THEN
            CALL allocate_any(blk_size, blk_values_3)
            CALL dims_tensor(tensor, tensor_dims)
            CALL enumerate_block_elements(blk_size, blk_offset, tensor_dims, blk_3=blk_values_3)
            CALL dbcsr_t_put_block(tensor, blk_index, blk_size, blk_values_3)
            DEALLOCATE (blk_values_3)
         ENDIF
         IF (dbcsr_t_ndims(tensor) == 4) THEN
            CALL allocate_any(blk_size, blk_values_4)
            CALL dims_tensor(tensor, tensor_dims)
            CALL enumerate_block_elements(blk_size, blk_offset, tensor_dims, blk_4=blk_values_4)
            CALL dbcsr_t_put_block(tensor, blk_index, blk_size, blk_values_4)
            DEALLOCATE (blk_values_4)
         ENDIF
      ENDDO
      CALL dbcsr_t_iterator_stop(iterator)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Enumerate tensor entries in block
!> \param blk_size size of block
!> \param blk_offset block offset (indices of first element)
!> \param tensor_size global tensor sizes
!> \param blk_2 ...
!> \param blk_3 ...
!> \param blk_4 ...
! **************************************************************************************************
   SUBROUTINE enumerate_block_elements(blk_size, blk_offset, tensor_size, blk_2, blk_3, blk_4)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size, blk_offset, tensor_size
      REAL(KIND=real_8), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: blk_2
      REAL(KIND=real_8), DIMENSION(:, :, :), &
         INTENT(OUT), OPTIONAL                           :: blk_3
      REAL(KIND=real_8), DIMENSION(:, :, :, :), &
         INTENT(OUT), OPTIONAL                           :: blk_4

      INTEGER                                            :: i_1, i_2, i_3, i_4, ndim
      INTEGER, DIMENSION(SIZE(blk_size))                 :: arr_ind, tens_ind

      ndim = SIZE(tensor_size)

      IF (ndim == 2) THEN
         DO i_2 = 1, blk_size(2)
         DO i_1 = 1, blk_size(1)
            arr_ind(:) = [i_1, i_2]
            tens_ind(:) = arr_ind(:)+blk_offset(:)-1
            blk_2(arr_ind(1), arr_ind(2)) = combine_index(tens_ind, tensor_size)
         ENDDO
         ENDDO
      ENDIF
      IF (ndim == 3) THEN
         DO i_3 = 1, blk_size(3)
         DO i_2 = 1, blk_size(2)
         DO i_1 = 1, blk_size(1)
            arr_ind(:) = [i_1, i_2, i_3]
            tens_ind(:) = arr_ind(:)+blk_offset(:)-1
            blk_3(arr_ind(1), arr_ind(2), arr_ind(3)) = combine_index(tens_ind, tensor_size)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
      IF (ndim == 4) THEN
         DO i_4 = 1, blk_size(4)
         DO i_3 = 1, blk_size(3)
         DO i_2 = 1, blk_size(2)
         DO i_1 = 1, blk_size(1)
            arr_ind(:) = [i_1, i_2, i_3, i_4]
            tens_ind(:) = arr_ind(:)+blk_offset(:)-1
            blk_4(arr_ind(1), arr_ind(2), arr_ind(3), arr_ind(4)) = combine_index(tens_ind, tensor_size)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Write all tensor blocks
!> \param tensor ...
!> \param io_unit_master for global output
!> \param io_unit_all for local output
! **************************************************************************************************
   SUBROUTINE write_tensor_blocks(tensor, io_unit_master, io_unit_all)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, INTENT(IN)                                :: io_unit_master, io_unit_all

      INTEGER                                            :: blk, proc
      TYPE(dbcsr_t_iterator_type)                        :: iterator
      REAL(KIND=real_8), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: blk_values_4
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :, :) :: blk_values_3
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :)    :: blk_values_2
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor))          :: blk_index, blk_size

      IF (io_unit_master > 0) THEN
         WRITE (io_unit_master, '(T7,A)') "(block index) @ process: (array index) value"
      ENDIF
      CALL dbcsr_t_iterator_start(iterator, tensor)
      DO WHILE (dbcsr_t_iterator_blocks_left(iterator))
         CALL dbcsr_t_iterator_next_block(iterator, blk_index, blk, blk_size=blk_size)
         CALL dbcsr_t_get_stored_coordinates(tensor, blk_index, proc)
         IF (dbcsr_t_ndims(tensor) == 2) THEN
            CALL allocate_any(blk_size, blk_values_2)
            CALL dbcsr_t_get_block(tensor, blk_index, blk_size, blk_values_2)
            CALL write_block(tensor%name, blk_size, blk_index, proc, io_unit_all, blk_values_2=blk_values_2)
            DEALLOCATE (blk_values_2)
         ENDIF
         IF (dbcsr_t_ndims(tensor) == 3) THEN
            CALL allocate_any(blk_size, blk_values_3)
            CALL dbcsr_t_get_block(tensor, blk_index, blk_size, blk_values_3)
            CALL write_block(tensor%name, blk_size, blk_index, proc, io_unit_all, blk_values_3=blk_values_3)
            DEALLOCATE (blk_values_3)
         ENDIF
         IF (dbcsr_t_ndims(tensor) == 4) THEN
            CALL allocate_any(blk_size, blk_values_4)
            CALL dbcsr_t_get_block(tensor, blk_index, blk_size, blk_values_4)
            CALL write_block(tensor%name, blk_size, blk_index, proc, io_unit_all, blk_values_4=blk_values_4)
            DEALLOCATE (blk_values_4)
         ENDIF
      ENDDO
      CALL dbcsr_t_iterator_stop(iterator)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Write a tensor block
!> \param name tensor name
!> \param blk_size block size
!> \param blk_index block index
!> \param proc which process am I
!> \param unit_nr unit number
!> \param blk_values_2 block values
!> \param blk_values_3 ...
!> \param blk_values_4 ...
! **************************************************************************************************
   SUBROUTINE write_block(name, blk_size, blk_index, proc, unit_nr, &
                          blk_values_2, blk_values_3, blk_values_4)

      CHARACTER(LEN=*), INTENT(IN)                       :: name
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size, blk_index
      INTEGER, INTENT(IN)                                :: proc, unit_nr
      REAL(KIND=real_8), &
         DIMENSION(blk_size(1), blk_size(2)), &
         INTENT(IN), OPTIONAL                            :: blk_values_2
      REAL(KIND=real_8), DIMENSION(blk_size(1), blk_size&
         (2), blk_size(3)), INTENT(IN), OPTIONAL         :: blk_values_3
      REAL(KIND=real_8), DIMENSION(blk_size(1), blk_size&
         (2), blk_size(3), blk_size(4)), INTENT(IN), &
         OPTIONAL                                        :: blk_values_4

      INTEGER                                            :: i_1, i_2, i_3, i_4, ndim

      ndim = SIZE(blk_size)

      IF (unit_nr > 0) THEN
         IF (ndim == 2) THEN
            DO i_2 = 1, blk_size(2)
            DO i_1 = 1, blk_size(1)
               WRITE (unit_nr, '(T7,A,T16,A,2I3,1X,A,1X,I3,A,1X,A,2I3,1X,A,1X,I7)') &
                  TRIM(name), "(", blk_index, ") @", proc, ':', &
                  "(", i_1, i_2, ")", &
                  INT(blk_values_2(i_1, i_2))
            ENDDO
            ENDDO
         ENDIF
         IF (ndim == 3) THEN
            DO i_3 = 1, blk_size(3)
            DO i_2 = 1, blk_size(2)
            DO i_1 = 1, blk_size(1)
               WRITE (unit_nr, '(T7,A,T16,A,3I3,1X,A,1X,I3,A,1X,A,3I3,1X,A,1X,I7)') &
                  TRIM(name), "(", blk_index, ") @", proc, ':', &
                  "(", i_1, i_2, i_3, ")", &
                  INT(blk_values_3(i_1, i_2, i_3))
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         IF (ndim == 4) THEN
            DO i_4 = 1, blk_size(4)
            DO i_3 = 1, blk_size(3)
            DO i_2 = 1, blk_size(2)
            DO i_1 = 1, blk_size(1)
               WRITE (unit_nr, '(T7,A,T16,A,4I3,1X,A,1X,I3,A,1X,A,4I3,1X,A,1X,I7)') &
                  TRIM(name), "(", blk_index, ") @", proc, ':', &
                  "(", i_1, i_2, i_3, i_4, ")", &
                  INT(blk_values_4(i_1, i_2, i_3, i_4))
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
      ENDIF
   END SUBROUTINE

END MODULE
