!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief General methods for testing DBT tensors.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_test
   #:include "dbt_macros.fypp"
   #:set maxdim = maxrank
   #:set ndims = range(2,maxdim+1)

   USE dbt_tas_base, ONLY: dbt_tas_info
   USE dbt_tas_util, ONLY: generate_larnv_seed
   USE dbt_methods, ONLY: &
      dbt_copy, dbt_get_block, dbt_iterator_type, dbt_iterator_blocks_left, &
      dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, &
      dbt_reserve_blocks, dbt_get_stored_coordinates, dbt_put_block, &
      dbt_contract, dbt_inverse_order
   USE dbt_block, ONLY: block_nd
   USE dbt_types, ONLY: &
      dbt_create, dbt_destroy, dbt_type, dbt_distribution_type, &
      dbt_distribution_destroy, &
      dims_tensor, ndims_tensor, dbt_distribution_new, dbt_nd_mp_comm, &
      mp_environ_pgrid, dbt_pgrid_type, dbt_pgrid_create, dbt_pgrid_destroy, dbt_get_info, &
      dbt_default_distvec
   USE dbt_io, ONLY: &
      dbt_write_blocks, dbt_write_block_indices
   USE kinds, ONLY: dp, default_string_length, int_8, dp
   USE message_passing, ONLY: mp_environ, &
                              mp_bcast, &
                              mp_comm_free, &
                              mp_comm_free, &
                              mp_sum, &
                              mp_cart_create
   USE dbt_allocate_wrap, ONLY: allocate_any
   USE dbt_index, ONLY: &
      combine_tensor_index, get_2d_indices_tensor, dbt_get_mapping_info
   USE dbt_tas_test, ONLY: dbt_tas_checksum

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_test'

   PUBLIC :: &
      dbt_setup_test_tensor, &
      dbt_contract_test, &
      dbt_test_formats, &
      dbt_checksum, &
      dbt_reset_randmat_seed

   INTERFACE dist_sparse_tensor_to_repl_dense_array
      #:for ndim in ndims
         MODULE PROCEDURE dist_sparse_tensor_to_repl_dense_${ndim}$d_array
      #:endfor
   END INTERFACE

   INTEGER, SAVE :: randmat_counter = 0
   INTEGER, PARAMETER, PRIVATE :: rand_seed_init = 12341313

CONTAINS

! **************************************************************************************************
!> \brief check if two (arbitrarily mapped and distributed) tensors are equal.
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_equal(tensor1, tensor2)
      TYPE(dbt_type), INTENT(INOUT)          :: tensor1, tensor2
      LOGICAL                                    :: dbt_equal

      TYPE(dbt_type)                         :: tensor2_tmp
      TYPE(dbt_iterator_type)                :: iter
      TYPE(block_nd)                             :: blk_data1, blk_data2
      INTEGER, DIMENSION(ndims_tensor(tensor1)) :: blk_size, ind_nd
      LOGICAL :: found

      ! create a copy of tensor2 that has exact same data format as tensor1
      CALL dbt_create(tensor1, tensor2_tmp)

      CALL dbt_reserve_blocks(tensor1, tensor2_tmp)
      CALL dbt_copy(tensor2, tensor2_tmp)

      dbt_equal = .TRUE.

!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor1,tensor2_tmp,dbt_equal) &
!$OMP PRIVATE(iter,ind_nd, blk_size,blk_data1,blk_data2,found)
      CALL dbt_iterator_start(iter, tensor1)

      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size)
         CALL dbt_get_block(tensor1, ind_nd, blk_data1, found)
         CPASSERT(found)
         CALL dbt_get_block(tensor2_tmp, ind_nd, blk_data2, found)
         CPASSERT(found)

         IF (.NOT. blocks_equal(blk_data1, blk_data2)) THEN
!$OMP CRITICAL
            dbt_equal = .FALSE.
!$OMP END CRITICAL
         END IF
      END DO

      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbt_destroy(tensor2_tmp)
   END FUNCTION

! **************************************************************************************************
!> \brief check if two blocks are equal
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION blocks_equal(block1, block2)
      TYPE(block_nd), INTENT(IN) :: block1, block2
      LOGICAL                    :: blocks_equal

      blocks_equal = MAXVAL(ABS(block1%blk - block2%blk)) .LT. 1.0E-12_dp

   END FUNCTION

! **************************************************************************************************
!> \brief Compute factorial
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION factorial(n)
      INTEGER, INTENT(IN) :: n
      INTEGER             :: k
      INTEGER             :: factorial
      factorial = PRODUCT((/(k, k=1, n)/))
   END FUNCTION

! **************************************************************************************************
!> \brief Compute all permutations p of (1, 2, ..., n)
!> \author Patrick Seewald
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
      RECURSIVE SUBROUTINE perm(i)
         INTEGER, INTENT(IN) :: i
         INTEGER :: j, t
         IF (i == n) THEN
            p(:, c) = pp(:)
            c = c + 1
         ELSE
            DO j = i, n
               t = pp(i)
               pp(i) = pp(j)
               pp(j) = t
               call perm(i + 1)
               t = pp(i)
               pp(i) = pp(j)
               pp(j) = t
            END DO
         END IF
      END SUBROUTINE
   END SUBROUTINE

! **************************************************************************************************
!> \brief Test equivalence of all tensor formats, using a random distribution.
!> \param blk_size_i block sizes along respective dimension
!> \param blk_ind_i index along respective dimension of non-zero blocks
!> \param ndims tensor rank
!> \param unit_nr output unit, needs to be a valid unit number on all mpi ranks
!> \param verbose if .TRUE., print all tensor blocks
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_test_formats(ndims, mp_comm, unit_nr, verbose, &
                               ${varlist("blk_size")}$, &
                               ${varlist("blk_ind")}$)
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ${varlist("blk_size")}$
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ${varlist("blk_ind")}$
      INTEGER, INTENT(IN)                         :: ndims
      INTEGER, INTENT(IN)                         :: unit_nr
      LOGICAL, INTENT(IN)                         :: verbose
      INTEGER, INTENT(IN)                         :: mp_comm
      TYPE(dbt_distribution_type)             :: dist1, dist2
      TYPE(dbt_type)                          :: tensor1, tensor2
      INTEGER                                     :: isep, iblk
      INTEGER, DIMENSION(:), ALLOCATABLE          :: ${varlist("dist1")}$, &
                                                     ${varlist("dist2")}$
      INTEGER                                     :: nblks, imap
      INTEGER, DIMENSION(ndims)                   :: pdims, myploc
      LOGICAL                                     :: eql
      INTEGER                                     :: iperm, idist, icount
      INTEGER, DIMENSION(:), ALLOCATABLE          :: map1, map2, map1_ref, map2_ref
      INTEGER, DIMENSION(ndims, factorial(ndims)) :: perm
      INTEGER                                     :: io_unit
      INTEGER                                     :: mynode, numnodes
      TYPE(dbt_pgrid_type)                    :: comm_nd
      CHARACTER(LEN=default_string_length)        :: tensor_name

      ! Process grid
      pdims(:) = 0
      CALL dbt_pgrid_create(mp_comm, pdims, comm_nd)
      CALL mp_environ(numnodes, mynode, mp_comm)

      io_unit = 0
      IF (mynode .EQ. 0) io_unit = unit_nr

      CALL permute(ndims, perm)
      ALLOCATE (map1_ref, source=perm(1:ndims/2, 1))
      ALLOCATE (map2_ref, source=perm(ndims/2 + 1:ndims, 1))

      IF (io_unit > 0) THEN
         WRITE (io_unit, *)
         WRITE (io_unit, '(A)') repeat("-", 80)
         WRITE (io_unit, '(A,1X,I1)') "Testing matrix representations of tensor rank", ndims
         WRITE (io_unit, '(A)') repeat("-", 80)
         WRITE (io_unit, '(A)') "Block sizes:"

         #:for dim in range(1, maxdim+1)
            IF (ndims >= ${dim}$) THEN
               WRITE (io_unit, '(T4,A,1X,I1,A,1X)', advance='no') 'Dim', ${dim}$, ':'
               DO iblk = 1, SIZE(blk_size_${dim}$)
                  WRITE (io_unit, '(I2,1X)', advance='no') blk_size_${dim}$ (iblk)
               END DO
               WRITE (io_unit, *)
            END IF
         #:endfor

         WRITE (io_unit, '(A)') "Non-zero blocks:"
         DO iblk = 1, SIZE(blk_ind_1)
            #:for ndim in ndims
               IF (ndims == ${ndim}$) THEN
                  WRITE (io_unit, '(T4,A, I3, A, ${ndim}$I3, 1X, A)') &
                     'Block', iblk, ': (', ${varlist("blk_ind", nmax=ndim, suffix='(iblk)')}$, ')'
               END IF
            #:endfor
         END DO

         WRITE (io_unit, *)
         WRITE (io_unit, '(A,1X)', advance='no') "Reference map:"
         WRITE (io_unit, '(A1,1X)', advance='no') "("
         DO imap = 1, SIZE(map1_ref)
            WRITE (io_unit, '(I1,1X)', advance='no') map1_ref(imap)
         END DO
         WRITE (io_unit, '(A1,1X)', advance='no') "|"
         DO imap = 1, SIZE(map2_ref)
            WRITE (io_unit, '(I1,1X)', advance='no') map2_ref(imap)
         END DO
         WRITE (io_unit, '(A1)') ")"

      END IF

      icount = 0
      DO iperm = 1, factorial(ndims)
         DO isep = 1, ndims - 1
            icount = icount + 1

            ALLOCATE (map1, source=perm(1:isep, iperm))
            ALLOCATE (map2, source=perm(isep + 1:ndims, iperm))

            CALL mp_environ(numnodes, mynode, mp_comm)
            CALL mp_environ_pgrid(comm_nd, pdims, myploc)

            #:for dim in range(1, maxdim+1)
               IF (${dim}$ <= ndims) THEN
                  nblks = SIZE(blk_size_${dim}$)
                  ALLOCATE (dist1_${dim}$ (nblks))
                  ALLOCATE (dist2_${dim}$ (nblks))
                  CALL dbt_default_distvec(nblks, pdims(${dim}$), blk_size_${dim}$, dist1_${dim}$)
                  CALL dbt_default_distvec(nblks, pdims(${dim}$), blk_size_${dim}$, dist2_${dim}$)
               END IF
            #:endfor

            WRITE (tensor_name, '(A,1X,I3,1X)') "Test", icount

            IF (io_unit > 0) THEN
               WRITE (io_unit, *)
               WRITE (io_unit, '(A,A,1X)', advance='no') TRIM(tensor_name), ':'
               WRITE (io_unit, '(A1,1X)', advance='no') "("
               DO imap = 1, SIZE(map1)
                  WRITE (io_unit, '(I1,1X)', advance='no') map1(imap)
               END DO
               WRITE (io_unit, '(A1,1X)', advance='no') "|"
               DO imap = 1, SIZE(map2)
                  WRITE (io_unit, '(I1,1X)', advance='no') map2(imap)
               END DO
               WRITE (io_unit, '(A1)') ")"

               WRITE (io_unit, '(T4,A)') "Reference distribution:"
               #:for dim in range(1, maxdim+1)
                  IF (${dim}$ <= ndims) THEN
                     WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec ${dim}$:"
                     DO idist = 1, SIZE(dist2_${dim}$)
                        WRITE (io_unit, '(I2,1X)', advance='no') dist2_${dim}$ (idist)
                     END DO
                     WRITE (io_unit, *)
                  END IF
               #:endfor

               WRITE (io_unit, '(T4,A)') "Test distribution:"
               #:for dim in range(1, maxdim+1)
                  IF (${dim}$ <= ndims) THEN
                     WRITE (io_unit, '(T7,A,1X)', advance='no') "Dist vec ${dim}$:"
                     DO idist = 1, SIZE(dist2_${dim}$)
                        WRITE (io_unit, '(I2,1X)', advance='no') dist1_${dim}$ (idist)
                     END DO
                     WRITE (io_unit, *)
                  END IF
               #:endfor
            END IF

            #:for ndim in ndims
               IF (ndims == ${ndim}$) THEN
                  CALL dbt_distribution_new(dist2, comm_nd, ${varlist("dist2", nmax=ndim)}$)
                  CALL dbt_create(tensor2, "Ref", dist2, map1_ref, map2_ref, &
                                  ${varlist("blk_size", nmax=ndim)}$)
                  CALL dbt_setup_test_tensor(tensor2, comm_nd%mp_comm_2d, .TRUE., ${varlist("blk_ind", nmax=ndim)}$)
               END IF
            #:endfor

            IF (verbose) CALL dbt_write_blocks(tensor2, io_unit, unit_nr)

            #:for ndim in ndims
               IF (ndims == ${ndim}$) THEN
                  CALL dbt_distribution_new(dist1, comm_nd, ${varlist("dist1", nmax=ndim)}$)
                  CALL dbt_create(tensor1, tensor_name, dist1, map1, map2, &
                                  ${varlist("blk_size", nmax=ndim)}$)
                  CALL dbt_setup_test_tensor(tensor1, comm_nd%mp_comm_2d, .TRUE., ${varlist("blk_ind", nmax=ndim)}$)
               END IF
            #:endfor

            IF (verbose) CALL dbt_write_blocks(tensor1, io_unit, unit_nr)

            eql = dbt_equal(tensor1, tensor2)

            IF (.NOT. eql) THEN
               IF (io_unit > 0) WRITE (io_unit, '(A,1X,A)') TRIM(tensor_name), 'Test failed!'
               CPABORT('')
            ELSE
               IF (io_unit > 0) WRITE (io_unit, '(A,1X,A)') TRIM(tensor_name), 'Test passed!'
            END IF
            DEALLOCATE (map1, map2)

            CALL dbt_destroy(tensor1)
            CALL dbt_distribution_destroy(dist1)

            CALL dbt_destroy(tensor2)
            CALL dbt_distribution_destroy(dist2)

            #:for dim in range(1, maxdim+1)
               IF (${dim}$ <= ndims) THEN
                  DEALLOCATE (dist1_${dim}$, dist2_${dim}$)
               END IF
            #:endfor

         END DO
      END DO
      CALL dbt_pgrid_destroy(comm_nd)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Allocate and fill test tensor - entries are enumerated by their index s.t. they only depend
!>        on global properties of the tensor but not on distribution, matrix representation, etc.
!> \param mp_comm communicator
!> \param blk_ind_i index along respective dimension of non-zero blocks
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_setup_test_tensor(tensor, mp_comm, enumerate, ${varlist("blk_ind")}$)
      TYPE(dbt_type), INTENT(INOUT)                  :: tensor
      INTEGER, INTENT(IN)                                :: mp_comm
      LOGICAL, INTENT(IN)                                :: enumerate
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: ${varlist("blk_ind")}$
      INTEGER                                            :: numnodes, mynode

      INTEGER                                            :: i, ib, my_nblks_alloc, nblks_alloc, proc, nze
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ${varlist("my_blk_ind")}$
      INTEGER, DIMENSION(ndims_tensor(tensor))          :: blk_index, blk_offset, blk_size, &
                                                           tensor_dims
      INTEGER, DIMENSION(:, :), ALLOCATABLE               :: ind_nd
      #:for ndim in ndims
         REAL(KIND=dp), ALLOCATABLE, &
            DIMENSION(${shape_colon(ndim)}$)                :: blk_values_${ndim}$
      #:endfor
      TYPE(dbt_iterator_type)                        :: iterator
      INTEGER, DIMENSION(4)                              :: iseed
      INTEGER, DIMENSION(2)                              :: blk_index_2d, nblks_2d

      nblks_alloc = SIZE(blk_ind_1)
      CALL mp_environ(numnodes, mynode, mp_comm)

      IF (.NOT. enumerate) THEN
         CPASSERT(randmat_counter .NE. 0)

         randmat_counter = randmat_counter + 1
      END IF

      ALLOCATE (ind_nd(nblks_alloc, ndims_tensor(tensor)))

!$OMP PARALLEL DEFAULT(NONE) SHARED(ind_nd,nblks_alloc,tensor,mynode,enumerate,randmat_counter) &
!$OMP SHARED(${varlist("blk_ind")}$) &
!$OMP PRIVATE(my_nblks_alloc,ib,proc,i,iterator,blk_offset,blk_size,blk_index) &
!$OMP PRIVATE(blk_index_2d,nblks_2d,iseed,nze,tensor_dims) &
!$OMP PRIVATE(${varlist("blk_values", nmin=2)}$) &
!$OMP PRIVATE(${varlist("my_blk_ind")}$)

      my_nblks_alloc = 0
!$OMP DO
      DO ib = 1, nblks_alloc
         #:for ndim in ndims
            IF (ndims_tensor(tensor) == ${ndim}$) THEN
               ind_nd(ib, :) = [${varlist("blk_ind", nmax=ndim, suffix="(ib)")}$]
            END IF
         #:endfor
         CALL dbt_get_stored_coordinates(tensor, ind_nd(ib, :), proc)
         IF (proc == mynode) THEN
            my_nblks_alloc = my_nblks_alloc + 1
         END IF
      END DO
!$OMP END DO

      #:for dim in range(1, maxdim+1)
         IF (ndims_tensor(tensor) >= ${dim}$) THEN
            ALLOCATE (my_blk_ind_${dim}$ (my_nblks_alloc))
         END IF
      #:endfor

      i = 0
!$OMP DO
      DO ib = 1, nblks_alloc
         CALL dbt_get_stored_coordinates(tensor, ind_nd(ib, :), proc)
         IF (proc == mynode) THEN
            i = i + 1
            #:for dim in range(1, maxdim+1)
               IF (ndims_tensor(tensor) >= ${dim}$) THEN
                  my_blk_ind_${dim}$ (i) = blk_ind_${dim}$ (ib)
               END IF
            #:endfor
         END IF
      END DO
!$OMP END DO

      #:for ndim in ndims
         IF (ndims_tensor(tensor) == ${ndim}$) THEN
            CALL dbt_reserve_blocks(tensor, ${varlist("my_blk_ind", nmax=ndim)}$)
         END IF
      #:endfor

      CALL dbt_iterator_start(iterator, tensor)
      DO WHILE (dbt_iterator_blocks_left(iterator))
         CALL dbt_iterator_next_block(iterator, blk_index, blk_size=blk_size, blk_offset=blk_offset)

         IF (.NOT. enumerate) THEN
            blk_index_2d = INT(get_2d_indices_tensor(tensor%nd_index_blk, blk_index))
            CALL dbt_get_mapping_info(tensor%nd_index_blk, dims_2d=nblks_2d)
            iseed = generate_larnv_seed(blk_index_2d(1), nblks_2d(1), blk_index_2d(2), nblks_2d(2), randmat_counter)
            nze = PRODUCT(blk_size)
         END IF

         #:for ndim in ndims
            IF (ndims_tensor(tensor) == ${ndim}$) THEN
               CALL allocate_any(blk_values_${ndim}$, shape_spec=blk_size)
               CALL dims_tensor(tensor, tensor_dims)
               IF (enumerate) THEN
                  CALL enumerate_block_elements(blk_size, blk_offset, tensor_dims, blk_${ndim}$=blk_values_${ndim}$)
               ELSE
                  CALL dlarnv(1, iseed, nze, blk_values_${ndim}$)
               END IF
               CALL dbt_put_block(tensor, blk_index, blk_size, blk_values_${ndim}$)
               DEALLOCATE (blk_values_${ndim}$)
            END IF
         #:endfor
      END DO
      CALL dbt_iterator_stop(iterator)
!$OMP END PARALLEL

   END SUBROUTINE

! **************************************************************************************************
!> \brief Enumerate tensor entries in block
!> \param blk_2 block values for 2 dimensions
!> \param blk_3 block values for 3 dimensions
!> \param blk_size size of block
!> \param blk_offset block offset (indices of first element)
!> \param tensor_size global tensor sizes
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE enumerate_block_elements(blk_size, blk_offset, tensor_size, ${varlist("blk", nmin=2)}$)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size, blk_offset, tensor_size
      #:for ndim in ndims
         REAL(KIND=dp), DIMENSION(${shape_colon(ndim)}$), &
            OPTIONAL, INTENT(OUT)                           :: blk_${ndim}$
      #:endfor
      INTEGER                                            :: ndim
      INTEGER, DIMENSION(SIZE(blk_size))                 :: arr_ind, tens_ind
      INTEGER                                            :: ${varlist("i")}$

      ndim = SIZE(tensor_size)

      #:for ndim in ndims
         IF (ndim == ${ndim}$) THEN
            #:for idim in range(ndim,0,-1)
               DO i_${idim}$ = 1, blk_size(${idim}$)
                  #:endfor
                  arr_ind(:) = [${varlist("i", nmax=ndim)}$]
                  tens_ind(:) = arr_ind(:) + blk_offset(:) - 1
                  blk_${ndim}$ (${arrlist("arr_ind", nmax=ndim)}$) = combine_tensor_index(tens_ind, tensor_size)
                  #:for idim in range(ndim,0,-1)
                     END DO
                  #:endfor
               END IF
            #:endfor

         END SUBROUTINE

         #:for ndim in ndims
! **************************************************************************************************
!> \brief Transform a distributed sparse tensor to a replicated dense array. This is only useful for
!>        testing tensor contraction by matrix multiplication of dense arrays.
!> \author Patrick Seewald
! **************************************************************************************************
            SUBROUTINE dist_sparse_tensor_to_repl_dense_${ndim}$d_array(tensor, array)

               TYPE(dbt_type), INTENT(INOUT)                          :: tensor
               REAL(dp), ALLOCATABLE, DIMENSION(${shape_colon(ndim)}$), &
                  INTENT(OUT)                                             :: array
               REAL(dp), ALLOCATABLE, DIMENSION(${shape_colon(ndim)}$)   :: block
               INTEGER, DIMENSION(ndims_tensor(tensor))                  :: dims_nd, ind_nd, blk_size, blk_offset
               TYPE(dbt_iterator_type)                                     :: iterator
               INTEGER                                                    :: idim
               INTEGER, DIMENSION(ndims_tensor(tensor))                   :: blk_start, blk_end
               LOGICAL                                                    :: found

               CPASSERT(ndims_tensor(tensor) .EQ. ${ndim}$)
               CALL dbt_get_info(tensor, nfull_total=dims_nd)
               CALL allocate_any(array, shape_spec=dims_nd)
               array(${shape_colon(ndim)}$) = 0.0_dp

!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor,array) &
!$OMP PRIVATE(iterator,ind_nd,blk_size,blk_offset,block,found,idim,blk_start,blk_end)
               CALL dbt_iterator_start(iterator, tensor)
               DO WHILE (dbt_iterator_blocks_left(iterator))
                  CALL dbt_iterator_next_block(iterator, ind_nd, blk_size=blk_size, blk_offset=blk_offset)
                  CALL dbt_get_block(tensor, ind_nd, block, found)
                  CPASSERT(found)

                  DO idim = 1, ndims_tensor(tensor)
                     blk_start(idim) = blk_offset(idim)
                     blk_end(idim) = blk_offset(idim) + blk_size(idim) - 1
                  END DO
                  array(${", ".join(["blk_start("+str(idim)+"):blk_end("+str(idim)+")" for idim in range(1, ndim + 1)])}$) = &
                     block(${shape_colon(ndim)}$)

                  DEALLOCATE (block)
               END DO
               CALL dbt_iterator_stop(iterator)
!$OMP END PARALLEL
               CALL mp_sum(array, tensor%pgrid%mp_comm_2d)

            END SUBROUTINE
         #:endfor

! **************************************************************************************************
!> \brief test tensor contraction
!> \note for testing/debugging, simply replace a call to dbt_contract with a call to this routine
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_contract_test(alpha, tensor_1, tensor_2, beta, tensor_3, &
                                      contract_1, notcontract_1, &
                                      contract_2, notcontract_2, &
                                      map_1, map_2, &
                                      unit_nr, &
                                      bounds_1, bounds_2, bounds_3, &
                                      log_verbose, write_int)

            REAL(dp), INTENT(IN) :: alpha
            TYPE(dbt_type), INTENT(INOUT)    :: tensor_1, tensor_2, tensor_3
            REAL(dp), INTENT(IN) :: beta
            INTEGER, DIMENSION(:), INTENT(IN)    :: contract_1, contract_2, &
                                                    notcontract_1, notcontract_2, &
                                                    map_1, map_2
            INTEGER, INTENT(IN)                  :: unit_nr
            INTEGER, DIMENSION(2, SIZE(contract_1)), &
               OPTIONAL                          :: bounds_1
            INTEGER, DIMENSION(2, SIZE(notcontract_1)), &
               OPTIONAL                          :: bounds_2
            INTEGER, DIMENSION(2, SIZE(notcontract_2)), &
               OPTIONAL                          :: bounds_3
            LOGICAL, INTENT(IN), OPTIONAL        :: log_verbose
            LOGICAL, INTENT(IN), OPTIONAL        :: write_int
            INTEGER                              :: io_unit, mynode, numnodes, mp_comm
            INTEGER, DIMENSION(:), ALLOCATABLE   :: size_1, size_2, size_3, &
                                                    order_t1, order_t2, order_t3
            INTEGER, DIMENSION(2, ndims_tensor(tensor_1)) :: bounds_t1
            INTEGER, DIMENSION(2, ndims_tensor(tensor_2)) :: bounds_t2

            #:for ndim in ndims
               REAL(KIND=dp), ALLOCATABLE, &
                  DIMENSION(${shape_colon(ndim)}$) :: array_1_${ndim}$d, &
                                                      array_2_${ndim}$d, &
                                                      array_3_${ndim}$d, &
                                                      array_1_${ndim}$d_full, &
                                                      array_2_${ndim}$d_full, &
                                                      array_3_0_${ndim}$d, &
                                                      array_1_rs${ndim}$d, &
                                                      array_2_rs${ndim}$d, &
                                                      array_3_rs${ndim}$d, &
                                                      array_3_0_rs${ndim}$d
            #:endfor
            REAL(KIND=dp), ALLOCATABLE, &
               DIMENSION(:, :)                   :: array_1_mm, &
                                                    array_2_mm, &
                                                    array_3_mm, &
                                                    array_3_test_mm
            LOGICAL                             :: eql, notzero
            LOGICAL, PARAMETER                  :: debug = .FALSE.
            REAL(KIND=dp)                   :: cs_1, cs_2, cs_3, eql_diff
            LOGICAL                             :: do_crop_1, do_crop_2

            mp_comm = tensor_1%pgrid%mp_comm_2d
            CALL mp_environ(numnodes, mynode, mp_comm)
            io_unit = -1
            IF (mynode .EQ. 0) io_unit = unit_nr

            cs_1 = dbt_checksum(tensor_1)
            cs_2 = dbt_checksum(tensor_2)
            cs_3 = dbt_checksum(tensor_3)

            IF (io_unit > 0) THEN
               WRITE (io_unit, *)
               WRITE (io_unit, '(A)') repeat("-", 80)
               WRITE (io_unit, '(A,1X,A,1X,A,1X,A,1X,A,1X,A)') "Testing tensor contraction", &
                  TRIM(tensor_1%name), "x", TRIM(tensor_2%name), "=", TRIM(tensor_3%name)
               WRITE (io_unit, '(A)') repeat("-", 80)
            END IF

            IF (debug) THEN
               IF (io_unit > 0) THEN
                  WRITE (io_unit, "(A, E9.2)") "checksum ", TRIM(tensor_1%name), cs_1
                  WRITE (io_unit, "(A, E9.2)") "checksum ", TRIM(tensor_2%name), cs_2
                  WRITE (io_unit, "(A, E9.2)") "checksum ", TRIM(tensor_3%name), cs_3
               END IF
            END IF

            IF (debug) THEN
               CALL dbt_write_block_indices(tensor_1, io_unit, unit_nr)
               CALL dbt_write_blocks(tensor_1, io_unit, unit_nr, write_int)
            END IF

            SELECT CASE (ndims_tensor(tensor_3))
               #:for ndim in ndims
                  CASE (${ndim}$)
                  CALL dist_sparse_tensor_to_repl_dense_array(tensor_3, array_3_0_${ndim}$d)
               #:endfor
            END SELECT

            CALL dbt_contract(alpha, tensor_1, tensor_2, beta, tensor_3, &
                              contract_1, notcontract_1, &
                              contract_2, notcontract_2, &
                              map_1, map_2, &
                              bounds_1=bounds_1, bounds_2=bounds_2, bounds_3=bounds_3, &
                              filter_eps=1.0E-12_dp, &
                              unit_nr=io_unit, log_verbose=log_verbose)

            cs_3 = dbt_checksum(tensor_3)

            IF (debug) THEN
               IF (io_unit > 0) THEN
                  WRITE (io_unit, "(A, E9.2)") "checksum ", TRIM(tensor_3%name), cs_3
               END IF
            END IF

            do_crop_1 = .FALSE.; do_crop_2 = .FALSE.!; do_crop_3 = .FALSE.

            ! crop tensor as first step
            bounds_t1(1, :) = 1
            CALL dbt_get_info(tensor_1, nfull_total=bounds_t1(2, :))

            bounds_t2(1, :) = 1
            CALL dbt_get_info(tensor_2, nfull_total=bounds_t2(2, :))

            IF (PRESENT(bounds_1)) THEN
               bounds_t1(:, contract_1) = bounds_1
               do_crop_1 = .TRUE.
               bounds_t2(:, contract_2) = bounds_1
               do_crop_2 = .TRUE.
            END IF

            IF (PRESENT(bounds_2)) THEN
               bounds_t1(:, notcontract_1) = bounds_2
               do_crop_1 = .TRUE.
            END IF

            IF (PRESENT(bounds_3)) THEN
               bounds_t2(:, notcontract_2) = bounds_3
               do_crop_2 = .TRUE.
            END IF

            ! Convert tensors to simple multidimensional arrays
            #:for i in range(1,4)
               SELECT CASE (ndims_tensor(tensor_${i}$))
                  #:for ndim in ndims
                     CASE (${ndim}$)
                     #:if i < 3
                        CALL dist_sparse_tensor_to_repl_dense_array(tensor_${i}$, array_${i}$_${ndim}$d_full)
                        CALL allocate_any(array_${i}$_${ndim}$d, shape_spec=SHAPE(array_${i}$_${ndim}$d_full))
                        array_${i}$_${ndim}$d = 0.0_dp
         array_${i}$_${ndim}$d(${", ".join(["bounds_t" + str(i) + "(1, " + str(idim) + "):bounds_t" + str(i) + "(2, " + str(idim) + ")" for idim in range(1, ndim+1)])}$) = &
         array_${i}$_${ndim}$d_full(${", ".join(["bounds_t" + str(i) + "(1, " + str(idim) + "):bounds_t" + str(i) + "(2, " + str(idim) + ")" for idim in range(1, ndim+1)])}$)
                     #:else
                        CALL dist_sparse_tensor_to_repl_dense_array(tensor_${i}$, array_${i}$_${ndim}$d)
                     #:endif

                  #:endfor
               END SELECT
            #:endfor

            ! Get array sizes

            #:for i in range(1,4)
               SELECT CASE (ndims_tensor(tensor_${i}$))
                  #:for ndim in ndims
                     CASE (${ndim}$)
                     ALLOCATE (size_${i}$, source=SHAPE(array_${i}$_${ndim}$d))

                  #:endfor
               END SELECT
            #:endfor

            #:for i in range(1,4)
               ALLOCATE (order_t${i}$ (ndims_tensor(tensor_${i}$)))
            #:endfor

            ASSOCIATE (map_t1_1 => notcontract_1, map_t1_2 => contract_1, &
                       map_t2_1 => notcontract_2, map_t2_2 => contract_2, &
                       map_t3_1 => map_1, map_t3_2 => map_2)

               #:for i in range(1,4)
                  order_t${i}$ (:) = dbt_inverse_order([map_t${i}$_1, map_t${i}$_2])

                  SELECT CASE (ndims_tensor(tensor_${i}$))
                     #:for ndim in ndims
                        CASE (${ndim}$)
                        CALL allocate_any(array_${i}$_rs${ndim}$d, source=array_${i}$_${ndim}$d, order=order_t${i}$)
                        CALL allocate_any(array_${i}$_mm, sizes_2d(size_${i}$, map_t${i}$_1, map_t${i}$_2))
                        array_${i}$_mm(:, :) = RESHAPE(array_${i}$_rs${ndim}$d, SHAPE(array_${i}$_mm))
                     #:endfor
                  END SELECT
               #:endfor

               SELECT CASE (ndims_tensor(tensor_3))
                  #:for ndim in ndims
                     CASE (${ndim}$)
                     CALL allocate_any(array_3_0_rs${ndim}$d, source=array_3_0_${ndim}$d, order=order_t3)
                     CALL allocate_any(array_3_test_mm, sizes_2d(size_3, map_t3_1, map_t3_2))
                     array_3_test_mm(:, :) = RESHAPE(array_3_0_rs${ndim}$d, SHAPE(array_3_mm))
                  #:endfor
               END SELECT

               array_3_test_mm(:, :) = beta*array_3_test_mm(:, :) + alpha*MATMUL(array_1_mm, transpose(array_2_mm))

            END ASSOCIATE

            eql_diff = MAXVAL(ABS(array_3_test_mm(:, :) - array_3_mm(:, :)))
            notzero = MAXVAL(ABS(array_3_test_mm(:, :))) .GT. 1.0E-12_dp

            eql = eql_diff .LT. 1.0E-11_dp

            IF (.NOT. eql .OR. .NOT. notzero) THEN
               IF (io_unit > 0) WRITE (io_unit, *) 'Test failed!', eql_diff
               CPABORT('')
            ELSE
               IF (io_unit > 0) WRITE (io_unit, *) 'Test passed!', eql_diff
            END IF

         END SUBROUTINE

! **************************************************************************************************
!> \brief mapped sizes in 2d
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION sizes_2d(nd_sizes, map1, map2)
            INTEGER, DIMENSION(:), INTENT(IN) :: nd_sizes, map1, map2
            INTEGER, DIMENSION(2)             :: sizes_2d
            sizes_2d(1) = PRODUCT(nd_sizes(map1))
            sizes_2d(2) = PRODUCT(nd_sizes(map2))
         END FUNCTION

! **************************************************************************************************
!> \brief checksum of a tensor consistent with block_checksum
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION dbt_checksum(tensor)
            TYPE(dbt_type), INTENT(IN) :: tensor
            REAL(KIND=dp) :: dbt_checksum
            dbt_checksum = dbt_tas_checksum(tensor%matrix_rep)
         END FUNCTION

! **************************************************************************************************
!> \brief Reset the seed used for generating random matrices to default value
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_reset_randmat_seed()
            randmat_counter = rand_seed_init
         END SUBROUTINE

      END MODULE
