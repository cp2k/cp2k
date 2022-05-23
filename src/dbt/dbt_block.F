!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods to operate on n-dimensional tensor blocks.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_block

   #:include "dbt_macros.fypp"
   #:set maxdim = maxrank
   #:set ndims = range(2,maxdim+1)

   USE OMP_LIB, ONLY: omp_get_thread_num, omp_get_num_threads
   USE dbcsr_api, ONLY: &
      dbcsr_type, dbcsr_release, &
      dbcsr_iterator_type, dbcsr_iterator_start, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
      dbcsr_has_symmetry, dbcsr_desymmetrize, dbcsr_get_num_blocks, dbcsr_iterator_stop, &
      dbcsr_reserve_blocks, dbcsr_finalize
   USE dbt_allocate_wrap, ONLY: &
      allocate_any
   USE dbt_tas_types, ONLY: &
      dbt_tas_iterator
   USE dbt_tas_base, ONLY: &
      dbt_tas_iterator_next_block, dbt_tas_iterator_blocks_left, dbt_tas_iterator_start, &
      dbt_tas_iterator_stop, dbt_tas_get_block_p, dbt_tas_put_block, dbt_tas_reserve_blocks, &
      dbt_tas_iterator_num_blocks
   USE kinds, ONLY: dp, int_8, dp
   USE dbt_index, ONLY: &
      nd_to_2d_mapping, ndims_mapping, get_nd_indices_tensor, destroy_nd_to_2d_mapping, get_2d_indices_tensor, &
      create_nd_to_2d_mapping
   USE dbt_array_list_methods, ONLY: &
      array_list, get_array_elements, destroy_array_list, sizes_of_arrays, create_array_list, &
      get_arrays
   USE dbt_types, ONLY: &
      dbt_type, ndims_tensor, dbt_blk_sizes, dbt_get_num_blocks, &
      dbt_finalize, ndims_matrix_row, ndims_matrix_column

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_block'

   PUBLIC :: &
      block_nd, &
      create_block, &
      dbt_get_block, &
      dbt_iterator_num_blocks, &
      dbt_iterator_blocks_left, &
      dbt_iterator_next_block, &
      dbt_iterator_start, &
      dbt_iterator_stop, &
      dbt_iterator_type, &
      dbt_put_block, &
      dbt_reserve_blocks, &
      destroy_block, &
      checker_tr, &
      ndims_iterator

   TYPE dbt_iterator_type
      TYPE(dbt_tas_iterator)      :: iter
      TYPE(dbt_type), POINTER     :: tensor => NULL()
   END TYPE dbt_iterator_type

   TYPE block_nd
      INTEGER, DIMENSION(:), ALLOCATABLE   :: sizes
      REAL(dp), DIMENSION(:), ALLOCATABLE :: blk
   END TYPE

   INTERFACE create_block
      MODULE PROCEDURE create_block_data
      MODULE PROCEDURE create_block_nodata
   END INTERFACE

   INTERFACE dbt_put_block
      #:for ndim in ndims
         MODULE PROCEDURE dbt_put_${ndim}$d_block
      #:endfor
      MODULE PROCEDURE dbt_put_anyd_block
   END INTERFACE

   INTERFACE dbt_get_block
      #:for ndim in ndims
         MODULE PROCEDURE dbt_get_${ndim}$d_block
         MODULE PROCEDURE dbt_allocate_and_get_${ndim}$d_block
      #:endfor
      MODULE PROCEDURE dbt_get_anyd_block
   END INTERFACE

   INTERFACE dbt_reserve_blocks
      MODULE PROCEDURE dbt_reserve_blocks_index
      MODULE PROCEDURE dbt_reserve_blocks_index_array
      MODULE PROCEDURE dbt_reserve_blocks_template
      MODULE PROCEDURE dbt_reserve_blocks_tensor_to_matrix
      MODULE PROCEDURE dbt_reserve_blocks_matrix_to_tensor
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief block size
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION block_size(block)
      TYPE(block_nd), INTENT(IN)         :: block
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_size

      ALLOCATE (block_size, source=block%sizes)
   END FUNCTION

! **************************************************************************************************
!> \brief Generalization of block_iterator_start for tensors.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_iterator_start(iterator, tensor)
      TYPE(dbt_iterator_type), INTENT(OUT)           :: iterator
      TYPE(dbt_type), INTENT(IN), TARGET             :: tensor

      CPASSERT(tensor%valid)
      CALL dbt_tas_iterator_start(iterator%iter, tensor%matrix_rep)
      iterator%tensor => tensor
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of block_iterator_stop for tensors.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_iterator_stop(iterator)
      TYPE(dbt_iterator_type), INTENT(INOUT) :: iterator

      CALL dbt_tas_iterator_stop(iterator%iter)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Number of dimensions.
!> \note specification function below must be defined before it is used in
!>       the source due to a bug in the IBM XL Fortran compiler (compilation fails)
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION ndims_iterator(iterator)
      TYPE(dbt_iterator_type), INTENT(IN) :: iterator
      INTEGER                                 :: ndims_iterator

      ndims_iterator = iterator%tensor%nd_index%ndim_nd
   END FUNCTION

! **************************************************************************************************
!> \brief iterate over nd blocks of an nd rank tensor, index only (blocks must be retrieved by
!>        calling dbt_get_block on tensor).
!> \param ind_nd nd index of block
!> \param blk_size blk size in each dimension
!> \param blk_offset blk offset in each dimension
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_iterator_next_block(iterator, ind_nd, blk_size, blk_offset)
      !!
      TYPE(dbt_iterator_type), INTENT(INOUT)     :: iterator
      INTEGER, DIMENSION(ndims_iterator(iterator)), &
         INTENT(OUT)                                 :: ind_nd
      INTEGER, DIMENSION(ndims_iterator(iterator)), &
         INTENT(OUT), OPTIONAL                       :: blk_size, blk_offset

      INTEGER(KIND=int_8), DIMENSION(2)              :: ind_2d

      CALL dbt_tas_iterator_next_block(iterator%iter, ind_2d(1), ind_2d(2))

      ind_nd(:) = get_nd_indices_tensor(iterator%tensor%nd_index_blk, ind_2d)
      IF (PRESENT(blk_size)) blk_size(:) = get_array_elements(iterator%tensor%blk_sizes, ind_nd)
      ! note: blk_offset needs to be determined by tensor metadata, can not be derived from 2d row/col
      ! offset since block index mapping is not consistent with element index mapping
      IF (PRESENT(blk_offset)) blk_offset(:) = get_array_elements(iterator%tensor%blk_offsets, ind_nd)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of block_iterator_num_blocks for tensors.
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION dbt_iterator_num_blocks(iterator)
      TYPE(dbt_iterator_type), INTENT(IN) :: iterator
      INTEGER                             :: dbt_iterator_num_blocks

      dbt_iterator_num_blocks = dbt_tas_iterator_num_blocks(iterator%iter)

   END FUNCTION

! **************************************************************************************************
!> \brief Generalization of block_iterator_blocks_left for tensors.
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION dbt_iterator_blocks_left(iterator)
      TYPE(dbt_iterator_type), INTENT(IN) :: iterator
      LOGICAL                                 :: dbt_iterator_blocks_left

      dbt_iterator_blocks_left = dbt_tas_iterator_blocks_left(iterator%iter)

   END FUNCTION

! **************************************************************************************************
!> \brief reserve blocks from indices as array object
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_reserve_blocks_index_array(tensor, blk_ind)
      TYPE(dbt_type), INTENT(INOUT)   :: tensor
      INTEGER, DIMENSION(:, :), INTENT(IN) :: blk_ind
      INTEGER                             :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_reserve_blocks_index_array'

      CALL timeset(routineN, handle)
      #:for ndim in ndims
         IF (ndims_tensor(tensor) == ${ndim}$) THEN
            CALL dbt_reserve_blocks(tensor, ${arrlist("blk_ind", nmax=ndim, ndim_pre=1)}$)
         END IF
      #:endfor
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief reserve tensor blocks using block indices
!> \param blk_ind index of blocks to reserve in each dimension
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_reserve_blocks_index(tensor, ${varlist("blk_ind")}$)
      TYPE(dbt_type), INTENT(INOUT)           :: tensor
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ${varlist("blk_ind")}$
      INTEGER                                     :: iblk, nblk, handle
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)          :: cols, rows
      INTEGER(KIND=int_8), DIMENSION(2)                       :: ind_2d
      TYPE(array_list)                            :: blks
      INTEGER, DIMENSION(ndims_tensor(tensor))   :: iblk_nd, ind_nd, nblk_tmp
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_reserve_blocks_index'

      CALL timeset(routineN, handle)
      CPASSERT(tensor%valid)

      CALL create_array_list(blks, ndims_tensor(tensor), &
                             ${varlist("blk_ind")}$)

      nblk_tmp(:) = sizes_of_arrays(blks)
      nblk = nblk_tmp(1)
      ALLOCATE (cols(nblk), rows(nblk))
      DO iblk = 1, nblk
         iblk_nd(:) = iblk
         ind_nd(:) = get_array_elements(blks, iblk_nd)
         ind_2d(:) = get_2d_indices_tensor(tensor%nd_index_blk, ind_nd)
         rows(iblk) = ind_2d(1); cols(iblk) = ind_2d(2)
      END DO

      CALL dbt_tas_reserve_blocks(tensor%matrix_rep, rows=rows, columns=cols)
      CALL dbt_finalize(tensor)
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief reserve tensor blocks using template
!> \param tensor_in template tensor
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_reserve_blocks_template(tensor_in, tensor_out)
      TYPE(dbt_type), INTENT(IN)           :: tensor_in
      TYPE(dbt_type), INTENT(INOUT)        :: tensor_out

      CHARACTER(LEN=*), PARAMETER          :: routineN = 'dbt_reserve_blocks_template'

      TYPE(dbt_iterator_type)              :: iter
      INTEGER                              :: handle, nblk, iblk
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: blk_ind

      CALL timeset(routineN, handle)

!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor_in,tensor_out) &
!$OMP PRIVATE(iter,nblk,iblk,blk_ind)
      CALL dbt_iterator_start(iter, tensor_in)
      nblk = dbt_iterator_num_blocks(iter)
      ALLOCATE (blk_ind(nblk, ndims_tensor(tensor_in)))
      DO iblk = 1, nblk
         CALL dbt_iterator_next_block(iter, ind_nd=blk_ind(iblk, :))
      END DO
      CPASSERT(.NOT. dbt_iterator_blocks_left(iter))
      CALL dbt_iterator_stop(iter)

      CALL dbt_reserve_blocks(tensor_out, blk_ind)
!$OMP END PARALLEL

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief reserve tensor blocks using matrix template
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_reserve_blocks_matrix_to_tensor(matrix_in, tensor_out)
      TYPE(dbcsr_type), TARGET, INTENT(IN)    :: matrix_in
      TYPE(dbt_type), INTENT(INOUT)  :: tensor_out
      TYPE(dbcsr_type), POINTER               :: matrix_in_desym

      INTEGER                            :: blk, iblk, nblk, nblk_per_thread, a, b
      INTEGER, ALLOCATABLE, DIMENSION(:) :: blk_ind_1, blk_ind_2
      INTEGER, DIMENSION(2)              :: ind_2d
      TYPE(dbcsr_iterator_type)          :: iter
      INTEGER                            :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_reserve_blocks_matrix_to_tensor'

      CALL timeset(routineN, handle)

      IF (dbcsr_has_symmetry(matrix_in)) THEN
         ALLOCATE (matrix_in_desym)
         CALL dbcsr_desymmetrize(matrix_in, matrix_in_desym)
      ELSE
         matrix_in_desym => matrix_in
      END IF

      nblk = dbcsr_get_num_blocks(matrix_in_desym)
      ALLOCATE (blk_ind_1(nblk), blk_ind_2(nblk))
      CALL dbcsr_iterator_start(iter, matrix_in_desym)
      DO iblk = 1, nblk
         CALL dbcsr_iterator_next_block(iter, ind_2d(1), ind_2d(2), blk)
         blk_ind_1(iblk) = ind_2d(1); blk_ind_2(iblk) = ind_2d(2)
      END DO
      CALL dbcsr_iterator_stop(iter)

!TODO: Parallelize creation of block list.
!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor_out,nblk,blk_ind_1,blk_ind_2) &
!$OMP PRIVATE(nblk_per_thread,a,b)
      nblk_per_thread = nblk/omp_get_num_threads() + 1
      a = omp_get_thread_num()*nblk_per_thread + 1
      b = MIN(a + nblk_per_thread, nblk)
      CALL dbt_reserve_blocks(tensor_out, blk_ind_1(a:b), blk_ind_2(a:b))
!$OMP END PARALLEL

      IF (dbcsr_has_symmetry(matrix_in)) THEN
         CALL dbcsr_release(matrix_in_desym)
         DEALLOCATE (matrix_in_desym)
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief reserve matrix blocks using tensor template
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_reserve_blocks_tensor_to_matrix(tensor_in, matrix_out)
      TYPE(dbt_type), INTENT(IN)        :: tensor_in
      TYPE(dbcsr_type), INTENT(INOUT)            :: matrix_out
      TYPE(dbt_iterator_type)           :: iter
      INTEGER, ALLOCATABLE, DIMENSION(:) :: blk_ind_1, blk_ind_2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_reserve_blocks_tensor_to_matrix'
      INTEGER :: handle, iblk, nblk
      INTEGER, DIMENSION(2)              :: ind_2d

      CALL timeset(routineN, handle)

      nblk = dbt_get_num_blocks(tensor_in)
      ALLOCATE (blk_ind_1(nblk), blk_ind_2(nblk))
      iblk = 0

!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor_in,matrix_out,iblk,blk_ind_1,blk_ind_2) &
!$OMP PRIVATE(iter,ind_2d)
      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_2d)
         IF (dbcsr_has_symmetry(matrix_out)) THEN
            IF (checker_tr(ind_2d(1), ind_2d(2))) CYCLE
            IF (ind_2d(1) > ind_2d(2)) CALL swap(ind_2d(1), ind_2d(2))
         END IF
!$OMP CRITICAL
         iblk = iblk + 1
         blk_ind_1(iblk) = ind_2d(1)
         blk_ind_2(iblk) = ind_2d(2)
!$OMP END CRITICAL
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_reserve_blocks(matrix_out, blk_ind_1(:iblk), blk_ind_2(:iblk))
      CALL dbcsr_finalize(matrix_out)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Swaps two integers
!> \author Patrick Seewald
! **************************************************************************************************
   ELEMENTAL SUBROUTINE swap(a, b)
      INTEGER, INTENT(INOUT)                             :: a, b
      INTEGER                                            :: tmp

      tmp = a
      a = b
      b = tmp
   END SUBROUTINE swap

! **************************************************************************************************
!> \brief Create block from array, array can be n-dimensional.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE create_block_data(block, sizes, array)
      TYPE(block_nd), INTENT(OUT)                       :: block
      INTEGER, DIMENSION(:), INTENT(IN)                 :: sizes
      REAL(dp), DIMENSION(PRODUCT(sizes)), INTENT(IN) :: array

      ALLOCATE (block%sizes, source=sizes)
      ALLOCATE (block%blk, source=array)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create and allocate block, but no data.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE create_block_nodata(block, sizes)
      INTEGER, INTENT(IN), DIMENSION(:)       :: sizes
      TYPE(block_nd), INTENT(OUT) :: block
      ALLOCATE (block%sizes, source=sizes)
      ALLOCATE (block%blk(PRODUCT(sizes)))
   END SUBROUTINE

! **************************************************************************************************
!> \brief
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE destroy_block(block)
      TYPE(block_nd), INTENT(INOUT) :: block
      DEALLOCATE (block%blk)
      DEALLOCATE (block%sizes)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Determines whether a transpose must be applied
!> \param row The absolute matrix row.
!> \param column The absolute matrix column
!> \param
!> \param
!> \param
!> \param
!> \param
!> \param
!> \author Patrick Seewald
! **************************************************************************************************
   ELEMENTAL FUNCTION checker_tr(row, column) RESULT(transpose)
      INTEGER, INTENT(IN)                                :: row, column
      LOGICAL                                            :: transpose

      transpose = BTEST(column + row, 0) .EQV. column .GE. row

   END FUNCTION checker_tr

! **************************************************************************************************
!> \brief Generic implementation of dbt_put_block, template for datatype
!> \param block block to put
!> \param summation whether block should be summed to existing block
!> \param ind block index
!> \param
!> \param
!> \param
!> \param
!> \param
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_put_anyd_block(tensor, ind, block, summation)
      TYPE(block_nd), INTENT(IN)       :: block
      TYPE(dbt_type), INTENT(INOUT)            :: tensor
      LOGICAL, INTENT(IN), OPTIONAL                :: summation
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(IN)                                :: ind

      SELECT CASE (ndims_tensor(tensor))
         #:for ndim in ndims
            CASE (${ndim}$)
            CALL dbt_put_${ndim}$d_block(tensor, ind, block%sizes, block%blk, summation=summation)
         #:endfor
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generic implementation of dbt_get_block (arbitrary tensor rank)
!> \param block block to get
!> \param found whether block was found
!> \param ind block index
!> \param
!> \param
!> \param
!> \param
!> \param
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_get_anyd_block(tensor, ind, block, found)
      TYPE(block_nd), INTENT(OUT)                  :: block
      LOGICAL, INTENT(OUT)                         :: found
      TYPE(dbt_type), INTENT(INOUT)            :: tensor
      INTEGER, DIMENSION(ndims_tensor(tensor)), &
         INTENT(IN)                                :: ind
      INTEGER, DIMENSION(ndims_tensor(tensor))    :: blk_size
      REAL(dp), DIMENSION(:), ALLOCATABLE         :: block_arr

      CALL dbt_blk_sizes(tensor, ind, blk_size)
      ALLOCATE (block_arr(PRODUCT(blk_size)))

      SELECT CASE (ndims_tensor(tensor))
         #:for ndim in ndims
            CASE (${ndim}$)
            CALL dbt_get_${ndim}$d_block(tensor, ind, blk_size, block_arr, found)
         #:endfor
      END SELECT
      CALL create_block(block, blk_size, block_arr)
   END SUBROUTINE

   #:for ndim in ndims
! **************************************************************************************************
!> \brief Template for dbt_put_block.
!> \param ind block index
!> \param sizes block size
!> \param block block to put
!> \param summation whether block should be summed to existing block
!> \param
!> \param
!> \param
!> \param
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE dbt_put_${ndim}$d_block(tensor, ind, sizes, block, summation)
         TYPE(dbt_type), INTENT(INOUT)                     :: tensor
         INTEGER, DIMENSION(${ndim}$), INTENT(IN) :: ind
         INTEGER, DIMENSION(${ndim}$), INTENT(IN) :: sizes
         REAL(dp), DIMENSION(${arrlist("sizes", nmax=ndim)}$), &
            INTENT(IN), TARGET                                 :: block
         LOGICAL, INTENT(IN), OPTIONAL                         :: summation
         INTEGER(KIND=int_8), DIMENSION(2)                     :: ind_2d
         INTEGER, DIMENSION(2)                                 :: shape_2d
         REAL(dp), POINTER, DIMENSION(:, :)                   :: block_2d
         INTEGER, DIMENSION(${ndim}$)                          :: shape_nd
         LOGICAL :: found, new_block
         REAL(dp), DIMENSION(${arrlist("sizes", nmax=ndim)}$) :: block_check

         LOGICAL, PARAMETER :: debug = .FALSE.
         INTEGER :: i

         new_block = .FALSE.

         IF (debug) THEN
            CALL dbt_get_block(tensor, ind, sizes, block_check, found=found)
            CPASSERT(found)
         END IF

         ASSOCIATE (map_nd => tensor%nd_index_blk%map_nd, &
                    map1_2d => tensor%nd_index_blk%map1_2d, &
                    map2_2d => tensor%nd_index_blk%map2_2d)

            shape_2d = [PRODUCT(sizes(map1_2d)), PRODUCT(sizes(map2_2d))]

            IF (ALL([map1_2d, map2_2d] == (/(i, i=1, ${ndim}$)/))) THEN
               ! to avoid costly reshape can do pointer bounds remapping as long as arrays are equivalent in memory
               block_2d(1:shape_2d(1), 1:shape_2d(2)) => block(${shape_colon(ndim)}$)
            ELSE
               ! need reshape due to rank reordering
               ALLOCATE (block_2d(shape_2d(1), shape_2d(2)))
               new_block = .TRUE.
               shape_nd(map_nd) = sizes
               block_2d(:, :) = RESHAPE(RESHAPE(block, SHAPE=shape_nd, order=map_nd), SHAPE=shape_2d)
            END IF

            ind_2d(:) = get_2d_indices_tensor(tensor%nd_index_blk, ind)

         END ASSOCIATE

         CALL dbt_tas_put_block(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d, summation=summation)

         IF (new_block) DEALLOCATE (block_2d)

      END SUBROUTINE
   #:endfor

   #:for ndim in ndims
! **************************************************************************************************
!> \brief allocate and get block
!> \param ind block index
!> \param block block to get
!> \param found whether block was found
!> \param
!> \param
!> \param
!> \param
!> \param
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE dbt_allocate_and_get_${ndim}$d_block(tensor, ind, block, found)
         TYPE(dbt_type), INTENT(INOUT)                     :: tensor
         INTEGER, DIMENSION(${ndim}$), INTENT(IN)  :: ind
         REAL(dp), DIMENSION(${shape_colon(ndim)}$), &
            ALLOCATABLE, INTENT(OUT)                           :: block
         LOGICAL, INTENT(OUT)                                  :: found
         INTEGER, DIMENSION(${ndim}$)                          :: blk_size

         CALL dbt_blk_sizes(tensor, ind, blk_size)
         CALL allocate_any(block, shape_spec=blk_size)
         CALL dbt_get_${ndim}$d_block(tensor, ind, blk_size, block, found)

      END SUBROUTINE
   #:endfor

   #:for ndim in ndims
! **************************************************************************************************
!> \brief Template for dbt_get_block.
!> \param ind block index
!> \param sizes block size
!> \param block block to get
!> \param found whether block was found
!> \author Patrick Seewald
! **************************************************************************************************
      SUBROUTINE dbt_get_${ndim}$d_block(tensor, ind, sizes, block, found)
         TYPE(dbt_type), INTENT(INOUT)                     :: tensor
         INTEGER, DIMENSION(${ndim}$), INTENT(IN) :: ind
         INTEGER, DIMENSION(${ndim}$), INTENT(IN) :: sizes
         REAL(dp), DIMENSION(${arrlist("sizes", nmax=ndim)}$), &
            INTENT(OUT)                                        :: block
         LOGICAL, INTENT(OUT)                                  :: found

         INTEGER(KIND=int_8), DIMENSION(2)                     :: ind_2d
         REAL(dp), DIMENSION(:, :), POINTER, CONTIGUOUS       :: block_2d_ptr
         INTEGER                                               :: i
         REAL(dp), DIMENSION(${shape_colon(ndim)}$), POINTER  :: block_ptr

         NULLIFY (block_2d_ptr)

         ind_2d(:) = get_2d_indices_tensor(tensor%nd_index_blk, ind)

         ASSOCIATE (map1_2d => tensor%nd_index_blk%map1_2d, &
                    map2_2d => tensor%nd_index_blk%map2_2d)

            CALL dbt_tas_get_block_p(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d_ptr)
            found = ASSOCIATED(block_2d_ptr)

            IF (found) THEN
               IF (ALL([map1_2d, map2_2d] == (/(i, i=1, ${ndim}$)/))) THEN
                  ! to avoid costly reshape can do pointer bounds remapping as long as arrays are equivalent in memory
                  block_ptr(${shape_explicit('block', ndim)}$) => block_2d_ptr(:, :)
                  block(${shape_colon(ndim)}$) = block_ptr(${shape_colon(ndim)}$)
               ELSE
                  ! need reshape due to rank reordering
                  block(${shape_colon(ndim)}$) = RESHAPE(block_2d_ptr, SHAPE=SHAPE(block), ORDER=[map1_2d, map2_2d])
               END IF
            END IF

         END ASSOCIATE

      END SUBROUTINE
   #:endfor

END MODULE
