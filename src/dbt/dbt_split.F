!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to split blocks and to convert between tensors with different block sizes.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_split

   #:include "dbt_macros.fypp"
   #:set maxdim = maxrank
   #:set ndims = range(2,maxdim+1)

   USE dbt_allocate_wrap, ONLY: allocate_any
   USE dbt_array_list_methods, ONLY: get_ith_array
   USE dbt_block, ONLY: dbt_iterator_type, &
                        dbt_get_block, &
                        dbt_put_block, &
                        dbt_iterator_start, &
                        dbt_iterator_blocks_left, &
                        dbt_iterator_stop, &
                        dbt_iterator_next_block, &
                        dbt_reserve_blocks, &
                        dbt_reserved_block_indices
   USE dbt_index, ONLY: dbt_get_mapping_info, &
                        dbt_inverse_order
   USE dbt_types, ONLY: dbt_create, &
                        dbt_type, &
                        ndims_tensor, &
                        dbt_distribution_type, &
                        dbt_distribution, &
                        dbt_distribution_destroy, &
                        dbt_distribution_new_expert, &
                        dbt_clear, &
                        dbt_finalize, &
                        dbt_get_num_blocks, &
                        dbt_blk_offsets, &
                        dbt_blk_sizes, &
                        ndims_matrix_row, &
                        ndims_matrix_column, &
                        dbt_filter, &
                        dbt_copy_contraction_storage
   USE kinds, ONLY: dp, dp

#include "../base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_split'
   PUBLIC :: &
      dbt_make_compatible_blocks, &
      dbt_split_blocks, &
      dbt_split_blocks_generic, &
      dbt_split_copyback, &
      dbt_crop

CONTAINS

! **************************************************************************************************
!> \brief Split tensor blocks into smaller blocks
!> \param tensor_in Input tensor
!> \param tensor_out Output tensor (splitted blocks)
!> \param blk_size_i block sizes for each of the tensor dimensions
!> \param nodata don't copy data from tensor_in to tensor_out
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_split_blocks_generic(tensor_in, tensor_out, ${varlist("blk_size")}$, nodata)
      TYPE(dbt_type), INTENT(INOUT)               :: tensor_in
      TYPE(dbt_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL     :: ${varlist("blk_size")}$
      LOGICAL, INTENT(IN), OPTIONAL                   :: nodata

      TYPE(dbt_distribution_type)                 :: dist_old, dist_split
      TYPE(dbt_iterator_type)                     :: iter
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("nd_dist_split")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("nd_blk_size_split")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("index_split_offset")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("inblock_offset")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("blk_nsplit")}$
      INTEGER                                         :: ${varlist("split_blk")}$
      INTEGER :: idim, i, isplit_sum, nsplit, handle, splitsum, bcount
      INTEGER, DIMENSION(:, :), ALLOCATABLE           :: blks_to_allocate
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist_d, blk_size_d, blk_size_split_d, dist_split_d
      INTEGER, DIMENSION(ndims_matrix_row(tensor_in)) :: map1_2d
      INTEGER, DIMENSION(ndims_matrix_column(tensor_in)) :: map2_2d
      INTEGER, DIMENSION(ndims_tensor(tensor_in)) :: blk_index, blk_size, blk_offset, &
                                                     blk_shape
      INTEGER, DIMENSION(${maxdim}$) :: bi_split, inblock_offset
      LOGICAL :: found

      #:for ndim in ndims
         REAL(dp), DIMENSION(${shape_colon(n=ndim)}$), ALLOCATABLE :: block_${ndim}$d
      #:endfor
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_split_blocks_generic'

      CALL timeset(routineN, handle)

      dist_old = dbt_distribution(tensor_in)

      DO idim = 1, ndims_tensor(tensor_in)
         CALL get_ith_array(dist_old%nd_dist, idim, dist_d)
         CALL get_ith_array(tensor_in%blk_sizes, idim, blk_size_d)

         #:for idim in range(1, maxdim+1)
            IF (idim == ${idim}$) THEN
               ! split block index offset for each normal block index
               ALLOCATE (index_split_offset_${idim}$ (SIZE(dist_d)))
               ! how many split blocks for each normal block index
               ALLOCATE (blk_nsplit_${idim}$ (SIZE(dist_d)))
               ! data offset of split blocks w.r.t. corresponding normal block
               ALLOCATE (inblock_offset_${idim}$ (SIZE(blk_size_${idim}$)))
               ALLOCATE (blk_size_split_d, source=blk_size_${idim}$)
            END IF
         #:endfor

         ! distribution vector for split blocks
         ALLOCATE (dist_split_d(SIZE(blk_size_split_d)))

         isplit_sum = 0 ! counting splits
         DO i = 1, SIZE(blk_size_d)
            nsplit = 0 ! number of splits for current normal block
            splitsum = 0 ! summing split block sizes for current normal block
            DO WHILE (splitsum < blk_size_d(i))
               nsplit = nsplit + 1
               isplit_sum = isplit_sum + 1
               #:for idim in range(1, maxdim+1)
                  IF (idim == ${idim}$) inblock_offset_${idim}$ (isplit_sum) = splitsum
               #:endfor
               dist_split_d(isplit_sum) = dist_d(i)
               splitsum = splitsum + blk_size_split_d(isplit_sum)
            END DO
            CPASSERT(splitsum == blk_size_d(i))
            #:for idim in range(1, maxdim+1)
               IF (idim == ${idim}$) THEN
                  blk_nsplit_${idim}$ (i) = nsplit
                  index_split_offset_${idim}$ (i) = isplit_sum - nsplit
               END IF
            #:endfor
         END DO

         #:for idim in range(1, maxdim+1)
            IF (idim == ${idim}$) THEN
               ALLOCATE (nd_dist_split_${idim}$, source=dist_split_d)
               ALLOCATE (nd_blk_size_split_${idim}$, source=blk_size_split_d)
            END IF
         #:endfor
         DEALLOCATE (dist_split_d)
         DEALLOCATE (blk_size_split_d)

      END DO

      CALL dbt_get_mapping_info(tensor_in%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)

      #:for ndim in ndims
         IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
            CALL dbt_distribution_new_expert(dist_split, tensor_in%pgrid, map1_2d, map2_2d, &
                                             ${varlist("nd_dist_split", nmax=ndim)}$)
            CALL dbt_create(tensor_out, tensor_in%name, dist_split, map1_2d, map2_2d, &
                            ${varlist("nd_blk_size_split", nmax=ndim)}$)
         END IF
      #:endfor

      CALL dbt_distribution_destroy(dist_split)

      IF (PRESENT(nodata)) THEN
         IF (nodata) THEN
            CALL timestop(handle)
            RETURN
         END IF
      END IF

      CALL dbt_iterator_start(iter, tensor_in)

      bcount = 0
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size)
         #:for ndim in ndims
            IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
               #:for idim in range(1,ndim+1)
                  DO split_blk_${idim}$ = 1, blk_nsplit_${idim}$ (blk_index(${idim}$))
                     #:endfor
                     bcount = bcount + 1
                     #:for idim in range(1,ndim+1)
                        END DO
                     #:endfor
                  END IF
               #:endfor
            END DO
            CALL dbt_iterator_stop(iter)

            ALLOCATE (blks_to_allocate(bcount, ndims_tensor(tensor_in)))

            CALL dbt_iterator_start(iter, tensor_in)

            bcount = 0
            DO WHILE (dbt_iterator_blocks_left(iter))
               CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)

               #:for ndim in ndims
                  IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
                     #:for idim in range(1,ndim+1)
                        DO split_blk_${idim}$ = 1, blk_nsplit_${idim}$ (blk_index(${idim}$))
                           bi_split(${idim}$) = index_split_offset_${idim}$ (blk_index(${idim}$)) + split_blk_${idim}$
                           #:endfor
                           bcount = bcount + 1
                           blks_to_allocate(bcount, :) = bi_split(1:ndims_tensor(tensor_in))
                           #:for idim in range(1,ndim+1)
                              END DO
                           #:endfor
                        END IF
                     #:endfor
                  END DO

                  CALL dbt_iterator_stop(iter)

                  CALL dbt_reserve_blocks(tensor_out, blks_to_allocate)

                  CALL dbt_iterator_start(iter, tensor_in)

                  DO WHILE (dbt_iterator_blocks_left(iter))
                     CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)
                     #:for ndim in ndims
                        IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
                           CALL dbt_get_block(tensor_in, blk_index, block_${ndim}$d, found)
                           CPASSERT(found)
                        END IF
                     #:endfor
                     #:for ndim in ndims
                        IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
                           #:for idim in range(1,ndim+1)
                              DO split_blk_${idim}$ = 1, blk_nsplit_${idim}$ (blk_index(${idim}$))
                                 ! split block index
                                 bi_split(${idim}$) = index_split_offset_${idim}$ (blk_index(${idim}$)) + split_blk_${idim}$
                                 blk_shape(${idim}$) = blk_size_${idim}$ (bi_split(${idim}$))
                                 #:endfor

                                 #:for idim in range(1,ndim+1)
                                    inblock_offset(${idim}$) = inblock_offset_${idim}$ (bi_split(${idim}$))
                                 #:endfor
                                 CALL dbt_put_block(tensor_out, bi_split(1:${ndim}$), &
                                                    blk_shape, &
                                                    block_${ndim}$d( &
                                      ${", ".join(["inblock_offset("+str(idim)+") + 1:inblock_offset("+str(idim)+") + blk_shape("+str(idim)+")" for idim in range(1, ndim+1)])}$))

                                 #:for idim in range(1,ndim+1)
                                    END DO
                                 #:endfor

                                 DEALLOCATE (block_${ndim}$d)
                              END IF
                           #:endfor
                        END DO
                        CALL dbt_iterator_stop(iter)

                        CALL dbt_finalize(tensor_out)

                        ! remove blocks that are exactly 0, these can occur if a cropping operation was performed before splitting
                        CALL dbt_filter(tensor_out, TINY(0.0_dp))

                        CALL timestop(handle)

                     END SUBROUTINE

! **************************************************************************************************
!> \brief Split tensor blocks into smaller blocks of maximum size PRODUCT(block_sizes).
!> \param tensor_in Input tensor
!> \param tensor_out Output tensor (split blocks)
!> \param block_sizes block sizes for each of the tensor dimensions
!> \param nodata don't copy data from tensor_in to tensor_out
!> \author Patrick Seewald
! **************************************************************************************************
                     SUBROUTINE dbt_split_blocks(tensor_in, tensor_out, block_sizes, nodata)

                        TYPE(dbt_type), INTENT(INOUT)               :: tensor_in
                        TYPE(dbt_type), INTENT(OUT)                 :: tensor_out
                        INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
                           INTENT(IN)                                   :: block_sizes
                        LOGICAL, INTENT(IN), OPTIONAL                   :: nodata

                        INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("nd_blk_size_split")}$
                        INTEGER :: idim, i, isplit_sum, blk_remainder, nsplit, isplit
                        INTEGER, DIMENSION(:), ALLOCATABLE :: blk_size_d, blk_size_split_d

                        DO idim = 1, ndims_tensor(tensor_in)
                           CALL get_ith_array(tensor_in%blk_sizes, idim, blk_size_d)

                           isplit_sum = 0
                           DO i = 1, SIZE(blk_size_d)
                              nsplit = (blk_size_d(i) + block_sizes(idim) - 1)/block_sizes(idim)
                              isplit_sum = isplit_sum + nsplit
                           END DO

                           ALLOCATE (blk_size_split_d(isplit_sum))

                           isplit_sum = 0
                           DO i = 1, SIZE(blk_size_d)
                              nsplit = (blk_size_d(i) + block_sizes(idim) - 1)/block_sizes(idim)
                              blk_remainder = blk_size_d(i)
                              DO isplit = 1, nsplit
                                 isplit_sum = isplit_sum + 1
                                 blk_size_split_d(isplit_sum) = MIN(block_sizes(idim), blk_remainder)
                                 blk_remainder = blk_remainder - block_sizes(idim)
                              END DO

                           END DO

                           #:for idim in range(1, maxdim+1)
                              IF (idim == ${idim}$) THEN
                                 ALLOCATE (nd_blk_size_split_${idim}$, source=blk_size_split_d)
                              END IF
                           #:endfor
                           DEALLOCATE (blk_size_split_d)
                        END DO

                        #:for ndim in ndims
                           IF (ndims_tensor(tensor_in) == ${ndim}$) CALL dbt_split_blocks_generic(tensor_in, tensor_out, &
                                                                                      ${varlist("nd_blk_size_split", nmax=ndim)}$, &
                                                                                                  nodata=nodata)
                        #:endfor

                     END SUBROUTINE

! **************************************************************************************************
!> \brief Copy tensor with split blocks to tensor with original block sizes.
!> \param tensor_split_in tensor with smaller blocks
!> \param tensor_out original tensor
!> \author Patrick Seewald
! **************************************************************************************************
                     SUBROUTINE dbt_split_copyback(tensor_split_in, tensor_out, summation)
                        TYPE(dbt_type), INTENT(INOUT)               :: tensor_split_in
                        TYPE(dbt_type), INTENT(INOUT)               :: tensor_out
                        LOGICAL, INTENT(IN), OPTIONAL                   :: summation
                        INTEGER, DIMENSION(:), ALLOCATABLE              :: first_split_d, last_split_d
                        INTEGER, DIMENSION(:), ALLOCATABLE              :: blk_size_split_d, blk_size_d
                        INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("last_split")}$, &
                                                                           ${varlist("first_split")}$, &
                                                                           ${varlist("split")}$
                     INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("inblock_offset")}$, ${varlist("blk_size_split")}$
                        INTEGER, DIMENSION(:, :), ALLOCATABLE            :: blks_to_allocate
                        INTEGER                                         :: idim, iblk, bcount
                        INTEGER                                         :: ${varlist("iblk")}$, isplit_sum, splitsum, nblk
                        TYPE(dbt_iterator_type)                     :: iter
                        INTEGER, DIMENSION(ndims_tensor(tensor_out)) :: blk_index, blk_size, blk_offset, blk_shape, blk_index_n
                        LOGICAL                                         :: found

                        INTEGER, DIMENSION(${maxdim}$)                  :: inblock_offset
                        INTEGER                                            :: handle
                        CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_split_copyback'
                        #:for ndim in ndims
                           REAL(dp), DIMENSION(${shape_colon(n=ndim)}$), ALLOCATABLE :: block_${ndim}$d
                           REAL(dp), DIMENSION(${shape_colon(n=ndim)}$), ALLOCATABLE :: block_split_${ndim}$d
                        #:endfor

                        CALL timeset(routineN, handle)
                        CPASSERT(tensor_out%valid)
                        IF (PRESENT(summation)) THEN
                           IF (.NOT. summation) CALL dbt_clear(tensor_out)
                        ELSE
                           CALL dbt_clear(tensor_out)
                        END IF

                        DO idim = 1, ndims_tensor(tensor_split_in)
                           CALL get_ith_array(tensor_split_in%blk_sizes, idim, blk_size_split_d)
                           CALL get_ith_array(tensor_out%blk_sizes, idim, blk_size_d)

                           #:for idim in range(1, maxdim+1)
                              IF (idim == ${idim}$) THEN
                                 ! data offset of split blocks w.r.t. corresponding normal block
                                 ALLOCATE (inblock_offset_${idim}$ (SIZE(blk_size_split_d)))
                                 ! normal block index for each split block
                                 ALLOCATE (split_${idim}$ (SIZE(blk_size_split_d)))
                              END IF
                           #:endfor

                           ALLOCATE (last_split_d(SIZE(blk_size_d)))
                           ALLOCATE (first_split_d(SIZE(blk_size_d)))
                           first_split_d(1) = 1
                           isplit_sum = 0
                           DO iblk = 1, SIZE(blk_size_d)
                              splitsum = 0
                              IF (iblk .GT. 1) first_split_d(iblk) = last_split_d(iblk - 1) + 1
                              DO WHILE (splitsum < blk_size_d(iblk))
                                 isplit_sum = isplit_sum + 1
                                 #:for idim in range(1, maxdim+1)
                                    IF (idim == ${idim}$) THEN
                                       inblock_offset_${idim}$ (isplit_sum) = splitsum
                                       split_${idim}$ (isplit_sum) = iblk
                                    END IF
                                 #:endfor
                                 splitsum = splitsum + blk_size_split_d(isplit_sum)
                              END DO
                              CPASSERT(splitsum == blk_size_d(iblk))
                              last_split_d(iblk) = isplit_sum
                           END DO
                           #:for idim in range(1, maxdim+1)
                              IF (idim == ${idim}$) THEN
                                 ALLOCATE (first_split_${idim}$, source=first_split_d)
                                 ALLOCATE (last_split_${idim}$, source=last_split_d)
                                 ALLOCATE (blk_size_split_${idim}$, source=blk_size_split_d)
                              END IF
                           #:endfor
                           DEALLOCATE (first_split_d, last_split_d)
                           DEALLOCATE (blk_size_split_d, blk_size_d)
                        END DO

                        nblk = dbt_get_num_blocks(tensor_split_in)
                        ALLOCATE (blks_to_allocate(nblk, ndims_tensor(tensor_split_in)))
                        CALL dbt_iterator_start(iter, tensor_split_in)
                        bcount = 0
                        DO WHILE (dbt_iterator_blocks_left(iter))
                           CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size)
                           #:for ndim in ndims
                              IF (ndims_tensor(tensor_out) == ${ndim}$) THEN
                                 #:for idim in range(1,ndim+1)
                                    blk_index_n(${idim}$) = split_${idim}$ (blk_index(${idim}$))
                                 #:endfor
                              END IF
                           #:endfor
                           blks_to_allocate(bcount + 1, :) = blk_index_n
                           bcount = bcount + 1
                        END DO
                        CALL dbt_iterator_stop(iter)
                        CALL dbt_reserve_blocks(tensor_out, blks_to_allocate)

                        CALL dbt_iterator_start(iter, tensor_out)
                        DO WHILE (dbt_iterator_blocks_left(iter))
                           CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)

                           #:for ndim in ndims
                              IF (ndims_tensor(tensor_out) == ${ndim}$) THEN
                                 CALL allocate_any(block_${ndim}$d, blk_size)
                                 block_${ndim}$d = 0.0_dp
                                 #:for idim in range(1,ndim+1)
                            DO iblk_${idim}$ = first_split_${idim}$ (blk_index(${idim}$)), last_split_${idim}$ (blk_index(${idim}$))
                                       #:endfor
                                       #:for idim in range(1,ndim+1)
                                          inblock_offset(${idim}$) = inblock_offset_${idim}$ (iblk_${idim}$)
                                       #:endfor

                            CALL dbt_get_block(tensor_split_in, [${", ".join(["iblk_"+str(idim) for idim in range(1, ndim+1)])}$], &
                                                          block_split_${ndim}$d, found)
                                       IF (found) THEN
                                          blk_shape(1:${ndim}$) = SHAPE(block_split_${ndim}$d)
                                          block_${ndim}$d( &
                     ${", ".join(["inblock_offset("+str(idim)+") + 1:inblock_offset("+str(idim)+") + blk_shape("+str(idim)+")" for idim in range(1, ndim+1)])}$) = &
                                             block_split_${ndim}$d
                                       END IF

                                       #:for idim in range(1,ndim+1)
                                          END DO
                                       #:endfor
                                       CALL dbt_put_block(tensor_out, blk_index, blk_size, block_${ndim}$d, summation=summation)
                                       DEALLOCATE (block_${ndim}$d)
                                    END IF
                                 #:endfor
                              END DO
                              CALL dbt_iterator_stop(iter)

                              CALL timestop(handle)

                           END SUBROUTINE

! **************************************************************************************************
!> \brief split two tensors with same total sizes but different block sizes such that they have
!>        equal block sizes
!> \param move_data memory optimization: transfer data s.t. tensor1 and tensor2 may be empty on return
!> \param tensor1_split tensor 1 with split blocks
!> \param tensor2_split tensor 2 with split blocks
!> \param nodata1 don't copy data of tensor 1
!> \param nodata2 don't copy data of tensor 2
!> \param
!> \param
!> \param
!> \param
!> \author Patrick Seewald
! **************************************************************************************************
                           SUBROUTINE dbt_make_compatible_blocks(tensor1, tensor2, tensor1_split, tensor2_split, &
                                                                 order, nodata1, nodata2, move_data)
                              TYPE(dbt_type), INTENT(INOUT) :: tensor1, tensor2
                              TYPE(dbt_type), INTENT(OUT)   :: tensor1_split, tensor2_split
                              INTEGER, DIMENSION(ndims_tensor(tensor1)), &
                                 INTENT(IN), OPTIONAL                        :: order
                              LOGICAL, INTENT(IN), OPTIONAL     :: nodata1, nodata2, move_data
                          INTEGER, DIMENSION(:), ALLOCATABLE  :: ${varlist("blk_size_split_1")}$, ${varlist("blk_size_split_2")}$, &
                                                                     blk_size_d_1, blk_size_d_2, blk_size_d_split
                              INTEGER :: size_sum_1, size_sum_2, size_sum, bind_1, bind_2, isplit, bs, idim, i
                              LOGICAL :: move_prv, nodata1_prv, nodata2_prv
                              INTEGER, DIMENSION(ndims_tensor(tensor1)) :: order_prv

                              IF (PRESENT(move_data)) THEN
                                 move_prv = move_data
                              ELSE
                                 move_prv = .FALSE.
                              END IF

                              IF (PRESENT(nodata1)) THEN
                                 nodata1_prv = nodata1
                              ELSE
                                 nodata1_prv = .FALSE.
                              END IF
                              IF (PRESENT(nodata2)) THEN
                                 nodata2_prv = nodata2
                              ELSE
                                 nodata2_prv = .FALSE.
                              END IF

                              IF (PRESENT(order)) THEN
                                 order_prv(:) = dbt_inverse_order(order)
                              ELSE
                                 order_prv(:) = (/(i, i=1, ndims_tensor(tensor1))/)
                              END IF

                              DO idim = 1, ndims_tensor(tensor2)
                                 CALL get_ith_array(tensor1%blk_sizes, order_prv(idim), blk_size_d_1)
                                 CALL get_ith_array(tensor2%blk_sizes, idim, blk_size_d_2)
                                 ALLOCATE (blk_size_d_split(SIZE(blk_size_d_1) + SIZE(blk_size_d_2)))
                                 size_sum_1 = 0
                                 size_sum_2 = 0
                                 size_sum = 0
                                 bind_1 = 0
                                 bind_2 = 0
                                 isplit = 0

                                 DO WHILE (bind_1 < SIZE(blk_size_d_1) .AND. bind_2 < SIZE(blk_size_d_2))
                                    IF (blk_size_d_1(bind_1 + 1) < blk_size_d_2(bind_2 + 1)) THEN
                                       bind_1 = bind_1 + 1
                                       bs = blk_size_d_1(bind_1)
                                       blk_size_d_2(bind_2 + 1) = blk_size_d_2(bind_2 + 1) - bs
                                       size_sum = size_sum + bs
                                       isplit = isplit + 1
                                       blk_size_d_split(isplit) = bs
                                    ELSEIF (blk_size_d_1(bind_1 + 1) > blk_size_d_2(bind_2 + 1)) THEN
                                       bind_2 = bind_2 + 1
                                       bs = blk_size_d_2(bind_2)
                                       blk_size_d_1(bind_1 + 1) = blk_size_d_1(bind_1 + 1) - bs
                                       size_sum = size_sum + bs
                                       isplit = isplit + 1
                                       blk_size_d_split(isplit) = bs
                                    ELSE
                                       bind_1 = bind_1 + 1
                                       bind_2 = bind_2 + 1
                                       bs = blk_size_d_1(bind_1)
                                       size_sum = size_sum + bs
                                       isplit = isplit + 1
                                       blk_size_d_split(isplit) = bs
                                    END IF
                                 END DO

                                 IF (bind_1 < SIZE(blk_size_d_1)) THEN
                                    bind_1 = bind_1 + 1
                                    bs = blk_size_d_1(bind_1)
                                    size_sum = size_sum + bs
                                    isplit = isplit + 1
                                    blk_size_d_split(isplit) = bs
                                 END IF

                                 IF (bind_2 < SIZE(blk_size_d_2)) THEN
                                    bind_2 = bind_2 + 1
                                    bs = blk_size_d_2(bind_2)
                                    size_sum = size_sum + bs
                                    isplit = isplit + 1
                                    blk_size_d_split(isplit) = bs
                                 END IF

                                 #:for idim in range(1, maxdim+1)
                                    IF (order_prv(idim) == ${idim}$) THEN
                                       ALLOCATE (blk_size_split_1_${idim}$, source=blk_size_d_split(:isplit))
                                    END IF
                                 #:endfor

                                 #:for idim in range(1, maxdim+1)
                                    IF (idim == ${idim}$) THEN
                                       ALLOCATE (blk_size_split_2_${idim}$, source=blk_size_d_split(:isplit))
                                    END IF
                                 #:endfor

                                 DEALLOCATE (blk_size_d_split, blk_size_d_1, blk_size_d_2)
                              END DO

                              #:for ndim in ndims
                                 IF (ndims_tensor(tensor1) == ${ndim}$) THEN
                   CALL dbt_split_blocks_generic(tensor1, tensor1_split, ${varlist("blk_size_split_1", nmax=ndim)}$, nodata=nodata1)
                                    IF (move_prv .AND. .NOT. nodata1_prv) CALL dbt_clear(tensor1)
                   CALL dbt_split_blocks_generic(tensor2, tensor2_split, ${varlist("blk_size_split_2", nmax=ndim)}$, nodata=nodata2)
                                    IF (move_prv .AND. .NOT. nodata2_prv) CALL dbt_clear(tensor2)
                                 END IF
                              #:endfor

                           END SUBROUTINE

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
                           SUBROUTINE dbt_crop(tensor_in, tensor_out, bounds, move_data)
                              TYPE(dbt_type), INTENT(INOUT) :: tensor_in
                              TYPE(dbt_type), INTENT(OUT) :: tensor_out
                              INTEGER, DIMENSION(2, ndims_tensor(tensor_in)), INTENT(IN) :: bounds
                              LOGICAL, INTENT(IN), OPTIONAL :: move_data
                              INTEGER, DIMENSION(2, ndims_tensor(tensor_in)) :: blk_bounds
                              TYPE(dbt_iterator_type)                     :: iter
                              INTEGER, DIMENSION(ndims_tensor(tensor_in)) :: blk_index, blk_size, blk_offset
                              LOGICAL :: found, move_data_prv
                              INTEGER :: idim, iblk, iblk_all, nblk
                              INTEGER, DIMENSION(:, :), ALLOCATABLE :: blk_ind, blk_ind_tmp
                              #:for ndim in ndims
                                 REAL(dp), DIMENSION(${shape_colon(n=ndim)}$), ALLOCATABLE :: block_${ndim}$d, block_put_${ndim}$d
                              #:endfor

                              IF (PRESENT(move_data)) THEN
                                 move_data_prv = move_data
                              ELSE
                                 move_data_prv = .FALSE.
                              END IF

                              CALL dbt_create(tensor_in, tensor_out)

                              ! reserve blocks inside bounds
                              ALLOCATE (blk_ind(dbt_get_num_blocks(tensor_in), ndims_tensor(tensor_in)))
                              CALL dbt_reserved_block_indices(tensor_in, blk_ind)
                              nblk = dbt_get_num_blocks(tensor_in)
                              ALLOCATE (blk_ind_tmp(dbt_get_num_blocks(tensor_in), ndims_tensor(tensor_in)))
                              blk_ind_tmp(:, :) = 0
                              iblk = 0
                              blk_loop: DO iblk_all = 1, nblk
                                 CALL dbt_blk_offsets(tensor_in, blk_ind(iblk_all, :), blk_offset)
                                 CALL dbt_blk_sizes(tensor_in, blk_ind(iblk_all, :), blk_size)
                                 DO idim = 1, ndims_tensor(tensor_in)
                                    IF (bounds(1, idim) > blk_offset(idim) - 1 + blk_size(idim)) CYCLE blk_loop
                                    IF (bounds(2, idim) < blk_offset(idim)) CYCLE blk_loop
                                 END DO
                                 iblk = iblk + 1
                                 blk_ind_tmp(iblk, :) = blk_ind(iblk_all, :)
                              END DO blk_loop

                              DEALLOCATE (blk_ind)
                              ALLOCATE (blk_ind(iblk, ndims_tensor(tensor_in)))
                              blk_ind(:, :) = blk_ind_tmp(:iblk, :)

                              CALL dbt_reserve_blocks(tensor_out, blk_ind)

                              ! copy blocks
                              CALL dbt_iterator_start(iter, tensor_out)
                              iter_loop: DO WHILE (dbt_iterator_blocks_left(iter))
                                 CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size, blk_offset=blk_offset)

                                 DO idim = 1, ndims_tensor(tensor_in)
                                    blk_bounds(1, idim) = MAX(bounds(1, idim) - blk_offset(idim) + 1, 1)
                                    blk_bounds(2, idim) = MIN(bounds(2, idim) - blk_offset(idim) + 1, blk_size(idim))
                                 END DO

                                 #:for ndim in ndims
                                    IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
                                       CALL dbt_get_block(tensor_in, blk_index, block_${ndim}$d, found)

                             ALLOCATE (block_put_${ndim}$d(${", ".join(["blk_size("+str(idim)+")" for idim in range(1, ndim+1)])}$))
                                       block_put_${ndim}$d = 0.0_dp
   block_put_${ndim}$d(${", ".join(["blk_bounds(1, "+str(idim)+"):blk_bounds(2,"+str(idim)+")" for idim in range(1, ndim+1)])}$) = &
           block_${ndim}$d(${", ".join(["blk_bounds(1, "+str(idim)+"):blk_bounds(2,"+str(idim)+")" for idim in range(1, ndim+1)])}$)
                                       CALL dbt_put_block(tensor_out, blk_index, blk_size, block_put_${ndim}$d)
                                       DEALLOCATE (block_${ndim}$d)
                                       DEALLOCATE (block_put_${ndim}$d)
                                    END IF
                                 #:endfor
                              END DO iter_loop
                              CALL dbt_iterator_stop(iter)
                              CALL dbt_finalize(tensor_out)

                              IF (move_data_prv) CALL dbt_clear(tensor_in)

                              ! transfer data for batched contraction
                              CALL dbt_copy_contraction_storage(tensor_in, tensor_out)

                           END SUBROUTINE

                        END MODULE
