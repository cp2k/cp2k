!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR tensor framework for block-sparse tensor contraction.
!>
!> Representation of n-rank (3 <= n <= 6) tensors as DBCSR matrices.
!> Support for arbitrary redistribution between different representations.
!> DBCSR routines are generalized to n dimensions by light-weight wrapper routines.
!> \todo  implement tensor contraction using DBCSR MM.
!>        implement communication optimal scheme for multiplication of non-square matrices.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsrt
   USE allocate_wrap,                   ONLY: allocate_any
   USE array_list_methods,              ONLY: array_list,&
                                              create_array_list,&
                                              destroy_array_list,&
                                              get_array_elements,&
                                              get_arrays,&
                                              sizes_of_arrays
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_distribution, dbcsr_distribution_mp, dbcsr_get_block_p, &
        dbcsr_get_stored_coordinates, dbcsr_iterator, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_mp_group, &
        dbcsr_put_block, dbcsr_reserve_blocks, dbcsr_type_complex_4, dbcsr_type_complex_8, &
        dbcsr_type_real_4, dbcsr_type_real_8
   USE dbcsrt_block,                    ONLY: &
        block_buffer_add_anyd_block, block_buffer_blocks_left, block_buffer_create, &
        block_buffer_get_next_anyd_block, block_buffer_type, block_nd, communicate_buffer, &
        create_block, destroy_block, ndims_buffer, reshape_2d_to_nd_block, reshape_nd_to_2d_block
   USE dbcsrt_block_c,                  ONLY: block_nd_c
   USE dbcsrt_block_d,                  ONLY: block_nd_d
   USE dbcsrt_block_s,                  ONLY: block_nd_s
   USE dbcsrt_block_z,                  ONLY: block_nd_z
   USE dbcsrt_index,                    ONLY: get_2d_indices,&
                                              get_nd_indices
   USE dbcsrt_index_types,              ONLY: create_nd_to_2d_mapping,&
                                              destroy_nd_to_2d_mapping,&
                                              get_mapping_info,&
                                              nd_to_2d_mapping,&
                                              ndims_mapping
   USE dbcsrt_types,                    ONLY: dbcsrt_blk_sizes,&
                                              dbcsrt_create,&
                                              dbcsrt_get_data_type,&
                                              dbcsrt_type,&
                                              ndims_tensor
   USE kinds,                           ONLY: real_4,&
                                              real_8
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_environ
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsrt'

   PUBLIC :: &
      dbcsrt_contract, &
      dbcsrt_copy, &
      dbcsrt_get_block, &
      dbcsrt_get_stored_coordinates, &
      dbcsrt_iterator, &
      dbcsrt_iterator_blocks_left, &
      dbcsrt_iterator_next_block, &
      dbcsrt_iterator_start, &
      dbcsrt_iterator_stop, &
      dbcsrt_ndims, &
      dbcsrt_permute_index, &
      dbcsrt_put_block, &
      dbcsrt_reserve_blocks

   TYPE dbcsrt_iterator
      TYPE(dbcsr_iterator)             :: iter
      TYPE(nd_to_2d_mapping)           :: nd_index_blk
      TYPE(nd_to_2d_mapping)           :: nd_index
      TYPE(array_list)                 :: blk_sizes, blk_offsets
   END TYPE dbcsrt_iterator

   INTERFACE dbcsrt_put_block
      MODULE PROCEDURE dbcsrt_put_3d_block_d
      MODULE PROCEDURE dbcsrt_put_3d_block_s
      MODULE PROCEDURE dbcsrt_put_3d_block_z
      MODULE PROCEDURE dbcsrt_put_3d_block_c
      MODULE PROCEDURE dbcsrt_put_4d_block_d
      MODULE PROCEDURE dbcsrt_put_4d_block_s
      MODULE PROCEDURE dbcsrt_put_4d_block_z
      MODULE PROCEDURE dbcsrt_put_4d_block_c
      MODULE PROCEDURE dbcsrt_put_5d_block_d
      MODULE PROCEDURE dbcsrt_put_5d_block_s
      MODULE PROCEDURE dbcsrt_put_5d_block_z
      MODULE PROCEDURE dbcsrt_put_5d_block_c
      MODULE PROCEDURE dbcsrt_put_6d_block_d
      MODULE PROCEDURE dbcsrt_put_6d_block_s
      MODULE PROCEDURE dbcsrt_put_6d_block_z
      MODULE PROCEDURE dbcsrt_put_6d_block_c
      MODULE PROCEDURE dbcsrt_put_anyd_block
   END INTERFACE

   INTERFACE dbcsrt_get_block
      MODULE PROCEDURE dbcsrt_get_3d_block_d
      MODULE PROCEDURE dbcsrt_get_3d_block_s
      MODULE PROCEDURE dbcsrt_get_3d_block_z
      MODULE PROCEDURE dbcsrt_get_3d_block_c
      MODULE PROCEDURE dbcsrt_get_4d_block_d
      MODULE PROCEDURE dbcsrt_get_4d_block_s
      MODULE PROCEDURE dbcsrt_get_4d_block_z
      MODULE PROCEDURE dbcsrt_get_4d_block_c
      MODULE PROCEDURE dbcsrt_get_5d_block_d
      MODULE PROCEDURE dbcsrt_get_5d_block_s
      MODULE PROCEDURE dbcsrt_get_5d_block_z
      MODULE PROCEDURE dbcsrt_get_5d_block_c
      MODULE PROCEDURE dbcsrt_get_6d_block_d
      MODULE PROCEDURE dbcsrt_get_6d_block_s
      MODULE PROCEDURE dbcsrt_get_6d_block_z
      MODULE PROCEDURE dbcsrt_get_6d_block_c
      MODULE PROCEDURE dbcsrt_get_anyd_block
   END INTERFACE

   INTERFACE dbcsrt_reserve_blocks
      MODULE PROCEDURE dbcsrt_reserve_blocks_index
      MODULE PROCEDURE dbcsrt_reserve_blocks_template
   END INTERFACE

   INTERFACE dbcsrt_ndims
      MODULE PROCEDURE ndims_tensor
      MODULE PROCEDURE ndims_iterator
      MODULE PROCEDURE ndims_mapping
      MODULE PROCEDURE ndims_buffer
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Number of non-zero blocks
!> \param tensor ...
!> \retval nblocks_reserved ...
! **************************************************************************************************
   FUNCTION nblocks_reserved(tensor)
      TYPE(dbcsrt_type), INTENT(IN)                      :: tensor
      INTEGER                                            :: nblocks_reserved

      INTEGER                                            :: blk, blk_p
      TYPE(dbcsrt_iterator)                              :: iterator
      INTEGER, DIMENSION(dbcsrt_ndims(tensor))           :: blk_offset, blk_size, ind_nd

      CPASSERT(tensor%valid)

      CALL dbcsrt_iterator_start(iterator, tensor)
      nblocks_reserved = 0
      DO WHILE (dbcsrt_iterator_blocks_left(iterator))
         CALL dbcsrt_iterator_next_block(iterator, ind_nd, blk, blk_p, blk_size, blk_offset)
         nblocks_reserved = nblocks_reserved+1
      ENDDO
      CALL dbcsrt_iterator_stop(iterator)
   END FUNCTION

! **************************************************************************************************
!> \brief indices of non-zero blocks
!> \param tensor ...
!> \param blk_ind_1 ...
!> \param blk_ind_2 ...
!> \param blk_ind_3 ...
!> \param blk_ind_4 ...
!> \param blk_ind_5 ...
!> \param blk_ind_6 ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_reserved_block_indices(tensor, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, &
                                            blk_ind_5, blk_ind_6)
      TYPE(dbcsrt_type), INTENT(IN)                      :: tensor
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: blk_ind_1, blk_ind_2, blk_ind_3, &
                                                            blk_ind_4, blk_ind_5, blk_ind_6

      INTEGER                                            :: blk, iblk, nblk
      TYPE(dbcsrt_iterator)                              :: iterator
      INTEGER, DIMENSION(dbcsrt_ndims(tensor))           :: ind_nd

      CPASSERT(tensor%valid)

      nblk = nblocks_reserved(tensor)

      IF (dbcsrt_ndims(tensor) >= 3) ALLOCATE (blk_ind_1(nblk), blk_ind_2(nblk), blk_ind_3(nblk))
      IF (dbcsrt_ndims(tensor) >= 4) ALLOCATE (blk_ind_4(nblk))
      IF (dbcsrt_ndims(tensor) >= 5) ALLOCATE (blk_ind_5(nblk))
      IF (dbcsrt_ndims(tensor) >= 6) ALLOCATE (blk_ind_6(nblk))

      CALL dbcsrt_iterator_start(iterator, tensor)
      DO iblk = 1, nblk
         CALL dbcsrt_iterator_next_block(iterator, ind_nd, blk)
         IF (dbcsrt_ndims(tensor) >= 3) THEN
            blk_ind_1(iblk) = ind_nd(1)
            blk_ind_2(iblk) = ind_nd(2)
            blk_ind_3(iblk) = ind_nd(3)
         ENDIF
         IF (dbcsrt_ndims(tensor) >= 4) blk_ind_4(iblk) = ind_nd(4)
         IF (dbcsrt_ndims(tensor) >= 5) blk_ind_5(iblk) = ind_nd(5)
         IF (dbcsrt_ndims(tensor) >= 6) blk_ind_6(iblk) = ind_nd(6)
      ENDDO
      CALL dbcsrt_iterator_stop(iterator)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_iterator_start for tensors.
!> \param iterator ...
!> \param tensor ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_iterator_start(iterator, tensor)
      TYPE(dbcsrt_iterator), INTENT(OUT)                 :: iterator
      TYPE(dbcsrt_type), INTENT(IN)                      :: tensor

      CPASSERT(tensor%valid)

      CALL dbcsr_iterator_start(iterator%iter, tensor%matrix_rep)
      iterator%nd_index_blk = tensor%nd_index_blk
      iterator%nd_index = tensor%nd_index
      iterator%blk_sizes = tensor%blk_sizes
      iterator%blk_offsets = tensor%blk_offsets
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_iterator_stop for tensors.
!> \param iterator ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_iterator_stop(iterator)
      TYPE(dbcsrt_iterator), INTENT(INOUT)               :: iterator

      CALL dbcsr_iterator_stop(iterator%iter)
      CALL destroy_nd_to_2d_mapping(iterator%nd_index)
      CALL destroy_nd_to_2d_mapping(iterator%nd_index_blk)
      CALL destroy_array_list(iterator%blk_sizes)
      CALL destroy_array_list(iterator%blk_offsets)
   END SUBROUTINE

! **************************************************************************************************
!> \brief iterate over nd blocks of an nd rank tensor
!>        index only
!> \param iterator ...
!> \param ind_nd nd index of block
!> \param blk is this needed?
!> \param blk_p is this needed?
!> \param blk_size blk size in each dimension
!> \param blk_offset blk offset in each dimension
! **************************************************************************************************
   SUBROUTINE dbcsrt_iterator_next_block(iterator, ind_nd, blk, blk_p, blk_size, blk_offset)
      TYPE(dbcsrt_iterator), INTENT(INOUT)               :: iterator
      INTEGER, DIMENSION(dbcsrt_ndims(iterator)), &
         INTENT(OUT)                                     :: ind_nd
      INTEGER, INTENT(OUT)                               :: blk
      INTEGER, INTENT(OUT), OPTIONAL                     :: blk_p
      INTEGER, DIMENSION(dbcsrt_ndims(iterator)), &
         INTENT(OUT), OPTIONAL                           :: blk_size, blk_offset

      INTEGER, DIMENSION(2)                              :: ind_2d

      CALL dbcsr_iterator_next_block(iterator%iter, ind_2d(1), ind_2d(2), blk, blk_p=blk_p)

      ind_nd(:) = get_nd_indices(iterator%nd_index_blk, ind_2d)
      IF (PRESENT(blk_size)) blk_size(:) = get_array_elements(iterator%blk_sizes, ind_nd)
      ! note: blk_offset needs to be determined by tensor metadata, can not be derived from 2d row/col
      ! offset since block index mapping is not consistent with element index mapping
      IF (PRESENT(blk_offset)) blk_offset(:) = get_array_elements(iterator%blk_offsets, ind_nd)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_iterator_blocks_left for tensors.
!> \param iterator ...
!> \retval dbcsrt_iterator_blocks_left ...
! **************************************************************************************************
   FUNCTION dbcsrt_iterator_blocks_left(iterator)
      TYPE(dbcsrt_iterator), INTENT(IN)                  :: iterator
      LOGICAL                                            :: dbcsrt_iterator_blocks_left

      dbcsrt_iterator_blocks_left = dbcsr_iterator_blocks_left(iterator%iter)
   END FUNCTION

! **************************************************************************************************
!> \brief Number of dimensions.
!> \param iterator ...
!> \retval ndims_iterator ...
! **************************************************************************************************
   PURE FUNCTION ndims_iterator(iterator)
      TYPE(dbcsrt_iterator), INTENT(IN)                  :: iterator
      INTEGER                                            :: ndims_iterator

      ndims_iterator = iterator%nd_index%ndim_nd
   END FUNCTION

! **************************************************************************************************
!> \brief reserve tensor blocks using block indices
!> \param tensor ...
!> \param blk_ind_1 index of blocks to reserve in each dimension
!> \param blk_ind_2 ...
!> \param blk_ind_3 ...
!> \param blk_ind_4 ...
!> \param blk_ind_5 ...
!> \param blk_ind_6 ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_reserve_blocks_index(tensor, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, blk_ind_5, blk_ind_6)
      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN), &
         OPTIONAL                                        :: blk_ind_1, blk_ind_2, blk_ind_3, &
                                                            blk_ind_4, blk_ind_5, blk_ind_6

      INTEGER                                            :: iblk, nblk
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cols, rows
      INTEGER, DIMENSION(2)                              :: ind_2d
      TYPE(array_list)                                   :: blks
      INTEGER, DIMENSION(dbcsrt_ndims(tensor))           :: iblk_nd, ind_nd, nblk_tmp

      CPASSERT(tensor%valid)

      CALL create_array_list(blks, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, blk_ind_5, blk_ind_6, &
                             dbcsrt_ndims(tensor))
      nblk_tmp(:) = sizes_of_arrays(blks)
      nblk = nblk_tmp(1)
      ALLOCATE (cols(nblk), rows(nblk))
      DO iblk = 1, nblk
         iblk_nd(:) = iblk
         ind_nd(:) = get_array_elements(blks, iblk_nd)
         ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind_nd)
         rows(iblk) = ind_2d(1); cols(iblk) = ind_2d(2)
      ENDDO

      CALL dbcsr_reserve_blocks(tensor%matrix_rep, rows=rows, columns=cols)
   END SUBROUTINE

! **************************************************************************************************
!> \brief reserve tensor blocks using template
!> \param tensor_in template tensor
!> \param tensor_out ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_reserve_blocks_template(tensor_in, tensor_out)
      TYPE(dbcsrt_type), INTENT(IN)                      :: tensor_in
      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor_out

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_ind_1, blk_ind_2, blk_ind_3, &
                                                            blk_ind_4, blk_ind_5, blk_ind_6

      CALL dbcsrt_reserved_block_indices(tensor_in, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, blk_ind_5, blk_ind_6)
      CALL dbcsrt_reserve_blocks(tensor_out, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, blk_ind_5, blk_ind_6)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_get_stored_coordinates for tensors.
!> \param tensor ...
!> \param ind_nd ...
!> \param processor ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_get_stored_coordinates(tensor, ind_nd, processor)
      TYPE(dbcsrt_type), INTENT(IN)                      :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), &
         INTENT(INOUT)                                   :: ind_nd
      INTEGER, INTENT(OUT), OPTIONAL                     :: processor

      INTEGER, DIMENSION(2)                              :: ind_2d

      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind_nd)
      CALL dbcsr_get_stored_coordinates(tensor%matrix_rep, ind_2d(1), ind_2d(2), processor)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Copy tensor data.
!>        Redistributes tensor data according to distributions of target and source tensor.
!>        Permutes tensor index according to `order` argument (if present).
!>
!> Source and target tensor formats are arbitrary as long as the following requirements are met:
!> * source and target tensors have the same number of blocks in each dimension, the same block sizes
!>   and the same block sparsity. If `order` argument is present, this must be the case after index
!>   permutation.
!> OR
!> * target tensor is not yet created, in this case an exact copy of source tensor is returned.
!> \param tensor_in Source
!> \param tensor_out Target
!> \param order Permutation of target tensor index. Exact same convention as order argument of RESHAPE
!>        intrinsic
! **************************************************************************************************
   SUBROUTINE dbcsrt_copy(tensor_in, tensor_out, order)
      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor_in, tensor_out
      INTEGER, DIMENSION(dbcsrt_ndims(tensor_in)), &
         INTENT(IN), OPTIONAL                            :: order

      TYPE(dbcsrt_type)                                  :: tensor_tmp

      IF (PRESENT(order)) THEN
         CALL dbcsrt_permute_index(tensor_in, tensor_tmp, order)
         CALL dbcsrt_copy_low(tensor_tmp, tensor_out)
      ELSE
         CALL dbcsrt_copy_low(tensor_in, tensor_out)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create new tensor by reordering index, data is copied exactly
!> \param tensor_in Source
!> \param tensor_out Target
!> \param order ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_permute_index(tensor_in, tensor_out, order)
      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor_in, tensor_out
      INTEGER, DIMENSION(dbcsrt_ndims(tensor_in)), &
         INTENT(IN)                                      :: order

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: data1, data2, data3, data4, data5, data6
      TYPE(nd_to_2d_mapping)                             :: nd_index_blk_rs, nd_index_rs

      CALL permute_index(tensor_in%nd_index, nd_index_rs, order)
      CALL permute_index(tensor_in%nd_index_blk, nd_index_blk_rs, order)

      CALL dbcsr_create(tensor_out%matrix_rep, tensor_in%matrix_rep)

      CALL dbcsr_copy(tensor_out%matrix_rep, tensor_in%matrix_rep)
      tensor_out%nd_index = nd_index_rs
      tensor_out%nd_index_blk = nd_index_blk_rs

      CALL get_arrays(tensor_in%blk_sizes, data1, data2, data3, data4, data5, data6, order)
      CALL create_array_list(tensor_out%blk_sizes, data1, data2, data3, data4, data5, data6, dbcsrt_ndims(tensor_in))
   END SUBROUTINE

! **************************************************************************************************
!> \brief reorder tensor index (no data)
!> \param map_in ...
!> \param map_out ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE permute_index(map_in, map_out, order)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map_in
      TYPE(nd_to_2d_mapping), INTENT(OUT)                :: map_out
      INTEGER, DIMENSION(dbcsrt_ndims(map_in)), &
         INTENT(IN)                                      :: order

      INTEGER                                            :: idim_nd, ndim_nd
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: map1_2d, map1_2d_reorder, map2_2d, &
                                                            map2_2d_reorder
      INTEGER, DIMENSION(dbcsrt_ndims(map_in))           :: dims_nd, dims_reorder

      CALL get_mapping_info(map_in, ndim_nd, dims_nd=dims_nd, map1_2d=map1_2d, map2_2d=map2_2d)

      CALL allocate_any(SHAPE(map1_2d), map1_2d_reorder)
      CALL allocate_any(SHAPE(map2_2d), map2_2d_reorder)
      DO idim_nd = 1, ndim_nd
         dims_reorder(order(idim_nd)) = dims_nd(idim_nd)
         map1_2d_reorder(order(idim_nd)) = map1_2d(idim_nd)
         map2_2d_reorder(order(idim_nd)) = map2_2d(idim_nd)
      ENDDO

      CALL create_nd_to_2d_mapping(map_out, dims_reorder, map1_2d_reorder, map2_2d_reorder)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generic implementation of dbcsrt_get_block (arbitrary tensor rank and arbitrary datatype)
!> \param tensor ...
!> \param ind block index
!> \param block block to get
! **************************************************************************************************
   SUBROUTINE dbcsrt_get_anyd_block(tensor, ind, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), &
         INTENT(IN)                                      :: ind
      TYPE(block_nd), INTENT(OUT)                        :: block

      SELECT CASE (dbcsrt_get_data_type (tensor))
      CASE (dbcsr_type_real_8)
         CALL dbcsrt_get_anyd_block_d(tensor, ind, block)
      CASE (dbcsr_type_real_4)
         CALL dbcsrt_get_anyd_block_s(tensor, ind, block)
      CASE (dbcsr_type_complex_8)
         CALL dbcsrt_get_anyd_block_z(tensor, ind, block)
      CASE (dbcsr_type_complex_4)
         CALL dbcsrt_get_anyd_block_c(tensor, ind, block)
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generic implementation of dbcsrt_put_block (arbitrary tensor rank and arbitrary datatype)
!> \param tensor ...
!> \param ind block index
!> \param block block to put
! **************************************************************************************************
   SUBROUTINE dbcsrt_put_anyd_block(tensor, ind, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), &
         INTENT(IN)                                      :: ind
      TYPE(block_nd), INTENT(IN)                         :: block

      SELECT CASE (block%data_type)
      CASE (dbcsr_type_real_8)
         CALL dbcsrt_put_anyd_block_d(tensor, ind, block%r_dp)
      CASE (dbcsr_type_real_4)
         CALL dbcsrt_put_anyd_block_s(tensor, ind, block%r_sp)
      CASE (dbcsr_type_complex_8)
         CALL dbcsrt_put_anyd_block_z(tensor, ind, block%c_dp)
      CASE (dbcsr_type_complex_4)
         CALL dbcsrt_put_anyd_block_c(tensor, ind, block%c_sp)
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief Helper of dbcsrt_copy: copy data after index permutation has been done.
!> \param tensor_in ...
!> \param tensor_out ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_copy_low(tensor_in, tensor_out)

      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor_in, tensor_out

      INTEGER                                            :: blk, iproc, mp_comm, mynode, ndata, &
                                                            numnodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: req_array
      TYPE(dbcsrt_iterator)                              :: iter
      TYPE(block_nd)                                     :: blk_data
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      INTEGER, DIMENSION(dbcsrt_ndims(tensor_in))        :: blk_size, ind_nd, index

      IF (tensor_out%valid) THEN

         mp_comm = dbcsr_mp_group(dbcsr_distribution_mp(dbcsr_distribution(tensor_in%matrix_rep)))
         CALL mp_environ(numnodes, mynode, mp_comm)
         ALLOCATE (buffer_send(0:numnodes-1))
         ALLOCATE (buffer_recv(0:numnodes-1))
         ALLOCATE (num_blocks_recv(0:numnodes-1))
         ALLOCATE (num_blocks_send(0:numnodes-1))
         ALLOCATE (num_entries_recv(0:numnodes-1))
         ALLOCATE (num_entries_send(0:numnodes-1))
         ALLOCATE (num_rec(0:2*numnodes-1))
         ALLOCATE (num_send(0:2*numnodes-1))
         num_send(:) = 0
         ALLOCATE (req_array(1:numnodes, 4))
         CALL dbcsrt_iterator_start(iter, tensor_in)
         DO WHILE (dbcsrt_iterator_blocks_left(iter))
            CALL dbcsrt_iterator_next_block(iter, ind_nd, blk, blk_size=blk_size)
            CALL dbcsrt_get_stored_coordinates(tensor_out, ind_nd, iproc)
            num_send(2*iproc) = num_send(2*iproc)+PRODUCT(blk_size)
            num_send(2*iproc+1) = num_send(2*iproc+1)+1
         ENDDO
         CALL dbcsrt_iterator_stop(iter)
         CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
         DO iproc = 0, numnodes-1
            num_entries_recv(iproc) = num_rec(2*iproc)
            num_blocks_recv(iproc) = num_rec(2*iproc+1)
            num_entries_send(iproc) = num_send(2*iproc)
            num_blocks_send(iproc) = num_send(2*iproc+1)

            CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc), &
                                     dbcsrt_get_data_type(tensor_in), dbcsrt_ndims(tensor_in))
            CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc), &
                                     dbcsrt_get_data_type(tensor_in), dbcsrt_ndims(tensor_in))
         ENDDO
         CALL dbcsrt_iterator_start(iter, tensor_in)
         DO WHILE (dbcsrt_iterator_blocks_left(iter))
            CALL dbcsrt_iterator_next_block(iter, ind_nd, blk, blk_size=blk_size)
            CALL dbcsrt_get_anyd_block(tensor_in, ind_nd, blk_data)
            CALL dbcsrt_get_stored_coordinates(tensor_out, ind_nd, iproc)
            CALL block_buffer_add_anyd_block(buffer_send(iproc), PRODUCT(blk_size), ind_nd, blk_data)
            CALL destroy_block(blk_data)
         ENDDO
         CALL dbcsrt_iterator_stop(iter)

         CALL communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)

         DO iproc = 0, numnodes-1
            ! First, we need to get the index to create block
            DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
               CALL block_buffer_get_next_anyd_block(buffer_recv(iproc), ndata, index)
               CALL dbcsrt_blk_sizes(tensor_in, index, blk_size)
               ! create block
               CALL create_block(blk_data, blk_size, dbcsrt_get_data_type(tensor_in))
               ! get actual block data
               CALL block_buffer_get_next_anyd_block(buffer_recv(iproc), ndata, index, blk_data)
               CALL dbcsrt_put_anyd_block(tensor_out, index, blk_data)
               CALL destroy_block(blk_data)
            ENDDO
         ENDDO
      ELSE
         CALL dbcsrt_create(tensor_in, tensor_out)
         CALL dbcsrt_reserve_blocks(tensor_in, tensor_out)
         CALL dbcsr_copy(tensor_in%matrix_rep, tensor_out%matrix_rep)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief Contract tensors by multiplying matrix representations. NYI.
!> \param tensor_1 first tensor
!> \param tensor_2 second tensor
!> \param contract_1 indices of tensor_1 to contract
!> \param contract_2 indices of tensor_2 to contract
!> \param tensor_prod ...
! **************************************************************************************************
   SUBROUTINE dbcsrt_contract(tensor_1, tensor_2, contract_1, contract_2, tensor_prod)
      TYPE(dbcsrt_type), INTENT(IN)                      :: tensor_1, tensor_2
      INTEGER, DIMENSION(:), INTENT(IN)                  :: contract_1, contract_2
      TYPE(dbcsrt_type), INTENT(INOUT)                   :: tensor_prod

      CPASSERT(tensor_1%valid)
      CPASSERT(tensor_2%valid)
      CPASSERT(SIZE(contract_1) .EQ. SIZE(contract_2))

      MARK_USED(tensor_1)
      MARK_USED(tensor_2)
      MARK_USED(tensor_prod)
      MARK_USED(contract_1)
      MARK_USED(contract_2)
      CPABORT("not yet implemented")
   END SUBROUTINE

#include "dbcsrt_inc_d.f90"
#include "dbcsrt_inc_s.f90"
#include "dbcsrt_inc_z.f90"
#include "dbcsrt_inc_c.f90"
END MODULE dbcsrt
