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
MODULE dbcsr_tensor
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
   USE dbcsr_tensor_block,              ONLY: &
        block_buffer_add_anyd_block, block_buffer_blocks_left, block_buffer_create, &
        block_buffer_get_next_anyd_block, block_buffer_type, block_nd, communicate_buffer, &
        create_block, destroy_block, ndims_buffer, reshape_2d_to_nd_block, reshape_nd_to_2d_block
   USE dbcsr_tensor_block_c,            ONLY: block_nd_c
   USE dbcsr_tensor_block_d,            ONLY: block_nd_d
   USE dbcsr_tensor_block_s,            ONLY: block_nd_s
   USE dbcsr_tensor_block_z,            ONLY: block_nd_z
   USE dbcsr_tensor_index,              ONLY: get_2d_indices,&
                                              get_nd_indices
   USE dbcsr_tensor_index_types,        ONLY: create_nd_to_2d_mapping,&
                                              destroy_nd_to_2d_mapping,&
                                              get_mapping_info,&
                                              nd_to_2d_mapping,&
                                              ndims_mapping
   USE dbcsr_tensor_types,              ONLY: dbcsr_t_blk_sizes,&
                                              dbcsr_t_create,&
                                              dbcsr_t_get_data_type,&
                                              dbcsr_t_type,&
                                              ndims_tensor
   USE kinds,                           ONLY: real_4,&
                                              real_8
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_environ
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor'

   PUBLIC :: &
      dbcsr_t_contract, &
      dbcsr_t_copy, &
      dbcsr_t_get_block, &
      dbcsr_t_get_stored_coordinates, &
      dbcsr_t_iterator, &
      dbcsr_t_iterator_blocks_left, &
      dbcsr_t_iterator_next_block, &
      dbcsr_t_iterator_start, &
      dbcsr_t_iterator_stop, &
      dbcsr_t_ndims, &
      dbcsr_t_permute_index, &
      dbcsr_t_put_block, &
      dbcsr_t_reserve_blocks

   TYPE dbcsr_t_iterator
      TYPE(dbcsr_iterator)             :: iter
      TYPE(nd_to_2d_mapping)           :: nd_index_blk
      TYPE(nd_to_2d_mapping)           :: nd_index
      TYPE(array_list)                 :: blk_sizes, blk_offsets
   END TYPE dbcsr_t_iterator

   INTERFACE dbcsr_t_put_block
      MODULE PROCEDURE dbcsr_t_put_3d_block_d
      MODULE PROCEDURE dbcsr_t_put_3d_block_s
      MODULE PROCEDURE dbcsr_t_put_3d_block_z
      MODULE PROCEDURE dbcsr_t_put_3d_block_c
      MODULE PROCEDURE dbcsr_t_put_4d_block_d
      MODULE PROCEDURE dbcsr_t_put_4d_block_s
      MODULE PROCEDURE dbcsr_t_put_4d_block_z
      MODULE PROCEDURE dbcsr_t_put_4d_block_c
      MODULE PROCEDURE dbcsr_t_put_5d_block_d
      MODULE PROCEDURE dbcsr_t_put_5d_block_s
      MODULE PROCEDURE dbcsr_t_put_5d_block_z
      MODULE PROCEDURE dbcsr_t_put_5d_block_c
      MODULE PROCEDURE dbcsr_t_put_6d_block_d
      MODULE PROCEDURE dbcsr_t_put_6d_block_s
      MODULE PROCEDURE dbcsr_t_put_6d_block_z
      MODULE PROCEDURE dbcsr_t_put_6d_block_c
      MODULE PROCEDURE dbcsr_t_put_anyd_block
   END INTERFACE

   INTERFACE dbcsr_t_get_block
      MODULE PROCEDURE dbcsr_t_get_3d_block_d
      MODULE PROCEDURE dbcsr_t_get_3d_block_s
      MODULE PROCEDURE dbcsr_t_get_3d_block_z
      MODULE PROCEDURE dbcsr_t_get_3d_block_c
      MODULE PROCEDURE dbcsr_t_get_4d_block_d
      MODULE PROCEDURE dbcsr_t_get_4d_block_s
      MODULE PROCEDURE dbcsr_t_get_4d_block_z
      MODULE PROCEDURE dbcsr_t_get_4d_block_c
      MODULE PROCEDURE dbcsr_t_get_5d_block_d
      MODULE PROCEDURE dbcsr_t_get_5d_block_s
      MODULE PROCEDURE dbcsr_t_get_5d_block_z
      MODULE PROCEDURE dbcsr_t_get_5d_block_c
      MODULE PROCEDURE dbcsr_t_get_6d_block_d
      MODULE PROCEDURE dbcsr_t_get_6d_block_s
      MODULE PROCEDURE dbcsr_t_get_6d_block_z
      MODULE PROCEDURE dbcsr_t_get_6d_block_c
      MODULE PROCEDURE dbcsr_t_get_anyd_block
   END INTERFACE

   INTERFACE dbcsr_t_reserve_blocks
      MODULE PROCEDURE dbcsr_t_reserve_blocks_index
      MODULE PROCEDURE dbcsr_t_reserve_blocks_template
   END INTERFACE

   INTERFACE dbcsr_t_ndims
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
      TYPE(dbcsr_t_type), INTENT(IN)                     :: tensor
      INTEGER                                            :: nblocks_reserved

      INTEGER                                            :: blk, blk_p
      TYPE(dbcsr_t_iterator)                             :: iterator
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor))          :: blk_offset, blk_size, ind_nd

      CPASSERT(tensor%valid)

      CALL dbcsr_t_iterator_start(iterator, tensor)
      nblocks_reserved = 0
      DO WHILE (dbcsr_t_iterator_blocks_left(iterator))
         CALL dbcsr_t_iterator_next_block(iterator, ind_nd, blk, blk_p, blk_size, blk_offset)
         nblocks_reserved = nblocks_reserved+1
      ENDDO
      CALL dbcsr_t_iterator_stop(iterator)
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
   SUBROUTINE dbcsr_t_reserved_block_indices(tensor, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, &
                                             blk_ind_5, blk_ind_6)
      TYPE(dbcsr_t_type), INTENT(IN)                     :: tensor
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: blk_ind_1, blk_ind_2, blk_ind_3, &
                                                            blk_ind_4, blk_ind_5, blk_ind_6

      INTEGER                                            :: blk, iblk, nblk
      TYPE(dbcsr_t_iterator)                             :: iterator
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor))          :: ind_nd

      CPASSERT(tensor%valid)

      nblk = nblocks_reserved(tensor)

      IF (dbcsr_t_ndims(tensor) >= 3) ALLOCATE (blk_ind_1(nblk), blk_ind_2(nblk), blk_ind_3(nblk))
      IF (dbcsr_t_ndims(tensor) >= 4) ALLOCATE (blk_ind_4(nblk))
      IF (dbcsr_t_ndims(tensor) >= 5) ALLOCATE (blk_ind_5(nblk))
      IF (dbcsr_t_ndims(tensor) >= 6) ALLOCATE (blk_ind_6(nblk))

      CALL dbcsr_t_iterator_start(iterator, tensor)
      DO iblk = 1, nblk
         CALL dbcsr_t_iterator_next_block(iterator, ind_nd, blk)
         IF (dbcsr_t_ndims(tensor) >= 3) THEN
            blk_ind_1(iblk) = ind_nd(1)
            blk_ind_2(iblk) = ind_nd(2)
            blk_ind_3(iblk) = ind_nd(3)
         ENDIF
         IF (dbcsr_t_ndims(tensor) >= 4) blk_ind_4(iblk) = ind_nd(4)
         IF (dbcsr_t_ndims(tensor) >= 5) blk_ind_5(iblk) = ind_nd(5)
         IF (dbcsr_t_ndims(tensor) >= 6) blk_ind_6(iblk) = ind_nd(6)
      ENDDO
      CALL dbcsr_t_iterator_stop(iterator)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_iterator_start for tensors.
!> \param iterator ...
!> \param tensor ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_iterator_start(iterator, tensor)
      TYPE(dbcsr_t_iterator), INTENT(OUT)                :: iterator
      TYPE(dbcsr_t_type), INTENT(IN)                     :: tensor

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
   SUBROUTINE dbcsr_t_iterator_stop(iterator)
      TYPE(dbcsr_t_iterator), INTENT(INOUT)              :: iterator

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
   SUBROUTINE dbcsr_t_iterator_next_block(iterator, ind_nd, blk, blk_p, blk_size, blk_offset)
      TYPE(dbcsr_t_iterator), INTENT(INOUT)              :: iterator
      INTEGER, DIMENSION(iterator%nd_index%ndim_nd), &
         INTENT(OUT)                                     :: ind_nd
      INTEGER, INTENT(OUT)                               :: blk
      INTEGER, INTENT(OUT), OPTIONAL                     :: blk_p
      INTEGER, DIMENSION(iterator%nd_index%ndim_nd), &
         INTENT(OUT), OPTIONAL                           :: blk_size, blk_offset

      INTEGER, DIMENSION(2)                              :: ind_2d

!     Note: specification function below not accepted by IBM XL compiler (compilation fails).
!     TODO: Uncomment once this has been fixed
!     INTEGER, DIMENSION(dbcsr_t_ndims(iterator)), &
!        INTENT(OUT)                                     :: ind_nd
!     INTEGER, DIMENSION(dbcsr_t_ndims(iterator)), &
!        INTENT(OUT), OPTIONAL                           :: blk_size, blk_offset

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
!> \retval dbcsr_t_iterator_blocks_left ...
! **************************************************************************************************
   FUNCTION dbcsr_t_iterator_blocks_left(iterator)
      TYPE(dbcsr_t_iterator), INTENT(IN)                 :: iterator
      LOGICAL                                            :: dbcsr_t_iterator_blocks_left

      dbcsr_t_iterator_blocks_left = dbcsr_iterator_blocks_left(iterator%iter)
   END FUNCTION

! **************************************************************************************************
!> \brief Number of dimensions.
!> \param iterator ...
!> \retval ndims_iterator ...
! **************************************************************************************************
   PURE FUNCTION ndims_iterator(iterator)
      TYPE(dbcsr_t_iterator), INTENT(IN)                 :: iterator
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
   SUBROUTINE dbcsr_t_reserve_blocks_index(tensor, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, blk_ind_5, blk_ind_6)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN), &
         OPTIONAL                                        :: blk_ind_1, blk_ind_2, blk_ind_3, &
                                                            blk_ind_4, blk_ind_5, blk_ind_6

      INTEGER                                            :: iblk, nblk
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cols, rows
      INTEGER, DIMENSION(2)                              :: ind_2d
      TYPE(array_list)                                   :: blks
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor))          :: iblk_nd, ind_nd, nblk_tmp

      CPASSERT(tensor%valid)

      CALL create_array_list(blks, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, blk_ind_5, blk_ind_6, &
                             dbcsr_t_ndims(tensor))
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
   SUBROUTINE dbcsr_t_reserve_blocks_template(tensor_in, tensor_out)
      TYPE(dbcsr_t_type), INTENT(IN)                     :: tensor_in
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor_out

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_ind_1, blk_ind_2, blk_ind_3, &
                                                            blk_ind_4, blk_ind_5, blk_ind_6

      CALL dbcsr_t_reserved_block_indices(tensor_in, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, blk_ind_5, blk_ind_6)
      CALL dbcsr_t_reserve_blocks(tensor_out, blk_ind_1, blk_ind_2, blk_ind_3, blk_ind_4, blk_ind_5, blk_ind_6)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of dbcsr_get_stored_coordinates for tensors.
!> \param tensor ...
!> \param ind_nd ...
!> \param processor ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_get_stored_coordinates(tensor, ind_nd, processor)
      TYPE(dbcsr_t_type), INTENT(IN)                     :: tensor
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor)), &
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
   SUBROUTINE dbcsr_t_copy(tensor_in, tensor_out, order)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor_in, tensor_out
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)), &
         INTENT(IN), OPTIONAL                            :: order

      TYPE(dbcsr_t_type)                                 :: tensor_tmp

      IF (PRESENT(order)) THEN
         CALL dbcsr_t_permute_index(tensor_in, tensor_tmp, order)
         CALL dbcsr_t_copy_low(tensor_tmp, tensor_out)
      ELSE
         CALL dbcsr_t_copy_low(tensor_in, tensor_out)
      ENDIF
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create new tensor by reordering index, data is copied exactly
!> \param tensor_in Source
!> \param tensor_out Target
!> \param order ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_permute_index(tensor_in, tensor_out, order)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor_in, tensor_out
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)), &
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
      CALL create_array_list(tensor_out%blk_sizes, data1, data2, data3, data4, data5, data6, dbcsr_t_ndims(tensor_in))
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
      INTEGER, DIMENSION(dbcsr_t_ndims(map_in)), &
         INTENT(IN)                                      :: order

      INTEGER                                            :: idim_nd, ndim_nd
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: map1_2d, map1_2d_reorder, map2_2d, &
                                                            map2_2d_reorder
      INTEGER, DIMENSION(dbcsr_t_ndims(map_in))          :: dims_nd, dims_reorder

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
!> \brief Generic implementation of dbcsr_t_get_block (arbitrary tensor rank and arbitrary datatype)
!> \param tensor ...
!> \param ind block index
!> \param block block to get
! **************************************************************************************************
   SUBROUTINE dbcsr_t_get_anyd_block(tensor, ind, block)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor)), &
         INTENT(IN)                                      :: ind
      TYPE(block_nd), INTENT(OUT)                        :: block

      SELECT CASE (dbcsr_t_get_data_type (tensor))
      CASE (dbcsr_type_real_8)
         CALL dbcsr_t_get_anyd_block_d(tensor, ind, block)
      CASE (dbcsr_type_real_4)
         CALL dbcsr_t_get_anyd_block_s(tensor, ind, block)
      CASE (dbcsr_type_complex_8)
         CALL dbcsr_t_get_anyd_block_z(tensor, ind, block)
      CASE (dbcsr_type_complex_4)
         CALL dbcsr_t_get_anyd_block_c(tensor, ind, block)
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generic implementation of dbcsr_t_put_block (arbitrary tensor rank and arbitrary datatype)
!> \param tensor ...
!> \param ind block index
!> \param block block to put
! **************************************************************************************************
   SUBROUTINE dbcsr_t_put_anyd_block(tensor, ind, block)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor)), &
         INTENT(IN)                                      :: ind
      TYPE(block_nd), INTENT(IN)                         :: block

      SELECT CASE (block%data_type)
      CASE (dbcsr_type_real_8)
         CALL dbcsr_t_put_anyd_block_d(tensor, ind, block%r_dp)
      CASE (dbcsr_type_real_4)
         CALL dbcsr_t_put_anyd_block_s(tensor, ind, block%r_sp)
      CASE (dbcsr_type_complex_8)
         CALL dbcsr_t_put_anyd_block_z(tensor, ind, block%c_dp)
      CASE (dbcsr_type_complex_4)
         CALL dbcsr_t_put_anyd_block_c(tensor, ind, block%c_sp)
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief Helper of dbcsr_t_copy: copy data after index permutation has been done.
!> \param tensor_in ...
!> \param tensor_out ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_copy_low(tensor_in, tensor_out)

      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor_in, tensor_out

      INTEGER                                            :: blk, iproc, mp_comm, mynode, ndata, &
                                                            numnodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: req_array
      TYPE(dbcsr_t_iterator)                             :: iter
      TYPE(block_nd)                                     :: blk_data
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in))       :: blk_size, ind_nd, index

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
         CALL dbcsr_t_iterator_start(iter, tensor_in)
         DO WHILE (dbcsr_t_iterator_blocks_left(iter))
            CALL dbcsr_t_iterator_next_block(iter, ind_nd, blk, blk_size=blk_size)
            CALL dbcsr_t_get_stored_coordinates(tensor_out, ind_nd, iproc)
            num_send(2*iproc) = num_send(2*iproc)+PRODUCT(blk_size)
            num_send(2*iproc+1) = num_send(2*iproc+1)+1
         ENDDO
         CALL dbcsr_t_iterator_stop(iter)
         CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
         DO iproc = 0, numnodes-1
            num_entries_recv(iproc) = num_rec(2*iproc)
            num_blocks_recv(iproc) = num_rec(2*iproc+1)
            num_entries_send(iproc) = num_send(2*iproc)
            num_blocks_send(iproc) = num_send(2*iproc+1)

            CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc), &
                                     dbcsr_t_get_data_type(tensor_in), dbcsr_t_ndims(tensor_in))
            CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc), &
                                     dbcsr_t_get_data_type(tensor_in), dbcsr_t_ndims(tensor_in))
         ENDDO
         CALL dbcsr_t_iterator_start(iter, tensor_in)
         DO WHILE (dbcsr_t_iterator_blocks_left(iter))
            CALL dbcsr_t_iterator_next_block(iter, ind_nd, blk, blk_size=blk_size)
            CALL dbcsr_t_get_anyd_block(tensor_in, ind_nd, blk_data)
            CALL dbcsr_t_get_stored_coordinates(tensor_out, ind_nd, iproc)
            CALL block_buffer_add_anyd_block(buffer_send(iproc), ind_nd, blk_data)
            CALL destroy_block(blk_data)
         ENDDO
         CALL dbcsr_t_iterator_stop(iter)

         CALL communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)

         DO iproc = 0, numnodes-1
            ! First, we need to get the index to create block
            DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
               CALL block_buffer_get_next_anyd_block(buffer_recv(iproc), ndata, index)
               CALL dbcsr_t_blk_sizes(tensor_in, index, blk_size)
               ! create block
               CALL create_block(blk_data, blk_size, dbcsr_t_get_data_type(tensor_in))
               ! get actual block data
               CALL block_buffer_get_next_anyd_block(buffer_recv(iproc), ndata, index, blk_data)
               CALL dbcsr_t_put_anyd_block(tensor_out, index, blk_data)
               CALL destroy_block(blk_data)
            ENDDO
         ENDDO
      ELSE
         CALL dbcsr_t_create(tensor_in, tensor_out)
         CALL dbcsr_t_reserve_blocks(tensor_in, tensor_out)
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
   SUBROUTINE dbcsr_t_contract(tensor_1, tensor_2, contract_1, contract_2, tensor_prod)
      TYPE(dbcsr_t_type), INTENT(IN)                     :: tensor_1, tensor_2
      INTEGER, DIMENSION(:), INTENT(IN)                  :: contract_1, contract_2
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor_prod

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

#include "dbcsr_tensor_inc_d.f90"
#include "dbcsr_tensor_inc_s.f90"
#include "dbcsr_tensor_inc_z.f90"
#include "dbcsr_tensor_inc_c.f90"
END MODULE
