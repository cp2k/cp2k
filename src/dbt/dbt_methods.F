!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBT tensor framework for block-sparse tensor contraction.
!>        Representation of n-rank tensors as DBT tall-and-skinny matrices.
!>        Support for arbitrary redistribution between different representations.
!>        Support for arbitrary tensor contractions
!> \todo implement checks and error messages
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_methods
   #:include "dbt_macros.fypp"
   #:set maxdim = maxrank
   #:set ndims = range(2,maxdim+1)

   USE dbcsr_api, ONLY: &
      dbcsr_type, dbcsr_release, &
      dbcsr_iterator_type, dbcsr_iterator_start, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
      dbcsr_has_symmetry, dbcsr_desymmetrize, dbcsr_put_block, dbcsr_clear, dbcsr_iterator_stop
   USE dbt_allocate_wrap, ONLY: &
      allocate_any
   USE dbt_array_list_methods, ONLY: &
      get_arrays, reorder_arrays, get_ith_array, array_list, array_sublist, check_equal, array_eq_i, &
      create_array_list, destroy_array_list, sizes_of_arrays
   USE dbm_api, ONLY: &
      dbm_clear
   USE dbt_tas_types, ONLY: &
      dbt_tas_split_info
   USE dbt_tas_base, ONLY: &
      dbt_tas_copy, dbt_tas_finalize, dbt_tas_get_info, dbt_tas_info
   USE dbt_tas_mm, ONLY: &
      dbt_tas_multiply, dbt_tas_batched_mm_init, dbt_tas_batched_mm_finalize, dbt_tas_result_index, &
      dbt_tas_batched_mm_complete, dbt_tas_set_batched_state
   USE dbt_block, ONLY: &
      dbt_iterator_type, dbt_get_block, dbt_put_block, dbt_iterator_start, &
      dbt_iterator_blocks_left, dbt_iterator_stop, dbt_iterator_next_block, &
      ndims_iterator, dbt_reserve_blocks, block_nd, destroy_block, checker_tr
   USE dbt_index, ONLY: &
      dbt_get_mapping_info, nd_to_2d_mapping, dbt_inverse_order, permute_index, get_nd_indices_tensor, &
      ndims_mapping_row, ndims_mapping_column, ndims_mapping
   USE dbt_types, ONLY: &
      dbt_create, dbt_type, ndims_tensor, dims_tensor, &
      dbt_distribution_type, dbt_distribution, dbt_nd_mp_comm, dbt_destroy, &
      dbt_distribution_destroy, dbt_distribution_new_expert, dbt_get_stored_coordinates, &
      blk_dims_tensor, dbt_hold, dbt_pgrid_type, mp_environ_pgrid, dbt_filter, &
      dbt_clear, dbt_finalize, dbt_get_num_blocks, dbt_scale, &
      dbt_get_num_blocks_total, dbt_get_info, ndims_matrix_row, ndims_matrix_column, &
      dbt_max_nblks_local, dbt_default_distvec, dbt_contraction_storage, dbt_nblks_total, &
      dbt_distribution_new, dbt_copy_contraction_storage, dbt_pgrid_destroy
   USE kinds, ONLY: &
      dp, default_string_length, int_8, dp
   USE message_passing, ONLY: &
      mp_environ, mp_max, mp_sum, mp_comm_free, mp_cart_create, mp_sync
   USE util, ONLY: &
      sort
   USE dbt_reshape_ops, ONLY: &
      dbt_reshape
   USE dbt_tas_split, ONLY: &
      dbt_tas_mp_comm, rowsplit, colsplit, dbt_tas_info_hold, dbt_tas_release_info, default_nsplit_accept_ratio, &
      default_pdims_accept_ratio, dbt_tas_create_split
   USE dbt_split, ONLY: &
      dbt_split_copyback, dbt_make_compatible_blocks, dbt_crop
   USE dbt_io, ONLY: &
      dbt_write_tensor_info, dbt_write_tensor_dist, prep_output_unit, dbt_write_split_info

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_methods'

   PUBLIC :: &
      dbt_contract, &
      dbt_copy, &
      dbt_get_block, &
      dbt_get_stored_coordinates, &
      dbt_inverse_order, &
      dbt_iterator_blocks_left, &
      dbt_iterator_next_block, &
      dbt_iterator_start, &
      dbt_iterator_stop, &
      dbt_iterator_type, &
      dbt_put_block, &
      dbt_reserve_blocks, &
      dbt_copy_matrix_to_tensor, &
      dbt_copy_tensor_to_matrix, &
      dbt_contract_index, &
      dbt_batched_contract_init, &
      dbt_batched_contract_finalize

CONTAINS

! **************************************************************************************************
!> \brief Copy tensor data.
!>        Redistributes tensor data according to distributions of target and source tensor.
!>        Permutes tensor index according to `order` argument (if present).
!>        Source and target tensor formats are arbitrary as long as the following requirements are met:
!>        * source and target tensors have the same rank and the same sizes in each dimension in terms
!>          of tensor elements (block sizes don't need to be the same).
!>          If `order` argument is present, sizes must match after index permutation.
!>        OR
!>        * target tensor is not yet created, in this case an exact copy of source tensor is returned.
!> \param tensor_in Source
!> \param tensor_out Target
!> \param order Permutation of target tensor index.
!>              Exact same convention as order argument of RESHAPE intrinsic.
!> \param bounds crop tensor data: start and end index for each tensor dimension
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_copy(tensor_in, tensor_out, order, summation, bounds, move_data, unit_nr)
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_in, tensor_out
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(IN), OPTIONAL                        :: order
      LOGICAL, INTENT(IN), OPTIONAL                  :: summation, move_data
      INTEGER, DIMENSION(2, ndims_tensor(tensor_in)), &
         INTENT(IN), OPTIONAL                        :: bounds
      INTEGER, INTENT(IN), OPTIONAL                  :: unit_nr
      INTEGER :: handle

      CALL mp_sync(tensor_in%pgrid%mp_comm_2d)
      CALL timeset("dbt_total", handle)

      ! make sure that it is safe to use dbt_copy during a batched contraction
      CALL dbt_tas_batched_mm_complete(tensor_in%matrix_rep, warn=.TRUE.)
      CALL dbt_tas_batched_mm_complete(tensor_out%matrix_rep, warn=.TRUE.)

      CALL dbt_copy_expert(tensor_in, tensor_out, order, summation, bounds, move_data, unit_nr)
      CALL mp_sync(tensor_in%pgrid%mp_comm_2d)
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief expert routine for copying a tensor. For internal use only.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_copy_expert(tensor_in, tensor_out, order, summation, bounds, move_data, unit_nr)
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_in, tensor_out
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(IN), OPTIONAL                        :: order
      LOGICAL, INTENT(IN), OPTIONAL                  :: summation, move_data
      INTEGER, DIMENSION(2, ndims_tensor(tensor_in)), &
         INTENT(IN), OPTIONAL                        :: bounds
      INTEGER, INTENT(IN), OPTIONAL                  :: unit_nr

      TYPE(dbt_type), POINTER                    :: in_tmp_1, in_tmp_2, &
                                                    in_tmp_3, out_tmp_1
      INTEGER                                        :: handle, unit_nr_prv
      INTEGER, DIMENSION(:), ALLOCATABLE             :: map1_in_1, map1_in_2, map2_in_1, map2_in_2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_copy'
      LOGICAL                                        :: dist_compatible_tas, dist_compatible_tensor, &
                                                        summation_prv, new_in_1, new_in_2, &
                                                        new_in_3, new_out_1, block_compatible, &
                                                        move_prv
      TYPE(array_list)                               :: blk_sizes_in

      CALL timeset(routineN, handle)

      CPASSERT(tensor_out%valid)

      unit_nr_prv = prep_output_unit(unit_nr)

      IF (PRESENT(move_data)) THEN
         move_prv = move_data
      ELSE
         move_prv = .FALSE.
      END IF

      dist_compatible_tas = .FALSE.
      dist_compatible_tensor = .FALSE.
      block_compatible = .FALSE.
      new_in_1 = .FALSE.
      new_in_2 = .FALSE.
      new_in_3 = .FALSE.
      new_out_1 = .FALSE.

      IF (PRESENT(summation)) THEN
         summation_prv = summation
      ELSE
         summation_prv = .FALSE.
      END IF

      IF (PRESENT(bounds)) THEN
         ALLOCATE (in_tmp_1)
         CALL dbt_crop(tensor_in, in_tmp_1, bounds=bounds, move_data=move_prv)
         new_in_1 = .TRUE.
         move_prv = .TRUE.
      ELSE
         in_tmp_1 => tensor_in
      END IF

      IF (PRESENT(order)) THEN
         CALL reorder_arrays(in_tmp_1%blk_sizes, blk_sizes_in, order=order)
         block_compatible = check_equal(blk_sizes_in, tensor_out%blk_sizes)
      ELSE
         block_compatible = check_equal(in_tmp_1%blk_sizes, tensor_out%blk_sizes)
      END IF

      IF (.NOT. block_compatible) THEN
         ALLOCATE (in_tmp_2, out_tmp_1)
         CALL dbt_make_compatible_blocks(in_tmp_1, tensor_out, in_tmp_2, out_tmp_1, order=order, &
                                         nodata2=.NOT. summation_prv, move_data=move_prv)
         new_in_2 = .TRUE.; new_out_1 = .TRUE.
         move_prv = .TRUE.
      ELSE
         in_tmp_2 => in_tmp_1
         out_tmp_1 => tensor_out
      END IF

      IF (PRESENT(order)) THEN
         ALLOCATE (in_tmp_3)
         CALL dbt_permute_index(in_tmp_2, in_tmp_3, order)
         new_in_3 = .TRUE.
      ELSE
         in_tmp_3 => in_tmp_2
      END IF

      ALLOCATE (map1_in_1(ndims_matrix_row(in_tmp_3)))
      ALLOCATE (map1_in_2(ndims_matrix_column(in_tmp_3)))
      CALL dbt_get_mapping_info(in_tmp_3%nd_index, map1_2d=map1_in_1, map2_2d=map1_in_2)

      ALLOCATE (map2_in_1(ndims_matrix_row(out_tmp_1)))
      ALLOCATE (map2_in_2(ndims_matrix_column(out_tmp_1)))
      CALL dbt_get_mapping_info(out_tmp_1%nd_index, map1_2d=map2_in_1, map2_2d=map2_in_2)

      IF (.NOT. PRESENT(order)) THEN
         IF (array_eq_i(map1_in_1, map2_in_1) .AND. array_eq_i(map1_in_2, map2_in_2)) THEN
            dist_compatible_tas = check_equal(in_tmp_3%nd_dist, out_tmp_1%nd_dist)
         ELSEIF (array_eq_i([map1_in_1, map1_in_2], [map2_in_1, map2_in_2])) THEN
            dist_compatible_tensor = check_equal(in_tmp_3%nd_dist, out_tmp_1%nd_dist)
         END IF
      END IF

      IF (dist_compatible_tas) THEN
         CALL dbt_tas_copy(out_tmp_1%matrix_rep, in_tmp_3%matrix_rep, summation)
         IF (move_prv) CALL dbt_clear(in_tmp_3)
      ELSEIF (dist_compatible_tensor) THEN
         CALL dbt_copy_nocomm(in_tmp_3, out_tmp_1, summation)
         IF (move_prv) CALL dbt_clear(in_tmp_3)
      ELSE
         CALL dbt_reshape(in_tmp_3, out_tmp_1, summation, move_data=move_prv)
      END IF

      IF (new_in_1) THEN
         CALL dbt_destroy(in_tmp_1)
         DEALLOCATE (in_tmp_1)
      END IF

      IF (new_in_2) THEN
         CALL dbt_destroy(in_tmp_2)
         DEALLOCATE (in_tmp_2)
      END IF

      IF (new_in_3) THEN
         CALL dbt_destroy(in_tmp_3)
         DEALLOCATE (in_tmp_3)
      END IF

      IF (new_out_1) THEN
         IF (unit_nr_prv /= 0) THEN
            CALL dbt_write_tensor_dist(out_tmp_1, unit_nr)
         END IF
         CALL dbt_split_copyback(out_tmp_1, tensor_out, summation)
         CALL dbt_destroy(out_tmp_1)
         DEALLOCATE (out_tmp_1)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief copy without communication, requires that both tensors have same process grid and distribution
!> \param summation Whether to sum matrices b = a + b
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_copy_nocomm(tensor_in, tensor_out, summation)
      TYPE(dbt_type), INTENT(INOUT) :: tensor_in
      TYPE(dbt_type), INTENT(INOUT) :: tensor_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation
      TYPE(dbt_iterator_type) :: iter
      INTEGER, DIMENSION(ndims_tensor(tensor_in))  :: ind_nd
      TYPE(block_nd) :: blk_data
      LOGICAL :: found

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_copy_nocomm'
      INTEGER :: handle

      CALL timeset(routineN, handle)
      CPASSERT(tensor_out%valid)

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL dbt_clear(tensor_out)
      ELSE
         CALL dbt_clear(tensor_out)
      END IF

      CALL dbt_reserve_blocks(tensor_in, tensor_out)

      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd)
         CALL dbt_get_block(tensor_in, ind_nd, blk_data, found)
         CPASSERT(found)
         CALL dbt_put_block(tensor_out, ind_nd, blk_data, summation=summation)
         CALL destroy_block(blk_data)
      END DO
      CALL dbt_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief copy matrix to tensor.
!> \param summation tensor_out = tensor_out + matrix_in
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_copy_matrix_to_tensor(matrix_in, tensor_out, summation)
      TYPE(dbcsr_type), TARGET, INTENT(IN)               :: matrix_in
      TYPE(dbt_type), INTENT(INOUT)             :: tensor_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation
      TYPE(dbcsr_type), POINTER                          :: matrix_in_desym

      INTEGER, DIMENSION(2)                              :: ind_2d
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)    :: block_arr
      REAL(KIND=dp), DIMENSION(:, :), POINTER        :: block
      TYPE(dbcsr_iterator_type)                          :: iter
      LOGICAL                                            :: tr

      INTEGER                                            :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_copy_matrix_to_tensor'

      CALL timeset(routineN, handle)
      CPASSERT(tensor_out%valid)

      NULLIFY (block)

      IF (dbcsr_has_symmetry(matrix_in)) THEN
         ALLOCATE (matrix_in_desym)
         CALL dbcsr_desymmetrize(matrix_in, matrix_in_desym)
      ELSE
         matrix_in_desym => matrix_in
      END IF

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL dbt_clear(tensor_out)
      ELSE
         CALL dbt_clear(tensor_out)
      END IF

      CALL dbt_reserve_blocks(matrix_in_desym, tensor_out)

      CALL dbcsr_iterator_start(iter, matrix_in_desym)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ind_2d(1), ind_2d(2), block, tr)
         CALL allocate_any(block_arr, source=block)
         CALL dbt_put_block(tensor_out, ind_2d, SHAPE(block_arr), block_arr, summation=summation)
         DEALLOCATE (block_arr)
      END DO
      CALL dbcsr_iterator_stop(iter)

      IF (dbcsr_has_symmetry(matrix_in)) THEN
         CALL dbcsr_release(matrix_in_desym)
         DEALLOCATE (matrix_in_desym)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief copy tensor to matrix
!> \param summation matrix_out = matrix_out + tensor_in
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_copy_tensor_to_matrix(tensor_in, matrix_out, summation)
      TYPE(dbt_type), INTENT(INOUT)      :: tensor_in
      TYPE(dbcsr_type), INTENT(INOUT)             :: matrix_out
      LOGICAL, INTENT(IN), OPTIONAL          :: summation
      TYPE(dbt_iterator_type)            :: iter
      INTEGER                                :: handle
      INTEGER, DIMENSION(2)                  :: ind_2d
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: block
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_copy_tensor_to_matrix'
      LOGICAL :: found

      CALL timeset(routineN, handle)

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL dbcsr_clear(matrix_out)
      ELSE
         CALL dbcsr_clear(matrix_out)
      END IF

      CALL dbt_reserve_blocks(tensor_in, matrix_out)

      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_2d)
         IF (dbcsr_has_symmetry(matrix_out) .AND. checker_tr(ind_2d(1), ind_2d(2))) CYCLE

         CALL dbt_get_block(tensor_in, ind_2d, block, found)
         CPASSERT(found)

         IF (dbcsr_has_symmetry(matrix_out) .AND. ind_2d(1) > ind_2d(2)) THEN
            CALL dbcsr_put_block(matrix_out, ind_2d(2), ind_2d(1), TRANSPOSE(block), summation=summation)
         ELSE
            CALL dbcsr_put_block(matrix_out, ind_2d(1), ind_2d(2), block, summation=summation)
         END IF
         DEALLOCATE (block)
      END DO
      CALL dbt_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Contract tensors by multiplying matrix representations.
!>        tensor_3(map_1, map_2) := alpha * tensor_1(notcontract_1, contract_1)
!>        * tensor_2(contract_2, notcontract_2)
!>        + beta * tensor_3(map_1, map_2)
!>
!> \note
!>      note 1: block sizes of the corresponding indices need to be the same in all tensors.
!>
!>      note 2: for best performance the tensors should have been created in matrix layouts
!>      compatible with the contraction, e.g. tensor_1 should have been created with either
!>      map1_2d == contract_1 and map2_2d == notcontract_1 or map1_2d == notcontract_1 and
!>      map2_2d == contract_1 (the same with tensor_2 and contract_2 / notcontract_2 and with
!>      tensor_3 and map_1 / map_2).
!>      Furthermore the two largest tensors involved in the contraction should map both to either
!>      tall or short matrices: the largest matrix dimension should be "on the same side"
!>      and should have identical distribution (which is always the case if the distributions were
!>      obtained with dbt_default_distvec).
!>
!>      note 3: if the same tensor occurs in multiple contractions, a different tensor object should
!>      be created for each contraction and the data should be copied between the tensors by use of
!>      dbt_copy. If the same tensor object is used in multiple contractions,
!>       matrix layouts are not compatible for all contractions (see note 2).
!>
!>      note 4: automatic optimizations are enabled by using the feature of batched contraction, see
!>      dbt_batched_contract_init, dbt_batched_contract_finalize.
!>      The arguments bounds_1, bounds_2, bounds_3 give the index ranges of the batches.
!>
!> \param tensor_1 first tensor (in)
!> \param tensor_2 second tensor (in)
!> \param contract_1 indices of tensor_1 to contract
!> \param contract_2 indices of tensor_2 to contract (1:1 with contract_1)
!> \param map_1 which indices of tensor_3 map to non-contracted indices of tensor_1 (1:1 with notcontract_1)
!> \param map_2 which indices of tensor_3 map to non-contracted indices of tensor_2 (1:1 with notcontract_2)
!> \param notcontract_1 indices of tensor_1 not to contract
!> \param notcontract_2 indices of tensor_2 not to contract
!> \param tensor_3 contracted tensor (out)
!> \param bounds_1 bounds corresponding to contract_1 AKA contract_2:
!>                 start and end index of an index range over which to contract.
!>                 For use in batched contraction.
!> \param bounds_2 bounds corresponding to notcontract_1: start and end index of an index range.
!>                 For use in batched contraction.
!> \param bounds_3 bounds corresponding to notcontract_2: start and end index of an index range.
!>                 For use in batched contraction.
!> \param optimize_dist Whether distribution should be optimized internally. In the current
!>                      implementation this guarantees optimal parameters only for dense matrices.
!> \param pgrid_opt_1 Optionally return optimal process grid for tensor_1.
!>                    This can be used to choose optimal process grids for subsequent tensor
!>                    contractions with tensors of similar shape and sparsity. Under some conditions,
!>                    pgrid_opt_1 can not be returned, in this case the pointer is not associated.
!> \param pgrid_opt_2 Optionally return optimal process grid for tensor_2.
!> \param pgrid_opt_3 Optionally return optimal process grid for tensor_3.
!> \param filter_eps As in DBM mm
!> \param flop As in DBM mm
!> \param move_data memory optimization: transfer data such that tensor_1 and tensor_2 are empty on return
!> \param retain_sparsity enforce the sparsity pattern of the existing tensor_3; default is no
!> \param unit_nr output unit for logging
!>                       set it to -1 on ranks that should not write (and any valid unit number on
!>                       ranks that should write output) if 0 on ALL ranks, no output is written
!> \param log_verbose verbose logging (for testing only)
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_contract(alpha, tensor_1, tensor_2, beta, tensor_3, &
                           contract_1, notcontract_1, &
                           contract_2, notcontract_2, &
                           map_1, map_2, &
                           bounds_1, bounds_2, bounds_3, &
                           optimize_dist, pgrid_opt_1, pgrid_opt_2, pgrid_opt_3, &
                           filter_eps, flop, move_data, retain_sparsity, unit_nr, log_verbose)
      REAL(dp), INTENT(IN)            :: alpha
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_1
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_2
      REAL(dp), INTENT(IN)            :: beta
      INTEGER, DIMENSION(:), INTENT(IN)              :: contract_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: contract_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: map_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: map_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: notcontract_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: notcontract_2
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_3
      INTEGER, DIMENSION(2, SIZE(contract_1)), &
         INTENT(IN), OPTIONAL                        :: bounds_1
      INTEGER, DIMENSION(2, SIZE(notcontract_1)), &
         INTENT(IN), OPTIONAL                        :: bounds_2
      INTEGER, DIMENSION(2, SIZE(notcontract_2)), &
         INTENT(IN), OPTIONAL                        :: bounds_3
      LOGICAL, INTENT(IN), OPTIONAL                  :: optimize_dist
      TYPE(dbt_pgrid_type), INTENT(OUT), &
         POINTER, OPTIONAL                           :: pgrid_opt_1
      TYPE(dbt_pgrid_type), INTENT(OUT), &
         POINTER, OPTIONAL                           :: pgrid_opt_2
      TYPE(dbt_pgrid_type), INTENT(OUT), &
         POINTER, OPTIONAL                           :: pgrid_opt_3
      REAL(KIND=dp), INTENT(IN), OPTIONAL        :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL     :: flop
      LOGICAL, INTENT(IN), OPTIONAL                  :: move_data
      LOGICAL, INTENT(IN), OPTIONAL                  :: retain_sparsity
      INTEGER, OPTIONAL, INTENT(IN)                  :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                  :: log_verbose

      INTEGER                     :: handle

      CALL mp_sync(tensor_1%pgrid%mp_comm_2d)
      CALL timeset("dbt_total", handle)
      CALL dbt_contract_expert(alpha, tensor_1, tensor_2, beta, tensor_3, &
                               contract_1, notcontract_1, &
                               contract_2, notcontract_2, &
                               map_1, map_2, &
                               bounds_1=bounds_1, &
                               bounds_2=bounds_2, &
                               bounds_3=bounds_3, &
                               optimize_dist=optimize_dist, &
                               pgrid_opt_1=pgrid_opt_1, &
                               pgrid_opt_2=pgrid_opt_2, &
                               pgrid_opt_3=pgrid_opt_3, &
                               filter_eps=filter_eps, &
                               flop=flop, &
                               move_data=move_data, &
                               retain_sparsity=retain_sparsity, &
                               unit_nr=unit_nr, &
                               log_verbose=log_verbose)
      CALL mp_sync(tensor_1%pgrid%mp_comm_2d)
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief expert routine for tensor contraction. For internal use only.
!> \param nblks_local number of local blocks on this MPI rank
!> \param result_index get indices of non-zero tensor blocks for tensor_3 without actually
!>                     performing contraction this is an estimate based on block norm multiplication
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_contract_expert(alpha, tensor_1, tensor_2, beta, tensor_3, &
                                  contract_1, notcontract_1, &
                                  contract_2, notcontract_2, &
                                  map_1, map_2, &
                                  bounds_1, bounds_2, bounds_3, &
                                  optimize_dist, pgrid_opt_1, pgrid_opt_2, pgrid_opt_3, &
                                  filter_eps, flop, move_data, retain_sparsity, &
                                  nblks_local, result_index, unit_nr, log_verbose)
      REAL(dp), INTENT(IN)            :: alpha
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_1
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_2
      REAL(dp), INTENT(IN)            :: beta
      INTEGER, DIMENSION(:), INTENT(IN)              :: contract_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: contract_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: map_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: map_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: notcontract_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: notcontract_2
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_3
      INTEGER, DIMENSION(2, SIZE(contract_1)), &
         INTENT(IN), OPTIONAL                        :: bounds_1
      INTEGER, DIMENSION(2, SIZE(notcontract_1)), &
         INTENT(IN), OPTIONAL                        :: bounds_2
      INTEGER, DIMENSION(2, SIZE(notcontract_2)), &
         INTENT(IN), OPTIONAL                        :: bounds_3
      LOGICAL, INTENT(IN), OPTIONAL                  :: optimize_dist
      TYPE(dbt_pgrid_type), INTENT(OUT), &
         POINTER, OPTIONAL                           :: pgrid_opt_1
      TYPE(dbt_pgrid_type), INTENT(OUT), &
         POINTER, OPTIONAL                           :: pgrid_opt_2
      TYPE(dbt_pgrid_type), INTENT(OUT), &
         POINTER, OPTIONAL                           :: pgrid_opt_3
      REAL(KIND=dp), INTENT(IN), OPTIONAL        :: filter_eps
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL     :: flop
      LOGICAL, INTENT(IN), OPTIONAL                  :: move_data
      LOGICAL, INTENT(IN), OPTIONAL                  :: retain_sparsity
      INTEGER, INTENT(OUT), OPTIONAL                 :: nblks_local
      INTEGER, DIMENSION(dbt_max_nblks_local(tensor_3), ndims_tensor(tensor_3)), &
         OPTIONAL, INTENT(OUT)                       :: result_index
      INTEGER, OPTIONAL, INTENT(IN)                  :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                  :: log_verbose

      TYPE(dbt_type), POINTER                    :: tensor_contr_1, tensor_contr_2, tensor_contr_3
      TYPE(dbt_type), TARGET                     :: tensor_algn_1, tensor_algn_2, tensor_algn_3
      TYPE(dbt_type), POINTER                    :: tensor_crop_1, tensor_crop_2
      TYPE(dbt_type), POINTER                    :: tensor_small, tensor_large

      INTEGER(int_8), DIMENSION(:, :), ALLOCATABLE  :: result_index_2d
      LOGICAL                                        :: assert_stmt, tensors_remapped
      INTEGER                                        :: max_mm_dim, max_tensor, mp_comm, &
                                                        iblk, nblk, unit_nr_prv, ref_tensor, mp_comm_opt, &
                                                        handle
      INTEGER, DIMENSION(SIZE(contract_1))           :: contract_1_mod
      INTEGER, DIMENSION(SIZE(notcontract_1))        :: notcontract_1_mod
      INTEGER, DIMENSION(SIZE(contract_2))           :: contract_2_mod
      INTEGER, DIMENSION(SIZE(notcontract_2))        :: notcontract_2_mod
      INTEGER, DIMENSION(SIZE(map_1))                :: map_1_mod
      INTEGER, DIMENSION(SIZE(map_2))                :: map_2_mod
      LOGICAL                                        :: trans_1, trans_2, trans_3
      LOGICAL                                        :: new_1, new_2, new_3, move_data_1, move_data_2
      INTEGER                                        :: ndims1, ndims2, ndims3
      INTEGER                                        :: occ_1, occ_2
      INTEGER, DIMENSION(:), ALLOCATABLE             :: dims1, dims2, dims3

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_contract'
      CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE    :: indchar1, indchar2, indchar3, indchar1_mod, &
                                                        indchar2_mod, indchar3_mod
      CHARACTER(LEN=1), DIMENSION(15) :: alph = &
                                         ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o']
      INTEGER, DIMENSION(2, ndims_tensor(tensor_1)) :: bounds_t1
      INTEGER, DIMENSION(2, ndims_tensor(tensor_2)) :: bounds_t2
      LOGICAL                                        :: do_crop_1, do_crop_2, do_write_3, nodata_3, do_batched, pgrid_changed, &
                                                        pgrid_changed_any, do_change_pgrid(2)
      TYPE(dbt_tas_split_info)                     :: split_opt, split, split_opt_avg
      INTEGER, DIMENSION(2) :: pdims_2d_opt, pdims_2d, pcoord_2d, pdims_sub, pdims_sub_opt
      LOGICAL, DIMENSION(2) :: periods_2d
      REAL(dp) :: pdim_ratio, pdim_ratio_opt

      NULLIFY (tensor_contr_1, tensor_contr_2, tensor_contr_3, tensor_crop_1, tensor_crop_2, &
               tensor_small)

      CALL timeset(routineN, handle)

      CPASSERT(tensor_1%valid)
      CPASSERT(tensor_2%valid)
      CPASSERT(tensor_3%valid)

      assert_stmt = SIZE(contract_1) .EQ. SIZE(contract_2)
      CPASSERT(assert_stmt)

      assert_stmt = SIZE(map_1) .EQ. SIZE(notcontract_1)
      CPASSERT(assert_stmt)

      assert_stmt = SIZE(map_2) .EQ. SIZE(notcontract_2)
      CPASSERT(assert_stmt)

      assert_stmt = SIZE(notcontract_1) + SIZE(contract_1) .EQ. ndims_tensor(tensor_1)
      CPASSERT(assert_stmt)

      assert_stmt = SIZE(notcontract_2) + SIZE(contract_2) .EQ. ndims_tensor(tensor_2)
      CPASSERT(assert_stmt)

      assert_stmt = SIZE(map_1) + SIZE(map_2) .EQ. ndims_tensor(tensor_3)
      CPASSERT(assert_stmt)

      unit_nr_prv = prep_output_unit(unit_nr)

      IF (PRESENT(flop)) flop = 0
      IF (PRESENT(result_index)) result_index = 0
      IF (PRESENT(nblks_local)) nblks_local = 0

      IF (PRESENT(move_data)) THEN
         move_data_1 = move_data
         move_data_2 = move_data
      ELSE
         move_data_1 = .FALSE.
         move_data_2 = .FALSE.
      END IF

      nodata_3 = .TRUE.
      IF (PRESENT(retain_sparsity)) THEN
         IF (retain_sparsity) nodata_3 = .FALSE.
      END IF

      CALL dbt_map_bounds_to_tensors(tensor_1, tensor_2, &
                                     contract_1, notcontract_1, &
                                     contract_2, notcontract_2, &
                                     bounds_t1, bounds_t2, &
                                     bounds_1=bounds_1, bounds_2=bounds_2, bounds_3=bounds_3, &
                                     do_crop_1=do_crop_1, do_crop_2=do_crop_2)

      IF (do_crop_1) THEN
         ALLOCATE (tensor_crop_1)
         CALL dbt_crop(tensor_1, tensor_crop_1, bounds_t1, move_data=move_data_1)
         move_data_1 = .TRUE.
      ELSE
         tensor_crop_1 => tensor_1
      END IF

      IF (do_crop_2) THEN
         ALLOCATE (tensor_crop_2)
         CALL dbt_crop(tensor_2, tensor_crop_2, bounds_t2, move_data=move_data_2)
         move_data_2 = .TRUE.
      ELSE
         tensor_crop_2 => tensor_2
      END IF

      ! shortcut for empty tensors
      ! this is needed to avoid unnecessary work in case user contracts different portions of a
      ! tensor consecutively to save memory
      mp_comm = tensor_crop_1%pgrid%mp_comm_2d
      occ_1 = dbt_get_num_blocks(tensor_crop_1)
      CALL mp_max(occ_1, mp_comm)
      occ_2 = dbt_get_num_blocks(tensor_crop_2)
      CALL mp_max(occ_2, mp_comm)

      IF (occ_1 == 0 .OR. occ_2 == 0) THEN
         CALL dbt_scale(tensor_3, beta)
         IF (do_crop_1) THEN
            CALL dbt_destroy(tensor_crop_1)
            DEALLOCATE (tensor_crop_1)
         END IF
         IF (do_crop_2) THEN
            CALL dbt_destroy(tensor_crop_2)
            DEALLOCATE (tensor_crop_2)
         END IF

         CALL timestop(handle)
         RETURN
      END IF

      IF (unit_nr_prv /= 0) THEN
         IF (unit_nr_prv > 0) THEN
            WRITE (unit_nr_prv, '(A)') repeat("-", 80)
            WRITE (unit_nr_prv, '(A,1X,A,1X,A,1X,A,1X,A,1X,A)') "DBT TENSOR CONTRACTION:", &
               TRIM(tensor_crop_1%name), 'x', TRIM(tensor_crop_2%name), '=', TRIM(tensor_3%name)
            WRITE (unit_nr_prv, '(A)') repeat("-", 80)
         END IF
         CALL dbt_write_tensor_info(tensor_crop_1, unit_nr_prv, full_info=log_verbose)
         CALL dbt_write_tensor_dist(tensor_crop_1, unit_nr_prv)
         CALL dbt_write_tensor_info(tensor_crop_2, unit_nr_prv, full_info=log_verbose)
         CALL dbt_write_tensor_dist(tensor_crop_2, unit_nr_prv)
      END IF

      ! align tensor index with data, tensor data is not modified
      ndims1 = ndims_tensor(tensor_crop_1)
      ndims2 = ndims_tensor(tensor_crop_2)
      ndims3 = ndims_tensor(tensor_3)
      ALLOCATE (indchar1(ndims1), indchar1_mod(ndims1))
      ALLOCATE (indchar2(ndims2), indchar2_mod(ndims2))
      ALLOCATE (indchar3(ndims3), indchar3_mod(ndims3))

      ! labeling tensor index with letters

      indchar1([notcontract_1, contract_1]) = alph(1:ndims1) ! arb. choice
      indchar2(notcontract_2) = alph(ndims1 + 1:ndims1 + SIZE(notcontract_2)) ! arb. choice
      indchar2(contract_2) = indchar1(contract_1)
      indchar3(map_1) = indchar1(notcontract_1)
      indchar3(map_2) = indchar2(notcontract_2)

      IF (unit_nr_prv /= 0) CALL dbt_print_contraction_index(tensor_crop_1, indchar1, &
                                                             tensor_crop_2, indchar2, &
                                                             tensor_3, indchar3, unit_nr_prv)
      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, '(T2,A)') "aligning tensor index with data"
      END IF

      CALL align_tensor(tensor_crop_1, contract_1, notcontract_1, &
                        tensor_algn_1, contract_1_mod, notcontract_1_mod, indchar1, indchar1_mod)

      CALL align_tensor(tensor_crop_2, contract_2, notcontract_2, &
                        tensor_algn_2, contract_2_mod, notcontract_2_mod, indchar2, indchar2_mod)

      CALL align_tensor(tensor_3, map_1, map_2, &
                        tensor_algn_3, map_1_mod, map_2_mod, indchar3, indchar3_mod)

      IF (unit_nr_prv /= 0) CALL dbt_print_contraction_index(tensor_algn_1, indchar1_mod, &
                                                             tensor_algn_2, indchar2_mod, &
                                                             tensor_algn_3, indchar3_mod, unit_nr_prv)

      ALLOCATE (dims1(ndims1))
      ALLOCATE (dims2(ndims2))
      ALLOCATE (dims3(ndims3))

      ! ideally we should consider block sizes and occupancy to measure tensor sizes but current solution should work for most
      ! cases and is more elegant. Note that we can not easily consider occupancy since it is unknown for result tensor
      CALL blk_dims_tensor(tensor_crop_1, dims1)
      CALL blk_dims_tensor(tensor_crop_2, dims2)
      CALL blk_dims_tensor(tensor_3, dims3)

      max_mm_dim = MAXLOC([PRODUCT(INT(dims1(notcontract_1), int_8)), &
                           PRODUCT(INT(dims1(contract_1), int_8)), &
                           PRODUCT(INT(dims2(notcontract_2), int_8))], DIM=1)
      max_tensor = MAXLOC([PRODUCT(INT(dims1, int_8)), PRODUCT(INT(dims2, int_8)), PRODUCT(INT(dims3, int_8))], DIM=1)
      SELECT CASE (max_mm_dim)
      CASE (1)
         IF (unit_nr_prv > 0) THEN
            WRITE (unit_nr_prv, '(T2,A)') "large tensors: 1, 3; small tensor: 2"
            WRITE (unit_nr_prv, '(T2,A)') "sorting contraction indices"
         END IF
         CALL index_linked_sort(contract_1_mod, contract_2_mod)
         CALL index_linked_sort(map_2_mod, notcontract_2_mod)
         SELECT CASE (max_tensor)
         CASE (1)
            CALL index_linked_sort(notcontract_1_mod, map_1_mod)
         CASE (3)
            CALL index_linked_sort(map_1_mod, notcontract_1_mod)
         CASE DEFAULT
            CPABORT("should not happen")
         END SELECT

         CALL reshape_mm_compatible(tensor_algn_1, tensor_algn_3, tensor_contr_1, tensor_contr_3, &
                                    contract_1_mod, notcontract_1_mod, map_2_mod, map_1_mod, &
                                    trans_1, trans_3, new_1, new_3, ref_tensor, nodata2=nodata_3, optimize_dist=optimize_dist, &
                                    move_data_1=move_data_1, unit_nr=unit_nr_prv)

         CALL reshape_mm_small(tensor_algn_2, contract_2_mod, notcontract_2_mod, tensor_contr_2, trans_2, &
                               new_2, move_data=move_data_2, unit_nr=unit_nr_prv)

         SELECT CASE (ref_tensor)
         CASE (1)
            tensor_large => tensor_contr_1
         CASE (2)
            tensor_large => tensor_contr_3
         END SELECT
         tensor_small => tensor_contr_2

      CASE (2)
         IF (unit_nr_prv > 0) THEN
            WRITE (unit_nr_prv, '(T2,A)') "large tensors: 1, 2; small tensor: 3"
            WRITE (unit_nr_prv, '(T2,A)') "sorting contraction indices"
         END IF

         CALL index_linked_sort(notcontract_1_mod, map_1_mod)
         CALL index_linked_sort(notcontract_2_mod, map_2_mod)
         SELECT CASE (max_tensor)
         CASE (1)
            CALL index_linked_sort(contract_1_mod, contract_2_mod)
         CASE (2)
            CALL index_linked_sort(contract_2_mod, contract_1_mod)
         CASE DEFAULT
            CPABORT("should not happen")
         END SELECT

         CALL reshape_mm_compatible(tensor_algn_1, tensor_algn_2, tensor_contr_1, tensor_contr_2, &
                                    notcontract_1_mod, contract_1_mod, notcontract_2_mod, contract_2_mod, &
                                    trans_1, trans_2, new_1, new_2, ref_tensor, optimize_dist=optimize_dist, &
                                    move_data_1=move_data_1, move_data_2=move_data_2, unit_nr=unit_nr_prv)
         trans_1 = .NOT. trans_1

         CALL reshape_mm_small(tensor_algn_3, map_1_mod, map_2_mod, tensor_contr_3, trans_3, &
                               new_3, nodata=nodata_3, unit_nr=unit_nr_prv)

         SELECT CASE (ref_tensor)
         CASE (1)
            tensor_large => tensor_contr_1
         CASE (2)
            tensor_large => tensor_contr_2
         END SELECT
         tensor_small => tensor_contr_3

      CASE (3)
         IF (unit_nr_prv > 0) THEN
            WRITE (unit_nr_prv, '(T2,A)') "large tensors: 2, 3; small tensor: 1"
            WRITE (unit_nr_prv, '(T2,A)') "sorting contraction indices"
         END IF
         CALL index_linked_sort(map_1_mod, notcontract_1_mod)
         CALL index_linked_sort(contract_2_mod, contract_1_mod)
         SELECT CASE (max_tensor)
         CASE (2)
            CALL index_linked_sort(notcontract_2_mod, map_2_mod)
         CASE (3)
            CALL index_linked_sort(map_2_mod, notcontract_2_mod)
         CASE DEFAULT
            CPABORT("should not happen")
         END SELECT

         CALL reshape_mm_compatible(tensor_algn_2, tensor_algn_3, tensor_contr_2, tensor_contr_3, &
                                    contract_2_mod, notcontract_2_mod, map_1_mod, map_2_mod, &
                                    trans_2, trans_3, new_2, new_3, ref_tensor, nodata2=nodata_3, optimize_dist=optimize_dist, &
                                    move_data_1=move_data_2, unit_nr=unit_nr_prv)

         trans_2 = .NOT. trans_2
         trans_3 = .NOT. trans_3

         CALL reshape_mm_small(tensor_algn_1, notcontract_1_mod, contract_1_mod, tensor_contr_1, &
                               trans_1, new_1, move_data=move_data_1, unit_nr=unit_nr_prv)

         SELECT CASE (ref_tensor)
         CASE (1)
            tensor_large => tensor_contr_2
         CASE (2)
            tensor_large => tensor_contr_3
         END SELECT
         tensor_small => tensor_contr_1

      END SELECT

      IF (unit_nr_prv /= 0) CALL dbt_print_contraction_index(tensor_contr_1, indchar1_mod, &
                                                             tensor_contr_2, indchar2_mod, &
                                                             tensor_contr_3, indchar3_mod, unit_nr_prv)
      IF (unit_nr_prv /= 0) THEN
         IF (new_1) CALL dbt_write_tensor_info(tensor_contr_1, unit_nr_prv, full_info=log_verbose)
         IF (new_1) CALL dbt_write_tensor_dist(tensor_contr_1, unit_nr_prv)
         IF (new_2) CALL dbt_write_tensor_info(tensor_contr_2, unit_nr_prv, full_info=log_verbose)
         IF (new_2) CALL dbt_write_tensor_dist(tensor_contr_2, unit_nr_prv)
      END IF

      IF (.NOT. PRESENT(result_index)) THEN
         CALL dbt_tas_multiply(trans_1, trans_2, trans_3, alpha, &
                               tensor_contr_1%matrix_rep, tensor_contr_2%matrix_rep, &
                               beta, &
                               tensor_contr_3%matrix_rep, filter_eps=filter_eps, flop=flop, &
                               unit_nr=unit_nr_prv, log_verbose=log_verbose, &
                               split_opt=split_opt, &
                               move_data_a=move_data_1, move_data_b=move_data_2, retain_sparsity=retain_sparsity)
      ELSE

         CALL dbt_tas_result_index(trans_1, trans_2, trans_3, tensor_contr_1%matrix_rep, tensor_contr_2%matrix_rep, &
                                   tensor_contr_3%matrix_rep, filter_eps=filter_eps, blk_ind=result_index_2d)

         nblk = SIZE(result_index_2d, 1)
         IF (PRESENT(nblks_local)) nblks_local = nblk
         IF (SIZE(result_index, 1) < nblk) THEN
            CALL cp_abort(__LOCATION__, &
        "allocated size of `result_index` is too small. This error occurs due to a high load imbalance of distributed tensor data.")
         END IF

         DO iblk = 1, nblk
            result_index(iblk, :) = get_nd_indices_tensor(tensor_contr_3%nd_index_blk, result_index_2d(iblk, :))
         END DO

         IF (new_1) THEN
            CALL dbt_destroy(tensor_contr_1)
            DEALLOCATE (tensor_contr_1)
         END IF
         IF (new_2) THEN
            CALL dbt_destroy(tensor_contr_2)
            DEALLOCATE (tensor_contr_2)
         END IF
         IF (new_3) THEN
            CALL dbt_destroy(tensor_contr_3)
            DEALLOCATE (tensor_contr_3)
         END IF
         IF (do_crop_1) THEN
            CALL dbt_destroy(tensor_crop_1)
            DEALLOCATE (tensor_crop_1)
         END IF
         IF (do_crop_2) THEN
            CALL dbt_destroy(tensor_crop_2)
            DEALLOCATE (tensor_crop_2)
         END IF

         CALL dbt_destroy(tensor_algn_1)
         CALL dbt_destroy(tensor_algn_2)
         CALL dbt_destroy(tensor_algn_3)

         CALL timestop(handle)
         RETURN
      END IF

      IF (PRESENT(pgrid_opt_1)) THEN
         IF (.NOT. new_1) THEN
            ALLOCATE (pgrid_opt_1)
            pgrid_opt_1 = opt_pgrid(tensor_1, split_opt)
         END IF
      END IF

      IF (PRESENT(pgrid_opt_2)) THEN
         IF (.NOT. new_2) THEN
            ALLOCATE (pgrid_opt_2)
            pgrid_opt_2 = opt_pgrid(tensor_2, split_opt)
         END IF
      END IF

      IF (PRESENT(pgrid_opt_3)) THEN
         IF (.NOT. new_3) THEN
            ALLOCATE (pgrid_opt_3)
            pgrid_opt_3 = opt_pgrid(tensor_3, split_opt)
         END IF
      END IF

      do_batched = tensor_small%matrix_rep%do_batched > 0

      tensors_remapped = .FALSE.
      IF (new_1 .OR. new_2 .OR. new_3) tensors_remapped = .TRUE.

      IF (tensors_remapped .AND. do_batched) THEN
         CALL cp_warn(__LOCATION__, &
                      "Internal process grid optimization disabled because tensors are not in contraction-compatible format")
      END IF

      CALL mp_environ(tensor_large%pgrid%mp_comm_2d, 2, pdims_2d, pcoord_2d, periods_2d)

      ! optimize process grid during batched contraction
      do_change_pgrid(:) = .FALSE.
      IF ((.NOT. tensors_remapped) .AND. do_batched) THEN
         ASSOCIATE (storage => tensor_small%contraction_storage)
            CPASSERT(storage%static)
            split = dbt_tas_info(tensor_large%matrix_rep)
            do_change_pgrid(:) = &
               update_contraction_storage(storage, split_opt, split)

            IF (ANY(do_change_pgrid)) THEN
               mp_comm_opt = dbt_tas_mp_comm(tensor_small%pgrid%mp_comm_2d, split_opt%split_rowcol, NINT(storage%nsplit_avg))
               CALL dbt_tas_create_split(split_opt_avg, mp_comm_opt, split_opt%split_rowcol, &
                                         NINT(storage%nsplit_avg), own_comm=.TRUE.)
               CALL mp_environ(split_opt_avg%mp_comm, 2, pdims_2d_opt, pcoord_2d, periods_2d)
            END IF

         END ASSOCIATE

         IF (do_change_pgrid(1) .AND. .NOT. do_change_pgrid(2)) THEN
            ! check if new grid has better subgrid, if not there is no need to change process grid
            CALL mp_environ(split_opt_avg%mp_comm_group, 2, pdims_sub_opt, pcoord_2d, periods_2d)
            CALL mp_environ(split%mp_comm_group, 2, pdims_sub, pcoord_2d, periods_2d)

            pdim_ratio = MAXVAL(REAL(pdims_sub, dp))/MINVAL(pdims_sub)
            pdim_ratio_opt = MAXVAL(REAL(pdims_sub_opt, dp))/MINVAL(pdims_sub_opt)
            IF (pdim_ratio/pdim_ratio_opt <= default_pdims_accept_ratio**2) THEN
               do_change_pgrid(1) = .FALSE.
               CALL dbt_tas_release_info(split_opt_avg)
            END IF
         END IF
      END IF

      IF (unit_nr_prv /= 0) THEN
         do_write_3 = .TRUE.
         IF (tensor_contr_3%matrix_rep%do_batched > 0) THEN
            IF (tensor_contr_3%matrix_rep%mm_storage%batched_out) do_write_3 = .FALSE.
         END IF
         IF (do_write_3) THEN
            CALL dbt_write_tensor_info(tensor_contr_3, unit_nr_prv, full_info=log_verbose)
            CALL dbt_write_tensor_dist(tensor_contr_3, unit_nr_prv)
         END IF
      END IF

      IF (new_3) THEN
         ! need redistribute if we created new tensor for tensor 3
         CALL dbt_scale(tensor_algn_3, beta)
         CALL dbt_copy_expert(tensor_contr_3, tensor_algn_3, summation=.TRUE., move_data=.TRUE.)
         IF (PRESENT(filter_eps)) CALL dbt_filter(tensor_algn_3, filter_eps)
         ! tensor_3 automatically has correct data because tensor_algn_3 contains a matrix
         ! pointer to data of tensor_3
      END IF

      ! transfer contraction storage
      CALL dbt_copy_contraction_storage(tensor_contr_1, tensor_1)
      CALL dbt_copy_contraction_storage(tensor_contr_2, tensor_2)
      CALL dbt_copy_contraction_storage(tensor_contr_3, tensor_3)

      IF (unit_nr_prv /= 0) THEN
         IF (new_3 .AND. do_write_3) CALL dbt_write_tensor_info(tensor_3, unit_nr_prv, full_info=log_verbose)
         IF (new_3 .AND. do_write_3) CALL dbt_write_tensor_dist(tensor_3, unit_nr_prv)
      END IF

      CALL dbt_destroy(tensor_algn_1)
      CALL dbt_destroy(tensor_algn_2)
      CALL dbt_destroy(tensor_algn_3)

      IF (do_crop_1) THEN
         CALL dbt_destroy(tensor_crop_1)
         DEALLOCATE (tensor_crop_1)
      END IF

      IF (do_crop_2) THEN
         CALL dbt_destroy(tensor_crop_2)
         DEALLOCATE (tensor_crop_2)
      END IF

      IF (new_1) THEN
         CALL dbt_destroy(tensor_contr_1)
         DEALLOCATE (tensor_contr_1)
      END IF
      IF (new_2) THEN
         CALL dbt_destroy(tensor_contr_2)
         DEALLOCATE (tensor_contr_2)
      END IF
      IF (new_3) THEN
         CALL dbt_destroy(tensor_contr_3)
         DEALLOCATE (tensor_contr_3)
      END IF

      IF (PRESENT(move_data)) THEN
         IF (move_data) THEN
            CALL dbt_clear(tensor_1)
            CALL dbt_clear(tensor_2)
         END IF
      END IF

      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, '(A)') repeat("-", 80)
         WRITE (unit_nr_prv, '(A)') "TENSOR CONTRACTION DONE"
         WRITE (unit_nr_prv, '(A)') repeat("-", 80)
      END IF

      IF (ANY(do_change_pgrid)) THEN
         pgrid_changed_any = .FALSE.
         SELECT CASE (max_mm_dim)
         CASE (1)
            IF (ALLOCATED(tensor_1%contraction_storage) .AND. ALLOCATED(tensor_3%contraction_storage)) THEN
               CALL dbt_change_pgrid_2d(tensor_1, tensor_1%pgrid%mp_comm_2d, pdims=pdims_2d_opt, &
                                        nsplit=split_opt_avg%ngroup, dimsplit=split_opt_avg%split_rowcol, &
                                        pgrid_changed=pgrid_changed, &
                                        unit_nr=unit_nr_prv)
               IF (pgrid_changed) pgrid_changed_any = .TRUE.
               CALL dbt_change_pgrid_2d(tensor_3, tensor_3%pgrid%mp_comm_2d, pdims=pdims_2d_opt, &
                                        nsplit=split_opt_avg%ngroup, dimsplit=split_opt_avg%split_rowcol, &
                                        pgrid_changed=pgrid_changed, &
                                        unit_nr=unit_nr_prv)
               IF (pgrid_changed) pgrid_changed_any = .TRUE.
            END IF
            IF (pgrid_changed_any) THEN
               IF (tensor_2%matrix_rep%do_batched == 3) THEN
                  ! set flag that process grid has been optimized to make sure that no grid optimizations are done
                  ! in TAS multiply algorithm
                  CALL dbt_tas_batched_mm_complete(tensor_2%matrix_rep)
               END IF
            END IF
         CASE (2)
            IF (ALLOCATED(tensor_1%contraction_storage) .AND. ALLOCATED(tensor_2%contraction_storage)) THEN
               CALL dbt_change_pgrid_2d(tensor_1, tensor_1%pgrid%mp_comm_2d, pdims=pdims_2d_opt, &
                                        nsplit=split_opt_avg%ngroup, dimsplit=split_opt_avg%split_rowcol, &
                                        pgrid_changed=pgrid_changed, &
                                        unit_nr=unit_nr_prv)
               IF (pgrid_changed) pgrid_changed_any = .TRUE.
               CALL dbt_change_pgrid_2d(tensor_2, tensor_2%pgrid%mp_comm_2d, pdims=pdims_2d_opt, &
                                        nsplit=split_opt_avg%ngroup, dimsplit=split_opt_avg%split_rowcol, &
                                        pgrid_changed=pgrid_changed, &
                                        unit_nr=unit_nr_prv)
               IF (pgrid_changed) pgrid_changed_any = .TRUE.
            END IF
            IF (pgrid_changed_any) THEN
               IF (tensor_3%matrix_rep%do_batched == 3) THEN
                  CALL dbt_tas_batched_mm_complete(tensor_3%matrix_rep)
               END IF
            END IF
         CASE (3)
            IF (ALLOCATED(tensor_2%contraction_storage) .AND. ALLOCATED(tensor_3%contraction_storage)) THEN
               CALL dbt_change_pgrid_2d(tensor_2, tensor_2%pgrid%mp_comm_2d, pdims=pdims_2d_opt, &
                                        nsplit=split_opt_avg%ngroup, dimsplit=split_opt_avg%split_rowcol, &
                                        pgrid_changed=pgrid_changed, &
                                        unit_nr=unit_nr_prv)
               IF (pgrid_changed) pgrid_changed_any = .TRUE.
               CALL dbt_change_pgrid_2d(tensor_3, tensor_3%pgrid%mp_comm_2d, pdims=pdims_2d_opt, &
                                        nsplit=split_opt_avg%ngroup, dimsplit=split_opt_avg%split_rowcol, &
                                        pgrid_changed=pgrid_changed, &
                                        unit_nr=unit_nr_prv)
               IF (pgrid_changed) pgrid_changed_any = .TRUE.
            END IF
            IF (pgrid_changed_any) THEN
               IF (tensor_1%matrix_rep%do_batched == 3) THEN
                  CALL dbt_tas_batched_mm_complete(tensor_1%matrix_rep)
               END IF
            END IF
         END SELECT
         CALL dbt_tas_release_info(split_opt_avg)
      END IF

      IF ((.NOT. tensors_remapped) .AND. do_batched) THEN
         ! freeze TAS process grids if tensor grids were optimized
         CALL dbt_tas_set_batched_state(tensor_1%matrix_rep, opt_grid=.TRUE.)
         CALL dbt_tas_set_batched_state(tensor_2%matrix_rep, opt_grid=.TRUE.)
         CALL dbt_tas_set_batched_state(tensor_3%matrix_rep, opt_grid=.TRUE.)
      END IF

      CALL dbt_tas_release_info(split_opt)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief align tensor index with data
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE align_tensor(tensor_in, contract_in, notcontract_in, &
                           tensor_out, contract_out, notcontract_out, indp_in, indp_out)
      TYPE(dbt_type), INTENT(INOUT)               :: tensor_in
      INTEGER, DIMENSION(:), INTENT(IN)            :: contract_in, notcontract_in
      TYPE(dbt_type), INTENT(OUT)              :: tensor_out
      INTEGER, DIMENSION(SIZE(contract_in)), &
         INTENT(OUT)                               :: contract_out
      INTEGER, DIMENSION(SIZE(notcontract_in)), &
         INTENT(OUT)                               :: notcontract_out
      CHARACTER(LEN=1), DIMENSION(ndims_tensor(tensor_in)), INTENT(IN) :: indp_in
      CHARACTER(LEN=1), DIMENSION(ndims_tensor(tensor_in)), INTENT(OUT) :: indp_out
      INTEGER, DIMENSION(ndims_tensor(tensor_in)) :: align

      CALL dbt_align_index(tensor_in, tensor_out, order=align)
      contract_out = align(contract_in)
      notcontract_out = align(notcontract_in)
      indp_out(align) = indp_in

   END SUBROUTINE

! **************************************************************************************************
!> \brief Prepare tensor for contraction: redistribute to a 2d format which can be contracted by
!>        matrix multiplication. This routine reshapes the two largest of the three tensors.
!>        Redistribution is avoided if tensors already in a consistent layout.
!> \param ind1_free indices of tensor 1 that are "free" (not linked to any index of tensor 2)
!> \param ind1_linked indices of tensor 1 that are linked to indices of tensor 2
!>                    1:1 correspondence with ind1_linked
!> \param trans1 transpose flag of matrix rep. tensor 1
!> \param trans2 transpose flag of matrix rep. tensor 2
!> \param new1 whether a new tensor 1 was created
!> \param new2 whether a new tensor 2 was created
!> \param nodata1 don't copy data of tensor 1
!> \param nodata2 don't copy data of tensor 2
!> \param move_data_1 memory optimization: transfer data s.t. tensor1 may be empty on return
!> \param move_data_2 memory optimization: transfer data s.t. tensor2 may be empty on return
!> \param optimize_dist experimental: optimize distribution
!> \param unit_nr output unit
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE reshape_mm_compatible(tensor1, tensor2, tensor1_out, tensor2_out, ind1_free, ind1_linked, &
                                    ind2_free, ind2_linked, trans1, trans2, new1, new2, ref_tensor, &
                                    nodata1, nodata2, move_data_1, &
                                    move_data_2, optimize_dist, unit_nr)
      TYPE(dbt_type), TARGET, INTENT(INOUT)   :: tensor1
      TYPE(dbt_type), TARGET, INTENT(INOUT)   :: tensor2
      TYPE(dbt_type), POINTER, INTENT(OUT)    :: tensor1_out, tensor2_out
      INTEGER, DIMENSION(:), INTENT(IN)           :: ind1_free, ind2_free
      INTEGER, DIMENSION(:), INTENT(IN)           :: ind1_linked, ind2_linked
      LOGICAL, INTENT(OUT)                        :: trans1, trans2
      LOGICAL, INTENT(OUT)                        :: new1, new2
      INTEGER, INTENT(OUT) :: ref_tensor
      LOGICAL, INTENT(IN), OPTIONAL               :: nodata1, nodata2
      LOGICAL, INTENT(INOUT), OPTIONAL            :: move_data_1, move_data_2
      LOGICAL, INTENT(IN), OPTIONAL               :: optimize_dist
      INTEGER, INTENT(IN), OPTIONAL               :: unit_nr
      INTEGER                                     :: compat1, compat1_old, compat2, compat2_old, &
                                                     comm_2d, unit_nr_prv
      TYPE(array_list)                            :: dist_list
      INTEGER, DIMENSION(:), ALLOCATABLE          :: mp_dims
      TYPE(dbt_distribution_type)             :: dist_in
      INTEGER(KIND=int_8)                         :: nblkrows, nblkcols
      LOGICAL                                     :: optimize_dist_prv
      INTEGER, DIMENSION(ndims_tensor(tensor1)) :: dims1
      INTEGER, DIMENSION(ndims_tensor(tensor2)) :: dims2

      NULLIFY (tensor1_out, tensor2_out)

      unit_nr_prv = prep_output_unit(unit_nr)

      CALL blk_dims_tensor(tensor1, dims1)
      CALL blk_dims_tensor(tensor2, dims2)

      IF (PRODUCT(int(dims1, int_8)) .GE. PRODUCT(int(dims2, int_8))) THEN
         ref_tensor = 1
      ELSE
         ref_tensor = 2
      END IF

      IF (PRESENT(optimize_dist)) THEN
         optimize_dist_prv = optimize_dist
      ELSE
         optimize_dist_prv = .FALSE.
      END IF

      compat1 = compat_map(tensor1%nd_index, ind1_linked)
      compat2 = compat_map(tensor2%nd_index, ind2_linked)
      compat1_old = compat1
      compat2_old = compat2

      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, '(T2,A,1X,A,A,1X)', advance='no') "compatibility of", TRIM(tensor1%name), ":"
         SELECT CASE (compat1)
         CASE (0)
            WRITE (unit_nr_prv, '(A)') "Not compatible"
         CASE (1)
            WRITE (unit_nr_prv, '(A)') "Normal"
         CASE (2)
            WRITE (unit_nr_prv, '(A)') "Transposed"
         END SELECT
         WRITE (unit_nr_prv, '(T2,A,1X,A,A,1X)', advance='no') "compatibility of", TRIM(tensor2%name), ":"
         SELECT CASE (compat2)
         CASE (0)
            WRITE (unit_nr_prv, '(A)') "Not compatible"
         CASE (1)
            WRITE (unit_nr_prv, '(A)') "Normal"
         CASE (2)
            WRITE (unit_nr_prv, '(A)') "Transposed"
         END SELECT
      END IF

      new1 = .FALSE.
      new2 = .FALSE.

      IF (compat1 == 0 .OR. optimize_dist_prv) THEN
         new1 = .TRUE.
      END IF

      IF (compat2 == 0 .OR. optimize_dist_prv) THEN
         new2 = .TRUE.
      END IF

      IF (ref_tensor == 1) THEN ! tensor 1 is reference and tensor 2 is reshaped compatible with tensor 1
         IF (compat1 == 0 .OR. optimize_dist_prv) THEN ! tensor 1 is not contraction compatible --> reshape
            IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A)') "Redistribution of", TRIM(tensor1%name)
            nblkrows = PRODUCT(INT(dims1(ind1_linked), KIND=int_8))
            nblkcols = PRODUCT(INT(dims1(ind1_free), KIND=int_8))
            comm_2d = dbt_tas_mp_comm(tensor1%pgrid%mp_comm_2d, nblkrows, nblkcols)
            ALLOCATE (tensor1_out)
            CALL dbt_remap(tensor1, ind1_linked, ind1_free, tensor1_out, comm_2d=comm_2d, &
                           nodata=nodata1, move_data=move_data_1)
            CALL mp_comm_free(comm_2d)
            compat1 = 1
         ELSE
            IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A)') "No redistribution of", TRIM(tensor1%name)
            tensor1_out => tensor1
         END IF
         IF (compat2 == 0 .OR. optimize_dist_prv) THEN ! tensor 2 is not contraction compatible --> reshape
            IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A,1X,A,1X,A)') "Redistribution of", &
               TRIM(tensor2%name), "compatible with", TRIM(tensor1%name)
            dist_in = dbt_distribution(tensor1_out)
            dist_list = array_sublist(dist_in%nd_dist, ind1_linked)
            IF (compat1 == 1) THEN ! linked index is first 2d dimension
               ! get distribution of linked index, tensor 2 must adopt this distribution
               ! get grid dimensions of linked index
               ALLOCATE (mp_dims(ndims_mapping_row(dist_in%pgrid%nd_index_grid)))
               CALL dbt_get_mapping_info(dist_in%pgrid%nd_index_grid, dims1_2d=mp_dims)
               ALLOCATE (tensor2_out)
               CALL dbt_remap(tensor2, ind2_linked, ind2_free, tensor2_out, comm_2d=dist_in%pgrid%mp_comm_2d, &
                              dist1=dist_list, mp_dims_1=mp_dims, nodata=nodata2, move_data=move_data_2)
            ELSEIF (compat1 == 2) THEN ! linked index is second 2d dimension
               ! get distribution of linked index, tensor 2 must adopt this distribution
               ! get grid dimensions of linked index
               ALLOCATE (mp_dims(ndims_mapping_column(dist_in%pgrid%nd_index_grid)))
               CALL dbt_get_mapping_info(dist_in%pgrid%nd_index_grid, dims2_2d=mp_dims)
               ALLOCATE (tensor2_out)
               CALL dbt_remap(tensor2, ind2_free, ind2_linked, tensor2_out, comm_2d=dist_in%pgrid%mp_comm_2d, &
                              dist2=dist_list, mp_dims_2=mp_dims, nodata=nodata2, move_data=move_data_2)
            ELSE
               CPABORT("should not happen")
            END IF
            compat2 = compat1
         ELSE
            IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A)') "No redistribution of", TRIM(tensor2%name)
            tensor2_out => tensor2
         END IF
      ELSE ! tensor 2 is reference and tensor 1 is reshaped compatible with tensor 2
         IF (compat2 == 0 .OR. optimize_dist_prv) THEN ! tensor 2 is not contraction compatible --> reshape
            IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A)') "Redistribution of", TRIM(tensor2%name)
            nblkrows = PRODUCT(INT(dims2(ind2_linked), KIND=int_8))
            nblkcols = PRODUCT(INT(dims2(ind2_free), KIND=int_8))
            comm_2d = dbt_tas_mp_comm(tensor2%pgrid%mp_comm_2d, nblkrows, nblkcols)
            ALLOCATE (tensor2_out)
            CALL dbt_remap(tensor2, ind2_linked, ind2_free, tensor2_out, nodata=nodata2, move_data=move_data_2)
            CALL mp_comm_free(comm_2d)
            compat2 = 1
         ELSE
            IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A)') "No redistribution of", TRIM(tensor2%name)
            tensor2_out => tensor2
         END IF
         IF (compat1 == 0 .OR. optimize_dist_prv) THEN ! tensor 1 is not contraction compatible --> reshape
            IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A,1X,A,1X,A)') "Redistribution of", TRIM(tensor1%name), &
               "compatible with", TRIM(tensor2%name)
            dist_in = dbt_distribution(tensor2_out)
            dist_list = array_sublist(dist_in%nd_dist, ind2_linked)
            IF (compat2 == 1) THEN
               ALLOCATE (mp_dims(ndims_mapping_row(dist_in%pgrid%nd_index_grid)))
               CALL dbt_get_mapping_info(dist_in%pgrid%nd_index_grid, dims1_2d=mp_dims)
               ALLOCATE (tensor1_out)
               CALL dbt_remap(tensor1, ind1_linked, ind1_free, tensor1_out, comm_2d=dist_in%pgrid%mp_comm_2d, &
                              dist1=dist_list, mp_dims_1=mp_dims, nodata=nodata1, move_data=move_data_1)
            ELSEIF (compat2 == 2) THEN
               ALLOCATE (mp_dims(ndims_mapping_column(dist_in%pgrid%nd_index_grid)))
               CALL dbt_get_mapping_info(dist_in%pgrid%nd_index_grid, dims2_2d=mp_dims)
               ALLOCATE (tensor1_out)
               CALL dbt_remap(tensor1, ind1_free, ind1_linked, tensor1_out, comm_2d=dist_in%pgrid%mp_comm_2d, &
                              dist2=dist_list, mp_dims_2=mp_dims, nodata=nodata1, move_data=move_data_1)
            ELSE
               CPABORT("should not happen")
            END IF
            compat1 = compat2
         ELSE
            IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A)') "No redistribution of", TRIM(tensor1%name)
            tensor1_out => tensor1
         END IF
      END IF

      SELECT CASE (compat1)
      CASE (1)
         trans1 = .FALSE.
      CASE (2)
         trans1 = .TRUE.
      CASE DEFAULT
         CPABORT("should not happen")
      END SELECT

      SELECT CASE (compat2)
      CASE (1)
         trans2 = .FALSE.
      CASE (2)
         trans2 = .TRUE.
      CASE DEFAULT
         CPABORT("should not happen")
      END SELECT

      IF (unit_nr_prv > 0) THEN
         IF (compat1 .NE. compat1_old) THEN
            WRITE (unit_nr_prv, '(T2,A,1X,A,A,1X)', advance='no') "compatibility of", TRIM(tensor1_out%name), ":"
            SELECT CASE (compat1)
            CASE (0)
               WRITE (unit_nr_prv, '(A)') "Not compatible"
            CASE (1)
               WRITE (unit_nr_prv, '(A)') "Normal"
            CASE (2)
               WRITE (unit_nr_prv, '(A)') "Transposed"
            END SELECT
         END IF
         IF (compat2 .NE. compat2_old) THEN
            WRITE (unit_nr_prv, '(T2,A,1X,A,A,1X)', advance='no') "compatibility of", TRIM(tensor2_out%name), ":"
            SELECT CASE (compat2)
            CASE (0)
               WRITE (unit_nr_prv, '(A)') "Not compatible"
            CASE (1)
               WRITE (unit_nr_prv, '(A)') "Normal"
            CASE (2)
               WRITE (unit_nr_prv, '(A)') "Transposed"
            END SELECT
         END IF
      END IF

      IF (new1 .AND. PRESENT(move_data_1)) move_data_1 = .TRUE.
      IF (new2 .AND. PRESENT(move_data_2)) move_data_2 = .TRUE.

   END SUBROUTINE

! **************************************************************************************************
!> \brief Prepare tensor for contraction: redistribute to a 2d format which can be contracted by
!>        matrix multiplication. This routine reshapes the smallest of the three tensors.
!> \param ind1 index that should be mapped to first matrix dimension
!> \param ind2 index that should be mapped to second matrix dimension
!> \param trans transpose flag of matrix rep.
!> \param new whether a new tensor was created for tensor_out
!> \param nodata don't copy tensor data
!> \param move_data memory optimization: transfer data s.t. tensor_in may be empty on return
!> \param unit_nr output unit
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE reshape_mm_small(tensor_in, ind1, ind2, tensor_out, trans, new, nodata, move_data, unit_nr)
      TYPE(dbt_type), TARGET, INTENT(INOUT)   :: tensor_in
      INTEGER, DIMENSION(:), INTENT(IN)           :: ind1, ind2
      TYPE(dbt_type), POINTER, INTENT(OUT)    :: tensor_out
      LOGICAL, INTENT(OUT)                        :: trans
      LOGICAL, INTENT(OUT)                        :: new
      LOGICAL, INTENT(IN), OPTIONAL               :: nodata, move_data
      INTEGER, INTENT(IN), OPTIONAL               :: unit_nr
      INTEGER                                     :: compat1, compat2, compat1_old, compat2_old, unit_nr_prv
      LOGICAL                                     :: nodata_prv

      NULLIFY (tensor_out)
      IF (PRESENT(nodata)) THEN
         nodata_prv = nodata
      ELSE
         nodata_prv = .FALSE.
      END IF

      unit_nr_prv = prep_output_unit(unit_nr)

      new = .FALSE.
      compat1 = compat_map(tensor_in%nd_index, ind1)
      compat2 = compat_map(tensor_in%nd_index, ind2)
      compat1_old = compat1; compat2_old = compat2
      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, '(T2,A,1X,A,A,1X)', advance='no') "compatibility of", TRIM(tensor_in%name), ":"
         IF (compat1 == 1 .AND. compat2 == 2) THEN
            WRITE (unit_nr_prv, '(A)') "Normal"
         ELSEIF (compat1 == 2 .AND. compat2 == 1) THEN
            WRITE (unit_nr_prv, '(A)') "Transposed"
         ELSE
            WRITE (unit_nr_prv, '(A)') "Not compatible"
         END IF
      END IF
      IF (compat1 == 0 .or. compat2 == 0) THEN ! index mapping not compatible with contract index

         IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A)') "Redistribution of", TRIM(tensor_in%name)

         ALLOCATE (tensor_out)
         CALL dbt_remap(tensor_in, ind1, ind2, tensor_out, nodata=nodata, move_data=move_data)
         CALL dbt_copy_contraction_storage(tensor_in, tensor_out)
         compat1 = 1
         compat2 = 2
         new = .TRUE.
      ELSE
         IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T2,A,1X,A)') "No redistribution of", TRIM(tensor_in%name)
         tensor_out => tensor_in
      END IF

      IF (compat1 == 1 .AND. compat2 == 2) THEN
         trans = .FALSE.
      ELSEIF (compat1 == 2 .AND. compat2 == 1) THEN
         trans = .TRUE.
      ELSE
         CPABORT("this should not happen")
      END IF

      IF (unit_nr_prv > 0) THEN
         IF (compat1_old .NE. compat1 .OR. compat2_old .NE. compat2) THEN
            WRITE (unit_nr_prv, '(T2,A,1X,A,A,1X)', advance='no') "compatibility of", TRIM(tensor_out%name), ":"
            IF (compat1 == 1 .AND. compat2 == 2) THEN
               WRITE (unit_nr_prv, '(A)') "Normal"
            ELSEIF (compat1 == 2 .AND. compat2 == 1) THEN
               WRITE (unit_nr_prv, '(A)') "Transposed"
            ELSE
               WRITE (unit_nr_prv, '(A)') "Not compatible"
            END IF
         END IF
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief update contraction storage that keeps track of process grids during a batched contraction
!>        and decide if tensor process grid needs to be optimized
!> \param split_opt optimized TAS process grid
!> \param split current TAS process grid
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION update_contraction_storage(storage, split_opt, split) RESULT(do_change_pgrid)
      TYPE(dbt_contraction_storage), INTENT(INOUT) :: storage
      TYPE(dbt_tas_split_info), INTENT(IN)           :: split_opt
      TYPE(dbt_tas_split_info), INTENT(IN)           :: split
      INTEGER, DIMENSION(2) :: pdims_opt, coor, pdims, pdims_sub
      LOGICAL, DIMENSION(2) :: periods
      LOGICAL, DIMENSION(2) :: do_change_pgrid
      REAL(kind=dp) :: change_criterion, pdims_ratio
      INTEGER :: nsplit_opt, nsplit

      CPASSERT(ALLOCATED(split_opt%ngroup_opt))
      nsplit_opt = split_opt%ngroup_opt
      nsplit = split%ngroup

      CALL mp_environ(split_opt%mp_comm, 2, pdims_opt, coor, periods)
      CALL mp_environ(split%mp_comm, 2, pdims, coor, periods)

      storage%ibatch = storage%ibatch + 1

      storage%nsplit_avg = (storage%nsplit_avg*REAL(storage%ibatch - 1, dp) + REAL(nsplit_opt, dp)) &
                           /REAL(storage%ibatch, dp)

      SELECT CASE (split_opt%split_rowcol)
      CASE (rowsplit)
         pdims_ratio = REAL(pdims(1), dp)/pdims(2)
      CASE (colsplit)
         pdims_ratio = REAL(pdims(2), dp)/pdims(1)
      END SELECT

      do_change_pgrid(:) = .FALSE.

      ! check for process grid dimensions
      CALL mp_environ(split%mp_comm_group, 2, pdims_sub, coor, periods)
      change_criterion = MAXVAL(REAL(pdims_sub, dp))/MINVAL(pdims_sub)
      IF (change_criterion > default_pdims_accept_ratio**2) do_change_pgrid(1) = .TRUE.

      ! check for split factor
      change_criterion = MAX(REAL(nsplit, dp)/storage%nsplit_avg, REAL(storage%nsplit_avg, dp)/nsplit)
      IF (change_criterion > default_nsplit_accept_ratio) do_change_pgrid(2) = .TRUE.

   END FUNCTION

! **************************************************************************************************
!> \brief Check if 2d index is compatible with tensor index
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION compat_map(nd_index, compat_ind)
      TYPE(nd_to_2d_mapping), INTENT(IN) :: nd_index
      INTEGER, DIMENSION(:), INTENT(IN)  :: compat_ind
      INTEGER, DIMENSION(ndims_mapping_row(nd_index)) :: map1
      INTEGER, DIMENSION(ndims_mapping_column(nd_index)) :: map2
      INTEGER                            :: compat_map

      CALL dbt_get_mapping_info(nd_index, map1_2d=map1, map2_2d=map2)

      compat_map = 0
      IF (array_eq_i(map1, compat_ind)) THEN
         compat_map = 1
      ELSEIF (array_eq_i(map2, compat_ind)) THEN
         compat_map = 2
      END IF

   END FUNCTION

! **************************************************************************************************
!> \brief
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE index_linked_sort(ind_ref, ind_dep)
      INTEGER, DIMENSION(:), INTENT(INOUT) :: ind_ref, ind_dep
      INTEGER, DIMENSION(SIZE(ind_ref))    :: sort_indices

      CALL sort(ind_ref, SIZE(ind_ref), sort_indices)
      ind_dep(:) = ind_dep(sort_indices)

   END SUBROUTINE

! **************************************************************************************************
!> \brief
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION opt_pgrid(tensor, tas_split_info)
      TYPE(dbt_type), INTENT(IN) :: tensor
      TYPE(dbt_tas_split_info), INTENT(IN) :: tas_split_info
      INTEGER, DIMENSION(ndims_matrix_row(tensor)) :: map1
      INTEGER, DIMENSION(ndims_matrix_column(tensor)) :: map2
      TYPE(dbt_pgrid_type) :: opt_pgrid
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: dims

      CALL dbt_get_mapping_info(tensor%pgrid%nd_index_grid, map1_2d=map1, map2_2d=map2)
      CALL blk_dims_tensor(tensor, dims)
      opt_pgrid = dbt_nd_mp_comm(tas_split_info%mp_comm, map1, map2, tdims=dims)

      ALLOCATE (opt_pgrid%tas_split_info, SOURCE=tas_split_info)
      CALL dbt_tas_info_hold(opt_pgrid%tas_split_info)
   END FUNCTION

! **************************************************************************************************
!> \brief Copy tensor to tensor with modified index mapping
!> \param map1_2d new index mapping
!> \param map2_2d new index mapping
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_remap(tensor_in, map1_2d, map2_2d, tensor_out, comm_2d, dist1, dist2, &
                        mp_dims_1, mp_dims_2, name, nodata, move_data)
      TYPE(dbt_type), INTENT(INOUT)      :: tensor_in
      INTEGER, DIMENSION(:), INTENT(IN)      :: map1_2d, map2_2d
      TYPE(dbt_type), INTENT(OUT)        :: tensor_out
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: name
      LOGICAL, INTENT(IN), OPTIONAL          :: nodata, move_data
      INTEGER, INTENT(IN), OPTIONAL          :: comm_2d
      TYPE(array_list), INTENT(IN), OPTIONAL :: dist1, dist2
      INTEGER, DIMENSION(SIZE(map1_2d)), OPTIONAL :: mp_dims_1
      INTEGER, DIMENSION(SIZE(map2_2d)), OPTIONAL :: mp_dims_2
      CHARACTER(len=default_string_length)   :: name_tmp
      INTEGER, DIMENSION(:), ALLOCATABLE     :: ${varlist("blk_sizes")}$, &
                                                ${varlist("nd_dist")}$
      TYPE(dbt_distribution_type)        :: dist
      INTEGER                                :: comm_2d_prv, handle, i
      INTEGER, DIMENSION(ndims_tensor(tensor_in)) :: pdims, myploc
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_remap'
      LOGICAL                               :: nodata_prv
      TYPE(dbt_pgrid_type)              :: comm_nd

      CALL timeset(routineN, handle)

      IF (PRESENT(name)) THEN
         name_tmp = name
      ELSE
         name_tmp = tensor_in%name
      END IF
      IF (PRESENT(dist1)) THEN
         CPASSERT(PRESENT(mp_dims_1))
      END IF

      IF (PRESENT(dist2)) THEN
         CPASSERT(PRESENT(mp_dims_2))
      END IF

      IF (PRESENT(comm_2d)) THEN
         comm_2d_prv = comm_2d
      ELSE
         comm_2d_prv = tensor_in%pgrid%mp_comm_2d
      END IF

      comm_nd = dbt_nd_mp_comm(comm_2d_prv, map1_2d, map2_2d, dims1_nd=mp_dims_1, dims2_nd=mp_dims_2)
      CALL mp_environ_pgrid(comm_nd, pdims, myploc)

      #:for ndim in ndims
         IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
            CALL get_arrays(tensor_in%blk_sizes, ${varlist("blk_sizes", nmax=ndim)}$)
         END IF
      #:endfor

      #:for ndim in ndims
         IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
            #:for idim in range(1, ndim+1)
               IF (PRESENT(dist1)) THEN
                  IF (ANY(map1_2d == ${idim}$)) THEN
                     i = MINLOC(map1_2d, dim=1, mask=map1_2d == ${idim}$) ! i is location of idim in map1_2d
                     CALL get_ith_array(dist1, i, nd_dist_${idim}$)
                  END IF
               END IF

               IF (PRESENT(dist2)) THEN
                  IF (ANY(map2_2d == ${idim}$)) THEN
                     i = MINLOC(map2_2d, dim=1, mask=map2_2d == ${idim}$) ! i is location of idim in map2_2d
                     CALL get_ith_array(dist2, i, nd_dist_${idim}$)
                  END IF
               END IF

               IF (.NOT. ALLOCATED(nd_dist_${idim}$)) THEN
                  ALLOCATE (nd_dist_${idim}$ (SIZE(blk_sizes_${idim}$)))
                  CALL dbt_default_distvec(SIZE(blk_sizes_${idim}$), pdims(${idim}$), blk_sizes_${idim}$, nd_dist_${idim}$)
               END IF
            #:endfor
            CALL dbt_distribution_new_expert(dist, comm_nd, map1_2d, map2_2d, &
                                             ${varlist("nd_dist", nmax=ndim)}$, own_comm=.TRUE.)
         END IF
      #:endfor

      #:for ndim in ndims
         IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
            CALL dbt_create(tensor_out, name_tmp, dist, map1_2d, map2_2d, &
                            ${varlist("blk_sizes", nmax=ndim)}$)
         END IF
      #:endfor

      IF (PRESENT(nodata)) THEN
         nodata_prv = nodata
      ELSE
         nodata_prv = .FALSE.
      END IF

      IF (.NOT. nodata_prv) CALL dbt_copy_expert(tensor_in, tensor_out, move_data=move_data)
      CALL dbt_distribution_destroy(dist)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Align index with data
!> \param order permutation resulting from alignment
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_align_index(tensor_in, tensor_out, order)
      TYPE(dbt_type), INTENT(INOUT)               :: tensor_in
      TYPE(dbt_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(ndims_matrix_row(tensor_in)) :: map1_2d
      INTEGER, DIMENSION(ndims_matrix_column(tensor_in)) :: map2_2d
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(OUT), OPTIONAL                        :: order
      INTEGER, DIMENSION(ndims_tensor(tensor_in))     :: order_prv
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_align_index'
      INTEGER                                         :: handle

      CALL timeset(routineN, handle)

      CALL dbt_get_mapping_info(tensor_in%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      order_prv = dbt_inverse_order([map1_2d, map2_2d])
      CALL dbt_permute_index(tensor_in, tensor_out, order=order_prv)

      IF (PRESENT(order)) order = order_prv

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create new tensor by reordering index, data is copied exactly (shallow copy)
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_permute_index(tensor_in, tensor_out, order)
      TYPE(dbt_type), INTENT(INOUT)                  :: tensor_in
      TYPE(dbt_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(IN)                                   :: order

      TYPE(nd_to_2d_mapping)                          :: nd_index_blk_rs, nd_index_rs
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_permute_index'
      INTEGER                                         :: handle
      INTEGER                                         :: ndims

      CALL timeset(routineN, handle)

      ndims = ndims_tensor(tensor_in)

      CALL permute_index(tensor_in%nd_index, nd_index_rs, order)
      CALL permute_index(tensor_in%nd_index_blk, nd_index_blk_rs, order)
      CALL permute_index(tensor_in%pgrid%nd_index_grid, tensor_out%pgrid%nd_index_grid, order)

      tensor_out%matrix_rep => tensor_in%matrix_rep
      tensor_out%owns_matrix = .FALSE.

      tensor_out%nd_index = nd_index_rs
      tensor_out%nd_index_blk = nd_index_blk_rs
      tensor_out%pgrid%mp_comm_2d = tensor_in%pgrid%mp_comm_2d
      IF (ALLOCATED(tensor_in%pgrid%tas_split_info)) THEN
         ALLOCATE (tensor_out%pgrid%tas_split_info, SOURCE=tensor_in%pgrid%tas_split_info)
      END IF
      tensor_out%refcount => tensor_in%refcount
      CALL dbt_hold(tensor_out)

      CALL reorder_arrays(tensor_in%blk_sizes, tensor_out%blk_sizes, order)
      CALL reorder_arrays(tensor_in%blk_offsets, tensor_out%blk_offsets, order)
      CALL reorder_arrays(tensor_in%nd_dist, tensor_out%nd_dist, order)
      CALL reorder_arrays(tensor_in%blks_local, tensor_out%blks_local, order)
      ALLOCATE (tensor_out%nblks_local(ndims))
      ALLOCATE (tensor_out%nfull_local(ndims))
      tensor_out%nblks_local(order) = tensor_in%nblks_local(:)
      tensor_out%nfull_local(order) = tensor_in%nfull_local(:)
      tensor_out%name = tensor_in%name
      tensor_out%valid = .TRUE.

      IF (ALLOCATED(tensor_in%contraction_storage)) THEN
         ALLOCATE (tensor_out%contraction_storage, SOURCE=tensor_in%contraction_storage)
         CALL destroy_array_list(tensor_out%contraction_storage%batch_ranges)
         CALL reorder_arrays(tensor_in%contraction_storage%batch_ranges, tensor_out%contraction_storage%batch_ranges, order)
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief get indices of non-zero tensor blocks for contraction result without actually
!>        performing contraction.
!>        this is an estimate based on block norm multiplication.
!>        See documentation of dbt_contract.
!> \param nblks_local number of local blocks on this MPI rank
!> \param result_index indices of local non-zero tensor blocks for tensor_3
!>        only the elements result_index(:nblks_local, :) are relevant (all others are set to 0)
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_contract_index(alpha, tensor_1, tensor_2, beta, tensor_3, &
                                 contract_1, notcontract_1, &
                                 contract_2, notcontract_2, &
                                 map_1, map_2, &
                                 bounds_1, bounds_2, bounds_3, &
                                 filter_eps, &
                                 nblks_local, result_index)
      REAL(dp), INTENT(IN)            :: alpha
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_1
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_2
      REAL(dp), INTENT(IN)            :: beta
      INTEGER, DIMENSION(:), INTENT(IN)              :: contract_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: contract_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: map_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: map_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: notcontract_1
      INTEGER, DIMENSION(:), INTENT(IN)              :: notcontract_2
      TYPE(dbt_type), INTENT(INOUT), TARGET      :: tensor_3
      INTEGER, DIMENSION(2, SIZE(contract_1)), &
         INTENT(IN), OPTIONAL                        :: bounds_1
      INTEGER, DIMENSION(2, SIZE(notcontract_1)), &
         INTENT(IN), OPTIONAL                        :: bounds_2
      INTEGER, DIMENSION(2, SIZE(notcontract_2)), &
         INTENT(IN), OPTIONAL                        :: bounds_3
      REAL(KIND=dp), INTENT(IN), OPTIONAL        :: filter_eps
      INTEGER, INTENT(OUT)                           :: nblks_local
      INTEGER, DIMENSION(dbt_max_nblks_local(tensor_3), ndims_tensor(tensor_3)), &
         INTENT(OUT)                                 :: result_index

      CALL dbt_contract_expert(alpha, tensor_1, tensor_2, beta, tensor_3, &
                               contract_1, notcontract_1, &
                               contract_2, notcontract_2, &
                               map_1, map_2, &
                               bounds_1=bounds_1, &
                               bounds_2=bounds_2, &
                               bounds_3=bounds_3, &
                               filter_eps=filter_eps, &
                               nblks_local=nblks_local, &
                               result_index=result_index)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Map contraction bounds to bounds referring to tensor indices
!>        see dbt_contract for docu of dummy arguments
!> \param bounds_t1 bounds mapped to tensor_1
!> \param bounds_t2 bounds mapped to tensor_2
!> \param do_crop_1 whether tensor 1 should be cropped
!> \param do_crop_2 whether tensor 2 should be cropped
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_map_bounds_to_tensors(tensor_1, tensor_2, &
                                        contract_1, notcontract_1, &
                                        contract_2, notcontract_2, &
                                        bounds_t1, bounds_t2, &
                                        bounds_1, bounds_2, bounds_3, &
                                        do_crop_1, do_crop_2)

      TYPE(dbt_type), INTENT(IN)      :: tensor_1, tensor_2
      INTEGER, DIMENSION(:), INTENT(IN)   :: contract_1, contract_2, &
                                             notcontract_1, notcontract_2
      INTEGER, DIMENSION(2, ndims_tensor(tensor_1)), &
         INTENT(OUT)                                 :: bounds_t1
      INTEGER, DIMENSION(2, ndims_tensor(tensor_2)), &
         INTENT(OUT)                                 :: bounds_t2
      INTEGER, DIMENSION(2, SIZE(contract_1)), &
         INTENT(IN), OPTIONAL                        :: bounds_1
      INTEGER, DIMENSION(2, SIZE(notcontract_1)), &
         INTENT(IN), OPTIONAL                        :: bounds_2
      INTEGER, DIMENSION(2, SIZE(notcontract_2)), &
         INTENT(IN), OPTIONAL                        :: bounds_3
      LOGICAL, INTENT(OUT), OPTIONAL                 :: do_crop_1, do_crop_2
      LOGICAL, DIMENSION(2)                          :: do_crop

      do_crop = .FALSE.

      bounds_t1(1, :) = 1
      CALL dbt_get_info(tensor_1, nfull_total=bounds_t1(2, :))

      bounds_t2(1, :) = 1
      CALL dbt_get_info(tensor_2, nfull_total=bounds_t2(2, :))

      IF (PRESENT(bounds_1)) THEN
         bounds_t1(:, contract_1) = bounds_1
         do_crop(1) = .TRUE.
         bounds_t2(:, contract_2) = bounds_1
         do_crop(2) = .TRUE.
      END IF

      IF (PRESENT(bounds_2)) THEN
         bounds_t1(:, notcontract_1) = bounds_2
         do_crop(1) = .TRUE.
      END IF

      IF (PRESENT(bounds_3)) THEN
         bounds_t2(:, notcontract_2) = bounds_3
         do_crop(2) = .TRUE.
      END IF

      IF (PRESENT(do_crop_1)) do_crop_1 = do_crop(1)
      IF (PRESENT(do_crop_2)) do_crop_2 = do_crop(2)

   END SUBROUTINE

! **************************************************************************************************
!> \brief print tensor contraction indices in a human readable way
!> \param indchar1 characters printed for index of tensor 1
!> \param indchar2 characters printed for index of tensor 2
!> \param indchar3 characters printed for index of tensor 3
!> \param unit_nr output unit
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_print_contraction_index(tensor_1, indchar1, tensor_2, indchar2, tensor_3, indchar3, unit_nr)
      TYPE(dbt_type), INTENT(IN) :: tensor_1, tensor_2, tensor_3
      CHARACTER(LEN=1), DIMENSION(ndims_tensor(tensor_1)), INTENT(IN) :: indchar1
      CHARACTER(LEN=1), DIMENSION(ndims_tensor(tensor_2)), INTENT(IN) :: indchar2
      CHARACTER(LEN=1), DIMENSION(ndims_tensor(tensor_3)), INTENT(IN) :: indchar3
      INTEGER, INTENT(IN) :: unit_nr
      INTEGER, DIMENSION(ndims_matrix_row(tensor_1)) :: map11
      INTEGER, DIMENSION(ndims_matrix_column(tensor_1)) :: map12
      INTEGER, DIMENSION(ndims_matrix_row(tensor_2)) :: map21
      INTEGER, DIMENSION(ndims_matrix_column(tensor_2)) :: map22
      INTEGER, DIMENSION(ndims_matrix_row(tensor_3)) :: map31
      INTEGER, DIMENSION(ndims_matrix_column(tensor_3)) :: map32
      INTEGER :: ichar1, ichar2, ichar3, unit_nr_prv

      unit_nr_prv = prep_output_unit(unit_nr)

      IF (unit_nr_prv /= 0) THEN
         CALL dbt_get_mapping_info(tensor_1%nd_index_blk, map1_2d=map11, map2_2d=map12)
         CALL dbt_get_mapping_info(tensor_2%nd_index_blk, map1_2d=map21, map2_2d=map22)
         CALL dbt_get_mapping_info(tensor_3%nd_index_blk, map1_2d=map31, map2_2d=map32)
      END IF

      IF (unit_nr_prv > 0) THEN
         WRITE (unit_nr_prv, '(T2,A)') "INDEX INFO"
         WRITE (unit_nr_prv, '(T15,A)', advance='no') "tensor index: ("
         DO ichar1 = 1, SIZE(indchar1)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar1(ichar1)
         END DO
         WRITE (unit_nr_prv, '(A)', advance='no') ") x ("
         DO ichar2 = 1, SIZE(indchar2)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar2(ichar2)
         END DO
         WRITE (unit_nr_prv, '(A)', advance='no') ") = ("
         DO ichar3 = 1, SIZE(indchar3)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar3(ichar3)
         END DO
         WRITE (unit_nr_prv, '(A)') ")"

         WRITE (unit_nr_prv, '(T15,A)', advance='no') "matrix index: ("
         DO ichar1 = 1, SIZE(map11)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar1(map11(ichar1))
         END DO
         WRITE (unit_nr_prv, '(A1)', advance='no') "|"
         DO ichar1 = 1, SIZE(map12)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar1(map12(ichar1))
         END DO
         WRITE (unit_nr_prv, '(A)', advance='no') ") x ("
         DO ichar2 = 1, SIZE(map21)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar2(map21(ichar2))
         END DO
         WRITE (unit_nr_prv, '(A1)', advance='no') "|"
         DO ichar2 = 1, SIZE(map22)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar2(map22(ichar2))
         END DO
         WRITE (unit_nr_prv, '(A)', advance='no') ") = ("
         DO ichar3 = 1, SIZE(map31)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar3(map31(ichar3))
         END DO
         WRITE (unit_nr_prv, '(A1)', advance='no') "|"
         DO ichar3 = 1, SIZE(map32)
            WRITE (unit_nr_prv, '(A1)', advance='no') indchar3(map32(ichar3))
         END DO
         WRITE (unit_nr_prv, '(A)') ")"
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Initialize batched contraction for this tensor.
!>
!>        Explanation: A batched contraction is a contraction performed in several consecutive steps
!>        by specification of bounds in dbt_contract. This can be used to reduce memory by
!>        a large factor. The routines dbt_batched_contract_init and
!>        dbt_batched_contract_finalize should be called to define the scope of a batched
!>        contraction as this enables important optimizations (adapting communication scheme to
!>        batches and adapting process grid to multiplication algorithm). The routines
!>        dbt_batched_contract_init and dbt_batched_contract_finalize must be
!>        called before the first and after the last contraction step on all 3 tensors.
!>
!>        Requirements:
!>        - the tensors are in a compatible matrix layout (see documentation of
!>          `dbt_contract`, note 2 & 3). If they are not, process grid optimizations are
!>          disabled and a warning is issued.
!>        - within the scope of a batched contraction, it is not allowed to access or change tensor
!>          data except by calling the routines dbt_contract & dbt_copy.
!>        - the bounds affecting indices of the smallest tensor must not change in the course of a
!>          batched contraction (todo: get rid of this requirement).
!>
!>        Side effects:
!>        - the parallel layout (process grid and distribution) of all tensors may change. In order
!>          to disable the process grid optimization including this side effect, call this routine
!>          only on the smallest of the 3 tensors.
!>
!> \note
!>        Note 1: for an example of batched contraction see `examples/dbt_example.F`.
!>        (todo: the example is outdated and should be updated).
!>
!>        Note 2: it is meaningful to use this feature if the contraction consists of one batch only
!>        but if multiple contractions involving the same 3 tensors are performed
!>        (batched_contract_init and batched_contract_finalize must then be called before/after each
!>        contraction call). The process grid is then optimized after the first contraction
!>        and future contraction may profit from this optimization.
!>
!> \param batch_range_i refers to the ith tensor dimension and contains all block indices starting
!>                      a new range. The size should be the number of ranges plus one, the last
!>                      element being the block index plus one of the last block in the last range.
!>                      For internal load balancing optimizations, optionally specify the index
!>                      ranges of batched contraction.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_batched_contract_init(tensor, ${varlist("batch_range")}$)
      TYPE(dbt_type), INTENT(INOUT) :: tensor
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)        :: ${varlist("batch_range")}$
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: tdims
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: ${varlist("batch_range_prv")}$
      LOGICAL :: static_range

      CALL dbt_get_info(tensor, nblks_total=tdims)

      static_range = .TRUE.
      #:for idim in range(1, maxdim+1)
         IF (ndims_tensor(tensor) >= ${idim}$) THEN
            IF (PRESENT(batch_range_${idim}$)) THEN
               ALLOCATE (batch_range_prv_${idim}$, source=batch_range_${idim}$)
               static_range = .FALSE.
            ELSE
               ALLOCATE (batch_range_prv_${idim}$ (2))
               batch_range_prv_${idim}$ (1) = 1
               batch_range_prv_${idim}$ (2) = tdims(${idim}$) + 1
            END IF
         END IF
      #:endfor

      ALLOCATE (tensor%contraction_storage)
      tensor%contraction_storage%static = static_range
      IF (static_range) THEN
         CALL dbt_tas_batched_mm_init(tensor%matrix_rep)
      END IF
      tensor%contraction_storage%nsplit_avg = 0.0_dp
      tensor%contraction_storage%ibatch = 0

      #:for ndim in range(1, maxdim+1)
         IF (ndims_tensor(tensor) == ${ndim}$) THEN
            CALL create_array_list(tensor%contraction_storage%batch_ranges, ${ndim}$, &
                                   ${varlist("batch_range_prv", nmax=ndim)}$)
         END IF
      #:endfor

   END SUBROUTINE

! **************************************************************************************************
!> \brief finalize batched contraction. This performs all communication that has been postponed in
!>         the contraction calls.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_batched_contract_finalize(tensor, unit_nr)
      TYPE(dbt_type), INTENT(INOUT) :: tensor
      INTEGER, INTENT(IN), OPTIONAL :: unit_nr
      LOGICAL :: do_write
      INTEGER :: unit_nr_prv, handle

      CALL mp_sync(tensor%pgrid%mp_comm_2d)
      CALL timeset("dbt_total", handle)
      unit_nr_prv = prep_output_unit(unit_nr)

      do_write = .FALSE.

      IF (tensor%contraction_storage%static) THEN
         IF (tensor%matrix_rep%do_batched > 0) THEN
            IF (tensor%matrix_rep%mm_storage%batched_out) do_write = .TRUE.
         END IF
         CALL dbt_tas_batched_mm_finalize(tensor%matrix_rep)
      END IF

      IF (do_write .AND. unit_nr_prv /= 0) THEN
         IF (unit_nr_prv > 0) THEN
            WRITE (unit_nr_prv, "(T2,A)") &
               "FINALIZING BATCHED PROCESSING OF MATMUL"
         END IF
         CALL dbt_write_tensor_info(tensor, unit_nr_prv)
         CALL dbt_write_tensor_dist(tensor, unit_nr_prv)
      END IF

      CALL destroy_array_list(tensor%contraction_storage%batch_ranges)
      DEALLOCATE (tensor%contraction_storage)
      CALL mp_sync(tensor%pgrid%mp_comm_2d)
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief change the process grid of a tensor
!> \param nodata optionally don't copy the tensor data (then tensor is empty on returned)
!> \param batch_range_i refers to the ith tensor dimension and contains all block indices starting
!>                      a new range. The size should be the number of ranges plus one, the last
!>                      element being the block index plus one of the last block in the last range.
!>                      For internal load balancing optimizations, optionally specify the index
!>                      ranges of batched contraction.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_change_pgrid(tensor, pgrid, ${varlist("batch_range")}$, &
                               nodata, pgrid_changed, unit_nr)
      TYPE(dbt_type), INTENT(INOUT)                  :: tensor
      TYPE(dbt_pgrid_type), INTENT(IN)               :: pgrid
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)        :: ${varlist("batch_range")}$
      !!
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata
      LOGICAL, INTENT(OUT), OPTIONAL                     :: pgrid_changed
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_change_pgrid'
      CHARACTER(default_string_length)                   :: name
      INTEGER                                            :: handle
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ${varlist("bs")}$, &
                                                            ${varlist("dist")}$
      INTEGER, DIMENSION(ndims_tensor(tensor))           :: pcoord, pcoord_ref, pdims, pdims_ref, &
                                                            tdims
      TYPE(dbt_type)                                 :: t_tmp
      TYPE(dbt_distribution_type)                    :: dist
      INTEGER, DIMENSION(ndims_matrix_row(tensor)) :: map1
      INTEGER, &
         DIMENSION(ndims_matrix_column(tensor))    :: map2
      LOGICAL, DIMENSION(ndims_tensor(tensor))             :: mem_aware
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: nbatch
      INTEGER :: ind1, ind2, batch_size, ibatch

      IF (PRESENT(pgrid_changed)) pgrid_changed = .FALSE.
      CALL mp_environ_pgrid(pgrid, pdims, pcoord)
      CALL mp_environ_pgrid(tensor%pgrid, pdims_ref, pcoord_ref)

      IF (ALL(pdims == pdims_ref)) THEN
         IF (ALLOCATED(pgrid%tas_split_info) .AND. ALLOCATED(tensor%pgrid%tas_split_info)) THEN
            IF (pgrid%tas_split_info%ngroup == tensor%pgrid%tas_split_info%ngroup) THEN
               RETURN
            END IF
         END IF
      END IF

      CALL timeset(routineN, handle)

      #:for idim in range(1, maxdim+1)
         IF (ndims_tensor(tensor) >= ${idim}$) THEN
            mem_aware(${idim}$) = PRESENT(batch_range_${idim}$)
            IF (mem_aware(${idim}$)) nbatch(${idim}$) = SIZE(batch_range_${idim}$) - 1
         END IF
      #:endfor

      CALL dbt_get_info(tensor, nblks_total=tdims, name=name)

      #:for idim in range(1, maxdim+1)
         IF (ndims_tensor(tensor) >= ${idim}$) THEN
            ALLOCATE (bs_${idim}$ (dbt_nblks_total(tensor, ${idim}$)))
            CALL get_ith_array(tensor%blk_sizes, ${idim}$, bs_${idim}$)
            ALLOCATE (dist_${idim}$ (tdims(${idim}$)))
            dist_${idim}$ = 0
            IF (mem_aware(${idim}$)) THEN
               DO ibatch = 1, nbatch(${idim}$)
                  ind1 = batch_range_${idim}$ (ibatch)
                  ind2 = batch_range_${idim}$ (ibatch + 1) - 1
                  batch_size = ind2 - ind1 + 1
                  CALL dbt_default_distvec(batch_size, pdims(${idim}$), &
                                           bs_${idim}$ (ind1:ind2), dist_${idim}$ (ind1:ind2))
               END DO
            ELSE
               CALL dbt_default_distvec(tdims(${idim}$), pdims(${idim}$), bs_${idim}$, dist_${idim}$)
            END IF
         END IF
      #:endfor

      CALL dbt_get_mapping_info(tensor%nd_index_blk, map1_2d=map1, map2_2d=map2)
      #:for ndim in ndims
         IF (ndims_tensor(tensor) == ${ndim}$) THEN
            CALL dbt_distribution_new(dist, pgrid, ${varlist("dist", nmax=ndim)}$)
            CALL dbt_create(t_tmp, name, dist, map1, map2, ${varlist("bs", nmax=ndim)}$)
         END IF
      #:endfor
      CALL dbt_distribution_destroy(dist)

      IF (PRESENT(nodata)) THEN
         IF (.NOT. nodata) CALL dbt_copy_expert(tensor, t_tmp, move_data=.TRUE.)
      ELSE
         CALL dbt_copy_expert(tensor, t_tmp, move_data=.TRUE.)
      END IF

      CALL dbt_copy_contraction_storage(tensor, t_tmp)

      CALL dbt_destroy(tensor)
      tensor = t_tmp

      IF (PRESENT(unit_nr)) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, "(T2,A,1X,A)") "OPTIMIZED PGRID INFO FOR", TRIM(tensor%name)
            WRITE (unit_nr, "(T4,A,1X,3I6)") "process grid dimensions:", pdims
            CALL dbt_write_split_info(pgrid, unit_nr)
         END IF
      END IF

      IF (PRESENT(pgrid_changed)) pgrid_changed = .TRUE.

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief map tensor to a new 2d process grid for the matrix representation.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_change_pgrid_2d(tensor, mp_comm, pdims, nodata, nsplit, dimsplit, pgrid_changed, unit_nr)
      TYPE(dbt_type), INTENT(INOUT)                  :: tensor
      INTEGER, INTENT(IN)               :: mp_comm
      INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: pdims
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata
      INTEGER, INTENT(IN), OPTIONAL :: nsplit, dimsplit
      LOGICAL, INTENT(OUT), OPTIONAL :: pgrid_changed
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr
      INTEGER, DIMENSION(ndims_matrix_row(tensor)) :: map1
      INTEGER, DIMENSION(ndims_matrix_column(tensor)) :: map2
      INTEGER, DIMENSION(ndims_tensor(tensor)) :: dims, nbatches
      TYPE(dbt_pgrid_type) :: pgrid
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("batch_range")}$
      INTEGER, DIMENSION(:), ALLOCATABLE :: array
      INTEGER :: idim

      CALL dbt_get_mapping_info(tensor%pgrid%nd_index_grid, map1_2d=map1, map2_2d=map2)
      CALL blk_dims_tensor(tensor, dims)

      IF (ALLOCATED(tensor%contraction_storage)) THEN
         ASSOCIATE (batch_ranges => tensor%contraction_storage%batch_ranges)
            nbatches = sizes_of_arrays(tensor%contraction_storage%batch_ranges) - 1
            ! for good load balancing the process grid dimensions should be chosen adapted to the
            ! tensor dimenions. For batched contraction the tensor dimensions should be divided by
            ! the number of batches (number of index ranges).
            DO idim = 1, ndims_tensor(tensor)
               CALL get_ith_array(tensor%contraction_storage%batch_ranges, idim, array)
               dims(idim) = array(nbatches(idim) + 1) - array(1)
               DEALLOCATE (array)
               dims(idim) = dims(idim)/nbatches(idim)
               IF (dims(idim) <= 0) dims(idim) = 1
            END DO
         END ASSOCIATE
      END IF

      pgrid = dbt_nd_mp_comm(mp_comm, map1, map2, pdims_2d=pdims, tdims=dims, nsplit=nsplit, dimsplit=dimsplit)
      IF (ALLOCATED(tensor%contraction_storage)) THEN
         #:for ndim in range(1, maxdim+1)
            IF (ndims_tensor(tensor) == ${ndim}$) THEN
               CALL get_arrays(tensor%contraction_storage%batch_ranges, ${varlist("batch_range", nmax=ndim)}$)
               CALL dbt_change_pgrid(tensor, pgrid, ${varlist("batch_range", nmax=ndim)}$, &
                                     nodata=nodata, pgrid_changed=pgrid_changed, unit_nr=unit_nr)
            END IF
         #:endfor
      ELSE
         CALL dbt_change_pgrid(tensor, pgrid, nodata=nodata, pgrid_changed=pgrid_changed, unit_nr=unit_nr)
      END IF
      CALL dbt_pgrid_destroy(pgrid)

   END SUBROUTINE

END MODULE
