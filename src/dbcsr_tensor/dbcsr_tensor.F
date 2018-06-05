!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBCSR tensor framework for block-sparse tensor contraction.
!>
!> Representation of n-rank tensors as DBCSR matrices.
!> Support for arbitrary redistribution between different representations.
!> Support for arbitrary tensor contractions
!> DBCSR routines are generalized to n dimensions by light-weight wrapper routines.
!> \todo implement communication optimal scheme for multiplication of non-square matrices.
!>       implement checks and error messages
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor

#:include "dbcsr_tensor.fypp"
#:set maxdim = maxrank
#:set ndims = range(2,maxdim+1)

   USE allocate_wrap,                   ONLY: allocate_any
   USE array_list_methods,              ONLY: get_arrays,&
                                              reorder_arrays,&
                                              get_ith_array
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_type, dbcsr_create, dbcsr_iterator_type, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, &
        dbcsr_finalize, dbcsr_transpose, dbcsr_no_transpose, dbcsr_get_data_type, dbcsr_multiply, &
        dbcsr_get_info, ${uselist(dtype_float_param)}$
   USE dbcsr_tensor_block,              ONLY: dbcsr_t_iterator_type, &
                                              dbcsr_t_get_block, &
                                              dbcsr_t_put_block, &
                                              dbcsr_t_iterator_start, &
                                              dbcsr_t_iterator_blocks_left, &
                                              dbcsr_t_iterator_stop, &
                                              dbcsr_t_iterator_next_block, &
                                              ndims_iterator,&
                                              dbcsr_t_reserve_blocks, &
                                              block_nd
   USE dbcsr_tensor_index,              ONLY: get_mapping_info,&
                                              nd_to_2d_mapping,&
                                              dbcsr_t_inverse_order,&
                                              permute_index
   USE dbcsr_tensor_types,              ONLY: dbcsr_t_create,&
                                              dbcsr_t_get_data_type,&
                                              dbcsr_t_type,&
                                              ndims_tensor,&
                                              dims_tensor,&
                                              dbcsr_t_distribution_type,&
                                              dbcsr_t_distribution,&
                                              dbcsr_t_nd_mp_comm,&
                                              dbcsr_t_destroy,&
                                              dbcsr_t_distribution_destroy,&
                                              dbcsr_t_distribution_new,&
                                              dbcsr_t_get_stored_coordinates
   USE kinds,                           ONLY: ${uselist(dtype_float_prec)}$,&
                                              default_string_length
   USE message_passing,                 ONLY: mp_environ,&
                                              mp_comm_free
   USE util,                            ONLY: sort
   USE dbcsr_tensor_reshape,            ONLY: dbcsr_t_reshape

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor'


   PUBLIC :: &
      dbcsr_t_contract, &
      dbcsr_t_copy, &
      dbcsr_t_dims, &
      dbcsr_t_get_block, &
      dbcsr_t_get_stored_coordinates, &
      dbcsr_t_inverse_order, &
      dbcsr_t_iterator_blocks_left, &
      dbcsr_t_iterator_next_block, &
      dbcsr_t_iterator_start, &
      dbcsr_t_iterator_stop, &
      dbcsr_t_iterator_type, &
      dbcsr_t_ndims, &
      dbcsr_t_put_block, &
      dbcsr_t_reserve_blocks, &
      dbcsr_t_split_blocks

   INTERFACE dbcsr_t_ndims
      MODULE PROCEDURE ndims_tensor
   END INTERFACE

   INTERFACE dbcsr_t_dims
      MODULE PROCEDURE dims_tensor
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Copy tensor data.
!>        Redistributes tensor data according to distributions of target and source tensor.
!>        Permutes tensor index according to `order` argument (if present).
!>
!> Source and target tensor formats are arbitrary as long as the following requirements are met:
!> * source and target tensors have the same number of blocks in each dimension and the same block sizes.
!>   If `order` argument is present, this must be the case after index permutation.
!> OR
!> * target tensor is not yet created, in this case an exact copy of source tensor is returned.
!> \param tensor_in Source
!> \param tensor_out Target
!> \param order Permutation of target tensor index. Exact same convention as order argument of RESHAPE
!>              intrinsic
! **************************************************************************************************
   SUBROUTINE dbcsr_t_copy(tensor_in, tensor_out, order)
      TYPE(dbcsr_t_type), INTENT(INOUT)              :: tensor_in, tensor_out
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(IN), OPTIONAL                        :: order

      TYPE(dbcsr_t_type)                             :: tensor_tmp
      INTEGER                                        :: handle

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_copy', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)

      IF (PRESENT(order)) THEN
         CALL dbcsr_t_permute_index(tensor_in, tensor_tmp, order)
         CALL dbcsr_t_reshape(tensor_tmp, tensor_out)
         CALL dbcsr_t_destroy(tensor_tmp)
      ELSE
         CALL dbcsr_t_reshape(tensor_in, tensor_out)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief copy matrix to tensor.
!> \param matrix_in ...
!> \param tensor_out ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_copy_matrix(matrix_in, tensor_out)
      TYPE(dbcsr_type), INTENT(IN)                        :: matrix_in
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor_out

      INTEGER, DIMENSION(2)                              :: ind_2d
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :)    :: block_arr
      REAL(KIND=real_8), DIMENSION(:, :), POINTER        :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      INTEGER                                            :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_copy_matrix', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)
      CPASSERT(tensor_out%valid)

      CALL dbcsr_t_reserve_blocks(matrix_in, tensor_out)

      CALL dbcsr_iterator_start(iter, matrix_in)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, ind_2d(1), ind_2d(2), block)
         CALL allocate_any(block_arr, source=block)
         CALL dbcsr_t_put_block(tensor_out, ind_2d, SHAPE(block_arr), block_arr)
         DEALLOCATE (block_arr)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Contract tensors by multiplying matrix representations.
!> \param tensor_1 first tensor (in)
!> \param tensor_2 second tensor (in)
!> \param contract_1 indices of tensor_1 to contract
!> \param notcontract_1 indices of tensor_1 not to contract
!> \param contract_2 indices of tensor_2 to contract (1:1 with contract_1)
!> \param notcontract_2 indices of tensor_2 not to contract
!> \param map_1 which indices of tensor_3 map to non-contracted indices of tensor_1 (1:1 with notcontract_1)
!> \param map_2 which indices of tensor_3 map to non-contracted indices of tensor_2 (1:1 with notcontract_2)
!> \param tensor_3 contracted tensor (out)
!> \todo Generate input arguments (contract_*, notcontract_*, map_*) from Einstein notation
! **************************************************************************************************
   SUBROUTINE dbcsr_t_contract(tensor_1, tensor_2, tensor_3, &
                               contract_1, notcontract_1, &
                               contract_2, notcontract_2, &
                               map_1, map_2)
      TYPE(dbcsr_t_type), INTENT(IN)                 :: tensor_1, tensor_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: contract_1, contract_2, map_1, map_2
      INTEGER, DIMENSION(:), INTENT(IN)              :: notcontract_1, notcontract_2
      TYPE(dbcsr_t_type), INTENT(INOUT)              :: tensor_3
      TYPE(dbcsr_t_type)                             :: tensor_contr_1, tensor_contr_2, tensor_tmp_1, tensor_tmp_2
      LOGICAL                                        :: assert_stmt
      INTEGER                                        :: i, data_type, mp_comm_2d, mp_comm_nd, handle
      INTEGER, DIMENSION(2)                          :: ndims_2d_3
      INTEGER, DIMENSION(:), ALLOCATABLE             :: pdims_3, myploc, sort_indices
      INTEGER, DIMENSION(SIZE(contract_1))           :: contract_1_mod
      INTEGER, DIMENSION(SIZE(notcontract_1))        :: notcontract_1_mod
      INTEGER, DIMENSION(SIZE(contract_2))           :: contract_2_mod
      INTEGER, DIMENSION(SIZE(notcontract_2))        :: notcontract_2_mod
      INTEGER, DIMENSION(SIZE(map_1))                :: map_1_mod, map_1_tmp
      INTEGER, DIMENSION(SIZE(map_2))                :: map_2_mod, map_2_tmp
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_3))    :: order
      LOGICAL, DIMENSION(:), ALLOCATABLE             :: periods
      CHARACTER(LEN=1)                               :: trans_1, trans_2
      INTEGER, DIMENSION(:), ALLOCATABLE             :: ${varlist("blk_sizes_3")}$
      INTEGER, DIMENSION(:), ALLOCATABLE             :: ${varlist("dist_3")}$
      TYPE(dbcsr_t_distribution_type)                :: dist_3
      TYPE(dbcsr_t_type)                             :: tensor_3_tmp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_contract', &
         routineP = moduleN//':'//routineN

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

      assert_stmt = SIZE(notcontract_1) + SIZE(contract_1) .EQ. dbcsr_t_ndims(tensor_1)
      CPASSERT(assert_stmt)

      assert_stmt = SIZE(notcontract_2) + SIZE(contract_2) .EQ. dbcsr_t_ndims(tensor_2)
      CPASSERT(assert_stmt)

      assert_stmt = SIZE(map_1) + SIZE(map_2) .EQ. dbcsr_t_ndims(tensor_3)
      CPASSERT(assert_stmt)

      assert_stmt = dbcsr_t_get_data_type(tensor_1) .EQ. dbcsr_t_get_data_type(tensor_2)
      CPASSERT(assert_stmt)

      ndims_2d_3(1) = SIZE(notcontract_1)
      ndims_2d_3(2) = SIZE(notcontract_2)
      data_type = dbcsr_t_get_data_type(tensor_1)

      CALL align_tensor(tensor_1, contract_1, notcontract_1, &
                        tensor_tmp_1, contract_1_mod, notcontract_1_mod)

      CALL align_tensor(tensor_2, contract_2, notcontract_2, &
                        tensor_tmp_2, contract_2_mod, notcontract_2_mod)

      CALL allocate_any(sort_indices, shape_spec=SHAPE(contract_1_mod))
      CALL sort(contract_1_mod, SIZE(contract_1_mod), sort_indices)
      contract_2_mod(:) = contract_2_mod(sort_indices)
      DEALLOCATE(sort_indices)

      CALL allocate_any(sort_indices, shape_spec=SHAPE(notcontract_1_mod))
      CALL sort(notcontract_1_mod, SIZE(notcontract_1_mod), sort_indices)
      map_1_mod(:) = map_1(sort_indices)
      DEALLOCATE(sort_indices)

      CALL allocate_any(sort_indices, shape_spec=SHAPE(notcontract_2_mod))
      CALL sort(notcontract_2_mod, SIZE(notcontract_2_mod), sort_indices)
      map_2_mod(:) = map_2(sort_indices)
      DEALLOCATE(sort_indices)

      CALL prep_tensor_lhs(tensor_tmp_1, contract_1_mod, notcontract_1_mod, 2, tensor_contr_1, trans_1)
      CALL prep_tensor_lhs(tensor_tmp_2, contract_2_mod, notcontract_2_mod, 1, tensor_contr_2, trans_2)

      CALL dbcsr_t_destroy(tensor_tmp_1)
      CALL dbcsr_t_destroy(tensor_tmp_2)

      ! create temporary tensor according to format of tensor 1, tensor 2
      map_1_tmp(:) = (/(i, i=1,ndims_2d_3(1))/)
      map_2_tmp(:) = (/(i, i=ndims_2d_3(1)+1,ndims_2d_3(1)+ndims_2d_3(2))/)

      CALL dbcsr_get_info(tensor_contr_1%matrix_rep, group=mp_comm_2d)
      mp_comm_nd = dbcsr_t_nd_mp_comm(mp_comm_2d, map_1_tmp, map_2_tmp)

      ALLOCATE(pdims_3(SUM(ndims_2d_3)))
      ALLOCATE(myploc(SUM(ndims_2d_3)))
      ALLOCATE(periods(SUM(ndims_2d_3)))
      CALL mp_environ(mp_comm_nd, SUM(ndims_2d_3), pdims_3, myploc, periods)

#:for dim in range(1, maxdim+1)
      IF (${dim}$ <= ndims_2d_3(1)) THEN
         CALL get_ith_array(tensor_contr_1%blk_sizes, notcontract_1_mod(${dim}$), blk_sizes_3_${dim}$)
      ELSEIF (${dim}$ <= SUM(ndims_2d_3)) THEN
         CALL get_ith_array(tensor_contr_2%blk_sizes, notcontract_2_mod(${dim}$-ndims_2d_3(1)), blk_sizes_3_${dim}$)
      ENDIF
#:endfor
      ! copy to tensor_3
#:for dim in range(1, maxdim+1)
      IF (${dim}$ <= SUM(ndims_2d_3)) THEN
         CALL cyclic_dist(dist_3_${dim}$, SIZE(blk_sizes_3_${dim}$), pdims_3(${dim}$))
      ENDIF
#:endfor

#:for ndim in ndims
      IF (SUM(ndims_2d_3) == ${ndim}$) THEN
         CALL dbcsr_t_distribution_new(dist_3, mp_comm_nd, mp_comm_2d, &
                                       map_1_tmp, map_2_tmp, &
                                       ${varlist("dist_3", nmax=ndim)}$)
         CALL dbcsr_t_create(tensor_3_tmp, "temporary contraction result", dist_3, &
                             map_1_tmp, map_2_tmp, &
                             data_type, &
                             ${varlist("blk_sizes_3", nmax=ndim)}$)
      ENDIF
#:endfor

      CALL dbcsr_multiply(trans_1, trans_2, 1.0_real_8, &
                          tensor_contr_1%matrix_rep, tensor_contr_2%matrix_rep, 0.0_real_8, &
                          tensor_3_tmp%matrix_rep)

      order(:) = [map_1_mod, map_2_mod]
      CALL dbcsr_t_copy(tensor_3_tmp, tensor_3, order = order)

      CALL dbcsr_t_destroy(tensor_contr_1)
      CALL dbcsr_t_destroy(tensor_contr_2)
      CALL dbcsr_t_destroy(tensor_3_tmp)
      CALL dbcsr_t_distribution_destroy(dist_3)
      CALL mp_comm_free(mp_comm_nd)

      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief align tensor index with data
!> \param tensor_in ...
!> \param contract_in ...
!> \param notcontract_in ...
!> \param tensor_out ...
!> \param contract_out ...
!> \param notcontract_out ...
! **************************************************************************************************
      SUBROUTINE align_tensor(tensor_in, contract_in, notcontract_in, &
                              tensor_out, contract_out, notcontract_out)
         TYPE(dbcsr_t_type), INTENT(IN)               :: tensor_in
         INTEGER, DIMENSION(:), INTENT(IN)            :: contract_in, notcontract_in
         TYPE(dbcsr_t_type), INTENT(OUT)              :: tensor_out
         INTEGER, DIMENSION(SIZE(contract_in)), &
            INTENT(OUT)                               :: contract_out
         INTEGER, DIMENSION(SIZE(notcontract_in)), &
            INTENT(OUT)                               :: notcontract_out
         INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)) :: align

         CALL dbcsr_t_align_index(tensor_in, tensor_out, order=align)
         contract_out = align(contract_in)
         notcontract_out = align(notcontract_in)

      END SUBROUTINE

! **************************************************************************************************
!> \brief Prepare tensor for contraction: redistribute to a 2d format which can be contracted by
!>        matrix multiplication.
!> \param tensor_in ...
!> \param contract ...
!> \param notcontract ...
!> \param contract_index ...
!> \param tensor_out ...
!> \param trans ...
! **************************************************************************************************
      SUBROUTINE prep_tensor_lhs(tensor_in, contract, notcontract, contract_index, tensor_out, trans)
         TYPE(dbcsr_t_type), INTENT(INOUT)           :: tensor_in
         INTEGER, DIMENSION(:), INTENT(IN)           :: contract, notcontract
         INTEGER, INTENT(IN)                         :: contract_index
         TYPE(dbcsr_t_type), INTENT(OUT)             :: tensor_out
         CHARACTER(LEN=1), INTENT(OUT)               :: trans
         CHARACTER(LEN=1)                            :: compat_1_trans, compat_2_trans
         INTEGER                                     :: compat

         compat = compat_map(tensor_in%nd_index, contract)
         IF (compat == 0) THEN ! index mapping not compatible with contract index
            CALL dbcsr_t_remap(tensor_in, contract, notcontract, tensor_out)
            compat = compat_map(tensor_out%nd_index, contract)
         ELSE
            CALL dbcsr_t_copy(tensor_in, tensor_out)
         ENDIF

         CPASSERT(compat > 0)

         SELECT CASE(contract_index)
         CASE(1)
            compat_1_trans = dbcsr_no_transpose
            compat_2_trans = dbcsr_transpose
         CASE(2)
            compat_1_trans = dbcsr_transpose
            compat_2_trans = dbcsr_no_transpose
         CASE DEFAULT
            CPABORT('')
         END SELECT

         SELECT CASE(compat)
         CASE(1)
            trans = compat_1_trans
         CASE(2)
            trans = compat_2_trans
         CASE DEFAULT
            CPABORT('')
         END SELECT

      END SUBROUTINE

! **************************************************************************************************
!> \brief Check if contraction index is compatible with tensor index
!> \param nd_index ...
!> \param contract ...
! **************************************************************************************************
      FUNCTION compat_map(nd_index, contract)
         TYPE(nd_to_2d_mapping), INTENT(IN) :: nd_index
         INTEGER, DIMENSION(:), INTENT(IN)  :: contract
         INTEGER, DIMENSION(:), ALLOCATABLE :: map1, map2
         INTEGER                            :: compat_map

         CALL get_mapping_info(nd_index, map1_2d=map1, map2_2d=map2)

         compat_map = 0
         IF(array_eq_i(map1, contract)) THEN
            compat_map = 1
         ELSEIF(array_eq_i(map2, contract)) THEN
            compat_map = 2
         ENDIF

      END FUNCTION

! **************************************************************************************************
!> \brief Check if 2 arrays are equal
!> \param arr1 ...
!> \param arr2
! **************************************************************************************************
      PURE FUNCTION array_eq_i(arr1, arr2)
         INTEGER, INTENT(IN), DIMENSION(:) :: arr1
         INTEGER, INTENT(IN), DIMENSION(:) :: arr2
         LOGICAL                           :: array_eq_i

         array_eq_i = .FALSE.
         IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)

      END FUNCTION

   END SUBROUTINE

! **************************************************************************************************
!> \brief set up cyclic distribution. This is used for all internally created temporary tensors.
!> \param dist_array ...
!> \param dist_size ...
!> \param nbins ...
! **************************************************************************************************
   SUBROUTINE cyclic_dist(dist_array, dist_size, nbins)
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: dist_array
      INTEGER, INTENT(in)                             :: dist_size, nbins
      INTEGER                                         :: i

      CALL allocate_any(dist_array, source=(/(MODULO(nbins-i, nbins), i=1, dist_size)/))

   END SUBROUTINE cyclic_dist

! **************************************************************************************************
!> \brief Copy tensor to tensor with modified index mapping
!> \param tensor_in ...
!> \param map1_2d new index mapping
!> \param map2_2d new index mapping
!> \param tensor_out ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_remap(tensor_in, map1_2d, map2_2d, tensor_out, name)
      TYPE(dbcsr_t_type), INTENT(INOUT)      :: tensor_in
      INTEGER, DIMENSION(:), INTENT(IN)      :: map1_2d, map2_2d
      TYPE(dbcsr_t_type), INTENT(OUT)        :: tensor_out
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: name
      CHARACTER(len=default_string_length)   :: name_tmp
      INTEGER, DIMENSION(:), ALLOCATABLE     :: ${varlist("blk_sizes")}$,&
                                                ${varlist("nd_dist")}$
      TYPE(dbcsr_t_distribution_type)        :: dist
      INTEGER                                :: comm_2d, comm_nd, handle
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)) :: pdims, myploc
      LOGICAL, DIMENSION(dbcsr_t_ndims(tensor_in)) :: periods
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_remap', &
         routineP = moduleN//':'//routineN

      CALL timeset(routineN, handle)

      IF (PRESENT(name)) THEN
         name_tmp = name
      ELSE
         name_tmp = tensor_in%name
      ENDIF

      CALL dbcsr_get_info(tensor_in%matrix_rep, group=comm_2d)
      comm_nd = dbcsr_t_nd_mp_comm(comm_2d, map1_2d, map2_2d)

      CALL mp_environ(comm_nd, dbcsr_t_ndims(tensor_in), pdims, myploc, periods)

#:for ndim in ndims
      IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
         CALL get_arrays(tensor_in%blk_sizes, ${varlist("blk_sizes", nmax=ndim)}$)
      ENDIF
#:endfor

      ! Note: it is not possible to derive distribution from tensor_in, since process grid has changed
#:for ndim in ndims
      IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
#:for idim in range(1, ndim+1)
         CALL cyclic_dist(nd_dist_${idim}$, SIZE(blk_sizes_${idim}$), pdims(${idim}$))
#:endfor
         CALL dbcsr_t_distribution_new(dist, comm_nd, comm_2d, map1_2d, map2_2d, &
                                       ${varlist("nd_dist", nmax=ndim)}$)
      ENDIF
#:endfor

#:for ndim in ndims
      IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
         CALL dbcsr_t_create(tensor_out, name_tmp, dist, &
                             map1_2d, map2_2d, dbcsr_get_data_type(tensor_in%matrix_rep),&
                             ${varlist("blk_sizes", nmax=ndim)}$)
         tensor_out%comm_prv = .TRUE. ! since this tensor is created internally (not by user)
      ENDIF
#:endfor

      CALL dbcsr_t_copy(tensor_in, tensor_out)
      CALL dbcsr_t_distribution_destroy(dist)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Align index with data
!> \param tensor_in ...
!> \param tensor_out ...
!> \param order permutation resulting from alignment
! **************************************************************************************************
   SUBROUTINE dbcsr_t_align_index(tensor_in, tensor_out, order)
      TYPE(dbcsr_t_type), INTENT(IN)               :: tensor_in
      TYPE(dbcsr_t_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(:), ALLOCATABLE              :: map1_2d, map2_2d
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(OUT), OPTIONAL                        :: order
      INTEGER, DIMENSION(ndims_tensor(tensor_in))     :: order_prv
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_align_index', &
         routineP = moduleN//':'//routineN
      INTEGER                                         :: handle

      CALL timeset(routineN, handle)

      CALL get_mapping_info(tensor_in%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      order_prv = dbcsr_t_inverse_order([map1_2d, map2_2d])
      CALL dbcsr_t_permute_index(tensor_in, tensor_out, order=order_prv)

      IF(PRESENT(order)) order = order_prv

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create new tensor by reordering index, data is copied exactly
!> \param tensor_in ...
!> \param tensor_out ...
!> \param order ...
! **************************************************************************************************
   SUBROUTINE dbcsr_t_permute_index(tensor_in, tensor_out, order)
      TYPE(dbcsr_t_type), INTENT(IN)                  :: tensor_in
      TYPE(dbcsr_t_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(ndims_tensor(tensor_in)), &
         INTENT(IN)                                   :: order

      TYPE(nd_to_2d_mapping)                          :: nd_index_blk_rs, nd_index_rs
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbcsr_t_permute_index', &
         routineP = moduleN//':'//routineN
      INTEGER                                         :: handle

      CALL timeset(routineN, handle)

      CALL permute_index(tensor_in%nd_index, nd_index_rs, order)
      CALL permute_index(tensor_in%nd_index_blk, nd_index_blk_rs, order)

      CALL dbcsr_create(tensor_out%matrix_rep, template=tensor_in%matrix_rep)
      CALL dbcsr_finalize(tensor_out%matrix_rep)

      CALL dbcsr_copy(tensor_out%matrix_rep, tensor_in%matrix_rep, shallow_data=.TRUE.)
      tensor_out%nd_index = nd_index_rs
      tensor_out%nd_index_blk = nd_index_blk_rs

      CALL reorder_arrays(tensor_in%blk_sizes, tensor_out%blk_sizes, order)
      CALL reorder_arrays(tensor_in%blk_offsets, tensor_out%blk_offsets, order)
      CALL reorder_arrays(tensor_in%nd_dist, tensor_out%nd_dist, order)
      tensor_out%name = tensor_in%name
      tensor_out%valid = .TRUE.

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Split tensor blocks into smaller blocks of maximum size PRODUCT(block_sizes).
!> \param tensor_in Input tensor
!> \param tensor_out Output tensor (splitted blocks)
!> \param block_sizes block sizes for each of the tensor dimensions
! **************************************************************************************************
   SUBROUTINE dbcsr_t_split_blocks(tensor_in, tensor_out, block_sizes)
      TYPE(dbcsr_t_type), INTENT(INOUT)               :: tensor_in
      TYPE(dbcsr_t_type), INTENT(OUT)                 :: tensor_out
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)), &
         INTENT(IN)                                   :: block_sizes

      TYPE(dbcsr_t_distribution_type)                 :: dist_old, dist_split
      TYPE(dbcsr_t_iterator_type)                     :: iter
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("nd_dist_split")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("nd_blk_size_split")}$
      INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("index_split_offset")}$
      INTEGER                                         :: ${varlist("split_blk")}$
      INTEGER :: idim, i, isplit_sum, blk_remainder, comm_2d, blk, nsplit, isplit
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist_d, blk_size_d, blk_size_split_d, dist_split_d, &
         map1_2d, map2_2d
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor_in)) :: blk_index, blk_size, blk_offset, nsplit_blk, &
         isplit_sum_blk, blk_remainder_nd, blk_shape
      INTEGER, DIMENSION(${maxdim}$) :: split_index

#:for dparam, dtype, dsuffix in dtype_float_list
#:for ndim in ndims
      ${dtype}$, DIMENSION(${shape_colon(n=ndim)}$), ALLOCATABLE :: block_${dsuffix}$_${ndim}$d
#:endfor
#:endfor

      dist_old = dbcsr_t_distribution(tensor_in)

      DO idim = 1, ndims_tensor(tensor_in)
         CALL get_ith_array(dist_old%nd_dist, idim, dist_d)
         CALL get_ith_array(tensor_in%blk_sizes, idim, blk_size_d)

#:for idim in range(1, maxdim+1)
         IF (idim == ${idim}$) THEN
            ALLOCATE(index_split_offset_${idim}$(SIZE(dist_d)))
         ENDIF
#:endfor

         isplit_sum=0
         DO i = 1, SIZE(dist_d)
            nsplit = (blk_size_d(i)+block_sizes(idim)-1)/block_sizes(idim)
            isplit_sum = isplit_sum+nsplit
         ENDDO
         ALLOCATE(dist_split_d(isplit_sum))
         ALLOCATE(blk_size_split_d(isplit_sum))

         isplit_sum = 0
         DO i = 1, SIZE(dist_d)
            nsplit = (blk_size_d(i)+block_sizes(idim)-1)/block_sizes(idim)
            blk_remainder = blk_size_d(i)
#:for idim in range(1, maxdim+1)
            IF(idim == ${idim}$) THEN
               index_split_offset_${idim}$(i) = isplit_sum
            ENDIF
#:endfor
            DO isplit = 1, nsplit
               isplit_sum = isplit_sum + 1
               dist_split_d(isplit_sum) = dist_d(i)
               blk_size_split_d(isplit_sum) = MIN(block_sizes(idim), blk_remainder)
               blk_remainder = blk_remainder - block_sizes(idim)
            ENDDO

         ENDDO

#:for idim in range(1, maxdim+1)
         IF(idim == ${idim}$) THEN
            CALL allocate_any(nd_dist_split_${idim}$, source=dist_split_d)
            CALL allocate_any(nd_blk_size_split_${idim}$, source=blk_size_split_d)
         ENDIF
#:endfor
         DEALLOCATE(dist_split_d)
         DEALLOCATE(blk_size_split_d)

      ENDDO

      CALL dbcsr_get_info(tensor_in%matrix_rep, group=comm_2d)
      CALL get_mapping_info(tensor_in%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)

#:for ndim in ndims
      IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
         CALL dbcsr_t_distribution_new(dist_split, tensor_in%comm_nd, comm_2d, map1_2d, map2_2d, &
                                       ${varlist("nd_dist_split", nmax=ndim)}$)
         CALL dbcsr_t_create(tensor_out, tensor_in%name, dist_split, map1_2d, map2_2d, &
                             dbcsr_t_get_data_type(tensor_in), ${varlist("nd_blk_size_split", nmax=ndim)}$)
      ENDIF
#:endfor

      CALL dbcsr_t_distribution_destroy(dist_split)

      CALL dbcsr_t_iterator_start(iter, tensor_in)

      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, blk_index, blk, blk_size=blk_size, blk_offset=blk_offset)
#:for dparam, dtype, dsuffix in dtype_float_list
         IF (dbcsr_t_get_data_type(tensor_in) == ${dparam}$) THEN
#:for ndim in ndims
            IF (dbcsr_t_ndims(tensor_in) == ${ndim}$) THEN
               CALL dbcsr_t_get_block(tensor_in, blk_index, block_${dsuffix}$_${ndim}$d)
            ENDIF
#:endfor
         ENDIF
#:endfor

         nsplit_blk = (blk_size + block_sizes - 1)/block_sizes

#:for ndim in ndims
         IF (dbcsr_t_ndims(tensor_in) == ${ndim}$) THEN
#:for idim in range(1,ndim+1)

            blk_remainder_nd(${idim}$) = blk_size(${idim}$)
            isplit_sum_blk(${idim}$) = 0
            DO split_blk_${idim}$ = 1, nsplit_blk(${idim}$)
               isplit_sum_blk(${idim}$) = isplit_sum_blk(${idim}$) + 1

               blk_shape(${idim}$) = MIN(block_sizes(${idim}$), blk_remainder_nd(${idim}$))
               blk_remainder_nd(${idim}$) = blk_remainder_nd(${idim}$) - block_sizes(${idim}$)
#:endfor

               split_index(1:${ndim}$) = [${", ".join(["index_split_offset_"+str(idim)+"(blk_index("+str(idim)+")) + isplit_sum_blk("+str(idim)+")" for idim in range(1, ndim+1)])}$]

#:for dparam, dtype, dsuffix in dtype_float_list

               IF (dbcsr_t_get_data_type(tensor_in) == ${dparam}$) THEN

                  CALL dbcsr_t_put_block(tensor_out, split_index(1:${ndim}$), &
                                         blk_shape, &
                                         block_${dsuffix}$_${ndim}$d( &
                                         ${", ".join(["(isplit_sum_blk("+str(idim)+") - 1)*block_sizes("+str(idim)+") + 1:MIN(isplit_sum_blk("+str(idim)+")*block_sizes("+str(idim)+"), blk_size("+str(idim)+"))" for idim in range(1, ndim+1)])}$))
               ENDIF
#:endfor

#:for idim in range(1,ndim+1)
            ENDDO
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
            IF (dbcsr_t_get_data_type(tensor_in) == ${dparam}$) THEN
               DEALLOCATE(block_${dsuffix}$_${ndim}$d)
            ENDIF
#:endfor
         ENDIF
#:endfor
      ENDDO
      CALL dbcsr_t_iterator_stop(iter)

      CALL dbcsr_finalize(tensor_out%matrix_rep)

   END SUBROUTINE

END MODULE
