!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief DBT tensor framework for block-sparse tensor contraction: Types and create/destroy routines.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_types
   #:include "dbt_macros.fypp"
   #:set maxdim = maxrank
   #:set ndims = range(2,maxdim+1)

   USE dbcsr_api, ONLY: dbcsr_type, dbcsr_get_info, dbcsr_distribution_type, dbcsr_distribution_get
   USE dbt_array_list_methods, ONLY: &
      array_list, array_offsets, create_array_list, destroy_array_list, get_array_elements, &
      sizes_of_arrays, sum_of_arrays, array_sublist, get_arrays, get_ith_array, array_eq_i
   USE dbm_api, ONLY: &
      dbm_distribution_obj, dbm_type
   USE kinds, ONLY: dp, dp, default_string_length
   USE dbt_tas_base, ONLY: &
      dbt_tas_create, dbt_tas_distribution_new, &
      dbt_tas_distribution_destroy, dbt_tas_finalize, dbt_tas_get_info, &
      dbt_tas_destroy, dbt_tas_get_stored_coordinates, dbt_tas_filter, &
      dbt_tas_get_num_blocks, dbt_tas_get_num_blocks_total, dbt_tas_get_nze, &
      dbt_tas_get_nze_total, dbt_tas_clear
   USE dbt_tas_types, ONLY: &
      dbt_tas_type, dbt_tas_distribution_type, dbt_tas_split_info, dbt_tas_mm_storage
   USE dbt_tas_mm, ONLY: dbt_tas_set_batched_state
   USE dbt_index, ONLY: &
      get_2d_indices_tensor, get_nd_indices_pgrid, create_nd_to_2d_mapping, destroy_nd_to_2d_mapping, &
      dbt_get_mapping_info, nd_to_2d_mapping, split_tensor_index, combine_tensor_index, combine_pgrid_index, &
      split_pgrid_index, ndims_mapping, ndims_mapping_row, ndims_mapping_column
   USE dbt_tas_split, ONLY: &
      dbt_tas_create_split_rows_or_cols, dbt_tas_release_info, dbt_tas_info_hold, &
      dbt_tas_create_split, dbt_tas_get_split_info, dbt_tas_set_strict_split
   USE kinds, ONLY: default_string_length, int_8, dp
   USE message_passing, ONLY: &
      mp_cart_create, mp_cart_rank, mp_environ, mp_dims_create, mp_comm_free, mp_comm_dup, mp_sum, mp_max
   USE dbt_tas_global, ONLY: dbt_tas_distribution, dbt_tas_rowcol_data, dbt_tas_default_distvec
   USE dbt_allocate_wrap, ONLY: allocate_any
   USE dbm_api, ONLY: dbm_scale
   USE util, ONLY: sort
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_types'

   PUBLIC  :: &
      blk_dims_tensor, &
      dbt_blk_offsets, &
      dbt_blk_sizes, &
      dbt_clear, &
      dbt_create, &
      dbt_destroy, &
      dbt_distribution, &
      dbt_distribution_destroy, &
      dbt_distribution_new, &
      dbt_distribution_new_expert, &
      dbt_distribution_type, &
      dbt_filter, &
      dbt_finalize, &
      dbt_get_info, &
      dbt_get_num_blocks, &
      dbt_get_num_blocks_total, &
      dbt_get_nze, &
      dbt_get_nze_total, &
      dbt_get_stored_coordinates, &
      dbt_hold, &
      dbt_mp_dims_create, &
      dbt_nd_mp_comm, &
      dbt_nd_mp_free, &
      dbt_pgrid_change_dims, &
      dbt_pgrid_create, &
      dbt_pgrid_create_expert, &
      dbt_pgrid_destroy, &
      dbt_pgrid_type, &
      dbt_pgrid_set_strict_split, &
      dbt_scale, &
      dbt_type, &
      dims_tensor, &
      mp_environ_pgrid, &
      ndims_tensor, &
      ndims_matrix_row, &
      ndims_matrix_column, &
      dbt_nblks_local, &
      dbt_nblks_total, &
      dbt_blk_size, &
      dbt_max_nblks_local, &
      dbt_default_distvec, &
      dbt_contraction_storage, &
      dbt_copy_contraction_storage

   TYPE dbt_pgrid_type
      TYPE(nd_to_2d_mapping)                  :: nd_index_grid
      INTEGER                                 :: mp_comm_2d
      TYPE(dbt_tas_split_info), ALLOCATABLE :: tas_split_info
      INTEGER                                 :: nproc
   END TYPE

   TYPE dbt_contraction_storage
      REAL(dp) :: nsplit_avg
      INTEGER :: ibatch
      TYPE(array_list) :: batch_ranges
      LOGICAL :: static
   END TYPE

   TYPE dbt_type
      TYPE(dbt_tas_type), POINTER        :: matrix_rep => NULL()
      TYPE(nd_to_2d_mapping)               :: nd_index_blk
      TYPE(nd_to_2d_mapping)               :: nd_index
      TYPE(array_list)                     :: blk_sizes
      TYPE(array_list)                     :: blk_offsets
      TYPE(array_list)                     :: nd_dist
      TYPE(dbt_pgrid_type)             :: pgrid
      TYPE(array_list)                     :: blks_local
      INTEGER, DIMENSION(:), ALLOCATABLE   :: nblks_local
      INTEGER, DIMENSION(:), ALLOCATABLE   :: nfull_local
      LOGICAL                              :: valid = .FALSE.
      LOGICAL                              :: owns_matrix = .FALSE.
      CHARACTER(LEN=default_string_length) :: name
      ! lightweight reference counting for communicators:
      INTEGER, POINTER                     :: refcount => NULL()
      TYPE(dbt_contraction_storage), ALLOCATABLE :: contraction_storage
   END TYPE dbt_type

   TYPE dbt_distribution_type
      TYPE(dbt_tas_distribution_type) :: dist
      TYPE(dbt_pgrid_type)      :: pgrid
      TYPE(array_list)              :: nd_dist
      ! lightweight reference counting for communicators:
      INTEGER, POINTER :: refcount => NULL()
   END TYPE

   ! tas matrix distribution function object for one matrix index
   TYPE, EXTENDS(dbt_tas_distribution) :: dbt_tas_dist_t
      ! tensor dimensions only for this matrix dimension:
      INTEGER, DIMENSION(:), ALLOCATABLE :: dims
      ! grid dimensions only for this matrix dimension:
      INTEGER, DIMENSION(:), ALLOCATABLE :: dims_grid
      ! dist only for tensor dimensions belonging to this matrix dimension:
      TYPE(array_list)       :: nd_dist
   CONTAINS
      ! map matrix index to process grid:
      PROCEDURE :: dist => tas_dist_t
      ! map process grid to matrix index:
      PROCEDURE :: rowcols => tas_rowcols_t
   END TYPE

   ! block size object for one matrix index
   TYPE, EXTENDS(dbt_tas_rowcol_data) :: dbt_tas_blk_size_t
      ! tensor dimensions only for this matrix dimension:
      INTEGER, DIMENSION(:), ALLOCATABLE :: dims
      ! block size only for this matrix dimension:
      TYPE(array_list) :: blk_size
   CONTAINS
      PROCEDURE :: data => tas_blk_size_t
   END TYPE

   INTERFACE dbt_create
      MODULE PROCEDURE dbt_create_new
      MODULE PROCEDURE dbt_create_template
      MODULE PROCEDURE dbt_create_matrix
   END INTERFACE

   INTERFACE dbt_tas_dist_t
      MODULE PROCEDURE new_dbt_tas_dist_t
   END INTERFACE

   INTERFACE dbt_tas_blk_size_t
      MODULE PROCEDURE new_dbt_tas_blk_size_t
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Create distribution object for one matrix dimension
!> \param nd_dist arrays for distribution vectors along all dimensions
!> \param map_blks tensor to matrix mapping object for blocks
!> \param map_grid tensor to matrix mapping object for process grid
!> \param which_dim for which dimension (1 or 2) distribution should be created
!> \return distribution object
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION new_dbt_tas_dist_t(nd_dist, map_blks, map_grid, which_dim)
      TYPE(array_list), INTENT(IN)       :: nd_dist
      TYPE(nd_to_2d_mapping), INTENT(IN) :: map_blks, map_grid
      INTEGER, INTENT(IN)                :: which_dim

      TYPE(dbt_tas_dist_t)               :: new_dbt_tas_dist_t
      INTEGER, DIMENSION(2)              :: grid_dims
      INTEGER(KIND=int_8), DIMENSION(2)  :: matrix_dims
      INTEGER, DIMENSION(:), ALLOCATABLE :: index_map

      IF (which_dim == 1) THEN
         ALLOCATE (new_dbt_tas_dist_t%dims(ndims_mapping_row(map_blks)))
         ALLOCATE (index_map(ndims_mapping_row(map_blks)))
         CALL dbt_get_mapping_info(map_blks, &
                                   dims_2d_i8=matrix_dims, &
                                   map1_2d=index_map, &
                                   dims1_2d=new_dbt_tas_dist_t%dims)
         ALLOCATE (new_dbt_tas_dist_t%dims_grid(ndims_mapping_row(map_grid)))
         CALL dbt_get_mapping_info(map_grid, &
                                   dims_2d=grid_dims, &
                                   dims1_2d=new_dbt_tas_dist_t%dims_grid)
      ELSEIF (which_dim == 2) THEN
         ALLOCATE (new_dbt_tas_dist_t%dims(ndims_mapping_column(map_blks)))
         ALLOCATE (index_map(ndims_mapping_column(map_blks)))
         CALL dbt_get_mapping_info(map_blks, &
                                   dims_2d_i8=matrix_dims, &
                                   map2_2d=index_map, &
                                   dims2_2d=new_dbt_tas_dist_t%dims)
         ALLOCATE (new_dbt_tas_dist_t%dims_grid(ndims_mapping_column(map_grid)))
         CALL dbt_get_mapping_info(map_grid, &
                                   dims_2d=grid_dims, &
                                   dims2_2d=new_dbt_tas_dist_t%dims_grid)
      ELSE
         CPABORT("Unknown value for which_dim")
      END IF

      new_dbt_tas_dist_t%nd_dist = array_sublist(nd_dist, index_map)
      new_dbt_tas_dist_t%nprowcol = grid_dims(which_dim)
      new_dbt_tas_dist_t%nmrowcol = matrix_dims(which_dim)
   END FUNCTION

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION tas_dist_t(t, rowcol)
      CLASS(dbt_tas_dist_t), INTENT(IN) :: t
      INTEGER(KIND=int_8), INTENT(IN) :: rowcol
      INTEGER, DIMENSION(${maxrank}$) :: ind_blk
      INTEGER, DIMENSION(${maxrank}$) :: dist_blk
      INTEGER :: tas_dist_t

      ind_blk(:SIZE(t%dims)) = split_tensor_index(rowcol, t%dims)
      dist_blk(:SIZE(t%dims)) = get_array_elements(t%nd_dist, ind_blk(:SIZE(t%dims)))
      tas_dist_t = combine_pgrid_index(dist_blk(:SIZE(t%dims)), t%dims_grid)
   END FUNCTION

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION tas_rowcols_t(t, dist)
      CLASS(dbt_tas_dist_t), INTENT(IN) :: t
      INTEGER, INTENT(IN) :: dist
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE :: tas_rowcols_t
      INTEGER, DIMENSION(${maxrank}$) :: dist_blk
      INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("dist")}$, ${varlist("blks")}$, blks_tmp, nd_ind
      INTEGER :: ${varlist("i")}$, i, iblk, iblk_count, nblks
      INTEGER(KIND=int_8) :: nrowcols
      TYPE(array_list) :: blks

      dist_blk(:SIZE(t%dims)) = split_pgrid_index(dist, t%dims_grid)

      #:for ndim in range(1, maxdim+1)
         IF (SIZE(t%dims) == ${ndim}$) THEN
            CALL get_arrays(t%nd_dist, ${varlist("dist", nmax=ndim)}$)
         END IF
      #:endfor

      #:for idim in range(1, maxdim+1)
         IF (SIZE(t%dims) .GE. ${idim}$) THEN
            nblks = SIZE(dist_${idim}$)
            ALLOCATE (blks_tmp(nblks))
            iblk_count = 0
            DO iblk = 1, nblks
               IF (dist_${idim}$ (iblk) == dist_blk(${idim}$)) THEN
                  iblk_count = iblk_count + 1
                  blks_tmp(iblk_count) = iblk
               END IF
            END DO
            ALLOCATE (blks_${idim}$ (iblk_count))
            blks_${idim}$ (:) = blks_tmp(:iblk_count)
            DEALLOCATE (blks_tmp)
         END IF
      #:endfor

      #:for ndim in range(1, maxdim+1)
         IF (SIZE(t%dims) == ${ndim}$) THEN
            CALL create_array_list(blks, ${ndim}$, ${varlist("blks", nmax=ndim)}$)
         END IF
      #:endfor

      nrowcols = PRODUCT(INT(sizes_of_arrays(blks), int_8))
      ALLOCATE (tas_rowcols_t(nrowcols))

      #:for ndim in range(1, maxdim+1)
         IF (SIZE(t%dims) == ${ndim}$) THEN
            ALLOCATE (nd_ind(${ndim}$))
            i = 0
            #:for idim in range(1,ndim+1)
               DO i_${idim}$ = 1, SIZE(blks_${idim}$)
                  #:endfor
                  i = i + 1

                  nd_ind(:) = get_array_elements(blks, [${varlist("i", nmax=ndim)}$])
                  tas_rowcols_t(i) = combine_tensor_index(nd_ind, t%dims)
                  #:for idim in range(1,ndim+1)
                     END DO
                  #:endfor
               END IF
            #:endfor

         END FUNCTION

! **************************************************************************************************
!> \brief Create block size object for one matrix dimension
!> \param blk_size arrays for block sizes along all dimensions
!> \param map_blks tensor to matrix mapping object for blocks
!> \param which_dim for which dimension (1 or 2) distribution should be created
!> \return block size object
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION new_dbt_tas_blk_size_t(blk_size, map_blks, which_dim)
            TYPE(array_list), INTENT(IN)                   :: blk_size
            TYPE(nd_to_2d_mapping), INTENT(IN)             :: map_blks
            INTEGER, INTENT(IN) :: which_dim
            INTEGER(KIND=int_8), DIMENSION(2) :: matrix_dims
            INTEGER, DIMENSION(:), ALLOCATABLE :: index_map
            TYPE(dbt_tas_blk_size_t) :: new_dbt_tas_blk_size_t

            IF (which_dim == 1) THEN
               ALLOCATE (index_map(ndims_mapping_row(map_blks)))
               ALLOCATE (new_dbt_tas_blk_size_t%dims(ndims_mapping_row(map_blks)))
               CALL dbt_get_mapping_info(map_blks, &
                                         dims_2d_i8=matrix_dims, &
                                         map1_2d=index_map, &
                                         dims1_2d=new_dbt_tas_blk_size_t%dims)
            ELSEIF (which_dim == 2) THEN
               ALLOCATE (index_map(ndims_mapping_column(map_blks)))
               ALLOCATE (new_dbt_tas_blk_size_t%dims(ndims_mapping_column(map_blks)))
               CALL dbt_get_mapping_info(map_blks, &
                                         dims_2d_i8=matrix_dims, &
                                         map2_2d=index_map, &
                                         dims2_2d=new_dbt_tas_blk_size_t%dims)
            ELSE
               CPABORT("Unknown value for which_dim")
            END IF

            new_dbt_tas_blk_size_t%blk_size = array_sublist(blk_size, index_map)
            new_dbt_tas_blk_size_t%nmrowcol = matrix_dims(which_dim)

            new_dbt_tas_blk_size_t%nfullrowcol = PRODUCT(INT(sum_of_arrays(new_dbt_tas_blk_size_t%blk_size), &
                                                             KIND=int_8))
         END FUNCTION

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION tas_blk_size_t(t, rowcol)
            CLASS(dbt_tas_blk_size_t), INTENT(IN) :: t
            INTEGER(KIND=int_8), INTENT(IN) :: rowcol
            INTEGER :: tas_blk_size_t
            INTEGER, DIMENSION(SIZE(t%dims)) :: ind_blk
            INTEGER, DIMENSION(SIZE(t%dims)) :: blk_size

            ind_blk(:) = split_tensor_index(rowcol, t%dims)
            blk_size(:) = get_array_elements(t%blk_size, ind_blk)
            tas_blk_size_t = PRODUCT(blk_size)

         END FUNCTION

! **************************************************************************************************
!> \brief Create a default nd process topology that is consistent with a given 2d topology.
!>        Purpose: a nd tensor defined on the returned process grid can be represented as a DBM
!>        matrix with the given 2d topology.
!>        This is needed to enable contraction of 2 tensors (must have the same 2d process grid).
!> \param comm_2d communicator with 2-dimensional topology
!> \param map1_2d which nd-indices map to first matrix index and in which order
!> \param map2_2d which nd-indices map to second matrix index and in which order
!> \param dims_nd nd dimensions
!> \param pdims_2d if comm_2d does not have a cartesian topology associated, can input dimensions
!>                 with pdims_2d
!> \param tdims tensor block dimensions. If present, process grid dimensions are created such that
!>              good load balancing is ensured even if some of the tensor dimensions are small
!>              (i.e. on the same order or smaller than nproc**(1/ndim) where ndim is the tensor rank)
!> \return with nd cartesian grid
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION dbt_nd_mp_comm(comm_2d, map1_2d, map2_2d, dims_nd, dims1_nd, dims2_nd, pdims_2d, tdims, &
                                 nsplit, dimsplit)
            INTEGER, INTENT(IN)                               :: comm_2d
            INTEGER, DIMENSION(:), INTENT(IN)                 :: map1_2d, map2_2d
            INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)), &
               INTENT(IN), OPTIONAL                           :: dims_nd
            INTEGER, DIMENSION(SIZE(map1_2d)), INTENT(IN), OPTIONAL :: dims1_nd
            INTEGER, DIMENSION(SIZE(map2_2d)), INTENT(IN), OPTIONAL :: dims2_nd
            INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL           :: pdims_2d
            INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)), &
               INTENT(IN), OPTIONAL                           :: tdims
            INTEGER, INTENT(IN), OPTIONAL :: nsplit, dimsplit
            INTEGER                                           :: ndim1, ndim2
            INTEGER                                           :: numtask
            INTEGER, DIMENSION(2)                             :: dims_2d, task_coor

            INTEGER, DIMENSION(SIZE(map1_2d)) :: dims1_nd_prv
            INTEGER, DIMENSION(SIZE(map2_2d)) :: dims2_nd_prv
            INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: dims_nd_prv
            INTEGER                                           :: handle
            CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_nd_mp_comm'
            TYPE(dbt_pgrid_type)                          :: dbt_nd_mp_comm

            CALL timeset(routineN, handle)

            ndim1 = SIZE(map1_2d); ndim2 = SIZE(map2_2d)

            IF (PRESENT(pdims_2d)) THEN
               dims_2d(:) = pdims_2d
            ELSE
               CALL mp_environ(numtask, dims_2d, task_coor, comm_2d)
            END IF

            IF (.NOT. PRESENT(dims_nd)) THEN
               dims1_nd_prv = 0; dims2_nd_prv = 0
               IF (PRESENT(dims1_nd)) THEN
                  dims1_nd_prv(:) = dims1_nd
               ELSE

                  IF (PRESENT(tdims)) THEN
                     CALL dbt_mp_dims_create(dims_2d(1), dims1_nd_prv, tdims(map1_2d))
                  ELSE
                     CALL mp_dims_create(dims_2d(1), dims1_nd_prv)
                  END IF
               END IF

               IF (PRESENT(dims2_nd)) THEN
                  dims2_nd_prv(:) = dims2_nd
               ELSE
                  IF (PRESENT(tdims)) THEN
                     CALL dbt_mp_dims_create(dims_2d(2), dims2_nd_prv, tdims(map2_2d))
                  ELSE
                     CALL mp_dims_create(dims_2d(2), dims2_nd_prv)
                  END IF
               END IF
               dims_nd_prv(map1_2d) = dims1_nd_prv
               dims_nd_prv(map2_2d) = dims2_nd_prv
            ELSE
               CPASSERT(PRODUCT(dims_nd(map1_2d)) == dims_2d(1))
               CPASSERT(PRODUCT(dims_nd(map2_2d)) == dims_2d(2))
               dims_nd_prv = dims_nd
            END IF

            CALL dbt_pgrid_create_expert(comm_2d, dims_nd_prv, dbt_nd_mp_comm, &
                                         tensor_dims=tdims, map1_2d=map1_2d, map2_2d=map2_2d, nsplit=nsplit, dimsplit=dimsplit)

            CALL timestop(handle)

         END FUNCTION

! **************************************************************************************************
!> \brief Create process grid dimensions corresponding to one dimension of the matrix representation
!>        of a tensor, imposing that no process grid dimension is greater than the corresponding
!>        tensor dimension.
!> \param nodes Total number of nodes available for this matrix dimension
!> \param dims process grid dimension corresponding to tensor_dims
!> \param tensor_dims tensor dimensions
!> \param lb_ratio load imbalance acceptance factor
!> \author Patrick Seewald
! **************************************************************************************************
         RECURSIVE SUBROUTINE dbt_mp_dims_create(nodes, dims, tensor_dims, lb_ratio)
            INTEGER, INTENT(IN) :: nodes
            INTEGER, DIMENSION(:), INTENT(INOUT) :: dims
            INTEGER, DIMENSION(:), INTENT(IN) :: tensor_dims
            REAL(dp), INTENT(IN), OPTIONAL :: lb_ratio

            INTEGER, DIMENSION(:), ALLOCATABLE :: tensor_dims_sorted, sort_indices, dims_store
            REAL(dp), DIMENSION(:), ALLOCATABLE :: sort_key
            INTEGER :: pdims_rem, idim, pdim
            REAL(dp) :: lb_ratio_prv

            IF (PRESENT(lb_ratio)) THEN
               lb_ratio_prv = lb_ratio
            ELSE
               lb_ratio_prv = 0.1_dp
            END IF

            ALLOCATE (dims_store, source=dims)

            ! get default process grid dimensions
            IF (any(dims == 0)) THEN
               CALL mp_dims_create(nodes, dims)
            END IF

            ! sort dimensions such that problematic grid dimensions (those who should be corrected) come first
            ALLOCATE (sort_key(SIZE(tensor_dims)))
            sort_key(:) = REAL(tensor_dims, dp)/dims

            ALLOCATE (tensor_dims_sorted, source=tensor_dims)
            ALLOCATE (sort_indices(SIZE(sort_key)))
            CALL sort(sort_key, SIZE(sort_key), sort_indices)
            tensor_dims_sorted(:) = tensor_dims_sorted(sort_indices)
            dims(:) = dims(sort_indices)

            ! remaining number of nodes
            pdims_rem = nodes

            DO idim = 1, SIZE(tensor_dims_sorted)
               IF (.NOT. accept_pdims_loadbalancing(pdims_rem, dims(idim), tensor_dims_sorted(idim), lb_ratio_prv)) THEN
                  pdim = tensor_dims_sorted(idim)
                  DO WHILE (.NOT. accept_pdims_loadbalancing(pdims_rem, pdim, tensor_dims_sorted(idim), lb_ratio_prv))
                     pdim = pdim - 1
                  END DO
                  dims(idim) = pdim
                  pdims_rem = pdims_rem/dims(idim)

                  IF (idim .NE. SIZE(tensor_dims_sorted)) THEN
                     dims(idim + 1:) = 0
                     CALL mp_dims_create(pdims_rem, dims(idim + 1:))
                  ELSEIF (lb_ratio_prv < 0.5_dp) THEN
                     ! resort to a less strict load imbalance factor
                     dims(:) = dims_store
                     CALL dbt_mp_dims_create(nodes, dims, tensor_dims, 0.5_dp)
                     RETURN
                  ELSE
                     ! resort to default process grid dimensions
                     dims(:) = dims_store
                     CALL mp_dims_create(nodes, dims)
                     RETURN
                  END IF

               ELSE
                  pdims_rem = pdims_rem/dims(idim)
               END IF
            END DO

            dims(sort_indices) = dims

         END SUBROUTINE

! **************************************************************************************************
!> \brief load balancing criterion whether to accept process grid dimension based on total number of
!>        cores and tensor dimension
!> \param pdims_avail available process grid dimensions (total number of cores)
!> \param pdim process grid dimension to test
!> \param tdim tensor dimension corresponding to pdim
!> \param lb_ratio load imbalance acceptance factor
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION accept_pdims_loadbalancing(pdims_avail, pdim, tdim, lb_ratio)
            INTEGER, INTENT(IN) :: pdims_avail
            INTEGER, INTENT(IN) :: pdim
            INTEGER, INTENT(IN) :: tdim
            REAL(dp), INTENT(IN) :: lb_ratio
            LOGICAL :: accept_pdims_loadbalancing

            accept_pdims_loadbalancing = .FALSE.
            IF (MOD(pdims_avail, pdim) == 0) THEN
               IF (REAL(tdim, dp)*lb_ratio < REAL(pdim, dp)) THEN
                  IF (MOD(tdim, pdim) == 0) accept_pdims_loadbalancing = .TRUE.
               ELSE
                  accept_pdims_loadbalancing = .TRUE.
               END IF
            END IF

         END FUNCTION

! **************************************************************************************************
!> \brief Release the MPI communicator.
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_nd_mp_free(mp_comm)
            INTEGER, INTENT(INOUT)                               :: mp_comm

            CALL mp_comm_free(mp_comm)
         END SUBROUTINE dbt_nd_mp_free

! **************************************************************************************************
!> \brief Create a tensor distribution.
!> \param pgrid process grid
!> \param nd_dist_i distribution vectors for all tensor dimensions
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_distribution_new(dist, pgrid, ${varlist("nd_dist")}$)
            TYPE(dbt_distribution_type), INTENT(OUT)    :: dist
            TYPE(dbt_pgrid_type), INTENT(IN)            :: pgrid
            INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL     :: ${varlist("nd_dist")}$
            INTEGER, DIMENSION(ndims_mapping_row(pgrid%nd_index_grid)) :: map1_2d
            INTEGER, DIMENSION(ndims_mapping_column(pgrid%nd_index_grid)) :: map2_2d
            INTEGER :: ndims

            CALL dbt_get_mapping_info(pgrid%nd_index_grid, map1_2d=map1_2d, map2_2d=map2_2d, ndim_nd=ndims)

            CALL dbt_distribution_new_expert(dist, pgrid, map1_2d, map2_2d, ${varlist("nd_dist")}$)

         END SUBROUTINE

! **************************************************************************************************
!> \brief Create a tensor distribution.
!> \param pgrid process grid
!> \param map1_2d which nd-indices map to first matrix index and in which order
!> \param map2_2d which nd-indices map to second matrix index and in which order
!> \param own_comm whether distribution should own communicator
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_distribution_new_expert(dist, pgrid, map1_2d, map2_2d, ${varlist("nd_dist")}$, own_comm)
            TYPE(dbt_distribution_type), INTENT(OUT)    :: dist
            TYPE(dbt_pgrid_type), INTENT(IN)            :: pgrid
            INTEGER, DIMENSION(:), INTENT(IN)               :: map1_2d
            INTEGER, DIMENSION(:), INTENT(IN)               :: map2_2d
            INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL     :: ${varlist("nd_dist")}$
            LOGICAL, INTENT(IN), OPTIONAL                   :: own_comm
            INTEGER                                         :: ndims, comm_2d
            INTEGER, DIMENSION(2)                           :: pdims_2d_check, &
                                                               pdims_2d, task_coor_2d
            INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: dims, nblks_nd, task_coor
            LOGICAL, DIMENSION(2)                           :: periods_2d
            TYPE(array_list)                                :: nd_dist
            TYPE(nd_to_2d_mapping)                          :: map_blks, map_grid
            INTEGER                                         :: handle
            TYPE(dbt_tas_dist_t)                          :: row_dist_obj, col_dist_obj
            TYPE(dbt_pgrid_type)                        :: pgrid_prv
            LOGICAL                                         :: need_pgrid_remap
            INTEGER, DIMENSION(ndims_mapping_row(pgrid%nd_index_grid)) :: map1_2d_check
            INTEGER, DIMENSION(ndims_mapping_column(pgrid%nd_index_grid)) :: map2_2d_check
            CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_distribution_new'

            CALL timeset(routineN, handle)
            ndims = SIZE(map1_2d) + SIZE(map2_2d)
            CPASSERT(ndims .GE. 2 .AND. ndims .LE. ${maxdim}$)

            CALL create_array_list(nd_dist, ndims, ${varlist("nd_dist")}$)

            nblks_nd(:) = sizes_of_arrays(nd_dist)

            need_pgrid_remap = .TRUE.
            IF (PRESENT(own_comm)) THEN
               CALL dbt_get_mapping_info(pgrid%nd_index_grid, map1_2d=map1_2d_check, map2_2d=map2_2d_check)
               IF (own_comm) THEN
                  IF (.NOT. array_eq_i(map1_2d_check, map1_2d) .OR. .NOT. array_eq_i(map2_2d_check, map2_2d)) THEN
                     CPABORT("map1_2d / map2_2d are not consistent with pgrid")
                  END IF
                  pgrid_prv = pgrid
                  need_pgrid_remap = .FALSE.
               END IF
            END IF

            IF (need_pgrid_remap) CALL dbt_pgrid_remap(pgrid, map1_2d, map2_2d, pgrid_prv)

            ! check that 2d process topology is consistent with nd topology.
            CALL mp_environ_pgrid(pgrid_prv, dims, task_coor)

            ! process grid index mapping
            CALL create_nd_to_2d_mapping(map_grid, dims, map1_2d, map2_2d, base=0, col_major=.FALSE.)

            ! blk index mapping
            CALL create_nd_to_2d_mapping(map_blks, nblks_nd, map1_2d, map2_2d)

            row_dist_obj = dbt_tas_dist_t(nd_dist, map_blks, map_grid, 1)
            col_dist_obj = dbt_tas_dist_t(nd_dist, map_blks, map_grid, 2)

            CALL dbt_get_mapping_info(map_grid, dims_2d=pdims_2d)

            comm_2d = pgrid_prv%mp_comm_2d

            CALL mp_environ(comm_2d, 2, pdims_2d_check, task_coor_2d, periods_2d)
            IF (ANY(pdims_2d_check .NE. pdims_2d)) THEN
               CPABORT("inconsistent process grid dimensions")
            END IF

            IF (ALLOCATED(pgrid_prv%tas_split_info)) THEN
               CALL dbt_tas_distribution_new(dist%dist, comm_2d, row_dist_obj, col_dist_obj, split_info=pgrid_prv%tas_split_info)
            ELSE
               CALL dbt_tas_distribution_new(dist%dist, comm_2d, row_dist_obj, col_dist_obj)
               ALLOCATE (pgrid_prv%tas_split_info, SOURCE=dist%dist%info)
               CALL dbt_tas_info_hold(pgrid_prv%tas_split_info)
            END IF

            dist%nd_dist = nd_dist
            dist%pgrid = pgrid_prv

            ALLOCATE (dist%refcount)
            dist%refcount = 1
            CALL timestop(handle)

         CONTAINS
            PURE FUNCTION array_eq_i(arr1, arr2)
               INTEGER, INTENT(IN), DIMENSION(:) :: arr1
               INTEGER, INTENT(IN), DIMENSION(:) :: arr2
               LOGICAL                           :: array_eq_i

               array_eq_i = .FALSE.
               IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)

            END FUNCTION

         END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy tensor distribution
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_distribution_destroy(dist)
            TYPE(dbt_distribution_type), INTENT(INOUT) :: dist
            INTEGER                                   :: handle
            CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_distribution_destroy'
            LOGICAL :: abort

            CALL timeset(routineN, handle)
            CALL dbt_tas_distribution_destroy(dist%dist)
            CALL destroy_array_list(dist%nd_dist)

            abort = .FALSE.
            IF (.NOT. ASSOCIATED(dist%refcount)) THEN
               abort = .TRUE.
            ELSEIF (dist%refcount < 1) THEN
               abort = .TRUE.
            END IF

            IF (abort) THEN
               CPABORT("can not destroy non-existing tensor distribution")
            END IF

            dist%refcount = dist%refcount - 1

            IF (dist%refcount == 0) THEN
               CALL dbt_pgrid_destroy(dist%pgrid)
               DEALLOCATE (dist%refcount)
            ELSE
               CALL dbt_pgrid_destroy(dist%pgrid, keep_comm=.TRUE.)
            END IF

            CALL timestop(handle)
         END SUBROUTINE

! **************************************************************************************************
!> \brief reference counting for distribution
!>        (only needed for communicator handle that must be freed when no longer needed)
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_distribution_hold(dist)
            TYPE(dbt_distribution_type), INTENT(IN) :: dist
            INTEGER, POINTER                            :: ref

            IF (dist%refcount < 1) THEN
               CPABORT("can not hold non-existing tensor distribution")
            END IF
            ref => dist%refcount
            ref = ref + 1
         END SUBROUTINE

! **************************************************************************************************
!> \brief get distribution from tensor
!> \return distribution
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION dbt_distribution(tensor)
            TYPE(dbt_type), INTENT(IN)  :: tensor
            TYPE(dbt_distribution_type) :: dbt_distribution

            CALL dbt_tas_get_info(tensor%matrix_rep, distribution=dbt_distribution%dist)
            dbt_distribution%pgrid = tensor%pgrid
            dbt_distribution%nd_dist = tensor%nd_dist
            dbt_distribution%refcount => dbt_distribution%refcount
         END FUNCTION

! **************************************************************************************************
!> \brief create a tensor.
!>        For performance, the arguments map1_2d and map2_2d (controlling matrix representation of
!>        tensor) should be consistent with the the contraction to be performed (see documentation
!>        of dbt_contract).
!> \param map1_2d which nd-indices to map to first 2d index and in which order
!> \param map2_2d which nd-indices to map to first 2d index and in which order
!> \param blk_size_i blk sizes in each dimension
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_create_new(tensor, name, dist, map1_2d, map2_2d, &
                                   ${varlist("blk_size")}$)
            TYPE(dbt_type), INTENT(OUT)                   :: tensor
            CHARACTER(len=*), INTENT(IN)                      :: name
            TYPE(dbt_distribution_type), INTENT(INOUT)    :: dist
            INTEGER, DIMENSION(:), INTENT(IN)                 :: map1_2d
            INTEGER, DIMENSION(:), INTENT(IN)                 :: map2_2d
            INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL       :: ${varlist("blk_size")}$
            INTEGER                                           :: ndims
            INTEGER(KIND=int_8), DIMENSION(2)                             :: dims_2d
            INTEGER, DIMENSION(SIZE(map1_2d) + SIZE(map2_2d)) :: dims, pdims, task_coor
            TYPE(dbt_tas_blk_size_t)                        :: col_blk_size_obj, row_blk_size_obj
            TYPE(dbt_distribution_type)                   :: dist_new
            TYPE(array_list)                                  :: blk_size, blks_local
            TYPE(nd_to_2d_mapping)                            :: map
            INTEGER                                   :: handle
            CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_create_new'
            INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("blks_local")}$
            INTEGER, DIMENSION(:), ALLOCATABLE              :: ${varlist("dist")}$
            INTEGER                                         :: iblk_count, iblk
            INTEGER, DIMENSION(:), ALLOCATABLE              :: nblks_local, nfull_local

            CALL timeset(routineN, handle)
            ndims = SIZE(map1_2d) + SIZE(map2_2d)
            CALL create_array_list(blk_size, ndims, ${varlist("blk_size")}$)
            dims = sizes_of_arrays(blk_size)

            CALL create_nd_to_2d_mapping(map, dims, map1_2d, map2_2d)
            CALL dbt_get_mapping_info(map, dims_2d_i8=dims_2d)

            row_blk_size_obj = dbt_tas_blk_size_t(blk_size, map, 1)
            col_blk_size_obj = dbt_tas_blk_size_t(blk_size, map, 2)

            CALL dbt_distribution_remap(dist, map1_2d, map2_2d, dist_new)

            ALLOCATE (tensor%matrix_rep)
            CALL dbt_tas_create(matrix=tensor%matrix_rep, &
                                name=TRIM(name)//" matrix", &
                                dist=dist_new%dist, &
                                row_blk_size=row_blk_size_obj, &
                                col_blk_size=col_blk_size_obj)

            tensor%owns_matrix = .TRUE.

            tensor%nd_index_blk = map
            tensor%name = name

            CALL dbt_tas_finalize(tensor%matrix_rep)
            CALL destroy_nd_to_2d_mapping(map)

            ! map element-wise tensor index
            CALL create_nd_to_2d_mapping(map, sum_of_arrays(blk_size), map1_2d, map2_2d)
            tensor%nd_index = map
            tensor%blk_sizes = blk_size

            CALL mp_environ_pgrid(dist_new%pgrid, pdims, task_coor)

            #:for ndim in range(1, maxdim+1)
               IF (ndims == ${ndim}$) THEN
                  CALL get_arrays(dist_new%nd_dist, ${varlist("dist", nmax=ndim)}$)
               END IF
            #:endfor

            ALLOCATE (nblks_local(ndims))
            ALLOCATE (nfull_local(ndims))
            nfull_local(:) = 0
            #:for idim in range(1, maxdim+1)
               IF (ndims .GE. ${idim}$) THEN
                  nblks_local(${idim}$) = COUNT(dist_${idim}$ == task_coor(${idim}$))
                  ALLOCATE (blks_local_${idim}$ (nblks_local(${idim}$)))
                  iblk_count = 0
                  DO iblk = 1, SIZE(dist_${idim}$)
                     IF (dist_${idim}$ (iblk) == task_coor(${idim}$)) THEN
                        iblk_count = iblk_count + 1
                        blks_local_${idim}$ (iblk_count) = iblk
                        nfull_local(${idim}$) = nfull_local(${idim}$) + blk_size_${idim}$ (iblk)
                     END IF
                  END DO
               END IF
            #:endfor

            #:for ndim in range(1, maxdim+1)
               IF (ndims == ${ndim}$) THEN
                  CALL create_array_list(blks_local, ${ndim}$, ${varlist("blks_local", nmax=ndim)}$)
               END IF
            #:endfor

            ALLOCATE (tensor%nblks_local(ndims))
            ALLOCATE (tensor%nfull_local(ndims))
            tensor%nblks_local(:) = nblks_local
            tensor%nfull_local(:) = nfull_local

            tensor%blks_local = blks_local

            tensor%nd_dist = dist_new%nd_dist
            tensor%pgrid = dist_new%pgrid

            CALL dbt_distribution_hold(dist_new)
            tensor%refcount => dist_new%refcount
            CALL dbt_distribution_destroy(dist_new)

            CALL array_offsets(tensor%blk_sizes, tensor%blk_offsets)

            tensor%valid = .TRUE.
            CALL timestop(handle)
         END SUBROUTINE

! **************************************************************************************************
!> \brief reference counting for tensors
!>        (only needed for communicator handle that must be freed when no longer needed)
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_hold(tensor)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER, POINTER :: ref

            IF (tensor%refcount < 1) THEN
               CPABORT("can not hold non-existing tensor")
            END IF
            ref => tensor%refcount
            ref = ref + 1

         END SUBROUTINE

! **************************************************************************************************
!> \brief create a tensor from template
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_create_template(tensor_in, tensor, name, dist, map1_2d, map2_2d)
            TYPE(dbt_type), INTENT(INOUT)      :: tensor_in
            TYPE(dbt_type), INTENT(OUT)        :: tensor
            CHARACTER(len=*), INTENT(IN), OPTIONAL :: name
            TYPE(dbt_distribution_type), &
               INTENT(INOUT), OPTIONAL             :: dist
            INTEGER, DIMENSION(:), INTENT(IN), &
               OPTIONAL                            :: map1_2d, map2_2d
            INTEGER                                :: handle
            CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_create_template'
            INTEGER, DIMENSION(:), ALLOCATABLE     :: ${varlist("bsize")}$
            INTEGER, DIMENSION(:), ALLOCATABLE     :: map1_2d_prv, map2_2d_prv
            CHARACTER(len=default_string_length)   :: name_prv
            TYPE(dbt_distribution_type)        :: dist_prv

            CALL timeset(routineN, handle)

            IF (PRESENT(dist) .OR. PRESENT(map1_2d) .OR. PRESENT(map2_2d)) THEN
               ! need to create matrix representation from scratch
               IF (PRESENT(dist)) THEN
                  dist_prv = dist
               ELSE
                  dist_prv = dbt_distribution(tensor_in)
               END IF
               IF (PRESENT(map1_2d) .AND. PRESENT(map2_2d)) THEN
                  ALLOCATE (map1_2d_prv, source=map1_2d)
                  ALLOCATE (map2_2d_prv, source=map2_2d)
               ELSE
                  ALLOCATE (map1_2d_prv(ndims_matrix_row(tensor_in)))
                  ALLOCATE (map2_2d_prv(ndims_matrix_column(tensor_in)))
                  CALL dbt_get_mapping_info(tensor_in%nd_index_blk, map1_2d=map1_2d_prv, map2_2d=map2_2d_prv)
               END IF
               IF (PRESENT(name)) THEN
                  name_prv = name
               ELSE
                  name_prv = tensor_in%name
               END IF

               #:for ndim in range(1, maxdim+1)
                  IF (ndims_tensor(tensor_in) == ${ndim}$) THEN
                     CALL get_arrays(tensor_in%blk_sizes, ${varlist("bsize", nmax=ndim)}$)
                     CALL dbt_create(tensor, name_prv, dist_prv, map1_2d_prv, map2_2d_prv, &
                                     ${varlist("bsize", nmax=ndim)}$)
                  END IF
               #:endfor
            ELSE
               ! create matrix representation from template
               ALLOCATE (tensor%matrix_rep)
               IF (.NOT. PRESENT(name)) THEN
                  CALL dbt_tas_create(tensor_in%matrix_rep, tensor%matrix_rep, &
                                      name=TRIM(tensor_in%name)//" matrix")
               ELSE
                  CALL dbt_tas_create(tensor_in%matrix_rep, tensor%matrix_rep, name=TRIM(name)//" matrix")
               END IF
               tensor%owns_matrix = .TRUE.
               CALL dbt_tas_finalize(tensor%matrix_rep)

               tensor%nd_index_blk = tensor_in%nd_index_blk
               tensor%nd_index = tensor_in%nd_index
               tensor%blk_sizes = tensor_in%blk_sizes
               tensor%blk_offsets = tensor_in%blk_offsets
               tensor%nd_dist = tensor_in%nd_dist
               tensor%blks_local = tensor_in%blks_local
               ALLOCATE (tensor%nblks_local(ndims_tensor(tensor_in)))
               tensor%nblks_local(:) = tensor_in%nblks_local
               ALLOCATE (tensor%nfull_local(ndims_tensor(tensor_in)))
               tensor%nfull_local(:) = tensor_in%nfull_local
               tensor%pgrid = tensor_in%pgrid

               tensor%refcount => tensor_in%refcount
               CALL dbt_hold(tensor)

               tensor%valid = .TRUE.
               IF (PRESENT(name)) THEN
                  tensor%name = name
               ELSE
                  tensor%name = tensor_in%name
               END IF
            END IF
            CALL timestop(handle)
         END SUBROUTINE

! **************************************************************************************************
!> \brief Create 2-rank tensor from matrix.
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_create_matrix(matrix_in, tensor, order, name)
            TYPE(dbcsr_type), INTENT(IN)                :: matrix_in
            TYPE(dbt_type), INTENT(OUT)        :: tensor
            INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: order
            CHARACTER(len=*), INTENT(IN), OPTIONAL      :: name

            CHARACTER(len=default_string_length)        :: name_in
            INTEGER, DIMENSION(2)                       :: order_in
            INTEGER                                     :: comm_2d
            TYPE(dbcsr_distribution_type)               :: matrix_dist
            TYPE(dbt_distribution_type)        :: dist
            INTEGER, DIMENSION(:), POINTER              :: row_blk_size, col_blk_size
            INTEGER, DIMENSION(:), POINTER              :: col_dist, row_dist
            INTEGER                                   :: handle
            CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_create_matrix'
            TYPE(dbt_pgrid_type)                  :: comm_nd
            INTEGER, DIMENSION(2)                     :: pdims_2d

            CALL timeset(routineN, handle)

            NULLIFY (row_blk_size, col_blk_size, col_dist, row_dist)
            IF (PRESENT(name)) THEN
               name_in = name
            ELSE
               CALL dbcsr_get_info(matrix_in, name=name_in)
            END IF

            IF (PRESENT(order)) THEN
               order_in = order
            ELSE
               order_in = [1, 2]
            END IF

            CALL dbcsr_get_info(matrix_in, distribution=matrix_dist)
            CALL dbcsr_distribution_get(matrix_dist, group=comm_2d, row_dist=row_dist, col_dist=col_dist, &
                                        nprows=pdims_2d(1), npcols=pdims_2d(2))
            comm_nd = dbt_nd_mp_comm(comm_2d, [order_in(1)], [order_in(2)], pdims_2d=pdims_2d)

            CALL dbt_distribution_new_expert( &
               dist, &
               comm_nd, &
               [order_in(1)], [order_in(2)], &
               row_dist, col_dist, own_comm=.TRUE.)

            CALL dbcsr_get_info(matrix_in, row_blk_size=row_blk_size, col_blk_size=col_blk_size)

            CALL dbt_create_new(tensor, name_in, dist, &
                                [order_in(1)], [order_in(2)], &
                                row_blk_size, &
                                col_blk_size)

            CALL dbt_distribution_destroy(dist)
            CALL timestop(handle)
         END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy a tensor
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_destroy(tensor)
            TYPE(dbt_type), INTENT(INOUT) :: tensor
            INTEGER                                   :: handle
            CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_destroy'
            LOGICAL :: abort

            CALL timeset(routineN, handle)
            IF (tensor%owns_matrix) THEN
               CALL dbt_tas_destroy(tensor%matrix_rep)
               DEALLOCATE (tensor%matrix_rep)
            ELSE
               NULLIFY (tensor%matrix_rep)
            END IF
            tensor%owns_matrix = .FALSE.

            CALL destroy_nd_to_2d_mapping(tensor%nd_index_blk)
            CALL destroy_nd_to_2d_mapping(tensor%nd_index)
            !CALL destroy_nd_to_2d_mapping(tensor%nd_index_grid)
            CALL destroy_array_list(tensor%blk_sizes)
            CALL destroy_array_list(tensor%blk_offsets)
            CALL destroy_array_list(tensor%nd_dist)
            CALL destroy_array_list(tensor%blks_local)

            DEALLOCATE (tensor%nblks_local, tensor%nfull_local)

            abort = .FALSE.
            IF (.NOT. ASSOCIATED(tensor%refcount)) THEN
               abort = .TRUE.
            ELSEIF (tensor%refcount < 1) THEN
               abort = .TRUE.
            END IF

            IF (abort) THEN
               CPABORT("can not destroy non-existing tensor")
            END IF

            tensor%refcount = tensor%refcount - 1

            IF (tensor%refcount == 0) THEN
               CALL dbt_pgrid_destroy(tensor%pgrid)
               !CALL mp_comm_free(tensor%comm_2d)
               !CALL mp_comm_free(tensor%comm_nd)
               DEALLOCATE (tensor%refcount)
            ELSE
               CALL dbt_pgrid_destroy(tensor%pgrid, keep_comm=.TRUE.)
            END IF

            tensor%valid = .FALSE.
            tensor%name = ""
            CALL timestop(handle)
         END SUBROUTINE

! **************************************************************************************************
!> \brief tensor rank
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION ndims_tensor(tensor)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER                        :: ndims_tensor

            ndims_tensor = tensor%nd_index%ndim_nd
         END FUNCTION

! **************************************************************************************************
!> \brief tensor dimensions
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dims_tensor(tensor, dims)
            TYPE(dbt_type), INTENT(IN)              :: tensor
            INTEGER, DIMENSION(ndims_tensor(tensor)), &
               INTENT(OUT)                              :: dims

            CPASSERT(tensor%valid)
            dims = tensor%nd_index%dims_nd
         END SUBROUTINE

! **************************************************************************************************
!> \brief tensor block dimensions
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE blk_dims_tensor(tensor, dims)
            TYPE(dbt_type), INTENT(IN)              :: tensor
            INTEGER, DIMENSION(ndims_tensor(tensor)), &
               INTENT(OUT)                              :: dims

            CPASSERT(tensor%valid)
            dims = tensor%nd_index_blk%dims_nd
         END SUBROUTINE

! **************************************************************************************************
!> \brief Size of tensor block
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_blk_sizes(tensor, ind, blk_size)
            TYPE(dbt_type), INTENT(IN)              :: tensor
            INTEGER, DIMENSION(ndims_tensor(tensor)), &
               INTENT(IN)                               :: ind
            INTEGER, DIMENSION(ndims_tensor(tensor)), &
               INTENT(OUT)                              :: blk_size

            blk_size(:) = get_array_elements(tensor%blk_sizes, ind)
         END SUBROUTINE

! **************************************************************************************************
!> \brief offset of tensor block
!> \param ind block index
!> \param blk_offset block offset
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_blk_offsets(tensor, ind, blk_offset)
            TYPE(dbt_type), INTENT(IN)              :: tensor
            INTEGER, DIMENSION(ndims_tensor(tensor)), &
               INTENT(IN)                               :: ind
            INTEGER, DIMENSION(ndims_tensor(tensor)), &
               INTENT(OUT)                              :: blk_offset

            CPASSERT(tensor%valid)
            blk_offset(:) = get_array_elements(tensor%blk_offsets, ind)
         END SUBROUTINE

! **************************************************************************************************
!> \brief Generalization of block_get_stored_coordinates for tensors.
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_get_stored_coordinates(tensor, ind_nd, processor)
            TYPE(dbt_type), INTENT(IN)               :: tensor
            INTEGER, DIMENSION(ndims_tensor(tensor)), &
               INTENT(IN)                                :: ind_nd
            INTEGER, INTENT(OUT)                         :: processor

            INTEGER(KIND=int_8), DIMENSION(2)                        :: ind_2d

            ind_2d(:) = get_2d_indices_tensor(tensor%nd_index_blk, ind_nd)
            CALL dbt_tas_get_stored_coordinates(tensor%matrix_rep, ind_2d(1), ind_2d(2), processor)
         END SUBROUTINE

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_pgrid_create(mp_comm, dims, pgrid, tensor_dims)
            INTEGER, INTENT(IN) :: mp_comm
            INTEGER, DIMENSION(:), INTENT(INOUT) :: dims
            TYPE(dbt_pgrid_type), INTENT(OUT) :: pgrid
            INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: tensor_dims
            INTEGER, DIMENSION(:), ALLOCATABLE :: map1_2d, map2_2d
            INTEGER :: i, ndims

            ndims = SIZE(dims)

            ALLOCATE (map1_2d(ndims/2))
            ALLOCATE (map2_2d(ndims - ndims/2))
            map1_2d(:) = (/(i, i=1, SIZE(map1_2d))/)
            map2_2d(:) = (/(i, i=SIZE(map1_2d) + 1, SIZE(map1_2d) + SIZE(map2_2d))/)

            CALL dbt_pgrid_create_expert(mp_comm, dims, pgrid, map1_2d, map2_2d, tensor_dims)

         END SUBROUTINE

! **************************************************************************************************
!> \brief Create an n-dimensional process grid.
!>        We can not use a n-dimensional MPI cartesian grid for tensors since the mapping between
!>        n-dim. and 2-dim. index allows for an arbitrary reordering of tensor index. Therefore we
!>        can not use n-dim. MPI Cartesian grid because it may not be consistent with the respective
!>        2d grid. The 2d Cartesian MPI grid is the reference grid (since tensor data is stored as
!>        DBM matrix) and this routine creates an object that is a n-dim. interface to this grid.
!>        map1_2d and map2_2d don't need to be specified (correctly), grid may be redefined in
!>        dbt_distribution_new. Note that pgrid is equivalent to a MPI cartesian grid only
!>        if map1_2d and map2_2d don't reorder indices (which is the case if
!>        [map1_2d, map2_2d] == [1, 2, ..., ndims]). Otherwise the mapping of grid coordinates to
!>        processes depends on the ordering of the indices and is not equivalent to a MPI cartesian
!>        grid.
!> \param mp_comm simple MPI Communicator
!> \param dims grid dimensions - if entries are 0, dimensions are chosen automatically.
!> \param pgrid n-dimensional grid object
!> \param map1_2d which nd-indices map to first matrix index and in which order
!> \param map2_2d which nd-indices map to first matrix index and in which order
!> \param tensor_dims tensor block dimensions. If present, process grid dimensions are created such
!>                    that good load balancing is ensured even if some of the tensor dimensions are
!>                    small (i.e. on the same order or smaller than nproc**(1/ndim) where ndim is
!>                    the tensor rank)
!> \param nsplit impose a constant split factor
!> \param dimsplit which matrix dimension to split
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_pgrid_create_expert(mp_comm, dims, pgrid, map1_2d, map2_2d, tensor_dims, nsplit, dimsplit)
            INTEGER, INTENT(IN) :: mp_comm
            INTEGER, DIMENSION(:), INTENT(INOUT) :: dims
            TYPE(dbt_pgrid_type), INTENT(OUT) :: pgrid
            INTEGER, DIMENSION(:), INTENT(IN) :: map1_2d, map2_2d
            INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: tensor_dims
            INTEGER, INTENT(IN), OPTIONAL :: nsplit, dimsplit
            INTEGER :: nproc, iproc, ndims, handle
            INTEGER, DIMENSION(2) :: pdims_2d, pos
            TYPE(dbt_tas_split_info) :: info

            CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_pgrid_create'

            CALL timeset(routineN, handle)

            ndims = SIZE(dims)

            CALL mp_environ(nproc, iproc, mp_comm)
            IF (ANY(dims == 0)) THEN
               IF (.NOT. PRESENT(tensor_dims)) THEN
                  CALL mp_dims_create(nproc, dims)
               ELSE
                  CALL dbt_mp_dims_create(nproc, dims, tensor_dims)
               END IF
            END IF
            CALL create_nd_to_2d_mapping(pgrid%nd_index_grid, dims, map1_2d, map2_2d, base=0, col_major=.FALSE.)
            CALL dbt_get_mapping_info(pgrid%nd_index_grid, dims_2d=pdims_2d)
            CALL mp_cart_create(mp_comm, 2, pdims_2d, pos, pgrid%mp_comm_2d)

            IF (PRESENT(nsplit)) THEN
               CPASSERT(PRESENT(dimsplit))
               CALL dbt_tas_create_split(info, pgrid%mp_comm_2d, dimsplit, nsplit, opt_nsplit=.FALSE.)
               ALLOCATE (pgrid%tas_split_info, SOURCE=info)
            END IF

            ! store number of MPI ranks because we need it for PURE function dbt_max_nblks_local
            pgrid%nproc = nproc

            CALL timestop(handle)
         END SUBROUTINE

! **************************************************************************************************
!> \brief destroy process grid
!> \param keep_comm  if .TRUE. communicator is not freed
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_pgrid_destroy(pgrid, keep_comm)
            TYPE(dbt_pgrid_type), INTENT(INOUT) :: pgrid
            LOGICAL, INTENT(IN), OPTIONAL           :: keep_comm
            LOGICAL :: keep_comm_prv
            IF (PRESENT(keep_comm)) THEN
               keep_comm_prv = keep_comm
            ELSE
               keep_comm_prv = .FALSE.
            END IF
            IF (.NOT. keep_comm_prv) CALL mp_comm_free(pgrid%mp_comm_2d)
            CALL destroy_nd_to_2d_mapping(pgrid%nd_index_grid)
            IF (ALLOCATED(pgrid%tas_split_info) .AND. .NOT. keep_comm_prv) THEN
               CALL dbt_tas_release_info(pgrid%tas_split_info)
               DEALLOCATE (pgrid%tas_split_info)
            END IF
         END SUBROUTINE

! **************************************************************************************************
!> \brief freeze current split factor such that it is never changed during contraction
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_pgrid_set_strict_split(pgrid)
            TYPE(dbt_pgrid_type), INTENT(INOUT) :: pgrid
            IF (ALLOCATED(pgrid%tas_split_info)) CALL dbt_tas_set_strict_split(pgrid%tas_split_info)
         END SUBROUTINE

! **************************************************************************************************
!> \brief remap a process grid (needed when mapping between tensor and matrix index is changed)
!> \param map1_2d new mapping
!> \param map2_2d new mapping
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_pgrid_remap(pgrid_in, map1_2d, map2_2d, pgrid_out)
            TYPE(dbt_pgrid_type), INTENT(IN) :: pgrid_in
            INTEGER, DIMENSION(:), INTENT(IN) :: map1_2d, map2_2d
            TYPE(dbt_pgrid_type), INTENT(OUT) :: pgrid_out
            INTEGER, DIMENSION(:), ALLOCATABLE :: dims
            INTEGER, DIMENSION(ndims_mapping_row(pgrid_in%nd_index_grid)) :: map1_2d_old
            INTEGER, DIMENSION(ndims_mapping_column(pgrid_in%nd_index_grid)) :: map2_2d_old

            ALLOCATE (dims(SIZE(map1_2d) + SIZE(map2_2d)))
            CALL dbt_get_mapping_info(pgrid_in%nd_index_grid, dims_nd=dims, map1_2d=map1_2d_old, map2_2d=map2_2d_old)
            CALL dbt_pgrid_create_expert(pgrid_in%mp_comm_2d, dims, pgrid_out, map1_2d=map1_2d, map2_2d=map2_2d)
            IF (array_eq_i(map1_2d_old, map1_2d) .AND. array_eq_i(map2_2d_old, map2_2d)) THEN
               IF (ALLOCATED(pgrid_in%tas_split_info)) THEN
                  ALLOCATE (pgrid_out%tas_split_info, SOURCE=pgrid_in%tas_split_info)
                  CALL dbt_tas_info_hold(pgrid_out%tas_split_info)
               END IF
            END IF
         END SUBROUTINE

! **************************************************************************************************
!> \brief change dimensions of an existing process grid.
!> \param pgrid process grid to be changed
!> \param pdims new process grid dimensions, should all be set > 0
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_pgrid_change_dims(pgrid, pdims)
            TYPE(dbt_pgrid_type), INTENT(INOUT) :: pgrid
            INTEGER, DIMENSION(:), INTENT(INOUT)    :: pdims
            TYPE(dbt_pgrid_type)                :: pgrid_tmp
            INTEGER                                 :: nsplit, dimsplit
            INTEGER, DIMENSION(ndims_mapping_row(pgrid%nd_index_grid)) :: map1_2d
            INTEGER, DIMENSION(ndims_mapping_column(pgrid%nd_index_grid)) :: map2_2d
            TYPe(nd_to_2d_mapping)                  :: nd_index_grid
            INTEGER, DIMENSION(2)                   :: pdims_2d

            CPASSERT(ALL(pdims > 0))
            CALL dbt_tas_get_split_info(pgrid%tas_split_info, nsplit=nsplit, split_rowcol=dimsplit)
            CALL dbt_get_mapping_info(pgrid%nd_index_grid, map1_2d=map1_2d, map2_2d=map2_2d)
            CALL create_nd_to_2d_mapping(nd_index_grid, pdims, map1_2d, map2_2d, base=0, col_major=.FALSE.)
            CALL dbt_get_mapping_info(nd_index_grid, dims_2d=pdims_2d)
            IF (MOD(pdims_2d(dimsplit), nsplit) == 0) THEN
               CALL dbt_pgrid_create_expert(pgrid%mp_comm_2d, pdims, pgrid_tmp, map1_2d=map1_2d, map2_2d=map2_2d, &
                                            nsplit=nsplit, dimsplit=dimsplit)
            ELSE
               CALL dbt_pgrid_create_expert(pgrid%mp_comm_2d, pdims, pgrid_tmp, map1_2d=map1_2d, map2_2d=map2_2d)
            END IF
            CALL dbt_pgrid_destroy(pgrid)
            pgrid = pgrid_tmp
         END SUBROUTINE

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_distribution_remap(dist_in, map1_2d, map2_2d, dist_out)
            TYPE(dbt_distribution_type), INTENT(IN)    :: dist_in
            INTEGER, DIMENSION(:), INTENT(IN) :: map1_2d, map2_2d
            TYPE(dbt_distribution_type), INTENT(OUT)    :: dist_out
            INTEGER, DIMENSION(:), ALLOCATABLE :: ${varlist("dist")}$
            INTEGER :: ndims
            ndims = SIZE(map1_2d) + SIZE(map2_2d)
            #:for ndim in range(1, maxdim+1)
               IF (ndims == ${ndim}$) THEN
                  CALL get_arrays(dist_in%nd_dist, ${varlist("dist", nmax=ndim)}$)
                  CALL dbt_distribution_new_expert(dist_out, dist_in%pgrid, map1_2d, map2_2d, ${varlist("dist", nmax=ndim)}$)
               END IF
            #:endfor
         END SUBROUTINE

! **************************************************************************************************
!> \brief as mp_environ but for special pgrid type
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE mp_environ_pgrid(pgrid, dims, task_coor)
            TYPE(dbt_pgrid_type), INTENT(IN) :: pgrid
            INTEGER, DIMENSION(ndims_mapping(pgrid%nd_index_grid)), INTENT(OUT) :: dims
            INTEGER, DIMENSION(ndims_mapping(pgrid%nd_index_grid)), INTENT(OUT) :: task_coor
            INTEGER, DIMENSION(2)                                          :: dims_2d, task_coor_2d
            INTEGER :: nproc

            CALL mp_environ(nproc, dims_2d, task_coor_2d, pgrid%mp_comm_2d)
            CALL mp_environ(nproc, dims_2d, task_coor_2d, pgrid%mp_comm_2d)
            CALL dbt_get_mapping_info(pgrid%nd_index_grid, dims_nd=dims)
            task_coor = get_nd_indices_pgrid(pgrid%nd_index_grid, task_coor_2d)
         END SUBROUTINE

! **************************************************************************************************
!> \brief As block_filter
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_filter(tensor, eps)
            TYPE(dbt_type), INTENT(INOUT)    :: tensor
            REAL(dp), INTENT(IN)                :: eps

            CALL dbt_tas_filter(tensor%matrix_rep, eps)

         END SUBROUTINE

! **************************************************************************************************
!> \brief As block_get_info but for tensors
!> \param nblks_total number of blocks along each dimension
!> \param nfull_total number of elements along each dimension
!> \param nblks_local local number of blocks along each dimension
!> \param nfull_local local number of elements along each dimension
!> \param my_ploc process coordinates in process grid
!> \param pdims process grid dimensions
!> \param blks_local_${idim}$ local blocks along dimension ${idim}$
!> \param proc_dist_${idim}$ distribution along dimension ${idim}$
!> \param blk_size_${idim}$ block sizes along dimension ${idim}$
!> \param blk_offset_${idim}$ block offsets along dimension ${idim}$
!> \param distribution distribution object
!> \param name name of tensor
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_get_info(tensor, nblks_total, &
                                 nfull_total, &
                                 nblks_local, &
                                 nfull_local, &
                                 pdims, &
                                 my_ploc, &
                                 ${varlist("blks_local")}$, &
                                 ${varlist("proc_dist")}$, &
                                 ${varlist("blk_size")}$, &
                                 ${varlist("blk_offset")}$, &
                                 distribution, &
                                 name)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(ndims_tensor(tensor)) :: nblks_total
            INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(ndims_tensor(tensor)) :: nfull_total
            INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(ndims_tensor(tensor)) :: nblks_local
            INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(ndims_tensor(tensor)) :: nfull_local
            INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(ndims_tensor(tensor)) :: my_ploc
            INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(ndims_tensor(tensor)) :: pdims
            #:for idim in range(1, maxdim+1)
               INTEGER, DIMENSION(dbt_nblks_local(tensor, ${idim}$)), INTENT(OUT), OPTIONAL :: blks_local_${idim}$
               INTEGER, DIMENSION(dbt_nblks_total(tensor, ${idim}$)), INTENT(OUT), OPTIONAL :: proc_dist_${idim}$
               INTEGER, DIMENSION(dbt_nblks_total(tensor, ${idim}$)), INTENT(OUT), OPTIONAL :: blk_size_${idim}$
               INTEGER, DIMENSION(dbt_nblks_total(tensor, ${idim}$)), INTENT(OUT), OPTIONAL :: blk_offset_${idim}$
            #:endfor
            TYPE(dbt_distribution_type), INTENT(OUT), OPTIONAL    :: distribution
            CHARACTER(len=*), INTENT(OUT), OPTIONAL                   :: name
            INTEGER, DIMENSION(ndims_tensor(tensor))                  :: pdims_tmp, my_ploc_tmp

            IF (PRESENT(nblks_total)) CALL dbt_get_mapping_info(tensor%nd_index_blk, dims_nd=nblks_total)
            IF (PRESENT(nfull_total)) CALL dbt_get_mapping_info(tensor%nd_index, dims_nd=nfull_total)
            IF (PRESENT(nblks_local)) nblks_local(:) = tensor%nblks_local
            IF (PRESENT(nfull_local)) nfull_local(:) = tensor%nfull_local

            IF (PRESENT(my_ploc) .OR. PRESENT(pdims)) CALL mp_environ_pgrid(tensor%pgrid, pdims_tmp, my_ploc_tmp)
            IF (PRESENT(my_ploc)) my_ploc = my_ploc_tmp
            IF (PRESENT(pdims)) pdims = pdims_tmp

            #:for idim in range(1, maxdim+1)
               IF (${idim}$ <= ndims_tensor(tensor)) THEN
                  IF (PRESENT(blks_local_${idim}$)) CALL get_ith_array(tensor%blks_local, ${idim}$, &
                                                                       dbt_nblks_local(tensor, ${idim}$), &
                                                                       blks_local_${idim}$)
                  IF (PRESENT(proc_dist_${idim}$)) CALL get_ith_array(tensor%nd_dist, ${idim}$, &
                                                                      dbt_nblks_total(tensor, ${idim}$), &
                                                                      proc_dist_${idim}$)
                  IF (PRESENT(blk_size_${idim}$)) CALL get_ith_array(tensor%blk_sizes, ${idim}$, &
                                                                     dbt_nblks_total(tensor, ${idim}$), &
                                                                     blk_size_${idim}$)
                  IF (PRESENT(blk_offset_${idim}$)) CALL get_ith_array(tensor%blk_offsets, ${idim}$, &
                                                                       dbt_nblks_total(tensor, ${idim}$), &
                                                                       blk_offset_${idim}$)
               END IF
            #:endfor

            IF (PRESENT(distribution)) distribution = dbt_distribution(tensor)
            IF (PRESENT(name)) name = tensor%name

         END SUBROUTINE

! **************************************************************************************************
!> \brief As block_get_num_blocks: get number of local blocks
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION dbt_get_num_blocks(tensor) RESULT(num_blocks)
            TYPE(dbt_type), INTENT(IN)    :: tensor
            INTEGER                           :: num_blocks
            num_blocks = dbt_tas_get_num_blocks(tensor%matrix_rep)
         END FUNCTION

! **************************************************************************************************
!> \brief Get total number of blocks
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION dbt_get_num_blocks_total(tensor) RESULT(num_blocks)
            TYPE(dbt_type), INTENT(IN)    :: tensor
            INTEGER(KIND=int_8)               :: num_blocks
            num_blocks = dbt_tas_get_num_blocks_total(tensor%matrix_rep)
         END FUNCTION

! **************************************************************************************************
!> \brief Clear tensor (s.t. it does not contain any blocks)
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_clear(tensor)
            TYPE(dbt_type), INTENT(INOUT) :: tensor

            CALL dbt_tas_clear(tensor%matrix_rep)
         END SUBROUTINE

! **************************************************************************************************
!> \brief Finalize tensor, as block_finalize. This should be taken care of internally in DBT
!>        tensors, there should not be any need to call this routine outside of DBT tensors.
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_finalize(tensor)
            TYPE(dbt_type), INTENT(INOUT) :: tensor
            CALL dbt_tas_finalize(tensor%matrix_rep)
         END SUBROUTINE

! **************************************************************************************************
!> \brief as block_scale
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_scale(tensor, alpha)
            TYPE(dbt_type), INTENT(INOUT) :: tensor
            REAL(dp), INTENT(IN) :: alpha
            CALL dbm_scale(tensor%matrix_rep%matrix, alpha)
         END SUBROUTINE

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION dbt_get_nze(tensor)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER                        :: dbt_get_nze
            dbt_get_nze = dbt_tas_get_nze(tensor%matrix_rep)
         END FUNCTION

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         FUNCTION dbt_get_nze_total(tensor)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER(KIND=int_8)            :: dbt_get_nze_total
            dbt_get_nze_total = dbt_tas_get_nze_total(tensor%matrix_rep)
         END FUNCTION

! **************************************************************************************************
!> \brief local number of blocks along dimension idim
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION dbt_nblks_local(tensor, idim)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER, INTENT(IN) :: idim
            INTEGER :: dbt_nblks_local

            IF (idim > ndims_tensor(tensor)) THEN
               dbt_nblks_local = 0
            ELSE
               dbt_nblks_local = tensor%nblks_local(idim)
            END IF

         END FUNCTION

! **************************************************************************************************
!> \brief total numbers of blocks along dimension idim
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION dbt_nblks_total(tensor, idim)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER, INTENT(IN) :: idim
            INTEGER :: dbt_nblks_total

            IF (idim > ndims_tensor(tensor)) THEN
               dbt_nblks_total = 0
            ELSE
               dbt_nblks_total = tensor%nd_index_blk%dims_nd(idim)
            END IF
         END FUNCTION

! **************************************************************************************************
!> \brief block size of block with index ind along dimension idim
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION dbt_blk_size(tensor, ind, idim)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER, DIMENSION(ndims_tensor(tensor)), &
               INTENT(IN) :: ind
            INTEGER, INTENT(IN) :: idim
            INTEGER, DIMENSION(ndims_tensor(tensor)) :: blk_size
            INTEGER :: dbt_blk_size

            IF (idim > ndims_tensor(tensor)) THEN
               dbt_blk_size = 0
            ELSE
               blk_size(:) = get_array_elements(tensor%blk_sizes, ind)
               dbt_blk_size = blk_size(idim)
            END IF
         END FUNCTION

! **************************************************************************************************
!> \brief how many tensor dimensions are mapped to matrix row
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION ndims_matrix_row(tensor)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER(int_8) :: ndims_matrix_row

            ndims_matrix_row = ndims_mapping_row(tensor%nd_index_blk)

         END FUNCTION

! **************************************************************************************************
!> \brief how many tensor dimensions are mapped to matrix column
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION ndims_matrix_column(tensor)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER(int_8) :: ndims_matrix_column

            ndims_matrix_column = ndims_mapping_column(tensor%nd_index_blk)
         END FUNCTION

! **************************************************************************************************
!> \brief returns an estimate of maximum number of local blocks in tensor
!>        (irrespective of the actual number of currently present blocks)
!>        this estimate is based on the following assumption: tensor data is dense and
!>        load balancing is within a factor of 2
!> \author Patrick Seewald
! **************************************************************************************************
         PURE FUNCTION dbt_max_nblks_local(tensor) RESULT(blk_count)
            TYPE(dbt_type), INTENT(IN) :: tensor
            INTEGER :: blk_count, nproc
            INTEGER, DIMENSION(ndims_tensor(tensor)) :: bdims
            INTEGER(int_8) :: blk_count_total
            INTEGER, PARAMETER :: max_load_imbalance = 2

            CALL dbt_get_mapping_info(tensor%nd_index_blk, dims_nd=bdims)

            blk_count_total = PRODUCT(INT(bdims, int_8))

            ! can not call an MPI routine due to PURE
            !CALL mp_environ(nproc, myproc, tensor%pgrid%mp_comm_2d)
            nproc = tensor%pgrid%nproc

            blk_count = INT(blk_count_total/nproc*max_load_imbalance)

         END FUNCTION

! **************************************************************************************************
!> \brief get a load-balanced and randomized distribution along one tensor dimension
!> \param nblk number of blocks (along one tensor dimension)
!> \param nproc number of processes (along one process grid dimension)
!> \param blk_size block sizes
!> \param dist distribution
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_default_distvec(nblk, nproc, blk_size, dist)
            INTEGER, INTENT(IN)                                :: nblk
            INTEGER, INTENT(IN)                                :: nproc
            INTEGER, DIMENSION(nblk), INTENT(IN)                :: blk_size
            INTEGER, DIMENSION(nblk), INTENT(OUT)               :: dist

            CALL dbt_tas_default_distvec(nblk, nproc, blk_size, dist)
         END SUBROUTINE

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
         SUBROUTINE dbt_copy_contraction_storage(tensor_in, tensor_out)
            TYPE(dbt_type), INTENT(IN) :: tensor_in
            TYPE(dbt_type), INTENT(INOUT) :: tensor_out
            TYPE(dbt_contraction_storage), ALLOCATABLE :: tensor_storage_tmp
            TYPE(dbt_tas_mm_storage), ALLOCATABLE :: tas_storage_tmp

            IF (tensor_in%matrix_rep%do_batched > 0) THEN
               ALLOCATE (tas_storage_tmp, SOURCE=tensor_in%matrix_rep%mm_storage)
               ! transfer data for batched contraction
               IF (ALLOCATED(tensor_out%matrix_rep%mm_storage)) DEALLOCATE (tensor_out%matrix_rep%mm_storage)
               CALL move_alloc(tas_storage_tmp, tensor_out%matrix_rep%mm_storage)
            END IF
            CALL dbt_tas_set_batched_state(tensor_out%matrix_rep, state=tensor_in%matrix_rep%do_batched, &
                                           opt_grid=tensor_in%matrix_rep%has_opt_pgrid)
            IF (ALLOCATED(tensor_in%contraction_storage)) THEN
               ALLOCATE (tensor_storage_tmp, SOURCE=tensor_in%contraction_storage)
            END IF
            IF (ALLOCATED(tensor_out%contraction_storage)) DEALLOCATE (tensor_out%contraction_storage)
            IF (ALLOCATED(tensor_storage_tmp)) CALL move_alloc(tensor_storage_tmp, tensor_out%contraction_storage)

         END SUBROUTINE

      END MODULE
