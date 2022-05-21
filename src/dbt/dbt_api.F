!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief This is the start of a dbt_api, all publically needed functions
!>        are exported here. The others remain private to the library.
!>        Currently, this is the CP2K used set.
!>        Ultimately, a reduced subset and well defined api will remain,
!>        possibly grouped in to standard and expert api.
!>        Currently, this is work in progress.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_api

   USE dbt_block,                       ONLY: dbt_iterator_blocks_left,&
                                              dbt_iterator_next_block,&
                                              dbt_iterator_num_blocks,&
                                              dbt_iterator_start,&
                                              dbt_iterator_stop,&
                                              dbt_iterator_type
   USE dbt_index,                       ONLY: dbt_get_mapping_info
   USE dbt_io,                          ONLY: dbt_write_blocks,&
                                              dbt_write_split_info,&
                                              dbt_write_tensor_dist,&
                                              dbt_write_tensor_info
   USE dbt_methods,                     ONLY: &
        dbt_batched_contract_finalize, dbt_batched_contract_init, dbt_contract, dbt_copy, &
        dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, dbt_get_block, &
        dbt_get_stored_coordinates, dbt_put_block, dbt_reserve_blocks
   USE dbt_split,                       ONLY: dbt_split_blocks
   USE dbt_test,                        ONLY: dbt_checksum,&
                                              dbt_contract_test
   USE dbt_types,                       ONLY: &
        dbt_blk_size, dbt_blk_sizes, dbt_clear, dbt_create, dbt_default_distvec, dbt_destroy, &
        dbt_dims => dims_tensor, dbt_distribution_destroy, dbt_distribution_new, &
        dbt_distribution_type, dbt_filter, dbt_finalize, dbt_get_info, dbt_get_num_blocks, &
        dbt_get_num_blocks_total, dbt_get_nze, dbt_get_nze_total, dbt_max_nblks_local, &
        dbt_mp_dims_create, dbt_mp_environ_pgrid => mp_environ_pgrid, dbt_nblks_local, &
        dbt_nblks_total, dbt_nd_mp_comm, dbt_nd_mp_free, dbt_ndims => ndims_tensor, &
        dbt_ndims_matrix_column => ndims_matrix_column, dbt_ndims_matrix_row => ndims_matrix_row, &
        dbt_pgrid_change_dims, dbt_pgrid_create, dbt_pgrid_create_expert, dbt_pgrid_destroy, &
        dbt_pgrid_type, dbt_scale, dbt_type

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_api'

   PUBLIC :: dbt_contract, dbt_contract_test
   PUBLIC :: dbt_get_block
   PUBLIC :: dbt_get_stored_coordinates
   PUBLIC :: dbt_put_block
   PUBLIC :: dbt_reserve_blocks
   PUBLIC :: dbt_create
   PUBLIC :: dbt_destroy
   PUBLIC :: dbt_distribution_destroy
   PUBLIC :: dbt_distribution_new
   PUBLIC :: dbt_distribution_type
   PUBLIC :: dbt_nd_mp_comm
   PUBLIC :: dbt_nd_mp_free
   PUBLIC :: dbt_type
   PUBLIC :: dbt_iterator_next_block
   PUBLIC :: dbt_iterator_blocks_left
   PUBLIC :: dbt_iterator_stop
   PUBLIC :: dbt_iterator_start
   PUBLIC :: dbt_iterator_type
   PUBLIC :: dbt_iterator_num_blocks
   PUBLIC :: dbt_split_blocks
   PUBLIC :: dbt_pgrid_type
   PUBLIC :: dbt_pgrid_create
   PUBLIC :: dbt_pgrid_create_expert
   PUBLIC :: dbt_pgrid_destroy
   PUBLIC :: dbt_filter
   PUBLIC :: dbt_mp_environ_pgrid
   PUBLIC :: dbt_copy_matrix_to_tensor
   PUBLIC :: dbt_blk_sizes
   PUBLIC :: dbt_copy
   PUBLIC :: dbt_copy_tensor_to_matrix
   PUBLIC :: dbt_get_info
   PUBLIC :: dbt_checksum
   PUBLIC :: dbt_finalize
   PUBLIC :: dbt_scale
   PUBLIC :: dbt_get_num_blocks, dbt_get_num_blocks_total
   PUBLIC :: dbt_get_nze, dbt_get_nze_total
   PUBLIC :: dbt_clear
   PUBLIC :: dbt_get_mapping_info
   PUBLIC :: dbt_write_split_info
   PUBLIC :: dbt_write_blocks
   PUBLIC :: dbt_write_tensor_dist
   PUBLIC :: dbt_write_tensor_info
   PUBLIC :: dbt_mp_dims_create
   PUBLIC :: dbt_batched_contract_init
   PUBLIC :: dbt_batched_contract_finalize
   PUBLIC :: dbt_ndims
   PUBLIC :: dbt_dims
   PUBLIC :: dbt_pgrid_change_dims
   PUBLIC :: dbt_ndims_matrix_row
   PUBLIC :: dbt_ndims_matrix_column
   PUBLIC :: dbt_nblks_local
   PUBLIC :: dbt_nblks_total
   PUBLIC :: dbt_blk_size
   PUBLIC :: dbt_max_nblks_local
   PUBLIC :: dbt_default_distvec

END MODULE dbt_api
