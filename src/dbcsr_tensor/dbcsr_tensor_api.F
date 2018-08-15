!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   This is the start of a dbcsr_tensor_api, all publically needed functions
!>          are exported here. The others remain private to the library.
!>          Currently, this is the CP2K used set.
!>          Ultimately, a reduced subset and well defined api will remain,
!>          possibly grouped in to standard and expert api.
!>          Currently, this is work in progress.
!> \author  Alfio Lazzaro
! **************************************************************************************************
MODULE dbcsr_tensor_api
   USE dbcsr_tensor,                    ONLY: dbcsr_t_contract,&
                                              dbcsr_t_get_block,&
                                              dbcsr_t_get_stored_coordinates,&
                                              dbcsr_t_put_block,&
                                              dbcsr_t_reserve_blocks,&
                                              dbcsr_t_split_blocks
   USE dbcsr_tensor_block,              ONLY: dbcsr_t_iterator_blocks_left,&
                                              dbcsr_t_iterator_next_block,&
                                              dbcsr_t_iterator_start,&
                                              dbcsr_t_iterator_stop,&
                                              dbcsr_t_iterator_type
   USE dbcsr_tensor_types,              ONLY: dbcsr_t_create,&
                                              dbcsr_t_destroy,&
                                              dbcsr_t_distribution_destroy,&
                                              dbcsr_t_distribution_new,&
                                              dbcsr_t_distribution_type,&
                                              dbcsr_t_nd_mp_comm,&
                                              dbcsr_t_nd_mp_free,&
                                              dbcsr_t_type

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_api'

   PUBLIC :: dbcsr_t_contract
   PUBLIC :: dbcsr_t_get_block
   PUBLIC :: dbcsr_t_get_stored_coordinates
   PUBLIC :: dbcsr_t_put_block
   PUBLIC :: dbcsr_t_reserve_blocks
   PUBLIC :: dbcsr_t_create
   PUBLIC :: dbcsr_t_destroy
   PUBLIC :: dbcsr_t_distribution_destroy
   PUBLIC :: dbcsr_t_distribution_new
   PUBLIC :: dbcsr_t_distribution_type
   PUBLIC :: dbcsr_t_nd_mp_comm
   PUBLIC :: dbcsr_t_nd_mp_free
   PUBLIC :: dbcsr_t_type
   PUBLIC :: dbcsr_t_iterator_next_block
   PUBLIC :: dbcsr_t_iterator_blocks_left
   PUBLIC :: dbcsr_t_iterator_stop
   PUBLIC :: dbcsr_t_iterator_start
   PUBLIC :: dbcsr_t_iterator_type
   PUBLIC :: dbcsr_t_split_blocks

END MODULE dbcsr_tensor_api
