!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Parameter stacks and their queues
!> \author  Urban Borstnik
!> \date    2011-06-17
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-06-17
! *****************************************************************************
MODULE dbcsr_pq_types
  USE array_types,                     ONLY: array_i1d_obj
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type
  USE dbcsr_error_handling
  USE dbcsr_types,                     ONLY: dbcsr_data_obj

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_pq_types'

  INTEGER, PARAMETER :: dbcsr_ps_width = 7
  INTEGER, PARAMETER :: dbcsr_ps_fqx_width = 3

  !> \var dbcsr_ps_fqx_tmp_offset  Index in the xlate array describing 
  !>                               the location of a block in src_data.
  !>                               \see dbcsr_ps_fin_type
  !> \var dbcsr_ps_fqx_tgt_offset  Index in the xlate array describing 
  !>                               the location of a block in the target data
  !>                               area.
  !>                               \see dbcsr_ps_fin_type
  !> \var dbcsr_ps_fqx_len         Index in the xlate array describing
  !>                               the block length.
  !>                               \see dbcsr_ps_fin_type
  INTEGER, PARAMETER :: dbcsr_ps_fqx_tmp_offset = 1
  INTEGER, PARAMETER :: dbcsr_ps_fqx_tgt_offset = 2
  INTEGER, PARAMETER :: dbcsr_ps_fqx_len = 3

  INTEGER, PARAMETER :: dbcsr_ps_state_empty   = 0
  INTEGER, PARAMETER :: dbcsr_ps_state_filling = 1
  INTEGER, PARAMETER :: dbcsr_ps_state_filled  = 2
  INTEGER, PARAMETER :: dbcsr_ps_state_working = 3

  LOGICAL, PARAMETER :: careful_mod = .TRUE.

  PUBLIC :: dbcsr_ps_type, dbcsr_ps_obj
  PUBLIC :: dbcsr_ps_group_type, dbcsr_ps_group_obj
  PUBLIC :: dbcsr_ps_set_type
  PUBLIC :: dbcsr_pq_type
  PUBLIC :: dbcsr_ps_target_obj, dbcsr_ps_target_type

  PUBLIC :: dbcsr_ps_width
  PUBLIC :: dbcsr_ps_state_empty, dbcsr_ps_state_filling,&
            dbcsr_ps_state_filled, dbcsr_ps_state_working

  PUBLIC :: dbcsr_ps_fin_q_type, dbcsr_ps_fin_type

  PUBLIC :: dbcsr_ps_fqx_width
  PUBLIC :: dbcsr_ps_fqx_tmp_offset, dbcsr_ps_fqx_tgt_offset, dbcsr_ps_fqx_len

  ! Just for correctness
  PUBLIC :: int_p_type


!> \brief Used for finalizing data crunched by foreign cores.
!> \var src_data  Temporary data to be moved to the target area
!> \var xlate     Translation array.
!>                \see dbcsr_ps_fqx_tmp_offset
!>                \see dbcsr_ps_fqx_tgt_offset
!>                \see dbcsr_ps_fqx_len
  TYPE dbcsr_ps_fin_type
     TYPE(dbcsr_data_obj)  :: src_data
     TYPE(array_i1d_obj)   :: xlate
     INTEGER :: n
  END TYPE dbcsr_ps_fin_type

  TYPE dbcsr_ps_fin_q_type
     TYPE(dbcsr_ps_fin_type), DIMENSION(:), POINTER :: f
     INTEGER :: n
  END TYPE dbcsr_ps_fin_q_type


  !> \brief Product matrix meta information.
  !>
  !> Every product data area (i.e., each thread's) should have one of these
  !> objects describing its current status.
  !> \var has_c_data          Whether product data data area is valid
  !> \var has_cuda_ab_data    Whether CUDA data descriptors (for left and right
  !>                          matrices) are valid
  !> \var has_c_data          Whether CUDA product data descriptor is valid
  !> \var product_data_area   Product matrix data
  !> \var zero_first    First element of C matrix data to zero
  !> \var zero_last     Last element of C matrix data to zero
  !> \var last_c_blk    Last block of the C matrix
  !> \var fin_queue     Foreign work finalization queue
  TYPE dbcsr_ps_target_type
     TYPE(dbcsr_data_obj)               :: product_data_area
     TYPE(dbcsr_cuda_mem_type), POINTER :: product_data_cuda
     LOGICAL :: has_c_data, has_cuda_c_data
     INTEGER :: zero_first, zero_last
     INTEGER :: last_c_blk
     TYPE(dbcsr_ps_fin_q_type) :: fin_queue
     TYPE(dbcsr_cuda_mem_type), POINTER :: c_locks_dev, stack_state_dev
     TYPE(dbcsr_cuda_mem_type), POINTER :: params_dev
     INTEGER :: refcount
  END TYPE dbcsr_ps_target_type


  TYPE dbcsr_ps_target_obj
     TYPE(dbcsr_ps_target_type), POINTER :: t
  END TYPE dbcsr_ps_target_obj


  !> \brief Parameter stack
  !> \var parameters    The data stored on the stack
  !> \var own_data      Whether the parameters array is owned by this object
  !> \var start_p       Offset into master array
  !> \var mem_type      Memory type used for the parameters variable.
  !> \var stack_p       Last element of the stack
  !> \var driver        Who/what does calculations
  !> \var driver_desc   Driver-specific instruction
  !> \var has_ab_data         Whether data areas (for left and right matrices)
  !>                          are valid
  !> \var left_data_area      Left matrix data
  !> \var right_data_area     Right matrix data
  !> \var product_data_cuda   Product matrix data in CUDA address space
  !> \var left_data_cuda      Left matrix data in CUDA address space
  !> \var right_data_cuda     Right matrix data in CUDA address space
  !> \var state               State of the stack
  !> \var refcount            Reference counter
  !> \var m                   Common m shared by all stack members; otherwise 0
  !> \var n                   Common n shared by all stack members; otherwise 0
  !> \var k                   Common k shared by all stack members; otherwise 0
  !> \var max_m               Maximum m that can be found in this stack
  !> \var max_n               Maximum n that can be found in this stack
  !> \var max_k               Maximum n that can be found in this stack
  !> \var state               State of the stack (used for asynchrony).
  !> \var own_state           Whether the state is allocated and must be
  !>                          deallocated or whether it points into a
  !>                          common state array. \see
  !>                          dbcsr_ps_set_type
  TYPE dbcsr_ps_type
     INTEGER, DIMENSION(:), POINTER :: parameters
     LOGICAL :: own_data
     INTEGER :: start_p
     INTEGER :: mem_type
     INTEGER :: stack_p
     INTEGER :: driver, driver_desc
     TYPE(dbcsr_ps_target_obj) :: t
     LOGICAL :: has_ab_data, has_cuda_ab_data, has_target
     TYPE(dbcsr_data_obj) :: left_data_area, right_data_area
     TYPE(dbcsr_cuda_mem_type), POINTER :: left_data_cuda, right_data_cuda
     INTEGER, POINTER :: state
     LOGICAL :: own_state
     LOGICAL :: defined_mnk
     INTEGER :: m, n, k, max_m, max_n, max_k
     INTEGER :: refcount
  END TYPE dbcsr_ps_type

  TYPE dbcsr_ps_obj
     TYPE(dbcsr_ps_type), POINTER :: s
  END TYPE dbcsr_ps_obj

  !> \brief A group of stacks
  !>
  !> All stacks in the group should share the same memory area that is
  !> described by the master dbcsr_ps_obj.
  TYPE dbcsr_ps_group_type
     TYPE(dbcsr_ps_obj) :: master
     TYPE(dbcsr_ps_obj), DIMENSION(:), POINTER :: stacks
  END TYPE dbcsr_ps_group_type


  TYPE dbcsr_ps_group_obj
     TYPE(dbcsr_ps_group_type) :: g
  END TYPE dbcsr_ps_group_obj

  TYPE int_p_type
     INTEGER, DIMENSION(:), POINTER :: DATA
     INTEGER :: mem_type
  END TYPE int_p_type

  !> \brief A set of parameter stacks
  !> \var groups     Groups, one for each pair of buffers and memory regions.
  !> \par Meaninings of stacks index positions
  !> * Index 1: Buffering (i.e., fill/read separation)
  !> * Index 2: Grouping by memory region (i.e., independent-memory
  !>            accelerator boards)
  !> \var all_states   States for all stacks in this set.  The states
  !>                   are gathered into one array that can be
  !>                   allocated with pinned memory.  Individual
  !>                   stacks point into here.
  TYPE dbcsr_ps_set_type
     TYPE(dbcsr_ps_group_type), DIMENSION(:,:), POINTER :: groups
     TYPE(int_p_type), DIMENSION(:), POINTER            :: all_states
     INTEGER :: nbuffers, nmemregs
     INTEGER :: next_buffer, next_memreg
  END TYPE dbcsr_ps_set_type

  !> \brief Parameter queue
  !> \var stacks    The parameter stacks
  !> \var nstacs    The number of valid stacks
  TYPE dbcsr_pq_type
     TYPE(dbcsr_ps_obj), DIMENSION(:), POINTER :: stacks
     INTEGER :: nstacks
  END TYPE dbcsr_pq_type

END MODULE dbcsr_pq_types
