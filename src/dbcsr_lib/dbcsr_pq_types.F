!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group   !
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
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type
  USE dbcsr_error_handling
  USE dbcsr_types,                     ONLY: dbcsr_data_obj

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_pq_types'

  ! Definitions for the members of a parameter stack.
  INTEGER, PARAMETER, PUBLIC :: dbcsr_ps_width = 7
  INTEGER, PARAMETER, PUBLIC :: p_m = 1
  INTEGER, PARAMETER, PUBLIC :: p_n = 2
  INTEGER, PARAMETER, PUBLIC :: p_k = 3
  INTEGER, PARAMETER, PUBLIC :: p_a_first = 4
  INTEGER, PARAMETER, PUBLIC :: p_b_first = 5
  INTEGER, PARAMETER, PUBLIC :: p_c_first = 6
  INTEGER, PARAMETER, PUBLIC :: p_c_blk = 7


  ! Possible queue states.
  INTEGER, PARAMETER :: dbcsr_ps_state_empty   = 0
  INTEGER, PARAMETER :: dbcsr_ps_state_filling = 1
  INTEGER, PARAMETER :: dbcsr_ps_state_filled  = 2
  INTEGER, PARAMETER :: dbcsr_ps_state_queued  = 3
  INTEGER, PARAMETER :: dbcsr_ps_state_working = 4

  ! The maximum number of regions to use for region-based locking of
  ! target data areas.
  INTEGER, PARAMETER :: max_regions = 1024
  

  LOGICAL, PARAMETER :: careful_mod = .FALSE.

  PUBLIC :: dbcsr_ps_type, dbcsr_ps_obj
  PUBLIC :: dbcsr_ps_group_type, dbcsr_ps_group_obj
  PUBLIC :: dbcsr_ps_set_type
  PUBLIC :: dbcsr_ps_target_obj, dbcsr_ps_target_type

  PUBLIC :: dbcsr_pq_type
  PUBLIC :: dbcsr_pq_type_p
  PUBLIC :: dbcsr_pq_all_type

  PUBLIC :: dbcsr_ps_state_empty, dbcsr_ps_state_filling,&
            dbcsr_ps_state_filled,&
            dbcsr_ps_state_queued, dbcsr_ps_state_working

  ! Just for correctness
  PUBLIC :: int_p_type

  PUBLIC :: max_regions

  PUBLIC :: ps_type_global
  PUBLIC :: ps_host_type
  PUBLIC :: ps_host_global
  PUBLIC :: ps_host_parameter


  !> \brief Product matrix meta information.
  !>
  !> Every product data area (i.e., each thread's) should have one of these
  !> objects describing its current status.
  !> \var product_data_area   Product matrix data
  !> \var product_data_cuda   Product matrix data on accelerator
  !> \var has_c_data          Whether product data data area is valid
  !> \var has_cuda_c_data     Whether CUDA product data descriptor is valid
  !> \var has_cuda_ab_data    Whether CUDA data descriptors (for left and right
  !>                          matrices) are valid
  !> \var zero_first    First element of C matrix data to zero
  !> \var zero_last     Last element of C matrix data to zero
  !> \var last_c_blk    Last block of the C matrix
  !> \var c_locks_dev   Locks on the accelerator card
  !> \var stack_state_dev   Value of "empty stack" to copy off accelerator
  !>                        card, which signals the stack can be reused.
  !> \var params_dev        Parameters on the card
  !> \var owner             Thread owning this target
  !> \var lock_owner        Current lock owner
  !> \var refcount          Reference counter
  !> \var target_lock       Lock for this instance
  !> \var owner_wants_lock  Other threads should let owner take lock
  !> \var n_use             Number of threads holding region locks on this
  !>                        instance
  !> \var bit_shift         Bit shift used to calculate regions from
  !>                        data area offset
  !> \var region_locks      Region locks
  TYPE dbcsr_ps_target_type
     TYPE(dbcsr_data_obj)               :: product_data_area
     TYPE(dbcsr_cuda_mem_type), POINTER :: product_data_cuda
     LOGICAL :: has_c_data, has_cuda_c_data
     INTEGER :: zero_first, zero_last
     INTEGER :: last_c_blk
     TYPE(dbcsr_cuda_mem_type), POINTER :: c_locks_dev, stack_state_dev
     TYPE(dbcsr_cuda_mem_type), POINTER :: params_dev
     !$ INTEGER(KIND=omp_lock_kind) :: target_lock = 0
     !$ LOGICAL :: owner_wants_lock
     !$ INTEGER :: n_use
     !$ INTEGER :: bit_shift
     !$ INTEGER(KIND=omp_lock_kind), DIMENSION(:), POINTER :: region_locks
     INTEGER :: owner, lock_owner
     INTEGER :: refcount
  END TYPE dbcsr_ps_target_type


  !> \brief Object for dbcsr_ps_target_type
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
  !> \var defined_mnk         All stack entries have a common m, n, and k values
  !> \var size_bin            Which bin to put this stack into
  !> \var state               State of the stack (used for asynchrony).
  !> \var own_state           Whether the state is allocated and must be
  !>                          deallocated or whether it points into a
  !>                          common state array. \see
  !>                          dbcsr_ps_set_type
  TYPE dbcsr_ps_type
     INTEGER :: id, sid
     INTEGER, DIMENSION(:), POINTER :: parameters
     INTEGER :: size_hint
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
     INTEGER :: size_bin
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
     INTEGER :: nbuffers, nmemregs, group_size
     INTEGER :: next_buffer, next_memreg
  END TYPE dbcsr_ps_set_type

  ! For host parameter stack allocation on CUDA versions
   TYPE ps_type_global
        TYPE(int_p_type),DIMENSION(:), POINTER   :: ps_local
   END TYPE

   TYPE(ps_type_global), DIMENSION(:), POINTER  :: ps_host_global

   TYPE ps_host_type
       TYPE(int_p_type),DIMENSION(:), POINTER   :: mem_stack
   END TYPE ps_host_type

   TYPE(ps_host_type), DIMENSION(:), POINTER  :: ps_host_parameter



  !> \brief Queue of parameter stacks
  !>
  !> Not a real queue but a list.
  !> \var stacks       The parameter stacks
  !> \var nstacks_total  The number of valid stacks
  !> \var nstacks_l    Number of valid stack for each size bin
  !> \var done         Signals that the queue owner is done with indexing
  !> \var all_queues   Collection of queues from all threads
  !> \var flush_level  The approximate level at which the queue is getting full
  !> \var queue_lock   Lock used to access queue
  TYPE dbcsr_pq_type
     TYPE(dbcsr_ps_obj), DIMENSION(:,:), ALLOCATABLE :: stacks
     INTEGER :: nstacks_total, nbins
     INTEGER, DIMENSION(:), ALLOCATABLE :: nstacks
     !$ INTEGER(KIND=omp_lock_kind), POINTER :: queue_lock
     LOGICAL :: done
     TYPE(dbcsr_pq_all_type), POINTER :: all_queues
     INTEGER :: flush_level
     INTEGER :: handoff
     INTEGER :: nworking
  END TYPE dbcsr_pq_type

  TYPE dbcsr_pq_type_p
     TYPE(dbcsr_pq_type), POINTER :: pq
  END TYPE dbcsr_pq_type_p

  !> \brief Structure to hold all queues
  !> \var queues  Collection of all queues
  !> \var n_done  Number of threads that have finished indexing
  !> \var n_working  Number of threads that are still indexing
  !> \var working    Bitmap of threads showing working status
  TYPE dbcsr_pq_all_type
     TYPE(dbcsr_pq_type_p), DIMENSION(:), POINTER :: queues
     INTEGER :: n_done, n_working
     LOGICAL, DIMENSION(:), POINTER :: working
  END TYPE dbcsr_pq_all_type

END MODULE dbcsr_pq_types
