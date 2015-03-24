!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the message passing library MPI
!> \par History
!>      JGH (02-Jan-2001): New error handling
!>                         Performance tools
!>      JGH (14-Jan-2001): New routines mp_comm_compare, mp_cart_coords,
!>                                      mp_rank_compare, mp_alltoall
!>      JGH (06-Feb-2001): New routines mp_comm_free
!>      JGH (22-Mar-2001): New routines mp_comm_dup
!>      fawzi (04-NOV-2004): storable performance info (for f77 interface)
!>      Wrapper routine for mpi_gatherv added (22.12.2005,MK)
!>      JGH (13-Feb-2006): Flexibel precision
!>      JGH (15-Feb-2006): single precision mp_alltoall
!> \author JGH
! *****************************************************************************
MODULE message_passing
  USE kinds,                           ONLY: &
       default_string_length, dp, int_4, int_4_size, int_8, int_8_size, &
       real_4, real_4_size, real_8, real_8_size
  USE machine,                         ONLY: default_output_unit,&
                                             m_abort,&
                                             m_flush,&
                                             m_walltime
  USE ISO_C_BINDING,                   ONLY: C_PTR,C_LOC,C_F_POINTER
#if defined(__parallel) && ! defined(__HAS_NO_MPI_MOD)
  USE mpi  ! errors mean mpi installation and fortran compiler mismatch: see INSTALL (-D__HAS_NO_MPI_MOD)
#endif

  IMPLICIT NONE
  PRIVATE

  ! parameters that might be needed
#if defined(__parallel)
#if defined(__HAS_NO_MPI_MOD)
    INCLUDE "mpif.h"
#endif
  INTEGER, PARAMETER     :: MP_STD_REAL = MPI_DOUBLE_PRECISION
  INTEGER, PARAMETER     :: MP_STD_COMPLEX = MPI_DOUBLE_COMPLEX
  INTEGER, PARAMETER     :: MP_STD_HALF_REAL = MPI_REAL
  INTEGER, PARAMETER     :: MP_STD_HALF_COMPLEX = MPI_COMPLEX
#endif

#ifdef __parallel
  LOGICAL, PARAMETER :: cp2k_is_parallel=.TRUE.
  INTEGER, PARAMETER, PUBLIC :: mp_any_tag=MPI_ANY_TAG
  INTEGER, PARAMETER, PUBLIC :: mp_any_source=MPI_ANY_SOURCE
  INTEGER, PARAMETER, PUBLIC :: mp_comm_null=MPI_COMM_NULL
  INTEGER, PARAMETER, PUBLIC :: mp_request_null=MPI_REQUEST_NULL
  INTEGER, PARAMETER, PUBLIC :: mp_status_size=MPI_STATUS_SIZE
  INTEGER, PARAMETER, PUBLIC :: mp_address_kind= MPI_ADDRESS_KIND

  INTEGER, PARAMETER, PUBLIC :: file_offset=MPI_OFFSET_KIND
  INTEGER, PARAMETER, PUBLIC :: file_amode_create=MPI_MODE_CREATE
  INTEGER, PARAMETER, PUBLIC :: file_amode_rdonly=MPI_MODE_RDONLY
  INTEGER, PARAMETER, PUBLIC :: file_amode_wronly=MPI_MODE_WRONLY
  INTEGER, PARAMETER, PUBLIC :: file_amode_rdwr=MPI_MODE_RDWR
  INTEGER, PARAMETER, PUBLIC :: file_amode_excl=MPI_MODE_EXCL
#else
  LOGICAL, PARAMETER :: cp2k_is_parallel=.FALSE.
  INTEGER, PARAMETER, PUBLIC :: mp_any_tag=-1
  INTEGER, PARAMETER, PUBLIC :: mp_any_source=-2
  INTEGER, PARAMETER, PUBLIC :: mp_comm_null=-3
  INTEGER, PARAMETER, PUBLIC :: mp_request_null=-4
  INTEGER, PARAMETER, PUBLIC :: mp_status_size=-5

  INTEGER, PARAMETER, PUBLIC :: file_offset=int_8
  INTEGER, PARAMETER, PUBLIC :: file_amode_create=1
  INTEGER, PARAMETER, PUBLIC :: file_amode_rdonly=2
  INTEGER, PARAMETER, PUBLIC :: file_amode_wronly=4
  INTEGER, PARAMETER, PUBLIC :: file_amode_rdwr=8
  INTEGER, PARAMETER, PUBLIC :: file_amode_excl=64
  INTEGER, PARAMETER, PUBLIC :: mp_address_kind=int_8
#endif

  ! we need to fix this to a given number (crossing fingers)
  ! so that the serial code using Fortran stream IO and the MPI have the same sizes.
  INTEGER, PARAMETER, PUBLIC :: mpi_character_size=1
  INTEGER, PARAMETER, PUBLIC :: mpi_integer_size=4

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'message_passing'

#if defined(__parallel)
  ! internal reference counter used to debug communicator leaks
  INTEGER, PRIVATE, SAVE :: debug_comm_count
#endif

  ! init and error
  PUBLIC :: mp_world_init, mp_world_finalize
  PUBLIC :: mp_abort

  ! performance gathering
  PUBLIC :: mp_perf_env_type
  PUBLIC :: mp_perf_env_retain, mp_perf_env_release
  PUBLIC :: add_mp_perf_env, rm_mp_perf_env, get_mp_perf_env, describe_mp_perf_env

  ! informational / generation of sub comms
  PUBLIC :: mp_environ, mp_comm_compare, mp_cart_coords, mp_rank_compare
  PUBLIC :: mp_cart_create, mp_dims_create, mp_cart_rank, mp_cart_sub, mp_comm_free
  PUBLIC :: mp_comm_dup, mp_comm_split, mp_comm_split_direct
  PUBLIC :: cp2k_is_parallel
  PUBLIC :: mp_proc_name
  PUBLIC :: mp_probe

  ! message passing
  PUBLIC :: mp_bcast, mp_sum, mp_max, mp_maxloc, mp_minloc, mp_min, mp_sync
  PUBLIC :: mp_gather, mp_scatter, mp_alltoall, mp_sendrecv, mp_allgather
  PUBLIC :: mp_isend, mp_irecv
  PUBLIC :: mp_shift, mp_isendrecv, mp_wait, mp_waitall, mp_waitany, mp_testany
  PUBLIC :: mp_gatherv
  PUBLIC :: mp_send, mp_recv

  ! Memory management
  PUBLIC :: mp_allocate, mp_deallocate

  ! MPI re-ordering
  PUBLIC :: mp_reordering

  ! default communicators
  PUBLIC :: MPI_COMM_SELF, MPI_COMM_WORLD, MPI_COMM_NULL

  ! I/O
  PUBLIC :: mp_file_open, mp_file_close
  PUBLIC :: mp_file_write_at
  PUBLIC :: mp_file_write_at_all, mp_file_read_at_all
  PUBLIC :: mp_file_get_size

  ! some 'advanced types' currently only used for dbcsr
  PUBLIC :: mp_type_descriptor_type
  PUBLIC :: mp_type_make
  PUBLIC :: mp_type_size

  ! some benchmarking code
  PUBLIC :: mpi_perf_test

  ! one-sided communication
  PUBLIC :: mp_win_create, mp_win_free, mp_win_lock_all,&
            mp_win_unlock_all, mp_rget

#ifndef __parallel
  INTEGER, PARAMETER :: MPI_COMM_SELF=0,MPI_COMM_WORLD=0, MPI_COMM_NULL=-1
#endif

  ! assumed to be private

! Interface declarations for non-data-oriented subroutines.

  INTERFACE mp_environ
     MODULE PROCEDURE mp_environ_l, mp_environ_c, mp_environ_c2
  END INTERFACE

  INTERFACE mp_waitall
     MODULE PROCEDURE mp_waitall_1, mp_waitall_2
  END INTERFACE

  !
  ! interfaces to deal easily with scalars / vectors / matrice / ...
  ! of the different types (integers, doubles, logicals, characters)
  !
  INTERFACE mp_minloc
     MODULE PROCEDURE mp_minloc_iv,&
                      mp_minloc_lv,&
                      mp_minloc_rv,&
                      mp_minloc_dv
  END INTERFACE

  INTERFACE mp_maxloc
     MODULE PROCEDURE mp_maxloc_iv,&
                      mp_maxloc_lv,&
                      mp_maxloc_rv,&
                      mp_maxloc_dv
  END INTERFACE

  INTERFACE mp_shift
     MODULE PROCEDURE mp_shift_im, mp_shift_i,&
                      mp_shift_lm, mp_shift_l,&
                      mp_shift_rm, mp_shift_r,&
                      mp_shift_dm, mp_shift_d,&
                      mp_shift_cm, mp_shift_c,&
                      mp_shift_zm, mp_shift_z
  END INTERFACE

  INTERFACE mp_bcast
     MODULE PROCEDURE mp_bcast_i, mp_bcast_iv, mp_bcast_im, mp_bcast_i3,&
                      mp_bcast_l, mp_bcast_lv, mp_bcast_lm, mp_bcast_l3,&
                      mp_bcast_r, mp_bcast_rv, mp_bcast_rm, mp_bcast_r3,&
                      mp_bcast_d, mp_bcast_dv, mp_bcast_dm, mp_bcast_d3,&
                      mp_bcast_c, mp_bcast_cv, mp_bcast_cm, mp_bcast_c3,&
                      mp_bcast_z, mp_bcast_zv, mp_bcast_zm, mp_bcast_z3
     MODULE PROCEDURE mp_bcast_b, mp_bcast_bv
     MODULE PROCEDURE mp_bcast_av, mp_bcast_am
  END INTERFACE

  INTERFACE mp_sum
     MODULE PROCEDURE mp_sum_i, mp_sum_iv, mp_sum_im, mp_sum_im3,mp_sum_im4,mp_sum_im5,mp_sum_im6,&
                      mp_sum_l, mp_sum_lv, mp_sum_lm, mp_sum_lm3,mp_sum_lm4,mp_sum_lm5,mp_sum_lm6,&
                      mp_sum_r, mp_sum_rv, mp_sum_rm, mp_sum_rm3,mp_sum_rm4,mp_sum_rm5,mp_sum_rm6,&
                      mp_sum_d, mp_sum_dv, mp_sum_dm, mp_sum_dm3,mp_sum_dm4,mp_sum_dm5,mp_sum_dm6,&
                      mp_sum_c, mp_sum_cv, mp_sum_cm, mp_sum_cm3,mp_sum_cm4,mp_sum_cm5,mp_sum_cm6,&
                      mp_sum_z, mp_sum_zv, mp_sum_zm, mp_sum_zm3,mp_sum_zm4,mp_sum_zm5,mp_sum_zm6,&
                      mp_sum_root_iv, mp_sum_root_im,&
                      mp_sum_root_lv, mp_sum_root_lm,&
                      mp_sum_root_rv, mp_sum_root_rm,&
                      mp_sum_root_dv, mp_sum_root_dm,&
                      mp_sum_root_cv, mp_sum_root_cm,&
                      mp_sum_root_zv, mp_sum_root_zm
     MODULE PROCEDURE mp_sum_b
  END INTERFACE

  INTERFACE mp_max
     MODULE PROCEDURE mp_max_i, mp_max_iv,&
                      mp_max_l, mp_max_lv,&
                      mp_max_r, mp_max_rv,&
                      mp_max_d, mp_max_dv,&
                      mp_max_c, mp_max_cv,&
                      mp_max_z, mp_max_zv
  END INTERFACE

  INTERFACE mp_min
     MODULE PROCEDURE mp_min_i, mp_min_iv,&
                      mp_min_l, mp_min_lv,&
                      mp_min_r, mp_min_rv,&
                      mp_min_d, mp_min_dv,&
                      mp_min_c, mp_min_cv,&
                      mp_min_z, mp_min_zv
  END INTERFACE


  INTERFACE mp_gather
     MODULE PROCEDURE mp_gather_i, mp_gather_iv, mp_gather_im,&
                      mp_gather_l, mp_gather_lv, mp_gather_lm,&
                      mp_gather_r, mp_gather_rv, mp_gather_rm,&
                      mp_gather_d, mp_gather_dv, mp_gather_dm,&
                      mp_gather_c, mp_gather_cv, mp_gather_cm,&
                      mp_gather_z, mp_gather_zv, mp_gather_zm
  END INTERFACE

  INTERFACE mp_gatherv
    MODULE PROCEDURE mp_gatherv_iv,&
                     mp_gatherv_lv,&
                     mp_gatherv_rv,&
                     mp_gatherv_dv,&
                     mp_gatherv_cv,&
                     mp_gatherv_zv
  END INTERFACE

!> todo: move allgatherv to a separate declaration
  INTERFACE mp_allgather
     MODULE PROCEDURE &
          mp_allgather_i,&
          mp_allgather_i12, mp_allgather_i23, mp_allgather_i34,&
          mp_allgather_l,&
          mp_allgather_l12, mp_allgather_l23, mp_allgather_l34,&
          mp_allgather_r,&
          mp_allgather_r12, mp_allgather_r23, mp_allgather_r34,&
          mp_allgather_d,&
          mp_allgather_d12, mp_allgather_d23, mp_allgather_d34,&
          mp_allgather_c,&
          mp_allgather_c12, mp_allgather_c23, mp_allgather_c34,&
          mp_allgather_z,&
          mp_allgather_z12, mp_allgather_z23, mp_allgather_z34,&
          mp_allgatherv_iv,&
          mp_allgatherv_lv,&
          mp_allgatherv_rv,&
          mp_allgatherv_dv,&
          mp_allgatherv_cv,&
          mp_allgatherv_zv
  END INTERFACE

  INTERFACE mp_scatter
     MODULE PROCEDURE mp_scatter_iv,&
                      mp_scatter_lv,&
                      mp_scatter_rv,&
                      mp_scatter_dv,&
                      mp_scatter_cv,&
                      mp_scatter_zv
  END INTERFACE

  INTERFACE mp_sum_scatter
     MODULE PROCEDURE mp_sum_scatter_iv,&
                      mp_sum_scatter_lv,&
                      mp_sum_scatter_rv,&
                      mp_sum_scatter_dv,&
                      mp_sum_scatter_cv,&
                      mp_sum_scatter_zv
  END INTERFACE

  INTERFACE mp_alltoall
     MODULE PROCEDURE mp_alltoall_i, mp_alltoall_i22, mp_alltoall_i33,&
                      mp_alltoall_i44, mp_alltoall_i45, mp_alltoall_i34,&
                      mp_alltoall_i11v, mp_alltoall_i22v, mp_alltoall_i54,&
                      mp_alltoall_l, mp_alltoall_l22, mp_alltoall_l33,&
                      mp_alltoall_l44, mp_alltoall_l45, mp_alltoall_l34,&
                      mp_alltoall_l11v, mp_alltoall_l22v, mp_alltoall_l54,&
                      mp_alltoall_r, mp_alltoall_r22, mp_alltoall_r33,&
                      mp_alltoall_r44, mp_alltoall_r45, mp_alltoall_r34,&
                      mp_alltoall_r11v, mp_alltoall_r22v, mp_alltoall_r54,&
                      mp_alltoall_d, mp_alltoall_d22, mp_alltoall_d33,&
                      mp_alltoall_d44, mp_alltoall_d45, mp_alltoall_d34,&
                      mp_alltoall_d11v, mp_alltoall_d22v, mp_alltoall_d54,&
                      mp_alltoall_c, mp_alltoall_c22, mp_alltoall_c33,&
                      mp_alltoall_c44, mp_alltoall_c45, mp_alltoall_c34,&
                      mp_alltoall_c11v, mp_alltoall_c22v, mp_alltoall_c54,&
                      mp_alltoall_z, mp_alltoall_z22, mp_alltoall_z33,&
                      mp_alltoall_z44, mp_alltoall_z45, mp_alltoall_z34,&
                      mp_alltoall_z11v, mp_alltoall_z22v, mp_alltoall_z54
  END INTERFACE

  INTERFACE mp_send
     MODULE PROCEDURE mp_send_i,mp_send_iv,&
                      mp_send_l,mp_send_lv,&
                      mp_send_r,mp_send_rv,&
                      mp_send_d,mp_send_dv,&
                      mp_send_c,mp_send_cv,&
                      mp_send_z,mp_send_zv
  END INTERFACE

  INTERFACE mp_recv
     MODULE PROCEDURE mp_recv_i, mp_recv_iv,&
                      mp_recv_l, mp_recv_lv,&
                      mp_recv_r, mp_recv_rv,&
                      mp_recv_d, mp_recv_dv,&
                      mp_recv_c, mp_recv_cv,&
                      mp_recv_z, mp_recv_zv
  END INTERFACE

  INTERFACE mp_sendrecv
     MODULE PROCEDURE mp_sendrecv_iv, mp_sendrecv_im2, mp_sendrecv_im3,&
                      mp_sendrecv_lv, mp_sendrecv_lm2, mp_sendrecv_lm3,&
                      mp_sendrecv_rv, mp_sendrecv_rm2, mp_sendrecv_rm3,&
                      mp_sendrecv_dv, mp_sendrecv_dm2, mp_sendrecv_dm3,&
                      mp_sendrecv_cv, mp_sendrecv_cm2, mp_sendrecv_cm3,&
                      mp_sendrecv_zv, mp_sendrecv_zm2, mp_sendrecv_zm3
  END INTERFACE

  INTERFACE mp_isendrecv
     MODULE PROCEDURE mp_isendrecv_i, mp_isendrecv_im2, mp_isendrecv_iv,&
                      mp_isendrecv_l, mp_isendrecv_lm2, mp_isendrecv_lv,&
                      mp_isendrecv_r, mp_isendrecv_rm2, mp_isendrecv_rv,&
                      mp_isendrecv_d, mp_isendrecv_dm2, mp_isendrecv_dv,&
                      mp_isendrecv_c, mp_isendrecv_cm2, mp_isendrecv_cv,&
                      mp_isendrecv_z, mp_isendrecv_zm2, mp_isendrecv_zv
  END INTERFACE

  INTERFACE mp_isend
     MODULE PROCEDURE mp_isend_iv, mp_isend_im2, mp_isend_im3,&
                      mp_isend_lv, mp_isend_lm2, mp_isend_lm3,&
                      mp_isend_rv, mp_isend_rm2, mp_isend_rm3,&
                      mp_isend_dv, mp_isend_dm2, mp_isend_dm3,&
                      mp_isend_cv, mp_isend_cm2, mp_isend_cm3,&
                      mp_isend_zv, mp_isend_zm2, mp_isend_zm3
     MODULE PROCEDURE mp_isend_custom
  END INTERFACE

  INTERFACE mp_irecv
     MODULE PROCEDURE mp_irecv_iv, mp_irecv_im2, mp_irecv_im3,&
                      mp_irecv_lv, mp_irecv_lm2, mp_irecv_lm3,&
                      mp_irecv_rv, mp_irecv_rm2, mp_irecv_rm3,&
                      mp_irecv_dv, mp_irecv_dm2, mp_irecv_dm3,&
                      mp_irecv_cv, mp_irecv_cm2, mp_irecv_cm3,&
                      mp_irecv_zv, mp_irecv_zm2, mp_irecv_zm3
     MODULE PROCEDURE mp_irecv_custom
  END INTERFACE

  INTERFACE mp_win_create
     MODULE PROCEDURE mp_win_create_iv,&
                      mp_win_create_lv,&
                      mp_win_create_rv,&
                      mp_win_create_dv,&
                      mp_win_create_cv,&
                      mp_win_create_zv
  END INTERFACE

  INTERFACE mp_rget
     MODULE PROCEDURE mp_rget_iv,&
                      mp_rget_lv,&
                      mp_rget_rv,&
                      mp_rget_dv,&
                      mp_rget_cv,&
                      mp_rget_zv
  END INTERFACE

  INTERFACE mp_allocate
     MODULE PROCEDURE mp_allocate_i,&
                      mp_allocate_l,&
                      mp_allocate_r,&
                      mp_allocate_d,&
                      mp_allocate_c,&
                      mp_allocate_z
  END INTERFACE

  INTERFACE mp_deallocate
     MODULE PROCEDURE mp_deallocate_i,&
                      mp_deallocate_l,&
                      mp_deallocate_r,&
                      mp_deallocate_d,&
                      mp_deallocate_c,&
                      mp_deallocate_z
  END INTERFACE

  INTERFACE mp_type_make
     MODULE PROCEDURE mp_type_make_struct
     MODULE PROCEDURE mp_type_make_i, mp_type_make_l,&
                      mp_type_make_r, mp_type_make_d,&
                      mp_type_make_c, mp_type_make_z
  END INTERFACE

  INTERFACE mp_file_write_at
     MODULE PROCEDURE mp_file_write_at_ch, mp_file_write_at_chv,&
                      mp_file_write_at_i, mp_file_write_at_iv,&
                      mp_file_write_at_r, mp_file_write_at_rv,&
                      mp_file_write_at_d, mp_file_write_at_dv,&
                      mp_file_write_at_c, mp_file_write_at_cv,&
                      mp_file_write_at_z, mp_file_write_at_zv,&
                      mp_file_write_at_l, mp_file_write_at_lv
  END INTERFACE

  INTERFACE mp_file_write_at_all
     MODULE PROCEDURE mp_file_write_at_all_ch, mp_file_write_at_all_chv,&
                      mp_file_write_at_all_i, mp_file_write_at_all_iv,&
                      mp_file_write_at_all_l, mp_file_write_at_all_lv,&
                      mp_file_write_at_all_r, mp_file_write_at_all_rv,&
                      mp_file_write_at_all_d, mp_file_write_at_all_dv,&
                      mp_file_write_at_all_c, mp_file_write_at_all_cv,&
                      mp_file_write_at_all_z, mp_file_write_at_all_zv
  END INTERFACE

  INTERFACE mp_file_read_at
     MODULE PROCEDURE mp_file_read_at_ch, mp_file_read_at_chv,&
                      mp_file_read_at_i, mp_file_read_at_iv,&
                      mp_file_read_at_r, mp_file_read_at_rv,&
                      mp_file_read_at_d, mp_file_read_at_dv,&
                      mp_file_read_at_c, mp_file_read_at_cv,&
                      mp_file_read_at_z, mp_file_read_at_zv,&
                      mp_file_read_at_l, mp_file_read_at_lv
  END INTERFACE

  INTERFACE mp_file_read_at_all
     MODULE PROCEDURE mp_file_read_at_all_ch, mp_file_read_at_all_chv,&
                      mp_file_read_at_all_i, mp_file_read_at_all_iv,&
                      mp_file_read_at_all_l, mp_file_read_at_all_lv,&
                      mp_file_read_at_all_r, mp_file_read_at_all_rv,&
                      mp_file_read_at_all_d, mp_file_read_at_all_dv,&
                      mp_file_read_at_all_c, mp_file_read_at_all_cv,&
                      mp_file_read_at_all_z, mp_file_read_at_all_zv
  END INTERFACE

  INTERFACE mp_alloc_mem
     MODULE PROCEDURE mp_alloc_mem_i, mp_alloc_mem_l,&
          mp_alloc_mem_d, mp_alloc_mem_z,&
          mp_alloc_mem_r, mp_alloc_mem_c
  END INTERFACE

  INTERFACE mp_free_mem
     MODULE PROCEDURE mp_free_mem_i, mp_free_mem_l,&
          mp_free_mem_d, mp_free_mem_z,&
          mp_free_mem_r, mp_free_mem_c
  END INTERFACE

! Type declarations
  TYPE mp_indexing_meta_type
     INTEGER, DIMENSION(:), POINTER :: index, chunks
  END TYPE mp_indexing_meta_type

  TYPE mp_type_descriptor_type
     INTEGER :: type_handle
     INTEGER :: length
#if defined(__parallel)
     INTEGER(kind=mpi_address_kind) :: base
#endif
     INTEGER(kind=int_4),  DIMENSION(:), POINTER :: data_i
     INTEGER(kind=int_8),  DIMENSION(:), POINTER :: data_l
     REAL(kind=real_4),    DIMENSION(:), POINTER :: data_r
     REAL(kind=real_8),    DIMENSION(:), POINTER :: data_d
     COMPLEX(kind=real_4), DIMENSION(:), POINTER :: data_c
     COMPLEX(kind=real_8), DIMENSION(:), POINTER :: data_z
     TYPE(mp_type_descriptor_type), DIMENSION(:), POINTER :: subtype
     INTEGER :: vector_descriptor(2)
     LOGICAL :: has_indexing
     TYPE(mp_indexing_meta_type) :: index_descriptor
  END TYPE mp_type_descriptor_type




  ! type internally used to store message passing performance indicators
! *****************************************************************************
  TYPE mp_perf_type
    CHARACTER ( LEN = 20 ) :: name
    INTEGER :: count
    REAL (KIND=dp) :: msg_size
    REAL (KIND=dp) :: time
  END TYPE mp_perf_type

  INTEGER, PARAMETER :: MAX_PERF = 21

! *****************************************************************************
  TYPE mp_perf_env_type
     !private
     INTEGER :: ref_count, id_nr
     TYPE(mp_perf_type), DIMENSION(MAX_PERF) :: mp_perfs
  END TYPE mp_perf_env_type

! *****************************************************************************
  TYPE mp_perf_env_p_type
     TYPE(mp_perf_env_type), POINTER         :: mp_perf_env
  END TYPE mp_perf_env_p_type

  ! introduce a stack of mp_perfs, first index is the stack pointer, for convience is replacing
  INTEGER, PARAMETER :: max_stack_size = 10
  INTEGER            :: stack_pointer = 0
  ! target attribute needed as a hack around ifc 7.1 bug
  TYPE(mp_perf_env_p_type), DIMENSION(max_stack_size), TARGET :: mp_perf_stack

  CHARACTER(LEN=20), PARAMETER :: sname(MAX_PERF) =  &
   (/"MP_Group            ", "MP_Bcast            ", "MP_Allreduce        ", &
     "MP_Gather           ", "MP_Sync             ", "MP_Alltoall         ", &
     "MP_SendRecv         ", "MP_ISendRecv        ", "MP_Wait             ", &
     "MP_comm_split       ", "MP_ISend            ", "MP_IRecv            ", &
     "MP_Send             ", "MP_Recv             ", "MP_Memory           ", &
     "MP_Put              ", "MP_Get              ", "MP_Fence            ", &
     "MP_Win_Lock         ", "MP_Win_Create       ", "MP_Win_Free         "/)
#if defined(__parallel)
  REAL(KIND=dp) :: t_start, t_end
#endif

  ! we make some assumptions on the length of INTEGERS, REALS and LOGICALS
  INTEGER, PARAMETER :: intlen=BIT_SIZE ( 0 ) / 8
  INTEGER, PARAMETER :: reallen=8
  INTEGER, PARAMETER :: loglen=BIT_SIZE ( 0 ) / 8
  INTEGER, PARAMETER :: charlen=1
  INTEGER, SAVE, PRIVATE :: last_mp_perf_env_id=0

  ! external timing hooks
  ! this interface (with subroutines in it) musst to be defined right before
  ! the regular subroutines/functions - otherwise prettify.py will screw up.
  INTERFACE
    SUBROUTINE timeset_interface(routineN, handle)
    CHARACTER(LEN=*), INTENT(IN)             :: routineN
    INTEGER, INTENT(OUT)                     :: handle

    END SUBROUTINE timeset_interface
    SUBROUTINE timestop_interface(handle)
    INTEGER, INTENT(IN)                      :: handle

    END SUBROUTINE timestop_interface
  END INTERFACE

  ! assumed to be private...
  PROCEDURE(timeset_interface), POINTER, SAVE  :: mp_external_timeset  => NULL()
  PROCEDURE(timestop_interface), POINTER, SAVE :: mp_external_timestop => NULL()

CONTAINS

! *****************************************************************************
!> \brief initializes the system default communicator
!> \param mp_comm [output] : handle of the default communicator
!> \par History
!>      2.2004 created [Joost VandeVondele ]
!> \note
!>      should only be called once
! *****************************************************************************
  SUBROUTINE mp_world_init(mp_comm)
    INTEGER, INTENT(OUT)                     :: mp_comm
#if defined(__parallel)
    INTEGER                                  :: ierr
!$  INTEGER                                  :: provided_tsl
!$  LOGICAL                                  :: no_threading_support

    mp_external_timeset  => NULL()
    mp_external_timestop => NULL()

#if defined(__NO_MPI_THREAD_SUPPORT_CHECK)
    ! Hack that does not request or check MPI thread suppolt level.
    ! User asserts that the MPI library will work correctly with
    ! threads.
!
!$  no_threading_support = .TRUE.
#else
    ! Does the right thing when using OpenMP: requests that the MPI
    ! library supports funneled mode and verifies that the MPI library
    ! provides that support.
    !
    ! Developers: Only the master thread will ever make calls to the
    ! MPI library.
!
!$  no_threading_support = .FALSE.
#endif
!$  IF (no_threading_support) THEN
       CALL mpi_init ( ierr )
       IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_init @ mp_world_init" )
!$  ELSE
!$OMP MASTER
!$     CALL mpi_init_thread (MPI_THREAD_FUNNELED, provided_tsl, ierr)
!$     IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_init_thread @ mp_world_init" )
!$     IF (provided_tsl .LT. MPI_THREAD_FUNNELED) THEN
!$        CALL mp_stop (0, "MPI library does not support the requested level of threading (MPI_THREAD_FUNNELED).")
!$     ENDIF
!$OMP END MASTER
!$  ENDIF
    CALL mpi_errhandler_set ( MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_errhandler_set @ mp_world_init" )
    mp_comm = MPI_COMM_WORLD
    debug_comm_count = 1
#else
    mp_comm = 0
#endif
    CALL add_mp_perf_env()
  END SUBROUTINE mp_world_init

! *****************************************************************************
!> \brief re-create the system default communicator with a different MPI 
!>        rank order
!> \param mp_comm [output] : handle of the default communicator
!> \param mp_new_comm ...
!> \param ranks_order ...
!> \par History
!>      1.2012 created [ Christiane Pousa ]
!> \note
!>      should only be called once, st very begining of CP2K run
! *****************************************************************************
  SUBROUTINE mp_reordering(mp_comm,mp_new_comm,ranks_order)
    INTEGER, INTENT(IN)                      :: mp_comm
    INTEGER, INTENT(out)                     :: mp_new_comm
    INTEGER, DIMENSION(:), POINTER           :: ranks_order

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_reordering', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: newcomm, newgroup, oldgroup
#endif

    CALL mp_timeset(routineN,handle)
    ierr = 0
#if defined(__parallel)
    t_start = m_walltime()    

    CALL mpi_comm_group( mp_comm, oldgroup, ierr );
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_group @ mp_reordering" )
    CALL mpi_group_incl(oldgroup,SIZE(ranks_order),ranks_order,newgroup, ierr)    
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_group_incl @ mp_reordering" )

    CALL mpi_comm_create(mp_comm, newgroup, newcomm, ierr)
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_create @ mp_reordering" )

    CALL mpi_group_free(oldgroup, ierr)
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_group_free @ mp_reordering" )
    CALL mpi_group_free(newgroup, ierr)
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_group_free @ mp_reordering" )

    ! update the system default communicator
    mp_new_comm = newcomm
    debug_comm_count = debug_comm_count + 1

    t_end = m_walltime ( )
    CALL add_perf(perf_id=1,count=1,time=t_end-t_start)
#endif 
    CALL mp_timestop(handle)
  END SUBROUTINE mp_reordering

! *****************************************************************************
!> \brief finalizes the system default communicator
!> \par History
!>      2.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE mp_world_finalize()

#if defined(__parallel)
    INTEGER                                  :: ierr
    CALL mpi_barrier ( MPI_COMM_WORLD,ierr ) ! call mpi directly to avoid 0 stack pointer
    CALL rm_mp_perf_env()
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_barrier @ mp_world_finalize" )
    debug_comm_count = debug_comm_count - 1
    IF (debug_comm_count .NE. 0) THEN
       ! A bug, we're leaking or double-freeing communicators. Needs to be fixed where the leak happens.
       ! Memory leak checking might be helpful to locate the culprit 
       CALL mp_abort("mp_world_finalize: assert failed:  leaking communicators")
    ENDIF
    CALL mpi_finalize ( ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_finalize @ mp_world_finalize" )
#else
    CALL rm_mp_perf_env()
#endif

  END SUBROUTINE mp_world_finalize

! all the following routines should work for a given communicator, not MPI_WORLD

! *****************************************************************************
!> \brief start and stop the performance indicators
!>      for every call to start there has to be (exactly) one call to stop
!> \param perf_env ...
!> \par History
!>      2.2004 created [Joost VandeVondele]
!> \note
!>      can be used to measure performance of a sub-part of a program.
!>      timings measured here will not show up in the outer start/stops
!>      Doesn't need a fresh communicator
! *****************************************************************************
  SUBROUTINE add_mp_perf_env(perf_env)
    TYPE(mp_perf_env_type), OPTIONAL, &
      POINTER                                :: perf_env

    stack_pointer = stack_pointer + 1
    IF (stack_pointer > max_stack_size) THEN
       CALL mp_abort ( "stack_pointer too large : message_passing @ add_mp_perf_env" )
    ENDIF
    NULLIFY(mp_perf_stack(stack_pointer)%mp_perf_env)
    IF (PRESENT(perf_env)) THEN
       mp_perf_stack(stack_pointer)%mp_perf_env => perf_env
       IF (ASSOCIATED(perf_env)) CALL mp_perf_env_retain(perf_env)
    END IF
    IF (.NOT.ASSOCIATED(mp_perf_stack(stack_pointer)%mp_perf_env)) THEN
       CALL mp_perf_env_create(mp_perf_stack(stack_pointer)%mp_perf_env)
    END IF
  END SUBROUTINE add_mp_perf_env

! *****************************************************************************
!> \brief ...
!> \param perf_env ...
! *****************************************************************************
  SUBROUTINE mp_perf_env_create(perf_env)
    TYPE(mp_perf_env_type), OPTIONAL, &
      POINTER                                :: perf_env

    INTEGER                                  :: i, stat

    NULLIFY(perf_env)
    ALLOCATE(perf_env,stat=stat)
    IF (stat/=0) THEN
       CALL mp_abort ( "allocation failed in mp_perf_env_create")
    ENDIF
    last_mp_perf_env_id=last_mp_perf_env_id+1
    perf_env%id_nr=last_mp_perf_env_id
    perf_env%ref_count=1
    DO i = 1, MAX_PERF
       perf_env%mp_perfs(i) % name = sname ( i )
       perf_env%mp_perfs(i) % count = 0
       perf_env%mp_perfs(i) % msg_size = 0.0_dp
       perf_env%mp_perfs(i) % time = 0.0_dp
    END DO

  END SUBROUTINE mp_perf_env_create

! *****************************************************************************
!> \brief ...
!> \param perf_env ...
! *****************************************************************************
  SUBROUTINE mp_perf_env_release(perf_env)
    TYPE(mp_perf_env_type), POINTER          :: perf_env

    INTEGER                                  :: stat

    IF (ASSOCIATED(perf_env)) THEN
       IF (perf_env%ref_count<1) THEN
          CALL mp_abort(&
               "invalid ref_count: message_passing @ mp_perf_env_release")
       END IF
       perf_env%ref_count=perf_env%ref_count-1
       IF (perf_env%ref_count==0) THEN
          DEALLOCATE(perf_env,stat=stat)
          IF (stat/=0) THEN
             CALL mp_abort("deallocation error: message_passing @ mp_perf_env_release")
          END IF
       END IF
    END IF
    NULLIFY(perf_env)
  END SUBROUTINE mp_perf_env_release

! *****************************************************************************
!> \brief ...
!> \param perf_env ...
! *****************************************************************************
  SUBROUTINE mp_perf_env_retain(perf_env)
    TYPE(mp_perf_env_type), POINTER          :: perf_env

    IF (.NOT.ASSOCIATED(perf_env)) THEN
       CALL mp_abort("unassociated perf_env: message_passing @ mp_perf_env_retain")
    END IF
    IF (perf_env%ref_count<1) THEN
       CALL mp_abort("invalid ref_count: message_passing @ mp_perf_env_retain")
    END IF
    perf_env%ref_count=perf_env%ref_count+1
  END SUBROUTINE mp_perf_env_retain

!.. reports the performance counters for the MPI run
! *****************************************************************************
!> \brief ...
!> \param perf_env ...
!> \param iw ...
! *****************************************************************************
  SUBROUTINE mp_perf_env_describe ( perf_env, iw )
    TYPE(mp_perf_env_type), POINTER          :: perf_env
    INTEGER, INTENT(IN)                      :: iw

#if defined(__parallel)
    INTEGER                                  :: i
    REAL(KIND=dp)                            :: per, vol
#endif

    IF (.NOT.ASSOCIATED(perf_env)) THEN
       CALL mp_abort ("unassociated perf_env : message_passing @ mp_perf_env_describe" )
    ENDIF
    IF (perf_env%ref_count<1) THEN
       CALL mp_abort ("invalid perf_env%ref_count : message_passing @ mp_perf_env_describe" )
    ENDIF
#if defined(__parallel)
    IF ( iw > 0 ) THEN
       WRITE ( iw, '( /, 1X, 79("-") )' )
       WRITE ( iw, '( " -", 77X, "-" )' )
       WRITE ( iw, '( " -", 24X, A, 24X, "-" )' ) ' MESSAGE PASSING PERFORMANCE '
       WRITE ( iw, '( " -", 77X, "-" )' )
       WRITE ( iw, '( 1X, 79("-"), / )' )
       WRITE ( iw, '( A, A, A )' ) ' ROUTINE', '             CALLS ', &
            ' TOT TIME [s]  AVE VOLUME [Bytes]  PERFORMANCE [MB/s]'
       DO i = 1, MAX_PERF

          IF ( perf_env%mp_perfs( i ) % count > 0 ) THEN
             vol = perf_env%mp_perfs( i ) % msg_size / REAL ( perf_env%mp_perfs( i ) % count,KIND=dp)
             IF ( perf_env%mp_perfs( i ) % time > 0.0_dp ) THEN
                per = perf_env%mp_perfs( i ) % msg_size / perf_env%mp_perfs (i) % time * 1.e-6_dp
             ELSE
                per = 0.0_dp
             ENDIF
             IF ( vol < 1.0_dp ) THEN
                WRITE ( iw, '(1X,A15,T17,I10,T27,F14.3)' ) &
                     ADJUSTL ( perf_env%mp_perfs( i ) % name ), perf_env%mp_perfs( i ) % count, &
                     perf_env%mp_perfs( i ) % time
             ELSE
                WRITE ( iw, '(1X,A15,T17,I10,T27,F14.3,T50,F11.0,T69,F12.2)' ) &
                     ADJUSTL ( perf_env%mp_perfs( i ) % name ), perf_env%mp_perfs( i ) % count, &
                     perf_env%mp_perfs( i ) % time, vol, per
             END IF
          ENDIF

       END DO
       WRITE ( iw, '( 1X, 79("-"), / )' )
    END IF
#endif
  END SUBROUTINE mp_perf_env_describe

! *****************************************************************************
!> \brief ...
! *****************************************************************************
  SUBROUTINE rm_mp_perf_env ()
    IF (stack_pointer<1) THEN
       CALL mp_abort ( "no perf_env in the stack : message_passing @ rm_mp_perf_env" )
    ENDIF
    CALL mp_perf_env_release(mp_perf_stack(stack_pointer)%mp_perf_env)
    stack_pointer = stack_pointer - 1
  END SUBROUTINE rm_mp_perf_env

! *****************************************************************************
!> \brief ...
!> \retval res ...
! *****************************************************************************
  FUNCTION get_mp_perf_env () RESULT(res)
    TYPE(mp_perf_env_type), POINTER          :: res

    IF (stack_pointer<1) THEN
       CALL mp_abort ( "no perf_env in the stack : message_passing @ get_mp_perf_env" )
    ENDIF
    res => mp_perf_stack(stack_pointer)%mp_perf_env
  END FUNCTION get_mp_perf_env

! *****************************************************************************
!> \brief ...
!> \param scr ...
! *****************************************************************************
  SUBROUTINE describe_mp_perf_env(scr)
    INTEGER, INTENT(in)                      :: scr

    TYPE(mp_perf_env_type), POINTER          :: perf_env

    perf_env => get_mp_perf_env()
    CALL mp_perf_env_describe(perf_env, scr)
  END SUBROUTINE describe_mp_perf_env

! *****************************************************************************
!> \brief adds the performance informations of one call
!> \param perf_id ...
!> \param count ...
!> \param time ...
!> \param msg_size ...
!> \author fawzi
! *****************************************************************************
  SUBROUTINE add_perf(perf_id,count,time,msg_size)
    INTEGER, INTENT(in)                      :: perf_id
    INTEGER, INTENT(in), OPTIONAL            :: count
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: time
    INTEGER, INTENT(in), OPTIONAL            :: msg_size

#if defined(__parallel)
    TYPE(mp_perf_type), POINTER              :: mp_perf

    mp_perf => mp_perf_stack (stack_pointer)%mp_perf_env%mp_perfs( perf_id )
    IF (PRESENT(count)) THEN
       mp_perf%count = mp_perf%count + count
    END IF
    IF (PRESENT(time)) THEN
       mp_perf%time = mp_perf%time + time
    END IF
    IF (PRESENT(msg_size)) THEN
       mp_perf%msg_size = mp_perf%msg_size+REAL(msg_size,dp)
    END IF
#endif

  END SUBROUTINE add_perf

! *****************************************************************************
!> \brief globally stops all tasks, can optionally print a message.
!>       this is intended to be low level, most of CP2K should rather use cp_error_handling
!> \param message ...
! *****************************************************************************
  SUBROUTINE mp_abort (message)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: message

    INTEGER                                  :: ierr, numtask, taskid

    ierr = 0
#if defined(__parallel)
    CALL mp_environ ( numtask, taskid, MPI_COMM_WORLD )
#else
    numtask=1
    taskid=0
#endif
    IF (PRESENT(message)) THEN
       WRITE(default_output_unit,'(A)')         ' CP2K| '//TRIM(message)
    ENDIF
    WRITE(default_output_unit,'(A,I0)')   ' CP2K| Abnormal program termination, stopped by process number ',taskid
    CALL m_flush(default_output_unit)

#if defined(__parallel)
    CALL mpi_abort ( MPI_COMM_WORLD, 1, ierr )
#else
    CALL m_abort() ! uncomment if you want nice core dumps
#endif
    ! this routine never returns
    STOP 1
  END SUBROUTINE mp_abort

! *****************************************************************************
!> \brief stops *after an mpi error* translating the error code
!> \param ierr an error code * returned by an mpi call *
!> \param prg_code ...
!> \note
!>       this function is private to message_passing.F
! *****************************************************************************
  SUBROUTINE mp_stop ( ierr, prg_code )
    INTEGER, INTENT(IN)                      :: ierr
    CHARACTER(LEN=*)                         :: prg_code

#if defined(__parallel)
    INTEGER                                  :: istat, len
    CHARACTER(LEN=MPI_MAX_ERROR_STRING )     :: error_string
    CHARACTER(LEN=MPI_MAX_ERROR_STRING+512)  :: full_error
#else
    CHARACTER(LEN=512)                       :: full_error
#endif

#if defined(__parallel)
    CALL mpi_error_string ( ierr, error_string, len, istat )
    WRITE(full_error,'(A,I0,A)') ' MPI error ',ierr,' in '//TRIM(prg_code)//' : '//error_string(1:len)
#else
    WRITE(full_error,'(A,I0,A)') ' MPI error (!?) ',ierr,' in '//TRIM(prg_code)
#endif

    CALL mp_abort(full_error)

  END SUBROUTINE mp_stop

! *****************************************************************************
!> \brief synchronizes with a barrier a given group of mpi tasks
!> \param group mpi communicator
! *****************************************************************************
  SUBROUTINE mp_sync ( group )
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sync', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_barrier ( group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_barrier @ mp_sync" )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=5,count=1,time=t_end-t_start)
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_sync

! *****************************************************************************
!> \brief returns number of tasks and task id for a given mpi group
!>       simple and cartesian version.. recursive needed in case of failing mpi_comm_rank.
!> \param numtask ...
!> \param taskid ...
!> \param groupid mpi communicator
!> \note
!>         ..mp_world_setup is gone, use mp_environ instead (i.e. give a groupid explicitly)
! *****************************************************************************
  RECURSIVE SUBROUTINE mp_environ_l ( numtask, taskid, groupid )

    INTEGER, INTENT(OUT)                     :: numtask, taskid
    INTEGER, INTENT(IN)                      :: groupid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_l', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr    = 0
    CALL mp_timeset(routineN,handle)

    numtask = 1
    taskid  = 0
#if defined(__parallel)
    CALL mpi_comm_rank ( groupid, taskid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_rank @ mp_environ_l" )

    CALL mpi_comm_size ( groupid, numtask, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_environ_l" )
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_environ_l

! *****************************************************************************
!> \brief ...
!> \param numtask ...
!> \param dims ...
!> \param task_coor ...
!> \param groupid ...
! *****************************************************************************
  SUBROUTINE mp_environ_c ( numtask, dims, task_coor, groupid )

    INTEGER, INTENT(OUT)                     :: numtask, dims( 2 ), &
                                                task_coor( 2 )
    INTEGER, INTENT(IN)                      :: groupid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    LOGICAL, DIMENSION(2)                    :: periods
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)
    numtask = 1
    task_coor = 0
    dims = 1
#if defined(__parallel)
    CALL mpi_comm_size ( groupid, numtask, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_environ_c" )

    CALL mpi_cart_get ( groupid, 2, dims, periods, task_coor, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_environ_c" )
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_environ_c

! *****************************************************************************
!> \brief ...
!> \param comm ...
!> \param ndims ...
!> \param dims ...
!> \param task_coor ...
!> \param periods ...
! *****************************************************************************
  SUBROUTINE mp_environ_c2 ( comm, ndims, dims, task_coor, periods )

    INTEGER, INTENT(IN)                      :: comm, ndims
    INTEGER, INTENT(OUT)                     :: dims( ndims ), &
                                                task_coor( ndims )
    LOGICAL, INTENT(out)                     :: periods( ndims )

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_c2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr = 0
    CALL mp_timeset(routineN,handle)

    task_coor = 0
    dims = 1
    periods=.FALSE.
#if defined(__parallel)
    CALL mpi_cart_get ( comm, ndims, dims, periods, task_coor, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_environ_c" )
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_environ_c2

!..mp_cart_create
! *****************************************************************************
!> \brief ...
!> \param comm_old ...
!> \param ndims ...
!> \param dims ...
!> \param pos ...
!> \param comm_cart ...
! *****************************************************************************
  SUBROUTINE mp_cart_create ( comm_old, ndims, dims, pos, comm_cart )

    INTEGER, INTENT(IN)                      :: comm_old, ndims
    INTEGER, INTENT(INOUT)                   :: dims( : )
    INTEGER, INTENT(OUT)                     :: pos( : ), comm_cart

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, nodes
#if defined(__parallel)
    LOGICAL, DIMENSION(1:ndims)              :: period
    LOGICAL                                  :: reorder
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

    nodes = 0
    pos ( 1:ndims ) = -1
    comm_cart = comm_old
#if defined(__parallel)

    t_start = m_walltime ( )
    CALL mpi_comm_size ( comm_old, nodes, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_cart_create" )

    IF (ANY(dims == 0)) CALL mpi_dims_create(nodes,ndims,dims,ierr)
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_dims_create @ mp_cart_create" )

    ! FIX ME.  Quick hack to avoid problems with realspace grids for compilers
    ! like IBM that actually reorder the processors when creating the new
    ! communicator
    reorder = .FALSE.
    period = .TRUE.
    CALL mpi_cart_create ( comm_old, ndims, dims, period, reorder, comm_cart, &
         ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_create @ mp_cart_create" )

    IF (comm_cart /= MPI_COMM_NULL) THEN
       debug_comm_count = debug_comm_count + 1
       CALL mpi_cart_get ( comm_cart, ndims, dims, period, pos, ierr )
       IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_cart_create" )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=1,count=1)
    CALL add_perf(perf_id=2,time=t_end-t_start)
#else
    pos ( 1:ndims ) = 0
    dims = 1
    comm_cart = 0
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_cart_create

!..mp_cart_coords
! *****************************************************************************
!> \brief ...
!> \param comm ...
!> \param rank ...
!> \param coords ...
! *****************************************************************************
  SUBROUTINE mp_cart_coords ( comm, rank, coords)

    INTEGER, INTENT(IN)                      :: comm, rank
    INTEGER, DIMENSION(:), INTENT(OUT)       :: coords

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_coords', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, m

    ierr = 0
    CALL mp_timeset(routineN,handle)

    m = SIZE ( coords )
#if defined(__parallel)
    CALL mpi_cart_coords ( comm, rank, m, coords, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_coords @ mp_cart_coords" )
#else
    coords = 0
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_cart_coords

!..mp_comm_compare
! *****************************************************************************
!> \brief ...
!> \param comm1 ...
!> \param comm2 ...
!> \param res ...
! *****************************************************************************
  SUBROUTINE mp_comm_compare ( comm1, comm2, res)

    INTEGER, INTENT(IN)                      :: comm1, comm2
    INTEGER, INTENT(OUT)                     :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_compare', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, iout

    ierr = 0
    CALL mp_timeset(routineN,handle)

    iout = 0
    res = 0
#if defined(__parallel)
    CALL mpi_comm_compare ( comm1, comm2, iout, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_compare @ mp_comm_compare" )
    SELECT CASE(iout)
    CASE(MPI_IDENT)
       res = 0
    CASE(MPI_CONGRUENT)
       res = 1
    CASE(MPI_SIMILAR)
       res = 2
    CASE (MPI_UNEQUAL)
       res = 3
    CASE default
       res = 4
    END SELECT
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_comm_compare

!..mp_cart_sub
! *****************************************************************************
!> \brief ...
!> \param comm ...
!> \param rdim ...
!> \param sub_comm ...
! *****************************************************************************
  SUBROUTINE mp_cart_sub ( comm, rdim, sub_comm )

    INTEGER, INTENT(IN)                      :: comm
    LOGICAL, DIMENSION(:), INTENT(IN)        :: rdim
    INTEGER, INTENT(OUT)                     :: sub_comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_sub', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr = 0
    CALL mp_timeset(routineN,handle)

    sub_comm = 0
#if defined(__parallel)
    CALL mpi_cart_sub ( comm, rdim, sub_comm, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_sub @ mp_cart_sub" )
    debug_comm_count = debug_comm_count + 1
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_cart_sub

!..mp_comm_free
! *****************************************************************************
!> \brief ...
!> \param comm ...
! *****************************************************************************
  SUBROUTINE mp_comm_free ( comm )

    INTEGER, INTENT(INOUT)                   :: comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_free', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    CALL mpi_comm_free ( comm, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_free @ mp_comm_free" )
    debug_comm_count = debug_comm_count - 1
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_comm_free

!..mp_comm_dup
! *****************************************************************************
!> \brief ...
!> \param comm1 ...
!> \param comm2 ...
! *****************************************************************************
  SUBROUTINE mp_comm_dup ( comm1, comm2 )

    INTEGER, INTENT(IN)                      :: comm1
    INTEGER, INTENT(OUT)                     :: comm2

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_dup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    CALL mpi_comm_dup ( comm1, comm2, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_dup @ mp_comm_dup" )
    debug_comm_count = debug_comm_count + 1
#else
    comm2 = comm1
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_comm_dup

!..mp_rank_compare
! *****************************************************************************
!> \brief ...
!> \param comm1 ...
!> \param comm2 ...
!> \param rank ...
! *****************************************************************************
  SUBROUTINE mp_rank_compare ( comm1, comm2, rank )

    INTEGER, INTENT(IN)                      :: comm1, comm2
    INTEGER, DIMENSION(:), INTENT(OUT)       :: rank

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_rank_compare', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: g1, g2, i, n, n1, n2
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rin
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

    rank = 0
#if defined(__parallel)
    CALL mpi_comm_size ( comm1, n1, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_size @ mp_rank_compare" )
    CALL mpi_comm_size ( comm2, n2, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_size @ mp_rank_compare" )
    n = MAX ( n1, n2 )
    CALL mpi_comm_group ( comm1, g1, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_group @ mp_rank_compare" )
    CALL mpi_comm_group ( comm2, g2, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_group @ mp_rank_compare" )
    ALLOCATE ( rin ( 0 : n - 1 ), STAT = ierr )
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ mp_rank_compare" )
    DO i = 0, n-1
       rin ( i ) = i
    END DO
    CALL mpi_group_translate_ranks ( g1, n, rin, g2, rank, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, &
         "mpi_group_translate_rank @ mp_rank_compare" )
    CALL mpi_group_free(g1, ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "group_free @ mp_rank_compare" )
    CALL mpi_group_free(g2, ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "group_free @ mp_rank_compare" )
    DEALLOCATE ( rin, STAT = ierr )
    IF ( ierr /= 0 ) CALL mp_abort( "deallocate @ mp_rank_compare" )
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_rank_compare

!..mp_dims_create
! *****************************************************************************
!> \brief ...
!> \param nodes ...
!> \param dims ...
! *****************************************************************************
  SUBROUTINE mp_dims_create ( nodes, dims )

    INTEGER, INTENT(IN)                      :: nodes
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: dims

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_dims_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, ndim

    ierr = 0
    CALL mp_timeset(routineN,handle)

    ndim = SIZE(dims)
#if defined(__parallel)
    IF (ANY(dims == 0)) CALL mpi_dims_create(nodes,ndim,dims,ierr)
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_dims_create @ mp_dims_create" )
#else
    dims = 1
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_dims_create

!..mp_cart_rank
! *****************************************************************************
!> \brief ...
!> \param group ...
!> \param pos ...
!> \param rank ...
! *****************************************************************************
  SUBROUTINE mp_cart_rank ( group, pos, rank )
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:), INTENT(IN)        :: pos
    INTEGER, INTENT(OUT)                     :: rank

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_rank', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr = 0
    CALL mp_timeset(routineN,handle)


#if defined(__parallel)
    CALL mpi_cart_rank ( group, pos, rank, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_rank @ mp_cart_rank" )
#else
    rank = 0
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_cart_rank


! *****************************************************************************
!> \brief waits for completion of the given request
!> \param request ...
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
!> \note
!>      see isendrecv
! *****************************************************************************
  SUBROUTINE mp_wait(request)
    INTEGER, INTENT(inout)                   :: request

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_wait', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE (status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )

    CALL mpi_wait(request,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_wait @ mp_wait" )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=9,count=1,time=t_end-t_start)
    DEALLOCATE (status)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_wait

! *****************************************************************************
!> \brief waits for completion of the given requests
!> \param requests ...
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
!> \note
!>      see isendrecv
! *****************************************************************************
  SUBROUTINE mp_waitall_1(requests)
    INTEGER, DIMENSION(:), INTENT(inout)     :: requests

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitall_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: count
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: status
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    count = SIZE(requests)
    ALLOCATE (status(MPI_STATUS_SIZE,count))
    t_start = m_walltime ( )

    CALL mpi_waitall(count,requests,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_waitall @ mp_waitall_1" )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=9,count=1,time=t_end-t_start)
    DEALLOCATE (status)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_waitall_1

! *****************************************************************************
!> \brief waits for completion of the given requests
!> \param requests ...
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
! *****************************************************************************
  SUBROUTINE mp_waitall_2(requests)
    INTEGER, DIMENSION(:, :), INTENT(inout)  :: requests

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitall_2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: count
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: status
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)


#if defined(__parallel)
    count = SIZE(requests)
    ALLOCATE (status(MPI_STATUS_SIZE,count))
    t_start = m_walltime ( )

    CALL mpi_waitall_internal(count,requests,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_waitall @ mp_waitall_2" )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=9,count=1,time=t_end-t_start)
    DEALLOCATE (status)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_waitall_2

! *****************************************************************************
!> \brief wrapper needed to deal with interfaces as present in openmpi 1.8.1
!>        the issue is with the rank or requests
!> \param count ...
!> \param array_of_requests ...
!> \param array_of_statuses ...
!> \param ierr ...
!> \author Joost VandeVondele 
! *****************************************************************************
#if defined(__parallel)
  SUBROUTINE mpi_waitall_internal(count, array_of_requests, array_of_statuses, ierr)
    INTEGER, INTENT(in)                      :: count
    INTEGER, DIMENSION(count), INTENT(inout) :: array_of_requests
    INTEGER, DIMENSION(MPI_STATUS_SIZE, *), &
      INTENT(out)                            :: array_of_statuses
    INTEGER, INTENT(out)                     :: ierr

     CALL mpi_waitall(count,array_of_requests,array_of_statuses,ierr)

  END SUBROUTINE mpi_waitall_internal
#endif

! *****************************************************************************
!> \brief waits for completion of any of the given requests
!> \param requests ...
!> \param completed ...
!> \par History
!>      09.2008 created
!> \author Iain Bethune (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
! *****************************************************************************
  SUBROUTINE mp_waitany(requests, completed)
    INTEGER, DIMENSION(:), INTENT(inout)     :: requests
    INTEGER, INTENT(out)                     :: completed

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitany', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: count
    INTEGER                                  :: status(MPI_STATUS_SIZE)
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    count = SIZE(requests)
    t_start = m_walltime ( )

    CALL mpi_waitany(count,requests,completed,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_waitany @ mp_waitany" )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=9,count=1,time=t_end-t_start)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_waitany


! *****************************************************************************
!> \brief tests for completion of the given requests
!> \param requests ...
!> \param completed ...
!> \param flag ...
!> \par History
!>      08.2011 created
!> \author Iain Bethune
! *****************************************************************************
  SUBROUTINE mp_testany(requests, completed, flag)
    INTEGER, DIMENSION(:), INTENT(inout)  :: requests
    INTEGER, INTENT(out), OPTIONAL           :: completed
    LOGICAL, INTENT(out), OPTIONAL           :: flag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_testany', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
#if defined(__parallel)
    INTEGER                                  :: completed_l, count
    INTEGER                                  :: status(MPI_STATUS_SIZE)
    LOGICAL                                  :: flag_l
#endif

    ierr = 0

#if defined(__parallel)
    count = SIZE(requests)

    CALL mpi_testany_internal(count,requests,completed_l,flag_l,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_testany @ mp_testany" )

    IF (PRESENT(completed)) completed = completed_l
    IF (PRESENT(flag)) flag = flag_l

#endif
  END SUBROUTINE mp_testany
! *****************************************************************************
!> \brief wrapper needed to deal with interfaces as present in openmpi 1.8.1
!>        the issue is with the rank or requests
!> \param count ...
!> \param array_of_requests ...
!> \param index ...
!> \param flag ...
!> \param status ...
!> \param ierr ...
!> \author Joost VandeVondele 
! *****************************************************************************
#if defined(__parallel)
  SUBROUTINE mpi_testany_internal(count, array_of_requests, index, flag, status, ierr)
    INTEGER, INTENT(in)                      :: count
    INTEGER, DIMENSION(count), INTENT(inout) :: array_of_requests
    INTEGER, INTENT(out)                     :: index
    LOGICAL, INTENT(out)                     :: flag
    INTEGER, DIMENSION(MPI_STATUS_SIZE), &
      INTENT(out)                            :: status
    INTEGER, INTENT(out)                     :: ierr

  END SUBROUTINE mpi_testany_internal
#endif

! *****************************************************************************
!> \brief the direct way to split a communicator each color is a sub_comm,
!>        the rank order is accoring to the order in the orig comm
!> \param comm ...
!> \param sub_comm ...
!> \param color ...
!> \param key ...
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE mp_comm_split_direct(comm,sub_comm,color,key)
    INTEGER, INTENT(in)                      :: comm
    INTEGER, INTENT(OUT)                     :: sub_comm
    INTEGER, INTENT(in)                      :: color
    INTEGER, INTENT(in), OPTIONAL            :: key

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_split_direct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, my_key

    ierr = 0
    CALL mp_timeset(routineN,handle)

    my_key = 0
#if defined(__parallel)
    t_start = m_walltime ( )
    IF(PRESENT(key)) my_key = key
    CALL mpi_comm_split(comm,color,my_key,sub_comm,ierr)
    debug_comm_count = debug_comm_count + 1
    IF (ierr/=mpi_success) CALL mp_stop(ierr,routineN)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=10,count=1,time=t_end-t_start)
#else
    CALL mp_comm_dup(comm,sub_comm)
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_comm_split_direct
! *****************************************************************************
!> \brief splits the given communicator in group in subgroups trying to organize
!>      them in a way that the communication within each subgroup is
!>      efficent (but not necessarily the comunication between subgroups)
!> \param comm the mpi communicator that you want to split
!> \param sub_comm the communicator for the subgroup (created, needs to be freed later)
!> \param ngroups actual number of groups
!> \param group_distribution input  : allocated with array with the nprocs entries (0 .. nprocs-1)
!> \param subgroup_min_size the minimum size of the subgroup
!> \param n_subgroups the number of subgroups wanted
!> \param group_partition n_subgroups sized array containing the number of cpus wanted per group.
!>                         should match the total number of cpus (only used if present and associated) (0..ngroups-1)
!> \par History
!>      10.2003 created [fawzi]
!>      02.2004 modified [Joost VandeVondele]
!> \author Fawzi Mohamed
!> \note
!>      at least one of subgroup_min_size and n_subgroups is needed,
!>      the other default to the value needed to use most processors.
!>      if less cpus are present than needed for subgroup min size, n_subgroups,
!>      just one comm is created that contains all cpus
! *****************************************************************************
  SUBROUTINE mp_comm_split(comm,sub_comm,ngroups, group_distribution, &
       subgroup_min_size, n_subgroups, group_partition)
    INTEGER, INTENT(in)                      :: comm
    INTEGER, INTENT(out)                     :: sub_comm, ngroups
    INTEGER, DIMENSION(:), POINTER           :: group_distribution
    INTEGER, INTENT(in), OPTIONAL            :: subgroup_min_size, n_subgroups
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: group_partition

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_split', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, mepos, nnodes
#if defined(__parallel)
    INTEGER                                  :: color, i, j, k, &
                                                my_subgroup_min_size
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

    ! actual number of groups

    IF (.NOT. PRESENT(subgroup_min_size) .AND. .NOT. PRESENT(n_subgroups)) THEN
       CALL mp_abort(routineP//" missing arguments ")
    ENDIF
    IF (PRESENT(subgroup_min_size) .AND. PRESENT(n_subgroups)) THEN
       CALL mp_abort(routineP//" too many arguments ")
    ENDIF

    CALL mp_environ(nnodes,mepos,comm)

    IF (.NOT. ASSOCIATED(group_distribution)) THEN
       CALL mp_abort(routineP//" group_distribution not associated")
    ENDIF
    IF (LBOUND(group_distribution,1) .NE. 0 .OR. &
         UBOUND(group_distribution,1).NE.nnodes-1) THEN
       CALL mp_abort(routineP//" group_distribution wrong bounds")
    ENDIF

#if defined(__parallel)
    t_start = m_walltime ( )
    IF (PRESENT(subgroup_min_size)) THEN
       IF (subgroup_min_size<0 .OR. subgroup_min_size>nnodes) THEN
          CALL mp_abort(routineP//" subgroup_min_size too small or too large")
       ENDIF
       ngroups= nnodes / subgroup_min_size
       my_subgroup_min_size = subgroup_min_size
    ELSE  ! n_subgroups
       IF (n_subgroups<=0) THEN
          CALL mp_abort(routineP//" n_subgroups too small")
       ENDIF
       IF (nnodes/n_subgroups > 0) THEN ! we have a least one cpu per group
          ngroups = n_subgroups
       ELSE ! well, only one group then
          ngroups = 1
       ENDIF
       my_subgroup_min_size = nnodes / ngroups
    ENDIF
    DO i=0,nnodes-1
       group_distribution(i)=i / my_subgroup_min_size
       ! if part of the rest, join the last group
       IF ( group_distribution(i) >= ngroups ) group_distribution(i)=ngroups-1
    ENDDO
    ! even the user gave a partition, see if we can use it to overwrite this choice
    IF (PRESENT(group_partition)) THEN
       IF (ASSOCIATED(group_partition)) THEN
          IF (ALL(group_partition>0) .AND. SUM(group_partition).EQ.nnodes .AND. ngroups == SIZE(group_partition)) THEN
             k=0
             DO i=0,SIZE(group_partition)-1
                DO j=1,group_partition(i)
                   group_distribution(k)=i
                   k=k+1
                ENDDO
             ENDDO
          ELSE
             ! just ignore silently as we have reasonable defaults. Probably a warning would not be to bad
          ENDIF
       ENDIF
    ENDIF
    color=group_distribution(mepos)
    CALL mpi_comm_split(comm,color,0,sub_comm,ierr)
    debug_comm_count = debug_comm_count + 1
    IF (ierr/=mpi_success) CALL mp_stop(ierr,"in "//routineP//" split")

    t_end = m_walltime ( )
    CALL add_perf(perf_id=10,count=1,time=t_end-t_start)
#else
    CALL mp_comm_dup(comm,sub_comm)
    group_distribution(0)=0
    ngroups=1
#endif
    CALL mp_timestop(handle)

  END SUBROUTINE mp_comm_split

! *****************************************************************************
!> \brief Extract the processor name of a MPI process
!> \param host_name The MPI process host name
!> \author Christiane Ribeiro Pousa
!> \note
!> - 2011-11-23 [UB] Copied from old ma_message_passing and correctly set 
!>                   string length.
! *****************************************************************************
  SUBROUTINE mp_proc_name(host_name)
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_proc_name', &
         routineP = moduleN//':'//routineN
    CHARACTER(LEN=default_string_length)     :: host_name

#if defined(__parallel)
    INTEGER                                  :: ierr, len_name
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME)    :: hname
    CALL mpi_get_processor_name(hname,len_name,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_get_processor_name @ "//routineN )
    len_name = MIN (default_string_length, MPI_MAX_PROCESSOR_NAME)
    host_name(1:len_name) = hname(1:len_name)
#else
    host_name = "localhost"
#endif
  END SUBROUTINE mp_proc_name

! *****************************************************************************
!> \brief probes for an incomming message with any tag
!> \param[inout] source the source of the possible incomming message, 
!>        if MP_ANY_SOURCE it is a blocking one and return value is the source 
!>        of the next incomming message
!>        if source is a different value it is a non-blocking probe retuning 
!>        MP_ANY_SOURCE if there is no incomming message 
!> \param[in] comm the communicator 
!> \param[out] tag the tag of the incomming message
!> \author Mandes
! *****************************************************************************
  SUBROUTINE mp_probe(source, comm, tag)
    INTEGER                                  :: source
    INTEGER, INTENT(IN)                      :: comm
    INTEGER, INTENT(OUT)                     :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_probe', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER, DIMENSION(mp_status_size)       :: status_single
    LOGICAL                                  :: flag
#endif

!   ---------------------------------------------------------------------------

    CALL mp_timeset(routineN,handle)

    ierr = 0
#if defined(__parallel)
    IF(source .EQ. mp_any_source) THEN
      CALL mpi_probe(mp_any_source, mp_any_tag, comm, status_single, ierr)
      IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_probe @ mp_probe" )
      source = status_single(MPI_SOURCE)
      tag = status_single(MPI_TAG)
    ELSE
      flag = .FALSE.
      CALL mpi_iprobe(source, mp_any_tag, comm, flag, status_single, ierr)
      IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_iprobe @ mp_probe" )
      IF(flag.EQV..FALSE.) THEN
        source = mp_any_source
        tag = -1 !status_single(MPI_TAG) ! in case of flag==false status is undefined
      ELSE
        tag = status_single(MPI_TAG)
      END IF
    END IF
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_probe

! *****************************************************************************
! Here come the data routines with none of the standard data types.
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param msg ...
!> \param source ...
!> \param gid ...
! *****************************************************************************
  SUBROUTINE mp_bcast_b(msg,source,gid)
    LOGICAL                                  :: msg
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_LOGICAL,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ " // routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*loglen)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_bcast_b

! *****************************************************************************
!> \brief ...
!> \param msg ...
!> \param source ...
!> \param gid ...
! *****************************************************************************
  SUBROUTINE mp_bcast_bv(msg,source,gid)
    LOGICAL                                  :: msg(:)
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_bv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    CALL mp_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_LOGICAL,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*loglen)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_bcast_bv

! *****************************************************************************
!> \brief ...
!> \param msg ...
!> \param source ...
!> \param gid ...
! *****************************************************************************
  SUBROUTINE mp_bcast_av(msg,source,gid)
    CHARACTER(LEN=*)                         :: msg
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_av', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: i, msglen, numtask, taskid
    INTEGER, DIMENSION(:), ALLOCATABLE       :: imsg
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )

    CALL mp_environ (numtask, taskid, gid )
    IF (taskid==source) msglen = LEN_TRIM(msg)

    CALL mp_bcast(msglen, source, gid)
    ! this is a workaround to avoid problems on the T3E
    ! at the moment we have a data alignment error when trying to
    ! broadcats characters on the T3E (not always!)
    ! JH 19/3/99 on galileo
    ! CALL mpi_bcast(msg,msglen,MPI_CHARACTER,source,gid,ierr)
    ALLOCATE (imsg(1:msglen))
    DO i = 1, msglen
       imsg(i) = ICHAR(msg(i:i))
    END DO
    CALL mpi_bcast(imsg,msglen,MPI_INTEGER,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    msg = ""
    DO i = 1, msglen
       msg(i:i) = CHAR(imsg(i))
    END DO
    DEALLOCATE (imsg)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*charlen)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_bcast_av

! *****************************************************************************
!> \brief ...
!> \param msg ...
!> \param source ...
!> \param gid ...
! *****************************************************************************
  SUBROUTINE mp_bcast_am(msg,source,gid)
    CHARACTER(LEN=*)                         :: msg(:)
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_am', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: i, j, k, msglen, msgsiz, &
                                                numtask, taskid
    INTEGER, ALLOCATABLE                     :: imsg(:), imsglen(:)
#endif

    ierr = 0
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mp_environ (numtask, taskid, gid )
    msgsiz = SIZE(msg)
    ! Determine size of the minimum array of integers to bradcast the string
    ALLOCATE (imsglen(1:msgsiz))
    DO j = 1, msgsiz
       IF (taskid==source) imsglen(j) = LEN_TRIM(msg(j))
    END DO
    CALL mp_bcast(imsglen, source, gid)
    msglen = SUM(imsglen)
    ! this is a workaround to avoid problems on the T3E
    ! at the moment we have a data alignment error when trying to
    ! broadcats characters on the T3E (not always!)
    ! JH 19/3/99 on galileo
    ! CALL mpi_bcast(msg,msglen,MPI_CHARACTER,source,gid,ierr)
    ALLOCATE (imsg(1:msglen))
    k = 0
    DO j = 1, msgsiz
       DO i = 1, imsglen(j)
          k = k + 1
          imsg(k) = ICHAR(msg(j)(i:i))
       END DO
    END DO
    CALL mpi_bcast(imsg,msglen,MPI_INTEGER,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    msg = ""
    k = 0
    DO j = 1, msgsiz
       DO i = 1, imsglen(j)
          k = k + 1
          msg(j)(i:i) = CHAR(imsg(k))
       END DO
    END DO
    DEALLOCATE (imsg)
    DEALLOCATE (imsglen)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*charlen*msgsiz)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_bcast_am

! *****************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! *****************************************************************************
  SUBROUTINE mp_minloc_dv(msg,gid)
    REAL(kind=real_8), INTENT(INOUT)         :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_dv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
    REAL(kind=real_8), ALLOCATABLE           :: res( : )
#endif

    ierr = 0
    IF ("d" .EQ. "l" .AND. real_8 .EQ. int_8) THEN
       CALL mp_abort&
            ("Minimal location not available with long integers @ "//routineN)
    ENDIF
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_8_size)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_minloc_dv

! *****************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! *****************************************************************************
  SUBROUTINE mp_minloc_iv(msg,gid)
    INTEGER(KIND=int_4), INTENT(INOUT)       :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_iv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
    INTEGER(KIND=int_4), ALLOCATABLE         :: res( : )
#endif

    ierr = 0
    IF ("i" .EQ. "l" .AND. int_4 .EQ. int_8) THEN
       CALL mp_abort&
            ("Minimal location not available with long integers @ "//routineN)
    ENDIF
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2INTEGER,MPI_MINLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*int_4_size)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_minloc_iv

! *****************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! *****************************************************************************
  SUBROUTINE mp_minloc_lv(msg,gid)
    INTEGER(KIND=int_8), INTENT(INOUT)       :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_lv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
    INTEGER(KIND=int_8), ALLOCATABLE         :: res( : )
#endif

    ierr = 0
    IF ("l" .EQ. "l" .AND. int_8 .EQ. int_8) THEN
       CALL mp_abort&
            ("Minimal location not available with long integers @ "//routineN)
    ENDIF
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_INTEGER8,MPI_MINLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*int_8_size)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_minloc_lv

! *****************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! *****************************************************************************
  SUBROUTINE mp_minloc_rv(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)         :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
    REAL(kind=real_4), ALLOCATABLE           :: res( : )
#endif

    ierr = 0
    IF ("r" .EQ. "l" .AND. real_4 .EQ. int_8) THEN
       CALL mp_abort&
            ("Minimal location not available with long integers @ "//routineN)
    ENDIF
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2REAL,MPI_MINLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_minloc_rv

! *****************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! *****************************************************************************
  SUBROUTINE mp_maxloc_dv(msg,gid)
    REAL(kind=real_8), INTENT(INOUT)         :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_dv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
    REAL(kind=real_8), ALLOCATABLE           :: res( : )
#endif

    ierr = 0
    IF ("d" .EQ. "l" .AND. real_8 .EQ. int_8) THEN
       CALL mp_abort&
            ("Maximal location not available with long integers @ "//routineN)
    ENDIF
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res,STAT=ierr)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_8_size)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_maxloc_dv

! *****************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! *****************************************************************************
  SUBROUTINE mp_maxloc_iv(msg,gid)
    INTEGER(KIND=int_4), INTENT(INOUT)       :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_iv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
    INTEGER(KIND=int_4), ALLOCATABLE         :: res( : )
#endif

    ierr = 0
    IF ("i" .EQ. "l" .AND. int_4 .EQ. int_8) THEN
       CALL mp_abort&
            ("Maximal location not available with long integers @ "//routineN)
    ENDIF
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2INTEGER,MPI_MAXLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res,STAT=ierr)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*int_4_size)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_maxloc_iv

! *****************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! *****************************************************************************
  SUBROUTINE mp_maxloc_lv(msg,gid)
    INTEGER(KIND=int_8), INTENT(INOUT)       :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_lv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
    INTEGER(KIND=int_8), ALLOCATABLE         :: res( : )
#endif

    ierr = 0
    IF ("l" .EQ. "l" .AND. int_8 .EQ. int_8) THEN
       CALL mp_abort&
            ("Maximal location not available with long integers @ "//routineN)
    ENDIF
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_INTEGER8,MPI_MAXLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res,STAT=ierr)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*int_8_size)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_maxloc_lv

! *****************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! *****************************************************************************
  SUBROUTINE mp_maxloc_rv(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)         :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
#if defined(__parallel)
    INTEGER                                  :: msglen
    REAL(kind=real_4), ALLOCATABLE           :: res( : )
#endif

    ierr = 0
    IF ("r" .EQ. "l" .AND. real_4 .EQ. int_8) THEN
       CALL mp_abort&
            ("Maximal location not available with long integers @ "//routineN)
    ENDIF
    CALL mp_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2REAL,MPI_MAXLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res,STAT=ierr)
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_maxloc_rv

! *****************************************************************************
!> \brief Logical OR reduction
!> \param[in,out] msg         Datum to perform inclusive disjunction (input)
!>                            and resultant inclusive disjunction (output)
!> \param[in] gid             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce
! *****************************************************************************
  SUBROUTINE mp_sum_b(msg,gid)
    LOGICAL, INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen
#if defined(__parallel)
    LOGICAL                                  :: res
#endif

    CALL mp_timeset(routineN,handle)
    ierr = 0
    msglen = 1
#if defined(__parallel)
    CALL mpi_allreduce(msg,res,msglen,MPI_LOGICAL,MPI_LOR,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
#endif
    CALL mp_timestop(handle)
  END SUBROUTINE mp_sum_b

! *****************************************************************************
!> \brief Opens a file 
!> \param[in] groupid    message passing environment identifier
!> \param[out] fh        file handle (file storage unit)
!> \param[in] filepath   path to the file 
!> \param amode_status   access mode
!> \param info ...
!> \par MPI-I/O mapping  mpi_file_open 
!> \par STREAM-I/O mapping  OPEN
!> 
!> \param[in](optinal) info   info object 
!> \par History
!>      11.2012 created [Hossein Bani-Hashemian]  
! *****************************************************************************
  SUBROUTINE mp_file_open(groupid, fh, filepath, amode_status, info)
    INTEGER, INTENT(IN)                      :: groupid
    INTEGER, INTENT(OUT)                     :: fh
    CHARACTER(len=*), INTENT(IN)             :: filepath
    INTEGER, INTENT(IN)                      :: amode_status
    INTEGER, INTENT(IN), OPTIONAL            :: info

    INTEGER                                  :: ierr, istat
#if defined(__parallel)
    INTEGER                                  :: my_info
#else
    CHARACTER(LEN=10)                        :: fstatus
    LOGICAL                                  :: exists, is_open
#endif

    ierr = 0
    istat = 0
#if defined(__parallel)
    my_info = mpi_info_null
    IF (PRESENT(info)) my_info = info
    CALL mpi_file_open(groupid, filepath, amode_status, my_info, fh, ierr)
    CALL mpi_file_set_errhandler(fh, MPI_ERRORS_RETURN, ierr)
    IF (ierr .NE. 0) CALL mp_stop ( ierr, "mpi_file_set_errhandler @ mp_file_open" )
#else
    IF ((amode_status .EQ. file_amode_create) .OR. &
        (amode_status .EQ. file_amode_create+file_amode_wronly) .OR. &
        (amode_status .EQ. file_amode_create+file_amode_wronly+file_amode_excl) ) THEN
       fstatus = "UNKNOWN"
    ELSE
       fstatus = "OLD"
    END IF
    ! Get a new unit number
    DO fh=1,999
       INQUIRE (UNIT=fh,EXIST=exists,OPENED=is_open,IOSTAT=istat)
       IF (exists.AND.(.NOT.is_open).AND.(istat == 0)) EXIT
    END DO
    OPEN(UNIT=fh, FILE=filepath, STATUS=fstatus, ACCESS="STREAM")
#endif
  END SUBROUTINE mp_file_open

! *****************************************************************************
!> \brief Closes a file 
!> \param[in] fh   file handle (file storage unit)                   
!> \par MPI-I/O mapping   mpi_file_close 
!> \par STREAM-I/O mapping   CLOSE 
!>
!> \par History
!>      11.2012 created [Hossein Bani-Hashemian]
! *****************************************************************************
  SUBROUTINE mp_file_close(fh)
    INTEGER, INTENT(INOUT)                   :: fh

    INTEGER                                  :: ierr

    ierr = 0
#if defined(__parallel)
    CALL mpi_file_set_errhandler(fh, MPI_ERRORS_RETURN, ierr)
    CALL mpi_file_close(fh, ierr)
    IF (ierr .NE. 0) CALL mp_stop ( ierr, "mpi_file_set_errhandler @ mp_file_close" )
#else
    CLOSE(fh)
#endif
  END SUBROUTINE mp_file_close

! *****************************************************************************
!> \brief Returns the file size
!> \param[in] fh file handle (file storage unit)
!> \param[out] file_size  the file size                   
!> \par MPI-I/O mapping   mpi_file_get_size
!> \par STREAM-I/O mapping   INQUIRE 
!>
!> \par History
!>      12.2012 created [Hossein Bani-Hashemian]
! *****************************************************************************
  SUBROUTINE mp_file_get_size(fh,file_size)
    INTEGER, INTENT(IN)                      :: fh
    INTEGER(kind=file_offset), INTENT(OUT)   :: file_size

    INTEGER                                  :: ierr

    ierr = 0
#if defined(__parallel)
    CALL mpi_file_set_errhandler(fh, MPI_ERRORS_RETURN, ierr)
    CALL mpi_file_get_size(fh, file_size, ierr)
    IF (ierr .NE. 0) CALL mp_stop ( ierr, "mpi_file_set_errhandler @ mp_file_get_size" )
#else
    INQUIRE(UNIT=fh, SIZE=file_size)
#endif
  END SUBROUTINE mp_file_get_size

! *****************************************************************************
!> \brief (parallel) Blocking individual file write using explicit offsets
!>        (serial) Unformatted stream write
!> \param[in] fh     file handle (file storage unit) 
!> \param[in] offset file offset (position)
!> \param[in] msg    data to be writen to the file
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_write_at
!> \par STREAM-I/O mapping   WRITE   
!> \param[in](optional) msglen number of the elements of data                    
! *****************************************************************************
  SUBROUTINE mp_file_write_at_chv(fh, offset, msg, msglen)
    CHARACTER, INTENT(IN)                      :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_chv', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status

    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT(fh, offset, msg, msg_len, MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_chv @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_write_at_chv

! *****************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! *****************************************************************************
  SUBROUTINE mp_file_write_at_ch(fh, offset, msg)
    CHARACTER(LEN=*), INTENT(IN)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_ch', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status

    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT(fh, offset, msg, LEN(msg), MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_ch @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_write_at_ch

! *****************************************************************************
!> \brief (parallel) Blocking collective file write using explicit offsets
!>        (serial) Unformatted stream write
!> \param fh ...
!> \param offset ...
!> \param msg ...
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_write_at_all 
!> \par STREAM-I/O mapping   WRITE   
! *****************************************************************************
  SUBROUTINE mp_file_write_at_all_chv(fh, offset, msg, msglen)
    CHARACTER, INTENT(IN)                      :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_all_chv', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status

    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT_ALL(fh, offset, msg, msg_len, MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_all_chv @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_write_at_all_chv

! *****************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! *****************************************************************************
  SUBROUTINE mp_file_write_at_all_ch(fh, offset, msg)
    CHARACTER(LEN=*), INTENT(IN)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER(kind=file_offset), INTENT(IN)      :: offset
    
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_all_ch', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status

    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT_ALL(fh, offset, msg, LEN(msg), MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_all_ch @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_write_at_all_ch

! *****************************************************************************
!> \brief (parallel) Blocking individual file read using explicit offsets 
!>        (serial) Unformatted stream read
!> \param[in] fh     file handle (file storage unit)
!> \param[in] offset file offset (position)
!> \param[out] msg   data to be read from the file
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_read_at 
!> \par STREAM-I/O mapping   READ  
!> \param[in](optional) msglen  number of elements of data                      
! *****************************************************************************
  SUBROUTINE mp_file_read_at_chv(fh, offset, msg, msglen)
    CHARACTER, INTENT(OUT)                     :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_chv', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status

    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT(fh, offset, msg, msg_len, MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_chv @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_chv

! *****************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! *****************************************************************************
  SUBROUTINE mp_file_read_at_ch(fh, offset, msg)
    CHARACTER(LEN=*), INTENT(OUT)              :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_ch', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status

    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT(fh, offset, msg, LEN(msg), MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_ch @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_ch

! *****************************************************************************
!> \brief (parallel) Blocking collective file read using explicit offsets
!>        (serial) Unformatted stream read
!> \param fh ...
!> \param offset ...
!> \param msg ...
!> \param msglen ...
!> \par MPI-I/O mapping    mpi_file_read_at_all 
!> \par STREAM-I/O mapping   READ  
! *****************************************************************************
  SUBROUTINE mp_file_read_at_all_chv(fh, offset, msg, msglen)
    CHARACTER, INTENT(OUT)                     :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_all_chv', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status

    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT_ALL(fh, offset, msg, msg_len, MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_all_chv @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_all_chv

! *****************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! *****************************************************************************
  SUBROUTINE mp_file_read_at_all_ch(fh, offset, msg)
    CHARACTER(LEN=*), INTENT(OUT)              :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER(kind=file_offset), INTENT(IN)      :: offset

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_all_ch', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status

    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT_ALL(fh, offset, msg, LEN(msg), MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_all_ch @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_all_ch

! *****************************************************************************
!> \brief Returns the size of a data type in bytes
!> \param[in] type_descriptor  data type                   
!> \param[out] type_size       size of the data type
!> \par MPI mapping
!>      mpi_type_size 
!> 
! *****************************************************************************
  SUBROUTINE mp_type_size (type_descriptor, type_size)
    TYPE(mp_type_descriptor_type), &
      INTENT(IN)                             :: type_descriptor
    INTEGER, INTENT(OUT)                     :: type_size

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_type_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0
#if defined(__parallel)
    CALL MPI_TYPE_SIZE (type_descriptor%type_handle, type_size, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_type_size @ "//routineN)
#else 
    SELECT CASE (type_descriptor%type_handle)
    CASE (1)
         type_size = real_4_size
    CASE (3)
         type_size = real_8_size
    CASE (5)
         type_size = 2*real_4_size
    CASE (7)
         type_size = 2*real_8_size
    END SELECT
#endif
  END SUBROUTINE mp_type_size

! *****************************************************************************
!> \brief ...
!> \param subtypes ...
!> \param vector_descriptor ...
!> \param index_descriptor ...
!> \retval type_descriptor ...
! *****************************************************************************
  FUNCTION mp_type_make_struct (subtypes,&
       vector_descriptor, index_descriptor) &
       RESULT (type_descriptor)
    TYPE(mp_type_descriptor_type), &
      DIMENSION(:), INTENT(IN)               :: subtypes
    INTEGER, DIMENSION(2), INTENT(IN), &
      OPTIONAL                               :: vector_descriptor
    TYPE(mp_indexing_meta_type), &
      INTENT(IN), OPTIONAL                   :: index_descriptor
    TYPE(mp_type_descriptor_type)            :: type_descriptor

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_type_make_struct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, n
#if defined(__parallel)
    INTEGER(kind=mpi_address_kind), &
      ALLOCATABLE, DIMENSION(:)              :: displacements
#endif
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: lengths, old_types

    ierr = 0
    n = SIZE (subtypes)
    !type_descriptor%mpi_type_handle = MPI_DATATYPE_NULL
    type_descriptor%length = 1
#if defined(__parallel)
    CALL mpi_get_address (MPI_BOTTOM, type_descriptor%base, ierr)
    IF (ierr /= 0) CALL mp_abort("MPI_get_address @ "//routineN)
    ALLOCATE (displacements(n))
#endif
    type_descriptor%vector_descriptor(1:2) = 1
    type_descriptor%has_indexing = .FALSE.
    ALLOCATE (type_descriptor%subtype(n))
    type_descriptor%subtype(:) = subtypes(:)
    ALLOCATE (lengths(n), old_types(n))
    DO i = 1, SIZE (subtypes)
#if defined(__parallel)
       displacements(i) = subtypes(i)%base
#endif
       old_types(i) = subtypes(i)%type_handle
       lengths(i) = subtypes(i)%length
    ENDDO
#if defined(__parallel)
    CALL MPI_Type_create_struct(n,&
         lengths, displacements, old_types,&
         type_descriptor%type_handle, ierr)
    IF (ierr /= 0) CALL mp_abort("MPI_Type_create_struct @ "//routineN)
    CALL MPI_Type_commit (type_descriptor%type_handle, ierr)
    IF (ierr /= 0) CALL mp_abort("MPI_Type_commit @ "//routineN)
#endif
    IF (PRESENT (vector_descriptor) .OR. PRESENT (index_descriptor)) THEN
       CALL mp_abort (routineN//" Vectors and indices NYI")
    ENDIF
  END FUNCTION mp_type_make_struct

! *****************************************************************************
!> \brief ...
!> \param type_descriptor ...
! *****************************************************************************
  RECURSIVE SUBROUTINE mp_type_free (type_descriptor)
    TYPE(mp_type_descriptor_type), &
      INTENT(inout)                          :: type_descriptor

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_type_free', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr

    ierr = 0

    ! If the subtype is associated, then it's a user-defined data type.

    IF (ASSOCIATED (type_descriptor%subtype)) THEN
        DO i = 1, SIZE(type_descriptor%subtype)
          CALL mp_type_free (type_descriptor%subtype(i))
       ENDDO
       DEALLOCATE (type_descriptor%subtype)
#if defined(__parallel)
       CALL MPI_Type_free (type_descriptor%type_handle, ierr)
       IF (ierr /= 0) CALL mp_abort("MPI_Type_free @ "//routineN)
#endif
    ENDIF
  END SUBROUTINE mp_type_free

! *****************************************************************************
!> \brief Non-blocking send of custom type
!> \param msgin ...
!> \param dest ...
!> \param comm ...
!> \param request ...
!> \param tag ...
! *****************************************************************************
  SUBROUTINE mp_isend_custom(msgin,dest,comm,request,tag)
    TYPE(mp_type_descriptor_type), &
      INTENT(IN)                             :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_custom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, my_tag

    ierr = 0
    my_tag = 0

#if defined(__parallel)
    IF (PRESENT(tag)) my_tag=tag

    CALL mpi_isend(MPI_BOTTOM,1,msgin%type_handle,dest,my_tag,&
            comm,request,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )
#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
  END SUBROUTINE mp_isend_custom

! *****************************************************************************
!> \brief Non-blocking receive of vector data
!> \param msgout ...
!> \param source ...
!> \param comm ...
!> \param request ...
!> \param tag ...
! *****************************************************************************
  SUBROUTINE mp_irecv_custom(msgout,source,comm,request,tag)
    TYPE(mp_type_descriptor_type), &
      INTENT(INOUT)                          :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_custom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, my_tag

    ierr = 0
    my_tag = 0

#if defined(__parallel)
    IF (PRESENT(tag)) my_tag=tag

    CALL mpi_irecv(MPI_BOTTOM,1,msgout%type_handle,source, my_tag,&
         comm,request,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )
#else
    ierr=1
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
  END SUBROUTINE mp_irecv_custom

! *****************************************************************************
!> \brief Window free
!> \param win ...
! *****************************************************************************
  SUBROUTINE mp_win_free(win)
    INTEGER, INTENT(INOUT)                   :: win

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_free', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    t_start = m_walltime ( )
    
    CALL mpi_win_free(win,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_win_free @ "//routineN )
    
    t_end = m_walltime ( )
    CALL add_perf(perf_id=21,count=1,time=t_end-t_start)
#else
    CALL mp_abort( "mp_win_free called in non parallel case" )
#endif
  END SUBROUTINE mp_win_free

! *****************************************************************************
!> \brief Window lock
!> \param win ...
! *****************************************************************************
  SUBROUTINE mp_win_lock_all(win)
    INTEGER, INTENT(INOUT)                   :: win

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_lock_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    t_start = m_walltime ( )

#if __MPI_VERSION > 2
    CALL mpi_win_lock_all(MPI_MODE_NOCHECK,win,ierr)
#else
    CALL mp_abort("mp_win_lock_all requires MPI-3 standard")
#endif
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_win_lock_all @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=19,count=1,time=t_end-t_start)
#else
    CALL mp_abort( "mp_win_lock_all called in non parallel case" )
#endif
  END SUBROUTINE mp_win_lock_all

! *****************************************************************************
!> \brief Window lock
!> \param win ...
! *****************************************************************************
  SUBROUTINE mp_win_unlock_all(win)
    INTEGER, INTENT(INOUT)                   :: win

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_unlock_all', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    t_start = m_walltime ( )

#if __MPI_VERSION > 2
    CALL mpi_win_unlock_all(win,ierr)
#else
    CALL mp_abort("mp_win_unlock_all requires MPI-3 standard")
#endif
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_win_unlock_all @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=19,count=1,time=t_end-t_start)
#else
    CALL mp_abort( "mp_win_unlock_all called in non parallel case" )
#endif
  END SUBROUTINE mp_win_unlock_all

#include "message_passing_i.f90"
#include "message_passing_l.f90"
#include "message_passing_r.f90"
#include "message_passing_d.f90"
#include "message_passing_c.f90"
#include "message_passing_z.f90"


! *****************************************************************************
!> \brief Tests the MPI library
!> \param comm ...
!> \param npow ...
!> \par History
!>      JGH  6-Feb-2001 : Test and performance code
!> \author JGH  1-JAN-2001
!> \note
!>      quickly adapted benchmark code, will only work on an even number of CPUs.
!>      comm is the relevant, initialized communicator
! *****************************************************************************
  SUBROUTINE mpi_perf_test(comm,npow)

    INTEGER                                  :: comm, npow

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mpi_perf_test', &
      routineP = moduleN//':'//routineN

#if defined(__parallel)

    INTEGER :: I, ierr, ierror, itask, itests, J, jtask, left, nbufmax, &
      ncount, Ngrid, Nloc, nprocs, Ntot, partner, right, &
      status(MPI_STATUS_SIZE), taskid
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, sdispl
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: maxdiff, res, res2, res3, t1, &
                                                t2, t3, t4, t5
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buffer1, buffer2, buffer3, &
                                                lgrid, lgrid2, lgrid3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: grid, grid2, grid3, &
                                                send_timings, send_timings2
    REAL(KIND=dp), PARAMETER :: threshold=1.0E-8_dp


    ! set system sizes !
    ngrid= 10**npow

    CALL mpi_comm_rank(comm,taskid,ierror)
    CALL mpi_comm_size(comm,Nprocs,ierror)
    ionode=(taskid==0)
    IF (ionode) WRITE(*,*) "Running with ",nprocs
    IF (ionode) WRITE(*,*) "running messages with npow = ",npow
    IF (ionode) WRITE(*,*) "use MPI X in the input for larger (e.g. 6) of smaller (e.g. 3) messages"
    IF (MODULO(nprocs,2).NE.0) THEN
       WRITE(*,*) "Testing only with an even number of tasks"
       RETURN
    ENDIF

    ! equal loads
    Nloc=Ngrid/nprocs
    Ntot=Nprocs*Nloc
    nbufmax=10**npow
    !
    ALLOCATE(rcount(nprocs))
    ALLOCATE(scount(nprocs))
    ALLOCATE(sdispl(nprocs))
    ALLOCATE(rdispl(nprocs))
    ALLOCATE(buffer1(nbufmax))
    ALLOCATE(buffer2(nbufmax))
    ALLOCATE(buffer3(nbufmax))
    ALLOCATE(grid (Nloc,Nprocs))
    ALLOCATE(grid2(Nloc,Nprocs))
    ALLOCATE(grid3(Nloc,Nprocs))
    ALLOCATE(lgrid (Nloc))
    ALLOCATE(lgrid2(Nloc))
    ALLOCATE(lgrid3(Nloc))
    ALLOCATE(send_timings(0:nprocs-1,0:nprocs-1))
    ALLOCATE(send_timings2(0:nprocs-1,0:nprocs-1))
    buffer1=0.0_dp
    buffer2=0.0_dp
    buffer3=0.0_dp
    ! timings
    send_timings=0.0_dp
    send_timings2=0.0_dp
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ some in memory tests                   ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(*,*) "Testing in memory copies just 1 CPU "
    IF (ionode) WRITE(*,*) "  could tell something about the motherboard / cache / compiler "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          buffer2(1:ncount)=buffer1(1:ncount)
          t2=t2+MPI_WTIME()-t1 +threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(*,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ some in memory tests                   ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(*,*) "Testing in memory copies all cpus"
    IF (ionode) WRITE(*,*) "  is the memory bandwidth affected on an SMP machine ?"
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          buffer2(1:ncount)=buffer1(1:ncount)
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(*,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ first test point to point communication ---------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(*,*) "Testing truely point to point communication (i with j only)"
    IF (ionode) WRITE(*,*) "  is there some different connection between i j (e.g. shared memory comm)"
    ncount=10**npow
    IF (ionode) WRITE(*,*) "For messages of ",ncount*8," bytes"
    IF (ncount.GT.nbufmax) STOP
    DO itask=0,nprocs-1
       DO jtask=itask+1,nprocs-1
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          IF (taskid.EQ. itask) THEN
             CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, jtask, itask*jtask, comm, ierror)
          ENDIF
          IF (taskid.EQ. jtask) THEN
             CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, itask, itask*jtask, comm, status, ierror)
          ENDIF
          send_timings(itask,jtask)=MPI_WTIME()-t1+threshold
       ENDDO
    ENDDO
    send_timings2(:,:)=send_timings
    CALL MPI_REDUCE(send_timings2, send_timings, nprocs**2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
    IF (ionode) THEN
       DO itask=0,nprocs-1
          DO jtask=itask+1,nprocs-1
             WRITE(*,'(I4,I4,F12.4,A)') itask,jtask,ncount*8.0E-6_dp/send_timings(itask,jtask)," Mb/s"
          ENDDO
       ENDDO
    ENDIF
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ second test point to point communication -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(*,*) "Testing all nearby point to point communication (0,1)(2,3)..."
    IF (ionode) WRITE(*,*) "    these could / should all be on the same shared memory node "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          IF (MODULO(taskid,2)==0) THEN
             CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid+1, 0 , comm, ierror)
          ELSE
             CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid-1, 0 , comm, status, ierror)
          ENDIF
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(*,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    CALL MPI_BARRIER(comm,ierror)
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ third test point to point communication -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(*,*) "Testing all far point to point communication (0,nprocs/2),(1,nprocs/2+1),.."
    IF (ionode) WRITE(*,*) "    these could all be going over the network, and stress it a lot"
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          ! first half with partner
          IF (taskid .LT. nprocs/2) THEN
             CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid+nprocs/2, 0 , comm, ierror)
          ELSE
             CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid-nprocs/2, 0 , comm, status, ierror)
          ENDIF
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(*,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test root to all broadcast               -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(*,*) "Testing root to all broadcast "
    IF (ionode) WRITE(*,*) "    using trees at least ? "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL  MPI_BCAST(buffer1, ncount, MPI_DOUBLE_PRECISION, 0, comm, ierror)
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(*,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test mp_sum like behavior                -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(*,*) "Test global summation (mp_sum / mpi_allreduce) "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL  MPI_ALLREDUCE(buffer1,buffer2,ncount,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(*,'(I9,A,F12.4,A)') 8*ncount," Bytes ",(3**(npow-i))*ncount*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO
    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ test all to all communication            -------------------
    ! -------------------------------------------------------------------------------------------
    CALL MPI_BARRIER(comm,ierror)
    IF (ionode) WRITE(*,*) "Test all to all communication (mpi_alltoallv)"
    IF (ionode) WRITE(*,*) "    mpi/network getting confused ? "
    DO i=1,npow
       ncount=10**i
       t2=0.0E0_dp
       IF (ncount.GT.nbufmax) STOP
       scount=ncount/nprocs
       rcount=ncount/nprocs
       DO j=1,nprocs
          sdispl(j)=(j-1)*(ncount/nprocs)
          rdispl(j)=(j-1)*(ncount/nprocs)
       ENDDO
       DO j=1,3**(npow-i)
          CALL MPI_BARRIER(comm,ierror)
          t1=MPI_WTIME()
          CALL mpi_alltoallv ( buffer1, scount, sdispl, MPI_DOUBLE_PRECISION, &
               buffer2, rcount, rdispl, MPI_DOUBLE_PRECISION, comm, ierr )
          t2=t2+MPI_WTIME()-t1+threshold
       ENDDO
       CALL MPI_REDUCE(t2,t1,1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierror)
       IF (ionode) THEN
          WRITE(*,'(I9,A,F12.4,A)') 8*(ncount/nprocs)*nprocs," Bytes ",(3**(npow-i))*(ncount/nprocs)*nprocs*8.0E-6_dp/t1," Mb/s"
       ENDIF
    ENDDO

    ! -------------------------------------------------------------------------------------------
    ! ------------------------------ other stuff                            ---------------------
    ! -------------------------------------------------------------------------------------------
    IF (ionode) WRITE(*,*) " Clean tests completed "
    IF (ionode) WRITE(*,*) " Testing MPI_REDUCE scatter"
    rcount=Nloc
    DO itests=1,3
       IF (ionode) WRITE(*,*) "------------------------------- test ",itests," ------------------------"
       ! *** reference ***
       DO j=1,Nprocs
          DO i=1,Nloc
             grid(i,j)=MODULO(i*j*taskid,itests)
          ENDDO
       ENDDO
       t1=MPI_WTIME()
       CALL MPI_REDUCE_SCATTER(grid, lgrid, rcount, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
       t2=MPI_WTIME()-t1+threshold
       CALL mpi_allreduce(t2,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) WRITE(*,*) "MPI_REDUCE_SCATTER    ",res
       ! *** simple shift ***
       DO j=1,Nprocs
          DO i=1,Nloc
             grid2(i,j)=MODULO(i*j*taskid,itests)
          ENDDO
       ENDDO
       left =MODULO(taskid-1,Nprocs)
       right=MODULO(taskid+1,Nprocs)
       t3=MPI_WTIME()
       lgrid2=0.0E0_dp
       DO i=1,Nprocs
          lgrid2(:)=lgrid2+grid(:,MODULO(taskid-i,Nprocs)+1)
          IF (i.EQ.nprocs) EXIT
          CALL MPI_SENDRECV_REPLACE(lgrid2,nloc,MPI_DOUBLE_PRECISION,right,0,left,0,comm,status,ierr)
       ENDDO
       t4=MPI_WTIME()-t3+threshold
       CALL mpi_allreduce(t4,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       maxdiff=MAXVAL(ABS(lgrid2-lgrid))
       CALL mpi_allreduce(maxdiff,res2,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) WRITE(*,*) "MPI_SENDRECV_REPLACE  ",res,res2
       ! *** involved shift ****
       IF (MODULO(nprocs,2)/=0) STOP
       DO j=1,Nprocs
          DO i=1,Nloc
             grid3(i,j)=MODULO(i*j*taskid,itests)
          ENDDO
       ENDDO
       t3=MPI_WTIME()
       ! first sum the grid in pairs (0,1),(2,3) should be within an LPAR and fast XXXXXXXXX
       ! 0 will only need parts 0,2,4,... correctly summed
       ! 1 will only need parts 1,3,5,... correctly summed
       ! *** could nicely be generalised ****
       IF (MODULO(taskid,2)==0) THEN
          partner=taskid+1
          DO i=1,Nprocs,2 ! sum the full grid with the partner
             CALL MPI_SENDRECV(grid3(1,i+1),nloc,MPI_DOUBLE_PRECISION,partner,17, &
                  lgrid3,nloc,MPI_DOUBLE_PRECISION,partner,19,comm,status,ierr)
             grid3(:,i)=grid3(:,i)+lgrid3(:)
          ENDDO
       ELSE
          partner=taskid-1
          DO i=1,Nprocs,2
             CALL MPI_SENDRECV(grid3(1,i),nloc,MPI_DOUBLE_PRECISION,partner,19, &
                  lgrid3,nloc,MPI_DOUBLE_PRECISION,partner,17,comm,status,ierr)
             grid3(:,i+1)=grid3(:,i+1)+lgrid3(:)
          ENDDO
       ENDIF
       t4=MPI_WTIME()-t3+threshold
       ! now send a given buffer from 1 to 3 to 5 .. adding the right part of the data
       ! since we've summed an lgrid does only need to pass by even or odd tasks
       left =MODULO(taskid-2,Nprocs)
       right=MODULO(taskid+2,Nprocs)
       t3=MPI_WTIME()
       lgrid3=0.0E0_dp
       DO i=1,Nprocs,2
          lgrid3(:)=lgrid3+grid3(:,MODULO(taskid-i-1,Nprocs)+1)
          IF (i.EQ.nprocs-1) EXIT
          CALL MPI_SENDRECV_REPLACE(lgrid3,nloc,MPI_DOUBLE_PRECISION,right,0,left,0,comm,status,ierr)
       ENDDO
       t5=MPI_WTIME()-t3+threshold
       CALL mpi_allreduce(t4,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       CALL mpi_allreduce(t5,res2,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       maxdiff=MAXVAL(ABS(lgrid3-lgrid))
       CALL mpi_allreduce(maxdiff,res3,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm, ierr)
       IF (ionode) WRITE(*,*) "INVOLVED SHIFT        ",res+res2,"(",res,",",res2,")",res3
    ENDDO
    DEALLOCATE(rcount)
    DEALLOCATE(scount)
    DEALLOCATE(sdispl)
    DEALLOCATE(rdispl)
    DEALLOCATE(buffer1)
    DEALLOCATE(buffer2)
    DEALLOCATE(buffer3)
    DEALLOCATE(grid )
    DEALLOCATE(grid2)
    DEALLOCATE(grid3)
    DEALLOCATE(lgrid )
    DEALLOCATE(lgrid2)
    DEALLOCATE(lgrid3)
    DEALLOCATE(send_timings)
    DEALLOCATE(send_timings2)
#else
    WRITE(*,*) "No MPI tests for a serial program"
#endif
  END SUBROUTINE mpi_perf_test

! *****************************************************************************
!> \brief Starts a timer region
!> \param routineN ...
!> \param handle ...
! *****************************************************************************
  SUBROUTINE mp_timeset(routineN, handle)
    CHARACTER(len=*), INTENT(IN)             :: routineN
    INTEGER, INTENT(OUT)                     :: handle

    IF(ASSOCIATED(mp_external_timeset))&
       CALL mp_external_timeset(routineN, handle)
  END SUBROUTINE mp_timeset

! *****************************************************************************
!> \brief Ends a timer region
!> \param handle ...
! *****************************************************************************
  SUBROUTINE mp_timestop(handle)
    INTEGER, INTENT(IN)                      :: handle

    IF(ASSOCIATED(mp_external_timestop))&
       CALL mp_external_timestop(handle)
  END SUBROUTINE mp_timestop

END MODULE message_passing
