!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Interface to the message passing library MPI
!> \par History
!>     taken form message_passing.F
!> \author vw
! *****************************************************************************

MODULE dbcsr_message_passing
  USE dbcsr_c_mpi_calls,               ONLY: mp_alloc_mem,&
                                             mp_free_mem
  USE kinds,                           ONLY: int_4,&
                                             int_8,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size
  USE machine,                         ONLY: default_output_unit,&
                                             m_abort,&
                                             m_flush

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

#if defined(__parallel)
  INCLUDE "mpif.h"
#endif


#ifdef __parallel
  INTEGER, PARAMETER, PUBLIC :: mp_any_tag=MPI_ANY_TAG
  INTEGER, PARAMETER, PUBLIC :: mp_any_source=MPI_ANY_SOURCE
  INTEGER, PARAMETER, PUBLIC :: mp_comm_null=MPI_COMM_NULL
  INTEGER, PARAMETER, PUBLIC :: mp_request_null=MPI_REQUEST_NULL

  INTEGER, PARAMETER, PUBLIC :: file_offset=MPI_OFFSET_KIND
  INTEGER, PARAMETER, PUBLIC :: file_amode_create=MPI_MODE_CREATE
  INTEGER, PARAMETER, PUBLIC :: file_amode_rdonly=MPI_MODE_RDONLY
  INTEGER, PARAMETER, PUBLIC :: file_amode_wronly=MPI_MODE_WRONLY
  INTEGER, PARAMETER, PUBLIC :: file_amode_rdwr=MPI_MODE_RDWR
  INTEGER, PARAMETER, PUBLIC :: file_amode_excl=MPI_MODE_EXCL
#else
  INTEGER, PARAMETER, PUBLIC :: mp_any_tag=-1
  INTEGER, PARAMETER, PUBLIC :: mp_any_source=-2
  INTEGER, PARAMETER, PUBLIC :: mp_comm_null=-3
  INTEGER, PARAMETER, PUBLIC :: mp_request_null=-4

  INTEGER, PARAMETER, PUBLIC :: file_offset=int_8
  INTEGER, PARAMETER, PUBLIC :: file_amode_create=1
  INTEGER, PARAMETER, PUBLIC :: file_amode_rdonly=2
  INTEGER, PARAMETER, PUBLIC :: file_amode_wronly=4
  INTEGER, PARAMETER, PUBLIC :: file_amode_rdwr=8
  INTEGER, PARAMETER, PUBLIC :: file_amode_excl=64
#endif

  ! we need to fix this to a given number (crossing fingers) so that the serial code using Fortran stream IO and the MPI have the same sizes.
  INTEGER, PARAMETER, PUBLIC :: mpi_character_size=1
  INTEGER, PARAMETER, PUBLIC :: mpi_integer_size=4

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_message_passing'

  ! parameters that might be needed

  ! init and error
  PUBLIC :: mp_world_init, mp_world_finalize
  PUBLIC :: mp_abort

  ! informational / generation of sub comms
  PUBLIC :: mp_environ, mp_comm_compare, mp_cart_coords, mp_rank_compare
  PUBLIC :: mp_cart_create, mp_dims_create, mp_cart_rank, mp_cart_sub, mp_comm_free
  PUBLIC :: mp_comm_dup, mp_comm_split
  PUBLIC :: mp_group

  ! message passing
  PUBLIC :: mp_bcast, mp_sum, dmp_max, mp_maxloc, mp_minloc, mp_min, mp_sync
  PUBLIC :: mp_gather, mp_scatter, mp_alltoall, mp_sendrecv, mp_allgather
  PUBLIC :: mp_cart_shift, mp_isend, mp_irecv
  PUBLIC :: mp_sum_scatter, mp_shift, mp_isendrecv, mp_wait, mp_waitall, mp_waitany, mp_testany
  PUBLIC :: mp_gatherv
  PUBLIC :: mp_send, mp_recv

  ! RMA
  PUBLIC :: mp_window_type
  PUBLIC :: mp_window_create, mp_window_destroy
  PUBLIC :: mp_window_lock, mp_window_unlock
  PUBLIC :: mp_window_fence
  PUBLIC :: mp_rma_get, mp_rma_put

  ! Memory management
  PUBLIC :: mp_allocate, mp_deallocate

  ! Advanced Types
  PUBLIC :: mp_indexing_meta_type, mp_type_descriptor_type
  PUBLIC :: mp_type_make, mp_type_free
  PUBLIC :: mp_type_size

  ! default communicators
  PUBLIC :: MPI_COMM_SELF, MPI_COMM_WORLD, MPI_COMM_NULL

  ! I/O
  PUBLIC :: mp_file_open, mp_file_close
  PUBLIC :: mp_file_write_at, mp_file_read_at
  PUBLIC :: mp_file_write_at_all, mp_file_read_at_all
  PUBLIC :: mp_file_get_size

#ifndef __parallel
  INTEGER, PARAMETER :: MPI_COMM_SELF=0, MPI_COMM_WORLD=0, MPI_COMM_NULL=-1
#endif

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
     MODULE PROCEDURE mp_bcast_b!, mp_bcast_bv
!     MODULE PROCEDURE mp_bcast_av, mp_bcast_am
     MODULE PROCEDURE mp_bcast_charv
  END INTERFACE

  INTERFACE mp_sum
     MODULE PROCEDURE mp_sum_i, mp_sum_iv, mp_sum_im, mp_sum_im3,&
                      mp_sum_l, mp_sum_lv, mp_sum_lm, mp_sum_lm3,&
                      mp_sum_r, mp_sum_rv, mp_sum_rm, mp_sum_rm3,&
                      mp_sum_d, mp_sum_dv, mp_sum_dm, mp_sum_dm3,&
                      mp_sum_c, mp_sum_cv, mp_sum_cm, mp_sum_cm3,&
                      mp_sum_z, mp_sum_zv, mp_sum_zm, mp_sum_zm3,&
                      mp_sum_root_iv, mp_sum_root_im,&
                      mp_sum_root_lv, mp_sum_root_lm,&
                      mp_sum_root_rv, mp_sum_root_rm,&
                      mp_sum_root_dv, mp_sum_root_dm,&
                      mp_sum_root_cv, mp_sum_root_cm,&
                      mp_sum_root_zv, mp_sum_root_zm
     MODULE PROCEDURE mp_sum_b
  END INTERFACE

  ! dmp_max, this fixes a bug in g95
  ! -> ambigous interface for mp_max
  !    in module dbcsr_message_passing and wrongly
  !    exported from e.g. module timings
  INTERFACE dmp_max
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
     MODULE PROCEDURE mp_gather_i, mp_gather_iv,&
                      mp_gather_l, mp_gather_lv,&
                      mp_gather_r, mp_gather_rv,&
                      mp_gather_d, mp_gather_dv,&
                      mp_gather_c, mp_gather_cv,&
                      mp_gather_z, mp_gather_zv
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
     MODULE PROCEDURE mp_isendrecv_im2, mp_isendrecv_iv,&
                      mp_isendrecv_lm2, mp_isendrecv_lv,&
                      mp_isendrecv_rm2, mp_isendrecv_rv,&
                      mp_isendrecv_dm2, mp_isendrecv_dv,&
                      mp_isendrecv_cm2, mp_isendrecv_cv,&
                      mp_isendrecv_zm2, mp_isendrecv_zv
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

  INTERFACE mp_window_create
     MODULE PROCEDURE mp_win_create_i,&
                      mp_win_create_l,&
                      mp_win_create_r,&
                      mp_win_create_d,&
                      mp_win_create_c,&
                      mp_win_create_z
  END INTERFACE

  INTERFACE mp_rma_get
     MODULE PROCEDURE mp_rma_get_i,&
                      mp_rma_get_l,&
                      mp_rma_get_r,&
                      mp_rma_get_d,&
                      mp_rma_get_c,&
                      mp_rma_get_z
  END INTERFACE

  INTERFACE mp_rma_put
     MODULE PROCEDURE mp_rma_put_i,&
                      mp_rma_put_l,&
                      mp_rma_put_r,&
                      mp_rma_put_d,&
                      mp_rma_put_c,&
                      mp_rma_put_z
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
                      mp_file_write_at_z, mp_file_write_at_zv
  END INTERFACE

  INTERFACE mp_file_write_at_all
     MODULE PROCEDURE mp_file_write_at_all_ch, mp_file_write_at_all_chv,&
                      mp_file_write_at_all_i, mp_file_write_at_all_iv,&
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
                      mp_file_read_at_z, mp_file_read_at_zv
  END INTERFACE

  INTERFACE mp_file_read_at_all
     MODULE PROCEDURE mp_file_read_at_all_ch, mp_file_read_at_all_chv,&
                      mp_file_read_at_all_i, mp_file_read_at_all_iv,&
                      mp_file_read_at_all_r, mp_file_read_at_all_rv,&
                      mp_file_read_at_all_d, mp_file_read_at_all_dv,&
                      mp_file_read_at_all_c, mp_file_read_at_all_cv,&
                      mp_file_read_at_all_z, mp_file_read_at_all_zv
  END INTERFACE
! Type declarations

! *****************************************************************************
!> \brief Windows used for MPI RMA operations
!> \var id          MPI window identifier
!> \var info        MPI window info
! *****************************************************************************
  TYPE mp_window_type
#if defined(__parallel)
     INTEGER        :: id
#else
     INTEGER(KIND=int_4), DIMENSION(:), POINTER       :: src_i
     INTEGER(KIND=int_8), DIMENSION(:), POINTER       :: src_l
     REAL(KIND=real_4), DIMENSION(:), POINTER         :: src_r
     REAL(KIND=real_8), DIMENSION(:), POINTER         :: src_d
     COMPLEX(KIND=real_4), DIMENSION(:), POINTER      :: src_c
     COMPLEX(KIND=real_8), DIMENSION(:), POINTER      :: src_z
#endif
  END TYPE mp_window_type


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



CONTAINS

! *****************************************************************************
!> \brief initializes the system default communicator
!> \param mp_comm [output] : handle of the default communicator
!> \note
!>      should only be called once
!> \par History
!>      2.2004 created [Joost VandeVondele ]
! *****************************************************************************
  SUBROUTINE mp_world_init(mp_comm)
    INTEGER, INTENT(OUT)                     :: mp_comm
#if defined(__parallel)
    INTEGER                                  :: ierr
!$  INTEGER                                  :: provided_tsl

!   For now we call mpi_init instead of mpi_init_thread
!$  IF (.TRUE.) THEN
       CALL mpi_init ( ierr )
       IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_init @ mp_world_init" )
!$  ELSE
!$OMP MASTER
!$     CALL mpi_init_thread (MPI_THREAD_MULTIPLE, provided_tsl, ierr)
!$     IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_init_thread @ mp_world_init" )
!$     IF (provided_tsl .LT. MPI_THREAD_FUNNELED) THEN
!$        CALL mp_stop (0, "Inadequate level of thread support is provided.")
!$     ENDIF
!$OMP END MASTER
!$  ENDIF
    CALL mpi_errhandler_set ( MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_errhandler_set @ mp_world_init" )
    mp_comm = MPI_COMM_WORLD
#else
    mp_comm = 0
#endif
  END SUBROUTINE mp_world_init

! *****************************************************************************
!> \brief finalizes the system default communicator
!> \par History
!>      2.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE mp_world_finalize()
#if defined(__parallel)
    INTEGER                                  :: ierr
    CALL mpi_barrier ( MPI_COMM_WORLD,ierr ) ! call mpi directly to avoid 0 stack pointer
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_barrier @ mp_world_finalize" )
    CALL mpi_finalize ( ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_finalize @ mp_world_finalize" )
#endif
  END SUBROUTINE mp_world_finalize

! all the following routines should work for a given communicator, not MPI_WORLD

! *****************************************************************************
!> \brief globally stops all tasks, can optionally print a message.
!>       this is intended to be low level, most of CP2K should rather use cp_error_handling
! *****************************************************************************
  SUBROUTINE mp_abort (message)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: message

    INTEGER                                  :: ierr, numtask, taskid

    WRITE (*,*) message

#if defined(__parallel)
    CALL mp_environ ( numtask, taskid, MPI_COMM_WORLD )
#else
    numtask=1
    taskid=0
#endif
    IF (PRESENT(message)) THEN
       WRITE(default_output_unit,'(A)') ' libdbcsr| '//TRIM(message)
    ENDIF
    WRITE(default_output_unit,'(A,I0)') ' libdbcsr| Abnormal program termination, stopped by process number ',taskid
    CALL m_flush(default_output_unit)

#if defined(__parallel)
    CALL mpi_abort ( MPI_COMM_WORLD, 1, ierr )
#else
    CALL m_abort()
#endif
    ! this routine never returns
    STOP 1
  END SUBROUTINE mp_abort

! *****************************************************************************
!> \brief stops *after an mpi error* translating the error code
!> \param ierr an error code * returned by an mpi call *
!> \note
!>       this function is private to message_passing.F
! *****************************************************************************
  SUBROUTINE mp_stop ( ierr, prg_code )
    INTEGER, INTENT(IN)                      :: ierr
    CHARACTER(LEN=*)                         :: prg_code

    INTEGER                                  :: istat, len
#if defined(__parallel)
    CHARACTER ( LEN = MPI_MAX_ERROR_STRING ) :: error_string
    CHARACTER(LEN=MPI_MAX_ERROR_STRING+512)  :: full_error
#else
    CHARACTER ( LEN = 512 )                  :: error_string
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

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_barrier ( group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_barrier @ mp_sync" )
#endif

  END SUBROUTINE mp_sync

! *****************************************************************************
!> \brief returns number of tasks and task id for a given mpi group
!>       simple and cartesian version
!> \param groupid mpi communicator
!> \note
!>         ..mp_world_setup is gone, use mp_environ instead (i.e. give a groupid explicitly)
! *****************************************************************************
  SUBROUTINE mp_environ_l ( numtask, taskid, groupid )

    INTEGER, INTENT(OUT)                     :: numtask, taskid
    INTEGER, INTENT(IN)                      :: groupid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_l', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr    = 0
    numtask = 1
    taskid  = 0
#if defined(__parallel)
    CALL mpi_comm_rank ( groupid, taskid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_rank @ mp_environ_l" )

    CALL mpi_comm_size ( groupid, numtask, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_environ_l" )
#endif

  END SUBROUTINE mp_environ_l

! *****************************************************************************
  SUBROUTINE mp_environ_c ( numtask, dims, task_coor, groupid )

    INTEGER, INTENT(OUT)                     :: numtask, dims( 2 ), &
                                                task_coor( 2 )
    INTEGER, INTENT(IN)                      :: groupid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: periods( 2 )

    ierr = 0
    numtask = 1
    task_coor = 0
    dims = 1
#if defined(__parallel)
    CALL mpi_comm_size ( groupid, numtask, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_environ_c" )

    CALL mpi_cart_get ( groupid, 2, dims, periods, task_coor, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_environ_c" )
#endif

  END SUBROUTINE mp_environ_c

! *****************************************************************************
  SUBROUTINE mp_environ_c2 ( comm, ndims, dims, task_coor, periods )

    INTEGER, INTENT(IN)                      :: comm, ndims
    INTEGER, INTENT(OUT)                     :: dims( ndims ), &
                                                task_coor( ndims )
    LOGICAL, INTENT(out)                     :: periods( ndims )

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_c2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

    task_coor = 0
    dims = 1
    periods=.FALSE.
#if defined(__parallel)
    CALL mpi_cart_get ( comm, ndims, dims, periods, task_coor, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_environ_c" )
#endif

  END SUBROUTINE mp_environ_c2

!..mp_cart_create
! *****************************************************************************
  SUBROUTINE mp_cart_create ( comm_old, ndims, dims, pos, comm_cart )

    INTEGER, INTENT(IN)                      :: comm_old, ndims
    INTEGER, INTENT(INOUT)                   :: dims( : )
    INTEGER, INTENT(OUT)                     :: pos( : ), comm_cart

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, nodes
    LOGICAL                                  :: period(1:ndims), reorder

    ierr = 0
    pos ( 1:ndims ) = -1
    comm_cart = comm_old
#if defined(__parallel)

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
       CALL mpi_cart_get ( comm_cart, ndims, dims, period, pos, ierr )
       IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_cart_create" )
    END IF
#else
    pos ( 1:ndims ) = 0
    dims = 1
    comm_cart = 0
#endif

  END SUBROUTINE mp_cart_create

!..mp_cart_coords
! *****************************************************************************
  SUBROUTINE mp_cart_coords ( comm, rank, coords)

    INTEGER, INTENT(IN)                      :: comm, rank
    INTEGER, DIMENSION(:), INTENT(OUT)       :: coords

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_coords', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, m

    ierr = 0

    m = SIZE ( coords )
#if defined(__parallel)
    CALL mpi_cart_coords ( comm, rank, m, coords, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_coords @ mp_cart_coords" )
#else
    coords = 0
#endif

  END SUBROUTINE mp_cart_coords

!..mp_cart_shift
! *****************************************************************************
  SUBROUTINE mp_cart_shift ( comm, dir, disp, source, dest )

    INTEGER, INTENT(IN)                      :: comm, dir, disp
    INTEGER, INTENT(OUT)                     :: source, dest

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_shift', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_cart_shift ( comm, dir, disp, source, dest, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_shift @ mp_cart_shift" )
#else
    source = 0
    dest = 0
#endif

  END SUBROUTINE mp_cart_shift

!..mp_comm_compare
! *****************************************************************************
  SUBROUTINE mp_comm_compare ( comm1, comm2, RESULT)

    INTEGER, INTENT(IN)                      :: comm1, comm2
    INTEGER, INTENT(OUT)                     :: RESULT

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_compare', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, iout

    ierr = 0

    RESULT = 0
#if defined(__parallel)
    CALL mpi_comm_compare ( comm1, comm2, iout, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_compare @ mp_comm_compare" )
    SELECT CASE(iout)
    CASE(MPI_IDENT)
       RESULT=0
    CASE(MPI_CONGRUENT)
       RESULT=1
    CASE(MPI_SIMILAR)
       RESULT=2
    CASE (MPI_UNEQUAL)
       RESULT=3
    CASE default
       RESULT=4
    END SELECT
#endif

  END SUBROUTINE mp_comm_compare

!..mp_cart_sub
! *****************************************************************************
  SUBROUTINE mp_cart_sub ( comm, rdim, sub_comm )

    INTEGER, INTENT(IN)                      :: comm
    LOGICAL, DIMENSION(:), INTENT(IN)        :: rdim
    INTEGER, INTENT(OUT)                     :: sub_comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_sub', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

    sub_comm = 0
#if defined(__parallel)
    CALL mpi_cart_sub ( comm, rdim, sub_comm, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_sub @ mp_cart_sub" )
#endif

  END SUBROUTINE mp_cart_sub

!..mp_comm_free
! *****************************************************************************
  SUBROUTINE mp_comm_free ( comm )

    INTEGER, INTENT(IN)                      :: comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_free', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_comm_free ( comm, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_free @ mp_comm_free" )
#endif

  END SUBROUTINE mp_comm_free

!..mp_comm_dup
! *****************************************************************************
  SUBROUTINE mp_comm_dup ( comm1, comm2 )

    INTEGER, INTENT(IN)                      :: comm1
    INTEGER, INTENT(OUT)                     :: comm2

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_dup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_comm_dup ( comm1, comm2, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_dup @ mp_comm_dup" )
#else
    comm2=comm1
#endif

  END SUBROUTINE mp_comm_dup

!..mp_rank_compare
! *****************************************************************************
  SUBROUTINE mp_rank_compare ( comm1, comm2, rank )

    INTEGER, INTENT(IN)                      :: comm1, comm2
    INTEGER, DIMENSION(:), INTENT(OUT)       :: rank

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_rank_compare', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: g1, g2, i, ierr, n, n1, n2
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rin

    ierr = 0
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

  END SUBROUTINE mp_rank_compare

!..mp_dims_create
! *****************************************************************************
  SUBROUTINE mp_dims_create ( nodes, dims )

    INTEGER, INTENT(IN)                      :: nodes
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: dims

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_dims_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, ndim

    ierr = 0

#if defined(__parallel)
    ndim = SIZE ( dims )
    IF (ANY(dims == 0)) CALL mpi_dims_create(nodes,ndim,dims,ierr)
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_dims_create @ mp_dims_create" )
#else
    dims = 1
#endif

  END SUBROUTINE mp_dims_create

!..mp_cart_rank
! *****************************************************************************
  SUBROUTINE mp_cart_rank ( group, pos, rank )
    INTEGER, INTENT(IN)                      :: group
    INTEGER, DIMENSION(:), INTENT(IN)        :: pos
    INTEGER, INTENT(OUT)                     :: rank

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_rank', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_cart_rank ( group, pos, rank, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_rank @ mp_cart_rank" )
#else
    rank = 0
#endif

  END SUBROUTINE mp_cart_rank


! *****************************************************************************
!> \brief waits for completion of the given request
!> \note
!>      see isendrecv
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
! *****************************************************************************
  SUBROUTINE mp_wait(request)
    INTEGER, INTENT(inout)                   :: request

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_wait', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))

    CALL mpi_wait(request,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_wait @ mp_wait" )

    DEALLOCATE(status)
#endif
  END SUBROUTINE mp_wait

! *****************************************************************************
!> \brief waits for completion of the given requests
!> \note
!>      see isendrecv
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
! *****************************************************************************
  SUBROUTINE mp_waitall_1(requests)
    INTEGER, DIMENSION(:), INTENT(inout)     :: requests

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitall_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: count, ierr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: status

    ierr = 0

#if defined(__parallel)
    count=SIZE(requests)
    ALLOCATE(status(MPI_STATUS_SIZE,count))

    CALL mpi_waitall(count,requests,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_waitall @ mp_waitall_1" )

    DEALLOCATE(status)
#endif
  END SUBROUTINE mp_waitall_1

! *****************************************************************************
!> \brief waits for completion of the given requests
!> \note
!>      see isendrecv
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
! *****************************************************************************
  SUBROUTINE mp_waitall_2(requests)
    INTEGER, DIMENSION(:, :), INTENT(inout)  :: requests

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitall_2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: count, ierr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: status

    ierr = 0

#if defined(__parallel)
    count=SIZE(requests)
    ALLOCATE(status(MPI_STATUS_SIZE,count))

    CALL mpi_waitall(count,requests,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_waitall @ mp_waitall_2" )

    DEALLOCATE(status)
#endif
  END SUBROUTINE mp_waitall_2

! *****************************************************************************
!> \brief waits for completion of any of the given requests
!> \par History
!>      09.2008 created
!> \author Iain Bethune (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
! *****************************************************************************
  SUBROUTINE mp_waitany(requests, completed)
    INTEGER, DIMENSION(:), INTENT(inout)     :: requests
    INTEGER, INTENT(out)                     :: completed

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitany', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: count, ierr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: status

    ierr = 0

#if defined(__parallel)
    count=SIZE(requests)
    ALLOCATE(status(MPI_STATUS_SIZE,count))

    CALL mpi_waitany(count,requests,completed,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_waitany @ mp_waitany" )

    DEALLOCATE(status)
#endif
  END SUBROUTINE mp_waitany

! *****************************************************************************
!> \brief tests for completion of the given requests
!> \par History
!>      08.2011 created
!> \author Iain Bethune
! *****************************************************************************
  SUBROUTINE mp_testany(requests, completed, flag)
    INTEGER, DIMENSION(:, :), INTENT(inout)  :: requests
    INTEGER, INTENT(out), OPTIONAL           :: completed
    LOGICAL, INTENT(out), OPTIONAL           :: flag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_testany', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: completed_l, count, ierr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: status
    LOGICAL                                  :: flag_l

    ierr = 0

#if defined(__parallel)
    count=SIZE(requests)
    ALLOCATE(status(MPI_STATUS_SIZE,count))

    CALL mpi_testany(count,requests,completed_l,flag_l,status,ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_testany @ mp_testany" )

    IF (PRESENT(completed)) completed = completed_l
    IF (PRESENT(flag)) flag = flag_l

    DEALLOCATE(status)
#endif
  END SUBROUTINE mp_testany


! *****************************************************************************
!> \brief splits the given communicator in group in subgroups trying to organize
!>      them in a way that the communication within each subgroup is
!>      efficent (but not necessarily the comunication between subgroups)
!> \param comm the mpi communicator that you want to split
!> \param sub_comm the communicator for the subgroup (created, needs to be freed later)
!> \param ngroups actual number of groups
!> \param group_distribution input  : allocated with array with the nprocs entries (0 .. nprocs-1)
!> \param output the number of the group the proc belongs to (0..ngroups-1)
!> \param subgroup_min_size the minimum size of the subgroup
!> \param n_subgroups the number of subgroups wanted
!> \param group_partition n_subgroups sized array containing the number of cpus wanted per group.
!>                         should match the total number of cpus (only used if present and associated) (0..ngroups-1)
!> \note
!>      at least one of subgroup_min_size and n_subgroups is needed,
!>      the other default to the value needed to use most processors.
!>      if less cpus are present than needed for subgroup min size, n_subgroups,
!>      just one comm is created that contains all cpus
!> \par History
!>      10.2003 created [fawzi]
!>      02.2004 modified [Joost VandeVondele]
!> \author Fawzi Mohamed
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

    INTEGER                                  :: color, i, ierr, j, k, mepos, &
                                                my_subgroup_min_size, nnodes

    ierr = 0

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
    IF (PRESENT(subgroup_min_size)) THEN
       IF (subgroup_min_size<0) THEN
          CALL mp_abort(routineP//" subgroup_min_size too small")
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
    IF (ierr/=mpi_success) CALL mp_stop(ierr,"in "//routineP//" split")

#else
    CALL mp_comm_dup(comm,sub_comm)
    group_distribution(0)=0
    ngroups=1
#endif

  END SUBROUTINE mp_comm_split

! a bunch of currently unused features
!..mp_group
! *****************************************************************************
  SUBROUTINE mp_group(group_list,group_size,base_group,groupid)

    INTEGER, INTENT(IN)                      :: group_list( : ), group_size, &
                                                base_group
    INTEGER, INTENT(OUT)                     :: groupid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_group', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: base, ierr, newgroup

    ierr = 0

    ierr = 0
    groupid = base_group
#if defined(__parallel)
    CALL mpi_comm_group ( base_group, base, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_group @ mp_group" )

    CALL mpi_group_incl ( base, group_size, group_list, newgroup, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_group_incl @ mp_group" )

    CALL mpi_comm_create ( base_group, newgroup, groupid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_create @ mp_group" )
#endif

  END SUBROUTINE mp_group



! *****************************************************************************
!> \brief Destroys an MPI RMA window.
!> \author UB
!> \param[in] window       window to free and destroy
! *****************************************************************************
  SUBROUTINE mp_window_destroy(window)
    TYPE(mp_window_type), INTENT(INOUT) :: window
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_window_destroy', &
         routineP = moduleN//':'//routineN
    INTEGER                   :: ierr
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    ierr = 0
    CALL MPI_WIN_FREE(window, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_free @ mp_window_destroy")
#else
    NULLIFY(window%src_d, window%src_i)
#endif
  END SUBROUTINE mp_window_destroy


! *****************************************************************************
!> \brief Initiates an RMA access
!> \author UB
!> \param window       RMA access window
!> \param[in] remote   process to which access is initiated (target)
!> \param[in] shared   (optional) access is shared (default is exclusive)
!> \param[in] nocheck  (optional) don't check for others' locks (default is
!>                     to be careful and check)
! *****************************************************************************
  SUBROUTINE mp_window_lock(window, remote, shared, nocheck)
    TYPE(mp_window_type), INTENT(INOUT) :: window
    INTEGER, INTENT(IN)      :: remote
    LOGICAL, INTENT(IN), OPTIONAL       :: shared, nocheck
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_window_lock', &
         routineP = moduleN//':'//routineN
    INTEGER                   :: ierr
    INTEGER                   :: exclusivity, assertion
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    ierr = 0
    exclusivity = MPI_LOCK_EXCLUSIVE
    IF (PRESENT (shared)) THEN
       IF (shared) exclusivity = MPI_LOCK_SHARED
    ENDIF
    assertion = 0
    IF (PRESENT (nocheck)) THEN
       IF (nocheck) assertion = MPI_MODE_NOCHECK
    ENDIF
    CALL MPI_WIN_LOCK (exclusivity, remote, assertion, window%id, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_lock @ mp_window_lock")
#endif
  END SUBROUTINE mp_window_lock

! *****************************************************************************
!> \brief Terminates an RMA access
!> \author UB
!> \param[in] window   RMA access window
!> \param[in] remote   process to which access is initiated (target)
! *****************************************************************************
  SUBROUTINE mp_window_unlock(window, remote)
    TYPE(mp_window_type), INTENT(INOUT) :: window
    INTEGER, INTENT(IN)      :: remote
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_window_unlock', &
         routineP = moduleN//':'//routineN
    INTEGER                   :: ierr
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    ierr = 0
    CALL MPI_WIN_UNLOCK (remote, window%id, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_unlock @ mp_window_unlock")
#endif
  END SUBROUTINE mp_window_unlock


! *****************************************************************************
!> \brief Creates an RMA window epoch
!> \author UB
!> \param[in] window       RMA access window
!> \param[in] nothing_stored  (optional) asserts that there were no local
!>                            window updates since last sync (default is
!>                            .FALSE.)
!> \param[in] no_put_follows  (optional) asserts that nobody will put or
!>                            accumulate data in local data of window until
!>                            next fence (default is .FALSE.)
!> \param[in] no_pre_access   (optional) asserts that no RMA access calls
!>                            were made. All procesess must specify the same
!>                            flag (default is .FALSE.)
!> \param[in] no_post_access  (optional) asserts that no RMA access will be
!>                            issued. All processess must specify the same
!>                            flag (default is .FALSE.)
! *****************************************************************************
  SUBROUTINE mp_window_fence(window, nothing_stored, no_put_follows,&
       no_pre_access, no_post_access)
    TYPE(mp_window_type), INTENT(INOUT) :: window
    LOGICAL, INTENT(IN), OPTIONAL       :: nothing_stored, no_put_follows,&
                                           no_pre_access, no_post_access
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_window_unlock', &
         routineP = moduleN//':'//routineN
    INTEGER                   :: ierr, assertion
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    ierr = 0
    assertion = 0
    IF (PRESENT (nothing_stored)) THEN
       IF (nothing_stored) assertion = assertion + MPI_MODE_NOSTORE
    ENDIF
    IF (PRESENT (no_put_follows)) THEN
       IF (no_put_follows) assertion = assertion + MPI_MODE_NOPUT
    ENDIF
    IF (PRESENT (no_pre_access)) THEN
       IF (no_pre_access) assertion = assertion + MPI_MODE_NOPRECEDE
    ENDIF
    IF (PRESENT (no_post_access)) THEN
       IF (no_post_access) assertion = assertion + MPI_MODE_NOSUCCEED
    ENDIF
    CALL MPI_WIN_FENCE (assertion, window%id, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_fence @ mp_window_fence")
#endif
  END SUBROUTINE mp_window_fence


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

!   ---------------------------------------------------------------------------

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


  RECURSIVE SUBROUTINE mp_type_free (type_descriptor)
    TYPE(mp_type_descriptor_type), &
      INTENT(inout)                          :: type_descriptor

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_type_free', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr

!   ---------------------------------------------------------------------------
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
!> \brief Broadcasts a datum to all processes.
!> \par MPI mapping
!>      mpi_bcast
!> \param[in] msg             Datum to broadcast
!> \param[in] source          Processes which broadcasts
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_bcast_b(msg,source,gid)
    LOGICAL                                  :: msg
    INTEGER, INTENT(IN)                      :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
    msglen = 1
#if defined(__parallel)
    CALL mpi_bcast(msg,msglen,MPI_LOGICAL,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
#endif
  END SUBROUTINE mp_bcast_b

! *****************************************************************************
!> \brief Broadcasts a datum to all processes.
!> \par MPI mapping
!>      mpi_bcast
!> \param[in] msg             Datum to broadcast
!> \param[in] source          Processes which broadcasts
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_bcast_charv(msg,source,gid)
    CHARACTER(LEN=*), DIMENSION(:)           :: msg
    INTEGER, INTENT(IN)                      :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_charv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
    msglen = SIZE(msg)*LEN(msg)
#if defined(__parallel)
    CALL mpi_bcast(msg,msglen,MPI_CHARACTER,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
#endif
  END SUBROUTINE mp_bcast_charv


! *****************************************************************************
!> \brief Logical OR reduction
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Datum to perform inclusive disjunction (input)
!>                            and resultant inclusive disjunction (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sum_b(msg,gid)
    LOGICAL, INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    LOGICAL                                  :: res

    ierr = 0
    msglen = 1
#if defined(__parallel)
    CALL mpi_allreduce(msg,res,msglen,MPI_LOGICAL,MPI_LOR,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
#endif
  END SUBROUTINE mp_sum_b


! *****************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_minloc_dv(msg,gid)
    REAL(kind=real_8), INTENT(INOUT)         :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_dv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    REAL(kind=real_8), ALLOCATABLE           :: res( : )

    ierr = 0
    IF ("d" .EQ. "l" .AND. real_8 .EQ. int_8) THEN
       CALL mp_abort&
            ("Minimal location not available with long integers @ "//routineN)
    ENDIF
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
#endif
  END SUBROUTINE mp_minloc_dv

! *****************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_minloc_iv(msg,gid)
    INTEGER(KIND=int_4), INTENT(INOUT)       :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_iv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    INTEGER(KIND=int_4), ALLOCATABLE         :: res( : )

    ierr = 0
    IF ("i" .EQ. "l" .AND. int_4 .EQ. int_8) THEN
       CALL mp_abort&
            ("Minimal location not available with long integers @ "//routineN)
    ENDIF
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2INTEGER,MPI_MINLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
#endif
  END SUBROUTINE mp_minloc_iv
! *****************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_minloc_lv(msg,gid)
    INTEGER(KIND=int_8), INTENT(INOUT)       :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_lv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    INTEGER(KIND=int_8), ALLOCATABLE         :: res( : )

    ierr = 0
    IF ("l" .EQ. "l" .AND. int_8 .EQ. int_8) THEN
       CALL mp_abort&
            ("Minimal location not available with long integers @ "//routineN)
    ENDIF
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_INTEGER8,MPI_MINLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
#endif
  END SUBROUTINE mp_minloc_lv

! *****************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_minloc_rv(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)         :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    REAL(kind=real_4), ALLOCATABLE           :: res( : )

    ierr = 0
    IF ("r" .EQ. "l" .AND. real_4 .EQ. int_8) THEN
       CALL mp_abort&
            ("Minimal location not available with long integers @ "//routineN)
    ENDIF
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2REAL,MPI_MINLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
#endif
  END SUBROUTINE mp_minloc_rv


! *****************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_maxloc_dv(msg,gid)
    REAL(kind=real_8), INTENT(INOUT)         :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_dv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    REAL(kind=real_8), ALLOCATABLE           :: res( : )

    ierr = 0
    IF ("d" .EQ. "l" .AND. real_8 .EQ. int_8) THEN
       CALL mp_abort&
            ("Maximal location not available with long integers @ "//routineN)
    ENDIF
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res,STAT=ierr)
#endif
  END SUBROUTINE mp_maxloc_dv


! *****************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_maxloc_iv(msg,gid)
    INTEGER(KIND=int_4), INTENT(INOUT)       :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_iv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    INTEGER(KIND=int_4), ALLOCATABLE         :: res( : )

    ierr = 0
    IF ("i" .EQ. "l" .AND. int_4 .EQ. int_8) THEN
       CALL mp_abort&
            ("Maximal location not available with long integers @ "//routineN)
    ENDIF
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2INTEGER,MPI_MAXLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res,STAT=ierr)
#endif
  END SUBROUTINE mp_maxloc_iv

! *****************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_maxloc_lv(msg,gid)
    INTEGER(KIND=int_8), INTENT(INOUT)       :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_lv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    INTEGER(KIND=int_8), ALLOCATABLE         :: res( : )

    ierr = 0
    IF ("l" .EQ. "l" .AND. int_8 .EQ. int_8) THEN
       CALL mp_abort&
            ("Maximal location not available with long integers @ "//routineN)
    ENDIF
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_INTEGER8,MPI_MAXLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res,STAT=ierr)
#endif
  END SUBROUTINE mp_maxloc_lv

! *****************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_maxloc_rv(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)         :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    REAL(kind=real_4), ALLOCATABLE           :: res( : )

    ierr = 0
    IF ("r" .EQ. "l" .AND. real_4 .EQ. int_8) THEN
       CALL mp_abort&
            ("Maximal location not available with long integers @ "//routineN)
    ENDIF
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,1,MPI_2REAL,MPI_MAXLOC,gid, ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res,STAT=ierr)
#endif
  END SUBROUTINE mp_maxloc_rv




! *****************************************************************************
!> \brief Non-blocking send of custom type
!> \sa mp_isendrecv_iv
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

#if defined(__parallel)
    my_tag = 0
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
!> \sa mp_isendrecv_iv
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

#if defined(__parallel)
    my_tag = 0
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
!> \brief Returns the size of a data type in bytes
!> \par MPI mapping
!>      mpi_type_size 
!> 
!> \param[in] type_descriptor  data type                   
!> \param[out] type_size       size of the data type
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
!> \brief Opens a file 
!> \par MPI-I/O mapping  mpi_file_open 
!> \par STREAM-I/O mapping  OPEN
!> 
!> \param[in] groupid    message passing environment identifier
!> \param[in] filepath   path to the file 
!> \param[in] amode      access mode
!> \param[out] fh        file handle (file storage unit)
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

    CHARACTER(len=10)                        :: fstatus
    INTEGER                                  :: ierr, my_info

    ierr = 0
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
    fh = 1  ! UNIT=1
    OPEN(UNIT=fh, FILE=filepath, STATUS=fstatus, ACCESS="STREAM")
#endif
  END SUBROUTINE mp_file_open

! *****************************************************************************
!> \brief Closes a file 
!> \par MPI-I/O mapping   mpi_file_close 
!> \par STREAM-I/O mapping   CLOSE 
!>
!> \param[in] fh   file handle (file storage unit)                   
!> \par History
!>      11.2012 created [Hossein Bani-Hashemian]
! *****************************************************************************
  SUBROUTINE mp_file_close(fh)
    INTEGER, INTENT(IN)                      :: fh

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
!> \par MPI-I/O mapping   mpi_file_get_size
!> \par STREAM-I/O mapping   INQUIRE 
!>
!> \param[in] fh file handle (file storage unit)
!> \param[out] file_size  the file size                   
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
!> \par MPI-I/O mapping   mpi_file_write_at
!> \par STREAM-I/O mapping   WRITE   
!> \param[in] fh     file handle (file storage unit) 
!> \param[in] offset file offset (position)
!> \param[in] msg    data to be writen to the file
!> \param[in](optional) msglen number of the elements of data                    
! *****************************************************************************
  SUBROUTINE mp_file_write_at_chv(fh, offset, msg, msglen)
    CHARACTER, INTENT(IN)                      :: msg(:)
    INTEGER, INTENT(IN)                        :: fh 
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_chv', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
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

  SUBROUTINE mp_file_write_at_ch(fh, offset, msg)
    CHARACTER(LEN=*), INTENT(IN)               :: msg
    INTEGER, INTENT(IN)                        :: fh 
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_ch', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
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
!> \par MPI-I/O mapping   mpi_file_write_at_all 
!> \par STREAM-I/O mapping   WRITE   
! *****************************************************************************
  SUBROUTINE mp_file_write_at_all_chv(fh, offset, msg, msglen)
    CHARACTER, INTENT(IN)                      :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_all_chv', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
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

  SUBROUTINE mp_file_write_at_all_ch(fh, offset, msg)
    CHARACTER(LEN=*), INTENT(IN)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_all_ch', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
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
!> \par MPI-I/O mapping   mpi_file_read_at 
!> \par STREAM-I/O mapping   READ  
!> \param[in] fh     file handle (file storage unit)
!> \param[in] offset file offset (position)
!> \param[out] msg   data to be read from the file
!> \param[in](optional) msglen  number of elements of data                      
! *****************************************************************************
  SUBROUTINE mp_file_read_at_chv(fh, offset, msg, msglen)
    CHARACTER, INTENT(OUT)                     :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_chv', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
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

  SUBROUTINE mp_file_read_at_ch(fh, offset, msg)
    CHARACTER(LEN=*), INTENT(OUT)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_ch', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
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
!> \par MPI-I/O mapping    mpi_file_read_at_all 
!> \par STREAM-I/O mapping   READ  
! *****************************************************************************
  SUBROUTINE mp_file_read_at_all_chv(fh, offset, msg, msglen)
    CHARACTER, INTENT(OUT)                     :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_all_chv', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
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

  SUBROUTINE mp_file_read_at_all_ch(fh, offset, msg)
    CHARACTER(LEN=*), INTENT(OUT)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_all_ch', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT_ALL(fh, offset, msg, LEN(msg), MPI_CHARACTER, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_all_ch @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_all_ch

#include "dbcsr_message_passing_i.f90"
#include "dbcsr_message_passing_l.f90"
#include "dbcsr_message_passing_r.f90"
#include "dbcsr_message_passing_d.f90"
#include "dbcsr_message_passing_c.f90"
#include "dbcsr_message_passing_z.f90"

END MODULE dbcsr_message_passing
