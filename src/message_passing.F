!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/message_passing [1.0] *
!!
!!   NAME
!!     message_passing
!!
!!   FUNCTION
!!     Interface to the message passing library MPI
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (02-Jan-2001): New error handling
!!                        Performance tools
!!     JGH (14-Jan-2001): New routines mp_comm_compare, mp_cart_coords, 
!!                                     mp_rank_compare, mp_transpose
!!     JGH (06-Feb-2001): New routines mp_comm_free
!!     JGH (22-Mar-2001): New routines mp_comm_dup
!!     JVdV   ( 05-2002): mp_allgather,mp_range
!!
!!   SOURCE
!******************************************************************************

MODULE message_passing

  USE kinds, ONLY : dbl
  USE machine, ONLY : m_cputime
  USE parallel_include

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mp_start, mp_end, mp_world_setup, mp_group, mp_cart_create, &
       mp_bcast, mp_stop, mp_sum, mp_max, mp_min, mp_sync, mp_environ, &
       mp_gather, mp_comm_compare, mp_cart_coords, mp_rank_compare, &
       mp_transpose, mp_dims_create, mp_cart_rank, mp_cart_sub, mp_comm_free, &
       mp_comm_dup, mp_cart_shift, mp_sendrecv, mp_range, mp_allgather, &
       mp_sum_scatter, mp_abort, mp_shift

  INTERFACE mp_shift
     MODULE PROCEDURE mp_shift_rm
  END INTERFACE
  
  INTERFACE mp_bcast
     MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, mp_bcast_z, &
          mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
          mp_bcast_cm, mp_bcast_im, mp_bcast_rm3
  END INTERFACE

  INTERFACE mp_sum
     MODULE PROCEDURE mp_sum_i1, mp_sum_r1, mp_sum_c1, mp_sum_iv, &
          mp_sum_rv, mp_sum_cv, mp_sum_im, mp_sum_rm, mp_sum_cm, &
          mp_sum_im3, mp_sum_rm3, mp_sum_cm3, mp_sum_root_rv, &
          mp_sum_root_rm, mp_sum_root_rm3
  END INTERFACE

  INTERFACE mp_max
     MODULE PROCEDURE mp_max_i, mp_max_iv, mp_max_r, mp_max_rv
  END INTERFACE

  INTERFACE mp_min
     MODULE PROCEDURE mp_min_i, mp_min_iv, mp_min_r, mp_min_rv
  END INTERFACE

  INTERFACE mp_gather
     MODULE PROCEDURE mp_gather_i, mp_gather_iv, mp_gather_r, mp_gather_rv
  END INTERFACE

  INTERFACE mp_allgather
     MODULE PROCEDURE mp_allgather_rm3, mp_allgatherv_rv
  END INTERFACE

  INTERFACE mp_sum_scatter
     MODULE PROCEDURE mp_sum_scatter_rv
  END INTERFACE

  INTERFACE mp_environ
     MODULE PROCEDURE mp_environ_l, mp_environ_c
  END INTERFACE

  INTERFACE mp_transpose
     MODULE PROCEDURE mp_transpose_c22v, mp_transpose_i
  END INTERFACE

  INTERFACE mp_sendrecv
     MODULE PROCEDURE mp_sendrecv_rm3
  END INTERFACE

  TYPE mp_perf_type
    CHARACTER ( LEN = 20 ) :: name
    INTEGER :: count
    REAL ( dbl ) :: msg_size
    REAL ( dbl ) :: time
  END TYPE mp_perf_type

  INTEGER, PARAMETER :: MAX_PERF = 7
  TYPE ( mp_perf_type ), DIMENSION ( MAX_PERF ) :: mp_perf
  CHARACTER ( LEN = 20 ), PARAMETER :: sname ( MAX_PERF ) =  &
   (/"MP_Group            ", "MP_Bcast            ", "MP_Allreduce        ", &
     "MP_Gather           ", "MP_Sync             ", "MP_Transpose        ", &
     "MP_SendRecv         "/)
  REAL ( dbl ) :: t_start, t_end
  INTEGER :: intlen, reallen, charlen, loglen

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

!..mp_start
SUBROUTINE mp_start()

  IMPLICIT NONE

! Locals
  INTEGER :: ierr, i

  ierr = 0
#if defined(__parallel)
  CALL mpi_init ( ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_init @ mp_start" )
  CALL mpi_errhandler_set ( MPI_COMM_WORLD, mpi_errors_return, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_errhandler_set @ mp_start" )
  DO i = 1, MAX_PERF
    mp_perf (i) % name = sname ( i )
    mp_perf (i) % count = 0
    mp_perf (i) % msg_size = 0._dbl
    mp_perf (i) % time = 0._dbl
  END DO
  ! we make some assumptions on the length of INTEGERS, REALS and LOGICALS
  intlen  = BIT_SIZE ( 0 ) / 8
  reallen = 8
  loglen  = intlen
  charlen = 1
#endif

END SUBROUTINE mp_start

!******************************************************************************

!..mp_end
SUBROUTINE mp_end ( scr )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ), OPTIONAL :: scr
! Locals
  INTEGER :: ierr, i, iw, numtask, taskid, groupid
  REAL ( dbl ) :: vol, per

  ierr = 0
#if defined(__parallel)
  CALL mp_world_setup ( numtask, taskid, groupid )
  CALL mp_sync ( MPI_COMM_WORLD )
  CALL mpi_finalize ( ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_finalize @ mp_end" )
  IF ( PRESENT ( scr ) ) THEN
    iw = scr
  ELSE
    iw = 6
  END IF
  IF ( taskid == 0 ) THEN
    WRITE ( iw, '( /, 1X, 79("-") )' )
    WRITE ( iw, '( " -", 77X, "-" )' )
    WRITE ( iw, '( " -", 24X, A, 24X, "-" )' ) ' MESSAGE PASSING PERFORMANCE '
    WRITE ( iw, '( " -", 77X, "-" )' )
    WRITE ( iw, '( 1X, 79("-"), / )' )
    WRITE ( iw, '( A, A, A )' ) ' ROUTINE', '             CALLS ', &
            ' TOT TIME [s]  AVE VOLUME [Bytes]  PERFORMANCE [MB/s]'
    DO i = 1, MAX_PERF

      IF ( mp_perf ( i ) % count > 0 ) THEN
        vol = mp_perf ( i ) % msg_size / REAL ( mp_perf ( i ) % count, dbl )
        IF ( mp_perf ( i ) % time > 0._dbl ) THEN
          per = mp_perf ( i ) % msg_size / mp_perf (i) % time * 1.e-6_dbl
        ELSE
          per = 0._dbl
        ENDIF
        IF ( vol < 1._dbl ) THEN
          WRITE ( iw, '(1X,A15,T17,I10,T27,F14.3)' ) &
            ADJUSTL ( mp_perf ( i ) % name ), mp_perf ( i ) % count, &
            mp_perf ( i ) % time
        ELSE
          WRITE ( iw, '(1X,A15,T17,I10,T27,F14.3,T50,F11.0,T69,F12.2)' ) &
            ADJUSTL ( mp_perf ( i ) % name ), mp_perf ( i ) % count, &
            mp_perf ( i ) % time, vol, per
        END IF
      ENDIF

    END DO
    WRITE ( iw, '( 1X, 79("-"), / )' )
  END IF
#endif

END SUBROUTINE mp_end

!******************************************************************************

!..mp_stop
SUBROUTINE mp_stop ( ierr, prg_code )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: ierr
  CHARACTER ( LEN=* ) :: prg_code

! Locals
  INTEGER :: numtask, taskid, groupid, len
#if defined(__parallel)
  CHARACTER ( LEN = MPI_MAX_ERROR_STRING ) :: error_string
#endif

  CALL mp_world_setup ( numtask, taskid, groupid )

  WRITE(*,'(/,A,T71,I10)') ' CP2K| Stopped by processor number',taskid
  WRITE(*,'(A,T31,A50)') ' CP2K| ', ADJUSTR ( prg_code )
  WRITE(*,'(A,T71,I10,/)') ' CP2K| Error number was ',ierr

#if defined(__parallel)
  CALL mpi_error_string ( ierr, error_string, len )
  WRITE(*,'(A,A)') ' CP2K| ', ADJUSTR ( error_string )
  CALL mpi_abort ( MPI_COMM_WORLD, prg_code )
#endif

  STOP "mp_stop"


END SUBROUTINE mp_stop

!******************************************************************************

!..mp_abort
SUBROUTINE mp_abort ( )

! Locals
  INTEGER :: numtask, taskid, groupid

  CALL mp_world_setup ( numtask, taskid, groupid )

  WRITE(*,'(/,A,T71,I10)') ' CP2K| Stopped by process number',taskid
  WRITE(*,'(A,/)')         ' CP2K| Abnormal program termination'

#if defined(__parallel)
  CALL mpi_abort ( MPI_COMM_WORLD, 0 )
#else
  STOP
#endif

END SUBROUTINE mp_abort

!******************************************************************************

!..mp_sync
SUBROUTINE mp_sync ( group )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: group

! Locals
  INTEGER :: ierr

#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_barrier ( group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_barrier @ mp_sync" )
  mp_perf ( 5 ) % count = mp_perf ( 5 ) % count + 1
  t_end = m_cputime ( )
  mp_perf ( 5 ) % time = mp_perf ( 5 ) % time + ( t_end - t_start )
#endif

END SUBROUTINE mp_sync

!******************************************************************************

!..mp_world_setup
SUBROUTINE mp_world_setup ( numtask, taskid, groupid )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( OUT ) :: numtask, taskid, groupid

! Locals
  INTEGER :: ierr

  ierr = 0
  numtask = 1
  taskid = 0
  groupid = 0
#if defined(__parallel)
  CALL mpi_comm_rank ( MPI_COMM_WORLD, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_rank @ mp_world_setup" )

  CALL mpi_comm_size ( MPI_COMM_WORLD, numtask, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_world_setup" )
  groupid = MPI_COMM_WORLD
#endif

END SUBROUTINE mp_world_setup

!******************************************************************************

!..mp_environ
SUBROUTINE mp_environ_l ( numtask, taskid, groupid )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: groupid
  INTEGER, INTENT ( OUT ) :: numtask, taskid

! Locals
  INTEGER :: ierr

  ierr = 0
  numtask = 1
  taskid = 0
#if defined(__parallel)
  CALL mpi_comm_rank ( groupid, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_rank @ mp_environ_l" )

  CALL mpi_comm_size ( groupid, numtask, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_environ_l" )
#endif

END SUBROUTINE mp_environ_l

SUBROUTINE mp_environ_c ( numtask, dims, task_coor, groupid )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: groupid
  INTEGER, INTENT ( OUT ) :: numtask, task_coor ( 2 ), dims ( 2 )

! Locals
  INTEGER :: ierr, rank
  LOGICAL :: periods ( 2 )

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

!******************************************************************************

!..mp_group
SUBROUTINE mp_group(group_list,group_size,base_group,groupid)

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: group_list ( : ), group_size, base_group
  INTEGER, INTENT ( OUT ) :: groupid

! Locals
  INTEGER :: base, newgroup, ierr

  ierr = 0
  groupid = base_group
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_comm_group ( base_group, base, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_group @ mp_group" )

  CALL mpi_group_incl ( base, group_size, group_list, newgroup, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_group_incl @ mp_group" )

  CALL mpi_comm_create ( base_group, newgroup, groupid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_create @ mp_group" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  t_end = m_cputime ( )
  mp_perf ( 1 ) % time = mp_perf ( 1 ) % time + ( t_end - t_start )
#endif

END SUBROUTINE mp_group

!******************************************************************************

!..mp_cart_create
SUBROUTINE mp_cart_create ( comm_old, ndims, dims, pos, comm_cart )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: comm_old, ndims
  INTEGER, INTENT ( INOUT ) :: dims ( : )
  INTEGER, INTENT ( OUT ) :: pos ( : ), comm_cart

! Locals
  INTEGER :: ierr, nodes
  LOGICAL :: period(1:ndims), reorder

  ierr = 0
  pos ( 1:ndims ) = 1
  comm_cart = comm_old
#if defined(__parallel)

  t_start = m_cputime ( )
  CALL mpi_comm_size ( comm_old, nodes, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_cart_create" )

  CALL mpi_dims_create ( nodes, ndims, dims, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_dims_create @ mp_cart_create" )

  reorder = .TRUE.
  period = .TRUE.
  CALL mpi_cart_create ( comm_old, ndims, dims, period, reorder, comm_cart, &
       ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_create @ mp_cart_create" )

  CALL mpi_cart_get ( comm_cart, ndims, dims, period, pos, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_get @ mp_cart_create" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  t_end = m_cputime ( )
  mp_perf ( 1 ) % time = mp_perf ( 1 ) % time + ( t_end - t_start )
#else
  dims = 1
  comm_cart = 0
#endif
  
END SUBROUTINE mp_cart_create

!******************************************************************************

!..mp_cart_coords 
SUBROUTINE mp_cart_coords ( comm, rank, coords)

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: comm, rank
  INTEGER, DIMENSION ( : ), INTENT ( OUT ) :: coords

! Locals
  INTEGER :: ierr, m

  ierr = 0
  m = SIZE ( coords )
#if defined(__parallel)
  CALL mpi_cart_coords ( comm, rank, m, coords, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_coords @ mp_cart_coords" )
#else
  coords = 0
#endif

END SUBROUTINE mp_cart_coords

!******************************************************************************

!..mp_cart_shift
SUBROUTINE mp_cart_shift ( comm, dir, disp, source, dest )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: comm, dir, disp
  INTEGER, INTENT ( OUT ) :: source, dest

! Locals
  INTEGER :: ierr

  ierr = 0
#if defined(__parallel)
  CALL mpi_cart_shift ( comm, dir, disp, source, dest, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_shift @ mp_cart_shift" )
#else
  source = 0
  dest = 0
#endif

END SUBROUTINE mp_cart_shift

!******************************************************************************

!..mp_comm_compare
SUBROUTINE mp_comm_compare ( comm1, comm2, result)

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: comm1, comm2
  INTEGER, INTENT ( OUT ) :: result

! Locals
  INTEGER :: ierr, iout

  ierr = 0
  result = 0
#if defined(__parallel)
  CALL mpi_comm_compare ( comm1, comm2, iout, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_compare @ mp_comm_compare" )
  IF ( iout == MPI_UNEQUAL ) result = 1
#endif

END SUBROUTINE mp_comm_compare

!******************************************************************************

!..mp_cart_sub
SUBROUTINE mp_cart_sub ( comm, rdim, sub_comm )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: comm
  LOGICAL, DIMENSION ( : ), INTENT ( IN ) :: rdim
  INTEGER, INTENT ( OUT ) :: sub_comm

! Locals
  INTEGER :: ierr

  ierr = 0
#if defined(__parallel)
  CALL mpi_cart_sub ( comm, rdim, sub_comm, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_cart_sub @ mp_cart_sub" )
#endif

END SUBROUTINE mp_cart_sub

!******************************************************************************

!..mp_comm_free
SUBROUTINE mp_comm_free ( comm )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: comm

! Locals
  INTEGER :: ierr

  ierr = 0
#if defined(__parallel)
  CALL mpi_comm_free ( comm, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_free @ mp_comm_free" )
#endif

END SUBROUTINE mp_comm_free

!******************************************************************************

!..mp_comm_dup
SUBROUTINE mp_comm_dup ( comm1, comm2 )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: comm1
  INTEGER, INTENT ( OUT ) :: comm2

! Locals
  INTEGER :: ierr

  ierr = 0
#if defined(__parallel)
  CALL mpi_comm_dup ( comm1, comm2, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_dup @ mp_comm_dup" )
#endif

END SUBROUTINE mp_comm_dup

!******************************************************************************

!..mp_rank_compare
SUBROUTINE mp_rank_compare ( comm1, comm2, rank )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: comm1, comm2
  INTEGER, DIMENSION ( : ), INTENT ( OUT ) :: rank

! Locals
  INTEGER :: ierr, n1, n2, n, g1, g2, i
  INTEGER, DIMENSION ( : ), ALLOCATABLE :: rin

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
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_rank_compare" )
  DO i = 0, n-1
    rin ( i ) = i
  END DO
  CALL mpi_group_translate_ranks ( g1, n, rin, g2, rank, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, &
      "mpi_group_translate_rank @ mp_rank_compare" )
  DEALLOCATE ( rin, STAT = ierr )
  IF ( ierr /= 0 ) CALL mp_stop( 0, "deallocate @ mp_rank_compare" )
#endif

END SUBROUTINE mp_rank_compare

!******************************************************************************

!..mp_dims_create
SUBROUTINE mp_dims_create ( nodes, dims )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: nodes
  INTEGER, DIMENSION ( : ), INTENT ( INOUT ) :: dims

! Locals
  INTEGER :: ndim, ierr

#if defined(__parallel)
  ndim = SIZE ( dims )
  CALL mpi_dims_create ( nodes, ndim, dims ,ierr)
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_dims_create @ mp_dims_create" )
#else
  dims = 1
#endif

END SUBROUTINE mp_dims_create

!******************************************************************************

!******************************************************************************
! this routine will shift around the data in msg
! i.e. msg will be moved from rank to rank+displ_in (in a circular way)
! displ_in will be 1 by default (others not tested)
! the message array needs to be the same size on all processes
!*****************************************************************************
SUBROUTINE mp_shift_rm( msg, group, displ_in)

  IMPLICIT NONE
  INTEGER, INTENT ( IN )               :: group
  REAL ( dbl ) , INTENT(INOUT)         :: msg ( :, : )
  INTEGER, INTENT ( IN ), OPTIONAL     :: displ_in
! 
  integer, dimension(:), allocatable ::  status
  INTEGER displ,myrank,nprocs,left,right,ierror,msglen,tag

#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  call mpi_comm_rank(group,myrank,ierror)
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_rank @ mp_shift" )
  call mpi_comm_size(group,nprocs,ierror)
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_size @ mp_shift" )
  IF (PRESENT(displ_in)) THEN
     displ=displ_in 
  ELSE
     displ=1
  ENDIF
  right=MOD(nprocs+myrank+displ,nprocs) 
  left =MOD(nprocs+myrank-displ,nprocs) 
  tag=17
  msglen = size(msg)
  t_start = m_cputime ( )
  CALL mpi_sendrecv_replace(msg,msglen,MPI_DOUBLE_PRECISION,right,tag,left,tag,group,status(1),ierror)
  t_end = m_cputime ( )
  IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_sendrecv_replace @ mp_shift" )
  mp_perf ( 7 ) % count = mp_perf ( 7 ) % count + 1
  mp_perf ( 7 ) % msg_size = mp_perf ( 7 ) % msg_size + msglen * reallen
  mp_perf ( 7 ) % time = mp_perf ( 7 ) % time + ( t_end - t_start )
  DEALLOCATE(status)
#endif

END SUBROUTINE mp_shift_rm

!******************************************************************************

!..mp_cart_rank
SUBROUTINE mp_cart_rank ( group, pos, rank )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ) :: group
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: pos
  INTEGER, INTENT ( OUT ) :: rank

! Locals
  INTEGER :: ierr

#if defined(__parallel)
  CALL mpi_cart_rank ( group, pos, rank, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_cart_rank @ mp_cart_rank" )
#else
  rank = 1
#endif

END SUBROUTINE mp_cart_rank
!******************************************************************************

!..mp_transpose
SUBROUTINE mp_transpose_c22v ( sb, scount, sdispl, rb, rcount, rdispl, group )

  IMPLICIT NONE

! Arguments
  COMPLEX ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: sb
  COMPLEX ( dbl ), DIMENSION ( :, : ), INTENT ( OUT ) :: rb
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: scount, sdispl,rcount, rdispl
  INTEGER, INTENT ( IN ) :: group

! Locals
  INTEGER :: ierr, msglen

  ierr = 0
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_alltoallv ( sb, scount, sdispl, MPI_DOUBLE_COMPLEX, &
                       rb, rcount, rdispl, MPI_DOUBLE_COMPLEX, group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ mp_transpose_c22v" )
  mp_perf ( 6 ) % count = mp_perf ( 6 ) % count + 1
  msglen = SUM ( scount ) + SUM ( rcount )
  mp_perf ( 6 ) % msg_size = mp_perf ( 6 ) % msg_size + msglen * 2 * reallen
  t_end = m_cputime ( )
  mp_perf ( 6 ) % time = mp_perf ( 6 ) % time + ( t_end - t_start )
#endif

END SUBROUTINE mp_transpose_c22v

SUBROUTINE mp_transpose_i ( sb, rb, count, group )

  IMPLICIT NONE

! Arguments
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: sb
  INTEGER, DIMENSION ( : ), INTENT ( OUT ) :: rb
  INTEGER, INTENT ( IN ) :: count, group

! Locals
  INTEGER :: ierr, msglen, np

  ierr = 0
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_alltoall ( sb, count, MPI_INTEGER, &
                      rb, count, MPI_INTEGER, group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ mp_transpose_i" )
  mp_perf ( 6 ) % count = mp_perf ( 6 ) % count + 1
  CALL mpi_comm_size ( group, np, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ mp_transpose_i" )
  msglen = 2 * count * np
  mp_perf ( 6 ) % msg_size = mp_perf ( 6 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 6 ) % time = mp_perf ( 6 ) % time + ( t_end - t_start )
#endif

END SUBROUTINE mp_transpose_i

!******************************************************************************

!..mp_bcast
SUBROUTINE mp_bcast_i1(msg,source,gid)
  IMPLICIT NONE
  INTEGER :: msg
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_INTEGER,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_i1" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_i1
SUBROUTINE mp_bcast_iv(msg,source,gid)
  IMPLICIT NONE
  INTEGER :: msg ( : )
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = size(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_INTEGER,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_iv" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_iv
SUBROUTINE mp_bcast_im(msg,source,gid)
  IMPLICIT NONE
  INTEGER :: msg ( :, : )
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = size(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_INTEGER,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_im" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_im
SUBROUTINE mp_bcast_r1(msg,source,gid)
  IMPLICIT NONE
  REAL ( dbl ) :: msg
  INTEGER :: msglen, source, gid
  INTEGER :: ierr

  ierr = 0
  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_PRECISION,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_r1" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_r1
SUBROUTINE mp_bcast_rv(msg,source,gid)
  IMPLICIT NONE
  REAL ( dbl ) :: msg ( : )
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = size(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_PRECISION,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_rv" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_rv
SUBROUTINE mp_bcast_rm(msg,source,gid)
  IMPLICIT NONE
  REAL ( dbl ) :: msg ( :, : )
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = size(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_PRECISION,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_rm" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_rm
SUBROUTINE mp_bcast_rm3(msg,source,gid)
  IMPLICIT NONE
  REAL ( dbl ) :: msg ( :, :, : )
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = size(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_PRECISION,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_rm" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_rm3
SUBROUTINE mp_bcast_c1(msg,source,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ) :: msg
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_c1" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * reallen*2
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_c1
SUBROUTINE mp_bcast_cv(msg,source,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ) :: msg ( : )
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = size(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_cv" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * reallen*2
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_cv
SUBROUTINE mp_bcast_cm(msg,source,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ) :: msg ( :, : )
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = size(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_DOUBLE_COMPLEX,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_cm")
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * reallen*2
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_cm
SUBROUTINE mp_bcast_l(msg,source,gid)
  IMPLICIT NONE
  LOGICAL :: msg
  INTEGER :: source, gid
  INTEGER :: msglen, ierr

  ierr = 0
  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_bcast(msg,msglen,MPI_LOGICAL,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_l" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * loglen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_l
SUBROUTINE mp_bcast_z(msg,source,gid)
  IMPLICIT NONE
  CHARACTER ( LEN = * ) :: msg
  INTEGER :: source, gid
  INTEGER :: msglen, ierr, i
  INTEGER, ALLOCATABLE :: imsg ( : )

  ierr = 0
  msglen = len(msg)
#if defined(__parallel)
  t_start = m_cputime ( )
! this is a workaround to avoid problems on the T3E
! at the moment we have a data alignment error when trying to
! broadcats characters on the T3E (not always!)
! JH 19/3/99 on galileo
! CALL mpi_bcast(msg,msglen,MPI_CHARACTER,source,gid,ierr)
  ALLOCATE (imsg(1:msglen))
  DO i = 1, msglen
     imsg(i) = ichar(msg(i:i))
  END DO
  CALL mpi_bcast(imsg,msglen,MPI_INTEGER,source,gid,ierr)
  DO i = 1, msglen
     msg(i:i) = char(imsg(i))
  END DO
  DEALLOCATE (imsg)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ mp_bcast_z" )
  mp_perf ( 2 ) % count = mp_perf ( 2 ) % count + 1
  mp_perf ( 2 ) % msg_size = mp_perf ( 2 ) % msg_size + msglen * charlen
  t_end = m_cputime ( )
  mp_perf ( 2 ) % time = mp_perf ( 2 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_bcast_z

!******************************************************************************

!..mp_sum
SUBROUTINE mp_sum_i1(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, res, ierr

  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_SUM,gid,ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_i1" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_i1
SUBROUTINE mp_sum_iv(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER, ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_iv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_SUM,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_iv" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_iv
SUBROUTINE mp_sum_im(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, ierr
  INTEGER, ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  m1 = size(msg,1)
  m2 = size(msg,2)
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_im" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_SUM,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_im" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_im
SUBROUTINE mp_sum_im3(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( :, :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, m3, ierr
  INTEGER, ALLOCATABLE :: res ( :, :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  m1 = size(msg,1)
  m2 = size(msg,2)
  m3 = size(msg,3)
  ALLOCATE (res(m1,m2,m3),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_im3" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_SUM,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_im3" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_im3
SUBROUTINE mp_sum_r1(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ) :: res

  msglen = 1
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_r1" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_r1
SUBROUTINE mp_sum_rv(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_rv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_rv" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_rv
SUBROUTINE mp_sum_root_rv(msg,root,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: root,gid
  INTEGER :: msglen, m1, ierr, taskid
  REAL ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_comm_rank ( gid, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ mp_sum_root_rm" )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  ALLOCATE (res(m1),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_root_rv" )
  CALL mpi_reduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,&
       root,gid,ierr)
  IF ( taskid == root ) THEN
    msg = res
  END IF
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ mp_sum_root_rv" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_root_rv
SUBROUTINE mp_sum_rm(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  m1 = size(msg,1)
  m2 = size(msg,2)
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_rm" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_rm" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_rm
SUBROUTINE mp_sum_root_rm(msg,root,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: root,gid
  INTEGER :: msglen, m1, m2, ierr, taskid
  REAL ( dbl ), ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_comm_rank ( gid, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ mp_sum_root_rm" )
  msglen = SIZE(msg)
  m1 = SIZE(msg,1)
  m2 = SIZE(msg,2)
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_root_rm" )
  CALL mpi_reduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,&
       root,gid,ierr)
  IF ( taskid == root ) THEN
    msg = res
  END IF
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ mp_sum_root_rm" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_root_rm
SUBROUTINE mp_sum_rm3(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( :, :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, m3, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( :, :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  m1 = size(msg,1)
  m2 = size(msg,2)
  m3 = size(msg,3)
  ALLOCATE (res(m1,m2,m3),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_rm3" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_rm3" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_rm3
SUBROUTINE mp_sum_root_rm3(msg,root,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( :, :, : )
  INTEGER, INTENT ( IN ) :: root,gid
  INTEGER :: msglen, m1, m2, m3, ierr, taskid
  REAL ( dbl ), ALLOCATABLE :: res ( :, :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_comm_rank ( gid, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ mp_sum_root_rm3" )
  msglen = size(msg)
  m1 = size(msg,1)
  m2 = size(msg,2)
  m3 = size(msg,3)
  ALLOCATE (res(m1,m2,m3),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_root_rm3" )
  CALL mpi_reduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,&
       root,gid,ierr)
  IF ( taskid == root ) THEN
    msg = res
  END IF
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ mp_sum_root_rm3" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_root_rm3
SUBROUTINE mp_sum_c1(msg,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  COMPLEX ( dbl ) :: res

  msglen = 2
#if defined(__parallel)
  t_start = m_cputime ( )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_c1" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen*2
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_c1
SUBROUTINE mp_sum_cv(msg,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  COMPLEX ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_cv" )
  CALL mpi_allreduce(msg,res,2*msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_cv" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen*2
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_cv
SUBROUTINE mp_sum_cm(msg,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ), INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, ierr
  COMPLEX ( dbl ), ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 2*size(msg)
  m1 = size(msg,1)
  m2 = size(msg,2)
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_cm" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_cm" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen*2
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_cm
SUBROUTINE mp_sum_cm3(msg,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ), INTENT ( INOUT ) :: msg ( :, :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, m3, ierr
  COMPLEX ( dbl ), ALLOCATABLE :: res ( :, :, : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 2*size(msg)
  m1 = size(msg,1)
  m2 = size(msg,2)
  m3 = size(msg,3)
  ALLOCATE (res(m1,m2,m3),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_sum_cm3" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_SUM,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_sum_cm3" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen*2
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_sum_cm3

!******************************************************************************

SUBROUTINE mp_max_i(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER :: res
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_MAX,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_max_i" )
  msg = res
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_max_i
SUBROUTINE mp_max_iv(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER, ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_max_iv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_MAX,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_max_iv" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_max_iv
SUBROUTINE mp_max_r(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ) :: res
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_MAX,gid, &
       ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_max_r" )
  msg = res
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_max_r
SUBROUTINE mp_max_rv(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_max_rv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_MAX,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_max_rv" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_max_rv

!******************************************************************************

SUBROUTINE mp_min_i(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER :: res
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_MIN,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_i" )
  msg = res
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_min_i
SUBROUTINE mp_min_iv(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER, ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_min_iv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_INTEGER,MPI_MIN,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_iv" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_min_iv
SUBROUTINE mp_min_r(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ) :: res
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_MIN,gid, &
       ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_r" )
  msg = res
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_min_r
SUBROUTINE mp_min_rv(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop( 0, "allocate @ mp_min_rv" )
  CALL mpi_allreduce(msg,res,msglen,MPI_DOUBLE_PRECISION,MPI_MIN,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ mp_min_rv" )
  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#endif
END SUBROUTINE mp_min_rv

!******************************************************************************

SUBROUTINE mp_gather_i(msg,msg_gather,root,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: msg
  INTEGER, INTENT ( OUT ) :: msg_gather ( : )
  INTEGER, INTENT ( IN ) :: root
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_gather(msg,msglen,MPI_INTEGER,msg_gather,&
                  msglen,MPI_INTEGER,root,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ mp_gather_i" )
  mp_perf ( 4 ) % count = mp_perf ( 4 ) % count + 1
  mp_perf ( 4 ) % msg_size = mp_perf ( 4 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 4 ) % time = mp_perf ( 4 ) % time + ( t_end - t_start )
#else
  msg_gather = msg
#endif
END SUBROUTINE mp_gather_i

SUBROUTINE mp_gather_iv(msg,msg_gather,root,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: msg ( : )
  INTEGER, INTENT ( OUT ) :: msg_gather ( : )
  INTEGER, INTENT ( IN ) :: root
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  CALL mpi_gather(msg,msglen,MPI_INTEGER,msg_gather,&
                  msglen,MPI_INTEGER,root,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ mp_gather_iv" )
  mp_perf ( 4 ) % count = mp_perf ( 4 ) % count + 1
  mp_perf ( 4 ) % msg_size = mp_perf ( 4 ) % msg_size + msglen * intlen
  t_end = m_cputime ( )
  mp_perf ( 4 ) % time = mp_perf ( 4 ) % time + ( t_end - t_start )
#else
  msg_gather = msg
#endif
END SUBROUTINE mp_gather_iv

SUBROUTINE mp_gather_r(msg,msg_gather,root,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( IN ) :: msg
  REAL ( dbl ), INTENT ( OUT ) :: msg_gather ( : )
  INTEGER, INTENT ( IN ) :: root
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = 1
  CALL mpi_gather(msg,msglen,MPI_DOUBLE_PRECISION,msg_gather,&
                  msglen,MPI_DOUBLE_PRECISION,root,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ mp_gather_r" )
  mp_perf ( 4 ) % count = mp_perf ( 4 ) % count + 1
  mp_perf ( 4 ) % msg_size = mp_perf ( 4 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 4 ) % time = mp_perf ( 4 ) % time + ( t_end - t_start )
#else
  msg_gather = msg
#endif
END SUBROUTINE mp_gather_r

SUBROUTINE mp_gather_rv(msg,msg_gather,root,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( IN ) :: msg ( : )
  REAL ( dbl ), INTENT ( OUT ) :: msg_gather ( : )
  INTEGER, INTENT ( IN ) :: root
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
#if defined(__parallel)
  t_start = m_cputime ( )
  msglen = size(msg)
  CALL mpi_gather(msg,msglen,MPI_DOUBLE_PRECISION,msg_gather,&
                  msglen,MPI_DOUBLE_PRECISION,root,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ mp_gather_rv" )
  mp_perf ( 4 ) % count = mp_perf ( 4 ) % count + 1
  mp_perf ( 4 ) % msg_size = mp_perf ( 4 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 4 ) % time = mp_perf ( 4 ) % time + ( t_end - t_start )
#else
  msg_gather = msg
#endif
END SUBROUTINE mp_gather_rv

!**********************************************************************************
! JVdV computes +- equal ranges within 1..elements 
! each pe is supposed to calculate local_offset+1 .. local_offset+local_planes
! so that allgather can do the job of broadcast the data 
! the optional arguments describe the full destribution
!**********************************************************************************

SUBROUTINE mp_range(mepos,num_pe,elements,local_offset,local_planes,offsets,planes)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: elements
  INTEGER, INTENT ( IN ) :: mepos
  INTEGER, INTENT ( IN ) :: num_pe
  INTEGER, INTENT ( OUT ) :: local_offset
  INTEGER, INTENT ( OUT ) :: local_planes
  INTEGER, INTENT ( OUT ), OPTIONAL :: offsets ( : )
  INTEGER, INTENT ( OUT ), OPTIONAL :: planes ( : )

  INTEGER :: equal_planes,i

  if (MODULO(elements,num_pe).eq.0) then
       equal_planes=elements/num_pe
  else
       equal_planes=1+elements/num_pe
  endif
  do i=1,num_pe
       IF (present(offsets)) THEN
          offsets(I)=(I-1)*equal_planes
       END IF
       IF (present(planes)) THEN
          planes(I)=MAX(0,MIN(elements-(I-1)*equal_planes,equal_planes))
       ENDIF
  enddo
  local_offset=mepos*equal_planes
  local_planes=MAX(0,MIN(elements-local_offset,equal_planes))
END SUBROUTINE mp_range

!********************************************************************************
!
! gathers all the data in msg.
! that is, if each PE computes msg( : , : , local_offset+1 .. local_offset+local+planes)
! mp_allgather can be used to distribute the data so that every PE has the full msg
! This is significantly faster then summing the data
!
!********************************************************************************

SUBROUTINE mp_allgather_rm3(mepos,num_pe,elements,msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( : , : , : )
  INTEGER, INTENT ( IN ) :: gid,num_pe,mepos,elements

  INTEGER :: ierr,plane_size,local_offset,local_planes
  INTEGER, POINTER, DIMENSION(:) :: offsets,planes

#if defined(__parallel)
  allocate(offsets(num_pe),planes(num_pe))
  call mp_range(mepos,num_pe,elements,local_offset,local_planes, &
                  offsets=offsets,planes=planes)

  plane_size=size(msg,1)*size(msg,2)
  
  offsets(:)=offsets(:)*plane_size
  planes(:)=planes(:)*plane_size
  CALL MPI_ALLGATHERV(msg(1,1,local_offset+1),local_planes*plane_size, &
                        MPI_DOUBLE_PRECISION,msg(1,1,1),planes,offsets,&
                        MPI_DOUBLE_PRECISION,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ mp_allgather_rm3" )

  deallocate(offsets,planes)
#else
  ierr = 0
#endif
END SUBROUTINE mp_allgather_rm3

! wrapper for MPI_ALLGATHERV subroutine
! one dimensional real arrays
SUBROUTINE mp_allgatherv_rv(msgout,msgin,rcount,rdispl,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( IN ) :: msgout ( : )
  REAL ( dbl ), INTENT ( OUT ) :: msgin ( : )
  INTEGER, INTENT ( IN ) :: rcount ( : ), rdispl ( : )
  INTEGER, INTENT ( IN ) :: gid

  INTEGER :: ierr, scount

#if defined(__parallel)
  scount = SIZE ( msgout )
  CALL MPI_ALLGATHERV(msgout, scount, MPI_DOUBLE_PRECISION, msgin, rcount, &
                      rdispl, MPI_DOUBLE_PRECISION, gid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgatherv @ mp_allgatherv_rv" )
#else
  msgin = msgout
  ierr = 0
#endif
END SUBROUTINE mp_allgatherv_rv

!******************************************************************************
! wrapper for MPI_REDUCE_SCATTER subroutine
! one dimensional real arrays
SUBROUTINE mp_sum_scatter_rv(msgout,msgin,rcount,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( IN ) :: msgout ( : )
  REAL ( dbl ), INTENT ( OUT ) :: msgin ( : )
  INTEGER, INTENT ( IN ) :: rcount ( : )
  INTEGER, INTENT ( IN ) :: gid

  INTEGER :: ierr

#if defined(__parallel)
  t_start = m_cputime ( )
  CALL MPI_REDUCE_SCATTER(msgout, msgin, rcount, MPI_DOUBLE_PRECISION, MPI_SUM, &
                          gid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce_scatter @ mp_sum_scatter_rv" )

  mp_perf ( 3 ) % count = mp_perf ( 3 ) % count + 1
  mp_perf ( 3 ) % msg_size = mp_perf ( 3 ) % msg_size + rcount(1) * reallen*2
  t_end = m_cputime ( )
  mp_perf ( 3 ) % time = mp_perf ( 3 ) % time + ( t_end - t_start )
#else
  msgin = msgout
  ierr = 0
#endif
END SUBROUTINE mp_sum_scatter_rv

!******************************************************************************

SUBROUTINE mp_sendrecv_rm3(msgin,dest,msgout,source,comm)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( IN ) :: msgin ( :, :, : )
  REAL ( dbl ), INTENT ( OUT ) :: msgout ( :, :, : )
  INTEGER, INTENT ( IN ) :: dest, source, comm
  INTEGER :: msglen, ierr, send_tag, recv_tag
  integer, dimension(:), allocatable ::  status
#if defined(__parallel)
  ALLOCATE(status(MPI_STATUS_SIZE))
  t_start = m_cputime ( )
  msglen = SIZE(msgin)
  send_tag = 0 ! cannot think of something better here, this might be dangerous
  recv_tag = 0 ! cannot think of something better here, this might be dangerous
  CALL mpi_sendrecv(msgin,msglen,MPI_DOUBLE_PRECISION,dest,send_tag,msgout,&
       msglen,MPI_DOUBLE_PRECISION,source,recv_tag,comm,status(1),ierr)
  ! we do not check the status
  IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ mp_send_recv_rm3" )
  mp_perf ( 7 ) % count = mp_perf ( 7 ) % count + 1
  mp_perf ( 7 ) % msg_size = mp_perf ( 7 ) % msg_size + msglen * reallen
  t_end = m_cputime ( )
  mp_perf ( 7 ) % time = mp_perf ( 7 ) % time + ( t_end - t_start )
  DEALLOCATE(status)
#else
  msgout = msgin
#endif
END SUBROUTINE mp_sendrecv_rm3

!******************************************************************************

END MODULE message_passing

!******************************************************************************
