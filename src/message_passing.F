!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/message_passing [1.0] *
!!
!!   NAME
!!     message_passing
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE message_passing
  
  USE kinds, ONLY : dbl
  USE parallel_include
  
  PRIVATE
  PUBLIC :: mp_start, mp_end, mp_environ, mp_group, mp_cart_create, &
       mp_bcast, mp_stop, mp_sum, mp_max, mp_min, mp_sync
  
  INTERFACE mp_bcast
     MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, mp_bcast_z, &
          mp_bcast_iv, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
          mp_bcast_cm, mp_bcast_im
  END INTERFACE
  
  INTERFACE mp_sum
     MODULE PROCEDURE mp_sum_i1, mp_sum_r1, mp_sum_c1, mp_sum_iv, &
          mp_sum_rv, mp_sum_cv, mp_sum_rm, mp_sum_cm
  END INTERFACE
  
  INTERFACE mp_max
     MODULE PROCEDURE mp_max_i, mp_max_r, mp_max_rv
  END INTERFACE
  
  INTERFACE mp_min
     MODULE PROCEDURE mp_min_i, mp_min_r, mp_min_rv
  END INTERFACE
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
  
!..mp_start
SUBROUTINE mp_start()
  
  IMPLICIT NONE
  
! Locals
  INTEGER :: ierr
  
  ierr = 0
#if defined(__parallel)
  CALL mpi_init ( ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8000 )
#endif
  
END SUBROUTINE mp_start

!******************************************************************************

!..mp_end
SUBROUTINE mp_end()
  
  IMPLICIT NONE
  
! Locals
  INTEGER :: ierr
  
  ierr = 0
#if defined(__parallel)
  CALL mpi_finalize ( ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8999 )
#endif
  
END SUBROUTINE mp_end

!******************************************************************************

!..mp_stop
SUBROUTINE mp_stop ( code )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: code
  
! Locals
  INTEGER :: numtask, taskid, groupid
  
#if defined(__parallel)
  CALL mpi_abort ( mpi_comm_world, code )
#endif
  
  CALL mp_environ ( numtask, taskid, groupid )
 
  WRITE(*,'(/,A,T71,I10,/)') ' CP2K| Stopped by processor number',taskid

  STOP "mp_stop"

  
END SUBROUTINE mp_stop

!******************************************************************************

!..mp_sync
SUBROUTINE mp_sync ( group )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: group
  
! Locals
  INTEGER :: ierr
  
#if defined(__parallel)
  CALL mpi_barrier ( group, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8400 )
#endif
  
END SUBROUTINE mp_sync

!******************************************************************************

!..mp_environ
SUBROUTINE mp_environ ( numtask, taskid, groupid )
  
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
  CALL mpi_comm_rank ( mpi_comm_world, taskid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8001 )
  
  CALL mpi_comm_size ( mpi_comm_world, numtask, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8002 )
  groupid = mpi_comm_world
#endif
  
END SUBROUTINE mp_environ

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
  CALL mpi_comm_group ( base_group, base, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8010 )
  
  CALL mpi_group_incl ( base, group_size, group_list, newgroup, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8011 )
  
  CALL mpi_comm_create ( base_group, newgroup, groupid, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8012 )
#endif
  
END SUBROUTINE mp_group

!******************************************************************************

!..mp_cart_create
SUBROUTINE mp_cart_create ( comm_old, ndims, dims, pos, comm_cart )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: comm_old, ndims
  INTEGER, INTENT ( OUT ) :: dims ( : ), pos ( : ), comm_cart
  
! Locals
  INTEGER :: ierr, nodes
  LOGICAL :: period(1:ndims), reorder
  
  ierr = 0
  dims ( 1:ndims ) = 1
  pos ( 1:ndims ) = 1
  comm_cart = comm_old
#if defined(__parallel)
  dims ( 1:ndims ) = 0
  
  CALL mpi_comm_size ( comm_old, nodes, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8020 )
  
  CALL mpi_dims_create ( nodes, ndims, dims, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8021 )
  
  reorder = .TRUE.
  period = .TRUE.
  CALL mpi_cart_create ( comm_old, ndims, dims, period, reorder, comm_cart, &
       ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8022 )
  
  CALL mpi_cart_get ( comm_cart, ndims, dims, period, pos, ierr )
  IF ( ierr /= 0 ) CALL mp_stop ( 8023 )
#endif
  
END SUBROUTINE mp_cart_create

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
  CALL mpi_bcast(msg,msglen,mpi_integer,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8101)
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
  CALL mpi_bcast(msg,msglen,mpi_integer,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8102)
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
  CALL mpi_bcast(msg,msglen,mpi_integer,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8102)
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
  CALL mpi_bcast(msg,msglen,mpi_double_precision,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8111)
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
  CALL mpi_bcast(msg,msglen,mpi_double_precision,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8112)
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
  CALL mpi_bcast(msg,msglen,mpi_double_precision,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8113)
#endif
END SUBROUTINE mp_bcast_rm
SUBROUTINE mp_bcast_c1(msg,source,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ) :: msg
  INTEGER :: source, gid
  INTEGER :: msglen, ierr
  
  ierr = 0
  msglen = 1
#if defined(__parallel)
  CALL mpi_bcast(msg,msglen,mpi_double_complex,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8121)
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
  CALL mpi_bcast(msg,msglen,mpi_double_complex,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8122)
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
  CALL mpi_bcast(msg,msglen,mpi_double_complex,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8123)
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
  CALL mpi_bcast(msg,msglen,mpi_logical,source,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8130)
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
! this is a workaround to avoid problems on the T3E
! at the moment we have a data alignment error when trying to
! broadcats characters on the T3E (not always!)
! JH 19/3/99 on galileo
! CALL mpi_bcast(msg,msglen,mpi_character,source,gid,ierr)
  ALLOCATE (imsg(1:msglen))
  DO i = 1, msglen
     imsg(i) = ichar(msg(i:i))
  END DO
  CALL mpi_bcast(imsg,msglen,mpi_integer,source,gid,ierr)
  DO i = 1, msglen
     msg(i:i) = char(imsg(i))
  END DO
  DEALLOCATE (imsg)
  IF ( ierr /= 0 ) CALL mp_stop(8140)
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
  CALL mpi_allreduce(msg,res,msglen,mpi_integer,mpi_sum,gid,ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop(8200)
#endif
END SUBROUTINE mp_sum_i1
SUBROUTINE mp_sum_iv(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER, ALLOCATABLE :: res ( : )
#if defined(__parallel)
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8201)
  CALL mpi_allreduce(msg,res,msglen,mpi_integer,mpi_sum,gid,ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop(8202)
#endif
END SUBROUTINE mp_sum_iv
SUBROUTINE mp_sum_r1(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ) :: res
  
  msglen = 1
#if defined(__parallel)
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,gid, &
       ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop(8203)
#endif
END SUBROUTINE mp_sum_r1
SUBROUTINE mp_sum_rv(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8204)
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop(8205)
#endif
END SUBROUTINE mp_sum_rv
SUBROUTINE mp_sum_rm(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  msglen = size(msg)
  m1 = size(msg(:,1))
  m2 = size(msg(1,:))
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8204)
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop(8205)
#endif
END SUBROUTINE mp_sum_rm
SUBROUTINE mp_sum_c1(msg,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  COMPLEX ( dbl ) :: res
  
  msglen = 2
#if defined(__parallel)
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,gid, &
       ierr)
  msg = res
  IF ( ierr /= 0 ) CALL mp_stop(8205)
#endif
END SUBROUTINE mp_sum_c1
SUBROUTINE mp_sum_cv(msg,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  COMPLEX ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8206)
  CALL mpi_allreduce(msg,res,2*msglen,mpi_double_precision,mpi_sum,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop(8207)
#endif
END SUBROUTINE mp_sum_cv
SUBROUTINE mp_sum_cm(msg,gid)
  IMPLICIT NONE
  COMPLEX ( dbl ), INTENT ( INOUT ) :: msg ( :, : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, m1, m2, ierr
  COMPLEX ( dbl ), ALLOCATABLE :: res ( :, : )
#if defined(__parallel)
  msglen = 2*size(msg)
  m1 = size(msg(:,1))
  m2 = size(msg(1,:))
  ALLOCATE (res(m1,m2),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8204)
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_sum,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop(8208)
#endif
END SUBROUTINE mp_sum_cm

!******************************************************************************

SUBROUTINE mp_max_i(msg,gid)
  IMPLICIT NONE
  INTEGER, INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  INTEGER :: res
#if defined(__parallel)
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,mpi_integer,mpi_max,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8300)
  msg = res
#endif
END SUBROUTINE mp_max_i
SUBROUTINE mp_max_r(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ) :: res
#if defined(__parallel)
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_max,gid, &
       ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8301)
  msg = res
#endif
END SUBROUTINE mp_max_r
SUBROUTINE mp_max_rv(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8302)
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_max,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop(8303)
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
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,mpi_integer,mpi_min,gid,ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8310)
  msg = res
#endif
END SUBROUTINE mp_min_i
SUBROUTINE mp_min_r(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ) :: res
#if defined(__parallel)
  msglen = 1
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_min,gid, &
       ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8311)
  msg = res
#endif
END SUBROUTINE mp_min_r
SUBROUTINE mp_min_rv(msg,gid)
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: msg ( : )
  INTEGER, INTENT ( IN ) :: gid
  INTEGER :: msglen, ierr
  REAL ( dbl ), ALLOCATABLE :: res ( : )
#if defined(__parallel)
  msglen = size(msg)
  ALLOCATE (res(1:msglen),STAT=ierr)
  IF ( ierr /= 0 ) CALL mp_stop(8312)
  CALL mpi_allreduce(msg,res,msglen,mpi_double_precision,mpi_min,gid, &
       ierr)
  msg = res
  DEALLOCATE (res)
  IF ( ierr /= 0 ) CALL mp_stop(8313)
#endif
END SUBROUTINE mp_min_rv

!******************************************************************************

END MODULE message_passing
