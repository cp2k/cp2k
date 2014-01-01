!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Shift around the data in msg
!> \par Example
!>      msg will be moved from rank to rank+displ_in (in a circular way)
!> \par Limitations
!>      * displ_in will be 1 by default (others not tested)
!>      * the message array needs to be the same size on all processes
!> \param[in,out] msg         Rank-2 data to shift
!> \param[in] group           message passing environment identifier
!> \param[in] displ_in        displacements (?)
! *****************************************************************************
  SUBROUTINE mp_shift_cm( msg, group, displ_in)

    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(IN), OPTIONAL            :: displ_in

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_shift_cm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: displ, ierror, left, &
                                                msglen, myrank, nprocs, &
                                                right, tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierror = 0

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL mpi_comm_rank(group,myrank,ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_rank @ "//routineN )
    CALL mpi_comm_size(group,nprocs,ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_size @ "//routineN )
    IF (PRESENT(displ_in)) THEN
       displ=displ_in
    ELSE
       displ=1
    ENDIF
    right=MODULO(myrank+displ,nprocs)
    left =MODULO(myrank-displ,nprocs)
    tag=17
    msglen = SIZE(msg)
    CALL mpi_sendrecv_replace(msg,msglen,MPI_COMPLEX,right,tag,left,tag, &
         group,status(1),ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_sendrecv_replace @ "//routineN )
    DEALLOCATE(status)
#endif

  END SUBROUTINE mp_shift_cm


! *****************************************************************************
!> \brief Shift around the data in msg
!> \par Example
!>      msg will be moved from rank to rank+displ_in (in a circular way)
!> \par Limitations
!>      * displ_in will be 1 by default (others not tested)
!>      * the message array needs to be the same size on all processes
!> \param[in,out] msg         Data to shift
!> \param[in] group           message passing environment identifier
!> \param[in] displ_in        displacements (?)
! *****************************************************************************
  SUBROUTINE mp_shift_c( msg, group, displ_in)

    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(IN), OPTIONAL            :: displ_in

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_shift_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: displ, ierror, left, &
                                                msglen, myrank, nprocs, &
                                                right, tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierror = 0

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL mpi_comm_rank(group,myrank,ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_rank @ "//routineN )
    CALL mpi_comm_size(group,nprocs,ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_comm_size @ "//routineN )
    IF (PRESENT(displ_in)) THEN
       displ=displ_in
    ELSE
       displ=1
    ENDIF
    right=MODULO(myrank+displ,nprocs)
    left =MODULO(myrank-displ,nprocs)
    tag=19
    msglen = SIZE(msg)
    CALL mpi_sendrecv_replace(msg,msglen,MPI_COMPLEX,right,tag,left,&
         tag,group,status(1),ierror)
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_sendrecv_replace @ "//routineN )
    DEALLOCATE(status)
#endif

  END SUBROUTINE mp_shift_c

! *****************************************************************************
!> \brief All-to-all data exchange, rank-1 data of different sizes
!> \par MPI mapping
!>      mpi_alltoallv
!> \par Array sizes
!>      The scount, rcount, and the sdispl and rdispl arrays have a
!>      size equal to the number of processes.
!> \par Offsets
!>      Values in sdispl and rdispl start with 0.
!> \param[in] sb              Data to send
!> \param[in] scount          Data counts for data sent to other processes
!> \param[in] sdispl          Respective data offsets for data sent to process
!> \param[in,out] rb          Buffer into which to receive data
!> \param[in] rcount          Data counts for data received from other
!>                            processes
!> \param[in] rdispl          Respective data offsets for data received from
!>                            other processes
!> \param[in] group           Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_alltoall_c11v ( sb, scount, sdispl, rb, rcount, rdispl, group )

    COMPLEX(kind=real_4), DIMENSION(:), INTENT(IN)        :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    COMPLEX(kind=real_4), DIMENSION(:), INTENT(INOUT)     :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c11v', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr,i

    ierr = 0
#if defined(__parallel)
    CALL mpi_alltoallv ( sb, scount, sdispl, MPI_COMPLEX, &
         rb, rcount, rdispl, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ "//routineN )
#else
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(rcount,rdispl,sdispl,rb,sb)
    DO i=1,rcount(1)
       rb(rdispl(1)+i)=sb(sdispl(1)+i)
    ENDDO
#endif

  END SUBROUTINE mp_alltoall_c11v

! *****************************************************************************
!> \brief All-to-all data exchange, rank-2 data of different sizes
!> \par MPI mapping
!>      mpi_alltoallv
!> \sa mp_alltoall_c11v
! *****************************************************************************
  SUBROUTINE mp_alltoall_c22v ( sb, scount, sdispl, rb, rcount, rdispl, group )

    COMPLEX(kind=real_4), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    COMPLEX(kind=real_4), DIMENSION(:, :), &
      INTENT(INOUT)                          :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c22v', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_alltoallv ( sb, scount, sdispl, MPI_COMPLEX, &
         rb, rcount, rdispl, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ "//routineN )
#else
    rb=sb
#endif

  END SUBROUTINE mp_alltoall_c22v

! *****************************************************************************
!> \brief All-to-all data exchange, rank 1 arrays, equal sizes
!> \par Index meaning
!> \par The first two indices specify the data while the last index counts
!>      the processes
!> \par Sizes of ranks
!>      All processes have the same data size.
!> \par MPI mapping
!>      mpi_alltoall
!> \param[in] sb    array with data to send
!> \param[out] rb   array into which data is received
!> \param[in] count  number of elements to send/receive (product of the
!>                   extents of the first two dimensions)
!> \param[in] group           Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_alltoall_c ( sb, rb, count, group )

    COMPLEX(kind=real_4), DIMENSION(:), INTENT(IN)        :: sb
    COMPLEX(kind=real_4), DIMENSION(:), INTENT(OUT)       :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_alltoall ( sb, count, MPI_COMPLEX, &
         rb, count, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
#else
    rb=sb
#endif

  END SUBROUTINE mp_alltoall_c

! *****************************************************************************
!> \brief All-to-all data exchange, rank-2 arrays, equal sizes
!> \sa mp_alltoall_c
! *****************************************************************************
  SUBROUTINE mp_alltoall_c22 ( sb, rb, count, group )

    COMPLEX(kind=real_4), DIMENSION(:, :), INTENT(IN)     :: sb
    COMPLEX(kind=real_4), DIMENSION(:, :), INTENT(OUT)    :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c22', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_alltoall ( sb, count, MPI_COMPLEX, &
         rb, count, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
#else
    rb=sb
#endif

  END SUBROUTINE mp_alltoall_c22

! *****************************************************************************
!> \brief All-to-all data exchange, rank-3 data with equal sizes
!> \sa mp_alltoall_c
! *****************************************************************************
  SUBROUTINE mp_alltoall_c33 ( sb, rb, count, group )

    COMPLEX(kind=real_4), DIMENSION(:, :, :), INTENT(IN)  :: sb
    COMPLEX(kind=real_4), DIMENSION(:, :, :), INTENT(OUT) :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c33', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_alltoall ( sb, count, MPI_COMPLEX, &
         rb, count, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
#else
    rb=sb
#endif

  END SUBROUTINE mp_alltoall_c33


! *****************************************************************************
!> \brief All-to-all data exchange, rank 4 data, equal sizes
!> \sa mp_alltoall_c
! *****************************************************************************
  SUBROUTINE mp_alltoall_c44 ( sb, rb, count, group )

    COMPLEX(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: sb
    COMPLEX(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c44', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_alltoall ( sb, count, MPI_COMPLEX, &
         rb, count, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
#else
    rb=sb
#endif

  END SUBROUTINE mp_alltoall_c44

! *****************************************************************************
!> \brief All-to-all data exchange, rank-4 data to rank-5 data
!> \note User must ensure size consistency.
!> \sa mp_alltoall_c
! *****************************************************************************
  SUBROUTINE mp_alltoall_c45 ( sb, rb, count, group )

    COMPLEX(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: sb
    COMPLEX(kind=real_4), &
      DIMENSION(:, :, :, :, :), INTENT(OUT)  :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c45', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_alltoall ( sb, count, MPI_COMPLEX, &
         rb, count, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
#endif

  END SUBROUTINE mp_alltoall_c45

! *****************************************************************************
!> \brief All-to-all data exchange, rank-3 data to rank-4 data
!> \note User must ensure size consistency.
!> \sa mp_alltoall_c
! *****************************************************************************
  SUBROUTINE mp_alltoall_c34 ( sb, rb, count, group )

    COMPLEX(kind=real_4), DIMENSION(:, :, :), &
      INTENT(IN)                             :: sb
    COMPLEX(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c34', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_alltoall ( sb, count, MPI_COMPLEX, &
         rb, count, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
#endif

  END SUBROUTINE mp_alltoall_c34

! *****************************************************************************
!> \brief All-to-all data exchange, rank-5 data to rank-4 data
!> \note User must ensure size consistency.
!> \sa mp_alltoall_c
! *****************************************************************************
  SUBROUTINE mp_alltoall_c54 ( sb, rb, count, group )

    COMPLEX(kind=real_4), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: sb
    COMPLEX(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_c54', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL mpi_alltoall ( sb, count, MPI_COMPLEX, &
         rb, count, MPI_COMPLEX, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
#endif

  END SUBROUTINE mp_alltoall_c54


! *****************************************************************************
!> \brief Send one datum to another process
!> \par MPI mapping
!>      mpi_send
!> \param[in] msg             Dum to send
!> \param[in] dest            Destination process
!> \param[in] tag             Transfer identifier
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_send_c(msg,dest,tag,gid)
    COMPLEX(kind=real_4)                                  :: msg
    INTEGER                                  :: dest, tag, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_send_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
    msglen = 1
#if defined(__parallel)
    CALL mpi_send(msg,msglen,MPI_COMPLEX,dest,tag,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_send @ "//routineN )
#endif
  END SUBROUTINE mp_send_c



! *****************************************************************************
!> \brief Send rank-1 data to another process
!> \param[in] msg             Rank-1 data to send
!> \sa mp_send_c
! *****************************************************************************
  SUBROUTINE mp_send_cv(msg,dest,tag,gid)
    COMPLEX(kind=real_4)                                  :: msg( : )
    INTEGER                                  :: dest, tag, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_send_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
    msglen = SIZE(msg)
#if defined(__parallel)
    CALL mpi_send(msg,msglen,MPI_COMPLEX,dest,tag,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_send @ "//routineN )
#endif
  END SUBROUTINE mp_send_cv


! *****************************************************************************
!> \brief Receive one datum from another process
!> \par MPI mapping
!>      mpi_send
!> \param[in,out] msg         Place received data into this variable
!> \param[in,out] source      Process to receieve from
!> \param[in,out] tag         Transfer identifier
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_recv_c(msg,source,tag,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(INOUT)                   :: source, tag
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_recv_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
    msglen = 1
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL mpi_recv(msg,msglen,MPI_COMPLEX,source,tag,gid,status,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_recv @ "//routineN )
    source = status(MPI_SOURCE)
    tag = status(MPI_TAG)
    DEALLOCATE(status)
#endif
  END SUBROUTINE mp_recv_c


! *****************************************************************************
!> \brief Receive rank-1 data from another process
!> \param[in,out] msg         Place receieved data into this rank-1 array
!> \sa mp_recv_c
! *****************************************************************************
  SUBROUTINE mp_recv_cv(msg,source,tag,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(INOUT)                   :: source, tag
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_recv_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
    msglen = SIZE(msg)
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL mpi_recv(msg,msglen,MPI_COMPLEX,source,tag,gid,status,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_recv @ "//routineN )
    source = status(MPI_SOURCE)
    tag = status(MPI_TAG)
    DEALLOCATE(status)
#endif
  END SUBROUTINE mp_recv_cv

! *****************************************************************************
!> \brief Broadcasts a datum to all processes.
!> \par MPI mapping
!>      mpi_bcast
!> \param[in] msg             Datum to broadcast
!> \param[in] source          Processes which broadcasts
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_bcast_c(msg,source,gid)
    COMPLEX(kind=real_4)                                  :: msg
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
    msglen = 1
#if defined(__parallel)
    CALL mpi_bcast(msg,msglen,MPI_COMPLEX,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
#endif
  END SUBROUTINE mp_bcast_c

! *****************************************************************************
!> \brief Broadcasts rank-1 data to all processes
!> \param[in] msg             Data to broadcast
!> \sa mp_bcast_c1
! *****************************************************************************
  SUBROUTINE mp_bcast_cv(msg,source,gid)
    COMPLEX(kind=real_4)                                  :: msg( : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
    msglen = SIZE(msg)
#if defined(__parallel)
    CALL mpi_bcast(msg,msglen,MPI_COMPLEX,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
#endif
  END SUBROUTINE mp_bcast_cv

! *****************************************************************************
!> \brief Broadcasts rank-2 data to all processes
!> \param[in] msg             Data to broadcast
!> \sa mp_bcast_c1
! *****************************************************************************
  SUBROUTINE mp_bcast_cm(msg,source,gid)
    COMPLEX(kind=real_4)                                  :: msg( :, : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_im', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
    msglen = SIZE(msg)
#if defined(__parallel)
    CALL mpi_bcast(msg,msglen,MPI_COMPLEX,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
#endif
  END SUBROUTINE mp_bcast_cm

! *****************************************************************************
!> \brief Broadcasts rank-2 data to all processes
!> \param[in] msg             Data to broadcast
!> \sa mp_bcast_c1
! *****************************************************************************
  SUBROUTINE mp_bcast_c3(msg,source,gid)
    COMPLEX(kind=real_4)                                  :: msg( :, :, : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_c3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
    msglen = SIZE(msg)
#if defined(__parallel)
    CALL mpi_bcast(msg,msglen,MPI_COMPLEX,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
#endif
  END SUBROUTINE mp_bcast_c3


! *****************************************************************************
!> \brief Sums a datum from all processes with result left on all processes.
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Datum to sum (input) and result (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sum_c(msg,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    COMPLEX(kind=real_4)                                  :: res

    ierr = 0
    msglen = 1
#if defined(__parallel)
    CALL mpi_allreduce(msg,res,msglen,MPI_COMPLEX,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
#endif
  END SUBROUTINE mp_sum_c


! *****************************************************************************
!> \brief Element-wise sum of a rank-1 array on all processes.
!> \sa mp_sum_c
!> \param[in,out] msg         Vector to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_cv(msg,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    COMPLEX(kind=real_4), ALLOCATABLE                     :: res( : )

    ierr = 0
#if defined(__parallel)
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,msglen,MPI_COMPLEX,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
    END IF
#endif
  END SUBROUTINE mp_sum_cv


! *****************************************************************************
!> \brief Element-wise sum of a rank-2 array on all processes.
!> \sa mp_sum_c
!> \param[in] msg             Matrix to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_cm(msg,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_cm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, m1, m2, msglen
    COMPLEX(kind=real_4), ALLOCATABLE                     :: res( :, : )

    ierr = 0
#if defined(__parallel)
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    m1 = SIZE(msg,1)
    m2 = SIZE(msg,2)
    ALLOCATE (res(m1,m2),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,msglen,MPI_COMPLEX,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
    END IF
#endif
  END SUBROUTINE mp_sum_cm

! *****************************************************************************
!> \brief Element-wise sum of a rank-3 array on all processes.
!> \sa mp_sum_c
!> \param[in] msg             Arary to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_cm3(msg,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_cm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, m1, m2, m3, &
                                                msglen
    COMPLEX(kind=real_4), ALLOCATABLE                     :: res( :, :, : )

    ierr = 0
#if defined(__parallel)
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    m1 = SIZE(msg,1)
    m2 = SIZE(msg,2)
    m3 = SIZE(msg,3)
    ALLOCATE (res(m1,m2,m3),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,msglen,MPI_COMPLEX,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
    END IF
#endif
  END SUBROUTINE mp_sum_cm3


! *****************************************************************************
!> \brief Element-wise sum of data from all processes with result left only on
!>        one.
!> \par MPI mapping
!>      mpi_reduce
!> \param[in,out] msg         Vector to sum (input) and (only on process root)
!>                            result (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sum_root_cv(msg,root,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_root_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, m1, msglen, &
                                                taskid
    COMPLEX(kind=real_4), ALLOCATABLE                     :: res( : )

    ierr = 0
#if defined(__parallel)
    CALL mpi_comm_rank ( gid, taskid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ "//routineN )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    m1 = SIZE(msg,1)
    ALLOCATE (res(m1),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_reduce(msg,res,msglen,MPI_COMPLEX,MPI_SUM,&
         root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ "//routineN )
    IF ( taskid == root ) THEN
       msg = res
    END IF
    DEALLOCATE (res)
    END IF
#endif
  END SUBROUTINE mp_sum_root_cv

! *****************************************************************************
!> \brief Element-wise sum of data from all processes with result left only on
!>        one.
!> \param[in,out] msg         Matrix to sum (input) and (only on process root)
!>                            result (output)
!> \sa mp_sum_root_cv
! *****************************************************************************
  SUBROUTINE mp_sum_root_cm(msg,root,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_root_rm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, m1, m2, msglen, &
                                                taskid
    COMPLEX(kind=real_4), ALLOCATABLE                     :: res( :, : )

    ierr = 0
#if defined(__parallel)
    CALL mpi_comm_rank ( gid, taskid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ "//routineN )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    m1 = SIZE(msg,1)
    m2 = SIZE(msg,2)
    ALLOCATE (res(m1,m2),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_reduce(msg,res,msglen,MPI_COMPLEX,MPI_SUM,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ "//routineN )
    IF ( taskid == root ) THEN
       msg = res
    END IF
    DEALLOCATE (res)
    END IF
#endif
  END SUBROUTINE mp_sum_root_cm


! *****************************************************************************
!> \brief Finds the maximum of a datum with the result left on all processes.
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Find maximum among these data (input) and
!>                            maximum (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_max_c(msg,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_max_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    COMPLEX(kind=real_4)                                  :: res

    ierr = 0
#if defined(__parallel)
    msglen = 1
    CALL mpi_allreduce(msg,res,msglen,MPI_COMPLEX,MPI_MAX,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
#endif
  END SUBROUTINE mp_max_c

! *****************************************************************************
!> \brief Finds the element-wise maximum of a vector with the result left on
!>        all processes.
!> \param[in,out] msg         Find maximum among these data (input) and
!>                            maximum (output)
!> \sa mp_max_c
! *****************************************************************************
  SUBROUTINE mp_max_cv(msg,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_max_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    COMPLEX(kind=real_4), ALLOCATABLE                     :: res( : )

    ierr = 0
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,msglen,MPI_COMPLEX,MPI_MAX,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
#endif
  END SUBROUTINE mp_max_cv


! *****************************************************************************
!> \brief Finds the minimum of a datum with the result left on all processes.
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Find minimum among these data (input) and
!>                            maximum (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_min_c(msg,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_min_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    COMPLEX(kind=real_4)                                  :: res

    ierr = 0
#if defined(__parallel)
    msglen = 1
    CALL mpi_allreduce(msg,res,msglen,MPI_COMPLEX,MPI_MIN,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
#endif
  END SUBROUTINE mp_min_c

! *****************************************************************************
!> \brief Finds the element-wise minimum of vector with the result left on
!>        all processes.
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Find minimum among these data (input) and
!>                            maximum (output)
!> \sa mp_min_c
! *****************************************************************************
  SUBROUTINE mp_min_cv(msg,gid)
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_min_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen
    COMPLEX(kind=real_4), ALLOCATABLE                     :: res( : )

    ierr = 0
#if defined(__parallel)
    msglen = SIZE(msg)
    ALLOCATE (res(1:msglen),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(msg,res,msglen,MPI_COMPLEX,MPI_MIN,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    msg = res
    DEALLOCATE (res)
#endif
  END SUBROUTINE mp_min_cv


! *****************************************************************************
!> \brief Scatters data from one processes to all others
!> \par MPI mapping
!>      mpi_scatter
!> \param[in] msg_scatter     Data to scatter (for root process)
!> \param[out] msg            Received data
!> \param[in] root            Process which scatters data
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_scatter_cv(msg_scatter,msg,root,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msg_scatter(:)
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msg( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_scatter_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
#if defined(__parallel)
    msglen = SIZE(msg)
    CALL mpi_scatter(msg_scatter,msglen,MPI_COMPLEX,msg,&
         msglen,MPI_COMPLEX,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_scatter @ "//routineN )
#else
    msg = msg_scatter
#endif
  END SUBROUTINE mp_scatter_cv

! *****************************************************************************
!> \brief Gathers a datum from all processes to one
!> \par MPI mapping
!>      mpi_gather
!> \param[in] msg             Datum to send to root
!> \param[out] msg_gather     Received data (on root)
!> \param[in] root            Process which gathers the data
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_gather_c(msg,msg_gather,root,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msg
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msg_gather( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gather_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
#if defined(__parallel)
    msglen = 1
    CALL mpi_gather(msg,msglen,MPI_COMPLEX,msg_gather,&
         msglen,MPI_COMPLEX,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ "//routineN )
#else
    msg_gather = msg
#endif
  END SUBROUTINE mp_gather_c

! *****************************************************************************
!> \brief Gathers data from all processes to one
!> \par Data length
!>      All data (msg) is equal-sized
!> \par MPI mapping
!>      mpi_gather
!> \param[in] msg             Datum to send to root
!> \sa mp_gather_c
! *****************************************************************************
  SUBROUTINE mp_gather_cv(msg,msg_gather,root,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msg( : )
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msg_gather( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gather_c]v', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen

    ierr = 0
#if defined(__parallel)
    msglen = SIZE(msg)
    CALL mpi_gather(msg,msglen,MPI_COMPLEX,msg_gather,&
         msglen,MPI_COMPLEX,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ "//routineN )
#else
    msg_gather = msg
#endif
  END SUBROUTINE mp_gather_cv


! *****************************************************************************
!> \brief Gathers data from all processes to one.
!> \par Data length
!>      Data can have different lengths
!> \par Offsets
!>      Offsets start at 0
!> \par MPI mapping
!>      mpi_gather
!> \param[in] sendbuf         Data to send to root
!> \param[out] recvbuf        Received data (on root)
!> \param[in] recvcounts      Sizes of data received from processes
!> \param[in] displs          Offsets of data received from processes
!> \param[in] root            Process which gathers the data
!> \param[in] comm            Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_gatherv_cv(sendbuf,recvbuf,recvcounts,displs,root,comm)

    COMPLEX(kind=real_4), DIMENSION(:), INTENT(IN)        :: sendbuf
    COMPLEX(kind=real_4), DIMENSION(:), INTENT(OUT)       :: recvbuf
    INTEGER, DIMENSION(:), INTENT(IN)        :: recvcounts, displs
    INTEGER, INTENT(IN)                      :: root, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gatherv_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, sendcount

    ierr = 0
    sendcount = SIZE(sendbuf)
#if defined(__parallel)
    CALL mpi_gatherv(sendbuf,sendcount,MPI_COMPLEX,&
         recvbuf,recvcounts,displs,MPI_COMPLEX,&
         root,comm,ierr)
    IF (ierr /= 0) CALL mp_stop(ierr,"mpi_gatherv @ "//routineN)
#else
    recvbuf(1+displs(1):) = sendbuf
#endif
  END SUBROUTINE mp_gatherv_cv


! *****************************************************************************
!> \brief Gathers a datum from all processes and all processes receive the
!>        same data
!> \par Data size
!>      All processes send equal-sized data
!> \par MPI mapping
!>      mpi_allgather
!> \param[in] msgout          Datum to send
!> \param[out] msgin          Received data
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_allgather_c(msgout,msgin,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgout
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, rcount, scount

    ierr = 0

#if defined(__parallel)
    scount = 1
    rcount = 1
    CALL MPI_ALLGATHER(msgout, scount, MPI_COMPLEX, &
                       msgin , rcount, MPI_COMPLEX, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin = msgout
#endif
  END SUBROUTINE mp_allgather_c

! *****************************************************************************
!> \brief Gathers vector data from all processes and all processes receive the
!>        same data
!> \par Data size
!>      All processes send equal-sized data
!> \par Ranks
!>      The last rank counts the processes
!> \par MPI mapping
!>      mpi_allgather
!> \param[in] msgout          Rank-1 data to send
!> \param[out] msgin          Received data
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_allgather_c12(msgout, msgin,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgout(:)
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgin(:, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_c12', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, rcount, scount

    ierr = 0

#if defined(__parallel)
    scount = SIZE (msgout(:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_COMPLEX, &
                       msgin , rcount, MPI_COMPLEX, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,1) = msgout(:)
#endif
  END SUBROUTINE mp_allgather_c12

! *****************************************************************************
!> \brief Gathers matrix data from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Rank-2 data to send
!> \sa mp_allgather_c12
! *****************************************************************************
  SUBROUTINE mp_allgather_c23(msgout, msgin,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgout(:,:)
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgin(:, :, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_c23', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, rcount, scount

    ierr = 0

#if defined(__parallel)
    scount = SIZE (msgout(:,:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_COMPLEX, &
                       msgin , rcount, MPI_COMPLEX, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,:,1) = msgout(:,:)
#endif
  END SUBROUTINE mp_allgather_c23

! *****************************************************************************
!> \brief Gathers rank-3 data from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Rank-3 data to send
!> \sa mp_allgather_c12
! *****************************************************************************
  SUBROUTINE mp_allgather_c34(msgout, msgin,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgout(:,:, :)
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgin(:, :, :, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_c34', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, rcount, scount

    ierr = 0

#if defined(__parallel)
    scount = SIZE (msgout(:,:,:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_COMPLEX, &
                       msgin , rcount, MPI_COMPLEX, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,:,:,1) = msgout(:,:,:)
#endif
  END SUBROUTINE mp_allgather_c34

! *****************************************************************************
!> \brief Gathers vector data from all processes and all processes receive the
!>        same data
!> \par Data size
!>      Processes can send different-sized data
!> \par Ranks
!>      The last rank counts the processes
!> \par Offsets
!>      Offsets are from 0
!> \par MPI mapping
!>      mpi_allgather
!> \param[in] msgout          Rank-1 data to send
!> \param[out] msgin          Received data
!> \param[in] rcount          Size of sent data for every process
!> \param[in] rdispl          Offset of sent data for every process
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_allgatherv_cv(msgout,msgin,rcount,rdispl,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgout( : )
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: rcount( : ), rdispl( : ), gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgatherv_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, scount

    ierr = 0

#if defined(__parallel)
    scount = SIZE ( msgout )
    CALL MPI_ALLGATHERV(msgout, scount, MPI_COMPLEX, msgin, rcount, &
                        rdispl, MPI_COMPLEX, gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgatherv @ "//routineN )
#else
    msgin = msgout
#endif
  END SUBROUTINE mp_allgatherv_cv


! *****************************************************************************
!> \brief Sums a vector and partitions the result among processes
!> \param[in] msgout          Data to sum
!> \param[out] msgin          Received portion of summed data
!> \param[in] rcount          Partition sizes of the summed data for
!>                            every process
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sum_scatter_cv(msgout,msgin,rcount,gid)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgout( : )
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: rcount( : ), gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_scatter_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr

    ierr = 0

#if defined(__parallel)
    CALL MPI_REDUCE_SCATTER(msgout, msgin, rcount, MPI_COMPLEX, MPI_SUM, &
         gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce_scatter @ "//routineN )

#else
    msgin = msgout
#endif
  END SUBROUTINE mp_sum_scatter_cv


! *****************************************************************************
!> \brief Sends and receives vector data
!> \param[in] msgin           Data to send
!> \param[in] dest            Process to send data to
!> \param[out] msgout         Received data
!> \param[in] source          Process from which to receive
!> \param[in] comm            Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sendrecv_cv(msgin,dest,msgout,source,comm)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgin( : )
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgout( : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen_in, &
                                                msglen_out, recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    msglen_in = SIZE(msgin)
    msglen_out = SIZE(msgout)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_COMPLEX,dest,send_tag,msgout,&
         msglen_out,MPI_COMPLEX,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
  END SUBROUTINE mp_sendrecv_cv


! *****************************************************************************
!> \brief Sends and receives matrix data
!> \sa mp_sendrecv_cv
! *****************************************************************************
  SUBROUTINE mp_sendrecv_cm2(msgin,dest,msgout,source,comm)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgin( :, : )
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgout( :, : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_cm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen_in, &
                                                msglen_out, recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    msglen_in = SIZE(msgin,1)*SIZE(msgin,2)
    msglen_out = SIZE(msgout,1)*SIZE(msgout,2)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_COMPLEX,dest,send_tag,msgout,&
         msglen_out,MPI_COMPLEX,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
  END SUBROUTINE mp_sendrecv_cm2


! *****************************************************************************
!> \brief Sends and receives rank-3 data
!> \sa mp_sendrecv_cv
! *****************************************************************************
  SUBROUTINE mp_sendrecv_cm3(msgin,dest,msgout,source,comm)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msgin( :, :, : )
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msgout( :, :, : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_cm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, msglen_in, &
                                                msglen_out, recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    msglen_in = SIZE(msgin)
    msglen_out = SIZE(msgout)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_COMPLEX,dest,send_tag,msgout,&
         msglen_out,MPI_COMPLEX,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
  END SUBROUTINE mp_sendrecv_cm3

! *****************************************************************************
!> \brief Non-blocking send and receieve of a vector
!> \par Implementation
!>      Calls mpi_isend and mpi_irecv.
!> \note
!>      The arguments must be pointers to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \par History
!>      11.2004 created [Joost VandeVondele]
!> \param[in] msgin           Vector data to send
!> \param[in] dest            Which process to send to
!> \param[out] msgout         Receive data into this pointer
!> \param[in] source          Process to receive from
!> \param[in] comm            Message passing environment identifier
!> \param[out] send_request   Request handle for the send
!> \param[out] recv_request   Request handle for the receive
!> \param[in] tag             (optional) tag to differentiate requests
! *****************************************************************************
  SUBROUTINE mp_isendrecv_cv(msgin,dest,msgout,source,comm,send_request,&
       recv_request,tag)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: msgin
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: send_request, recv_request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isendrecv_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, lower1, msglen, &
                                                my_tag
    COMPLEX(kind=real_4)                                  :: foo

    ierr = 0
#if defined(__parallel)
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       CALL mpi_irecv(msgout(lower1),msglen,MPI_COMPLEX,source, my_tag,&
            comm,recv_request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_COMPLEX,source, my_tag,&
            comm,recv_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    msglen = SIZE(msgin,1)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       CALL mpi_isend(msgin(lower1),msglen,MPI_COMPLEX,dest,my_tag,&
            comm,send_request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_COMPLEX,dest,my_tag,&
            comm,send_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )
#else
    send_request=0
    recv_request=0
    msgout = msgin
#endif
  END SUBROUTINE mp_isendrecv_cv


! *****************************************************************************
!> \brief Non-blocking send and receieve of a matrix
!> \note
!>      The arguments must be pointers to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \par History
!>      08.2003 created [f&j]
!> \sa mp_isendrecv_cv
! *****************************************************************************
  SUBROUTINE mp_isendrecv_cm2(msgin,dest,msgout,source,comm,send_request,&
       recv_request,tag)
    COMPLEX(kind=real_4), DIMENSION(:, :), POINTER        :: msgin
    INTEGER, INTENT(IN)                      :: dest
    COMPLEX(kind=real_4), DIMENSION(:, :), POINTER        :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: send_request, recv_request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isendrecv_cm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, lower1, lower2, &
                                                msglen, my_tag
    COMPLEX(kind=real_4)                                  :: foo

    ierr = 0

#if defined(__parallel)
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       CALL mpi_irecv(msgout(lower1,lower2),msglen,MPI_COMPLEX,source, my_tag,&
            comm,recv_request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_COMPLEX,source, my_tag,&
            comm,recv_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    msglen = SIZE(msgin,1)*SIZE(msgin,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       CALL mpi_isend(msgin(lower1,lower2),msglen,MPI_COMPLEX,dest,my_tag,&
            comm,send_request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_COMPLEX,dest,my_tag,&
            comm,send_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ mp_isendrecv_rm2" )
#else
    send_request=0
    recv_request=0
    msgout = msgin
#endif
  END SUBROUTINE mp_isendrecv_cm2



! *****************************************************************************
!> \brief Non-blocking send of vector data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \par History
!>      08.2003 created [f&j]
!> \sa mp_isendrecv_cv
! *****************************************************************************
  SUBROUTINE mp_isend_cv(msgin,dest,comm,request,tag)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, lower1, msglen, &
                                                my_tag
    COMPLEX(kind=real_4)                                  :: foo(1)

    ierr = 0

#if defined(__parallel)
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       CALL mpi_isend(msgin(lower1),msglen,MPI_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
  END SUBROUTINE mp_isend_cv


! *****************************************************************************
!> \brief Non-blocking send of matrix data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \author fawzi
!> \par History
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_cv
!> \sa mp_isend_cv
! *****************************************************************************
  SUBROUTINE mp_isend_cm2(msgin,dest,comm,request,tag)
    COMPLEX(kind=real_4), DIMENSION(:, :), POINTER  :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_cm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, lower1, lower2, &
                                                msglen, my_tag
    COMPLEX(kind=real_4)                                  :: foo(1)

    ierr = 0

#if defined(__parallel)
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin,1)*SIZE(msgin,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       CALL mpi_isend(msgin(lower1,lower2),msglen,MPI_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
  END SUBROUTINE mp_isend_cm2


! *****************************************************************************
!> \brief Non-blocking send of rank-3 data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \author fawzi
!> \par History
!>      9.2008 added _rm3 subroutine [Iain Bethune] (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_cv
!> \sa mp_isend_cv
! *****************************************************************************
  SUBROUTINE mp_isend_cm3(msgin,dest,comm,request,tag)
    COMPLEX(kind=real_4), DIMENSION(:, :, :), &
      POINTER                                :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_cm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, lower1, lower2, &
                                                lower3, msglen, my_tag
    COMPLEX(kind=real_4)                                  :: foo(1)

    ierr = 0

#if defined(__parallel)
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin,1)*SIZE(msgin,2)*SIZE(msgin,3)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       lower3=LBOUND(msgin,3)
       CALL mpi_isend(msgin(lower1,lower2,lower3),msglen,MPI_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_COMPLEX,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
  END SUBROUTINE mp_isend_cm3


! *****************************************************************************
!> \brief Non-blocking receive of vector data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \par History
!>      08.2003 created [f&j]
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_cv
! *****************************************************************************
  SUBROUTINE mp_irecv_cv(msgout,source,comm,request,tag)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_cv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, lower1, msglen, &
                                                my_tag
    COMPLEX(kind=real_4)                                  :: foo(1)

    ierr = 0

#if defined(__parallel)
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       CALL mpi_irecv(msgout(lower1),msglen,MPI_COMPLEX,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_COMPLEX,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
  END SUBROUTINE mp_irecv_cv


! *****************************************************************************
!> \brief Non-blocking receive of matrix data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \author fawzi
!> \par History
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_cv
!> \sa mp_irecv_cv
! *****************************************************************************
  SUBROUTINE mp_irecv_cm2(msgout,source,comm,request,tag)
    COMPLEX(kind=real_4), DIMENSION(:, :), POINTER        :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_cm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, lower1, lower2, &
                                                msglen, my_tag
    COMPLEX(kind=real_4)                                  :: foo(1)

    ierr = 0

#if defined(__parallel)
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       CALL mpi_irecv(msgout(lower1,lower2),msglen,MPI_COMPLEX,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_COMPLEX,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
  END SUBROUTINE mp_irecv_cm2


! *****************************************************************************
!> \brief Non-blocking send of rank-3 data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \author fawzi
!> \par History
!>      9.2008 added _rm3 subroutine [Iain Bethune] (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_cv
!> \sa mp_irecv_cv
! *****************************************************************************
  SUBROUTINE mp_irecv_cm3(msgout,source,comm,request,tag)
    COMPLEX(kind=real_4), DIMENSION(:, :, :), &
      POINTER                                :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_cm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, lower1, lower2, &
                                                lower3, msglen, my_tag
    COMPLEX(kind=real_4)                                  :: foo(1)

    ierr = 0

#if defined(__parallel)
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)*SIZE(msgout,3)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       lower3=LBOUND(msgout,3)
       CALL mpi_irecv(msgout(lower1,lower2,lower3),msglen,MPI_COMPLEX,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_COMPLEX,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_ircv @ "//routineN )

#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
  END SUBROUTINE mp_irecv_cm3

! *****************************************************************************
!> \brief Creates an MPI RMA window.
!> \author UB
!> \param[out] window      created window id
!> \param[in] range        window contents
!> \param[in] len          (optional) window size
!> \param[in] gid          message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_win_create_c(window, range, len, gid)
    TYPE(mp_window_type), INTENT(OUT)    :: window
    COMPLEX(kind=real_4), DIMENSION(:), POINTER       :: range
    INTEGER, INTENT(IN)                  :: gid
    INTEGER, INTENT(IN), OPTIONAL        :: len
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_create_c', &
         routineP = moduleN//':'//routineN
    INTEGER                              :: ierr
    INTEGER                              :: data_type_size
#if defined(__parallel)
    ! requires that your MPI installation matches your Fortran compiler... be aware
    ! e.g. on a 64 bit machine, NAG will require MPI_ADDRESS_KIND=4, and gfortran MPI_ADDRESS_KIND=8
    INTEGER(KIND=MPI_ADDRESS_KIND)       :: winsize
#endif
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    ierr = 0
    IF (PRESENT (len)) THEN
       winsize = len
    ELSE
       winsize = SIZE (range)
    ENDIF
    CALL MPI_TYPE_SIZE(MPI_COMPLEX, data_type_size, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_type_size @ "//routineN)
    winsize = winsize*data_type_size
    CALL MPI_WIN_CREATE(range, winsize, data_type_size,&
      MPI_INFO_NULL, gid, window%id, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_create @ "//routineN)
#else
    NULLIFY(window%src_r, window%src_d)
    NULLIFY(window%src_c, window%src_z)
    NULLIFY(window%src_i, window%src_l)
    window%src_c => range
#endif

  END SUBROUTINE mp_win_create_c


! *****************************************************************************
!> \brief Fetch access using RMA
!> \author UB
!> \param[in] window       RMA access window
!> \param[in] remote       remote process
!> \param[in] offset       offset in remote window (starting at 0)
!> \param[in] len          amount of data to fetch
!> \param[out] storage     integer storage of fetched data
! *****************************************************************************
  SUBROUTINE mp_rma_get_c(window, remote, offset, len, storage)
    TYPE(mp_window_type), INTENT(IN) :: window
    INTEGER, INTENT(IN)              :: remote
    INTEGER, INTENT(IN)              :: offset
    COMPLEX(kind=real_4), DIMENSION(:), POINTER   :: storage
    INTEGER, INTENT(IN)              :: len
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_rma_get_c', &
         routineP = moduleN//':'//routineN
    INTEGER                   :: ierr
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    ierr = 0
    CALL MPI_GET(storage, len, MPI_COMPLEX, remote,&
         INT (offset, MPI_ADDRESS_KIND), len, MPI_COMPLEX, window%id, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_get @ "//routineN)
#else
    storage(LBOUND(storage,1):LBOUND(storage,1)+len-1) &
         = window%src_c(LBOUND(window%src_i,1)+offset&
         : LBOUND(window%src_i,1)+offset+len-1)
#endif
  END SUBROUTINE mp_rma_get_c


! *****************************************************************************
!> \brief Store access using RMA
!> \author UB
!> \param[in] window       RMA access window
!> \param[in] storage      integer data to store
!> \param[in] len          amount of data to store
!> \param[in] remote       remote process
!> \param[in] offset       offset in remote window (starting at 0)
! *****************************************************************************
  SUBROUTINE mp_rma_put_c(window, storage, len, remote, offset)
    TYPE(mp_window_type), INTENT(IN) :: window
    INTEGER, INTENT(IN)              :: remote
    INTEGER, INTENT(IN)              :: offset
    COMPLEX(kind=real_4), DIMENSION(:), POINTER   :: storage
    INTEGER, INTENT(IN)              :: len
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_rma_put_c', &
         routineP = moduleN//':'//routineN
    INTEGER                        :: ierr
#if !defined (__parallel)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: win_data
#endif
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    ierr = 0
    CALL MPI_PUT(storage, len, MPI_COMPLEX, remote,&
         INT (offset, MPI_ADDRESS_KIND), len, MPI_COMPLEX, window, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_put @ "//routineN)
#else
    win_data => window%src_c
    win_data(LBOUND(window%src_i,1)+offset &
         : LBOUND(window%src_i,1)+offset+len-1) =&
    storage(LBOUND(storage,1):LBOUND(storage,1)+len-1)
#endif
  END SUBROUTINE mp_rma_put_c


! *****************************************************************************
!> \brief Allocates special parallel memory
!> \author UB
!> \param[in]  data      pointer to integer array to allocate
!> \param[in]  len       number of integers to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_allocate_c(DATA, len, stat)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER      :: DATA
    INTEGER, INTENT(IN)                 :: len
    INTEGER, INTENT(OUT), OPTIONAL      :: stat
    INTEGER                   :: ierr
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allocate_c', &
         routineP = moduleN//':'//routineN
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    !$ IF (OMP_GET_THREAD_NUM() .NE. 0) call mp_stop(ierr, routineN)
    ierr = 0
    NULLIFY(DATA)
    CALL mp_alloc_mem(DATA, len, stat=ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "mpi_alloc_mem @ "//routineN)
    ENDIF
#else
    ALLOCATE(DATA(len), stat=ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "ALLOCATE @ "//routineN)
    ENDIF
#endif
  END SUBROUTINE mp_allocate_c


! *****************************************************************************
!> \brief Deallocates special parallel memory
!> \author UB
!> \param[in] data         pointer to special memory to deallocate
! *****************************************************************************
  SUBROUTINE mp_deallocate_c(DATA, stat)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER      :: DATA
    INTEGER, INTENT(OUT), OPTIONAL      :: stat
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_deallocate_c', &
         routineP = moduleN//':'//routineN
    INTEGER                   :: ierr
!   ---------------------------------------------------------------------------
#if defined(__parallel)
    !$ IF (OMP_GET_THREAD_NUM() .NE. 0) call mp_stop(ierr, routineN)
    ierr = 0
    CALL mp_free_mem(DATA, ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "mpi_free_mem @ "//routineN)
    ENDIF
    NULLIFY(DATA)
#else
    DEALLOCATE(DATA, stat=ierr)
    IF (PRESENT (stat)) THEN
       stat=ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "DEALLOCATE @ "//routineN)
    ENDIF
    NULLIFY(DATA)
#endif
  END SUBROUTINE mp_deallocate_c

  FUNCTION mp_type_make_c (ptr,&
       vector_descriptor, index_descriptor) &
       RESULT (type_descriptor)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER                    :: ptr
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL       :: vector_descriptor
    TYPE(mp_indexing_meta_type), INTENT(IN), OPTIONAL :: index_descriptor
    TYPE(mp_type_descriptor_type)                     :: type_descriptor

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_make_type_c', &
         routineP = moduleN//':'//routineN
    INTEGER :: ierr
!   ---------------------------------------------------------------------------
    NULLIFY (type_descriptor%subtype)
    type_descriptor%length = SIZE (ptr)
#if defined(__parallel)
    type_descriptor%type_handle = MPI_COMPLEX
    CALL MPI_Get_address (ptr, type_descriptor%base, ierr)
    IF (ierr /= 0) CALL mp_abort("MPI_Get_address @ "//routineN)
#else
    type_descriptor%type_handle = 5
#endif
    type_descriptor%vector_descriptor(1:2) = 1
    type_descriptor%has_indexing = .FALSE.
    type_descriptor%data_c => ptr
    IF (PRESENT (vector_descriptor) .OR. PRESENT (index_descriptor)) THEN
       CALL mp_abort (routineN//": Vectors and indices NYI")
    ENDIF
  END FUNCTION mp_type_make_c

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
  SUBROUTINE mp_file_write_at_cv(fh, offset, msg, msglen)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msg(:)
    INTEGER, INTENT(IN)                        :: fh 
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_cv', &
                                   routineP = moduleN//':'//routineN

    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT(fh, offset, msg, msg_len, MPI_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_cv @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg(1:msg_len)
#endif
  END SUBROUTINE mp_file_write_at_cv

  SUBROUTINE mp_file_write_at_c(fh, offset, msg)
    COMPLEX(kind=real_4), INTENT(IN)               :: msg
    INTEGER, INTENT(IN)                        :: fh 
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_c', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT(fh, offset, msg, 1, MPI_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_c @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_write_at_c
! *****************************************************************************
!> \brief (parallel) Blocking collective file write using explicit offsets
!>        (serial) Unformatted stream write
!> \par MPI-I/O mapping   mpi_file_write_at_all 
!> \par STREAM-I/O mapping   WRITE   
! *****************************************************************************
  SUBROUTINE mp_file_write_at_all_cv(fh, offset, msg, msglen)
    COMPLEX(kind=real_4), INTENT(IN)                      :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_all_cv', &
                                   routineP = moduleN//':'//routineN

    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT_ALL(fh, offset, msg, msg_len, MPI_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_all_cv @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg(1:msg_len)
#endif
  END SUBROUTINE mp_file_write_at_all_cv

  SUBROUTINE mp_file_write_at_all_c(fh, offset, msg)
    COMPLEX(kind=real_4), INTENT(IN)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_at_all_c', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_WRITE_AT_ALL(fh, offset, msg, 1, MPI_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_write_at_all_c @ "//routineN)
    DEALLOCATE(status)
#else
    WRITE(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_write_at_all_c
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
  SUBROUTINE mp_file_read_at_cv(fh, offset, msg, msglen)
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_cv', &
                                   routineP = moduleN//':'//routineN

    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT(fh, offset, msg, msg_len, MPI_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_cv @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg(1:msg_len)
#endif
  END SUBROUTINE mp_file_read_at_cv

  SUBROUTINE mp_file_read_at_c(fh, offset, msg)
    COMPLEX(kind=real_4), INTENT(OUT)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_c', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT(fh, offset, msg, 1, MPI_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_c @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_c
! *****************************************************************************
!> \brief (parallel) Blocking collective file read using explicit offsets
!>        (serial) Unformatted stream read
!> \par MPI-I/O mapping    mpi_file_read_at_all 
!> \par STREAM-I/O mapping   READ  
! *****************************************************************************
  SUBROUTINE mp_file_read_at_all_cv(fh, offset, msg, msglen)
    COMPLEX(kind=real_4), INTENT(OUT)                     :: msg(:)
    INTEGER, INTENT(IN)                        :: fh
    INTEGER, INTENT(IN), OPTIONAL              :: msglen
    INTEGER                                    :: ierr, msg_len
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_all_cv', &
                                   routineP = moduleN//':'//routineN

    msg_len = SIZE(msg)
    IF (PRESENT(msglen)) msg_len = msglen
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT_ALL(fh, offset, msg, msg_len, MPI_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_all_cv @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg(1:msg_len)
#endif
  END SUBROUTINE mp_file_read_at_all_cv

  SUBROUTINE mp_file_read_at_all_c(fh, offset, msg)
    COMPLEX(kind=real_4), INTENT(OUT)               :: msg
    INTEGER, INTENT(IN)                        :: fh
    INTEGER                                    :: ierr
    INTEGER, ALLOCATABLE, DIMENSION(:)         :: status
    INTEGER(kind=file_offset), INTENT(IN)        :: offset
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_at_all_c', &
                                   routineP = moduleN//':'//routineN

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    CALL MPI_FILE_READ_AT_ALL(fh, offset, msg, 1, MPI_COMPLEX, status, ierr)
    IF (ierr .NE. 0) CALL mp_abort("mpi_file_read_at_all_c @ "//routineN)
    DEALLOCATE(status)
#else
    READ(UNIT=fh, POS=offset+1) msg
#endif
  END SUBROUTINE mp_file_read_at_all_c
