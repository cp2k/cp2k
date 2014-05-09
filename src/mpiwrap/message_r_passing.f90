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
  SUBROUTINE mp_shift_rm( msg, group, displ_in)

    REAL(kind=real_4), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(IN), OPTIONAL            :: displ_in

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_shift_rm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: displ, handle, ierror, left, &
                                                msglen, myrank, nprocs, &
                                                right, tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierror = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


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
    t_start = m_walltime ( )
    CALL mpi_sendrecv_replace(msg,msglen,MPI_REAL,right,tag,left,tag, &
         group,status(1),ierror)
    t_end = m_walltime ( )
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_sendrecv_replace @ "//routineN )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
    DEALLOCATE(status)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_shift_rm


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
  SUBROUTINE mp_shift_r( msg, group, displ_in)

    REAL(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(IN), OPTIONAL            :: displ_in

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_shift_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: displ, handle, ierror, left, &
                                                msglen, myrank, nprocs, &
                                                right, tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierror = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


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
    t_start = m_walltime ( )
    CALL mpi_sendrecv_replace(msg,msglen,MPI_REAL,right,tag,left,&
         tag,group,status(1),ierror)
    t_end = m_walltime ( )
    IF ( ierror /= 0 ) CALL mp_stop ( ierror, "mpi_sendrecv_replace @ "//routineN )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
    DEALLOCATE(status)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_shift_r

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
  SUBROUTINE mp_alltoall_r11v ( sb, scount, sdispl, rb, rcount, rdispl, group )

    REAL(kind=real_4), DIMENSION(:), INTENT(IN)        :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    REAL(kind=real_4), DIMENSION(:), INTENT(INOUT)     :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r11v', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen, i

    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    ierr = 0
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoallv ( sb, scount, sdispl, MPI_REAL, &
         rb, rcount, rdispl, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ "//routineN )
    t_end = m_walltime ( )
    msglen = SUM ( scount ) + SUM ( rcount )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) SHARED(rcount,rdispl,sdispl,rb,sb)
    DO i=1,rcount(1)
       rb(rdispl(1)+i)=sb(sdispl(1)+i)
    ENDDO
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r11v

! *****************************************************************************
!> \brief All-to-all data exchange, rank-2 data of different sizes
!> \par MPI mapping
!>      mpi_alltoallv
!> \sa mp_alltoall_r11v
! *****************************************************************************
  SUBROUTINE mp_alltoall_r22v ( sb, scount, sdispl, rb, rcount, rdispl, group )

    REAL(kind=real_4), DIMENSION(:, :), &
      INTENT(IN)                             :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    REAL(kind=real_4), DIMENSION(:, :), &
      INTENT(INOUT)                          :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    INTEGER, INTENT(IN)                      :: group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r22v', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoallv ( sb, scount, sdispl, MPI_REAL, &
         rb, rcount, rdispl, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoallv @ "//routineN )
    msglen = SUM ( scount ) + SUM ( rcount )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*2*real_4_size)
#else
    rb=sb
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r22v

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
  SUBROUTINE mp_alltoall_r ( sb, rb, count, group )

    REAL(kind=real_4), DIMENSION(:), INTENT(IN)        :: sb
    REAL(kind=real_4), DIMENSION(:), INTENT(OUT)       :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen, np

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_REAL, &
         rb, count, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    rb=sb
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r

! *****************************************************************************
!> \brief All-to-all data exchange, rank-2 arrays, equal sizes
!> \sa mp_alltoall_r
! *****************************************************************************
  SUBROUTINE mp_alltoall_r22 ( sb, rb, count, group )

    REAL(kind=real_4), DIMENSION(:, :), INTENT(IN)     :: sb
    REAL(kind=real_4), DIMENSION(:, :), INTENT(OUT)    :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r22', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen, np

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_REAL, &
         rb, count, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * SIZE(sb) * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    rb=sb
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r22

! *****************************************************************************
!> \brief All-to-all data exchange, rank-3 data with equal sizes
!> \sa mp_alltoall_r
! *****************************************************************************
  SUBROUTINE mp_alltoall_r33 ( sb, rb, count, group )

    REAL(kind=real_4), DIMENSION(:, :, :), INTENT(IN)  :: sb
    REAL(kind=real_4), DIMENSION(:, :, :), INTENT(OUT) :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r33', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen, np

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_REAL, &
         rb, count, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    rb=sb
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r33


! *****************************************************************************
!> \brief All-to-all data exchange, rank 4 data, equal sizes
!> \sa mp_alltoall_r
! *****************************************************************************
  SUBROUTINE mp_alltoall_r44 ( sb, rb, count, group )

    REAL(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: sb
    REAL(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r44', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen, np

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_REAL, &
         rb, count, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    rb=sb
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r44

! *****************************************************************************
!> \brief All-to-all data exchange, rank-4 data to rank-5 data
!> \note User must ensure size consistency.
!> \sa mp_alltoall_r
! *****************************************************************************
  SUBROUTINE mp_alltoall_r45 ( sb, rb, count, group )

    REAL(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: sb
    REAL(kind=real_4), &
      DIMENSION(:, :, :, :, :), INTENT(OUT)  :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r45', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen, np

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_REAL, &
         rb, count, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r45

! *****************************************************************************
!> \brief All-to-all data exchange, rank-3 data to rank-4 data
!> \note User must ensure size consistency.
!> \sa mp_alltoall_r
! *****************************************************************************
  SUBROUTINE mp_alltoall_r34 ( sb, rb, count, group )

    REAL(kind=real_4), DIMENSION(:, :, :), &
      INTENT(IN)                             :: sb
    REAL(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r34', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen, np

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_REAL, &
         rb, count, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r34

! *****************************************************************************
!> \brief All-to-all data exchange, rank-5 data to rank-4 data
!> \note User must ensure size consistency.
!> \sa mp_alltoall_r
! *****************************************************************************
  SUBROUTINE mp_alltoall_r54 ( sb, rb, count, group )

    REAL(kind=real_4), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: sb
    REAL(kind=real_4), DIMENSION(:, :, :, :), &
      INTENT(OUT)                            :: rb
    INTEGER, INTENT(IN)                      :: count, group

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_alltoall_r54', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen, np

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_alltoall ( sb, count, MPI_REAL, &
         rb, count, MPI_REAL, group, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_alltoall @ "//routineN )
    CALL mpi_comm_size ( group, np, ierr )
    IF ( ierr /= 0 ) CALL mp_stop ( ierr, "mpi_comm_size @ "//routineN )
    msglen = 2 * count * np
    t_end = m_walltime ( )
    CALL add_perf(perf_id=6,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_alltoall_r54


! *****************************************************************************
!> \brief Send one datum to another process
!> \par MPI mapping
!>      mpi_send
!> \param[in] msg             Dum to send
!> \param[in] dest            Destination process
!> \param[in] tag             Transfer identifier
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_send_r(msg,dest,tag,gid)
    REAL(kind=real_4)                                  :: msg
    INTEGER                                  :: dest, tag, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_send_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_send(msg,msglen,MPI_REAL,dest,tag,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_send @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=13,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_send_r



! *****************************************************************************
!> \brief Send rank-1 data to another process
!> \param[in] msg             Rank-1 data to send
!> \sa mp_send_r
! *****************************************************************************
  SUBROUTINE mp_send_rv(msg,dest,tag,gid)
    REAL(kind=real_4)                                  :: msg( : )
    INTEGER                                  :: dest, tag, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_send_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_send(msg,msglen,MPI_REAL,dest,tag,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_send @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=13,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_send_rv


! *****************************************************************************
!> \brief Receive one datum from another process
!> \par MPI mapping
!>      mpi_send
!> \param[in,out] msg         Place received data into this variable
!> \param[in,out] source      Process to receieve from
!> \param[in,out] tag         Transfer identifier
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_recv_r(msg,source,tag,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(INOUT)                   :: source, tag
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_recv_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    CALL mpi_recv(msg,msglen,MPI_REAL,source,tag,gid,status,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_recv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=14,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
    source = status(MPI_SOURCE)
    tag = status(MPI_TAG)
    DEALLOCATE(status)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_recv_r


! *****************************************************************************
!> \brief Receive rank-1 data from another process
!> \param[in,out] msg         Place receieved data into this rank-1 array
!> \sa mp_recv_r
! *****************************************************************************
  SUBROUTINE mp_recv_rv(msg,source,tag,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(INOUT)                   :: source, tag
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_recv_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    CALL mpi_recv(msg,msglen,MPI_REAL,source,tag,gid,status,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_recv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=14,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
    source = status(MPI_SOURCE)
    tag = status(MPI_TAG)
    DEALLOCATE(status)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_recv_rv

! *****************************************************************************
!> \brief Broadcasts a datum to all processes.
!> \par MPI mapping
!>      mpi_bcast
!> \param[in] msg             Datum to broadcast
!> \param[in] source          Processes which broadcasts
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_bcast_r(msg,source,gid)
    REAL(kind=real_4)                                  :: msg
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_REAL,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_bcast_r

! *****************************************************************************
!> \brief Broadcasts rank-1 data to all processes
!> \param[in] msg             Data to broadcast
!> \sa mp_bcast_r1
! *****************************************************************************
  SUBROUTINE mp_bcast_rv(msg,source,gid)
    REAL(kind=real_4)                                  :: msg( : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_REAL,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_bcast_rv

! *****************************************************************************
!> \brief Broadcasts rank-2 data to all processes
!> \param[in] msg             Data to broadcast
!> \sa mp_bcast_r1
! *****************************************************************************
  SUBROUTINE mp_bcast_rm(msg,source,gid)
    REAL(kind=real_4)                                  :: msg( :, : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_im', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_REAL,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_bcast_rm

! *****************************************************************************
!> \brief Broadcasts rank-3 data to all processes
!> \param[in] msg             Data to broadcast
!> \sa mp_bcast_r1
! *****************************************************************************
  SUBROUTINE mp_bcast_r3(msg,source,gid)
    REAL(kind=real_4)                                  :: msg( :, :, : )
    INTEGER                                  :: source, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = SIZE(msg)
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_bcast(msg,msglen,MPI_REAL,source,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_bcast @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=2,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_bcast_r3


! *****************************************************************************
!> \brief Sums a datum from all processes with result left on all processes.
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Datum to sum (input) and result (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sum_r(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    msglen = 1
#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_r


! *****************************************************************************
!> \brief Element-wise sum of a rank-1 array on all processes.
!> \sa mp_sum_r
!> \param[in,out] msg         Vector to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_rv(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_rv


! *****************************************************************************
!> \brief Element-wise sum of a rank-2 array on all processes.
!> \sa mp_sum_r
!> \param[in] msg             Matrix to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_rm(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_rm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, m1, m2, msglen, step
    INTEGER, PARAMETER :: max_msg=2**25

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    ! chunk up the call so that message sizes are limited, to avoid overflows in mpich triggered in large rpa calcs
    step=MAX(1,SIZE(msg,2)/MAX(1,SIZE(msg)/max_msg))
    DO m1=LBOUND(msg,2),UBOUND(msg,2), step
       msglen = SIZE(msg,1)*(MIN(UBOUND(msg,2),m1+step-1)-m1+1)
       IF (msglen>0) THEN
          CALL mpi_allreduce(MPI_IN_PLACE,msg(LBOUND(msg,1),m1),msglen,MPI_REAL,MPI_SUM,gid,ierr)
          IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
       END IF
    ENDDO
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_rm

! *****************************************************************************
!> \brief Element-wise sum of a rank-3 array on all processes.
!> \sa mp_sum_r
!> \param[in] msg             Arary to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_rm3(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_rm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, m1, m2, m3, &
                                                msglen
    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_rm3

! *****************************************************************************
!> \brief Element-wise sum of a rank-4 array on all processes.
!> \sa mp_sum_r
!> \param[in] msg             Arary to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_rm4(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( :, :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_rm4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, m1, m2, m3, m4, &
                                                msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_rm4

! *****************************************************************************
!> \brief Element-wise sum of a rank-5 array on all processes.
!> \sa mp_sum_r
!> \param[in] msg             Arary to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_rm5(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( :, :, :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_rm5', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, m1, m2, m3, m4, m5, &
                                                msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_rm5

! *****************************************************************************
!> \brief Element-wise sum of a rank-6 array on all processes.
!> \sa mp_sum_r
!> \param[in] msg             Arary to sum and result
! *****************************************************************************
  SUBROUTINE mp_sum_rm6(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( :, :, :, :, :, : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_rm6', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, m1, m2, m3, m4, m5, m6, &
                                                msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_SUM,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_rm6


! *****************************************************************************
!> \brief Element-wise sum of data from all processes with result left only on
!>        one.
!> \par MPI mapping
!>      mpi_reduce
!> \param[in,out] msg         Vector to sum (input) and (only on process root)
!>                            result (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sum_root_rv(msg,root,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_root_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, m1, msglen, &
                                                taskid
    REAL(kind=real_4), ALLOCATABLE                     :: res( : )

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_comm_rank ( gid, taskid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ "//routineN )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    m1 = SIZE(msg,1)
    ALLOCATE (res(m1),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_reduce(msg,res,msglen,MPI_REAL,MPI_SUM,&
         root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ "//routineN )
    IF ( taskid == root ) THEN
       msg = res
    END IF
    DEALLOCATE (res)
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_root_rv

! *****************************************************************************
!> \brief Element-wise sum of data from all processes with result left only on
!>        one.
!> \param[in,out] msg         Matrix to sum (input) and (only on process root)
!>                            result (output)
!> \sa mp_sum_root_rv
! *****************************************************************************
  SUBROUTINE mp_sum_root_rm(msg,root,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( :, : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_root_rm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, m1, m2, msglen, &
                                                taskid
    REAL(kind=real_4), ALLOCATABLE                     :: res( :, : )

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    CALL mpi_comm_rank ( gid, taskid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_comm_rank @ "//routineN )
    msglen = SIZE(msg)
    IF (msglen>0) THEN
    m1 = SIZE(msg,1)
    m2 = SIZE(msg,2)
    ALLOCATE (res(m1,m2),STAT=ierr)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_reduce(msg,res,msglen,MPI_REAL,MPI_SUM,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce @ "//routineN )
    IF ( taskid == root ) THEN
       msg = res
    END IF
    DEALLOCATE (res)
    END IF
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_root_rm


! *****************************************************************************
!> \brief Finds the maximum of a datum with the result left on all processes.
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Find maximum among these data (input) and
!>                            maximum (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_max_r(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_max_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = 1
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_MAX,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_max_r

! *****************************************************************************
!> \brief Finds the element-wise maximum of a vector with the result left on
!>        all processes.
!> \param[in,out] msg         Find maximum among these data (input) and
!>                            maximum (output)
!> \sa mp_max_r
! *****************************************************************************
  SUBROUTINE mp_max_rv(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_max_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_MAX,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_max_rv


! *****************************************************************************
!> \brief Finds the minimum of a datum with the result left on all processes.
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Find minimum among these data (input) and
!>                            maximum (output)
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_min_r(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_min_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = 1
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_MIN,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_min_r

! *****************************************************************************
!> \brief Finds the element-wise minimum of vector with the result left on
!>        all processes.
!> \par MPI mapping
!>      mpi_allreduce
!> \param[in,out] msg         Find minimum among these data (input) and
!>                            maximum (output)
!> \sa mp_min_r
! *****************************************************************************
  SUBROUTINE mp_min_rv(msg,gid)
    REAL(kind=real_4), INTENT(INOUT)                   :: msg( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_min_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    IF ( ierr /= 0 ) CALL mp_abort( "allocate @ "//routineN )
    CALL mpi_allreduce(MPI_IN_PLACE,msg,msglen,MPI_REAL,MPI_MIN,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allreduce @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_min_rv


! *****************************************************************************
!> \brief Scatters data from one processes to all others
!> \par MPI mapping
!>      mpi_scatter
!> \param[in] msg_scatter     Data to scatter (for root process)
!> \param[out] msg            Received data
!> \param[in] root            Process which scatters data
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_scatter_rv(msg_scatter,msg,root,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msg_scatter(:)
    REAL(kind=real_4), INTENT(OUT)                     :: msg( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_scatter_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    CALL mpi_scatter(msg_scatter,msglen,MPI_REAL,msg,&
         msglen,MPI_REAL,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_scatter @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    msg = msg_scatter
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_scatter_rv

! *****************************************************************************
!> \brief Gathers a datum from all processes to one
!> \par MPI mapping
!>      mpi_gather
!> \param[in] msg             Datum to send to root
!> \param[out] msg_gather     Received data (on root)
!> \param[in] root            Process which gathers the data
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_gather_r(msg,msg_gather,root,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msg
    REAL(kind=real_4), INTENT(OUT)                     :: msg_gather( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gather_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = 1
    CALL mpi_gather(msg,msglen,MPI_REAL,msg_gather,&
         msglen,MPI_REAL,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    msg_gather = msg
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_gather_r

! *****************************************************************************
!> \brief Gathers data from all processes to one
!> \par Data length
!>      All data (msg) is equal-sized
!> \par MPI mapping
!>      mpi_gather
!> \param[in] msg             Datum to send to root
!> \sa mp_gather_r
! *****************************************************************************
  SUBROUTINE mp_gather_rv(msg,msg_gather,root,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msg( : )
    REAL(kind=real_4), INTENT(OUT)                     :: msg_gather( : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gather_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    CALL mpi_gather(msg,msglen,MPI_REAL,msg_gather,&
         msglen,MPI_REAL,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    msg_gather = msg
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_gather_rv

! *****************************************************************************
!> \brief Gathers data from all processes to one
!> \par Data length
!>      All data (msg) is equal-sized
!> \par MPI mapping
!>      mpi_gather
!> \param[in] msg             Datum to send to root
!> \sa mp_gather_r
! *****************************************************************************
  SUBROUTINE mp_gather_rm(msg,msg_gather,root,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msg( :, : )
    REAL(kind=real_4), INTENT(OUT)                     :: msg_gather( :, : )
    INTEGER, INTENT(IN)                      :: root, gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gather_rm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    msglen = SIZE(msg)
    CALL mpi_gather(msg,msglen,MPI_REAL,msg_gather,&
         msglen,MPI_REAL,root,gid,ierr)
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_gather @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=4,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    msg_gather = msg
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_gather_rm

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
  SUBROUTINE mp_gatherv_rv(sendbuf,recvbuf,recvcounts,displs,root,comm)

    REAL(kind=real_4), DIMENSION(:), INTENT(IN)        :: sendbuf
    REAL(kind=real_4), DIMENSION(:), INTENT(OUT)       :: recvbuf
    INTEGER, DIMENSION(:), INTENT(IN)        :: recvcounts, displs
    INTEGER, INTENT(IN)                      :: root, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_gatherv_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, sendcount

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

    sendcount = SIZE(sendbuf)
#if defined(__parallel)
    t_start = m_walltime()
    CALL mpi_gatherv(sendbuf,sendcount,MPI_REAL,&
         recvbuf,recvcounts,displs,MPI_REAL,&
         root,comm,ierr)
    IF (ierr /= 0) CALL mp_stop(ierr,"mpi_gatherv @ "//routineN)
    t_end = m_walltime()
    CALL add_perf(perf_id=4,&
         count=1,&
         time=t_end-t_start,&
         msg_size=sendcount*real_4_size)
#else
    recvbuf(1+displs(1):) = sendbuf
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_gatherv_rv


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
  SUBROUTINE mp_allgather_r(msgout,msgin,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msgout
    REAL(kind=real_4), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_r', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, rcount, scount

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    scount = 1
    rcount = 1
    CALL MPI_ALLGATHER(msgout, scount, MPI_REAL, &
                       msgin , rcount, MPI_REAL, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin = msgout
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_allgather_r

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
  SUBROUTINE mp_allgather_r12(msgout, msgin,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msgout(:)
    REAL(kind=real_4), INTENT(OUT)                     :: msgin(:, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_r12', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, rcount, scount

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    scount = SIZE (msgout(:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_REAL, &
                       msgin , rcount, MPI_REAL, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,1) = msgout(:)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_allgather_r12

! *****************************************************************************
!> \brief Gathers matrix data from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Rank-2 data to send
!> \sa mp_allgather_r12
! *****************************************************************************
  SUBROUTINE mp_allgather_r23(msgout, msgin,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msgout(:,:)
    REAL(kind=real_4), INTENT(OUT)                     :: msgin(:, :, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_r23', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, rcount, scount

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    scount = SIZE (msgout(:,:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_REAL, &
                       msgin , rcount, MPI_REAL, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,:,1) = msgout(:,:)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_allgather_r23

! *****************************************************************************
!> \brief Gathers rank-3 data from all processes and all processes receive the
!>        same data
!> \param[in] msgout          Rank-3 data to send
!> \sa mp_allgather_r12
! *****************************************************************************
  SUBROUTINE mp_allgather_r34(msgout, msgin,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msgout(:,:, :)
    REAL(kind=real_4), INTENT(OUT)                     :: msgin(:, :, :, :)
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgather_r34', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, rcount, scount

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    scount = SIZE (msgout(:,:,:))
    rcount = scount
    CALL MPI_ALLGATHER(msgout, scount, MPI_REAL, &
                       msgin , rcount, MPI_REAL, &
                       gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgather @ "//routineN )
#else
    msgin(:,:,:,1) = msgout(:,:,:)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_allgather_r34

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
  SUBROUTINE mp_allgatherv_rv(msgout,msgin,rcount,rdispl,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msgout( : )
    REAL(kind=real_4), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: rcount( : ), rdispl( : ), gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allgatherv_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, scount

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    scount = SIZE ( msgout )
    CALL MPI_ALLGATHERV(msgout, scount, MPI_REAL, msgin, rcount, &
                        rdispl, MPI_REAL, gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_allgatherv @ "//routineN )
#else
    msgin = msgout
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_allgatherv_rv


! *****************************************************************************
!> \brief Sums a vector and partitions the result among processes
!> \param[in] msgout          Data to sum
!> \param[out] msgin          Received portion of summed data
!> \param[in] rcount          Partition sizes of the summed data for
!>                            every process
!> \param[in] gid             Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sum_scatter_rv(msgout,msgin,rcount,gid)
    REAL(kind=real_4), INTENT(IN)                      :: msgout( : )
    REAL(kind=real_4), INTENT(OUT)                     :: msgin( : )
    INTEGER, INTENT(IN)                      :: rcount( : ), gid

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_scatter_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    CALL MPI_REDUCE_SCATTER(msgout, msgin, rcount, MPI_REAL, MPI_SUM, &
         gid, ierr )
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_reduce_scatter @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=3,count=1,time=t_end-t_start,&
         msg_size=rcount(1)*2*real_4_size)
#else
    msgin = msgout
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sum_scatter_rv


! *****************************************************************************
!> \brief Sends and receives vector data
!> \param[in] msgin           Data to send
!> \param[in] dest            Process to send data to
!> \param[out] msgout         Received data
!> \param[in] source          Process from which to receive
!> \param[in] comm            Message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_sendrecv_rv(msgin,dest,msgout,source,comm)
    REAL(kind=real_4), INTENT(IN)                      :: msgin( : )
    INTEGER, INTENT(IN)                      :: dest
    REAL(kind=real_4), INTENT(OUT)                     :: msgout( : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen_in, &
                                                msglen_out, recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    msglen_in = SIZE(msgin)
    msglen_out = SIZE(msgout)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_REAL,dest,send_tag,msgout,&
         msglen_out,MPI_REAL,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,&
         msg_size=(msglen_in+msglen_out)*real_4_size/2)
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sendrecv_rv


! *****************************************************************************
!> \brief Sends and receives matrix data
!> \sa mp_sendrecv_rv
! *****************************************************************************
  SUBROUTINE mp_sendrecv_rm2(msgin,dest,msgout,source,comm)
    REAL(kind=real_4), INTENT(IN)                      :: msgin( :, : )
    INTEGER, INTENT(IN)                      :: dest
    REAL(kind=real_4), INTENT(OUT)                     :: msgout( :, : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_rm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen_in, &
                                                msglen_out, recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    msglen_in = SIZE(msgin,1)*SIZE(msgin,2)
    msglen_out = SIZE(msgout,1)*SIZE(msgout,2)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_REAL,dest,send_tag,msgout,&
         msglen_out,MPI_REAL,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,&
         msg_size=(msglen_in+msglen_out)*real_4_size/2)
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sendrecv_rm2


! *****************************************************************************
!> \brief Sends and receives rank-3 data
!> \sa mp_sendrecv_rv
! *****************************************************************************
  SUBROUTINE mp_sendrecv_rm3(msgin,dest,msgout,source,comm)
    REAL(kind=real_4), INTENT(IN)                      :: msgin( :, :, : )
    INTEGER, INTENT(IN)                      :: dest
    REAL(kind=real_4), INTENT(OUT)                     :: msgout( :, :, : )
    INTEGER, INTENT(IN)                      :: source, comm

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_sendrecv_rm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, msglen_in, &
                                                msglen_out, recv_tag, send_tag
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: status

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    ALLOCATE(status(MPI_STATUS_SIZE))
    t_start = m_walltime ( )
    msglen_in = SIZE(msgin)
    msglen_out = SIZE(msgout)
    send_tag = 0 ! cannot think of something better here, this might be dangerous
    recv_tag = 0 ! cannot think of something better here, this might be dangerous
    CALL mpi_sendrecv(msgin,msglen_in,MPI_REAL,dest,send_tag,msgout,&
         msglen_out,MPI_REAL,source,recv_tag,comm,status(1),ierr)
    ! we do not check the status
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_sendrecv @ "//routineN )
    t_end = m_walltime ( )
    CALL add_perf(perf_id=7,count=1,time=t_end-t_start,&
         msg_size=(msglen_in+msglen_out)*real_4_size/2)
    DEALLOCATE(status)
#else
    msgout = msgin
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_sendrecv_rm3

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
  SUBROUTINE mp_isendrecv_rv(msgin,dest,msgout,source,comm,send_request,&
       recv_request,tag)
    REAL(kind=real_4), DIMENSION(:), POINTER           :: msgin
    INTEGER, INTENT(IN)                      :: dest
    REAL(kind=real_4), DIMENSION(:), POINTER           :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: send_request, recv_request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isendrecv_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, lower1, msglen, &
                                                my_tag
    REAL(kind=real_4)                                  :: foo

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       CALL mpi_irecv(msgout(lower1),msglen,MPI_REAL,source, my_tag,&
            comm,recv_request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_REAL,source, my_tag,&
            comm,recv_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    msglen = SIZE(msgin,1)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       CALL mpi_isend(msgin(lower1),msglen,MPI_REAL,dest,my_tag,&
            comm,send_request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_REAL,dest,my_tag,&
            comm,send_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

    msglen = (msglen+SIZE(msgout,1))/2.0_dp
    t_end = m_walltime ( )
    CALL add_perf(perf_id=8,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    send_request=0
    recv_request=0
    msgout = msgin
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_isendrecv_rv


! *****************************************************************************
!> \brief Non-blocking send and receieve of a matrix
!> \note
!>      The arguments must be pointers to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \par History
!>      08.2003 created [f&j]
!> \sa mp_isendrecv_rv
! *****************************************************************************
  SUBROUTINE mp_isendrecv_rm2(msgin,dest,msgout,source,comm,send_request,&
       recv_request,tag)
    REAL(kind=real_4), DIMENSION(:, :), POINTER        :: msgin
    INTEGER, INTENT(IN)                      :: dest
    REAL(kind=real_4), DIMENSION(:, :), POINTER        :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: send_request, recv_request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isendrecv_rm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, lower1, lower2, &
                                                msglen, my_tag
    REAL(kind=real_4)                                  :: foo

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       CALL mpi_irecv(msgout(lower1,lower2),msglen,MPI_REAL,source, my_tag,&
            comm,recv_request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_REAL,source, my_tag,&
            comm,recv_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    msglen = SIZE(msgin,1)*SIZE(msgin,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       CALL mpi_isend(msgin(lower1,lower2),msglen,MPI_REAL,dest,my_tag,&
            comm,send_request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_REAL,dest,my_tag,&
            comm,send_request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ mp_isendrecv_rm2" )

    msglen = (msglen+SIZE(msgout,1)*SIZE(msgout,2))/2.0_dp
    t_end = m_walltime ( )
    CALL add_perf(perf_id=8,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    send_request=0
    recv_request=0
    msgout = msgin
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_isendrecv_rm2



! *****************************************************************************
!> \brief Non-blocking send of vector data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \par History
!>      08.2003 created [f&j]
!> \sa mp_isendrecv_rv
! *****************************************************************************
  SUBROUTINE mp_isend_rv(msgin,dest,comm,request,tag)
    REAL(kind=real_4), DIMENSION(:), POINTER           :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, lower1, msglen, &
                                                my_tag
    REAL(kind=real_4)                                  :: foo(1)

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       CALL mpi_isend(msgin(lower1),msglen,MPI_REAL,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_REAL,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=11,count=1,time=t_end-t_start,msg_size=2*msglen*real_4_size)
#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_isend_rv


! *****************************************************************************
!> \brief Non-blocking send of matrix data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \author fawzi
!> \par History
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_rv
!> \sa mp_isend_rv
! *****************************************************************************
  SUBROUTINE mp_isend_rm2(msgin,dest,comm,request,tag)
    REAL(kind=real_4), DIMENSION(:, :), POINTER  :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_rm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, lower1, lower2, &
                                                msglen, my_tag
    REAL(kind=real_4)                                  :: foo(1)

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin,1)*SIZE(msgin,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       CALL mpi_isend(msgin(lower1,lower2),msglen,MPI_REAL,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_REAL,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=11,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_isend_rm2


! *****************************************************************************
!> \brief Non-blocking send of rank-3 data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \author fawzi
!> \par History
!>      9.2008 added _rm3 subroutine [Iain Bethune] (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_rv
!> \sa mp_isend_rv
! *****************************************************************************
  SUBROUTINE mp_isend_rm3(msgin,dest,comm,request,tag)
    REAL(kind=real_4), DIMENSION(:, :, :), &
      POINTER                                :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_rm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, lower1, lower2, &
                                                lower3, msglen, my_tag
    REAL(kind=real_4)                                  :: foo(1)

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgin,1)*SIZE(msgin,2)*SIZE(msgin,3)
    IF (msglen>0) THEN
       lower1=LBOUND(msgin,1)
       lower2=LBOUND(msgin,2)
       lower3=LBOUND(msgin,3)
       CALL mpi_isend(msgin(lower1,lower2,lower3),msglen,MPI_REAL,dest,my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_isend(foo,msglen,MPI_REAL,dest,my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_isend @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=11,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    ierr=1
    CALL mp_stop( ierr, "mp_isend called in non parallel case" )
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_isend_rm3


! *****************************************************************************
!> \brief Non-blocking receive of vector data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \par History
!>      08.2003 created [f&j]
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_rv
! *****************************************************************************
  SUBROUTINE mp_irecv_rv(msgout,source,comm,request,tag)
    REAL(kind=real_4), DIMENSION(:), POINTER           :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, lower1, msglen, &
                                                my_tag
    REAL(kind=real_4)                                  :: foo(1)

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       CALL mpi_irecv(msgout(lower1),msglen,MPI_REAL,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_REAL,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=12,count=1,time=t_end-t_start,msg_size=2*msglen*real_4_size)
#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_irecv_rv


! *****************************************************************************
!> \brief Non-blocking receive of matrix data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \author fawzi
!> \par History
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_rv
!> \sa mp_irecv_rv
! *****************************************************************************
  SUBROUTINE mp_irecv_rm2(msgout,source,comm,request,tag)
    REAL(kind=real_4), DIMENSION(:, :), POINTER        :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_rm2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, lower1, lower2, &
                                                msglen, my_tag
    REAL(kind=real_4)                                  :: foo(1)

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       CALL mpi_irecv(msgout(lower1,lower2),msglen,MPI_REAL,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_REAL,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_irecv @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=12,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_irecv_rm2


! *****************************************************************************
!> \brief Non-blocking send of rank-3 data
!> \note
!>      The argument must be a pointer to be sure that we do not get
!>      temporaries. They must point to contiguous memory.
!> \author fawzi
!> \par History
!>      9.2008 added _rm3 subroutine [Iain Bethune] (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
!>      2009-11-25 [UB] Made type-generic for templates
!> \sa mp_isendrecv_rv
!> \sa mp_irecv_rv
! *****************************************************************************
  SUBROUTINE mp_irecv_rm3(msgout,source,comm,request,tag)
    REAL(kind=real_4), DIMENSION(:, :, :), &
      POINTER                                :: msgout
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(out)                     :: request
    INTEGER, INTENT(in), OPTIONAL            :: tag

    CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_rm3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, lower1, lower2, &
                                                lower3, msglen, my_tag
    REAL(kind=real_4)                                  :: foo(1)

    ierr = 0
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)


#if defined(__parallel)
    t_start = m_walltime ( )
    my_tag = 0
    IF (PRESENT(tag)) my_tag=tag

    msglen = SIZE(msgout,1)*SIZE(msgout,2)*SIZE(msgout,3)
    IF (msglen>0) THEN
       lower1=LBOUND(msgout,1)
       lower2=LBOUND(msgout,2)
       lower3=LBOUND(msgout,3)
       CALL mpi_irecv(msgout(lower1,lower2,lower3),msglen,MPI_REAL,source, my_tag,&
            comm,request,ierr)
    ELSE
       CALL mpi_irecv(foo,msglen,MPI_REAL,source, my_tag,&
            comm,request,ierr)
    END IF
    IF ( ierr /= 0 ) CALL mp_stop( ierr, "mpi_ircv @ "//routineN )

    t_end = m_walltime ( )
    CALL add_perf(perf_id=12,count=1,time=t_end-t_start,msg_size=msglen*real_4_size)
#else
    CALL mp_abort( "mp_irecv called in non parallel case" )
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_irecv_rm3

! *****************************************************************************
!> \brief Creates an MPI RMA window.
!> \author UB
!> \param[out] window      created window id
!> \param[in] range        window contents
!> \param[in] len          (optional) window size
!> \param[in] gid          message passing environment identifier
! *****************************************************************************
  SUBROUTINE mp_win_create_r(window, range, len, gid)
    TYPE(mp_window_type), INTENT(OUT)    :: window
    REAL(kind=real_4), DIMENSION(:), POINTER       :: range
    INTEGER, INTENT(IN)                  :: gid
    INTEGER, INTENT(IN), OPTIONAL        :: len
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_create_r', &
         routineP = moduleN//':'//routineN
    INTEGER                              :: ierr, handle
    INTEGER                              :: data_type_size
#if defined(__parallel)
    ! requires that your MPI installation matches your Fortran compiler... be aware
    ! e.g. on a 64 bit machine, NAG will require MPI_ADDRESS_KIND=4, and gfortran MPI_ADDRESS_KIND=8
    INTEGER(KIND=MPI_ADDRESS_KIND)       :: winsize
#endif
!   ---------------------------------------------------------------------------
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime()
    ierr = 0
    IF (PRESENT (len)) THEN
       winsize = len
    ELSE
       winsize = SIZE (range)
    ENDIF
    CALL MPI_TYPE_SIZE(MPI_REAL, data_type_size, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_type_size @ "//routineN)
    winsize = winsize*data_type_size
    CALL MPI_WIN_CREATE(range, winsize, data_type_size,&
      MPI_INFO_NULL, gid, window%id, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_create @ "//routineN)
    t_end = m_walltime()
    CALL add_perf(perf_id=20, time=t_end-t_start)
#else
    NULLIFY(window%src_r, window%src_d)
    NULLIFY(window%src_c, window%src_z)
    NULLIFY(window%src_i, window%src_l)
    window%src_r => range
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)

  END SUBROUTINE mp_win_create_r


! *****************************************************************************
!> \brief Fetch access using RMA
!> \author UB
!> \param[in] window       RMA access window
!> \param[in] remote       remote process
!> \param[in] offset       offset in remote window (starting at 0)
!> \param[in] len          amount of data to fetch
!> \param[out] storage     integer storage of fetched data
! *****************************************************************************
  SUBROUTINE mp_rma_get_r(window, remote, offset, len, storage)
    TYPE(mp_window_type), INTENT(IN) :: window
    INTEGER, INTENT(IN)              :: remote
    INTEGER, INTENT(IN)              :: offset
    REAL(kind=real_4), DIMENSION(:), POINTER   :: storage
    INTEGER, INTENT(IN)              :: len
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_rma_get_r', &
         routineP = moduleN//':'//routineN
    INTEGER                   :: ierr, handle
!   ---------------------------------------------------------------------------
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime()
    ierr = 0
    CALL MPI_GET(storage, len, MPI_REAL, remote,&
         INT (offset, MPI_ADDRESS_KIND), len, MPI_REAL, window%id, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_get @ "//routineN)
    t_end = m_walltime()
    CALL add_perf(perf_id=17, count=1, time=t_end-t_start,&
         msg_size=len*real_4_size)
#else
    storage(LBOUND(storage,1):LBOUND(storage,1)+len-1) &
         = window%src_r(LBOUND(window%src_i,1)+offset&
         : LBOUND(window%src_i,1)+offset+len-1)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_rma_get_r


! *****************************************************************************
!> \brief Store access using RMA
!> \author UB
!> \param[in] window       RMA access window
!> \param[in] storage      integer data to store
!> \param[in] len          amount of data to store
!> \param[in] remote       remote process
!> \param[in] offset       offset in remote window (starting at 0)
! *****************************************************************************
  SUBROUTINE mp_rma_put_r(window, storage, len, remote, offset)
    TYPE(mp_window_type), INTENT(IN) :: window
    INTEGER, INTENT(IN)              :: remote
    INTEGER, INTENT(IN)              :: offset
    REAL(kind=real_4), DIMENSION(:), POINTER   :: storage
    INTEGER, INTENT(IN)              :: len
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_rma_put_r', &
         routineP = moduleN//':'//routineN
    INTEGER                        :: ierr, handle
    REAL(kind=real_4), DIMENSION(:), POINTER :: win_data
!   ---------------------------------------------------------------------------
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime()
    ierr = 0
    CALL MPI_PUT(storage, len, MPI_REAL, remote,&
         INT (offset, MPI_ADDRESS_KIND), len, MPI_REAL, window, ierr)
    IF (ierr /= 0) CALL mp_stop(ierr, "mpi_put @ "//routineN)
    t_end = m_walltime()
    CALL add_perf(perf_id=16, count=1, time=t_end-t_start,&
         msg_size=real_4_size*len)
#else
    win_data => window%src_r
    win_data(LBOUND(window%src_i,1)+offset &
         : LBOUND(window%src_i,1)+offset+len-1) =&
    storage(LBOUND(storage,1):LBOUND(storage,1)+len-1)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_rma_put_r


! *****************************************************************************
!> \brief Allocates special parallel memory
!> \author UB
!> \param[in]  data      pointer to integer array to allocate
!> \param[in]  len       number of integers to allocate
!> \param[out] stat      (optional) allocation status result
! *****************************************************************************
  SUBROUTINE mp_allocate_r(DATA, len, stat)
    REAL(kind=real_4), DIMENSION(:), POINTER      :: DATA
    INTEGER, INTENT(IN)                 :: len
    INTEGER, INTENT(OUT), OPTIONAL      :: stat
    INTEGER                   :: ierr, handle
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_allocate_r', &
         routineP = moduleN//':'//routineN
!   ---------------------------------------------------------------------------
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime()
    ierr = 0
    NULLIFY(DATA)
    CALL mp_alloc_mem(DATA, len, stat=ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "mpi_alloc_mem @ "//routineN)
    ENDIF
    t_end = m_walltime()
    CALL add_perf(perf_id=15, count=1, time=t_end-t_start)
#else
    ALLOCATE(DATA(len), stat=ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "ALLOCATE @ "//routineN)
    ENDIF
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_allocate_r


! *****************************************************************************
!> \brief Deallocates special parallel memory
!> \author UB
!> \param[in] data         pointer to special memory to deallocate
! *****************************************************************************
  SUBROUTINE mp_deallocate_r(DATA, stat)
    REAL(kind=real_4), DIMENSION(:), POINTER      :: DATA
    INTEGER, INTENT(OUT), OPTIONAL      :: stat
    CHARACTER(len=*), PARAMETER :: routineN = 'mp_deallocate_r', &
         routineP = moduleN//':'//routineN
    INTEGER                   :: ierr, handle
!   ---------------------------------------------------------------------------
    IF (ASSOCIATED(mp_external_timeset)) CALL mp_external_timeset(routineN,handle)

#if defined(__parallel)
    t_start = m_walltime()
    ierr = 0
    CALL mp_free_mem(DATA, ierr)
    IF (PRESENT (stat)) THEN
       stat = ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "mpi_free_mem @ "//routineN)
    ENDIF
    NULLIFY(DATA)
    t_end = m_walltime()
    CALL add_perf(perf_id=15, count=1, time=t_end-t_start)
#else
    DEALLOCATE(DATA, stat=ierr)
    IF (PRESENT (stat)) THEN
       stat=ierr
    ELSE
       IF (ierr /= 0) CALL mp_stop(ierr, "DEALLOCATE @ "//routineN)
    ENDIF
    NULLIFY(DATA)
#endif
    IF (ASSOCIATED(mp_external_timestop)) CALL mp_external_timestop(handle)
  END SUBROUTINE mp_deallocate_r
