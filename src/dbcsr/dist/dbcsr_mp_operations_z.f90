!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Row/column and global all-to-all
!> \param sb ...
!> \param scount ...
!> \param sdispl ...
!> \param rb ...
!> \param rcount ...
!> \param rdispl ...
!> \param[in] mp_env         MP Environment
!> \param[in] most_ptp       (optional) Use point-to-point for row/column;
!>                           default is no
!> \param[in] remainder_ptp  (optional) Use point-to-point for remaining;
!>                           default is no
!> \param[in] no_hybrid      (optional) Use regular global collective; default
!>                           is no
!> \par Communicator selection
!>      Uses row and column communicators for row/column
!>      sends. Remaining sends are performed using the global
!>      communicator.  Point-to-point isend/irecv are used if ptp is
!>      set, otherwise a alltoall collective call is issued.
!>      see mp_alltoall 
! *****************************************************************************
  SUBROUTINE hybrid_alltoall_z1 (sb, scount, sdispl,&
       rb, rcount, rdispl, mp_env, most_ptp, remainder_ptp, no_hybrid)
    COMPLEX(kind=real_8), DIMENSION(:), &
      INTENT(in), TARGET                     :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    COMPLEX(kind=real_8), DIMENSION(:), &
      INTENT(INOUT), TARGET                  :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL, INTENT(in), OPTIONAL            :: most_ptp, remainder_ptp,&
                                                no_hybrid

    CHARACTER(len=*), PARAMETER :: routineN = 'hybrid_alltoall_z', &
      routineP = moduleN//':'//routineN

    INTEGER :: all_group, mynode, mypcol, myprow, nall_rr, nall_sr, ncol_rr, &
      ncol_sr, npcols, nprows, nrow_rr, nrow_sr, numnodes, dst, src,&
      prow, pcol, send_cnt, recv_cnt, tag, grp, i
    INTEGER, ALLOCATABLE, DIMENSION(:) :: all_rr, all_sr, col_rr, col_sr, &
      new_rcount, new_rdispl, new_scount, new_sdispl, row_rr, row_sr
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: most_collective, &
                                                remainder_collective, no_h
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: send_data_p, recv_data_p
    TYPE(dbcsr_mp_obj)                       :: mpe

    !CALL dbcsr_assert (mp_env%mp%subgroups_defined, dbcsr_warning_level,&
    !     dbcsr_caller_error, routineN, "Row/col communicators undefined.")
    IF (.NOT. dbcsr_mp_has_subgroups (mp_env)) THEN
       mpe = mp_env
       CALL dbcsr_mp_grid_setup (mpe)
    ENDIF
    most_collective = .TRUE.
    remainder_collective = .TRUE.
    no_h = .FALSE.
    IF (PRESENT (most_ptp)) most_collective = .NOT. most_ptp
    IF (PRESENT (remainder_ptp)) remainder_collective = .NOT. remainder_ptp
    IF (PRESENT (no_hybrid)) no_h = no_hybrid
    all_group = dbcsr_mp_group (mp_env)
    ! Don't use subcommunicators if they're not defined.
    no_h = no_h .OR. .NOT. dbcsr_mp_has_subgroups (mp_env) .OR. .NOT. has_MPI
    subgrouped: IF (mp_env%mp%subgroups_defined .AND. .NOT. no_h) THEN
       mynode = dbcsr_mp_mynode (mp_env)
       numnodes = dbcsr_mp_numnodes (mp_env)
       nprows = dbcsr_mp_nprows (mp_env)
       npcols = dbcsr_mp_npcols (mp_env)
       myprow = dbcsr_mp_myprow (mp_env)
       mypcol = dbcsr_mp_mypcol (mp_env)
       pgrid => dbcsr_mp_pgrid (mp_env)
       ALLOCATE (row_sr(0:npcols-1)) ; nrow_sr = 0
       ALLOCATE (row_rr(0:npcols-1)) ; nrow_rr = 0
       ALLOCATE (col_sr(0:nprows-1)) ; ncol_sr = 0
       ALLOCATE (col_rr(0:nprows-1)) ; ncol_rr = 0
       ALLOCATE (all_sr(0:numnodes-1)) ; nall_sr = 0
       ALLOCATE (all_rr(0:numnodes-1)) ; nall_rr = 0
       ALLOCATE (new_scount(numnodes), new_rcount(numnodes))
       ALLOCATE (new_sdispl(numnodes), new_rdispl(numnodes))
       IF (.NOT.remainder_collective) THEN
          CALL remainder_point_to_point ()
       ENDIF
       IF (.NOT.most_collective) THEN
          CALL most_point_to_point ()
       ELSE
          CALL most_alltoall ()
       ENDIF
       IF (remainder_collective) THEN
          CALL remainder_alltoall ()
       ENDIF
       ! Wait for all issued sends and receives.
       IF (.NOT.most_collective) THEN
          CALL mp_waitall (row_sr(0:nrow_sr-1))
          CALL mp_waitall (col_sr(0:ncol_sr-1))
          CALL mp_waitall (row_rr(0:nrow_rr-1))
          CALL mp_waitall (col_rr(0:ncol_rr-1))
       END IF
       IF (.NOT.remainder_collective) THEN
          CALL mp_waitall (all_sr(1:nall_sr))
          CALL mp_waitall (all_rr(1:nall_rr))
       ENDIF
    ELSE
       CALL mp_alltoall (sb, scount, sdispl,&
            rb, rcount, rdispl,&
            all_group)
    ENDIF subgrouped
  CONTAINS
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE most_alltoall()
      FORALL (pcol = 0 : npcols-1)
         new_scount(1+pcol) = scount(1+pgrid(myprow, pcol))
         new_rcount(1+pcol) = rcount(1+pgrid(myprow, pcol))
         new_sdispl(1+pcol) = sdispl(1+pgrid(myprow, pcol))
         new_rdispl(1+pcol) = rdispl(1+pgrid(myprow, pcol))
      END FORALL
      CALL mp_alltoall (sb, new_scount(1:npcols), new_sdispl(1:npcols),&
           rb, new_rcount(1:npcols), new_rdispl(1:npcols),&
           dbcsr_mp_my_row_group (mp_env))
      FORALL (prow = 0 : nprows-1)
         new_scount(1+prow) = scount(1+pgrid(prow, mypcol))
         new_rcount(1+prow) = rcount(1+pgrid(prow, mypcol))
         new_sdispl(1+prow) = sdispl(1+pgrid(prow, mypcol))
         new_rdispl(1+prow) = rdispl(1+pgrid(prow, mypcol))
      END FORALL
      CALL mp_alltoall (sb, new_scount(1:nprows), new_sdispl(1:nprows),&
           rb, new_rcount(1:nprows), new_rdispl(1:nprows),&
           dbcsr_mp_my_col_group (mp_env))
    END SUBROUTINE most_alltoall
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE most_point_to_point ()
      ! Go through my prow and exchange.
      DO i = 0, npcols - 1
         pcol = MOD (mypcol+i, npcols)
         grp = dbcsr_mp_my_row_group (mp_env)
         !
         dst = dbcsr_mp_get_process (mp_env, myprow, pcol)
         send_cnt = scount(dst+1)
         IF (send_cnt .GT. 0) THEN
            send_data_p => sb( 1+sdispl(dst+1) : 1+sdispl(dst+1)+send_cnt-1 )
            tag = 4*mypcol
            CALL mp_isend (send_data_p, pcol, grp, row_sr(nrow_sr), tag)
            nrow_sr = nrow_sr+1
         ENDIF
         !
         pcol = MODULO (mypcol-i, npcols)
         src = dbcsr_mp_get_process (mp_env, myprow, pcol)
         recv_cnt = rcount(src+1)
         IF (recv_cnt .GT. 0) THEN
            recv_data_p => rb( 1+rdispl(src+1) : 1+rdispl(src+1)+recv_cnt-1 )
            tag = 4*pcol
            CALL mp_irecv (recv_data_p, pcol, grp, row_rr(nrow_rr), tag)
            nrow_rr = nrow_rr+1
         ENDIF
      ENDDO
      ! go through my pcol and exchange
      DO i = 0, nprows - 1
         prow = MOD (myprow+i, nprows)
         grp = dbcsr_mp_my_col_group (mp_env)
         !
         dst = dbcsr_mp_get_process (mp_env, prow, mypcol)
         send_cnt = scount(dst+1)
         IF (send_cnt .GT. 0) THEN
            send_data_p => sb( 1+sdispl(dst+1) : 1+sdispl(dst+1)+send_cnt-1 )
            tag = 4*myprow+1
            CALL mp_isend (send_data_p, prow, grp, col_sr(ncol_sr), tag)
            ncol_sr = ncol_sr + 1
         ENDIF
         !
         prow = MODULO (myprow-i, nprows)
         src = dbcsr_mp_get_process (mp_env, prow, mypcol)
         recv_cnt = rcount(src+1)
         IF (recv_cnt .GT. 0) THEN
            recv_data_p => rb( 1+rdispl(src+1) : 1+rdispl(src+1)+recv_cnt-1 )
            tag = 4*prow+1
            CALL mp_irecv (recv_data_p, prow, grp, col_rr(ncol_rr), tag)
            ncol_rr = ncol_rr + 1
         ENDIF
      ENDDO
    END SUBROUTINE most_point_to_point
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE remainder_alltoall ()
      new_scount(:) = scount(:)
      new_rcount(:) = rcount(:)
      FORALL (prow = 0:nprows-1)
         new_scount(1+pgrid(prow, mypcol)) = 0
         new_rcount(1+pgrid(prow, mypcol)) = 0
      END FORALL
      FORALL (pcol = 0:npcols-1)
         new_scount(1+pgrid(myprow, pcol)) = 0
         new_rcount(1+pgrid(myprow, pcol)) = 0
      END FORALL
      CALL mp_alltoall (sb, new_scount, sdispl,&
           rb, new_rcount, rdispl, all_group)
    END SUBROUTINE remainder_alltoall
! *****************************************************************************
!> \brief ...
! *****************************************************************************
    SUBROUTINE remainder_point_to_point()
    INTEGER                                  :: col, row

      DO row = 0, nprows-1
         prow = MOD(row+myprow, nprows)
         IF (prow .EQ. myprow) CYCLE
         DO col = 0, npcols-1
            pcol = MOD (col+mypcol, npcols)
            IF (pcol .EQ. mypcol) CYCLE
            dst = dbcsr_mp_get_process (mp_env, prow, pcol)
            send_cnt = scount(dst+1)
            IF (send_cnt .GT. 0) THEN
               send_data_p => sb( 1+sdispl(dst+1) : 1+sdispl(dst+1)+send_cnt-1 )
               tag = 4*mynode+2
               CALL mp_isend (send_data_p, dst, all_group, all_sr(nall_sr+1), tag)
               nall_sr = nall_sr + 1
            ENDIF
            !
            src = dbcsr_mp_get_process (mp_env, prow, pcol)
            recv_cnt = rcount(src+1)
            IF (recv_cnt .GT. 0) THEN
               recv_data_p => rb( 1+rdispl(src+1) : 1+rdispl(src+1)+recv_cnt-1 )
               tag = 4*src+2
               CALL mp_irecv (recv_data_p, src, all_group, all_rr(nall_rr+1), tag)
               nall_rr = nall_rr+1
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE remainder_point_to_point
  END SUBROUTINE hybrid_alltoall_z1
