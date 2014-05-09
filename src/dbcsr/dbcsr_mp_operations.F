!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Wrappers to message passing calls.
!> \author  Urban Borstnik
!> \date    2010-02-18
!> \version 0.9
!>
!> <b>Modification history:</b>
!> - 2010-02-18 Moved from dbcsr_util
! *****************************************************************************
MODULE dbcsr_mp_operations

  USE dbcsr_config,                    ONLY: has_MPI
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_type
  USE dbcsr_error_handling
  USE dbcsr_methods,                   ONLY: &
       dbcsr_mp_get_process, dbcsr_mp_grid_setup, dbcsr_mp_group, &
       dbcsr_mp_has_subgroups, dbcsr_mp_my_col_group, dbcsr_mp_my_row_group, &
       dbcsr_mp_mynode, dbcsr_mp_mypcol, dbcsr_mp_myprow, dbcsr_mp_npcols, &
       dbcsr_mp_nprows, dbcsr_mp_numnodes, dbcsr_mp_pgrid
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_mp_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE message_passing,                 ONLY: &
       mp_allgather, mp_alltoall, mp_irecv, mp_isend, mp_recv, mp_send, &
       mp_type_descriptor_type, mp_type_make, mp_waitall

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mp_operations'

  ! MP routines
  PUBLIC :: dbcsr_hybrid_alltoall
  PUBLIC :: hybrid_alltoall_s1, hybrid_alltoall_d1,&
            hybrid_alltoall_c1, hybrid_alltoall_z1,&
            hybrid_alltoall_i1, hybrid_alltoall_any
  PUBLIC :: dbcsr_allgatherv
  PUBLIC :: dbcsr_send_any, dbcsr_recv_any
  PUBLIC :: dbcsr_isend_any, dbcsr_irecv_any
  ! Type helpers
  PUBLIC :: dbcsr_mp_type_from_anytype

  INTERFACE dbcsr_hybrid_alltoall
     MODULE PROCEDURE hybrid_alltoall_s1, hybrid_alltoall_d1,&
                      hybrid_alltoall_c1, hybrid_alltoall_z1
     MODULE PROCEDURE hybrid_alltoall_i1
     MODULE PROCEDURE hybrid_alltoall_any
  END INTERFACE

CONTAINS


  SUBROUTINE hybrid_alltoall_any (sb, scount, sdispl,&
       rb, rcount, rdispl, mp_env, most_ptp, remainder_ptp, no_hybrid,&
       error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL, INTENT(in), OPTIONAL            :: most_ptp, remainder_ptp, &
                                                no_hybrid
    TYPE(dbcsr_error_type)                   :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hybrid_alltoall_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)

    SELECT CASE (dbcsr_data_get_type (sb))
    CASE (dbcsr_type_real_4)
       CALL hybrid_alltoall_s1 (sb%d%r_sp, scount, sdispl,&
       rb%d%r_sp, rcount, rdispl, mp_env,&
       most_ptp, remainder_ptp, no_hybrid)
    CASE (dbcsr_type_real_8)
       CALL hybrid_alltoall_d1 (sb%d%r_dp, scount, sdispl,&
       rb%d%r_dp, rcount, rdispl, mp_env,&
       most_ptp, remainder_ptp, no_hybrid)
    CASE (dbcsr_type_complex_4)
       CALL hybrid_alltoall_c1 (sb%d%c_sp, scount, sdispl,&
       rb%d%c_sp, rcount, rdispl, mp_env,&
       most_ptp, remainder_ptp, no_hybrid)
    CASE (dbcsr_type_complex_8)
       CALL hybrid_alltoall_z1 (sb%d%c_dp, scount, sdispl,&
       rb%d%c_dp, rcount, rdispl, mp_env,&
       most_ptp, remainder_ptp, no_hybrid)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_internal_error,&
            routineN, "Invalid data type",__LINE__,error)
    END SELECT

    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE hybrid_alltoall_any


! *****************************************************************************
!> \brief Row/column and global all-to-all
!> \par Communicator selection
!>      Uses row and column communicators for row/column
!>      sends. Remaining sends are performed using the global
!>      communicator.  Point-to-point isend/irecv are used if ptp is
!>      set, otherwise a alltoall collective call is issued.
!> \param[in] mp_env         MP Environment
!> \param[in] most_ptp       (optional) Use point-to-point for row/column;
!>                           default is no
!> \param[in] remainder_ptp  (optional) Use point-to-point for remaining;
!>                           default is no
!> \param[in] no_hybrid      (optional) Use regular global collective; default
!>                           is no
!> \see mp_alltoall
! *****************************************************************************
  SUBROUTINE hybrid_alltoall_i1 (sb, scount, sdispl,&
       rb, rcount, rdispl, mp_env, most_ptp, remainder_ptp, no_hybrid)
    INTEGER, DIMENSION(:), INTENT(in), &
      TARGET                                 :: sb
    INTEGER, DIMENSION(:), INTENT(IN)        :: scount, sdispl
    INTEGER, DIMENSION(:), &
      INTENT(INOUT), TARGET                  :: rb
    INTEGER, DIMENSION(:), INTENT(IN)        :: rcount, rdispl
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL, INTENT(IN), OPTIONAL            :: most_ptp, remainder_ptp,&
                                                no_hybrid

    CHARACTER(len=*), PARAMETER :: routineN = 'hybrid_alltoall_i1', &
      routineP = moduleN//':'//routineN

    INTEGER :: all_group, mynode, mypcol, myprow, nall_rr, nall_sr, ncol_rr, &
      ncol_sr, npcols, nprows, nrow_rr, nrow_sr, numnodes, dst, src,&
      prow, pcol, send_cnt, recv_cnt, tag, grp, i
    INTEGER, ALLOCATABLE, DIMENSION(:) :: all_rr, all_sr, col_rr, col_sr, &
      new_rcount, new_rdispl, new_scount, new_sdispl, row_rr, row_sr
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid
    LOGICAL                                  :: most_collective, &
                                                remainder_collective, no_h
    INTEGER, DIMENSION(:), POINTER           :: send_data_p, recv_data_p
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
    SUBROUTINE most_point_to_point ()
      ! Go through my prow and exchange.
      DO i = 0, npcols - 1
         pcol = MOD (mypcol+i, npcols)
         grp = dbcsr_mp_my_row_group (mp_env)
         !
         dst = dbcsr_mp_get_process (mp_env, myprow, pcol)
         send_cnt = scount(dst+1)
         send_data_p => sb( 1+sdispl(dst+1) : 1+sdispl(dst+1)+send_cnt-1 )
         tag = 4*mypcol
         IF (send_cnt .GT. 0) THEN
            CALL mp_isend (send_data_p, pcol, grp, row_sr(nrow_sr), tag)
            nrow_sr = nrow_sr+1
         ENDIF
         !
         pcol = MODULO (mypcol-i, npcols)
         src = dbcsr_mp_get_process (mp_env, myprow, pcol)
         recv_cnt = rcount(src+1)
         recv_data_p => rb( 1+rdispl(src+1) : 1+rdispl(src+1)+recv_cnt-1 )
         tag = 4*pcol
         IF (recv_cnt .GT. 0) THEN
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
               tag = 4*mynode+2
               send_data_p => sb( 1+sdispl(dst+1) : 1+sdispl(dst+1)+send_cnt-1 )
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
  END SUBROUTINE hybrid_alltoall_i1

! *****************************************************************************
!> \brief Creates an MPI combined type from the given anytype.
!> \param[in] data_area       Data area of any type
!> \result mp_type            Type descriptor
! *****************************************************************************
  FUNCTION dbcsr_mp_type_from_anytype(data_area) RESULT (mp_type)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: data_area
    TYPE(mp_type_descriptor_type)            :: mp_type

    SELECT CASE (data_area%d%data_type)
    CASE (dbcsr_type_real_4)
       mp_type = mp_type_make (data_area%d%r_sp)
    CASE (dbcsr_type_real_8)
       mp_type = mp_type_make (data_area%d%r_dp)
    CASE (dbcsr_type_complex_4)
       mp_type = mp_type_make (data_area%d%c_sp)
    CASE (dbcsr_type_complex_8)
       mp_type = mp_type_make (data_area%d%c_dp)
    END SELECT
  END FUNCTION dbcsr_mp_type_from_anytype


! *****************************************************************************
!> \brief Blocking send of encapsulated data.
!> \sa mp_send_iv
! *****************************************************************************
  SUBROUTINE dbcsr_send_any(msgin,dest,tag,comm, error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: msgin
    INTEGER, INTENT(IN)                      :: dest, tag, comm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_send_any', &
      routineP = moduleN//':'//routineN

    SELECT CASE (dbcsr_data_get_type (msgin))
    CASE (dbcsr_type_real_4)
       CALL mp_send (msgin%d%r_sp, dest, tag, comm)
    CASE (dbcsr_type_real_8)
       CALL mp_send (msgin%d%r_dp, dest, tag, comm)
    CASE (dbcsr_type_complex_4)
       CALL mp_send (msgin%d%c_sp, dest, tag, comm)
    CASE (dbcsr_type_complex_8)
       CALL mp_send (msgin%d%c_dp, dest, tag, comm)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Incorrect data type", __LINE__, error=error)
    END SELECT
  END SUBROUTINE dbcsr_send_any

! *****************************************************************************
!> \brief Non-blocking send of encapsulated data.
!> \sa mp_isend_iv
! *****************************************************************************
  SUBROUTINE dbcsr_isend_any(msgin,dest,comm,request,tag, error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: msgin
    INTEGER, INTENT(IN)                      :: dest, comm
    INTEGER, INTENT(OUT)                     :: request
    INTEGER, INTENT(IN), OPTIONAL            :: tag
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_isend_any', &
      routineP = moduleN//':'//routineN

    SELECT CASE (dbcsr_data_get_type (msgin))
    CASE (dbcsr_type_real_4)
       CALL mp_isend (msgin%d%r_sp, dest, comm, request, tag)
    CASE (dbcsr_type_real_8)
       CALL mp_isend (msgin%d%r_dp, dest, comm, request, tag)
    CASE (dbcsr_type_complex_4)
       CALL mp_isend (msgin%d%c_sp, dest, comm, request, tag)
    CASE (dbcsr_type_complex_8)
       CALL mp_isend (msgin%d%c_dp, dest, comm, request, tag)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Incorrect data type", __LINE__, error=error)
    END SELECT
  END SUBROUTINE dbcsr_isend_any

! *****************************************************************************
!> \brief Blocking recv of encapsulated data.
!> \sa mp_recv_iv
! *****************************************************************************
  SUBROUTINE dbcsr_recv_any(msgin,source,tag, comm, error)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: msgin
    INTEGER, INTENT(IN)                      :: source, tag, comm
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_recv_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: io_comm, io_source, io_tag

    io_source = source ; io_comm = comm; io_tag = tag

    SELECT CASE (dbcsr_data_get_type (msgin))
    CASE (dbcsr_type_real_4)
       CALL mp_recv (msgin%d%r_sp, io_source, io_tag, io_comm)
    CASE (dbcsr_type_real_8)
       CALL mp_recv (msgin%d%r_dp, io_source, io_tag, io_comm)
    CASE (dbcsr_type_complex_4)
       CALL mp_recv (msgin%d%c_sp, io_source, io_tag, io_comm)
    CASE (dbcsr_type_complex_8)
       CALL mp_recv (msgin%d%c_dp, io_source, io_tag, io_comm)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Incorrect data type", __LINE__, error=error)
    END SELECT
  END SUBROUTINE dbcsr_recv_any

! *****************************************************************************
!> \brief Non-blocking recv of encapsulated data.
!> \sa mp_irecv_iv
! *****************************************************************************
  SUBROUTINE dbcsr_irecv_any(msgin,source,comm,request,tag, error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: msgin
    INTEGER, INTENT(IN)                      :: source, comm
    INTEGER, INTENT(OUT)                     :: request
    INTEGER, INTENT(IN), OPTIONAL            :: tag
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_irecv_any', &
      routineP = moduleN//':'//routineN

    SELECT CASE (dbcsr_data_get_type (msgin))
    CASE (dbcsr_type_real_4)
       CALL mp_irecv (msgin%d%r_sp, source, comm, request, tag)
    CASE (dbcsr_type_real_8)
       CALL mp_irecv (msgin%d%r_dp, source, comm, request, tag)
    CASE (dbcsr_type_complex_4)
       CALL mp_irecv (msgin%d%c_sp, source, comm, request, tag)
    CASE (dbcsr_type_complex_8)
       CALL mp_irecv (msgin%d%c_dp, source, comm, request, tag)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_caller_error,&
            routineN, "Incorrect data type", __LINE__, error=error)
    END SELECT
  END SUBROUTINE dbcsr_irecv_any


! *****************************************************************************
!> \brief Allgather of encapsulated data
!> \sa mp_allgatherv_dv
! *****************************************************************************
  SUBROUTINE dbcsr_allgatherv(send_data, recv_data, recv_count, recv_displ, gid)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: send_data
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: recv_data
    INTEGER, DIMENSION(:), INTENT(IN)        :: recv_count, recv_displ
    INTEGER, INTENT(IN)                      :: gid

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_allgatherv', &
      routineP = moduleN//':'//routineN

    TYPE(dbcsr_error_type)                   :: error

    CALL dbcsr_assert (dbcsr_data_get_type (send_data), "EQ",&
         dbcsr_data_get_type (recv_data), dbcsr_fatal_level,&
         dbcsr_wrong_args_error, routineN, "Data type mismatch", __LINE__,&
         error=error)
    SELECT CASE (dbcsr_data_get_type(send_data))
    CASE (dbcsr_type_real_4)
       CALL mp_allgather(send_data%d%r_sp, recv_data%d%r_sp,&
       recv_count, recv_displ, gid)
    CASE (dbcsr_type_real_8)
       CALL mp_allgather(send_data%d%r_dp, recv_data%d%r_dp,&
       recv_count, recv_displ, gid)
    CASE (dbcsr_type_complex_4)
       CALL mp_allgather(send_data%d%c_sp, recv_data%d%c_sp,&
       recv_count, recv_displ, gid)
    CASE (dbcsr_type_complex_8)
       CALL mp_allgather(send_data%d%c_dp, recv_data%d%c_dp,&
       recv_count, recv_displ, gid)
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            routineN, "Invalid data type", __LINE__, error=error)
    END SELECT
  END SUBROUTINE dbcsr_allgatherv


#include "dbcsr_mp_operations_d.f90"
#include "dbcsr_mp_operations_z.f90"
#include "dbcsr_mp_operations_s.f90"
#include "dbcsr_mp_operations_c.f90"

END MODULE dbcsr_mp_operations
