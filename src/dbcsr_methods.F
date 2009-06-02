!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Base methods on DBCSR data structures
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.95
!>
!> <b>Modification history:</b>
!> - Created
! *****************************************************************************
MODULE dbcsr_methods
  USE array_types,                     ONLY: array_hold,&
                                             array_i1d_obj,&
                                             array_release,&
                                             array_size,&
                                             array_data,&
                                             array_new
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_mp_obj
  USE distribution_2d_types,           ONLY: distribution_2d_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_methods'

  PUBLIC :: dbcsr_mp_new, dbcsr_mp_hold, dbcsr_mp_release
  PUBLIC :: dbcsr_mp_pgrid, dbcsr_mp_numnodes, dbcsr_mp_mynode, dbcsr_mp_group,&
            dbcsr_mp_new_transposed, dbcsr_mp_nprows, dbcsr_mp_npcols,&
            dbcsr_mp_myprow, dbcsr_mp_mypcol
  PUBLIC :: dbcsr_distribution_new, dbcsr_distribution_hold,&
            dbcsr_distribution_release, dbcsr_distribution_init
  PUBLIC :: dbcsr_distribution_mp, dbcsr_distribution_processor,&
            dbcsr_distribution_nrows, dbcsr_distribution_ncols,&
            dbcsr_distribution_row_dist, dbcsr_distribution_col_dist,&
            dbcsr_distribution_nlocal_rows, dbcsr_distribution_nlocal_cols,&
            dbcsr_distribution_local_rows, dbcsr_distribution_local_cols

  

  PRIVATE

#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)


CONTAINS

! *****************************************************************************
!> \brief Creates new process grid
!> \param[out] mp_env         multiprocessor environment
!> \param[in] pgrid           process grid
!> \param[in] mynode          my processor number
!> \param[in] numnodes         total number of processors (processes)
! *****************************************************************************
  SUBROUTINE dbcsr_mp_new(mp_env, pgrid, mp_group, mynode, numnodes, myprow,&
       mypcol)
    TYPE(dbcsr_mp_obj), INTENT(OUT)          :: mp_env
    INTEGER, DIMENSION(0:, 0:), INTENT(IN)   :: pgrid
    INTEGER, INTENT(IN)                      :: mp_group, mynode
    INTEGER, INTENT(IN), OPTIONAL            :: numnodes, myprow, mypcol

    INTEGER                                  :: pcol, prow

!   ---------------------------------------------------------------------------

    ALLOCATE(mp_env%mp)
    mp_env%mp%refcount = 1
    ALLOCATE (mp_env%mp%pgrid (0:SIZE(pgrid, 1)-1, 0:SIZE(pgrid, 2)-1 ))
    mp_env%mp%pgrid(:,:) = pgrid(:,:)
    mp_env%mp%mynode = mynode
    mp_env%mp%mp_group = mp_group
    IF (PRESENT (numnodes)) THEN
       mp_env%mp%numnodes = numnodes
    ELSE
       mp_env%mp%numnodes = SIZE (pgrid)
    ENDIF
    IF (PRESENT (myprow) .AND. PRESENT (mypcol)) THEN
       mp_env%mp%myprow = myprow
       mp_env%mp%mypcol = mypcol
    ELSE
       mp_env%mp%myprow = -33777
       mp_env%mp%mypcol = -33777
       column_loop: DO pcol = LBOUND (pgrid, 2), UBOUND (pgrid, 2)
          row_loop: DO prow = LBOUND (pgrid, 1), UBOUND (pgrid, 1)
             test_position: IF (pgrid (prow, pcol) .EQ. mynode) THEN
                mp_env%mp%myprow = prow
                mp_env%mp%mypcol = pcol
                EXIT column_loop
             ENDIF test_position
          ENDDO row_loop
       ENDDO column_loop
    ENDIF
  END SUBROUTINE dbcsr_mp_new

! *****************************************************************************
!> \brief Marks another use of the mp_env
!> \param[in,out] mp_env      multiprocessor environment
! *****************************************************************************
  PURE SUBROUTINE dbcsr_mp_hold(mp_env)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env

!   ---------------------------------------------------------------------------

    mp_env%mp%refcount = mp_env%mp%refcount+1
  END SUBROUTINE dbcsr_mp_hold

! *****************************************************************************
!> \brief Releases and potentially destrops an mp_env
!> \param[in,out] mp_env         multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_mp_release(mp_env)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (mp_env%mp)) THEN
       mp_env%mp%refcount = mp_env%mp%refcount - 1
       IF (mp_env%mp%refcount .LE. 0) THEN
          DEALLOCATE (mp_env%mp%pgrid)
          DEALLOCATE (mp_env%mp)
          NULLIFY (mp_env%mp)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_mp_release

  FUNCTION dbcsr_mp_pgrid(mp_env) RESULT (pgrid)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    pgrid => mp_env%mp%pgrid
  END FUNCTION dbcsr_mp_pgrid
  PURE FUNCTION dbcsr_mp_numnodes(mp_env) RESULT (numnodes)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: numnodes

    numnodes = mp_env%mp%numnodes
  END FUNCTION dbcsr_mp_numnodes
  PURE FUNCTION dbcsr_mp_mynode(mp_env) RESULT (mynode)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mynode

    mynode = mp_env%mp%mynode
  END FUNCTION dbcsr_mp_mynode
  PURE FUNCTION dbcsr_mp_group(mp_env) RESULT (mp_group)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mp_group

    mp_group = mp_env%mp%mp_group
  END FUNCTION dbcsr_mp_group
  PURE FUNCTION dbcsr_mp_nprows(mp_env) RESULT (nprows)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: nprows

    nprows = SIZE (mp_env%mp%pgrid, 1)
  END FUNCTION dbcsr_mp_nprows
  PURE FUNCTION dbcsr_mp_npcols(mp_env) RESULT (npcols)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: npcols

    npcols = SIZE (mp_env%mp%pgrid, 2)
  END FUNCTION dbcsr_mp_npcols
  PURE FUNCTION dbcsr_mp_myprow(mp_env) RESULT (myprow)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: myprow

    myprow = mp_env%mp%myprow
  END FUNCTION dbcsr_mp_myprow
  PURE FUNCTION dbcsr_mp_mypcol(mp_env) RESULT (mypcol)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mypcol

    mypcol = mp_env%mp%mypcol
  END FUNCTION dbcsr_mp_mypcol


! *****************************************************************************
!> \brief Transposes a multiprocessor environment
!> \param[out] mp_t           transposed multiprocessor environment
!> \param[in] mp              original multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_mp_new_transposed(mp_t, mp)
    TYPE(dbcsr_mp_obj), INTENT(OUT)          :: mp_t
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp

!   ---------------------------------------------------------------------------

    CALL dbcsr_mp_new (mp_t, TRANSPOSE (dbcsr_mp_pgrid (mp)),&
         dbcsr_mp_group (mp),&
         dbcsr_mp_mynode (mp), dbcsr_mp_numnodes (mp),&
         dbcsr_mp_mypcol (mp), dbcsr_mp_myprow (mp))
  END SUBROUTINE dbcsr_mp_new_transposed


! *****************************************************************************
!> \brief Creates new distribution
!> \param[out] dist           distribution
!> \param[in] mp_env          multiprocessing environment
!> \param[in] row_dist, col_dist        row and column distributions
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_new(dist, mp_env, row_dist, col_dist,&
       local_rows, local_cols)
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    TYPE(array_i1d_obj), INTENT(IN)          :: row_dist, col_dist
    TYPE(array_i1d_obj), INTENT(IN), OPTIONAL     :: local_rows, local_cols
    TYPE(distribution_2d_type), POINTER      :: dist2d
    INTEGER                                  :: i, seq, mypcoor
    INTEGER, DIMENSION(:), POINTER           :: dd, ld

!   ---------------------------------------------------------------------------

    ALLOCATE (dist%d)
    dist%d%refcount = 1
    dist%d%row_dist = row_dist
    CALL array_hold (dist%d%row_dist)
    dist%d%col_dist = col_dist
    CALL array_hold (dist%d%col_dist)
    dist%d%mp_env = mp_env
    CALL dbcsr_mp_hold (dist%d%mp_env)
    IF (PRESENT (local_rows)) THEN
       dist%d%local_rows = local_rows
       CALL array_hold (dist%d%local_rows)
    ELSE
       dd => array_data (row_dist)
       mypcoor = dbcsr_mp_myprow (mp_env)
       i = COUNT (dd.EQ.mypcoor)
       NULLIFY (ld)
       ALLOCATE (ld(i))
       seq = 1
       DO i = 1, array_size (row_dist)
          IF (dd(i) .EQ. mypcoor) THEN
             ld(seq) = i
             seq = seq+1
          ENDIF
       ENDDO
       CALL array_new (dist%d%local_rows, ld, gift=.TRUE.)
    ENDIF
    IF (PRESENT (local_cols)) THEN
       dist%d%local_cols = local_cols
       CALL array_hold (dist%d%local_cols)
    ELSE
       dd => array_data (col_dist)
       mypcoor = dbcsr_mp_mypcol (mp_env)
       i = COUNT (dd.EQ.mypcoor)
       NULLIFY (ld)
       ALLOCATE (ld(i))
       seq = 1
       DO i = 1, array_size (col_dist)
          IF (dd(i) .EQ. mypcoor) THEN
             ld(seq) = i
             seq = seq+1
          ENDIF
       ENDDO
       CALL array_new (dist%d%local_cols, ld, gift=.TRUE.)
    ENDIF
    !WRITE(*,*)' row_dist', array_data(dist%d%row_dist)
    !WRITE(*,*)' col_dist', array_data(dist%d%col_dist)
    !WRITE(*,*)' loc_rows', array_data(dist%d%local_rows)
    !WRITE(*,*)' loc_cols', array_data(dist%d%local_cols)
  END SUBROUTINE dbcsr_distribution_new

! *****************************************************************************
!> \brief Marks another use of the distribution
!> \param[in,out] mp_env      multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_hold(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

!   ---------------------------------------------------------------------------

    dist%d%refcount = dist%d%refcount + 1
  END SUBROUTINE dbcsr_distribution_hold

! *****************************************************************************
!> \brief Releases and potentially destrops a distribution
!> \param[in,out] mp_env         multiprocessor environment
! *****************************************************************************
  SUBROUTINE dbcsr_distribution_release(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT)                          :: dist

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (dist%d)) THEN
       dist%d%refcount = dist%d%refcount - 1
       IF (dist%d%refcount .EQ. 0) THEN
          CALL array_release (dist%d%row_dist)
          CALL array_release (dist%d%col_dist)
          CALL array_release (dist%d%local_rows)
          CALL array_release (dist%d%local_cols)
          CALL dbcsr_mp_release (dist%d%mp_env)
          DEALLOCATE (dist%d)
          CALL dbcsr_distribution_init (dist)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_distribution_release

  SUBROUTINE dbcsr_distribution_init (dist)
    TYPE (dbcsr_distribution_obj), INTENT(OUT) :: dist
    NULLIFY (dist%d)
  END SUBROUTINE dbcsr_distribution_init

  FUNCTION dbcsr_distribution_mp(dist) RESULT (mp_env)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(dbcsr_mp_obj)                       :: mp_env

!   ---------------------------------------------------------------------------

    mp_env = dist%d%mp_env
  END FUNCTION dbcsr_distribution_mp
  PURE FUNCTION dbcsr_distribution_nrows(dist) RESULT (nrows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nrows
    nrows = array_size (dist%d%row_dist)
  END FUNCTION dbcsr_distribution_nrows
  PURE FUNCTION dbcsr_distribution_ncols(dist) RESULT (ncols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: ncols
    ncols = array_size (dist%d%col_dist)
  END FUNCTION dbcsr_distribution_ncols
  FUNCTION dbcsr_distribution_row_dist(dist) RESULT (row_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: row_dist

!   ---------------------------------------------------------------------------

    row_dist = dist%d%row_dist
  END FUNCTION dbcsr_distribution_row_dist
  FUNCTION dbcsr_distribution_col_dist(dist) RESULT (col_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: col_dist

!   ---------------------------------------------------------------------------

    col_dist = dist%d%col_dist
  END FUNCTION dbcsr_distribution_col_dist

  PURE FUNCTION dbcsr_distribution_nlocal_rows(dist) RESULT (nlocalrows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nlocalrows
    nlocalrows = array_size (dist%d%local_rows)
  END FUNCTION dbcsr_distribution_nlocal_rows
  PURE FUNCTION dbcsr_distribution_nlocal_cols(dist) RESULT (nlocalcols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER                                  :: nlocalcols
    nlocalcols = array_size (dist%d%local_cols)
  END FUNCTION dbcsr_distribution_nlocal_cols
  FUNCTION dbcsr_distribution_local_rows(dist) RESULT (local_rows)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: local_rows
    local_rows = dist%d%local_rows
  END FUNCTION dbcsr_distribution_local_rows
  FUNCTION dbcsr_distribution_local_cols(dist) RESULT (local_cols)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: local_cols
    local_cols = dist%d%local_cols
  END FUNCTION dbcsr_distribution_local_cols
  !
  PURE FUNCTION dbcsr_distribution_processor(dist, row, col)&
       RESULT (processor)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    INTEGER, INTENT(IN)                      :: row, col
    INTEGER                                  :: processor
    processor = dist%d%mp_env%mp%pgrid(dist%d%row_dist%low%data(row),&
         dist%d%col_dist%low%data(col))
  END FUNCTION dbcsr_distribution_processor
END MODULE dbcsr_methods
