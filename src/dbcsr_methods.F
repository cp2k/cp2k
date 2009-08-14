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
  USE array_types,                     ONLY: array_data,&
                                             array_hold,&
                                             array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release,&
                                             array_size
  USE btree_I8_k_cp2d_v,               ONLY: btree_destroy_c => btree_delete,&
                                             btree_new_c => btree_new
  USE btree_I8_k_dp2d_v,               ONLY: btree_destroy_d => btree_delete,&
                                             btree_new_d => btree_new
  USE btree_I8_k_sp2d_v,               ONLY: btree_destroy_s => btree_delete,&
                                             btree_new_s => btree_new
  USE btree_I8_k_zp2d_v,               ONLY: btree_destroy_z => btree_delete,&
                                             btree_new_z => btree_new
  USE dbcsr_types,                     ONLY: &
       dbcsr_1d_set_obj, dbcsr_2d_set_obj, dbcsr_data_obj, &
       dbcsr_distribution_obj, dbcsr_iterator, dbcsr_mp_obj, &
       dbcsr_mutable_obj, dbcsr_obj, dbcsr_type, dbcsr_type_complex_4, &
       dbcsr_type_complex_8, dbcsr_type_real_4, dbcsr_type_real_8
  USE distribution_2d_types,           ONLY: distribution_2d_release
  USE f77_blas
  USE f77_blas_generic,                ONLY: bl_copy
  USE kinds,                           ONLY: default_string_length,&
                                             real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_allocate,&
                                             mp_deallocate

  !$ USE OMP_LIB
#include "cp_common_uses.h"
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_methods'

  PUBLIC :: dbcsr_hold, dbcsr_release, dbcsr_destroy
  PUBLIC :: dbcsr_mp_new, dbcsr_mp_hold, dbcsr_mp_release
  PUBLIC :: dbcsr_mp_pgrid, dbcsr_mp_numnodes, dbcsr_mp_mynode, dbcsr_mp_group,&
            dbcsr_mp_new_transposed, dbcsr_mp_nprows, dbcsr_mp_npcols,&
            dbcsr_mp_myprow, dbcsr_mp_mypcol, dbcsr_mp_pgrid_equal
  PUBLIC :: dbcsr_distribution_new, dbcsr_distribution_hold,&
            dbcsr_distribution_release, dbcsr_distribution_init
  PUBLIC :: dbcsr_distribution_mp, dbcsr_distribution_processor,&
            dbcsr_distribution_nrows, dbcsr_distribution_ncols,&
            dbcsr_distribution_row_dist, dbcsr_distribution_col_dist,&
            dbcsr_distribution_nlocal_rows, dbcsr_distribution_nlocal_cols,&
            dbcsr_distribution_local_rows, dbcsr_distribution_local_cols,&
            dbcsr_distribution_thread_dist,&
            dbcsr_distribution_make_threads
  PUBLIC :: dbcsr_get_info, dbcsr_distribution,&
            dbcsr_get_matrix_type, dbcsr_get_data_type,&
            dbcsr_row_block_sizes, dbcsr_col_block_sizes,&
            dbcsr_nblkrows_total, dbcsr_nblkcols_total, dbcsr_nfullrows_total,&
            dbcsr_nfullcols_total, dbcsr_nblkcols_local, dbcsr_nblkrows_local,&
            dbcsr_nfullrows_local, dbcsr_nfullcols_local,&
            dbcsr_max_row_size, dbcsr_max_col_size,&
            dbcsr_uses_special_memory
  PUBLIC :: dbcsr_data_init, dbcsr_get_data, dbcsr_data_new, dbcsr_data_hold,&
            dbcsr_data_release, dbcsr_data_copyall, dbcsr_get_data_size,&
            dbcsr_data_set_pointer, dbcsr_data_area
  PUBLIC :: dbcsr_get_num_blocks
  PUBLIC :: dbcsr_get_data_p, dbcsr_get_data_p_s, dbcsr_get_data_p_c,&
            dbcsr_get_data_p_d, dbcsr_get_data_p_z

  PUBLIC :: dbcsr_iterator_start, dbcsr_iterator_next_block,&
            dbcsr_iterator_stop, dbcsr_iterator_blocks_left

  PUBLIC :: dbcsr_allocate_matrix_set, dbcsr_set_new, dbcsr_set_put, dbcsr_set_get,&
            dbcsr_set_destroy, dbcsr_set_hold, dbcsr_set_release

  PUBLIC :: dbcsr_mutable_init, dbcsr_mutable_new, dbcsr_mutable_destroy,&
            dbcsr_mutable_release, dbcsr_mutable_hold,&
            dbcsr_mutable_instantiated

  PRIVATE

#define DBG IF (dbg) WRITE(*,*)routineN//" ",
#define DBGV IF (bcsr_verbose) WRITE(*,*)routineN//" ",
#define DBGI IF (info) WRITE(*,*)routineN//" ",
#define DEBUG_HEADER        LOGICAL :: dbg, info
#define DEBUG_BODY        dbg = .FALSE. ; CALL dbcsr_set_debug(dbg, info=info)
#define DEBUG_BODYY       dbg = .TRUE. ; CALL dbcsr_set_debug(dbg, info=info)

  INTERFACE dbcsr_get_data
     MODULE PROCEDURE get_data_s, get_data_d, get_data_c, get_data_z,&
          get_data_m_s, get_data_m_d, get_data_m_c, get_data_m_z
  END INTERFACE

  INTERFACE dbcsr_data_set_pointer
     MODULE PROCEDURE set_data_p_s, set_data_p_d, set_data_p_c, set_data_p_z
  END INTERFACE

  INTERFACE dbcsr_get_data_p
     MODULE PROCEDURE dbcsr_get_data_c_s, dbcsr_get_data_c_c,&
          dbcsr_get_data_c_d, dbcsr_get_data_c_z
  END INTERFACE

  INTERFACE dbcsr_iterator_next_block
     MODULE PROCEDURE iterator_next_block_d, iterator_next_block_index,&
                      iterator_next_1d_block_d
  END INTERFACE

  ! For the sets

  INTERFACE dbcsr_allocate_matrix_set
     MODULE PROCEDURE set_init_1d, set_init_2d
  END INTERFACE

  INTERFACE dbcsr_set_new
     MODULE PROCEDURE set_new_1d, set_new_2d
  END INTERFACE

  INTERFACE dbcsr_set_put
     MODULE PROCEDURE set_put_1d, set_put_2d
  END INTERFACE

  INTERFACE dbcsr_set_get
     MODULE PROCEDURE set_get_1d, set_get_2d
  END INTERFACE

  INTERFACE dbcsr_set_destroy
     MODULE PROCEDURE set_destroy_1d, set_destroy_2d
  END INTERFACE

  INTERFACE dbcsr_set_hold
     MODULE PROCEDURE set_hold_1d, set_hold_2d
  END INTERFACE

  INTERFACE dbcsr_set_release
     MODULE PROCEDURE set_release_1d, set_release_2d
  END INTERFACE

CONTAINS


! *****************************************************************************
!> \brief Registers another reference for a DBCSR matrix
!> \param[in,out] matrix    DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_hold (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_hold', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (matrix%m),&
         cp_failure_level, cp_caller_error, routineN,&
         "Matrix not instantiated")
    matrix%m%refcount = matrix%m%refcount + 1
  END SUBROUTINE dbcsr_hold


! *****************************************************************************
!> \brief Releases a reference for a DBCSR matrix
!>
!> If there are no references left, the matrix is destroyed.
!> \param[in,out] set    DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_release (matrix)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_release', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (matrix%m),&
         cp_warning_level, cp_caller_error, routineN,&
         "Matrix not instantiated")
    IF (ASSOCIATED (matrix%m)) THEN
       !matrix%m%refcount = matrix%m%refcount - 1
       !IF (matrix%m%refcount .EQ. 0) THEN
          CALL dbcsr_destroy (matrix)
       !ENDIF
    ENDIF
  END SUBROUTINE dbcsr_release

! *****************************************************************************
!> \brief Deallocates and destroys a matrix.
!> \param[in,out] matrix      matrix
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_destroy(matrix, error)
    TYPE(dbcsr_obj), INTENT(INOUT)           :: matrix
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_destroy', &
      routineP = moduleN//':'//routineN

    TYPE(cp_error_type)                      :: my_error

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (matrix%m), cp_warning_level, cp_caller_error,&
         routineN, "Can not destroy nonexisting matrix object.")
    IF (ASSOCIATED (matrix%m)) THEN
       matrix%m%refcount = matrix%m%refcount - 1
       IF (matrix%m%refcount .EQ. 0) THEN
          IF (matrix%m%special_memory) THEN
             CALL MP_DEALLOCATE(matrix%m%index)
          ELSE
             DEALLOCATE(matrix%m%index)
          ENDIF
          !IF (matrix%m%special_memory) THEN
          !   CALL MP_DEALLOCATE(matrix%m%data)
          !ELSE
          !   DEALLOCATE(matrix%m%data)
          !ENDIF
          CALL dbcsr_data_release (matrix%m%data_area)
          CALL array_release (matrix%m%row_blk_size)
          CALL array_release (matrix%m%col_blk_size)
          CALL dbcsr_distribution_release(matrix%m%dist)
          matrix%m%valid = .FALSE.
          matrix%m%initialized = .TRUE.
          CALL cp_error_init(my_error)
          IF (ASSOCIATED (matrix%m%dist2d)) &
               CALL distribution_2d_release(matrix%m%dist2d,my_error)
          DEALLOCATE (matrix%m)
       ENDIF
    ENDIF
  END SUBROUTINE dbcsr_destroy

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
  FUNCTION dbcsr_mp_pgrid_equal(mp_env1, mp_env2) RESULT (equal_pgrid)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env1, mp_env2
    LOGICAL                                  :: equal_pgrid

    INTEGER, DIMENSION(:, :), POINTER        :: pgrid1, pgrid2

    IF (dbcsr_mp_nprows (mp_env1) .EQ. dbcsr_mp_nprows (mp_env2)&
         .AND. dbcsr_mp_nprows (mp_env1) .EQ. dbcsr_mp_nprows (mp_env2)) THEN
       pgrid1 => dbcsr_mp_pgrid (mp_env1)
       pgrid2 => dbcsr_mp_pgrid (mp_env2)
       IF (ALL (pgrid1 .EQ. pgrid2)) THEN
          equal_pgrid = .TRUE.
       ELSE
          equal_pgrid = .FALSE.
       ENDIF
    ELSE
       equal_pgrid = .FALSE.
    ENDIF
  END FUNCTION dbcsr_mp_pgrid_equal

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
    TYPE(array_i1d_obj), INTENT(IN), &
      OPTIONAL                               :: local_rows, local_cols

    INTEGER                                  :: i, mypcoor, seq
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
    CALL array_nullify (dist%d%thread_dist)
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
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT)                            :: dist

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

  FUNCTION dbcsr_distribution_thread_dist(dist) RESULT (thread_dist)
    TYPE(dbcsr_distribution_obj), INTENT(IN) :: dist
    TYPE(array_i1d_obj)                      :: thread_dist

!   ---------------------------------------------------------------------------

    thread_dist = dist%d%thread_dist
  END FUNCTION dbcsr_distribution_thread_dist

  SUBROUTINE dbcsr_distribution_make_threads(dist)
    TYPE(dbcsr_distribution_obj), &
      INTENT(INOUT), VOLATILE                :: dist

    INTEGER                                  :: first, last, nrows, nthreads, &
                                                rpt, t
    INTEGER, DIMENSION(:), POINTER           :: td

!   ---------------------------------------------------------------------------

    nthreads = 1
    !$ nthreads = OMP_GET_NUM_THREADS ()
    nrows = dbcsr_distribution_nlocal_rows (dist)
    ALLOCATE (td(nrows))
    rpt = CEILING (REAL(nrows) / REAL(nthreads))
    DO t = 0, nthreads-1
       first = MIN (t*rpt+1, nrows)
       last = MIN ((t+1)*rpt, nrows)
       td(first:last) = t
    END DO
    CALL array_new (dist%d%thread_dist, td, gift=.TRUE.)
  END SUBROUTINE dbcsr_distribution_make_threads
    

! Pertaining to the dbcsr matrix.

  FUNCTION dbcsr_nblkrows_total(matrix) RESULT (nblkrows_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkrows_total

    nblkrows_total = matrix%m%nblkrows_total
  END FUNCTION dbcsr_nblkrows_total

  FUNCTION dbcsr_nblkcols_total(matrix) RESULT (nblkcols_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkcols_total

    nblkcols_total = matrix%m%nblkcols_total
  END FUNCTION dbcsr_nblkcols_total
  FUNCTION dbcsr_nfullrows_total(matrix) RESULT (nfullrows_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullrows_total

    nfullrows_total = matrix%m%nfullrows_total
  END FUNCTION dbcsr_nfullrows_total
  FUNCTION dbcsr_nfullcols_total(matrix) RESULT (nfullcols_total)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullcols_total

    nfullcols_total = matrix%m%nfullcols_total
  END FUNCTION dbcsr_nfullcols_total
  FUNCTION dbcsr_nblkrows_local(matrix) RESULT (nblkrows_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkrows_local

    nblkrows_local = matrix%m%nblkrows_local
  END FUNCTION dbcsr_nblkrows_local
  FUNCTION dbcsr_nblkcols_local(matrix) RESULT (nblkcols_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nblkcols_local

    nblkcols_local = matrix%m%nblkcols_local
  END FUNCTION dbcsr_nblkcols_local
  FUNCTION dbcsr_nfullrows_local(matrix) RESULT (nfullrows_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullrows_local

    nfullrows_local = matrix%m%nfullrows_local
  END FUNCTION dbcsr_nfullrows_local
  FUNCTION dbcsr_nfullcols_local(matrix) RESULT (nfullcols_local)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: nfullcols_local

    nfullcols_local = matrix%m%nfullcols_local
  END FUNCTION dbcsr_nfullcols_local
  FUNCTION dbcsr_max_row_size(matrix) RESULT (max_row_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: max_row_size

    max_row_size = matrix%m%max_rbs
  END FUNCTION dbcsr_max_row_size
  FUNCTION dbcsr_max_col_size(matrix) RESULT (max_col_size)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: max_col_size

    max_col_size = matrix%m%max_cbs
  END FUNCTION dbcsr_max_col_size

  FUNCTION dbcsr_distribution (matrix) RESULT (distribution)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(dbcsr_distribution_obj)             :: distribution

    distribution = matrix%m%dist
  END FUNCTION dbcsr_distribution

  FUNCTION dbcsr_name (matrix) RESULT (name)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    CHARACTER(len=default_string_length)     :: name

    name = matrix%m%name
  END FUNCTION dbcsr_name

  FUNCTION dbcsr_row_block_sizes (matrix) RESULT (row_blk_sizes)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(array_i1d_obj)                      :: row_blk_sizes

    row_blk_sizes = matrix%m%row_blk_size
  END FUNCTION dbcsr_row_block_sizes

  FUNCTION dbcsr_col_block_sizes (matrix) RESULT (col_blk_sizes)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    TYPE(array_i1d_obj)                      :: col_blk_sizes

    col_blk_sizes = matrix%m%col_blk_size
  END FUNCTION dbcsr_col_block_sizes
! *****************************************************************************
!> \brief Gets information about a matrix
!> \param[in] matrix          matrix to query
!> \param[out] nblkrows_total
!> \param[out] nblkcols_total
!> \param[out] nfullrows_total
!> \param[out] nfullcols_total
!> \param[out] nblkrows_local
!> \param[out] nblkrows_local
!> \param[out] nfullrows_local
!> \param[out] nfullrows_local
!> \param[out] my_prow, my_pcol
!> \param[out] local_rows, local_cols
!> \param[out] proc_row_dist, proc_col_dist
!> \param[out] row_blk_size
!> \param[out] col_blk_size
!> \param[out] distribution   the data distribution of the matrix
!> \param[out] name           matrix name
!> \param[out] data_area      data_area
!> \param[out] matrix_type    matrix type (regular, symmetric, see
!>                            dbcsr_types.F for values)
!> \param[out] data_type      data type (single/double precision real/complex)
!> \param[in,out] error       cp2k error
! *****************************************************************************
  SUBROUTINE dbcsr_get_info(matrix, nblkrows_total, nblkcols_total,&
       nfullrows_total, nfullcols_total,&
       nblkrows_local, nblkcols_local,&
       nfullrows_local, nfullcols_local,&
       my_prow, my_pcol,&
       local_rows, local_cols, proc_row_dist, proc_col_dist,&
       row_blk_size, col_blk_size, distribution, name, data_area,&
       matrix_type, data_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER, INTENT(OUT), OPTIONAL :: nblkrows_total, nblkcols_total, &
      nfullrows_total, nfullcols_total, nblkrows_local, nblkcols_local, &
      nfullrows_local, nfullcols_local, my_prow, my_pcol
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: local_rows, local_cols, &
                                                proc_row_dist, proc_col_dist
    TYPE(array_i1d_obj), INTENT(OUT), &
      OPTIONAL                               :: row_blk_size, col_blk_size
    TYPE(dbcsr_distribution_obj), &
      INTENT(OUT), OPTIONAL                  :: distribution
    CHARACTER(len=*), INTENT(OUT), OPTIONAL  :: name
    TYPE(dbcsr_data_obj), INTENT(OUT), &
      OPTIONAL                               :: data_area
    CHARACTER, OPTIONAL                      :: matrix_type, data_type

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_info', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (matrix%m), cp_fatal_level, cp_caller_error,&
         routineP,"Nonexisting matrix")
    CALL cp_assert (matrix%m%valid, cp_warning_level, cp_caller_error,&
         routineP,"Invalid matrix")
    IF (PRESENT (nblkrows_total)) nblkrows_total = matrix%m%nblkrows_total
    IF (PRESENT (nblkcols_total)) nblkcols_total = matrix%m%nblkcols_total
    IF (PRESENT (nfullrows_total)) nfullrows_total = matrix%m%nfullrows_total
    IF (PRESENT (nfullcols_total)) nfullcols_total = matrix%m%nfullcols_total
    IF (PRESENT (nblkrows_local)) nblkrows_local = matrix%m%nblkrows_local
    IF (PRESENT (nblkcols_local)) nblkcols_local = matrix%m%nblkcols_local
    IF (PRESENT (nfullrows_local)) nfullrows_local = matrix%m%nfullrows_local
    IF (PRESENT (nfullcols_local)) nfullcols_local = matrix%m%nfullcols_local
    IF (PRESENT (row_blk_size)) row_blk_size = matrix%m%row_blk_size
    IF (PRESENT (col_blk_size)) col_blk_size = matrix%m%col_blk_size
    IF (PRESENT (distribution)) distribution = matrix%m%dist
    IF (PRESENT (name)) name = matrix%m%name
    IF (PRESENT (data_area)) data_area = matrix%m%data_area
    IF (PRESENT (matrix_type)) matrix_type = matrix%m%matrix_type
    IF (PRESENT (data_type)) data_type = matrix%m%data_type
    IF (PRESENT (local_rows)) &
         local_rows => array_data (dbcsr_distribution_local_rows (matrix%m%dist))
    IF (PRESENT (local_cols)) &
         local_cols => array_data (dbcsr_distribution_local_cols (matrix%m%dist))
    IF (PRESENT (proc_row_dist)) &
         proc_row_dist => array_data (dbcsr_distribution_row_dist (matrix%m%dist))
    IF (PRESENT (proc_col_dist)) &
         proc_col_dist => array_data (dbcsr_distribution_col_dist (matrix%m%dist))
    IF (PRESENT (my_prow)) &
       my_prow = dbcsr_mp_myprow (dbcsr_distribution_mp (matrix%m%dist))
    IF (PRESENT (my_pcol)) &
       my_pcol = dbcsr_mp_mypcol (dbcsr_distribution_mp (matrix%m%dist))
  END SUBROUTINE dbcsr_get_info


! *****************************************************************************
!> \brief Initializes a data area
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_init (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    NULLIFY (area%d)
  END SUBROUTINE dbcsr_data_init
! *****************************************************************************
!> \brief Initializes a data area and all the actual data pointers
!> \param[inout] area         data area
!> \param[in] data_type       select data type to use
!> \param[in] data_size       (optional) allocate this much data
!> \param[in] special_memory  (optional) whether to use specially-allocated
!>                            memory, default is False
! *****************************************************************************
  SUBROUTINE dbcsr_data_new (area, data_type, data_size, special_memory)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    CHARACTER, INTENT(IN)                    :: data_type
    INTEGER, INTENT(IN), OPTIONAL            :: data_size
    LOGICAL, INTENT(IN), OPTIONAL            :: special_memory

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_new', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (.NOT. ASSOCIATED (area%d)) THEN
       ALLOCATE (area%d)
    ENDIF
    NULLIFY (area%d%r_sp)
    NULLIFY (area%d%r_dp)
    NULLIFY (area%d%c_sp)
    NULLIFY (area%d%c_dp)
    area%d%refcount = 1
    IF (PRESENT (special_memory)) THEN
       area%d%special = special_memory
    ELSE
       area%d%special = .FALSE.
    ENDIF
    IF (PRESENT (data_size)) THEN
       IF (area%d%special) THEN
          SELECT CASE (data_type)
          CASE (dbcsr_type_real_4)
             CALL mp_allocate (area%d%r_sp, data_size)
          CASE (dbcsr_type_real_8)
             CALL mp_allocate (area%d%r_dp, data_size)
          CASE (dbcsr_type_complex_4)
             CALL mp_allocate (area%d%c_sp, data_size)
          CASE (dbcsr_type_complex_8)
             CALL mp_allocate (area%d%c_dp, data_size)
          CASE default
             CALL cp_assert (.FALSE., cp_failure_level, cp_caller_error,&
                  routineN, "Invalid matrix type.")
          END SELECT
       ELSE
          SELECT CASE (data_type)
          CASE (dbcsr_type_real_8)
             ALLOCATE (area%d%r_dp(data_size))
          CASE (dbcsr_type_real_4)
             ALLOCATE (area%d%r_sp(data_size))
          CASE (dbcsr_type_complex_8)
             ALLOCATE (area%d%c_dp(data_size))
          CASE (dbcsr_type_complex_4)
             ALLOCATE (area%d%c_sp(data_size))
          CASE default
             CALL cp_assert (.FALSE., cp_failure_level, cp_caller_error,&
                  routineN, "Invalid matrix type.")
          END SELECT
       ENDIF
    ENDIF
    area%d%data_type = data_type
  END SUBROUTINE dbcsr_data_new


! *****************************************************************************
!> \brief Removes a reference and/or clears the data area.
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_release (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_release', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (area%d), cp_caller_error, cp_warning_level,&
         routineN, "Can't release empty data area.")
    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    CALL cp_assert (area%d%refcount .GT. 0, cp_caller_error, cp_warning_level,&
         routineN, "Can't release unreferenced data area.")
    IF (area%d%refcount .GT. 0) THEN
       area%d%refcount = area%d%refcount - 1
    ELSE
       RETURN
    ENDIF
    ! If we're releasing the last reference, then free the memory.
    IF (area%d%refcount .EQ. 0) THEN
       IF (area%d%special) THEN
          SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_4)
             CALL mp_deallocate (area%d%r_sp)
             NULLIFY (area%d%r_sp)
          CASE (dbcsr_type_real_8)
             CALL mp_deallocate (area%d%r_dp)
             NULLIFY (area%d%r_dp)
          CASE (dbcsr_type_complex_4)
             CALL mp_deallocate (area%d%c_sp)
             NULLIFY (area%d%c_sp)
          CASE (dbcsr_type_complex_8)
             CALL mp_deallocate (area%d%c_dp)
             NULLIFY (area%d%c_dp)
          CASE default
             CALL cp_assert (.FALSE., cp_failure_level, cp_caller_error,&
                  routineN, "Invalid data type.")
          END SELECT
       ELSE
          SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_8)
             IF (ASSOCIATED (area%d%r_dp)) &
                  DEALLOCATE (area%d%r_dp)
             NULLIFY (area%d%r_dp)
          CASE (dbcsr_type_real_4)
             IF (ASSOCIATED (area%d%r_sp)) &
                  DEALLOCATE (area%d%r_sp)
             NULLIFY (area%d%r_sp)
          CASE (dbcsr_type_complex_8)
             IF (ASSOCIATED (area%d%c_dp)) &
                  DEALLOCATE (area%d%c_dp)
             NULLIFY (area%d%c_dp)
          CASE (dbcsr_type_complex_4)
             IF (ASSOCIATED (area%d%c_sp)) &
                  DEALLOCATE (area%d%c_sp)
             NULLIFY (area%d%c_sp)
          CASE default
             CALL cp_assert (.FALSE., cp_failure_level, cp_caller_error,&
                  routineN, "Invalid data type.")
          END SELECT
       ENDIF
       DEALLOCATE (area%d)
       NULLIFY (area%d)
    ENDIF
  END SUBROUTINE dbcsr_data_release


! *****************************************************************************
!> \brief Registers another use of the data area
!> \param[inout] area         data area
! *****************************************************************************
  SUBROUTINE dbcsr_data_hold (area)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_hold', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (area%d), cp_caller_error, cp_warning_level,&
         routineN, "Can't hold an empty data area.")
    CALL cp_assert (area%d%refcount .GT. 0, cp_caller_error, cp_warning_level,&
         routineN, "Should not hold an area with zero references.")
    IF (.NOT. ASSOCIATED (area%d)) THEN
       RETURN
    ENDIF
    area%d%refcount = area%d%refcount + 1
  END SUBROUTINE dbcsr_data_hold


! *****************************************************************************
!> \brief Returns the single precision real data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_s (area, DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(KIND=real_4), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%r_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_s
! *****************************************************************************
!> \brief Returns the single precision real data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_m_s (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(KIND=real_4), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (matrix%m%data_area%d)) THEN
       DATA => matrix%m%data_area%d%r_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_m_s

! *****************************************************************************
!> \brief Returns the single precision real data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_s (area) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(KIND=real_4), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%r_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_s

! *****************************************************************************
!> \brief Returns the double precision real data
!> \param[in] area       data area
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_d (area, DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(KIND=real_8), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%r_dp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_d
! *****************************************************************************
!> \brief Returns the double precision real data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_m_d (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    REAL(KIND=real_8), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (matrix%m%data_area%d)) THEN
       DATA => matrix%m%data_area%d%r_dp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_m_d

! *****************************************************************************
!> \brief Returns the double precision real data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_d (area) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(KIND=real_8), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%r_dp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_d

! *****************************************************************************
!> \brief Returns the single precision complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_c (area, DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_c
! *****************************************************************************
!> \brief Returns the single precision complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_m_c (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: DATA

    IF (ASSOCIATED (matrix%m%data_area%d)) THEN
       DATA => matrix%m%data_area%d%c_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_m_c

! *****************************************************************************
!> \brief Returns the single precision complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_c (area) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_c

! *****************************************************************************
!> \brief Returns the double precision complex data
!> \param[in] area       data area
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_z (area, DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_dp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_z
! *****************************************************************************
!> \brief Returns the double precision complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  SUBROUTINE get_data_m_z (matrix, DATA)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: DATA

    IF (ASSOCIATED (matrix%m%data_area%d)) THEN
       DATA => matrix%m%data_area%d%c_dp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END SUBROUTINE get_data_m_z

! *****************************************************************************
!> \brief Returns the double precision complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_p_z (area) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_dp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_p_z

! *****************************************************************************
!> \brief Returns the single precision real data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_s (area, coersion) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(KIND=real_4), INTENT(IN)            :: coersion
    REAL(KIND=real_4), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%r_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_s

! *****************************************************************************
!> \brief Returns the double precision real data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_d (area, coersion) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    REAL(KIND=real_8), INTENT(IN)            :: coersion
    REAL(KIND=real_8), DIMENSION(:), POINTER :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%r_dp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_d

! *****************************************************************************
!> \brief Returns the single precision complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_c (area, coersion) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(KIND=real_4), INTENT(IN)         :: coersion
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_sp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_c

! *****************************************************************************
!> \brief Returns the double precision complex data
!> \param[in] area       data area
!> \param[in] coersion   force datatype
!> \param[out] data      pointer to data
! *****************************************************************************
  FUNCTION dbcsr_get_data_c_z (area, coersion) RESULT (DATA)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    COMPLEX(KIND=real_8), INTENT(IN)         :: coersion
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: DATA

    IF (ASSOCIATED (area%d)) THEN
       DATA => area%d%c_dp
    ELSE
       NULLIFY (DATA)
    ENDIF
  END FUNCTION dbcsr_get_data_c_z


! *****************************************************************************
!> \brief Returns the data area
!> \param matrix     matrix from which to get data
!> \param area       data area
! *****************************************************************************
  FUNCTION dbcsr_data_area (matrix) RESULT (data_area)
    TYPE(dbcsr_type), INTENT(IN)             :: matrix
    TYPE(dbcsr_data_obj)                     :: data_area

    data_area = matrix%data_area
  END FUNCTION dbcsr_data_area

! *****************************************************************************
!> \brief Returns the matrix type
!> \param matrix              query this matrix
!> \param matrix_type         matrix_type (see dbcsr_types.F for possible
!>                            values)
! *****************************************************************************
  PURE FUNCTION dbcsr_get_matrix_type (matrix) RESULT (matrix_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    CHARACTER                                :: matrix_type

    matrix_type = matrix%m%matrix_type
  END FUNCTION dbcsr_get_matrix_type

! *****************************************************************************
!> \brief Returns the data type stored in the matrix
!> \param matrix              query this matrix
!> \param data_type           data_type (see dbcsr_types.F for possible
!>                            values)
! *****************************************************************************
  PURE FUNCTION dbcsr_get_data_type (matrix) RESULT (data_type)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    CHARACTER                                :: data_type

    data_type = matrix%m%data_type
  END FUNCTION dbcsr_get_data_type

! *****************************************************************************
!> \brief Returns the matrix uses specially-allocated memory
!> \param[in] matrix          query this matrix
!> \param[out] uses_special   whether the matrix uses specially allocated
!>                            memory
! *****************************************************************************
  PURE FUNCTION dbcsr_uses_special_memory (matrix) RESULT (uses_special)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    LOGICAL                                  :: uses_special

    uses_special = matrix%m%special_memory
  END FUNCTION dbcsr_uses_special_memory

! *****************************************************************************
!> \brief Returns the allocated data size
!> \param[in] area       data area
!> \param data_size      size of data
! *****************************************************************************
  FUNCTION dbcsr_get_data_size (area) RESULT (data_size)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: area
    INTEGER                                  :: data_size

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_data_size', &
      routineP = moduleN//':'//routineN

    data_size = 0
    IF (ASSOCIATED (area%d)) THEN
       SELECT CASE (area%d%data_type)
          CASE (dbcsr_type_real_8)
             IF (ASSOCIATED (area%d%r_dp))&
                  data_size = SIZE (area%d%r_dp)
          CASE (dbcsr_type_real_4)
             IF (ASSOCIATED (area%d%r_sp))&
                  data_size = SIZE (area%d%r_sp)
          CASE (dbcsr_type_complex_8)
             IF (ASSOCIATED (area%d%c_dp))&
                  data_size = SIZE (area%d%c_dp)
          CASE (dbcsr_type_complex_4)
             IF (ASSOCIATED (area%d%c_sp))&
                  data_size = SIZE (area%d%c_sp)
          CASE default
             CALL cp_assert (.FALSE., cp_failure_level, cp_caller_error,&
                  routineN, "Incorrect data type")
          END SELECT
    ELSE
       CALL cp_assert (.FALSE., cp_warning_level, cp_caller_error, routineN,&
            "Uninitialized data area")
       data_size = 0
    ENDIF
  END FUNCTION dbcsr_get_data_size

! *****************************************************************************
!> \brief Copies a data area, deep by default.
!> \param[inout] target_area     target data area
!> \param[in]    source_area     source data area
!> \param[in]    shallow      shallow copy (default is deep)
! *****************************************************************************
  SUBROUTINE dbcsr_data_copyall (target_area, source_area, shallow)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: target_area
    TYPE(dbcsr_data_obj), INTENT(IN)         :: source_area
    LOGICAL, INTENT(IN), OPTIONAL            :: shallow

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_data_copyall', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: t_data_c
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: t_data_z
    LOGICAL                                  :: shallow_copy
    REAL(KIND=real_4), DIMENSION(:), POINTER :: t_data_r
    REAL(KIND=real_8), DIMENSION(:), POINTER :: t_data_d

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (source_area%d), cp_warning_level,&
         cp_caller_error, routineN, "Attempt to copy unassigned data")
    CALL cp_assert (source_area%d%refcount .GT. 0, cp_warning_level,&
         cp_caller_error, routineN, "Attempt to copy unheld data")
    shallow_copy = .FALSE.
    IF (PRESENT (shallow)) shallow_copy = shallow
    IF (shallow_copy) THEN
       target_area = source_area
       CALL dbcsr_data_hold (target_area)
    ELSE
       IF (.NOT. ASSOCIATED (source_area%d)) THEN
          RETURN
       ENDIF
       SELECT CASE (source_area%d%data_type)
          CASE (dbcsr_type_real_8)
             CALL dbcsr_get_data (target_area, t_data_d)
             CALL bl_copy (dbcsr_get_data_size (source_area),&
                  dbcsr_get_data_p_d (source_area), 1,&
                  t_data_d, 1)
          CASE (dbcsr_type_real_4)
             CALL dbcsr_get_data (target_area, t_data_r)
             CALL bl_copy (dbcsr_get_data_size (source_area),&
                  dbcsr_get_data_p_s (source_area), 1,&
                  t_data_r, 1)
          CASE (dbcsr_type_complex_8)
             CALL dbcsr_get_data (target_area, t_data_z)
             CALL bl_copy (dbcsr_get_data_size (source_area),&
                  dbcsr_get_data_p_z (source_area), 1,&
                  t_data_z, 1)
          CASE (dbcsr_type_complex_4)
             CALL dbcsr_get_data (target_area, t_data_c)
             CALL bl_copy (dbcsr_get_data_size (source_area),&
                  dbcsr_get_data_p_c (source_area), 1,&
                  t_data_c, 1)
          CASE default
             CALL cp_assert (.FALSE., cp_failure_level, cp_caller_error, &
                  routineN, "Invalid matrix type")
          END SELECT
    ENDIF
    ! If we're releasing the last reference, then free the memory.
  END SUBROUTINE dbcsr_data_copyall

! *****************************************************************************
!> \brief Sets a data pointer.
!>
!> Assumes that no memory will be lost when repointing the pointer in the data
!> area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_s (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    REAL(KIND=real_4), DIMENSION(:), POINTER :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_s', &
      routineP = moduleN//':'//routineN

    area%d%r_sp => p
  END SUBROUTINE set_data_p_s

! *****************************************************************************
!> \brief Sets a data pointer.
!>
!> Assumes that no memory will be lost when repointing the pointer in the data
!> area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_d (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    REAL(KIND=real_8), DIMENSION(:), POINTER :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_d', &
      routineP = moduleN//':'//routineN

    area%d%r_dp => p
  END SUBROUTINE set_data_p_d

! *****************************************************************************
!> \brief Sets a data pointer.
!>
!> Assumes that no memory will be lost when repointing the pointer in the data
!> area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_c (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(KIND=real_4), DIMENSION(:), &
      POINTER                                :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_c', &
      routineP = moduleN//':'//routineN

    area%d%c_sp => p
  END SUBROUTINE set_data_p_c

! *****************************************************************************
!> \brief Sets a data pointer.
!>
!> Assumes that no memory will be lost when repointing the pointer in the data
!> area and that the area is initialized.
!> \param[inout] area     target data area
!> \param[in]    p        source data pointer
! *****************************************************************************
  SUBROUTINE set_data_p_z (area, p)
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: area
    COMPLEX(KIND=real_8), DIMENSION(:), &
      POINTER                                :: p

    CHARACTER(len=*), PARAMETER :: routineN = 'set_data_p_z', &
      routineP = moduleN//':'//routineN

    area%d%c_dp => p
  END SUBROUTINE set_data_p_z


! *****************************************************************************
!> \brief Returns the number of blocks in the matrix
!> \param matrix     matrix from which to get data
!> \param area       data area
! *****************************************************************************
  PURE FUNCTION dbcsr_get_num_blocks (matrix) RESULT (num_blocks)
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix
    INTEGER                                  :: num_blocks

    num_blocks = matrix%m%nblks
  END FUNCTION dbcsr_get_num_blocks


! *****************************************************************************
! Iterator functions
! *****************************************************************************


! *****************************************************************************
!> \brief Sets up an iterator
!> \param[out] iterator   the iterator
!> \param[in] matrix      DBCSR matrix
!> \param[in] contiguous_pointers  (optional) Whether returned pointers need
!>                                 to be contiguous; default is FALSE.
!>
!> Contiguous pointers may incur reallocation penalties but enable quick
!> passing of arrays to routines with unspecified interfaces (i.e., direct
!> calls to BLACS or MPI).
! *****************************************************************************
  SUBROUTINE dbcsr_iterator_start (iterator, matrix, contiguous_pointers)
    TYPE(dbcsr_iterator), INTENT(OUT)        :: iterator
    TYPE(dbcsr_obj)                          :: matrix
    LOGICAL, INTENT(IN), OPTIONAL            :: contiguous_pointers

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_iterator_start', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: max_c_size, max_r_size, &
                                                maxsize

!   ---------------------------------------------------------------------------

    IF (PRESENT (contiguous_pointers)) THEN
       iterator%contiguous_pointers = contiguous_pointers
    ELSE
       iterator%contiguous_pointers = .FALSE.
    ENDIF
    iterator%matrix = matrix
    iterator%row = 0
    iterator%pos = 0
    iterator%rbs => array_data (matrix%m%row_blk_size)
    iterator%cbs => array_data (matrix%m%col_blk_size)
    iterator%row_size = 0
    max_r_size = 0
    max_c_size = 0
    IF (matrix%m%nblks .GE. 1) THEN
       iterator%row = 1
       iterator%pos = 1
       CALL find_first_valid_block (iterator%pos, matrix%m%nblks,&
            iterator%matrix%m%blk_p)
       CALL find_proper_row (iterator%pos, iterator%row,&
            iterator%matrix%m%nblkrows_total, iterator%matrix%m%row_p)
       IF (iterator%row .GT. 0) THEN
          iterator%row_size = iterator%rbs(iterator%row)
          max_r_size = matrix%m%max_rbs
          max_c_size = matrix%m%max_cbs
       ENDIF
    ENDIF
    NULLIFY (iterator%buffer_r_sp, iterator%buffer_r_dp,&
         iterator%buffer_c_sp, iterator%buffer_c_dp)
    maxsize = max_r_size * max_c_size
    SELECT CASE (matrix%m%data_area%d%data_type)
    CASE (dbcsr_type_real_4)
       ALLOCATE (iterator%buffer_r_sp (max_r_size, max_c_size))
    CASE (dbcsr_type_real_8)
       ALLOCATE (iterator%buffer_r_dp (max_r_size, max_c_size))
    CASE (dbcsr_type_complex_4)
       ALLOCATE (iterator%buffer_c_sp (max_r_size, max_c_size))
    CASE (dbcsr_type_complex_8)
       ALLOCATE (iterator%buffer_c_dp (max_r_size, max_c_size))
    CASE default
       CALL cp_assert(.FALSE., cp_warning_level, cp_internal_error,&
            routineN, "Invalid data type.")
    END SELECT
  END SUBROUTINE dbcsr_iterator_start

! *****************************************************************************
!> \brief Stops up an iterator
!> \param[out] iterator   the iterator
!> \param[in] matrix      DBCSR matrix
! *****************************************************************************
  SUBROUTINE dbcsr_iterator_stop (iterator)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_iterator_stop', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    iterator%row = 0
    iterator%pos = 0
    SELECT CASE (iterator%matrix%m%data_area%d%data_type)
    CASE (dbcsr_type_real_4)
       CALL cp_assert (ASSOCIATED (iterator%buffer_r_sp), cp_warning_level,&
            cp_internal_error, routineN, "Array should have been allocated")
       IF (ASSOCIATED (iterator%buffer_r_sp)) THEN
          DEALLOCATE (iterator%buffer_r_sp)
       ENDIF
    CASE (dbcsr_type_real_8)
       CALL cp_assert (ASSOCIATED (iterator%buffer_r_dp), cp_warning_level,&
            cp_internal_error, routineN, "Array should have been allocated")
       IF (ASSOCIATED (iterator%buffer_r_dp)) THEN
          DEALLOCATE (iterator%buffer_r_dp)
       ENDIF
    CASE (dbcsr_type_complex_4)
       CALL cp_assert (ASSOCIATED (iterator%buffer_c_sp), cp_warning_level,&
            cp_internal_error, routineN, "Array should have been allocated")
       IF (ASSOCIATED (iterator%buffer_c_sp)) THEN
          DEALLOCATE (iterator%buffer_c_sp)
       ENDIF
    CASE (dbcsr_type_complex_8)
       CALL cp_assert (ASSOCIATED (iterator%buffer_c_dp), cp_warning_level,&
            cp_internal_error, routineN, "Array should have been allocated")
       IF (ASSOCIATED (iterator%buffer_c_dp)) THEN
          DEALLOCATE (iterator%buffer_c_dp)
       ENDIF
    CASE default
       CALL cp_assert(.FALSE., cp_warning_level, cp_internal_error,&
            routineN, "Invalid data type.")
    END SELECT
    !CALL dbcsr_init (iterator%matrix)
    NULLIFY (iterator%buffer_r_sp, iterator%buffer_r_dp,&
         iterator%buffer_c_sp, iterator%buffer_c_dp, iterator%rbs, iterator%cbs)
  END SUBROUTINE dbcsr_iterator_stop


! *****************************************************************************
!> \brief Finds the first valid block, inclusive from the current position.
!>
!> If there is no valid block, pos is set to 0
!> \param[in,out] pos        input: current position; output: next valid
!>                           position or 0
!> \param[in] maxpos         maximal allowed position
!> \param[in] blk_p          block pointers, used to check validity
! *****************************************************************************
  PURE SUBROUTINE find_first_valid_block (pos, maxpos, blk_p)
    INTEGER, INTENT(INOUT)                   :: pos
    INTEGER, INTENT(IN)                      :: maxpos
    INTEGER, DIMENSION(1:maxpos), INTENT(IN) :: blk_p

    CHARACTER(len=*), PARAMETER :: routineN = 'find_first_valid_block', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    DO WHILE (pos .LE. maxpos)
       IF (blk_p(pos) .EQ. 0) THEN
          !WRITE(*,*)routineN//" Skipping deleted block."
          pos = pos+1
       ELSE
          EXIT
       ENDIF
    ENDDO
    IF (pos .GT. maxpos) pos = 0
  END SUBROUTINE find_first_valid_block

! *****************************************************************************
!> \brief Finds the first valid block, inclusive from the current position.
!>
!> If there is no valid block, pos is set to 0
!> \param[in pos             current position
!> \param[in,out] row        input: current row; output: the row corresponding
!>                           to the position
!> \param[in] maxrows        maxmimum row
!> \param[in] row_p          row pointers
! *****************************************************************************
  PURE SUBROUTINE find_proper_row (pos, row, maxrows, row_p)
    INTEGER, INTENT(IN)                      :: pos
    INTEGER, INTENT(INOUT)                   :: row
    INTEGER, INTENT(IN)                      :: maxrows
    INTEGER, DIMENSION(1:maxrows+1), &
      INTENT(IN)                             :: row_p

    CHARACTER(len=*), PARAMETER :: routineN = 'find_proper_row', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (pos .GT. 0) THEN
       DO WHILE (row_p(row+1) .LT. pos)
          row = row+1
          IF (row .GT. maxrows) THEN
             row = 0
             EXIT
          ENDIF
       ENDDO
    ELSE
       row = 0
    ENDIF
  END SUBROUTINE find_proper_row


! *****************************************************************************
!> \brief Gets the next data block, double precision real
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] block         pointer to the data block
!> \param[out] block_number  (optional) block number
! *****************************************************************************
  SUBROUTINE iterator_next_block_d (iterator, row, column, block, block_number)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(KIND=real_8), DIMENSION(:, :), &
      POINTER                                :: block
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_block_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, csize, nze, rsize
    REAL(KIND=real_8), DIMENSION(:), POINTER :: lin_blk_p, rbp

!   ---------------------------------------------------------------------------
! If we're pointing to a valid block, return that block.

    IF (iterator%pos .LE. iterator%matrix%m%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%matrix%m%col_i(iterator%pos)
       !IF (PRESENT (block)) THEN
          blk_p = iterator%matrix%m%blk_p(iterator%pos)
          rsize = iterator%row_size
          csize = iterator%cbs(column)
          nze = rsize * csize
          rbp => dbcsr_get_data_p (iterator%matrix%m%data_area, 0.0_real_8)
#if defined(__PTR_RANK_REMAP)
          block(1:rsize, 1:csize) => rbp(blk_p:blk_p+nze-1)
#elif !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
          lin_blk_p => rbp(blk_p:blk_p+nze-1)
          CALL pointer_d_rank_remap2 (block, rsize, csize, lin_blk_p)
#else
          IF (iterator%contiguous_pointers) THEN
             IF (SIZE (iterator%buffer_r_dp, 1) .NE. rsize) THEN
                DEALLOCATE (iterator%buffer_r_dp)
                ALLOCATE (iterator%buffer_r_dp(&
                     rsize, iterator%matrix%m%max_cbs))
             ENDIF
          ENDIF
          iterator%buffer_r_dp =&
               RESHAPE (rbp(blk_p:blk_p+nze-1), (/rsize, csize/))
          block => iterator%buffer_r_dp(1:rsize, 1:csize)
#endif /* __PTR_RANK_REMAP || !__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION */
       !ENDIF
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       iterator%pos = iterator%pos+1
       CALL find_first_valid_block (iterator%pos, iterator%matrix%m%nblks,&
            iterator%matrix%m%blk_p)
       CALL find_proper_row (iterator%pos, iterator%row,&
            iterator%matrix%m%nblkrows_total, iterator%matrix%m%row_p)
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
    ELSE
       row = 0
       column = 0
       !IF (PRESENT (block)) NULLIFY (block)
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_block_d

! *****************************************************************************
!> \brief Gets the next data block, double precision real
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] block         pointer to the data block
!> \param[out] block_number  (optional) block number
! *****************************************************************************
  SUBROUTINE iterator_next_1d_block_d (iterator, row, column, block, block_number)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column
    REAL(KIND=real_8), DIMENSION(:), POINTER :: block
    INTEGER, INTENT(OUT), OPTIONAL           :: block_number

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_1d_block_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk_p, csize, nze, rsize
    REAL(KIND=real_8), DIMENSION(:), POINTER :: rbp

!   ---------------------------------------------------------------------------
! If we're pointing to a valid block, return that block.

    IF (iterator%pos .LE. iterator%matrix%m%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%matrix%m%col_i(iterator%pos)
       !IF (PRESENT (block)) THEN
          blk_p = iterator%matrix%m%blk_p(iterator%pos)
          rsize = iterator%row_size
          csize = iterator%cbs(column)
          nze = rsize * csize
          rbp => dbcsr_get_data_p (iterator%matrix%m%data_area, 0.0_real_8)
          block => rbp(blk_p:blk_p+nze-1)
       !ENDIF
       IF (PRESENT (block_number)) block_number = iterator%pos
       ! Move to the next non-deleted position.
       iterator%pos = iterator%pos+1
       CALL find_first_valid_block (iterator%pos, iterator%matrix%m%nblks,&
            iterator%matrix%m%blk_p)
       CALL find_proper_row (iterator%pos, iterator%row,&
            iterator%matrix%m%nblkrows_total, iterator%matrix%m%row_p)
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
    ELSE
       row = 0
       column = 0
       !IF (PRESENT (block)) NULLIFY (block)
       IF (PRESENT (block_number)) block_number = 0
    ENDIF
  END SUBROUTINE iterator_next_1d_block_d


! *****************************************************************************
!> \brief Gets the index information of the next block, no data.
!> \param[in,out] iterator   the iterator
!> \param[out] row           row of the data block
!> \param[out] column        column of the data block
!> \param[out] blk           block number
!> \param[out] blk_p         (optional) index into block data array
! *****************************************************************************
  SUBROUTINE iterator_next_block_index (iterator, row, column, blk, blk_p)
    TYPE(dbcsr_iterator), INTENT(INOUT)      :: iterator
    INTEGER, INTENT(OUT)                     :: row, column, blk
    INTEGER, INTENT(OUT), OPTIONAL           :: blk_p

    CHARACTER(len=*), PARAMETER :: routineN = 'iterator_next_block_index', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (iterator%pos .LE. iterator%matrix%m%nblks&
         .AND. iterator%pos .NE. 0) THEN
       row = iterator%row
       column = iterator%matrix%m%col_i(iterator%pos)
       blk = iterator%pos
       IF (PRESENT (blk_p)) blk_p = iterator%matrix%m%blk_p(iterator%pos)
       iterator%pos = iterator%pos+1
       CALL find_first_valid_block (iterator%pos, iterator%matrix%m%nblks,&
            iterator%matrix%m%blk_p)
       CALL find_proper_row (iterator%pos, iterator%row,&
            iterator%matrix%m%nblkrows_total, iterator%matrix%m%row_p)
       IF (iterator%row .GT. 0) iterator%row_size = iterator%rbs(iterator%row)
    ELSE
       row = 0
       column = 0
    ENDIF
    
  END SUBROUTINE iterator_next_block_index


! *****************************************************************************
!> \brief Returns whether there any blocks left in the iterator.
!> \param[in] iterator       the iterator
!> \param[in]
! *****************************************************************************
  FUNCTION dbcsr_iterator_blocks_left (iterator) RESULT (blocks_left)
    TYPE(dbcsr_iterator), INTENT(IN)         :: iterator
    LOGICAL                                  :: blocks_left

    blocks_left = iterator%pos .NE. 0
  END FUNCTION dbcsr_iterator_blocks_left


! *****************************************************************************
! Sets
! *****************************************************************************


! *****************************************************************************
!> \brief Initializes a 1-d set of DBCSR matrices
!> \param[out] set    the set of matrices
!> \param[in] n       (optional) number of matrices to be included
! *****************************************************************************
  SUBROUTINE set_init_1d (set, n)
    TYPE(dbcsr_1d_set_obj), INTENT(OUT)      :: set
    INTEGER, INTENT(IN), OPTIONAL            :: n

    NULLIFY (set%s, set%refcount)
    IF (PRESENT (n)) CALL dbcsr_set_new (set, n)
  END SUBROUTINE set_init_1d

! *****************************************************************************
!> \brief Initializes a 2-d set of DBCSR matrices
!> \param[out] set    the set of matrices
!> \param[in] n       (optional) number of matrices to be included
! *****************************************************************************
  SUBROUTINE set_init_2d (set, n)
    TYPE(dbcsr_2d_set_obj), INTENT(OUT)      :: set
    INTEGER, DIMENSION(2), INTENT(IN), &
      OPTIONAL                               :: n

    NULLIFY (set%s, set%refcount)
    IF (PRESENT (n)) CALL dbcsr_set_new (set, n)
  END SUBROUTINE set_init_2d

! *****************************************************************************
!> \brief Allocates a 1-d set of DBCSR matrices
!> \param[in,out] set  the set of matrices
!> \param[in] n        number of matrices to be included
! *****************************************************************************
  SUBROUTINE set_new_1d (set, n)
    TYPE(dbcsr_1d_set_obj), INTENT(INOUT)    :: set
    INTEGER, INTENT(IN)                      :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'set_new_1d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (.NOT. ASSOCIATED (set%s), cp_warning_level,&
         cp_caller_error, routineN,&
         "Set was already allocated, memory leakage will occur")
    ALLOCATE (set%s(n))
    ALLOCATE (set%refcount)
    set%refcount = 1
  END SUBROUTINE set_new_1d

! *****************************************************************************
!> \brief Allocates a 2-d set of DBCSR matrices
!> \param[in,out] set  the set of matrices
!> \param[in] n        number of matrices to be included
! *****************************************************************************
  SUBROUTINE set_new_2d (set, n)
    TYPE(dbcsr_2d_set_obj), INTENT(INOUT)    :: set
    INTEGER, DIMENSION(2), INTENT(IN)        :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'set_new_2d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (.NOT. ASSOCIATED (set%s), cp_warning_level,&
         cp_caller_error, routineN,&
         "Set was already allocated, memory leakage will occur")
    ALLOCATE (set%s(n(1), n(2)))
    ALLOCATE (set%refcount)
    set%refcount = 1
  END SUBROUTINE set_new_2d

! *****************************************************************************
!> \brief Inserts a DBCSR matrix into the set
!>
!> The matrix is inserted into the set and its reference is registered.
!> \param[in,out] set    the set of matrices
!> \param[in] position   array position into which matrix is inserted
!> \param[in] matrix     matrices to be included
! *****************************************************************************
  SUBROUTINE set_put_1d (set, position, matrix)
    TYPE(dbcsr_1d_set_obj), INTENT(INOUT)    :: set
    INTEGER, INTENT(IN)                      :: position
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'set_put_1d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%s), cp_fatal_level, cp_caller_error,&
         routineN, "Can not insert into unallocated set")
    CALL cp_assert(position .GE. LBOUND (set%s, 1)&
         .AND. position .LE. UBOUND (set%s, 1),&
         cp_failure_level, cp_caller_error, routineN,&
         "Insertion position out of bounds.")
    set%s(position) = matrix
    CALL dbcsr_hold (set%s(position))
  END SUBROUTINE set_put_1d

! *****************************************************************************
!> \brief Inserts a DBCSR matrix into the set
!>
!> The matrix is inserted into the set and its reference is registered.
!> \param[in,out] set the set of matrices
!> \param[in] position   array position into which matrix is inserted
!> \param[in] matrix  matrices to be included
! *****************************************************************************
  SUBROUTINE set_put_2d (set, position, matrix)
    TYPE(dbcsr_2d_set_obj), INTENT(INOUT)    :: set
    INTEGER, DIMENSION(2), INTENT(IN)        :: position
    TYPE(dbcsr_obj), INTENT(IN)              :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'set_put_2d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%s), cp_fatal_level, cp_caller_error,&
         routineN, "Can not insert into unallocated set")
    CALL cp_assert(position(1) .GE. LBOUND (set%s, 1)&
         .AND. position(1) .LE. UBOUND (set%s, 1)&
         .AND. position(2) .GE. LBOUND (set%s, 2)&
         .AND. position(2) .LE. UBOUND (set%s, 2),&
         cp_failure_level, cp_caller_error, routineN,&
         "Insertion position out of bounds.")
    set%s(position(1), position(2)) = matrix
    CALL dbcsr_hold (set%s(position(1), position(2)))
  END SUBROUTINE set_put_2d

! *****************************************************************************
!> \brief Gets a DBCSR matrix from the set
!> \param[in] set     the set of matrices
!> \param[in] position   array position into which matrix is inserted
!> \param matrix      matrix (no new reference is registered)
! *****************************************************************************
  FUNCTION set_get_1d (set, position) RESULT (matrix)
    TYPE(dbcsr_1d_set_obj), INTENT(IN)       :: set
    INTEGER, INTENT(IN)                      :: position
    TYPE(dbcsr_obj)                          :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'set_get_1d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%s), cp_fatal_level, cp_caller_error,&
         routineN, "Can not read from an unallocated set")
    CALL cp_assert (position .GE. LBOUND (set%s,1)&
         .AND. position .LE. UBOUND (set%s,1),&
         cp_failure_level, cp_caller_error, routineN,&
         "Read position out of bounds.")
    matrix = set%s(position)
  END FUNCTION set_get_1d

! *****************************************************************************
!> \brief Gets a DBCSR matrix from the set
!> \param[in] set     the set of matrices
!> \param[in] position   array position into which matrix is inserted
!> \param matrix      matrix (no new reference is registered)
! *****************************************************************************
  FUNCTION set_get_2d (set, position) RESULT (matrix)
    TYPE(dbcsr_2d_set_obj), INTENT(IN)       :: set
    INTEGER, DIMENSION(2), INTENT(IN)        :: position
    TYPE(dbcsr_obj)                          :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'set_get_2d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%s), cp_fatal_level, cp_caller_error,&
         routineN, "Can not read from an unallocated set")
    CALL cp_assert(position(1) .GE. LBOUND (set%s, 1)&
         .AND. position(1) .LE. UBOUND (set%s, 1)&
         .AND. position(2) .GE. LBOUND (set%s, 2)&
         .AND. position(2) .LE. UBOUND (set%s, 2),&
         cp_failure_level, cp_caller_error, routineN,&
         "Read position out of bounds.")
    matrix = set%s(position(1), position(2))
  END FUNCTION set_get_2d

! *****************************************************************************
!> \brief Destroys a 1-d set and releases all held resources
!> \param[in,out] set    the set of matrices
! *****************************************************************************
  SUBROUTINE set_destroy_1d (set)
    TYPE(dbcsr_1d_set_obj), INTENT(INOUT)    :: set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_destroy_1d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%s), cp_fatal_level, cp_caller_error,&
         routineN, "Can not destroy an unallocated set")
    CALL cp_assert (ASSOCIATED (set%refcount),&
         cp_failure_level, cp_caller_error, routineN,&
         "Reference count does not exist.")
    CALL cp_assert (set%refcount .EQ. 0, cp_warning_level, cp_caller_error,&
         routineN, "References are still held")
    DO i = LBOUND (set%s, 1), UBOUND (set%s, 1)
       CALL dbcsr_release (set%s(i))
    ENDDO
    DEALLOCATE (set%s)
    DEALLOCATE (set%refcount)
    CALL dbcsr_allocate_matrix_set (set)
  END SUBROUTINE set_destroy_1d

! *****************************************************************************
!> \brief Destroys a 2-d set and releases all held resources
!> \param[in,out] set    the set of matrices
! *****************************************************************************
  SUBROUTINE set_destroy_2d (set)
    TYPE(dbcsr_2d_set_obj), INTENT(INOUT)    :: set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_destroy_2d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%s), cp_fatal_level, cp_caller_error,&
         routineN, "Can not destroy an unallocated set")
    CALL cp_assert (ASSOCIATED (set%refcount),&
         cp_failure_level, cp_caller_error, routineN,&
         "Reference count does not exist.")
    CALL cp_assert (set%refcount .EQ. 0, cp_warning_level, cp_caller_error,&
         routineN, "References are still held")
    DO i = LBOUND (set%s, 1), UBOUND (set%s, 1)
       DO j = LBOUND (set%s, 2), UBOUND (set%s, 2)
          CALL dbcsr_release (set%s(i, j))
       ENDDO
    ENDDO
    DEALLOCATE (set%s)
    DEALLOCATE (set%refcount)
    CALL dbcsr_allocate_matrix_set (set)
  END SUBROUTINE set_destroy_2d


! *****************************************************************************
!> \brief Registers another reference for a 1-d set
!> \param[in,out] set    the set of matrices
! *****************************************************************************
  SUBROUTINE set_hold_1d (set)
    TYPE(dbcsr_1d_set_obj), INTENT(INOUT)    :: set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_hold_1d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%refcount),&
         cp_failure_level, cp_caller_error, routineN,&
         "Reference count does not exist.")
    set%refcount = set%refcount + 1
  END SUBROUTINE set_hold_1d

! *****************************************************************************
!> \brief Registers another reference for a 2-d set
!> \param[in,out] set    the set of matrices
! *****************************************************************************
  SUBROUTINE set_hold_2d (set)
    TYPE(dbcsr_2d_set_obj), INTENT(INOUT)    :: set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_hold_2d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%refcount),&
         cp_failure_level, cp_caller_error, routineN,&
         "Reference count does not exist.")
    set%refcount = set%refcount + 1
  END SUBROUTINE set_hold_2d


! *****************************************************************************
!> \brief Releases a reference for a 1-d set
!>
!> If there are no references left, the set is destroyed.
!> \param[in,out] set    the set of matrices
! *****************************************************************************
  SUBROUTINE set_release_1d (set)
    TYPE(dbcsr_1d_set_obj), INTENT(INOUT)    :: set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_release_1d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%refcount),&
         cp_failure_level, cp_caller_error, routineN,&
         "Reference count does not exist.")
    set%refcount = set%refcount - 1
    IF (set%refcount .EQ. 0) THEN
       CALL dbcsr_set_destroy (set)
    ENDIF
  END SUBROUTINE set_release_1d

! *****************************************************************************
!> \brief Releases a reference for a 2-d set
!>
!> If there are no references left, the set is destroyed.
!> \param[in,out] set    the set of matrices
! *****************************************************************************
  SUBROUTINE set_release_2d (set)
    TYPE(dbcsr_2d_set_obj), INTENT(INOUT)    :: set

    CHARACTER(len=*), PARAMETER :: routineN = 'set_release_2d', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (set%refcount),&
         cp_failure_level, cp_caller_error, routineN,&
         "Reference count does not exist.")
    set%refcount = set%refcount - 1
    IF (set%refcount .EQ. 0) THEN
       CALL dbcsr_set_destroy (set)
    ENDIF
  END SUBROUTINE set_release_2d

#if !defined(__NO_ASSUMED_SIZE_NOCOPY_ASSUMPTION)
! *****************************************************************************
!> \brief Sets a rank-2 pointer to rank-1 data using ugly hacks.
!>
!> Merge this with the same routine in dbcsr_util* upon refactoring.
! *****************************************************************************
  SUBROUTINE pointer_d_rank_remap2 (r2p, d1, d2, r1p)
    REAL(KIND=real_8), DIMENSION(:, :), &
      POINTER                                :: r2p
    INTEGER, INTENT(IN)                      :: d1, d2
    REAL(KIND=real_8), DIMENSION(d1, *), &
      TARGET                                 :: r1p

    r2p => r1p(1:d1, 1:d2)
  END SUBROUTINE pointer_d_rank_remap2
#endif


! *****************************************************************************
! Mutable data
! *****************************************************************************

! *****************************************************************************
!> \brief Initializes a mutable data type
!> \param[out] mutable        mutable data
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_init (mutable)
    TYPE(dbcsr_mutable_obj), INTENT(OUT)     :: mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_init', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    NULLIFY (mutable%m)
  END SUBROUTINE dbcsr_mutable_init

! *****************************************************************************
!> \brief Destroys a mutable data type
!> \param[in,out] mutable     mutable data
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_destroy (mutable)
    TYPE(dbcsr_mutable_obj), INTENT(INOUT)   :: mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_destroy', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED (mutable%m)) THEN
       !CALL cp_assert (mutable%m%refcount .EQ. 0, cp_warning_level,&
       !     cp_caller_error, routineN, "Destroying with non-0 reference count")
       CALL btree_destroy_s (mutable%m%btree_s)
       CALL btree_destroy_d (mutable%m%btree_d)
       CALL btree_destroy_c (mutable%m%btree_c)
       CALL btree_destroy_z (mutable%m%btree_z)
       DEALLOCATE (mutable%m)
    ENDIF
    NULLIFY (mutable%m)
  END SUBROUTINE dbcsr_mutable_destroy


! *****************************************************************************
!> \brief Registers another reference to the mutable data type
!> \param[in,out] mutable     mutable data
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_hold (mutable)
    TYPE(dbcsr_mutable_obj), INTENT(INOUT)   :: mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_hold', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (mutable%m), cp_fatal_level,&
         cp_caller_error, routineN, "Mutable data area not instantiated")
    mutable%m%refcount = mutable%m%refcount + 1
  END SUBROUTINE dbcsr_mutable_hold

! *****************************************************************************
!> \brief Deregisters a reference to the mutable data type
!>
!> The object is destroy when there is no reference to it left.
!> \param[in,out] mutable     mutable data
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_release (mutable)
    TYPE(dbcsr_mutable_obj), INTENT(INOUT)   :: mutable

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_release', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (ASSOCIATED (mutable%m), cp_fatal_level,&
         cp_caller_error, routineN, "Mutable data area not instantiated")
    mutable%m%refcount = mutable%m%refcount - 1
    IF (mutable%m%refcount .EQ. 0) THEN
       CALL dbcsr_mutable_destroy (mutable)
    ENDIF
  END SUBROUTINE dbcsr_mutable_release

! *****************************************************************************
!> \brief Creates a new mutable instance.
!>
!> \param[in,out] mutable     mutable data
!> \param[in] data_type       data type to be stored here (see dbcsr_types for
!>                            possibilities)
! *****************************************************************************
  SUBROUTINE dbcsr_mutable_new (mutable, data_type)
    TYPE(dbcsr_mutable_obj), INTENT(INOUT)   :: mutable
    CHARACTER, INTENT(IN)                    :: data_type

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mutable_new', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL cp_assert (.NOT.ASSOCIATED (mutable%m), cp_fatal_level,&
         cp_caller_error, routineN, "Mutable data area already instantiated")
    CALL cp_assert (data_type.EQ.dbcsr_type_real_4&
         .OR. data_type.EQ.dbcsr_type_real_8&
         .OR. data_type.EQ.dbcsr_type_complex_4&
         .OR. data_type.EQ.dbcsr_type_complex_8, cp_fatal_level,&
         cp_wrong_args_error, routineN, "Invalid data type")
    ALLOCATE (mutable%m)
    mutable%m%data_type = data_type
    CALL btree_new_s (mutable%m%btree_s)
    CALL btree_new_d (mutable%m%btree_d)
    CALL btree_new_c (mutable%m%btree_c)
    CALL btree_new_z (mutable%m%btree_z)
  END SUBROUTINE dbcsr_mutable_new

! *****************************************************************************
!> \brief Deregisters a reference to the mutable data type
!>
!> The object is destroy when there is no reference to it left.
!> \param[in] mutable         mutable data
!> \param[out] instantiated   whether the object is instantiated
! *****************************************************************************
  PURE FUNCTION dbcsr_mutable_instantiated (mutable) RESULT (instantiated)
    TYPE(dbcsr_mutable_obj), INTENT(IN)      :: mutable
    LOGICAL                                  :: instantiated

!   ---------------------------------------------------------------------------

    instantiated = ASSOCIATED (mutable%m)
  END FUNCTION dbcsr_mutable_instantiated

END MODULE dbcsr_methods
