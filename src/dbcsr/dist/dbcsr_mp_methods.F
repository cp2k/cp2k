!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Routines related to message passing
!> \author  Ole Schuett
!> \date    10-2014
! *****************************************************************************
MODULE dbcsr_mp_methods
  USE dbcsr_config,                    ONLY: use_subcommunicators
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error
  USE dbcsr_methods,                   ONLY: dbcsr_mp_grid_remove,&
                                             dbcsr_mp_release
  USE dbcsr_types,                     ONLY: dbcsr_mp_obj
  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_cart_sub,&
                                             mp_comm_free

  !$ USE OMP_LIB
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mp_methods'

  INTEGER, PUBLIC, SAVE :: dbcsr_matrix_counter = 111111

  PUBLIC :: dbcsr_mp_new, dbcsr_mp_hold, dbcsr_mp_release,&
            dbcsr_mp_pgrid, dbcsr_mp_numnodes, dbcsr_mp_mynode, dbcsr_mp_group,&
            dbcsr_mp_new_transposed, dbcsr_mp_nprows, dbcsr_mp_npcols,&
            dbcsr_mp_myprow, dbcsr_mp_mypcol,&
            dbcsr_mp_my_row_group, dbcsr_mp_my_col_group,&
            dbcsr_mp_has_subgroups, dbcsr_mp_get_process,&
            dbcsr_mp_grid_setup, dbcsr_mp_grid_remove,&
            dbcsr_mp_init, dbcsr_mp_active


CONTAINS

! *****************************************************************************
!> \brief Initializes a new process grid
!> \param mp_env ...
! *****************************************************************************
  SUBROUTINE dbcsr_mp_init (mp_env)
    TYPE(dbcsr_mp_obj), INTENT(OUT)          :: mp_env

    NULLIFY (mp_env%mp)
  END SUBROUTINE dbcsr_mp_init

! *****************************************************************************
!> \brief Checks whether this process is part of the message passing environment
!> \param mp_env ...
!> \retval active ...
! *****************************************************************************
  FUNCTION dbcsr_mp_active (mp_env) RESULT (active)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL                                  :: active

    active = ASSOCIATED (mp_env%mp)
  END FUNCTION dbcsr_mp_active


! *****************************************************************************
!> \brief Creates new process grid
!> \param[out] mp_env         multiprocessor environment
!> \param[in] pgrid           process grid
!> \param mp_group ...
!> \param[in] mynode          my processor number
!> \param[in] numnodes        total number of processors (processes)
!> \param myprow ...
!> \param mypcol ...
!> \par History
!>      UB (2010-02-04)       Duplicates own communicator and sets up
!>                            cartesian grid
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
     mp_env%mp%subgroups_defined = .FALSE.
!!!! KG workaround in place.
!!!! This will be the replacement:
!    ALLOCATE(mp_env%mp)
!    mp_env%mp%refcount = 1
!    ndims = 2
!    dims(1:2) = (/ SIZE (pgrid, 1), SIZE (pgrid, 2) /)
!    CALL mp_cart_create (mp_group, ndims,&
!         dims, my_pos,&
!         mp_env%mp%mp_group)
!    CALL mp_environ (mp_env%mp%numnodes, mp_env%mp%mynode, mp_env%mp%mp_group)
!    mp_env%mp%myprow = my_pos(1)
!    mp_env%mp%mypcol = my_pos(2)
!    ALLOCATE (mp_env%mp%pgrid (0:SIZE(pgrid, 1)-1, 0:SIZE(pgrid, 2)-1 ))
!    column_loop: DO pcol = 0, SIZE (mp_env%mp%pgrid, 2)-1
!       row_loop: DO prow = 0, SIZE (mp_env%mp%pgrid, 1)-1
!          CALL mp_cart_rank (mp_env%mp%mp_group, (/ prow, pcol /),&
!               mp_env%mp%pgrid (prow, pcol))
!       ENDDO row_loop
!    ENDDO column_loop
!    mp_env%mp%subgroups_defined = .FALSE.
  END SUBROUTINE dbcsr_mp_new

! *****************************************************************************
!> \brief Sets up MPI cartesian process grid
!> \param[in,out] mp_env      multiprocessor environment
!> \param[in] force           (optional) force creation of subcommunicators
! *****************************************************************************
  SUBROUTINE dbcsr_mp_grid_setup(mp_env, force)
    TYPE(dbcsr_mp_obj), INTENT(INOUT)        :: mp_env
    LOGICAL, INTENT(IN), OPTIONAL            :: force

    INTEGER                                  :: ndims, tmp_group
    INTEGER, DIMENSION(2)                    :: dims, my_pos
    LOGICAL                                  :: should_make
    LOGICAL, DIMENSION(2)                    :: remain
    TYPE(dbcsr_error_type)                   :: error

!   ---------------------------------------------------------------------------
! Intel MPI Workaround
!CALL dbcsr_assert (..NOT. mp_env%mp%subgroups_defined,&
!     dbcsr_fatal_level, dbcsr_internal_error, "mp_grid_setup",&
!     "Subcommunicators should be turned off")

    IF (PRESENT (force)) THEN
       should_make = use_subcommunicators .OR. force
    ELSE
       should_make = use_subcommunicators
    ENDIF

    IF (.NOT. mp_env%mp%subgroups_defined .AND. should_make) THEN
       ! KG workaround.
       ! This will be deleted (replaced by code in mp_new).
       ndims = 2
       dims(1:2) = (/ SIZE (mp_env%mp%pgrid, 1), SIZE (mp_env%mp%pgrid, 2) /)
       CALL mp_cart_create (mp_env%mp%mp_group, ndims,&
            dims, my_pos,&
            tmp_group)
       CALL dbcsr_assert (my_pos(1) .EQ. mp_env%mp%myprow,&
            dbcsr_fatal_level, dbcsr_internal_error, "mp_grid_setup",&
            "Got different MPI process grid",__LINE__,error)
       CALL dbcsr_assert (my_pos(2) .EQ. mp_env%mp%mypcol,&
            dbcsr_fatal_level, dbcsr_internal_error, "mp_grid_setup",&
            "Got different MPI process grid",__LINE__,error)
       !
       remain = (/.FALSE., .TRUE./)
       CALL mp_cart_sub (tmp_group, remain, mp_env%mp%prow_group)
       remain = (/.TRUE., .FALSE./)
       CALL mp_cart_sub (tmp_group, remain, mp_env%mp%pcol_group)
       CALL mp_comm_free (tmp_group)
       mp_env%mp%subgroups_defined = .TRUE.
    ENDIF
  END SUBROUTINE dbcsr_mp_grid_setup

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
!> \brief ...
!> \param mp_env ...
!> \param prow ...
!> \param pcol ...
!> \retval process ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_get_process(mp_env, prow, pcol) RESULT (process)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, INTENT(IN)                      :: prow, pcol
    INTEGER                                  :: process

    process = mp_env%mp%pgrid(prow, pcol)
  END FUNCTION dbcsr_mp_get_process

! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \param node ...
!> \param prow ...
!> \param pcol ...
! *****************************************************************************
  PURE SUBROUTINE dbcsr_mp_get_coordinates(mp_env, node, prow, pcol)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, INTENT(IN)                      :: node
    INTEGER, INTENT(OUT)                     :: prow, pcol

    INTEGER                                  :: tmp

    pcol = MOD (node, SIZE(mp_env%mp%pgrid,2))
    prow = node / SIZE(mp_env%mp%pgrid,2)
    IF (mp_env%mp%pgrid(prow, pcol) .NE. node) THEN
       tmp = pcol ; pcol = prow ; prow = tmp
    ENDIF
  END SUBROUTINE dbcsr_mp_get_coordinates

! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval pgrid ...
! *****************************************************************************
  FUNCTION dbcsr_mp_pgrid(mp_env) RESULT (pgrid)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER, DIMENSION(:, :), POINTER        :: pgrid

    pgrid => mp_env%mp%pgrid
  END FUNCTION dbcsr_mp_pgrid
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval numnodes ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_numnodes(mp_env) RESULT (numnodes)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: numnodes

    numnodes = mp_env%mp%numnodes
  END FUNCTION dbcsr_mp_numnodes
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval mynode ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_mynode(mp_env) RESULT (mynode)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mynode

    mynode = mp_env%mp%mynode
  END FUNCTION dbcsr_mp_mynode
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval mp_group ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_group(mp_env) RESULT (mp_group)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mp_group

    mp_group = mp_env%mp%mp_group
  END FUNCTION dbcsr_mp_group
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval nprows ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_nprows(mp_env) RESULT (nprows)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: nprows

    nprows = SIZE (mp_env%mp%pgrid, 1)
  END FUNCTION dbcsr_mp_nprows
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval npcols ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_npcols(mp_env) RESULT (npcols)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: npcols

    npcols = SIZE (mp_env%mp%pgrid, 2)
  END FUNCTION dbcsr_mp_npcols
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval myprow ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_myprow(mp_env) RESULT (myprow)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: myprow

    myprow = mp_env%mp%myprow
  END FUNCTION dbcsr_mp_myprow
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval mypcol ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_mypcol(mp_env) RESULT (mypcol)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: mypcol

    mypcol = mp_env%mp%mypcol
  END FUNCTION dbcsr_mp_mypcol
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval has_subgroups ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_has_subgroups(mp_env) RESULT (has_subgroups)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    LOGICAL                                  :: has_subgroups

    has_subgroups = mp_env%mp%subgroups_defined
  END FUNCTION dbcsr_mp_has_subgroups
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval row_group ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_my_row_group(mp_env) RESULT (row_group)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: row_group

    row_group = mp_env%mp%prow_group
  END FUNCTION dbcsr_mp_my_row_group
! *****************************************************************************
!> \brief ...
!> \param mp_env ...
!> \retval col_group ...
! *****************************************************************************
  PURE FUNCTION dbcsr_mp_my_col_group(mp_env) RESULT (col_group)
    TYPE(dbcsr_mp_obj), INTENT(IN)           :: mp_env
    INTEGER                                  :: col_group

    col_group = mp_env%mp%pcol_group
  END FUNCTION dbcsr_mp_my_col_group

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

END MODULE dbcsr_mp_methods
