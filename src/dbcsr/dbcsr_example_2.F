!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR example 2
!>          This example shows how to set a dbcsr matrix
!> \author  V. Weber
!> \date    2010
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
PROGRAM dbcsr_example_2
  USE dbcsr_api,                       ONLY: &
       array_data, array_i1d_obj, array_new, array_nullify, array_release, &
       dbcsr_col_block_sizes, dbcsr_create, dbcsr_distribution, &
       dbcsr_distribution_mp, dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_release, dbcsr_error_type, dbcsr_finalize, &
       dbcsr_finalize_lib, dbcsr_get_stored_coordinates, dbcsr_init, &
       dbcsr_init_lib, dbcsr_max_col_size, dbcsr_max_row_size, &
       dbcsr_mp_mynode, dbcsr_mp_new, dbcsr_mp_obj, dbcsr_mp_release, &
       dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_obj, dbcsr_print, &
       dbcsr_put_block, dbcsr_release, dbcsr_row_block_sizes, &
       dbcsr_set_default_config, dbcsr_type_no_symmetry, dbcsr_type_real_8
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init

  !$ USE OMP_LIB

  IMPLICIT NONE

  TYPE(dbcsr_obj)                          :: matrix_a
  TYPE(dbcsr_error_type)                   :: error

  TYPE(array_i1d_obj)                      :: col_blk_sizes, row_blk_sizes
  INTEGER, DIMENSION(:), POINTER           :: rbs, cbs
  INTEGER                                  :: mp_comm, group, numnodes, mynode, &
       prow, pcol, nblkrows_total, nblkcols_total, node_holds_blk, max_nze, nze, &
       row, col, row_s, col_s, io_unit
  INTEGER, DIMENSION(2)                    :: npdims, myploc, coord
  INTEGER, DIMENSION(:,:), POINTER         :: pgrid
  TYPE(array_i1d_obj)                      :: col_dist, row_dist
  TYPE(dbcsr_distribution_obj)             :: dist
  TYPE(dbcsr_mp_obj)                       :: mp_env
  REAL(real_8), DIMENSION(:), ALLOCATABLE  :: values
  LOGICAL                                  :: tr

  !***************************************************************************************

  CALL mp_world_init(mp_comm)
  npdims(:) = 0
  CALL mp_cart_create (mp_comm, 2, npdims, myploc, group)
  CALL mp_environ (numnodes, mynode, group)
  ALLOCATE (pgrid(0:npdims(1)-1, 0:npdims(2)-1))
  DO prow = 0, npdims(1)-1
     DO pcol = 0, npdims(2)-1
        coord = (/ prow, pcol /)
        CALL mp_cart_rank (group, coord, pgrid(prow, pcol))
     ENDDO
  ENDDO
  CALL dbcsr_mp_new (mp_env, pgrid, group, mynode, numnodes,&
       myprow=myploc(1), mypcol=myploc(2))
  WRITE(*,*) 'mynode ',mynode,' numnodes',numnodes
  DEALLOCATE(pgrid)

  !***************************************************************************************
  !
  ! initialize libdbcsr
  CALL dbcsr_init_lib (group=mp_comm, error=error)
  CALL dbcsr_set_default_config (error=error)

  !
  ! set standard output parameters
  io_unit = 0
  IF (mynode.EQ.0) io_unit = default_output_unit

  !
  ! the matrix will contain nblkrows_total row blocks and nblkcols_total column blocks
  nblkrows_total = 4
  nblkcols_total = 4

  !
  ! set the block size for each row and column
  ALLOCATE (rbs(nblkrows_total), cbs(nblkcols_total))
  rbs(:) = 2
  cbs(:) = 2

  !
  ! set up the block size arrays
  CALL array_nullify (row_blk_sizes)
  CALL array_nullify (col_blk_sizes)
  CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
  CALL array_new (col_blk_sizes, cbs, gift=.TRUE.)

  !
  ! set the row and column distributions (here the distribution is set randomly)
  CALL random_dist (row_dist, nblkrows_total, npdims(1))
  CALL random_dist (col_dist, nblkcols_total, npdims(2))

  !
  ! set the dbcsr distribution object
  CALL dbcsr_distribution_new (dist, mp_env, row_dist, col_dist)

  !
  ! initialize the dbcsr matrix
  CALL dbcsr_init (matrix_a)

  !
  ! create the dbcsr matrix, i.e. a double precision non symmetric matrix
  ! with nblkrows_total x nblkcols_total blocks and
  ! sizes "sum(row_blk_sizes)" x "sum(col_blk_sizes)", distributed as
  ! specified by the dist object
  CALL dbcsr_create (matrix=matrix_a,&
       name="this is my matrix a",&
       dist=dist, &
       matrix_type=dbcsr_type_no_symmetry,&
       row_blk_size=row_blk_sizes,&
       col_blk_size=col_blk_sizes,&
       data_type=dbcsr_type_real_8,&
       error=error)

  !
  ! get the block sizes from the matrix
  rbs => array_data (dbcsr_row_block_sizes (matrix_a))
  cbs => array_data (dbcsr_col_block_sizes (matrix_a))

  !
  ! get the node id from the matrix
  mynode = dbcsr_mp_mynode (dbcsr_distribution_mp (dbcsr_distribution (matrix_a)))

  !
  ! get the maximum block size of the matrix
  max_nze = dbcsr_max_row_size (matrix_a) * dbcsr_max_col_size (matrix_a)

  !
  ! allocate a 1d buffer that is needed to put a block
  ! into the matrix (2d buffer can also be used)
  ALLOCATE (values (max_nze))

  !
  ! loop over the blocks, build a tridiagonal matrix
  DO row = 1, dbcsr_nblkrows_total (matrix_a)
     DO col = MAX(row-1,1),MIN(row+1,dbcsr_nblkcols_total (matrix_a))
        !
        ! get the node id that holds this (row, col) block
        tr = .FALSE.
        row_s=row ; col_s=col
        CALL dbcsr_get_stored_coordinates (matrix_a, row_s, col_s, tr, node_holds_blk)
        !
        ! put the block on the right node
        IF(node_holds_blk .EQ. mynode) THEN
           !
           ! get the size of the block
           nze = rbs(row_s) * cbs(col_s)
           !
           ! fill the matrix with the random block
           CALL RANDOM_NUMBER (values (1:nze))
           CALL dbcsr_put_block (matrix_a, row_s, col_s, values(1:nze))
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE (values)

  !
  ! finalize the dbcsr matrix
  CALL dbcsr_finalize(matrix_a, error=error)

  !
  ! print the matrix
  CALL dbcsr_print(matrix_a, error=error)

  !
  ! release the matrix
  CALL dbcsr_release (matrix_a)

  CALL dbcsr_distribution_release (dist)
  CALL dbcsr_mp_release (mp_env)
  CALL array_release (row_dist)
  CALL array_release (col_dist)
  CALL array_release (row_blk_sizes)
  CALL array_release (col_blk_sizes)

  !***************************************************************************************
  !
  ! finalize libdbcsr
  CALL dbcsr_finalize_lib (mp_comm, io_unit, error)

  !
  ! finalize mpi
  CALL mp_comm_free(group)
  CALL mp_world_finalize()

  !***************************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param dist_array ...
!> \param dist_size ...
!> \param nbins ...
! *****************************************************************************
  SUBROUTINE random_dist (dist_array, dist_size, nbins)
    TYPE(array_i1d_obj), INTENT(out)         :: dist_array
    INTEGER, INTENT(in)                      :: dist_size, nbins

    INTEGER                                  :: i
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: grid_dist

    ALLOCATE (grid_dist(dist_size))
    CALL array_nullify (dist_array)

    FORALL (i = 1 : dist_size)
       grid_dist(i) = MODULO (nbins-i, nbins)
    END FORALL

    CALL array_new (dist_array, grid_dist, lb=1)
    DEALLOCATE (grid_dist)

  END SUBROUTINE random_dist

END PROGRAM dbcsr_example_2
