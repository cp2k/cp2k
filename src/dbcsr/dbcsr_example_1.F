!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR example 1
!>          This example shows how to create a dbcsr matrix
!> \author  V. Weber
!> \date    2010
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
PROGRAM dbcsr_example_1
  USE dbcsr_api,                       ONLY: &
       array_i1d_obj, array_new, array_nullify, array_release, dbcsr_create, &
       dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_release, dbcsr_error_type, dbcsr_finalize, &
       dbcsr_finalize_lib, dbcsr_init, dbcsr_init_lib, dbcsr_mp_new, &
       dbcsr_mp_obj, dbcsr_mp_release, dbcsr_obj, dbcsr_print, dbcsr_release, &
       dbcsr_set_default_config, dbcsr_type_no_symmetry, dbcsr_type_real_8
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
       prow, pcol, nblkrows_total, nblkcols_total, io_unit
  INTEGER, DIMENSION(2)                    :: npdims, myploc, coord
  INTEGER, DIMENSION(:,:), POINTER         :: pgrid
  TYPE(array_i1d_obj)                      :: col_dist, row_dist
  TYPE(dbcsr_distribution_obj)             :: dist
  TYPE(dbcsr_mp_obj)                       :: mp_env

  !***************************************************************************************

  !
  ! initialize mpi
  CALL mp_world_init(mp_comm)

  !
  !
  npdims(:) = 0
  CALL mp_cart_create (mp_comm, 2, npdims, myploc, group)


  !
  !
  CALL mp_environ (numnodes, mynode, group)

  !
  !
  ALLOCATE (pgrid(0:npdims(1)-1, 0:npdims(2)-1))
  DO prow = 0, npdims(1)-1
     DO pcol = 0, npdims(2)-1
        coord = (/ prow, pcol /)
        CALL mp_cart_rank (group, coord, pgrid(prow, pcol))
     ENDDO
  ENDDO

  !
  ! setup the mp enviroment
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
  nblkcols_total = 3


  !
  ! set the block size for each row and column
  ALLOCATE (rbs(nblkrows_total), cbs(nblkcols_total))
  rbs(:) = 2
  cbs(:) = 3

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
  ! finalize the dbcsr matrix
  CALL dbcsr_finalize(matrix_a, error=error)

  !
  ! print the *empty* matrix
  CALL dbcsr_print(matrix_a, error=error)

  !
  ! release the matrix
  CALL dbcsr_release (matrix_a)

  !
  ! release the distribution
  CALL dbcsr_distribution_release (dist)

  !
  ! release the message passing environement
  CALL dbcsr_mp_release (mp_env)

  !
  ! release the block distributions
  CALL array_release (row_dist)
  CALL array_release (col_dist)

  !
  ! release the block size arrays
  CALL array_release (row_blk_sizes)
  CALL array_release (col_blk_sizes)


  !***************************************************************************************
  !
  ! finalize libdbcsr
  CALL dbcsr_finalize_lib (mp_comm, io_unit, error)

  ! free comm
  CALL mp_comm_free(group)

  !
  ! finalize mpi
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

END PROGRAM dbcsr_example_1
