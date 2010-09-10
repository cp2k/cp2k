!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Tests for DBCSR operations
!> \author  VW
!> \date    2010
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
PROGRAM dbcsr_test_driver

  USE array_types,                     ONLY: array_i1d_obj
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_kinds,                     ONLY: default_string_length,&
                                             real_8
  USE dbcsr_message_passing,           ONLY: mp_bcast,&
                                             mp_cart_create,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init
  USE dbcsr_test_methods,              ONLY: dbcsr_test_read_args
  USE dbcsr_test_multiply,             ONLY: dbcsr_test_multiplies
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_mp_obj

  !$ USE OMP_LIB

  IMPLICIT NONE


  TYPE(array_i1d_obj)                      :: col_blk_sizes, row_blk_sizes
  INTEGER, DIMENSION(:), POINTER           :: rbs, cbs, nproc, bs_m, bs_n, bs_k
  INTEGER                                  :: mp_comm, group, numnodes, mynode, &
       prow, pcol, nblkrows_total, nblkcols_total, ionode, narg, matrix_sizes(3), error_handler,&
       limits(6)
  INTEGER, DIMENSION(2)                    :: npdims, myploc
  INTEGER, DIMENSION(:,:), POINTER         :: pgrid
  REAL(real_8)                             :: sparsities(3), beta(2), alpha(2)
  TYPE(array_i1d_obj)                      :: col_dist, row_dist
  TYPE(dbcsr_distribution_obj)             :: dist
  TYPE(dbcsr_mp_obj)                       :: mp_env
  TYPE(dbcsr_error_type)                   :: error
  CHARACTER(len=default_string_length)     :: operation, args(100)
  CHARACTER :: symmetries(3)
  LOGICAL :: reatain_sparsity
  
  CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check_multiply'


  !***************************************************************************************

  CALL dbcsr_error_set(routineN, error_handler, error)

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
  ! set parameters
  ionode = 0
  IF(mynode.EQ.0) ionode = 6

  !
  ! read args
  IF(mynode.eq.0) CALL dbcsr_test_read_args(narg, args)
  CALL mp_bcast(narg, 0, group)
  CALL mp_bcast(args, 0, group)
  WRITE(*,*) 'narg',narg

  CALL dbcsr_assert( narg.GE.1 , dbcsr_fatal_level, dbcsr_wrong_args_error, &
       routineN, "nargs not correct", __LINE__, error )
  !
  ! Select the operation
  SELECT CASE(args(1))
     !CASE()
  CASE('dbcsr_multiply')
     CALL dbcsr_test_multiplies( group, ionode, nproc, narg, args, error=error)
  CASE DEFAULT
     CALL dbcsr_assert( .FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
          routineN, "operation not found", __LINE__, error )
  END SELECT

  !
  ! finalize mpi
  CALL mp_world_finalize()

  CALL dbcsr_error_stop(error_handler, error)

END PROGRAM dbcsr_test_driver
