!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Prerformance for DBCSR operations
!> \author  VW
!> \date    2010
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010
! *****************************************************************************
PROGRAM dbcsr_performance_driver

  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_kinds,                     ONLY: default_string_length
  USE dbcsr_message_passing,           ONLY: mp_bcast,&
                                             mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init
  USE dbcsr_methods,                   ONLY: dbcsr_mp_new,&
                                             dbcsr_mp_release
  USE dbcsr_performance_multiply,      ONLY: dbcsr_perf_multiply
  USE dbcsr_test_methods,              ONLY: dbcsr_test_read_args
  USE dbcsr_types,                     ONLY: dbcsr_mp_obj

  !$ USE OMP_LIB

  IMPLICIT NONE


  INTEGER                                  :: mp_comm, group, numnodes, mynode, &
       prow, pcol, ionode, narg, error_handler
  INTEGER, DIMENSION(2)                    :: npdims, myploc
  INTEGER, DIMENSION(:,:), POINTER         :: pgrid
  TYPE(dbcsr_mp_obj)                       :: mp_env
  TYPE(dbcsr_error_type)                   :: error
  CHARACTER(len=default_string_length)     :: operation, args(100)


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
  !
  ALLOCATE (pgrid(0:npdims(1)-1, 0:npdims(2)-1))
  DO prow = 0, npdims(1)-1
     DO pcol = 0, npdims(2)-1
        CALL mp_cart_rank (group, (/ prow, pcol /), pgrid(prow, pcol))
     ENDDO
  ENDDO

  !
  ! setup the mp enviroment
  CALL dbcsr_mp_new (mp_env, pgrid, group, mynode, numnodes,&
       myprow=myploc(1), mypcol=myploc(2))
  DEALLOCATE(pgrid)

  !
  ! set parameters
  ionode = 0
  IF(mynode.EQ.0) ionode = 6

  !
  ! read args
  IF(mynode.eq.0) CALL dbcsr_test_read_args(narg, args)
  CALL mp_bcast(narg, 0, group)
  CALL mp_bcast(args, 0, group)

  CALL dbcsr_assert( narg, "GE", 1 , dbcsr_fatal_level, dbcsr_wrong_args_error, &
       routineN, "nargs not correct", __LINE__, error )
  !
  ! Select the operation
  SELECT CASE(args(1))
     !CASE()
  CASE('dbcsr_multiply')
     CALL dbcsr_perf_multiply( group, mp_env, npdims, ionode, narg, args, error=error)
  CASE DEFAULT
     CALL dbcsr_assert( .FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
          routineN, "operation not found", __LINE__, error )
  END SELECT

  !
  ! clean mp enviroment
  CALL dbcsr_mp_release (mp_env)

  !
  ! finalize mpi
  CALL mp_world_finalize()


  CALL dbcsr_error_stop(error_handler, error)

END PROGRAM dbcsr_performance_driver
