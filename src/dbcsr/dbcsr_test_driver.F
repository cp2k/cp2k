!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE dbcsr_config,                    ONLY: dbcsr_set_default_config
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_wrong_args_error
  USE dbcsr_methods,                   ONLY: dbcsr_mp_new,&
                                             dbcsr_mp_release
  USE dbcsr_operations,                ONLY: dbcsr_finalize_lib,&
                                             dbcsr_init_lib
  USE dbcsr_test_add,                  ONLY: dbcsr_test_adds
  USE dbcsr_test_methods,              ONLY: dbcsr_test_read_args
  USE dbcsr_test_multiply,             ONLY: dbcsr_test_multiplies
  USE dbcsr_types,                     ONLY: dbcsr_mp_obj
  USE kinds,                           ONLY: default_string_length
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init

  !$ USE OMP_LIB

  IMPLICIT NONE


  INTEGER                                  :: mp_comm, group, numnodes, mynode, &
       prow, pcol, io_unit, narg, error_handler
  INTEGER, DIMENSION(2)                    :: npdims, myploc
  INTEGER, DIMENSION(:,:), POINTER         :: pgrid
  TYPE(dbcsr_mp_obj)                       :: mp_env
  TYPE(dbcsr_error_type)                   :: error
  CHARACTER(len=default_string_length)     :: args(100)


  CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_check_multiply'


  !***************************************************************************************

  !
  ! initialize libdbcsr errors
  CALL dbcsr_error_set (routineN, error_handler, error)

  !
  ! initialize mpi
  CALL mp_world_init (mp_comm)

  ! setup the mp environment
  npdims(:) = 0
  CALL mp_cart_create (mp_comm, 2, npdims, myploc, group)
  CALL mp_environ (numnodes, mynode, group)
  ALLOCATE(pgrid(0:npdims(1)-1,0:npdims(2)-1))
  DO prow = 0, npdims(1)-1
     DO pcol = 0, npdims(2)-1
        CALL mp_cart_rank (group, (/ prow, pcol /), pgrid(prow, pcol))
     ENDDO
  ENDDO
  CALL dbcsr_mp_new (mp_env, pgrid, group, mynode, numnodes,&
       myprow=myploc(1), mypcol=myploc(2))
  DEALLOCATE(pgrid)

  !
  ! set standard output parameters
  io_unit = 0
  IF (mynode.EQ.0) io_unit = default_output_unit

  !
  ! read and distribute input args
  IF (mynode.eq.0) CALL dbcsr_test_read_args (narg, args)
  CALL mp_bcast (narg, 0, group)
  CALL mp_bcast (args, 0, group)
  CALL dbcsr_assert (narg.GE.1, dbcsr_fatal_level, dbcsr_wrong_args_error, &
     routineN, "nargs not correct", __LINE__, error)

  !
  ! initialize libdbcsr
  CALL dbcsr_init_lib (group=mp_comm, error=error)
  CALL dbcsr_set_default_config (error=error)

  !
  ! select the operation
  SELECT CASE(args(1))
    CASE('dbcsr_add')
       CALL dbcsr_test_adds (group, mp_env, npdims, io_unit, narg, args, error=error)
    CASE('dbcsr_multiply')
       CALL dbcsr_test_multiplies (group, mp_env, npdims, io_unit, narg, args, error=error)
    CASE DEFAULT
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error, &
          routineN, "operation not found", __LINE__, error)
  END SELECT

  !
  ! finalize libdbcsr
  CALL dbcsr_finalize_lib (mp_comm, io_unit, error)

  !
  ! clean mp enviroment
  CALL dbcsr_mp_release (mp_env)

  !
  ! finalize mpi
  CALL mp_comm_free(group)
  CALL mp_world_finalize ()

  !
  ! finalize libdbcsr errors
  CALL dbcsr_error_stop (error_handler, error)

END PROGRAM dbcsr_test_driver
