!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Tests for DBCSR operations
!> \author  VW
!> \date    2010
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2010
!> - 2015 Rewrote into standalone selftest (Ole Schuett)
! *****************************************************************************
PROGRAM dbcsr_selftest

  USE acc_device,                      ONLY: acc_get_ndevices,&
                                             acc_set_active_device
  USE dbcsr_config,                    ONLY: dbcsr_set_default_config
  USE dbcsr_error_handling,            ONLY: dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type
  USE dbcsr_lib,                       ONLY: dbcsr_finalize_lib,&
                                             dbcsr_init_lib
  USE dbcsr_mp_methods,                ONLY: dbcsr_mp_new,&
                                             dbcsr_mp_release
  USE dbcsr_test_add,                  ONLY: dbcsr_test_adds
  USE dbcsr_test_multiply,             ONLY: dbcsr_test_multiplies
  USE dbcsr_types,                     ONLY: dbcsr_mp_obj
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: default_output_unit
  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_cart_rank,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init

  !$ USE OMP_LIB

  IMPLICIT NONE


  INTEGER                                  :: mp_comm, group, numnodes, mynode, &
       prow, pcol, io_unit, error_handler
  INTEGER, DIMENSION(2)                    :: npdims, myploc
  INTEGER, DIMENSION(:,:), POINTER         :: pgrid
  TYPE(dbcsr_mp_obj)                       :: mp_env
  TYPE(dbcsr_error_type)                   :: error


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
  ! initialize libdbcsr
  IF (acc_get_ndevices() > 0) &
     CALL acc_set_active_device(MOD(mynode, acc_get_ndevices()))
  CALL dbcsr_init_lib (group=mp_comm, error=error)
  CALL dbcsr_set_default_config (error=error)


  ! run tests

  ! add -----------------------------------------------------------------------

  CALL dbcsr_test_adds("add_1",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 25/),&
       sparsities=(/0.7_dp, 0.5_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 1.0_dp, dp), beta=CMPLX(2.0_dp, 2.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2, 1, 3/),&
       limits=(/1, 50, 1, 25/), error=error)

  CALL dbcsr_test_adds("add_1",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50/),&
       sparsities=(/0.4_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(3.0_dp, 2.0_dp, dp), beta=CMPLX(4.0_dp, 0.5_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/),&
       limits=(/1, 50, 1, 50/), error=error)

  ! multiply ------------------------------------------------------------------

  CALL dbcsr_test_multiplies("multiply_ALPHA",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(-3.0_dp, -4.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 4/), bs_n=(/1, 4/), bs_k=(/1, 4/),&
       limits=(/2, 6, 3, 7, 6, 7/), error=error)

  CALL dbcsr_test_multiplies("multiply_BETA",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(3.0_dp, -2.0_dp, dp),&
       bs_m=(/1, 4/), bs_n=(/1, 4/), bs_k=(/1, 4/),&
       limits=(/2, 6, 3, 7, 6, 7/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_COL_1",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 1, 20, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_COL_2",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 9, 18, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_COL_3",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 9, 18, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_COL_4",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/25, 50, 75/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 25, 9, 18, 1, 75/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_K_1",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 1, 50, 1, 20/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_K_2",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 1, 50, 9, 18/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_K_3",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 1, 50, 9, 18/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_K_4",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/25, 50, 75/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 25, 1, 50, 9, 18/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_MIX_1",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/9, 18, 11, 20, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_MIX_2",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 9, 10, 11, 20/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_MIX_3",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/9, 20, 1, 50, 11, 18/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_MIX_4",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/11, 20, 11, 20, 13, 18/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_MIX_5",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/11, 20, 11, 20, 13, 18/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_MIX_6",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/25, 50, 75/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/11, 20, 11, 20, 13, 18/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_MIX_7",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/25, 50, 75/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 1.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2, 1, 2/), bs_n=(/1, 3, 1, 3, 1, 2/), bs_k=(/1, 0/),&
       limits=(/11, 20, 11, 20, 6, 10/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_ROW_1",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 20, 1, 50, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_ROW_2",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/9, 18, 1, 50, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_ROW_3",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/9, 18, 1, 50, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_LIMITS_ROW_4",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/25, 50, 75/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/9, 18, 1, 50, 1, 75/), error=error)

  CALL dbcsr_test_multiplies("multiply_RT",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 1, 50, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_SQ",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/50, 50, 50/),&
       sparsities=(/0.0_dp, 0.0_dp, 0.0_dp/), retain_sparsity=.FALSE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 50, 1, 50, 1, 50/), error=error)

  CALL dbcsr_test_multiplies("multiply_ALPHA",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(-3.0_dp, -4.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 4/), bs_n=(/1, 4/), bs_k=(/1, 4/),&
       limits=(/2, 6, 3, 7, 6, 7/), error=error)

  ! multiply-ghost ------------------------------------------------------------

  CALL dbcsr_test_multiplies("ub2",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 4/), bs_n=(/1, 4/), bs_k=(/1, 4/),&
       limits=(/2, 6, 3, 7, 6, 7/), error=error)

  CALL dbcsr_test_multiplies("ub-k-ghost",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 4/), bs_n=(/1, 4, 1, 4/), bs_k=(/1, 0/),&
       limits=(/2, 6, 3, 7, 2, 7/), error=error)

  CALL dbcsr_test_multiplies("ub-m-ghost",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 4/), bs_n=(/1, 0/), bs_k=(/1, 4/),&
       limits=(/2, 6, 3, 7, 2, 7/), error=error)

  CALL dbcsr_test_multiplies("ub-mnk-ghost",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 4, 1, 0/), bs_n=(/1, 4, 1, 0/), bs_k=(/1, 4/),&
       limits=(/2, 6, 3, 7, 2, 7/), error=error)

  CALL dbcsr_test_multiplies("ub-n-ghost",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 4, 1, 4/), bs_n=(/1, 0/), bs_k=(/1, 4/),&
       limits=(/2, 6, 3, 7, 2, 7/), error=error)

  CALL dbcsr_test_multiplies("ub",&
       group, mp_env, npdims, io_unit, matrix_sizes=(/20, 20, 20/),&
       sparsities=(/0.5_dp, 0.5_dp, 0.5_dp/), retain_sparsity=.TRUE.,&
       alpha=CMPLX(1.0_dp, 0.0_dp, dp), beta=CMPLX(0.0_dp, 0.0_dp, dp),&
       bs_m=(/1, 2/), bs_n=(/1, 2/), bs_k=(/1, 2/),&
       limits=(/1, 20, 1, 20, 9, 18/), error=error)


  ! end of test cases ---------------------------------------------------------

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

END PROGRAM dbcsr_selftest
