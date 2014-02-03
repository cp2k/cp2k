!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief pw_types
!> \author CJM
! *****************************************************************************
MODULE ewald_pw_types
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dg_types,                        ONLY: dg_create,&
                                             dg_release,&
                                             dg_retain,&
                                             dg_type
  USE dgs,                             ONLY: dg_pme_grid_setup
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE f77_blas
  USE input_constants,                 ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: MPI_COMM_SELF
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_release,&
                                             pw_grid_setup
  USE pw_poisson_methods,              ONLY: pw_poisson_set
  USE pw_poisson_types,                ONLY: pw_poisson_create,&
                                             pw_poisson_release,&
                                             pw_poisson_retain,&
                                             pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE qs_util,                         ONLY: exp_radius
  USE realspace_grid_types,            ONLY: &
       init_input_type, realspace_grid_desc_type, realspace_grid_input_type, &
       realspace_grid_type, rs_grid_create, rs_grid_create_descriptor, &
       rs_grid_print, rs_grid_release, rs_grid_release_descriptor, &
       rs_grid_retain_descriptor
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_pw_types'
  INTEGER, PRIVATE, SAVE :: last_ewald_pw_id_nr=0
  PUBLIC :: ewald_pw_type, ewald_pw_release, &
            ewald_pw_retain, ewald_pw_create, &
            ewald_pw_get, ewald_pw_set

! *****************************************************************************
  TYPE ewald_pw_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE ( pw_pool_type ), POINTER       :: pw_small_pool
     TYPE ( pw_pool_type ), POINTER       :: pw_big_pool
     TYPE ( realspace_grid_desc_type ), POINTER    :: rs_desc
     TYPE ( pw_poisson_type ), POINTER    :: poisson_env
     TYPE ( dg_type ), POINTER            :: dg
  END TYPE ewald_pw_type

CONTAINS

! *****************************************************************************
!> \brief retains the structure ewald_pw_type
! *****************************************************************************
  SUBROUTINE ewald_pw_retain(ewald_pw,error)
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(ewald_pw),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(ewald_pw%ref_count>0,cp_failure_level,routineP,error,failure)
       ewald_pw%ref_count=ewald_pw%ref_count+1
    END IF
  END SUBROUTINE ewald_pw_retain

! *****************************************************************************
!> \brief creates the structure ewald_pw_type
! *****************************************************************************
  SUBROUTINE ewald_pw_create ( ewald_pw, ewald_env, cell, cell_ref, print_section, error )
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(dg_type), POINTER                   :: dg

    NULLIFY ( dg )
    failure=.FALSE.
    ALLOCATE ( ewald_pw, stat=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY ( ewald_pw % pw_big_pool )
       NULLIFY ( ewald_pw % pw_small_pool )
       NULLIFY ( ewald_pw % rs_desc )
       NULLIFY ( ewald_pw % poisson_env )
       CALL dg_create ( dg, error )
       ewald_pw % dg => dg
       ewald_pw%ref_count=1
       last_ewald_pw_id_nr=last_ewald_pw_id_nr+1
       ewald_pw%id_nr=last_ewald_pw_id_nr
       CALL ewald_pw_init ( ewald_pw, ewald_env, cell, cell_ref, print_section,&
            error=error )
    END IF
  END SUBROUTINE ewald_pw_create

!****f* ewald_pw_types/ewald_pw_release [1.0] *

! *****************************************************************************
!> \brief releases the memory used by the ewald_pw
! *****************************************************************************
  SUBROUTINE ewald_pw_release ( ewald_pw, error )
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.
    IF (ASSOCIATED(ewald_pw)) THEN
       CPPrecondition(ewald_pw%ref_count>0,cp_failure_level,routineP,error,failure)
       ewald_pw%ref_count=ewald_pw%ref_count-1
       IF (ewald_pw%ref_count<1) THEN
          CALL pw_pool_release ( ewald_pw % pw_small_pool, error = error )
          CALL pw_pool_release ( ewald_pw % pw_big_pool, error = error )
          CALL rs_grid_release_descriptor ( ewald_pw % rs_desc, error = error )
          CALL pw_poisson_release ( ewald_pw % poisson_env, error = error )
          CALL dg_release ( ewald_pw % dg, error )
          DEALLOCATE (  ewald_pw , stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY ( ewald_pw )
    CALL timestop(handle)
  END SUBROUTINE ewald_pw_release

!****** ewald_pw_types/ewald_pw_init [1.0] *

! *****************************************************************************
!> \par History
!>      JGH (12-Jan-2001): Added SPME part
!>      JGH (15-Mar-2001): Work newly distributed between initialize, setup,
!>                         and force routine
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_pw_init ( ewald_pw, ewald_env, cell, cell_ref, print_section, error )
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bo(2,3), ewald_type, gmax(3), &
                                                handle, npts_s( 3 ), ns_max, &
                                                o_spline, output_unit, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, alphasq, &
                                                cutoff_radius, epsilon, norm
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_grid_type), POINTER              :: pw_big_grid, pw_small_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_input_type)          :: input_settings
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(section_vals_type), POINTER         :: poisson_section, &
                                                rs_grid_section

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY ( pw_big_grid )
    NULLIFY ( pw_small_grid, poisson_section )

    CPPrecondition(ASSOCIATED ( ewald_pw ),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED ( ewald_env ),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED ( cell ),cp_failure_level,routineP,error,failure)
    CPPrecondition(ewald_pw%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( .NOT. failure ) THEN
       CALL ewald_env_get ( ewald_env = ewald_env,&
                            para_env = para_env, &
                            gmax = gmax, alpha = alpha,&
                            ns_max = ns_max,  &
                            ewald_type = ewald_type, &
                            o_spline=o_spline,&
                            poisson_section=poisson_section,&
                            epsilon = epsilon ,error=error)

       rs_grid_section => section_vals_get_subs_vals(poisson_section,"EWALD%RS_GRID",error=error)

       SELECT CASE ( ewald_type )
       CASE ( do_ewald_ewald )
          ! set up Classic EWALD sum
          logger => cp_error_get_logger ( error )
          output_unit = cp_print_key_unit_nr(logger,print_section,"",extension=".Log",error=error)
          CALL pw_grid_create( pw_big_grid, MPI_COMM_SELF ,error=error)

          IF ( ANY(gmax == 2 * ( gmax / 2 )) ) THEN
             CALL stop_program(routineN,moduleN,__LINE__,"gmax has to be odd.")
          END IF
          bo ( 1, : ) = -gmax / 2
          bo ( 2, : ) = +gmax / 2
          CALL pw_grid_setup ( cell_ref%hmat, pw_big_grid, grid_span=HALFSPACE, bounds=bo, spherical=.TRUE., &
                                fft_usage=.FALSE., iounit=output_unit, error=error)
          NULLIFY ( pw_pool )
          CALL pw_pool_create ( pw_pool, pw_grid = pw_big_grid ,error=error)
          ewald_pw % pw_big_pool => pw_pool
          CALL pw_pool_retain ( ewald_pw % pw_big_pool ,error=error)
          CALL pw_pool_release ( pw_pool ,error=error)
          CALL pw_grid_release ( pw_big_grid, error )
          CALL cp_print_key_finished_output(output_unit,logger,print_section,"",error=error)

       CASE ( do_ewald_pme )
          ! set up Particle-Mesh EWALD sum
          logger => cp_error_get_logger ( error )
          output_unit = cp_print_key_unit_nr(logger,print_section,"",extension=".Log",error=error)
          IF (.NOT.ASSOCIATED(ewald_pw%poisson_env)) THEN
             CALL pw_poisson_create(ewald_pw%poisson_env,error=error)
          END IF
          CALL pw_grid_create( pw_small_grid, MPI_COMM_SELF ,error=error)
          CALL pw_grid_create( pw_big_grid, para_env%group ,error=error)
          IF ( ns_max == 2 * ( ns_max / 2 ) ) THEN
             CALL stop_program(routineN,moduleN,__LINE__,"ns_max has to be odd.")
          END IF
          npts_s ( : ) = ns_max
          ! compute cut-off radius
          alphasq = alpha ** 2
          norm = ( 2.0_dp * alphasq / pi ) ** ( 1.5_dp )
          cutoff_radius = exp_radius ( 0, 2.0_dp * alphasq, epsilon , norm )

          CALL dg_pme_grid_setup ( cell_ref%hmat, npts_s, cutoff_radius,   &
               pw_small_grid, pw_big_grid, rs_dims=(/para_env%num_pe,1/),&
               iounit=output_unit, fft_usage=.TRUE., error=error )
          ! Write some useful info
          IF ( output_unit > 0 ) THEN
             WRITE ( output_unit, '( A,T71,E10.4 )' ) &
                  ' EWALD| Gaussian tolerance (effective) ', epsilon
             WRITE ( output_unit, '( A,T63,3I6 )' ) &
                  ' EWALD| Small box grid ', pw_small_grid % npts
             WRITE ( output_unit, '( A,T63,3I6 )' ) &
                  ' EWALD| Full box grid ', pw_big_grid % npts
          END IF

          ! pw pools initialized
          NULLIFY ( pw_pool )
          CALL pw_pool_create( pw_pool, pw_grid = pw_big_grid ,error=error)
          ewald_pw % pw_big_pool => pw_pool
          CALL pw_pool_retain ( ewald_pw % pw_big_pool ,error=error)
          CALL pw_pool_release ( pw_pool ,error=error)

          NULLIFY ( pw_pool )
          CALL pw_pool_create( pw_pool, pw_grid = pw_small_grid ,error=error)
          ewald_pw%pw_small_pool => pw_pool
          CALL pw_pool_retain ( ewald_pw % pw_small_pool ,error=error)
          CALL pw_pool_release ( pw_pool ,error=error)

          NULLIFY ( rs_desc )
          CALL init_input_type(input_settings,nsmax=MAXVAL ( pw_small_grid % npts(1:3) ),&
                               rs_grid_section=rs_grid_section,ilevel=1, &
                               higher_grid_layout=(/-1,-1,-1/),error=error)
          CALL rs_grid_create_descriptor(rs_desc,pw_big_grid, input_settings, error)

          CALL rs_grid_create(rs, rs_desc, error)
          CALL rs_grid_print(rs,output_unit,error=error)
          CALL rs_grid_release(rs, error)

          CALL cp_print_key_finished_output(output_unit,logger,print_section,"",error=error)

          ewald_pw%rs_desc => rs_desc

          CALL rs_grid_retain_descriptor ( ewald_pw % rs_desc, error)
          CALL rs_grid_release_descriptor ( rs_desc, error)

          CALL pw_grid_release ( pw_small_grid, error )
          CALL pw_grid_release ( pw_big_grid, error )

       CASE ( do_ewald_spme )
          ! set up the Smooth-Particle-Mesh EWALD sum
          logger => cp_error_get_logger ( error )
          output_unit = cp_print_key_unit_nr(logger,print_section,"",extension=".Log",error=error)
          IF (.NOT.ASSOCIATED(ewald_pw%poisson_env)) THEN
             CALL pw_poisson_create(ewald_pw%poisson_env,error=error)
          END IF
          CALL pw_grid_create( pw_big_grid, para_env%group ,error=error)
          npts_s = gmax
          CALL pw_grid_setup ( cell_ref%hmat, pw_big_grid, grid_span=HALFSPACE, npts=npts_s, spherical=.TRUE.,&
               rs_dims=(/para_env%num_pe,1/), iounit=output_unit, fft_usage=.TRUE., error=error)

          ! pw pools initialized
          NULLIFY ( pw_pool )
          CALL pw_pool_create(pw_pool, pw_grid=pw_big_grid,error=error)
          ewald_pw%pw_big_pool => pw_pool
          CALL pw_pool_retain ( ewald_pw % pw_big_pool ,error=error)
          CALL pw_pool_release ( pw_pool ,error=error)

          NULLIFY ( rs_desc )
          CALL init_input_type(input_settings,nsmax=o_spline,&
                               rs_grid_section=rs_grid_section,ilevel=1, &
                               higher_grid_layout=(/-1,-1,-1/),error=error)
          CALL rs_grid_create_descriptor(rs_desc, pw_big_grid, input_settings, error)

          CALL rs_grid_create(rs, rs_desc, error)
          CALL rs_grid_print(rs,output_unit,error=error)
          CALL rs_grid_release(rs, error)
          CALL cp_print_key_finished_output(output_unit,logger,print_section,"",error=error)

          ewald_pw%rs_desc => rs_desc

          CALL rs_grid_retain_descriptor ( ewald_pw % rs_desc, error)
          CALL rs_grid_release_descriptor ( rs_desc, error)

          CALL pw_grid_release ( pw_big_grid, error )
       CASE ( do_ewald_none )
          ! No EWALD sums..
       CASE default
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       ! Poisson Environment
       IF (ASSOCIATED(ewald_pw%poisson_env)) THEN
          ALLOCATE(pw_pools(1),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          pw_pools(1)%pool => ewald_pw%pw_big_pool
          CALL pw_poisson_set(ewald_pw%poisson_env,cell_hmat=cell%hmat,parameters=poisson_section,&
               use_level=1,pw_pools=pw_pools,error=error)
          DEALLOCATE(pw_pools,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    ENDIF
    CALL timestop(handle)
  END SUBROUTINE ewald_pw_init

!****** ewald_pw_types/ewald_pw_get [1.0] *

! *****************************************************************************
!> \brief get the ewald_pw environment to the correct program.
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_pw_get (ewald_pw, pw_big_pool, pw_small_pool, rs_desc, poisson_env, dg )

    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_big_pool, pw_small_pool
    TYPE(realspace_grid_desc_type), &
      OPTIONAL, POINTER                      :: rs_desc
    TYPE(pw_poisson_type), OPTIONAL, POINTER :: poisson_env
    TYPE(dg_type), OPTIONAL, POINTER         :: dg

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_pw_get', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( poisson_env ) ) poisson_env => ewald_pw % poisson_env
    IF ( PRESENT ( pw_big_pool ) )  pw_big_pool => ewald_pw % pw_big_pool
    IF ( PRESENT ( pw_small_pool ) )  pw_small_pool => ewald_pw % pw_small_pool
    IF ( PRESENT ( rs_desc ) )  rs_desc => ewald_pw % rs_desc
    IF ( PRESENT ( dg ) )  dg => ewald_pw % dg

  END SUBROUTINE ewald_pw_get

!****** ewald_pw_types/ewald_pw_set [1.0] *

! *****************************************************************************
!> \brief set the ewald_pw environment to the correct program.
!> \author CJM
! *****************************************************************************
  SUBROUTINE ewald_pw_set (ewald_pw, pw_big_pool, pw_small_pool, rs_desc, dg,&
       poisson_env, error )

    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_big_pool, pw_small_pool
    TYPE(realspace_grid_desc_type), &
      OPTIONAL, POINTER                      :: rs_desc
    TYPE(dg_type), OPTIONAL, POINTER         :: dg
    TYPE(pw_poisson_type), OPTIONAL, POINTER :: poisson_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_pw_set', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( pw_big_pool ) )  THEN
        CALL pw_pool_retain ( pw_big_pool ,error=error)
        CALL pw_pool_release ( ewald_pw % pw_big_pool ,error=error)
        ewald_pw % pw_big_pool => pw_big_pool
    ENDIF
    IF ( PRESENT ( pw_small_pool ) )  THEN
        CALL pw_pool_retain ( pw_small_pool ,error=error)
        CALL pw_pool_release ( ewald_pw % pw_small_pool ,error=error)
        ewald_pw % pw_small_pool => pw_small_pool
    ENDIF
    IF ( PRESENT ( rs_desc ) )  THEN
        CALL rs_grid_retain_descriptor (rs_desc, error=error)
        CALL rs_grid_release_descriptor (ewald_pw % rs_desc, error=error)
        ewald_pw % rs_desc => rs_desc
    ENDIF
    IF ( PRESENT ( dg ) )  THEN
        CALL dg_retain ( dg ,error=error)
        CALL dg_release ( ewald_pw % dg, error = error )
        ewald_pw % dg => dg
    ENDIF
    IF (PRESENT(poisson_env)) THEN
       IF (ASSOCIATED(poisson_env)) &
            CALL pw_poisson_retain(poisson_env,error=error)
       CALL pw_poisson_release(ewald_pw%poisson_env,error=error)
       ewald_pw%poisson_env => poisson_env
    END IF

  END SUBROUTINE ewald_pw_set

END MODULE ewald_pw_types
