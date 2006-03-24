!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
!!****** cp2k/ewald_pw_types [1.0] *
!!
!!   NAME
!!     ewald_pw_types
!!
!!   FUNCTION
!!     pw_types    
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE ewald_pw_types
  USE cell_types,                      ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
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
  USE input_constants,                 ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: MPI_COMM_SELF
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_create,&
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


  TYPE ewald_pw_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE ( pw_pool_type ), POINTER       :: pw_small_pool
     TYPE ( pw_pool_type ), POINTER       :: pw_big_pool
     TYPE ( cell_type ),    POINTER       :: cell, cell_ref
     TYPE ( pw_poisson_type ), POINTER    :: poisson_env
     TYPE ( dg_type ), POINTER            :: dg
  END TYPE ewald_pw_type
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!

!!*** **********************************************************************
!!****f* ewald_pw_types/ewald_pw_retain [1.0] *
!!
!!   NAME
!!     ewald_pw_retain
!!
!!   FUNCTION
!!     retains the structure ewald_pw_type
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
SUBROUTINE ewald_pw_retain(ewald_pw,error)
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
!!*** **********************************************************************
!!****f* ewald_pw_types/ewald_pw_create [1.0] *
!!
!!   NAME
!!     ewald_pw_create
!!
!!   FUNCTION
!!     creates the structure ewald_pw_type
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
SUBROUTINE ewald_pw_create ( ewald_pw, cell, cell_ref, ewald_env, force_env_section, error )
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
     NULLIFY ( ewald_pw % poisson_env )
     NULLIFY ( ewald_pw % cell ) 
     NULLIFY ( ewald_pw % cell_ref ) 
     CALL dg_create ( dg, error )
     ewald_pw % dg => dg 
     ewald_pw%ref_count=1
     last_ewald_pw_id_nr=last_ewald_pw_id_nr+1
     ewald_pw%id_nr=last_ewald_pw_id_nr
     CALL ewald_pw_init ( ewald_pw, cell, cell_ref, ewald_env, force_env_section,&
          error=error )
  END IF
END SUBROUTINE ewald_pw_create
!***************************************************************************
!****f* ewald_pw_types/ewald_pw_release [1.0] *
!!
!!   NAME
!!     release_ewald_pw
!!
!!   FUNCTION
!!     releases the memory used by the ewald_pw
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE ewald_pw_release ( ewald_pw, error )               


    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

   CALL timeset(routineN,"I","",handle)
   failure=.FALSE.  
   IF (ASSOCIATED(ewald_pw)) THEN
     CPPrecondition(ewald_pw%ref_count>0,cp_failure_level,routineP,error,failure)
     ewald_pw%ref_count=ewald_pw%ref_count-1
     IF (ewald_pw%ref_count<1) THEN
      CALL pw_pool_release ( ewald_pw % pw_small_pool, error = error )
      CALL pw_pool_release ( ewald_pw % pw_big_pool, error = error )
      CALL pw_poisson_release ( ewald_pw % poisson_env, error = error )
      CALL cell_release ( ewald_pw % cell, error )
      CALL cell_release ( ewald_pw % cell_ref, error )
      CALL dg_release ( ewald_pw % dg, error )
      DEALLOCATE (  ewald_pw , stat = stat )
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END IF
  NULLIFY ( ewald_pw )
  CALL timestop(0.0_dp,handle)

  END SUBROUTINE ewald_pw_release

!******************************************************************************
!****** ewald_pw_types/ewald_pw_init [1.0] *
!!
!!   NAME
!!     ewald_pw_init
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (12-Jan-2001): Added SPME part
!!     JGH (15-Mar-2001): Work newly distributed between initialize, setup,
!!                        and force routine
!!
!!   SOURCE
!!*****
!******************************************************************************
SUBROUTINE ewald_pw_init ( ewald_pw, cell, cell_ref, ewald_env, print_section, error )
  

    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_pw_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ewald_type, gmax, handle, &
                                                npts_s( 3 ), ns_max, &
                                                o_spline, output_unit, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, alphasq, cutoff, &
                                                cutoff_radius, ecut, epsilon, &
                                                i, norm
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_grid_type), POINTER              :: pw_big_grid, pw_small_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(section_vals_type), POINTER         :: poisson_section

!------------------------------------------------------------------------------

  CALL timeset("ewald_pw_init","I","",handle)
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
                         epsilon = epsilon )
    SELECT CASE ( ewald_type )
    CASE ( do_ewald_ewald )
       CALL pw_grid_create( pw_big_grid, MPI_COMM_SELF )
       ! set up EWALD "Classic"
       IF ( gmax == 2 * ( gmax / 2 ) ) THEN
          CALL stop_program ( "initialize_ewalds", "gmax has to be odd" )
       END IF
       pw_big_grid % bounds ( 1, : ) = -gmax / 2
       pw_big_grid % bounds ( 2, : ) = +gmax / 2
       npts_s = (/ gmax, gmax, gmax /)
       pw_big_grid % grid_span = HALFSPACE

       CALL pw_find_cutoff ( npts_s, cell_ref, cutoff )
       ecut = 0.5_dp * cutoff * cutoff
      ! CALL pw_find_cutoff ( pw_big_grid % bounds, cell_ref % h_inv, cutoff )
      ! ecut = cutoff 

       CALL pw_grid_setup ( cell_ref, pw_big_grid, print_section, cutoff = ecut, &
                           fft_usage = .FALSE. )
       NULLIFY ( pw_pool )
       CALL pw_pool_create ( pw_pool, pw_grid = pw_big_grid )
       ewald_pw % pw_big_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_big_pool )
       CALL pw_pool_release ( pw_pool )
       CALL pw_grid_release ( pw_big_grid, error )
    CASE ( do_ewald_pme )
       IF (.NOT.ASSOCIATED(ewald_pw%poisson_env)) THEN
          CALL pw_poisson_create(ewald_pw%poisson_env,error=error)
       END IF
       CALL pw_grid_create( pw_small_grid, MPI_COMM_SELF )
       CALL pw_grid_create( pw_big_grid, para_env%group )
       ! set up PME, the double grid method
       IF ( ns_max == 2 * ( ns_max / 2 ) ) THEN
          CALL stop_program ( "initialize_ewalds", "ns_max has to be odd" )
       END IF
       npts_s ( : ) = ns_max
! compute cut-off radius
       alphasq = alpha ** 2
       norm = ( 2.0_dp * alphasq / pi ) ** ( 1.5_dp )
       cutoff_radius = exp_radius ( 0, 2.0_dp * alphasq, epsilon , norm )

       CALL dg_pme_grid_setup ( cell_ref, npts_s, cutoff_radius,   &
                                pw_small_grid, pw_big_grid,   & 
                                para_env,rs_dims=(/para_env%num_pe,1/),&
                                print_section=print_section, &
                                error = error )
! Write some useful info
       logger => cp_error_get_logger ( error )
       output_unit = cp_print_key_unit_nr(logger,print_section,"",&
                                       extension=".Log",error=error)
       IF ( output_unit > 0 ) THEN
         WRITE ( output_unit, '( A,T71,E10.4 )' ) &
            ' EWALD| Gaussian tolerance (effective) ', epsilon
         WRITE ( output_unit, '( A,T63,3I6 )' ) &
            ' EWALD| Small box grid ', pw_small_grid % npts
         WRITE ( output_unit, '( A,T63,3I6 )' ) &
            ' EWALD| Full box grid ', pw_big_grid % npts
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,print_section,&
                                         "",error=error)  

     ! pw pools initialized
       NULLIFY ( pw_pool )
       CALL pw_pool_create( pw_pool, pw_grid = pw_big_grid )     
       ewald_pw % pw_big_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_big_pool )
       CALL pw_pool_release ( pw_pool )

       NULLIFY ( pw_pool )
       CALL pw_pool_create( pw_pool, pw_grid = pw_small_grid )     
       ewald_pw%pw_small_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_small_pool )
       CALL pw_pool_release ( pw_pool )
       CALL pw_grid_release ( pw_small_grid, error )
       CALL pw_grid_release ( pw_big_grid, error )
    CASE ( do_ewald_spme )
       IF (.NOT.ASSOCIATED(ewald_pw%poisson_env)) THEN
          CALL pw_poisson_create(ewald_pw%poisson_env,error=error)
       END IF
       CALL pw_grid_create( pw_big_grid, para_env%group )
     ! set up SPME
       npts_s = (/ gmax, gmax, gmax /)
       pw_big_grid % grid_span = HALFSPACE
       pw_big_grid % bounds ( 1, : ) = -gmax/2
       pw_big_grid % bounds ( 2, : ) = -gmax/2 + gmax - 1

       CALL pw_find_cutoff ( npts_s, cell_ref, cutoff )
       ecut = 0.5_dp * cutoff * cutoff
     !  CALL pw_find_cutoff ( pw_big_grid % bounds, cell_ref % h_inv, cutoff )
     !  ecut = cutoff 
       CALL pw_grid_setup( cell_ref, pw_big_grid, print_section, cutoff = ecut, &
            fft_usage = .TRUE., rs_dims=(/para_env%num_pe,1/))

     ! pw pools initialized
       NULLIFY (pw_pool)
       CALL pw_pool_create(pw_pool, pw_grid=pw_big_grid)     
       ewald_pw%pw_big_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_big_pool )
       CALL pw_pool_release ( pw_pool )
       CALL pw_grid_release ( pw_big_grid, error )
    CASE ( do_ewald_none )
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    IF (ASSOCIATED(ewald_pw%poisson_env)) THEN
       ALLOCATE(pw_pools(1),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       pw_pools(1)%pool => ewald_pw%pw_big_pool
       CALL pw_poisson_set(ewald_pw%poisson_env,cell=cell,parameters=poisson_section,&
            use_level=1,pw_pools=pw_pools,error=error)
!FM       CALL pw_poisson_rebuild(ewald_pw%poisson_env,error=error)
       DEALLOCATE(pw_pools,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL cell_retain ( cell )
    CALL cell_retain ( cell_ref )
    CALL cell_release ( ewald_pw % cell )
    CALL cell_release ( ewald_pw % cell_ref )
    ewald_pw % cell => cell
    ewald_pw % cell_ref => cell_ref
  ENDIF
  CALL timestop(0.0_dp,handle)

END SUBROUTINE ewald_pw_init

! *****************************************************************************
  SUBROUTINE ewald_pw_get (ewald_pw, cell, cell_ref, pw_big_pool,  &
                            pw_small_pool, poisson_env, dg )
!
!   Purpose: get the ewald_pw environment to the correct program.
!
! ***************************************************************************

    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_big_pool, pw_small_pool
    TYPE(pw_poisson_type), OPTIONAL, POINTER :: poisson_env
    TYPE(dg_type), OPTIONAL, POINTER         :: dg

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_pw_get', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( poisson_env ) ) poisson_env => ewald_pw % poisson_env
    IF ( PRESENT ( pw_big_pool ) )  pw_big_pool => ewald_pw % pw_big_pool
    IF ( PRESENT ( pw_small_pool ) )  pw_small_pool => ewald_pw % pw_small_pool
    IF ( PRESENT ( dg ) )  dg => ewald_pw % dg
    IF ( PRESENT ( cell ) )  cell => ewald_pw % cell
    IF ( PRESENT ( cell_ref ) )  cell_ref => ewald_pw % cell_ref
    
  END SUBROUTINE ewald_pw_get 
         
! *****************************************************************************
  SUBROUTINE ewald_pw_set (ewald_pw, cell, cell_ref, pw_big_pool,  &
                           pw_small_pool, dg, poisson_env, error )

!!   Purpose: Set the ewald_pw environment to the correct program.

!   ***************************************************************************

    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    TYPE(pw_pool_type), OPTIONAL, POINTER    :: pw_big_pool, pw_small_pool
    TYPE(dg_type), OPTIONAL, POINTER         :: dg
    TYPE(pw_poisson_type), OPTIONAL, POINTER :: poisson_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_pw_set', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( pw_big_pool ) )  THEN 
        CALL pw_pool_retain ( pw_big_pool )
        CALL pw_pool_release ( ewald_pw % pw_big_pool )
        ewald_pw % pw_big_pool => pw_big_pool
    ENDIF
    IF ( PRESENT ( pw_small_pool ) )  THEN 
        CALL pw_pool_retain ( pw_small_pool )
        CALL pw_pool_release ( ewald_pw % pw_small_pool )
        ewald_pw % pw_small_pool => pw_small_pool
    ENDIF
    IF ( PRESENT ( dg ) )  THEN 
        CALL dg_retain ( dg )
        CALL dg_release ( ewald_pw % dg, error = error )
        ewald_pw % dg => dg
    ENDIF
    IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(ewald_pw%cell,error=error)
       ewald_pw%cell => cell
    END IF
    IF (PRESENT(cell_ref)) THEN
       CALL cell_retain(cell_ref, error=error)
       CALL cell_release(ewald_pw%cell_ref,error=error)
       ewald_pw%cell_ref => cell_ref
    END IF
    IF (PRESENT(poisson_env)) THEN
       IF (ASSOCIATED(poisson_env)) &
            CALL pw_poisson_retain(poisson_env,error=error)
       CALL pw_poisson_release(ewald_pw%poisson_env,error=error)
       ewald_pw%poisson_env => poisson_env
    END IF
    
  END SUBROUTINE ewald_pw_set 
!******************************************************************************

END MODULE ewald_pw_types

!******************************************************************************
