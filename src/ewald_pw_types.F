!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C ) 2002 - 2003  CP2K developers group                         !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE message_passing, only: MPI_COMM_SELF
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dgs,                             ONLY: dg_pme_grid_setup
  USE dg_types,                        ONLY: dg_type, &
                                             dg_create, &
                                             dg_retain, &
                                             dg_release
  USE ewald_environment_types,         ONLY: ewald_environment_type, &
                                             ewald_env_get
  USE kinds,                           ONLY: dp, &
                                             default_string_length
  USE mathconstants,                   ONLY: zero, pi
  USE pw_grid_types,                   ONLY: pw_grid_type, &
                                             HALFSPACE
  USE pw_grids,                        ONLY: pw_grid_release, &
                                             pw_grid_setup, &
                                             pw_grid_construct, &
                                             pw_find_cutoff
  USE pw_pool_types,                   ONLY: pw_pool_release,&
                                             pw_pool_type, &
                                             pw_pool_create, &
                                             pw_pool_retain
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE qs_util,                         ONLY: exp_radius
  USE simulation_cell,                 ONLY: cell_type, &
                                             cell_retain, &
                                             cell_release
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE green_types,                     ONLY: greens_function_type,&
                                             pw_green_release,&
                                             pw_green_retain
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='ewald_pw_types'
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
     TYPE(greens_function_type), POINTER  :: green     
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
  TYPE(ewald_pw_type), POINTER :: ewald_pw
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
 
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='ewald_pw_retain',&
        routineP=moduleN//':'//routineN
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
SUBROUTINE ewald_pw_create ( ewald_pw, green, cell, cell_ref, ewald_env, error )
  TYPE(ewald_pw_type), POINTER :: ewald_pw
  TYPE(greens_function_type), POINTER  :: green
  TYPE ( cell_type ), POINTER :: cell, cell_ref
  TYPE(ewald_environment_type), POINTER :: ewald_env
  TYPE(cp_error_type), OPTIONAL, INTENT ( inout ) :: error

! Locals 
  TYPE ( dg_type ), POINTER :: dg
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='ewald_pw_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  NULLIFY ( dg )
  failure=.FALSE.
  ALLOCATE ( ewald_pw, stat=stat ) 
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY ( ewald_pw % pw_big_pool ) 
     NULLIFY ( ewald_pw % pw_small_pool ) 
     NULLIFY ( ewald_pw % green )
     NULLIFY ( ewald_pw % cell ) 
     NULLIFY ( ewald_pw % cell_ref ) 
     CALL dg_create ( dg, error )
     ewald_pw % dg => dg 
     ewald_pw%ref_count=1
     last_ewald_pw_id_nr=last_ewald_pw_id_nr+1
     ewald_pw%id_nr=last_ewald_pw_id_nr
     CALL ewald_pw_init ( ewald_pw, green, cell, cell_ref, ewald_env, error )
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

  IMPLICIT NONE

!arguments
    TYPE(ewald_pw_type), POINTER       :: ewald_pw
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)            :: error

!Locals
    CHARACTER(len=*), PARAMETER :: routineN='ewald_pw_release',&
       routineP=moduleN//':'//routineN

    LOGICAL :: failure
    INTEGER :: stat,handle

!   ---------------------------------------------------------------------------
 
   CALL timeset(routineN,"I","",handle)
   failure=.FALSE.  
   IF (ASSOCIATED(ewald_pw)) THEN
     CPPrecondition(ewald_pw%ref_count>0,cp_failure_level,routineP,error,failure)
     ewald_pw%ref_count=ewald_pw%ref_count-1
     IF (ewald_pw%ref_count<1) THEN
      CALL pw_pool_release ( ewald_pw % pw_small_pool, error = error )
      CALL pw_pool_release ( ewald_pw % pw_big_pool, error = error )
      CALL pw_green_release( ewald_pw % green, error=error)
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

SUBROUTINE ewald_pw_init ( ewald_pw, green, cell, cell_ref, ewald_env, error )
  
  IMPLICIT NONE

! Arguments
  TYPE ( ewald_environment_type ), POINTER :: ewald_env
  TYPE ( ewald_pw_type ), POINTER :: ewald_pw
  TYPE(greens_function_type), POINTER                   :: green
  TYPE ( cell_type ), POINTER :: cell, cell_ref
  TYPE ( cp_error_type ), OPTIONAL, INTENT ( inout ) :: error
  
! Locals
  INTEGER :: ns_max, o_spline, n, isos
  REAL (KIND=dp) :: alpha, epsilon
  TYPE ( pw_grid_type ), POINTER :: pw_big_grid
  TYPE ( pw_grid_type ), POINTER :: pw_small_grid
  INTEGER ::  iat, jat, gmax, npts_s ( 3 ), iw, output_unit, handle,stat
  REAL (KIND=dp) :: i, cutoff, ecut, alphasq, norm, cutoff_radius
  TYPE ( cp_para_env_type ), POINTER :: para_env
  LOGICAL :: pw_grid_information, failure, ionode
  TYPE ( cp_logger_type ), POINTER :: logger
  TYPE ( pw_pool_type ), POINTER :: pw_pool
  INTEGER :: ewald_type
  CHARACTER(len=*), PARAMETER :: routineN='ewald_pw_init',&
        routineP=moduleN//':'//routineN

!------------------------------------------------------------------------------
  CALL timeset("ewald_pw_init","I","",handle)
  failure=.FALSE.

  NULLIFY ( pw_big_grid )
  NULLIFY ( pw_small_grid )
  logger => cp_error_get_logger ( error )
  pw_grid_information=logger%print_keys%pw_grid_information
  ionode=(logger%para_env%mepos==logger%para_env%source)
  output_unit=-1
  IF (ionode) output_unit=cp_logger_get_default_unit_nr(logger)

  CPPrecondition(ASSOCIATED ( ewald_pw ),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED ( ewald_env ),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED ( green ),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED ( cell ),cp_failure_level,routineP,error,failure)
  CPPrecondition(ewald_pw%ref_count>0,cp_failure_level,routineP,error,failure)
  IF ( .NOT. failure ) THEN
    CALL ewald_env_get ( ewald_env = ewald_env,&
                         para_env = para_env, &  
                         gmax = gmax, alpha = alpha,&
                         ns_max = ns_max,  &
                         ewald_type = ewald_type, &
                         o_spline=o_spline,&
                         epsilon = epsilon )
    SELECT CASE ( ewald_type )
    CASE ( do_ewald_ewald )
       ALLOCATE(pw_big_grid,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pw_grid_construct( pw_big_grid, MPI_COMM_SELF )
       ! set up EWALD "Classic"
       IF ( gmax == 2 * ( gmax / 2 ) ) THEN
          IF ( ionode ) &
            CALL stop_program ( "initialize_ewalds", "gmax has to be odd" )
       END IF
       pw_big_grid % bounds ( 1, : ) = -gmax / 2
       pw_big_grid % bounds ( 2, : ) = +gmax / 2
       npts_s = (/ gmax, gmax, gmax /)
       pw_big_grid % grid_span = HALFSPACE

       CALL pw_find_cutoff ( npts_s, cell_ref, cutoff )

       ecut = 0.5_dp * cutoff * cutoff
       CALL pw_grid_setup ( cell_ref, pw_big_grid, cutoff = ecut, &
                           info = output_unit, fft_usage = .FALSE. )
       NULLIFY ( pw_pool )
       CALL pw_pool_create ( pw_pool, pw_grid = pw_big_grid )
       ewald_pw % pw_big_pool => pw_pool
       CALL pw_pool_retain ( ewald_pw % pw_big_pool )
       CALL pw_pool_release ( pw_pool )
       CALL pw_grid_release ( pw_small_grid, error )
    CASE ( do_ewald_pme )
       ALLOCATE(pw_big_grid,pw_small_grid,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pw_grid_construct( pw_small_grid, MPI_COMM_SELF )
       CALL pw_grid_construct( pw_big_grid, para_env%group )
       ! set up PME, the double grid method
       IF ( ns_max == 2 * ( ns_max / 2 ) ) THEN
          IF ( ionode ) &
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
                                error = error )
! Write some useful info
       IF ( ionode ) THEN
         WRITE ( output_unit, '( A,T71,E10.4 )' ) &
            ' EWALD| Gaussian tolerance (effective) ', epsilon
         WRITE ( output_unit, '( A,T63,3I6 )' ) &
            ' EWALD| Small box grid ', pw_small_grid % npts
         WRITE ( output_unit, '( A,T63,3I6 )' ) &
            ' EWALD| Full box grid ', pw_big_grid % npts
       END IF

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
       ALLOCATE(pw_big_grid,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pw_grid_construct( pw_big_grid, para_env%group )
     ! set up SPME
       npts_s = (/ gmax, gmax, gmax /)
       CALL pw_find_cutoff ( npts_s, cell_ref, cutoff )
       pw_big_grid % grid_span = HALFSPACE
       pw_big_grid % bounds ( 1, : ) = -gmax/2
       pw_big_grid % bounds ( 2, : ) = -gmax/2 + gmax - 1

       ecut = 0.5_dp * cutoff * cutoff

       CALL pw_grid_setup( cell_ref, pw_big_grid, cutoff = ecut, &
            info = output_unit, &
            fft_usage = .TRUE., rs_dims=(/para_env%num_pe,1/))

       n = o_spline
       ALLOCATE ( green % p3m_coeff ( -(n-1):n-1, 0:n-1 ), STAT = isos )
       IF ( isos /= 0 ) CALL stop_memory ( 'spme','coeff', 2*n*n )
       CALL spme_coeff_calculate ( n,  green % p3m_coeff )
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
    CALL cell_retain ( cell )
    CALL cell_retain ( cell_ref )
    CALL pw_green_retain ( green )
    CALL cell_release ( ewald_pw % cell )
    CALL cell_release ( ewald_pw % cell_ref )
    CALL pw_green_release ( ewald_pw % green )
    ewald_pw % cell => cell
    ewald_pw % cell_ref => cell_ref
    ewald_pw % green => green
  ENDIF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE ewald_pw_init

! *****************************************************************************
  SUBROUTINE ewald_pw_get (ewald_pw, cell, cell_ref, pw_big_pool,  &
                            pw_small_pool, green, dg )
!
!   Purpose: get the ewald_pw environment to the correct program.
!
! ***************************************************************************

    TYPE ( ewald_pw_type ), POINTER          :: ewald_pw
    TYPE ( pw_pool_type ), POINTER, OPTIONAL :: pw_big_pool
    TYPE ( pw_pool_type ), POINTER, OPTIONAL :: pw_small_pool
    TYPE ( dg_type ), POINTER, OPTIONAL      :: dg
    TYPE(greens_function_type), POINTER, OPTIONAL :: green
    TYPE ( cell_type ), POINTER, OPTIONAL    :: cell
    TYPE ( cell_type ), POINTER, OPTIONAL    :: cell_ref

    CHARACTER(LEN=*), PARAMETER :: routine_name = "ewald_pw_set", &
      routineP = moduleN//':'//routine_name

    IF ( PRESENT ( green ) ) green => ewald_pw % green
    IF ( PRESENT ( pw_big_pool ) )  pw_big_pool => ewald_pw % pw_big_pool
    IF ( PRESENT ( pw_small_pool ) )  pw_small_pool => ewald_pw % pw_small_pool
    IF ( PRESENT ( dg ) )  dg => ewald_pw % dg
    IF ( PRESENT ( cell ) )  cell => ewald_pw % cell
    IF ( PRESENT ( cell_ref ) )  cell_ref => ewald_pw % cell_ref
    
  END SUBROUTINE ewald_pw_get 
         
! *****************************************************************************
  SUBROUTINE ewald_pw_set (ewald_pw, cell, cell_ref, pw_big_pool,  &
                           pw_small_pool, dg, green, error )

!!   Purpose: Set the ewald_pw environment to the correct program.

!   ***************************************************************************

    TYPE ( ewald_pw_type ), POINTER          :: ewald_pw
    TYPE ( pw_pool_type ), POINTER, OPTIONAL :: pw_big_pool
    TYPE ( pw_pool_type ), POINTER, OPTIONAL :: pw_small_pool
    TYPE(greens_function_type), POINTER, OPTIONAL :: green
    TYPE ( dg_type ), POINTER, OPTIONAL      :: dg
    TYPE ( cell_type ), POINTER, OPTIONAL    :: cell
    TYPE ( cell_type ), POINTER, OPTIONAL    :: cell_ref
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "ewald_pw_set", &
      routineP = moduleN//':'//routine_name

    IF ( PRESENT (green ) ) THEN
       CALL pw_green_retain ( green )
       CALL pw_green_release ( ewald_pw % green )
       ewald_pw % green => green
    END IF
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
    
  END SUBROUTINE ewald_pw_set 
!******************************************************************************
!!****** ewald_pw_types/spme_coeff_calculate [1.0] *
!!
!!   NAME
!!     spme_coeff_calculate
!!
!!   FUNCTION
!!     Calculates the coefficients for the charge assignment function
!!
!!   AUTHOR
!!     DG (29-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE spme_coeff_calculate ( n, coeff )


    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp), &
      DIMENSION(-(n-1):n-1, 0:n-1), &
      INTENT(OUT)                            :: coeff

    INTEGER                                  :: i, j, l, m
    REAL(KIND=dp)                                :: b
    REAL(KIND=dp), DIMENSION(n, -n:n, 0:n-1)     :: a

!------------------------------------------------------------------------------

  a = 0.0_dp
  a ( 1, 0, 0 ) = 1.0_dp

  DO i = 2, n
     m = i-1
     DO j = -m, m, 2
        DO l = 0, m-1
           b = ( a ( m, j-1, l ) + &
               REAL ( (-1) ** l,KIND=dp) * a ( m, j+1, l ) ) / &
               REAL ( ( l + 1 ) * 2 ** ( l + 1 ) ,KIND=dp)
           a ( i, j, 0 ) = a ( i, j, 0 ) + b
        END DO
        DO l = 0, m-1
           a ( i, j, l+1 ) = ( a ( m, j+1, l ) - &
                               a ( m, j-1, l ) ) / REAL ( l + 1,KIND=dp)
        END DO
     END DO
  END DO

  coeff = 0.0_dp
  DO i = 0, n-1
    DO j = -(n-1), n-1, 2
      coeff ( j, i ) = a ( n, j, i )
    END DO
  END DO

END SUBROUTINE spme_coeff_calculate

!******************************************************************************

END MODULE ewald_pw_types

!******************************************************************************
