!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000,2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/green_types [1.0] *
!!
!!   NAME
!!     green_types
!!
!!   FUNCTION
!!     Routines to calculate the Green's function for the Poisson equation '''
!!
!!   AUTHOR
!!     JGH (19-Dec-2000) 
!! 
!!   MODIFICATION HISTORY
!!     JGH (9-Mar-2001) : include influence_function into greens_function_type
!!                        add cell volume as indicator for updates
!!     JGH (30-Mar-2001) : Added B-spline routines
!!
!!   SOURCE
!******************************************************************************

MODULE green_types
  USE bessel_lib,                      ONLY: bessj0,&
                                             bessj1,&
                                             bessk0,&
                                             bessk1
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             twopi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_release
  USE pw_types,                        ONLY: REALDATA1D,&
                                             RECIPROCALSPACE,&
                                             COMPLEXDATA1D,&
                                             pw_deallocate,&
                                             pw_release,&
                                             pw_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE input_constants,                 ONLY: periodic3d,& 
                                             analytic2d,& 
                                             analytic1d,& 
                                             analytic0d,& 
                                             mt2d,&       
                                             mt1d,&       
                                             mt0d,&       
                                             hockney2d,&
                                             hockney1d,&
                                             hockney0d 
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_green_fn,&
            greens_function_type,&
            pw_green_create,&
            pw_green_retain,&
            pw_green_release,&
            influence_factor
  
  TYPE greens_function_type
    INTEGER :: method
    INTEGER :: special_dimension
    INTEGER :: id_nr,  ref_count
    REAL (KIND=dp) :: radius
    REAL (KIND=dp) :: MT_alpha
    REAL (KIND=dp) :: MT_rel_cutoff
    REAL (KIND=dp) :: slab_size
    REAL (KIND=dp) :: alpha
    REAL (KIND=dp) :: volume
    LOGICAL :: p3m
    INTEGER :: p3m_order
    REAL (KIND=dp) :: p3m_alpha
    REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: p3m_coeff
    REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: p3m_bm2
    TYPE ( pw_type ), POINTER :: influence_function
    TYPE ( pw_type ), POINTER :: screen_function
    TYPE ( pw_type ), POINTER :: p3m_charge
    TYPE ( pw_grid_type), POINTER :: super_ref_pw_grid
  END TYPE greens_function_type


  INTEGER, SAVE, PRIVATE :: last_greens_fn_id_nr=0
  CHARACTER(len=*),PARAMETER,PRIVATE :: moduleN='green_types'
!!*****
!******************************************************************************

CONTAINS

  !******************************************************************************
  !!****** green_types/pw_green_retain [1.0] *
  !!
  !!   NAME
  !!     pw_green_retain
  !!
  !!   FUNCTION
  !!     retains the type 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pw_green_retain(gftype,error)
    ! Arguments
    TYPE(greens_function_type),pointer       :: gftype
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error    
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'pw_green_retain', &
         routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure
    
    failure=.FALSE.
 
    CPPrecondition(ASSOCIATED(gftype),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(gftype%ref_count>0,cp_failure_level,routineP,error,failure)
       gftype%ref_count=gftype%ref_count+1
    END IF
  END SUBROUTINE pw_green_retain

  !******************************************************************************
  !!****** green_types/pw_green_release [1.0] *
  !!
  !!   NAME
  !!     pw_green_release
  !!
  !!   FUNCTION
  !!     destroys the type (deallocates data)
  !!
  !!   AUTHOR
  !!     Joost VandeVondele
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !******************************************************************************

  SUBROUTINE pw_green_release ( gftype, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(greens_function_type),pointer       :: gftype
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error    
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'pw_green_release', &
      routineP = moduleN//':'//routineN    
    LOGICAL :: failure 
    INTEGER :: stat
    
    failure = .FALSE.
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(gftype)) THEN
          CPPrecondition(gftype%ref_count>0,cp_failure_level,routineP,error,failure)
          gftype%ref_count=gftype%ref_count-1
          IF (gftype%ref_count<1) THEN
             IF (ASSOCIATED(gftype%influence_function)) THEN
                CALL pw_deallocate(gftype%influence_function)
                DEALLOCATE(gftype%influence_function, stat=stat)
                CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
             ENDIF             
             CALL pw_release(gftype%screen_function)
             CALL pw_grid_release(gftype%super_ref_pw_grid)
             IF (gftype % p3m) THEN
                IF (ASSOCIATED(gftype % p3m_charge)) &
                     CALL pw_deallocate(gftype % p3m_charge)
                IF (ASSOCIATED(gftype % p3m_charge)) &
                     DEALLOCATE ( gftype % p3m_charge, stat=stat )
                CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
                IF (ASSOCIATED(gftype % p3m_bm2)) &
                     DEALLOCATE ( gftype % p3m_bm2 , stat=stat)
                CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
                IF (ASSOCIATED(gftype % p3m_coeff)) &
                     DEALLOCATE ( gftype % p3m_coeff , stat=stat)
                CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
             ENDIF
             DEALLOCATE(gftype, stat=stat)
             CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
          END IF
       END IF
    END IF
  END SUBROUTINE pw_green_release
  
  !!****f* green_fn/pw_green_create [1.0] *
  !!
  !!   NAME
  !!     pw_green_fn_init
  !!
  !!   FUNCTION
  !!     initializes the green function
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - green_fn: the object to initialize
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     fawzi
  !!     Teodoro Laino
  !!
  !!*** **********************************************************************
  SUBROUTINE pw_green_create(green,error)
    IMPLICIT NONE
    ! Arguments
    TYPE(greens_function_type), POINTER  :: green
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'pw_green_fn_init', &
         routineP = moduleN//':'//routineN
    LOGICAL :: failure
    INTEGER :: stat

    failure = .FALSE.
    IF (.NOT.failure) THEN       
       ALLOCATE(green, stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
       green%p3m=.FALSE.
       green%special_dimension = 0
       green%radius = 0.0_dp
       green%slab_size = 0.0_dp
       green%alpha = 0.0_dp
       green%volume = 0.0_dp
       green%method=PERIODIC3D
       last_greens_fn_id_nr = last_greens_fn_id_nr+1
       green % id_nr = last_greens_fn_id_nr
       green%ref_count=1
       green%MT_alpha=1.0_dp
       green%MT_rel_cutoff=1.0_dp
       green%p3m=.FALSE.
       green%p3m_order=0
       green%p3m_alpha=0.0_dp
       
       NULLIFY (green%influence_function,green%p3m_charge)
       NULLIFY (green%p3m_coeff,green%p3m_bm2)
       NULLIFY (green%screen_function)
       NULLIFY (green%super_ref_pw_grid)
    END IF
  END SUBROUTINE pw_green_create
  !***************************************************************************
  
  !!****** green_types/pw_green_fn [1.0] *
  !!
  !!   NAME
  !!     pw_green_fn
  !!
  !!   FUNCTION
  !!     Calculates the Green's function in reciprocal space  '''
  !!
  !!   AUTHOR
  !!     JGH (19-Dec-2000)
  !!
  !!   MODIFICATION HISTORY
  !!     JGH (9-Mar-2001) : update check through volume
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE pw_green_fn ( gftype )
    IMPLICIT NONE
    ! Arguments
    TYPE(greens_function_type), POINTER  :: gftype
    ! Local Variables
    INTEGER                                  :: ig, iz, method, nz
    REAL(KIND=dp)                            :: alpha, alpha2, g2, g3d, gg, &
         gxy, j0g, j1g, k0g, k1g, rlength, zlength
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_type), POINTER                   :: gf
    
    method = gftype % method
    gf   => gftype % influence_function
    grid => gftype % influence_function % pw_grid
    SELECT CASE ( method )
    CASE ( PERIODIC3D, ANALYTIC2D, ANALYTIC1D, ANALYTIC0D, MT2D, MT1D, MT0D )
       IF ( gf % in_space /= RECIPROCALSPACE .OR. &
            gf % in_use /= COMPLEXDATA1D ) &
            CALL stop_program ( "green_types", "gf wrong space or DATA TYPE" )
    CASE ( HOCKNEY2D, HOCKNEY1D, HOCKNEY0D )
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    CASE DEFAULT
       CALL stop_program ( "pw_green_fn", " Method not implemented ")
    END SELECT
    
    ! check if influence function is consistent with current box volume
    IF ( gftype % volume /= grid % vol ) THEN       
       gftype % volume = grid % vol         
       SELECT CASE ( method )            
       CASE ( PERIODIC3D )            
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             g2 = grid % gsq ( ig )
             gf % cc ( ig ) = fourpi / g2
          END DO
          IF ( grid % have_g0 ) gf % cc ( 1 ) = 0.0_dp
          
       CASE ( ANALYTIC2D )
          
          iz = gftype % special_dimension ! iz is the direction with NO PBC
          zlength = gftype % slab_size    ! zlength is the thickness of the cell
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             nz = grid % g_hat ( iz, ig )
             g2 = grid % gsq ( ig )
             g3d = fourpi / g2
             gg = 0.5_dp * SQRT ( g2 )
             gf % cc ( ig ) = g3d * ( 1.0_dp - (-1.0_dp)**nz * EXP ( - gg * zlength ) )
          END DO
          IF ( grid % have_g0 ) gf % cc ( 1 ) = 0.0_dp
          
       CASE ( ANALYTIC1D )
          
          ! iz is the direction of the PBC ( can be 1,2,3 -> x,y,z )
          iz = gftype % special_dimension
          ! rlength is the radius of the tube
          rlength = gftype % radius
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             g2 = grid % gsq ( ig )
             g3d = fourpi / g2
             gxy = SQRT ( g2 - grid % g(iz,ig) * grid % g(iz,ig) )
             j0g = bessj0 ( rlength * gxy )
             j1g = bessj1 ( rlength * gxy )
             k0g = bessk0 ( rlength * grid % g(iz,ig) )
             k1g = bessk1 ( rlength * grid % g(iz,ig) )
             gf % cc ( ig ) = g3d * ( 1.0_dp - rlength * &
                  ( gxy * j1g * k0g - grid % g(iz,ig) * j0g * k1g ) )
          END DO
          IF ( grid % have_g0 ) gf % cc ( 1 ) = 0.0_dp
          
       CASE ( ANALYTIC0D )
          
          rlength = gftype % radius   ! rlength is the radius of the sphere
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             g2 = grid % gsq ( ig )
             gg = 0.5_dp * SQRT ( g2 )
             g3d = fourpi / g2
             gf % cc ( ig ) = g3d * ( 1.0_dp - COS ( rlength * gg ) )
          END DO
          IF ( grid % have_g0 ) &
               gf % cc ( 1 ) = 0.5_dp * fourpi * rlength * rlength
          
       CASE ( MT2D, MT1D, MT0D )
          alpha   = gftype % MT_alpha   ! alpha is the convergence tuning parameter
          alpha2  = alpha * alpha
          DO ig = grid % first_gne0, grid % ngpts_cut_local
             g2 = grid % gsq ( ig ) 
             g3d = fourpi / g2
             gf%cc ( ig ) = g3d + gftype%screen_function%cc(ig) 
          END DO
          IF ( grid % have_g0 ) &
               gf%cc(1) = gftype%screen_function%cc(1) 
       CASE ( HOCKNEY2D )
          CALL stop_program ( "pw_green_fn", " Method not implemented ")
       CASE ( HOCKNEY1D )
          CALL stop_program ( "pw_green_fn", " Method not implemented ")
       CASE ( HOCKNEY0D )
          CALL stop_program ( "pw_green_fn", " Method not implemented ")
       CASE DEFAULT
          CALL stop_program ( "pw_green_fn", " Method not implemented ")
       END SELECT
       
    END IF
    
  END SUBROUTINE pw_green_fn
  
  !!****** green_types/influence_factor [1.0] *
  !!
  !!   NAME
  !!     influence_factor
  !!
  !!   FUNCTION
  !!     Calculates the influence_factor for the
  !!     SPME Green's function in reciprocal space'''
  !!
  !!   AUTHOR
  !!     DH (29-Mar-2001)
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !******************************************************************************    
  SUBROUTINE influence_factor ( gftype, error )
    IMPLICIT NONE
    ! Arguments
    TYPE(greens_function_type), POINTER  :: gftype
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    COMPLEX(KIND=dp)                             :: b_m, exp_m, sum_m
    INTEGER                                  :: dim, ierr, j, k, l, n, pt
    INTEGER, DIMENSION(3)                    :: npts
    INTEGER, DIMENSION(:), POINTER           :: lb, ub
    REAL(KIND=dp)                                :: l_arg, prod_arg, val
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: m_assign
    
    n = gftype % p3m_order
    
    ! calculate the assignment function values
    
    lb => gftype % influence_function % pw_grid % bounds (1, : )
    ub => gftype % influence_function % pw_grid % bounds (2, : )
    IF (ASSOCIATED(gftype % p3m_bm2)) THEN
       IF (LBOUND(gftype % p3m_bm2,2)/=MINVAL(lb).OR.&
            UBOUND(gftype % p3m_bm2,2)/=MAXVAL(ub)) THEN
          DEALLOCATE(gftype % p3m_bm2,stat=ierr)
          IF ( ierr /= 0 ) CALL stop_memory ( 'influence_factor, deallocating', &
               'gftype % p3m_bm2', 3*SIZE(gftype % p3m_bm2,2) )
       END IF
    END IF
    IF (.NOT.ASSOCIATED(gftype % p3m_bm2)) THEN
       ALLOCATE ( gftype % p3m_bm2 ( 3, MINVAL(lb):MAXVAL(ub) ), STAT = ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( 'influence_factor', &
            'gftype % p3m_bm2', 3*(MAXVAL(ub)-MINVAL(lb)+1) )
    END IF
    
    ALLOCATE ( m_assign ( 0:n-2 ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( 'influence_factor', 'm_assign', n-1 )
    m_assign = 0.0_dp
    DO k = 0, n-2
       j = -(n-1) + 2 * k
       DO l = 0, n-1
          l_arg = 0.5_dp ** l
          prod_arg = gftype % p3m_coeff ( j, l ) * l_arg
          m_assign ( k ) =  m_assign ( k ) + prod_arg
       END DO
    END DO
    
    ! calculate the absolute b values
    
    npts ( : ) = ub ( : ) - lb ( : ) + 1
    DO dim = 1, 3
       DO pt = lb (dim), ub (dim)
          val = twopi * ( REAL ( pt,KIND=dp) / REAL ( npts ( dim ),KIND=dp) )
          exp_m = CMPLX ( COS ( val ), -SIN ( val ),KIND=dp)
          sum_m = CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
          DO k = 0, n-2
             sum_m  =  sum_m + m_assign ( k ) * exp_m ** k
          END DO
          b_m = exp_m ** ( n - 1 ) / sum_m
          gftype % p3m_bm2 ( dim, pt ) = SQRT ( REAL ( b_m * CONJG ( b_m ),KIND=dp) )
       END DO
    END DO
    
    DEALLOCATE ( m_assign, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( 'influence_factor', 'm_assign', 0 )    
  END SUBROUTINE influence_factor
  
END MODULE green_types

!******************************************************************************
