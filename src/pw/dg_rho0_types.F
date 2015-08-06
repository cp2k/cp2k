!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
MODULE dg_rho0_types

  
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_zero
  USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE pw_types,                        ONLY: REALDATA3D,&
                                             pw_create,&
                                             pw_p_type,&
                                             pw_release
  USE termination,                     ONLY: stop_program
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: dg_rho0_type, dg_rho0_init, dg_rho0_set, dg_rho0_get, &
           dg_rho0_create, dg_rho0_retain, dg_rho0_release

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dg_rho0_types'
  INTEGER, PRIVATE, SAVE :: last_dg_rho0_id_nr=0

! *****************************************************************************
!> \brief   Type for Gaussian Densities
!>              type = type of gaussian (PME)
!>              grid = grid number
!>              gcc = Gaussian contraction coefficient
!>              zet = Gaussian exponent
! *****************************************************************************
  TYPE dg_rho0_type
     INTEGER :: ref_count, id_nr
     INTEGER :: TYPE
     INTEGER :: grid
     INTEGER :: kind
     REAL (KIND=dp) :: cutoff_radius
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: gcc
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: zet
     TYPE ( pw_p_type ) :: density
  END TYPE dg_rho0_type

CONTAINS

! *****************************************************************************
!> \brief   Set the dg_rho0_type
!> \param dg_rho0 ...
!> \param TYPE ...
!> \param grid ...
!> \param kind ...
!> \param cutoff_radius ...
!> \param gcc ...
!> \param zet ...
!> \param density ...
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dg_rho0_set ( dg_rho0, TYPE, grid, kind, cutoff_radius, &
       gcc, zet, density )
    INTEGER, OPTIONAL                        :: TYPE
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    INTEGER, OPTIONAL                        :: grid, kind
    REAL(KIND=dp), OPTIONAL                  :: cutoff_radius
    REAL(KIND=dp), OPTIONAL, POINTER         :: gcc( : ), zet( : )
    TYPE(pw_p_type), OPTIONAL                :: density

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_rho0_set', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( grid ) ) dg_rho0 % grid = grid
    IF ( PRESENT ( kind ) ) dg_rho0 % kind = kind
    IF ( PRESENT ( density ) ) dg_rho0 % density = density
    IF ( PRESENT ( gcc ) ) dg_rho0 % gcc => gcc
    IF ( PRESENT ( zet ) ) dg_rho0 % zet => zet
    IF ( PRESENT ( TYPE ) ) dg_rho0 % type = TYPE
    IF ( PRESENT ( cutoff_radius ) ) dg_rho0 % cutoff_radius = cutoff_radius

  END SUBROUTINE dg_rho0_set

! *****************************************************************************
!> \brief  Get the dg_rho0_type
!> \param dg_rho0 ...
!> \param cutoff_radius ...
!> \param TYPE ...
!> \param grid ...
!> \param kind ...
!> \param gcc ...
!> \param zet ...
!> \param density ...
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dg_rho0_get ( dg_rho0, cutoff_radius, TYPE, grid, kind, gcc, zet, density )
    INTEGER, OPTIONAL                        :: TYPE, kind, grid
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    REAL(KIND=dp), OPTIONAL                  :: cutoff_radius
    REAL(KIND=dp), OPTIONAL, POINTER         :: gcc( : ), zet( : )
    TYPE(pw_p_type), OPTIONAL, POINTER       :: density

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_rho0_get', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( grid ) ) grid = dg_rho0 % grid
    IF ( PRESENT ( kind ) ) kind = dg_rho0 % kind
    IF ( PRESENT ( density ) ) density = dg_rho0 % density
    IF ( PRESENT ( gcc ) ) gcc => dg_rho0 % gcc
    IF ( PRESENT ( zet ) ) zet => dg_rho0 % zet
    IF ( PRESENT ( TYPE ) ) TYPE = dg_rho0 % type
    IF ( PRESENT ( cutoff_radius ) ) cutoff_radius = dg_rho0 % cutoff_radius

  END SUBROUTINE dg_rho0_get

! *****************************************************************************
!> \brief   create the dg_rho0 structure
!> \param dg_rho0 ...
!> \param error ...
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dg_rho0_create ( dg_rho0, error )
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_rho0_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    ALLOCATE ( dg_rho0, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    NULLIFY ( dg_rho0 % gcc )
    NULLIFY ( dg_rho0 % zet )
    dg_rho0 % cutoff_radius = 0.0_dp
    dg_rho0 % grid = 0
    dg_rho0 % kind = 0
    dg_rho0 % type = do_ewald_none
    last_dg_rho0_id_nr=last_dg_rho0_id_nr+1
    dg_rho0%id_nr=last_dg_rho0_id_nr
    dg_rho0%ref_count=1
    NULLIFY ( dg_rho0%density%pw )

  END SUBROUTINE dg_rho0_create

! *****************************************************************************
!> \brief retains the given dg_rho0_type
!> \param dg_rho0 the dg_rho0_type to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! *****************************************************************************
  SUBROUTINE dg_rho0_retain ( dg_rho0, error)
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dg_rho0_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(dg_rho0),cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(dg_rho0%ref_count>0,cp_failure_level,routineP,error)
    dg_rho0%ref_count=dg_rho0%ref_count+1
  END SUBROUTINE dg_rho0_retain

! *****************************************************************************
!> \brief releases the given dg_rho0_type
!> \param dg_rho0 the dg_rho0_type to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! *****************************************************************************
  SUBROUTINE dg_rho0_release(dg_rho0, error)
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dg_rho0_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(dg_rho0)) THEN
       CPPreconditionNoFail(dg_rho0%ref_count>0,cp_failure_level,routineP,error)
       dg_rho0%ref_count=dg_rho0%ref_count-1
       IF (dg_rho0%ref_count==0) THEN
          IF ( ASSOCIATED ( dg_rho0 % gcc ) ) THEN
             DEALLOCATE ( dg_rho0 % gcc, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF ( ASSOCIATED ( dg_rho0 % zet ) ) THEN
             DEALLOCATE ( dg_rho0 % zet, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          CALL pw_release ( dg_rho0 % density % pw ,error=error)
          NULLIFY ( dg_rho0 % gcc )
          NULLIFY ( dg_rho0 % zet )
          DEALLOCATE (  dg_rho0 , stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(dg_rho0)
  END SUBROUTINE dg_rho0_release

! *****************************************************************************
!> \brief ...
!> \param dg_rho0 ...
!> \param pw_grid ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dg_rho0_init ( dg_rho0, pw_grid, error )
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_rho0_init', &
      routineP = moduleN//':'//routineN

    CALL pw_release ( dg_rho0 % density % pw ,error=error)
    SELECT CASE ( dg_rho0 % type )
    CASE ( do_ewald_ewald )
       CALL pw_create ( dg_rho0 % density % pw, pw_grid, REALDATA3D ,error=error)
       CALL dg_rho0_pme_gauss ( dg_rho0 % density, dg_rho0 % zet ( 1 ),error=error)
    CASE ( do_ewald_pme )
       CALL pw_create ( dg_rho0 % density % pw, pw_grid, REALDATA3D ,error=error)
       CALL dg_rho0_pme_gauss ( dg_rho0 % density, dg_rho0 % zet ( 1 ),error=error)
    CASE ( do_ewald_spme )
       CALL stop_program(routineN,moduleN,__LINE__,'SPME type not implemented')
    END SELECT

  END SUBROUTINE dg_rho0_init

! *****************************************************************************
!> \brief ...
!> \param dg_rho0 ...
!> \param alpha ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dg_rho0_pme_gauss ( dg_rho0, alpha, error)

    TYPE(pw_p_type), INTENT(INOUT)           :: dg_rho0
    REAL(KIND=dp), INTENT(IN)                :: alpha
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_rho0_pme_gauss', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: IMPOSSIBLE = 10000

    INTEGER                                  :: gpt, l0, ln, lp, m0, mn, mp, &
                                                n0, nn, np
    INTEGER, DIMENSION(:), POINTER           :: ghat
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(KIND=dp)                            :: const, e_gsq
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    TYPE(pw_grid_type), POINTER              :: pw_grid

    const = 1.0_dp / ( 8.0_dp * alpha ** 2 )

    pw_grid => dg_rho0 % pw % pw_grid
    bds => pw_grid % bounds

    IF ( -bds ( 1, 1 ) == bds ( 2, 1 ) ) THEN
       l0 = IMPOSSIBLE
    ELSE
       l0 = bds ( 1, 1 )
    END IF

    IF ( -bds ( 1, 2 ) == bds ( 2, 2 ) ) THEN
       m0 = IMPOSSIBLE
    ELSE
       m0 = bds ( 1, 2 )
    END IF

    IF ( -bds ( 1, 3 ) == bds ( 2, 3 ) ) THEN
       n0 = IMPOSSIBLE
    ELSE
       n0 = bds ( 1, 3 )
    END IF

    CALL pw_zero ( dg_rho0%pw,error=error)

    rho0 => dg_rho0 % pw % cr3d

    DO gpt = 1, pw_grid % ngpts_cut_local
       ghat => pw_grid % g_hat ( :, gpt )

       lp = pw_grid % mapl % pos ( ghat ( 1 ) )
       ln = pw_grid % mapl % neg ( ghat ( 1 ) )
       mp = pw_grid % mapm % pos ( ghat ( 2 ) )
       mn = pw_grid % mapm % neg ( ghat ( 2 ) )
       np = pw_grid % mapn % pos ( ghat ( 3 ) )
       nn = pw_grid % mapn % neg ( ghat ( 3 ) )

       e_gsq = EXP ( -const * pw_grid % gsq ( gpt ) ) / pw_grid % vol

       !*apsi
       lp = lp + bds ( 1, 1 )
       mp = mp + bds ( 1, 2 )
       np = np + bds ( 1, 3 )
       ln = ln + bds ( 1, 1 )
       mn = mn + bds ( 1, 2 )
       nn = nn + bds ( 1, 3 )

       rho0 ( lp, mp, np ) = e_gsq
       rho0 ( ln, mn, nn ) = e_gsq

       IF ( ghat ( 1 ) == l0 .OR. ghat ( 2 ) == m0 .OR. ghat ( 3 ) == n0 ) THEN
          rho0 ( lp, mp, np ) = 0.0_dp
          rho0 ( ln, mn, nn ) = 0.0_dp
       END IF

    END DO

  END SUBROUTINE dg_rho0_pme_gauss

END MODULE dg_rho0_types
