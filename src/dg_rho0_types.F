!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/dg_rho0_types [1.0] *
!!
!!   NAME
!!     dg_rho0_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE dg_rho0_types

  USE coefficient_types,               ONLY: PW_REALDATA3D,&
                                             coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_init,&
                                             coeff_type,&
                                             coeff_zero
  USE input_constants,                 ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: dg_rho0_type, dg_rho0_init, dg_rho0_set, dg_rho0_get, &
           dg_rho0_create, dg_rho0_retain, dg_rho0_release

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dg_rho0_types'
  INTEGER, PRIVATE, SAVE :: last_dg_rho0_id_nr=0

  ! type = type of gaussian (PME)
  ! grid = grid number
  ! gcc = Gaussian contraction coefficient
  ! zet = Gaussian exponent
  TYPE dg_rho0_type
     INTEGER :: ref_count, id_nr
     INTEGER :: TYPE
     INTEGER :: grid
     INTEGER :: kind
     REAL (KIND=dp) :: cutoff_radius
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: gcc
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: zet
     TYPE ( coeff_type ) :: density
  END TYPE dg_rho0_type

!!*****
!******************************************************************************

CONTAINS
!-----------------------------------------------------------------------------!
 SUBROUTINE dg_rho0_set ( dg_rho0, TYPE, grid, kind, cutoff_radius, & 
                          gcc, zet, density ) 

!   Purpose: Set the dg_rho0_type

!   ***************************************************************************

    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    INTEGER, OPTIONAL                        :: TYPE, grid, kind
    REAL(KIND=dp), OPTIONAL                  :: cutoff_radius
    REAL(KIND=dp), OPTIONAL, POINTER         :: gcc( : ), zet( : )
    TYPE(coeff_type), OPTIONAL               :: density

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_rho0_set', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( grid ) ) dg_rho0 % grid = grid
    IF ( PRESENT ( kind ) ) dg_rho0 % kind = kind
    IF ( PRESENT ( density ) ) dg_rho0 % density = density
    IF ( PRESENT ( gcc ) ) dg_rho0 % gcc => gcc
    IF ( PRESENT ( zet ) ) dg_rho0 % zet => zet
    IF ( PRESENT ( TYPE ) ) dg_rho0 % type = TYPE
    IF ( PRESENT ( cutoff_radius ) ) dg_rho0 % cutoff_radius = cutoff_radius

  END SUBROUTINE dg_rho0_set
!-----------------------------------------------------------------------------!
 SUBROUTINE dg_rho0_get ( dg_rho0, cutoff_radius, TYPE, grid, kind, gcc, zet, density ) 

!   Purpose: Get the dg_rho0_type

!   ***************************************************************************

    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    INTEGER, OPTIONAL                        :: kind, grid, TYPE
    REAL(KIND=dp), OPTIONAL                  :: cutoff_radius
    REAL(KIND=dp), OPTIONAL, POINTER         :: gcc( : ), zet( : )
    TYPE(coeff_type), OPTIONAL, POINTER      :: density

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_rho0_get', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( grid ) ) grid = dg_rho0 % grid
    IF ( PRESENT ( kind ) ) kind = dg_rho0 % kind
    IF ( PRESENT ( density ) ) density = dg_rho0 % density
    IF ( PRESENT ( gcc ) ) gcc => dg_rho0 % gcc
    IF ( PRESENT ( zet ) ) zet => dg_rho0 % zet
    IF ( PRESENT ( TYPE ) ) TYPE = dg_rho0 % type
    IF ( PRESENT ( cutoff_radius ) ) cutoff_radius = dg_rho0 % cutoff_radius

  END SUBROUTINE dg_rho0_get

! *****************************************************************************

  SUBROUTINE dg_rho0_create ( dg_rho0, error )

!   Purpose: create the dg_rho0 structure

!   ***************************************************************************

    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_rho0_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

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
    CALL coeff_init(dg_rho0%density,error=error)

  END SUBROUTINE dg_rho0_create
!***************************************************************************

!!****f* dg_rho0_types/dg_rho0_retain [1.0] *
!!
!!   NAME
!!     dg_rho0_retain
!!
!!   SYNOPSIS
!!     Subroutine dg_rho0_retain(dg, error)
!!       Type(dg_rho0_type), Pointer:: dg
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine dg_rho0_retain
!!
!!   FUNCTION
!!     retains the given dg_rho0_type
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - dg_rho0_type: the dg_rho0_type to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dg_rho0_retain ( dg_rho0, error)
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dg_rho0_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(dg_rho0),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(dg_rho0%ref_count>0,cp_failure_level,routineP,error)
     dg_rho0%ref_count=dg_rho0%ref_count+1
  END IF
END SUBROUTINE dg_rho0_retain
!***************************************************************************
!!****f* dg_rho0_types/dg_rho0_release [1.0] *
!!
!!   NAME
!!     dg_rho0_release
!!
!!   SYNOPSIS
!!     Subroutine dg_rho0_release(dg, error)
!!       Type(dg_rho0_type), Pointer:: dg
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine dg_rho0_release
!!
!!   FUNCTION
!!     releases the given dg_rho0_type
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - dg: the dg_rho0_type to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE dg_rho0_release(dg_rho0, error)
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
        CALL coeff_deallocate ( dg_rho0 % density )
        NULLIFY ( dg_rho0 % gcc )
        NULLIFY ( dg_rho0 % zet )
        DEALLOCATE (  dg_rho0 , stat = stat )
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(dg_rho0)
END SUBROUTINE dg_rho0_release

!******************************************************************************

SUBROUTINE dg_rho0_init ( dg_rho0, pw_grid )
    TYPE(dg_rho0_type), POINTER              :: dg_rho0
    TYPE(pw_grid_type), POINTER              :: pw_grid

  CALL coeff_deallocate ( dg_rho0 % density )
  SELECT CASE ( dg_rho0 % type )
  CASE ( do_ewald_ewald )
    CALL coeff_allocate ( dg_rho0 % density, pw_grid, PW_REALDATA3D )
    CALL dg_rho0_pme_gauss ( dg_rho0 % density, dg_rho0 % zet ( 1 ) )
  CASE ( do_ewald_pme )
    CALL coeff_allocate ( dg_rho0 % density, pw_grid, PW_REALDATA3D )
    CALL dg_rho0_pme_gauss ( dg_rho0 % density, dg_rho0 % zet ( 1 ) )
  CASE ( do_ewald_spme )
    CALL stop_program ( 'dg_rho0_init', 'SPME type not implemented' )
  END SELECT

END SUBROUTINE dg_rho0_init

!******************************************************************************

SUBROUTINE dg_rho0_pme_gauss ( dg_rho0, alpha )


    TYPE(coeff_type), INTENT(INOUT)          :: dg_rho0
    REAL(KIND=dp), INTENT(IN)                :: alpha

    INTEGER, PARAMETER                       :: IMPOSSIBLE = 10000

    INTEGER                                  :: gpt, l0, ln, lp, m0, mn, mp, &
                                                n0, nn, np
    INTEGER, DIMENSION(:), POINTER           :: ghat
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(KIND=dp)                            :: const, e_gsq
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    TYPE(pw_grid_type), POINTER              :: pw_grid

!------------------------------------------------------------------------------

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

  CALL coeff_zero ( dg_rho0 )

  rho0 => dg_rho0 % pw % cr3d

  DO gpt = 1, pw_grid % ngpts_cut
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
!******************************************************************************

END MODULE dg_rho0_types
