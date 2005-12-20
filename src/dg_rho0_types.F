!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE coefficient_types,               ONLY: coeff_type, &
                                             coeff_allocate, &
                                             coeff_deallocate, &
                                             coeff_zero, &
                                             coeff_init,&
                                             PW_REALDATA3D, &
                                             PW_COMPLEXDATA3D 
  USE kinds, ONLY: dp, default_string_length
  USE mathconstants, ONLY : pi
  USE pw_grid_types, ONLY : pw_grid_type
  USE termination, ONLY : stop_program
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme  

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: dg_rho0_type, dg_rho0_init, dg_rho0_set, dg_rho0_get, &
           dg_rho0_create, dg_rho0_retain, dg_rho0_release

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='dg_rho0_types'
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
 SUBROUTINE dg_rho0_set ( dg_rho0, type, grid, kind, cutoff_radius, & 
                          gcc, zet, density ) 

    IMPLICIT NONE
!   Purpose: Set the dg_rho0_type

!   ***************************************************************************

    TYPE ( dg_rho0_type ), POINTER :: dg_rho0
    INTEGER, OPTIONAL :: grid, kind
    REAL (KIND=dp), POINTER, OPTIONAL :: gcc ( : ), zet ( : )
    REAL (KIND=dp), OPTIONAL :: cutoff_radius
    TYPE ( coeff_type ), OPTIONAL :: density
    INTEGER, OPTIONAL :: TYPE
    
!   ---------------------------------------------------------------------------
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "dg_rho0_set",&
         routineP=moduleN//':'//routine_name

    IF ( PRESENT ( grid ) ) dg_rho0 % grid = grid
    IF ( PRESENT ( kind ) ) dg_rho0 % kind = kind
    IF ( PRESENT ( density ) ) dg_rho0 % density = density
    IF ( PRESENT ( gcc ) ) dg_rho0 % gcc => gcc
    IF ( PRESENT ( zet ) ) dg_rho0 % zet => zet
    IF ( PRESENT ( type ) ) dg_rho0 % type = type
    IF ( PRESENT ( cutoff_radius ) ) dg_rho0 % cutoff_radius = cutoff_radius

  END SUBROUTINE dg_rho0_set
!-----------------------------------------------------------------------------!
 SUBROUTINE dg_rho0_get ( dg_rho0, cutoff_radius, type, grid, kind, gcc, zet, density ) 

    IMPLICIT NONE
!   Purpose: Get the dg_rho0_type

!   ***************************************************************************

    TYPE ( dg_rho0_type ), POINTER :: dg_rho0
    INTEGER, OPTIONAL :: grid, kind
    REAL (KIND=dp), POINTER, OPTIONAL :: gcc ( : ), zet ( : )
    REAL (KIND=dp), OPTIONAL :: cutoff_radius
    TYPE ( coeff_type ), POINTER, OPTIONAL :: density
    INTEGER, OPTIONAL :: TYPE
    
!   ---------------------------------------------------------------------------
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "dg_rho0_get",&
         routineP=moduleN//':'//routine_name

    IF ( PRESENT ( grid ) ) grid = dg_rho0 % grid
    IF ( PRESENT ( kind ) ) kind = dg_rho0 % kind
    IF ( PRESENT ( density ) ) density = dg_rho0 % density
    IF ( PRESENT ( gcc ) ) gcc => dg_rho0 % gcc
    IF ( PRESENT ( zet ) ) zet => dg_rho0 % zet
    IF ( PRESENT ( type ) ) type = dg_rho0 % type
    IF ( PRESENT ( cutoff_radius ) ) cutoff_radius = dg_rho0 % cutoff_radius

  END SUBROUTINE dg_rho0_get

! *****************************************************************************

  SUBROUTINE dg_rho0_create ( dg_rho0, error )

!   Purpose: create the dg_rho0 structure

!   ***************************************************************************

    TYPE( dg_rho0_type ), POINTER :: dg_rho0
    TYPE ( cp_error_type ), INTENT ( INOUT ), OPTIONAL :: error
!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "dg_rho0_create",&
         routineP=moduleN//':'//routine_name
    INTEGER :: istat

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
    call coeff_init(dg_rho0%density,error=error)

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
  TYPE(dg_rho0_type), POINTER :: dg_rho0
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='dg_rho0_retain',&
        routineP=moduleN//':'//routineN

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
  TYPE(dg_rho0_type), POINTER :: dg_rho0
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='dg_rho0_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i,istat

  failure=.FALSE.
  
  IF (ASSOCIATED(dg_rho0)) THEN
     CPPreconditionNoFail(dg_rho0%ref_count>0,cp_failure_level,routineP,error)
     dg_rho0%ref_count=dg_rho0%ref_count-1
     IF (dg_rho0%ref_count==0) THEN
        IF ( ASSOCIATED ( dg_rho0 % gcc ) ) THEN 
           DEALLOCATE ( dg_rho0 % gcc, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF ( ASSOCIATED ( dg_rho0 % zet ) ) then
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
  IMPLICIT NONE
  ! Arguments
  TYPE ( dg_rho0_type ), POINTER :: dg_rho0
  TYPE ( pw_grid_type ), POINTER :: pw_grid

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

  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ), INTENT ( INOUT ) :: dg_rho0
  REAL (KIND=dp), INTENT ( IN ) :: alpha

! Locals
  TYPE ( pw_grid_type ), POINTER :: pw_grid
  INTEGER :: l0, m0, n0
  INTEGER, PARAMETER :: IMPOSSIBLE = 10000
  REAL (KIND=dp) :: const, e_gsq
  INTEGER :: lp, mp, np, ln, mn, nn, gpt
  INTEGER, DIMENSION ( :, : ), POINTER :: bds
  INTEGER, DIMENSION ( : ), POINTER :: ghat
  REAL (KIND=dp), DIMENSION ( :, :, : ), POINTER :: rho0

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
