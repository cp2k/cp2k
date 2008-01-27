
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/scp_rspw_methods [1.0] *
!!
!!   NAME
!!     scp_rspw_methods
!!
!!   FUNCTION
!!     set the SCP real-space plane-wave environment
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     CJM
!!
!!   SOURCE
!******************************************************************************

MODULE scp_rspw_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_util,                         ONLY: exp_radius
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE scp_rspw_types,                  ONLY: scp_rspw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_rspw_methods'
  PUBLIC ::  scp_qs_rspw_rebuild, scp_rspw_create, maxradius_scp_init

!-----------------------------------------------------------------------------!
CONTAINS
!-----------------------------------------------------------------------------!

!!*** **********************************************************************
SUBROUTINE scp_rspw_create(rspw,error)
    TYPE(scp_rspw_type), POINTER             :: rspw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_rspw_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(rspw, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY( rspw%qs_pw_env )
     rspw%ref_count=1
  END IF
END SUBROUTINE scp_rspw_create

!!****f* scp_density_calculate/maxradius_scp_init [1.0] *
!!
!!   NAME
!!     maxradius_scp_init
!!
!!   FUNCTION
!!     Initialize the maxradius to be used in cube file 
!!     initialization. Can be called from qs pw_env_methods
!!     or wherever needed (e.g. Fist, or scp_rspw routines )
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     04.2007 
!!
!!*** **********************************************************************  
SUBROUTINE maxradius_scp_init ( scp_env, atomic_kind_set, gridlevel_info, &
                              radius, eps_gvg_rspace, eps_rho_rspace, error)
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: radius
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: eps_gvg_rspace, eps_rho_rspace
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'maxradius_scp_init', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: safety_factor = 1.2

    INTEGER                                  :: igrid_level, ikind, ipgf, &
                                                iset, ishell, la, &
                                                lgrid_level, nkind, nseta
    INTEGER, DIMENSION(:), POINTER           :: npgfa, nshella
    INTEGER, DIMENSION(:, :), POINTER        :: lshella
    REAL(KIND=dp)                            :: maxradius, my_eps_gvg_rspace, &
                                                my_eps_rho_rspace, zetp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set

    IF ( PRESENT ( eps_gvg_rspace ) ) my_eps_gvg_rspace = eps_gvg_rspace
    IF ( PRESENT ( eps_rho_rspace ) ) my_eps_rho_rspace = eps_rho_rspace
    IF ( ( .NOT. PRESENT ( eps_gvg_rspace ) ).AND. (.NOT. PRESENT ( eps_rho_rspace ) ) ) THEN
      CALL get_scp_env ( scp_env, eps_rho_rspace = my_eps_rho_rspace, &
                         eps_gvg_rspace=my_eps_gvg_rspace, error = error )
    END IF
! Getting maxradius to initialize cubfiles for the rho_pol density defined
! in the AUX_BASIS_SET
      DO igrid_level=1, gridlevel_info % ngrid_levels
        maxradius=0.0_dp
        nkind = SIZE ( atomic_kind_set )
        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          
          CALL get_atomic_kind ( atomic_kind=atomic_kind,&
                                 orb_basis_set=aux_basis_set )
          CALL get_gto_basis_set ( gto_basis_set=aux_basis_set, &
                                  npgf=npgfa, nset=nseta, zet=zeta, &
                                  l=lshella,nshell=nshella)
           DO iset=1,nseta
              DO ipgf=1,npgfa(iset)
                 DO ishell=1,nshella(iset)
                   zetp = zeta(ipgf,iset)
                   la = lshella(ishell,iset)
                   lgrid_level = gaussian_gridlevel(gridlevel_info,zetp)
                   IF (lgrid_level .EQ. igrid_level) THEN
                     !density
                     maxradius=MAX(maxradius,exp_radius( la, zetp, &
                                  eps_rho_rspace, 1.0_dp))
                     !potential
                     maxradius=MAX(maxradius,exp_radius( la, zetp, &
                                  eps_gvg_rspace, 1.0_dp))
                   ENDIF
                 END DO
              END DO
           END DO
         END DO
        ! safety first, not very efficient, nor general
        ! one could possibly decide that this *is* the maximum allowed radius
         maxradius = maxradius * safety_factor
         IF ( maxradius > radius ( igrid_level ) ) radius(igrid_level)=maxradius
      END DO
  END SUBROUTINE maxradius_scp_init

!!*** **********************************************************************
!!**** scp_rspw_methods/scp_qs_rspw_rebuild        [1.0] *
!!
!!   NAME
!!     scp_qs_rspw_rebuild
!!
!!   FUNCTION
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - rspw: real-space/plane-wave env to be initialized
!!     - qs_env: the qs environment
!!     - error: error
!!
!!   AUTHOR
!!    CJM
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  SUBROUTINE scp_qs_rspw_rebuild( rspw, qs_env, error )

    TYPE(scp_rspw_type), POINTER             :: rspw
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_qs_rspw_rebuild', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(pw_env_type), POINTER               :: pw_env

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,'I','',handle)

    failure = .FALSE.
    NULLIFY ( pw_env )
    NULLIFY ( atomic_kind_set )

    CPPrecondition(ASSOCIATED ( rspw ),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED ( qs_env ),cp_failure_level,routineP,error,failure)
    CPPrecondition(rspw%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
!   Get the current pw_env from QS
       CALL get_qs_env ( qs_env = qs_env, pw_env = pw_env, &
                         atomic_kind_set = atomic_kind_set, &
                         error = error )
       CALL pw_env_retain ( pw_env, error ) 
       CALL pw_env_release ( rspw % qs_pw_env, error=error ) 
       rspw % qs_pw_env => pw_env
       CALL pw_env_release ( pw_env, error=error )
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE scp_qs_rspw_rebuild

! **************************************************************************

END MODULE scp_rspw_methods

!******************************************************************************
