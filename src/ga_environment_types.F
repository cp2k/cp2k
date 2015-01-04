!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module that collects GA functionality for semi-empirical
!>
!> \author CJM and Bruce Palmer (05.2011)
!> \par History
!>
! *****************************************************************************
MODULE ga_environment_types
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ga_environment_types'

#ifdef _USE_GA
  INTEGER, PRIVATE, SAVE :: last_ga_env_id_nr=0
#endif

  ! Data types
! *****************************************************************************
  TYPE ga_environment_type
    INTEGER :: se_dims ( 3 ), npairs, ncells
    INTEGER ::  ref_count, id_nr
    INTEGER ::  g_cnt, g_cellnum, g_offset, g_atoms, g_ks, g_d_info, g_pair_list
    INTEGER ::  g_count, pdx, pdy, pdz
    REAL ( dp ) :: rcoul
    REAL ( dp ) :: time_put, time_get, time_acc
    INTEGER :: iput, iget, iacc, ielm_get, ielm_acc, ielm_put
  END TYPE  ga_environment_type

! Public types

  PUBLIC :: ga_environment_type

! Public subroutine

  PUBLIC :: ga_env_release, ga_env_retain, init_ga_env, get_ga_env

CONTAINS
! *****************************************************************************
!> \brief Specific GA tool interfaces for the semi-empirical code
!> \param ga_env ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE init_ga_env ( ga_env, error )
    TYPE(ga_environment_type), POINTER       :: ga_env
    TYPE(cp_error_type), INTENT(inout)       :: error
#ifdef _USE_GA
    INTEGER :: stat
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'init_ga_env', &
      routineP = moduleN//':'//routineN


    IF (.NOT. ASSOCIATED ( ga_env ) ) THEN
      ALLOCATE (ga_env,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ga_env%pdx = 0
      ga_env%pdy = 0
      ga_env%pdz = 0
      ga_env%g_pair_list = 0
      ga_env%npairs = 0
      ga_env%ncells = 0
      ga_env%rcoul = 0.0_dp
      ga_env%g_cnt =0
      ga_env%g_cellnum=0
      ga_env%g_offset=0
      ga_env%g_atoms=0
      ga_env%g_ks=0
      ga_env%g_d_info=0
      ga_env%se_dims = 0
      ga_env%ref_count=1
      last_ga_env_id_nr=last_ga_env_id_nr+1
      ga_env%id_nr=last_ga_env_id_nr
    ENDIF
#endif
  END SUBROUTINE init_ga_env
! *****************************************************************************
!> \brief Specific GA tool interfaces for the semi-empirical code
!> \param ga_env ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE ga_env_release ( ga_env, error )
    TYPE ( ga_environment_type ), POINTER :: ga_env
    TYPE(cp_error_type), INTENT(inout)       :: error

#ifdef _USE_GA
    CHARACTER(len=*), PARAMETER :: routineN = 'ga_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
  failure = .FALSE.
  IF (ASSOCIATED(ga_env)) THEN
     CPPrecondition(ga_env%ref_count>0,cp_failure_level,routineP,error,failure)
     ga_env%ref_count=ga_env%ref_count-1
     IF (ga_env%ref_count<1) THEN
! deallocate
       DEALLOCATE(ga_env, stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(ga_env)
#endif
  END SUBROUTINE ga_env_release
! *****************************************************************************
!> \brief Specific GA tool interfaces for the semi-empirical code
!> \param ga_env ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE ga_env_retain ( ga_env, error )
    TYPE ( ga_environment_type ), POINTER :: ga_env
    TYPE(cp_error_type), INTENT(inout)       :: error


#ifdef _USE_GA
   CHARACTER(len=*), PARAMETER :: routineN = 'ga_env_retain', &
     routineP = moduleN//':'//routineN

   LOGICAL                                  :: failure

   failure=.FALSE.

   CPPrecondition(ASSOCIATED(ga_env),cp_failure_level,routineP,error,failure)
   IF (.NOT. failure) THEN
      CPPrecondition(ga_env%ref_count>0,cp_failure_level,routineP,error,failure)
      ga_env%ref_count=ga_env%ref_count+1
   END IF
#endif

  END SUBROUTINE ga_env_retain
! *****************************************************************************
!> \brief Specific GA tool interfaces for the semi-empirical code
!> \param ga_env ...
!> \param g_cnt ...
!> \param g_count ...
!> \param g_cellnum ...
!> \param g_offset ...
!> \param g_atoms ...
!> \param g_pair_list ...
!> \param g_ks ...
!> \param g_d_info ...
!> \param ncells ...
!> \param npairs ...
!> \param se_dims ...
!> \param rcoul ...
!> \param pdx ...
!> \param pdy ...
!> \param pdz ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_ga_env ( ga_env, g_cnt, g_count, g_cellnum, g_offset, g_atoms, &
                          g_pair_list, g_ks, g_d_info, ncells, npairs, se_dims, &
                          rcoul, pdx, pdy, pdz, error )
    TYPE ( ga_environment_type ), POINTER :: ga_env
    INTEGER, OPTIONAL :: g_count, g_cnt, g_cellnum, g_offset, g_atoms, g_ks, g_d_info
    INTEGER, OPTIONAL :: g_pair_list, ncells, npairs, pdx, pdy, pdz
    REAL ( dp ), OPTIONAL :: rcoul
    INTEGER, INTENT ( INOUT ), OPTIONAL :: se_dims ( 3 )
    TYPE(cp_error_type), INTENT(inout)       :: error
#ifdef _USE_GA
    INTEGER :: stat
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'get_ga_env', &
      routineP = moduleN//':'//routineN

    failure=.FALSE.
    CPPrecondition(ga_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(g_atoms)) g_atoms = ga_env % g_atoms
    IF (PRESENT(g_pair_list)) g_pair_list = ga_env % g_pair_list
    IF (PRESENT(g_cellnum)) g_cellnum = ga_env % g_cellnum
    IF (PRESENT(g_count)) g_count = ga_env % g_count
    IF (PRESENT(g_cnt)) g_cnt = ga_env % g_cnt
    IF (PRESENT(g_ks)) g_ks = ga_env % g_ks
    IF (PRESENT(g_d_info)) g_d_info = ga_env % g_d_info
    IF (PRESENT(se_dims)) se_dims = ga_env % se_dims
    IF (PRESENT(npairs)) npairs = ga_env % npairs
    IF (PRESENT(ncells)) ncells = ga_env % ncells
    IF (PRESENT(rcoul)) rcoul = ga_env % rcoul
    IF (PRESENT(pdx)) pdx = ga_env % pdx
    IF (PRESENT(pdy)) pdy = ga_env % pdy
    IF (PRESENT(rcoul)) pdz = ga_env % pdz
#endif
  END SUBROUTINE get_ga_env

END MODULE ga_environment_types
