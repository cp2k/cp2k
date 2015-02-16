!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types containing essential information for running implicit (iterative)
!>        Poisson solver
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! *****************************************************************************
MODULE ps_implicit_types

  USE dct,                             ONLY: dct_type
  USE dielectric_types,                ONLY: dielectric_parameters,&
                                             dielectric_release,&
                                             dielectric_type
  USE dirichlet_bc_types,              ONLY: dirichlet_bc_p_type,&
                                             dirichlet_bc_parameters,&
                                             dirichlet_boundary_region_release
  USE kinds,                           ONLY: dp
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_release,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_implicit_types'

  TYPE :: ps_implicit_parameters
     INTEGER                           :: boundary_condition
     LOGICAL                           :: zero_initial_guess
     INTEGER                           :: max_iter
     REAL(KIND=dp)                     :: tol
     REAL(KIND=dp)                     :: omega
     TYPE(dielectric_parameters)       :: dielectric_params
     TYPE(dirichlet_bc_parameters)     :: dbc_params
  END TYPE ps_implicit_parameters

  TYPE :: ps_implicit_type
     TYPE(ps_implicit_parameters)      :: params
     TYPE(dielectric_type), POINTER    :: dielectric
     TYPE (dct_type)                   :: dct_env
     TYPE(pw_type), POINTER            :: initial_guess  
     REAL(dp), DIMENSION(:), POINTER   :: initial_lambda => NULL()
     REAL(dp)                          :: ehartree
     INTEGER                           :: times_called 
     TYPE(dirichlet_bc_p_type), & 
       DIMENSION(:), POINTER           :: gates => NULL()
     REAL(dp), DIMENSION(:,:), POINTER :: QS    => NULL()
     REAL(dp), DIMENSION(:,:), POINTER :: Rinv  => NULL()
     REAL(dp), DIMENSION(:,:), POINTER :: B   => NULL()
     REAL(dp), DIMENSION(:,:), POINTER :: Bt  => NULL()
     REAL(dp), DIMENSION(:), POINTER   :: v_D => NULL()
     INTEGER, DIMENSION(:), POINTER    :: idx_3dto1d => NULL()
  END TYPE ps_implicit_type

  PUBLIC ps_implicit_type, ps_implicit_parameters
  PUBLIC ps_implicit_release
 
  INTEGER, PARAMETER, PUBLIC           :: PERIODIC_BC       = 0,   &
                                          MIXED_BC          = 1,   &
                                          MIXED_PERIODIC_BC = 2,   &
                                          NEUMANN_BC        = 3
CONTAINS

! *****************************************************************************
!> \brief   Deallocates ps_implicit
!> \param implicit_env the implicit_env to be deallocated
!> \param pw_pool pool of plane-wave grid  
!> \param error cp2k error
!>
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE ps_implicit_release(implicit_env, pw_pool, error)
    TYPE(ps_implicit_type), POINTER          :: implicit_env
    TYPE(pw_pool_type), INTENT(IN), &
      OPTIONAL, POINTER                      :: pw_pool
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ps_implicit_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: can_give_back

    CALL timeset(routineN,handle)

    IF (ASSOCIATED(implicit_env)) THEN
       can_give_back = PRESENT(pw_pool)
       IF (can_give_back) can_give_back = ASSOCIATED(pw_pool)
       IF (can_give_back) THEN
          CALL pw_pool_give_back_pw(pw_pool, implicit_env%initial_guess, &
                                    accept_non_compatible=.TRUE., error=error)
          CALL dirichlet_boundary_region_release(implicit_env%gates, &
                                                 pw_pool=pw_pool, error=error)
       ELSE
          CALL pw_release(implicit_env%initial_guess, error=error)
          CALL dirichlet_boundary_region_release(implicit_env%gates, error=error)
       END IF

       DEALLOCATE(implicit_env%initial_lambda)
       DEALLOCATE(implicit_env%B)
       DEALLOCATE(implicit_env%Bt)
       DEALLOCATE(implicit_env%QS)
       DEALLOCATE(implicit_env%Rinv)
       DEALLOCATE(implicit_env%v_D)
       DEALLOCATE(implicit_env%idx_3dto1d)
    
       CALL dielectric_release(implicit_env%dielectric,pw_pool,error)

       DEALLOCATE(implicit_env)
    END IF

    CALL timestop(handle)

  END SUBROUTINE ps_implicit_release

END MODULE ps_implicit_types
