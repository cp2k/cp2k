!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types containing essential information for running implicit (iterative)
!>        Poisson solver
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Hossein Bani-Hashemian
! **************************************************************************************************
MODULE ps_implicit_types

   USE dct,                             ONLY: dct_type,&
                                              dct_type_release
   USE dielectric_types,                ONLY: dielectric_release,&
                                              dielectric_type
   USE dirichlet_bc_types,              ONLY: dbc_release,&
                                              dirichlet_bc_p_type
   USE kinds,                           ONLY: dp
   USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_release,&
                                              pw_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_implicit_types'

   TYPE :: ps_implicit_parameters
      INTEGER                           :: boundary_condition
      LOGICAL                           :: zero_initial_guess
      INTEGER                           :: max_iter
      REAL(KIND=dp)                     :: tol
      REAL(KIND=dp)                     :: omega
      INTEGER                           :: neumann_directions
   END TYPE ps_implicit_parameters

   TYPE :: ps_implicit_type
      LOGICAL                           :: do_dbc_cube
      TYPE(dielectric_type), POINTER    :: dielectric
      TYPE(dct_type)                    :: dct_env
      TYPE(pw_type), POINTER            :: initial_guess
      TYPE(pw_type), POINTER            :: v_eps
      TYPE(pw_type), POINTER            :: cstr_charge
      REAL(dp), DIMENSION(:), POINTER   :: initial_lambda => NULL()
      REAL(dp)                          :: ehartree
      REAL(dp)                          :: electric_enthalpy
      INTEGER                           :: times_called
      TYPE(dirichlet_bc_p_type), &
         DIMENSION(:), POINTER           :: contacts => NULL()
      REAL(dp), DIMENSION(:, :), POINTER :: QS => NULL()
      REAL(dp), DIMENSION(:, :), POINTER :: Rinv => NULL()
      REAL(dp), DIMENSION(:, :), POINTER :: B => NULL()
      REAL(dp), DIMENSION(:, :), POINTER :: Bt => NULL()
      REAL(dp), DIMENSION(:), POINTER   :: v_D => NULL()
      REAL(dp), DIMENSION(:), POINTER   :: osc_frac => NULL()
      REAL(dp), DIMENSION(:), POINTER   :: frequency => NULL()
      REAL(dp), DIMENSION(:), POINTER   :: phase => NULL()
      INTEGER, DIMENSION(:), POINTER    :: idx_1dto3d => NULL()
   END TYPE ps_implicit_type

   PUBLIC ps_implicit_type, ps_implicit_parameters
   PUBLIC ps_implicit_release

   INTEGER, PARAMETER, PUBLIC           :: PERIODIC_BC = 0, &
                                           MIXED_BC = 1, &
                                           MIXED_PERIODIC_BC = 2, &
                                           NEUMANN_BC = 3
CONTAINS

! **************************************************************************************************
!> \brief   Deallocates ps_implicit
!> \param ps_implicit_env the implicit_env to be deallocated
!> \param ps_implicit_params implicit env parameters
!> \param pw_pool pool of plane-wave grid
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE ps_implicit_release(ps_implicit_env, ps_implicit_params, pw_pool)
      TYPE(ps_implicit_type), POINTER                    :: ps_implicit_env
      TYPE(ps_implicit_parameters), INTENT(INOUT)        :: ps_implicit_params
      TYPE(pw_pool_type), INTENT(IN), OPTIONAL, POINTER  :: pw_pool

      CHARACTER(len=*), PARAMETER :: routineN = 'ps_implicit_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: boundary_condition, handle
      LOGICAL                                            :: can_give_back, do_dbc_cube

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(ps_implicit_env)) THEN
         can_give_back = PRESENT(pw_pool)
         IF (can_give_back) can_give_back = ASSOCIATED(pw_pool)
         do_dbc_cube = ps_implicit_env%do_dbc_cube

         IF (can_give_back) THEN
            CALL pw_pool_give_back_pw(pw_pool, ps_implicit_env%initial_guess, &
                                      accept_non_compatible=.TRUE.)
            CALL pw_pool_give_back_pw(pw_pool, ps_implicit_env%v_eps, &
                                      accept_non_compatible=.TRUE.)
            CALL pw_pool_give_back_pw(pw_pool, ps_implicit_env%cstr_charge, &
                                      accept_non_compatible=.TRUE.)
            CALL dbc_release(ps_implicit_env%contacts, do_dbc_cube, pw_pool=pw_pool)
         ELSE
            CALL pw_release(ps_implicit_env%initial_guess)
            CALL pw_release(ps_implicit_env%v_eps)
            CALL pw_release(ps_implicit_env%cstr_charge)
            CALL dbc_release(ps_implicit_env%contacts, do_dbc_cube)
         END IF

         DEALLOCATE (ps_implicit_env%initial_lambda)
         DEALLOCATE (ps_implicit_env%B)
         DEALLOCATE (ps_implicit_env%Bt)
         DEALLOCATE (ps_implicit_env%QS)
         DEALLOCATE (ps_implicit_env%Rinv)
         DEALLOCATE (ps_implicit_env%v_D)
         DEALLOCATE (ps_implicit_env%osc_frac)
         DEALLOCATE (ps_implicit_env%frequency)
         DEALLOCATE (ps_implicit_env%phase)
         DEALLOCATE (ps_implicit_env%idx_1dto3d)

         CALL dielectric_release(ps_implicit_env%dielectric, pw_pool)
         boundary_condition = ps_implicit_params%boundary_condition
         IF (boundary_condition .EQ. MIXED_BC .OR. boundary_condition .EQ. NEUMANN_BC) THEN
            CALL dct_type_release(ps_implicit_env%dct_env)
         END IF

         DEALLOCATE (ps_implicit_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE ps_implicit_release

END MODULE ps_implicit_types
