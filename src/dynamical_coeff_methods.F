!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief collection of routines handling the integrators for dynamical_coeff_type
!> \author Teo
! *****************************************************************************
MODULE dynamical_coeff_methods

  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_p_type,&
                                             dyn_coeff_set_type,&
                                             dyn_coeff_type,&
                                             get_dyn_coeff_set
  USE f77_blas
  USE input_constants,                 ONLY: dyn_coeff_mass_atom,&
                                             dyn_coeff_mass_coeff,&
                                             dyn_coeff_mass_global,&
                                             dyn_coeff_mass_kind
  USE kinds,                           ONLY: dp
  USE thermostat_methods,              ONLY: apply_thermostat_coeff
  USE thermostat_types,                ONLY: thermostat_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dynamical_coeff_methods'

  PUBLIC ::dynamical_coeff_verlet1, dynamical_coeff_verlet2

CONTAINS

! *****************************************************************************
!> \brief velocity verlet 1st part for dynamical_coeff_types
!>      it handles both cases of NVE or NPT over coeff_types
!> \par History
!>      Split from the integrator
!> \author Teodoro Laino - Zurich University - 02.2008
! *****************************************************************************
  SUBROUTINE dynamical_coeff_verlet1 ( dyn_coeff_set, thermostat_coeff, para_env, dt, error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(thermostat_type), POINTER           :: thermostat_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp)                            :: dt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dynamical_coeff_verlet1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icoef, iparticle_kind, &
                                                iparticle_local, mass_storage
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dmcoef, global_mass
    TYPE(dyn_coeff_p_type), DIMENSION(:), &
      POINTER                                :: coeffs_of_kind
    TYPE(dyn_coeff_type), POINTER            :: coeffs

    failure = .FALSE. 

    NULLIFY(coeffs_of_kind,coeffs)
    IF(ASSOCIATED(dyn_coeff_set)) THEN
       CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set,&
            coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
            mass_storage=mass_storage, error=error)
       CALL apply_thermostat_coeff(thermostat_coeff, dyn_coeff_set, para_env%group, error)
       IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
       DO iparticle_kind=1,SIZE(coeffs_of_kind)
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
             coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
             IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
             DO iparticle_local=1,coeffs%n_els
                DO icoef=1,coeffs%ncoef_atom
                   SELECT CASE (mass_storage)
                   CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind)
                      coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                           dmcoef*coeffs%forces(iparticle_local,icoef)
                   CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                      coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                           coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                   END SELECT
                   coeffs%pos(iparticle_local,icoef)=coeffs%pos(iparticle_local,icoef)+ &
                        dt*coeffs%vel(iparticle_local,icoef)
                END DO
             END DO
          END IF
       END DO
    END IF

  END SUBROUTINE dynamical_coeff_verlet1

! *****************************************************************************
!> \brief velocity verlet 2nd part for dynamical_coeff_types
!>      it handles both cases of NVE or NPT over coeff_types
!> \par History
!>      Split from the integrator
!> \author Teodoro Laino - Zurich University - 02.2008
! *****************************************************************************
  SUBROUTINE dynamical_coeff_verlet2 ( dyn_coeff_set, thermostat_coeff, para_env, dt, error)

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(thermostat_type), POINTER           :: thermostat_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp)                            :: dt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dynamical_coeff_verlet2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icoef, iparticle_kind, &
                                                iparticle_local, mass_storage
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dmcoef, global_mass
    TYPE(dyn_coeff_p_type), DIMENSION(:), &
      POINTER                                :: coeffs_of_kind
    TYPE(dyn_coeff_type), POINTER            :: coeffs

    failure = .FALSE. 
    NULLIFY(coeffs_of_kind,coeffs)
    IF(ASSOCIATED(dyn_coeff_set)) THEN
       CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set,&
            coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
            mass_storage=mass_storage,error=error)
       IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
       DO iparticle_kind=1,SIZE(coeffs_of_kind)
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
             coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
             IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
             DO iparticle_local=1,coeffs%n_els
                DO icoef=1,coeffs%ncoef_atom
                   SELECT CASE (mass_storage)
                   CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind)
                      coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                           dmcoef*coeffs%forces(iparticle_local,icoef)
                   CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                      coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                           coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                   END SELECT
                END DO
             END DO
          END IF
       END DO
       CALL apply_thermostat_coeff(thermostat_coeff, dyn_coeff_set, para_env%group, error)
    END IF

  END SUBROUTINE dynamical_coeff_verlet2

END MODULE dynamical_coeff_methods

