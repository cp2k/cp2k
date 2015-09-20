!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Callback used by global geometry optimization schemes
!> \author Ole Schuett
! *****************************************************************************
MODULE glbopt_callback
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type,&
                                             pack_subsys_particles
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE kinds,                           ONLY: dp
  USE md_ener_types,                   ONLY: md_ener_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE mdctrl_types,                    ONLY: glbopt_mdctrl_data_type
#include "../common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_callback'

 PUBLIC :: glbopt_md_callback


 CONTAINS


! *****************************************************************************
!> \brief Callback used to hook into the main MD-loop.
!>        It recognizes and counts bumps in the potential energy.
!>        When MD_BUMPS_MAX is reached, the MD simulations is stoped.
!> \param mdctrl_data ...
!> \param md_env ...
!> \param should_stop ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE glbopt_md_callback(mdctrl_data, md_env, should_stop)
    TYPE(glbopt_mdctrl_data_type), POINTER   :: mdctrl_data
    TYPE(md_environment_type), POINTER       :: md_env
    LOGICAL, INTENT(inout)                   :: should_stop

    CHARACTER(len=*), PARAMETER :: routineN = 'glbopt_md_callback', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iw, n_atoms
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, passed_minimum
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: positions
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(md_ener_type), POINTER              :: md_ener

    failure = .FALSE.
    CPASSERT(ASSOCIATED(mdctrl_data))
    CPASSERT(ASSOCIATED(md_env))

    iw = mdctrl_data%output_unit

    ! add new potential energy value to history
    NULLIFY(md_ener, itimes)
    CALL get_md_env(md_env=md_env, md_ener=md_ener, itimes=itimes, force_env=force_env)
    mdctrl_data%itimes = itimes

    mdctrl_data%epot_history(:) = EOSHIFT(mdctrl_data%epot_history, shift=-1)
    mdctrl_data%epot_history(1) = md_ener%epot

    ! check if we passed a minimum
    passed_minimum = .TRUE.
    DO i=1, mdctrl_data%bump_steps_upwards
      IF(mdctrl_data%epot_history(i) <= mdctrl_data%epot_history(i+1)) &
         passed_minimum = .FALSE.
    END DO

    DO i=mdctrl_data%bump_steps_upwards+1, mdctrl_data%bump_steps_upwards+mdctrl_data%bump_steps_downwards
      IF(mdctrl_data%epot_history(i) >= mdctrl_data%epot_history(i+1)) &
         passed_minimum = .FALSE.
    END DO


    ! count the passed bumps and stop md_run when md_bumps_max is reached.
    IF(passed_minimum) &
       mdctrl_data%md_bump_counter = mdctrl_data%md_bump_counter + 1

    IF(mdctrl_data%md_bump_counter >= mdctrl_data%md_bumps_max) THEN
       should_stop = .TRUE.
       IF(iw>0) WRITE (iw,"(A)") " GLBOPT| Stopping MD because of MD_BUMPS_MAX."
    END IF

    CALL force_env_get(force_env, subsys=subsys)
    CALL cp_subsys_get(subsys, natom=n_atoms)
    ALLOCATE(positions(3*n_atoms))
    CALL pack_subsys_particles(subsys, r=positions)

  END SUBROUTINE glbopt_md_callback

END MODULE glbopt_callback

