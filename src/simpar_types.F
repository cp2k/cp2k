!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Type for storing MD parameters
!> \author CJM
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
! *****************************************************************************
MODULE simpar_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
!> \brief Simulation parameter type for molecular dynamics
!> \par History
!>         created [CJM]
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
! *****************************************************************************
  TYPE simpar_type
     INTEGER        :: nsteps
     REAL (KIND=dp) :: dt
     REAL (KIND=dp) :: dt_fact
     REAL (KIND=dp) :: dr_tol
     REAL (KIND=dp) :: dsc_tol
     REAL (KIND=dp) :: temp_ext
     REAL (KIND=dp) :: temp_baro_ext
     REAL (KIND=dp) :: temp_baro
     REAL (KIND=dp) :: temp_tol
     REAL (KIND=dp) :: temp_baro_tol
     REAL (KIND=dp) :: p_ext
     REAL (KIND=dp) :: cmass
     REAL (KIND=dp) :: cmass_nph
     REAL (KIND=dp) :: v0
     REAL (KIND=dp) :: e0
     REAL (KIND=dp) :: v_shock
     REAL (KIND=dp) :: p0
     REAL (KIND=dp) :: f_annealing
     REAL (KIND=dp) :: f_annealing_cell
     REAL (KIND=dp) :: gamma_nph
     INTEGER        :: ensemble
     LOGICAL        :: constraint
     LOGICAL        :: annealing
     LOGICAL        :: annealing_cell
     LOGICAL        :: dump_lm
     LOGICAL        :: angvel_zero
     LOGICAL        :: variable_dt
     INTEGER        :: nfree, nfree_rot_transl
     INTEGER        :: info_constraint
     INTEGER        :: lagrange_multipliers
     REAL (KIND=dp) :: tau_cell
     ! Constraints Parameters
     REAL (KIND=dp) :: shake_tol, roll_tol
     ! Langevin Parameters
     REAL (KIND=dp) :: gamma
     REAL (KIND=dp) :: noisy_gamma
     REAL (KIND=dp) :: shadow_gamma
     REAL (KIND=dp) :: var_w
     ! RESPA Parameters
     LOGICAL        :: multi_time_switch, do_respa
     INTEGER        :: n_time_steps
     ! SHELL parameters
     REAL (KIND=dp) :: temp_sh_ext
     REAL (KIND=dp) :: temp_sh_tol
     LOGICAL        :: temperature_per_kind
     LOGICAL        :: scale_temperature_per_kind
     LOGICAL        :: do_thermal_region
     ! ADIABATIC parameters
     REAL (KIND=dp) :: temp_slow
     REAL (KIND=dp) :: temp_fast
     REAL (KIND=dp) :: temp_tol_fast, temp_tol_slow
     INTEGER :: n_resp_fast
     ! Velocity softening Parameters
     INTEGER        :: soften_nsteps
     REAL (KIND=dp) :: soften_alpha
     REAL (KIND=dp) :: soften_delta
  END TYPE simpar_type


  PUBLIC :: simpar_type,&
            create_simpar_type,&
            release_simpar_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'simpar_types'

CONTAINS
! *****************************************************************************
!> \brief Creates the simulation parameters type
!> \param simpar ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE create_simpar_type(simpar, error)
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_simpar_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    ALLOCATE(simpar, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE create_simpar_type

! *****************************************************************************
!> \brief Releases the simulation parameters type
!> \param simpar ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE release_simpar_type(simpar, error)
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_simpar_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    DEALLOCATE(simpar, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE release_simpar_type

END MODULE simpar_types
