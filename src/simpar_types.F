!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Type for storing MD parameters
!> \author CJM
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
! **************************************************************************************************
MODULE simpar_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief Simulation parameter type for molecular dynamics
!> \par History
!>         created [CJM]
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
!>         reorganization of the original routines/modules
! **************************************************************************************************
   TYPE simpar_type
      INTEGER       :: nsteps = -1
      INTEGER       :: max_steps = -1
      REAL(KIND=dp) :: dt = -1.0_dp
      REAL(KIND=dp) :: dt_fact = -1.0_dp
      REAL(KIND=dp) :: dr_tol = -1.0_dp
      REAL(KIND=dp) :: dsc_tol = -1.0_dp
      REAL(KIND=dp) :: temp_ext = -1.0_dp
      REAL(KIND=dp) :: temp_baro_ext = -1.0_dp
      REAL(KIND=dp) :: temp_baro = -1.0_dp
      REAL(KIND=dp) :: temp_tol = -1.0_dp
      REAL(KIND=dp) :: temp_baro_tol = -1.0_dp
      REAL(KIND=dp) :: p_ext = -1.0_dp
      REAL(KIND=dp) :: cmass = -1.0_dp
      REAL(KIND=dp) :: cmass_nph = -1.0_dp
      REAL(KIND=dp) :: v0 = -1.0_dp
      REAL(KIND=dp) :: e0 = -1.0_dp
      REAL(KIND=dp) :: v_shock = -1.0_dp
      REAL(KIND=dp) :: p0 = -1.0_dp
      REAL(KIND=dp) :: f_annealing = -1.0_dp
      REAL(KIND=dp) :: f_annealing_cell = -1.0_dp
      REAL(KIND=dp) :: f_temperature_annealing = -1.0_dp
      REAL(KIND=dp) :: gamma_nph = -1.0_dp
      INTEGER        :: ensemble = -1
      LOGICAL        :: constraint = .FALSE.
      LOGICAL        :: annealing = .FALSE.
      LOGICAL        :: annealing_cell = .FALSE.
      LOGICAL        :: temperature_annealing = .FALSE.
      LOGICAL        :: dump_lm = .FALSE.
      LOGICAL        :: angvel_zero = .FALSE.
      LOGICAL        :: variable_dt = .FALSE.
      INTEGER        :: nfree = -1, nfree_rot_transl = -1
      INTEGER        :: info_constraint = -1
      INTEGER        :: lagrange_multipliers = -1
      REAL(KIND=dp) :: tau_cell = -1.0_dp
      ! Constraints Parameters
      REAL(KIND=dp) :: shake_tol = -1.0_dp, roll_tol = -1.0_dp
      ! Langevin Parameters
      REAL(KIND=dp) :: gamma = -1.0_dp
      REAL(KIND=dp) :: noisy_gamma = -1.0_dp
      REAL(KIND=dp) :: shadow_gamma = -1.0_dp
      REAL(KIND=dp) :: var_w = -1.0_dp
      ! RESPA Parameters
      LOGICAL        :: multi_time_switch = .FALSE., do_respa = .FALSE.
      INTEGER        :: n_time_steps = -1
      ! SHELL parameters
      REAL(KIND=dp) :: temp_sh_ext = -1.0_dp
      REAL(KIND=dp) :: temp_sh_tol = -1.0_dp
      LOGICAL        :: temperature_per_kind = .FALSE.
      LOGICAL        :: scale_temperature_per_kind = .FALSE.
      LOGICAL        :: do_thermal_region = .FALSE.
      ! ADIABATIC parameters
      REAL(KIND=dp) :: temp_slow = -1.0_dp
      REAL(KIND=dp) :: temp_fast = -1.0_dp
      REAL(KIND=dp) :: temp_tol_fast = -1.0_dp, temp_tol_slow = -1.0_dp
      INTEGER :: n_resp_fast = -1
      ! Velocity softening Parameters
      INTEGER        :: soften_nsteps = -1
      REAL(KIND=dp) :: soften_alpha = -1.0_dp
      REAL(KIND=dp) :: soften_delta = -1.0_dp
      ! MD initialisation method
      INTEGER       :: initialization_method = -1
   END TYPE simpar_type

   PUBLIC :: simpar_type, &
             create_simpar_type, &
             release_simpar_type
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'simpar_types'

CONTAINS
! **************************************************************************************************
!> \brief Creates the simulation parameters type
!> \param simpar ...
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE create_simpar_type(simpar)
      TYPE(simpar_type), POINTER                         :: simpar

      CPASSERT(.NOT. ASSOCIATED(simpar))
      ALLOCATE (simpar)
   END SUBROUTINE create_simpar_type

! **************************************************************************************************
!> \brief Releases the simulation parameters type
!> \param simpar ...
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE release_simpar_type(simpar)
      TYPE(simpar_type), POINTER                         :: simpar

      CPASSERT(ASSOCIATED(simpar))
      DEALLOCATE (simpar)
   END SUBROUTINE release_simpar_type

END MODULE simpar_types
