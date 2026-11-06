!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief State of the fictitious electronic degree of freedom (electron number q, its
!>        momentum, potentiostat and charge-mixing settings) used in constant-potential MD.
!> \author Nanchen DONGFANG
! **************************************************************************************************
MODULE qs_electrons_types
   USE kinds,                           ONLY: dp
   USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                              next_rng_seed,&
                                              rng_stream_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: qs_electrons_type, &
             qs_electrons_init, &
             qs_electrons_dealloc

   TYPE qs_electrons_type
      REAL(KIND=dp)            :: q = 0.0_dp                  ! electronic coordinate (total electron number)
      REAL(KIND=dp)            :: p = 0.0_dp                  ! conjugate momentum
      REAL(KIND=dp)            :: mass = 1.0_dp                ! fictitious electronic mass
      REAL(KIND=dp)            :: neutral_electrons = 0.0_dp
      REAL(KIND=dp)            :: target_workfunction = 0.0_dp
      REAL(KIND=dp)            :: target_fermi_level = 0.0_dp  ! eV, used when control_fermi is TRUE
      ! Control mode. FALSE = work-function
      ! control (drive Phi -> TARGET_WORKFUNCTION), TRUE = Fermi-level control
      ! (drive e_F -> TARGET_FERMI_LEVEL). The two are NOT equivalent: they differ
      ! by the vacuum level, and d(V_vac)/dq is not zero for a charged slab.
      LOGICAL                  :: control_fermi = .FALSE.
      ! Frozen mixing weight. When lambda_is_fixed is TRUE the mixing weight is
      ! held at lambda_fixed, q(t) is NOT propagated, and no endpoint is ever
      ! rebuilt: the run becomes ordinary MD on a fixed linear combination of the
      ! two integer-charge endpoints. Diagnostic use -- it isolates the mixing
      ! from the potentiostat, so a run with LAMBDA set contains no
      ! constant-potential dynamics at all.
      REAL(KIND=dp)            :: lambda_fixed = 0.0_dp
      LOGICAL                  :: lambda_is_fixed = .FALSE.
      ! Boundary-crossing policy. FALSE = tear down and re-run qs_init on the
      ! stale endpoint. TRUE = mutate the existing endpoint's electron count
      ! in place, never re-entering qs_init.
      LOGICAL                  :: charge_mixing_inplace = .FALSE.
      ! Report-only mode, selected by "&ELECTRONS 0" (section present, parameter
      ! false). The .cpcont file is written for an ordinary MD run -- in
      ! particular the Fermi level and work function are computed and logged
      ! every step -- but q is never propagated: no kick, no drift, no
      ! potentiostat, no charge mixing, no endpoint rebuilds. q stays at the
      ! system's actual electron count for the whole run and
      ! target_workfunction is 0, so the Legendre term vanishes and
      ! Grand_canonical_potential reduces to E_pot.
      LOGICAL                  :: report_only = .FALSE.
      REAL(KIND=dp)            :: tau_potentiostat = 0.0_dp    ! relaxation time, plain number in femtoseconds
      REAL(KIND=dp)            :: temp_potentiostat = 0.0_dp
      REAL(KIND=dp)            :: thermostat_energy = 0.0_dp   ! accumulated heat exchanged, for conserved-quantity bookkeeping
      REAL(KIND=dp)            :: current_workfunction = 0.0_dp
      REAL(KIND=dp)            :: current_fermi_level_eV = 0.0_dp
      INTEGER                  :: potential_average_axis = 3
      INTEGER                  :: cpcont_unit = -1
      REAL(KIND=dp)            :: potential_average_min = 0.0_dp
      REAL(KIND=dp)            :: potential_average_max = 0.0_dp
      INTEGER                  :: charge_mixing_span = 1       ! S: N_lo=S*floor(q/S), N_hi=N_lo+S
      INTEGER :: charge_mixing_current_n_lo = -HUGE(0)
      INTEGER :: charge_mixing_current_n_hi = -HUGE(0)
      LOGICAL                  :: charge_mixing_enabled = .FALSE.
      TYPE(rng_stream_type)    :: gaussian_rng_stream = rng_stream_type()
   END TYPE qs_electrons_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_electrons_types'

CONTAINS

! **************************************************************************************************
!> \brief Initialize the electronic coordinate/potentiostat state for constant-potential MD.
!> \param electrons ...
!> \param q_init ...
!> \param p_init ...
!> \param mass ...
!> \param neutral_electrons ...
!> \param target_workfunction ...
!> \param tau_potentiostat  relaxation time, a plain number in femtoseconds (converted to
!>        "number of thermostat calls" inside qs_electrons_thermostat_step)
!> \param temp_potentiostat ...
!> \param potential_average_axis axis (1=x, 2=y, 3=z) for the planar potential average
!> \param potential_average_min lower bound of the averaging slab, Angstrom
!> \param potential_average_max upper bound of the averaging slab, Angstrom
!> \param charge_mixing_span S: endpoint separation in electrons; presence of this
!>        argument is what enables charge mixing
!> \param target_fermi_level target Fermi level in eV (absolute, on CP2K's internal scale)
!> \param control_fermi TRUE for Fermi-level control, FALSE for work-function control
!> \param lambda_fixed frozen mixing weight, honoured only when lambda_is_fixed
!> \param lambda_is_fixed TRUE when LAMBDA was given: q is not propagated
!> \param charge_mixing_inplace TRUE to mutate an endpoint's electron count in place
!>        rather than re-running qs_init on it
!> \param report_only ...
! **************************************************************************************************
   SUBROUTINE qs_electrons_init(electrons, q_init, p_init, mass, neutral_electrons, &
                                target_workfunction, tau_potentiostat, temp_potentiostat, &
                                potential_average_axis, potential_average_min, potential_average_max, &
                                charge_mixing_span, target_fermi_level, control_fermi, &
                                lambda_fixed, lambda_is_fixed, charge_mixing_inplace, &
                                report_only)
      TYPE(qs_electrons_type), POINTER                   :: electrons
      REAL(KIND=dp), INTENT(IN)                          :: q_init, p_init, mass, neutral_electrons, &
                                                            target_workfunction, tau_potentiostat, &
                                                            temp_potentiostat
      INTEGER, INTENT(IN)                                :: potential_average_axis
      REAL(KIND=dp), INTENT(IN)                          :: potential_average_min, &
                                                            potential_average_max
      INTEGER, INTENT(IN), OPTIONAL                      :: charge_mixing_span
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: target_fermi_level
      LOGICAL, INTENT(IN), OPTIONAL                      :: control_fermi
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: lambda_fixed
      LOGICAL, INTENT(IN), OPTIONAL                      :: lambda_is_fixed, charge_mixing_inplace, &
                                                            report_only

      CPASSERT(.NOT. ASSOCIATED(electrons))
      ALLOCATE (electrons)
      electrons%q = q_init
      electrons%p = p_init
      electrons%mass = mass
      electrons%neutral_electrons = neutral_electrons
      electrons%target_workfunction = target_workfunction
      IF (PRESENT(target_fermi_level)) electrons%target_fermi_level = target_fermi_level
      IF (PRESENT(control_fermi)) electrons%control_fermi = control_fermi
      IF (PRESENT(lambda_fixed)) electrons%lambda_fixed = lambda_fixed
      IF (PRESENT(lambda_is_fixed)) electrons%lambda_is_fixed = lambda_is_fixed
      IF (PRESENT(charge_mixing_inplace)) electrons%charge_mixing_inplace = charge_mixing_inplace
      IF (PRESENT(report_only)) electrons%report_only = report_only
      electrons%tau_potentiostat = tau_potentiostat
      electrons%temp_potentiostat = temp_potentiostat
      electrons%thermostat_energy = 0.0_dp
      electrons%potential_average_axis = potential_average_axis
      electrons%potential_average_min = potential_average_min
      electrons%potential_average_max = potential_average_max
      ! Initialise the Gaussian random-number stream used by the CSVR potentiostat.
      ! This MUST be done explicitly, exactly as csvr_thermo_create does for the
      ! ionic thermostats.  Left at its default initialisation the component has
      ! distribution_type = UNIFORM and bg = cg = ig = 0, and rn32() then returns
      ! the constant m1*norm = 0.99999999976716947 on every call without ever
      ! advancing its state.  The linear noise term in rescaling_factor,
      !     2*rr*SQRT(kk*sigma/ndeg*(1-factor)*factor),
      ! consequently stops averaging to zero and the thermostat turns into a
      ! systematic energy source for the electronic degree of freedom.
      electrons%gaussian_rng_stream = rng_stream_type( &
                                      ! rng_name_length is 40 -- keep this string shorter
                                      name="Wiener process for e-potentiostat", &
                                      distribution_type=GAUSSIAN, &
                                      extended_precision=.TRUE., &
                                      seed=next_rng_seed())

      IF (PRESENT(charge_mixing_span)) THEN
         electrons%charge_mixing_span = charge_mixing_span
         electrons%charge_mixing_enabled = .TRUE.
         ! Left at the -HUGE(0) sentinel on purpose.  These record which
         ! endpoints the two FORCE_EVAL blocks are ACTUALLY built for, which is
         ! not knowable here -- md_run calls charge_mixing_sync_endpoint_state
         ! to read it back off the sub_force_envs before the first force
         ! evaluation.  Seeding them from q_init instead (the original code)
         ! ASSERTED that the input already agreed, so the first
         ! charge_mixing_maybe_rebuild saw no divergence and never corrected
         ! anything: an input whose two blocks were written high-charge-first
         ! mixed backwards, silently, for the entire run.
      END IF

   END SUBROUTINE qs_electrons_init

! **************************************************************************************************
!> \brief Release the electronic state.
!> \param electrons ...
! **************************************************************************************************
   SUBROUTINE qs_electrons_dealloc(electrons)
      TYPE(qs_electrons_type), POINTER                   :: electrons

      IF (ASSOCIATED(electrons)) DEALLOCATE (electrons)

   END SUBROUTINE qs_electrons_dealloc

END MODULE qs_electrons_types
