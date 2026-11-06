!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Dual integer-electron-endpoint charge mixing for constant-potential MD,
!>        built on top of CP2K's stock MIXED force_eval
!>        (METHOD MIXED / MIXING_TYPE LINEAR_COMBINATION). Energy, forces and
!>        virial are already interpolated by stock force_env_methods.F's
!>        mixed_energy_forces using the LAMBDA value read from the input
!>        structure at every call; this module (a) pushes the dynamically
!>        tracked lambda into that input structure each step, (b) mixes
!>        Fermi level / workfunction, which stock MIXED does not touch, and
!>        (c) rebuilds a sub_force_env at a new integer CHARGE whenever q(t)
!>        drifts outside the currently-built bracket ("boundary crossing").
!> \author Nanchen DONGFANG
! **************************************************************************************************
MODULE qs_electrons_charge_mixing
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE force_env_methods,               ONLY: force_env_create
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_release,&
                                              force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_electrons_types,              ONLY: qs_electrons_type
   USE qs_electrons_workfunction,       ONLY: qs_electrons_compute_workfunction
   USE qs_environment,                  ONLY: qs_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_create,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type,&
                                              set_mo_set
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_set,&
                                              qs_subsys_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_electrons_charge_mixing'

   PUBLIC :: qs_electrons_charge_mixing_lambda, &
             push_lambda_to_mixed_env, &
             compute_mixed_workfunction, &
             charge_mixing_maybe_rebuild, &
             charge_mixing_sync_endpoint_state

CONTAINS

! **************************************************************************************************
!> \brief Given the current electronic coordinate q and the integer span S,
!>        compute the bracketing integers N_lo/N_hi and the mixing weight
!>        lambda = (q - N_lo)/S.
!> \param electrons ...
!> \param lambda ...
!> \param n_lo ...
!> \param n_hi ...
! **************************************************************************************************
   SUBROUTINE qs_electrons_charge_mixing_lambda(electrons, lambda, n_lo, n_hi)
      TYPE(qs_electrons_type), POINTER                   :: electrons
      REAL(KIND=dp), INTENT(OUT)                         :: lambda
      INTEGER, INTENT(OUT), OPTIONAL                     :: n_lo, n_hi

      INTEGER                                            :: my_n_hi, my_n_lo, span

      span = MAX(1, electrons%charge_mixing_span)
      my_n_lo = span*FLOOR(electrons%q/REAL(span, KIND=dp))
      my_n_hi = my_n_lo + span
      IF (electrons%lambda_is_fixed) THEN
         ! LAMBDA given in the input: hold the mixing weight there. n_lo/n_hi are
         ! still reported from q for the caller's benefit, but nothing rebuilds --
         ! charge_mixing_maybe_rebuild returns early in this mode -- so the two
         ! endpoints stay at whatever CHARGE the input built them with.
         lambda = electrons%lambda_fixed
      ELSE
         lambda = (electrons%q - REAL(my_n_lo, KIND=dp))/REAL(span, KIND=dp)
      END IF

      IF (PRESENT(n_lo)) n_lo = my_n_lo
      IF (PRESENT(n_hi)) n_hi = my_n_hi

   END SUBROUTINE qs_electrons_charge_mixing_lambda

! **************************************************************************************************
!> \brief Push a dynamically computed lambda into the MIXED force_eval's
!>        input structure. force_env_methods.F::mixed_energy_forces reads
!>        MIXED%LINEAR%LAMBDA fresh from the input structure every call, so
!>        this alone is sufficient to make energy/force/virial mixing track
!>        our own extended-variable dynamics with zero changes to stock CP2K.
!>
!>        IMPORTANT sign convention: CP2K's own mixed_energy_forces computes
!>        E = LAMBDA*E(force_eval#1) + (1-LAMBDA)*E(force_eval#2). In our
!>        input, force_eval#1 is the N_lo endpoint and force_eval#2 is the
!>        N_hi endpoint, but our own `lambda` argument here uses the opposite
!>        (standard interpolation) convention, lambda = (q-N_lo)/S, i.e.
!>        lambda=0 at N_lo and lambda=1 at N_hi. So CP2K's LAMBDA must be set
!>        to (1-lambda), NOT lambda directly, or the mixing comes out
!>        backwards (mostly N_hi when q is near N_lo).
!> \param force_env the top-level MIXED force_env (force_env%in_use must be
!>        use_mixed_force; force_env%force_env_section must contain a MIXED
!>        subsection with MIXING_TYPE LINEAR_COMBINATION)
!> \param lambda our convention: 0 at N_lo, 1 at N_hi
! **************************************************************************************************
   SUBROUTINE push_lambda_to_mixed_env(force_env, lambda)
      TYPE(force_env_type), POINTER                      :: force_env
      REAL(KIND=dp), INTENT(IN)                          :: lambda

      TYPE(section_vals_type), POINTER                   :: force_env_section, mixed_section

      NULLIFY (force_env_section, mixed_section)
      CALL force_env_get(force_env, force_env_section=force_env_section)
      mixed_section => section_vals_get_subs_vals(force_env_section, "MIXED")
      CALL section_vals_val_set(mixed_section, "LINEAR%LAMBDA", r_val=1.0_dp - lambda)

   END SUBROUTINE push_lambda_to_mixed_env

! **************************************************************************************************
!> \brief Fermi level / workfunction mixing across the two integer-charge
!>        endpoints. Stock MIXED only combines energy/force/virial/dipole
!>        (mixed_environment_utils.F::mixed_map_forces), so these two
!>        quantities -- needed every step by our own O-B-A-B-O potentiostat
!>        to decide how to move q(t) -- are interpolated here, mirroring
!>        FFChargeMixTwoSockets._mix_results's fermi/workfunction lines.
!> \param force_env the top-level MIXED force_env; force_env%sub_force_env(1)
!>        and (2) must each be an ordinary QS force_env
!> \param axis ...
!> \param wf_min ...
!> \param wf_max ...
!> \param q current electronic coordinate (total electron count)
!> \param span charge_mixing_span (S)
!> \param fermi_level_eV ...
!> \param workfunction_eV ...
!> \param lambda_in optional explicit mixing weight; callers pass the value from
!>        qs_electrons_charge_mixing_lambda so a FIXED lambda is honoured here too
! **************************************************************************************************
   SUBROUTINE compute_mixed_workfunction(force_env, axis, wf_min, wf_max, q, span, &
                                         fermi_level_eV, workfunction_eV, lambda_in)
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER, INTENT(IN)                                :: axis
      REAL(KIND=dp), INTENT(IN)                          :: wf_min, wf_max, q
      INTEGER, INTENT(IN)                                :: span
      REAL(KIND=dp), INTENT(OUT)                         :: fermi_level_eV, workfunction_eV
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: lambda_in

      INTEGER                                            :: n_lo
      REAL(KIND=dp)                                      :: fermi_hi, fermi_lo, lambda, wf_hi, wf_lo

! Optional explicit mixing weight. Callers pass the value returned by
! qs_electrons_charge_mixing_lambda so that a FIXED lambda is honoured here
! too; without it this routine would re-derive lambda from q and the mixed
! work function would disagree with the mixed energy.

      CPASSERT(ASSOCIATED(force_env%sub_force_env))
      CPASSERT(SIZE(force_env%sub_force_env) == 2)

      IF (PRESENT(lambda_in)) THEN
         lambda = lambda_in
      ELSE
         n_lo = span*FLOOR(q/REAL(span, KIND=dp))
         lambda = (q - REAL(n_lo, KIND=dp))/REAL(span, KIND=dp)
      END IF

      ! With MIXED%NGROUPS 2 (our default, for concurrent execution), CP2K
      ! splits ranks across the two sub-force_evals, so a given rank only has
      ! ONE of the two sub_force_env's built locally (f77_interface.F::
      ! create_force_env CYCLEs past the other for that rank's group).
      !
      ! Every rank within a given group calls qs_electrons_compute_workfunction
      ! (needed in case any internal fallback tier ever becomes a collective
      ! operation across that group, e.g. an OT calculate_subspace_eigenvalues
      ! path), but the Fermi level / workfunction result is IDENTICAL/replicated
      ! across every rank in that group -- so before the cross-group reduction,
      ! only that group's own "source" rank is allowed to keep its value; every
      ! other rank in the SAME group zeroes it out first. This exactly mirrors
      ! how force_env_methods.F::mixed_energy_forces combines energies(:)
      ! across groups (only para_env%is_source() assigns energies(iforce_eval)
      ! before the shared CALL force_env%para_env%sum(energies)). Without this
      ! is_source() gate, every rank in a group would contribute the SAME
      ! value to the sum, inflating the result by that group's rank count --
      ! which is what produced the ~12x-too-large workfunction seen earlier.
      fermi_lo = 0.0_dp
      wf_lo = 0.0_dp
      fermi_hi = 0.0_dp
      wf_hi = 0.0_dp

      IF (ASSOCIATED(force_env%sub_force_env(1)%force_env)) THEN
         CALL qs_electrons_compute_workfunction(force_env%sub_force_env(1)%force_env, axis, wf_min, wf_max, &
                                                fermi_lo, wf_lo)
         IF (.NOT. force_env%sub_force_env(1)%force_env%para_env%is_source()) THEN
            fermi_lo = 0.0_dp
            wf_lo = 0.0_dp
         END IF
      END IF
      IF (ASSOCIATED(force_env%sub_force_env(2)%force_env)) THEN
         CALL qs_electrons_compute_workfunction(force_env%sub_force_env(2)%force_env, axis, wf_min, wf_max, &
                                                fermi_hi, wf_hi)
         IF (.NOT. force_env%sub_force_env(2)%force_env%para_env%is_source()) THEN
            fermi_hi = 0.0_dp
            wf_hi = 0.0_dp
         END IF
      END IF

      CALL force_env%para_env%sum(fermi_lo)
      CALL force_env%para_env%sum(wf_lo)
      CALL force_env%para_env%sum(fermi_hi)
      CALL force_env%para_env%sum(wf_hi)

      fermi_level_eV = (1.0_dp - lambda)*fermi_lo + lambda*fermi_hi
      workfunction_eV = (1.0_dp - lambda)*wf_lo + lambda*wf_hi

   END SUBROUTINE compute_mixed_workfunction

! **************************************************************************************************
!> \brief Check whether q(t) has drifted outside the currently-built
!>        [N_lo,N_hi] bracket and, if so, rebuild whichever endpoint(s) are
!>        now stale at the new integer electron count. Must be called every
!>        step, right after qs_electrons_drift and before
!>        push_lambda_to_mixed_env, so LAMBDA is always pushed against a
!>        bracket that actually matches the physically-built sub_force_env
!>        pair.
!> \param force_env the top-level MIXED force_env
!> \param electrons ...
! **************************************************************************************************
   SUBROUTINE charge_mixing_maybe_rebuild(force_env, electrons)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(qs_electrons_type), POINTER                   :: electrons

      INTEGER                                            :: n_hi, n_lo
      LOGICAL                                            :: moved
      REAL(KIND=dp)                                      :: lambda_unused

      ! A frozen mixing weight means the endpoints must stay exactly as the input
      ! built them; rebuilding to bracket q would silently change what LAMBDA is
      ! a weight BETWEEN.
      IF (electrons%lambda_is_fixed) RETURN

      ! n_lo = S*FLOOR(q/S), n_hi = n_lo + S, slot 1 <- n_lo, slot 2 <- n_hi.
      ! The endpoints are only updated when the bracket actually changes, by
      ! comparing against the tracked pair, because updating an endpoint mutates
      ! the SCF state and is not free.
      CALL qs_electrons_charge_mixing_lambda(electrons, lambda_unused, n_lo, n_hi)
      moved = (n_lo /= electrons%charge_mixing_current_n_lo) .OR. &
              (n_hi /= electrons%charge_mixing_current_n_hi)

      IF (n_lo /= electrons%charge_mixing_current_n_lo) THEN
         IF (electrons%charge_mixing_inplace) THEN
            CALL charge_mixing_set_endpoint_charge(force_env, 1, n_lo, electrons)
         ELSE
            CALL rebuild_charge_mixing_endpoint(force_env, 1, n_lo, electrons)
         END IF
         electrons%charge_mixing_current_n_lo = n_lo
      END IF
      IF (n_hi /= electrons%charge_mixing_current_n_hi) THEN
         IF (electrons%charge_mixing_inplace) THEN
            CALL charge_mixing_set_endpoint_charge(force_env, 2, n_hi, electrons)
         ELSE
            CALL rebuild_charge_mixing_endpoint(force_env, 2, n_hi, electrons)
         END IF
         electrons%charge_mixing_current_n_hi = n_hi
      END IF

      ! Read the counts back and abort if the endpoints did not actually take
      ! the change, rather than silently mixing the wrong states.
      CALL charge_mixing_verify_endpoints(force_env, electrons, moved)

   END SUBROUTINE charge_mixing_maybe_rebuild

! **************************************************************************************************
!> \brief Record which integer electron counts the two FORCE_EVAL blocks are
!>        ACTUALLY built for, by interrogating them, so that the very first
!>        charge_mixing_maybe_rebuild can tell whether the input agrees with q.
!>
!>        This closes the "startup endpoint mismatch" hole.  The port assumes
!>        throughout that sub_force_env(1) carries N_lo and sub_force_env(2)
!>        carries N_hi -- push_lambda_to_mixed_env relies on it, since CP2K
!>        evaluates E = LAMBDA*E(1) + (1-LAMBDA)*E(2) and we therefore push
!>        LAMBDA = 1 - lambda.  Nothing in stock CP2K enforces that ordering:
!>        an input whose two blocks are written high-charge-first would mix
!>        BACKWARDS, silently, for the whole run, and the only symptom would be
!>        a potentiostat that pushes q the wrong way.
!>
!>        Reading the counts back here and letting maybe_rebuild correct them
!>        removes the dependence on the order of the input blocks.
!>
!>        Reading rather than blindly forcing a rebuild matters: seeding the
!>        tracking with -HUGE(0) would also fix the ordering, but it would force
!>        two full qs_init rebuilds at step 0 of every run, including the vast
!>        majority whose input is already correct.
!> \param force_env the top-level MIXED force_env
!> \param electrons ...
! **************************************************************************************************
   SUBROUTINE charge_mixing_sync_endpoint_state(force_env, electrons)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(qs_electrons_type), POINTER                   :: electrons

      INTEGER, DIMENSION(2)                              :: n_have

      IF (.NOT. electrons%charge_mixing_enabled) RETURN

      CALL charge_mixing_read_endpoint_counts(force_env, n_have)
      electrons%charge_mixing_current_n_lo = n_have(1)
      electrons%charge_mixing_current_n_hi = n_have(2)

   END SUBROUTINE charge_mixing_sync_endpoint_state

! **************************************************************************************************
!> \brief Read back the integer electron count each of the two FORCE_EVAL
!>        endpoints is currently built for.
!>
!>        Under NGROUPS > 1 a rank sees only its OWN group's sub_force_env and
!>        the other pointer is null, so each slot starts at a sentinel that
!>        loses the MAX reduction and is contributed by whichever ranks hold it.
!>        Electron counts are positive, so MAX is safe, and every rank leaves
!>        with the same pair -- which matters because the caller may CPABORT on
!>        the result and that decision has to be collective.
!> \param force_env the top-level MIXED force_env
!> \param n_have (1) = slot 1's count, (2) = slot 2's count
! **************************************************************************************************
   SUBROUTINE charge_mixing_read_endpoint_counts(force_env, n_have)
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER, DIMENSION(2), INTENT(OUT)                 :: n_have

      INTEGER                                            :: idx
      TYPE(force_env_type), POINTER                      :: sub
      TYPE(qs_subsys_type), POINTER                      :: subsys

      n_have = -HUGE(0)
      IF (.NOT. ASSOCIATED(force_env%sub_force_env)) RETURN

      DO idx = 1, 2
         NULLIFY (sub, subsys)
         IF (idx > SIZE(force_env%sub_force_env)) CYCLE
         sub => force_env%sub_force_env(idx)%force_env
         IF (.NOT. ASSOCIATED(sub)) CYCLE
         IF (.NOT. ASSOCIATED(sub%qs_env)) CYCLE
         CALL get_qs_env(sub%qs_env, subsys=subsys)
         IF (.NOT. ASSOCIATED(subsys)) CYCLE
         CALL qs_subsys_get(subsys, nelectron_total=n_have(idx))
      END DO

      CALL force_env%para_env%max(n_have)

   END SUBROUTINE charge_mixing_read_endpoint_counts

! **************************************************************************************************
!> \brief Assert that the two endpoints really are built for the electron counts
!>        we think they are, and log the bracket whenever it moves.
!>
!>        The in-place path mutates dft_control%charge and
!>        the subsys electron counts directly, so a silent failure to take
!>        effect would look exactly like a converged run -- just one mixing two
!>        copies of the same integer-charge state, which makes E independent of
!>        lambda and leaves q with no restoring force at all.  Checking is cheap
!>        and turns that failure mode into an abort instead of a plausible
!>        trajectory.
!> \param force_env ...
!> \param electrons ...
!> \param moved whether the bracket changed on this step (controls logging)
! **************************************************************************************************
   SUBROUTINE charge_mixing_verify_endpoints(force_env, electrons, moved)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(qs_electrons_type), POINTER                   :: electrons
      LOGICAL, INTENT(IN)                                :: moved

      CHARACTER(LEN=256)                                 :: msg
      INTEGER                                            :: iw
      INTEGER, DIMENSION(2)                              :: n_have

      CALL charge_mixing_read_endpoint_counts(force_env, n_have)

      IF (n_have(1) /= electrons%charge_mixing_current_n_lo .OR. &
          n_have(2) /= electrons%charge_mixing_current_n_hi) THEN
         WRITE (msg, '(A,I0,A,I0,A,I0,A,I0,A)') &
            "CHARGE MIXING: endpoints disagree with the requested bracket. Wanted (", &
            electrons%charge_mixing_current_n_lo, ",", electrons%charge_mixing_current_n_hi, &
            ") but the FORCE_EVALs hold (", n_have(1), ",", n_have(2), ")."
         CPABORT(TRIM(msg))
      END IF

      ! A collapsed bracket makes E independent of lambda: the mixed energy is
      ! then a single integer-charge result, dPhi/dq is zero, and q drifts under
      ! no restoring force whatsoever.  It cannot arise from the arithmetic
      ! above, but it is the one corruption that would masquerade as a healthy
      ! trajectory, so it is worth naming explicitly.
      IF (n_have(2) - n_have(1) /= MAX(1, electrons%charge_mixing_span)) THEN
         WRITE (msg, '(A,I0,A,I0,A,I0,A)') &
            "CHARGE MIXING: endpoint separation is ", n_have(2) - n_have(1), &
            " but CHARGE_MIXING_SPAN is ", MAX(1, electrons%charge_mixing_span), &
            ". The two FORCE_EVALs no longer bracket q."
         CPABORT(TRIM(msg))
      END IF

      IF (moved) THEN
         iw = cp_logger_get_default_io_unit()
         IF (iw > 0) THEN
            WRITE (iw, '(T2,A,I0,A,I0,A,F14.8)') &
               "CHARGE MIXING| bracket now [", n_have(1), ",", n_have(2), "]  q = ", electrons%q
         END IF
      END IF

   END SUBROUTINE charge_mixing_verify_endpoints

! **************************************************************************************************
!> \brief Tear down and rebuild force_env%sub_force_env(idx) at a new integer
!>        electron count, preserving the CURRENT (dynamics-evolved) geometry
!>        and cell by handing the old subsys straight to qs_init. This is the
!>        fallback path used when CHARGE_MIXING_INPLACE is off.
!> \param force_env the top-level MIXED force_env
!> \param idx 1 (the N_lo endpoint) or 2 (the N_hi endpoint)
!> \param new_n_electron the new integer electron count for this endpoint
!> \param electrons used only for neutral_electrons, to convert electron
!>        count back into CP2K's own CHARGE = neutral - n_electron convention
! **************************************************************************************************
   SUBROUTINE rebuild_charge_mixing_endpoint(force_env, idx, new_n_electron, electrons)
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER, INTENT(IN)                                :: idx, new_n_electron
      TYPE(qs_electrons_type), POINTER                   :: electrons

      INTEGER                                            :: my_group, new_charge
      LOGICAL                                            :: subsys_explicit
      TYPE(cp_logger_type), POINTER                      :: sub_logger
      TYPE(cp_subsys_type), POINTER                      :: old_subsys
      TYPE(force_env_type), POINTER                      :: new_sub, old_sub
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_environment_type), POINTER                 :: new_qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section, force_env_section, &
                                                            mixed_subsys_section, root_section, &
                                                            subsys_section

      NULLIFY (new_sub, old_sub, dft_section, force_env_section, subsys_section, &
               para_env, globenv, root_section, new_qs_env, old_subsys, sub_logger, &
               mixed_subsys_section)

      old_sub => force_env%sub_force_env(idx)%force_env
      IF (.NOT. ASSOCIATED(old_sub)) RETURN   ! not this rank's group under NGROUPS>1

      ! Install this group's own logger for the duration of the rebuild, exactly as
      ! mixed_energy_forces does before it touches a sub_force_env.  This is NOT
      ! cosmetic: distribute_molecules_1d (distribution_methods.F:102) takes NO
      ! para_env argument at all.  It reads its communicator straight off the
      ! default logger --
      !
      !    logger => cp_get_default_logger()
      !    ASSOCIATE (group => logger%para_env, mype => logger%para_env%mepos + 1, &
      !               npe   => logger%para_env%num_pe)
      !       ...
      !       nbins = npe
      !
      ! -- deals the atoms into nbins bins and lets a rank claim one only when
      ! ipe == mype.  At startup mixed_energy_forces has already installed the
      ! group logger, so the atoms are dealt across the GROUP.  A rebuild driven
      ! from the MD integrator runs outside mixed_energy_forces, where the default
      ! logger is the GLOBAL one: the atoms are then dealt across every rank in the
      ! job while the DFT-D3 energy is later summed over this group only, so every
      ! atom binned to a rank outside the group is silently dropped.  Measured:
      ! TOT_vdw 14112 -> 7042 pairs, rank-0 local particles 2 -> 1, at NGROUPS 2
      ! with a 16-rank group inside a 32-rank job.
      !
      ! distribution_2d is immune because distribute_molecules_2d is handed its
      ! blacs_env explicitly (rows2d/cols2d were verified unchanged at 9/9), and
      ! sab_vdw/sab_cn are the only neighbour lists affected because
      ! operator_type="PP" -> otype 4 is the only branch of build_neighbor_lists
      ! that reads atom2d%list_1d, which is built from local_particles.
      !
      ! The same wrong-logger para_env is what produced MPI_ERR_TRUNCATE in
      ! mp_sum_l at the first crossing: qs_init logs collectively, and outside
      ! mixed_energy_forces the DEFAULT logger's para_env spans every rank.
      IF (ASSOCIATED(force_env%mixed_env)) THEN
         IF (ASSOCIATED(force_env%mixed_env%sub_logger) .AND. &
             ASSOCIATED(force_env%mixed_env%group_distribution)) THEN
            my_group = force_env%mixed_env%group_distribution(force_env%para_env%mepos)
            sub_logger => force_env%mixed_env%sub_logger(my_group + 1)%p
         END IF
      END IF
      IF (ASSOCIATED(sub_logger)) CALL cp_add_default_logger(sub_logger)

      new_charge = NINT(electrons%neutral_electrons) - new_n_electron

      para_env => old_sub%para_env
      globenv => old_sub%globenv
      root_section => old_sub%root_section
      force_env_section => old_sub%force_env_section
      CALL force_env_get(old_sub, subsys=old_subsys)

      ! Mutate CHARGE in place on the (shared, reference-counted) input
      ! section before rebuilding -- qs_init reads DFT%CHARGE from
      ! force_env_section exactly as it did on the original, once-at-startup
      ! build, so this alone determines the new endpoint's electron count.
      dft_section => section_vals_get_subs_vals(force_env_section, "DFT")
      CALL section_vals_val_set(dft_section, "CHARGE", i_val=new_charge)

      ! Which &SUBSYS does qs_init read the kinds and basis sets from?
      !
      ! In the ORIGINAL input each sub-FORCE_EVAL carries its own &SUBSYS, so
      ! taking the child of this endpoint's force_env_section works. In a RESTART
      ! file it does not: CP2K's restart writer emits &SUBSYS exactly once, on the
      ! force_eval that owns the particles, and leaves the sub-FORCE_EVALs without
      ! one. CP2K's own startup path copes with that, which is why a restarted run
      ! reaches the first SCF happily -- but this rebuild then reads an empty
      ! section and aborts, first in cell_methods.F ("cell volume is too small")
      ! and, once the cell is inherited from the live environment, in
      ! qs_kind_types.F ("No &KIND section ... to associate to the atomic kind").
      !
      ! So: use this endpoint's own &SUBSYS when it is explicit, and otherwise
      ! fall back to the MIXED force_env's. The two describe the same atoms; only
      ! DFT%CHARGE differs between endpoints, and that lives in the &DFT section
      ! mutated above, not here.
      subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      CALL section_vals_get(subsys_section, explicit=subsys_explicit)
      IF (.NOT. subsys_explicit) THEN
         IF (ASSOCIATED(force_env%force_env_section)) THEN
            mixed_subsys_section => section_vals_get_subs_vals(force_env%force_env_section, "SUBSYS")
            CALL section_vals_get(mixed_subsys_section, explicit=subsys_explicit)
            IF (subsys_explicit) subsys_section => mixed_subsys_section
         END IF
      END IF

      ALLOCATE (new_qs_env)
      CALL qs_env_create(new_qs_env, globenv)
      ! cp_subsys=old_subsys hands the current MD-evolved positions, velocities and
      ! cell (including the reference cell) straight to the new environment, so
      ! neither positions nor the cell are re-read from the input file.
      CALL qs_init(new_qs_env, para_env, root_section=root_section, globenv=globenv, &
                   cp_subsys=old_subsys, force_env_section=force_env_section, &
                   subsys_section=subsys_section, use_motion_section=.FALSE.)

      ! force_env_create's signature IS confirmed (force_env_methods.F line
      ! 728): it ALLOCATEs new_sub itself (must be passed unallocated, which
      ! it is here), sets ref_count=1, and -- critically -- RETAINS every
      ! shared object it is handed before storing it:
      !   force_env%globenv => globenv;        CALL globenv_retain(...)
      !   force_env%root_section => root_section; CALL section_vals_retain(...)
      !   force_env%para_env => para_env;      CALL force_env%para_env%retain()
      !   CALL section_vals_retain(force_env_section); force_env%force_env_section => ...
      ! So handing it the SAME para_env/globenv/root_section/force_env_section
      ! that old_sub already owns is safe: new_sub's copies are retained
      ! (ref_count bumped) here, BEFORE old_sub is released below, so
      ! old_sub's release cannot pull any of them to zero out from under
      ! new_sub. qs_env is NOT retained (force_env%qs_env => qs_env is a
      ! direct, unshared assignment) -- correct here since new_qs_env was
      ! freshly allocated solely for this new_sub, not shared with old_sub.
      CALL force_env_create(new_sub, root_section=root_section, para_env=para_env, &
                            globenv=globenv, qs_env=new_qs_env, force_env_section=force_env_section)

      ! force_env_release's signature IS confirmed (force_env_types.F line
! **************************************************************************************************
!> \brief ...
!> \param force_env ...
! **************************************************************************************************
      ! 206): RECURSIVE SUBROUTINE force_env_release(force_env) with
      ! force_env as TYPE(force_env_type), POINTER -- exactly as called here.
      ! It is refcounted (decrements ref_count, only actually deallocates at
      ! zero) and nullifies the pointer on return, so old_sub itself does not
      ! need a manual NULLIFY afterward.
      ! ---------------------------------------------------------------------
      ! REMOVED: the sab_vdw / sab_cn handover from the old dispersion_env to
      ! the new one.  It was tested on-cluster and changed nothing (E_disp
      ! stayed at -0.16998 Ha), and the reason is structural, not empirical:
      ! build_neighbor_lists opens with CALL release_neighbor_list_sets(ab_list)
      ! (qs_neighbor_lists.F:1148), and build_qs_neighbor_lists calls it
      ! unconditionally for both sab_vdw and sab_cn.  Whatever is handed over
      ! here is therefore freed and regenerated at the new endpoint's first SCF.
      ! The handover was dead code.  Any fix must act at the construction site,
      ! not here.  The dispersion_env parameter dump that used to sit here is
      ! also gone: it answered its question (every parameter is identical
      ! between the two builds).  The live diagnostic is now the @@ NLDIAG /
      ! @@ NLDIST / @@ NLTOTAL block in qs_neighbor_lists.F, which reports pair
      ! counts and the atomic distribution for both the startup build and the
      ! rebuild in the same job.
      ! ---------------------------------------------------------------------

      CALL force_env_release(old_sub)
      force_env%sub_force_env(idx)%force_env => new_sub

      IF (ASSOCIATED(sub_logger)) CALL cp_rm_default_logger()

   END SUBROUTINE rebuild_charge_mixing_endpoint

! **************************************************************************************************
!> \brief Change the electron count of an existing endpoint IN PLACE, without
!>        re-entering qs_init.
!>
!>        Both endpoints persist across boundary crossings and only their
!>        electron count changes. Everything qs_init constructs that does NOT depend on the
!>        electron count -- dispersion_env and its neighbour lists, pw_env, the
!>        atomic distributions, the basis, the wavefunction -- simply persists.
!>
!>        qs_environment.F derives the whole electron-count-dependent state from
!>        one line (line 1549):
!>
!>            nelectron = nelectron - dft_control%charge
!>
!>        and from it sets nelectron_spin / maxocc / n_mo, storing the counts on
!>        the subsys and the mo_sets. Reproducing just that chain is sufficient.
!>
!>        The allocated MO count (mo_set%nmo) is deliberately NOT touched: nmo
!>        is n_occupied + ADDED_MOS as qs_init sized it for the INPUT charge,
!>        and the occupations are recomputed from the Fermi level at every SCF,
!>        so nothing stale survives. Each crossing moves an endpoint by
!>        CHARGE_MIXING_SPAN electrons, so ADDED_MOS has to cover the furthest
!>        endpoint the run reaches; that is checked below rather than assumed.
!>
!>        OT is refused. Under OT every one of the nmo columns is fully occupied
!>        (no ADDED_MOS, no fractional occupation), so changing nelectron without
!>        resizing mo_coeff leaves the density at the OLD electron count while
!>        the subsys reports the new one. Use CHARGE_MIXING_INPLACE F with OT.
!>
!>        Restricted (nspins=1): the whole count goes in one channel; with
!>        smearing it need not be even. Unrestricted: the multiplicity is
!>        re-derived exactly as qs_init does it (qs_environment.F ~1551). The
!>        value dft_control%multiplicity holds is NOT reusable: qs_init has
!>        already overwritten the input's MULTIPLICITY 0 ("automatic") with 1 or
!>        2 for the INPUT parity, so reusing it after an odd span crossing gave
!>        (N+m-1)/2 + (N-m+1)/2 = N-1 electrons -- silently, since
!>        nelectron_total was still set to N. MULTIPLICITY is therefore re-read
!>        from the input: 0 is resolved again for the new parity, an explicit
!>        value is kept and must be compatible, as qs_init requires.
!> \param force_env the top-level MIXED force_env
!> \param idx which sub_force_env (1 = N_lo slot, 2 = N_hi slot)
!> \param new_n_electron the endpoint's new integer electron count
!> \param electrons ...
! **************************************************************************************************
   SUBROUTINE charge_mixing_set_endpoint_charge(force_env, idx, new_n_electron, electrons)
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER, INTENT(IN)                                :: idx, new_n_electron
      TYPE(qs_electrons_type), POINTER                   :: electrons

      CHARACTER(LEN=256)                                 :: msg
      INTEGER                                            :: input_multiplicity, ispin, multiplicity, &
                                                            new_charge, nmo
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      REAL(KIND=dp)                                      :: maxocc
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_type), POINTER                      :: sub
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, sub_section

      NULLIFY (sub, dft_control, mos, subsys, scf_control, dft_section, sub_section)

      sub => force_env%sub_force_env(idx)%force_env
      IF (.NOT. ASSOCIATED(sub)) RETURN     ! not this rank's group under NGROUPS>1
      IF (.NOT. ASSOCIATED(sub%qs_env)) RETURN

      new_charge = NINT(electrons%neutral_electrons) - new_n_electron

      CALL get_qs_env(sub%qs_env, dft_control=dft_control, mos=mos, subsys=subsys, &
                      scf_control=scf_control)

      IF (ASSOCIATED(scf_control)) THEN
         IF (scf_control%use_ot) THEN
            msg = "CHARGE_MIXING_INPLACE cannot change the electron count under OT; "// &
                  "use CHARGE_MIXING_INPLACE F"
            CPABORT(TRIM(msg))
         END IF
      END IF

      IF (dft_control%nspins == 1) THEN
         multiplicity = dft_control%multiplicity
         nelectron_spin(1) = new_n_electron
         nelectron_spin(2) = 0
      ELSE
         CALL force_env_get(sub, force_env_section=sub_section)
         dft_section => section_vals_get_subs_vals(sub_section, "DFT")
         CALL section_vals_val_get(dft_section, "MULTIPLICITY", i_val=input_multiplicity)
         multiplicity = input_multiplicity
         IF (multiplicity == 0) THEN
            IF (MODULO(new_n_electron, 2) == 0) THEN
               multiplicity = 1
            ELSE
               multiplicity = 2
            END IF
         END IF
         IF (MODULO(new_n_electron + multiplicity - 1, 2) /= 0) THEN
            WRITE (msg, '(A,I0,A,I0,A)') "CHARGE MIXING: MULTIPLICITY ", multiplicity, &
               " is incompatible with ", new_n_electron, &
               " electrons. Use MULTIPLICITY 0 (automatic) or an even CHARGE_MIXING_SPAN."
            CPABORT(TRIM(msg))
         END IF
         nelectron_spin(1) = (new_n_electron + multiplicity - 1)/2
         nelectron_spin(2) = (new_n_electron - multiplicity + 1)/2
         IF (nelectron_spin(2) < 0) then
           CPABORT("CHARGE MIXING: too few electrons for this MULTIPLICITY")
         end if
      END IF

      ! Every orbital the new occupation needs must already be allocated: nmo is
      ! never resized here. Without this an endpoint pushed past n_occ+ADDED_MOS
      ! would quietly hold fewer electrons than nelectron_total claims.
      IF (ASSOCIATED(mos)) THEN
         IF (SIZE(mos) /= dft_control%nspins) THEN
            msg = "CHARGE MIXING: mo_set count differs from nspins"
            CPABORT(TRIM(msg))
         END IF
         DO ispin = 1, SIZE(mos)
            CALL get_mo_set(mos(ispin), nmo=nmo, maxocc=maxocc)
            IF (REAL(nelectron_spin(ispin), KIND=dp) > maxocc*REAL(nmo, KIND=dp) + 1.0E-8_dp) THEN
               WRITE (msg, '(A,I0,A,I0,A,I0,A)') "CHARGE MIXING: endpoint at ", new_n_electron, &
                  " electrons needs more orbitals in spin ", ispin, " than the ", nmo, &
                  " allocated. Increase ADDED_MOS or use CHARGE_MIXING_INPLACE F."
               CPABORT(TRIM(msg))
            END IF
         END DO
      END IF

      dft_control%charge = REAL(new_charge, KIND=dp)
      dft_control%multiplicity = multiplicity

      CALL qs_subsys_set(subsys, nelectron_total=new_n_electron, &
                         nelectron_spin=nelectron_spin)

      IF (ASSOCIATED(mos)) THEN
         DO ispin = 1, SIZE(mos)
            CALL set_mo_set(mos(ispin), nelectron=nelectron_spin(ispin), &
                            n_el_f=REAL(nelectron_spin(ispin), KIND=dp))
         END DO
      END IF

   END SUBROUTINE charge_mixing_set_endpoint_charge

END MODULE qs_electrons_charge_mixing
