!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Methods for sampling helium variables
!> \author Lukasz Walewski
!> \date   2009-06-10
! *****************************************************************************
MODULE helium_sampling

  USE cp_external_control,             ONLY: external_control
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_rm_iter_level
  USE global_types,                    ONLY: global_environment_type
  USE helium_common,                   ONLY: &
       helium_boxmean_3d, helium_calc_plength, helium_calc_rdf, &
       helium_calc_rho, helium_calc_wnumber, helium_eval_expansion, &
       helium_pbc, helium_rotate, helium_spline, &
       helium_update_transition_matrix
  USE helium_interactions,             ONLY: helium_bead_solute_e,&
                                             helium_calc_energy,&
                                             helium_solute_e_f
  USE helium_io,                       ONLY: &
       helium_write_accepts, helium_write_coordinates, helium_write_energy, &
       helium_write_force, helium_write_line, helium_write_perm, &
       helium_write_plength, helium_write_rdf, helium_write_rho, &
       helium_write_sdensity, helium_write_wnumber
  USE helium_types,                    ONLY: e_id_total,&
                                             helium_solvent_type
  USE input_cp2k_restarts,             ONLY: write_restart
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_walltime
  USE message_passing,                 ONLY: mp_sum
  USE parallel_rng_types,              ONLY: next_random_number
  USE pint_types,                      ONLY: pint_env_type
  USE splines_types,                   ONLY: spline_data_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_sampling'

  PUBLIC :: helium_do_run
  PUBLIC :: helium_sample
  PUBLIC :: helium_step

  CONTAINS

  ! ***************************************************************************
  !> \brief  Performs MC simulation for helium (only)
  !> \author Lukasz Walewski
  !> \date   2009-07-14
  !> \note   This routine gets called only when HELIUM_ONLY is set to .TRUE.,
  !>         so do not put any property calculations here!
  ! ***************************************************************************
  SUBROUTINE helium_do_run(helium, globenv, error)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_do_run', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: step, tot_steps
    LOGICAL                                  :: failure, should_stop
    TYPE(pint_env_type), POINTER             :: pint_env

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    IF (failure) THEN
      CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

    NULLIFY(pint_env)

    NULLIFY(helium%logger)
    helium%logger => cp_error_get_logger(error)

    ! create iteration level
    ! The iteration level 'MD' can be misleading since this is a pure MC
    ! simulation, however the helium code is meant to be run mainly from PINT.
    ! We keep it like that for compatibility with PINT+HELIUM simulations.
    CALL cp_add_iter_level(helium%logger%iter_info,"MD",error=error)
    CALL cp_iterate(helium%logger%iter_info,iter_nr=helium%first_step,error=error)
    tot_steps = helium%first_step
    DO step = 1, helium%num_steps

      tot_steps = tot_steps + 1
      CALL cp_iterate(helium%logger%iter_info,last=(step==helium%num_steps),iter_nr=tot_steps,error=error)
      helium%current_step = tot_steps

      CALL helium_step(helium,pint_env,error)

      ! call write_restart here to avoid interference with PINT write_restart
      CALL write_restart(root_section=helium%input,helium_env=helium,error=error)

      ! exit from the main loop if soft exit has been requested
      CALL external_control(should_stop,"MD",globenv=globenv,error=error)
      IF (should_stop) EXIT

    END DO

    ! remove iteration level
    CALL cp_rm_iter_level(helium%logger%iter_info,"MD",error=error)

    RETURN
  END SUBROUTINE helium_do_run

  ! ***************************************************************************
  !> \brief  Sample the helium phase space
  !> \author Lukasz Walewski
  !> \date   2009-10-27
  !> \descr  Samples helium variable space according to multilevel Metropolis
  !>         MC scheme, calculates the forces exerted by helium solvent on the
  !>         solute and stores them in helium%force_avrg array. The forces are
  !>         averaged over outer MC loop.
  !> \note   The implicit assumption is that we simulate solute _with_ solvent
  !>         most of the time, so for performance reasons I do not check that
  !>         everywhere I should. This leads to some redundancy in the case of
  !>         helium-only calculations.
  ! ***************************************************************************
  SUBROUTINE helium_sample( helium, pint_env, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_sample', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: irot, nslices
    LOGICAL                                  :: failure, solute_present
    REAL(kind=dp)                            :: inv_iter_rot, inv_num_pe, &
                                                rnd, rtmp
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    solute_present = helium%solute_present
    helium%force_avrg(:,:) = 0.0_dp
    helium%energy_avrg(:) = 0.0_dp
    helium%energy_inst(:) = 0.0_dp
    helium%wnumber_avrg(:) = 0.0_dp
    helium%sdensity_avrg = 0.0_dp
    helium%plength_avrg(:) = 0.0_dp
    helium%rdf_avrg(:) = 0.0_dp
    helium%num_accepted(:,:) = 0.0_dp

    ! helium parallelization: each processor gets different RN stream and
    ! runs independent helium simulation, the properties and forces are
    ! averaged over parallel helium environments once per step.

    ! helium sampling (outer MC loop)
    inv_iter_rot = 1.0_dp / REAL(helium%iter_rot,dp)
    DO irot = 1, helium%iter_rot

      ! rotate helium beads in imaginary time at random, store current
      ! 'rotation state' in helium%relrot wich is within (0, helium%beads-1)
      ! (this is needed to sample different fragments of the permutation
      ! paths in try_permutations)
      rnd = next_random_number(helium%rng_stream_uniform,error=error)
      nslices = rnd * helium%beads
      CALL helium_rotate(helium,nslices,error)

      CALL helium_try_permutations(helium,pint_env,error)

      ! calculate instantaneous forces on the solute and He properties
      IF (solute_present) CALL helium_solute_e_f(pint_env, helium, rtmp)
      CALL helium_calc_energy(helium,pint_env)
      CALL helium_calc_wnumber( helium )
      CALL helium_calc_plength( helium, error )
      CALL helium_calc_rdf( helium, error )
      IF (helium%rho_present) CALL helium_calc_rho( helium, error )

      ! instantaneous force output according to HELIUM%PRINT%FORCES_INST
      ! Warning: file I/O here may cost A LOT of cpu time!
      ! switched off here to save cpu
      !CALL helium_write_force_inst( helium, error )

      ! collect instantaneous values for averaging
      helium%force_avrg(:,:) = &
        helium%force_avrg(:,:) + helium%force_inst(:,:)
      helium%energy_avrg(:) = helium%energy_avrg(:) + helium%energy_inst(:)
      helium%wnumber_avrg(:) = helium%wnumber_avrg(:) + helium%wnumber_inst(:)
      helium%sdensity_avrg = helium%sdensity_avrg + helium%sdensity_inst
      helium%plength_avrg(:) = helium%plength_avrg(:) + helium%plength_inst(:)
      helium%rdf_avrg(:) = helium%rdf_avrg(:) + helium%rdf_inst(:)
      IF (helium%rho_present) THEN
        helium%rho_avrg(:,:,:,:) = helium%rho_avrg(:,:,:,:) + &
                                   helium%rho_inst(:,:,:,:) * inv_iter_rot
      END IF

    END DO

    ! restore the original alignment of beads in imaginary time
    ! (this is useful to make a single bead's position easy to follow
    ! in the trajectory, otherwise it's index would change every step)
    CALL helium_rotate(helium,-helium%relrot, error)

    ! actually average the forces and properties over the outer MC loop
    helium%force_avrg(:,:) = helium%force_avrg(:,:) * inv_iter_rot
    helium%energy_avrg(:) = helium%energy_avrg(:) * inv_iter_rot
    helium%wnumber_avrg(:) = helium%wnumber_avrg(:) * inv_iter_rot
    helium%sdensity_avrg = helium%sdensity_avrg * inv_iter_rot
    helium%plength_avrg(:) = helium%plength_avrg(:) * inv_iter_rot
    helium%rdf_avrg(:) = helium%rdf_avrg(:) * inv_iter_rot

    ! average over helium environments sitting at different processors
    ! WARNING: do not average winding number, etc. over He environments!
    inv_num_pe = 1.0_dp / REAL(logger%para_env%num_pe,dp)
    IF (solute_present) CALL mp_sum(helium%force_avrg,logger%para_env%group)
    CALL mp_sum(helium%energy_avrg,logger%para_env%group)
    CALL mp_sum(helium%sdensity_avrg,logger%para_env%group)
    CALL mp_sum(helium%plength_avrg,logger%para_env%group)
    CALL mp_sum(helium%rdf_avrg,logger%para_env%group)
    CALL mp_sum(helium%num_accepted,logger%para_env%group)
    helium%force_avrg(:,:) = helium%force_avrg(:,:) * inv_num_pe
    helium%energy_avrg(:) = helium%energy_avrg(:) * inv_num_pe
    helium%sdensity_avrg = helium%sdensity_avrg * inv_num_pe
    helium%plength_avrg(:) = helium%plength_avrg(:) * inv_num_pe
    helium%rdf_avrg(:) = helium%rdf_avrg(:) * inv_num_pe
    helium%num_accepted(:,:) = helium%num_accepted(:,:) * inv_num_pe

    RETURN
  END SUBROUTINE helium_sample

  ! ***************************************************************************
  !> \brief  Perform MC step for helium
  !> \author Lukasz Walewski
  !> \date   2009-06-12
  ! ***************************************************************************
  SUBROUTINE helium_step( helium, pint_env, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_step', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: msgstr, stmp, time_unit
    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: time_start, time_stop, &
                                                time_used

    CALL timeset(routineN,handle)
    time_start = m_walltime()

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    ! perform the actual phase space sampling
    CALL helium_sample( helium, pint_env, error )

    ! write out the averaged properties
    CALL helium_write_energy( helium, error )
    CALL helium_write_sdensity( helium, error )
    CALL helium_write_wnumber( helium, error )
    CALL helium_write_plength( helium, error )
    CALL helium_write_rdf( helium, error )
    IF (helium%rho_present) CALL helium_write_rho( helium, error )
    CALL helium_write_perm( helium, error )
    CALL helium_write_accepts( helium, error )
    CALL helium_write_coordinates( helium, error=error )
    CALL helium_write_force( helium, error )

    time_stop = m_walltime()
    time_used = time_stop - time_start
    time_unit = "sec"
    IF ( time_used .GE. 60.0_dp ) THEN
      time_used = time_used / 60.0_dp
      time_unit = "min"
      IF ( time_used .GE. 60.0_dp ) THEN
        time_used = time_used / 60.0_dp
        time_unit = "hours"
      END IF
    END IF
    msgstr = "MC step"
    stmp = ""
    WRITE(stmp, *) helium%current_step
    msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(stmp)) // " of"
    stmp = ""
    WRITE(stmp, *) helium%last_step
    msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(stmp)) // " in"
    stmp = ""
    WRITE(stmp, '(F20.1)') time_used
    msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(stmp))
    msgstr = TRIM(ADJUSTL(msgstr)) // " " // TRIM(ADJUSTL(time_unit)) // "."
    CALL helium_write_line(TRIM(msgstr), error)

    ! print out the total energy - for regtest evaluation
    stmp = ""
    WRITE(stmp, *) helium%energy_avrg(e_id_total)
    msgstr = "Total energy = " // TRIM(ADJUSTL(stmp))
    CALL helium_write_line(TRIM(msgstr), error)

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_step

  ! ***************************************************************************
  !> \brief  ...
  !> \param  pint_env  - path integral environment
  !> \param  error     - error handling (see module cp_error_handling)
  !> \author hforbert
  !> \par    History
  !>        2010-06-17 ratio for m-value added (m-sampling related) [lwalewski]
  ! ***************************************************************************
  SUBROUTINE helium_try_permutations(helium, pint_env, error)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_try_permutations', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ncycles, ni, nselected, res
    REAL(KIND=dp)                            :: r, rnd, x

   IF (helium%maxcycle>1) CALL helium_update_transition_matrix(helium,error)
     helium%work(:,:,:)=helium%pos(:,:,:)

   ! the inner MC loop (without rotation in imaginary time)
   DO ni = 1, helium%iter_norot

      ! set the probability threshold for m_value: 1/(1+(m-1)/helium%m_ratio)
      r = 1.0_dp / ( 1.0_dp + (helium%maxcycle-1)/helium%m_ratio )

      ! draw permutation length for this trial from the distribution of choice
      x = next_random_number(helium%rng_stream_uniform,error=error)
      IF ( x .LT. r ) THEN
        ncycles = helium%m_value
      ELSE
        DO
          x = next_random_number(helium%rng_stream_uniform,error=error)
          ncycles = INT(helium%maxcycle*x)+1
          IF ( ncycles .NE. helium%m_value ) EXIT
        END DO
      END IF

      IF (ncycles<1) ncycles = 1
      IF (ncycles>helium%maxcycle) ncycles = helium%maxcycle
      helium%num_accepted(1,ncycles) = helium%num_accepted(1,ncycles) + 1

      ! check, if permutation of this length can be constructed
      IF (ncycles == 1) THEN
         rnd = next_random_number(helium%rng_stream_uniform,error=error)
         helium%ptable(1)=1+INT(rnd*helium%atoms)
         helium%ptable(2)=-1
         helium%pweight=0.0_dp
         nselected = 1
      ELSE
         nselected = helium_select_permutation(helium,ncycles,error)
      END IF

      IF (nselected /= 0) THEN
         ! the permutation was successfully selected - actually sample this permutation
         res = helium_slice_metro_cyclic(helium, pint_env, ncycles, error)
      END IF

    END DO

    RETURN
  END SUBROUTINE helium_try_permutations

  ! ***************************************************************************
  !> \brief ...
  !> \return 1 - if the MC move was accepted, 0 - otherwise
  !> \author hforbert
  ! ***************************************************************************
  FUNCTION helium_slice_metro_cyclic(helium,pint_env,n,error) RESULT(res)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env
    INTEGER, INTENT(IN)                      :: n
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: res

    INTEGER                                  :: c, hbeads, i, j, k, l, level, &
                                                m, nb, pk1, pk2, stride
    INTEGER, DIMENSION(:), POINTER           :: p, perm
    LOGICAL                                  :: nperiodic
    REAL(KIND=dp)                            :: cell_size, ds, dtk, e1, e2, &
                                                pds, prev_ds, r, sigma, tmp, x
    REAL(KIND=dp), DIMENSION(3)              :: rm1, rm2, tmp1, tmp2, x1
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pos, work
    TYPE(spline_data_p_type), &
      DIMENSION(:, :), POINTER               :: uij

! trial permutation implicit in p
! since we (momentarily) only use cyclic permutations:
! n = 1 : no permutation, sample p[0] anew
! n = 2 : p[0] -> p[1] -> p[0]
! n = 3 : p[0] -> p[1] -> p[2] -> p[0]
! n = 4 : p[0] -> p[1] -> p[2] -> p[3] -> p[0]

  m = helium%bisection
  p => helium%ptable
  prev_ds = helium%pweight

  helium%num_accepted(2,n) = helium%num_accepted(2,n) + 1
  level = 1
  res = 0

IF (m==0) RETURN
   hbeads = helium%beads
   nb = helium%atoms
   pos => helium%pos
   work => helium%work
   perm => helium%permutation
   uij => helium%uij
   cell_size = (0.5_dp*helium%cell_size)**2
   nperiodic = .NOT.helium%periodic

!TODO: use CP2k error/failure mechanisms

   pds = prev_ds
   i = hbeads - m + 1

   ! sanity checks:
   IF (i<1) THEN
      PRINT *,"slice_metro_cyclic: i<1 test failed: ",i
      STOP
   END IF
   j = 1
   k = m
   DO
      IF (k<2) EXIT
      j = j * 2
      k = k/2
   END DO
   IF (j /= m) THEN
      PRINT *,"slice_metro_cyclic: m not a power of 2! ",m
      STOP
   END IF
   IF (m < 2) THEN
      PRINT *,"slice_metro_cyclic: m less than 2! ",m
      STOP
   END IF

!   work(:,:,:) = pos(:,:,:)

   stride = m
   DO
      IF (stride <= 2) EXIT
      ! calc new trial positions
      ! trial action: 0.5*stride*endpointapprox
      sigma = SQRT(0.5_dp*helium%hb2m*(stride/2)*helium%tau)
      dtk = 0.0_dp
      ds = 0.0_dp

      j = i+stride/2
      DO
         IF (j>hbeads-stride/2) EXIT
         pk1 = j-stride/2
         pk2 = j+stride/2
         ! calculate log(T(s)):
         DO k = 1, n
            CALL helium_boxmean_3d( helium, pos(:,p(k),pk1), pos(:,p(k),pk2), x1 )
            tmp1(:) = x1(:) - pos(:,p(k),j)
            CALL helium_pbc( helium, tmp1 )
            tmp1(:) = tmp1(:) / sigma
            dtk = dtk - 0.5_dp * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
         END DO
         ! calculate log(T(sprime)) and sprime itself
         DO k = 1, n
            CALL helium_boxmean_3d( helium, work(:,p(k),pk1), work(:,p(k),pk2), tmp1 )
            DO c = 1, 3
               x = next_random_number(rng_stream=helium%rng_stream_gaussian,&
                   variance=1.0_dp,error=error)
               x = sigma * x
               tmp1(c) = tmp1(c) + x
               tmp2(c) = x
            END DO
            CALL helium_pbc( helium, tmp1 )
            CALL helium_pbc( helium, tmp2 )
            work(:,p(k),j) = tmp1(:)
            tmp2(:) = tmp2(:) / sigma
            dtk = dtk + 0.5_dp * ( tmp2(1)*tmp2(1) + tmp2(2)*tmp2(2) + tmp2(3)*tmp2(3) )
         END DO
         j = j + stride
      END DO

      j = hbeads - stride/2 + 1
      pk1 = j - stride/2
      DO k = 1, n
         CALL helium_boxmean_3d( helium, pos(:,p(k),pk1), pos(:,perm(p(k)),1), x1 )
         tmp1(:) = x1(:) - pos(:,p(k),j)
         CALL helium_pbc( helium, tmp1 )
         tmp1(:) = tmp1(:) / sigma
         dtk = dtk - 0.5_dp * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
      END DO
      DO k = 1, n
         CALL helium_boxmean_3d( helium, work(:,p(k),pk1), work(:,perm(p(1+MOD(k,n))),1), tmp1 )
         DO c = 1, 3
            x = next_random_number(rng_stream=helium%rng_stream_gaussian,&
                variance=1.0_dp,error=error)
            x = sigma * x
            tmp1(c) = tmp1(c) + x
            tmp2(c) = x
         END DO
         CALL helium_pbc( helium, tmp1 )
         CALL helium_pbc( helium, tmp2 )
         work(:,p(k),j) = tmp1(:)
         tmp2(:) = tmp2(:) / sigma
         dtk = dtk + 0.5_dp * ( tmp2(1)*tmp2(1) + tmp2(2)*tmp2(2) + tmp2(3)*tmp2(3) )
      END DO
      ! ok we got the new positions
      ! calculate action_k(s)-action_k(sprime)
      x = 1.0_dp/(helium%tau*helium%hb2m*stride)
      j = i
      DO
         IF (j>hbeads-stride/2) EXIT
         pk1 = j+stride/2
         DO k = 1, n
            tmp1(:) = pos(:,p(k),j) - pos(:,p(k),pk1)
            CALL helium_pbc( helium, tmp1 )
            ds = ds + x * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
            tmp1(:) = work(:,p(k),j) - work(:,p(k),pk1)
            CALL helium_pbc( helium, tmp1 )
            ds = ds - x * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
            ! interaction change
            IF (helium%solute_present) THEN
               CALL helium_bead_solute_e(pint_env, helium, p(k), pk1, energy=e1)
               CALL helium_bead_solute_e(pint_env, helium, p(k), pk1, work(:,p(k),pk1), e2)
               ds=ds + (stride/2)*(e1-e2) * helium%tau
            END IF
            DO l = 1, nb
               IF (l /= p(k)) THEN
                  tmp1(:) = pos(:,p(k),pk1) - pos(:,l,pk1)
                  CALL helium_pbc( helium, tmp1 )
                  r = tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3)
                  IF ( (r < cell_size) .OR. nperiodic) THEN
                     r = SQRT(r)
                     ds = ds + (stride/2)*helium_spline(uij(1,1)%spline_data,r)
                  END IF
                  tmp1(:) = work(:,p(k),pk1) - work(:,l,pk1)
                  CALL helium_pbc( helium, tmp1 )
                  r = tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3)
                  IF ( (r < cell_size) .OR. nperiodic) THEN
                     r = SQRT(r)
                     ds = ds - (stride/2)*helium_spline(uij(1,1)%spline_data,r)
                  END IF
               END IF
            END DO
            ! counted p[k], p[m] twice. subtract those again
            IF (k < n) THEN
               DO l = k+1, n
                  tmp1(:) = pos(:,p(k),pk1) - pos(:,p(l),pk1)
                  CALL helium_pbc( helium, tmp1 )
                  r = tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3)
                  IF ( (r < cell_size) .OR. nperiodic) THEN
                     r = SQRT(r)
                     ds = ds - (stride/2)*helium_spline(uij(1,1)%spline_data,r)
                  END IF
                  tmp1(:) = work(:,p(k),pk1) - work(:,p(l),pk1)
                  CALL helium_pbc( helium, tmp1 )
                  r = tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3)
                  IF ( (r < cell_size) .OR. nperiodic) THEN
                     r = SQRT(r)
                     ds = ds + (stride/2)*helium_spline(uij(1,1)%spline_data,r)
                  END IF
               END DO
            END IF
         END DO
         j = j + stride/2
      END DO
      ! last link
      pk1 = hbeads - stride/2 + 1
      DO k = 1, n
         tmp1(:) = pos(:,p(k),pk1) - pos(:,perm(p(k)),1)
         CALL helium_pbc( helium, tmp1 )
         ds = ds + x * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
         tmp1(:) = work(:,p(k),pk1) - work(:,perm(p(1+MOD(k,n))),1)
         CALL helium_pbc( helium, tmp1 )
         ds = ds - x * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
      END DO
      ! ok now accept or reject:
      tmp = next_random_number(helium%rng_stream_uniform,error=error)
!      IF ((dtk+ds-pds < 0.0_dp).AND.(EXP(dtk+ds-pds)<tmp)) THEN
      IF (dtk+ds-pds < 0.0_dp) THEN
         IF (EXP(dtk+ds-pds)<tmp) THEN
            res = 0
DO c = i, hbeads
DO k = 1, n
work(:,p(k),c)=pos(:,p(k),c)
END DO
END DO
            RETURN
         END IF
      END IF
      ! accepted. go on to the next level
      helium%num_accepted(level+2,n) = helium%num_accepted(level+2,n) + 1
      level = level+1
      pds = ds
      stride = stride / 2
   END DO
   ! we are on the lowest level now
   ! calc new trial positions
   ! trial action: endpointapprox for T, full action otherwise
   sigma = SQRT(0.5_dp*helium%hb2m*helium%tau)
   dtk = 0.0_dp
   ds = 0.0_dp
   DO j = i+1, hbeads-1, 2
      pk1 = j-1
      pk2 = j+1
      ! calculate log(T(s)):
      DO k = 1, n
         CALL helium_boxmean_3d( helium, pos(:,p(k),pk1), pos(:,p(k),pk2), x1 )
         tmp1(:) = x1(:) - pos(:,p(k),j)
         CALL helium_pbc( helium, tmp1 )
         tmp1(:) = tmp1(:) / sigma
         dtk = dtk - 0.5_dp * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
      END DO
      ! calculate log(T(sprime)) and sprime itself
      DO k = 1, n
         CALL helium_boxmean_3d( helium, work(:,p(k),pk1), work(:,p(k),pk2), tmp1 )
         DO c = 1, 3
            x = next_random_number(rng_stream=helium%rng_stream_gaussian,&
                variance=1.0_dp,error=error)
            x = sigma * x
            tmp1(c) = tmp1(c) + x
            tmp2(c) = x
         END DO
         CALL helium_pbc( helium, tmp1 )
         CALL helium_pbc( helium, tmp2 )
         work(:,p(k),j) = tmp1(:)
         tmp2(:) = tmp2(:) / sigma
         dtk = dtk + 0.5_dp * ( tmp2(1)*tmp2(1) + tmp2(2)*tmp2(2) + tmp2(3)*tmp2(3) )
      END DO
   END DO
   j = hbeads
   pk1 = j-1
   DO k = 1, n
      CALL helium_boxmean_3d( helium, pos(:,p(k),pk1), pos(:,perm(p(k)),1), x1 )
      tmp1(:) = x1(:) - pos(:,p(k),j)
      CALL helium_pbc( helium, tmp1 )
      tmp1(:) = tmp1(:) / sigma
      dtk = dtk - 0.5_dp * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
   END DO
   DO k = 1, n
     CALL helium_boxmean_3d( helium, work(:,p(k),pk1), work(:,perm(p(1+MOD(k,n))),1), tmp1 )
     DO c = 1, 3
        x = next_random_number(rng_stream=helium%rng_stream_gaussian,&
            variance=1.0_dp,error=error)
        x = sigma * x
        tmp1(c) = tmp1(c) + x
        tmp2(c) = x
     END DO
     CALL helium_pbc( helium, tmp1 )
     CALL helium_pbc( helium, tmp2 )
     work(:,p(k),j) = tmp1(:)
     tmp2 = tmp2 / sigma
     dtk = dtk + 0.5_dp * ( tmp2(1)*tmp2(1) + tmp2(2)*tmp2(2) + tmp2(3)*tmp2(3) )
   END DO
   ! ok we got the new positions.
   ! calculate action_k(s)-action_k(sprime)
! interaction change
!TODO interaction ONLY here? or some simplified 12-6 in the upper part?
   IF (helium%solute_present) THEN
     DO j = i , hbeads
       DO k = 1, n
         CALL helium_bead_solute_e(pint_env, helium, p(k), j, energy=e1)
         CALL helium_bead_solute_e(pint_env, helium, p(k), j, work(:,p(k),j), e2)
         ds = ds + (e1-e2) * helium%tau
       END DO
     END DO
   END IF
   x = 1.0_dp/(helium%tau*helium%hb2m*stride)
   DO j = i, hbeads-1
      pk1 = j+1
      DO k = 1, n
         tmp1(:) = pos(:,p(k),j) - pos(:,p(k),pk1)
         CALL helium_pbc( helium, tmp1 )
         ds = ds + x * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
         tmp1(:) = work(:,p(k),j) - work(:,p(k),pk1)
         CALL helium_pbc( helium, tmp1 )
         ds = ds - x * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
         DO l = 1, nb
            IF (l /= p(k)) THEN
               rm1(:)=pos(:,p(k),j)-pos(:,l,j)
               rm2(:)=pos(:,p(k),pk1)-pos(:,l,pk1)
               ds = ds + helium_eval_expansion(helium,rm1,rm2,uij,1)
               rm1(:)=work(:,p(k),j)-work(:,l,j)
               rm2(:)=work(:,p(k),pk1)-work(:,l,pk1)
               ds = ds - helium_eval_expansion(helium,rm1,rm2,uij,1)
            END IF
         END DO
         ! counted p[k], p[m] twice. subtract those again
         IF (k < n) THEN
            DO l = k+1, n
               rm1(:)=pos(:,p(k),j)-pos(:,p(l),j)
               rm2(:)=pos(:,p(k),pk1)-pos(:,p(l),pk1)
               ds = ds - helium_eval_expansion(helium,rm1,rm2,uij,1)
               rm1(:)=work(:,p(k),j)-work(:,p(l),j)
               rm2(:)=work(:,p(k),pk1)-work(:,p(l),pk1)
               ds = ds + helium_eval_expansion(helium,rm1,rm2,uij,1)
            END DO
         END IF
      END DO
   END DO
   j = hbeads
   pk1 = 1
   DO k = 1, n
      tmp1(:) = pos(:,p(k),j) - pos(:,perm(p(k)),1)
      CALL helium_pbc( helium, tmp1 )
      ds = ds + x * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
      tmp1(:) = work(:,p(k),j) - work(:,perm(p(1+MOD(k,n))),1)
      CALL helium_pbc( helium, tmp1 )
      ds = ds - x * ( tmp1(1)*tmp1(1) + tmp1(2)*tmp1(2) + tmp1(3)*tmp1(3) )
      DO l = 1, nb
         IF (l /= p(k)) THEN
            rm1(:)=pos(:,p(k),j)-pos(:,l,j)
            rm2(:)=pos(:,perm(p(k)),1)-pos(:,perm(l),1)
            ds = ds + helium_eval_expansion(helium,rm1,rm2,uij,1)
         END IF
      END DO
      ! counted p[k], p[m] twice. subtract those again
      IF (k < n) THEN
         DO l = k+1, n
            rm1(:)=pos(:,p(k),j)- pos(:,p(l),j)
            rm2(:)=pos(:,perm(p(k)),pk1)-pos(:,perm(p(l)),pk1)
            ds = ds - helium_eval_expansion(helium,rm1,rm2,uij,1)
         END DO
      END IF
   END DO
   IF (n>1) THEN
      !k,c,l
      c = perm(p(1))
      DO k = 1, n - 1
         perm(p(k)) = perm(p(k+1))
      END DO
      perm(p(n)) = c
   END IF
   DO k = 1, n
      DO l = 1, nb
         IF (l /= p(k)) THEN
            rm1(:)=work(:,p(k),j)-work(:,l,j)
            rm2(:)=work(:,perm(p(k)),1)-work(:,perm(l),1)
            ds = ds - helium_eval_expansion(helium,rm1,rm2,uij,1)
         END IF
      END DO
      ! counted p[k], p[m] twice. subtract those again
      IF (k < n) THEN
         DO l = k+1, n
            rm1(:)=work(:,p(k),j)-work(:,p(l),j)
            rm2(:)=work(:,perm(p(k)),1)-work(:,perm(p(l)),1)
            ds = ds + helium_eval_expansion(helium,rm1,rm2,uij,1)
         END DO
      END IF
   END DO
   ! ok now accept or reject:
   tmp = next_random_number(helium%rng_stream_uniform,error=error)
!   IF ((dtk+ds-pds<0.0_dp).AND.(EXP(dtk+ds-pds)<tmp)) THEN
   IF (dtk+ds-pds<0.0_dp) THEN
      IF (EXP(dtk+ds-pds)<tmp) THEN
         res = 0
DO c = i, hbeads
DO k = 1, n
work(:,p(k),c)=pos(:,p(k),c)
END DO
END DO
IF (n>1) THEN
   c = perm(p(n))
   DO k = n - 1, 1, -1
      perm(p(k+1)) = perm(p(k))
   END DO
   perm(p(1)) = c
END IF
         RETURN
      END IF
   END IF
   ! accepted.
   ! copy trial over to the real thing
DO c = i, hbeads
DO k = 1, n
pos(:,p(k),c)=work(:,p(k),c)
END DO
END DO
DO k = 1, n
helium%iperm(perm(p(k))) = p(k)
END DO
!   pos(:,:,:) = work(:,:,:)

    helium%num_accepted(level+2,n) = helium%num_accepted(level+2,n) + 1
    res = 1

    RETURN
  END FUNCTION helium_slice_metro_cyclic

  ! ***************************************************************************
  !> \brief  ...
  !> \author hforbert
  !> \return 1 - if the permutation was successfully selected, 0 - otherwise
  ! ***************************************************************************
  FUNCTION helium_select_permutation(helium,len,error) RESULT(res)
    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: len
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: res

    INTEGER                                  :: i, j, k, n
    INTEGER, DIMENSION(:), POINTER           :: iperm, p, perm
    INTEGER, DIMENSION(:, :), POINTER        :: nmatrix
    REAL(KIND=dp)                            :: rnd, s1, s2, t
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: ipmatrix, pmatrix, tmatrix

    s1 = 0.0_dp
    s2 = 0.0_dp
    res = 0
    n = helium%atoms
    tmatrix => helium%tmatrix
    pmatrix => helium%pmatrix
    ipmatrix => helium%ipmatrix
    perm => helium%permutation
    iperm => helium%iperm
    p => helium%ptable
    nmatrix => helium%nmatrix

    p(len+1)=-1
    rnd = next_random_number(helium%rng_stream_uniform,error=error)
    p(1) = INT(n*rnd)+1
    DO k = 1, len - 1
      t = next_random_number(helium%rng_stream_uniform,error=error)
      ! find the corresponding path to connect to
      ! using the precalculated optimal decision tree:
      i = n-1
      DO
        IF (tmatrix(p(k),i) > t) THEN
          i = nmatrix(p(k),2*i-1)
        ELSE
          i = nmatrix(p(k),2*i)
        END IF
        IF (i<0) EXIT
      END DO
      i = -i
      ! which particle was it previously connected to?
      p(k+1) = iperm(i)
      ! is it unique? quit if it was already part of the permutation
      DO j = 1, k
        IF (p(j) == p(k+1)) RETURN
      END DO
      ! acummulate the needed values for the final
      ! accept/reject step:
      s1 = s1 + ipmatrix(p(k),i)
      s2 = s2 + ipmatrix(p(k),perm(p(k)))
    END DO
    ! close the permutation loop:
    s1 = s1 + ipmatrix(p(len),perm(p(1)))
    s2 = s2 + ipmatrix(p(len),perm(p(len)))
    ! final accept/reject:
    rnd = next_random_number(helium%rng_stream_uniform,error=error)
    t = s1*rnd
    IF (t > s2) RETURN
    ! ok, we have accepted the permutation
    ! calculate the action bias for the subsequent resampling
    ! of the paths:
    s1 = pmatrix(p(len),perm(p(1)))-pmatrix(p(len),perm(p(len)))
    DO k = 1, len-1
      s1 = s1 + pmatrix(p(k),perm(p(k+1))) - pmatrix(p(k),perm(p(k)))
    END DO
    helium%pweight = s1
    res = 1
    RETURN
  END FUNCTION helium_select_permutation

END MODULE helium_sampling
