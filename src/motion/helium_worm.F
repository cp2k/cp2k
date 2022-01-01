!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods dealing with the canonical worm alogrithm
!> \author Felix Uhl [fuhl]
!> \date   2018-10-10
! **************************************************************************************************
MODULE helium_worm

   USE helium_common,                   ONLY: helium_calc_plength,&
                                              helium_eval_chain,&
                                              helium_eval_expansion,&
                                              helium_pbc,&
                                              helium_update_coord_system
   USE helium_interactions,             ONLY: helium_bead_solute_e_f,&
                                              helium_calc_energy,&
                                              helium_solute_e_f
   USE helium_io,                       ONLY: helium_write_line
   USE helium_types,                    ONLY: helium_solvent_type
   USE input_constants,                 ONLY: helium_forces_average,&
                                              helium_forces_last
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE pint_types,                      ONLY: pint_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_worm'

   PUBLIC :: helium_sample_worm

CONTAINS

! **************************************************************************************************
!> \brief Main worm sampling routine
!> \param helium ...
!> \param pint_env ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE helium_sample_worm(helium, pint_env)

      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      TYPE(pint_env_type), INTENT(INOUT), POINTER        :: pint_env

      CHARACTER(len=default_string_length)               :: stmp
      INTEGER :: ac, iatom, ibead, icrawl, iMC, imove, ncentracc, ncentratt, ncloseacc, ncloseatt, &
         ncrawlbwdacc, ncrawlbwdatt, ncrawlfwdacc, ncrawlfwdatt, ncycle, nMC, nmoveheadacc, &
         nmoveheadatt, nmovetailacc, nmovetailatt, nopenacc, nopenatt, npswapacc, nstagacc, &
         nstagatt, nstat, nswapacc, nswapatt, ntot, staging_l
      REAL(KIND=dp)                                      :: rtmp

      CALL helium_update_coord_system(helium, pint_env)

      ncentratt = 0
      ncentracc = 0
      nstagatt = 0
      nstagacc = 0
      nopenatt = 0
      nopenacc = 0
      ncloseatt = 0
      ncloseacc = 0
      nmoveheadatt = 0
      nmoveheadacc = 0
      nmovetailatt = 0
      nmovetailacc = 0
      ncrawlfwdatt = 0
      ncrawlfwdacc = 0
      ncrawlbwdatt = 0
      ncrawlbwdacc = 0
      nswapatt = 0
      npswapacc = 0
      nswapacc = 0
      nstat = 0
      ntot = 0
      helium%proarea%inst(:) = 0.d0
      helium%prarea2%inst(:) = 0.d0
      helium%wnumber%inst(:) = 0.d0
      helium%wnmber2%inst(:) = 0.d0
      helium%mominer%inst(:) = 0.d0
      IF (helium%solute_present) helium%force_avrg(:, :) = 0.0d0
      helium%energy_avrg(:) = 0.0d0
      helium%plength_avrg(:) = 0.0d0
      IF (helium%worm_max_open_cycles > 0) THEN
         helium%savepos(:, :, :) = helium%pos(:, :, :)
         helium%saveiperm(:) = helium%iperm(:)
         helium%savepermutation(:) = helium%permutation(:)
      END IF

      nMC = helium%iter_rot
      ncycle = 0
      IF (helium%worm_allow_open) THEN
         DO ! Exit criterion at the end of the loop
            DO iMC = 1, nMC
               imove = helium%rng_stream_uniform%next(1, helium%worm_all_limit)
               IF (helium%worm_is_closed) THEN
                  IF ((imove >= helium%worm_centroid_min) .AND. (imove <= helium%worm_centroid_max)) THEN
                     ! centroid move
                     iatom = helium%rng_stream_uniform%next(1, helium%atoms)
                     CALL worm_centroid_move(pint_env, helium, &
                                             iatom, helium%worm_centroid_drmax, ac)
                     ncentratt = ncentratt + 1
                     ncentracc = ncentracc + ac
                     ! Note: weights for open and centroid move are taken from open sampling
                     ! staging is adjusted to conserve these weights
                  ELSE IF ((imove >= helium%worm_centroid_max + 1) .AND. (imove <= helium%worm_open_close_min - 1)) THEN
                     ! staging move
                     iatom = helium%rng_stream_uniform%next(1, helium%atoms)
                     ibead = helium%rng_stream_uniform%next(1, helium%beads)
                     staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                     CALL worm_staging_move(pint_env, helium, &
                                            iatom, ibead, staging_l, ac)
                     nstagatt = nstagatt + 1
                     nstagacc = nstagacc + ac
                  ELSE IF ((imove >= helium%worm_open_close_min) .AND. (imove <= helium%worm_open_close_max)) THEN
                     ! attempt opening of worm
                     iatom = helium%rng_stream_uniform%next(1, helium%atoms)
                     ibead = helium%rng_stream_uniform%next(1, helium%beads)
                     staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                     CALL worm_open_move(pint_env, helium, &
                                         iatom, ibead, staging_l, ac)
                     nopenatt = nopenatt + 1
                     nopenacc = nopenacc + ac
                  ELSE
                     ! this must not occur
                     CPABORT("Undefined move selected in helium worm sampling!")
                  END IF
               ELSE ! worm is open
                  IF ((imove >= helium%worm_centroid_min) .AND. (imove <= helium%worm_centroid_max)) THEN
                     ! centroid move
                     iatom = helium%rng_stream_uniform%next(1, helium%atoms)
                     CALL worm_centroid_move(pint_env, helium, &
                                             iatom, helium%worm_centroid_drmax, ac)
                     ncentratt = ncentratt + 1
                     ncentracc = ncentracc + ac
                  ELSE IF ((imove >= helium%worm_staging_min) .AND. (imove <= helium%worm_staging_max)) THEN
                     ! staging move
                     iatom = helium%rng_stream_uniform%next(1, helium%atoms)
                     ibead = helium%rng_stream_uniform%next(1, helium%beads)
                     staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                     CALL worm_staging_move(pint_env, helium, &
                                            iatom, ibead, staging_l, ac)
                     nstagatt = nstagatt + 1
                     nstagacc = nstagacc + ac
                  ELSE IF ((imove >= helium%worm_fcrawl_min) .AND. (imove <= helium%worm_fcrawl_max)) THEN
                     ! crawl forward
                     DO icrawl = 1, helium%worm_repeat_crawl
                        staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                        CALL worm_crawl_move_forward(pint_env, helium, &
                                                     staging_l, ac)
                        ncrawlfwdatt = ncrawlfwdatt + 1
                        ncrawlfwdacc = ncrawlfwdacc + ac
                     END DO
                  ELSE IF ((imove >= helium%worm_bcrawl_min) .AND. (imove <= helium%worm_bcrawl_max)) THEN
                     ! crawl backward
                     DO icrawl = 1, helium%worm_repeat_crawl
                        staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                        CALL worm_crawl_move_backward(pint_env, helium, &
                                                      staging_l, ac)
                        ncrawlbwdatt = ncrawlbwdatt + 1
                        ncrawlbwdacc = ncrawlbwdacc + ac
                     END DO
                  ELSE IF ((imove >= helium%worm_head_min) .AND. (imove <= helium%worm_head_max)) THEN
                     ! move head
                     staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                     CALL worm_head_move(pint_env, helium, &
                                         staging_l, ac)
                     nmoveheadatt = nmoveheadatt + 1
                     nmoveheadacc = nmoveheadacc + ac
                  ELSE IF ((imove >= helium%worm_tail_min) .AND. (imove <= helium%worm_tail_max)) THEN
                     ! move tail
                     staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                     CALL worm_tail_move(pint_env, helium, &
                                         staging_l, ac)
                     nmovetailatt = nmovetailatt + 1
                     nmovetailacc = nmovetailacc + ac
                  ELSE IF ((imove >= helium%worm_swap_min) .AND. (imove <= helium%worm_swap_max)) THEN
                     staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                     CALL worm_swap_move(pint_env, helium, &
                                         helium%atoms, staging_l, ac)
                     npswapacc = npswapacc + ac
                     nswapacc = nswapacc + ac
                     nswapatt = nswapatt + 1
                  ELSE IF ((imove >= helium%worm_open_close_min) .AND. (imove <= helium%worm_open_close_max)) THEN
                     ! attempt closing of worm
                     staging_l = helium%rng_stream_uniform%next(2, helium%worm_staging_l)
                     CALL worm_close_move(pint_env, helium, &
                                          staging_l, ac)
                     ncloseatt = ncloseatt + 1
                     ncloseacc = ncloseacc + ac
                  ELSE
                     ! this must not occur
                     CPABORT("Undefined move selected in helium worm sampling!")
                  END IF
               END IF

               ! Accumulate statistics if we are in the Z-sector:
               IF (helium%worm_is_closed) THEN
                  nstat = nstat + 1
                  IF (helium%solute_present) THEN
                     IF (helium%get_helium_forces == helium_forces_average) THEN
                        !TODO needs proper averaging!
                        CALL helium_solute_e_f(pint_env, helium, rtmp)
                        helium%force_avrg(:, :) = helium%force_avrg(:, :) + helium%force_inst(:, :)
                     END IF
                  END IF
               END IF
               ntot = ntot + 1
            END DO ! MC loop

            IF (helium%worm_is_closed) THEN
               EXIT
            ELSE
               ncycle = ncycle + 1
               ! reset position and permutation and start MC cycle again
               ! if stuck in open position for too long
               IF (helium%worm_max_open_cycles > 0 .AND. ncycle > helium%worm_max_open_cycles) THEN
                  nMC = helium%iter_rot
                  ncycle = 0
                  helium%pos = helium%savepos
                  helium%permutation = helium%savepermutation
                  helium%iperm = helium%saveiperm
                  CPWARN("Trapped in open state, reset helium to closed state from last MD step.")
               ELSE
                  nMC = MAX(helium%iter_rot/10, 10)
               END IF
            END IF
         END DO !attempts loop
      ELSE ! only closed configurations allowed
         DO iMC = 1, nMC
            imove = helium%rng_stream_uniform%next(1, helium%worm_all_limit)

            IF ((imove >= helium%worm_centroid_min) .AND. (imove <= helium%worm_centroid_max)) THEN
               ! centroid move
               iatom = helium%rng_stream_uniform%next(1, helium%atoms)
               CALL worm_centroid_move(pint_env, helium, &
                                       iatom, helium%worm_centroid_drmax, ac)
               ncentratt = ncentratt + 1
               ncentracc = ncentracc + ac
            ELSE IF ((imove >= helium%worm_staging_min) .AND. (imove <= helium%worm_staging_max)) THEN
               ! staging move
               iatom = helium%rng_stream_uniform%next(1, helium%atoms)
               ibead = helium%rng_stream_uniform%next(1, helium%beads)
               CALL worm_staging_move(pint_env, helium, &
                                      iatom, ibead, helium%worm_staging_l, ac)
               nstagatt = nstagatt + 1
               nstagacc = nstagacc + ac
            ELSE
               ! this must not occur
               CPABORT("Undefined move selected in helium worm sampling!")
            END IF

            ! Accumulate statistics if we are in closed configurations (which we always are)
            nstat = nstat + 1
            ntot = ntot + 1
            IF (helium%solute_present) THEN
               IF (helium%get_helium_forces == helium_forces_average) THEN
                  ! TODO: needs proper averaging
                  CALL helium_solute_e_f(pint_env, helium, rtmp)
                  helium%force_avrg(:, :) = helium%force_avrg(:, :) + helium%force_inst(:, :)
               END IF
            END IF
         END DO ! MC loop
      END IF

      ! Save naccepted and ntot
      helium%num_accepted(1, 1) = ncentracc + nstagacc + nopenacc + ncloseacc + nswapacc + &
                                  nmoveheadacc + nmovetailacc + ncrawlfwdacc + ncrawlbwdacc
      helium%num_accepted(2, 1) = ncentratt + nstagatt + nopenatt + ncloseatt + nswapatt + &
                                  nmoveheadatt + nmovetailatt + ncrawlfwdatt + ncrawlbwdatt
      helium%worm_nstat = nstat

      ! Calculate energy and permutation path length
      ! Multiply times helium%worm_nstat for consistent averaging in helium_sampling
      CALL helium_calc_energy(helium, pint_env)
      helium%energy_avrg(:) = helium%energy_inst(:)*REAL(helium%worm_nstat, dp)

      CALL helium_calc_plength(helium)
      helium%plength_avrg(:) = helium%plength_inst(:)*REAL(helium%worm_nstat, dp)

      ! Calculate last_force
      IF (helium%solute_present) THEN
         IF (helium%get_helium_forces == helium_forces_last) THEN
            CALL helium_solute_e_f(pint_env, helium, rtmp)
            helium%force_avrg(:, :) = helium%force_avrg(:, :) + helium%force_inst(:, :)
         END IF
      END IF

      IF (helium%worm_show_statistics) THEN
         WRITE (stmp, *) '--------------------------------------------------'
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Opening: ', &
            REAL(nopenacc, dp)/REAL(MAX(1, nopenatt), dp), &
            nopenacc, nopenatt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Closing: ', &
            REAL(ncloseacc, dp)/REAL(MAX(1, ncloseatt), dp), &
            ncloseacc, ncloseatt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Move Head: ', &
            REAL(nmoveheadacc, dp)/REAL(MAX(1, nmoveheadatt), dp), &
            nmoveheadacc, nmoveheadatt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Move Tail: ', &
            REAL(nmovetailacc, dp)/REAL(MAX(1, nmovetailatt), dp), &
            nmovetailacc, nmovetailatt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Crawl FWD: ', &
            REAL(ncrawlfwdacc, dp)/REAL(MAX(1, ncrawlfwdatt), dp), &
            ncrawlfwdacc, ncrawlfwdatt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Crawl BWD: ', &
            REAL(ncrawlbwdacc, dp)/REAL(MAX(1, ncrawlbwdatt), dp), &
            ncrawlbwdacc, ncrawlbwdatt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Staging: ', &
            REAL(nstagacc, dp)/REAL(MAX(1, nstagatt), dp), &
            nstagacc, nstagatt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Centroid: ', &
            REAL(ncentracc, dp)/REAL(MAX(1, ncentratt), dp), &
            ncentracc, ncentratt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Select: ', &
            REAL(npswapacc, dp)/REAL(MAX(1, nswapatt), dp), &
            npswapacc, nswapatt
         CALL helium_write_line(stmp)
         WRITE (stmp, '(A10,F15.5,2I10)') 'Swapping: ', &
            REAL(nswapacc, dp)/REAL(MAX(1, nswapatt), dp), &
            nswapacc, nswapatt
         CALL helium_write_line(stmp)
         WRITE (stmp, *) "Open State Probability:   ", REAL(ntot - nstat, dp)/REAL(MAX(1, ntot), dp), ntot - nstat, ntot
         CALL helium_write_line(stmp)
         WRITE (stmp, *) "Closed State Probability: ", REAL(nstat, dp)/REAL(MAX(1, ntot), dp), nstat, ntot
         CALL helium_write_line(stmp)
      END IF

      !CALL center_pos(helium)

      RETURN
   END SUBROUTINE helium_sample_worm

! **************************************************************************************************
!> \brief Centroid Move (accounting for exchanged particles)
!> \param pint_env ...
!> \param helium ...
!> \param iatom ...
!> \param drmax ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_centroid_move(pint_env, helium, iatom, drmax, ac)

      TYPE(pint_env_type), INTENT(INOUT), POINTER        :: pint_env
      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      INTEGER, INTENT(IN)                                :: iatom
      REAL(KIND=dp), INTENT(IN)                          :: drmax
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: ia, ib, ic, jatom
      LOGICAL                                            :: should_reject, worm_in_moved_cycle
      REAL(KIND=dp)                                      :: rtmp, rtmpo, sdiff, snew, sold
      REAL(KIND=dp), DIMENSION(3)                        :: dr, dro, new_com, old_com

      DO ic = 1, 3
         rtmp = helium%rng_stream_uniform%next()
         dr(ic) = (2.0_dp*rtmp - 1.0_dp)*drmax
      END DO

      IF (helium%worm_is_closed) THEN
         worm_in_moved_cycle = .FALSE.
         ! Perform move for first atom
         DO ib = 1, helium%beads
            helium%work(:, iatom, ib) = helium%work(:, iatom, ib) + dr(:)
         END DO
         ! move along permutation cycle
         jatom = helium%permutation(iatom)
         DO WHILE (jatom /= iatom)
            DO ib = 1, helium%beads
               helium%work(:, jatom, ib) = helium%work(:, jatom, ib) + dr(:)
            END DO
            ! next atom in chain
            jatom = helium%permutation(jatom)
         END DO
      ELSE ! worm is open
         worm_in_moved_cycle = (helium%worm_atom_idx == iatom)
         ! while moving, check if worm is in moved cycle
         ! Perform move for first atom
         DO ib = 1, helium%beads
            helium%work(:, iatom, ib) = helium%work(:, iatom, ib) + dr(:)
         END DO
         ! move along permutation cycle
         jatom = helium%permutation(iatom)
         DO WHILE (jatom /= iatom)
            DO ib = 1, helium%beads
               helium%work(:, jatom, ib) = helium%work(:, jatom, ib) + dr(:)
            END DO
            worm_in_moved_cycle = worm_in_moved_cycle .OR. (helium%worm_atom_idx == jatom)
            ! next atom in chain
            jatom = helium%permutation(jatom)
         END DO
         ! if atom contains had bead move that as well
         IF (worm_in_moved_cycle) THEN
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:) + dr(:)
         END IF
      END IF

      sold = worm_centroid_move_action(helium, helium%pos, iatom, &
                                       helium%worm_xtra_bead, worm_in_moved_cycle)

      snew = worm_centroid_move_action(helium, helium%work, iatom, &
                                       helium%worm_xtra_bead_work, worm_in_moved_cycle)

      IF (helium%solute_present) THEN
         sold = sold + worm_centroid_move_inter_action(pint_env, helium, helium%pos, iatom, &
                                                       helium%worm_xtra_bead, worm_in_moved_cycle)
         snew = snew + worm_centroid_move_inter_action(pint_env, helium, helium%work, iatom, &
                                                       helium%worm_xtra_bead_work, worm_in_moved_cycle)
      END IF

      ! Metropolis:
      sdiff = sold - snew
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            DO ib = 1, helium%beads
               helium%work(:, iatom, ib) = helium%pos(:, iatom, ib)
            END DO
            ! move along permutation cycle
            jatom = helium%permutation(iatom)
            DO WHILE (jatom /= iatom)
               DO ib = 1, helium%beads
                  helium%work(:, jatom, ib) = helium%pos(:, jatom, ib)
               END DO
               ! next atom in chain
               jatom = helium%permutation(jatom)
            END DO
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! restore original coordinates
            ! write back only changed atoms
            DO ib = 1, helium%beads
               helium%work(:, iatom, ib) = helium%pos(:, iatom, ib)
            END DO
            ! move along permutation cycle
            jatom = helium%permutation(iatom)
            DO WHILE (jatom /= iatom)
               DO ib = 1, helium%beads
                  helium%work(:, jatom, ib) = helium%pos(:, jatom, ib)
               END DO
               ! next atom in chain
               jatom = helium%permutation(jatom)
            END DO
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      ! write changed coordinates to position array
      DO ib = 1, helium%beads
         helium%pos(:, iatom, ib) = helium%work(:, iatom, ib)
      END DO
      ! move along permutation cycle
      jatom = helium%permutation(iatom)
      DO WHILE (jatom /= iatom)
         DO ib = 1, helium%beads
            helium%pos(:, jatom, ib) = helium%work(:, jatom, ib)
         END DO
         ! next atom in chain
         jatom = helium%permutation(jatom)
      END DO
      helium%worm_xtra_bead(:) = helium%worm_xtra_bead_work(:)

      RETURN

   END SUBROUTINE worm_centroid_move

! **************************************************************************************************
!> \brief Action for centroid move
!> \param helium ...
!> \param pos ...
!> \param iatom ...
!> \param xtrapos ...
!> \param winc ...
!> \return ...
!> \author fuhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION worm_centroid_move_action(helium, pos, iatom, xtrapos, winc) &
      RESULT(partaction)

      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: pos
      INTEGER, INTENT(IN)                                :: iatom
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xtrapos
      LOGICAL, INTENT(IN)                                :: winc

      INTEGER                                            :: ia, ib, jatom, katom, opatom, patom, &
                                                            wbead
      LOGICAL                                            :: incycle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work2, work3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), DIMENSION(3)                        :: r, rp

      ALLOCATE (work(3, helium%beads + 1), work2(helium%beads + 1), work3(SIZE(helium%uoffdiag, 1) + 1))

      ! compute action difference
      ! action before the move from pos array
      partaction = 0.0_dp

      ! pair of first atom with non-in-cycle atoms
      jatom = iatom
      DO ia = 1, helium%atoms
         IF (ia == jatom) CYCLE
         incycle = .FALSE.
         katom = helium%permutation(jatom)
         DO WHILE (katom /= jatom)
            IF (katom == ia) THEN
               incycle = .TRUE.
               EXIT
            END IF
            katom = helium%permutation(katom)
         END DO
         IF (incycle) CYCLE
         ! if not in cycle, compute pair action
         DO ib = 1, helium%beads
            work(:, ib) = pos(:, ia, ib) - pos(:, jatom, ib)
         END DO
         work(:, helium%beads + 1) = pos(:, helium%permutation(ia), 1) - &
                                     pos(:, helium%permutation(jatom), 1)
         partaction = partaction + helium_eval_chain(helium, work, helium%beads + 1, work2, work3)
      END DO
      ! all other cycle atoms with non-in-cycle atoms
      jatom = helium%permutation(iatom)
      DO WHILE (jatom /= iatom)
         DO ia = 1, helium%atoms
            IF (ia == jatom) CYCLE
            incycle = .FALSE.
            katom = helium%permutation(jatom)
            DO WHILE (katom /= jatom)
               IF (katom == ia) THEN
                  incycle = .TRUE.
                  EXIT
               END IF
               katom = helium%permutation(katom)
            END DO
            IF (incycle) CYCLE
            ! if not in cycle, compute pair action
            DO ib = 1, helium%beads
               work(:, ib) = pos(:, ia, ib) - pos(:, jatom, ib)
            END DO
            work(:, helium%beads + 1) = pos(:, helium%permutation(ia), 1) - &
                                        pos(:, helium%permutation(jatom), 1)
            partaction = partaction + helium_eval_chain(helium, work, helium%beads + 1, work2, work3)
         END DO
         jatom = helium%permutation(jatom)
      END DO
      ! correct pair action for open worm configurations
      IF (.NOT. helium%worm_is_closed) THEN
         wbead = helium%worm_bead_idx
         IF (winc) THEN
            IF (helium%worm_bead_idx == 1) THEN
               ! patom is the atom in front of the lone head bead
               patom = helium%iperm(helium%worm_atom_idx)
               ! go through all atoms, not in the cycle, and correct pair action
               DO ia = 1, helium%atoms
                  IF (ia == helium%worm_atom_idx) CYCLE
                  incycle = .FALSE.
                  katom = helium%permutation(helium%worm_atom_idx)
                  DO WHILE (katom /= helium%worm_atom_idx)
                     IF (katom == ia) THEN
                        incycle = .TRUE.
                        EXIT
                     END IF
                     katom = helium%permutation(katom)
                  END DO
                  IF (incycle) CYCLE
                  ! find other previous atom
                  ! opatom is the atom in front of the first bead of the second atom
                  opatom = helium%iperm(ia)
                  ! if not in cycle, compute pair action
                  ! subtract pair action for closed link
                  r(:) = pos(:, helium%worm_atom_idx, 1) - pos(:, ia, 1)
                  rp(:) = pos(:, patom, helium%beads) - pos(:, opatom, helium%beads)
                  partaction = partaction - helium_eval_expansion(helium, r, rp, work3)
                  ! and add corrected extra link
                  ! rp stays the same
                  r(:) = xtrapos(:) - pos(:, ia, 1)
                  partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
               END DO
            ELSE ! bead index /= 1
               ! atom index is constant
               ! go through all atoms, not in the cycle, and correct pair action
               DO ia = 1, helium%atoms
                  IF (ia == helium%worm_atom_idx) CYCLE
                  incycle = .FALSE.
                  katom = helium%permutation(helium%worm_atom_idx)
                  DO WHILE (katom /= helium%worm_atom_idx)
                     IF (katom == ia) THEN
                        incycle = .TRUE.
                        EXIT
                     END IF
                     katom = helium%permutation(katom)
                  END DO
                  IF (incycle) CYCLE
                  ! if not in cycle, compute pair action
                  ! subtract pair action for closed link
                  r(:) = pos(:, helium%worm_atom_idx, wbead) - pos(:, ia, wbead)
                  rp(:) = pos(:, helium%worm_atom_idx, wbead - 1) - pos(:, ia, wbead - 1)
                  partaction = partaction - helium_eval_expansion(helium, r, rp, work3)
                  ! and add corrected extra link
                  ! rp stays the same
                  r(:) = xtrapos(:) - pos(:, ia, wbead)
                  partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
               END DO
            END IF
         ELSE ! worm is not the moved cycle
            IF (helium%worm_bead_idx == 1) THEN
               ! patom is the atom in front of the lone head bead
               patom = helium%iperm(helium%worm_atom_idx)
               opatom = helium%iperm(iatom)
               !correct action contribution for first atom in moved cycle
               r(:) = pos(:, helium%worm_atom_idx, 1) - pos(:, iatom, 1)
               rp(:) = pos(:, patom, helium%beads) - pos(:, opatom, helium%beads)
               partaction = partaction - helium_eval_expansion(helium, r, rp, work3)
               ! and add corrected extra link
               ! rp stays the same
               r(:) = xtrapos(:) - pos(:, iatom, 1)
               partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
               ! go through all other atoms, not in the exchange cycle, and correct pair action
               ia = helium%permutation(iatom)
               DO WHILE (ia /= iatom)
                  opatom = helium%iperm(ia)
                  ! subtract pair action for closed link
                  r(:) = pos(:, helium%worm_atom_idx, 1) - pos(:, ia, 1)
                  rp(:) = pos(:, patom, helium%beads) - pos(:, opatom, helium%beads)
                  partaction = partaction - helium_eval_expansion(helium, r, rp, work3)
                  ! and add corrected extra link
                  ! rp stays the same
                  r(:) = xtrapos(:) - pos(:, ia, 1)
                  partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
                  ia = helium%permutation(ia)
               END DO
            ELSE ! bead index /= 1
               ! patom is the atom in front of the lone head bead
               !correct action contribution for first atom in moved cycle
               r(:) = pos(:, helium%worm_atom_idx, wbead) - pos(:, iatom, wbead)
               rp(:) = pos(:, helium%worm_atom_idx, wbead - 1) - pos(:, iatom, wbead - 1)
               partaction = partaction - helium_eval_expansion(helium, r, rp, work3)
               ! and add corrected extra link
               ! rp stays the same
               r(:) = xtrapos(:) - pos(:, iatom, wbead)
               partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
               ! go through all other atoms, not in the exchange cycle, and correct pair action
               ia = helium%permutation(iatom)
               DO WHILE (ia /= iatom)
                  ! subtract pair action for closed link
                  r(:) = pos(:, helium%worm_atom_idx, wbead) - pos(:, ia, wbead)
                  rp(:) = pos(:, helium%worm_atom_idx, wbead - 1) - pos(:, ia, wbead - 1)
                  partaction = partaction - helium_eval_expansion(helium, r, rp, work3)
                  ! and add corrected extra link
                  ! rp stays the same
                  r(:) = xtrapos(:) - pos(:, ia, wbead)
                  partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
                  ia = helium%permutation(ia)
               END DO
            END IF
         END IF
      END IF
      DEALLOCATE (work, work2, work3)

   END FUNCTION worm_centroid_move_action

! **************************************************************************************************
!> \brief Centroid interaction
!> \param pint_env ...
!> \param helium ...
!> \param pos ...
!> \param iatom ...
!> \param xtrapos ...
!> \param winc ...
!> \return ...
!> \author fuhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION worm_centroid_move_inter_action(pint_env, helium, pos, iatom, xtrapos, winc) &
      RESULT(partaction)

      TYPE(pint_env_type), INTENT(IN), POINTER           :: pint_env
      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: pos
      INTEGER, INTENT(IN)                                :: iatom
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xtrapos
      LOGICAL, INTENT(IN)                                :: winc

      INTEGER                                            :: jatom, jbead
      REAL(KIND=dp)                                      :: energy

      partaction = 0.0_dp
      ! spcial treatment if worm is in moved cycle
      IF (winc) THEN
         ! first atom by hand
         jatom = iatom
         ! if it is worm atom it gets special treatment
         IF (jatom == helium%worm_atom_idx) THEN
            ! up to worm intersection
            DO jbead = 1, helium%worm_bead_idx - 1
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           jatom, jbead, pos(:, jatom, jbead), energy=energy)
               partaction = partaction + energy
            END DO
            ! head and tail each with 1/2 weight
            jbead = helium%worm_bead_idx
            ! tail
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        jatom, jbead, pos(:, jatom, jbead), energy=energy)
            partaction = partaction + 0.5_dp*energy
            ! head
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        jatom, jbead, xtrapos, energy=energy)
            partaction = partaction + 0.5_dp*energy
            ! rest of ring polymer
            DO jbead = helium%worm_bead_idx + 1, helium%beads
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           jatom, jbead, pos(:, jatom, jbead), energy=energy)
               partaction = partaction + energy
            END DO
         ELSE
            DO jbead = 1, helium%beads
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           jatom, jbead, pos(:, jatom, jbead), energy=energy)
               partaction = partaction + energy
            END DO
         END IF
         jatom = helium%permutation(iatom)
         DO WHILE (jatom /= iatom)
            IF (jatom == helium%worm_atom_idx) THEN
               ! up to worm intersection
               DO jbead = 1, helium%worm_bead_idx - 1
                  CALL helium_bead_solute_e_f(pint_env, helium, &
                                              jatom, jbead, pos(:, jatom, jbead), energy=energy)
                  partaction = partaction + energy
               END DO
               ! head and tail each with 1/2 weight
               jbead = helium%worm_bead_idx
               ! tail
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           jatom, jbead, pos(:, jatom, jbead), energy=energy)
               partaction = partaction + 0.5_dp*energy
               ! head
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           jatom, jbead, xtrapos, energy=energy)
               partaction = partaction + 0.5_dp*energy
               ! rest of ring polymer
               DO jbead = helium%worm_bead_idx + 1, helium%beads
                  CALL helium_bead_solute_e_f(pint_env, helium, &
                                              jatom, jbead, pos(:, jatom, jbead), energy=energy)
                  partaction = partaction + energy
               END DO
            ELSE
               DO jbead = 1, helium%beads
                  CALL helium_bead_solute_e_f(pint_env, helium, &
                                              jatom, jbead, pos(:, jatom, jbead), energy=energy)
                  partaction = partaction + energy
               END DO
            END IF
            jatom = helium%permutation(jatom)
         END DO
      ELSE ! worm not in moved cycle
         ! first atom by hand
         jatom = iatom
         ! if it is worm atom it gets special treatment
         DO jbead = 1, helium%beads
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        jatom, jbead, pos(:, jatom, jbead), energy=energy)
            partaction = partaction + energy
         END DO
         jatom = helium%permutation(iatom)
         DO WHILE (jatom /= iatom)
            DO jbead = 1, helium%beads
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           jatom, jbead, pos(:, jatom, jbead), energy=energy)
               partaction = partaction + energy
            END DO
            jatom = helium%permutation(jatom)
         END DO
      END IF

      partaction = partaction*helium%tau

   END FUNCTION worm_centroid_move_inter_action

! **************************************************************************************************
!> \brief General path construct based on staging moves
!> \param helium ...
!> \param ri ...
!> \param rj ...
!> \param l ...
!> \param new_path ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE path_construct(helium, ri, rj, l, new_path)

      !constructs a path of length l between the positions ri and rj
      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ri, rj
      INTEGER                                            :: l
      REAL(KIND=dp), DIMENSION(3, l), INTENT(OUT)        :: new_path

      INTEGER                                            :: idim, istage
      REAL(KIND=dp)                                      :: imass, invstagemass, rk, weight
      REAL(KIND=dp), DIMENSION(3)                        :: re, rs

      imass = 1.0_dp/helium%he_mass_au
      ! dealing with periodicity
      rs(:) = ri(:)
      re(:) = rj(:) - rs(:)
      CALL helium_pbc(helium, re)
      re(:) = re(:) + rs(:)

      ! first construction by hand
      ! reusable weight factor 1/(l+1)
      rk = REAL(l, dp)
      weight = 1.0_dp/(rk + 1.0_dp)
      ! staging mass needed for modified variance
      invstagemass = rk*weight*imass
      ! proposing new positions
      DO idim = 1, 3
         new_path(idim, 1) = helium%rng_stream_gaussian%next(variance=helium%tau*invstagemass)
      END DO
      new_path(:, 1) = new_path(:, 1) + weight*(re(:) + rk*rs(:))

      DO istage = 2, l
         ! reusable weight factor 1/(k+1)
         rk = REAL(l - istage + 1, dp)
         weight = 1.0_dp/(rk + 1.0_dp)
         ! staging mass needed for modified variance
         invstagemass = rk*weight*imass
         ! proposing new positions
         DO idim = 1, 3
            new_path(idim, istage) = helium%rng_stream_gaussian%next(variance=helium%tau*invstagemass)
         END DO
         new_path(:, istage) = new_path(:, istage) + weight*(rk*new_path(:, istage - 1) + re(:))
      END DO

      RETURN

   END SUBROUTINE path_construct

! **************************************************************************************************
!> \brief Staging move
!> \param pint_env ...
!> \param helium ...
!> \param startatom ...
!> \param startbead ...
!> \param l ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_staging_move(pint_env, helium, startatom, startbead, l, ac)

      TYPE(pint_env_type), INTENT(INOUT), POINTER        :: pint_env
      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      INTEGER, INTENT(IN)                                :: startatom, startbead, l
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: endatom, endbead, ia, ib, ibead, jbead
      LOGICAL                                            :: should_reject, worm_overlap
      REAL(KIND=dp)                                      :: rtmp, rtmpo, sdiff, snew, sold
      REAL(KIND=dp), DIMENSION(3)                        :: dr, dro, new_com, old_com
      REAL(KIND=dp), DIMENSION(3, l)                     :: newsection

      ac = 0
      endbead = startbead + l + 1
      ! Check if the imaginary time section belongs to two atoms
      IF (endbead > helium%beads) THEN
         endatom = helium%permutation(startatom)
         endbead = endbead - helium%beads
      ELSE
         endatom = startatom
      END IF

      ! check if the imaginary time section overlaps with the worm opening
      IF (helium%worm_is_closed) THEN
         worm_overlap = .FALSE.
      ELSE
         ! if it does give it special treatment during action evaluation
         worm_overlap = ((startbead < endbead) .AND. &
                         (helium%worm_bead_idx > startbead) .AND. &
                         (helium%worm_bead_idx <= endbead)) &
                        .OR. &
                        ((startbead > endbead) .AND. &
                         ((helium%worm_bead_idx > startbead) .OR. &
                          (helium%worm_bead_idx <= endbead)))
         IF (worm_overlap) THEN
            ! if in addition the worm end is IN the reconstructed section reject immediately
            IF (((helium%worm_atom_idx == startatom) .AND. (helium%worm_bead_idx >= startbead)) .OR. &
                ((helium%worm_atom_idx == endatom) .AND. (helium%worm_bead_idx <= endbead))) THEN
               ac = 0
               RETURN
            END IF
         END IF
      END IF
      ! action before the move
      IF (worm_overlap) THEN
         sold = worm_path_action_worm_corrected(helium, helium%pos, &
                                                startatom, startbead, endatom, endbead, &
                                                helium%worm_xtra_bead, helium%worm_atom_idx, helium%worm_bead_idx)
      ELSE
         sold = worm_path_action(helium, helium%pos, &
                                 startatom, startbead, endatom, endbead)
      END IF

      IF (helium%solute_present) THEN
         ! no special head treatment needed, because a swap can't go over
         ! the worm gap and due to primitive coupling no cross bead terms are considered
         sold = sold + worm_path_inter_action(pint_env, helium, helium%pos, &
                                              startatom, startbead, endatom, endbead)
      END IF

      ! construct a new path connecting the start and endbead
      CALL path_construct(helium, &
                          helium%pos(:, startatom, startbead), &
                          helium%pos(:, endatom, endbead), l, &
                          newsection)

      ! write new path segment to work array
      ! first the part that is guaranteed to fit on the coorinates of startatom
      jbead = 1
      DO ibead = startbead + 1, MIN(helium%beads, startbead + l)
         helium%work(:, startatom, ibead) = newsection(:, jbead)
         jbead = jbead + 1
      END DO
      ! transfer the rest of the beads to coordinates of endatom if necessary
      IF (helium%beads < startbead + l) THEN
         DO ibead = 1, endbead - 1
            helium%work(:, endatom, ibead) = newsection(:, jbead)
            jbead = jbead + 1
         END DO
      END IF

      ! action after the move
      IF (worm_overlap) THEN
         snew = worm_path_action_worm_corrected(helium, helium%work, &
                                                startatom, startbead, endatom, endbead, &
                                                helium%worm_xtra_bead_work, helium%worm_atom_idx_work, helium%worm_bead_idx_work)
      ELSE
         snew = worm_path_action(helium, helium%work, &
                                 startatom, startbead, endatom, endbead)
      END IF

      IF (helium%solute_present) THEN
         ! no special head treatment needed, because a swap can't go over
         ! the worm gap and due to primitive coupling no cross bead terms are considered
         snew = snew + worm_path_inter_action(pint_env, helium, helium%work, &
                                              startatom, startbead, endatom, endbead)
      END IF

      ! Metropolis:
      sdiff = sold - snew
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            jbead = 1
            DO ibead = startbead + 1, MIN(helium%beads, startbead + l)
               helium%work(:, startatom, ibead) = helium%pos(:, startatom, ibead)
               jbead = jbead + 1
            END DO
            ! transfer the rest of the beads to coordinates of endatom if necessary
            IF (helium%beads < startbead + l) THEN
               DO ibead = 1, endbead - 1
                  helium%work(:, endatom, ibead) = helium%pos(:, endatom, ibead)
                  jbead = jbead + 1
               END DO
            END IF
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! restore original coordinates
            ! write back only changed atoms
            jbead = 1
            DO ibead = startbead + 1, MIN(helium%beads, startbead + l)
               helium%work(:, startatom, ibead) = helium%pos(:, startatom, ibead)
               jbead = jbead + 1
            END DO
            ! transfer the rest of the beads to coordinates of endatom if necessary
            IF (helium%beads < startbead + l) THEN
               DO ibead = 1, endbead - 1
                  helium%work(:, endatom, ibead) = helium%pos(:, endatom, ibead)
                  jbead = jbead + 1
               END DO
            END IF
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      ! write changed coordinates to position array
      jbead = 1
      DO ibead = startbead + 1, MIN(helium%beads, startbead + l)
         helium%pos(:, startatom, ibead) = helium%work(:, startatom, ibead)
         jbead = jbead + 1
      END DO
      ! transfer the rest of the beads to coordinates of endatom if necessary
      IF (helium%beads < startbead + l) THEN
         DO ibead = 1, endbead - 1
            helium%pos(:, endatom, ibead) = helium%work(:, endatom, ibead)
            jbead = jbead + 1
         END DO
      END IF

   END SUBROUTINE worm_staging_move

! **************************************************************************************************
!> \brief Open move to off-diagonal elements of the density matrix, allows to sample permutations
!> \param pint_env ...
!> \param helium ...
!> \param endatom ...
!> \param endbead ...
!> \param l ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_open_move(pint_env, helium, endatom, endbead, l, ac)

      TYPE(pint_env_type), INTENT(INOUT), POINTER        :: pint_env
      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      INTEGER, INTENT(IN)                                :: endatom, endbead, l
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: ia, ib, idim, kbead, startatom, startbead
      LOGICAL                                            :: should_reject
      REAL(KIND=dp)                                      :: distsq, mass, rtmp, rtmpo, sdiff, snew, &
                                                            sold, xr
      REAL(KIND=dp), DIMENSION(3)                        :: distvec, dr, dro, new_com, old_com

      mass = helium%he_mass_au

      ! get index of the atom and bead, where the resampling of the head begins
      IF (l < endbead) THEN
         ! startbead belongs to the same atom
         startatom = endatom
         startbead = endbead - l
      ELSE
         ! startbead belongs to a different atom
         ! find previous atom (assuming l < nbeads)
         startatom = helium%iperm(endatom)
         startbead = endbead + helium%beads - l
      END IF
      sold = worm_path_action(helium, helium%pos, &
                              startatom, startbead, endatom, endbead)

      IF (helium%solute_present) THEN
         ! yes this is correct, as the bead, that splits into tail and head only changes half
         ! therefore only half of its action needs to be considred
         ! and is cheated in here by passing it as head bead
         sold = sold + worm_path_inter_action_head(pint_env, helium, helium%pos, &
                                                   startatom, startbead, &
                                                   helium%pos(:, endatom, endbead), endatom, endbead)
      END IF

      helium%worm_is_closed = .FALSE.
      helium%worm_atom_idx_work = endatom
      helium%worm_bead_idx_work = endbead

      ! alternative grow with consecutive gaussians
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         ! gro head from startbead
         DO kbead = startbead + 1, endbead - 1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, startatom, kbead) = helium%work(idim, startatom, kbead - 1) + xr
            END DO
         END DO
         ! last grow head bead
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, startatom, endbead - 1) + xr
         END DO
      ELSE IF (endbead /= 1) THEN
         ! is distributed among two atoms
         ! grow from startbead
         DO kbead = startbead + 1, helium%beads
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, startatom, kbead) = helium%work(idim, startatom, kbead - 1) + xr
            END DO
         END DO
         ! bead one of endatom relative to last on startatom
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%work(idim, endatom, 1) = helium%work(idim, startatom, helium%beads) + xr
         END DO
         ! everything on endatom
         DO kbead = 2, endbead - 1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, endatom, kbead) = helium%work(idim, endatom, kbead - 1) + xr
            END DO
         END DO
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, endatom, endbead - 1) + xr
         END DO
      ELSE ! imagtimewrap and headbead = 1
         ! is distributed among two atoms
         ! grow from startbead
         DO kbead = startbead + 1, helium%beads
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, startatom, kbead) = helium%work(idim, startatom, kbead - 1) + xr
            END DO
         END DO
         ! bead one of endatom relative to last on startatom
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, startatom, helium%beads) + xr
         END DO
      END IF

      snew = worm_path_action_worm_corrected(helium, helium%work, &
                                             startatom, startbead, endatom, endbead, &
                                             helium%worm_xtra_bead_work, helium%worm_atom_idx_work, helium%worm_bead_idx_work)

      IF (helium%solute_present) THEN
         snew = snew + worm_path_inter_action_head(pint_env, helium, helium%work, &
                                                   startatom, startbead, &
                                                   helium%worm_xtra_bead_work, helium%worm_atom_idx_work, helium%worm_bead_idx_work)
      END IF

      ! Metropolis:
      ! first compute ln of free density matrix
      distvec(:) = helium%pos(:, startatom, startbead) - helium%pos(:, endatom, endbead)
      CALL helium_pbc(helium, distvec)
      distsq = DOT_PRODUCT(distvec, distvec)
      ! action difference
      sdiff = sold - snew
      ! modify action difference due to extra bead
      sdiff = sdiff + distsq/(2.0_dp*helium%hb2m*REAL(l, dp)*helium%tau)
      sdiff = sdiff + 1.5_dp*LOG(REAL(l, dp)*helium%tau)
      sdiff = sdiff + helium%worm_ln_openclose_scale
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            ! transfer the new coordinates to work array
            IF (startbead < endbead) THEN
               ! everything belongs to the same atom
               DO kbead = startbead + 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the head bead
               DO kbead = startbead + 1, helium%beads
                  helium%work(:, startatom, kbead) = helium%pos(:, startatom, kbead)
               END DO
               ! transfer to atom containing the head bead
               DO kbead = 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            END IF
            helium%worm_is_closed = .TRUE.
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! restore original coordinates
            ! write back only changed atoms
            IF (startbead < endbead) THEN
               ! everything belongs to the same atom
               DO kbead = startbead + 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the head bead
               DO kbead = startbead + 1, helium%beads
                  helium%work(:, startatom, kbead) = helium%pos(:, startatom, kbead)
               END DO
               ! transfer to atom containing the head bead
               DO kbead = 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            END IF
            helium%worm_is_closed = .TRUE.
            !helium%worm_atom_idx_work = helium%worm_atom_idx
            !helium%worm_bead_idx_work = helium%worm_bead_idx
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      ! write changed coordinates to position array
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         DO kbead = startbead + 1, endbead - 1
            helium%pos(:, endatom, kbead) = helium%work(:, endatom, kbead)
         END DO
      ELSE
         ! is distributed among two atoms
         ! transfer to atom not containing the head bead
         DO kbead = startbead + 1, helium%beads
            helium%pos(:, startatom, kbead) = helium%work(:, startatom, kbead)
         END DO
         ! transfer to atom containing the head bead
         DO kbead = 1, endbead - 1
            helium%pos(:, endatom, kbead) = helium%work(:, endatom, kbead)
         END DO
      END IF
      helium%worm_xtra_bead(:) = helium%worm_xtra_bead_work(:)
      helium%worm_atom_idx = helium%worm_atom_idx_work
      helium%worm_bead_idx = helium%worm_bead_idx_work

      RETURN

   END SUBROUTINE worm_open_move

! **************************************************************************************************
!> \brief Close move back to diagonal elements of density matrix, permutation fixed in closed state
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_close_move(pint_env, helium, l, ac)

      TYPE(pint_env_type), INTENT(INOUT), POINTER        :: pint_env
      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      INTEGER, INTENT(IN)                                :: l
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: endatom, endbead, ia, ib, jbead, kbead, &
                                                            startatom, startbead
      LOGICAL                                            :: should_reject
      REAL(KIND=dp)                                      :: distsq, mass, rtmp, rtmpo, sdiff, snew, &
                                                            sold
      REAL(KIND=dp), DIMENSION(3)                        :: distvec, dr, dro, new_com, old_com
      REAL(KIND=dp), DIMENSION(3, l-1)                   :: newsection

      mass = helium%he_mass_au

      endatom = helium%worm_atom_idx
      endbead = helium%worm_bead_idx
      ! get index of the atom and bead, where the resampling of the head begins
      IF (l < endbead) THEN
         ! startbead belongs to the same atom
         startatom = endatom
         startbead = endbead - l
      ELSE
         ! startbead belongs to a different atom
         ! find previous atom (assuming l < nbeads)
         startatom = helium%iperm(endatom)
         startbead = endbead + helium%beads - l
      END IF
      sold = worm_path_action_worm_corrected(helium, helium%pos, &
                                             startatom, startbead, endatom, endbead, &
                                             helium%worm_xtra_bead, helium%worm_atom_idx, helium%worm_bead_idx)

      IF (helium%solute_present) THEN
         sold = sold + worm_path_inter_action_head(pint_env, helium, helium%pos, &
                                                   startatom, startbead, &
                                                   helium%worm_xtra_bead, helium%worm_atom_idx, helium%worm_bead_idx)
      END IF

      ! close between head and tail
      ! only l-1 beads need to be reconstructed
      CALL path_construct(helium, &
                          helium%pos(:, startatom, startbead), &
                          helium%pos(:, endatom, endbead), l - 1, &
                          newsection)

      ! transfer the new coordinates to work array
      jbead = 1
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         DO kbead = startbead + 1, endbead - 1
            helium%work(:, endatom, kbead) = newsection(:, jbead)
            jbead = jbead + 1
         END DO
      ELSE
         ! is distributed among two atoms
         ! transfer to atom not containing the head bead
         DO kbead = startbead + 1, helium%beads
            helium%work(:, startatom, kbead) = newsection(:, jbead)
            jbead = jbead + 1
         END DO
         ! transfer to atom containing the head bead
         DO kbead = 1, endbead - 1
            helium%work(:, endatom, kbead) = newsection(:, jbead)
            jbead = jbead + 1
         END DO
      END IF

      helium%worm_is_closed = .TRUE.

      snew = worm_path_action(helium, helium%work, &
                              startatom, startbead, endatom, endbead)

      IF (helium%solute_present) THEN
         ! yes this is correct, as the bead, that was split into tail and head only changes half
         ! therefore only half of its action needs to be considred
         ! and is cheated in here by passing it as head bead
         snew = snew + worm_path_inter_action_head(pint_env, helium, helium%work, &
                                                   startatom, startbead, &
                                                   helium%work(:, endatom, endbead), endatom, endbead)
      END IF

      ! Metropolis:
      ! first compute ln of free density matrix
      distvec(:) = helium%pos(:, startatom, startbead) - helium%pos(:, endatom, endbead)
      CALL helium_pbc(helium, distvec)
      distsq = DOT_PRODUCT(distvec, distvec)
      ! action difference
      sdiff = sold - snew
      ! modify action difference due to extra bead
      sdiff = sdiff - distsq/(2.0_dp*helium%hb2m*REAL(l, dp)*helium%tau)
      sdiff = sdiff - 1.5_dp*LOG(REAL(l, dp)*helium%tau)
      sdiff = sdiff - helium%worm_ln_openclose_scale
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            ! transfer the new coordinates to work array
            IF (startbead < endbead) THEN
               ! everything belongs to the same atom
               DO kbead = startbead + 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the head bead
               DO kbead = startbead + 1, helium%beads
                  helium%work(:, startatom, kbead) = helium%pos(:, startatom, kbead)
               END DO
               ! transfer to atom containing the head bead
               DO kbead = 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            END IF
            helium%worm_is_closed = .FALSE.
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! restore original coordinates
            ! write back only changed atoms
            IF (startbead < endbead) THEN
               ! everything belongs to the same atom
               DO kbead = startbead + 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the head bead
               DO kbead = startbead + 1, helium%beads
                  helium%work(:, startatom, kbead) = helium%pos(:, startatom, kbead)
               END DO
               ! transfer to atom containing the head bead
               DO kbead = 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            END IF
            helium%worm_is_closed = .FALSE.
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      ! write changed coordinates to position array
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         DO kbead = startbead + 1, endbead - 1
            helium%pos(:, endatom, kbead) = helium%work(:, endatom, kbead)
         END DO
      ELSE
         ! is distributed among two atoms
         ! transfer to atom not containing the head bead
         DO kbead = startbead + 1, helium%beads
            helium%pos(:, startatom, kbead) = helium%work(:, startatom, kbead)
         END DO
         ! transfer to atom containing the head bead
         DO kbead = 1, endbead - 1
            helium%pos(:, endatom, kbead) = helium%work(:, endatom, kbead)
         END DO
      END IF

      RETURN

   END SUBROUTINE worm_close_move

! **************************************************************************************************
!> \brief Move head in open state
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_head_move(pint_env, helium, l, ac)

      TYPE(pint_env_type), INTENT(INOUT), POINTER        :: pint_env
      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      INTEGER, INTENT(IN)                                :: l
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: endatom, endbead, ia, ib, idim, kbead, &
                                                            startatom, startbead
      LOGICAL                                            :: should_reject
      REAL(KIND=dp)                                      :: mass, rtmp, rtmpo, sdiff, snew, sold, xr
      REAL(KIND=dp), DIMENSION(3)                        :: dr, dro, new_com, old_com

      mass = helium%he_mass_au

      ! get index of the atom and bead, where the resampling of the head begins
      endatom = helium%worm_atom_idx
      endbead = helium%worm_bead_idx
      IF (l < endbead) THEN
         ! startbead belongs to the same atom
         startatom = endatom
         startbead = endbead - l
      ELSE
         ! startbead belongs to a different atom
         ! find previous atom (assuming l < nbeads)
         startatom = helium%iperm(endatom)
         startbead = endbead + helium%beads - l
      END IF

      sold = worm_path_action_worm_corrected(helium, helium%pos, &
                                             startatom, startbead, endatom, endbead, &
                                             helium%worm_xtra_bead, helium%worm_atom_idx, helium%worm_bead_idx)

      IF (helium%solute_present) THEN
         sold = sold + worm_path_inter_action_head(pint_env, helium, helium%pos, &
                                                   startatom, startbead, &
                                                   helium%worm_xtra_bead, helium%worm_atom_idx, helium%worm_bead_idx)
      END IF

      ! alternative grow with consecutive gaussians
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         ! gro head from startbead
         DO kbead = startbead + 1, endbead - 1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, startatom, kbead) = helium%work(idim, startatom, kbead - 1) + xr
            END DO
         END DO
         ! last grow head bead
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, startatom, endbead - 1) + xr
         END DO
      ELSE IF (endbead /= 1) THEN
         ! is distributed among two atoms
         ! grow from startbead
         DO kbead = startbead + 1, helium%beads
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, startatom, kbead) = helium%work(idim, startatom, kbead - 1) + xr
            END DO
         END DO
         ! bead one of endatom relative to last on startatom
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%work(idim, endatom, 1) = helium%work(idim, startatom, helium%beads) + xr
         END DO
         ! everything on endatom
         DO kbead = 2, endbead - 1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, endatom, kbead) = helium%work(idim, endatom, kbead - 1) + xr
            END DO
         END DO
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, endatom, endbead - 1) + xr
         END DO
      ELSE ! imagtimewrap and headbead = 1
         ! is distributed among two atoms
         ! grow from startbead
         DO kbead = startbead + 1, helium%beads
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, startatom, kbead) = helium%work(idim, startatom, kbead - 1) + xr
            END DO
         END DO
         ! bead one of endatom relative to last on startatom
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, startatom, helium%beads) + xr
         END DO
      END IF

      snew = worm_path_action_worm_corrected(helium, helium%work, &
                                             startatom, startbead, endatom, endbead, &
                                             helium%worm_xtra_bead_work, helium%worm_atom_idx_work, helium%worm_bead_idx_work)

      IF (helium%solute_present) THEN
         snew = snew + worm_path_inter_action_head(pint_env, helium, helium%work, &
                                                   startatom, startbead, &
                                                   helium%worm_xtra_bead_work, helium%worm_atom_idx_work, helium%worm_bead_idx_work)
      END IF

      ! Metropolis:
      ! action difference
      sdiff = sold - snew
      ! modify action difference due to extra bead
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            ! transfer the new coordinates to work array
            IF (startbead < endbead) THEN
               ! everything belongs to the same atom
               DO kbead = startbead + 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the head bead
               DO kbead = startbead + 1, helium%beads
                  helium%work(:, startatom, kbead) = helium%pos(:, startatom, kbead)
               END DO
               ! transfer to atom containing the head bead
               DO kbead = 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            END IF
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! restore original coordinates
            ! write back only changed atoms
            IF (startbead < endbead) THEN
               ! everything belongs to the same atom
               DO kbead = startbead + 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the head bead
               DO kbead = startbead + 1, helium%beads
                  helium%work(:, startatom, kbead) = helium%pos(:, startatom, kbead)
               END DO
               ! transfer to atom containing the head bead
               DO kbead = 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            END IF
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      ! write changed coordinates to position array
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         DO kbead = startbead + 1, endbead - 1
            helium%pos(:, endatom, kbead) = helium%work(:, endatom, kbead)
         END DO
      ELSE
         ! is distributed among two atoms
         ! transfer to atom not containing the head bead
         DO kbead = startbead + 1, helium%beads
            helium%pos(:, startatom, kbead) = helium%work(:, startatom, kbead)
         END DO
         ! transfer to atom containing the head bead
         DO kbead = 1, endbead - 1
            helium%pos(:, endatom, kbead) = helium%work(:, endatom, kbead)
         END DO
      END IF
      helium%worm_xtra_bead(:) = helium%worm_xtra_bead_work(:)

      RETURN

   END SUBROUTINE worm_head_move

! **************************************************************************************************
!> \brief Move tail in open state
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_tail_move(pint_env, helium, l, ac)

      TYPE(pint_env_type), INTENT(INOUT), POINTER        :: pint_env
      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      INTEGER, INTENT(IN)                                :: l
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: endatom, endbead, ia, ib, idim, kbead, &
                                                            startatom, startbead
      LOGICAL                                            :: should_reject
      REAL(KIND=dp)                                      :: mass, rtmp, rtmpo, sdiff, snew, sold, xr
      REAL(KIND=dp), DIMENSION(3)                        :: dr, dro, new_com, old_com

      mass = helium%he_mass_au

      ! get index of the atom and bead, where the resampling of the tail ends
      startatom = helium%worm_atom_idx
      startbead = helium%worm_bead_idx
      endbead = startbead + l

      IF (endbead <= helium%beads) THEN
         ! endbead belongs to the same atom
         endatom = startatom
      ELSE
         ! endbead belongs to a different atom
         ! find next atom (assuming l < nbeads)
         endatom = helium%permutation(startatom)
         endbead = endbead - helium%beads
      END IF

      !yes this is correct, as the head does not play any role here
      sold = worm_path_action(helium, helium%pos, &
                              startatom, startbead, endatom, endbead)

      IF (helium%solute_present) THEN
         sold = sold + worm_path_inter_action_tail(pint_env, helium, helium%pos, &
                                                   endatom, endbead, &
                                                   helium%worm_atom_idx, helium%worm_bead_idx)
      END IF

      ! alternative grow with consecutive gaussians
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         ! gro tail from endbead to startbead (confusing eh?)
         DO kbead = endbead - 1, startbead, -1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, startatom, kbead) = helium%work(idim, startatom, kbead + 1) + xr
            END DO
         END DO
      ELSE
         ! is distributed among two atoms
         ! grow from endbead
         DO kbead = endbead - 1, 1, -1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, endatom, kbead) = helium%work(idim, endatom, kbead + 1) + xr
            END DO
         END DO

         ! over imaginary time boundary
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%work(idim, startatom, helium%beads) = helium%work(idim, endatom, 1) + xr
         END DO

         ! rest on startatom
         DO kbead = helium%beads - 1, startbead, -1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, startatom, kbead) = helium%work(idim, startatom, kbead + 1) + xr
            END DO
         END DO
      END IF

      !yes this is correct, as the head does not play any role here
      snew = worm_path_action(helium, helium%work, &
                              startatom, startbead, endatom, endbead)

      IF (helium%solute_present) THEN
         snew = snew + worm_path_inter_action_tail(pint_env, helium, helium%work, &
                                                   endatom, endbead, &
                                                   helium%worm_atom_idx_work, helium%worm_bead_idx_work)
      END IF

      ! Metropolis:
      ! action difference
      sdiff = sold - snew
      ! modify action difference due to extra bead
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            ! transfer the new coordinates to work array
            IF (startbead < endbead) THEN
               ! everything belongs to the same atom
               DO kbead = startbead, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the tail bead
               DO kbead = startbead, helium%beads
                  helium%work(:, startatom, kbead) = helium%pos(:, startatom, kbead)
               END DO
               ! transfer to atom containing the tail bead
               DO kbead = 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            END IF
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! restore original coordinates
            ! write back only changed atoms
            IF (startbead < endbead) THEN
               ! everything belongs to the same atom
               DO kbead = startbead, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the tail bead
               DO kbead = startbead, helium%beads
                  helium%work(:, startatom, kbead) = helium%pos(:, startatom, kbead)
               END DO
               ! transfer to atom containing the tail bead
               DO kbead = 1, endbead - 1
                  helium%work(:, endatom, kbead) = helium%pos(:, endatom, kbead)
               END DO
            END IF
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      ! write changed coordinates to position array
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         DO kbead = startbead, endbead - 1
            helium%pos(:, endatom, kbead) = helium%work(:, endatom, kbead)
         END DO
      ELSE
         ! is distributed among two atoms
         ! transfer to atom not containing the tail bead
         DO kbead = startbead, helium%beads
            helium%pos(:, startatom, kbead) = helium%work(:, startatom, kbead)
         END DO
         ! transfer to atom containing the tail bead
         DO kbead = 1, endbead - 1
            helium%pos(:, endatom, kbead) = helium%work(:, endatom, kbead)
         END DO
      END IF

      RETURN

   END SUBROUTINE worm_tail_move

! **************************************************************************************************
!> \brief Crawl forward in open state, can avoid being trapped in open state
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_crawl_move_forward(pint_env, helium, l, ac)

      TYPE(pint_env_type), INTENT(IN), POINTER           :: pint_env
      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      INTEGER, INTENT(IN)                                :: l
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: ia, ib, idim, kbead
      LOGICAL                                            :: should_reject
      REAL(KIND=dp)                                      :: mass, newtailpot, oldheadpot, &
                                                            oldtailpot, rtmp, rtmpo, sdiff, snew, &
                                                            sold, xr
      REAL(KIND=dp), DIMENSION(3)                        :: dr, dro, new_com, old_com

      mass = helium%he_mass_au

      ! determine position of new head in imaginary time
      helium%worm_bead_idx_work = helium%worm_bead_idx + l
      IF (helium%worm_bead_idx_work > helium%beads) THEN
         helium%worm_bead_idx_work = helium%worm_bead_idx_work - helium%beads
         helium%worm_atom_idx_work = helium%permutation(helium%worm_atom_idx)
      ELSE
         helium%worm_atom_idx_work = helium%worm_atom_idx
      END IF

      ! compute action before move
      ! head is not involved in pair action
      sold = worm_path_action(helium, helium%pos, &
                              helium%worm_atom_idx, helium%worm_bead_idx, &
                              helium%worm_atom_idx_work, helium%worm_bead_idx_work)
      IF (helium%solute_present) THEN
         !this will leave out the old and new tail bead
         ! due to efficiency reasons they are treated separately
         sold = sold + worm_path_inter_action(pint_env, helium, helium%pos, &
                                              helium%worm_atom_idx, helium%worm_bead_idx, &
                                              helium%worm_atom_idx_work, helium%worm_bead_idx_work)

         ! compute old/new head/tail interactions
         ! old tail
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     helium%worm_atom_idx, helium%worm_bead_idx, &
                                     helium%pos(:, helium%worm_atom_idx, helium%worm_bead_idx), &
                                     energy=oldtailpot)
         oldtailpot = oldtailpot*helium%tau

         ! new tail
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     helium%worm_atom_idx_work, helium%worm_bead_idx_work, &
                                     helium%pos(:, helium%worm_atom_idx_work, helium%worm_bead_idx_work), &
                                     energy=newtailpot)
         newtailpot = newtailpot*helium%tau

         ! old head
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     helium%worm_atom_idx, helium%worm_bead_idx, &
                                     helium%worm_xtra_bead, &
                                     energy=oldheadpot)
         oldheadpot = oldheadpot*helium%tau

         ! new head is not known yet

         sold = sold + 0.5_dp*(oldtailpot + oldheadpot) + newtailpot
      END IF

      ! copy over old head position to working array and grow from there
      helium%work(:, helium%worm_atom_idx, helium%worm_bead_idx) = helium%worm_xtra_bead

      ! grow head part with consecutive gaussians
      IF (helium%worm_bead_idx < helium%worm_bead_idx_work) THEN
         ! everything belongs to the same atom
         ! gro head from startbead
         DO kbead = helium%worm_bead_idx + 1, helium%worm_bead_idx_work - 1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, helium%worm_atom_idx, kbead) = helium%work(idim, helium%worm_atom_idx, kbead - 1) + xr
            END DO
         END DO
         ! last grow head bead
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, helium%worm_atom_idx, helium%worm_bead_idx_work - 1) + xr
         END DO
      ELSE IF (helium%worm_bead_idx_work /= 1) THEN
         ! is distributed among two atoms
         ! grow from startbead
         DO kbead = helium%worm_bead_idx + 1, helium%beads
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, helium%worm_atom_idx, kbead) = helium%work(idim, helium%worm_atom_idx, kbead - 1) + xr
            END DO
         END DO
         ! bead one of endatom relative to last on helium%worm_atom_idx
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%work(idim, helium%worm_atom_idx_work, 1) = helium%work(idim, helium%worm_atom_idx, helium%beads) + xr
         END DO
         ! everything on endatom
         DO kbead = 2, helium%worm_bead_idx_work - 1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, helium%worm_atom_idx_work, kbead) = helium%work(idim, helium%worm_atom_idx_work, kbead - 1) + xr
            END DO
         END DO
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, helium%worm_atom_idx_work, helium%worm_bead_idx_work - 1) + xr
         END DO
      ELSE ! imagtimewrap and headbead = 1
         ! is distributed among two atoms
         ! grow from startbead
         DO kbead = helium%worm_bead_idx + 1, helium%beads
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, helium%worm_atom_idx, kbead) = helium%work(idim, helium%worm_atom_idx, kbead - 1) + xr
            END DO
         END DO
         ! bead one of endatom relative to last on helium%worm_atom_idx
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%worm_xtra_bead_work(idim) = helium%work(idim, helium%worm_atom_idx, helium%beads) + xr
         END DO
      END IF

      snew = worm_path_action_worm_corrected(helium, helium%work, &
                                             helium%worm_atom_idx, helium%worm_bead_idx, &
                                             helium%worm_atom_idx_work, helium%worm_bead_idx_work, &
                                             helium%worm_xtra_bead_work, helium%worm_atom_idx_work, helium%worm_bead_idx_work)

      IF (helium%solute_present) THEN
         snew = snew + worm_path_inter_action_head(pint_env, helium, helium%work, &
                                                   helium%worm_atom_idx, helium%worm_bead_idx, &
                                                   helium%worm_xtra_bead_work, helium%worm_atom_idx_work, helium%worm_bead_idx_work)
         snew = snew + 0.5_dp*newtailpot + oldheadpot
      END IF

      ! Metropolis:
      ! action difference
      sdiff = sold - snew
      ! modify action difference due to extra bead
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            ! transfer the new coordinates to work array
            IF (helium%worm_bead_idx < helium%worm_bead_idx_work) THEN
               ! everything belongs to the same atom
               DO kbead = helium%worm_bead_idx, helium%worm_bead_idx_work - 1
                  helium%work(:, helium%worm_atom_idx_work, kbead) = helium%pos(:, helium%worm_atom_idx_work, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the head bead
               DO kbead = helium%worm_bead_idx, helium%beads
                  helium%work(:, helium%worm_atom_idx, kbead) = helium%pos(:, helium%worm_atom_idx, kbead)
               END DO
               ! transfer to atom containing the head bead
               DO kbead = 1, helium%worm_bead_idx_work - 1
                  helium%work(:, helium%worm_atom_idx_work, kbead) = helium%pos(:, helium%worm_atom_idx_work, kbead)
               END DO
            END IF
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            helium%worm_bead_idx_work = helium%worm_bead_idx
            helium%worm_atom_idx_work = helium%worm_atom_idx
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! restore original coordinates
            ! write back only changed atoms
            IF (helium%worm_bead_idx < helium%worm_bead_idx_work) THEN
               ! everything belongs to the same atom
               DO kbead = helium%worm_bead_idx, helium%worm_bead_idx_work - 1
                  helium%work(:, helium%worm_atom_idx_work, kbead) = helium%pos(:, helium%worm_atom_idx_work, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the head bead
               DO kbead = helium%worm_bead_idx, helium%beads
                  helium%work(:, helium%worm_atom_idx, kbead) = helium%pos(:, helium%worm_atom_idx, kbead)
               END DO
               ! transfer to atom containing the head bead
               DO kbead = 1, helium%worm_bead_idx_work - 1
                  helium%work(:, helium%worm_atom_idx_work, kbead) = helium%pos(:, helium%worm_atom_idx_work, kbead)
               END DO
            END IF
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            helium%worm_bead_idx_work = helium%worm_bead_idx
            helium%worm_atom_idx_work = helium%worm_atom_idx
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      ! write changed coordinates to position array
      IF (helium%worm_bead_idx < helium%worm_bead_idx_work) THEN
         ! everything belongs to the same atom
         DO kbead = helium%worm_bead_idx, helium%worm_bead_idx_work - 1
            helium%pos(:, helium%worm_atom_idx_work, kbead) = helium%work(:, helium%worm_atom_idx_work, kbead)
         END DO
      ELSE
         ! is distributed among two atoms
         ! transfer to atom not containing the head bead
         DO kbead = helium%worm_bead_idx, helium%beads
            helium%pos(:, helium%worm_atom_idx, kbead) = helium%work(:, helium%worm_atom_idx, kbead)
         END DO
         ! transfer to atom containing the head bead
         DO kbead = 1, helium%worm_bead_idx_work - 1
            helium%pos(:, helium%worm_atom_idx_work, kbead) = helium%work(:, helium%worm_atom_idx_work, kbead)
         END DO
      END IF
      helium%worm_xtra_bead(:) = helium%worm_xtra_bead_work(:)
      helium%worm_bead_idx = helium%worm_bead_idx_work
      helium%worm_atom_idx = helium%worm_atom_idx_work

      RETURN

   END SUBROUTINE worm_crawl_move_forward

! **************************************************************************************************
!> \brief Crawl backwards in open state, can avoid being trapped in open state
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_crawl_move_backward(pint_env, helium, l, ac)

      TYPE(pint_env_type), INTENT(IN), POINTER           :: pint_env
      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      INTEGER, INTENT(IN)                                :: l
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: ia, ib, idim, kbead
      LOGICAL                                            :: should_reject
      REAL(KIND=dp)                                      :: mass, newheadpot, oldheadpot, &
                                                            oldtailpot, rtmp, rtmpo, sdiff, snew, &
                                                            sold, xr
      REAL(KIND=dp), DIMENSION(3)                        :: dr, dro, new_com, old_com

      mass = helium%he_mass_au

      ! determine position of new head in imaginary time
      helium%worm_bead_idx_work = helium%worm_bead_idx - l
      IF (helium%worm_bead_idx_work < 1) THEN
         helium%worm_bead_idx_work = helium%worm_bead_idx_work + helium%beads
         helium%worm_atom_idx_work = helium%iperm(helium%worm_atom_idx)
      ELSE
         helium%worm_atom_idx_work = helium%worm_atom_idx
      END IF

      ! compute action before move
      ! head is not involved in pair action
      sold = worm_path_action_worm_corrected(helium, helium%work, &
                                             helium%worm_atom_idx_work, helium%worm_bead_idx_work, &
                                             helium%worm_atom_idx, helium%worm_bead_idx, &
                                             helium%worm_xtra_bead, helium%worm_atom_idx, helium%worm_bead_idx)

      IF (helium%solute_present) THEN
         !this will leave out the old and new tail bead
         ! due to efficiency reasons they are treated separately
         sold = sold + worm_path_inter_action(pint_env, helium, helium%pos, &
                                              helium%worm_atom_idx_work, helium%worm_bead_idx_work, &
                                              helium%worm_atom_idx, helium%worm_bead_idx)

         ! compute old/new head/tail interactions
         ! old tail
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     helium%worm_atom_idx, helium%worm_bead_idx, &
                                     helium%pos(:, helium%worm_atom_idx, helium%worm_bead_idx), &
                                     energy=oldtailpot)
         oldtailpot = oldtailpot*helium%tau

         ! old head
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     helium%worm_atom_idx, helium%worm_bead_idx, &
                                     helium%worm_xtra_bead, &
                                     energy=oldheadpot)
         oldheadpot = oldheadpot*helium%tau

         ! new head
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     helium%worm_atom_idx_work, helium%worm_bead_idx_work, &
                                     helium%pos(:, helium%worm_atom_idx_work, helium%worm_bead_idx_work), &
                                     energy=newheadpot)
         newheadpot = newheadpot*helium%tau

         ! new tail not known yet

         sold = sold + 0.5_dp*(oldtailpot + oldheadpot) + newheadpot
      END IF

      ! copy position to the head bead
      helium%worm_xtra_bead_work = helium%pos(:, helium%worm_atom_idx_work, helium%worm_bead_idx_work)

      ! alternative grow with consecutive gaussians
      IF (helium%worm_bead_idx_work < helium%worm_bead_idx) THEN
         ! everything belongs to the same atom
         ! gro tail from endbead to startbead (confusing eh?)
         DO kbead = helium%worm_bead_idx - 1, helium%worm_bead_idx_work, -1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, helium%worm_atom_idx, kbead) = helium%work(idim, helium%worm_atom_idx, kbead + 1) + xr
            END DO

         END DO
      ELSE
         ! is distributed among two atoms
         ! grow from endbead
         DO kbead = helium%worm_bead_idx - 1, 1, -1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, helium%worm_atom_idx, kbead) = helium%work(idim, helium%worm_atom_idx, kbead + 1) + xr
            END DO
         END DO

         ! over imaginary time boundary
         DO idim = 1, 3
            xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
            helium%work(idim, helium%worm_atom_idx_work, helium%beads) = helium%work(idim, helium%worm_atom_idx, 1) + xr
         END DO

         ! rest on startatom
         DO kbead = helium%beads - 1, helium%worm_bead_idx_work, -1
            DO idim = 1, 3
               xr = helium%rng_stream_gaussian%next(variance=helium%hb2m*helium%tau)
               helium%work(idim, helium%worm_atom_idx_work, kbead) = helium%work(idim, helium%worm_atom_idx_work, kbead + 1) + xr
            END DO
         END DO
      END IF

      snew = worm_path_action(helium, helium%work, &
                              helium%worm_atom_idx_work, helium%worm_bead_idx_work, &
                              helium%worm_atom_idx, helium%worm_bead_idx)

      IF (helium%solute_present) THEN
         snew = snew + worm_path_inter_action_tail(pint_env, helium, helium%work, &
                                                   helium%worm_atom_idx, helium%worm_bead_idx, &
                                                   helium%worm_atom_idx_work, helium%worm_bead_idx_work)
         snew = snew + 0.5_dp*newheadpot + oldtailpot
      END IF

      ! Metropolis:
      ! action difference
      sdiff = sold - snew
      ! modify action difference due to extra bead
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            ! transfer the new coordinates to work array
            IF (helium%worm_bead_idx_work < helium%worm_bead_idx) THEN
               ! everything belongs to the same atom
               DO kbead = helium%worm_bead_idx_work, helium%worm_bead_idx - 1
                  helium%work(:, helium%worm_atom_idx, kbead) = helium%pos(:, helium%worm_atom_idx, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the tail bead
               DO kbead = helium%worm_bead_idx_work, helium%beads
                  helium%work(:, helium%worm_atom_idx_work, kbead) = helium%pos(:, helium%worm_atom_idx_work, kbead)
               END DO
               ! transfer to atom containing the tail bead
               DO kbead = 1, helium%worm_bead_idx - 1
                  helium%work(:, helium%worm_atom_idx, kbead) = helium%pos(:, helium%worm_atom_idx, kbead)
               END DO
            END IF
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            helium%worm_bead_idx_work = helium%worm_bead_idx
            helium%worm_atom_idx_work = helium%worm_atom_idx
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! restore original coordinates
            ! write back only changed atoms
            ! write back only changed atoms
            ! transfer the new coordinates to work array
            IF (helium%worm_bead_idx_work < helium%worm_bead_idx) THEN
               ! everything belongs to the same atom
               DO kbead = helium%worm_bead_idx_work, helium%worm_bead_idx - 1
                  helium%work(:, helium%worm_atom_idx, kbead) = helium%pos(:, helium%worm_atom_idx, kbead)
               END DO
            ELSE
               ! is distributed among two atoms
               ! transfer to atom not containing the tail bead
               DO kbead = helium%worm_bead_idx_work, helium%beads
                  helium%work(:, helium%worm_atom_idx_work, kbead) = helium%pos(:, helium%worm_atom_idx_work, kbead)
               END DO
               ! transfer to atom containing the tail bead
               DO kbead = 1, helium%worm_bead_idx - 1
                  helium%work(:, helium%worm_atom_idx, kbead) = helium%pos(:, helium%worm_atom_idx, kbead)
               END DO
            END IF
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            helium%worm_bead_idx_work = helium%worm_bead_idx
            helium%worm_atom_idx_work = helium%worm_atom_idx
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      ! write changed coordinates to position array
      ! write back only changed atoms
      IF (helium%worm_bead_idx_work < helium%worm_bead_idx) THEN
         ! everything belongs to the same atom
         DO kbead = helium%worm_bead_idx_work, helium%worm_bead_idx - 1
            helium%pos(:, helium%worm_atom_idx, kbead) = helium%work(:, helium%worm_atom_idx, kbead)
         END DO
      ELSE
         ! is distributed among two atoms
         ! transfer to atom not containing the tail bead
         DO kbead = helium%worm_bead_idx_work, helium%beads
            helium%pos(:, helium%worm_atom_idx_work, kbead) = helium%work(:, helium%worm_atom_idx_work, kbead)
         END DO
         ! transfer to atom containing the tail bead
         DO kbead = 1, helium%worm_bead_idx - 1
            helium%pos(:, helium%worm_atom_idx, kbead) = helium%work(:, helium%worm_atom_idx, kbead)
         END DO
      END IF
      helium%worm_xtra_bead(:) = helium%worm_xtra_bead_work(:)
      helium%worm_bead_idx = helium%worm_bead_idx_work
      helium%worm_atom_idx = helium%worm_atom_idx_work

      RETURN

   END SUBROUTINE worm_crawl_move_backward

! **************************************************************************************************
!> \brief Free particle density matrix
!> \param helium ...
!> \param startpos ...
!> \param endpos ...
!> \param l ...
!> \return ...
!> \author fuhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION free_density_matrix(helium, startpos, endpos, l) RESULT(rho)

      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: startpos, endpos
      INTEGER, INTENT(IN)                                :: l

      REAL(KIND=dp)                                      :: distsq, prefac
      REAL(KIND=dp), DIMENSION(3)                        :: dvec

      dvec(:) = startpos(:) - endpos(:)
      CALL helium_pbc(helium, dvec)
      distsq = DOT_PRODUCT(dvec, dvec)

      prefac = 0.5_dp/(helium%hb2m*REAL(l, dp)*helium%tau)
      rho = -prefac*distsq
      rho = EXP(rho)
      rho = rho*(SQRT(prefac/pi))**3

      RETURN

   END FUNCTION free_density_matrix

! **************************************************************************************************
!> \brief Swap move in open state to change permutation state
!> \param pint_env ...
!> \param helium ...
!> \param natoms ...
!> \param l ...
!> \param ac ...
!> \author fuhl
! **************************************************************************************************
   SUBROUTINE worm_swap_move(pint_env, helium, natoms, l, ac)

      TYPE(pint_env_type), INTENT(INOUT), POINTER        :: pint_env
      TYPE(helium_solvent_type), INTENT(INOUT), POINTER  :: helium
      INTEGER, INTENT(IN)                                :: natoms, l
      INTEGER, INTENT(OUT)                               :: ac

      INTEGER                                            :: bendatom, bstartatom, endbead, &
                                                            excludeatom, fendatom, fstartatom, ia, &
                                                            iatom, ib, ibead, jbead, kbead, &
                                                            startbead, tmpi
      LOGICAL                                            :: should_reject
      REAL(KIND=dp)                                      :: backwarddensmatsum, forwarddensmatsum, &
                                                            mass, newheadpotential, &
                                                            oldheadpotential, rtmp, rtmpo, sdiff, &
                                                            snew, sold
      REAL(KIND=dp), DIMENSION(3)                        :: dr, dro, new_com, old_com
      REAL(KIND=dp), DIMENSION(3, l)                     :: newsection
      REAL(KIND=dp), DIMENSION(natoms)                   :: backwarddensmat, forwarddensmat

      ! first the endbead of the reconstruction is needed
      startbead = helium%worm_bead_idx
      endbead = helium%worm_bead_idx + l + 1
      fstartatom = helium%worm_atom_idx
      excludeatom = fstartatom
      ! compute the atomwise probabilities to be the worms swap partner
      ! Check if the imaginary time section belongs to two atoms
      IF (endbead > helium%beads) THEN
         endbead = endbead - helium%beads
         ! exclude atom is the one not to connect to because it will result in an unnatural state
         excludeatom = helium%permutation(excludeatom)
      END IF
      DO iatom = 1, excludeatom - 1
         forwarddensmat(iatom) = free_density_matrix(helium, helium%worm_xtra_bead(:), &
                                                     helium%pos(:, iatom, endbead), l + 1)
      END DO
      forwarddensmat(excludeatom) = 0.0_dp
      DO iatom = excludeatom + 1, natoms
         forwarddensmat(iatom) = free_density_matrix(helium, helium%worm_xtra_bead(:), &
                                                     helium%pos(:, iatom, endbead), l + 1)
      END DO
      forwarddensmatsum = SUM(forwarddensmat)

      ! Select an atom with its corresponding probability
      rtmp = helium%rng_stream_uniform%next()*forwarddensmatsum
      fendatom = 1
      DO WHILE (rtmp >= forwarddensmat(fendatom))
         rtmp = rtmp - forwarddensmat(fendatom)
         fendatom = fendatom + 1
      END DO
      ! just for numerical safety
      fendatom = MIN(fendatom, natoms)

      ! ensure detailed balance
      IF (endbead > startbead) THEN
         bstartatom = fendatom
      ELSE
         bstartatom = helium%iperm(fendatom)
      END IF
      bendatom = fendatom

      DO iatom = 1, excludeatom - 1
         backwarddensmat(iatom) = free_density_matrix(helium, &
                                                      helium%pos(:, bstartatom, startbead), &
                                                      helium%pos(:, iatom, endbead), l + 1)
      END DO
      backwarddensmat(excludeatom) = 0.0_dp
      DO iatom = excludeatom + 1, natoms
         backwarddensmat(iatom) = free_density_matrix(helium, &
                                                      helium%pos(:, bstartatom, startbead), &
                                                      helium%pos(:, iatom, endbead), l + 1)
      END DO

      backwarddensmatsum = SUM(backwarddensmat)

      mass = helium%he_mass_au

      !compute action before the move
      sold = worm_path_action(helium, helium%pos, &
                              bstartatom, startbead, fendatom, endbead)

      IF (helium%solute_present) THEN
         ! no special head treatment needed, as it will change due to swapping
         ! the worm gap and due to primitive coupling no cross bead terms are considered
         sold = sold + worm_path_inter_action(pint_env, helium, helium%pos, &
                                              bstartatom, startbead, fendatom, endbead)
         ! compute potential of old and new head here (only once, as later only a rescaling is necessary)
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     fstartatom, startbead, helium%worm_xtra_bead, &
                                     energy=oldheadpotential)
         oldheadpotential = oldheadpotential*helium%tau

         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     bstartatom, startbead, helium%pos(:, bstartatom, startbead), &
                                     energy=newheadpotential)
         newheadpotential = newheadpotential*helium%tau

         sold = sold + 0.5_dp*oldheadpotential + newheadpotential
      END IF

      ! construct a new path connecting the start and endbead
      ! need to be the old coordinates due to reordering of the work coordinates
      CALL path_construct(helium, &
                          helium%worm_xtra_bead(:), &
                          helium%pos(:, fendatom, endbead), l, &
                          newsection)

      !write new path segment to work array
      !first the part that is guaranteed to fit on the coorinates of startatom of the
      jbead = 1
      IF (startbead < endbead) THEN
         ! everything belongs to the same atom
         DO kbead = startbead + 1, endbead - 1
            helium%work(:, fstartatom, kbead) = newsection(:, jbead)
            jbead = jbead + 1
         END DO
      ELSE
         ! is distributed among two atoms
         ! transfer to atom not containing the head bead
         DO kbead = startbead + 1, helium%beads
            helium%work(:, fstartatom, kbead) = newsection(:, jbead)
            jbead = jbead + 1
         END DO
         ! rest to the second atom
         DO ibead = 1, endbead - 1
            helium%work(:, fendatom, ibead) = newsection(:, jbead)
            jbead = jbead + 1
         END DO
      END IF

      !exchange coordinates head coordinate first
      helium%work(:, helium%worm_atom_idx, helium%worm_bead_idx) = helium%worm_xtra_bead(:)
      helium%worm_xtra_bead_work(:) = helium%pos(:, bstartatom, startbead)

      ! move coordinates from former worm atom tail to new worm atom tail
      DO ib = startbead, helium%beads
         helium%work(:, bstartatom, ib) = helium%pos(:, fstartatom, ib)
      END DO
      ! need to copy the rest of pselatom to former worm atom
      IF (endbead > startbead) THEN
         DO ib = endbead, helium%beads
            helium%work(:, fstartatom, ib) = helium%pos(:, bstartatom, ib)
         END DO
      END IF

      !update permtable
      tmpi = helium%permutation(bstartatom)
      helium%permutation(bstartatom) = helium%permutation(fstartatom)
      helium%permutation(fstartatom) = tmpi
      !update inverse permtable
      DO ib = 1, SIZE(helium%permutation)
         helium%iperm(helium%permutation(ib)) = ib
      END DO
      helium%worm_atom_idx_work = bstartatom
      ! action after the move

      IF (endbead > startbead) THEN
         snew = worm_path_action(helium, helium%work, &
                                 fstartatom, startbead, fstartatom, endbead)
         IF (helium%solute_present) THEN
            ! no special head treatment needed, because a swap can't go over
            ! the worm gap and due to primitive coupling no cross bead terms are considered
            snew = snew + worm_path_inter_action(pint_env, helium, helium%work, &
                                                 fstartatom, startbead, fstartatom, endbead)

            ! add the previously computed old and new head actions
            snew = snew + oldheadpotential + 0.5_dp*newheadpotential
         END IF
      ELSE
         snew = worm_path_action(helium, helium%work, &
                                 fstartatom, startbead, helium%permutation(fstartatom), endbead)
         IF (helium%solute_present) THEN
            ! no special head treatment needed, because a swap can't go over
            ! the worm gap and due to primitive coupling no cross bead terms are considered
            snew = snew + worm_path_inter_action(pint_env, helium, helium%work, &
                                                 fstartatom, startbead, helium%permutation(fstartatom), endbead)

            ! add the previously computed old and new head actions
            snew = snew + oldheadpotential + 0.5_dp*newheadpotential
         END IF
      END IF

      ! Metropolis:
      sdiff = sold - snew
      sdiff = sdiff + LOG(forwarddensmatsum/backwarddensmatsum)
      IF (sdiff < 0) THEN
         should_reject = .FALSE.
         IF (sdiff < -100.0_dp) THEN ! To protect from exponential underflow
            should_reject = .TRUE.
         ELSE
            rtmp = helium%rng_stream_uniform%next()
            IF (EXP(sdiff) < rtmp) THEN
               should_reject = .TRUE.
            END IF
         END IF
         IF (should_reject) THEN
            ! rejected !
            ! write back only changed atoms
            DO kbead = startbead, helium%beads
               helium%work(:, bstartatom, kbead) = helium%pos(:, bstartatom, kbead)
            END DO
            DO kbead = startbead, helium%beads
               helium%work(:, fstartatom, kbead) = helium%pos(:, fstartatom, kbead)
            END DO
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            helium%worm_atom_idx_work = helium%worm_atom_idx
            IF (startbead > endbead) THEN
               DO kbead = 1, endbead
                  helium%work(:, fendatom, kbead) = helium%pos(:, fendatom, kbead)
               END DO
            END IF
            ! reset permtable
            tmpi = helium%permutation(bstartatom)
            helium%permutation(bstartatom) = helium%permutation(fstartatom)
            helium%permutation(fstartatom) = tmpi
            !update inverse permtable
            DO ib = 1, SIZE(helium%permutation)
               helium%iperm(helium%permutation(ib)) = ib
            END DO
            ac = 0
            RETURN
         END IF
      END IF

      ! for now accepted
      ! rejection of the whole move if any centroid is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      ! AND ist not moving towards the center
      IF (.NOT. helium%periodic) THEN
         IF (helium%solute_present) THEN
            new_com(:) = helium%center(:)
            old_com(:) = new_com(:)
         ELSE
            new_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  new_com(:) = new_com(:) + helium%work(:, ia, ib)
               END DO
            END DO
            new_com(:) = new_com(:)/(REAL(helium%atoms*helium%beads, dp))
            ! also compute the old center of mass (optimization potential)
            old_com(:) = 0.0_dp
            DO ia = 1, helium%atoms
               DO ib = 1, helium%beads
                  old_com(:) = old_com(:) + helium%pos(:, ia, ib)
               END DO
            END DO
            old_com(:) = old_com(:)/(REAL(helium%atoms*helium%beads, dp))
         END IF
         should_reject = .FALSE.
         atom: DO ia = 1, helium%atoms
            dr(:) = 0.0_dp
            DO ib = 1, helium%beads
               dr(:) = dr(:) + helium%work(:, ia, ib) - new_com(:)
            END DO
            dr(:) = dr(:)/REAL(helium%beads, dp)
            rtmp = DOT_PRODUCT(dr, dr)
            IF (rtmp >= helium%droplet_radius**2) THEN
               dro(:) = 0.0_dp
               DO ib = 1, helium%beads
                  dro(:) = dro(:) + helium%pos(:, ia, ib) - old_com(:)
               END DO
               dro(:) = dro(:)/REAL(helium%beads, dp)
               rtmpo = DOT_PRODUCT(dro, dro)
               ! only reject if it moves away from COM outside the droplet radius
               IF (rtmpo < rtmp) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT atom
               END IF
            END IF
         END DO atom
         IF (should_reject) THEN
            ! write back only changed atoms
            DO kbead = startbead, helium%beads
               helium%work(:, bstartatom, kbead) = helium%pos(:, bstartatom, kbead)
            END DO
            DO kbead = startbead, helium%beads
               helium%work(:, fstartatom, kbead) = helium%pos(:, fstartatom, kbead)
            END DO
            helium%worm_xtra_bead_work(:) = helium%worm_xtra_bead(:)
            helium%worm_atom_idx_work = helium%worm_atom_idx
            IF (startbead > endbead) THEN
               DO kbead = 1, endbead
                  helium%work(:, fendatom, kbead) = helium%pos(:, fendatom, kbead)
               END DO
            END IF

            ! reset permtable
            tmpi = helium%permutation(bstartatom)
            helium%permutation(bstartatom) = helium%permutation(fstartatom)
            helium%permutation(fstartatom) = tmpi
            !update inverse permtable
            DO ib = 1, SIZE(helium%permutation)
               helium%iperm(helium%permutation(ib)) = ib
            END DO
            ac = 0
            RETURN
         END IF
      END IF

      ! finally accepted
      ac = 1
      DO kbead = startbead, helium%beads
         helium%pos(:, bstartatom, kbead) = helium%work(:, bstartatom, kbead)
      END DO
      DO kbead = startbead, helium%beads
         helium%pos(:, fstartatom, kbead) = helium%work(:, fstartatom, kbead)
      END DO
      helium%worm_xtra_bead(:) = helium%worm_xtra_bead_work(:)
      helium%worm_atom_idx = helium%worm_atom_idx_work
      IF (startbead > endbead) THEN
         DO kbead = 1, endbead
            helium%pos(:, fendatom, kbead) = helium%work(:, fendatom, kbead)
         END DO
      END IF

      RETURN

   END SUBROUTINE worm_swap_move

! **************************************************************************************************
!> \brief Action along path
!> \param helium ...
!> \param pos ...
!> \param startatom ...
!> \param startbead ...
!> \param endatom ...
!> \param endbead ...
!> \return ...
!> \author fuhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION worm_path_action(helium, pos, &
                                           startatom, startbead, endatom, endbead) RESULT(partaction)

      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: pos
      INTEGER, INTENT(IN)                                :: startatom, startbead, endatom, endbead

      INTEGER                                            :: iatom, ibead
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work2, work3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work

      ALLOCATE (work(3, helium%beads + 1), work2(helium%beads + 1), work3(SIZE(helium%uoffdiag, 1) + 1))
      partaction = 0.0_dp
      ! do action in two ways, depending
      ! if coordinates are not wrapping
      IF (startbead < endbead) THEN
         ! helium pair action
         ! every atom, with the one (or two) who got a resampling
         DO iatom = 1, helium%atoms
            ! avoid self interaction
            IF (iatom == startatom) CYCLE
            ! first the section up to the worm gap
            ! two less, because we need to work on the worm intersection separately
            DO ibead = startbead, endbead
               work(:, ibead + 1 - startbead) = pos(:, iatom, ibead) - pos(:, startatom, ibead)
            END DO
            partaction = partaction + helium_eval_chain(helium, work, endbead - startbead + 1, work2, work3)
         END DO
      ELSE !(startbead > endbead)
         ! helium pair action
         ! every atom, with the one (or two) who got a resampling
         DO iatom = 1, helium%atoms
            ! avoid self interaction
            IF (iatom == startatom) CYCLE
            ! first the section up to the worm gap
            ! two less, because we need to work on the worm intersection separately
            DO ibead = startbead, helium%beads
               work(:, ibead + 1 - startbead) = pos(:, iatom, ibead) - pos(:, startatom, ibead)
            END DO
            ! wrapping bond
            work(:, helium%beads + 2 - startbead) = pos(:, helium%permutation(iatom), 1) - pos(:, endatom, 1)
            partaction = partaction + helium_eval_chain(helium, work, helium%beads - startbead + 2, work2, work3)
         END DO

         ! ending atom
         DO iatom = 1, helium%atoms
            ! avoid self interaction
            IF (iatom == endatom) CYCLE
            !from first to endbead
            IF (endbead > 1) THEN
               DO ibead = 1, endbead
                  work(:, ibead) = pos(:, iatom, ibead) - pos(:, endatom, ibead)
               END DO
               partaction = partaction + helium_eval_chain(helium, work, endbead, work2, work3)
            END IF
         END DO

      END IF
      DEALLOCATE (work, work2, work3)

      RETURN

   END FUNCTION worm_path_action

! **************************************************************************************************
!> \brief Corrected path action for worm
!> \param helium ...
!> \param pos ...
!> \param startatom ...
!> \param startbead ...
!> \param endatom ...
!> \param endbead ...
!> \param xtrapos ...
!> \param worm_atom_idx ...
!> \param worm_bead_idx ...
!> \return ...
!> \author fuhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION worm_path_action_worm_corrected(helium, pos, &
                                   startatom, startbead, endatom, endbead, xtrapos, worm_atom_idx, worm_bead_idx) RESULT(partaction)

      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: pos
      INTEGER, INTENT(IN)                                :: startatom, startbead, endatom, endbead
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xtrapos
      INTEGER, INTENT(IN)                                :: worm_atom_idx, worm_bead_idx

      INTEGER                                            :: iatom, ibead
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: work2, work3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), DIMENSION(3)                        :: r, rp

      ALLOCATE (work(3, helium%beads + 1), work2(helium%beads + 1), work3(SIZE(helium%uoffdiag, 1) + 1))
      partaction = 0.0_dp
      ! do action in two ways, depending
      ! if coordinates are not wrapping
      IF (startbead < endbead) THEN
         ! helium pair action
         ! every atom, with the one (or two) who got a resampling
         DO iatom = 1, helium%atoms
            ! avoid self interaction
            IF (iatom == startatom) CYCLE
            ! first the section up to the worm gap
            ! two less, because we need to work on the worm intersection separately
            IF (worm_bead_idx - startbead > 1) THEN
               DO ibead = startbead, worm_bead_idx - 1
                  work(:, ibead + 1 - startbead) = pos(:, iatom, ibead) - pos(:, startatom, ibead)
               END DO
               partaction = partaction + helium_eval_chain(helium, work, worm_bead_idx - startbead, work2, work3)
            END IF

            IF (endbead > worm_bead_idx) THEN
               DO ibead = worm_bead_idx, endbead
                  work(:, ibead + 1 - worm_bead_idx) = pos(:, iatom, ibead) - pos(:, startatom, ibead)
               END DO
               partaction = partaction + helium_eval_chain(helium, work, endbead - worm_bead_idx + 1, work2, work3)
            END IF
         END DO

         IF (worm_atom_idx /= startatom) THEN
            DO iatom = 1, helium%atoms
               IF (iatom == startatom) CYCLE
               IF (iatom == worm_atom_idx) CYCLE
               r(:) = pos(:, iatom, worm_bead_idx - 1) - pos(:, startatom, worm_bead_idx - 1)
               rp(:) = pos(:, iatom, worm_bead_idx) - pos(:, startatom, worm_bead_idx)
               partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
            END DO
            ! add pair action with worm
            r(:) = pos(:, startatom, worm_bead_idx - 1) - pos(:, worm_atom_idx, worm_bead_idx - 1)
            rp(:) = pos(:, startatom, worm_bead_idx) - xtrapos(:)
            partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
         ELSE
            ! worm intersection
            DO iatom = 1, helium%atoms
               IF (iatom == startatom) CYCLE
               r(:) = pos(:, iatom, worm_bead_idx - 1) - pos(:, startatom, worm_bead_idx - 1)
               rp(:) = pos(:, iatom, worm_bead_idx) - xtrapos(:)
               partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
            END DO
         END IF
      ELSE !(startbead > endbead)
         ! section wraps around the atom
         ! two cases: worm gap is before or after wrap
         IF (worm_bead_idx > startbead) THEN
            ! action terms up to worm gap on starting atom
            DO iatom = 1, helium%atoms
               ! avoid self interaction
               IF (iatom == startatom) CYCLE
               ! first the section up to the worm gap
               ! two less, because we need to work on the worm intersection separately
               IF (worm_bead_idx - startbead > 1) THEN
                  DO ibead = startbead, worm_bead_idx - 1
                     work(:, ibead + 1 - startbead) = pos(:, iatom, ibead) - pos(:, startatom, ibead)
                  END DO
                  partaction = partaction + helium_eval_chain(helium, work, worm_bead_idx - startbead, work2, work3)
               END IF

               ! up to the wrapping border
               DO ibead = worm_bead_idx, helium%beads
                  work(:, ibead + 1 - worm_bead_idx) = pos(:, iatom, ibead) - pos(:, startatom, ibead)
               END DO
               ! wrapping bond
               work(:, helium%beads - worm_bead_idx + 2) = pos(:, helium%permutation(iatom), 1) - &
                                                           pos(:, helium%permutation(startatom), 1)
               partaction = partaction + helium_eval_chain(helium, work, helium%beads - worm_bead_idx + 2, work2, work3)

            END DO

            ! ending atom
            DO iatom = 1, helium%atoms
               ! avoid self interaction
               IF (iatom == endatom) CYCLE
               !from first to endbead
               IF (endbead > 1) THEN
                  DO ibead = 1, endbead
                     work(:, ibead) = pos(:, iatom, ibead) - pos(:, endatom, ibead)
                  END DO
                  partaction = partaction + helium_eval_chain(helium, work, endbead, work2, work3)
               END IF
            END DO

            IF (worm_atom_idx /= startatom) THEN
               DO iatom = 1, helium%atoms
                  IF (iatom == startatom) CYCLE
                  IF (iatom == worm_atom_idx) CYCLE
                  r(:) = pos(:, iatom, worm_bead_idx - 1) - pos(:, startatom, worm_bead_idx - 1)
                  rp(:) = pos(:, iatom, worm_bead_idx) - pos(:, startatom, worm_bead_idx)
                  partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
               END DO
               ! add pair action with worm
               r(:) = pos(:, startatom, worm_bead_idx - 1) - pos(:, worm_atom_idx, worm_bead_idx - 1)
               rp(:) = pos(:, startatom, worm_bead_idx) - xtrapos(:)
               partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
            ELSE
               ! worm intersection
               DO iatom = 1, helium%atoms
                  IF (iatom == startatom) CYCLE
                  r(:) = pos(:, iatom, worm_bead_idx - 1) - pos(:, startatom, worm_bead_idx - 1)
                  rp(:) = pos(:, iatom, worm_bead_idx) - xtrapos(:)
                  partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
               END DO
            END IF
         ELSE !(worm_bead_idx < endbead)
            IF (worm_bead_idx /= 1) THEN
               DO iatom = 1, helium%atoms
                  ! avoid self interaction
                  IF (iatom == startatom) CYCLE
                  ! first the section up to the end of the atom
                  ! one less, because we need to work on the wrapping separately
                  DO ibead = startbead, helium%beads
                     work(:, ibead + 1 - startbead) = pos(:, iatom, ibead) - pos(:, startatom, ibead)
                  END DO
                  ! wrapping bond
             work(:, helium%beads - startbead + 2) = pos(:, helium%permutation(iatom), 1) - pos(:, helium%permutation(startatom), 1)
                  partaction = partaction + helium_eval_chain(helium, work, helium%beads - startbead + 2, work2, work3)
               END DO

               ! ending atom
               DO iatom = 1, helium%atoms
                  ! avoid self interaction
                  IF (iatom == endatom) CYCLE
                  !from first to two before the worm gap
                  IF (worm_bead_idx > 2) THEN
                     DO ibead = 1, worm_bead_idx - 1
                        work(:, ibead) = pos(:, iatom, ibead) - pos(:, endatom, ibead)
                     END DO
                     partaction = partaction + helium_eval_chain(helium, work, worm_bead_idx - 1, work2, work3)
                  END IF

                  IF (endbead > worm_bead_idx) THEN
                     DO ibead = worm_bead_idx, endbead
                        work(:, ibead - worm_bead_idx + 1) = pos(:, iatom, ibead) - pos(:, endatom, ibead)
                     END DO
                     partaction = partaction + helium_eval_chain(helium, work, endbead - worm_bead_idx + 1, work2, work3)
                  END IF
               END DO

               IF (worm_atom_idx /= endatom) THEN
                  DO iatom = 1, helium%atoms
                     IF (iatom == endatom) CYCLE
                     IF (iatom == worm_atom_idx) CYCLE
                     r(:) = pos(:, iatom, worm_bead_idx - 1) - pos(:, endatom, worm_bead_idx - 1)
                     rp(:) = pos(:, iatom, worm_bead_idx) - pos(:, endatom, worm_bead_idx)
                     partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
                  END DO
                  ! add pair action with worm
                  r(:) = pos(:, endatom, worm_bead_idx - 1) - pos(:, worm_atom_idx, worm_bead_idx - 1)
                  rp(:) = pos(:, endatom, worm_bead_idx) - xtrapos(:)
                  partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
               ELSE
                  ! worm intersection
                  DO iatom = 1, helium%atoms
                     IF (iatom == endatom) CYCLE
                     r(:) = pos(:, iatom, worm_bead_idx - 1) - pos(:, endatom, worm_bead_idx - 1)
                     rp(:) = pos(:, iatom, worm_bead_idx) - xtrapos(:)
                     partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
                  END DO
               END IF
            ELSE !(worm_bead_idx == 1)
               ! special treatment if first bead is worm bead
               DO iatom = 1, helium%atoms
                  ! avoid self interaction
                  IF (iatom == startatom) CYCLE
                  ! first the section up to the end of the atom
                  ! one less, because we need to work on the wrapping separately
                  IF (helium%beads > startbead) THEN
                     DO ibead = startbead, helium%beads
                        work(:, ibead - startbead + 1) = pos(:, iatom, ibead) - pos(:, startatom, ibead)
                     END DO
                     partaction = partaction + helium_eval_chain(helium, work, helium%beads - startbead + 1, work2, work3)
                  END IF
               END DO

               ! ending atom
               DO iatom = 1, helium%atoms
                  IF (endbead > 1) THEN
                     DO ibead = 1, endbead
                        work(:, ibead) = pos(:, iatom, ibead) - pos(:, endatom, ibead)
                     END DO
                     partaction = partaction + helium_eval_chain(helium, work, endbead, work2, work3)
                  END IF
               END DO

               IF (worm_atom_idx /= endatom) THEN
                  DO iatom = 1, helium%atoms
                     IF (iatom == endatom) CYCLE
                     IF (iatom == worm_atom_idx) CYCLE
                     r(:) = pos(:, helium%iperm(iatom), helium%beads) - pos(:, startatom, helium%beads)
                     rp(:) = pos(:, iatom, 1) - pos(:, endatom, 1)
                     partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
                  END DO
                  ! add pair action with worm
                  r(:) = pos(:, startatom, helium%beads) - pos(:, helium%iperm(worm_atom_idx), helium%beads)
                  rp(:) = pos(:, endatom, 1) - xtrapos(:)
                  partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
               ELSE
                  ! worm intersection
                  DO iatom = 1, helium%atoms
                     IF (iatom == endatom) CYCLE
                     r(:) = pos(:, helium%iperm(iatom), helium%beads) - pos(:, startatom, helium%beads)
                     rp(:) = pos(:, iatom, 1) - xtrapos(:)
                     partaction = partaction + helium_eval_expansion(helium, r, rp, work3)
                  END DO
               END IF
            END IF
         END IF
      END IF
      DEALLOCATE (work, work2, work3)

      RETURN

   END FUNCTION worm_path_action_worm_corrected

! **************************************************************************************************
!> \brief Path interaction
!> \param pint_env ...
!> \param helium ...
!> \param pos ...
!> \param startatom ...
!> \param startbead ...
!> \param endatom ...
!> \param endbead ...
!> \return ...
!> \author fuhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION worm_path_inter_action(pint_env, helium, pos, &
                                                 startatom, startbead, endatom, endbead) RESULT(partaction)

      TYPE(pint_env_type), INTENT(IN), POINTER           :: pint_env
      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: pos
      INTEGER, INTENT(IN)                                :: startatom, startbead, endatom, endbead

      INTEGER                                            :: ibead
      REAL(KIND=dp)                                      :: energy

      partaction = 0.0_dp

      ! helium interaction with the solute
      ! do action in two ways, depending
      ! if coordinates are not wrapping
      IF (startbead < endbead) THEN

         ! interaction is only beadwise due to primitive coupling
         ! startatom and endatom are the same
         DO ibead = startbead + 1, endbead - 1
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        startatom, ibead, pos(:, startatom, ibead), energy=energy)
            partaction = partaction + energy
         END DO

      ELSE !(startbead > endbead)

         ! interaction is only beadwise due to primitive coupling
         ! startatom and endatom are different
         DO ibead = startbead + 1, helium%beads
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        startatom, ibead, pos(:, startatom, ibead), energy=energy)
            partaction = partaction + energy
         END DO
         ! second atom after imaginary time wrap
         DO ibead = 1, endbead - 1
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        endatom, ibead, pos(:, endatom, ibead), energy=energy)
            partaction = partaction + energy
         END DO
      END IF

      partaction = partaction*helium%tau

      RETURN

   END FUNCTION worm_path_inter_action

! **************************************************************************************************
!> \brief Path interaction for head
!> \param pint_env ...
!> \param helium ...
!> \param pos ...
!> \param startatom ...
!> \param startbead ...
!> \param xtrapos ...
!> \param worm_atom_idx ...
!> \param worm_bead_idx ...
!> \return ...
!> \author fuhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION worm_path_inter_action_head(pint_env, helium, pos, &
                                                     startatom, startbead, xtrapos, worm_atom_idx, worm_bead_idx) RESULT(partaction)

      TYPE(pint_env_type), INTENT(IN), POINTER           :: pint_env
      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: pos
      INTEGER, INTENT(IN)                                :: startatom, startbead
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xtrapos
      INTEGER, INTENT(IN)                                :: worm_atom_idx, worm_bead_idx

      INTEGER                                            :: ibead
      REAL(KIND=dp)                                      :: energy

      partaction = 0.0_dp

      ! helium interaction with the solute
      ! if coordinates are not wrapping
      IF (startbead < worm_bead_idx) THEN
         DO ibead = startbead + 1, worm_bead_idx - 1
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        startatom, ibead, pos(:, startatom, ibead), energy=energy)
            partaction = partaction + energy
         END DO
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     startatom, ibead, xtrapos, energy=energy)
         partaction = partaction + 0.5_dp*energy
      ELSE !(startbead > worm_bead_idx)
         DO ibead = startbead + 1, helium%beads
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        startatom, ibead, pos(:, startatom, ibead), energy=energy)
            partaction = partaction + energy
         END DO
         DO ibead = 1, worm_bead_idx - 1
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        worm_atom_idx, ibead, pos(:, worm_atom_idx, ibead), energy=energy)
            partaction = partaction + energy
         END DO
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     worm_atom_idx, ibead, xtrapos, energy=energy)
         partaction = partaction + 0.5_dp*energy
      END IF

      partaction = partaction*helium%tau

      RETURN

   END FUNCTION worm_path_inter_action_head

! **************************************************************************************************
!> \brief Path interaction for tail
!> \param pint_env ...
!> \param helium ...
!> \param pos ...
!> \param endatom ...
!> \param endbead ...
!> \param worm_atom_idx ...
!> \param worm_bead_idx ...
!> \return ...
!> \author fuhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION worm_path_inter_action_tail(pint_env, helium, pos, &
                                                      endatom, endbead, worm_atom_idx, worm_bead_idx) RESULT(partaction)

      TYPE(pint_env_type), INTENT(IN), POINTER           :: pint_env
      TYPE(helium_solvent_type), INTENT(IN), POINTER     :: helium
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         POINTER                                         :: pos
      INTEGER, INTENT(IN)                                :: endatom, endbead, worm_atom_idx, &
                                                            worm_bead_idx

      INTEGER                                            :: ibead
      REAL(KIND=dp)                                      :: energy

      partaction = 0.0_dp

      ! helium interaction with the solute
      ! if coordinates are not wrapping
      IF (worm_bead_idx < endbead) THEN
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     worm_atom_idx, worm_bead_idx, pos(:, worm_atom_idx, worm_bead_idx), energy=energy)
         partaction = partaction + 0.5_dp*energy
         DO ibead = worm_bead_idx + 1, endbead - 1
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        endatom, ibead, pos(:, endatom, ibead), energy=energy)
            partaction = partaction + energy
         END DO
      ELSE !(worm_bead_idx > endbead)
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     worm_atom_idx, worm_bead_idx, pos(:, worm_atom_idx, worm_bead_idx), energy=energy)
         partaction = partaction + 0.5_dp*energy
         DO ibead = worm_bead_idx + 1, helium%beads
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        worm_atom_idx, ibead, pos(:, worm_atom_idx, ibead), energy=energy)
            partaction = partaction + energy
         END DO
         DO ibead = 1, endbead - 1
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        endatom, ibead, pos(:, endatom, ibead), energy=energy)
            partaction = partaction + energy
         END DO
      END IF

      partaction = partaction*helium%tau

      RETURN

   END FUNCTION worm_path_inter_action_tail

END MODULE helium_worm
