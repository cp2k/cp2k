!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Worm algorithm
!> \author Ricardo
!> \date   2014-10-15
! *****************************************************************************
MODULE helium_worm

   USE helium_common,                   ONLY: helium_calc_plength,&
                                              helium_com,&
                                              helium_cycle_number,&
                                              helium_eval_expansion,&
                                              helium_pbc,&
                                              helium_total_moment_of_inertia,&
                                              helium_total_projected_area,&
                                              helium_total_winding_number
   USE helium_interactions,             ONLY: helium_bead_solute_e_f,&
                                              helium_calc_energy,&
                                              helium_solute_e_f
   USE helium_io,                       ONLY: helium_write_line
   USE helium_types,                    ONLY: he_mass,&
                                              helium_solvent_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE parallel_rng_types,              ONLY: next_random_number
   USE physcon,                         ONLY: massunit
   USE pint_types,                      ONLY: pint_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_worm'

   PUBLIC :: helium_sample_worm

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param pint_env ...
!> \param nMC ...
!> \param nstat ...
! **************************************************************************************************
   SUBROUTINE helium_sample_worm(helium, pint_env, nMC, nstat)
      TYPE(helium_solvent_type), POINTER                 :: helium
      TYPE(pint_env_type), POINTER                       :: pint_env
      INTEGER, INTENT(IN)                                :: nMC
      INTEGER                                            :: nstat

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_sample_worm', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: stmp
      INTEGER :: ac, ib, iMC, imoveG, ip, irep, itry, iupdate, l, natt, ncentracc, ncentratt, &
         ncloseacc, ncloseatt, nG, nmoveheadacc, nmoveheadatt, nmovetailacc, nmovetailatt, &
         nopenacc, nopenatt, nstagacc, nstagatt, nswapacc, nswapatt, ntot, nwrong, nZ
      REAL(KIND=dp)                                      :: rtmp

      helium%center(:) = helium_com(helium)

      nwrong = 0
      natt = 25
      DO itry = 1, natt !We give it several attempts to sample the Z-sector
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
         nswapatt = 0
         nswapacc = 0
         nstat = 0
         nZ = 0
         nG = 0
         ntot = 0
         helium%proarea%inst(:) = 0.d0
         helium%prarea2%inst(:) = 0.d0
         helium%wnumber%inst(:) = 0.d0
         helium%wnmber2%inst(:) = 0.d0
         helium%mominer%inst(:) = 0.d0
         IF (helium%solute_present) helium%force_avrg(:, :) = 0.0d0
         helium%energy_avrg(:) = 0.0d0
         helium%plength_avrg(:) = 0.0d0

         DO iMC = 1, nMC

            !centroid move:
            IF (MOD(iMC, helium%worm_centroid_move_freq) == 0) THEN
               DO ip = 1, helium%atoms
                  CALL centroid_move(pint_env, helium, ip, helium%worm_centroid_move_drmax, ac)
                  ntot = ntot+1
                  nG = nG+helium%isector
                  ncentratt = ncentratt+1
                  IF (ac == 1) THEN
                     ncentracc = ncentracc+1
                     helium%center(:) = helium_com(helium)
                  END IF
               END DO
            END IF

            !staging move:
            DO irep = 1, helium%worm_staging_move_rep
               DO ip = 1, helium%atoms
                  ib = next_random_number(helium%rng_stream_uniform, 1, helium%beads)
                  CALL staging(pint_env, helium, ip, ib, helium%worm_staging_move_l, ac)
                  ntot = ntot+1
                  nG = nG+helium%isector
                  nstagatt = nstagatt+1
                  IF (ac == 1) THEN
                     nstagacc = nstagacc+1
                     helium%center(:) = helium_com(helium)
                  END IF
               END DO
            END DO

            !G-sector moves:
            IF (helium%worm_gsector_move .EQV. .TRUE.) THEN
               !open vs close:
               iupdate = next_random_number(helium%rng_stream_uniform, 0, 1) ! page 129 of Rota's thesis
               IF (iupdate == 0) THEN
                  IF (helium%isector == 0) THEN
                     !open move
                     ip = next_random_number(helium%rng_stream_uniform, 1, helium%atoms)
                     l = next_random_number(helium%rng_stream_uniform, 1, helium%worm_open_close_move_lmax)
                     CALL OPEN (pint_env, helium, ip, l, helium%worm_open_close_move_C, ac)
                     ntot = ntot+1
                     nG = nG+helium%isector
                     nopenatt = nopenatt+1
                     IF (ac == 1) THEN
                        nopenacc = nopenacc+1
                        helium%center(:) = helium_com(helium)
                     END IF
                  END IF
               ELSE
                  IF (helium%isector == 1) THEN
                     !close move:
                     l = next_random_number(helium%rng_stream_uniform, 1, helium%worm_open_close_move_lmax)
                     CALL CLOSE (pint_env, helium, l, helium%worm_open_close_move_C, ac)
                     ntot = ntot+1
                     nG = nG+helium%isector
                     ncloseatt = ncloseatt+1
                     IF (ac == 1) THEN
                        ncloseacc = ncloseacc+1
                        helium%center(:) = helium_com(helium)
                     END IF
                  END IF
               END IF
               IF (helium%isector == 1) THEN
                  DO imoveG = 1, helium%worm_gsector_move_rep
                     !head nove:
                     l = next_random_number(helium%rng_stream_uniform, 1, &
                                            helium%worm_head_tail_move_lmax)
                     CALL movehead(pint_env, helium, l, ac)
                     ntot = ntot+1
                     nG = nG+helium%isector
                     nmoveheadatt = nmoveheadatt+1
                     IF (ac == 1) THEN
                        nmoveheadacc = nmoveheadacc+1
                        helium%center(:) = helium_com(helium)
                     END IF
                     !tail move:
                     l = next_random_number(helium%rng_stream_uniform, 1, &
                                            helium%worm_head_tail_move_lmax)
                     CALL movetail(pint_env, helium, l, ac)
                     ntot = ntot+1
                     nG = nG+helium%isector
                     nmovetailatt = nmovetailatt+1
                     IF (ac == 1) THEN
                        nmovetailacc = nmovetailacc+1
                        helium%center(:) = helium_com(helium)
                     END IF
                     !swap move:
                     l = next_random_number(helium%rng_stream_uniform, 1, &
                                            helium%worm_swap_move_lmax)
                     CALL swap(pint_env, helium, l, ac)
                     ntot = ntot+1
                     nG = nG+helium%isector
                     nswapatt = nswapatt+1
                     IF (ac == 1) THEN
                        nswapacc = nswapacc+1
                        helium%center(:) = helium_com(helium)
                     END IF
                  END DO
               END IF
            END IF

            ! Accumulate statistics if we are in the Z-sector:
            IF (helium%isector == 0) THEN
               nstat = nstat+1

               CALL helium_calc_energy(helium, pint_env)
               helium%energy_avrg(:) = helium%energy_avrg(:)+helium%energy_inst(:)

               IF (helium%solute_present) THEN
                  CALL helium_solute_e_f(pint_env, helium, rtmp)
                  helium%force_avrg(:, :) = helium%force_avrg(:, :)+helium%force_inst(:, :)
               END IF

               IF (helium%worm_gsector_move .EQV. .TRUE.) THEN
!TODO_WORM        CALL helium_calc_wnumber( helium )
                  CALL helium_calc_plength(helium)
                  helium%plength_avrg(:) = helium%plength_avrg(:)+helium%plength_inst(:)
                  helium%proarea%inst(:) = helium_total_projected_area(helium)
                  helium%prarea2%inst(:) = helium%proarea%inst(:)*helium%proarea%inst(:)
                  helium%wnumber%inst(:) = helium_total_winding_number(helium)
                  helium%wnmber2%inst(:) = helium%wnumber%inst(:)*helium%wnumber%inst(:)
                  helium%mominer%inst(:) = helium_total_moment_of_inertia(helium)
               END IF

               ! instantaneous force output according to HELIUM%PRINT%FORCES_INST
               ! Warning: file I/O here may cost A LOT of cpu time!
               ! switched off here to save cpu
               !CALL helium_write_force_inst( helium )

            END IF

         END DO ! MC loop

         ! if most of the whole MC loop is spent in the G-sector (>90%), repeat loop.
         ! If after natt times we are still stuck in the G-sector, stop and kindly
         ! suggest changing parameters.
         IF (nG > 0.9*ntot) THEN
            nwrong = nwrong+1
            IF (nwrong >= natt) THEN
               stmp = ""
               WRITE (stmp, *) "Failed to sample Z sector completely!"
               CALL helium_write_line(stmp)
               stmp = ""
               WRITE (stmp, *) "Try changing C parameter or increase N_OUTER"
               CALL helium_write_line(stmp)
               CPABORT("Stopping.")
            END IF
         ELSE
            EXIT
         END IF

      END DO !attempts loop

      ! STATISTICS
      IF (helium%worm_write .EQV. .TRUE.) THEN
         stmp = ""
         WRITE (stmp, *) '--------------------------------------------------'
         CALL helium_write_line(stmp)
         stmp = ""
         WRITE (stmp, *) 'Statistics after ', nwrong, ' failed attempts'
         CALL helium_write_line(stmp)
         IF (ncentratt > 0) THEN
            stmp = ""
            WRITE (stmp, '(a10,f15.5,2i10)') 'Centroid', REAL(ncentracc, KIND=dp)/REAL(ncentratt, KIND=dp), ncentracc, ncentratt
            CALL helium_write_line(stmp)
         END IF
         IF (nstagatt > 0) THEN
            stmp = ""
            WRITE (stmp, '(a10,f15.5,2i10)') 'Staging', REAL(nstagacc, KIND=dp)/REAL(nstagatt, KIND=dp), nstagacc, nstagatt
            CALL helium_write_line(stmp)
         END IF
         IF (nopenatt > 0) THEN
            stmp = ""
            WRITE (stmp, '(a10,f15.5,2i10)') 'Open', REAL(nopenacc, KIND=dp)/REAL(nopenatt, KIND=dp), nopenacc, nopenatt
            CALL helium_write_line(stmp)
         END IF
         IF (ncloseatt > 0) THEN
            stmp = ""
            WRITE (stmp, '(a10,f15.5,2i10)') 'Close', REAL(ncloseacc, KIND=dp)/REAL(ncloseatt, KIND=dp), ncloseacc, ncloseatt
            CALL helium_write_line(stmp)
         END IF
         IF (nmoveheadatt > 0) THEN
            stmp = ""
            WRITE (stmp, '(a10,f15.5,2i10)') 'Move_head', REAL(nmoveheadacc, KIND=dp)/REAL(nmoveheadatt, KIND=dp), &
               nmoveheadacc, nmoveheadatt
            CALL helium_write_line(stmp)
         END IF
         IF (nmovetailatt > 0) THEN
            stmp = ""
            WRITE (stmp, '(a10,f15.5,2i10)') 'Move_tail', REAL(nmovetailacc, KIND=dp)/REAL(nmovetailatt, KIND=dp), &
               nmovetailacc, nmovetailatt
            CALL helium_write_line(stmp)
         END IF
         IF (nswapatt > 0) THEN
            stmp = ""
            WRITE (stmp, '(a10,f15.5,2i10)') 'Swap', REAL(nswapacc, KIND=dp)/REAL(nswapatt, KIND=dp), nswapacc, nswapatt
            CALL helium_write_line(stmp)
         END IF

         nZ = ntot-nG

         stmp = ""
         WRITE (stmp, '(a10,f5.1,a3,f5.1)') 'Z/G: ', REAL(nZ, KIND=dp)*100.d0/REAL(ntot, KIND=dp), &
            ' / ', REAL(nG, KIND=dp)*100.d0/REAL(ntot, KIND=dp)
         CALL helium_write_line(stmp)
         stmp = ""
         WRITE (stmp, '(a20,i5,a3,i5)') 'Z-sector statistics:', nstat, ' / ', nMC
         CALL helium_write_line(stmp)
      END IF

      CALL center_pos(helium)
      RETURN
   END SUBROUTINE helium_sample_worm

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \param ip ...
!> \param ibi ...
!> \param l ...
!> \param ac ...
! **************************************************************************************************
   SUBROUTINE staging(pint_env, helium, ip, ibi, l, ac)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: ip, ibi, l
      INTEGER, INTENT(out)                               :: ac

      INTEGER                                            :: ib, ibb, ibf, ipp, nbeads
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp)                                      :: det, dU, mass, U, U1, U2, xrr

      perm => helium%permutation
      nbeads = helium%beads
      mass = he_mass*massunit
      ac = 0
      ibf = ibi+l
      ! Careful with cyclicity in open paths:
      IF (helium%isector == 1 .AND. ip == helium%iw .AND. ibi+l > nbeads) RETURN

      ! Compute the potential action of the moving beads before any trial move
      U1 = 0.d0
      DO ib = ibi, ibf-1 !only moving beads
         ipp = ip
         ibb = ib
         IF (ibb > nbeads) THEN
            ipp = perm(ip)
            ibb = ibb-nbeads
         END IF
         CALL compute_potential_action_link(pint_env, helium, ipp, ibb, U)
         U1 = U1+U
      END DO

      ! Staging move:
      CALL staging_move(helium, ip, ibi, l)
      ! cyclic symmetry when we are in the Z-sector, or in the G-sector and we are not dealing with the worm:
      IF (helium%isector == 0 .OR. (helium%isector == 1 .AND. ip .NE. helium%iw)) THEN
         helium%work(:, ip, nbeads+1) = helium%work(:, perm(ip), 1)
      END IF

      ! if no PBC apply confinement:
      IF (.NOT. helium%periodic) THEN
         IF (.NOT. indroplet(helium, ip, ibi, l)) THEN
            ! reject move
            helium%work = helium%pos
            ac = 0
            RETURN
         END IF
      END IF

      ! Compute the potential action after the move:
      U2 = 0.d0
      DO ib = ibi, ibf-1
         ipp = ip
         ibb = ib
         IF (ibb > nbeads) THEN
            ipp = perm(ip)
            ibb = ibb-nbeads
         END IF
         CALL compute_potential_action_link(pint_env, helium, ipp, ibb, U)
         U2 = U2+U
      END DO

      ! ---Metropolis algorithm---
      dU = U2-U1

      IF (dU < 0.d0) THEN
         !accept move
         helium%pos = helium%work
         ac = 1
         RETURN
      END IF
      dU = MIN(MAX(-200.0_dp, dU), 200.0_dp) ! guard against overflow
      det = dexp(-dU)
      xrr = next_random_number(helium%rng_stream_uniform)
      IF (det > xrr) THEN
         ! accept move
         helium%pos = helium%work
         ac = 1
         RETURN
      ELSE
         ! reject move
         helium%work = helium%pos
         ac = 0
         RETURN
      END IF
   END SUBROUTINE staging

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param ip ...
!> \param ibi ...
!> \param l ...
! **************************************************************************************************
   SUBROUTINE staging_move(helium, ip, ibi, l)
      ! Updates l-1 positions, from pos(:,ibi+1,ip) to pos(:,ibi+l-1,ip)
      ! No cyclicity symmetry is taken into account here!
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: ip, ibi, l

      INTEGER                                            :: ibf, ibk, ibk1, ibkf, id, ipk, ipk1, &
                                                            ipkf, k, nbeads
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp)                                      :: dr(3), mass, mk, rf(3), ri(3), rp(3), &
                                                            xlk, xlk1, xr(3)

      perm => helium%permutation
      mass = he_mass*massunit

      ibf = ibi+l
      nbeads = helium%beads
      ibkf = ibf
      ipkf = ip
      IF (ibkf > nbeads+1) THEN
         ibkf = ibkf-nbeads
         ipkf = perm(ip)
      ENDIF ! Final point(does not depend on k)
      rf(:) = helium%work(:, ipkf, ibkf)

      DO k = 0, l-2
         xlk = REAL(l-k, KIND=dp)
         xlk1 = REAL(l-k-1, KIND=dp)
         mk = mass*xlk/xlk1
         ibk = ibi+k
         ipk = ip
         IF (ibk > nbeads) THEN
            ibk = ibk-nbeads
            ipk = perm(ip)
         ENDIF
         ibk1 = ibi+k+1
         ipk1 = ip
         IF (ibk1 > nbeads) THEN
            ibk1 = ibk1-nbeads
            ipk1 = perm(ip)
         ENDIF

         ri(:) = helium%work(:, ipk, ibk)
         dr(:) = ri(:)-rf(:)
         CALL helium_pbc(helium, dr)
         ri(:) = rf(:)+dr(:)
         rp(:) = (rf(:)+xlk1*ri(:))/xlk
         CALL helium_pbc(helium, rp)

         DO id = 1, 3
            xr(id) = next_random_number(rng_stream=helium%rng_stream_gaussian, variance=helium%tau/mk)
         END DO

         CALL helium_pbc(helium, xr)
         rp(:) = rp(:)+xr(:)
         CALL helium_pbc(helium, rp)
         helium%work(:, ipk1, ibk1) = rp(:)
      END DO
   END SUBROUTINE staging_move

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param r0 ...
!> \param rl ...
!> \param l ...
!> \param r ...
! **************************************************************************************************
   SUBROUTINE staging_reconstruct(helium, r0, rl, l, r)
      ! Reconstruct a segment of length (l-1) between the fixed points r1 and rl
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      REAL(KIND=dp), INTENT(in)                          :: r0(3), rl(3)
      INTEGER, INTENT(in)                                :: l
      REAL(KIND=dp), INTENT(out)                         :: r(0:l, 3)

      INTEGER                                            :: id, k
      REAL(KIND=dp)                                      :: dr(3), mass, mk, ri(3), rp(3), xlk, &
                                                            xlk1, xr(3)

      mass = he_mass*massunit
      ! Fixed points
      r(0, :) = r0(:)
      r(l, :) = rl(:)

      DO k = 0, l-2
         xlk = REAL(l-k, KIND=dp)
         xlk1 = REAL(l-k-1, KIND=dp)
         mk = mass*xlk/xlk1
         ri(:) = r(k, :)
         dr(:) = ri(:)-rl(:)
         CALL helium_pbc(helium, dr)
         ri(:) = rl(:)+dr(:)
         rp(:) = (rl(:)+xlk1*ri(:))/xlk
         CALL helium_pbc(helium, rp)
         DO id = 1, 3
            xr(id) = next_random_number(rng_stream=helium%rng_stream_gaussian, variance=helium%tau/mk)
         END DO
         CALL helium_pbc(helium, xr)
         rp(:) = rp(:)+xr(:)
         CALL helium_pbc(helium, rp)
         r(k+1, :) = rp(:)
      END DO

   END SUBROUTINE staging_reconstruct

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \param ip ...
!> \param l ...
!> \param C ...
!> \param ac ...
! **************************************************************************************************
   SUBROUTINE OPEN (pint_env, helium, ip, l, C, ac)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: ip, l
      REAL(KIND=dp), INTENT(in)                          :: C
      INTEGER, INTENT(out)                               :: ac

      INTEGER                                            :: ib, id, nbeads
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp)                                      :: dU, expdU, head(3), mass, Popen, rho, &
                                                            rnew(3), tail(3), U, U1, U2, xr, xrr

      perm => helium%permutation
      nbeads = helium%beads
      mass = he_mass*massunit

      ! Computes the potential action before the move:
      U1 = 0.d0
      DO ib = nbeads+1-l, nbeads
         CALL compute_potential_action_link(pint_env, helium, ip, ib, U)
         U1 = U1+U
      END DO

      ! Sample the proposal of the head of the worm:
      DO id = 1, 3
         xr = next_random_number(rng_stream=helium%rng_stream_gaussian, &
                                 variance=helium%tau*l/mass)
         rnew(id) = helium%work(id, ip, nbeads+1-l)+xr
      END DO

      CALL helium_pbc(helium, rnew)
      helium%work(:, ip, nbeads+1) = rnew(:)

      head(:) = helium%work(:, ip, helium%beads+1)
      tail(:) = helium%work(:, perm(ip), 1)
      IF (distance2(head, tail)/(4.d0*helium%beads*helium%tau*helium%hb2m) > 4) THEN
         !reject move
         helium%work = helium%pos
         ac = 0
         RETURN
      END IF

      ! Perform a staging move from time slices
      ! nbeads+1-l (fixed) and the new head at nbeads+1 (fixed):
      CALL staging_move(helium, ip, nbeads+1-l, l)

      ! if no PBC apply confinement:
      IF (.NOT. helium%periodic) THEN
      IF (.NOT. indroplet(helium, ip, nbeads+1-l, l)) THEN
         ! reject move
         helium%work = helium%pos
         ac = 0
         RETURN
      END IF
      END IF

      ! Compute the potential action after the move:
      U2 = 0.d0
      DO ib = nbeads+1-l, nbeads
         CALL compute_potential_action_link(pint_env, helium, ip, ib, U)
         U2 = U2+U
      END DO

      ! Metropolis:
      dU = U2-U1
      dU = MIN(MAX(-200.0_dp, dU), 200.0_dp) ! guard against overflow
      expdU = dexp(-dU)
      ! Free density matrix between the two time slices:
      rho = free_density_matrix(helium, ip, nbeads+1-l, ip, nbeads+1)
      Popen = C*expdU/rho
      xrr = next_random_number(helium%rng_stream_uniform)

      IF (Popen > xrr) THEN !accept move
         helium%pos = helium%work
         helium%isector = 1 !change to the G-sector
         helium%iw = ip !define the worm particle
         helium%it = perm(helium%iw)
         ac = 1
      ELSE !reject move
         helium%work = helium%pos
         ac = 0
      END IF

   END SUBROUTINE OPEN

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param C ...
!> \param ac ...
! **************************************************************************************************
   SUBROUTINE CLOSE (pint_env, helium, l, C, ac)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: l
      REAL(KIND=dp), INTENT(in)                          :: C
      INTEGER, INTENT(out)                               :: ac

      INTEGER                                            :: ib, nbeads
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp)                                      :: dr(3), dU, expdU, head(3), mass, Pclose, &
                                                            rho, rnew(3), tail(3), U, U1, U2, xrr

      perm => helium%permutation
      nbeads = helium%beads
      mass = he_mass*massunit

      head(:) = helium%work(:, helium%iw, helium%beads+1)
      tail(:) = helium%work(:, helium%it, 1)
      IF (distance2(head, tail)/(4.d0*helium%beads*helium%tau*helium%hb2m) > 4) THEN
         !reject move
         ac = 0
         RETURN
      END IF

      ! Compute the potential action before the move:
      U1 = 0.d0
      DO ib = nbeads+1-l, nbeads
         CALL compute_potential_action_link(pint_env, helium, helium%iw, ib, U)
         U1 = U1+U
      END DO

      ! Try closing the path:
      dr(:) = helium%work(:, perm(helium%iw), 1)-helium%work(:, helium%iw, nbeads+1)
      CALL helium_pbc(helium, dr)
      rnew(:) = helium%work(:, helium%iw, nbeads+1)+dr(:)
      CALL helium_pbc(helium, rnew)
      helium%work(:, helium%iw, nbeads+1) = rnew(:)

      CALL staging_move(helium, helium%iw, nbeads+1-l, l)

      ! if no PBC apply confinement:
      IF (.NOT. helium%periodic) THEN
      IF (.NOT. indroplet(helium, helium%iw, nbeads+1-l, l)) THEN
         ! reject move
         helium%work = helium%pos
         ac = 0
         RETURN
      END IF
      END IF

      ! Compute the potential action after the move:
      U2 = 0.d0
      DO ib = nbeads+1-l, nbeads
         CALL compute_potential_action_link(pint_env, helium, helium%iw, ib, U)
         U2 = U2+U
      END DO

      ! Metropolis:
      dU = U2-U1
      dU = MIN(MAX(-200.0_dp, dU), 200.0_dp) ! guard against overflow
      expdU = dexp(-dU)
      ! Free density matrix between the two time slices:
      rho = free_density_matrix(helium, helium%iw, nbeads+1, helium%iw, nbeads+1-l)
      Pclose = rho*expdU/C
      xrr = next_random_number(helium%rng_stream_uniform)
      IF (Pclose > xrr) THEN !accept move
         helium%pos = helium%work
         helium%isector = 0
         helium%iw = 0
         helium%it = 0
         ac = 1
      ELSE !reject move
         helium%work = helium%pos
         ac = 0
      END IF
   END SUBROUTINE CLOSE

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param ac ...
! **************************************************************************************************
   SUBROUTINE movehead(pint_env, helium, l, ac)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: l
      INTEGER, INTENT(out)                               :: ac

      INTEGER                                            :: ib, id, nbeads
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp)                                      :: dU, expdU, mass, Pmovehead, rnew(3), U, &
                                                            U1, U2, xr

      perm => helium%permutation
      nbeads = helium%beads
      mass = he_mass*massunit

      ! Compute the potential action before the move:
      U1 = 0.d0
      DO ib = nbeads+1-l, nbeads
         CALL compute_potential_action_link(pint_env, helium, helium%iw, ib, U)
         U1 = U1+U
      END DO

      ! Sample the proposal of the new head of the worm:
      DO id = 1, 3
         xr = next_random_number(rng_stream=helium%rng_stream_gaussian, &
                                 variance=helium%tau*l/mass)
         rnew(id) = helium%work(id, helium%iw, nbeads+1-l)+xr
      END DO
      CALL helium_pbc(helium, rnew)
      helium%work(:, helium%iw, nbeads+1) = rnew(:)

      CALL staging_move(helium, helium%iw, nbeads+1-l, l)

      ! if no PBC apply confinement:
      IF (.NOT. helium%periodic) THEN
      IF (.NOT. indroplet(helium, helium%iw, nbeads+1-l, l)) THEN
         ! reject move
         helium%work = helium%pos
         ac = 0
         RETURN
      END IF
      END IF

      ! Compute the potential action after the move:
      U2 = 0.d0
      DO ib = nbeads+1-l, nbeads
         CALL compute_potential_action_link(pint_env, helium, helium%iw, ib, U)
         U2 = U2+U
      END DO

      ! Metropolis:
      IF (U2 < U1) THEN !accept move
         helium%pos = helium%work
         ac = 1
         RETURN
      END IF
      dU = U2-U1
      dU = MIN(MAX(-200.0_dp, dU), 200.0_dp) ! guard against overflow
      expdU = dexp(-dU)
      Pmovehead = expdU
      xr = next_random_number(helium%rng_stream_uniform)
      IF (Pmovehead > xr) THEN !accept move
         helium%pos = helium%work
         ac = 1
      ELSE !reject move
         helium%work = helium%pos
         ac = 0
      END IF
   END SUBROUTINE movehead

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param ac ...
! **************************************************************************************************
   SUBROUTINE movetail(pint_env, helium, l, ac)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: l
      INTEGER, INTENT(out)                               :: ac

      INTEGER                                            :: ib, id, nbeads
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp)                                      :: dU, expdU, mass, Pmovetail, rnew(3), U, &
                                                            U1, U2, xr

      perm => helium%permutation
      nbeads = helium%beads
      mass = he_mass*massunit
      !it=perm(helium%iw) !tail of the worm

      ! Compute the potential action before the move:
      U1 = 0.d0
      DO ib = 1, l
         CALL compute_potential_action_link(pint_env, helium, helium%it, ib, U)
         U1 = U1+U
      END DO

      ! Sample the proposal of the new tail of the worm:
      DO id = 1, 3
         xr = next_random_number(rng_stream=helium%rng_stream_gaussian, &
                                 variance=helium%tau*l/mass)
         rnew(id) = helium%work(id, helium%it, l+1)+xr
      END DO

      CALL helium_pbc(helium, rnew)
      helium%work(:, helium%it, 1) = rnew(:)

      CALL staging_move(helium, helium%it, 1, l)

      ! if no PBC apply confinement:
      IF (.NOT. helium%periodic) THEN
      IF (.NOT. indroplet(helium, helium%it, 1, l)) THEN
         ! reject move
         helium%work = helium%pos
         ac = 0
         RETURN
      END IF
      END IF

      ! Compute the potential action after the move:
      U2 = 0.d0
      DO ib = 1, l
         CALL compute_potential_action_link(pint_env, helium, helium%it, ib, U)
         U2 = U2+U
      END DO

      ! Metropolis:
      IF (U2 < U1) THEN !accept move
         helium%pos = helium%work
         ac = 1
         RETURN
      END IF
      dU = U2-U1
      dU = MIN(MAX(-200.0_dp, dU), 200.0_dp) ! guard against overflow
      expdU = dexp(-dU)
      Pmovetail = expdU
      xr = next_random_number(helium%rng_stream_uniform)
      IF (Pmovetail > xr) THEN !accept move
         helium%pos = helium%work
         ac = 1
      ELSE !reject move
         helium%work = helium%pos
         ac = 0
      END IF
   END SUBROUTINE movetail

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \param l ...
!> \param ac ...
! **************************************************************************************************
   SUBROUTINE swap(pint_env, helium, l, ac)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: l
      INTEGER, INTENT(out)                               :: ac

      INTEGER                                            :: ia, ib, ik, ip, natoms, nbeads
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp)                                      :: dr(3), dU, expdU, mass, Pswap1, Pswap2, &
                                                            rnew(3), sumk, sumw, U, U1, U2, xrr, &
                                                            xx, Pcum(helium%atoms), P(helium%atoms)

      perm => helium%permutation
      natoms = helium%atoms
      nbeads = helium%beads
      mass = he_mass*massunit
      ac = 0

      !Compute probabilities:
      sumw = 0.d0
      Pcum = 0.d0
      DO ia = 1, natoms
         xx = free_density_matrix(helium, ia, l+1, helium%iw, nbeads+1, l)
         P(ia) = xx
         sumw = sumw+xx
         Pcum(ia) = sumw
      END DO
      !Normalilze cumulative probability:
      Pcum = Pcum/sumw

      !Particle ik is chosen with probability P(ik)/sumw (tower sampling):
      xrr = next_random_number(helium%rng_stream_uniform)
      ik = 0
      DO ip = 1, natoms
         IF (Pcum(ip) > xrr) THEN
            ik = ip
            EXIT
         END IF
      END DO

      !First test: ik cannot be the particle containing the tail of the worm.
      !it=perm(helium%iw) !tail of the worm
      IF (ik == helium%it) RETURN

      !Second test: detailed balance.
      sumk = 0.d0
      DO ia = 1, natoms
         sumk = sumk+free_density_matrix(helium, ia, l+1, ik, 1)
      END DO
      Pswap1 = sumw/sumk
      xrr = next_random_number(helium%rng_stream_uniform)
      IF (Pswap1 < xrr) THEN
         RETURN ! reject move
      END IF

      ! Once the preliminary tests have been passed, we proceed to perform the move

      ! Compute the potential action before the move:
      U1 = 0.d0
      DO ib = 1, l
         CALL compute_potential_action_link(pint_env, helium, ik, ib, U)
         U1 = U1+U
      END DO

      ! Connect bead 1 of particle ik with the head of the worm:
      dr(:) = helium%work(:, helium%iw, nbeads+1)-helium%work(:, ik, 1)
      CALL helium_pbc(helium, dr)
      rnew(:) = helium%work(:, ik, 1)+dr(:)
      CALL helium_pbc(helium, rnew)
      helium%work(:, ik, 1) = rnew(:)
      !helium%work(:,ik,1)=helium%work(:,helium%iw,nbeads+1)
      CALL staging_move(helium, ik, 1, l)

      ! if no PBC apply confinement:
      IF (.NOT. helium%periodic) THEN
      IF (.NOT. indroplet(helium, ik, 1, l)) THEN
         ! reject move
         helium%work = helium%pos
         ac = 0
         RETURN
      END IF
      END IF

      ! Compute the potential action after the move:
      U2 = 0.d0
      DO ib = 1, l
         CALL compute_potential_action_link(pint_env, helium, ik, ib, U)
         U2 = U2+U
      END DO

      ! Metropolis:
      IF (U2 < U1) THEN !accept move
         helium%pos = helium%work
         !change permutation table:
         ia = perminv(helium, ik)
         perm(ia) = helium%it
         perm(helium%iw) = ik
         helium%iw = ia
         helium%it = perm(helium%iw)
         ac = 1
         RETURN
      END IF
      dU = U2-U1
      dU = MIN(MAX(-200.0_dp, dU), 200.0_dp) ! guard against overflow
      expdU = dexp(-dU)
      Pswap2 = expdU
      xrr = next_random_number(helium%rng_stream_uniform)

      IF (Pswap2 > xrr) THEN !accept move
         helium%pos = helium%work
         !change permutation table:
         ia = perminv(helium, ik)
         perm(ia) = helium%it
         perm(helium%iw) = ik
         helium%iw = ia
         helium%it = perm(helium%iw)
         ac = 1
      ELSE !reject move
         helium%work = helium%pos
         ac = 0
      END IF

   END SUBROUTINE swap

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \param iatom ...
!> \param drmax ...
!> \param ac ...
! **************************************************************************************************
   SUBROUTINE centroid_move(pint_env, helium, iatom, drmax, ac)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: iatom
      REAL(KIND=dp), INTENT(in)                          :: drmax
      INTEGER, INTENT(out)                               :: ac

      INTEGER                                            :: ib, icycle, id, ip, nbeads
      REAL(KIND=dp)                                      :: dr(3), dU, expdU, U, U1, U2, xl, xr

      nbeads = helium%beads

      ! Identify which cycle the atom iatom belongs to:
      icycle = helium_cycle_number(helium, iatom, helium%permutation)

      ! Compute the potential action before the move:
      U1 = 0.d0
      DO ip = 1, helium%atoms
         IF (helium_cycle_number(helium, ip, helium%permutation) == icycle) THEN
            DO ib = 1, nbeads
               CALL compute_potential_action_link(pint_env, helium, iatom, ib, U)
               U1 = U1+U
            END DO
         END IF
      END DO

      DO id = 1, 3
         xl = next_random_number(helium%rng_stream_uniform)
         dr(id) = xl*drmax
      END DO

      ! Perform move:
      DO ip = 1, helium%atoms
         IF (helium_cycle_number(helium, ip, helium%permutation) == icycle) THEN
            CALL centroid_displacement(helium, iatom, dr)
         END IF
      END DO

      DO ip = 1, helium%atoms
         IF (helium_cycle_number(helium, ip, helium%permutation) == icycle) THEN
            ! if no PBC apply confinement:
            IF (.NOT. helium%periodic) THEN
               IF (.NOT. indroplet(helium, ip, 1, helium%beads)) THEN
                  ! reject move
                  helium%work = helium%pos
                  ac = 0
                  RETURN
               END IF
            END IF
         END IF
      END DO

      ! Compute the potential action after the move:
      U2 = 0.d0
      DO ip = 1, helium%atoms
         IF (helium_cycle_number(helium, ip, helium%permutation) == icycle) THEN
            DO ib = 1, nbeads
               CALL compute_potential_action_link(pint_env, helium, iatom, ib, U)
               U2 = U2+U
            END DO
         END IF
      END DO

      ! Metropolis:
      dU = U2-U1
      IF (dU < 0.d0) THEN
         !accept move
         helium%pos = helium%work
         ac = 1
      END IF
      dU = MIN(MAX(-200.0_dp, dU), 200.0_dp) ! guard against overflow
      expdU = dexp(-dU)
      xr = next_random_number(helium%rng_stream_uniform)
      IF (expdU > xr) THEN
         !accept move
         helium%pos = helium%work
         ac = 1
      ELSE
         !reject move:
         helium%work = helium%pos
         ac = 0
      END IF
   END SUBROUTINE centroid_move

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \param iparticle ...
!> \param ibead ...
!> \param U ...
! **************************************************************************************************
   SUBROUTINE compute_potential_action_link(pint_env, helium, iparticle, ibead, U)
      ! Compute the potential action of
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: iparticle, ibead
      REAL(kind=dp), INTENT(out)                         :: U

      INTEGER                                            :: ip
      REAL(kind=dp)                                      :: V1, V2
      REAL(kind=dp), DIMENSION(3)                        :: r, rp

      U = 0.d0
      DO ip = 1, helium%atoms
         IF (ip .NE. iparticle) THEN
            r(:) = helium%work(:, ip, ibead)-helium%work(:, iparticle, ibead)
            rp(:) = helium%work(:, ip, ibead+1)-helium%work(:, iparticle, ibead+1)
            !pair product action for He-He:
            U = U+helium_eval_expansion(helium, r, rp, helium%uij, 1)
         END IF
      END DO

      !Careful with double counting the action in
      !time slices 1 and M+1 when using the primitive approximation:
      IF (helium%solute_present) THEN
         CALL helium_bead_solute_e_f(pint_env, helium, &
                                     iparticle, ibead, helium%work(:, iparticle, ibead), V1)
         IF (ibead < helium%beads) THEN
            ! note that only the bead index and the explicit position
            ! are relevant to the helium_bead_solute_e_f subroutine.
            ! The particle index is not.
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        iparticle, ibead+1, helium%work(:, iparticle, ibead+1), V2)
         ELSE IF (ibead == helium%beads) THEN
            CALL helium_bead_solute_e_f(pint_env, helium, &
                                        helium%permutation(iparticle), 1, &
                                        helium%work(:, iparticle, ibead+1), V2)
         END IF
         !primitive action (symmetrized) for He-solute:
         U = U+0.5_dp*(V1+V2)*helium%tau
      END IF
   END SUBROUTINE compute_potential_action_link

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param ipi ...
!> \param ibi ...
!> \param ipf ...
!> \param ibf ...
!> \param length ...
!> \return ...
! **************************************************************************************************
   FUNCTION free_density_matrix(helium, ipi, ibi, ipf, ibf, length)
      !In the swap move, l=iabs(ibf-ibi) does not hold,
      !therefore it must be passed as an argument
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: ipi, ibi, ipf, ibf
      INTEGER, INTENT(in), OPTIONAL                      :: length
      REAL(KIND=dp)                                      :: free_density_matrix

      INTEGER                                            :: id, l
      REAL(KIND=dp)                                      :: d2, dx, factor, factorpi, vec(3)

      IF (PRESENT(length)) THEN
         l = length
      ELSE
         l = iabs(ibf-ibi)
      END IF
      factor = 1.d0/(2.d0*helium%hb2m*l*helium%tau)
      factorpi = factor/pi
      ! Compute distance squared d2:
      vec(:) = helium%work(:, ipi, ibi)-helium%work(:, ipf, ibf)
      CALL helium_pbc(helium, vec)
      d2 = 0.d0
      DO id = 1, 3
         dx = vec(id)
         d2 = d2+dx*dx
      END DO
      free_density_matrix = factorpi*dsqrt(factorpi)*dexp(-factor*d2)
   END FUNCTION free_density_matrix

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param iatom ...
!> \return ...
! **************************************************************************************************
   FUNCTION perminv(helium, iatom)
      ! Inverse permutation function
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: iatom
      INTEGER                                            :: perminv

      INTEGER                                            :: ip

      !perm(ia)=iparticle ---> perminv(iparticle)=ia
      DO ip = 1, helium%atoms
         IF (helium%permutation(ip) == iatom) THEN
            perminv = ip
            EXIT
         END IF
      END DO
   END FUNCTION perminv

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param iatom ...
!> \param rc ...
! **************************************************************************************************
   SUBROUTINE centroid(helium, iatom, rc)
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: iatom
      REAL(KIND=dp), INTENT(out)                         :: rc(3)

      INTEGER                                            :: ib, id
      REAL(KIND=dp)                                      :: x

      DO id = 1, 3
         x = 0.d0
         DO ib = 1, helium%beads
            x = x+helium%work(id, iatom, ib)
         END DO
         rc(id) = x/REAL(helium%beads, KIND=dp)
      END DO
   END SUBROUTINE centroid

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param iatom ...
!> \param dr ...
! **************************************************************************************************
   SUBROUTINE centroid_displacement(helium, iatom, dr)
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: iatom
      REAL(KIND=dp), INTENT(in)                          :: dr(3)

      INTEGER                                            :: ib

      DO ib = 1, helium%beads+1
         helium%work(:, iatom, ib) = helium%work(:, iatom, ib)+dr(:)
      END DO
   END SUBROUTINE centroid_displacement

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param ip ...
!> \param ibi ...
!> \param l ...
!> \return ...
! **************************************************************************************************
   FUNCTION indroplet(helium, ip, ibi, l)
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      INTEGER, INTENT(in)                                :: ip, ibi, l
      LOGICAL                                            :: indroplet

      INTEGER                                            :: ib, ibb, ipp
      REAL(KIND=dp)                                      :: d2, rcutoff2

      rcutoff2 = helium%droplet_radius**2
      indroplet = .TRUE.
      ! Compute the center of mass of helium%pos (before move)
      IF (helium%solute_present) THEN
         helium%center(:) = 0.d0
      END IF
      ! Check distances at helium%work (after move)
      DO ib = ibi, ibi+l
         ibb = ib
         ipp = ip
         IF (ib > helium%beads+1) THEN
            ibb = ibb-helium%beads
            ipp = helium%permutation(ip)
         END IF
         d2 = distance2(helium%work(:, ipp, ibb), helium%center(:))
         IF (d2 .GE. rcutoff2) THEN
            indroplet = .FALSE.
            RETURN
         END IF
      END DO
   END FUNCTION indroplet

! **************************************************************************************************
!> \brief ...
!> \param v1 ...
!> \param v2 ...
!> \return ...
! **************************************************************************************************
   FUNCTION distance2(v1, v2)
      REAL(KIND=dp), INTENT(in)                          :: v1(3), v2(3)
      REAL(KIND=dp)                                      :: distance2

      INTEGER                                            :: i

      distance2 = 0.d0
      DO i = 1, 3
         distance2 = distance2+(v1(i)-v2(i))**2
      END DO
   END FUNCTION distance2

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \return ...
! **************************************************************************************************
   FUNCTION unit_rand_vec(helium)
      ! Unitary random vector:
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium
      REAL(KIND=dp)                                      :: unit_rand_vec(3)

      REAL(KIND=dp)                                      :: phi, theta, u, v

      u = next_random_number(helium%rng_stream_uniform)
      v = next_random_number(helium%rng_stream_uniform)
      theta = 2.d0*pi*u
      phi = dacos(2.d0*v-1.d0)
      unit_rand_vec(1) = dsin(phi)*dcos(theta)
      unit_rand_vec(2) = dsin(phi)*dsin(theta)
      unit_rand_vec(3) = dcos(phi)
   END FUNCTION unit_rand_vec

! **************************************************************************************************
!> \brief ...
!> \param helium ...
! **************************************************************************************************
   SUBROUTINE center_pos(helium)
      TYPE(helium_solvent_type), INTENT(in), POINTER     :: helium

      INTEGER                                            :: ib, ip

      DO ib = 1, helium%beads
      DO ip = 1, helium%atoms
         helium%pos(:, ip, ib) = helium%pos(:, ip, ib)-helium%center(:)
      END DO
      END DO
      helium%work = helium%pos
      helium%center(:) = 0.d0
   END SUBROUTINE center_pos

END MODULE helium_worm
