!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods that handle helium-solvent and helium-helium interactions
!> \author Lukasz Walewski
!> \date   2009-06-10
! **************************************************************************************************
MODULE helium_interactions

   USE helium_common,                   ONLY: helium_eval_expansion,&
                                              helium_pbc
   USE helium_types,                    ONLY: e_id_interact,&
                                              e_id_kinetic,&
                                              e_id_potential,&
                                              e_id_thermo,&
                                              e_id_total,&
                                              e_id_virial,&
                                              helium_solvent_p_type,&
                                              helium_solvent_type
   USE input_constants,                 ONLY: helium_sampling_worm,&
                                              helium_solute_intpot_ch5,&
                                              helium_solute_intpot_cubetable,&
                                              helium_solute_intpot_mwater,&
                                              helium_solute_intpot_none,&
                                              helium_solute_intpot_seltable
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE physcon,                         ONLY: angstrom,&
                                              kelvin
   USE pint_types,                      ONLY: pint_env_type
   USE splines_types,                   ONLY: spline_data_p_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_interactions'

   PUBLIC :: helium_calc_energy
   PUBLIC :: helium_total_link_action
   PUBLIC :: helium_total_pair_action
   PUBLIC :: helium_total_pair_actionw
   PUBLIC :: helium_total_inter_action
   PUBLIC :: helium_total_inter_actionw
   PUBLIC :: helium_solute_e_f
   PUBLIC :: helium_bead_solute_e_f
   PUBLIC :: helium_intpot_scan

CONTAINS

! ***************************************************************************
!> \brief  Calculate the helium energy (including helium-solute interaction)
!> \param    helium     helium environment
!> \param    pint_env   path integral environment
!> \par History
!>         2009-06 moved I/O out from here [lwalewski]
!> \author hforbert
! **************************************************************************************************
   SUBROUTINE helium_calc_energy(helium, pint_env)
      TYPE(helium_solvent_type), POINTER                 :: helium
      TYPE(pint_env_type), POINTER                       :: pint_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_calc_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: b, bead, c, i, j, n
      INTEGER, DIMENSION(:), POINTER                     :: perm
      LOGICAL                                            :: nperiodic
      REAL(KIND=dp)                                      :: a, cell_size, en, interac, kin, pot, &
                                                            rmax, rmin, vkin
      REAL(KIND=dp), DIMENSION(3)                        :: r, rp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: pos
      TYPE(spline_data_p_type), DIMENSION(:, :), POINTER :: eij

      pos => helium%pos
      perm => helium%permutation
      eij => helium%eij
      cell_size = 0.5_dp*helium%cell_size
      nperiodic = .NOT. helium%periodic
      n = helium%atoms
      b = helium%beads
      en = 0.0_dp
      pot = 0.0_dp
      rmin = 1.0e20_dp
      rmax = 0.0_dp
      DO i = 1, n-1
         DO j = i+1, n
            rp(:) = pos(:, i, 1)-pos(:, j, 1)
            CALL helium_pbc(helium, rp)
            DO bead = 2, b
               a = 0.0_dp
               DO c = 1, 3
                  r(c) = rp(c)
                  a = a+r(c)**2
                  rp(c) = pos(c, i, bead)-pos(c, j, bead)
               END DO
               CALL helium_pbc(helium, rp)
               en = en+helium_eval_expansion(helium, r, rp, eij, 0)
               a = SQRT(a)
               IF (a < rmin) rmin = a
               IF (a > rmax) rmax = a
               IF ((a < cell_size) .OR. nperiodic) THEN
                  pot = pot+helium_vij(a)
               END IF
            END DO
            a = 0.0_dp
            DO c = 1, 3
               r(c) = rp(c)
               a = a+r(c)**2
               rp(c) = pos(c, perm(i), 1)-pos(c, perm(j), 1)
            END DO
            CALL helium_pbc(helium, rp)
            en = en+helium_eval_expansion(helium, r, rp, eij, 0)
            a = SQRT(a)
            IF (a < rmin) rmin = a
            IF (a > rmax) rmax = a
            IF ((a < cell_size) .OR. nperiodic) THEN
               pot = pot+helium_vij(a)
            END IF
         END DO
      END DO
      pot = pot/b
      en = en/b

      ! helium-solute interaction energy (all beads of all particles)
      interac = 0.0_dp
      IF (helium%solute_present) THEN
         CALL helium_solute_e(pint_env, helium, interac)
      END IF
      interac = interac/b

!TODO:
      vkin = 0.0_dp
!   vkin = helium_virial_energy(helium)

      kin = 0.0_dp
      DO i = 1, n
         r(:) = pos(:, i, b)-pos(:, perm(i), 1)
         CALL helium_pbc(helium, r)
         kin = kin+r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
         DO bead = 2, b
            r(:) = pos(:, i, bead-1)-pos(:, i, bead)
            CALL helium_pbc(helium, r)
            kin = kin+r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
         END DO
      END DO
      kin = 1.5_dp*n/helium%tau-0.5*kin/(b*helium%tau**2*helium%hb2m)

! TODO: move printing somwhere else ?
!   print *,"POT = ",(pot/n+helium%e_corr)*kelvin,"K"
!   print *,"INTERAC = ",interac*kelvin,"K"
!   print *,"RMIN= ",rmin*angstrom,"A"
!   print *,"RMAX= ",rmax*angstrom,"A"
!   print *,"EVIRIAL not valid!"
!   print *,"ETHERMO= ",((en+kin)/n+helium%e_corr)*kelvin,"K"
!   print *,"ECORR= ",helium%e_corr*kelvin,"K"
!!   kin = helium_total_action(helium)
!!   print *,"ACTION= ",kin
!   print *,"WINDING#= ",helium_calc_winding(helium)

      helium%energy_inst(e_id_potential) = pot/n+helium%e_corr
      helium%energy_inst(e_id_kinetic) = (en-pot+kin)/n
      helium%energy_inst(e_id_interact) = interac
      helium%energy_inst(e_id_thermo) = (en+kin)/n+helium%e_corr
      helium%energy_inst(e_id_virial) = 0.0_dp !(en+vkin)/n+helium%e_corr
      helium%energy_inst(e_id_total) = (en+vkin)/n+helium%e_corr

      RETURN
   END SUBROUTINE helium_calc_energy

! ***************************************************************************
!> \brief  Computes the total harmonic link action of the helium
!> \param helium ...
!> \return ...
!> \date   2016-05-03
!> \author Felix Uhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION helium_total_link_action(helium) RESULT(linkaction)

      TYPE(helium_solvent_type), POINTER                 :: helium

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_total_link_action', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, ibead
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp), DIMENSION(3)                        :: r

      perm => helium%permutation
      linkaction = 0.0_dp

      ! Harmonic Link action
      ! (r(m-1) - r(m))**2/(4*lambda*tau)
      DO ibead = 1, helium%beads-1
         DO iatom = 1, helium%atoms
            r(:) = helium%pos(:, iatom, ibead)-helium%pos(:, iatom, ibead+1)
            CALL helium_pbc(helium, r)
            linkaction = linkaction+(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
         END DO
      END DO
      DO iatom = 1, helium%atoms
         ! choose last bead connection according to permutation table
         r(:) = helium%pos(:, iatom, helium%beads)-helium%pos(:, perm(iatom), 1)
         CALL helium_pbc(helium, r)
         linkaction = linkaction+(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
      END DO
      linkaction = linkaction/(2.0_dp*helium%tau*helium%hb2m)

      RETURN

   END FUNCTION helium_total_link_action

! ***************************************************************************
!> \brief  Computes the total pair action of the helium
!> \param helium ...
!> \return ...
!> \date   2016-05-03
!> \author Felix Uhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION helium_total_pair_action(helium) RESULT(pairaction)

      TYPE(helium_solvent_type), POINTER                 :: helium

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_total_pair_action', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, ibead, jatom, opatom, patom
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp), DIMENSION(3)                        :: r, rp

      perm => helium%permutation
      pairaction = 0.0_dp

      ! He-He pair action
      DO ibead = 1, helium%beads-1
         DO iatom = 1, helium%atoms-1
            DO jatom = iatom+1, helium%atoms
!WRITE(*,*) "@FUg total     ", iatom ,"--",jatom,":", ibead
!WRITE(*,*) "@FUg total     ", iatom ,"--",jatom,":", ibead + 1
               r(:) = helium%pos(:, iatom, ibead)-helium%pos(:, jatom, ibead)
               rp(:) = helium%pos(:, iatom, ibead+1)-helium%pos(:, jatom, ibead+1)
!WRITE(*,*) "@FUg", helium_eval_expansion(helium, r, rp, helium%uij, 1)
               pairaction = pairaction+helium_eval_expansion(helium, r, rp, helium%uij, 1)
            END DO
         END DO
      END DO
      !Ensure right permutation for pair action of last and first beads.
      DO iatom = 1, helium%atoms-1
         DO jatom = iatom+1, helium%atoms
!WRITE(*,*) "@FUg total     ", iatom ,"--",jatom,":", helium%beads
!WRITE(*,*) "@FUg", helium%pos(:, iatom, helium%beads)
!WRITE(*,*) "@FUg", helium%pos(:, jatom, helium%beads)
            r(:) = helium%pos(:, iatom, helium%beads)-helium%pos(:, jatom, helium%beads)
!WRITE(*,*) "@FUg", r
!WRITE(*,*) "@FUg total     ", perm(iatom) ,"--",perm(jatom),":", 1
!WRITE(*,*) "@FUg", helium%pos(:, perm(iatom), 1)
!WRITE(*,*) "@FUg", helium%pos(:, perm(jatom), 1)
            rp(:) = helium%pos(:, perm(iatom), 1)-helium%pos(:, perm(jatom), 1)
!WRITE(*,*) "@FUg", rp
!WRITE(*,*) "@FUg", helium_eval_expansion(helium, r, rp, helium%uij, 1)
            pairaction = pairaction+helium_eval_expansion(helium, r, rp, helium%uij, 1)
         END DO
      END DO

      ! correct for open worm configurations
      IF (.NOT. helium%worm_is_closed) THEN
!WRITE(*,*) "@FUg total correcting action for open worm"
         ! special treatment if double bead is first bead
         iatom = helium%worm_atom_idx
         IF (helium%worm_bead_idx == 1) THEN
            ! patom is the atom in front of the lone head bead
            patom = helium%iperm(iatom)
            ! go through all atoms
            DO jatom = 1, helium%atoms
               IF (jatom == helium%worm_atom_idx) CYCLE
               opatom = helium%iperm(jatom)
!WRITE(*,*) "@FUg total rm  ", iatom ,"--",jatom,":", 1, helium%beads
               ! substract pair action for closed link
               r(:) = helium%pos(:, iatom, 1)-helium%pos(:, jatom, 1)
               rp(:) = helium%pos(:, patom, helium%beads)-helium%pos(:, opatom, helium%beads)
               pairaction = pairaction-helium_eval_expansion(helium, r, rp, helium%uij, 1)
!WRITE(*,*) "@FUg ", r(:)
!WRITE(*,*) "@FUg ", rp(:)
!WRITE(*,*) "@FUg - ", helium_eval_expansion(helium, r, rp, helium%uij, 1)
               ! and add corrected extra link
               ! rp stays the same
!WRITE(*,*) "@FUg total add ", iatom ,"--",jatom,":", "        xtra", helium%beads
               r(:) = helium%worm_xtra_bead(:)-helium%pos(:, jatom, 1)
               pairaction = pairaction+helium_eval_expansion(helium, r, rp, helium%uij, 1)
!WRITE(*,*) "@FUg ", r(:)
!WRITE(*,*) "@FUg + ", helium_eval_expansion(helium, r, rp, helium%uij, 1)
            END DO
         ELSE
            ! bead stays constant
            ibead = helium%worm_bead_idx
            ! go through all atoms
            DO jatom = 1, helium%atoms
               IF (jatom == helium%worm_atom_idx) CYCLE
!WRITE(*,*) "@FUg total rm  ", iatom ,"--",jatom,":", ibead, ibead-1
               ! substract pair action for closed link
               r(:) = helium%pos(:, iatom, ibead)-helium%pos(:, jatom, ibead)
               rp(:) = helium%pos(:, iatom, ibead-1)-helium%pos(:, jatom, ibead-1)
               pairaction = pairaction-helium_eval_expansion(helium, r, rp, helium%uij, 1)
!WRITE(*,*) "@FUg - ", helium_eval_expansion(helium, r, rp, helium%uij, 1)
               ! and add corrected extra link
               ! rp stays the same
!WRITE(*,*) "@FUg total add ", iatom,"--",jatom,":", "        xtra", ibead-1
               r(:) = helium%worm_xtra_bead(:)-helium%pos(:, jatom, ibead)
               pairaction = pairaction+helium_eval_expansion(helium, r, rp, helium%uij, 1)
!WRITE(*,*) "@FUg + ", helium_eval_expansion(helium, r, rp, helium%uij, 1)
            END DO
         END IF
      END IF

      RETURN
   END FUNCTION helium_total_pair_action

! ***************************************************************************
!> \brief  Computes the total pair action of the helium
!> \param helium ...
!> \return ...
!> \date   2016-05-03
!> \author Felix Uhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION helium_total_pair_actionw(helium) RESULT(pairaction)

      TYPE(helium_solvent_type), POINTER                 :: helium

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_total_pair_actionw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, ibead, jatom, opatom, patom
      INTEGER, DIMENSION(:), POINTER                     :: perm
      REAL(KIND=dp), DIMENSION(3)                        :: r, rp

      perm => helium%permutation
      pairaction = 0.0_dp

      ! He-He pair action
      DO ibead = 1, helium%beads-1
         DO iatom = 1, helium%atoms-1
            DO jatom = iatom+1, helium%atoms
!WRITE(*,*) "@FUg total     ", iatom ,"--",jatom,":", ibead
!WRITE(*,*) "@FUg total     ", iatom ,"--",jatom,":", ibead + 1
               r(:) = helium%work(:, iatom, ibead)-helium%work(:, jatom, ibead)
               rp(:) = helium%work(:, iatom, ibead+1)-helium%work(:, jatom, ibead+1)
!WRITE(*,*) "@FUg", helium_eval_expansion(helium, r, rp, helium%uij, 1)
               pairaction = pairaction+helium_eval_expansion(helium, r, rp, helium%uij, 1)
            END DO
         END DO
      END DO
      !Ensure right permutation for pair action of last and first beads.
      DO iatom = 1, helium%atoms-1
         DO jatom = iatom+1, helium%atoms
!WRITE(*,*) "@FUg total     ", iatom ,"--",jatom,":", helium%beads
!WRITE(*,*) "@FUg", helium%work(:, iatom, helium%beads)
!WRITE(*,*) "@FUg", helium%work(:, jatom, helium%beads)
            r(:) = helium%work(:, iatom, helium%beads)-helium%work(:, jatom, helium%beads)
!WRITE(*,*) "@FUg", r
!WRITE(*,*) "@FUg total     ", perm(iatom) ,"--",perm(jatom),":", 1
!WRITE(*,*) "@FUg", helium%work(:, perm(iatom), 1)
!WRITE(*,*) "@FUg", helium%work(:, perm(jatom), 1)
            rp(:) = helium%work(:, perm(iatom), 1)-helium%work(:, perm(jatom), 1)
!WRITE(*,*) "@FUg", rp
!WRITE(*,*) "@FUg", helium_eval_expansion(helium, r, rp, helium%uij, 1)
            pairaction = pairaction+helium_eval_expansion(helium, r, rp, helium%uij, 1)
         END DO
      END DO

      ! correct for open worm configurations
      IF (.NOT. helium%worm_is_closed) THEN
!WRITE(*,*) "@FUg total correcting actionw for open worm"
         ! special treatment if double bead is first bead
         iatom = helium%worm_atom_idx_work
         IF (helium%worm_bead_idx_work == 1) THEN
            ! patom is the atom in front of the lone head bead
            patom = helium%iperm(iatom)
            ! go through all atoms
            DO jatom = 1, helium%atoms
               IF (jatom == helium%worm_atom_idx_work) CYCLE
               opatom = helium%iperm(jatom)
!WRITE(*,*) "@FUg total rm  ", iatom ,"--",jatom,":", 1, helium%beads
               ! substract pair action for closed link
               r(:) = helium%work(:, iatom, 1)-helium%work(:, jatom, 1)
               rp(:) = helium%work(:, patom, helium%beads)-helium%work(:, opatom, helium%beads)
               pairaction = pairaction-helium_eval_expansion(helium, r, rp, helium%uij, 1)
!WRITE(*,*) "@FUg ", r(:)
!WRITE(*,*) "@FUg ", rp(:)
!WRITE(*,*) "@FUg - ", helium_eval_expansion(helium, r, rp, helium%uij, 1)
               ! and add corrected extra link
               ! rp stays the same
!WRITE(*,*) "@FUg total add ", iatom ,"--",jatom,":", "        xtra", helium%beads
               r(:) = helium%worm_xtra_bead_work(:)-helium%work(:, jatom, 1)
               pairaction = pairaction+helium_eval_expansion(helium, r, rp, helium%uij, 1)
!WRITE(*,*) "@FUg ", r(:)
!WRITE(*,*) "@FUg + ", helium_eval_expansion(helium, r, rp, helium%uij, 1)
            END DO
         ELSE
            ! bead stays constant
            ibead = helium%worm_bead_idx_work
            ! go through all atoms
            DO jatom = 1, helium%atoms
               IF (jatom == helium%worm_atom_idx_work) CYCLE
!WRITE(*,*) "@FUg total rm  ", iatom ,"--",jatom,":", ibead, ibead-1
               ! substract pair action for closed link
               r(:) = helium%work(:, iatom, ibead)-helium%work(:, jatom, ibead)
               rp(:) = helium%work(:, iatom, ibead-1)-helium%work(:, jatom, ibead-1)
               pairaction = pairaction-helium_eval_expansion(helium, r, rp, helium%uij, 1)
!WRITE(*,*) "@FUg - ", helium_eval_expansion(helium, r, rp, helium%uij, 1)
               ! and add corrected extra link
               ! rp stays the same
!WRITE(*,*) "@FUg total add ", iatom ,"--",jatom,":", "        xtra", ibead-1
               r(:) = helium%worm_xtra_bead_work(:)-helium%work(:, jatom, ibead)
               pairaction = pairaction+helium_eval_expansion(helium, r, rp, helium%uij, 1)
!WRITE(*,*) "@FUg + ",  helium_eval_expansion(helium, r, rp, helium%uij, 1)
            END DO
         END IF
      END IF

      RETURN
   END FUNCTION helium_total_pair_actionw

! ***************************************************************************
!> \brief  Computes the total interaction of the helium with the solute
!> \param pint_env ...
!> \param helium ...
!> \return ...
!> \date   2016-05-03
!> \author Felix Uhl
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION helium_total_inter_action(pint_env, helium) RESULT(interaction)

      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), POINTER                 :: helium

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_total_inter_action', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, ibead
      REAL(KIND=dp)                                      :: e

      interaction = 0.0_dp

      ! InterAction with solute
      IF (helium%solute_present) THEN
         DO ibead = 1, helium%beads
            DO iatom = 1, helium%atoms

               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           iatom, ibead, helium%pos(:, iatom, ibead), e)
!WRITE(*,*) "@FUg +1.0 ",  ibead, iatom, e
               interaction = interaction+e
            END DO
         END DO
         IF (helium%sampling_method == helium_sampling_worm) THEN
            IF (.NOT. helium%worm_is_closed) THEN
               ! subtract half of tail bead interaction again
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           helium%worm_atom_idx, helium%worm_bead_idx, &
                                           helium%pos(:, helium%worm_atom_idx, helium%worm_bead_idx), e)
!WRITE(*,*) "@FUg -0.5 ",  helium%worm_atom_idx, helium%worm_bead_idx, e
               interaction = interaction-0.5_dp*e
               ! add half of head bead interaction
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           helium%worm_atom_idx, helium%worm_bead_idx, &
                                           helium%worm_xtra_bead, e)
!WRITE(*,*) "@FUg +0.5 ",  "head ", e
               interaction = interaction+0.5_dp*e
            END IF
         END IF
      END IF

      interaction = interaction*helium%tau

      RETURN
   END FUNCTION helium_total_inter_action

! **************************************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium ...
!> \return ...
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION helium_total_inter_actionw(pint_env, helium) RESULT(interaction)

      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), POINTER                 :: helium

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_total_inter_actionw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, ibead
      REAL(KIND=dp)                                      :: e

      interaction = 0.0_dp

      ! InterAction with solute
      IF (helium%solute_present) THEN
         DO ibead = 1, helium%beads
            DO iatom = 1, helium%atoms
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           iatom, ibead, helium%work(:, iatom, ibead), e)
!WRITE(*,*) "@FUg +1.0 ",  ibead, iatom, e
               interaction = interaction+e
            END DO
         END DO
         IF (helium%sampling_method == helium_sampling_worm) THEN
            IF (.NOT. helium%worm_is_closed) THEN
               ! subtract half of tail bead interaction again
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           helium%worm_atom_idx_work, helium%worm_bead_idx_work, &
                                           helium%work(:, helium%worm_atom_idx_work, helium%worm_bead_idx_work), e)
!WRITE(*,*) "@FUg -0.5 ",  helium%worm_atom_idx, helium%worm_bead_idx, e
               interaction = interaction-0.5_dp*e
               ! add half of head bead interaction
               CALL helium_bead_solute_e_f(pint_env, helium, &
                                           helium%worm_atom_idx_work, helium%worm_bead_idx_work, &
                                           helium%worm_xtra_bead_work, e)
!WRITE(*,*) "@FUg +0.5 ",  "head ", e
               interaction = interaction+0.5_dp*e
            END IF
         END IF
      END IF

      interaction = interaction*helium%tau

      RETURN
   END FUNCTION helium_total_inter_actionw

! ***************************************************************************
!> \brief Calculate general helium-solute interaction energy (and forces)
!>        between one helium bead and the corresponding solute time slice.
!> \param pint_env           path integral environment
!> \param helium ...
!> \param helium_part_index  helium particle index
!> \param helium_slice_index helium time slice index
!> \param helium_r_opt       explicit helium bead coordinates (optional)
!> \param energy             calculated energy
!> \param force              calculated force (if requested)
!> \par History
!>         2019-09 Added multiple-time striding in imag. time [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_bead_solute_e_f(pint_env, helium, helium_part_index, &
                                     helium_slice_index, helium_r_opt, energy, force)

      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), POINTER                 :: helium
      INTEGER, INTENT(IN)                                :: helium_part_index, helium_slice_index
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: helium_r_opt
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL, POINTER                               :: force

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_bead_solute_e_f', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: hbeads, hi, qi, stride
      REAL(KIND=dp), DIMENSION(3)                        :: helium_r
      REAL(KIND=dp), DIMENSION(:), POINTER               :: my_force

      hbeads = helium%beads
      ! helium bead index that is invariant wrt the rotations
      hi = MOD(helium_slice_index-1+hbeads+helium%relrot, hbeads)+1
      ! solute bead index that belongs to hi helium index
      qi = ((hi-1)*pint_env%p)/hbeads+1

      ! coordinates of the helium bead
      IF (PRESENT(helium_r_opt)) THEN
         helium_r(:) = helium_r_opt(:)
      ELSE
         helium_r(:) = helium%pos(:, helium_part_index, helium_slice_index)
      END IF

      SELECT CASE (helium%solute_interaction)

      CASE (helium_solute_intpot_mwater)
         IF (PRESENT(force)) THEN
            force(:, :) = 0.0_dp
            my_force => force(qi, :)
            CALL helium_intpot_model_water( &
               pint_env%x(qi, :), &
               helium, &
               helium_r, &
               energy, &
               my_force &
               )
         ELSE
            CALL helium_intpot_model_water( &
               pint_env%x(qi, :), &
               helium, &
               helium_r, &
               energy &
               )
         END IF

      CASE (helium_solute_intpot_ch5)
         IF (PRESENT(force)) THEN
            force(:, :) = 0.0_dp
            my_force => force(qi, :)
            CALL helium_intpot_ch5_kuchenbaecker_e_f( &
               pint_env%x(qi, :), &
               helium, &
               helium_r, &
               energy, &
               my_force &
               )
         ELSE
            CALL helium_intpot_ch5_kuchenbaecker_e_f( &
               pint_env%x(qi, :), &
               helium, &
               helium_r, &
               energy &
               )

         END IF

      CASE (helium_solute_intpot_seltable)
         IF (PRESENT(force)) THEN
            ! not implemented for this potential, thus
            force(:, :) = 0.0_dp
         END IF
         CALL helium_intpot_seltable_e( &
            pint_env%x(qi, :), &
            helium, &
            helium_r, &
            energy &
            )

      CASE (helium_solute_intpot_cubetable)
         IF (PRESENT(force)) THEN
            ! not implemented for this potential, thus
            force(:, :) = 0.0_dp
         END IF
         CALL helium_intpot_cubetable( &
            pint_env%x(qi, :), &
            helium, &
            helium_r, &
            energy &
            )

      CASE (helium_solute_intpot_none)
         energy = 0.0_dp
         IF (PRESENT(force)) THEN
            force(:, :) = 0.0_dp
         END IF

      CASE DEFAULT

      END SELECT

      ! Account for Imaginary time striding in forces:
      IF (PRESENT(force)) THEN
         IF (hbeads < pint_env%p) THEN
            stride = pint_env%p/hbeads
            force = force*REAL(stride, dp)
         END IF
      END IF

      RETURN
   END SUBROUTINE helium_bead_solute_e_f

! ***************************************************************************
!> \brief Calculate total helium-solute interaction energy and forces.
!> \param   pint_env   path integral environment
!> \param helium ...
!> \param   energy     calculated interaction energy
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_solute_e_f(pint_env, helium, energy)

      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), POINTER                 :: helium
      REAL(KIND=dp), INTENT(OUT)                         :: energy

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_solute_e_f', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, jb, jc
      REAL(KIND=dp)                                      :: my_energy
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: force

      NULLIFY (force)
      force => helium%force_inst

      energy = 0.0_dp
      force(:, :) = 0.0_dp

      ! calculate the total interaction energy and gradients between the
      ! solute and the helium, sum over all beads of all He particles
      DO ia = 1, helium%atoms
         DO ib = 1, helium%beads
            CALL helium_bead_solute_e_f(pint_env, helium, ia, ib, &
                                        energy=my_energy, force=helium%rtmp_p_ndim_2d)
            energy = energy+my_energy
            DO jb = 1, pint_env%p
               DO jc = 1, pint_env%ndim
                  force(jb, jc) = force(jb, jc)+helium%rtmp_p_ndim_2d(jb, jc)
               END DO
            END DO
         END DO
      END DO

      RETURN
   END SUBROUTINE helium_solute_e_f

! ***************************************************************************
!> \brief Calculate total helium-solute interaction energy.
!> \param   pint_env   path integral environment
!> \param helium ...
!> \param   energy     calculated interaction energy
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_solute_e(pint_env, helium, energy)

      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_type), POINTER                 :: helium
      REAL(KIND=dp), INTENT(OUT)                         :: energy

      INTEGER                                            :: ia, ib
      REAL(KIND=dp)                                      :: my_energy

      energy = 0.0_dp

      DO ia = 1, helium%atoms
         DO ib = 1, helium%beads
            CALL helium_bead_solute_e_f(pint_env, helium, ia, ib, energy=my_energy)
            energy = energy+my_energy
         END DO
      END DO

      RETURN
   END SUBROUTINE helium_solute_e

   ! ***************************************************************************
   !> \brief Calculate helium-solute interaction energy and forces between one helium bead
   !>        and the corresponding solute time slice asuming CH5+ solute.
   !> \param solute_x  solute positions ARR(3*NATOMS)
   !> \param helium    only needed for helium_pbc call at the moment
   !> \param helium_x  helium bead position ARR(3)
   !> \param energy    calculated interaction energy
   !> \note  Based on original subroutine from Harald's standalone program
   !>        with induction.
   !> \author Felix Uhl
   ! ***************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param solute_x ...
!> \param helium ...
!> \param helium_x ...
!> \param energy ...
!> \param force ...
! **************************************************************************************************
   SUBROUTINE helium_intpot_ch5_kuchenbaecker_e_f(solute_x, helium, helium_x, energy, force)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: solute_x
      TYPE(helium_solvent_type), POINTER                 :: helium
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: helium_x
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), &
         OPTIONAL, POINTER                               :: force

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_intpot_ch5_kuchenbaecker_e_f', &
         routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER :: alpha2h = 10.0000883993_dp, alphac = 537.2365703187_dp, &
         alphah = 11.1227066949_dp, alphawx = 0.1393170339_dp, beta2h = 4.9999902914_dp, &
         betac = 2.6414399416_dp, betah = 3.0453141136_dp, betawx = 1.8943211280_dp, &
         Csixc = 25.1486400697_dp, Csixh = -20.7731766063_dp, Csixwx = 0.9756843988_dp, &
         Csixx = 6.1486734828_dp, deltac = 2.1660761212_dp, deltadc = 5.3180112328_dp, &
         deltadh = 1.0000000000_dp, deltadx = 8.0000000000_dp, deltah = 4.7880909104_dp, &
         deltaqc = 1.0000000001_dp, deltaqh = 1.0000000000_dp, deltaqx = 1.6120189633_dp, &
         deltawx = 3.8778221297_dp
      REAL(KIND=dp), PARAMETER :: deltax = 0.0062370139_dp, dp_pol_He = 1.3824_dp, &
         findu = 2.0050419163_dp, fvdw = 0.8134562093_dp, fwvdw = 1.6108570993_dp, &
         qc = 4.6431264647_dp, qh = -0.6771456067_dp, qp_pol_He = 2.538_dp*3.0_dp/2.0_dp, &
         qx = -0.0514796862_dp

      INTEGER                                            :: i, ig, j, k, m, n, num_carbon, &
                                                            num_hydrogen
      REAL(KIND=dp) :: ddenom, ddenomsq, Eindd, Eindq, Evdw, expbetar, expbetarsq, qdenom, &
         qdenomsq, Qxx, Qxy, Qxz, Qyy, Qyz, Qzz, r, rcube, rfif, rfour, rsqr, vdenom, vdenomsq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: F, rvekt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dQxx, dQxy, dQxz, dQyy, dQyz, dQzz, &
                                                            gradient, indusites, solute, vdwsites, &
                                                            vdwwsites
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Fprime

! Induction Parameters
! Helium polarizabilities
! Induction charges
! dipole field
! quadrupole field
! extrasite position
!VDW Parameters
! pre exponential factor
! exponential factor
! C6-VDW factor
! denominator shift
! extrasite position
! angular extrasite position
! TODO avoid hard wired array size!
!Field
!Dipole Field gradient (term index, diff index, atom index)
!Quadrupole (Fieldgradient)
!Quadrupole gradients
! r, r**2, r**3, r**4, r**5

      num_carbon = 1
      num_hydrogen = 5

      energy = 0.0_dp

      ! Allocate all arrays for energy and optionally the gradient ones
      ALLOCATE (F(3))
      F(:) = 0.0_dp
      ALLOCATE (indusites(3, 5))
      indusites(:, :) = 0.0_dp
      ALLOCATE (vdwsites(3, 5))
      vdwsites(:, :) = 0.0_dp
      ALLOCATE (vdwwsites(3, 10))
      vdwwsites(:, :) = 0.0_dp
      ALLOCATE (solute(3, 6))
      ALLOCATE (rvekt(3))

      Qxx = 0.0_dp
      Qyy = 0.0_dp
      Qzz = 0.0_dp
      Qxy = 0.0_dp
      Qxz = 0.0_dp
      Qyz = 0.0_dp

      Evdw = 0.0_dp
      Eindd = 0.0_dp
      Eindq = 0.0_dp

      IF (PRESENT(force)) THEN
         ALLOCATE (Fprime(3, 3, 6))
         Fprime(:, :, :) = 0.0_dp
         ALLOCATE (dQxx(3, 6))
         dQxx(:, :) = 0.0_dp
         ALLOCATE (dQyy(3, 6))
         dQyy(:, :) = 0.0_dp
         ALLOCATE (dQzz(3, 6))
         dQzz(:, :) = 0.0_dp
         ALLOCATE (dQxy(3, 6))
         dQxy(:, :) = 0.0_dp
         ALLOCATE (dQxz(3, 6))
         dQxz(:, :) = 0.0_dp
         ALLOCATE (dQyz(3, 6))
         dQyz(:, :) = 0.0_dp
         ALLOCATE (gradient(3, 6))
         gradient(:, :) = 0.0_dp
      END IF

      ig = 1
      DO j = 1, helium%solute_atoms
         IF ("C " == helium%solute_element(j)) THEN
            DO i = 1, 3
               solute(i, ig) = solute_x(3*(j-1)+i)
            END DO
            ig = ig+1
         END IF
      END DO
      DO j = 1, helium%solute_atoms
         IF ("H " == helium%solute_element(j)) THEN
            DO i = 1, 3
               solute(i, ig) = solute_x(3*(j-1)+i)
            END DO
            ig = ig+1
         END IF
      END DO

      ! calculate induction sites
      DO j = 2, 6
         DO i = 1, 3
            indusites(i, j-1) = solute(i, 1)+findu*(solute(i, j)-solute(i, 1))
         END DO
      END DO

      ! calculate vdw sites
      DO j = 2, 6
         DO i = 1, 3
            vdwsites(i, j-1) = solute(i, 1)+fvdw*(solute(i, j)-solute(i, 1))
         END DO
      END DO

      ! calculate vdw angular sites
      ! calculate angular vdw sites
      ig = 0
      DO k = 2, 5
         DO j = k+1, 6
            ig = ig+1
            DO i = 1, 3
               vdwwsites(i, ig) = solute(i, 1)+0.5_dp*fwvdw*(solute(i, k)+solute(i, j)-2.0_dp*solute(i, 1))
            END DO
         END DO
      END DO

      ! CARBON
      DO j = 1, 3
         rvekt(j) = solute(j, 1)-helium_x(j)
      END DO
      CALL helium_pbc(helium, rvekt)
      rsqr = rvekt(1)*rvekt(1)+rvekt(2)*rvekt(2)+rvekt(3)*rvekt(3)
      r = SQRT(rsqr)
      rcube = r*rsqr
      rfour = rsqr*rsqr
      rfif = rsqr*rcube
      ddenom = rcube+deltadc**3
      qdenom = rfif+deltaqc**5
      vdenom = rcube*rcube+deltac**6
      expbetar = EXP(-betac*r)
      expbetarsq = EXP(-1.5_dp*rsqr)
      !vdwenergy
      Evdw = Evdw+alphac*expbetar
      Evdw = Evdw+Csixc/vdenom
      Evdw = Evdw+1.0d6*expbetarsq
      !dipole induction terms
      DO j = 1, 3
         F(j) = F(j)+qc*rvekt(j)/ddenom
      END DO
      !quadrupole induction terms
      Qxx = Qxx+qc*(3.0_dp*rvekt(1)*rvekt(1)-rsqr)/qdenom
      Qyy = Qyy+qc*(3.0_dp*rvekt(2)*rvekt(2)-rsqr)/qdenom
      Qzz = Qzz+qc*(3.0_dp*rvekt(3)*rvekt(3)-rsqr)/qdenom
      Qxy = Qxy+qc*(3.0_dp*rvekt(1)*rvekt(2))/qdenom
      Qxz = Qxz+qc*(3.0_dp*rvekt(1)*rvekt(3))/qdenom
      Qyz = Qyz+qc*(3.0_dp*rvekt(2)*rvekt(3))/qdenom
      IF (PRESENT(force)) THEN
         vdenomsq = vdenom*vdenom
         DO k = 1, 3
            gradient(k, 1) = gradient(k, 1)-alphac*betac*expbetar*rvekt(k)/r
            gradient(k, 1) = gradient(k, 1)-6.0_dp*Csixc*rfour*rvekt(k)/vdenomsq
            gradient(k, 1) = gradient(k, 1)-3.0_dp*1.0d6*expbetarsq*rvekt(k)
         END DO
         !dipole gradients
         ddenomsq = ddenom*ddenom
         DO i = 1, 3
            Fprime(i, i, 1) = Fprime(i, i, 1)+qc*(ddenom-3.0_dp*rvekt(i)**2*r)/ddenomsq
         END DO
         !off diagnoal elements
         DO i = 1, 2
            DO j = i+1, 3
               Fprime(j, i, 1) = Fprime(j, i, 1)+qc*(-3.0_dp*rvekt(i)*rvekt(j)*r)/ddenomsq
            END DO
         END DO
         !quadrupole gradients
         qdenomsq = qdenom*qdenom
         dQxx(1, 1) = dQxx(1, 1)+qc*(4.0_dp*rvekt(1)*qdenom- &
                                     15.0_dp*rvekt(1)**3*rcube+ &
                                     5.0_dp*rvekt(1)*rfif)/qdenomsq
         dQxx(2, 1) = dQxx(2, 1)+qc*(-2.0_dp*rvekt(2)*qdenom- &
                                     15.0_dp*rvekt(2)*rvekt(1)**2*rcube+ &
                                     5.0_dp*rvekt(2)*rfif)/qdenomsq
         dQxx(3, 1) = dQxx(3, 1)+qc*(-2.0_dp*rvekt(3)*qdenom- &
                                     15.0_dp*rvekt(3)*rvekt(1)**2*rcube+ &
                                     5.0_dp*rvekt(3)*rfif)/qdenomsq

         dQyy(1, 1) = dQyy(1, 1)+qc*(-2.0_dp*rvekt(1)*qdenom- &
                                     15.0_dp*rvekt(1)*rvekt(2)**2*rcube+ &
                                     5.0_dp*rvekt(1)*rfif)/qdenomsq
         dQyy(2, 1) = dQyy(2, 1)+qc*(4.0_dp*rvekt(2)*qdenom- &
                                     15.0_dp*rvekt(2)**3*rcube+ &
                                     5.0_dp*rvekt(2)*rfif)/qdenomsq
         dQyy(3, 1) = dQyy(3, 1)+qc*(-2.0_dp*rvekt(3)*qdenom- &
                                     15.0_dp*rvekt(3)*rvekt(2)**2*rcube+ &
                                     5.0_dp*rvekt(3)*rfif)/qdenomsq

         dQzz(1, 1) = dQzz(1, 1)+qc*(-2.0_dp*rvekt(1)*qdenom- &
                                     15.0_dp*rvekt(1)*rvekt(3)**2*rcube+ &
                                     5.0_dp*rvekt(1)*rfif)/qdenomsq
         dQzz(2, 1) = dQzz(2, 1)+qc*(-2.0_dp*rvekt(2)*qdenom- &
                                     15.0_dp*rvekt(2)*rvekt(3)**2*rcube+ &
                                     5.0_dp*rvekt(2)*rfif)/qdenomsq
         dQzz(3, 1) = dQzz(3, 1)+qc*(4.0_dp*rvekt(3)*qdenom- &
                                     15.0_dp*rvekt(3)**3*rcube+ &
                                     5.0_dp*rvekt(3)*rfif)/qdenomsq

         dQxy(1, 1) = dQxy(1, 1)+qc*(3.0_dp*rvekt(2)*qdenom- &
                                     15.0_dp*rvekt(2)*rvekt(1)**2*rcube)/qdenomsq
         dQxy(2, 1) = dQxy(2, 1)+qc*(3.0_dp*rvekt(1)*qdenom- &
                                     15.0_dp*rvekt(1)*rvekt(2)**2*rcube)/qdenomsq
         dQxy(3, 1) = dQxy(3, 1)-qc*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq

         dQxz(1, 1) = dQxz(1, 1)+qc*(3.0_dp*rvekt(3)*qdenom- &
                                     15.0_dp*rvekt(3)*rvekt(1)**2*rcube)/qdenomsq
         dQxz(2, 1) = dQxz(2, 1)-qc*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq
         dQxz(3, 1) = dQxz(3, 1)+qc*(3.0_dp*rvekt(1)*qdenom- &
                                     15.0_dp*rvekt(1)*rvekt(3)**2*rcube)/qdenomsq

         dQyz(1, 1) = dQyz(1, 1)-qc*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq
         dQyz(2, 1) = dQyz(2, 1)+qc*(3.0_dp*rvekt(3)*qdenom- &
                                     15.0_dp*rvekt(3)*rvekt(2)**2*rcube)/qdenomsq
         dQyz(3, 1) = dQyz(3, 1)+qc*(3.0_dp*rvekt(2)*qdenom- &
                                     15.0_dp*rvekt(2)*rvekt(3)**2*rcube)/qdenomsq
      END IF

      ! HYDROGEN
      DO i = 2, 6
         DO j = 1, 3
            rvekt(j) = solute(j, i)-helium_x(j)
         END DO
         CALL helium_pbc(helium, rvekt)
         rsqr = rvekt(1)*rvekt(1)+rvekt(2)*rvekt(2)+rvekt(3)*rvekt(3)
         r = SQRT(rsqr)
         rcube = r*rsqr
         rfour = rsqr*rsqr
         rfif = rsqr*rcube
         ddenom = rcube+deltadh**3
         qdenom = rfif+deltaqh**5
         vdenom = rcube*rcube+deltah**6
         expbetar = EXP(-betah*r)
         expbetarsq = EXP(-beta2h*rsqr)
         !vdw energy
         Evdw = Evdw+alphah*expbetar
         Evdw = Evdw+Csixh/vdenom
         Evdw = Evdw+alpha2h*expbetarsq
         !dipole induction terms
         DO j = 1, 3
            F(j) = F(j)+qh*rvekt(j)/ddenom
         END DO
         !quadrupole induction terms
         Qxx = Qxx+qh*(3.0_dp*rvekt(1)*rvekt(1)-rsqr)/qdenom
         Qyy = Qyy+qh*(3.0_dp*rvekt(2)*rvekt(2)-rsqr)/qdenom
         Qzz = Qzz+qh*(3.0_dp*rvekt(3)*rvekt(3)-rsqr)/qdenom
         Qxy = Qxy+qh*(3.0_dp*rvekt(1)*rvekt(2))/qdenom
         Qxz = Qxz+qh*(3.0_dp*rvekt(1)*rvekt(3))/qdenom
         Qyz = Qyz+qh*(3.0_dp*rvekt(2)*rvekt(3))/qdenom
         IF (PRESENT(force)) THEN
            vdenomsq = vdenom*vdenom
            DO k = 1, 3
               gradient(k, i) = gradient(k, i)-alphah*betah*expbetar*rvekt(k)/r
               gradient(k, i) = gradient(k, i)-6.0_dp*Csixh*rfour*rvekt(k)/vdenomsq
               gradient(k, i) = gradient(k, i)-2.0_dp*alpha2h*beta2h*expbetarsq*rvekt(k)/r
            END DO
            ! dipole field gradient with respect to ith hydrogen coordinate
            ddenomsq = ddenom*ddenom
            DO j = 1, 3
               Fprime(j, j, i) = Fprime(j, j, i)+qh*(ddenom-3.0_dp*rvekt(j)**2*r)/ddenomsq
            END DO
            !off diagnoal elements
            DO j = 1, 2
               DO k = j+1, 3
                  Fprime(k, j, i) = Fprime(k, j, i)+qh*(-3.0_dp*rvekt(j)*rvekt(k)*r)/ddenomsq
               END DO
            END DO
            ! Quadrupole field gradient
            qdenomsq = qdenom*qdenom
            dQxx(1, i) = dQxx(1, i)+qh*(4.0_dp*rvekt(1)*qdenom- &
                                        15.0_dp*rvekt(1)**3*rcube+ &
                                        5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQxx(2, i) = dQxx(2, i)+qh*(-2.0_dp*rvekt(2)*qdenom- &
                                        15.0_dp*rvekt(2)*rvekt(1)**2*rcube+ &
                                        5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQxx(3, i) = dQxx(3, i)+qh*(-2.0_dp*rvekt(3)*qdenom- &
                                        15.0_dp*rvekt(3)*rvekt(1)**2*rcube+ &
                                        5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQyy(1, i) = dQyy(1, i)+qh*(-2.0_dp*rvekt(1)*qdenom- &
                                        15.0_dp*rvekt(1)*rvekt(2)**2*rcube+ &
                                        5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQyy(2, i) = dQyy(2, i)+qh*(4.0_dp*rvekt(2)*qdenom- &
                                        15.0_dp*rvekt(2)**3*rcube+ &
                                        5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQyy(3, i) = dQyy(3, i)+qh*(-2.0_dp*rvekt(3)*qdenom- &
                                        15.0_dp*rvekt(3)*rvekt(2)**2*rcube+ &
                                        5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQzz(1, i) = dQzz(1, i)+qh*(-2.0_dp*rvekt(1)*qdenom- &
                                        15.0_dp*rvekt(1)*rvekt(3)**2*rcube+ &
                                        5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQzz(2, i) = dQzz(2, i)+qh*(-2.0_dp*rvekt(2)*qdenom- &
                                        15.0_dp*rvekt(2)*rvekt(3)**2*rcube+ &
                                        5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQzz(3, i) = dQzz(3, i)+qh*(4.0_dp*rvekt(3)*qdenom- &
                                        15.0_dp*rvekt(3)**3*rcube+ &
                                        5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQxy(1, i) = dQxy(1, i)+qh*(3.0_dp*rvekt(2)*qdenom- &
                                        15.0_dp*rvekt(2)*rvekt(1)**2*rcube)/qdenomsq
            dQxy(2, i) = dQxy(2, i)+qh*(3.0_dp*rvekt(1)*qdenom- &
                                        15.0_dp*rvekt(1)*rvekt(2)**2*rcube)/qdenomsq
            dQxy(3, i) = dQxy(3, i)-qh*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq

            dQxz(1, i) = dQxz(1, i)+qh*(3.0_dp*rvekt(3)*qdenom- &
                                        15.0_dp*rvekt(3)*rvekt(1)**2*rcube)/qdenomsq
            dQxz(2, i) = dQxz(2, i)-qh*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq
            dQxz(3, i) = dQxz(3, i)+qh*(3.0_dp*rvekt(1)*qdenom- &
                                        15.0_dp*rvekt(1)*rvekt(3)**2*rcube)/qdenomsq

            dQyz(1, i) = dQyz(1, i)-qh*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq
            dQyz(2, i) = dQyz(2, i)+qh*(3.0_dp*rvekt(3)*qdenom- &
                                        15.0_dp*rvekt(3)*rvekt(2)**2*rcube)/qdenomsq
            dQyz(3, i) = dQyz(3, i)+qh*(3.0_dp*rvekt(2)*qdenom- &
                                        15.0_dp*rvekt(2)*rvekt(3)**2*rcube)/qdenomsq
         END IF
      END DO

      ! EXTRASITES
      DO i = 1, 5
         !VDW
         DO j = 1, 3
            rvekt(j) = vdwsites(j, i)-helium_x(j)
         END DO
         CALL helium_pbc(helium, rvekt)
         rsqr = rvekt(1)*rvekt(1)+rvekt(2)*rvekt(2)+rvekt(3)*rvekt(3)
         r = SQRT(rsqr)
         rcube = r*rsqr
         rfour = rsqr*rsqr
         rfif = rsqr*rcube
         vdenom = rcube*rcube+deltax**6
         !vdw energy
         Evdw = Evdw+Csixx/vdenom
         IF (PRESENT(force)) THEN
            vdenomsq = vdenom*vdenom
            DO k = 1, 3
               ! Carbon gradient via chain rule
               gradient(k, 1) = gradient(k, 1)-(1.0_dp-fvdw)*6.0_dp*Csixx*rfour*rvekt(k)/vdenomsq
               ! Hydrogen gradient via chain rule
               gradient(k, i+1) = gradient(k, i+1)-fvdw*6.0_dp*Csixx*rfour*rvekt(k)/vdenomsq
            END DO
         END IF
         !INDUCTION
         DO j = 1, 3
            rvekt(j) = indusites(j, i)-helium_x(j)
         END DO
         CALL helium_pbc(helium, rvekt)
         rsqr = rvekt(1)*rvekt(1)+rvekt(2)*rvekt(2)+rvekt(3)*rvekt(3)
         r = SQRT(rsqr)
         rcube = r*rsqr
         rfif = rsqr*rcube
         ddenom = rcube+deltadx**3
         qdenom = rfif+deltaqx**5
         ! dipole induction terms
         DO j = 1, 3
            F(j) = F(j)+qx*rvekt(j)/ddenom
         END DO
         !quadrupole induction terms
         Qxx = Qxx+qx*(3.0_dp*rvekt(1)*rvekt(1)-rsqr)/qdenom
         Qyy = Qyy+qx*(3.0_dp*rvekt(2)*rvekt(2)-rsqr)/qdenom
         Qzz = Qzz+qx*(3.0_dp*rvekt(3)*rvekt(3)-rsqr)/qdenom
         Qxy = Qxy+qx*(3.0_dp*rvekt(1)*rvekt(2))/qdenom
         Qxz = Qxz+qx*(3.0_dp*rvekt(1)*rvekt(3))/qdenom
         Qyz = Qyz+qx*(3.0_dp*rvekt(2)*rvekt(3))/qdenom
         IF (PRESENT(force)) THEN
            ddenomsq = ddenom*ddenom
            !dipole gradient with respect to carbon coordinate
            DO j = 1, 3
               Fprime(j, j, 1) = Fprime(j, j, 1)+qx*(1.0_dp-findu)*(ddenom-3.0_dp*rvekt(j)**2*r)/ddenomsq
            END DO
            !off diagnoal elements
            DO j = 1, 2
               DO k = j+1, 3
                  Fprime(k, j, 1) = Fprime(k, j, 1)+qx*(1.0_dp-findu)*(-3.0_dp*rvekt(j)*rvekt(k)*r)/ddenomsq
               END DO
            END DO
            ! dipole gradient with respect to ith hydrogen coordinate
            DO j = 1, 3
               Fprime(j, j, i+1) = Fprime(j, j, i+1)+qx*findu*(ddenom-3.0_dp*rvekt(j)**2*r)/ddenomsq
            END DO
            !off diagnoal elements
            DO j = 1, 2
               DO k = j+1, 3
                  Fprime(k, j, i+1) = Fprime(k, j, i+1)+qx*findu*(-3.0_dp*rvekt(j)*rvekt(k)*r)/ddenomsq
               END DO
            END DO
            ! Quadrupole field gradient
            qdenomsq = qdenom*qdenom
            ! gradient with respect to the carbon atom
            dQxx(1, 1) = dQxx(1, 1)+qx*(1.0_dp-findu)*(4.0_dp*rvekt(1)*qdenom- &
                                                       15.0_dp*rvekt(1)**3*rcube+ &
                                                       5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQxx(2, 1) = dQxx(2, 1)+qx*(1.0_dp-findu)*(-2.0_dp*rvekt(2)*qdenom- &
                                                       15.0_dp*rvekt(2)*rvekt(1)**2*rcube+ &
                                                       5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQxx(3, 1) = dQxx(3, 1)+qx*(1.0_dp-findu)*(-2.0_dp*rvekt(3)*qdenom- &
                                                       15.0_dp*rvekt(3)*rvekt(1)**2*rcube+ &
                                                       5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQyy(1, 1) = dQyy(1, 1)+qx*(1.0_dp-findu)*(-2.0_dp*rvekt(1)*qdenom- &
                                                       15.0_dp*rvekt(1)*rvekt(2)**2*rcube+ &
                                                       5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQyy(2, 1) = dQyy(2, 1)+qx*(1.0_dp-findu)*(4.0_dp*rvekt(2)*qdenom- &
                                                       15.0_dp*rvekt(2)**3*rcube+ &
                                                       5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQyy(3, 1) = dQyy(3, 1)+qx*(1.0_dp-findu)*(-2.0_dp*rvekt(3)*qdenom- &
                                                       15.0_dp*rvekt(3)*rvekt(2)**2*rcube+ &
                                                       5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQzz(1, 1) = dQzz(1, 1)+qx*(1.0_dp-findu)*(-2.0_dp*rvekt(1)*qdenom- &
                                                       15.0_dp*rvekt(1)*rvekt(3)**2*rcube+ &
                                                       5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQzz(2, 1) = dQzz(2, 1)+qx*(1.0_dp-findu)*(-2.0_dp*rvekt(2)*qdenom- &
                                                       15.0_dp*rvekt(2)*rvekt(3)**2*rcube+ &
                                                       5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQzz(3, 1) = dQzz(3, 1)+qx*(1.0_dp-findu)*(4.0_dp*rvekt(3)*qdenom- &
                                                       15.0_dp*rvekt(3)**3*rcube+ &
                                                       5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQxy(1, 1) = dQxy(1, 1)+qx*(1.0_dp-findu)*(3.0_dp*rvekt(2)*qdenom- &
                                                       15.0_dp*rvekt(2)*rvekt(1)**2*rcube)/qdenomsq
            dQxy(2, 1) = dQxy(2, 1)+qx*(1.0_dp-findu)*(3.0_dp*rvekt(1)*qdenom- &
                                                       15.0_dp*rvekt(1)*rvekt(2)**2*rcube)/qdenomsq
            dQxy(3, 1) = dQxy(3, 1)-qx*(1.0_dp-findu)*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq

            dQxz(1, 1) = dQxz(1, 1)+qx*(1.0_dp-findu)*(3.0_dp*rvekt(3)*qdenom- &
                                                       15.0_dp*rvekt(3)*rvekt(1)**2*rcube)/qdenomsq
            dQxz(2, 1) = dQxz(2, 1)-qx*(1.0_dp-findu)*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq
            dQxz(3, 1) = dQxz(3, 1)+qx*(1.0_dp-findu)*(3.0_dp*rvekt(1)*qdenom- &
                                                       15.0_dp*rvekt(1)*rvekt(3)**2*rcube)/qdenomsq

            dQyz(1, 1) = dQyz(1, 1)-qx*(1.0_dp-findu)*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq
            dQyz(2, 1) = dQyz(2, 1)+qx*(1.0_dp-findu)*(3.0_dp*rvekt(3)*qdenom- &
                                                       15.0_dp*rvekt(3)*rvekt(2)**2*rcube)/qdenomsq
            dQyz(3, 1) = dQyz(3, 1)+qx*(1.0_dp-findu)*(3.0_dp*rvekt(2)*qdenom- &
                                                       15.0_dp*rvekt(2)*rvekt(3)**2*rcube)/qdenomsq

            ! gradient with respect to the hydrogen atom
            dQxx(1, i+1) = dQxx(1, i+1)+qx*findu*(4.0_dp*rvekt(1)*qdenom- &
                                                  15.0_dp*rvekt(1)**3*rcube+ &
                                                  5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQxx(2, i+1) = dQxx(2, i+1)+qx*findu*(-2.0_dp*rvekt(2)*qdenom- &
                                                  15.0_dp*rvekt(2)*rvekt(1)**2*rcube+ &
                                                  5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQxx(3, i+1) = dQxx(3, i+1)+qx*findu*(-2.0_dp*rvekt(3)*qdenom- &
                                                  15.0_dp*rvekt(3)*rvekt(1)**2*rcube+ &
                                                  5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQyy(1, i+1) = dQyy(1, i+1)+qx*findu*(-2.0_dp*rvekt(1)*qdenom- &
                                                  15.0_dp*rvekt(1)*rvekt(2)**2*rcube+ &
                                                  5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQyy(2, i+1) = dQyy(2, i+1)+qx*findu*(4.0_dp*rvekt(2)*qdenom- &
                                                  15.0_dp*rvekt(2)**3*rcube+ &
                                                  5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQyy(3, i+1) = dQyy(3, i+1)+qx*findu*(-2.0_dp*rvekt(3)*qdenom- &
                                                  15.0_dp*rvekt(3)*rvekt(2)**2*rcube+ &
                                                  5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQzz(1, i+1) = dQzz(1, i+1)+qx*findu*(-2.0_dp*rvekt(1)*qdenom- &
                                                  15.0_dp*rvekt(1)*rvekt(3)**2*rcube+ &
                                                  5.0_dp*rvekt(1)*rfif)/qdenomsq
            dQzz(2, i+1) = dQzz(2, i+1)+qx*findu*(-2.0_dp*rvekt(2)*qdenom- &
                                                  15.0_dp*rvekt(2)*rvekt(3)**2*rcube+ &
                                                  5.0_dp*rvekt(2)*rfif)/qdenomsq
            dQzz(3, i+1) = dQzz(3, i+1)+qx*findu*(4.0_dp*rvekt(3)*qdenom- &
                                                  15.0_dp*rvekt(3)**3*rcube+ &
                                                  5.0_dp*rvekt(3)*rfif)/qdenomsq

            dQxy(1, i+1) = dQxy(1, i+1)+qx*findu*(3.0_dp*rvekt(2)*qdenom- &
                                                  15.0_dp*rvekt(2)*rvekt(1)**2*rcube)/qdenomsq
            dQxy(2, i+1) = dQxy(2, i+1)+qx*findu*(3.0_dp*rvekt(1)*qdenom- &
                                                  15.0_dp*rvekt(1)*rvekt(2)**2*rcube)/qdenomsq
            dQxy(3, i+1) = dQxy(3, i+1)-qx*findu*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq

            dQxz(1, i+1) = dQxz(1, i+1)+qx*findu*(3.0_dp*rvekt(3)*qdenom- &
                                                  15.0_dp*rvekt(3)*rvekt(1)**2*rcube)/qdenomsq
            dQxz(2, i+1) = dQxz(2, i+1)-qx*findu*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq
            dQxz(3, i+1) = dQxz(3, i+1)+qx*findu*(3.0_dp*rvekt(1)*qdenom- &
                                                  15.0_dp*rvekt(1)*rvekt(3)**2*rcube)/qdenomsq

            dQyz(1, i+1) = dQyz(1, i+1)-qx*findu*(15.0_dp*rvekt(1)*rvekt(2)*rvekt(3)*rcube)/qdenomsq
            dQyz(2, i+1) = dQyz(2, i+1)+qx*findu*(3.0_dp*rvekt(3)*qdenom- &
                                                  15.0_dp*rvekt(3)*rvekt(2)**2*rcube)/qdenomsq
            dQyz(3, i+1) = dQyz(3, i+1)+qx*findu*(3.0_dp*rvekt(2)*qdenom- &
                                                  15.0_dp*rvekt(2)*rvekt(3)**2*rcube)/qdenomsq
         END IF
      END DO

      ! ANGULAR EXTRASITES
      i = 0
      DO n = 2, 5
         DO m = n+1, 6
            i = i+1
            !VDW
            DO j = 1, 3
               rvekt(j) = vdwwsites(j, i)-helium_x(j)
            END DO
            CALL helium_pbc(helium, rvekt)
            rsqr = rvekt(1)*rvekt(1)+rvekt(2)*rvekt(2)+rvekt(3)*rvekt(3)
            r = SQRT(rsqr)
            rcube = r*rsqr
            rfour = rsqr*rsqr
            rfif = rsqr*rcube
            vdenom = rcube*rcube+deltawx**6
            expbetar = EXP(-betawx*r)
            !vdw energy
            Evdw = Evdw+alphawx*expbetar
            Evdw = Evdw+Csixwx/vdenom
            IF (PRESENT(force)) THEN
               vdenomsq = vdenom*vdenom
               DO k = 1, 3
                  ! Carbon gradient via chain rule
                  gradient(k, 1) = gradient(k, 1)-(1.0_dp-fwvdw)*alphawx*betawx*expbetar*rvekt(k)/r
                  gradient(k, 1) = gradient(k, 1)-(1.0_dp-fwvdw)*6.0_dp*Csixwx*rfour*rvekt(k)/vdenomsq
                  ! nth Hydrogen gradient via chain rule
                  gradient(k, n) = gradient(k, n)-0.5_dp*fwvdw*alphawx*betawx*expbetar*rvekt(k)/r
                  gradient(k, n) = gradient(k, n)-0.5_dp*fwvdw*6.0_dp*Csixwx*rfour*rvekt(k)/vdenomsq
                  ! mth Hydrogen gradient via chain rule
                  gradient(k, m) = gradient(k, m)-0.5_dp*fwvdw*alphawx*betawx*expbetar*rvekt(k)/r
                  gradient(k, m) = gradient(k, m)-0.5_dp*fwvdw*6.0_dp*Csixwx*rfour*rvekt(k)/vdenomsq
               END DO
            END IF
         END DO
      END DO

      ! Dipole energy induced by field
      Eindd = Eindd-0.5_dp*dp_pol_He*(F(1)**2+F(2)**2+F(3)**2)
      !quadrupole induced by field
      ! Grad(alpha,aplha)*Grad(alpha,alpha)
      Eindq = 2.0_dp/3.0_dp*(Qxx*Qxx+Qyy*Qyy+Qzz*Qzz)
      ! Grad(alpha,alpha)*Grad(beta,beta)
      Eindq = Eindq-2.0_dp/3.0_dp*(Qxx*Qyy+Qxx*Qzz+Qyy*Qzz)
      ! Grad(alpha,beta)*Grad(beta,alpha) + Grad(alpha,beta)*Grad(alpha,beta)
      Eindq = Eindq+2.0_dp/1.0_dp*(Qxy*Qxy+Qxz*Qxz+Qyz*Qyz)
      Eindq = -1.0_dp/6.0_dp*qp_pol_He*Eindq
      energy = Evdw+Eindd+Eindq

      IF (PRESENT(force)) THEN
         ! symetrerize dipole gradient tensor
         DO i = 1, 6
            DO j = 1, 2
               DO k = j+1, 3
                  Fprime(j, k, i) = Fprime(k, j, i)
               END DO
            END DO
         END DO
         !gradient of dipole field
         DO j = 1, 6
            DO i = 1, 3
               gradient(i, j) = gradient(i, j)-dp_pol_He*(F(1)*Fprime(1, i, j)+ &
                                                          F(2)*Fprime(2, i, j)+ &
                                                          F(3)*Fprime(3, i, j))
            END DO
         END DO
         !gradient of quadrupole field
         DO j = 1, 6
            DO i = 1, 3
               gradient(i, j) = gradient(i, j)-1.0_dp/6.0_dp*qp_pol_He* &
                                (4.0_dp/3.0_dp*(Qxx*dQxx(i, j)+Qyy*dQyy(i, j)+Qzz*dQzz(i, j)) &
                                 -2.0_dp/3.0_dp*(Qyy*dQxx(i, j)+Qxx*dQyy(i, j)+ &
                                                 Qzz*dQxx(i, j)+Qxx*dQzz(i, j)+ &
                                                 Qyy*dQzz(i, j)+Qzz*dQyy(i, j))+ &
                                 4.0_dp/1.0_dp*(Qxy*dQxy(i, j)+Qxz*dQxz(i, j)+Qyz*dQyz(i, j)))
            END DO
         END DO

         ig = 1
         DO j = 1, helium%solute_atoms
            IF ("C " == helium%solute_element(j)) THEN
               DO i = 1, 3
                  force(3*(j-1)+i) = -gradient(i, j)
               END DO
               ig = ig+1
            END IF
         END DO
         DO j = 1, helium%solute_atoms
            IF ("H " == helium%solute_element(j)) THEN
               DO i = 1, 3
                  force(3*(j-1)+i) = -gradient(i, j)
               END DO
               ig = ig+1
            END IF
         END DO

         DEALLOCATE (Fprime)
         DEALLOCATE (dQxx)
         DEALLOCATE (dQyy)
         DEALLOCATE (dQzz)
         DEALLOCATE (dQxy)
         DEALLOCATE (dQxz)
         DEALLOCATE (dQyz)
         DEALLOCATE (gradient)
      END IF

      DEALLOCATE (F)
      DEALLOCATE (indusites)
      DEALLOCATE (vdwsites)
      DEALLOCATE (vdwwsites)

!IF (PRESENT(force)) THEN
!WRITE(*,*) "@FU: C ", (solute(i,1),i=1,3)
!DO j = 2,6
!   WRITE(*,*) "@FU: H ", (solute(i,j),i=1,3)
!END DO
!DO i = 1, 5
!      WRITE(*,*) "I ", indusites(1,i), indusites(2,i), indusites(3,i)
!   END DO
!   DO i = 1, 5
!      WRITE(*,*) "V ", vdwsites(1,i), vdwsites(2,i), vdwsites(3,i)
!   END DO
!   DO i = 1, 10
!      WRITE(*,*) "W ", vdwwsites(1,i), vdwwsites(2,i), vdwwsites(3,i)
!   END DO
!WRITE(*,*) "@FU: He Pos xyz: ", (helium_x(i),i=1,3)
!WRITE(*,*) "@FU: Eindd: ", Eindd
!WRITE(*,*) "@FU: Eindq: ", Eindq
!WRITE(*,*) "@FU: Eind:  ", Eindd + Eindq
!WRITE(*,*) "@FU: Evdw:  ", Evdw
!WRITE(*,*) "@FU: Eges:  ", Eindd + Eindq + Evdw
!DO i = 1, 6
!      WRITE(*,'(4X,3(1X,G16.9))')   (gradient(j,i), j = 1,3)
!END DO
!
!STOP
!END IF
      RETURN

   END SUBROUTINE helium_intpot_ch5_kuchenbaecker_e_f

! **************************************************************************************************
!> \brief ...
!> \param solute_x ...
!> \param helium ...
!> \param helium_x ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE helium_intpot_seltable_e(solute_x, helium, helium_x, energy)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: solute_x
      TYPE(helium_solvent_type), POINTER                 :: helium
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: helium_x
      REAL(KIND=dp), INTENT(OUT)                         :: energy

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_intpot_seltable_e', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: nr = 41
      REAL(KIND=dp), PARAMETER :: phimax = 360.0_dp, phimin = 0.0_dp, phistep = 7.5_dp, &
         pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164_dp, &
         thetamax = 180.0_dp, thetamin = 0.0_dp, thetastep = 7.5_dp
      INTEGER, PARAMETER :: ntheta = INT((thetamax-thetamin)/thetastep+0.5_dp), &
         nphi = INT((phimax-phimin)/phistep+0.5_dp)

      INTEGER                                            :: i, j, phiindx, readerr, rindx, selindx, &
                                                            thetaindx
      INTEGER, SAVE                                      :: nval = 0
      LOGICAL, SAVE                                      :: readpotential = .FALSE.
      REAL(KIND=dp)                                      :: phi, r, theta
      REAL(KIND=dp), DIMENSION(3)                        :: dvek
      REAL(KIND=dp), POINTER, SAVE                       :: tabpot(:, :)

      !!INCLUDE './heliumpottable.dat'
      IF (.NOT. readpotential) THEN
         readpotential = .TRUE.
         readerr = 0
         OPEN (137, FILE="./heliumpottable.dat", FORM="formatted", STATUS="old")
         DO WHILE (readerr == 0)
            READ (137, *, iostat=readerr) energy
            IF (readerr == 0) nval = nval+1
         END DO
         CLOSE (137)
         WRITE (*, *) "Allocating for ", nval, "points"
         ALLOCATE (tabpot(4, nval))
         WRITE (*, *) "reading potential"
         OPEN (137, FILE="./heliumpottable.dat", FORM="formatted", STATUS="old")
         DO j = 1, nval
            READ (137, *, iostat=readerr) (tabpot(i, j), i=1, 4)
            IF (readerr /= 0) WRITE (*, *) "FUHL READ ERROR!!!"
            IF (readerr /= 0) STOP
         END DO
         CLOSE (137)
         WRITE (*, *) "read potential"
      END IF

      DO j = 1, 3
         dvek(j) = helium_x(j)
      END DO
      CALL helium_pbc(helium, dvek)
      energy = 0.0_dp
      r = SQRT(dvek(1)**2+dvek(2)**2+dvek(3)**2)
      IF (r < 1.0_dp) THEN
         energy = 10.0_dp**6
         RETURN
      END IF
      theta = dvek(3)/r
      IF (ABS(theta) > REAL(1.0-EPSILON(1.0), dp)) THEN
         phi = 0.0_dp
      ELSE
         phi = ATAN2(dvek(2), dvek(1))*180.0_dp/pi
         IF (phi < phimin) phi = phimax+phi
      END IF
      theta = ACOS(theta)*180.0_dp/pi

      ! get theta index
      IF (theta < thetamin+thetastep/2.0_dp) THEN
         thetaindx = 0
      ELSE IF (theta > thetamax-thetastep/2.0_dp) THEN
         thetaindx = ntheta
      ELSE
         thetaindx = INT((theta-thetamin)/thetastep+0.5_dp)
      END IF

      ! get phi index
      IF (phi > phimax-phistep/2.0_dp) THEN
         phiindx = 1
      ELSE IF (thetaindx == ntheta) THEN
         phiindx = 0
      ELSE
         phiindx = INT((phi-phimin)/phistep+0.5_dp)
      END IF

      ! get r index
      !rindx = MIN(MAX(INT((r - 2.5)/0.25 + 0.5_dp) + 1 , 1),nr)
      r = r*0.529177249_dp
      IF (r < 2.5_dp) THEN
         rindx = 1
      ELSE IF ((r >= 2.5_dp) .AND. (r < 4.9_dp)) THEN
         rindx = MIN(MAX(INT((r-2.5)/0.1_dp+0.5_dp)+1, 2), 24)
      ELSE IF ((r >= 4.9_dp) .AND. (r < 5.25_dp)) THEN
         rindx = 25
      ELSE IF ((r >= 5.25_dp) .AND. (r < 13.0_dp)) THEN
         rindx = INT((r-5_dp)/0.5_dp+0.5_dp)+25
      ELSE
         rindx = 41
      END IF

      ! get final index
      selindx = 0
      IF (thetaindx == 0) THEN
         selindx = rindx
      ELSE IF (thetaindx == INT(thetamax/thetastep+0.5_dp)) THEN
         selindx = rindx+nr+(ntheta-1)*nphi*nr
      ELSE
         selindx = nr+(thetaindx-1)*nphi*nr+phiindx*nr+rindx
      END IF

      energy = tabpot(4, selindx)

      RETURN

   END SUBROUTINE helium_intpot_seltable_e
! ***************************************************************************
!> \brief Calculate l-th Legendre polynomial P_{l}(x) at point x
!> \param x ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   FUNCTION Pl(x, n)
      REAL(KIND=dp), INTENT(IN)                          :: x
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp)                                      :: Pl

      INTEGER                                            :: k
      REAL(KIND=dp)                                      :: pln(0:n)

      pln(0) = 1.0_dp
      pln(1) = x

      IF (n <= 1) THEN
         Pl = pln(n)
      ELSE
         DO k = 1, n-1
            pln(k+1) = ((2.0*k+1.0)*x*pln(k)-REAL(k, dp)*pln(k-1))/(REAL(k+1, dp))
         END DO
         Pl = pln(n)
      END IF
      RETURN
   END FUNCTION Pl

! ***************************************************************************
!> \brief  Scan the helium-solute interaction energy within the periodic cell
!> \param pint_env ...
!> \param helium_env ...
!> \date   2014-01-22
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_intpot_scan(pint_env, helium_env)

      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_intpot_scan', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ic, ix, iy, iz, k, nbin
      LOGICAL                                            :: wrapped
      REAL(KIND=dp)                                      :: delr, my_en, ox, oy, oz
      REAL(kind=dp), DIMENSION(3)                        :: pbc1, pbc2, pos

      CALL timeset(routineN, handle)

      ! Perform scan only on ionode, since this is only used to output the intpot
      IF (pint_env%logger%para_env%ionode) THEN
         ! Assume ionode always to have at least one helium_env
         k = 1
         helium_env(k)%helium%rho_inst(1, :, :, :) = 0.0_dp
         nbin = helium_env(k)%helium%rho_nbin
         delr = helium_env(k)%helium%rho_delr
         ox = helium_env(k)%helium%center(1)-helium_env(k)%helium%rho_maxr/2.0_dp
         oy = helium_env(k)%helium%center(2)-helium_env(k)%helium%rho_maxr/2.0_dp
         oz = helium_env(k)%helium%center(3)-helium_env(k)%helium%rho_maxr/2.0_dp

         DO ix = 1, nbin
            DO iy = 1, nbin
               DO iz = 1, nbin

                  ! put the probe in the center of the current voxel
                  pos(:) = (/ox+(ix-0.5_dp)*delr, oy+(iy-0.5_dp)*delr, oz+(iz-0.5_dp)*delr/)

                  ! calc interaction energy for the current probe position
                  helium_env(k)%helium%pos(:, 1, 1) = pos(:)
                  CALL helium_bead_solute_e_f(pint_env, helium_env(k)%helium, 1, 1, energy=my_en)

                  ! check if the probe fits within the unit cell
                  pbc1(:) = pos(:)-helium_env(k)%helium%center
                  pbc2(:) = pbc1(:)
                  CALL helium_pbc(helium_env(k)%helium, pbc2)
                  wrapped = .FALSE.
                  DO ic = 1, 3
                     IF (ABS(pbc1(ic)-pbc2(ic)) .GT. 10.0_dp*EPSILON(0.0_dp)) THEN
                        wrapped = .TRUE.
                     END IF
                  END DO

                  ! set the interaction energy value
                  IF (wrapped) THEN
                     helium_env(k)%helium%rho_inst(1, ix, iy, iz) = 0.0_dp
                  ELSE
                     helium_env(k)%helium%rho_inst(1, ix, iy, iz) = my_en
                  END IF

               END DO
            END DO
         END DO
      END IF

      CALL timestop(handle)
      RETURN
   END SUBROUTINE helium_intpot_scan

! ***************************************************************************
!> \brief Calculate model helium-solute interaction energy and forces
!>        between one helium bead and the corresponding solute time
!>        slice asuming water solute.
!> \param solute_x  solute positions ARR(3*NATOMS)
!>        to global atom indices
!> \param helium    only needed for helium_pbc call at the moment
!> \param helium_x  helium bead position ARR(3)
!> \param energy    calculated interaction energy
!> \param force ...
!> \author Felix Uhl
! **************************************************************************************************
   SUBROUTINE helium_intpot_model_water(solute_x, helium, helium_x, energy, force)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: solute_x
      TYPE(helium_solvent_type), POINTER                 :: helium
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: helium_x
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL, POINTER                               :: force

      CHARACTER(LEN=*), PARAMETER :: routineN = 'helium_intpot_model_water', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ig
      REAL(KIND=dp)                                      :: d, d2, dd, ep, eps, s1, s2, sig
      REAL(KIND=dp), DIMENSION(3)                        :: dr, solute_r

      energy = 0.0_dp
      IF (PRESENT(force)) THEN
         force(:) = 0.0_dp
      END IF

      sig = 2.69_dp ! 1.4 Angstrom
      eps = 60.61e-6_dp ! 19 K
      s1 = 0.0_dp
      DO i = 1, SIZE(helium%solute_element)
         IF (helium%solute_element(i) == "H ") THEN
            ig = i-1
            solute_r(1) = solute_x(3*ig+1)
            solute_r(2) = solute_x(3*ig+2)
            solute_r(3) = solute_x(3*ig+3)
            dr(:) = solute_r(:)-helium_x(:)
            CALL helium_pbc(helium, dr)
            d2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            d = SQRT(d2)
            dd = (sig/d)**6
            ep = 4.0_dp*eps*dd*(dd-1.0_dp)
            s1 = s1+ep
            s2 = 24.0_dp*eps*dd*(2.0_dp*dd-1.0_dp)/d2
            IF (PRESENT(force)) THEN
               force(3*ig+1) = force(3*ig+1)+s2*dr(1)
               force(3*ig+2) = force(3*ig+2)+s2*dr(2)
               force(3*ig+3) = force(3*ig+3)+s2*dr(3)
            END IF
         END IF
      END DO ! i = 1, num_hydrogen
      energy = energy+s1

      sig = 5.01_dp ! 2.6 Angstrom
      eps = 104.5e-6_dp ! 33 K
      s1 = 0.0_dp
      DO i = 1, SIZE(helium%solute_element)
         IF (helium%solute_element(i) == "O ") THEN
            ig = i-1
            solute_r(1) = solute_x(3*ig+1)
            solute_r(2) = solute_x(3*ig+2)
            solute_r(3) = solute_x(3*ig+3)
            dr(:) = solute_r(:)-helium_x(:)
            CALL helium_pbc(helium, dr)
            d2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
            d = SQRT(d2)
            dd = (sig/d)**6
            ep = 4.0_dp*eps*dd*(dd-1.0_dp)
            s1 = s1+ep
            s2 = 24.0_dp*eps*dd*(2.0_dp*dd-1.0_dp)/d2
            IF (PRESENT(force)) THEN
               force(3*ig+1) = force(3*ig+1)+s2*dr(1)
               force(3*ig+2) = force(3*ig+2)+s2*dr(2)
               force(3*ig+3) = force(3*ig+3)+s2*dr(3)
            END IF
         END IF
      END DO ! i = 1, num_chlorine
      energy = energy+s1

      RETURN

   END SUBROUTINE helium_intpot_model_water

! **************************************************************************************************
!> \brief ...
!> \param solute_x ...
!> \param helium ...
!> \param helium_x ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE helium_intpot_cubetable(solute_x, helium, helium_x, energy)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: solute_x
      TYPE(helium_solvent_type), POINTER                 :: helium
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: helium_x
      REAL(KIND=dp), INTENT(OUT)                         :: energy

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_intpot_cubetable', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: dummy
      INTEGER                                            :: e_ind, handle, i, ix, iy, iz, readerr
      INTEGER, DIMENSION(3)                              :: indices
      INTEGER, DIMENSION(3), SAVE                        :: cubenvox
      INTEGER, SAVE                                      :: cubenat
      LOGICAL, SAVE                                      :: readpot = .FALSE.
      REAL(KIND=dp)                                      :: thresh
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         SAVE                                            :: pottable
      REAL(KIND=dp), DIMENSION(3)                        :: dvek
      REAL(KIND=dp), DIMENSION(3), SAVE                  :: cubeorig
      REAL(KIND=dp), DIMENSION(3, 3), SAVE               :: cubevec

      CALL timeset(routineN, handle)

      energy = 0.0_dp
      indices(:) = 0

      ! Read cube file:
      IF (.NOT. readpot) THEN
         readpot = .TRUE.
         cubenat = 0
         cubenvox(:) = 0
         cubeorig(:) = 0.0_dp
         cubevec(:, :) = 0.0_dp
         OPEN (137, FILE="./heliumpot.cube", FORM="formatted", STATUS="old")
         ! Read header:
         READ (137, *, iostat=readerr) dummy
         READ (137, *, iostat=readerr) dummy
         READ (137, *, iostat=readerr) cubenat, cubeorig(:)
         DO i = 1, 3
            READ (137, *, iostat=readerr) cubenvox(i), cubevec(:, i)
         END DO
         DO i = 1, cubenat
            READ (137, *, iostat=readerr) dummy
         END DO

         !WRITE(*,*) cubenat, cubeorig
         !WRITE(*,*) cubenvox, cubeorig

         !ndim=cubenvox(1)*cubenvox(2)*cubenvox(3)
         ALLOCATE (pottable(cubenvox(1), cubenvox(2), cubenvox(3)))
         pottable(:, :, :) = 0.0_dp

         i = 1
         ! Read pot:
         DO ix = 1, cubenvox(1)
            DO iy = 1, cubenvox(2)
               DO iz = 1, cubenvox(3)
                  READ (137, *, iostat=readerr) pottable(ix, iy, iz)
                  IF (readerr /= 0) WRITE (*, *) "READ ERROR, PROVIDE CUBE IN RIGHT FORMAT!!!"
                  IF (readerr /= 0) STOP
                  i = i+1
               END DO
            END DO
         END DO
         CLOSE (137)
      END IF

      DO i = 1, 3
         dvek(i) = helium_x(i)
      END DO
      CALL helium_pbc(helium, dvek)

      ! Get indices:
      indices(:) = 0
      !WRITE(*,*) dvek
      thresh = 0.00001_dp
      DO i = 1, 3
         IF (dvek(i)-cubeorig(i) <= thresh) THEN
            indices(i) = 1
            ! Assume cube cubes only!!
         ELSE IF (dvek(i)-(cubeorig(i)+(cubenvox(i)-1)*cubevec(i, i)) >= thresh) THEN
            indices(i) = cubenvox(i)
         ELSE
            indices(i) = CEILING((dvek(i)-cubeorig(i))/cubevec(i, i)+0.5_dp)
         END IF
      END DO
      !WRITE(*,*) indices

      ! calculate energy:
      e_ind = indices(3)+(indices(2)-1)*cubenvox(3)+(indices(1)-1)*cubenvox(3)*cubenvox(2)
      !WRITE(*,*) pottable(73,18,63)
      energy = pottable(indices(1), indices(2), indices(3))

      ! Deallocate:
      !DEALLOCATE (pottable)

      CALL timestop(handle)

      RETURN

   END SUBROUTINE helium_intpot_cubetable
! ***************************************************************************
!> \brief Helium-helium pair interaction potential.
!> \param r ...
!> \return ...
! **************************************************************************************************
   FUNCTION helium_vij(r) RESULT(vij)

      REAL(kind=dp), INTENT(IN)                          :: r
      REAL(kind=dp)                                      :: vij

      REAL(kind=dp)                                      :: f, x, x2

      x = angstrom*r/2.9673_dp
      IF (x < 1.241314_dp) THEN
         x2 = 1.241314_dp/x-1.0_dp
         f = EXP(-x2*x2)
      ELSE
         f = 1.0_dp
      END IF
      x2 = 1.0_dp/(x*x)
      vij = 10.8_dp/kelvin*(544850.4_dp*EXP(-13.353384_dp*x)-f* &
                            ((0.1781_dp*x2+0.4253785_dp)*x2+1.3732412_dp)*x2*x2*x2)
      RETURN
   END FUNCTION helium_vij

#if 0

   ! this block is currently turned off

! ***************************************************************************
!> \brief Helium-helium pair interaction potential's derivative.
!> \param r ...
!> \return ...
! **************************************************************************************************
   FUNCTION helium_d_vij(r) RESULT(dvij)

      REAL(kind=dp), INTENT(IN)                          :: r
      REAL(kind=dp)                                      :: dvij

      REAL(kind=dp)                                      :: f, fp, x, x2, y

      x = angstrom*r/2.9673_dp
      x = r/2.9673_dp
      x2 = 1.0_dp/(x*x)
      IF (x < 1.241314_dp) THEN
         y = 1.241314_dp/x-1.0_dp
         f = EXP(-y*y)
         fp = 2.0_dp*1.241314_dp*f*y* &
              ((0.1781_dp*x2+0.4253785_dp)*x2+1.3732412_dp)*x2*x2*x2*x2
      ELSE
         f = 1.0_dp
         fp = 0.0_dp
      END IF

      dvij = angstrom*(10.8_dp/2.9673_dp)*( &
             (-13.353384_dp*544850.4_dp)*EXP(-13.353384_dp*x)-fp+ &
             f*(((10.0_dp*0.1781_dp)*x2+(8.0_dp*0.4253785_dp))*x2+(6.0_dp*1.3732412_dp))* &
             x2*x2*x2/x)/(r*kelvin)
      RETURN
   END FUNCTION helium_d_vij

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param n ...
!> \param i ...
!> \return ...
! **************************************************************************************************
   FUNCTION helium_atom_action(helium, n, i) RESULT(res)

      TYPE(helium_solvent_type), POINTER                 :: helium
      INTEGER, INTENT(IN)                                :: n, i
      REAL(KIND=dp)                                      :: res

      INTEGER                                            :: c, j
      REAL(KIND=dp)                                      :: r(3), rp(3), s, t

      s = 0.0_dp
      t = 0.0_dp
      IF (n < helium%beads) THEN
         DO c = 1, 3
            r(c) = helium%pos(c, i, n)-helium%pos(c, i, n+1)
         END DO
         CALL helium_pbc(helium, r)
         t = r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
         DO j = 1, i-1
            DO c = 1, 3
               r(c) = helium%pos(c, i, n)-helium%pos(c, j, n)
               rp(c) = helium%pos(c, i, n+1)-helium%pos(c, j, n+1)
            END DO
            CALL helium_pbc(helium, r)
            CALL helium_pbc(helium, rp)
            s = s+helium_eval_expansion(helium, r, rp, helium%uij, 1)
         END DO
         DO j = i+1, helium%atoms
            DO c = 1, 3
               r(c) = helium%pos(c, i, n)-helium%pos(c, j, n)
               rp(c) = helium%pos(c, i, n+1)-helium%pos(c, j, n+1)
            END DO
            CALL helium_pbc(helium, r)
            CALL helium_pbc(helium, rp)
            s = s+helium_eval_expansion(helium, r, rp, helium%uij, 1)
         END DO
      ELSE
         DO c = 1, 3
            r(c) = helium%pos(c, i, n)-helium%pos(c, helium%permutation(i), 1)
         END DO
         CALL helium_pbc(helium, r)
         t = r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
         DO j = 1, i-1
            DO c = 1, 3
               r(c) = helium%pos(c, i, n)-helium%pos(c, j, n)
               rp(c) = helium%pos(c, helium%permutation(i), 1)-helium%pos(c, helium%permutation(j), 1)
            END DO
            CALL helium_pbc(helium, r)
            CALL helium_pbc(helium, rp)
            s = s+helium_eval_expansion(helium, r, rp, helium%uij, 1)
         END DO
         DO j = i+1, helium%atoms
            DO c = 1, 3
               r(c) = helium%pos(c, i, n)-helium%pos(c, j, n)
               rp(c) = helium%pos(c, helium%permutation(i), 1)-helium%pos(c, helium%permutation(j), 1)
            END DO
            CALL helium_pbc(helium, r)
            CALL helium_pbc(helium, rp)
            s = s+helium_eval_expansion(helium, r, rp, helium%uij, 1)
         END DO
      END IF
      t = t/(2.0_dp*helium%tau*helium%hb2m)
      s = s*0.5_dp
      res = s+t
      RETURN

   END FUNCTION helium_atom_action

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   FUNCTION helium_link_action(helium, n) RESULT(res)

      TYPE(helium_solvent_type), POINTER                 :: helium
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp)                                      :: res

      INTEGER                                            :: c, i, j
      REAL(KIND=dp)                                      :: r(3), rp(3), s, t

      s = 0.0_dp
      t = 0.0_dp
      IF (n < helium%beads) THEN
         DO i = 1, helium%atoms
            DO c = 1, 3
               r(c) = helium%pos(c, i, n)-helium%pos(c, i, n+1)
            END DO
            CALL helium_pbc(helium, r)
            t = t+r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
            DO j = 1, i-1
               DO c = 1, 3
                  r(c) = helium%pos(c, i, n)-helium%pos(c, j, n)
                  rp(c) = helium%pos(c, i, n+1)-helium%pos(c, j, n+1)
               END DO
               CALL helium_pbc(helium, r)
               CALL helium_pbc(helium, rp)
               s = s+helium_eval_expansion(helium, r, rp, helium%uij, 1)
            END DO
         END DO
      ELSE
         DO i = 1, helium%atoms
            DO c = 1, 3
               r(c) = helium%pos(c, i, n)-helium%pos(c, helium%permutation(i), 1)
            END DO
            CALL helium_pbc(helium, r)
            t = t+r(1)*r(1)+r(2)*r(2)+r(3)*r(3)
            DO j = 1, i-1
               DO c = 1, 3
                  r(c) = helium%pos(c, i, n)-helium%pos(c, j, n)
                  rp(c) = helium%pos(c, helium%permutation(i), 1)-helium%pos(c, helium%permutation(j), 1)
               END DO
               CALL helium_pbc(helium, r)
               CALL helium_pbc(helium, rp)
               s = s+helium_eval_expansion(helium, r, rp, helium%uij, 1)
            END DO
         END DO
      END IF
      t = t/(2.0_dp*helium%tau*helium%hb2m)
      res = s+t
      RETURN

   END FUNCTION helium_link_action

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \return ...
! **************************************************************************************************
   FUNCTION helium_total_action(helium) RESULT(res)

      TYPE(helium_solvent_type), POINTER                 :: helium
      REAL(KIND=dp)                                      :: res

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: s

      s = 0.0_dp
      DO i = 1, helium%beads
         s = s+helium_link_action(helium, i)
      END DO
      res = s
      RETURN

   END FUNCTION helium_total_action

! **************************************************************************************************
!> \brief ...
!> \param helium ...
!> \param part ...
!> \param ref_bead ...
!> \param delta_bead ...
!> \param d ...
! **************************************************************************************************
   SUBROUTINE helium_delta_pos(helium, part, ref_bead, delta_bead, d)

      TYPE(helium_solvent_type), POINTER                 :: helium
      INTEGER, INTENT(IN)                                :: part, ref_bead, delta_bead
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: d

      INTEGER                                            :: b, bead, db, nbead, np, p
      REAL(KIND=dp), DIMENSION(3)                        :: r

      b = helium%beads

      d(:) = 0.0_dp
      IF (delta_bead > 0) THEN
         bead = ref_bead
         p = part
         db = delta_bead
         DO
            IF (db < 1) EXIT
            nbead = bead+1
            np = p
            IF (nbead > b) THEN
               nbead = nbead-b
               np = helium%permutation(np)
            END IF
            r(:) = helium%pos(:, p, bead)-helium%pos(:, np, nbead)
            CALL helium_pbc(helium, r)
            d(:) = d(:)+r(:)
            bead = nbead
            p = np
            db = db-1
         END DO
      ELSEIF (delta_bead < 0) THEN
         bead = ref_bead
         p = part
         db = delta_bead
         DO
            IF (db >= 0) EXIT
            nbead = bead-1
            np = p
            IF (nbead < 1) THEN
               nbead = nbead+b
               np = helium%iperm(np)
            END IF
            r(:) = helium%pos(:, p, bead)-helium%pos(:, np, nbead)
            CALL helium_pbc(helium, r)
            d(:) = d(:)+r(:)
            bead = nbead
            p = np
            db = db+1
         END DO
      END IF
      RETURN
   END SUBROUTINE helium_delta_pos

#endif

END MODULE helium_interactions
