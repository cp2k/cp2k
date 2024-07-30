!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of charge equilibration in xTB
!> \author JGH
! **************************************************************************************************
MODULE xtb_eeq
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE eeq_method,                      ONLY: eeq_solver
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE particle_types,                  ONLY: particle_type
   USE qs_dispersion_pairpot,           ONLY: dcnum_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE spme,                            ONLY: spme_forces
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_eeq'

   PUBLIC :: xtb_eeq_calculation, xtb_eeq_forces

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param cnumbers ...
!> \param eeq_energy ...
!> \param lambda ...
! **************************************************************************************************
   SUBROUTINE xtb_eeq_calculation(qs_env, charges, cnumbers, eeq_energy, lambda)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cnumbers
      REAL(KIND=dp), INTENT(INOUT)                       :: eeq_energy, lambda

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_eeq_calculation'

      INTEGER                                            :: handle, iatom, ikind, iunit, jkind, &
                                                            natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      LOGICAL                                            :: defined, do_ewald, do_sparse
      REAL(KIND=dp)                                      :: ala, alb, gama, kappa, scn, totalcharge, &
                                                            xi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chia, gam
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      iunit = cp_logger_get_default_unit_nr()

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      atprop=atprop, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)

      xtb_control => dft_control%qs_control%xtb_control

      totalcharge = dft_control%charge

      IF (atprop%energy) THEN
         CALL atprop_array_init(atprop%atecoul, natom)
      END IF

      ! gamma[a,b]
      ALLOCATE (gab(nkind, nkind), gam(nkind))
      gab = 0.0_dp
      gam = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
         IF (.NOT. defined) CYCLE
         CALL get_xtb_atom_param(xtb_atom_a, alpg=ala, eta=gama)
         gam(ikind) = gama
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
            IF (.NOT. defined) CYCLE
            CALL get_xtb_atom_param(xtb_atom_b, alpg=alb)
            !
            gab(ikind, jkind) = SQRT(1._dp/(ala*ala + alb*alb))
            !
         END DO
      END DO

      ! Chi[a,a]
      ALLOCATE (chia(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, xi=xi, kappa0=kappa)
         !
         scn = SQRT(cnumbers(iatom)) + 1.0e-14_dp
         chia(iatom) = xi - kappa*scn
         !
      END DO

      do_ewald = xtb_control%do_ewald
      do_sparse = .FALSE.
      IF (natom > 100000 .AND. do_ewald) do_sparse = .TRUE.

      IF (do_ewald) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL eeq_solver(qs_env, charges, lambda, eeq_energy, &
                         particle_set, kind_of, cell, chia, gam, gab, &
                         totalcharge=totalcharge, ewald=do_ewald, sparse=do_sparse, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw, iounit=iunit)
      ELSE
         CALL eeq_solver(qs_env, charges, lambda, eeq_energy, &
                         particle_set, kind_of, cell, chia, gam, gab, &
                         totalcharge=totalcharge, iounit=iunit)
      END IF

      DEALLOCATE (gab, gam, chia)

      CALL timestop(handle)

   END SUBROUTINE xtb_eeq_calculation

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param dcharges ...
!> \param cnumbers ...
!> \param dcnum ...
! **************************************************************************************************
   SUBROUTINE xtb_eeq_forces(qs_env, charges, dcharges, cnumbers, dcnum)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges, dcharges, cnumbers
      TYPE(dcnum_type), DIMENSION(:), INTENT(IN)         :: dcnum

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xtb_eeq_forces'

      INTEGER                                            :: atom_a, atom_b, atom_c, handle, i, ia, &
                                                            iatom, ikind, iunit, jatom, jkind, &
                                                            katom, kkind, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL                                            :: defined, do_ewald, do_sparse, use_virial
      REAL(KIND=dp)                                      :: ala, alb, alpha, ctot, dr, dr2, drk, &
                                                            elag, gam2, gama, grc, kappa, qlam, &
                                                            qq, scn, totalcharge, xi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: gam
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: epforce, gab
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, ri, rij, rik, rj
      REAL(KIND=dp), DIMENSION(:), POINTER               :: chrgx, dchia, qlag
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      iunit = cp_logger_get_default_unit_nr()

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      atprop=atprop, &
                      force=force, &
                      virial=virial, &
                      cell=cell, &
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      xtb_control => dft_control%qs_control%xtb_control

      totalcharge = dft_control%charge

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      ! gamma[a,b]
      ALLOCATE (gab(nkind, nkind), gam(nkind))
      gab = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
         IF (.NOT. defined) CYCLE
         CALL get_xtb_atom_param(xtb_atom_a, alpg=ala, eta=gama)
         gam(ikind) = gama
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
            IF (.NOT. defined) CYCLE
            CALL get_xtb_atom_param(xtb_atom_b, alpg=alb)
            !
            gab(ikind, jkind) = SQRT(1._dp/(ala*ala + alb*alb))
            !
         END DO
      END DO

      ALLOCATE (qlag(natom))

      do_ewald = xtb_control%do_ewald
      do_sparse = .FALSE.
      IF (natom > 100000 .AND. do_ewald) do_sparse = .TRUE.

      IF (do_ewald) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL eeq_solver(qs_env, qlag, qlam, elag, &
                         particle_set, kind_of, cell, dcharges, gam, gab, &
                         ewald=do_ewald, sparse=do_sparse, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw, iounit=iunit)
         qlag = -qlag
         qlam = -qlam
      ELSE
         CALL eeq_solver(qs_env, qlag, qlam, elag, &
                         particle_set, kind_of, cell, -dcharges, gam, gab, iounit=iunit)
      END IF

      ALLOCATE (chrgx(natom), dchia(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, xi=xi, kappa0=kappa, eta=gama)
         !
         scn = SQRT(cnumbers(iatom)) + 1.0e-14_dp
         ctot = 0.5_dp*(charges(iatom) - qlag(iatom))
         dchia(iatom) = -ctot*kappa/scn
      END DO

      ! Forces from q*X
      CALL get_qs_env(qs_env=qs_env, &
                      local_particles=local_particles)
      DO ikind = 1, nkind
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            atom_a = atom_of_kind(iatom)
            DO i = 1, dcnum(iatom)%neighbors
               katom = dcnum(iatom)%nlist(i)
               kkind = kind_of(katom)
               atom_c = atom_of_kind(katom)
               rik = dcnum(iatom)%rik(:, i)
               drk = SQRT(SUM(rik(:)**2))
               IF (drk > 1.e-3_dp) THEN
                  fdik(:) = dchia(iatom)*dcnum(iatom)%dvals(i)*rik(:)/drk
                  force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) - fdik(:)
                  force(kkind)%rho_elec(:, atom_c) = force(kkind)%rho_elec(:, atom_c) + fdik(:)
                  IF (use_virial) THEN
                     CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                     IF (atprop%stress) THEN
                        CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
                        CALL virial_pair_force(atprop%atstress(:, :, katom), -0.5_dp, fdik, rik)
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO
      ! Forces from (0.5*q+l)*dA/dR*q
      IF (do_ewald) THEN
         CALL ewald_env_get(ewald_env, alpha=alpha)
      END IF
      DO ikind = 1, nkind
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            atom_a = atom_of_kind(iatom)
            ri(1:3) = particle_set(iatom)%r(1:3)
            DO jatom = 1, natom
               IF (iatom == jatom) CYCLE
               jkind = kind_of(jatom)
               atom_b = atom_of_kind(jatom)
               qq = (0.5_dp*charges(iatom) - qlag(iatom))*charges(jatom)
               rj(1:3) = particle_set(jatom)%r(1:3)
               rij(1:3) = ri(1:3) - rj(1:3)
               rij = pbc(rij, cell)
               dr2 = SUM(rij**2)
               dr = SQRT(dr2)
               gama = gab(ikind, jkind)
               gam2 = gama*gama
               grc = 2._dp*gama*EXP(-gam2*dr2)*oorootpi/dr - erf(gama*dr)/dr2
               IF (do_ewald) THEN
                  grc = grc - 2._dp*alpha*EXP(-alpha**2*dr2)*oorootpi/dr + erf(alpha*dr)/dr2
               END IF
               fdik(:) = qq*grc*rij(:)/dr
               force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) + fdik(:)
               force(jkind)%rho_elec(:, atom_b) = force(jkind)%rho_elec(:, atom_b) - fdik(:)
               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rij)
                  IF (atprop%stress) THEN
                     CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rij)
                     CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp, fdik, rij)
                  END IF
               END IF
            END DO
         END DO
      END DO
      ! Forces from Ewald potential: (q+l)*A*q
      IF (do_ewald) THEN
         ALLOCATE (epforce(3, natom))
         epforce = 0.0_dp
         dchia = -charges + qlag
         chrgx = charges
         CALL spme_forces(ewald_env, ewald_pw, cell, particle_set, chrgx, &
                          particle_set, dchia, epforce)
         dchia = charges
         chrgx = qlag
         CALL spme_forces(ewald_env, ewald_pw, cell, particle_set, chrgx, &
                          particle_set, dchia, epforce)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            i = atom_of_kind(iatom)
            force(ikind)%rho_elec(:, i) = force(ikind)%rho_elec(:, i) + epforce(:, iatom)
            IF (use_virial) THEN
               fdik(1:3) = epforce(1:3, iatom)
               CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
               END IF
            END IF
         END DO
         DEALLOCATE (epforce)
      END IF

      DEALLOCATE (gab, chrgx, dchia, qlag)

      CALL timestop(handle)

   END SUBROUTINE xtb_eeq_forces

END MODULE xtb_eeq
