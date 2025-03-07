!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_unit_nr
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE eeq_input,                       ONLY: eeq_solver_type
   USE eeq_method,                      ONLY: eeq_efield_energy,&
                                              eeq_efield_force_loc,&
                                              eeq_efield_force_periodic,&
                                              eeq_efield_pot,&
                                              eeq_solver
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_dispersion_cnum,              ONLY: dcnum_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE spme,                            ONLY: spme_forces,&
                                              spme_virial
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_eeq'

   PUBLIC :: xtb_eeq_calculation, xtb_eeq_forces, xtb_eeq_lagrange

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param cnumbers ...
!> \param eeq_sparam ...
!> \param eeq_energy ...
!> \param ef_energy ...
!> \param lambda ...
! **************************************************************************************************
   SUBROUTINE xtb_eeq_calculation(qs_env, charges, cnumbers, &
                                  eeq_sparam, eeq_energy, ef_energy, lambda)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: charges
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cnumbers
      TYPE(eeq_solver_type), INTENT(IN)                  :: eeq_sparam
      REAL(KIND=dp), INTENT(INOUT)                       :: eeq_energy, ef_energy, lambda

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_eeq_calculation'

      INTEGER                                            :: enshift_type, handle, iatom, ikind, &
                                                            iunit, jkind, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      LOGICAL                                            :: defined, do_ewald
      REAL(KIND=dp)                                      :: ala, alb, cn, esg, gama, kappa, scn, &
                                                            sgamma, totalcharge, xi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: chia, efr, gam
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(mp_para_env_type), POINTER                    :: para_env
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
      enshift_type = xtb_control%enshift_type
      IF (enshift_type == 0) THEN
         enshift_type = 2
         IF (ALL(cell%perd == 0)) enshift_type = 1
      END IF
      sgamma = 8.0_dp ! see D4 for periodic systems paper
      esg = 1.0_dp + EXP(sgamma)
      ALLOCATE (chia(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, xi=xi, kappa0=kappa)
         !
         IF (enshift_type == 1) THEN
            scn = SQRT(cnumbers(iatom)) + 1.0e-14_dp
         ELSE IF (enshift_type == 2) THEN
            cn = cnumbers(iatom)/esg
            scn = LOG(esg/(esg - cnumbers(iatom)))
         ELSE
            CPABORT("Unknown enshift_type")
         END IF
         chia(iatom) = xi - kappa*scn
         !
      END DO

      ef_energy = 0.0_dp
      IF (dft_control%apply_period_efield .OR. dft_control%apply_efield .OR. &
          dft_control%apply_efield_field) THEN
         ALLOCATE (efr(natom))
         efr(1:natom) = 0.0_dp
         CALL eeq_efield_pot(qs_env, efr)
         chia(1:natom) = chia(1:natom) + efr(1:natom)
      END IF

      do_ewald = xtb_control%do_ewald

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      IF (do_ewald) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL eeq_solver(charges, lambda, eeq_energy, &
                         particle_set, kind_of, cell, chia, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, &
                         totalcharge=totalcharge, ewald=do_ewald, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw, iounit=iunit)
      ELSE
         CALL eeq_solver(charges, lambda, eeq_energy, &
                         particle_set, kind_of, cell, chia, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, &
                         totalcharge=totalcharge, iounit=iunit)
      END IF

      IF (dft_control%apply_period_efield .OR. dft_control%apply_efield .OR. &
          dft_control%apply_efield_field) THEN
         CALL eeq_efield_energy(qs_env, charges, ef_energy)
         eeq_energy = eeq_energy - SUM(charges*efr)
         DEALLOCATE (efr)
      END IF

      DEALLOCATE (gab, gam, chia)

      CALL timestop(handle)

   END SUBROUTINE xtb_eeq_calculation

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param charges ...
!> \param dcharges ...
!> \param qlagrange ...
!> \param cnumbers ...
!> \param dcnum ...
!> \param eeq_sparam ...
! **************************************************************************************************
   SUBROUTINE xtb_eeq_forces(qs_env, charges, dcharges, qlagrange, cnumbers, dcnum, eeq_sparam)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: charges, dcharges
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: qlagrange
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cnumbers
      TYPE(dcnum_type), DIMENSION(:), INTENT(IN)         :: dcnum
      TYPE(eeq_solver_type), INTENT(IN)                  :: eeq_sparam

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xtb_eeq_forces'

      INTEGER                                            :: atom_a, atom_b, atom_c, enshift_type, &
                                                            handle, i, ia, iatom, ikind, iunit, &
                                                            jatom, jkind, katom, kkind, natom, &
                                                            nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL                                            :: defined, do_ewald, use_virial
      REAL(KIND=dp)                                      :: ala, alb, alpha, cn, ctot, dr, dr2, drk, &
                                                            elag, esg, fe, gam2, gama, grc, kappa, &
                                                            qlam, qq, qq1, qq2, rcut, scn, sgamma, &
                                                            totalcharge, xi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: gam
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: epforce, gab
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, ri, rij, rik, rj
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pvir
      REAL(KIND=dp), DIMENSION(:), POINTER               :: chrgx, dchia, qlag
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_tbe
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

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      IF (do_ewald) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL eeq_solver(qlag, qlam, elag, &
                         particle_set, kind_of, cell, -dcharges, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, &
                         ewald=do_ewald, ewald_env=ewald_env, ewald_pw=ewald_pw, iounit=iunit)
      ELSE
         CALL eeq_solver(qlag, qlam, elag, &
                         particle_set, kind_of, cell, -dcharges, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, iounit=iunit)
      END IF

      enshift_type = xtb_control%enshift_type
      IF (enshift_type == 0) THEN
         enshift_type = 2
         IF (ALL(cell%perd == 0)) enshift_type = 1
      END IF
      sgamma = 8.0_dp ! see D4 for periodic systems paper
      esg = 1.0_dp + EXP(sgamma)
      ALLOCATE (chrgx(natom), dchia(natom))
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, xi=xi, kappa0=kappa)
         !
         ctot = 0.5_dp*(charges(iatom) - qlag(iatom))
         IF (enshift_type == 1) THEN
            scn = SQRT(cnumbers(iatom)) + 1.0e-14_dp
            dchia(iatom) = -ctot*kappa/scn
         ELSE IF (enshift_type == 2) THEN
            cn = cnumbers(iatom)
            scn = 1.0_dp/(esg - cn)
            dchia(iatom) = -ctot*kappa*scn
         ELSE
            CPABORT("Unknown enshift_type")
         END IF
      END DO

      ! Efield
      IF (dft_control%apply_period_efield) THEN
         CALL eeq_efield_force_periodic(qs_env, charges, qlag)
      ELSE IF (dft_control%apply_efield) THEN
         CALL eeq_efield_force_loc(qs_env, charges, qlag)
      ELSE IF (dft_control%apply_efield_field) THEN
         CPABORT("apply field")
      END IF

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
                  END IF
               END IF
            END DO
         END DO
      END DO

      ! Forces from (0.5*q+l)*dA/dR*q
      IF (do_ewald) THEN
         CALL get_qs_env(qs_env, sab_tbe=sab_tbe)
         CALL ewald_env_get(ewald_env, alpha=alpha, rcut=rcut)
         rcut = 1.0_dp*rcut
         CALL neighbor_list_iterator_create(nl_iterator, sab_tbe)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rij)
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            !
            dr2 = SUM(rij**2)
            dr = SQRT(dr2)
            IF (dr > rcut .OR. dr < 1.E-6_dp) CYCLE
            fe = 1.0_dp
            IF (iatom == jatom) fe = 0.5_dp
            qq = (0.5_dp*charges(iatom) - qlag(iatom))*charges(jatom)
            gama = gab(ikind, jkind)
            gam2 = gama*gama
            grc = 2._dp*gama*EXP(-gam2*dr2)*oorootpi/dr - erf(gama*dr)/dr2 &
                  - 2._dp*alpha*EXP(-alpha**2*dr2)*oorootpi/dr + erf(alpha*dr)/dr2
            qq1 = (0.5_dp*charges(iatom) - qlag(iatom))*charges(jatom)
            qq2 = (0.5_dp*charges(jatom) - qlag(jatom))*charges(iatom)
            fdik(:) = -qq1*grc*rij(:)/dr
            force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) + fdik(:)
            force(jkind)%rho_elec(:, atom_b) = force(jkind)%rho_elec(:, atom_b) - fdik(:)
            IF (use_virial) THEN
               CALL virial_pair_force(virial%pv_virial, fe, fdik, rij)
            END IF
            fdik(:) = qq2*grc*rij(:)/dr
            force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) - fdik(:)
            force(jkind)%rho_elec(:, atom_b) = force(jkind)%rho_elec(:, atom_b) + fdik(:)
            IF (use_virial) THEN
               CALL virial_pair_force(virial%pv_virial, -fe, fdik, rij)
            END IF
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)
      ELSE
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
                  fdik(:) = qq*grc*rij(:)/dr
                  force(ikind)%rho_elec(:, atom_a) = force(ikind)%rho_elec(:, atom_a) + fdik(:)
                  force(jkind)%rho_elec(:, atom_b) = force(jkind)%rho_elec(:, atom_b) - fdik(:)
               END DO
            END DO
         END DO
      END IF

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
         END DO
         DEALLOCATE (epforce)

         ! virial
         IF (use_virial) THEN
            chrgx = charges - qlag
            CALL spme_virial(ewald_env, ewald_pw, particle_set, cell, chrgx, pvir)
            virial%pv_virial = virial%pv_virial + pvir
            chrgx = qlag
            CALL spme_virial(ewald_env, ewald_pw, particle_set, cell, chrgx, pvir)
            virial%pv_virial = virial%pv_virial - pvir
         END IF
      END IF

      ! return Lagrange multipliers
      qlagrange(1:natom) = qlag(1:natom)

      DEALLOCATE (gab, chrgx, dchia, qlag)

      CALL timestop(handle)

   END SUBROUTINE xtb_eeq_forces

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dcharges ...
!> \param qlagrange ...
!> \param eeq_sparam ...
! **************************************************************************************************
   SUBROUTINE xtb_eeq_lagrange(qs_env, dcharges, qlagrange, eeq_sparam)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: dcharges
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: qlagrange
      TYPE(eeq_solver_type), INTENT(IN)                  :: eeq_sparam

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xtb_eeq_lagrange'

      INTEGER                                            :: handle, ikind, iunit, jkind, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL                                            :: defined, do_ewald
      REAL(KIND=dp)                                      :: ala, alb, elag, gama, qlam, totalcharge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: gam
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: qlag
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(mp_para_env_type), POINTER                    :: para_env
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
                      dft_control=dft_control)
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)

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

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
      IF (do_ewald) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL eeq_solver(qlag, qlam, elag, &
                         particle_set, kind_of, cell, -dcharges, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, &
                         ewald=do_ewald, ewald_env=ewald_env, ewald_pw=ewald_pw, iounit=iunit)
      ELSE
         CALL eeq_solver(qlag, qlam, elag, &
                         particle_set, kind_of, cell, -dcharges, gam, gab, &
                         para_env, blacs_env, dft_control, eeq_sparam, iounit=iunit)
      END IF

      ! return Lagrange multipliers
      qlagrange(1:natom) = qlag(1:natom)

      DEALLOCATE (gab, qlag)

      CALL timestop(handle)

   END SUBROUTINE xtb_eeq_lagrange

END MODULE xtb_eeq
