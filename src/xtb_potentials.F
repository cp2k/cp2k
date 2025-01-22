!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief xTB (repulsive) pair potentials
!>        Reference: Stefan Grimme, Christoph Bannwarth, Philip Shushkov
!>                   JCTC 13, 1989-2009, (2017)
!>                   DOI: 10.1021/acs.jctc.7b00118
!> \author JGH
! **************************************************************************************************
MODULE xtb_potentials
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE fparser,                         ONLY: evalfd,&
                                              finalizef
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pair_potential,                  ONLY: init_genpot
   USE pair_potential_types,            ONLY: not_initialized,&
                                              pair_potential_p_type,&
                                              pair_potential_pp_create,&
                                              pair_potential_pp_release,&
                                              pair_potential_pp_type,&
                                              pair_potential_single_clean,&
                                              pair_potential_single_copy,&
                                              pair_potential_single_type
   USE pair_potential_util,             ONLY: ener_pot
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
   USE string_utilities,                ONLY: compress,&
                                              uppercase
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   TYPE neighbor_atoms_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: coord => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER        :: rab => NULL()
      INTEGER, DIMENSION(:), POINTER              :: katom => NULL()
   END TYPE neighbor_atoms_type

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_potentials'

   PUBLIC :: xtb_pp_radius, repulsive_potential, srb_potential
   PUBLIC :: nonbonded_correction, xb_interaction
   PUBLIC :: neighbor_atoms_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param erep ...
!> \param kf ...
!> \param enscale ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE repulsive_potential(qs_env, erep, kf, enscale, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(INOUT)                            :: erep
      REAL(dp), INTENT(IN)                               :: kf, enscale
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'repulsive_potential'

      INTEGER                                            :: atom_a, atom_b, handle, iatom, ikind, &
                                                            jatom, jkind, za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: defined, use_virial
      REAL(KIND=dp)                                      :: alphaa, alphab, den2, den4, derepij, dr, &
                                                            ena, enb, ens, erepij, f1, sal, &
                                                            zneffa, zneffb
      REAL(KIND=dp), DIMENSION(3)                        :: force_rr, rij
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_xtb_pp
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b

      CALL timeset(routineN, handle)

      erep = 0._dp

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      atprop=atprop, &
                      sab_xtb_pp=sab_xtb_pp)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)

      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, virial=virial, force=force)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      CALL neighbor_list_iterator_create(nl_iterator, sab_xtb_pp)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=za, xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
         IF (.NOT. defined) CYCLE
         CALL get_qs_kind(qs_kind_set(jkind), zatom=zb, xtb_parameter=xtb_atom_b)
         CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
         IF (.NOT. defined) CYCLE

         dr = SQRT(SUM(rij(:)**2))
         ! repulsive potential
         IF (dr > 0.001_dp) THEN

            ! atomic parameters
            CALL get_xtb_atom_param(xtb_atom_a, en=ena, alpha=alphaa, zneff=zneffa)
            CALL get_xtb_atom_param(xtb_atom_b, en=enb, alpha=alphab, zneff=zneffb)

            ! scaling (not in papers! but in code)
            den2 = (ena - enb)**2
            den4 = den2*den2
            sal = SQRT(alphaa*alphab)
            ens = 1.0_dp + (0.01_dp*den2 + 0.01_dp*den4)*enscale

            erepij = zneffa*zneffb/dr*EXP(-ens*sal*dr**kf)
            erep = erep + erepij
            IF (atprop%energy) THEN
               atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*erepij
               atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*erepij
            END IF
            IF (calculate_forces .AND. (iatom /= jatom .OR. dr > 0.001_dp)) THEN
               derepij = -(1.0_dp/dr + ens*sal*kf*dr**(kf - 1.0_dp))*erepij
               force_rr(1) = derepij*rij(1)/dr
               force_rr(2) = derepij*rij(2)/dr
               force_rr(3) = derepij*rij(3)/dr
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               force(ikind)%repulsive(:, atom_a) = force(ikind)%repulsive(:, atom_a) - force_rr(:)
               force(jkind)%repulsive(:, atom_b) = force(jkind)%repulsive(:, atom_b) + force_rr(:)
               IF (use_virial) THEN
                  f1 = 1.0_dp
                  IF (iatom == jatom) f1 = 0.5_dp
                  CALL virial_pair_force(virial%pv_virial, -f1, force_rr, rij)
               END IF
            END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL timestop(handle)

   END SUBROUTINE repulsive_potential

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param esrb ...
!> \param calculate_forces ...
!> \param xtb_control ...
!> \param cnumbers ...
!> \param dcnum ...
! **************************************************************************************************
   SUBROUTINE srb_potential(qs_env, esrb, calculate_forces, xtb_control, cnumbers, dcnum)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(INOUT)                            :: esrb
      LOGICAL, INTENT(IN)                                :: calculate_forces
      TYPE(xtb_control_type), POINTER                    :: xtb_control
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cnumbers
      TYPE(dcnum_type), DIMENSION(:), INTENT(IN)         :: dcnum

      CHARACTER(len=*), PARAMETER                        :: routineN = 'srb_potential'
      REAL(KIND=dp), DIMENSION(5:9), PARAMETER :: &
         cnfac = (/0.05646607_dp, 0.10514203_dp, 0.09753494_dp, 0.30470380_dp, 0.23261783_dp/), &
         ensrb = (/2.20568300_dp, 2.49640820_dp, 2.81007174_dp, 4.51078438_dp, 4.67476223_dp/), &
         r0srb = (/1.35974851_dp, 0.98310699_dp, 0.98423007_dp, 0.76716063_dp, 1.06139799_dp/)

      INTEGER                                            :: atom_a, atom_b, atom_c, handle, i, &
                                                            iatom, ikind, jatom, jkind, katom, &
                                                            kkind, za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: defined, use_virial
      REAL(KIND=dp)                                      :: c1srb, c2srb, den1, den2, desrbij, dr, &
                                                            dr0, drk, enta, entb, esrbij, etasrb, &
                                                            f1, fhua, fhub, gscal, ksrb, rra0, &
                                                            rrb0, shift
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, force_rr, rij, rik
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_xtb_pp
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b

      CALL timeset(routineN, handle)

      esrb = 0._dp

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      atprop=atprop, &
                      sab_xtb_pp=sab_xtb_pp)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind, &
                               kind_of=kind_of)

      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, virial=virial, force=force)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      ! SRB parameters
      ksrb = xtb_control%ksrb
      etasrb = xtb_control%esrb
      c1srb = xtb_control%c1srb*0.01_dp
      c2srb = xtb_control%c2srb*0.01_dp
      gscal = xtb_control%gscal
      shift = xtb_control%shift

      CALL neighbor_list_iterator_create(nl_iterator, sab_xtb_pp)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, z=za, electronegativity=enta, defined=defined)
         IF (.NOT. defined) CYCLE
         CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
         CALL get_xtb_atom_param(xtb_atom_b, z=zb, electronegativity=entb, defined=defined)
         IF (.NOT. defined) CYCLE

         dr = SQRT(SUM(rij(:)**2))
         ! short-ranged correction term
         IF (dr > 0.001_dp) THEN
         IF (za >= 5 .AND. za <= 9 .AND. zb >= 5 .AND. zb <= 9 .AND. za /= zb) THEN
            rra0 = r0srb(za) + cnfac(za)*cnumbers(iatom) + shift
            rrb0 = r0srb(zb) + cnfac(zb)*cnumbers(jatom) + shift
            den1 = ABS(ensrb(za) - ensrb(zb))
            dr0 = (rra0 + rrb0)*(1._dp - c1srb*den1 - c2srb*den1*den1)
            den2 = (enta - entb)**2
            esrbij = ksrb*EXP(-etasrb*(1._dp + gscal*den2)*(dr - dr0)**2)
            esrb = esrb + esrbij
            IF (atprop%energy) THEN
               atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*esrbij
               atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*esrbij
            END IF
            IF (calculate_forces) THEN
               desrbij = 2.0_dp*esrbij*(-etasrb*(1._dp + gscal*den2)*(dr - dr0))
               force_rr(1) = desrbij*rij(1)/dr
               force_rr(2) = desrbij*rij(2)/dr
               force_rr(3) = desrbij*rij(3)/dr
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               force(ikind)%repulsive(:, atom_a) = force(ikind)%repulsive(:, atom_a) - force_rr(:)
               force(jkind)%repulsive(:, atom_b) = force(jkind)%repulsive(:, atom_b) + force_rr(:)
               IF (use_virial) THEN
                  f1 = 1.0_dp
                  IF (iatom == jatom) f1 = 0.5_dp
                  CALL virial_pair_force(virial%pv_virial, -f1, force_rr, rij)
               END IF
               ! coordination number derivatives
               ! iatom
               fhua = -desrbij*cnfac(za)*(1._dp - c1srb*den1 - c2srb*den1*den1)
               DO i = 1, dcnum(iatom)%neighbors
                  katom = dcnum(iatom)%nlist(i)
                  kkind = kind_of(katom)
                  atom_c = atom_of_kind(katom)
                  rik = dcnum(iatom)%rik(:, i)
                  drk = SQRT(SUM(rik(:)**2))
                  IF (drk > 1.e-3_dp) THEN
                     fdik(:) = fhua*dcnum(iatom)%dvals(i)*rik(:)/drk
                     force(ikind)%repulsive(:, atom_a) = force(ikind)%repulsive(:, atom_a) - fdik(:)
                     force(kkind)%repulsive(:, atom_c) = force(kkind)%repulsive(:, atom_c) + fdik(:)
                     IF (use_virial) THEN
                        CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                     END IF
                  END IF
               END DO
               ! jatom
               fhub = -desrbij*cnfac(zb)*(1._dp - c1srb*den1 - c2srb*den1*den1)
               DO i = 1, dcnum(jatom)%neighbors
                  katom = dcnum(jatom)%nlist(i)
                  kkind = kind_of(katom)
                  atom_c = atom_of_kind(katom)
                  rik = dcnum(jatom)%rik(:, i)
                  drk = SQRT(SUM(rik(:)**2))
                  IF (drk > 1.e-3_dp) THEN
                     fdik(:) = fhub*dcnum(jatom)%dvals(i)*rik(:)/drk
                     force(jkind)%repulsive(:, atom_b) = force(jkind)%repulsive(:, atom_b) - fdik(:)
                     force(kkind)%repulsive(:, atom_c) = force(kkind)%repulsive(:, atom_c) + fdik(:)
                     IF (use_virial) THEN
                        CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                     END IF
                  END IF
               END DO
            END IF
         END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL timestop(handle)

   END SUBROUTINE srb_potential

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param ppradius ...
!> \param eps_pair ...
!> \param kfparam ...
! **************************************************************************************************
   SUBROUTINE xtb_pp_radius(qs_kind_set, ppradius, eps_pair, kfparam)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: ppradius
      REAL(KIND=dp), INTENT(IN)                          :: eps_pair, kfparam

      INTEGER                                            :: ikind, ir, jkind, nkind
      LOGICAL                                            :: defa, defb
      REAL(KIND=dp)                                      :: alphaa, alphab, erep, rab, rab0, rcova, &
                                                            rcovb, saa, zneffa, zneffb
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b

      ppradius = 0.0_dp
      nkind = SIZE(ppradius, 1)
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, rcov=rcova, alpha=alphaa, zneff=zneffa, defined=defa)
         IF (.NOT. defa) CYCLE
         DO jkind = ikind, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, rcov=rcovb, alpha=alphab, zneff=zneffb, defined=defb)
            IF (.NOT. defb) CYCLE
            rab = 0.0_dp
            DO ir = 1, 24
               rab = rab + 1.0_dp
               saa = SQRT(alphaa*alphab)
               erep = zneffa*zneffb/rab*EXP(-saa*rab**kfparam)
               IF (erep < eps_pair) EXIT
            END DO
            rab0 = rcova + rcovb
            rab = MAX(rab, rab0 + 2.0_dp)
            ppradius(ikind, jkind) = rab
            ppradius(jkind, ikind) = ppradius(ikind, jkind)
         END DO
      END DO

   END SUBROUTINE xtb_pp_radius

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param exb ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE xb_interaction(qs_env, exb, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(INOUT)                       :: exb
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'xb_interaction'

      INTEGER                                            :: atom_a, atom_b, handle, iatom, ikind, &
                                                            jatom, jkind, na, natom, nkind, zat
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: defined, use_virial
      REAL(KIND=dp)                                      :: dr, kx2, kxr, rcova, rcovb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: kx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rcab
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_atoms_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: neighbor_atoms
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_xb, sab_xtb_pp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      sab_xb=sab_xb, &
                      sab_xtb_pp=sab_xtb_pp)

      nkind = SIZE(atomic_kind_set)
      xtb_control => dft_control%qs_control%xtb_control

      ! global parameters
      kxr = xtb_control%kxr
      kx2 = xtb_control%kx2

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, virial=virial, force=force)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      ! list of neighbor atoms for XB term
      ALLOCATE (neighbor_atoms(nkind))
      DO ikind = 1, nkind
         NULLIFY (neighbor_atoms(ikind)%coord)
         NULLIFY (neighbor_atoms(ikind)%rab)
         NULLIFY (neighbor_atoms(ikind)%katom)
         CALL get_atomic_kind(atomic_kind_set(ikind), z=zat, natom=na)
         IF (zat == 17 .OR. zat == 35 .OR. zat == 53 .OR. zat == 85) THEN
            ALLOCATE (neighbor_atoms(ikind)%coord(3, na))
            neighbor_atoms(ikind)%coord(1:3, 1:na) = 0.0_dp
            ALLOCATE (neighbor_atoms(ikind)%rab(na))
            neighbor_atoms(ikind)%rab(1:na) = HUGE(0.0_dp)
            ALLOCATE (neighbor_atoms(ikind)%katom(na))
            neighbor_atoms(ikind)%katom(1:na) = 0
         END IF
      END DO
      ! kx parameters
      ALLOCATE (kx(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, kx=kx(ikind))
      END DO
      !
      ALLOCATE (rcab(nkind, nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, rcov=rcova)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, rcov=rcovb)
            rcab(ikind, jkind) = kxr*(rcova + rcovb)
         END DO
      END DO

      CALL neighbor_list_iterator_create(nl_iterator, sab_xtb_pp)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined)
         IF (.NOT. defined) CYCLE
         CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
         CALL get_xtb_atom_param(xtb_atom_b, defined=defined)
         IF (.NOT. defined) CYCLE

         dr = SQRT(SUM(rij(:)**2))

         ! neighbor atom for XB term
         IF (dr > 1.e-3_dp) THEN
            IF (ASSOCIATED(neighbor_atoms(ikind)%rab)) THEN
               atom_a = atom_of_kind(iatom)
               IF (dr < neighbor_atoms(ikind)%rab(atom_a)) THEN
                  neighbor_atoms(ikind)%rab(atom_a) = dr
                  neighbor_atoms(ikind)%coord(1:3, atom_a) = rij(1:3)
                  neighbor_atoms(ikind)%katom(atom_a) = jatom
               END IF
            END IF
            IF (ASSOCIATED(neighbor_atoms(jkind)%rab)) THEN
               atom_b = atom_of_kind(jatom)
               IF (dr < neighbor_atoms(jkind)%rab(atom_b)) THEN
                  neighbor_atoms(jkind)%rab(atom_b) = dr
                  neighbor_atoms(jkind)%coord(1:3, atom_b) = -rij(1:3)
                  neighbor_atoms(jkind)%katom(atom_b) = iatom
               END IF
            END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      exb = 0.0_dp
      CALL xb_neighbors(neighbor_atoms, para_env)
      CALL xb_energy(exb, neighbor_atoms, atom_of_kind, kind_of, sab_xb, kx, kx2, rcab, &
                     calculate_forces, use_virial, force, virial, atprop)

      DO ikind = 1, nkind
         IF (ASSOCIATED(neighbor_atoms(ikind)%coord)) THEN
            DEALLOCATE (neighbor_atoms(ikind)%coord)
         END IF
         IF (ASSOCIATED(neighbor_atoms(ikind)%rab)) THEN
            DEALLOCATE (neighbor_atoms(ikind)%rab)
         END IF
         IF (ASSOCIATED(neighbor_atoms(ikind)%katom)) THEN
            DEALLOCATE (neighbor_atoms(ikind)%katom)
         END IF
      END DO
      DEALLOCATE (neighbor_atoms)
      DEALLOCATE (kx, rcab)

      CALL timestop(handle)

   END SUBROUTINE xb_interaction

! **************************************************************************************************
!> \brief  Distributes the neighbor atom information to all processors
!>
!> \param neighbor_atoms ...
!> \param para_env ...
!> \par History
!>       1.2019 JGH
!> \version 1.1
! **************************************************************************************************
   SUBROUTINE xb_neighbors(neighbor_atoms, para_env)
      TYPE(neighbor_atoms_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: neighbor_atoms
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: iatom, ikind, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: matom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dmloc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: coord

      nkind = SIZE(neighbor_atoms)
      DO ikind = 1, nkind
         IF (ASSOCIATED(neighbor_atoms(ikind)%rab)) THEN
            natom = SIZE(neighbor_atoms(ikind)%rab)
            ALLOCATE (dmloc(2*natom), matom(natom), coord(3, natom))
            dmloc = 0.0_dp
            DO iatom = 1, natom
               dmloc(2*iatom - 1) = neighbor_atoms(ikind)%rab(iatom)
               dmloc(2*iatom) = REAL(para_env%mepos, KIND=dp)
            END DO
            CALL para_env%minloc(dmloc)
            coord = 0.0_dp
            matom = 0
            DO iatom = 1, natom
               neighbor_atoms(ikind)%rab(iatom) = dmloc(2*iatom - 1)
               IF (NINT(dmloc(2*iatom)) == para_env%mepos) THEN
                  coord(1:3, iatom) = neighbor_atoms(ikind)%coord(1:3, iatom)
                  matom(iatom) = neighbor_atoms(ikind)%katom(iatom)
               END IF
            END DO
            CALL para_env%sum(coord)
            neighbor_atoms(ikind)%coord(1:3, :) = coord(1:3, :)
            CALL para_env%sum(matom)
            neighbor_atoms(ikind)%katom(:) = matom(:)
            DEALLOCATE (dmloc, matom, coord)
         END IF
      END DO

   END SUBROUTINE xb_neighbors

! **************************************************************************************************
!> \brief  Computes a correction for nonbonded interactions based on a generic potential
!>
!> \param enonbonded energy contribution
!> \param force ...
!> \param qs_env ...
!> \param xtb_control ...
!> \param sab_xtb_nonbond ...
!> \param atomic_kind_set ...
!> \param calculate_forces ...
!> \param use_virial ...
!> \param virial ...
!> \param atprop ...
!> \param atom_of_kind ..
!> \par History
!>      12.2018 JGH
!> \version 1.1
! **************************************************************************************************
   SUBROUTINE nonbonded_correction(enonbonded, force, qs_env, xtb_control, sab_xtb_nonbond, &
                                   atomic_kind_set, calculate_forces, use_virial, virial, atprop, atom_of_kind)
      REAL(dp), INTENT(INOUT)                            :: enonbonded
      TYPE(qs_force_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: force
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(xtb_control_type), POINTER                    :: xtb_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(IN), POINTER                             :: sab_xtb_nonbond
      TYPE(atomic_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: atomic_kind_set
      LOGICAL, INTENT(IN)                                :: calculate_forces, use_virial
      TYPE(virial_type), INTENT(IN), POINTER             :: virial
      TYPE(atprop_type), INTENT(IN), POINTER             :: atprop
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind

      CHARACTER(len=*), PARAMETER :: routineN = 'nonbonded_correction'

      CHARACTER(LEN=default_string_length)               :: def_error, this_error
      INTEGER                                            :: atom_i, atom_j, handle, iatom, ikind, &
                                                            jatom, jkind, kk, ntype
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: do_ewald
      REAL(KIND=dp)                                      :: dedf, dr, dx, energy_cutoff, err, fval, &
                                                            lerr, rcut
      REAL(KIND=dp), DIMENSION(3)                        :: fij, rij
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(pair_potential_p_type), POINTER               :: nonbonded
      TYPE(pair_potential_pp_type), POINTER              :: potparm
      TYPE(pair_potential_single_type), POINTER          :: pot
      TYPE(section_vals_type), POINTER                   :: nonbonded_section

      CALL timeset(routineN, handle)

      NULLIFY (nonbonded)
      NULLIFY (potparm)
      NULLIFY (ewald_env)
      nonbonded => xtb_control%nonbonded
      do_ewald = xtb_control%do_ewald
      CALL get_qs_env(qs_env=qs_env, ewald_env=ewald_env)

      ntype = SIZE(atomic_kind_set)
      CALL pair_potential_pp_create(potparm, ntype)
      !Assign input and potential info to potparm_nonbond
      CALL force_field_pack_nonbond_pot_correction(atomic_kind_set, nonbonded, potparm, ewald_env, do_ewald)
      !Initialize genetic potential
      CALL init_genpot(potparm, ntype)

      NULLIFY (pot)
      enonbonded = 0._dp
      energy_cutoff = 0._dp

      CALL neighbor_list_iterator_create(nl_iterator, sab_xtb_nonbond)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         pot => potparm%pot(ikind, jkind)%pot
         dr = SQRT(rij(1)**2 + rij(2)**2 + rij(3)**2)
         rcut = SQRT(pot%rcutsq)
         IF (dr <= rcut .AND. dr > 1.E-3_dp) THEN
            fval = 1.0_dp
            IF (ikind == jkind) fval = 0.5_dp
            ! splines not implemented
            enonbonded = enonbonded + fval*ener_pot(pot, dr, energy_cutoff)
            IF (atprop%energy) THEN
               atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*fval*ener_pot(pot, dr, energy_cutoff)
               atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*fval*ener_pot(pot, dr, energy_cutoff)
            END IF
         END IF

         IF (calculate_forces) THEN

            kk = SIZE(pot%type)
            IF (kk /= 1) THEN
               CALL cp_warn(__LOCATION__, "Generic potential with type > 1 not implemented.")
               CPABORT("pot type")
            END IF
            ! rmin and rmax and rcut
            IF ((pot%set(kk)%rmin /= not_initialized) .AND. (dr < pot%set(kk)%rmin)) CYCLE
            ! An upper boundary for the potential definition was defined
            IF ((pot%set(kk)%rmax /= not_initialized) .AND. (dr >= pot%set(kk)%rmax)) CYCLE
            ! If within limits let's compute the potential...
            IF (dr <= rcut .AND. dr > 1.E-3_dp) THEN

               NULLIFY (nonbonded_section)
               nonbonded_section => section_vals_get_subs_vals(qs_env%input, "DFT%QS%xTB%NONBONDED")
               CALL section_vals_val_get(nonbonded_section, "DX", r_val=dx)
               CALL section_vals_val_get(nonbonded_section, "ERROR_LIMIT", r_val=lerr)

               dedf = fval*evalfd(pot%set(kk)%gp%myid, 1, pot%set(kk)%gp%values, dx, err)
               IF (ABS(err) > lerr) THEN
                  WRITE (this_error, "(A,G12.6,A)") "(", err, ")"
                  WRITE (def_error, "(A,G12.6,A)") "(", lerr, ")"
                  CALL compress(this_error, .TRUE.)
                  CALL compress(def_error, .TRUE.)
                  CALL cp_warn(__LOCATION__, &
                               'ASSERTION (cond) failed at line '//cp_to_string(__LINE__)// &
                               ' Error '//TRIM(this_error)//' in computing numerical derivatives larger then'// &
                               TRIM(def_error)//' .')
               END IF

               atom_i = atom_of_kind(iatom)
               atom_j = atom_of_kind(jatom)

               fij(1:3) = dedf*rij(1:3)/pot%set(kk)%gp%values(1)

               force(ikind)%repulsive(:, atom_i) = force(ikind)%repulsive(:, atom_i) - fij(:)
               force(jkind)%repulsive(:, atom_j) = force(jkind)%repulsive(:, atom_j) + fij(:)

               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, -1._dp, fij, rij)
               END IF

            END IF
         END IF
         NULLIFY (pot)
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CALL finalizef()

      IF (ASSOCIATED(potparm)) THEN
         CALL pair_potential_pp_release(potparm)
      END IF

      CALL timestop(handle)

   END SUBROUTINE nonbonded_correction

! **************************************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param nonbonded ...
!> \param potparm ...
!> \param ewald_env ...
!> \param do_ewald ...
! **************************************************************************************************
   SUBROUTINE force_field_pack_nonbond_pot_correction(atomic_kind_set, nonbonded, potparm, ewald_env, do_ewald)

      ! routine based on force_field_pack_nonbond
      TYPE(atomic_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: atomic_kind_set
      TYPE(pair_potential_p_type), INTENT(IN), POINTER   :: nonbonded
      TYPE(pair_potential_pp_type), INTENT(INOUT), &
         POINTER                                         :: potparm
      TYPE(ewald_environment_type), INTENT(IN), POINTER  :: ewald_env
      LOGICAL, INTENT(IN)                                :: do_ewald

      CHARACTER(LEN=default_string_length)               :: name_atm_a, name_atm_a_local, &
                                                            name_atm_b, name_atm_b_local
      INTEGER                                            :: ikind, ingp, iw, jkind
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: ewald_rcut
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(pair_potential_single_type), POINTER          :: pot

      NULLIFY (pot, logger)

      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)

      DO ikind = 1, SIZE(atomic_kind_set)
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, name=name_atm_a_local)
         DO jkind = ikind, SIZE(atomic_kind_set)
            atomic_kind => atomic_kind_set(jkind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, name=name_atm_b_local)
            found = .FALSE.
            name_atm_a = name_atm_a_local
            name_atm_b = name_atm_b_local
            CALL uppercase(name_atm_a)
            CALL uppercase(name_atm_b)
            pot => potparm%pot(ikind, jkind)%pot
            IF (ASSOCIATED(nonbonded)) THEN
               DO ingp = 1, SIZE(nonbonded%pot)
                  IF ((TRIM(nonbonded%pot(ingp)%pot%at1) == "*") .OR. &
                      (TRIM(nonbonded%pot(ingp)%pot%at2) == "*")) CYCLE

                  !IF (iw > 0) WRITE (iw, *) "TESTING ", TRIM(name_atm_a), TRIM(name_atm_b), &
                  !   " with ", TRIM(nonbonded%pot(ingp)%pot%at1), &
                  !   TRIM(nonbonded%pot(ingp)%pot%at2)
                  IF ((((name_atm_a) == (nonbonded%pot(ingp)%pot%at1)) .AND. &
                       ((name_atm_b) == (nonbonded%pot(ingp)%pot%at2))) .OR. &
                      (((name_atm_b) == (nonbonded%pot(ingp)%pot%at1)) .AND. &
                       ((name_atm_a) == (nonbonded%pot(ingp)%pot%at2)))) THEN
                     CALL pair_potential_single_copy(nonbonded%pot(ingp)%pot, pot)
                     ! multiple potential not implemented, simply overwriting
                     IF (found) &
                        CALL cp_warn(__LOCATION__, &
                                     "Multiple NONBONDED declaration: "//TRIM(name_atm_a)// &
                                     " and "//TRIM(name_atm_b)//" OVERWRITING! ")
                     !IF (iw > 0) WRITE (iw, *) "    FOUND ", TRIM(name_atm_a), " ", TRIM(name_atm_b)
                     found = .TRUE.
                  END IF
               END DO
            END IF
            IF (.NOT. found) THEN
               CALL pair_potential_single_clean(pot)
               !IF (iw > 0) WRITE (iw, *) " NOTFOUND ", TRIM(name_atm_a), " ", TRIM(name_atm_b)
            END IF
         END DO !jkind
      END DO !ikind

      ! Cutoff is defined always as the maximum between the FF and Ewald
      IF (do_ewald) THEN
         CALL ewald_env_get(ewald_env, rcut=ewald_rcut)
         pot%rcutsq = MAX(pot%rcutsq, ewald_rcut*ewald_rcut)
         !IF (iw > 0) WRITE (iw, *) " RCUT     ", SQRT(pot%rcutsq), ewald_rcut
      END IF

   END SUBROUTINE force_field_pack_nonbond_pot_correction

! **************************************************************************************************
!> \brief  Computes the interaction term between Br/I/At and donor atoms
!>
!> \param exb ...
!> \param neighbor_atoms ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param sab_xb ...
!> \param kx ...
!> \param kx2 ...
!> \param rcab ...
!> \param calculate_forces ...
!> \param use_virial ...
!> \param force ...
!> \param virial ...
!> \param atprop ...
!> \par History
!>      12.2018 JGH
!> \version 1.1
! **************************************************************************************************
   SUBROUTINE xb_energy(exb, neighbor_atoms, atom_of_kind, kind_of, sab_xb, kx, kx2, rcab, &
                        calculate_forces, use_virial, force, virial, atprop)
      REAL(dp), INTENT(INOUT)                            :: exb
      TYPE(neighbor_atoms_type), DIMENSION(:), &
         INTENT(IN)                                      :: neighbor_atoms
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_xb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kx
      REAL(dp), INTENT(IN)                               :: kx2
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: rcab
      LOGICAL, INTENT(IN)                                :: calculate_forces, use_virial
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial
      TYPE(atprop_type), POINTER                         :: atprop

      INTEGER                                            :: atom_a, atom_b, atom_c, iatom, ikind, &
                                                            jatom, jkind, katom, kkind
      INTEGER, DIMENSION(3)                              :: cell
      REAL(KIND=dp)                                      :: alp, aterm, cosa, daterm, ddab, ddax, &
                                                            ddbx, ddr, ddr12, ddr6, deval, dr, &
                                                            drab, drax, drbx, eval, xy
      REAL(KIND=dp), DIMENSION(3)                        :: fia, fij, fja, ria, rij, rja
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      ! exonent in angular term
      alp = 6.0_dp
      ! loop over all atom pairs
      CALL neighbor_list_iterator_create(nl_iterator, sab_xb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         ! ikind, iatom : Halogen
         ! jkind, jatom : Donor
         atom_a = atom_of_kind(iatom)
         katom = neighbor_atoms(ikind)%katom(atom_a)
         IF (katom == 0) CYCLE
         dr = SQRT(rij(1)**2 + rij(2)**2 + rij(3)**2)
         ddr = rcab(ikind, jkind)/dr
         ddr6 = ddr**6
         ddr12 = ddr6*ddr6
         eval = kx(ikind)*(ddr12 - kx2*ddr6)/(1.0_dp + ddr12)
         ! angle
         ria(1:3) = neighbor_atoms(ikind)%coord(1:3, atom_a)
         rja(1:3) = rij(1:3) - ria(1:3)
         drax = ria(1)**2 + ria(2)**2 + ria(3)**2
         drbx = dr*dr
         drab = rja(1)**2 + rja(2)**2 + rja(3)**2
         xy = SQRT(drbx*drax)
         ! cos angle B-X-A
         cosa = (drbx + drax - drab)/xy
         aterm = (0.5_dp - 0.25_dp*cosa)**alp
         !
         exb = exb + aterm*eval
         IF (atprop%energy) THEN
            atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*aterm*eval
            atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*aterm*eval
         END IF
         !
         IF (calculate_forces) THEN
            kkind = kind_of(katom)
            atom_b = atom_of_kind(jatom)
            atom_c = atom_of_kind(katom)
            !
            deval = 6.0_dp*kx(ikind)*ddr6*(kx2*ddr12 + 2.0_dp*ddr6 - kx2)/(1.0_dp + ddr12)**2
            deval = -rcab(ikind, jkind)*deval/(dr*dr)/ddr
            fij(1:3) = aterm*deval*rij(1:3)/dr
            force(ikind)%repulsive(:, atom_a) = force(ikind)%repulsive(:, atom_a) - fij(:)
            force(jkind)%repulsive(:, atom_b) = force(jkind)%repulsive(:, atom_b) + fij(:)
            IF (use_virial) THEN
               CALL virial_pair_force(virial%pv_virial, -1._dp, fij, rij)
            END IF
            !
            fij(1:3) = 0.0_dp
            fia(1:3) = 0.0_dp
            fja(1:3) = 0.0_dp
            daterm = -0.25_dp*alp*(0.5_dp - 0.25_dp*cosa)**(alp - 1.0_dp)
            ddbx = 0.5_dp*(drab - drax + drbx)/xy/drbx
            ddax = 0.5_dp*(drab + drax - drbx)/xy/drax
            ddab = -1._dp/xy
            fij(1:3) = 2.0_dp*daterm*ddbx*rij(1:3)*eval
            fia(1:3) = 2.0_dp*daterm*ddax*ria(1:3)*eval
            fja(1:3) = 2.0_dp*daterm*ddab*rja(1:3)*eval
            force(ikind)%repulsive(:, atom_a) = force(ikind)%repulsive(:, atom_a) - fij(:) - fia(:)
            force(jkind)%repulsive(:, atom_b) = force(jkind)%repulsive(:, atom_b) + fij(:) + fja(:)
            force(kkind)%repulsive(:, atom_c) = force(kkind)%repulsive(:, atom_c) + fia(:) - fja(:)
            IF (use_virial) THEN
               CALL virial_pair_force(virial%pv_virial, -1._dp, fij, rij)
               CALL virial_pair_force(virial%pv_virial, -1._dp, fia, ria)
               CALL virial_pair_force(virial%pv_virial, -1._dp, fja, rja)
            END IF
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

   END SUBROUTINE xb_energy

END MODULE xtb_potentials
