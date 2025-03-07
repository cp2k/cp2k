!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of D2 dispersion
!> \author JGH
! **************************************************************************************************
MODULE qs_dispersion_d2

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: bohr,&
                                              kjmol
   USE qs_dispersion_types,             ONLY: qs_atom_dispersion_type,&
                                              qs_dispersion_type
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
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_d2'

   PUBLIC :: calculate_dispersion_d2_pairpot, dftd2_param

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param evdw ...
!> \param calculate_forces ...
!> \param atevdw ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_d2_pairpot(qs_env, dispersion_env, evdw, calculate_forces, atevdw)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), INTENT(OUT)                         :: evdw
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: atevdw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_d2_pairpot'

      INTEGER                                            :: atom_a, atom_b, handle, iatom, ikind, &
                                                            jatom, jkind, mepos, natom, nkind, &
                                                            num_pe, za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, atomnumber, kind_of
      LOGICAL                                            :: atenergy, atex, floating_a, ghost_a, &
                                                            use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: dodisp, floating, ghost
      REAL(KIND=dp)                                      :: c6, dd, devdw, dfdmp, dr, er, fac, fdmp, &
                                                            rcc, rcut, s6, xp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c6d2, radd2
      REAL(KIND=dp), DIMENSION(3)                        :: fdij, rij
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_virial_thread
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atener
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_vdw
      TYPE(qs_atom_dispersion_type), POINTER             :: disp_a
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      evdw = 0._dp

      NULLIFY (atomic_kind_set, qs_kind_set, sab_vdw)

      CALL get_qs_env(qs_env=qs_env, nkind=nkind, natom=natom, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, cell=cell, virial=virial, atprop=atprop)

      ! atomic energy and stress arrays
      atenergy = atprop%energy
      IF (atenergy) THEN
         CALL atprop_array_init(atprop%atevdw, natom)
         atener => atprop%atevdw
      END IF
      ! external atomic energy
      atex = .FALSE.
      IF (PRESENT(atevdw)) THEN
         atex = .TRUE.
      END IF

      NULLIFY (force)
      CALL get_qs_env(qs_env=qs_env, force=force)
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      pv_virial_thread(:, :) = 0._dp

      ALLOCATE (dodisp(nkind), ghost(nkind), floating(nkind), atomnumber(nkind), c6d2(nkind), radd2(nkind))
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
         CALL get_qs_kind(qs_kind_set(ikind), dispersion=disp_a, ghost=ghost_a, floating=floating_a)
         dodisp(ikind) = disp_a%defined
         ghost(ikind) = ghost_a
         floating(ikind) = floating_a
         atomnumber(ikind) = za
         c6d2(ikind) = disp_a%c6
         radd2(ikind) = disp_a%vdw_radii
      END DO

      rcut = 2._dp*dispersion_env%rc_disp
      s6 = dispersion_env%scaling
      dd = dispersion_env%exp_pre

      sab_vdw => dispersion_env%sab_vdw
      num_pe = 1
      CALL neighbor_list_iterator_create(nl_iterator, sab_vdw, nthread=num_pe)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom, r=rij)

         IF (ghost(ikind) .OR. ghost(jkind) .OR. floating(ikind) .OR. floating(jkind)) CYCLE

         IF (.NOT. (dodisp(ikind) .AND. dodisp(jkind))) CYCLE

         za = atomnumber(ikind)
         zb = atomnumber(jkind)
         ! vdW potential
         dr = SQRT(SUM(rij(:)**2))
         IF (dr <= rcut) THEN
            fac = 1._dp
            IF (iatom == jatom) fac = 0.5_dp
            IF (dr > 0.001_dp) THEN
               c6 = SQRT(c6d2(ikind)*c6d2(jkind))
               rcc = radd2(ikind) + radd2(jkind)
               er = EXP(-dd*(dr/rcc - 1._dp))
               fdmp = 1._dp/(1._dp + er)
               xp = s6*c6/dr**6
               evdw = evdw - xp*fdmp*fac
               IF (calculate_forces) THEN
                  dfdmp = dd/rcc*er*fdmp*fdmp
                  devdw = -xp*(-6._dp*fdmp/dr + dfdmp)
                  fdij(:) = devdw*rij(:)/dr*fac
                  atom_a = atom_of_kind(iatom)
                  atom_b = atom_of_kind(jatom)
                  force(ikind)%dispersion(:, atom_a) = force(ikind)%dispersion(:, atom_a) - fdij(:)
                  force(jkind)%dispersion(:, atom_b) = force(jkind)%dispersion(:, atom_b) + fdij(:)
                  IF (use_virial) THEN
                     CALL virial_pair_force(pv_virial_thread, -1._dp, fdij, rij)
                  END IF
               END IF
               IF (atenergy) THEN
                  atener(iatom) = atener(iatom) - 0.5_dp*xp*fdmp*fac
                  atener(jatom) = atener(jatom) - 0.5_dp*xp*fdmp*fac
               END IF
               IF (atex) THEN
                  atevdw(iatom) = atevdw(iatom) - 0.5_dp*xp*fdmp*fac
                  atevdw(jatom) = atevdw(jatom) - 0.5_dp*xp*fdmp*fac
               END IF
            END IF
         END IF

      END DO

      virial%pv_virial = virial%pv_virial + pv_virial_thread

      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (dodisp, ghost, floating, atomnumber, radd2, c6d2)

      CALL timestop(handle)

   END SUBROUTINE calculate_dispersion_d2_pairpot

! **************************************************************************************************
!> \brief ...
!> \param z ...
!> \param c6 ...
!> \param r ...
!> \param found ...
! **************************************************************************************************
   SUBROUTINE dftd2_param(z, c6, r, found)

      INTEGER, INTENT(in)                                :: z
      REAL(KIND=dp), INTENT(inout)                       :: c6, r
      LOGICAL, INTENT(inout)                             :: found

      REAL(KIND=dp), DIMENSION(54), PARAMETER :: c6val = (/0.14_dp, 0.08_dp, 1.61_dp, 1.61_dp, &
         3.13_dp, 1.75_dp, 1.23_dp, 0.70_dp, 0.75_dp, 0.63_dp, 5.71_dp, 5.71_dp, 10.79_dp, 9.23_dp,&
         7.84_dp, 5.57_dp, 5.07_dp, 4.61_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, &
         10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 10.80_dp, 16.99_dp, 17.10_dp, &
         16.37_dp, 12.64_dp, 12.47_dp, 12.01_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, &
         24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 24.67_dp, 37.32_dp, 38.71_dp, &
         38.44_dp, 31.74_dp, 31.50_dp, 29.99_dp/)
      REAL(KIND=dp), DIMENSION(54), PARAMETER :: rval = (/1.001_dp, 1.012_dp, 0.825_dp, 1.408_dp, &
         1.485_dp, 1.452_dp, 1.397_dp, 1.342_dp, 1.287_dp, 1.243_dp, 1.144_dp, 1.364_dp, 1.639_dp, &
         1.716_dp, 1.705_dp, 1.683_dp, 1.639_dp, 1.595_dp, 1.485_dp, 1.474_dp, 1.562_dp, 1.562_dp, &
         1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.562_dp, 1.650_dp, &
         1.727_dp, 1.760_dp, 1.771_dp, 1.749_dp, 1.727_dp, 1.628_dp, 1.606_dp, 1.639_dp, 1.639_dp, &
         1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.639_dp, 1.672_dp, &
         1.804_dp, 1.881_dp, 1.892_dp, 1.892_dp, 1.881_dp/)

!
! GRIMME DISPERSION PARAMETERS
! Stefan Grimme, Semiempirical GGA-Type Density Functional Constructed
!                with a Long-Range Dispersion Correction, J. Comp. Chem. 27: 1787-1799 (2006)
! doi:10.1002/jcc.20495
!
! Conversion factor [Jnm^6mol^-1] -> [a.u.] : 17.34527758021901
! Conversion factor [A] -> [a.u.] : 1.889726132885643
!
! C6 values in [Jnm^6/mol]
! vdW radii [A]

      IF (z > 0 .AND. z <= 54) THEN
         found = .TRUE.
         c6 = c6val(z)*1000._dp*bohr**6/kjmol
         r = rval(z)*bohr
      ELSE
         found = .FALSE.
      END IF

   END SUBROUTINE dftd2_param

! **************************************************************************************************

END MODULE qs_dispersion_d2
