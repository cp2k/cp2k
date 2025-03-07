!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of D3 dispersion
!> \author JGH
! **************************************************************************************************
MODULE qs_dispersion_d3

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc,&
                                              plane_distance
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE input_constants,                 ONLY: vdw_pairpot_dftd3,&
                                              vdw_pairpot_dftd3bj
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: kcalmol
   USE qs_dispersion_cnum,              ONLY: d3_cnumber,&
                                              dcnum_distribute,&
                                              dcnum_type,&
                                              exclude_d3_kind_pair
   USE qs_dispersion_types,             ONLY: dftd3_pp,&
                                              qs_atom_dispersion_type,&
                                              qs_dispersion_type
   USE qs_dispersion_utils,             ONLY: cellhash
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_d3'

   PUBLIC :: calculate_dispersion_d3_pairpot, dftd3_c6_param

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param evdw ...
!> \param calculate_forces ...
!> \param unit_nr ...
!> \param atevdw ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_d3_pairpot(qs_env, dispersion_env, evdw, calculate_forces, &
                                              unit_nr, atevdw)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), INTENT(OUT)                         :: evdw
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: atevdw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_d3_pairpot'

      INTEGER :: atom_a, atom_b, atom_c, atom_d, handle, hashb, hashc, i, ia, iat, iatom, icx, &
         icy, icz, idmp, ikind, ilist, imol, jatom, jkind, katom, kkind, kstart, latom, lkind, &
         max_elem, maxc, mepos, na, natom, nb, nc, nkind, num_pe, za, zb, zc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, atomnumber, kind_of
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c, ncell, periodic
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL :: atenergy, atex, debugall, domol, exclude_pair, floating_a, floating_b, &
         floating_c, ghost_a, ghost_b, ghost_c, is000, use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: dodisp, exclude
      REAL(KIND=dp) :: a1, a2, alp6, alp8, ang, c6, c8, c9, cc6ab, cc6bc, cc6ca, cnum, dc6a, dc6b, &
         dc8a, dc8b, dcc6aba, dcc6abb, dcc6bcb, dcc6bcc, dcc6caa, dcc6cac, de6, de8, de91, de921, &
         de922, dea, dfdab6, dfdab8, dfdabc, dr, drk, e6, e6tot, e8, e8tot, e9, e9tot, elrc6, &
         elrc8, elrc9, eps_cn, esrb, f0ab, fac, fac0, fdab6, fdab8, fdabc, gsrb, kgc8, nab, nabc, &
         r0, r0ab, r2ab, r2abc, r2bc, r2ca, r6, r8, rabc, rc2, rcc, rcut, rs6, rs8, s1, s2, s3, &
         s6, s8, s8i, s9, srbe, ssrb, t1srb, t2srb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: atom2mol, c6d2, cnkind, cnumbers, &
                                                            cnumfix, radd2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rcpbc
      REAL(KIND=dp), DIMENSION(3)                        :: fdij, fdik, ra, rab, rb, rb0, rbc, rc, &
                                                            rc0, rca, rij, rik, sab_max
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_virial_thread
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atener
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_vdw
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_atom_dispersion_type), POINTER             :: disp_a, disp_b, disp_c
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      evdw = 0._dp

      NULLIFY (atomic_kind_set, qs_kind_set, sab_vdw)

      CALL get_qs_env(qs_env=qs_env, nkind=nkind, natom=natom, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, cell=cell, virial=virial, para_env=para_env, atprop=atprop)

      debugall = dispersion_env%verbose

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

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)

      NULLIFY (force)
      CALL get_qs_env(qs_env=qs_env, force=force)
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      pv_virial_thread(:, :) = 0._dp

      ALLOCATE (dodisp(nkind), exclude(nkind), atomnumber(nkind), c6d2(nkind), radd2(nkind))
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
         CALL get_qs_kind(qs_kind_set(ikind), dispersion=disp_a, ghost=ghost_a, floating=floating_a)
         dodisp(ikind) = disp_a%defined
         exclude(ikind) = ghost_a .OR. floating_a
         atomnumber(ikind) = za
         c6d2(ikind) = disp_a%c6
         radd2(ikind) = disp_a%vdw_radii
      END DO

      ALLOCATE (rcpbc(3, natom))
      DO iatom = 1, natom
         rcpbc(:, iatom) = pbc(particle_set(iatom)%r(:), cell)
      END DO

      rcut = 2._dp*dispersion_env%rc_disp
      rc2 = rcut*rcut

      maxc = SIZE(dispersion_env%c6ab, 3)
      max_elem = SIZE(dispersion_env%c6ab, 1)
      alp6 = dispersion_env%alp
      alp8 = alp6 + 2._dp
      s6 = dispersion_env%s6
      s8 = dispersion_env%s8
      s9 = dispersion_env%s6
      rs6 = dispersion_env%sr6
      rs8 = 1._dp
      a1 = dispersion_env%a1
      a2 = dispersion_env%a2
      eps_cn = dispersion_env%eps_cn
      e6tot = 0._dp
      e8tot = 0._dp
      e9tot = 0._dp
      esrb = 0._dp
      domol = dispersion_env%domol
      ! molecule correction
      kgc8 = dispersion_env%kgc8
      IF (domol) THEN
         CALL get_qs_env(qs_env=qs_env, molecule_set=molecule_set)
         ALLOCATE (atom2mol(natom))
         DO imol = 1, SIZE(molecule_set)
            DO iat = molecule_set(imol)%first_atom, molecule_set(imol)%last_atom
               atom2mol(iat) = imol
            END DO
         END DO
      END IF
      ! damping type
      idmp = 0
      IF (dispersion_env%pp_type == vdw_pairpot_dftd3) THEN
         idmp = 1
      ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
         idmp = 2
      END IF
      ! SRB parameters
      ssrb = dispersion_env%srb_params(1)
      gsrb = dispersion_env%srb_params(2)
      t1srb = dispersion_env%srb_params(3)
      t2srb = dispersion_env%srb_params(4)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) " Scaling parameter (s6) ", s6
         WRITE (unit_nr, *) " Scaling parameter (s8) ", s8
         IF (dispersion_env%pp_type == vdw_pairpot_dftd3) THEN
            WRITE (unit_nr, *) " Zero Damping parameter (sr6)", rs6
            WRITE (unit_nr, *) " Zero Damping parameter (sr8)", rs8
         ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
            WRITE (unit_nr, *) " BJ Damping parameter (a1) ", a1
            WRITE (unit_nr, *) " BJ Damping parameter (a2) ", a2
         END IF
         WRITE (unit_nr, *) " Cutoff coordination numbers", eps_cn
         IF (dispersion_env%lrc) THEN
            WRITE (unit_nr, *) " Apply a long range correction"
         END IF
         IF (dispersion_env%srb) THEN
            WRITE (unit_nr, *) " Apply a short range bond correction"
            WRITE (unit_nr, *) " SRB parameters (s,g,t1,t2) ", ssrb, gsrb, t1srb, t2srb
         END IF
         IF (domol) THEN
            WRITE (unit_nr, *) " Inter-molecule scaling parameter (s8) ", kgc8
         END IF
      END IF
      ! Calculate coordination numbers
      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)
      ALLOCATE (cnumbers(natom))
      cnumbers = 0._dp

      IF (calculate_forces .OR. debugall) THEN
         ALLOCATE (dcnum(natom))
         dcnum(:)%neighbors = 0
         DO iatom = 1, natom
            ALLOCATE (dcnum(iatom)%nlist(10), dcnum(iatom)%dvals(10), dcnum(iatom)%rik(3, 10))
         END DO
      ELSE
         ALLOCATE (dcnum(1))
      END IF

      CALL d3_cnumber(qs_env, dispersion_env, cnumbers, dcnum, exclude, atomnumber, &
                      calculate_forces, 1)

      CALL para_env%sum(cnumbers)
      ! for parallel runs we have to update dcnum on all processors
      IF (calculate_forces .OR. debugall) THEN
         CALL dcnum_distribute(dcnum, para_env)
         IF (unit_nr > 0 .AND. SIZE(dcnum) > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, *) "  ATOM       Coordination   Neighbors"
            DO i = 1, natom
               WRITE (unit_nr, "(I6,F20.10,I12)") i, cnumbers(i), dcnum(i)%neighbors
            END DO
            WRITE (unit_nr, *)
         END IF
      END IF

      nab = 0._dp
      nabc = 0._dp
      IF (dispersion_env%doabc) THEN
         rcc = 2._dp*dispersion_env%rc_disp
         CALL get_cell(cell=cell, periodic=periodic)
         sab_max(1) = rcc/plane_distance(1, 0, 0, cell)
         sab_max(2) = rcc/plane_distance(0, 1, 0, cell)
         sab_max(3) = rcc/plane_distance(0, 0, 1, cell)
         ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) " Calculate C9 Terms"
            WRITE (unit_nr, "(A,T20,I3,A,I3)") "  Search in cells ", -ncell(1), ":", ncell(1)
            WRITE (unit_nr, "(T20,I3,A,I3)") - ncell(2), ":", ncell(2)
            WRITE (unit_nr, "(T20,I3,A,I3)") - ncell(3), ":", ncell(3)
            WRITE (unit_nr, *)
         END IF
         IF (dispersion_env%c9cnst) THEN
            IF (unit_nr > 0) WRITE (unit_nr, *) " Use reference coordination numbers for C9 term"
            ALLOCATE (cnumfix(natom))
            cnumfix = 0._dp
            ! first use the default values
            DO iatom = 1, natom
               ikind = kind_of(iatom)
               CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
               cnumfix(iatom) = dispersion_env%cn(za)
            END DO
            ! now check for changes from default
            IF (ASSOCIATED(dispersion_env%cnkind)) THEN
               DO i = 1, SIZE(dispersion_env%cnkind)
                  ikind = dispersion_env%cnkind(i)%kind
                  cnum = dispersion_env%cnkind(i)%cnum
                  CPASSERT(ikind <= nkind)
                  CPASSERT(ikind > 0)
                  CALL get_atomic_kind(atomic_kind_set(ikind), natom=na, atom_list=atom_list)
                  DO ia = 1, na
                     iatom = atom_list(ia)
                     cnumfix(iatom) = cnum
                  END DO
               END DO
            END IF
            IF (ASSOCIATED(dispersion_env%cnlist)) THEN
               DO i = 1, SIZE(dispersion_env%cnlist)
                  DO ilist = 1, dispersion_env%cnlist(i)%natom
                     iatom = dispersion_env%cnlist(i)%atom(ilist)
                     cnumfix(iatom) = dispersion_env%cnlist(i)%cnum
                  END DO
               END DO
            END IF
            IF (unit_nr > 0) THEN
               DO i = 1, natom
                  IF (ABS(cnumbers(i) - cnumfix(i)) > 0.5_dp) THEN
                     WRITE (unit_nr, "(A,T20,A,I6,T41,2F10.3)") "  Difference in CN ", "Atom:", &
                        i, cnumbers(i), cnumfix(i)
                  END IF
               END DO
               WRITE (unit_nr, *)
            END IF
         END IF
      END IF

      sab_vdw => dispersion_env%sab_vdw

      num_pe = 1
      CALL neighbor_list_iterator_create(nl_iterator, sab_vdw, nthread=num_pe)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom, r=rij)

         IF (exclude(ikind) .OR. exclude(jkind)) CYCLE

         IF (.NOT. (dodisp(ikind) .AND. dodisp(jkind))) CYCLE

         za = atomnumber(ikind)
         zb = atomnumber(jkind)
         ! vdW potential
         dr = SQRT(SUM(rij(:)**2))
         IF (dr <= rcut) THEN
            nab = nab + 1._dp
            fac = 1._dp
            IF (iatom == jatom) fac = 0.5_dp
            IF (disp_a%type == dftd3_pp .AND. dr > 0.001_dp) THEN
               IF (dispersion_env%nd3_exclude_pair > 0) THEN
                  CALL exclude_d3_kind_pair(dispersion_env%d3_exclude_pair, ikind, jkind, &
                                            exclude=exclude_pair)
                  IF (exclude_pair) CYCLE
               END IF
               ! C6 coefficient
               CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, za, zb, &
                          cnumbers(iatom), cnumbers(jatom), dispersion_env%k3, c6, dc6a, dc6b)
               c8 = 3.0d0*c6*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb)
               dc8a = 3.0d0*dc6a*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb)
               dc8b = 3.0d0*dc6b*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb)
               r6 = dr**6
               r8 = r6*dr*dr
               s8i = s8
               IF (domol) THEN
                  IF (atom2mol(iatom) /= atom2mol(jatom)) THEN
                     s8i = kgc8
                  END IF
               END IF
               ! damping
               IF (idmp == 1) THEN
                  ! zero
                  CALL damping_d3(dr, dispersion_env%r0ab(za, zb), rs6, alp6, rcut, fdab6, dfdab6)
                  CALL damping_d3(dr, dispersion_env%r0ab(za, zb), rs8, alp8, rcut, fdab8, dfdab8)
                  e6 = s6*fac*c6*fdab6/r6
                  e8 = s8i*fac*c8*fdab8/r8
               ELSE IF (idmp == 2) THEN
                  ! BJ
                  r0ab = SQRT(3.0d0*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb))
                  f0ab = a1*r0ab + a2
                  fdab6 = 1.0_dp/(r6 + f0ab**6)
                  fdab8 = 1.0_dp/(r8 + f0ab**8)
                  e6 = s6*fac*c6*fdab6
                  e8 = s8i*fac*c8*fdab8
               ELSE
                  CPABORT("Unknown DFT-D3 damping function:")
               END IF
               IF (dispersion_env%srb .AND. dr .LT. 30.0d0) THEN
                  srbe = ssrb*(REAL((za*zb), KIND=dp))**t1srb*EXP(-gsrb*dr*dispersion_env%r0ab(za, zb)**t2srb)
                  esrb = esrb + srbe
                  evdw = evdw - srbe
               ELSE
                  srbe = 0.0_dp
               END IF
               evdw = evdw - e6 - e8
               e6tot = e6tot - e6
               e8tot = e8tot - e8
               IF (atenergy) THEN
                  atener(iatom) = atener(iatom) - 0.5_dp*(e6 + e8 + srbe)
                  atener(jatom) = atener(jatom) - 0.5_dp*(e6 + e8 + srbe)
               END IF
               IF (atex) THEN
                  atevdw(iatom) = atevdw(iatom) - 0.5_dp*(e6 + e8 + srbe)
                  atevdw(jatom) = atevdw(jatom) - 0.5_dp*(e6 + e8 + srbe)
               END IF
               IF (calculate_forces) THEN
                  ! damping
                  IF (idmp == 1) THEN
                     ! zero
                     de6 = -s6*c6/r6*(dfdab6 - 6._dp*fdab6/dr)
                     de8 = -s8i*c8/r8*(dfdab8 - 8._dp*fdab8/dr)
                  ELSE IF (idmp == 2) THEN
                     ! BJ
                     de6 = s6*c6*fdab6*fdab6*6.0_dp*dr**5
                     de8 = s8i*c8*fdab8*fdab8*8.0_dp*dr**7
                  ELSE
                     CPABORT("Unknown DFT-D3 damping function:")
                  END IF
                  fdij(:) = (de6 + de8)*rij(:)/dr*fac
                  IF (dispersion_env%srb .AND. dr .LT. 30.0d0) THEN
                     fdij(:) = fdij(:) + srbe*gsrb*dispersion_env%r0ab(za, zb)**t2srb*rij(:)/dr
                  END IF
                  atom_a = atom_of_kind(iatom)
                  atom_b = atom_of_kind(jatom)
                  force(ikind)%dispersion(:, atom_a) = force(ikind)%dispersion(:, atom_a) - fdij(:)
                  force(jkind)%dispersion(:, atom_b) = force(jkind)%dispersion(:, atom_b) + fdij(:)
                  IF (use_virial) THEN
                     CALL virial_pair_force(pv_virial_thread, -1._dp, fdij, rij)
                  END IF
                  ! forces from the r-dependence of the coordination numbers
                  IF (idmp == 1) THEN
                     ! zero
                     dc6a = -s6*fac*dc6a*fdab6/r6
                     dc6b = -s6*fac*dc6b*fdab6/r6
                     dc8a = -s8i*fac*dc8a*fdab8/r8
                     dc8b = -s8i*fac*dc8b*fdab8/r8
                  ELSE IF (idmp == 2) THEN
                     ! BJ
                     dc6a = -s6*fac*dc6a*fdab6
                     dc6b = -s6*fac*dc6b*fdab6
                     dc8a = -s8i*fac*dc8a*fdab8
                     dc8b = -s8i*fac*dc8b*fdab8
                  ELSE
                     CPABORT("Unknown DFT-D3 damping function:")
                  END IF
                  DO i = 1, dcnum(iatom)%neighbors
                     katom = dcnum(iatom)%nlist(i)
                     kkind = kind_of(katom)
                     rik = dcnum(iatom)%rik(:, i)
                     drk = SQRT(SUM(rik(:)**2))
                     fdik(:) = (dc6a + dc8a)*dcnum(iatom)%dvals(i)*rik(:)/drk
                     atom_c = atom_of_kind(katom)
                     force(ikind)%dispersion(:, atom_a) = force(ikind)%dispersion(:, atom_a) - fdik(:)
                     force(kkind)%dispersion(:, atom_c) = force(kkind)%dispersion(:, atom_c) + fdik(:)
                     IF (use_virial) THEN
                        CALL virial_pair_force(pv_virial_thread, -1._dp, fdik, rik)
                     END IF
                  END DO
                  DO i = 1, dcnum(jatom)%neighbors
                     katom = dcnum(jatom)%nlist(i)
                     kkind = kind_of(katom)
                     rik = dcnum(jatom)%rik(:, i)
                     drk = SQRT(SUM(rik(:)**2))
                     fdik(:) = (dc6b + dc8b)*dcnum(jatom)%dvals(i)*rik(:)/drk
                     atom_c = atom_of_kind(katom)
                     force(jkind)%dispersion(:, atom_b) = force(jkind)%dispersion(:, atom_b) - fdik(:)
                     force(kkind)%dispersion(:, atom_c) = force(kkind)%dispersion(:, atom_c) + fdik(:)
                     IF (use_virial) THEN
                        CALL virial_pair_force(pv_virial_thread, -1._dp, fdik, rik)
                     END IF
                  END DO
               END IF
               IF (dispersion_env%doabc) THEN
                  CALL get_iterator_info(nl_iterator, cell=cell_b)
                  hashb = cellhash(cell_b, ncell)
                  is000 = (ALL(cell_b == 0))
                  rb0(:) = MATMUL(cell%hmat, cell_b)
                  ra(:) = pbc(particle_set(iatom)%r(:), cell)
                  rb(:) = pbc(particle_set(jatom)%r(:), cell) + rb0
                  r2ab = SUM((ra - rb)**2)
                  DO icx = -ncell(1), ncell(1)
                     DO icy = -ncell(2), ncell(2)
                        DO icz = -ncell(3), ncell(3)
                           cell_c(1) = icx
                           cell_c(2) = icy
                           cell_c(3) = icz
                           hashc = cellhash(cell_c, ncell)
                           IF (is000 .AND. (ALL(cell_c == 0))) THEN
                              ! CASE 1: all atoms in (000), use only ordered triples
                              kstart = MAX(jatom + 1, iatom + 1)
                              fac0 = 1.0_dp
                           ELSE IF (is000) THEN
                              ! CASE 2: AB in (000), C in other cell
                              !         This case covers also all instances with BC in same
                              !         cell not (000)
                              kstart = 1
                              fac0 = 1.0_dp
                           ELSE
                              ! These are case 2 again, cycle
                              IF (hashc == hashb) CYCLE
                              IF (ALL(cell_c == 0)) CYCLE
                              ! CASE 3: A in (000) and B and C in different cells
                              kstart = 1
                              fac0 = 1.0_dp/3.0_dp
                           END IF
                           rc0 = MATMUL(cell%hmat, cell_c)
                           DO katom = kstart, natom
                              kkind = kind_of(katom)
                              IF (exclude(kkind) .OR. .NOT. dodisp(kkind)) CYCLE
                              rc(:) = rcpbc(:, katom) + rc0(:)
                              r2bc = SUM((rb - rc)**2)
                              IF (r2bc >= rc2) CYCLE
                              r2ca = SUM((rc - ra)**2)
                              IF (r2ca >= rc2) CYCLE
                              r2abc = r2ab*r2bc*r2ca
                              IF (r2abc <= 0.001_dp) CYCLE
                              IF (dispersion_env%nd3_exclude_pair > 0) THEN
                                 CALL exclude_d3_kind_pair(dispersion_env%d3_exclude_pair, ikind, jkind, &
                                                           kkind, exclude_pair)
                                 IF (exclude_pair) CYCLE
                              END IF
                              ! this is an empirical scaling
                              IF (r2abc <= 0.01*rc2*rc2*rc2) THEN
                                 kkind = kind_of(katom)
                                 atom_c = atom_of_kind(katom)
                                 zc = atomnumber(kkind)
                                 ! avoid double counting!
                                 fac = 1._dp
                                 IF (iatom == jatom .OR. iatom == katom .OR. jatom == katom) fac = 0.5_dp
                                 IF (iatom == jatom .AND. iatom == katom) fac = 1._dp/3._dp
                                 fac = fac*fac0
                                 nabc = nabc + 1._dp
                                 IF (dispersion_env%c9cnst) THEN
                                    CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, za, zb, &
                                               cnumfix(iatom), cnumfix(jatom), dispersion_env%k3, cc6ab, dcc6aba, dcc6abb)
                                    CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, zb, zc, &
                                               cnumfix(jatom), cnumfix(katom), dispersion_env%k3, cc6bc, dcc6bcb, dcc6bcc)
                                    CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, zc, za, &
                                               cnumfix(katom), cnumfix(iatom), dispersion_env%k3, cc6ca, dcc6cac, dcc6caa)
                                 ELSE
                                    CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, za, zb, &
                                               cnumbers(iatom), cnumbers(jatom), dispersion_env%k3, cc6ab, dcc6aba, dcc6abb)
                                    CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, zb, zc, &
                                               cnumbers(jatom), cnumbers(katom), dispersion_env%k3, cc6bc, dcc6bcb, dcc6bcc)
                                    CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, zc, za, &
                                               cnumbers(katom), cnumbers(iatom), dispersion_env%k3, cc6ca, dcc6cac, dcc6caa)
                                 END IF
                                 c9 = -SQRT(cc6ab*cc6bc*cc6ca)
                                 rabc = r2abc**(1._dp/6._dp)
                                 r0 = (dispersion_env%r0ab(za, zb)*dispersion_env%r0ab(zb, zc)* &
                                       dispersion_env%r0ab(zc, za))**(1._dp/3._dp)
                                 ! bug fixed 3.10.2017
                                 ! correct value from alp6=14 to 16 as used in original paper
                                 CALL damping_d3(rabc, r0, 4._dp/3._dp, 16.0_dp, rcut, fdabc, dfdabc)
                                 s1 = r2ab + r2bc - r2ca
                                 s2 = r2ab + r2ca - r2bc
                                 s3 = r2ca + r2bc - r2ab
                                 ang = 0.375_dp*s1*s2*s3/r2abc + 1.0_dp

                                 e9 = s9*fac*fdabc*c9*ang/r2abc**1.50d0
                                 evdw = evdw - e9
                                 e9tot = e9tot - e9
                                 IF (atenergy) THEN
                                    atener(iatom) = atener(iatom) - e9/3._dp
                                    atener(jatom) = atener(jatom) - e9/3._dp
                                    atener(katom) = atener(katom) - e9/3._dp
                                 END IF
                                 IF (atex) THEN
                                    atevdw(iatom) = atevdw(iatom) - e9/3._dp
                                    atevdw(jatom) = atevdw(jatom) - e9/3._dp
                                    atevdw(katom) = atevdw(katom) - e9/3._dp
                                 END IF

                                 IF (calculate_forces) THEN
                                    rab = rb - ra; rbc = rc - rb; rca = ra - rc
                                    de91 = s9*fac*dfdabc*c9*ang/r2abc**1.50d0
                                    de91 = -de91/3._dp*rabc + 3._dp*e9
                                    dea = s9*fac*fdabc*c9/r2abc**2.50d0*0.75_dp
                                    fdij(:) = de91*rab(:)/r2ab
                                    fdij(:) = fdij(:) + dea*s1*s2*s3*rab(:)/r2ab
                                    fdij(:) = fdij(:) - dea*(s2*s3 + s1*s3 - s1*s2)*rab(:)
                                    force(ikind)%dispersion(:, atom_a) = force(ikind)%dispersion(:, atom_a) - fdij(:)
                                    force(jkind)%dispersion(:, atom_b) = force(jkind)%dispersion(:, atom_b) + fdij(:)
                                    IF (use_virial) THEN
                                       CALL virial_pair_force(pv_virial_thread, -1._dp, fdij, rab)
                                    END IF
                                    fdij(:) = de91*rbc(:)/r2bc
                                    fdij(:) = fdij(:) + dea*s1*s2*s3*rbc(:)/r2bc
                                    fdij(:) = fdij(:) - dea*(s2*s3 - s1*s3 + s1*s2)*rbc(:)
                                    force(jkind)%dispersion(:, atom_b) = force(jkind)%dispersion(:, atom_b) - fdij(:)
                                    force(kkind)%dispersion(:, atom_c) = force(kkind)%dispersion(:, atom_c) + fdij(:)
                                    IF (use_virial) THEN
                                       CALL virial_pair_force(pv_virial_thread, -1._dp, fdij, rbc)
                                    END IF
                                    fdij(:) = de91*rca(:)/r2ca
                                    fdij(:) = fdij(:) + dea*s1*s2*s3*rca(:)/r2ca
                                    fdij(:) = fdij(:) - dea*(-s2*s3 + s1*s3 + s1*s2)*rca(:)
                                    force(kkind)%dispersion(:, atom_c) = force(kkind)%dispersion(:, atom_c) - fdij(:)
                                    force(ikind)%dispersion(:, atom_a) = force(ikind)%dispersion(:, atom_a) + fdij(:)
                                    IF (use_virial) THEN
                                       CALL virial_pair_force(pv_virial_thread, -1._dp, fdij, rca)
                                    END IF

                                    IF (.NOT. dispersion_env%c9cnst) THEN
                                       ! forces from the r-dependence of the coordination numbers
                                       ! atomic stress not implemented

                                       de91 = 0.5_dp*e9/cc6ab
                                       de921 = de91*dcc6aba
                                       de922 = de91*dcc6abb
                                       DO i = 1, dcnum(iatom)%neighbors
                                          latom = dcnum(iatom)%nlist(i)
                                          lkind = kind_of(latom)
                                          rik(1) = dcnum(iatom)%rik(1, i)
                                          rik(2) = dcnum(iatom)%rik(2, i)
                                          rik(3) = dcnum(iatom)%rik(3, i)
                                          drk = SQRT(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                                          fdik(:) = -de921*dcnum(iatom)%dvals(i)*rik(:)/drk
                                          atom_d = atom_of_kind(latom)
                                          force(ikind)%dispersion(:, atom_a) = force(ikind)%dispersion(:, atom_a) - fdik(:)
                                          force(lkind)%dispersion(:, atom_d) = force(lkind)%dispersion(:, atom_d) + fdik(:)
                                          IF (use_virial) THEN
                                             CALL virial_pair_force(pv_virial_thread, -1._dp, fdik, rik)
                                          END IF
                                       END DO
                                       DO i = 1, dcnum(jatom)%neighbors
                                          latom = dcnum(jatom)%nlist(i)
                                          lkind = kind_of(latom)
                                          rik(1) = dcnum(jatom)%rik(1, i)
                                          rik(2) = dcnum(jatom)%rik(2, i)
                                          rik(3) = dcnum(jatom)%rik(3, i)
                                          drk = SQRT(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                                          fdik(:) = -de922*dcnum(jatom)%dvals(i)*rik(:)/drk
                                          atom_d = atom_of_kind(latom)
                                          force(jkind)%dispersion(:, atom_b) = force(jkind)%dispersion(:, atom_b) - fdik(:)
                                          force(lkind)%dispersion(:, atom_d) = force(lkind)%dispersion(:, atom_d) + fdik(:)
                                          IF (use_virial) THEN
                                             CALL virial_pair_force(pv_virial_thread, -1._dp, fdik, rik)
                                          END IF
                                       END DO

                                       de91 = 0.5_dp*e9/cc6bc
                                       de921 = de91*dcc6bcb
                                       de922 = de91*dcc6bcc
                                       DO i = 1, dcnum(jatom)%neighbors
                                          latom = dcnum(jatom)%nlist(i)
                                          lkind = kind_of(latom)
                                          rik(1) = dcnum(jatom)%rik(1, i)
                                          rik(2) = dcnum(jatom)%rik(2, i)
                                          rik(3) = dcnum(jatom)%rik(3, i)
                                          drk = SQRT(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                                          fdik(:) = -de921*dcnum(jatom)%dvals(i)*rik(:)/drk
                                          atom_d = atom_of_kind(latom)
                                          force(jkind)%dispersion(:, atom_b) = force(jkind)%dispersion(:, atom_b) - fdik(:)
                                          force(lkind)%dispersion(:, atom_d) = force(lkind)%dispersion(:, atom_d) + fdik(:)
                                          IF (use_virial) THEN
                                             CALL virial_pair_force(pv_virial_thread, -1._dp, fdik, rik)
                                          END IF
                                       END DO
                                       DO i = 1, dcnum(katom)%neighbors
                                          latom = dcnum(katom)%nlist(i)
                                          lkind = kind_of(latom)
                                          rik(1) = dcnum(katom)%rik(1, i)
                                          rik(2) = dcnum(katom)%rik(2, i)
                                          rik(3) = dcnum(katom)%rik(3, i)
                                          drk = SQRT(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                                          fdik(:) = -de922*dcnum(katom)%dvals(i)*rik(:)/drk
                                          atom_d = atom_of_kind(latom)
                                          force(kkind)%dispersion(:, atom_c) = force(kkind)%dispersion(:, atom_c) - fdik(:)
                                          force(lkind)%dispersion(:, atom_d) = force(lkind)%dispersion(:, atom_d) + fdik(:)
                                          IF (use_virial) THEN
                                             CALL virial_pair_force(pv_virial_thread, -1._dp, fdik, rik)
                                          END IF
                                       END DO

                                       de91 = 0.5_dp*e9/cc6ca
                                       de921 = de91*dcc6cac
                                       de922 = de91*dcc6caa
                                       DO i = 1, dcnum(katom)%neighbors
                                          latom = dcnum(katom)%nlist(i)
                                          lkind = kind_of(latom)
                                          rik(1) = dcnum(katom)%rik(1, i)
                                          rik(2) = dcnum(katom)%rik(2, i)
                                          rik(3) = dcnum(katom)%rik(3, i)
                                          drk = SQRT(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                                          fdik(:) = -de921*dcnum(katom)%dvals(i)*rik(:)/drk
                                          atom_d = atom_of_kind(latom)
                                          force(kkind)%dispersion(:, atom_c) = force(kkind)%dispersion(:, atom_c) - fdik(:)
                                          force(lkind)%dispersion(:, atom_d) = force(lkind)%dispersion(:, atom_d) + fdik(:)
                                          IF (use_virial) THEN
                                             CALL virial_pair_force(pv_virial_thread, -1._dp, fdik, rik)
                                          END IF
                                       END DO
                                       DO i = 1, dcnum(iatom)%neighbors
                                          latom = dcnum(iatom)%nlist(i)
                                          lkind = kind_of(latom)
                                          rik(1) = dcnum(iatom)%rik(1, i)
                                          rik(2) = dcnum(iatom)%rik(2, i)
                                          rik(3) = dcnum(iatom)%rik(3, i)
                                          drk = SQRT(rik(1)*rik(1) + rik(2)*rik(2) + rik(3)*rik(3))
                                          fdik(:) = -de922*dcnum(iatom)%dvals(i)*rik(:)/drk
                                          atom_d = atom_of_kind(latom)
                                          force(ikind)%dispersion(:, atom_a) = force(ikind)%dispersion(:, atom_a) - fdik(:)
                                          force(lkind)%dispersion(:, atom_d) = force(lkind)%dispersion(:, atom_d) + fdik(:)
                                          IF (use_virial) THEN
                                             CALL virial_pair_force(pv_virial_thread, -1._dp, fdik, rik)
                                          END IF
                                       END DO
                                    END IF

                                 END IF

                              END IF
                           END DO
                        END DO
                     END DO
                  END DO

               END IF
            END IF
         END IF
      END DO

      virial%pv_virial = virial%pv_virial + pv_virial_thread

      CALL neighbor_list_iterator_release(nl_iterator)

      elrc6 = 0._dp
      elrc8 = 0._dp
      elrc9 = 0._dp
      ! Long range correction (atomic contributions not implemented)
      IF (dispersion_env%lrc) THEN
         ALLOCATE (cnkind(nkind))
         cnkind = 0._dp
         ! first use the default values
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
            cnkind(ikind) = dispersion_env%cn(za)
         END DO
         ! now check for changes from default
         IF (ASSOCIATED(dispersion_env%cnkind)) THEN
            DO i = 1, SIZE(dispersion_env%cnkind)
               ikind = dispersion_env%cnkind(i)%kind
               cnkind(ikind) = dispersion_env%cnkind(i)%cnum
            END DO
         END IF
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=na, z=za)
            CALL get_qs_kind(qs_kind_set(ikind), dispersion=disp_a, ghost=ghost_a, floating=floating_a)
            IF (.NOT. disp_a%defined .OR. ghost_a .OR. floating_a) CYCLE
            DO jkind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(jkind), natom=nb, z=zb)
               CALL get_qs_kind(qs_kind_set(jkind), dispersion=disp_b, ghost=ghost_b, floating=floating_b)
               IF (.NOT. disp_b%defined .OR. ghost_b .OR. floating_b) CYCLE
               IF (dispersion_env%nd3_exclude_pair > 0) THEN
                  CALL exclude_d3_kind_pair(dispersion_env%d3_exclude_pair, ikind, jkind, &
                                            exclude=exclude_pair)
                  IF (exclude_pair) CYCLE
               END IF
               CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, za, zb, &
                          cnkind(ikind), cnkind(jkind), dispersion_env%k3, cc6ab, dcc6aba, dcc6abb)
               elrc6 = elrc6 - s6*twopi*REAL(na*nb, KIND=dp)*cc6ab/(3._dp*rcut**3*cell%deth)
               c8 = 3.0d0*cc6ab*dispersion_env%r2r4(za)*dispersion_env%r2r4(zb)
               elrc8 = elrc8 - s8*twopi*REAL(na*nb, KIND=dp)*c8/(5._dp*rcut**5*cell%deth)
               IF (dispersion_env%doabc) THEN
                  DO kkind = 1, nkind
                     CALL get_atomic_kind(atomic_kind_set(kkind), natom=nc, z=zc)
                     CALL get_qs_kind(qs_kind_set(kkind), dispersion=disp_c, ghost=ghost_c, floating=floating_c)
                     IF (.NOT. disp_c%defined .OR. ghost_c .OR. floating_c) CYCLE
                     IF (dispersion_env%nd3_exclude_pair > 0) THEN
                        CALL exclude_d3_kind_pair(dispersion_env%d3_exclude_pair, ikind, jkind, kkind, &
                                                  exclude_pair)
                        IF (exclude_pair) CYCLE
                     END IF
                     CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, za, zb, &
                                cnkind(ikind), cnkind(jkind), dispersion_env%k3, cc6ab, dcc6aba, dcc6abb)
                     CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, zc, za, &
                                cnkind(kkind), cnkind(ikind), dispersion_env%k3, cc6ca, dcc6aba, dcc6abb)
                     CALL getc6(maxc, max_elem, dispersion_env%c6ab, dispersion_env%maxci, zb, zc, &
                                cnkind(jkind), cnkind(kkind), dispersion_env%k3, cc6bc, dcc6aba, dcc6abb)
                     c9 = -SQRT(cc6ab*cc6bc*cc6ca)
                     elrc9 = elrc9 - s9*64._dp*twopi*REAL(na*nb*nc, KIND=dp)*c9/(6._dp*rcut**3*cell%deth**2)
                  END DO
               END IF
            END DO
         END DO
         IF (use_virial) THEN
            IF (para_env%is_source()) THEN
               DO i = 1, 3
                  virial%pv_virial(i, i) = virial%pv_virial(i, i) + (elrc6 + elrc8 + 2._dp*elrc9)
               END DO
            END IF
         END IF
         DEALLOCATE (cnkind)
      END IF

      DEALLOCATE (cnumbers)
      IF (dispersion_env%doabc .AND. dispersion_env%c9cnst) THEN
         DEALLOCATE (cnumfix)
      END IF
      IF (calculate_forces .OR. debugall) THEN
         DO iatom = 1, natom
            DEALLOCATE (dcnum(iatom)%nlist, dcnum(iatom)%dvals, dcnum(iatom)%rik)
         END DO
         DEALLOCATE (dcnum)
      ELSE
         DEALLOCATE (dcnum)
      END IF

      ! set dispersion energy
      IF (para_env%is_source()) THEN
         evdw = evdw + (elrc6 + elrc8 + elrc9)
      END IF

      CALL para_env%sum(e6tot)
      CALL para_env%sum(e8tot)
      CALL para_env%sum(e9tot)
      CALL para_env%sum(nab)
      CALL para_env%sum(nabc)
      e6tot = e6tot + elrc6
      e8tot = e8tot + elrc8
      e9tot = e9tot + elrc9
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(A,F20.0)") "  E6 vdW terms              :", nab
         WRITE (unit_nr, *) " E6 vdW energy [au/kcal]   :", e6tot, e6tot*kcalmol
         WRITE (unit_nr, *) " E8 vdW energy [au/kcal]   :", e8tot, e8tot*kcalmol
         WRITE (unit_nr, *) " %E8 on total vdW energy   :", e8tot/evdw*100._dp
         WRITE (unit_nr, "(A,F20.0)") "  E9 vdW terms              :", nabc
         WRITE (unit_nr, *) " E9 vdW energy [au/kcal]   :", e9tot, e9tot*kcalmol
         WRITE (unit_nr, *) " %E9 on total vdW energy   :", e9tot/evdw*100._dp
         IF (dispersion_env%lrc) THEN
            WRITE (unit_nr, *) " E LRC C6 [au/kcal]        :", elrc6, elrc6*kcalmol
            WRITE (unit_nr, *) " E LRC C8 [au/kcal]        :", elrc8, elrc8*kcalmol
            WRITE (unit_nr, *) " E LRC C9 [au/kcal]        :", elrc9, elrc9*kcalmol
         END IF
      END IF

      DEALLOCATE (dodisp, exclude, atomnumber, rcpbc, radd2, c6d2)

      IF (domol) THEN
         DEALLOCATE (atom2mol)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_dispersion_d3_pairpot

! **************************************************************************************************
!> \brief ...
!> \param c6ab ...
!> \param maxci ...
!> \param filename ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE dftd3_c6_param(c6ab, maxci, filename, para_env)

      REAL(KIND=dp), DIMENSION(:, :, :, :, :)            :: c6ab
      INTEGER, DIMENSION(:)                              :: maxci
      CHARACTER(LEN=*)                                   :: filename
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: funit, iadr, iat, jadr, jat, kk, nl, &
                                                            nlines, nn
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: pars

      IF (para_env%is_source()) THEN
         ! Read the DFT-D3 C6AB parameters from file "filename"
         CALL open_file(file_name=filename, unit_number=funit, file_form="FORMATTED")
         READ (funit, *) nl, nlines
      END IF
      CALL para_env%bcast(nl)
      CALL para_env%bcast(nlines)
      ALLOCATE (pars(nl))
      IF (para_env%is_source()) THEN
         READ (funit, *) pars(1:nl)
         CALL close_file(unit_number=funit)
      END IF
      CALL para_env%bcast(pars)

      ! Store C6AB coefficients in an array
      c6ab = -1._dp
      maxci = 0
      kk = 1
      DO nn = 1, nlines
         iat = NINT(pars(kk + 1))
         jat = NINT(pars(kk + 2))
         CALL limit(iat, jat, iadr, jadr)
         maxci(iat) = MAX(maxci(iat), iadr)
         maxci(jat) = MAX(maxci(jat), jadr)
         c6ab(iat, jat, iadr, jadr, 1) = pars(kk)
         c6ab(iat, jat, iadr, jadr, 2) = pars(kk + 3)
         c6ab(iat, jat, iadr, jadr, 3) = pars(kk + 4)

         c6ab(jat, iat, jadr, iadr, 1) = pars(kk)
         c6ab(jat, iat, jadr, iadr, 2) = pars(kk + 4)
         c6ab(jat, iat, jadr, iadr, 3) = pars(kk + 3)
         kk = (nn*5) + 1
      END DO

      DEALLOCATE (pars)

   END SUBROUTINE dftd3_c6_param

! **************************************************************************************************
!> \brief ...
!> \param iat ...
!> \param jat ...
!> \param iadr ...
!> \param jadr ...
! **************************************************************************************************
   SUBROUTINE limit(iat, jat, iadr, jadr)
      INTEGER                                            :: iat, jat, iadr, jadr

      INTEGER                                            :: i

      iadr = 1
      jadr = 1
      i = 100
      DO WHILE (iat .GT. 100)
         iat = iat - 100
         iadr = iadr + 1
      END DO

      i = 100
      DO WHILE (jat .GT. 100)
         jat = jat - 100
         jadr = jadr + 1
      END DO
   END SUBROUTINE limit

! **************************************************************************************************
!> \brief ...
!> \param rab ...
!> \param rcutab ...
!> \param srn ...
!> \param alpn ...
!> \param rcut ...
!> \param fdab ...
!> \param dfdab ...
! **************************************************************************************************
   SUBROUTINE damping_d3(rab, rcutab, srn, alpn, rcut, fdab, dfdab)

      REAL(KIND=dp), INTENT(IN)                          :: rab, rcutab, srn, alpn, rcut
      REAL(KIND=dp), INTENT(OUT)                         :: fdab, dfdab

      REAL(KIND=dp)                                      :: a, b, c, d, dd, dfab, dfcc, dz, fab, &
                                                            fcc, rl, rr, ru, z, zz

      rl = rcut - 1._dp
      ru = rcut
      IF (rab >= ru) THEN
         fcc = 0._dp
         dfcc = 0._dp
      ELSEIF (rab <= rl) THEN
         fcc = 1._dp
         dfcc = 0._dp
      ELSE
         z = rab*rab - rl*rl
         dz = 2._dp*rab
         zz = z*z*z
         d = ru*ru - rl*rl
         dd = d*d*d
         a = -10._dp/dd
         b = 15._dp/(dd*d)
         c = -6._dp/(dd*d*d)
         fcc = 1._dp + zz*(a + b*z + c*z*z)
         dfcc = zz*dz/z*(3._dp*a + 4._dp*b*z + 5._dp*c*z*z)
      END IF

      rr = 6._dp*(rab/(srn*rcutab))**(-alpn)
      fab = 1._dp/(1._dp + rr)
      dfab = fab*fab*rr*alpn/rab
      fdab = fab*fcc
      dfdab = dfab*fcc + fab*dfcc

   END SUBROUTINE damping_d3

! **************************************************************************************************
!> \brief ...
!> \param maxc ...
!> \param max_elem ...
!> \param c6ab ...
!> \param mxc ...
!> \param iat ...
!> \param jat ...
!> \param nci ...
!> \param ncj ...
!> \param k3 ...
!> \param c6 ...
!> \param dc6a ...
!> \param dc6b ...
! **************************************************************************************************
   SUBROUTINE getc6(maxc, max_elem, c6ab, mxc, iat, jat, nci, ncj, k3, c6, dc6a, dc6b)

      INTEGER, INTENT(IN)                                :: maxc, max_elem
      REAL(KIND=dp), INTENT(IN) :: c6ab(max_elem, max_elem, maxc, maxc, 3)
      INTEGER, INTENT(IN)                                :: mxc(max_elem), iat, jat
      REAL(KIND=dp), INTENT(IN)                          :: nci, ncj, k3
      REAL(KIND=dp), INTENT(OUT)                         :: c6, dc6a, dc6b

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: c6mem, cn1, cn2, csum, dtmpa, dtmpb, &
                                                            dwa, dwb, dza, dzb, r, rsave, rsum, &
                                                            tmp1

! the exponential is sensitive to numerics
! when nci or ncj is much larger than cn1/cn2

      c6mem = -1.0e+99_dp
      rsave = 1.0e+99_dp
      rsum = 0.0_dp
      csum = 0.0_dp
      dza = 0.0_dp
      dzb = 0.0_dp
      dwa = 0.0_dp
      dwb = 0.0_dp
      c6 = 0.0_dp
      DO i = 1, mxc(iat)
         DO j = 1, mxc(jat)
            c6 = c6ab(iat, jat, i, j, 1)
            IF (c6 > 0.0_dp) THEN
               cn1 = c6ab(iat, jat, i, j, 2)
               cn2 = c6ab(iat, jat, i, j, 3)
               ! distance
               r = (cn1 - nci)**2 + (cn2 - ncj)**2
               IF (r < rsave) THEN
                  rsave = r
                  c6mem = c6
               END IF
               tmp1 = EXP(k3*r)
               dtmpa = -2.0_dp*k3*(cn1 - nci)*tmp1
               dtmpb = -2.0_dp*k3*(cn2 - ncj)*tmp1
               rsum = rsum + tmp1
               csum = csum + tmp1*c6
               dza = dza + dtmpa*c6
               dwa = dwa + dtmpa
               dzb = dzb + dtmpb*c6
               dwb = dwb + dtmpb
            END IF
         END DO
      END DO

      IF (c6 == 0.0_dp) c6mem = 0.0_dp

      IF (rsum > 1.0e-66_dp) THEN
         c6 = csum/rsum
         dc6a = (dza - c6*dwa)/rsum
         dc6b = (dzb - c6*dwb)/rsum
      ELSE
         c6 = c6mem
         dc6a = 0._dp
         dc6b = 0._dp
      END IF

   END SUBROUTINE getc6

END MODULE qs_dispersion_d3
