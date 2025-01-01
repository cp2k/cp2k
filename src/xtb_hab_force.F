!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of xTB Hamiltonian derivative
!>        Reference: Stefan Grimme, Christoph Bannwarth, Philip Shushkov
!>                   JCTC 13, 1989-2009, (2017)
!>                   DOI: 10.1021/acs.jctc.7b00118
!> \author JGH
! **************************************************************************************************
MODULE xtb_hab_force
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_overlap,                      ONLY: overlap_ab
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_dispersion_cnum,              ONLY: cnumber_init,&
                                              cnumber_release,&
                                              dcnum_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                              get_memory_usage
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: create_sab_matrix
   USE xtb_hcore,                       ONLY: gfn1_huckel,&
                                              gfn1_kpair
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_hab_force'

   PUBLIC :: build_xtb_hab_force

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_matrix ...
! **************************************************************************************************
   SUBROUTINE build_xtb_hab_force(qs_env, p_matrix)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: p_matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_xtb_hab_force'

      INTEGER :: atom_a, atom_b, atom_c, handle, i, iatom, ic, icol, ikind, img, ir, irow, iset, &
         j, jatom, jkind, jset, katom, kkind, la, lb, ldsab, maxder, n1, n2, na, natom, natorb_a, &
         natorb_b, nb, ncoa, ncob, nderivatives, nimg, nkind, nsa, nsb, nseta, nsetb, sgfa, sgfb, &
         za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(25)                             :: laoa, laob, naoa, naob
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: defined, diagblock, found, use_virial
      REAL(KIND=dp)                                      :: dfp, dhij, dr, drk, drx, f0, fhua, fhub, &
                                                            fhud, foab, hij, rcova, rcovab, rcovb, &
                                                            rrab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cnumbers
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dfblock, dhuckel, huckel, owork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: kijab
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, fdika, fdikb, force_ab, rij, rik
      REAL(KIND=dp), DIMENSION(5)                        :: dpia, dpib, kpolya, kpolyb, pia, pib
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fblock, pblock, rpgfa, rpgfb, sblock, &
                                                            scon_a, scon_b, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(block_p_type), DIMENSION(2:4)                 :: dsblocks
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_s
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      NULLIFY (matrix_h, matrix_s, atomic_kind_set, qs_kind_set, sab_orb)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      sab_orb=sab_orb)

      CPASSERT(dft_control%qs_control%xtb_control%gfn_type == 1)

      nkind = SIZE(atomic_kind_set)
      xtb_control => dft_control%qs_control%xtb_control
      nimg = dft_control%nimages
      nderivatives = 1
      maxder = ncoset(nderivatives)

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      NULLIFY (force)
      CALL get_qs_env(qs_env=qs_env, force=force)
      use_virial = .FALSE.
      CPASSERT(nimg == 1)

      ! set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ! allocate overlap matrix
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
      CALL dbcsr_allocate_matrix_set(matrix_s, maxder, nimg)
      CALL create_sab_matrix(ks_env, matrix_s, "xTB OVERLAP MATRIX", basis_set_list, basis_set_list, &
                             sab_orb, .TRUE.)
      ! initialize H matrix
      CALL dbcsr_allocate_matrix_set(matrix_h, 1, nimg)
      DO img = 1, nimg
         ALLOCATE (matrix_h(1, img)%matrix)
         CALL dbcsr_create(matrix_h(1, img)%matrix, template=matrix_s(1, 1)%matrix, &
                           name="HAMILTONIAN MATRIX")
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_h(1, img)%matrix, sab_orb)
      END DO

      ! Calculate coordination numbers
      ! needed for effective atomic energy levels (Eq. 12)
      ! code taken from D3 dispersion energy
      CALL cnumber_init(qs_env, cnumbers, dcnum, 1, .TRUE.)

      ! Calculate Huckel parameters
      CALL gfn1_huckel(qs_env, cnumbers, huckel, dhuckel, .TRUE.)

      ! Calculate KAB parameters and electronegativity correction
      CALL gfn1_kpair(qs_env, kijab)

      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined, natorb=natorb_a)
         IF (.NOT. defined .OR. natorb_a < 1) CYCLE
         CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
         CALL get_xtb_atom_param(xtb_atom_b, defined=defined, natorb=natorb_b)
         IF (.NOT. defined .OR. natorb_b < 1) CYCLE

         dr = SQRT(SUM(rij(:)**2))

         ! atomic parameters
         CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa, lao=laoa, rcov=rcova, &
                                 nshell=nsa, kpoly=kpolya)
         CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob, lao=laob, rcov=rcovb, &
                                 nshell=nsb, kpoly=kpolyb)

         ic = 1
         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)
         NULLIFY (sblock, fblock)
         CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=matrix_h(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=fblock, found=found)
         CPASSERT(found)

         NULLIFY (pblock)
         CALL dbcsr_get_block_p(matrix=p_matrix, &
                                row=irow, col=icol, block=pblock, found=found)
         CPASSERT(ASSOCIATED(pblock))
         DO i = 2, 4
            NULLIFY (dsblocks(i)%block)
            CALL dbcsr_get_block_p(matrix=matrix_s(i, ic)%matrix, &
                                   row=irow, col=icol, BLOCK=dsblocks(i)%block, found=found)
            CPASSERT(found)
         END DO

         ! overlap
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         scon_a => basis_set_a%scon
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         scon_b => basis_set_b%scon
         zetb => basis_set_b%zet

         ldsab = get_memory_usage(qs_kind_set, "ORB", "ORB")
         ALLOCATE (oint(ldsab, ldsab, maxder), owork(ldsab, ldsab))
         ALLOCATE (sint(natorb_a, natorb_b, maxder))
         sint = 0.0_dp

         DO iset = 1, nseta
            ncoa = npgfa(iset)*ncoset(la_max(iset))
            n1 = npgfa(iset)*(ncoset(la_max(iset)) - ncoset(la_min(iset) - 1))
            sgfa = first_sgfa(1, iset)
            DO jset = 1, nsetb
               IF (set_radius_a(iset) + set_radius_b(jset) < dr) CYCLE
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               n2 = npgfb(jset)*(ncoset(lb_max(jset)) - ncoset(lb_min(jset) - 1))
               sgfb = first_sgfb(1, jset)
               CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                               lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                               rij, sab=oint(:, :, 1), dab=oint(:, :, 2:4))
               ! Contraction
               DO i = 1, 4
                  CALL contraction(oint(:, :, i), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                   cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
                  CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, i), sgfa, sgfb, trans=.FALSE.)
               END DO
            END DO
         END DO
         ! update S matrix
         IF (iatom <= jatom) THEN
            sblock(:, :) = sblock(:, :) + sint(:, :, 1)
         ELSE
            sblock(:, :) = sblock(:, :) + TRANSPOSE(sint(:, :, 1))
         END IF
         DO i = 2, 4
            IF (iatom <= jatom) THEN
               dsblocks(i)%block(:, :) = dsblocks(i)%block(:, :) + sint(:, :, i)
            ELSE
               dsblocks(i)%block(:, :) = dsblocks(i)%block(:, :) - TRANSPOSE(sint(:, :, i))
            END IF
         END DO

         ! Calculate Pi = Pia * Pib (Eq. 11)
         rcovab = rcova + rcovb
         rrab = SQRT(dr/rcovab)
         pia(1:nsa) = 1._dp + kpolya(1:nsa)*rrab
         pib(1:nsb) = 1._dp + kpolyb(1:nsb)*rrab
         IF (dr > 1.e-6_dp) THEN
            drx = 0.5_dp/rrab/rcovab
         ELSE
            drx = 0.0_dp
         END IF
         dpia(1:nsa) = drx*kpolya(1:nsa)
         dpib(1:nsb) = drx*kpolyb(1:nsb)

         ! diagonal block
         diagblock = .FALSE.
         IF (iatom == jatom .AND. dr < 0.001_dp) diagblock = .TRUE.
         !
         ! Eq. 10
         !
         IF (diagblock) THEN
            DO i = 1, natorb_a
               na = naoa(i)
               fblock(i, i) = fblock(i, i) + huckel(na, iatom)
            END DO
         ELSE
            DO j = 1, natorb_b
               nb = naob(j)
               DO i = 1, natorb_a
                  na = naoa(i)
                  hij = 0.5_dp*(huckel(na, iatom) + huckel(nb, jatom))*pia(na)*pib(nb)
                  IF (iatom <= jatom) THEN
                     fblock(i, j) = fblock(i, j) + hij*sint(i, j, 1)*kijab(i, j, ikind, jkind)
                  ELSE
                     fblock(j, i) = fblock(j, i) + hij*sint(i, j, 1)*kijab(i, j, ikind, jkind)
                  END IF
               END DO
            END DO
         END IF

         f0 = 1.0_dp
         IF (irow == iatom) f0 = -1.0_dp
         ! Derivative wrt coordination number
         fhua = 0.0_dp
         fhub = 0.0_dp
         fhud = 0.0_dp
         IF (diagblock) THEN
            DO i = 1, natorb_a
               la = laoa(i)
               na = naoa(i)
               fhud = fhud + pblock(i, i)*dhuckel(na, iatom)
            END DO
         ELSE
            DO j = 1, natorb_b
               lb = laob(j)
               nb = naob(j)
               DO i = 1, natorb_a
                  la = laoa(i)
                  na = naoa(i)
                  hij = 0.5_dp*pia(na)*pib(nb)
                  IF (iatom <= jatom) THEN
                     fhua = fhua + hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)*pblock(i, j)*dhuckel(na, iatom)
                     fhub = fhub + hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)*pblock(i, j)*dhuckel(nb, jatom)
                  ELSE
                     fhua = fhua + hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)*pblock(j, i)*dhuckel(na, iatom)
                     fhub = fhub + hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)*pblock(j, i)*dhuckel(nb, jatom)
                  END IF
               END DO
            END DO
            IF (iatom /= jatom) THEN
               fhua = 2.0_dp*fhua
               fhub = 2.0_dp*fhub
            END IF
         END IF
         ! iatom
         atom_a = atom_of_kind(iatom)
         DO i = 1, dcnum(iatom)%neighbors
            katom = dcnum(iatom)%nlist(i)
            kkind = kind_of(katom)
            atom_c = atom_of_kind(katom)
            rik = dcnum(iatom)%rik(:, i)
            drk = SQRT(SUM(rik(:)**2))
            IF (drk > 1.e-3_dp) THEN
               fdika(:) = fhua*dcnum(iatom)%dvals(i)*rik(:)/drk
               force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a) - fdika(:)
               force(kkind)%all_potential(:, atom_c) = force(kkind)%all_potential(:, atom_c) + fdika(:)
               fdikb(:) = fhud*dcnum(iatom)%dvals(i)*rik(:)/drk
               force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a) - fdikb(:)
               force(kkind)%all_potential(:, atom_c) = force(kkind)%all_potential(:, atom_c) + fdikb(:)
            END IF
         END DO
         ! jatom
         atom_b = atom_of_kind(jatom)
         DO i = 1, dcnum(jatom)%neighbors
            katom = dcnum(jatom)%nlist(i)
            kkind = kind_of(katom)
            atom_c = atom_of_kind(katom)
            rik = dcnum(jatom)%rik(:, i)
            drk = SQRT(SUM(rik(:)**2))
            IF (drk > 1.e-3_dp) THEN
               fdik(:) = fhub*dcnum(jatom)%dvals(i)*rik(:)/drk
               force(jkind)%all_potential(:, atom_b) = force(jkind)%all_potential(:, atom_b) - fdik(:)
               force(kkind)%all_potential(:, atom_c) = force(kkind)%all_potential(:, atom_c) + fdik(:)
            END IF
         END DO
         IF (diagblock) THEN
            force_ab = 0._dp
         ELSE
            ! force from R dendent Huckel element
            n1 = SIZE(fblock, 1)
            n2 = SIZE(fblock, 2)
            ALLOCATE (dfblock(n1, n2))
            dfblock = 0.0_dp
            DO j = 1, natorb_b
               lb = laob(j)
               nb = naob(j)
               DO i = 1, natorb_a
                  la = laoa(i)
                  na = naoa(i)
                  dhij = 0.5_dp*(huckel(na, iatom) + huckel(nb, jatom))*(dpia(na)*pib(nb) + pia(na)*dpib(nb))
                  IF (iatom <= jatom) THEN
                     dfblock(i, j) = dfblock(i, j) + dhij*sint(i, j, 1)*kijab(i, j, ikind, jkind)
                  ELSE
                     dfblock(j, i) = dfblock(j, i) + dhij*sint(i, j, 1)*kijab(i, j, ikind, jkind)
                  END IF
               END DO
            END DO
            dfp = f0*SUM(dfblock(:, :)*pblock(:, :))
            DO ir = 1, 3
               foab = 2.0_dp*dfp*rij(ir)/dr
               ! force from overlap matrix contribution to H
               DO j = 1, natorb_b
                  lb = laob(j)
                  nb = naob(j)
                  DO i = 1, natorb_a
                     la = laoa(i)
                     na = naoa(i)
                     hij = 0.5_dp*(huckel(na, iatom) + huckel(nb, jatom))*pia(na)*pib(nb)
                     IF (iatom <= jatom) THEN
                        foab = foab + 2.0_dp*hij*sint(i, j, ir + 1)*pblock(i, j)*kijab(i, j, ikind, jkind)
                     ELSE
                        foab = foab - 2.0_dp*hij*sint(i, j, ir + 1)*pblock(j, i)*kijab(i, j, ikind, jkind)
                     END IF
                  END DO
               END DO
               force_ab(ir) = foab
            END DO
            DEALLOCATE (dfblock)
         END IF

         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)
         IF (irow == iatom) force_ab = -force_ab
         force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a) - force_ab(:)
         force(jkind)%all_potential(:, atom_b) = force(jkind)%all_potential(:, atom_b) + force_ab(:)

         DEALLOCATE (oint, owork, sint)

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(matrix_h, 1)
         DO img = 1, nimg
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
      END DO
      CALL dbcsr_deallocate_matrix_set(matrix_s)
      CALL dbcsr_deallocate_matrix_set(matrix_h)

      ! deallocate coordination numbers
      CALL cnumber_release(cnumbers, dcnum, .TRUE.)

      ! deallocate Huckel parameters
      DEALLOCATE (huckel, dhuckel)
      ! deallocate KAB parameters
      DEALLOCATE (kijab)

      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_xtb_hab_force

END MODULE xtb_hab_force

