!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of charge response in xTB (EEQ only)
!>        Reference: Stefan Grimme, Christoph Bannwarth, Philip Shushkov
!>                   JCTC 13, 1989-2009, (2017)
!>                   DOI: 10.1021/acs.jctc.7b00118
!> \author JGH
! **************************************************************************************************
MODULE xtb_qresp
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_overlap,                      ONLY: overlap_ab
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE eeq_input,                       ONLY: eeq_solver_type
   USE input_constants,                 ONLY: vdw_pairpot_dftd4
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_dispersion_cnum,              ONLY: cnumber_init,&
                                              cnumber_release,&
                                              dcnum_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                              get_memory_usage
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE xtb_eeq,                         ONLY: xtb_eeq_calculation,&
                                              xtb_eeq_lagrange
   USE xtb_hcore,                       ONLY: gfn0_huckel,&
                                              gfn0_kpair
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_qresp'

   PUBLIC :: build_xtb_qresp

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param qresp ...
! **************************************************************************************************
   SUBROUTINE build_xtb_qresp(qs_env, qresp)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: qresp

      INTEGER                                            :: gfn_type
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      gfn_type = dft_control%qs_control%xtb_control%gfn_type

      qresp = 0.0_dp
      SELECT CASE (gfn_type)
      CASE (0)
         CALL build_gfn0_qresp(qs_env, qresp)
      CASE (1)
         CPABORT("QRESP: gfn_type = 1 not available")
      CASE (2)
         CPABORT("QRESP: gfn_type = 2 not available")
      CASE DEFAULT
         CPABORT("QRESP: Unknown gfn_type")
      END SELECT

   END SUBROUTINE build_xtb_qresp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param qresp ...
! **************************************************************************************************
   SUBROUTINE build_gfn0_qresp(qs_env, qresp)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: qresp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_gfn0_qresp'

      INTEGER :: atom_a, atom_b, handle, i, iatom, ic, icol, ikind, img, irow, iset, j, jatom, &
         jkind, jset, ldsab, lmaxa, lmaxb, maxder, n1, n2, na, natom, natorb_a, natorb_b, nb, &
         ncoa, ncob, nderivatives, nimg, nkind, nsa, nsb, nseta, nsetb, sgfa, sgfb, za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(25)                             :: laoa, laob, naoa, naob
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, diagblock, do_nonbonded, found
      REAL(KIND=dp)                                      :: dr, drx, eeq_energy, ef_energy, etaa, &
                                                            etab, f2, fqa, fqb, hij, qlambda, &
                                                            rcova, rcovab, rcovb, rrab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: charges, cnumbers, dcharges
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dhuckel, dqhuckel, huckel, owork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: kijab
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(5)                        :: hena, henb, kpolya, kpolyb, pia, pib
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock, rpgfa, rpgfb, scon_a, scon_b, &
                                                            zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_w
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(eeq_solver_type)                              :: eeq_sparam
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      NULLIFY (matrix_p, atomic_kind_set, qs_kind_set, sab_orb, ks_env)
      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb)

      nkind = SIZE(atomic_kind_set)
      xtb_control => dft_control%qs_control%xtb_control
      do_nonbonded = xtb_control%do_nonbonded
      nimg = dft_control%nimages
      nderivatives = 0
      maxder = ncoset(nderivatives)

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      NULLIFY (rho, matrix_w)
      CALL get_qs_env(qs_env=qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

      IF (SIZE(matrix_p, 1) == 2) THEN
         DO img = 1, nimg
            CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            CALL dbcsr_add(matrix_w(1, img)%matrix, matrix_w(2, img)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         END DO
      END IF

      NULLIFY (cell_to_index)
      IF (nimg > 1) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      ! set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ! Calculate coordination numbers
      ! needed for effective atomic energy levels
      ! code taken from D3 dispersion energy
      CALL cnumber_init(qs_env, cnumbers, dcnum, 2, .TRUE.)

      ALLOCATE (charges(natom), dcharges(natom))
      charges = 0.0_dp
      CALL xtb_eeq_calculation(qs_env, charges, cnumbers, eeq_sparam, eeq_energy, ef_energy, qlambda)
      dcharges = qlambda/REAL(para_env%num_pe, KIND=dp)

      CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env)
      IF (dispersion_env%pp_type == vdw_pairpot_dftd4 .AND. dispersion_env%ext_charges) THEN
         IF (ASSOCIATED(dispersion_env%dcharges)) THEN
            dcharges(1:natom) = dcharges(1:natom) + dispersion_env%dcharges(1:natom)
         ELSE
            CPWARN("gfn0-xTB variational dipole DFTD4 contribution missing")
         END IF
      END IF

      ! Calculate Huckel parameters
      CALL gfn0_huckel(qs_env, cnumbers, charges, huckel, dhuckel, dqhuckel, .TRUE.)

      ! Calculate KAB parameters and electronegativity correction
      CALL gfn0_kpair(qs_env, kijab)

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
         CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa, lao=laoa, rcov=rcova, eta=etaa, &
                                 lmax=lmaxa, nshell=nsa, kpoly=kpolya, hen=hena)
         CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob, lao=laob, rcov=rcovb, eta=etab, &
                                 lmax=lmaxb, nshell=nsb, kpoly=kpolyb, hen=henb)

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cell(1), cell(2), cell(3))
            CPASSERT(ic > 0)
         END IF

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         NULLIFY (pblock)
         CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                row=irow, col=icol, block=pblock, found=found)
         CPASSERT(ASSOCIATED(pblock))

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
                               rij, sab=oint(:, :, 1))
               ! Contraction
               CALL contraction(oint(:, :, 1), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
               CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, 1), sgfa, sgfb, trans=.FALSE.)
            END DO
         END DO

         ! Calculate Pi = Pia * Pib (Eq. 11)
         rcovab = rcova + rcovb
         rrab = SQRT(dr/rcovab)
         pia(1:nsa) = 1._dp + kpolya(1:nsa)*rrab
         pib(1:nsb) = 1._dp + kpolyb(1:nsb)*rrab

         ! diagonal block
         diagblock = .FALSE.
         IF (iatom == jatom .AND. dr < 0.001_dp) diagblock = .TRUE.
         !
         f2 = 1.0_dp
         IF (iatom /= jatom) f2 = 2.0_dp
         ! Derivative wrt coordination number
         fqa = 0.0_dp
         fqb = 0.0_dp
         IF (diagblock) THEN
            DO i = 1, natorb_a
               na = naoa(i)
               fqa = fqa + pblock(i, i)*dqhuckel(na, iatom)
            END DO
            dcharges(iatom) = dcharges(iatom) + fqa
         ELSE
            DO j = 1, natorb_b
               nb = naob(j)
               DO i = 1, natorb_a
                  na = naoa(i)
                  hij = 0.5_dp*pia(na)*pib(nb)
                  drx = f2*hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)
                  IF (iatom <= jatom) THEN
                     fqa = fqa + drx*pblock(i, j)*dqhuckel(na, iatom)
                     fqb = fqb + drx*pblock(i, j)*dqhuckel(nb, jatom)
                  ELSE
                     fqa = fqa + drx*pblock(j, i)*dqhuckel(na, iatom)
                     fqb = fqb + drx*pblock(j, i)*dqhuckel(nb, jatom)
                  END IF
               END DO
            END DO
            dcharges(iatom) = dcharges(iatom) + fqa
            dcharges(jatom) = dcharges(jatom) + fqb
         END IF

         DEALLOCATE (oint, owork, sint)

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! EEQ response
      CALL para_env%sum(dcharges)
      CALL xtb_eeq_lagrange(qs_env, dcharges, qresp, eeq_sparam)

      ! deallocate coordination numbers
      CALL cnumber_release(cnumbers, dcnum, .TRUE.)

      ! deallocate Huckel parameters
      DEALLOCATE (huckel, dhuckel, dqhuckel)
      ! deallocate KAB parameters
      DEALLOCATE (kijab)

      ! deallocate charges
      DEALLOCATE (charges, dcharges)

      DEALLOCATE (basis_set_list)
      IF (SIZE(matrix_p, 1) == 2) THEN
         DO img = 1, nimg
            CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, alpha_scalar=1.0_dp, &
                           beta_scalar=-1.0_dp)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_gfn0_qresp

END MODULE xtb_qresp

