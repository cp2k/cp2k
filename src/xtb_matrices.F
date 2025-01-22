!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Overlap and Hamiltonian matrices in xTB
!>        Reference: Stefan Grimme, Christoph Bannwarth, Philip Shushkov
!>                   JCTC 13, 1989-2009, (2017)
!>                   DOI: 10.1021/acs.jctc.7b00118
!> \author JGH
! **************************************************************************************************
MODULE xtb_matrices
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_overlap,                      ONLY: overlap_ab
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE eeq_input,                       ONLY: eeq_solver_type
   USE input_constants,                 ONLY: vdw_pairpot_dftd4
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_condnum,                      ONLY: overlap_condnum
   USE qs_dispersion_cnum,              ONLY: cnumber_init,&
                                              cnumber_release,&
                                              dcnum_type
   USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                              get_memory_usage
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: create_sab_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_eeq,                         ONLY: xtb_eeq_calculation,&
                                              xtb_eeq_forces
   USE xtb_hcore,                       ONLY: gfn0_huckel,&
                                              gfn0_kpair,&
                                              gfn1_huckel,&
                                              gfn1_kpair
   USE xtb_potentials,                  ONLY: nonbonded_correction,&
                                              repulsive_potential,&
                                              srb_potential,&
                                              xb_interaction
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_matrices'

   PUBLIC :: build_xtb_matrices

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_xtb_matrices(qs_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: gfn_type
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      gfn_type = dft_control%qs_control%xtb_control%gfn_type

      SELECT CASE (gfn_type)
      CASE (0)
         CALL build_gfn0_xtb_matrices(qs_env, calculate_forces)
      CASE (1)
         CALL build_gfn1_xtb_matrices(qs_env, calculate_forces)
      CASE (2)
         CPABORT("gfn_type = 2 not yet available")
      CASE DEFAULT
         CPABORT("Unknown gfn_type")
      END SELECT

   END SUBROUTINE build_xtb_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_gfn0_xtb_matrices(qs_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_gfn0_xtb_matrices'

      INTEGER :: atom_a, atom_b, atom_c, handle, i, iatom, ic, icol, ikind, img, ir, irow, iset, &
         j, jatom, jkind, jset, katom, kkind, la, lb, ldsab, lmaxa, lmaxb, maxder, n1, n2, na, &
         natom, natorb_a, natorb_b, nb, ncoa, ncob, nderivatives, nimg, nkind, nsa, nsb, nseta, &
         nsetb, sgfa, sgfb, za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(25)                             :: laoa, laob, naoa, naob
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, diagblock, do_nonbonded, found, &
                                                            use_virial
      REAL(KIND=dp) :: dfp, dhij, dr, drk, drx, eeq_energy, ef_energy, enonbonded, enscale, erep, &
         esrb, etaa, etab, f0, f1, f2, fhua, fhub, fhud, foab, fqa, fqb, hij, kf, qlambda, rcova, &
         rcovab, rcovb, rrab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: charges, cnumbers, dcharges, qlagrange
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dfblock, dhuckel, dqhuckel, huckel, owork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: kijab
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, fdika, fdikb, force_ab, rij, rik
      REAL(KIND=dp), DIMENSION(5)                        :: dpia, dpib, hena, henb, kpolya, kpolyb, &
                                                            pia, pib
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eeq_q, set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fblock, pblock, rpgfa, rpgfb, sblock, &
                                                            scon_a, scon_b, wblock, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(block_p_type), DIMENSION(2:4)                 :: dsblocks
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, matrix_s, matrix_w
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
         POINTER                                         :: sab_orb, sab_xtb_nonbond
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      NULLIFY (logger, virial, atprop)
      logger => cp_get_default_logger()

      NULLIFY (matrix_h, matrix_s, matrix_p, matrix_w, atomic_kind_set, &
               qs_kind_set, sab_orb, ks_env)
      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      matrix_s_kp=matrix_s, &
                      para_env=para_env, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb)

      nkind = SIZE(atomic_kind_set)
      xtb_control => dft_control%qs_control%xtb_control
      do_nonbonded = xtb_control%do_nonbonded
      nimg = dft_control%nimages
      nderivatives = 0
      IF (calculate_forces) nderivatives = 1
      IF (dft_control%tddfpt2_control%enabled) nderivatives = 1
      maxder = ncoset(nderivatives)

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      IF (calculate_forces) THEN
         NULLIFY (rho, force, matrix_w)
         CALL get_qs_env(qs_env=qs_env, &
                         rho=rho, matrix_w_kp=matrix_w, &
                         virial=virial, force=force)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_add(matrix_w(1, img)%matrix, matrix_w(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF
      ! atomic energy decomposition
      IF (atprop%energy) THEN
         CALL atprop_array_init(atprop%atecc, natom)
      END IF

      NULLIFY (cell_to_index)
      IF (nimg > 1) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      ! set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ! allocate overlap matrix
      CALL dbcsr_allocate_matrix_set(matrix_s, maxder, nimg)
      CALL create_sab_matrix(ks_env, matrix_s, "xTB OVERLAP MATRIX", basis_set_list, basis_set_list, &
                             sab_orb, .TRUE.)
      CALL set_ks_env(ks_env, matrix_s_kp=matrix_s)

      ! initialize H matrix
      CALL dbcsr_allocate_matrix_set(matrix_h, 1, nimg)
      DO img = 1, nimg
         ALLOCATE (matrix_h(1, img)%matrix)
         CALL dbcsr_create(matrix_h(1, img)%matrix, template=matrix_s(1, 1)%matrix, &
                           name="HAMILTONIAN MATRIX")
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_h(1, img)%matrix, sab_orb)
      END DO
      CALL set_ks_env(ks_env, matrix_h_kp=matrix_h)

      ! Calculate coordination numbers
      ! needed for effective atomic energy levels
      ! code taken from D3 dispersion energy
      CALL cnumber_init(qs_env, cnumbers, dcnum, 2, calculate_forces)

      ALLOCATE (charges(natom))
      charges = 0.0_dp
      CALL xtb_eeq_calculation(qs_env, charges, cnumbers, eeq_sparam, eeq_energy, ef_energy, qlambda)
      IF (calculate_forces) THEN
         ALLOCATE (dcharges(natom))
         dcharges = qlambda/REAL(para_env%num_pe, KIND=dp)
      END IF
      energy%eeq = eeq_energy
      energy%efield = ef_energy

      CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env)
      ! prepare charges (needed for D4)
      IF (dispersion_env%pp_type == vdw_pairpot_dftd4) THEN
         dispersion_env%ext_charges = .TRUE.
         IF (ASSOCIATED(dispersion_env%charges)) DEALLOCATE (dispersion_env%charges)
         ALLOCATE (dispersion_env%charges(natom))
         dispersion_env%charges = charges
         IF (calculate_forces) THEN
            IF (ASSOCIATED(dispersion_env%dcharges)) DEALLOCATE (dispersion_env%dcharges)
            ALLOCATE (dispersion_env%dcharges(natom))
            dispersion_env%dcharges = 0.0_dp
         END IF
      END IF
      CALL calculate_dispersion_pairpot(qs_env, dispersion_env, &
                                        energy%dispersion, calculate_forces)
      IF (calculate_forces) THEN
         IF (dispersion_env%pp_type == vdw_pairpot_dftd4 .AND. dispersion_env%ext_charges) THEN
            dcharges(1:natom) = dcharges(1:natom) + dispersion_env%dcharges(1:natom)
         END IF
      END IF

      ! Calculate Huckel parameters
      CALL gfn0_huckel(qs_env, cnumbers, charges, huckel, dhuckel, dqhuckel, calculate_forces)

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
         NULLIFY (sblock, fblock)
         CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=matrix_h(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=fblock, found=found)
         CPASSERT(found)

         IF (calculate_forces) THEN
            NULLIFY (pblock)
            CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                   row=irow, col=icol, block=pblock, found=found)
            CPASSERT(ASSOCIATED(pblock))
            NULLIFY (wblock)
            CALL dbcsr_get_block_p(matrix=matrix_w(1, ic)%matrix, &
                                   row=irow, col=icol, block=wblock, found=found)
            CPASSERT(ASSOCIATED(wblock))
            DO i = 2, 4
               NULLIFY (dsblocks(i)%block)
               CALL dbcsr_get_block_p(matrix=matrix_s(i, ic)%matrix, &
                                      row=irow, col=icol, BLOCK=dsblocks(i)%block, found=found)
               CPASSERT(found)
            END DO
         END IF

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
               IF (calculate_forces) THEN
                  CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                  rij, sab=oint(:, :, 1), dab=oint(:, :, 2:4))
               ELSE
                  CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                  rij, sab=oint(:, :, 1))
               END IF
               ! Contraction
               CALL contraction(oint(:, :, 1), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
               CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, 1), sgfa, sgfb, trans=.FALSE.)
               IF (calculate_forces) THEN
                  DO i = 2, 4
                     CALL contraction(oint(:, :, i), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                      cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
                     CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, i), sgfa, sgfb, trans=.FALSE.)
                  END DO
               END IF
            END DO
         END DO
         ! forces W matrix
         IF (calculate_forces) THEN
            DO i = 1, 3
               IF (iatom <= jatom) THEN
                  force_ab(i) = SUM(sint(:, :, i + 1)*wblock(:, :))
               ELSE
                  force_ab(i) = SUM(sint(:, :, i + 1)*TRANSPOSE(wblock(:, :)))
               END IF
            END DO
            f1 = 2.0_dp
            force(ikind)%overlap(:, atom_a) = force(ikind)%overlap(:, atom_a) - f1*force_ab(:)
            force(jkind)%overlap(:, atom_b) = force(jkind)%overlap(:, atom_b) + f1*force_ab(:)
            IF (use_virial .AND. dr > 1.e-3_dp) THEN
               IF (iatom == jatom) f1 = 1.0_dp
               CALL virial_pair_force(virial%pv_virial, -f1, force_ab, rij)
            END IF
         END IF
         ! update S matrix
         IF (iatom <= jatom) THEN
            sblock(:, :) = sblock(:, :) + sint(:, :, 1)
         ELSE
            sblock(:, :) = sblock(:, :) + TRANSPOSE(sint(:, :, 1))
         END IF
         IF (calculate_forces) THEN
            DO i = 2, 4
               IF (iatom <= jatom) THEN
                  dsblocks(i)%block(:, :) = dsblocks(i)%block(:, :) + sint(:, :, i)
               ELSE
                  dsblocks(i)%block(:, :) = dsblocks(i)%block(:, :) - TRANSPOSE(sint(:, :, i))
               END IF
            END DO
         END IF

         ! Calculate Pi = Pia * Pib (Eq. 11)
         rcovab = rcova + rcovb
         rrab = SQRT(dr/rcovab)
         pia(1:nsa) = 1._dp + kpolya(1:nsa)*rrab
         pib(1:nsb) = 1._dp + kpolyb(1:nsb)*rrab
         IF (calculate_forces) THEN
            IF (dr > 1.e-6_dp) THEN
               drx = 0.5_dp/rrab/rcovab
            ELSE
               drx = 0.0_dp
            END IF
            dpia(1:nsa) = drx*kpolya(1:nsa)
            dpib(1:nsb) = drx*kpolyb(1:nsb)
         END IF

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
         IF (calculate_forces) THEN
            f0 = 1.0_dp
            IF (irow == iatom) f0 = -1.0_dp
            f2 = 1.0_dp
            IF (iatom /= jatom) f2 = 2.0_dp
            ! Derivative wrt coordination number
            fhua = 0.0_dp
            fhub = 0.0_dp
            fhud = 0.0_dp
            fqa = 0.0_dp
            fqb = 0.0_dp
            IF (diagblock) THEN
               DO i = 1, natorb_a
                  la = laoa(i)
                  na = naoa(i)
                  fhud = fhud + pblock(i, i)*dhuckel(na, iatom)
                  fqa = fqa + pblock(i, i)*dqhuckel(na, iatom)
               END DO
               dcharges(iatom) = dcharges(iatom) + fqa
            ELSE
               DO j = 1, natorb_b
                  lb = laob(j)
                  nb = naob(j)
                  DO i = 1, natorb_a
                     la = laoa(i)
                     na = naoa(i)
                     hij = 0.5_dp*pia(na)*pib(nb)
                     drx = f2*hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)
                     IF (iatom <= jatom) THEN
                        fhua = fhua + drx*pblock(i, j)*dhuckel(na, iatom)
                        fhub = fhub + drx*pblock(i, j)*dhuckel(nb, jatom)
                        fqa = fqa + drx*pblock(i, j)*dqhuckel(na, iatom)
                        fqb = fqb + drx*pblock(i, j)*dqhuckel(nb, jatom)
                     ELSE
                        fhua = fhua + drx*pblock(j, i)*dhuckel(na, iatom)
                        fhub = fhub + drx*pblock(j, i)*dhuckel(nb, jatom)
                        fqa = fqa + drx*pblock(j, i)*dqhuckel(na, iatom)
                        fqb = fqb + drx*pblock(j, i)*dqhuckel(nb, jatom)
                     END IF
                  END DO
               END DO
               dcharges(iatom) = dcharges(iatom) + fqa
               dcharges(jatom) = dcharges(jatom) + fqb
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
                  IF (use_virial) THEN
                     fdik = fdika + fdikb
                     CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                  END IF
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
                  IF (use_virial) THEN
                     CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                  END IF
               END IF
            END DO
            ! force from R dendent Huckel element: Pia*Pib
            IF (diagblock) THEN
               force_ab = 0._dp
            ELSE
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
         END IF

         IF (calculate_forces) THEN
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            IF (irow == iatom) force_ab = -force_ab
            force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a) - force_ab(:)
            force(jkind)%all_potential(:, atom_b) = force(jkind)%all_potential(:, atom_b) + force_ab(:)
            IF (use_virial) THEN
               f1 = 1.0_dp
               IF (iatom == jatom) f1 = 0.5_dp
               CALL virial_pair_force(virial%pv_virial, -f1, force_ab, rij)
            END IF
         END IF

         DEALLOCATE (oint, owork, sint)

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(matrix_h, 1)
         DO img = 1, nimg
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
      END DO

      ! EEQ forces (response and direct)
      IF (calculate_forces) THEN
         CALL para_env%sum(dcharges)
         ALLOCATE (qlagrange(natom))
         CALL xtb_eeq_forces(qs_env, charges, dcharges, qlagrange, cnumbers, dcnum, eeq_sparam)
      END IF

      kf = xtb_control%kf
      enscale = xtb_control%enscale
      erep = 0.0_dp
      CALL repulsive_potential(qs_env, erep, kf, enscale, calculate_forces)

      esrb = 0.0_dp
      CALL srb_potential(qs_env, esrb, calculate_forces, xtb_control, cnumbers, dcnum)

      enonbonded = 0.0_dp
      IF (do_nonbonded) THEN
         ! nonbonded interactions
         NULLIFY (sab_xtb_nonbond)
         CALL get_qs_env(qs_env=qs_env, sab_xtb_nonbond=sab_xtb_nonbond)
         CALL nonbonded_correction(enonbonded, force, qs_env, xtb_control, sab_xtb_nonbond, &
                                   atomic_kind_set, calculate_forces, use_virial, virial, atprop, atom_of_kind)
      END IF

      ! set repulsive energy
      erep = erep + esrb + enonbonded
      IF (do_nonbonded) THEN
         CALL para_env%sum(enonbonded)
         energy%xtb_nonbonded = enonbonded
      END IF
      CALL para_env%sum(esrb)
      energy%srb = esrb
      CALL para_env%sum(erep)
      energy%repulsive = erep

      ! save EEQ charges
      NULLIFY (eeq_q)
      CALL get_qs_env(qs_env, eeq=eeq_q)
      IF (ASSOCIATED(eeq_q)) THEN
         CPASSERT(SIZE(eeq_q) == natom)
      ELSE
         ALLOCATE (eeq_q(natom))
         eeq_q(1:natom) = charges(1:natom)
      END IF
      CALL set_qs_env(qs_env, eeq=eeq_q)

      ! deallocate coordination numbers
      CALL cnumber_release(cnumbers, dcnum, calculate_forces)

      ! deallocate Huckel parameters
      DEALLOCATE (huckel)
      IF (calculate_forces) THEN
         DEALLOCATE (dhuckel, dqhuckel)
      END IF
      ! deallocate KAB parameters
      DEALLOCATE (kijab)

      ! deallocate charges
      DEALLOCATE (charges)
      IF (calculate_forces) THEN
         DEALLOCATE (dcharges, qlagrange)
      END IF

      ! AO matrix outputs
      CALL ao_matrix_output(qs_env, matrix_h, matrix_s, calculate_forces)

      DEALLOCATE (basis_set_list)
      IF (calculate_forces) THEN
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, alpha_scalar=1.0_dp, &
                              beta_scalar=-1.0_dp)
            END DO
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_gfn0_xtb_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_gfn1_xtb_matrices(qs_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_gfn1_xtb_matrices'

      INTEGER :: atom_a, atom_b, atom_c, handle, i, iatom, ic, icol, ikind, img, ir, irow, iset, &
         j, jatom, jkind, jset, katom, kkind, la, lb, ldsab, lmaxa, lmaxb, maxder, n1, n2, na, &
         natom, natorb_a, natorb_b, nb, ncoa, ncob, nderivatives, nimg, nkind, nsa, nsb, nseta, &
         nsetb, sgfa, sgfb, za, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(25)                             :: laoa, laob, naoa, naob
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, diagblock, do_nonbonded, found, &
                                                            use_virial, xb_inter
      REAL(KIND=dp)                                      :: dfp, dhij, dr, drk, drx, enonbonded, &
                                                            enscale, erep, etaa, etab, exb, f0, &
                                                            f1, fhua, fhub, fhud, foab, hij, kf, &
                                                            rcova, rcovab, rcovb, rrab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cnumbers
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dfblock, dhuckel, huckel, owork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: kijab
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, fdika, fdikb, force_ab, rij, rik
      REAL(KIND=dp), DIMENSION(5)                        :: dpia, dpib, kpolya, kpolyb, pia, pib
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fblock, pblock, rpgfa, rpgfb, sblock, &
                                                            scon_a, scon_b, wblock, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(block_p_type), DIMENSION(2:4)                 :: dsblocks
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, matrix_s, matrix_w
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sab_xtb_nonbond
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b
      TYPE(xtb_control_type), POINTER                    :: xtb_control

      CALL timeset(routineN, handle)

      NULLIFY (logger, virial, atprop)
      logger => cp_get_default_logger()

      NULLIFY (matrix_h, matrix_s, matrix_p, matrix_w, atomic_kind_set, &
               qs_kind_set, sab_orb, ks_env)

      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      matrix_s_kp=matrix_s, &
                      para_env=para_env, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb)

      nkind = SIZE(atomic_kind_set)
      xtb_control => dft_control%qs_control%xtb_control
      xb_inter = xtb_control%xb_interaction
      do_nonbonded = xtb_control%do_nonbonded
      nimg = dft_control%nimages
      nderivatives = 0
      IF (calculate_forces) nderivatives = 1
      IF (dft_control%tddfpt2_control%enabled) nderivatives = 1
      maxder = ncoset(nderivatives)

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      IF (calculate_forces) THEN
         NULLIFY (rho, force, matrix_w)
         CALL get_qs_env(qs_env=qs_env, &
                         rho=rho, matrix_w_kp=matrix_w, &
                         virial=virial, force=force)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_add(matrix_w(1, img)%matrix, matrix_w(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF
      ! atomic energy decomposition
      IF (atprop%energy) THEN
         CALL atprop_array_init(atprop%atecc, natom)
      END IF

      NULLIFY (cell_to_index)
      IF (nimg > 1) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      ! set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ! allocate overlap matrix
      CALL dbcsr_allocate_matrix_set(matrix_s, maxder, nimg)
      CALL create_sab_matrix(ks_env, matrix_s, "xTB OVERLAP MATRIX", basis_set_list, basis_set_list, &
                             sab_orb, .TRUE.)
      CALL set_ks_env(ks_env, matrix_s_kp=matrix_s)

      ! initialize H matrix
      CALL dbcsr_allocate_matrix_set(matrix_h, 1, nimg)
      DO img = 1, nimg
         ALLOCATE (matrix_h(1, img)%matrix)
         CALL dbcsr_create(matrix_h(1, img)%matrix, template=matrix_s(1, 1)%matrix, &
                           name="HAMILTONIAN MATRIX")
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_h(1, img)%matrix, sab_orb)
      END DO
      CALL set_ks_env(ks_env, matrix_h_kp=matrix_h)

      ! Calculate coordination numbers
      ! needed for effective atomic energy levels (Eq. 12)
      ! code taken from D3 dispersion energy
      CALL cnumber_init(qs_env, cnumbers, dcnum, 1, calculate_forces)

      ! vdW Potential
      CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env)
      CALL calculate_dispersion_pairpot(qs_env, dispersion_env, &
                                        energy%dispersion, calculate_forces)

      ! Calculate Huckel parameters
      CALL gfn1_huckel(qs_env, cnumbers, huckel, dhuckel, calculate_forces)

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
         CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa, lao=laoa, rcov=rcova, eta=etaa, &
                                 lmax=lmaxa, nshell=nsa, kpoly=kpolya)
         CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob, lao=laob, rcov=rcovb, eta=etab, &
                                 lmax=lmaxb, nshell=nsb, kpoly=kpolyb)

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cell(1), cell(2), cell(3))
            CPASSERT(ic > 0)
         END IF

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)
         NULLIFY (sblock, fblock)
         CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=matrix_h(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=fblock, found=found)
         CPASSERT(found)

         IF (calculate_forces) THEN
            NULLIFY (pblock)
            CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                   row=irow, col=icol, block=pblock, found=found)
            CPASSERT(found)
            NULLIFY (wblock)
            CALL dbcsr_get_block_p(matrix=matrix_w(1, ic)%matrix, &
                                   row=irow, col=icol, block=wblock, found=found)
            CPASSERT(found)
            DO i = 2, 4
               NULLIFY (dsblocks(i)%block)
               CALL dbcsr_get_block_p(matrix=matrix_s(i, ic)%matrix, &
                                      row=irow, col=icol, BLOCK=dsblocks(i)%block, found=found)
               CPASSERT(found)
            END DO
         END IF

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
               IF (calculate_forces) THEN
                  CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                  rij, sab=oint(:, :, 1), dab=oint(:, :, 2:4))
               ELSE
                  CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                  rij, sab=oint(:, :, 1))
               END IF
               ! Contraction
               CALL contraction(oint(:, :, 1), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
               CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, 1), sgfa, sgfb, trans=.FALSE.)
               IF (calculate_forces) THEN
                  DO i = 2, 4
                     CALL contraction(oint(:, :, i), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                      cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
                     CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, i), sgfa, sgfb, trans=.FALSE.)
                  END DO
               END IF
            END DO
         END DO
         ! forces W matrix
         IF (calculate_forces) THEN
            DO i = 1, 3
               IF (iatom <= jatom) THEN
                  force_ab(i) = SUM(sint(:, :, i + 1)*wblock(:, :))
               ELSE
                  force_ab(i) = SUM(sint(:, :, i + 1)*TRANSPOSE(wblock(:, :)))
               END IF
            END DO
            f1 = 2.0_dp
            force(ikind)%overlap(:, atom_a) = force(ikind)%overlap(:, atom_a) - f1*force_ab(:)
            force(jkind)%overlap(:, atom_b) = force(jkind)%overlap(:, atom_b) + f1*force_ab(:)
            IF (use_virial .AND. dr > 1.e-3_dp) THEN
               IF (iatom == jatom) f1 = 1.0_dp
               CALL virial_pair_force(virial%pv_virial, -f1, force_ab, rij)
            END IF
         END IF
         ! update S matrix
         IF (iatom <= jatom) THEN
            sblock(:, :) = sblock(:, :) + sint(:, :, 1)
         ELSE
            sblock(:, :) = sblock(:, :) + TRANSPOSE(sint(:, :, 1))
         END IF
         IF (calculate_forces) THEN
            DO i = 2, 4
               IF (iatom <= jatom) THEN
                  dsblocks(i)%block(:, :) = dsblocks(i)%block(:, :) + sint(:, :, i)
               ELSE
                  dsblocks(i)%block(:, :) = dsblocks(i)%block(:, :) - TRANSPOSE(sint(:, :, i))
               END IF
            END DO
         END IF

         ! Calculate Pi = Pia * Pib (Eq. 11)
         rcovab = rcova + rcovb
         rrab = SQRT(dr/rcovab)
         pia(1:nsa) = 1._dp + kpolya(1:nsa)*rrab
         pib(1:nsb) = 1._dp + kpolyb(1:nsb)*rrab
         IF (calculate_forces) THEN
            IF (dr > 1.e-6_dp) THEN
               drx = 0.5_dp/rrab/rcovab
            ELSE
               drx = 0.0_dp
            END IF
            dpia(1:nsa) = drx*kpolya(1:nsa)
            dpib(1:nsb) = drx*kpolyb(1:nsb)
         END IF

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
         IF (calculate_forces) THEN
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
                  IF (use_virial) THEN
                     fdik = fdika + fdikb
                     CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                  END IF
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
                  IF (use_virial) THEN
                     CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
                  END IF
               END IF
            END DO
            ! force from R dendent Huckel element: Pia*Pib
            IF (diagblock) THEN
               force_ab = 0._dp
            ELSE
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
         END IF

         IF (calculate_forces) THEN
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            IF (irow == iatom) force_ab = -force_ab
            force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a) - force_ab(:)
            force(jkind)%all_potential(:, atom_b) = force(jkind)%all_potential(:, atom_b) + force_ab(:)
            IF (use_virial) THEN
               f1 = 1.0_dp
               IF (iatom == jatom) f1 = 0.5_dp
               CALL virial_pair_force(virial%pv_virial, -f1, force_ab, rij)
            END IF
         END IF

         DEALLOCATE (oint, owork, sint)

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(matrix_h, 1)
         DO img = 1, nimg
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
      END DO

      kf = xtb_control%kf
      enscale = xtb_control%enscale
      erep = 0.0_dp
      CALL repulsive_potential(qs_env, erep, kf, enscale, calculate_forces)

      exb = 0.0_dp
      IF (xb_inter) THEN
         CALL xb_interaction(qs_env, exb, calculate_forces)
      END IF

      enonbonded = 0.0_dp
      IF (do_nonbonded) THEN
         ! nonbonded interactions
         NULLIFY (sab_xtb_nonbond)
         CALL get_qs_env(qs_env=qs_env, sab_xtb_nonbond=sab_xtb_nonbond)
         CALL nonbonded_correction(enonbonded, force, qs_env, xtb_control, sab_xtb_nonbond, &
                                   atomic_kind_set, calculate_forces, use_virial, virial, atprop, atom_of_kind)
      END IF

      ! set repulsive energy
      erep = erep + exb + enonbonded
      IF (xb_inter) THEN
         CALL para_env%sum(exb)
         energy%xtb_xb_inter = exb
      END IF
      IF (do_nonbonded) THEN
         CALL para_env%sum(enonbonded)
         energy%xtb_nonbonded = enonbonded
      END IF
      CALL para_env%sum(erep)
      energy%repulsive = erep

      ! deallocate coordination numbers
      CALL cnumber_release(cnumbers, dcnum, calculate_forces)

      ! deallocate Huckel parameters
      DEALLOCATE (huckel)
      IF (calculate_forces) THEN
         DEALLOCATE (dhuckel)
      END IF
      ! deallocate KAB parameters
      DEALLOCATE (kijab)

      ! AO matrix outputs
      CALL ao_matrix_output(qs_env, matrix_h, matrix_s, calculate_forces)

      DEALLOCATE (basis_set_list)
      IF (calculate_forces) THEN
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, alpha_scalar=1.0_dp, &
                              beta_scalar=-1.0_dp)
            END DO
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_gfn1_xtb_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE ao_matrix_output(qs_env, matrix_h, matrix_s, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_s
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: after, i, img, iw, nimg
      LOGICAL                                            :: norml1, norml2, omit_headers, use_arnoldi
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env

      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env, para_env=para_env)
      nimg = SIZE(matrix_h, 2)
      CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           qs_env%input, "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN", &
                                   extension=".Log")
         CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
         after = MIN(MAX(after, 1), 16)
         DO img = 1, nimg
            CALL cp_dbcsr_write_sparse_matrix(matrix_h(1, img)%matrix, 4, after, qs_env, para_env, &
                                              output_unit=iw, omit_headers=omit_headers)
         END DO
         CALL cp_print_key_finished_output(iw, logger, qs_env%input, "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           qs_env%input, "DFT%PRINT%AO_MATRICES/OVERLAP"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%AO_MATRICES/OVERLAP", &
                                   extension=".Log")
         CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
         after = MIN(MAX(after, 1), 16)
         DO img = 1, nimg
            CALL cp_dbcsr_write_sparse_matrix(matrix_s(1, img)%matrix, 4, after, qs_env, para_env, &
                                              output_unit=iw, omit_headers=omit_headers)
            IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                 qs_env%input, "DFT%PRINT%AO_MATRICES/DERIVATIVES"), cp_p_file)) THEN
               DO i = 2, SIZE(matrix_s, 1)
                  CALL cp_dbcsr_write_sparse_matrix(matrix_s(i, img)%matrix, 4, after, qs_env, para_env, &
                                                    output_unit=iw, omit_headers=omit_headers)
               END DO
            END IF
         END DO
         CALL cp_print_key_finished_output(iw, logger, qs_env%input, "DFT%PRINT%AO_MATRICES/OVERLAP")
      END IF

      ! *** Overlap condition number
      IF (.NOT. calculate_forces) THEN
         IF (cp_print_key_should_output(logger%iter_info, qs_env%input, &
                                        "DFT%PRINT%OVERLAP_CONDITION") .NE. 0) THEN
            iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%OVERLAP_CONDITION", &
                                      extension=".Log")
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%1-NORM", l_val=norml1)
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%DIAGONALIZATION", l_val=norml2)
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%ARNOLDI", l_val=use_arnoldi)
            CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env)
            CALL overlap_condnum(matrix_s, condnum, iw, norml1, norml2, use_arnoldi, blacs_env)
         END IF
      END IF

   END SUBROUTINE ao_matrix_output

END MODULE xtb_matrices

