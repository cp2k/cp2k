!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_type
   USE efield_tb_methods,               ONLY: efield_tb_matrix
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_minloc,&
                                              mp_sum
   USE mulliken,                        ONLY: ao_charges
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_charge_mixing,                ONLY: charge_mixing
   USE qs_condnum,                      ONLY: overlap_condnum
   USE qs_dispersion_pairpot,           ONLY: d3_cnumber,&
                                              dcnum_distribute,&
                                              dcnum_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                              get_memory_usage
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: create_sab_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_coulomb,                     ONLY: build_xtb_coulomb
   USE xtb_parameters,                  ONLY: xtb_set_kab
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   TYPE neighbor_atoms_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: coord
      REAL(KIND=dp), DIMENSION(:), POINTER        :: rab
      INTEGER, DIMENSION(:), POINTER              :: katom
   END TYPE neighbor_atoms_type

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_matrices'

   PUBLIC :: build_xtb_matrices, build_xtb_ks_matrix

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_xtb_matrices(qs_env, para_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_xtb_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER :: after, atom_a, atom_b, atom_c, handle, i, iatom, ic, icol, ikind, img, ir, irow, &
         iset, iw, j, jatom, jkind, jset, katom, kkind, la, lb, ldsab, lmaxa, lmaxb, maxder, maxs, &
         n1, n2, na, natom, natorb_a, natorb_b, nb, ncoa, ncob, nderivatives, nimg, nkind, nsa, &
         nsb, nseta, nsetb, nshell, sgfa, sgfb, za, zat, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, atomnumber, kind_of
      INTEGER, DIMENSION(25)                             :: laoa, laob, lval, naoa, naob
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, diagblock, floating_a, found, &
                                                            ghost_a, norml1, norml2, omit_headers, &
                                                            use_arnoldi, use_virial, xb_inter
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: floating, ghost
      REAL(KIND=dp) :: alphaa, alphab, derepij, dfp, dhij, dr, drk, drx, ena, enb, erep, erepij, &
         etaa, etab, exb, f0, f1, fen, fhua, fhub, fhud, foab, hij, k2sh, kab, kcnd, kcnp, kcns, &
         kd, ken, kf, kg, kia, kjb, kp, ks, ksp, kx2, kxr, rcova, rcovab, rcovb, rrab, zneffa, &
         zneffb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cnumbers, kx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dfblock, huckel, kcnlk, owork, rcab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: kijab
      REAL(KIND=dp), DIMENSION(0:3)                      :: kl
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, fdika, fdikb, force_ab, force_rr, &
                                                            rij, rik
      REAL(KIND=dp), DIMENSION(5)                        :: dpia, dpib, hena, henb, kappaa, kappab, &
                                                            kpolya, kpolyb, pia, pib
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fblock, pblock, rpgfa, rpgfb, sblock, &
                                                            scon_a, scon_b, wblock, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(block_p_type), DIMENSION(2:4)                 :: dsblocks
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, matrix_s, matrix_w
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_atoms_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: neighbor_atoms
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sab_xb
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
                      atprop=atprop, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb)

      nkind = SIZE(atomic_kind_set)
      xtb_control => dft_control%qs_control%xtb_control
      xb_inter = xtb_control%xb_interaction
      nimg = dft_control%nimages
      nderivatives = 0
      IF (calculate_forces) nderivatives = 1
      IF (dft_control%tddfpt2_control%enabled) nderivatives = 1
      maxder = ncoset(nderivatives)

      ! global parameters (Table 2 from Ref.)
      ks = xtb_control%ks
      kp = xtb_control%kp
      kd = xtb_control%kd
      ksp = xtb_control%ksp
      k2sh = xtb_control%k2sh
      kg = xtb_control%kg
      kf = xtb_control%kf
      kcns = xtb_control%kcns
      kcnp = xtb_control%kcnp
      kcnd = xtb_control%kcnd
      ken = xtb_control%ken
      kxr = xtb_control%kxr
      kx2 = xtb_control%kx2

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)
      ALLOCATE (atom_of_kind(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxsgf=maxs, basis_type="ORB")

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
      ALLOCATE (cnumbers(natom))
      cnumbers = 0._dp
      IF (calculate_forces) THEN
         ALLOCATE (dcnum(natom))
         dcnum(:)%neighbors = 0
         DO iatom = 1, natom
            ALLOCATE (dcnum(iatom)%nlist(10), dcnum(iatom)%dvals(10), dcnum(iatom)%rik(3, 10))
         END DO
      ELSE
         ALLOCATE (dcnum(1))
      END IF

      ALLOCATE (ghost(nkind), floating(nkind), atomnumber(nkind))
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
         CALL get_qs_kind(qs_kind_set(ikind), ghost=ghost_a, floating=floating_a)
         ghost(ikind) = ghost_a
         floating(ikind) = floating_a
         atomnumber(ikind) = za
      END DO
      CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env)
      CALL d3_cnumber(qs_env, dispersion_env, cnumbers, dcnum, ghost, floating, atomnumber, &
                      calculate_forces, .FALSE.)
      DEALLOCATE (ghost, floating, atomnumber)
      CALL mp_sum(cnumbers, para_env%group)
      IF (calculate_forces) THEN
         CALL dcnum_distribute(dcnum, para_env)
      END IF

      ! Calculate Huckel parameters
      ! Eq 12
      ! huckel(nshell,natom)
      ALLOCATE (kcnlk(0:3, nkind))
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
         IF (metal(za)) THEN
            kcnlk(0:3, ikind) = 0.0_dp
         ELSEIF (early3d(za)) THEN
            kcnlk(0, ikind) = kcns
            kcnlk(1, ikind) = kcnp
            kcnlk(2, ikind) = 0.005_dp
            kcnlk(3, ikind) = 0.0_dp
         ELSE
            kcnlk(0, ikind) = kcns
            kcnlk(1, ikind) = kcnp
            kcnlk(2, ikind) = kcnd
            kcnlk(3, ikind) = 0.0_dp
         END IF
      END DO
      ALLOCATE (huckel(5, natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, nshell=nshell, lval=lval, hen=hena)
         huckel(:, iatom) = 0.0_dp
         DO i = 1, nshell
            huckel(i, iatom) = hena(i)*(1._dp + kcnlk(lval(i), ikind)*cnumbers(iatom))
         END DO
      END DO

      ! Calculate KAB parameters and Huckel parameters and electronegativity correction
      kl(0) = ks
      kl(1) = kp
      kl(2) = kd
      kl(3) = 0.0_dp
      ALLOCATE (kijab(maxs, maxs, nkind, nkind))
      kijab = 0.0_dp
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined, natorb=natorb_a)
         IF (.NOT. defined .OR. natorb_a < 1) CYCLE
         CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa, lao=laoa, electronegativity=ena)
         DO jkind = 1, nkind
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_b, defined=defined, natorb=natorb_b)
            IF (.NOT. defined .OR. natorb_b < 1) CYCLE
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob, lao=laob, electronegativity=enb)
            ! get Fen = (1+ken*deltaEN^2)
            fen = 1.0_dp + ken*(ena - enb)**2
            ! Kab
            kab = xtb_set_kab(za, zb, xtb_control)
            DO j = 1, natorb_b
               lb = laob(j)
               nb = naob(j)
               DO i = 1, natorb_a
                  la = laoa(i)
                  na = naoa(i)
                  kia = kl(la)
                  kjb = kl(lb)
                  IF (zb == 1 .AND. nb == 2) kjb = k2sh
                  IF (za == 1 .AND. na == 2) kia = k2sh
                  IF ((zb == 1 .AND. nb == 2) .OR. (za == 1 .AND. na == 2)) THEN
                     kijab(i, j, ikind, jkind) = 0.5_dp*(kia + kjb)
                  ELSE
                     IF ((la == 0 .AND. lb == 1) .OR. (la == 1 .AND. lb == 0)) THEN
                        kijab(i, j, ikind, jkind) = ksp*kab*fen
                     ELSE
                        kijab(i, j, ikind, jkind) = 0.5_dp*(kia + kjb)*kab*fen
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO

      IF (xb_inter) THEN
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
         !
         CALL get_qs_env(qs_env=qs_env, sab_xb=sab_xb)
      END IF

      ! initialize repulsion energy
      erep = 0._dp

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

         ! neighbor atom for XB term
         IF (xb_inter .AND. (dr > 1.e-3_dp)) THEN
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

         ! atomic parameters
         CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa, lao=laoa, rcov=rcova, eta=etaa, &
                                 lmax=lmaxa, nshell=nsa, alpha=alphaa, zneff=zneffa, kpoly=kpolya, &
                                 kappa=kappaa, hen=hena)
         CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob, lao=laob, rcov=rcovb, eta=etab, &
                                 lmax=lmaxb, nshell=nsb, alpha=alphab, zneff=zneffb, kpoly=kpolyb, &
                                 kappa=kappab, hen=henb)

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
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp*f1, force_ab, rij)
                  CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp*f1, force_ab, rij)
               END IF
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
         DO i = 1, nsa
            pia(i) = 1._dp + kpolya(i)*rrab
         END DO
         DO i = 1, nsb
            pib(i) = 1._dp + kpolyb(i)*rrab
         END DO
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
                  fhud = fhud + pblock(i, i)*kcnlk(la, ikind)*hena(na)
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
                        fhua = fhua + hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)*pblock(i, j)*kcnlk(la, ikind)*hena(na)
                        fhub = fhub + hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)*pblock(i, j)*kcnlk(lb, jkind)*henb(nb)
                     ELSE
                        fhua = fhua + hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)*pblock(j, i)*kcnlk(la, ikind)*hena(na)
                        fhub = fhub + hij*kijab(i, j, ikind, jkind)*sint(i, j, 1)*pblock(j, i)*kcnlk(lb, jkind)*henb(nb)
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
                     IF (atprop%stress) THEN
                        CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
                        CALL virial_pair_force(atprop%atstress(:, :, katom), -0.5_dp, fdik, rik)
                     END IF
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
                     IF (atprop%stress) THEN
                        CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
                        CALL virial_pair_force(atprop%atstress(:, :, katom), -0.5_dp, fdik, rik)
                     END IF
                  END IF
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
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp*f1, force_ab, rij)
                  CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp*f1, force_ab, rij)
               END IF
            END IF
         END IF

         DEALLOCATE (oint, owork, sint)

         ! repulsive potential
         IF (dr > 0.001_dp) THEN
            erepij = zneffa*zneffb/dr*EXP(-SQRT(alphaa*alphab)*dr**kf)
            erep = erep + erepij
            IF (atprop%energy) THEN
               atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*erepij
               atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*erepij
            END IF
            IF (calculate_forces .AND. (iatom /= jatom .OR. dr > 0.001_dp)) THEN
               derepij = -(1.0_dp/dr + SQRT(alphaa*alphab)*kf*dr**(kf - 1.0_dp))*erepij
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
                  IF (atprop%stress) THEN
                     CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp*f1, force_rr, rij)
                     CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp*f1, force_rr, rij)
                  END IF
               END IF
            END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(matrix_h, 1)
         DO img = 1, nimg
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
      ENDDO

      exb = 0.0_dp
      IF (xb_inter) THEN
         CALL xb_neighbors(neighbor_atoms, para_env)
         CALL xb_interaction(exb, neighbor_atoms, atom_of_kind, kind_of, sab_xb, kx, kx2, rcab, &
                             calculate_forces, use_virial, force, virial, atprop)
      END IF
      ! set repulsive energy
      erep = erep + exb
      CALL mp_sum(erep, para_env%group)
      energy%repulsive = erep

      ! deallocate coordination numbers
      DEALLOCATE (cnumbers)
      IF (calculate_forces) THEN
         DO iatom = 1, natom
            DEALLOCATE (dcnum(iatom)%nlist, dcnum(iatom)%dvals, dcnum(iatom)%rik)
         END DO
      END IF
      DEALLOCATE (dcnum)

      ! deallocate Huckel parameters
      DEALLOCATE (huckel)
      ! deallocate KAB parameters
      DEALLOCATE (kijab)

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
            CALL overlap_condnum(matrix_s, iw, norml1, norml2, use_arnoldi, blacs_env)
         END IF
      END IF

      DEALLOCATE (basis_set_list)
      DEALLOCATE (kind_of)
      DEALLOCATE (atom_of_kind)
      IF (calculate_forces) THEN
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, alpha_scalar=1.0_dp, &
                              beta_scalar=-1.0_dp)
            END DO
         END IF
      END IF

      IF (xb_inter) THEN
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
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_xtb_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param ext_ks_matrix ...
! **************************************************************************************************
   SUBROUTINE build_xtb_ks_matrix(qs_env, calculate_forces, just_energy, ext_ks_matrix)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: ext_ks_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'build_xtb_ks_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, iatom, ikind, img, &
                                                            iounit, is, ispin, na, natom, natorb, &
                                                            nimg, nkind, ns, nsgf, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      LOGICAL                                            :: pass_check
      REAL(KIND=dp)                                      :: achrg, chmax, pc_ener, qmmm_el
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mcharge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, charges
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p1, mo_derivs, p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix, matrix_h, matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: mo_coeff, s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: scf_section
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, logger, scf_section, matrix_p, particle_set, ks_env, &
               ks_matrix, rho, energy)
      CPASSERT(ASSOCIATED(qs_env))

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      para_env=para_env, &
                      ks_env=ks_env, &
                      matrix_ks_kp=ks_matrix, &
                      rho=rho, &
                      energy=energy)

      IF (PRESENT(ext_ks_matrix)) THEN
         ! remap pointer to allow for non-kpoint external ks matrix
         ! ext_ks_matrix is used in linear response code
         ns = SIZE(ext_ks_matrix)
         ks_matrix(1:ns, 1:1) => ext_ks_matrix(1:ns)
      END IF

      energy%hartree = 0.0_dp
      energy%qmmm_el = 0.0_dp

      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")
      nspins = dft_control%nspins
      nimg = dft_control%nimages
      CPASSERT(ASSOCIATED(matrix_h))
      CPASSERT(ASSOCIATED(rho))
      CPASSERT(SIZE(ks_matrix) > 0)

      DO ispin = 1, nspins
         DO img = 1, nimg
            ! copy the core matrix into the fock matrix
            CALL dbcsr_copy(ks_matrix(ispin, img)%matrix, matrix_h(1, img)%matrix)
         END DO
      END DO

      IF (dft_control%qs_control%xtb_control%coulomb_interaction) THEN
         ! Mulliken charges
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, matrix_s_kp=matrix_s)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         natom = SIZE(particle_set)
         ALLOCATE (mcharge(natom), charges(natom, 5))
         charges = 0.0_dp
         nkind = SIZE(atomic_kind_set)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
         ALLOCATE (aocg(nsgf, natom))
         aocg = 0.0_dp
         IF (nimg > 1) THEN
            CALL ao_charges(matrix_p, matrix_s, aocg, para_env)
         ELSE
            p_matrix => matrix_p(:, 1)
            s_matrix => matrix_s(1, 1)%matrix
            CALL ao_charges(p_matrix, s_matrix, aocg, para_env)
         END IF
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, occupation=occ)
            DO iatom = 1, na
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               charges(atom_a, :) = REAL(occ(:), KIND=dp)
               DO is = 1, natorb
                  ns = lao(is) + 1
                  charges(atom_a, ns) = charges(atom_a, ns) - aocg(is, atom_a)
               END DO
            END DO
         END DO
         DEALLOCATE (aocg)

         ! charge mixing
         IF (dft_control%qs_control%do_ls_scf) THEN
            !
         ELSE
            CALL get_qs_env(qs_env=qs_env, scf_env=scf_env)
            CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                               charges, para_env, scf_env%iter_count)
         END IF

         DO iatom = 1, natom
            mcharge(iatom) = SUM(charges(iatom, :))
         END DO

         CALL build_xtb_coulomb(qs_env, ks_matrix, rho, charges, mcharge, energy, &
                                calculate_forces, just_energy)

         CALL efield_tb_matrix(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)

         IF (dft_control%qs_control%xtb_control%check_atomic_charges) THEN
            pass_check = .TRUE.
            DO ikind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
               CALL get_xtb_atom_param(xtb_kind, chmax=chmax)
               DO iatom = 1, na
                  atom_a = atomic_kind_set(ikind)%atom_list(iatom)
                  achrg = mcharge(atom_a)
                  IF (ABS(achrg) > chmax) THEN
                     IF (iounit > 0) THEN
                        WRITE (iounit, "(A,A,I3,I6,A,F4.2,A,F6.2)") " Charge outside chemical range:", &
                           "  Kind Atom=", ikind, atom_a, "   Limit=", chmax, "  Charge=", achrg
                     END IF
                     pass_check = .FALSE.
                  END IF
               END DO
            END DO
            IF (.NOT. pass_check) THEN
               CALL cp_warn(__LOCATION__, "Atomic charges outside chemical range were detected."// &
                            " Switch-off CHECK_ATOMIC_CHARGES keyword in the &xTB section"// &
                            " if you want to force to continue the calculation.")
               CPABORT("xTB Charges")
            END IF
         END IF

         DEALLOCATE (mcharge, charges)
      END IF

      IF (qs_env%qmmm) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         DO ispin = 1, nspins
            ! If QM/MM sumup the 1el Hamiltonian
            CALL dbcsr_add(ks_matrix(ispin, 1)%matrix, qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           1.0_dp, 1.0_dp)
            CALL qs_rho_get(rho, rho_ao=matrix_p1)
            ! Compute QM/MM Energy
            CALL dbcsr_dot(qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           matrix_p1(ispin)%matrix, qmmm_el)
            energy%qmmm_el = energy%qmmm_el + qmmm_el
         END DO
         pc_ener = qs_env%ks_qmmm_env%pc_ener
         energy%qmmm_el = energy%qmmm_el + pc_ener
      END IF

      energy%total = energy%core + energy%hartree + energy%efield + energy%qmmm_el + &
                     energy%repulsive + energy%dispersion + energy%dftb3

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Charge fluctuation energy:                     ", energy%hartree, &
            "London dispersion energy:                      ", energy%dispersion
         IF (ABS(energy%efield) > 1.e-12_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Electric field interaction energy:             ", energy%efield
         END IF
         WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "DFTB3 3rd Order Energy Correction              ", energy%dftb3
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
         END IF
      END IF
      CALL cp_print_key_finished_output(iounit, logger, scf_section, &
                                        "PRINT%DETAILED_ENERGY")
      ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C
      IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         CALL get_qs_env(qs_env, mo_derivs=mo_derivs, mos=mo_array)
         DO ispin = 1, SIZE(mo_derivs)
            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, mo_coeff_b=mo_coeff)
            IF (.NOT. mo_array(ispin)%mo_set%use_mo_coeff_b) THEN
               CPABORT("")
            ENDIF
            CALL dbcsr_multiply('n', 'n', 1.0_dp, ks_matrix(ispin, 1)%matrix, mo_coeff, &
                                0.0_dp, mo_derivs(ispin)%matrix)
         ENDDO
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE build_xtb_ks_matrix

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
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'xb_neighbors', routineP = moduleN//':'//routineN

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
            CALL mp_minloc(dmloc, para_env%group)
            coord = 0.0_dp
            matom = 0
            DO iatom = 1, natom
               neighbor_atoms(ikind)%rab(iatom) = dmloc(2*iatom - 1)
               IF (NINT(dmloc(2*iatom)) == para_env%mepos) THEN
                  coord(1:3, iatom) = neighbor_atoms(ikind)%coord(1:3, iatom)
                  matom(iatom) = neighbor_atoms(ikind)%katom(iatom)
               END IF
            END DO
            CALL mp_sum(coord, para_env%group)
            neighbor_atoms(ikind)%coord(1:3, :) = coord(1:3, :)
            CALL mp_sum(matom, para_env%group)
            neighbor_atoms(ikind)%katom(:) = matom(:)
            DEALLOCATE (dmloc, matom, coord)
         END IF
      END DO

   END SUBROUTINE xb_neighbors
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
   SUBROUTINE xb_interaction(exb, neighbor_atoms, atom_of_kind, kind_of, sab_xb, kx, kx2, rcab, &
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

      CHARACTER(len=*), PARAMETER :: routineN = 'xb_interaction', routineP = moduleN//':'//routineN

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
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fij, rij)
                  CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp, fij, rij)
               END IF
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
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fij, rij)
                  CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp, fij, rij)
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fia, ria)
                  CALL virial_pair_force(atprop%atstress(:, :, katom), -0.5_dp, fia, ria)
                  CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp, fja, rja)
                  CALL virial_pair_force(atprop%atstress(:, :, katom), -0.5_dp, fja, rja)
               END IF
            END IF
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

   END SUBROUTINE xb_interaction

! **************************************************************************************************
!> \brief ...
!> \param z ...
!> \return ...
! **************************************************************************************************
   FUNCTION metal(z) RESULT(ismetal)
      INTEGER                                            :: z
      LOGICAL                                            :: ismetal

      SELECT CASE (z)
      CASE DEFAULT
         ismetal = .TRUE.
      CASE (1:2, 6:10, 14:18, 32:36, 50:54, 82:86)
         ismetal = .FALSE.
      END SELECT

   END FUNCTION metal

! **************************************************************************************************
!> \brief ...
!> \param z ...
!> \return ...
! **************************************************************************************************
   FUNCTION early3d(z) RESULT(isearly3d)
      INTEGER                                            :: z
      LOGICAL                                            :: isearly3d

      isearly3d = .FALSE.
      IF (z >= 21 .AND. z <= 24) isearly3d = .TRUE.

   END FUNCTION early3d

END MODULE xtb_matrices

