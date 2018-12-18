!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Overlap and Hamiltonian matrices in xTB
!>        Reference: Stefan Grimme, Christoph Bannwarth, Philip Shushkov
!>                   JCTC 13, 1989-2009, (2017)
!>                   DOI: 10.1021/acs.jctc.7b00118
!> \author JGH
! **************************************************************************************************
MODULE xtb_matrices
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE block_p_types,                   ONLY: block_p_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              xtb_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_type, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_multiply, dbcsr_p_type, dbcsr_type, dbcsr_type_antisymmetric, &
        dbcsr_type_symmetric
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_sum
   USE mulliken,                        ONLY: ao_charges
   USE particle_types,                  ONLY: particle_type
   USE qs_dispersion_pairpot,           ONLY: d3_cnumber,&
                                              dcnum_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
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
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_coulomb,                     ONLY: build_xtb_coulomb
   USE xtb_parameters,                  ONLY: xtb_set_kab
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   ! screening for gamma function
   REAL(dp), PARAMETER                    :: tol_gamma = 1.e-4_dp
   ! small real number
   REAL(dp), PARAMETER                    :: rtiny = 1.e-10_dp
   !
   TYPE neighbor_atoms_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: coord
      REAL(KIND=dp), DIMENSION(:), POINTER       :: rab
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
         iw, j, jatom, jkind, katom, kkind, la, lb, lmaxa, lmaxb, n1, n2, na, natom, natorb_a, &
         natorb_b, nb, nderivatives, nimg, nkind, nmat, nsa, nsb, nshell, za, zat, zb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, atomnumber, kind_of
      INTEGER, DIMENSION(25)                             :: laoa, laob, lval, naoa, naob
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, diagblock, floating_a, found, &
                                                            ghost_a, h2scase, omit_headers, &
                                                            use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: floating, ghost
      REAL(KIND=dp) :: alphaa, alphab, derepij, dfij, dfp, dhij, dr, drk, drx, ena, enb, erep, &
         erepij, etaa, etab, exb, f0, f1, fen, fhua, fhub, fhud, fij, foab, hij, k2sh, kab, kcnd, &
         kcnp, kcns, kd, ken, kf, kg, kia, kjb, kp, ks, ksp, kx2, kxr, rcova, rcovab, rcovb, rcut, &
         rcuta, rcutb, rrab, xgammaa, xgammab, zneffa, zneffb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cnumbers
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dfblock, dgamma, huckel, kabset
      REAL(KIND=dp), DIMENSION(0:3)                      :: kcnl, kl
      REAL(KIND=dp), DIMENSION(3)                        :: fdik, force_ab, force_rr, rij, rik
      REAL(KIND=dp), DIMENSION(5)                        :: dpia, dpib, hena, henb, kappaa, kappab, &
                                                            kpolya, kpolyb, pia, pib
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fblock, gblock, pblock, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(block_p_type), DIMENSION(2:4)                 :: dgblocks, dsblocks
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, matrix_s, matrix_w
      TYPE(dcnum_type), ALLOCATABLE, DIMENSION(:)        :: dcnum
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_atoms_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: neighbor_atoms
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
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

      NULLIFY (matrix_h, matrix_s, matrix_p, matrix_w, gamma_matrix, atomic_kind_set, &
               qs_kind_set, sab_orb, ks_env)

      CALL get_qs_env(qs_env=qs_env, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      matrix_s_kp=matrix_s, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      ks_env=ks_env)

      xtb_control => dft_control%qs_control%xtb_control
      nimg = dft_control%nimages
      ! Allocate the overlap and Hamiltonian matrix
      CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb)
      nderivatives = 0
      IF (calculate_forces) nderivatives = 1

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

      ! gamma matrix allocation
      IF (calculate_forces) THEN
         nmat = 4
      ELSE
         nmat = 1
      END IF
      CALL get_qs_env(qs_env, gamma_matrix=gamma_matrix)
      CALL setup_gamma(qs_env, nmat, gamma_matrix, sab_orb)

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

      ! overlap matrix
      IF (calculate_forces) THEN
         CALL build_overlap_matrix(ks_env, nderivative=nderivatives, matrixkp_s=matrix_s, &
                                   matrix_name="xTB OVERLAP MATRIX", &
                                   basis_type_a="ORB", &
                                   basis_type_b="ORB", &
                                   sab_nl=sab_orb, calculate_forces=.TRUE., &
                                   matrixkp_p=matrix_w)
      ELSE
         CALL build_overlap_matrix(ks_env, nderivative=nderivatives, matrixkp_s=matrix_s, &
                                   matrix_name="xTB OVERLAP MATRIX", &
                                   basis_type_a="ORB", &
                                   basis_type_b="ORB", &
                                   sab_nl=sab_orb)
      END IF
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

      nkind = SIZE(atomic_kind_set)
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

      ! Calculate Huckel parameters
      ! Eq 12
      ! huckel(nshell,natom)
      kcnl(0) = kcns
      kcnl(1) = kcnp
      kcnl(2) = kcnd
      kcnl(3) = 0.0_dp
      kl(0) = ks
      kl(1) = kp
      kl(2) = kd
      kl(3) = 0.0_dp
      ALLOCATE (huckel(5, natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, nshell=nshell, lval=lval, hen=hena)
         huckel(:, iatom) = 0.0_dp
         DO i = 1, nshell
            huckel(i, iatom) = hena(i)*(1._dp+kcnl(lval(i))*cnumbers(iatom))
         END DO
      END DO

      ! Calculate KAB parameters
      ! precalculate the KAB for all kind pairs (Table 2)
      ALLOCATE (kabset(nkind, nkind))
      kabset = 1.0_dp
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), z=za)
         DO jkind = ikind+1, nkind
            CALL get_atomic_kind(atomic_kind_set(jkind), z=zb)
            kabset(ikind, jkind) = xtb_set_kab(za, zb)
            kabset(ikind, jkind) = kabset(jkind, ikind)
         END DO
      ENDDO

      ! list of neighbor atoms for XB term
      ALLOCATE (neighbor_atoms(nkind))
      DO ikind = 1, nkind
         NULLIFY (neighbor_atoms(ikind)%coord)
         NULLIFY (neighbor_atoms(ikind)%rab)
         CALL get_atomic_kind(atomic_kind_set(ikind), z=zat, natom=na)
         IF (zat == 17 .OR. zat == 35 .OR. zat == 53 .OR. zat == 85) THEN
            ALLOCATE (neighbor_atoms(ikind)%coord(3, na))
            neighbor_atoms(ikind)%coord(1:3, 1:na) = 0.0_dp
            ALLOCATE (neighbor_atoms(ikind)%rab(na))
            neighbor_atoms(ikind)%rab(1:na) = HUGE(0.0_dp)
         END IF
      END DO

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
         IF (ASSOCIATED(neighbor_atoms(ikind)%rab)) THEN
            atom_a = atom_of_kind(iatom)
            IF (dr < neighbor_atoms(ikind)%rab(atom_a)) THEN
               neighbor_atoms(ikind)%rab(atom_a) = dr
               neighbor_atoms(ikind)%coord(1:3, atom_a) = rij(1:3)
            END IF
         END IF
         IF (ASSOCIATED(neighbor_atoms(jkind)%rab)) THEN
            atom_b = atom_of_kind(jatom)
            IF (dr < neighbor_atoms(jkind)%rab(atom_b)) THEN
               neighbor_atoms(jkind)%rab(atom_b) = dr
               neighbor_atoms(jkind)%coord(1:3, atom_b) = -rij(1:3)
            END IF
         END IF

         ! atomic parameters
         CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa, lao=laoa, rcov=rcova, eta=etaa, xgamma=xgammaa, &
                                 lmax=lmaxa, nshell=nsa, alpha=alphaa, zneff=zneffa, kpoly=kpolya, &
                                 kappa=kappaa, hen=hena, electronegativity=ena)
         CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob, lao=laob, rcov=rcovb, eta=etab, xgamma=xgammab, &
                                 lmax=lmaxb, nshell=nsb, alpha=alphab, zneff=zneffb, kpoly=kpolyb, &
                                 kappa=kappab, hen=henb, electronegativity=enb)

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cell(1), cell(2), cell(3))
            CPASSERT(ic > 0)
         END IF

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)
         NULLIFY (sblock, fblock, gblock)
         CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=matrix_h(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=fblock, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=gamma_matrix(1)%matrix, &
                                row=irow, col=icol, BLOCK=gblock, found=found)
         CPASSERT(found)

         IF (calculate_forces) THEN
            NULLIFY (pblock)
            CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                   row=irow, col=icol, block=pblock, found=found)
            CPASSERT(ASSOCIATED(pblock))
            DO i = 2, 4
               NULLIFY (dsblocks(i)%block)
               CALL dbcsr_get_block_p(matrix=matrix_s(i, ic)%matrix, &
                                      row=irow, col=icol, BLOCK=dsblocks(i)%block, found=found)
               CPASSERT(found)
               NULLIFY (dgblocks(i)%block)
               CALL dbcsr_get_block_p(matrix=gamma_matrix(i)%matrix, &
                                      row=irow, col=icol, BLOCK=dgblocks(i)%block, found=found)
               CPASSERT(found)
            END DO
         END IF

         ! Calculate Pi = Pia * Pib (Eq. 11)
         rcovab = rcova+rcovb
         rrab = SQRT(dr/rcovab)
         DO i = 1, nsa
            pia(i) = 1._dp+kpolya(i)*rrab
         END DO
         DO i = 1, nsb
            pib(i) = 1._dp+kpolyb(i)*rrab
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
         ! get KAB
         kab = kabset(ikind, jkind)
         ! get Fen = (1+ken*deltaEN^2)
         fen = 1.0_dp+ken*(ena-enb)**2
         !
         DO j = 1, natorb_b
            lb = laob(j)
            nb = naob(j)
            DO i = 1, natorb_a
               la = laoa(i)
               na = naoa(i)
               IF (diagblock .AND. i == j) THEN
                  fblock(i, i) = fblock(i, i)+huckel(na, iatom)
               ELSE
                  kia = kl(la)
                  kjb = kl(lb)
                  h2scase = .FALSE.
                  IF (zb == 1 .AND. nb == 2) THEN
                     kjb = k2sh
                     h2scase = .TRUE.
                  END IF
                  IF (za == 1 .AND. na == 2) THEN
                     kia = k2sh
                     h2scase = .TRUE.
                  END IF
                  hij = 0.5_dp*(huckel(na, iatom)+huckel(nb, jatom))*pia(na)*pib(nb)
                  IF (h2scase) THEN
                     fij = 0.5_dp*(kia+kjb)*hij
                  ELSE
                     IF ((la == 0 .AND. lb == 1) .OR. (la == 1 .AND. lb == 0)) THEN
                        fij = ksp*hij*kab*fen
                     ELSE
                        fij = 0.5_dp*(kia+kjb)*hij*kab*fen
                     END IF
                  END IF
                  IF (iatom <= jatom) THEN
                     fblock(i, j) = fblock(i, j)+fij*sblock(i, j)
                  ELSE
                     fblock(j, i) = fblock(j, i)+fij*sblock(j, i)
                  END IF
               END IF
            END DO
         END DO
         IF (calculate_forces) THEN
            f0 = 1.0_dp
            IF (irow == iatom) f0 = -1.0_dp
            ! Derivative wrt coordination number
            fhua = 0.0_dp
            fhub = 0.0_dp
            fhud = 0.0_dp
            DO j = 1, natorb_b
               lb = laob(j)
               nb = naob(j)
               DO i = 1, natorb_a
                  la = laoa(i)
                  na = naoa(i)
                  IF (diagblock .AND. i == j) THEN
                     fhud = fhud+pblock(i, i)*kcnl(la)*hena(na)
                  ELSE
                     kia = kl(la)
                     kjb = kl(lb)
                     h2scase = .FALSE.
                     IF (zb == 1 .AND. nb == 2) THEN
                        kjb = k2sh
                        h2scase = .TRUE.
                     END IF
                     IF (za == 1 .AND. na == 2) THEN
                        kia = k2sh
                        h2scase = .TRUE.
                     END IF
                     hij = 0.5_dp*pia(na)*pib(nb)
                     IF (h2scase) THEN
                        fij = 0.5_dp*(kia+kjb)*hij
                     ELSE
                        IF ((la == 0 .AND. lb == 1) .OR. (la == 1 .AND. lb == 0)) THEN
                           fij = ksp*hij*kab*fen
                        ELSE
                           fij = 0.5_dp*(kia+kjb)*hij*kab*fen
                        END IF
                     END IF
                     IF (iatom <= jatom) THEN
                        fhua = fhua+fij*sblock(i, j)*pblock(i, j)*kcnl(la)*hena(na)
                        fhub = fhub+fij*sblock(i, j)*pblock(i, j)*kcnl(lb)*henb(nb)
                     ELSE
                        fhua = fhua+fij*sblock(j, i)*pblock(j, i)*kcnl(la)*hena(na)
                        fhub = fhub+fij*sblock(j, i)*pblock(j, i)*kcnl(lb)*henb(nb)
                     END IF
                  END IF
               END DO
            END DO
            IF (.NOT. diagblock) THEN
               fhua = 2.0_dp*fhua
               fhub = 2.0_dp*fhub
            END IF
            ! iatom
            atom_a = atom_of_kind(iatom)
            DO i = 1, dcnum(iatom)%neighbors
               katom = dcnum(iatom)%nlist(i)
               kkind = kind_of(katom)
               atom_c = atom_of_kind(katom)
               rik = dcnum(iatom)%rik(:, i)
               drk = SQRT(SUM(rik(:)**2))
               fdik(:) = fhua*dcnum(iatom)%dvals(i)*rik(:)/drk
               force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a)-fdik(:)
               force(kkind)%all_potential(:, atom_c) = force(kkind)%all_potential(:, atom_c)+fdik(:)
               fdik(:) = fhud*dcnum(iatom)%dvals(i)*rik(:)/drk
               force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a)-fdik(:)
               force(kkind)%all_potential(:, atom_c) = force(kkind)%all_potential(:, atom_c)+fdik(:)
               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
               END IF
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
                  CALL virial_pair_force(atprop%atstress(:, :, katom), -0.5_dp, fdik, rik)
               END IF
            END DO
            ! iatom
            atom_b = atom_of_kind(jatom)
            DO i = 1, dcnum(jatom)%neighbors
               katom = dcnum(jatom)%nlist(i)
               kkind = kind_of(katom)
               atom_c = atom_of_kind(katom)
               rik = dcnum(jatom)%rik(:, i)
               drk = SQRT(SUM(rik(:)**2))
               fdik(:) = fhub*dcnum(jatom)%dvals(i)*rik(:)/drk
               force(jkind)%all_potential(:, atom_b) = force(jkind)%all_potential(:, atom_b)-fdik(:)
               force(kkind)%all_potential(:, atom_c) = force(kkind)%all_potential(:, atom_c)+fdik(:)
               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, -1._dp, fdik, rik)
               END IF
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, fdik, rik)
                  CALL virial_pair_force(atprop%atstress(:, :, katom), -0.5_dp, fdik, rik)
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
                     kia = kl(la)
                     kjb = kl(lb)
                     h2scase = .FALSE.
                     IF (zb == 1 .AND. nb == 2) THEN
                        kjb = k2sh
                        h2scase = .TRUE.
                     END IF
                     IF (za == 1 .AND. na == 2) THEN
                        kia = k2sh
                        h2scase = .TRUE.
                     END IF
                     dhij = 0.5_dp*(huckel(na, iatom)+huckel(nb, jatom))*(dpia(na)*pib(nb)+pia(na)*dpib(nb))
                     IF (h2scase) THEN
                        dfij = 0.5_dp*(kia+kjb)*dhij
                     ELSE
                        IF ((la == 0 .AND. lb == 1) .OR. (la == 1 .AND. lb == 0)) THEN
                           dfij = ksp*dhij*kab*fen
                        ELSE
                           dfij = 0.5_dp*(kia+kjb)*dhij*kab*fen
                        END IF
                     END IF
                     IF (iatom <= jatom) THEN
                        dfblock(i, j) = dfblock(i, j)+dfij*sblock(i, j)
                     ELSE
                        dfblock(j, i) = dfblock(j, i)+dfij*sblock(j, i)
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
                        kia = kl(la)
                        kjb = kl(lb)
                        h2scase = .FALSE.
                        IF (zb == 1 .AND. nb == 2) THEN
                           kjb = k2sh
                           h2scase = .TRUE.
                        END IF
                        IF (za == 1 .AND. na == 2) THEN
                           kia = k2sh
                           h2scase = .TRUE.
                        END IF
                        hij = 0.5_dp*(huckel(na, iatom)+huckel(nb, jatom))*pia(na)*pib(nb)
                        IF (h2scase) THEN
                           fij = 0.5_dp*(kia+kjb)*hij
                        ELSE
                           IF ((la == 0 .AND. lb == 1) .OR. (la == 1 .AND. lb == 0)) THEN
                              fij = ksp*hij*kab*fen
                           ELSE
                              fij = 0.5_dp*(kia+kjb)*hij*kab*fen
                           END IF
                        END IF
                        IF (iatom <= jatom) THEN
                           foab = foab+2.0_dp*fij*dsblocks(ir+1)%block(i, j)*pblock(i, j)
                        ELSE
                           foab = foab+2.0_dp*fij*dsblocks(ir+1)%block(j, i)*pblock(j, i)
                        END IF
                     END DO
                  END DO
                  force_ab(ir) = foab
               END DO
               DEALLOCATE (dfblock)
            END IF
            IF (use_virial) THEN
               CALL virial_pair_force(virial%pv_virial, -f0, force_ab, rij)
               IF (atprop%stress) THEN
                  f1 = 0.5_dp*f0
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -f1, force_ab, rij)
                  CALL virial_pair_force(atprop%atstress(:, :, jatom), -f1, force_ab, rij)
               END IF
            END IF
         END IF

         IF (calculate_forces) THEN
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)
            IF (irow == iatom) force_ab = -force_ab
            force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a)-force_ab(:)
            force(jkind)%all_potential(:, atom_b) = force(jkind)%all_potential(:, atom_b)+force_ab(:)
         END IF

         ! gamma matrix
         CALL get_xtb_atom_param(xtb_atom_a, rcut=rcuta)
         CALL get_xtb_atom_param(xtb_atom_b, rcut=rcutb)
         rcut = rcuta+rcutb
         IF (irow == iatom) THEN
            CALL gamma_rab_sr(gblock, dr, lmaxa+1, kappaa, etaa, lmaxb+1, kappab, etab, kg, rcut)
         ELSE
            CALL gamma_rab_sr(gblock, dr, lmaxb+1, kappab, etab, lmaxa+1, kappaa, etaa, kg, rcut)
         END IF
         IF (calculate_forces .AND. (iatom /= jatom .OR. dr > 0.001_dp)) THEN
            n1 = SIZE(gblock, 1)
            n2 = SIZE(gblock, 2)
            ALLOCATE (dgamma(n1, n2))
            IF (irow == iatom) THEN
               CALL dgamma_rab_sr(dgamma, dr, lmaxa+1, kappaa, etaa, lmaxb+1, kappab, etab, kg, rcut)
            ELSE
               CALL dgamma_rab_sr(dgamma, dr, lmaxb+1, kappab, etab, lmaxa+1, kappaa, etaa, kg, rcut)
            END IF
            DO i = 1, 3
               CPASSERT(ASSOCIATED(dgblocks(i+1)%block))
               IF (irow == iatom) THEN
                  dgblocks(i+1)%block(1:n1, 1:n2) = dgblocks(i+1)%block(1:n1, 1:n2)+dgamma(1:n1, 1:n2)*rij(i)/dr
               ELSE
                  dgblocks(i+1)%block(1:n1, 1:n2) = dgblocks(i+1)%block(1:n1, 1:n2)-dgamma(1:n1, 1:n2)*rij(i)/dr
               END IF
            END DO
            DEALLOCATE (dgamma)
         END IF

         ! repulsive potential
         IF (dr > 0.001_dp) THEN
            erepij = zneffa*zneffb/dr*EXP(-SQRT(alphaa*alphab)*dr**kf)
            erep = erep+erepij
            IF (atprop%energy) THEN
               atprop%atecc(iatom) = atprop%atecc(iatom)+0.5_dp*erepij
               atprop%atecc(jatom) = atprop%atecc(jatom)+0.5_dp*erepij
            END IF
            IF (calculate_forces .AND. (iatom /= jatom .OR. dr > 0.001_dp)) THEN
               derepij = -(1.0_dp/dr+SQRT(alphaa*alphab)*kf*dr**(kf-1.0_dp))*erepij
               force_rr(1) = derepij*rij(1)/dr
               force_rr(2) = derepij*rij(2)/dr
               force_rr(3) = derepij*rij(3)/dr
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               force(ikind)%repulsive(:, atom_a) = &
                  force(ikind)%repulsive(:, atom_a)-force_rr(:)
               force(jkind)%repulsive(:, atom_b) = &
                  force(jkind)%repulsive(:, atom_b)+force_rr(:)
               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, -1._dp, force_rr, rij)
                  IF (atprop%stress) THEN
                     CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp, force_rr, rij)
                     CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp, force_rr, rij)
                  END IF
               END IF
            END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(gamma_matrix)
         CALL dbcsr_finalize(gamma_matrix(i)%matrix)
      ENDDO
      CALL set_ks_env(ks_env, gamma_matrix=gamma_matrix)
      DO i = 1, SIZE(matrix_h, 1)
         DO img = 1, nimg
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
         END DO
      ENDDO

!qxtb
      exb = 0.0_dp
!   calculate halogen correction term here
!     CALL xb_interaction(exb,neighbor_atoms,sab_xb,kx,kx2,kxr)
      erep = erep+exb
!qxtb
      ! set repulsive energy
      CALL mp_sum(erep, para_env%group)
      energy%repulsive = erep

      ! deallocate coordination numbers
      DEALLOCATE (cnumbers)
      IF (calculate_forces) THEN
         DO iatom = 1, natom
            DEALLOCATE (dcnum(iatom)%nlist, dcnum(iatom)%dvals, dcnum(iatom)%rik)
         END DO
         DEALLOCATE (dcnum)
      ELSE
         DEALLOCATE (dcnum)
      END IF

      ! deallocate Huckel parameters
      DEALLOCATE (huckel)
      ! deallocate KAB parameters
      DEALLOCATE (kabset)

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

         CALL cp_print_key_finished_output(iw, logger, qs_env%input, &
                                           "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN")
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

         CALL cp_print_key_finished_output(iw, logger, qs_env%input, &
                                           "DFT%PRINT%AO_MATRICES/OVERLAP")
      END IF

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

      DO ikind = 1, nkind
         IF (ASSOCIATED(neighbor_atoms(ikind)%coord)) THEN
            DEALLOCATE (neighbor_atoms(ikind)%coord)
         END IF
         IF (ASSOCIATED(neighbor_atoms(ikind)%rab)) THEN
            DEALLOCATE (neighbor_atoms(ikind)%rab)
         END IF
      END DO
      DEALLOCATE (neighbor_atoms)

      CALL timestop(handle)

   END SUBROUTINE build_xtb_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE build_xtb_ks_matrix(qs_env, calculate_forces, just_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'build_xtb_ks_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, iatom, ikind, img, &
                                                            iounit, is, ispin, natom, natorb, &
                                                            nkind, ns, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: aocharge, mcharge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: charges
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix, matrix_h, matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: scf_section
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)
      NULLIFY (dft_control, logger, scf_section, matrix_p, particle_set, ks_env, &
               ks_matrix, rho, energy)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(qs_env))

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

      energy%hartree = 0.0_dp

      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")
      nspins = dft_control%nspins
      CPASSERT(ASSOCIATED(matrix_h))
      CPASSERT(ASSOCIATED(rho))
      CPASSERT(SIZE(ks_matrix) > 0)

      DO ispin = 1, nspins
         DO img = 1, SIZE(ks_matrix, 2)
            ! copy the core matrix into the fock matrix
            CALL dbcsr_copy(ks_matrix(ispin, img)%matrix, matrix_h(1, img)%matrix)
         END DO
      END DO

      ! Mulliken charges
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, matrix_s_kp=matrix_s)
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      natom = SIZE(particle_set)
      ALLOCATE (mcharge(natom), charges(natom, 5))
      charges = 0.0_dp
      nkind = SIZE(atomic_kind_set)
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, occupation=occ)
         ALLOCATE (aocharge(natorb))
         DO iatom = 1, natom
            atom_a = atomic_kind_set(ikind)%atom_list(iatom)
            CALL ao_charges(matrix_p, matrix_s, aocharge, atom_a, para_env)
            charges(atom_a, :) = REAL(occ(:), KIND=dp)
            DO is = 1, natorb
               ns = lao(is)+1
               charges(atom_a, ns) = charges(atom_a, ns)-aocharge(is)
            END DO
            mcharge(atom_a) = SUM(charges(atom_a, :))
         END DO
         DEALLOCATE (aocharge)
      END DO

      CALL build_xtb_coulomb(qs_env, ks_matrix, rho, charges, mcharge, energy, &
                             calculate_forces, just_energy)

      DEALLOCATE (mcharge, charges)

      IF (qs_env%qmmm) THEN
         CPABORT("QMMM with xTB not implemented")
      END IF

      energy%total = energy%core+energy%hartree+ &
                     energy%repulsive+energy%dispersion+energy%dftb3

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Charge fluctuation energy:                     ", energy%hartree, &
            "London dispersion energy:                      ", energy%dispersion
         WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "DFTB3 3rd Order Energy Correction              ", energy%dftb3
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
!> \brief  Computes the short-range gamma parameter from
!>         Nataga-Mishimoto-Ohno-Klopman formula for xTB
!>         WARNING: The xTB function (gamma - 1/r) has still an l-dependent longrange (1/r^3)
!>                  behaviour. We use a cutoff function to smoothly remove this part.
!>                  However, this will change energies and effect final results.
!>
!> \param gmat ...
!> \param rab ...
!> \param nla ...
!> \param kappaa ...
!> \param etaa ...
!> \param nlb ...
!> \param kappab ...
!> \param etab ...
!> \param kg ...
!> \param rcut ...
!> \par History
!>      10.2018 JGH
!> \version 1.1
! **************************************************************************************************
   SUBROUTINE gamma_rab_sr(gmat, rab, nla, kappaa, etaa, nlb, kappab, etab, kg, rcut)
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: gmat
      REAL(dp), INTENT(IN)                               :: rab
      INTEGER, INTENT(IN)                                :: nla
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kappaa
      REAL(dp), INTENT(IN)                               :: etaa
      INTEGER, INTENT(IN)                                :: nlb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kappab
      REAL(dp), INTENT(IN)                               :: etab, kg, rcut

      CHARACTER(len=*), PARAMETER :: routineN = 'gamma_rab_sr', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: fcut, r, rk, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eta

      ALLOCATE (eta(nla, nlb))

      DO j = 1, nlb
         DO i = 1, nla
            eta(i, j) = 1._dp/(etaa*(1._dp+kappaa(i)))+1._dp/(etab*(1._dp+kappab(j)))
            eta(i, j) = 2._dp/eta(i, j)
         END DO
      END DO

      IF (rab < 1.e-6) THEN
         ! on site terms
         gmat(:, :) = gmat(:, :)+eta(:, :)
      ELSEIF (rab > rcut) THEN
         ! do nothing
      ELSE
         rk = rab**kg
         eta = eta**(-kg)
         IF (rab < rcut-rsmooth) THEN
            fcut = 1.0_dp
         ELSE
            r = rab-(rcut-rsmooth)
            x = r/rsmooth
            fcut = -10._dp*x**5+15._dp*x**4-10._dp*x**3+1._dp
         END IF
         gmat(:, :) = gmat(:, :)+fcut*(1._dp/(rk+eta(:, :)))**(1._dp/kg)
         gmat(:, :) = gmat(:, :)-fcut*1.0_dp/rab
      END IF

      DEALLOCATE (eta)

   END SUBROUTINE gamma_rab_sr

! **************************************************************************************************
!> \brief  Computes the derivative of the short-range gamma parameter from
!>         Nataga-Mishimoto-Ohno-Klopman formula for xTB
!>         WARNING: The xTB function (gamma - 1/r) has still an l-dependent longrange (1/r^3)
!>                  behaviour. We use a cutoff function to smoothly remove this part.
!>                  However, this will change energies and effect final results.
!>
!> \param dgmat ...
!> \param rab ...
!> \param nla ...
!> \param kappaa ...
!> \param etaa ...
!> \param nlb ...
!> \param kappab ...
!> \param etab ...
!> \param kg ...
!> \param rcut ...
!> \par History
!>      10.2018 JGH
!> \version 1.1
! **************************************************************************************************
   SUBROUTINE dgamma_rab_sr(dgmat, rab, nla, kappaa, etaa, nlb, kappab, etab, kg, rcut)
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: dgmat
      REAL(dp), INTENT(IN)                               :: rab
      INTEGER, INTENT(IN)                                :: nla
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kappaa
      REAL(dp), INTENT(IN)                               :: etaa
      INTEGER, INTENT(IN)                                :: nlb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kappab
      REAL(dp), INTENT(IN)                               :: etab, kg, rcut

      CHARACTER(len=*), PARAMETER :: routineN = 'dgamma_rab_sr', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: dfcut, fcut, r, rk, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eta

      ALLOCATE (eta(nla, nlb))

      DO j = 1, nlb
         DO i = 1, nla
            eta(i, j) = 1._dp/(etaa*(1._dp+kappaa(i)))+1._dp/(etab*(1._dp+kappab(j)))
            eta(i, j) = 2._dp/eta(i, j)
         END DO
      END DO

      IF (rab < 1.e-6) THEN
         ! on site terms
         dgmat(:, :) = 0.0_dp
      ELSEIF (rab > rcut) THEN
         dgmat(:, :) = 0.0_dp
      ELSE
         eta = eta**(-kg)
         rk = rab**kg
         IF (rab < rcut-rsmooth) THEN
            fcut = 1.0_dp
            dfcut = 0.0_dp
         ELSE
            r = rab-(rcut-rsmooth)
            x = r/rsmooth
            fcut = -10._dp*x**5+15._dp*x**4-10._dp*x**3+1._dp
            dfcut = -50._dp*x**4+60._dp*x**3-30._dp*x**2
            dfcut = dfcut/rsmooth
         END IF
         dgmat(:, :) = dfcut*(1._dp/(rk+eta(:, :)))**(1._dp/kg)
         dgmat(:, :) = dgmat(:, :)-dfcut/rab+fcut/rab**2
         dgmat(:, :) = dgmat(:, :)-fcut/(rk+eta(:, :))*(1._dp/(rk+eta(:, :)))**(1._dp/kg)*rk/rab
      END IF

      DEALLOCATE (eta)

   END SUBROUTINE dgamma_rab_sr

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nmat ...
!> \param gammat ...
!> \param sab_nl ...
! **************************************************************************************************
   SUBROUTINE setup_gamma(qs_env, nmat, gammat, sab_nl)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nmat
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gammat
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_gamma', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iatom, ikind, natom, &
                                                            neighbor_list_id, ns
      INTEGER, DIMENSION(:), POINTER                     :: kind_of, row_blk_sizes
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(xtb_atom_type), POINTER                       :: xtb_parameter

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      neighbor_list_id=neighbor_list_id, &
                      dbcsr_dist=dbcsr_dist)

      natom = SIZE(particle_set)
      ALLOCATE (row_blk_sizes(natom))
      ALLOCATE (kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      DO iatom = 1, natom
         ikind = kind_of(iatom)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, xtb_parameter=xtb_parameter)
         CALL get_xtb_atom_param(xtb_parameter, lmax=ns)
         row_blk_sizes(iatom) = ns+1
      ENDDO

      CALL dbcsr_allocate_matrix_set(gammat, nmat)
      ALLOCATE (gammat(1)%matrix)
      CALL dbcsr_create(matrix=gammat(1)%matrix, &
                        name="GAMMA MATRIX", &
                        dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                        row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                        nze=0, mutable_work=.TRUE.)

      DO i = 2, nmat
         ALLOCATE (gammat(i)%matrix)
         CALL dbcsr_create(matrix=gammat(i)%matrix, &
                           name="DERIVATIVE GAMMA MATRIX", &
                           dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0, mutable_work=.TRUE.)
      END DO

      DEALLOCATE (row_blk_sizes, kind_of)

      ! setup the matrices using the neighbor list
      DO i = 1, nmat
         CALL cp_dbcsr_alloc_block_from_nbl(gammat(i)%matrix, sab_nl)
      END DO

   END SUBROUTINE setup_gamma

END MODULE xtb_matrices

