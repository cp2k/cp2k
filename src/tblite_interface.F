!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief interface to tblite
!> \author JVP
!> \history creation 09.2024
! **************************************************************************************************

MODULE tblite_interface

#if defined(__TBLITE)
   USE mctc_env, ONLY: error_type
   USE mctc_io, ONLY: structure_type, new
   USE mctc_io_symbols, ONLY: symbol_to_number
   USE tblite_basis_type, ONLY: get_cutoff
   USE tblite_container, ONLY: container_cache
   USE tblite_integral_multipole, ONLY: multipole_cgto, multipole_grad_cgto, maxl, msao
   USE tblite_scf_info, ONLY: scf_info, atom_resolved, shell_resolved, &
                              orbital_resolved, not_used
   USE tblite_scf_potential, ONLY: potential_type, new_potential
   USE tblite_wavefunction_type, ONLY: wavefunction_type, new_wavefunction
   USE tblite_xtb_calculator, ONLY: xtb_calculator, new_xtb_calculator
   USE tblite_xtb_gfn1, ONLY: new_gfn1_calculator
   USE tblite_xtb_gfn2, ONLY: new_gfn2_calculator
   USE tblite_xtb_h0, ONLY: get_selfenergy, get_hamiltonian, get_occupation, &
                            get_hamiltonian_gradient, tb_hamiltonian
   USE tblite_xtb_ipea1, ONLY: new_ipea1_calculator
#endif
   USE ai_contraction, ONLY: block_add, &
                             contraction
   USE ai_overlap, ONLY: overlap_ab
   USE atomic_kind_types, ONLY: atomic_kind_type, get_atomic_kind, get_atomic_kind_set
   USE atprop_types, ONLY: atprop_type
   USE basis_set_types, ONLY: gto_basis_set_type, gto_basis_set_p_type, &
                              & allocate_gto_basis_set, write_gto_basis_set, process_gto_basis
   USE cell_types, ONLY: cell_type, get_cell
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_control_types, ONLY: dft_control_type
   USE cp_dbcsr_api, ONLY: dbcsr_type, dbcsr_p_type, dbcsr_create, dbcsr_add, &
                           dbcsr_get_block_p, dbcsr_finalize, &
                           dbcsr_iterator_type, dbcsr_iterator_blocks_left, &
                           dbcsr_iterator_start, dbcsr_iterator_stop, &
                           dbcsr_iterator_next_block
   USE cp_dbcsr_contrib, ONLY: dbcsr_print
   USE cp_dbcsr_cp2k_link, ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations, ONLY: dbcsr_allocate_matrix_set
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_type, cp_logger_get_default_io_unit
   USE cp_output_handling, ONLY: cp_print_key_should_output, &
                                 cp_print_key_unit_nr, cp_print_key_finished_output
   USE input_constants, ONLY: gfn1xtb, gfn2xtb, ipea1xtb
   USE input_section_types, ONLY: section_vals_val_get
   USE kinds, ONLY: dp, default_string_length
   USE kpoint_types, ONLY: get_kpoint_info, kpoint_type
   USE memory_utilities, ONLY: reallocate
   USE message_passing, ONLY: mp_para_env_type
   USE mulliken, ONLY: ao_charges
   USE orbital_pointers, ONLY: ncoset
   USE particle_types, ONLY: particle_type
   USE qs_charge_mixing, ONLY: charge_mixing
   USE qs_condnum, ONLY: overlap_condnum
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, qs_environment_type
   USE qs_force_types, ONLY: qs_force_type
   USE qs_integral_utils, ONLY: basis_set_list_setup, get_memory_usage
   USE qs_kind_types, ONLY: get_qs_kind, qs_kind_type, get_qs_kind_set
   USE qs_ks_types, ONLY: qs_ks_env_type, set_ks_env
   USE qs_neighbor_list_types, ONLY: neighbor_list_iterator_create, neighbor_list_iterate, &
                                     get_iterator_info, neighbor_list_set_p_type, &
                                     neighbor_list_iterator_p_type, neighbor_list_iterator_release
   USE qs_overlap, ONLY: create_sab_matrix
   USE qs_rho_types, ONLY: qs_rho_get, qs_rho_type
   USE qs_scf_types, ONLY: qs_scf_env_type
   USE input_section_types, ONLY: section_vals_get_subs_vals, section_vals_type
   USE string_utilities, ONLY: integer_to_string
   USE tblite_types, ONLY: tblite_type, deallocate_tblite_type, allocate_tblite_type
   USE virial_types, ONLY: virial_type
   USE xtb_types, ONLY: get_xtb_atom_param, xtb_atom_type
   USE xtb_types, ONLY: xtb_atom_type

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_interface'

   INTEGER, PARAMETER                                 :: dip_n = 3
   INTEGER, PARAMETER                                 :: quad_n = 6

   PUBLIC :: tb_set_calculator, tb_init_geometry, tb_init_wf
   PUBLIC :: tb_get_basis, build_tblite_matrices
   PUBLIC :: tb_get_energy, tb_update_charges, tb_ham_add_coulomb
   PUBLIC :: tb_get_multipole
   PUBLIC :: tb_derive_dH_diag, tb_derive_dH_off

CONTAINS

! **************************************************************************************************
!> \brief intialize geometry objects ...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_geometry(qs_env, tb)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_init_geometry'

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER                                            :: iatom, natom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle, ikind
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL, DIMENSION(3)                              :: lperiod

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, qs_kind_set=qs_kind_set)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom))
      CALL allocate_tblite_type(tb)
      ALLOCATE (tb%el_num(natom))
      tb%el_num = -9
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=tb%el_num(iatom))
         IF (tb%el_num(iatom) < 1 .OR. tb%el_num(iatom) > 85) THEN
            CPABORT("only elements 1-85 are supported by tblite")
         END IF
      END DO

      !get information about cell / lattice
      CALL get_cell(cell=cell, periodic=periodic)
      lperiod(1) = periodic(1) == 1
      lperiod(2) = periodic(2) == 1
      lperiod(3) = periodic(3) == 1

      !prepare for the call to the dispersion function
      CALL new(tb%mol, tb%el_num, xyz, lattice=cell%hmat, periodic=lperiod)

      DEALLOCATE (xyz)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_geometry

! **************************************************************************************************
!> \brief updating coordinates...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_update_geometry(qs_env, tb)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_update_geometry'

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER                                            :: iatom, natom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom))
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
      END DO
      tb%mol%xyz(:, :) = xyz

      DEALLOCATE (xyz)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_update_geometry

! **************************************************************************************************
!> \brief initialize wavefunction ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_wf(tb)

      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE)

      INTEGER, PARAMETER                                 :: nSpin = 1 !number of spin channels

      TYPE(scf_info)                                     :: info

      info = tb%calc%variable_info()
      IF (info%charge > shell_resolved) CPABORT("tblite: no support for orbital resolved charge")
      IF (info%dipole > atom_resolved) CPABORT("tblite: no support for shell resolved dipole moment")
      IF (info%quadrupole > atom_resolved) &
         CPABORT("tblite: no support shell resolved quadrupole moment")

      CALL new_wavefunction(tb%wfn, tb%mol%nat, tb%calc%bas%nsh, tb%calc%bas%nao, nSpin, 0.0_dp)

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)

      !allocate quantities later required
      ALLOCATE (tb%e_hal(tb%mol%nat), tb%e_rep(tb%mol%nat), tb%e_disp(tb%mol%nat))
      ALLOCATE (tb%e_scd(tb%mol%nat), tb%e_es(tb%mol%nat))
      ALLOCATE (tb%selfenergy(tb%calc%bas%nsh))
      IF (ALLOCATED(tb%calc%ncoord)) ALLOCATE (tb%cn(tb%mol%nat))

#else
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_wf

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param typ ...
! **************************************************************************************************
   SUBROUTINE tb_set_calculator(tb, typ)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER                                            :: typ

#if defined(__TBLITE)

      TYPE(error_type), ALLOCATABLE                      :: error

      SELECT CASE (typ)
      CASE default
         CPABORT("Unknown xtb type")
      CASE (gfn1xtb)
         CALL new_gfn1_calculator(tb%calc, tb%mol, error)
      CASE (gfn2xtb)
         CALL new_gfn2_calculator(tb%calc, tb%mol, error)
      CASE (ipea1xtb)
         CALL new_ipea1_calculator(tb%calc, tb%mol, error)
      END SELECT

#else
      MARK_USED(tb)
      MARK_USED(typ)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_set_calculator

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_init_ham(qs_env, tb, para_env)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env

#if defined(__TBLITE)

      TYPE(container_cache)                              :: hcache, rcache

      tb%e_hal = 0.0_dp
      tb%e_rep = 0.0_dp
      tb%e_disp = 0.0_dp
      IF (ALLOCATED(tb%grad)) THEN
         tb%grad = 0.0_dp
         CALL tb_zero_force(qs_env)
      END IF
      tb%sigma = 0.0_dp

      IF (ALLOCATED(tb%calc%halogen)) THEN
         CALL tb%calc%halogen%update(tb%mol, hcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal, &
            & tb%grad, tb%sigma)
            CALL tb_grad2force(qs_env, tb, para_env, 0)
         ELSE
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%repulsion)) THEN
         CALL tb%calc%repulsion%update(tb%mol, rcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep, &
            & tb%grad, tb%sigma)
            CALL tb_grad2force(qs_env, tb, para_env, 1)
         ELSE
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%dispersion)) THEN
         CALL tb%calc%dispersion%update(tb%mol, tb%dcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp, &
            & tb%grad, tb%sigma)
            CALL tb_grad2force(qs_env, tb, para_env, 2)
         ELSE
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp)
         END IF
      END IF

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)
      IF (ALLOCATED(tb%calc%coulomb)) THEN
         CALL tb%calc%coulomb%update(tb%mol, tb%cache)
      END IF

      IF (ALLOCATED(tb%grad)) THEN
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            CALL tb%calc%ncoord%get_cn(tb%mol, tb%cn, tb%dcndr, tb%dcndL)
         END IF
         CALL get_selfenergy(tb%calc%h0, tb%mol%id, tb%calc%bas%ish_at, &
         &  tb%calc%bas%nsh_id, cn=tb%cn, selfenergy=tb%selfenergy, dsedcn=tb%dsedcn)
      ELSE
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            CALL tb%calc%ncoord%get_cn(tb%mol, tb%cn)
         END IF
         CALL get_selfenergy(tb%calc%h0, tb%mol%id, tb%calc%bas%ish_at, &
         &  tb%calc%bas%nsh_id, cn=tb%cn, selfenergy=tb%selfenergy, dsedcn=tb%dsedcn)
      END IF

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(para_env)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_ham

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE tb_get_energy(qs_env, tb, energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(qs_energy_type), POINTER                      :: energy

#if defined(__TBLITE)

      INTEGER                                            :: iounit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: scf_section

      NULLIFY (scf_section, logger)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")

      energy%repulsive = SUM(tb%e_rep)
      energy%el_stat = SUM(tb%e_es)
      energy%dispersion = SUM(tb%e_disp)
      energy%dispersion_sc = SUM(tb%e_scd)
      energy%xtb_xb_inter = SUM(tb%e_hal)

      energy%total = energy%core + energy%repulsive + energy%el_stat + energy%dispersion &
                     + energy%dispersion_sc + energy%xtb_xb_inter + energy%kTS &
                     + energy%efield + energy%qmmm_el

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Electrostatic energy:                          ", energy%el_stat, &
            "Self-consistent dispersion energy:             ", energy%dispersion_sc, &
            "Non-self consistent dispersion energy:         ", energy%dispersion
         WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "Correction for halogen bonding:                ", energy%xtb_xb_inter
         IF (ABS(energy%efield) > 1.e-12_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Electric field interaction energy:             ", energy%efield
         END IF
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
         END IF
      END IF
      CALL cp_print_key_finished_output(iounit, logger, scf_section, &
                                        "PRINT%DETAILED_ENERGY")

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(energy)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_energy

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param gto_basis_set ...
!> \param element_symbol ...
!> \param param ...
!> \param occ ...
! **************************************************************************************************
   SUBROUTINE tb_get_basis(tb, gto_basis_set, element_symbol, param, occ)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis_set
      CHARACTER(len=2), INTENT(IN)                       :: element_symbol
      TYPE(xtb_atom_type), POINTER                       :: param
      INTEGER, DIMENSION(5), INTENT(out)                 :: occ

#if defined(__TBLITE)

      CHARACTER(LEN=default_string_length)               :: sng
      INTEGER                                            :: ang, i_type, id_atom, ind_ao, ipgf, iSH, &
                                                            ishell, ityp, maxl, mprim, natorb, &
                                                            nset, nshell
      LOGICAL                                            :: do_ortho

      CALL allocate_gto_basis_set(gto_basis_set)

      !identifying element in the bas data
      CALL symbol_to_number(i_type, element_symbol)
      DO id_atom = 1, tb%mol%nat
         IF (i_type == tb%el_num(id_atom)) EXIT
      END DO
      param%z = i_type
      param%symbol = element_symbol
      param%defined = .TRUE.
      ityp = tb%mol%id(id_atom)

      !getting size information
      nset = tb%calc%bas%nsh_id(ityp)
      nshell = 1
      mprim = 0
      DO ishell = 1, nset
         mprim = MAX(mprim, tb%calc%bas%cgto(ishell, ityp)%nprim)
      END DO
      param%nshell = nset
      natorb = 0

      !write basis set information
      CALL integer_to_string(mprim, sng)
      gto_basis_set%name = element_symbol//"_STO-"//TRIM(sng)//"G"
      gto_basis_set%nset = nset
      CALL reallocate(gto_basis_set%lmax, 1, nset)
      CALL reallocate(gto_basis_set%lmin, 1, nset)
      CALL reallocate(gto_basis_set%npgf, 1, nset)
      CALL reallocate(gto_basis_set%nshell, 1, nset)
      CALL reallocate(gto_basis_set%n, 1, 1, 1, nset)
      CALL reallocate(gto_basis_set%l, 1, 1, 1, nset)
      CALL reallocate(gto_basis_set%zet, 1, mprim, 1, nset)
      CALL reallocate(gto_basis_set%gcc, 1, mprim, 1, 1, 1, nset)

      ind_ao = 0
      maxl = 0
      DO ishell = 1, nset
         ang = tb%calc%bas%cgto(ishell, ityp)%ang
         natorb = natorb + (2*ang + 1)
         param%lval(ishell) = ang
         maxl = MAX(ang, maxl)
         gto_basis_set%lmax(ishell) = ang
         gto_basis_set%lmin(ishell) = ang
         gto_basis_set%npgf(ishell) = tb%calc%bas%cgto(ishell, ityp)%nprim
         gto_basis_set%nshell(ishell) = nshell
         gto_basis_set%n(1, ishell) = ang + 1
         gto_basis_set%l(1, ishell) = ang
         DO ipgf = 1, gto_basis_set%npgf(ishell)
            gto_basis_set%gcc(ipgf, 1, ishell) = tb%calc%bas%cgto(ishell, ityp)%coeff(ipgf)
            gto_basis_set%zet(ipgf, ishell) = tb%calc%bas%cgto(ishell, ityp)%alpha(ipgf)
         END DO
         DO ipgf = 1, (2*ang + 1)
            ind_ao = ind_ao + 1
            param%lao(ind_ao) = ang
            param%nao(ind_ao) = ishell
         END DO
      END DO

      do_ortho = .FALSE.
      CALL process_gto_basis(gto_basis_set, do_ortho, nset, maxl)

      !setting additional values in parameter
      param%rcut = get_cutoff(tb%calc%bas)
      param%natorb = natorb
      param%lmax = maxl !max angular momentum

      !getting occupation
      occ = 0
      DO iSh = 1, MIN(tb%calc%bas%nsh_at(id_atom), 5)
         occ(iSh) = NINT(tb%calc%h0%refocc(iSh, ityp))
         param%occupation(iSh) = occ(iSh)
      END DO
      param%zeff = SUM(occ) !effective core charge

      !set normalization process
      gto_basis_set%norm_type = 3

#else
      MARK_USED(tb)
      MARK_USED(gto_basis_set)
      MARK_USED(element_symbol)
      MARK_USED(param)
      MARK_USED(occ)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_basis

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_tblite_matrices(qs_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_tblite_matrices'

      INTEGER                                            :: handle, maxder, nderivatives, nimg, img, nkind, i, ic, iw, &
                                                            iatom, jatom, ikind, jkind, iset, jset, n1, n2, icol, irow, &
                                                            ia, ib, sgfa, sgfb, atom_a, atom_b, &
                                                            ldsab, nseta, nsetb, natorb_a, natorb_b, sgfa0
      LOGICAL                                            :: found, norml1, norml2, use_arnoldi
      REAL(KIND=dp)                                      :: dr, rr
      INTEGER, DIMENSION(3)                              :: cell
      REAL(KIND=dp)                                      :: hij, shpoly
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: owork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint, hint
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min
      INTEGER, DIMENSION(:), POINTER                     :: npgfa, npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, zeta, zetb, scon_a, scon_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sblock, fblock

      TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
      TYPE(atprop_type), POINTER                            :: atprop
      TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
      TYPE(cp_logger_type), POINTER                         :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER          :: matrix_h, matrix_s, matrix_p, matrix_w
      TYPE(dft_control_type), POINTER                       :: dft_control
      TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
      TYPE(gto_basis_set_type), POINTER                     :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER     :: basis_set_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                              :: nl_iterator
      TYPE(mp_para_env_type), POINTER                       :: para_env
      TYPE(qs_energy_type), POINTER                         :: energy
      TYPE(qs_ks_env_type), POINTER                         :: ks_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set
      TYPE(qs_rho_type), POINTER                            :: rho
      TYPE(tblite_type), POINTER                            :: tb
      TYPE(tb_hamiltonian), POINTER                         :: h0
      TYPE(virial_type), POINTER                            :: virial

      CALL timeset(routineN, handle)

      NULLIFY (ks_env, energy, atomic_kind_set, qs_kind_set)
      NULLIFY (matrix_h, matrix_s, atprop, dft_control)
      NULLIFY (sab_orb, rho, tb)

      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, para_env=para_env, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      matrix_s_kp=matrix_s, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb, &
                      rho=rho, tb_tblite=tb)
      h0 => tb%calc%h0

      !update geometry (required for debug / geometry optimization)
      CALL tb_update_geometry(qs_env, tb)

      nkind = SIZE(atomic_kind_set)
      nderivatives = 0
      IF (calculate_forces) THEN
         nderivatives = 1
         IF (ALLOCATED(tb%grad)) DEALLOCATE (tb%grad)
         ALLOCATE (tb%grad(3, tb%mol%nat))
         IF (ALLOCATED(tb%dsedcn)) DEALLOCATE (tb%dsedcn)
         ALLOCATE (tb%dsedcn(tb%calc%bas%nsh))
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            IF (ALLOCATED(tb%dcndr)) DEALLOCATE (tb%dcndr)
            ALLOCATE (tb%dcndr(3, tb%mol%nat, tb%mol%nat))
            IF (ALLOCATED(tb%dcndL)) DEALLOCATE (tb%dcndL)
            ALLOCATE (tb%dcndL(3, 3, tb%mol%nat))
         END IF
      END IF
      maxder = ncoset(nderivatives)
      nimg = dft_control%nimages

      !intialise hamiltonian
      CALL tb_init_ham(qs_env, tb, para_env)

      ! get density matrtix
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

      ! set up matrices for force calculations
      IF (calculate_forces) THEN
         NULLIFY (force, matrix_w, virial)
         CALL get_qs_env(qs_env=qs_env, &
                         matrix_w_kp=matrix_w, &
                         virial=virial, force=force)

         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_add(matrix_w(1, img)%matrix, matrix_w(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
         tb%use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)

      ! set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ! allocate overlap matrix
      CALL dbcsr_allocate_matrix_set(matrix_s, maxder, nimg)
      CALL create_sab_matrix(ks_env, matrix_s, "OVERLAP MATRIX", basis_set_list, basis_set_list, &
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

      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)
         IF (icol == jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         dr = NORM2(rij(:))

         ic = 1
         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         NULLIFY (fblock)
         CALL dbcsr_get_block_p(matrix=matrix_h(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=fblock, found=found)
         CPASSERT(found)

         ! --------- Overlap
         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         atom_a = atom_of_kind(icol)
         atom_b = atom_of_kind(irow)
         ! basis a
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
         ! basis b
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
         natorb_a = 0
         DO iset = 1, nseta
            natorb_a = natorb_a + (2*basis_set_a%l(1, iset) + 1)
         END DO
         natorb_b = 0
         DO iset = 1, nsetb
            natorb_b = natorb_b + (2*basis_set_b%l(1, iset) + 1)
         END DO
         ALLOCATE (sint(natorb_a, natorb_b, maxder))
         sint = 0.0_dp
         ALLOCATE (hint(natorb_a, natorb_b, maxder))
         hint = 0.0_dp

         !----------------- overlap integrals
         DO iset = 1, nseta
            n1 = npgfa(iset)*(ncoset(la_max(iset)) - ncoset(la_min(iset) - 1))
            sgfa = first_sgfa(1, iset)
            DO jset = 1, nsetb
               IF (set_radius_a(iset) + set_radius_b(jset) < dr) CYCLE
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

         ! update S matrix
         IF (icol <= irow) THEN
            sblock(:, :) = sblock(:, :) + sint(:, :, 1)
         ELSE
            sblock(:, :) = sblock(:, :) + TRANSPOSE(sint(:, :, 1))
         END IF

         ! --------- Hamiltonian
         IF (icol == irow .AND. dr < 0.001_dp) THEN
            !get diagonal F matrix from selfenergy
            n1 = tb%calc%bas%ish_at(icol)
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               hij = tb%selfenergy(n1 + iset)
               DO ia = sgfa, sgfa + nsgfa(iset) - 1
                  hint(ia, ia, 1) = hij
               END DO
            END DO
         ELSE
            !get off-diagonal F matrix
            rr = SQRT(dr/(h0%rad(jkind) + h0%rad(ikind)))
            n1 = tb%calc%bas%ish_at(icol)
            sgfa0 = 1
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               sgfa0 = sgfa0 + tb%calc%bas%cgto(iset, tb%mol%id(icol))%nprim
               n2 = tb%calc%bas%ish_at(irow)
               DO jset = 1, nsetb
                  sgfb = first_sgfb(1, jset)
                  shpoly = (1.0_dp + h0%shpoly(iset, ikind)*rr) &
                           *(1.0_dp + h0%shpoly(jset, jkind)*rr)
                  hij = 0.5_dp*(tb%selfenergy(n1 + iset) + tb%selfenergy(n2 + jset)) &
                        *h0%hscale(iset, jset, ikind, jkind)*shpoly
                  DO ia = sgfa, sgfa + nsgfa(iset) - 1
                     DO ib = sgfb, sgfb + nsgfb(jset) - 1
                        hint(ia, ib, 1) = hij*sint(ia, ib, 1)
                     END DO
                  END DO
               END DO
            END DO
         END IF

         ! update F matrix
         IF (icol <= irow) THEN
            fblock(:, :) = fblock(:, :) + hint(:, :, 1)
         ELSE
            fblock(:, :) = fblock(:, :) + TRANSPOSE(hint(:, :, 1))
         END IF

         DEALLOCATE (oint, owork, sint, hint)

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO img = 1, nimg
         DO i = 1, SIZE(matrix_s, 1)
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
         DO i = 1, SIZE(matrix_h, 1)
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
         END DO
      END DO

      !compute multipole moments for gfn2
      IF (dft_control%qs_control%xtb_control%tblite_method == gfn2xtb) &
         CALL tb_get_multipole(qs_env, tb)

      ! output overlap information
      NULLIFY (logger)
      logger => cp_get_default_logger()
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

      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(calculate_forces)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE build_tblite_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dft_control ...
!> \param tb ...
!> \param calculate_forces ...
!> \param use_rho ...
! **************************************************************************************************
   SUBROUTINE tb_update_charges(qs_env, dft_control, tb, calculate_forces, use_rho)

      TYPE(qs_environment_type), POINTER                                 :: qs_env
      TYPE(dft_control_type), POINTER                                    :: dft_control
      TYPE(tblite_type), POINTER                                         :: tb
      LOGICAL, INTENT(IN)                                                :: calculate_forces
      LOGICAL, INTENT(IN)                                                :: use_rho

#if defined(__TBLITE)

      INTEGER                                            :: iatom, ikind, is, ns, atom_a, ii, im
      INTEGER                                            :: nimg, nkind, nsgf, natorb, na
      INTEGER                                            :: n_atom, max_orb, max_shell
      LOGICAL                                            :: do_dipole, do_quadrupole
      REAL(KIND=dp)                                      :: norm
      INTEGER, DIMENSION(5)                              :: occ
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(25)                             :: nao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ch_atom, ch_shell, ch_ref, ch_orb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, ao_dip, ao_quad

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      ! also compute multipoles needed by GFN2
      do_dipole = .FALSE.
      do_quadrupole = .FALSE.

      ! compute mulliken charges required for charge update
      NULLIFY (particle_set, qs_kind_set, atomic_kind_set)
      CALL get_qs_env(qs_env=qs_env, scf_env=scf_env, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, matrix_s_kp=matrix_s, rho=rho, para_env=para_env)
      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         IF (ASSOCIATED(tb%dipbra)) do_dipole = .TRUE.
         IF (ASSOCIATED(tb%quadbra)) do_quadrupole = .TRUE.
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF
      n_atom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      nimg = dft_control%nimages
      CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
      ALLOCATE (aocg(nsgf, n_atom))
      aocg = 0.0_dp
      IF (do_dipole) ALLOCATE (ao_dip(n_atom, dip_n))
      IF (do_quadrupole) ALLOCATE (ao_quad(n_atom, quad_n))
      max_orb = 0
      max_shell = 0
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, natorb=natorb)
         max_orb = MAX(max_orb, natorb)
      END DO
      DO is = 1, n_atom
         max_shell = MAX(max_shell, tb%calc%bas%nsh_at(is))
      END DO
      ALLOCATE (ch_atom(n_atom, 1), ch_shell(n_atom, max_shell))
      ALLOCATE (ch_orb(max_orb, n_atom), ch_ref(max_orb, n_atom))
      ch_atom = 0.0_dp
      ch_shell = 0.0_dp
      ch_orb = 0.0_dp
      ch_ref = 0.0_dp
      IF (nimg > 1) THEN
         CALL ao_charges(matrix_p, matrix_s, aocg, para_env)
         IF (do_dipole .OR. do_quadrupole) THEN
            CPABORT("missing contraction with density matrix for multiple k-points")
         END IF
      ELSE
         NULLIFY (p_matrix, s_matrix)
         p_matrix => matrix_p(:, 1)
         s_matrix => matrix_s(1, 1)%matrix
         CALL ao_charges(p_matrix, s_matrix, aocg, para_env)
         IF (do_dipole) THEN
            DO im = 1, dip_n
               CALL contract_dens(p_matrix, tb%dipbra(im)%matrix, tb%dipket(im)%matrix, ao_dip(:, im), para_env)
            END DO
         END IF
         IF (do_quadrupole) THEN
            DO im = 1, quad_n
               CALL contract_dens(p_matrix, tb%quadbra(im)%matrix, tb%quadket(im)%matrix, ao_quad(:, im), para_env)
            END DO
         END IF
      END IF
      NULLIFY (xtb_kind)
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, nao=nao, occupation=occ)
         DO iatom = 1, na
            atom_a = atomic_kind_set(ikind)%atom_list(iatom)
            DO is = 1, natorb
               ns = lao(is) + 1
               norm = 2*lao(is) + 1
               ch_ref(is, atom_a) = tb%calc%h0%refocc(nao(is), ikind)/norm
               ch_orb(is, atom_a) = aocg(is, atom_a) - ch_ref(is, atom_a)
               ch_shell(atom_a, ns) = ch_orb(is, atom_a) + ch_shell(atom_a, ns)
            END DO
            ch_atom(atom_a, 1) = SUM(ch_orb(:, atom_a))
         END DO
      END DO
      DEALLOCATE (aocg)

      ! charge mixing
      IF (dft_control%qs_control%do_ls_scf) THEN
         !
      ELSE
         CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                            ch_shell, para_env, scf_env%iter_count)
      END IF

      !setting new wave function
      CALL tb%pot%reset
      tb%e_es = 0.0_dp
      tb%e_scd = 0.0_dp
      DO iatom = 1, n_atom
         ii = tb%calc%bas%ish_at(iatom)
         DO is = 1, tb%calc%bas%nsh_at(iatom)
            tb%wfn%qsh(ii + is, 1) = -ch_shell(iatom, is)
         END DO
         tb%wfn%qat(iatom, 1) = -ch_atom(iatom, 1)
      END DO

      IF (do_dipole) THEN
         DO iatom = 1, n_atom
            DO im = 1, dip_n
               tb%wfn%dpat(im, iatom, 1) = -ao_dip(iatom, im)
            END DO
         END DO
         DEALLOCATE (ao_dip)
      END IF
      IF (do_quadrupole) THEN
         DO iatom = 1, n_atom
            DO im = 1, quad_n
               tb%wfn%qpat(im, iatom, 1) = -ao_quad(iatom, im)
            END DO
         END DO
         DEALLOCATE (ao_quad)
      END IF

      IF (ALLOCATED(tb%calc%coulomb)) THEN
         CALL tb%calc%coulomb%get_potential(tb%mol, tb%cache, tb%wfn, tb%pot)
         CALL tb%calc%coulomb%get_energy(tb%mol, tb%cache, tb%wfn, tb%e_es)
      END IF
      IF (ALLOCATED(tb%calc%dispersion)) THEN
         CALL tb%calc%dispersion%get_potential(tb%mol, tb%dcache, tb%wfn, tb%pot)
         CALL tb%calc%dispersion%get_energy(tb%mol, tb%dcache, tb%wfn, tb%e_scd)
      END IF

      IF (calculate_forces) THEN
         IF (ALLOCATED(tb%calc%coulomb)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%coulomb%get_gradient(tb%mol, tb%cache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_grad2force(qs_env, tb, para_env, 3)
         END IF

         IF (ALLOCATED(tb%calc%dispersion)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%dispersion%get_gradient(tb%mol, tb%dcache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_grad2force(qs_env, tb, para_env, 2)
         END IF
      END IF

      DEALLOCATE (ch_atom, ch_shell, ch_orb, ch_ref)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(dft_control)
      MARK_USED(calculate_forces)
      MARK_USED(use_rho)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_update_charges

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE tb_ham_add_coulomb(qs_env, tb, dft_control)

      TYPE(qs_environment_type), POINTER                       :: qs_env
      TYPE(tblite_type), POINTER                               :: tb
      TYPE(dft_control_type), POINTER                          :: dft_control

#if defined(__TBLITE)

      INTEGER                                            :: ikind, jkind, iatom, jatom, icol, irow
      INTEGER                                            :: ic, is, nimg, ni, nj, i, j
      INTEGER                                            :: la, lb, za, zb
      LOGICAL                                            :: found
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(25)                             :: naoa, naob
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, sum_shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ashift, bshift
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ksblock, sblock
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dip_ket1, dip_ket2, dip_ket3, &
                                                            dip_bra1, dip_bra2, dip_bra3
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: quad_ket1, quad_ket2, quad_ket3, &
                                                            quad_ket4, quad_ket5, quad_ket6, &
                                                            quad_bra1, quad_bra2, quad_bra3, &
                                                            quad_bra4, quad_bra5, quad_bra6

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b

      nimg = dft_control%nimages

      NULLIFY (matrix_s, ks_matrix, n_list, qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, sab_orb=n_list, matrix_s_kp=matrix_s, matrix_ks_kp=ks_matrix, qs_kind_set=qs_kind_set)

      !creating sum of shell lists
      ALLOCATE (sum_shell(tb%mol%nat))
      i = 0
      DO j = 1, tb%mol%nat
         sum_shell(j) = i
         i = i + tb%calc%bas%nsh_at(j)
      END DO

      IF (nimg == 1) THEN
         ! no k-points; all matrices have been transformed to periodic bsf
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  kind_of=kind_of)
         CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

            ikind = kind_of(irow)
            jkind = kind_of(icol)

            ! atomic parameters
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa)
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob)

            ni = SIZE(sblock, 1)
            ALLOCATE (ashift(ni, ni))
            ashift = 0.0_dp
            DO i = 1, ni
               la = naoa(i) + sum_shell(irow)
               ashift(i, i) = tb%pot%vsh(la, 1)
            END DO

            nj = SIZE(sblock, 2)
            ALLOCATE (bshift(nj, nj))
            bshift = 0.0_dp
            DO j = 1, nj
               lb = naob(j) + sum_shell(icol)
               bshift(j, j) = tb%pot%vsh(lb, 1)
            END DO

            DO is = 1, SIZE(ks_matrix, 1)
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock - 0.5_dp*(MATMUL(ashift, sblock) &
                                           + MATMUL(sblock, bshift))
               ksblock = ksblock - 0.5_dp*(tb%pot%vat(irow, 1) &
                                           + tb%pot%vat(icol, 1))*sblock
            END DO
            DEALLOCATE (ashift, bshift)
         END DO
         CALL dbcsr_iterator_stop(iter)

         IF (ASSOCIATED(tb%dipbra)) THEN
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

               NULLIFY (dip_bra1, dip_bra2, dip_bra3)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra3, found=found)
               CPASSERT(found)
               NULLIFY (dip_ket1, dip_ket2, dip_ket3)
               CALL dbcsr_get_block_p(matrix=tb%dipket(1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket3, found=found)
               CPASSERT(found)

               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock - 0.5_dp*(dip_ket1*tb%pot%vdp(1, irow, 1) &
                                              + dip_ket2*tb%pot%vdp(2, irow, 1) &
                                              + dip_ket3*tb%pot%vdp(3, irow, 1) &
                                              + dip_bra1*tb%pot%vdp(1, icol, 1) &
                                              + dip_bra2*tb%pot%vdp(2, icol, 1) &
                                              + dip_bra3*tb%pot%vdp(3, icol, 1))
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
         END IF

         IF (ASSOCIATED(tb%quadbra)) THEN
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

               NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra6, found=found)
               CPASSERT(found)

               NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
               CALL dbcsr_get_block_p(matrix=tb%quadket(1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket6, found=found)
               CPASSERT(found)

               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)

                  ksblock = ksblock - 0.5_dp*(quad_ket1*tb%pot%vqp(1, irow, 1) &
                                              + quad_ket2*tb%pot%vqp(2, irow, 1) &
                                              + quad_ket3*tb%pot%vqp(3, irow, 1) &
                                              + quad_ket4*tb%pot%vqp(4, irow, 1) &
                                              + quad_ket5*tb%pot%vqp(5, irow, 1) &
                                              + quad_ket6*tb%pot%vqp(6, irow, 1) &
                                              + quad_bra1*tb%pot%vqp(1, icol, 1) &
                                              + quad_bra2*tb%pot%vqp(2, icol, 1) &
                                              + quad_bra3*tb%pot%vqp(3, icol, 1) &
                                              + quad_bra4*tb%pot%vqp(4, icol, 1) &
                                              + quad_bra5*tb%pot%vqp(5, icol, 1) &
                                              + quad_bra6*tb%pot%vqp(6, icol, 1))
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
         END IF

      ELSE
         CPABORT("GFN methods with k-points not tested")
         NULLIFY (kpoints)
         CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
         NULLIFY (cell_to_index)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)

         CALL neighbor_list_iterator_create(nl_iterator, n_list)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rij, cell=cellind)

            icol = MAX(iatom, jatom)
            irow = MIN(iatom, jatom)

            IF (icol == jatom) THEN
               rij = -rij
               i = ikind
               ikind = jkind
               jkind = i
            END IF

            ic = cell_to_index(cellind(1), cellind(2), cellind(3))
            CPASSERT(ic > 0)

            NULLIFY (sblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                   row=irow, col=icol, block=sblock, found=found)
            CPASSERT(found)

            ! atomic parameters
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa)
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob)

            ni = SIZE(sblock, 1)
            ALLOCATE (ashift(ni, ni))
            ashift = 0.0_dp
            DO i = 1, ni
               la = naoa(i) + sum_shell(irow)
               ashift(i, i) = tb%pot%vsh(la, 1)
            END DO

            nj = SIZE(sblock, 2)
            ALLOCATE (bshift(nj, nj))
            bshift = 0.0_dp
            DO j = 1, nj
               lb = naob(j) + sum_shell(icol)
               bshift(j, j) = tb%pot%vsh(lb, 1)
            END DO

            DO is = 1, SIZE(ks_matrix, 1)
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock - 0.5_dp*(MATMUL(ashift, sblock) &
                                           + MATMUL(sblock, bshift))
               ksblock = ksblock - 0.5_dp*(tb%pot%vat(irow, 1) &
                                           + tb%pot%vat(icol, 1))*sblock
            END DO
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         IF (ASSOCIATED(tb%dipbra)) THEN
            CALL neighbor_list_iterator_create(nl_iterator, n_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, &
                                      iatom=iatom, jatom=jatom, cell=cellind)
               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)

               NULLIFY (dip_bra1, dip_bra2, dip_bra3)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra3, found=found)
               CPASSERT(found)
               NULLIFY (dip_ket1, dip_ket2, dip_ket3)
               CALL dbcsr_get_block_p(matrix=tb%dipket(1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket3, found=found)
               CPASSERT(found)

               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock - 0.5_dp*(dip_ket1*tb%pot%vdp(1, irow, 1) &
                                              + dip_ket2*tb%pot%vdp(2, irow, 1) &
                                              + dip_ket3*tb%pot%vdp(3, irow, 1) &
                                              + dip_bra1*tb%pot%vdp(1, icol, 1) &
                                              + dip_bra2*tb%pot%vdp(2, icol, 1) &
                                              + dip_bra3*tb%pot%vdp(3, icol, 1))
               END DO
            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

         IF (ASSOCIATED(tb%quadbra)) THEN
            CALL neighbor_list_iterator_create(nl_iterator, n_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, &
                                      iatom=iatom, jatom=jatom, cell=cellind)
               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)

               NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra6, found=found)
               CPASSERT(found)

               NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
               CALL dbcsr_get_block_p(matrix=tb%quadket(1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket6, found=found)
               CPASSERT(found)

               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)

                  ksblock = ksblock - 0.5_dp*(quad_ket1*tb%pot%vqp(1, irow, 1) &
                                              + quad_ket2*tb%pot%vqp(2, irow, 1) &
                                              + quad_ket3*tb%pot%vqp(3, irow, 1) &
                                              + quad_ket4*tb%pot%vqp(4, irow, 1) &
                                              + quad_ket5*tb%pot%vqp(5, irow, 1) &
                                              + quad_ket6*tb%pot%vqp(6, irow, 1) &
                                              + quad_bra1*tb%pot%vqp(1, icol, 1) &
                                              + quad_bra2*tb%pot%vqp(2, icol, 1) &
                                              + quad_bra3*tb%pot%vqp(3, icol, 1) &
                                              + quad_bra4*tb%pot%vqp(4, icol, 1) &
                                              + quad_bra5*tb%pot%vqp(5, icol, 1) &
                                              + quad_bra6*tb%pot%vqp(6, icol, 1))
               END DO
            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

      END IF

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(dft_control)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_ham_add_coulomb

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_get_multipole(qs_env, tb)

      TYPE(qs_environment_type), POINTER                       :: qs_env
      TYPE(tblite_type), POINTER                               :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tb_get_multipole'

      INTEGER                                     :: ikind, jkind, iatom, jatom, icol, irow, iset, jset, ityp, jtyp
      INTEGER                                     :: nkind, natom, handle, nimg, i, inda, indb
      INTEGER                                     :: atom_a, atom_b, nseta, nsetb, ia, ib, ij
      LOGICAL                                     :: found
      REAL(KIND=dp)                               :: r2
      REAL(KIND=dp), DIMENSION(3)                 :: rij
      INTEGER, DIMENSION(:), POINTER              :: la_max, lb_max
      INTEGER, DIMENSION(:), POINTER              :: nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER           :: first_sgfa, first_sgfb
      INTEGER, ALLOCATABLE                        :: atom_of_kind(:)
      REAL(KIND=dp), ALLOCATABLE                  :: stmp(:)
      REAL(KIND=dp), ALLOCATABLE                  :: dtmp(:, :), qtmp(:, :), dtmpj(:, :), qtmpj(:, :)
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: dip_ket1, dip_ket2, dip_ket3, &
                                                     dip_bra1, dip_bra2, dip_bra3
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: quad_ket1, quad_ket2, quad_ket3, &
                                                     quad_ket4, quad_ket5, quad_ket6, &
                                                     quad_bra1, quad_bra2, quad_bra3, &
                                                     quad_bra4, quad_bra5, quad_bra6

      TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                       :: dft_control
      TYPE(gto_basis_set_type), POINTER                     :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER     :: basis_set_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                              :: nl_iterator
      TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      NULLIFY (atomic_kind_set, qs_kind_set, sab_orb, particle_set)
      NULLIFY (dft_control, matrix_s)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      sab_orb=sab_orb, &
                      particle_set=particle_set, &
                      dft_control=dft_control, &
                      matrix_s_kp=matrix_s)
      natom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      nimg = dft_control%nimages

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)

      !set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ALLOCATE (stmp(msao(tb%calc%bas%maxl)**2))
      ALLOCATE (dtmp(dip_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (qtmp(quad_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (dtmpj(dip_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (qtmpj(quad_n, msao(tb%calc%bas%maxl)**2))

      ! allocate dipole/quadrupole moment matrix elemnts
      CALL dbcsr_allocate_matrix_set(tb%dipbra, dip_n)
      CALL dbcsr_allocate_matrix_set(tb%dipket, dip_n)
      CALL dbcsr_allocate_matrix_set(tb%quadbra, quad_n)
      CALL dbcsr_allocate_matrix_set(tb%quadket, quad_n)
      DO i = 1, dip_n
         ALLOCATE (tb%dipbra(i)%matrix)
         ALLOCATE (tb%dipket(i)%matrix)
         CALL dbcsr_create(tb%dipbra(i)%matrix, template=matrix_s(1, 1)%matrix, &
                           name="DIPOLE BRAMATRIX")
         CALL dbcsr_create(tb%dipket(i)%matrix, template=matrix_s(1, 1)%matrix, &
                           name="DIPOLE KETMATRIX")
         CALL cp_dbcsr_alloc_block_from_nbl(tb%dipbra(i)%matrix, sab_orb)
         CALL cp_dbcsr_alloc_block_from_nbl(tb%dipket(i)%matrix, sab_orb)
      END DO
      DO i = 1, quad_n
         ALLOCATE (tb%quadbra(i)%matrix)
         ALLOCATE (tb%quadket(i)%matrix)
         CALL dbcsr_create(tb%quadbra(i)%matrix, template=matrix_s(1, 1)%matrix, &
                           name="QUADRUPOLE BRAMATRIX")
         CALL dbcsr_create(tb%quadket(i)%matrix, template=matrix_s(1, 1)%matrix, &
                           name="QUADRUPOLE KETMATRIX")
         CALL cp_dbcsr_alloc_block_from_nbl(tb%quadbra(i)%matrix, sab_orb)
         CALL cp_dbcsr_alloc_block_from_nbl(tb%quadket(i)%matrix, sab_orb)
      END DO

      !loop over all atom pairs with a non-zero overlap (sab_orb)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)

         r2 = NORM2(rij(:))**2

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         IF (icol == jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         ityp = tb%mol%id(icol)
         jtyp = tb%mol%id(irow)

         NULLIFY (dip_bra1, dip_bra2, dip_bra3)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(1)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(2)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(3)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra3, found=found)
         CPASSERT(found)

         NULLIFY (dip_ket1, dip_ket2, dip_ket3)
         CALL dbcsr_get_block_p(matrix=tb%dipket(1)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipket(2)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipket(3)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket3, found=found)
         CPASSERT(found)

         NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(1)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(2)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(3)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra3, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(4)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra4, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(5)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra5, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(6)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra6, found=found)
         CPASSERT(found)

         NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
         CALL dbcsr_get_block_p(matrix=tb%quadket(1)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(2)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(3)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket3, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(4)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket4, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(5)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket5, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(6)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket6, found=found)
         CPASSERT(found)

         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         atom_a = atom_of_kind(icol)
         atom_b = atom_of_kind(irow)
         ! basis a
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         ! basis b
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set

         ! --------- Hamiltonian
         IF (icol == irow) THEN
            DO iset = 1, nseta
               DO jset = 1, nsetb
                  CALL multipole_cgto(tb%calc%bas%cgto(jset, ityp), tb%calc%bas%cgto(iset, ityp), &
                       & r2, -rij, tb%calc%bas%intcut, stmp, dtmp, qtmp)

                  DO inda = 1, nsgfa(iset)
                     ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                     DO indb = 1, nsgfb(jset)
                        ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb
                        ij = indb + nsgfb(jset)*(inda - 1)

                        dip_ket1(ib, ia) = dtmp(1, ij)
                        dip_ket2(ib, ia) = dtmp(2, ij)
                        dip_ket3(ib, ia) = dtmp(3, ij)

                        quad_ket1(ib, ia) = qtmp(1, ij)
                        quad_ket2(ib, ia) = qtmp(2, ij)
                        quad_ket3(ib, ia) = qtmp(3, ij)
                        quad_ket4(ib, ia) = qtmp(4, ij)
                        quad_ket5(ib, ia) = qtmp(5, ij)
                        quad_ket6(ib, ia) = qtmp(6, ij)

                        dip_bra1(ib, ia) = dtmp(1, ij)
                        dip_bra2(ib, ia) = dtmp(2, ij)
                        dip_bra3(ib, ia) = dtmp(3, ij)

                        quad_bra1(ib, ia) = qtmp(1, ij)
                        quad_bra2(ib, ia) = qtmp(2, ij)
                        quad_bra3(ib, ia) = qtmp(3, ij)
                        quad_bra4(ib, ia) = qtmp(4, ij)
                        quad_bra5(ib, ia) = qtmp(5, ij)
                        quad_bra6(ib, ia) = qtmp(6, ij)
                     END DO
                  END DO
               END DO
            END DO
         ELSE
            DO iset = 1, nseta
               DO jset = 1, nsetb
                  CALL multipole_cgto(tb%calc%bas%cgto(jset, jtyp), tb%calc%bas%cgto(iset, ityp), &
                      & r2, -rij, tb%calc%bas%intcut, stmp, dtmp, qtmp)
                  CALL multipole_cgto(tb%calc%bas%cgto(iset, ityp), tb%calc%bas%cgto(jset, jtyp), &
                      & r2, rij, tb%calc%bas%intcut, stmp, dtmpj, qtmpj)

                  DO inda = 1, nsgfa(iset)
                     ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                     DO indb = 1, nsgfb(jset)
                        ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb

                        ij = indb + nsgfb(jset)*(inda - 1)

                        dip_bra1(ib, ia) = dtmp(1, ij)
                        dip_bra2(ib, ia) = dtmp(2, ij)
                        dip_bra3(ib, ia) = dtmp(3, ij)

                        quad_bra1(ib, ia) = qtmp(1, ij)
                        quad_bra2(ib, ia) = qtmp(2, ij)
                        quad_bra3(ib, ia) = qtmp(3, ij)
                        quad_bra4(ib, ia) = qtmp(4, ij)
                        quad_bra5(ib, ia) = qtmp(5, ij)
                        quad_bra6(ib, ia) = qtmp(6, ij)

                        ij = inda + nsgfa(iset)*(indb - 1)

                        dip_ket1(ib, ia) = dtmpj(1, ij)
                        dip_ket2(ib, ia) = dtmpj(2, ij)
                        dip_ket3(ib, ia) = dtmpj(3, ij)

                        quad_ket1(ib, ia) = qtmpj(1, ij)
                        quad_ket2(ib, ia) = qtmpj(2, ij)
                        quad_ket3(ib, ia) = qtmpj(3, ij)
                        quad_ket4(ib, ia) = qtmpj(4, ij)
                        quad_ket5(ib, ia) = qtmpj(5, ij)
                        quad_ket6(ib, ia) = qtmpj(6, ij)
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, dip_n
         CALL dbcsr_finalize(tb%dipbra(i)%matrix)
         CALL dbcsr_finalize(tb%dipket(i)%matrix)
      END DO
      DO i = 1, quad_n
         CALL dbcsr_finalize(tb%quadbra(i)%matrix)
         CALL dbcsr_finalize(tb%quadket(i)%matrix)
      END DO

      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_multipole

! **************************************************************************************************
!> \brief compute the mulliken properties (AO resolved)
!> \param p_matrix ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param output ...
!> \param para_env ...
!> \par History
!>      adapted from ao_charges_2
!> \note
! **************************************************************************************************
   SUBROUTINE contract_dens(p_matrix, bra_mat, ket_mat, output, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_type), POINTER                          :: bra_mat, ket_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'contract_dens'

      INTEGER                                            :: handle, i, iblock_col, iblock_row, &
                                                            ispin, j, nspin
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: bra, ket, p_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      nspin = SIZE(p_matrix)
      output = 0.0_dp
      DO ispin = 1, nspin
         CALL dbcsr_iterator_start(iter, bra_mat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            NULLIFY (p_block, bra, ket)

            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, bra)
            CALL dbcsr_get_block_p(matrix=p_matrix(ispin)%matrix, &
                                   row=iblock_row, col=iblock_col, BLOCK=p_block, found=found)
            IF (.NOT. found) CYCLE
            CALL dbcsr_get_block_p(matrix=ket_mat, &
                                   row=iblock_row, col=iblock_col, BLOCK=ket, found=found)
            IF (.NOT. found) CPABORT("missing block")

            IF (.NOT. (ASSOCIATED(bra) .AND. ASSOCIATED(p_block))) CYCLE
            DO j = 1, SIZE(p_block, 1)
               DO i = 1, SIZE(p_block, 2)
                  output(iblock_row) = output(iblock_row) + p_block(j, i)*ket(j, i)
               END DO
            END DO
            IF (iblock_col /= iblock_row) THEN
               DO j = 1, SIZE(p_block, 1)
                  DO i = 1, SIZE(p_block, 2)
                     output(iblock_col) = output(iblock_col) + p_block(j, i)*bra(j, i)
                  END DO
               END DO
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO

      CALL para_env%sum(output)
      CALL timestop(handle)

   END SUBROUTINE contract_dens

! **************************************************************************************************
!> \brief save gradient to force
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
!> \param ityp ...
!> \note
! **************************************************************************************************
   SUBROUTINE tb_grad2force(qs_env, tb, para_env, ityp)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env
      INTEGER                                            :: ityp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tb_grad2force'

      INTEGER                                            :: atoma, handle, iatom, ikind, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      CALL timeset(routineN, handle)

      NULLIFY (force, atomic_kind_set)
      CALL get_qs_env(qs_env=qs_env, force=force, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      natom = SIZE(particle_set)

      SELECT CASE (ityp)
      CASE DEFAULT
         CPABORT("unknown force type")
      CASE (0)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%all_potential(:, atoma) = &
               force(ikind)%all_potential(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (1)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%repulsive(:, atoma) = &
               force(ikind)%repulsive(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (2)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%dispersion(:, atoma) = &
               force(ikind)%dispersion(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (3)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%rho_elec(:, atoma) = &
               force(ikind)%rho_elec(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (4)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%overlap(:, atoma) = &
               force(ikind)%overlap(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (5)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%efield(:, atoma) = &
               force(ikind)%efield(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE tb_grad2force

! **************************************************************************************************
!> \brief set gradient to zero
!> \param qs_env ...
!> \note
! **************************************************************************************************
   SUBROUTINE tb_zero_force(qs_env)

      TYPE(qs_environment_type)                          :: qs_env

      INTEGER                                            :: iatom, ikind, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      NULLIFY (force, atomic_kind_set)
      CALL get_qs_env(qs_env=qs_env, force=force, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               kind_of=kind_of)

      natom = SIZE(particle_set)

      DO iatom = 1, natom
         ikind = kind_of(iatom)
         force(ikind)%all_potential = 0.0_dp
         force(ikind)%repulsive = 0.0_dp
         force(ikind)%dispersion = 0.0_dp
         force(ikind)%rho_elec = 0.0_dp
         force(ikind)%overlap = 0.0_dp
         force(ikind)%efield = 0.0_dp
      END DO

   END SUBROUTINE tb_zero_force

! **************************************************************************************************
!> \brief compute the derivative of the energy
!> \param qs_env ...
!> \param use_rho ...
! **************************************************************************************************
   SUBROUTINE tb_derive_dH_diag(qs_env, use_rho)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: use_rho

#if defined(__TBLITE)
      INTEGER                                            :: i, iatom, ic, ikind, ind, is, ish, &
                                                            jatom, jkind
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tblite_type), POINTER                         :: tb

      ! compute mulliken charges required for charge update
      NULLIFY (scf_env, rho, tb, sab_orb, para_env, kpoints)
      CALL get_qs_env(qs_env=qs_env, scf_env=scf_env, rho=rho, tb_tblite=tb, &
                      sab_orb=sab_orb, para_env=para_env, kpoints=kpoints)
      NULLIFY (cell_to_index)
      CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)

      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF

      ALLOCATE (dE(tb%mol%nat))

      dE = 0.0_dp
      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, cell=cellind)

         IF (iatom /= jatom) CYCLE

         is = tb%calc%bas%ish_at(iatom)

         ic = cell_to_index(cellind(1), cellind(2), cellind(3))
         CPASSERT(ic > 0)
         CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                row=iatom, col=jatom, BLOCK=pblock, found=found)

         IF (.NOT. found) CPABORT("block not found")

         ind = 0
         DO ish = 1, tb%calc%bas%nsh_id(ikind)
            DO i = 1, msao(tb%calc%bas%cgto(ish, ikind)%ang)
               ind = ind + 1
               dE(iatom) = dE(iatom) + tb%dsedcn(is + ish)*pblock(ind, ind)
            END DO
         END DO
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CALL para_env%sum(dE)

      tb%grad = 0.0_dp
      CALL tb_add_grad(tb%grad, tb%dcndr, dE, tb%mol%nat)
      IF (tb%use_virial) CALL tb_add_sig(tb%sigma, tb%dcndL, dE, tb%mol%nat)
      CALL tb_grad2force(qs_env, tb, para_env, 4)
      DEALLOCATE (dE)

#else
      MARK_USED(qs_env)
      MARK_USED(use_rho)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_derive_dH_diag

! **************************************************************************************************
!> \brief compute the derivative of the energy
!> \param qs_env ...
!> \param use_rho ...
! **************************************************************************************************
   SUBROUTINE tb_derive_dH_off(qs_env, use_rho)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: use_rho

#if defined(__TBLITE)
      INTEGER                                            :: i, ij, iatom, ic, icol, ikind, ni, nj, nkind, nel, &
                                                            sgfi, sgfj, ityp, jatom, jkind, jrow, jtyp, iset, jset, &
                                                            nseti, nsetj, ia, ib, inda, indb
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:), POINTER                     :: nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: r2, itemp, jtemp, rr, hij, shpoly, dshpoly, idHdc, jdHdc, &
                                                            scal, hp, i_a_shift, j_a_shift, ishift, jshift
      REAL(KIND=dp), DIMENSION(3)                        :: rij, dgrad
      REAL(KIND=dp), DIMENSION(3, 3)                      :: hsigma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE, t_ov, idip, jdip, iquad, jquad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: t_dip, t_quad, t_d_ov
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: t_i_dip, t_i_quad, t_j_dip, t_j_quad
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock, wblock

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_w
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER            :: basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tb_hamiltonian), POINTER                      :: h0
      TYPE(tblite_type), POINTER                         :: tb

      ! compute mulliken charges required for charge update
      NULLIFY (scf_env, rho, tb, sab_orb, para_env, kpoints, matrix_w)
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      scf_env=scf_env, &
                      rho=rho, &
                      tb_tblite=tb, &
                      sab_orb=sab_orb, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set, &
                      kpoints=kpoints, &
                      matrix_w_kp=matrix_w)
      NULLIFY (cell_to_index)
      CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)

      h0 => tb%calc%h0

      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF

      ! set up basis set lists
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ALLOCATE (dE(tb%mol%nat))

      nel = msao(tb%calc%bas%maxl)**2
      ALLOCATE (t_ov(nel))
      ALLOCATE (t_d_ov(3, nel))
      ALLOCATE (t_dip(dip_n, nel))
      ALLOCATE (t_i_dip(3, dip_n, nel), t_j_dip(3, dip_n, nel))
      ALLOCATE (t_quad(quad_n, nel))
      ALLOCATE (t_i_quad(3, quad_n, nel), t_j_quad(3, quad_n, nel))

      ALLOCATE (idip(dip_n), jdip(dip_n))
      ALLOCATE (iquad(quad_n), jquad(quad_n))

      dE = 0.0_dp
      tb%grad = 0.0_dp
      hsigma = 0.0_dp
      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cellind)

         IF (iatom == jatom) CYCLE

         icol = MAX(iatom, jatom)
         jrow = MIN(iatom, jatom)

         IF (icol == jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         ityp = tb%mol%id(icol)
         jtyp = tb%mol%id(jrow)

         r2 = NORM2(rij(:))
         rr = SQRT(r2/(h0%rad(ikind) + h0%rad(jkind)))
         r2 = r2**2

         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         first_sgfa => basis_set_a%first_sgf
         nsgfa => basis_set_a%nsgf_set
         nseti = basis_set_a%nset
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         first_sgfb => basis_set_b%first_sgf
         nsgfb => basis_set_b%nsgf_set
         nsetj = basis_set_b%nset

         ic = cell_to_index(cellind(1), cellind(2), cellind(3))
         CPASSERT(ic > 0)
         NULLIFY (pblock)
         CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                row=jrow, col=icol, BLOCK=pblock, found=found)
         IF (.NOT. found) CPABORT("pblock not found")

         NULLIFY (wblock)
         CALL dbcsr_get_block_p(matrix=matrix_w(1, ic)%matrix, &
                                row=jrow, col=icol, block=wblock, found=found)
         CPASSERT(found)

         i_a_shift = tb%pot%vat(icol, 1)
         j_a_shift = tb%pot%vat(jrow, 1)
         idip(:) = tb%pot%vdp(:, icol, 1)
         jdip(:) = tb%pot%vdp(:, jrow, 1)
         iquad(:) = tb%pot%vqp(:, icol, 1)
         jquad(:) = tb%pot%vqp(:, jrow, 1)

         ni = tb%calc%bas%ish_at(icol)
         DO iset = 1, nseti
            sgfi = first_sgfa(1, iset)
            ishift = i_a_shift + tb%pot%vsh(ni + iset, 1)
            nj = tb%calc%bas%ish_at(jrow)
            DO jset = 1, nsetj
               sgfj = first_sgfb(1, jset)
               jshift = j_a_shift + tb%pot%vsh(nj + jset, 1)

               !get integrals and derivatives
               CALL multipole_grad_cgto(tb%calc%bas%cgto(iset, ityp), tb%calc%bas%cgto(jset, jtyp), &
                  & r2, rij, tb%calc%bas%intcut, t_ov, t_dip, t_quad, t_d_ov, t_i_dip, t_i_quad, &
                  & t_j_dip, t_j_quad)

               shpoly = (1.0_dp + h0%shpoly(iset, ikind)*rr) &
                        *(1.0_dp + h0%shpoly(jset, jkind)*rr)
               dshpoly = ((1.0_dp + h0%shpoly(iset, ikind)*rr)*h0%shpoly(jset, jkind)*rr &
                          + (1.0_dp + h0%shpoly(jset, jkind)*rr)*h0%shpoly(iset, ikind)*rr) &
                         *0.5_dp/r2
               scal = h0%hscale(iset, jset, ikind, jkind)
               hij = 0.5_dp*(tb%selfenergy(ni + iset) + tb%selfenergy(nj + jset))*scal

               idHdc = tb%dsedcn(ni + iset)*shpoly*scal
               jdHdc = tb%dsedcn(nj + jset)*shpoly*scal

               itemp = 0.0_dp
               jtemp = 0.0_dp
               dgrad = 0.0_dp
               DO inda = 1, nsgfa(iset)
                  ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                  DO indb = 1, nsgfb(jset)
                     ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb

                     ij = inda + nsgfa(iset)*(indb - 1)

                     itemp = itemp + idHdc*pblock(ib, ia)*t_ov(ij)
                     jtemp = jtemp + jdHdc*pblock(ib, ia)*t_ov(ij)

                     hp = 2*hij*pblock(ib, ia)
                     dgrad(:) = dgrad(:) &
                                - (ishift + jshift)*pblock(ib, ia)*t_d_ov(:, ij) &
                                - 2*wblock(ib, ia)*t_d_ov(:, ij) &
                                + hp*shpoly*t_d_ov(:, ij) &
                                + hp*dshpoly*t_ov(ij)*rij &
                                - pblock(ib, ia)*( &
                                MATMUL(t_i_dip(:, :, ij), idip) &
                                + MATMUL(t_j_dip(:, :, ij), jdip) &
                                + MATMUL(t_i_quad(:, :, ij), iquad) &
                                + MATMUL(t_j_quad(:, :, ij), jquad))

                  END DO
               END DO
               dE(icol) = dE(icol) + itemp
               dE(jrow) = dE(jrow) + jtemp
               tb%grad(:, icol) = tb%grad(:, icol) - dgrad
               tb%grad(:, jrow) = tb%grad(:, jrow) + dgrad
               IF ((tb%use_virial) .AND. (icol == jrow)) THEN
                  DO ia = 1, 3
                     DO ib = 1, 3
                        hsigma(ia, ib) = hsigma(ia, ib) + 0.25_dp*(rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                     END DO
                  END DO
               ELSE
                  DO ia = 1, 3
                     DO ib = 1, 3
                        hsigma(ia, ib) = hsigma(ia, ib) + 0.50_dp*(rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL para_env%sum(hsigma)
      tb%sigma = tb%sigma + hsigma

      CALL para_env%sum(dE)
      CALL para_env%sum(tb%grad)
      CALL tb_add_grad(tb%grad, tb%dcndr, dE, tb%mol%nat)
      IF (tb%use_virial) CALL tb_add_sig(tb%sigma, tb%dcndL, dE, tb%mol%nat)
      CALL tb_grad2force(qs_env, tb, para_env, 4)

      DEALLOCATE (dE)
      DEALLOCATE (basis_set_list)

      DEALLOCATE (t_ov, t_d_ov)
      DEALLOCATE (t_dip, t_i_dip, t_j_dip)
      DEALLOCATE (t_quad, t_i_quad, t_j_quad)
      DEALLOCATE (idip, jdip, iquad, jquad)

      IF (tb%use_virial) CALL tb_add_stress(qs_env, tb, para_env)

#else
      MARK_USED(qs_env)
      MARK_USED(use_rho)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_derive_dH_off

! **************************************************************************************************
!> \brief save stress tensor
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_add_stress(qs_env, tb, para_env)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env

      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (virial)
      CALL get_qs_env(qs_env=qs_env, virial=virial)

      virial%pv_virial = virial%pv_virial - tb%sigma/para_env%num_pe

   END SUBROUTINE tb_add_stress

! **************************************************************************************************
!> \brief add contrib. to gradient
!> \param grad ...
!> \param deriv ...
!> \param dE ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE tb_add_grad(grad, deriv, dE, natom)

      REAL(KIND=dp), DIMENSION(:, :)                     :: grad
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: deriv
      REAL(KIND=dp), DIMENSION(:)                        :: dE
      INTEGER                                            :: natom

      INTEGER                                            :: i, j

      DO i = 1, natom
         DO j = 1, natom
            grad(:, i) = grad(:, i) + deriv(:, i, j)*dE(j)
         END DO
      END DO

   END SUBROUTINE tb_add_grad

! **************************************************************************************************
!> \brief add contrib. to sigma
!> \param sig ...
!> \param deriv ...
!> \param dE ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE tb_add_sig(sig, deriv, dE, natom)

      REAL(KIND=dp), DIMENSION(:, :)                     :: sig
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: deriv
      REAL(KIND=dp), DIMENSION(:)                        :: dE
      INTEGER                                            :: natom

      INTEGER                                            :: i, j

      DO i = 1, 3
         DO j = 1, natom
            sig(:, i) = sig(:, i) + deriv(:, i, j)*dE(j)
         END DO
      END DO

   END SUBROUTINE tb_add_sig

END MODULE tblite_interface

