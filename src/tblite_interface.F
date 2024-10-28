!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief interface to tblite [Johann Pototschnig]
!> \History creation 09.2024
! **************************************************************************************************
MODULE tblite_interface

#if defined(__TBLITE)
   USE mctc_env, ONLY: error_type
   USE mctc_io, ONLY: structure_type, new
   USE mctc_io_symbols, ONLY: symbol_to_number
   USE tblite_xtb_calculator, ONLY: xtb_calculator, new_xtb_calculator
   USE tblite_xtb_gfn1, ONLY: new_gfn1_calculator
   USE tblite_xtb_gfn2, ONLY: new_gfn2_calculator
   USE tblite_xtb_ipea1, ONLY: new_ipea1_calculator
   USE tblite_container, ONLY: container_cache
   USE tblite_wavefunction_type, ONLY: wavefunction_type, new_wavefunction
   USE tblite_scf_info, ONLY: scf_info, atom_resolved, shell_resolved, orbital_resolved, &
      & not_used
   USE tblite_scf_potential, ONLY: potential_type, new_potential
   USE tblite_xtb_h0, ONLY: get_selfenergy, get_hamiltonian, get_occupation, &
      & get_hamiltonian_gradient, tb_hamiltonian
   USE tblite_basis_type, ONLY: get_cutoff
#endif

   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE cell_types, ONLY: cell_type, get_cell
   USE particle_types, ONLY: particle_type
   USE kinds, ONLY: dp, default_string_length
   USE atomic_kind_types, ONLY: get_atomic_kind
   USE input_constants, ONLY: gfn1xtb, gfn2xtb, ipea1xtb
   USE basis_set_types, ONLY: gto_basis_set_type, allocate_gto_basis_set, write_gto_basis_set, process_gto_basis
   USE tblite_types, ONLY: tblite_type, deallocate_tblite_type, allocate_tblite_type
   USE memory_utilities, ONLY: reallocate
   USE qs_kind_types, ONLY: get_qs_kind, qs_kind_type
   USE string_utilities, ONLY: integer_to_string

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_interface'

   PUBLIC :: tb_init_geometry, tb_init_ham, tb_init_wf
   PUBLIC :: tb_set_calculator, tb_get_basis

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_geometry(qs_env, tb)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_init_geometry'

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER                                            :: iatom, natom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle, ikind
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL, DIMENSION(3)                              :: lperiod

#if defined(__TBLITE)

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
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE tb_init_geometry

   ! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param do_grad ...
! **************************************************************************************************
   SUBROUTINE tb_init_wf(tb, do_grad)

      TYPE(tblite_type), POINTER                         :: tb
      LOGICAL, INTENT(in)                                :: do_grad

      REAL(KIND=dp), PARAMETER                           :: K2Eh = 3.16681156346E-06
      INTEGER                                            :: nSpin = 1 !number of spin channels

#if defined(__TBLITE)
      TYPE(scf_info)                                     :: info

      info = tb%calc%variable_info()
      IF (info%charge > shell_resolved) CPABORT("tblite: no support for orbital resolved charge")
      IF (info%dipole > atom_resolved) CPABORT("tblite: no support for shell resolved dipole moment")
      IF (info%quadrupole > atom_resolved) CPABORT("tblite: no support shell resolved quadrup[ole moment")

      CALL new_wavefunction(tb%wfn, tb%mol%nat, tb%calc%bas%nsh, tb%calc%bas%nao, nSpin, 0.0_dp*K2Eh)

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)

      !allocate quantities later required
      ALLOCATE (tb%e_hal(tb%mol%nat), tb%e_rep(tb%mol%nat), tb%e_disp(tb%mol%nat))
      ALLOCATE (tb%e_scd(tb%mol%nat), tb%e_es(tb%mol%nat))
      ALLOCATE (tb%selfenergy(tb%calc%bas%nsh))

      IF (ALLOCATED(tb%calc%ncoord)) ALLOCATE (tb%cn(tb%mol%nat))

      IF (do_grad) ALLOCATE (tb%grad(3, tb%mol%nat))

      IF (ALLOCATED(tb%grad)) THEN
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            ALLOCATE (tb%dcndr(3, tb%mol%nat, tb%mol%nat), tb%dcndL(3, 3, tb%mol%nat))
         END IF
         ALLOCATE (tb%dsedcn(tb%calc%bas%nsh))
      END IF

#else
      CPABORT("Build without TBLITE")
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
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE tb_set_calculator

! **************************************************************************************************
!> \brief ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_ham(tb)

      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE)
      TYPE(container_cache)                              :: hcache, rcache

      tb%e_hal = 0.0_dp
      tb%e_rep = 0.0_dp
      tb%e_hal = 0.0_dp
      IF (ALLOCATED(tb%grad)) THEN
         tb%grad = 0.0_dp
         tb%sigma = 0.0_dp
      END IF

      IF (ALLOCATED(tb%calc%halogen)) THEN
         CALL tb%calc%halogen%update(tb%mol, hcache)
         IF (ALLOCATED(tb%grad)) THEN
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal, &
            & tb%grad, tb%sigma)
         ELSE
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%repulsion)) THEN
         CALL tb%calc%repulsion%update(tb%mol, rcache)
         IF (ALLOCATED(tb%grad)) THEN
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep, &
            & tb%grad, tb%sigma)
         ELSE
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%dispersion)) THEN
         CALL tb%calc%dispersion%update(tb%mol, tb%dcache)
         IF (ALLOCATED(tb%grad)) THEN
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp, &
            & tb%grad, tb%sigma)
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
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE tb_init_ham

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param energies ...
! **************************************************************************************************
   SUBROUTINE get_energy(tb, energies)

      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: energies

#if defined(__TBLITE)

      energies(:) = tb%e_hal + tb%e_rep + tb%e_disp + tb%e_scd + tb%e_es

#else
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE get_energy

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param ham ...
! **************************************************************************************************
   SUBROUTINE tb_get_hamiltonian(tb, ham)

      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), INTENT(OUT)                         :: ham(:)

#if defined(__TBLITE)
      ham = 0.0_dp

#else
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE tb_get_hamiltonian

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param gto_basis_set ...
!> \param element_symbol ...
!> \param r_cutoff ...
!> \param occ ...
! **************************************************************************************************
   SUBROUTINE tb_get_basis(tb, gto_basis_set, element_symbol, r_cutoff, occ)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis_set
      CHARACTER(len=2), INTENT(IN)                       :: element_symbol
      REAL(dp), INTENT(out)                              :: r_cutoff
      INTEGER, DIMENSION(5), INTENT(out)                 :: occ

      CHARACTER(LEN=default_string_length)               :: sng
      INTEGER                                            :: ang, i_type, id_atom, iId, ipgf, iSH, &
                                                            ishell, maxl, mprim, nset, nshell
      LOGICAL                                            :: do_ortho

      CALL allocate_gto_basis_set(gto_basis_set)

#if defined(__TBLITE)

      !identifying element in the bas data
      CALL symbol_to_number(i_type, element_symbol)
      WRITE (*, *) "element: ", i_type, element_symbol
      DO id_atom = 1, tb%mol%nat
         IF (i_type == tb%el_num(id_atom)) EXIT
      END DO

      !getting size information
      nset = tb%calc%bas%nsh_at(id_atom)
      nshell = 1
      mprim = 0
      DO ishell = 1, nset
         mprim = MAX(mprim, tb%calc%bas%cgto(ishell, id_atom)%nprim)
      END DO

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

      maxl = 0
      DO ishell = 1, nset
         ang = tb%calc%bas%cgto(ishell, id_atom)%ang
         maxl = MAX(ang, maxl)
         gto_basis_set%lmax(ishell) = ang
         gto_basis_set%lmin(ishell) = ang
         gto_basis_set%npgf(ishell) = tb%calc%bas%cgto(ishell, id_atom)%nprim
         gto_basis_set%nshell(ishell) = nshell
         gto_basis_set%n(1, ishell) = ang + 1
         gto_basis_set%l(1, ishell) = ang
         DO ipgf = 1, gto_basis_set%npgf(ishell)
            gto_basis_set%gcc(ipgf, 1, ishell) = tb%calc%bas%cgto(ishell, id_atom)%coeff(ipgf)
            gto_basis_set%zet(ipgf, ishell) = tb%calc%bas%cgto(ishell, id_atom)%alpha(ipgf)
         END DO
      END DO

      do_ortho = .FALSE.
      CALL process_gto_basis(gto_basis_set, do_ortho, nset, maxl)

      !getting cutoff
      r_cutoff = get_cutoff(tb%calc%bas)

      !getting occupation
      occ = 0
      iId = tb%mol%id(id_atom)
      DO iSh = 1, MIN(tb%calc%bas%nsh_at(id_atom), 5)
         occ(iSh) = NINT(tb%calc%h0%refocc(iSh, iId))
      END DO

      !set normalization process
      gto_basis_set%norm_type = 3

#else
      CPABORT("Build without TBLITE")
#endif

   END SUBROUTINE tb_get_basis

END MODULE tblite_interface

