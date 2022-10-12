!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The module to read/write QCSchema JSON files for interfacing CP2K with other programs
!> \par History
!>      10.2022 created [SB]
!> \author SB
! **************************************************************************************************
MODULE qcschema

   USE atomic_kind_types,              ONLY: get_atomic_kind
   USE basis_set_types,                ONLY: gto_basis_set_type
   USE cp2k_info,                      ONLY: cp2k_version
   USE cp_control_types,               ONLY: dft_control_type
   USE cp_fm_types,                    ONLY: cp_fm_type, cp_fm_get_element
   USE cp_log_handling,                ONLY: cp_get_default_logger, &
                                             cp_logger_get_default_io_unit, &
                                             cp_logger_type
   USE dbcsr_api,                      ONLY: dbcsr_csr_type
   USE input_section_types,            ONLY: section_vals_get, &
                                             section_vals_get_subs_vals, &
                                             section_vals_type
   USE kinds,                          ONLY: dp, default_string_length, &
                                             default_path_length, int_8, int_4
   USE mp2_types,                      ONLY: mp2_type
   USE message_passing,                ONLY: mp_comm_type
   USE particle_types,                 ONLY: particle_type
   USE periodic_table,                 ONLY: get_ptable_info
   USE qs_active_space_types,          ONLY: active_space_type, csr_idx_to_combined, &
                                             eri_type, get_irange_csr, csr_idx_from_combined
   USE qs_force_types,                 ONLY: qs_force_type
   USE qs_ks_types,                    ONLY: qs_ks_env_type
   USE qs_kind_types,                  ONLY: qs_kind_type, get_qs_kind
   USE qs_mo_types,                    ONLY: mo_set_p_type, get_mo_set
   USE qs_scf_types,                   ONLY: qs_scf_env_type
   USE qs_energy_types,                ONLY: qs_energy_type
   USE qs_environment_types,           ONLY: get_qs_env, &
                                             qs_environment_type
   USE hdf5,                           ONLY: h5open_f, h5close_f, &
                                             h5fcreate_f, h5fopen_f, h5fclose_f, &
                                             h5screate_f, h5screate_simple_f, h5sclose_f, &
                                             h5sget_simple_extent_npoints_f, H5S_SCALAR_F, &
                                             h5gcreate_f, h5gopen_f, h5gclose_f, &
                                             h5acreate_f, h5awrite_f, h5aread_f, h5aopen_f, h5aclose_f, &
                                             h5dcreate_f, h5dopen_f, h5dwrite_f, h5dread_f, h5dclose_f, &
                                             h5dget_space_f, &
                                             h5tset_size_f, h5tset_cset_f, h5tcopy_f, h5tset_strpad_f, &
                                             h5tcreate_f, h5tenum_insert_f, h5tclose_f, &
                                             H5T_STRING, H5T_CSET_UTF8_F, H5T_ENUM_F, &
                                             H5T_STR_NULLPAD_F, H5T_NATIVE_INTEGER, &
                                             H5T_NATIVE_DOUBLE, HID_T, H5F_ACC_TRUNC_F, &
                                             SIZE_T, H5T_C_S1, HSIZE_T, H5F_ACC_RDONLY_F
   USE iso_c_binding,                  ONLY: c_ptr, C_LOC, c_int, C_F_POINTER

#include "./base/base_uses.f90"


   IMPLICIT NONE

   PRIVATE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'qcschema'

   PUBLIC :: qcschema_type
   PUBLIC :: qcschema_env_create, qcschema_env_release, qcschema_to_hdf5, &
             read_pmat_from_hdf5, qcschema_update_fock, read_return_energy_from_hdf5


! **************************************************************************************************
!> \brief A derived type to store the program information that generated the QCSchema file.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#provenance
! **************************************************************************************************
   TYPE qcschema_provenance
      CHARACTER(LEN=default_string_length) :: creator ! The name of the creator of this object
      CHARACTER(LEN=default_string_length) :: version ! The version of the creator of this object
      CHARACTER(LEN=default_string_length) :: routine ! The routine that was used to create this object
   END TYPE qcschema_provenance


! **************************************************************************************************
!> \brief A derived type to store the topological information of the physical system.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#topology
! **************************************************************************************************
   TYPE qcschema_topology
      CHARACTER(LEN=default_string_length)         :: name        ! of the molecule
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE  :: symbols     ! of the atoms
      REAL(KIND=dp),  DIMENSION(:), ALLOCATABLE    :: geometry    ! row major, in bohr
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: masses
      INTEGER, DIMENSION(:), ALLOCATABLE           :: atomic_numbers
      INTEGER                                      :: molecular_charge
      INTEGER                                      :: molecular_multiplicity
      CHARACTER(LEN=default_string_length)         :: schema_name
      INTEGER                                      :: schema_version
      TYPE(qcschema_provenance)                    :: provenance
   END TYPE qcschema_topology


! **************************************************************************************************
!> \brief A derived type to store the information of a single electron shell in a basis set.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L43
! **************************************************************************************************
   TYPE qcschema_electron_shell
      ! The angular momenta of this electron shell as a list of integers
      INTEGER, DIMENSION(:), POINTER               :: angular_momentum => Null()
      ! The type of this shell: spherical or cartesian
      CHARACTER(LEN=9)                             :: harmonic_type
      ! The exponents of this contracted shell. The official spec stores these values as strings
      REAL(KIND=dp), DIMENSION(:), POINTER         :: exponents => Null()
      ! The general contraction coefficients of this contracted shell
      REAL(KIND=dp), DIMENSION(:,:), POINTER       :: coefficients => Null()
   END TYPE qcschema_electron_shell


! **************************************************************************************************
!> \brief A derived type to store the information of an ECP in a basis set.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L90
! **************************************************************************************************
   TYPE qcschema_ecp
      ! The type of this potential
      CHARACTER(LEN=default_string_length)         :: ecp_type
      ! The angular momenta of this potential as a list of integers
      INTEGER, DIMENSION(:), POINTER               :: angular_momentum => Null()
      ! The exponents of the r terms
      INTEGER, DIMENSION(:), POINTER               :: r_exponents => Null()
      ! The exponents of the Gaussian terms
      REAL(KIND=dp), DIMENSION(:), POINTER         :: gaussian_exponents => Null()
      ! The general contraction coefficients of this potential
      REAL(KIND=dp), DIMENSION(:,:), POINTER       :: coefficients => Null()
   END TYPE qcschema_ecp


! **************************************************************************************************
!> \brief A derived type to store the information of a single atom/center in the basis.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L146
! **************************************************************************************************
   TYPE qcschema_center_basis
      ! The list of electronic shells for this element
      TYPE(qcschema_electron_shell), DIMENSION(:), POINTER :: electron_shells => Null()
      ! The list of effective core potentials for this element
      TYPE(qcschema_ecp), DIMENSION(:), POINTER            :: ecp_potentials => Null()
      ! The number of electrons replaced by an ECP
      INTEGER                                              :: ecp_electrons
   END TYPE qcschema_center_basis


! **************************************************************************************************
!> \brief A derived type to store the information of the basis set used in the calculation.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_basis.html#basis-set-schema
! **************************************************************************************************
   TYPE qcschema_basis_set
      ! The name of the basis set
      CHARACTER(LEN=default_string_length) :: name
      ! A dictionary mapping the keys provided by `atom_map` to their basis center data
      TYPE(qcschema_center_basis), DIMENSION(:), POINTER :: center_data => Null()
      ! The list of atomic kinds, indicating the keys used to store the basis in `center_data`
      ! Not clear if this will be of the length of the basis set size, or rather just one
      ! entry for atomic kind. E.g. only one entry for hydrogen even though there might be
      ! many hydrogen atoms in the molecule. If this is the case, then we really need a
      ! hash table for `center_data`
      CHARACTER(LEN=2), DIMENSION(:), POINTER            :: atom_map => Null()
      ! The version of this specific schema
      INTEGER                                            :: schema_version
      ! The name of this schema. This value is expected to be `qcschema_basis`
      CHARACTER(LEN=default_string_length)               :: schema_name
      ! A description of this basis set
      CHARACTER(LEN=default_string_length)               :: description
   END TYPE qcschema_basis_set


! **************************************************************************************************
!> \brief A derived type to store any additional computed wavefunction properties.
!>        Matrix quantities are stored as flat, column-major arrays.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_wf.html#wavefunction-schema
! **************************************************************************************************
   TYPE qcschema_wavefunction

      ! The name of the method used to obtain the wf
      CHARACTER(LEN=default_string_length)         :: method

      ! The basis set used during the computation
      TYPE(qcschema_basis_set)                     :: basis_set

      ! SCF quantities in AO or MO basis
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_orbitals_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_orbitals_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eigenvalues_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eigenvalues_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_occupations_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_occupations_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_density_mo_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_density_mo_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_fock_mo_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_fock_mo_b

      ! Electron repulsion integrals
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri_mo_aa
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri_mo_ab
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: scf_eri_mo_bb

      ! Quantities with localized orbitals. All `nmo` orbitals are included,
      ! even if only a subset were localized
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: localized_orbitals_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: localized_orbitals_b
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: localized_fock_a
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: localized_fock_b

      ! Whether the computation used restricted spin orbitals
      LOGICAL :: restricted

   END TYPE qcschema_wavefunction


! **************************************************************************************************
!> \brief A derived type to store the computed properties of the original calculation.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_props.html#properties-schema
! **************************************************************************************************
   TYPE qcschema_properties

      REAL(KIND=dp) :: return_energy = 0.0_dp

      INTEGER :: calcinfo_nbasis = 0 ! AO basis size
      INTEGER :: calcinfo_nmo    = 0 ! MO basis size
      INTEGER :: calcinfo_nalpha = 0 ! # of alpha electrons
      INTEGER :: calcinfo_nbeta  = 0 ! # of beta electrons
      INTEGER :: calcinfo_natom  = 0

      ! SCF results
      INTEGER :: scf_iterations = 0
      REAL(KIND=dp) :: scf_one_electron_energy          = 0.0_dp
      REAL(KIND=dp) :: scf_two_electron_energy          = 0.0_dp
      REAL(KIND=dp) :: nuclear_repulsion_energy         = 0.0_dp
      REAL(KIND=dp) :: scf_vv10_energy                  = 0.0_dp
      REAL(KIND=dp) :: scf_xc_energy                    = 0.0_dp
      REAL(KIND=dp) :: scf_dispersion_correction_energy = 0.0_dp
      REAL(KIND=dp) :: scf_total_energy                 = 0.0_dp
      ! the dipole moment is calculated on the fly and not stored
      REAL(KIND=dp), DIMENSION(3) :: scf_dipole_moment  = 0.0_dp

      ! MP2 results
      REAL(KIND=dp) :: mp2_same_spin_correlation_energy     = 0.0_dp
      REAL(KIND=dp) :: mp2_opposite_spin_correlation_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_singles_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_doubles_energy = 0.0_dp
      ! these are the only two that are saved
      REAL(KIND=dp) :: mp2_correlation_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_total_energy       = 0.0_dp

      ! internal flags to know the type of calculation
      LOGICAL :: mp2 = .FALSE.

   END TYPE qcschema_properties


! **************************************************************************************************
!> \brief The full QCSchema output type.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#output-components
! **************************************************************************************************
   TYPE qcschema_type
      TYPE(qcschema_topology)     :: topology
      TYPE(qcschema_provenance)   :: provenance
      TYPE(qcschema_properties)   :: properties
      TYPE(qcschema_wavefunction) :: wavefunction
      TYPE(qcschema_basis_set)    :: basis
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: return_result
      CHARACTER(LEN=default_string_length) :: driver
      LOGICAL :: success
   END TYPE qcschema_type


   CONTAINS


! **************************************************************************************************
!> \brief Create and initialize a qcschema object from a quickstep environment
!> \param qcschema_env the qcschema environment to populate
!> \param qs_env the qs environment with all the info of the computation
! **************************************************************************************************
   SUBROUTINE qcschema_env_create(qcschema_env, qs_env)
      TYPE(qcschema_type), INTENT(INOUT)                 :: qcschema_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input, hfx_sections
      INTEGER :: handle, output_unit, natoms, nspins, nmo, nao, nel, &
                 nalpha, nbeta, i_glb_start, i_glb_end, &
                 atomic_number, iatom, ikind
      CHARACTER(LEN=2) :: atomic_symbol
      CHARACTER(LEN=default_string_length) :: method, basis_set_name
      REAL(KIND=dp)    :: dispersion, one_el_en, two_el_en, mass
      LOGICAL          :: do_hfx
      CHARACTER(LEN=*), PARAMETER :: routineN = 'qcschema_env_create'

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      ! reset everything
      CALL qcschema_env_release(qcschema_env)

      ! collect environment info
      IF (ASSOCIATED(qs_env)) THEN
         CALL get_qs_env(qs_env, ks_env=ks_env, energy=energy, &
                         dft_control=dft_control, force=force, &
                         particle_set=particle_set, &
                         scf_env=scf_env, mp2_env=mp2_env, &
                         input=input, qs_kind_set=kind_set, &
                         active_space=active_space_env)
      ELSE
         CPABORT("QS environment not associated, QCSchema interface quitting")
      END IF

      ! we need the AS environemnt to get all the SCF data
      IF (.NOT. ASSOCIATED(active_space_env)) THEN
         CPABORT("Active space environment not associated, QCSchema interface quitting")
      END IF

      !========================================================================================!
      ! *** QCSchema provenance ***
      !========================================================================================!

      qcschema_env%provenance%creator = 'CP2K'
      qcschema_env%provenance%version = cp2k_version
      qcschema_env%provenance%routine = routineN


      !========================================================================================!
      ! *** QCSchema topology ***
      !========================================================================================!

      qcschema_env%topology%schema_name = 'qcschema'
      qcschema_env%topology%schema_version = 3

      natoms = SIZE(particle_set)

      ALLOCATE (qcschema_env%topology%geometry(3*natoms))
      ALLOCATE (qcschema_env%topology%symbols(natoms))
      ALLOCATE (qcschema_env%topology%atomic_numbers(natoms))
      ALLOCATE (qcschema_env%topology%masses(natoms))

      DO iatom = 1, natoms
         ! set the geometry as a flat array
         qcschema_env%topology%geometry((iatom-1)*3+1:(iatom)*3) = particle_set(iatom)%r(1:3)

         ! set the atomic symbols
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=atomic_symbol)
         qcschema_env%topology%symbols(iatom) = atomic_symbol

         ! set the atomic numbers and masses
         CALL get_ptable_info(atomic_symbol, number=atomic_number, amass=mass)
         qcschema_env%topology%atomic_numbers(iatom) = atomic_number
         qcschema_env%topology%masses(iatom) = mass
      END DO

      qcschema_env%topology%molecular_charge = dft_control%charge
      qcschema_env%topology%molecular_multiplicity = dft_control%multiplicity


      !========================================================================================!
      ! *** QCSchema properties ***
      !========================================================================================!

      nspins = active_space_env%nspins

      nao = active_space_env%mos_active(1)%nao
      nmo = SIZE(active_space_env%active_orbitals, dim=1)
      nel = active_space_env%nelec_active

      ! TODO: I think this info should be on the AS env
      ! IF (MOD(nel,2) == 1) THEN
      !    nalpha = (nel+1)/2
      !    nbeta  = (nel-1)/2
      ! ELSE
      !    nalpha = nel/2
      !    nbeta  = nel/2
      ! END IF
      if (nspins == 1) then
         nalpha = active_space_env%nelec_active_spinwise(1)/2
         nbeta  = nalpha
      else
         nalpha = active_space_env%nelec_active_spinwise(1)
         nbeta  = active_space_env%nelec_active_spinwise(2)
      end if

      qcschema_env%properties%calcinfo_natom  = natoms
      qcschema_env%properties%calcinfo_nbasis = nao
      qcschema_env%properties%calcinfo_nmo    = nmo
      qcschema_env%properties%calcinfo_nalpha = nalpha
      qcschema_env%properties%calcinfo_nbeta  = nbeta

      ! energy results
      qcschema_env%properties%return_energy = energy%total
      qcschema_env%properties%scf_total_energy = energy%total
      ! TODO: really get the nuclear repulsion...
      qcschema_env%properties%nuclear_repulsion_energy = active_space_env%energy_inactive
      ! SCF info
      qcschema_env%properties%scf_iterations = scf_env%iter_count
      ! one-electron energy is the sum of all core terms
      one_el_en = energy%core_overlap + energy%core_self + energy%core
      qcschema_env%properties%scf_two_electron_energy = one_el_en
      ! two-electron energy is the sum of hartree and exact exchange (if there)
      two_el_en = energy%hartree + energy%ex + energy%hartree_1c
      qcschema_env%properties%scf_one_electron_energy = two_el_en
      ! xc energy
      qcschema_env%properties%scf_xc_energy = &
      energy%exc + energy%exc_aux_fit + energy%exc1 + energy%exc1_aux_fit
      ! dispersion energy
      dispersion = energy%dispersion + energy%gcp
      qcschema_env%properties%scf_dispersion_correction_energy = dispersion

      ! Some methods of CP2K are not supported by QCSchema, let's warn the user
      IF (dft_control%smear) CPABORT('WARNING: smearing not supported in QCSchema')
      IF (dft_control%dft_plus_u) CPABORT('WARNING: DFT+U not supported in QCSchema')
      IF (dft_control%do_sccs) CPABORT('WARNING: SCCS not supported in QCSchema')
      IF (qs_env%qmmm) CPABORT('WARNING: QM/MM not supported in QCSchema')
      IF (dft_control%qs_control%mulliken_restraint) &
         CPABORT('WARNING: Mulliken restrains not supported in QCSchema')
      IF (dft_control%qs_control%semi_empirical) &
         CPABORT('WARNING: semi_empirical methods not supported in QCSchema')
      IF (dft_control%qs_control%dftb) CPABORT('WARNING: DFTB not supported in QCSchema')
      IF (dft_control%qs_control%xtb) CPABORT('WARNING: xTB not supported in QCSchema')

      ! MP2 info
      IF (ASSOCIATED(qs_env%mp2_env)) THEN
         qcschema_env%properties%mp2 = .TRUE.
         ! this info is computed on the fly, but not stored!
         ! qcschema_env%properties%mp2_same_spin_correlation_energy
         ! qcschema_env%properties%mp2_opposite_spin_correlation_energy

         qcschema_env%properties%mp2_correlation_energy = energy%mp2
         qcschema_env%properties%mp2_total_energy = energy%total

         ! update the scf energy
         qcschema_env%properties%scf_total_energy = energy%total - energy%mp2
      END IF


      !========================================================================================!
      ! *** QCSchema wavefunction ***
      !========================================================================================!

      IF (nspins == 1) THEN
         qcschema_env%wavefunction%restricted = .TRUE.
      ELSE
         qcschema_env%wavefunction%restricted = .FALSE.
      END IF

      ! get the right indices to extract the MO information
      i_glb_start = active_space_env%active_orbitals(1, 1)
      i_glb_end   = active_space_env%active_orbitals(nmo, 1)

      ! alpha MO energies
      ALLOCATE (qcschema_env%wavefunction%scf_eigenvalues_a(nmo))
      qcschema_env%wavefunction%scf_eigenvalues_a(1:nmo) = &
         active_space_env%mos_active(1)%eigenvalues(i_glb_start:i_glb_end)

      ! alpha MO occupations
      ALLOCATE (qcschema_env%wavefunction%scf_occupations_a(nmo))
      qcschema_env%wavefunction%scf_occupations_a(1:nmo) = &
         active_space_env%mos_active(1)%occupation_numbers(i_glb_start:i_glb_end)

      ALLOCATE (qcschema_env%wavefunction%scf_fock_mo_a(nmo*nmo))
      ! alpha Fock matrix
      CALL subspace_symm_matrix_to_array(active_space_env%fock_sub(1), &
         qcschema_env%wavefunction%scf_fock_mo_a, i_glb_start, i_glb_end)

      ALLOCATE (qcschema_env%wavefunction%scf_density_mo_a(nmo*nmo))
      ! alpha density matrix
      CALL subspace_symm_matrix_to_array(active_space_env%p_active(1), &
         qcschema_env%wavefunction%scf_density_mo_a, i_glb_start, i_glb_end)

      ALLOCATE (qcschema_env%wavefunction%scf_orbitals_a(nao*nmo))
      ! alpha MOs coefficients
      CALL subspace_full_matrix_to_array(active_space_env%mos_active(1)%mo_coeff, &
         qcschema_env%wavefunction%scf_orbitals_a, 1, i_glb_start, nao, nmo)

      IF (nspins == 2) THEN
         i_glb_start = active_space_env%active_orbitals(1, 2)
         i_glb_end   = active_space_env%active_orbitals(nmo, 2)

         ! beta MO energies
         ALLOCATE (qcschema_env%wavefunction%scf_eigenvalues_b(nmo))
         qcschema_env%wavefunction%scf_eigenvalues_b(1:nmo) = &
            active_space_env%mos_active(2)%eigenvalues(i_glb_start:i_glb_end)

         ! beta MO occupations
         ALLOCATE (qcschema_env%wavefunction%scf_occupations_b(nmo))
         qcschema_env%wavefunction%scf_occupations_b(1:nmo) = &
            active_space_env%mos_active(2)%occupation_numbers(i_glb_start:i_glb_end)

         ALLOCATE (qcschema_env%wavefunction%scf_fock_mo_b(nmo*nmo))
         ! beta Fock matrix
         CALL subspace_symm_matrix_to_array(active_space_env%fock_sub(2), &
            qcschema_env%wavefunction%scf_fock_mo_b, i_glb_start, i_glb_end)

         ALLOCATE (qcschema_env%wavefunction%scf_density_mo_b(nmo*nmo))
         ! beta density matrix
         CALL subspace_symm_matrix_to_array(active_space_env%p_active(2), &
            qcschema_env%wavefunction%scf_density_mo_b, i_glb_start, i_glb_end)

         ALLOCATE (qcschema_env%wavefunction%scf_orbitals_b(nao*nmo))
         ! beta MOs coefficients
         CALL subspace_full_matrix_to_array(active_space_env%mos_active(2)%mo_coeff, &
            qcschema_env%wavefunction%scf_orbitals_b, 1, i_glb_start, nao, nmo)
      END IF

      ! get the alpha-alpha eri
      ALLOCATE (qcschema_env%wavefunction%scf_eri_mo_aa(nmo**4))
      CALL eri_to_array(active_space_env%eri, qcschema_env%wavefunction%scf_eri_mo_aa, &
                        active_space_env%active_orbitals, 1, 1)

      IF (nspins == 2) THEN
         ! get the alpha-beta eri
         ALLOCATE (qcschema_env%wavefunction%scf_eri_mo_ab(nmo**4))
         CALL eri_to_array(active_space_env%eri, qcschema_env%wavefunction%scf_eri_mo_ab, &
                           active_space_env%active_orbitals, 1, 2)

         ! get the beta-beta eri
         ALLOCATE (qcschema_env%wavefunction%scf_eri_mo_bb(nmo**4))
         CALL eri_to_array(active_space_env%eri, qcschema_env%wavefunction%scf_eri_mo_bb, &
                           active_space_env%active_orbitals, 2, 2)
      END IF

      !========================================================================================!
      ! *** QCSchema model ***
      !========================================================================================!

      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(kind_set(ikind), basis_set=basis_set)

         basis_set_name = basis_set%name

         ! make sure that we do not run a mixed basis set
         IF (iatom > 1) THEN
            CPASSERT(basis_set_name == basis_set%name)
         END IF
      END DO
      qcschema_env%wavefunction%basis_set%name = basis_set_name

      ! figure out which method was used for the calculation
      IF (dft_control%uks) THEN
         method = 'U'
      ELSE IF (dft_control%roks) THEN
         method = 'RO'
      ELSE
         method = 'R'
      END IF

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx)

      IF (do_hfx) THEN
         method = TRIM(method) // 'HF'
      ELSE IF (qcschema_env%properties%mp2) THEN
         method = TRIM(method) // 'MP2'
      ELSE
         method = TRIM(method) // 'KS'
      END IF

      qcschema_env%wavefunction%method = TRIM(method)


      !========================================================================================!
      ! *** QCSchema root ***
      !========================================================================================!

      ! driver
      IF (ASSOCIATED(force)) THEN
         qcschema_env%driver = 'gradient'
      ELSE
         qcschema_env%driver = 'energy'
      END IF

      ! success
      ! TODO: how to check if the calculation was succesful?
      qcschema_env%success = .TRUE.

      ! return result
      IF (qcschema_env%success) THEN
         IF (qcschema_env%driver == 'energy') THEN
            ALLOCATE (qcschema_env%return_result(1))
            qcschema_env%return_result(1) = energy%total
         ELSE
            ALLOCATE (qcschema_env%return_result(3*SIZE(particle_set)))
            ! TODO: populate with forces!!
            qcschema_env%return_result = 0.0_dp
         END IF
      ELSE
         CPABORT("The calculation to build the AS is unsuccessful")
      END IF


      CALL timestop(handle)

   END SUBROUTINE qcschema_env_create


! **************************************************************************************************
!> \brief Releases the allocated memory of a qcschema environment
!> \param qcschema_env the qcschema environment to release
! **************************************************************************************************
   SUBROUTINE qcschema_env_release(qcschema_env)
      TYPE(qcschema_type), INTENT(INOUT)      :: qcschema_env

      IF (ALLOCATED(qcschema_env%return_result)) THEN
         DEALLOCATE (qcschema_env%return_result)
      END IF

      IF (ALLOCATED(qcschema_env%topology%atomic_numbers)) THEN
         DEALLOCATE (qcschema_env%topology%atomic_numbers)
      END IF

      IF (ALLOCATED(qcschema_env%topology%masses)) THEN
         DEALLOCATE (qcschema_env%topology%masses)
      END IF

      IF (ALLOCATED(qcschema_env%topology%geometry)) THEN
         DEALLOCATE (qcschema_env%topology%geometry)
      END IF

      IF (ALLOCATED(qcschema_env%topology%symbols)) THEN
         DEALLOCATE (qcschema_env%topology%symbols)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_density_mo_a)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_density_mo_a)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_density_mo_b)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_density_mo_b)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_fock_mo_a)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_fock_mo_a)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_fock_mo_b)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_fock_mo_b)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_orbitals_a)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_orbitals_a)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_orbitals_b)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_orbitals_b)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_eigenvalues_a)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_eigenvalues_a)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_eigenvalues_b)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_eigenvalues_b)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_occupations_a)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_occupations_a)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_occupations_b)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_occupations_b)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_eri)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_eri)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_eri_mo_aa)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_eri_mo_aa)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_eri_mo_bb)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_eri_mo_bb)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%scf_eri_mo_ab)) THEN
         DEALLOCATE (qcschema_env%wavefunction%scf_eri_mo_ab)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%localized_orbitals_a)) THEN
         DEALLOCATE (qcschema_env%wavefunction%localized_orbitals_a)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%localized_orbitals_b)) THEN
         DEALLOCATE (qcschema_env%wavefunction%localized_orbitals_b)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%localized_fock_a)) THEN
         DEALLOCATE (qcschema_env%wavefunction%localized_fock_a)
      END IF

      IF (ALLOCATED(qcschema_env%wavefunction%localized_fock_b)) THEN
         DEALLOCATE (qcschema_env%wavefunction%localized_fock_b)
      END IF

   END SUBROUTINE qcschema_env_release


! **************************************************************************************************
!> \brief Writes a qcschema object to an hdf5 file
!> \param qcschema_env the qcschema environment to write to file
! **************************************************************************************************
   SUBROUTINE qcschema_to_hdf5(qcschema_env)
      TYPE(qcschema_type), INTENT(IN)           :: qcschema_env
      TYPE(cp_logger_type), POINTER             :: logger
      INTEGER                                   :: error, output_unit
      INTEGER(KIND=int_8)                       :: nresult
      INTEGER(KIND=hid_t)                       :: file_id, group_id
      CHARACTER(LEN=default_path_length)        :: filename

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)


      ! initialize HDF5 Fortran API
      CALL h5open_f(error)
      IF (error < 0) CPABORT('ERROR: initialization of HDF5 interface failed')


      ! create qcschema hdf5 file
      filename = TRIM(logger%iter_info%project_name) // '.qcschema.hdf5'
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 file failed')


      ! !===========================================================================!
      ! *** Root group ***
      ! !===========================================================================!
      ! driver
      CALL h5awrite_fixlen_string(file_id, 'driver', TRIM(qcschema_env%driver))
      ! return result
      nresult = SIZE(qcschema_env%return_result)
      IF (SIZE(qcschema_env%return_result) == 1) THEN
         CALL h5awrite_double_scalar(file_id, 'return_result', qcschema_env%return_result(1))
      ELSE
         CALL h5awrite_double_simple(file_id, 'return_result', qcschema_env%return_result)
      END IF
      ! schema name
      CALL h5awrite_fixlen_string(file_id, 'schema_name', TRIM(qcschema_env%topology%schema_name))
      ! schema version
      CALL h5awrite_integer_scalar(file_id, 'schema_version', qcschema_env%topology%schema_version)
      ! success
      CALL h5awrite_boolean(file_id, 'success', qcschema_env%success)


      !========================================================================================!
      ! *** QCSchema provenance ***
      !========================================================================================!
      ! create the provenance group
      CALL h5gcreate_f(file_id, 'provenance', group_id, error)
      ! populate provenance
      CALL h5awrite_fixlen_string(group_id, 'creator', TRIM(qcschema_env%provenance%creator))
      CALL h5awrite_fixlen_string(group_id, 'routine', TRIM(qcschema_env%provenance%routine))
      CALL h5awrite_fixlen_string(group_id, 'version', TRIM(qcschema_env%provenance%version))
      ! close provenance group
      CALL h5gclose_f(group_id, error)


      !========================================================================================!
      ! *** QCSchema molecule ***
      !========================================================================================!
      ! create the molecule group
      CALL h5gcreate_f(file_id, 'molecule', group_id, error)
      ! populate molecule
      CALL h5awrite_double_simple(group_id, 'geometry', qcschema_env%topology%geometry)
      CALL h5awrite_integer_simple(group_id, 'atomic_numbers', qcschema_env%topology%atomic_numbers)
      CALL h5awrite_double_simple(group_id, 'masses', qcschema_env%topology%masses)
      CALL h5awrite_integer_scalar(group_id, 'molecular_charge', qcschema_env%topology%molecular_charge)
      CALL h5awrite_integer_scalar(group_id, 'molecular_multiplicity', qcschema_env%topology%molecular_multiplicity)
      CALL h5awrite_string_simple(group_id, 'symbols', qcschema_env%topology%symbols)

      CALL h5awrite_fixlen_string(group_id, 'schema_name', 'qcschema_molecule')
      CALL h5awrite_integer_scalar(group_id, 'schema_version', 2)
      ! close molecule group
      CALL h5gclose_f(group_id, error)


      !========================================================================================!
      ! *** QCSchema properties ***
      !========================================================================================!
      ! create the properties group
      CALL h5gcreate_f(file_id, 'properties', group_id, error)
      ! populate properties
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_natom',  qcschema_env%properties%calcinfo_natom )
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nbasis', qcschema_env%properties%calcinfo_nbasis)
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nmo',    qcschema_env%properties%calcinfo_nmo   )
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nalpha', qcschema_env%properties%calcinfo_nalpha)
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nbeta',  qcschema_env%properties%calcinfo_nbeta )

      ! CALL h5dwrite_double_simple(group_id, 'scf_dipole_moment', &
      !                             qcschema_env%properties%scf_dipole_moment)

      ! energies, scf, mp2, ...
      CALL h5awrite_double_scalar(group_id, 'return_energy', qcschema_env%properties%return_energy)
      CALL h5awrite_double_scalar(group_id, 'scf_total_energy', qcschema_env%properties%scf_total_energy)
      CALL h5awrite_double_scalar(group_id, 'nuclear_repulsion_energy', &
                                  qcschema_env%properties%nuclear_repulsion_energy)

      IF (qcschema_env%properties%scf_iterations /= 0) THEN
         CALL h5awrite_integer_scalar(group_id, 'scf_iterations', qcschema_env%properties%scf_iterations)
      END IF

      IF (qcschema_env%properties%scf_one_electron_energy /= 0.0_dp) THEN
         CALL h5awrite_double_scalar(group_id, 'scf_one_electron_energy', &
                                     qcschema_env%properties%scf_one_electron_energy)
      END IF

      IF (qcschema_env%properties%scf_two_electron_energy /= 0.0_dp) THEN
         CALL h5awrite_double_scalar(group_id, 'scf_two_electron_energy', &
                                     qcschema_env%properties%scf_two_electron_energy)
      END IF

      IF (qcschema_env%properties%scf_xc_energy /= 0.0_dp) THEN
         CALL h5awrite_double_scalar(group_id, 'scf_xc_energy', &
                                     qcschema_env%properties%scf_xc_energy)
      END IF

      IF (qcschema_env%properties%scf_dispersion_correction_energy /= 0.0_dp) THEN
         CALL h5awrite_double_scalar(group_id, 'scf_dispersion_correction_energy', &
                                     qcschema_env%properties%scf_dispersion_correction_energy)
      END IF

      IF (qcschema_env%properties%mp2) THEN
         CALL h5awrite_double_scalar(group_id, 'mp2_correlation_energy', &
                                     qcschema_env%properties%mp2_correlation_energy)
      END IF

      ! close properties group
      CALL h5gclose_f(group_id, error)


      !========================================================================================!
      ! *** QCSchema wavefunction ***
      !========================================================================================!
      ! create the wavefunction group
      CALL h5gcreate_f(file_id, 'wavefunction', group_id, error)

      CALL h5awrite_fixlen_string(group_id, 'basis', TRIM(qcschema_env%wavefunction%basis_set%name))

      CALL h5dwrite_double_simple(group_id, 'scf_orbitals_a', &
                                  qcschema_env%wavefunction%scf_orbitals_a)

      CALL h5dwrite_double_simple(group_id, 'scf_eigenvalues_a', &
                                  qcschema_env%wavefunction%scf_eigenvalues_a)

      CALL h5dwrite_double_simple(group_id, 'scf_occupations_a', &
                                  qcschema_env%wavefunction%scf_occupations_a)

      CALL h5dwrite_double_simple(group_id, 'scf_fock_mo_a', &
                                  qcschema_env%wavefunction%scf_fock_mo_a)

      CALL h5dwrite_double_simple(group_id, 'scf_density_mo_a', &
                                  qcschema_env%wavefunction%scf_density_mo_a)

      CALL h5dwrite_double_simple(group_id, 'scf_eri_mo_aa', &
                                  qcschema_env%wavefunction%scf_eri_mo_aa)

      IF (.NOT. qcschema_env%wavefunction%restricted) THEN
         CALL h5dwrite_double_simple(group_id, 'scf_orbitals_b', &
                                     qcschema_env%wavefunction%scf_orbitals_b)

         CALL h5dwrite_double_simple(group_id, 'scf_eigenvalues_b', &
                                     qcschema_env%wavefunction%scf_eigenvalues_b)

         CALL h5dwrite_double_simple(group_id, 'scf_occupations_b', &
                                     qcschema_env%wavefunction%scf_occupations_b)

         CALL h5dwrite_double_simple(group_id, 'scf_fock_mo_b', &
                                     qcschema_env%wavefunction%scf_fock_mo_b)

         CALL h5dwrite_double_simple(group_id, 'scf_density_mo_b', &
                                     qcschema_env%wavefunction%scf_density_mo_b)

         CALL h5dwrite_double_simple(group_id, 'scf_eri_mo_bb', &
                                     qcschema_env%wavefunction%scf_eri_mo_bb)

         CALL h5dwrite_double_simple(group_id, 'scf_eri_mo_ba', &
                                     qcschema_env%wavefunction%scf_eri_mo_ab)

      END IF

      ! close wavefunction group
      CALL h5gclose_f(group_id, error)


      !========================================================================================!
      ! *** QCSchema model ***
      !========================================================================================!
      ! create the model group
      CALL h5gcreate_f(file_id, 'model', group_id, error)
      CALL h5awrite_fixlen_string(group_id, 'basis', TRIM(qcschema_env%wavefunction%basis_set%name))
      CALL h5awrite_fixlen_string(group_id, 'method', TRIM(qcschema_env%wavefunction%method))
      ! close model group
      CALL h5gclose_f(group_id, error)

      ! create the keywords group
      CALL h5gcreate_f(file_id, 'keywords', group_id, error)
      ! close keywords group
      CALL h5gclose_f(group_id, error)


      CALL h5fclose_f(file_id, error)

      CALL h5close_f(error)

      ! WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") 'qcschema hdf5 file written!'

   END SUBROUTINE qcschema_to_hdf5


! **************************************************************************************************
!> \brief Reads the electron density from a qcschema hdf5 file
!> \param filename the path to the qcschema hdf5 file
!> \param qcschema_env the qcschema environment onto which it writes the density
! **************************************************************************************************
   SUBROUTINE read_pmat_from_hdf5(filename, qcschema_env)
      CHARACTER(LEN=default_path_length), INTENT(IN)  :: filename
      TYPE(qcschema_type), INTENT(INOUT)              :: qcschema_env
      INTEGER                                         :: error, nmo
      INTEGER(KIND=hid_t)                             :: file_id, group_id


      ! initialize HDF5 Fortran API
      CALL h5open_f(error)
      IF (error < 0) CPABORT('ERROR: initialization of HDF5 interface failed')

      ! open qcschema hdf5 file
      CALL h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)
      IF (error < 0) CPABORT('ERROR: opening of HDF5 file failed')

      ! open the wave function group
      CALL h5gopen_f(file_id, 'wavefunction', group_id, error)
      IF (error < 0) CPABORT('ERROR: opening of "wavefunction" group failed')

      ! allocate the space for the array containing the density
      nmo = qcschema_env%properties%calcinfo_nmo
      IF (.NOT. ALLOCATED (qcschema_env%wavefunction%scf_density_mo_a)) THEN
         ALLOCATE (qcschema_env%wavefunction%scf_density_mo_a(nmo*nmo))
      END IF

      ! read the alpha density
      CALL h5dread_double_simple(group_id, 'scf_density_mo_a', qcschema_env%wavefunction%scf_density_mo_a)

      IF (.NOT. qcschema_env%wavefunction%restricted) THEN
         IF (.NOT. ALLOCATED (qcschema_env%wavefunction%scf_density_mo_b)) THEN
            ALLOCATE (qcschema_env%wavefunction%scf_density_mo_b(nmo*nmo))
         END IF
         ! read the beta density
         CALL h5dread_double_simple(group_id, 'scf_density_mo_b', qcschema_env%wavefunction%scf_density_mo_b)
      END IF

      ! close everything
      CALL h5gclose_f(group_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   END SUBROUTINE read_pmat_from_hdf5


! **************************************************************************************************
!> \brief Reads the return energy from a qcschema hdf5 file
!> \param filename the path to the qcschema hdf5 file
!> \param qcschema_env the qcschema environment onto which it writes the energy
! **************************************************************************************************
   SUBROUTINE read_return_energy_from_hdf5(filename, qcschema_env)
      CHARACTER(LEN=default_path_length), INTENT(IN)  :: filename
      TYPE(qcschema_type), INTENT(INOUT)              :: qcschema_env
      INTEGER                                         :: error
      INTEGER(KIND=hid_t)                             :: file_id, group_id


      ! initialize HDF5 Fortran API
      CALL h5open_f(error)
      IF (error < 0) CPABORT('ERROR: initialization of HDF5 interface failed')

      ! open qcschema hdf5 file
      CALL h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)
      IF (error < 0) CPABORT('ERROR: opening of HDF5 file failed')

      ! open the properties group
      CALL h5gopen_f(file_id, 'properties', group_id, error)
      IF (error < 0) CPABORT('ERROR: opening of "properties" group failed')

      ! read the alpha density
      CALL h5aread_double_scalar(group_id, 'return_energy', qcschema_env%properties%return_energy)

      ! close everything
      CALL h5gclose_f(group_id, error)
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)

   END SUBROUTINE read_return_energy_from_hdf5


! **************************************************************************************************
!> \brief Updates the Fock matrix and the inactive energy in a qcschema object
!> \param qcschema_env the qcschema environment
!> \param active_space_env the active space environment with the updated data
! **************************************************************************************************
   SUBROUTINE qcschema_update_fock(qcschema_env, active_space_env)
      TYPE(qcschema_type), INTENT(INOUT)              :: qcschema_env
      TYPE(active_space_type), INTENT(IN), POINTER    :: active_space_env
      INTEGER                                         :: nmo, i_glb_start, i_glb_end

      nmo = qcschema_env%properties%calcinfo_nmo

      ! get the right indices to extract the MO information
      i_glb_start = active_space_env%active_orbitals(1, 1)
      i_glb_end   = active_space_env%active_orbitals(nmo, 1)

      ! alpha Fock matrix
      CALL subspace_symm_matrix_to_array(active_space_env%fock_sub(1), &
         qcschema_env%wavefunction%scf_fock_mo_a, i_glb_start, i_glb_end)

      ! beta Fock matrix
      IF (active_space_env%nspins == 2) THEN
         i_glb_start = active_space_env%active_orbitals(1, 2)
         i_glb_end   = active_space_env%active_orbitals(nmo, 2)

         CALL subspace_symm_matrix_to_array(active_space_env%fock_sub(2), &
            qcschema_env%wavefunction%scf_fock_mo_b, i_glb_start, i_glb_end)
      END IF

      ! update inactive energy
      qcschema_env%properties%nuclear_repulsion_energy = active_space_env%energy_inactive

   END SUBROUTINE qcschema_update_fock


! **************************************************************************************************
!> \brief Copy a (square portion) of a `cp_fm_type` matrix to a standard 1D Fortran array
!> \param source_matrix the matrix from where the data is taken
!> \param target_array the array were the data is copied to
!> \param start_idx the index from which the submatrix starts
!> \param end_idx the index at which the submatrix ends
! **************************************************************************************************
   SUBROUTINE subspace_symm_matrix_to_array(source_matrix, target_array, start_idx, end_idx)
      INTEGER, INTENT(IN)                                :: start_idx, end_idx
      TYPE(cp_fm_type), INTENT(IN)                       :: source_matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: target_array

      INTEGER                                            :: i, j, i_sub, j_sub, norb
      REAL(KIND=dp)                                      :: mval

      CPASSERT(end_idx >= start_idx)
      CPASSERT(start_idx > 0)

      norb = end_idx - start_idx + 1

      CPASSERT(SIZE(target_array) == norb*norb)

      DO j = 1, norb
         j_sub = start_idx + j - 1
         DO i = j, norb
            i_sub = start_idx + i - 1
            CALL cp_fm_get_element(source_matrix, i_sub, j_sub, mval)
            target_array(i + (j-1)*norb) = mval
            target_array(j + (i-1)*norb) = mval
         END DO
      END DO
   END SUBROUTINE subspace_symm_matrix_to_array


! **************************************************************************************************
!> \brief Copy an arbitrary portion of a `cp_fm_type` matrix to a standard 1D Fortran array
!> \param matrix the matrix from where the data is taken
!> \param array the array were the data is copied to
!> \param nrows the number of rows to copy
!> \param ncols the number of cols to copy
!> \param first_row the row index from which the submatrix starts
!> \param first_col the col index at which the submatrix starts
! **************************************************************************************************
   SUBROUTINE subspace_full_matrix_to_array(matrix, array, first_row, first_col, nrows, ncols)
      INTEGER, INTENT(IN)                                :: first_row, first_col, nrows, ncols
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: array

      INTEGER                                            :: i, j, i_sub, j_sub
      REAL(KIND=dp)                                      :: mval

      CPASSERT(nrows > 0)
      CPASSERT(ncols > 0)
      CPASSERT(first_row > 0)
      CPASSERT(first_col > 0)

      CPASSERT(SIZE(array) == nrows*ncols)

      DO i = 1, nrows
         i_sub = first_row + i - 1
         DO j = 1, ncols
            j_sub = first_col + j - 1
            CALL cp_fm_get_element(matrix, i_sub, j_sub, mval)
            array(j + (i-1)*ncols) = mval
         END DO
      END DO
   END SUBROUTINE subspace_full_matrix_to_array


! **************************************************************************************************
!> \brief Copy the eri tensor for spins isp1 and isp2 to a standard 1D Fortran array
!> \param eri the eri tensor in dbcsr format
!> \param array the 1D Fortran array where the eri are copied to
!> \param isp1 the spin of the bra
!> \param isp2 the spin of the ket
! **************************************************************************************************
   SUBROUTINE eri_to_array(eri_env, array, active_orbitals, spin1, spin2)
      TYPE(eri_type), INTENT(IN)                         :: eri_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: array
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: active_orbitals
      INTEGER, INTENT(IN)                                :: spin1, spin2
      TYPE(dbcsr_csr_type), POINTER                      :: eri
      INTEGER                                            :: i1, i2, i3, i4, &
                                                            i, j, k, l, &
                                                            ijkl, jikl, ijlk, jilk, &
                                                            klij, lkij, klji, lkji, &
                                                            i12, i12l, i34, i34l, &
                                                            nmo_active, nmo_total, &
                                                            nindex, irptr, nmo_inactive
      INTEGER, DIMENSION(2)                              :: irange
      REAL(KIND=dp)                                      :: erival
      TYPE(mp_comm_type)                                 :: mp_group

      nmo_active = SIZE(active_orbitals, 1)
      nmo_total = eri_env%norb
      nmo_inactive = nmo_total - nmo_active
      nindex = (nmo_total*(nmo_total + 1))/2
      if (spin1 == 1 .and. spin2 == 1) then
         eri => eri_env%eri(1)%csr_mat
      else if ((spin1 == 1 .and. spin2 == 2) .or. (spin1 == 2 .and. spin2 == 1)) then
         eri => eri_env%eri(2)%csr_mat
      else
         eri => eri_env%eri(3)%csr_mat
      end if

      CALL mp_group%set_handle(eri%mp_group)
      irange = get_irange_csr(nindex, mp_group)

      array = 0.0_dp

      DO i = 1, nmo_active
         i1 = active_orbitals(i, spin1)
         DO j = i, nmo_active
            i2 = active_orbitals(j, spin2)
            i12 = csr_idx_to_combined(i1, i2, nmo_total)
            IF (i12 >= irange(1) .AND. i12 <= irange(2)) THEN
               i12l = i12 - irange(1) + 1
               irptr = eri%rowptr_local(i12l) - 1
               DO i34l = 1, eri%nzerow_local(i12l)
                  i34 = eri%colind_local(irptr + i34l)
                  CALL csr_idx_from_combined(i34, nmo_total, i3, i4)
                  k = i3 - nmo_inactive
                  l = i4 - nmo_inactive
                  erival = eri%nzval_local%r_dp(irptr + i34l)

                  ! 8-fold permutational symmetry
                  ijkl = i + (j - 1)*nmo_active + (k - 1)*nmo_active**2 + (l - 1)*nmo_active**3
                  jikl = j + (i - 1)*nmo_active + (k - 1)*nmo_active**2 + (l - 1)*nmo_active**3
                  ijlk = i + (j - 1)*nmo_active + (l - 1)*nmo_active**2 + (k - 1)*nmo_active**3
                  jilk = j + (i - 1)*nmo_active + (l - 1)*nmo_active**2 + (k - 1)*nmo_active**3
                  array(ijkl) = erival
                  array(jikl) = erival
                  array(ijlk) = erival
                  array(jilk) = erival
                  IF (spin1 == spin2) THEN
                     klij = k + (l - 1)*nmo_active + (i - 1)*nmo_active**2 + (j - 1)*nmo_active**3
                     lkij = l + (k - 1)*nmo_active + (i - 1)*nmo_active**2 + (j - 1)*nmo_active**3
                     klji = k + (l - 1)*nmo_active + (j - 1)*nmo_active**2 + (i - 1)*nmo_active**3
                     lkji = l + (k - 1)*nmo_active + (j - 1)*nmo_active**2 + (i - 1)*nmo_active**3
                     array(klij) = erival
                     array(lkij) = erival
                     array(klji) = erival
                     array(lkji) = erival
                  END IF
               END DO
            END IF
         END DO
      END DO

   END SUBROUTINE eri_to_array


!===================================================================================================!
! ***** HDF5 API WRAPPERS *****
!===================================================================================================!


! **************************************************************************************************
!> \brief Write a variable-length string attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the string to write
! **************************************************************************************************
   SUBROUTINE h5awrite_varlen_string(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      CHARACTER(LEN=*), INTENT(IN), TARGET   :: attr_data
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      TYPE(c_ptr), TARGET                    :: in_between_ptr
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! create a variable-length string type
      CALL h5tcopy_f(H5T_STRING, type_id, error)
      CAll h5tset_cset_f(type_id, H5T_CSET_UTF8_F, error)
      CALL h5tset_strpad_f(type_id, H5T_STR_NULLPAD_F, error)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! weird in-between pointer needed for variable-length
      ! string to a scalar dataspace
      in_between_ptr = C_LOC(attr_data)
      ! the actual pointer to be passed
      buffer = C_LOC(in_between_ptr)

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 datatype failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_varlen_string


! **************************************************************************************************
!> \brief Write a fixed-length string attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the string to write
! **************************************************************************************************
   SUBROUTINE h5awrite_fixlen_string(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      CHARACTER(LEN=*), INTENT(IN), TARGET   :: attr_data
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! create a fixed-length string datatype
      CALL h5tcopy_f(H5T_C_S1, type_id, error)
      CAll h5tset_cset_f(type_id, H5T_CSET_UTF8_F, error)
      CALL h5tset_size_f(type_id, LEN(attr_data, size_t), error)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data)

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 datatype failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_fixlen_string


! **************************************************************************************************
!> \brief Write a boolean attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the logical to write (.true. or .false.)
! **************************************************************************************************
   SUBROUTINE h5awrite_boolean(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      LOGICAL, INTENT(IN)                    :: attr_data
      INTEGER, TARGET                        :: attr_data_to_int
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      INTEGER(KIND=size_t)                   :: enum_size = 1 ! 8-bit integers in enum bool_type
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! create the datatype
      CALL h5tcreate_f(H5T_ENUM_F, enum_size, type_id, error)
      CALL h5tenum_insert_f(type_id, "FALSE", 0, error)
      CALL h5tenum_insert_f(type_id,  "TRUE", 1, error)

      IF (attr_data) THEN
         attr_data_to_int = 1
      ELSE
         attr_data_to_int = 0
      END IF
      ! the C pointer to the actual data
      buffer = C_LOC(attr_data_to_int)

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 datatype failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_boolean


! **************************************************************************************************
!> \brief Write a (scalar) integer attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the integer to write
! **************************************************************************************************
   SUBROUTINE h5awrite_integer_scalar(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      INTEGER, INTENT(IN), TARGET            :: attr_data
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error, output_unit

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 dataspace failed'
         RETURN
      ENDIF

      ! the C pointer to the actual data
      buffer = C_LOC(attr_data)

      ! set the type of data
      type_id = H5T_NATIVE_INTEGER

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: creation of HDF5 attribute failed'
         RETURN
      ENDIF

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: writing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 attribute failed'
         RETURN
      ENDIF

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T5,A,/)") &
         ' ERROR: closing of HDF5 dataspace failed'
         RETURN
      ENDIF

   END SUBROUTINE h5awrite_integer_scalar


! **************************************************************************************************
!> \brief Write a (scalar) double precision attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the double to write
! **************************************************************************************************
   SUBROUTINE h5awrite_double_scalar(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)        :: loc_id
      CHARACTER(LEN=*), INTENT(IN)           :: attr_name
      REAL(KIND=dp), INTENT(IN), TARGET      :: attr_data
      INTEGER(KIND=hid_t)                    :: space_id, attr_id, type_id
      TYPE(c_ptr)                            :: buffer
      INTEGER                                :: error

      ! create a scalar dataspace
      CALL h5screate_f(H5S_SCALAR_F, space_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 dataspace failed')

      ! the C pointer to the actual data
      buffer = C_LOC(attr_data)

      ! set the type of data
      type_id = H5T_NATIVE_DOUBLE

      ! create the attribute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 attribute failed')

      ! write the string attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: writing of HDF5 attribute failed')

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 attribute failed')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 dataspace failed')

   END SUBROUTINE h5awrite_double_scalar


! **************************************************************************************************
!> \brief Write an array of fixed-length string attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the array of strings
! **************************************************************************************************
   SUBROUTINE h5awrite_string_simple(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                    :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                       :: attr_name
      CHARACTER(LEN=*), INTENT(IN), DIMENSION(:), TARGET :: attr_data
      INTEGER(KIND=hid_t)                                :: space_id, attr_id, type_id
      INTEGER                                            :: error
      TYPE(c_ptr)                                        :: buffer
      INTEGER(KIND=hsize_t), DIMENSION(2)                :: dims

      dims(1) = LEN(attr_data(1), kind=hsize_t) ! length of a string entry
      dims(2) = SIZE(attr_data, kind=hsize_t)   ! length of array of strings

      ! create a fixed-length string datatype
      CALL h5tcopy_f(H5T_C_S1, type_id, error)
      CAll h5tset_cset_f(type_id, H5T_CSET_UTF8_F, error)
      CALL h5tset_size_f(type_id, dims(1), error)

      ! create a simple dataspace
      CALL h5screate_simple_f(1, dims(2:2), space_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 dataspace failed')

      ! create the atrtibute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 attribute failed')

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data(1))

      ! write the string array attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: writing of HDF5 attribute failed')

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 attribute failed')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 dataspace failed')

      ! close datatype
      CALL h5tclose_f(type_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 datatype failed')

   END SUBROUTINE h5awrite_string_simple


! **************************************************************************************************
!> \brief Write an array of doubles attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the array of doubles
! **************************************************************************************************
   SUBROUTINE h5awrite_double_simple(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                 :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                    :: attr_name
      REAL(KIND=dp), INTENT(IN), DIMENSION(:), TARGET :: attr_data
      INTEGER(KIND=hid_t)                             :: space_id, attr_id, type_id
      INTEGER                                         :: error
      TYPE(c_ptr)                                     :: buffer
      INTEGER(KIND=hsize_t), DIMENSION(1)             :: dims

      dims(1) = SIZE(attr_data, kind=hsize_t)   ! length of array of strings

      ! set the type of data
      type_id = H5T_NATIVE_DOUBLE

      ! create a simple dataspace
      CALL h5screate_simple_f(1, dims, space_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 dataspace failed')

      ! create the atrtibute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 attribute failed')

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data(1))

      ! write the string array attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: writing of HDF5 attribute failed')

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 attribute failed')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 dataspace failed')

   END SUBROUTINE h5awrite_double_simple


! **************************************************************************************************
!> \brief Write an array of integers attribute
!> \param loc_id either file id or group id
!> \param attr_name the name of the attribute
!> \param attr_data the attribute data, i.e. the array of integers
! **************************************************************************************************
   SUBROUTINE h5awrite_integer_simple(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                 :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                    :: attr_name
      INTEGER, INTENT(IN), DIMENSION(:), TARGET       :: attr_data
      INTEGER(KIND=hid_t)                             :: space_id, attr_id, type_id
      INTEGER                                         :: error
      TYPE(c_ptr)                                     :: buffer
      INTEGER(KIND=hsize_t), DIMENSION(1)             :: dims

      dims(1) = SIZE(attr_data, kind=hsize_t)   ! length of array of strings

      ! set the type of data
      type_id = H5T_NATIVE_INTEGER

      ! create a simple dataspace
      CALL h5screate_simple_f(1, dims, space_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 dataspace failed')

      ! create the atrtibute
      CALL h5acreate_f(loc_id, attr_name, type_id, space_id, attr_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 attribute failed')

      ! the actual pointer to be passed
      buffer = C_LOC(attr_data(1))

      ! write the string array attribute to file
      CALL h5awrite_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: writing of HDF5 attribute failed')

      ! close attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 attribute failed')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 dataspace failed')

   END SUBROUTINE h5awrite_integer_simple


! **************************************************************************************************
!> \brief Write a dataset containing an array of doubles
!> \param loc_id either file id or group id
!> \param dset_name the name of the dataset
!> \param dset_data the dataset data, i.e. the array of doubles
! **************************************************************************************************
   SUBROUTINE h5dwrite_double_simple(loc_id, dset_name, dset_data)
      INTEGER(KIND=hid_t), INTENT(IN)                 :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                    :: dset_name
      REAL(KIND=dp), INTENT(IN), DIMENSION(:), TARGET :: dset_data
      INTEGER(KIND=hid_t)                             :: space_id, dset_id, type_id
      INTEGER                                         :: error
      TYPE(c_ptr)                                     :: buffer
      INTEGER(KIND=hsize_t), DIMENSION(1)             :: dims

      dims(1) = SIZE(dset_data, kind=hsize_t)   ! length of array

      ! set the type of data
      type_id = H5T_NATIVE_DOUBLE

      ! create a simple dataspace
      CALL h5screate_simple_f(1, dims, space_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 dataspace failed')

      ! create the dataset
      CALL h5dcreate_f(loc_id, dset_name, type_id, space_id, dset_id, error)
      IF (error < 0) CPABORT('ERROR: creation of HDF5 dataset failed')

      ! the actual pointer to be passed
      buffer = C_LOC(dset_data(1))

      ! write the string array attribute to file
      CALL h5dwrite_f(dset_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: writing of HDF5 dataset failed')

      ! close dataset
      CALL h5dclose_f(dset_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 dataset failed')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: closing of HDF5 dataspace failed')

   END SUBROUTINE h5dwrite_double_simple


! **************************************************************************************************
!> \brief Read a dataset containing an array of doubles
!> \param loc_id either file id or group id
!> \param dset_name the name of the dataset
!> \param dset_data where the read dataset data will be written
! **************************************************************************************************
   SUBROUTINE h5dread_double_simple(loc_id, dset_name, dset_data)
      INTEGER(KIND=hid_t), INTENT(IN)                       :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                          :: dset_name
      REAL(KIND=dp), INTENT(OUT), DIMENSION(:)              :: dset_data
      INTEGER(KIND=hid_t)                                   :: space_id, dset_id, type_id, npoints
      INTEGER                                               :: error
      INTEGER(KIND=hsize_t), DIMENSION(1)                   :: dims

      dims(1) = SIZE(dset_data, kind=hsize_t)   ! length of array

      ! set the type of data
      type_id = H5T_NATIVE_DOUBLE

      ! open the dataset
      CALL h5dopen_f(loc_id, dset_name, dset_id, error)
      IF (error < 0) CPABORT('ERROR: opening HDF5 dataset failed')

      ! get information on the dataspace
      CALL h5dget_space_f(dset_id, space_id, error)
      IF (error < 0) CPABORT('ERROR: fetching HDF5 dataspace info failed')

      ! get dataspace dims
      CALL h5sget_simple_extent_npoints_f(space_id, npoints, error)
      IF (error < 0) CPABORT('ERROR: fetching HDF5 dataspace dimension failed')

      ! Why doesn't this pass?
      ! CPASSERT(dims(1) == INT(npoints*npoints, kind=hsize_t))

      ! read the data
      CALL h5dread_f(dset_id, type_id, dset_data, dims, error)
      IF (error < 0) CPABORT('ERROR: reading HDF5 dataset failed')

      ! close dataset
      CALL h5dclose_f(dset_id, error)
      IF (error < 0) CPABORT('ERROR: closing HDF5 dataset failed')

      ! close dataspace
      CALL h5sclose_f(space_id, error)
      IF (error < 0) CPABORT('ERROR: closing HDF5 dataspace failed')

   END SUBROUTINE h5dread_double_simple


! **************************************************************************************************
!> \brief Read an attribute containing a scalar double
!> \param loc_id either file id or group id
!> \param dset_name the name of the dataset
!> \param dset_data where the read dataset data will be written
! **************************************************************************************************
   SUBROUTINE h5aread_double_scalar(loc_id, attr_name, attr_data)
      INTEGER(KIND=hid_t), INTENT(IN)                       :: loc_id
      CHARACTER(LEN=*), INTENT(IN)                          :: attr_name
      REAL(KIND=dp), INTENT(OUT), TARGET                    :: attr_data
      INTEGER(KIND=hid_t)                                   :: attr_id, type_id
      TYPE(c_ptr)                                           :: buffer
      INTEGER                                               :: error

      ! set the type of data
      type_id = H5T_NATIVE_DOUBLE

      ! open the attribute
      CALL h5aopen_f(loc_id, attr_name, attr_id, error)
      IF (error < 0) CPABORT('ERROR: opening HDF5 attribute failed')

      buffer = C_LOC(attr_data)
      ! read the data
      CALL h5aread_f(attr_id, type_id, buffer, error)
      IF (error < 0) CPABORT('ERROR: reading HDF5 attribute failed')

      ! close the attribute
      CALL h5aclose_f(attr_id, error)
      IF (error < 0) CPABORT('ERROR: closing HDF5 attribute failed')

   END SUBROUTINE h5aread_double_scalar


END MODULE qcschema
