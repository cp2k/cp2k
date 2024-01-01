!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The module to read/write QCSchema HDF5 files for interfacing CP2K with other programs
!> \par History
!>      10.2022 created [SB]
!> \author Stefano Battaglia
! **************************************************************************************************
MODULE qcschema

   USE atomic_kind_types, ONLY: get_atomic_kind
   USE basis_set_types, ONLY: gto_basis_set_type
   USE cp2k_info, ONLY: cp2k_version
   USE cp_control_types, ONLY: dft_control_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
#ifdef __HDF5
   USE hdf5_wrapper, ONLY: &
      h5aread_double_scalar, h5awrite_boolean, h5awrite_double_scalar, h5awrite_double_simple, &
      h5awrite_fixlen_string, h5awrite_integer_scalar, h5awrite_integer_simple, &
      h5awrite_string_simple, h5close, h5dread_double_simple, h5dwrite_double_simple, h5fclose, &
      h5fcreate, h5fopen, h5gclose, h5gcreate, h5gopen, h5open, hdf5_id
#endif
   USE input_section_types, ONLY: section_vals_get, &
                                  section_vals_get_subs_vals, &
                                  section_vals_type
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp, &
                    int_8
   USE mp2_types, ONLY: mp2_type
   USE particle_types, ONLY: particle_type
   USE periodic_table, ONLY: get_ptable_info
   USE qs_active_space_types, ONLY: active_space_type
   USE qs_active_space_utils, ONLY: eri_to_array, &
                                    subspace_matrix_to_array
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_force_types, ONLY: qs_force_type
   USE qs_kind_types, ONLY: get_qs_kind, &
                            qs_kind_type
   USE qs_ks_types, ONLY: qs_ks_env_type
   USE qs_scf_types, ONLY: qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qcschema'

   PUBLIC :: qcschema_type
   PUBLIC :: qcschema_env_create, qcschema_env_release, qcschema_to_hdf5

! **************************************************************************************************
!> \brief A derived type to store the program information that generated the QCSchema file.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#provenance
! **************************************************************************************************
   TYPE qcschema_provenance
      CHARACTER(LEN=default_string_length) :: creator = "" ! The name of the creator of this object
      CHARACTER(LEN=default_string_length) :: version = "" ! The version of the creator of this object
      CHARACTER(LEN=default_string_length) :: routine = "" ! The routine that was used to create this object
   END TYPE qcschema_provenance

! **************************************************************************************************
!> \brief A derived type to store the topological information of the physical system.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#topology
! **************************************************************************************************
   TYPE qcschema_topology
      CHARACTER(LEN=default_string_length)         :: name = "" ! of the molecule
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE  :: symbols   ! of the atoms
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: geometry  ! row major, in bohr
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: masses
      INTEGER, DIMENSION(:), ALLOCATABLE           :: atomic_numbers
      INTEGER                                      :: molecular_charge = 0
      INTEGER                                      :: molecular_multiplicity = 1
      CHARACTER(LEN=default_string_length)         :: schema_name = ""
      INTEGER                                      :: schema_version = 0
      TYPE(qcschema_provenance)                    :: provenance = qcschema_provenance()
   END TYPE qcschema_topology

! **************************************************************************************************
!> \brief A derived type to store the information of a single electron shell in a basis set.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L43
! **************************************************************************************************
   TYPE qcschema_electron_shell
      ! The angular momenta of this electron shell as a list of integers
      INTEGER, DIMENSION(:), POINTER               :: angular_momentum => NULL()
      ! The type of this shell: spherical or cartesian
      CHARACTER(LEN=9)                             :: harmonic_type = ""
      ! The exponents of this contracted shell. The official spec stores these values as strings
      REAL(KIND=dp), DIMENSION(:), POINTER         :: exponents => NULL()
      ! The general contraction coefficients of this contracted shell
      REAL(KIND=dp), DIMENSION(:, :), POINTER      :: coefficients => NULL()
   END TYPE qcschema_electron_shell

! **************************************************************************************************
!> \brief A derived type to store the information of an ECP in a basis set.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L90
! **************************************************************************************************
   TYPE qcschema_ecp
      ! The type of this potential
      CHARACTER(LEN=default_string_length)         :: ecp_type = ""
      ! The angular momenta of this potential as a list of integers
      INTEGER, DIMENSION(:), POINTER               :: angular_momentum => NULL()
      ! The exponents of the r terms
      INTEGER, DIMENSION(:), POINTER               :: r_exponents => NULL()
      ! The exponents of the Gaussian terms
      REAL(KIND=dp), DIMENSION(:), POINTER         :: gaussian_exponents => NULL()
      ! The general contraction coefficients of this potential
      REAL(KIND=dp), DIMENSION(:, :), POINTER      :: coefficients => NULL()
   END TYPE qcschema_ecp

! **************************************************************************************************
!> \brief A derived type to store the information of a single atom/center in the basis.
!>        For more information refer to:
!>        https://github.com/MolSSI/QCSchema/blob/1d5ff3baa5/qcschema/dev/definitions.py#L146
! **************************************************************************************************
   TYPE qcschema_center_basis
      ! The list of electronic shells for this element
      TYPE(qcschema_electron_shell), DIMENSION(:), POINTER :: electron_shells => NULL()
      ! The list of effective core potentials for this element
      TYPE(qcschema_ecp), DIMENSION(:), POINTER            :: ecp_potentials => NULL()
      ! The number of electrons replaced by an ECP
      INTEGER                                              :: ecp_electrons = 0
   END TYPE qcschema_center_basis

! **************************************************************************************************
!> \brief A derived type to store the information of the basis set used in the calculation.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_basis.html#basis-set-schema
! **************************************************************************************************
   TYPE qcschema_basis_set
      ! The name of the basis set
      CHARACTER(LEN=default_string_length) :: name = ""
      ! A dictionary mapping the keys provided by `atom_map` to their basis center data
      TYPE(qcschema_center_basis), DIMENSION(:), POINTER :: center_data => NULL()
      ! The list of atomic kinds, indicating the keys used to store the basis in `center_data`
      ! Not clear if this will be of the length of the basis set size, or rather just one
      ! entry for atomic kind. E.g. only one entry for hydrogen even though there might be
      ! many hydrogen atoms in the molecule. If this is the case, then we really need a
      ! hash table for `center_data`
      CHARACTER(LEN=2), DIMENSION(:), POINTER            :: atom_map => NULL()
      ! The version of this specific schema
      INTEGER                                            :: schema_version = -1
      ! The name of this schema. This value is expected to be `qcschema_basis`
      CHARACTER(LEN=default_string_length)               :: schema_name = ""
      ! A description of this basis set
      CHARACTER(LEN=default_string_length)               :: description = ""
   END TYPE qcschema_basis_set

! **************************************************************************************************
!> \brief A derived type to store any additional computed wavefunction properties.
!>        Matrix quantities are stored as flat, column-major arrays.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_wf.html#wavefunction-schema
! **************************************************************************************************
   TYPE qcschema_wavefunction

      ! The name of the method used to obtain the wf
      CHARACTER(LEN=default_string_length)         :: method = ""

      ! The basis set used during the computation
      TYPE(qcschema_basis_set)                     :: basis_set = qcschema_basis_set()

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
      LOGICAL :: restricted = .FALSE.

   END TYPE qcschema_wavefunction

! **************************************************************************************************
!> \brief A derived type to store the computed properties of the original calculation.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/auto_props.html#properties-schema
! **************************************************************************************************
   TYPE qcschema_properties

      REAL(KIND=dp) :: return_energy = 0.0_dp

      INTEGER :: calcinfo_nbasis = 0 ! AO basis size
      INTEGER :: calcinfo_nmo = 0 ! MO basis size
      INTEGER :: calcinfo_nalpha = 0 ! # of alpha electrons
      INTEGER :: calcinfo_nbeta = 0 ! # of beta electrons
      INTEGER :: calcinfo_natom = 0

      ! SCF results
      INTEGER :: scf_iterations = 0
      REAL(KIND=dp) :: scf_one_electron_energy = 0.0_dp
      REAL(KIND=dp) :: scf_two_electron_energy = 0.0_dp
      REAL(KIND=dp) :: nuclear_repulsion_energy = 0.0_dp
      REAL(KIND=dp) :: scf_vv10_energy = 0.0_dp
      REAL(KIND=dp) :: scf_xc_energy = 0.0_dp
      REAL(KIND=dp) :: scf_dispersion_correction_energy = 0.0_dp
      REAL(KIND=dp) :: scf_total_energy = 0.0_dp
      ! the dipole moment is calculated on the fly and not stored
      REAL(KIND=dp), DIMENSION(3) :: scf_dipole_moment = 0.0_dp

      ! MP2 results
      REAL(KIND=dp) :: mp2_same_spin_correlation_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_opposite_spin_correlation_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_singles_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_doubles_energy = 0.0_dp
      ! these are the only two that are saved
      REAL(KIND=dp) :: mp2_correlation_energy = 0.0_dp
      REAL(KIND=dp) :: mp2_total_energy = 0.0_dp

      ! internal flags to know the type of calculation
      LOGICAL :: mp2 = .FALSE.

   END TYPE qcschema_properties

! **************************************************************************************************
!> \brief The full QCSchema output type.
!>        For more information refer to:
!>        https://molssi-qc-schema.readthedocs.io/en/latest/spec_components.html#output-components
! **************************************************************************************************
   TYPE qcschema_type
      TYPE(qcschema_topology)     :: topology = qcschema_topology()
      TYPE(qcschema_provenance)   :: provenance = qcschema_provenance()
      TYPE(qcschema_properties)   :: properties = qcschema_properties()
      TYPE(qcschema_wavefunction) :: wavefunction = qcschema_wavefunction()
      TYPE(qcschema_basis_set)    :: basis = qcschema_basis_set()
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: return_result
      CHARACTER(LEN=default_string_length) :: driver = ""
      LOGICAL :: success = .FALSE.
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qcschema_env_create'

      CHARACTER(LEN=2)                                   :: atomic_symbol
      CHARACTER(LEN=default_string_length)               :: basis_set_name, method
      INTEGER                                            :: atomic_number, handle, i, i_glb, iatom, &
                                                            ikind, nalpha, nao, natoms, nbeta, &
                                                            nel, nmo, nspins, output_unit
      LOGICAL                                            :: do_hfx
      REAL(KIND=dp)                                      :: dispersion, mass, one_el_en, two_el_en
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

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
         qcschema_env%topology%geometry((iatom - 1)*3 + 1:(iatom)*3) = particle_set(iatom)%r(1:3)

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
      nmo = active_space_env%nmo_active
      nel = active_space_env%nelec_active

      IF (nspins == 1) THEN
         nalpha = active_space_env%nelec_active/2
         nbeta = nalpha
      ELSE
         nalpha = (active_space_env%nelec_active + active_space_env%multiplicity - 1)/2
         nbeta = (active_space_env%nelec_active - active_space_env%multiplicity + 1)/2
      END IF

      qcschema_env%properties%calcinfo_natom = natoms
      qcschema_env%properties%calcinfo_nbasis = nao
      qcschema_env%properties%calcinfo_nmo = nmo
      qcschema_env%properties%calcinfo_nalpha = nalpha
      qcschema_env%properties%calcinfo_nbeta = nbeta

      ! energy results
      qcschema_env%properties%return_energy = energy%total
      qcschema_env%properties%scf_total_energy = energy%total
      ! here we abuse the nuclear repulsion energy to store the inactive energy
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

      ! alpha MO energies
      ALLOCATE (qcschema_env%wavefunction%scf_eigenvalues_a(nmo))
      DO i = 1, nmo
         i_glb = active_space_env%active_orbitals(i, 1)
         qcschema_env%wavefunction%scf_eigenvalues_a(i) = &
            active_space_env%mos_active(1)%eigenvalues(i_glb)
      END DO

      ! alpha MO occupations
      ALLOCATE (qcschema_env%wavefunction%scf_occupations_a(nmo))
      DO i = 1, nmo
         i_glb = active_space_env%active_orbitals(i, 1)
         qcschema_env%wavefunction%scf_occupations_a(i) = &
            active_space_env%mos_active(1)%occupation_numbers(i_glb)
      END DO

      ! alpha Fock matrix
      ALLOCATE (qcschema_env%wavefunction%scf_fock_mo_a(nmo*nmo))
      CALL subspace_matrix_to_array(active_space_env%fock_sub(1), &
                                    qcschema_env%wavefunction%scf_fock_mo_a, &
                                    active_space_env%active_orbitals(:, 1), &
                                    active_space_env%active_orbitals(:, 1))

      ! alpha density matrix
      ALLOCATE (qcschema_env%wavefunction%scf_density_mo_a(nmo*nmo))
      CALL subspace_matrix_to_array(active_space_env%p_active(1), &
                                    qcschema_env%wavefunction%scf_density_mo_a, &
                                    active_space_env%active_orbitals(:, 1), &
                                    active_space_env%active_orbitals(:, 1))

      ! alpha MOs coefficients
      ALLOCATE (qcschema_env%wavefunction%scf_orbitals_a(nao*nmo))
      CALL subspace_matrix_to_array(active_space_env%mos_active(1)%mo_coeff, &
                                    qcschema_env%wavefunction%scf_orbitals_a, &
                                    (/(i, i=1, nao)/), active_space_env%active_orbitals(:, 1))

      IF (nspins == 2) THEN
         ! beta MO energies
         ALLOCATE (qcschema_env%wavefunction%scf_eigenvalues_b(nmo))
         DO i = 1, nmo
            i_glb = active_space_env%active_orbitals(i, 2)
            qcschema_env%wavefunction%scf_eigenvalues_b(i) = &
               active_space_env%mos_active(2)%eigenvalues(i_glb)
         END DO

         ! beta MO occupations
         ALLOCATE (qcschema_env%wavefunction%scf_occupations_b(nmo))
         DO i = 1, nmo
            i_glb = active_space_env%active_orbitals(i, 2)
            qcschema_env%wavefunction%scf_occupations_b(i) = &
               active_space_env%mos_active(2)%occupation_numbers(i_glb)
         END DO

         ! beta Fock matrix
         ALLOCATE (qcschema_env%wavefunction%scf_fock_mo_b(nmo*nmo))
         CALL subspace_matrix_to_array(active_space_env%fock_sub(2), &
                                       qcschema_env%wavefunction%scf_fock_mo_b, &
                                       active_space_env%active_orbitals(:, 2), &
                                       active_space_env%active_orbitals(:, 2))

         ! beta density matrix
         ALLOCATE (qcschema_env%wavefunction%scf_density_mo_b(nmo*nmo))
         CALL subspace_matrix_to_array(active_space_env%p_active(2), &
                                       qcschema_env%wavefunction%scf_density_mo_b, &
                                       active_space_env%active_orbitals(:, 2), &
                                       active_space_env%active_orbitals(:, 2))

         ! beta MOs coefficients
         ALLOCATE (qcschema_env%wavefunction%scf_orbitals_b(nao*nmo))
         CALL subspace_matrix_to_array(active_space_env%mos_active(2)%mo_coeff, &
                                       qcschema_env%wavefunction%scf_orbitals_b, &
                                       (/(i, i=1, nao)/), active_space_env%active_orbitals(:, 2))
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
         method = TRIM(method)//'HF'
      ELSE IF (qcschema_env%properties%mp2) THEN
         method = TRIM(method)//'MP2'
      ELSE
         method = TRIM(method)//'KS'
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
      TYPE(qcschema_type), INTENT(INOUT)                 :: qcschema_env

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
!> \brief Updates the Fock matrix and the inactive energy in a qcschema object
!> \param qcschema_env the qcschema environment
!> \param active_space_env the active space environment with the updated data
! **************************************************************************************************
   SUBROUTINE qcschema_update_fock(qcschema_env, active_space_env)
      TYPE(qcschema_type), INTENT(INOUT)                 :: qcschema_env
      TYPE(active_space_type), INTENT(IN), POINTER       :: active_space_env

      ! alpha Fock matrix
      CALL subspace_matrix_to_array(active_space_env%fock_sub(1), &
                                    qcschema_env%wavefunction%scf_fock_mo_a, &
                                    active_space_env%active_orbitals(:, 1), &
                                    active_space_env%active_orbitals(:, 1))

      ! beta Fock matrix
      IF (active_space_env%nspins == 2) THEN
         CALL subspace_matrix_to_array(active_space_env%fock_sub(2), &
                                       qcschema_env%wavefunction%scf_fock_mo_b, &
                                       active_space_env%active_orbitals(:, 2), &
                                       active_space_env%active_orbitals(:, 2))
      END IF

      ! update inactive energy
      qcschema_env%properties%nuclear_repulsion_energy = active_space_env%energy_inactive

   END SUBROUTINE qcschema_update_fock

! **************************************************************************************************
!> \brief Writes a qcschema object to an hdf5 file
!> \param qcschema_env the qcschema environment to write to file
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE qcschema_to_hdf5(qcschema_env, filename)
      TYPE(qcschema_type), INTENT(IN)                    :: qcschema_env
      CHARACTER(LEN=default_path_length), INTENT(IN)     :: filename
#ifndef __HDF5
      CPABORT("CP2K was compiled without the HDF5 library")
      MARK_USED(filename)
      MARK_USED(qcschema_env)
#else
      INTEGER                                            :: output_unit
      INTEGER(KIND=hdf5_id)                              :: file_id, group_id
      INTEGER(KIND=int_8)                                :: nresult
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      ! initialize HDF5 Fortran API
      CALL h5open()

      ! create qcschema hdf5 file
      ! filename = TRIM(logger%iter_info%project_name) // 'hdf5'
      CALL h5fcreate(TRIM(filename), file_id)

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
      CALL h5gcreate(file_id, 'provenance', group_id)
      ! populate provenance
      CALL h5awrite_fixlen_string(group_id, 'creator', TRIM(qcschema_env%provenance%creator))
      CALL h5awrite_fixlen_string(group_id, 'routine', TRIM(qcschema_env%provenance%routine))
      CALL h5awrite_fixlen_string(group_id, 'version', TRIM(qcschema_env%provenance%version))
      ! close provenance group
      CALL h5gclose(group_id)

      !========================================================================================!
      ! *** QCSchema molecule ***
      !========================================================================================!
      ! create the molecule group
      CALL h5gcreate(file_id, 'molecule', group_id)
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
      CALL h5gclose(group_id)

      !========================================================================================!
      ! *** QCSchema properties ***
      !========================================================================================!
      ! create the properties group
      CALL h5gcreate(file_id, 'properties', group_id)
      ! populate properties
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_natom', qcschema_env%properties%calcinfo_natom)
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nbasis', qcschema_env%properties%calcinfo_nbasis)
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nmo', qcschema_env%properties%calcinfo_nmo)
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nalpha', qcschema_env%properties%calcinfo_nalpha)
      CALL h5awrite_integer_scalar(group_id, 'calcinfo_nbeta', qcschema_env%properties%calcinfo_nbeta)

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
      CALL h5gclose(group_id)

      !========================================================================================!
      ! *** QCSchema wavefunction ***
      !========================================================================================!
      ! create the wavefunction group
      CALL h5gcreate(file_id, 'wavefunction', group_id)

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
      CALL h5gclose(group_id)

      !========================================================================================!
      ! *** QCSchema model ***
      !========================================================================================!
      ! create the model group
      CALL h5gcreate(file_id, 'model', group_id)
      CALL h5awrite_fixlen_string(group_id, 'basis', TRIM(qcschema_env%wavefunction%basis_set%name))
      CALL h5awrite_fixlen_string(group_id, 'method', TRIM(qcschema_env%wavefunction%method))
      ! close model group
      CALL h5gclose(group_id)

      ! create the keywords group
      CALL h5gcreate(file_id, 'keywords', group_id)
      ! close keywords group
      CALL h5gclose(group_id)

      CALL h5fclose(file_id)
      CALL h5close()
#endif

   END SUBROUTINE qcschema_to_hdf5

#ifdef __HDF5
! **************************************************************************************************
!> \brief Reads the electron density from a qcschema hdf5 file
!> \param filename the path to the qcschema hdf5 file
!> \param qcschema_env the qcschema environment onto which it writes the density
! **************************************************************************************************
   SUBROUTINE read_pmat_from_hdf5(filename, qcschema_env)
      CHARACTER(LEN=default_path_length), INTENT(IN)     :: filename
      TYPE(qcschema_type), INTENT(INOUT)                 :: qcschema_env

      INTEGER                                            :: nmo
      INTEGER(KIND=hdf5_id)                              :: file_id, group_id

      ! initialize HDF5 Fortran API
      CALL h5open()

      ! open qcschema hdf5 file
      CALL h5fopen(TRIM(filename), file_id)

      ! open the wave function group
      CALL h5gopen(file_id, 'wavefunction', group_id)

      ! allocate the space for the array containing the density
      nmo = qcschema_env%properties%calcinfo_nmo
      IF (.NOT. ALLOCATED(qcschema_env%wavefunction%scf_density_mo_a)) THEN
         ALLOCATE (qcschema_env%wavefunction%scf_density_mo_a(nmo*nmo))
      END IF

      ! read the alpha density
      CALL h5dread_double_simple(group_id, 'scf_density_mo_a', qcschema_env%wavefunction%scf_density_mo_a)

      IF (.NOT. qcschema_env%wavefunction%restricted) THEN
         IF (.NOT. ALLOCATED(qcschema_env%wavefunction%scf_density_mo_b)) THEN
            ALLOCATE (qcschema_env%wavefunction%scf_density_mo_b(nmo*nmo))
         END IF
         ! read the beta density
         CALL h5dread_double_simple(group_id, 'scf_density_mo_b', qcschema_env%wavefunction%scf_density_mo_b)
      END IF

      ! close everything
      CALL h5gclose(group_id)
      CALL h5fclose(file_id)
      CALL h5close()

   END SUBROUTINE read_pmat_from_hdf5

! **************************************************************************************************
!> \brief Reads the return energy from a qcschema hdf5 file
!> \param filename the path to the qcschema hdf5 file
!> \param qcschema_env the qcschema environment onto which it writes the energy
! **************************************************************************************************
   SUBROUTINE read_return_energy_from_hdf5(filename, qcschema_env)
      CHARACTER(LEN=default_path_length), INTENT(IN)     :: filename
      TYPE(qcschema_type), INTENT(INOUT)                 :: qcschema_env

      INTEGER(KIND=hdf5_id)                              :: file_id, group_id

      ! initialize HDF5 Fortran API
      CALL h5open()

      ! open qcschema hdf5 file
      CALL h5fopen(TRIM(filename), file_id)

      ! open the properties group
      CALL h5gopen(file_id, 'properties', group_id)

      ! read the return energy
      CALL h5aread_double_scalar(group_id, 'return_energy', qcschema_env%properties%return_energy)

      ! close everything
      CALL h5gclose(group_id)
      CALL h5fclose(file_id)
      CALL h5close()

   END SUBROUTINE read_return_energy_from_hdf5

! **************************************************************************************************
!> \brief Reads the active space energy from a qcschema file and stores it in active_space_env
!> \param active_space_env ...
!> \param qcschema_env ...
!> \author Stefano Battaglia
! **************************************************************************************************
   SUBROUTINE read_active_energy_from_hdf5(active_space_env, qcschema_env)
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(qcschema_type)                                :: qcschema_env

      CHARACTER(LEN=default_path_length)                 :: qcschema_filename

      ! File name
      qcschema_filename = active_space_env%qcschema_filename
      ! read active space energy
      CALL read_return_energy_from_hdf5(qcschema_filename, qcschema_env)

      active_space_env%energy_active = qcschema_env%properties%return_energy
      active_space_env%energy_total = active_space_env%energy_inactive + active_space_env%energy_active

   END SUBROUTINE read_active_energy_from_hdf5
#endif

END MODULE qcschema
