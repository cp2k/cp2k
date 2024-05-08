!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The module to read/write TREX IO files for interfacing CP2K with other programs
!> \par History
!>      05.2024 created [SB]
!> \author Stefano Battaglia
! **************************************************************************************************
MODULE trexio_utils

   USE atomic_kind_types, ONLY: get_atomic_kind
   USE basis_set_types, ONLY: gto_basis_set_type
   USE cp2k_info, ONLY: cp2k_version
   USE cp_control_types, ONLY: dft_control_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE input_section_types, ONLY: section_vals_get, &
                                  section_vals_val_get, &
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
   USE trexio, ONLY: trexio_open, trexio_close, &
                     TREXIO_HDF5, TREXIO_SUCCESS, &
                     trexio_string_of_error, trexio_t, trexio_exit_code, &
                     trexio_write_metadata_code_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'trexio_utils'

   PUBLIC :: write_trexio

   CONTAINS

! **************************************************************************************************
!> \brief Write a trexio file
!> \param qs_env the qs environment with all the info of the computation
!> \param filename the name of the file to write
! **************************************************************************************************
   SUBROUTINE write_trexio(qs_env, trexio_section)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: trexio_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_trexio'

      INTEGER                                            :: handle, output_unit
      ! CHARACTER(LEN=2)                                   :: atomic_symbol
      CHARACTER(len=default_path_length)                 :: filename
      ! CHARACTER(LEN=default_string_length)               :: basis_set_name, method
                                                            ! atomic_number, handle, i, i_glb, iatom, &
                                                            ! ikind, nalpha, nao, natoms, nbeta, &
                                                            ! nel, nmo, nspins, output_unit
      INTEGER(trexio_t)                                  :: f        ! The TREXIO file handle
      INTEGER(trexio_exit_code)                          :: rc       ! TREXIO return code
      ! LOGICAL                                            :: do_hfx
      ! REAL(KIND=dp)                                      :: dispersion, mass, one_el_en, two_el_en
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      ! TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(qs_env))

      ! collect environment info
      CALL get_qs_env(qs_env, ks_env=ks_env, energy=energy, &
                      dft_control=dft_control, force=force, &
                      particle_set=particle_set, &
                      scf_env=scf_env, mp2_env=mp2_env, &
                      input=input, qs_kind_set=kind_set)


      ! get filename
      CALL section_vals_val_get(trexio_section, "FILENAME", c_val=filename)
      WRITE (output_unit, '(T3,A)') TRIM(filename)

      ! open trexio file
      f = trexio_open('h2o-TREXIO.h5', 'w', TREXIO_HDF5, rc)
      CALL trexio_error(rc)



      !========================================================================================!
      ! Metadata
      !========================================================================================!

      rc = trexio_write_metadata_code_num(f, 1)
      CALL trexio_error(rc)







      ! natoms = SIZE(particle_set)

      ! ALLOCATE (qcschema_env%topology%geometry(3*natoms))
      ! ALLOCATE (qcschema_env%topology%symbols(natoms))
      ! ALLOCATE (qcschema_env%topology%atomic_numbers(natoms))
      ! ALLOCATE (qcschema_env%topology%masses(natoms))

      ! DO iatom = 1, natoms
      !    ! set the geometry as a flat array
      !    qcschema_env%topology%geometry((iatom - 1)*3 + 1:(iatom)*3) = particle_set(iatom)%r(1:3)

      !    ! set the atomic symbols
      !    CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=atomic_symbol)
      !    qcschema_env%topology%symbols(iatom) = atomic_symbol

      !    ! set the atomic numbers and masses
      !    CALL get_ptable_info(atomic_symbol, number=atomic_number, amass=mass)
      !    qcschema_env%topology%atomic_numbers(iatom) = atomic_number
      !    qcschema_env%topology%masses(iatom) = mass
      ! END DO

      ! qcschema_env%topology%molecular_charge = dft_control%charge
      ! qcschema_env%topology%molecular_multiplicity = dft_control%multiplicity

      ! !========================================================================================!
      ! ! *** QCSchema properties ***
      ! !========================================================================================!

      ! nspins = active_space_env%nspins

      ! nao = active_space_env%mos_active(1)%nao
      ! nmo = active_space_env%nmo_active
      ! nel = active_space_env%nelec_active

      ! IF (nspins == 1) THEN
      !    nalpha = active_space_env%nelec_active/2
      !    nbeta = nalpha
      ! ELSE
      !    nalpha = (active_space_env%nelec_active + active_space_env%multiplicity - 1)/2
      !    nbeta = (active_space_env%nelec_active - active_space_env%multiplicity + 1)/2
      ! END IF

      ! qcschema_env%properties%calcinfo_natom = natoms
      ! qcschema_env%properties%calcinfo_nbasis = nao
      ! qcschema_env%properties%calcinfo_nmo = nmo
      ! qcschema_env%properties%calcinfo_nalpha = nalpha
      ! qcschema_env%properties%calcinfo_nbeta = nbeta

      ! ! energy results
      ! qcschema_env%properties%return_energy = energy%total
      ! qcschema_env%properties%scf_total_energy = energy%total
      ! ! here we abuse the nuclear repulsion energy to store the inactive energy
      ! qcschema_env%properties%nuclear_repulsion_energy = active_space_env%energy_inactive
      ! ! SCF info
      ! qcschema_env%properties%scf_iterations = scf_env%iter_count
      ! ! one-electron energy is the sum of all core terms
      ! one_el_en = energy%core_overlap + energy%core_self + energy%core
      ! qcschema_env%properties%scf_two_electron_energy = one_el_en
      ! ! two-electron energy is the sum of hartree and exact exchange (if there)
      ! two_el_en = energy%hartree + energy%ex + energy%hartree_1c
      ! qcschema_env%properties%scf_one_electron_energy = two_el_en
      ! ! xc energy
      ! qcschema_env%properties%scf_xc_energy = &
      !    energy%exc + energy%exc_aux_fit + energy%exc1 + energy%exc1_aux_fit
      ! ! dispersion energy
      ! dispersion = energy%dispersion + energy%gcp
      ! qcschema_env%properties%scf_dispersion_correction_energy = dispersion

      ! ! Some methods of CP2K are not supported by QCSchema, let's warn the user
      ! IF (dft_control%smear) CPABORT('WARNING: smearing not supported in QCSchema')
      ! IF (dft_control%dft_plus_u) CPABORT('WARNING: DFT+U not supported in QCSchema')
      ! IF (dft_control%do_sccs) CPABORT('WARNING: SCCS not supported in QCSchema')
      ! IF (qs_env%qmmm) CPABORT('WARNING: QM/MM not supported in QCSchema')
      ! IF (dft_control%qs_control%mulliken_restraint) &
      !    CPABORT('WARNING: Mulliken restrains not supported in QCSchema')
      ! IF (dft_control%qs_control%semi_empirical) &
      !    CPABORT('WARNING: semi_empirical methods not supported in QCSchema')
      ! IF (dft_control%qs_control%dftb) CPABORT('WARNING: DFTB not supported in QCSchema')
      ! IF (dft_control%qs_control%xtb) CPABORT('WARNING: xTB not supported in QCSchema')

      ! ! MP2 info
      ! IF (ASSOCIATED(qs_env%mp2_env)) THEN
      !    qcschema_env%properties%mp2 = .TRUE.
      !    ! this info is computed on the fly, but not stored!
      !    ! qcschema_env%properties%mp2_same_spin_correlation_energy
      !    ! qcschema_env%properties%mp2_opposite_spin_correlation_energy

      !    qcschema_env%properties%mp2_correlation_energy = energy%mp2
      !    qcschema_env%properties%mp2_total_energy = energy%total

      !    ! update the scf energy
      !    qcschema_env%properties%scf_total_energy = energy%total - energy%mp2
      ! END IF

      ! !========================================================================================!
      ! ! *** QCSchema wavefunction ***
      ! !========================================================================================!

      ! IF (nspins == 1) THEN
      !    qcschema_env%wavefunction%restricted = .TRUE.
      ! ELSE
      !    qcschema_env%wavefunction%restricted = .FALSE.
      ! END IF

      ! ! alpha MO energies
      ! ALLOCATE (qcschema_env%wavefunction%scf_eigenvalues_a(nmo))
      ! DO i = 1, nmo
      !    i_glb = active_space_env%active_orbitals(i, 1)
      !    qcschema_env%wavefunction%scf_eigenvalues_a(i) = &
      !       active_space_env%mos_active(1)%eigenvalues(i_glb)
      ! END DO

      ! ! alpha MO occupations
      ! ALLOCATE (qcschema_env%wavefunction%scf_occupations_a(nmo))
      ! DO i = 1, nmo
      !    i_glb = active_space_env%active_orbitals(i, 1)
      !    qcschema_env%wavefunction%scf_occupations_a(i) = &
      !       active_space_env%mos_active(1)%occupation_numbers(i_glb)
      ! END DO

      ! ! alpha Fock matrix
      ! ALLOCATE (qcschema_env%wavefunction%scf_fock_mo_a(nmo*nmo))
      ! CALL subspace_matrix_to_array(active_space_env%fock_sub(1), &
      !                               qcschema_env%wavefunction%scf_fock_mo_a, &
      !                               active_space_env%active_orbitals(:, 1), &
      !                               active_space_env%active_orbitals(:, 1))

      ! ! alpha density matrix
      ! ALLOCATE (qcschema_env%wavefunction%scf_density_mo_a(nmo*nmo))
      ! CALL subspace_matrix_to_array(active_space_env%p_active(1), &
      !                               qcschema_env%wavefunction%scf_density_mo_a, &
      !                               active_space_env%active_orbitals(:, 1), &
      !                               active_space_env%active_orbitals(:, 1))

      ! ! alpha MOs coefficients
      ! ALLOCATE (qcschema_env%wavefunction%scf_orbitals_a(nao*nmo))
      ! CALL subspace_matrix_to_array(active_space_env%mos_active(1)%mo_coeff, &
      !                               qcschema_env%wavefunction%scf_orbitals_a, &
      !                               (/(i, i=1, nao)/), active_space_env%active_orbitals(:, 1))

      ! IF (nspins == 2) THEN
      !    ! beta MO energies
      !    ALLOCATE (qcschema_env%wavefunction%scf_eigenvalues_b(nmo))
      !    DO i = 1, nmo
      !       i_glb = active_space_env%active_orbitals(i, 2)
      !       qcschema_env%wavefunction%scf_eigenvalues_b(i) = &
      !          active_space_env%mos_active(2)%eigenvalues(i_glb)
      !    END DO

      !    ! beta MO occupations
      !    ALLOCATE (qcschema_env%wavefunction%scf_occupations_b(nmo))
      !    DO i = 1, nmo
      !       i_glb = active_space_env%active_orbitals(i, 2)
      !       qcschema_env%wavefunction%scf_occupations_b(i) = &
      !          active_space_env%mos_active(2)%occupation_numbers(i_glb)
      !    END DO

      !    ! beta Fock matrix
      !    ALLOCATE (qcschema_env%wavefunction%scf_fock_mo_b(nmo*nmo))
      !    CALL subspace_matrix_to_array(active_space_env%fock_sub(2), &
      !                                  qcschema_env%wavefunction%scf_fock_mo_b, &
      !                                  active_space_env%active_orbitals(:, 2), &
      !                                  active_space_env%active_orbitals(:, 2))

      !    ! beta density matrix
      !    ALLOCATE (qcschema_env%wavefunction%scf_density_mo_b(nmo*nmo))
      !    CALL subspace_matrix_to_array(active_space_env%p_active(2), &
      !                                  qcschema_env%wavefunction%scf_density_mo_b, &
      !                                  active_space_env%active_orbitals(:, 2), &
      !                                  active_space_env%active_orbitals(:, 2))

      !    ! beta MOs coefficients
      !    ALLOCATE (qcschema_env%wavefunction%scf_orbitals_b(nao*nmo))
      !    CALL subspace_matrix_to_array(active_space_env%mos_active(2)%mo_coeff, &
      !                                  qcschema_env%wavefunction%scf_orbitals_b, &
      !                                  (/(i, i=1, nao)/), active_space_env%active_orbitals(:, 2))
      ! END IF

      ! ! get the alpha-alpha eri
      ! ALLOCATE (qcschema_env%wavefunction%scf_eri_mo_aa(nmo**4))
      ! CALL eri_to_array(active_space_env%eri, qcschema_env%wavefunction%scf_eri_mo_aa, &
      !                   active_space_env%active_orbitals, 1, 1)

      ! IF (nspins == 2) THEN
      !    ! get the alpha-beta eri
      !    ALLOCATE (qcschema_env%wavefunction%scf_eri_mo_ab(nmo**4))
      !    CALL eri_to_array(active_space_env%eri, qcschema_env%wavefunction%scf_eri_mo_ab, &
      !                      active_space_env%active_orbitals, 1, 2)

      !    ! get the beta-beta eri
      !    ALLOCATE (qcschema_env%wavefunction%scf_eri_mo_bb(nmo**4))
      !    CALL eri_to_array(active_space_env%eri, qcschema_env%wavefunction%scf_eri_mo_bb, &
      !                      active_space_env%active_orbitals, 2, 2)
      ! END IF

      ! !========================================================================================!
      ! ! *** QCSchema model ***
      ! !========================================================================================!

      ! DO iatom = 1, natoms
      !    CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
      !    CALL get_qs_kind(kind_set(ikind), basis_set=basis_set)

      !    basis_set_name = basis_set%name

      !    ! make sure that we do not run a mixed basis set
      !    IF (iatom > 1) THEN
      !       CPASSERT(basis_set_name == basis_set%name)
      !    END IF
      ! END DO
      ! qcschema_env%wavefunction%basis_set%name = basis_set_name

      ! ! figure out which method was used for the calculation
      ! IF (dft_control%uks) THEN
      !    method = 'U'
      ! ELSE IF (dft_control%roks) THEN
      !    method = 'RO'
      ! ELSE
      !    method = 'R'
      ! END IF

      ! hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      ! CALL section_vals_get(hfx_sections, explicit=do_hfx)

      ! IF (do_hfx) THEN
      !    method = TRIM(method)//'HF'
      ! ELSE IF (qcschema_env%properties%mp2) THEN
      !    method = TRIM(method)//'MP2'
      ! ELSE
      !    method = TRIM(method)//'KS'
      ! END IF

      ! qcschema_env%wavefunction%method = TRIM(method)

      ! Close the TREXIO file
      rc = trexio_close(f)
      CALL trexio_error(rc)

      CALL timestop(handle)

   END SUBROUTINE write_trexio


   SUBROUTINE trexio_error(rc)
      INTEGER(trexio_exit_code), INTENT(IN) :: rc
      CHARACTER(LEN=128) :: err_msg

      IF (rc /= TREXIO_SUCCESS) THEN
         CALL trexio_string_of_error(rc, err_msg)
         CPABORT('TREXIO Error: '//TRIM(err_msg))
      END IF

   END SUBROUTINE trexio_error

END MODULE trexio_utils
