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

#ifdef __TREXIO
   USE atomic_kind_types, ONLY: get_atomic_kind
   USE basis_set_types, ONLY: gto_basis_set_type
   USE cell_types, ONLY: cell_type
   USE cp2k_info, ONLY: cp2k_version
   USE cp_control_types, ONLY: dft_control_type
   USE cp_files, ONLY: close_file, file_exists, open_file
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
   USE kpoint_types, ONLY: kpoint_type, get_kpoint_info
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
                     trexio_write_metadata_code, trexio_write_metadata_code_num, &
                     trexio_write_nucleus_coord, trexio_write_nucleus_num, &
                     trexio_write_nucleus_charge, trexio_write_nucleus_label, &
                     trexio_write_nucleus_repulsion, trexio_write_cell_a, &
                     trexio_write_cell_b, trexio_write_cell_c, &
                     trexio_write_cell_g_a, trexio_write_cell_g_b, &
                     trexio_write_cell_g_c, trexio_write_cell_two_pi, &
                     trexio_write_pbc_periodic, trexio_write_pbc_k_point_num, &
                     trexio_write_pbc_k_point, trexio_write_pbc_k_point_weight, &
                     trexio_write_electron_num, trexio_write_electron_up_num, &
                     trexio_write_electron_dn_num
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'trexio_utils'

   PUBLIC :: write_trexio

   CONTAINS

#ifdef __TREXIO
! **************************************************************************************************
!> \brief Write a trexio file
!> \param qs_env the qs environment with all the info of the computation
!> \param filename the name of the file to write
! **************************************************************************************************
   SUBROUTINE write_trexio(qs_env, trexio_section)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: trexio_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_trexio'

      INTEGER                                            :: handle, output_unit, unit_trexio
      ! CHARACTER(LEN=2)                                   :: atomic_symbol
      CHARACTER(len=default_path_length)                 :: filename
      ! CHARACTER(LEN=default_string_length)               :: basis_set_name, method
                                                            ! atomic_number, handle, i, i_glb, iatom, &
                                                            ! ikind, nalpha, nao, natoms, nbeta, &
                                                            ! nel, nmo, nspins, output_unit
      INTEGER(trexio_t)                                  :: f        ! The TREXIO file handle
      INTEGER(trexio_exit_code)                          :: rc       ! TREXIO return code
      LOGICAL                                            :: explicit, do_kpoints
      REAL(KIND=dp)                                      :: e_nn
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      ! TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: input

      INTEGER                                            :: iatom, natoms, z, periodic, nkp, nel_tot, nspins
      INTEGER, DIMENSION(2)                              :: nel_spin
      INTEGER, DIMENSION(3)                              :: nkp_grid
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: coord  ! row major, in bohr
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: charge ! of the atoms
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE        :: label  ! of the elements
      CHARACTER(LEN=2)                                   :: atomic_symbol

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(qs_env))

      ! collect environment info
      CALL get_qs_env(qs_env, ks_env=ks_env, energy=energy, &
                      dft_control=dft_control, force=force, &
                      particle_set=particle_set, cell=cell, &
                      scf_env=scf_env, mp2_env=mp2_env, &
                      do_kpoints=do_kpoints, kpoints=kpoints, &
                      input=input, qs_kind_set=kind_set, &
                      nelectron_total=nel_tot, nelectron_spin=nel_spin)

      ! get filename
      CALL section_vals_val_get(trexio_section, "FILENAME", c_val=filename, explicit=explicit)
      IF (.NOT. explicit) THEN
         filename = TRIM(logger%iter_info%project_name)//'-TREXIO.h5'
      ELSE
         filename = TRIM(filename)//'.h5'
      END IF

      ! inquire whether a file with the same name already exists, if yes, delete it
      IF (file_exists(filename)) THEN
         CALL open_file(filename, unit_number=unit_trexio)
         CALL close_file(unit_number=unit_trexio, file_status="DELETE")
      END IF

      ! open trexio file
      f = trexio_open(filename, 'w', TREXIO_HDF5, rc)
      CALL trexio_error(rc)

      !========================================================================================!
      ! Metadata group
      !========================================================================================!

      rc = trexio_write_metadata_code_num(f, 1)
      CALL trexio_error(rc)

      rc = trexio_write_metadata_code(f, cp2k_version, LEN_TRIM(cp2k_version)+1)
      CALL trexio_error(rc)

      !========================================================================================!
      ! Nucleus group
      !========================================================================================!

      natoms = SIZE(particle_set)
      rc = trexio_write_nucleus_num(f, natoms)
      CALL trexio_error(rc)

      ALLOCATE (coord(3,natoms))
      ALLOCATE (label(natoms))
      ALLOCATE (charge(natoms))
      DO iatom = 1, natoms
         coord(iatom,:) = particle_set(iatom)%r(1:3)

         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=atomic_symbol, z=z)
         label(iatom) = atomic_symbol
         charge(iatom) = REAL(z, kind=dp)
      END DO

      rc = trexio_write_nucleus_coord(f, coord)
      CALL trexio_error(rc)
      DEALLOCATE (coord)

      rc = trexio_write_nucleus_charge(f, charge)
      CALL trexio_error(rc)
      DEALLOCATE (charge)

      rc = trexio_write_nucleus_label(f, label, 3)
      CALL trexio_error(rc)
      DEALLOCATE (label)

      ! nuclear repulsion energy well-defined for molecules only
      IF (SUM(cell%perd) == 0) THEN
         CALL nuclear_repulsion_energy(particle_set, e_nn)
         rc = trexio_write_nucleus_repulsion(f, e_nn)
         CALL trexio_error(rc)
      END IF

      !========================================================================================!
      ! Cell group
      !========================================================================================!

      rc = trexio_write_cell_a(f, cell%hmat(:,1))
      CALL trexio_error(rc)

      rc = trexio_write_cell_b(f, cell%hmat(:,2))
      CALL trexio_error(rc)

      rc = trexio_write_cell_c(f, cell%hmat(:,3))
      CALL trexio_error(rc)

      rc = trexio_write_cell_g_a(f, cell%h_inv(:,1))
      CALL trexio_error(rc)

      rc = trexio_write_cell_g_b(f, cell%h_inv(:,2))
      CALL trexio_error(rc)

      rc = trexio_write_cell_g_c(f, cell%h_inv(:,3))
      CALL trexio_error(rc)

      rc = trexio_write_cell_two_pi(f, 0)
      CALL trexio_error(rc)

      !========================================================================================!
      ! PBC group
      !========================================================================================!

      periodic = 0
      IF (SUM(cell%perd) /= 0) periodic = 1
      rc = trexio_write_pbc_periodic(f, periodic)
      CALL trexio_error(rc)

      IF (do_kpoints) THEN
         CALL get_kpoint_info(kpoints, nkp=nkp, nkp_grid=nkp_grid, wkp=wkp)

         rc = trexio_write_pbc_k_point_num(f, nkp)
         CALL trexio_error(rc)

         rc = trexio_write_pbc_k_point(f, REAL(nkp_grid, KIND=dp))
         CALL trexio_error(rc)

         rc = trexio_write_pbc_k_point_weight(f, wkp)
         CALL trexio_error(rc)
      END IF

      !========================================================================================!
      ! Electron group
      !========================================================================================!

      CALL get_qs_env(qs_env, nelectron_total=nel_tot)

      rc = trexio_write_electron_num(f, nel_tot)
      CALL trexio_error(rc)

      nspins = dft_control%nspins
      IF (nspins == 1) THEN
         ! it is a spin-restricted calculation and we need to split the electrons manually,
         ! because in CP2K they are all otherwise weirdly stored in nelectron_spin(1)
         nel_spin(1) = nel_tot/2
         nel_spin(2) = nel_tot/2
      ELSE
         ! for UKS/ROKS, the two spin channels are populated correctly and according to
         ! the multiplicity
         CALL get_qs_env(qs_env, nelectron_spin=nel_spin)
      END IF
      rc = trexio_write_electron_up_num(f, nel_spin(1))
      CALL trexio_error(rc)
      rc = trexio_write_electron_dn_num(f, nel_spin(2))
      CALL trexio_error(rc)


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


      !========================================================================================!
      ! *** QCSchema model ***
      !========================================================================================!

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


! **************************************************************************************************
!> \brief Handles TREXIO errors
!> \param rc the TREXIO return code
! **************************************************************************************************
   SUBROUTINE trexio_error(rc)
      INTEGER(trexio_exit_code), INTENT(IN) :: rc
      CHARACTER(LEN=128) :: err_msg

      IF (rc /= TREXIO_SUCCESS) THEN
         CALL trexio_string_of_error(rc, err_msg)
         CPABORT('TREXIO Error: '//TRIM(err_msg))
      END IF

   END SUBROUTINE trexio_error


! **************************************************************************************************
!> \brief Computes the nuclear repulsion energy of a molecular system
!> \param particle_set the set of particles in the system
!> \param e_nn the nuclear repulsion energy
! **************************************************************************************************
   SUBROUTINE nuclear_repulsion_energy(particle_set, e_nn)
      TYPE(particle_type), DIMENSION(:), POINTER, INTENT(IN) :: particle_set
      REAL(KIND=dp), INTENT(OUT) :: e_nn

      INTEGER :: i, j, z_i, z_j
      REAL(KIND=dp) :: r_ij

      e_nn = 0.0_dp
      DO i = 1, SIZE(particle_set)
         DO j = i+1, SIZE(particle_set)
            r_ij = NORM2(particle_set(i)%r - particle_set(j)%r)
            CALL get_atomic_kind(particle_set(i)%atomic_kind, z=z_i)
            CALL get_atomic_kind(particle_set(j)%atomic_kind, z=z_j)
            e_nn = e_nn + REAL(z_i, KIND=dp)*REAL(z_j, KIND=dp)/r_ij
         END DO
      END DO

   END SUBROUTINE nuclear_repulsion_energy

#endif

END MODULE trexio_utils
