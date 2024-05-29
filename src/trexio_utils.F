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
   USE basis_set_types, ONLY: gto_basis_set_type, get_gto_basis_set
   USE cell_types, ONLY: cell_type
   USE cp2k_info, ONLY: cp2k_version
   USE cp_control_types, ONLY: dft_control_type
   USE cp_files, ONLY: close_file, file_exists, open_file
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE external_potential_types, ONLY: sgp_potential_type, get_potential
   USE input_section_types, ONLY: section_vals_type, &
                                  section_vals_val_get
   !                                section_vals_get_subs_vals, &
   !                                section_vals_get
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp
   USE kpoint_types, ONLY: kpoint_type, get_kpoint_info
   USE particle_types, ONLY: particle_type
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   ! USE qs_force_types, ONLY: qs_force_type
   USE qs_kind_types, ONLY: get_qs_kind, get_qs_kind_set, &
                            qs_kind_type
   ! USE qs_ks_types, ONLY: qs_ks_env_type
   ! USE qs_scf_types, ONLY: qs_scf_env_type
   USE trexio, ONLY: trexio_open, trexio_close, &
                     TREXIO_HDF5, TREXIO_SUCCESS, &
                     trexio_string_of_error, trexio_t, trexio_exit_code, &
                     trexio_write_metadata_code, trexio_write_metadata_code_num, &
                     trexio_write_nucleus_coord, trexio_write_nucleus_num, &
                     trexio_write_nucleus_charge, trexio_write_nucleus_label, &
                     trexio_write_nucleus_repulsion, &
                     trexio_write_cell_a, trexio_write_cell_b, trexio_write_cell_c, &
                     trexio_write_cell_g_a, trexio_write_cell_g_b, &
                     trexio_write_cell_g_c, trexio_write_cell_two_pi, &
                     trexio_write_pbc_periodic, trexio_write_pbc_k_point_num, &
                     trexio_write_pbc_k_point, trexio_write_pbc_k_point_weight, &
                     trexio_write_electron_num, trexio_write_electron_up_num, &
                     trexio_write_electron_dn_num, &
                     trexio_write_state_num, trexio_write_state_id, &
                     trexio_write_state_energy, &
                     trexio_write_basis_type, trexio_write_basis_prim_num, &
                     trexio_write_basis_shell_num, trexio_write_basis_nucleus_index, &
                     trexio_write_basis_shell_ang_mom, trexio_write_basis_shell_factor, &
                     trexio_write_basis_r_power, trexio_write_basis_shell_index, &
                     trexio_write_basis_exponent, trexio_write_basis_coefficient, &
                     trexio_write_basis_prim_factor, &
                     trexio_write_ecp_z_core, trexio_write_ecp_max_ang_mom_plus_1, &
                     trexio_write_ecp_num, trexio_write_ecp_ang_mom, &
                     trexio_write_ecp_nucleus_index, trexio_write_ecp_exponent, &
                     trexio_write_ecp_coefficient, trexio_write_ecp_power
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
      CHARACTER(len=default_path_length)                 :: filename
      INTEGER(trexio_t)                                  :: f        ! The TREXIO file handle
      INTEGER(trexio_exit_code)                          :: rc       ! TREXIO return code
      LOGICAL                                            :: explicit, do_kpoints, ecp_semi_local, &
                                                            ecp_local, sgp_potential_present
      REAL(KIND=dp)                                      :: e_nn, zeff
      TYPE(cell_type), POINTER                           :: cell => Null()
      TYPE(cp_logger_type), POINTER                      :: logger => Null()
      TYPE(dft_control_type), POINTER                    :: dft_control => Null()
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(kpoint_type), POINTER                         :: kpoints => Null()
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set => Null()
      TYPE(qs_energy_type), POINTER                      :: energy => Null()
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set => Null()
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential => Null()
      ! TYPE(qs_ks_env_type), POINTER                      :: ks_env => Null()
      ! TYPE(qs_scf_env_type), POINTER                     :: scf_env => Null()

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE        :: label
      INTEGER                                            :: iatom, natoms, periodic, nkp, nel_tot, &
                                                            nspins, ikind, ishell_loc, ishell, &
                                                            shell_num, prim_num, nset, iset, ipgf, z, &
                                                            sl_lmax, ecp_num, nloc, nsemiloc, sl_l, iecp
      INTEGER, DIMENSION(2)                              :: nel_spin
      INTEGER, DIMENSION(3)                              :: nkp_grid
      INTEGER, DIMENSION(0:10)                           :: npot
      INTEGER, DIMENSION(:),    ALLOCATABLE              :: nucleus_index, shell_ang_mom, r_power, &
                                                            shell_index, z_core, max_ang_mom_plus_1, &
                                                            ang_mom, powers
      INTEGER, DIMENSION(:),    POINTER                  :: nshell => Null(), npgf => Null()
      INTEGER, DIMENSION(:, :), POINTER                  :: l => Null()
      REAL(KIND=dp), DIMENSION(:),   ALLOCATABLE         :: charge, shell_factor, exponents, coefficients, &
                                                            prim_factor
      REAL(KIND=dp), DIMENSION(:),   POINTER             :: wkp => Null()
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: coord
      REAL(KIND=dp), DIMENSION(:,:), POINTER             :: zetas => Null()
      REAL(KIND=dp), DIMENSION(:,:,:), POINTER           :: gcc => Null()

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(qs_env))

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

      !========================================================================================!
      ! Open the TREXIO file
      !========================================================================================!
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
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set, qs_kind_set=kind_set, natom=natoms)

      rc = trexio_write_nucleus_num(f, natoms)
      CALL trexio_error(rc)

      ALLOCATE (coord(3,natoms))
      ALLOCATE (label(natoms))
      ALLOCATE (charge(natoms))
      DO iatom = 1, natoms
         ! store the coordinates
         coord(iatom,:) = particle_set(iatom)%r(1:3)
         ! figure out the element symbol and to which kind_set entry this atomic_kind corresponds to
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=element_symbol, kind_number=ikind)
         ! store the element symbol
         label(iatom) = element_symbol
         ! get and store the effective nuclear charge of this kind_type (ikind)
         CALL get_qs_kind(kind_set(ikind), zeff=zeff)
         charge(iatom) = zeff
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
         CALL nuclear_repulsion_energy(particle_set, kind_set, e_nn)
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
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints, kpoints=kpoints)

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
      CALL get_qs_env(qs_env, dft_control=dft_control, nelectron_total=nel_tot)

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

      !========================================================================================!
      ! State group
      !========================================================================================!
      CALL get_qs_env(qs_env, energy=energy)

      rc = trexio_write_state_num(f, 1)
      CALL trexio_error(rc)

      rc = trexio_write_state_id(f, 1)
      CALL trexio_error(rc)

      rc = trexio_write_state_energy(f, energy%total)
      CALL trexio_error(rc)

      !========================================================================================!
      ! Basis group
      !========================================================================================!
      CALL get_qs_kind_set(kind_set, nshell=shell_num, npgf_seg=prim_num)

      rc = trexio_write_basis_type(f, 'Gaussian', LEN_TRIM('Gaussian')+1)
      CALL trexio_error(rc)

      rc = trexio_write_basis_shell_num(f, shell_num)
      CALL trexio_error(rc)

      rc = trexio_write_basis_prim_num(f, prim_num)
      CALL trexio_error(rc)

      ! one-to-one mapping between shells and ...
      ALLOCATE (nucleus_index(shell_num)) ! ...atomic indices
      ALLOCATE (shell_ang_mom(shell_num)) ! ...angular momenta
      ALLOCATE (shell_index(prim_num))    ! ...indices of primitive functions
      ALLOCATE (exponents(prim_num))      ! ...primitive exponents
      ALLOCATE (coefficients(prim_num))   ! ...contraction coefficients

      ishell = 0
      ipgf = 0
      DO iatom = 1, natoms
         ! get the qs_kind (index position in kind_set) for this atom (atomic_kind)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         ! get the primary (orbital) basis set associated to this qs_kind
         CALL get_qs_kind(kind_set(ikind), basis_set=basis_set, basis_type="ORB")
         ! get the info from the basis set
         CALL get_gto_basis_set(basis_set, &
                                nset=nset, &
                                nshell=nshell, &
                                npgf=npgf, &
                                zet=zetas, &
                                gcc=gcc, &
                                l=l)

         DO iset = 1, nset
            DO ishell_loc = 1, nshell(iset)
               ishell = ishell + 1

               ! nucleus_index array
               nucleus_index(ishell) = iatom

               ! shell_ang_mom array
               shell_ang_mom(ishell) = l(ishell_loc,iset)

               ! shell_index array
               shell_index(ipgf+1:ipgf+npgf(iset)) = ishell

               ! exponents array
               exponents(ipgf+1:ipgf+npgf(iset)) = zetas(1:npgf(iset), iset)

               ! coefficients array
               coefficients(ipgf+1:ipgf+npgf(iset)) = gcc(1:npgf(iset), ishell_loc, iset)

               ipgf = ipgf + npgf(iset)
            END DO
         END DO
      END DO
      ! just a failsafe check
      CPASSERT(ishell == shell_num)
      CPASSERT(ipgf == prim_num)

      rc = trexio_write_basis_nucleus_index(f, nucleus_index)
      CALL trexio_error(rc)
      DEALLOCATE (nucleus_index)

      rc = trexio_write_basis_shell_ang_mom(f, shell_ang_mom)
      CALL trexio_error(rc)
      DEALLOCATE (shell_ang_mom)

      ! Cartesian Gaussians are not normalized in CP2K
      ALLOCATE (shell_factor(shell_num))  ! 1-to-1 map bw shells and normalization factors
      shell_factor(:) = 1.0_dp
      rc = trexio_write_basis_shell_factor(f, shell_factor)
      CALL trexio_error(rc)
      DEALLOCATE (shell_factor)

      ! This is always 0 for Gaussian basis sets
      ALLOCATE (r_power(shell_num))       ! 1-to-1 map bw shells radial function powers
      r_power(:) = 0
      rc = trexio_write_basis_r_power(f, r_power)
      CALL trexio_error(rc)
      DEALLOCATE (r_power)

      rc = trexio_write_basis_shell_index(f, shell_index)
      CALL trexio_error(rc)
      DEALLOCATE (shell_index)

      rc = trexio_write_basis_exponent(f, exponents)
      CALL trexio_error(rc)
      DEALLOCATE (exponents)

      rc = trexio_write_basis_coefficient(f, coefficients)
      CALL trexio_error(rc)
      DEALLOCATE (coefficients)

      ! Primitive functions are not normalized in CP2K
      ALLOCATE (prim_factor(prim_num))
      prim_factor(:) = 1.0_dp
      rc = trexio_write_basis_prim_factor(f, prim_factor)
      CALL trexio_error(rc)
      DEALLOCATE (prim_factor)

      !========================================================================================!
      ! ECP group
      !========================================================================================!
      CALL get_qs_kind_set(kind_set, sgp_potential_present=sgp_potential_present)

      ! figure out whether we actually have ECP potentials
      IF (sgp_potential_present) THEN
         ecp_num = 0
         DO iatom = 1, natoms
            ! get the qs_kind (index position in kind_set) for this atom (atomic_kind)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
            ! get the the sgp_potential associated to this qs_kind
            CALL get_qs_kind(kind_set(ikind), sgp_potential=sgp_potential)

            ! get the info on the potential
            IF (ASSOCIATED(sgp_potential)) THEN
               CALL get_potential(potential=sgp_potential, ecp_local=ecp_local, ecp_semi_local=ecp_semi_local)
               IF (ecp_local) THEN
                  ! get number of local terms
                  CALL get_potential(potential=sgp_potential, nloc=nloc)
                  ecp_num = ecp_num + nloc
               END IF
               IF (ecp_semi_local) THEN
                  ! get number of semilocal terms
                  CALL get_potential(potential=sgp_potential, npot=npot)
                  ecp_num = ecp_num + SUM(npot)
               END IF
            END IF
         END DO
      END IF

      ! if we have ECP potentials, populate the ECP group
      IF (ecp_num > 0) THEN
         ALLOCATE (z_core(natoms))
         ALLOCATE (max_ang_mom_plus_1(natoms))
         max_ang_mom_plus_1(:) = 0

         ALLOCATE (ang_mom(ecp_num))
         ALLOCATE (nucleus_index(ecp_num))
         ALLOCATE (exponents(ecp_num))
         ALLOCATE (coefficients(ecp_num))
         ALLOCATE (powers(ecp_num))

         iecp = 0
         DO iatom = 1, natoms
            ! get the qs_kind (index position in kind_set) for this atom (atomic_kind)
            CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind, z=z)
            ! get the the sgp_potential associated to this qs_kind
            CALL get_qs_kind(kind_set(ikind), sgp_potential=sgp_potential, zeff=zeff)

            ! number of core electrons removed by the ECP
            z_core(iatom) = z - INT(zeff)

            ! get the info on the potential
            IF (ASSOCIATED(sgp_potential)) THEN
               CALL get_potential(potential=sgp_potential, ecp_local=ecp_local, ecp_semi_local=ecp_semi_local)

               ! deal with the local part
               IF (ecp_local) THEN
                  CALL get_potential(potential=sgp_potential, nloc=nloc, sl_lmax=sl_lmax)
                  ang_mom(iecp+1:iecp+nloc) = sl_lmax+1
                  nucleus_index(iecp+1:iecp+nloc) = iatom
                  exponents(iecp+1:iecp+nloc) = sgp_potential%bloc(1:nloc)
                  coefficients(iecp+1:iecp+nloc) = sgp_potential%aloc(1:nloc)
                  powers(iecp+1:iecp+nloc) = sgp_potential%nrloc(1:nloc) - 2
                  iecp = iecp + nloc
               END IF

               ! deal with the semilocal part
               IF (ecp_semi_local) THEN
                  CALL get_potential(potential=sgp_potential, npot=npot, sl_lmax=sl_lmax)
                  max_ang_mom_plus_1(iatom) = sl_lmax + 1

                  DO sl_l = 0, sl_lmax
                     nsemiloc = npot(sl_l)
                     ang_mom(iecp+1:iecp+nsemiloc) = sl_l
                     nucleus_index(iecp+1:iecp+nsemiloc) = iatom
                     exponents(iecp+1:iecp+nsemiloc) = sgp_potential%bpot(1:nsemiloc, sl_l)
                     coefficients(iecp+1:iecp+nsemiloc) = sgp_potential%apot(1:nsemiloc, sl_l)
                     powers(iecp+1:iecp+nsemiloc) = sgp_potential%nrpot(1:nsemiloc, sl_l) - 2
                     iecp = iecp + nsemiloc
                  END DO
               END IF
            END IF
         END DO

         ! fail-safe check
         CPASSERT(iecp == ecp_num)

         rc = trexio_write_ecp_num(f, ecp_num)
         CALL trexio_error(rc)

         rc = trexio_write_ecp_z_core(f, z_core)
         CALL trexio_error(rc)
         DEALLOCATE (z_core)

         rc = trexio_write_ecp_max_ang_mom_plus_1(f, max_ang_mom_plus_1)
         CALL trexio_error(rc)
         DEALLOCATE (max_ang_mom_plus_1)

         rc = trexio_write_ecp_ang_mom(f, ang_mom)
         CALL trexio_error(rc)
         DEALLOCATE (ang_mom)

         rc = trexio_write_ecp_nucleus_index(f, nucleus_index)
         CALL trexio_error(rc)
         DEALLOCATE (nucleus_index)

         rc = trexio_write_ecp_exponent(f, exponents)
         CALL trexio_error(rc)
         DEALLOCATE (exponents)

         rc = trexio_write_ecp_coefficient(f, coefficients)
         CALL trexio_error(rc)
         DEALLOCATE (coefficients)

         rc = trexio_write_ecp_power(f, powers)
         CALL trexio_error(rc)
         DEALLOCATE (powers)
      END IF

      !========================================================================================!
      ! Grid group
      !========================================================================================!
      ! TODO


      !========================================================================================!
      ! AO group
      !========================================================================================!


      !========================================================================================!
      ! MO group
      !========================================================================================!


      !========================================================================================!
      ! Close the TREXIO file
      !========================================================================================!
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
   SUBROUTINE nuclear_repulsion_energy(particle_set, kind_set, e_nn)
      TYPE(particle_type), DIMENSION(:), POINTER, INTENT(IN) :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER, INTENT(IN) :: kind_set
      REAL(KIND=dp), INTENT(OUT) :: e_nn

      INTEGER :: i, j, ikind, jkind, natoms
      REAL(KIND=dp) :: r_ij, zeff_i, zeff_j

      natoms = SIZE(particle_set)
      e_nn = 0.0_dp
      DO i = 1, natoms
         CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(kind_set(ikind), zeff=zeff_i)
         DO j = i+1, natoms
            r_ij = NORM2(particle_set(i)%r - particle_set(j)%r)

            CALL get_atomic_kind(particle_set(j)%atomic_kind, kind_number=jkind)
            CALL get_qs_kind(kind_set(jkind), zeff=zeff_j)

            e_nn = e_nn + zeff_i*zeff_j/r_ij
         END DO
      END DO

   END SUBROUTINE nuclear_repulsion_energy

#endif

END MODULE trexio_utils
