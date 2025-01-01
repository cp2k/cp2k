!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
   USE basis_set_types, ONLY: gto_basis_set_type, get_gto_basis_set
   USE cell_types, ONLY: cell_type
   USE cp2k_info, ONLY: cp2k_version
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_control_types, ONLY: dft_control_type
   USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm
   USE cp_files, ONLY: close_file, file_exists, open_file
   USE cp_fm_types, ONLY: cp_fm_get_info, cp_fm_type, cp_fm_create, cp_fm_set_all, &
                          cp_fm_get_submatrix, cp_fm_to_fm_submat_general, cp_fm_release
   USE cp_fm_struct, ONLY: cp_fm_struct_create, &
                           cp_fm_struct_release, &
                           cp_fm_struct_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE external_potential_types, ONLY: sgp_potential_type, get_potential
   USE input_section_types, ONLY: section_vals_type, section_vals_get, &
                                  section_vals_val_get
   USE kinds, ONLY: default_path_length, dp
   USE kpoint_types, ONLY: kpoint_type, kpoint_env_p_type, &
                           get_kpoint_info, get_kpoint_env
   USE mathconstants, ONLY: fourpi, pi
   USE message_passing, ONLY: mp_para_env_type
   USE orbital_pointers, ONLY: nco, nso
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types, ONLY: particle_type
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_kind_types, ONLY: get_qs_kind, get_qs_kind_set, &
                            qs_kind_type
   USE qs_mo_types, ONLY: mo_set_type, get_mo_set
#ifdef __TREXIO
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
                     trexio_write_ecp_coefficient, trexio_write_ecp_power, &
                     trexio_write_ao_cartesian, trexio_write_ao_num, &
                     trexio_write_ao_shell, trexio_write_ao_normalization, &
                     trexio_write_mo_num, trexio_write_mo_energy, &
                     trexio_write_mo_occupation, trexio_write_mo_spin, &
                     trexio_write_mo_class, trexio_write_mo_coefficient, &
                     trexio_write_mo_coefficient_im, trexio_write_mo_k_point, &
                     trexio_write_mo_type
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'trexio_utils'

   PUBLIC :: write_trexio

CONTAINS

! **************************************************************************************************
!> \brief Write a trexio file
!> \param qs_env the qs environment with all the info of the computation
!> \param trexio_section the section with the trexio info
! **************************************************************************************************
   SUBROUTINE write_trexio(qs_env, trexio_section)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: trexio_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_trexio'

      INTEGER                                            :: handle

#ifdef __TREXIO
      INTEGER                                            :: output_unit, unit_trexio
      CHARACTER(len=default_path_length)                 :: filename
      INTEGER(trexio_t)                                  :: f        ! The TREXIO file handle
      INTEGER(trexio_exit_code)                          :: rc       ! TREXIO return code
      LOGICAL                                            :: explicit, do_kpoints, ecp_semi_local, &
                                                            ecp_local, sgp_potential_present, ionode, &
                                                            use_real_wfn, save_cartesian
      REAL(KIND=dp)                                      :: e_nn, zeff, expzet, prefac, zeta, gcca
      TYPE(cell_type), POINTER                           :: cell => Null()
      TYPE(cp_logger_type), POINTER                      :: logger => Null()
      TYPE(dft_control_type), POINTER                    :: dft_control => Null()
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(kpoint_type), POINTER                         :: kpoints => Null()
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set => Null()
      TYPE(qs_energy_type), POINTER                      :: energy => Null()
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set => Null()
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential => Null()
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos => Null()
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos_kp => Null()
      TYPE(kpoint_env_p_type), DIMENSION(:), POINTER     :: kp_env => Null()
      TYPE(mp_para_env_type), POINTER                    :: para_env => Null(), para_env_inter_kp => Null()
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env => Null()
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct => Null()
      TYPE(cp_fm_type)                                   :: fm_mo_coeff, fm_dummy, fm_mo_coeff_im

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE        :: label
      INTEGER                                            :: iatom, natoms, periodic, nkp, nel_tot, &
                                                            nspins, ikind, ishell_loc, ishell, &
                                                            shell_num, prim_num, nset, iset, ipgf, z, &
                                                            sl_lmax, ecp_num, nloc, nsemiloc, sl_l, iecp, &
                                                            igf, icgf, ncgf, ngf_shell, lshell, ao_num, nmo, &
                                                            mo_num, ispin, ikp, imo, ikp_loc, nsgf, &
                                                            i, k, l, m
      INTEGER, DIMENSION(2)                              :: nel_spin, kp_range, nmo_spin
      INTEGER, DIMENSION(3)                              :: nkp_grid
      INTEGER, DIMENSION(0:10)                           :: npot
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nucleus_index, shell_ang_mom, r_power, &
                                                            shell_index, z_core, max_ang_mom_plus_1, &
                                                            ang_mom, powers, ao_shell, mo_spin, mo_kpoint
      INTEGER, DIMENSION(:), POINTER                     :: nshell => Null(), npgf => Null()
      INTEGER, DIMENSION(:, :), POINTER                  :: l_shell_set => Null()
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: charge, shell_factor, exponents, coefficients, &
                                                            prim_factor, ao_normalization, mo_energy, &
                                                            mo_occupation
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp => Null(), norm_cgf => Null()
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE        :: coord, mo_coefficient, mo_coefficient_im, &
                                                            mos_sgf, diag_nsgf, diag_ncgf, temp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: zetas => Null()
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc => Null()
#endif

      CALL timeset(routineN, handle)

#ifdef __TREXIO
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

      CALL get_qs_env(qs_env, para_env=para_env)
      ionode = para_env%is_source()

      ! inquire whether a file with the same name already exists, if yes, delete it
      IF (ionode) THEN
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

         rc = trexio_write_metadata_code(f, cp2k_version, LEN_TRIM(cp2k_version) + 1)
         CALL trexio_error(rc)

         !========================================================================================!
         ! Nucleus group
         !========================================================================================!
         CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set, qs_kind_set=kind_set, natom=natoms)

         rc = trexio_write_nucleus_num(f, natoms)
         CALL trexio_error(rc)

         ALLOCATE (coord(3, natoms))
         ALLOCATE (label(natoms))
         ALLOCATE (charge(natoms))
         DO iatom = 1, natoms
            ! store the coordinates
            coord(:, iatom) = particle_set(iatom)%r(1:3)
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
         rc = trexio_write_cell_a(f, cell%hmat(:, 1))
         CALL trexio_error(rc)

         rc = trexio_write_cell_b(f, cell%hmat(:, 2))
         CALL trexio_error(rc)

         rc = trexio_write_cell_c(f, cell%hmat(:, 3))
         CALL trexio_error(rc)

         rc = trexio_write_cell_g_a(f, cell%h_inv(:, 1))
         CALL trexio_error(rc)

         rc = trexio_write_cell_g_b(f, cell%h_inv(:, 2))
         CALL trexio_error(rc)

         rc = trexio_write_cell_g_c(f, cell%h_inv(:, 3))
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

      END IF ! ionode

      !========================================================================================!
      ! Basis group
      !========================================================================================!
      CALL get_qs_env(qs_env, qs_kind_set=kind_set, natom=natoms, particle_set=particle_set)
      CALL get_qs_kind_set(kind_set, nshell=shell_num, npgf_seg=prim_num)

      IF (ionode) THEN
         rc = trexio_write_basis_type(f, 'Gaussian', LEN_TRIM('Gaussian') + 1)
         CALL trexio_error(rc)

         rc = trexio_write_basis_shell_num(f, shell_num)
         CALL trexio_error(rc)

         rc = trexio_write_basis_prim_num(f, prim_num)
         CALL trexio_error(rc)
      END IF ! ionode

      ! one-to-one mapping between shells and ...
      ALLOCATE (nucleus_index(shell_num)) ! ...atomic indices
      ALLOCATE (shell_ang_mom(shell_num)) ! ...angular momenta
      ALLOCATE (shell_index(prim_num))    ! ...indices of primitive functions
      ALLOCATE (exponents(prim_num))      ! ...primitive exponents
      ALLOCATE (coefficients(prim_num))   ! ...contraction coefficients
      ALLOCATE (prim_factor(prim_num))    ! ...primitive normalization factors

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
                                l=l_shell_set)

         DO iset = 1, nset
            DO ishell_loc = 1, nshell(iset)
               ishell = ishell + 1

               ! nucleus_index array
               nucleus_index(ishell) = iatom

               ! shell_ang_mom array
               l = l_shell_set(ishell_loc, iset)
               shell_ang_mom(ishell) = l

               ! shell_index array
               shell_index(ipgf + 1:ipgf + npgf(iset)) = ishell

               ! exponents array
               exponents(ipgf + 1:ipgf + npgf(iset)) = zetas(1:npgf(iset), iset)

               ! compute on the fly the normalization factor as in normalise_gcc_orb
               ! and recover the original contraction coefficients to store them separately
               expzet = 0.25_dp*REAL(2*l + 3, dp)
               prefac = 2.0_dp**l*(2.0_dp/pi)**0.75_dp
               DO i = 1, npgf(iset)
                  gcca = gcc(i, ishell_loc, iset)
                  zeta = zetas(i, iset)

                  ! primitives normalization factors array
                  prim_factor(i + ipgf) = prefac*zeta**expzet

                  ! contractio coefficients array
                  coefficients(i + ipgf) = gcca/prim_factor(i + ipgf)
               END DO

               ipgf = ipgf + npgf(iset)
            END DO
         END DO
      END DO
      ! just a failsafe check
      CPASSERT(ishell == shell_num)
      CPASSERT(ipgf == prim_num)

      IF (ionode) THEN
         rc = trexio_write_basis_nucleus_index(f, nucleus_index)
         CALL trexio_error(rc)

         rc = trexio_write_basis_shell_ang_mom(f, shell_ang_mom)
         CALL trexio_error(rc)

         ! Normalization factors are shoved in the AO group
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

         rc = trexio_write_basis_exponent(f, exponents)
         CALL trexio_error(rc)

         rc = trexio_write_basis_coefficient(f, coefficients)
         CALL trexio_error(rc)

         ! Normalization factors are shoved in the AO group
         rc = trexio_write_basis_prim_factor(f, prim_factor)
         CALL trexio_error(rc)
      END IF

      DEALLOCATE (nucleus_index)
      DEALLOCATE (shell_index)
      DEALLOCATE (exponents)
      DEALLOCATE (coefficients)
      DEALLOCATE (prim_factor)
      ! shell_ang_mom is needed in the MO group, so will be deallocated there

      !========================================================================================!
      ! ECP group
      !========================================================================================!
      IF (ionode) THEN
         CALL get_qs_kind_set(kind_set, sgp_potential_present=sgp_potential_present)

         ! figure out whether we actually have ECP potentials
         ecp_num = 0
         IF (sgp_potential_present) THEN
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
                     ang_mom(iecp + 1:iecp + nloc) = sl_lmax + 1
                     nucleus_index(iecp + 1:iecp + nloc) = iatom
                     exponents(iecp + 1:iecp + nloc) = sgp_potential%bloc(1:nloc)
                     coefficients(iecp + 1:iecp + nloc) = sgp_potential%aloc(1:nloc)
                     powers(iecp + 1:iecp + nloc) = sgp_potential%nrloc(1:nloc) - 2
                     iecp = iecp + nloc
                  END IF

                  ! deal with the semilocal part
                  IF (ecp_semi_local) THEN
                     CALL get_potential(potential=sgp_potential, npot=npot, sl_lmax=sl_lmax)
                     max_ang_mom_plus_1(iatom) = sl_lmax + 1

                     DO sl_l = 0, sl_lmax
                        nsemiloc = npot(sl_l)
                        ang_mom(iecp + 1:iecp + nsemiloc) = sl_l
                        nucleus_index(iecp + 1:iecp + nsemiloc) = iatom
                        exponents(iecp + 1:iecp + nsemiloc) = sgp_potential%bpot(1:nsemiloc, sl_l)
                        coefficients(iecp + 1:iecp + nsemiloc) = sgp_potential%apot(1:nsemiloc, sl_l)
                        powers(iecp + 1:iecp + nsemiloc) = sgp_potential%nrpot(1:nsemiloc, sl_l) - 2
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

      END IF ! ionode

      !========================================================================================!
      ! Grid group
      !========================================================================================!
      ! TODO

      !========================================================================================!
      ! AO group
      !========================================================================================!
      CALL get_qs_env(qs_env, qs_kind_set=kind_set)
      CALL get_qs_kind_set(kind_set, ncgf=ncgf, nsgf=nsgf)

      CALL section_vals_val_get(trexio_section, "CARTESIAN", l_val=save_cartesian)
      IF (save_cartesian) THEN
         ao_num = ncgf
      ELSE
         ao_num = nsgf
      END IF

      IF (ionode) THEN
         IF (save_cartesian) THEN
            rc = trexio_write_ao_cartesian(f, 1)
         ELSE
            rc = trexio_write_ao_cartesian(f, 0)
         END IF
         CALL trexio_error(rc)

         rc = trexio_write_ao_num(f, ao_num)
         CALL trexio_error(rc)
      END IF

      ! one-to-one mapping between AOs and ...
      ALLOCATE (ao_shell(ao_num))         ! ..shells
      ALLOCATE (ao_normalization(ao_num)) ! ..normalization factors

      ! we need to be consistent with the basis group on the shell indices
      ishell = 0  ! global shell index
      igf = 0     ! global AO index
      DO iatom = 1, natoms
         ! get the qs_kind (index position in kind_set) for this atom (atomic_kind)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         ! get the primary (orbital) basis set associated to this qs_kind
         CALL get_qs_kind(kind_set(ikind), basis_set=basis_set, basis_type="ORB")
         ! get the info from the basis set
         CALL get_gto_basis_set(basis_set, &
                                nset=nset, &
                                nshell=nshell, &
                                norm_cgf=norm_cgf, &
                                ncgf=ncgf, &
                                nsgf=nsgf, &
                                l=l_shell_set)

         icgf = 0
         DO iset = 1, nset
            DO ishell_loc = 1, nshell(iset)
               ! global shell index
               ishell = ishell + 1
               ! angular momentum l of this shell
               lshell = l_shell_set(ishell_loc, iset)

               ! number of AOs in this shell
               IF (save_cartesian) THEN
                  ngf_shell = nco(lshell)
               ELSE
                  ngf_shell = nso(lshell)
               END IF

               ! one-to-one mapping between AOs and shells
               ao_shell(igf + 1:igf + ngf_shell) = ishell

               ! one-to-one mapping between AOs and normalization factors
               IF (save_cartesian) THEN
                  ao_normalization(igf + 1:igf + ngf_shell) = norm_cgf(icgf + 1:icgf + ngf_shell)
               ELSE
                  ! allocate some temporary arrays
                  ALLOCATE (diag_ncgf(nco(lshell), nco(lshell)))
                  ALLOCATE (diag_nsgf(nso(lshell), nso(lshell)))
                  ALLOCATE (temp(nso(lshell), nco(lshell)))
                  diag_ncgf = 0.0_dp
                  diag_nsgf = 0.0_dp
                  temp = 0.0_dp

                  DO i = 1, nco(lshell)
                     diag_ncgf(i, i) = norm_cgf(icgf + i)
                  END DO

                  ! transform the normalization factors from Cartesian to solid harmonics
                  temp(:, :) = MATMUL(orbtramat(lshell)%c2s, diag_ncgf)
                  diag_nsgf(:, :) = MATMUL(temp, TRANSPOSE(orbtramat(lshell)%s2c))
                  DO i = 1, nso(lshell)
                     ao_normalization(igf + i) = diag_nsgf(i, i)
                  END DO

                  DEALLOCATE (diag_ncgf)
                  DEALLOCATE (diag_nsgf)
                  DEALLOCATE (temp)
               END IF

               igf = igf + ngf_shell
               icgf = icgf + nco(lshell)
            END DO
         END DO
         ! just a failsafe check
         CPASSERT(icgf == ncgf)
      END DO

      IF (ionode) THEN
         rc = trexio_write_ao_shell(f, ao_shell)
         CALL trexio_error(rc)

         rc = trexio_write_ao_normalization(f, ao_normalization)
         CALL trexio_error(rc)
      END IF

      DEALLOCATE (ao_shell)
      DEALLOCATE (ao_normalization)

      !========================================================================================!
      ! MO group
      !========================================================================================!
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints, kpoints=kpoints, dft_control=dft_control, &
                      particle_set=particle_set, qs_kind_set=kind_set, blacs_env=blacs_env)
      nspins = dft_control%nspins
      CALL get_qs_kind_set(kind_set, nsgf=nsgf)
      nmo_spin = 0

      ! figure out that total number of MOs
      mo_num = 0
      IF (do_kpoints) THEN
         CALL get_kpoint_info(kpoints, kp_env=kp_env, nkp=nkp, use_real_wfn=use_real_wfn)
         CALL get_kpoint_env(kp_env(1)%kpoint_env, mos=mos_kp)
         DO ispin = 1, nspins
            CALL get_mo_set(mos_kp(1, ispin), nmo=nmo)
            nmo_spin(ispin) = nmo
         END DO
         mo_num = nkp*SUM(nmo_spin)

         ! we create a distributed fm matrix to gather the MOs from everywhere (in sph basis)
         CALL cp_fm_struct_create(fm_struct, para_env=para_env, context=blacs_env, &
                                  nrow_global=nsgf, ncol_global=mo_num)
         CALL cp_fm_create(fm_mo_coeff, fm_struct)
         CALL cp_fm_set_all(fm_mo_coeff, 0.0_dp)
         IF (.NOT. use_real_wfn) THEN
            CALL cp_fm_create(fm_mo_coeff_im, fm_struct)
            CALL cp_fm_set_all(fm_mo_coeff_im, 0.0_dp)
         END IF
         CALL cp_fm_struct_release(fm_struct)
      ELSE
         CALL get_qs_env(qs_env, mos=mos)
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), nmo=nmo)
            nmo_spin(ispin) = nmo
         END DO
         mo_num = SUM(nmo_spin)
      END IF

      ! allocate all the arrays
      ALLOCATE (mo_coefficient(ao_num, mo_num))
      mo_coefficient(:, :) = 0.0_dp
      ALLOCATE (mo_energy(mo_num))
      mo_energy(:) = 0.0_dp
      ALLOCATE (mo_occupation(mo_num))
      mo_occupation(:) = 0.0_dp
      ALLOCATE (mo_spin(mo_num))
      mo_spin(:) = 0
      ! extra arrays for kpoints
      IF (do_kpoints) THEN
         ALLOCATE (mo_coefficient_im(ao_num, mo_num))
         mo_coefficient_im(:, :) = 0.0_dp
         ALLOCATE (mo_kpoint(mo_num))
         mo_kpoint(:) = 0
      END IF

      ! in case of kpoints, we do this in 2 steps:
      ! 1. we gather the MOs of each kpt and pipe them into a single large distributed fm matrix;
      ! 2. we possibly transform the MOs of each kpt to Cartesian AOs and write them in the single large local array;
      IF (do_kpoints) THEN
         CALL get_kpoint_info(kpoints, kp_env=kp_env, nkp=nkp, kp_range=kp_range)

         DO ispin = 1, nspins
            DO ikp = 1, nkp
               nmo = nmo_spin(ispin)
               ! global index to store the MOs
               imo = (ikp - 1)*nmo + (ispin - 1)*nmo_spin(1)*nkp

               ! do I have this kpoint on this rank?
               IF (ikp >= kp_range(1) .AND. ikp <= kp_range(2)) THEN
                  ikp_loc = ikp - kp_range(1) + 1
                  ! get the mo set for this kpoint
                  CALL get_kpoint_env(kp_env(ikp_loc)%kpoint_env, mos=mos_kp)

                  ! if MOs are stored with dbcsr, copy them to fm
                  IF (mos_kp(1, ispin)%use_mo_coeff_b) THEN
                     CALL copy_dbcsr_to_fm(mos_kp(1, ispin)%mo_coeff_b, mos_kp(1, ispin)%mo_coeff)
                  END IF
                  ! copy real part of MO coefficients to large distributed fm matrix
                  CALL cp_fm_to_fm_submat_general(mos_kp(1, ispin)%mo_coeff, fm_mo_coeff, &
                                                  nsgf, nmo, 1, 1, 1, imo + 1, blacs_env)

                  ! copy MO energies to local arrays
                  mo_energy(imo + 1:imo + nmo) = mos_kp(1, ispin)%eigenvalues(1:nmo)

                  ! copy MO occupations to local arrays
                  mo_occupation(imo + 1:imo + nmo) = mos_kp(1, ispin)%occupation_numbers(1:nmo)

                  ! same for the imaginary part of MO coefficients
                  IF (.NOT. use_real_wfn) THEN
                     IF (mos_kp(2, ispin)%use_mo_coeff_b) THEN
                        CALL copy_dbcsr_to_fm(mos_kp(2, ispin)%mo_coeff_b, mos_kp(2, ispin)%mo_coeff)
                     END IF
                     CALL cp_fm_to_fm_submat_general(mos_kp(2, ispin)%mo_coeff, fm_mo_coeff_im, &
                                                     nsgf, nmo, 1, 1, 1, imo + 1, blacs_env)
                  END IF
               ELSE
                  ! call with a dummy fm for receiving the data
                  CALL cp_fm_to_fm_submat_general(fm_dummy, fm_mo_coeff, &
                                                  nsgf, nmo, 1, 1, 1, imo + 1, blacs_env)
                  IF (.NOT. use_real_wfn) THEN
                     CALL cp_fm_to_fm_submat_general(fm_dummy, fm_mo_coeff_im, &
                                                     nsgf, nmo, 1, 1, 1, imo + 1, blacs_env)
                  END IF
               END IF
            END DO
         END DO
      END IF

      ! reduce MO energies and occupations to the master node
      IF (do_kpoints) THEN
         CALL get_kpoint_info(kpoints, para_env_inter_kp=para_env_inter_kp)
         CALL para_env_inter_kp%sum(mo_energy)
         CALL para_env_inter_kp%sum(mo_occupation)
      END IF

      ! second step: here we actually put everything in the local arrays for writing to trexio
      DO ispin = 1, nspins
         ! get number of MOs for this spin
         nmo = nmo_spin(ispin)
         ! allocate local temp array to transform the MOs of each kpoint/spin
         ALLOCATE (mos_sgf(nsgf, nmo))

         IF (do_kpoints) THEN
            DO ikp = 1, nkp
               ! global index to store the MOs
               imo = (ikp - 1)*nmo + (ispin - 1)*nmo_spin(1)*nkp

               ! store kpoint index
               mo_kpoint(imo + 1:imo + nmo) = ikp
               ! store the MO spins
               mo_spin(imo + 1:imo + nmo) = ispin - 1

               ! transform and store the MO coefficients
               CALL cp_fm_get_submatrix(fm_mo_coeff, mos_sgf, 1, imo + 1, nsgf, nmo)
               IF (save_cartesian) THEN
                  CALL spherical_to_cartesian_mo(mos_sgf, particle_set, kind_set, mo_coefficient(:, imo + 1:imo + nmo))
               ELSE
                  ! we have to reorder the MOs since CP2K and TREXIO have different conventions
                  ! from m = -l, -l+1, ..., 0, ..., l-1, l   of CP2K
                  ! to   m =  0, +1, -1, +2, -2, ..., +l, -l of TREXIO
                  i = 0
                  DO ishell = 1, shell_num
                     l = shell_ang_mom(ishell)
                     DO k = 1, 2*l + 1
                        ! map from running index k: 1...2l+1 to magnetic quantum number m in TREXIO convention
                        m = (-1)**k*FLOOR(REAL(k, KIND=dp)/2.0_dp)
                        mo_coefficient(i + k, imo + 1:imo + nmo) = mos_sgf(i + l + m + 1, :)
                     END DO
                     i = i + 2*l + 1
                  END DO
                  CPASSERT(i == nsgf)
               END IF

               ! we have to do it for the imaginary part as well
               IF (.NOT. use_real_wfn) THEN
                  CALL cp_fm_get_submatrix(fm_mo_coeff_im, mos_sgf, 1, imo + 1, nsgf, nmo)
                  IF (save_cartesian) THEN
                     CALL spherical_to_cartesian_mo(mos_sgf, particle_set, kind_set, mo_coefficient_im(:, imo + 1:imo + nmo))
                  ELSE
                     ! we have to reorder the MOs since CP2K and TREXIO have different conventions
                     ! from m = -l, -l+1, ..., 0, ..., l-1, l   of CP2K
                     ! to   m =  0, +1, -1, +2, -2, ..., +l, -l of TREXIO
                     i = 0
                     DO ishell = 1, shell_num
                        l = shell_ang_mom(ishell)
                        DO k = 1, 2*l + 1
                           ! map from running index k: 1...2l+1 to magnetic quantum number m in TREXIO convention
                           m = (-1)**k*FLOOR(REAL(k, KIND=dp)/2.0_dp)
                           mo_coefficient_im(i + k, imo + 1:imo + nmo) = mos_sgf(i + l + m + 1, :)
                        END DO
                        i = i + 2*l + 1
                     END DO
                     CPASSERT(i == nsgf)
                  END IF
               END IF
            END DO
         ELSE ! no k-points
            ! global index to store the MOs
            imo = (ispin - 1)*nmo_spin(1)
            ! store the MO energies
            mo_energy(imo + 1:imo + nmo) = mos(ispin)%eigenvalues
            ! store the MO occupations
            mo_occupation(imo + 1:imo + nmo) = mos(ispin)%occupation_numbers
            ! store the MO spins
            mo_spin(imo + 1:imo + nmo) = ispin - 1

            ! check if we are using the dbcsr mo_coeff and copy them to fm if needed
            IF (mos(ispin)%use_mo_coeff_b) CALL copy_dbcsr_to_fm(mos(ispin)%mo_coeff_b, mos(ispin)%mo_coeff)

            ! allocate a normal fortran array to store the spherical MO coefficients
            CALL cp_fm_get_submatrix(mos(ispin)%mo_coeff, mos_sgf)

            IF (save_cartesian) THEN
               CALL spherical_to_cartesian_mo(mos_sgf, particle_set, kind_set, mo_coefficient(:, imo + 1:imo + nmo))
            ELSE
               ! we have to reorder the MOs since CP2K and TREXIO have different conventions
               ! from m = -l, -l+1, ..., 0, ..., l-1, l   of CP2K
               ! to   m =  0, +1, -1, +2, -2, ..., +l, -l of TREXIO
               i = 0
               DO ishell = 1, shell_num
                  l = shell_ang_mom(ishell)
                  DO k = 1, 2*l + 1
                     ! map from running index k: 1...2l+1 to magnetic quantum number m in TREXIO convention
                     m = (-1)**k*FLOOR(REAL(k, KIND=dp)/2.0_dp)
                     mo_coefficient(i + k, imo + 1:imo + nmo) = mos_sgf(i + l + m + 1, :)
                  END DO
                  i = i + 2*l + 1
               END DO
               CPASSERT(i == nsgf)
            END IF
         END IF

         DEALLOCATE (mos_sgf)
      END DO

      IF (ionode) THEN
         rc = trexio_write_mo_type(f, 'Canonical', LEN_TRIM('Canonical') + 1)
         CALL trexio_error(rc)

         rc = trexio_write_mo_num(f, mo_num)
         CALL trexio_error(rc)

         rc = trexio_write_mo_coefficient(f, mo_coefficient)
         CALL trexio_error(rc)

         rc = trexio_write_mo_energy(f, mo_energy)
         CALL trexio_error(rc)

         rc = trexio_write_mo_occupation(f, mo_occupation)
         CALL trexio_error(rc)

         rc = trexio_write_mo_spin(f, mo_spin)
         CALL trexio_error(rc)

         IF (do_kpoints) THEN
            rc = trexio_write_mo_coefficient_im(f, mo_coefficient_im)
            CALL trexio_error(rc)

            rc = trexio_write_mo_k_point(f, mo_kpoint)
            CALL trexio_error(rc)
         END IF
      END IF

      DEALLOCATE (shell_ang_mom)
      DEALLOCATE (mo_coefficient)
      DEALLOCATE (mo_energy)
      DEALLOCATE (mo_occupation)
      DEALLOCATE (mo_spin)
      IF (do_kpoints) THEN
         DEALLOCATE (mo_coefficient_im)
         DEALLOCATE (mo_kpoint)
         CALL cp_fm_release(fm_mo_coeff)
         CALL cp_fm_release(fm_mo_coeff_im)
      END IF

      !========================================================================================!
      ! RDM group
      !========================================================================================!
      !TODO

      !========================================================================================!
      ! Close the TREXIO file
      !========================================================================================!
      IF (ionode) THEN
         rc = trexio_close(f)
         CALL trexio_error(rc)
      END IF
#else
      MARK_USED(qs_env)
      MARK_USED(trexio_section)
      CPWARN('TREXIO support has not been enabled in this build.')
#endif
      CALL timestop(handle)

   END SUBROUTINE write_trexio

#ifdef __TREXIO
! **************************************************************************************************
!> \brief Handles TREXIO errors
!> \param rc the TREXIO return code
! **************************************************************************************************
   SUBROUTINE trexio_error(rc)
      INTEGER(trexio_exit_code), INTENT(IN)              :: rc

      CHARACTER(LEN=128)                                 :: err_msg

      IF (rc /= TREXIO_SUCCESS) THEN
         CALL trexio_string_of_error(rc, err_msg)
         CPABORT('TREXIO Error: '//TRIM(err_msg))
      END IF

   END SUBROUTINE trexio_error

! **************************************************************************************************
!> \brief Computes the nuclear repulsion energy of a molecular system
!> \param particle_set the set of particles in the system
!> \param kind_set the set of qs_kinds in the system
!> \param e_nn the nuclear repulsion energy
! **************************************************************************************************
   SUBROUTINE nuclear_repulsion_energy(particle_set, kind_set, e_nn)
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: kind_set
      REAL(KIND=dp), INTENT(OUT)                         :: e_nn

      INTEGER                                            :: i, ikind, j, jkind, natoms
      REAL(KIND=dp)                                      :: r_ij, zeff_i, zeff_j

      natoms = SIZE(particle_set)
      e_nn = 0.0_dp
      DO i = 1, natoms
         CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(kind_set(ikind), zeff=zeff_i)
         DO j = i + 1, natoms
            r_ij = NORM2(particle_set(i)%r - particle_set(j)%r)

            CALL get_atomic_kind(particle_set(j)%atomic_kind, kind_number=jkind)
            CALL get_qs_kind(kind_set(jkind), zeff=zeff_j)

            e_nn = e_nn + zeff_i*zeff_j/r_ij
         END DO
      END DO

   END SUBROUTINE nuclear_repulsion_energy

! **************************************************************************************************
!> \brief Computes a spherical to cartesian MO transformation (solid harmonics in reality)
!> \param mos_sgf the MO coefficients in spherical AO basis
!> \param particle_set the set of particles in the system
!> \param qs_kind_set the set of qs_kinds in the system
!> \param mos_cgf the transformed MO coefficients in Cartesian AO basis
! **************************************************************************************************
   SUBROUTINE spherical_to_cartesian_mo(mos_sgf, particle_set, qs_kind_set, mos_cgf)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: mos_sgf
      TYPE(particle_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: mos_cgf

      INTEGER                                            :: iatom, icgf, ikind, iset, isgf, ishell, &
                                                            lshell, ncgf, nmo, nset, nsgf
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: l
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      CALL get_qs_kind_set(qs_kind_set, ncgf=ncgf, nsgf=nsgf)

      mos_cgf = 0.0_dp
      nmo = SIZE(mos_sgf, 2)

      ! Transform spherical MOs to Cartesian MOs
      icgf = 1
      isgf = 1
      DO iatom = 1, SIZE(particle_set)
         NULLIFY (orb_basis_set)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)

         IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                   nset=nset, &
                                   nshell=nshell, &
                                   l=l)
            DO iset = 1, nset
               DO ishell = 1, nshell(iset)
                  lshell = l(ishell, iset)
                  CALL dgemm("T", "N", nco(lshell), nmo, nso(lshell), 1.0_dp, &
                             orbtramat(lshell)%c2s, nso(lshell), &
                             mos_sgf(isgf, 1), nsgf, 0.0_dp, &
                             mos_cgf(icgf, 1), ncgf)
                  icgf = icgf + nco(lshell)
                  isgf = isgf + nso(lshell)
               END DO
            END DO
         ELSE
            ! assume atom without basis set
            CPABORT("Unknown basis set type")
         END IF
      END DO ! iatom

   END SUBROUTINE spherical_to_cartesian_mo
#endif

END MODULE trexio_utils
