!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Environment for NEGF based quantum transport calculations
! **************************************************************************************************
MODULE negf_env_types
   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_p_type,&
                                              force_env_type,&
                                              use_qs_force
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_env,&
                                              get_kpoint_info,&
                                              kpoint_env_p_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_sum
   USE negf_atom_map,                   ONLY: negf_atom_map_type,&
                                              negf_map_atomic_indices
   USE negf_control_types,              ONLY: negf_control_contact_type,&
                                              negf_control_type
   USE negf_matrix_utils,               ONLY: invert_cell_to_index,&
                                              negf_copy_contact_matrix,&
                                              negf_copy_sym_dbcsr_to_fm_submat,&
                                              number_of_atomic_orbitals
   USE negf_qs_env_types,               ONLY: negf_qs_env_type
   USE negf_subgroup_types,             ONLY: negf_subgroup_env_type
   USE negf_vectors,                    ONLY: contact_direction_vector,&
                                              projection_on_direction_vector
   USE particle_types,                  ONLY: particle_type
   USE qs_density_mixing_types,         ONLY: mixing_storage_create,&
                                              mixing_storage_release,&
                                              mixing_storage_type
   USE qs_energy,                       ONLY: qs_energies
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_env_types'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: negf_env_type, negf_env_contact_type
   PUBLIC :: negf_env_create, negf_env_release

! **************************************************************************************************
!> \brief  Contact-specific NEGF environment.
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_env_contact_type
      REAL(kind=dp), DIMENSION(3)                        :: direction_vector, origin
      !> index if the axis in a contact unit cell which coinsides with the transport (contact) direction
      !> 0 (undefined), 1 (x), 2 (y), 3 (z)
      INTEGER                                            :: direction_axis
      !> atoms belonging to a primary contact unit cell
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist_cell0
      !> atoms belonging to a secondary contact unit cell (will be removed one day ...)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist_cell1
      !> list of equivalent atoms in an appropriate contact force environment
      TYPE(negf_atom_map_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: atom_map_cell0
      !> energy of the HOMO
      REAL(kind=dp)                                      :: homo_energy
      !> diagonal (h_00) and off-diagonal (h_01) blocks of the contact Kohn-Sham matrix ([number_of_spins]).
      !> The matrix h_01 is of the shape [nao_cell0 x nao_cell1]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: h_00, h_01
      !> diagonal and off-diagonal blocks of the density matrix
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: rho_00, rho_01
      !> diagonal and off-diagonal blocks of the overlap matrix
      TYPE(cp_fm_type), POINTER                          :: s_00 => null(), s_01 => null()
   END TYPE negf_env_contact_type

! **************************************************************************************************
!> \brief  NEGF environment.
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_env_type
      !> contact-specific NEGF environments
      TYPE(negf_env_contact_type), ALLOCATABLE, &
         DIMENSION(:)                                     :: contacts
      !> Kohn-Sham matrix of the scattering region
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)       :: h_s
      !> Kohn-Sham matrix of the scattering region -- contact interface ([nspins, ncontacts])
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)    :: h_sc
      !> overlap matrix of the scattering region
      TYPE(cp_fm_type), POINTER                           :: s_s => null()
      !> overlap matrix of the scattering region -- contact interface for every contact ([ncontacts])
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)       :: s_sc
      !> structure needed for density mixing
      TYPE(mixing_storage_type), POINTER                  :: mixing_storage
      !> density mixing method
      INTEGER                                             :: mixing_method
   END TYPE negf_env_type

! **************************************************************************************************
!> \brief  Allocatable list of the type 'negf_atom_map_type'.
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE negf_atom_map_contact_type
      TYPE(negf_atom_map_type), ALLOCATABLE, DIMENSION(:) :: atom_map
   END TYPE negf_atom_map_contact_type

CONTAINS

! **************************************************************************************************
!> \brief Create a new NEGF environment and compute the relevant Kohn-Sham matrices.
!> \param negf_env            NEGF environment to create
!> \param sub_env             NEGF parallel (sub)group environment
!> \param negf_control        NEGF control
!> \param force_env           the primary force environment
!> \param negf_mixing_section pointer to a mixing section within the NEGF input section
!> \param log_unit            output unit number
!> \par History
!>   * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_env_create(negf_env, sub_env, negf_control, force_env, negf_mixing_section, log_unit)
      TYPE(negf_env_type), INTENT(inout)                 :: negf_env
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(section_vals_type), POINTER                   :: negf_mixing_section
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_create', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: contact_str, force_env_str, &
                                                            n_force_env_str
      INTEGER                                            :: handle, icontact, in_use, ispin, &
                                                            n_force_env, ncontacts, nspins
      LOGICAL                                            :: do_kpoints
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_p_type), DIMENSION(:), POINTER      :: sub_force_env
      TYPE(negf_atom_map_contact_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: map_contact
      TYPE(negf_qs_env_type), POINTER                    :: negf_qs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env, qs_env_contact
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(qs_subsys_type), POINTER                      :: subsys, subsys_contact

      CALL timeset(routineN, handle)

      ! ensure we have Quickstep enabled for all force_env
      NULLIFY (sub_force_env)
      CALL force_env_get(force_env, in_use=in_use, qs_env=qs_env, sub_force_env=sub_force_env)

      IF (ASSOCIATED(sub_force_env)) THEN
         n_force_env = SIZE(sub_force_env)
      ELSE
         n_force_env = 0
      END IF

      IF (in_use == use_qs_force) THEN
         DO icontact = 1, n_force_env
            CALL force_env_get(sub_force_env(icontact)%force_env, in_use=in_use)
            IF (in_use /= use_qs_force) EXIT
         END DO
      END IF

      IF (in_use /= use_qs_force) THEN
         CPABORT("Quickstep is required for NEGF run.")
      END IF

      ! check that all mentioned FORCE_EVAL sections are actually present
      ncontacts = SIZE(negf_control%contacts)

      DO icontact = 1, ncontacts
         IF (negf_control%contacts(icontact)%force_env_index > n_force_env .OR. &
             negf_control%contacts(icontact)%force_env_index <= 0) THEN
            WRITE (contact_str, '(I11)') icontact
            WRITE (force_env_str, '(I11)') negf_control%contacts(icontact)%force_env_index
            WRITE (n_force_env_str, '(I11)') n_force_env

            CALL cp_abort(__LOCATION__, &
                          "Contact number "//TRIM(ADJUSTL(contact_str))//" is linked with the FORCE_EVAL section number "// &
                          TRIM(ADJUSTL(force_env_str))//", however only "//TRIM(ADJUSTL(n_force_env_str))// &
                          " FORCE_EVAL sections have been found. Note that FORCE_EVAL sections are enumerated from 0"// &
                          " and that the primary (0-th) section must contain all the atoms.")
         END IF
      END DO

      ! create basic matrices and neighbour lists for the primary force_env,
      ! so we know how matrix elements are actually distributed across CPUs.
      CALL qs_energies_init(qs_env, calc_forces=.FALSE.)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, do_kpoints=do_kpoints, &
                      matrix_ks_kp=matrix_ks_kp, negf_env=negf_qs_env, subsys=subsys)

      IF (do_kpoints) THEN
         CPABORT("k-points are currently not supported for device FORCE_EVAL")
      END IF

      ! stage 1: map the atoms between the device force_env and all contact force_env-s
      ALLOCATE (negf_env%contacts(ncontacts))
      ALLOCATE (map_contact(ncontacts))

      DO icontact = 1, ncontacts
         CALL force_env_get(sub_force_env(negf_control%contacts(icontact)%force_env_index)%force_env, qs_env=qs_env_contact)
         CALL get_qs_env(qs_env_contact, subsys=subsys_contact)

         CALL negf_env_contact_init_maps(contact_env=negf_env%contacts(icontact), &
                                         contact_control=negf_control%contacts(icontact), &
                                         atom_map=map_contact(icontact)%atom_map, &
                                         eps_geometry=negf_control%eps_geometry, &
                                         subsys_device=subsys, &
                                         subsys_contact=subsys_contact)

         IF (negf_env%contacts(icontact)%direction_axis == 0) THEN
            WRITE (contact_str, '(I11)') icontact
            WRITE (force_env_str, '(I11)') negf_control%contacts(icontact)%force_env_index
            CALL cp_abort(__LOCATION__, &
                          "One lattice vector of the contact unit cell (FORCE_EVAL section "// &
                          TRIM(ADJUSTL(force_env_str))//") must be parallel to the direction of the contact "// &
                          TRIM(ADJUSTL(contact_str))//".")
         END IF
      END DO

      ! obtain a Konh-Sham and overlap matrices for each contact
      ALLOCATE (negf_qs_env%contacts(ncontacts))

      ! The following code is kept commented until we find a way to obtain matrix elements in a consistent way
      ! ***********************************************************************************************************************
      !
      ! stage 2: obtain relevant Kohn-Sham matrix blocks for each contact
      DO icontact = 1, ncontacts
         IF (log_unit > 0) &
            WRITE (log_unit, '(/,T2,A,T70,I11)') "NEGF| Construct the Kohn-Sham matrix for the bulk contact ", icontact

         CALL force_env_get(sub_force_env(negf_control%contacts(icontact)%force_env_index)%force_env, qs_env=qs_env_contact)

         CALL qs_energies(qs_env_contact, consistent_energies=.FALSE., calc_forces=.FALSE.)

         CALL negf_env_contact_init_matrices(contact_env=negf_env%contacts(icontact), &
                                             sub_env=sub_env, qs_env_contact=qs_env_contact)

         !
         !  CALL negf_get_bulk_contact_matrix(negf_control_aux%contacts(icontact)%ks_block_map, &
         !                                    negf_control_aux%contacts(icontact)%ks_block_storage,&
         !                                    matrix_ks_kp, negf_control%contacts(icontact)%atomlist_bulk, &
         !                                    map_contact(icontact)%atom_map, &
         !                                    matrix_ks_kp_contact, cell_to_index, blacs_env%para_env)

         negf_qs_env%contacts(icontact)%direction_vector = negf_env%contacts(icontact)%direction_vector
         negf_qs_env%contacts(icontact)%origin_scattering = negf_env%contacts(icontact)%origin
         negf_qs_env%contacts(icontact)%v_external = negf_control%contacts(icontact)%v_external
      END DO
      !negf_control_aux%update_ks_matrix = .TRUE.
      ! ***********************************************************************************************************************

      ! obtain an initial KS-matrix for the scattering region

      CALL qs_energies(qs_env, consistent_energies=.FALSE., calc_forces=.FALSE.)
      CALL negf_env_device_init_matrices(negf_env, negf_control, sub_env, qs_env)

      ! stage2 (alternative): extract the appropriate parts of the Kohn-Sham matrix which
      !                       related to the corresponding contact
      ! Will be removed later
      DO icontact = 1, ncontacts
         CALL get_qs_env(qs_env, blacs_env=blacs_env, matrix_s_kp=matrix_s_kp, matrix_ks_kp=matrix_ks_kp)
         nspins = SIZE(negf_env%contacts(1)%h_00)

         CALL get_qs_env(qs_env, rho=rho_struct)
         CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)

         DO ispin = 1, nspins
            CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_ks_kp(ispin, 1)%matrix, &
                                                  fm=negf_env%contacts(icontact)%h_00(ispin)%matrix, &
                                                  atomlist_row=negf_env%contacts(icontact)%atomlist_cell0, &
                                                  atomlist_col=negf_env%contacts(icontact)%atomlist_cell0, &
                                                  subsys=subsys, mpi_comm_global=blacs_env%para_env%group, &
                                                  do_upper_diag=.TRUE., do_lower=.TRUE.)

            CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_ks_kp(ispin, 1)%matrix, &
                                                  fm=negf_env%contacts(icontact)%h_01(ispin)%matrix, &
                                                  atomlist_row=negf_env%contacts(icontact)%atomlist_cell0, &
                                                  atomlist_col=negf_env%contacts(icontact)%atomlist_cell1, &
                                                  subsys=subsys, mpi_comm_global=blacs_env%para_env%group, &
                                                  do_upper_diag=.TRUE., do_lower=.TRUE.)
         END DO

         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_s_kp(1, 1)%matrix, &
                                               fm=negf_env%contacts(icontact)%s_00, &
                                               atomlist_row=negf_env%contacts(icontact)%atomlist_cell0, &
                                               atomlist_col=negf_env%contacts(icontact)%atomlist_cell0, &
                                               subsys=subsys, mpi_comm_global=blacs_env%para_env%group, &
                                               do_upper_diag=.TRUE., do_lower=.TRUE.)

         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_s_kp(1, 1)%matrix, &
                                               fm=negf_env%contacts(icontact)%s_01, &
                                               atomlist_row=negf_env%contacts(icontact)%atomlist_cell0, &
                                               atomlist_col=negf_env%contacts(icontact)%atomlist_cell1, &
                                               subsys=subsys, mpi_comm_global=blacs_env%para_env%group, &
                                               do_upper_diag=.TRUE., do_lower=.TRUE.)
      END DO

      ! electron density mixing;
      ! the input section below should be consistent with the subroutine create_negf_section()
      NULLIFY (negf_env%mixing_storage)
      CALL section_vals_val_get(negf_mixing_section, "METHOD", i_val=negf_env%mixing_method)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      CALL mixing_storage_create(negf_env%mixing_storage, negf_mixing_section, &
                                 negf_env%mixing_method, dft_control%qs_control%cutoff)

      CALL timestop(handle)
   END SUBROUTINE negf_env_create

! **************************************************************************************************
!> \brief Establish mapping between the primary and the contact force environments
!> \param contact_env         NEGF environment for the given contact (modified on exit)
!> \param contact_control     NEGF control
!> \param atom_map            atomic map
!> \param eps_geometry        accuracy in mapping atoms between different force environments
!> \param subsys_device       QuickStep subsystem of the device force environment
!> \param subsys_contact      QuickStep subsystem of the contact force environment
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_env_contact_init_maps(contact_env, contact_control, atom_map, &
                                         eps_geometry, subsys_device, subsys_contact)
      TYPE(negf_env_contact_type), INTENT(inout)         :: contact_env
      TYPE(negf_control_contact_type), INTENT(in)        :: contact_control
      TYPE(negf_atom_map_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(inout)                     :: atom_map
      REAL(kind=dp), INTENT(in)                          :: eps_geometry
      TYPE(qs_subsys_type), POINTER                      :: subsys_device, subsys_contact

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_env_contact_init_maps', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, natoms

      CALL timeset(routineN, handle)

      CALL contact_direction_vector(contact_env%origin, &
                                    contact_env%direction_vector, &
                                    contact_control%atomlist_screening, &
                                    contact_control%atomlist_bulk, &
                                    subsys_device)

      contact_env%direction_axis = contact_direction_axis(contact_env%direction_vector, subsys_contact, eps_geometry)

      IF (contact_env%direction_axis > 0) THEN
         natoms = SIZE(contact_control%atomlist_bulk)
         ALLOCATE (atom_map(natoms))

         ! map atom listed in 'contact_control%atomlist_bulk' to the corresponding atom/cell replica from the contact force_env
         CALL negf_map_atomic_indices(atom_map=atom_map, &
                                      atom_list=contact_control%atomlist_bulk, &
                                      subsys_device=subsys_device, &
                                      subsys_contact=subsys_contact, &
                                      eps_geometry=eps_geometry)

         ! list atoms from 'contact_control%atomlist_bulk' which belong to
         ! the primary unit cell of the bulk region for the given contact
         CALL list_atoms_in_bulk_primary_unit_cell(atomlist_cell0=contact_env%atomlist_cell0, &
                                                   atom_map_cell0=contact_env%atom_map_cell0, &
                                                   atomlist_bulk=contact_control%atomlist_bulk, &
                                                   atom_map=atom_map, &
                                                   origin=contact_env%origin, &
                                                   direction_vector=contact_env%direction_vector, &
                                                   direction_axis=contact_env%direction_axis, &
                                                   subsys_device=subsys_device)

         ! currently we extract contacts' off-diagonal matrix blocks using matrices from the device region,
         ! so we need to know which atoms belong to the "second replica" of the contact unit cell.
         ! Currently we can extract the corresponding blocks from the contact force_env as well,
         ! however each matrix elements obtained this way have a random phase with respect to the device's
         ! force_env. This is not a problem for multiple force environment DFT calculations but it is for NEGF.
         CALL list_atoms_in_bulk_secondary_unit_cell(atomlist_cell1=contact_env%atomlist_cell1, &
                                                     atomlist_bulk=contact_control%atomlist_bulk, &
                                                     atom_map=atom_map, &
                                                     origin=contact_env%origin, &
                                                     direction_vector=contact_env%direction_vector, &
                                                     direction_axis=contact_env%direction_axis, &
                                                     subsys_device=subsys_device)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_env_contact_init_maps

! **************************************************************************************************
!> \brief Extract relevant matrix blocks for the given contact.
!> \param contact_env         NEGF environment for the contact (modified on exit)
!> \param sub_env             NEGF parallel (sub)group environment
!> \param qs_env_contact      QuickStep environment for the contact force environment
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_env_contact_init_matrices(contact_env, sub_env, qs_env_contact)
      TYPE(negf_env_contact_type), INTENT(inout)         :: contact_env
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      TYPE(qs_environment_type), POINTER                 :: qs_env_contact

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_env_contact_init_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ispin, nao, natoms, &
                                                            nimages, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env_contact, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      matrix_ks_kp=matrix_ks_kp, &
                      matrix_s_kp=matrix_s_kp, &
                      para_env=para_env, &
                      subsys=subsys)
      CALL get_qs_env(qs_env_contact, rho=rho_struct)
      CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)

      CALL negf_homo_energy_estimate(contact_env%homo_energy, qs_env_contact)

      natoms = SIZE(contact_env%atomlist_cell0)
      ALLOCATE (atom_list(natoms))
      DO iatom = 1, natoms
         atom_list(iatom) = contact_env%atom_map_cell0(iatom)%iatom

         ! with no k-points there is one-to-one correspondence between the primary unit cell
         ! of the contact force_env and the first contact unit cell of the device force_env
         IF (SUM(ABS(contact_env%atom_map_cell0(iatom)%cell(:))) > 0) THEN
            CPABORT("NEGF K-points are not currently supported")
         END IF
      END DO

      nspins = dft_control%nspins
      nimages = dft_control%nimages

      IF (do_kpoints) THEN
         CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
      ELSE
         ALLOCATE (cell_to_index(0:0, 0:0, 0:0))
         cell_to_index(0, 0, 0) = 1
      END IF

      ALLOCATE (index_to_cell(3, nimages))
      CALL invert_cell_to_index(cell_to_index, nimages, index_to_cell)
      IF (.NOT. do_kpoints) DEALLOCATE (cell_to_index)

      NULLIFY (fm_struct)
      nao = number_of_atomic_orbitals(subsys, atom_list)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=sub_env%blacs_env)

      ! ++ create matrices: s_00, s_01
      NULLIFY (contact_env%s_00, contact_env%s_01)
      CALL cp_fm_create(contact_env%s_00, fm_struct)
      CALL cp_fm_create(contact_env%s_01, fm_struct)

      ! ++ create matrices: h_00, h_01, rho_00, rho_01
      ALLOCATE (contact_env%h_00(nspins), contact_env%h_01(nspins))
      ALLOCATE (contact_env%rho_00(nspins), contact_env%rho_01(nspins))
      DO ispin = 1, nspins
         NULLIFY (contact_env%h_00(ispin)%matrix)
         NULLIFY (contact_env%h_01(ispin)%matrix)
         NULLIFY (contact_env%rho_00(ispin)%matrix)
         NULLIFY (contact_env%rho_01(ispin)%matrix)

         CALL cp_fm_create(contact_env%h_00(ispin)%matrix, fm_struct)
         CALL cp_fm_create(contact_env%h_01(ispin)%matrix, fm_struct)
         CALL cp_fm_create(contact_env%rho_00(ispin)%matrix, fm_struct)
         CALL cp_fm_create(contact_env%rho_01(ispin)%matrix, fm_struct)
      END DO

      CALL cp_fm_struct_release(fm_struct)

      ! extract matrices: s_00, s_01
      CALL negf_copy_contact_matrix(fm_cell0=contact_env%s_00, fm_cell1=contact_env%s_01, &
                                    direction_axis=contact_env%direction_axis, &
                                    matrix_kp=matrix_s_kp(1, :), &
                                    index_to_cell=index_to_cell, atom_list=atom_list, &
                                    subsys=subsys, mpi_comm_global=para_env%group)

      ! extract matrices: h_00, h_01, rho_00, rho_01
      DO ispin = 1, nspins
         CALL negf_copy_contact_matrix(fm_cell0=contact_env%h_00(ispin)%matrix, &
                                       fm_cell1=contact_env%h_01(ispin)%matrix, &
                                       direction_axis=contact_env%direction_axis, &
                                       matrix_kp=matrix_ks_kp(ispin, :), &
                                       index_to_cell=index_to_cell, atom_list=atom_list, &
                                       subsys=subsys, mpi_comm_global=para_env%group)

         CALL negf_copy_contact_matrix(fm_cell0=contact_env%rho_00(ispin)%matrix, &
                                       fm_cell1=contact_env%rho_01(ispin)%matrix, &
                                       direction_axis=contact_env%direction_axis, &
                                       matrix_kp=rho_ao_kp(ispin, :), &
                                       index_to_cell=index_to_cell, atom_list=atom_list, &
                                       subsys=subsys, mpi_comm_global=para_env%group)
      END DO

      DEALLOCATE (index_to_cell)
      DEALLOCATE (atom_list)

      CALL timestop(handle)
   END SUBROUTINE negf_env_contact_init_matrices

! **************************************************************************************************
!> \brief Extract relevant matrix blocks for the scattering region as well as
!>        all the scattering -- contact interface regions.
!> \param negf_env            NEGF environment (modified on exit)
!> \param negf_control        NEGF control
!> \param sub_env             NEGF parallel (sub)group environment
!> \param qs_env              Primary QuickStep environment
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_env_device_init_matrices(negf_env, negf_control, sub_env, qs_env)
      TYPE(negf_env_type), INTENT(inout)                 :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_env_device_init_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icontact, ispin, nao_c, nao_s, &
                                                            ncontacts, nspins
      LOGICAL                                            :: do_kpoints
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL timeset(routineN, handle)

      IF (ALLOCATED(negf_control%atomlist_S_screening)) THEN
         CALL get_qs_env(qs_env, &
                         dft_control=dft_control, &
                         do_kpoints=do_kpoints, &
                         matrix_ks_kp=matrix_ks_kp, &
                         matrix_s_kp=matrix_s_kp, &
                         para_env=para_env, &
                         subsys=subsys)
         IF (do_kpoints) THEN
            CALL cp_abort(__LOCATION__, &
                          "K-points in device region have not been implemented yet.")
         END IF

         ncontacts = SIZE(negf_control%contacts)
         nspins = dft_control%nspins

         NULLIFY (fm_struct)
         nao_s = number_of_atomic_orbitals(subsys, negf_control%atomlist_S_screening)

         ! ++ create matrices: h_s, s_s
         NULLIFY (negf_env%s_s, fm_struct)
         ALLOCATE (negf_env%h_s(nspins))

         CALL cp_fm_struct_create(fm_struct, nrow_global=nao_s, ncol_global=nao_s, context=sub_env%blacs_env)
         CALL cp_fm_create(negf_env%s_s, fm_struct)
         DO ispin = 1, nspins
            NULLIFY (negf_env%h_s(ispin)%matrix)
            CALL cp_fm_create(negf_env%h_s(ispin)%matrix, fm_struct)
         END DO
         CALL cp_fm_struct_release(fm_struct)

         ! ++ create matrices: h_sc, s_sc
         ALLOCATE (negf_env%h_sc(nspins, ncontacts), negf_env%s_sc(ncontacts))
         DO icontact = 1, ncontacts
            nao_c = number_of_atomic_orbitals(subsys, negf_env%contacts(icontact)%atomlist_cell0)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nao_s, ncol_global=nao_c, context=sub_env%blacs_env)

            NULLIFY (negf_env%s_sc(icontact)%matrix)
            CALL cp_fm_create(negf_env%s_sc(icontact)%matrix, fm_struct)

            DO ispin = 1, nspins
               NULLIFY (negf_env%h_sc(ispin, icontact)%matrix)
               CALL cp_fm_create(negf_env%h_sc(ispin, icontact)%matrix, fm_struct)
            END DO

            CALL cp_fm_struct_release(fm_struct)
         END DO

         ! extract matrices: h_s, s_s
         DO ispin = 1, nspins
            CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_ks_kp(ispin, 1)%matrix, &
                                                  fm=negf_env%h_s(ispin)%matrix, &
                                                  atomlist_row=negf_control%atomlist_S_screening, &
                                                  atomlist_col=negf_control%atomlist_S_screening, &
                                                  subsys=subsys, mpi_comm_global=para_env%group, &
                                                  do_upper_diag=.TRUE., do_lower=.TRUE.)
         END DO

         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_s_kp(1, 1)%matrix, &
                                               fm=negf_env%s_s, &
                                               atomlist_row=negf_control%atomlist_S_screening, &
                                               atomlist_col=negf_control%atomlist_S_screening, &
                                               subsys=subsys, mpi_comm_global=para_env%group, &
                                               do_upper_diag=.TRUE., do_lower=.TRUE.)

         ! extract matrices: h_sc, s_sc
         DO icontact = 1, ncontacts
            DO ispin = 1, nspins
               CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_ks_kp(ispin, 1)%matrix, &
                                                     fm=negf_env%h_sc(ispin, icontact)%matrix, &
                                                     atomlist_row=negf_control%atomlist_S_screening, &
                                                     atomlist_col=negf_env%contacts(icontact)%atomlist_cell0, &
                                                     subsys=subsys, mpi_comm_global=para_env%group, &
                                                     do_upper_diag=.TRUE., do_lower=.TRUE.)
            END DO

            CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_s_kp(1, 1)%matrix, &
                                                  fm=negf_env%s_sc(icontact)%matrix, &
                                                  atomlist_row=negf_control%atomlist_S_screening, &
                                                  atomlist_col=negf_env%contacts(icontact)%atomlist_cell0, &
                                                  subsys=subsys, mpi_comm_global=para_env%group, &
                                                  do_upper_diag=.TRUE., do_lower=.TRUE.)
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_env_device_init_matrices

! **************************************************************************************************
!> \brief Detect the axis towards secondary unit cell.
!> \param direction_vector    direction vector
!> \param subsys_contact      QuickStep subsystem of the contact force environment
!> \param eps_geometry        accuracy in mapping atoms between different force environments
!> \retval direction_axis     direction axis: 0 (undefined), 1 (x), 2(y), 3 (z)
!> \par History
!>   * 08.2017 created [Sergey Chulkov]
! **************************************************************************************************
   FUNCTION contact_direction_axis(direction_vector, subsys_contact, eps_geometry) RESULT(direction_axis)
      REAL(kind=dp), DIMENSION(3), INTENT(in)            :: direction_vector
      TYPE(qs_subsys_type), POINTER                      :: subsys_contact
      REAL(kind=dp), INTENT(in)                          :: eps_geometry
      INTEGER                                            :: direction_axis

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contact_direction_axis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, naxes
      REAL(kind=dp), DIMENSION(3)                        :: scaled
      TYPE(cell_type), POINTER                           :: cell

      CALL qs_subsys_get(subsys_contact, cell=cell)
      CALL real_to_scaled(scaled, direction_vector, cell)

      naxes = 0

      DO i = 1, 3
         IF (ABS(scaled(i)) > eps_geometry) THEN
            direction_axis = i
            naxes = naxes+1
         END IF
      END DO

      ! direction_vector is not parallel to one of the unit cell's axis
      IF (naxes > 1) direction_axis = 0
   END FUNCTION contact_direction_axis

! **************************************************************************************************
!> \brief Estimate energy of the highest spin-alpha occupied molecular orbital.
!> \param homo_energy  HOMO energy (initialised on exit)
!> \param qs_env       QuickStep environment
!> \par History
!>   * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_homo_energy_estimate(homo_energy, qs_env)
      REAL(kind=dp), INTENT(out)                         :: homo_energy
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_homo_energy_estimate', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: gamma_point = 1

      INTEGER                                            :: handle, homo, ikpgr, ikpoint, imo, &
                                                            ispin, kplocal, nmo, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: do_kpoints
      REAL(kind=dp)                                      :: my_homo_energy
      REAL(kind=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_kp
      TYPE(kpoint_env_p_type), DIMENSION(:), POINTER     :: kp_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(mo_set_p_type), DIMENSION(:, :), POINTER      :: mos_kp

      CALL timeset(routineN, handle)
      my_homo_energy = 0.0_dp

      CALL get_qs_env(qs_env, para_env=para_env, mos=mos, kpoints=kpoints, do_kpoints=do_kpoints)

      IF (do_kpoints) THEN
         CALL get_kpoint_info(kpoints, kp_env=kp_env, kp_range=kp_range, para_env_kp=para_env_kp)

         ! looking for a processor that holds the gamma point
         IF (para_env_kp%mepos == 0 .AND. kp_range(1) <= gamma_point .AND. kp_range(2) >= gamma_point) THEN
            kplocal = kp_range(2)-kp_range(1)+1

            DO ikpgr = 1, kplocal
               CALL get_kpoint_env(kp_env(ikpgr)%kpoint_env, nkpoint=ikpoint, mos=mos_kp)

               IF (ikpoint == gamma_point) THEN
                  ! mos_kp(component, spin), where component = 1 (real), or 2 (imaginary)
                  CALL get_mo_set(mos_kp(1, 1)%mo_set, homo=homo, eigenvalues=eigenvalues) ! mu=fermi_level

                  my_homo_energy = eigenvalues(homo)
                  EXIT
               END IF
            END DO
         END IF

         CALL mp_sum(my_homo_energy, para_env%group)
      ELSE
         ! Hamiltonian of the bulk contact region has been computed without k-points.
         ! Try to obtain the HOMO energy assuming there is no OT. We probably should abort here
         ! as we do need a second replica of the bulk contact unit cell along transport
         ! direction anyway which is not available without k-points.

         CPWARN("It is strongly advised to use k-points within all contact FORCE_EVAL-s")

         nspins = SIZE(mos)

         spin_loop: DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin)%mo_set, homo=homo, nmo=nmo, eigenvalues=eigenvalues)

            DO imo = nmo, 1, -1
               IF (eigenvalues(imo) /= 0.0_dp) EXIT spin_loop
            END DO
         END DO spin_loop

         IF (imo == 0) THEN
            CPABORT("Orbital transformation (OT) for contact FORCE_EVAL-s is not supported")
         END IF

         my_homo_energy = eigenvalues(homo)
      END IF

      homo_energy = my_homo_energy
      CALL timestop(handle)
   END SUBROUTINE negf_homo_energy_estimate

! **************************************************************************************************
!> \brief List atoms from the contact's primary unit cell.
!> \param atomlist_cell0    list of atoms belonging to the contact's primary unit cell
!>                          (allocate and initialised on exit)
!> \param atom_map_cell0    atomic map of atoms from 'atomlist_cell0' (allocate and initialised on exit)
!> \param atomlist_bulk     list of atoms belonging to the bulk contact region
!> \param atom_map          atomic map of atoms from 'atomlist_bulk'
!> \param origin            origin of the contact
!> \param direction_vector  direction vector of the contact
!> \param direction_axis    axis towards secondary unit cell
!> \param subsys_device     QuickStep subsystem of the device force environment
!> \par History
!>   * 08.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE list_atoms_in_bulk_primary_unit_cell(atomlist_cell0, atom_map_cell0, atomlist_bulk, atom_map, &
                                                   origin, direction_vector, direction_axis, subsys_device)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(inout)  :: atomlist_cell0
      TYPE(negf_atom_map_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(inout)                     :: atom_map_cell0
      INTEGER, DIMENSION(:), INTENT(in)                  :: atomlist_bulk
      TYPE(negf_atom_map_type), DIMENSION(:), INTENT(in) :: atom_map
      REAL(kind=dp), DIMENSION(3), INTENT(in)            :: origin, direction_vector
      INTEGER, INTENT(in)                                :: direction_axis
      TYPE(qs_subsys_type), POINTER                      :: subsys_device

      CHARACTER(LEN=*), PARAMETER :: routineN = 'list_atoms_in_bulk_primary_unit_cell', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_min, dir_axis_min, handle, iatom, &
                                                            natoms_bulk, natoms_cell0
      REAL(kind=dp)                                      :: proj, proj_min
      REAL(kind=dp), DIMENSION(3)                        :: vector
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      CALL qs_subsys_get(subsys_device, particle_set=particle_set)

      natoms_bulk = SIZE(atomlist_bulk)
      CPASSERT(SIZE(atom_map, 1) == natoms_bulk)

      ! looking for the nearest atom from the scattering region
      proj_min = 1.0_dp
      atom_min = 1
      DO iatom = 1, natoms_bulk
         vector = particle_set(atomlist_bulk(iatom))%r-origin
         proj = projection_on_direction_vector(vector, direction_vector)

         IF (proj < proj_min) THEN
            proj_min = proj
            atom_min = iatom
         END IF
      END DO

      dir_axis_min = atom_map(atom_min)%cell(direction_axis)

      natoms_cell0 = 0
      DO iatom = 1, natoms_bulk
         IF (atom_map(iatom)%cell(direction_axis) == dir_axis_min) &
            natoms_cell0 = natoms_cell0+1
      END DO

      ALLOCATE (atomlist_cell0(natoms_cell0))
      ALLOCATE (atom_map_cell0(natoms_cell0))

      natoms_cell0 = 0
      DO iatom = 1, natoms_bulk
         IF (atom_map(iatom)%cell(direction_axis) == dir_axis_min) THEN
            natoms_cell0 = natoms_cell0+1
            atomlist_cell0(natoms_cell0) = atomlist_bulk(iatom)
            atom_map_cell0(natoms_cell0) = atom_map(iatom)
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE list_atoms_in_bulk_primary_unit_cell

! **************************************************************************************************
!> \brief List atoms from the contact's secondary unit cell.
!> \param atomlist_cell1    list of atoms belonging to the contact's secondary unit cell
!>                          (allocate and initialised on exit)
!> \param atomlist_bulk     list of atoms belonging to the bulk contact region
!> \param atom_map          atomic map of atoms from 'atomlist_bulk'
!> \param origin            origin of the contact
!> \param direction_vector  direction vector of the contact
!> \param direction_axis    axis towards the secondary unit cell
!> \param subsys_device     QuickStep subsystem of the device force environment
!> \par History
!>   * 11.2017 created [Sergey Chulkov]
!> \note Cloned from list_atoms_in_bulk_primary_unit_cell. Will be removed once we can managed to
!>        maintain consistency between real-space matrices from different force_eval sections.
! **************************************************************************************************
   SUBROUTINE list_atoms_in_bulk_secondary_unit_cell(atomlist_cell1, atomlist_bulk, atom_map, &
                                                     origin, direction_vector, direction_axis, subsys_device)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(inout)  :: atomlist_cell1
      INTEGER, DIMENSION(:), INTENT(in)                  :: atomlist_bulk
      TYPE(negf_atom_map_type), DIMENSION(:), INTENT(in) :: atom_map
      REAL(kind=dp), DIMENSION(3), INTENT(in)            :: origin, direction_vector
      INTEGER, INTENT(in)                                :: direction_axis
      TYPE(qs_subsys_type), POINTER                      :: subsys_device

      CHARACTER(LEN=*), PARAMETER :: routineN = 'list_atoms_in_bulk_secondary_unit_cell', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_max, atom_min, dir_axis_min, &
                                                            handle, iatom, natoms_bulk, &
                                                            natoms_cell1, offset
      REAL(kind=dp)                                      :: proj, proj_max, proj_min
      REAL(kind=dp), DIMENSION(3)                        :: vector
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      CALL qs_subsys_get(subsys_device, particle_set=particle_set)

      natoms_bulk = SIZE(atomlist_bulk)
      CPASSERT(SIZE(atom_map, 1) == natoms_bulk)

      ! looking for the nearest atom from the scattering region
      proj_min = 1.0_dp
      atom_min = 1
      proj_max = 0.0_dp
      atom_max = 1
      DO iatom = 1, natoms_bulk
         vector = particle_set(atomlist_bulk(iatom))%r-origin
         proj = projection_on_direction_vector(vector, direction_vector)

         IF (proj < proj_min) THEN
            proj_min = proj
            atom_min = iatom
         END IF

         IF (proj > proj_max) THEN
            proj_max = proj
            atom_max = iatom
         END IF
      END DO

      dir_axis_min = atom_map(atom_min)%cell(direction_axis)

      IF (atom_max < atom_min) THEN
         offset = -1
      ELSE
         offset = 1
      END IF

      natoms_cell1 = 0
      DO iatom = 1, natoms_bulk
         IF (atom_map(iatom)%cell(direction_axis) == dir_axis_min+offset) &
            natoms_cell1 = natoms_cell1+1
      END DO

      ALLOCATE (atomlist_cell1(natoms_cell1))

      natoms_cell1 = 0
      DO iatom = 1, natoms_bulk
         IF (atom_map(iatom)%cell(direction_axis) == dir_axis_min+offset) THEN
            natoms_cell1 = natoms_cell1+1
            atomlist_cell1(natoms_cell1) = atomlist_bulk(iatom)
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE list_atoms_in_bulk_secondary_unit_cell

! **************************************************************************************************
!> \brief Release a NEGF environment variable.
!> \param negf_env  NEGF environment to release
!> \par History
!>   * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_env_release(negf_env)
      TYPE(negf_env_type), INTENT(inout)                 :: negf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icontact, ispin, nspins

      CALL timeset(routineN, handle)

      IF (ALLOCATED(negf_env%contacts)) THEN
         DO icontact = SIZE(negf_env%contacts), 1, -1
            CALL negf_env_contact_release(negf_env%contacts(icontact))
         END DO

         DEALLOCATE (negf_env%contacts)
      END IF

      ! h_s
      IF (ALLOCATED(negf_env%h_s)) THEN
         DO ispin = SIZE(negf_env%h_s), 1, -1
            IF (ASSOCIATED(negf_env%h_s(ispin)%matrix)) &
               CALL cp_fm_release(negf_env%h_s(ispin)%matrix)
         END DO
         DEALLOCATE (negf_env%h_s)
      END IF

      ! h_sc
      IF (ALLOCATED(negf_env%h_sc)) THEN
         nspins = SIZE(negf_env%h_sc, 1)
         DO icontact = SIZE(negf_env%h_sc, 2), 1, -1
            DO ispin = nspins, 1, -1
               IF (ASSOCIATED(negf_env%h_sc(ispin, icontact)%matrix)) &
                  CALL cp_fm_release(negf_env%h_sc(ispin, icontact)%matrix)
            END DO
         END DO
         DEALLOCATE (negf_env%h_sc)
      END IF

      ! s_s
      IF (ASSOCIATED(negf_env%s_s)) &
         CALL cp_fm_release(negf_env%s_s)

      ! s_sc
      IF (ALLOCATED(negf_env%s_sc)) THEN
         DO icontact = SIZE(negf_env%s_sc), 1, -1
            IF (ASSOCIATED(negf_env%s_sc(icontact)%matrix)) &
               CALL cp_fm_release(negf_env%s_sc(icontact)%matrix)
         END DO
         DEALLOCATE (negf_env%s_sc)
      END IF

      ! density mixing
      CALL mixing_storage_release(negf_env%mixing_storage)

      CALL timestop(handle)
   END SUBROUTINE negf_env_release

! **************************************************************************************************
!> \brief Release a NEGF contact environment variable.
!> \param contact_env  NEGF contact environment to release
! **************************************************************************************************
   SUBROUTINE negf_env_contact_release(contact_env)
      TYPE(negf_env_contact_type), INTENT(inout)         :: contact_env

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_contact_release', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      ! h_00
      IF (ALLOCATED(contact_env%h_00)) THEN
         DO ispin = SIZE(contact_env%h_00), 1, -1
            IF (ASSOCIATED(contact_env%h_00(ispin)%matrix)) &
               CALL cp_fm_release(contact_env%h_00(ispin)%matrix)
         END DO
         DEALLOCATE (contact_env%h_00)
      END IF

      ! h_01
      IF (ALLOCATED(contact_env%h_01)) THEN
         DO ispin = SIZE(contact_env%h_01), 1, -1
            IF (ASSOCIATED(contact_env%h_01(ispin)%matrix)) &
               CALL cp_fm_release(contact_env%h_01(ispin)%matrix)
         END DO
         DEALLOCATE (contact_env%h_01)
      END IF

      ! rho_00
      IF (ALLOCATED(contact_env%rho_00)) THEN
         DO ispin = SIZE(contact_env%rho_00), 1, -1
            IF (ASSOCIATED(contact_env%rho_00(ispin)%matrix)) &
               CALL cp_fm_release(contact_env%rho_00(ispin)%matrix)
         END DO
         DEALLOCATE (contact_env%rho_00)
      END IF

      ! rho_01
      IF (ALLOCATED(contact_env%rho_01)) THEN
         DO ispin = SIZE(contact_env%rho_01), 1, -1
            IF (ASSOCIATED(contact_env%rho_01(ispin)%matrix)) &
               CALL cp_fm_release(contact_env%rho_01(ispin)%matrix)
         END DO
         DEALLOCATE (contact_env%rho_01)
      END IF

      ! s_00
      IF (ASSOCIATED(contact_env%s_00)) &
         CALL cp_fm_release(contact_env%s_00)

      ! s_01
      IF (ASSOCIATED(contact_env%s_01)) &
         CALL cp_fm_release(contact_env%s_01)

      IF (ALLOCATED(contact_env%atomlist_cell0)) DEALLOCATE (contact_env%atomlist_cell0)
      IF (ALLOCATED(contact_env%atomlist_cell1)) DEALLOCATE (contact_env%atomlist_cell1)
      IF (ALLOCATED(contact_env%atom_map_cell0)) DEALLOCATE (contact_env%atom_map_cell0)

      CALL timestop(handle)
   END SUBROUTINE negf_env_contact_release
END MODULE negf_env_types
