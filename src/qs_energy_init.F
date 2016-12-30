!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility subroutine for qs energy calculation
!> \par History
!>      11.2016 split out from qs_energy_utils
!> \author MK (29.10.2002)
! **************************************************************************************************
MODULE qs_energy_init
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_copy,&
                                              dbcsr_init,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE efield_utils,                    ONLY: calculate_ecore_efield
   USE input_constants,                 ONLY: kg_tnadd_atomic,&
                                              kg_tnadd_embed
   USE input_section_types,             ONLY: section_vals_type
   USE kg_environment,                  ONLY: kg_build_neighborlist,&
                                              kg_build_subsets
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_type
   USE lri_environment_methods,         ONLY: build_lri_matrices
   USE lri_environment_types,           ONLY: lri_environment_type
   USE molecule_types_new,              ONLY: molecule_of_atom,&
                                              molecule_type
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
   USE qs_dftb_dispersion,              ONLY: calculate_dftb_dispersion
   USE qs_dftb_matrices,                ONLY: build_dftb_matrices
   USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_external_density,             ONLY: external_read_density
   USE qs_external_potential,           ONLY: external_c_potential,&
                                              external_e_potential
   USE qs_ks_methods,                   ONLY: qs_ks_allocate_basics
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
   USE qs_update_s_mstruct,             ONLY: qs_env_update_s_mstruct
   USE scptb_core_interactions,         ONLY: scptb_core_interaction
   USE scptb_core_matrix,               ONLY: build_scptb_core_matrix
   USE se_core_core,                    ONLY: se_core_core_interaction
   USE se_core_matrix,                  ONLY: build_se_core_matrix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_init'

   PUBLIC :: qs_energies_init

CONTAINS

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Driver routine for the initial
!>        setup and calculations for a qs energy calculation
!> \param qs_env ...
!> \param calc_forces ...
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE qs_energies_init(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_init', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: img, ispin, nimg, nspin
      LOGICAL                                            :: has_unit_metric, molecule_only
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_w_mp2
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_w
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      NULLIFY (ks_env, matrix_w, matrix_w_mp2, matrix_s, dft_control)

      CALL qs_energies_init_kg(qs_env, molecule_only)
      CALL qs_energies_init_hamiltonians(qs_env, calc_forces, molecule_only)
      CALL qs_ks_allocate_basics(qs_env)
      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)

      ! if need forces allocate energy weighted density matrices
      IF (calc_forces .AND. .NOT. has_unit_metric) THEN
         CALL get_qs_env(qs_env, &
                         ks_env=ks_env, &
                         matrix_s_kp=matrix_s, &
                         dft_control=dft_control)
         nspin = dft_control%nspins
         nimg = dft_control%nimages
         matrix => matrix_s(1, 1)%matrix
         CALL dbcsr_allocate_matrix_set(matrix_w, nspin, nimg)
         DO ispin = 1, nspin
            DO img = 1, nimg
               ALLOCATE (matrix_w(ispin, img)%matrix)
               CALL dbcsr_init(matrix_w(ispin, img)%matrix)
               CALL dbcsr_copy(matrix_w(ispin, img)%matrix, matrix, name="W MATRIX")
               CALL dbcsr_set(matrix_w(ispin, img)%matrix, 0.0_dp)
            END DO
         END DO
         CALL set_ks_env(ks_env, matrix_w_kp=matrix_w)
         IF (ASSOCIATED(qs_env%mp2_env)) THEN
            ! create the MP2 energy weighted density matrix
            CALL dbcsr_allocate_matrix_set(matrix_w_mp2, dft_control%nspins)
            DO ispin = 1, nspin
               ALLOCATE (matrix_w_mp2(ispin)%matrix)
               CALL dbcsr_init(matrix_w_mp2(ispin)%matrix)
               CALL dbcsr_copy(matrix_w_mp2(ispin)%matrix, matrix, &
                               name="W MATRIX MP2")
               CALL dbcsr_set(matrix_w_mp2(ispin)%matrix, 0.0_dp)
            END DO
            CALL set_ks_env(ks_env, matrix_w_mp2=matrix_w_mp2)
         END IF
      ENDIF

   END SUBROUTINE qs_energies_init

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Puts initialization of the Kim-Gordon
!>        settings into separate subroutine
!> \param qs_env ...
!> \param molecule_only ...
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE qs_energies_init_kg(qs_env, molecule_only)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: molecule_only

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_init_kg', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, isubset, natom
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set

      CALL timeset(routineN, handle)
      NULLIFY (dft_control, para_env)
      molecule_only = .FALSE.

      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env)
      IF (dft_control%qs_control%do_kg) THEN

         ! create neighbor lists with molecular blocks
         molecule_only = .TRUE.

         ! get the set of molecules
         CALL get_qs_env(qs_env=qs_env, molecule_set=molecule_set, natom=natom)

         qs_env%kg_env%natom = natom

         ! store set of molecules in kg_env
         qs_env%kg_env%molecule_set => molecule_set

         ! build the (new) full neighborlist
         CALL kg_build_neighborlist(qs_env, sab_orb=qs_env%kg_env%sab_orb_full)

         IF (.NOT. ALLOCATED(qs_env%kg_env%atom_to_molecule)) THEN
            ALLOCATE (qs_env%kg_env%atom_to_molecule(natom))

            ! get the mapping from atoms to molecules
            CALL molecule_of_atom(molecule_set, atom_to_mol=qs_env%kg_env%atom_to_molecule)
         END IF

         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN
            ! allocate the subset list
            IF (.NOT. ASSOCIATED(qs_env%kg_env%subset_of_mol)) ALLOCATE (qs_env%kg_env%subset_of_mol(SIZE(molecule_set)))

            CALL kg_build_subsets(qs_env%kg_env, para_env)

            DO isubset = 1, qs_env%kg_env%nsubsets
               ! build the (new) molecular neighborlist of the current subset
               CALL kg_build_neighborlist(qs_env, sab_orb=qs_env%kg_env%subset(isubset)%sab_orb, molecular=.TRUE., &
                                          subset_of_mol=qs_env%kg_env%subset_of_mol, current_subset=isubset)
            END DO
         ELSE IF (qs_env%kg_env%tnadd_method == kg_tnadd_atomic) THEN
            ! build the A-C list for the nonadditive kinetic energy potential
            CALL kg_build_neighborlist(qs_env, sac_kin=qs_env%kg_env%sac_kin)
         END IF

      END IF
      CALL timestop(handle)

   END SUBROUTINE qs_energies_init_kg

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Moves computation of the different
!>        core hamiltonians into separate subroutine
!> \param qs_env        QS environment
!> \param calc_forces   Calculate forces
!> \param molecule_only restrict neighbor list to molecules
!> \par History
!>      05.2013 created [Florian Schiffmann]
!>      08.2014 Kpoints [JGH]
! **************************************************************************************************

   SUBROUTINE qs_energies_init_hamiltonians(qs_env, calc_forces, molecule_only)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces
      LOGICAL                                            :: molecule_only

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_init_hamiltonians', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_kpoints
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      input=input, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      kpoints=kpoints, &
                      do_kpoints=do_kpoints)

      ! create neighbor lists for standard use in QS
      CALL build_qs_neighbor_lists(qs_env, para_env, molecular=molecule_only, &
                                   force_env_section=input)

      ! calculate cell index for k-point calculations
      IF (do_kpoints) THEN
         CALL get_qs_env(qs_env, sab_kp=sab_nl)
         CALL kpoint_init_cell_index(kpoints, sab_nl, para_env, dft_control)
      ENDIF

      dft_control%qs_control%becke_control%need_pot = .TRUE.

      ! Calculate the overlap and the core Hamiltonian integral matrix
      IF (dft_control%qs_control%semi_empirical) THEN
         CALL build_se_core_matrix(qs_env=qs_env, para_env=para_env, &
                                   calculate_forces=.FALSE.)
         CALL qs_env_update_s_mstruct(qs_env)
         CALL se_core_core_interaction(qs_env, para_env, calculate_forces=.FALSE.)
         CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env, energy=energy)
         CALL calculate_dispersion_pairpot(qs_env, dispersion_env, energy%dispersion, calc_forces)
      ELSEIF (dft_control%qs_control%dftb) THEN
         CALL build_dftb_matrices(qs_env=qs_env, para_env=para_env, &
                                  calculate_forces=.FALSE.)
         CALL calculate_dftb_dispersion(qs_env=qs_env, para_env=para_env, &
                                        calculate_forces=.FALSE.)
         CALL qs_env_update_s_mstruct(qs_env)
      ELSEIF (dft_control%qs_control%scptb) THEN
         CALL build_scptb_core_matrix(qs_env=qs_env, calculate_forces=.FALSE.)
         CALL qs_env_update_s_mstruct(qs_env)
         CALL scptb_core_interaction(qs_env, calculate_forces=.FALSE.)
         CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env, energy=energy)
         CALL calculate_dispersion_pairpot(qs_env, dispersion_env, energy%dispersion, calc_forces)
      ELSE
         CALL build_core_hamiltonian_matrix(qs_env=qs_env, calculate_forces=.FALSE.)
         CALL qs_env_update_s_mstruct(qs_env)
         CALL calculate_ecore_self(qs_env)
         CALL calculate_ecore_efield(qs_env, calculate_forces=.FALSE.)
         CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces=.FALSE.)
         !swap external_e_potential before external_c_potential, to ensure
         !that external potential on grid is loaded before calculating energy of cores
         CALL external_e_potential(qs_env)
         IF (.NOT. dft_control%qs_control%gapw) THEN
            CALL external_c_potential(qs_env, calculate_forces=.FALSE.)
         END IF
         ! LRIGPW  matrices
         IF (dft_control%qs_control%lrigpw) THEN
            CALL get_qs_env(qs_env=qs_env, lri_env=lri_env)
            CALL build_lri_matrices(lri_env, qs_env, calc_forces)
         ENDIF

         ! ZMP addition to read external density
         CALL external_read_density(qs_env)

         ! Add possible pair potential dispersion energy - Evaluate first so we can print
         ! energy info at the end of the SCF
         CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env, energy=energy)
         CALL calculate_dispersion_pairpot(qs_env, dispersion_env, energy%dispersion, calc_forces)

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_energies_init_hamiltonians

END MODULE qs_energy_init
