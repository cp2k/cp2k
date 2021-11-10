!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility subroutine for qs energy calculation
!> \par History
!>      11.2016 split out from qs_energy_utils
!> \author MK (29.10.2002)
! **************************************************************************************************
MODULE qs_energy_init
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE efield_utils,                    ONLY: calculate_ecore_efield
   USE input_constants,                 ONLY: kg_tnadd_atomic,&
                                              kg_tnadd_embed,&
                                              kg_tnadd_embed_ri,&
                                              kg_tnadd_none
   USE input_section_types,             ONLY: section_vals_type
   USE kg_environment,                  ONLY: kg_build_neighborlist,&
                                              kg_build_subsets
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_type
   USE lri_environment_methods,         ONLY: build_lri_matrices,&
                                              calculate_lri_integrals
   USE lri_environment_types,           ONLY: lri_environment_type
   USE molecule_types,                  ONLY: molecule_of_atom,&
                                              molecule_type
   USE optimize_embedding_potential,    ONLY: given_embed_pot
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
   USE qs_gcp_method,                   ONLY: calculate_gcp_pairpot
   USE qs_gcp_types,                    ONLY: qs_gcp_type
   USE qs_ks_methods,                   ONLY: qs_ks_allocate_basics
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
   USE qs_update_s_mstruct,             ONLY: qs_env_update_s_mstruct
   USE ri_environment_methods,          ONLY: build_ri_matrices
   USE se_core_core,                    ONLY: se_core_core_interaction
   USE se_core_matrix,                  ONLY: build_se_core_matrix
   USE xtb_matrices,                    ONLY: build_xtb_matrices
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

      INTEGER                                            :: img, ispin, nimg, nspin
      LOGICAL                                            :: has_unit_metric, molecule_only
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_w
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      NULLIFY (ks_env, matrix_w, matrix_s, dft_control)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%qs_control%do_kg) THEN
         molecule_only = .TRUE.
         CALL qs_energies_init_kg(qs_env)
      ELSE
         molecule_only = .FALSE.
      END IF
      CALL qs_energies_init_hamiltonians(qs_env, calc_forces, molecule_only)
      CALL qs_ks_allocate_basics(qs_env)

      ! if need forces allocate energy weighted density matrices
      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)
      IF (calc_forces .AND. .NOT. has_unit_metric) THEN
         CALL get_qs_env(qs_env, &
                         ks_env=ks_env, &
                         matrix_s_kp=matrix_s)
         nspin = dft_control%nspins
         nimg = dft_control%nimages
         matrix => matrix_s(1, 1)%matrix
         CALL dbcsr_allocate_matrix_set(matrix_w, nspin, nimg)
         DO ispin = 1, nspin
            DO img = 1, nimg
               ALLOCATE (matrix_w(ispin, img)%matrix)
               CALL dbcsr_copy(matrix_w(ispin, img)%matrix, matrix, name="W MATRIX")
               CALL dbcsr_set(matrix_w(ispin, img)%matrix, 0.0_dp)
            END DO
         END DO
         CALL set_ks_env(ks_env, matrix_w_kp=matrix_w)
      END IF

   END SUBROUTINE qs_energies_init

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Puts initialization of the Kim-Gordon
!>        settings into separate subroutine
!> \param qs_env ...
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE qs_energies_init_kg(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_init_kg'

      INTEGER                                            :: handle, isubset, natom
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: soo_list, soo_list1

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env)
      CPASSERT(dft_control%qs_control%do_kg)

      kg_env => qs_env%kg_env

      ! get the set of molecules
      CALL get_qs_env(qs_env=qs_env, molecule_set=molecule_set, natom=natom)
      kg_env%natom = natom
      ! store set of molecules in kg_env
      kg_env%molecule_set => molecule_set
      ! build the (new) full neighborlist
      CALL kg_build_neighborlist(qs_env, sab_orb=kg_env%sab_orb_full)

      IF (.NOT. ALLOCATED(kg_env%atom_to_molecule)) THEN
         ALLOCATE (kg_env%atom_to_molecule(natom))
         ! get the mapping from atoms to molecules
         CALL molecule_of_atom(molecule_set, atom_to_mol=kg_env%atom_to_molecule)
      END IF

      SELECT CASE (kg_env%tnadd_method)
      CASE (kg_tnadd_embed)
         ! allocate the subset list
         IF (.NOT. ASSOCIATED(kg_env%subset_of_mol)) THEN
            ALLOCATE (kg_env%subset_of_mol(SIZE(molecule_set)))
         END IF
         !
         CALL kg_build_subsets(kg_env, para_env)
         !
         DO isubset = 1, kg_env%nsubsets
            ! build the (new) molecular neighborlist of the current subset
            CALL kg_build_neighborlist(qs_env, sab_orb=kg_env%subset(isubset)%sab_orb, molecular=.TRUE., &
                                       subset_of_mol=kg_env%subset_of_mol, current_subset=isubset)
         END DO
      CASE (kg_tnadd_embed_ri)
!deb should be deleted as soon as atomic grids work
         ! allocate the subset list
         IF (.NOT. ASSOCIATED(kg_env%subset_of_mol)) THEN
            ALLOCATE (kg_env%subset_of_mol(SIZE(molecule_set)))
         END IF
         !
         CALL kg_build_subsets(kg_env, para_env)
         !
         DO isubset = 1, kg_env%nsubsets
            ! build the (new) molecular neighborlist of the current subset
            CALL kg_build_neighborlist(qs_env, sab_orb=kg_env%subset(isubset)%sab_orb, molecular=.TRUE., &
                                       subset_of_mol=kg_env%subset_of_mol, current_subset=isubset)
         END DO
!deb
         ! LRI neighborlist
         NULLIFY (soo_list)
         CALL kg_build_neighborlist(qs_env, sab_orb=soo_list, molecular=.TRUE.)
         kg_env%lri_env%soo_list => soo_list
         CALL calculate_lri_integrals(kg_env%lri_env, qs_env)
         IF (qs_env%energy_correction) THEN
            NULLIFY (soo_list1)
            CALL kg_build_neighborlist(qs_env, sab_orb=soo_list1, molecular=.TRUE.)
            kg_env%lri_env1%soo_list => soo_list1
            CALL calculate_lri_integrals(kg_env%lri_env1, qs_env)
         END IF

         ! Atomic grids
      CASE (kg_tnadd_atomic)
         ! build the A-C list for the nonadditive kinetic energy potential
         CALL kg_build_neighborlist(qs_env, sac_kin=kg_env%sac_kin)
      CASE (kg_tnadd_none)
         ! nothing to do
      CASE DEFAULT
         CPABORT("KG:TNADD METHOD")
      END SELECT

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

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_init_hamiltonians'

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
      TYPE(qs_gcp_type), POINTER                         :: gcp_env
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
      END IF
      IF (dft_control%qs_control%cdft) THEN
         IF (.NOT. (dft_control%qs_control%cdft_control%external_control)) &
            dft_control%qs_control%cdft_control%need_pot = .TRUE.
         IF (ASSOCIATED(dft_control%qs_control%cdft_control%group)) THEN
            ! In case CDFT weight function was built beforehand (in mixed force_eval)
            IF (ASSOCIATED(dft_control%qs_control%cdft_control%group(1)%weight%pw)) &
               dft_control%qs_control%cdft_control%need_pot = .FALSE.
         END IF
      END IF

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
      ELSEIF (dft_control%qs_control%xtb) THEN
         CALL build_xtb_matrices(qs_env=qs_env, para_env=para_env, &
                                 calculate_forces=.FALSE.)
         CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env, energy=energy)
         CALL calculate_dispersion_pairpot(qs_env, dispersion_env, energy%dispersion, calc_forces)
         CALL qs_env_update_s_mstruct(qs_env)
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
         ! LRIGPW/RIGPW  matrices
         IF (dft_control%qs_control%lrigpw) THEN
            CALL get_qs_env(qs_env=qs_env, lri_env=lri_env)
            CALL build_lri_matrices(lri_env, qs_env)
         ELSE IF (dft_control%qs_control%rigpw) THEN
            CALL get_qs_env(qs_env=qs_env, lri_env=lri_env)
            CALL build_ri_matrices(lri_env, qs_env, calculate_forces=.FALSE.)
         END IF

         ! ZMP addition to read external density
         CALL external_read_density(qs_env)

         ! Add possible pair potential dispersion energy - Evaluate first so we can print
         ! energy info at the end of the SCF
         CALL get_qs_env(qs_env=qs_env, dispersion_env=dispersion_env, energy=energy)
         CALL calculate_dispersion_pairpot(qs_env, dispersion_env, energy%dispersion, calc_forces)
         ! Add possible pair potential gCP energy - Evaluate first so we can print
         ! energy info at the end of the SCF
         CALL get_qs_env(qs_env=qs_env, gcp_env=gcp_env, energy=energy)
         IF (ASSOCIATED(gcp_env)) THEN
            CALL calculate_gcp_pairpot(qs_env, gcp_env, energy%gcp, calc_forces)
         END IF

      END IF
      ! Embedding potential
      IF (dft_control%qs_control%dfet_embedded) THEN
         dft_control%apply_embed_pot = .TRUE.
         CALL given_embed_pot(qs_env)
      END IF
      ! Matrix embedding potential
      IF (dft_control%qs_control%dmfet_embedded) THEN
         dft_control%apply_dmfet_pot = .TRUE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_energies_init_hamiltonians

END MODULE qs_energy_init
