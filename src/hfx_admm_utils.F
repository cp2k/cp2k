!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for hfx and admm methods
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!>     Made GAPW compatible 12.2019 (A. Bussy)
!> \author MI
! **************************************************************************************************
MODULE hfx_admm_utils
   USE admm_dm_types,                   ONLY: admm_dm_create
   USE admm_methods,                    ONLY: scale_dm
   USE admm_types,                      ONLY: admm_env_create,&
                                              admm_gapw_type,&
                                              admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: copy_gto_basis_set,&
                                              get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              plane_distance
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_m_by_n_from_row_template,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE external_potential_types,        ONLY: copy_potential
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_pw_methods,                  ONLY: pw_hfx
   USE hfx_ri,                          ONLY: hfx_ri_update_forces,&
                                              hfx_ri_update_ks
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: &
        do_admm_aux_exch_func_bee, do_admm_aux_exch_func_bee_libxc, do_admm_aux_exch_func_default, &
        do_admm_aux_exch_func_default_libxc, do_admm_aux_exch_func_none, &
        do_admm_aux_exch_func_opt, do_admm_aux_exch_func_opt_libxc, do_admm_aux_exch_func_pbex, &
        do_admm_aux_exch_func_pbex_libxc, do_admm_aux_exch_func_sx_libxc, do_potential_coulomb, &
        do_potential_long, do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_short, &
        do_potential_truncated, xc_funct_no_shortcut, xc_none
   USE input_section_types,             ONLY: section_vals_duplicate,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_remove_values,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                              paw_proj_set_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_interactions,                 ONLY: init_interaction_radii
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              init_gapw_basis_set,&
                                              init_gapw_nlcc,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_set_create
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup,&
                                              write_neighbor_lists
   USE qs_oce_methods,                  ONLY: build_oce_matrices
   USE qs_oce_types,                    ONLY: allocate_oce_set,&
                                              create_oce_set
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_atom_methods,             ONLY: init_rho_atom
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_type
   USE rt_propagation_types,            ONLY: rt_prop_type
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list
   USE virial_types,                    ONLY: virial_type
   USE xc_adiabatic_utils,              ONLY: rescale_xc_potential
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! *** Public subroutines ***
   PUBLIC :: hfx_ks_matrix, hfx_admm_init, create_admm_xc_section, tddft_hfx_matrix

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_admm_utils'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE hfx_admm_init(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_admm_init'

      INTEGER                                            :: handle, ispin, n_rep_hf, nao_aux_fit, &
                                                            natoms, nelectron, nmo
      LOGICAL                                            :: s_mstruct_changed, use_virial
      REAL(dp)                                           :: maxocc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: aux_fit_fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_aux_fit
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos, mos_aux_fit
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input, xc_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (admm_env, hfx_sections, mos, mos_aux_fit, para_env, virial, &
               mo_coeff_aux_fit, xc_section, rho, ks_env, dft_control, input, &
               qs_kind_set, mo_coeff_b, aux_fit_fm_struct, blacs_env)

      CALL get_qs_env(qs_env, &
                      mos=mos, &
                      admm_env=admm_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      s_mstruct_changed=s_mstruct_changed, &
                      rho=rho, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      input=input, &
                      virial=virial)

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      IF (n_rep_hf > 1) &
         CPABORT("ADMM can handle only one HF section.")

      IF (.NOT. ASSOCIATED(admm_env)) THEN
         ! setup admm environment
         CALL get_qs_env(qs_env, input=input, natom=natoms, qs_kind_set=qs_kind_set)
         CALL get_qs_kind_set(qs_kind_set, nsgf=nao_aux_fit, basis_type="AUX_FIT")
         CALL admm_env_create(admm_env, dft_control%admm_control, mos, para_env, natoms, nao_aux_fit)
         CALL set_qs_env(qs_env, admm_env=admm_env)
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
         CALL create_admm_xc_section(x_data=qs_env%x_data, xc_section=xc_section, &
                                     admm_env=admm_env)

         ! Initialize the GAPW data types
         IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) &
            CALL init_admm_gapw(qs_env)

         ! ADMM neighbor lists and overlap matrices
         CALL admm_init_hamiltonians(admm_env, qs_env)

         !The aux_fit task list and densities
         ALLOCATE (admm_env%rho_aux_fit)
         CALL qs_rho_create(admm_env%rho_aux_fit)
         ALLOCATE (admm_env%rho_aux_fit_buffer)
         CALL qs_rho_create(admm_env%rho_aux_fit_buffer)
         CALL admm_update_s_mstruct(admm_env, qs_env)
         IF (admm_env%do_gapw) CALL update_admm_gapw(qs_env)

         !The ADMM KS matrices
         CALL admm_alloc_ks_matrices(admm_env, qs_env)

         !The aux_fit MOs and derivatives
         ALLOCATE (mos_aux_fit(dft_control%nspins))
         DO ispin = 1, dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin), nmo=nmo, nelectron=nelectron, maxocc=maxocc)
            CALL allocate_mo_set(mo_set=mos_aux_fit(ispin), &
                                 nao=nao_aux_fit, &
                                 nmo=nmo, &
                                 nelectron=nelectron, &
                                 n_el_f=REAL(nelectron, dp), &
                                 maxocc=maxocc, &
                                 flexible_electron_count=dft_control%relax_multiplicity)
         END DO
         admm_env%mos_aux_fit => mos_aux_fit

         DO ispin = 1, dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin), nmo=nmo)
            CALL cp_fm_struct_create(aux_fit_fm_struct, context=blacs_env, para_env=para_env, &
                                     nrow_global=nao_aux_fit, ncol_global=nmo)
            CALL get_mo_set(mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit, mo_coeff_b=mo_coeff_b)
            IF (.NOT. ASSOCIATED(mo_coeff_aux_fit)) THEN
               CALL init_mo_set(mos_aux_fit(ispin), fm_struct=aux_fit_fm_struct, &
                                name="qs_env%mo_aux_fit"//TRIM(ADJUSTL(cp_to_string(ispin))))
            END IF
            CALL cp_fm_struct_release(aux_fit_fm_struct)

            IF (.NOT. ASSOCIATED(mo_coeff_b)) THEN
               CALL cp_fm_get_info(mos_aux_fit(ispin)%mo_coeff, ncol_global=nmo)
               CALL dbcsr_init_p(mos_aux_fit(ispin)%mo_coeff_b)
               CALL cp_dbcsr_m_by_n_from_row_template(mos_aux_fit(ispin)%mo_coeff_b, &
                                                      template=admm_env%matrix_s_aux_fit(1)%matrix, &
                                                      n=nmo, sym=dbcsr_type_no_symmetry)
            END IF
         END DO

         IF (qs_env%requires_mo_derivs) THEN
            ALLOCATE (admm_env%mo_derivs_aux_fit(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               CALL get_mo_set(admm_env%mos_aux_fit(ispin), mo_coeff=mo_coeff_aux_fit)
               CALL cp_fm_create(admm_env%mo_derivs_aux_fit(ispin), mo_coeff_aux_fit%matrix_struct)
            END DO
         END IF

      ELSE IF (s_mstruct_changed) THEN
         CALL admm_init_hamiltonians(admm_env, qs_env)
         CALL admm_update_s_mstruct(admm_env, qs_env)
         CALL admm_alloc_ks_matrices(admm_env, qs_env)
         IF (admm_env%do_gapw) CALL update_admm_gapw(qs_env)
      END IF

      IF (admm_env%do_gapw .AND. dft_control%do_admm_dm) THEN
         CPABORT("GAPW ADMM not implemented for MCWEENY or NONE_DM purification.")
      END IF

      !ADMMS and ADMMP stress tensors only available for close-shell systesms, because virial cannot
      !be scaled by gsi spin component wise
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (use_virial .AND. admm_env%do_admms .AND. dft_control%nspins == 2) THEN
         CPABORT("ADMMS stress tensor is only available for closed-shell systems")
      END IF
      IF (use_virial .AND. admm_env%do_admmp .AND. dft_control%nspins == 2) THEN
         CPABORT("ADMMP stress tensor is only available for closed-shell systems")
      END IF

      IF (dft_control%do_admm_dm .AND. .NOT. ASSOCIATED(admm_env%admm_dm)) THEN
         CALL admm_dm_create(admm_env%admm_dm, dft_control%admm_control, nspins=dft_control%nspins, natoms=natoms)
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_admm_init

! **************************************************************************************************
!> \brief Sets up the admm_gapw env
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE init_admm_gapw(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: ikind, nkind
      TYPE(admm_gapw_type), POINTER                      :: admm_gapw_env
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: aux_fit_basis, aux_fit_soft_basis, &
                                                            orb_basis, soft_basis
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: admm_kind_set, qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input

      NULLIFY (admm_kind_set, aux_fit_basis, atomic_kind_set, aux_fit_soft_basis, &
               dft_control, input, orb_basis, para_env, qs_kind_set, soft_basis)

      CALL get_qs_env(qs_env, admm_env=admm_env, &
                      atomic_kind_set=atomic_kind_set, &
                      dft_control=dft_control, &
                      input=input, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set)

      admm_env%do_gapw = .TRUE.
      ALLOCATE (admm_env%admm_gapw_env)
      admm_gapw_env => admm_env%admm_gapw_env
      NULLIFY (admm_gapw_env%local_rho_set)
      NULLIFY (admm_gapw_env%admm_kind_set)
      NULLIFY (admm_gapw_env%task_list)

      !Create a new kind set for the ADMM stuff (paw_proj soft AUX_FIT basis, etc)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (admm_gapw_env%admm_kind_set(nkind))
      admm_kind_set => admm_gapw_env%admm_kind_set

      !In this new kind set, we want the AUX_FIT basis to be known as ORB, such that GAPW routines work
      DO ikind = 1, nkind
         !copying over simple data  of interest from qs_kind_set
         admm_kind_set(ikind)%name = qs_kind_set(ikind)%name
         admm_kind_set(ikind)%element_symbol = qs_kind_set(ikind)%element_symbol
         admm_kind_set(ikind)%natom = qs_kind_set(ikind)%natom
         admm_kind_set(ikind)%hard_radius = qs_kind_set(ikind)%hard_radius
         admm_kind_set(ikind)%max_rad_local = qs_kind_set(ikind)%max_rad_local
         admm_kind_set(ikind)%gpw_type_forced = qs_kind_set(ikind)%gpw_type_forced
         admm_kind_set(ikind)%ngrid_rad = qs_kind_set(ikind)%ngrid_rad
         admm_kind_set(ikind)%ngrid_ang = qs_kind_set(ikind)%ngrid_ang

         !copying potentials of interest from qs_kind_set
         IF (ASSOCIATED(qs_kind_set(ikind)%all_potential)) THEN
            CALL copy_potential(qs_kind_set(ikind)%all_potential, admm_kind_set(ikind)%all_potential)
         END IF
         IF (ASSOCIATED(qs_kind_set(ikind)%gth_potential)) THEN
            CALL copy_potential(qs_kind_set(ikind)%gth_potential, admm_kind_set(ikind)%gth_potential)
         END IF
         IF (ASSOCIATED(qs_kind_set(ikind)%sgp_potential)) THEN
            CALL copy_potential(qs_kind_set(ikind)%sgp_potential, admm_kind_set(ikind)%sgp_potential)
         END IF

         NULLIFY (orb_basis)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=aux_fit_basis, basis_type="AUX_FIT")
         CALL copy_gto_basis_set(aux_fit_basis, orb_basis)
         CALL add_basis_set_to_container(admm_kind_set(ikind)%basis_sets, orb_basis, "ORB")
      END DO

      !Create the corresponding soft basis set (and projectors)
      CALL init_gapw_basis_set(admm_kind_set, dft_control%qs_control, input, &
                               modify_qs_control=.FALSE.)

      !Make sure the basis and the projectors are well initialized
      CALL init_interaction_radii(dft_control%qs_control, admm_kind_set)

      !We also init the atomic grids and harmonics
      CALL local_rho_set_create(admm_gapw_env%local_rho_set)
      CALL init_rho_atom(admm_gapw_env%local_rho_set%rho_atom_set, &
                         atomic_kind_set, admm_kind_set, dft_control, para_env)

      !Make sure that any NLCC potential is well initialized
      CALL init_gapw_nlcc(admm_kind_set)

      !Need to have access to the soft AUX_FIT basis from the qs_env => add it to the qs_kinds
      DO ikind = 1, nkind
         NULLIFY (aux_fit_soft_basis)
         CALL get_qs_kind(admm_kind_set(ikind), basis_set=soft_basis, basis_type="ORB_SOFT")
         CALL copy_gto_basis_set(soft_basis, aux_fit_soft_basis)
         CALL add_basis_set_to_container(qs_kind_set(ikind)%basis_sets, aux_fit_soft_basis, "AUX_FIT_SOFT")
      END DO

   END SUBROUTINE init_admm_gapw

! **************************************************************************************************
!> \brief Builds the ADMM nmeighbor lists and overlap matrix on the model of qs_energies_init_hamiltonians()
!> \param admm_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE admm_init_hamiltonians(admm_env, qs_env)

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'admm_init_hamiltonians'

      INTEGER                                            :: handle, ikind, nkind
      LOGICAL                                            :: mic, molecule_only
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: aux_fit_present, orb_present
      REAL(dp)                                           :: pdist, subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: aux_fit_radius, orb_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(gto_basis_set_type), POINTER                  :: aux_fit_basis_set, orb_basis_set
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: neighbor_list_section

      NULLIFY (particle_set, cell, kpoints, distribution_1d, distribution_2d, molecule_set, &
               atomic_kind_set, dft_control, neighbor_list_section, aux_fit_basis_set, orb_basis_set, &
               ks_env, para_env, qs_kind_set)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, particle_set=particle_set, cell=cell, kpoints=kpoints, &
                      local_particles=distribution_1d, distribution_2d=distribution_2d, &
                      molecule_set=molecule_set, atomic_kind_set=atomic_kind_set, &
                      dft_control=dft_control, para_env=para_env, qs_kind_set=qs_kind_set)
      ALLOCATE (orb_present(nkind), aux_fit_present(nkind))
      ALLOCATE (orb_radius(nkind), aux_fit_radius(nkind), pair_radius(nkind, nkind))
      aux_fit_radius(:) = 0.0_dp

      molecule_only = .FALSE.
      IF (dft_control%qs_control%do_kg) molecule_only = .TRUE.
      mic = molecule_only
      IF (kpoints%nkp > 0) THEN
         mic = .FALSE.
      ELSEIF (dft_control%qs_control%semi_empirical) THEN
         mic = .TRUE.
      END IF

      pdist = dft_control%qs_control%pairlist_radius

      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)
      neighbor_list_section => section_vals_get_subs_vals(qs_env%input, "DFT%PRINT%NEIGHBOR_LISTS")

      ALLOCATE (atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, molecule_only, particle_set=particle_set)

      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
         IF (ASSOCIATED(orb_basis_set)) THEN
            orb_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, kind_radius=orb_radius(ikind))
         ELSE
            orb_present(ikind) = .FALSE.
         END IF

         CALL get_qs_kind(qs_kind_set(ikind), basis_set=aux_fit_basis_set, basis_type="AUX_FIT")
         IF (ASSOCIATED(aux_fit_basis_set)) THEN
            aux_fit_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=aux_fit_basis_set, kind_radius=aux_fit_radius(ikind))
         ELSE
            aux_fit_present(ikind) = .FALSE.
         END IF
      END DO

      IF (pdist < 0.0_dp) THEN
         pdist = MAX(plane_distance(1, 0, 0, cell), &
                     plane_distance(0, 1, 0, cell), &
                     plane_distance(0, 0, 1, cell))
      END IF
      CALL pair_radius_setup(aux_fit_present, aux_fit_present, aux_fit_radius, aux_fit_radius, pair_radius, pdist)
      CALL build_neighbor_lists(admm_env%sab_aux_fit, particle_set, atom2d, cell, pair_radius, &
                                mic=mic, molecular=molecule_only, subcells=subcells, nlname="sab_aux_fit")
      CALL build_neighbor_lists(admm_env%sab_aux_fit_asymm, particle_set, atom2d, cell, pair_radius, &
                                mic=mic, symmetric=.FALSE., molecular=molecule_only, subcells=subcells, &
                                nlname="sab_aux_fit_asymm")
      CALL pair_radius_setup(aux_fit_present, orb_present, aux_fit_radius, orb_radius, pair_radius)
      CALL build_neighbor_lists(admm_env%sab_aux_fit_vs_orb, particle_set, atom2d, cell, pair_radius, &
                                mic=mic, symmetric=.FALSE., molecular=molecule_only, subcells=subcells, &
                                nlname="sab_aux_fit_vs_orb")

      CALL write_neighbor_lists(admm_env%sab_aux_fit, particle_set, cell, para_env, neighbor_list_section, &
                                "/SAB_AUX_FIT", "sab_aux_fit", "AUX_FIT_ORBITAL AUX_FIT_ORBITAL")
      CALL write_neighbor_lists(admm_env%sab_aux_fit_vs_orb, particle_set, cell, para_env, neighbor_list_section, &
                                "/SAB_AUX_FIT_VS_ORB", "sab_aux_fit_vs_orb", "ORBITAL AUX_FIT_ORBITAL")

      CALL atom2d_cleanup(atom2d)

      !The ADMM overlap matrices (initially in qs_core_hamiltonian.F)
      CALL get_qs_env(qs_env, ks_env=ks_env)

      IF (ASSOCIATED(admm_env%matrix_s_aux_fit)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_s_aux_fit)
      CALL build_overlap_matrix(ks_env, matrix_s=admm_env%matrix_s_aux_fit, &
                                matrix_name="AUX_FIT_OVERLAP", &
                                basis_type_a="AUX_FIT", &
                                basis_type_b="AUX_FIT", &
                                sab_nl=admm_env%sab_aux_fit)
      IF (ASSOCIATED(admm_env%matrix_s_aux_fit_vs_orb)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_s_aux_fit_vs_orb)
      CALL build_overlap_matrix(ks_env, matrix_s=admm_env%matrix_s_aux_fit_vs_orb, &
                                matrix_name="MIXED_OVERLAP", &
                                basis_type_a="AUX_FIT", &
                                basis_type_b="ORB", &
                                sab_nl=admm_env%sab_aux_fit_vs_orb)

      CALL timestop(handle)

   END SUBROUTINE admm_init_hamiltonians

! **************************************************************************************************
!> \brief Updates the ADMM task_list and density based on the model of qs_env_update_s_mstruct()
!> \param admm_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE admm_update_s_mstruct(admm_env, qs_env)

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'admm_update_s_mstruct'

      INTEGER                                            :: handle
      LOGICAL                                            :: skip_load_balance_distributed
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      NULLIFY (ks_env, dft_control)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control)

      !The aux_fit task_list
      skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
      IF (ASSOCIATED(admm_env%task_list_aux_fit)) CALL deallocate_task_list(admm_env%task_list_aux_fit)
      CALL allocate_task_list(admm_env%task_list_aux_fit)
      CALL generate_qs_task_list(ks_env, admm_env%task_list_aux_fit, &
                                 reorder_rs_grid_ranks=.FALSE., soft_valid=.FALSE., basis_type="AUX_FIT", &
                                 skip_load_balance_distributed=skip_load_balance_distributed, &
                                 sab_orb_external=admm_env%sab_aux_fit)

      !The aux_fit densities
      CALL qs_rho_rebuild(admm_env%rho_aux_fit, qs_env=qs_env, admm=.TRUE.)
      CALL qs_rho_rebuild(admm_env%rho_aux_fit_buffer, qs_env=qs_env, admm=.TRUE.)

      CALL timestop(handle)

   END SUBROUTINE admm_update_s_mstruct

! **************************************************************************************************
!> \brief Update the admm_gapw_env internals to the current qs_env (i.e. atomic positions)
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE update_admm_gapw(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'update_admm_gapw'

      INTEGER                                            :: handle, ikind, nkind
      LOGICAL                                            :: paw_atom
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: aux_present, oce_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: aux_radius, oce_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(admm_gapw_type), POINTER                      :: admm_gapw_env
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(gto_basis_set_type), POINTER                  :: aux_fit_basis
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sap_oce
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(paw_proj_set_type), POINTER                   :: paw_proj
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: admm_kind_set, qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      NULLIFY (ks_env, qs_kind_set, admm_kind_set, aux_fit_basis, cell, distribution_1d)
      NULLIFY (distribution_2d, paw_proj, particle_set, molecule_set, admm_env, admm_gapw_env)
      NULLIFY (dft_control, atomic_kind_set, sap_oce)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, ks_env=ks_env, qs_kind_set=qs_kind_set, admm_env=admm_env, &
                      dft_control=dft_control)
      admm_gapw_env => admm_env%admm_gapw_env
      admm_kind_set => admm_gapw_env%admm_kind_set
      nkind = SIZE(qs_kind_set)

      !Update the task lisft for the AUX_FIT_SOFT basis
      IF (ASSOCIATED(admm_gapw_env%task_list)) CALL deallocate_task_list(admm_gapw_env%task_list)
      CALL allocate_task_list(admm_gapw_env%task_list)

      !note: we set soft_valid to .FALSE. want to use AUX_FIT_SOFT and not the normal ORB SOFT basis
      CALL generate_qs_task_list(ks_env, admm_gapw_env%task_list, reorder_rs_grid_ranks=.FALSE., &
                                 soft_valid=.FALSE., basis_type="AUX_FIT_SOFT", &
                                 skip_load_balance_distributed=dft_control%qs_control%skip_load_balance_distributed, &
                                 sab_orb_external=admm_env%sab_aux_fit)

      !Update the precomputed oce integrals
      !a sap_oce neighbor list is required => build it here
      ALLOCATE (aux_present(nkind), oce_present(nkind))
      aux_present = .FALSE.; oce_present = .FALSE.
      ALLOCATE (aux_radius(nkind), oce_radius(nkind))
      aux_radius = 0.0_dp; oce_radius = 0.0_dp

      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=aux_fit_basis, basis_type="AUX_FIT")
         IF (ASSOCIATED(aux_fit_basis)) THEN
            aux_present(ikind) = .TRUE.
            CALL get_gto_basis_set(aux_fit_basis, kind_radius=aux_radius(ikind))
         END IF

         !note: get oce info from admm_kind_set
         CALL get_qs_kind(admm_kind_set(ikind), paw_atom=paw_atom, paw_proj_set=paw_proj)
         IF (paw_atom) THEN
            oce_present(ikind) = .TRUE.
            CALL get_paw_proj_set(paw_proj, rcprj=oce_radius(ikind))
         END IF
      END DO

      ALLOCATE (pair_radius(nkind, nkind))
      pair_radius = 0.0_dp
      CALL pair_radius_setup(aux_present, oce_present, aux_radius, oce_radius, pair_radius)

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      distribution_2d=distribution_2d, local_particles=distribution_1d, &
                      particle_set=particle_set, molecule_set=molecule_set)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ALLOCATE (atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set)
      CALL build_neighbor_lists(sap_oce, particle_set, atom2d, cell, pair_radius, &
                                subcells=subcells, operator_type="ABBA", nlname="AUX_PAW-PRJ")
      CALL atom2d_cleanup(atom2d)

      !actually compute the oce matrices
      CALL create_oce_set(admm_gapw_env%oce)
      CALL allocate_oce_set(admm_gapw_env%oce, nkind)

      !always compute the derivative, cheap anyways
      CALL build_oce_matrices(admm_gapw_env%oce%intac, calculate_forces=.TRUE., nder=1, &
                              qs_kind_set=admm_kind_set, particle_set=particle_set, &
                              sap_oce=sap_oce, eps_fit=dft_control%qs_control%gapw_control%eps_fit)

      CALL release_neighbor_list_sets(sap_oce)

      CALL timestop(handle)

   END SUBROUTINE update_admm_gapw

! **************************************************************************************************
!> \brief Allocates the various ADMM KS matrices
!> \param admm_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE admm_alloc_ks_matrices(admm_env, qs_env)

      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'admm_alloc_ks_matrices'

      INTEGER                                            :: handle, ispin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (dft_control, matrix_s_aux_fit)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      matrix_s_aux_fit => admm_env%matrix_s_aux_fit
      IF (ASSOCIATED(admm_env%matrix_ks_aux_fit)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_ks_aux_fit)
      CALL dbcsr_allocate_matrix_set(admm_env%matrix_ks_aux_fit, dft_control%nspins)
      DO ispin = 1, dft_control%nspins
         ALLOCATE (admm_env%matrix_ks_aux_fit(ispin)%matrix)
         CALL dbcsr_create(admm_env%matrix_ks_aux_fit(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix, &
                           name="KOHN-SHAM_MATRIX for ADMM")
         CALL cp_dbcsr_alloc_block_from_nbl(admm_env%matrix_ks_aux_fit(ispin)%matrix, admm_env%sab_aux_fit)
         CALL dbcsr_set(admm_env%matrix_ks_aux_fit(ispin)%matrix, 0.0_dp)
      END DO

      IF (ASSOCIATED(admm_env%matrix_ks_aux_fit_dft)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_ks_aux_fit_dft)
      CALL dbcsr_allocate_matrix_set(admm_env%matrix_ks_aux_fit_dft, dft_control%nspins)
      DO ispin = 1, dft_control%nspins
         ALLOCATE (admm_env%matrix_ks_aux_fit_dft(ispin)%matrix)
         CALL dbcsr_create(admm_env%matrix_ks_aux_fit_dft(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix, &
                           name="KOHN-SHAM_MATRIX for ADMM")
         CALL cp_dbcsr_alloc_block_from_nbl(admm_env%matrix_ks_aux_fit_dft(ispin)%matrix, admm_env%sab_aux_fit)
         CALL dbcsr_set(admm_env%matrix_ks_aux_fit_dft(ispin)%matrix, 0.0_dp)
      END DO

      IF (ASSOCIATED(admm_env%matrix_ks_aux_fit_hfx)) CALL dbcsr_deallocate_matrix_set(admm_env%matrix_ks_aux_fit_hfx)
      CALL dbcsr_allocate_matrix_set(admm_env%matrix_ks_aux_fit_hfx, dft_control%nspins)
      DO ispin = 1, dft_control%nspins
         ALLOCATE (admm_env%matrix_ks_aux_fit_hfx(ispin)%matrix)
         CALL dbcsr_create(admm_env%matrix_ks_aux_fit_hfx(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix, &
                           name="KOHN-SHAM_MATRIX for ADMM")
         CALL cp_dbcsr_alloc_block_from_nbl(admm_env%matrix_ks_aux_fit_hfx(ispin)%matrix, admm_env%sab_aux_fit)
         CALL dbcsr_set(admm_env%matrix_ks_aux_fit_hfx(ispin)%matrix, 0.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE admm_alloc_ks_matrices

! **************************************************************************************************
!> \brief Add the hfx contributions to the Hamiltonian
!>
!> \param qs_env ...
!> \param matrix_ks ...
!> \param rho ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param v_rspace_new ...
!> \param v_tau_rspace ...
!> \par History
!>     refactoring 03-2011 [MI]
! **************************************************************************************************

   SUBROUTINE hfx_ks_matrix(qs_env, matrix_ks, rho, energy, calculate_forces, &
                            just_energy, v_rspace_new, v_tau_rspace)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy
      TYPE(pw_type), DIMENSION(:), POINTER               :: v_rspace_new, v_tau_rspace

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_ks_matrix'

      INTEGER                                            :: handle, img, irep, ispin, mspin, &
                                                            n_rep_hf, nimages, ns, nspins
      LOGICAL                                            :: distribute_fock_matrix, &
                                                            do_adiabatic_rescaling, &
                                                            hfx_treat_lsd_in_core, &
                                                            s_mstruct_changed, use_virial
      REAL(dp)                                           :: eh1, ehfx, ehfxrt, eold
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: hf_energy
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: matrix_ks_1d, matrix_ks_aux_fit, &
         matrix_ks_aux_fit_hfx, matrix_ks_aux_fit_im, matrix_ks_im, rho_ao_1d, rho_ao_resp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_h_im, matrix_ks_orb, &
                                                            rho_ao_orb
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo_array
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_orb
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(section_vals_type), POINTER                   :: adiabatic_rescaling_section, &
                                                            hfx_sections, input
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool, dft_control, hfx_sections, input, &
               para_env, poisson_env, pw_env, virial, matrix_ks_im, &
               matrix_ks_orb, rho_ao_orb, matrix_h, matrix_h_im, matrix_ks_aux_fit, &
               matrix_ks_aux_fit_im, matrix_ks_aux_fit_hfx)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      matrix_h_kp=matrix_h, &
                      matrix_h_im_kp=matrix_h_im, &
                      para_env=para_env, &
                      pw_env=pw_env, &
                      virial=virial, &
                      matrix_ks_im=matrix_ks_im, &
                      s_mstruct_changed=s_mstruct_changed, &
                      x_data=x_data)

      IF (dft_control%do_admm) THEN
         CALL get_admm_env(qs_env%admm_env, mos_aux_fit=mo_array, matrix_ks_aux_fit=matrix_ks_aux_fit, &
                           matrix_ks_aux_fit_im=matrix_ks_aux_fit_im, matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx)
      ELSE
         CALL get_qs_env(qs_env=qs_env, mos=mo_array)
      END IF

      nspins = dft_control%nspins
      nimages = dft_control%nimages

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      IF (use_virial .AND. calculate_forces) virial%pv_fock_4c = 0.0_dp

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                i_rep_section=1)
      adiabatic_rescaling_section => section_vals_get_subs_vals(input, "DFT%XC%ADIABATIC_RESCALING")
      CALL section_vals_get(adiabatic_rescaling_section, explicit=do_adiabatic_rescaling)

      ! *** Initialize the auxiliary ks matrix to zero if required
      IF (dft_control%do_admm) THEN
         DO ispin = 1, nspins
            CALL dbcsr_set(matrix_ks_aux_fit(ispin)%matrix, 0.0_dp)
         END DO
      END IF
      DO ispin = 1, nspins
         DO img = 1, nimages
            CALL dbcsr_set(matrix_ks(ispin, img)%matrix, 0.0_dp)
         END DO
      END DO

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)

      ALLOCATE (hf_energy(n_rep_hf))

      eold = 0.0_dp

      DO irep = 1, n_rep_hf
         ! Remember: Vhfx is added, energy is calclulated from total Vhfx,
         ! so energy of last iteration is correct

         IF (do_adiabatic_rescaling .AND. hfx_treat_lsd_in_core) &
            CPABORT("HFX_TREAT_LSD_IN_CORE not implemented for adiabatically rescaled hybrids")
         ! everything is calculated with adiabatic rescaling but the potential is not added in a first step
         distribute_fock_matrix = .NOT. do_adiabatic_rescaling

         mspin = 1
         IF (hfx_treat_lsd_in_core) mspin = nspins

         ! fetch the correct matrices for normal HFX or ADMM
         IF (dft_control%do_admm) THEN
            CALL get_admm_env(qs_env%admm_env, matrix_ks_aux_fit=matrix_ks_1d, rho_aux_fit=rho_orb)
            ns = SIZE(matrix_ks_1d)
            matrix_ks_orb(1:ns, 1:1) => matrix_ks_1d(1:ns)
         ELSE
            CALL get_qs_env(qs_env=qs_env, matrix_ks_kp=matrix_ks_orb, rho=rho_orb)
         END IF
         CALL qs_rho_get(rho_struct=rho_orb, rho_ao_kp=rho_ao_orb)
         ! Finally the real hfx calulation
         ehfx = 0.0_dp

         IF (x_data(irep, 1)%do_hfx_ri) THEN
            CALL hfx_ri_update_ks(qs_env, x_data(irep, 1)%ri_data, matrix_ks_orb, ehfx, &
                                  mo_array, rho_ao_orb, &
                                  s_mstruct_changed, nspins, &
                                  x_data(irep, 1)%general_parameter%fraction)
            IF (dft_control%do_admm) THEN
               !for ADMMS, we need the exchange matrix k(d) for both spins
               DO ispin = 1, nspins
                  CALL dbcsr_copy(matrix_ks_aux_fit_hfx(ispin)%matrix, matrix_ks_orb(ispin, 1)%matrix, &
                                  name="HF exch. part of matrix_ks_aux_fit for ADMMS")
               END DO
            END IF

         ELSE

            DO ispin = 1, mspin
               CALL integrate_four_center(qs_env, x_data, matrix_ks_orb, eh1, rho_ao_orb, hfx_sections, &
                                          para_env, s_mstruct_changed, irep, distribute_fock_matrix, &
                                          ispin=ispin)
               ehfx = ehfx + eh1
            END DO
         END IF

         IF (calculate_forces .AND. .NOT. do_adiabatic_rescaling) THEN
            !Scale auxiliary density matrix for ADMMP (see Merlot2014) with gsi(ispin) to scale force
            IF (dft_control%do_admm) THEN
               CALL scale_dm(qs_env, rho_ao_orb, scale_back=.FALSE.)
            END IF
            NULLIFY (rho_ao_resp)

            IF (x_data(irep, 1)%do_hfx_ri) THEN

               CALL hfx_ri_update_forces(qs_env, x_data(irep, 1)%ri_data, nspins, &
                                         x_data(irep, 1)%general_parameter%fraction, &
                                         rho_ao=rho_ao_orb, mos=mo_array, &
                                         rho_ao_resp=rho_ao_resp, &
                                         use_virial=use_virial)

            ELSE

               CALL derivatives_four_center(qs_env, rho_ao_orb, rho_ao_resp, hfx_sections, &
                                            para_env, irep, use_virial)

            END IF

            !Scale auxiliary density matrix for ADMMP back with 1/gsi(ispin)
            IF (dft_control%do_admm) THEN
               CALL scale_dm(qs_env, rho_ao_orb, scale_back=.TRUE.)
            END IF
         END IF

         !! If required, the calculation of the forces will be done later with adiabatic rescaling
         IF (do_adiabatic_rescaling) hf_energy(irep) = ehfx

         ! special case RTP/EMD we have a full complex density and HFX has a contribution from the imaginary part
         ehfxrt = 0.0_dp
         IF (qs_env%run_rtp) THEN

            CALL get_qs_env(qs_env=qs_env, rtp=rtp)
            DO ispin = 1, nspins
               CALL dbcsr_set(matrix_ks_im(ispin)%matrix, 0.0_dp)
            END DO
            IF (dft_control%do_admm) THEN
               ! matrix_ks_orb => matrix_ks_aux_fit_im
               ns = SIZE(matrix_ks_aux_fit_im)
               matrix_ks_orb(1:ns, 1:1) => matrix_ks_aux_fit_im(1:ns)
               DO ispin = 1, nspins
                  CALL dbcsr_set(matrix_ks_aux_fit_im(ispin)%matrix, 0.0_dp)
               END DO
            ELSE
               ! matrix_ks_orb => matrix_ks_im
               ns = SIZE(matrix_ks_im)
               matrix_ks_orb(1:ns, 1:1) => matrix_ks_im(1:ns)
            END IF

            CALL qs_rho_get(rho_orb, rho_ao_im=rho_ao_1d)
            ns = SIZE(rho_ao_1d)
            rho_ao_orb(1:ns, 1:1) => rho_ao_1d(1:ns)

            ehfxrt = 0.0_dp

            IF (x_data(irep, 1)%do_hfx_ri) THEN
               CALL hfx_ri_update_ks(qs_env, x_data(irep, 1)%ri_data, matrix_ks_orb, ehfx, &
                                     mo_array, rho_ao_orb, &
                                     .FALSE., nspins, &
                                     x_data(irep, 1)%general_parameter%fraction)
               IF (dft_control%do_admm) THEN
                  !for ADMMS, we need the exchange matrix k(d) for both spins
                  DO ispin = 1, nspins
                     CALL dbcsr_copy(matrix_ks_aux_fit_hfx(ispin)%matrix, matrix_ks_orb(ispin, 1)%matrix, &
                                     name="HF exch. part of matrix_ks_aux_fit for ADMMS")
                  END DO
               END IF

            ELSE
               DO ispin = 1, mspin
                  CALL integrate_four_center(qs_env, x_data, matrix_ks_orb, eh1, rho_ao_orb, hfx_sections, &
                                             para_env, .FALSE., irep, distribute_fock_matrix, &
                                             ispin=ispin)
                  ehfxrt = ehfxrt + eh1
               END DO
            END IF

            IF (calculate_forces .AND. .NOT. do_adiabatic_rescaling) THEN
               NULLIFY (rho_ao_resp)

               IF (x_data(irep, 1)%do_hfx_ri) THEN

                  CALL hfx_ri_update_forces(qs_env, x_data(irep, 1)%ri_data, nspins, &
                                            x_data(irep, 1)%general_parameter%fraction, &
                                            rho_ao=rho_ao_orb, mos=mo_array, &
                                            use_virial=use_virial)

               ELSE
                  CALL derivatives_four_center(qs_env, rho_ao_orb, rho_ao_resp, hfx_sections, &
                                               para_env, irep, use_virial)
               END IF
            END IF

            !! If required, the calculation of the forces will be done later with adiabatic rescaling
            IF (do_adiabatic_rescaling) hf_energy(irep) = ehfx + ehfxrt

            IF (dft_control%rtp_control%velocity_gauge) THEN
               CPASSERT(ASSOCIATED(matrix_h_im))
               DO ispin = 1, nspins
                  CALL dbcsr_add(matrix_ks_im(ispin)%matrix, matrix_h_im(1, 1)%matrix, &
                                 1.0_dp, 1.0_dp)
               END DO
            END IF

         END IF

         IF (.NOT. qs_env%run_rtp) THEN
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                            poisson_env=poisson_env)
            eh1 = ehfx - eold
            CALL pw_hfx(qs_env, eh1, hfx_sections, poisson_env, auxbas_pw_pool, irep)
            eold = ehfx
         END IF

      END DO

      ! *** Set the total HFX energy
      energy%ex = ehfx + ehfxrt

      ! *** Add Core-Hamiltonian-Matrix ***
      DO ispin = 1, nspins
         DO img = 1, nimages
            CALL dbcsr_add(matrix_ks(ispin, img)%matrix, matrix_h(1, img)%matrix, &
                           1.0_dp, 1.0_dp)
         END DO
      END DO
      IF (use_virial .AND. calculate_forces) THEN
         virial%pv_exx = virial%pv_exx - virial%pv_fock_4c
         virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
         virial%pv_calculate = .FALSE.
      END IF

      !! If we perform adiabatic rescaling we are now able to rescale the xc-potential
      IF (do_adiabatic_rescaling) THEN
         CALL rescale_xc_potential(qs_env, matrix_ks, rho, energy, v_rspace_new, v_tau_rspace, &
                                   hf_energy, just_energy, calculate_forces, use_virial)
      END IF ! do_adiabatic_rescaling

      !update the hfx aux_fit matrixIF (dft_control%do_admm) THEN
      IF (dft_control%do_admm) THEN
         DO ispin = 1, nspins
            CALL dbcsr_add(matrix_ks_aux_fit_hfx(ispin)%matrix, matrix_ks_aux_fit(ispin)%matrix, &
                           0.0_dp, 1.0_dp)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_ks_matrix

! **************************************************************************************************
!> \brief This routine modifies the xc section depending on the potential type
!>        used for the HF exchange and the resulting correction term. Currently
!>        three types of corrections are implemented:
!>
!>        coulomb:     Ex,hf = Ex,hf' + (PBEx-PBEx')
!>        shortrange:  Ex,hf = Ex,hf' + (XWPBEX-XWPBEX')
!>        truncated:   Ex,hf = Ex,hf' + ( (XWPBEX0-PBE_HOLE_TC_LR) -(XWPBEX0-PBE_HOLE_TC_LR)' )
!>
!>        with ' denoting the auxiliary basis set and
!>
!>          PBEx:           PBE exchange functional
!>          XWPBEX:         PBE exchange hole for short-range potential (erfc(omega*r)/r)
!>          XWPBEX0:        PBE exchange hole for standard coulomb potential
!>          PBE_HOLE_TC_LR: PBE exchange hole for longrange truncated coulomb potential
!>
!>        Above explanation is correct for the deafult case. If a specific functional is requested
!>        for the correction term (cfun), we get
!>        Ex,hf = Ex,hf' + (cfun-cfun')
!>        for all cases of operators.
!>
!> \param x_data ...
!> \param xc_section the original xc_section
!> \param admm_env the ADMM environment
!> \par History
!>      12.2009 created [Manuel Guidon]
!>      05.2021 simplify for case of no correction [JGH]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE create_admm_xc_section(x_data, xc_section, admm_env)
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(admm_type), POINTER                           :: admm_env

      LOGICAL, PARAMETER                                 :: debug_functional = .FALSE.
#if defined (__LIBXC)
      REAL(KIND=dp), PARAMETER :: x_factor_c = 0.930525736349100025_dp
#endif

      CHARACTER(LEN=20)                                  :: name_x_func
      INTEGER                                            :: hfx_potential_type, ifun, iounit, nfun
      LOGICAL                                            :: funct_found
      REAL(dp)                                           :: cutoff_radius, hfx_fraction, omega, &
                                                            scale_coulomb, scale_longrange, scale_x
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: xc_fun, xc_fun_section

      logger => cp_get_default_logger()
      NULLIFY (admm_env%xc_section_aux, admm_env%xc_section_primary)

      !! ** Duplicate existing xc-section
      CALL section_vals_duplicate(xc_section, admm_env%xc_section_aux)
      CALL section_vals_duplicate(xc_section, admm_env%xc_section_primary)
      !** Now modify the auxiliary basis
      !** First remove all functionals
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux, "XC_FUNCTIONAL")

      !* Overwrite possible shortcut
      CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                i_val=xc_funct_no_shortcut)

      !** Get number of Functionals in the list
      ifun = 0
      nfun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         nfun = nfun + 1
      END DO

      ifun = 0
      DO ifun = 1, nfun
         xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=1)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         CALL section_vals_remove_values(xc_fun)
      END DO

      IF (ASSOCIATED(x_data)) THEN
         hfx_potential_type = x_data(1, 1)%potential_parameter%potential_type
         hfx_fraction = x_data(1, 1)%general_parameter%fraction
      ELSE
         CPWARN("ADMM requested without a DFT%XC%HF section. It will be ignored for the SCF.")
         admm_env%aux_exch_func = do_admm_aux_exch_func_none
      END IF

      !in case of no admm exchange corr., no auxiliary exchange functional needed
      IF (admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                   i_val=xc_none)
         hfx_fraction = 0.0_dp
      ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_default) THEN
         ! default PBE Functional
         !! ** Add functionals evaluated with auxiliary basis
         SELECT CASE (hfx_potential_type)
         CASE (do_potential_coulomb)
            CALL section_vals_val_set(xc_fun_section, "PBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_X", &
                                      r_val=-hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_C", &
                                      r_val=0.0_dp)
         CASE (do_potential_short)
            omega = x_data(1, 1)%potential_parameter%omega
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=-hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=0.0_dp)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                      r_val=omega)
         CASE (do_potential_truncated)
            cutoff_radius = x_data(1, 1)%potential_parameter%cutoff_radius
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                      r_val=hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                      r_val=cutoff_radius)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=0.0_dp)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=-hfx_fraction)
         CASE (do_potential_long)
            omega = x_data(1, 1)%potential_parameter%omega
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=-hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                      r_val=omega)
         CASE (do_potential_mix_cl)
            omega = x_data(1, 1)%potential_parameter%omega
            scale_coulomb = x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = x_data(1, 1)%potential_parameter%scale_longrange
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=hfx_fraction*scale_longrange)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=-hfx_fraction*(scale_longrange + scale_coulomb))
            CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                      r_val=omega)
         CASE (do_potential_mix_cl_trunc)
            omega = x_data(1, 1)%potential_parameter%omega
            cutoff_radius = x_data(1, 1)%potential_parameter%cutoff_radius
            scale_coulomb = x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = x_data(1, 1)%potential_parameter%scale_longrange
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                      r_val=hfx_fraction*(scale_longrange + scale_coulomb))
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                      r_val=cutoff_radius)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                      r_val=hfx_fraction*scale_longrange)
            CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                      r_val=-hfx_fraction*(scale_longrange + scale_coulomb))
            CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                      r_val=omega)
         CASE DEFAULT
            CPABORT("Unknown potential operator!")
         END SELECT

         !** Now modify the functionals for the primary basis
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary, "XC_FUNCTIONAL")
         !* Overwrite possible shortcut
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                   i_val=xc_funct_no_shortcut)

         SELECT CASE (hfx_potential_type)
         CASE (do_potential_coulomb)
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "PBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "PBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_X", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_C", &
                                         r_val=0.0_dp)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "PBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "PBE%SCALE_X", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_short)
            omega = x_data(1, 1)%potential_parameter%omega
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=0.0_dp)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_long)
            omega = x_data(1, 1)%potential_parameter%omega
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=-hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            END IF
         CASE (do_potential_truncated)
            cutoff_radius = x_data(1, 1)%potential_parameter%cutoff_radius
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "PBE_HOLE_T_C_LR") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=-hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=0.0_dp)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_mix_cl_trunc)
            cutoff_radius = x_data(1, 1)%potential_parameter%cutoff_radius
            omega = x_data(1, 1)%potential_parameter%omega
            scale_coulomb = x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = x_data(1, 1)%potential_parameter%scale_longrange
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "PBE_HOLE_T_C_LR") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=-hfx_fraction*(scale_coulomb + scale_longrange))
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=hfx_fraction*(scale_coulomb + scale_longrange))
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=-hfx_fraction*scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*scale_longrange
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            END IF
         CASE (do_potential_mix_cl)
            omega = x_data(1, 1)%potential_parameter%omega
            scale_coulomb = x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = x_data(1, 1)%potential_parameter%scale_longrange
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "XWPBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "XWPBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=hfx_fraction*(scale_coulomb + scale_longrange))
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=-hfx_fraction*scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X0", &
                                         r_val=scale_x)

               CALL section_vals_val_get(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*scale_longrange
               CALL section_vals_val_set(xc_fun_section, "XWPBE%SCALE_X", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "XWPBE%OMEGA", &
                                         r_val=omega)
            END IF
         END SELECT
      ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_default_libxc) THEN
         ! default PBE Functional
         !! ** Add functionals evaluated with auxiliary basis
#if defined (__LIBXC)
         SELECT CASE (hfx_potential_type)
         CASE (do_potential_coulomb)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                      r_val=-hfx_fraction)
         CASE (do_potential_short)
            omega = x_data(1, 1)%potential_parameter%omega
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                      r_val=-hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                      r_val=omega)
         CASE (do_potential_truncated)
            cutoff_radius = x_data(1, 1)%potential_parameter%cutoff_radius
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                      r_val=hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                      r_val=cutoff_radius)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                      r_val=-hfx_fraction)
         CASE (do_potential_long)
            omega = x_data(1, 1)%potential_parameter%omega
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                      r_val=hfx_fraction)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                      r_val=omega)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                      r_val=-hfx_fraction)
         CASE (do_potential_mix_cl)
            omega = x_data(1, 1)%potential_parameter%omega
            scale_coulomb = x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = x_data(1, 1)%potential_parameter%scale_longrange
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                      r_val=hfx_fraction*scale_longrange)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                      r_val=omega)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                      r_val=-hfx_fraction*(scale_longrange + scale_coulomb))
         CASE (do_potential_mix_cl_trunc)
            omega = x_data(1, 1)%potential_parameter%omega
            cutoff_radius = x_data(1, 1)%potential_parameter%cutoff_radius
            scale_coulomb = x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = x_data(1, 1)%potential_parameter%scale_longrange
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                      r_val=hfx_fraction*(scale_longrange + scale_coulomb))
            CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                      r_val=cutoff_radius)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                      r_val=hfx_fraction*scale_longrange)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                      r_val=omega)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                      r_val=-hfx_fraction*(scale_longrange + scale_coulomb))
         CASE DEFAULT
            CPABORT("Unknown potential operator!")
         END SELECT

         !** Now modify the functionals for the primary basis
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary, "XC_FUNCTIONAL")
         !* Overwrite possible shortcut
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                   i_val=xc_funct_no_shortcut)

         SELECT CASE (hfx_potential_type)
         CASE (do_potential_coulomb)
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_PBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=hfx_fraction)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_short)
            omega = x_data(1, 1)%potential_parameter%omega
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_WPBEH") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                         r_val=omega)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_long)
            omega = x_data(1, 1)%potential_parameter%omega
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_WPBEH") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=-hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                         r_val=omega)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                         r_val=omega)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_PBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=hfx_fraction)
            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_truncated)
            cutoff_radius = x_data(1, 1)%potential_parameter%cutoff_radius
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "PBE_HOLE_T_C_LR") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=-hfx_fraction)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_PBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=hfx_fraction)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_mix_cl_trunc)
            cutoff_radius = x_data(1, 1)%potential_parameter%cutoff_radius
            omega = x_data(1, 1)%potential_parameter%omega
            scale_coulomb = x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = x_data(1, 1)%potential_parameter%scale_longrange
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "PBE_HOLE_T_C_LR") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=-hfx_fraction*(scale_coulomb + scale_longrange))
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%SCALE_X", &
                                         r_val=scale_x)
               CALL section_vals_val_set(xc_fun_section, "PBE_HOLE_T_C_LR%CUTOFF_RADIUS", &
                                         r_val=cutoff_radius)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_WPBEH") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=-hfx_fraction*scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                         r_val=omega)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*scale_longrange
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                         r_val=omega)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_PBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=hfx_fraction*(scale_coulomb + scale_longrange))
            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
            END IF
         CASE (do_potential_mix_cl)
            omega = x_data(1, 1)%potential_parameter%omega
            scale_coulomb = x_data(1, 1)%potential_parameter%scale_coulomb
            scale_longrange = x_data(1, 1)%potential_parameter%scale_longrange
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_WPBEH") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=-hfx_fraction*scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                         r_val=omega)

            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x - hfx_fraction*scale_longrange
               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%SCALE", &
                                         r_val=scale_x)

               CALL section_vals_val_set(xc_fun_section, "GGA_X_WPBEH%_OMEGA", &
                                         r_val=omega)
            END IF
            ifun = 0
            funct_found = .FALSE.
            DO
               ifun = ifun + 1
               xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
               IF (.NOT. ASSOCIATED(xc_fun)) EXIT
               IF (xc_fun%section%name == "GGA_X_PBE") THEN
                  funct_found = .TRUE.
               END IF
            END DO
            IF (.NOT. funct_found) THEN
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%_SECTION_PARAMETERS_", &
                                         l_val=.TRUE.)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=hfx_fraction*(scale_coulomb + scale_longrange))
            ELSE
               CALL section_vals_val_get(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
               scale_x = scale_x + hfx_fraction*(scale_coulomb + scale_longrange)
               CALL section_vals_val_set(xc_fun_section, "GGA_X_PBE%SCALE", &
                                         r_val=scale_x)
            END IF
         END SELECT
#else
         CALL cp_abort(__LOCATION__, "In order use a LibXC-based ADMM "// &
                       "exchange correction functionals, you have to compile and link against LibXC!")
#endif

         ! PBEX (always bare form), OPTX and Becke88 functional
      ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex .OR. &
               admm_env%aux_exch_func == do_admm_aux_exch_func_opt .OR. &
               admm_env%aux_exch_func == do_admm_aux_exch_func_bee) THEN
         IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex) THEN
            name_x_func = 'PBE'
         ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt) THEN
            name_x_func = 'OPTX'
         ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_bee) THEN
            name_x_func = 'BECKE88'
         END IF
         !primary basis
         CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_SECTION_PARAMETERS_", &
                                   l_val=.TRUE.)
         CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_X", &
                                   r_val=-hfx_fraction)

         IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex) THEN
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_C", r_val=0.0_dp)
         END IF

         IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt) THEN
            IF (admm_env%aux_exch_func_param) THEN
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%A1", &
                                         r_val=admm_env%aux_x_param(1))
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%A2", &
                                         r_val=admm_env%aux_x_param(2))
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%GAMMA", &
                                         r_val=admm_env%aux_x_param(3))
            END IF
         END IF

         !** Now modify the functionals for the primary basis
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary, "XC_FUNCTIONAL")
         !* Overwrite possible L")
         !* Overwrite possible shortcut
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                   i_val=xc_funct_no_shortcut)

         ifun = 0
         funct_found = .FALSE.
         DO
            ifun = ifun + 1
            xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
            IF (.NOT. ASSOCIATED(xc_fun)) EXIT
            IF (xc_fun%section%name == TRIM(name_x_func)) THEN
               funct_found = .TRUE.
            END IF
         END DO
         IF (.NOT. funct_found) THEN
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_X", &
                                      r_val=hfx_fraction)
            IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex) THEN
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_C", &
                                         r_val=0.0_dp)
            ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt) THEN
               IF (admm_env%aux_exch_func_param) THEN
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%A1", &
                                            r_val=admm_env%aux_x_param(1))
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%A2", &
                                            r_val=admm_env%aux_x_param(2))
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%GAMMA", &
                                            r_val=admm_env%aux_x_param(3))
               END IF
            END IF

         ELSE
            CALL section_vals_val_get(xc_fun_section, TRIM(name_x_func)//"%SCALE_X", &
                                      r_val=scale_x)
            scale_x = scale_x + hfx_fraction
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE_X", &
                                      r_val=scale_x)
            IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt) THEN
               CPASSERT(.NOT. admm_env%aux_exch_func_param)
            END IF
         END IF

      ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex_libxc .OR. &
               admm_env%aux_exch_func == do_admm_aux_exch_func_opt_libxc .OR. &
               admm_env%aux_exch_func == do_admm_aux_exch_func_sx_libxc .OR. &
               admm_env%aux_exch_func == do_admm_aux_exch_func_bee_libxc) THEN
#if defined(__LIBXC)
         IF (admm_env%aux_exch_func == do_admm_aux_exch_func_pbex_libxc) THEN
            name_x_func = 'GGA_X_PBE'
         ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt_libxc) THEN
            name_x_func = 'GGA_X_OPTX'
         ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_bee_libxc) THEN
            name_x_func = 'GGA_X_B88'
         ELSE IF (admm_env%aux_exch_func == do_admm_aux_exch_func_sx_libxc) THEN
            name_x_func = 'LDA_X'
         END IF
         !primary basis
         CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_SECTION_PARAMETERS_", &
                                   l_val=.TRUE.)
         CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE", &
                                   r_val=-hfx_fraction)

         IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt_libxc) THEN
            IF (admm_env%aux_exch_func_param) THEN
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_A", &
                                         r_val=admm_env%aux_x_param(1))
               ! LibXC rescales the second parameter of the OPTX functional (see documentation there)
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_B", &
                                         r_val=admm_env%aux_x_param(2)/x_factor_c)
               CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_GAMMA", &
                                         r_val=admm_env%aux_x_param(3))
            END IF
         END IF

         !** Now modify the functionals for the primary basis
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary, "XC_FUNCTIONAL")
         !* Overwrite possible L")
         !* Overwrite possible shortcut
         CALL section_vals_val_set(xc_fun_section, "_SECTION_PARAMETERS_", &
                                   i_val=xc_funct_no_shortcut)

         ifun = 0
         funct_found = .FALSE.
         DO
            ifun = ifun + 1
            xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
            IF (.NOT. ASSOCIATED(xc_fun)) EXIT
            IF (xc_fun%section%name == TRIM(name_x_func)) THEN
               funct_found = .TRUE.
            END IF
         END DO
         IF (.NOT. funct_found) THEN
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_SECTION_PARAMETERS_", &
                                      l_val=.TRUE.)
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE", &
                                      r_val=hfx_fraction)
            IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt_libxc) THEN
               IF (admm_env%aux_exch_func_param) THEN
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_A", &
                                            r_val=admm_env%aux_x_param(1))
                  ! LibXC rescales the second parameter of the OPTX functional (see documentation there)
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_B", &
                                            r_val=admm_env%aux_x_param(2)/x_factor_c)
                  CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%_GAMMA", &
                                            r_val=admm_env%aux_x_param(3))
               END IF
            END IF

         ELSE
            CALL section_vals_val_get(xc_fun_section, TRIM(name_x_func)//"%SCALE", &
                                      r_val=scale_x)
            scale_x = scale_x + hfx_fraction
            CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func)//"%SCALE", &
                                      r_val=scale_x)
            IF (admm_env%aux_exch_func == do_admm_aux_exch_func_opt_libxc) THEN
               CPASSERT(.NOT. admm_env%aux_exch_func_param)
            END IF
         END IF
#else
         CALL cp_abort(__LOCATION__, "In order use a LibXC-based ADMM "// &
                       "exchange correction functionals, you have to compile and link against LibXC!")
#endif

      ELSE
         CPABORT("Unknown exchange correction functional!")
      END IF

      IF (debug_functional) THEN
         iounit = cp_logger_get_default_io_unit(logger)
         IF (iounit > 0) THEN
            WRITE (iounit, "(A)") " ADMM Primary Basis Set Functional"
         END IF
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary, "XC_FUNCTIONAL")
         ifun = 0
         funct_found = .FALSE.
         DO
            ifun = ifun + 1
            xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
            IF (.NOT. ASSOCIATED(xc_fun)) EXIT

            scale_x = -1000.0_dp
            IF (xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
               CALL section_vals_val_get(xc_fun, "SCALE_X", r_val=scale_x)
            END IF
            IF (xc_fun%section%name == "XWPBE") THEN
               CALL section_vals_val_get(xc_fun, "SCALE_X0", r_val=hfx_fraction)
               IF (iounit > 0) THEN
                  WRITE (iounit, "(T5,A,T25,2F10.3)") TRIM(xc_fun%section%name), scale_x, hfx_fraction
               END IF
            ELSE
               IF (iounit > 0) THEN
                  WRITE (iounit, "(T5,A,T25,F10.3)") TRIM(xc_fun%section%name), scale_x
               END IF
            END IF
         END DO

         IF (iounit > 0) THEN
            WRITE (iounit, "(A)") " Auxiliary Basis Set Functional"
         END IF
         xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux, "XC_FUNCTIONAL")
         ifun = 0
         funct_found = .FALSE.
         DO
            ifun = ifun + 1
            xc_fun => section_vals_get_subs_vals2(xc_fun_section, i_section=ifun)
            IF (.NOT. ASSOCIATED(xc_fun)) EXIT
            scale_x = -1000.0_dp
            IF (xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
               CALL section_vals_val_get(xc_fun, "SCALE_X", r_val=scale_x)
            END IF
            IF (xc_fun%section%name == "XWPBE") THEN
               CALL section_vals_val_get(xc_fun, "SCALE_X0", r_val=hfx_fraction)
               IF (iounit > 0) THEN
                  WRITE (iounit, "(T5,A,T25,2F10.3)") TRIM(xc_fun%section%name), scale_x, hfx_fraction
               END IF
            ELSE
               IF (iounit > 0) THEN
                  WRITE (iounit, "(T5,A,T25,F10.3)") TRIM(xc_fun%section%name), scale_x
               END IF
            END IF
         END DO
      END IF

   END SUBROUTINE create_admm_xc_section

! **************************************************************************************************
!> \brief Add the hfx contributions to the Hamiltonian
!>
!> \param matrix_ks Kohn-Sham matrix (updated on exit)
!> \param rho_ao    electron density expressed in terms of atomic orbitals
!> \param qs_env    Quickstep environment
!> \param update_energy whether to update energy (default: yes)
!> \param recalc_integrals whether to recalculate integrals (default: value of HF%TREAT_LSD_IN_CORE)
!> \param external_hfx_sections ...
!> \param external_x_data ...
!> \note
!>     Simplified version of subroutine hfx_ks_matrix()
! **************************************************************************************************
   SUBROUTINE tddft_hfx_matrix(matrix_ks, rho_ao, qs_env, update_energy, recalc_integrals, &
                               external_hfx_sections, external_x_data)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         TARGET                                          :: matrix_ks, rho_ao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: update_energy, recalc_integrals
      TYPE(section_vals_type), OPTIONAL, POINTER         :: external_hfx_sections
      TYPE(hfx_type), DIMENSION(:, :), OPTIONAL, TARGET  :: external_x_data

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_hfx_matrix'

      INTEGER                                            :: handle, irep, ispin, mspin, n_rep_hf, &
                                                            nspins
      LOGICAL                                            :: distribute_fock_matrix, &
                                                            hfx_treat_lsd_in_core, &
                                                            my_update_energy, s_mstruct_changed
      REAL(KIND=dp)                                      :: eh1, ehfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, hfx_sections, input, para_env, matrix_ks_kp, rho_ao_kp)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      energy=energy, &
                      input=input, &
                      para_env=para_env, &
                      s_mstruct_changed=s_mstruct_changed, &
                      x_data=x_data)

      ! This should probably be the HF section from the TDDFPT XC section!
      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")

      IF (PRESENT(external_hfx_sections)) hfx_sections => external_hfx_sections
      IF (PRESENT(external_x_data)) x_data => external_x_data

      my_update_energy = .TRUE.
      IF (PRESENT(update_energy)) my_update_energy = update_energy

      IF (PRESENT(recalc_integrals)) s_mstruct_changed = recalc_integrals

      CPASSERT(dft_control%nimages == 1)
      nspins = dft_control%nspins

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                i_rep_section=1)

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
      distribute_fock_matrix = .TRUE.

      mspin = 1
      IF (hfx_treat_lsd_in_core) mspin = nspins

      matrix_ks_kp(1:nspins, 1:1) => matrix_ks(1:nspins)
      rho_ao_kp(1:nspins, 1:1) => rho_ao(1:nspins)

      DO irep = 1, n_rep_hf
         ! the real hfx calulation
         ehfx = 0.0_dp

         IF (x_data(irep, 1)%do_hfx_ri) THEN
            CALL hfx_ri_update_ks(qs_env, x_data(irep, 1)%ri_data, matrix_ks_kp, ehfx, &
                                  rho_ao=rho_ao_kp, geometry_did_change=s_mstruct_changed, &
                                  nspins=nspins, hf_fraction=x_data(irep, 1)%general_parameter%fraction)

         ELSE
            DO ispin = 1, mspin
               CALL integrate_four_center(qs_env, x_data, matrix_ks_kp, eh1, rho_ao_kp, hfx_sections, para_env, &
                                          s_mstruct_changed, irep, distribute_fock_matrix, ispin=ispin)
               ehfx = ehfx + eh1
            END DO
         END IF
      END DO
      IF (my_update_energy) energy%ex = ehfx

      CALL timestop(handle)
   END SUBROUTINE tddft_hfx_matrix

END MODULE hfx_admm_utils
