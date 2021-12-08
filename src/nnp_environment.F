!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods dealing with Neural Network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_environment

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Behler2007,&
                                              Behler2011,&
                                              Schran2020a,&
                                              Schran2020b,&
                                              cite_reference
   USE cell_methods,                    ONLY: read_cell,&
                                              write_cell
   USE cell_types,                      ONLY: cell_release,&
                                              cell_type,&
                                              get_cell
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_parser_methods,               ONLY: parser_read_line,&
                                              parser_search_string
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release,&
                                              parser_reset
   USE cp_subsys_methods,               ONLY: cp_subsys_create
   USE cp_subsys_types,                 ONLY: cp_subsys_release,&
                                              cp_subsys_set,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                              distribution_1d_type
   USE distribution_methods,            ONLY: distribute_molecules_1d
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                              write_molecule_kind_set
   USE molecule_types,                  ONLY: molecule_type
   USE nnp_acsf,                        ONLY: nnp_init_acsf_groups,&
                                              nnp_sort_acsf,&
                                              nnp_sort_ele,&
                                              nnp_write_acsf
   USE nnp_environment_types,           ONLY: &
        nnp_actfnct_cos, nnp_actfnct_exp, nnp_actfnct_gaus, nnp_actfnct_invsig, nnp_actfnct_lin, &
        nnp_actfnct_quad, nnp_actfnct_sig, nnp_actfnct_softplus, nnp_actfnct_tanh, nnp_env_set, &
        nnp_type
   USE nnp_model,                       ONLY: nnp_write_arc
   USE particle_methods,                ONLY: write_fist_particle_coordinates,&
                                              write_particle_distances,&
                                              write_structure_data
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_environment'

   PUBLIC :: nnp_init

CONTAINS

! **************************************************************************************************
!> \brief Read and initialize all the information for neural network potentials
!> \param nnp_env ...
!> \param root_section ...
!> \param para_env ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_init(nnp_env, root_section, para_env, force_env_section, subsys_section, &
                       use_motion_section)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: root_section
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(section_vals_type), INTENT(INOUT), POINTER    :: force_env_section, subsys_section
      LOGICAL, INTENT(IN)                                :: use_motion_section

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: explicit, use_ref_cell
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: cell_section, nnp_section

      CALL timeset(routineN, handle)
      CALL cite_reference(Behler2007)
      CALL cite_reference(Behler2011)
      CALL cite_reference(Schran2020a)
      CALL cite_reference(Schran2020b)

      CPASSERT(ASSOCIATED(nnp_env))

      NULLIFY (cell_section, nnp_section, cell, cell_ref, subsys)

      IF (.NOT. ASSOCIATED(subsys_section)) THEN
         subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      END IF
      cell_section => section_vals_get_subs_vals(subsys_section, "CELL")
      nnp_section => section_vals_get_subs_vals(force_env_section, "NNP")
      CALL section_vals_get(nnp_section, explicit=explicit)
      IF (.NOT. explicit) THEN
         CPWARN("NNP section not explicitly stated. Using default file names.")
      END IF

      CALL nnp_env_set(nnp_env=nnp_env, nnp_input=nnp_section, &
                       force_env_input=force_env_section)

      CALL read_cell(cell=cell, cell_ref=cell_ref, use_ref_cell=use_ref_cell, cell_section=cell_section, &
                     para_env=para_env)
      CALL get_cell(cell=cell, abc=abc)
      CALL write_cell(cell=cell, subsys_section=subsys_section)

      CALL cp_subsys_create(subsys, para_env, root_section, &
                            force_env_section=force_env_section, subsys_section=subsys_section, &
                            use_motion_section=use_motion_section)

      CALL nnp_init_subsys(nnp_env=nnp_env, subsys=subsys, cell=cell, &
                           cell_ref=cell_ref, use_ref_cell=use_ref_cell, &
                           subsys_section=subsys_section)

      CALL cell_release(cell)
      CALL cell_release(cell_ref)
      CALL cp_subsys_release(subsys)

      CALL timestop(handle)

   END SUBROUTINE nnp_init

! **************************************************************************************************
!> \brief Read and initialize all the information for neural network potentials
!> \param nnp_env ...
!> \param subsys ...
!> \param cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param subsys_section ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_init_subsys(nnp_env, subsys, cell, cell_ref, use_ref_cell, subsys_section)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
      TYPE(cp_subsys_type), INTENT(IN), POINTER          :: subsys
      TYPE(cell_type), INTENT(INOUT), POINTER            :: cell, cell_ref
      LOGICAL, INTENT(IN)                                :: use_ref_cell
      TYPE(section_vals_type), INTENT(IN), POINTER       :: subsys_section

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_init_subsys', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, natom
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, molecule_kind_set, particle_set, molecule_set, &
               local_molecules, local_particles)

      particle_set => subsys%particles%els
      atomic_kind_set => subsys%atomic_kinds%els
      molecule_kind_set => subsys%molecule_kinds%els
      molecule_set => subsys%molecules%els

      !Print the molecule kind set
      CALL write_molecule_kind_set(molecule_kind_set, subsys_section)

      !Print the atomic coordinates
      CALL write_fist_particle_coordinates(particle_set, subsys_section, &
                                           charges=Null())
      CALL write_particle_distances(particle_set, cell=cell, &
                                    subsys_section=subsys_section)
      CALL write_structure_data(particle_set, cell=cell, &
                                input_section=subsys_section)

      !Distribute molecules and atoms using the new data structures
      CALL distribute_molecules_1d(atomic_kind_set=atomic_kind_set, &
                                   particle_set=particle_set, &
                                   local_particles=local_particles, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   local_molecules=local_molecules, &
                                   force_env_section=nnp_env%force_env_input)

      natom = SIZE(particle_set)

      ALLOCATE (nnp_env%nnp_forces(3, natom))

      nnp_env%nnp_forces(:, :) = 0.0_dp

      nnp_env%nnp_potential_energy = 0.0_dp

      ! Set up arrays for calculation:
      nnp_env%num_atoms = natom
      ALLOCATE (nnp_env%ele_ind(natom))
      ALLOCATE (nnp_env%nuc_atoms(natom))
      ALLOCATE (nnp_env%coord(3, natom))
      ALLOCATE (nnp_env%atoms(natom))
      ALLOCATE (nnp_env%sort(natom))
      ALLOCATE (nnp_env%sort_inv(natom))

      CALL cp_subsys_set(subsys, cell=cell)

      CALL nnp_env_set(nnp_env=nnp_env, subsys=subsys, &
                       cell_ref=cell_ref, use_ref_cell=use_ref_cell, &
                       local_molecules=local_molecules, &
                       local_particles=local_particles)

      CALL distribution_1d_release(local_particles)
      CALL distribution_1d_release(local_molecules)

      CALL nnp_init_model(nnp_env=nnp_env)

      CALL timestop(handle)

   END SUBROUTINE nnp_init_subsys

! **************************************************************************************************
!> \brief Initialize the Neural Network Potential
!> \param nnp_env ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_init_model(nnp_env)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_init_model', routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: def_str_len = 256, &
                                                            default_path_length = 256

      CHARACTER(len=1), ALLOCATABLE, DIMENSION(:)        :: cactfnct
      CHARACTER(len=2)                                   :: ele
      CHARACTER(len=def_str_len)                         :: dummy, line
      CHARACTER(len=default_path_length)                 :: base_name, file_name
      INTEGER                                            :: handle, i, i_com, io, iweight, j, k, l, &
                                                            n_weight, nele, nuc_ele, symfnct_type, &
                                                            unit_nr
      LOGICAL                                            :: at_end, atom_e_found, explicit, first, &
                                                            found
      REAL(KIND=dp)                                      :: energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      REAL(KIND=dp), DIMENSION(7)                        :: test_array
      REAL(KIND=dp), DIMENSION(:), POINTER               :: work
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_parser_type), POINTER                      :: parser
      TYPE(section_vals_type), POINTER                   :: bias_section, model_section

      CALL timeset(routineN, handle)

      NULLIFY (parser, logger)

      logger => cp_get_default_logger()

      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) "NNP| Neural Network Potential Force Environment"
      END IF

      model_section => section_vals_get_subs_vals(nnp_env%nnp_input, "MODEL")
      CALL section_vals_get(model_section, n_repetition=nnp_env%n_committee)
      ALLOCATE (nnp_env%atomic_energy(nnp_env%num_atoms, nnp_env%n_committee))
      ALLOCATE (nnp_env%committee_energy(nnp_env%n_committee))
      ALLOCATE (nnp_env%myforce(3, nnp_env%num_atoms, nnp_env%n_committee))
      ALLOCATE (nnp_env%committee_forces(3, nnp_env%num_atoms, nnp_env%n_committee))
      ALLOCATE (nnp_env%committee_stress(3, 3, nnp_env%n_committee))

      CALL section_vals_val_get(nnp_env%nnp_input, "NNP_INPUT_FILE_NAME", c_val=file_name)
      CALL parser_create(parser, file_name, para_env=logger%para_env)

      ! read number of elements and cut_type and check for scale and center
      nnp_env%scale_acsf = .FALSE.
      nnp_env%scale_sigma_acsf = .FALSE.
      ! Defaults for scale min and max:
      nnp_env%scmin = 0.0_dp
      nnp_env%scmax = 1.0_dp
      nnp_env%center_acsf = .FALSE.
      nnp_env%normnodes = .FALSE.
      nnp_env%n_hlayer = 0

      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
         WRITE (unit_nr, *) "NNP| Reading NNP input from file: ", TRIM(file_name)
      END IF

      CALL parser_search_string(parser, "number_of_elements", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, nnp_env%n_ele
      ELSE
         CPABORT("NNP| number of elements missing in NNP_INPUT_FILE")
      END IF

      CALL parser_search_string(parser, "scale_symmetry_functions_sigma", .TRUE., found, &
                                search_from_begin_of_file=.TRUE.)
      nnp_env%scale_sigma_acsf = found

      CALL parser_search_string(parser, "scale_symmetry_functions", .TRUE., found, &
                                search_from_begin_of_file=.TRUE.)
      nnp_env%scale_acsf = found

      ! Test if there are two keywords of this:
      CALL parser_search_string(parser, "scale_symmetry_functions", .TRUE., found)
      IF (found .AND. nnp_env%scale_sigma_acsf) THEN
         CPWARN('Two scaling keywords in the input, we will ignore sigma scaling in this case')
         nnp_env%scale_sigma_acsf = .FALSE.
      ELSE IF (.NOT. found .AND. nnp_env%scale_sigma_acsf) THEN
         nnp_env%scale_acsf = .FALSE.
      END IF

      CALL parser_search_string(parser, "scale_min_short_atomic", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) READ (line, *) dummy, nnp_env%scmin

      CALL parser_search_string(parser, "scale_max_short_atomic", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) READ (line, *) dummy, nnp_env%scmax

      CALL parser_search_string(parser, "center_symmetry_functions", .TRUE., found, &
                                search_from_begin_of_file=.TRUE.)
      nnp_env%center_acsf = found
      ! n2p2 overwrites sigma scaling, if centering is requested:
      IF (nnp_env%scale_sigma_acsf .AND. nnp_env%center_acsf) THEN
         nnp_env%scale_sigma_acsf = .FALSE.
      END IF

      CALL parser_search_string(parser, "normalize_nodes", .TRUE., found, &
                                search_from_begin_of_file=.TRUE.)
      nnp_env%normnodes = found

      CALL parser_search_string(parser, "cutoff_type", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, nnp_env%cut_type
      ELSE
         CPABORT("NNP| no cutoff type specified in NNP_INPUT_FILE")
      END IF

      CALL parser_search_string(parser, "global_hidden_layers_short", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, nnp_env%n_hlayer
      ELSE
         CPABORT("NNP| number of hidden layers missing in NNP_INPUT_FILE")
      END IF
      nnp_env%n_layer = nnp_env%n_hlayer + 2

      nele = nnp_env%n_ele
      ALLOCATE (nnp_env%rad(nele))
      ALLOCATE (nnp_env%ang(nele))
      ALLOCATE (nnp_env%n_rad(nele))
      ALLOCATE (nnp_env%n_ang(nele))
      ALLOCATE (nnp_env%actfnct(nnp_env%n_hlayer + 1))
      ALLOCATE (cactfnct(nnp_env%n_hlayer + 1))
      ALLOCATE (nnp_env%ele(nele))
      ALLOCATE (nnp_env%nuc_ele(nele))
      ALLOCATE (nnp_env%arc(nele))
      DO i = 1, nele
         ALLOCATE (nnp_env%arc(i)%layer(nnp_env%n_layer))
         ALLOCATE (nnp_env%arc(i)%n_nodes(nnp_env%n_layer))
      END DO
      ALLOCATE (nnp_env%n_hnodes(nnp_env%n_hlayer))
      ALLOCATE (nnp_env%atom_energies(nele))
      nnp_env%atom_energies = 0.0_dp

      ! read elements, broadcast and sort
      CALL parser_reset(parser)
      DO
         CALL parser_search_string(parser, "elements", .TRUE., found, line)
         IF (found) THEN
            READ (line, *) dummy
            IF (TRIM(ADJUSTL(dummy)) == "elements") THEN
               READ (line, *) dummy, nnp_env%ele(:)
               CALL nnp_sort_ele(nnp_env%ele, nnp_env%nuc_ele)
               EXIT
            END IF
         ELSE
            CPABORT("NNP| elements not specified in NNP_INPUT_FILE")
         END IF
      END DO

      CALL parser_search_string(parser, "remove_atom_energies", .TRUE., atom_e_found, &
                                search_from_begin_of_file=.TRUE.)

      IF (atom_e_found) THEN
         CALL parser_reset(parser)
         i = 0
         DO
            CALL parser_search_string(parser, "atom_energy", .TRUE., found, line)
            IF (found) THEN
               READ (line, *) dummy, ele, energy
               DO j = 1, nele
                  IF (nnp_env%ele(j) == TRIM(ele)) THEN
                     i = i + 1
                     nnp_env%atom_energies(j) = energy
                  END IF
               END DO
               IF (i == nele) EXIT
            ELSE
               CPABORT("NNP| atom energies are not specified")
            END IF
         END DO
      END IF

      CALL parser_search_string(parser, "global_nodes_short", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, nnp_env%n_hnodes(:)
      ELSE
         CPABORT("NNP| global_nodes_short not specified in NNP_INPUT_FILE")
      END IF

      CALL parser_search_string(parser, "global_activation_short", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, cactfnct(:)
      ELSE
         CPABORT("NNP| global_activation_short not specified in NNP_INPUT_FILE")
      END IF

      DO i = 1, nnp_env%n_hlayer + 1
         SELECT CASE (cactfnct(i))
         CASE ("t")
            nnp_env%actfnct(i) = nnp_actfnct_tanh
         CASE ("g")
            nnp_env%actfnct(i) = nnp_actfnct_gaus
         CASE ("l")
            nnp_env%actfnct(i) = nnp_actfnct_lin
         CASE ("c")
            nnp_env%actfnct(i) = nnp_actfnct_cos
         CASE ("s")
            nnp_env%actfnct(i) = nnp_actfnct_sig
         CASE ("S")
            nnp_env%actfnct(i) = nnp_actfnct_invsig
         CASE ("e")
            nnp_env%actfnct(i) = nnp_actfnct_exp
         CASE ("p")
            nnp_env%actfnct(i) = nnp_actfnct_softplus
         CASE ("h")
            nnp_env%actfnct(i) = nnp_actfnct_quad
         CASE DEFAULT
            CPABORT("NNP| Activation function unkown")
         END SELECT
      END DO

      ! determine n_rad and n_ang
      DO i = 1, nele
         nnp_env%n_rad(i) = 0
         nnp_env%n_ang(i) = 0
      END DO

      ! count symfunctions
      CALL parser_reset(parser)
      first = .TRUE.
      DO
         CALL parser_search_string(parser, "symfunction_short", .TRUE., found, line)
         IF (found) THEN
            READ (line, *) dummy, ele, symfnct_type
            DO i = 1, nele
               IF (TRIM(ele) .EQ. nnp_env%ele(i)) THEN
                  IF (symfnct_type .EQ. 2) THEN
                     nnp_env%n_rad(i) = nnp_env%n_rad(i) + 1
                  ELSE IF (symfnct_type .EQ. 3) THEN
                     nnp_env%n_ang(i) = nnp_env%n_ang(i) + 1
                  ELSE
                     CPABORT("NNP| Symmetry function type not supported")
                  END IF
               END IF
            END DO
            first = .FALSE.
         ELSE
            IF (first) CPABORT("NNP| no symfunction_short specified in NNP_INPUT_FILE")
            ! no additional symfnct found
            EXIT
         END IF
      END DO

      DO i = 1, nele
         ALLOCATE (nnp_env%rad(i)%y(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%funccut(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%eta(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%rs(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%loc_min(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%loc_max(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%loc_av(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%sigma(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%ele(nnp_env%n_rad(i)))
         ALLOCATE (nnp_env%rad(i)%nuc_ele(nnp_env%n_rad(i)))
         nnp_env%rad(i)%funccut = 0.0_dp
         nnp_env%rad(i)%eta = 0.0_dp
         nnp_env%rad(i)%rs = 0.0_dp
         nnp_env%rad(i)%ele = 'X'
         nnp_env%rad(i)%nuc_ele = 0

         ALLOCATE (nnp_env%ang(i)%y(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%funccut(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%eta(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%zeta(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%prefzeta(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%lam(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%loc_min(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%loc_max(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%loc_av(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%sigma(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%ele1(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%ele2(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%nuc_ele1(nnp_env%n_ang(i)))
         ALLOCATE (nnp_env%ang(i)%nuc_ele2(nnp_env%n_ang(i)))
         nnp_env%ang(i)%funccut = 0.0_dp
         nnp_env%ang(i)%eta = 0.0_dp
         nnp_env%ang(i)%zeta = 0.0_dp
         nnp_env%ang(i)%prefzeta = 1.0_dp
         nnp_env%ang(i)%lam = 0.0_dp
         nnp_env%ang(i)%ele1 = 'X'
         nnp_env%ang(i)%ele2 = 'X'
         nnp_env%ang(i)%nuc_ele1 = 0
         nnp_env%ang(i)%nuc_ele2 = 0

         ! set number of nodes
         nnp_env%arc(i)%n_nodes(1) = nnp_env%n_rad(i) + nnp_env%n_ang(i)
         nnp_env%arc(i)%n_nodes(2:nnp_env%n_layer - 1) = nnp_env%n_hnodes
         nnp_env%arc(i)%n_nodes(nnp_env%n_layer) = 1
         DO j = 1, nnp_env%n_layer
            ALLOCATE (nnp_env%arc(i)%layer(j)%node(nnp_env%arc(i)%n_nodes(j)))
            ALLOCATE (nnp_env%arc(i)%layer(j)%node_grad(nnp_env%arc(i)%n_nodes(j)))
            ALLOCATE (nnp_env%arc(i)%layer(j)%tmp_der(nnp_env%arc(i)%n_nodes(1), nnp_env%arc(i)%n_nodes(j)))
         END DO
      END DO

      ! read, bcast and sort symfnct parameters
      DO i = 1, nele
         nnp_env%n_rad(i) = 0
         nnp_env%n_ang(i) = 0
      END DO
      CALL parser_reset(parser)
      first = .TRUE.
      nnp_env%max_cut = 0.0_dp
      DO
         CALL parser_search_string(parser, "symfunction_short", .TRUE., found, line)
         IF (found) THEN
            READ (line, *) dummy, ele, symfnct_type
            DO i = 1, nele
               IF (TRIM(ele) .EQ. nnp_env%ele(i)) THEN
                  IF (symfnct_type .EQ. 2) THEN
                     nnp_env%n_rad(i) = nnp_env%n_rad(i) + 1
                     READ (line, *) dummy, ele, symfnct_type, &
                        nnp_env%rad(i)%ele(nnp_env%n_rad(i)), &
                        nnp_env%rad(i)%eta(nnp_env%n_rad(i)), &
                        nnp_env%rad(i)%rs(nnp_env%n_rad(i)), &
                        nnp_env%rad(i)%funccut(nnp_env%n_rad(i))
                     IF (nnp_env%max_cut < nnp_env%rad(i)%funccut(nnp_env%n_rad(i))) THEN
                        nnp_env%max_cut = nnp_env%rad(i)%funccut(nnp_env%n_rad(i))
                     END IF
                  ELSE IF (symfnct_type .EQ. 3) THEN
                     nnp_env%n_ang(i) = nnp_env%n_ang(i) + 1
                     READ (line, *) dummy, ele, symfnct_type, &
                        nnp_env%ang(i)%ele1(nnp_env%n_ang(i)), &
                        nnp_env%ang(i)%ele2(nnp_env%n_ang(i)), &
                        nnp_env%ang(i)%eta(nnp_env%n_ang(i)), &
                        nnp_env%ang(i)%lam(nnp_env%n_ang(i)), &
                        nnp_env%ang(i)%zeta(nnp_env%n_ang(i)), &
                        nnp_env%ang(i)%funccut(nnp_env%n_ang(i))
                     nnp_env%ang(i)%prefzeta(nnp_env%n_ang(i)) = &
                        2.0_dp**(1.0_dp - nnp_env%ang(i)%zeta(nnp_env%n_ang(i)))
                     IF (nnp_env%max_cut < nnp_env%ang(i)%funccut(nnp_env%n_ang(i))) THEN
                        nnp_env%max_cut = nnp_env%ang(i)%funccut(nnp_env%n_ang(i))
                     END IF
                  ELSE
                     CPABORT("NNP| Symmetry function type not supported")
                  END IF
               END IF
            END DO
            first = .FALSE.
         ELSE
            IF (first) CPABORT("NNP| no symfunction_short specified in NNP_INPUT_FILE")
            ! no additional symfnct found
            EXIT
         END IF
      END DO

      DO i = 1, nele
         DO j = 1, nnp_env%n_rad(i)
            CALL get_ptable_info(nnp_env%rad(i)%ele(j), number=nnp_env%rad(i)%nuc_ele(j))
         END DO
         DO j = 1, nnp_env%n_ang(i)
            CALL get_ptable_info(nnp_env%ang(i)%ele1(j), number=nnp_env%ang(i)%nuc_ele1(j))
            CALL get_ptable_info(nnp_env%ang(i)%ele2(j), number=nnp_env%ang(i)%nuc_ele2(j))
            ! sort ele1 and ele2
            IF (nnp_env%ang(i)%nuc_ele1(j) .GT. nnp_env%ang(i)%nuc_ele2(j)) THEN
               ele = nnp_env%ang(i)%ele1(j)
               nnp_env%ang(i)%ele1(j) = nnp_env%ang(i)%ele2(j)
               nnp_env%ang(i)%ele2(j) = ele
               nuc_ele = nnp_env%ang(i)%nuc_ele1(j)
               nnp_env%ang(i)%nuc_ele1(j) = nnp_env%ang(i)%nuc_ele2(j)
               nnp_env%ang(i)%nuc_ele2(j) = nuc_ele
            END IF
         END DO
      END DO
      ! Done with input.nn file
      CALL parser_release(parser)

      ! sort symmetry functions and output information
      CALL nnp_sort_acsf(nnp_env)
      CALL nnp_write_acsf(nnp_env, logger%para_env)
      CALL nnp_write_arc(nnp_env, logger%para_env)

      ! read scaling information from file
      IF (nnp_env%scale_acsf .OR. nnp_env%center_acsf .OR. nnp_env%scale_sigma_acsf) THEN
         IF (logger%para_env%ionode) THEN
            WRITE (unit_nr, *) "NNP| Reading scaling information from file: ", TRIM(file_name)
         END IF
         CALL section_vals_val_get(nnp_env%nnp_input, "SCALE_FILE_NAME", &
                                   c_val=file_name)
         CALL parser_create(parser, file_name, para_env=logger%para_env)

         ! Get number of elements in scaling file
         CALL parser_read_line(parser, 1)
         k = 0
         DO WHILE (k < 7)
            READ (parser%input_line, *, IOSTAT=io) test_array(1:k)
            IF (io == -1) EXIT
            k = k + 1
         END DO
         k = k - 1

         IF (k == 5 .AND. nnp_env%scale_sigma_acsf) THEN
            CPABORT("Sigma scaling requested, but scaling.data does not contain sigma.")
         END IF

         CALL parser_reset(parser)
         DO i = 1, nnp_env%n_ele
            DO j = 1, nnp_env%n_rad(i)
               CALL parser_read_line(parser, 1)
               IF (nnp_env%scale_sigma_acsf) THEN
                  READ (parser%input_line, *) dummy, dummy, &
                     nnp_env%rad(i)%loc_min(j), &
                     nnp_env%rad(i)%loc_max(j), &
                     nnp_env%rad(i)%loc_av(j), &
                     nnp_env%rad(i)%sigma(j)
               ELSE
                  READ (parser%input_line, *) dummy, dummy, &
                     nnp_env%rad(i)%loc_min(j), &
                     nnp_env%rad(i)%loc_max(j), &
                     nnp_env%rad(i)%loc_av(j)
               END IF
            END DO
            DO j = 1, nnp_env%n_ang(i)
               CALL parser_read_line(parser, 1)
               IF (nnp_env%scale_sigma_acsf) THEN
                  READ (parser%input_line, *) dummy, dummy, &
                     nnp_env%ang(i)%loc_min(j), &
                     nnp_env%ang(i)%loc_max(j), &
                     nnp_env%ang(i)%loc_av(j), &
                     nnp_env%ang(i)%sigma(j)
               ELSE
                  READ (parser%input_line, *) dummy, dummy, &
                     nnp_env%ang(i)%loc_min(j), &
                     nnp_env%ang(i)%loc_max(j), &
                     nnp_env%ang(i)%loc_av(j)
               END IF
            END DO
         END DO
         CALL parser_release(parser)
      END IF

      CALL nnp_init_acsf_groups(nnp_env)

      ! read weights from file
      DO i = 1, nnp_env%n_ele
         DO j = 2, nnp_env%n_layer
            ALLOCATE (nnp_env%arc(i)%layer(j)%weights(nnp_env%arc(i)%n_nodes(j - 1), &
                                                      nnp_env%arc(i)%n_nodes(j), nnp_env%n_committee))
            ALLOCATE (nnp_env%arc(i)%layer(j)%bweights(nnp_env%arc(i)%n_nodes(j), nnp_env%n_committee))
         END DO
      END DO
      DO i_com = 1, nnp_env%n_committee
         CALL section_vals_val_get(model_section, "WEIGHTS", c_val=base_name, i_rep_section=i_com)
         IF (logger%para_env%ionode) THEN
            WRITE (unit_nr, *) "NNP| Initializing weights for model: ", i_com
         END IF
         DO i = 1, nnp_env%n_ele
            WRITE (file_name, '(A,I0.3,A)') TRIM(base_name)//".", nnp_env%nuc_ele(i), ".data"
            IF (logger%para_env%ionode) THEN
               WRITE (unit_nr, *) "NNP| Reading weights from file: ", TRIM(file_name)
            END IF
            CALL parser_create(parser, file_name, para_env=logger%para_env)
            n_weight = 0
            DO WHILE (.TRUE.)
               CALL parser_read_line(parser, 1, at_end)
               IF (at_end) EXIT
               n_weight = n_weight + 1
            END DO

            ALLOCATE (weights(n_weight))

            CALL parser_reset(parser)
            DO j = 1, n_weight
               CALL parser_read_line(parser, 1)
               READ (parser%input_line, *) weights(j)
            END DO
            CALL parser_release(parser)

            ! sort weights into corresponding arrays
            iweight = 0
            DO j = 2, nnp_env%n_layer
               DO k = 1, nnp_env%arc(i)%n_nodes(j - 1)
                  DO l = 1, nnp_env%arc(i)%n_nodes(j)
                     iweight = iweight + 1
                     nnp_env%arc(i)%layer(j)%weights(k, l, i_com) = weights(iweight)
                  END DO
               END DO

               DO k = 1, nnp_env%arc(i)%n_nodes(j)
                  iweight = iweight + 1
                  nnp_env%arc(i)%layer(j)%bweights(k, i_com) = weights(iweight)
               END DO
            END DO

            DEALLOCATE (weights)
         END DO
      END DO

      !Initialize extrapolation counter
      nnp_env%expol = 0

      ! Bias the standard deviation of committee disagreement
      bias_section => section_vals_get_subs_vals(nnp_env%nnp_input, "BIAS")
      CALL section_vals_get(bias_section, explicit=explicit)
      nnp_env%bias = .FALSE.
      IF (explicit) THEN
         IF (nnp_env%n_committee > 1) THEN
            IF (logger%para_env%ionode) THEN
               WRITE (unit_nr, *) "NNP| Biasing of committee disagreement enabled"
            END IF
            nnp_env%bias = .TRUE.
            ALLOCATE (nnp_env%bias_forces(3, nnp_env%num_atoms))
            ALLOCATE (nnp_env%bias_e_avrg(nnp_env%n_committee))
            CALL section_vals_val_get(bias_section, "SIGMA_0", r_val=nnp_env%bias_sigma0)
            CALL section_vals_val_get(bias_section, "K_B", r_val=nnp_env%bias_kb)
            nnp_env%bias_e_avrg(:) = 0.0_dp
            CALL section_vals_val_get(bias_section, "ALIGN_NNP_ENERGIES", explicit=explicit)
            nnp_env%bias_align = explicit
            IF (explicit) THEN
               NULLIFY (work)
               CALL section_vals_val_get(bias_section, "ALIGN_NNP_ENERGIES", r_vals=work)
               IF (SIZE(work) .NE. nnp_env%n_committee) THEN
                  CPABORT("ALIGN_NNP_ENERGIES size mismatch wrt committee size.")
               END IF
               nnp_env%bias_e_avrg(:) = work
               IF (logger%para_env%ionode) THEN
                  WRITE (unit_nr, *) "NNP| Biasing is aligned by shifting the energy prediction of the C-NNP members"
               END IF
            END IF
         ELSE
            CPWARN("NNP committee size is 1, BIAS section is ignored.")
         END IF
      END IF

      IF (logger%para_env%ionode) THEN
         WRITE (unit_nr, *) "NNP| NNP force environment initialized"
      END IF

      CALL timestop(handle)

   END SUBROUTINE nnp_init_model

END MODULE nnp_environment
