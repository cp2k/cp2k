!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods dealing with Neural Network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_environment

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Behler2007,&
                                              Behler2011,&
                                              Schran2020a,&
                                              Schran2020b,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_parser_methods,               ONLY: parser_read_line,&
                                              parser_search_string
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release,&
                                              parser_reset
   USE cp_subsys_methods,               ONLY: cp_subsys_create
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                              distribution_1d_type
   USE distribution_methods,            ONLY: distribute_molecules_1d
   USE input_constants,                 ONLY: nnp_gpu_prec_double,&
                                              nnp_gpu_prec_mixed,&
                                              nnp_use_gpu_auto,&
                                              nnp_use_gpu_off,&
                                              nnp_use_gpu_on
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
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
   USE nnp_gpu,                         ONLY: nnp_gpu_available,&
                                              nnp_gpu_init_device,&
                                              nnp_gpu_model_capacity_ok,&
                                              nnp_gpu_upload_descriptor_statics,&
                                              nnp_gpu_upload_network
   USE nnp_model,                       ONLY: nnp_write_arc
   USE particle_methods,                ONLY: write_fist_particle_coordinates,&
                                              write_particle_distances,&
                                              write_structure_data
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_environment'

   PUBLIC :: nnp_init
   PUBLIC :: nnp_init_model

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
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(section_vals_type), INTENT(INOUT), POINTER    :: force_env_section, subsys_section
      LOGICAL, INTENT(IN)                                :: use_motion_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_init'

      INTEGER                                            :: handle
      LOGICAL                                            :: explicit
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: nnp_section

      CALL timeset(routineN, handle)
      CALL cite_reference(Behler2007)
      CALL cite_reference(Behler2011)
      CALL cite_reference(Schran2020a)
      CALL cite_reference(Schran2020b)

      CPASSERT(ASSOCIATED(nnp_env))

      NULLIFY (nnp_section, subsys)

      IF (.NOT. ASSOCIATED(subsys_section)) THEN
         subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      END IF
      nnp_section => section_vals_get_subs_vals(force_env_section, "NNP")
      CALL section_vals_get(nnp_section, explicit=explicit)
      IF (.NOT. explicit) THEN
         CPWARN("NNP section not explicitly stated. Using default file names.")
      END IF

      CALL nnp_env_set(nnp_env=nnp_env, nnp_input=nnp_section, &
                       force_env_input=force_env_section)

      CALL cp_subsys_create(subsys, para_env, root_section, &
                            force_env_section=force_env_section, subsys_section=subsys_section, &
                            use_motion_section=use_motion_section)

      CALL nnp_init_subsys(nnp_env=nnp_env, subsys=subsys, &
                           subsys_section=subsys_section)

      CALL timestop(handle)

   END SUBROUTINE nnp_init

! **************************************************************************************************
!> \brief Read and initialize all the information for neural network potentials
!> \param nnp_env ...
!> \param subsys ...
!> \param subsys_section ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_init_subsys(nnp_env, subsys, subsys_section)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
      TYPE(cp_subsys_type), INTENT(IN), POINTER          :: subsys
      TYPE(section_vals_type), INTENT(IN), POINTER       :: subsys_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_init_subsys'

      INTEGER                                            :: handle, natom
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: my_cell, my_cell_ref
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, molecule_kind_set, my_cell, my_cell_ref, &
               particle_set, molecule_set, local_molecules, local_particles)

      particle_set => subsys%particles%els
      atomic_kind_set => subsys%atomic_kinds%els
      molecule_kind_set => subsys%molecule_kinds%els
      molecule_set => subsys%molecules%els
      my_cell => subsys%cell
      my_cell_ref => subsys%cell_ref

      !Print the molecule kind set
      CALL write_molecule_kind_set(molecule_kind_set, subsys_section)

      !Set up cell
      CALL nnp_env_set(nnp_env=nnp_env, subsys=subsys, &
                       cell=my_cell, cell_ref=my_cell_ref, &
                       use_ref_cell=subsys%use_ref_cell)

      !Print the atomic coordinates
      CALL write_fist_particle_coordinates(particle_set, subsys_section)
      CALL write_particle_distances(particle_set, cell=my_cell, &
                                    subsys_section=subsys_section)
      CALL write_structure_data(particle_set, cell=my_cell, &
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

      CALL nnp_env_set(nnp_env=nnp_env, &
                       local_molecules=local_molecules, &
                       local_particles=local_particles)

      CALL distribution_1d_release(local_particles)
      CALL distribution_1d_release(local_molecules)

      CALL nnp_init_model(nnp_env=nnp_env, printtag="NNP")

      CALL timestop(handle)

   END SUBROUTINE nnp_init_subsys

! **************************************************************************************************
!> \brief Initialize the Neural Network Potential
!> \param nnp_env ...
!> \param printtag ...
!> \param gpu_capable .FALSE. for callers whose force evaluation has no device
!>                    path, such as the path-integral helium-solute coupling.
!>                    The GPU keywords are then not consulted at all and the
!>                    environment stays on the CPU, so it can neither claim the
!>                    per-process device nor make a &FORCE_EVAL NNP lose it.
!>                    Defaults to .TRUE.
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_init_model(nnp_env, printtag, gpu_capable)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
      CHARACTER(LEN=*), INTENT(IN)                       :: printtag
      LOGICAL, INTENT(IN), OPTIONAL                      :: gpu_capable

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_init_model'
      INTEGER, PARAMETER                                 :: def_str_len = 256

      CHARACTER(len=1), ALLOCATABLE, DIMENSION(:)        :: cactfnct
      CHARACTER(len=2)                                   :: ele
      CHARACTER(len=def_str_len)                         :: dummy, line
      CHARACTER(len=default_path_length)                 :: file_name
      INTEGER                                            :: handle, i, io, j, k, nele, nuc_ele, &
                                                            symfnct_type, unit_nr
      LOGICAL                                            :: atom_e_found, explicit, first, found, &
                                                            my_gpu_capable
      REAL(KIND=dp)                                      :: energy
      REAL(KIND=dp), DIMENSION(7)                        :: test_array
      REAL(KIND=dp), DIMENSION(:), POINTER               :: work
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_parser_type)                               :: parser
      TYPE(section_vals_type), POINTER                   :: bias_section, model_section

      CALL timeset(routineN, handle)

      NULLIFY (logger)

      logger => cp_get_default_logger()

      unit_nr = cp_logger_get_default_unit_nr(logger)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) TRIM(printtag)//"| Neural Network Potential Force Environment"
      END IF

      model_section => section_vals_get_subs_vals(nnp_env%nnp_input, "MODEL")
      CALL section_vals_get(model_section, n_repetition=nnp_env%n_committee)
      ALLOCATE (nnp_env%atomic_energy(nnp_env%num_atoms, nnp_env%n_committee))
      ALLOCATE (nnp_env%committee_energy(nnp_env%n_committee))
      ALLOCATE (nnp_env%myforce(3, nnp_env%num_atoms, nnp_env%n_committee))
      ALLOCATE (nnp_env%committee_forces(3, nnp_env%num_atoms, nnp_env%n_committee))
      ALLOCATE (nnp_env%committee_stress(3, 3, nnp_env%n_committee))

      CALL section_vals_val_get(nnp_env%nnp_input, "NNP_INPUT_FILE_NAME", c_val=file_name)
      CALL section_vals_val_get(nnp_env%nnp_input, "RAD_SPLINE_N", i_val=nnp_env%rad_spline_n)
      IF (nnp_env%rad_spline_n < 2) THEN
         CPABORT("NNP| RAD_SPLINE_N must be >= 2.")
      END IF
      CALL section_vals_val_get(nnp_env%nnp_input, "VERLET_SKIN", r_val=nnp_env%verlet_skin)

      my_gpu_capable = .TRUE.
      IF (PRESENT(gpu_capable)) my_gpu_capable = gpu_capable

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

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) TRIM(printtag)//"| Reading NNP input from file: ", TRIM(file_name)
      END IF

      CALL parser_search_string(parser, "number_of_elements", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, nnp_env%n_ele
      ELSE
         CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                       "| number of elements missing in NNP_INPUT_FILE")
      END IF

      CALL parser_search_string(parser, "scale_symmetry_functions_sigma", .TRUE., found, &
                                search_from_begin_of_file=.TRUE.)
      nnp_env%scale_sigma_acsf = found

      CALL parser_search_string(parser, "scale_symmetry_functions", .TRUE., found, &
                                search_from_begin_of_file=.TRUE.)
      nnp_env%scale_acsf = found

      ! parser_search_string matches substrings, so a search for
      ! "scale_symmetry_functions" also fires on "..._sigma". Re-search from the
      ! current parser position to confirm the bare keyword appears on its own line.
      CALL parser_search_string(parser, "scale_symmetry_functions", .TRUE., found)
      IF (found .AND. nnp_env%scale_sigma_acsf) THEN
         CPWARN("Ignoring sigma ACSF scaling; both keywords were set.")
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
      ! Print warning if centering and scaling is requested:
      IF (nnp_env%center_acsf .AND. nnp_env%scale_acsf) THEN
         IF ((ABS(nnp_env%scmin) > EPSILON(0.0_dp)*1.0E+4_dp) .OR. (ABS(nnp_env%scmax - 1.0_dp) > EPSILON(0.0_dp)*1.0E+4_dp)) THEN
            CALL cp_warn(__LOCATION__, &
                         "Centering and scaling of symmetry functions requested while scale_min_short_atomic != 0 and/or "// &
                         "scale_max_short_atomic != 1. Make sure that scaling and centering of symmetry functions in CP2K "// &
                         "is consistent with your training code. "// &
                         "In CP2K: G* = (G - ave(G)) / (max(G) - min(G)) * (Smax - Smin) + Smin")
         END IF
      END IF

      CALL parser_search_string(parser, "normalize_nodes", .TRUE., found, &
                                search_from_begin_of_file=.TRUE.)
      nnp_env%normnodes = found

      CALL parser_search_string(parser, "cutoff_type", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, nnp_env%cut_type
      ELSE
         CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                       "| no cutoff type specified in NNP_INPUT_FILE")
      END IF

      CALL parser_search_string(parser, "global_hidden_layers_short", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, nnp_env%n_hlayer
      ELSE
         CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                       "| number of hidden layers missing in NNP_INPUT_FILE")
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
            CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                          "| elements not specified in NNP_INPUT_FILE")
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
               CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                             "| atom energies are not specified")
            END IF
         END DO
      END IF

      CALL parser_search_string(parser, "global_nodes_short", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, nnp_env%n_hnodes(:)
      ELSE
         CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                       "NNP| global_nodes_short not specified in NNP_INPUT_FILE")
      END IF

      CALL parser_search_string(parser, "global_activation_short", .TRUE., found, line, &
                                search_from_begin_of_file=.TRUE.)
      IF (found) THEN
         READ (line, *) dummy, cactfnct(:)
      ELSE
         CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                       "| global_activation_short not specified in NNP_INPUT_FILE")
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
            CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                          "| Activation function unkown")
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
               IF (TRIM(ele) == nnp_env%ele(i)) THEN
                  IF (symfnct_type == 2) THEN
                     nnp_env%n_rad(i) = nnp_env%n_rad(i) + 1
                  ELSE IF (symfnct_type == 3) THEN
                     nnp_env%n_ang(i) = nnp_env%n_ang(i) + 1
                  ELSE
                     CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                                   "| Symmetry function type not supported")
                  END IF
               END IF
            END DO
            first = .FALSE.
         ELSE
            IF (first) CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                                     "| no symfunction_short specified in NNP_INPUT_FILE")
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
               IF (TRIM(ele) == nnp_env%ele(i)) THEN
                  IF (symfnct_type == 2) THEN
                     nnp_env%n_rad(i) = nnp_env%n_rad(i) + 1
                     READ (line, *) dummy, ele, symfnct_type, &
                        nnp_env%rad(i)%ele(nnp_env%n_rad(i)), &
                        nnp_env%rad(i)%eta(nnp_env%n_rad(i)), &
                        nnp_env%rad(i)%rs(nnp_env%n_rad(i)), &
                        nnp_env%rad(i)%funccut(nnp_env%n_rad(i))
                     IF (nnp_env%max_cut < nnp_env%rad(i)%funccut(nnp_env%n_rad(i))) THEN
                        nnp_env%max_cut = nnp_env%rad(i)%funccut(nnp_env%n_rad(i))
                     END IF
                  ELSE IF (symfnct_type == 3) THEN
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
                     CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                                   "| Symmetry function type not supported")
                  END IF
               END IF
            END DO
            first = .FALSE.
         ELSE
            IF (first) CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                                     "| no symfunction_short specified in NNP_INPUT_FILE")
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
            IF (nnp_env%ang(i)%nuc_ele1(j) > nnp_env%ang(i)%nuc_ele2(j)) THEN
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
      CALL nnp_write_acsf(nnp_env, logger%para_env, printtag)
      CALL nnp_write_arc(nnp_env, logger%para_env, printtag)

      ! read scaling information from file
      IF (nnp_env%scale_acsf .OR. nnp_env%center_acsf .OR. nnp_env%scale_sigma_acsf) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) TRIM(printtag)//"| Reading scaling information from file: ", TRIM(file_name)
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

         ! Reject degenerate scaling.data at load time. The (loc_max - loc_min)
         ! and sigma denominators are used unguarded inside nnp_scale_acsf;
         ! a zero range or zero sigma there produces silent NaN forces.
         IF (nnp_env%scale_acsf) THEN
            DO i = 1, nnp_env%n_ele
               DO j = 1, nnp_env%n_rad(i)
                  IF (nnp_env%rad(i)%loc_max(j) <= nnp_env%rad(i)%loc_min(j)) THEN
                     WRITE (line, '(A,I0,A,I0,A,2(1X,ES13.6))') &
                        "scaling.data: radial sf range non-positive for element ", i, &
                        " sf ", j, " (loc_min, loc_max) =", &
                        nnp_env%rad(i)%loc_min(j), nnp_env%rad(i)%loc_max(j)
                     CPABORT(TRIM(line))
                  END IF
               END DO
               DO j = 1, nnp_env%n_ang(i)
                  IF (nnp_env%ang(i)%loc_max(j) <= nnp_env%ang(i)%loc_min(j)) THEN
                     WRITE (line, '(A,I0,A,I0,A,2(1X,ES13.6))') &
                        "scaling.data: angular sf range non-positive for element ", i, &
                        " sf ", j, " (loc_min, loc_max) =", &
                        nnp_env%ang(i)%loc_min(j), nnp_env%ang(i)%loc_max(j)
                     CPABORT(TRIM(line))
                  END IF
               END DO
            END DO
         END IF
         IF (nnp_env%scale_sigma_acsf) THEN
            DO i = 1, nnp_env%n_ele
               DO j = 1, nnp_env%n_rad(i)
                  IF (nnp_env%rad(i)%sigma(j) <= 0.0_dp) THEN
                     WRITE (line, '(A,I0,A,I0,A,1X,ES13.6)') &
                        "scaling.data: radial sf sigma non-positive for element ", i, &
                        " sf ", j, " sigma =", nnp_env%rad(i)%sigma(j)
                     CPABORT(TRIM(line))
                  END IF
               END DO
               DO j = 1, nnp_env%n_ang(i)
                  IF (nnp_env%ang(i)%sigma(j) <= 0.0_dp) THEN
                     WRITE (line, '(A,I0,A,I0,A,1X,ES13.6)') &
                        "scaling.data: angular sf sigma non-positive for element ", i, &
                        " sf ", j, " sigma =", nnp_env%ang(i)%sigma(j)
                     CPABORT(TRIM(line))
                  END IF
               END DO
            END DO
         END IF
      END IF

      CALL nnp_init_acsf_groups(nnp_env)

      ! Allocate per-layer weight tables, then load each committee member
      ! from disk. Allocation lives at the call site because n_committee /
      ! n_layer / n_nodes are owned here; the loader only fills them.
      DO i = 1, nnp_env%n_ele
         DO j = 2, nnp_env%n_layer
            ALLOCATE (nnp_env%arc(i)%layer(j)%weights(nnp_env%arc(i)%n_nodes(j - 1), &
                                                      nnp_env%arc(i)%n_nodes(j), nnp_env%n_committee))
            ALLOCATE (nnp_env%arc(i)%layer(j)%bweights(nnp_env%arc(i)%n_nodes(j), nnp_env%n_committee))
         END DO
      END DO
      CALL nnp_read_committee_weights(nnp_env, model_section, printtag)

      ! Resolve USE_GPU only now that the model definition is fully read: the
      ! resolution sizes the model against the compiled-in device-table
      ! capacities, which needs the group tables populated above.
      IF (my_gpu_capable) THEN
         CALL nnp_resolve_gpu(nnp_env, printtag, unit_nr, logger)
      ELSE
         ! This caller has no device force path, and its section carries no GPU
         ! keywords to read.
         nnp_env%use_gpu = nnp_use_gpu_off
         nnp_env%gpu_precision = nnp_gpu_prec_double
      END IF

      ! Upload the model-static descriptor tables (per-element scaling
      ! parameters and per-(group, element) symmetry-function maps) to the
      ! device once the model definition is fully read; the per-step descriptor
      ! kernels read them from device memory.
      IF (nnp_env%use_gpu == nnp_use_gpu_on) THEN
         CALL nnp_gpu_upload_descriptor_statics(nnp_env)
         CALL nnp_gpu_upload_network(nnp_env)
      END IF

      nnp_env%expol = 0

      ! Bias the standard deviation of committee disagreement
      NULLIFY (bias_section)
      explicit = .FALSE.
      !HELIUM NNP does not currently define a bias term
      bias_section => section_vals_get_subs_vals(nnp_env%nnp_input, "BIAS", can_return_null=.TRUE.)
      IF (ASSOCIATED(bias_section)) CALL section_vals_get(bias_section, explicit=explicit)
      nnp_env%bias = .FALSE.
      IF (explicit) THEN
         IF (nnp_env%n_committee > 1) THEN
            IF (unit_nr > 0) THEN
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
               IF (SIZE(work) /= nnp_env%n_committee) THEN
                  CPABORT("ALIGN_NNP_ENERGIES size mismatch wrt committee size.")
               END IF
               nnp_env%bias_e_avrg(:) = work
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) TRIM(printtag)//"| Biasing is aligned by shifting the energy prediction of the C-NNP members"
               END IF
            END IF
         ELSE
            CPWARN("NNP committee size is 1, BIAS section is ignored.")
         END IF
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) TRIM(printtag)//"| NNP force environment initialized"
      END IF

      CALL timestop(handle)

   END SUBROUTINE nnp_init_model

! **************************************************************************************************
!> \brief Resolve &NNP / USE_GPU and GPU_PRECISION once, before any force code
!>        runs, and claim the device when the GPU path is selected.
!>
!> AUTO probes for an offload device and otherwise falls back to the CPU path;
!> builds without NNP GPU support always report no device. The resolved value is
!> either nnp_use_gpu_on or nnp_use_gpu_off, so the force evaluation never has
!> to interpret AUTO.
!>
!> \param nnp_env  NNP environment; use_gpu and gpu_precision are set here.
!> \param printtag Log-line prefix.
!> \param unit_nr  Output unit, or a non-positive value when there is none.
!> \param logger   Default logger, for its parallel environment.
!> \author Dhruv Sharma (ds2173@cam.ac.uk)
! **************************************************************************************************
   SUBROUTINE nnp_resolve_gpu(nnp_env, printtag, unit_nr, logger)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
      CHARACTER(LEN=*), INTENT(IN)                       :: printtag
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      INTEGER, PARAMETER :: gpu_auto_min_atoms_per_rank = 1000

      CHARACTER(LEN=2*default_string_length)             :: cap_reason
      INTEGER                                            :: n_gpu_ranks, use_gpu_input
      LOGICAL                                            :: gpu_accepted, gpu_everywhere

      CALL section_vals_val_get(nnp_env%nnp_input, "USE_GPU", i_val=use_gpu_input)

      ! The outcome has to be the same on every rank. A rank-dependent answer
      ! would either hang (some ranks abort below while the rest go on into the
      ! collectives) or run a silent mix of CPU and GPU arithmetic, so count the
      ! ranks that see a device and decide from the total.
      n_gpu_ranks = MERGE(1, 0, nnp_gpu_available())
      CALL logger%para_env%sum(n_gpu_ranks)
      gpu_everywhere = (n_gpu_ranks == logger%para_env%num_pe)

      SELECT CASE (use_gpu_input)
      CASE (nnp_use_gpu_off)
         nnp_env%use_gpu = nnp_use_gpu_off
      CASE (nnp_use_gpu_on)
         IF (n_gpu_ranks == 0) THEN
            CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                          "| USE_GPU ON requested, but this build has no NNP GPU "// &
                          "support or no device is available on this node.")
         ELSE IF (.NOT. gpu_everywhere) THEN
            CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                          "| USE_GPU ON requested, but only some of the MPI ranks see "// &
                          "a device. Give every rank one, or set USE_GPU OFF.")
         END IF
         nnp_env%use_gpu = nnp_use_gpu_on
      CASE (nnp_use_gpu_auto)
         ! Below roughly gpu_auto_min_atoms_per_rank rank-local atoms the kernel
         ! launch and transfer overhead outweighs the device throughput, so AUTO
         ! keeps such systems on the CPU path rather than making an existing
         ! input slower, and its trajectory irreproducible, the moment it is run
         ! on an accelerated build. USE_GPU ON overrides the threshold.
         IF (gpu_everywhere .AND. &
             nnp_env%num_atoms >= gpu_auto_min_atoms_per_rank*logger%para_env%num_pe) THEN
            nnp_env%use_gpu = nnp_use_gpu_on
         ELSE
            nnp_env%use_gpu = nnp_use_gpu_off
         END IF
      CASE DEFAULT
         CPABORT("Unknown USE_GPU value.")
      END SELECT

      CALL section_vals_val_get(nnp_env%nnp_input, "GPU_PRECISION", &
                                i_val=nnp_env%gpu_precision)

      ! A model larger than the compiled-in device-table capacities cannot take
      ! the GPU path. Sized here, before the device is claimed and before the
      ! banner prints, so an oversized model falls back cleanly under AUTO and
      ! aborts with the limit named under ON.
      IF (nnp_env%use_gpu == nnp_use_gpu_on) THEN
         IF (.NOT. nnp_gpu_model_capacity_ok(nnp_env, cap_reason)) THEN
            nnp_env%use_gpu = nnp_use_gpu_off
            IF (use_gpu_input == nnp_use_gpu_on) THEN
               CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                             "| USE_GPU ON requested, but "//TRIM(cap_reason)// &
                             ". Set USE_GPU OFF, or raise the limit in "// &
                             "src/nnp/gpu/nnp_gpu_state.h and rebuild.")
            END IF
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) TRIM(printtag)//"| GPU offload: declined, "// &
                  TRIM(cap_reason)//"; using the CPU path"
            END IF
         END IF
      END IF

      IF (nnp_env%use_gpu /= nnp_use_gpu_on) RETURN

      ! Claim the device and push the precision selector; any failure aborts
      ! here rather than mid-MD. The claim is declined
      ! when another NNP environment in this process already holds the device,
      ! which is a hard error only if the user asked for the GPU explicitly.
      CALL nnp_gpu_init_device(nnp_env%gpu_precision, gpu_accepted)
      IF (.NOT. gpu_accepted) THEN
         nnp_env%use_gpu = nnp_use_gpu_off
         IF (use_gpu_input == nnp_use_gpu_on) THEN
            CALL cp_abort(__LOCATION__, TRIM(printtag)// &
                          "| USE_GPU ON requested, but another NNP environment in this "// &
                          "process already uses the GPU, and only one may. Set USE_GPU "// &
                          "OFF on this section.")
         END IF
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) TRIM(printtag)//"| GPU offload: declined, another NNP "// &
               "environment holds the device; using the CPU path"
         END IF
      ELSE IF (unit_nr > 0) THEN
         IF (nnp_env%gpu_precision == nnp_gpu_prec_mixed) THEN
            WRITE (unit_nr, *) TRIM(printtag)//"| GPU offload: enabled (MIXED precision)"
         ELSE
            WRITE (unit_nr, *) TRIM(printtag)//"| GPU offload: enabled (DOUBLE precision)"
         END IF
      END IF

   END SUBROUTINE nnp_resolve_gpu

! **************************************************************************************************
!> \brief Read committee weights from disk into pre-allocated arc layer tables.
!>
!> One file per (committee member, element), named <base>.<nuc>.data. The file
!> is a flat list of doubles; the order on disk is layer-major over (i->j weights
!> then biases) so a single counter walks every value into its slot.
!>
!> \param nnp_env       NNP environment with arc()%layer()%weights/bweights already allocated.
!> \param model_section &MODEL section block; supplies the per-member WEIGHTS base path.
!> \param printtag      Log-line prefix.
! **************************************************************************************************
   SUBROUTINE nnp_read_committee_weights(nnp_env, model_section, printtag)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
      TYPE(section_vals_type), POINTER                   :: model_section
      CHARACTER(LEN=*), INTENT(IN)                       :: printtag

      CHARACTER(len=default_path_length)                 :: base_name, file_name
      INTEGER                                            :: i, i_com, iweight, j, k, l, n_weight, &
                                                            unit_nr
      LOGICAL                                            :: at_end
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_parser_type)                               :: parser

      NULLIFY (logger)
      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_unit_nr(logger)

      DO i_com = 1, nnp_env%n_committee
         CALL section_vals_val_get(model_section, "WEIGHTS", c_val=base_name, i_rep_section=i_com)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) TRIM(printtag)//"| Initializing weights for model: ", i_com
         END IF
         DO i = 1, nnp_env%n_ele
            WRITE (file_name, '(A,I0.3,A)') TRIM(base_name)//".", nnp_env%nuc_ele(i), ".data"
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) TRIM(printtag)//"| Reading weights from file: ", TRIM(file_name)
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

   END SUBROUTINE nnp_read_committee_weights

END MODULE nnp_environment
