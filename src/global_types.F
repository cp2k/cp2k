!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****h* cp2k/global_types [1.0] *
!!
!!   NAME
!!     global_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     - print keys, basis_set_file name and potential_file_name added to the
!!       global type (27.02.2001, MK)
!!     - JGH (28.11.2001) : added pp_library_path to type
!!     - Merged with MODULE print_keys (17.01.2002, MK)
!!
!!   SOURCE
!******************************************************************************

MODULE global_types

  USE cp_para_types,                   ONLY: cp_para_env_type
  IMPLICIT NONE

  PRIVATE

! *** Public parameters ***

  PUBLIC :: DEBUG,HIGH,LOW,MEDIUM,SILENT
  PUBLIC :: BLACS_GRID_SQUARE, BLACS_GRID_ROW, BLACS_GRID_COL
  PUBLIC :: WANNIER_JACOBI, WANNIER_CRAZY, WANNIER_DIRECT, WANNIER_NONE

! *** Public data types

  PUBLIC :: global_environment_type,&
            print_key_type

! *** Public subroutines ***

  PUBLIC :: init_print_keys

! *** Print level parameters ***

  INTEGER, PARAMETER :: SILENT = 0,&
                        LOW = 1,&
                        MEDIUM = 2,&
                        HIGH = 3,&
                        DEBUG = 4

  INTEGER, PARAMETER :: BLACS_GRID_SQUARE=1, &
                        BLACS_GRID_ROW=2, &
                        BLACS_GRID_COL=3

  INTEGER, PARAMETER :: WANNIER_JACOBI=1, WANNIER_CRAZY=2, WANNIER_DIRECT=3, WANNIER_NONE=4

!!*****
!!****s* global_types/print_key_type [1.0] *
!!
!!   NAME
!!     print_key_type
!!
!!   FUNCTION
!!     
!!   NOTES
!!     -
!!
!!   SOURCE
!!
  TYPE print_key_type
    INTEGER :: level,precision,ref_count
    LOGICAL :: angles,&
               atomic_coordinates,&
               basic_data_types,&
               basis_set,&
               blacs_info,&
               cartesian_om,&
               cell_parameters,&
               charge_density_matrix,&
               core_hamiltonian_matrix,&
               core_charge_radii,&
               cputime,&
               density_matrix,&
               derivatives,&
               dft_control_parameters,&
               diis_information,&
               distribution,&
               each_scf_step,&
               forces,&
               hartree_matrix,&
               interatomic_distances,&
               kind_information,&
               kind_radii,&
               kinetic_energy_matrix,&
               kohn_sham_matrix,&
               memory,&
               mo_eigenvalues,&
               mo_eigenvectors,&
               mo_occupation_numbers
    LOGICAL :: optical_conductivity,&
               ortho_matrix,&
               overlap_matrix,&
               overlap_matrix_magnitude, &
               density_matrix_magnitude, &
               pgf_radii,&
               physcon,&
               potential,&
               potential_energy_matrix,&
               ppl_radii,&
               ppnl_radii,&
               program_banner,&
               program_run_information,&
               projectors,&
               pw_grid_information,&
               sab_orb_neighbor_lists,&
               sac_ppl_neighbor_lists,&
               sac_ppnl_neighbor_lists,&
               sbc_ppnl_neighbor_lists,&
               scf,&
               scf_energies,&
               set_radii,&
               spherical_harmonics,&
               title,&
               total_densities,&
               total_numbers
    LOGICAL :: w_matrix,&
               e_density_cube, &
               v_hartree_cube, &
               lumos, &
               homos, &
               cubes, &
               mulliken_analysis, &
               mo_orthonormality, &
               molecular_states
    LOGICAL :: ep_matrixes, ep_rho1_cube, ep_operator
    LOGICAL :: timecumul
    INTEGER    nlumos,nhomos
    LOGICAL :: localization
    INTEGER :: localization_method
    LOGICAL :: HPM
  END TYPE print_key_type

!!*****
!!****s* global_types/global_environment_type [1.0] *
!!
!!   NAME
!!     global_environment_type
!!
!!   FUNCTION
!!     
!!   NOTES
!!     -
!!
!!   SOURCE
!!
  TYPE global_environment_type
     TYPE(cp_para_env_type), POINTER  :: para_env
     TYPE(print_key_type)             :: print
     CHARACTER(LEN=200)               :: basis_set_file_name,&
                                         coord_file_name,&
                                         default_fft_library,&
                                         error_file_name,&
                                         input_file_name,&
                                         output_file_name,&
                                         potential_file_name,&
                                         pp_library_path,&
                                         program_name,&
                                         project_name,&
                                         restart_file_name
     CHARACTER(LEN=60) :: run_type
     INTEGER :: group        !! group handle for all processors
     INTEGER :: idum         !! random number seed
     INTEGER :: mepos        !! position number for local processor
     INTEGER :: num_pe       !! total number of processors
     INTEGER :: print_level  !! how much information is written
     INTEGER :: scr          !! output unit
     INTEGER :: source       !! processor id for global output node
     LOGICAL :: ionode       !! this processor performs global output
     INTEGER :: blacs_grid_layout !! will store the user preference for the blacs grid
  END TYPE global_environment_type

! *****************************************************************************

CONTAINS

!!*****
!!****f* global_types/init_print_keys [1.0] *
!!
!!   NAME
!!     init_print_keys
!!
!!   FUNCTION
!!     
!!   NOTES
!!     -
!!
!!   SOURCE
!!
  SUBROUTINE init_print_keys(print_key,level)

    TYPE(print_key_type) :: print_key
    INTEGER, INTENT(IN)  :: level

!   *** Local parameters ***

    LOGICAL, PARAMETER :: f = .FALSE.,&
                          t = .TRUE.

!   ---------------------------------------------------------------------------

    print_key%ref_count=1
    print_key%angles = key_value(level,(/f,f,t/))
    print_key%atomic_coordinates = key_value(level,(/t,t,t/))
    print_key%basic_data_types = key_value(level,(/f,f,f/))
    print_key%basis_set = key_value(level,(/f,t,t/))
    print_key%blacs_info = key_value(level,(/f,t,t/))
    print_key%cartesian_om = key_value(level,(/f,f,f/))
    print_key%cell_parameters = key_value(level,(/t,t,t/))
    print_key%charge_density_matrix = key_value(level,(/f,f,t/))
    print_key%core_hamiltonian_matrix = key_value(level,(/f,f,t/))
    print_key%core_charge_radii = key_value(level,(/f,t,t/))
    print_key%cputime = key_value(level,(/f,f,f/))
    print_key%dft_control_parameters = key_value(level,(/t,t,t/))
    print_key%density_matrix = key_value(level,(/f,f,t/))
    print_key%derivatives = key_value(level,(/f,f,f/))
    print_key%diis_information = key_value(level,(/f,f,t/))
    print_key%distribution = key_value(level,(/f,t,t/))
    print_key%each_scf_step = key_value(level,(/f,f,t/))
    print_key%ep_matrixes = key_value(level,(/f,f,t/))
    print_key%ep_rho1_cube = key_value(level,(/f,f,t/))
    print_key%ep_operator = key_value(level,(/f,f,t/))
    print_key%forces = key_value(level,(/f,t,t/))
    print_key%hartree_matrix = key_value(level,(/f,f,t/))
    print_key%interatomic_distances = key_value(level,(/f,f,t/))
    print_key%kind_information = key_value(level,(/f,t,t/))
    print_key%kind_radii = key_value(level,(/f,t,t/))
    print_key%kinetic_energy_matrix = key_value(level,(/f,f,t/))
    print_key%kohn_sham_matrix = key_value(level,(/f,f,t/))
    print_key%memory = key_value(level,(/f,f,t/))
    print_key%mo_eigenvalues = key_value(level,(/f,f,t/))
    print_key%mo_eigenvectors = key_value(level,(/f,f,t/))
    print_key%mo_occupation_numbers = key_value(level,(/f,f,t/))
    print_key%ortho_matrix = key_value(level,(/f,f,t/))
    print_key%overlap_matrix = key_value(level,(/f,f,t/))
    print_key%pgf_radii = key_value(level,(/f,f,t/))
    print_key%physcon = key_value(level,(/f,t,t/))
    print_key%potential = key_value(level,(/f,f,t/))
    print_key%ppl_radii = key_value(level,(/f,t,t/))
    print_key%ppnl_radii = key_value(level,(/f,f,t/))
    print_key%program_banner = key_value(level,(/t,t,t/))
    print_key%program_run_information = key_value(level,(/t,t,t/))
    print_key%projectors = key_value(level,(/f,t,t/))
    print_key%potential_energy_matrix = key_value(level,(/f,f,f/))
    print_key%pw_grid_information = key_value(level,(/f,t,t/))
    print_key%sab_orb_neighbor_lists = key_value(level,(/f,f,t/))
    print_key%sac_ppl_neighbor_lists = key_value(level,(/f,f,t/))
    print_key%sac_ppnl_neighbor_lists = key_value(level,(/f,f,t/))
    print_key%sbc_ppnl_neighbor_lists = key_value(level,(/f,f,t/))
    print_key%scf = key_value(level,(/t,t,t/))
    print_key%scf_energies = key_value(level,(/f,t,t/))
    print_key%set_radii = key_value(level,(/f,t,t/))
    print_key%spherical_harmonics = key_value(level,(/f,f,t/))
    print_key%title = key_value(level,(/t,t,t/))
    print_key%total_densities = key_value(level,(/f,t,t/))
    print_key%total_numbers = key_value(level,(/t,t,t/))
    print_key%w_matrix = key_value(level,(/f,f,f/))
    print_key%density_matrix_magnitude = key_value(level,(/f,f,f/))
    print_key%overlap_matrix_magnitude = key_value(level,(/f,f,f/))
    print_key%e_density_cube = key_value(level,(/f,f,f/))
    print_key%v_hartree_cube = key_value(level,(/f,f,f/))
    print_key%cubes = key_value(level,(/f,f,f/))
    print_key%lumos = key_value(level,(/f,f,f/))
    print_key%mo_orthonormality = key_value(level,(/f,t,t/))
    print_key%optical_conductivity = key_value(level,(/f,f,f/))
    print_key%homos = key_value(level,(/f,f,f/))
    print_key%nlumos = 4
    print_key%nhomos = 4
    print_key%timecumul = key_value(level,(/f,f,f/))
    print_key%localization = key_value(level,(/f,f,f/))
    print_key%localization_method = WANNIER_CRAZY
    print_key%mulliken_analysis = key_value(level,(/f,t,t/))
    print_key%molecular_states = key_value(level,(/f,f,f/))
    print_key%HPM = key_value(level,(/f,f,f/))

  END SUBROUTINE init_print_keys

! *****************************************************************************

  FUNCTION key_value(level,level_values) RESULT(value)

    INTEGER, INTENT(IN)               :: level
    LOGICAL, DIMENSION(3), INTENT(IN) :: level_values

    LOGICAL :: value

!   ---------------------------------------------------------------------------

    SELECT CASE (level)
    CASE (:SILENT)
      value = .FALSE.
    CASE (LOW:HIGH)
      value = level_values(level)
    CASE (DEBUG:)
      value = .TRUE.
    END SELECT

  END FUNCTION key_value

! *****************************************************************************

END MODULE global_types
!!*****
