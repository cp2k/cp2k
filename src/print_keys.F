!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/print_keys [1.0] *
!!
!!   NAME
!!     print_keys
!!
!!   FUNCTION
!!     Declaration of the print key type and the initialization of variables of
!!     that type.
!!
!!   AUTHOR
!!     Matthias Krack (27.02.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE print_keys

  IMPLICIT NONE

  PRIVATE

  TYPE print_key_type
    INTEGER :: level
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
               density_matrix,&
               diis_information,&
               distribution,&
               each_scf_step,&
               hartree_matrix,&
               interaction_thresholds,&
               interatomic_distances,&
               kind_radii,&
               kinetic_energy_matrix,&
               kohn_sham_matrix,&
               memory,&
               method_specifications,&
               mo_eigenvalues,&
               mo_eigenvectors,&
               mo_occupation_numbers
    LOGICAL :: nl_overlap,&
               nl_verfc,&
               nl_vppl,&
               nl_vppnl,&
               ortho_matrix,&
               overlap_matrix,&
               pgf_radii,&
               physcon,&
               potential,&
               ppl_radii,&
               ppnl_radii,&
               program_banner,&
               program_run_information,&
               projectors,&
               pw_grid_information,&
               scf,&
               scf_energies,&
               set_radii,&
               spherical_harmonics,&
               timing_information,&
               title,&
               total_densities,&
               total_numbers
  END TYPE print_key_type

  INTEGER, PARAMETER :: SILENT = 0,&
                        LOW = 1,&
                        MEDIUM = 2,&
                        HIGH = 3,&
                        DEBUG = 4

! *** Public data types ***

  PUBLIC :: print_key_type

! *** Public variables ***

  PUBLIC :: DEBUG,HIGH,LOW,MEDIUM,SILENT

! *** Public subroutines ***

  PUBLIC :: init_print_keys

! *****************************************************************************

CONTAINS

! *****************************************************************************
!!****** cp2k/init_print_keys [1.0] *
!!
!!   NAME
!!     init_print_keys
!!
!!   FUNCTION
!!     Initialization of the print keys for the defined print level.
!!
!!   AUTHOR
!!     Matthias Krack (27.02.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE init_print_keys(print_key,level)

    TYPE(print_key_type), POINTER :: print_key
    INTEGER, INTENT(IN)           :: level

!   *** Local parameters ***

    LOGICAL, PARAMETER :: f = .FALSE.,&
                          t = .TRUE.

!   ---------------------------------------------------------------------------

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
    print_key%density_matrix = key_value(level,(/f,f,t/))
    print_key%diis_information = key_value(level,(/f,f,t/))
    print_key%distribution = key_value(level,(/f,t,t/))
    print_key%each_scf_step = key_value(level,(/f,f,t/))
    print_key%hartree_matrix = key_value(level,(/f,t,t/))
    print_key%interaction_thresholds = key_value(level,(/f,t,t/))
    print_key%interatomic_distances = key_value(level,(/f,f,t/))
    print_key%kind_radii = key_value(level,(/f,t,t/))
    print_key%kinetic_energy_matrix = key_value(level,(/f,f,t/))
    print_key%kohn_sham_matrix = key_value(level,(/f,f,t/))
    print_key%memory = key_value(level,(/f,t,t/))
    print_key%method_specifications = key_value(level,(/t,t,t/))
    print_key%mo_eigenvalues = key_value(level,(/f,t,t/))
    print_key%mo_eigenvectors = key_value(level,(/f,f,t/))
    print_key%mo_occupation_numbers = key_value(level,(/f,t,t/))
    print_key%nl_overlap = key_value(level,(/f,f,t/))
    print_key%nl_verfc = key_value(level,(/f,f,t/))
    print_key%nl_vppl = key_value(level,(/f,f,t/))
    print_key%nl_vppnl = key_value(level,(/f,f,t/))
    print_key%ortho_matrix = key_value(level,(/f,f,t/))
    print_key%overlap_matrix = key_value(level,(/f,f,t/))
    print_key%pgf_radii = key_value(level,(/f,f,t/))
    print_key%physcon = key_value(level,(/f,f,t/))
    print_key%potential = key_value(level,(/f,f,t/))
    print_key%ppl_radii = key_value(level,(/f,t,t/))
    print_key%ppnl_radii = key_value(level,(/f,f,t/))
    print_key%program_banner = key_value(level,(/t,t,t/))
    print_key%program_run_information = key_value(level,(/t,t,t/))
    print_key%projectors = key_value(level,(/f,t,t/))
    print_key%pw_grid_information = key_value(level,(/f,t,t/))
    print_key%scf = key_value(level,(/t,t,t/))
    print_key%scf_energies = key_value(level,(/f,t,t/))
    print_key%set_radii = key_value(level,(/f,f,t/))
    print_key%spherical_harmonics = key_value(level,(/f,f,t/))
    print_key%timing_information = key_value(level,(/t,t,t/))
    print_key%title = key_value(level,(/t,t,t/))
    print_key%total_densities = key_value(level,(/f,f,t/))
    print_key%total_numbers = key_value(level,(/t,t,t/))

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

END MODULE print_keys

!!*****
!******************************************************************************
