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
               cartesian_om,&
               cell_parameters,&
               charge_density_matrix,&
               core_hamiltonian_matrix,&
               distribution,&
               erf_ppl_radii,&
               erfc_radii,&
               exp_ppl_radii,&
               hartree_matrix,&
               interaction_thresholds,&
               interatomic_distances,&
               kind_radii,&
               kinetic_energy_matrix,&
               kohn_sham_matrix,&
               memory,&
               method_specifications,&
               mo_eigenvalues,&
               mo_eigenvectors
    LOGICAL :: nl_overlap,&
               nl_verfc,&
               nl_verf_ppl,&
               nl_vexp_ppl,&
               nl_vprj_ppnl,&
               ortho_matrix,&
               overlap_matrix,&
               pgf_radii,&
               physcon,&
               potential,&
               prj_ppnl_radii,&
               program_banner,&
               program_run_information,&
               projectors,&
               set_radii,&
               spherical_harmonics,&
               timing_information,&
               title,&
               total_numbers
  END TYPE print_key_type

! *** Public data types ***

  PUBLIC :: print_key_type

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
    print_key%cartesian_om = key_value(level,(/f,f,f/))
    print_key%cell_parameters = key_value(level,(/t,t,t/))
    print_key%charge_density_matrix = key_value(level,(/f,f,t/))
    print_key%core_hamiltonian_matrix = key_value(level,(/f,f,t/))
    print_key%distribution = key_value(level,(/f,t,t/))
    print_key%erf_ppl_radii = key_value(level,(/f,t,t/))
    print_key%erfc_radii = key_value(level,(/f,t,t/))
    print_key%exp_ppl_radii = key_value(level,(/f,t,t/))
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
    print_key%nl_overlap = key_value(level,(/f,f,t/))
    print_key%nl_verfc = key_value(level,(/f,f,t/))
    print_key%nl_verf_ppl = key_value(level,(/f,f,t/))
    print_key%nl_vexp_ppl = key_value(level,(/f,f,t/))
    print_key%nl_vprj_ppnl = key_value(level,(/f,f,t/))
    print_key%ortho_matrix = key_value(level,(/f,f,t/))
    print_key%overlap_matrix = key_value(level,(/f,f,t/))
    print_key%pgf_radii = key_value(level,(/f,f,t/))
    print_key%physcon = key_value(level,(/f,f,t/))
    print_key%potential = key_value(level,(/f,f,t/))
    print_key%prj_ppnl_radii = key_value(level,(/f,f,t/))
    print_key%program_banner = key_value(level,(/t,t,t/))
    print_key%program_run_information = key_value(level,(/t,t,t/))
    print_key%projectors = key_value(level,(/f,t,t/))
    print_key%set_radii = key_value(level,(/f,t,t/))
    print_key%spherical_harmonics = key_value(level,(/f,f,t/))
    print_key%timing_information = key_value(level,(/t,t,t/))
    print_key%title = key_value(level,(/t,t,t/))
    print_key%total_numbers = key_value(level,(/t,t,t/))

  END SUBROUTINE init_print_keys

! *****************************************************************************

  FUNCTION key_value(level,level_values) RESULT(value)

    INTEGER, INTENT(IN)               :: level
    LOGICAL, DIMENSION(3), INTENT(IN) :: level_values

    LOGICAL :: value

!   ---------------------------------------------------------------------------

    SELECT CASE (level)
    CASE (:0)
      value = .FALSE.
    CASE (1:3)
      value = level_values(level)
    CASE (4:)
      value = .TRUE.
    END SELECT

  END FUNCTION key_value

! *****************************************************************************

END MODULE print_keys

!!*****
!******************************************************************************
