!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_environment_types [1.0] *
!!
!!   NAME
!!     kg_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gt  SEPT-24-02 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_environment_types

  USE atomic_kind_types, ONLY: atomic_kind_type
  USE coefficient_types, ONLY : coeff_type
  USE dft_types, ONLY : dft_control_type
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE fist_environment_types, ONLY : pimd_environment_type
  USE kinds, ONLY : dbl
  USE linklist_types, ONLY : linklist_internal_data_type
  USE md, ONLY : virial_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE pair_potential, ONLY : potentialparm_type
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type

  PRIVATE

  character(len=*), private, parameter :: moduleN='kg_environment_types'

  TYPE kg_thermodynamic_type
     TYPE ( virial_type ), POINTER ::  virial
     REAL ( dbl ) :: kin, pot, e_gspace, e_self, e_neut, e_xc, e_hartree,&
                     e_thomas_fermi
  END TYPE kg_thermodynamic_type

  TYPE kg_replica_environment_type
    TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
    TYPE ( cell_type ), POINTER :: box
    TYPE ( cell_type ), POINTER :: box_ref
    TYPE ( coeff_type ), POINTER  :: coef_pos ( : )
    TYPE ( coeff_type ), POINTER  :: coef_vel ( : )
    TYPE ( coeff_type ), POINTER  :: coef_force ( : )
    TYPE ( linklist_internal_data_type), POINTER :: ll_data ( : )
    TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
    TYPE ( particle_type ), POINTER :: part ( : )
    TYPE ( particle_node_type ), POINTER :: pnode ( : )
    TYPE ( pimd_environment_type ), POINTER :: pimd_env
  END TYPE kg_replica_environment_type

  TYPE kg_environment_type
    LOGICAL, POINTER :: pimd
    TYPE (dft_control_type), POINTER :: dft_control
    TYPE ( ewald_parameters_type ), POINTER :: ewald_param
    TYPE ( kg_replica_environment_type ), POINTER :: rep_env ( : )
    TYPE ( kg_thermodynamic_type ), POINTER :: thermo
    TYPE ( potentialparm_type ), POINTER :: potparm ( :, : )
  END TYPE kg_environment_type

! *** Public data types ***

  PUBLIC :: kg_environment_type,        & 
            kg_replica_environment_type,     & 
            kg_thermodynamic_type

! *** Public subroutines ***

  PUBLIC :: get_kg_env,&
            init_kg_env,&
            set_kg_env

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_kg_env( kg_env, pimd, dft_control, ewald_param, &
                         rep_env, potparm, thermo )
    IMPLICIT NONE

!   Purpose: Get the KG environment.

!   ***************************************************************************

    TYPE ( kg_environment_type ), INTENT(IN)              :: kg_env
    LOGICAL, OPTIONAL, POINTER                              :: pimd
    TYPE (dft_control_type), OPTIONAL, POINTER              :: dft_control
    TYPE ( ewald_parameters_type ), OPTIONAL, POINTER       :: ewald_param
    TYPE ( kg_replica_environment_type ), OPTIONAL, POINTER :: rep_env ( : )
    TYPE ( potentialparm_type ), OPTIONAL, POINTER          :: potparm ( :, : )
    TYPE ( kg_thermodynamic_type ), OPTIONAL, POINTER       :: thermo

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( pimd ) ) pimd => kg_env % pimd
    IF ( PRESENT ( dft_control ) ) dft_control => kg_env % dft_control
    IF ( PRESENT ( ewald_param ) ) ewald_param => kg_env % ewald_param
    IF ( PRESENT ( potparm ) ) potparm => kg_env % potparm
    IF ( PRESENT ( thermo ) ) thermo => kg_env % thermo
    IF ( PRESENT ( rep_env ) ) rep_env => kg_env % rep_env

  END SUBROUTINE get_kg_env

! *****************************************************************************

  SUBROUTINE init_kg_env ( kg_env )

!   Purpose: Initialise the KG environment.

!   ***************************************************************************

    TYPE(kg_environment_type), INTENT(OUT) :: kg_env

!   ---------------------------------------------------------------------------
    NULLIFY (kg_env%pimd)
    NULLIFY (kg_env%dft_control)
    NULLIFY (kg_env%ewald_param)
    NULLIFY (kg_env%potparm)
    NULLIFY (kg_env%thermo)
    NULLIFY (kg_env%rep_env)

  END SUBROUTINE init_kg_env

! *****************************************************************************

  SUBROUTINE set_kg_env ( kg_env, pimd, dft_control, ewald_param, &
                         rep_env, potparm, thermo )

    IMPLICIT NONE
!   Purpose: Set the kg environment.

    TYPE ( kg_environment_type ), INTENT(out)             :: kg_env
    LOGICAL, OPTIONAL, POINTER                              :: pimd
    TYPE (dft_control_type), OPTIONAL, POINTER              :: dft_control
    TYPE ( ewald_parameters_type ), OPTIONAL, POINTER       :: ewald_param
    TYPE ( potentialparm_type ), OPTIONAL, POINTER          :: potparm ( :, : )
    TYPE ( kg_thermodynamic_type ), OPTIONAL, POINTER     :: thermo
    TYPE ( kg_replica_environment_type ), OPTIONAL, POINTER    :: rep_env ( : )

    IF ( PRESENT ( pimd ) ) kg_env % pimd => pimd
    IF ( PRESENT ( dft_control ) ) kg_env % dft_control => dft_control
    IF ( PRESENT ( ewald_param ) ) kg_env % ewald_param => ewald_param
    IF ( PRESENT ( potparm ) ) kg_env % potparm => potparm
    IF ( PRESENT ( thermo ) ) kg_env % thermo => thermo
    IF ( PRESENT ( rep_env ) ) kg_env % rep_env => rep_env

!***************************************************************************
  END SUBROUTINE set_kg_env
!***************************************************************************

END MODULE kg_environment_types
