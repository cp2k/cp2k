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
  USE kinds, ONLY : dbl
  USE linklist_types, ONLY : linklist_internal_data_type
  USE md, ONLY : virial_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type

  PRIVATE

  character(len=*), private, parameter :: moduleN='kg_environment_types'

  TYPE kg_thermodynamic_type
     TYPE ( virial_type ), POINTER ::  virial
     REAL ( dbl ) :: kin, pot, e_rspace_core, e_self, e_neut, e_xc, e_hartree,&
                     e_thomas_fermi
  END TYPE kg_thermodynamic_type

  TYPE kg_environment_type
    TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
    TYPE ( cell_type ), POINTER :: box
    TYPE ( cell_type ), POINTER :: box_ref
    TYPE ( coeff_type ), POINTER  :: coef_pos ( : )
    TYPE ( coeff_type ), POINTER  :: coef_vel ( : )
    TYPE ( coeff_type ), POINTER  :: coef_force ( : )
    TYPE (dft_control_type), POINTER :: dft_control
    TYPE ( linklist_internal_data_type), POINTER :: ll_data ( : )
    TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
    TYPE ( particle_type ), POINTER :: part ( : )
    TYPE ( particle_node_type ), POINTER :: pnode ( : )
    TYPE ( kg_thermodynamic_type ), POINTER :: thermo
  END TYPE kg_environment_type

! *** Public data types ***

  PUBLIC :: kg_environment_type,        & 
            kg_thermodynamic_type

! *** Public subroutines ***

  PUBLIC :: get_kg_env,&
            init_kg_env,&
            set_kg_env

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_kg_env( kg_env, atomic_kind_set, box, box_ref, coef_pos, &
                         coef_vel, coef_force, dft_control, ll_data, &
                         molecule, part, pnode, thermo )
    IMPLICIT NONE

!   Purpose: Get the KG environment.

!   ***************************************************************************

    TYPE ( kg_environment_type ), TARGET              :: kg_env
    TYPE ( atomic_kind_type ), OPTIONAL, POINTER :: atomic_kind_set ( : )
    TYPE ( cell_type ), OPTIONAL, POINTER :: box
    TYPE ( cell_type ), OPTIONAL, POINTER :: box_ref
    TYPE ( coeff_type ), OPTIONAL, POINTER  :: coef_pos ( : )
    TYPE ( coeff_type ), OPTIONAL, POINTER  :: coef_vel ( : )
    TYPE ( coeff_type ), OPTIONAL, POINTER  :: coef_force ( : )
    TYPE ( dft_control_type ), OPTIONAL, POINTER     :: dft_control
    TYPE ( linklist_internal_data_type), OPTIONAL, POINTER :: ll_data ( : )
    TYPE ( molecule_structure_type ), OPTIONAL, POINTER :: molecule ( : )
    TYPE ( particle_type ), OPTIONAL, POINTER :: part ( : )
    TYPE ( particle_node_type ), OPTIONAL, POINTER :: pnode ( : )
    TYPE ( kg_thermodynamic_type ), OPTIONAL, POINTER       :: thermo

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( atomic_kind_set ) ) atomic_kind_set => kg_env % atomic_kind_set
    IF ( PRESENT ( box ) ) box => kg_env % box
    IF ( PRESENT ( box_ref ) ) box_ref  => kg_env % box_ref
    IF ( PRESENT ( coef_pos ) ) coef_pos  => kg_env % coef_pos
    IF ( PRESENT ( coef_vel ) ) coef_vel  => kg_env % coef_vel
    IF ( PRESENT ( coef_force ) ) coef_force  => kg_env % coef_force
    IF ( PRESENT ( dft_control ) ) dft_control => kg_env % dft_control
    IF ( PRESENT ( ll_data ) ) ll_data => kg_env % ll_data
    IF ( PRESENT ( molecule ) ) molecule => kg_env % molecule
    IF ( PRESENT ( part ) ) part => kg_env % part
    IF ( PRESENT ( pnode ) ) pnode => kg_env % pnode
    IF ( PRESENT ( thermo ) ) thermo => kg_env % thermo

  END SUBROUTINE get_kg_env

! *****************************************************************************

  SUBROUTINE init_kg_env ( kg_env )

!   Purpose: Initialise the KG environment.

!   ***************************************************************************

    TYPE(kg_environment_type), POINTER :: kg_env

!   ---------------------------------------------------------------------------
    NULLIFY (kg_env%atomic_kind_set)
    NULLIFY (kg_env%box)
    NULLIFY (kg_env%box_ref)
    NULLIFY (kg_env%coef_pos)
    NULLIFY (kg_env%coef_vel)
    NULLIFY (kg_env%coef_force)
    NULLIFY (kg_env%dft_control)
    NULLIFY (kg_env%ll_data)
    NULLIFY (kg_env%molecule)
    NULLIFY (kg_env%part)
    NULLIFY (kg_env%pnode)
    NULLIFY (kg_env%thermo)

  END SUBROUTINE init_kg_env

! *****************************************************************************

  SUBROUTINE set_kg_env( kg_env, atomic_kind_set, box, box_ref, coef_pos, &
                         coef_vel, coef_force, dft_control, ll_data, &
                         molecule, part, pnode, thermo )
    IMPLICIT NONE

!   Purpose: set the KG environment.

!   ***************************************************************************

    TYPE ( kg_environment_type ),  POINTER  :: kg_env
    TYPE ( atomic_kind_type ), OPTIONAL, TARGET :: atomic_kind_set ( : )
    TYPE ( cell_type ), OPTIONAL, TARGET :: box
    TYPE ( cell_type ), OPTIONAL, TARGET :: box_ref
    TYPE ( coeff_type ), OPTIONAL, TARGET  :: coef_pos ( : )
    TYPE ( coeff_type ), OPTIONAL, TARGET  :: coef_vel ( : )
    TYPE ( coeff_type ), OPTIONAL, TARGET  :: coef_force ( : )
    TYPE (dft_control_type), OPTIONAL, TARGET     :: dft_control
    TYPE ( linklist_internal_data_type), OPTIONAL, TARGET :: ll_data ( : )
    TYPE ( molecule_structure_type ), OPTIONAL, TARGET :: molecule ( : )
    TYPE ( particle_type ), OPTIONAL, TARGET :: part ( : )
    TYPE ( particle_node_type ), OPTIONAL, TARGET :: pnode ( : )
    TYPE ( kg_thermodynamic_type ), OPTIONAL, TARGET  :: thermo

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( atomic_kind_set ) ) kg_env % atomic_kind_set => atomic_kind_set
    IF ( PRESENT ( box ) ) kg_env % box =>  box
    IF ( PRESENT ( box_ref ) ) kg_env % box_ref  =>  box_ref
    IF ( PRESENT ( coef_pos ) ) kg_env % coef_pos  => coef_pos
    IF ( PRESENT ( coef_vel ) ) kg_env % coef_vel  =>  coef_vel
    IF ( PRESENT ( coef_force ) ) kg_env % coef_force  =>  coef_force
    IF ( PRESENT ( dft_control ) ) kg_env % dft_control => dft_control
    IF ( PRESENT ( ll_data ) ) kg_env % ll_data => ll_data
    IF ( PRESENT ( molecule ) ) kg_env % molecule => molecule
    IF ( PRESENT ( part ) ) kg_env % part => part
    IF ( PRESENT ( pnode ) ) kg_env % pnode => pnode
    IF ( PRESENT ( thermo ) ) kg_env % thermo => thermo

!***************************************************************************
  END SUBROUTINE set_kg_env
!***************************************************************************

END MODULE kg_environment_types

!***************************************************************************
