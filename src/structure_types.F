!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/structure_types [1.0] *
!!
!!   NAME
!!     structure_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     Harald Forbert (Dec-2000): Add ll-data to structure_type
!!     CJM 14-April-2001 : added coeff_type to structure_type
!!     CJM 12-SEPT-2002 : major rewrite
!!
!!   SOURCE
!******************************************************************************

MODULE structure_types

  USE coefficient_types, ONLY : coeff_type
  USE fist_environment_types, ONLY : replica_environment_type, &
                                     pimd_environment_type
  USE kinds, ONLY : dbl
  USE kg_environment_types, ONLY : kg_environment_type
  USE linklist_types, ONLY : linklist_internal_data_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE particle_types, ONLY : particle_type
  USE qs_environment_types, ONLY : qs_environment_type
  USE simulation_cell, ONLY : cell_type

  IMPLICIT NONE

  PRIVATE

  TYPE structure_type
     TYPE ( coeff_type ), POINTER  :: coef_pos ( : )        
     TYPE ( coeff_type ), POINTER  :: coef_vel ( : )        
     TYPE ( coeff_type ), POINTER  :: coef_force ( : )      
     TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
     TYPE ( particle_type ), POINTER :: part ( : )
     TYPE ( particle_node_type ), POINTER :: pnode ( : )
     TYPE ( pimd_environment_type ), POINTER :: pimd_env
  END TYPE structure_type

! *** Public subroutines and data types ***
  PUBLIC :: structure_type, set_structure_type, &
            get_structure_type, init_structure_type

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_structure_type( struc, part, pnode, molecule, &
                                pimd_env, coef_pos, coef_vel, coef_force )

    IMPLICIT NONE
!   Purpose: Get the structure type

!   ***************************************************************************

    TYPE ( structure_type ), INTENT ( IN ) :: struc
    TYPE ( particle_type ), POINTER, OPTIONAL :: part ( : )
    TYPE ( particle_node_type ), POINTER, OPTIONAL :: pnode ( : )
    TYPE ( molecule_structure_type ), POINTER, OPTIONAL :: molecule ( : )
    TYPE ( pimd_environment_type ), POINTER, OPTIONAL :: pimd_env
    TYPE ( coeff_type ), POINTER, OPTIONAL  :: coef_pos ( : )        
    TYPE ( coeff_type ), POINTER, OPTIONAL  :: coef_vel ( : )        
    TYPE ( coeff_type ), POINTER, OPTIONAL  :: coef_force ( : )      
!   ---------------------------------------------------------------------------

    IF ( PRESENT ( pnode ) ) pnode => struc % pnode
    IF ( PRESENT ( part ) ) part => struc % part
    IF ( PRESENT ( molecule ) ) molecule => struc % molecule
    IF ( PRESENT ( coef_pos ) ) coef_pos => struc % coef_pos
    IF ( PRESENT ( coef_vel ) ) coef_vel => struc % coef_vel
    IF ( PRESENT ( coef_force ) ) coef_force => struc % coef_force
    IF ( PRESENT ( pimd_env ) ) pimd_env => struc % pimd_env

  END SUBROUTINE get_structure_type

! *****************************************************************************

  SUBROUTINE init_structure_type ( struc )

!   Purpose: Initialise the FIST environment.

!   ***************************************************************************

    TYPE(structure_type), INTENT(OUT) :: struc

!   ---------------------------------------------------------------------------

    NULLIFY (struc%pnode)
    NULLIFY (struc%part)
    NULLIFY (struc%molecule)
    NULLIFY (struc%pimd_env)
    NULLIFY (struc%coef_pos)
    NULLIFY (struc%coef_vel)
    NULLIFY (struc%coef_force)

  END SUBROUTINE init_structure_type

! *****************************************************************************

  SUBROUTINE set_structure_type( struc, kg_env, rep_env, qs_env )

!   Purpose: Set the structure environment to the correct program.


!   ***************************************************************************
    TYPE ( structure_type ), INTENT ( OUT ) :: struc 
    TYPE ( kg_environment_type ), INTENT ( IN ), OPTIONAL    :: kg_env
    TYPE ( qs_environment_type ), INTENT ( IN ), OPTIONAL    :: qs_env
    TYPE ( replica_environment_type ), INTENT ( IN ), OPTIONAL  :: rep_env 

    IF ( PRESENT ( rep_env ) ) THEN 
      struc % pnode => rep_env % pnode
      struc % part => rep_env % part
      struc % molecule => rep_env % molecule
      struc % pimd_env => rep_env % pimd_env
    ENDIF
    IF ( PRESENT ( kg_env ) ) THEN 
      struc % pnode => kg_env % pnode
      struc % part => kg_env % part
      struc % molecule => kg_env % molecule
      struc % coef_pos => kg_env % coef_pos
      struc % coef_vel => kg_env % coef_vel
      struc % coef_force => kg_env % coef_force
    ENDIF
    IF ( PRESENT ( qs_env ) ) THEN 
    ENDIF

!   ---------------------------------------------------------------------------

  END SUBROUTINE set_structure_type

!***************************************************************************

END MODULE structure_types

!******************************************************************************
