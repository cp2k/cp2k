!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/initialize_extended_types [1.0] *
!!
!!   NAME
!!     initialize_extended_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE initialize_extended_types
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type
  USE molecule_types, ONLY : molecule_structure_type, molecule_type
  USE message_passing, ONLY : mp_sum
  USE nhc, ONLY : nhc_control_particles, nhc_control_barostat
  USE nose, ONLY : initial_nhc, lnhc_parameters_type, initial_npt, &
                   npt_info_type
  USE simulation_cell, ONLY : cell_type
  
  IMPLICIT NONE
  
  PUBLIC :: initialize_extended_system
  PRIVATE
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** initialize_extended_types/initialize_extended_system [1.0] *
!!
!!   NAME
!!     initialize_extended_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_extended_system ( box, simpar, molecule, mol_setup, &
     globenv, nhc_part, nhc_baro, npt_info )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( INOUT ) :: mol_setup
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( lnhc_parameters_type ),  INTENT ( INOUT ) :: nhc_part
  TYPE ( lnhc_parameters_type ),  INTENT ( INOUT ) :: nhc_baro
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  
!------------------------------------------------------------------------------
  
! fire up the thermostats, if not NVE
  IF ( simpar % ensemble /= 'NVE' ) THEN
     CALL nhc_control_particles ( simpar, molecule, mol_setup, nhc_part, globenv )
     
     IF ( simpar % read_type /= 'ALL' ) THEN
        
! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
        CALL initial_nhc ( nhc_part, simpar % tau_nhc,  &
	                   simpar % temp_ext, globenv )
     END IF
     
  END IF

  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' )  THEN
     CALL nhc_control_barostat ( simpar, nhc_baro, globenv )
     IF ( simpar % read_type /= 'ALL' ) THEN
        
! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
        CALL initial_nhc ( nhc_baro, simpar % tau_nhc,  &
	                   simpar % temp_ext, globenv )
     END IF
  ENDIF
! compute the total number of degrees of freedom for temperature
  
! Initialize simpar%nfree to dim (because at this point no
! assumptions have been made about thermostatting options
! and COM is ALWAYS zeroed in this code.  If there are
! multiple thermostats etc...simpar%nfree will be set
! to zero in nhc_control.f)
  simpar % nfree = 6 - SUM ( box % perd ( 1:3 ) )
  
  CALL get_nfree ( simpar, molecule, globenv )
  
  IF ( simpar % ensemble ( 1 : 3 ) == 'NPT' ) CALL initial_npt ( npt_info, &
       simpar % tau_cell, simpar % temp_ext, box % deth, simpar % nfree, &
       simpar % ensemble, globenv )
  
END SUBROUTINE initialize_extended_system

!!*****
!******************************************************************************
!!****** initialize_extended_types/get_nfree [1.0] *
!!
!!   NAME
!!     get_nfree
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     29-Nov-00 (JGH) correct counting of DOF if constraints are off
!!     18-Dec-00 (JGH) bug fix in initialization of nconst
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_nfree ( simpar, molecule, globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( molecule_structure_type ), INTENT ( IN ), DIMENSION ( : ) :: molecule
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  
! Locals
  INTEGER :: imol, nconst, natoms, nfree, iw
  
!------------------------------------------------------------------------------
  
  natoms = 0
  DO imol = 1, SIZE ( molecule )
     natoms = natoms + molecule ( imol ) % natoms_mol
  END DO
  nconst = 0
  IF ( simpar % constraint ) THEN
    DO imol = 1, SIZE ( molecule )
       nconst = nconst + molecule ( imol ) % ndistconst_mol
       nconst = nconst + 3 * molecule ( imol ) % ng3x3const_mol
    END DO
  END IF
  
#if defined(__parallel)
  CALL mp_sum(nconst,globenv % group)
  CALL mp_sum(natoms,globenv % group)
#endif
  
! definining simpar%nfree.  If multiple thermostats present
! it is initialized to 0.  It one thermostat or no thermostats
! it is initialized to 3
  nfree = 3 * natoms - nconst - simpar % nfree
  IF ( globenv % ionode .AND. globenv % print_level > 2 ) THEN
     iw = globenv % scr
     WRITE ( iw, '( /, A )' ) ' Calculation of degrees of freedom'
     WRITE ( iw, '( T48, A, T71, I10 )' ) '       Number of atoms:', natoms
     WRITE ( iw, '( T48, A, T71, I10 )' ) ' Number of constraints:', nconst
     WRITE ( iw, '( T34, A, T71, I10 )' ) &
          ' Invariants(translation + rotations):', simpar % nfree
     WRITE ( iw, '( T48, A, T71, I10 )' ) '    Degrees of freedom:', nfree
     WRITE ( iw, '()' )
  END IF
  simpar % nfree = nfree
  
END SUBROUTINE get_nfree

!!*****
!******************************************************************************

END MODULE initialize_extended_types

!******************************************************************************
