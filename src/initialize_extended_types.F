!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE initialize_extended_types
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type
  USE molecule_types, ONLY : molecule_structure_type, molecule_type
  USE message_passing, ONLY : mp_sum
  USE nhc, ONLY : nhc_control
  USE nose, ONLY : initial_nhc, extended_parameters_type, initial_npt
  USE simulation_cell, ONLY : cell_type
  
  IMPLICIT NONE
  
  PUBLIC :: initialize_extended_type
  PRIVATE
  
CONTAINS

!******************************************************************************

SUBROUTINE initialize_extended_type ( box, simpar, molecule, mol_setup, nhcp, &
     globenv )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( INOUT ) :: mol_setup
  TYPE ( extended_parameters_type ),  INTENT ( INOUT ) :: nhcp
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  
!------------------------------------------------------------------------------
  
! fire up the thermostats, if not NVE
  IF ( simpar % ensemble /= 'NVE' ) THEN
     CALL nhc_control ( simpar, molecule, mol_setup, nhcp, globenv )
     
     IF ( simpar % read_type /= 'ALL' ) THEN
        
! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
        CALL initial_nhc ( nhcp, simpar % tau_nhc, simpar % temp_ext, globenv )
     END IF
     
  END IF
! compute the total number of degrees of freedom for temperature
  
! Initialize simpar%nfree to dim (because at this point no
! assumptions have been made about thermostatting options
! and COM is ALWAYS zeroed in this code.  If there are
! multiple thermostats etc...simpar%nfree will be set
! to zero in nhc_control.f)
  simpar % nfree = 6 - SUM ( box % perd ( 1:3 ) )
  
  CALL get_nfree ( simpar, molecule, globenv )
  
  IF ( simpar % ensemble ( 1:3 ) == 'NPT' ) CALL initial_npt ( nhcp, &
       simpar % tau_cell, simpar % temp_ext, box % deth, simpar % nfree, &
       simpar % ensemble, globenv )
  
END SUBROUTINE initialize_extended_type

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
  
  nconst = 0
  natoms = 0
  DO imol = 1, SIZE ( molecule )
     nconst = nconst + molecule ( imol ) % ndistconst_mol
     nconst = nconst + 3 * molecule ( imol ) % ng3x3const_mol
     natoms = natoms + molecule ( imol ) % natoms_mol
  END DO
  
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

!******************************************************************************

END MODULE initialize_extended_types
