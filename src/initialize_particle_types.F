!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE initialize_particle_types
  
  USE atoms_input, ONLY : system_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type
  USE molecule_types, ONLY : molecule_type
  USE particle_types, ONLY : particle_type, particle_prop_type
  USE stop_program, ONLY : stop_prg
  USE string_utilities, ONLY : str_search
  
  IMPLICIT NONE
  
  PUBLIC :: initialize_particle_type
  PRIVATE
  
CONTAINS

!******************************************************************************

SUBROUTINE initialize_particle_type ( atom_names, simpar, mol_setup, &
     ainp, pstat, part )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
  TYPE ( particle_prop_type ), INTENT ( IN ), DIMENSION ( : ), TARGET :: pstat
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( system_type ), INTENT ( INOUT ) :: ainp
  CHARACTER ( len = * ), DIMENSION ( : ), INTENT ( IN ) :: atom_names
  
!locals
  INTEGER :: natoms,iat,nmol_type,i,j,k,n,iw,natom_types
  
!------------------------------------------------------------------------------
  
  natoms = SIZE ( ainp % c, 2 )
  DO iat = 1, natoms
     part(iat) % iatom = iat
     part(iat) % r(1) = ainp % c(1,iat)
     part(iat) % r(2) = ainp % c(2,iat)
     part(iat) % r(3) = ainp % c(3,iat)
  END DO
  
  IF ( simpar%read_type=='POSVEL') THEN
     DO iat = 1, natoms
        part(iat) %v(1) = ainp%v(1,iat)
        part(iat) %v(2) = ainp%v(2,iat)
        part(iat) %v(3) = ainp%v(3,iat)
     END DO
  END IF
  
  iat = 0
  nmol_type = SIZE(mol_setup)
  natom_types = SIZE(atom_names)
  DO i = 1, nmol_type
     DO k = 1, mol_setup(i) % num_mol
        DO j = 1, mol_setup(i) % molpar % natom
           iat = iat + 1
           n = str_search ( atom_names, natom_types, &
                mol_setup ( i ) % molpar % aname ( j ) )
           part ( iat ) % prop => pstat ( n )
           
!!!!           IF ( globenv % ionode .AND. globenv % print_level > 9 ) &
!!!!                WRITE ( globenv % scr, '( A, I10, T60, A, T71, I10 )' ) &
!!!!                ' DEBUG|  Particle:', iat, ' Type:', n
        END DO
     END DO
  END DO
  
END SUBROUTINE initialize_particle_type

!******************************************************************************

END MODULE initialize_particle_types
