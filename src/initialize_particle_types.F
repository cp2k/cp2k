!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/initialize_particle_types [1.0] *
!!
!!   NAME
!!     initialize_particle_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     cjm, 20-feb-2001:  only assigns velocites for init and pos
!!     restart types
!!   SOURCE
!******************************************************************************

MODULE initialize_particle_types

  USE atoms_input,                     ONLY: system_type
  USE kinds,                           ONLY: dp
  USE md,                              ONLY: simulation_parameters_type
  USE molecule_types,                  ONLY: molecule_type
  USE particle_types,                  ONLY: particle_prop_type,&
                                             particle_type
  USE string_utilities,                ONLY: str_search

  IMPLICIT NONE

  PUBLIC :: initialize_particle_type
  PRIVATE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** initialize_particle_types/initialize_particle_type [1.0] *
!!
!!   NAME
!!     initialize_particle_type
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_particle_type ( atom_names, simpar, mol_setup, &
     ainp, pstat, part )


    CHARACTER(len=*), DIMENSION(:), &
      INTENT(IN)                             :: atom_names
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(molecule_type), DIMENSION(:), &
      INTENT(IN)                             :: mol_setup
    TYPE(system_type), INTENT(INOUT)         :: ainp
    TYPE(particle_prop_type), DIMENSION(:), &
      INTENT(IN), TARGET                     :: pstat
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part

    INTEGER                                  :: i, iat, j, k, n, natom_types, &
                                                natoms, nmol_type

!------------------------------------------------------------------------------

  natoms = SIZE ( ainp % c, 2 )
  DO iat = 1, natoms
     part(iat) % iatom = iat
     part(iat) % r(1) = ainp % c(1,iat)
     part(iat) % r(2) = ainp % c(2,iat)
     part(iat) % r(3) = ainp % c(3,iat)
  END DO

  part ( : ) % v ( 1 ) =  0.0_dp
  part ( : ) % v ( 2 ) =  0.0_dp
  part ( : ) % v ( 3 ) =  0.0_dp

  IF ( ( simpar%read_type=='POSVEL' ) .OR. ( simpar%read_type=='ALL' ) ) THEN
     DO iat = 1, natoms
        part ( iat ) % v ( 1 ) = ainp % v ( 1, iat )
        part ( iat ) % v ( 2 ) = ainp % v ( 2, iat )
        part ( iat ) % v ( 3 ) = ainp % v ( 3, iat )
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

! initialize polarization variables

  part ( : ) % abase = 0
  part ( : ) % kind = 0
  part ( : ) % grid_id = 0
  part ( : ) % first_cgf = 0
  part ( : ) % last_cgf = 0

 DO iat = 1, natoms
   NULLIFY ( part ( iat ) % coef_list )
 END DO

END SUBROUTINE initialize_particle_type

!!*****
!******************************************************************************

END MODULE initialize_particle_types

!******************************************************************************
