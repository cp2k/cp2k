!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/initialize_molecule_types [1.0] *
!!
!!   NAME
!!     initialize_molecule_types
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

MODULE initialize_molecule_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error,&
                                             cp_unimplemented_error_nr
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_create,&
                                             distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE linklists,                       ONLY: atom_link_list,&
                                             bend_link_list,&
                                             bond_link_list,&
                                             distconst_link_list,&
                                             exclusion,&
                                             g3x3_link_list,&
                                             g4x6_link_list,&
                                             torsion_link_list
  USE md,                              ONLY: simulation_parameters_type
  USE molecule_types,                  ONLY: intra_parameters_type,&
                                             molecule_structure_type,&
                                             molecule_type,&
                                             particle_node_type
  USE particle_types,                  ONLY: particle_type
  USE util,                            ONLY: get_limit
  IMPLICIT NONE

  PUBLIC :: initialize_molecule_type, mol_local_particles_create
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="initialize_molecule_type"
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** initialize_molecule_types/initialize_molecule_type [1.0] *
!!
!!   NAME
!!     initialize_molecule_types
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

SUBROUTINE initialize_molecule_type ( mol_setup, pnode, part, molecule,&
                                      globenv, intra_param )

  IMPLICIT NONE
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( IN ) :: mol_setup
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ), TARGET :: part
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE ( intra_parameters_type ), INTENT ( IN ), OPTIONAL :: intra_param

! Locals
  INTEGER :: imol, i, j, k, atombase, iat, nlim ( 2 ), iw, nmol_type
  INTEGER, DIMENSION(:), ALLOCATABLE :: i_att

!------------------------------------------------------------------------------

  iw = globenv % scr
  IF ( globenv % num_pe == 1 ) THEN
     imol = 0
     atombase = 1
     nmol_type = SIZE ( mol_setup )
     DO i = 1, nmol_type
        DO k = 1, mol_setup ( i ) %num_mol
           imol = imol + 1
           molecule(imol) % mol_name = mol_setup(i) % molpar% name
           molecule(imol) % imol = imol
           molecule(imol) % atombase = atombase
           molecule(imol) % moltype = i
           molecule(imol) % natoms_mol = mol_setup(i) % molpar% natom
           molecule(imol) % nbonds_mol = mol_setup(i) % molpar% nbonds
           molecule(imol) % nbends_mol = mol_setup(i) % molpar% nbends
           molecule(imol) % ntorsions_mol = mol_setup(i) % molpar % ntorsions
           molecule(imol) % ndistconst_mol = mol_setup(i) % molpar % ndcon
           molecule(imol) % ng3x3const_mol = mol_setup(i) % molpar % n3x3con
           molecule(imol) % ng4x6const_mol = mol_setup(i) % molpar % n4x6con
           molecule(imol) % nstates = 0
           NULLIFY ( molecule ( imol ) % states )
           IF (globenv%ionode .AND. globenv%print_level>9) WRITE (iw, &
                '(A,I10,T40,A,I5,A,T61,A20 )' ) ' Molecule:', imol, &
                ' Type:', i, '     Name:', &
                ADJUSTR ( mol_setup(i) % molpar % name)
           atombase = atombase + mol_setup(i) % molpar % natom
        END DO
     END DO
     DO iat = 1, SIZE(pnode)
        pnode(iat) %p => part(iat)
     END DO

  ELSE

     imol = 0
     atombase = 1
     nmol_type = SIZE ( mol_setup )
     DO i = 1, nmol_type
        nlim = get_limit ( mol_setup ( i ) % num_mol, globenv % num_pe, &
             globenv % mepos )
        j = SUM( mol_setup(1:i) % num_mol )
        DO k = nlim(1), nlim(2)
           imol = imol + 1
           molecule(imol) % mol_name = mol_setup(i) % molpar%name
           molecule(imol) % imol = j + k
           molecule(imol) % atombase = atombase + &
                (k-1)*mol_setup(i) % molpar % natom
           molecule(imol) % moltype = i
           molecule(imol) % natoms_mol = mol_setup(i) % molpar%natom
           molecule(imol) % nbonds_mol = mol_setup(i) % molpar%nbonds
           molecule(imol) % nbends_mol = mol_setup(i) % molpar%nbends
           molecule(imol) % ntorsions_mol = mol_setup(i) % molpar%ntorsions
           molecule(imol) % ndistconst_mol = mol_setup(i) % molpar%ndcon
           molecule(imol) % ng3x3const_mol = mol_setup(i) % molpar%n3x3con
           molecule(imol) % ng4x6const_mol = mol_setup(i) % molpar%n4x6con
           molecule(imol) % nstates = 0
           NULLIFY ( molecule ( imol ) % states )
           IF (globenv%ionode .AND. globenv%print_level>9) WRITE (iw, &
                '(A,I10,T40,A,I5,A,T61,A20 )' ) ' Molecule:', imol, &
                ' Type:', i, '     Name:', ADJUSTR ( mol_setup(i) % molpar%name)
        END DO
        atombase = atombase &
             + mol_setup(i) % num_mol * mol_setup(i) % molpar % natom
     END DO

     k = 0
     DO imol = 1, SIZE(molecule)
        atombase = molecule(imol) % atombase - 1
        DO j = 1, molecule(imol) % natoms_mol
           k = k + 1
           pnode ( k ) % p => part ( atombase + j )
        END DO
     END DO
  END IF

!..atom link list
  CALL atom_link_list(molecule,pnode)

  IF ( PRESENT ( intra_param ) ) THEN
!..initialize bonds
    CALL bond_link_list(pnode,molecule,mol_setup,intra_param % bond_param)

!..initialize bends
    CALL bend_link_list(pnode,molecule,mol_setup,intra_param % bend_param)

!..initialize torsions
    CALL torsion_link_list(pnode,molecule,mol_setup,intra_param % torsion_param)

!..initialize distance constraints
    CALL distconst_link_list(pnode,molecule,mol_setup, &
       intra_param % constraint_distance)

!..initialize g3x3 constraints
   CALL g3x3_link_list(pnode,molecule,mol_setup, &
        intra_param % constraint_distance)

!..initialize g4x6 constraints
   CALL g4x6_link_list(pnode,molecule,mol_setup, &
        intra_param % constraint_distance)

!..create exclusion list
   CALL exclusion(molecule,pnode)
  ENDIF

END SUBROUTINE initialize_molecule_type

!!*****
!******************************************************************************

!!****f* initialize_molecule_type/mol_local_particles_create [1.0] *
!!
!!   NAME
!!     mol_local_particles_create
!!
!!   SYNOPSIS
!!     Subroutine mol_local_particles_create(local_particles,&
!!         atomic_kind_set, particles, local_molecules, para_env, error)
!!       Type(distribution_1d_type), Pointer:: local_particles
!!       Type(atomic_kind_type), Dimension(:), Pointer:: atomic_kind_set
!!       Type(molecule_structure_type), Dimension(:), Pointer::&
!!         local_molecules
!!       Type(particle_type), Dimension(:), Pointer:: particles
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine mol_local_particles_create
!!
!!   FUNCTION
!!     create a local_particles structure from a molecule structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - local_particles: the object to create
!!     - atomic_kind_set: tha atomic kind set object
!!     - local_molecule: structure of the molecules that are local to this
!!       processor
!!     - particles: array with all the particles
!!     - para_env: parallel einvironment on which the distribution is based
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     06.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE mol_local_particles_create(local_particles, atomic_kind_set,&
     particles, local_molecules, para_env, error)
  TYPE(distribution_1d_type), POINTER :: local_particles
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(molecule_structure_type), DIMENSION(:), POINTER :: local_molecules
  TYPE(particle_type), DIMENSION(:), POINTER :: particles
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  TYPE(atomic_kind_type), POINTER :: atomic_kind
  CHARACTER(len=*), PARAMETER :: routineN='mol_local_particles_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, ikind, imol,j, atombase
  INTEGER, DIMENSION(:), ALLOCATABLE :: multiindex
  INTEGER, DIMENSION(:), POINTER :: atom_list

  failure=.FALSE.
  NULLIFY(atom_list)
  
  CPPrecondition(.NOT.ASSOCIATED(local_particles),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(local_molecules),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(multiindex(SIZE(atomic_kind_set)),stat=stat)

     IF (para_env%num_pe==1) THEN
        DO ikind=1,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind,natom=multiindex(ikind))
        END DO

        CALL distribution_1d_create(local_particles,para_env=para_env,n_el=multiindex,&
             listbased_distribution=.FALSE.,error=error)

        DO ikind=1,SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind,atom_list=atom_list)
          local_particles%list(ikind)%array=atom_list
        END DO

     ELSE
        ! not so efficent at the moment

        multiindex=0
        !$omp parallel do private(imol,j,atombase,ikind)
        DO imol = 1, SIZE(local_molecules)
           atombase = local_molecules(imol) % atombase - 1
           DO j = 1, local_molecules(imol) % natoms_mol
              CALL get_atomic_kind(particles( atombase + j )%atomic_kind,&
                   kind_number=ikind)
              multiindex(ikind)=multiindex(ikind)+1
           END DO
        END DO

        CALL distribution_1d_create(local_particles,para_env=para_env,&
             n_el=multiindex,&
             listbased_distribution=.FALSE.,error=error)

        multiindex=0
        !$omp parallel do private(imol,j,atombase,ikind)
        DO imol = 1, SIZE(local_molecules)
           atombase = local_molecules(imol) % atombase - 1
           DO j = 1, local_molecules(imol) % natoms_mol
              CALL get_atomic_kind(particles( atombase + j )%atomic_kind,&
                   kind_number=ikind)
              multiindex(ikind)=multiindex(ikind)+1
              local_particles%list(ikind)%array(multiindex(ikind))=atombase+j
           END DO
        END DO
        
     END IF

     DEALLOCATE(multiindex,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

  END IF
END SUBROUTINE mol_local_particles_create
!***************************************************************************

END MODULE initialize_molecule_types

!******************************************************************************
