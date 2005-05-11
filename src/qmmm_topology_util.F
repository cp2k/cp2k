!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qmmm_topology_util [1.0] *
!!
!!   NAME
!!     qmmm_topology_util
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************


MODULE qmmm_topology_util
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE molecule_types_new,              ONLY: molecule_type,&
                                             get_molecule
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE topology_types,                  ONLY: topology_parameters_type
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_utilities,                ONLY: compress
  IMPLICIT NONE
  PRIVATE
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "qmmm_topology_util"

  PUBLIC :: qmmm_coordinate_control,&
            qmmm_connectivity_control

CONTAINS
  !!*************************************************************************
  !!****** qmmm_topology_util/qmmm_coordinate_control [1.0] *
  !!
  !!   NAME
  !!     qmmm_coordinate_control
  !!
  !!   FUNCTION
  !!     Modifies the atom_info%label_atmname
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     11.2004 created [tlaino]
  !!
  !***************************************************************************
  SUBROUTINE qmmm_coordinate_control(topology, globenv, qmmm_env, error)
    IMPLICIT NONE
    
    TYPE(topology_parameters_type), &
         INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
         INTENT(IN)                             :: globenv
    TYPE(qmmm_env_mm_type), pointer             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_coordinate_control', &
         routineP = moduleN//':'//routineN
    CHARACTER (LEN=default_string_length)    :: prefix_lnk
    INTEGER                                  :: handle, iw, iatm
    INTEGER                                  :: ifirst, ilast, myind
    
    failure = .FALSE.
    iw = globenv%scr
    CALL write_checkpoint_information("entering "//routineN,globenv)
    CALL timeset(routineN,'I','',handle)
    IF((globenv%ionode).AND.(globenv%print%level>0))THEN
       WRITE(iw,*) "  Entering qmmm_coordinate_control"
    END IF
    !
    ! setting ilast and ifirst for QM molecule
    !
    CPPostcondition(SIZE(qmmm_env%qm_atom_index) /= 0,cp_failure_level,routineP,error,failure)
    DO iatm = 1, SIZE(qmmm_env%qm_atom_index)
       topology%atom_info%label_atmname(qmmm_env%qm_atom_index(iatm)) = "_QM_"//&
            TRIM(topology%atom_info%label_atmname(qmmm_env%qm_atom_index(iatm)))
       topology%atom_info%label_resname(qmmm_env%qm_atom_index(iatm)) = "_QM_"//&
            TRIM(topology%atom_info%label_resname(qmmm_env%qm_atom_index(iatm)))
    END DO
    !
    ! Modify type for MM link atoms 
    !
    IF (ASSOCIATED(qmmm_env%mm_link_atoms)) THEN
       DO iatm = 1, SIZE(qmmm_env%mm_link_atoms)
          prefix_lnk = "_LNK000"
          WRITE(prefix_lnk(5:),'(I20)')iatm
          CALL compress(prefix_lnk,.TRUE.)
          topology%atom_info%label_atmname(qmmm_env%mm_link_atoms(iatm)) = TRIM(prefix_lnk)//&
               TRIM(topology%atom_info%label_atmname(qmmm_env%mm_link_atoms(iatm)))
          topology%atom_info%label_resname(qmmm_env%mm_link_atoms(iatm)) = TRIM(prefix_lnk)//&
               TRIM(topology%atom_info%label_resname(qmmm_env%mm_link_atoms(iatm))) 
       END DO
    END IF
    !
    IF((globenv%ionode).AND.(globenv%print%level>0))THEN
       WRITE(iw,*) "  Exiting  qmmm_coordinate_control"
    END IF
    
    CALL timestop(0.0_dp,handle)
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END SUBROUTINE qmmm_coordinate_control
  

  !!*************************************************************************
  !!****** qmmm_topology_util/qmmm_connectivity_control [1.0] *
  !!
  !!   NAME
  !!     qmmm_connectivity_control
  !!
  !!   FUNCTION
  !!     Set up the connectivity for QM/MM calculations
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !***************************************************************************
  SUBROUTINE qmmm_connectivity_control(topology, globenv, molecule_set,&
       qmmm_env, error)
    IMPLICIT NONE
    
    TYPE(topology_parameters_type), &
         INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
         INTENT(IN)                             :: globenv
    TYPE(molecule_type), DIMENSION(:), &
         POINTER                                :: molecule_set
    TYPE(qmmm_env_mm_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    
    LOGICAL :: failure, detected_link
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_connectivity_control', &
         routineP = moduleN//':'//routineN
    INTEGER                                  :: stat, iw, imolecule, handle, i
    INTEGER :: first_atom,last_atom, qm_mol_num, natom
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, qm_molecule_index

    NULLIFY(qm_atom_index, qm_molecule_index, molecule, molecule_kind)
    failure = .FALSE.
    detected_link = .FALSE.
    iw = globenv%scr
    CALL write_checkpoint_information("entering "//routineN,globenv)
    CALL timeset(routineN,'I','',handle)
    IF (.NOT.failure) THEN
       IF((globenv%ionode).AND.(globenv%print%level>0))THEN
          WRITE(iw,*) "  Entering qmmm_connectivity_control"
       END IF
       qm_mol_num    = 0
       qm_atom_index => qmmm_env%qm_atom_index
       DO imolecule=1,SIZE(molecule_set)
          IF ((globenv%ionode).AND.(globenv%print%level>1))&
               WRITE(iw,*)"Entering molecule number ::",imolecule
          molecule => molecule_set(imolecule)
          CALL get_molecule(molecule,molecule_kind=molecule_kind,&
               first_atom=first_atom,last_atom=last_atom)
          IF (ANY(qm_atom_index >= first_atom .AND. qm_atom_index <= last_atom))&
               qm_mol_num =qm_mol_num+1
       END DO
       !
       ALLOCATE(qm_molecule_index(qm_mol_num), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       qm_mol_num    = 0
       DO imolecule=1,SIZE(molecule_set)
          IF ((globenv%ionode).AND.(globenv%print%level>1))&
               WRITE(iw,*)"Entering molecule number ::",imolecule
          molecule => molecule_set(imolecule)
          CALL get_molecule(molecule,molecule_kind=molecule_kind,&
               first_atom=first_atom,last_atom=last_atom)
          natom = last_atom - first_atom + 1
          IF (ANY(qm_atom_index >= first_atom .AND. qm_atom_index <= last_atom)) THEN
             qm_mol_num  = qm_mol_num    + 1
             !
             ! Check if all atoms of the molecule are QM or if a QM/MM link scheme
             ! need to be used...
             !
             detected_link = .FALSE.
             DO i = first_atom, last_atom
                IF (.NOT.ANY(qm_atom_index == i)) detected_link = .TRUE.
             END DO
             IF (detected_link) THEN
                IF (globenv%ionode) THEN
                    WRITE(iw,fmt='(A)',ADVANCE="NO") " QM/MM link detected..."
                ENDIF
                IF (.NOT.qmmm_env%qmmm_link) THEN 
                   WRITE(iw,fmt='(A)')" Missing LINK section in input file!!"
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END IF
             END IF
             qm_molecule_index(qm_mol_num) = imolecule
          END IF
       END DO
       IF (ASSOCIATED(qmmm_env%qm_molecule_index)) DEALLOCATE(qmmm_env%qm_molecule_index)
       NULLIFY(qmmm_env%qm_molecule_index)
       qmmm_env%qm_molecule_index => qm_molecule_index
       IF((globenv%ionode).AND.(globenv%print%level>0))THEN
          WRITE(iw,*) "    QM molecule index ::",qm_molecule_index
          WRITE(iw,*) "  Exiting qmmm_connectivity_control"
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
    CALL write_checkpoint_information("leaving "//routineN,globenv)
  END SUBROUTINE qmmm_connectivity_control


END MODULE qmmm_topology_util
