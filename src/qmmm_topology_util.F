!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author teo
! *****************************************************************************
MODULE qmmm_topology_util
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_topology_util'

  PUBLIC :: qmmm_coordinate_control,&
            qmmm_connectivity_control

CONTAINS

! *****************************************************************************
!> \brief Modifies the atom_info%id_atmname
!> \param topology ...
!> \param qmmm_env ...
!> \param subsys_section ...
!> \param error ...
!> \par History
!>      11.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_coordinate_control(topology, qmmm_env, subsys_section, error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_coordinate_control', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: prefix_lnk
    INTEGER                                  :: handle, iatm, iw
    LOGICAL                                  :: failure, qmmm_index_in_range
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)
    IF (iw>0) WRITE(iw,*) "  Entering qmmm_coordinate_control"
    !
    ! setting ilast and ifirst for QM molecule
    !
    CPPostcondition(SIZE(qmmm_env%qm_atom_index) /= 0,cp_failure_level,routineP,error,failure)
    qmmm_index_in_range=(MAXVAL(qmmm_env%qm_atom_index)<=SIZE(topology%atom_info%id_atmname)) &
                        .AND.(MINVAL(qmmm_env%qm_atom_index)>0)
    CPPostcondition(qmmm_index_in_range,cp_failure_level,routineP,error,failure)
    DO iatm = 1, SIZE(qmmm_env%qm_atom_index)
       topology%atom_info%id_atmname(qmmm_env%qm_atom_index(iatm)) = str2id(s2s("_QM_"//&
            TRIM(id2str(topology%atom_info%id_atmname(qmmm_env%qm_atom_index(iatm))))))
       topology%atom_info%id_resname(qmmm_env%qm_atom_index(iatm)) = str2id(s2s("_QM_"//&
            TRIM(id2str(topology%atom_info%id_resname(qmmm_env%qm_atom_index(iatm))))))
    END DO
    !
    ! Modify type for MM link atoms
    !
    IF (ASSOCIATED(qmmm_env%mm_link_atoms)) THEN
       DO iatm = 1, SIZE(qmmm_env%mm_link_atoms)
          prefix_lnk = "_LNK000"
          WRITE(prefix_lnk(5:),'(I20)')iatm
          CALL compress(prefix_lnk,.TRUE.)
          topology%atom_info%id_atmname(qmmm_env%mm_link_atoms(iatm)) = str2id(s2s(TRIM(prefix_lnk)//&
               TRIM(id2str(topology%atom_info%id_atmname(qmmm_env%mm_link_atoms(iatm))))))
          topology%atom_info%id_resname(qmmm_env%mm_link_atoms(iatm)) = str2id(s2s(TRIM(prefix_lnk)//&
               TRIM(id2str(topology%atom_info%id_resname(qmmm_env%mm_link_atoms(iatm))))))
       END DO
    END IF
    !
    IF (iw>0) WRITE(iw,*) "  Exiting  qmmm_coordinate_control"
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    CALL timestop(handle)
  END SUBROUTINE qmmm_coordinate_control

! *****************************************************************************
!> \brief Set up the connectivity for QM/MM calculations
!> \param molecule_set ...
!> \param qmmm_env ...
!> \param subsys_section ...
!> \param error ...
!> \par History
!>      12.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_connectivity_control(molecule_set,&
       qmmm_env, subsys_section, error)

    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_connectivity_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: first_atom, handle, i, &
                                                imolecule, iw, last_atom, &
                                                natom, output_unit, &
                                                qm_mol_num, stat
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, &
                                                qm_molecule_index
    LOGICAL                                  :: detected_link, failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

    NULLIFY(qm_atom_index, qm_molecule_index, molecule, molecule_kind)
    failure = .FALSE.
    detected_link = .FALSE.
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/UTIL_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)
    IF (.NOT.failure) THEN
       qm_mol_num    = 0
       qm_atom_index => qmmm_env%qm_atom_index
       DO imolecule=1,SIZE(molecule_set)
          IF (iw>0) WRITE(iw,*)"Entering molecule number ::",imolecule
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
          IF (iw>0) WRITE(iw,*)"Entering molecule number ::",imolecule
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
                IF (iw>0) WRITE(iw,fmt='(A)',ADVANCE="NO") " QM/MM link detected..."
                IF (.NOT.qmmm_env%qmmm_link) THEN
                   IF (iw>0) WRITE(iw,fmt='(A)')" Missing LINK section in input file!!"
                   WRITE(output_unit,'(T2,"QMMM_CONNECTIVITY|",A)')&
                        " ERROR in the QM/MM connectivity. A QM/MM LINK was detected but",&
                        " no LINK section was provided in the Input file!",&
                        " This very probably can be identified as an error in the specified QM",&
                        " indexes or in a missing LINK section. Check your structure!"
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END IF
             END IF
             qm_molecule_index(qm_mol_num) = imolecule
          END IF
       END DO
       IF (ASSOCIATED(qmmm_env%qm_molecule_index)) DEALLOCATE(qmmm_env%qm_molecule_index)
       NULLIFY(qmmm_env%qm_molecule_index)
       qmmm_env%qm_molecule_index => qm_molecule_index
       IF(iw>0) WRITE(iw,*) "    QM molecule index ::",qm_molecule_index
    END IF
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
                                     "PRINT%TOPOLOGY_INFO/UTIL_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE qmmm_connectivity_control

END MODULE qmmm_topology_util
