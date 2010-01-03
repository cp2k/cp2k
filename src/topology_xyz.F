!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE topology_xyz
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_xyz'

  PRIVATE
  PUBLIC :: read_coordinate_xyz

CONTAINS

! *****************************************************************************
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE read_coordinate_xyz (topology,para_env,subsys_section,error)
    TYPE(topology_parameters_type)           :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_xyz', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 1000

    CHARACTER(LEN=default_string_length)     :: my_default_index
    INTEGER                                  :: frame, handle, iw, j, natom
    LOGICAL                                  :: failure, my_end
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    NULLIFY(parser, logger)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/XYZ_INFO",&
         extension=".subsysLog",error=error)
    CALL timeset(routineN,handle)

    atom_info => topology%atom_info

    IF(iw>0) WRITE(iw,*) "    Reading in XYZ file ",TRIM(topology%coord_file_name)
    CALL parser_create(parser,topology%coord_file_name, para_env=para_env, &
                       parse_white_lines=.TRUE.,error=error)

    ! Element is assigned on the basis of the atm_name
    topology%aa_element = .TRUE.

    natom = 0
    frame = 0
    CALL parser_get_next_line(parser,1,error=error)
    Frames: DO
       ! Atom numbers
       CALL parser_get_object  (parser,natom,error=error)
       frame = frame + 1
       IF (frame == 1) THEN
          CALL reallocate(atom_info%label_molname, 1,    natom)
          CALL reallocate(atom_info%label_resname, 1,    natom)
          CALL reallocate(atom_info%label_resid,   1,    natom)
          CALL reallocate(atom_info%label_atmname, 1,    natom)
          CALL reallocate(atom_info%r,             1,3,1,natom)
          CALL reallocate(atom_info%atm_mass,      1,    natom)
          CALL reallocate(atom_info%atm_charge,    1,    natom)
          CALL reallocate(atom_info%occup,         1,    natom)
          CALL reallocate(atom_info%beta,          1,    natom)
          CALL reallocate(atom_info%element,       1,    natom)
       ELSE IF(natom > SIZE(atom_info%label_atmname)) THEN
          CALL stop_program (routineP,&
               "Atom number differs in different frames!")
       END IF
       ! Dummy line
       CALL parser_get_next_line(parser,2,error=error)
       DO j = 1, natom
          ! Atom coordinates
          READ(parser%input_line,*)atom_info%label_atmname(j),&
                                   atom_info%r(1,j),&
                                   atom_info%r(2,j),&
                                   atom_info%r(3,j)
          ! For default, set atom name to residue name to molecule name
          WRITE(my_default_index,'(I0)')j   
          atom_info%label_molname(j)   =  TRIM(atom_info%label_atmname(j))//TRIM(my_default_index)
          atom_info%label_resname(j)   =  atom_info%label_molname(j)
          atom_info%label_resid(j)     =  1
          atom_info%element(j)         =  atom_info%label_atmname(j)
          atom_info%atm_mass(j)        =  HUGE(0.0_dp)
          atom_info%atm_charge(j)      = -HUGE(0.0_dp)
          IF(iw>0) WRITE(iw,'(T2,"XYZ_INFO| ",A4,3F8.3,A4,I4)') &
               TRIM(atom_info%label_atmname(j)),&
               atom_info%r(1,j),&
               atom_info%r(2,j),&
               atom_info%r(3,j),&
               TRIM(atom_info%label_resname(j)),&
               atom_info%label_resid(j)
          atom_info%r(1,j) = cp_unit_to_cp2k(atom_info%r(1,j),"angstrom",error=error)
          atom_info%r(2,j) = cp_unit_to_cp2k(atom_info%r(2,j),"angstrom",error=error)
          atom_info%r(3,j) = cp_unit_to_cp2k(atom_info%r(3,j),"angstrom",error=error)
          ! If there's a white line or end of file exit.. otherwise read other available
          ! snapshots
          CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
          my_end = my_end .OR. (LEN_TRIM(parser%input_line)==0)
          IF (my_end) THEN 
             CALL cp_assert(j==natom,cp_fatal_level,cp_assertion_failed,routineP,&
                  "Number of lines in XYZ format not equal to the number of atoms."//&
                  " Error in XYZ format. Very probably the line with title is missing or is empty."//&
                  " Please check the XYZ file and rerun your job!"//&
CPSourceFileRef,&
             only_ionode=.TRUE.)
             EXIT Frames
          END IF
       END DO
    END DO Frames
    CALL parser_release(parser,error=error)

    topology%natoms            = natom
    topology%molname_generated = .TRUE.
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/XYZ_INFO",error=error)

  END SUBROUTINE read_coordinate_xyz

END MODULE topology_xyz
