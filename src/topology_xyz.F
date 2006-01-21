!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/topology_xyz [1.0] *
!!
!!   NAME
!!     topology_xyz
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE topology_xyz

  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
       default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE physcon,                         ONLY: bohr
  USE cp_files,                        ONLY: close_file,&
       open_file
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
       timestop
  USE topology_types,                  ONLY: topology_parameters_type,&
       atom_info_type
  USE cp_parser,                       ONLY: parser_release,&
                                             parser_get_object,&
                                             parser_search_string,&
                                             parser_create,&
                                             parser_test_next_token,&
                                             cp_parser_type,&
                                             parser_get_next_line
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_finished_output  
  USE input_section_types,             ONLY: section_vals_type

  

#include "cp_common_uses.h"  


  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_xyz"

  PRIVATE
  PUBLIC :: read_coordinate_xyz

CONTAINS

  !******************************************************************************
  !!****** topology_xyz/read_coordinate_xyz [1.0] *
  !!
  !!   NAME
  !!     read_coordinate_xyz
  !!
  !!   FUNCTION
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!
  !!   OPTIONS
  !!
  !!   INPUTS
  !!*** *************************************************************************
  ! (1) The first line contains the number of atoms in the molecule.
  ! (2) The second line contains the name or title of the molecule. It can be blank if desired.
  ! (3) The rest of the file is composed of one line for each atom in the molecule. 
  !     Each atom is designated by the atomic number and then the x, y, and z coordinate value.
  !     Separate each value by at least 1 space; there is no strict format for these lines.
  SUBROUTINE read_coordinate_xyz (topology,globenv,subsys_section,error)
    IMPLICIT NONE
    ! Arguments
    TYPE(topology_parameters_type)              :: topology
    TYPE(global_environment_type), &
         pointer, OPTIONAL                   :: globenv
    TYPE(section_vals_type), POINTER :: subsys_section

    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    TYPE (atom_info_type),POINTER               :: atom_info
    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_xyz', &
         routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                          :: nblock_atom = 1000

    CHARACTER(LEN=default_string_length)        :: root_mol_name
    INTEGER                                     :: handle, iw, natom, i, j, inum_mol
    INTEGER                                     :: frame
    LOGICAL                                     ::  my_end
    TYPE(cp_parser_type), POINTER               :: parser    
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(parser, logger)
    logger => cp_error_get_logger(error)
    CALL write_checkpoint_information("Entering "//routineN,globenv,error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/XYZ_INFO",&
         extension=".subsysLog",error=error)    
    CALL timeset(routineN,'I','',handle)
    
    atom_info => topology%atom_info

    IF(iw>0) WRITE(iw,*) "    Reading in XYZ file ",TRIM(topology%coord_file_name)
    CALL parser_create(parser,topology%coord_file_name,&
         para_env=globenv%para_env,error=error)      

    natom         = 0
    inum_mol      = 1
    WRITE(root_mol_name,'(A3,I7)')"MOL",inum_mol
    DO i =  1, 10
       IF (root_mol_name(i:i) == " ") root_mol_name(i:i) = "0" 
    END DO
    frame = 0
    CALL parser_get_next_line(parser,1,error=error)
    Frames: DO 
       ! Atom numbers
       CALL parser_get_object  (parser,natom)
       topology%natoms = natom
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
          atom_info%label_molname(j) = TRIM(root_mol_name)
          CALL parser_get_object  (parser,atom_info%label_atmname(j))
          CALL uppercase(atom_info%label_atmname(j))
          ! For default, set atom name to residue name to molecule name
          atom_info%label_molname(j)   = TRIM(atom_info%label_molname(j))//"_"//&
                                         TRIM(atom_info%label_atmname(j))
          atom_info%label_resname(j)   = atom_info%label_molname(j)
          atom_info%label_resid(j)     = 1
          CALL parser_get_object  (parser,atom_info%r(1,j))
          CALL parser_get_object  (parser,atom_info%r(2,j))
          CALL parser_get_object  (parser,atom_info%r(3,j))
          atom_info%element(j)         = atom_info%label_atmname(j)
          atom_info%atm_mass(j)        = 0.0_dp
          atom_info%atm_charge(j)      = -10000.0_dp
          IF(iw>0) WRITE(iw,'(5X,A4,3F8.3,A4,I4)') &
               TRIM(atom_info%label_atmname(j)),&
               atom_info%r(1,j),&
               atom_info%r(2,j),&
               atom_info%r(3,j),&
               TRIM(atom_info%label_resname(j)),&
               atom_info%label_resid(j)
          CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(1,j))
          CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(2,j))
          CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(3,j))
          CALL parser_get_next_line(parser,1,at_end=my_end,error=error)
          IF (my_end) EXIT Frames
       END DO
    END DO Frames
    CALL parser_release(parser,error=error)

    IF(.NOT.topology%para_res) atom_info%label_resid(:) = 1
    topology%natoms = natom
    CALL timestop(0.0_dp,handle)
    CALL write_checkpoint_information("exiting "//routineN,globenv,error)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO/XYZ_INFO",error=error)

  END SUBROUTINE read_coordinate_xyz

END MODULE topology_xyz
