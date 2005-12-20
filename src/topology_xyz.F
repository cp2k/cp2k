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
  SUBROUTINE read_coordinate_xyz (topology,globenv)
    IMPLICIT NONE
    ! Arguments
    TYPE(topology_parameters_type)              :: topology
    TYPE(global_environment_type), &
         pointer, OPTIONAL                   :: globenv
    ! Local Variables
    TYPE (atom_info_type),POINTER               :: atom_info
    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_xyz', &
         routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                          :: nblock_atom = 1000

    CHARACTER(LEN=default_string_length)        :: line, record, root_mol_name
    INTEGER                                     :: handle, istat, iw, natom, i, j, k, m, inum_mol
    INTEGER                                     :: file_unit, frame
    LOGICAL                                     :: PM1, PM2, PM3, PM4

    PM1 = .FALSE.
    PM2 = .FALSE.
    PM3 = .FALSE.
    PM4 = .FALSE.
    IF(PRESENT(globenv)) THEN
       CALL write_checkpoint_information("entering "//routineN,globenv)
       iw = globenv%scr
       IF((globenv%ionode).AND.(globenv%print%level>0)) PM1= .TRUE.
       IF((globenv%ionode).AND.(globenv%print%level>1)) PM2= .TRUE.
       IF((globenv%ionode).AND.(globenv%print%level>2)) PM3= .TRUE.
       IF((globenv%ionode).AND.(globenv%print%level>3)) PM4= .TRUE.
    END IF
    
    IF(PM1) WRITE(iw,*) "  Entering read_coordinate_xyz"

    CALL write_checkpoint_information("entering "//routineN,globenv)
    CALL timeset(routineN,'I','',handle)
    
    atom_info => topology%atom_info

    IF(PM1) WRITE(iw,*) "    Reading in XYZ file ",TRIM(topology%coord_file_name)
    CALL open_file(file_name=topology%coord_file_name,&
                   file_action="READ",&
                   unit_number=file_unit)

    natom         = 0
    inum_mol      = 1
    WRITE(root_mol_name,'(A3,I7)')"MOL",inum_mol
    DO i =  1, 10
       if (root_mol_name(i:i) == " ") root_mol_name(i:i) = "0" 
    END DO
    frame = 0
    Frames: DO 
       ! Atom numbers
       READ(file_unit,*,IOSTAT=istat)natom
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
       line = ''; READ(file_unit,'(A80)',IOSTAT=istat) line 
       DO j = 1, natom
          line = ''
          READ(file_unit,'(A80)',IOSTAT=istat) line 
          IF(istat /= 0) EXIT Frames
          ! Atom coordinates
          atom_info%label_molname(j) = TRIM(root_mol_name)
          k = 1
          DO WHILE (line(k:k) == " ")
             k = k + 1
          END DO
          m = k
          DO WHILE (line(m:m) /= " ")
             m = m + 1
          END DO
          atom_info%label_atmname(j)   = line(k:m-1)
          CALL uppercase(atom_info%label_atmname(j))
          ! For default, set atom name to residue name to molecule name
          atom_info%label_molname(j)   = TRIM(atom_info%label_molname(j))//"_"//&
                                         TRIM(atom_info%label_atmname(j))
          atom_info%label_resname(j)   = atom_info%label_molname(j)
          atom_info%label_resid(j)     = 1
          READ(line(m:),*,IOSTAT=istat) atom_info%r(1:3,j)
          atom_info%element(j)         = atom_info%label_atmname(j)
          atom_info%atm_mass(j)        = 0.0_dp
          atom_info%atm_charge(j)      = -10000.0_dp
          IF(PM2) WRITE(iw,*) j,line
          IF(PM2) WRITE(iw,'(5X,A4,3F8.3,A4,I4)') &
               TRIM(atom_info%label_atmname(j)),&
               atom_info%r(1,j),&
               atom_info%r(2,j),&
               atom_info%r(3,j),&
               TRIM(atom_info%label_resname(j)),&
               atom_info%label_resid(j)
          CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(1,j))
          CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(2,j))
          CALL convert_to_cp2k_units("ANGSTROM",length=atom_info%r(3,j))
       END DO
    END DO Frames
    CALL close_file(unit_number=file_unit)

    IF(.NOT.topology%para_res) atom_info%label_resid(:) = 1
    topology%natoms = natom

    IF(PM1) WRITE(iw,*) "  Exiting read_coordinate_xyz"

    CALL timestop(0.0_dp,handle)
    IF(PRESENT(globenv)) THEN
       CALL write_checkpoint_information("exiting "//routineN,globenv)
    END IF
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
  END SUBROUTINE read_coordinate_xyz

END MODULE topology_xyz
