!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/topology_control [1.0] *
!!
!!   NAME
!!     topology_control
!!
!!   FUNCTION
!!     Control for reading in different topologies and coordinates
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE topology

  USE atomic_kind_types,                ONLY: atomic_kind_type
  USE atoms_input,                      ONLY: read_atoms_input
  USE global_types,                     ONLY: global_environment_type
  USE kinds,                            ONLY: dbl
  USE molecule_types_new,               ONLY: allocate_molecule_set,&
                                              molecule_type,&
                                              set_molecule_set
  USE molecule_kind_types,              ONLY: allocate_molecule_kind_set,&
                                              atom_type,&
                                              get_molecule_kind,&
                                              molecule_kind_type,&
                                              read_molecule_kind_set,&
                                              set_molecule_kind
  USE particle_types,                   ONLY: particle_type,&
                                              read_particle_set
  USE simulation_cell,                  ONLY: cell_type,get_cell,&
                                              read_cell,write_cell
  USE termination,                      ONLY: stop_program,&
                                              stop_memory
  USE topology_input,                   ONLY: read_topology_section
  USE topology_gromos,                  ONLY: read_topology_gromos
  USE topology_pdb,                     ONLY: read_coordinate_pdb
  USE topology_psf,                     ONLY: read_topology_psf
  USE topology_types,                   ONLY: topology_parameters_type
  USE topology_util,                    ONLY: topology_coordinate_pack,&
                                              topology_generate_bond,&
                                              topology_generate_bend,&
                                              topology_generate_dihe,&
                                              topology_generate_molecule,&
                                              topology_connectivity_pack

  PRIVATE
  PUBLIC :: topology_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology/topology_control [1.0] *
!!
!!   NAME
!!     topology_control
!!
!!   FUNCTION
!!     1. Read in the topology section in the input file if any
!!     2. Read in the coordinates
!!     3. Read in or generate the molecular connectivity
!!     4. Pack everything into the molecular types
!!     5. Pack everything into the atomic types
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE topology_control (atomic_kind_set,particle_set,&
                             molecule_kind_set,molecule_set,&
                             globenv)

  IMPLICIT NONE

! Arguments
  TYPE ( atomic_kind_type ), DIMENSION(:), POINTER  :: atomic_kind_set
  TYPE ( particle_type ), DIMENSION(:), POINTER     :: particle_set
  TYPE ( molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
  TYPE ( molecule_type), DIMENSION(:), POINTER      :: molecule_set
  TYPE ( global_environment_type ), INTENT ( IN )   :: globenv


! Locals
  TYPE ( topology_parameters_type ) :: topology
  INTEGER :: iw

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  iw = globenv%scr

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Get the cell info
  !-----------------------------------------------------------------------------
  NULLIFY(topology%cell)
  CALL read_cell(topology%cell,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Read in the topology section in the input file if any
  !-----------------------------------------------------------------------------
  CALL read_topology_section(topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Read in the coordinates
  !-----------------------------------------------------------------------------
  CALL coordinate_control(particle_set,atomic_kind_set,&
                          molecule_kind_set,molecule_set,&
                          topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Read in or generate the molecular connectivity
  !-----------------------------------------------------------------------------
  CALL connectivity_control(topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Pack everything into the molecular types
  !-----------------------------------------------------------------------------
  CALL topology_connectivity_pack(molecule_kind_set,molecule_set,&
                                  topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Pack everything into the atomic types
  !-----------------------------------------------------------------------------
  CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
                                molecule_kind_set,molecule_set,&
                                topology,globenv)



  !CALL stop_program ("topology_control","program check")


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE topology_control


!******************************************************************************
!!****** topology/connectivity_control [1.0] *
!!
!!   NAME
!!     connectivity_control
!!
!!   FUNCTION
!!     1. If reading in from external file, make sure its there first
!!     2. Generate the connectivity if no information to be read in
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE connectivity_control(topology,globenv)

  IMPLICIT NONE

! Arguments
  TYPE ( topology_parameters_type ), INTENT (INOUT)   :: topology
  TYPE ( global_environment_type ), INTENT ( IN )  :: globenv

! Locals
  LOGICAL found
  INTEGER :: iw, ierr
  INTEGER :: natom,nbond

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  iw = globenv%scr


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. If reading in from external file, make sure its there first
  !-----------------------------------------------------------------------------
  IF(topology%connectivity) THEN
    !------------------------------------------------------
    !Make sure the file really do exist before proceding
    IF(globenv%ionode) THEN
      INQUIRE(FILE=topology%conn_file_name,EXIST=found)
      IF(found) THEN
        WRITE(iw,*) 'Using connectivity file ',topology%conn_file_name
      ELSE
        CALL stop_program ("connectivity_control","connectivity file missing")
      END IF
    END IF
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Read in or generate the connectivity information any way we can 
  !-----------------------------------------------------------------------------
  SELECT CASE (topology%conn_type)
  CASE ("PSF")
    CALL read_topology_psf(topology,globenv)
  CASE ("GROMOS")
    !WRITE(iw,*) 'Using topology file ',topology%conn_file_name
    !CALL read_topology_gromos(molecule_set,topology,globenv)
    CALL stop_program ("connectivity_control",&
                       "gromos not implemented yet...")
  CASE DEFAULT
    !Will use topology generate to set up the connectivity information
    CALL topology_generate_bond(topology)
    CALL topology_generate_bend(topology)
    CALL topology_generate_dihe(topology)
    CALL topology_generate_molecule(topology)
  END SELECT


END SUBROUTINE connectivity_control




!******************************************************************************
!!****** topology/coordinate_control [1.0] *
!!
!!   NAME
!!     coordinate_control
!!
!!   FUNCTION
!!     1. If reading in from external file, make sure its there first
!!     2. Read in the coordinates from the corresponding locations
!!
!!   AUTHOR
!!     IKUO 08.11.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE coordinate_control(particle_set,atomic_kind_set,&
                              molecule_kind_set,molecule_set,&
                              topology,globenv)

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ),POINTER :: molecule_set
  TYPE ( molecule_kind_type ), DIMENSION ( : ),POINTER :: molecule_kind_set
  TYPE ( topology_parameters_type ), INTENT (INOUT) :: topology
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

! Locals
  LOGICAL found
  CHARACTER ( LEN = 20 ) :: unit_type
  CHARACTER ( LEN = 4 ) :: fext
  INTEGER :: iw, ierr
  INTEGER :: i, natom
  INTEGER :: input_select

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  iw = globenv%scr

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. If reading in from external file, make sure its there first
  !-----------------------------------------------------------------------------
  IF(topology%coordinate) THEN
    !------------------------------------------------------
    !Make sure the file really do exist before proceding
    IF(globenv%ionode) THEN
      INQUIRE(FILE=topology%coord_file_name,EXIST=found)
      IF(found) THEN
      ELSE
        CALL stop_program ("coordinate_control","coordinate file missing")
      END IF
    END IF
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Read in the coordinates from the corresponding locations
  !-----------------------------------------------------------------------------
  SELECT CASE (topology%coord_type)
  CASE ("PDB")
    CALL read_coordinate_pdb (topology,globenv)
  CASE ("CHM")
    !CALL read_coordinate_chm (r,label_resname,label_atmname,&
    !                          topology,globenv)
    CALL stop_program ("coordinate_control",&
                       "CHM coordinate not implemented yet...")
  CASE ("XYZ")
    !CALL read_coordinate_xyz (r,label_resname,label_atmname,&
    !                          topology,globenv)
    CALL stop_program ("coordinate_control",&
                       "XYZ coordinate not implemented yet...")
  CASE DEFAULT
    !CALL read_particle_set (particle_set,atomic_kind_set,molecule_set,&
    !                        molecule_kind_set,topology%cell,globenv)
    CALL read_atoms_input ( topology, globenv )
  END SELECT


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE coordinate_control

END MODULE topology

!******************************************************************************
