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

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atoms_input,                     ONLY: read_atoms_input
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE force_fields,                    ONLY: read_force_field_amber,&
                                             read_force_field_charmm,&
                                             read_force_field_section
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE pair_potential_types,            ONLY: pair_potential_type
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: read_cell
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_input,                  ONLY: read_topology_section
  USE topology_pdb,                    ONLY: read_coordinate_pdb,&
                                             write_coordinate_pdb
  USE topology_psf,                    ONLY: read_topology_psf,&
                                             write_topology_psf
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: topology_connectivity_pack,&
                                             topology_constraint_pack,&
                                             topology_coordinate_pack,&
                                             topology_force_field_pack,&
                                             topology_generate_bend,&
                                             topology_generate_bond,&
                                             topology_generate_dihe,&
                                             topology_generate_molecule,&
                                             topology_reorder_atoms

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology"

  PRIVATE
  PUBLIC :: topology_control, connectivity_control, force_field_control

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
!!
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE topology_control (atomic_kind_set,particle_set,&
                             molecule_kind_set,molecule_set,&
                             globenv,potparm)


    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(pair_potential_type), OPTIONAL, &
      POINTER                                :: potparm( :, : )

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    TYPE(topology_parameters_type)           :: topology

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  iw = globenv%scr

  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)

  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "Entering topology_control"
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. Initialize the topology structure type
  !-----------------------------------------------------------------------------
  CALL init_topology(topology)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Get the cell info
  !-----------------------------------------------------------------------------
  CALL read_cell(topology%cell,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Read in the topology section in the input file if any
  !-----------------------------------------------------------------------------
  CALL read_topology_section(topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Read in the coordinates
  !-----------------------------------------------------------------------------
  CALL coordinate_control(topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Read in or generate the molecular connectivity
  !-----------------------------------------------------------------------------
  CALL connectivity_control(topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Read in force field informations
  !-----------------------------------------------------------------------------
  CALL force_field_control(topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Pack everything into the molecular types
  !-----------------------------------------------------------------------------
  CALL topology_connectivity_pack(molecule_kind_set,molecule_set,&
                                  topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Pack everything into the atomic types
  !-----------------------------------------------------------------------------
  CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
                                molecule_kind_set,molecule_set,&
                                topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 9. Pack everything into the atomic types
  !-----------------------------------------------------------------------------
  IF((topology%ff).AND.(PRESENT(potparm))) THEN
      CALL topology_force_field_pack(particle_set,atomic_kind_set,&
                                     molecule_kind_set,molecule_set,&
                                     topology,globenv,potparm)
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 10. Deal with the constraint stuff if requested
  !-----------------------------------------------------------------------------
  IF(topology%constraint) THEN
    CALL topology_constraint_pack(molecule_kind_set,molecule_set,&
                                  topology,globenv)
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 11. Dump the topology informations
  !-----------------------------------------------------------------------------
  IF(topology%dump_topology) THEN
    CALL write_topology_psf(topology,globenv)
    CALL write_coordinate_pdb(topology,globenv)
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 12. Cleanup the topology structure type
  !-----------------------------------------------------------------------------
  CALL deallocate_topology(topology)



  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "Exiting  topology_control"
  END IF

  CALL timestop(0.0_dbl,handle)
  CALL write_checkpoint_information("leaving "//routineN,globenv)
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


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'connectivity_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: found

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  iw = globenv%scr
  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)

  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Entering connectivity_control"
  END IF


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
    CALL topology_generate_bond(topology,globenv)
    CALL topology_generate_bend(topology,globenv)
    CALL topology_generate_dihe(topology,globenv)
    CALL topology_generate_molecule(topology,globenv)
    IF(topology%reorder_atom) THEN
      CALL topology_reorder_atoms(topology,globenv)
      CALL topology_generate_bond(topology,globenv)
      CALL topology_generate_bend(topology,globenv)
      CALL topology_generate_dihe(topology,globenv)
      CALL topology_generate_molecule(topology,globenv)
    END IF
  END SELECT


  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Exiting  connectivity_control"
  END IF

  CALL timestop(0.0_dbl,handle)
  CALL write_checkpoint_information("leaving "//routineN,globenv)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
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
SUBROUTINE coordinate_control(topology,globenv)


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'coordinate_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: found

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  iw = globenv%scr
  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)
  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Entering coordinate_control"
  END IF


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


  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Exiting  coordinate_control"
  END IF

  CALL timestop(0.0_dbl,handle)
  CALL write_checkpoint_information("leaving "//routineN,globenv)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE coordinate_control




!******************************************************************************
!!****** topology/force_field_control [1.0] *
!!
!!   NAME
!!     force_field_control
!!
!!   FUNCTION
!!     1. If reading in from external file, make sure its there first
!!     2. Read in the force_field from the corresponding locations
!!
!!   AUTHOR
!!     IKUO 08.11.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE force_field_control(topology,globenv)


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'force_field_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: found

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  iw = globenv%scr
  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)
  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Entering force_field_control"
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 1. If reading in from external file, make sure its there first
  !-----------------------------------------------------------------------------
  IF((topology%ff).AND.(topology%ff_type/="ON")) THEN
    !------------------------------------------------------
    !Make sure the file really do exist before proceding
    IF(globenv%ionode) THEN
      INQUIRE(FILE=topology%ff_file_name,EXIST=found)
      IF(found) THEN
      ELSE
        CALL stop_program ("force_field_control","force field file missing")
      END IF
    END IF
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Read in the force field from the corresponding locations
  !-----------------------------------------------------------------------------
  IF((topology%ff).AND.(topology%ff_type/="ON")) THEN
    SELECT CASE (topology%ff_type)
    CASE ("CHM")
      CALL stop_program ("force_field_control",&
                         "force field type not implemented")
      CALL read_force_field_charmm(topology,globenv)
    CASE ("AMB")
      CALL stop_program ("force_field_control",&
                         "force field type not implemented")
      CALL read_force_field_amber(topology,globenv)
    CASE DEFAULT
      CALL stop_program ("force_field_control",&
                         "force field type not implemented")
    END SELECT
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Parse the input for force field section and append any info
  !-----------------------------------------------------------------------------
  CALL read_force_field_section(topology,globenv)



  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Exiting  force_field_control"
  END IF

  CALL timestop(0.0_dbl,handle)
  CALL write_checkpoint_information("leaving "//routineN,globenv)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE force_field_control

END MODULE topology

!******************************************************************************
