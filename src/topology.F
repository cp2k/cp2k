!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
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
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE topology
  USE cp_error_handling,               ONLY: cp_error_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atoms_input,                     ONLY: read_atoms_input
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: read_cell
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_input,                  ONLY: read_constraints_section,&
                                             read_topology_section
  USE topology_pdb,                    ONLY: read_coordinate_pdb,&
                                             write_coordinate_pdb
  USE topology_psf,                    ONLY: read_topology_psf,&
                                             write_topology_psf,&
                                             read_topology_mc
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: topology_connectivity_pack,&
                                             topology_constraint_pack,&
                                             topology_coordinate_pack,&
                                             topology_set_atm_mass,&
                                             topology_generate_bend,&
                                             topology_generate_bond,&
                                             topology_generate_dihe,&
                                             topology_generate_impr,&
                                             topology_generate_molecule,&
                                             topology_generate_onfo,&
                                             topology_generate_ub,&
                                             topology_reorder_atoms
  USE qmmm_topology_util,              ONLY: qmmm_coordinate_control,&
                                             qmmm_connectivity_control
  USE qmmm_types,                      only: qmmm_env_mm_type

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology"

  PRIVATE

  PUBLIC :: topology_control, &
            connectivity_control, &
            coordinate_control

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
                             globenv, qmmm, qmmm_env, error)

  TYPE(atomic_kind_type), DIMENSION(:), POINTER    :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER       :: particle_set
  TYPE(molecule_kind_type), DIMENSION(:), POINTER  :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), POINTER       :: molecule_set
  TYPE(global_environment_type), INTENT(IN)        :: globenv
  TYPE(qmmm_env_mm_type), POINTER, OPTIONAL        :: qmmm_env
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  LOGICAL, OPTIONAL :: qmmm

  CHARACTER(len=*), PARAMETER :: routineN = 'topology_control', &
      routineP = moduleN//':'//routineN

  INTEGER                                          :: handle, iw
  TYPE(topology_parameters_type)                   :: topology

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
  CALL read_cell ( topology%cell, topology%cell_ref, globenv )


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Read in the topology section in the input file if any
  !-----------------------------------------------------------------------------
  CALL read_topology_section(topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. Read in the constraints section
  !-----------------------------------------------------------------------------
  CALL read_constraints_section(topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Read in the coordinates
  !-----------------------------------------------------------------------------
  CALL coordinate_control(topology,globenv)

  !-----------------------------------------------------------------------------
  ! Teodoro 11.2004:
  !
  ! The Idea behind the implementation of the QM/MM in Fist is to  define
  ! independently from the topology e/o PDB input file new atom types for
  ! QM  atoms. 
  ! We then modify the force_field_pack routine so to nullify all qm internal
  ! interactions, keep the proper QM/MM internal fields, non-bonded, etc..
  !
  ! Two main routines:
  !    - qmmm_coordinate_control
  !    - qmmm_ff_precond_only_qm (use in force_fields.F)
  !-----------------------------------------------------------------------------
  ! 5.bis. Check if QM/MM and perform the required modifications
  !-----------------------------------------------------------------------------
  IF (PRESENT(qmmm)) THEN
     IF (qmmm) THEN
        IF (.NOT.PRESENT(qmmm_env)) &
             CALL stop_program("Error in qmmm_env in Fist","Illdefined subroutine!")     
        CALL qmmm_coordinate_control(topology, globenv, qmmm_env, error)
     END IF
  END IF
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Read in or generate the molecular connectivity
  !-----------------------------------------------------------------------------
  CALL connectivity_control(topology,globenv,qmmm, qmmm_env, error)

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Pack everything into the molecular types
  !-----------------------------------------------------------------------------
  CALL topology_connectivity_pack(molecule_kind_set,molecule_set,&
                                  topology,globenv)
  !-----------------------------------------------------------------------------
  ! Teo 12.2004:
  ! 
  ! This part takes care of the molecule in which QM atoms were defined.
  ! Preliminary setup for QM/MM link region
  !-----------------------------------------------------------------------------
  ! 7.bis. Set up the QM/MM linkage (if any)
  !-----------------------------------------------------------------------------
  IF (PRESENT(qmmm)) THEN
     IF (qmmm) THEN
        IF (.NOT.PRESENT(qmmm_env)) &
             CALL stop_program("Error in qmmm_env in Fist","Illdefined subroutine!")     
        CALL qmmm_connectivity_control(topology, globenv, molecule_set, &
             qmmm_env, error)
     END IF
  END IF
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Pack everything into the atomic types
  !-----------------------------------------------------------------------------
  CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
                                molecule_kind_set,molecule_set,&
                                topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 9. Deal with the constraint stuff if requested
  !-----------------------------------------------------------------------------
  CALL topology_constraint_pack(molecule_kind_set,molecule_set,&
                                topology,globenv)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 10. Dump the topology informations
  !-----------------------------------------------------------------------------
  IF(topology%dump_topology) THEN
    CALL write_topology_psf(topology,globenv)
    CALL write_coordinate_pdb(topology,globenv)
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 11. Cleanup the topology structure type
  !-----------------------------------------------------------------------------
  CALL deallocate_topology(topology)


  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "Exiting  topology_control"
  END IF

  CALL timestop(0.0_dp,handle)
  CALL write_checkpoint_information("leaving "//routineN,globenv)


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
SUBROUTINE connectivity_control(topology,globenv,qmmm,qmmm_env,error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv
    TYPE(qmmm_env_mm_type), POINTER, OPTIONAL :: qmmm_env
    LOGICAL, INTENT(in), OPTIONAL            :: qmmm
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error    

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
        CALL stop_program ("connectivity_control",&
           "connectivity file missing")
      END IF
    END IF
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Read in or generate the connectivity information any way we can 
  !-----------------------------------------------------------------------------
  SELECT CASE (topology%conn_type)
  CASE ("OHH_MC")
 ! this is a special case, and it's being use for the moment in MC simulations
 ! of water
    CALL read_topology_mc(topology,globenv)
  CASE ("PSF")
    CALL read_topology_psf(topology,globenv)
    !
    IF (PRESENT(qmmm)) THEN
       IF (qmmm) THEN
          IF (.NOT.PRESENT(qmmm_env)) &
               CALL stop_program("Error in qmmm_env in Fist","Illdefined subroutine!")
          CALL qmmm_coordinate_control(topology, globenv, qmmm_env, error)
       END IF
    END IF
    !
    CALL topology_generate_molecule(topology,globenv)
  CASE ("GROMOS")
    !WRITE(iw,*) 'Using topology file ',topology%conn_file_name
    !CALL read_topology_gromos(molecule_set,topology,globenv)
    CALL stop_program ("connectivity_control",&
                       "gromos not implemented yet...")
  CASE ("GENERATE")
    CALL topology_set_atm_mass(topology,globenv)
    CALL topology_generate_bond(topology,globenv)
    CALL topology_generate_bend(topology,globenv)
    CALL topology_generate_ub(topology,globenv)
    CALL topology_generate_dihe(topology,globenv)
    CALL topology_generate_impr(topology,globenv)
    CALL topology_generate_onfo(topology,globenv)
    CALL topology_generate_molecule(topology,globenv)
    IF(topology%reorder_atom) THEN
      CALL topology_reorder_atoms(topology,globenv)
      CALL topology_set_atm_mass(topology,globenv)
      CALL topology_generate_bond(topology,globenv)
      CALL topology_generate_bend(topology,globenv)
      CALL topology_generate_ub(topology,globenv)
      CALL topology_generate_dihe(topology,globenv)
      CALL topology_generate_impr(topology,globenv)
      CALL topology_generate_onfo(topology,globenv)
      CALL topology_generate_molecule(topology,globenv)
    END IF
  CASE ("OFF")
    CALL topology_set_atm_mass(topology,globenv)
    CALL topology_generate_bend(topology,globenv)
    CALL topology_generate_ub(topology,globenv)
    CALL topology_generate_dihe(topology,globenv)
    CALL topology_generate_impr(topology,globenv)
    CALL topology_generate_onfo(topology,globenv)
    CALL topology_generate_molecule(topology,globenv)
  CASE DEFAULT
    CALL stop_program ("connectivity_control",&
                       "no default")
  END SELECT


  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Exiting  connectivity_control"
  END IF

  CALL timestop(0.0_dp,handle)
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

  CALL timestop(0.0_dp,handle)
  CALL write_checkpoint_information("leaving "//routineN,globenv)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE coordinate_control

END MODULE topology

!******************************************************************************
