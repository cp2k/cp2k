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
!!     Control for reading in different topologies
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE topology_control

  USE atomic_kind_types,                ONLY : atomic_kind_type
  USE global_types,                     ONLY : global_environment_type
  USE input_types,                      ONLY : setup_parameters_type
  USE kinds,                            ONLY : dbl
  USE message_passing,                  ONLY : mp_bcast
  USE molecule_types,                   ONLY : molecule_type
  USE particle_types,                   ONLY : particle_type
  USE termination,                      ONLY : stop_program, &
                                               stop_memory
  USE topology_gromos,                  ONLY : read_topology_gromos
  USE topology_pdb,                     ONLY : read_coordinate_pdb
  USE topology_psf,                     ONLY : read_topology_psf

  PRIVATE
  PUBLIC :: read_coordinate_control,read_topology_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology/read_topology_control [1.0] *
!!
!!   NAME
!!     read_topology_control
!!
!!   FUNCTION
!!     none
!!
!!   AUTHOR
!!     IKUO 08.01.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE read_topology_control (mol_setup,setup,globenv)

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), POINTER :: mol_setup
  TYPE ( setup_parameters_type ), INTENT (IN) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  LOGICAL found
  INTEGER :: iw, ierr

  iw = globenv%scr

  IF(setup%topology) THEN
    !Make sure the topology file is there
    IF(globenv%ionode) THEN
      INQUIRE(FILE=setup%topology_file_name,EXIST=found)
      IF(found) THEN
        WRITE(iw,*) 'Using topology file ',setup%topology_file_name
      ELSE
        CALL stop_program ("read_topology_control ","topology file missing")
      END IF
    END IF

    !Determine which type of topology file to readin
    SELECT CASE (setup%topology_type)
    CASE ('PSF')
      CALL read_topology_psf(mol_setup,setup,globenv)
!   CASE ('GROMOS')
!     CALL read_topology_gromos(mol_setup,setup,globenv)
!   CASE ('GENERATE')
!     !Later on, in case of reading something in like a PDB
!     CALL topology_generate_bonds()
!     CALL topology_generate_molecule()
!     CALL topology_generate_bends()
!     CALL topology_generate_dihe()
    CASE DEFAULT
      CALL stop_program (" read_topology_file ","Unknown type of topology file")
    END SELECT
  END IF

! send the outcome of the inquiry to all cpus
! CALL mp_bcast(topo_file % active, globenv%source, globenv%group)
! CALL mp_bcast(topo_file % type, globenv%source, globenv%group)

END SUBROUTINE read_topology_control

!******************************************************************************
!!****** topology/read_coordinate_control [1.0] *
!!
!!   NAME
!!     read_coordinate_control
!!
!!   FUNCTION
!!     none
!!
!!   AUTHOR
!!     IKUO 08.11.2003
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE read_coordinate_control (particle_set,atomic_kind_set,mol_setup,&
                                    atom_names,setup,globenv)

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT(IN) :: mol_setup
  TYPE ( setup_parameters_type ), INTENT (IN) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
  CHARACTER ( len = * ), DIMENSION ( : ), INTENT ( IN ) :: atom_names

! Locals
  LOGICAL found
  INTEGER :: iw, ierr

  iw = globenv%scr

  IF(setup%coordinate) THEN
    !Make sure the coordinate file is there
    IF(globenv%ionode) THEN
      INQUIRE(FILE=setup%coord_file_name,EXIST=found)
      IF(found) THEN
        WRITE(iw,*) 'Using coordinate file ',setup%coord_file_name
      ELSE
        CALL stop_program ("read_coordinate_control ","coordinate file missing")
      END IF
    END IF

    !Determine which type of topology file to readin
    SELECT CASE (setup%coordinate_type)
    CASE ('PDB')
      CALL read_coordinate_pdb (particle_set,atomic_kind_set,mol_setup,&
                                atom_names,setup,globenv)
      !IF((setup%topology==.TRUE).AND.(setup%topology_type="GENERATE"))THEN
      !  WRITE(iw,*) "GENERATE TOPOLOGY FROM COORDINATES"
      !  CALL topology_generate_bonds()
      !  CALL topology_generate_molecule()
      !  CALL topology_generate_bends()
      !  CALL topology_generate_dihe()
      !END IF
    CASE ('CHM')
      !CALL read_coordinate_chm(mol_setup,setup,globenv)
      !IF((setup%topology==.TRUE).AND.(setup%topology_type="GENERATE"))THEN
      !  WRITE(iw,*) "GENERATE TOPOLOGY FROM COORDINATES"
      !  CALL topology_generate_bonds()
      !  CALL topology_generate_molecule()
      !  CALL topology_generate_bends()
      !  CALL topology_generate_dihe()
      !END IF
    CASE DEFAULT
      CALL stop_program (" read_coordinate_file ",&
                         "Unknown type of coordinate file")
    END SELECT
  END IF

! send the outcome of the inquiry to all cpus
! CALL mp_bcast(topo_file % active, globenv%source, globenv%group)
! CALL mp_bcast(topo_file % type, globenv%source, globenv%group)

END SUBROUTINE read_coordinate_control

END MODULE topology_control

!******************************************************************************
