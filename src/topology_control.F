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
  USE topology_util,                    ONLY : topology_coord_pack,&
                                               topology_generate_bond,&
                                               topology_generate_bend,&
                                               topology_generate_dihe,&
                                               topology_generate_molecule,&
                                               topology_pack

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
  CHARACTER (LEN=20), POINTER :: label_resname(:),label_atmname(:)
  REAL(dbl), DIMENSION(:,:), POINTER :: r
  LOGICAL found
  INTEGER :: iw, ierr
  INTEGER :: natom,nbond
  INTEGER, POINTER :: map_mol_typ(:) !map_mol_typ(iatom)= mol_typ number
  INTEGER, POINTER :: map_mol_num(:) !map_mol_num(iatom)= mol number
  INTEGER, POINTER :: bond_a(:),bond_b(:)
  INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)
  INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)
  REAL (dbl), POINTER :: atm_charge(:)
  REAL (dbl), POINTER :: atm_mass(:)

  iw = globenv%scr

  IF(setup%topology) THEN

    !Determine which type of topology file to readin
    SELECT CASE (setup%topology_type)
    CASE ('PSF')
      !-----------------------------------------------------------------------
      !Make sure the topology file is there
      IF(globenv%ionode) THEN
        INQUIRE(FILE=setup%topology_file_name,EXIST=found)
        IF(found) THEN
          WRITE(iw,*) 'Using topology file ',setup%topology_file_name
          CALL read_topology_psf(mol_setup,setup,globenv)
        ELSE
          CALL stop_program ("read_topology_control","topology file missing")
        END IF
      END IF
      !-----------------------------------------------------------------------
!   CASE ('GROMOS')
!     !-----------------------------------------------------------------------
!     !Make sure the topology file is there
!     IF(globenv%ionode) THEN
!       INQUIRE(FILE=setup%topology_file_name,EXIST=found)
!       IF(found) THEN
!         WRITE(iw,*) 'Using topology file ',setup%topology_file_name
!         CALL read_topology_gromos(mol_setup,setup,globenv)
!       ELSE
!         CALL stop_program ("read_topology_control","topology file missing")
!       END IF
!     END IF
!     !-----------------------------------------------------------------------
    CASE ('GENERATE')
      !-----------------------------------------------------------------------
      !Read in the coordinates first PDB
      IF(setup%coordinate) THEN
        SELECT CASE (setup%coordinate_type)
        CASE ("PDB")
          CALL read_coordinate_pdb (r,label_resname,label_atmname,&
                                    setup,globenv)
        CASE ("CHM")
          !CALL read_coordinate_chm (r,label_resname,label_atmname,&
          !                          setup,globenv)
          CALL stop_program (" read_topology_control ",&
                             "chm coordinate readin not implemented yet")
        CASE DEFAULT
          CALL stop_program (" read_topology_control ",&
                             "Unknown type of coordinate file")
        END SELECT
      ELSE
        CALL stop_program ("read_topology_control",&
                           "setup%coordinate is .FALSE.")
      END IF

      CALL topology_generate_bond(setup%cell,r,label_atmname,bond_a,bond_b,&
                                  atm_mass,atm_charge)
      CALL topology_generate_bend(bond_a,bond_b,theta_a,theta_b,theta_c)
      CALL topology_generate_dihe(theta_a,theta_b,theta_c,&
                                  phi_a,phi_b,phi_c,phi_d)
      natom = SIZE(r(1,:))
      nbond = SIZE(bond_a)
      ALLOCATE (map_mol_typ(natom),STAT=ierr)
        IF(ierr/=0) CALL stop_memory &
          ('read_topology_psf','map_mol_typ',natom)
      ALLOCATE (map_mol_num(natom),STAT=ierr)
        IF(ierr/=0) CALL stop_memory &
          ('read_topology_psf','map_mol_num',natom)
      CALL topology_generate_molecule(natom,label_resname,nbond,&
                                      bond_a,bond_b,&
                                      map_mol_typ,map_mol_num)
      CALL topology_pack(bond_a,bond_b,theta_a,theta_b,theta_c,&
                         phi_a,phi_b,phi_c,phi_d,&
                         label_atmname,label_resname,atm_mass,atm_charge,&
                         map_mol_typ,map_mol_num,&
                         mol_setup,setup,globenv)
      !CALL stop_program ("read_topology_control ","checkpointfile")
      !-----------------------------------------------------------------------
    CASE DEFAULT
      CALL stop_program ("read_topology_control ",&
                         "Unknown type of topology file")
    END SELECT
  END IF

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
                                    setup,globenv)

  IMPLICIT NONE

! Arguments
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT(IN) :: mol_setup
  TYPE ( setup_parameters_type ), INTENT (IN) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
  TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

! Locals
  LOGICAL found
  CHARACTER (LEN=20), POINTER :: label_resname(:),label_atmname(:)
  REAL(dbl), DIMENSION(:,:), POINTER :: r
  INTEGER :: iw, ierr
  INTEGER :: i, natom

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
      CALL read_coordinate_pdb (r,label_resname,label_atmname,&
                                setup,globenv)
      CALL topology_coord_pack (r,label_atmname,particle_set,&
                                atomic_kind_set,mol_setup,setup,globenv)
    CASE ('CHM')
      !CALL read_coordinate_chm (r,label_resname,label_atmname,&
      !                          setup,globenv)
      !CALL topology_coord_pack (r,label_atmname,particle_set,&
      !                          atomic_kind_set,mol_setup,setup,globenv)
      CALL stop_program (" read_coordinate_control ",&
                         "Not implemented yet...")
    CASE DEFAULT
      CALL stop_program (" read_coordinate_control ",&
                         "Unknown type of coordinate file")
    END SELECT
  END IF

END SUBROUTINE read_coordinate_control

END MODULE topology_control

!******************************************************************************
