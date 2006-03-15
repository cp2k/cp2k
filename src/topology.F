!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/topology [1.0] *
!!
!!   NAME
!!     topology
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
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atoms_input,                     ONLY: read_atoms_input
  USE cell_types,                      ONLY: read_cell
  USE colvar_types,                    ONLY: colvar_p_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_conn_generate,&
                                             do_conn_off,&
                                             do_conn_psf,&
                                             do_conn_mol_psf,&
                                             do_coord_pdb,&
                                             do_coord_xyz
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_topology_util,              ONLY: qmmm_connectivity_control,&
                                             qmmm_coordinate_control
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_input,                  ONLY: read_constraints_section,&
                                             read_topology_section
  USE topology_pdb,                    ONLY: read_coordinate_pdb,&
                                             write_coordinate_pdb
  USE topology_psf,                    ONLY: read_topology_psf,&
                                             write_topology_psf
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: &
       check_subsys_element, topology_connectivity_pack, &
       topology_constraint_pack, topology_coordinate_pack, &
       topology_generate_bend, topology_generate_bond, &
       topology_generate_dihe, topology_generate_impr, &
       topology_generate_molecule, topology_generate_onfo, &
       topology_generate_ub, topology_molecules_check, &
       topology_reorder_atoms, topology_set_atm_mass
  USE topology_xyz,                    ONLY: read_coordinate_xyz
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology'

  PRIVATE

! *** Public parameters ***
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
       molecule_kind_set,molecule_set,colvar_p,globenv, qmmm, qmmm_env, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(colvar_p_type), DIMENSION(:), POINTER :: colvar_p
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure, my_qmmm
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: cell_section, &
                                                constraint_section, &
                                                subsys_section, &
                                                topology_section
    TYPE(topology_parameters_type)           :: topology

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,'I','',handle)

    NULLIFY(cell_section,subsys_section)

    subsys_section   =>  section_vals_get_subs_vals(globenv%input_file,&
         "FORCE_EVAL%SUBSYS",error=error)
    cell_section => section_vals_get_subs_vals(subsys_section,&
         "CELL",error=error)
    constraint_section =>  section_vals_get_subs_vals(subsys_section,&
         "CONSTRAINT",error=error)
    topology_section   =>  section_vals_get_subs_vals(subsys_section,&
         "TOPOLOGY",error=error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO",&
         extension=".mmLog",error=error)
    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Initialize the topology structure type
    !-----------------------------------------------------------------------------
    CALL init_topology(topology)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Get the cell info
    !-----------------------------------------------------------------------------
    CALL read_cell ( topology%cell, topology%cell_ref, globenv, error=error,&
         cell_section=cell_section)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Read in the topology section in the input file if any
    !-----------------------------------------------------------------------------
    CALL read_topology_section(topology,topology_section,globenv,error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Read in the constraints section
    !-----------------------------------------------------------------------------
    CALL read_constraints_section(topology,globenv,colvar_p,constraint_section,error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Read in the coordinates
    !-----------------------------------------------------------------------------
    CALL coordinate_control(topology,globenv,subsys_section,error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 6. Read in or generate the molecular connectivity
    !-----------------------------------------------------------------------------
    CALL connectivity_control(topology,globenv,my_qmmm, qmmm_env, subsys_section,error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 7. Pack everything into the molecular types
    !-----------------------------------------------------------------------------
    CALL topology_connectivity_pack(molecule_kind_set,molecule_set,&
         topology,globenv,subsys_section,error)
    !-----------------------------------------------------------------------------
    ! Teo 12.2004:
    ! 
    ! This part takes care of the molecule in which QM atoms were defined.
    ! Preliminary setup for QM/MM link region
    !-----------------------------------------------------------------------------
    ! 7.bis. Set up the QM/MM linkage (if any)
    !-----------------------------------------------------------------------------
    IF (my_qmmm) THEN
       CALL qmmm_connectivity_control(topology, globenv, molecule_set, &
            qmmm_env, error)
    END IF
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 8. Pack everything into the atomic types
    !-----------------------------------------------------------------------------
    IF (my_qmmm) THEN
       CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
            molecule_kind_set,molecule_set,&
            topology,globenv, my_qmmm, qmmm_env, subsys_section,error=error)
    ELSE 
       CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
            molecule_kind_set,molecule_set,&
            topology,globenv,subsys_section=subsys_section,error=error)
    END IF
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 9. Deal with the constraint stuff if requested
    !-----------------------------------------------------------------------------
    IF (my_qmmm) THEN
       CALL topology_constraint_pack(molecule_kind_set,molecule_set,&
            topology,qmmm_env,particle_set,globenv,subsys_section,error)
    ELSE
       CALL topology_constraint_pack(molecule_kind_set,molecule_set,&
            topology,particle_set=particle_set, globenv=globenv,&
            subsys_section=subsys_section,error=error)
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 10. Dump the topology informations
    !-----------------------------------------------------------------------------
    IF(topology%dump_topology) THEN
       CALL write_topology_psf(topology,globenv,subsys_section,error)
       CALL write_coordinate_pdb(topology,globenv,subsys_section,error)
    END IF


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 11. Cleanup the topology structure type
    !-----------------------------------------------------------------------------
    CALL deallocate_topology(topology, error)
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO",error=error)  
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
  SUBROUTINE connectivity_control(topology,globenv,qmmm,qmmm_env,subsys_section,error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    LOGICAL, INTENT(in), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'connectivity_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw
    INTEGER                                  :: i_rep, n_rep,imol,nmol
    LOGICAL                                  :: failure, found, my_qmmm
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(len=default_string_length)     :: filename

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO",&
         extension=".mmLog",error=error)
    CALL timeset(routineN,'I','',handle)

    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. If reading in from external file, make sure its there first
    !-----------------------------------------------------------------------------
!    IF(topology%connectivity) THEN
!       !------------------------------------------------------
!       !Make sure the file really do exist before proceding
!       IF(globenv%ionode) THEN
!
!          INQUIRE(FILE=topology%conn_file_name,EXIST=found)
!          IF(found) THEN
!             IF(iw>0) WRITE(iw,*) 'Using connectivity file ',topology%conn_file_name
!          ELSE
!
!             section => section_vals_get_subs_vals(globenv%input_file,&
!                        "FORCE_EVAL%SUBSYS%TOPOLOGY%MOL_SET",error=error)
!             CALL section_vals_get(section,n_repetition=n_rep,error=error)
!             DO i_rep=1,n_rep
!                CALL section_vals_val_get(section,"CONN_FILE",&
!                          c_val=filename,i_rep_section=i_rep,error=error)
!                INQUIRE(FILE=file_name,EXIST=found)
!                IF(.NOT.found) THEN
!                  CALL stop_program ("connectivity_control","connectivity file missing")
!             END DO
!          END IF
!       END IF
!    END IF


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Read in or generate the connectivity information any way we can 
    !-----------------------------------------------------------------------------
    SELECT CASE (topology%conn_type)
    CASE (do_conn_mol_psf)
       section => section_vals_get_subs_vals(globenv%input_file,&
                  "FORCE_EVAL%SUBSYS%TOPOLOGY%MOL_SET",error=error)
       CALL section_vals_get(section,n_repetition=n_rep,error=error)
       DO i_rep=1,n_rep
         CALL section_vals_val_get(section,"NMOL",i_val=nmol,&
                                    i_rep_section=i_rep,error=error)
         CALL section_vals_val_get(section,"CONN_FILE",c_val=filename,&
                                    i_rep_section=i_rep,error=error)
         DO imol=1,nmol
           CALL read_topology_psf(filename,topology,globenv,subsys_section,error)
         END DO
       END DO
       CALL topology_generate_molecule(topology,my_qmmm,qmmm_env,globenv,subsys_section,error)
       IF(topology%molecules_check) &
            CALL topology_molecules_check(topology,globenv,subsys_section,error)
    CASE (do_conn_psf)
       CALL read_topology_psf(topology%conn_file_name,topology,globenv,subsys_section,error)
       CALL topology_generate_molecule(topology,my_qmmm,qmmm_env,globenv,subsys_section,error)
       IF(topology%molecules_check) &
            CALL topology_molecules_check(topology,globenv,subsys_section,error)
    CASE (do_conn_generate)
       CALL topology_set_atm_mass(topology,globenv,subsys_section,error)
       CALL topology_generate_bond(topology,globenv,subsys_section,error)
       CALL topology_generate_bend(topology,globenv,subsys_section,error)
       CALL topology_generate_ub(topology,globenv,subsys_section,error)
       CALL topology_generate_dihe(topology,globenv,subsys_section,error)
       CALL topology_generate_impr(topology,globenv,subsys_section,error)
       CALL topology_generate_onfo(topology,globenv,subsys_section,error)
       CALL topology_generate_molecule(topology,my_qmmm,qmmm_env,globenv,subsys_section,error)
       IF(topology%molecules_check) &
            CALL topology_molecules_check(topology,globenv,subsys_section,error)
       IF(topology%reorder_atom) THEN
          CALL topology_reorder_atoms(topology,globenv,subsys_section,error)
          CALL topology_set_atm_mass(topology,globenv,subsys_section,error)
          CALL topology_generate_bond(topology,globenv,subsys_section,error)
          CALL topology_generate_bend(topology,globenv,subsys_section,error)
          CALL topology_generate_ub(topology,globenv,subsys_section,error)
          CALL topology_generate_dihe(topology,globenv,subsys_section,error)
          CALL topology_generate_impr(topology,globenv,subsys_section,error)
          CALL topology_generate_onfo(topology,globenv,subsys_section,error)
          CALL topology_generate_molecule(topology,my_qmmm,qmmm_env,globenv,subsys_section,error)
       END IF
    CASE (do_conn_off)
       CALL topology_set_atm_mass(topology,globenv,subsys_section,error)
       CALL topology_generate_bond(topology,globenv,subsys_section,error) 
       CALL topology_generate_bend(topology,globenv,subsys_section,error)
       CALL topology_generate_ub(topology,globenv,subsys_section,error)
       CALL topology_generate_dihe(topology,globenv,subsys_section,error)
       CALL topology_generate_impr(topology,globenv,subsys_section,error)
       CALL topology_generate_onfo(topology,globenv,subsys_section,error)
       CALL topology_generate_molecule(topology,my_qmmm,qmmm_env,globenv,subsys_section,error)
    CASE DEFAULT
       CALL stop_program ("connectivity_control",&
            "no default")
    END SELECT

    IF (my_qmmm) THEN
       CALL qmmm_coordinate_control(topology, globenv, qmmm_env,subsys_section, error)
    END IF
    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO",error=error)

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
  SUBROUTINE coordinate_control(topology,globenv,subsys_section,error)


    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coordinate_control', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_element
    INTEGER                                  :: handle, i, id, iw
    LOGICAL                                  :: failure, found
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO",&
         extension=".mmLog",error=error)
    CALL timeset(routineN,'I','',handle)


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
    CASE (do_coord_pdb)
       CALL read_coordinate_pdb (topology,globenv,subsys_section,error)
       CALL read_atoms_input ( topology, globenv, overwrite=.TRUE., subsys_section=subsys_section,error=error)
    CASE (do_coord_xyz)
       CALL read_coordinate_xyz (topology,globenv,subsys_section,error)
       CALL read_atoms_input ( topology, globenv, overwrite=.TRUE., subsys_section=subsys_section,error=error)
    CASE DEFAULT
       CALL read_atoms_input ( topology, globenv,  subsys_section=subsys_section, error=error)
    END SELECT

    ! if the e.g. the &COORD block is empty, we don't proceed
    IF (topology%natoms.LE.0) THEN
       CALL stop_program ("coordinate_control","No atomic coordinates have been found")
    ENDIF

    ! Fix element name for QS
    DO i = 1, topology%natoms
       id = INDEX(topology%atom_info%element(i),"_")-1
       IF ( id == -1 ) id = LEN_TRIM(topology%atom_info%element(i))
       CALL check_subsys_element(topology%atom_info%element(i)(1:id),my_element,globenv,subsys_section,error)     
       topology%atom_info%element(i)=my_element 
    END DO


    CALL timestop(0.0_dp,handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO",error=error)
  END SUBROUTINE coordinate_control

END MODULE topology

!******************************************************************************
