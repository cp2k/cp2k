!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Control for reading in different topologies and coordinates
!> \par History
!>      none
! *****************************************************************************
MODULE topology
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atoms_input,                     ONLY: read_atoms_input
  USE cell_types,                      ONLY: cell_create,&
                                             cell_retain,&
                                             cell_type,&
                                             read_cell,&
                                             write_cell
  USE colvar_types,                    ONLY: colvar_p_type,&
                                             post_process_colvar
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE exclusion_types,                 ONLY: exclusion_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_conn_amb7, do_conn_g87, do_conn_g96, do_conn_generate, &
       do_conn_mol_set, do_conn_off, do_conn_psf, do_conn_psf_u, &
       do_conn_user, do_coord_cif, do_coord_cp2k, do_coord_crd, do_coord_g96, &
       do_coord_off, do_coord_pdb, do_coord_xtl, do_coord_xyz
  USE input_cp2k_binary_restarts,      ONLY: read_binary_coordinates
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine
  USE mm_mapping_library,              ONLY: create_ff_map,&
                                             destroy_ff_map
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_topology_util,              ONLY: qmmm_connectivity_control,&
                                             qmmm_coordinate_control
  USE qmmm_types,                      ONLY: qmmm_env_mm_type
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_amber,                  ONLY: read_connectivity_amber,&
                                             read_coordinate_crd
  USE topology_cif,                    ONLY: read_coordinate_cif
  USE topology_connectivity_util,      ONLY: topology_conn_multiple,&
                                             topology_connectivity_pack
  USE topology_constraint_util,        ONLY: topology_constraint_pack
  USE topology_coordinate_util,        ONLY: topology_coordinate_pack
  USE topology_cp2k,                   ONLY: read_coordinate_cp2k
  USE topology_generate_util,          ONLY: topology_generate_bend,&
                                             topology_generate_bond,&
                                             topology_generate_dihe,&
                                             topology_generate_impr,&
                                             topology_generate_molecule,&
                                             topology_generate_onfo,&
                                             topology_generate_ub
  USE topology_gromos,                 ONLY: read_coordinate_g96,&
                                             read_topology_gromos
  USE topology_input,                  ONLY: read_constraints_section,&
                                             read_topology_section
  USE topology_multiple_unit_cell,     ONLY: topology_muc
  USE topology_pdb,                    ONLY: read_coordinate_pdb,&
                                             write_coordinate_pdb
  USE topology_psf,                    ONLY: idm_psf,&
                                             psf_post_process,&
                                             read_topology_psf,&
                                             write_topology_psf
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             pre_read_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: check_subsys_element,&
                                             topology_molecules_check,&
                                             topology_reorder_atoms,&
                                             topology_set_atm_mass
  USE topology_xtl,                    ONLY: read_coordinate_xtl
  USE topology_xyz,                    ONLY: read_coordinate_xyz
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology'

  PRIVATE

! *** Public parameters ***
  PUBLIC :: topology_control, &
            connectivity_control, &
            coordinate_control

CONTAINS

! *****************************************************************************
  SUBROUTINE topology_control(atomic_kind_set, particle_set, molecule_kind_set,&
       molecule_set, colvar_p, gci, root_section, para_env, qmmm, qmmm_env,&
       force_env_section, subsys_section, use_motion_section, exclusions, error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(exclusion_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: exclusions
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'topology_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, iw2
    LOGICAL                                  :: binary_coord_read, failure, &
                                                my_qmmm
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: cell_section, &
                                                constraint_section, &
                                                topology_section
    TYPE(topology_parameters_type)           :: topology

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL timeset(routineN,handle)
    NULLIFY(cell_section, constraint_section, topology_section)

    cell_section => section_vals_get_subs_vals(subsys_section,"CELL",error=error)
    IF (use_motion_section) THEN
       constraint_section =>  section_vals_get_subs_vals(root_section,"MOTION%CONSTRAINT",error=error)
    END IF
    topology_section   =>  section_vals_get_subs_vals(subsys_section,"TOPOLOGY",error=error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO",&
         extension=".mmLog",error=error)
    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm

    ! 1. Initialize the topology structure type
    CALL init_topology(topology, error)

    ! 2. Get the cell info
    CALL read_cell (topology%cell, topology%cell_ref, cell_section=cell_section,&
                    para_env=para_env, error=error)
    CALL write_cell(topology%cell, subsys_section, label="CELL_TOP|", error=error)
    CALL setup_cell_muc(topology%cell_muc, topology%cell, subsys_section, error)

    ! 3. Read in the topology section in the input file if any
    CALL read_topology_section(topology,topology_section,error)

    ! 4. Read in the constraints section
    CALL read_constraints_section(topology,colvar_p,constraint_section,error)

    ! 5. Read in the coordinates
    CALL read_binary_coordinates(topology,root_section,para_env,subsys_section,&
                                 binary_coord_read,error)
    IF (.NOT.binary_coord_read) THEN
       CALL coordinate_control(topology,root_section,para_env,subsys_section,error)
    END IF

    ! 6. Read in or generate the molecular connectivity
    CALL connectivity_control(topology,para_env,my_qmmm, qmmm_env,subsys_section,&
         force_env_section,error)

    ! 7. Pack everything into the molecular types
    CALL topology_connectivity_pack(molecule_kind_set,molecule_set,&
         topology,subsys_section,error)

    ! 8. Set up the QM/MM linkage (if any)
    !    This part takes care of the molecule in which QM atoms were defined.
    !    Preliminary setup for QM/MM link region
    IF (my_qmmm) THEN
       CALL qmmm_connectivity_control(molecule_set, qmmm_env, subsys_section, error)
    END IF

    ! 9. Pack everything into the atomic types
    IF (my_qmmm) THEN
       CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
            molecule_kind_set,molecule_set,&
            topology, my_qmmm, qmmm_env, subsys_section,&
            force_env_section=force_env_section, exclusions=exclusions,&
            error=error)
    ELSE
       CALL topology_coordinate_pack(particle_set,atomic_kind_set,&
            molecule_kind_set,molecule_set,&
            topology,subsys_section=subsys_section,&
            force_env_section=force_env_section, exclusions=exclusions, &
            error=error)
    END IF

    !10. Post-Process colvar definitions (if needed)
    CALL topology_post_proc_colvar(colvar_p, particle_set, error)

    !11. Deal with the constraint stuff if requested
    IF (my_qmmm) THEN
       CALL topology_constraint_pack(molecule_kind_set,molecule_set,&
            topology,qmmm_env,particle_set,root_section,subsys_section,&
            gci,error)
    ELSE
       CALL topology_constraint_pack(molecule_kind_set,molecule_set,&
            topology,particle_set=particle_set,input_file=root_section, &
            subsys_section=subsys_section,gci=gci,error=error)
    END IF

    !12. Dump the topology informations
    iw2 = cp_print_key_unit_nr(logger,subsys_section,"TOPOLOGY%DUMP_PDB",&
          file_status="REPLACE",extension=".pdb",error=error)
    IF (iw2>0) THEN
       CALL write_coordinate_pdb(iw2,topology,subsys_section,error)
    END IF
    CALL cp_print_key_finished_output(iw2,logger,subsys_section,&
         "TOPOLOGY%DUMP_PDB",error=error)
    iw2 = cp_print_key_unit_nr(logger,subsys_section,"TOPOLOGY%DUMP_PSF",&
          file_status="REPLACE",extension=".psf",error=error)
    IF (iw2>0) THEN
       CALL write_topology_psf(iw2,topology,subsys_section,force_env_section,error)
    END IF
    CALL cp_print_key_finished_output(iw2,logger,subsys_section,&
         "TOPOLOGY%DUMP_PSF",error=error)

    !13. Cleanup the topology structure type
    CALL deallocate_topology(topology, error)
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO",error=error)
  END SUBROUTINE topology_control

! *****************************************************************************
!> \brief 1. If reading in from external file, make sure its there first
!>      2. Generate the connectivity if no information to be read in
!> \par History
!>      none
!> \author IKUO 08.01.2003
! *****************************************************************************
  SUBROUTINE connectivity_control(topology,para_env,qmmm,qmmm_env,subsys_section,&
       force_env_section, error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(in), OPTIONAL            :: qmmm
    TYPE(qmmm_env_mm_type), OPTIONAL, &
      POINTER                                :: qmmm_env
    TYPE(section_vals_type), POINTER         :: subsys_section, &
                                                force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'connectivity_control', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: map0 = ICHAR("0"), &
                                                map9 = ICHAR("9")

    CHARACTER(len=default_string_length)     :: element0, my_element
    CHARACTER(len=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: elements
    INTEGER                                  :: handle, handle2, i, id, itmp, &
                                                iw, j, k, stat
    LOGICAL                                  :: check, failure, my_qmmm, &
                                                use_mm_map_first
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO",&
         extension=".mmLog",error=error)
    CALL timeset(routineN,handle)

    my_qmmm = .FALSE.
    IF (PRESENT(qmmm).AND.PRESENT(qmmm_env)) my_qmmm = qmmm

    ! 1. Read in the connectivity information (if this is the case)
    SELECT CASE (topology%conn_type)
    CASE (do_conn_generate, do_conn_off, do_conn_user)
       ! Do nothing for the time being.. after we check element and proceed with the workflow..
    CASE DEFAULT
       ! Prepare arrays
       CALL pre_read_topology(topology,error)

       ! Read connectivity from file
       CALL read_topology_conn(topology, topology%conn_type, topology%conn_file_name, &
                               para_env, subsys_section, error)

       ! Post process of PSF and AMBER information
       SELECT CASE (topology%conn_type)
       CASE (do_conn_mol_set, do_conn_psf, do_conn_psf_u, do_conn_amb7)
          CALL psf_post_process(topology,subsys_section,error)
       END SELECT
    END SELECT

    ! 2. In case element was autoassigned let's keep up2date the element name
    !    with the atom_name
    IF (topology%aa_element) THEN
       check = SIZE(topology%atom_info%id_element)==SIZE(topology%atom_info%id_atmname)
       CPPostcondition(check,cp_failure_level,routineP,error,failure)
       topology%atom_info%id_element = topology%atom_info%id_atmname
    END IF

    ! 3. Check for the element name..
    CALL timeset(routineN//"_check_element_name",handle2)
    ! Fix element name

    ! we will only translate names if we actually have a connectivity file given
    SELECT CASE(topology%conn_type)
    CASE(do_conn_mol_set,do_conn_psf, do_conn_psf_u, do_conn_g96, do_conn_g87, do_conn_amb7)
        use_mm_map_first=.TRUE.
    CASE DEFAULT
        use_mm_map_first=.FALSE.
    END SELECT
    CALL create_ff_map("AMBER",error)
    CALL create_ff_map("CHARMM",error)
    CALL create_ff_map("GROMOS",error)

    ALLOCATE(elements(SIZE(topology%atom_info%id_element)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i=1,SIZE(elements)
       elements(i) = id2str(topology%atom_info%id_element(i))
    ENDDO

    DO i = 1, topology%natoms
       IF (elements(i)=="__DEF__") CYCLE
       ! If present an underscore let's skip all that over the underscore
       id = INDEX(elements(i),"_")-1
       IF (id==-1) id = LEN_TRIM(elements(i))
       ! Many atomic kind have been defined as ELEMENT+LETTER+NUMBER
       ! the number at the end can vary arbitrarily..
       ! Let's check all ELEMENT+LETTER skipping the number.. we should
       ! identify in any case the element
       DO j = id, 1, -1
          itmp = ICHAR(elements(i)(j:j))
          IF ((itmp<map0).OR.(itmp>map9)) EXIT
       END DO
       element0 = elements(i)(1:j)
       ! ALWAYS check for elements..
       CALL check_subsys_element(element0,id2str(topology%atom_info%id_atmname(i)),my_element,&
            subsys_section,use_mm_map_first,error)
       ! Earn time fixing same element labels for same atoms
       element0 = elements(i)
       DO k = i, topology%natoms
          IF (element0==id2str(topology%atom_info%id_element(k))) THEN
             topology%atom_info%id_element(k)=str2id(s2s(my_element))
             elements(k)="__DEF__"
          END IF
       END DO
    END DO
    DEALLOCATE(elements,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL destroy_ff_map("GROMOS",error)
    CALL destroy_ff_map("CHARMM",error)
    CALL destroy_ff_map("AMBER",error)
    CALL timestop(handle2)

    ! 4. Generate the connectivity information otherwise
    SELECT CASE (topology%conn_type)
    CASE (do_conn_generate)
       CALL topology_set_atm_mass(topology,subsys_section,error)
       CALL topology_generate_bond(topology,para_env,subsys_section,error)
       IF(topology%reorder_atom) THEN
          ! If we generate connectivity we can save memory reordering the molecules
          ! in this case once a first connectivity has been created we match according
          ! molecule names provided in the PDB and reorder the connectivity according to that.
          CALL topology_reorder_atoms(topology,qmmm,qmmm_env,subsys_section,&
                                      force_env_section,error=error)
          CALL topology_set_atm_mass(topology,subsys_section,error)
          CALL topology_generate_bond(topology,para_env,subsys_section,error)
       END IF
       CALL topology_generate_bend(topology,subsys_section,error)
       CALL topology_generate_ub(topology,subsys_section,error)
       CALL topology_generate_dihe(topology,subsys_section,error)
       CALL topology_generate_impr(topology,subsys_section,error)
       CALL topology_generate_onfo(topology,subsys_section,error)
    CASE (do_conn_off, do_conn_user)
       CALL topology_set_atm_mass(topology,subsys_section,error)
       CALL topology_generate_bend(topology,subsys_section,error)
       CALL topology_generate_ub(topology,subsys_section,error)
       CALL topology_generate_dihe(topology,subsys_section,error)
       CALL topology_generate_impr(topology,subsys_section,error)
       CALL topology_generate_onfo(topology,subsys_section,error)
    END SELECT

    ! 5. Handle multiple unit_cell - Update atoms_info
    CALL topology_muc(topology, subsys_section, error=error)

    ! 6. Handle multiple unit_cell - Update conn_info
    CALL topology_conn_multiple(topology, subsys_section, error)

    ! 7. Generate Molecules
    CALL topology_generate_molecule(topology,my_qmmm,qmmm_env,subsys_section,error)
    IF(topology%molecules_check) CALL topology_molecules_check(topology,subsys_section,error)

    ! 8. Modify for QM/MM
    IF (my_qmmm) THEN
       CALL qmmm_coordinate_control(topology, qmmm_env,subsys_section, error)
    END IF
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO",error=error)

  END SUBROUTINE connectivity_control

! *****************************************************************************
!> \brief Reads connectivity from file
!> \author Teodoro Laino [tlaino] - 10.2009
! *****************************************************************************
  RECURSIVE SUBROUTINE read_topology_conn(topology, conn_type, conn_file_name, para_env,&
       subsys_section, error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    INTEGER, INTENT(IN)                      :: conn_type
    CHARACTER(LEN=default_path_length), &
      INTENT(IN)                             :: conn_file_name
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_topology_conn', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename
    INTEGER                                  :: i_rep, imol, loc_conn_type, &
                                                n_rep, nmol
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: section

    failure = .FALSE.
    NULLIFY (section)

    SELECT CASE (conn_type)
    CASE (do_conn_mol_set)
       section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%MOL_SET",error=error)
       section => section_vals_get_subs_vals(section,"MOLECULE",error=error)
       CALL section_vals_get(section,n_repetition=n_rep,error=error)
       DO i_rep=1,n_rep
          CALL section_vals_val_get(section,"NMOL",i_val=nmol,i_rep_section=i_rep,error=error)
          CALL section_vals_val_get(section,"CONN_FILE_NAME",c_val=filename,i_rep_section=i_rep,error=error)
          CALL section_vals_val_get(section,"CONN_FILE_FORMAT",i_val=loc_conn_type,i_rep_section=i_rep,error=error)

          SELECT CASE (loc_conn_type)
          CASE (do_conn_psf, do_conn_psf_u, do_conn_g96, do_conn_g87, do_conn_amb7)
             DO imol=1,nmol
                CALL read_topology_conn(topology,loc_conn_type,filename,para_env,subsys_section,&
                     error)
             END DO
          CASE DEFAULT
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "MOL_SET feature implemented only for PSF/UPSF, G87/G96 and AMBER "//&
                  "connectivity type."//&
CPSourceFileRef,&
                  error=error,failure=failure)
          END SELECT
       END DO
       CALL cp_assert(SIZE(topology%atom_info%id_molname)==topology%natoms,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "Number of atoms in connectivity control is larger than the "//&
            "number of atoms in coordinate control. check coordinates and "//&
            "connectivity. "//&
 CPSourceFileRef,&
            error=error,failure=failure)

       ! Merge defined structures
       section => section_vals_get_subs_vals(subsys_section,"TOPOLOGY%MOL_SET%MERGE_MOLECULES",error=error)
       CALL idm_psf(topology,section,subsys_section,error)

    CASE (do_conn_g96,do_conn_g87)
       CALL read_topology_gromos(conn_file_name,topology,para_env,subsys_section,error)
    CASE (do_conn_psf, do_conn_psf_u)
       CALL read_topology_psf(conn_file_name,topology,para_env,subsys_section,conn_type, error)
    CASE (do_conn_amb7)
       CALL read_connectivity_amber(conn_file_name,topology,para_env,subsys_section,error)
    END SELECT

  END SUBROUTINE read_topology_conn

! *****************************************************************************
!> \brief 1. If reading in from external file, make sure its there first
!>      2. Read in the coordinates from the corresponding locations
!> \par History
!>      - Teodoro Laino [tlaino] - University of Zurich 10.2008
!>        adding support for AMBER coordinates
!> \author IKUO 08.11.2003
! *****************************************************************************
  SUBROUTINE coordinate_control(topology,root_section,para_env,subsys_section,error)

    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coordinate_control', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: message
    INTEGER                                  :: handle, handle2, istat, iw
    LOGICAL                                  :: failure, found, save_mem
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: global_section

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO",&
         extension=".mmLog",error=error)
    CALL timeset(routineN,handle)

    NULLIFY(global_section)
    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. If reading in from external file, make sure its there first
    !-----------------------------------------------------------------------------
    IF (topology%coordinate) THEN
       INQUIRE (FILE=topology%coord_file_name,EXIST=found,IOSTAT=istat)
       IF (istat /= 0) THEN
          WRITE (UNIT=message,FMT="(A,I0,A)")&
            "An error occurred inquiring the file <"//&
            TRIM(topology%coord_file_name)//"> (IOSTAT = ",istat,")"
          CALL stop_program(routineN,moduleN,__LINE__,TRIM(message))
       END IF
       IF (.NOT.found) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Coordinate file <"//TRIM(topology%coord_file_name)//&
                            "> not found.",para_env)
       END IF
    END IF
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Read in the coordinates from the corresponding locations
    !-----------------------------------------------------------------------------
    CALL timeset(routineN//"_READ_COORDINATE",handle2)
    SELECT CASE (topology%coord_type)
    CASE (do_coord_off)
       ! Do nothing.. we will parse later from the &COORD section..
    CASE (do_coord_g96)
       CALL read_coordinate_g96 (topology,para_env,subsys_section,error)
    CASE (do_coord_crd)
       CALL read_coordinate_crd (topology,para_env,subsys_section,error)
    CASE (do_coord_pdb)
       CALL read_coordinate_pdb (topology,para_env,subsys_section,error)
    CASE (do_coord_xyz)
       CALL read_coordinate_xyz (topology,para_env,subsys_section,error)
    CASE (do_coord_cif)
       CALL read_coordinate_cif (topology,para_env,subsys_section,error)
    CASE (do_coord_xtl)
       CALL read_coordinate_xtl (topology,para_env,subsys_section,error)
    CASE (do_coord_cp2k)
       CALL read_coordinate_cp2k(topology,para_env,subsys_section,error)
    CASE DEFAULT
       ! We should never reach this point..
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    ! Parse &COORD section and in case overwrite
    IF (topology%coord_type /= do_coord_cp2k) THEN
       CALL read_atoms_input(topology,overwrite=(topology%coord_type /= do_coord_off),&
                             subsys_section=subsys_section,save_mem=save_mem,error=error)
    END IF
    CALL section_vals_val_set(subsys_section,"TOPOLOGY%NUMBER_OF_ATOMS",&
                              i_val=topology%natoms,error=error)
    CALL timestop(handle2)
    ! Check on atom numbers
    CALL cp_assert(topology%natoms>0,cp_fatal_level,cp_assertion_failed,routineP,&
         "No atomic coordinates have been found! "//&
CPSourceFileRef,&
         error,failure)
    CALL timestop(handle)
    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
         "PRINT%TOPOLOGY_INFO",error=error)
  END SUBROUTINE coordinate_control

! *****************************************************************************
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - 07.2007
! *****************************************************************************
  SUBROUTINE topology_post_proc_colvar(colvar_p, particle_set, error)

    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'topology_post_proc_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (ASSOCIATED(colvar_p)) THEN
       DO i = 1, SIZE(colvar_p)
          CALL post_process_colvar(colvar_p(i)%colvar, particle_set, error)
       END DO
    END IF
  END SUBROUTINE topology_post_proc_colvar

! *****************************************************************************
!> \brief  Setup the cell used for handling properly the multiple_unit_cell option
!> \author Teodoro Laino [tlaino] - 06.2009
! *****************************************************************************
  SUBROUTINE setup_cell_muc(cell_muc, cell, subsys_section, error)

    TYPE(cell_type), POINTER                 :: cell_muc, cell
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_cell_muc', &
      routineP = moduleN//':'//routineN

    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat_ref

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(cell_muc),cp_failure_level,routineP,error,failure)

    CALL section_vals_val_get(subsys_section, "CELL%MULTIPLE_UNIT_CELL",&
                              i_vals=multiple_unit_cell, error=error)
    IF (ANY(multiple_unit_cell/=1)) THEN
       ! Restore the original cell
       hmat_ref(:,1) = cell%hmat(:,1)/multiple_unit_cell(1)
       hmat_ref(:,2) = cell%hmat(:,2)/multiple_unit_cell(2)
       hmat_ref(:,3) = cell%hmat(:,3)/multiple_unit_cell(3)

       ! Create the MUC cell
       CALL cell_create(cell_muc,hmat=hmat_ref,periodic=cell%perd,error=error)
       CALL write_cell(cell_muc,subsys_section,label="CELL_UC|",error=error)
    ELSE
       ! If a multiple_unit_cell was not requested just point to the original cell
       CALL cell_retain(cell, error)
       cell_muc => cell
    END IF

  END SUBROUTINE setup_cell_muc

END MODULE topology

