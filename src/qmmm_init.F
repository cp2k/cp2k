!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_init [1.0] *
!!
!!   NAME
!!     qmmm_init
!!
!!   FUNCTION
!!     Initialize a QM/MM calculation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     5.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_init
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             read_atomic_kind_set, &
                                             set_atomic_kind
  USE cp_output_handling, ONLY: cp_should_output
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_assertion_failed
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_logger_type,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_p_type,&
                                             cp_subsystem_type
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE force_fields,                    ONLY: force_field_control
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_create,&
                                             qmmm_env_qm_release,&
                                             qmmm_env_qm_type,&
                                             qmmm_env_mm_type,&
                                             qmmm_env_mm_create,&
                                             qmmm_env_mm_release,&
                                             qmmm_links_type,&
                                             add_set_type,&
                                             create_add_set_type,&
                                             add_set_release
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type
  USE qmmm_gaussian_init,              ONLY: qmmm_gaussian_initialize
  USE qmmm_pw_grid,                    ONLY: qmmm_pw_grid_init
  USE qmmm_elpot,                      ONLY: qmmm_potential_init
  USE qs_environment,                  ONLY: qs_init_subsys
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_env_release,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_parser,                       ONLY: finish_parser,&
                                             open_file,&
                                             read_object,&
                                             search_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE simulation_cell,                 ONLY: cell_retain,&
                                             cell_type,&
                                             init_qmmm_cell,&
                                             cell_clone,&
                                             cell_release,&
                                             write_cell
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: connectivity_control,&
                                             coordinate_control
  USE topology_pdb,                    ONLY: write_coordinate_pdb
  USE topology_psf,                    ONLY: write_topology_psf
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: topology_connectivity_pack,&
                                             topology_constraint_pack,&
                                             topology_coordinate_pack
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             allocate_potential,&
                                             set_potential,&
                                             get_potential
  USE fist_main,                       ONLY: fist_create_force_env
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_describe_html, section_describe_index_html, section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  USE input_cp2k, ONLY: RADIUS_QMMM_DEFAULT, do_qmmm_link_imomm, do_qmmm_link_gho,&
                        do_qmmm_gauss
  USE fist_environment_types,  ONLY: set_fist_env
  
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_init'

  PUBLIC :: qmmm_create_force_env
!***
!****************************************************************************
CONTAINS

!!****f* qmmm_init/qmmm_create_force_env [1.0] *
!!
!!   NAME
!!     qmmm_create_force_env
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE  qmmm_create_force_env(force_env, globenv, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_create_force_env', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: qm_atom_type
    INTEGER                                  :: handle, isubf, stat
    INTEGER                                  :: iw
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_atom_index, qm_subsys_index
    INTEGER, DIMENSION(:), POINTER           :: mm_link_atoms
    LOGICAL                                  :: failure, qmmm_link
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_mm
    TYPE(cp_subsystem_type), POINTER         :: subsys_qm
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: sub_force_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env_mm
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp), DIMENSION(:), POINTER          :: mm_atom_chrg, mm_el_pot_radius
    REAL(KIND=dp)                            :: eps_mm_rspace
    INTEGER                                  :: qmmm_coupl_type
    TYPE(section_vals_type), POINTER         :: force_env_section, qmmm_section
    TYPE(qmmm_links_type), pointer           :: qmmm_links
    REAL(KIND=dp),  DIMENSION(:), POINTER    :: mm_link_scale_factor, fist_scale_charge_link
    LOGICAL                                  :: move_mm_charges, add_mm_charges
    TYPE(add_set_type), POINTER              :: added_charges

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    NULLIFY(sub_force_env, qm_atom_index, mm_atom_index, qm_atom_type,   &
         qmmm_env_qm, subsys_mm,subsys_qm, mm_cell, qm_cell_small,       &
         qs_env, para_env, mm_atom_chrg,mm_el_pot_radius, qmmm_env_mm,   &
         mm_link_atoms, mm_link_scale_factor, qmmm_links, added_charges, &
         fist_scale_charge_link)
    para_env => globenv%para_env
    iw = globenv%scr
    IF (.NOT. failure) THEN
       ! Input section...
       force_env_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL",error=error)
       qmmm_section => section_vals_get_subs_vals(force_env_section,"QMMM",error=error)

       ALLOCATE(sub_force_env(2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO isubf=1,2
          NULLIFY(sub_force_env(isubf)%force_env)
       END DO
       ! Create QM/MM Environments..
       CALL qmmm_env_qm_create(qmmm_env_qm,error=error)
       CALL qmmm_env_mm_create(qmmm_env_mm,error=error)

       ! Set up QM/MM Options
       CALL setup_qmmm_vars_mm(qmmm_section,&
                               qmmm_env_mm,&
                               qm_atom_index,&
                               mm_link_atoms,&
                               mm_link_scale_factor,&
                               fist_scale_charge_link,&
                               qmmm_coupl_type,&
                               qmmm_link,&
                               error=error)

       qmmm_env_mm%qm_atom_index          => qm_atom_index       
       qmmm_env_mm%mm_link_atoms          => mm_link_atoms
       qmmm_env_mm%mm_link_scale_factor   => mm_link_scale_factor
       qmmm_env_mm%fist_scale_charge_link => fist_scale_charge_link
       qmmm_env_mm%qmmm_coupl_type        =  qmmm_coupl_type
       qmmm_env_mm%qmmm_link              =  qmmm_link
       ! Center the qm subsys into the qm box
       IF (globenv%run_type == "DEBUG") qmmm_env_qm%center_qm_subsys = .FALSE.

       ! First Initialize Fist...
       globenv%program_name="FIST"       
       CALL fist_create_force_env(sub_force_env(1)%force_env, &
            globenv, qmmm=.TRUE., qmmm_env=qmmm_env_mm, error=error)

       CALL force_env_get(sub_force_env(1)%force_env,subsys=subsys_mm,&
            cell=mm_cell,error=error) 
       CPPrecondition(SIZE(subsys_mm)==1,cp_failure_level,routineP,error,failure)

       IF (.NOT.ASSOCIATED(qm_cell_small)) THEN
          ALLOCATE (qm_cell_small,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL cell_clone(mm_cell, qm_cell_small,.true.)

       ! Set up QM/MM Options       
       CALL setup_qmmm_vars_qm (qmmm_section,&
                                qmmm_env_qm,&
                                subsys_mm,&
                                qm_atom_type,&
                                qm_atom_index,&
                                mm_atom_index,&
                                qm_subsys_index,&
                                qm_cell_small,&
                                qmmm_coupl_type,&
                                eps_mm_rspace,&
                                qmmm_link,&
                                error)
                               
       qmmm_env_qm%qm_atom_index   => qm_atom_index
       qmmm_env_qm%mm_atom_index   => mm_atom_index
       qmmm_env_qm%qm_subsys_index => qm_subsys_index
       qmmm_env_qm%eps_mm_rspace   = eps_mm_rspace
       qmmm_env_qm%qmmm_coupl_type = qmmm_coupl_type
       qmmm_env_qm%qmmm_link       = qmmm_link
       qmmm_env_qm%num_qm_atoms = size(qm_atom_index)
       qmmm_env_qm%num_mm_atoms = size(mm_atom_index)
       WRITE(iw,'(A,I6)')" Number of QM atoms: ",qmmm_env_qm%num_qm_atoms 
       WRITE(iw,'(A,I6)')" Number of MM atoms: ",qmmm_env_qm%num_mm_atoms
       WRITE(iw,'(A)')" QM cell ::"
       CALL write_cell(qm_cell_small, globenv)
       !
       ! Assign charges and mm_el_pot_radius from fist_topology
       !
       ALLOCATE(mm_atom_chrg(SIZE(mm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mm_el_pot_radius(SIZE(mm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       mm_atom_chrg     = 0.0_dp
       mm_el_pot_radius = 0.0_dp
       
       CALL assign_mm_charges_and_radius(subsys=subsys_mm,&
                                         mm_atom_chrg=mm_atom_chrg,&
                                         mm_el_pot_radius=mm_el_pot_radius,&
                                         mm_atom_index=mm_atom_index,&
                                         mm_link_atoms=mm_link_atoms,&
                                         mm_link_scale_factor=mm_link_scale_factor,&
                                         error=error)

       qmmm_env_qm%mm_atom_chrg     => mm_atom_chrg
       qmmm_env_qm%mm_el_pot_radius => mm_el_pot_radius

       IF (qmmm_link) THEN 
          CALL setup_qmmm_links(qmmm_section, qmmm_links, error=error)
          qmmm_env_qm%qmmm_links    => qmmm_links

          CALL print_qmmm_links(globenv, qmmm_section, qmmm_links, error)

          CALL add_set_release  (qmmm_env_qm%added_charges, error=error)
          CALL move_or_add_atoms(qmmm_section, move_mm_charges, add_mm_charges, &
                                 mm_atom_chrg, mm_el_pot_radius, added_charges, &
                                 mm_atom_index, qm_cell_small, error)
          qmmm_env_qm%move_mm_charges =  move_mm_charges
          qmmm_env_qm%add_mm_charges  =  add_mm_charges
          qmmm_env_qm%added_charges   => added_charges
       END IF
       
       CALL print_qmmm_charges(mm_atom_index, mm_atom_chrg, mm_el_pot_radius,&
                               qmmm_env_qm%added_charges, globenv, qmmm_section,&
                               error)

       globenv%program_name="QS"
       CALL create_small_subsystem(subsys_qm,&
            big_subsys=subsys_mm,small_para_env=para_env,&
            small_cell=qm_cell_small,sub_atom_index=qm_atom_index,&
            sub_subsys_index=qm_subsys_index,&
            sub_atom_kind_name=qm_atom_type, globenv=globenv,error=error)
       CALL read_atomic_kind_set(subsys_qm%atomic_kinds%els,globenv)
     
       CALL qs_env_create(qs_env, para_env=para_env, error=error)

       CALL set_qs_env(qs_env,input=force_env_section,error=error)

       CALL qs_init_subsys(qs_env,subsystem=subsys_qm,cell=qm_cell_small,&
            cell_ref=qm_cell_small, use_ref_cell=.FALSE.,qmmm=.TRUE.,&
            globenv=globenv,&
            error=error)

       CALL cp_subsys_release(subsys_qm,error=error)
       CALL force_env_create(sub_force_env(2)%force_env,&
            qs_env=qs_env,globenv=globenv, error=error)

       !
       ! Initialize MM Potential fitted with Gaussian
       !
       CALL qmmm_init_gaussian_type (qmmm_env_qm = qmmm_env_qm,&
                                     qm_cell_small=qm_cell_small,&
                                     mm_cell=mm_cell,&
                                     globenv=globenv,&
                                     qs_env=qs_env,&
                                     added_charges=qmmm_env_qm%added_charges,&
                                     error=error)
       !
       ! Initialize the MM potential stored on vector
       !
       CALL qmmm_init_potential(qmmm_env_qm=qmmm_env_qm,&
                                qm_cell_small=qm_cell_small,&
                                mm_cell=mm_cell,&
                                globenv=globenv,&
                                qs_env=qs_env,&
                                added_charges=qmmm_env_qm%added_charges,&
                                error=error)
       !
       ! Initialize the qmmm_pw_grid
       !
       CALL qmmm_pw_grid_init(qmmm_env=qmmm_env_qm,&
                              pw_env=qs_env%pw_env,&
                              cell=qm_cell_small,&
                              error=error)

       CALL cell_release(qm_cell_small, error)                       
       CALL qs_env_release(qs_env,error=error)

       globenv%program_name="QMMM"
       CALL force_env_create(force_env, sub_force_env=sub_force_env,&
            qmmm_env=qmmm_env_qm,globenv=globenv,error=error)

       CALL qmmm_env_qm_release(qmmm_env_qm, error)
       CALL qmmm_env_mm_release(qmmm_env_mm, error)
       DEALLOCATE(qm_atom_type,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_create_force_env
!***************************************************************************

!!****f* qmmm_init/create_small_subsystem [1.0] *
!!
!!   NAME
!!     create_small_subsystem
!!
!!   FUNCTION
!!     updates the molecule information of the given subsys
!!
!!   NOTES
!!     not really ready to be used with different para_envs for the small
!!     and big part
!!
!!   INPUTS
!!     - small_subsys: the subsystem to create
!!     - big_subsys: the superset of small_subsys
!!     - small_para_env: the parallel environment for the new (small)
!!       subsys
!!     - sub_atom_index: indexes of the atoms that should be in small_subsys
!!     - sub_subsys_index: indexes of the subsys of big_subsys from which
!!       the corresponding atom should be taken
!!     - globenv: place where to find the input file
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE create_small_subsystem(small_subsystem,big_subsys,small_cell,&
       small_para_env,sub_atom_index,sub_subsys_index,sub_atom_kind_name,&
       globenv,error)
    TYPE(cp_subsystem_type), POINTER         :: small_subsystem
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: big_subsys
    TYPE(cell_type), POINTER                 :: small_cell
    TYPE(cp_para_env_type), POINTER          :: small_para_env
    INTEGER, DIMENSION(:), INTENT(in)        :: sub_atom_index, &
                                                sub_subsys_index
    CHARACTER(len=default_string_length), &
      DIMENSION(:), INTENT(in)               :: sub_atom_kind_name
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_small_subsystem', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, isubsys, nat, stat, id_
    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(mol_kind_new_list_type), POINTER    :: mol_kinds
    TYPE(mol_new_list_type), POINTER         :: mols
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(topology_parameters_type)           :: topology

    failure=.FALSE.
    NULLIFY(mol_kinds,mols,particles,atomic_kinds, atomic_kind_set, particle_set, &
         molecule_kind_set,molecule_set,particles,atomic_kinds)

    CPPrecondition(.NOT.ASSOCIATED(small_subsystem),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(big_subsys),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CPAssert(big_subsys(1)%subsys%para_env%group==small_para_env%group,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 1. Initialize the topology structure type
       !-----------------------------------------------------------------------------
       CALL init_topology(topology)


       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 2. Get the cell info
       !-----------------------------------------------------------------------------
       topology%cell=> small_cell
       CALL cell_retain(small_cell,error=error)

       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 3. Read in the topology section in the input file if any
       !-----------------------------------------------------------------------------

       ! initialize atom coords from the bigger system
       nat=SIZE(sub_atom_index)
       topology%natoms=nat
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%r),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%label_atmname),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%label_molname),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%label_resname),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%atm_mass),cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.ASSOCIATED(topology%atom_info%atm_charge),cp_failure_level,routineP,error,failure)
       ALLOCATE(topology%atom_info%r(3,nat),topology%atom_info%label_atmname(nat),&
            topology%atom_info%label_molname(nat),topology%atom_info%label_resname(nat),&
            topology%atom_info%element(nat),&
            topology%atom_info%atm_mass(nat), topology%atom_info%atm_charge(nat),stat=stat)
       isubsys=-1
       DO iat=1,nat
          IF (isubsys/=sub_subsys_index(iat)) THEN
             isubsys=sub_subsys_index(iat)
             CALL cp_subsys_get(big_subsys(isubsys)%subsys,&
                  particles=particles, error=error)
          END IF
          topology%atom_info%r(:,iat)=particles%els(sub_atom_index(iat))%r
          topology%atom_info%label_atmname(iat)=sub_atom_kind_name(iat)
          topology%atom_info%label_molname(iat)=topology%atom_info%label_atmname(iat)
          topology%atom_info%label_resname(iat)=topology%atom_info%label_atmname(iat)
          id_ = INDEX(topology%atom_info%label_atmname(iat),"_")-1
          IF ( id_ == -1 ) id_ = LEN_TRIM(topology%atom_info%label_atmname(iat))
          topology%atom_info%element(iat)=topology%atom_info%label_atmname(iat)(1:id_)
          topology%atom_info%atm_mass(iat)=0._dp
          topology%atom_info%atm_charge(iat)=0._dp
       END DO
       topology%bondparm_factor=0._dp


       !-----------------------------------------------------------------------------
       !-----------------------------------------------------------------------------
       ! 5. Read in or generate the molecular connectivity
       !-----------------------------------------------------------------------------
       CALL connectivity_control(topology,globenv)


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


       CALL cp_subsys_create(small_subsystem,para_env=small_para_env,error=error)
       CALL particle_list_create(particles,els_ptr=particle_set,error=error)
       CALL atomic_kind_list_create(atomic_kinds,els_ptr=atomic_kind_set,error=error)
       CALL mol_new_list_create(mols,els_ptr=molecule_set,error=error)
       CALL mol_kind_new_list_create(mol_kinds,els_ptr=molecule_kind_set,error=error)
       CALL cp_subsys_set(small_subsystem,particles=particles,atomic_kinds=atomic_kinds,&
            molecules_new=mols,molecule_kinds_new=mol_kinds,error=error)
       CALL particle_list_release(particles,error=error)
       CALL atomic_kind_list_release(atomic_kinds,error=error)
       CALL mol_new_list_release(mols,error=error)
       CALL mol_kind_new_list_release(mol_kinds,error=error)
    END IF
  END SUBROUTINE create_small_subsystem
!***************************************************************************

!!****f* qmmm_init/assign_mm_charges_and_radius [1.0] *
!!
!!   NAME
!!     assign_mm_charges_and_radius
!!
!!   FUNCTION
!!     Assigns charges and radius to evaluate the MM electrostatic potential
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     - subsys: the subsystem containing the MM charges
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*** **********************************************************************
 SUBROUTINE assign_mm_charges_and_radius(subsys, mm_atom_chrg,mm_el_pot_radius,&
                                         mm_atom_index, mm_link_atoms, mm_link_scale_factor,&
                                         error)
   IMPLICIT NONE
   ! Arguments
   TYPE(cp_subsystem_p_type), &
        DIMENSION(:), POINTER                  :: subsys
   REAL(dp), DIMENSION(:), POINTER             :: mm_atom_chrg, mm_el_pot_radius,mm_link_scale_factor
   INTEGER, DIMENSION(:), POINTER              :: mm_atom_index, mm_link_atoms
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error
   ! Local Variables
   TYPE(atomic_kind_type), POINTER             :: my_kind
    TYPE(elp_potential_type), POINTER          :: my_potential
   TYPE(particle_list_type), POINTER           :: particles
   TYPE(particle_type), DIMENSION(:), &
        POINTER                                :: particle_set
   INTEGER  :: I, IndMM, ilink, iw
   REAL(dp) :: qi, ri

   NULLIFY(particle_set, my_kind)
   CALL cp_subsys_get(subsys=subsys(1)%subsys,&
                      particles=particles)
   particle_set => particles%els

   DO I = 1, SIZE(mm_atom_index)
      IndMM = mm_atom_index(I)
      my_kind => particle_set(IndMM)%atomic_kind
      CALL get_atomic_kind(atomic_kind=my_kind, elp_potential=my_potential)
      CALL get_potential  (potential=my_potential,& 
                           qeff=qi,&
                           core_charge_radius=ri)
      mm_atom_chrg(I)     = qi
      mm_el_pot_radius(I) = ri
   END DO

   DO ilink = 1, SIZE(mm_link_atoms)
      DO i = 1, SIZE(mm_atom_index)
         IF (mm_atom_index(i) == mm_link_atoms(ilink)) EXIT
      END DO
      IndMM = mm_atom_index(I)
      mm_atom_chrg(i) = mm_atom_chrg(i) * mm_link_scale_factor(ilink)
   END DO

 END SUBROUTINE assign_mm_charges_and_radius

!!****f* qmmm_init/print_qmmm_charges [1.0] *
!!
!!   NAME
!!     print_qmmm_charges
!!
!!   FUNCTION
!!     Print info on charges generating the qmmm potential..
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     - subsys: the subsystem containing the MM charges
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2005 created [tlaino]
!!
!!*** **********************************************************************
 SUBROUTINE print_qmmm_charges(mm_atom_index, mm_atom_chrg, mm_el_pot_radius, added_charges,&
                               globenv, qmmm_section, error)
   IMPLICIT NONE
   ! Arguments
   REAL(dp), DIMENSION(:), POINTER             :: mm_atom_chrg, mm_el_pot_radius
   INTEGER, DIMENSION(:), POINTER              :: mm_atom_index
   TYPE(section_vals_type), POINTER            :: qmmm_section
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error
    TYPE(global_environment_type), POINTER     :: globenv
    TYPE(add_set_type), POINTER                :: added_charges
    ! Local Variables
    INTEGER                                    :: IndMM, I, iw, ind1, ind2
    REAL(KIND=dp)                              :: qi, ri, qtot
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(section_vals_type), POINTER :: print_key

    qtot = 0.0_dp
    iw = globenv%scr
    logger => cp_error_get_logger(error)
    print_key => section_vals_get_subs_vals(qmmm_section,"qmmm_charges",error=error)
    IF (cp_should_output(logger%iter_info,print_key,error=error).AND.&
         logger%para_env%mepos==logger%para_env%source) THEN 
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
       WRITE(iw,FMT='(/5X,A)')"MM    POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL"
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
       DO I = 1, SIZE(mm_atom_index)
          IndMM = mm_atom_index(I)
          qi = mm_atom_chrg(I)    
          qtot = qtot + qi 
          ri = mm_el_pot_radius(I) 
          WRITE(iw,'(7X,A9,I5,A8,F12.6,A8,F12.6)')' MM ATOM:',IndMM,' RADIUS:',ri,&
                                                  ' CHARGE:',qi
       END DO
       IF (added_charges%num_mm_atoms /= 0) THEN
          WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
          WRITE(iw,'(/5X,A)')"ADDED POINT CHARGES GENERATING THE QM/MM ELECTROSTATIC POTENTIAL"
          WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
          DO I = 1, SIZE(added_charges%mm_atom_index)
             IndMM = added_charges%mm_atom_index(I)
             qi = added_charges%mm_atom_chrg(I)     
             qtot = qtot + qi
             ri = added_charges%mm_el_pot_radius(I) 
             ind1 = added_charges%add_env(I)%Index1
             ind2 = added_charges%add_env(I)%Index2
             WRITE(iw,'(7X,A9,I5,A8,F12.6,A8,F12.6,I5,I5)')'MM POINT:',IndMM,' RADIUS:',ri,&
                                                           ' CHARGE:',qi,ind1,ind2
          END DO

       END IF
       WRITE(iw,FMT="(/,T2,A)") REPEAT("-",73)
       WRITE(iw,'(/,T36,A,F12.6)')' TOTAL CHARGE:',qtot
       WRITE(iw,FMT="(/,T2,A,/)") REPEAT("-",73)
    END IF
  END SUBROUTINE print_qmmm_charges

!!****f* qmmm_init/print_qmmm_links [1.0] *
!!
!!   NAME
!!     print_qmmm_links
!!
!!   FUNCTION
!!     Print info on qm/mm links
!!
!!   NOTES
!!     
!!
!!   INPUTS
!!     - subsys: the subsystem containing the MM charges
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     01.2005 created [tlaino]
!!
!!*** **********************************************************************
 SUBROUTINE print_qmmm_links(globenv, qmmm_section, qmmm_links, error)
   IMPLICIT NONE
   ! Arguments
   TYPE(section_vals_type), POINTER            :: qmmm_section
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error
    TYPE(global_environment_type), POINTER     :: globenv
    TYPE(qmmm_links_type), POINTER             :: qmmm_links
    ! Local Variables
    INTEGER                                    :: mm_index, qm_index, iw, i
    REAL(KIND=dp)                              :: alpha
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(section_vals_type), POINTER :: print_key

    iw = globenv%scr
    logger => cp_error_get_logger(error)
    print_key => section_vals_get_subs_vals(qmmm_section,"qmmm_link_info",error=error)
    IF (cp_should_output(logger%iter_info,print_key,error=error).AND.&
         logger%para_env%mepos==logger%para_env%source) THEN 
       WRITE(iw,FMT="(/,T2, A)") REPEAT("-",73)
       WRITE(iw,FMT="(/,T31,A)") " QM/MM LINKS "
       WRITE(iw,FMT="(/,T2, A)") REPEAT("-",73)
       IF (ASSOCIATED(qmmm_links)) THEN
          DO I = 1, SIZE(qmmm_links%imomm)
             qm_index = qmmm_links%imomm(I)%link%qm_index
             mm_index = qmmm_links%imomm(I)%link%mm_index
             alpha    = qmmm_links%imomm(I)%link%alpha
             WRITE(iw,FMT="(T2,A,T20,A9,I6,1X,A9,I6,T55,A6,F12.6)")"TYPE: IMOMM",&
                  "QM INDEX:",qm_index,"MM INDEX:",mm_index,"ALPHA:",alpha
          END DO
       END IF
       WRITE(iw,FMT="(/,T2,A,/)") REPEAT("-",73)
    END IF
  END SUBROUTINE print_qmmm_links

  !!****f* qmmm_init/qmmm_init_gaussian_type [1.0] *
  !!
  !!   NAME
  !!     qmmm_init_gaussian_type
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     1.2005 created [tlaino]
  !!
  !!*************************************************************************
 SUBROUTINE qmmm_init_gaussian_type(qmmm_env_qm, qm_cell_small, mm_cell, globenv,&
                                    qs_env, added_charges, error )
   IMPLICIT NONE
   ! Arguments
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(add_set_type), POINTER              :: added_charges
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
   ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_init_gaussian_type',&
         routineP = moduleN//':'//routineN
    REAL(KIND=dp), DIMENSION(:), POINTER :: maxradius, maxradius2
    INTEGER :: i, stat
    LOGICAL :: failure

    failure=.FALSE.
    stat=0
    NULLIFY(maxradius, maxradius2)

    CALL qmmm_gaussian_initialize(qmmm_gaussian_fns=qmmm_env_qm%pgfs,&
                                  qm_cell_small=qm_cell_small,&
                                  mm_cell=mm_cell,&
                                  globenv=globenv,&
                                  pw_env=qs_env%pw_env,&
                                  mm_el_pot_radius=qmmm_env_qm%mm_el_pot_radius,&
                                  qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                                  eps_mm_rspace=qmmm_env_qm%eps_mm_rspace,&
                                  maxradius=maxradius,&
                                  error=error)

    IF (qmmm_env_qm%move_mm_charges.OR.qmmm_env_qm%add_mm_charges) THEN
       CALL qmmm_gaussian_initialize(qmmm_gaussian_fns=added_charges%pgfs,&
                                     qm_cell_small=qm_cell_small,&
                                     mm_cell=mm_cell,&
                                     globenv=globenv,&
                                     pw_env=qs_env%pw_env,&
                                     mm_el_pot_radius=added_charges%mm_el_pot_radius,&
                                     qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                                     eps_mm_rspace=qmmm_env_qm%eps_mm_rspace,&
                                     maxradius=maxradius2,&
                                     error=error)

       IF (qmmm_env_qm%qmmm_coupl_type == do_qmmm_gauss) THEN
          DO i = 1, SIZE(maxradius)
             maxradius(i) = MAX( maxradius(i), maxradius2(i))
          END DO
       END IF
       IF (ASSOCIATED(maxradius2)) DEALLOCATE(maxradius2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    qmmm_env_qm%maxradius => maxradius

 END SUBROUTINE qmmm_init_gaussian_type

  !!****f* qmmm_init/qmmm_init_potential [1.0] *
  !!
  !!   NAME
  !!     qmmm_init_potential
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     1.2005 created [tlaino]
  !!
  !!*************************************************************************
 SUBROUTINE qmmm_init_potential(qmmm_env_qm, qm_cell_small, mm_cell, globenv, qs_env,&
                                added_charges, error)
   IMPLICIT NONE
   ! Arguments
   TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small
   TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
   TYPE(global_environment_type), POINTER   :: globenv
   TYPE(qs_environment_type), POINTER       :: qs_env
   TYPE(add_set_type), POINTER              :: added_charges
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error   
   ! Local Variables
   
   CALL qmmm_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                            mm_el_pot_radius=qmmm_env_qm%mm_el_pot_radius,&
                            potentials=qmmm_env_qm%potentials,&
                            pgfs=qmmm_env_qm%pgfs,&
                            qm_cell_small=qm_cell_small,&
                            mm_cell=mm_cell,&
                            globenv=globenv,&
                            pw_env=qs_env%pw_env,&
                            error=error)

   IF (qmmm_env_qm%move_mm_charges.OR.qmmm_env_qm%add_mm_charges) THEN

      CALL qmmm_potential_init(qmmm_coupl_type=qmmm_env_qm%qmmm_coupl_type,&
                               mm_el_pot_radius=added_charges%mm_el_pot_radius,&
                               potentials=added_charges%potentials,&
                               pgfs=added_charges%pgfs,&
                               qm_cell_small=qm_cell_small,&
                               mm_cell=mm_cell,&
                               globenv=globenv,&
                               pw_env=qs_env%pw_env,&
                               error=error)      
   END IF

 END SUBROUTINE qmmm_init_potential

  !!****f* qmmm_init/setup_qmmm_vars_qm [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_vars_qm
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     11.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE setup_qmmm_vars_qm( qmmm_section, qmmm_env, subsys_mm, qm_atom_type,&
                              qm_atom_index, mm_atom_index, qm_subsys_index,&
                              qm_cell_small, qmmm_coupl_type, eps_mm_rspace, qmmm_link,&
                              error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_env_qm_type), POINTER             :: qmmm_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_mm
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: qm_atom_type
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_atom_index, qm_subsys_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    INTEGER, INTENT(OUT)                     :: qmmm_coupl_type
    REAL(KIND=dp), INTENT(OUT)               :: eps_mm_rspace
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(OUT)                     :: qmmm_link
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_vars_qm',&
         routineP = moduleN//':'//routineN
    INTEGER :: i, j, icount, size_mm_system, stat, ikind, ikindr, nkind, n_rep, iatm
    LOGICAL :: is_mm
    INTEGER                                  :: qm_subsys_type
    LOGICAL                                  :: failure
    CHARACTER(len=default_string_length)     :: qm_atom_kind, mm_atom_kind
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type),POINTER           :: atomic_kind
    CHARACTER(LEN=default_string_length)     :: atmname
    REAL(KIND=dp)                                :: tmp_radius
    TYPE(elp_potential_type), POINTER        :: elp_potential
    TYPE(section_vals_type), POINTER         :: mm_kinds
    LOGICAL :: explicit
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cell_par
    INTEGER, DIMENSION(:), POINTER           :: mm_link_atoms

    NULLIFY(mm_link_atoms)
    qmmm_link = .FALSE.
    failure   = .FALSE.
    CALL section_vals_get(qmmm_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(qmmm_section,"E_COUPL",i_val=qmmm_coupl_type,error=error)
       CALL section_vals_val_get(qmmm_section,"EPS_MM_RSPACE",r_val=eps_mm_rspace,error=error)
       !
       ! Initialization of arrays and core_charge_radius...
       !
       tmp_radius=0.0_dp
       CALL cp_subsys_get(subsys=subsys_mm(1)%subsys,&
                          atomic_kinds=atomic_kinds,&
                          error=error)
       DO Ikind = 1, SIZE(atomic_kinds%els)
          atomic_kind => atomic_kinds%els(Ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               elp_potential=elp_potential)
          CALL set_potential  (potential=elp_potential,&
                               core_charge_radius=tmp_radius)
          CALL set_atomic_kind(atomic_kind=atomic_kind,&
                               elp_potential=elp_potential)
       END DO
       CALL setup_qm_atom_list(qmmm_section=qmmm_section,&
                               qm_atom_index=qm_atom_index,&
                               qm_atom_type=qm_atom_type,&
                               mm_link_atoms=mm_link_atoms,&
                               qm_subsys_index=qm_subsys_index,&
                               qmmm_link=qmmm_link,&
                               error=error)
       !
       ! MM_KINDS
       !
       mm_kinds => section_vals_get_subs_vals(qmmm_section,"MM_KIND",error=error)
       CALL section_vals_get(mm_kinds,explicit=explicit,n_repetition=nkind,error=error)
       !
       ! Default
       !
       tmp_radius = RADIUS_QMMM_DEFAULT * qm_cell_small%unit_of_length
       Set_Radius_Pot_0: DO IkindR = 1, SIZE(atomic_kinds%els)
          atomic_kind => atomic_kinds%els(IkindR)
          CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
               elp_potential=elp_potential)
          CALL set_potential  (potential=elp_potential,&
               core_charge_radius=tmp_radius)
          CALL set_atomic_kind(atomic_kind=atomic_kind,&
               elp_potential=elp_potential)
       END DO Set_Radius_Pot_0
       !
       ! If present overwrite the kind specified in input file...
       !
       IF (explicit) THEN
          DO ikind=1,nkind
             CALL section_vals_val_get(mm_kinds,"_SECTION_PARAMETERS_",i_rep_section=ikind,&
                  c_val=mm_atom_kind,error=error)
             CALL section_vals_val_get(mm_kinds,"RADIUS",i_rep_section=ikind,r_val=tmp_radius,error=error)
             tmp_radius = tmp_radius * qm_cell_small%unit_of_length
             Set_Radius_Pot_1: DO IkindR = 1, SIZE(atomic_kinds%els)
                atomic_kind => atomic_kinds%els(IkindR)
                CALL get_atomic_kind(atomic_kind=atomic_kind,name=atmname)
                IF (TRIM(mm_atom_kind)==atmname) THEN
                   CALL get_atomic_kind(atomic_kind=atomic_kind,&
                        elp_potential=elp_potential)
                   CALL set_potential  (potential=elp_potential,&
                        core_charge_radius=tmp_radius)
                   CALL set_atomic_kind(atomic_kind=atomic_kind,&
                        elp_potential=elp_potential)
                END IF
             END DO Set_Radius_Pot_1
          END DO
       END IF
       !
       ! QM Cell Info
       !
       CALL section_vals_val_get(qmmm_section,"A_QM",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
           CALL section_vals_val_get(qmmm_section,"A_QM",r_vals=cell_par,error=error)
          qm_cell_small%hmat(:,1) = cell_par(:)
       END IF
       CALL section_vals_val_get(qmmm_section,"B_QM",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(qmmm_section,"B_QM",r_vals=cell_par,error=error)
          qm_cell_small%hmat(:,2) = cell_par(:)
       END IF
       CALL section_vals_val_get(qmmm_section,"C_QM",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(qmmm_section,"C_QM",r_vals=cell_par,error=error)
          qm_cell_small%hmat(:,3) = cell_par(:)
       END IF
       CALL section_vals_val_get(qmmm_section,"ABC_QM",n_rep_val=n_rep,error=error)
       IF (n_rep /= 0) THEN
          CALL section_vals_val_get(qmmm_section,"ABC_QM",r_vals=cell_par,error=error)
          qm_cell_small%hmat(1,1) = cell_par(1)
          qm_cell_small%hmat(2,2) = cell_par(2)
          qm_cell_small%hmat(3,3) = cell_par(3)
       END IF
       qm_cell_small%hmat=qm_cell_small%hmat(:,:)*qm_cell_small%unit_of_length
    ELSE
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "QMMM section not present in input file!"//&
          CPSourceFileRef,&
          error,failure)
    ENDIF
    !
    ! Initialize cell parameters..
    !
    CALL init_qmmm_cell(qm_cell_small)
    !
    ! Build MM atoms list
    !
    size_mm_system = SIZE(subsys_mm(1)%subsys%particles%els) - SIZE(qm_atom_index)
    IF (qmmm_link) size_mm_system = size_mm_system + SIZE(mm_link_atoms)
    ALLOCATE(mm_atom_index(size_mm_system))
    icount = 0
    DO i = 1, SIZE(subsys_mm(1)%subsys%particles%els)
       is_mm = .TRUE.
       IF ( ANY(qm_atom_index == i)) is_mm =.FALSE.
       IF (ANY( mm_link_atoms == i ).AND.qmmm_link) is_mm =.TRUE.
       IF (is_mm) THEN
          icount = icount + 1
          mm_atom_index(icount) = i
       END IF
    END DO
    CPPostcondition(icount==size_mm_system,cp_failure_level,routineP,error,failure)
    IF (ASSOCIATED(mm_link_atoms)) THEN 
       DEALLOCATE(mm_link_atoms, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE setup_qmmm_vars_qm


  !!****f* qmmm_init/setup_qmmm_vars_mm [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_vars_mm
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE setup_qmmm_vars_mm( qmmm_section, qmmm_env, qm_atom_index,  &
                                 mm_link_atoms, mm_link_scale_factor, &
                                 fist_scale_charge_link, qmmm_coupl_type,&
                                 qmmm_link, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(qmmm_env_mm_type), POINTER          :: qmmm_env
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, mm_link_atoms
    INTEGER, INTENT(OUT)                     :: qmmm_coupl_type
    LOGICAL, INTENT(OUT)                     :: qmmm_link
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_link_scale_factor,fist_scale_charge_link
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_vars_mm',&
         routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure
    LOGICAL :: explicit
    
    qmmm_link = .FALSE.
    CALL section_vals_get(qmmm_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(qmmm_section,"E_COUPL",i_val=qmmm_coupl_type,error=error)
       CALL setup_qm_atom_list(qmmm_section,qm_atom_index=qm_atom_index,qmmm_link=qmmm_link,&
            mm_link_atoms=mm_link_atoms, mm_link_scale_factor = mm_link_scale_factor,&
            fist_scale_charge_link=fist_scale_charge_link, error=error)
    END IF
  END SUBROUTINE setup_qmmm_vars_mm
     
  !!****f* qmmm_init/setup_qm_atom_list [1.0] *
  !!
  !!   NAME
  !!     setup_qm_atom_list
  !!
  !!   FUNCTION
  !!
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE setup_qm_atom_list(qmmm_section,qm_atom_index,qm_atom_type,&
       qm_subsys_index, mm_link_atoms, mm_link_scale_factor, qmmm_link, &
       fist_scale_charge_link, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER, OPTIONAL        :: qm_atom_type
    INTEGER, DIMENSION(:), POINTER, OPTIONAL :: qm_atom_index, qm_subsys_index
    INTEGER, DIMENSION(:), POINTER, OPTIONAL :: mm_link_atoms
    REAL(KIND=dp), DIMENSION(:), POINTER, OPTIONAL     :: mm_link_scale_factor
    REAL(KIND=dp), DIMENSION(:), POINTER, OPTIONAL     :: fist_scale_charge_link
    LOGICAL, OPTIONAL, INTENT(OUT) :: qmmm_link
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qm_atom_list',&
         routineP = moduleN//':'//routineN
    LOGICAL :: explicit, failure
    INTEGER :: num_qm_atom_tot, nkind, ikind, nlinks, stat, mm_index
    INTEGER :: qm_subsys_type
    INTEGER, DIMENSION(:), POINTER           :: mm_indexes
    TYPE(section_vals_type), POINTER         :: qm_kinds, qmmm_links
    CHARACTER(len=default_string_length)     :: qm_atom_kind
    REAL(KIND=dp) :: scale_f

    failure         = .FALSE.
    num_qm_atom_tot = 0
    nlinks          = 0
    stat            = 0
    !
    ! QM_KINDS
    !
    qm_kinds => section_vals_get_subs_vals(qmmm_section,"QM_KIND",error=error)
    CALL section_vals_get(qm_kinds,n_repetition=nkind,error=error)
    DO ikind=1,nkind
       CALL section_vals_val_get(qm_kinds,"MM_INDEX",i_rep_section=ikind,i_vals=mm_indexes,&
            error=error)
       num_qm_atom_tot = num_qm_atom_tot + SIZE(mm_indexes)             
    END DO
    !
    ! QM/MM LINKS
    !
    qmmm_links => section_vals_get_subs_vals(qmmm_section,"LINK",error=error)
    CALL section_vals_get(qmmm_links,explicit=explicit,error=error)
    IF (explicit) THEN
       qmmm_link = .TRUE.
       CALL section_vals_get(qmmm_links,n_repetition=nlinks,error=error)
       num_qm_atom_tot = num_qm_atom_tot + nlinks
    END IF
    IF (PRESENT(mm_link_scale_factor)) ALLOCATE(mm_link_scale_factor(nlinks), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(fist_scale_charge_link)) ALLOCATE(fist_scale_charge_link(nlinks), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(mm_link_atoms))   ALLOCATE(mm_link_atoms(nlinks),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(qm_atom_index))   ALLOCATE(qm_atom_index(num_qm_atom_tot),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(qm_atom_type))    ALLOCATE(qm_atom_type(num_qm_atom_tot),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(qm_subsys_index)) ALLOCATE(qm_subsys_index(num_qm_atom_tot),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(qm_atom_index))   qm_atom_index    =   0
    IF (PRESENT(qm_atom_type))    qm_atom_type     = " "
    IF (PRESENT(qm_subsys_index)) qm_subsys_index  =   0
    num_qm_atom_tot = 1
    DO ikind=1,nkind
       CALL section_vals_val_get(qm_kinds,"MM_INDEX",i_rep_section=ikind,i_vals=mm_indexes,&
            error=error)
       IF (PRESENT(qm_atom_index)) THEN
          qm_atom_index(num_qm_atom_tot:num_qm_atom_tot + SIZE(mm_indexes)-1)   = mm_indexes(:)
       END IF
       IF (PRESENT(qm_atom_type)) THEN
          CALL section_vals_val_get(qm_kinds,"_SECTION_PARAMETERS_",i_rep_section=ikind,&
               c_val=qm_atom_kind,error=error)
          qm_atom_type (num_qm_atom_tot:num_qm_atom_tot + SIZE(mm_indexes)-1)   = qm_atom_kind
       END IF
       IF (PRESENT(qm_subsys_index)) THEN
          CALL section_vals_val_get(qm_kinds,"SUBSYS",i_rep_section=ikind,i_val=qm_subsys_type,&
               error=error)
          qm_subsys_index(num_qm_atom_tot:num_qm_atom_tot + SIZE(mm_indexes)-1) = qm_subsys_type
       END IF
       num_qm_atom_tot = num_qm_atom_tot + SIZE(mm_indexes)
    END DO
    IF (PRESENT(mm_link_scale_factor))     mm_link_scale_factor = 0.0_dp
    IF (PRESENT(fist_scale_charge_link)) fist_scale_charge_link = 0.0_dp
    IF (PRESENT(mm_link_atoms))            mm_link_atoms        = 0
    IF (explicit) THEN
       DO ikind = 1, nlinks
          IF (PRESENT(qm_subsys_index)) THEN
             CALL section_vals_val_get(qmmm_links,"SUBSYS",i_rep_section=ikind,&
                  i_val=qm_subsys_type,error=error)
             qm_subsys_index(num_qm_atom_tot:num_qm_atom_tot) = qm_subsys_type
          END IF
          IF (PRESENT(qm_atom_type)) THEN
             qm_atom_type (num_qm_atom_tot:num_qm_atom_tot)   = "H_LINK"
          END IF
          IF (PRESENT(qm_atom_index)) THEN
             CALL section_vals_val_get(qmmm_links,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,&
                  error=error)
             qm_atom_index(num_qm_atom_tot:num_qm_atom_tot )  = mm_index
             num_qm_atom_tot = num_qm_atom_tot + 1
          END IF
          IF (PRESENT(mm_link_atoms)) THEN
             CALL section_vals_val_get(qmmm_links,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,&
                  error=error)
             mm_link_atoms (ikind) = mm_index
          END IF
          IF (PRESENT(mm_link_scale_factor)) THEN
             CALL section_vals_val_get(qmmm_links,"QMMM_SCALE_FACTOR",i_rep_section=ikind,r_val=scale_f,&
                  error=error)
             mm_link_scale_factor(ikind) = scale_f
          END IF
          IF (PRESENT(fist_scale_charge_link)) THEN
             CALL section_vals_val_get(qmmm_links,"FIST_SCALE_FACTOR",i_rep_section=ikind,r_val=scale_f,&
                  error=error)
             fist_scale_charge_link(ikind) = scale_f
          END IF
       END DO
    END IF
    CPPostcondition(num_qm_atom_tot-1==SIZE(qm_atom_index),cp_failure_level,routineP,error,failure)

  END SUBROUTINE setup_qm_atom_list


  !!****f* qmmm_init/setup_qmmm_links [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_links
  !!
  !!   FUNCTION
  !!     this routine sets up all variables to treat qmmm links
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE  setup_qmmm_links(qmmm_section, qmmm_links, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section, qmmm_link_section
    TYPE(qmmm_links_type), POINTER           :: qmmm_links
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local variables
    CHARACTER(len=*), PARAMETER :: routineN = 'setup_qmmm_links',&
         routineP = moduleN//':'//routineN
    INTEGER  :: ikind, nlinks, n_imomm, n_gho, link_type, qm_index, mm_index, stat
    LOGICAL  :: failure, charge_scale
    REAL(KIND=dp) :: alpha, scale_factor

    failure = .FALSE.
    n_imomm = 0
    n_gho   = 0
    qmmm_link_section => section_vals_get_subs_vals(qmmm_section,"LINK",error=error)
    CALL section_vals_get(qmmm_link_section,n_repetition=nlinks,error=error)    
    CPPostcondition(nlinks /= 0,cp_failure_level,routineP,error,failure)
    DO ikind= 1, nlinks
       CALL section_vals_val_get(qmmm_link_section,"LINK_TYPE",i_rep_section=ikind,i_val=link_type,error=error)
       IF (link_type ==  do_qmmm_link_imomm) n_imomm = n_imomm + 1
       IF (link_type ==  do_qmmm_link_gho  ) n_gho   = n_gho   + 1
    END DO
    CPPostcondition(n_imomm + n_gho /= 0,cp_failure_level,routineP,error,failure)
    ALLOCATE(qmmm_links, stat=stat)
    CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
    IF ( n_imomm /= 0 ) THEN
       ALLOCATE(qmmm_links%imomm(n_imomm), stat=stat)
       CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
       DO ikind = 1, n_imomm
          NULLIFY(qmmm_links%imomm(ikind)%link)
          ALLOCATE(qmmm_links%imomm(ikind)%link,stat=stat)
          CPPostcondition(stat == 0,cp_failure_level,routineP,error,failure)
       END DO
       n_imomm = 0
       DO ikind = 1, nlinks
          CALL section_vals_val_get(qmmm_link_section,"LINK_TYPE",i_rep_section=ikind,i_val=link_type,error=error)
          IF (link_type ==  do_qmmm_link_imomm) THEN
             n_imomm = n_imomm + 1
             CALL section_vals_val_get(qmmm_link_section,"QM_INDEX",i_rep_section=ikind,i_val=qm_index,error=error)
             CALL section_vals_val_get(qmmm_link_section,"MM_INDEX",i_rep_section=ikind,i_val=mm_index,error=error)
             CALL section_vals_val_get(qmmm_link_section,"ALPHA_IMOMM",i_rep_section=ikind,r_val=alpha,error=error)
             qmmm_links%imomm(n_imomm)%link%qm_index = qm_index
             qmmm_links%imomm(n_imomm)%link%mm_index = mm_index
             qmmm_links%imomm(n_imomm)%link%alpha    = alpha
          END IF
       END DO
    END IF
    IF ( n_gho /= 0) THEN
       ! not yet implemented
       ! still to define : type, implementation into QS
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure) 
    END IF
  END SUBROUTINE setup_qmmm_links


  !!****f* qmmm_init/setup_qmmm_links [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_links
  !!
  !!   FUNCTION
  !!     this routine sets up all variables to treat qmmm links
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
  SUBROUTINE move_or_add_atoms(qmmm_section, move_mm_charges, add_mm_charges, &
                               mm_atom_chrg, mm_el_pot_radius, added_charges, &
                               mm_atom_index, qm_cell_small, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(section_vals_type), POINTER         :: qmmm_section
    LOGICAL, INTENT(OUT)                     :: move_mm_charges, add_mm_charges
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_atom_chrg, mm_el_pot_radius
    TYPE(add_set_type), POINTER              :: added_charges
    INTEGER, POINTER, DIMENSION(:)           :: mm_atom_index
    TYPE(cell_type), POINTER                 :: qm_cell_small
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                            :: error
    ! Local variables
    CHARACTER(len=*), PARAMETER :: routineN = 'move_or_add_atoms',&
         routineP = moduleN//':'//routineN
    LOGICAL :: failure, explicit
    INTEGER :: nlinks, n_moves, n_adds, icount, ikind, n_move_tot, n_add_tot
    INTEGER :: Index1, Index2, i_add, ind1
    TYPE(section_vals_type), POINTER         :: qmmm_link_section, move_section,&
                                                add_section
    REAL(KIND=dp) :: radius, alpha, charge

    failure = .FALSE.
    explicit = .FALSE.
    move_mm_charges = .FALSE.
    add_mm_charges  = .FALSE.
    NULLIFY( qmmm_link_section, move_section, add_section)
    IF (.NOT.failure) THEN
       qmmm_link_section => section_vals_get_subs_vals(qmmm_section,"LINK",error=error)
       CALL section_vals_get(qmmm_link_section,n_repetition=nlinks,error=error)
       CPPostcondition(nlinks /= 0,cp_failure_level,routineP,error,failure)
       icount     = 0
       n_move_tot = 0
       n_add_tot  = 0
       DO ikind= 1, nlinks
          move_section => section_vals_get_subs_vals(qmmm_link_section,"MOVE_MM_CHARGE",&
               i_rep_section=ikind, error=error)
          CALL section_vals_get(move_section,n_repetition=n_moves,error=error)
          add_section => section_vals_get_subs_vals(qmmm_link_section,"ADD_MM_CHARGE",&
               i_rep_section=ikind, error=error)
          CALL section_vals_get(add_section,n_repetition=n_adds,error=error)
          n_move_tot = n_move_tot + n_moves
          n_add_tot  = n_add_tot  + n_adds
       END DO
       icount = n_move_tot + n_add_tot
       IF (n_add_tot   /= 0)  add_mm_charges  = .TRUE.
       IF (n_move_tot  /= 0) move_mm_charges  = .TRUE.
       !
       ! create add_set_type
       !
       CALL create_add_set_type(added_charges, ndim=icount, error=error)
       !
       ! Fill in structures
       !
       icount = 0
       DO ikind= 1, nlinks
          move_section => section_vals_get_subs_vals(qmmm_link_section,"MOVE_MM_CHARGE",&
                                                     i_rep_section=ikind, error=error)
          CALL section_vals_get(move_section, explicit=explicit, n_repetition=n_moves, error=error)
          !
          ! Moving charge atoms
          !
          IF (explicit) THEN
             DO i_add = 1, n_moves
                icount = icount + 1
                CALL section_vals_val_get(move_section,"ATOM_INDEX_1",i_val=Index1,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(move_section,"ATOM_INDEX_2",i_val=Index2,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(move_section,"ALPHA",r_val=alpha,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(move_section,"RADIUS",r_val=radius,i_rep_section=i_add,error=error)
                
                radius = radius * qm_cell_small%unit_of_length
                CALL set_add_set_type(added_charges, icount, Index1, Index2, alpha, radius,&
                                      mm_atom_chrg=mm_atom_chrg, mm_el_pot_radius=mm_el_pot_radius,&
                                      mm_atom_index=mm_atom_index, move=n_moves, Ind1=ind1,&
                                      error=error)
             END DO
             mm_atom_chrg(ind1) = 0.0_dp
          END IF
          
          add_section => section_vals_get_subs_vals(qmmm_link_section,"ADD_MM_CHARGE",&
                                                    i_rep_section=ikind, error=error)
          CALL section_vals_get(add_section, explicit=explicit, n_repetition=n_adds, error=error)
          !
          ! Adding charge atoms
          !
          IF (explicit) THEN
             DO i_add = 1, n_adds
                icount = icount + 1
                CALL section_vals_val_get(add_section,"ATOM_INDEX_1",i_val=Index1,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(add_section,"ATOM_INDEX_2",i_val=Index2,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(add_section,"ALPHA",r_val=alpha,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(add_section,"RADIUS",r_val=radius,i_rep_section=i_add,error=error)
                CALL section_vals_val_get(add_section,"CHARGE",r_val=charge,i_rep_section=i_add,error=error)

                radius = radius * qm_cell_small%unit_of_length
                CALL set_add_set_type(added_charges, icount, Index1, Index2, alpha, radius, charge,&
                                      mm_atom_chrg=mm_atom_chrg, mm_el_pot_radius=mm_el_pot_radius,&
                                      mm_atom_index=mm_atom_index, error=error)
             END DO
          END IF
       END DO
    END IF

   END SUBROUTINE move_or_add_atoms

  !!****f* qmmm_init/setup_qmmm_links [1.0] *
  !!
  !!   NAME
  !!     setup_qmmm_links
  !!
  !!   FUNCTION
  !!     this routine sets up all variables of the add_set_type type
  !!
  !!   NOTES
  !!
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     12.2004 created [tlaino]
  !!
  !!*************************************************************************
   SUBROUTINE set_add_set_type(added_charges, icount, Index1, Index2, alpha, radius, charge,&
        mm_atom_chrg, mm_el_pot_radius, mm_atom_index, move, ind1, error)
     IMPLICIT NONE
     ! Arguments
     REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_atom_chrg, mm_el_pot_radius
     TYPE(add_set_type), POINTER              :: added_charges
     INTEGER, POINTER, DIMENSION(:)           :: mm_atom_index
     INTEGER, INTENT(IN)  :: Index1, Index2, icount
     INTEGER, INTENT(OUT),OPTIONAL  :: ind1
     REAL(KIND=dp), INTENT(IN) :: alpha, radius
     REAL(KIND=dp), INTENT(IN), OPTIONAL :: charge
     INTEGER, OPTIONAL, INTENT(in) :: move
     TYPE(cp_error_type), INTENT(inout), &
          OPTIONAL                            :: error
     ! Local Variables
     CHARACTER(len=*), PARAMETER :: routineN = 'set_add_set_type',&
          routineP = moduleN//':'//routineN
     LOGICAL :: failure
     INTEGER :: i, my_move
     REAL(KIND=dp) :: my_charge, my_radius
     
     failure = .FALSE.
     my_move = 0
     my_radius = radius
     IF (PRESENT(charge)) my_charge = charge
     IF (PRESENT(move))   my_move   = move
     IF (.NOT.failure) THEN
        i = 1
        GetId: DO WHILE (i<=SIZE(mm_atom_index))
           IF (Index1 == mm_atom_index(i)) EXIT GetId
           i = i+1
        END DO GetId
        IF (PRESENT(ind1)) ind1 = i
        CPPostcondition(i<=SIZE(mm_atom_index),cp_failure_level,routineP,error,failure)
        IF (.NOT.PRESENT(charge)) my_charge = mm_atom_chrg(i)/REAL(my_move,KIND=dp)
        IF (my_radius == 0.0_dp) my_radius = mm_el_pot_radius(i)

        added_charges%add_env(icount)%Index1   = Index1
        added_charges%add_env(icount)%Index2   = Index2
        added_charges%add_env(icount)%alpha    = alpha 
        added_charges%mm_atom_index(icount)    = icount
        added_charges%mm_atom_chrg(icount)     = my_charge
        added_charges%mm_el_pot_radius(icount) = my_radius
     END IF
   END SUBROUTINE set_add_set_type

END MODULE qmmm_init
