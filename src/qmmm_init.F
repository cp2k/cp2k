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
                                             read_atomic_kind_set, &
                                             get_atomic_kind
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
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
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dbl
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
  USE qmmm_types,                      ONLY: qmmm_env_create,&
                                             qmmm_env_release,&
                                             qmmm_env_type
  USE qs_environment,                  ONLY: qs_init_subsys
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_env_release,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_retain,&
                                             cell_type,&
                                             init_qmmm_cell,&
                                             cell_clone,&
                                             cell_release
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: connectivity_control,&
                                             force_field_control,&
                                             coordinate_control
  USE topology_pdb,                    ONLY: write_coordinate_pdb
  USE topology_psf,                    ONLY: write_topology_psf
  USE topology_types,                  ONLY: deallocate_topology,&
                                             init_topology,&
                                             topology_parameters_type
  USE topology_util,                   ONLY: topology_connectivity_pack,&
                                             topology_constraint_pack,&
                                             topology_coordinate_pack
  USE qs_parser,                       ONLY: finish_parser,&
                                             open_file,&
                                             read_object,&
                                             search_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE external_potential_types,        ONLY: elp_potential_type,&
                                             get_potential
  USE fist_main,                       ONLY: fist_create_force_env
  USE qmmm_input,                      ONLY: qmmm_namelist_read
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type
  USE qmmm_gaussian_init,              ONLY: qmmm_gaussian_initialize
  USE qmmm_pw_grid,                    ONLY: qmmm_pw_grid_init
  USE qmmm_elpot,                      ONLY: qmmm_potential_init

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
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: mm_cell, qm_cell_small
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_mm
    TYPE(cp_subsystem_type), POINTER         :: subsys_qm
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: sub_force_env
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(dp), DIMENSION(:), POINTER          :: mm_atom_chrg, mm_el_pot_radius
    REAL(KIND=dp)                            :: eps_mm_rspace
    CHARACTER(len=default_string_length)     :: qmmm_coupl_type

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    NULLIFY(sub_force_env, qm_atom_index, mm_atom_index, qm_atom_type,qmmm_env, subsys_mm,&
         subsys_qm, mm_cell, qm_cell_small, qs_env, para_env, mm_atom_chrg,&
         mm_el_pot_radius )
    para_env => globenv%para_env
    iw = globenv%scr
    IF (.NOT. failure) THEN
       ALLOCATE(sub_force_env(2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO isubf=1,2
          NULLIFY(sub_force_env(isubf)%force_env)
       END DO

       ! Create QM/MM Environment..
       CALL qmmm_env_create(qmmm_env,error=error)

       ! First Initialize Fist...
       globenv%program_name="FIST"       
       CALL fist_create_force_env(sub_force_env(1)%force_env, &
            globenv, error=error)
 
       CALL force_env_get(sub_force_env(1)%force_env,subsys=subsys_mm,&
            cell=mm_cell,error=error) 
       CPPrecondition(SIZE(subsys_mm)==1,cp_failure_level,routineP,error,failure)

       IF (.NOT.ASSOCIATED(qm_cell_small)) THEN
          ALLOCATE (qm_cell_small,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL cell_clone(mm_cell, qm_cell_small,.true.)

       ! Set up QM/MM Options
       CALL qmmm_namelist_read(globenv,&
                               qmmm_env,&
                               subsys_mm,&
                               qm_atom_type,&
                               qm_atom_index,&
                               mm_atom_index,&
                               qm_subsys_index,&
                               qm_cell_small,&
                               qmmm_coupl_type,&
                               eps_mm_rspace,&
                               error)
                               
       qmmm_env%qm_atom_index   => qm_atom_index
       qmmm_env%mm_atom_index   => mm_atom_index
       qmmm_env%qm_subsys_index => qm_subsys_index
       qmmm_env%eps_mm_rspace   = eps_mm_rspace
       qmmm_env%qmmm_coupl_type = qmmm_coupl_type
       qmmm_env%num_qm_atoms = size(qm_atom_index)
       qmmm_env%num_mm_atoms = size(mm_atom_index)
       Write(iw,'(A,I6)')"Number of QM atoms: ",qmmm_env%num_qm_atoms 
       Write(iw,'(A,I6)')"Number of MM atoms: ",qmmm_env%num_mm_atoms
       !
       ! Assign charges and mm_el_pot_radius from fist_topology
       !
       ALLOCATE(mm_atom_chrg(SIZE(mm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mm_el_pot_radius(SIZE(mm_atom_index)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       mm_atom_chrg     = 0.0_dbl
       mm_el_pot_radius = 0.0_dbl
       
       CALL assign_mm_charges_and_radius(subsys=subsys_mm,&
                                         mm_atom_chrg=mm_atom_chrg,&
                                         mm_el_pot_radius=mm_el_pot_radius,&
                                         mm_atom_index=mm_atom_index,&
                                         error=error)

       qmmm_env%mm_atom_chrg     => mm_atom_chrg
       qmmm_env%mm_el_pot_radius => mm_el_pot_radius


       globenv%program_name="QS"
       CALL create_small_subsystem(subsys_qm,&
            big_subsys=subsys_mm,small_para_env=para_env,&
            small_cell=qm_cell_small,sub_atom_index=qm_atom_index,&
            sub_subsys_index=qm_subsys_index,&
            sub_atom_kind_name=qm_atom_type, globenv=globenv,error=error)
       CALL read_atomic_kind_set(subsys_qm%atomic_kinds%els,globenv)
     
       CALL qs_env_create(qs_env, para_env=para_env, error=error)

       CALL qs_init_subsys(qs_env,subsystem=subsys_qm,cell=qm_cell_small,&
            cell_ref=qm_cell_small, use_ref_cell=.FALSE.,qmmm=.true.,&
            globenv=globenv,&
            error=error)

       CALL cp_subsys_release(subsys_qm,error=error)
       CALL force_env_create(sub_force_env(2)%force_env,&
            qs_env=qs_env,globenv=globenv, error=error)

       !
       ! Initialize MM Potential fitted with Gaussian
       !
       CALL qmmm_gaussian_initialize(qmmm_env=qmmm_env,&
                                     qm_cell_small=qm_cell_small,&
                                     mm_cell=mm_cell,&
                                     globenv=globenv,&
                                     pw_env=qs_env%pw_env,&
                                     error=error)
       !
       ! Initialize the MM potential stored on vector
       !
       CALL qmmm_potential_init(qmmm_env=qmmm_env,&
                                qm_cell_small=qm_cell_small,&
                                mm_cell=mm_cell,&
                                globenv=globenv,&
                                pw_env=qs_env%pw_env,&
                                error=error)
       !
       ! Initialize the qmmm_pw_grid
       !
       CALL qmmm_pw_grid_init(qmmm_env=qmmm_env,&
                              pw_env=qs_env%pw_env,&
                              cell=qm_cell_small,&
                              error=error)

       CALL cell_release(qm_cell_small, error)                       
       CALL qs_env_release(qs_env,error=error)

       globenv%program_name="QMMM"
       CALL force_env_create(force_env, sub_force_env=sub_force_env,&
            qmmm_env=qmmm_env,globenv=globenv,error=error)

       CALL qmmm_env_release(qmmm_env, error)
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

    INTEGER                                  :: iat, isubsys, nat, stat
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

    CPPrecondition(.not.ASSOCIATED(small_subsystem),cp_failure_level,routineP,error,failure)
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
       !CALL read_topology_section(topology,globenv)

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
                                         mm_atom_index, error)
   IMPLICIT NONE
   ! Arguments
   TYPE(cp_subsystem_p_type), &
        DIMENSION(:), POINTER                  :: subsys
   REAL(dp), DIMENSION(:), POINTER             :: mm_atom_chrg, mm_el_pot_radius
   INTEGER, DIMENSION(:), POINTER              :: mm_atom_index
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error
   ! Local Variables
   TYPE(atomic_kind_type), POINTER             :: my_kind
    TYPE(elp_potential_type), POINTER          :: my_potential
   TYPE(particle_list_type), POINTER           :: particles
   TYPE(particle_type), DIMENSION(:), &
        POINTER                                :: particle_set
   INTEGER  :: I, IndMM
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
      WRITE(*,'(A,I5,A,F12.6,A,F12.6)')'MM ATOM:',IndMM,' RADIUS: ',ri,&
                                       ' CHARGE: ',qi
   END DO
   
 END SUBROUTINE assign_mm_charges_and_radius

END MODULE qmmm_init
