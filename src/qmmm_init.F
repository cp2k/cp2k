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
!!     performs a QM/MM calculation
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
                                             read_atomic_kind_set
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
  USE qmmm_types,                      ONLY: qmmm_env_create,&
                                             qmmm_env_type
  USE qs_environment,                  ONLY: qs_init_subsys
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_env_release,&
                                             qs_environment_type
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE simulation_cell,                 ONLY: cell_retain,&
                                             cell_type,&
                                             init_qmmm_cell
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: connectivity_control,&
                                             force_field_control
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
    CHARACTER(len=default_string_length)     :: keyword, qm_atom_kind
    INTEGER                                  :: handle, isubf, stat,string_length
    INTEGER                                  :: num_qm_atoms, i
    INTEGER                                  :: qm_num_atom_kind, qm_num_atom_tot
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index, qm_subsys_index
    LOGICAL                                  :: failure, section_found, keyword_found
    TYPE(cell_type), POINTER                 :: qm_cell, qm_cell_small
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys_mm
    TYPE(cp_subsystem_type), POINTER         :: subsys_qm
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: sub_force_env
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    CHARACTER(len=3)                         :: test_result

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    NULLIFY(sub_force_env, qm_atom_index, qm_atom_type,qmmm_env, subsys_mm,&
         subsys_qm, qm_cell, qm_cell_small, qs_env, para_env)
    para_env => globenv%para_env
    IF (.NOT. failure) THEN
       ALLOCATE(sub_force_env(2),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO isubf=1,2
          NULLIFY(sub_force_env(isubf)%force_env)
       END DO
       globenv%program_name="QS"
       
       ! quando sara' pronto fist bisogna modificare questa chiamata...
       CALL quickstep_create_force_env(sub_force_env(1)%force_env, &
            globenv, error=error)

       CALL qmmm_env_create(qmmm_env,error=error)
 
       CALL force_env_get(sub_force_env(1)%force_env,subsys=subsys_mm,&
            cell=qm_cell,error=error) 
       CPPrecondition(SIZE(subsys_mm)==1,cp_failure_level,routineP,error,failure)

       IF (.NOT.ASSOCIATED(qm_cell_small)) THEN
          ALLOCATE (qm_cell_small,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       qm_cell_small = qm_cell

       CALL qmmm_namelist_read
       qmmm_env%qm_atom_index => qm_atom_index
       qmmm_env%qm_subsys_index => qm_subsys_index

       CALL create_small_subsystem(subsys_qm,&
            big_subsys=subsys_mm,small_para_env=para_env,&
            small_cell=qm_cell_small,sub_atom_index=qm_atom_index,&
            sub_subsys_index=qm_subsys_index,&
            sub_atom_kind_name=qm_atom_type, globenv=globenv,error=error)
       CALL read_atomic_kind_set(subsys_qm%atomic_kinds%els,globenv)
     
       CALL qs_env_create(qs_env, para_env=para_env, error=error)

       CALL qs_init_subsys(qs_env,subsystem=subsys_qm,cell=qm_cell_small,&
            cell_ref=qm_cell_small, use_ref_cell=.FALSE.,globenv=globenv,&
            error=error)

       CALL cp_subsys_release(subsys_qm,error=error)
       CALL force_env_create(sub_force_env(2)%force_env,&
            qs_env=qs_env,globenv=globenv,error=error)

       CALL qs_env_release(qs_env,error=error)

       globenv%program_name="QMMM"
       CALL force_env_create(force_env, sub_force_env=sub_force_env,&
            qmmm_env=qmmm_env,globenv=globenv,error=error)

       NULLIFY(qm_atom_index,qm_subsys_index)
       DEALLOCATE(qm_atom_type,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    CALL timestop(0.0_dp,handle)
    !!
    CONTAINS
      !!****f* qmmm_init/qmmm_create_force_env/qmmm_namelist_read [1.0] *
      !!
      !!   NAME
      !!     qmmm_namelist_read
      !!
      !!   FUNCTION
      !!     Read "QMMM" Namelist Arguments
      !!     Arguments:
      !!     &QMMM
      !!        QM_ATOMS Number_of_QM_atoms
      !!        Atom_LABEL Number_of_this_atomic_Kind
      !!        Index1 Index2 ... IndexN
      !!        Atom_LABEL Number_of_this_atomic_Kind
      !!        Index1 Index2 ... IndexL
      !!   
      !!        ABC_QM x y z
      !!     &END
      !!   NOTES
      !!     -
      !!
      !!   INPUTS
      !!     - error: variable to control error logging, stopping,... 
      !!       see module cp_error_handling 
      !!
      !!   AUTHOR
      !!     Teodoro Laino
      !!
      !!   MODIFICATION HISTORY
      !!     05.2004 created [tlaino]
      !!
      !!*************************************************************************
      SUBROUTINE qmmm_namelist_read
        Implicit None
        CHARACTER(len=*), PARAMETER :: subroutineN = 'qmmm_namelist_read',&
             subroutineP1 = moduleN//':'//routineN//":"//subroutineN

        keyword_found   = .false.
        !
        ! Check for &QMMM Namelist
        !
        CALL start_parser(file_name=globenv%input_file_name,&
             globenv=globenv,&
             start_section_label="QMMM",&
             section_found=section_found)
       
        If (section_found) then
           !
           ! Read Info from QMMM Namelist
           !
           DO
              test_result = test_object(newline=.TRUE.)
              IF (test_result == "EOS") THEN
                 EXIT
              ELSE IF (test_result == "EOF") THEN
                 CALL stop_parser(subroutineP1,modulen,__LINE__,"EOF")
              END IF
              CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
              SELECT CASE(TRIM(keyword))
              CASE("QMATOMS","QM_ATOMS")
                 qm_num_atom_tot = 0
                 CALL read_object(num_qm_atoms)
                 ALLOCATE(qm_atom_index(num_qm_atoms),& 
                          qm_atom_type(num_qm_atoms),&
                          qm_subsys_index(num_qm_atoms),&
                          stat=stat)
                 CPPostcondition(stat==0,cp_failure_level,subroutineP1,error,failure)
                 ! Read the QM Atomic Infos..
                 do while (qm_num_atom_tot /= num_qm_atoms)
                    CALL read_object(qm_atom_kind,newline=.TRUE.)
                    CALL read_object(qm_num_atom_kind)
                    if (qm_num_atom_kind /= 0) then
                       i = qm_num_atom_tot + 1
                       CALL read_object(qm_atom_index(i),newline=.TRUE.)
                       qm_atom_type(i) = qm_atom_kind
                       do i = qm_num_atom_tot + 2, qm_num_atom_kind + qm_num_atom_tot
                          CALL read_object(qm_atom_index(i))
                          qm_atom_type(i) = qm_atom_kind
                       End do
                    End If
                    qm_num_atom_tot = qm_num_atom_tot + qm_num_atom_kind
                 End do
                 qm_subsys_index=1
                 keyword_found = .true.
              CASE ("A_QM","B_QM","C_QM","ABC_QM")
                 SELECT CASE(TRIM(keyword))
                 CASE ("A_QM")
                    CALL read_object(qm_cell_small%hmat(1,1))
                    CALL read_object(qm_cell_small%hmat(2,1))
                    CALL read_object(qm_cell_small%hmat(3,1))
                 CASE ("ABC_QM")
                    CALL read_object(qm_cell_small%hmat(1,1))
                    CALL read_object(qm_cell_small%hmat(2,2))
                    CALL read_object(qm_cell_small%hmat(3,3))
                 CASE ("B_QM")
                    CALL read_object(qm_cell_small%hmat(1,2))
                    CALL read_object(qm_cell_small%hmat(2,2))
                    CALL read_object(qm_cell_small%hmat(3,2))
                 CASE ("C_QM")
                    CALL read_object(qm_cell_small%hmat(1,3))
                    CALL read_object(qm_cell_small%hmat(2,3))
                    CALL read_object(qm_cell_small%hmat(3,3))
                 END SELECT
                 call init_qmmm_cell(qm_cell_small)
              CASE("&END")
                 EXIT 
              CASE DEFAULT
                 CALL stop_parser(subroutineP1,"Invalid keyword found")
              END SELECT              
           END DO
           CPPostcondition(keyword_found,cp_failure_level,subroutineP1,error,failure)
        else
           !
           ! I keep this part for the regtest...
           ! 
           ALLOCATE(qm_atom_index(3), qm_atom_type(3),qm_subsys_index(3), stat=stat)
           CPPostcondition(stat==0,cp_failure_level,subroutineP1,error,failure)
           qm_atom_index=(/ 1,2,3 /)
           qm_atom_type=(/ "O","H","H"/)
           qm_subsys_index=1
        endif
        CALL finish_parser()
      END SUBROUTINE qmmm_namelist_read
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
       ! 4. Read in the coordinates
       !-----------------------------------------------------------------------------
       !CALL coordinate_control(topology,globenv)


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
!FM  IF((topology%ff).AND.(PRESENT(potparm))) THEN
!FM      CALL topology_force_field_pack(particle_set,atomic_kind_set,&
!FM                                     molecule_kind_set,molecule_set,&
!FM                                     topology,globenv,potparm)
!FM  END IF


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


END MODULE qmmm_init
