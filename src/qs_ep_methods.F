!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_ep_methods [1.0] *
!!
!!   NAME
!!     qs_ep_methods
!!
!!   FUNCTION
!!     in this module there are the routines for the calculation of 
!!     psi and energy without scf, but perturbatively as in PRL 82,22;
!!     26 Nov 2001,"Accurate Total Energies without Self-Consistency",
!!     David M. Benoit, Daniel Sebastiani, and Michele Parrinello 
!!
!!   NOTES
!!     At the moment much is hard-coded
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ep_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set,&
                                             set_atomic_kind
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_get_mstruct
  USE cp_fm_struct,                    ONLY: cp_fm_struct_get,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dft_types,                       ONLY: dft_control_type
  USE global_types,                    ONLY: DEBUG,&
                                             global_environment_type
  USE greens_fn,                       ONLY: pw_green_fn_init
  USE kinds,                           ONLY: wp=>dp
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_core,&
                                             calculate_rho_elec,&
                                             calculate_total_rho,&
                                             calculate_wavefunction
  USE qs_core_energies,                ONLY: calculate_ecore,&
                                             calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_rebuild_aomo_pools,&
                                             qs_env_rebuild_rho,&
                                             qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             init_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ep_types,                     ONLY: qs_ep_env_type,&
                                             qs_ep_release,&
                                             qs_ep_retain
  USE qs_ks_methods,                   ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_p_type,&
                                             qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             calculate_density_matrix,&
                                             calculate_subspace_eigenvalues,&
                                             correct_mo_eigenvalues,&
                                             get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             read_mo_set,&
                                             set_mo_occupation,&
                                             write_mo_set
  USE qs_neighbor_lists,               ONLY: build_neighbor_lists
  USE qs_overlap,                      ONLY: distribute_overlap_matrix
  USE qs_p_sparse_psi,                 ONLY: p_proj_create,&
                                             p_proj_release,&
                                             p_proj_retain,&
                                             p_proj_transfer_blacs_to_f,&
                                             qs_p_projection_p_type,&
                                             qs_p_projection_type
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_functionals_tests,            ONLY: xc_functional_test
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ep_methods'
  INTEGER, PRIVATE, SAVE :: last_ep_env_id=0

  PUBLIC :: qs_ep_calculation

!!***
!****************************************************************************
CONTAINS

!!****f* qs_ep_methods/qs_ep_create [1.0] *
!!
!!   NAME
!!     qs_ep_create
!!
!!   SYNOPSIS
!!     Subroutine qs_ep_calculation(qs_env, globenv, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), Intent (INOUT), Target:: globenv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ep_calculation
!!
!!   FUNCTION
!!     creates the environement to perform an ep perturbation
!!
!!   NOTES
!!     Don't look at the code, it's in development!
!!     Whole initialization of qs_env should be retought (also for dynamics)
!!
!!   INPUTS
!!     - ep_env: the ep_env env to create
!!     - main_qs_env: the qs_env of the full system
!!     - projections: the projections into the subsystems
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_ep_create(ep_env,main_qs_env,projections, error)
    TYPE(qs_ep_env_type), POINTER :: ep_env
    TYPE(qs_environment_type), INTENT(inout),TARGET :: main_qs_env
    TYPE(qs_p_projection_p_type), DIMENSION(:), INTENT(in) :: projections
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='qs_ep_create',&
         routineP=moduleN//':'//routineN
!FM    INTEGER :: ikind, isub, iatom, natoms, natoms_kind, ispin, stat,&
!FM         handle, nat1, nat2, n_ao, n_mo, nelectron, next_cgf, next_sgf
!FM    INTEGER, DIMENSION(:), POINTER :: atom_list, atom_list1, atom_list2
!FM    TYPE(dft_control_type), POINTER :: dft_control
!FM    TYPE(scf_control_type), POINTER :: scf_control
!FM    TYPE(cp_para_env_type), POINTER :: para_env
!FM    TYPE(mo_set_p_type), DIMENSION(:), POINTER :: c
!FM    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: k,p,h,s
!FM    TYPE(cp_logger_type), POINTER :: logger
!FM    TYPE(cp_full_matrix_type), POINTER :: eigenvectors
!FM    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
!FM    TYPE(qs_charges_type), POINTER :: qs_charges
!FM    TYPE(qs_p_projection_p_type), DIMENSION(:), POINTER :: sub_proj
!FM    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
!FM    TYPE(qs_ks_env_type), POINTER :: ks_env
!FM    TYPE(coeff_type), POINTER   :: rho_core
!FM    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
!FM    TYPE(pw_env_type), POINTER :: pw_env
!FM    TYPE(cell_type), POINTER :: cell

!FM    failure=.FALSE.
!FM    NULLIFY(atom_list, atom_list1, atom_list2,dft_control,scf_control,&
!FM         para_env,c,k,p,h,s,logger,eigenvectors,atomic_kind_set,qs_charges,&
!FM         sub_proj, particle_set, ks_env, rho_core, auxbas_pw_pool, pw_env)
!FM
!FM    allocate(ep_env, stat=stat)
!FM    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM    if (.not. failure) then
!FM       ep_env%ref_count=1
!FM       last_ep_env_id=last_ep_env_id+1
!FM       ep_env%id_nr=last_ep_env_id
!FM       allocate(ep_env%sub_qs_env(size(projections)),&
!FM            ep_env%sub_proj(size(projections)),stat=stat)
!FM       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM       
!FM    end if
!FM    DO isub_qs=1,2
!FM       CALL init_qs_env(sub_qs_env(isub_qs),para_env=para_env,&
!FM            blacs_env=qs_env%blacs_env,error=error)
!FM    END DO
!FM    natoms=SIZE(qs_env%particle_set)
!FM    ALLOCATE(sub_qs_env(1)%particle_set(natoms/2),&
!FM         sub_qs_env(2)%particle_set((natoms+1)/2),&
!FM         stat=stat)
!FM    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM
!FM    sub_qs_env(1)%particle_set = qs_env%particle_set(1:(natoms/2))
!FM    sub_qs_env(2)%particle_set = qs_env%particle_set((natoms/2+1):natoms)
!FM
!FM    ALLOCATE(sub_qs_env(1)%atomic_kind_set(SIZE(qs_env%atomic_kind_set)),&
!FM         sub_qs_env(2)%atomic_kind_set(SIZE(qs_env%atomic_kind_set)),&
!FM         stat=stat)
!FM    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM
!FM    DO isub_qs=1,2
!FM       sub_qs_env(isub_qs)%atomic_kind_set = qs_env%atomic_kind_set
!FM       ALLOCATE(sub_qs_env(isub_qs)%cell, stat=stat)
!FM       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM       sub_qs_env(isub_qs)%cell = qs_env%cell
!FM       sub_qs_env(isub_qs)%dft_control => qs_env%dft_control
!FM       sub_qs_env(isub_qs)%scf_control => qs_env%scf_control
!FM       CALL scf_c_retain(qs_env%scf_control, error=error)
!FM       sub_qs_env(isub_qs)%pw_env => qs_env%pw_env
!FM       CALL pw_env_retain(qs_env%pw_env, error=error)
!FM    END DO
!FM
!FM    ! atomic_kind_set setup
!FM    DO ikind=1,SIZE(qs_env%atomic_kind_set)
!FM       CALL get_atomic_kind(qs_env%atomic_kind_set(ikind),&
!FM            atom_list=atom_list)
!FM
!FM       natoms_kind=0
!FM       DO iatom=1,SIZE(atom_list)
!FM          IF (atom_list(iatom) <= natoms/2) natoms_kind=natoms_kind+1
!FM       END DO
!FM
!FM       ALLOCATE(atom_list1(natoms_kind),&
!FM            atom_list2(SIZE(atom_list)-natoms_kind),&
!FM            stat=stat)
!FM       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM
!FM       nat1=0; nat2=0
!FM       DO iatom=1,SIZE(atom_list)
!FM          IF (atom_list(iatom) <= natoms/2) THEN
!FM             nat1=nat1+1
!FM             atom_list1(nat1)=atom_list(iatom)
!FM             sub_qs_env(1)%particle_set(atom_list1(nat1))%atomic_kind &
!FM                  => sub_qs_env(1)%atomic_kind_set(ikind)
!FM          ELSE
!FM             nat2=nat2+1
!FM             atom_list2(nat2)=atom_list(iatom)-natoms/2
!FM             sub_qs_env(2)%particle_set(atom_list2(nat2))%atomic_kind &
!FM                  => sub_qs_env(2)%atomic_kind_set(ikind)
!FM          END IF
!FM       END DO
!FM
!FM       CALL set_atomic_kind(sub_qs_env(1)%atomic_kind_set(ikind),&
!FM            atom_list=atom_list1)
!FM       CALL set_atomic_kind(sub_qs_env(2)%atomic_kind_set(ikind),&
!FM            atom_list=atom_list2)
!FM    END DO
!FM
!FM    DO isub_qs=1,2
!FM       particle_set => sub_qs_env(isub_qs)%particle_set
!FM       next_cgf=1
!FM       next_sgf=1
!FM       DO iatom=1,SIZE(particle_set)
!FM          particle_set(iatom)%last_cgf=next_cgf+particle_set(iatom)%last_cgf-&
!FM               particle_set(iatom)%first_cgf
!FM          particle_set(iatom)%first_cgf=next_cgf
!FM          next_cgf=particle_set(iatom)%last_cgf+1
!FM          particle_set(iatom)%last_sgf=next_sgf+particle_set(iatom)%last_sgf-&
!FM               particle_set(iatom)%first_sgf
!FM          particle_set(iatom)%first_sgf=next_sgf
!FM          next_sgf=particle_set(iatom)%last_sgf+1
!FM          !FM cgf_aux not initialized! (at the Moment unused)
!FM          particle_set(iatom)%iatom=iatom
!FM       END DO
!FM    END DO
!FM
!FM    DO isub_qs=1,2
!FM       !   *** allocate mo_types ***
!FM       CALL get_qs_env(sub_qs_env(isub_qs), atomic_kind_set=atomic_kind_set, &
!FM            dft_control=dft_control, scf_control=scf_control)
!FM       CALL get_atomic_kind_set(atomic_kind_set, nsgf=n_ao, nelectron=nelectron)
!FM       nelectron = nelectron-dft_control%charge
!FM       IF (globenv%print%mo_eigenvalues.OR.&
!FM            (scf_control%level_shift /= 0.0_wp).OR.&
!FM            (scf_control%smear /= 0.0_wp).OR.&
!FM            (scf_control%eps_jacobi /= 0.0_wp)) THEN
!FM          n_mo = n_ao
!FM       ELSE
!FM          n_mo = nelectron/2
!FM       END IF
!FM       ALLOCATE(c(dft_control%nspins),stat=stat)
!FM       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM       IF (.NOT.failure) THEN
!FM          DO ispin=1,dft_control%nspins
!FM             NULLIFY(c(ispin)%mo_set)
!FM             IF (dft_control%nspins==1) THEN
!FM                CALL allocate_mo_set(c(ispin)%mo_set,&
!FM                     nao=n_ao,nmo=n_mo,nelectron=nelectron,maxocc=2.0_wp,&
!FM                     error=error)
!FM             ELSE
!FM                CALL allocate_mo_set(c(ispin)%mo_set,&
!FM                     nao=n_ao,nmo=n_mo,nelectron=nelectron,maxocc=1.0_wp,&
!FM                     error=error)
!FM             END IF
!FM          END DO
!FM       END IF
!FM       CALL set_qs_env(sub_qs_env(isub_qs), c=c, error=error)
!FM
!FM       !   *** Allocate the data structure for Quickstep energies ***
!FM       CALL allocate_qs_energy(sub_qs_env(isub_qs)%energy)
!FM       ! allocate qs_charges
!FM       CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
!FM       CALL set_qs_env(sub_qs_env(isub_qs), qs_charges=qs_charges,error=error)
!FM       CALL qs_charges_release(qs_charges, error=error)
!FM       ! rebuilds ao_mo_fm_pools
!FM       CALL qs_env_rebuild_aomo_pools(sub_qs_env(isub_qs), error=error)
!FM       ! init green's function
!FM       CALL get_qs_env(sub_qs_env(isub_qs), cell=cell, pw_env=pw_env,&
!FM            rho_core=rho_core, error=error)
!FM       CALL pw_green_fn_init(cell%green,&
!FM            pw_env%pw_pools(pw_env%auxbas_grid)%pool%pw_grid)
!FM       ! reallocate rho_core
!FM       IF (ASSOCIATED(rho_core)) THEN
!FM          CALL coeff_deallocate(rho_core)
!FM          DEALLOCATE(rho_core,stat=stat)
!FM          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
!FM       ELSE
!FM          ALLOCATE(rho_core,stat=stat)
!FM          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!FM       END IF
!FM       IF (.NOT.failure) THEN
!FM          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
!FM          CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=rho_core, &
!FM               use_data=COMPLEXDATA1D, error=error)
!FM          rho_core%pw%in_space=RECIPROCALSPACE
!FM          CALL set_qs_env(sub_qs_env(isub_qs), rho_core=rho_core,error=error)
!FM       END IF
!FM       
!FM    END DO
!FM
!FM    CALL timestop(0.0_wp,handle)

  END SUBROUTINE qs_ep_create
!***************************************************************************

!!****f* qs_ep_methods/qs_ep_calculation [1.0] *
!!
!!   NAME
!!     qs_ep_calculation
!!
!!   SYNOPSIS
!!     Subroutine qs_ep_calculation(qs_env, globenv, error)
!!       Type(qs_environment_type), Intent (INOUT), Target:: qs_env
!!       Type(global_environment_type), Intent (INOUT), Target:: globenv
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_ep_calculation
!!
!!   FUNCTION
!!     Perform an ep calculation
!!
!!   NOTES
!!     Don't look at the code, it's in development!
!!     Whole initialization of qs_env should be retought (also for dynamics)
!!
!!   INPUTS
!!     - qs_env: the qs_env of the full system as returned after qs_init
!!     - globenv: globenv, needed for some initializations
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ep_calculation(qs_env,globenv, error)
  TYPE(qs_environment_type), INTENT(inout),TARGET :: qs_env
  TYPE(global_environment_type), INTENT(inout), TARGET :: globenv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_ep_calculation',&
        routineP=moduleN//':'//routineN
  INTEGER :: isub, iatom, natoms, handle, stat
  TYPE(qs_p_projection_p_type), DIMENSION(2) :: sub_proj
  TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
  TYPE(qs_ep_env_type), POINTER :: ep_env
  TYPE(cp_logger_type), POINTER :: logger
  INTEGER, DIMENSION(:), ALLOCATABLE :: atoms

  CALL timeset(routineN,"I","",handle)
  
  failure=.FALSE.
  NULLIFY(particle_set,ep_env)
  logger => cp_error_get_logger(error)

  CALL cp_log(logger,cp_note_level,routineP,&
       "Finishing initialization of qs_env")
  CALL qs_ep_finish_qs_init(qs_env, globenv, error)
  CALL cp_error_check(error,failure)

  ! define projections (i.e. subsystems)
  IF (.NOT. failure) THEN
     CALL get_qs_env(qs_env,particle_set=particle_set)
     natoms=SIZE(particle_set)

     NULLIFY(sub_proj(1)%projection, sub_proj(2)%projection)

     ALLOCATE(atoms(natoms/2),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO iatom=1,natoms/2
        atoms(iatom)=iatom
     END DO
     CALL p_proj_create(sub_proj(1)%projection, atoms=atoms,&
          particle_set=particle_set, error=error)
     DEALLOCATE(atoms,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

     ALLOCATE(atoms((natoms+1)/2),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO iatom=1,(natoms+1)/2
        atoms(iatom)=iatom+natoms/2
     END DO
     CALL p_proj_create(sub_proj(1)%projection, atoms=atoms,&
          particle_set=particle_set, error=error)
     DEALLOCATE(atoms,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF

  ! create ep_env
  IF (.not.failure) THEN
     CALL qs_ep_create(ep_env,main_qs_env=qs_env, projections=sub_proj,&
          error=error)
  END IF

  ! energy in subsystem with diag
  DO isub=1,2
     CALL cp_log(logger,cp_note_level,routineP,&
          "Diagonalization of subsystem "//ADJUSTL(cp_to_string(isub)))
!     CALL qs_energies(sub_qs_env(isub),globenv)
  END DO
  
  ! perturbative approach for complete system
  CALL cp_log(logger,cp_note_level,routineP,&
       "Starting perturbation")

  CALL timestop(0.0_wp, handle)
END SUBROUTINE qs_ep_calculation
!***************************************************************************

!!****f* qs_ep_methods/qs_ep_finish_qs_init [1.0] *
!!
!!   NAME
!!     qs_ep_finish_qs_init
!!
!!   FUNCTION
!!     finish the initialization of the qs environement
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the qs_environment to fully initialize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ep_finish_qs_init(qs_env,globenv,error)
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(global_environment_type), INTENT(in), TARGET :: globenv
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_ep_finish_qs_init',&
       routineP=moduleN//':'//routineN
  INTEGER :: handle, stat, ispin
  TYPE(qs_ks_env_type), POINTER :: ks_env
  TYPE(mo_set_p_type), DIMENSION(:), POINTER :: c
  TYPE(cp_full_matrix_type), POINTER :: eigenvectors
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(dft_control_type), POINTER :: dft_control
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: k,p,h,s
 
  failure=.FALSE.
  NULLIFY(ks_env, c, eigenvectors, para_env, dft_control,k,p,h,s)

  CALL xc_functional_test ( globenv )
  CALL timeset(routineN,"I","",handle)

  CALL build_neighbor_lists(qs_env,globenv)

  CALL distribute_overlap_matrix(qs_env,globenv)
  para_env => qs_env%para_env
  CALL qs_env_update_s_mstruct(qs_env, error=error)

  !   *** Calculate the overlap and the core Hamiltonian integral matrix ***

  CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
       globenv=globenv,&
       calculate_forces=.FALSE.)

  !   *** Updates rho structure ***

  CALL qs_env_rebuild_rho(qs_env=qs_env, error=error)

  !   *** init the ks_env

  CALL get_qs_env(qs_env, ks_env=ks_env)
  IF (.not.ASSOCIATED(ks_env)) THEN
     CALL qs_ks_create(ks_env,qs_env=qs_env)
     CALL set_qs_env(qs_env, ks_env=ks_env)
  END IF

  !   *** Initializes the MOs ***
  CALL get_qs_env(qs_env,c=c)
  CPPrecondition(ASSOCIATED(c),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     DO ispin=1,SIZE(c)
        CALL get_mo_set(c(ispin)%mo_set,eigenvectors=eigenvectors)
        IF (.NOT.ASSOCIATED(eigenvectors)) THEN
           CALL init_mo_set(c(ispin)%mo_set, &
                qs_env%ao_mo_fm_pools(ispin)%pool,&
                error=error)
        END IF
     END DO
  END IF

  !   *** Allocate k and p and put them in the QS environment ***
  CALL get_qs_env(qs_env,dft_control=dft_control,s=s)
  CALL allocate_matrix_set(k,dft_control%nspins)
  CALL allocate_matrix_set(p,dft_control%nspins)
  DO ispin=1,dft_control%nspins
     NULLIFY (k(ispin)%matrix, p(ispin)%matrix)
     CALL replicate_matrix_structure(s(1)%matrix,k(ispin)%matrix,&
          "KOHN-SHAM MATRIX-"//ADJUSTL(cp_to_string(ispin)))
     CALL replicate_matrix_structure(s(1)%matrix,p(ispin)%matrix,&
          "DENSITY MATRIX-"//ADJUSTL(cp_to_string(ispin)))
  END DO
  CALL set_qs_env(qs_env=qs_env,&
       k=k,&
       p=p)

!FM!   *** Redistribute the core Hamiltonian matrix in ***
!FM!   *** order to eliminate redundant atomic blocks ***
!FM
!FM    CALL copy_sparse_to_blacs_matrix(h(1)%matrix,scf_work1)
!FM    IF ((globenv%num_pe > 1).AND.gth_potential_present) THEN
!FM      CALL deallocate_matrix(h(1)%matrix)
!FM      CALL replicate_matrix_structure(s(1)%matrix,h(1)%matrix,&
!FM                                      "CORE HAMILTONIAN MATRIX")
!FM      CALL copy_blacs_to_sparse_matrix(scf_work1,h(1)%matrix)
!FM      CALL set_qs_env(qs_env,h=h)
!FM    END IF

  CALL calculate_ecore_self(qs_env)
  CALL calculate_ecore_overlap(qs_env,globenv,.FALSE.)

END SUBROUTINE qs_ep_finish_qs_init
!***************************************************************************

END MODULE qs_ep_methods
