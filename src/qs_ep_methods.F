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
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
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
                                             qs_env_create,&
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
    INTEGER :: ikind, isub, iatom, natoms, ispin, stat, isub_atom,&
         handle, nat, n_ao, n_mo, nelectron, next_cgf, next_sgf, nsub
    INTEGER, DIMENSION(:), POINTER :: atom_list, sub_atom_list
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(scf_control_type), POINTER :: scf_control
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    TYPE(mo_set_p_type), DIMENSION(:), POINTER :: sub_c
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(qs_charges_type), POINTER :: sub_qs_charges
    TYPE(qs_p_projection_type), POINTER :: proj
    TYPE(qs_environment_type), POINTER :: sub_qs
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set, p_set
    TYPE(coeff_type), POINTER   :: sub_rho_core
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(cell_type), POINTER :: sub_cell
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(atom_list, sub_atom_list, dft_control,scf_control,&
         para_env, blacs_env, sub_c,atomic_kind_set,sub_qs_charges,&
         proj, sub_qs, particle_set, p_set, sub_rho_core, auxbas_pw_pool,&
         pw_env, sub_cell, logger)

    ALLOCATE(ep_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       nsub=SIZE(projections)
       CALL get_qs_env(main_qs_env, dft_control=dft_control, &
            particle_set=particle_set,&
            atomic_kind_set=atomic_kind_set,&
            scf_control=scf_control,&
            pw_env=pw_env,&
            para_env=para_env,&
            blacs_env=blacs_env, error=error)

       ep_env%ref_count=1
       last_ep_env_id=last_ep_env_id+1
       ep_env%id_nr=last_ep_env_id
       ep_env%main_qs_env => main_qs_env
       ALLOCATE(ep_env%sub_qs_env(nsub),&
            ep_env%sub_proj(nsub),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       natoms=SIZE(particle_set)
       DO isub=1,nsub
          CALL qs_env_create(ep_env%sub_qs_env(isub)%qs_env,&
               para_env=para_env,&
               blacs_env=blacs_env,error=error)

          ep_env%sub_proj(isub)%projection => projections(isub)%projection
          CALL p_proj_retain(ep_env%sub_proj(isub)%projection, error=error)

          sub_qs => ep_env%sub_qs_env(isub)%qs_env
          proj => ep_env%sub_proj(isub)%projection
          ALLOCATE(sub_qs%particle_set(SIZE(proj%atoms)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO

       DO isub=1,nsub
          sub_qs => ep_env%sub_qs_env(isub)%qs_env
          proj => ep_env%sub_proj(isub)%projection
          p_set => sub_qs%particle_set

          DO iatom=1,SIZE(proj%atoms)
             p_set(iatom) = particle_set(proj%atoms(iatom))
          END DO

          ALLOCATE(sub_qs%atomic_kind_set(SIZE(atomic_kind_set)),&
               stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

          sub_qs%atomic_kind_set = atomic_kind_set
          ALLOCATE(sub_qs%cell, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          sub_qs%cell = ep_env%main_qs_env%cell
          sub_qs%dft_control => dft_control
          sub_qs%scf_control => scf_control
          CALL scf_c_retain(scf_control, error=error)
          sub_qs%pw_env => pw_env
          CALL pw_env_retain(pw_env, error=error)

          ! atomic_kind_set setup
          DO ikind=1,SIZE(atomic_kind_set)
             CALL get_atomic_kind(atomic_kind_set(ikind),&
                  atom_list=atom_list)

             ALLOCATE(sub_atom_list(SIZE(atom_list)),&
                  stat=stat)! alloc outside the do?
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

             nat=0
             DO iatom=1,SIZE(atom_list)
                DO isub_atom=1,SIZE(proj%atoms)
                   IF (atom_list(iatom) == proj%atoms(isub_atom)) THEN
                      nat=nat+1
                      sub_atom_list(nat)=isub_atom

                      sub_qs%particle_set(isub_atom)%atomic_kind &
                           => sub_qs%atomic_kind_set(ikind)
                      EXIT
                   END IF
                END DO
             END DO

             CALL set_atomic_kind(sub_qs%atomic_kind_set(ikind),&
                  atom_list=sub_atom_list(1:nat))
             DEALLOCATE(sub_atom_list, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END DO

          next_cgf=1
          next_sgf=1
          DO iatom=1,SIZE(p_set)
             p_set(iatom)%last_cgf=next_cgf+p_set(iatom)%last_cgf-&
                  p_set(iatom)%first_cgf
             p_set(iatom)%first_cgf=next_cgf
             next_cgf=p_set(iatom)%last_cgf+1
             p_set(iatom)%last_sgf=next_sgf+p_set(iatom)%last_sgf-&
                  p_set(iatom)%first_sgf
             p_set(iatom)%first_sgf=next_sgf
             next_sgf=p_set(iatom)%last_sgf+1
!FM cgf_aux not initialized! (at the moment unused)
             p_set(iatom)%iatom=iatom
          END DO

          !   *** allocate mo_types ***
          CALL get_atomic_kind_set(sub_qs%atomic_kind_set,&
               nsgf=n_ao, nelectron=nelectron)
          CPPrecondition(dft_control%charge==0,cp_failure_level,routineP,error,failure)
          nelectron = nelectron-dft_control%charge
          IF (logger%print_keys%mo_eigenvalues.OR.&
               (scf_control%level_shift /= 0.0_wp).OR.&
               (scf_control%smear /= 0.0_wp).OR.&
               (scf_control%eps_jacobi /= 0.0_wp)) THEN
             n_mo = n_ao
          ELSE
             n_mo = nelectron/2
          END IF
          ALLOCATE(sub_c(dft_control%nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             DO ispin=1,dft_control%nspins
                NULLIFY(sub_c(ispin)%mo_set)
                IF (dft_control%nspins==1) THEN
                   CALL allocate_mo_set(sub_c(ispin)%mo_set,&
                        nao=n_ao,nmo=n_mo,nelectron=nelectron,maxocc=2.0_wp,&
                        error=error)
                ELSE
                   CALL allocate_mo_set(sub_c(ispin)%mo_set,&
                        nao=n_ao,nmo=n_mo,nelectron=nelectron,maxocc=1.0_wp,&
                        error=error)
                END IF
             END DO
          END IF
          CALL set_qs_env(sub_qs, c=sub_c, error=error)

          !   *** Allocate the data structure for Quickstep energies ***
          CALL allocate_qs_energy(sub_qs%energy)
          ! allocate qs_charges
          CALL qs_charges_create(sub_qs_charges,nspins=dft_control%nspins,error=error)
          CALL set_qs_env(sub_qs, qs_charges=sub_qs_charges,error=error)
          CALL qs_charges_release(sub_qs_charges, error=error)
          ! rebuilds ao_mo_fm_pools
          CALL qs_env_rebuild_aomo_pools(sub_qs, error=error)
          ! init green's function
          CALL get_qs_env(sub_qs, cell=sub_cell, pw_env=pw_env,&
               rho_core=sub_rho_core, error=error)
          CALL pw_green_fn_init(sub_cell%green,&
               pw_env%pw_pools(pw_env%auxbas_grid)%pool%pw_grid)
          ! reallocate rho_core
          IF (ASSOCIATED(sub_rho_core)) THEN
             CALL coeff_deallocate(sub_rho_core)
             DEALLOCATE(sub_rho_core,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          ELSE
             ALLOCATE(sub_rho_core,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (.NOT.failure) THEN
             CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
             CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=sub_rho_core, &
                  use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
             CALL set_qs_env(sub_qs, rho_core=sub_rho_core,error=error)
          END IF

       END DO
    END IF

    CALL timestop(0.0_wp,handle)

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
  TYPE(qs_environment_type), POINTER :: sub_qs
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
     CALL p_proj_create(sub_proj(2)%projection, atoms=atoms,&
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
  DO isub=1,SIZE(ep_env%sub_qs_env)
     CALL cp_log(logger,cp_note_level,routineP,&
          "Diagonalization of subsystem "//ADJUSTL(cp_to_string(isub)))
     sub_qs => ep_env%sub_qs_env(isub)%qs_env
     CALL qs_energies(sub_qs,globenv)
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

  CALL timestop(0.0_wp,handle)
END SUBROUTINE qs_ep_finish_qs_init
!***************************************************************************

END MODULE qs_ep_methods
