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
  USE qs_build_KS_matrix,              ONLY: qs_ks_create,&
                                             qs_ks_did_change,&
                                             qs_ks_update_qs_env
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
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ep_methods'

  PUBLIC :: qs_ep_calculation

!!***
!****************************************************************************
CONTAINS

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
    TYPE(qs_environment_type), DIMENSION(1:2) :: sub_qs_env
    TYPE(qs_ks_env_type), DIMENSION(1:2) :: sub_ks_env
    INTEGER :: ikind, isub_qs, iatom, natoms, natoms_kind, ispin, stat,&
         handle, nat1, nat2, n_ao, n_mo, nelectron, next_cgf, next_sgf
    INTEGER, DIMENSION(:), POINTER :: atom_list, atom_list1, atom_list2
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(scf_control_type), POINTER :: scf_control
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), POINTER :: c
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: k,p,h,s
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(cp_full_matrix_type), POINTER :: eigenvectors
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(qs_charges_type), POINTER :: qs_charges
    TYPE(qs_p_projection_p_type), DIMENSION(:), POINTER :: sub_proj
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    TYPE(qs_ks_env_type), POINTER :: ks_env
    TYPE(coeff_type), POINTER   :: rho_core
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(cell_type), POINTER :: cell

    failure=.FALSE.
    NULLIFY(atom_list, atom_list1, atom_list2,dft_control,scf_control,&
         para_env,c,k,p,h,s,logger,eigenvectors,atomic_kind_set,qs_charges,&
         sub_proj, particle_set, ks_env, rho_core, auxbas_pw_pool, pw_env)
    logger => cp_error_get_logger(error)

    CALL timeset(routineN,"I","",handle)

    ! finish pre initialization of qs_env
    CALL cp_log(logger,cp_note_level,routineP,&
         "Finishing initialization of qs_env")
    CALL build_neighbor_lists(qs_env,globenv)

    CALL distribute_overlap_matrix(qs_env,globenv)
    para_env => qs_env%para_env ! cp_para_env_from_globenv(globenv)
    CALL qs_env_update_s_mstruct(qs_env)

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


    ! split up qs_env
    CALL cp_log(logger,cp_note_level,routineP,&
         "Began splitting in subsystems ")
    DO isub_qs=1,2
       CALL init_qs_env(sub_qs_env(isub_qs),para_env=para_env,&
            blacs_env=qs_env%blacs_env,error=error)
    END DO
    natoms=SIZE(qs_env%particle_set)
    ALLOCATE(sub_qs_env(1)%particle_set(natoms/2),&
         sub_qs_env(2)%particle_set((natoms+1)/2),&
         stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    sub_qs_env(1)%particle_set = qs_env%particle_set(1:(natoms/2))
    sub_qs_env(2)%particle_set = qs_env%particle_set((natoms/2+1):natoms)

    ALLOCATE(sub_qs_env(1)%atomic_kind_set(SIZE(qs_env%atomic_kind_set)),&
         sub_qs_env(2)%atomic_kind_set(SIZE(qs_env%atomic_kind_set)),&
         stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO isub_qs=1,2
       sub_qs_env(isub_qs)%atomic_kind_set = qs_env%atomic_kind_set
       ALLOCATE(sub_qs_env(isub_qs)%cell, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       sub_qs_env(isub_qs)%cell = qs_env%cell
       sub_qs_env(isub_qs)%dft_control => qs_env%dft_control
       sub_qs_env(isub_qs)%scf_control => qs_env%scf_control
       CALL scf_c_retain(qs_env%scf_control, error=error)
       sub_qs_env(isub_qs)%pw_env => qs_env%pw_env
       CALL pw_env_retain(qs_env%pw_env, error=error)
    END DO

    ! atomic_kind_set setup
    CALL cp_log(logger,cp_note_level,routineP,&
         "atomic kind setup")
    DO ikind=1,SIZE(qs_env%atomic_kind_set)
       CALL get_atomic_kind(qs_env%atomic_kind_set(ikind),&
            atom_list=atom_list)

       natoms_kind=0
       DO iatom=1,SIZE(atom_list)
          IF (atom_list(iatom) <= natoms/2) natoms_kind=natoms_kind+1
       END DO

       ALLOCATE(atom_list1(natoms_kind),&
            atom_list2(SIZE(atom_list)-natoms_kind),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       nat1=0; nat2=0
       DO iatom=1,SIZE(atom_list)
          IF (atom_list(iatom) <= natoms/2) THEN
             nat1=nat1+1
             atom_list1(nat1)=atom_list(iatom)
             sub_qs_env(1)%particle_set(atom_list1(nat1))%atomic_kind &
                  => sub_qs_env(1)%atomic_kind_set(ikind)
          ELSE
             nat2=nat2+1
             atom_list2(nat2)=atom_list(iatom)-natoms/2
             sub_qs_env(2)%particle_set(atom_list2(nat2))%atomic_kind &
                  => sub_qs_env(2)%atomic_kind_set(ikind)
          END IF
       END DO

       CALL set_atomic_kind(sub_qs_env(1)%atomic_kind_set(ikind),&
            atom_list=atom_list1)
       CALL set_atomic_kind(sub_qs_env(2)%atomic_kind_set(ikind),&
            atom_list=atom_list2)
    END DO

    DO isub_qs=1,2
       particle_set => sub_qs_env(isub_qs)%particle_set
       next_cgf=1
       next_sgf=1
       DO iatom=1,SIZE(particle_set)
          particle_set(iatom)%last_cgf=next_cgf+particle_set(iatom)%last_cgf-&
               particle_set(iatom)%first_cgf
          particle_set(iatom)%first_cgf=next_cgf
          next_cgf=particle_set(iatom)%last_cgf+1
          particle_set(iatom)%last_sgf=next_sgf+particle_set(iatom)%last_sgf-&
               particle_set(iatom)%first_sgf
          particle_set(iatom)%first_sgf=next_sgf
          next_sgf=particle_set(iatom)%last_sgf+1
          !FM cgf_aux not initialized! (at the Moment unused)
          particle_set(iatom)%iatom=iatom
       END DO
    END DO

    DO isub_qs=1,2
       !   *** allocate mo_types ***
       CALL get_qs_env(sub_qs_env(isub_qs), atomic_kind_set=atomic_kind_set, &
            dft_control=dft_control, scf_control=scf_control)
       CALL get_atomic_kind_set(atomic_kind_set, nsgf=n_ao, nelectron=nelectron)
       nelectron = nelectron-dft_control%charge
       IF (globenv%print%mo_eigenvalues.OR.&
            (scf_control%level_shift /= 0.0_wp).OR.&
            (scf_control%smear /= 0.0_wp).OR.&
            (scf_control%eps_jacobi /= 0.0_wp)) THEN
          n_mo = n_ao
       ELSE
          n_mo = nelectron/2
       END IF
       ALLOCATE(c(dft_control%nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          DO ispin=1,dft_control%nspins
             NULLIFY(c(ispin)%mo_set)
             IF (dft_control%nspins==1) THEN
                CALL allocate_mo_set(c(ispin)%mo_set,&
                     nao=n_ao,nmo=n_mo,nelectron=nelectron,maxocc=2.0_wp,&
                     error=error)
             ELSE
                CALL allocate_mo_set(c(ispin)%mo_set,&
                     nao=n_ao,nmo=n_mo,nelectron=nelectron,maxocc=1.0_wp,&
                     error=error)
             END IF
          END DO
       END IF
       CALL set_qs_env(sub_qs_env(isub_qs), c=c, error=error)

       !   *** Allocate the data structure for Quickstep energies ***
       CALL allocate_qs_energy(sub_qs_env(isub_qs)%energy)
       ! allocate qs_charges
       CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
       CALL set_qs_env(sub_qs_env(isub_qs), qs_charges=qs_charges,error=error)
       CALL qs_charges_release(qs_charges, error=error)
       ! rebuilds ao_mo_fm_pools
       CALL qs_env_rebuild_aomo_pools(sub_qs_env(isub_qs), error=error)
       ! init green's function
       CALL get_qs_env(sub_qs_env(isub_qs), cell=cell, pw_env=pw_env,&
            rho_core=rho_core, error=error)
       CALL pw_green_fn_init(cell%green,&
            pw_env%pw_pools(pw_env%auxbas_grid)%pool%pw_grid)
       ! reallocate rho_core
       IF (ASSOCIATED(rho_core)) THEN
          CALL coeff_deallocate(rho_core)
          DEALLOCATE(rho_core,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ELSE
          ALLOCATE(rho_core,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
          CALL pw_pool_init_coeff(pool=auxbas_pw_pool, coeff=rho_core, &
               use_data=COMPLEXDATA1D, error=error)
          rho_core%pw%in_space=RECIPROCALSPACE
          CALL set_qs_env(sub_qs_env(isub_qs), rho_core=rho_core,error=error)
       END IF
       
    END DO

    ! energy in subsystem with diag
    DO isub_qs=1,2
       CALL cp_log(logger,cp_note_level,routineP,&
            "Diagonalization of subsystem "//ADJUSTL(cp_to_string(isub_qs)))
       CALL qs_energies(sub_qs_env(isub_qs),globenv)
    END DO

    ! perturbative approach for complete system
    CALL cp_log(logger,cp_note_level,routineP,&
         "Starting perturbation")

!    ALLOCATE(atom_list1(natoms/2),atom_list2((natoms+1)/2),stat=stat)
!FM    call p_proj_create(sub_proj(isub_qs)%projection, &
!FM         atoms=, proj_indexes, particle_set, error)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE qs_ep_calculation
!***************************************************************************

END MODULE qs_ep_methods
