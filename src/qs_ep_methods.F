
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
!!     psi and energy without scf, but perturbatively as in PRL 87(22);226401;
!!     26 Nov 2001,"Accurate Total Energies without Self-Consistency",
!!     David M. Benoit, Daniel Sebastiani, and Michele Parrinello 
!!
!!   NOTES
!!     Initialization procedure of qs_env not satisfactory
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
!MK  USE qs_overlap,                      ONLY: distribute_overlap_matrix
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_check,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_type,&
                                             cp_fm_write
  USE cp_fm_vect,                      ONLY: cp_fm_vect_set_all,&
                                             cp_fm_vect_write
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE green_methods,                   ONLY: pw_green_fn_rebuild
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE particle_types,                  ONLY: particle_type
  USE preconditioner,                  ONLY: apply_preconditioner,&
                                             init_preconditioner,&
                                             make_preconditioner
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_retain,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_type
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_create,&
                                             qs_environment_p_type,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ep_types,                     ONLY: ep_env_get,&
                                             qs_ep_env_type
  USE qs_kpp1_env_types,               ONLY: kpp1_write
  USE qs_ks_methods,                   ONLY: qs_ks_create
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             mpools_rebuild_fm_pools,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_psi0_changed,&
                                             p_op_l1,&
                                             p_op_l2_fawzi,&
                                             p_postortho,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: p_env_write,&
                                             qs_p_env_type
  USE qs_p_sparse_psi,                 ONLY: p_proj_create,&
                                             p_proj_retain,&
                                             p_proj_transfer_blacs_to_f,&
                                             qs_p_projection_p_type,&
                                             qs_p_projection_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE scf_control_types,               ONLY: scf_c_retain,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             cell_clone
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE green_types,                     ONLY: greens_function_type

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_ep_methods'
  INTEGER, PRIVATE, SAVE :: last_ep_env_id=0

  PUBLIC :: qs_ep_calculation, ep_env_create, ep_env_finish_qs_init

!!***
!****************************************************************************

CONTAINS

!!****f* qs_ep_methods/qs_ep_calculation [1.0] *
!!
!!   NAME
!!     qs_ep_calculation
!!
!!   FUNCTION
!!     Perform an ep calculation
!!
!!   NOTES
!!     still in active development, ugly code
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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ep_calculation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: file_unit, handle, iatom, &
                                                ispin, isub, lfomo, natoms, &
                                                nmo, nspins, stat, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atoms
    LOGICAL                                  :: failure, print_ep_matrixes
    REAL(KIND=dp)                            :: e_corr, e_corr_s, maxocc
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_environment_type), POINTER       :: sub_qs
    TYPE(qs_ep_env_type), POINTER            :: ep_env
    TYPE(qs_p_projection_p_type), &
      DIMENSION(2)                           :: sub_proj

  CALL timeset(routineN,"I","",handle)
  
  failure=.FALSE.
  NULLIFY(particle_set,ep_env,sub_qs,dft_control, energy, mos)
  logger => cp_error_get_logger(error)
  unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
  print_ep_matrixes=logger%print_keys%ep_matrixes

  CALL cp_log(logger,cp_note_level,routineP,&
       "Finishing initialization of qs_env")
  CALL ep_env_finish_qs_init(qs_env, globenv, error)
  CALL cp_error_check(error,failure)

  ! define projections (i.e. subsystems)
  CALL cp_log(logger,cp_note_level,routineP,&
       "Defining projections (i.e. subsystmes)")
  IF (.NOT. failure) THEN
     CALL get_qs_env(qs_env,particle_set=particle_set, &
          dft_control=dft_control, mos=mos, error=error)
     natoms=SIZE(particle_set)
     nspins=dft_control%nspins

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
     CALL cp_log(logger,cp_note_level,routineP,&
          "Creating ep_env")
     CALL ep_env_create(ep_env,main_qs_env=qs_env, projections=sub_proj,&
          error=error)
     CALL cp_log(logger,cp_note_level,routineP,&
          "Creating ep_env_created")
     CALL cp_error_check(error,failure)
  END IF

  IF (print_ep_matrixes) THEN
     IF (.not.failure) THEN
!FM        CALL cp_sm_output(logger=logger, &
!FM             outputName="ep_s", &
!FM             fromWhere=routineP, iter=0, &
!FM             matrix=ep_env%main_qs_env%matrix_s(1)%matrix,&
!FM             para_env=logger%para_env,error=error)
     END IF
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

  CALL ep_env_transfer_psi0(ep_env, error=error)

  CALL p_env_psi0_changed(ep_env%main_p_env, qs_env=ep_env%main_qs_env,&
       Hrho_psi0d=ep_env%m_pi_Hrho_psi0d, error=error)

  CALL p_postortho(p_env=ep_env%main_p_env,qs_env=ep_env%main_qs_env,&
          v=ep_env%m_pi_Hrho_psi0d, &
          n_cols=ep_env%main_p_env%n_mo,&
          error=error)
  DO ispin=1,SIZE(mos)
     CALL get_mo_set(mos(ispin)%mo_set,lfomo=lfomo, &
          nmo=nmo, maxocc=maxocc)
     IF (lfomo>nmo) THEN
        CALL cp_fm_scale_and_add(alpha=-maxocc,&
             matrix_a=ep_env%m_pi_Hrho_psi0d(ispin)%matrix,error=error)
     ELSE
        CALL cp_unimplemented_error(fromWhere=routineP,&
             message="symmetrized onesided smearing to do",&
             error=error)
     END IF
  END DO

  IF (print_ep_matrixes) THEN
     CALL cp_fm_vect_write(matrixes=ep_env%m_pi_Hrho_psi0d, &
          unit_nr=cp_logger_get_default_unit_nr(logger), &
          long_description=.TRUE.,&
          local=.FALSE., error=error)
     
     CALL p_env_write(ep_env%main_p_env,unit_nr=unit_nr,&
          long_description=debug_this_module, local=.FALSE.,&
          error=error)

!FM     DO ispin=1,nspins
!FM        CALL cp_sm_output(logger=logger, &
!FM             outputName="ep_ks-"//cp_to_string(ispin), &
!FM             fromWhere=routineP, iter=0, &
!FM             matrix=ep_env%main_qs_env%matrix_ks(ispin)%matrix,&
!FM             para_env=ep_env%main_qs_env%para_env,error=error)
!FM     END DO
  END IF

  CALL get_qs_env(ep_env%main_qs_env, energy=energy, error=error)
  WRITE (UNIT=unit_nr,FMT="(/,(T3,A,T60,F20.10))")&
       "Overlap energy of the core charge distribution:",energy%core_overlap,&
       "Self energy of the core charge distribution:   ",energy%core_self,&
       "Core Hamiltonian energy:                       ",energy%core,&
       "Hartree energy:                                ",energy%hartree,&
       "Exchange-correlation energy:                   ",energy%exc,&
       "Total energy:                                  ",energy%total

  IF (logger%print_keys%ep_operator) THEN
     CALL p_env_write_ep_matrix(p_env=ep_env%main_p_env, &
          qs_env=ep_env%main_qs_env, &
          name="ep_op", error=error)
  END IF

  CALL stupid_solve(ep_env,eps_r=1.0e-5_dp,error=error)

  ! calculate the energy correction (half from bra part, half from ket part)
  e_corr=0.0_dp
  DO ispin=1,nspins
     CALL cp_fm_trace(matrix_a=ep_env%m_pi_Hrho_psi0d(ispin)%matrix,&
          matrix_b=ep_env%psi1(ispin)%matrix,trace=e_corr_s,error=error)
     e_corr=e_corr-e_corr_s
  END DO
  e_corr=e_corr
  WRITE(unit=unit_nr,fmt="(a,es18.10)")" ep energy correction = ", e_corr

  IF (print_ep_matrixes) THEN
     CALL cp_fm_vect_write(matrixes=ep_env%psi1, &
          unit_nr=cp_logger_get_default_unit_nr(logger), &
          long_description=.TRUE.,&
          local=.FALSE., error=error)
  END IF
  
  IF (logger%print_keys%ep_rho1_cube) THEN
     CALL open_file(file_name="rho1.cube",file_status="UNKNOWN",&
          file_action="WRITE",unit_number=file_unit)
     CALL rs_pw_to_cube ( ep_env%main_p_env%rho1%rho_r(1)%pw, iunit=file_unit,&
          ionode=logger%para_env%mepos==logger%para_env%source,&
          title="rho1")
     CALL close_file(file_unit)
  END IF

  IF (logger%print_keys%cubes) THEN
     CALL open_file(file_name="rho0.cube",file_status="UNKNOWN",&
          file_action="WRITE",unit_number=file_unit)
     CALL rs_pw_to_cube ( ep_env%main_qs_env%rho%rho_r(1)%pw, iunit=file_unit,&
          ionode=logger%para_env%mepos==logger%para_env%source,&
          title="rho0")
     CALL close_file(file_unit)
  END IF

  IF (print_ep_matrixes) THEN
     CALL kpp1_write(kpp1_env=ep_env%main_p_env%kpp1_env,&
          unit_nr=unit_nr, long_description=.TRUE., local=.FALSE.,&
          error=error)
     CALL p_env_write(ep_env%main_p_env,unit_nr=unit_nr,&
          long_description=debug_this_module, local=.FALSE.,&
          error=error)
  END IF

  CALL timestop(0.0_dp, handle)
END SUBROUTINE qs_ep_calculation
!***************************************************************************

!!****f* qs_ep_methods/ep_env_create [1.0] *
!!
!!   NAME
!!     ep_env_create
!!
!!   FUNCTION
!!     creates the environement to perform an ep perturbation
!!
!!   NOTES
!!     Could take more advantage of what is calculated in main_qs_env
!!     (for example the S matrix)
!!     
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
  SUBROUTINE ep_env_create(ep_env,main_qs_env,projections, error)
    TYPE(qs_ep_env_type), POINTER            :: ep_env
    TYPE(qs_environment_type), POINTER       :: main_qs_env
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), INTENT(in)               :: projections
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, ikind, ispin, &
                                                isub, isub_atom, n_ao, n_mo, &
                                                nat, natoms, nelectron, nsub, &
                                                stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list, sub_atom_list
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: sub_cell
    TYPE(greens_function_type), POINTER      :: green
    TYPE(coeff_type), POINTER                :: sub_rho_core
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi0
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: sub_mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: p_set, particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_charges_type), POINTER           :: sub_qs_charges
    TYPE(qs_environment_type), POINTER       :: sub_qs
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(qs_p_projection_type), POINTER      :: proj
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(atom_list, sub_atom_list, dft_control,scf_control,&
         para_env, blacs_env, sub_mos,atomic_kind_set,sub_qs_charges,&
         proj, sub_qs, particle_set, p_set, sub_rho_core, auxbas_pw_pool,&
         pw_env, sub_cell, psi0, ao_mo_fm_pools, mpools, green)

    ALLOCATE(ep_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       NULLIFY(ep_env%sub_qs_env, ep_env%sub_proj, ep_env%main_qs_env,&
            ep_env%main_p_env, ep_env%m_pi_Hrho_psi0d, ep_env%psi1,&
            ep_env%precond)
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
               error=error)

          ep_env%sub_proj(isub)%projection => projections(isub)%projection
          CALL p_proj_retain(ep_env%sub_proj(isub)%projection, error=error)

          sub_qs => ep_env%sub_qs_env(isub)%qs_env
          proj => ep_env%sub_proj(isub)%projection
          !ALLOCATE(sub_qs%particle_set(SIZE(proj%atoms)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO

       DO isub=1,nsub
          sub_qs => ep_env%sub_qs_env(isub)%qs_env
          proj => ep_env%sub_proj(isub)%projection
!          p_set => sub_qs%particle_set

          DO iatom=1,SIZE(proj%atoms)
             p_set(iatom) = particle_set(proj%atoms(iatom))
          END DO

!          ALLOCATE(sub_qs%atomic_kind_set(SIZE(atomic_kind_set)),&
!               stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

!          sub_qs%atomic_kind_set = atomic_kind_set
          ALLOCATE(sub_qs%cell, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL cell_clone(ep_env%main_qs_env%cell,sub_qs%cell, error=error)
          sub_qs%dft_control => dft_control
          sub_qs%scf_control => scf_control
          CALL scf_c_retain(scf_control, error=error)
          sub_qs%pw_env => pw_env
          CALL pw_env_retain(pw_env, error=error)

          ! atomic_kind_set setup
          DO ikind=1,SIZE(atomic_kind_set)
             atomic_kind => atomic_kind_set(ikind)
             CALL get_atomic_kind(atomic_kind,atom_list=atom_list)

             ALLOCATE(sub_atom_list(SIZE(atom_list)),&
                  stat=stat)! alloc outside the do?
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

             nat=0
             DO iatom=1,SIZE(atom_list)
                DO isub_atom=1,SIZE(proj%atoms)
                   IF (atom_list(iatom) == proj%atoms(isub_atom)) THEN
                      nat=nat+1
                      sub_atom_list(nat)=isub_atom

!                      sub_qs%particle_set(isub_atom)%atomic_kind &
!                           => sub_qs%atomic_kind_set(ikind)
                      EXIT
                   END IF
                END DO
             END DO

!             CALL set_atomic_kind(sub_qs%atomic_kind_set(ikind),&
!                  atom_list=sub_atom_list(1:nat))
             DEALLOCATE(sub_atom_list, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END DO

!MK this components are obsolete, since they are now available
!MK should be retrieved on the fly by calling get_particle_set
!MK          next_cgf=1
!MK          next_sgf=1
!MK          DO iatom=1,SIZE(p_set)
!MK             p_set(iatom)%last_cgf=next_cgf+p_set(iatom)%last_cgf-&
!MK                  p_set(iatom)%first_cgf
!MK             p_set(iatom)%first_cgf=next_cgf
!MK             next_cgf=p_set(iatom)%last_cgf+1
!MK             p_set(iatom)%last_sgf=next_sgf+p_set(iatom)%last_sgf-&
!MK                  p_set(iatom)%first_sgf
!MK             p_set(iatom)%first_sgf=next_sgf
!MK             next_sgf=p_set(iatom)%last_sgf+1
!FM cgf_aux not initialized! (at the moment unused)
!MK             p_set(iatom)%iatom=iatom
!MK          END DO

          !   *** allocate mo_types ***
!          CALL get_atomic_kind_set(sub_qs%atomic_kind_set,&
!               nsgf=n_ao, nelectron=nelectron)
          CPPrecondition(dft_control%charge==0,cp_failure_level,routineP,error,failure)
          nelectron = nelectron-dft_control%charge
          IF (logger%print_keys%mo_eigenvalues.OR.&
               (scf_control%level_shift /= 0.0_dp).OR.&
               (scf_control%smear /= 0.0_dp).OR.&
               (scf_control%eps_jacobi /= 0.0_dp)) THEN
             n_mo = n_ao
          ELSE
             n_mo = nelectron/2
          END IF
          ALLOCATE(sub_mos(dft_control%nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             DO ispin=1,dft_control%nspins
                NULLIFY(sub_mos(ispin)%mo_set)
                IF (dft_control%nspins==1) THEN
                   CALL allocate_mo_set(sub_mos(ispin)%mo_set,&
                        nao=n_ao,nmo=n_mo,nelectron=nelectron,maxocc=2.0_dp,&
                        error=error)
                ELSE
                   CALL allocate_mo_set(sub_mos(ispin)%mo_set,&
                        nao=n_ao,nmo=n_mo,nelectron=nelectron,maxocc=1.0_dp,&
                        error=error)
                END IF
             END DO
          END IF
          CALL set_qs_env(sub_qs, mos=sub_mos, error=error)

          !   *** Allocate the data structure for Quickstep energies ***
          CALL allocate_qs_energy(sub_qs%energy)
          ! allocate qs_charges
          CALL qs_charges_create(sub_qs_charges,nspins=dft_control%nspins,error=error)
          CALL set_qs_env(sub_qs, qs_charges=sub_qs_charges,error=error)
          CALL qs_charges_release(sub_qs_charges, error=error)
          ! rebuilds ao_mo_fm_pools
          CALL mpools_rebuild_fm_pools(sub_qs%mpools, mos=sub_qs%mos, &
               blacs_env=sub_qs%blacs_env, para_env=sub_qs%para_env,&
               nrow_block=scf_control%nrow_block,&
               ncol_block=scf_control%ncol_block,error=error)
          ! init green's function
          CALL get_qs_env(sub_qs, cell=sub_cell, green=green, pw_env=pw_env,&
               rho_core=sub_rho_core, error=error)
          CALL pw_green_fn_rebuild(green,&
               pw_env%pw_pools(pw_env%auxbas_grid)%pool)
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
       
       ! create the perturbation environment
       CALL p_env_create(ep_env%main_p_env, qs_env=ep_env%main_qs_env,&
            orthogonal_orbitals=.FALSE.,error=error)

       ! alloc m_pi_Hrho_psi0d
       CALL get_qs_env(ep_env%main_qs_env, mpools=mpools,&
            error=error)
       CALL mpools_get(mpools,ao_mo_fm_pools=ao_mo_fm_pools, error=error)
       
       CALL fm_pools_create_fm_vect(ao_mo_fm_pools,ep_env%m_pi_Hrho_psi0d,&
            name="ep_env"//TRIM(ADJUSTL(cp_to_string(ep_env%id_nr)))//&
            "%m_pi_Hrho_psi0d",error=error)
       ! alloc psi1
       CALL fm_pools_create_fm_vect(ao_mo_fm_pools,ep_env%psi1,&
            name="ep_env"//TRIM(ADJUSTL(cp_to_string(ep_env%id_nr)))//"%psi1",&
            error=error)

       ALLOCATE(psi0(dft_control%nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          DO ispin=1,dft_control%nspins
             NULLIFY(psi0(ispin)%matrix)
             CALL get_mo_set(ep_env%main_qs_env%mos(ispin)%mo_set,&
                  mo_coeff=psi0(ispin)%matrix)
          END DO
       END IF
       ! initializes the preconditioner
       ALLOCATE(ep_env%precond(dft_control%nspins),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ispin=1
       NULLIFY(ep_env%precond(ispin)%preconditioner)
       ALLOCATE(ep_env%precond(ispin)%preconditioner,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CALL init_preconditioner(ep_env%precond(ispin)%preconditioner,&
            para_env=para_env,&
            blacs_env=blacs_env)
       CALL make_preconditioner(ep_env%precond(ispin)%preconditioner,&
            "FULL_KINETIC",&
            ep_env%main_qs_env%matrix_ks(ispin)%matrix,&
            ep_env%main_qs_env%matrix_s(1)%matrix,&
            ep_env%main_qs_env%kinetic(1)%matrix,&
            psi0(ispin)%matrix,0.2_dp)
       ! same precond for both spins
       DO ispin=2,dft_control%nspins
          ep_env%precond(ispin)%preconditioner => &
               ep_env%precond(1)%preconditioner
       END DO
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE ep_env_create
!***************************************************************************

!!****f* qs_ep_methods/ep_env_finish_qs_init [1.0] *
!!
!!   NAME
!!     ep_env_finish_qs_init
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
SUBROUTINE ep_env_finish_qs_init(qs_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_finish_qs_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s

  failure=.FALSE.
  NULLIFY(ks_env, mos, mo_coeff, dft_control,matrix_ks,matrix_s)

  CALL timeset(routineN,"I","",handle)

  CALL build_qs_neighbor_lists(qs_env,globenv)

! *** Calculate the overlap and the core Hamiltonian integral matrix ***

  CALL build_core_hamiltonian_matrix(qs_env=qs_env,&
                                     globenv=globenv,&
                                     calculate_forces=.FALSE.)
  CALL qs_env_update_s_mstruct(qs_env, error=error)

!  *** init the ks_env

  CALL get_qs_env(qs_env, ks_env=ks_env)
  IF (.not.ASSOCIATED(ks_env)) THEN
     CALL qs_ks_create(ks_env,qs_env=qs_env)
     CALL set_qs_env(qs_env, ks_env=ks_env)
  END IF

  !   *** Initializes the MOs ***
  CALL get_qs_env(qs_env,mos=mos)
  CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     DO ispin=1,SIZE(mos)
        CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
        IF (.NOT.ASSOCIATED(mo_coeff)) THEN
           CALL init_mo_set(mos(ispin)%mo_set, &
                qs_env%mpools%ao_mo_fm_pools(ispin)%pool,&
                name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                halfhole=-1,error=error)
        END IF
     END DO
  END IF

  !   *** Allocate k and put it in the QS environment ***
  CALL get_qs_env(qs_env,dft_control=dft_control,matrix_s=matrix_s)
  CALL allocate_matrix_set(matrix_ks,dft_control%nspins)
  DO ispin=1,dft_control%nspins
     NULLIFY (matrix_ks(ispin)%matrix)
     CALL replicate_matrix_structure(matrix_s(1)%matrix,matrix_ks(ispin)%matrix,&
          "KOHN-SHAM MATRIX-"//ADJUSTL(cp_to_string(ispin)))
  END DO
  CALL set_qs_env(qs_env=qs_env,&
       matrix_ks=matrix_ks)

  CALL calculate_ecore_self(qs_env)
  CALL calculate_ecore_overlap(qs_env,globenv,.FALSE.)

  CALL timestop(0.0_dp,handle)
END SUBROUTINE ep_env_finish_qs_init
!***************************************************************************

!!****f* qs_ep_methods/ep_env_transfer_psi0 [1.0] *
!!
!!   NAME
!!     ep_env_transfer_psi0
!!
!!   FUNCTION
!!     transfers psi0 from the sub_qs_env to the main_qs_env
!!     (i.e. sub_qs_env%c are injected into main_qs_env%c)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
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
SUBROUTINE ep_env_transfer_psi0(ep_env, error)
    TYPE(qs_ep_env_type), POINTER            :: ep_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_transfer_psi0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, isub, nmo, &
                                                nspins, start_col, sub_nmo, &
                                                unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_fm_type), POINTER                :: orbitals, sub_orbitals
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, sub_mos
    TYPE(qs_environment_p_type), &
      DIMENSION(:), POINTER                  :: sub_qs_env
    TYPE(qs_environment_type), POINTER       :: main_qs_env
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), POINTER                  :: sub_proj

  CALL timeset(routineN,"I","",handle)

  failure=.FALSE.
  NULLIFY(sub_orbitals, sub_mos, main_qs_env, sub_qs_env, sub_proj, logger)
  logger => cp_error_get_logger(error)
  unit_nr=cp_logger_get_default_unit_nr(logger, local=.FALSE.)
  
  CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
  CALL ep_env_get(ep_env,sub_qs_env=sub_qs_env,&
       sub_proj=sub_proj,main_qs_env=main_qs_env,error=error)
  CALL cp_error_check(error,failure)
  CPPrecondition(ASSOCIATED(sub_qs_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(sub_proj),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(main_qs_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     nspins=main_qs_env%dft_control%nspins
     CALL get_qs_env(main_qs_env,mos=mos,error=error)
     DO ispin=1,nspins ! transfer also occupation_numbers ?
        start_col=1
        CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=orbitals,nmo=nmo)
        CALL cp_fm_set_all(orbitals,0.0_dp,error=error)

        DO isub=1,SIZE(sub_qs_env)
           CALL get_qs_env(sub_qs_env(isub)%qs_env,mos=sub_mos,error=error)
           CALL get_mo_set(sub_mos(ispin)%mo_set,mo_coeff=sub_orbitals,&
                nmo=sub_nmo)
           CPInvariant(start_col+sub_nmo<nmo+2,cp_failure_level,routineP,error,failure)
           IF (debug_this_module .AND. logger%print_keys%ep_matrixes) THEN
              WRITE (unit=unit_nr,fmt="(a,i4,a,i4,a,i4)")&
                   "Printing psi of sub_env ",isub,&
                   " spin=",ispin," nmo=",sub_nmo
              CALL cp_fm_write(sub_orbitals,unit_nr,&
                   long_description=.TRUE.,&
                   local=.FALSE.,error=error)
           END IF

           CALL p_proj_transfer_blacs_to_f(p_proj=sub_proj(isub)%projection,&
                min_m=sub_orbitals, full_m=orbitals,&
                start_col_full=start_col, ncol=sub_nmo, error=error)

           start_col=start_col+sub_nmo
        END DO
        CPPostcondition(start_col==nmo+1,cp_failure_level,routineP,error,failure)

        IF (logger%print_keys%ep_matrixes) THEN
           WRITE (unit=unit_nr,fmt="(a,i4,a,i4)")&
                "Printing psi of main_qs_env spin=",ispin," nmo=",nmo
           CALL cp_fm_write(orbitals,unit_nr,long_description=.TRUE.,&
                local=.FALSE.,error=error)
        END IF
     END DO
  END IF
END SUBROUTINE ep_env_transfer_psi0
!***************************************************************************

!!****f* qs_ep_methods/stupid_solve [1.0] *
!!
!!   NAME
!!     stupid_solve
!!
!!   FUNCTION
!!     solves a linear system with conjugated gradient
!!
!!   NOTES
!!     The linear solver should go in nother module, this is a quick and 
!!     dirty implementation.
!!     A: postortho (p_op_l1+p_op_l2) preortho
!!     x: ep_env%psi1
!!     ugly scaling of p (to avoid cutoff os small values on the grids)
!!
!!   ARGUMENTS
!!     - ep_env: the sytem to solve
!!     - eps_r: wanted residual error
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE stupid_solve(ep_env,eps_r,error)
    TYPE(qs_ep_env_type), POINTER            :: ep_env
    REAL(KIND=dp), INTENT(in)                :: eps_r
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'stupid_solve', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ispin, j, nspins
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: alpha, Ap_p, Ap_p_s, beta, p_ortho_norm, p_ortho_norm_s, &
      r_norm2, r_norm2_s, r_z, r_z_new, r_z_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: Ap, p, p_ortho, r, z
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

  failure=.FALSE.
  logger => cp_error_get_logger(error)
  NULLIFY(p,p_ortho,Ap,r,z,ao_mo_fm_pools,dft_control,matrix_s,mpools)

  CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL get_qs_env(ep_env%main_qs_env, &
          dft_control=dft_control,matrix_s=matrix_s,mpools=mpools,error=error)
     CALL mpools_get(mpools,ao_mo_fm_pools=ao_mo_fm_pools, error=error)
     nspins=dft_control%nspins

     CALL fm_pools_create_fm_vect(ao_mo_fm_pools,p,name="ep_CG_p",&
          error=error)
     CALL fm_pools_create_fm_vect(ao_mo_fm_pools,p_ortho,&
          name="ep_CG_p_ortho",&
          error=error)
     CALL fm_pools_create_fm_vect(ao_mo_fm_pools,Ap,name="ep_CG_Ap",&
          error=error)
     CALL fm_pools_create_fm_vect(ao_mo_fm_pools,r,name="ep_CG_r",&
          error=error)
     CALL fm_pools_create_fm_vect(ao_mo_fm_pools,z,name="ep_CG_z",&
          error=error)
     CALL cp_error_check(error,failure)
  END IF
  IF (.not.failure) THEN
     DO i=1,20
        ! orthogonalize x wrt. to psi0
        CALL p_preortho(p_env=ep_env%main_p_env, &
             qs_env=ep_env%main_qs_env, v=ep_env%psi1,&
             n_cols=ep_env%main_p_env%n_mo, error=error)
        ! localize x
        CALL ep_env_localize_matrix(ep_env, ep_env%psi1, error=error)
        ! r=b-A x
        CALL p_op_ep(p_env=ep_env%main_p_env, qs_env=ep_env%main_qs_env,&
             v=ep_env%psi1, res=r, error=error)
        DO ispin=1,nspins
           CALL cp_fm_scale_and_add(alpha=-1.0_dp,matrix_a=r(ispin)%matrix,&
                beta=1.0_dp,matrix_b=ep_env%m_pi_Hrho_psi0d(ispin)%matrix,&
                error=error)
        END DO
        ! localize r
        CALL ep_env_localize_matrix(ep_env, r, error=error)
!FM        CALL cp_fm_vect_write(matrixes=r, &
!FM             unit_nr=cp_logger_get_default_unit_nr(logger),&
!FM             long_description=.TRUE., error=error) !pippo
        
        ! check convergence (with ||r||)
        r_norm2=0.0_dp
        DO ispin=1,nspins
           CALL cp_fm_trace(matrix_a=r(ispin)%matrix,&
                matrix_b=r(ispin)%matrix,trace=r_norm2_s,error=error)
           r_norm2=r_norm2+r_norm2_s
        END DO
        IF (r_norm2<eps_r**2) EXIT
        
        ! z=M^-1 r
        DO ispin=1,nspins
           CALL apply_preconditioner(ep_env%precond(ispin)%preconditioner,&
                matrix_in=r(ispin)%matrix,&
                matrix_out=z(ispin)%matrix)
        END DO
        ! localize z
        CALL ep_env_localize_matrix(ep_env, z, error=error)

        ! p=z
        DO ispin=1,nspins
           CALL cp_fm_scale_and_add(alpha=0.0_dp,matrix_a=p(ispin)%matrix,beta=1.0_dp,&
                matrix_b=z(ispin)%matrix,error=error)
        END DO
        r_z=0.0_dp
        DO ispin=1,nspins
           CALL cp_fm_trace(matrix_a=r(ispin)%matrix,&
                matrix_b=z(ispin)%matrix,trace=r_z_s,error=error)
           r_z=r_z+r_z_s
        END DO

        DO j=1,10
           CALL cp_log(logger,cp_note_level,routineP,"i="//cp_to_string(i)//&
                "j="//cp_to_string(j)//"r_norm2="//cp_to_string(r_norm2))
           
           ! put the component of p orthogonal to psi0 in p_ortho
           DO ispin=1,nspins
              CALL cp_fm_scale_and_add(alpha=0.0_dp, matrix_a=p_ortho(ispin)%matrix,&
                   beta=1.0_dp,&
                   matrix_b=p(ispin)%matrix,error=error)
           END DO
           CALL p_preortho(p_env=ep_env%main_p_env,&
             qs_env=ep_env%main_qs_env, v=p_ortho,&
             n_cols=ep_env%main_p_env%n_mo, error=error)
           ! localize p_ortho
           CALL ep_env_localize_matrix(ep_env, p_ortho, error=error)

           ! calc norm of p_ortho
           p_ortho_norm=0.0_dp
           DO ispin=1,nspins
              CALL cp_sm_fm_multiply(sparse_matrix=matrix_s(1)%matrix,&
                   v_in=p_ortho(ispin)%matrix,v_out=Ap(ispin)%matrix,&
                   ncol=ep_env%main_p_env%n_mo(ispin),&
                   error=error)
              CALL cp_fm_trace(matrix_a=p_ortho(ispin)%matrix,&
                   matrix_b=Ap(ispin)%matrix,&
                   trace=p_ortho_norm_s,error=error)
              p_ortho_norm=p_ortho_norm+p_ortho_norm_s
           END DO
           p_ortho_norm=SQRT(p_ortho_norm)
!FM           CALL cp_fm_vect_write(matrixes=p_ortho, &
!FM                unit_nr=cp_logger_get_default_unit_nr(logger),&
!FM                long_description=.TRUE., error=error) !pippo

           ! rescale p_ortho
           CPPreconditionNoFail(p_ortho_norm>eps_r*0.5,cp_warning_level,routineP,error)
           DO ispin=1,nspins
              CALL cp_fm_scale_and_add(alpha=1.0_dp/p_ortho_norm,&
                   matrix_a=p_ortho(ispin)%matrix,error=error)
           END DO

           ! Ap=A p_ortho=p_ortho_norm A p
           CALL p_op_ep(p_env=ep_env%main_p_env, qs_env=ep_env%main_qs_env,&
                v=p_ortho, res=Ap,error=error)
           ! localize Ap
           CALL ep_env_localize_matrix(ep_env, Ap, error=error)

!FM           CALL cp_fm_vect_write(matrixes=Ap, &
!FM                unit_nr=cp_logger_get_default_unit_nr(logger),&
!FM                long_description=.TRUE., error=error) !pippo
           
           Ap_p=0.0_dp
           DO ispin=1,nspins
              CALL cp_fm_trace(matrix_a=Ap(ispin)%matrix,&
                   matrix_b=p(ispin)%matrix,trace=Ap_p_s,error=error)
              Ap_p=Ap_p+Ap_p_s
           END DO
           
!FM           ! check if A is almost indefinite (in this case Conjugate 
!FM           ! Residual would be better)
!FM           CALL cp_assert(p_ortho_norm*Ap_p>r_norm2*MIN(eps_r,1.0e-6_dp),&
!FM                cp_warning_level,cp_assertion_failed,routineP,&
!FM                "A is almost indefinite (Conjugate Residual would be better)"//&
!FM                CPSourceFileRef,&
!FM                error=error)
           
           ! this alpha is the alpha in alg. desc. times p_ortho_norm
           alpha=r_z/Ap_p

           ! x=x+alpha p_ortho
           DO ispin=1,nspins
              CALL cp_fm_scale_and_add(alpha=1.0_dp,matrix_a=ep_env%psi1(ispin)%matrix,&
                   beta=alpha,matrix_b=p_ortho(ispin)%matrix,&
                   error=error)
           END DO
!FM        CALL cp_fm_vect_write(matrixes=ep_env%psi1, &
!FM             unit_nr=cp_logger_get_default_unit_nr(logger),&
!FM             long_description=.TRUE., error=error) !pippo
           
           ! r=r- alpha Ap
           DO ispin=1,nspins
              CALL cp_fm_scale_and_add(alpha=1.0_dp,matrix_a=r(ispin)%matrix,&
                   beta=-alpha,matrix_b=Ap(ispin)%matrix,&
                   error=error)
           END DO
           
           r_norm2=0.0_dp
           DO ispin=1,nspins
              CALL cp_fm_trace(matrix_a=r(ispin)%matrix,&
                   matrix_b=r(ispin)%matrix,trace=r_norm2_s,error=error)
              r_norm2=r_norm2+r_norm2_s
           END DO
           IF (r_norm2<eps_r**2) EXIT
           
           ! z=M^-1 r
           DO ispin=1,nspins
              CALL apply_preconditioner(ep_env%precond(ispin)%preconditioner,&
                   matrix_in=r(ispin)%matrix,&
                   matrix_out=z(ispin)%matrix)
           END DO
           ! localize z
           CALL ep_env_localize_matrix(ep_env, z, error=error)

           r_z_new=0.0_dp
           DO ispin=1,nspins
              CALL cp_fm_trace(matrix_a=r(ispin)%matrix,&
                   matrix_b=z(ispin)%matrix,trace=r_z_s,error=error)
              r_z_new=r_z_new+r_z_s
           END DO

           beta= r_z_new/r_z
           r_z=r_z_new
           
           ! p=beta p+z
           DO ispin=1,nspins
              CALL cp_fm_scale_and_add(alpha=beta,matrix_a=p(ispin)%matrix,&
                   beta=1.0_dp,matrix_b=z(ispin)%matrix,&
                   error=error)
           END DO
        END DO
        CALL cp_log(logger,cp_note_level,routineP,"i="//cp_to_string(i)//&
             "j=final r_norm2="//cp_to_string(r_norm2))
     END DO
     
     CALL cp_log(logger,cp_note_level,routineP,&
          "solve finished, r_norm2="//cp_to_string(r_norm2))
     CPPostcondition(r_norm2<eps_r**2,cp_warning_level,routineP,error,failure)
     
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,p,error=error)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,Ap,error=error)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,r,error=error)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,z,error=error)
  END IF
           
END SUBROUTINE stupid_solve
!***************************************************************************

!!****f* qs_ep_methods/p_op_ep [1.0] *
!!
!!   NAME
!!     p_op_ep
!!
!!   FUNCTION
!!     evaluates the ep operator (p_op_l1+p_op_l2) on the given vector
!!
!!   NOTES
!!     no preorthogonalization, thus v must be orthogonal to the psi0
!!
!!   ARGUMENTS
!!     - p_env : perturbation calculation environment
!!     - qs_env: the qs_env that is perturbed by this p_env
!!     - v     : the matrix to operate on (must be orthogonal to the psi0)
!!     - res   : the result
!!     - error : error handling object (optional)
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE p_op_ep(p_env, qs_env, v, res, error)
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(in)                             :: v
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(inout)                          :: res
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_op_ep', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, lfomo, nmo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: maxocc, scale_factor
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

  CALL timeset(routineN,"I","",handle)

  failure=.FALSE.
  NULLIFY(mos)

  CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL get_qs_env(qs_env,mos=mos,error=error)
     CALL p_op_l1(p_env=p_env, qs_env=qs_env,&
          v=v, res=res,error=error)
     CALL get_mo_set(mos(1)%mo_set,lfomo=lfomo, nmo=nmo, maxocc=maxocc)
     IF (lfomo>nmo) THEN
        scale_factor=maxocc
        DO ispin=2,SIZE(mos)
           CALL get_mo_set(mos(ispin)%mo_set,lfomo=lfomo, &
                nmo=nmo, maxocc=maxocc)
           CPAssert(scale_factor==maxocc,cp_failure_level,routineP,error,failure)
        END DO
        CALL p_op_l2_fawzi(p_env=p_env, qs_env=qs_env,&
             v=v, res=res,&
             alpha=scale_factor,beta=scale_factor,error=error)
     ELSE
        CALL cp_unimplemented_error(fromWhere=routineP,&
             message="symmetrized onesided smearing to do",&
             error=error)
     END IF
     CALL p_postortho(p_env=p_env, qs_env=qs_env, v=res, n_cols=p_env%n_mo,&
          error=error)
  END IF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE p_op_ep
!***************************************************************************

!!****f* qs_ep_methods/p_env_write_ep_matrix [1.0] *
!!
!!   NAME
!!     p_env_write_ep_matrix
!!
!!   FUNCTION
!!     writes out the whole ep interaction matrix (as defined by p_op_ep)
!!     one column per file
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - p_env: the p_env where to evaluate the interaction matrix
!!     - qs_env: the qs_env that is perturbed by the p_env
!!     - name: name used for the matrixes outputted
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE p_env_write_ep_matrix(p_env, qs_env, name, error)
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    CHARACTER(len=*), INTENT(in)             :: name
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_write_ep_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iao, imo, ispin, ispin2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: v_norm, v_norm_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: res, v
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

  failure=.FALSE.
  logger => cp_error_get_logger(error)
  NULLIFY(ao_mo_fm_pools, dft_control,matrix_s)
  
  CPPrecondition(ASSOCIATED(p_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(p_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL get_qs_env(qs_env,&
          dft_control=dft_control,matrix_s=matrix_s,error=error)
     CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
          error=error)
     CALL fm_pools_create_fm_vect(ao_mo_fm_pools,v,error=error)
     CALL fm_pools_create_fm_vect(ao_mo_fm_pools,res,name=name,&
          error=error)
     
     DO ispin=1,dft_control%nspins
        DO imo=1,p_env%n_mo(ispin)
           DO iao=1,p_env%n_ao(ispin)
              CALL cp_fm_vect_set_all(v,0.0_dp,error=error)
              CALL cp_fm_set_element(v(ispin)%matrix,iao,imo,1.0_dp,&
                   error=error)

              CALL p_preortho(p_env=p_env, qs_env=qs_env, v=v,&
                   n_cols=p_env%n_mo, error=error)
              ! scale v
              v_norm=0.0_dp
              DO ispin2=1,dft_control%nspins
                 CALL cp_sm_fm_multiply(sparse_matrix=matrix_s(1)%matrix,&
                      v_in=v(ispin2)%matrix,v_out=res(ispin2)%matrix,&
                      ncol=p_env%n_mo(ispin2), &
                      error=error)
                 CALL cp_fm_trace(matrix_a=v(ispin2)%matrix, &
                      matrix_b=res(ispin2)%matrix, trace=v_norm_s, error=error)
                 v_norm=v_norm+v_norm_s
              END DO
              v_norm=SQRT(v_norm)
              CPPreconditionNoFail(v_norm>1.0e-8,cp_warning_level,routineP,error)
              DO ispin2=1,dft_control%nspins
                 CALL cp_fm_scale_and_add(1.0_dp/v_norm,v(ispin2)%matrix,error=error)
              END DO

              CALL p_op_ep(p_env=p_env, qs_env=qs_env, v=v,&
                   res=res, error=error)

              DO ispin2=1,dft_control%nspins
                 CALL cp_fm_scale_and_add(v_norm,res(ispin2)%matrix,error=error)
              END DO
              
              CALL cp_fm_vect_write(matrixes=res, &
                   unit_nr=cp_logger_get_default_unit_nr(logger),&
                   long_description=.TRUE., error=error)
              CALL m_flush(cp_logger_get_default_unit_nr(logger))
           END DO
        END DO
     END DO

     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, v, error=error)
     CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools, res, error=error)
  END IF
END SUBROUTINE p_env_write_ep_matrix
!***************************************************************************

!!****f* qs_ep_methods/ep_env_localize_matrix [1.0] *
!!
!!   NAME
!!     ep_env_localize_matrix
!!
!!   FUNCTION
!!     localizes the wavefunction matrix, restraining it the the subsystem
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE ep_env_localize_matrix(ep_env, wf_coeff, error)
    TYPE(qs_ep_env_type), POINTER            :: ep_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: wf_coeff
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_localize_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, isub, nmo, &
                                                nspins, start_col, sub_nmo, &
                                                unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: sub_mos
    TYPE(qs_environment_p_type), &
      DIMENSION(:), POINTER                  :: sub_qs_env
    TYPE(qs_environment_type), POINTER       :: main_qs_env
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), POINTER                  :: sub_proj

  CALL timeset(routineN,"I","",handle)

  failure=.FALSE.
  NULLIFY(main_qs_env,sub_qs_env, sub_proj, sub_mos)
  logger => cp_error_get_logger(error)
  unit_nr=cp_logger_get_default_unit_nr(logger, local=.FALSE.)
  
  CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
  CALL ep_env_get(ep_env,sub_qs_env=sub_qs_env,&
       sub_proj=sub_proj,main_qs_env=main_qs_env,error=error)
  CALL cp_error_check(error,failure)
  CPPrecondition(ASSOCIATED(sub_qs_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(sub_proj),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(main_qs_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     nspins=main_qs_env%dft_control%nspins
     DO ispin=1,nspins
        start_col=1
        CALL cp_fm_get_info(wf_coeff(ispin)%matrix,ncol_global=nmo)

        DO isub=1,SIZE(sub_qs_env)
           CALL get_qs_env(sub_qs_env(isub)%qs_env,mos=sub_mos,error=error)
           CALL get_mo_set(sub_mos(ispin)%mo_set, nmo=sub_nmo)
           CPInvariant(start_col+sub_nmo<nmo+2,cp_failure_level,routineP,error,failure)

!FM           CALL p_proj_restrain_f(p_proj=sub_proj(isub)%projection,&
!FM                full_m=wf_coeff(ispin)%matrix,&
!FM                start_col=start_col, ncol=sub_nmo, error=error)
           
           start_col=start_col+sub_nmo
        END DO
        CPPostcondition(start_col==nmo+1,cp_failure_level,routineP,error,failure)
     END DO
  END IF
END SUBROUTINE ep_env_localize_matrix
!***************************************************************************

END MODULE qs_ep_methods
