!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief in this module there are the routines for the calculation of
!>      psi and energy without scf, but perturbatively as in PRL 87(22);226401;
!>      26 Nov 2001,"Accurate Total Energies without Self-Consistency",
!>      David M. Benoit, Daniel Sebastiani, and Michele Parrinello
!> \par History
!>      11.2002 created [fawzi]
!>      06.2006 ressurected, renamed qs_ep_methods => ep_methods
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE ep_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_scale,&
                                             cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_symm,&
                                             cp_fm_trace
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_give_back_fm,&
                                             fm_pools_create_fm_vect,&
                                             fm_pools_give_back_fm_vect
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type,&
                                             cp_fm_write
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc,&
                                             cp_fm_vect_set_all,&
                                             cp_fm_vect_write
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: &
       cp_add_iter_level, cp_iter_string, cp_iterate, cp_p_file, cp_p_store, &
       cp_print_key_finished_output, cp_print_key_log, &
       cp_print_key_should_output, cp_print_key_unit_nr, cp_rm_iter_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE ep_qs_methods,                   ONLY: ep_qs_create
  USE ep_qs_types,                     ONLY: ep_qs_release,&
                                             ep_qs_type
  USE ep_types,                        ONLY: ep_energy_type,&
                                             ep_env_p_type,&
                                             ep_env_type,&
                                             ep_force_type
  USE f77_interface,                   ONLY: calc_energy,&
                                             calc_energy_force,&
                                             create_force_env,&
                                             destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type,&
                                             get_natom
  USE force_env_types,                 ONLY: force_env_get
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_release,&
                                             globenv_retain
  USE input_constants,                 ONLY: do_ep,&
                                             do_qs,&
                                             ot_precond_full_kinetic,&
                                             ot_precond_solver_default,&
                                             wfi_use_prev_wf_method_nr
  USE input_cp2k_read,                 ONLY: empty_initial_variables,&
                                             read_input
  USE input_section_types,             ONLY: &
       section_get_ivals, section_type, section_vals_get_subs_vals, &
       section_vals_release, section_vals_retain, section_vals_type, &
       section_vals_val_get, section_vals_val_set, section_vals_write
  USE iso_c_binding
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE preconditioner,                  ONLY: apply_preconditioner,&
                                             make_preconditioner
  USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             preconditioner_p_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env,&
                                             qs_env_update_s_mstruct
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_env_retain,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: write_forces
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type,&
                                             zero_qs_force
  USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                             qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             set_ks_env
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE qs_p_env_methods,                ONLY: p_env_create,&
                                             p_env_did_change,&
                                             p_env_psi0_changed,&
                                             p_op_l1,&
                                             p_op_l2_fawzi,&
                                             p_postortho,&
                                             p_preortho
  USE qs_p_env_types,                  ONLY: p_env_release,&
                                             qs_p_env_type
  USE qs_p_sparse_psi,                 ONLY: p_proj_create,&
                                             p_proj_release,&
                                             qs_p_projection_p_type,&
                                             qs_p_projection_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_set,&
                                             qs_subsys_type
  USE qs_wf_history_methods,           ONLY: wfi_change_memory_depth,&
                                             wfi_create,&
                                             wfs_duplicate_snapshot
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             qs_wf_snapshot_type,&
                                             wfi_get_snapshot,&
                                             wfi_release,&
                                             wfs_release
  USE replica_methods,                 ONLY: rep_env_calc_e_f,&
                                             rep_env_create
  USE replica_types,                   ONLY: rep_env_local_index,&
                                             rep_env_release,&
                                             replica_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ep_methods'
  INTEGER, SAVE :: last_ep_env_id=0, last_force_id=0
  LOGICAL, SAVE, PRIVATE :: module_initialized=.FALSE.

  TYPE(ep_env_p_type), DIMENSION(:), POINTER, PRIVATE :: ep_envs

  PUBLIC :: ep_env_create

  PUBLIC :: ep_env_retain, ep_env_release
  PUBLIC :: ep_envs_get_ep_env, ep_env_calc_e_f
CONTAINS

! *****************************************************************************
!> \brief writes out the qs energies
!> \param qs_env ...
!> \param output_unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE print_qs_energies(qs_env,output_unit,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_qs_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho

    failure=.FALSE.

    CALL timeset(routineN,handle)
    IF (output_unit>0) THEN
       NULLIFY(rho,energy,qs_charges,dft_control)
       CALL get_qs_env(qs_env,rho=rho,energy=energy,dft_control=dft_control,&
            qs_charges=qs_charges,error=error)
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T61,F20.10))")&
            "Total electronic density (r-space): ",&
            accurate_sum(rho%tot_rho_r),&
            "Total core charge density (r-space):",&
            qs_charges%total_rho_core_rspace
       WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
            "Total charge density (r-space):     ",&
            accurate_sum(rho%tot_rho_r)+&
            qs_charges%total_rho_core_rspace,&
            "Total charge density (g-space):     ",qs_charges%total_rho_gspace
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T56,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
       IF (energy%e_hartree /= 0.0_dp) &
            WRITE (UNIT=output_unit,FMT="(T3,A,/,T3,A,T56,F25.14)")&
            "Coulomb Electron-Electron Interaction Energy ",&
            "- Already included in the total Hartree term ",energy%e_hartree
       IF (dft_control%qs_control%mulliken_restraint) THEN
          WRITE (UNIT=output_unit,FMT="(T3,A,T61,F20.10)")&
               "Mulliken restraint energy: ",energy%mulliken
       END IF
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T56,F25.14))")&
            "Total energy:                                  ",energy%total
       CALL m_flush(output_unit)
    END IF
    CALL timestop(handle)
  END SUBROUTINE print_qs_energies

! *****************************************************************************
!> \brief calculates c0 for a single molecule
!> \param ep_env the ep_env with all the
!> \param mol_pos ...
!> \param new_c ...
!> \param imol ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Clotilde Cucinotta
! *****************************************************************************
  SUBROUTINE calc_c0_tilde(ep_env,mol_pos,new_c,imol,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    REAL(kind=dp), DIMENSION(:), INTENT(in)  :: mol_pos
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      INTENT(inout)                          :: new_c
    INTEGER, INTENT(in)                      :: imol
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: icol, irow, ispin
    REAL(dp), DIMENSION(:), POINTER          :: coeff_standard1d

! legge dall'input

    CALL section_vals_val_get(ep_env%root_section,"FORCE_EVAL%EP%START_COEFFS",&
         r_vals=coeff_standard1d,error=error)

    ! trova rotazione usando la posizione mol_pos (coordinate xyz xyz xyz)
    ! che porta la molecola nella posizione standard

    ! ruota gli orbitali standard per avere gli orbitali della molecola
    ! attuale in new_c
    DO ispin=1,ep_env%nspins
       DO icol=1,ep_env%sub_nmo(ispin)
          DO irow=1,ep_env%sub_nao(ispin)
             new_c(ispin)%array(icol,irow)=0.0_dp
          END DO
       END DO
    END DO
  END SUBROUTINE calc_c0_tilde


! *****************************************************************************
!> \brief creates the environement to perform an ep perturbation
!> \param ep_env the ep_env env to create
!> \param root_section ...
!> \param input_declaration ...
!> \param para_env ...
!> \param globenv ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      Could take more advantage of what is calculated in main_qs_env
!>      (for example the S matrix)
! *****************************************************************************
  SUBROUTINE ep_env_create(ep_env,root_section,input_declaration,para_env,globenv, error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(.NOT.ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       ALLOCATE(ep_env,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

       last_ep_env_id=last_ep_env_id+1
       ep_env%id_nr=last_ep_env_id
       ep_env%ref_count=1
       ep_env%f_env_id=-1
       ep_env%nspins=-1
       ep_env%nat=-1
       ep_env%nat_per_mol=-1
       ep_env%nmol=-1
       NULLIFY(ep_env%mol_envs,ep_env%sub_proj,ep_env%main_qs_env,&
            ep_env%main_p_env,ep_env%sub_p_env,ep_env%m_pi_Hrho_psi0d,&
            ep_env%psi1, ep_env%precond,ep_env%sub_nmo,ep_env%sub_nao,&
            ep_env%full_nmo,ep_env%full_nao,ep_env%input,ep_env%at2sub,&
            ep_env%force,ep_env%root_section,ep_env%base_C0)
       CALL ep_energy_zero(ep_env%energy,error=error)
       ep_env%root_section => root_section
       CALL section_vals_retain(root_section,error)
       ep_env%globenv => globenv
       CALL globenv_retain(globenv,error=error)
       ep_env%para_env=>para_env
       CALL cp_para_env_retain(ep_env%para_env,error=error)

       CALL ep_envs_add_ep_env(ep_env,error=error)
       CALL ep_env_init_low(ep_env%id_nr,input_declaration,ierr)

       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE ep_env_create


! *****************************************************************************
!> \brief creates the environement to perform an ep perturbation
!> \param ep_env_id ...
!> \param input_declaration ...
!> \param ierr ...
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      Could take more advantage of what is calculated in main_qs_env
!>      (for example the S matrix)
! *****************************************************************************
  SUBROUTINE ep_env_init_low(ep_env_id,input_declaration,ierr)
    INTEGER, INTENT(in)                      :: ep_env_id
    TYPE(section_type), POINTER              :: input_declaration
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_init_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: comp_path
    CHARACTER(LEN=default_string_length)     :: project_name
    INTEGER                                  :: at_per_mol, handle, i, iat, &
                                                imol, ispin, nat, nmol, stat, &
                                                unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atoms
    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_error_type)                      :: error, new_error
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(ep_qs_type), POINTER                :: ep_qs_env
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: full_mos, sub_mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_environment_type), POINTER       :: sub_qs_env
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), POINTER                  :: sub_proj
    TYPE(section_vals_type), POINTER         :: comp_input, input

!    TYPE(qs_p_projection_p_type), &
!      DIMENSION(:), INTENT(in)               :: projections
!    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL timeset(routineN,handle)
    failure=.FALSE.
    ep_env => ep_envs_get_ep_env(ep_env_id)
    NULLIFY(particle_set, qs_kind_set)
    CALL cp_error_init(error)
    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       input => ep_env%root_section
       para_env => ep_env%para_env
       logger => cp_error_get_logger(error)
       CALL cp_log(logger, level=cp_note_level, fromWhere=routineP ,&
            message=routineP//" creating main_qs_env", local=.FALSE.)

       CALL cp_print_key_log(logger,input,"FORCE_EVAL%EP%PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//" creating main_qs_env",error=error)
       ! create_main_qs_env
       CALL section_vals_val_set(input,"FORCE_EVAL%METHOD",&
            i_val=do_qs,error=error)
       CALL section_vals_val_get(input,"GLOBAL%PROJECT_NAME",&
            c_val=project_name,error=error)
       CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
            c_val=TRIM(project_name)//"_mainQS",error=error)
       CALL open_file(file_name=TRIM(project_name)//"_mainQS.inp",&
            file_status="UNKNOWN",file_form="FORMATTED",file_action="WRITE",&
            unit_number=unit_nr)
       CALL section_vals_write(input,unit_nr=unit_nr,hide_root=.TRUE.,&
            error=error)
       CALL close_file(unit_nr)
       CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",&
            c_val=project_name,error=error)
       CALL section_vals_val_set(input,"FORCE_EVAL%METHOD",&
            i_val=do_ep,error=error)

       CALL create_force_env(ep_env%f_env_id,input_declaration,&
            TRIM(project_name)//"_mainQS.inp",&
            TRIM(project_name)//"_mainQS.out",para_env%group,ierr=ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       globenv => ep_env%globenv
       NULLIFY(f_env)
       CALL f_env_add_defaults(f_env_id=ep_env%f_env_id,f_env=f_env,&
            new_error=new_error, failure=failure)
       CALL force_env_get(f_env%force_env,qs_env=ep_env%main_qs_env,&
            globenv=globenv,error=new_error)
       CALL qs_env_retain(ep_env%main_qs_env,error=error)
       CALL cp_print_key_log(logger,input,"FORCE_EVAL%EP%PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//" created main_qs_env "//&
            TRIM(ADJUSTR(cp_to_string(ep_env%main_qs_env%id_nr))),error=error)

       !-       ! reduce memory usage, get rid of an input structure
       !-       CALL section_vals_retain(root_section,error=error)
       !-       CALL section_vals_release(ep_env%root_section,error=error)
       !-       ep_env%root_section => root_section
       !-       input => ep_env%root_section
       CALL section_vals_release(ep_env%input,error=error)
       ep_env%input => section_vals_get_subs_vals(ep_env%root_section,&
            "FORCE_EVAL%EP",error=error)
       CALL section_vals_retain(ep_env%input,error=error)
       !-       CALL globenv_retain(globenv,error=error)
       !-       CALL globenv_release(ep_env%globenv,error=error)
       !-       ep_env%globenv => globenv

       ! *** initial initializations (atom position dependent)
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//" initial setup of main_qs_env",&
            error=error)

       CALL ep_env_finish_qs_init(ep_env%main_qs_env, globenv, error)

       CALL get_qs_env(ep_env%main_qs_env,dft_control=dft_control,&
            blacs_env=blacs_env,para_env=para_env,error=error)

       CALL cp_para_env_retain(para_env,error=error)
       CALL cp_para_env_release(ep_env%para_env,error=error)
       ep_env%para_env => para_env

       ep_env%nspins=dft_control%nspins
       ALLOCATE(ep_env%full_nao(ep_env%nspins),ep_env%full_nmo(ep_env%nspins),&
            stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       NULLIFY(full_mos)
       CALL get_qs_env(ep_env%main_qs_env,mos=full_mos, particle_set=particle_set,&
                       qs_kind_set=qs_kind_set, mpools=mpools,error=new_error)
       DO ispin=1,ep_env%nspins
          CALL get_mo_set(full_mos(ispin)%mo_set,nao=ep_env%full_nao(ispin),&
               nmo=ep_env%full_nmo(ispin))
       END DO

       ! add ep_qs_env
       NULLIFY(ep_qs_env)
       CALL ep_qs_create(ep_qs_env,ep_env%main_qs_env,error=error)
!       CALL set_qs_env(ep_env%main_qs_env,ep_qs_env=ep_qs_env,error=error)
       CALL ep_qs_release(ep_qs_env,error=error)
       !
       CALL f_env_rm_defaults(f_env,new_error,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " initial setup of main_qs_env finished",error=error)

       ! create mol_envs
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " creating mol_envs",error=error)
       CALL get_natom(ep_env%f_env_id,nat,ierr)
       CALL section_vals_val_get(input,"FORCE_EVAL%EP%AT_PER_MOL",&
            i_val=at_per_mol,error=error)
       nmol=nat/at_per_mol
       CPPrecondition(MODULO(nat,at_per_mol)==0,cp_failure_level,routineP,error,failure)
       ep_env%nmol=nmol
       ep_env%nat=nat
       ep_env%nat_per_mol=at_per_mol
       CALL ep_force_create(ep_env%force,nat=nat,error=error)
       CALL section_vals_val_get(input,"FORCE_EVAL%EP%COMP_INPUT",&
            c_val=comp_path,error=error)
       comp_input => read_input(input_declaration, comp_path, &
                                initial_variables=empty_initial_variables, &
                                para_env=para_env,error=error)
       CALL rep_env_create(ep_env%mol_envs, para_env=para_env, input=comp_input,&
            input_declaration=input_declaration,nrep=nmol, prep=1,keep_wf_history=.TRUE.,error=error)
       CALL cp_assert(ASSOCIATED(ep_env%mol_envs),cp_failure_level,&
            cp_assertion_failed,routineP,"ep not implemented for a number of"//&
            " processors bigger than the number of molecules",error,failure)
       CALL section_vals_release(comp_input,error=error)
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " mol_envs created",error=error)

       !sub_nao,sub_nmo
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " mol_envs initial setup",error=error)
       ALLOCATE(ep_env%sub_nao(ep_env%nspins),ep_env%sub_nmo(ep_env%nspins),&
            stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       ep_env%sub_nao=0
       ep_env%sub_nmo=0
       IF (ASSOCIATED(ep_env%mol_envs)) THEN
          NULLIFY(sub_mos,sub_qs_env)
          CALL f_env_add_defaults(f_env_id=ep_env%mol_envs%f_env_id,f_env=f_env,&
               new_error=new_error, failure=failure)
          CALL force_env_get(f_env%force_env,qs_env=sub_qs_env,error=new_error)
          CALL get_qs_env(sub_qs_env,mos=sub_mos,error=new_error)
          DO ispin=1,ep_env%nspins
             CALL get_mo_set(sub_mos(ispin)%mo_set,nao=ep_env%sub_nao(ispin),&
                  nmo=ep_env%sub_nmo(ispin))
          END DO

          CALL f_env_rm_defaults(f_env,new_error,ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL mp_max(ep_env%sub_nmo,para_env%group)
       CALL mp_max(ep_env%sub_nao,para_env%group)
       ! dummy force calc to ensure that everything is allocated for force
       ! evaluation
       CALL calc_energy_force(ep_env%mol_envs%f_env_id,.TRUE.,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       ! guarantee wf & overlap in wf_history
       DO imol=1,SIZE(ep_env%mol_envs%local_rep_indices)
          ep_env%mol_envs%wf_history(imol)%wf_history%store_wf=.TRUE.
          ep_env%mol_envs%wf_history(imol)%wf_history%store_overlap=.TRUE.
          CALL wfi_change_memory_depth(&
               ep_env%mol_envs%wf_history(imol)%wf_history,&
               MAX(1,ep_env%mol_envs%wf_history(imol)%wf_history%memory_depth),&
               error=error)
       END DO
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " mol_envs setup fininshed",error=error)
       ! create the perturbation environment
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " creating sub_p_env",error=error)
       CALL p_env_create(ep_env%sub_p_env, qs_env=sub_qs_env,&
            orthogonal_orbitals=.TRUE.,error=error)
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " sub_p_env created",error=error)

       ! create sub_proj and at2sub
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " creating projections",error=error)
       ALLOCATE(sub_proj(nmol),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       ALLOCATE(ep_env%at2sub(nat),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       ALLOCATE(atoms(at_per_mol),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       iat=0
       DO imol=1,nmol
          DO i=1,at_per_mol
             iat=iat+1
             atoms(i)=iat
             ep_env%at2sub(iat)=imol
          END DO
          CALL p_proj_create(sub_proj(imol)%projection, atoms=atoms,&
               particle_set=particle_set,qs_kind_set=qs_kind_set, error=error)
       END DO
       DEALLOCATE(atoms,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ep_env%sub_proj => sub_proj
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " created projections",error=error)

       ! create the perturbation environment
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " creating main_p_env",error=error)
       CALL p_env_create(ep_env%main_p_env, qs_env=ep_env%main_qs_env,&
            orthogonal_orbitals=.FALSE.,error=error)
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " created main_p_env",error=error)

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

       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " ep initial setup finished",error=error)
    END IF
    CALL cp_error_check(error,failure)
    IF (failure) THEN
       ierr=1
    ELSE
       ierr=0
    END IF
    CALL cp_error_dealloc_ref(error)
    CALL timestop(handle)
  END SUBROUTINE ep_env_init_low

! *****************************************************************************
!> \brief finish the initialization of the qs environement
!> \param qs_env the qs_environment to fully initialize
!> \param globenv global environment finish the the initialization
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE ep_env_finish_qs_init(qs_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_finish_qs_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    failure=.FALSE.
    NULLIFY(ks_env, mos, mo_coeff, dft_control,matrix_ks,matrix_s)

    CALL timeset(routineN,handle)
    para_env=>qs_env%para_env

    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=qs_env%input,error=error)

    ! *** Calculate the overlap and the core Hamiltonian integral matrix ***

    CALL build_core_hamiltonian_matrix(qs_env=qs_env,calculate_forces=.FALSE.,error=error)
    CALL qs_env_update_s_mstruct(qs_env, error=error)

    !   *** Initializes the MOs ***
    CALL get_qs_env(qs_env,mos=mos,error=error)
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO ispin=1,SIZE(mos)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
          IF (.NOT.ASSOCIATED(mo_coeff)) THEN
             CALL init_mo_set(mos(ispin)%mo_set, &
                  qs_env%mpools%ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
       END DO
    END IF

    !   *** Allocate k and put it in the QS environment ***
    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    dft_control=dft_control,&
                    matrix_s=matrix_s,&
                    error=error)

    CALL cp_dbcsr_allocate_matrix_set(matrix_ks,dft_control%nspins,error=error)
    DO ispin=1,dft_control%nspins
       ALLOCATE (matrix_ks(ispin)%matrix)
       CALL cp_dbcsr_init(matrix_ks(ispin)%matrix,error=error)
       CALL cp_dbcsr_copy(matrix_ks(ispin)%matrix,matrix_s(1)%matrix,&
            "KOHN-SHAM MATRIX-"//ADJUSTL(cp_to_string(ispin)),error=error)
    END DO

    CALL set_ks_env(ks_env,matrix_ks=matrix_ks,error=error)

    CALL calculate_ecore_self(qs_env,error=error)
    CALL calculate_ecore_overlap(qs_env,para_env,.FALSE.,error=error)

    CALL timestop(handle)
  END SUBROUTINE ep_env_finish_qs_init

! *****************************************************************************
!> \brief transfers psi0 from the sub_qs_env to the main_qs_env
!>      (i.e. sub_qs_env%c are injected into main_qs_env%c)
!> \param ep_env ep environment in which the transfer should be performed
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE ep_env_transfer_psi0(ep_env, error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_transfer_psi0', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: max_blocksize = 100

    INTEGER :: blocksize, blocksize2, full_nmo, handle, i, icol, ierr, ispin, &
      isub, local_idx, my_start_col_full, my_start_col_min, ncol_min, nmo, &
      nrow_min, nspins, stat, sub_nmo, unit_nr
    LOGICAL                                  :: failure, print_matrixes
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sub_mo_matrix, tmp_full
    TYPE(cp_error_type)                      :: new_error, suberror
    TYPE(cp_fm_type), POINTER                :: orbitals, sub_orbitals
    TYPE(cp_logger_type), POINTER            :: logger, new_logger
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_environment_type), POINTER       :: main_qs_env
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), POINTER                  :: sub_proj
    TYPE(qs_p_projection_type), POINTER      :: proj
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot
    TYPE(replica_env_type), POINTER          :: rep_env

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(sub_orbitals, main_qs_env, sub_proj, logger, orbitals)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    CALL ep_env_get(ep_env,&
         sub_proj=sub_proj,main_qs_env=main_qs_env,error=error)
    CALL cp_error_check(error,failure)
    CPPrecondition(ASSOCIATED(sub_proj),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(main_qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       nspins=main_qs_env%dft_control%nspins
       CALL get_qs_env(main_qs_env,mos=mos,error=error)
       rep_env => ep_env%mol_envs

       DO ispin=1,nspins ! transfer also occupation_numbers ?
          my_start_col_full=1
          sub_nmo=ep_env%sub_nmo(ispin)
          full_nmo=ep_env%full_nmo(ispin)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=orbitals,nmo=nmo)
          CALL cp_fm_set_all(orbitals,0.0_dp,error=error)
          blocksize=MIN(max_blocksize,sub_nmo)
          ALLOCATE(sub_mo_matrix(blocksize,ep_env%sub_nao(ispin)),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          ALLOCATE(tmp_full(blocksize,ep_env%full_nao(ispin)),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          print_matrixes=BTEST(cp_print_key_should_output(logger%iter_info,&
               ep_env%input,"PRINT%EP_MATRIXES/PSI0_BLOCKS",error=error),cp_p_file)
          DO isub=1,ep_env%nmol
             ! does not work if there are "empty" processors in replica env
             CPAssert(ASSOCIATED(rep_env),cp_failure_level,routineP,error,failure)

             my_start_col_min=1
             proj => ep_env%sub_proj(isub)%projection
             blocksize2=blocksize
             DO icol=0,sub_nmo-1,blocksize
                IF (icol+blocksize>sub_nmo) blocksize2=sub_nmo-icol
                IF (rep_env%rep_is_local(isub)) THEN
                   CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
                        new_error=new_error, failure=failure)
                   local_idx=rep_env_local_index(rep_env,isub,error=error)
                   CPAssert(local_idx>0,cp_failure_level,routineP,error,failure)
                   new_logger => cp_error_get_logger(new_error)
                   snapshot => wfi_get_snapshot(&
                        rep_env%wf_history(local_idx)%wf_history,&
                        1,error=new_error)
                   CPAssert(ASSOCIATED(snapshot),cp_failure_level,routineP,new_error,failure)
                   CPAssert(ASSOCIATED(snapshot%wf),cp_failure_level,routineP,new_error,failure)
                   sub_orbitals => snapshot%wf(ispin)%matrix
                   new_logger%iter_info%iteration(1)=isub
                   IF (print_matrixes) THEN
                      CALL cp_error_init(suberror,template_error=new_error)
                      unit_nr=cp_print_key_unit_nr(new_logger,ep_env%input,&
                           "PRINT%EP_MATRIXES/PSI0_BLOCKS",extension=".psi0_block",&
                           ignore_should_output=.TRUE.,error=new_error)
                      CALL cp_fm_write(sub_orbitals, unit_nr,&
                           long_description=.TRUE., error=suberror)
                      CALL cp_error_reset(suberror)
                      CALL cp_print_key_finished_output(unit_nr,new_logger,ep_env%input,&
                           "PRINT%EP_MATRIXES/PSI0_BLOCKS",&
                           ignore_should_output=.TRUE.,error=new_error)
                      CALL cp_error_dealloc_ref(suberror)
                   END IF
                   CALL cp_fm_get_info(sub_orbitals, nrow_global=nrow_min, ncol_global=ncol_min,&
                        error=error)
                   CPPostcondition(nrow_min==ep_env%sub_nao(ispin),cp_failure_level,routineP,error,failure)
                   CPPostcondition(ncol_min==sub_nmo,cp_failure_level,routineP,error,failure)
                   CALL cp_assert(blocksize2+icol<=sub_nmo+1,&
                        cp_failure_level,cp_assertion_failed,routineP,&
                        "icol+blocksize2>sub_nmo"//&
                        CPSourceFileRef,&
                        error,failure)
                   CALL cp_fm_get_submatrix(sub_orbitals,sub_mo_matrix,transpose=.TRUE.,&
                        start_col=icol+1, n_cols=blocksize2,error=new_error)
                   CALL f_env_rm_defaults(f_env,new_error,ierr)
                   CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
                END IF
                CALL mp_bcast(sub_mo_matrix,rep_env%inter_rep_rank(rep_env%replica_owner(isub)),&
                     rep_env%para_env_inter_rep%group)

                ! be smarter and use only a part of the full width?
                ! i.e.size(tmp_full,2) could easily reduced to (min(proj_indexes):max(proj_indexes))
                CALL dcopy(SIZE(tmp_full,1)*SIZE(tmp_full,2),0.0_dp,0,tmp_full(1,1),1)
                DO i=1,SIZE(proj%proj_indexes)
                   tmp_full(1:blocksize2,proj%proj_indexes(i))=sub_mo_matrix(1:blocksize2,i)
                END DO

                CALL cp_assert(icol+my_start_col_full+blocksize2<=ep_env%full_nmo(ispin)+1,&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     'icol+my_start_col_full+blocksize2>ep_env%full_nmo(ispin)+1 '//&
                     CPSourceFileRef,&
                     error,failure)
                CALL cp_fm_set_submatrix(orbitals,tmp_full,&
                     start_col=icol+my_start_col_full, n_cols=blocksize2, &
                     transpose=.TRUE.,error=error)
             END DO
             my_start_col_full=my_start_col_full+sub_nmo
             CALL cp_assert(my_start_col_full<=full_nmo+1,&
                  cp_failure_level,cp_assertion_failed,routineP,&
                  "my_start_col_full>full_nmo+1, "//&
                  CPSourceFileRef,&
                  error,failure)

          END DO
          CALL cp_assert(my_start_col_full==ep_env%full_nmo(ispin)+1,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "my_start_col_full/=ep_env%full_nmo(ispin)+1, "//&
               CPSourceFileRef,&
               error,failure)

          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               ep_env%input,"PRINT%EP_MATRIXES/PSI0",error=error),cp_p_file)) THEN
             CALL cp_error_init(suberror,template_error=error)
             unit_nr=cp_print_key_unit_nr(logger,ep_env%input,&
                  "PRINT%EP_MATRIXES/PSI0",extension=".psi0",&
                  error=error)
             CALL cp_fm_write(orbitals, unit_nr,&
                  long_description=.TRUE., error=suberror)
             CALL cp_error_reset(suberror)
             CALL cp_print_key_finished_output(unit_nr,logger,ep_env%input,&
                  "PRINT%EP_MATRIXES/PSI0",&
                  error=error)
             CALL cp_error_dealloc_ref(suberror)
          END IF
       END DO

    END IF
    CALL timestop(handle)
  END SUBROUTINE ep_env_transfer_psi0

! *****************************************************************************
!> \brief solves a linear system with conjugated gradient
!> \param ep_env the sytem to solve
!> \param eps_r wanted residual error
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author fawzi
!> \note
!>      The linear solver should go in another module, this is a quick and
!>      dirty implementation.
!>      A: postortho (p_op_l1+p_op_l2) preortho
!>      x: ep_env%psi1
!>      ugly scaling of p (to avoid cutoff of small values on the grids)
! *****************************************************************************
  SUBROUTINE stupid_solve(ep_env,eps_r,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    REAL(KIND=dp), INTENT(in)                :: eps_r
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'stupid_solve', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ispin, j, nspins, &
                                                output_unit
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: alpha, Ap_p, Ap_p_s, beta, p_ortho_norm, p_ortho_norm_s, &
      r_norm2, r_norm2_s, r_z, r_z_new, r_z_s
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: Ap, p, p_ortho, r, z
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_matrix_pools_type), POINTER      :: mpools

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
       CALL cp_add_iter_level(logger%iter_info,level_name="EP_LIN_SOLVER",error=error)
       CALL cp_iterate(logger%iter_info,error=error)
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
                  matrix_out=z(ispin)%matrix,error=error)
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
             CALL cp_iterate(logger%iter_info,error=error)
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
                CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,&
                     p_ortho(ispin)%matrix,Ap(ispin)%matrix,&
                     ncol=ep_env%main_p_env%n_mo(ispin),&
                     error=error)
                CALL cp_fm_trace(matrix_a=p_ortho(ispin)%matrix,&
                     matrix_b=Ap(ispin)%matrix,&
                     trace=p_ortho_norm_s,error=error)
                p_ortho_norm=p_ortho_norm+p_ortho_norm_s
             END DO
             p_ortho_norm=SQRT(p_ortho_norm)

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

             Ap_p=0.0_dp
             DO ispin=1,nspins
                CALL cp_fm_trace(matrix_a=Ap(ispin)%matrix,&
                     matrix_b=p(ispin)%matrix,trace=Ap_p_s,error=error)
                Ap_p=Ap_p+Ap_p_s
             END DO

             ! this alpha is the alpha in alg. desc. times p_ortho_norm
             alpha=r_z/Ap_p

             ! x=x+alpha p_ortho
             DO ispin=1,nspins
                CALL cp_fm_scale_and_add(alpha=1.0_dp,matrix_a=ep_env%psi1(ispin)%matrix,&
                     beta=alpha,matrix_b=p_ortho(ispin)%matrix,&
                     error=error)
             END DO

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
                     matrix_out=z(ispin)%matrix,error=error)
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

          output_unit=cp_print_key_unit_nr(logger,ep_env%input,"PRINT%RUN_INFO/LIN_SOLV",&
               extension=".epLog",error=error)
          IF (output_unit>0) THEN
             WRITE(output_unit,"(a,a,i6,a,f20.15)")routineP,"i=",i,&
                  "j=final r_norm2=",r_norm2
             CALL cp_print_key_finished_output(output_unit,logger,ep_env%input,&
                  "PRINT%RUN_INFO/LIN_SOLV",&
                  error=error)
          END IF
       END DO

       CALL cp_rm_iter_level(logger%iter_info,level_name="EP_LIN_SOLVER",error=error)
       output_unit=cp_print_key_unit_nr(logger,ep_env%input,"PRINT%RUN_INFO/LIN_SOLV",&
            extension=".epLog",error=error)
       IF (output_unit>0) THEN
          WRITE(output_unit,"(a,f20.15)") "solve finished, r_norm2=",r_norm2
          CALL cp_print_key_finished_output(output_unit,logger,ep_env%input,&
               "PRINT%RUN_INFO/LIN_SOLV",&
               error=error)
       END IF
       CPPostcondition(r_norm2<eps_r**2,cp_warning_level,routineP,error,failure)

       CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,p,error=error)
       CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,Ap,error=error)
       CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,r,error=error)
       CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,z,error=error)
       CALL fm_pools_give_back_fm_vect(ao_mo_fm_pools,p_ortho,error=error)
    END IF

  END SUBROUTINE stupid_solve

! *****************************************************************************
!> \brief evaluates the ep operator (p_op_l1+p_op_l2) on the given vector
!> \param p_env perturbation calculation environment
!> \param qs_env the qs_env that is perturbed by this p_env
!> \param v the matrix to operate on (must be orthogonal to the psi0)
!> \param res the result
!> \param error error handling object (optional)
!> \par History
!>      12.2002 created [fawzi]
!> \author fawzi
!> \note
!>      no preorthogonalization, thus v must be orthogonal to the psi0
! *****************************************************************************
  SUBROUTINE p_op_ep(p_env, qs_env, v, res, error)
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(in)                             :: v
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(inout)                          :: res
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_op_ep', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, lfomo, nmo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: maxocc, scale_factor
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)

  END SUBROUTINE p_op_ep

! *****************************************************************************
!> \brief writes out the whole ep interaction matrix (as defined by p_op_ep)
!>      one column per file
!> \param p_env the p_env where to evaluate the interaction matrix
!> \param qs_env the qs_env that is perturbed by the p_env
!> \param name name used for the matrixes outputted
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE p_env_write_ep_matrix(p_env, qs_env, name, error)
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    CHARACTER(len=*), INTENT(in)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'p_env_write_ep_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iao, imo, ispin, ispin2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: v_norm, v_norm_s
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: res, v
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control

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
                   CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,&
                        v(ispin2)%matrix,res(ispin2)%matrix,&
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

! *****************************************************************************
!> \brief localizes the wavefunction matrix, restraining it the the subsystem
!> \param ep_env ...
!> \param wf_coeff ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_env_localize_matrix(ep_env, wf_coeff, error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: wf_coeff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_localize_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, isub, nmo, &
                                                nspins, start_col
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: main_qs_env
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), POINTER                  :: sub_proj

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(main_qs_env, sub_proj)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    CALL ep_env_get(ep_env,&
         sub_proj=sub_proj,main_qs_env=main_qs_env,error=error)
    CALL cp_error_check(error,failure)
    CPPrecondition(ASSOCIATED(sub_proj),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(main_qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       nspins=main_qs_env%dft_control%nspins
       DO ispin=1,nspins
          start_col=1
          CALL cp_fm_get_info(wf_coeff(ispin)%matrix,ncol_global=nmo,error=error)

          DO isub=1,SIZE(sub_proj)
             CPInvariant(start_col+ep_env%sub_nmo(ispin)<nmo+2,cp_failure_level,routineP,error,failure)
             start_col=start_col+ep_env%sub_nmo(ispin)
          END DO
          CPPostcondition(start_col==nmo+1,cp_failure_level,routineP,error,failure)
       END DO
    END IF

    CALL timestop(handle)
  END SUBROUTINE ep_env_localize_matrix

! *****************************************************************************
!> \brief calculates the energy and forces using the ep method
!> \param ep_env_id the id of the ep_env that should be updated
!> \param calc_f if non zero then calculates also the forces
!> \param ierr will be nonzero if there is an error
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_env_calc_e_f_low(ep_env_id, calc_f, ierr)
    INTEGER, INTENT(in)                      :: ep_env_id, calc_f
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_calc_e_f_low', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(ep_env_type), POINTER               :: ep_env

    failure=.FALSE.

    CALL cp_error_init(error)
    ep_env => ep_envs_get_ep_env(ep_env_id)
    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL ep_env_calc_e0(ep_env,error=error)
    END IF
    CALL cp_error_check(error,failure)
    IF (failure) THEN
       ierr=1
    ELSE
       ierr=0
    END IF
    CALL cp_error_dealloc_ref(error)
  END SUBROUTINE ep_env_calc_e_f_low

! *****************************************************************************
!> \brief checks the allocation status of the needed matrixes, and if necessary
!>      allocate them
!> \param qs_env the qs_env, the scf_env lives in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2006 taken from scf_env_check_i_alloc [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE qs_check_i_alloc(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_check_i_alloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, ispin
    LOGICAL                                  :: failure, &
                                                gth_potential_present, &
                                                uniform_occupation
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env

!                                                my_transition_potential, &

    CALL timeset(routineN,handle)

    NULLIFY(matrix_ks, ao_mo_fm_pools, matrix_s, dft_control, mos, &
         ks_env)
    NULLIFY(qs_kind_set, mo_coeff)
    logger => cp_error_get_logger(error)
    failure=.FALSE.

!    my_transition_potential = .FALSE.
    uniform_occupation = .TRUE.

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         mos=mos,&
         matrix_ks=matrix_ks,&
         ks_env=ks_env,&
         qs_kind_set=qs_kind_set,&
         matrix_s=matrix_s,error=error)
    CALL get_qs_kind_set(qs_kind_set, gth_potential_present=gth_potential_present)
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
         error=error)


    !   *** finish initialization of the MOs ***
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO ispin=1,SIZE(mos)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
          IF (.NOT.ASSOCIATED(mo_coeff)) THEN
             CALL init_mo_set(mos(ispin)%mo_set,&
                  ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
!          IF(my_transition_potential .AND. ispin==1) THEN
!             CALL get_mo_set(mos(ispin)%mo_set,&
!                  occupation_numbers=occupation_numbers)
!             occupation_numbers(dft_control%xas_estate) = &
!                  dft_control%xas_control%occ_estate
!             occupation_numbers(homo) = dft_control%xas_control%occ_homo
!             uniform_occupation = .FALSE.
!
!             CALL set_mo_set(mos(ispin)%mo_set,&
!                  uniform_occupation=uniform_occupation,error=error)
!          END IF
       END DO
    END IF

    !   *** get the mo_derivs OK if needed ***
    IF (qs_env%requires_mo_derivs) THEN
       ! this is never used ...
       !CALL get_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
       !IF (.NOT.ASSOCIATED(mo_derivs)) THEN
       !   IF (dft_control%restricted) THEN ! right now, there might be more mos than needed derivs
       !      ALLOCATE(mo_derivs(1))
       !      CALL get_mo_set(mos(1)%mo_set,mo_coeff=mo_coeff)
       !      CALL cp_fm_create(mo_derivs(1)%matrix,mo_coeff%matrix_struct,error=error)
       !   ELSE
       !      ALLOCATE(mo_derivs(dft_control%nspins))
       !      DO ispin=1,dft_control%nspins
       !         CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
       !         CALL cp_fm_create(mo_derivs(ispin)%matrix,mo_coeff%matrix_struct,error=error)
       !      ENDDO
       !   ENDIF
       !   CALL set_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
       !ENDIF
    ELSE
       ! nothing should be done
    ENDIF

    !   *** Allocate matrix_ks and put it in the QS environment ***

    IF (.not.ASSOCIATED(matrix_ks)) THEN
       CALL cp_dbcsr_allocate_matrix_set( matrix_ks, dft_control%nspins, error )
       DO ispin=1,SIZE(matrix_ks)
          ALLOCATE(matrix_ks(ispin)%matrix)
          CALL cp_dbcsr_init(matrix_ks(ispin)%matrix,error=error)
          CALL cp_dbcsr_copy(matrix_ks(ispin)%matrix, matrix_s(ispin)%matrix,&
               name="EP"//TRIM(cp_iter_string(logger%iter_info,error=error))//&
               "KOHN-SHAM_MATRIX", error=error)
       ENDDO

       CALL set_ks_env(ks_env,matrix_ks=matrix_ks,error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_check_i_alloc

! *****************************************************************************
!> \brief ...
!> \param ep_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_env_calc_e0(ep_env,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_calc_e0', &
      routineP = moduleN//':'//routineN

    INTEGER :: flag, handle, i, iat, iatom, idir, ii, ikind, imol, ispin, iw, &
      lfomo, natom, nkind, nmo, nspin, output_unit, stat, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                natom_of_kind
    INTEGER, DIMENSION(:), POINTER           :: at
    LOGICAL                                  :: do_rotations, e0_only, failure
    REAL(dp)                                 :: e_corr_s, eps_lin_solv, maxocc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_2d_r_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: rotatedC0
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: kinetic, matrix_ks, matrix_s, &
                                                matrix_w
    TYPE(cp_error_type)                      :: suberror
    TYPE(cp_fm_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: psi0
    TYPE(cp_fm_type), POINTER                :: orbitals
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: print_section, qs_env_input

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY (logger)
    logger => cp_error_get_logger(error)

    NULLIFY (atomic_kind_set, dft_control, force, ks_env, mos, matrix_ks, &
         matrix_s, matrix_w, particle_set, rho, scf_control,para_env,&
         blacs_env,orbitals, kinetic)

    globenv => ep_env%globenv
    qs_env => ep_env%main_qs_env
    para_env=>qs_env%para_env

    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " started",error=error)
    CALL cp_error_init(suberror,template_error=error)

    CALL qs_env_rebuild_pw_env(qs_env,error)

    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         dft_control=dft_control,&
         force=force,&
         particle_set=particle_set,&
         scf_control=scf_control,&
         ks_env=ks_env,&
         input=qs_env_input,&
         subsys=subsys,&
         energy=energy,&
         mos=mos,&
         rho=rho,&
         error=error)


    CALL qs_subsys_get(subsys,particles=particles,error=error)


    ! ensure storage of electronic kinetic energy
    CALL section_vals_val_get(qs_env_input,"DFT%PRINT%AO_MATRICES%__CONTROL_VAL",&
         i_val=flag,error=error)
    flag=IBSET(flag,cp_p_store)
    CALL section_vals_val_set(qs_env_input,"DFT%PRINT%AO_MATRICES%__CONTROL_VAL",&
         i_val=flag,error=error)
    CALL section_vals_val_set(qs_env_input,"DFT%PRINT%AO_MATRICES%KINETIC_ENERGY",&
         l_val=.TRUE.,error=error)

    natom = SIZE(particle_set)

    ! zero out the forces
    DO iatom=1,natom
       particle_set(iatom)%f=0.0_dp
    END DO
    ! zero energy
    CALL ep_energy_zero(ep_env%energy,error=error)

    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

    ALLOCATE (kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         atom_of_kind=atom_of_kind,&
         kind_of=kind_of)

    IF (.NOT.ASSOCIATED(force)) THEN
       !   *** Allocate the force data structure ***
       nkind = SIZE(atomic_kind_set)
       ALLOCATE (natom_of_kind(nkind),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            natom_of_kind=natom_of_kind)
       CALL allocate_qs_force(force,natom_of_kind,error)
       DEALLOCATE (natom_of_kind,STAT=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CALL qs_subsys_set(subsys,force=force,error=error)
    END IF
    CALL zero_qs_force(force,error)
    CALL ep_force_zero(ep_env%force,error=error)

    ! !!!!! qs energy
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " qs_energy setup in main_qs_env",error=error)
    CALL build_qs_neighbor_lists(qs_env,para_env,force_env_section=qs_env%input,error=error)

    ! *** Calculate the overlap and the core Hamiltonian integral matrix ***
    CALL build_core_hamiltonian_matrix(qs_env=qs_env,calculate_forces=.FALSE.,error=error)
    CALL qs_env_update_s_mstruct(qs_env,error=error)
    CALL qs_check_i_alloc(qs_env,error=error)
    CALL calculate_ecore_self(qs_env,error=error)
    CALL calculate_ecore_overlap(qs_env, para_env, &
         calculate_forces=.FALSE.,error=error)

    ! *** calculate psi0
    ! ** update positions of replicas
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " calculating psi0",error=error)
    rep_env => ep_env%mol_envs
    CALL get_qs_env(ep_env%main_qs_env,particle_set=particle_set,&
         para_env=para_env,blacs_env=blacs_env,mos=mos,error=error)
    DO imol=1,ep_env%nmol
       at => ep_env%sub_proj(imol)%projection%atoms
       ii=0
       DO iat=1,SIZE(at)
          DO idir=1,3
             ii=ii+1
             rep_env%r(ii,imol)=particle_set(at(iat))%r(idir)
          END DO
       END DO
    END DO

    CALL section_vals_val_get(ep_env%root_section,"FORCE_EVAL%EP%ROTATE",&
         l_val=do_rotations,error=error)
    IF (.NOT.do_rotations) THEN
       ! ** calc e and psi00
       IF (ASSOCIATED(rep_env)) THEN
          CALL rep_env_calc_e_f(rep_env,calc_f=.FALSE.,error=error)
       END IF
       ep_env%energy%e_no_int=SUM(rep_env%f(rep_env%ndim+1,:))
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " calculated psi0 locally, e_no_int="//&
            TRIM(ADJUSTL(cp_to_string(ep_env%energy%e_no_int))),error=error)
       ! ** transfer psi0
       CALL ep_env_transfer_psi0(ep_env,error=error)
    ELSE
       ! matrix for rotated C
       ALLOCATE(rotatedC0(ep_env%nspins),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       DO ispin=1,ep_env%nspins
          ALLOCATE(rotatedC0(ispin)%array(ep_env%sub_nmo(ispin),ep_env%sub_nao(ispin)),stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=orbitals,nmo=nmo)
          CALL cp_fm_set_all(orbitals,0.0_dp,error=error)
       END DO
       DO imol=1,ep_env%nmol
          IF (rep_env%rep_is_local(imol)) THEN
             ! this should return the C0 for a molecule
             CALL calc_c0_tilde(ep_env,&
                  mol_pos=rep_env%r(:,imol),new_c=rotatedC0,&
                  imol=imol,error=error)
          END IF
          DO ispin=1,ep_env%nspins
             ! transfer to matrix in the full systems (psi0)
             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=orbitals,nmo=nmo)
             CALL mp_bcast(rotatedC0(ispin)%array, rep_env%replica_owner(imol),&
                  rep_env%para_env_inter_rep%group)
             CALL cp_fm_set_submatrix(orbitals,rotatedC0(ispin)%array,&
                  start_col=(imol-1)*ep_env%sub_nmo(ispin)+1,&
                  start_row=(imol-1)*ep_env%sub_nao(ispin)+1,&
                  transpose=.TRUE.,error=error)
          END DO
       END DO
       DO ispin=1,ep_env%nspins
          DEALLOCATE(rotatedC0(ispin)%array,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END DO
       DEALLOCATE(rotatedC0,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " transferred psi0 to main system",error=error)

    ! update main_p_env & psi0d
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " update main_p_env",error=error)
    CALL p_env_did_change(ep_env%sub_p_env, s_struct_changed=.TRUE., &
         grid_changed=.FALSE.,error=error)
    CALL p_env_psi0_changed(ep_env%main_p_env, qs_env=ep_env%main_qs_env,&
         Hrho_psi0d=ep_env%m_pi_Hrho_psi0d, error=error)
    CALL p_postortho(p_env=ep_env%main_p_env,qs_env=ep_env%main_qs_env,&
         v=ep_env%m_pi_Hrho_psi0d, &
         n_cols=ep_env%main_p_env%n_mo,&
         error=error)

    ! ** update qs_env%rho
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " update qs_env%rho and ks matrix",error=error)
    CALL get_qs_env(qs_env, rho=rho, mos=mos, error=error)
    DO ispin=1,ep_env%nspins
       CALL cp_dbcsr_set(rho%rho_ao(ispin)%matrix,0._dp,error=error)
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho%rho_ao(ispin)%matrix,&
            matrix_v=mos(ispin)%mo_set%mo_coeff,&
            matrix_g=ep_env%main_p_env%psi0d(ispin)%matrix,&
            ncol=ep_env%full_nmo(ispin),&
            error=error)
       CALL get_mo_set(mos(ispin)%mo_set, maxocc=maxocc)
       CALL cp_dbcsr_scale(rho%rho_ao(ispin)%matrix,alpha_scalar=maxocc,error=error)
    END DO
    CALL qs_rho_update_rho(rho, qs_env=qs_env, error=error)
    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)
    CALL qs_ks_update_qs_env(ks_env,qs_env,calculate_forces=.FALSE.,&
         just_energy=.FALSE.,error=error)
    ep_env%energy%e0=energy%total
    unit_nr=cp_print_key_unit_nr(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",error=error)
    CALL print_qs_energies(qs_env,unit_nr,error)
    IF (unit_nr>0) CALL cp_print_key_finished_output(unit_nr,&
         logger,ep_env%input,"PRINT%RUN_INFO", error=error)

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
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " finished qs_enegy setup",error=error)
    ! !!!! end qs_energy

    ! init preconditioner
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " init proeconditioner main_p_env",error=error)
    ALLOCATE(psi0(dft_control%nspins),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO ispin=1,dft_control%nspins
          NULLIFY(psi0(ispin)%matrix)
          CALL get_mo_set(ep_env%main_qs_env%mos(ispin)%mo_set,&
               mo_coeff=psi0(ispin)%matrix)
       END DO
       ALLOCATE(ep_env%precond(dft_control%nspins),stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ispin=1
       NULLIFY(ep_env%precond(ispin)%preconditioner)
       ALLOCATE(ep_env%precond(ispin)%preconditioner,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CALL init_preconditioner(ep_env%precond(ispin)%preconditioner,&
            para_env=para_env,&
            blacs_env=blacs_env,error=error)
       CALL get_qs_env(qs_env,&
                       matrix_ks=matrix_ks,&
                       matrix_s=matrix_s,&
                       kinetic=kinetic,&
                       error=error)
       CALL make_preconditioner(ep_env%precond(ispin)%preconditioner,&
            ot_precond_full_kinetic,&
            ot_precond_solver_default,&
            matrix_ks(ispin)%matrix,&
            matrix_s(1)%matrix,&
            kinetic(1)%matrix,&
            mos(ispin)%mo_set,0.2_dp,error=error)
       ! same precond for both spins
       DO ispin=2,dft_control%nspins
          ep_env%precond(ispin)%preconditioner => &
               ep_env%precond(1)%preconditioner
       END DO
    END IF
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " preconditioner setup main_p_env finished",error=error)

    CALL section_vals_val_get(ep_env%root_section,"FORCE_EVAL%EP%E0_ONLY",&
         l_val=e0_only,error=error)
    IF (.NOT. e0_only) THEN
       ! *** perturb to do
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " perturbation start",error=error)
       CALL section_vals_val_get(ep_env%root_section,"FORCE_EVAL%EP%EPS_LIN_SOLV",&
            r_val=eps_lin_solv,error=error)
       CALL stupid_solve(ep_env,eps_r=eps_lin_solv,error=error)

       ep_env%energy%e1=0.0_dp
       DO ispin=1,ep_env%nspins
          CALL cp_fm_trace(matrix_a=ep_env%m_pi_Hrho_psi0d(ispin)%matrix,&
               matrix_b=ep_env%psi1(ispin)%matrix,trace=e_corr_s,error=error)
          ep_env%energy%e1=ep_env%energy%e1-e_corr_s
       END DO
       CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
            extension=".epLog",message=routineP//&
            " perturbation end",error=error)
       ! *** end perturb
    END IF

    ep_env%energy%e_tot=ep_env%energy%e0+ep_env%energy%e1
    output_unit=cp_print_key_unit_nr(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",error=error)
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(T3,A,T56,F25.14)")&
            "no interaction energy:",ep_env%energy%e_no_int
       WRITE (UNIT=output_unit,FMT="(T3,A,T56,F25.14)")&
            "psi0 energy (e0):",ep_env%energy%e0
       WRITE (UNIT=output_unit,FMT="(T3,A,T56,F25.14)")&
            "correction energy (e1):",ep_env%energy%e1
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T56,F25.14))")&
            "Total energy:                                  ",ep_env%energy%e_tot
       CALL cp_print_key_finished_output(unit_nr,logger,ep_env%input,&
            "PRINT%RUN_INFO",&
            error=error)
    END IF

    output_unit=cp_print_key_unit_nr(logger,ep_env%input,"PRINT%ENERGY",&
         extension=".epEnergy",error=error)
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(4F25.14)")&
            ep_env%energy%e_no_int,ep_env%energy%e0,ep_env%energy%e1,ep_env%energy%e_tot
       CALL cp_print_key_finished_output(output_unit,logger,ep_env%input,&
            "PRINT%ENERGY", error=error)
    END IF

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         ep_env%input,"PRINT%EP_MATRIXES/PSI1",error=error),cp_p_file)) THEN
       unit_nr=cp_print_key_unit_nr(logger,ep_env%input,&
            "PRINT%EP_MATRIXES/PSI1",extension=".psi1",&
            error=error)
       CALL cp_fm_vect_write(ep_env%psi1, unit_nr,&
            long_description=.TRUE., error=suberror)
       CALL cp_error_reset(suberror)
       CALL cp_print_key_finished_output(unit_nr,logger,ep_env%input,&
            "PRINT%EP_MATRIXES/PSI1",&
            error=error)
    END IF
    IF (ASSOCIATED(ep_env%main_p_env%rho1)) THEN
       IF (ep_env%main_p_env%rho1%rho_r_valid) THEN
          IF (BTEST(cp_print_key_should_output(logger%iter_info,ep_env%input,&
               "PRINT%EP_RHO_CUBE/RHO1",error=error),cp_p_file) ) THEN
             DO ispin=1,dft_control%nspins
                unit_nr=cp_print_key_unit_nr(logger,ep_env%input,"PRINT%EP_RHO_CUBE/RHO1",&
                     extension=".cube",middle_name="rho1_"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                     log_filename=.FALSE.,error=error)
                CALL cp_pw_to_cube ( ep_env%main_p_env%rho1%rho_r(1)%pw, unit_nr=unit_nr,&
                     title="rho1_"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                     particles=particles,&
                     stride=section_get_ivals(ep_env%input,"PRINT%EP_RHO_CUBE%STRIDE",error=error),&
                     error=suberror)
                CALL cp_error_reset(suberror)
                CALL cp_print_key_finished_output(unit_nr,logger,ep_env%input,&
                     "PRINT%EP_RHO_CUBE/RHO1",error=error)
             END DO
          END IF
       END IF
    END IF

    IF (BTEST(cp_print_key_should_output(logger%iter_info,ep_env%input,&
         "PRINT%EP_RHO_CUBE/RHO0",error=error),cp_p_file) ) THEN
       DO ispin=1,dft_control%nspins
          unit_nr=cp_print_key_unit_nr(logger,ep_env%input,"PRINT%EP_RHO_CUBE/RHO0",&
               extension=".cube",middle_name="rho0_"//TRIM(ADJUSTL(cp_to_string(ispin))),&
               log_filename=.FALSE.,error=error)
          CALL cp_pw_to_cube (rho%rho_r(1)%pw, unit_nr=unit_nr,&
               title="rho0_"//TRIM(ADJUSTL(cp_to_string(ispin))),&
               particles=particles,&
               stride=section_get_ivals(ep_env%input,"PRINT%EP_RHO_CUBE%STRIDE",error=error),&
               error=suberror)
          CALL cp_error_reset(suberror)
          CALL cp_print_key_finished_output(unit_nr,logger,ep_env%input,&
               "PRINT%EP_RHO_CUBE/RHO0",error=error)
       END DO
    END IF

    ! forces
    ! Build W matrix (from the orthonormality constraint lambda C^T S C)
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " starting forces",error=error)

    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " calculate dc/dr forces",error=error)
    ! commented out (fawzi    24-Apr-07)
    !    CALL ep_calc_f_sub(ep_env,error=error)
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " finished calculating forces dc/dr",error=error)

    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " creating w_matrix main_qs_env",error=error)
    nspin = SIZE(mos)
    CALL get_qs_env(qs_env, matrix_s=matrix_s, error=error)
    CALL cp_dbcsr_allocate_matrix_set(matrix_w,nspin,error=error)
    DO ispin=1,nspin
       ALLOCATE(matrix_w(ispin)%matrix)
       CALL cp_dbcsr_init(matrix_w(ispin)%matrix,error=error)
       CALL cp_dbcsr_copy(matrix_w(ispin)%matrix,matrix_s(1)%matrix,&
            "W MATRIX"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)
    END DO

    CALL ep_calc_w_matrix_full(ep_env,matrix_w,error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX",&
            extension=".Log",error=error)
!FM       DO ispin=1,nspin
!FM          CALL write_sparse_matrix(matrix_w(ispin)%matrix,4,6,qs_env,globenv,output_unit=iw)
!FM       END DO
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/W_MATRIX", error=error)
    END IF

    CALL set_ks_env(ks_env,matrix_w=matrix_w,error=error)

    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " w_matrix created and added",error=error)

    !*** compute core forces (also overwrites matrix_w) ***
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " compute forces in main_qs_env",error=error)
    CALL build_core_hamiltonian_matrix(qs_env=qs_env,calculate_forces=.TRUE.,error=error)
    CALL calculate_ecore_self(qs_env,error=error)
    CALL calculate_ecore_overlap(qs_env, para_env, &
         calculate_forces=.TRUE.,&
         error=error)

    ! *** compute grid-based forces ***
    CALL qs_ks_update_qs_env(ks_env=ks_env,&
         qs_env=qs_env,&
         calculate_forces=.TRUE.,error=error)

    !  *** replicate forces ***
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " summing up forces",error=error)
    DO ikind=1,SIZE(force)
       CALL mp_sum(force(ikind)%overlap,para_env%group)
       CALL mp_sum(force(ikind)%kinetic,para_env%group)
       CALL mp_sum(force(ikind)%gth_ppl,para_env%group)
       CALL mp_sum(force(ikind)%gth_nlcc,para_env%group)
       CALL mp_sum(force(ikind)%gth_ppnl,para_env%group)
       CALL mp_sum(force(ikind)%all_potential,para_env%group)
       CALL mp_sum(force(ikind)%core_overlap,para_env%group)
       CALL mp_sum(force(ikind)%rho_core,para_env%group)
       CALL mp_sum(force(ikind)%rho_elec,para_env%group)
       CALL mp_sum(force(ikind)%vhxc_atom,para_env%group)
       CALL mp_sum(force(ikind)%g0s_Vh_elec,para_env%group)
       CALL mp_sum(force(ikind)%fock_4c,para_env%group)
       force(ikind)%total(:,:) = force(ikind)%total(:,:) +&
            force(ikind)%core_overlap(:,:) +&
            force(ikind)%gth_ppl(:,:) +&
            force(ikind)%gth_nlcc(:,:) +&
            force(ikind)%gth_ppnl(:,:) +&
            force(ikind)%all_potential(:,:) +&
            force(ikind)%kinetic(:,:) +&
            force(ikind)%overlap(:,:) +&
            force(ikind)%rho_core(:,:) +&
            force(ikind)%rho_elec(:,:) +&
            force(ikind)%vhxc_atom(:,:) +&
            force(ikind)%g0s_Vh_elec(:,:) +&
            force(ikind)%fock_4c(:,:)
    END DO

    DO iatom=1,natom

       i = atom_of_kind(iatom)
       ikind = kind_of(iatom)
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       ! the force is - dE/dR, what is called force is actually the gradient
       ! Things should have the right name
       ! The minus sign below is a hack
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       force(ikind)%other(1:3,i)=-particle_set(iatom)%f(1:3)+force(ikind)%ch_pulay(1:3,i)+ep_env%force%f0_internal(1:3,iatom)
       force(ikind)%total(1:3,i)=force(ikind)%total(1:3,i)+force(ikind)%other(1:3,i)
       particle_set(iatom)%f = -force(ikind)%total(1:3,i)
    END DO
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " finished forces",error=error)

    output_unit = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%DERIVATIVES",&
         extension=".Log",error=error)
    print_section => section_vals_get_subs_vals(qs_env%input,"DFT%PRINT%DERIVATIVES",error=error)
    CALL write_forces(force,atomic_kind_set,0,output_unit=output_unit,&
         print_section=print_section,error=error)

    CALL cp_print_key_finished_output(output_unit,logger,qs_env%input,&
         "DFT%PRINT%DERIVATIVES",error=error)

    ! destroy preconditioner
    IF (ASSOCIATED(ep_env%precond)) THEN
       !CASE(ot_precond_full_all,ot_precond_full_single, ot_precond_full_single_inverse)
       ! these depend on the ks matrix
       DO ispin=1,SIZE(ep_env%precond)
          CALL destroy_preconditioner(ep_env%precond(ispin)%preconditioner,error=error)
          DEALLOCATE(ep_env%precond(ispin)%preconditioner)
       ENDDO
       DEALLOCATE(ep_env%precond,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       !CASE(ot_precond_none,ot_precond_full_kinetic,ot_precond_s_inverse, &
       !     ot_precond_sparse_diag) ! these are 'independent'
    END IF

    CALL cp_dbcsr_deallocate_matrix_set( matrix_w, error=error )

    CALL set_ks_env(ks_env,matrix_w=matrix_w,error=error)

    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

    DEALLOCATE (kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

    CALL cp_error_dealloc_ref(suberror)
    CALL cp_print_key_log(logger,ep_env%input,"PRINT%RUN_INFO",&
         extension=".epLog",message=routineP//&
         " finished",error=error)

    CALL timestop(handle)
  END SUBROUTINE ep_env_calc_e0

! *****************************************************************************
!> \brief Calculate the W matrix for the full system
!> \param ep_env the ep environment
!> \param w_matrix sparse matrix vector that will contain the w matrix
!> \param error controls error handling
! *****************************************************************************
  SUBROUTINE ep_calc_w_matrix_full(ep_env,w_matrix,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: w_matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_calc_w_matrix_full', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, ispin
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_fm_pool_type), POINTER           :: maxao_maxmo_fm_pool
    TYPE(cp_fm_type), POINTER                :: psi0, weighted_vectors
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(weighted_vectors,f_env,maxao_maxmo_fm_pool,mos)

    CALL f_env_add_defaults(ep_env%f_env_id,f_env,new_error,failure)
    CALL mpools_get(ep_env%main_qs_env%mpools,&
         maxao_maxmo_fm_pool=maxao_maxmo_fm_pool,&
         error=error)
    CALL fm_pool_create_fm(maxao_maxmo_fm_pool,weighted_vectors,"weighted_vectors",&
         error=error)
    CALL get_qs_env(ep_env%main_qs_env,mos=mos,error=new_error)
    DO ispin=1,ep_env%nspins
       CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=psi0)
       CALL cp_fm_symm('R', 'U', ep_env%full_nao(ispin), ep_env%full_nmo(ispin), &
            -REAL(3-ep_env%nspins,dp), &
            ep_env%main_p_env%m_epsilon(ispin)%matrix,psi0, 0.0_dp,weighted_vectors,error=error)
       CALL cp_dbcsr_set(w_matrix(ispin)%matrix,0.0_dp,error=error)
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=w_matrix(ispin)%matrix,&
            matrix_v=psi0,&
            matrix_g=weighted_vectors,&
            ncol=ep_env%full_nmo(ispin),error=error)
    END DO
    CALL fm_pool_give_back_fm(maxao_maxmo_fm_pool,weighted_vectors,error=error)
    CALL f_env_rm_defaults(f_env,new_error,ierr)
    CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)

  END SUBROUTINE ep_calc_w_matrix_full

! *****************************************************************************
!> \brief calculates the derivative of dc with respect to an atomic coord
!>      with finite differences
!> \param ep_env the ep_environment
!> \param dc_dr the derivative
!> \param imol the subsystem you are interested in
!> \param iat atom to derive
!> \param idir direction to derive
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      for debugging purposes only, non symmetric derivative
! *****************************************************************************
  SUBROUTINE ep_calc_dc_dr_fdiff(ep_env,dc_dr,imol,iat,idir,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_fm_p_type), DIMENSION(:)         :: dc_dr
    INTEGER, INTENT(in)                      :: imol, iat, idir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_calc_dc_dr_fdiff', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: dr = 0.001

    INTEGER                                  :: handle, ierr, ilocal, ipos, &
                                                ispin, ndims
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: e_pot, orig
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: mo_mo_fm_pools
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_mo_m
    TYPE(cp_error_type)                      :: new_error
    REAL(kind=dp), &
      DIMENSION(ep_env%mol_envs%ndim)        :: pos
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(qs_environment_type), POINTER       :: sub_qs_env
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(qs_wf_history_type), POINTER        :: tmp_wf_history
    TYPE(qs_wf_snapshot_type), POINTER       :: snapshot, tmp_snap
    TYPE(replica_env_type), POINTER          :: rep_env

    failure=.FALSE.

    CALL timeset(routineN,handle)
    NULLIFY(snapshot,tmp_snap,mo_mo_m)
    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       rep_env => ep_env%mol_envs
       CALL f_env_add_defaults(rep_env%f_env_id,f_env,new_error,failure)
       CALL force_env_get(f_env%force_env,qs_env=sub_qs_env,&
            error=new_error)
       CALL get_qs_env(sub_qs_env,mpools=mpools,error=new_error)
       CALL mpools_get(mpools,mo_mo_fm_pools=mo_mo_fm_pools,error=new_error)

       ! ** create tmp wf history
       ilocal=rep_env_local_index(rep_env,imol,error=error)
       snapshot => wfi_get_snapshot(&
            rep_env%wf_history(ilocal)%wf_history,&
            1,error=new_error)
       CPAssert(ASSOCIATED(snapshot),cp_failure_level,routineP,new_error,failure)
       CPAssert(ASSOCIATED(snapshot%wf),cp_failure_level,routineP,new_error,failure)
       CALL wfi_create(tmp_wf_history,&
            interpolation_method_nr=wfi_use_prev_wf_method_nr,&
            extrapolation_order = 1,&
            has_unit_metric = .FALSE.,&
            error=new_error)
       CALL wfi_change_memory_depth(tmp_wf_history,&
            MAX(1,tmp_wf_history%memory_depth),error=error)
       tmp_wf_history%store_wf=.TRUE.

       CALL wfs_duplicate_snapshot(input_snapshot=snapshot,&
            output_snapshot=tmp_snap,&
            qs_env=sub_qs_env, error=new_error)
       tmp_wf_history%last_state_index=MODULO(tmp_wf_history%snapshot_count,&
            tmp_wf_history%memory_depth)+1
       CALL wfs_release(tmp_wf_history%past_states(tmp_wf_history%last_state_index)%snapshot,&
            error=new_error)
       tmp_wf_history%past_states(tmp_wf_history%last_state_index)%snapshot => &
            tmp_snap
       NULLIFY(tmp_snap)
       CALL set_qs_env(sub_qs_env,wf_history=tmp_wf_history,error=new_error)

       ndims=ep_env%mol_envs%ndim
       pos=rep_env%r(:,imol)
       IF (.not.failure) THEN
          ipos=3*(iat-1)+idir
          orig=pos(ipos)
          pos(ipos)=orig+dr
          CALL calc_energy(env_id=ep_env%mol_envs%f_env_id,pos=pos,n_el=ndims,&
               e_pot=e_pot,ierr=ierr)

       END IF
       IF (.NOT.failure) THEN
          tmp_snap => wfi_get_snapshot(&
               tmp_wf_history,&
               1,error=new_error)
          CPAssert(ASSOCIATED(tmp_snap),cp_failure_level,routineP,new_error,failure)
          CPAssert(ASSOCIATED(tmp_snap%wf),cp_failure_level,routineP,new_error,failure)
          CALL fm_pools_create_fm_vect(mo_mo_fm_pools,mo_mo_m,error=new_error)
          DO ispin=1,ep_env%nspins
                CALL cp_dbcsr_sm_fm_multiply(snapshot%overlap,&
                     tmp_snap%wf(ispin)%matrix,dc_dr(ispin)%matrix,&
                     ncol=ep_env%sub_nmo(ispin),error=new_error)
                CALL cp_gemm(transa='T',transb='N',m=ep_env%sub_nmo(ispin),&
                     n=ep_env%sub_nmo(ispin),k=ep_env%sub_nao(ispin),alpha=1.0_dp,&
                     matrix_a=snapshot%wf(ispin)%matrix,&
                     matrix_b=dc_dr(ispin)%matrix,beta=0._dp,&
                     matrix_c=mo_mo_m(ispin)%matrix,error=new_error)
                CALL cp_fm_scale_and_add(matrix_a=dc_dr(ispin)%matrix,&
                     matrix_b=tmp_snap%wf(ispin)%matrix,alpha=0.0_dp,beta=1.0_dp,&
                     error=new_error)
                CALL cp_gemm(transa='N',transb='N',m=ep_env%sub_nao(ispin),&
                     n=ep_env%sub_nmo(ispin),k=ep_env%sub_nmo(ispin),alpha=-1.0_dp/dr,&
                     matrix_a=snapshot%wf(ispin)%matrix,&
                     matrix_b=mo_mo_m(ispin)%matrix,beta=1._dp/dr,&
                     matrix_c=dc_dr(ispin)%matrix,error=new_error)
          END DO
          CALL fm_pools_give_back_fm_vect(mo_mo_fm_pools,mo_mo_m,error=new_error)
       END IF
       CALL wfi_release(tmp_wf_history,error=new_error)

       CALL f_env_rm_defaults(f_env,new_error,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE ep_calc_dc_dr_fdiff

! *****************************************************************************
!> \brief zeroes out the energies
!> \param ep_energy the energy object to zero out
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_energy_zero(ep_energy,error)
    TYPE(ep_energy_type), INTENT(out)        :: ep_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_energy_zero', &
      routineP = moduleN//':'//routineN

    ep_energy%e_no_int=0._dp
    ep_energy%e0=0._dp
    ep_energy%e1=0._dp
    ep_energy%e_tot=0._dp
  END SUBROUTINE ep_energy_zero

! *****************************************************************************
!> \brief retains an ep environment (see doc/ReferenceCounting.html)
!> \param ep_env the environment to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE ep_env_retain(ep_env, error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ep_env%ref_count=ep_env%ref_count+1
    END IF
  END SUBROUTINE ep_env_retain

! *****************************************************************************
!> \brief releases the given ep_env (see doc/ReferenceCounting.html)
!> \param ep_env the ep_env to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE ep_env_release(ep_env, error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(ep_env)) THEN
       CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
       ep_env%ref_count=ep_env%ref_count-1
       IF (ep_env%ref_count<1) THEN
          CALL rep_env_release(ep_env%mol_envs,error=error)
          IF (ASSOCIATED(ep_env%sub_proj)) THEN
             DO i=1,SIZE(ep_env%sub_proj)
                CALL p_proj_release(ep_env%sub_proj(i)%projection,error=error)
             END DO
             DEALLOCATE(ep_env%sub_proj,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          CALL qs_env_release(ep_env%main_qs_env,error=error)
          CALL p_env_release(ep_env%main_p_env,error=error)
          CALL p_env_release(ep_env%sub_p_env,error=error)
          CALL cp_fm_vect_dealloc(ep_env%m_pi_Hrho_psi0d,error=error)
          CALL cp_fm_vect_dealloc(ep_env%psi1,error=error)
          CALL destroy_force_env(ep_env%f_env_id,ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
          IF (ASSOCIATED(ep_env%precond)) THEN
             DO i=1,1 !SIZE(ep_env%precond) ! same precond for all spins
                CALL destroy_preconditioner(ep_env%precond(i)%preconditioner,error=error)
             END DO
             DEALLOCATE(ep_env%precond,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(ep_env%sub_nmo)) THEN
             DEALLOCATE(ep_env%sub_nmo,ep_env%sub_nao,ep_env%full_nmo,ep_env%full_nao,&
                  stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(ep_env%at2sub)) THEN
             DEALLOCATE(ep_env%at2sub,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          CALL globenv_release(ep_env%globenv,error=error)
          CALL cp_para_env_release(ep_env%para_env,error=error)
          CALL section_vals_release(ep_env%input,error=error)
          CALL section_vals_release(ep_env%root_section,error=error)
          CALL ep_envs_rm_ep_env(ep_env,error=error)
          CALL ep_force_release(ep_env%force,error=error)
          IF (ASSOCIATED(ep_env%base_C0)) THEN
             DO i=1,SIZE(ep_env%base_C0)
                DEALLOCATE(ep_env%base_C0(i)%array,stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             END DO
             DEALLOCATE(ep_env%base_C0,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          DEALLOCATE(ep_env, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(ep_env)
  END SUBROUTINE ep_env_release

! *****************************************************************************
!> \brief returns the various attributes of the ep_env
!> \param ep_env the ep_env to get the info from
!> \param id_nr ...
!> \param mol_envs ...
!> \param sub_proj ...
!> \param main_qs_env ...
!> \param main_p_env ...
!> \param sub_p_env ...
!> \param preconditioner ...
!> \param m_pi_Hrho_psi0d ...
!> \param psi1 ...
!> \param f_env_id ...
!> \param globenv ...
!> \param at2sub ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For the other arguments see the attributes of ep_env_type
!> \par History
!>      12.2002 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_env_get(ep_env,id_nr, mol_envs, sub_proj, main_qs_env,&
       main_p_env,sub_p_env, preconditioner, m_pi_Hrho_psi0d, psi1, f_env_id,&
       globenv,at2sub,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    INTEGER, INTENT(out), OPTIONAL           :: id_nr
    TYPE(replica_env_type), OPTIONAL, &
      POINTER                                :: mol_envs
    TYPE(qs_p_projection_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sub_proj
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: main_qs_env
    TYPE(qs_p_env_type), OPTIONAL, POINTER   :: main_p_env, sub_p_env
    TYPE(preconditioner_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: preconditioner
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: m_pi_Hrho_psi0d, psi1
    INTEGER, INTENT(out), OPTIONAL           :: f_env_id
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: at2sub
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(ep_env),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CPPrecondition(ep_env%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT. failure) THEN
       IF (PRESENT(id_nr)) id_nr=ep_env%id_nr
       IF (PRESENT(mol_envs)) mol_envs => ep_env%mol_envs
       IF (PRESENT(sub_proj)) sub_proj => ep_env%sub_proj
       IF (PRESENT(main_qs_env)) main_qs_env => ep_env%main_qs_env
       IF (PRESENT(main_p_env)) main_p_env => ep_env%main_p_env
       IF (PRESENT(sub_p_env)) sub_p_env => ep_env%sub_p_env
       IF (PRESENT(m_pi_Hrho_psi0d)) m_pi_Hrho_psi0d => ep_env%m_pi_Hrho_psi0d
       IF (PRESENT(psi1)) psi1 => ep_env%psi1
!FM     IF (PRESENT(nspins)) THEN
!FM        CPPrecondition(ASSOCIATED(ep_env%main_qs_env),cp_failure_level,routineP,error,failure)
!FM        IF (.not.failure) THEN
!FM           nspins=ep_env%main_qs_env%dft_control%nspins
!FM        END IF
!FM     END IF
       IF (PRESENT(preconditioner)) THEN
          preconditioner => ep_env%precond
       END IF
       IF (PRESENT(globenv)) globenv => ep_env%globenv
       IF (PRESENT(f_env_id)) f_env_id=ep_env%f_env_id
       IF (PRESENT(at2sub)) at2sub => ep_env%at2sub
    END IF
  END SUBROUTINE ep_env_get

! *****************************************************************************
!> \brief returns the replica environment with the given id_nr
!> \param id_nr the id_nr of the requested ep_envs
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION ep_envs_get_ep_env(id_nr) RESULT(res)
    INTEGER, INTENT(in)                      :: id_nr
    TYPE(ep_env_type), POINTER               :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_envs_get_ep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    NULLIFY(res)
    IF (module_initialized) THEN
       IF (ASSOCIATED(ep_envs)) THEN
          DO i=1,SIZE(ep_envs)
             IF (ep_envs(i)%ep_env%id_nr==id_nr) THEN
                res => ep_envs(i)%ep_env
                EXIT
             END IF
          END DO
       END IF
    END IF
  END FUNCTION ep_envs_get_ep_env

! *****************************************************************************
!> \brief adds the given ep_env to the list of controlled ep_envs.
!> \param ep_env the ep_env to add
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_envs_add_ep_env(ep_env,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_envs_add_ep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(ep_env_p_type), DIMENSION(:), &
      POINTER                                :: new_ep_envs
    TYPE(ep_env_type), POINTER               :: ep_env2

    failure=.FALSE.

    IF (ASSOCIATED(ep_env)) THEN
       ep_env2 => ep_envs_get_ep_env(ep_env%id_nr)
       IF (.NOT.ASSOCIATED(ep_env2)) THEN
          IF (module_initialized) THEN
             IF (ASSOCIATED(ep_envs)) THEN
                ALLOCATE(new_ep_envs(SIZE(ep_envs)+1),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                DO i=1,SIZE(ep_envs)
                   new_ep_envs(i)%ep_env => ep_envs(i)%ep_env
                END DO
                DEALLOCATE(ep_envs,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
                ep_envs => new_ep_envs
             ELSE
                ALLOCATE(ep_envs(1),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
             END IF
          ELSE
             ALLOCATE(ep_envs(1),stat=stat)
             CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          END IF
          ep_envs(SIZE(ep_envs))%ep_env => ep_env
          module_initialized=.TRUE.
       END IF
    END IF
  END SUBROUTINE ep_envs_add_ep_env

! *****************************************************************************
!> \brief removes the given ep_env to the list of controlled ep_envs.
!> \param ep_env the ep_env to remove
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_envs_rm_ep_env(ep_env,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_envs_rm_ep_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ii, stat
    LOGICAL                                  :: failure
    TYPE(ep_env_p_type), DIMENSION(:), &
      POINTER                                :: new_ep_envs

    failure=.FALSE.

    IF (ASSOCIATED(ep_env)) THEN
       CPPrecondition(module_initialized,cp_failure_level,routineP,error,failure)
       ALLOCATE(new_ep_envs(SIZE(ep_envs)-1),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       ii=0
       DO i=1,SIZE(ep_envs)
          IF (ep_envs(i)%ep_env%id_nr/=ep_env%id_nr) THEN
             ii=ii+1
             new_ep_envs(ii)%ep_env => ep_envs(i)%ep_env
          END IF
       END DO
       CPPostcondition(ii==SIZE(new_ep_envs),cp_failure_level,routineP,error,failure)
       DEALLOCATE(ep_envs,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       ep_envs => new_ep_envs
       IF (SIZE(ep_envs)==0) THEN
          DEALLOCATE(ep_envs,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE ep_envs_rm_ep_env

! *****************************************************************************
!> \brief updates forces and energy in the given ep_env
!> \param ep_env ...
!> \param calc_f ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_env_calc_e_f(ep_env,calc_f,error)
    TYPE(ep_env_type), POINTER               :: ep_env
    LOGICAL, INTENT(in)                      :: calc_f
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_env_calc_e_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, my_calc_f
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.failure) THEN
       IF (calc_f) THEN
          my_calc_f=1
       ELSE
          my_calc_f=0
       END IF
       CALL ep_env_calc_e_f_low(ep_env%id_nr,my_calc_f,ierr)
       CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE ep_env_calc_e_f

! *****************************************************************************
!> \brief creates an ep_force_type
!> \param force the type to create
!> \param nat the number of atoms
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_force_create(force,nat,error)
    TYPE(ep_force_type), POINTER             :: force
    INTEGER, INTENT(in)                      :: nat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_force_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(force),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       ALLOCATE(force,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       last_force_id=last_force_id+1
       force%id_nr=last_force_id
       force%ref_count=1
       ALLOCATE(force%f0_internal(3,nat),stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
  END SUBROUTINE ep_force_create

! *****************************************************************************
!> \brief zeroes the forces
!> \param force the force type to zero out
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_force_zero(force,error)
    TYPE(ep_force_type), POINTER             :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_force_zero', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(force),cp_failure_level,routineP,error,failure)
    CPPrecondition(force%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       force%f0_internal=0._dp
    END IF
  END SUBROUTINE ep_force_zero

! *****************************************************************************
!> \brief retains the given force type
!> \param force the force type to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_force_retain(force,error)
    TYPE(ep_force_type), POINTER             :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_force_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(force),cp_failure_level,routineP,error,failure)
    CPPrecondition(force%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       force%ref_count=force%ref_count+1
    END IF
  END SUBROUTINE ep_force_retain

! *****************************************************************************
!> \brief retains the given force type
!> \param force the force type to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE ep_force_release(force,error)
    TYPE(ep_force_type), POINTER             :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ep_force_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(force)) THEN
       CPPrecondition(force%ref_count>0,cp_failure_level,routineP,error,failure)
       force%ref_count=force%ref_count-1
       IF (force%ref_count==0) THEN
          DEALLOCATE(force%f0_internal,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(force,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(force)
  END SUBROUTINE ep_force_release

! *****************************************************************************
! The following routines are never called, but introduce a dependence on the
! old SM code.
!> \brief transfers a full sparse matrix to local matrix
!> \param ep_env the ep_environment
!> \param full_sm the full matrix to copy
!> \param local_sm the targt local matrix
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
! SUBROUTINE ep_full_to_local_sm(ep_env, full_sm, local_sm, ispin, error)
!   TYPE(ep_env_type), POINTER               :: ep_env
!   TYPE(real_matrix_type), POINTER          :: full_sm
!   TYPE(real_matrix_p_type), &
!     DIMENSION(:, :), POINTER               :: local_sm
!   INTEGER                                  :: ispin
!   TYPE(cp_error_type), INTENT(inout)       :: error
!   CHARACTER(len=*), PARAMETER :: routineN = 'ep_full_to_local_sm', &
!     routineP = moduleN//':'//routineN
!   INTEGER                                  :: disp, i, iblock_col, &
!                                               iblock_row, ii, ip, irep, j, &
!                                               lcol, lrow, num_pe, stat, &
!                                               target_p
!   INTEGER, ALLOCATABLE, DIMENSION(:)       :: pos_att, rcv_sizes, rdisp, &
!                                               sdisp, send_sizes
!   INTEGER, DIMENSION(:), POINTER           :: atoms, first_col, first_row, &
!                                               last_col, last_row
!   LOGICAL                                  :: failure
!   REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: rcv_buf, send_buf
!   REAL(kind=dp), DIMENSION(:, :), POINTER  :: m_block
!   TYPE(real_block_node_type), POINTER      :: block_node
!   TYPE(replica_env_type), POINTER          :: rep_env
!   failure=.FALSE.
!   NULLIFY(m_block)
!   CPPrecondition(ASSOCIATED(full_sm),cp_failure_level,routineP,error,failure)
!   CPPrecondition(ASSOCIATED(local_sm),cp_failure_level,routineP,error,failure)
!   IF (.NOT.failure) THEN
!      CPPrecondition(SIZE(local_sm,1)<=ispin,cp_failure_level,routineP,error,failure)
!   END IF
!   IF (.NOT. failure) THEN
!      rep_env => ep_env%mol_envs
!      num_pe=rep_env%para_env_inter_rep%num_pe
!      ALLOCATE(send_sizes(0:num_pe-1),stat=stat)
!      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
!      send_sizes=0
!      DO iblock_row=1,full_sm%nblock_row
!         block_node => first_block_node(full_sm,iblock_row)
!         DO WHILE (ASSOCIATED(block_node))
!            CALL get_block_node(block_node=block_node,&
!                 block_col=iblock_col, block=m_block)
!            IF (ep_env%at2sub(iblock_col)==ep_env%at2sub(iblock_row)) THEN
!               target_p=rep_env%inter_rep_rank(&
!                    rep_env%replica_owner(ep_env%at2sub(iblock_col)))
!               send_sizes(target_p)=send_sizes(target_p)+SIZE(m_block,1)*SIZE(m_block,2)+2
!            END IF
!            block_node => next_block_node(block_node)
!         END DO
!      END DO
!      ALLOCATE(sdisp(0:num_pe), stat=stat)
!      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
!      disp=0
!      sdisp(0)=0
!      DO ip=1,num_pe
!         disp=disp+send_sizes(ip-1)
!         sdisp(ip)=disp
!      END DO
!      ALLOCATE(send_buf(disp),stat=stat)
!      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
!      ALLOCATE(rcv_sizes(0:num_pe-1),stat=stat)
!      CALL mp_alltoall(send_sizes,rcv_sizes,1,num_pe)
!      ALLOCATE(rdisp(0:num_pe),stat=stat)
!      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
!      disp=0
!      rdisp(0)=0
!      DO ip=1,num_pe
!         disp=disp+rcv_sizes(ip-1)
!         rdisp(ip)=disp
!      END DO
!      ALLOCATE(rcv_buf(disp),stat=stat)
!      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
!      ALLOCATE(pos_att(0:num_pe-1),stat=stat)
!      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
!      pos_att=sdisp(0:num_pe-1)
!      DO iblock_row=1,full_sm%nblock_row
!         block_node => first_block_node(full_sm,iblock_row)
!         DO WHILE (ASSOCIATED(block_node))
!            CALL get_block_node(block_node=block_node,&
!                 block_col=iblock_col,&
!                 BLOCK=m_block)
!            IF (ep_env%at2sub(iblock_col)==ep_env%at2sub(iblock_row)) THEN
!               target_p=rep_env%inter_rep_rank(&
!                    rep_env%replica_owner(ep_env%at2sub(iblock_col)))
!               send_buf(pos_att(target_p)+1)=REAL(iblock_row,dp)
!               send_buf(pos_att(target_p)+2)=REAL(iblock_col,dp)
!               ii=pos_att(target_p)+2
!               DO j=1,SIZE(m_block,2)
!                  DO i=1,SIZE(m_block,1)
!                     ii=ii+1
!                     send_buf(ii)=m_block(i,j)
!                  END DO
!               END DO
!               pos_att(target_p)=ii
!            END IF
!            block_node => next_block_node(block_node)
!         END DO
!      END DO
!      CPPostcondition(ALL(sdisp(1:)==pos_att),cp_failure_level,routineP,error,failure)
!      CALL mp_alltoall ( send_buf, send_sizes, sdisp,&
!           rcv_buf, rcv_sizes, rdisp, rep_env%para_env_inter_rep%group )
!      CALL get_matrix_info(full_sm,&
!           first_row=first_row,last_row=last_row,&
!           first_col=first_col,last_col=last_col)
!      pos_att=rdisp(0:num_pe-1)
!      DO ip=0,num_pe-1
!         DO
!            IF (pos_att(ip)>=rdisp(ip+1)) EXIT
!            iblock_row=INT(rcv_buf(pos_att(ip)+1))
!            iblock_col=INT(rcv_buf(pos_att(ip)+2))
!FM             PRINT *,"ip",ip," pos_att(ip)",pos_att(ip),&
!FM                  " rdisp(ip+1)",rdisp(ip+1)," iblock_col",iblock_col,&
!FM                  " iblock_row",iblock_row
!            irep=ep_env%at2sub(iblock_row)
!            CPAssert(ep_env%at2sub(iblock_col)==irep,cp_failure_level,routineP,error,failure)
!            ! allow for different basis in full vs sub system, i.e. use proj?
!            ALLOCATE(m_block(last_row(iblock_row)-first_row(iblock_row)+1,&
!                 last_col(iblock_col)-first_col(iblock_col)+1),stat=stat)
!            ii=pos_att(ip)+2
!            DO j=1,SIZE(m_block,2)
!               DO i=1,SIZE(m_block,1)
!                  ii=ii+1
!                  m_block(i,j)=rcv_buf(ii)
!               END DO
!            END DO
!            pos_att(ip)=ii
!            atoms => ep_env%sub_proj(irep)%projection%atoms
!            lrow=iblock_row-atoms(1)+1
!            lcol=iblock_col-atoms(1)+1
!            CPAssert(lrow>0,cp_failure_level,routineP,error,failure)
!            CPAssert(lrow<=SIZE(atoms),cp_failure_level,routineP,error,failure)
!            CPAssert(lcol>0,cp_failure_level,routineP,error,failure)
!            CPAssert(lcol<=SIZE(atoms),cp_failure_level,routineP,error,failure)
!            CPAssert(atoms(lrow)==iblock_row,cp_failure_level,routineP,error,failure)
!            CPAssert(atoms(lcol)==iblock_col,cp_failure_level,routineP,error,failure)
!            CALL add_block_node(matrix=local_sm(ispin,&
!                 rep_env_local_index(rep_env,irep,error=error))%matrix,&
!                 block_row=lrow,&
!                 block_col=lcol,&
!                 BLOCK=m_block,error=error)
!            DEALLOCATE(m_block,stat=stat)
!            CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
!         END DO
!      END DO
!   END IF
! END SUBROUTINE ep_full_to_local_sm



































































END MODULE ep_methods
