!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for all ALMO-based SCF methods
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin 
! *****************************************************************************
MODULE almo_scf
  USE almo_scf_aux2_methods,           ONLY: construct_object01,&
                                             copy_object01_gen,&
                                             init_object01_gen,&
                                             op2_object01_gen,&
                                             op3_object01,&
                                             release_object01_gen
  USE almo_scf_aux2_types,             ONLY: object01_type,&
                                             select_row
  USE almo_scf_aux_types,              ONLY: extrapolate_object00,&
                                             init_object00_gen,&
                                             object00_type,&
                                             push_object00,&
                                             release_object00
  USE almo_scf_methods,                ONLY: &
       almo_scf_get_hp_blk_and_tv_blk, almo_scf_get_hp_xx_and_tv_xx, &
       almo_scf_get_t_blk, almo_scf_loc_hp_blk, almo_scf_loc_hp_xx, &
       almo_scf_ortho_blk, almo_scf_p_get_t, almo_scf_p_get_t_blk, &
       get_group_complex, get_group_distr, get_group_inv, get_group_rdown, &
       get_group_sqrt, invert_blk_once, operations_on_sets
  USE almo_scf_special,                ONLY: almo_level2_spec1_0,&
                                             almo_level2_spec2_1,&
                                             almo_levelX_spec2_0,&
                                             almo_levelX_spec3_0,&
                                             almo_levelX_spec5_0,&
                                             almo_levelX_spec6_0
  USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                             almo_mat_dim_occ,&
                                             almo_mat_dim_virt,&
                                             almo_mat_dim_virt_disc,&
                                             almo_mat_dim_virt_full,&
                                             almo_max_cutoff_multiplier,&
                                             almo_objectM1_type
  USE bibliography,                    ONLY: Khaliullin2013,&
                                             cite_reference
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                             cp_dbcsr_cholesky_invert
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_function_of_elements, cp_dbcsr_get_info, &
       cp_dbcsr_get_stored_coordinates, cp_dbcsr_hadamard_product, &
       cp_dbcsr_init, cp_dbcsr_multiply, cp_dbcsr_nblkcols_total, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_norm, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_reserve_block2d, cp_dbcsr_scale, &
       cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_type, cp_dbcsr_work_create, &
       dbcsr_distribution_mp, dbcsr_func_artanh, dbcsr_func_dtanh, &
       dbcsr_func_inverse, dbcsr_func_tanh, dbcsr_mp_mynode, &
       dbcsr_norm_maxabsnorm, dbcsr_type_no_symmetry, dbcsr_type_symmetric
  USE cp_external_control,             ONLY: external_control
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE input_constants,                 ONLY: &
       almo_constraint_distance, almo_deloc_full_scf, almo_deloc_none, &
       almo_deloc_qscf, almo_deloc_qx, almo_deloc_x, &
       almo_deloc_x_then_full_scf, almo_deloc_xk, almo_domain_layout_atomic, &
       almo_domain_layout_molecular, almo_mat_distr_atomic, &
       almo_mat_distr_molecular, almo_scf_diag, almo_scf_dm_sign, &
       almo_scf_pcg, cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, &
       cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, &
       cg_zero, do_bondparm_vdw, tensor_orthogonal, virt_full, virt_minimal, &
       virt_number, virt_occ_size
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                             matrix_sqrt_Newton_Schulz
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE molecule_types_new,              ONLY: get_molecule_set_info,&
                                             molecule_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf'

  PUBLIC :: almo_entry_scf

  LOGICAL, PARAMETER :: debug_mode = .FALSE.
  LOGICAL, PARAMETER :: safe_mode = .FALSE.

CONTAINS

! *****************************************************************************
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_entry_scf(qs_env, calc_forces, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: calc_forces
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_entry_scf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(almo_objectM1_type)                 :: main_objectM1

    CALL timeset(routineN,handle)

    CALL cite_reference(Khaliullin2013)

    ! no forces so far 
    CALL cp_assert(.NOT.calc_forces,cp_failure_level,cp_assertion_failed,routineP,&
           " No forces/gradients with almo.... yet(?!) "//&
                CPSourceFileRef, only_ionode=.TRUE.)

    CALL almo_level1_exp1_0(qs_env,main_objectM1,error)
    CALL almo_level1_exp2_1(qs_env,main_objectM1,error)
    CALL almo_level1_exp3_2(qs_env,main_objectM1,error)
    CALL almo_level1_exp4_1(qs_env,main_objectM1,error)

    CALL timestop(handle)

  END SUBROUTINE almo_entry_scf

! *****************************************************************************
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level1_exp1_0(qs_env,main_objectM1,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level1_exp1_0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ao, handle, i, iao, idomain, &
                                                ispin, naos, natoms, &
                                                ndomains, nmols, nspins, &
                                                unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
    failure = .FALSE.

    ! define the output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    CALL get_qs_env(qs_env,&
      nelectron_total=main_objectM1%nelectrons_total,&
      matrix_s=matrix_s,&
      dft_control=dft_control,&
      molecule_set=molecule_set,&
      input=input,&
      has_unit_metric=main_objectM1%orthogonal_basis,&
      para_env=main_objectM1%para_env,&
      blacs_env=main_objectM1%blacs_env,&
      nelectron_spin=main_objectM1%nelectrons_spin,&
      error=error)
    CALL cp_para_env_retain(main_objectM1%para_env,error)
    CALL cp_blacs_env_retain(main_objectM1%blacs_env, error)

    main_objectM1%nspins=dft_control%nspins
    main_objectM1%natoms=cp_dbcsr_nblkrows_total(matrix_s(1)%matrix)
    main_objectM1%nmolecules=SIZE(molecule_set)
    CALL cp_dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=naos )
    main_objectM1%naos=naos

    nspins=main_objectM1%nspins
    nmols=main_objectM1%nmolecules
    natoms=main_objectM1%natoms
    
    CALL almo_level2_aux3_p1(input,main_objectM1,unit_nr,error)
   
    IF (main_objectM1%domain_layout_mos==almo_domain_layout_molecular) THEN
       main_objectM1%ndomains=main_objectM1%nmolecules
    ELSE
       main_objectM1%ndomains=main_objectM1%natoms
    ENDIF

    ndomains=main_objectM1%ndomains
    ALLOCATE(main_objectM1%domain_index_of_atom(natoms))
    ALLOCATE(main_objectM1%domain_index_of_ao(naos))
    ALLOCATE(main_objectM1%first_atom_of_domain(ndomains))
    ALLOCATE(main_objectM1%nbasis_of_domain(ndomains))
    ALLOCATE(main_objectM1%nocc_of_domain(ndomains,nspins))
    ALLOCATE(main_objectM1%nvirt_full_of_domain(ndomains,nspins))
    ALLOCATE(main_objectM1%nvirt_of_domain(ndomains,nspins))
    ALLOCATE(main_objectM1%nvirt_disc_of_domain(ndomains,nspins))
    ALLOCATE(main_objectM1%mu_of_domain(ndomains,nspins))
    ALLOCATE(main_objectM1%cpu_of_domain(ndomains))
    ALLOCATE(main_objectM1%charge_of_domain(ndomains))

    IF (main_objectM1%domain_layout_mos==almo_domain_layout_molecular) THEN
       CALL get_molecule_set_info(molecule_set,&
               atom_to_mol=main_objectM1%domain_index_of_atom,&
               mol_to_first_atom=main_objectM1%first_atom_of_domain,&
               mol_to_nelectrons=main_objectM1%nocc_of_domain,&
               mol_to_nbasis=main_objectM1%nbasis_of_domain,&
               error=error)
       DO ispin=1,nspins
          main_objectM1%nvirt_full_of_domain(:,ispin)=&
             main_objectM1%nbasis_of_domain(:)-&
             main_objectM1%nocc_of_domain(:,ispin)
          SELECT CASE (main_objectM1%deloc_truncate_virt)
          CASE (virt_full)
             main_objectM1%nvirt_of_domain(:,ispin)=&
                main_objectM1%nvirt_full_of_domain(:,ispin)
             main_objectM1%nvirt_disc_of_domain(:,ispin)=0
          CASE (virt_number)
             DO idomain=1,ndomains
                main_objectM1%nvirt_of_domain(idomain,ispin)=&
                   MIN(main_objectM1%deloc_virt_per_domain,&
                   main_objectM1%nvirt_full_of_domain(idomain,ispin))
                main_objectM1%nvirt_disc_of_domain(idomain,ispin)=&
                   main_objectM1%nvirt_full_of_domain(idomain,ispin)-&
                   main_objectM1%nvirt_of_domain(idomain,ispin)
             ENDDO
          CASE (virt_occ_size)
             DO idomain=1,ndomains
                main_objectM1%nvirt_of_domain(idomain,ispin)=&
                   MIN(main_objectM1%nocc_of_domain(idomain,ispin),&
                   main_objectM1%nvirt_full_of_domain(idomain,ispin))
                main_objectM1%nvirt_disc_of_domain(idomain,ispin)=&
                   main_objectM1%nvirt_full_of_domain(idomain,ispin)-&
                   main_objectM1%nvirt_of_domain(idomain,ispin)
             ENDDO
          CASE DEFAULT
             CPErrorMessage(cp_failure_level,routineP,"illegal method for virtual space truncation",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       ENDDO 
    ELSE 
       main_objectM1%domain_index_of_atom(:)=(/(i, i=1,natoms)/)
    ENDIF
    
    ao=1
    DO idomain=1,ndomains
      DO iao=1,main_objectM1%nbasis_of_domain(idomain)
         main_objectM1%domain_index_of_ao(ao)=idomain
         ao=ao+1
      ENDDO
    ENDDO

    main_objectM1%mu_of_domain(:,:)=main_objectM1%mu

    IF (main_objectM1%mat_distr_aos==almo_mat_distr_atomic) THEN
     ALLOCATE(main_objectM1%domain_index_of_ao_block(natoms))
     main_objectM1%domain_index_of_ao_block(:)=&
        main_objectM1%domain_index_of_atom(:)
    ELSE IF (main_objectM1%mat_distr_aos==almo_mat_distr_molecular) THEN
     ALLOCATE(main_objectM1%domain_index_of_ao_block(nmols))
     main_objectM1%domain_index_of_ao_block(:)=(/(i, i=1,nmols)/)
    ENDIF
    IF (main_objectM1%mat_distr_mos==almo_mat_distr_atomic) THEN
     ALLOCATE(main_objectM1%domain_index_of_mo_block(natoms))
     main_objectM1%domain_index_of_mo_block(:)=&
        main_objectM1%domain_index_of_atom(:)
    ELSE IF (main_objectM1%mat_distr_mos==almo_mat_distr_molecular) THEN
     ALLOCATE(main_objectM1%domain_index_of_mo_block(nmols))
     main_objectM1%domain_index_of_mo_block(:)=(/(i, i=1,nmols)/)
    ENDIF

       main_objectM1%need_previous_ks=.TRUE.

       main_objectM1%need_virtuals=.TRUE.
       main_objectM1%need_orbital_energies=.TRUE.

    CALL almo_level3_aux1_0(main_objectM1,unit_nr,error)

    CALL almo_level2_aux1(main_objectM1,matrix_s(1)%matrix,error)
    main_objectM1%s_inv_done=.FALSE.
    main_objectM1%s_sqrt_done=.FALSE.
    CALL almo_level2_aux0_0(matrix_s(1)%matrix,main_objectM1,error)

    CALL almo_level2_spec1_0(qs_env,main_objectM1,error)
    CALL get_group_distr(main_objectM1,error)

    ALLOCATE(main_objectM1%concept_p(ndomains,nspins))
    CALL init_object01_gen(main_objectM1%concept_p,error)
   
    ALLOCATE(main_objectM1%concept_hfh_xx(ndomains,nspins))
    CALL init_object01_gen(main_objectM1%concept_hfh_xx,error)
   
    ALLOCATE(main_objectM1%concept_soi(ndomains,nspins))
    CALL init_object01_gen(main_objectM1%concept_soi,error)
    ALLOCATE(main_objectM1%concept_ss_inv(ndomains,nspins))
    CALL init_object01_gen(main_objectM1%concept_ss_inv,error)
    ALLOCATE(main_objectM1%concept_ss(ndomains,nspins))
    CALL init_object01_gen(main_objectM1%concept_ss,error)
    ALLOCATE(main_objectM1%concept_enter(ndomains,nspins))
    CALL init_object01_gen(main_objectM1%concept_enter,error=error)
    ALLOCATE(main_objectM1%concept_gst(ndomains,nspins))
    CALL init_object01_gen(main_objectM1%concept_gst,error=error)
    ALLOCATE(main_objectM1%concept_dpp(ndomains,nspins))
    CALL init_object01_gen(main_objectM1%concept_dpp,error=error)

    CALL almo_level2_spec2_1(qs_env,main_objectM1,error)

    CALL timestop(handle)

  END SUBROUTINE almo_level1_exp1_0

! *****************************************************************************
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level2_aux3_p1(input,main_objectM1,unit_nr,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level2_aux3_p1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: almo_scf_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    almo_scf_section => section_vals_get_subs_vals(input,"DFT%ALMO_SCF",error=error)

    ! read user input
    CALL section_vals_val_get(almo_scf_section,"EPS_FILTER",&
                              r_val=main_objectM1%eps_filter,error=error)
    CALL section_vals_val_get(almo_scf_section,"BLOCKED_EPS_ITER",&
                              r_val=main_objectM1%eps_scf_bd,error=error)
    CALL section_vals_val_get(almo_scf_section,"BLOCKED_MAX_ITER",&
                              i_val=main_objectM1%max_scf_bd,error=error)
    CALL section_vals_val_get(almo_scf_section,"BLOCKED_N_DIIS",&
                              i_val=main_objectM1%counter04_bd,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_EPS_ITER",&
                              r_val=main_objectM1%eps_scf_q,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_EPS_LIN_SEARCH",&
                              r_val=main_objectM1%eps_lin_search,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_MAX_ITER",&
                              i_val=main_objectM1%max_scf_q,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_CONJUGATOR",&
                              i_val=main_objectM1%scf_conjugator,&
                              error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_METHOD",&
                              i_val=main_objectM1%deloc_method,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_ALGORITHM",&
                              i_val=main_objectM1%almo_update_algorithm_q,&
                              error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_R_CUTOFF_FACTOR",&
                              r_val=main_objectM1%quencher_r0_factor,&
                              error=error)

    main_objectM1%domain_layout_aos=almo_domain_layout_molecular
    main_objectM1%domain_layout_mos=almo_domain_layout_molecular
    main_objectM1%mat_distr_aos=almo_mat_distr_molecular 
    main_objectM1%mat_distr_mos=almo_mat_distr_molecular

    main_objectM1%outer_max_scf_q=0
    main_objectM1%constraint_type=almo_constraint_distance
    main_objectM1%mu=-0.1_dp
    main_objectM1%fixed_mu=.FALSE.
    main_objectM1%almo_update_algorithm=almo_scf_diag
    main_objectM1%eps_prev_guess=main_objectM1%eps_filter/10.0_dp
    main_objectM1%mixing_fraction=0.45_dp
    
    main_objectM1%deloc_cayley_tensor_type=tensor_orthogonal
    main_objectM1%deloc_cayley_conjugator=cg_hager_zhang
    main_objectM1%deloc_cayley_max_iter=100
    main_objectM1%deloc_use_occ_orbs=.TRUE.
    main_objectM1%deloc_cayley_use_virt_orbs=.FALSE.
    main_objectM1%deloc_cayley_linear=.FALSE.
    main_objectM1%deloc_cayley_eps_convergence=1.0E-6_dp
    main_objectM1%deloc_cayley_occ_precond=.TRUE.
    main_objectM1%deloc_cayley_vir_precond=.TRUE.
    main_objectM1%deloc_truncate_virt=virt_full
    main_objectM1%deloc_virt_per_domain=-1
    
    main_objectM1%opt_k_eps_convergence=1.0E-5_dp
    main_objectM1%opt_k_max_iter=100
    main_objectM1%opt_k_outer_max_iter=1
    main_objectM1%opt_k_trial_step_size=0.05_dp
    main_objectM1%opt_k_conjugator=cg_hager_zhang
    main_objectM1%opt_k_trial_step_size_multiplier=1.05_dp
    main_objectM1%opt_k_conj_iter_start=0
    main_objectM1%opt_k_prec_iter_start=0
    main_objectM1%opt_k_conj_iter_freq=10000000
    main_objectM1%opt_k_prec_iter_freq=10000000
    
    main_objectM1%quencher_radius_type=do_bondparm_vdw
    main_objectM1%quencher_r1_factor=main_objectM1%quencher_r0_factor
    
    main_objectM1%quencher_s0=1.0E-4_dp
    main_objectM1%quencher_s1=1.0E-6_dp

    main_objectM1%envelope_amplitude=1.0_dp

    main_objectM1%logical01=.FALSE. 
    main_objectM1%logical02=.TRUE. 
    main_objectM1%logical03=.TRUE. 
    main_objectM1%logical04=.TRUE. 
    main_objectM1%logical05=.FALSE.

    main_objectM1%real01=main_objectM1%eps_filter/10.0_dp
    main_objectM1%real02=0.0_dp 
    main_objectM1%real03=0.0_dp
    main_objectM1%real04=0.5_dp
    main_objectM1%real05=0.1_dp

    main_objectM1%integer01=10
    main_objectM1%integer02=4 
    main_objectM1%integer03=0
    main_objectM1%integer04=0
    main_objectM1%integer05=0
          
    IF (main_objectM1%almo_update_algorithm_q.eq.almo_scf_diag .AND. &
        main_objectM1%deloc_method.ne.almo_deloc_qx) THEN
      CPErrorMessage(cp_failure_level,routineP,"D-DIAG is implemented only for X_R",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (main_objectM1%deloc_truncate_virt.EQ.virt_number .AND. &
        main_objectM1%deloc_virt_per_domain.LE.0) THEN
      CPErrorMessage(cp_failure_level,routineP,"specify a positive number of virtual orbitals",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (main_objectM1%deloc_truncate_virt.EQ.virt_minimal) THEN
      CPErrorMessage(cp_failure_level,routineP,"VIRT TRUNCATION TO MINIMAL BASIS IS NIY",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (main_objectM1%domain_layout_mos.NE.almo_domain_layout_molecular) THEN
      CPErrorMessage(cp_failure_level,routineP,"use MOLECULAR domains",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (main_objectM1%domain_layout_aos.NE.almo_domain_layout_molecular) THEN
      CPErrorMessage(cp_failure_level,routineP,"use MOLECULAR domains",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (main_objectM1%mat_distr_mos.NE.almo_mat_distr_molecular) THEN
      CPErrorMessage(cp_failure_level,routineP,"use MOLECULAR distr for MOs",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (main_objectM1%mat_distr_aos==almo_mat_distr_molecular .AND. &
        main_objectM1%domain_layout_aos==almo_domain_layout_atomic) THEN 
       CPErrorMessage(cp_failure_level,routineP,"AO blocks cannot be larger than domains",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (main_objectM1%mat_distr_mos==almo_mat_distr_molecular .AND. &
        main_objectM1%domain_layout_mos==almo_domain_layout_atomic) THEN 
       CPErrorMessage(cp_failure_level,routineP,"MO blocks cannot be larger than domains",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (main_objectM1%quencher_r1_factor.gt.almo_max_cutoff_multiplier) THEN
       CPErrorMessage(cp_failure_level,routineP,"DELOCALIZE_R_CUTOFF_FACTOR is larger than almo_max_cutoff_multiplier",error)
       CPErrorMessage(cp_failure_level,routineP,"increase hard-coded almo_max_cutoff_multiplier",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_level2_aux3_p1

! *****************************************************************************
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level3_aux1_0(main_objectM1,unit_nr,error)
    
    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level3_aux1_0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'()')
       WRITE(unit_nr,'(T2,A,A,A)') REPEAT("-",32)," ALMO SETTINGS ",REPEAT("-",32)
    
       WRITE(unit_nr,'(T2,A,T48,E33.3)') "eps_filter:",main_objectM1%eps_filter
       WRITE(unit_nr,'(T2,A,T48,E33.3)') "blocked_eps_iter:",main_objectM1%eps_scf_bd
       WRITE(unit_nr,'(T2,A,T48,I33)') "blocked_max_iter:",main_objectM1%max_scf_bd
       WRITE(unit_nr,'(T2,A,T48,I33)') "blocked_n_diis:",main_objectM1%counter04_bd
       
       SELECT CASE(main_objectM1%deloc_method)
       CASE(almo_deloc_none)
           WRITE(unit_nr,'(T2,A,T48,A33)') "delocalization:","NONE"
       CASE(almo_deloc_x)
           WRITE(unit_nr,'(T2,A,T48,A33)') "delocalization:","X"
       CASE(almo_deloc_full_scf)
           WRITE(unit_nr,'(T2,A,T48,A33)') "delocalization:","SCF"
       CASE(almo_deloc_qx)
           WRITE(unit_nr,'(T2,A,T48,A33)') "delocalization:","X_R"
       CASE(almo_deloc_qscf)
           WRITE(unit_nr,'(T2,A,T48,A33)') "delocalization:","SCF_R"
       CASE(almo_deloc_x_then_full_scf)
           WRITE(unit_nr,'(T2,A,T48,A33)') "delocalization:","X_THEN_SCF"
       END SELECT
       
       IF (main_objectM1%deloc_method.ne.almo_deloc_none) THEN
          
          SELECT CASE(main_objectM1%almo_update_algorithm_q)
          CASE(almo_scf_diag)
              WRITE(unit_nr,'(T2,A,T48,A33)') "algorithm:","Blocked diagonalization"
          CASE(almo_scf_pcg)
              WRITE(unit_nr,'(T2,A,T48,A33)') "algorithm:","PCG"
          END SELECT
          
          WRITE(unit_nr,'(T2,A,T48,F33.5)') "delocalize_r_cutoff_factor:",main_objectM1%quencher_r0_factor
          
          IF (main_objectM1%almo_update_algorithm_q.eq.almo_scf_pcg) THEN
             WRITE(unit_nr,'(T2,A,T48,E33.3)') "delocalize_eps_iter:",main_objectM1%eps_scf_q
             WRITE(unit_nr,'(T2,A,T48,I33)') "delocalize_max_iter:",main_objectM1%max_scf_q
             WRITE(unit_nr,'(T2,A,T48,E33.3)') "delocalize_eps_lin_search:",main_objectM1%eps_lin_search
             SELECT CASE(main_objectM1%scf_conjugator)
             CASE(cg_zero)
                 WRITE(unit_nr,'(T2,A,T48,A33)') "conjugation:","Steepest Descent"
             CASE(cg_polak_ribiere)
                 WRITE(unit_nr,'(T2,A,T48,A33)') "conjugation:","Polak-Ribiere"
             CASE(cg_fletcher_reeves)
                 WRITE(unit_nr,'(T2,A,T48,A33)') "conjugation:","Fletcher-Reeves"
             CASE(cg_hestenes_stiefel)
                 WRITE(unit_nr,'(T2,A,T48,A33)') "conjugation:","Hestenes-Stiefel"
             CASE(cg_fletcher)
                 WRITE(unit_nr,'(T2,A,T48,A33)') "conjugation:","Fletcher"
             CASE(cg_liu_storey)
                 WRITE(unit_nr,'(T2,A,T48,A33)') "conjugation:","Liu-Storey"
             CASE(cg_dai_yuan)
                 WRITE(unit_nr,'(T2,A,T48,A33)') "conjugation:","Dai-Yuan"
             CASE(cg_hager_zhang)
                 WRITE(unit_nr,'(T2,A,T48,A33)') "conjugation:","Hager-Zhang"
             END SELECT
          ENDIF
          
       ENDIF
 
       WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
       WRITE(unit_nr,'()')

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_level3_aux1_0

! *****************************************************************************
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level2_aux0_0(matrix_s,main_objectM1,error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s
    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level2_aux0_0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    IF (main_objectM1%orthogonal_basis) THEN
       CALL cp_dbcsr_set(main_objectM1%matrix_s(1),0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(main_objectM1%matrix_s(1),1.0_dp,error=error)
       CALL cp_dbcsr_set(main_objectM1%so_b(1),0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(main_objectM1%so_b(1),1.0_dp,error=error)
    ELSE
       CALL almo_levelX_spec5_0(matrix_s,main_objectM1%matrix_s(1),&
                              main_objectM1,.FALSE.,error=error)
       CALL almo_levelX_spec5_0(matrix_s,main_objectM1%so_b(1),&
                              main_objectM1,.TRUE.,error=error)
    ENDIF
    
    CALL cp_dbcsr_filter(main_objectM1%matrix_s(1),main_objectM1%eps_filter,error=error)
    CALL cp_dbcsr_filter(main_objectM1%so_b(1),main_objectM1%eps_filter,error=error)

    IF (main_objectM1%almo_update_algorithm.eq.almo_scf_diag) THEN
       CALL matrix_sqrt_Newton_Schulz(main_objectM1%so_bs(1),&
               main_objectM1%so_bsi(1),&
               main_objectM1%so_b(1),&
               main_objectM1%eps_filter, 3, 1.0E-4_dp, 40, error=error)
    ELSE IF (main_objectM1%almo_update_algorithm.eq.almo_scf_dm_sign) THEN
       CALL invert_Hotelling(main_objectM1%so_bi(1),&
               main_objectM1%so_b(1),&
               threshold=main_objectM1%eps_filter,&
               error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_level2_aux0_0

! *****************************************************************************
!> \par History
!>       2011.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level1_exp2_1(qs_env,main_objectM1,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level1_exp2_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    IF (main_objectM1%almo_update_algorithm.eq.almo_scf_pcg) THEN
       CALL almo_level2_exp3_1(qs_env,main_objectM1,&
               quench_t=main_objectM1%quench_t_blk,&
               matrix_t_in=main_objectM1%enter_b,&
               matrix_t_out=main_objectM1%enter_b,&
               assume_t0_q0x=.FALSE.,perturbation_only=.FALSE.,&
               special_case=0,error=error)
       CALL almo_scf_ortho_blk(main_objectM1,error)
    ELSE
       CALL almo_level2_exp1(qs_env,main_objectM1,error)
    ENDIF
    
    DO ispin=1,main_objectM1%nspins
       CALL cp_dbcsr_copy(main_objectM1%hfh_con(ispin),&
                 main_objectM1%hfh(ispin),error=error)
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE almo_level1_exp2_1

! *****************************************************************************
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level2_exp1(qs_env,main_objectM1,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level2_exp1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iscf, ispin, nspin, &
                                                unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: local_nocc_of_domain
    LOGICAL                                  :: converged, failure, &
                                                prepare_to_exit, should_stop, &
                                                use_object00, &
                                                use_prev_as_guess
    REAL(KIND=dp) :: energy_diff, energy_new, energy_old, error_norm, &
      error_norm_ispin, t1, t2, true_mixing_fraction
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: local_mu
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: matrix_mixing_old_blk
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(object00_type), ALLOCATABLE, &
      DIMENSION(:)                           :: almo_object00

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    use_object00=.TRUE.
    use_prev_as_guess=.FALSE.

    nspin=main_objectM1%nspins
    ALLOCATE(local_mu(main_objectM1%ndomains))
    ALLOCATE(local_nocc_of_domain(main_objectM1%ndomains))

    ALLOCATE(matrix_mixing_old_blk(nspin))
    ALLOCATE(almo_object00(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(matrix_mixing_old_blk(ispin),error=error)
       CALL cp_dbcsr_create(matrix_mixing_old_blk(ispin),&
              template=main_objectM1%hfh_b(ispin),error=error)
       CALL init_object00_gen(object00_env=almo_object00(ispin),&
              sample_err=main_objectM1%hfh_b(ispin),&
              sample_var=main_objectM1%so_b(1),&
              error_type=1,&
              max_length=main_objectM1%counter04_bd,&
              error=error)
    ENDDO

    energy_old=0.0_dp
    iscf=0
    prepare_to_exit=.FALSE.
    true_mixing_fraction=0.0_dp

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'(T2,A,A,A)') REPEAT("-",20), &
          " Optimization of block-diagonal ALMOs ", REPEAT("-",21)
       WRITE(unit_nr,*)
       WRITE(unit_nr,'(T2,A13,A6,A23,A14,A14,A9)') "Method","Iter",&
               "Total Energy","Change","Convergence","Time"
       WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
    ENDIF

    t1 = m_walltime()
    DO 

      iscf=iscf+1

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
      DO ispin=1,nspin
         CALL almo_levelX_spec5_0(matrix_ks(ispin)%matrix,&
                 main_objectM1%hfh(ispin),&
                 main_objectM1,.FALSE.,error=error)
         CALL almo_levelX_spec5_0(matrix_ks(ispin)%matrix,&
                 main_objectM1%hfh_b(ispin),&
                 main_objectM1,.TRUE.,error=error)
         CALL cp_dbcsr_filter(main_objectM1%hfh(ispin),&
                  main_objectM1%eps_filter,error=error)
      ENDDO

      CALL almo_scf_loc_hp_blk(main_objectM1,error=error)

      IF (use_object00) THEN
         DO ispin=1,nspin
            CALL push_object00(object00_env=almo_object00(ispin),&
                    var=main_objectM1%hfh_b(ispin),&
                    err=main_objectM1%gst_b(ispin),&
                    error=error)
         ENDDO
      ENDIF
 
      DO ispin=1,nspin
         CALL cp_dbcsr_norm(main_objectM1%gst_b(ispin),&
                 dbcsr_norm_maxabsnorm,&
                 norm_scalar=error_norm_ispin, error=error)
         IF (ispin.eq.1) error_norm=error_norm_ispin
         IF (ispin.gt.1 .AND. error_norm_ispin.gt.error_norm) &
            error_norm=error_norm_ispin
      ENDDO
      
      IF (error_norm.lt.main_objectM1%eps_prev_guess) THEN
         use_prev_as_guess = .TRUE.
      ELSE
         use_prev_as_guess = .FALSE.
      ENDIF 
      
      converged=.TRUE.
      IF (error_norm.gt.main_objectM1%eps_scf_bd) converged=.FALSE.
      CALL external_control(should_stop,"SCF",&
              start_time=qs_env%start_time,&
              target_time=qs_env%target_time,error=error)
      IF (should_stop .OR. iscf>=main_objectM1%max_scf_bd .OR. converged)  THEN
         prepare_to_exit=.TRUE.
      ENDIF
      
      IF (.NOT.prepare_to_exit) THEN
         
         IF (iscf.ne.1) THEN
            IF (use_object00) THEN
               DO ispin=1,nspin
                  CALL extrapolate_object00(object00_env=almo_object00(ispin),&
                          extr_var=main_objectM1%hfh_b(ispin),&
                          error=error)
               ENDDO
            ELSE 
               true_mixing_fraction=main_objectM1%mixing_fraction
               DO ispin=1,nspin
                  CALL cp_dbcsr_add(main_objectM1%hfh_b(ispin),&
                                    matrix_mixing_old_blk(ispin),& 
                                    true_mixing_fraction,&
                                    1.0_dp-true_mixing_fraction,&
                                    error=error)
               END DO
            ENDIF
         ENDIF
         DO ispin=1,nspin
            CALL cp_dbcsr_copy(matrix_mixing_old_blk(ispin),&
                   main_objectM1%hfh_b(ispin), error=error)
         ENDDO
   
         SELECT CASE (main_objectM1%almo_update_algorithm) 
         CASE (almo_scf_diag) 
   
            CALL almo_scf_get_hp_blk_and_tv_blk(main_objectM1,error) 
   
         CASE (almo_scf_dm_sign)
   
            DO ispin=1,nspin
      
               local_nocc_of_domain(:)=main_objectM1%nocc_of_domain(:,ispin)
               local_mu(:)=main_objectM1%mu_of_domain(:,ispin)
               IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "CVS only: density_matrix_sign has not been updated in SVN"
               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
               main_objectM1%mu_of_domain(:,ispin)=local_mu(:)
      
            ENDDO
      
            CALL almo_scf_get_t_blk(main_objectM1,error)
            CALL almo_scf_ortho_blk(main_objectM1,error=error)
   
         END SELECT

         CALL almo_scf_p_get_t_blk(main_objectM1,&
                 use_sigma_inv_guess=use_prev_as_guess,error=error)
         
         CALL almo_levelX_spec3_0(qs_env,main_objectM1,energy_new,error)
   
      ENDIF 

      energy_diff=energy_new-energy_old
      energy_old=energy_new
      main_objectM1%almo_scf_energy=energy_new

      t2 = m_walltime()
      IF (unit_nr>0) THEN
         WRITE(unit_nr,'(T2,A13,I6,F23.10,E14.5,F14.9,F9.2)') "ALMO SCF DIIS",&
               iscf,&
               energy_new,energy_diff,error_norm, t2-t1
      ENDIF
      t1 = m_walltime()
 
      IF (prepare_to_exit) EXIT

    ENDDO 

    IF (.NOT.converged)  THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "BLOCK-DIAGONAL ALMO SCF not converged! "
    ENDIF

    DO ispin=1,nspin
       CALL cp_dbcsr_release(matrix_mixing_old_blk(ispin),error=error)
       CALL release_object00(object00_env=almo_object00(ispin),error=error)
    ENDDO
    DEALLOCATE(almo_object00)
    DEALLOCATE(matrix_mixing_old_blk)
    DEALLOCATE(local_mu)
    DEALLOCATE(local_nocc_of_domain)
 
    CALL timestop(handle)

  END SUBROUTINE almo_level2_exp1

! *****************************************************************************
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level2_exp2_4(qs_env,main_objectM1,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level2_exp2_4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iscf, ispin, nspin, &
                                                unit_nr
    LOGICAL                                  :: converged, prepare_to_exit, &
                                                should_stop
    REAL(KIND=dp) :: denergy_tot, energy_diff, energy_new, energy_old, &
      error_norm, error_norm_0, spin_factor, t1, t2
    REAL(KIND=dp), DIMENSION(2)              :: denergy_spin
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type)                      :: matrix_p_almo_scf_converged
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(object00_type), ALLOCATABLE, &
      DIMENSION(:)                           :: almo_object00
    TYPE(object01_type), ALLOCATABLE, &
      DIMENSION(:, :)                        :: obj01_mixing_old_blk

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=main_objectM1%nspins
    IF (nspin == 1) THEN
     spin_factor = 2.0_dp
    ELSE
     spin_factor = 1.0_dp
    ENDIF
   
    ispin=1 
    CALL get_group_sqrt(&
            matrix_s=main_objectM1%matrix_s(1),&
            subm_s_sqrt=main_objectM1%concept_ss(:,ispin),&
            subm_s_sqrt_inv=main_objectM1%concept_ss_inv(:,ispin),&
            dpattern=main_objectM1%quench_t(ispin),&
            map=main_objectM1%domain_map(ispin),&
            node_of_domain=main_objectM1%cpu_of_domain,&
            error=error)
    DO ispin=1,nspin
       CALL construct_object01(&
               matrix=main_objectM1%quench_t(ispin),&
               object01=main_objectM1%concept_enter(:,ispin),&
               distr_pattern=main_objectM1%quench_t(ispin),&
               domain_map=main_objectM1%domain_map(ispin),&
               node_of_domain=main_objectM1%cpu_of_domain,&
               job_type=select_row,error=error)
    ENDDO

    ALLOCATE(obj01_mixing_old_blk(main_objectM1%ndomains,nspin))
    CALL init_object01_gen(obj01_mixing_old_blk,error=error)
    ALLOCATE(almo_object00(nspin))
    

    DO ispin=1,nspin
       CALL init_object00_gen(object00_env=almo_object00(ispin),&
              sample_err=main_objectM1%concept_ss(:,ispin),&
              sample_var=main_objectM1%concept_ss_inv(:,ispin),&
              error_type=1,&
              max_length=main_objectM1%counter04_bd,&
              error=error)
    ENDDO

    denergy_tot=0.0_dp
    energy_old=0.0_dp
    iscf=0
    prepare_to_exit=.FALSE.

    t1 = m_walltime()
    DO 

      iscf=iscf+1

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
      DO ispin=1,nspin
         CALL almo_levelX_spec5_0(matrix_ks(ispin)%matrix,&
                 main_objectM1%hfh(ispin),&
                 main_objectM1,.FALSE.,error=error)
         CALL almo_levelX_spec5_0(matrix_ks(ispin)%matrix,&
                 main_objectM1%hfh_b(ispin),&
                 main_objectM1,.TRUE.,error=error)
         CALL cp_dbcsr_filter(main_objectM1%hfh(ispin),&
                  main_objectM1%eps_filter,error=error)
      ENDDO

      CALL almo_scf_loc_hp_xx(main_objectM1,error=error)
      
      DO ispin=1,nspin
         CALL push_object00(object00_env=almo_object00(ispin),&
                 d_var=main_objectM1%concept_hfh_xx(:,ispin),&
                 d_err=main_objectM1%concept_gst(:,ispin),&
                 error=error)
      ENDDO
 
      converged=.TRUE.
      DO ispin=1,nspin
         CALL cp_dbcsr_norm(main_objectM1%gst_xx(ispin),&
                 dbcsr_norm_maxabsnorm,&
                 norm_scalar=error_norm, error=error)
         CALL op3_object01(main_objectM1%concept_gst(:,ispin),&
                 norm=error_norm_0,&
                 error=error)
         IF (error_norm.gt.main_objectM1%eps_scf_q) THEN
            converged=.FALSE.
            EXIT 
         ENDIF
      ENDDO
      CALL external_control(should_stop,"SCF",&
              start_time=qs_env%start_time,&
              target_time=qs_env%target_time,error=error)
      IF (should_stop .OR. iscf>=main_objectM1%max_scf_q .OR. converged)  THEN
         prepare_to_exit=.TRUE.
      ENDIF

      IF (.NOT.prepare_to_exit) THEN
         
         IF (iscf.ne.1) THEN
            IF (.FALSE.) THEN 
               DO ispin=1,nspin
                  CALL op2_object01_gen(&
                          main_objectM1%mixing_fraction,&
                          main_objectM1%concept_hfh_xx(:,ispin),&
                          1.0_dp-main_objectM1%mixing_fraction,&
                          obj01_mixing_old_blk(:,ispin),&
                          'N',error)
               END DO
            ELSE
               DO ispin=1,nspin
                  CALL extrapolate_object00(object00_env=almo_object00(ispin),&
                          d_extr_var=main_objectM1%concept_hfh_xx(:,ispin),&
                          error=error)
               ENDDO
            ENDIF
         ENDIF
         DO ispin=1,nspin
            CALL copy_object01_gen(&
                   main_objectM1%concept_hfh_xx(:,ispin),&
                   obj01_mixing_old_blk(:,ispin),&
                   copy_data=.TRUE., error=error)
         ENDDO
   
         CALL almo_scf_get_hp_xx_and_tv_xx(main_objectM1,error) 

         DO ispin=1,nspin
            
            IF (iscf.eq.1) THEN
               CALL cp_dbcsr_init(matrix_p_almo_scf_converged, error=error)
               CALL cp_dbcsr_create(matrix_p_almo_scf_converged,&
                       template=main_objectM1%dpp(ispin),&
                       error=error)
               CALL cp_dbcsr_copy(matrix_p_almo_scf_converged,&
                       main_objectM1%dpp(ispin),error=error)
            ENDIF

            CALL almo_scf_p_get_t(&
                    t=main_objectM1%enter(ispin),&
                    p=main_objectM1%dpp(ispin),&
                    eps_filter=main_objectM1%eps_filter,&
                    orthog_orbs=.FALSE.,&
                    s=main_objectM1%matrix_s(1),&
                    sigma=main_objectM1%rem(ispin),&
                    sigma_inv=main_objectM1%rem_i(ispin),&
                    use_guess=.TRUE.,&
                    error=error)
            CALL cp_dbcsr_scale(main_objectM1%dpp(ispin),spin_factor,&
                    error=error)

            IF (iscf.eq.1) THEN
              
               CALL cp_dbcsr_add(matrix_p_almo_scf_converged,&
                                main_objectM1%dpp(ispin),-1.0_dp,1.0_dp,&
                                error=error)
               CALL cp_dbcsr_trace(main_objectM1%hfh_con(ispin),&
                                   matrix_p_almo_scf_converged,&
                                   denergy_spin(ispin),error=error)
            
               CALL cp_dbcsr_release(matrix_p_almo_scf_converged, error=error)
               
               denergy_tot=denergy_tot+denergy_spin(ispin)
  
            ENDIF

         ENDDO
         
         IF (iscf.eq.1) THEN
            IF (unit_nr>0) THEN
               WRITE(unit_nr,*)
               WRITE(unit_nr,'(T2,A35,F25.10)') "ENERGY OF BLOCK-DIAGONAL ALMOs:",&
                  main_objectM1%almo_scf_energy
               WRITE(unit_nr,'(T2,A35,F25.10)') "ENERGY LOWERING:",&
                  denergy_tot
               WRITE(unit_nr,'(T2,A35,F25.10)') "CORRECTED ENERGY:",&
                  main_objectM1%almo_scf_energy+denergy_tot
               WRITE(unit_nr,*)
            ENDIF
            CALL almo_levelX_spec2_0(qs_env,&
                    main_objectM1%almo_scf_energy+denergy_tot,&
                    error)
         ENDIF

         IF (.NOT.main_objectM1%perturbative_delocalization) THEN
            CALL almo_levelX_spec3_0(qs_env,main_objectM1,energy_new,error)
         ENDIF
   
      ENDIF 

      IF (main_objectM1%perturbative_delocalization) THEN
         
         converged=.TRUE.
         prepare_to_exit=.TRUE.

      ELSE

         energy_diff=energy_new-energy_old
         energy_old=energy_new
         main_objectM1%almo_scf_energy=energy_new

         t2 = m_walltime()
         IF (unit_nr>0) THEN
            WRITE(unit_nr,'(T2,A,I6,F20.9,E11.3,E11.3,E11.3,F8.2)') "ALMO SCF",&
                  iscf,&
                  energy_new,energy_diff,error_norm,error_norm_0, t2-t1
         ENDIF
         t1 = m_walltime()

      ENDIF
 
      IF (prepare_to_exit) EXIT

    ENDDO 

    IF (.NOT.converged)  THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "ALMO SCF quenched eigensolver not converged! "
    ENDIF

    DO ispin=1,nspin
       CALL release_object01_gen(obj01_mixing_old_blk(:,ispin),error=error)
       CALL release_object00(object00_env=almo_object00(ispin),error=error)
    ENDDO
    DEALLOCATE(almo_object00)
    DEALLOCATE(obj01_mixing_old_blk)
 
    CALL timestop(handle)

  END SUBROUTINE almo_level2_exp2_4

! *****************************************************************************
!> \brief All optimization procedures described in JCTC, (2013)
!> \par History
!>       2011.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level2_exp3_1(qs_env,main_objectM1,quench_t,&
                matrix_t_in,matrix_t_out,assume_t0_q0x,perturbation_only,&
                special_case,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: quench_t, matrix_t_in, &
                                                matrix_t_out
    LOGICAL, INTENT(IN)                      :: assume_t0_q0x, &
                                                perturbation_only
    INTEGER, INTENT(IN), OPTIONAL            :: special_case
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level2_exp3_1', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: iunit = 78

    CHARACTER(LEN=20)                        :: iter_type
    INTEGER :: cg_iteration, fixed_line_search_niter, handle, ispin, &
      iteration, line_search_iteration, max_iter, my_special_case, ndomains, &
      outer_iteration, outer_max_iter, prec_type, precond_domain_projector, &
      unit_nr
    LOGICAL :: converged, do_md, failure, first_md_iteration, just_started, &
      line_search, md_in_theta_space, optimize_theta, outer_prepare_to_exit, &
      prepare_to_exit, reset_conjugator, skip_grad, use_guess, &
      use_preconditioner
    REAL(kind=dp) :: appr_sec_der, beta, denom, e0, e1, energy_diff, &
      energy_new, energy_old, eps_skip_gradients, g0, g1, grad_norm, &
      grad_norm_frob, kappa, kin_energy, line_search_error, &
      next_step_size_guess, numer, prec_sf_mixing_s, spin_factor, step_size, &
      t1, t2, t_norm, tau, time_step
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type) :: FTsiginv, fvo_0, grad, m_theta, m_tmp_nn_1, &
      m_tmp_no_1, m_tmp_no_2, m_tmp_no_3, m_tmp_oo_1, matrix_p_0, &
      matrix_sigma_0, matrix_sigma_inv_0, matrix_t_0, prec_oo, prec_oo_inv, &
      prec_vv, prev_grad, prev_minus_prec_grad, prev_step, siginvTFTsiginv, &
      ST, step, STsiginv_0, velocity
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(object01_type), ALLOCATABLE, &
      DIMENSION(:)                           :: domain_r_down

    CALL timeset(routineN,handle)

    my_special_case=-1
    IF (PRESENT(special_case)) my_special_case=special_case

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    IF (unit_nr>0) THEN
       WRITE(unit_nr,*)
       WRITE(unit_nr,'(T2,A,A,A)') REPEAT("-",28), &
          " Optimization of ALMOs ", REPEAT("-",28)
       WRITE(unit_nr,*)
       WRITE(unit_nr,'(T2,A13,A6,A23,A14,A14,A9)') "Method","Iter",&
               "Objective Function","Change","Convergence","Time"
       WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
    ENDIF

    do_md=main_objectM1%logical01
    use_preconditioner=main_objectM1%logical04
    optimize_theta=main_objectM1%logical05
    prec_sf_mixing_s=main_objectM1%real04
    prec_type=main_objectM1%integer02
    eps_skip_gradients=main_objectM1%real01
    fixed_line_search_niter=0 
    
    IF (main_objectM1%nspins == 1) THEN
       spin_factor = 2.0_dp
    ELSE
       spin_factor = 1.0_dp
    ENDIF
    
IF (main_objectM1%nspins.gt.1) THEN
   CPErrorMessage(cp_failure_level,routineP,"MODIFY CODE FOR UNRESTRICTED",error)
   CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
ENDIF

    DO ispin=1,main_objectM1%nspins

       CALL cp_dbcsr_init(m_theta,error=error)
       CALL cp_dbcsr_init(prec_vv,error=error)
       CALL cp_dbcsr_init(fvo_0,error=error)
       CALL cp_dbcsr_init(STsiginv_0,error=error)
       CALL cp_dbcsr_init(m_tmp_no_1,error=error)
       CALL cp_dbcsr_init(m_tmp_no_2,error=error)
       CALL cp_dbcsr_init(m_tmp_no_3,error=error)
       CALL cp_dbcsr_init(ST,error=error)
       CALL cp_dbcsr_init(FTsiginv,error=error)
       CALL cp_dbcsr_init(m_tmp_oo_1,error=error)
       CALL cp_dbcsr_init(m_tmp_nn_1,error=error)
       CALL cp_dbcsr_init(siginvTFTsiginv,error=error)
       CALL cp_dbcsr_init(prec_oo,error=error)
       CALL cp_dbcsr_init(prec_oo_inv,error=error)
       CALL cp_dbcsr_init(prev_grad,error=error)
       CALL cp_dbcsr_init(prev_step,error=error)
       CALL cp_dbcsr_init(grad,error=error)
       CALL cp_dbcsr_init(step,error=error)
       CALL cp_dbcsr_init(prev_minus_prec_grad,error=error)
       CALL cp_dbcsr_create(prec_vv,&
               template=main_objectM1%hfh(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(prec_oo,&
               template=main_objectM1%rem(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(prec_oo_inv,&
               template=main_objectM1%rem(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(m_tmp_oo_1,&
               template=main_objectM1%rem(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(siginvTFTsiginv,&
               template=main_objectM1%rem(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(STsiginv_0,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(fvo_0,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(m_tmp_no_1,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(m_tmp_no_2,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(m_tmp_no_3,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(FTsiginv,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(ST,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(m_theta,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(prev_grad,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(grad,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(prev_step,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(step,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(prev_minus_prec_grad,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       
       ndomains = cp_dbcsr_nblkcols_total(quench_t(ispin))
       ALLOCATE(domain_r_down(ndomains))
       CALL init_object01_gen(domain_r_down,error)

       CALL cp_dbcsr_init(matrix_t_0,error=error)
       CALL cp_dbcsr_init(matrix_sigma_inv_0,error=error)
       CALL cp_dbcsr_init(matrix_sigma_0,error=error)
       CALL cp_dbcsr_init(matrix_p_0,error=error)
       CALL cp_dbcsr_create(matrix_t_0,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(matrix_sigma_inv_0,&
               template=main_objectM1%rem_i(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(matrix_sigma_0,&
               template=main_objectM1%rem_i(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(matrix_p_0,&
               template=main_objectM1%hfh(ispin),&
               error=error)
       CALL cp_dbcsr_copy(matrix_t_0,matrix_t_in(ispin),&
               error=error)

       CALL cp_dbcsr_set(step,0.0_dp,error=error)

       md_in_theta_space=.FALSE. 
       IF (do_md) THEN
          CALL cp_dbcsr_init(velocity,error=error)
          CALL cp_dbcsr_create(velocity,&
                  template=matrix_t_out(ispin),error=error)
          CALL cp_dbcsr_copy(velocity,quench_t(ispin),error=error)
          CALL cp_dbcsr_set(velocity,0.0_dp,error=error)
          CALL cp_dbcsr_copy(prev_step,quench_t(ispin),error=error)
          CALL cp_dbcsr_set(prev_step,0.0_dp,error=error)
          time_step=main_objectM1%real05
       ENDIF
       
       IF (assume_t0_q0x) THEN
          CALL cp_dbcsr_set(m_theta,0.0_dp,error=error)
       ELSE
          IF (optimize_theta) THEN
             CALL cp_dbcsr_norm(matrix_t_0,&
                     dbcsr_norm_maxabsnorm, norm_scalar=grad_norm, error=error)
             IF (grad_norm.gt.main_objectM1%envelope_amplitude) THEN
                CPErrorMessage(cp_failure_level,routineP,"Max norm of the initial guess is too large",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF
             CALL cp_dbcsr_copy(m_theta,matrix_t_0,error=error)
             CALL cp_dbcsr_function_of_elements(m_theta,&
                     !func=dbcsr_func_asin,&
                     func=dbcsr_func_artanh,&
                     a0=0.0_dp,&
                     a1=1.0_dp/main_objectM1%envelope_amplitude,&
                     error=error)
             CALL cp_dbcsr_scale(m_theta,main_objectM1%envelope_amplitude,&
                     error=error)
          ELSE
             CALL cp_dbcsr_copy(m_theta,matrix_t_0,error=error)
             CALL cp_dbcsr_norm(m_theta,&
                     dbcsr_norm_maxabsnorm, norm_scalar=grad_norm, error=error)
          ENDIF
       ENDIF

       IF (my_special_case.eq.-1) THEN
          CALL get_group_inv(&
                 matrix_s=main_objectM1%matrix_s(1),&
                 subm_s_inv=main_objectM1%concept_soi(:,ispin),&
                 dpattern=quench_t(ispin),&
                 map=main_objectM1%domain_map(ispin),&
                 node_of_domain=main_objectM1%cpu_of_domain,&
                 error=error)
       ENDIF

       outer_max_iter=main_objectM1%outer_max_scf_q
       outer_prepare_to_exit=.FALSE.
       outer_iteration=0
       grad_norm=0.0_dp
       grad_norm_frob=0.0_dp
       use_guess=.FALSE.

       DO

          max_iter=main_objectM1%max_scf_q
          prepare_to_exit=.FALSE.
          line_search=.FALSE.
          converged=.FALSE.
          iteration=0
          cg_iteration=0
          line_search_iteration=0
          energy_new=0.0_dp
          energy_old=0.0_dp
          line_search_error=0.0_dp
          t1 = m_walltime()

          DO

             just_started=(iteration.eq.0).AND.(outer_iteration.eq.0)

             IF (iteration.eq.main_objectM1%integer01.AND.do_md) THEN
                CALL cp_dbcsr_set(velocity,0.0_dp,error=error)
                CALL cp_dbcsr_set(prev_step,0.0_dp,error=error)
                md_in_theta_space=.TRUE.
                first_md_iteration=.TRUE.
             ENDIF

             IF (assume_t0_q0x.AND.just_started) THEN
                   CALL cp_dbcsr_set(matrix_t_out(ispin),0.0_dp,error=error)
             ENDIF
             IF (optimize_theta) THEN
                CALL cp_dbcsr_norm(m_theta,&
                        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
             ENDIF
             IF (optimize_theta) THEN
                CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
                CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
                        !func=dbcsr_func_sin,&
                        func=dbcsr_func_tanh,&
                        a0=0.0_dp,&
                        a1=1.0_dp/main_objectM1%envelope_amplitude,&
                        error=error)
                CALL cp_dbcsr_scale(m_tmp_no_1,&
                        main_objectM1%envelope_amplitude,error=error)
             ELSE
                CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
             ENDIF
             CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
                     quench_t(ispin),&
                     matrix_t_out(ispin),error=error)
             CALL cp_dbcsr_norm(matrix_t_out(ispin),&
                     dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)

             IF (assume_t0_q0x.AND.(.NOT.just_started)) THEN
                IF (my_special_case.eq.1) THEN
                   CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                           STsiginv_0,&
                           matrix_t_out(ispin),&
                           0.0_dp,m_tmp_oo_1,&
                           filter_eps=main_objectM1%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                           matrix_t_0,&
                           m_tmp_oo_1,&
                           1.0_dp,matrix_t_out(ispin),&
                           filter_eps=main_objectM1%eps_filter,&
                           error=error)
                ELSE IF (my_special_case.eq.0) THEN
                   CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                ELSE
                   IF (prec_type.eq.4) THEN
                      CALL operations_on_sets(&
                              matrix_in=matrix_t_out(ispin),&
                              matrix_out=m_tmp_no_1,&
                              operator1=domain_r_down(:),&
                              operator2=main_objectM1%concept_soi(:,ispin),&
                              dpattern=quench_t(ispin),&
                              map=main_objectM1%domain_map(ispin),&
                              node_of_domain=main_objectM1%cpu_of_domain,&
                              my_action=1,&
                              filter_eps=main_objectM1%eps_filter,&
                              use_trimmer=.FALSE.,&
                              error=error)
                   ELSE
                      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   ENDIF 
                   CALL cp_dbcsr_copy(matrix_t_out(ispin),&
                           m_tmp_no_1,error=error)
                ENDIF 
                CALL cp_dbcsr_norm(matrix_t_out(ispin),&
                        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
             ENDIF 

             IF (assume_t0_q0x) THEN
                CALL cp_dbcsr_add(matrix_t_out(ispin),&
                        matrix_t_0,1.0_dp,1.0_dp,&
                        error=error)
                CALL cp_dbcsr_norm(matrix_t_out(ispin),&
                        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
             ENDIF

             CALL cp_dbcsr_filter(matrix_t_out(ispin),&
                     eps=main_objectM1%eps_filter,&
                     error=error)

             CALL almo_scf_p_get_t(&
                     t=matrix_t_out(ispin),&
                     p=main_objectM1%dpp(ispin),&
                     eps_filter=main_objectM1%eps_filter,&
                     orthog_orbs=.FALSE.,&
                     s=main_objectM1%matrix_s(1),&
                     sigma=main_objectM1%rem(ispin),&
                     sigma_inv=main_objectM1%rem_i(ispin),&
                     use_guess=use_guess,&
                     error=error)
             CALL cp_dbcsr_scale(main_objectM1%dpp(ispin),&
                     spin_factor,error=error)

             IF ( .NOT.(perturbation_only.AND.(.NOT.just_started)) ) THEN
                CALL almo_levelX_spec3_0(qs_env,main_objectM1,energy_new,error)
                CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
                CALL almo_levelX_spec5_0(matrix_ks(ispin)%matrix,&
                        main_objectM1%hfh(ispin),&
                        main_objectM1,.FALSE.,error=error)
                CALL cp_dbcsr_filter(main_objectM1%hfh(ispin),&
                        main_objectM1%eps_filter,error=error)
             ENDIF

             IF (just_started) THEN
                CALL cp_dbcsr_copy(matrix_p_0,main_objectM1%dpp(ispin),&
                        error=error)
                CALL cp_dbcsr_desymmetrize(main_objectM1%rem(ispin),&
                        matrix_sigma_0,&
                        error=error)
                CALL cp_dbcsr_desymmetrize(main_objectM1%rem_i(ispin),&
                        matrix_sigma_inv_0,&
                        error=error)
             ENDIF

IF (my_special_case.eq.-1.AND.prec_type.eq.4) THEN 
             IF (assume_t0_q0x.AND.just_started) THEN
                CALL get_group_rdown(&
                        matrix_t=matrix_t_0,&
                        matrix_sigma_inv=matrix_sigma_inv_0,&
                        matrix_s=main_objectM1%matrix_s(1),&
                        subm_r_down=domain_r_down(:),&
                        dpattern=quench_t(ispin),&
                        map=main_objectM1%domain_map(ispin),&
                        node_of_domain=main_objectM1%cpu_of_domain,&
                        filter_eps=main_objectM1%eps_filter,&
                        error=error)
             ENDIF 
ENDIF

             IF (perturbation_only) THEN
                CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                        main_objectM1%hfh(ispin),&
                        matrix_t_out(ispin),&
                        0.0_dp,m_tmp_no_1,&
                        filter_eps=main_objectM1%eps_filter,&
                        error=error)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                        m_tmp_no_1,&
                        main_objectM1%rem_i(ispin),&
                        0.0_dp,FTsiginv,&
                        filter_eps=main_objectM1%eps_filter,&
                        error=error)
                CALL cp_dbcsr_trace(matrix_t_out(ispin),&
                        FTsiginv,energy_new,"T","N",&
                        error=error)
                energy_new=energy_new*spin_factor
             ENDIF

             IF (line_search_iteration.eq.0.AND.iteration.ne.0) &
                CALL cp_dbcsr_copy(prev_grad,grad,error=error)
             
             skip_grad = ( iteration.gt.0 .AND. &
                           fixed_line_search_niter.ne.0 .AND. &
                           line_search_iteration.ne.fixed_line_search_niter )
                           
             IF (.NOT.skip_grad) THEN
                
                IF (.NOT.perturbation_only) THEN
                   CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           main_objectM1%hfh(ispin),&
                           matrix_t_out(ispin),&
                           0.0_dp,m_tmp_no_1,&
                           filter_eps=main_objectM1%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_multiply("N","N",spin_factor,&
                           m_tmp_no_1,&
                           main_objectM1%rem_i(ispin),&
                           0.0_dp,FTsiginv,&
                           filter_eps=main_objectM1%eps_filter,&
                           error=error)
                ENDIF
                CALL cp_dbcsr_copy(m_tmp_no_2,quench_t(ispin),&
                        error=error)
                CALL cp_dbcsr_copy(m_tmp_no_2,&
                        FTsiginv,keep_sparsity=.TRUE.,error=error)
                CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                        matrix_t_out(ispin),&
                        FTsiginv,&
                        0.0_dp,m_tmp_oo_1,&
                        filter_eps=main_objectM1%eps_filter,&
                        error=error)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                        main_objectM1%rem_i(ispin),&
                        m_tmp_oo_1,&
                        0.0_dp,siginvTFTsiginv,&
                        filter_eps=main_objectM1%eps_filter,&
                        error=error)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                        main_objectM1%matrix_s(1),&
                        matrix_t_out(ispin),&
                        0.0_dp,ST,&
                        filter_eps=main_objectM1%eps_filter,&
                        error=error)
                IF (assume_t0_q0x.AND.just_started&
                   .AND.special_case.eq.1) THEN
                   CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           ST,&
                           main_objectM1%rem_i(ispin),&
                           0.0_dp,STsiginv_0,&
                           filter_eps=main_objectM1%eps_filter,&
                           error=error)
                ENDIF
                CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                        ST,&
                        siginvTFTsiginv,&
                        1.0_dp,m_tmp_no_2,&
                        retain_sparsity=.TRUE.,&
                        error=error)
                CALL cp_dbcsr_scale(m_tmp_no_2,&
                        2.0_dp*spin_factor,&
                        error=error)
                CALL cp_dbcsr_filter(m_tmp_no_2,&
                        eps=main_objectM1%eps_filter,&
                        error=error)

                IF (perturbation_only.AND.just_started) THEN
                   CALL cp_dbcsr_copy(fvo_0,m_tmp_no_2,error=error)
                   CALL cp_dbcsr_scale(fvo_0,&
                           0.5_dp,error=error)
                ENDIF

                CALL cp_dbcsr_norm(m_tmp_no_2,&
                        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)

                IF (assume_t0_q0x) THEN
                   IF (my_special_case.eq.1) THEN
                      CALL cp_dbcsr_copy(grad,m_tmp_no_2,error=error)
                      CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                              matrix_t_0,&
                              grad,&
                              0.0_dp,m_tmp_oo_1,&
                              filter_eps=main_objectM1%eps_filter,&
                              error=error)
                      CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                              STsiginv_0,&
                              m_tmp_oo_1,&
                              1.0_dp,grad,&
                              filter_eps=main_objectM1%eps_filter,&
                              error=error)
                   ELSE IF (my_special_case.eq.0) THEN
                      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   ELSE
                      IF (prec_type.eq.4) THEN
                         CALL operations_on_sets(&
                                 matrix_in=m_tmp_no_2,&
                                 matrix_out=grad,&
                                 operator2=domain_r_down(:),&
                                 operator1=main_objectM1%concept_soi(:,ispin),&
                                 dpattern=quench_t(ispin),&
                                 map=main_objectM1%domain_map(ispin),&
                                 node_of_domain=main_objectM1%cpu_of_domain,&
                                 my_action=1,&
                                 filter_eps=main_objectM1%eps_filter,&
                                 use_trimmer=.FALSE.,&
                                 error=error)
                      ELSE
                         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                      ENDIF
                   ENDIF 
                   CALL cp_dbcsr_copy(m_tmp_no_2,grad,error=error)
                ENDIF

                IF (optimize_theta) THEN
                   CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
                   CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
                           func=dbcsr_func_dtanh,&
                           a0=0.0_dp,&
                           a1=1.0_dp/main_objectM1%envelope_amplitude,&
                           error=error)
                   CALL cp_dbcsr_scale(m_tmp_no_1,&
                           main_objectM1%envelope_amplitude,&
                           error=error)
                   CALL cp_dbcsr_set(m_tmp_no_3,0.0_dp,error=error)
                   CALL cp_dbcsr_filter(m_tmp_no_3,&
                           eps=main_objectM1%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_hadamard_product(m_tmp_no_2,&
                           m_tmp_no_1,&
                           m_tmp_no_3,&
                           b_assume_value=1.0_dp,&
                           error=error)
                   CALL cp_dbcsr_hadamard_product(m_tmp_no_3,&
                           quench_t(ispin),&
                           grad,&
                           error=error)
                ELSE 
                   CALL cp_dbcsr_hadamard_product(m_tmp_no_2,&
                           quench_t(ispin),&
                           grad,&
                           error=error)
                ENDIF
                CALL cp_dbcsr_filter(grad,eps=main_objectM1%eps_filter,&
                        error=error)

             ENDIF
             
             grad_norm_frob=cp_dbcsr_frobenius_norm(grad)
             CALL cp_dbcsr_norm(grad, dbcsr_norm_maxabsnorm,&
                     norm_scalar=grad_norm, error=error)
             converged=(grad_norm.lt.main_objectM1%eps_scf_q)
             IF (converged.OR.(iteration.ge.max_iter)) THEN
                prepare_to_exit=.TRUE.
             ENDIF
             IF (grad_norm.lt.main_objectM1%eps_prev_guess) THEN
                use_guess=.TRUE.
             ENDIF

             IF (md_in_theta_space) THEN

                IF (.NOT.first_md_iteration) THEN
                   CALL cp_dbcsr_copy(prev_step,step,error=error)
                ENDIF
                CALL cp_dbcsr_copy(step,grad,error=error)
                CALL cp_dbcsr_scale(step,-1.0_dp,error=error)
   
                IF (.NOT.first_md_iteration) THEN
                   CALL cp_dbcsr_add(velocity,&
                           step,1.0_dp,0.5_dp*time_step,error=error)
                   CALL cp_dbcsr_add(velocity,&
                           prev_step,1.0_dp,0.5_dp*time_step,error=error)
                ENDIF
                kin_energy=cp_dbcsr_frobenius_norm(velocity)
                kin_energy=0.5_dp*kin_energy*kin_energy
   
                CALL cp_dbcsr_add(m_theta,&
                        velocity,1.0_dp,time_step,error=error)
                CALL cp_dbcsr_add(m_theta,&
                        step,1.0_dp,0.5_dp*time_step*time_step,error=error)
                
                iter_type="MD"

                t2 = m_walltime()
                IF (unit_nr>0) THEN
                   WRITE(unit_nr,'(T2,A,A2,I5,F16.7,F17.9,F17.9,F17.9,E12.3,F10.3)') &
                           "ALMO SCF ",iter_type,iteration,time_step*iteration,&
                           energy_new,kin_energy,energy_new+kin_energy,grad_norm,&
                           t2-t1
                ENDIF
                t1 = m_walltime()
                
                IF (first_md_iteration) THEN
                   first_md_iteration=.FALSE.
                ENDIF

             ELSE

                IF (.NOT.prepare_to_exit) THEN
      
                   IF (iteration.ne.0) THEN

                      IF (fixed_line_search_niter.eq.0) THEN

                         CALL cp_dbcsr_trace(grad,step,line_search_error,&
                                 "T","N",error=error)
                         CALL cp_dbcsr_trace(grad,grad,denom,"T","N",&
                                 error=error)
                         line_search_error=line_search_error/SQRT(denom)
                         CALL cp_dbcsr_trace(step,step,denom,"T","N",&
                                 error=error)
                         line_search_error=line_search_error/SQRT(denom)
                         IF (ABS(line_search_error).gt.main_objectM1%eps_lin_search) THEN
                            line_search=.TRUE.
                            line_search_iteration=line_search_iteration+1
                         ELSE
                            line_search=.FALSE.
                            line_search_iteration=0
                            IF (grad_norm.lt.eps_skip_gradients) THEN
                               fixed_line_search_niter=ABS(main_objectM1%integer04)
                            ENDIF
                         ENDIF
                      
                      ELSE
                         
                         IF (.NOT.line_search) THEN
                            line_search=.TRUE.
                            line_search_iteration=line_search_iteration+1
                         ELSE
                            IF (line_search_iteration.eq.fixed_line_search_niter) THEN
                               line_search=.FALSE.
                               line_search_iteration=0
                            ELSE
                               line_search_iteration=line_search_iteration+1
                            ENDIF
                         ENDIF
                      
                      ENDIF 
                   ENDIF
   
                   IF (line_search) THEN
                         energy_diff=0.0_dp
                   ELSE
                         energy_diff=energy_new-energy_old
                         energy_old=energy_new
                   ENDIF
   
                   IF (.NOT.line_search) THEN
   
                      cg_iteration=cg_iteration+1 
                      
                      IF ( (just_started .AND. perturbation_only) .OR. &
                           (iteration.eq.0 .AND. (.NOT.perturbation_only)) ) THEN
      
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                                 ST,&
                                 main_objectM1%rem_i(ispin),&
                                 0.0_dp,m_tmp_no_3,&
                                 filter_eps=main_objectM1%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_create(m_tmp_nn_1,&
                                 template=main_objectM1%matrix_s(1),&
                                 matrix_type=dbcsr_type_no_symmetry,error=error) 
                         CALL cp_dbcsr_desymmetrize(main_objectM1%matrix_s(1),&
                                 m_tmp_nn_1,error=error)
                         IF (my_special_case.eq.1) THEN
                         ELSE
                            CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
                                    ST,&
                                    m_tmp_no_3,&
                                    1.0_dp,m_tmp_nn_1,&
                                    filter_eps=main_objectM1%eps_filter,&
                                    error=error)
                         ENDIF
                         CALL cp_dbcsr_create(prec_vv,&
                                 template=main_objectM1%hfh(ispin),&
                                 matrix_type=dbcsr_type_no_symmetry,error=error) 
                         CALL cp_dbcsr_desymmetrize(main_objectM1%hfh(ispin),&
                                 prec_vv,error=error)
                         CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
                                 FTsiginv,&
                                 ST,&
                                 1.0_dp,prec_vv,&
                                 filter_eps=main_objectM1%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
                                 ST,&
                                 FTsiginv,&
                                 1.0_dp,prec_vv,&
                                 filter_eps=main_objectM1%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                                 ST,&
                                 siginvTFTsiginv,&
                                 0.0_dp,m_tmp_no_3,&
                                 filter_eps=main_objectM1%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_multiply("N","T",1.0_dp,&
                                 m_tmp_no_3,&
                                 ST,&
                                 1.0_dp,prec_vv,&
                                 filter_eps=main_objectM1%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_add(prec_vv,m_tmp_nn_1,&
                                 1.0_dp-prec_sf_mixing_s,&
                                 prec_sf_mixing_s,error=error)
                         CALL cp_dbcsr_scale(prec_vv,2.0_dp*spin_factor,error=error)
                         CALL cp_dbcsr_copy(m_tmp_nn_1,prec_vv,error=error)
                         
                         IF (my_special_case.eq.0) THEN 
                            CALL invert_blk_once(matrix_in=m_tmp_nn_1,&
                                    matrix_out=prec_vv,error=error)
                         ELSE IF (my_special_case.eq.1) THEN 
                            CALL cp_dbcsr_cholesky_decompose(prec_vv,&
                                    para_env=main_objectM1%para_env,&
                                    blacs_env=main_objectM1%blacs_env,error=error)
                            CALL cp_dbcsr_cholesky_invert(prec_vv,&
                                    para_env=main_objectM1%para_env,&
                                    blacs_env=main_objectM1%blacs_env,&
                                    upper_to_full=.TRUE.,error=error)
                            CALL cp_dbcsr_filter(prec_vv,&
                                    eps=main_objectM1%eps_filter,&
                                    error=error)
                         ELSE
                            IF (assume_t0_q0x) THEN
                               precond_domain_projector=-1
                            ELSE
                               precond_domain_projector=0
                            ENDIF
                            IF (prec_type.eq.4) THEN
                               CALL get_group_complex(&
                                  matrix_main=m_tmp_nn_1,&
                                  subm_s_inv=main_objectM1%concept_soi(:,ispin),&
                                  subm_r_down=domain_r_down(:),&
                                  matrix_trimmer=quench_t(ispin),&
                                  dpattern=quench_t(ispin),&
                                  map=main_objectM1%domain_map(ispin),&
                                  node_of_domain=main_objectM1%cpu_of_domain,&
                                  preconditioner=main_objectM1%concept_p(:,ispin),&
                                  use_trimmer=.FALSE.,&
                                  my_action=precond_domain_projector,&
                                  error=error)
                            ENDIF 
                         ENDIF

                      ENDIF
      
                      CALL cp_dbcsr_copy(prev_step,step,error=error)
      
                      IF (use_preconditioner) THEN
                         
                         SELECT CASE (prec_type)
                         CASE (4)
              
                            IF (my_special_case.eq.0.OR.my_special_case.eq.1) THEN
                               
                               CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                                       prec_vv,&
                                       grad,&
                                       0.0_dp,step,&
                                       filter_eps=main_objectM1%eps_filter,&
                                       error=error)

                            ELSE
                               
                               IF (optimize_theta) THEN
                                  CPErrorMessage(cp_failure_level,routineP,"theta is NYI",error)
                                  CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                               ENDIF
   
                               CALL operations_on_sets(&
                                       matrix_in=grad,&
                                       matrix_out=step,&
                                       operator1=main_objectM1%concept_p(:,ispin),&
                                       dpattern=quench_t(ispin),&
                                       map=main_objectM1%domain_map(ispin),&
                                       node_of_domain=main_objectM1%cpu_of_domain,&
                                       my_action=0,&
                                       filter_eps=main_objectM1%eps_filter,&
                                       error=error)
                               CALL cp_dbcsr_scale(step,-1.0_dp,error=error)
      
                               CALL cp_dbcsr_copy(m_tmp_no_3,&
                                       quench_t(ispin),&
                                       error=error)
                               CALL cp_dbcsr_function_of_elements(m_tmp_no_3,&
                                       func=dbcsr_func_inverse,&
                                       a0=0.0_dp,&
                                       a1=1.0_dp,&
                                       error=error)
                               CALL cp_dbcsr_copy(m_tmp_no_2,step,error=error)
                               CALL cp_dbcsr_hadamard_product(&
                                       m_tmp_no_2,&
                                       m_tmp_no_3,&
                                       step,&
                                       error=error)
                               CALL cp_dbcsr_copy(m_tmp_no_3,quench_t(ispin),error=error)
                                                           
                            ENDIF
               
                         CASE DEFAULT
                            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         END SELECT 

                      ELSE

                         CALL cp_dbcsr_copy(step,grad,error=error)
                         CALL cp_dbcsr_scale(step,-1.0_dp,error=error)

                      ENDIF
   
                      IF (iteration.eq.0) THEN
                         reset_conjugator=.TRUE.
                      ENDIF
   
                      IF (.NOT.reset_conjugator) THEN
     
                         SELECT CASE (main_objectM1%scf_conjugator)
                         CASE (cg_hestenes_stiefel)
                            CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
                            CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,&
                                    1.0_dp,-1.0_dp,error=error)
                            CALL cp_dbcsr_trace(m_tmp_no_1,step,numer,&
                                    "T","N",error=error)
                            CALL cp_dbcsr_trace(m_tmp_no_1,prev_step,denom,&
                                    "T","N",error=error)
                            beta=-1.0_dp*numer/denom
                         CASE (cg_fletcher_reeves)
                            CALL cp_dbcsr_trace(grad,step,numer,"T","N",error=error)
                            CALL cp_dbcsr_trace(prev_grad,prev_minus_prec_grad,denom,"T","N",error=error)
                            beta=numer/denom
                         CASE (cg_polak_ribiere)
                            CALL cp_dbcsr_trace(prev_grad,prev_minus_prec_grad,denom,"T","N",error=error)
                            CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
                            CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,1.0_dp,-1.0_dp,error=error)
                            CALL cp_dbcsr_trace(m_tmp_no_1,step,numer,"T","N",error=error)
                            beta=numer/denom
                         CASE (cg_fletcher)
                            CALL cp_dbcsr_trace(grad,step,numer,"T","N",error=error)
                            CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N",error=error)
                            beta=numer/denom
                         CASE (cg_liu_storey)
                            CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N",error=error)
                            CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
                            CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,1.0_dp,-1.0_dp,error=error)
                            CALL cp_dbcsr_trace(m_tmp_no_1,step,numer,"T","N",error=error)
                            beta=numer/denom
                         CASE (cg_dai_yuan)
                            CALL cp_dbcsr_trace(grad,step,numer,"T","N",error=error)
                            CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
                            CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,1.0_dp,-1.0_dp,error=error)
                            CALL cp_dbcsr_trace(m_tmp_no_1,prev_step,denom,"T","N",error=error)
                            beta=-1.0_dp*numer/denom
                         CASE (cg_hager_zhang)
                            CALL cp_dbcsr_copy(m_tmp_no_1,grad,error=error)
                            CALL cp_dbcsr_add(m_tmp_no_1,prev_grad,1.0_dp,-1.0_dp,error=error)
                            CALL cp_dbcsr_trace(m_tmp_no_1,prev_step,denom,"T","N",error=error)
                            CALL cp_dbcsr_trace(m_tmp_no_1,prev_minus_prec_grad,numer,"T","N",error=error)
                            kappa=-2.0_dp*numer/denom
                            CALL cp_dbcsr_trace(m_tmp_no_1,step,numer,"T","N",error=error)
                            tau=-1.0_dp*numer/denom
                            CALL cp_dbcsr_trace(prev_step,grad,numer,"T","N",error=error)
                            beta=tau-kappa*numer/denom
                         CASE (cg_zero)
                            beta=0.0_dp
                         CASE DEFAULT
                            CPErrorMessage(cp_failure_level,routineP,"illegal conjugator",error)
                            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         END SELECT
      
                         IF (beta.lt.0.0_dp) THEN
                            IF (unit_nr>0) THEN
                               WRITE(unit_nr,*) "Beta is negative, ", beta
                            ENDIF
                            reset_conjugator=.TRUE.
                         ENDIF
      
                      ENDIF
      
                      IF (reset_conjugator) THEN 
      
                         beta=0.0_dp
                         IF (unit_nr>0 .AND. (.NOT.just_started)) THEN
                            WRITE(unit_nr,*) "(Re)-setting conjugator to zero"
                         ENDIF
                         reset_conjugator=.FALSE.
   
                      ENDIF
                   
                      CALL cp_dbcsr_copy(prev_minus_prec_grad,step,error=error)
      
                      CALL cp_dbcsr_add(step,prev_step,1.0_dp,beta,error=error)
   
                   ENDIF
      
                   IF (.NOT.line_search) THEN
                      e0=energy_new
                      CALL cp_dbcsr_trace(grad,step,g0,"T","N",error=error)
                      IF (iteration.eq.0) THEN
                         step_size=main_objectM1%real05 
                      ELSE
                         IF (next_step_size_guess.le.0.0_dp) THEN
                            step_size=main_objectM1%real05
                         ELSE
                            step_size=next_step_size_guess*1.05_dp
                         ENDIF
                      ENDIF
                      next_step_size_guess=step_size
                   ELSE
                      IF (fixed_line_search_niter.eq.0) THEN
                         e1=energy_new
                         CALL cp_dbcsr_trace(grad,step,g1,"T","N",error=error)
                         appr_sec_der=(g1-g0)/step_size
                         step_size=-g1/appr_sec_der
                         e0=e1
                         g0=g1
                      ELSE
                         e1=energy_new
                         appr_sec_der=2.0*( (e1-e0)/step_size - g0 )/step_size
                         g1=appr_sec_der*step_size + g0
                         IF (unit_nr>0) THEN
                            WRITE(unit_nr,'(A2,7F12.5)') &
                                    "EG",e0,e1,g0,g1,appr_sec_der,step_size,-g1/appr_sec_der
                         ENDIF
                         step_size=-g1/appr_sec_der
                         e0=e1
                         g0=g1
                      ENDIF
                      next_step_size_guess=next_step_size_guess+step_size
                   ENDIF
   
                   CALL cp_dbcsr_add(m_theta,step,1.0_dp,step_size,error=error)
   
                ENDIF 
              
                IF (line_search) THEN
                   iter_type="LS"
                ELSE
                   iter_type="CG"
                ENDIF
   
                t2 = m_walltime()
                IF (unit_nr>0) THEN
                   iter_type=TRIM("ALMO SCF "//iter_type)
                   WRITE(unit_nr,'(T2,A13,I6,F23.10,E14.5,F14.9,F9.2)') &
                           iter_type,iteration,&
                           energy_new,energy_diff,grad_norm,&
                           t2-t1
                ENDIF

                t1 = m_walltime()
   
   
             ENDIF 


             iteration=iteration+1
             IF (prepare_to_exit) EXIT

          ENDDO 

          IF (converged.OR.(outer_iteration.ge.outer_max_iter)) THEN
             outer_prepare_to_exit=.TRUE.
          ENDIF
          
          outer_iteration=outer_iteration+1
          IF (outer_prepare_to_exit) EXIT

       ENDDO 
       
       IF (converged)  THEN
             
          IF (perturbation_only) THEN
            
             CALL cp_dbcsr_add(matrix_t_0,matrix_t_out(ispin),&
                     -1.0_dp,1.0_dp,&
                     error=error)

             CALL cp_dbcsr_trace(matrix_t_0,&
                     fvo_0,energy_new,"T","N",&
                     error=error)
             IF (unit_nr>0) THEN
                WRITE(unit_nr,*)
                WRITE(unit_nr,'(T2,A35,F25.10)') "ENERGY OF BLOCK-DIAGONAL ALMOs:",&
                   main_objectM1%almo_scf_energy
                WRITE(unit_nr,'(T2,A35,F25.10)') "ENERGY LOWERING:",&
                   energy_new
                WRITE(unit_nr,'(T2,A35,F25.10)') "CORRECTED ENERGY:",&
                   main_objectM1%almo_scf_energy+energy_new
                WRITE(unit_nr,*)
             ENDIF
             CALL almo_levelX_spec2_0(qs_env,&
                     main_objectM1%almo_scf_energy+energy_new,&
                     error)

          ENDIF

       ENDIF 

       IF (md_in_theta_space) THEN
          CALL cp_dbcsr_release(velocity,error=error)
       ENDIF
       CALL cp_dbcsr_release(m_theta,error=error)
       CALL cp_dbcsr_release(prec_vv,error=error)
       CALL cp_dbcsr_release(prec_oo,error=error)
       CALL cp_dbcsr_release(prec_oo_inv,error=error)
       CALL cp_dbcsr_release(m_tmp_no_1,error=error)
       CALL cp_dbcsr_release(fvo_0,error=error)
       CALL cp_dbcsr_release(STsiginv_0,error=error)
       CALL cp_dbcsr_release(m_tmp_no_2,error=error)
       CALL cp_dbcsr_release(m_tmp_no_3,error=error)
       CALL cp_dbcsr_release(m_tmp_oo_1,error=error)
       CALL cp_dbcsr_release(ST,error=error)
       CALL cp_dbcsr_release(FTsiginv,error=error)
       CALL cp_dbcsr_release(siginvTFTsiginv,error=error)
       CALL cp_dbcsr_release(m_tmp_nn_1,error=error)
       CALL cp_dbcsr_release(prev_grad,error=error)
       CALL cp_dbcsr_release(prev_step,error=error)
       CALL cp_dbcsr_release(grad,error=error)
       CALL cp_dbcsr_release(step,error=error)
       CALL cp_dbcsr_release(prev_minus_prec_grad,error=error)
       CALL cp_dbcsr_release(matrix_p_0, error=error)
       CALL cp_dbcsr_release(matrix_t_0, error=error)
       CALL cp_dbcsr_release(matrix_sigma_0,error=error)
       CALL cp_dbcsr_release(matrix_sigma_inv_0,error=error)
    
       IF (.NOT.converged)  THEN
          IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "Optimization not converged! "
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
    
       DEALLOCATE(domain_r_down)

    ENDDO 
    
    CALL timestop(handle)

  END SUBROUTINE almo_level2_exp3_1

! *****************************************************************************
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level1_exp3_2(qs_env,main_objectM1,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level1_exp3_2', &
      routineP = moduleN//':'//routineN

    INTEGER :: col, handle, hold, iblock_col, iblock_row, ispin, mynode, &
      nblkcols_tot, nblkrows_tot, row, unit_nr
    LOGICAL                                  :: failure, tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: no_quench
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    SELECT CASE (main_objectM1%deloc_method)
    CASE (almo_deloc_x,almo_deloc_full_scf,almo_deloc_x_then_full_scf)
       ALLOCATE(no_quench(main_objectM1%nspins))
       CALL cp_dbcsr_init(no_quench(1),error=error)
       CALL cp_dbcsr_create(no_quench(1),&
               template=main_objectM1%enter(1),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(no_quench(1))))
       CALL cp_dbcsr_work_create(no_quench(1),&
               work_mutable=.TRUE., error=error)
       nblkrows_tot = cp_dbcsr_nblkrows_total(no_quench(1))
       nblkcols_tot = cp_dbcsr_nblkcols_total(no_quench(1))
       DO row = 1, nblkrows_tot
          DO col = 1, nblkcols_tot
             tr = .FALSE.
             iblock_row = row
             iblock_col = col
             CALL cp_dbcsr_get_stored_coordinates(no_quench(1),&
                     iblock_row, iblock_col, tr, hold)
             IF (hold.EQ.mynode) THEN
                NULLIFY (p_new_block)
                CALL cp_dbcsr_reserve_block2d(no_quench(1),&
                        iblock_row, iblock_col, p_new_block)
                CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
                p_new_block(:,:) = 1.0_dp
             ENDIF
          ENDDO
       ENDDO
       CALL cp_dbcsr_finalize(no_quench(1),error=error)
       IF (main_objectM1%nspins.gt.1) THEN
          DO ispin=2,main_objectM1%nspins
             CALL cp_dbcsr_init(no_quench(ispin),error=error)
             CALL cp_dbcsr_create(no_quench(ispin),&
                     template=main_objectM1%enter(1),&
                     matrix_type=dbcsr_type_no_symmetry,error=error)
             CALL cp_dbcsr_copy(no_quench(ispin),no_quench(1),error=error)
          ENDDO
       ENDIF
    END SELECT

    SELECT CASE (main_objectM1%deloc_method)
    CASE (almo_deloc_none,almo_deloc_full_scf)
       DO ispin=1,main_objectM1%nspins
          CALL cp_dbcsr_copy(main_objectM1%enter(ispin),&
                  main_objectM1%enter_b(ispin),error=error)
       ENDDO
    CASE (almo_deloc_x,almo_deloc_xk,almo_deloc_x_then_full_scf)
       CALL almo_level2_exp3_1(qs_env,main_objectM1,&
               no_quench,main_objectM1%enter_b,&
               main_objectM1%enter,&
               assume_t0_q0x=.TRUE.,perturbation_only=.TRUE.,&
               special_case=1,error=error)
    CASE (almo_deloc_qx)
       main_objectM1%perturbative_delocalization=.TRUE.
       IF (main_objectM1%almo_update_algorithm_q.eq.almo_scf_diag) THEN
          DO ispin=1,main_objectM1%nspins
             CALL cp_dbcsr_copy(main_objectM1%enter(ispin),&
                     main_objectM1%enter_b(ispin),error=error)
          ENDDO
          CALL almo_level2_exp2_4(qs_env,main_objectM1,error)
       ELSE IF (main_objectM1%almo_update_algorithm_q.eq.almo_scf_pcg) THEN
          CALL almo_level2_exp3_1(qs_env,main_objectM1,&
                  main_objectM1%quench_t,main_objectM1%enter_b,&
                  main_objectM1%enter,&
                  assume_t0_q0x=.TRUE.,perturbation_only=.TRUE.,&
                  special_case=-1,error=error)
       ENDIF
    CASE (almo_deloc_qscf)
       main_objectM1%perturbative_delocalization=.FALSE.
       IF (main_objectM1%almo_update_algorithm_q.eq.almo_scf_diag) THEN
          DO ispin=1,main_objectM1%nspins
             CALL cp_dbcsr_copy(main_objectM1%enter(ispin),&
                     main_objectM1%enter_b(ispin),error=error)
          ENDDO
          CALL almo_level2_exp2_4(qs_env,main_objectM1,error)
       ELSE IF (main_objectM1%almo_update_algorithm_q.eq.almo_scf_pcg) THEN
          CALL almo_level2_exp3_1(qs_env,main_objectM1,&
                  main_objectM1%quench_t,main_objectM1%enter_b,&
                  main_objectM1%enter,&
                  assume_t0_q0x=.TRUE.,perturbation_only=.FALSE.,&
                  special_case=-1,error=error)
       ENDIF
    CASE DEFAULT
       CPErrorMessage(cp_failure_level,routineP,"illegal post scf method",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    
    SELECT CASE (main_objectM1%deloc_method)
    CASE (almo_deloc_full_scf,almo_deloc_x_then_full_scf)
       IF (main_objectM1%deloc_truncate_virt.ne.virt_full) THEN
          CPErrorMessage(cp_failure_level,routineP,"full scf is NYI for truncated virtual space",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       CALL almo_level2_exp3_1(qs_env,main_objectM1,&
               no_quench,main_objectM1%enter,&
               main_objectM1%enter,&
               assume_t0_q0x=.FALSE.,perturbation_only=.FALSE.,&
               special_case=1,error=error)
    END SELECT
    
    SELECT CASE (main_objectM1%deloc_method)
    CASE (almo_deloc_x,almo_deloc_full_scf,almo_deloc_x_then_full_scf)
       DO ispin=1, main_objectM1%nspins
          CALL cp_dbcsr_release(no_quench(ispin),error=error)
       ENDDO
       DEALLOCATE(no_quench)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE almo_level1_exp3_2

! *****************************************************************************
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level2_aux1(main_objectM1,matrix_s0,error)
    
    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s0
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level2_aux1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins

    CALL timeset(routineN,handle)
   
    nspins = main_objectM1%nspins

    CALL almo_levelX_spec6_0(matrix_new=main_objectM1%matrix_s(1),&
         matrix_qs=matrix_s0,&
         main_objectM1=main_objectM1,&
         name_new="S",&
         size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
         symmetry_new=dbcsr_type_symmetric,&
         spin_key=0,&
         init_domains=.FALSE.,&
         error=error)
    CALL almo_levelX_spec6_0(matrix_new=main_objectM1%so_b(1),&
         matrix_qs=matrix_s0,&
         main_objectM1=main_objectM1,&
         name_new="S_BLK",&
         size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
         symmetry_new=dbcsr_type_symmetric,&
         spin_key=0,&
         init_domains=.TRUE.,&
         error=error)
    IF (main_objectM1%almo_update_algorithm.eq.almo_scf_diag) THEN
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%so_bsi(1),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="S_BLK_SQRT_INV",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=0,&
            init_domains=.TRUE.,&
            error=error)
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%so_bs(1),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="S_BLK_SQRT",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=0,&
            init_domains=.TRUE.,&
            error=error)
    ELSE IF (main_objectM1%almo_update_algorithm.eq.almo_scf_dm_sign) THEN
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%so_bi(1),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="S_BLK_INV",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=0,&
            init_domains=.TRUE.,&
            error=error)
    ENDIF

    ALLOCATE(main_objectM1%enter_b(nspins))
    ALLOCATE(main_objectM1%quench_t_blk(nspins))
    ALLOCATE(main_objectM1%gst_b(nspins))
    ALLOCATE(main_objectM1%gst_xx(nspins))
    ALLOCATE(main_objectM1%rem(nspins))
    ALLOCATE(main_objectM1%rem_i(nspins))
    ALLOCATE(main_objectM1%rem_s(nspins))
    ALLOCATE(main_objectM1%rem_si(nspins))
    ALLOCATE(main_objectM1%rem_b(nspins))
    ALLOCATE(main_objectM1%enter(nspins))
    ALLOCATE(main_objectM1%enter_tr(nspins))
    DO ispin=1,nspins
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%quench_t_blk(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="Q_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%enter_b(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="T_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%gst_b(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="ERR_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%gst_xx(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="ERR_XX",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%enter_tr(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="T_TR",&
            size_keys=(/almo_mat_dim_occ,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%rem(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="SIG",&
            size_keys=(/almo_mat_dim_occ,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%rem_b(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="SIG_BLK",&
            size_keys=(/almo_mat_dim_occ,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       CALL almo_levelX_spec6_0(&
            matrix_new=main_objectM1%rem_i(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="SIGINV",&
            size_keys=(/almo_mat_dim_occ,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       CALL almo_levelX_spec6_0(&
            matrix_new=main_objectM1%enter(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="T",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       CALL cp_dbcsr_init(main_objectM1%rem_s(ispin),error=error)
       CALL cp_dbcsr_init(main_objectM1%rem_si(ispin),error=error)
       CALL cp_dbcsr_create(main_objectM1%rem_s(ispin),&
                            template=main_objectM1%rem(ispin),&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL cp_dbcsr_create(main_objectM1%rem_si(ispin),&
                            template=main_objectM1%rem(ispin),&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
    ENDDO
   
    IF (main_objectM1%need_virtuals) THEN
       ALLOCATE(main_objectM1%v_b(nspins))
       ALLOCATE(main_objectM1%v_fb(nspins))
       ALLOCATE(main_objectM1%vib(nspins))
       ALLOCATE(main_objectM1%vo(nspins))
       ALLOCATE(main_objectM1%xia(nspins))
       ALLOCATE(main_objectM1%ov(nspins))
       ALLOCATE(main_objectM1%ov_full(nspins))
       ALLOCATE(main_objectM1%rem_vv(nspins))
       ALLOCATE(main_objectM1%rem_vv_b(nspins))
       ALLOCATE(main_objectM1%rem_vv_s(nspins))
       ALLOCATE(main_objectM1%rem_vv_si(nspins))
       ALLOCATE(main_objectM1%vv_fb(nspins))
                   
       IF (main_objectM1%deloc_truncate_virt.ne.virt_full) THEN
          ALLOCATE(main_objectM1%kx_b(nspins))
          ALLOCATE(main_objectM1%kx_bo(nspins))
          ALLOCATE(main_objectM1%k_t(nspins))
          ALLOCATE(main_objectM1%v_d(nspins))
          ALLOCATE(main_objectM1%v_db(nspins))
          ALLOCATE(main_objectM1%ov_disc(nspins))
          ALLOCATE(main_objectM1%vv_db(nspins))
          ALLOCATE(main_objectM1%vv_d(nspins))
          ALLOCATE(main_objectM1%opt_k_t_dd(nspins))
          ALLOCATE(main_objectM1%opt_k_t_rr(nspins))
          ALLOCATE(main_objectM1%opt_k_denom(nspins))
       ENDIF

       DO ispin=1,nspins
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%v_fb(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="V_FULL_BLK",&
               size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt_full/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%v_b(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="V_BLK",&
               size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%vib(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="V",&
               size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%ov_full(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="OV_FULL",&
               size_keys=(/almo_mat_dim_occ,almo_mat_dim_virt_full/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%ov(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="OV",&
               size_keys=(/almo_mat_dim_occ,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%vo(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="VO",&
               size_keys=(/almo_mat_dim_virt,almo_mat_dim_occ/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%xia(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="VO",&
               size_keys=(/almo_mat_dim_virt,almo_mat_dim_occ/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%rem_vv(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="SIG_VV",&
               size_keys=(/almo_mat_dim_virt,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_symmetric,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%vv_fb(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="VV_FULL_BLK",&
               size_keys=(/almo_mat_dim_virt_full,almo_mat_dim_virt_full/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.TRUE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%rem_vv_b(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="SIG_VV_BLK",&
               size_keys=(/almo_mat_dim_virt,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_symmetric,&
               spin_key=ispin,&
               init_domains=.TRUE.,&
               error=error)
          CALL cp_dbcsr_init(main_objectM1%rem_vv_s(ispin),error=error)
          CALL cp_dbcsr_init(main_objectM1%rem_vv_si(ispin),error=error)
          CALL cp_dbcsr_create(main_objectM1%rem_vv_s(ispin),&
                  template=main_objectM1%rem_vv(ispin),&
                  matrix_type=dbcsr_type_no_symmetry, error=error)
          CALL cp_dbcsr_create(main_objectM1%rem_vv_si(ispin),&
                  template=main_objectM1%rem_vv(ispin),&
                  matrix_type=dbcsr_type_no_symmetry, error=error)
          
          IF (main_objectM1%deloc_truncate_virt.ne.virt_full) THEN
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%opt_k_t_rr(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="OPT_K_U_RR",&
                  size_keys=(/almo_mat_dim_virt,almo_mat_dim_virt/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%vv_d(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="VV_DISC",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_symmetric,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%opt_k_t_dd(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="OPT_K_U_DD",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%vv_db(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="VV_DISC_BLK",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_symmetric,&
                  spin_key=ispin,&
                  init_domains=.TRUE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%kx_b(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="K_BLK",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.TRUE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%kx_bo(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="K_BLK_1",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.TRUE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%opt_k_denom(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="OPT_K_DENOM",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%k_t(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="K_TR",&
                  size_keys=(/almo_mat_dim_virt,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%v_db(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="V_DISC_BLK",&
                  size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%v_d(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="V_DISC",&
                  size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL almo_levelX_spec6_0(matrix_new=main_objectM1%ov_disc(ispin),&
                  matrix_qs=matrix_s0,&
                  main_objectM1=main_objectM1,&
                  name_new="OV_DISC",&
                  size_keys=(/almo_mat_dim_occ,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
        
          ENDIF 

       ENDDO 
    ENDIF
    
    IF (main_objectM1%need_orbital_energies) THEN
       ALLOCATE(main_objectM1%eoo(nspins))
       ALLOCATE(main_objectM1%evv(nspins))
       DO ispin=1,nspins
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%eoo(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="E_OCC",&
               size_keys=(/almo_mat_dim_occ,almo_mat_dim_occ/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL almo_levelX_spec6_0(matrix_new=main_objectM1%evv(ispin),&
               matrix_qs=matrix_s0,&
               main_objectM1=main_objectM1,&
               name_new="E_VIRT",&
               size_keys=(/almo_mat_dim_virt_full,almo_mat_dim_virt_full/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
       ENDDO
    ENDIF

    ALLOCATE(main_objectM1%dpp(nspins))
    ALLOCATE(main_objectM1%dpp_b(nspins))
    ALLOCATE(main_objectM1%hfh(nspins))
    ALLOCATE(main_objectM1%hfh_b(nspins))
    IF (main_objectM1%need_previous_ks)&
       ALLOCATE(main_objectM1%hfh_con(nspins))
    DO ispin=1,nspins
       CALL cp_dbcsr_init(main_objectM1%dpp(ispin),error=error)
       CALL cp_dbcsr_create(main_objectM1%dpp(ispin),&
                            template=main_objectM1%matrix_s(1),&
                            matrix_type=dbcsr_type_symmetric, error=error)
       CALL cp_dbcsr_init(main_objectM1%hfh(ispin),error=error)
       CALL cp_dbcsr_create(main_objectM1%hfh(ispin),&
                            template=main_objectM1%matrix_s(1),&
                            matrix_type=dbcsr_type_symmetric, error=error)
       IF (main_objectM1%need_previous_ks) THEN
          CALL cp_dbcsr_init(main_objectM1%hfh_con(ispin),&
                          error=error)
          CALL cp_dbcsr_create(main_objectM1%hfh_con(ispin),&
                            template=main_objectM1%matrix_s(1),&
                            matrix_type=dbcsr_type_symmetric, error=error)
       ENDIF
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%dpp_b(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="P_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       CALL almo_levelX_spec6_0(matrix_new=main_objectM1%hfh_b(ispin),&
            matrix_qs=matrix_s0,&
            main_objectM1=main_objectM1,&
            name_new="KS_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
    ENDDO
    
    CALL timestop(handle)
  
  END SUBROUTINE almo_level2_aux1 

! *****************************************************************************
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_level1_exp4_1(qs_env,main_objectM1,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_objectM1_type)                 :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_level1_exp4_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    CALL cp_dbcsr_release(main_objectM1%matrix_s(1),error=error)
    CALL cp_dbcsr_release(main_objectM1%so_b(1),error=error)
    IF (main_objectM1%almo_update_algorithm.eq.almo_scf_diag) THEN
       CALL cp_dbcsr_release(main_objectM1%so_bsi(1),error=error)
       CALL cp_dbcsr_release(main_objectM1%so_bs(1),error=error)
    ELSE IF (main_objectM1%almo_update_algorithm.eq.almo_scf_dm_sign) THEN
       CALL cp_dbcsr_release(main_objectM1%so_bi(1),error=error)
    ENDIF
    DO ispin=1,main_objectM1%nspins
       CALL cp_dbcsr_release(main_objectM1%quench_t(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%quench_t_blk(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%enter_b(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%gst_b(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%gst_xx(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%enter_tr(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%rem(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%rem_b(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%rem_i(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%enter(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%rem_s(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%rem_si(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%dpp(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%hfh(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%dpp_b(ispin),error=error)
       CALL cp_dbcsr_release(main_objectM1%hfh_b(ispin),error=error)
       IF (main_objectM1%need_previous_ks) THEN
          CALL cp_dbcsr_release(main_objectM1%hfh_con(ispin),&
                  error=error)
       ENDIF
       IF (main_objectM1%need_virtuals) THEN
          CALL cp_dbcsr_release(main_objectM1%v_b(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%v_fb(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%vib(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%vo(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%xia(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%ov(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%ov_full(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%rem_vv(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%rem_vv_b(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%rem_vv_s(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%rem_vv_si(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%vv_fb(ispin),error=error)
          IF (main_objectM1%deloc_truncate_virt.ne.virt_full) THEN
             CALL cp_dbcsr_release(main_objectM1%k_t(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%kx_b(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%kx_bo(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%v_d(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%v_db(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%ov_disc(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%vv_db(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%vv_d(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%opt_k_t_dd(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%opt_k_t_rr(ispin),error=error)
             CALL cp_dbcsr_release(main_objectM1%opt_k_denom(ispin),error=error)
          ENDIF
       ENDIF
       IF (main_objectM1%need_orbital_energies) THEN
          CALL cp_dbcsr_release(main_objectM1%eoo(ispin),error=error)
          CALL cp_dbcsr_release(main_objectM1%evv(ispin),error=error)
       ENDIF
    ENDDO

    DEALLOCATE(main_objectM1%dpp)
    DEALLOCATE(main_objectM1%dpp_b)
    DEALLOCATE(main_objectM1%hfh)
    DEALLOCATE(main_objectM1%hfh_b)
    DEALLOCATE(main_objectM1%enter_b)
    DEALLOCATE(main_objectM1%gst_b)
    DEALLOCATE(main_objectM1%gst_xx)
    DEALLOCATE(main_objectM1%enter)
    DEALLOCATE(main_objectM1%enter_tr)
    DEALLOCATE(main_objectM1%rem)
    DEALLOCATE(main_objectM1%rem_b)
    DEALLOCATE(main_objectM1%rem_s)
    DEALLOCATE(main_objectM1%rem_si)
    DEALLOCATE(main_objectM1%rem_i)
    DEALLOCATE(main_objectM1%quench_t)
    DEALLOCATE(main_objectM1%quench_t_blk)
    IF (main_objectM1%need_virtuals) THEN
       DEALLOCATE(main_objectM1%v_b)
       DEALLOCATE(main_objectM1%v_fb)
       DEALLOCATE(main_objectM1%vib)
       DEALLOCATE(main_objectM1%vo)
       DEALLOCATE(main_objectM1%xia)
       DEALLOCATE(main_objectM1%ov)
       DEALLOCATE(main_objectM1%ov_full)
       DEALLOCATE(main_objectM1%rem_vv)
       DEALLOCATE(main_objectM1%rem_vv_b)
       DEALLOCATE(main_objectM1%rem_vv_s)
       DEALLOCATE(main_objectM1%rem_vv_si)
       DEALLOCATE(main_objectM1%vv_fb)
       IF (main_objectM1%deloc_truncate_virt.ne.virt_full) THEN
          DEALLOCATE(main_objectM1%k_t)
          DEALLOCATE(main_objectM1%kx_b)
          DEALLOCATE(main_objectM1%v_d)
          DEALLOCATE(main_objectM1%v_db)
          DEALLOCATE(main_objectM1%ov_disc)
          DEALLOCATE(main_objectM1%vv_db)
          DEALLOCATE(main_objectM1%vv_d)
          DEALLOCATE(main_objectM1%kx_bo)
          DEALLOCATE(main_objectM1%opt_k_t_dd)
          DEALLOCATE(main_objectM1%opt_k_t_rr)
          DEALLOCATE(main_objectM1%opt_k_denom)
       ENDIF
    ENDIF
    IF (main_objectM1%need_previous_ks) THEN
       DEALLOCATE(main_objectM1%hfh_con)
    ENDIF
    IF (main_objectM1%need_orbital_energies) THEN
       DEALLOCATE(main_objectM1%eoo)
       DEALLOCATE(main_objectM1%evv)
    ENDIF

    DO ispin=1,main_objectM1%nspins
       CALL release_object01_gen(&
              main_objectM1%concept_p(:,ispin),&
              error=error)
       CALL release_object01_gen(main_objectM1%concept_soi(:,ispin),&
              error=error)
       CALL release_object01_gen(main_objectM1%concept_ss_inv(:,ispin),&
              error=error)
       CALL release_object01_gen(main_objectM1%concept_ss(:,ispin),&
              error=error)
       CALL release_object01_gen(main_objectM1%concept_hfh_xx(:,ispin),&
              error=error)
       CALL release_object01_gen(main_objectM1%concept_enter(:,ispin),&
              error=error)
       CALL release_object01_gen(main_objectM1%concept_gst(:,ispin),&
              error=error)
       CALL release_object01_gen(main_objectM1%concept_dpp(:,ispin),&
              error=error)
    ENDDO
    DEALLOCATE(main_objectM1%concept_p)
    DEALLOCATE(main_objectM1%concept_soi)
    DEALLOCATE(main_objectM1%concept_ss_inv)
    DEALLOCATE(main_objectM1%concept_ss)
    DEALLOCATE(main_objectM1%concept_hfh_xx)
    DEALLOCATE(main_objectM1%concept_enter)
    DEALLOCATE(main_objectM1%concept_gst)
    DEALLOCATE(main_objectM1%concept_dpp)
    DO ispin=1,main_objectM1%nspins
       DEALLOCATE(main_objectM1%domain_map(ispin)%pairs)
       DEALLOCATE(main_objectM1%domain_map(ispin)%index1)
    ENDDO
    DEALLOCATE(main_objectM1%domain_map)
    DEALLOCATE(main_objectM1%domain_index_of_ao)
    DEALLOCATE(main_objectM1%domain_index_of_atom)
    DEALLOCATE(main_objectM1%first_atom_of_domain)
    DEALLOCATE(main_objectM1%nbasis_of_domain)
    DEALLOCATE(main_objectM1%nocc_of_domain)
    DEALLOCATE(main_objectM1%nvirt_full_of_domain)
    DEALLOCATE(main_objectM1%nvirt_of_domain)
    DEALLOCATE(main_objectM1%nvirt_disc_of_domain)
    DEALLOCATE(main_objectM1%mu_of_domain)
    DEALLOCATE(main_objectM1%cpu_of_domain)
DEALLOCATE(main_objectM1%charge_of_domain)
    
    DEALLOCATE(main_objectM1%domain_index_of_ao_block)
    DEALLOCATE(main_objectM1%domain_index_of_mo_block)

    CALL cp_para_env_release(main_objectM1%para_env,error)
    CALL cp_blacs_env_release(main_objectM1%blacs_env,error)

    CALL timestop(handle)

  END SUBROUTINE almo_level1_exp4_1

END MODULE almo_scf

