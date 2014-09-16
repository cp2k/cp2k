!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for all ALMO-based SCF methods
!>        RZK-warning marks unresolved issues
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin 
! *****************************************************************************
MODULE almo_scf
  USE almo_scf_diis_types,             ONLY: almo_scf_diis_extrapolate,&
                                             almo_scf_diis_init,&
                                             almo_scf_diis_push,&
                                             almo_scf_diis_release,&
                                             almo_scf_diis_type
  USE almo_scf_methods,                ONLY: &
       almo_scf_domain_operations, almo_scf_ks_blk_to_tv_blk, &
       almo_scf_ks_to_ks_blk, almo_scf_ks_to_ks_xx, almo_scf_ks_xx_to_tv_xx, &
       almo_scf_p_blk_to_t_blk, almo_scf_t_blk_to_p, &
       almo_scf_t_blk_to_t_blk_orthonormal, almo_scf_t_to_p, &
       apply_domain_operators, apply_projector, &
       construct_domain_preconditioner, construct_domain_r_down, &
       construct_domain_s_inv, construct_domain_s_sqrt, distribute_domains, &
       get_overlap, invert_diagonal_blk, newton_grad_to_step
!almo_scf_v_to_v_orthonormal_blk,&
!most_orthogonal_almos,&
  USE almo_scf_qs,                     ONLY: almo_scf_construct_quencher,&
                                             almo_scf_dm_to_ks,&
                                             almo_scf_init_qs,&
                                             almo_scf_update_ks_energy,&
                                             matrix_almo_create,&
                                             matrix_qs_to_almo
  USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                             almo_mat_dim_occ,&
                                             almo_mat_dim_virt,&
                                             almo_mat_dim_virt_disc,&
                                             almo_mat_dim_virt_full,&
                                             almo_max_cutoff_multiplier,&
                                             almo_scf_env_type
  USE bibliography,                    ONLY: Khaliullin2013,&
                                             cite_reference
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                             cp_dbcsr_cholesky_invert,&
                                             cp_dbcsr_cholesky_restore
  USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_frobenius_norm, &
       cp_dbcsr_function_of_elements, cp_dbcsr_get_diag, cp_dbcsr_get_info, &
       cp_dbcsr_get_stored_coordinates, cp_dbcsr_hadamard_product, &
       cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_nblkcols_total, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_norm, cp_dbcsr_p_type, &
       cp_dbcsr_print, cp_dbcsr_print_block_sum, cp_dbcsr_release, &
       cp_dbcsr_reserve_block2d, cp_dbcsr_scale, cp_dbcsr_set, &
       cp_dbcsr_set_diag, cp_dbcsr_trace, cp_dbcsr_triu, cp_dbcsr_type, &
       cp_dbcsr_work_create, dbcsr_distribution_mp, dbcsr_func_artanh, &
       dbcsr_func_dtanh, dbcsr_func_inverse, dbcsr_func_inverse_special, &
       dbcsr_func_tanh, dbcsr_mp_mynode, dbcsr_mp_numnodes, &
       dbcsr_norm_maxabsnorm, dbcsr_type_no_symmetry, dbcsr_type_symmetric
  USE cp_external_control,             ONLY: external_control
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE ct_methods,                      ONLY: analytic_line_search,&
                                             ct_step_execute,&
                                             diagonalize_diagonal_blocks
!USE recursive_matrix,                ONLY: sqrti_factorization
  USE ct_types,                        ONLY: ct_step_env_clean,&
                                             ct_step_env_get,&
                                             ct_step_env_init,&
                                             ct_step_env_set,&
                                             ct_step_env_type
  USE domain_submatrix_methods,        ONLY: add_submatrices,&
                                             construct_submatrices,&
                                             copy_submatrices,&
                                             init_submatrices,&
                                             maxnorm_submatrices,&
                                             release_submatrices
  USE domain_submatrix_types,          ONLY: domain_submatrix_type,&
                                             select_row
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
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
!USE dm_,               ONLY: density_matrix_sign
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
!> \brief The entry point into ALMO SCF routines
!> \param qs_env - pointer to the QS environment
!> \param calc_forces - calculate forces?
!> \param error - error
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
    TYPE(almo_scf_env_type)                  :: almo_scf_env

    CALL timeset(routineN,handle)

    CALL cite_reference(Khaliullin2013)

    ! no forces so far 
    CALL cp_assert(.NOT.calc_forces,cp_failure_level,cp_assertion_failed,routineP,&
           " No forces/gradients with ALMOs.... yet(!) "//&
                CPSourceFileRef, only_ionode=.TRUE.)

    ! initialize scf
    CALL almo_scf_init(qs_env,almo_scf_env,error)

    ! perform SCF for block diagonal ALMOs
    CALL almo_scf_main(qs_env,almo_scf_env,error)

    ! allow electron delocalization
    CALL almo_scf_delocalization(qs_env,almo_scf_env,error)

    ! return computed quantities to the qs_env
    !CALL almo_return_to_qs(qs_env,almo_scf_env,error)

    ! do post scf processing
    CALL almo_scf_clean_up(qs_env,almo_scf_env,error)

    CALL timestop(handle)

  END SUBROUTINE almo_entry_scf

! *****************************************************************************
!> \brief Initialization of the almo_scf_env_type.
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_init(qs_env,almo_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_init', &
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

    ! get info from the qs_env
    CALL get_qs_env(qs_env,&
      nelectron_total=almo_scf_env%nelectrons_total,&
      matrix_s=matrix_s,&
      dft_control=dft_control,&
      molecule_set=molecule_set,&
      input=input,&
      has_unit_metric=almo_scf_env%orthogonal_basis,&
      para_env=almo_scf_env%para_env,&
      blacs_env=almo_scf_env%blacs_env,&
      nelectron_spin=almo_scf_env%nelectrons_spin,&
      error=error)
    CALL cp_para_env_retain(almo_scf_env%para_env,error)
    CALL cp_blacs_env_retain(almo_scf_env%blacs_env, error)

    ! copy basic quantities
    almo_scf_env%nspins=dft_control%nspins
    almo_scf_env%natoms=cp_dbcsr_nblkrows_total(matrix_s(1)%matrix)
    almo_scf_env%nmolecules=SIZE(molecule_set)
    CALL cp_dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=naos )
    almo_scf_env%naos=naos

    ! convenient local varibales
    nspins=almo_scf_env%nspins
    nmols=almo_scf_env%nmolecules
    natoms=almo_scf_env%natoms
    
    ! parse the almo_scf section and set appropriate quantities
    CALL almo_scf_init_read_write_input(input,almo_scf_env,unit_nr,error)
   
    ! Define groups: either atomic or molecular
    IF (almo_scf_env%domain_layout_mos==almo_domain_layout_molecular) THEN
       almo_scf_env%ndomains=almo_scf_env%nmolecules
    ELSE
       almo_scf_env%ndomains=almo_scf_env%natoms
    ENDIF

    ! allocate domain descriptors
    ndomains=almo_scf_env%ndomains
    ALLOCATE(almo_scf_env%domain_index_of_atom(natoms))
    ALLOCATE(almo_scf_env%domain_index_of_ao(naos))
    ALLOCATE(almo_scf_env%first_atom_of_domain(ndomains))
    ALLOCATE(almo_scf_env%nbasis_of_domain(ndomains))
    ALLOCATE(almo_scf_env%nocc_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%nvirt_full_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%nvirt_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%nvirt_disc_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%mu_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%cpu_of_domain(ndomains))
    ! RZK-warning charges on domains are not yet implemented
    ALLOCATE(almo_scf_env%charge_of_domain(ndomains))

    ! fill out domain descriptors and group descriptors
    IF (almo_scf_env%domain_layout_mos==almo_domain_layout_molecular) THEN
       ! get domain info from molecule_set 
       CALL get_molecule_set_info(molecule_set,&
               atom_to_mol=almo_scf_env%domain_index_of_atom,&
               mol_to_first_atom=almo_scf_env%first_atom_of_domain,&
               mol_to_nelectrons=almo_scf_env%nocc_of_domain,&
               mol_to_nbasis=almo_scf_env%nbasis_of_domain,&
               error=error)
       ! RZK-warning adjust number of electrons on each molecule taking into account its charge
       !almo_scf_env%nocc_of_domain(:)
       DO ispin=1,nspins
          ! take care of the full virtual subspace
          almo_scf_env%nvirt_full_of_domain(:,ispin)=&
             almo_scf_env%nbasis_of_domain(:)-&
             almo_scf_env%nocc_of_domain(:,ispin)
          ! and the truncated virtual subspace
          SELECT CASE (almo_scf_env%deloc_truncate_virt)
          CASE (virt_full)
             almo_scf_env%nvirt_of_domain(:,ispin)=&
                almo_scf_env%nvirt_full_of_domain(:,ispin)
             almo_scf_env%nvirt_disc_of_domain(:,ispin)=0
          CASE (virt_number)
             DO idomain=1,ndomains
                almo_scf_env%nvirt_of_domain(idomain,ispin)=&
                   MIN(almo_scf_env%deloc_virt_per_domain,&
                   almo_scf_env%nvirt_full_of_domain(idomain,ispin))
                almo_scf_env%nvirt_disc_of_domain(idomain,ispin)=&
                   almo_scf_env%nvirt_full_of_domain(idomain,ispin)-&
                   almo_scf_env%nvirt_of_domain(idomain,ispin)
             ENDDO
          CASE (virt_occ_size)
             DO idomain=1,ndomains
                almo_scf_env%nvirt_of_domain(idomain,ispin)=&
                   MIN(almo_scf_env%nocc_of_domain(idomain,ispin),&
                   almo_scf_env%nvirt_full_of_domain(idomain,ispin))
                almo_scf_env%nvirt_disc_of_domain(idomain,ispin)=&
                   almo_scf_env%nvirt_full_of_domain(idomain,ispin)-&
                   almo_scf_env%nvirt_of_domain(idomain,ispin)
             ENDDO
          CASE DEFAULT
             CPErrorMessage(cp_failure_level,routineP,"illegal method for virtual space truncation",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       ENDDO ! spin
    ELSE ! domains are atomic
       ! RZK-warning do the same for atomic domains/groups 
       almo_scf_env%domain_index_of_atom(1:natoms)=(/(i, i=1,natoms)/)
    ENDIF
    
    ao=1
    DO idomain=1,ndomains
      DO iao=1,almo_scf_env%nbasis_of_domain(idomain)
         almo_scf_env%domain_index_of_ao(ao)=idomain
         ao=ao+1
      ENDDO
    ENDDO

    almo_scf_env%mu_of_domain(:,:)=almo_scf_env%mu

    ! build domain (i.e. layout) indices for distribution blocks 
    ! ao blocks
    IF (almo_scf_env%mat_distr_aos==almo_mat_distr_atomic) THEN
     ALLOCATE(almo_scf_env%domain_index_of_ao_block(natoms))
     almo_scf_env%domain_index_of_ao_block(:)=&
        almo_scf_env%domain_index_of_atom(:)
    ELSE IF (almo_scf_env%mat_distr_aos==almo_mat_distr_molecular) THEN
     ALLOCATE(almo_scf_env%domain_index_of_ao_block(nmols))
     ! if distr blocks are molecular then domain layout is also molecular
     almo_scf_env%domain_index_of_ao_block(:)=(/(i, i=1,nmols)/)
    ENDIF
    ! mo blocks
    IF (almo_scf_env%mat_distr_mos==almo_mat_distr_atomic) THEN
     ALLOCATE(almo_scf_env%domain_index_of_mo_block(natoms))
     almo_scf_env%domain_index_of_mo_block(:)=&
        almo_scf_env%domain_index_of_atom(:)
    ELSE IF (almo_scf_env%mat_distr_mos==almo_mat_distr_molecular) THEN
     ALLOCATE(almo_scf_env%domain_index_of_mo_block(nmols))
     ! if distr blocks are molecular then domain layout is also molecular
     almo_scf_env%domain_index_of_mo_block(:)=(/(i, i=1,nmols)/)
    ENDIF

    ! set all flags
    !almo_scf_env%need_previous_ks=.FALSE.
    !IF (almo_scf_env%deloc_method==almo_deloc_harris) THEN
       almo_scf_env%need_previous_ks=.TRUE.
    !ENDIF

    !almo_scf_env%need_virtuals=.FALSE.
    !almo_scf_env%need_orbital_energies=.FALSE.
    !IF (almo_scf_env%almo_update_algorithm==almo_scf_diag) THEN
       almo_scf_env%need_virtuals=.TRUE.
       almo_scf_env%need_orbital_energies=.TRUE.
    !ENDIF

    ! FINISH setting job parameters here, print out job info
    CALL almo_scf_print_job_info(almo_scf_env,unit_nr,error)

    ! create all matrices
    CALL almo_scf_env_create_matrices(almo_scf_env,matrix_s(1)%matrix,error)
    ! set up matrix S and all required functions of S
    almo_scf_env%s_inv_done=.FALSE.
    almo_scf_env%s_sqrt_done=.FALSE.
    CALL almo_scf_init_ao_overlap(matrix_s(1)%matrix,almo_scf_env,error)

    ! create the quencher (imposes sparsity template)
    CALL almo_scf_construct_quencher(qs_env,almo_scf_env,error)
    CALL distribute_domains(almo_scf_env,error)

    ! allocate and init the domain preconditioner
    ALLOCATE(almo_scf_env%domain_preconditioner(ndomains,nspins))
    CALL init_submatrices(almo_scf_env%domain_preconditioner,error)
   
    ! allocate and init projected KS for domains
    ALLOCATE(almo_scf_env%domain_ks_xx(ndomains,nspins))
    CALL init_submatrices(almo_scf_env%domain_ks_xx,error)
   
    ! init ao-overlap subblocks
    ALLOCATE(almo_scf_env%domain_s_inv(ndomains,nspins))
    CALL init_submatrices(almo_scf_env%domain_s_inv,error)
    ALLOCATE(almo_scf_env%domain_s_sqrt_inv(ndomains,nspins))
    CALL init_submatrices(almo_scf_env%domain_s_sqrt_inv,error)
    ALLOCATE(almo_scf_env%domain_s_sqrt(ndomains,nspins))
    CALL init_submatrices(almo_scf_env%domain_s_sqrt,error)
    ALLOCATE(almo_scf_env%domain_t(ndomains,nspins))
    CALL init_submatrices(almo_scf_env%domain_t,error=error)
    ALLOCATE(almo_scf_env%domain_err(ndomains,nspins))
    CALL init_submatrices(almo_scf_env%domain_err,error=error)
    ALLOCATE(almo_scf_env%domain_r_down_up(ndomains,nspins))
    CALL init_submatrices(almo_scf_env%domain_r_down_up,error=error)

    ! initialization of the QS settings with the ALMO flavor
    CALL almo_scf_init_qs(qs_env,almo_scf_env,error)

    CALL timestop(handle)

  END SUBROUTINE almo_scf_init

! *****************************************************************************
!> \brief Parses the ALMO input section
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_init_read_write_input(input,almo_scf_env,unit_nr,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'almo_scf_init_read_write_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: almo_scf_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    almo_scf_section => section_vals_get_subs_vals(input,"DFT%ALMO_SCF",error=error)
    !orthogonalization_opt_section => section_vals_get_subs_vals(input,&
    !   "DFT%ALMO_SCF%ORTHOGONALIZATION_OPT",error=error)
   
    ! RZK-warning the values of these keywords are hardcoded 
    ! but ideally they should also be read from the input file
    almo_scf_env%order_lanczos=3
    almo_scf_env%eps_lanczos=1.0E-4_dp
    almo_scf_env%max_iter_lanczos=40

    ! read user input
    CALL section_vals_val_get(almo_scf_section,"EPS_FILTER",&
                              r_val=almo_scf_env%eps_filter,error=error)
    CALL section_vals_val_get(almo_scf_section,"BLOCKED_EPS_ITER",&
                              r_val=almo_scf_env%eps_scf_bd,error=error)
    CALL section_vals_val_get(almo_scf_section,"BLOCKED_MAX_ITER",&
                              i_val=almo_scf_env%max_scf_bd,error=error)
    CALL section_vals_val_get(almo_scf_section,"BLOCKED_N_DIIS",&
                              i_val=almo_scf_env%ndiis_bd,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_EPS_ITER",&
                              r_val=almo_scf_env%eps_scf_q,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_EPS_LIN_SEARCH",&
                              r_val=almo_scf_env%eps_lin_search,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_MAX_ITER",&
                              i_val=almo_scf_env%max_scf_q,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_CONJUGATOR",&
                              i_val=almo_scf_env%scf_conjugator,&
                              error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_METHOD",&
                              i_val=almo_scf_env%deloc_method,error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_ALGORITHM",&
                              i_val=almo_scf_env%almo_update_algorithm_q,&
                              error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_R_CUTOFF_FACTOR",&
                              r_val=almo_scf_env%quencher_r0_factor,&
                              error=error)

    ! do not do EDA in the official release before the output files are written
    ! using the proper CP2K routines
    almo_scf_env%almo_eda=0

    ! set to some small value because the xALMO DIIS procedure does not converge
    ! and currently is not used
    almo_scf_env%ndiis_q=2

    !CALL section_vals_val_get(almo_scf_section,"OUTER_MAX_SCF_Q",&
    !                          i_val=almo_scf_env%outer_max_scf_q,error=error)
    !CALL section_vals_val_get(almo_scf_section,"DOMAIN_LAYOUT_AOS",&
    !                          i_val=almo_scf_env%domain_layout_aos,error=error)
    !CALL section_vals_val_get(almo_scf_section,"DOMAIN_LAYOUT_MOS",&
    !                          i_val=almo_scf_env%domain_layout_mos,error=error)
    !CALL section_vals_val_get(almo_scf_section,"MATRIX_CLUSTERING_AOS",&
    !                          i_val=almo_scf_env%mat_distr_aos,error=error)
    !CALL section_vals_val_get(almo_scf_section,"MATRIX_CLUSTERING_MOS",&
    !                          i_val=almo_scf_env%mat_distr_mos,error=error)
    !CALL section_vals_val_get(almo_scf_section,"CONSTRAINT_TYPE",&
    !                          i_val=almo_scf_env%constraint_type,error=error)
    !CALL section_vals_val_get(almo_scf_section,"MU",&
    !                          r_val=almo_scf_env%mu,error=error)
    !CALL section_vals_val_get(almo_scf_section,"FIXED_MU",&
    !                          l_val=almo_scf_env%fixed_mu,error=error)
    !CALL section_vals_val_get(almo_scf_section,"EPS_USE_PREV_AS_GUESS",&
    !                          r_val=almo_scf_env%eps_prev_guess,error=error)
    !CALL section_vals_val_get(almo_scf_section,"MIXING_FRACTION",&
    !                          r_val=almo_scf_env%mixing_fraction,error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_TENSOR_TYPE",&
    !                          i_val=almo_scf_env%deloc_cayley_tensor_type,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_CONJUGATOR",&
    !                          i_val=almo_scf_env%deloc_cayley_conjugator,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_MAX_ITER",&
    !                          i_val=almo_scf_env%deloc_cayley_max_iter,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_USE_OCC_ORBS",&
    !                          l_val=almo_scf_env%deloc_use_occ_orbs,error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_USE_VIRT_ORBS",&
    !                          l_val=almo_scf_env%deloc_cayley_use_virt_orbs,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_LINEAR",&
    !                          l_val=almo_scf_env%deloc_cayley_linear,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_EPS_CONVERGENCE",&
    !                          r_val=almo_scf_env%deloc_cayley_eps_convergence,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_OCC_PRECOND",&
    !                          l_val=almo_scf_env%deloc_cayley_occ_precond,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_CAYLEY_VIR_PRECOND",&
    !                          l_val=almo_scf_env%deloc_cayley_vir_precond,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"ALMO_UPDATE_ALGORITHM_BD",&
    !                          i_val=almo_scf_env%almo_update_algorithm,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_TRUNCATE_VIRTUALS",&
    !                          i_val=almo_scf_env%deloc_truncate_virt,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"DELOC_VIRT_PER_DOMAIN",&
    !                          i_val=almo_scf_env%deloc_virt_per_domain,&
    !                          error=error)
    !
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_EPS_CONVERGENCE",&
    !                          r_val=almo_scf_env%opt_k_eps_convergence,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_MAX_ITER",&
    !                          i_val=almo_scf_env%opt_k_max_iter,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_OUTER_MAX_ITER",&
    !                          i_val=almo_scf_env%opt_k_outer_max_iter,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_TRIAL_STEP_SIZE",&
    !                          r_val=almo_scf_env%opt_k_trial_step_size,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJUGATOR",&
    !                          i_val=almo_scf_env%opt_k_conjugator,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_TRIAL_STEP_SIZE_MULTIPLIER",&
    !                          r_val=almo_scf_env%opt_k_trial_step_size_multiplier,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJ_ITER_START",&
    !                          i_val=almo_scf_env%opt_k_conj_iter_start,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_PREC_ITER_START",&
    !                          i_val=almo_scf_env%opt_k_prec_iter_start,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_CONJ_ITER_FREQ_RESET",&
    !                          i_val=almo_scf_env%opt_k_conj_iter_freq,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"OPT_K_PREC_ITER_FREQ_UPDATE",&
    !                          i_val=almo_scf_env%opt_k_prec_iter_freq,&
    !                          error=error)
    !
    !CALL section_vals_val_get(almo_scf_section,"QUENCHER_RADIUS_TYPE",&
    !                          i_val=almo_scf_env%quencher_radius_type,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"QUENCHER_R0_FACTOR",&
    !                          r_val=almo_scf_env%quencher_r0_factor,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"QUENCHER_R1_FACTOR",&
    !                          r_val=almo_scf_env%quencher_r1_factor,&
    !                          error=error)
    !!CALL section_vals_val_get(almo_scf_section,"QUENCHER_R0_SHIFT",&
    !!                          r_val=almo_scf_env%quencher_r0_shift,&
    !!                          error=error)
    !!CALL section_vals_val_get(almo_scf_section,"QUENCHER_R1_SHIFT",&
    !!                          r_val=almo_scf_env%quencher_r1_shift,&
    !!                          error=error)
    !!almo_scf_env%quencher_r0_shift = cp_unit_to_cp2k(&
    !!   almo_scf_env%quencher_r0_shift,"angstrom",error=error)
    !!almo_scf_env%quencher_r1_shift = cp_unit_to_cp2k(&
    !!   almo_scf_env%quencher_r1_shift,"angstrom",error=error)
    !
    !CALL section_vals_val_get(almo_scf_section,"QUENCHER_AO_OVERLAP_0",&
    !                          r_val=almo_scf_env%quencher_s0,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"QUENCHER_AO_OVERLAP_1",&
    !                          r_val=almo_scf_env%quencher_s1,&
    !                          error=error)

    !CALL section_vals_val_get(almo_scf_section,"ENVELOPE_AMPLITUDE",&
    !                          r_val=almo_scf_env%envelope_amplitude,&
    !                          error=error)

    !CALL section_vals_val_get(almo_scf_section,"LOGICAL01",&
    !                          l_val=almo_scf_env%logical01,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"LOGICAL02",&
    !                          l_val=almo_scf_env%logical02,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"LOGICAL03",&
    !                          l_val=almo_scf_env%logical03,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"LOGICAL04",&
    !                          l_val=almo_scf_env%logical04,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"LOGICAL05",&
    !                          l_val=almo_scf_env%logical05,&
    !                          error=error)

    !CALL section_vals_val_get(almo_scf_section,"REAL01",&
    !                          r_val=almo_scf_env%real01,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"REAL02",&
    !                          r_val=almo_scf_env%real02,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"REAL03",&
    !                          r_val=almo_scf_env%real03,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"REAL04",&
    !                          r_val=almo_scf_env%real04,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"REAL05",&
    !                          r_val=almo_scf_env%real05,&
    !                          error=error)

    !CALL section_vals_val_get(almo_scf_section,"INTEGER01",&
    !                          i_val=almo_scf_env%integer01,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"INTEGER02",&
    !                          i_val=almo_scf_env%integer02,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"INTEGER03",&
    !                          i_val=almo_scf_env%integer03,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"INTEGER04",&
    !                          i_val=almo_scf_env%integer04,&
    !                          error=error)
    !CALL section_vals_val_get(almo_scf_section,"INTEGER05",&
    !                          i_val=almo_scf_env%integer05,&
    !                          error=error)

    !! read molecular charges 
    !CALL section_vals_val_get(almo_scf_section,"INT_LIST01", &
    !        n_rep_val=n_rep,error=error)
    !counter_i = 0
    !DO k = 1,n_rep
    !  CALL section_vals_val_get(almo_scf_section,"INT_LIST01",&
    !          i_rep_val=k,i_vals=tmplist, error=error)
    !   DO jj = 1,SIZE(tmplist)
    !      counter_i=counter_i+1
    !      almo_scf_env%charge_of_domain(counter_i)=tmplist(jj)
    !   ENDDO
    !ENDDO

    !! subsection ORTHOGONALIZATION_OPT
    !almo_scf_env%orthogonalization_opt%eps_filter=almo_scf_env%eps_filter
    !almo_scf_env%orthogonalization_opt%para_env => almo_scf_env%para_env
    !almo_scf_env%orthogonalization_opt%blacs_env => almo_scf_env%blacs_env
    !CALL section_vals_val_get(orthogonalization_opt_section,"PRECONDITIONER",&
    !                          i_val=almo_scf_env%orthogonalization_opt%preconditioner,&
    !                          error=error)
    !CALL section_vals_val_get(orthogonalization_opt_section,"CONJUGATOR",&
    !                          i_val=almo_scf_env%orthogonalization_opt%conjugator,&
    !                          error=error)
    !CALL section_vals_val_get(orthogonalization_opt_section,"MAX_ITER",&
    !                          i_val=almo_scf_env%orthogonalization_opt%max_iter,&
    !                          error=error)
    !CALL section_vals_val_get(orthogonalization_opt_section,"GRAD_TARGET",&
    !                          r_val=almo_scf_env%orthogonalization_opt%grad_target,&
    !                          error=error)
    !CALL section_vals_val_get(orthogonalization_opt_section,"LIN_SEARCH_GRAD_TARGET",&
    !                          r_val=almo_scf_env%orthogonalization_opt%lin_search_grad_target,&
    !                          error=error)
    !CALL section_vals_val_get(orthogonalization_opt_section,"STEP_SIZE_GUESS",&
    !                          r_val=almo_scf_env%orthogonalization_opt%step_size_guess,&
    !                          error=error)

    almo_scf_env%domain_layout_aos=almo_domain_layout_molecular
    almo_scf_env%domain_layout_mos=almo_domain_layout_molecular
    almo_scf_env%mat_distr_aos=almo_mat_distr_molecular 
    almo_scf_env%mat_distr_mos=almo_mat_distr_molecular

    almo_scf_env%outer_max_scf_q=0
    almo_scf_env%constraint_type=almo_constraint_distance
    almo_scf_env%mu=-0.1_dp
    almo_scf_env%fixed_mu=.FALSE.
    almo_scf_env%almo_update_algorithm=almo_scf_diag
    almo_scf_env%eps_prev_guess=almo_scf_env%eps_filter/10.0_dp
    almo_scf_env%mixing_fraction=0.45_dp
    
    almo_scf_env%deloc_cayley_tensor_type=tensor_orthogonal
    almo_scf_env%deloc_cayley_conjugator=cg_hager_zhang
    almo_scf_env%deloc_cayley_max_iter=100
    almo_scf_env%deloc_use_occ_orbs=.TRUE.
    almo_scf_env%deloc_cayley_use_virt_orbs=.FALSE.
    almo_scf_env%deloc_cayley_linear=.FALSE.
    almo_scf_env%deloc_cayley_eps_convergence=1.0E-6_dp
    almo_scf_env%deloc_cayley_occ_precond=.TRUE.
    almo_scf_env%deloc_cayley_vir_precond=.TRUE.
    almo_scf_env%deloc_truncate_virt=virt_full
    almo_scf_env%deloc_virt_per_domain=-1
    
    almo_scf_env%opt_k_eps_convergence=1.0E-5_dp
    almo_scf_env%opt_k_max_iter=100
    almo_scf_env%opt_k_outer_max_iter=1
    almo_scf_env%opt_k_trial_step_size=0.05_dp
    almo_scf_env%opt_k_conjugator=cg_hager_zhang
    almo_scf_env%opt_k_trial_step_size_multiplier=1.05_dp
    almo_scf_env%opt_k_conj_iter_start=0
    almo_scf_env%opt_k_prec_iter_start=0
    almo_scf_env%opt_k_conj_iter_freq=10000000
    almo_scf_env%opt_k_prec_iter_freq=10000000
    
    almo_scf_env%quencher_radius_type=do_bondparm_vdw
    almo_scf_env%quencher_r1_factor=almo_scf_env%quencher_r0_factor
    !almo_scf_env%quencher_r0_shift=0.0_dp
    !almo_scf_env%quencher_r1_shift=0.0_dp
    !almo_scf_env%quencher_r0_shift = cp_unit_to_cp2k(&
    !   almo_scf_env%quencher_r0_shift,"angstrom",error=error)
    !almo_scf_env%quencher_r1_shift = cp_unit_to_cp2k(&
    !   almo_scf_env%quencher_r1_shift,"angstrom",error=error)
    
    almo_scf_env%quencher_s0=1.0E-4_dp
    almo_scf_env%quencher_s1=1.0E-6_dp

    almo_scf_env%envelope_amplitude=1.0_dp

    almo_scf_env%logical01=.FALSE. ! md in eDOF space
    almo_scf_env%logical02=.TRUE.  ! not used
    almo_scf_env%logical03=.TRUE.  ! not used 
    almo_scf_env%logical04=.TRUE.  ! use preconditioner
    almo_scf_env%logical05=.FALSE. ! optimize theta

    almo_scf_env%real01=almo_scf_env%eps_filter/10.0_dp ! skip gradients
    almo_scf_env%real02=0.0_dp ! not used
    almo_scf_env%real03=0.0_dp ! not used
    almo_scf_env%real04=0.5_dp ! mixing s-f precond
    almo_scf_env%real05=0.1_dp ! step size guess

    almo_scf_env%integer01=10 ! start eDOF-md
    almo_scf_env%integer02=4  ! preconditioner type
    almo_scf_env%integer03=0  ! not used
    almo_scf_env%integer04=0  ! fixed number of line searches (no grad)
    almo_scf_env%integer05=0  ! not used
          
    IF (almo_scf_env%almo_update_algorithm_q.eq.almo_scf_diag .AND. &
        almo_scf_env%deloc_method.ne.almo_deloc_qx) THEN
      CPErrorMessage(cp_failure_level,routineP,"DOMAIN-DIAG is implemented only for X_R",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (almo_scf_env%deloc_truncate_virt.EQ.virt_number .AND. &
        almo_scf_env%deloc_virt_per_domain.LE.0) THEN
      CPErrorMessage(cp_failure_level,routineP,"specify a positive number of virtual orbitals",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (almo_scf_env%deloc_truncate_virt.EQ.virt_minimal) THEN
      CPErrorMessage(cp_failure_level,routineP,"VIRT TRUNCATION TO MINIMAL BASIS IS NIY",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (almo_scf_env%domain_layout_mos.NE.almo_domain_layout_molecular) THEN
      CPErrorMessage(cp_failure_level,routineP,"use MOLECULAR domains",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (almo_scf_env%domain_layout_aos.NE.almo_domain_layout_molecular) THEN
      CPErrorMessage(cp_failure_level,routineP,"use MOLECULAR domains",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (almo_scf_env%mat_distr_mos.NE.almo_mat_distr_molecular) THEN
      CPErrorMessage(cp_failure_level,routineP,"use MOLECULAR distr for MOs",error)
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (almo_scf_env%mat_distr_aos==almo_mat_distr_molecular .AND. &
        almo_scf_env%domain_layout_aos==almo_domain_layout_atomic) THEN 
       CPErrorMessage(cp_failure_level,routineP,"AO blocks cannot be larger than domains",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (almo_scf_env%mat_distr_mos==almo_mat_distr_molecular .AND. &
        almo_scf_env%domain_layout_mos==almo_domain_layout_atomic) THEN 
       CPErrorMessage(cp_failure_level,routineP,"MO blocks cannot be larger than domains",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (almo_scf_env%quencher_r1_factor.gt.almo_max_cutoff_multiplier) THEN
       CPErrorMessage(cp_failure_level,routineP,"DELOCALIZE_R_CUTOFF_FACTOR is larger than almo_max_cutoff_multiplier",error)
       CPErrorMessage(cp_failure_level,routineP,"increase hard-coded almo_max_cutoff_multiplier",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_scf_init_read_write_input

! *****************************************************************************
!> \brief Prints out a short summary about the ALMO SCF job
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_print_job_info(almo_scf_env,unit_nr,error)
    
    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_print_job_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'()')
       WRITE(unit_nr,'(T2,A,A,A)') REPEAT("-",32)," ALMO SETTINGS ",REPEAT("-",32)
    
       WRITE(unit_nr,'(T2,A,T48,E33.3)') "eps_filter:",almo_scf_env%eps_filter
       WRITE(unit_nr,'(T2,A,T48,E33.3)') "blocked_eps_iter:",almo_scf_env%eps_scf_bd
       WRITE(unit_nr,'(T2,A,T48,I33)') "blocked_max_iter:",almo_scf_env%max_scf_bd
       WRITE(unit_nr,'(T2,A,T48,I33)') "blocked_n_diis:",almo_scf_env%ndiis_bd
       
       SELECT CASE(almo_scf_env%deloc_method)
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
       
       IF (almo_scf_env%deloc_method.ne.almo_deloc_none) THEN
          
          SELECT CASE(almo_scf_env%almo_update_algorithm_q)
          CASE(almo_scf_diag)
              WRITE(unit_nr,'(T2,A,T48,A33)') "algorithm:","Domain diagonalization"
          CASE(almo_scf_pcg)
              WRITE(unit_nr,'(T2,A,T48,A33)') "algorithm:","PCG"
          END SELECT
          
          WRITE(unit_nr,'(T2,A,T48,F33.5)') "delocalize_r_cutoff_factor:",almo_scf_env%quencher_r0_factor
          WRITE(unit_nr,'(T2,A,T48,E33.3)') "delocalize_eps_iter:",almo_scf_env%eps_scf_q
          WRITE(unit_nr,'(T2,A,T48,I33)') "delocalize_max_iter:",almo_scf_env%max_scf_q
          
          IF (almo_scf_env%almo_update_algorithm_q.eq.almo_scf_pcg) THEN
             WRITE(unit_nr,'(T2,A,T48,E33.3)') "delocalize_eps_lin_search:",almo_scf_env%eps_lin_search
             SELECT CASE(almo_scf_env%scf_conjugator)
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
 
       !SELECT CASE(almo_scf_env%domain_layout_mos)
       !CASE(almo_domain_layout_orbital)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Delocalization constraints","ORBITAL"
       !CASE(almo_domain_layout_atomic)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Delocalization constraints","ATOMIC"
       !CASE(almo_domain_layout_molecular)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Delocalization constraints","MOLECULAR"
       !END SELECT

       !SELECT CASE(almo_scf_env%domain_layout_aos)
       !CASE(almo_domain_layout_atomic)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Basis function domains","ATOMIC"
       !CASE(almo_domain_layout_molecular)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Basis function domains","MOLECULAR"
       !END SELECT

       !SELECT CASE(almo_scf_env%mat_distr_aos)
       !CASE(almo_mat_distr_atomic)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Parallel distribution for AOs","ATOMIC"
       !CASE(almo_mat_distr_molecular)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Parallel distribution for AOs","MOLECULAR"
       !END SELECT

       !SELECT CASE(almo_scf_env%mat_distr_mos)
       !CASE(almo_mat_distr_atomic)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Parallel distribution for MOs","ATOMIC"
       !CASE(almo_mat_distr_molecular)
       !    WRITE(unit_nr,'(T2,A,T48,A33)') "Parallel distribution for MOs","MOLECULAR"
       !END SELECT

       !!! MOL INFO: WRITE(unit_nr,'(T2,A)') REPEAT("-",79)

       !!! MOL INFO: ! print fragment info
       !!! MOL INFO: WRITE(unit_nr,'(T2,A11,A17,A17,A17,A16)') "Fragment","Basis Set","Occupied MOs",&
       !!! MOL INFO:                                           "Active Virt","Discarded Virt"
       !!! MOL INFO: WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
       !!! MOL INFO: DO idomain=1,almo_scf_env%ndomains

       !!! MOL INFO:    WRITE(unit_nr,'(T2,I11,I17,I17,I17,I16)') &
       !!! MOL INFO:       idomain,almo_scf_env%nbasis_of_domain(idomain),&
       !!! MOL INFO:               SUM(almo_scf_env%nocc_of_domain(idomain,:)),&
       !!! MOL INFO:               SUM(almo_scf_env%nvirt_of_domain(idomain,:)),&
       !!! MOL INFO:               SUM(almo_scf_env%nvirt_disc_of_domain(idomain,:))

       !!! MOL INFO: ENDDO ! cycle over domains
       
       WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
       WRITE(unit_nr,'()')

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_scf_print_job_info

! *****************************************************************************
!> \brief Initializes the ALMO SCF copy of the AO overlap matrix 
!>        and all necessary functions (sqrt, inverse...)
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_init_ao_overlap(matrix_s,almo_scf_env,error)
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s
    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_init_ao_overlap', &
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

    ! make almo copy of S
    ! also copy S to S_blk (i.e. to S with the domain structure imposed)
    IF (almo_scf_env%orthogonal_basis) THEN
       CALL cp_dbcsr_set(almo_scf_env%matrix_s(1),0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(almo_scf_env%matrix_s(1),1.0_dp,error=error)
       CALL cp_dbcsr_set(almo_scf_env%matrix_s_blk(1),0.0_dp,error=error)
       CALL cp_dbcsr_add_on_diag(almo_scf_env%matrix_s_blk(1),1.0_dp,error=error)
    ELSE
       CALL matrix_qs_to_almo(matrix_s,almo_scf_env%matrix_s(1),&
                              almo_scf_env,.FALSE.,error=error)
       CALL matrix_qs_to_almo(matrix_s,almo_scf_env%matrix_s_blk(1),&
                              almo_scf_env,.TRUE.,error=error)
    ENDIF
    
    CALL cp_dbcsr_filter(almo_scf_env%matrix_s(1),almo_scf_env%eps_filter,error=error)
    CALL cp_dbcsr_filter(almo_scf_env%matrix_s_blk(1),almo_scf_env%eps_filter,error=error)

    IF (almo_scf_env%almo_update_algorithm.eq.almo_scf_diag) THEN
       CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_s_blk_sqrt(1),&
               almo_scf_env%matrix_s_blk_sqrt_inv(1),&
               almo_scf_env%matrix_s_blk(1),&
               threshold=almo_scf_env%eps_filter,&
               order=almo_scf_env%order_lanczos,&
               eps_lanczos=almo_scf_env%eps_lanczos,&
               max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
               error=error)
    ELSE IF (almo_scf_env%almo_update_algorithm.eq.almo_scf_dm_sign) THEN
       CALL invert_Hotelling(almo_scf_env%matrix_s_blk_inv(1),&
               almo_scf_env%matrix_s_blk(1),&
               threshold=almo_scf_env%eps_filter,&
               error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_scf_init_ao_overlap

! *****************************************************************************
!> \brief Selects the subroutine for the optimization of block-daigonal ALMOs.
!>        Keep it short and clean.
!> \par History
!>       2011.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_main(qs_env,almo_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_main', &
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

    IF (almo_scf_env%almo_update_algorithm.eq.almo_scf_pcg) THEN
       ! try this optimizer if the DIIS optimizer fails
       CALL almo_scf_xalmo_pcg(qs_env,almo_scf_env,&
               quench_t=almo_scf_env%quench_t_blk,&
               matrix_t_in=almo_scf_env%matrix_t_blk,&
               matrix_t_out=almo_scf_env%matrix_t_blk,&
               assume_t0_q0x=.FALSE.,perturbation_only=.FALSE.,&
               special_case=0,error=error)
       CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env,error)
    ELSE
       ! mixing/DIIS optimizer
       CALL almo_scf_block_diagonal(qs_env,almo_scf_env,error)
    ENDIF
    
    ! we might need a copy of the converged KS
    DO ispin=1,almo_scf_env%nspins
       CALL cp_dbcsr_copy(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                 almo_scf_env%matrix_ks(ispin),error=error)
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE almo_scf_main

! *****************************************************************************
!> \brief An SCF procedure that optimizes block-diagonal ALMOs using DIIS
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_block_diagonal(qs_env,almo_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_block_diagonal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iscf, ispin, nspin, &
                                                unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: local_nocc_of_domain
    LOGICAL                                  :: converged, failure, &
                                                prepare_to_exit, should_stop, &
                                                use_diis, use_prev_as_guess
    REAL(KIND=dp) :: energy_diff, energy_new, energy_old, error_norm, &
      error_norm_ispin, t1, t2, true_mixing_fraction
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: local_mu
    TYPE(almo_scf_diis_type), ALLOCATABLE, &
      DIMENSION(:)                           :: almo_diis
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: matrix_mixing_old_blk
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! use DIIS, it's superior to simple mixing
    use_diis=.TRUE.
    use_prev_as_guess=.FALSE.

    nspin=almo_scf_env%nspins
    ALLOCATE(local_mu(almo_scf_env%ndomains))
    ALLOCATE(local_nocc_of_domain(almo_scf_env%ndomains))

    ! init mixing matrices
    ALLOCATE(matrix_mixing_old_blk(nspin))
    ALLOCATE(almo_diis(nspin))
    DO ispin=1,nspin
       CALL cp_dbcsr_init(matrix_mixing_old_blk(ispin),error=error)
       CALL cp_dbcsr_create(matrix_mixing_old_blk(ispin),&
              template=almo_scf_env%matrix_ks_blk(ispin),error=error)
       CALL almo_scf_diis_init(diis_env=almo_diis(ispin),&
              sample_err=almo_scf_env%matrix_ks_blk(ispin),&
              sample_var=almo_scf_env%matrix_s_blk(1),&
              error_type=1,&
              max_length=almo_scf_env%ndiis_bd,&
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

    ! the real SCF loop
    t1 = m_walltime()
    DO 

      iscf=iscf+1

      ! get a copy of the current KS matrix
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
      DO ispin=1,nspin
         CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix,&
                 almo_scf_env%matrix_ks(ispin),&
                 almo_scf_env,.FALSE.,error=error)
         CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix,&
                 almo_scf_env%matrix_ks_blk(ispin),&
                 almo_scf_env,.TRUE.,error=error)
         CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin),&
                  almo_scf_env%eps_filter,error=error)
      ENDDO

      ! obtain projected KS matrix and the DIIS-error vector
      CALL almo_scf_ks_to_ks_blk(almo_scf_env,error=error)

      ! inform the DIIS handler about the new KS matrix and its error vector
      IF (use_diis) THEN
         DO ispin=1,nspin
            CALL almo_scf_diis_push(diis_env=almo_diis(ispin),&
                    var=almo_scf_env%matrix_ks_blk(ispin),&
                    err=almo_scf_env%matrix_err_blk(ispin),&
                    error=error)
         ENDDO
      ENDIF
 
      ! get error_norm: choose the largest of the two spins
      DO ispin=1,nspin
         !error_norm=cp_dbcsr_frobenius_norm(almo_scf_env%matrix_err_blk(ispin))
         CALL cp_dbcsr_norm(almo_scf_env%matrix_err_blk(ispin),&
                 dbcsr_norm_maxabsnorm,&
                 norm_scalar=error_norm_ispin, error=error)
         IF (ispin.eq.1) error_norm=error_norm_ispin
         IF (ispin.gt.1 .AND. error_norm_ispin.gt.error_norm) &
            error_norm=error_norm_ispin
      ENDDO
      
      IF (error_norm.lt.almo_scf_env%eps_prev_guess) THEN
         use_prev_as_guess = .TRUE.
      ELSE
         use_prev_as_guess = .FALSE.
      ENDIF 
      
      ! check convergence
      converged=.TRUE.
      IF (error_norm.gt.almo_scf_env%eps_scf_bd) converged=.FALSE.
      ! check other exit criteria: max SCF steps and timing
      CALL external_control(should_stop,"SCF",&
              start_time=qs_env%start_time,&
              target_time=qs_env%target_time,error=error)
      IF (should_stop .OR. iscf>=almo_scf_env%max_scf_bd .OR. converged)  THEN
         prepare_to_exit=.TRUE.
      ENDIF
      
      IF (.NOT.prepare_to_exit) THEN ! update the ALMOs and density matrix
         
         ! perform mixing of KS matrices
         IF (iscf.ne.1) THEN
            IF (use_diis) THEN ! use diis instead of mixing
               DO ispin=1,nspin
                  CALL almo_scf_diis_extrapolate(diis_env=almo_diis(ispin),&
                          extr_var=almo_scf_env%matrix_ks_blk(ispin),&
                          error=error)
               ENDDO
            ELSE ! use mixing
               true_mixing_fraction=almo_scf_env%mixing_fraction
               DO ispin=1,nspin
                  CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),&
                                    matrix_mixing_old_blk(ispin),& 
                                    true_mixing_fraction,&
                                    1.0_dp-true_mixing_fraction,&
                                    error=error)
               END DO
            ENDIF
         ENDIF
         ! save the new matrix for the future mixing
         DO ispin=1,nspin
            CALL cp_dbcsr_copy(matrix_mixing_old_blk(ispin),&
                   almo_scf_env%matrix_ks_blk(ispin), error=error)
         ENDDO
   
         ! obtain ALMOs from the new KS matrix
         SELECT CASE (almo_scf_env%almo_update_algorithm) 
         CASE (almo_scf_diag) 
   
            CALL almo_scf_ks_blk_to_tv_blk(almo_scf_env,error) 
   
         CASE (almo_scf_dm_sign)
   
            ! update the density matrix
            DO ispin=1,nspin
      
               local_nocc_of_domain(:)=almo_scf_env%nocc_of_domain(:,ispin)
               local_mu(:)=almo_scf_env%mu_of_domain(:,ispin)
               ! RZK UPDATE! the update algorithm is removed because
               ! RZK UPDATE! it requires updating core LS_SCF routines
               ! RZK UPDATE! (the code exists in the CVS version)
               IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "CVS only: density_matrix_sign has not been updated in SVN"
               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
               ! RZK UPDATE!  CALL density_matrix_sign(almo_scf_env%matrix_p_blk(ispin),&
               ! RZK UPDATE!          local_mu,&
               ! RZK UPDATE!          almo_scf_env%fixed_mu,&
               ! RZK UPDATE!          almo_scf_env%matrix_ks_blk(ispin),&
               ! RZK UPDATE!          !matrix_mixing_old_blk(ispin),&
               ! RZK UPDATE!          almo_scf_env%matrix_s_blk(1), &
               ! RZK UPDATE!          almo_scf_env%matrix_s_blk_inv(1), &
               ! RZK UPDATE!          local_nocc_of_domain,&
               ! RZK UPDATE!          almo_scf_env%eps_filter,&
               ! RZK UPDATE!          almo_scf_env%domain_index_of_ao,&
               ! RZK UPDATE!          error)
               almo_scf_env%mu_of_domain(:,ispin)=local_mu(:)
      
            ENDDO
      
            ! obtain ALMOs from matrix_p_blk: T_new = P_blk S_blk T_old
            CALL almo_scf_p_blk_to_t_blk(almo_scf_env,error)
            CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env,error=error)
   
         END SELECT

         ! obtain density matrix from ALMOs
         CALL almo_scf_t_blk_to_p(almo_scf_env,&
                 use_sigma_inv_guess=use_prev_as_guess,error=error)

         ! compute the new KS matrix and new energy
         CALL almo_scf_dm_to_ks(qs_env,almo_scf_env,energy_new,error)
   
      ENDIF ! prepare_to_exit

      energy_diff=energy_new-energy_old
      energy_old=energy_new
      almo_scf_env%almo_scf_energy=energy_new

      t2 = m_walltime()
      ! brief report on the current SCF loop
      IF (unit_nr>0) THEN
         WRITE(unit_nr,'(T2,A13,I6,F23.10,E14.5,F14.9,F9.2)') "ALMO SCF DIIS",&
               iscf,&
               energy_new,energy_diff,error_norm, t2-t1
      ENDIF
      t1 = m_walltime()
 
      IF (prepare_to_exit) EXIT

    ENDDO ! end scf cycle

    IF (.NOT.converged)  THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "SCF for block-diagonal ALMOs not converged! "
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    DO ispin=1,nspin
       CALL cp_dbcsr_release(matrix_mixing_old_blk(ispin),error=error)
       CALL almo_scf_diis_release(diis_env=almo_diis(ispin),error=error)
    ENDDO
    DEALLOCATE(almo_diis)
    DEALLOCATE(matrix_mixing_old_blk)
    DEALLOCATE(local_mu)
    DEALLOCATE(local_nocc_of_domain)
 
    CALL timestop(handle)

  END SUBROUTINE almo_scf_block_diagonal

! *****************************************************************************
!> \brief An eigensolver-based SCF to optimize extended ALMOs (i.e. ALMOs on
!>        overlapping domains)
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_xalmo_eigensolver(qs_env,almo_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_xalmo_eigensolver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iscf, ispin, nspin, &
                                                unit_nr
    LOGICAL                                  :: converged, failure, &
                                                prepare_to_exit, should_stop
    REAL(KIND=dp) :: denergy_tot, energy_diff, energy_new, energy_old, &
      error_norm, error_norm_0, spin_factor, t1, t2
    REAL(KIND=dp), DIMENSION(2)              :: denergy_spin
    TYPE(almo_scf_diis_type), ALLOCATABLE, &
      DIMENSION(:)                           :: almo_diis
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type)                      :: matrix_p_almo_scf_converged
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:, :)           :: submatrix_mixing_old_blk

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=almo_scf_env%nspins
    IF (nspin == 1) THEN
     spin_factor = 2.0_dp
    ELSE
     spin_factor = 1.0_dp
    ENDIF
   
    ! RZK-warning domain_s_sqrt and domain_s_sqrt_inv do not have spin
    ! components yet (may be used later)
    ispin=1 
    CALL construct_domain_s_sqrt(&
            matrix_s=almo_scf_env%matrix_s(1),&
            subm_s_sqrt=almo_scf_env%domain_s_sqrt(:,ispin),&
            subm_s_sqrt_inv=almo_scf_env%domain_s_sqrt_inv(:,ispin),&
            dpattern=almo_scf_env%quench_t(ispin),&
            map=almo_scf_env%domain_map(ispin),&
            node_of_domain=almo_scf_env%cpu_of_domain,&
            error=error)
    ! TRY: construct s_inv
    !CALL construct_domain_s_inv(&
    !       matrix_s=almo_scf_env%matrix_s(1),&
    !       subm_s_inv=almo_scf_env%domain_s_inv(:,ispin),&
    !       dpattern=almo_scf_env%quench_t(ispin),&
    !       map=almo_scf_env%domain_map(ispin),&
    !       node_of_domain=almo_scf_env%cpu_of_domain,&
    !       error=error)

    ! construct the domain template for the occupied orbitals
    DO ispin=1,nspin
       ! RZK-warning we need only the matrix structure, not data
       ! replace construct_submatrices with lighter procedure with
       ! no heavy communications
       CALL construct_submatrices(&
               matrix=almo_scf_env%quench_t(ispin),&
               submatrix=almo_scf_env%domain_t(:,ispin),&
               distr_pattern=almo_scf_env%quench_t(ispin),&
               domain_map=almo_scf_env%domain_map(ispin),&
               node_of_domain=almo_scf_env%cpu_of_domain,&
               job_type=select_row,error=error)
    ENDDO

    ! init mixing matrices
    ALLOCATE(submatrix_mixing_old_blk(almo_scf_env%ndomains,nspin))
    CALL init_submatrices(submatrix_mixing_old_blk,error=error)
    ALLOCATE(almo_diis(nspin))
    
    ! TRY: construct block-projector
    !ALLOCATE(submatrix_tmp(almo_scf_env%ndomains))
    !DO ispin=1,nspin
    !   CALL init_submatrices(submatrix_tmp,error=error)
    !   CALL construct_domain_r_down(&
    !           matrix_t=almo_scf_env%matrix_t_blk(ispin),&
    !           matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
    !           matrix_s=almo_scf_env%matrix_s(1),&
    !           subm_r_down=submatrix_tmp(:),&
    !           dpattern=almo_scf_env%quench_t(ispin),&
    !           map=almo_scf_env%domain_map(ispin),&
    !           node_of_domain=almo_scf_env%cpu_of_domain,&
    !           filter_eps=almo_scf_env%eps_filter,&
    !           error=error)
    !   CALL multiply_submatrices('N','N',1.0_dp,&
    !           submatrix_tmp(:),&
    !           almo_scf_env%domain_s_inv(:,1),0.0_dp,&
    !           almo_scf_env%domain_r_down_up(:,ispin),&
    !           error)
    !   CALL release_submatrices(submatrix_tmp,error=error)
    !ENDDO
    !DEALLOCATE(submatrix_tmp)

    DO ispin=1,nspin
       ! use s_sqrt since they are already properly constructed
       ! and have the same distributions as domain_err and domain_ks_xx
       CALL almo_scf_diis_init(diis_env=almo_diis(ispin),&
              sample_err=almo_scf_env%domain_s_sqrt(:,ispin),&
              sample_var=almo_scf_env%domain_s_sqrt_inv(:,ispin),&
              error_type=1,&
              max_length=almo_scf_env%ndiis_q,&
              error=error)
    ENDDO

    denergy_tot=0.0_dp
    energy_old=0.0_dp
    iscf=0
    prepare_to_exit=.FALSE.

    ! the SCF loop
    t1 = m_walltime()
    DO 

      iscf=iscf+1

      ! get a copy of the current KS matrix
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
      DO ispin=1,nspin
         CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix,&
                 almo_scf_env%matrix_ks(ispin),&
                 almo_scf_env,.FALSE.,error=error)
         CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix,&
                 almo_scf_env%matrix_ks_blk(ispin),&
                 almo_scf_env,.TRUE.,error=error)
         CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin),&
                  almo_scf_env%eps_filter,error=error)
      ENDDO

      ! obtain projected KS matrix and the DIIS-error vector
      CALL almo_scf_ks_to_ks_xx(almo_scf_env,error=error)
      
      ! inform the DIIS handler about the new KS matrix and its error vector
      DO ispin=1,nspin
         CALL almo_scf_diis_push(diis_env=almo_diis(ispin),&
                 d_var=almo_scf_env%domain_ks_xx(:,ispin),&
                 d_err=almo_scf_env%domain_err(:,ispin),&
                 error=error)
      ENDDO
 
      ! check convergence
      converged=.TRUE.
      DO ispin=1,nspin
         !error_norm=cp_dbcsr_frobenius_norm(almo_scf_env%matrix_err_blk(ispin))
         CALL cp_dbcsr_norm(almo_scf_env%matrix_err_xx(ispin),&
                 dbcsr_norm_maxabsnorm,&
                 norm_scalar=error_norm, error=error)
         CALL maxnorm_submatrices(almo_scf_env%domain_err(:,ispin),&
                 norm=error_norm_0,&
                 error=error)
         IF (error_norm.gt.almo_scf_env%eps_scf_q) THEN
            converged=.FALSE.
            EXIT ! no need to check the other spin
         ENDIF
      ENDDO
      ! check other exit criteria: max SCF steps and timing
      CALL external_control(should_stop,"SCF",&
              start_time=qs_env%start_time,&
              target_time=qs_env%target_time,error=error)
      IF (should_stop .OR. iscf>=almo_scf_env%max_scf_q .OR. converged)  THEN
         prepare_to_exit=.TRUE.
      ENDIF

      IF (.NOT.prepare_to_exit) THEN ! update the ALMOs and density matrix
         
         ! perform mixing of KS matrices
         IF (iscf.ne.1) THEN
            IF (.FALSE.) THEN ! use diis instead of mixing
               DO ispin=1,nspin
                  CALL add_submatrices(&
                          almo_scf_env%mixing_fraction,&
                          almo_scf_env%domain_ks_xx(:,ispin),&
                          1.0_dp-almo_scf_env%mixing_fraction,&
                          submatrix_mixing_old_blk(:,ispin),&
                          'N',error)
                  !CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),&
                  !        matrix_mixing_old_blk(ispin),& 
                  !        almo_scf_env%mixing_fraction,&
                  !        1.0_dp-almo_scf_env%mixing_fraction,&
                  !        error=error)
               END DO
            ELSE
               DO ispin=1,nspin
                  CALL almo_scf_diis_extrapolate(diis_env=almo_diis(ispin),&
                          d_extr_var=almo_scf_env%domain_ks_xx(:,ispin),&
                          error=error)
               ENDDO
            ENDIF
         ENDIF
         ! save the new matrix for the future mixing
         DO ispin=1,nspin
            CALL copy_submatrices(&
                   almo_scf_env%domain_ks_xx(:,ispin),&
                   submatrix_mixing_old_blk(:,ispin),&
                   copy_data=.TRUE., error=error)
         ENDDO
   
         ! obtain a new set of ALMOs from the updated KS matrix
         CALL almo_scf_ks_xx_to_tv_xx(almo_scf_env,error) 

         ! update the density matrix
         DO ispin=1,nspin
            
            ! save the initial density matrix (to get the perturbative energy lowering)
            IF (iscf.eq.1) THEN
               CALL cp_dbcsr_init(matrix_p_almo_scf_converged, error=error)
               CALL cp_dbcsr_create(matrix_p_almo_scf_converged,&
                       template=almo_scf_env%matrix_p(ispin),&
                       error=error)
               CALL cp_dbcsr_copy(matrix_p_almo_scf_converged,&
                       almo_scf_env%matrix_p(ispin),error=error)
            ENDIF

            ! update now 
            CALL almo_scf_t_to_p(&
                    t=almo_scf_env%matrix_t(ispin),&
                    p=almo_scf_env%matrix_p(ispin),&
                    eps_filter=almo_scf_env%eps_filter,&
                    orthog_orbs=.FALSE.,&
                    s=almo_scf_env%matrix_s(1),&
                    sigma=almo_scf_env%matrix_sigma(ispin),&
                    sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
                    use_guess=.TRUE.,&
                    error=error)
            CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin),spin_factor,&
                    error=error)

            ! obtain perturbative estimate (at no additional cost) 
            ! of the energy lowering relative to the block-diagonal ALMOs 
            IF (iscf.eq.1) THEN
              
               CALL cp_dbcsr_add(matrix_p_almo_scf_converged,&
                                almo_scf_env%matrix_p(ispin),-1.0_dp,1.0_dp,&
                                error=error)
               CALL cp_dbcsr_trace(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                                   matrix_p_almo_scf_converged,&
                                   denergy_spin(ispin),error=error)
            
               CALL cp_dbcsr_release(matrix_p_almo_scf_converged, error=error)
               
               denergy_tot=denergy_tot+denergy_spin(ispin)
  
               ! RZK-warning Energy correction can be evaluated using matrix_x
               ! as shown in the attempt below and in the PCG procedure.
               ! Using matrix_x allows immediate decomposition of the energy
               ! lowering into 2-body components for EDA. However, it does not
               ! work here because the diagonalization routine does not necessarily 
               ! produce orbitals with the same sign as the block-diagonal ALMOs
               ! Any fixes?!
 
               !CALL cp_dbcsr_init(matrix_x, error=error)
               !CALL cp_dbcsr_create(matrix_x,&
               !        template=almo_scf_env%matrix_t(ispin),&
               !        error=error)
               !
               !CALL cp_dbcsr_init(matrix_tmp_no, error=error)
               !CALL cp_dbcsr_create(matrix_tmp_no,&
               !        template=almo_scf_env%matrix_t(ispin),&
               !        error=error)
               !
               !CALL cp_dbcsr_copy(matrix_x,&
               !        almo_scf_env%matrix_t_blk(ispin),error=error)
               !CALL cp_dbcsr_add(matrix_x,almo_scf_env%matrix_t(ispin),&
               !        -1.0_dp,1.0_dp,&
               !        error=error)

               !CALL cp_dbcsr_trace(matrix_x,&
               !        almo_scf_env%matrix_err_xx(ispin),denergy,"T","N",&
               !        error=error)

               !denergy=denergy*spin_factor

               !IF (unit_nr>0) THEN
               !   WRITE(unit_nr,*) "_ENERGY-0: ", almo_scf_env%almo_scf_energy
               !   WRITE(unit_nr,*) "_ENERGY-D: ", denergy
               !   WRITE(unit_nr,*) "_ENERGY-F: ", almo_scf_env%almo_scf_energy+denergy
               !ENDIF
               !! RZK-warning update will not work since the energy is overwritten almost immediately
               !!CALL almo_scf_update_ks_energy(qs_env,&
               !!        almo_scf_env%almo_scf_energy+denergy,&
               !!        error)

               !! print out the results of the decomposition analysis
               !CALL cp_dbcsr_hadamard_product(matrix_x,&
               !        almo_scf_env%matrix_err_xx(ispin),&
               !        matrix_tmp_no,error=error)
               !CALL cp_dbcsr_scale(matrix_tmp_no,spin_factor,error=error)
               !CALL cp_dbcsr_filter(matrix_tmp_no,almo_scf_env%eps_filter,error=error)
               !
               !IF (unit_nr>0) THEN
               !   WRITE(unit_nr,*)
               !   WRITE(unit_nr,'(T2,A)') "DECOMPOSITION OF THE DELOCALIZATION ENERGY"
               !ENDIF

               !mynode=dbcsr_mp_mynode(dbcsr_distribution_mp(&
               !   cp_dbcsr_distribution(matrix_tmp_no)))
               !WRITE(mynodestr,'(I6.6)') mynode
               !mylogfile='EDA.'//TRIM(ADJUSTL(mynodestr))
               !OPEN (iunit,file=mylogfile,status='REPLACE')
               !CALL cp_dbcsr_print_block_sum(matrix_tmp_no,iunit,error=error)
               !CLOSE(iunit)
               !
               !CALL cp_dbcsr_release(matrix_tmp_no, error=error)
               !CALL cp_dbcsr_release(matrix_x, error=error)

            ENDIF ! iscf.eq.1

         ENDDO
         
         ! print out the energy lowering
         IF (iscf.eq.1) THEN
            IF (unit_nr>0) THEN
               WRITE(unit_nr,*)
               WRITE(unit_nr,'(T2,A35,F25.10)') "ENERGY OF BLOCK-DIAGONAL ALMOs:",&
                  almo_scf_env%almo_scf_energy
               WRITE(unit_nr,'(T2,A35,F25.10)') "ENERGY LOWERING:",&
                  denergy_tot
               WRITE(unit_nr,'(T2,A35,F25.10)') "CORRECTED ENERGY:",&
                  almo_scf_env%almo_scf_energy+denergy_tot
               WRITE(unit_nr,*)
            ENDIF
            CALL almo_scf_update_ks_energy(qs_env,&
                    almo_scf_env%almo_scf_energy+denergy_tot,&
                    error)
         ENDIF

         ! compute the new KS matrix and new energy
         IF (.NOT.almo_scf_env%perturbative_delocalization) THEN
            CALL almo_scf_dm_to_ks(qs_env,almo_scf_env,energy_new,error)
         ENDIF
   
      ENDIF ! prepare_to_exit

      IF (almo_scf_env%perturbative_delocalization) THEN
         
         ! exit after the first step if we do not need the SCF procedure
         converged=.TRUE.
         prepare_to_exit=.TRUE.

      ELSE ! not a perturbative treatment

         energy_diff=energy_new-energy_old
         energy_old=energy_new
         almo_scf_env%almo_scf_energy=energy_new

         t2 = m_walltime()
         ! brief report on the current SCF loop
         IF (unit_nr>0) THEN
            WRITE(unit_nr,'(T2,A,I6,F20.9,E11.3,E11.3,E11.3,F8.2)') "ALMO SCF",&
                  iscf,&
                  energy_new,energy_diff,error_norm,error_norm_0, t2-t1
         ENDIF
         t1 = m_walltime()

      ENDIF
 
      IF (prepare_to_exit) EXIT

    ENDDO ! end scf cycle

    IF (.NOT.converged)  THEN
      IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "SCF for ALMOs on overlapping domains not converged! "
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    DO ispin=1,nspin
       CALL release_submatrices(submatrix_mixing_old_blk(:,ispin),error=error)
       CALL almo_scf_diis_release(diis_env=almo_diis(ispin),error=error)
    ENDDO
    DEALLOCATE(almo_diis)
    DEALLOCATE(submatrix_mixing_old_blk)
 
    CALL timestop(handle)

  END SUBROUTINE almo_scf_xalmo_eigensolver

! *****************************************************************************
!> \brief Optimization of ALMOs using PCG-like minimizers
!> \param special_case - to reduce the overhead special cases are implemented 
!>                      -1. - no special case (i.e. xALMOs)
!>                       0. - block-diagonal ALMOs
!>                       1. - fully delocalized ALMOs
!> \par History
!>       2011.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_xalmo_pcg(qs_env,almo_scf_env,quench_t,&
                matrix_t_in,matrix_t_out,assume_t0_q0x,perturbation_only,&
                special_case,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: quench_t, matrix_t_in, &
                                                matrix_t_out
    LOGICAL, INTENT(IN)                      :: assume_t0_q0x, &
                                                perturbation_only
    INTEGER, INTENT(IN), OPTIONAL            :: special_case
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_xalmo_pcg', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: iunit = 78

    CHARACTER(LEN=100)                       :: mylogfile, mynodestr
    CHARACTER(LEN=20)                        :: iter_type
    INTEGER :: cg_iteration, dim0, fixed_line_search_niter, handle, ispin, &
      iteration, jj, line_search_iteration, max_iter, my_special_case, &
      mynode, ncores, ndomains, occ1, outer_iteration, outer_max_iter, &
      prec_type, precond_domain_projector, unit_nr, zero_neg_eiv
    LOGICAL :: converged, do_md, failure, first_md_iteration, just_started, &
      line_search, md_in_theta_space, optimize_theta, outer_prepare_to_exit, &
      prepare_to_exit, reset_conjugator, skip_grad, use_guess, &
      use_preconditioner
    REAL(kind=dp) :: appr_sec_der, beta, denom, e0, e1, energy_diff, &
      energy_new, energy_old, eps_skip_gradients, g0, g1, grad_norm, &
      grad_norm_frob, kappa, kin_energy, line_search_error, &
      next_step_size_guess, numer, prec_sf_mixing_s, shift_value, &
      spin_factor, step_size, t1, t2, t_norm, tau, time_step
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(cp_dbcsr_type) :: FTsiginv, fvo_0, grad, inv_eiv, m_theta, &
      m_tmp_nn_1, m_tmp_no_1, m_tmp_no_2, m_tmp_no_3, m_tmp_oo_1, matrix_p_0, &
      matrix_sigma_0, matrix_sigma_inv_0, matrix_t_0, prec_oo, prec_oo_inv, &
      prec_vv, prev_grad, prev_minus_prec_grad, prev_step, siginvTFTsiginv, &
      ST, step, STsiginv_0, velocity
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)              :: domain_r_down

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

    ! set local flags
    do_md=almo_scf_env%logical01
    !perturbation_only=almo_scf_env%logical02
    !use_projector=almo_scf_env%logical03
    use_preconditioner=almo_scf_env%logical04
    optimize_theta=almo_scf_env%logical05
    ! set local parameters
    prec_sf_mixing_s=almo_scf_env%real04
    prec_type=almo_scf_env%integer02
    eps_skip_gradients=almo_scf_env%real01
    fixed_line_search_niter=0 ! init to zero, change when eps is small enough
    !fixed_line_search_niter=abs(almo_scf_env%integer04)! init to zero, change when eps is small enough
    
    ncores = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
           cp_dbcsr_distribution(almo_scf_env%matrix_s(1))))
    
    ! Since it is extrememly difficult to converge the iterative
    ! procedure using T as an optimized variable, assume 
    ! T = T_0 + (1-R_0)*X and optimize X for both perturbative 
    ! and non-perturbative cases

    IF (almo_scf_env%nspins == 1) THEN
       spin_factor = 2.0_dp
    ELSE
       spin_factor = 1.0_dp
    ENDIF
    
    !!!!!! RZK-warning THIS PROCEDURE WILL WORK ONLY FOR CLOSED SHELL SYSTEMS
    !!!!!! TO ADAPT IT FOR UNRESTRICTED ORBITALS - UPDATE KS MATRIX WITH PARTIALLY
    !!!!!! OPTIMIZED ORBITALS - BOTH ALPNA AND BETA
    IF (almo_scf_env%nspins.gt.1) THEN
       CPErrorMessage(cp_failure_level,routineP,"UNRESTRICTED ALMO SCF IS NYI(!)",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    DO ispin=1,almo_scf_env%nspins

       ! init temporary storage
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
               template=almo_scf_env%matrix_ks(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(prec_oo,&
               template=almo_scf_env%matrix_sigma(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(prec_oo_inv,&
               template=almo_scf_env%matrix_sigma(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(m_tmp_oo_1,&
               template=almo_scf_env%matrix_sigma(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(siginvTFTsiginv,&
               template=almo_scf_env%matrix_sigma(ispin),&
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
       
       !ndomains = cp_dbcsr_nblkcols_total(quench_t(ispin))
       ndomains = almo_scf_env%ndomains 
       ALLOCATE(domain_r_down(ndomains))
       CALL init_submatrices(domain_r_down,error)

       ! create matrices to store the initial state
       CALL cp_dbcsr_init(matrix_t_0,error=error)
       CALL cp_dbcsr_init(matrix_sigma_inv_0,error=error)
       CALL cp_dbcsr_init(matrix_sigma_0,error=error)
       CALL cp_dbcsr_init(matrix_p_0,error=error)
       CALL cp_dbcsr_create(matrix_t_0,&
               template=matrix_t_out(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(matrix_sigma_inv_0,&
               template=almo_scf_env%matrix_sigma_inv(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(matrix_sigma_0,&
               template=almo_scf_env%matrix_sigma_inv(ispin),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_create(matrix_p_0,&
               template=almo_scf_env%matrix_ks(ispin),&
               error=error)
       CALL cp_dbcsr_copy(matrix_t_0,matrix_t_in(ispin),&
               error=error)

       CALL cp_dbcsr_set(step,0.0_dp,error=error)

       md_in_theta_space=.FALSE. ! turn on later after several minimization steps
       IF (do_md) THEN
          CALL cp_dbcsr_init(velocity,error=error)
          CALL cp_dbcsr_create(velocity,&
                  template=matrix_t_out(ispin),error=error)
          CALL cp_dbcsr_copy(velocity,quench_t(ispin),error=error)
          CALL cp_dbcsr_set(velocity,0.0_dp,error=error)
          CALL cp_dbcsr_copy(prev_step,quench_t(ispin),error=error)
          CALL cp_dbcsr_set(prev_step,0.0_dp,error=error)
          time_step=almo_scf_env%real05
       ENDIF
       
       ! create initial guess from the initial orbitals
       IF (assume_t0_q0x) THEN
          CALL cp_dbcsr_set(m_theta,0.0_dp,error=error)
       ELSE
          IF (optimize_theta) THEN
             ! check that all MO coefficients of the guess are less
             ! than the maximum allowed amplitude
             CALL cp_dbcsr_norm(matrix_t_0,&
                     dbcsr_norm_maxabsnorm, norm_scalar=grad_norm, error=error)
             !IF (unit_nr>0) THEN
             !   WRITE(unit_nr,*) "Maximum norm of the initial guess: ", grad_norm
             !   WRITE(unit_nr,*) "Maximum allowed amplitude: ", almo_scf_env%envelope_amplitude
             !ENDIF
             IF (grad_norm.gt.almo_scf_env%envelope_amplitude) THEN
                CPErrorMessage(cp_failure_level,routineP,"Max norm of the initial guess is too large",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF
             ! use artanh of block-diagonal ALMOs as an initial guess for independent variables
             CALL cp_dbcsr_copy(m_theta,matrix_t_0,error=error)
             CALL cp_dbcsr_function_of_elements(m_theta,&
                     !func=dbcsr_func_asin,&
                     func=dbcsr_func_artanh,&
                     a0=0.0_dp,&
                     a1=1.0_dp/almo_scf_env%envelope_amplitude,&
                     error=error)
             CALL cp_dbcsr_scale(m_theta,almo_scf_env%envelope_amplitude,&
                     error=error)
          ELSE
             ! simply copy MO coefficients to m_theta
             CALL cp_dbcsr_copy(m_theta,matrix_t_0,error=error)
             CALL cp_dbcsr_norm(m_theta,&
                     dbcsr_norm_maxabsnorm, norm_scalar=grad_norm, error=error)
             !IF (unit_nr>0) THEN
             !   WRITE(unit_nr,*) "Maximum norm of the initial guess: ", grad_norm
             !ENDIF
          ENDIF
       ENDIF

       ! invert S domains if necessary
       ! RZK-warning must be done outside the spin loop to save time
       IF (my_special_case.eq.-1) THEN
          CALL construct_domain_s_inv(&
                 matrix_s=almo_scf_env%matrix_s(1),&
                 subm_s_inv=almo_scf_env%domain_s_inv(:,ispin),&
                 dpattern=quench_t(ispin),&
                 map=almo_scf_env%domain_map(ispin),&
                 node_of_domain=almo_scf_env%cpu_of_domain,&
                 error=error)
       ENDIF

       ! start the outer SCF loop
       outer_max_iter=almo_scf_env%outer_max_scf_q
       outer_prepare_to_exit=.FALSE.
       outer_iteration=0
       grad_norm=0.0_dp
       grad_norm_frob=0.0_dp
       use_guess=.FALSE.

       DO

          ! start the inner SCF loop
          max_iter=almo_scf_env%max_scf_q
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

             ! switch to MD after several minimization steps
             IF (iteration.eq.almo_scf_env%integer01.AND.do_md) THEN
                CALL cp_dbcsr_set(velocity,0.0_dp,error=error)
                CALL cp_dbcsr_set(prev_step,0.0_dp,error=error)
                md_in_theta_space=.TRUE.
                first_md_iteration=.TRUE.
             ENDIF

             ! compute the MO coefficients from theta
             IF (assume_t0_q0x.AND.just_started) THEN
                   CALL cp_dbcsr_set(matrix_t_out(ispin),0.0_dp,error=error)
             ENDIF
             IF (optimize_theta) THEN
                CALL cp_dbcsr_norm(m_theta,&
                        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                !IF (unit_nr>0) THEN
                !   WRITE(unit_nr,*) "VAR1: |T|: ", t_norm
                !   !WRITE(unit_nr,*) "VAR1: SIN(|T|/A), COS(|T|/A), A*SIN(|T|/A): ",&
                !   !   SIN(t_norm/almo_scf_env%envelope_amplitude),&
                !   !   COS(t_norm/almo_scf_env%envelope_amplitude),&
                !   !   almo_scf_env%envelope_amplitude*SIN(t_norm/almo_scf_env%envelope_amplitude)
                !   WRITE(unit_nr,*) "VAR1: TANH(|T|/A), D(A*TANH(|T|/A)), A*TANH(|T|/A): ",&
                !      TANH(t_norm/almo_scf_env%envelope_amplitude),&
                !      1.0_dp-(TANH(t_norm/almo_scf_env%envelope_amplitude))**2,&
                !      almo_scf_env%envelope_amplitude*TANH(t_norm/almo_scf_env%envelope_amplitude)
                !ENDIF
             ENDIF
             IF (optimize_theta) THEN
                CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
                CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
                        !func=dbcsr_func_sin,&
                        func=dbcsr_func_tanh,&
                        a0=0.0_dp,&
                        a1=1.0_dp/almo_scf_env%envelope_amplitude,&
                        error=error)
                CALL cp_dbcsr_scale(m_tmp_no_1,&
                        almo_scf_env%envelope_amplitude,error=error)
             ELSE
                CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
             ENDIF
             CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
                     quench_t(ispin),&
                     matrix_t_out(ispin),error=error)
             CALL cp_dbcsr_norm(matrix_t_out(ispin),&
                     dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)

             !IF (optimize_theta) THEN
             !   IF (unit_nr>0) THEN
             !      WRITE(unit_nr,*) "VAR1: |B*A*TANH(T/A)|: ", t_norm
             !   ENDIF
             !ELSE
             !   IF (unit_nr>0) THEN
             !      WRITE(unit_nr,*) "VAR1: |B*X|: ", t_norm
             !   ENDIF
             !ENDIF

             ! project out R_0
             IF (assume_t0_q0x.AND.(.NOT.just_started)) THEN
                IF (my_special_case.eq.1) THEN
                   CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                           STsiginv_0,&
                           matrix_t_out(ispin),&
                           0.0_dp,m_tmp_oo_1,&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                           matrix_t_0,&
                           m_tmp_oo_1,&
                           1.0_dp,matrix_t_out(ispin),&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
                ELSE IF (my_special_case.eq.0) THEN
                   ! cannot use projector with special case 0 - block-daigonal
                   CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                ELSE
                   ! no special case
                   IF (prec_type.eq.4) THEN
                      CALL apply_domain_operators(&
                              matrix_in=matrix_t_out(ispin),&
                              matrix_out=m_tmp_no_1,&
                              operator1=domain_r_down(:),&
                              operator2=almo_scf_env%domain_s_inv(:,ispin),&
                              dpattern=quench_t(ispin),&
                              map=almo_scf_env%domain_map(ispin),&
                              node_of_domain=almo_scf_env%cpu_of_domain,&
                              my_action=1,&
                              filter_eps=almo_scf_env%eps_filter,&
                              !matrix_trimmer=,&
                              use_trimmer=.FALSE.,&
                              error=error)
                   ELSE
                      CALL almo_scf_domain_operations(&
                              matrix_in=matrix_t_out(ispin),&
                              matrix_out=m_tmp_no_1,&
                              quench_t=quench_t(ispin),&
                              my_action=2,&
                              matrix_s=almo_scf_env%matrix_s(1),&
                              matrix_sigma_inv=matrix_sigma_inv_0,&
                              matrix_t=matrix_t_0,&
                              error=error)
                   ENDIF ! prec_type
                   CALL cp_dbcsr_copy(matrix_t_out(ispin),&
                           m_tmp_no_1,error=error)
                ENDIF ! special case
                CALL cp_dbcsr_norm(matrix_t_out(ispin),&
                        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                !IF (unit_nr>0) THEN
                !   WRITE(unit_nr,*) "VAR1: |Q.X|: ", t_norm
                !ENDIF
             ENDIF ! end assume_t0_q0x

             IF (assume_t0_q0x) THEN
                CALL cp_dbcsr_add(matrix_t_out(ispin),&
                        matrix_t_0,1.0_dp,1.0_dp,&
                        error=error)
                CALL cp_dbcsr_norm(matrix_t_out(ispin),&
                        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                !IF (unit_nr>0) THEN
                !   WRITE(unit_nr,*) "VAR1: |T0+Q.X|: ", t_norm
                !ENDIF
             ENDIF

             !!! RZK-warning an attempt (so far unsuccesful) to improve convergence
             !!! by keeping ALMOs as orthogonal as possible
             !!! before partial orthogonalization save the original state
             !!CALL cp_dbcsr_copy(m_tmp_no_1,matrix_t_out(ispin),error=error) 
             !!! make ALMOs as orthogonal as possible
             !!CALL most_orthogonal_almos(ket=matrix_t_out(ispin),&
             !!        overlap=almo_scf_env%matrix_sigma(ispin),&
             !!        overlap_blk=almo_scf_env%matrix_sigma_blk(ispin),&
             !!        metric=almo_scf_env%matrix_s(1),&
             !!        options=almo_scf_env%orthogonalization_opt,&
             !!        error=error)
             !!! find the change after orthogonalization
             !!CALL cp_dbcsr_add(m_tmp_no_1,matrix_t_out(ispin),-1.0_dp,1.0_dp,error=error) 
             !!CALL cp_dbcsr_scale(m_tmp_no_1,1.0_dp/step_size,error=error) 
             !!CALL cp_dbcsr_print(m_tmp_no_1,error=error) 
             !!! add the change to the step
             !!IF (iteration.gt.0) THEN
             !!   CALL cp_dbcsr_add(step,m_tmp_no_1,1.0_dp,1.0_dp,error=error) 
             !!ENDIF

             !! RZK-warning simple orthogonalization for block-diagonal ALMOs
             !CALL orthogonalize_mos(ket=matrix_t_out(ispin),&
             !        overlap=almo_scf_env%matrix_sigma_blk(ispin),&
             !        metric=almo_scf_env%matrix_s(1),&
             !        retain_locality=.TRUE.,&
             !        only_normalize=.TRUE.,&
             !        eps_filter=almo_scf_env%eps_filter,&
             !        order_lanczos=almo_scf_env%order_lanczos,&
             !        eps_lanczos=almo_scf_env%eps_lanczos,&
             !        max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
             !        error=error)

             CALL cp_dbcsr_filter(matrix_t_out(ispin),&
                     eps=almo_scf_env%eps_filter,&
                     error=error)

             ! compute the density matrix
             CALL almo_scf_t_to_p(&
                     t=matrix_t_out(ispin),&
                     p=almo_scf_env%matrix_p(ispin),&
                     eps_filter=almo_scf_env%eps_filter,&
                     orthog_orbs=.FALSE.,&
                     s=almo_scf_env%matrix_s(1),&
                     sigma=almo_scf_env%matrix_sigma(ispin),&
                     sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
                     use_guess=use_guess,&
                     error=error)
             CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin),&
                     spin_factor,error=error)

             !! RZK-warning to debug lets see eigenvalues of the ALMO overlap
             !CALL cp_dbcsr_copy(almo_scf_env%matrix_sigma_blk(ispin),&
             !               almo_scf_env%matrix_sigma(ispin),&
             !               keep_sparsity=.TRUE.,&
             !               error=error)
             !CALL cp_dbcsr_init(u_sigma,error=error)
             !CALL cp_dbcsr_create(u_sigma,template=almo_scf_env%matrix_sigma(ispin),&
             !        matrix_type=dbcsr_type_no_symmetry,error=error)
             !CALL cp_dbcsr_get_info(almo_scf_env%matrix_sigma_blk(ispin), nfullrows_total=occ1 )
             !ALLOCATE(evals(occ1))
             !CALL cp_dbcsr_syevd(almo_scf_env%matrix_sigma_blk(ispin),u_sigma,evals,&
             !        almo_scf_env%para_env,almo_scf_env%blacs_env,error=error)
             !WRITE(*,*) 'SIGMA_BLCK: ', evals
             !DEALLOCATE(evals)
             !CALL cp_dbcsr_create(u_sigma,template=almo_scf_env%matrix_sigma(ispin),&
             !        matrix_type=dbcsr_type_no_symmetry,error=error)
             !CALL cp_dbcsr_get_info(almo_scf_env%matrix_sigma(ispin), nfullrows_total=occ1 )
             !ALLOCATE(evals(occ1))
             !CALL cp_dbcsr_get_diag(almo_scf_env%matrix_sigma(ispin),evals,error=error)
             !WRITE(*,*) 'TRACE: ', SUM(evals)
             !CALL cp_dbcsr_syevd(almo_scf_env%matrix_sigma(ispin),u_sigma,evals,&
             !        almo_scf_env%para_env,almo_scf_env%blacs_env,error=error)
             !WRITE(*,*) 'SIGMA_FULL: ', evals
             !WRITE(*,*) 'SUMEI: ', SUM(evals)
             !DEALLOCATE(evals)
             !CALL cp_dbcsr_release(u_sigma,error=error)

             ! update the KS matrix and energy if necessary
             IF ( .NOT.(perturbation_only.AND.(.NOT.just_started)) ) THEN
                !!! RZK-warning the KS matrix must be updated outside the spin loop
                !!! Now the code works only for restricted orbitals
                !IF (unit_nr>0) THEN
                !   WRITE(unit_nr,*) "....updating KS matrix...."
                !ENDIF
                CALL almo_scf_dm_to_ks(qs_env,almo_scf_env,energy_new,error)
                CALL get_qs_env(qs_env, matrix_ks=matrix_ks, error=error)
                CALL matrix_qs_to_almo(matrix_ks(ispin)%matrix,&
                        almo_scf_env%matrix_ks(ispin),&
                        almo_scf_env,.FALSE.,error=error)
                CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin),&
                        almo_scf_env%eps_filter,error=error)
             ENDIF

             ! save the initial state
             IF (just_started) THEN
                CALL cp_dbcsr_copy(matrix_p_0,almo_scf_env%matrix_p(ispin),&
                        error=error)
                CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma(ispin),&
                        matrix_sigma_0,&
                        error=error)
                CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_inv(ispin),&
                        matrix_sigma_inv_0,&
                        error=error)
             ENDIF

             IF (my_special_case.eq.-1.AND.prec_type.eq.4) THEN
                ! construct domain-projector
                IF (assume_t0_q0x.AND.just_started) THEN
                   CALL construct_domain_r_down(&
                           matrix_t=matrix_t_0,&
                           matrix_sigma_inv=matrix_sigma_inv_0,&
                           matrix_s=almo_scf_env%matrix_s(1),&
                           subm_r_down=domain_r_down(:),&
                           dpattern=quench_t(ispin),&
                           map=almo_scf_env%domain_map(ispin),&
                           node_of_domain=almo_scf_env%cpu_of_domain,&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
                   !CALL construct_domain_r_down(&
                   !        matrix_t=matrix_t_out(ispin),&
                   !        matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
                   !        matrix_s=almo_scf_env%matrix_s(1),&
                   !        subm_r_down=domain_r_down(:),&
                   !        dpattern=quench_t(ispin),&
                   !        map=almo_scf_env%domain_map(ispin),&
                   !        node_of_domain=almo_scf_env%cpu_of_domain,&
                   !        filter_eps=almo_scf_env%eps_filter,&
                   !        error=error)
                ENDIF ! assume_t0_q0x
             ENDIF ! debug condition

             IF (perturbation_only) THEN
                ! calculate objective function Tr(F_0 R)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                        almo_scf_env%matrix_ks(ispin),&
                        matrix_t_out(ispin),&
                        0.0_dp,m_tmp_no_1,&
                        filter_eps=almo_scf_env%eps_filter,&
                        error=error)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                        m_tmp_no_1,&
                        almo_scf_env%matrix_sigma_inv(ispin),&
                        0.0_dp,FTsiginv,&
                        filter_eps=almo_scf_env%eps_filter,&
                        error=error)
                CALL cp_dbcsr_trace(matrix_t_out(ispin),&
                        FTsiginv,energy_new,"T","N",&
                        error=error)
                energy_new=energy_new*spin_factor
             ENDIF

             ! save the previous gradient to compute beta
             ! do it only if the previous grad was computed 
             ! for .NOT.line_search
             IF (line_search_iteration.eq.0.AND.iteration.ne.0) &
                CALL cp_dbcsr_copy(prev_grad,grad,error=error)
             
             ! compute the energy gradient if necessary
             skip_grad = ( iteration.gt.0 .AND. &
                           fixed_line_search_niter.ne.0 .AND. &
                           line_search_iteration.ne.fixed_line_search_niter )
                           
             IF (.NOT.skip_grad) THEN
                
                !IF (unit_nr>0) THEN
                !   WRITE(unit_nr,*) "....computing gradient...."
                !ENDIF

                ! do d_E/d_theta first
                IF (.NOT.perturbation_only) THEN
                   CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           almo_scf_env%matrix_ks(ispin),&
                           matrix_t_out(ispin),&
                           0.0_dp,m_tmp_no_1,&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_multiply("N","N",spin_factor,&
                           m_tmp_no_1,&
                           almo_scf_env%matrix_sigma_inv(ispin),&
                           0.0_dp,FTsiginv,&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
                   !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                   !        m_tmp_no_1,&
                   !        almo_scf_env%matrix_sigma_inv(ispin),&
                   !        0.0_dp,m_tmp_no_2,&
                   !        retain_sparsity=.TRUE.,&
                   !        error=error)
                ENDIF
                CALL cp_dbcsr_copy(m_tmp_no_2,quench_t(ispin),&
                        error=error)
                CALL cp_dbcsr_copy(m_tmp_no_2,&
                        FTsiginv,keep_sparsity=.TRUE.,error=error)
                CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                        matrix_t_out(ispin),&
                        FTsiginv,&
                        0.0_dp,m_tmp_oo_1,&
                        filter_eps=almo_scf_env%eps_filter,&
                        error=error)
                !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                !        m_tmp_oo_1,&
                !        almo_scf_env%matrix_sigma_inv(ispin),&
                !        0.0_dp,m_tmp_oo_2,&
                !        filter_eps=almo_scf_env%eps_filter,&
                !        error=error)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                        almo_scf_env%matrix_sigma_inv(ispin),&
                        m_tmp_oo_1,&
                        0.0_dp,siginvTFTsiginv,&
                        filter_eps=almo_scf_env%eps_filter,&
                        error=error)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                        almo_scf_env%matrix_s(1),&
                        matrix_t_out(ispin),&
                        0.0_dp,ST,&
                        filter_eps=almo_scf_env%eps_filter,&
                        error=error)
                ! save S*T_0*siginv_0
                IF (assume_t0_q0x.AND.just_started&
                   .AND.special_case.eq.1) THEN
                   CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           ST,&
                           almo_scf_env%matrix_sigma_inv(ispin),&
                           0.0_dp,STsiginv_0,&
                           filter_eps=almo_scf_env%eps_filter,&
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
                        eps=almo_scf_env%eps_filter,&
                        error=error)

                IF (perturbation_only.AND.just_started) THEN
                   ! save the first gradient
                   ! it is equal to F_vo and necessary to compute
                   ! the correction to the energy
                   CALL cp_dbcsr_copy(fvo_0,m_tmp_no_2,error=error)
                   CALL cp_dbcsr_scale(fvo_0,&
                           0.5_dp,error=error)
                ENDIF

                ! a short print-out
                CALL cp_dbcsr_norm(m_tmp_no_2,&
                        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                !IF (unit_nr>0) THEN
                !   WRITE(unit_nr,*) "Maximum norm of dE/dT: ", t_norm
                !ENDIF

                ! project out the occupied space from the gradient
                IF (assume_t0_q0x) THEN
                   IF (my_special_case.eq.1) THEN
                      CALL cp_dbcsr_copy(grad,m_tmp_no_2,error=error)
                      CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                              matrix_t_0,&
                              grad,&
                              0.0_dp,m_tmp_oo_1,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                      CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                              STsiginv_0,&
                              m_tmp_oo_1,&
                              1.0_dp,grad,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                   ELSE IF (my_special_case.eq.0) THEN
                      ! should not be here - cannot project the zero-order space from itself
                      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   ELSE
                      ! no special case
                      IF (prec_type.eq.4) THEN
                         CALL apply_domain_operators(&
                                 matrix_in=m_tmp_no_2,&
                                 matrix_out=grad,&
                                 operator2=domain_r_down(:),&
                                 operator1=almo_scf_env%domain_s_inv(:,ispin),&
                                 dpattern=quench_t(ispin),&
                                 map=almo_scf_env%domain_map(ispin),&
                                 node_of_domain=almo_scf_env%cpu_of_domain,&
                                 my_action=1,&
                                 filter_eps=almo_scf_env%eps_filter,&
                                 !matrix_trimmer=,&
                                 use_trimmer=.FALSE.,&
                                 error=error)
                      ELSE
                         CALL almo_scf_domain_operations(&
                                 matrix_in=m_tmp_no_2,&
                                 matrix_out=grad,&
                                 quench_t=quench_t(ispin),&
                                 my_action=1,&
                                 matrix_s=almo_scf_env%matrix_s(1),&
                                 matrix_sigma_inv=matrix_sigma_inv_0,&
                                 matrix_t=matrix_t_0,&
                                 error=error)
                      ENDIF ! prec_type
                   ENDIF ! my_special_case
                   ! If projecting the current occupied space is necessary use:
                   !CALL almo_scf_domain_operations(&
                   !        matrix_in=m_tmp_no_2,&
                   !        matrix_out=grad,&
                   !        quench_t=quench_t(ispin),&
                   !        my_action=1,&
                   !        matrix_s=almo_scf_env%matrix_s(1),&
                   !        matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
                   !        matrix_t=matrix_t_out(ispin),&
                   !        error=error)
                   CALL cp_dbcsr_copy(m_tmp_no_2,grad,error=error)
                ENDIF

                ! transform d_E/d_T to d_E/d_theta
                IF (optimize_theta) THEN
                   CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
                   CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
                           !func=dbcsr_func_cos,&
                           func=dbcsr_func_dtanh,&
                           a0=0.0_dp,&
                           a1=1.0_dp/almo_scf_env%envelope_amplitude,&
                           error=error)
                   CALL cp_dbcsr_scale(m_tmp_no_1,&
                           almo_scf_env%envelope_amplitude,&
                           error=error)
                   CALL cp_dbcsr_set(m_tmp_no_3,0.0_dp,error=error)
                   CALL cp_dbcsr_filter(m_tmp_no_3,&
                           eps=almo_scf_env%eps_filter,&
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
                ELSE ! simply copy
                   CALL cp_dbcsr_hadamard_product(m_tmp_no_2,&
                           quench_t(ispin),&
                           grad,&
                           error=error)
                ENDIF
                CALL cp_dbcsr_filter(grad,eps=almo_scf_env%eps_filter,&
                        error=error)

                !INFO CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                !INFO         grad,&
                !INFO         matrix_t_out(ispin),&
                !INFO         0.0_dp,m_tmp_oo_1,&
                !INFO         filter_eps=almo_scf_env%eps_filter,&
                !INFO         error=error)
                !INFO CALL cp_dbcsr_norm(m_tmp_oo_1,&
                !INFO         dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                !INFO IF (unit_nr>0) THEN
                !INFO    WRITE(unit_nr,*) "Grad orthogonality error: ", t_norm
                !INFO ENDIF

             ENDIF ! skip_grad
             
             ! check convergence and other exit criteria
             grad_norm_frob=cp_dbcsr_frobenius_norm(grad)
             CALL cp_dbcsr_norm(grad, dbcsr_norm_maxabsnorm,&
                     norm_scalar=grad_norm, error=error)
             converged=(grad_norm.lt.almo_scf_env%eps_scf_q)
             IF (converged.OR.(iteration.ge.max_iter)) THEN
                prepare_to_exit=.TRUE.
             ENDIF
             IF (grad_norm.lt.almo_scf_env%eps_prev_guess) THEN
                use_guess=.TRUE.
             ENDIF

             IF (md_in_theta_space) THEN

                IF (.NOT.first_md_iteration) THEN
                   CALL cp_dbcsr_copy(prev_step,step,error=error)
                ENDIF
                CALL cp_dbcsr_copy(step,grad,error=error)
                CALL cp_dbcsr_scale(step,-1.0_dp,error=error)
   
                ! update velocities v(i) = v(i-1) + 0.5*dT*(a(i-1) + a(i))
                IF (.NOT.first_md_iteration) THEN
                   CALL cp_dbcsr_add(velocity,&
                           step,1.0_dp,0.5_dp*time_step,error=error)
                   CALL cp_dbcsr_add(velocity,&
                           prev_step,1.0_dp,0.5_dp*time_step,error=error)
                ENDIF
                kin_energy=cp_dbcsr_frobenius_norm(velocity)
                kin_energy=0.5_dp*kin_energy*kin_energy
   
                ! update positions theta(i) = theta(i-1) + dT*v(i-1) + 0.5*dT*dT*a(i-1)
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

             ELSE ! optimizization (not MD)

                IF (.NOT.prepare_to_exit) THEN
      
                   ! check the gradient along the step direction
                   IF (iteration.ne.0) THEN

                      IF (fixed_line_search_niter.eq.0) THEN

                         CALL cp_dbcsr_trace(grad,step,line_search_error,&
                                 "T","N",error=error)
                         ! normalize the result
                         !IF (unit_nr>0) THEN
                         !   WRITE(unit_nr,*) "Angle between step/grad: ", line_search_error
                         !ENDIF
                         CALL cp_dbcsr_trace(grad,grad,denom,"T","N",&
                                 error=error)
                         !IF (unit_nr>0) THEN
                         !   WRITE(unit_nr,*) "Frobenius norm of grad:  ", SQRT(denom)
                         !ENDIF
                         line_search_error=line_search_error/SQRT(denom)
                         CALL cp_dbcsr_trace(step,step,denom,"T","N",&
                                 error=error)
                         !IF (unit_nr>0) THEN
                         !   WRITE(unit_nr,*) "Frobenius norm of step:  ", SQRT(denom)
                         !ENDIF
                         line_search_error=line_search_error/SQRT(denom)
                         IF (ABS(line_search_error).gt.almo_scf_env%eps_lin_search) THEN
                            line_search=.TRUE.
                            line_search_iteration=line_search_iteration+1
                         ELSE
                            line_search=.FALSE.
                            line_search_iteration=0
                            IF (grad_norm.lt.eps_skip_gradients) THEN
                               fixed_line_search_niter=ABS(almo_scf_env%integer04)
                            ENDIF
                         ENDIF
                      
                      ELSE ! decision for fixed_line_search_niter
                         
                         IF (.NOT.line_search) THEN
                            line_search=.TRUE.
                            line_search_iteration=line_search_iteration+1
                         ELSE
                            IF (line_search_iteration.eq.fixed_line_search_niter) THEN
                               line_search=.FALSE.
                               line_search_iteration=0
                               line_search_iteration=line_search_iteration+1
                            ENDIF
                         ENDIF
                      
                      ENDIF ! fixed_line_search_niter fork
                   ENDIF
   
                   IF (line_search) THEN
                         energy_diff=0.0_dp
                   ELSE
                         energy_diff=energy_new-energy_old
                         energy_old=energy_new
                   ENDIF
   
                   ! update the step direction
                   IF (.NOT.line_search) THEN
   
                      !IF (unit_nr>0) THEN
                      !   WRITE(unit_nr,*) "....updating step direction...."
                      !ENDIF
                      
                      cg_iteration=cg_iteration+1 
                      
                      IF ( (just_started .AND. perturbation_only) .OR. &
                           (iteration.eq.0 .AND. (.NOT.perturbation_only)) ) THEN
      
                         ! compute the preconditioner
                         !IF (unit_nr>0) THEN
                         !   WRITE(unit_nr,*) "....computing preconditioner...."
                         !ENDIF

                         ! calculate (1-R)F(1-R) and S-SRS
                         ! RZK-warning take advantage: some elements will be removed by the quencher
                         ! RZK-warning S operations can be performed outside the spin loop to save time
                         ! IT IS REQUIRED THAT PRECONDITIONER DOES NOT BREAK THE LOCALITY!!!!
                         ! RZK-warning: further optimization is ABSOLUTELY NECESSARY
                         !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                         !        almo_scf_env%matrix_s(1),&
                         !        matrix_t_out(ispin),&
                         !        0.0_dp,m_tmp_no_1,&
                         !        filter_eps=almo_scf_env%eps_filter,&
                         !        error=error)
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                                 ST,&
                                 almo_scf_env%matrix_sigma_inv(ispin),&
                                 0.0_dp,m_tmp_no_3,&
                                 filter_eps=almo_scf_env%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_create(m_tmp_nn_1,&
                                 template=almo_scf_env%matrix_s(1),&
                                 matrix_type=dbcsr_type_no_symmetry,error=error) 
                         CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_s(1),&
                                 m_tmp_nn_1,error=error)
                         IF (my_special_case.eq.1) THEN
                            ! use S instead of S-SRS
                         ELSE
                            CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
                                    ST,&
                                    m_tmp_no_3,&
                                    1.0_dp,m_tmp_nn_1,&
                                    filter_eps=almo_scf_env%eps_filter,&
                                    error=error)
                         ENDIF
                         !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                         !        almo_scf_env%matrix_ks(ispin),&
                         !        matrix_t_out(ispin),&
                         !        0.0_dp,m_tmp_no_1,&
                         !        filter_eps=almo_scf_env%eps_filter,&
                         !        error=error)
                         ! re-create nn matrix because desymmetrize is buggy -
                         ! it will create multiple copies of blocks
                         CALL cp_dbcsr_create(prec_vv,&
                                 template=almo_scf_env%matrix_ks(ispin),&
                                 matrix_type=dbcsr_type_no_symmetry,error=error) 
                         CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_ks(ispin),&
                                 prec_vv,error=error)
                         CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
                                 FTsiginv,&
                                 ST,&
                                 1.0_dp,prec_vv,&
                                 filter_eps=almo_scf_env%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
                                 ST,&
                                 FTsiginv,&
                                 1.0_dp,prec_vv,&
                                 filter_eps=almo_scf_env%eps_filter,&
                                 error=error)
                         !CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                         !        matrix_t_out(ispin),&
                         !        m_tmp_no_1,&
                         !        0.0_dp,m_tmp_oo_1,&
                         !        filter_eps=almo_scf_env%eps_filter,&
                         !        error=error)
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                                 ST,&
                                 siginvTFTsiginv,&
                                 0.0_dp,m_tmp_no_3,&
                                 filter_eps=almo_scf_env%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_multiply("N","T",1.0_dp,&
                                 m_tmp_no_3,&
                                 ST,&
                                 1.0_dp,prec_vv,&
                                 filter_eps=almo_scf_env%eps_filter,&
                                 error=error)
                         CALL cp_dbcsr_add(prec_vv,m_tmp_nn_1,&
                                 1.0_dp-prec_sf_mixing_s,&
                                 prec_sf_mixing_s,error=error)
                         CALL cp_dbcsr_scale(prec_vv,2.0_dp*spin_factor,error=error)
                         CALL cp_dbcsr_copy(m_tmp_nn_1,prec_vv,error=error)
                         
                         ! invert using various algorithms
                         IF (my_special_case.eq.0) THEN ! non-overlapping diagonal blocks
                            CALL invert_diagonal_blk(matrix_in=m_tmp_nn_1,&
                                    matrix_out=prec_vv,error=error)
                         ELSE IF (my_special_case.eq.1) THEN ! the entire system is a block
                            ! invert using cholesky (works with S matrix, will not work with S-SRS matrix)
                            CALL cp_dbcsr_cholesky_decompose(prec_vv,&
                                    para_env=almo_scf_env%para_env,&
                                    blacs_env=almo_scf_env%blacs_env,error=error)
                            CALL cp_dbcsr_cholesky_invert(prec_vv,&
                                    para_env=almo_scf_env%para_env,&
                                    blacs_env=almo_scf_env%blacs_env,&
                                    upper_to_full=.TRUE.,error=error)
                            CALL cp_dbcsr_filter(prec_vv,&
                                    eps=almo_scf_env%eps_filter,&
                                    error=error)
                         ELSE
                            !!! use a sophisticated domain preconditioner
                            IF (assume_t0_q0x) THEN
                               precond_domain_projector=-1
                            ELSE
                               precond_domain_projector=0
                            ENDIF
                            ! for other experimental preconditioner types the inversion is
                            ! done together with applying the preconditioner
                            IF (prec_type.eq.4) THEN
                               CALL construct_domain_preconditioner(&
                                  matrix_main=m_tmp_nn_1,&
                                  subm_s_inv=almo_scf_env%domain_s_inv(:,ispin),&
                                  subm_r_down=domain_r_down(:),&
                                  matrix_trimmer=quench_t(ispin),&
                                  dpattern=quench_t(ispin),&
                                  map=almo_scf_env%domain_map(ispin),&
                                  node_of_domain=almo_scf_env%cpu_of_domain,&
                                  preconditioner=almo_scf_env%domain_preconditioner(:,ispin),&
                                  use_trimmer=.FALSE.,&
                                  my_action=precond_domain_projector,&
                                  error=error)
                            ENDIF ! prec type
                         ENDIF

                         ! invert using cholesky (works with S matrix, will not work with S-SRS matrix)
                         !!!CALL cp_dbcsr_cholesky_decompose(prec_vv,&
                         !!!        para_env=almo_scf_env%para_env,&
                         !!!        blacs_env=almo_scf_env%blacs_env,error=error)
                         !!!CALL cp_dbcsr_cholesky_invert(prec_vv,&
                         !!!        para_env=almo_scf_env%para_env,&
                         !!!        blacs_env=almo_scf_env%blacs_env,&
                         !!!        upper_to_full=.TRUE.,error=error)
                         !!!CALL cp_dbcsr_filter(prec_vv,&
                         !!!        eps=almo_scf_env%eps_filter,&
                         !!!        error=error)
      
                         ! re-create the matrix because desymmetrize is buggy -
                         ! it will create multiple copies of blocks
                         !!!DESYM!CALL cp_dbcsr_create(prec_vv,&
                         !!!DESYM!        template=almo_scf_env%matrix_s(1),&
                         !!!DESYM!        matrix_type=dbcsr_type_no_symmetry,error=error) 
                         !!!DESYM!CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_s(1),&
                         !!!DESYM!        prec_vv,error=error)
                         !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                         !        almo_scf_env%matrix_s(1),&
                         !        matrix_t_out(ispin),&
                         !        0.0_dp,m_tmp_no_1,&
                         !        filter_eps=almo_scf_env%eps_filter,&
                         !        error=error)
                         !CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                         !        m_tmp_no_1,&
                         !        almo_scf_env%matrix_sigma_inv(ispin),&
                         !        0.0_dp,m_tmp_no_3,&
                         !        filter_eps=almo_scf_env%eps_filter,&
                         !        error=error)
                         !CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
                         !        m_tmp_no_3,&
                         !        m_tmp_no_1,&
                         !        1.0_dp,prec_vv,&
                         !        filter_eps=almo_scf_env%eps_filter,&
                         !        error=error)
                         !CALL cp_dbcsr_add_on_diag(prec_vv,&
                         !        prec_sf_mixing_s,error=error)
   
                         !CALL cp_dbcsr_create(prec_oo,&
                         !        template=almo_scf_env%matrix_sigma(ispin),&
                         !        matrix_type=dbcsr_type_no_symmetry,error=error) 
                         !CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma(ispin),&
                         !        matrix_type=dbcsr_type_no_symmetry,error=error) 
                         !CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma(ispin),&
                         !        prec_oo,error=error)
                         !CALL cp_dbcsr_filter(prec_oo,&
                         !        eps=almo_scf_env%eps_filter,&
                         !        error=error)
   
                         !! invert using cholesky
                         !CALL cp_dbcsr_create(prec_oo_inv,&
                         !        template=prec_oo,&
                         !        matrix_type=dbcsr_type_no_symmetry,error=error) 
                         !CALL cp_dbcsr_desymmetrize(prec_oo,&
                         !        prec_oo_inv,error=error)
                         !CALL cp_dbcsr_cholesky_decompose(prec_oo_inv,&
                         !        para_env=almo_scf_env%para_env,&
                         !        blacs_env=almo_scf_env%blacs_env,error=error)
                         !CALL cp_dbcsr_cholesky_invert(prec_oo_inv,&
                         !        para_env=almo_scf_env%para_env,&
                         !        blacs_env=almo_scf_env%blacs_env,&
                         !        upper_to_full=.TRUE.,error=error)

                      ENDIF
      
                      ! save the previous step
                      CALL cp_dbcsr_copy(prev_step,step,error=error)
      
                      ! compute the new step (apply preconditioner if available)
                      IF (use_preconditioner) THEN
                         
                         !IF (unit_nr>0) THEN
                         !   WRITE(unit_nr,*) "....applying preconditioner...."
                         !ENDIF

                         SELECT CASE (prec_type)
                         CASE (1)
                            ! expensive Newton-Raphson step (the Hessian is still approximate)
                            ! RZK-warning THIS PREC HAS NOT BEEN IMPLEMENTED FOR THETA
                            IF (ncores.gt.1) THEN                            
                               CPErrorMessage(cp_failure_level,routineP,"serial code only",error)
                               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                            ENDIF
                            CALL newton_grad_to_step(&
                                    matrix_grad=m_tmp_no_2,&
                                    matrix_step=m_tmp_no_1,&
                                    matrix_s=almo_scf_env%matrix_s(1),&
                                    matrix_ks=almo_scf_env%matrix_ks(ispin),&
                                    matrix_t=matrix_t_out(ispin),&
                                    matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
                                    !matrix_ks=matrix_ks_0,&
                                    !matrix_t=matrix_t_0,&
                                    !matrix_sigma_inv=matrix_sigma_inv_0,&
                                    quench_t=quench_t(ispin),&
                                    spin_factor=spin_factor,&
                                    eps_filter=almo_scf_env%eps_filter,&
                                    error=error)
                         CASE (2)
                            IF (my_special_case.eq.0.OR.my_special_case.eq.1) THEN
                               !IF (unit_nr>0) THEN
                               !   WRITE(unit_nr,*) "....special case ",my_special_case,"...."
                               !ENDIF
                               CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                                       prec_vv,&
                                       grad,&
                                       0.0_dp,step,&
                                       filter_eps=almo_scf_env%eps_filter,&
                                       error=error)
                            ELSE
                               ! block-by-block inversion
                               ! create trimmer
                               IF (optimize_theta) THEN
                                  CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
                                  CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
                                          !func=dbcsr_func_cos,&
                                          func=dbcsr_func_dtanh,&
                                          a0=0.0_dp,&
                                          a1=1.0_dp/almo_scf_env%envelope_amplitude,&
                                          error=error)
                                  CALL cp_dbcsr_scale(m_tmp_no_1,&
                                          almo_scf_env%envelope_amplitude,&
                                          error=error)
                                  CALL cp_dbcsr_set(m_tmp_no_3,0.0_dp,error=error)
                                  CALL cp_dbcsr_filter(m_tmp_no_3,&
                                          eps=almo_scf_env%eps_filter,&
                                          error=error)
                                  CALL cp_dbcsr_hadamard_product(&
                                          quench_t(ispin),&
                                          m_tmp_no_1,&
                                          m_tmp_no_3,&
                                          b_assume_value=1.0_dp,&
                                          error=error)
                                  !!!! EXPERIMENT
                                  !CALL cp_dbcsr_copy(m_tmp_no_1,&
                                  !        m_tmp_no_3,&
                                  !        error=error)
                                  !CALL cp_dbcsr_copy(m_tmp_no_3,&
                                  !        quench_t_saved,&
                                  !        error=error)
                                  !!!! END EXPERIMENT
                               ELSE ! optimize MOs
                                  CALL cp_dbcsr_copy(m_tmp_no_3,&
                                           quench_t(ispin),&
                                           error=error)
                               ENDIF
   
                               !!!!WRITE(*,*) "TRIMMER: "
                               !!!!CALL cp_dbcsr_print(m_tmp_no_3,error=error)
                               !!!!WRITE(*,*) "GRAD: "
                               !!!!CALL cp_dbcsr_print(grad,error=error)
                               !!!!!!!!!! ADD_ON REMOVE LATER !!!!
                               IF (optimize_theta) THEN
                                  CALL almo_scf_domain_operations(&
                                          !!! RZK-warning: check what is given as matrix_in
                                          matrix_in=m_tmp_no_2,&
                                          matrix_out=step,&
                                          matrix_main=m_tmp_nn_1,&
                                          quench_t=quench_t(ispin),&
                                          my_action=-1,&
                                          matrix_s=almo_scf_env%matrix_s(1),&
                                          matrix_sigma=almo_scf_env%matrix_sigma(ispin),&
                                          matrix_sigma_inv=matrix_sigma_inv_0,&
                                          matrix_t=matrix_t_0,&
                                          !!! RZK-warning: check use_trimmer (TRUE?)
                                          use_trimmer=.FALSE.,&
                                          matrix_trimmer=m_tmp_no_3,&
                                          error=error)
   
                                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                  !!WRITE(*,*) "STEP: "
                                  !!CALL cp_dbcsr_print(step,error=error)
                                  !CALL almo_scf_domain_operations(&
                                  !        matrix_in=step,&
                                  !        matrix_out=m_tmp_no_3,&
                                  !        !matrix_main=m_tmp_nn_1,&
                                  !        quench_t=quench_t(ispin),&
                                  !        my_action=2,&
                                  !        matrix_s=almo_scf_env%matrix_s(1),&
                                  !        !matrix_sigma=almo_scf_env%matrix_sigma(ispin),&
                                  !        matrix_sigma_inv=matrix_sigma_inv_0,&
                                  !        matrix_t=matrix_t_0,&
                                  !        use_trimmer=.FALSE.,&
                                  !        error=error)
                                  !CALL cp_dbcsr_copy(step,m_tmp_no_3,error=error)
                                  !!!CALL cp_dbcsr_print(step,error=error)
   
                                  !! RZK-warning trimmer is ready - m_tmp_no_1, invert it now
                                  !! it is a dangerous operation - there are zeros there
                                  shift_value=1.0E-7_dp
                                  !CALL cp_dbcsr_print(m_tmp_no_1,error=error)
                                  CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
                                          func=dbcsr_func_inverse_special,&
                                          a0=shift_value,&
                                          error=error)
                                  !CALL cp_dbcsr_print(m_tmp_no_1,error=error)
                                  CALL cp_dbcsr_copy(m_tmp_no_3,&
                                          step,&
                                          error=error)
                                  CALL cp_dbcsr_set(step,0.0_dp,error=error)
                                  CALL cp_dbcsr_filter(step,&
                                          eps=almo_scf_env%eps_filter,&
                                          error=error)
                                  CALL cp_dbcsr_hadamard_product(&
                                          m_tmp_no_3,&
                                          m_tmp_no_1,&
                                          step,&
                                          b_assume_value=1.0_dp/shift_value,&
                                          error=error)
                                  
                                  ! evaluate the trimmer again!!!
                                  CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
                                  CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
                                          !func=dbcsr_func_cos,&
                                          func=dbcsr_func_dtanh,&
                                          a0=0.0_dp,&
                                          a1=1.0_dp/almo_scf_env%envelope_amplitude,&
                                          error=error)
                                  CALL cp_dbcsr_scale(m_tmp_no_1,&
                                          almo_scf_env%envelope_amplitude,&
                                          error=error)
                                  CALL cp_dbcsr_set(m_tmp_no_3,0.0_dp,error=error)
                                  CALL cp_dbcsr_filter(m_tmp_no_3,&
                                          eps=almo_scf_env%eps_filter,&
                                          error=error)
                                  CALL cp_dbcsr_hadamard_product(&
                                          quench_t(ispin),&
                                          m_tmp_no_1,&
                                          m_tmp_no_3,&
                                          b_assume_value=1.0_dp,&
                                          error=error)
                                  CALL almo_scf_domain_operations(&
                                          matrix_in=grad,&
                                          matrix_out=step,&
                                          matrix_main=m_tmp_nn_1,&
                                          quench_t=quench_t(ispin),&
                                          my_action=-2,&
                                          matrix_s=almo_scf_env%matrix_s(1),&
                                          matrix_sigma=almo_scf_env%matrix_sigma(ispin),&
                                          matrix_sigma_inv=matrix_sigma_inv_0,&
                                          matrix_t=matrix_t_0,&
                                          use_trimmer=.TRUE.,&
                                          matrix_trimmer=m_tmp_no_3,&
                                          error=error)
                                  
                                  CALL cp_dbcsr_print(step,error=error)             
   
                               ELSE ! non-theta optimization
                                  CALL almo_scf_domain_operations(&
                                          !!! RZK-warning: check what is given as matrix_in
                                          matrix_in=m_tmp_no_2,&
                                          matrix_out=step,&
                                          matrix_main=m_tmp_nn_1,&
                                          quench_t=quench_t(ispin),&
                                          my_action=-1,&
                                          matrix_s=almo_scf_env%matrix_s(1),&
                                          matrix_sigma=almo_scf_env%matrix_sigma(ispin),&
                                          matrix_sigma_inv=matrix_sigma_inv_0,&
                                          matrix_t=matrix_t_0,&
                                          !!! RZK-warning: check use_trimmer
                                          use_trimmer=.FALSE.,&
                                          matrix_trimmer=m_tmp_no_3,&
                                          error=error)
                                  !CALL almo_scf_domain_operations(&
                                  !        matrix_in=m_tmp_no_2,&
                                  !        matrix_out=step,&
                                  !        matrix_main=m_tmp_nn_1,&
                                  !        quench_t=quench_t(ispin),&
                                  !        my_action=-2,&
                                  !        matrix_s=almo_scf_env%matrix_s(1),&
                                  !        matrix_sigma=almo_scf_env%matrix_sigma(ispin),&
                                  !        matrix_sigma_inv=matrix_sigma_inv_0,&
                                  !        matrix_t=matrix_t_0,&
                                  !        use_trimmer=.FALSE.,&
                                  !        matrix_trimmer=m_tmp_no_3,&
                                  !        error=error)
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
                                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                  !!WRITE(*,*) "STEP: "
                                  !!CALL cp_dbcsr_print(step,error=error)
                                  !!CALL almo_scf_domain_operations(&
                                  !!        matrix_in=step,&
                                  !!        matrix_out=m_tmp_no_1,&
                                  !!        !matrix_main=m_tmp_nn_1,&
                                  !!        quench_t=quench_t(ispin),&
                                  !!        my_action=2,&
                                  !!        matrix_s=almo_scf_env%matrix_s(1),&
                                  !!        !matrix_sigma=almo_scf_env%matrix_sigma(ispin),&
                                  !!        matrix_sigma_inv=matrix_sigma_inv_0,&
                                  !!        matrix_t=matrix_t_0,&
                                  !!        use_trimmer=.FALSE.,&
                                  !!        !matrix_trimmer=m_tmp_no_3,&
                                  !!        error=error)
                                  !!CALL cp_dbcsr_copy(step,m_tmp_no_1,error=error)
                                  !!!CALL cp_dbcsr_print(step,error=error)
                                  
                                  CALL cp_dbcsr_copy(m_tmp_no_3,quench_t(ispin),error=error)
                                  CALL almo_scf_domain_operations(&
                                          matrix_in=grad,&
                                          matrix_out=step,&
                                          matrix_main=m_tmp_nn_1,&
                                          quench_t=quench_t(ispin),&
                                          my_action=-2,&
                                          matrix_s=almo_scf_env%matrix_s(1),&
                                          matrix_sigma=almo_scf_env%matrix_sigma(ispin),&
                                          matrix_sigma_inv=matrix_sigma_inv_0,&
                                          matrix_t=matrix_t_0,&
                                          use_trimmer=.TRUE.,&
                                          matrix_trimmer=m_tmp_no_3,&
                                          error=error)
                               ENDIF
                               
                               CALL cp_dbcsr_filter(step,almo_scf_env%eps_filter,error=error)
   
                               CALL cp_dbcsr_create(m_tmp_oo_1,&
                                       template=almo_scf_env%matrix_sigma_blk(ispin),&
                                       matrix_type=dbcsr_type_no_symmetry,error=error) 
                               CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_blk(ispin),m_tmp_oo_1,error=error)
                               CALL get_overlap(bra=matrix_t_out(ispin),&
                                       ket=step,&
                                       overlap=m_tmp_oo_1,&
                                       metric=almo_scf_env%matrix_s(1),&
                                       retain_overlap_sparsity=.TRUE.,&
                                       eps_filter=almo_scf_env%eps_filter,&
                                       error=error)
                               CALL cp_dbcsr_norm(m_tmp_oo_1,&
                                       dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                               IF (unit_nr>0) THEN
                                  WRITE(unit_nr,*) "Step block-orthogonality error: ", t_norm 
                               ENDIF

                            ENDIF ! no special cases

                         CASE (3)

                            ! RZK-warning THIS PREC HAS NOT BEEN IMPLEMENTED FOR THETA

                            ! inversion
                            CALL cp_dbcsr_get_info(m_tmp_nn_1, nfullrows_total=dim0 )
                            ALLOCATE(evals(dim0))
                            CALL cp_dbcsr_syevd(m_tmp_nn_1,prec_vv,evals,&
                                    almo_scf_env%para_env,almo_scf_env%blacs_env,error=error)
                            ! invert eigenvalues and use eigenvectors to compute the Hessian inverse
                            ! take special care of zero eigenvalues
                            zero_neg_eiv=0
                            CALL cp_dbcsr_get_info(almo_scf_env%matrix_sigma(ispin), nfullrows_total=occ1 )
                            DO jj=1, dim0
                               IF (jj.le.occ1) THEN
                                  evals(jj)=evals(jj)*0.0_dp
                                  zero_neg_eiv=zero_neg_eiv+1
                               ELSE
                                  evals(jj)=1.0_dp/evals(jj)
                               ENDIF
                            ENDDO
                            IF (unit_nr>0) THEN
                               WRITE(*,*) 'ZERO OR NEGATIVE EIGENVALUES: ', zero_neg_eiv
                            ENDIF
                            CALL cp_dbcsr_init(inv_eiv,error=error)
                            CALL cp_dbcsr_create(inv_eiv,&
                                 template=m_tmp_nn_1,&
                                 matrix_type=dbcsr_type_no_symmetry,error=error) 
                            CALL cp_dbcsr_add_on_diag(inv_eiv,1.0_dp,error=error)
                            CALL cp_dbcsr_set_diag(inv_eiv,evals,error=error)
                            CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                                    prec_vv,&
                                    inv_eiv,&
                                    0.0_dp,m_tmp_nn_1,&
                                    filter_eps=almo_scf_env%eps_filter,&
                                    error=error)
                            CALL cp_dbcsr_multiply("N","T",1.0_dp,&
                                    m_tmp_nn_1,&
                                    prec_vv,&
                                    0.0_dp,inv_eiv,&
                                    filter_eps=almo_scf_env%eps_filter,&
                                    error=error)
                            CALL cp_dbcsr_copy(prec_vv,inv_eiv,error=error)
                            CALL cp_dbcsr_release(inv_eiv,error=error)
                            DEALLOCATE(evals)
   
                            !!CALL cp_dbcsr_copy(step,&
                            !!        quench_t(ispin),&
                            !!        error=error)
                            !!CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                            !!        m_tmp_no_2,&
                            !!        !grad,& - this choice is worse
                            !!        prec_oo,&
                            !!        0.0_dp,step,&
                            !!        !retain_sparsity=.TRUE.,&
                            !!        filter_eps=almo_scf_env%eps_filter,&
                            !!        error=error)
                            CALL cp_dbcsr_copy(m_tmp_no_1,&
                                    quench_t(ispin),&
                                    error=error)
                            !!CALL cp_dbcsr_hadamard_product(&
                            !!        quench_t(ispin),&
                            !!        step,&
                            !!        m_tmp_no_1,&
                            !!        error=error)
                            CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                                    prec_vv,&
                                    m_tmp_no_2,&
                                    0.0_dp,m_tmp_no_1,&
                                    retain_sparsity=.TRUE.,&
                                    error=error)
   
                         CASE (4)
              
                            IF (my_special_case.eq.0.OR.my_special_case.eq.1) THEN
                               
                               !IF (unit_nr>0) THEN
                               !   WRITE(unit_nr,*) "....special case ",my_special_case,"...."
                               !ENDIF
                               CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
                                       prec_vv,&
                                       grad,&
                                       0.0_dp,step,&
                                       filter_eps=almo_scf_env%eps_filter,&
                                       error=error)

                            ELSE
                               
                               !!! RZK-warning Currently for non-theta only
                               IF (optimize_theta) THEN
                                  CPErrorMessage(cp_failure_level,routineP,"theta is NYI",error)
                                  CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                               ENDIF
   
                               CALL apply_domain_operators(&
                                       matrix_in=grad,&
                                       matrix_out=step,&
                                       operator1=almo_scf_env%domain_preconditioner(:,ispin),&
                                       !operator2=,&
                                       dpattern=quench_t(ispin),&
                                       map=almo_scf_env%domain_map(ispin),&
                                       node_of_domain=almo_scf_env%cpu_of_domain,&
                                       my_action=0,&
                                       filter_eps=almo_scf_env%eps_filter,&
                                       !matrix_trimmer=,&
                                       !use_trimmer=.FALSE.,&
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
                                                           
                               !CALL cp_dbcsr_create(m_tmp_oo_1,&
                               !        template=almo_scf_env%matrix_sigma_blk(ispin),&
                               !        matrix_type=dbcsr_type_no_symmetry,error=error) 
                               !CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_blk(ispin),m_tmp_oo_1,error=error)
                               !CALL get_overlap(bra=matrix_t_out(ispin),&
                               !        ket=step,&
                               !        overlap=m_tmp_oo_1,&
                               !        metric=almo_scf_env%matrix_s(1),&
                               !        retain_overlap_sparsity=.TRUE.,&
                               !        eps_filter=almo_scf_env%eps_filter,&
                               !        error=error)
                               !CALL cp_dbcsr_norm(m_tmp_oo_1,&
                               !        dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                               !IF (unit_nr>0) THEN
                               !   WRITE(unit_nr,*) "Step block-orthogonality error: ", t_norm 
                               !ENDIF

                            ENDIF ! special case
               
                         END SELECT ! preconditioner type fork

                      ELSE

                         !!!! NO PRECONDITIONER
                         CALL cp_dbcsr_copy(step,grad,error=error)
                         CALL cp_dbcsr_scale(step,-1.0_dp,error=error)

                      ENDIF
   
                      !INFO CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                      !INFO         almo_scf_env%matrix_s(1),&
                      !INFO         matrix_t_out(ispin),&
                      !INFO         0.0_dp,m_tmp_no_1,&
                      !INFO         filter_eps=almo_scf_env%eps_filter,&
                      !INFO         error=error)
                      !INFO CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                      !INFO         step,&
                      !INFO         m_tmp_no_1,&
                      !INFO         0.0_dp,m_tmp_oo_1,&
                      !INFO         filter_eps=almo_scf_env%eps_filter,&
                      !INFO         error=error)
                      !INFO CALL cp_dbcsr_norm(m_tmp_oo_1,&
                      !INFO         dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                      !INFO IF (unit_nr>0) THEN
                      !INFO    WRITE(unit_nr,*) "Dir orthogonality error: ", t_norm
                      !INFO ENDIF
   
                      ! check whether we need to reset conjugate directions
                      IF (iteration.eq.0) THEN
                         reset_conjugator=.TRUE.
                      ENDIF
   
                      ! compute the conjugation coefficient - beta
                      IF (.NOT.reset_conjugator) THEN
     
                         SELECT CASE (almo_scf_env%scf_conjugator)
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
                   
                      ! save the preconditioned gradient (useful for beta)
                      CALL cp_dbcsr_copy(prev_minus_prec_grad,step,error=error)
      
                      !IF (unit_nr>0) THEN
                      !   WRITE(unit_nr,*) "....final beta....", beta
                      !ENDIF
                      
                      ! conjugate the step direction
                      CALL cp_dbcsr_add(step,prev_step,1.0_dp,beta,error=error)
                      !CALL almo_scf_domain_operations(&
                      !        matrix_in=step,&
                      !        matrix_out=m_tmp_no_3,&
                      !        matrix_s=almo_scf_env%matrix_s(1),&
                      !        matrix_sigma=almo_scf_env%matrix_sigma(ispin),&
                      !        matrix_sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
                      !        matrix_t=matrix_t_out(ispin),&
                      !        quench_t=quench_t(ispin),&
                      !        my_action=2,&
                      !        error=error)
                      !CALL cp_dbcsr_copy(step,m_tmp_no_3,error=error)
                      !CALL cp_dbcsr_print(step,error=error)
   
                   ENDIF ! update the step direction 
      
                   ! estimate the step size 
                   IF (.NOT.line_search) THEN
                      e0=energy_new
                      CALL cp_dbcsr_trace(grad,step,g0,"T","N",error=error)
                      ! we just changed the direction and
                      ! we have only E and grad from the current step
                      ! it is not enouhg to compute step_size - just guess it
                      IF (iteration.eq.0) THEN
                         step_size=almo_scf_env%real05 ! wild guess
                      ELSE
                         IF (next_step_size_guess.le.0.0_dp) THEN
                            step_size=almo_scf_env%real05 ! wild guess
                         ELSE
                            ! take the last value 
                            step_size=next_step_size_guess*1.05_dp
                         ENDIF
                      ENDIF
                      next_step_size_guess=step_size
                   ELSE
                      IF (fixed_line_search_niter.eq.0) THEN
                         e1=energy_new
                         CALL cp_dbcsr_trace(grad,step,g1,"T","N",error=error)
                         ! we have accumulated some points along this direction
                         ! use only the most recent g0 (quadratic approximation)
                         appr_sec_der=(g1-g0)/step_size
                         !IF (unit_nr>0) THEN
                         !   WRITE(unit_nr,'(A2,7F12.5)') &
                         !           "EG",e0,e1,g0,g1,appr_sec_der,step_size,-g1/appr_sec_der
                         !ENDIF
                         step_size=-g1/appr_sec_der
                         e0=e1
                         g0=g1
                      ELSE
                         ! use e0, g0 and e1 to compute g1 and make a step
                         ! if the next iteration is also line_search
                         ! use e1 and the calculated g1 as e0 and g0
                         e1=energy_new
                         appr_sec_der=2.0*( (e1-e0)/step_size - g0 )/step_size
                         g1=appr_sec_der*step_size + g0
                         IF (unit_nr>0) THEN
                            WRITE(unit_nr,'(A2,7F12.5)') &
                                    "EG",e0,e1,g0,g1,appr_sec_der,step_size,-g1/appr_sec_der
                         ENDIF
                         !appr_sec_der=(g1-g0)/step_size
                         step_size=-g1/appr_sec_der
                         e0=e1
                         g0=g1
                      ENDIF
                      next_step_size_guess=next_step_size_guess+step_size
                   ENDIF
   
                   !INFO CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                   !INFO         almo_scf_env%matrix_s(1),&
                   !INFO         matrix_t_out(ispin),&
                   !INFO         0.0_dp,m_tmp_no_1,&
                   !INFO         filter_eps=almo_scf_env%eps_filter,&
                   !INFO         error=error)
                   !INFO CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                   !INFO         step,&
                   !INFO         m_tmp_no_1,&
                   !INFO         0.0_dp,m_tmp_oo_1,&
                   !INFO         filter_eps=almo_scf_env%eps_filter,&
                   !INFO         error=error)
                   !INFO CALL cp_dbcsr_norm(m_tmp_oo_1,&
                   !INFO         dbcsr_norm_maxabsnorm, norm_scalar=t_norm, error=error)
                   !INFO IF (unit_nr>0) THEN
                   !INFO    WRITE(unit_nr,*) "Step orthogonality error: ", t_norm
                   !INFO ENDIF
   
                   ! update theta
                   CALL cp_dbcsr_add(m_theta,step,1.0_dp,step_size,error=error)
   
                ENDIF ! not.prepare_to_exit
              
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
                   !WRITE(unit_nr,'(T2,A11,I6,F20.12,E12.3,E12.3,E12.3,F12.5,F10.3)') &
                   !        "ALMO SCF ",iter_type,iteration,&
                   !        energy_new,energy_diff,grad_norm,line_search_error,&
                   !        step_size,t2-t1
                ENDIF

                t1 = m_walltime()
   
   
             ENDIF ! MD in theta space


             iteration=iteration+1
             IF (prepare_to_exit) EXIT

          ENDDO ! inner SCF loop

          IF (converged.OR.(outer_iteration.ge.outer_max_iter)) THEN
             outer_prepare_to_exit=.TRUE.
          ENDIF
          
          outer_iteration=outer_iteration+1
          IF (outer_prepare_to_exit) EXIT

       ENDDO ! outer SCF loop
       
       ! post SCF-loop calculations
       IF (converged)  THEN
          
          ! obtain MO coefficients from final theta
          ! DANGER: if decide to uncomment make sure that Theta->T
          ! procedure is consistent with the procedure in the loop
          !IF (optimize_theta) THEN
          !   CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
          !   CALL cp_dbcsr_function_of_elements(m_tmp_no_1,&
          !           func=dbcsr_func_tanh,&
          !           a0=0.0_dp,&
          !           a1=1.0_dp/almo_scf_env%envelope_amplitude,&
          !           error=error)
          !   CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
          !           quench_t(ispin),&
          !           matrix_t_out(ispin),error=error)
          !   CALL cp_dbcsr_scale(matrix_t_out(ispin),&
          !           almo_scf_env%envelope_amplitude,error=error)
          !ELSE
          !   CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
          !   CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
          !           quench_t(ispin),&
          !           matrix_t_out(ispin),error=error)
          !ENDIF
          !IF (perturbation_only) THEN
          !   CALL cp_dbcsr_add(matrix_t_out(ispin),&
          !           matrix_t_0,1.0_dp,1.0_dp,&
          !           error=error)
          !ENDIF
          !CALL cp_dbcsr_filter(matrix_t_out(ispin),&
          !        eps=almo_scf_env%eps_filter,&
          !        error=error)
          !CALL cp_dbcsr_norm(matrix_t_out(ispin),&
          !        dbcsr_norm_maxabsnorm, norm_scalar=grad_norm, error=error)
          !IF (unit_nr>0) THEN
          !   WRITE(unit_nr,*) "Maximum norm of the ALMOs: ", grad_norm
          !ENDIF
             
          !!!!! experiment: bump final amplitudes and get T
          !!!!IF (.NOT.optimize_theta.AND.perturbation_only) THEN
          !!!!   CALL cp_dbcsr_copy(m_tmp_no_1,m_theta,error=error)
          !!!!   CALL cp_dbcsr_hadamard_product(m_tmp_no_1,&
          !!!!           quench_t_saved,&
          !!!!           matrix_t_out(ispin),error=error)
          !!!!   IF (use_projector) THEN
          !!!!      CALL almo_scf_domain_operations(&
          !!!!              matrix_in=matrix_t_out(ispin),&
          !!!!              matrix_out=m_tmp_no_1,&
          !!!!              quench_t=quench_t(ispin),&
          !!!!              my_action=2,&
          !!!!              matrix_s=almo_scf_env%matrix_s(1),&
          !!!!              matrix_sigma_inv=matrix_sigma_inv_0,&
          !!!!              matrix_t=matrix_t_0,&
          !!!!              error=error)
          !!!!      CALL cp_dbcsr_copy(matrix_t_out(ispin),&
          !!!!              m_tmp_no_1,error=error)
          !!!!   ENDIF
          !!!!   IF (perturbation_only) THEN
          !!!!      CALL cp_dbcsr_add(matrix_t_out(ispin),&
          !!!!              matrix_t_0,1.0_dp,1.0_dp,&
          !!!!              error=error)
          !!!!   ENDIF
          !!!!ENDIF
          !!!!
          !!!!IF (.NOT.optimize_theta) THEN
          !!!!   CALL cp_dbcsr_copy(quench_t(ispin),quench_t_saved,error=error)
          !!!!ENDIF

          IF (perturbation_only) THEN
            
             CALL cp_dbcsr_add(matrix_t_0,matrix_t_out(ispin),&
                     -1.0_dp,1.0_dp,&
                     error=error)

             CALL cp_dbcsr_trace(matrix_t_0,&
                     fvo_0,energy_new,"T","N",&
                     error=error)
             ! print out the energy lowering
             IF (unit_nr>0) THEN
                WRITE(unit_nr,*)
                WRITE(unit_nr,'(T2,A35,F25.10)') "ENERGY OF BLOCK-DIAGONAL ALMOs:",&
                   almo_scf_env%almo_scf_energy
                WRITE(unit_nr,'(T2,A35,F25.10)') "ENERGY LOWERING:",&
                   energy_new
                WRITE(unit_nr,'(T2,A35,F25.10)') "CORRECTED ENERGY:",&
                   almo_scf_env%almo_scf_energy+energy_new
                WRITE(unit_nr,*)
             ENDIF
             !IF (unit_nr>0) THEN
             !   WRITE(unit_nr,*) "_ENERGY-0: ", almo_scf_env%almo_scf_energy
             !   WRITE(unit_nr,*) "_ENERGY-D: ", energy_new
             !   WRITE(unit_nr,*) "_ENERGY-F: ", almo_scf_env%almo_scf_energy+energy_new
             !ENDIF
             CALL almo_scf_update_ks_energy(qs_env,&
                     almo_scf_env%almo_scf_energy+energy_new,&
                     error)

             ! similar method to evaluate the energy correction
             !CALL cp_dbcsr_add(matrix_p_0,&
             !                 almo_scf_env%matrix_p(ispin),-1.0_dp,1.0_dp,&
             !                 error=error)
             !CALL cp_dbcsr_trace(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
             !                    matrix_p_0,&
             !                    energy_new,error=error)
             !IF (unit_nr>0) THEN
             !   WRITE(unit_nr,*) "alt-ENERGY-D: ", energy_new
             !ENDIF
             !CALL almo_scf_update_ks_energy(qs_env,&
             !        almo_scf_env%almo_scf_energy+energy_new,&
             !        error)

             IF (almo_scf_env%almo_eda.gt.0) THEN

                ! print out the results of decomposition analysis
                CALL cp_dbcsr_hadamard_product(matrix_t_0,&
                        fvo_0,m_tmp_no_1,error=error)
                IF (unit_nr>0) THEN
                   WRITE(unit_nr,*)
                   WRITE(unit_nr,'(T2,A)') "DECOMPOSITION OF THE DELOCALIZATION ENERGY"
                ENDIF
                CALL cp_dbcsr_filter(m_tmp_no_1,almo_scf_env%eps_filter,error=error)

                mynode=dbcsr_mp_mynode(dbcsr_distribution_mp(&
                   cp_dbcsr_distribution(m_tmp_no_1)))
                WRITE(mynodestr,'(I6.6)') mynode
                mylogfile='EDA.'//TRIM(ADJUSTL(mynodestr))
                OPEN (iunit,file=mylogfile,status='REPLACE')
                
                CALL cp_dbcsr_print_block_sum(m_tmp_no_1,iunit,error=error)
                CLOSE(iunit)
                
                !IF (unit_nr>0) THEN
                !   WRITE(unit_nr,*)
                !   WRITE(unit_nr,'(T2,A)') "CHARGE DECOMPOSITION"
                !ENDIF
                !CALL cp_dbcsr_print_block_sum(cta_matrix(ispin),error=error)

             ENDIF ! do ALMO EDA

          ENDIF ! if perturbation only

       ENDIF ! if converged

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

    ENDDO ! ispin
    
    CALL timestop(handle)

  END SUBROUTINE almo_scf_xalmo_pcg

! *****************************************************************************
!> \brief selects various post scf routines
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_delocalization(qs_env,almo_scf_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_delocalization', &
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

    SELECT CASE (almo_scf_env%deloc_method)
    CASE (almo_deloc_x,almo_deloc_full_scf,almo_deloc_x_then_full_scf)

       ! RZK-warning hack into the quenched routine:
       ! create a quench matrix with all-all-all blocks 1.0
       ! it is a waste of memory but since matrices are distributed
       ! we can tolerate it for now
       ALLOCATE(no_quench(almo_scf_env%nspins))
       CALL cp_dbcsr_init(no_quench(1),error=error)
       CALL cp_dbcsr_create(no_quench(1),&
               template=almo_scf_env%matrix_t(1),&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(&
          cp_dbcsr_distribution(no_quench(1))))
       CALL cp_dbcsr_work_create(no_quench(1),&
               work_mutable=.TRUE., error=error)
       nblkrows_tot = cp_dbcsr_nblkrows_total(no_quench(1))
       nblkcols_tot = cp_dbcsr_nblkcols_total(no_quench(1))
       ! RZK-warning: is it a quadratic-scaling routine?
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
       !CALL cp_dbcsr_print(no_quench(1),error=error)
       IF (almo_scf_env%nspins.gt.1) THEN
          DO ispin=2,almo_scf_env%nspins
             CALL cp_dbcsr_init(no_quench(ispin),error=error)
             CALL cp_dbcsr_create(no_quench(ispin),&
                     template=almo_scf_env%matrix_t(1),&
                     matrix_type=dbcsr_type_no_symmetry,error=error)
             CALL cp_dbcsr_copy(no_quench(ispin),no_quench(1),error=error)
          ENDDO
       ENDIF

    END SELECT

    SELECT CASE (almo_scf_env%deloc_method)
    CASE (almo_deloc_none,almo_deloc_full_scf)

       DO ispin=1,almo_scf_env%nspins
          CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin),&
                  almo_scf_env%matrix_t_blk(ispin),error=error)
       ENDDO

    CASE (almo_deloc_x,almo_deloc_xk,almo_deloc_x_then_full_scf)

       !!!! RZK-warning a whole class of delocalization methods
       !!!! are commented out at the moment because some of their
       !!!! routines have not been thoroughly tested.
       !!!! Some of them are not to be working correctly

       !!!! if we have virtuals pre-optimize and truncate them
       !!!IF (almo_scf_env%need_virtuals) THEN
       !!!   SELECT CASE (almo_scf_env%deloc_truncate_virt)
       !!!   CASE (virt_full)
       !!!      ! simply copy virtual orbitals from matrix_v_full_blk to matrix_v_blk
       !!!      DO ispin=1,almo_scf_env%nspins
       !!!         CALL cp_dbcsr_copy(almo_scf_env%matrix_v_blk(ispin),&
       !!!                 almo_scf_env%matrix_v_full_blk(ispin),error=error)
       !!!      ENDDO
       !!!   CASE (virt_number,virt_occ_size)
       !!!      CALL split_v_blk(almo_scf_env,error)
       !!!      !CALL truncate_subspace_v_blk(qs_env,almo_scf_env,error)
       !!!   CASE DEFAULT
       !!!      CPErrorMessage(cp_failure_level,routineP,"illegal method for virtual space truncation",error)
       !!!      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       !!!   END SELECT
       !!!ENDIF
       !!!CALL harris_foulkes_correction(qs_env,almo_scf_env,error)

       CALL almo_scf_xalmo_pcg(qs_env,almo_scf_env,&
               no_quench,almo_scf_env%matrix_t_blk,&
               almo_scf_env%matrix_t,&
               assume_t0_q0x=.TRUE.,perturbation_only=.TRUE.,&
               special_case=1,error=error)

    CASE (almo_deloc_qx)

       almo_scf_env%perturbative_delocalization=.TRUE.
       IF (almo_scf_env%almo_update_algorithm_q.eq.almo_scf_diag) THEN
          DO ispin=1,almo_scf_env%nspins
             CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin),&
                     almo_scf_env%matrix_t_blk(ispin),error=error)
          ENDDO
          CALL almo_scf_xalmo_eigensolver(qs_env,almo_scf_env,error)
       ELSE IF (almo_scf_env%almo_update_algorithm_q.eq.almo_scf_pcg) THEN
          CALL almo_scf_xalmo_pcg(qs_env,almo_scf_env,&
                  almo_scf_env%quench_t,almo_scf_env%matrix_t_blk,&
                  almo_scf_env%matrix_t,&
                  assume_t0_q0x=.TRUE.,perturbation_only=.TRUE.,&
                  special_case=-1,error=error)
       ENDIF

    CASE (almo_deloc_qscf)

       almo_scf_env%perturbative_delocalization=.FALSE.
       IF (almo_scf_env%almo_update_algorithm_q.eq.almo_scf_diag) THEN
          DO ispin=1,almo_scf_env%nspins
             CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin),&
                     almo_scf_env%matrix_t_blk(ispin),error=error)
          ENDDO
          CALL almo_scf_xalmo_eigensolver(qs_env,almo_scf_env,error)
       ELSE IF (almo_scf_env%almo_update_algorithm_q.eq.almo_scf_pcg) THEN
          CALL almo_scf_xalmo_pcg(qs_env,almo_scf_env,&
                  almo_scf_env%quench_t,almo_scf_env%matrix_t_blk,&
                  almo_scf_env%matrix_t,&
                  assume_t0_q0x=.TRUE.,perturbation_only=.FALSE.,&
                  special_case=-1,error=error)
       ENDIF

    CASE DEFAULT
       CPErrorMessage(cp_failure_level,routineP,"illegal post scf method",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    
    SELECT CASE (almo_scf_env%deloc_method)
    CASE (almo_deloc_full_scf,almo_deloc_x_then_full_scf)
       IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
          CPErrorMessage(cp_failure_level,routineP,"full scf is NYI for truncated virtual space",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       CALL almo_scf_xalmo_pcg(qs_env,almo_scf_env,&
               no_quench,almo_scf_env%matrix_t,&
               almo_scf_env%matrix_t,&
               assume_t0_q0x=.FALSE.,perturbation_only=.FALSE.,&
               special_case=1,error=error)
    END SELECT
    
    ! clean up 
    SELECT CASE (almo_scf_env%deloc_method)
    CASE (almo_deloc_x,almo_deloc_full_scf,almo_deloc_x_then_full_scf)
       DO ispin=1, almo_scf_env%nspins
          CALL cp_dbcsr_release(no_quench(ispin),error=error)
       ENDDO
       DEALLOCATE(no_quench)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE almo_scf_delocalization

! *****************************************************************************
!> \brief Split the matrix of virtual orbitals into two:
!>        retained orbs and discarded 
!> \par History
!>       2011.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE split_v_blk(almo_scf_env,error)

    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'split_v_blk', &
      routineP = moduleN//':'//routineN

    INTEGER :: discarded_v, handle, iblock_col, iblock_col_size, iblock_row, &
      iblock_row_size, ispin, retained_v
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: data_p, p_new_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)

    DO ispin=1,almo_scf_env%nspins

       CALL cp_dbcsr_work_create(almo_scf_env%matrix_v_blk(ispin),&
               work_mutable=.TRUE.,error=error)
       CALL cp_dbcsr_work_create(almo_scf_env%matrix_v_disc_blk(ispin),&
               work_mutable=.TRUE.,error=error)

       CALL cp_dbcsr_iterator_start(iter,almo_scf_env%matrix_v_full_blk(ispin))

       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

          CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
                  row_size=iblock_row_size,col_size=iblock_col_size)

          IF (iblock_row.ne.iblock_col) THEN
             CPErrorMessage(cp_failure_level,routineP,"off-diagonal block found",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          ENDIF

          retained_v=almo_scf_env%nvirt_of_domain(iblock_col,ispin)
          discarded_v=almo_scf_env%nvirt_disc_of_domain(iblock_col,ispin)
          CPPrecondition(retained_v.gt.0,cp_failure_level,routineP,error,failure)
          CPPrecondition(discarded_v.gt.0,cp_failure_level,routineP,error,failure)

          NULLIFY (p_new_block)
          CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_v_disc_blk(ispin),&
                  iblock_row,iblock_col,p_new_block)
          CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
          CPPrecondition((retained_v+discarded_v.eq.iblock_col_size),cp_failure_level,routineP,error,failure)
          p_new_block(:,:) = data_p(:,(retained_v+1):iblock_col_size)

          NULLIFY (p_new_block)
          CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_v_blk(ispin),&
                  iblock_row,iblock_col,p_new_block)
          CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
          p_new_block(:,:) = data_p(:,1:retained_v)

       ENDDO ! iterator
       CALL cp_dbcsr_iterator_stop(iter)

       CALL cp_dbcsr_finalize(almo_scf_env%matrix_v_blk(ispin),error=error)       
       CALL cp_dbcsr_finalize(almo_scf_env%matrix_v_disc_blk(ispin),error=error)

    ENDDO ! ispin
    
    CALL timestop(handle)

  END SUBROUTINE split_v_blk

! *****************************************************************************
!> \brief various methods for calculating the Harris-Foulkes correction
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE harris_foulkes_correction(qs_env,almo_scf_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_foulkes_correction', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: cayley_transform = 1, &
                                                dm_ls_step = 2

    INTEGER :: algorithm_id, handle, handle1, handle2, handle3, handle4, &
      handle5, handle6, handle7, handle8, ispin, iteration, n, nmins, nspin, &
      opt_k_max_iter, outer_opt_k_iteration, outer_opt_k_max_iter, unit_nr
    INTEGER, DIMENSION(1)                    :: fake, nelectron_spin_real
    LOGICAL :: converged, failure, line_search, md_in_k_space, &
      outer_opt_k_prepare_to_exit, prepare_to_exit, reset_conjugator, &
      reset_step_size, use_cubic_approximation, use_quadratic_approximation
    REAL(KIND=dp) :: aa, bb, beta, conjugacy_error, &
      conjugacy_error_threshold, delta_obj_function, denom, &
      energy_correction_final, frob_matrix, frob_matrix_base, fun0, fun1, &
      gfun0, gfun1, grad_norm, grad_norm_frob, kappa, kin_energy, &
      line_search_error, line_search_error_threshold, num_threshold, numer, &
      obj_function, quadratic_approx_error, quadratic_approx_error_threshold, &
      safety_multiplier, spin_factor, step_size, step_size_quadratic_approx, &
      step_size_quadratic_approx2, t1, t1a, t1cholesky, t2, t2a, t2cholesky, &
      tau, time_step, x_opt_eps_adaptive, x_opt_eps_adaptive_factor
    REAL(KIND=dp), DIMENSION(1)              :: local_mu
    REAL(KIND=dp), DIMENSION(2)              :: energy_correction
    REAL(KIND=dp), DIMENSION(3)              :: minima
    !TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type) :: grad, k_vd_index_down, k_vr_index_down, &
      matrix_k_central, matrix_tmp1, matrix_tmp2, prec, prev_grad, &
      prev_minus_prec_grad, prev_step, sigma_oo_curr, sigma_oo_curr_inv, &
      sigma_vv_sqrt, sigma_vv_sqrt_guess, sigma_vv_sqrt_inv, &
      sigma_vv_sqrt_inv_guess, step, t_curr, tmp1_n_vr, tmp2_n_o, tmp3_vd_vr, &
      tmp4_o_vr, tmp_k_blk, vd_fixed, vd_index_sqrt, vd_index_sqrt_inv, &
      velocity, vr_fixed, vr_index_sqrt, vr_index_sqrt_inv
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: matrix_p_almo_scf_converged
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(ct_step_env_type)                   :: ct_step_env

!TYPE(cp_dbcsr_type) :: tmp11,tmp22,tmp33
!REAL(kind=dp)                            :: k_var1, k_var2
!TYPE(cp_dbcsr_type)                      :: fake_step 
!
!TYPE(cp_dbcsr_type)                      :: sigma_dr, sigma_dr2, sigma_rr, sigma_rr2
!TYPE(cp_dbcsr_type)                      :: fake_a,fake_b,fake_k0
!INTEGER                                  :: retained_v,discarded_v,i_row,j_col
!TYPE(cp_dbcsr_type) :: matrix_rst0, matrix_rst1, matrix_rst2, ss_vv
!REAL(KIND=dp)       :: filter_memorize, init_filter, occ_vv
!INTEGER             :: ppp

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    nspin=almo_scf_env%nspins
    energy_correction_final=0.0_dp
    IF (nspin.eq.1) THEN
       spin_factor = 2.0_dp
    ELSE
       spin_factor = 1.0_dp
    ENDIF


    IF (almo_scf_env%deloc_use_occ_orbs) THEN
       algorithm_id=cayley_transform
    ELSE 
       algorithm_id=dm_ls_step
    ENDIF

    t1 = m_walltime()

    SELECT CASE (algorithm_id)
    CASE (cayley_transform)

       ! rescale density matrix by spin factor
       ! so the orbitals and density are consistent with each other
       IF (almo_scf_env%nspins == 1) THEN
        CALL cp_dbcsr_scale(almo_scf_env%matrix_p(1),1.0_dp/spin_factor,error=error)
       ENDIF
   
       ! transform matrix_t not matrix_t_blk (we might need ALMOs later)
       DO ispin=1,nspin
          
          CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin),&
                   almo_scf_env%matrix_t_blk(ispin),error=error) 
  
          ! obtain orthogonalization matrices for ALMOs
          ! RZK-warning - remove this sqrt(sigma) and inv(sqrt(sigma))
          ! ideally ALMO scf should use sigma and sigma_inv in
          ! the tensor_up_down representation
   
          IF (unit_nr>0) THEN
             WRITE(unit_nr,*) "sqrt and inv(sqrt) of MO overlap matrix"
          ENDIF
          CALL cp_dbcsr_init(almo_scf_env%matrix_sigma_sqrt(ispin),error=error)
          CALL cp_dbcsr_init(almo_scf_env%matrix_sigma_sqrt_inv(ispin),error=error)
          CALL cp_dbcsr_create(almo_scf_env%matrix_sigma_sqrt(ispin),&
                               template=almo_scf_env%matrix_sigma(ispin),&
                               matrix_type=dbcsr_type_no_symmetry,error=error) 
          CALL cp_dbcsr_create(almo_scf_env%matrix_sigma_sqrt_inv(ispin),&
                               template=almo_scf_env%matrix_sigma(ispin),&
                               matrix_type=dbcsr_type_no_symmetry,error=error) 
      
          CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_sigma_sqrt(ispin),&
                                         almo_scf_env%matrix_sigma_sqrt_inv(ispin),&
                                         almo_scf_env%matrix_sigma(ispin),&
                                         threshold=almo_scf_env%eps_filter,&
                                         order=almo_scf_env%order_lanczos,&
                                         eps_lanczos=almo_scf_env%eps_lanczos,&
                                         max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
                                         error=error)
      
          IF (safe_mode) THEN
             CALL cp_dbcsr_init(matrix_tmp1,error=error)
             CALL cp_dbcsr_create(matrix_tmp1,template=almo_scf_env%matrix_sigma(ispin),&
                                  matrix_type=dbcsr_type_no_symmetry,error=error) 
             CALL cp_dbcsr_init(matrix_tmp2,error=error)
             CALL cp_dbcsr_create(matrix_tmp2,template=almo_scf_env%matrix_sigma(ispin),&
                                  matrix_type=dbcsr_type_no_symmetry,error=error) 
      
             CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_sigma_sqrt_inv(ispin),&
                                    almo_scf_env%matrix_sigma(ispin),&
                                    0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter,&
                                    error=error)
             CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_tmp1,&
                                    almo_scf_env%matrix_sigma_sqrt_inv(ispin),&
                                    0.0_dp,matrix_tmp2,filter_eps=almo_scf_env%eps_filter,&
                                    error=error)
      
             frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
             CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp,error=error)
             frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
             IF (unit_nr>0) THEN
                WRITE(unit_nr,*) "Error for (inv(sqrt(SIG))*SIG*inv(sqrt(SIG))-I)",frob_matrix/frob_matrix_base
             ENDIF
      
             CALL cp_dbcsr_release(matrix_tmp1,error=error) 
             CALL cp_dbcsr_release(matrix_tmp2,error=error) 
          ENDIF
       ENDDO

       IF (almo_scf_env%almo_update_algorithm.eq.almo_scf_diag) THEN
  
          DO ispin=1,nspin
   
             t1a = m_walltime()

             line_search_error_threshold=almo_scf_env%real01
             conjugacy_error_threshold=almo_scf_env%real02
             quadratic_approx_error_threshold=almo_scf_env%real03
             x_opt_eps_adaptive_factor=almo_scf_env%real04

             !! the outer loop for k optimization
             outer_opt_k_max_iter=almo_scf_env%opt_k_outer_max_iter
             outer_opt_k_prepare_to_exit=.FALSE.
             outer_opt_k_iteration=0
             grad_norm=0.0_dp
             grad_norm_frob=0.0_dp
             CALL cp_dbcsr_set(almo_scf_env%matrix_x(ispin),0.0_dp,error=error)
             IF (almo_scf_env%deloc_truncate_virt.eq.virt_full) outer_opt_k_max_iter=0
             
             DO

                ! obtain proper retained virtuals (1-R)|ALMO_vr>
                CALL apply_projector(psi_in=almo_scf_env%matrix_v_blk(ispin),&
                        psi_out=almo_scf_env%matrix_v(ispin),&
                        psi_projector=almo_scf_env%matrix_t_blk(ispin),&
                        metric=almo_scf_env%matrix_s(1),&
                        project_out=.TRUE.,&
                        psi_projector_orthogonal=.FALSE.,&
                        proj_in_template=almo_scf_env%matrix_ov(ispin),&
                        eps_filter=almo_scf_env%eps_filter,&
                        sig_inv_projector=almo_scf_env%matrix_sigma_inv(ispin),&
                        !sig_inv_template=almo_scf_env%matrix_sigma_inv(ispin),&
                        error=error)

                ! save initial retained virtuals 
                CALL cp_dbcsr_init(vr_fixed,error=error)
                CALL cp_dbcsr_create(vr_fixed,&
                        template=almo_scf_env%matrix_v(ispin),error=error)
                CALL cp_dbcsr_copy(vr_fixed,almo_scf_env%matrix_v(ispin),&
                        error=error)

                ! init matrices common for optimized and non-optimized virts
                CALL cp_dbcsr_init(sigma_vv_sqrt,error=error)
                CALL cp_dbcsr_init(sigma_vv_sqrt_inv,error=error)
                CALL cp_dbcsr_create(sigma_vv_sqrt,&
                        template=almo_scf_env%matrix_sigma_vv(ispin),&
                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                CALL cp_dbcsr_create(sigma_vv_sqrt_inv,&
                        template=almo_scf_env%matrix_sigma_vv(ispin),&
                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                CALL cp_dbcsr_init(sigma_vv_sqrt_inv_guess,error=error)
                CALL cp_dbcsr_init(sigma_vv_sqrt_guess,error=error)
                CALL cp_dbcsr_create(sigma_vv_sqrt_inv_guess,&
                        template=almo_scf_env%matrix_sigma_vv(ispin),&
                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                CALL cp_dbcsr_create(sigma_vv_sqrt_guess,&
                        template=almo_scf_env%matrix_sigma_vv(ispin),&
                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                CALL cp_dbcsr_set(sigma_vv_sqrt_guess,0.0_dp,error=error)
                CALL cp_dbcsr_add_on_diag(sigma_vv_sqrt_guess,1.0_dp,error=error)
                CALL cp_dbcsr_filter(sigma_vv_sqrt_guess,almo_scf_env%eps_filter,&
                        error=error)
                CALL cp_dbcsr_set(sigma_vv_sqrt_inv_guess,0.0_dp,error=error)
                CALL cp_dbcsr_add_on_diag(sigma_vv_sqrt_inv_guess,1.0_dp,error=error)
                CALL cp_dbcsr_filter(sigma_vv_sqrt_inv_guess,almo_scf_env%eps_filter,&
                        error=error)
      
                ! do things required to optimize virtuals
                IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
                
                   ! project retained virtuals out of discarded block-by-block
                   ! (1-Q^VR_ALMO)|ALMO_vd>
                   ! this is probably not necessary, do it just to be safe
                   !CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin),&
                   !        psi_out=almo_scf_env%matrix_v_disc(ispin),&
                   !        psi_projector=almo_scf_env%matrix_v_blk(ispin),&
                   !        metric=almo_scf_env%matrix_s_blk(1),&
                   !        project_out=.TRUE.,&
                   !        psi_projector_orthogonal=.FALSE.,&
                   !        proj_in_template=almo_scf_env%matrix_k_tr(ispin),&
                   !        eps_filter=almo_scf_env%eps_filter,&
                   !        sig_inv_template=almo_scf_env%matrix_sigma_vv(ispin),&
                   !        error=error)
                   !CALL cp_dbcsr_copy(almo_scf_env%matrix_v_disc_blk(ispin),&
                   !        almo_scf_env%matrix_v_disc(ispin),error=error)

                   ! construct discarded virtuals (1-R)|ALMO_vd>
                   CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin),&
                           psi_out=almo_scf_env%matrix_v_disc(ispin),&
                           psi_projector=almo_scf_env%matrix_t_blk(ispin),&
                           metric=almo_scf_env%matrix_s(1),&
                           project_out=.TRUE.,&
                           psi_projector_orthogonal=.FALSE.,&
                           proj_in_template=almo_scf_env%matrix_ov_disc(ispin),&
                           eps_filter=almo_scf_env%eps_filter,&
                           sig_inv_projector=almo_scf_env%matrix_sigma_inv(ispin),&
                           !sig_inv_template=almo_scf_env%matrix_sigma_inv(ispin),&
                           error=error)
                   
                   ! save initial discarded
                   CALL cp_dbcsr_init(vd_fixed,error=error)
                   CALL cp_dbcsr_create(vd_fixed,&
                           template=almo_scf_env%matrix_v_disc(ispin),error=error)
                   CALL cp_dbcsr_copy(vd_fixed,almo_scf_env%matrix_v_disc(ispin),&
                           error=error)

                   !! create the down metric in the retained k-subspace
                   CALL cp_dbcsr_init(k_vr_index_down,error=error)
                   CALL cp_dbcsr_create(k_vr_index_down,&
                           template=almo_scf_env%matrix_sigma_vv_blk(ispin),&
                           matrix_type=dbcsr_type_no_symmetry,error=error)
                   !CALL cp_dbcsr_copy(k_vr_index_down,&
                   !        almo_scf_env%matrix_sigma_vv_blk(ispin),error=error)

                   !CALL get_overlap(bra=almo_scf_env%matrix_v_blk(ispin),&
                   !        ket=almo_scf_env%matrix_v_blk(ispin),&
                   !        overlap=k_vr_index_down,&
                   !        metric=almo_scf_env%matrix_s_blk(1),&
                   !        retain_overlap_sparsity=.FALSE.,&
                   !        eps_filter=almo_scf_env%eps_filter,&
                   !        error=error)

                   !! create the up metric in the discarded k-subspace
                   CALL cp_dbcsr_init(k_vd_index_down,error=error)
                   CALL cp_dbcsr_create(k_vd_index_down,&
                           template=almo_scf_env%matrix_vv_disc_blk(ispin),&
                           matrix_type=dbcsr_type_no_symmetry,error=error)
                   !CALL cp_dbcsr_init(k_vd_index_up,error=error)
                   !CALL cp_dbcsr_create(k_vd_index_up,&
                   !        template=almo_scf_env%matrix_vv_disc_blk(ispin),&
                   !        matrix_type=dbcsr_type_no_symmetry,error=error)
                   !CALL cp_dbcsr_copy(k_vd_index_down,&
                   !        almo_scf_env%matrix_vv_disc_blk(ispin),error=error)

                   !CALL get_overlap(bra=almo_scf_env%matrix_v_disc_blk(ispin),&
                   !        ket=almo_scf_env%matrix_v_disc_blk(ispin),&
                   !        overlap=k_vd_index_down,&
                   !        metric=almo_scf_env%matrix_s_blk(1),&
                   !        retain_overlap_sparsity=.FALSE.,&
                   !        eps_filter=almo_scf_env%eps_filter,&
                   !        error=error)

                   !IF (unit_nr>0) THEN
                   !   WRITE(unit_nr,*) "Inverting blocked overlap matrix of discarded virtuals"
                   !ENDIF
                   !CALL invert_Hotelling(k_vd_index_up,&
                   !        k_vd_index_down,&
                   !        almo_scf_env%eps_filter,&
                   !        error)
                   !IF (safe_mode) THEN
                   !   CALL cp_dbcsr_init(matrix_tmp1,error=error)
                   !   CALL cp_dbcsr_create(matrix_tmp1,template=k_vd_index_down,&
                   !                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                   !   CALL cp_dbcsr_multiply("N","N",1.0_dp,k_vd_index_up,&
                   !                          k_vd_index_down,&
                   !                          0.0_dp, matrix_tmp1,&
                   !                          filter_eps=almo_scf_env%eps_filter,error=error)
                   !   frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
                   !   CALL cp_dbcsr_add_on_diag(matrix_tmp1,-1.0_dp,error=error)
                   !   frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
                   !   IF (unit_nr>0) THEN
                   !      WRITE(unit_nr,*) "Error for (inv(SIG)*SIG-I)",&
                   !            frob_matrix/frob_matrix_base
                   !   ENDIF
                   !   CALL cp_dbcsr_release(matrix_tmp1,error=error)
                   !ENDIF

                   ! init matrices necessary for optimization of truncated virts
                   ! init blocked gradient before setting K to zero
                   ! otherwise the block structure might be lost
                   CALL cp_dbcsr_init(grad,error=error)
                   CALL cp_dbcsr_create(grad,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error)
                   CALL cp_dbcsr_copy(grad,almo_scf_env%matrix_k_blk(ispin),error=error)

                   ! init MD in the k-space
                   md_in_k_space=almo_scf_env%logical01
                   IF (md_in_k_space) THEN
                      CALL cp_dbcsr_init(velocity,error=error)
                      CALL cp_dbcsr_create(velocity,&
                              template=almo_scf_env%matrix_k_blk(ispin),error=error)
                      CALL cp_dbcsr_copy(velocity,almo_scf_env%matrix_k_blk(ispin),error=error)
                      CALL cp_dbcsr_set(velocity,0.0_dp,error=error)
                      time_step=almo_scf_env%opt_k_trial_step_size
                   ENDIF

                   CALL cp_dbcsr_init(prev_step,error=error)
                   CALL cp_dbcsr_create(prev_step,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error)

                   CALL cp_dbcsr_init(prev_minus_prec_grad,error=error)
                   CALL cp_dbcsr_create(prev_minus_prec_grad,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error)

                   ! initialize diagonal blocks of the preconditioner to 1.0_dp
                   CALL cp_dbcsr_init(prec,error=error)
                   CALL cp_dbcsr_create(prec,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error)
                   CALL cp_dbcsr_copy(prec,almo_scf_env%matrix_k_blk(ispin),error=error)
                   CALL cp_dbcsr_set(prec,1.0_dp,error=error)
                   
                   ! generate initial K (extrapolate if previous values are available)
                   CALL cp_dbcsr_set(almo_scf_env%matrix_k_blk(ispin),0.0_dp,error=error)
                   ! matrix_k_central stores current k because matrix_k_blk is updated 
                   ! during linear search
                   CALL cp_dbcsr_init(matrix_k_central,error=error)
                   CALL cp_dbcsr_create(matrix_k_central,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error)
                   CALL cp_dbcsr_copy(matrix_k_central,&
                           almo_scf_env%matrix_k_blk(ispin),error=error)
                   CALL cp_dbcsr_init(tmp_k_blk,error=error)
                   CALL cp_dbcsr_create(tmp_k_blk,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error)
                   CALL cp_dbcsr_init(step,error=error)
                   CALL cp_dbcsr_create(step,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error)
                   CALL cp_dbcsr_set(step,0.0_dp,error=error)
                   CALL cp_dbcsr_init(t_curr,error=error)
                   CALL cp_dbcsr_create(t_curr,&
                           template=almo_scf_env%matrix_t(ispin),error=error)
                   CALL cp_dbcsr_init(sigma_oo_curr,error=error)
                   CALL cp_dbcsr_init(sigma_oo_curr_inv,error=error)
                   CALL cp_dbcsr_create(sigma_oo_curr,&
                           template=almo_scf_env%matrix_sigma(ispin),&
                           matrix_type=dbcsr_type_no_symmetry,error=error) 
                   CALL cp_dbcsr_create(sigma_oo_curr_inv,&
                           template=almo_scf_env%matrix_sigma(ispin),&
                           matrix_type=dbcsr_type_no_symmetry,error=error) 
                   CALL cp_dbcsr_init(tmp1_n_vr,error=error)
                   CALL cp_dbcsr_create(tmp1_n_vr,&
                           template=almo_scf_env%matrix_v(ispin),error=error)
                   CALL cp_dbcsr_init(tmp3_vd_vr,error=error)
                   CALL cp_dbcsr_create(tmp3_vd_vr,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error) 
                   CALL cp_dbcsr_init(tmp2_n_o,error=error)
                   CALL cp_dbcsr_create(tmp2_n_o,&
                           template=almo_scf_env%matrix_t(ispin),error=error)
                   CALL cp_dbcsr_init(tmp4_o_vr,error=error)
                   CALL cp_dbcsr_create(tmp4_o_vr,&
                           template=almo_scf_env%matrix_ov(ispin),error=error)
                   CALL cp_dbcsr_init(prev_grad,error=error)
                   CALL cp_dbcsr_create(prev_grad,&
                           template=almo_scf_env%matrix_k_blk(ispin),error=error)
                   CALL cp_dbcsr_set(prev_grad,0.0_dp,error=error)

                   !CALL cp_dbcsr_init(sigma_oo_guess,error=error)
                   !CALL cp_dbcsr_create(sigma_oo_guess,&
                   !        template=almo_scf_env%matrix_sigma(ispin),&
                   !        matrix_type=dbcsr_type_no_symmetry,error=error) 
                   !CALL cp_dbcsr_set(sigma_oo_guess,0.0_dp,error=error)
                   !CALL cp_dbcsr_add_on_diag(sigma_oo_guess,1.0_dp,error=error)
                   !CALL cp_dbcsr_filter(sigma_oo_guess,almo_scf_env%eps_filter,error=error)
                   !CALL cp_dbcsr_print(sigma_oo_guess,error=error)

                ENDIF ! done constructing discarded virtuals

                ! init variables 
                opt_k_max_iter=almo_scf_env%opt_k_max_iter
                iteration=0
                converged=.FALSE.
                prepare_to_exit=.FALSE.
                beta=0.0_dp
                line_search=.FALSE.
                obj_function=0.0_dp
                conjugacy_error=0.0_dp
                line_search_error=0.0_dp
                fun0=0.0_dp
                fun1=0.0_dp
                gfun0=0.0_dp
                gfun1=0.0_dp
                step_size_quadratic_approx=0.0_dp
                reset_step_size=.TRUE.
                IF (almo_scf_env%deloc_truncate_virt.eq.virt_full) opt_k_max_iter=0

                ! start cg iterations to optimize matrix_k_blk
                DO

                   CALL timeset('k_opt_vr',handle1)

                   IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN

                      ! construct k-excited virtuals
                      CALL cp_dbcsr_multiply("N","N",1.0_dp,vd_fixed,&
                              almo_scf_env%matrix_k_blk(ispin),&
                              0.0_dp,almo_scf_env%matrix_v(ispin),&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                      CALL cp_dbcsr_add(almo_scf_env%matrix_v(ispin),vr_fixed,&
                              +1.0_dp,+1.0_dp,error=error)
                   ENDIF
                
                   ! decompose the overlap matrix of the current retained orbitals
                   !IF (unit_nr>0) THEN
                   !   WRITE(unit_nr,*) "decompose the active VV overlap matrix"
                   !ENDIF
                   CALL get_overlap(bra=almo_scf_env%matrix_v(ispin),&
                           ket=almo_scf_env%matrix_v(ispin),&
                           overlap=almo_scf_env%matrix_sigma_vv(ispin),&
                           metric=almo_scf_env%matrix_s(1),&
                           retain_overlap_sparsity=.FALSE.,&
                           eps_filter=almo_scf_env%eps_filter,&
                           error=error)
                   ! use either cholesky or sqrt
                   !! RZK-warning: strangely, cholesky does not work with k-optimization
                   IF (almo_scf_env%deloc_truncate_virt.eq.virt_full) THEN
                      CALL timeset('cholesky',handle2)
                      t1cholesky = m_walltime()

                      ! re-create sigma_vv_sqrt because desymmetrize is buggy -
                      ! it will create multiple copies of blocks
                      CALL cp_dbcsr_create(sigma_vv_sqrt,&
                              template=almo_scf_env%matrix_sigma_vv(ispin),&
                              matrix_type=dbcsr_type_no_symmetry,error=error) 
                      CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_vv(ispin),&
                              sigma_vv_sqrt,error=error)
                      CALL cp_dbcsr_cholesky_decompose(sigma_vv_sqrt,&
                              para_env=almo_scf_env%para_env,&
                              blacs_env=almo_scf_env%blacs_env,error=error)
                      CALL cp_dbcsr_triu(sigma_vv_sqrt,error=error)
                      CALL cp_dbcsr_filter(sigma_vv_sqrt,almo_scf_env%eps_filter,error=error)
                      ! apply SOLVE to compute U^(-1) : U*U^(-1)=I
                      CALL cp_dbcsr_get_info(sigma_vv_sqrt, nfullrows_total=n )
                      CALL cp_dbcsr_init(matrix_tmp1,error=error)
                      CALL cp_dbcsr_create(matrix_tmp1,template=almo_scf_env%matrix_sigma_vv(ispin),&
                              matrix_type=dbcsr_type_no_symmetry,error=error)
                      CALL cp_dbcsr_set(matrix_tmp1,0.0_dp,error=error)
                      CALL cp_dbcsr_add_on_diag(matrix_tmp1,1.0_dp,error=error)
                      CALL cp_dbcsr_cholesky_restore(matrix_tmp1,n,sigma_vv_sqrt,&
                              sigma_vv_sqrt_inv,op="SOLVE",pos="RIGHT",&
                              para_env=almo_scf_env%para_env,&
                              blacs_env=almo_scf_env%blacs_env,error=error)
                      CALL cp_dbcsr_filter(sigma_vv_sqrt_inv,almo_scf_env%eps_filter,error=error)
                      CALL cp_dbcsr_release(matrix_tmp1,error=error)
                      IF (safe_mode) THEN
                         CALL cp_dbcsr_create(matrix_tmp1, template=almo_scf_env%matrix_sigma_vv(ispin),&
                             matrix_type=dbcsr_type_no_symmetry,error=error)
                         CALL cp_dbcsr_desymmetrize(almo_scf_env%matrix_sigma_vv(ispin),&
                                 matrix_tmp1,error=error)
                         CALL cp_dbcsr_multiply("T","N",1.0_dp,sigma_vv_sqrt,&
                                 sigma_vv_sqrt,&
                                 -1.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter,&
                                 error=error)
                         frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
                         CALL cp_dbcsr_add_on_diag(matrix_tmp1,1.0_dp,error=error)
                         frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
                         IF (unit_nr>0) THEN
                            WRITE(unit_nr,*) "Error for ( U^T * U - Sig )",&
                               frob_matrix/frob_matrix_base
                         ENDIF
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,sigma_vv_sqrt_inv,&
                                 sigma_vv_sqrt,&
                                 0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter,&
                                 error=error)
                         frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
                         CALL cp_dbcsr_add_on_diag(matrix_tmp1,-1.0_dp,error=error)
                         frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
                         IF (unit_nr>0) THEN
                            WRITE(unit_nr,*) "Error for ( inv(U) * U - I )",&
                               frob_matrix/frob_matrix_base
                         ENDIF
                         CALL cp_dbcsr_release(matrix_tmp1,error=error)
                      ENDIF ! safe_mode
                      t2cholesky = m_walltime()
                      IF (unit_nr>0) THEN
                         WRITE(unit_nr,*) "Cholesky+inverse wall-time: ",t2cholesky-t1cholesky
                      ENDIF
                      CALL timestop(handle2)
                   ELSE
                      CALL matrix_sqrt_Newton_Schulz(sigma_vv_sqrt,&
                              sigma_vv_sqrt_inv,&
                              almo_scf_env%matrix_sigma_vv(ispin),&
                              !matrix_sqrt_inv_guess=sigma_vv_sqrt_inv_guess,&
                              !matrix_sqrt_guess=sigma_vv_sqrt_guess,&
                              threshold=almo_scf_env%eps_filter,&
                              order=almo_scf_env%order_lanczos,&
                              eps_lanczos=almo_scf_env%eps_lanczos,&
                              max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
                              error=error)
                      CALL cp_dbcsr_copy(sigma_vv_sqrt_inv_guess,sigma_vv_sqrt_inv,error=error)
                      CALL cp_dbcsr_copy(sigma_vv_sqrt_guess,sigma_vv_sqrt,error=error)
                      IF (safe_mode) THEN
                         CALL cp_dbcsr_init(matrix_tmp1,error=error)
                         CALL cp_dbcsr_create(matrix_tmp1,template=almo_scf_env%matrix_sigma_vv(ispin),&
                                              matrix_type=dbcsr_type_no_symmetry,error=error) 
                         CALL cp_dbcsr_init(matrix_tmp2,error=error)
                         CALL cp_dbcsr_create(matrix_tmp2,template=almo_scf_env%matrix_sigma_vv(ispin),&
                                              matrix_type=dbcsr_type_no_symmetry,error=error) 
       
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,sigma_vv_sqrt_inv,&
                                                almo_scf_env%matrix_sigma_vv(ispin),&
                                                0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter,&
                                                error=error)
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_tmp1,&
                                                sigma_vv_sqrt_inv,&
                                                0.0_dp,matrix_tmp2,filter_eps=almo_scf_env%eps_filter,&
                                                error=error)
       
                         frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
                         CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp,error=error)
                         frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
                         IF (unit_nr>0) THEN
                            WRITE(unit_nr,*) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)",&
                               frob_matrix/frob_matrix_base
                         ENDIF
       
                         CALL cp_dbcsr_release(matrix_tmp1,error=error) 
                         CALL cp_dbcsr_release(matrix_tmp2,error=error) 
                      ENDIF
                   ENDIF   
                   CALL timestop(handle1)
    
                   ! compute excitation amplitudes (to the current set of retained virtuals)
                   ! set convergence criterion for x-optimization
                   IF ((iteration.eq.0).AND.(.NOT.line_search).AND.&
                      (outer_opt_k_iteration.eq.0)) THEN
                      x_opt_eps_adaptive=&
                         almo_scf_env%deloc_cayley_eps_convergence
                   ELSE
                      x_opt_eps_adaptive=&
                         MAX(ABS(almo_scf_env%deloc_cayley_eps_convergence),&
                         ABS(x_opt_eps_adaptive_factor*grad_norm))
                   ENDIF
                   CALL ct_step_env_init(ct_step_env,error)
                   CALL ct_step_env_set(ct_step_env,&
                           para_env=almo_scf_env%para_env,&
                           blacs_env=almo_scf_env%blacs_env,&
                           use_occ_orbs=.TRUE.,&
                           use_virt_orbs=.TRUE.,&
                           occ_orbs_orthogonal=.FALSE.,&
                           virt_orbs_orthogonal=.FALSE.,&
                           pp_preconditioner_full=almo_scf_env%deloc_cayley_occ_precond,&
                           qq_preconditioner_full=almo_scf_env%deloc_cayley_vir_precond,&
                           tensor_type=almo_scf_env%deloc_cayley_tensor_type,&
                           neglect_quadratic_term=almo_scf_env%deloc_cayley_linear,&
                           conjugator=almo_scf_env%deloc_cayley_conjugator,&
                           max_iter=almo_scf_env%deloc_cayley_max_iter,&
                           calculate_energy_corr=.TRUE.,&
                           update_p=.FALSE.,&
                           update_q=.FALSE.,&
                           eps_convergence=x_opt_eps_adaptive,&
                           eps_filter=almo_scf_env%eps_filter,&
                           !nspins=1,&
                           q_index_up=sigma_vv_sqrt_inv,&
                           q_index_down=sigma_vv_sqrt,&
                           p_index_up=almo_scf_env%matrix_sigma_sqrt_inv(ispin),&
                           p_index_down=almo_scf_env%matrix_sigma_sqrt(ispin),&
                           matrix_ks=almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                           matrix_t=almo_scf_env%matrix_t(ispin),&
                           matrix_qp_template=almo_scf_env%matrix_vo(ispin),&
                           matrix_pq_template=almo_scf_env%matrix_ov(ispin),&
                           matrix_v=almo_scf_env%matrix_v(ispin),&
                           matrix_x_guess=almo_scf_env%matrix_x(ispin),&
                           error=error)
                   ! perform calculations
                   CALL ct_step_execute(ct_step_env,error)
                   ! get the energy correction
                   CALL ct_step_env_get(ct_step_env,&
                           energy_correction=energy_correction(ispin),&
                           copy_matrix_x=almo_scf_env%matrix_x(ispin),&
                           error=error)
                   CALL ct_step_env_clean(ct_step_env,error)
                   ! RZK-warning matrix_x is being transformed 
                   ! back and forth between orth and up_down representations
                   energy_correction(1)=energy_correction(1)*spin_factor
   
                   IF (opt_k_max_iter.ne.0) THEN
   
                      CALL timeset('k_opt_t_curr',handle3)
                      
                      ! construct current occupied orbitals T_blk + V_r*X
                      CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                              almo_scf_env%matrix_v(ispin),&
                              almo_scf_env%matrix_x(ispin),&
                              0.0_dp,t_curr,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                      CALL cp_dbcsr_add(t_curr,almo_scf_env%matrix_t_blk(ispin),&
                              +1.0_dp,+1.0_dp,error=error)
   
                      ! calculate current occupied overlap
                      !IF (unit_nr>0) THEN
                      !   WRITE(unit_nr,*) "Inverting current occ overlap matrix"
                      !ENDIF
                      CALL get_overlap(bra=t_curr,&
                              ket=t_curr,&
                              overlap=sigma_oo_curr,&
                              metric=almo_scf_env%matrix_s(1),&
                              retain_overlap_sparsity=.FALSE.,&
                              eps_filter=almo_scf_env%eps_filter,&
                              error=error)
                      IF (iteration.eq.0) THEN
                         CALL invert_Hotelling(sigma_oo_curr_inv,&
                                 sigma_oo_curr,&
                                 threshold=almo_scf_env%eps_filter,&
                                 use_inv_as_guess=.FALSE.,&
                                 error=error)
                      ELSE
                         CALL invert_Hotelling(sigma_oo_curr_inv,&
                                 sigma_oo_curr,&
                                 threshold=almo_scf_env%eps_filter,&
                                 use_inv_as_guess=.TRUE.,&
                                 error=error)
                         !CALL cp_dbcsr_copy(sigma_oo_guess,sigma_oo_curr_inv,error=error)
                      ENDIF
                      IF (safe_mode) THEN
                         CALL cp_dbcsr_init(matrix_tmp1,error=error)
                         CALL cp_dbcsr_create(matrix_tmp1,template=sigma_oo_curr,&
                                              matrix_type=dbcsr_type_no_symmetry,error=error) 
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,sigma_oo_curr,&
                                                sigma_oo_curr_inv,&
                                                0.0_dp, matrix_tmp1,&
                                                filter_eps=almo_scf_env%eps_filter,error=error)
                         frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
                         CALL cp_dbcsr_add_on_diag(matrix_tmp1,-1.0_dp,error=error)
                         frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
                         !CALL cp_dbcsr_filter(matrix_tmp1,almo_scf_env%eps_filter,error=error)
                         !CALL cp_dbcsr_print(matrix_tmp1,error=error)
                         IF (unit_nr>0) THEN
                            WRITE(unit_nr,*) "Error for (SIG*inv(SIG)-I)",&
                                  frob_matrix/frob_matrix_base, frob_matrix_base
                         ENDIF
                         CALL cp_dbcsr_release(matrix_tmp1,error=error)
                      ENDIF
                      IF (safe_mode) THEN
                         CALL cp_dbcsr_init(matrix_tmp1,error=error)
                         CALL cp_dbcsr_create(matrix_tmp1,template=sigma_oo_curr,&
                                              matrix_type=dbcsr_type_no_symmetry,error=error) 
                         CALL cp_dbcsr_multiply("N","N",1.0_dp,sigma_oo_curr_inv,&
                                                sigma_oo_curr,&
                                                0.0_dp, matrix_tmp1,&
                                                filter_eps=almo_scf_env%eps_filter,error=error)
                         frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
                         CALL cp_dbcsr_add_on_diag(matrix_tmp1,-1.0_dp,error=error)
                         frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
                         !CALL cp_dbcsr_filter(matrix_tmp1,almo_scf_env%eps_filter,error=error)
                         !CALL cp_dbcsr_print(matrix_tmp1,error=error)
                         IF (unit_nr>0) THEN
                            WRITE(unit_nr,*) "Error for (inv(SIG)*SIG-I)",&
                                  frob_matrix/frob_matrix_base, frob_matrix_base
                         ENDIF
                         CALL cp_dbcsr_release(matrix_tmp1,error=error)
                      ENDIF
   
                      CALL timestop(handle3)
                      CALL timeset('k_opt_vd',handle4)
   
                      ! construct current discarded virtuals:
                      ! (1-R_curr)(1-Q^VR_curr)|ALMO_vd_basis> =
                      ! = (1-Q^VR_curr)|ALMO_vd_basis>
                      ! use sigma_vv_sqrt to store the inverse of the overlap
                      ! sigma_vv_inv is computed from sqrt/cholesky
                      CALL cp_dbcsr_multiply("N","T",1.0_dp,&
                              sigma_vv_sqrt_inv,&
                              sigma_vv_sqrt_inv,&
                              0.0_dp,sigma_vv_sqrt,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                      CALL apply_projector(psi_in=almo_scf_env%matrix_v_disc_blk(ispin),&
                              psi_out=almo_scf_env%matrix_v_disc(ispin),&
                              psi_projector=almo_scf_env%matrix_v(ispin),&
                              metric=almo_scf_env%matrix_s(1),&
                              project_out=.FALSE.,&
                              psi_projector_orthogonal=.FALSE.,&
                              proj_in_template=almo_scf_env%matrix_k_tr(ispin),&
                              eps_filter=almo_scf_env%eps_filter,&
                              sig_inv_projector=sigma_vv_sqrt,&
                              !sig_inv_template=almo_scf_env%matrix_sigma_vv(ispin),&
                              error=error)
                      CALL cp_dbcsr_add(almo_scf_env%matrix_v_disc(ispin),&
                              vd_fixed,-1.0_dp,+1.0_dp,error=error)
                      
                      CALL timestop(handle4)
                      CALL timeset('k_opt_grad',handle5)
   
                      ! evaluate the gradient from the assembled components
                      ! grad_xx = c0 [ (Vd_curr^tr)*F*T_curr*sigma_oo_curr_inv*(X^tr)]_xx
                      ! save previous gradient to calculate conjugation coef
                      IF (line_search) THEN
                         CALL cp_dbcsr_copy(prev_grad,grad,error=error)
                      ENDIF
                      CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                              almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                              t_curr,&
                              0.0_dp,tmp2_n_o,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                      CALL cp_dbcsr_multiply("N","T",1.0_dp,&
                              sigma_oo_curr_inv,&
                              almo_scf_env%matrix_x(ispin),&
                              0.0_dp,tmp4_o_vr,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                      CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                              tmp2_n_o,&
                              tmp4_o_vr,&
                              0.0_dp,tmp1_n_vr,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                      CALL cp_dbcsr_multiply("T","N",2.0_dp*spin_factor,&
                              almo_scf_env%matrix_v_disc(ispin),&
                              tmp1_n_vr,&
                              0.0_dp,grad,&
                              retain_sparsity=.TRUE.,&
                              !filter_eps=almo_scf_env%eps_filter,&
                              error=error)
                      ! keep tmp2_n_o for the next step
                      ! keep tmp4_o_vr for the preconditioner
               
                      ! check convergence and other exit criteria
                      grad_norm_frob=cp_dbcsr_frobenius_norm(grad)
                      CALL cp_dbcsr_norm(grad, dbcsr_norm_maxabsnorm, norm_scalar=grad_norm, error=error)
                      converged=(grad_norm.lt.almo_scf_env%opt_k_eps_convergence)
                      IF (converged.OR.(iteration.ge.opt_k_max_iter)) THEN
                         prepare_to_exit=.TRUE.
                      ENDIF
                      CALL timestop(handle5)
   
   
                      IF (.NOT.prepare_to_exit) THEN
   
                         CALL timeset('k_opt_energy',handle6)
   
                         ! compute "energy" c0*Tr[sig_inv_oo*t*F*t]
                         CALL cp_dbcsr_multiply("T","N",spin_factor,&
                                 t_curr,&
                                 tmp2_n_o,&
                                 0.0_dp,sigma_oo_curr,&
                                 filter_eps=almo_scf_env%eps_filter,&
                                 error=error)
                         delta_obj_function=fun0
                         CALL cp_dbcsr_trace(sigma_oo_curr_inv,&
                                 sigma_oo_curr,obj_function,"T","N",&
                                 error=error)
                         delta_obj_function=obj_function-delta_obj_function
                         IF (line_search) THEN
                            fun1=obj_function
                         ELSE
                            fun0=obj_function
                         ENDIF
   
                         CALL timestop(handle6)

                         ! update the step direction
                         IF (.NOT.line_search) THEN
                            
                            CALL timeset('k_opt_step',handle7)
   
                            IF ((.NOT.md_in_k_space).AND.&
                               (iteration.ge.MAX(0,almo_scf_env%opt_k_prec_iter_start).AND.&
                               MOD(iteration-almo_scf_env%opt_k_prec_iter_start,&
                               almo_scf_env%opt_k_prec_iter_freq).eq.0)) THEN
                               
                            !IF ((iteration.eq.0).AND.(.NOT.md_in_k_space)) THEN
   
                               ! compute the preconditioner
                               IF (unit_nr>0) THEN
                                  WRITE(unit_nr,*) "Computing preconditioner"
                               ENDIF
                               !CALL opt_k_create_preconditioner(prec,&
                               !        almo_scf_env%matrix_v_disc(ispin),&
                               !        almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                               !        almo_scf_env%matrix_x(ispin),&
                               !        tmp4_o_vr,&
                               !        almo_scf_env%matrix_s(1),&
                               !        grad,&
                               !        !almo_scf_env%matrix_v_disc_blk(ispin),&
                               !        vd_fixed,&
                               !        t_curr,&
                               !        k_vd_index_up,&
                               !        k_vr_index_down,&
                               !        tmp1_n_vr,&
                               !        spin_factor,&
                               !        almo_scf_env%eps_filter,error)
                               CALL opt_k_create_preconditioner_blk(almo_scf_env,&
                                       almo_scf_env%matrix_v_disc(ispin),&
                                       tmp4_o_vr,&
                                       t_curr,&
                                       ispin,&
                                       spin_factor,&
                                       error)
   
                            ENDIF
   
                            ! save the previous step
                            CALL cp_dbcsr_copy(prev_step,step,error=error)
   
                            ! compute the new step
                            CALL opt_k_apply_preconditioner_blk(almo_scf_env,&
                                    step,grad,ispin,error)
                            !CALL cp_dbcsr_hadamard_product(prec,grad,step,error=error)
                            CALL cp_dbcsr_scale(step,-1.0_dp,error=error)
   
                            ! check whether we need to reset conjugate directions
                            reset_conjugator=.FALSE.
                            ! first check if manual reset is active
                            IF (iteration.lt.MAX(almo_scf_env%opt_k_conj_iter_start,1).OR.&
                               MOD(iteration-almo_scf_env%opt_k_conj_iter_start,&
                               almo_scf_env%opt_k_conj_iter_freq).eq.0) THEN
   
                               reset_conjugator=.TRUE.
   
                            ELSE
   
                               ! check for the errors in the cg algorithm 
                               !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk,error=error)
                               !CALL cp_dbcsr_trace(grad,tmp_k_blk,numer,"T","N",error=error)
                               !CALL cp_dbcsr_trace(prev_grad,tmp_k_blk,denom,"T","N",error=error)
                               CALL cp_dbcsr_trace(grad,prev_minus_prec_grad,numer,"T","N",error=error)
                               CALL cp_dbcsr_trace(prev_grad,prev_minus_prec_grad,denom,"T","N",error=error)
                               conjugacy_error=numer/denom
   
                               IF (conjugacy_error.gt.MIN(0.5_dp,conjugacy_error_threshold)) THEN
                                  reset_conjugator=.TRUE.
                                  IF (unit_nr>0) THEN
                                     WRITE(unit_nr,*) "Lack of progress, conjugacy error is ", conjugacy_error
                                  ENDIF
                               ENDIF
   
                               ! check the gradient along the previous direction
                               IF ((iteration.ne.0).AND.(.NOT.reset_conjugator)) THEN
                                  CALL cp_dbcsr_trace(grad,prev_step,numer,"T","N",error=error)
                                  CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N",error=error)
                                  line_search_error=numer/denom
                                  IF (line_search_error.gt.line_search_error_threshold) THEN
                                     reset_conjugator=.TRUE.
                                     IF (unit_nr>0) THEN
                                        WRITE(unit_nr,*) "Bad line search, line search error is ", line_search_error
                                     ENDIF
                                  ENDIF
                               ENDIF
   
                            ENDIF
   
                            ! compute the conjugation coefficient - beta
                            IF (.NOT.reset_conjugator) THEN
   
                               SELECT CASE (almo_scf_env%opt_k_conjugator)
                               CASE (cg_hestenes_stiefel)
                                  CALL cp_dbcsr_copy(tmp_k_blk,grad,error=error)
                                  CALL cp_dbcsr_add(tmp_k_blk,prev_grad,1.0_dp,-1.0_dp,error=error)
                                  CALL cp_dbcsr_trace(tmp_k_blk,step,numer,"T","N",error=error)
                                  CALL cp_dbcsr_trace(tmp_k_blk,prev_step,denom,"T","N",error=error)
                                  beta=-1.0_dp*numer/denom
                               CASE (cg_fletcher_reeves)
                                  !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk,error=error)
                                  !CALL cp_dbcsr_trace(prev_grad,tmp_k_blk,denom,"T","N",error=error)
                                  !CALL cp_dbcsr_hadamard_product(prec,grad,tmp_k_blk,error=error)
                                  !CALL cp_dbcsr_trace(grad,tmp_k_blk,numer,"T","N",error=error)
                                  !beta=numer/denom
                                  CALL cp_dbcsr_trace(grad,step,numer,"T","N",error=error)
                                  CALL cp_dbcsr_trace(prev_grad,prev_minus_prec_grad,denom,"T","N",error=error)
                                  beta=numer/denom
                               CASE (cg_polak_ribiere)
                                  !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk,error=error)
                                  !CALL cp_dbcsr_trace(prev_grad,tmp_k_blk,denom,"T","N",error=error)
                                  !CALL cp_dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp,error=error)
                                  !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk,error=error)
                                  !CALL cp_dbcsr_trace(tmp_k_blk,grad,numer,"T","N",error=error)
                                  CALL cp_dbcsr_trace(prev_grad,prev_minus_prec_grad,denom,"T","N",error=error)
                                  CALL cp_dbcsr_copy(tmp_k_blk,grad,error=error)
                                  CALL cp_dbcsr_add(tmp_k_blk,prev_grad,1.0_dp,-1.0_dp,error=error)
                                  CALL cp_dbcsr_trace(tmp_k_blk,step,numer,"T","N",error=error)
                                  beta=numer/denom
                               CASE (cg_fletcher)
                                  !CALL cp_dbcsr_hadamard_product(prec,grad,tmp_k_blk,error=error)
                                  !CALL cp_dbcsr_trace(grad,tmp_k_blk,numer,"T","N",error=error)
                                  !CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N",error=error)
                                  !beta=-1.0_dp*numer/denom
                                  CALL cp_dbcsr_trace(grad,step,numer,"T","N",error=error)
                                  CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N",error=error)
                                  beta=numer/denom
                               CASE (cg_liu_storey)
                                  CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N",error=error)
                                  !CALL cp_dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp,error=error)
                                  !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk,error=error)
                                  !CALL cp_dbcsr_trace(tmp_k_blk,grad,numer,"T","N",error=error)
                                  CALL cp_dbcsr_copy(tmp_k_blk,grad,error=error)
                                  CALL cp_dbcsr_add(tmp_k_blk,prev_grad,1.0_dp,-1.0_dp,error=error)
                                  CALL cp_dbcsr_trace(tmp_k_blk,step,numer,"T","N",error=error)
                                  beta=numer/denom
                               CASE (cg_dai_yuan)
                                  !CALL cp_dbcsr_hadamard_product(prec,grad,tmp_k_blk,error=error)
                                  !CALL cp_dbcsr_trace(grad,tmp_k_blk,numer,"T","N",error=error)
                                  !CALL cp_dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp,error=error)
                                  !CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N",error=error)
                                  !beta=numer/denom
                                  CALL cp_dbcsr_trace(grad,step,numer,"T","N",error=error)
                                  CALL cp_dbcsr_copy(tmp_k_blk,grad,error=error)
                                  CALL cp_dbcsr_add(tmp_k_blk,prev_grad,1.0_dp,-1.0_dp,error=error)
                                  CALL cp_dbcsr_trace(tmp_k_blk,prev_step,denom,"T","N",error=error)
                                  beta=-1.0_dp*numer/denom
                               CASE (cg_hager_zhang)
                                  !CALL cp_dbcsr_add(prev_grad,grad,-1.0_dp,1.0_dp,error=error)
                                  !CALL cp_dbcsr_trace(prev_grad,prev_step,denom,"T","N",error=error)
                                  !CALL cp_dbcsr_hadamard_product(prec,prev_grad,tmp_k_blk,error=error)
                                  !CALL cp_dbcsr_trace(tmp_k_blk,prev_grad,numer,"T","N",error=error)
                                  !kappa=2.0_dp*numer/denom
                                  !CALL cp_dbcsr_trace(tmp_k_blk,grad,numer,"T","N",error=error)
                                  !tau=numer/denom
                                  !CALL cp_dbcsr_trace(prev_step,grad,numer,"T","N",error=error)
                                  !beta=tau-kappa*numer/denom
                                  CALL cp_dbcsr_copy(tmp_k_blk,grad,error=error)
                                  CALL cp_dbcsr_add(tmp_k_blk,prev_grad,1.0_dp,-1.0_dp,error=error)
                                  CALL cp_dbcsr_trace(tmp_k_blk,prev_step,denom,"T","N",error=error)
                                  CALL cp_dbcsr_trace(tmp_k_blk,prev_minus_prec_grad,numer,"T","N",error=error)
                                  kappa=-2.0_dp*numer/denom
                                  CALL cp_dbcsr_trace(tmp_k_blk,step,numer,"T","N",error=error)
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
   
                            IF (md_in_k_space) THEN
                               reset_conjugator=.TRUE.
                            ENDIF
   
                            IF (reset_conjugator) THEN 
   
                               beta=0.0_dp
                               !reset_step_size=.TRUE.
   
                               IF (unit_nr>0) THEN
                                  WRITE(unit_nr,*) "(Re)-setting conjugator to zero"
                               ENDIF
   
                            ENDIF
                      
                            ! save the preconditioned gradient
                            CALL cp_dbcsr_copy(prev_minus_prec_grad,step,error=error)
   
                            ! conjugate the step direction
                            CALL cp_dbcsr_add(step,prev_step,1.0_dp,beta,error=error)
   
                            CALL timestop(handle7)

                         ! update the step direction
                         ELSE ! step update
                            conjugacy_error=0.0_dp
                         ENDIF
   
                         ! compute the gradient with respect to the step size in the curr direction
                         IF (line_search) THEN
                            CALL cp_dbcsr_trace(grad,step,gfun1,"T","N",error=error)
                            line_search_error=gfun1/gfun0
                         ELSE
                            CALL cp_dbcsr_trace(grad,step,gfun0,"T","N",error=error)
                         ENDIF
   
                         ! make a step - update k 
                         IF (line_search) THEN
                           
                            ! check if the trial step provides enough numerical accuracy
                            safety_multiplier=1.0E+1_dp ! must be more than one
                            num_threshold=MAX(EPSILON(1.0_dp),&
                               safety_multiplier*(almo_scf_env%eps_filter**2)*almo_scf_env%ndomains)
                            IF (ABS(fun1-fun0-gfun0*step_size).lt.num_threshold) THEN
                               IF (unit_nr>0) THEN
                                  WRITE(unit_nr,'(T3,A,1X,E17.7)') &
                                     "Numerical accuracy is too low to observe non-linear behavior",&
                                     ABS(fun1-fun0-gfun0*step_size)
                                  WRITE(unit_nr,'(T3,A,1X,E17.7,A,1X,E12.3)') "Error computing ",&
                                     ABS(gfun0),&
                                     " is smaller than the threshold",num_threshold
                               ENDIF
                               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                            ENDIF
                            IF (ABS(gfun0).lt.num_threshold) THEN
                               IF (unit_nr>0) THEN
                                  WRITE(unit_nr,'(T3,A,1X,E17.7,A,1X,E12.3)') "Linear gradient",&
                                     ABS(gfun0),&
                                     " is smaller than the threshold",num_threshold
                               ENDIF
                               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                            ENDIF
                            
                            use_quadratic_approximation=.TRUE.
                            use_cubic_approximation=.FALSE.
   
                            ! find the minimum assuming quadratic form
                            ! use f0, f1, g0
                            step_size_quadratic_approx=-(gfun0*step_size*step_size)/(2.0_dp*(fun1-fun0-gfun0*step_size))
                            ! use f0, f1, g1
                            step_size_quadratic_approx2=-(fun1-fun0-step_size*gfun1/2.0_dp)/(gfun1-(fun1-fun0)/step_size)
                            
                            IF ((step_size_quadratic_approx.lt.0.0_dp).AND.&
                                (step_size_quadratic_approx2.lt.0.0_dp)) THEN
                               IF (unit_nr>0) THEN
                                  WRITE(unit_nr,'(T3,A,1X,E17.7,1X,E17.7,1X,A)') &
                                     "Quadratic approximation gives negative steps",&
                                     step_size_quadratic_approx,step_size_quadratic_approx2,&
                                     "trying cubic..."
                               ENDIF
                               use_cubic_approximation=.TRUE.
                               use_quadratic_approximation=.FALSE.
                            ELSE
                               IF (step_size_quadratic_approx.lt.0.0_dp) THEN
                                  step_size_quadratic_approx=step_size_quadratic_approx2
                               ENDIF
                               IF (step_size_quadratic_approx2.lt.0.0_dp) THEN
                                  step_size_quadratic_approx2=step_size_quadratic_approx
                               ENDIF
                            ENDIF
   
                            ! check accuracy of the quadratic approximation
                            IF (use_quadratic_approximation) THEN
                               quadratic_approx_error=ABS(step_size_quadratic_approx-&
                                  step_size_quadratic_approx2)/step_size_quadratic_approx
                               IF (quadratic_approx_error.gt.quadratic_approx_error_threshold) THEN
                                  IF (unit_nr>0) THEN
                                     WRITE(unit_nr,'(T3,A,1X,E17.7,1X,E17.7,1X,A)') "Quadratic approximation is poor",&
                                        step_size_quadratic_approx,step_size_quadratic_approx2,&
                                        "Try cubic approximation"
                                  ENDIF
                                  use_cubic_approximation=.TRUE.
                                  use_quadratic_approximation=.FALSE.
                               ENDIF
                            ENDIF
   
                            ! check if numerics is fine enough to capture the cubic form
                            IF (use_cubic_approximation) THEN
                               
                               ! if quadratic approximation is not accurate enough
                               ! try to find the minimum assuming cubic form
                               ! aa*x**3 + bb*x**2 + cc*x + dd = f(x)
                               bb = (-step_size*gfun1+3.0_dp*(fun1-fun0)-2.0_dp*step_size*gfun0)/(step_size*step_size)
                               aa = (gfun1-2.0_dp*step_size*bb-gfun0)/(3.0_dp*step_size*step_size)
                               
                               IF (ABS(gfun1-2.0_dp*step_size*bb-gfun0).lt.num_threshold) THEN
                                  IF (unit_nr>0) THEN
                                     WRITE(unit_nr,'(T3,A,1X,E17.7)') &
                                        "Numerical accuracy is too low to observe cubic behavior",&
                                        ABS(gfun1-2.0_dp*step_size*bb-gfun0)
                                  ENDIF
                                  use_cubic_approximation=.FALSE.
                                  use_quadratic_approximation=.TRUE.
                               ENDIF
                               IF (ABS(gfun1).lt.num_threshold) THEN
                                  IF (unit_nr>0) THEN
                                     WRITE(unit_nr,'(T3,A,1X,E17.7,A,1X,E12.3)') "Linear gradient",&
                                        ABS(gfun1),&
                                        " is smaller than the threshold",num_threshold
                                  ENDIF
                                  use_cubic_approximation=.FALSE.
                                  use_quadratic_approximation=.TRUE.
                               ENDIF
                            ENDIF
   
                            ! find the step assuming cubic approximation
                            IF (use_cubic_approximation) THEN
                               ! to obtain the minimum of the cubic function solve the quadratic equation
                               ! 0.0*x**3 + 3.0*aa*x**2 + 2.0*bb*x + cc = 0
                               CALL analytic_line_search(0.0_dp,3.0_dp*aa,2.0_dp*bb,gfun0,minima,nmins)
                               IF (nmins.lt.1) THEN
                                  IF (unit_nr>0) THEN
                                     WRITE(unit_nr,'(T3,A)')&
                                        "Cubic approximation gives zero soultions! Use quadratic approximation"
                                  ENDIF
                                  use_quadratic_approximation=.TRUE.
                                  use_cubic_approximation=.TRUE.
                               ELSE
                                  step_size=minima(1)
                                  IF (nmins.gt.1) THEN
                                     IF (unit_nr>0) THEN
                                        WRITE(unit_nr,'(T3,A)')&
                                           "More than one solution found! Use quadratic approximation"
                                     ENDIF
                                     use_quadratic_approximation=.TRUE.
                                     use_cubic_approximation=.TRUE.
                                  ENDIF
                               ENDIF
                            ENDIF
   
                            IF (use_quadratic_approximation) THEN ! use quadratic approximation
                               IF (unit_nr>0) THEN
                                  WRITE(unit_nr,'(T3,A)') "Use quadratic approximation"
                               ENDIF
                               step_size=(step_size_quadratic_approx+step_size_quadratic_approx2)*0.5_dp
                            ENDIF
   
                            ! one more check on the step size
                            IF (step_size.lt.0.0_dp) THEN
                               CPErrorMessage(cp_failure_level,routineP,"Negative step proposed",error)
                               CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
                            ENDIF
   
                            CALL cp_dbcsr_copy(almo_scf_env%matrix_k_blk(ispin),&
                                    matrix_k_central,error=error)
                            CALL cp_dbcsr_add(almo_scf_env%matrix_k_blk(ispin),&
                                    step,1.0_dp,step_size,error=error)
                            CALL cp_dbcsr_copy(matrix_k_central,&
                                    almo_scf_env%matrix_k_blk(ispin),error=error)
                            line_search=.FALSE.
   
                         ELSE
   
                            IF (md_in_k_space) THEN
   
                               ! update velocities v(i) = v(i-1) + 0.5*dT*(a(i-1) + a(i))
                               IF (iteration.ne.0) THEN
                                  CALL cp_dbcsr_add(velocity,&
                                          step,1.0_dp,0.5_dp*time_step,error=error)
                                  CALL cp_dbcsr_add(velocity,&
                                          prev_step,1.0_dp,0.5_dp*time_step,error=error)
                               ENDIF
                               kin_energy=cp_dbcsr_frobenius_norm(velocity)
                               kin_energy=0.5_dp*kin_energy*kin_energy
   
                               ! update positions k(i) = k(i-1) + dT*v(i-1) + 0.5*dT*dT*a(i-1)
                               CALL cp_dbcsr_add(almo_scf_env%matrix_k_blk(ispin),&
                                       velocity,1.0_dp,time_step,error=error)
                               CALL cp_dbcsr_add(almo_scf_env%matrix_k_blk(ispin),&
                                       step,1.0_dp,0.5_dp*time_step*time_step,error=error)
   
                            ELSE
   
                               IF (reset_step_size) THEN
                                  step_size=almo_scf_env%opt_k_trial_step_size
                                  reset_step_size=.FALSE.
                               ELSE
                                  step_size=step_size*almo_scf_env%opt_k_trial_step_size_multiplier 
                               ENDIF
                               CALL cp_dbcsr_copy(almo_scf_env%matrix_k_blk(ispin),&
                                       matrix_k_central,error=error)
                               CALL cp_dbcsr_add(almo_scf_env%matrix_k_blk(ispin),&
                                       step,1.0_dp,step_size,error=error)
                               line_search=.TRUE.
                            ENDIF
   
                         ENDIF
   
                      ENDIF ! .NOT.prepare_to_exit
   
                      ! print the status of the optimization
                      t2a = m_walltime()
                      IF (unit_nr>0) THEN
                         IF (md_in_k_space) THEN
                               WRITE(unit_nr,'(T6,A,1X,I5,1X,E12.3,E16.7,F15.9,F15.9,F15.9,E12.3,F15.9,F15.9,F8.3)') &
                                  "K iter CG",iteration,time_step,time_step*iteration,&
                                  energy_correction(ispin),obj_function,delta_obj_function,grad_norm,&
                                  kin_energy,kin_energy+obj_function,beta
                         ELSE
                            IF (line_search.OR.prepare_to_exit) THEN
                               WRITE(unit_nr,'(T6,A,1X,I3,1X,E12.3,F16.10,F16.10,E12.3,E12.3,E12.3,F8.3,F8.3,F10.3)') &
                                     "K iter CG",iteration,step_size,&
                                     energy_correction(ispin),delta_obj_function,grad_norm,&
                                     gfun0,line_search_error,beta,conjugacy_error,t2a-t1a
                                     !(flop1+flop2)/(1.0E6_dp*(t2-t1))
                            ELSE
                               WRITE(unit_nr,'(T6,A,1X,I3,1X,E12.3,F16.10,F16.10,E12.3,E12.3,E12.3,F8.3,F8.3,F10.3)') &
                                     "K iter LS",iteration,step_size,&
                                     energy_correction(ispin),delta_obj_function,grad_norm,&
                                     gfun1,line_search_error,beta,conjugacy_error,t2a-t1a
                                     !(flop1+flop2)/(1.0E6_dp*(t2-t1))
                            ENDIF
                         ENDIF
                         CALL m_flush(unit_nr)
                      ENDIF
                      t1a = m_walltime()
   
                   ELSE ! opt_k_max_iter .eq. 0
                      prepare_to_exit=.TRUE.
                   ENDIF ! opt_k_max_iter .ne. 0
   
                   IF (.NOT.line_search) iteration = iteration + 1
   
                   IF (prepare_to_exit) EXIT
   
                ENDDO ! end iterations on K
             
                IF (converged.OR.(outer_opt_k_iteration.ge.outer_opt_k_max_iter)) THEN
                   outer_opt_k_prepare_to_exit=.TRUE.
                ENDIF
             
                IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
                           
                   IF (unit_nr>0) THEN
                      WRITE(unit_nr,*) "Updating ALMO virtuals"
                   ENDIF
             
                   CALL timeset('k_opt_v0_update',handle8)
             
                   ! update retained ALMO virtuals to restart the cg iterations
                   CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           almo_scf_env%matrix_v_disc_blk(ispin),&
                           almo_scf_env%matrix_k_blk(ispin),&
                           0.0_dp,vr_fixed,&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_add(vr_fixed,almo_scf_env%matrix_v_blk(ispin),&
                           +1.0_dp,+1.0_dp,error=error)
                   
                   ! update discarded ALMO virtuals to restart the cg iterations
                   CALL cp_dbcsr_multiply("N","T",1.0_dp,&
                           almo_scf_env%matrix_v_blk(ispin),&
                           almo_scf_env%matrix_k_blk(ispin),&
                           0.0_dp,vd_fixed,&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_add(vd_fixed,almo_scf_env%matrix_v_disc_blk(ispin),&
                           -1.0_dp,+1.0_dp,error=error)
                   
                   ! orthogonalize new orbitals on fragments
                   CALL get_overlap(bra=vr_fixed,&
                           ket=vr_fixed,&
                           overlap=k_vr_index_down,&
                           metric=almo_scf_env%matrix_s_blk(1),&
                           retain_overlap_sparsity=.FALSE.,&
                           eps_filter=almo_scf_env%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_init(vr_index_sqrt_inv,error=error)
                   CALL cp_dbcsr_create(vr_index_sqrt_inv,template=k_vr_index_down,&
                                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                   CALL cp_dbcsr_init(vr_index_sqrt,error=error)
                   CALL cp_dbcsr_create(vr_index_sqrt,template=k_vr_index_down,&
                                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                   CALL matrix_sqrt_Newton_Schulz(vr_index_sqrt,&
                           vr_index_sqrt_inv,&
                           k_vr_index_down,&
                           threshold=almo_scf_env%eps_filter,&
                           order=almo_scf_env%order_lanczos,&
                           eps_lanczos=almo_scf_env%eps_lanczos,&
                           max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
                           error=error)
                   IF (safe_mode) THEN
                      CALL cp_dbcsr_init(matrix_tmp1,error=error)
                      CALL cp_dbcsr_create(matrix_tmp1,template=k_vr_index_down,&
                                           matrix_type=dbcsr_type_no_symmetry,error=error) 
                      CALL cp_dbcsr_init(matrix_tmp2,error=error)
                      CALL cp_dbcsr_create(matrix_tmp2,template=k_vr_index_down,&
                                           matrix_type=dbcsr_type_no_symmetry,error=error) 
                   
                      CALL cp_dbcsr_multiply("N","N",1.0_dp,vr_index_sqrt_inv,&
                                             k_vr_index_down,&
                                             0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter,&
                                             error=error)
                      CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_tmp1,&
                                             vr_index_sqrt_inv,&
                                             0.0_dp,matrix_tmp2,filter_eps=almo_scf_env%eps_filter,&
                                             error=error)
                   
                      frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
                      CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp,error=error)
                      frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
                      IF (unit_nr>0) THEN
                         WRITE(unit_nr,*) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)",&
                            frob_matrix/frob_matrix_base
                      ENDIF
                   
                      CALL cp_dbcsr_release(matrix_tmp1,error=error) 
                      CALL cp_dbcsr_release(matrix_tmp2,error=error) 
                   ENDIF
                   CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           vr_fixed,& 
                           vr_index_sqrt_inv,&
                           0.0_dp,almo_scf_env%matrix_v_blk(ispin),&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
             
                   CALL get_overlap(bra=vd_fixed,&
                           ket=vd_fixed,&
                           overlap=k_vd_index_down,&
                           metric=almo_scf_env%matrix_s_blk(1),&
                           retain_overlap_sparsity=.FALSE.,&
                           eps_filter=almo_scf_env%eps_filter,&
                           error=error)
                   CALL cp_dbcsr_init(vd_index_sqrt_inv,error=error)
                   CALL cp_dbcsr_create(vd_index_sqrt_inv,template=k_vd_index_down,&
                                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                   CALL cp_dbcsr_init(vd_index_sqrt,error=error)
                   CALL cp_dbcsr_create(vd_index_sqrt,template=k_vd_index_down,&
                                        matrix_type=dbcsr_type_no_symmetry,error=error) 
                   CALL matrix_sqrt_Newton_Schulz(vd_index_sqrt,&
                           vd_index_sqrt_inv,&
                           k_vd_index_down,&
                           threshold=almo_scf_env%eps_filter,&
                           order=almo_scf_env%order_lanczos,&
                           eps_lanczos=almo_scf_env%eps_lanczos,&
                           max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
                           error=error)
                   IF (safe_mode) THEN
                      CALL cp_dbcsr_init(matrix_tmp1,error=error)
                      CALL cp_dbcsr_create(matrix_tmp1,template=k_vd_index_down,&
                                           matrix_type=dbcsr_type_no_symmetry,error=error) 
                      CALL cp_dbcsr_init(matrix_tmp2,error=error)
                      CALL cp_dbcsr_create(matrix_tmp2,template=k_vd_index_down,&
                                           matrix_type=dbcsr_type_no_symmetry,error=error) 
                   
                      CALL cp_dbcsr_multiply("N","N",1.0_dp,vd_index_sqrt_inv,&
                                             k_vd_index_down,&
                                             0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter,&
                                             error=error)
                      CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_tmp1,&
                                             vd_index_sqrt_inv,&
                                             0.0_dp,matrix_tmp2,filter_eps=almo_scf_env%eps_filter,&
                                             error=error)
                   
                      frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
                      CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp,error=error)
                      frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
                      IF (unit_nr>0) THEN
                         WRITE(unit_nr,*) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)",&
                            frob_matrix/frob_matrix_base
                      ENDIF
                   
                      CALL cp_dbcsr_release(matrix_tmp1,error=error) 
                      CALL cp_dbcsr_release(matrix_tmp2,error=error) 
                   ENDIF
                   CALL cp_dbcsr_multiply("N","N",1.0_dp,&
                           vd_fixed,& 
                           vd_index_sqrt_inv,&
                           0.0_dp,almo_scf_env%matrix_v_disc_blk(ispin),&
                           filter_eps=almo_scf_env%eps_filter,&
                           error=error)
                  
                   CALL cp_dbcsr_release(vr_index_sqrt_inv,error=error)
                   CALL cp_dbcsr_release(vr_index_sqrt,error=error)
                   CALL cp_dbcsr_release(vd_index_sqrt_inv,error=error)
                   CALL cp_dbcsr_release(vd_index_sqrt,error=error)
             
                   CALL timestop(handle8)

                ENDIF ! ne.virt_full
             
                ! RZK-warning released outside the outer loop
                CALL cp_dbcsr_release(sigma_vv_sqrt,error=error)
                CALL cp_dbcsr_release(sigma_vv_sqrt_inv,error=error)
                IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
                   CALL cp_dbcsr_release(k_vr_index_down,error=error)
                   CALL cp_dbcsr_release(k_vd_index_down,error=error)
                   !CALL cp_dbcsr_release(k_vd_index_up,error=error)
                   CALL cp_dbcsr_release(matrix_k_central,error=error)
                   CALL cp_dbcsr_release(vr_fixed,error=error)
                   CALL cp_dbcsr_release(vd_fixed,error=error)
                   CALL cp_dbcsr_release(grad,error=error)
                   CALL cp_dbcsr_release(prec,error=error)
                   CALL cp_dbcsr_release(prev_grad,error=error)
                   CALL cp_dbcsr_release(tmp3_vd_vr,error=error)
                   CALL cp_dbcsr_release(tmp1_n_vr,error=error)
                   CALL cp_dbcsr_release(tmp_k_blk,error=error)
                   CALL cp_dbcsr_release(t_curr,error=error)
                   CALL cp_dbcsr_release(sigma_oo_curr,error=error)
                   CALL cp_dbcsr_release(sigma_oo_curr_inv,error=error)
                   CALL cp_dbcsr_release(step,error=error)
                   CALL cp_dbcsr_release(tmp2_n_o,error=error)
                   CALL cp_dbcsr_release(tmp4_o_vr,error=error)
                   CALL cp_dbcsr_release(prev_step,error=error)
                   CALL cp_dbcsr_release(prev_minus_prec_grad,error=error)
                   IF (md_in_k_space) THEN
                      CALL cp_dbcsr_release(velocity,error=error)
                   ENDIF
                
                ENDIF

                outer_opt_k_iteration=outer_opt_k_iteration+1
                IF (outer_opt_k_prepare_to_exit) EXIT

             ENDDO ! outer loop for k

          ENDDO ! ispin

          ! RZK-warning update mo orbitals

       ELSE ! virtual orbitals might not be available use projected AOs

          ! compute sqrt(S) and inv(sqrt(S))
          ! RZK-warning - remove this sqrt(S) and inv(sqrt(S))
          ! ideally ALMO scf should use sigma and sigma_inv in
          ! the tensor_up_down representation
          IF (.NOT.almo_scf_env%s_sqrt_done) THEN
   
             IF (unit_nr>0) THEN
                WRITE(unit_nr,*) "sqrt and inv(sqrt) of AO overlap matrix"
             ENDIF
             CALL cp_dbcsr_init(almo_scf_env%matrix_s_sqrt(1),error=error)
             CALL cp_dbcsr_init(almo_scf_env%matrix_s_sqrt_inv(1),error=error)
             CALL cp_dbcsr_create(almo_scf_env%matrix_s_sqrt(1),&
                                  template=almo_scf_env%matrix_s(1),&
                                  matrix_type=dbcsr_type_no_symmetry,error=error) 
             CALL cp_dbcsr_create(almo_scf_env%matrix_s_sqrt_inv(1),&
                                  template=almo_scf_env%matrix_s(1),&
                                  matrix_type=dbcsr_type_no_symmetry,error=error) 
      
             CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_s_sqrt(1),&
                                            almo_scf_env%matrix_s_sqrt_inv(1),&
                                            almo_scf_env%matrix_s(1),&
                                            threshold=almo_scf_env%eps_filter,&
                                            order=almo_scf_env%order_lanczos,&
                                            eps_lanczos=almo_scf_env%eps_lanczos,&
                                            max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
                                            error=error)
      
             IF (safe_mode) THEN
                CALL cp_dbcsr_init(matrix_tmp1,error=error)
                CALL cp_dbcsr_create(matrix_tmp1,template=almo_scf_env%matrix_s(1),&
                                     matrix_type=dbcsr_type_no_symmetry,error=error) 
                CALL cp_dbcsr_init(matrix_tmp2,error=error)
                CALL cp_dbcsr_create(matrix_tmp2,template=almo_scf_env%matrix_s(1),&
                                     matrix_type=dbcsr_type_no_symmetry,error=error) 
      
                CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_s_sqrt_inv(1),&
                                       almo_scf_env%matrix_s(1),&
                                       0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter,&
                                       error=error)
                CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_tmp1,almo_scf_env%matrix_s_sqrt_inv(1),&
                                       0.0_dp,matrix_tmp2,filter_eps=almo_scf_env%eps_filter,error=error)
      
                frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
                CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp,error=error)
                frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
                IF (unit_nr>0) THEN
                   WRITE(unit_nr,*) "Error for (inv(sqrt(S))*S*inv(sqrt(S))-I)",frob_matrix/frob_matrix_base
                ENDIF
      
                CALL cp_dbcsr_release(matrix_tmp1,error=error) 
                CALL cp_dbcsr_release(matrix_tmp2,error=error) 
             ENDIF
   
             almo_scf_env%s_sqrt_done=.TRUE.
             
          ENDIF
   
          DO ispin=1,nspin
          
             CALL ct_step_env_init(ct_step_env,error)
             CALL ct_step_env_set(ct_step_env,&
                para_env=almo_scf_env%para_env,&
                blacs_env=almo_scf_env%blacs_env,&
                use_occ_orbs=.TRUE.,&
                use_virt_orbs=almo_scf_env%deloc_cayley_use_virt_orbs,&
                occ_orbs_orthogonal=.FALSE.,&
                virt_orbs_orthogonal=almo_scf_env%orthogonal_basis,&
                tensor_type=almo_scf_env%deloc_cayley_tensor_type,&
                neglect_quadratic_term=almo_scf_env%deloc_cayley_linear,&
                calculate_energy_corr=.TRUE.,&
                update_p=.TRUE.,&
                update_q=.FALSE.,&
                pp_preconditioner_full=almo_scf_env%deloc_cayley_occ_precond,&
                qq_preconditioner_full=almo_scf_env%deloc_cayley_vir_precond,&
                eps_convergence=almo_scf_env%deloc_cayley_eps_convergence,&
                eps_filter=almo_scf_env%eps_filter,&
                !nspins=almo_scf_env%nspins,&
                q_index_up=almo_scf_env%matrix_s_sqrt_inv(1),&
                q_index_down=almo_scf_env%matrix_s_sqrt(1),&
                p_index_up=almo_scf_env%matrix_sigma_sqrt_inv(ispin),&
                p_index_down=almo_scf_env%matrix_sigma_sqrt(ispin),&
                matrix_ks=almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                matrix_p=almo_scf_env%matrix_p(ispin),&
                matrix_qp_template=almo_scf_env%matrix_t(ispin),&
                matrix_pq_template=almo_scf_env%matrix_t_tr(ispin),&
                matrix_t=almo_scf_env%matrix_t(ispin),&
                conjugator=almo_scf_env%deloc_cayley_conjugator,&
                max_iter=almo_scf_env%deloc_cayley_max_iter,&
                error=error)
             
             ! perform calculations
             CALL ct_step_execute(ct_step_env,error)
   
             ! for now we do not need the new set of orbitals
             ! just get the energy correction
             CALL ct_step_env_get(ct_step_env,&
                energy_correction=energy_correction(ispin),&
                !copy_da_energy_matrix=matrix_eda(ispin),&
                !copy_da_charge_matrix=matrix_cta(ispin),&
                error=error)
   
             CALL ct_step_env_clean(ct_step_env,error)

          ENDDO

          energy_correction(1)=energy_correction(1)*spin_factor
       
       ENDIF
       
       ! print the energy correction and exit
       DO ispin=1,nspin
   
          IF (unit_nr>0) THEN
             WRITE(unit_nr,*)
             WRITE(unit_nr,'(T2,A,I6,F20.9)') "ECORR",ispin,&
                     energy_correction(ispin)
             WRITE(unit_nr,*)
          ENDIF
          energy_correction_final=energy_correction_final+energy_correction(ispin)
          
          !!! print out the results of decomposition analysis
          !!IF (unit_nr>0) THEN
          !!   WRITE(unit_nr,*)
          !!   WRITE(unit_nr,'(T2,A)') "ENERGY DECOMPOSITION"
          !!ENDIF
          !!CALL cp_dbcsr_print_block_sum(eda_matrix(ispin),error=error)
          !!IF (unit_nr>0) THEN
          !!   WRITE(unit_nr,*)
          !!   WRITE(unit_nr,'(T2,A)') "CHARGE DECOMPOSITION"
          !!ENDIF
          !!CALL cp_dbcsr_print_block_sum(cta_matrix(ispin),error=error)

          ! obtain density matrix from updated MOs
          ! RZK-later sigma and sigma_inv are lost here
          CALL almo_scf_t_to_p(t=almo_scf_env%matrix_t(ispin),&
                               p=almo_scf_env%matrix_p(ispin),&
                               eps_filter=almo_scf_env%eps_filter,&
                               orthog_orbs=.FALSE.,&
                               s=almo_scf_env%matrix_s(1),&
                               sigma=almo_scf_env%matrix_sigma(ispin),&
                               sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
                               error=error)
          
          IF (almo_scf_env%nspins==1) &
             CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin),&
                    spin_factor,error=error)
       
       ENDDO
         

    CASE(dm_ls_step)
   
       ! compute the inverse of S
       IF (.NOT.almo_scf_env%s_inv_done) THEN
          IF (unit_nr>0) THEN
             WRITE(unit_nr,*) "Inverting AO overlap matrix"
          ENDIF
          CALL cp_dbcsr_init(almo_scf_env%matrix_s_inv(1),error=error)
          CALL cp_dbcsr_create(almo_scf_env%matrix_s_inv(1),&
                               template=almo_scf_env%matrix_s(1),&
                               matrix_type=dbcsr_type_no_symmetry,&
                               error=error) 
          IF (.NOT.almo_scf_env%s_sqrt_done) THEN
             CALL invert_Hotelling(almo_scf_env%matrix_s_inv(1),&
                                   almo_scf_env%matrix_s(1),&
                                   threshold=almo_scf_env%eps_filter,&
                                   error=error)
          ELSE
             CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_s_sqrt_inv(1),&
                                    almo_scf_env%matrix_s_sqrt_inv(1), &
                                    0.0_dp,almo_scf_env%matrix_s_inv(1),&
                                    filter_eps=almo_scf_env%eps_filter,error=error)
          ENDIF
   
          IF (safe_mode) THEN
             CALL cp_dbcsr_init(matrix_tmp1,error=error)
             CALL cp_dbcsr_create(matrix_tmp1,template=almo_scf_env%matrix_s(1),&
                                  matrix_type=dbcsr_type_no_symmetry,error=error) 
             CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_inv(1),&
                                    almo_scf_env%matrix_s(1),&
                                    0.0_dp, matrix_tmp1,&
                                    filter_eps=almo_scf_env%eps_filter,error=error)
             frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp1)
             CALL cp_dbcsr_add_on_diag(matrix_tmp1,-1.0_dp,error=error)
             frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp1)
             IF (unit_nr>0) THEN
                WRITE(unit_nr,*) "Error for (inv(S)*S-I)",&
                      frob_matrix/frob_matrix_base
             ENDIF
             CALL cp_dbcsr_release(matrix_tmp1,error=error)
          ENDIF
   
          almo_scf_env%s_inv_done=.TRUE.
   
       ENDIF
       
       DO ispin=1,nspin
          ! RZK-warning the preconditioner is very important
          !       IF (.FALSE.) THEN
          !           CALL apply_matrix_preconditioner(almo_scf_env%matrix_ks(ispin),&
          !                   "forward",almo_scf_env%matrix_s_blk_sqrt(1),&
          !                   almo_scf_env%matrix_s_blk_sqrt_inv(1),error)
          !       ENDIF
          !CALL cp_dbcsr_filter(almo_scf_env%matrix_ks(ispin),&
          !         almo_scf_env%eps_filter,error=error)
       ENDDO
   
       ALLOCATE(matrix_p_almo_scf_converged(nspin))
       DO ispin=1,nspin
          CALL cp_dbcsr_init(matrix_p_almo_scf_converged(ispin),error=error)
          CALL cp_dbcsr_create(matrix_p_almo_scf_converged(ispin),&
                 template=almo_scf_env%matrix_p(ispin),error=error)
          CALL cp_dbcsr_copy(matrix_p_almo_scf_converged(ispin),&
                 almo_scf_env%matrix_p(ispin),error=error)
       ENDDO
       
       ! update the density matrix
       DO ispin=1,nspin
   
          nelectron_spin_real(1)=almo_scf_env%nelectrons_spin(ispin)
          IF (almo_scf_env%nspins==1) &
                 nelectron_spin_real(1)=nelectron_spin_real(1)/2
   
          local_mu(1)=SUM(almo_scf_env%mu_of_domain(:,ispin))/almo_scf_env%ndomains
          fake(1)=123523
   
          ! RZK UPDATE! the update algorithm is removed because
          ! RZK UPDATE! it requires updating core LS_SCF routines
          ! RZK UPDATE! (the code exists in the CVS version)
          IF (unit_nr>0) WRITE(unit_nr,'(T2,A)') "CVS only: density_matrix_sign has not been updated in SVN"
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          ! RZK UPDATE!CALL density_matrix_sign(almo_scf_env%matrix_p(ispin),&
          ! RZK UPDATE!                     local_mu,&
          ! RZK UPDATE!                     almo_scf_env%fixed_mu,&
          ! RZK UPDATE!                     almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
          ! RZK UPDATE!                     almo_scf_env%matrix_s(1), &
          ! RZK UPDATE!                     almo_scf_env%matrix_s_inv(1), &
          ! RZK UPDATE!                     nelectron_spin_real,&
          ! RZK UPDATE!                     almo_scf_env%eps_filter,&
          ! RZK UPDATE!                     fake,&
          ! RZK UPDATE!                     error)
          almo_scf_env%mu=local_mu(1)
   
          !IF (almo_scf_env%has_s_preconditioner) THEN
          !    CALL apply_matrix_preconditioner(&
          !             almo_scf_env%matrix_p_blk(ispin),&
          !             "forward",almo_scf_env%matrix_s_blk_sqrt(1),&
          !             almo_scf_env%matrix_s_blk_sqrt_inv(1),error)
          !ENDIF
          !CALL cp_dbcsr_filter(almo_scf_env%matrix_p(ispin),&
          !        almo_scf_env%eps_filter,error=error)
   
          IF (almo_scf_env%nspins==1) &
             CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin),&
                    spin_factor,error=error)
   
          !CALL cp_dbcsr_trace(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
          !  almo_scf_env%matrix_p(ispin),&
          !  energy_correction(ispin),error=error)
          !IF (unit_nr>0) THEN
          !   WRITE(unit_nr,*)
          !   WRITE(unit_nr,'(T2,A,I6,F20.9)') "EFAKE",ispin,&
          !           energy_correction(ispin)
          !   WRITE(unit_nr,*)
          !ENDIF
          CALL cp_dbcsr_add(matrix_p_almo_scf_converged(ispin),&
                           almo_scf_env%matrix_p(ispin),-1.0_dp,1.0_dp,&
                           error=error)
          CALL cp_dbcsr_trace(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                              matrix_p_almo_scf_converged(ispin),&
                              energy_correction(ispin),error=error)
         
          energy_correction_final=energy_correction_final+energy_correction(ispin)
   
          IF (unit_nr>0) THEN
             WRITE(unit_nr,*)
             WRITE(unit_nr,'(T2,A,I6,F20.9)') "ECORR",ispin,&
                     energy_correction(ispin)
             WRITE(unit_nr,*)
          ENDIF
   
       ENDDO
       
       DO ispin=1,nspin
          CALL cp_dbcsr_release(matrix_p_almo_scf_converged(ispin),error=error)
       ENDDO
       DEALLOCATE(matrix_p_almo_scf_converged)
   
    END SELECT ! algorithm selection

    t2 = m_walltime()

    IF (unit_nr>0) THEN
       WRITE(unit_nr,*)
       WRITE(unit_nr,'(T2,A,F18.9,F18.9,F18.9,F12.6)') "ETOT",&
               almo_scf_env%almo_scf_energy,&
               energy_correction_final,&
               almo_scf_env%almo_scf_energy+energy_correction_final,&
               t2-t1
       WRITE(unit_nr,*)
    ENDIF

    CALL timestop(handle)
  
  END SUBROUTINE harris_foulkes_correction

! *****************************************************************************
!> \brief Computes a diagonal preconditioner for the cg optimization of k matrix
!> \par History
!>       2011.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE opt_k_create_preconditioner(prec,vd_prop,f,x,oo_inv_x_tr,s,grad,&
       vd_blk,t,template_vd_vd_blk,template_vr_vr_blk,template_n_vr,&
       spin_factor,eps_filter,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: prec
    TYPE(cp_dbcsr_type), INTENT(IN) :: vd_prop, f, x, oo_inv_x_tr, s, grad, &
      vd_blk, t, template_vd_vd_blk, template_vr_vr_blk, template_n_vr
    REAL(KIND=dp), INTENT(IN)                :: spin_factor, eps_filter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'opt_k_create_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, p_nrows, q_nrows
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: p_diagonal, q_diagonal
    TYPE(cp_dbcsr_type)                      :: pp_diag, qq_diag, t1, t2, &
                                                tmp, tmp1_n_vr, tmp2_n_vr, &
                                                tmp_n_vd, tmp_vd_vd_blk, &
                                                tmp_vr_vr_blk

! init diag blocks outside
! init diag blocks otside
!INTEGER                                  :: iblock_row, iblock_col,&
!                                            nblkrows_tot, nblkcols_tot
!REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
!INTEGER                                  :: mynode, hold, row, col

    CALL timeset(routineN,handle)

    ! initialize a matrix to 1.0
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=prec,error=error)
    ! in order to use cp_dbcsr_set matrix blocks must exist
    CALL cp_dbcsr_copy(tmp,prec,error=error)
    CALL cp_dbcsr_set(tmp,1.0_dp,error=error)

    ! compute qq = (Vd^tr)*F*Vd
    CALL cp_dbcsr_init(tmp_n_vd,error=error)
    CALL cp_dbcsr_create(tmp_n_vd,template=vd_prop,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,f,vd_prop,&
             0.0_dp,tmp_n_vd,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_init(tmp_vd_vd_blk,error=error)
    CALL cp_dbcsr_create(tmp_vd_vd_blk,&
            template=template_vd_vd_blk,error=error)
    CALL cp_dbcsr_copy(tmp_vd_vd_blk,template_vd_vd_blk,error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,vd_prop,tmp_n_vd,&
             0.0_dp,tmp_vd_vd_blk,&
             retain_sparsity=.TRUE.,&
             filter_eps=eps_filter,error=error)
    ! copy diagonal elements of the result into rows of a matrix
    CALL cp_dbcsr_get_info(tmp_vd_vd_blk, nfullrows_total=q_nrows )
    ALLOCATE(q_diagonal(q_nrows))
    CALL cp_dbcsr_get_diag(tmp_vd_vd_blk,q_diagonal,error=error)
    CALL cp_dbcsr_init(qq_diag,error=error)
    CALL cp_dbcsr_create(qq_diag,&
            template=template_vd_vd_blk,error=error)
    CALL cp_dbcsr_add_on_diag(qq_diag,1.0_dp,error=error)
    CALL cp_dbcsr_set_diag(qq_diag,q_diagonal,error=error)
    CALL cp_dbcsr_init(t1,error=error)
    CALL cp_dbcsr_create(t1,template=prec,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,qq_diag,tmp,&
             0.0_dp,t1,filter_eps=eps_filter,error=error)

    ! compute pp = X*sigma_oo_inv*X^tr
    CALL cp_dbcsr_init(tmp_vr_vr_blk,error=error)
    CALL cp_dbcsr_create(tmp_vr_vr_blk,template=template_vr_vr_blk,error=error)
    CALL cp_dbcsr_copy(tmp_vr_vr_blk,template_vr_vr_blk,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,x,oo_inv_x_tr,&
             0.0_dp,tmp_vr_vr_blk,&
             retain_sparsity=.TRUE.,&
             filter_eps=eps_filter,error=error)
    ! copy diagonal elements of the result into cols of a matrix
    CALL cp_dbcsr_get_info(tmp_vr_vr_blk,nfullrows_total=p_nrows)
    ALLOCATE(p_diagonal(p_nrows))
    CALL cp_dbcsr_get_diag(tmp_vr_vr_blk,p_diagonal,error=error)
    CALL cp_dbcsr_init(pp_diag,error=error)
    CALL cp_dbcsr_create(pp_diag,template=template_vr_vr_blk,error=error)
    CALL cp_dbcsr_add_on_diag(pp_diag,1.0_dp,error=error)
    CALL cp_dbcsr_set_diag(pp_diag,p_diagonal,error=error)
    CALL cp_dbcsr_set(tmp,1.0_dp,error=error)
    CALL cp_dbcsr_init(t2,error=error)
    CALL cp_dbcsr_create(t2,template=prec,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,tmp,pp_diag,&
             0.0_dp,t2,filter_eps=eps_filter,error=error)

    CALL cp_dbcsr_hadamard_product(t1,t2,prec,error=error)

    ! compute qq = (Vd^tr)*S*Vd
    CALL cp_dbcsr_multiply("N","N",1.0_dp,s,vd_prop,&
             0.0_dp,tmp_n_vd,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,vd_prop,tmp_n_vd,&
             0.0_dp,tmp_vd_vd_blk,&
             retain_sparsity=.TRUE.,&
             filter_eps=eps_filter,error=error)
    ! copy diagonal elements of the result into rows of a matrix
    CALL cp_dbcsr_get_diag(tmp_vd_vd_blk,q_diagonal,error=error)
    CALL cp_dbcsr_add_on_diag(qq_diag,1.0_dp,error=error)
    CALL cp_dbcsr_set_diag(qq_diag,q_diagonal,error=error)
    CALL cp_dbcsr_set(tmp,1.0_dp,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,qq_diag,tmp,&
             0.0_dp,t1,filter_eps=eps_filter,error=error)
    
    ! compute pp = X*sig_oo_inv*(T^tr)*F*T*sig_oo_inv*(X^tr)
    CALL cp_dbcsr_init(tmp1_n_vr,error=error)
    CALL cp_dbcsr_create(tmp1_n_vr,template=template_n_vr,error=error)
    CALL cp_dbcsr_init(tmp2_n_vr,error=error)
    CALL cp_dbcsr_create(tmp2_n_vr,template=template_n_vr,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,t,oo_inv_x_tr,&
             0.0_dp,tmp1_n_vr,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,f,tmp1_n_vr,&
             0.0_dp,tmp2_n_vr,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,tmp1_n_vr,tmp2_n_vr,&
             0.0_dp,tmp_vr_vr_blk,&
             retain_sparsity=.TRUE.,&
             filter_eps=eps_filter,error=error)
    ! copy diagonal elements of the result into cols of a matrix
    CALL cp_dbcsr_get_diag(tmp_vr_vr_blk,p_diagonal,error=error)
    CALL cp_dbcsr_add_on_diag(pp_diag,1.0_dp,error=error)
    CALL cp_dbcsr_set_diag(pp_diag,p_diagonal,error=error)
    CALL cp_dbcsr_set(tmp,1.0_dp,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,tmp,pp_diag,&
             0.0_dp,t2,filter_eps=eps_filter,error=error)

    CALL cp_dbcsr_hadamard_product(t1,t2,tmp,error=error)
    CALL cp_dbcsr_add(prec,tmp,1.0_dp,-1.0_dp,error=error)
    CALL cp_dbcsr_scale(prec,2.0_dp*spin_factor,error=error)

    ! compute qp = X*sig_oo_inv*(T^tr)*S*Vd
    CALL cp_dbcsr_multiply("N","N",1.0_dp,s,vd_blk,&
             0.0_dp,tmp_n_vd,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,tmp_n_vd,tmp1_n_vr,&
             0.0_dp,tmp,retain_sparsity=.TRUE.,&
             filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_hadamard_product(grad,tmp,t1,error=error)
    ! gradient already contains 2.0*spin_factor
    CALL cp_dbcsr_scale(t1,-2.0_dp,error=error)

    CALL cp_dbcsr_add(prec,t1,1.0_dp,1.0_dp,error=error)
    
    CALL cp_dbcsr_function_of_elements(prec,dbcsr_func_inverse,error=error)
    CALL cp_dbcsr_filter(prec,eps_filter,error=error)

    DEALLOCATE(q_diagonal)
    DEALLOCATE(p_diagonal)
    CALL cp_dbcsr_release(tmp,error=error)
    CALL cp_dbcsr_release(qq_diag,error=error)
    CALL cp_dbcsr_release(t1,error=error)
    CALL cp_dbcsr_release(pp_diag,error=error)
    CALL cp_dbcsr_release(t2,error=error)
    CALL cp_dbcsr_release(tmp_n_vd,error=error)
    CALL cp_dbcsr_release(tmp_vd_vd_blk,error=error)
    CALL cp_dbcsr_release(tmp_vr_vr_blk,error=error)
    CALL cp_dbcsr_release(tmp1_n_vr,error=error)
    CALL cp_dbcsr_release(tmp2_n_vr,error=error)

    CALL timestop(handle)

  END SUBROUTINE opt_k_create_preconditioner

! *****************************************************************************
!> \brief Computes a block-diagonal preconditioner for the optimization of 
!>        k matrix
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE opt_k_create_preconditioner_blk(almo_scf_env,vd_prop,oo_inv_x_tr,&
    t_curr,ispin,spin_factor,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_dbcsr_type), INTENT(IN)          :: vd_prop, oo_inv_x_tr, t_curr
    INTEGER, INTENT(IN)                      :: ispin
    REAL(KIND=dp), INTENT(IN)                :: spin_factor
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'opt_k_create_preconditioner_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: eps_filter
    TYPE(cp_dbcsr_type) :: opt_k_e_dd, opt_k_e_rr, s_dd_sqrt, s_rr_sqrt, t1, &
      tmp, tmp1_n_vr, tmp2_n_vr, tmp_n_vd, tmp_vd_vd_blk, tmp_vr_vr_blk

! matrices that has been computed outside the routine already

    CALL timeset(routineN,handle)
    
    eps_filter=almo_scf_env%eps_filter

    ! compute S_qq = (Vd^tr)*S*Vd
    CALL cp_dbcsr_init(tmp_n_vd,error=error)
    CALL cp_dbcsr_create(tmp_n_vd,template=almo_scf_env%matrix_v_disc(ispin),&
            error=error)
    CALL cp_dbcsr_init(tmp_vd_vd_blk,error=error)
    CALL cp_dbcsr_create(tmp_vd_vd_blk,&
            template=almo_scf_env%matrix_vv_disc_blk(ispin),&
            matrix_type=dbcsr_type_no_symmetry,error=error) 
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            almo_scf_env%matrix_s(1),&
            vd_prop,&
            0.0_dp,tmp_n_vd,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_copy(tmp_vd_vd_blk,&
            almo_scf_env%matrix_vv_disc_blk(ispin),error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,vd_prop,tmp_n_vd,&
            0.0_dp,tmp_vd_vd_blk,&
            retain_sparsity=.TRUE.,&
            error=error)

    CALL cp_dbcsr_init(s_dd_sqrt,error=error)
    CALL cp_dbcsr_create(s_dd_sqrt,&
            template=almo_scf_env%matrix_vv_disc_blk(ispin),&
            matrix_type=dbcsr_type_no_symmetry,error=error) 
    CALL matrix_sqrt_Newton_Schulz(s_dd_sqrt,&
            almo_scf_env%opt_k_t_dd(ispin),&
            tmp_vd_vd_blk,&
            threshold=eps_filter,&
            order=almo_scf_env%order_lanczos,&
            eps_lanczos=almo_scf_env%eps_lanczos,&
            max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
            error=error)
    
    ! compute F_qq = (Vd^tr)*F*Vd
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
            vd_prop,&
            0.0_dp,tmp_n_vd,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_copy(tmp_vd_vd_blk,&
            almo_scf_env%matrix_vv_disc_blk(ispin),error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,vd_prop,tmp_n_vd,&
            0.0_dp,tmp_vd_vd_blk,&
            retain_sparsity=.TRUE.,&
            error=error)
    CALL cp_dbcsr_release(tmp_n_vd,error=error)

    ! bring to the blocked-orthogonalized basis
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            tmp_vd_vd_blk,&
            almo_scf_env%opt_k_t_dd(ispin),&
            0.0_dp,s_dd_sqrt,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            almo_scf_env%opt_k_t_dd(ispin),&
            s_dd_sqrt,&
            0.0_dp,tmp_vd_vd_blk,filter_eps=eps_filter,error=error)

    ! diagonalize the matrix
    CALL cp_dbcsr_init(opt_k_e_dd,error=error)
    CALL cp_dbcsr_create(opt_k_e_dd,&
            template=almo_scf_env%matrix_vv_disc_blk(ispin),&
            error=error)
    CALL cp_dbcsr_release(s_dd_sqrt,error=error)
    CALL cp_dbcsr_init(s_dd_sqrt,error=error)
    CALL cp_dbcsr_create(s_dd_sqrt,&
            template=almo_scf_env%matrix_vv_disc_blk(ispin),&
            matrix_type=dbcsr_type_no_symmetry,error=error) 
    CALL diagonalize_diagonal_blocks(tmp_vd_vd_blk,&
            s_dd_sqrt,&
            opt_k_e_dd,&
            error=error)

    ! obtain the transformation matrix in the discarded subspace
    ! T = S^{-1/2}.U
    CALL cp_dbcsr_copy(tmp_vd_vd_blk,&
            almo_scf_env%opt_k_t_dd(ispin),error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            tmp_vd_vd_blk,&
            s_dd_sqrt,&
            0.0_dp,almo_scf_env%opt_k_t_dd(ispin),&
            filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_release(s_dd_sqrt,error=error)
    CALL cp_dbcsr_release(tmp_vd_vd_blk,error=error)

    ! copy diagonal elements of the result into rows of a matrix
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,&
            template=almo_scf_env%matrix_k_blk_ones(ispin),&
            error=error)
    CALL cp_dbcsr_copy(tmp,&
            almo_scf_env%matrix_k_blk_ones(ispin),&
            error=error)
    CALL cp_dbcsr_init(t1,error=error)
    CALL cp_dbcsr_create(t1,&
            template=almo_scf_env%matrix_k_blk_ones(ispin),&
            error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            opt_k_e_dd,tmp,&
            0.0_dp,t1,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_release(opt_k_e_dd,error=error)
    
    ! compute S_pp = X*sigma_oo_inv*X^tr
    CALL cp_dbcsr_init(tmp_vr_vr_blk,error=error)
    CALL cp_dbcsr_create(tmp_vr_vr_blk,&
            template=almo_scf_env%matrix_sigma_vv_blk(ispin),&
            matrix_type=dbcsr_type_no_symmetry,error=error) 
    CALL cp_dbcsr_copy(tmp_vr_vr_blk,&
            almo_scf_env%matrix_sigma_vv_blk(ispin),error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            almo_scf_env%matrix_x(ispin),&
            oo_inv_x_tr,&
            0.0_dp,tmp_vr_vr_blk,&
            retain_sparsity=.TRUE.,&
            error=error)
      
    ! obtain the orthogonalization matrix
    CALL cp_dbcsr_init(s_rr_sqrt,error=error)
    CALL cp_dbcsr_create(s_rr_sqrt,&
            template=almo_scf_env%matrix_sigma_vv_blk(ispin),&
            matrix_type=dbcsr_type_no_symmetry,error=error) 
    CALL matrix_sqrt_Newton_Schulz(s_rr_sqrt,&
            almo_scf_env%opt_k_t_rr(ispin),&
            tmp_vr_vr_blk,&
            threshold=eps_filter,&
            order=almo_scf_env%order_lanczos,&
            eps_lanczos=almo_scf_env%eps_lanczos,&
            max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
            error=error)

    ! compute F_pp = X*sig_oo_inv*(T^tr)*F*T*sig_oo_inv*(X^tr)
    CALL cp_dbcsr_init(tmp1_n_vr,error=error)
    CALL cp_dbcsr_create(tmp1_n_vr,&
            template=almo_scf_env%matrix_v(ispin),error=error)
    CALL cp_dbcsr_init(tmp2_n_vr,error=error)
    CALL cp_dbcsr_create(tmp2_n_vr,&
            template=almo_scf_env%matrix_v(ispin),error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,t_curr,oo_inv_x_tr,&
            0.0_dp,tmp1_n_vr,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
            tmp1_n_vr,&
            0.0_dp,tmp2_n_vr,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,tmp1_n_vr,tmp2_n_vr,&
            0.0_dp,tmp_vr_vr_blk,&
            retain_sparsity=.TRUE.,&
            error=error)
    CALL cp_dbcsr_release(tmp1_n_vr,error=error)
    CALL cp_dbcsr_release(tmp2_n_vr,error=error)
    
    ! bring to the blocked-orthogonalized basis
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            tmp_vr_vr_blk,&
            almo_scf_env%opt_k_t_rr(ispin),&
            0.0_dp,s_rr_sqrt,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            almo_scf_env%opt_k_t_rr(ispin),&
            s_rr_sqrt,&
            0.0_dp,tmp_vr_vr_blk,filter_eps=eps_filter,error=error)

    ! diagonalize the matrix
    CALL cp_dbcsr_init(opt_k_e_rr,error=error)
    CALL cp_dbcsr_create(opt_k_e_rr,&
            template=almo_scf_env%matrix_sigma_vv_blk(ispin),&
            error=error)
    CALL cp_dbcsr_release(s_rr_sqrt,error=error)
    CALL cp_dbcsr_init(s_rr_sqrt,error=error)
    CALL cp_dbcsr_create(s_rr_sqrt,&
            template=almo_scf_env%matrix_sigma_vv_blk(ispin),&
            matrix_type=dbcsr_type_no_symmetry,error=error) 
    CALL diagonalize_diagonal_blocks(tmp_vr_vr_blk,&
            s_rr_sqrt,&
            opt_k_e_rr,&
            error=error)

    ! obtain the transformation matrix in the retained subspace
    ! T = S^{-1/2}.U
    CALL cp_dbcsr_copy(tmp_vr_vr_blk,&
            almo_scf_env%opt_k_t_rr(ispin),error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            tmp_vr_vr_blk,&
            s_rr_sqrt,&
            0.0_dp,almo_scf_env%opt_k_t_rr(ispin),&
            filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_release(s_rr_sqrt,error=error)
    CALL cp_dbcsr_release(tmp_vr_vr_blk,error=error)

    ! copy diagonal elements of the result into cols of a matrix
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            tmp,opt_k_e_rr,&
            0.0_dp,almo_scf_env%opt_k_denom(ispin),&
            filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_release(opt_k_e_rr,error=error)
    CALL cp_dbcsr_release(tmp,error=error)

    ! form the denominator matrix
    CALL cp_dbcsr_add(almo_scf_env%opt_k_denom(ispin),t1,&
            -1.0_dp,1.0_dp,error=error)
    CALL cp_dbcsr_release(t1,error=error)
    CALL cp_dbcsr_scale(almo_scf_env%opt_k_denom(ispin),&
            2.0_dp*spin_factor,error=error)

    CALL cp_dbcsr_function_of_elements(almo_scf_env%opt_k_denom(ispin),&
            dbcsr_func_inverse,error=error)
    CALL cp_dbcsr_filter(almo_scf_env%opt_k_denom(ispin),&
            eps_filter,error=error)

    CALL timestop(handle)

  END SUBROUTINE opt_k_create_preconditioner_blk

! *****************************************************************************
!> \brief Applies a block-diagonal preconditioner for the optimization of 
!>        k matrix (preconditioner matrices must be calculated and stored 
!>        beforehand)
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE opt_k_apply_preconditioner_blk(almo_scf_env,step,grad,ispin,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_dbcsr_type), INTENT(OUT)         :: step
    TYPE(cp_dbcsr_type), INTENT(IN)          :: grad
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'opt_k_apply_preconditioner_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: eps_filter
    TYPE(cp_dbcsr_type)                      :: tmp_k

    CALL timeset(routineN,handle)
    
    eps_filter=almo_scf_env%eps_filter

    CALL cp_dbcsr_init(tmp_k,error=error)
    CALL cp_dbcsr_create(tmp_k,template=almo_scf_env%matrix_k_blk(ispin),&
            error=error)

    ! transform gradient to the correct "diagonal" basis
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            grad,almo_scf_env%opt_k_t_rr(ispin),&
            0.0_dp,tmp_k,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,&
            almo_scf_env%opt_k_t_dd(ispin),tmp_k,&
            0.0_dp,step,filter_eps=eps_filter,error=error)

    ! apply diagonal preconditioner
    CALL cp_dbcsr_hadamard_product(step,&
            almo_scf_env%opt_k_denom(ispin),tmp_k,error=error)

    ! back-transform the result to the initial basis
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            almo_scf_env%opt_k_t_dd(ispin),tmp_k,&
            0.0_dp,step,filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_multiply("N","T",1.0_dp,&
            step,almo_scf_env%opt_k_t_rr(ispin),&
            0.0_dp,tmp_k,filter_eps=eps_filter,error=error)

    CALL cp_dbcsr_copy(step,tmp_k,error=error)

    CALL cp_dbcsr_release(tmp_k,error=error)

    CALL timestop(handle)

  END SUBROUTINE opt_k_apply_preconditioner_blk

! *****************************************************************************
!> \brief create various matrices 
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_env_create_matrices(almo_scf_env,matrix_s0,error)
    
    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s0
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_env_create_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins

    CALL timeset(routineN,handle)
   
    nspins = almo_scf_env%nspins

    ! AO overlap matrix and its various functions
    CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s(1),&
         matrix_qs=matrix_s0,&
         almo_scf_env=almo_scf_env,&
         name_new="S",&
         size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
         symmetry_new=dbcsr_type_symmetric,&
         spin_key=0,&
         init_domains=.FALSE.,&
         error=error)
    CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s_blk(1),&
         matrix_qs=matrix_s0,&
         almo_scf_env=almo_scf_env,&
         name_new="S_BLK",&
         size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
         symmetry_new=dbcsr_type_symmetric,&
         spin_key=0,&
         init_domains=.TRUE.,&
         error=error)
    IF (almo_scf_env%almo_update_algorithm.eq.almo_scf_diag) THEN
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s_blk_sqrt_inv(1),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="S_BLK_SQRT_INV",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=0,&
            init_domains=.TRUE.,&
            error=error)
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s_blk_sqrt(1),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="S_BLK_SQRT",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=0,&
            init_domains=.TRUE.,&
            error=error)
    ELSE IF (almo_scf_env%almo_update_algorithm.eq.almo_scf_dm_sign) THEN
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s_blk_inv(1),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="S_BLK_INV",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=0,&
            init_domains=.TRUE.,&
            error=error)
    ENDIF

    ! MO coeff matrices and their derivatives
    ALLOCATE(almo_scf_env%matrix_t_blk(nspins))
    ALLOCATE(almo_scf_env%quench_t_blk(nspins))
    ALLOCATE(almo_scf_env%matrix_err_blk(nspins))
    ALLOCATE(almo_scf_env%matrix_err_xx(nspins))
    ALLOCATE(almo_scf_env%matrix_sigma(nspins))
    ALLOCATE(almo_scf_env%matrix_sigma_inv(nspins))
    ALLOCATE(almo_scf_env%matrix_sigma_sqrt(nspins))
    ALLOCATE(almo_scf_env%matrix_sigma_sqrt_inv(nspins))
    ALLOCATE(almo_scf_env%matrix_sigma_blk(nspins))
    ALLOCATE(almo_scf_env%matrix_t(nspins))
    ALLOCATE(almo_scf_env%matrix_t_tr(nspins))
    DO ispin=1,nspins
       ! create the blocked quencher
       CALL matrix_almo_create(matrix_new=almo_scf_env%quench_t_blk(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="Q_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       ! create ALMO coefficient matrix
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_t_blk(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="T_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       ! create the error matrix
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_err_blk(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="ERR_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       ! create the error matrix for the quenched ALMOs
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_err_xx(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="ERR_XX",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       ! create a matrix with dimensions of a transposed mo coefficient matrix
       ! it might be necessary to perform the correction step using cayley 
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_t_tr(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="T_TR",&
            size_keys=(/almo_mat_dim_occ,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       ! create mo overlap matrix
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="SIG",&
            size_keys=(/almo_mat_dim_occ,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       ! create blocked mo overlap matrix
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma_blk(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="SIG_BLK",&
            size_keys=(/almo_mat_dim_occ,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       ! create inverse mo overlap matrix
       CALL matrix_almo_create(&
            matrix_new=almo_scf_env%matrix_sigma_inv(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="SIGINV",&
            size_keys=(/almo_mat_dim_occ,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       ! create various templates that will be necessary later
       CALL matrix_almo_create(&
            matrix_new=almo_scf_env%matrix_t(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="T",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_occ/),&
            symmetry_new=dbcsr_type_no_symmetry,&
            spin_key=ispin,&
            init_domains=.FALSE.,&
            error=error)
       CALL cp_dbcsr_init(almo_scf_env%matrix_sigma_sqrt(ispin),error=error)
       CALL cp_dbcsr_init(almo_scf_env%matrix_sigma_sqrt_inv(ispin),error=error)
       CALL cp_dbcsr_create(almo_scf_env%matrix_sigma_sqrt(ispin),&
                            template=almo_scf_env%matrix_sigma(ispin),&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
       CALL cp_dbcsr_create(almo_scf_env%matrix_sigma_sqrt_inv(ispin),&
                            template=almo_scf_env%matrix_sigma(ispin),&
                            matrix_type=dbcsr_type_no_symmetry, error=error)
    ENDDO
   
    ! create virtual orbitals if necessary
    IF (almo_scf_env%need_virtuals) THEN
       ALLOCATE(almo_scf_env%matrix_v_blk(nspins))
       ALLOCATE(almo_scf_env%matrix_v_full_blk(nspins))
       ALLOCATE(almo_scf_env%matrix_v(nspins))
       ALLOCATE(almo_scf_env%matrix_vo(nspins))
       ALLOCATE(almo_scf_env%matrix_x(nspins))
       ALLOCATE(almo_scf_env%matrix_ov(nspins))
       ALLOCATE(almo_scf_env%matrix_ov_full(nspins))
       ALLOCATE(almo_scf_env%matrix_sigma_vv(nspins))
       ALLOCATE(almo_scf_env%matrix_sigma_vv_blk(nspins))
       ALLOCATE(almo_scf_env%matrix_sigma_vv_sqrt(nspins))
       ALLOCATE(almo_scf_env%matrix_sigma_vv_sqrt_inv(nspins))
       ALLOCATE(almo_scf_env%matrix_vv_full_blk(nspins))
                   
       IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
          ALLOCATE(almo_scf_env%matrix_k_blk(nspins))
          ALLOCATE(almo_scf_env%matrix_k_blk_ones(nspins))
          ALLOCATE(almo_scf_env%matrix_k_tr(nspins))
          ALLOCATE(almo_scf_env%matrix_v_disc(nspins))
          ALLOCATE(almo_scf_env%matrix_v_disc_blk(nspins))
          ALLOCATE(almo_scf_env%matrix_ov_disc(nspins))
          ALLOCATE(almo_scf_env%matrix_vv_disc_blk(nspins))
          ALLOCATE(almo_scf_env%matrix_vv_disc(nspins))
          ALLOCATE(almo_scf_env%opt_k_t_dd(nspins))
          ALLOCATE(almo_scf_env%opt_k_t_rr(nspins))
          ALLOCATE(almo_scf_env%opt_k_denom(nspins))
       ENDIF

       DO ispin=1,nspins
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v_full_blk(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="V_FULL_BLK",&
               size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt_full/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v_blk(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="V_BLK",&
               size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="V",&
               size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_ov_full(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="OV_FULL",&
               size_keys=(/almo_mat_dim_occ,almo_mat_dim_virt_full/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_ov(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="OV",&
               size_keys=(/almo_mat_dim_occ,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_vo(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="VO",&
               size_keys=(/almo_mat_dim_virt,almo_mat_dim_occ/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_x(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="VO",&
               size_keys=(/almo_mat_dim_virt,almo_mat_dim_occ/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma_vv(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="SIG_VV",&
               size_keys=(/almo_mat_dim_virt,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_symmetric,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_vv_full_blk(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="VV_FULL_BLK",&
               size_keys=(/almo_mat_dim_virt_full,almo_mat_dim_virt_full/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.TRUE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma_vv_blk(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="SIG_VV_BLK",&
               size_keys=(/almo_mat_dim_virt,almo_mat_dim_virt/),&
               symmetry_new=dbcsr_type_symmetric,&
               spin_key=ispin,&
               init_domains=.TRUE.,&
               error=error)
          CALL cp_dbcsr_init(almo_scf_env%matrix_sigma_vv_sqrt(ispin),error=error)
          CALL cp_dbcsr_init(almo_scf_env%matrix_sigma_vv_sqrt_inv(ispin),error=error)
          CALL cp_dbcsr_create(almo_scf_env%matrix_sigma_vv_sqrt(ispin),&
                  template=almo_scf_env%matrix_sigma_vv(ispin),&
                  matrix_type=dbcsr_type_no_symmetry, error=error)
          CALL cp_dbcsr_create(almo_scf_env%matrix_sigma_vv_sqrt_inv(ispin),&
                  template=almo_scf_env%matrix_sigma_vv(ispin),&
                  matrix_type=dbcsr_type_no_symmetry, error=error)
          
          IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
             CALL matrix_almo_create(matrix_new=almo_scf_env%opt_k_t_rr(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="OPT_K_U_RR",&
                  size_keys=(/almo_mat_dim_virt,almo_mat_dim_virt/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_vv_disc(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="VV_DISC",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_symmetric,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%opt_k_t_dd(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="OPT_K_U_DD",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_vv_disc_blk(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="VV_DISC_BLK",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_symmetric,&
                  spin_key=ispin,&
                  init_domains=.TRUE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_k_blk(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="K_BLK",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.TRUE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_k_blk_ones(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="K_BLK_1",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.TRUE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%opt_k_denom(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="OPT_K_DENOM",&
                  size_keys=(/almo_mat_dim_virt_disc,almo_mat_dim_virt/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_k_tr(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="K_TR",&
                  size_keys=(/almo_mat_dim_virt,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v_disc_blk(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="V_DISC_BLK",&
                  size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v_disc(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="V_DISC",&
                  size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
             CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_ov_disc(ispin),&
                  matrix_qs=matrix_s0,&
                  almo_scf_env=almo_scf_env,&
                  name_new="OV_DISC",&
                  size_keys=(/almo_mat_dim_occ,almo_mat_dim_virt_disc/),&
                  symmetry_new=dbcsr_type_no_symmetry,&
                  spin_key=ispin,&
                  init_domains=.FALSE.,&
                  error=error)
        
          ENDIF ! end need_discarded_virtuals

       ENDDO ! spin
    ENDIF
    
    ! create matrices of orbital energies if necessary
    IF (almo_scf_env%need_orbital_energies) THEN
       ALLOCATE(almo_scf_env%matrix_eoo(nspins))
       ALLOCATE(almo_scf_env%matrix_evv_full(nspins))
       DO ispin=1,nspins
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_eoo(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="E_OCC",&
               size_keys=(/almo_mat_dim_occ,almo_mat_dim_occ/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
          CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_evv_full(ispin),&
               matrix_qs=matrix_s0,&
               almo_scf_env=almo_scf_env,&
               name_new="E_VIRT",&
               size_keys=(/almo_mat_dim_virt_full,almo_mat_dim_virt_full/),&
               symmetry_new=dbcsr_type_no_symmetry,&
               spin_key=ispin,&
               init_domains=.FALSE.,&
               error=error)
       ENDDO
    ENDIF

    ! Density and KS matrices 
    ALLOCATE(almo_scf_env%matrix_p(nspins))
    ALLOCATE(almo_scf_env%matrix_p_blk(nspins))
    ALLOCATE(almo_scf_env%matrix_ks(nspins))
    ALLOCATE(almo_scf_env%matrix_ks_blk(nspins))
    IF (almo_scf_env%need_previous_ks)&
       ALLOCATE(almo_scf_env%matrix_ks_almo_scf_converged(nspins))
    DO ispin=1,nspins
       CALL cp_dbcsr_init(almo_scf_env%matrix_p(ispin),error=error)
       ! RZK-warning copy with symmery but remember that this might cause problems 
       CALL cp_dbcsr_create(almo_scf_env%matrix_p(ispin),&
                            template=almo_scf_env%matrix_s(1),&
                            matrix_type=dbcsr_type_symmetric, error=error)
       CALL cp_dbcsr_init(almo_scf_env%matrix_ks(ispin),error=error)
       CALL cp_dbcsr_create(almo_scf_env%matrix_ks(ispin),&
                            template=almo_scf_env%matrix_s(1),&
                            matrix_type=dbcsr_type_symmetric, error=error)
       IF (almo_scf_env%need_previous_ks) THEN
          CALL cp_dbcsr_init(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                          error=error)
          CALL cp_dbcsr_create(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                            template=almo_scf_env%matrix_s(1),&
                            matrix_type=dbcsr_type_symmetric, error=error)
       ENDIF
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_p_blk(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="P_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
       CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_ks_blk(ispin),&
            matrix_qs=matrix_s0,&
            almo_scf_env=almo_scf_env,&
            name_new="KS_BLK",&
            size_keys=(/almo_mat_dim_aobasis,almo_mat_dim_aobasis/),&
            symmetry_new=dbcsr_type_symmetric,&
            spin_key=ispin,&
            init_domains=.TRUE.,&
            error=error)
    ENDDO
    
    CALL timestop(handle)
  
  END SUBROUTINE almo_scf_env_create_matrices 

! *****************************************************************************
!> \brief clean up procedures for almo scf
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_clean_up(qs_env,almo_scf_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(almo_scf_env_type)                  :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_clean_up', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! release matrices
    CALL cp_dbcsr_release(almo_scf_env%matrix_s(1),error=error)
    CALL cp_dbcsr_release(almo_scf_env%matrix_s_blk(1),error=error)
    IF (almo_scf_env%almo_update_algorithm.eq.almo_scf_diag) THEN
       CALL cp_dbcsr_release(almo_scf_env%matrix_s_blk_sqrt_inv(1),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_s_blk_sqrt(1),error=error)
    ELSE IF (almo_scf_env%almo_update_algorithm.eq.almo_scf_dm_sign) THEN
       CALL cp_dbcsr_release(almo_scf_env%matrix_s_blk_inv(1),error=error)
    ENDIF
    DO ispin=1,almo_scf_env%nspins
       CALL cp_dbcsr_release(almo_scf_env%quench_t(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%quench_t_blk(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_t_blk(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_err_blk(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_err_xx(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_t_tr(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_sigma(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_sigma_blk(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_sigma_inv(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_t(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_sigma_sqrt(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_sigma_sqrt_inv(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_p(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_ks(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_p_blk(ispin),error=error)
       CALL cp_dbcsr_release(almo_scf_env%matrix_ks_blk(ispin),error=error)
       IF (almo_scf_env%need_previous_ks) THEN
          CALL cp_dbcsr_release(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                  error=error)
       ENDIF
       IF (almo_scf_env%need_virtuals) THEN
          CALL cp_dbcsr_release(almo_scf_env%matrix_v_blk(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_v_full_blk(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_v(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_vo(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_x(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_ov(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_ov_full(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_sigma_vv(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_sigma_vv_blk(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_sigma_vv_sqrt(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_sigma_vv_sqrt_inv(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_vv_full_blk(ispin),error=error)
          IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
             CALL cp_dbcsr_release(almo_scf_env%matrix_k_tr(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%matrix_k_blk(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%matrix_k_blk_ones(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%matrix_v_disc(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%matrix_v_disc_blk(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%matrix_ov_disc(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%matrix_vv_disc_blk(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%matrix_vv_disc(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%opt_k_t_dd(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%opt_k_t_rr(ispin),error=error)
             CALL cp_dbcsr_release(almo_scf_env%opt_k_denom(ispin),error=error)
          ENDIF
       ENDIF
       IF (almo_scf_env%need_orbital_energies) THEN
          CALL cp_dbcsr_release(almo_scf_env%matrix_eoo(ispin),error=error)
          CALL cp_dbcsr_release(almo_scf_env%matrix_evv_full(ispin),error=error)
       ENDIF
    ENDDO

    ! deallocate matrices
    DEALLOCATE(almo_scf_env%matrix_p)
    DEALLOCATE(almo_scf_env%matrix_p_blk)
    DEALLOCATE(almo_scf_env%matrix_ks)
    DEALLOCATE(almo_scf_env%matrix_ks_blk)
    DEALLOCATE(almo_scf_env%matrix_t_blk)
    DEALLOCATE(almo_scf_env%matrix_err_blk)
    DEALLOCATE(almo_scf_env%matrix_err_xx)
    DEALLOCATE(almo_scf_env%matrix_t)
    DEALLOCATE(almo_scf_env%matrix_t_tr)
    DEALLOCATE(almo_scf_env%matrix_sigma)
    DEALLOCATE(almo_scf_env%matrix_sigma_blk)
    DEALLOCATE(almo_scf_env%matrix_sigma_sqrt)
    DEALLOCATE(almo_scf_env%matrix_sigma_sqrt_inv)
    DEALLOCATE(almo_scf_env%matrix_sigma_inv)
    DEALLOCATE(almo_scf_env%quench_t)
    DEALLOCATE(almo_scf_env%quench_t_blk)
    IF (almo_scf_env%need_virtuals) THEN
       DEALLOCATE(almo_scf_env%matrix_v_blk)
       DEALLOCATE(almo_scf_env%matrix_v_full_blk)
       DEALLOCATE(almo_scf_env%matrix_v)
       DEALLOCATE(almo_scf_env%matrix_vo)
       DEALLOCATE(almo_scf_env%matrix_x)
       DEALLOCATE(almo_scf_env%matrix_ov)
       DEALLOCATE(almo_scf_env%matrix_ov_full)
       DEALLOCATE(almo_scf_env%matrix_sigma_vv)
       DEALLOCATE(almo_scf_env%matrix_sigma_vv_blk)
       DEALLOCATE(almo_scf_env%matrix_sigma_vv_sqrt)
       DEALLOCATE(almo_scf_env%matrix_sigma_vv_sqrt_inv)
       DEALLOCATE(almo_scf_env%matrix_vv_full_blk)
       IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
          DEALLOCATE(almo_scf_env%matrix_k_tr)
          DEALLOCATE(almo_scf_env%matrix_k_blk)
          DEALLOCATE(almo_scf_env%matrix_v_disc)
          DEALLOCATE(almo_scf_env%matrix_v_disc_blk)
          DEALLOCATE(almo_scf_env%matrix_ov_disc)
          DEALLOCATE(almo_scf_env%matrix_vv_disc_blk)
          DEALLOCATE(almo_scf_env%matrix_vv_disc)
          DEALLOCATE(almo_scf_env%matrix_k_blk_ones)
          DEALLOCATE(almo_scf_env%opt_k_t_dd)
          DEALLOCATE(almo_scf_env%opt_k_t_rr)
          DEALLOCATE(almo_scf_env%opt_k_denom)
       ENDIF
    ENDIF
    IF (almo_scf_env%need_previous_ks) THEN
       DEALLOCATE(almo_scf_env%matrix_ks_almo_scf_converged)
    ENDIF
    IF (almo_scf_env%need_orbital_energies) THEN
       DEALLOCATE(almo_scf_env%matrix_eoo)
       DEALLOCATE(almo_scf_env%matrix_evv_full)
    ENDIF

    ! clean up other variables
    DO ispin=1,almo_scf_env%nspins
       CALL release_submatrices(&
              almo_scf_env%domain_preconditioner(:,ispin),&
              error=error)
       CALL release_submatrices(almo_scf_env%domain_s_inv(:,ispin),&
              error=error)
       CALL release_submatrices(almo_scf_env%domain_s_sqrt_inv(:,ispin),&
              error=error)
       CALL release_submatrices(almo_scf_env%domain_s_sqrt(:,ispin),&
              error=error)
       CALL release_submatrices(almo_scf_env%domain_ks_xx(:,ispin),&
              error=error)
       CALL release_submatrices(almo_scf_env%domain_t(:,ispin),&
              error=error)
       CALL release_submatrices(almo_scf_env%domain_err(:,ispin),&
              error=error)
       CALL release_submatrices(almo_scf_env%domain_r_down_up(:,ispin),&
              error=error)
    ENDDO
    DEALLOCATE(almo_scf_env%domain_preconditioner)
    DEALLOCATE(almo_scf_env%domain_s_inv)
    DEALLOCATE(almo_scf_env%domain_s_sqrt_inv)
    DEALLOCATE(almo_scf_env%domain_s_sqrt)
    DEALLOCATE(almo_scf_env%domain_ks_xx)
    DEALLOCATE(almo_scf_env%domain_t)
    DEALLOCATE(almo_scf_env%domain_err)
    DEALLOCATE(almo_scf_env%domain_r_down_up)
    DO ispin=1,almo_scf_env%nspins
       DEALLOCATE(almo_scf_env%domain_map(ispin)%pairs)
       DEALLOCATE(almo_scf_env%domain_map(ispin)%index1)
    ENDDO
    DEALLOCATE(almo_scf_env%domain_map)
    DEALLOCATE(almo_scf_env%domain_index_of_ao)
    DEALLOCATE(almo_scf_env%domain_index_of_atom)
    DEALLOCATE(almo_scf_env%first_atom_of_domain)
    DEALLOCATE(almo_scf_env%nbasis_of_domain)
    DEALLOCATE(almo_scf_env%nocc_of_domain)
    DEALLOCATE(almo_scf_env%nvirt_full_of_domain)
    DEALLOCATE(almo_scf_env%nvirt_of_domain)
    DEALLOCATE(almo_scf_env%nvirt_disc_of_domain)
    DEALLOCATE(almo_scf_env%mu_of_domain)
    DEALLOCATE(almo_scf_env%cpu_of_domain)
    DEALLOCATE(almo_scf_env%charge_of_domain)
    
    DEALLOCATE(almo_scf_env%domain_index_of_ao_block)
    DEALLOCATE(almo_scf_env%domain_index_of_mo_block)

    CALL cp_para_env_release(almo_scf_env%para_env,error)
    CALL cp_blacs_env_release(almo_scf_env%blacs_env,error)

    CALL timestop(handle)

  END SUBROUTINE almo_scf_clean_up

!! *****************************************************************************
!!> \brief Reduce the number of virtual orbitals by rotating them within 
!!>        a domain. The rotation is such that minimizes the frobenius norm of
!!>        the Fov domain-blocks of the discarded virtuals
!!> \par History
!!>       2011.08 created [Rustam Z Khaliullin]
!!> \author Rustam Z Khaliullin
!! *****************************************************************************
!  SUBROUTINE truncate_subspace_v_blk(qs_env,almo_scf_env,error)
!
!    TYPE(qs_environment_type), POINTER       :: qs_env
!    TYPE(almo_scf_env_type)                  :: almo_scf_env
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'truncate_subspace_v_blk', &
!      routineP = moduleN//':'//routineN
!
!    INTEGER                                  :: handle, ispin, iblock_row, &
!                                                iblock_col, iblock_row_size, &
!                                                iblock_col_size, retained_v, &
!                                                iteration, line_search_step, &
!                                                unit_nr, line_search_step_last 
!    REAL(KIND=dp)                            :: t1, obj_function, grad_norm,&
!                                                c0, b0, a0, obj_function_new,&
!                                                t2, alpha, ff1, ff2, step1,&
!                                                step2,&
!                                                frob_matrix_base,&
!                                                frob_matrix
!    LOGICAL                                  :: safe_mode, converged, &
!                                                prepare_to_exit, failure
!    TYPE(cp_logger_type), POINTER            :: logger
!    TYPE(cp_dbcsr_type)                      :: Fon, Fov, Fov_filtered, &
!                                                temp1_oo, temp2_oo, Fov_original, &
!                                                temp0_ov, U_blk_tot, U_blk, &
!                                                grad_blk, step_blk, matrix_filter, &
!                                                v_full_new,v_full_tmp,&
!                                                matrix_sigma_vv_full,&
!                                                matrix_sigma_vv_full_sqrt,&
!                                                matrix_sigma_vv_full_sqrt_inv,&
!                                                matrix_tmp1,&
!                                                matrix_tmp2
!
!    REAL(kind=dp), DIMENSION(:, :), POINTER  :: data_p, p_new_block
!    TYPE(cp_dbcsr_iterator)                  :: iter
!
!
!REAL(kind=dp), DIMENSION(:), ALLOCATABLE     :: eigenvalues, WORK
!REAL(kind=dp), DIMENSION(:,:), ALLOCATABLE   :: data_copy, left_vectors, right_vectors
!INTEGER                                      :: LWORK, INFO
!TYPE(cp_dbcsr_type)                          :: temp_u_v_full_blk
!       
!    CALL timeset(routineN,handle)
!
!    safe_mode=.TRUE.
!
!    ! get a useful output_unit
!    logger => cp_error_get_logger(error)
!    IF (logger%para_env%mepos==logger%para_env%source) THEN
!       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
!    ELSE
!       unit_nr=-1
!    ENDIF
!
!    DO ispin=1,almo_scf_env%nspins
!       
!       t1 = m_walltime()
!       
!       !!!!!!!!!!!!!!!!!
!       ! 0. Orthogonalize virtuals
!       !    Unfortunately, we have to do it in the FULL V subspace :(
!       
!       CALL cp_dbcsr_init(v_full_new,error=error)
!       CALL cp_dbcsr_create(v_full_new,&
!               template=almo_scf_env%matrix_v_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!       
!       ! project the occupied subspace out
!       CALL almo_scf_p_out_from_v(almo_scf_env%matrix_v_full_blk(ispin),&
!              v_full_new,almo_scf_env%matrix_ov_full(ispin),&
!              ispin,almo_scf_env,error)
!         
!       ! init overlap and its functions
!       CALL cp_dbcsr_init(matrix_sigma_vv_full,error=error)
!       CALL cp_dbcsr_init(matrix_sigma_vv_full_sqrt,error=error)
!       CALL cp_dbcsr_init(matrix_sigma_vv_full_sqrt_inv,error=error)
!       CALL cp_dbcsr_create(matrix_sigma_vv_full,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,error=error) 
!       CALL cp_dbcsr_create(matrix_sigma_vv_full_sqrt,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,error=error) 
!       CALL cp_dbcsr_create(matrix_sigma_vv_full_sqrt_inv,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,error=error)
!
!       ! construct VV overlap
!       CALL almo_scf_mo_to_sigma(v_full_new,&
!               matrix_sigma_vv_full,&
!               almo_scf_env%matrix_s(1),&
!               almo_scf_env%eps_filter,error)
!
!       IF (unit_nr>0) THEN
!          WRITE(unit_nr,*) "sqrt and inv(sqrt) of the FULL virtual MO overlap"
!       ENDIF
!
!       ! construct orthogonalization matrices       
!       CALL matrix_sqrt_Newton_Schulz(matrix_sigma_vv_full_sqrt,&
!                                      matrix_sigma_vv_full_sqrt_inv,&
!                                      matrix_sigma_vv_full,&
!                                      threshold=almo_scf_env%eps_filter,&
!                                      order=almo_scf_env%order_lanczos,&
!                                      eps_lanczos=almo_scf_env%eps_lanczos,&
!                                      max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
!                                      error=error)
!       IF (safe_mode) THEN
!          CALL cp_dbcsr_init(matrix_tmp1,error=error)
!          CALL cp_dbcsr_create(matrix_tmp1,template=matrix_sigma_vv_full,&
!                               matrix_type=dbcsr_type_no_symmetry,error=error) 
!          CALL cp_dbcsr_init(matrix_tmp2,error=error)
!          CALL cp_dbcsr_create(matrix_tmp2,template=matrix_sigma_vv_full,&
!                               matrix_type=dbcsr_type_no_symmetry,error=error) 
!      
!          CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_sigma_vv_full_sqrt_inv,&
!                                 matrix_sigma_vv_full,&
!                                 0.0_dp,matrix_tmp1,filter_eps=almo_scf_env%eps_filter,&
!                                 error=error)
!          CALL cp_dbcsr_multiply("N","N",1.0_dp,matrix_tmp1,&
!                                 matrix_sigma_vv_full_sqrt_inv,&
!                                 0.0_dp,matrix_tmp2,filter_eps=almo_scf_env%eps_filter,&
!                                 error=error)
!      
!          frob_matrix_base=cp_dbcsr_frobenius_norm(matrix_tmp2)
!          CALL cp_dbcsr_add_on_diag(matrix_tmp2,-1.0_dp,error=error)
!          frob_matrix=cp_dbcsr_frobenius_norm(matrix_tmp2)
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,*) "Error for (inv(sqrt(SIGVV))*SIGVV*inv(sqrt(SIGVV))-I)",frob_matrix/frob_matrix_base
!          ENDIF
!      
!          CALL cp_dbcsr_release(matrix_tmp1,error=error) 
!          CALL cp_dbcsr_release(matrix_tmp2,error=error) 
!       ENDIF
!     
!       ! discard unnecessary overlap functions
!       CALL cp_dbcsr_release(matrix_sigma_vv_full,error=error)
!       CALL cp_dbcsr_release(matrix_sigma_vv_full_sqrt,error=error)
!
!! this can be re-written because we have (1-P)|v>
!
!       !!!!!!!!!!!!!!!!!!!
!       ! 1. Compute F_ov
!       CALL cp_dbcsr_init(Fon,error=error)
!       CALL cp_dbcsr_create(Fon,&
!               template=almo_scf_env%matrix_v_full_blk(ispin),&
!               error=error)
!       CALL cp_dbcsr_init(Fov,error=error)
!       CALL cp_dbcsr_create(Fov,&
!               template=almo_scf_env%matrix_ov_full(ispin),&
!               error=error)
!       CALL cp_dbcsr_init(Fov_filtered,error=error)
!       CALL cp_dbcsr_create(Fov_filtered,&
!               template=almo_scf_env%matrix_ov_full(ispin),&
!               error=error)
!       CALL cp_dbcsr_init(temp1_oo,error=error)
!       CALL cp_dbcsr_create(temp1_oo,&
!               template=almo_scf_env%matrix_sigma(ispin),&
!               !matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!       CALL cp_dbcsr_init(temp2_oo,error=error)
!       CALL cp_dbcsr_create(temp2_oo,&
!               template=almo_scf_env%matrix_sigma(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!
!       CALL cp_dbcsr_multiply("T","N",1.0_dp,almo_scf_env%matrix_t_blk(ispin),&
!               almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
!               0.0_dp,Fon,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_v_full_blk(ispin),&
!               0.0_dp,Fov,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_t_blk(ispin),&
!               0.0_dp,temp1_oo,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,temp1_oo,&
!               almo_scf_env%matrix_sigma_inv(ispin),&
!               0.0_dp,temp2_oo,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!       CALL cp_dbcsr_release(temp1_oo,error=error)
!
!       CALL cp_dbcsr_multiply("T","N",1.0_dp,almo_scf_env%matrix_t_blk(ispin),&
!               almo_scf_env%matrix_s(1),&
!               0.0_dp,Fon,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,Fon,&
!               almo_scf_env%matrix_v_full_blk(ispin),&
!               0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!       CALL cp_dbcsr_release(Fon,error=error)
!
!       CALL cp_dbcsr_multiply("N","N",-1.0_dp,temp2_oo,&
!               Fov_filtered,&
!               1.0_dp,Fov,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!       CALL cp_dbcsr_release(temp2_oo,error=error)
!
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_sigma_inv(ispin),&
!               Fov,0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!       
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,Fov_filtered,&
!               matrix_sigma_vv_full_sqrt_inv,&
!               0.0_dp,Fov,filter_eps=almo_scf_env%eps_filter,&
!               error=error)
!       !CALL cp_dbcsr_copy(Fov,Fov_filtered,error=error)
!CALL cp_dbcsr_print(Fov,error=error)
!       
!       IF (safe_mode) THEN
!          CALL cp_dbcsr_init(Fov_original,error=error)
!          CALL cp_dbcsr_create(Fov_original,template=Fov,error=error)
!          CALL cp_dbcsr_copy(Fov_original,Fov,error=error)
!       ENDIF
!
!!! remove diagonal blocks
!!CALL cp_dbcsr_iterator_start(iter,Fov)
!!DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
!!
!!   CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!!           row_size=iblock_row_size,col_size=iblock_col_size)
!!   
!!   IF (iblock_row.eq.iblock_col) data_p(:,:)=0.0_dp
!!
!!ENDDO
!!CALL cp_dbcsr_iterator_stop(iter)
!!CALL cp_dbcsr_finalize(Fov,error=error)       
!
!!! perform svd of blocks
!!!!! THIS ROUTINE WORKS ONLY ON ONE CPU AND ONLY FOR 2 MOLECULES !!!
!!CALL cp_dbcsr_init(temp_u_v_full_blk,error=error)
!!CALL cp_dbcsr_create(temp_u_v_full_blk,&
!!        template=almo_scf_env%matrix_vv_full_blk(ispin),&
!!        matrix_type=dbcsr_type_no_symmetry,&
!!        error=error)
!!CALL cp_dbcsr_work_create(temp_u_v_full_blk,&
!!        work_mutable=.TRUE.,error=error)
!!CALL cp_dbcsr_iterator_start(iter,Fov)
!!DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
!!
!!   CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!!           row_size=iblock_row_size,col_size=iblock_col_size)
!!
!!   IF (iblock_row.ne.iblock_col) THEN
!!
!!      ! Prepare data
!!      allocate(eigenvalues(min(iblock_row_size,iblock_col_size)))
!!      allocate(data_copy(iblock_row_size,iblock_col_size))
!!      allocate(left_vectors(iblock_row_size,iblock_row_size))
!!      allocate(right_vectors(iblock_col_size,iblock_col_size))
!!      data_copy(:,:)=data_p(:,:)
!!
!!      ! Query the optimal workspace for dgesvd
!!      LWORK = -1
!!      allocate(WORK(MAX(1,LWORK)))
!!      CALL DGESVD('N','A',iblock_row_size,iblock_col_size,data_copy,&
!!              iblock_row_size,eigenvalues,left_vectors,iblock_row_size,&
!!              right_vectors,iblock_col_size,WORK,LWORK,INFO)
!!      LWORK = INT(WORK( 1 ))
!!      deallocate(WORK)
!!
!!      ! Allocate the workspace and perform svd
!!      allocate(WORK(MAX(1,LWORK)))
!!      CALL DGESVD('N','A',iblock_row_size,iblock_col_size,data_copy,&
!!              iblock_row_size,eigenvalues,left_vectors,iblock_row_size,&
!!              right_vectors,iblock_col_size,WORK,LWORK,INFO)
!!      deallocate(WORK)
!!      IF( INFO.NE.0 ) THEN
!!         CPErrorMessage(cp_failure_level,routineP,"DGESVD failed",error)
!!         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!!      END IF
!!
!!      ! copy right singular vectors into a unitary matrix
!!      NULLIFY (p_new_block)
!!      CALL cp_dbcsr_reserve_block2d(temp_u_v_full_blk,iblock_col,iblock_col,p_new_block)
!!      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
!!      p_new_block(:,:) = right_vectors(:,:)
!!
!!      deallocate(eigenvalues)
!!      deallocate(data_copy)
!!      deallocate(left_vectors)
!!      deallocate(right_vectors)
!!
!!   ENDIF
!!ENDDO
!!CALL cp_dbcsr_iterator_stop(iter)
!!CALL cp_dbcsr_finalize(temp_u_v_full_blk,error=error)       
!!!CALL cp_dbcsr_print(temp_u_v_full_blk,error=error)       
!!CALL cp_dbcsr_multiply("N","T",1.0_dp,Fov,temp_u_v_full_blk,&
!!        0.0_dp,Fov_filtered,filter_eps=almo_scf_env%eps_filter,&
!!        error=error)
!!CALL cp_dbcsr_copy(Fov,Fov_filtered,error=error)
!!CALL cp_dbcsr_print(Fov,error=error)       
!       
!       !!!!!!!!!!!!!!!!!!!
!       ! 2. Initialize variables
!
!       ! temp space 
!       CALL cp_dbcsr_init(temp0_ov,error=error)
!       CALL cp_dbcsr_create(temp0_ov,&
!               template=almo_scf_env%matrix_ov_full(ispin),&
!               error=error)
!
!       ! current unitary matrix
!       CALL cp_dbcsr_init(U_blk,error=error)
!       CALL cp_dbcsr_create(U_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!       
!       ! unitary matrix accumulator
!       CALL cp_dbcsr_init(U_blk_tot,error=error)
!       CALL cp_dbcsr_create(U_blk_tot,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!       CALL cp_dbcsr_add_on_diag(U_blk_tot,1.0_dp,error=error)
!
!!CALL cp_dbcsr_add_on_diag(U_blk,1.0_dp,error=error)
!!CALL cp_dbcsr_multiply("N","T",1.0_dp,U_blk,temp_u_v_full_blk,&
!!        0.0_dp,U_blk_tot,filter_eps=almo_scf_env%eps_filter,&
!!        error=error)
!!CALL cp_dbcsr_release(temp_u_v_full_blk,error=error)
!
!       ! init gradient
!       CALL cp_dbcsr_init(grad_blk,error=error)
!       CALL cp_dbcsr_create(grad_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!       
!       ! init step matrix
!       CALL cp_dbcsr_init(step_blk,error=error)
!       CALL cp_dbcsr_create(step_blk,&
!               template=almo_scf_env%matrix_vv_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!
!       ! "retain discarded" filter (0.0 - retain, 1.0 - discard)
!       CALL cp_dbcsr_init(matrix_filter,error=error)
!       CALL cp_dbcsr_create(matrix_filter,&
!               template=almo_scf_env%matrix_ov_full(ispin),&
!               error=error)
!       ! copy Fov into the filter matrix temporarily
!       ! so we know which blocks contain significant elements
!       CALL cp_dbcsr_copy(matrix_filter,Fov,error=error)
!
!       ! fill out filter elements block-by-block
!       CALL cp_dbcsr_iterator_start(iter,matrix_filter)
!       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
!
!          CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!                  row_size=iblock_row_size,col_size=iblock_col_size)
!          
!          retained_v=almo_scf_env%nvirt_of_domain(iblock_col,ispin)
!
!          data_p(:,1:retained_v)=0.0_dp
!          data_p(:,(retained_v+1):iblock_col_size)=1.0_dp
!
!       ENDDO
!       CALL cp_dbcsr_iterator_stop(iter)
!       CALL cp_dbcsr_finalize(matrix_filter,error=error)       
!       
!       ! apply the filter
!       CALL cp_dbcsr_hadamard_product(Fov,matrix_filter,Fov_filtered,error=error)
!  
!       !!!!!!!!!!!!!!!!!!!!!
!       ! 3. start iterative minimization of the elements to be discarded
!       iteration=0
!       converged=.FALSE.
!       prepare_to_exit=.FALSE.
!       DO
!   
!          iteration=iteration+1
!   
!          !!!!!!!!!!!!!!!!!!!!!!!!!
!          ! 4. compute the gradient
!          CALL cp_dbcsr_set(grad_blk,0.0_dp,error)
!          ! create the diagonal blocks only
!          CALL cp_dbcsr_add_on_diag(grad_blk,1.0_dp,error=error)
!          
!          CALL cp_dbcsr_multiply("T","N",2.0_dp,Fov_filtered,Fov,&
!                  0.0_dp,grad_blk,retain_sparsity=.TRUE.,&
!                  filter_eps=almo_scf_env%eps_filter,error=error)
!          CALL cp_dbcsr_multiply("T","N",-2.0_dp,Fov,Fov_filtered,&
!                  1.0_dp,grad_blk,retain_sparsity=.TRUE.,&
!                  filter_eps=almo_scf_env%eps_filter,error=error)
!
!          !!!!!!!!!!!!!!!!!!!!!!!
!          ! 5. check convergence
!          obj_function = 0.5_dp*(cp_dbcsr_frobenius_norm(Fov_filtered))**2
!          grad_norm = cp_dbcsr_frobenius_norm(grad_blk)
!          converged=(grad_norm.lt.almo_scf_env%truncate_v_eps_convergence)
!          IF (converged.OR.(iteration.ge.almo_scf_env%truncate_v_max_iter)) THEN
!             prepare_to_exit=.TRUE.
!          ENDIF
!
!          IF (.NOT.prepare_to_exit) THEN
!
!             !!!!!!!!!!!!!!!!!!!!!!!
!             ! 6. perform steps in the direction of the gradient
!             !    a. first, perform a trial step to "see" the parameters
!             !       of the parabola along the gradient:
!             !       a0 * x^2 + b0 * x + c0
!             !    b. then perform the step to the bottom of the parabola
!
!             ! get c0
!             c0 = obj_function
!             ! get b0 <= d_f/d_alpha along grad
!             !!!CALL cp_dbcsr_multiply("N","N",4.0_dp,Fov,grad_blk,&
!             !!!        0.0_dp,temp0_ov,&
!             !!!        filter_eps=almo_scf_env%eps_filter,error=error)
!             !!!CALL cp_dbcsr_trace(Fov_filtered,temp0_ov,b0,"T","N",error=error)
!
!             alpha=almo_scf_env%truncate_v_trial_step_size
!           
!             line_search_step_last=3
!             DO line_search_step=1,line_search_step_last
!                CALL cp_dbcsr_copy(step_blk,grad_blk,error=error)
!                CALL cp_dbcsr_scale(step_blk,-1.0_dp*alpha,error=error)
!                CALL generator_to_unitary(step_blk,U_blk,&
!                        almo_scf_env%eps_filter,error)
!                CALL cp_dbcsr_multiply("N","N",1.0_dp,Fov,U_blk,0.0_dp,temp0_ov,&
!                        filter_eps=almo_scf_env%eps_filter,error=error)
!                CALL cp_dbcsr_hadamard_product(temp0_ov,matrix_filter,&
!                        Fov_filtered,error=error)
!               
!                obj_function_new = 0.5_dp*(cp_dbcsr_frobenius_norm(Fov_filtered))**2
!                IF (line_search_step.eq.1) THEN
!                   ff1 = obj_function_new
!                   step1 = alpha
!                ELSE IF (line_search_step.eq.2) THEN
!                   ff2 = obj_function_new
!                   step2 = alpha
!                ENDIF
!
!                IF (unit_nr>0.AND.(line_search_step.ne.line_search_step_last)) THEN
!                   WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.3,E12.3,E12.3,E12.3)') &
!                         "JOINT_SVD_lin",&
!                         iteration,&
!                         alpha,&
!                         obj_function,&
!                         obj_function_new,&
!                         obj_function_new-obj_function
!                ENDIF
!
!                IF (line_search_step.eq.1) THEN 
!                   alpha=2.0_dp*alpha
!                ENDIF
!                IF (line_search_step.eq.2) THEN 
!                   a0 = ((ff1-c0)/step1 - (ff2-c0)/step2) / (step1 - step2)
!                   b0 = (ff1-c0)/step1 - a0*step1
!                   ! step size in to the bottom of "the parabola" 
!                   alpha=-b0/(2.0_dp*a0)
!                   ! update the default step size
!                   almo_scf_env%truncate_v_trial_step_size=alpha
!                ENDIF
!                !!!IF (line_search_step.eq.1) THEN 
!                !!!   a0 = (obj_function_new - b0 * alpha - c0) / (alpha*alpha)
!                !!!   ! step size in to the bottom of "the parabola" 
!                !!!   alpha=-b0/(2.0_dp*a0)
!                !!!   !IF (alpha.gt.10.0_dp) alpha=10.0_dp
!                !!!ENDIF
!
!             ENDDO
!
!             ! update Fov and U_blk_tot (use grad_blk as tmp storage)
!             CALL cp_dbcsr_copy(Fov,temp0_ov,error=error)
!             CALL cp_dbcsr_multiply("N","N",1.0_dp,U_blk_tot,U_blk,&
!                     0.0_dp,grad_blk,&
!                     filter_eps=almo_scf_env%eps_filter,error=error)
!             CALL cp_dbcsr_copy(U_blk_tot,grad_blk,error=error)
!
!          ENDIF
!       
!          t2 = m_walltime()
!          
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,'(T6,A,1X,I3,1X,F10.3,E12.3,E12.3,E12.3,E12.3,F10.3)') &
!                   "JOINT_SVD_itr",&
!                   iteration,&
!                   alpha,&
!                   obj_function,&
!                   obj_function_new,&
!                   obj_function_new-obj_function,&
!                   grad_norm,&
!                   t2-t1
!                   !(flop1+flop2)/(1.0E6_dp*(t2-t1))
!             CALL m_flush(unit_nr)
!          ENDIF
!
!          t1 = m_walltime()
!
!          IF (prepare_to_exit) EXIT
!
!       ENDDO ! stop iterations
!
!       IF (safe_mode) THEN
!          CALL cp_dbcsr_multiply("N","N",1.0_dp,Fov_original,&
!                  U_blk_tot,0.0_dp,temp0_ov,&
!                  filter_eps=almo_scf_env%eps_filter,error=error)
!CALL cp_dbcsr_print(temp0_ov,error=error)
!          CALL cp_dbcsr_hadamard_product(temp0_ov,matrix_filter,&
!                  Fov_filtered,error=error)
!          obj_function_new = 0.5_dp*(cp_dbcsr_frobenius_norm(Fov_filtered))**2
!          
!          IF (unit_nr>0) THEN
!             WRITE(unit_nr,'(T6,A,1X,E12.3)') &
!                   "SANITY CHECK:",&
!                   obj_function_new
!             CALL m_flush(unit_nr)
!          ENDIF
!
!          CALL cp_dbcsr_release(Fov_original,error)
!       ENDIF
!       
!       CALL cp_dbcsr_release(temp0_ov,error=error)
!       CALL cp_dbcsr_release(U_blk,error=error)
!       CALL cp_dbcsr_release(grad_blk,error=error)
!       CALL cp_dbcsr_release(step_blk,error=error)
!       CALL cp_dbcsr_release(matrix_filter,error=error)
!       CALL cp_dbcsr_release(Fov,error=error)
!       CALL cp_dbcsr_release(Fov_filtered,error=error)
!
!       ! compute rotated virtual orbitals
!       CALL cp_dbcsr_init(v_full_tmp,error=error)
!       CALL cp_dbcsr_create(v_full_tmp,&
!               template=almo_scf_env%matrix_v_full_blk(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!               v_full_new,&
!               matrix_sigma_vv_full_sqrt_inv,0.0_dp,v_full_tmp,&
!               filter_eps=almo_scf_env%eps_filter,error=error)
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!               v_full_tmp,&
!               U_blk_tot,0.0_dp,v_full_new,&
!               filter_eps=almo_scf_env%eps_filter,error=error)
!       
!       CALL cp_dbcsr_release(matrix_sigma_vv_full_sqrt_inv,error=error)
!       CALL cp_dbcsr_release(v_full_tmp,error=error)
!       CALL cp_dbcsr_release(U_blk_tot,error=error)
!
!!!!! orthogonalized virtuals are not blocked       
!       ! copy new virtuals into the truncated matrix
!       !CALL cp_dbcsr_work_create(almo_scf_env%matrix_v_blk(ispin),&
!       CALL cp_dbcsr_work_create(almo_scf_env%matrix_v(ispin),&
!               work_mutable=.TRUE.,error=error)
!       CALL cp_dbcsr_iterator_start(iter,v_full_new)
!       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
!
!          CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,&
!                  row_size=iblock_row_size,col_size=iblock_col_size)
!          
!          retained_v=almo_scf_env%nvirt_of_domain(iblock_col,ispin)
!
!          NULLIFY (p_new_block)
!          !CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_v_blk(ispin),&
!          CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_v(ispin),&
!                  iblock_row,iblock_col,p_new_block)
!          CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
!          CPPrecondition(retained_v.gt.0,cp_failure_level,routineP,error,failure)
!          p_new_block(:,:) = data_p(:,1:retained_v)
!
!       ENDDO ! iterator
!       CALL cp_dbcsr_iterator_stop(iter)
!       !!CALL cp_dbcsr_finalize(almo_scf_env%matrix_v_blk(ispin),error=error)       
!       CALL cp_dbcsr_finalize(almo_scf_env%matrix_v(ispin),error=error)       
!       
!       CALL cp_dbcsr_release(v_full_new,error=error)
!
!    ENDDO ! ispin
!
!    CALL timestop(handle)
!  
!  END SUBROUTINE truncate_subspace_v_blk

END MODULE almo_scf

