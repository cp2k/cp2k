!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for all ALMO-based SCF methods
!>        RZK-warning marks unresolved issues
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin 
! *****************************************************************************
MODULE almo_scf
  USE almo_scf_methods,                ONLY: almo_scf_t_blk_to_t_blk_orthonormal,&
                                             distribute_domains
  USE almo_scf_optimizer,              ONLY: almo_scf_block_diagonal,&
                                             almo_scf_xalmo_eigensolver,&
                                             almo_scf_xalmo_pcg
  USE almo_scf_qs,                     ONLY: almo_scf_construct_quencher,&
                                             almo_scf_init_qs,&
                                             matrix_almo_create,&
                                             matrix_qs_to_almo
  USE almo_scf_types,                  ONLY: &
       almo_mat_dim_aobasis, almo_mat_dim_occ, almo_mat_dim_virt, &
       almo_mat_dim_virt_disc, almo_mat_dim_virt_full, &
       almo_max_cutoff_multiplier, almo_scf_env_type, optimizer_options_type, &
       print_optimizer_options
  USE bibliography,                    ONLY: Khaliullin2013,&
                                             cite_reference
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add_on_diag, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_distribution, cp_dbcsr_filter, cp_dbcsr_finalize, &
       cp_dbcsr_get_info, cp_dbcsr_get_stored_coordinates, cp_dbcsr_init, &
       cp_dbcsr_nblkcols_total, cp_dbcsr_nblkrows_total, cp_dbcsr_p_type, &
       cp_dbcsr_release, cp_dbcsr_reserve_block2d, cp_dbcsr_set, &
       cp_dbcsr_type, cp_dbcsr_work_create, dbcsr_distribution_mp, &
       dbcsr_mp_mynode, dbcsr_type_no_symmetry, dbcsr_type_symmetric
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE domain_submatrix_methods,        ONLY: init_submatrices,&
                                             release_submatrices
  USE input_constants,                 ONLY: &
       almo_constraint_distance, almo_deloc_none, almo_deloc_scf, &
       almo_deloc_x, almo_deloc_x_then_scf, almo_deloc_xalmo_1diag, &
       almo_deloc_xalmo_scf, almo_deloc_xalmo_x, almo_deloc_xk, &
       almo_domain_layout_atomic, almo_domain_layout_molecular, &
       almo_mat_distr_atomic, almo_mat_distr_molecular, almo_scf_diag, &
       almo_scf_dm_sign, almo_scf_pcg, cg_hager_zhang, do_bondparm_vdw, &
       optimizer_diis, optimizer_pcg, tensor_orthogonal, virt_full, &
       virt_minimal, virt_number, virt_occ_size, xalmo_case_block_diag, &
       xalmo_case_fully_deloc, xalmo_case_normal
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                             matrix_sqrt_Newton_Schulz
  USE kinds,                           ONLY: dp
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
    CALL almo_scf_clean_up(almo_scf_env,error)

    CALL timestop(handle)

  END SUBROUTINE almo_entry_scf

! *****************************************************************************
!> \brief Initialization of the almo_scf_env_type.
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param error ...
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

    INTEGER :: ao, handle, i, iao, idomain, ispin, multip, naos, natoms, &
      ndomains, nelec, nelec_a, nelec_b, nmols, nspins, unit_nr
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

    ! set optimizers' types
    almo_scf_env%opt_block_diag_diis%optimizer_type=optimizer_diis
    almo_scf_env%opt_block_diag_pcg%optimizer_type=optimizer_pcg
    almo_scf_env%opt_xalmo_diis%optimizer_type=optimizer_diis
    almo_scf_env%opt_xalmo_pcg%optimizer_type=optimizer_pcg
    
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
    CALL almo_scf_init_read_write_input(input,almo_scf_env,error)
   
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
    ALLOCATE(almo_scf_env%last_atom_of_domain(ndomains))
    ALLOCATE(almo_scf_env%nbasis_of_domain(ndomains))
    ALLOCATE(almo_scf_env%nocc_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%nvirt_full_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%nvirt_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%nvirt_disc_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%mu_of_domain(ndomains,nspins))
    ALLOCATE(almo_scf_env%cpu_of_domain(ndomains))
    ALLOCATE(almo_scf_env%charge_of_domain(ndomains))
    ALLOCATE(almo_scf_env%multiplicity_of_domain(ndomains))

    ! fill out domain descriptors and group descriptors
    IF (almo_scf_env%domain_layout_mos==almo_domain_layout_molecular) THEN
       ! get domain info from molecule_set 
       CALL get_molecule_set_info(molecule_set,&
               atom_to_mol=almo_scf_env%domain_index_of_atom,&
               mol_to_first_atom=almo_scf_env%first_atom_of_domain,&
               mol_to_last_atom=almo_scf_env%last_atom_of_domain,&
               mol_to_nelectrons=almo_scf_env%nocc_of_domain(1:ndomains,1),&
               mol_to_nbasis=almo_scf_env%nbasis_of_domain,&
               mol_to_charge=almo_scf_env%charge_of_domain,&
               mol_to_multiplicity=almo_scf_env%multiplicity_of_domain,&
               error=error)
       ! calculate number of alpha and beta occupied orbitals from
       ! the number of electrons and multiplicity of each molecule
       ! Na + Nb = Ne
       ! Na - Nb = Mult - 1 (assume Na > Nb as we do not have more info from get_molecule_set_info)
       DO idomain=1,ndomains
          nelec=almo_scf_env%nocc_of_domain(idomain,1)
          multip=almo_scf_env%multiplicity_of_domain(idomain)
          nelec_a=(nelec+multip-1)/2
          nelec_b=nelec-nelec_a
          almo_scf_env%nocc_of_domain(idomain,1)=nelec_a
          IF (nelec_a.ne.nelec_b) THEN 
             IF (nspins.eq.1) THEN
                WRITE(*,*) "Domain ", idomain, " out of ", ndomains, ". Electrons = ", nelec
                CPErrorMessage(cp_failure_level,routineP,"odd e- -- use unrestricted methods",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF
             almo_scf_env%nocc_of_domain(idomain,2)=nelec_b
             ! RZK-warning: open-shell procedures have not been tested yet
             ! Stop the program now
             CPErrorMessage(cp_failure_level,routineP,"Unrestricted ALMO methods are NYI",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          ENDIF
       ENDDO
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

    ! create all matrices
    CALL almo_scf_env_create_matrices(almo_scf_env,matrix_s(1)%matrix,error)
    ! set up matrix S and all required functions of S
    almo_scf_env%s_inv_done=.FALSE.
    almo_scf_env%s_sqrt_done=.FALSE.
    CALL almo_scf_init_ao_overlap(matrix_s(1)%matrix,almo_scf_env,error)

    ! create the quencher (imposes sparsity template)
    CALL almo_scf_construct_quencher(qs_env,almo_scf_env,error)
    CALL distribute_domains(almo_scf_env,error)

    ! FINISH setting job parameters here, print out job info
    CALL almo_scf_print_job_info(almo_scf_env,unit_nr,error)

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
!> \param input ...
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_init_read_write_input(input,almo_scf_env,error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'almo_scf_init_read_write_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER :: almo_opt_diis_section, &
      almo_opt_pcg_section, almo_scf_section, xalmo_opt_pcg_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    
    almo_scf_section => section_vals_get_subs_vals(input,"DFT%ALMO_SCF",error=error)
    almo_opt_diis_section => section_vals_get_subs_vals(almo_scf_section,&
                             "ALMO_OPTIMIZER_DIIS",error=error)
    almo_opt_pcg_section => section_vals_get_subs_vals(almo_scf_section,&
                            "ALMO_OPTIMIZER_PCG",error=error)
    xalmo_opt_pcg_section => section_vals_get_subs_vals(almo_scf_section,&
                             "XALMO_OPTIMIZER_PCG",error=error)
   
    ! RZK-warning the values of these keywords are hardcoded 
    ! but ideally they should also be read from the input file
    almo_scf_env%order_lanczos=3
    almo_scf_env%eps_lanczos=1.0E-4_dp
    almo_scf_env%max_iter_lanczos=40

    ! read user input
    ! common ALMO options
    CALL section_vals_val_get(almo_scf_section,"EPS_FILTER",&
                              r_val=almo_scf_env%eps_filter,error=error)
    CALL section_vals_val_get(almo_scf_section,"ALMO_SCF_GUESS",&
                              i_val=almo_scf_env%almo_scf_guess,&
                              error=error)
    CALL section_vals_val_get(almo_scf_section,"ALMO_ALGORITHM",&
                              i_val=almo_scf_env%almo_update_algorithm,&
                              error=error)
    CALL section_vals_val_get(almo_scf_section,"DELOCALIZE_METHOD",&
                              i_val=almo_scf_env%deloc_method,error=error)
    CALL section_vals_val_get(almo_scf_section,"XALMO_R_CUTOFF_FACTOR",&
                              r_val=almo_scf_env%quencher_r0_factor,&
                              error=error)

    ! optimizers
    CALL section_vals_val_get(almo_opt_diis_section,"EPS_ERROR",&
            r_val=almo_scf_env%opt_block_diag_diis%eps_error,&
            error=error)
    CALL section_vals_val_get(almo_opt_diis_section,"MAX_ITER",&
            i_val=almo_scf_env%opt_block_diag_diis%max_iter,error=error)
    CALL section_vals_val_get(almo_opt_diis_section,"N_DIIS",&
            i_val=almo_scf_env%opt_block_diag_diis%ndiis,error=error)

    CALL section_vals_val_get(almo_opt_pcg_section,"EPS_ERROR",&
            r_val=almo_scf_env%opt_block_diag_pcg%eps_error,&
            error=error)
    CALL section_vals_val_get(almo_opt_pcg_section,"MAX_ITER",&
            i_val=almo_scf_env%opt_block_diag_pcg%max_iter,error=error)
    CALL section_vals_val_get(almo_opt_pcg_section,"MAX_ITER_OUTER_LOOP",&
            i_val=almo_scf_env%opt_block_diag_pcg%max_iter_outer_loop,error=error)
    CALL section_vals_val_get(almo_opt_pcg_section,"LIN_SEARCH_EPS_ERROR",&
            r_val=almo_scf_env%opt_block_diag_pcg%lin_search_eps_error,error=error)
    CALL section_vals_val_get(almo_opt_pcg_section,"LIN_SEARCH_STEP_SIZE_GUESS",&
            r_val=almo_scf_env%opt_block_diag_pcg%lin_search_step_size_guess,error=error)
    CALL section_vals_val_get(almo_opt_pcg_section,"CONJUGATOR",&
            i_val=almo_scf_env%opt_block_diag_pcg%conjugator,error=error)
    CALL section_vals_val_get(almo_opt_pcg_section,"PRECONDITIONER",&
            i_val=almo_scf_env%opt_block_diag_pcg%preconditioner,error=error)
    
    CALL section_vals_val_get(xalmo_opt_pcg_section,"EPS_ERROR",&
            r_val=almo_scf_env%opt_xalmo_pcg%eps_error,&
            error=error)
    CALL section_vals_val_get(xalmo_opt_pcg_section,"MAX_ITER",&
            i_val=almo_scf_env%opt_xalmo_pcg%max_iter,error=error)
    CALL section_vals_val_get(xalmo_opt_pcg_section,"MAX_ITER_OUTER_LOOP",&
            i_val=almo_scf_env%opt_xalmo_pcg%max_iter_outer_loop,error=error)
    CALL section_vals_val_get(xalmo_opt_pcg_section,"LIN_SEARCH_EPS_ERROR",&
            r_val=almo_scf_env%opt_xalmo_pcg%lin_search_eps_error,error=error)
    CALL section_vals_val_get(xalmo_opt_pcg_section,"LIN_SEARCH_STEP_SIZE_GUESS",&
            r_val=almo_scf_env%opt_xalmo_pcg%lin_search_step_size_guess,error=error)
    CALL section_vals_val_get(xalmo_opt_pcg_section,"CONJUGATOR",&
            i_val=almo_scf_env%opt_xalmo_pcg%conjugator,error=error)
    CALL section_vals_val_get(xalmo_opt_pcg_section,"PRECONDITIONER",&
            i_val=almo_scf_env%opt_xalmo_pcg%preconditioner,error=error)
    
    ! do not do EDA in the official release until the output files are written
    ! using the proper CP2K routines
    almo_scf_env%almo_eda=0

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

    !! how to read lists 
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

    almo_scf_env%domain_layout_aos=almo_domain_layout_molecular
    almo_scf_env%domain_layout_mos=almo_domain_layout_molecular
    almo_scf_env%mat_distr_aos=almo_mat_distr_molecular 
    almo_scf_env%mat_distr_mos=almo_mat_distr_molecular

    almo_scf_env%constraint_type=almo_constraint_distance
    almo_scf_env%mu=-0.1_dp
    almo_scf_env%fixed_mu=.FALSE.
    almo_scf_env%eps_prev_guess=almo_scf_env%eps_filter/100.0_dp
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

    almo_scf_env%integer01=10 ! start eDOF-md
    almo_scf_env%integer02=4  ! preconditioner type
    almo_scf_env%integer03=0  ! not used
    almo_scf_env%integer04=0  ! fixed number of line searches (no grad)
    almo_scf_env%integer05=0  ! not used
       
    ! select options that are not directly provided 
    ! but are inderectly determined by the input
    IF (almo_scf_env%deloc_method.eq.almo_deloc_xalmo_1diag) THEN
        almo_scf_env%xalmo_update_algorithm=almo_scf_diag
    ELSE
        almo_scf_env%xalmo_update_algorithm=almo_scf_pcg
    ENDIF

    ! check for conflicts between options
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
       CPErrorMessage(cp_failure_level,routineP,"XALMO_R_CUTOFF_FACTOR is larger than almo_max_cutoff_multiplier",error)
       CPErrorMessage(cp_failure_level,routineP,"increase the hard-coded almo_max_cutoff_multiplier",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_scf_init_read_write_input

! *****************************************************************************
!> \brief Prints out a short summary about the ALMO SCF job
!> \param almo_scf_env ...
!> \param unit_nr ...
!> \param error ...
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

    CHARACTER(len=13)                        :: neig_string
    CHARACTER(len=33)                        :: deloc_method_string
    INTEGER                                  :: handle, idomain, index1_prev

    CALL timeset(routineN,handle)

    IF (unit_nr>0) THEN
       WRITE(unit_nr,'()')
       WRITE(unit_nr,'(T2,A,A,A)') REPEAT("-",32)," ALMO SETTINGS ",REPEAT("-",32)

       WRITE(unit_nr,'(T2,A,T48,E33.3)') "eps_filter:",almo_scf_env%eps_filter
   
       WRITE(unit_nr,'(T2,A)') "optimization of block-diagonal ALMOs:"
       SELECT CASE(almo_scf_env%almo_update_algorithm)
       CASE(almo_scf_diag)
          ! the DIIS algorith is the only choice for the diagonlaization-based algorithm
          CALL print_optimizer_options(almo_scf_env%opt_block_diag_diis,unit_nr,error=error)
       CASE(almo_scf_pcg)
          ! print out PCG options
          CALL print_optimizer_options(almo_scf_env%opt_block_diag_pcg,unit_nr,error=error)
       END SELECT
       
       SELECT CASE(almo_scf_env%deloc_method)
       CASE(almo_deloc_none)
          deloc_method_string="NONE"
       CASE(almo_deloc_x)
          deloc_method_string="FULL_X"
       CASE(almo_deloc_scf)
          deloc_method_string="FULL_SCF"
       CASE(almo_deloc_x_then_scf)
          deloc_method_string="FULL_X_THEN_SCF"
       CASE(almo_deloc_xalmo_1diag)
          deloc_method_string="XALMO_1DIAG"
       CASE(almo_deloc_xalmo_x)
          deloc_method_string="XALMO_X"
       CASE(almo_deloc_xalmo_scf)
          deloc_method_string="XALMO_SCF"
       END SELECT
       WRITE(unit_nr,'(T2,A,T48,A33)') "delocalization:",TRIM(deloc_method_string)
       
       IF (almo_scf_env%deloc_method.ne.almo_deloc_none) THEN
          
          SELECT CASE(almo_scf_env%deloc_method)
          CASE(almo_deloc_x,almo_deloc_scf,almo_deloc_x_then_scf)
             WRITE(unit_nr,'(T2,A,T48,A33)') "delocalization cutoff radius:",&
                                               "infinite"
              deloc_method_string="FULL_X_THEN_SCF"
          CASE(almo_deloc_xalmo_1diag,almo_deloc_xalmo_x,almo_deloc_xalmo_scf)
             WRITE(unit_nr,'(T2,A,T48,F33.5)') "XALMO cutoff radius:",&
                                               almo_scf_env%quencher_r0_factor
          END SELECT
          
          IF (almo_scf_env%deloc_method.eq.almo_deloc_xalmo_1diag) THEN
             ! print nothing because no actual optimization is done
          ELSE
             WRITE(unit_nr,'(T2,A)') "optimization of extended orbitals:"
             SELECT CASE(almo_scf_env%xalmo_update_algorithm)
             CASE(almo_scf_diag)
                    CALL print_optimizer_options(almo_scf_env%opt_xalmo_diis,unit_nr,error=error)
             CASE(almo_scf_pcg)
                 CALL print_optimizer_options(almo_scf_env%opt_xalmo_pcg,unit_nr,error=error)
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

       WRITE(unit_nr,'(T2,A)') REPEAT("-",79)

       IF (almo_scf_env%ndomains.le.200) THEN

          ! print fragment info
          WRITE(unit_nr,'(T2,A13,A13,A13,A13,A13,A13)') &
             "Fragment","Basis Set","Occupied","Virtual","Charge","Deloc Neig"!,"Discarded Virt"
          WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
          DO idomain=1,almo_scf_env%ndomains

             IF (idomain.eq.1) THEN
                index1_prev = 1
             ELSE
                index1_prev = almo_scf_env%domain_map(1)%index1(idomain-1)
             ENDIF
       
             SELECT CASE(almo_scf_env%deloc_method)
             CASE(almo_deloc_none)
                neig_string="NONE"
             CASE(almo_deloc_x,almo_deloc_scf,almo_deloc_x_then_scf)
                neig_string="ALL"
             CASE(almo_deloc_xalmo_1diag,almo_deloc_xalmo_x,almo_deloc_xalmo_scf)
                WRITE(neig_string,'(I13)') &
                   almo_scf_env%domain_map(1)%index1(idomain)-index1_prev-1 ! minus self
             CASE DEFAULT
                neig_string="N/A"
             END SELECT

             WRITE(unit_nr,'(T2,I13,I13,I13,I13,I13,A13)') &
                idomain,almo_scf_env%nbasis_of_domain(idomain),&
                        SUM(almo_scf_env%nocc_of_domain(idomain,:)),&
                        SUM(almo_scf_env%nvirt_of_domain(idomain,:)),&
                        !SUM(almo_scf_env%nvirt_disc_of_domain(idomain,:)),&
                        almo_scf_env%charge_of_domain(idomain),&
                        ADJUSTR(TRIM(neig_string))

          ENDDO ! cycle over domains
          
          SELECT CASE(almo_scf_env%deloc_method)
          CASE(almo_deloc_xalmo_1diag,almo_deloc_xalmo_x,almo_deloc_xalmo_scf)
             
             WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
             
             ! print fragment neigbors
             WRITE(unit_nr,'(T2,A78)') &
                "Neighbor lists (including self)"
             WRITE(unit_nr,'(T2,A)') REPEAT("-",79)
             DO idomain=1,almo_scf_env%ndomains

                IF (idomain.eq.1) THEN
                   index1_prev = 1
                ELSE
                   index1_prev = almo_scf_env%domain_map(1)%index1(idomain-1)
                ENDIF

                WRITE(unit_nr,'(T2,I13,":")') idomain
                WRITE(unit_nr,'(T15,5I13)') &
                   almo_scf_env%domain_map(1)%pairs&
                      (index1_prev:almo_scf_env%domain_map(1)%index1(idomain)-1,1) ! includes self

             ENDDO ! cycle over domains
          
          END SELECT
          
          WRITE(unit_nr,'(T2,A)') REPEAT("-",79)

       ENDIF 

       WRITE(unit_nr,'()')

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_scf_print_job_info

! *****************************************************************************
!> \brief Initializes the ALMO SCF copy of the AO overlap matrix 
!>        and all necessary functions (sqrt, inverse...)
!> \param matrix_s ...
!> \param almo_scf_env ...
!> \param error ...
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
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param error ...
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
       
       ! ALMO PCG optimizer as a special case of XALMO PCG
       CALL almo_scf_xalmo_pcg(qs_env=qs_env,&
                               almo_scf_env=almo_scf_env,&
                               optimizer=almo_scf_env%opt_block_diag_pcg,&
                               quench_t=almo_scf_env%quench_t_blk,&
                               matrix_t_in=almo_scf_env%matrix_t_blk,&
                               matrix_t_out=almo_scf_env%matrix_t_blk,&
                               assume_t0_q0x=.FALSE.,&
                               perturbation_only=.FALSE.,&
                               special_case=xalmo_case_block_diag,&
                               error=error)

       CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env,error)
    
    ELSE
       
       ! mixing/DIIS optimizer
       CALL almo_scf_block_diagonal(qs_env,almo_scf_env,&
               almo_scf_env%opt_block_diag_diis,error)
    
    ENDIF
    
    ! we might need a copy of the converged KS
    DO ispin=1,almo_scf_env%nspins
       CALL cp_dbcsr_copy(almo_scf_env%matrix_ks_almo_scf_converged(ispin),&
                 almo_scf_env%matrix_ks(ispin),error=error)
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE almo_scf_main

! *****************************************************************************
!> \brief selects various post scf routines
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param error ...
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
    LOGICAL                                  :: almo_experimental, failure, tr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_new_block
    TYPE(cp_dbcsr_type), ALLOCATABLE, &
      DIMENSION(:)                           :: no_quench
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(optimizer_options_type)             :: arbitrary_optimizer

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    ! create a local optimizer that handles XALMO DIIS
    ! the options of this optimizer are arbitrary because
    ! XALMO DIIS SCF does not converge for yet unknown reasons and 
    ! currently used in the code to get perturbative estimates only
    arbitrary_optimizer%optimizer_type=optimizer_diis
    arbitrary_optimizer%max_iter=3
    arbitrary_optimizer%eps_error=1.0E-6_dp
    arbitrary_optimizer%ndiis=2

    SELECT CASE (almo_scf_env%deloc_method)
    CASE (almo_deloc_x,almo_deloc_scf,almo_deloc_x_then_scf)

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
       ! As a matter of fact it is! But this block treats
       ! fully delocalized MOs. So it is unavoidable.
       ! C'est la vie
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
    CASE (almo_deloc_none,almo_deloc_scf)

       DO ispin=1,almo_scf_env%nspins
          CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin),&
                  almo_scf_env%matrix_t_blk(ispin),error=error)
       ENDDO

    CASE (almo_deloc_x,almo_deloc_xk,almo_deloc_x_then_scf)

       !!!! RZK-warning a whole class of delocalization methods
       !!!! are commented out at the moment because some of their
       !!!! routines have not been thoroughly tested.

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

       CALL almo_scf_xalmo_pcg(qs_env=qs_env,&
                               almo_scf_env=almo_scf_env,&
                               optimizer=almo_scf_env%opt_xalmo_pcg,&
                               quench_t=no_quench,&
                               matrix_t_in=almo_scf_env%matrix_t_blk,&
                               matrix_t_out=almo_scf_env%matrix_t,&
                               assume_t0_q0x=.TRUE.,&
                               perturbation_only=.TRUE.,&
                               special_case=xalmo_case_fully_deloc,&
                               error=error)

    CASE (almo_deloc_xalmo_1diag)

       IF (almo_scf_env%xalmo_update_algorithm.eq.almo_scf_diag) THEN

          almo_scf_env%perturbative_delocalization=.TRUE.
          DO ispin=1,almo_scf_env%nspins
             CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin),&
                     almo_scf_env%matrix_t_blk(ispin),error=error)
          ENDDO
          CALL almo_scf_xalmo_eigensolver(qs_env,almo_scf_env,&
                  arbitrary_optimizer,error)
       
       ELSE
          
          CPErrorMessage(cp_failure_level,routineP,"Other algorithms do not exist",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       
       ENDIF

    CASE (almo_deloc_xalmo_x)

       IF (almo_scf_env%xalmo_update_algorithm.eq.almo_scf_pcg) THEN

          CALL almo_scf_xalmo_pcg(qs_env=qs_env,&
                                  almo_scf_env=almo_scf_env,&
                                  optimizer=almo_scf_env%opt_xalmo_pcg,&
                                  quench_t=almo_scf_env%quench_t,&
                                  matrix_t_in=almo_scf_env%matrix_t_blk,&
                                  matrix_t_out=almo_scf_env%matrix_t,&
                                  assume_t0_q0x=.TRUE.,&
                                  perturbation_only=.TRUE.,&
                                  special_case=xalmo_case_normal,&
                                  error=error)
       
       ELSE

          CPErrorMessage(cp_failure_level,routineP,"Other algorithms do not exist",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       
       ENDIF

    CASE (almo_deloc_xalmo_scf)

       IF (almo_scf_env%xalmo_update_algorithm.eq.almo_scf_diag) THEN
          
          CPErrorMessage(cp_failure_level,routineP,"Should not be here: convergence will fail!",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          
          almo_scf_env%perturbative_delocalization=.FALSE.
          DO ispin=1,almo_scf_env%nspins
             CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin),&
                     almo_scf_env%matrix_t_blk(ispin),error=error)
          ENDDO
          CALL almo_scf_xalmo_eigensolver(qs_env,almo_scf_env,&
                  arbitrary_optimizer,error)

       ELSE IF (almo_scf_env%xalmo_update_algorithm.eq.almo_scf_pcg) THEN
          
          CALL almo_scf_xalmo_pcg(qs_env=qs_env,&
                                  almo_scf_env=almo_scf_env,&
                                  optimizer=almo_scf_env%opt_xalmo_pcg,&
                                  quench_t=almo_scf_env%quench_t,&
                                  matrix_t_in=almo_scf_env%matrix_t_blk,&
                                  matrix_t_out=almo_scf_env%matrix_t,&
                                  assume_t0_q0x=.TRUE.,&
                                  perturbation_only=.FALSE.,&
                                  special_case=xalmo_case_normal,&
                                  error=error)
       
          ! RZK-warning THIS IS A HACK TO GET ORBITAL ENERGIES
          almo_experimental=.FALSE.
          IF (almo_experimental) THEN
             almo_scf_env%perturbative_delocalization=.TRUE.
             !DO ispin=1,almo_scf_env%nspins
             !   CALL cp_dbcsr_copy(almo_scf_env%matrix_t(ispin),&
             !           almo_scf_env%matrix_t_blk(ispin),error=error)
             !ENDDO
             CALL almo_scf_xalmo_eigensolver(qs_env,almo_scf_env,&
                     arbitrary_optimizer,error)
          ENDIF ! experimental
       
       ENDIF

    CASE DEFAULT

       CPErrorMessage(cp_failure_level,routineP,"Illegal delocalization method",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    
    END SELECT
    
    SELECT CASE (almo_scf_env%deloc_method)
    CASE (almo_deloc_scf,almo_deloc_x_then_scf)
       
       IF (almo_scf_env%deloc_truncate_virt.ne.virt_full) THEN
          CPErrorMessage(cp_failure_level,routineP,"full scf is NYI for truncated virtual space",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
       
       CALL almo_scf_xalmo_pcg(qs_env=qs_env,&
                               almo_scf_env=almo_scf_env,&
                               optimizer=almo_scf_env%opt_xalmo_pcg,&
                               quench_t=no_quench,&
                               matrix_t_in=almo_scf_env%matrix_t,&
                               matrix_t_out=almo_scf_env%matrix_t,&
                               assume_t0_q0x=.FALSE.,&
                               perturbation_only=.FALSE.,&
                               special_case=xalmo_case_fully_deloc,&
                               error=error)
       
    END SELECT
    
    ! clean up 
    SELECT CASE (almo_scf_env%deloc_method)
    CASE (almo_deloc_x,almo_deloc_scf,almo_deloc_x_then_scf)
       DO ispin=1, almo_scf_env%nspins
          CALL cp_dbcsr_release(no_quench(ispin),error=error)
       ENDDO
       DEALLOCATE(no_quench)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE almo_scf_delocalization

! *****************************************************************************
!> \brief create various matrices 
!> \param almo_scf_env ...
!> \param matrix_s0 ...
!> \param error ...
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
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_clean_up(almo_scf_env,error)

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
    DEALLOCATE(almo_scf_env%last_atom_of_domain)
    DEALLOCATE(almo_scf_env%nbasis_of_domain)
    DEALLOCATE(almo_scf_env%nocc_of_domain)
    DEALLOCATE(almo_scf_env%nvirt_full_of_domain)
    DEALLOCATE(almo_scf_env%nvirt_of_domain)
    DEALLOCATE(almo_scf_env%nvirt_disc_of_domain)
    DEALLOCATE(almo_scf_env%mu_of_domain)
    DEALLOCATE(almo_scf_env%cpu_of_domain)
    DEALLOCATE(almo_scf_env%charge_of_domain)
    DEALLOCATE(almo_scf_env%multiplicity_of_domain)
    
    DEALLOCATE(almo_scf_env%domain_index_of_ao_block)
    DEALLOCATE(almo_scf_env%domain_index_of_mo_block)

    CALL cp_para_env_release(almo_scf_env%para_env,error)
    CALL cp_blacs_env_release(almo_scf_env%blacs_env,error)

    CALL timestop(handle)

  END SUBROUTINE almo_scf_clean_up

END MODULE almo_scf

