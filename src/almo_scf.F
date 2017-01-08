!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for all ALMO-based SCF methods
!>        RZK-warning marks unresolved issues
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE almo_scf
   USE almo_scf_methods,                ONLY: almo_scf_p_blk_to_t_blk,&
                                              almo_scf_t_blk_to_p,&
                                              almo_scf_t_blk_to_t_blk_orthonormal,&
                                              distribute_domains,&
                                              orthogonalize_mos
   USE almo_scf_optimizer,              ONLY: almo_scf_block_diagonal,&
                                              almo_scf_xalmo_eigensolver,&
                                              almo_scf_xalmo_pcg
   USE almo_scf_qs,                     ONLY: almo_scf_construct_quencher,&
                                              almo_scf_init_qs,&
                                              calculate_w_matrix_almo,&
                                              matrix_almo_create,&
                                              matrix_almo_to_qs,&
                                              matrix_qs_to_almo
   USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                              almo_mat_dim_occ,&
                                              almo_mat_dim_virt,&
                                              almo_mat_dim_virt_disc,&
                                              almo_mat_dim_virt_full,&
                                              almo_scf_env_type,&
                                              optimizer_options_type,&
                                              print_optimizer_options
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Khaliullin2013,&
                                              Kolafa2004,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_retain
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_retain
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_binary_read, dbcsr_checksum, dbcsr_copy, dbcsr_create, &
        dbcsr_distribution, dbcsr_distribution_mp, dbcsr_filter, dbcsr_finalize, dbcsr_get_info, &
        dbcsr_get_stored_coordinates, dbcsr_init, dbcsr_mp_mynode, dbcsr_multiply, &
        dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_block2d, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric, dbcsr_work_create
   USE domain_submatrix_methods,        ONLY: init_submatrices,&
                                              release_submatrices
   USE input_constants,                 ONLY: &
        almo_deloc_none, almo_deloc_scf, almo_deloc_x, almo_deloc_x_then_scf, &
        almo_deloc_xalmo_1diag, almo_deloc_xalmo_scf, almo_deloc_xalmo_x, almo_deloc_xk, &
        almo_domain_layout_molecular, almo_mat_distr_atomic, almo_mat_distr_molecular, &
        almo_scf_diag, almo_scf_dm_sign, almo_scf_pcg, atomic_guess, molecular_guess, &
        optimizer_diis, optimizer_pcg, restart_guess, virt_full, virt_number, virt_occ_size, &
        xalmo_case_block_diag, xalmo_case_fully_deloc, xalmo_case_normal
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathlib,                         ONLY: binomial
   USE molecule_types,                  ONLY: get_molecule_set_info,&
                                              molecule_type
   USE mscfg_types,                     ONLY: get_matrix_from_submatrices,&
                                              molecular_scf_guess_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_initial_guess,                ONLY: calculate_atomic_block_dm,&
                                              calculate_mopac_dm
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              scf_env_release
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf'

   PUBLIC :: almo_entry_scf

   LOGICAL, PARAMETER :: debug_mode = .FALSE.
   LOGICAL, PARAMETER :: safe_mode = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief The entry point into ALMO SCF routines
!> \param qs_env   pointer to the QS environment
!> \param calc_forces   calculate forces?
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_entry_scf(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_entry_scf', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(almo_scf_env_type), POINTER                   :: almo_scf_env

      CALL timeset(routineN, handle)

      CALL cite_reference(Khaliullin2013)

      ! get a pointer to the almo environment
      CALL get_qs_env(qs_env, almo_scf_env=almo_scf_env)

      ! initialize scf
      CALL almo_scf_init(qs_env, almo_scf_env, calc_forces)

      ! create the initial guess for ALMOs
      CALL almo_scf_initial_guess(qs_env, almo_scf_env)

      ! perform SCF for block diagonal ALMOs
      CALL almo_scf_main(qs_env, almo_scf_env)

      ! allow electron delocalization
      CALL almo_scf_delocalization(qs_env, almo_scf_env)

      ! electron correlation methods
      !CALL almo_correlation_main(qs_env,almo_scf_env)

      ! do post scf processing
      CALL almo_scf_post(qs_env, almo_scf_env)

      ! clean up the mess
      CALL almo_scf_clean_up(almo_scf_env)

      CALL timestop(handle)

   END SUBROUTINE almo_entry_scf

! **************************************************************************************************
!> \brief Initialization of the almo_scf_env_type.
!> \param qs_env ...
!> \param almo_scf_env ...
!> \param calc_forces ...
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_init(qs_env, almo_scf_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env
      LOGICAL                                            :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: ao, handle, i, iao, idomain, ispin, &
                                                            multip, naos, natoms, ndomains, nelec, &
                                                            nelec_a, nelec_b, nmols, nspins, &
                                                            unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ! define the output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! set optimizers' types
      almo_scf_env%opt_block_diag_diis%optimizer_type = optimizer_diis
      almo_scf_env%opt_block_diag_pcg%optimizer_type = optimizer_pcg
      almo_scf_env%opt_xalmo_diis%optimizer_type = optimizer_diis
      almo_scf_env%opt_xalmo_pcg%optimizer_type = optimizer_pcg

      ! get info from the qs_env
      CALL get_qs_env(qs_env, &
                      nelectron_total=almo_scf_env%nelectrons_total, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control, &
                      molecule_set=molecule_set, &
                      input=input, &
                      has_unit_metric=almo_scf_env%orthogonal_basis, &
                      para_env=almo_scf_env%para_env, &
                      blacs_env=almo_scf_env%blacs_env, &
                      nelectron_spin=almo_scf_env%nelectrons_spin)
      CALL cp_para_env_retain(almo_scf_env%para_env)
      CALL cp_blacs_env_retain(almo_scf_env%blacs_env)

      ! copy basic quantities
      almo_scf_env%nspins = dft_control%nspins
      almo_scf_env%natoms = dbcsr_nblkrows_total(matrix_s(1)%matrix)
      almo_scf_env%nmolecules = SIZE(molecule_set)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=naos)
      almo_scf_env%naos = naos

      ! convenient local varibales
      nspins = almo_scf_env%nspins
      nmols = almo_scf_env%nmolecules
      natoms = almo_scf_env%natoms

      ! parse the almo_scf section and set appropriate quantities
      !CALL almo_scf_init_read_write_input(input, almo_scf_env)

      ! Define groups: either atomic or molecular
      IF (almo_scf_env%domain_layout_mos == almo_domain_layout_molecular) THEN
         almo_scf_env%ndomains = almo_scf_env%nmolecules
      ELSE
         almo_scf_env%ndomains = almo_scf_env%natoms
      ENDIF

      ! allocate domain descriptors
      ndomains = almo_scf_env%ndomains
      ALLOCATE (almo_scf_env%domain_index_of_atom(natoms))
      ALLOCATE (almo_scf_env%domain_index_of_ao(naos))
      ALLOCATE (almo_scf_env%first_atom_of_domain(ndomains))
      ALLOCATE (almo_scf_env%last_atom_of_domain(ndomains))
      ALLOCATE (almo_scf_env%nbasis_of_domain(ndomains))
      ALLOCATE (almo_scf_env%nocc_of_domain(ndomains, nspins))
      ALLOCATE (almo_scf_env%nvirt_full_of_domain(ndomains, nspins))
      ALLOCATE (almo_scf_env%nvirt_of_domain(ndomains, nspins))
      ALLOCATE (almo_scf_env%nvirt_disc_of_domain(ndomains, nspins))
      ALLOCATE (almo_scf_env%mu_of_domain(ndomains, nspins))
      ALLOCATE (almo_scf_env%cpu_of_domain(ndomains))
      ALLOCATE (almo_scf_env%charge_of_domain(ndomains))
      ALLOCATE (almo_scf_env%multiplicity_of_domain(ndomains))

      ! fill out domain descriptors and group descriptors
      IF (almo_scf_env%domain_layout_mos == almo_domain_layout_molecular) THEN
         ! get domain info from molecule_set
         CALL get_molecule_set_info(molecule_set, &
                                    atom_to_mol=almo_scf_env%domain_index_of_atom, &
                                    mol_to_first_atom=almo_scf_env%first_atom_of_domain, &
                                    mol_to_last_atom=almo_scf_env%last_atom_of_domain, &
                                    mol_to_nelectrons=almo_scf_env%nocc_of_domain(1:ndomains, 1), &
                                    mol_to_nbasis=almo_scf_env%nbasis_of_domain, &
                                    mol_to_charge=almo_scf_env%charge_of_domain, &
                                    mol_to_multiplicity=almo_scf_env%multiplicity_of_domain)
         ! calculate number of alpha and beta occupied orbitals from
         ! the number of electrons and multiplicity of each molecule
         ! Na + Nb = Ne
         ! Na - Nb = Mult - 1 (assume Na > Nb as we do not have more info from get_molecule_set_info)
         DO idomain = 1, ndomains
            nelec = almo_scf_env%nocc_of_domain(idomain, 1)
            multip = almo_scf_env%multiplicity_of_domain(idomain)
            nelec_a = (nelec+multip-1)/2
            nelec_b = nelec-nelec_a
            almo_scf_env%nocc_of_domain(idomain, 1) = nelec_a
            IF (nelec_a .NE. nelec_b) THEN
               IF (nspins .EQ. 1) THEN
                  WRITE (*, *) "Domain ", idomain, " out of ", ndomains, ". Electrons = ", nelec
                  CPABORT("odd e- -- use unrestricted methods")
               ENDIF
               almo_scf_env%nocc_of_domain(idomain, 2) = nelec_b
               ! RZK-warning: open-shell procedures have not been tested yet
               ! Stop the program now
               CPABORT("Unrestricted ALMO methods are NYI")
            ENDIF
         ENDDO
         DO ispin = 1, nspins
            ! take care of the full virtual subspace
            almo_scf_env%nvirt_full_of_domain(:, ispin) = &
               almo_scf_env%nbasis_of_domain(:)- &
               almo_scf_env%nocc_of_domain(:, ispin)
            ! and the truncated virtual subspace
            SELECT CASE (almo_scf_env%deloc_truncate_virt)
            CASE (virt_full)
               almo_scf_env%nvirt_of_domain(:, ispin) = &
                  almo_scf_env%nvirt_full_of_domain(:, ispin)
               almo_scf_env%nvirt_disc_of_domain(:, ispin) = 0
            CASE (virt_number)
               DO idomain = 1, ndomains
                  almo_scf_env%nvirt_of_domain(idomain, ispin) = &
                     MIN(almo_scf_env%deloc_virt_per_domain, &
                         almo_scf_env%nvirt_full_of_domain(idomain, ispin))
                  almo_scf_env%nvirt_disc_of_domain(idomain, ispin) = &
                     almo_scf_env%nvirt_full_of_domain(idomain, ispin)- &
                     almo_scf_env%nvirt_of_domain(idomain, ispin)
               ENDDO
            CASE (virt_occ_size)
               DO idomain = 1, ndomains
                  almo_scf_env%nvirt_of_domain(idomain, ispin) = &
                     MIN(almo_scf_env%nocc_of_domain(idomain, ispin), &
                         almo_scf_env%nvirt_full_of_domain(idomain, ispin))
                  almo_scf_env%nvirt_disc_of_domain(idomain, ispin) = &
                     almo_scf_env%nvirt_full_of_domain(idomain, ispin)- &
                     almo_scf_env%nvirt_of_domain(idomain, ispin)
               ENDDO
            CASE DEFAULT
               CPABORT("illegal method for virtual space truncation")
            END SELECT
         ENDDO ! spin
      ELSE ! domains are atomic
         ! RZK-warning do the same for atomic domains/groups
         almo_scf_env%domain_index_of_atom(1:natoms) = (/(i, i=1, natoms)/)
      ENDIF

      ao = 1
      DO idomain = 1, ndomains
         DO iao = 1, almo_scf_env%nbasis_of_domain(idomain)
            almo_scf_env%domain_index_of_ao(ao) = idomain
            ao = ao+1
         ENDDO
      ENDDO

      almo_scf_env%mu_of_domain(:, :) = almo_scf_env%mu

      ! build domain (i.e. layout) indices for distribution blocks
      ! ao blocks
      IF (almo_scf_env%mat_distr_aos == almo_mat_distr_atomic) THEN
         ALLOCATE (almo_scf_env%domain_index_of_ao_block(natoms))
         almo_scf_env%domain_index_of_ao_block(:) = &
            almo_scf_env%domain_index_of_atom(:)
      ELSE IF (almo_scf_env%mat_distr_aos == almo_mat_distr_molecular) THEN
         ALLOCATE (almo_scf_env%domain_index_of_ao_block(nmols))
         ! if distr blocks are molecular then domain layout is also molecular
         almo_scf_env%domain_index_of_ao_block(:) = (/(i, i=1, nmols)/)
      ENDIF
      ! mo blocks
      IF (almo_scf_env%mat_distr_mos == almo_mat_distr_atomic) THEN
         ALLOCATE (almo_scf_env%domain_index_of_mo_block(natoms))
         almo_scf_env%domain_index_of_mo_block(:) = &
            almo_scf_env%domain_index_of_atom(:)
      ELSE IF (almo_scf_env%mat_distr_mos == almo_mat_distr_molecular) THEN
         ALLOCATE (almo_scf_env%domain_index_of_mo_block(nmols))
         ! if distr blocks are molecular then domain layout is also molecular
         almo_scf_env%domain_index_of_mo_block(:) = (/(i, i=1, nmols)/)
      ENDIF

      ! set all flags
      !almo_scf_env%need_previous_ks=.FALSE.
      !IF (almo_scf_env%deloc_method==almo_deloc_harris) THEN
      almo_scf_env%need_previous_ks = .TRUE.
      !ENDIF

      !almo_scf_env%need_virtuals=.FALSE.
      !almo_scf_env%need_orbital_energies=.FALSE.
      !IF (almo_scf_env%almo_update_algorithm==almo_scf_diag) THEN
      almo_scf_env%need_virtuals = .TRUE.
      almo_scf_env%need_orbital_energies = .TRUE.
      !ENDIF

      almo_scf_env%calc_forces = calc_forces
      IF (calc_forces) THEN
         IF (almo_scf_env%deloc_method == almo_deloc_x .OR. &
             almo_scf_env%deloc_method == almo_deloc_xalmo_x .OR. &
             almo_scf_env%deloc_method == almo_deloc_xalmo_1diag) &
            CPABORT(" Forces are not implemented for this method yet ")
      ENDIF

      almo_scf_env%calc_forces = calc_forces
      IF (calc_forces) THEN
         IF (almo_scf_env%deloc_method == almo_deloc_x .OR. &
             almo_scf_env%deloc_method == almo_deloc_xalmo_x .OR. &
             almo_scf_env%deloc_method == almo_deloc_xalmo_1diag) &
            CPABORT(" Forces are not implemented for this method yet ")
      ENDIF

      ! create all matrices
      CALL almo_scf_env_create_matrices(almo_scf_env, matrix_s(1)%matrix)
      ! set up matrix S and all required functions of S
      almo_scf_env%s_inv_done = .FALSE.
      almo_scf_env%s_sqrt_done = .FALSE.
      CALL almo_scf_init_ao_overlap(matrix_s(1)%matrix, almo_scf_env)

      ! create the quencher (imposes sparsity template)
      CALL almo_scf_construct_quencher(qs_env, almo_scf_env)
      CALL distribute_domains(almo_scf_env)

      ! FINISH setting job parameters here, print out job info
      CALL almo_scf_print_job_info(almo_scf_env, unit_nr)

      ! allocate and init the domain preconditioner
      ALLOCATE (almo_scf_env%domain_preconditioner(ndomains, nspins))
      CALL init_submatrices(almo_scf_env%domain_preconditioner)

      ! allocate and init projected KS for domains
      ALLOCATE (almo_scf_env%domain_ks_xx(ndomains, nspins))
      CALL init_submatrices(almo_scf_env%domain_ks_xx)

      ! init ao-overlap subblocks
      ALLOCATE (almo_scf_env%domain_s_inv(ndomains, nspins))
      CALL init_submatrices(almo_scf_env%domain_s_inv)
      ALLOCATE (almo_scf_env%domain_s_sqrt_inv(ndomains, nspins))
      CALL init_submatrices(almo_scf_env%domain_s_sqrt_inv)
      ALLOCATE (almo_scf_env%domain_s_sqrt(ndomains, nspins))
      CALL init_submatrices(almo_scf_env%domain_s_sqrt)
      ALLOCATE (almo_scf_env%domain_t(ndomains, nspins))
      CALL init_submatrices(almo_scf_env%domain_t)
      ALLOCATE (almo_scf_env%domain_err(ndomains, nspins))
      CALL init_submatrices(almo_scf_env%domain_err)
      ALLOCATE (almo_scf_env%domain_r_down_up(ndomains, nspins))
      CALL init_submatrices(almo_scf_env%domain_r_down_up)

      ! initialization of the QS settings with the ALMO flavor
      CALL almo_scf_init_qs(qs_env, almo_scf_env)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_init

! **************************************************************************************************
!> \brief create the scf initial guess for ALMOs
!> \param qs_env ...
!> \param almo_scf_env ...
!> \par History
!>       2016.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_initial_guess(qs_env, almo_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_initial_guess', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: handle, iaspc, ispin, istore, naspc, &
                                                            nspins, unit_nr
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      LOGICAL                                            :: aspc_guess, has_unit_metric
      REAL(KIND=dp)                                      :: alpha, cs_pos
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(molecular_scf_guess_env_type), POINTER        :: mscfg_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      NULLIFY (rho, rho_ao)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! get basic quantities from the qs_env
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      has_unit_metric=has_unit_metric, &
                      para_env=para_env, &
                      nelectron_spin=nelectron_spin, &
                      mscfg_env=mscfg_env, &
                      rho=rho)

      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CPASSERT(ASSOCIATED(mscfg_env))

      ! initial guess on the first simulation step is determined by almo_scf_env%almo_scf_guess
      ! the subsequent simulation steps are determined by extrapolation_order
      ! if extrapolation order is zero then again almo_scf_env%almo_scf_guess is used
      ! ... the number of stored history points will remain zero if extrapolation order is zero
      IF (almo_scf_env%almo_history%istore == 0) THEN
         aspc_guess = .FALSE.
      ELSE
         aspc_guess = .TRUE.
      ENDIF

      nspins = almo_scf_env%nspins

      ! create an initial guess
      IF (.NOT. aspc_guess) THEN

         SELECT CASE (almo_scf_env%almo_scf_guess)
         CASE (molecular_guess)

            DO ispin = 1, nspins

               ! the calculations on "isolated" molecules has already been done
               ! all we need to do is convert the MOs of molecules into
               ! the ALMO matrix taking into account different distributions
               CALL get_matrix_from_submatrices(mscfg_env, &
                                                almo_scf_env%matrix_t_blk(ispin), ispin)
               CALL dbcsr_filter(almo_scf_env%matrix_t_blk(ispin), &
                                 almo_scf_env%eps_filter)

            ENDDO

         CASE (atomic_guess)

            DO ispin = 1, nspins
               CALL dbcsr_set(rho_ao(ispin)%matrix, 0.0_dp)
               IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical) THEN
                  CALL calculate_mopac_dm(rho_ao(ispin)%matrix, &
                                          matrix_s(1)%matrix, has_unit_metric, &
                                          dft_control, particle_set, atomic_kind_set, qs_kind_set, &
                                          nspins, nelectron_spin(ispin), &
                                          para_env)
               ELSE
                  CALL calculate_atomic_block_dm(rho_ao(ispin)%matrix, &
                                                 matrix_s(1)%matrix, particle_set, atomic_kind_set, qs_kind_set, &
                                                 ispin, nspins, nelectron_spin(ispin), unit_nr)
               ENDIF

               ! copy the atomic-block dm into matrix_p_blk
               CALL matrix_qs_to_almo(rho_ao(ispin)%matrix, &
                                      almo_scf_env%matrix_p_blk(ispin), almo_scf_env, &
                                      .FALSE.)
               CALL dbcsr_filter(almo_scf_env%matrix_p_blk(ispin), &
                                 almo_scf_env%eps_filter)

            ENDDO ! ispin

            ! obtain orbitals from the density matrix
            ! (the current version of ALMO SCF needs orbitals)
            CALL almo_scf_p_blk_to_t_blk(almo_scf_env, ionic=.FALSE.)

         CASE (restart_guess)

            project_name = logger%iter_info%project_name

            DO ispin = 1, nspins
               WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_ALMO_SPIN_", ispin, "_RESTART.mo"
               CALL dbcsr_binary_read(file_name, distribution=dbcsr_distribution(almo_scf_env%matrix_t_blk(ispin)), &
                                      matrix_new=almo_scf_env%matrix_t_blk(ispin))
               cs_pos = dbcsr_checksum(almo_scf_env%matrix_t_blk(ispin), pos=.TRUE.)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,A,E20.8)') "Read restart ALMO "//TRIM(file_name)//" with checksum: ", cs_pos
               ENDIF
            ENDDO
         END SELECT

      ELSE !aspc_guess

         CALL cite_reference(Kolafa2004)

         naspc = MIN(almo_scf_env%almo_history%istore, almo_scf_env%almo_history%nstore)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, FMT="(/,T2,A,/,/,T3,A,I0)") &
               "Parameters for the always stable predictor-corrector (ASPC) method:", &
               "ASPC order: ", naspc
         END IF

         DO ispin = 1, nspins

            ! extrapolation
            DO iaspc = 1, naspc
               istore = MOD(almo_scf_env%almo_history%istore-iaspc, almo_scf_env%almo_history%nstore)+1
               alpha = (-1.0_dp)**(iaspc+1)*REAL(iaspc, KIND=dp)* &
                       binomial(2*naspc, naspc-iaspc)/binomial(2*naspc-2, naspc-1)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, FMT="(T3,A2,I0,A4,F10.6)") &
                     "B(", iaspc, ") = ", alpha
               END IF
               IF (iaspc == 1) THEN
                  CALL dbcsr_copy(almo_scf_env%matrix_t_blk(ispin), &
                                  almo_scf_env%almo_history%matrix_t(ispin), &
                                  keep_sparsity=.TRUE.)
                  CALL dbcsr_scale(almo_scf_env%matrix_t_blk(ispin), alpha)
               ELSE
                  CALL dbcsr_multiply("N", "N", alpha, &
                                      almo_scf_env%almo_history%matrix_p_up_down(ispin, istore), &
                                      almo_scf_env%almo_history%matrix_t(ispin), &
                                      1.0_dp, almo_scf_env%matrix_t_blk(ispin), &
                                      retain_sparsity=.TRUE.)
               ENDIF
            ENDDO !iaspc

         ENDDO !ispin

      ENDIF !aspc_guess?

      CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env)
      CALL almo_scf_t_blk_to_p(almo_scf_env, &
                               use_sigma_inv_guess=.FALSE.)

      DO ispin = 1, nspins
         CALL matrix_almo_to_qs(almo_scf_env%matrix_p(ispin), &
                                rho_ao(ispin)%matrix, &
                                almo_scf_env)
      ENDDO

      CALL qs_rho_update_rho(rho, qs_env=qs_env)
      CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
      CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., &
                               just_energy=.FALSE.)

      IF (unit_nr > 0) THEN
         IF (almo_scf_env%almo_scf_guess .EQ. molecular_guess) THEN
            WRITE (unit_nr, '(T2,A38,F40.10)') "Single-molecule energy:", &
               SUM(mscfg_env%energy_of_frag)
         ENDIF
         WRITE (unit_nr, '(T2,A38,F40.10)') "Energy of the initial guess:", energy%total
         WRITE (unit_nr, '()')
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE almo_scf_initial_guess

! **************************************************************************************************
!> \brief store a history of matrices for later use in almo_scf_initial_guess
!> \param almo_scf_env ...
!> \par History
!>       2016.11 created [Rustam Z Khaliullin]
!> \author Rustam Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_store_result(almo_scf_env)
      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_store_result', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istore, unit_nr
      LOGICAL :: delocalization_uses_extrapolation
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_no_tmp1, matrix_no_tmp2, &
                                                            matrix_no_tmp3, matrix_no_tmp4

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (almo_scf_env%almo_history%nstore > 0) THEN

         almo_scf_env%almo_history%istore = almo_scf_env%almo_history%istore+1

         DO ispin = 1, SIZE(almo_scf_env%matrix_t_blk)

            istore = MOD(almo_scf_env%almo_history%istore-1, almo_scf_env%almo_history%nstore)+1

            IF (almo_scf_env%almo_history%istore == 1) THEN
               CALL dbcsr_init(almo_scf_env%almo_history%matrix_t(ispin))
               CALL dbcsr_create(almo_scf_env%almo_history%matrix_t(ispin), &
                                 template=almo_scf_env%matrix_t_blk(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
            ENDIF
            CALL dbcsr_copy(almo_scf_env%almo_history%matrix_t(ispin), &
                            almo_scf_env%matrix_t_blk(ispin))

            IF (almo_scf_env%almo_history%istore <= almo_scf_env%almo_history%nstore) THEN
               CALL dbcsr_init(almo_scf_env%almo_history%matrix_p_up_down(ispin, istore))
               CALL dbcsr_create(almo_scf_env%almo_history%matrix_p_up_down(ispin, istore), &
                                 template=almo_scf_env%matrix_s(1), &
                                 matrix_type=dbcsr_type_no_symmetry)
            ENDIF

            CALL dbcsr_init(matrix_no_tmp1)
            CALL dbcsr_create(matrix_no_tmp1, template=almo_scf_env%matrix_t_blk(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_init(matrix_no_tmp2)
            CALL dbcsr_create(matrix_no_tmp2, template=almo_scf_env%matrix_t_blk(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            ! compute contra-covariant density matrix
            CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1), &
                                almo_scf_env%matrix_t_blk(ispin), &
                                0.0_dp, matrix_no_tmp1, &
                                filter_eps=almo_scf_env%eps_filter)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_no_tmp1, &
                                almo_scf_env%matrix_sigma_inv_0deloc(ispin), &
                                0.0_dp, matrix_no_tmp2, &
                                filter_eps=almo_scf_env%eps_filter)
            CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                almo_scf_env%matrix_t_blk(ispin), &
                                matrix_no_tmp2, &
                                0.0_dp, almo_scf_env%almo_history%matrix_p_up_down(ispin, istore), &
                                filter_eps=almo_scf_env%eps_filter)

            CALL dbcsr_release(matrix_no_tmp1)
            CALL dbcsr_release(matrix_no_tmp2)

         ENDDO

      ENDIF

      ! exrapolate xalmos?
      delocalization_uses_extrapolation = &
         .NOT. ((almo_scf_env%deloc_method .EQ. almo_deloc_none) .OR. &
                (almo_scf_env%deloc_method .EQ. almo_deloc_xalmo_1diag))
      IF (almo_scf_env%xalmo_history%nstore > 0 .AND. &
          delocalization_uses_extrapolation) THEN

         almo_scf_env%xalmo_history%istore = almo_scf_env%xalmo_history%istore+1

         DO ispin = 1, SIZE(almo_scf_env%matrix_t)

            istore = MOD(almo_scf_env%xalmo_history%istore-1, almo_scf_env%xalmo_history%nstore)+1

            IF (almo_scf_env%xalmo_history%istore == 1) THEN
               CALL dbcsr_init(almo_scf_env%xalmo_history%matrix_t(ispin))
               CALL dbcsr_create(almo_scf_env%xalmo_history%matrix_t(ispin), &
                                 template=almo_scf_env%matrix_t(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
            ENDIF
            CALL dbcsr_copy(almo_scf_env%xalmo_history%matrix_t(ispin), &
                            almo_scf_env%matrix_t(ispin))

            IF (almo_scf_env%xalmo_history%istore <= almo_scf_env%xalmo_history%nstore) THEN
               !CALL dbcsr_init(almo_scf_env%xalmo_history%matrix_x(ispin, istore))
               !CALL dbcsr_create(almo_scf_env%xalmo_history%matrix_x(ispin, istore), &
               !        template=almo_scf_env%matrix_t(ispin), &
               !        matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_init(almo_scf_env%xalmo_history%matrix_p_up_down(ispin, istore))
               CALL dbcsr_create(almo_scf_env%xalmo_history%matrix_p_up_down(ispin, istore), &
                                 template=almo_scf_env%matrix_s(1), &
                                 matrix_type=dbcsr_type_no_symmetry)
            ENDIF

            CALL dbcsr_init(matrix_no_tmp3)
            CALL dbcsr_create(matrix_no_tmp3, template=almo_scf_env%matrix_t(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_init(matrix_no_tmp4)
            CALL dbcsr_create(matrix_no_tmp4, template=almo_scf_env%matrix_t(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            ! compute contra-covariant density matrix
            CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1), &
                                almo_scf_env%matrix_t(ispin), &
                                0.0_dp, matrix_no_tmp3, &
                                filter_eps=almo_scf_env%eps_filter)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_no_tmp3, &
                                almo_scf_env%matrix_sigma_inv(ispin), &
                                0.0_dp, matrix_no_tmp4, &
                                filter_eps=almo_scf_env%eps_filter)
            CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                almo_scf_env%matrix_t(ispin), &
                                matrix_no_tmp4, &
                                0.0_dp, almo_scf_env%xalmo_history%matrix_p_up_down(ispin, istore), &
                                filter_eps=almo_scf_env%eps_filter)

            ! store the difference between t and t0
            !CALL dbcsr_copy(almo_scf_env%xalmo_history%matrix_x(ispin, istore),&
            !        almo_scf_env%matrix_t(ispin))
            !CALL dbcsr_add(almo_scf_env%xalmo_history%matrix_x(ispin, istore),&
            !        almo_scf_env%matrix_t_blk(ispin),1.0_dp,-1.0_dp)

            CALL dbcsr_release(matrix_no_tmp3)
            CALL dbcsr_release(matrix_no_tmp4)

         ENDDO

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE almo_scf_store_result

! **************************************************************************************************
!> \brief Prints out a short summary about the ALMO SCF job
!> \param almo_scf_env ...
!> \param unit_nr ...
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_print_job_info(almo_scf_env, unit_nr)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_print_job_info', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=13)                                  :: neig_string
      CHARACTER(len=33)                                  :: deloc_method_string
      INTEGER                                            :: handle, idomain, index1_prev

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A,A,A)') REPEAT("-", 32), " ALMO SETTINGS ", REPEAT("-", 32)

         WRITE (unit_nr, '(T2,A,T48,E33.3)') "eps_filter:", almo_scf_env%eps_filter

         WRITE (unit_nr, '(T2,A)') "optimization of block-diagonal ALMOs:"
         SELECT CASE (almo_scf_env%almo_update_algorithm)
         CASE (almo_scf_diag)
            ! the DIIS algorith is the only choice for the diagonlaization-based algorithm
            CALL print_optimizer_options(almo_scf_env%opt_block_diag_diis, unit_nr)
         CASE (almo_scf_pcg)
            ! print out PCG options
            CALL print_optimizer_options(almo_scf_env%opt_block_diag_pcg, unit_nr)
         END SELECT

         SELECT CASE (almo_scf_env%deloc_method)
         CASE (almo_deloc_none)
            deloc_method_string = "NONE"
         CASE (almo_deloc_x)
            deloc_method_string = "FULL_X"
         CASE (almo_deloc_scf)
            deloc_method_string = "FULL_SCF"
         CASE (almo_deloc_x_then_scf)
            deloc_method_string = "FULL_X_THEN_SCF"
         CASE (almo_deloc_xalmo_1diag)
            deloc_method_string = "XALMO_1DIAG"
         CASE (almo_deloc_xalmo_x)
            deloc_method_string = "XALMO_X"
         CASE (almo_deloc_xalmo_scf)
            deloc_method_string = "XALMO_SCF"
         END SELECT
         WRITE (unit_nr, '(T2,A,T48,A33)') "delocalization:", TRIM(deloc_method_string)

         IF (almo_scf_env%deloc_method .NE. almo_deloc_none) THEN

            SELECT CASE (almo_scf_env%deloc_method)
            CASE (almo_deloc_x, almo_deloc_scf, almo_deloc_x_then_scf)
               WRITE (unit_nr, '(T2,A,T48,A33)') "delocalization cutoff radius:", &
                  "infinite"
               deloc_method_string = "FULL_X_THEN_SCF"
            CASE (almo_deloc_xalmo_1diag, almo_deloc_xalmo_x, almo_deloc_xalmo_scf)
               WRITE (unit_nr, '(T2,A,T48,F33.5)') "XALMO cutoff radius:", &
                  almo_scf_env%quencher_r0_factor
            END SELECT

            IF (almo_scf_env%deloc_method .EQ. almo_deloc_xalmo_1diag) THEN
               ! print nothing because no actual optimization is done
            ELSE
               WRITE (unit_nr, '(T2,A)') "optimization of extended orbitals:"
               SELECT CASE (almo_scf_env%xalmo_update_algorithm)
               CASE (almo_scf_diag)
                  CALL print_optimizer_options(almo_scf_env%opt_xalmo_diis, unit_nr)
               CASE (almo_scf_pcg)
                  CALL print_optimizer_options(almo_scf_env%opt_xalmo_pcg, unit_nr)
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

         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

         IF (almo_scf_env%ndomains .LE. 200) THEN

            ! print fragment info
            WRITE (unit_nr, '(T2,A13,A13,A13,A13,A13,A13)') &
               "Fragment", "Basis Set", "Occupied", "Virtual", "Charge", "Deloc Neig" !,"Discarded Virt"
            WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
            DO idomain = 1, almo_scf_env%ndomains

               IF (idomain .EQ. 1) THEN
                  index1_prev = 1
               ELSE
                  index1_prev = almo_scf_env%domain_map(1)%index1(idomain-1)
               ENDIF

               SELECT CASE (almo_scf_env%deloc_method)
               CASE (almo_deloc_none)
                  neig_string = "NONE"
               CASE (almo_deloc_x, almo_deloc_scf, almo_deloc_x_then_scf)
                  neig_string = "ALL"
               CASE (almo_deloc_xalmo_1diag, almo_deloc_xalmo_x, almo_deloc_xalmo_scf)
                  WRITE (neig_string, '(I13)') &
                     almo_scf_env%domain_map(1)%index1(idomain)-index1_prev-1 ! minus self
               CASE DEFAULT
                  neig_string = "N/A"
               END SELECT

               WRITE (unit_nr, '(T2,I13,I13,I13,I13,I13,A13)') &
                  idomain, almo_scf_env%nbasis_of_domain(idomain), &
                  SUM(almo_scf_env%nocc_of_domain(idomain, :)), &
                  SUM(almo_scf_env%nvirt_of_domain(idomain, :)), &
                  !SUM(almo_scf_env%nvirt_disc_of_domain(idomain,:)),&
                  almo_scf_env%charge_of_domain(idomain), &
                  ADJUSTR(TRIM(neig_string))

            ENDDO ! cycle over domains

            SELECT CASE (almo_scf_env%deloc_method)
            CASE (almo_deloc_xalmo_1diag, almo_deloc_xalmo_x, almo_deloc_xalmo_scf)

               WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

               ! print fragment neigbors
               WRITE (unit_nr, '(T2,A78)') &
                  "Neighbor lists (including self)"
               WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
               DO idomain = 1, almo_scf_env%ndomains

                  IF (idomain .EQ. 1) THEN
                     index1_prev = 1
                  ELSE
                     index1_prev = almo_scf_env%domain_map(1)%index1(idomain-1)
                  ENDIF

                  WRITE (unit_nr, '(T2,I13,":")') idomain
                  WRITE (unit_nr, '(T15,5I13)') &
                     almo_scf_env%domain_map(1)%pairs &
                     (index1_prev:almo_scf_env%domain_map(1)%index1(idomain)-1, 1) ! includes self

               ENDDO ! cycle over domains

            END SELECT

            WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

         ENDIF

         WRITE (unit_nr, '()')

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE almo_scf_print_job_info

! **************************************************************************************************
!> \brief Initializes the ALMO SCF copy of the AO overlap matrix
!>        and all necessary functions (sqrt, inverse...)
!> \param matrix_s ...
!> \param almo_scf_env ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_init_ao_overlap(matrix_s, almo_scf_env)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s
      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_init_ao_overlap', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! make almo copy of S
      ! also copy S to S_blk (i.e. to S with the domain structure imposed)
      IF (almo_scf_env%orthogonal_basis) THEN
         CALL dbcsr_set(almo_scf_env%matrix_s(1), 0.0_dp)
         CALL dbcsr_add_on_diag(almo_scf_env%matrix_s(1), 1.0_dp)
         CALL dbcsr_set(almo_scf_env%matrix_s_blk(1), 0.0_dp)
         CALL dbcsr_add_on_diag(almo_scf_env%matrix_s_blk(1), 1.0_dp)
      ELSE
         CALL matrix_qs_to_almo(matrix_s, almo_scf_env%matrix_s(1), &
                                almo_scf_env, .FALSE.)
         CALL matrix_qs_to_almo(matrix_s, almo_scf_env%matrix_s_blk(1), &
                                almo_scf_env, .TRUE.)
      ENDIF

      CALL dbcsr_filter(almo_scf_env%matrix_s(1), almo_scf_env%eps_filter)
      CALL dbcsr_filter(almo_scf_env%matrix_s_blk(1), almo_scf_env%eps_filter)

      IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_diag) THEN
         CALL matrix_sqrt_Newton_Schulz(almo_scf_env%matrix_s_blk_sqrt(1), &
                                        almo_scf_env%matrix_s_blk_sqrt_inv(1), &
                                        almo_scf_env%matrix_s_blk(1), &
                                        threshold=almo_scf_env%eps_filter, &
                                        order=almo_scf_env%order_lanczos, &
                                        eps_lanczos=almo_scf_env%eps_lanczos, &
                                        max_iter_lanczos=almo_scf_env%max_iter_lanczos)
      ELSE IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_dm_sign) THEN
         CALL invert_Hotelling(almo_scf_env%matrix_s_blk_inv(1), &
                               almo_scf_env%matrix_s_blk(1), &
                               threshold=almo_scf_env%eps_filter)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE almo_scf_init_ao_overlap

! **************************************************************************************************
!> \brief Selects the subroutine for the optimization of block-daigonal ALMOs.
!>        Keep it short and clean.
!> \param qs_env ...
!> \param almo_scf_env ...
!> \par History
!>       2011.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_main(qs_env, almo_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_main', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_pcg) THEN

         ! ALMO PCG optimizer as a special case of XALMO PCG
         CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
                                 almo_scf_env=almo_scf_env, &
                                 optimizer=almo_scf_env%opt_block_diag_pcg, &
                                 quench_t=almo_scf_env%quench_t_blk, &
                                 matrix_t_in=almo_scf_env%matrix_t_blk, &
                                 matrix_t_out=almo_scf_env%matrix_t_blk, &
                                 assume_t0_q0x=.FALSE., &
                                 perturbation_only=.FALSE., &
                                 special_case=xalmo_case_block_diag)

         CALL almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env)

      ELSE

         ! mixing/DIIS optimizer
         CALL almo_scf_block_diagonal(qs_env, almo_scf_env, &
                                      almo_scf_env%opt_block_diag_diis)

      ENDIF

      ! we might need a copy of the converged KS and sigma_inv
      DO ispin = 1, almo_scf_env%nspins
         CALL dbcsr_copy(almo_scf_env%matrix_ks_0deloc(ispin), &
                         almo_scf_env%matrix_ks(ispin))
         CALL dbcsr_copy(almo_scf_env%matrix_sigma_inv_0deloc(ispin), &
                         almo_scf_env%matrix_sigma_inv(ispin))
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE almo_scf_main

! **************************************************************************************************
!> \brief selects various post scf routines
!> \param qs_env ...
!> \param almo_scf_env ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_delocalization(qs_env, almo_scf_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_delocalization', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, handle, hold, iblock_col, &
                                                            iblock_row, ispin, mynode, &
                                                            nblkcols_tot, nblkrows_tot, row, &
                                                            unit_nr
      LOGICAL                                            :: almo_experimental, tr
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_new_block
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: no_quench
      TYPE(optimizer_options_type)                       :: arbitrary_optimizer

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! create a local optimizer that handles XALMO DIIS
      ! the options of this optimizer are arbitrary because
      ! XALMO DIIS SCF does not converge for yet unknown reasons and
      ! currently used in the code to get perturbative estimates only
      arbitrary_optimizer%optimizer_type = optimizer_diis
      arbitrary_optimizer%max_iter = 3
      arbitrary_optimizer%eps_error = 1.0E-6_dp
      arbitrary_optimizer%ndiis = 2

      SELECT CASE (almo_scf_env%deloc_method)
      CASE (almo_deloc_x, almo_deloc_scf, almo_deloc_x_then_scf)

         ! RZK-warning hack into the quenched routine:
         ! create a quench matrix with all-all-all blocks 1.0
         ! it is a waste of memory but since matrices are distributed
         ! we can tolerate it for now
         ALLOCATE (no_quench(almo_scf_env%nspins))
         CALL dbcsr_init(no_quench(1))
         CALL dbcsr_create(no_quench(1), &
                           template=almo_scf_env%matrix_t(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         mynode = dbcsr_mp_mynode(dbcsr_distribution_mp( &
                                  dbcsr_distribution(no_quench(1))))
         CALL dbcsr_work_create(no_quench(1), &
                                work_mutable=.TRUE.)
         nblkrows_tot = dbcsr_nblkrows_total(no_quench(1))
         nblkcols_tot = dbcsr_nblkcols_total(no_quench(1))
         ! RZK-warning: is it a quadratic-scaling routine?
         ! As a matter of fact it is! But this block treats
         ! fully delocalized MOs. So it is unavoidable.
         ! C'est la vie
         DO row = 1, nblkrows_tot
            DO col = 1, nblkcols_tot
               tr = .FALSE.
               iblock_row = row
               iblock_col = col
               CALL dbcsr_get_stored_coordinates(no_quench(1), &
                                                 iblock_row, iblock_col, hold)
               IF (hold .EQ. mynode) THEN
                  NULLIFY (p_new_block)
                  CALL dbcsr_reserve_block2d(no_quench(1), &
                                             iblock_row, iblock_col, p_new_block)
                  CPASSERT(ASSOCIATED(p_new_block))
                  p_new_block(:, :) = 1.0_dp
               ENDIF
            ENDDO
         ENDDO
         CALL dbcsr_finalize(no_quench(1))
         IF (almo_scf_env%nspins .GT. 1) THEN
            DO ispin = 2, almo_scf_env%nspins
               CALL dbcsr_init(no_quench(ispin))
               CALL dbcsr_create(no_quench(ispin), &
                                 template=almo_scf_env%matrix_t(1), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_copy(no_quench(ispin), no_quench(1))
            ENDDO
         ENDIF

      END SELECT

      SELECT CASE (almo_scf_env%deloc_method)
      CASE (almo_deloc_none, almo_deloc_scf)

         DO ispin = 1, almo_scf_env%nspins
            CALL dbcsr_copy(almo_scf_env%matrix_t(ispin), &
                            almo_scf_env%matrix_t_blk(ispin))
         ENDDO

      CASE (almo_deloc_x, almo_deloc_xk, almo_deloc_x_then_scf)

         !!!! RZK-warning a whole class of delocalization methods
         !!!! are commented out at the moment because some of their
         !!!! routines have not been thoroughly tested.

         !!!! if we have virtuals pre-optimize and truncate them
         !!!IF (almo_scf_env%need_virtuals) THEN
         !!!   SELECT CASE (almo_scf_env%deloc_truncate_virt)
         !!!   CASE (virt_full)
         !!!      ! simply copy virtual orbitals from matrix_v_full_blk to matrix_v_blk
         !!!      DO ispin=1,almo_scf_env%nspins
         !!!         CALL dbcsr_copy(almo_scf_env%matrix_v_blk(ispin),&
         !!!                 almo_scf_env%matrix_v_full_blk(ispin))
         !!!      ENDDO
         !!!   CASE (virt_number,virt_occ_size)
         !!!      CALL split_v_blk(almo_scf_env)
         !!!      !CALL truncate_subspace_v_blk(qs_env,almo_scf_env)
         !!!   CASE DEFAULT
         !!!      CPErrorMessage(cp_failure_level,routineP,"illegal method for virtual space truncation")
         !!!      CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
         !!!   END SELECT
         !!!ENDIF
         !!!CALL harris_foulkes_correction(qs_env,almo_scf_env)

         CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
                                 almo_scf_env=almo_scf_env, &
                                 optimizer=almo_scf_env%opt_xalmo_pcg, &
                                 quench_t=no_quench, &
                                 matrix_t_in=almo_scf_env%matrix_t_blk, &
                                 matrix_t_out=almo_scf_env%matrix_t, &
                                 assume_t0_q0x=.TRUE., &
                                 perturbation_only=.TRUE., &
                                 special_case=xalmo_case_fully_deloc)

      CASE (almo_deloc_xalmo_1diag)

         IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_diag) THEN

            almo_scf_env%perturbative_delocalization = .TRUE.
            DO ispin = 1, almo_scf_env%nspins
               CALL dbcsr_copy(almo_scf_env%matrix_t(ispin), &
                               almo_scf_env%matrix_t_blk(ispin))
            ENDDO
            CALL almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, &
                                            arbitrary_optimizer)

         ELSE

            CPABORT("Other algorithms do not exist")

         ENDIF

      CASE (almo_deloc_xalmo_x)

         IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_pcg) THEN

            CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
                                    almo_scf_env=almo_scf_env, &
                                    optimizer=almo_scf_env%opt_xalmo_pcg, &
                                    quench_t=almo_scf_env%quench_t, &
                                    matrix_t_in=almo_scf_env%matrix_t_blk, &
                                    matrix_t_out=almo_scf_env%matrix_t, &
                                    assume_t0_q0x=.TRUE., &
                                    perturbation_only=.TRUE., &
                                    special_case=xalmo_case_normal)

         ELSE

            CPABORT("Other algorithms do not exist")

         ENDIF

      CASE (almo_deloc_xalmo_scf)

         IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_diag) THEN

            CPABORT("Should not be here: convergence will fail!")

            almo_scf_env%perturbative_delocalization = .FALSE.
            DO ispin = 1, almo_scf_env%nspins
               CALL dbcsr_copy(almo_scf_env%matrix_t(ispin), &
                               almo_scf_env%matrix_t_blk(ispin))
            ENDDO
            CALL almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, &
                                            arbitrary_optimizer)

         ELSE IF (almo_scf_env%xalmo_update_algorithm .EQ. almo_scf_pcg) THEN

            CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
                                    almo_scf_env=almo_scf_env, &
                                    optimizer=almo_scf_env%opt_xalmo_pcg, &
                                    quench_t=almo_scf_env%quench_t, &
                                    matrix_t_in=almo_scf_env%matrix_t_blk, &
                                    matrix_t_out=almo_scf_env%matrix_t, &
                                    assume_t0_q0x=.TRUE., &
                                    perturbation_only=.FALSE., &
                                    special_case=xalmo_case_normal)

            ! RZK-warning THIS IS A HACK TO GET ORBITAL ENERGIES
            almo_experimental = .FALSE.
            IF (almo_experimental) THEN
               almo_scf_env%perturbative_delocalization = .TRUE.
               !DO ispin=1,almo_scf_env%nspins
               !   CALL dbcsr_copy(almo_scf_env%matrix_t(ispin),&
               !           almo_scf_env%matrix_t_blk(ispin))
               !ENDDO
               CALL almo_scf_xalmo_eigensolver(qs_env, almo_scf_env, &
                                               arbitrary_optimizer)
            ENDIF ! experimental

         ENDIF

      CASE DEFAULT

         CPABORT("Illegal delocalization method")

      END SELECT

      SELECT CASE (almo_scf_env%deloc_method)
      CASE (almo_deloc_scf, almo_deloc_x_then_scf)

         IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN
            CPABORT("full scf is NYI for truncated virtual space")
         ENDIF

         CALL almo_scf_xalmo_pcg(qs_env=qs_env, &
                                 almo_scf_env=almo_scf_env, &
                                 optimizer=almo_scf_env%opt_xalmo_pcg, &
                                 quench_t=no_quench, &
                                 matrix_t_in=almo_scf_env%matrix_t, &
                                 matrix_t_out=almo_scf_env%matrix_t, &
                                 assume_t0_q0x=.FALSE., &
                                 perturbation_only=.FALSE., &
                                 special_case=xalmo_case_fully_deloc)

      END SELECT

      ! clean up
      SELECT CASE (almo_scf_env%deloc_method)
      CASE (almo_deloc_x, almo_deloc_scf, almo_deloc_x_then_scf)
         DO ispin = 1, almo_scf_env%nspins
            CALL dbcsr_release(no_quench(ispin))
         ENDDO
         DEALLOCATE (no_quench)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE almo_scf_delocalization

! *****************************************************************************
!> \brief after SCF we have the final density and KS matrices compute various
!>        post-scf quantities
!> \param qs_env ...
!> \param almo_scf_env ...
!> \par History
!>       2015.03 created  [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_post(qs_env, almo_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_post', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_w
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: matrix_t_orthogonal
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CALL timeset(routineN, handle)

      ! store the matrices for a next scf run
      CALL almo_scf_store_result(almo_scf_env)

      ! orthogonalize orbitals before returning them to QS
      ALLOCATE (matrix_t_orthogonal(almo_scf_env%nspins))

      DO ispin = 1, almo_scf_env%nspins

         CALL dbcsr_init(matrix_t_orthogonal(ispin))
         CALL dbcsr_create(matrix_t_orthogonal(ispin), &
                           template=almo_scf_env%matrix_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_copy(matrix_t_orthogonal(ispin), &
                         almo_scf_env%matrix_t(ispin))

         CALL orthogonalize_mos(ket=matrix_t_orthogonal(ispin), &
                                overlap=almo_scf_env%matrix_sigma(ispin), &
                                metric=almo_scf_env%matrix_s(1), &
                                retain_locality=.FALSE., &
                                only_normalize=.FALSE., &
                                eps_filter=almo_scf_env%eps_filter, &
                                order_lanczos=almo_scf_env%order_lanczos, &
                                eps_lanczos=almo_scf_env%eps_lanczos, &
                                max_iter_lanczos=almo_scf_env%max_iter_lanczos)
      ENDDO

      ! return orthogonalized orbitals to QS
      NULLIFY (mos, mo_coeff, scf_env)

      CALL get_qs_env(qs_env, mos=mos, scf_env=scf_env)

      DO ispin = 1, almo_scf_env%nspins
         ! Currently only fm version of mo_set is usable.
         ! First transform the matrix_t to fm version
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)

         CALL copy_dbcsr_to_fm(matrix_t_orthogonal(ispin), mo_coeff)

         CALL dbcsr_release(matrix_t_orthogonal(ispin))
      ENDDO
      DEALLOCATE (matrix_t_orthogonal)

      ! calculate post scf properties
      ! CALL almo_post_scf_compute_properties(qs_env, almo_scf_env)
      CALL almo_post_scf_compute_properties(qs_env)

      ! In normal scf procedure this is done by scf,
      ! but it's not called by ALMO
      ! best to release it somewhere else
      IF (ASSOCIATED(scf_env)) THEN
         CALL scf_env_release(scf_env)
      ENDIF

      ! compute the W matrix if associated
      IF (almo_scf_env%calc_forces) THEN
         CALL get_qs_env(qs_env, matrix_w=matrix_w)
         IF (ASSOCIATED(matrix_w)) THEN
            CALL calculate_w_matrix_almo(matrix_w, almo_scf_env)
         ELSE
            CPABORT("Matrix W is needed but not associated")
         ENDIF
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE almo_scf_post

! **************************************************************************************************
!> \brief create various matrices
!> \param almo_scf_env ...
!> \param matrix_s0 ...
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_env_create_matrices(almo_scf_env, matrix_s0)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s0

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_env_create_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins

      CALL timeset(routineN, handle)

      nspins = almo_scf_env%nspins

      ! AO overlap matrix and its various functions
      CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s(1), &
                              matrix_qs=matrix_s0, &
                              almo_scf_env=almo_scf_env, &
                              name_new="S", &
                              size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
                              symmetry_new=dbcsr_type_symmetric, &
                              spin_key=0, &
                              init_domains=.FALSE.)
      CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s_blk(1), &
                              matrix_qs=matrix_s0, &
                              almo_scf_env=almo_scf_env, &
                              name_new="S_BLK", &
                              size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
                              symmetry_new=dbcsr_type_symmetric, &
                              spin_key=0, &
                              init_domains=.TRUE.)
      IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_diag) THEN
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s_blk_sqrt_inv(1), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="S_BLK_SQRT_INV", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
                                 symmetry_new=dbcsr_type_symmetric, &
                                 spin_key=0, &
                                 init_domains=.TRUE.)
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s_blk_sqrt(1), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="S_BLK_SQRT", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
                                 symmetry_new=dbcsr_type_symmetric, &
                                 spin_key=0, &
                                 init_domains=.TRUE.)
      ELSE IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_dm_sign) THEN
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_s_blk_inv(1), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="S_BLK_INV", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
                                 symmetry_new=dbcsr_type_symmetric, &
                                 spin_key=0, &
                                 init_domains=.TRUE.)
      ENDIF

      ! MO coeff matrices and their derivatives
      ALLOCATE (almo_scf_env%matrix_t_blk(nspins))
      ALLOCATE (almo_scf_env%quench_t_blk(nspins))
      ALLOCATE (almo_scf_env%matrix_err_blk(nspins))
      ALLOCATE (almo_scf_env%matrix_err_xx(nspins))
      ALLOCATE (almo_scf_env%matrix_sigma(nspins))
      ALLOCATE (almo_scf_env%matrix_sigma_inv(nspins))
      ALLOCATE (almo_scf_env%matrix_sigma_sqrt(nspins))
      ALLOCATE (almo_scf_env%matrix_sigma_sqrt_inv(nspins))
      ALLOCATE (almo_scf_env%matrix_sigma_blk(nspins))
      ALLOCATE (almo_scf_env%matrix_sigma_inv_0deloc(nspins))
      ALLOCATE (almo_scf_env%matrix_t(nspins))
      ALLOCATE (almo_scf_env%matrix_t_tr(nspins))
      DO ispin = 1, nspins
         ! create the blocked quencher
         CALL matrix_almo_create(matrix_new=almo_scf_env%quench_t_blk(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="Q_BLK", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                                 symmetry_new=dbcsr_type_no_symmetry, &
                                 spin_key=ispin, &
                                 init_domains=.TRUE.)
         ! create ALMO coefficient matrix
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_t_blk(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="T_BLK", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                                 symmetry_new=dbcsr_type_no_symmetry, &
                                 spin_key=ispin, &
                                 init_domains=.TRUE.)
         ! create the error matrix
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_err_blk(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="ERR_BLK", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
                                 symmetry_new=dbcsr_type_no_symmetry, &
                                 spin_key=ispin, &
                                 init_domains=.TRUE.)
         ! create the error matrix for the quenched ALMOs
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_err_xx(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="ERR_XX", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                                 symmetry_new=dbcsr_type_no_symmetry, &
                                 spin_key=ispin, &
                                 init_domains=.FALSE.)
         ! create a matrix with dimensions of a transposed mo coefficient matrix
         ! it might be necessary to perform the correction step using cayley
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_t_tr(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="T_TR", &
                                 size_keys=(/almo_mat_dim_occ, almo_mat_dim_aobasis/), &
                                 symmetry_new=dbcsr_type_no_symmetry, &
                                 spin_key=ispin, &
                                 init_domains=.FALSE.)
         ! create mo overlap matrix
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="SIG", &
                                 size_keys=(/almo_mat_dim_occ, almo_mat_dim_occ/), &
                                 symmetry_new=dbcsr_type_symmetric, &
                                 spin_key=ispin, &
                                 init_domains=.FALSE.)
         ! create blocked mo overlap matrix
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma_blk(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="SIG_BLK", &
                                 size_keys=(/almo_mat_dim_occ, almo_mat_dim_occ/), &
                                 symmetry_new=dbcsr_type_symmetric, &
                                 spin_key=ispin, &
                                 init_domains=.TRUE.)
         ! create blocked inverse mo overlap matrix
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma_inv_0deloc(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="SIGINV_BLK", &
                                 size_keys=(/almo_mat_dim_occ, almo_mat_dim_occ/), &
                                 symmetry_new=dbcsr_type_symmetric, &
                                 spin_key=ispin, &
                                 init_domains=.TRUE.)
         ! create inverse mo overlap matrix
         CALL matrix_almo_create( &
            matrix_new=almo_scf_env%matrix_sigma_inv(ispin), &
            matrix_qs=matrix_s0, &
            almo_scf_env=almo_scf_env, &
            name_new="SIGINV", &
            size_keys=(/almo_mat_dim_occ, almo_mat_dim_occ/), &
            symmetry_new=dbcsr_type_symmetric, &
            spin_key=ispin, &
            init_domains=.FALSE.)
         ! create various templates that will be necessary later
         CALL matrix_almo_create( &
            matrix_new=almo_scf_env%matrix_t(ispin), &
            matrix_qs=matrix_s0, &
            almo_scf_env=almo_scf_env, &
            name_new="T", &
            size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
            symmetry_new=dbcsr_type_no_symmetry, &
            spin_key=ispin, &
            init_domains=.FALSE.)
         CALL dbcsr_init(almo_scf_env%matrix_sigma_sqrt(ispin))
         CALL dbcsr_init(almo_scf_env%matrix_sigma_sqrt_inv(ispin))
         CALL dbcsr_create(almo_scf_env%matrix_sigma_sqrt(ispin), &
                           template=almo_scf_env%matrix_sigma(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(almo_scf_env%matrix_sigma_sqrt_inv(ispin), &
                           template=almo_scf_env%matrix_sigma(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
      ENDDO

      ! create virtual orbitals if necessary
      IF (almo_scf_env%need_virtuals) THEN
         ALLOCATE (almo_scf_env%matrix_v_blk(nspins))
         ALLOCATE (almo_scf_env%matrix_v_full_blk(nspins))
         ALLOCATE (almo_scf_env%matrix_v(nspins))
         ALLOCATE (almo_scf_env%matrix_vo(nspins))
         ALLOCATE (almo_scf_env%matrix_x(nspins))
         ALLOCATE (almo_scf_env%matrix_ov(nspins))
         ALLOCATE (almo_scf_env%matrix_ov_full(nspins))
         ALLOCATE (almo_scf_env%matrix_sigma_vv(nspins))
         ALLOCATE (almo_scf_env%matrix_sigma_vv_blk(nspins))
         ALLOCATE (almo_scf_env%matrix_sigma_vv_sqrt(nspins))
         ALLOCATE (almo_scf_env%matrix_sigma_vv_sqrt_inv(nspins))
         ALLOCATE (almo_scf_env%matrix_vv_full_blk(nspins))

         IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN
            ALLOCATE (almo_scf_env%matrix_k_blk(nspins))
            ALLOCATE (almo_scf_env%matrix_k_blk_ones(nspins))
            ALLOCATE (almo_scf_env%matrix_k_tr(nspins))
            ALLOCATE (almo_scf_env%matrix_v_disc(nspins))
            ALLOCATE (almo_scf_env%matrix_v_disc_blk(nspins))
            ALLOCATE (almo_scf_env%matrix_ov_disc(nspins))
            ALLOCATE (almo_scf_env%matrix_vv_disc_blk(nspins))
            ALLOCATE (almo_scf_env%matrix_vv_disc(nspins))
            ALLOCATE (almo_scf_env%opt_k_t_dd(nspins))
            ALLOCATE (almo_scf_env%opt_k_t_rr(nspins))
            ALLOCATE (almo_scf_env%opt_k_denom(nspins))
         ENDIF

         DO ispin = 1, nspins
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v_full_blk(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="V_FULL_BLK", &
                                    size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_virt_full/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v_blk(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="V_BLK", &
                                    size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_virt/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="V", &
                                    size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_virt/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_ov_full(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="OV_FULL", &
                                    size_keys=(/almo_mat_dim_occ, almo_mat_dim_virt_full/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_ov(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="OV", &
                                    size_keys=(/almo_mat_dim_occ, almo_mat_dim_virt/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_vo(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="VO", &
                                    size_keys=(/almo_mat_dim_virt, almo_mat_dim_occ/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_x(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="VO", &
                                    size_keys=(/almo_mat_dim_virt, almo_mat_dim_occ/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma_vv(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="SIG_VV", &
                                    size_keys=(/almo_mat_dim_virt, almo_mat_dim_virt/), &
                                    symmetry_new=dbcsr_type_symmetric, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_vv_full_blk(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="VV_FULL_BLK", &
                                    size_keys=(/almo_mat_dim_virt_full, almo_mat_dim_virt_full/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.TRUE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_sigma_vv_blk(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="SIG_VV_BLK", &
                                    size_keys=(/almo_mat_dim_virt, almo_mat_dim_virt/), &
                                    symmetry_new=dbcsr_type_symmetric, &
                                    spin_key=ispin, &
                                    init_domains=.TRUE.)
            CALL dbcsr_init(almo_scf_env%matrix_sigma_vv_sqrt(ispin))
            CALL dbcsr_init(almo_scf_env%matrix_sigma_vv_sqrt_inv(ispin))
            CALL dbcsr_create(almo_scf_env%matrix_sigma_vv_sqrt(ispin), &
                              template=almo_scf_env%matrix_sigma_vv(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(almo_scf_env%matrix_sigma_vv_sqrt_inv(ispin), &
                              template=almo_scf_env%matrix_sigma_vv(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN
               CALL matrix_almo_create(matrix_new=almo_scf_env%opt_k_t_rr(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="OPT_K_U_RR", &
                                       size_keys=(/almo_mat_dim_virt, almo_mat_dim_virt/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.FALSE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_vv_disc(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="VV_DISC", &
                                       size_keys=(/almo_mat_dim_virt_disc, almo_mat_dim_virt_disc/), &
                                       symmetry_new=dbcsr_type_symmetric, &
                                       spin_key=ispin, &
                                       init_domains=.FALSE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%opt_k_t_dd(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="OPT_K_U_DD", &
                                       size_keys=(/almo_mat_dim_virt_disc, almo_mat_dim_virt_disc/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.FALSE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_vv_disc_blk(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="VV_DISC_BLK", &
                                       size_keys=(/almo_mat_dim_virt_disc, almo_mat_dim_virt_disc/), &
                                       symmetry_new=dbcsr_type_symmetric, &
                                       spin_key=ispin, &
                                       init_domains=.TRUE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_k_blk(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="K_BLK", &
                                       size_keys=(/almo_mat_dim_virt_disc, almo_mat_dim_virt/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.TRUE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_k_blk_ones(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="K_BLK_1", &
                                       size_keys=(/almo_mat_dim_virt_disc, almo_mat_dim_virt/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.TRUE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%opt_k_denom(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="OPT_K_DENOM", &
                                       size_keys=(/almo_mat_dim_virt_disc, almo_mat_dim_virt/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.FALSE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_k_tr(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="K_TR", &
                                       size_keys=(/almo_mat_dim_virt, almo_mat_dim_virt_disc/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.FALSE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v_disc_blk(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="V_DISC_BLK", &
                                       size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_virt_disc/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.FALSE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_v_disc(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="V_DISC", &
                                       size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_virt_disc/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.FALSE.)
               CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_ov_disc(ispin), &
                                       matrix_qs=matrix_s0, &
                                       almo_scf_env=almo_scf_env, &
                                       name_new="OV_DISC", &
                                       size_keys=(/almo_mat_dim_occ, almo_mat_dim_virt_disc/), &
                                       symmetry_new=dbcsr_type_no_symmetry, &
                                       spin_key=ispin, &
                                       init_domains=.FALSE.)

            ENDIF ! end need_discarded_virtuals

         ENDDO ! spin
      ENDIF

      ! create matrices of orbital energies if necessary
      IF (almo_scf_env%need_orbital_energies) THEN
         ALLOCATE (almo_scf_env%matrix_eoo(nspins))
         ALLOCATE (almo_scf_env%matrix_evv_full(nspins))
         DO ispin = 1, nspins
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_eoo(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="E_OCC", &
                                    size_keys=(/almo_mat_dim_occ, almo_mat_dim_occ/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
            CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_evv_full(ispin), &
                                    matrix_qs=matrix_s0, &
                                    almo_scf_env=almo_scf_env, &
                                    name_new="E_VIRT", &
                                    size_keys=(/almo_mat_dim_virt_full, almo_mat_dim_virt_full/), &
                                    symmetry_new=dbcsr_type_no_symmetry, &
                                    spin_key=ispin, &
                                    init_domains=.FALSE.)
         ENDDO
      ENDIF

      ! Density and KS matrices
      ALLOCATE (almo_scf_env%matrix_p(nspins))
      ALLOCATE (almo_scf_env%matrix_p_blk(nspins))
      ALLOCATE (almo_scf_env%matrix_ks(nspins))
      ALLOCATE (almo_scf_env%matrix_ks_blk(nspins))
      IF (almo_scf_env%need_previous_ks) &
         ALLOCATE (almo_scf_env%matrix_ks_0deloc(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_init(almo_scf_env%matrix_p(ispin))
         ! RZK-warning copy with symmery but remember that this might cause problems
         CALL dbcsr_create(almo_scf_env%matrix_p(ispin), &
                           template=almo_scf_env%matrix_s(1), &
                           matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_init(almo_scf_env%matrix_ks(ispin))
         CALL dbcsr_create(almo_scf_env%matrix_ks(ispin), &
                           template=almo_scf_env%matrix_s(1), &
                           matrix_type=dbcsr_type_symmetric)
         IF (almo_scf_env%need_previous_ks) THEN
            CALL dbcsr_init(almo_scf_env%matrix_ks_0deloc(ispin))
            CALL dbcsr_create(almo_scf_env%matrix_ks_0deloc(ispin), &
                              template=almo_scf_env%matrix_s(1), &
                              matrix_type=dbcsr_type_symmetric)
         ENDIF
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_p_blk(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="P_BLK", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
                                 symmetry_new=dbcsr_type_symmetric, &
                                 spin_key=ispin, &
                                 init_domains=.TRUE.)
         CALL matrix_almo_create(matrix_new=almo_scf_env%matrix_ks_blk(ispin), &
                                 matrix_qs=matrix_s0, &
                                 almo_scf_env=almo_scf_env, &
                                 name_new="KS_BLK", &
                                 size_keys=(/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
                                 symmetry_new=dbcsr_type_symmetric, &
                                 spin_key=ispin, &
                                 init_domains=.TRUE.)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE almo_scf_env_create_matrices

! **************************************************************************************************
!> \brief clean up procedures for almo scf
!> \param almo_scf_env ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_clean_up(almo_scf_env)

      TYPE(almo_scf_env_type)                            :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_clean_up', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%mepos == logger%para_env%source) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! release matrices
      CALL dbcsr_release(almo_scf_env%matrix_s(1))
      CALL dbcsr_release(almo_scf_env%matrix_s_blk(1))
      IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_diag) THEN
         CALL dbcsr_release(almo_scf_env%matrix_s_blk_sqrt_inv(1))
         CALL dbcsr_release(almo_scf_env%matrix_s_blk_sqrt(1))
      ELSE IF (almo_scf_env%almo_update_algorithm .EQ. almo_scf_dm_sign) THEN
         CALL dbcsr_release(almo_scf_env%matrix_s_blk_inv(1))
      ENDIF
      DO ispin = 1, almo_scf_env%nspins
         CALL dbcsr_release(almo_scf_env%quench_t(ispin))
         CALL dbcsr_release(almo_scf_env%quench_t_blk(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_t_blk(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_err_blk(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_err_xx(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_t_tr(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_sigma(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_sigma_blk(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_sigma_inv_0deloc(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_sigma_inv(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_t(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_sigma_sqrt(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_sigma_sqrt_inv(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_p(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_ks(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_p_blk(ispin))
         CALL dbcsr_release(almo_scf_env%matrix_ks_blk(ispin))
         IF (almo_scf_env%need_previous_ks) THEN
            CALL dbcsr_release(almo_scf_env%matrix_ks_0deloc(ispin))
         ENDIF
         IF (almo_scf_env%need_virtuals) THEN
            CALL dbcsr_release(almo_scf_env%matrix_v_blk(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_v_full_blk(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_v(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_vo(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_x(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_ov(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_ov_full(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_sigma_vv(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_sigma_vv_blk(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_sigma_vv_sqrt(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_sigma_vv_sqrt_inv(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_vv_full_blk(ispin))
            IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN
               CALL dbcsr_release(almo_scf_env%matrix_k_tr(ispin))
               CALL dbcsr_release(almo_scf_env%matrix_k_blk(ispin))
               CALL dbcsr_release(almo_scf_env%matrix_k_blk_ones(ispin))
               CALL dbcsr_release(almo_scf_env%matrix_v_disc(ispin))
               CALL dbcsr_release(almo_scf_env%matrix_v_disc_blk(ispin))
               CALL dbcsr_release(almo_scf_env%matrix_ov_disc(ispin))
               CALL dbcsr_release(almo_scf_env%matrix_vv_disc_blk(ispin))
               CALL dbcsr_release(almo_scf_env%matrix_vv_disc(ispin))
               CALL dbcsr_release(almo_scf_env%opt_k_t_dd(ispin))
               CALL dbcsr_release(almo_scf_env%opt_k_t_rr(ispin))
               CALL dbcsr_release(almo_scf_env%opt_k_denom(ispin))
            ENDIF
         ENDIF
         IF (almo_scf_env%need_orbital_energies) THEN
            CALL dbcsr_release(almo_scf_env%matrix_eoo(ispin))
            CALL dbcsr_release(almo_scf_env%matrix_evv_full(ispin))
         ENDIF
      ENDDO

      ! deallocate matrices
      DEALLOCATE (almo_scf_env%matrix_p)
      DEALLOCATE (almo_scf_env%matrix_p_blk)
      DEALLOCATE (almo_scf_env%matrix_ks)
      DEALLOCATE (almo_scf_env%matrix_ks_blk)
      DEALLOCATE (almo_scf_env%matrix_t_blk)
      DEALLOCATE (almo_scf_env%matrix_err_blk)
      DEALLOCATE (almo_scf_env%matrix_err_xx)
      DEALLOCATE (almo_scf_env%matrix_t)
      DEALLOCATE (almo_scf_env%matrix_t_tr)
      DEALLOCATE (almo_scf_env%matrix_sigma)
      DEALLOCATE (almo_scf_env%matrix_sigma_blk)
      DEALLOCATE (almo_scf_env%matrix_sigma_inv_0deloc)
      DEALLOCATE (almo_scf_env%matrix_sigma_sqrt)
      DEALLOCATE (almo_scf_env%matrix_sigma_sqrt_inv)
      DEALLOCATE (almo_scf_env%matrix_sigma_inv)
      DEALLOCATE (almo_scf_env%quench_t)
      DEALLOCATE (almo_scf_env%quench_t_blk)
      IF (almo_scf_env%need_virtuals) THEN
         DEALLOCATE (almo_scf_env%matrix_v_blk)
         DEALLOCATE (almo_scf_env%matrix_v_full_blk)
         DEALLOCATE (almo_scf_env%matrix_v)
         DEALLOCATE (almo_scf_env%matrix_vo)
         DEALLOCATE (almo_scf_env%matrix_x)
         DEALLOCATE (almo_scf_env%matrix_ov)
         DEALLOCATE (almo_scf_env%matrix_ov_full)
         DEALLOCATE (almo_scf_env%matrix_sigma_vv)
         DEALLOCATE (almo_scf_env%matrix_sigma_vv_blk)
         DEALLOCATE (almo_scf_env%matrix_sigma_vv_sqrt)
         DEALLOCATE (almo_scf_env%matrix_sigma_vv_sqrt_inv)
         DEALLOCATE (almo_scf_env%matrix_vv_full_blk)
         IF (almo_scf_env%deloc_truncate_virt .NE. virt_full) THEN
            DEALLOCATE (almo_scf_env%matrix_k_tr)
            DEALLOCATE (almo_scf_env%matrix_k_blk)
            DEALLOCATE (almo_scf_env%matrix_v_disc)
            DEALLOCATE (almo_scf_env%matrix_v_disc_blk)
            DEALLOCATE (almo_scf_env%matrix_ov_disc)
            DEALLOCATE (almo_scf_env%matrix_vv_disc_blk)
            DEALLOCATE (almo_scf_env%matrix_vv_disc)
            DEALLOCATE (almo_scf_env%matrix_k_blk_ones)
            DEALLOCATE (almo_scf_env%opt_k_t_dd)
            DEALLOCATE (almo_scf_env%opt_k_t_rr)
            DEALLOCATE (almo_scf_env%opt_k_denom)
         ENDIF
      ENDIF
      IF (almo_scf_env%need_previous_ks) THEN
         DEALLOCATE (almo_scf_env%matrix_ks_0deloc)
      ENDIF
      IF (almo_scf_env%need_orbital_energies) THEN
         DEALLOCATE (almo_scf_env%matrix_eoo)
         DEALLOCATE (almo_scf_env%matrix_evv_full)
      ENDIF

      ! clean up other variables
      DO ispin = 1, almo_scf_env%nspins
         CALL release_submatrices( &
            almo_scf_env%domain_preconditioner(:, ispin))
         CALL release_submatrices(almo_scf_env%domain_s_inv(:, ispin))
         CALL release_submatrices(almo_scf_env%domain_s_sqrt_inv(:, ispin))
         CALL release_submatrices(almo_scf_env%domain_s_sqrt(:, ispin))
         CALL release_submatrices(almo_scf_env%domain_ks_xx(:, ispin))
         CALL release_submatrices(almo_scf_env%domain_t(:, ispin))
         CALL release_submatrices(almo_scf_env%domain_err(:, ispin))
         CALL release_submatrices(almo_scf_env%domain_r_down_up(:, ispin))
      ENDDO
      DEALLOCATE (almo_scf_env%domain_preconditioner)
      DEALLOCATE (almo_scf_env%domain_s_inv)
      DEALLOCATE (almo_scf_env%domain_s_sqrt_inv)
      DEALLOCATE (almo_scf_env%domain_s_sqrt)
      DEALLOCATE (almo_scf_env%domain_ks_xx)
      DEALLOCATE (almo_scf_env%domain_t)
      DEALLOCATE (almo_scf_env%domain_err)
      DEALLOCATE (almo_scf_env%domain_r_down_up)
      DO ispin = 1, almo_scf_env%nspins
         DEALLOCATE (almo_scf_env%domain_map(ispin)%pairs)
         DEALLOCATE (almo_scf_env%domain_map(ispin)%index1)
      ENDDO
      DEALLOCATE (almo_scf_env%domain_map)
      DEALLOCATE (almo_scf_env%domain_index_of_ao)
      DEALLOCATE (almo_scf_env%domain_index_of_atom)
      DEALLOCATE (almo_scf_env%first_atom_of_domain)
      DEALLOCATE (almo_scf_env%last_atom_of_domain)
      DEALLOCATE (almo_scf_env%nbasis_of_domain)
      DEALLOCATE (almo_scf_env%nocc_of_domain)
      DEALLOCATE (almo_scf_env%nvirt_full_of_domain)
      DEALLOCATE (almo_scf_env%nvirt_of_domain)
      DEALLOCATE (almo_scf_env%nvirt_disc_of_domain)
      DEALLOCATE (almo_scf_env%mu_of_domain)
      DEALLOCATE (almo_scf_env%cpu_of_domain)
      DEALLOCATE (almo_scf_env%charge_of_domain)
      DEALLOCATE (almo_scf_env%multiplicity_of_domain)

      DEALLOCATE (almo_scf_env%domain_index_of_ao_block)
      DEALLOCATE (almo_scf_env%domain_index_of_mo_block)

      CALL cp_para_env_release(almo_scf_env%para_env)
      CALL cp_blacs_env_release(almo_scf_env%blacs_env)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_clean_up

! **************************************************************************************************
!> \brief Do post scf calculations with ALMO
!>        WARNING: ALMO post scf calculation may not work for certain quantities,
!>        like forces, since ALMO wave function is only 'partially' optimized
!> \param qs_env ...
!> \par History
!>       2016.12 created [Yifei Shi]
!> \author Yifei Shi
! **************************************************************************************************
   SUBROUTINE almo_post_scf_compute_properties(qs_env)
      ! SUBROUTINE almo_post_scf_compute_properties(qs_env, almo_scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_post_scf_compute_properties', &
         routineP = moduleN//':'//routineN

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: dft_section, input

!TYPE(almo_scf_env_type)                            :: almo_scf_env

      CALL get_qs_env(qs_env, scf_env=scf_env, input=input)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL qs_scf_compute_properties(qs_env, dft_section)

   END SUBROUTINE almo_post_scf_compute_properties

END MODULE almo_scf

