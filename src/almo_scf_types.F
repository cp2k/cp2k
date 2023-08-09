!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for all ALMO-based methods
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE almo_scf_types
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE dbcsr_api,                       ONLY: dbcsr_release,&
                                              dbcsr_type
   USE domain_submatrix_types,          ONLY: domain_map_type,&
                                              domain_submatrix_type
   USE input_constants,                 ONLY: &
        cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, cg_hager_zhang, cg_hestenes_stiefel, &
        cg_liu_storey, cg_polak_ribiere, cg_zero, optimizer_diis, optimizer_pcg, optimizer_trustr, &
        trustr_cauchy, trustr_dogleg, trustr_steihaug, xalmo_prec_domain, xalmo_prec_full, &
        xalmo_prec_zero
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_types'

   INTEGER, PARAMETER, PUBLIC   :: almo_mat_dim_aobasis = 1, &
                                   almo_mat_dim_occ = 2, &
                                   almo_mat_dim_virt = 3, &
                                   almo_mat_dim_virt_full = 4, &
                                   almo_mat_dim_domains = 5, &
                                   almo_mat_dim_virt_disc = 6
   REAL(KIND=dp), PARAMETER, PUBLIC :: almo_max_cutoff_multiplier = 2.2_dp

   PUBLIC :: almo_scf_env_type, optimizer_options_type, &
             print_optimizer_options, almo_scf_env_release, &
             almo_scf_history_type

   ! methods that add penalty terms to the energy functional
   TYPE penalty_type

      REAL(KIND=dp)      :: final_determinant = 0.0_dp, penalty_strength = 0.0_dp, &
                            determinant_tolerance = 0.0_dp, penalty_strength_dec_factor = 0.0_dp, &
                            compactification_filter_start = 0.0_dp
      INTEGER            :: operator_type = 0
      LOGICAL            :: virtual_nlmos = .FALSE.

   END TYPE penalty_type

   ! almo-based electronic structure analysis
   TYPE almo_analysis_type

      ! switch analysis on/off
      LOGICAL :: do_analysis = .FALSE.

      INTEGER :: frozen_mo_energy_term = 0

   END TYPE almo_analysis_type

   TYPE optimizer_options_type

      REAL(KIND=dp)  :: eps_error = 0.0_dp, &
                        eps_error_early = 0.0_dp, &
                        lin_search_eps_error = 0.0_dp, &
                        lin_search_step_size_guess = 0.0_dp, &
                        rho_do_not_update = 0.0_dp, &
                        model_grad_norm_ratio = 0.0_dp, &
                        initial_trust_radius = 0.0_dp, &
                        max_trust_radius = 0.0_dp, &
                        neglect_threshold = 0.0_dp

      INTEGER        :: optimizer_type = 0 ! diis, pcg, etc.
      TYPE(penalty_type)  :: opt_penalty = penalty_type()

      INTEGER        :: preconditioner = 0, & ! preconditioner type
                        conjugator = 0, & ! conjugator type
                        max_iter = 0, &
                        max_iter_early = 0, &
                        max_iter_outer_loop = 0, &
                        trustr_algorithm = 0, &
                        ndiis = 0 ! diis history length

      LOGICAL        :: early_stopping_on = .FALSE.

   END TYPE optimizer_options_type

   TYPE almo_scf_history_type
      INTEGER :: istore = 0, nstore = 0
      TYPE(dbcsr_type), DIMENSION(:, :), ALLOCATABLE :: matrix_p_up_down
      !TYPE(dbcsr_type), DIMENSION(:, :), ALLOCATABLE :: matrix_x
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_t
   END TYPE

   ! the structure contains general info about the system
   TYPE almo_scf_env_type

      TYPE(mp_para_env_type), POINTER  :: para_env => NULL()
      TYPE(cp_blacs_env_type), POINTER  :: blacs_env => NULL()

      INTEGER :: nspins = 0, nelectrons_total = 0, naos = 0
      INTEGER :: natoms = 0, nmolecules = 0
      INTEGER, DIMENSION(2) :: nelectrons_spin = 0

      ! Definitions:
      ! I.  Domain - a subset of basis functions (e.g. AOs),
      ! II. Group  - a subset of electrons delocalized within a domain.
      !
      ! The following variables specify the group-domain structure
      ! of the system. Several rules must be obeyed:
      ! 1. There should be no zero domains (i.e. domain contains at least one AO).
      ! 2. There should be no empty domains (i.e. all domains must be populated
      !     by at least one electron).
      ! 3. If two groups are localized within the same domain they are combined
      ! It follows that the number of domains is equal to the number of groups
      !
      ! Number of domains
      INTEGER :: ndomains = 0

      ! List of atoms, whose basis functions are included into the domain.
      ! It is assumed that:
      !   (a) basis functions are localized and atom-labeled,
      !   (b) basis functions are grouped into atomic sets (i.e. if a basis
      !       function on an atom is in domain A then all basis functions on
      !       this atom are in domain A)
      !TYPE(domain_list_type), DIMENSION(:), ALLOCATABLE   :: atom_list_of_domain
      ! List of basis functions included into the domain
      !TYPE(domain_list_type), DIMENSION(:), ALLOCATABLE   :: basis_list_of_domain

      ! Number of electrons of each spin for a given domain (second dim is spin).
      ! Note that some domains can be populated only with alpha or beta electrons.
      INTEGER, DIMENSION(:, :), ALLOCATABLE                :: nocc_of_domain
      ! Number of basis functions for a given domain
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: nbasis_of_domain
      ! Define number of virtuals for a given domain: nvirt = nbasis - nocc
      INTEGER, DIMENSION(:, :), ALLOCATABLE                :: nvirt_full_of_domain
      ! Define the dimension of truncated virtual subspace for a given domain:
      INTEGER, DIMENSION(:, :), ALLOCATABLE                :: nvirt_of_domain
      ! Define the dimension of discarded virtual subspace for a given domain:
      INTEGER, DIMENSION(:, :), ALLOCATABLE                :: nvirt_disc_of_domain
      ! Each domain has its own mu - "fermi" level
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE          :: mu_of_domain
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: first_atom_of_domain
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: last_atom_of_domain
      ! The following arrays are useful only with non-overlapping domains
      ! RZK-warning generalization is required
      INTEGER, DIMENSION(:), ALLOCATABLE        :: domain_index_of_ao
      INTEGER, DIMENSION(:), ALLOCATABLE        :: domain_index_of_atom

      ! Charge of domains
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: charge_of_domain
      ! Charge of domains
      INTEGER, DIMENSION(:), ALLOCATABLE                  :: multiplicity_of_domain

      ! The matrix contains information about the delocalization of
      ! alpha and beta electrons.
      ! Rows denote basis function, columns denote electrons.
      ! Non-zero (j,i) entry means that electron j can delocalize over
      ! basis function i. 0.0 means no delocalization
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE   :: quench_t
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE   :: quench_t_blk
      ! Local array for a compact description of quench_t
      TYPE(domain_map_type), DIMENSION(:), ALLOCATABLE :: domain_map

      ! Several special cases for the structure of the group-domain matrix:
      ! 1. The basis functions can be grouped into:
      !    a. molecular sets
      !    b. atomic sets
      ! 2. Electrons can be grouped into:
      !    a. molecular sets
      !    b. atomic sets
      INTEGER :: domain_layout_mos = 0, domain_layout_aos = 0
      ! ALMO  constraint type.
      INTEGER :: constraint_type = 0

      ! Desciptors of molecules
      !INTEGER, DIMENSION(:), ALLOCATABLE        :: molecule_index_of_atom
      !INTEGER, DIMENSION(:), ALLOCATABLE        :: first_atom_of_molecule
      !INTEGER, DIMENSION(:), ALLOCATABLE        :: nbasis_of_molecule
      !INTEGER, DIMENSION(:,:), ALLOCATABLE      :: nocc_of_molecule
      !INTEGER, DIMENSION(:,:), ALLOCATABLE      :: nvirt_of_molecule
      !REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE :: mu_of_molecule

      ! Descriptors of atoms
      !INTEGER, DIMENSION(:), ALLOCATABLE        :: nbasis_of_atom
      !INTEGER, DIMENSION(:,:), ALLOCATABLE      :: nocc_of_atom
      !INTEGER, DIMENSION(:,:), ALLOCATABLE      :: nvirt_of_atom
      !REAL(KIND=dp),DIMENSION(:,:), ALLOCATABLE :: mu_of_atom

      ! All AO and MO matrices are distributed for parallel computations.
      ! The following flags specify what constitues a block for a parallel
      ! distribution. Both AOs and MOs can be divided into atomic or
      ! molecular blocks. Domain blocks should be equal or larger than
      ! the distribution blocks (otherwise retain_sparsity does not work).
      ! Possible values: almo_mat_distr_atomic, almo_mat_distr_molecular
      INTEGER :: mat_distr_aos = 0, mat_distr_mos = 0
      ! Define mappping from a distribution block to a domain
      INTEGER, DIMENSION(:), ALLOCATABLE :: domain_index_of_ao_block
      INTEGER, DIMENSION(:), ALLOCATABLE :: domain_index_of_mo_block

      LOGICAL              :: need_previous_ks = .FALSE.
      LOGICAL              :: need_virtuals = .FALSE.
      LOGICAL              :: need_orbital_energies = .FALSE.
      LOGICAL              :: s_inv_done = .FALSE.
      LOGICAL              :: s_sqrt_done = .FALSE.
      REAL(KIND=dp)        :: almo_scf_energy = 0.0_dp
      LOGICAL              :: orthogonal_basis = .FALSE., fixed_mu = .FALSE.
      LOGICAL              :: return_orthogonalized_mos = .FALSE., construct_nlmos = .FALSE.

      !! Smearing control
      !! smear flag allow to retrieve eigenvalues in almo_scf with diag algorithm and create occupation-scaled ALMO orbitals
      LOGICAL               :: smear = .FALSE.
      !! store relevant smearing parameters
      REAL(KIND=dp)         :: smear_e_temp = 0.0_dp !! electronic temperature, required for Fermi-Dirac
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: kTS !! electronic entropy contribution of each spin system
      !! mo_energies(imo, ispin) stores the eigenvalue corresponding to the orbital imo with spin ispin
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: mo_energies
      !! since S-ALMO creates partially occupied orbitals, there is a need to store the real number of electron-pairs
      !! of each spin and for each fragment
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: real_ne_of_domain

      ! Controls for the SCF procedure
      REAL(KIND=dp)         :: eps_filter = 0.0_dp
      INTEGER               :: xalmo_trial_wf = 0
      INTEGER               :: almo_scf_guess = 0
      REAL(KIND=dp)         :: eps_prev_guess = 0.0_dp
      INTEGER               :: order_lanczos = 0
      REAL(KIND=dp)         :: matrix_iter_eps_error_factor = 0.0_dp
      REAL(KIND=dp)         :: eps_lanczos = 0.0_dp
      INTEGER               :: max_iter_lanczos = 0
      REAL(KIND=dp)         :: mixing_fraction = 0.0_dp
      REAL(KIND=dp)         :: mu = 0.0_dp
      ! SCF procedure for the block-diagonal ALMOs
      INTEGER               :: almo_update_algorithm = 0
      ! SCF procedure for the quenched ALMOs (xALMOs)
      INTEGER               :: xalmo_update_algorithm = 0
      ! mo overlap inversion algorithm
      INTEGER               :: sigma_inv_algorithm = 0

      ! Determinant of the ALMO overlap matrix
      REAL(KIND=dp)         :: overlap_determinant = 0.0_dp

      ! ALMO SCF delocalization control
      LOGICAL               :: perturbative_delocalization = .FALSE.
      INTEGER               :: quencher_radius_type = 0
      REAL(KIND=dp)         :: quencher_r0_factor = 0.0_dp, &
                               quencher_r1_factor = 0.0_dp, &
                               !quencher_r0_shift,&
                               !quencher_r1_shift,&
                               quencher_s0 = 0.0_dp, &
                               quencher_s1 = 0.0_dp, &
                               envelope_amplitude = 0.0_dp

      ! guess options
      ! This prevents a bug in GCC 8/9
      TYPE(almo_scf_history_type) :: almo_history = almo_scf_history_type(matrix_p_up_down=null(), matrix_t=null())
      TYPE(almo_scf_history_type) :: xalmo_history = almo_scf_history_type(matrix_p_up_down=null(), matrix_t=null())
      INTEGER :: almo_extrapolation_order = 0
      INTEGER :: xalmo_extrapolation_order = 0

      ! forces
      LOGICAL :: calc_forces = .FALSE.

      !!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! MATRICES !!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!

      ! AO overlap NxN
      TYPE(dbcsr_type), DIMENSION(1)   :: matrix_s
      TYPE(dbcsr_type), DIMENSION(1)   :: matrix_s_inv
      TYPE(dbcsr_type), DIMENSION(1)   :: matrix_s_sqrt
      TYPE(dbcsr_type), DIMENSION(1)   :: matrix_s_sqrt_inv
      ! block-diagonal AO overlap NxN
      TYPE(dbcsr_type), DIMENSION(1)   :: matrix_s_blk
      TYPE(dbcsr_type), DIMENSION(1)   :: matrix_s_blk_inv
      TYPE(dbcsr_type), DIMENSION(1)   :: matrix_s_blk_sqrt
      TYPE(dbcsr_type), DIMENSION(1)   :: matrix_s_blk_sqrt_inv

      ! occupied ALMO coeff NxOCC (alpha,beta - if necessary)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_t_blk
      ! occupied MO coeff NxOCC (alpha,beta - if necessary)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_t
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_t_tr
      ! MO overlap OCCxOCC and its inverse (alpha, beta - if necessary)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_sigma, &
                                                     matrix_sigma_inv, &
                                                     matrix_sigma_sqrt, &
                                                     matrix_sigma_sqrt_inv, &
                                                     matrix_sigma_blk, &
                                                     matrix_sigma_inv_0deloc

      ! error vector (alpha,beta - if necessary)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_err_blk
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_err_xx

      ! MO overlap VIRTxVIRT and its derivatives
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_sigma_vv, &
                                                     matrix_sigma_vv_blk, &
                                                     matrix_sigma_vv_sqrt, &
                                                     matrix_sigma_vv_sqrt_inv

      ! template of various VIRT x VIR matrices
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_vv_full_blk, &
                                                     matrix_vv_disc_blk, &
                                                     matrix_vv_disc

      ! VIRT-OCC MO overlap
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_vo, matrix_ov
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ov_full, &
                                                     matrix_ov_disc
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_x

      ! VIRT_DISC x VIRT_RETAINED
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_k_blk
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_k_tr
      ! matrix_k_blk_ones is blocked with all elements equal to 1.0
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_k_blk_ones

      ! virtual ALMO coeff NxV
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_v_blk, &
                                                     matrix_v, &
                                                     matrix_v_full_blk, &
                                                     matrix_v_disc, &
                                                     matrix_v_disc_blk

      ! kohn-sham matrix (alpha,beta - if necessary)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks
      ! the diff between ks_blk and ks_0deloc is that blk is a blocked matrix
      ! 0deloc stores the matrix that correponds to zero-delocalization state
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks_blk
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks_0deloc
      ! density NxN (alpha,beta - if necessary)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p_blk

      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_eoo
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_evv_full

      ! preconditioner for k-optimization
      ! RZK-warning: do they have to be stored?
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: opt_k_t_rr, &
                                                     opt_k_t_dd, &
                                                     opt_k_denom

      ! second dimension is spin
      TYPE(domain_submatrix_type), DIMENSION(:, :), ALLOCATABLE :: domain_preconditioner
      TYPE(domain_submatrix_type), DIMENSION(:, :), ALLOCATABLE :: domain_s_inv
      TYPE(domain_submatrix_type), DIMENSION(:, :), ALLOCATABLE :: domain_s_sqrt
      TYPE(domain_submatrix_type), DIMENSION(:, :), ALLOCATABLE :: domain_s_sqrt_inv
      TYPE(domain_submatrix_type), DIMENSION(:, :), ALLOCATABLE :: domain_ks_xx
      TYPE(domain_submatrix_type), DIMENSION(:, :), ALLOCATABLE :: domain_t
      TYPE(domain_submatrix_type), DIMENSION(:, :), ALLOCATABLE :: domain_err
      TYPE(domain_submatrix_type), DIMENSION(:, :), ALLOCATABLE :: domain_r_down_up

      INTEGER, DIMENSION(:), ALLOCATABLE                       :: cpu_of_domain

      ! Options for various subsection options collected neatly
      TYPE(almo_analysis_type)                       :: almo_analysis = almo_analysis_type()

      ! Options for various optimizers collected neatly
      TYPE(optimizer_options_type)                   :: opt_block_diag_diis = optimizer_options_type()
      TYPE(optimizer_options_type)                   :: opt_block_diag_pcg = optimizer_options_type()
      TYPE(optimizer_options_type)                   :: opt_xalmo_diis = optimizer_options_type()
      TYPE(optimizer_options_type)                   :: opt_xalmo_pcg = optimizer_options_type()
      TYPE(optimizer_options_type)                   :: opt_xalmo_trustr = optimizer_options_type()
      TYPE(optimizer_options_type)                   :: opt_nlmo_pcg = optimizer_options_type()
      TYPE(optimizer_options_type)                   :: opt_block_diag_trustr = optimizer_options_type()
      TYPE(optimizer_options_type)                   :: opt_xalmo_newton_pcg_solver = optimizer_options_type()
      TYPE(optimizer_options_type)                   :: opt_k_pcg = optimizer_options_type()

      ! keywords that control electron delocalization treatment
      ! RZK-warning: many of these varibles should be collected
      !  into an optimizer_options_type variable
      INTEGER         :: deloc_method = 0
      LOGICAL         :: deloc_use_occ_orbs = .FALSE.
      LOGICAL         :: deloc_cayley_use_virt_orbs = .FALSE.
      INTEGER         :: deloc_cayley_tensor_type = 0
      LOGICAL         :: deloc_cayley_linear = .FALSE.
      INTEGER         :: deloc_cayley_conjugator = 0
      REAL(KIND=dp)   :: deloc_cayley_eps_convergence = 0.0_dp
      INTEGER         :: deloc_cayley_max_iter = 0
      INTEGER         :: deloc_truncate_virt = 0
      INTEGER         :: deloc_virt_per_domain = 0
      LOGICAL         :: deloc_cayley_occ_precond = .FALSE.
      LOGICAL         :: deloc_cayley_vir_precond = .FALSE.

      !! keywords that control optimization of retained orbitals
      INTEGER         :: opt_k_conjugator = 0 !-> conjugartor
      REAL(KIND=dp)   :: opt_k_eps_convergence = 0.0_dp !-> eps_error
      REAL(KIND=dp)   :: opt_k_trial_step_size = 0.0_dp !-> lin_search_step_size_guess
      INTEGER         :: opt_k_max_iter = 0 !-> max_iter
      INTEGER         :: opt_k_outer_max_iter = 0 !-> max_iter for a separate 'outer' optimizer
      REAL(KIND=dp)   :: opt_k_trial_step_size_multiplier = 0.0_dp !-> ?
      INTEGER         :: opt_k_conj_iter_start = 0 !-> ?
      INTEGER         :: opt_k_prec_iter_start = 0 !-> ?
      INTEGER         :: opt_k_conj_iter_freq = 0 !-> ?
      INTEGER         :: opt_k_prec_iter_freq = 0 !-> ?

      ! development keywords
      INTEGER         :: integer01 = 0
      INTEGER         :: integer02 = 0
      INTEGER         :: integer03 = 0
      INTEGER         :: integer04 = 0
      INTEGER         :: integer05 = 0
      REAL(KIND=dp)   :: real01 = 0.0_dp
      REAL(KIND=dp)   :: real02 = 0.0_dp
      REAL(KIND=dp)   :: real03 = 0.0_dp
      REAL(KIND=dp)   :: real04 = 0.0_dp
      REAL(KIND=dp)   :: real05 = 0.0_dp
      LOGICAL         :: logical01 = .FALSE.
      LOGICAL         :: logical02 = .FALSE.
      LOGICAL         :: logical03 = .FALSE.
      LOGICAL         :: logical04 = .FALSE.
      LOGICAL         :: logical05 = .FALSE.

   END TYPE almo_scf_env_type

CONTAINS

! **************************************************************************************************
!> \brief Prints out the options of an optimizer
!> \param optimizer   options to print
!> \param unit_nr   output stream
!> \par History
!>       2014.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE print_optimizer_options(optimizer, unit_nr)

      TYPE(optimizer_options_type), INTENT(IN)           :: optimizer
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(33)                                      :: conj_string, prec_string, type_string

      IF (unit_nr .GT. 0) THEN

         SELECT CASE (optimizer%optimizer_type)
         CASE (optimizer_diis)
            type_string = "DIIS"
         CASE (optimizer_pcg)
            type_string = "PCG"
         CASE (optimizer_trustr)
            type_string = "TRUST REGION"
         END SELECT

         WRITE (unit_nr, '(T4,A,T48,A33)') "optimizer type:", TRIM(type_string)
         WRITE (unit_nr, '(T4,A,T48,I33)') "maximum iterations:", optimizer%max_iter
         WRITE (unit_nr, '(T4,A,T48,E33.3)') "target error:", optimizer%eps_error

         IF (optimizer%optimizer_type .EQ. optimizer_diis) THEN

            WRITE (unit_nr, '(T4,A,T48,I33)') "maximum DIIS history:", optimizer%ndiis

         END IF

         IF (optimizer%optimizer_type .EQ. optimizer_trustr .OR. &
             optimizer%optimizer_type .EQ. optimizer_pcg) THEN

            WRITE (unit_nr, '(T4,A,T48,I33)') "maximum outer loop iterations:", &
               optimizer%max_iter_outer_loop

            SELECT CASE (optimizer%preconditioner)
            CASE (xalmo_prec_zero)
               prec_string = "NONE"
            CASE (xalmo_prec_domain)
               prec_string = "0.5 KS + 0.5 S, DOMAINS"
            CASE (xalmo_prec_full)
               prec_string = "0.5 KS + 0.5 S, FULL"
            END SELECT
            WRITE (unit_nr, '(T4,A,T48,A33)') "preconditioner:", TRIM(prec_string)

            SELECT CASE (optimizer%conjugator)
            CASE (cg_zero)
               conj_string = "Steepest descent"
            CASE (cg_polak_ribiere)
               conj_string = "Polak-Ribiere"
            CASE (cg_fletcher_reeves)
               conj_string = "Fletcher-Reeves"
            CASE (cg_hestenes_stiefel)
               conj_string = "Hestenes-Stiefel"
            CASE (cg_fletcher)
               conj_string = "Fletcher"
            CASE (cg_liu_storey)
               conj_string = "Liu-Storey"
            CASE (cg_dai_yuan)
               conj_string = "Dai-Yuan"
            CASE (cg_hager_zhang)
               conj_string = "Hager-Zhang"
            END SELECT
            WRITE (unit_nr, '(T4,A,T48,A33)') "conjugator:", TRIM(conj_string)

         END IF

         IF (optimizer%optimizer_type .EQ. optimizer_pcg) THEN

            WRITE (unit_nr, '(T4,A,T48,E33.3)') "line search step size guess:", &
               optimizer%lin_search_step_size_guess
            WRITE (unit_nr, '(T4,A,T48,E33.3)') "line search target error:", &
               optimizer%lin_search_eps_error
            IF (optimizer%neglect_threshold .GT. 0.0_dp) THEN
               WRITE (unit_nr, '(T4,A,T48,E33.3)') "low-curvature threshold:", &
                  optimizer%neglect_threshold
            END IF

         END IF

         IF (optimizer%optimizer_type .EQ. optimizer_trustr) THEN

            SELECT CASE (optimizer%trustr_algorithm)
            CASE (trustr_steihaug)
               conj_string = "Steihaug's CG"
            CASE (trustr_cauchy)
               conj_string = "Cauchy point"
            CASE (trustr_dogleg)
               conj_string = "Dogleg"
            END SELECT
            WRITE (unit_nr, '(T4,A,T48,A33)') "Subproblem algorithm:", TRIM(conj_string)

            WRITE (unit_nr, '(T4,A,T48,E33.3)') "gradient decrease accepted:", &
               optimizer%model_grad_norm_ratio
            WRITE (unit_nr, '(T4,A,T48,E33.3)') "initial trust radius:", &
               optimizer%initial_trust_radius
            WRITE (unit_nr, '(T4,A,T48,E33.3)') "max trust radius:", &
               optimizer%max_trust_radius
            WRITE (unit_nr, '(T4,A,T48,E33.3)') "rho of no update lies between .0 and .25:", &
               optimizer%rho_do_not_update

         END IF

      END IF

   END SUBROUTINE print_optimizer_options

! **************************************************************************************************
!> \brief release the almo scf envirnoment
!> \param almo_scf_env ...
!> \par History
!>       2016.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_env_release(almo_scf_env)
      TYPE(almo_scf_env_type), POINTER                   :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'almo_scf_env_release'

      INTEGER                                            :: handle, ispin, istore

      CALL timeset(routineN, handle)

      ! delete history
      DO ispin = 1, SIZE(almo_scf_env%almo_history%matrix_t)
         DO istore = 1, MIN(almo_scf_env%almo_history%istore, almo_scf_env%almo_history%nstore)
            CALL dbcsr_release(almo_scf_env%almo_history%matrix_p_up_down(ispin, istore))
         END DO
         IF (almo_scf_env%almo_history%istore .GT. 0) &
            CALL dbcsr_release(almo_scf_env%almo_history%matrix_t(ispin))
      END DO
      DEALLOCATE (almo_scf_env%almo_history%matrix_p_up_down)
      DEALLOCATE (almo_scf_env%almo_history%matrix_t)
      ! delete xalmo history
      DO ispin = 1, SIZE(almo_scf_env%xalmo_history%matrix_t)
         DO istore = 1, MIN(almo_scf_env%xalmo_history%istore, almo_scf_env%xalmo_history%nstore)
            CALL dbcsr_release(almo_scf_env%xalmo_history%matrix_p_up_down(ispin, istore))
            !CALL dbcsr_release(almo_scf_env%xalmo_history%matrix_x(ispin, istore))
         END DO
         IF (almo_scf_env%xalmo_history%istore .GT. 0) &
            CALL dbcsr_release(almo_scf_env%xalmo_history%matrix_t(ispin))
      END DO
      DEALLOCATE (almo_scf_env%xalmo_history%matrix_p_up_down)
      !DEALLOCATE (almo_scf_env%xalmo_history%matrix_x)
      DEALLOCATE (almo_scf_env%xalmo_history%matrix_t)

      DEALLOCATE (almo_scf_env)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_env_release

END MODULE almo_scf_types

