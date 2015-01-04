!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types for all ALMO-based methods
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE almo_scf_types
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE domain_submatrix_types,          ONLY: domain_map_type,&
                                             domain_submatrix_type
  USE input_constants,                 ONLY: &
       cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, cg_hager_zhang, &
       cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero, &
       optimizer_diis, optimizer_pcg, prec_default, prec_zero
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_types'

  INTEGER, PARAMETER, PUBLIC   :: almo_mat_dim_aobasis=1, &
                                  almo_mat_dim_occ=2, &
                                  almo_mat_dim_virt=3, &
                                  almo_mat_dim_virt_full=4, &
                                  almo_mat_dim_domains=5, &
                                  almo_mat_dim_virt_disc=6
  REAL(KIND=dp), PARAMETER, PUBLIC :: almo_max_cutoff_multiplier=2.2_dp

  PUBLIC :: almo_scf_env_type, optimizer_options_type, &
            print_optimizer_options

  TYPE optimizer_options_type 
    
    ! SCF procedure for the block-diagonal ALMOs
    !REAL(KIND=dp)         :: eps_scf_bd -> eps_error
    !INTEGER               :: max_scf_bd -> max_iter
    !INTEGER               :: ndiis_bd -> ndiis
    
    ! SCF procedure for the quenched ALMOs (xALMOs)
    !REAL(KIND=dp)         :: eps_scf_q -> eps_error
    !REAL(KIND=dp)         :: eps_lin_search ->  eps_lin_search_error
    !INTEGER               :: max_scf_q -> max_iter
    !INTEGER               :: scf_conjugator -> conjugator
    !INTEGER               :: outer_max_scf_q -> max_iter for a separate optimizer
    !INTEGER               :: ndiis_q -> ndiis

    REAL(KIND=dp)  :: eps_error,&
                      lin_search_eps_error,&
                      lin_search_step_size_guess

    INTEGER        :: optimizer_type ! diis, pcg, etc.
    
    INTEGER        :: preconditioner,& ! preconditioner type
                      conjugator,& ! conjugator type
                      max_iter,&
                      max_iter_outer_loop,&
                      ndiis ! diis history length

  END TYPE optimizer_options_type


  ! the structure contains general info about the system
  TYPE almo_scf_env_type
    
    TYPE(cp_para_env_type),  POINTER  :: para_env
    TYPE(cp_blacs_env_type), POINTER  :: blacs_env

    INTEGER :: nspins, nelectrons_total, naos
    INTEGER :: natoms, nmolecules
    INTEGER, DIMENSION(2) :: nelectrons_spin

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
    INTEGER :: ndomains

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
    INTEGER, DIMENSION(:,:), ALLOCATABLE                :: nocc_of_domain
    ! Number of basis functions for a given domain
    INTEGER, DIMENSION(:), ALLOCATABLE                  :: nbasis_of_domain
    ! Define number of virtuals for a given domain: nvirt = nbasis - nocc
    INTEGER, DIMENSION(:,:), ALLOCATABLE                :: nvirt_full_of_domain
    ! Define the dimension of truncated virtual subspace for a given domain:
    INTEGER, DIMENSION(:,:), ALLOCATABLE                :: nvirt_of_domain
    ! Define the dimension of discarded virtual subspace for a given domain:
    INTEGER, DIMENSION(:,:), ALLOCATABLE                :: nvirt_disc_of_domain
    ! Each domain has its own mu - "fermi" level
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE          :: mu_of_domain
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
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE   :: quench_t 
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE   :: quench_t_blk
    ! Local array for a compact description of quench_t
    TYPE(domain_map_type), DIMENSION(:), ALLOCATABLE :: domain_map
    
    ! Several special cases for the structure of the group-domain matrix:
    ! 1. The basis functions can be grouped into:
    !    a. molecular sets
    !    b. atomic sets
    ! 2. Electrons can be grouped into:
    !    a. molecular sets
    !    b. atomic sets
    INTEGER :: domain_layout_mos, domain_layout_aos
    ! ALMO  constraint type. 
    INTEGER :: constraint_type

    ! decomposition analysis (energy da, charge da)
    INTEGER :: almo_eda

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
    INTEGER :: mat_distr_aos, mat_distr_mos
    ! Define mappping from a distribution block to a domain
    INTEGER, DIMENSION(:), ALLOCATABLE :: domain_index_of_ao_block
    INTEGER, DIMENSION(:), ALLOCATABLE :: domain_index_of_mo_block

    LOGICAL              :: need_previous_ks
    LOGICAL              :: need_virtuals
    LOGICAL              :: need_orbital_energies
    LOGICAL              :: s_inv_done
    LOGICAL              :: s_sqrt_done
    REAL(KIND=dp)        :: almo_scf_energy
    LOGICAL               :: orthogonal_basis, fixed_mu
  
    ! Controls for the SCF procedure
    REAL(KIND=dp)         :: eps_filter
    INTEGER               :: almo_scf_guess
    REAL(KIND=dp)         :: eps_prev_guess
    INTEGER               :: order_lanczos
    REAL(KIND=dp)         :: eps_lanczos
    INTEGER               :: max_iter_lanczos
    REAL(KIND=dp)         :: mixing_fraction
    REAL(KIND=dp)         :: mu
    ! SCF procedure for the block-diagonal ALMOs
    INTEGER               :: almo_update_algorithm
    ! SCF procedure for the quenched ALMOs (xALMOs)
    INTEGER               :: xalmo_update_algorithm

    ! ALMO SCF delocalization control
    LOGICAL               :: perturbative_delocalization
    INTEGER               :: quencher_radius_type
    REAL(KIND=dp)         :: quencher_r0_factor,&
                             quencher_r1_factor,&
                             !quencher_r0_shift,&
                             !quencher_r1_shift,&
                             quencher_s0,&
                             quencher_s1,&
                             envelope_amplitude

    !!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! MATRICES !!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!

    ! AO overlap NxN
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s_inv
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s_sqrt
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s_sqrt_inv
    ! block-diagonal AO overlap NxN
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s_blk
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s_blk_inv
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s_blk_sqrt
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s_blk_sqrt_inv
    
    ! occupied ALMO coeff NxOCC (alpha,beta - if necessary)
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_t_blk
    ! occupied MO coeff NxOCC (alpha,beta - if necessary)
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_t
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_t_tr
    ! MO overlap OCCxOCC and its inverse (alpha, beta - if necessary)
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_sigma,&
                                                      matrix_sigma_inv,&
                                                      matrix_sigma_sqrt,&
                                                      matrix_sigma_sqrt_inv,&
                                                      matrix_sigma_blk
    
    ! error vector (alpha,beta - if necessary)
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_err_blk
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_err_xx

    ! MO overlap VIRTxVIRT and its derivatives
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_sigma_vv,&
                                                      matrix_sigma_vv_blk,&
                                                      matrix_sigma_vv_sqrt,&
                                                      matrix_sigma_vv_sqrt_inv

    ! template of various VIRT x VIR matrices
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_vv_full_blk,&
                                                      matrix_vv_disc_blk,&
                                                      matrix_vv_disc
   
    ! VIRT-OCC MO overlap
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_vo, matrix_ov
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ov_full,&
                                                      matrix_ov_disc
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_x
    
    ! VIRT_DISC x VIRT_RETAINED 
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_k_blk
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_k_tr
    ! matrix_k_blk_ones is blocked with all elements equal to 1.0
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_k_blk_ones
    
    ! virtual ALMO coeff NxV 
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_v_blk,&
                                                      matrix_v,&
                                                      matrix_v_full_blk,&
                                                      matrix_v_disc,&
                                                      matrix_v_disc_blk
    
    ! kohn-sham matrix (alpha,beta - if necessary)
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_ks_blk
    ! density NxN (alpha,beta - if necessary)
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_p_blk

    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_eoo
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: matrix_evv_full
    TYPE(cp_dbcsr_type), DIMENSION(:), &
       ALLOCATABLE                  :: matrix_ks_almo_scf_converged

    ! preconditioner for k-optimization
    ! RZK-warning: do they have to be stored?
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: opt_k_t_rr,&
                                                      opt_k_t_dd,&
                                                      opt_k_denom

    ! second dimension is spin
    TYPE(domain_submatrix_type), DIMENSION(:,:), ALLOCATABLE :: domain_preconditioner
    TYPE(domain_submatrix_type), DIMENSION(:,:), ALLOCATABLE :: domain_s_inv
    TYPE(domain_submatrix_type), DIMENSION(:,:), ALLOCATABLE :: domain_s_sqrt
    TYPE(domain_submatrix_type), DIMENSION(:,:), ALLOCATABLE :: domain_s_sqrt_inv
    TYPE(domain_submatrix_type), DIMENSION(:,:), ALLOCATABLE :: domain_ks_xx
    TYPE(domain_submatrix_type), DIMENSION(:,:), ALLOCATABLE :: domain_t
    TYPE(domain_submatrix_type), DIMENSION(:,:), ALLOCATABLE :: domain_err
    TYPE(domain_submatrix_type), DIMENSION(:,:), ALLOCATABLE :: domain_r_down_up

    INTEGER, DIMENSION(:), ALLOCATABLE                       :: cpu_of_domain

    ! Options for various optimizers collected neatly
    TYPE(optimizer_options_type)                   :: opt_block_diag_diis
    TYPE(optimizer_options_type)                   :: opt_block_diag_pcg
    TYPE(optimizer_options_type)                   :: opt_xalmo_diis
    TYPE(optimizer_options_type)                   :: opt_xalmo_pcg
    TYPE(optimizer_options_type)                   :: opt_k_pcg
    
    ! keywords that control electron delocalization treatment
    ! RZK-warning: many of these varibles should be collected
    !  into an optimizer_options_type variable
    INTEGER         :: deloc_method
    LOGICAL         :: deloc_use_occ_orbs
    LOGICAL         :: deloc_cayley_use_virt_orbs
    INTEGER         :: deloc_cayley_tensor_type
    LOGICAL         :: deloc_cayley_linear
    INTEGER         :: deloc_cayley_conjugator
    REAL(KIND=dp)   :: deloc_cayley_eps_convergence
    INTEGER         :: deloc_cayley_max_iter
    INTEGER         :: deloc_truncate_virt
    INTEGER         :: deloc_virt_per_domain
    LOGICAL         :: deloc_cayley_occ_precond
    LOGICAL         :: deloc_cayley_vir_precond

    !! keywords that control optimization of retained orbitals
    INTEGER         :: opt_k_conjugator !-> conjugartor
    REAL(KIND=dp)   :: opt_k_eps_convergence !-> eps_error
    REAL(KIND=dp)   :: opt_k_trial_step_size !-> lin_search_step_size_guess
    INTEGER         :: opt_k_max_iter !-> max_iter
    INTEGER         :: opt_k_outer_max_iter !-> max_iter for a separate 'outer' optimizer
    REAL(KIND=dp)   :: opt_k_trial_step_size_multiplier !-> ?
    INTEGER         :: opt_k_conj_iter_start !-> ?
    INTEGER         :: opt_k_prec_iter_start !-> ?
    INTEGER         :: opt_k_conj_iter_freq !-> ?
    INTEGER         :: opt_k_prec_iter_freq !-> ?

    ! development keywords
    INTEGER         :: integer01
    INTEGER         :: integer02
    INTEGER         :: integer03
    INTEGER         :: integer04
    INTEGER         :: integer05
    REAL(KIND=dp)   :: real01 
    REAL(KIND=dp)   :: real02 
    REAL(KIND=dp)   :: real03 
    REAL(KIND=dp)   :: real04 
    REAL(KIND=dp)   :: real05
    LOGICAL         :: logical01
    LOGICAL         :: logical02
    LOGICAL         :: logical03
    LOGICAL         :: logical04
    LOGICAL         :: logical05
 
  END TYPE almo_scf_env_type

  CONTAINS

! *****************************************************************************
!> \brief Prints out the options of an optimizer
!> \param optimizer - options to print
!> \param unit_nr - output stream
!> \param error - error
!> \par History
!>       2014.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE print_optimizer_options(optimizer, unit_nr, error) 

    TYPE(optimizer_options_type), INTENT(IN) :: optimizer
    INTEGER, INTENT(IN)                      :: unit_nr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_optimizer_options', &
      routineP = moduleN//':'//routineN

    CHARACTER(33)                            :: conj_string, prec_string, &
                                                type_string

    IF (unit_nr.gt.0) THEN

       SELECT CASE(optimizer%optimizer_type)
       CASE(optimizer_diis)
          type_string="DIIS"
       CASE(optimizer_pcg)
          type_string="PCG"
       END SELECT

       WRITE(unit_nr,'(T4,A,T48,A33)') "optimizer type:",TRIM(type_string)
       WRITE(unit_nr,'(T4,A,T48,I33)') "maximum iterations:",optimizer%max_iter
       WRITE(unit_nr,'(T4,A,T48,E33.3)') "target error:",optimizer%eps_error

       SELECT CASE(optimizer%optimizer_type)
       CASE(optimizer_diis)
          
          WRITE(unit_nr,'(T4,A,T48,I33)') "maximum DIIS history:",optimizer%ndiis
       
       CASE(optimizer_pcg)
          
          WRITE(unit_nr,'(T4,A,T48,E33.3)') "line search step size guess:",&
                optimizer%lin_search_step_size_guess
          WRITE(unit_nr,'(T4,A,T48,E33.3)') "line search target error:",&
                optimizer%lin_search_eps_error
          WRITE(unit_nr,'(T4,A,T48,I33)') "maximum outer loop iterations:",&
                optimizer%max_iter_outer_loop
          
          SELECT CASE(optimizer%preconditioner)
          CASE(prec_zero)
              prec_string="NONE"
          CASE(prec_default)
              prec_string="0.5 H + 0.5 S"
          END SELECT
          WRITE(unit_nr,'(T4,A,T48,A33)') "preconditioner:",TRIM(prec_string)
          
          SELECT CASE(optimizer%conjugator)
          CASE(cg_zero)
              conj_string="Steepest descent"
          CASE(cg_polak_ribiere)
              conj_string="Polak-Ribiere"
          CASE(cg_fletcher_reeves)
              conj_string="Fletcher-Reeves"
          CASE(cg_hestenes_stiefel)
              conj_string="Hestenes-Stiefel"
          CASE(cg_fletcher)
              conj_string="Fletcher"
          CASE(cg_liu_storey)
              conj_string="Liu-Storey"
          CASE(cg_dai_yuan)
              conj_string="Dai-Yuan"
          CASE(cg_hager_zhang)
              conj_string="Hager-Zhang"
          END SELECT
          WRITE(unit_nr,'(T4,A,T48,A33)') "conjugator:",TRIM(conj_string)

       END SELECT !optimizer type

    ENDIF

  END SUBROUTINE print_optimizer_options

END MODULE almo_scf_types

