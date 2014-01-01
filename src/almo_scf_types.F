!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types for all ALMO-based methods
!> \par History
!>       2011.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE almo_scf_types
  USE almo_scf_aux2_types,             ONLY: object01_type,&
                                             object02_type
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

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

  PUBLIC :: almo_objectM1_type

  TYPE almo_objectM1_type
    
    TYPE(cp_para_env_type),  POINTER  :: para_env
    TYPE(cp_blacs_env_type), POINTER  :: blacs_env

    INTEGER :: nspins, nelectrons_total, naos
    INTEGER :: natoms, nmolecules
    INTEGER, DIMENSION(2) :: nelectrons_spin

    INTEGER :: ndomains

    INTEGER, DIMENSION(:,:), ALLOCATABLE                :: nocc_of_domain
    INTEGER, DIMENSION(:), ALLOCATABLE                  :: nbasis_of_domain
    INTEGER, DIMENSION(:,:), ALLOCATABLE                :: nvirt_full_of_domain
    INTEGER, DIMENSION(:,:), ALLOCATABLE                :: nvirt_of_domain
    INTEGER, DIMENSION(:,:), ALLOCATABLE                :: nvirt_disc_of_domain
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE          :: mu_of_domain
    INTEGER, DIMENSION(:), ALLOCATABLE                  :: first_atom_of_domain
    INTEGER, DIMENSION(:), ALLOCATABLE        :: domain_index_of_ao
    INTEGER, DIMENSION(:), ALLOCATABLE        :: domain_index_of_atom


INTEGER, DIMENSION(:), ALLOCATABLE                  :: charge_of_domain
   
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE   :: quench_t 
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE   :: quench_t_blk
    TYPE(object02_type), DIMENSION(:), ALLOCATABLE :: domain_map
    
    INTEGER :: domain_layout_mos, domain_layout_aos
    INTEGER :: constraint_type

    INTEGER :: mat_distr_aos, mat_distr_mos
    INTEGER, DIMENSION(:), ALLOCATABLE :: domain_index_of_ao_block
    INTEGER, DIMENSION(:), ALLOCATABLE :: domain_index_of_mo_block

    LOGICAL              :: need_previous_ks
    LOGICAL              :: need_virtuals
    LOGICAL              :: need_orbital_energies
    LOGICAL              :: s_inv_done
    LOGICAL              :: s_sqrt_done
    REAL(KIND=dp)        :: almo_scf_energy
    LOGICAL               :: orthogonal_basis, fixed_mu
  
    REAL(KIND=dp)         :: eps_filter
    REAL(KIND=dp)         :: eps_prev_guess

    REAL(KIND=dp)         :: eps_scf_bd
    REAL(KIND=dp)         :: mixing_fraction
    REAL(KIND=dp)         :: mu
    INTEGER               :: max_scf_bd
    INTEGER               :: counter04_bd
    INTEGER               :: almo_update_algorithm

    REAL(KIND=dp)         :: eps_scf_q
    REAL(KIND=dp)         :: eps_lin_search
    INTEGER               :: max_scf_q
    INTEGER               :: scf_conjugator
    INTEGER               :: outer_max_scf_q
    INTEGER               :: almo_update_algorithm_q

    LOGICAL               :: perturbative_delocalization
    INTEGER               :: quencher_radius_type
    REAL(KIND=dp)         :: quencher_r0_factor,&
                             quencher_r1_factor,&
                             quencher_s0,&
                             quencher_s1,&
                             envelope_amplitude

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

    INTEGER         :: opt_k_conjugator
    REAL(KIND=dp)   :: opt_k_eps_convergence
    REAL(KIND=dp)   :: opt_k_trial_step_size
    INTEGER         :: opt_k_max_iter
    INTEGER         :: opt_k_outer_max_iter
    REAL(KIND=dp)   :: opt_k_trial_step_size_multiplier
    INTEGER         :: opt_k_conj_iter_start
    INTEGER         :: opt_k_prec_iter_start
    INTEGER         :: opt_k_conj_iter_freq
    INTEGER         :: opt_k_prec_iter_freq
 
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
 
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: matrix_s
    
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: so_i
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: so_s
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: so_si
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: so_b
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: so_bi
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: so_bs
    TYPE(cp_dbcsr_type), DIMENSION(1)   :: so_bsi
    
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: enter_b
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: enter
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: enter_tr
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: rem,&
                                                      rem_i,&
                                                      rem_s,&
                                                      rem_si,&
                                                      rem_b
    
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: gst_b
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: gst_xx

    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: rem_vv,&
                                                      rem_vv_b,&
                                                      rem_vv_s,&
                                                      rem_vv_si

    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: vv_fb,&
                                                      vv_db,&
                                                      vv_d
   
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: vo, ov
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: ov_full,&
                                                      ov_disc
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: xia
    
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: kx_b
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: k_t
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: kx_bo
    
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: v_b,&
                                                      vib,&
                                                      v_fb,&
                                                      v_d,&
                                                      v_db
    
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: hfh
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: hfh_b
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: dpp
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: dpp_b

    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: eoo
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: evv
    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: hfh_con

    TYPE(cp_dbcsr_type), DIMENSION(:), ALLOCATABLE :: opt_k_t_rr,&
                                                      opt_k_t_dd,&
                                                      opt_k_denom

    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: concept_p
    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: concept_soi
    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: concept_ss
    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: concept_ss_inv
    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: concept_hfh_xx
    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: concept_enter
    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: concept_gst
    TYPE(object01_type), DIMENSION(:,:), ALLOCATABLE :: concept_dpp

    INTEGER, DIMENSION(:), ALLOCATABLE                       :: cpu_of_domain

  END TYPE almo_objectM1_type

END MODULE almo_scf_types

