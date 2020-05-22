!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types and set/get functions for HFX
!> \par History
!>      04.2008 created [Manuel Guidon]
!>      05.2019 Moved erfc_cutoff to common/mathlib (A. Bussy)
!> \author Manuel Guidon
! **************************************************************************************************
MODULE hfx_types
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE auto_basis,                      ONLY: create_ri_aux_basis_set
   USE basis_set_container_types,       ONLY: add_basis_set_to_container
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type,&
                                              init_orb_basis_set
   USE bibliography,                    ONLY: cite_reference,&
                                              guidon2008,&
                                              guidon2009
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              plane_distance,&
                                              scaled_to_real
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_files,                        ONLY: close_file,&
                                              file_exists,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_batched_contract_finalize, dbcsr_t_create, dbcsr_t_default_distvec, &
        dbcsr_t_destroy, dbcsr_t_distribution_destroy, dbcsr_t_distribution_new, &
        dbcsr_t_distribution_type, dbcsr_t_mp_dims_create, dbcsr_t_pgrid_create, &
        dbcsr_t_pgrid_destroy, dbcsr_t_pgrid_type, dbcsr_t_type
   USE hfx_helpers,                     ONLY: count_cells_perd,&
                                              next_image_cell_perd
   USE input_constants,                 ONLY: &
        do_hfx_auto_shells, do_potential_coulomb, do_potential_gaussian, do_potential_id, &
        do_potential_long, do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_mix_lg, &
        do_potential_short, do_potential_truncated, hfx_ri_do_2c_diag, hfx_ri_do_2c_iter
   USE input_cp2k_hfx,                  ONLY: ri_mo,&
                                              ri_pmat
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE libint_wrapper,                  ONLY: &
        cp_libint_cleanup_eri, cp_libint_cleanup_eri1, cp_libint_init_eri, cp_libint_init_eri1, &
        cp_libint_set_contrdepth, cp_libint_static_cleanup, cp_libint_static_init, cp_libint_t, &
        prim_data_f_size
   USE machine,                         ONLY: m_chdir,&
                                              m_getcwd
   USE mathlib,                         ONLY: erfc_cutoff
   USE message_passing,                 ONLY: mp_cart_create,&
                                              mp_environ
   USE orbital_pointers,                ONLY: nco,&
                                              ncoset,&
                                              nso
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_loc_types,                    ONLY: qs_loc_env_new_type
   USE qs_tensors_types,                ONLY: &
        create_2c_tensor, create_3c_tensor, create_tensor_batches, default_block_size, &
        distribution_3d_create, distribution_3d_destroy, distribution_3d_type, pgf_block_sizes, &
        split_block_sizes
   USE string_utilities,                ONLY: compress
   USE t_c_g0,                          ONLY: free_C0

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: hfx_type, hfx_create, hfx_release, &
             hfx_set_distr_energy, &
             hfx_set_distr_forces, &
             hfx_cell_type, hfx_distribution, &
             hfx_potential_type, hfx_screening_type, &
             hfx_memory_type, hfx_load_balance_type, hfx_general_type, &
             hfx_container_type, hfx_cache_type, &
             hfx_basis_type, parse_memory_section, &
             hfx_init_container, &
             hfx_basis_info_type, hfx_screen_coeff_type, &
             hfx_reset_memory_usage_counter, pair_list_type, pair_list_element_type, &
             pair_set_list_type, hfx_p_kind, hfx_2D_map, hfx_pgf_list, &
             hfx_pgf_product_list, hfx_block_range_type, &
             alloc_containers, dealloc_containers, hfx_task_list_type, init_t_c_g0_lmax, &
             hfx_create_neighbor_cells, hfx_create_basis_types, hfx_release_basis_types, &
             hfx_ri_type, hfx_compression_type

#define CACHE_SIZE 1024
#define BITS_MAX_VAL 6

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_types'
   INTEGER, PARAMETER, PUBLIC                 :: max_atom_block = 32
   INTEGER, PARAMETER, PUBLIC                 :: max_images = 27
   REAL(dp), PARAMETER, PUBLIC                :: log_zero = -1000.0_dp
   REAL(dp), PARAMETER, PUBLIC                :: powell_min_log = -20.0_dp
   REAL(KIND=dp), DIMENSION(0:10), &
      PARAMETER, PUBLIC                       :: mul_fact = (/1.0_dp, &
                                                              1.1781_dp, &
                                                              1.3333_dp, &
                                                              1.4726_dp, &
                                                              1.6000_dp, &
                                                              1.7181_dp, &
                                                              1.8286_dp, &
                                                              1.9328_dp, &
                                                              2.0317_dp, &
                                                              2.1261_dp, &
                                                              2.2165_dp/)

   INTEGER, SAVE                                         :: init_t_c_g0_lmax = -1

!***

! **************************************************************************************************
   TYPE hfx_potential_type
      INTEGER                                  :: potential_type !! 1/r/ erfc(wr)/r ...
      REAL(dp)                                 :: omega !! w
      REAL(dp)                                 :: scale_coulomb !! scaling factor for mixed potential
      REAL(dp)                                 :: scale_longrange !! scaling factor for mixed potential
      REAL(dp)                                 :: scale_gaussian !! scaling factor for mixed potential
      REAL(dp)                                 :: cutoff_radius !! cutoff radius if cutoff potential in use
      CHARACTER(default_path_length)           :: filename
   END TYPE

! **************************************************************************************************
   TYPE hfx_screening_type
      REAL(dp)                                 :: eps_schwarz !! threshold
      REAL(dp)                                 :: eps_schwarz_forces !! threshold
      LOGICAL                                  :: do_p_screening_forces !! screen on P^2 ?
      LOGICAL                                  :: do_initial_p_screening !! screen on initial guess?
   END TYPE

! **************************************************************************************************
   TYPE hfx_memory_type
      INTEGER                                  :: max_memory !! user def max memory MiB
      INTEGER(int_8)                           :: max_compression_counter !! corresponding number of reals
      INTEGER(int_8)                           :: final_comp_counter_energy
      LOGICAL                                  :: do_all_on_the_fly !! max mem == 0 ?
      REAL(dp)                                 :: eps_storage_scaling
      INTEGER                                  :: cache_size
      INTEGER                                  :: bits_max_val
      INTEGER                                  :: actual_memory_usage
      INTEGER                                  :: actual_memory_usage_disk
      INTEGER(int_8)                           :: max_compression_counter_disk
      LOGICAL                                  :: do_disk_storage
      CHARACTER(len=default_path_length)       :: storage_location
      INTEGER(int_8)                           :: ram_counter
      INTEGER(int_8)                           :: ram_counter_forces
      INTEGER(int_8)                           :: size_p_screen
      LOGICAL                                  :: treat_forces_in_core
      LOGICAL                                  :: recalc_forces
   END TYPE

! **************************************************************************************************
   TYPE hfx_periodic_type
      INTEGER                                  :: number_of_shells !! number of periodic image cells
      LOGICAL                                  :: do_periodic !! periodic ?
      INTEGER                                  :: perd(3) !! x,xy,xyz,...
      INTEGER                                  :: mode
      REAL(dp)                                 :: R_max_stress
      INTEGER                                  :: number_of_shells_from_input
   END TYPE

! **************************************************************************************************
   TYPE hfx_load_balance_type
      INTEGER                                  :: nbins
      INTEGER                                  :: block_size
      INTEGER                                  :: nblocks
      LOGICAL                                  :: rtp_redistribute
      LOGICAL                                  :: blocks_initialized
      LOGICAL                                  :: do_randomize
   END TYPE

! **************************************************************************************************
   TYPE hfx_general_type
      REAL(dp)                                 :: fraction !! for hybrids
      LOGICAL                                  :: treat_lsd_in_core
   END TYPE

! **************************************************************************************************
   TYPE hfx_cell_type
      REAL(dp)                                 :: cell(3)
      REAL(dp)                                 :: cell_r(3)
   END TYPE

! **************************************************************************************************
   TYPE hfx_distribution
      INTEGER(int_8)                           :: istart
      INTEGER(int_8)                           :: number_of_atom_quartets
      INTEGER(int_8)                           :: cost
      REAL(KIND=dp)                            :: time_first_scf
      REAL(KIND=dp)                            :: time_other_scf
      REAL(KIND=dp)                            :: time_forces
      INTEGER(int_8)                           :: ram_counter
   END TYPE

! **************************************************************************************************
   TYPE pair_list_element_type
      INTEGER, DIMENSION(2) :: pair
      INTEGER, DIMENSION(2) :: set_bounds
      INTEGER, DIMENSION(2) :: kind_pair
      REAL(KIND=dp)         :: r1(3), r2(3)
      REAL(KIND=dp)         :: dist2
   END TYPE

   ! **************************************************************************************************
   TYPE pair_set_list_type
      INTEGER, DIMENSION(2) :: pair
   END TYPE

! **************************************************************************************************
   TYPE pair_list_type
      TYPE(pair_list_element_type), DIMENSION(max_atom_block**2) :: elements
      INTEGER :: n_element
   END TYPE pair_list_type

! **************************************************************************************************
   TYPE hfx_cache_type
      INTEGER(int_8), DIMENSION(CACHE_SIZE)    :: DATA
      INTEGER                                  :: element_counter
   END TYPE

! **************************************************************************************************
   TYPE hfx_container_node
      TYPE(hfx_container_node), POINTER        :: next, prev
      INTEGER(int_8), DIMENSION(CACHE_SIZE)    :: DATA
   END TYPE

! **************************************************************************************************
   TYPE hfx_container_type
      TYPE(hfx_container_node), POINTER        :: first, current
      INTEGER                                  :: element_counter
      INTEGER(int_8)                           :: file_counter
      CHARACTER(LEN=5)                         :: desc
      INTEGER                                  :: unit
      CHARACTER(default_path_length)           :: filename
   END TYPE

! **************************************************************************************************
   TYPE hfx_basis_type
      INTEGER, DIMENSION(:), POINTER           :: lmax
      INTEGER, DIMENSION(:), POINTER           :: lmin
      INTEGER, DIMENSION(:), POINTER           :: npgf
      INTEGER                                  :: nset
      REAL(dp), DIMENSION(:, :), POINTER        :: zet
      INTEGER, DIMENSION(:), POINTER           :: nsgf
      INTEGER, DIMENSION(:, :), POINTER         :: first_sgf
      REAL(dp), DIMENSION(:, :), POINTER        :: sphi
      INTEGER                                  :: nsgf_total
      INTEGER, DIMENSION(:, :), POINTER         :: nl
      INTEGER, DIMENSION(:, :), POINTER         :: nsgfl
      INTEGER, DIMENSION(:), POINTER           :: nshell
      REAL(dp), DIMENSION(:, :, :, :), POINTER &
         :: sphi_ext
      REAL(dp), DIMENSION(:), POINTER          :: set_radius
      REAL(dp), DIMENSION(:, :), POINTER        :: pgf_radius
      REAL(dp)                                 :: kind_radius
   END TYPE

! **************************************************************************************************
   TYPE hfx_basis_info_type
      INTEGER                                  :: max_set
      INTEGER                                  :: max_sgf
      INTEGER                                  :: max_am
   END TYPE

! **************************************************************************************************
   TYPE hfx_screen_coeff_type
      REAL(dp)                                 :: x(2)
   END TYPE

! **************************************************************************************************
   TYPE hfx_p_kind
      REAL(dp), DIMENSION(:, :, :, :), POINTER    :: p_kind
   END TYPE

! **************************************************************************************************
   TYPE hfx_2D_map
      INTEGER, DIMENSION(:), POINTER           :: iatom_list
      INTEGER, DIMENSION(:), POINTER           :: jatom_list
   END TYPE

! **************************************************************************************************
   TYPE hfx_pgf_image
      REAL(dp)                                 :: ra(3), rb(3)
      REAL(dp)                                 :: rab2
      REAL(dp)                                 :: S1234
      REAL(dp)                                 :: P(3)
      REAL(dp)                                 :: R
      REAL(dp)                                 :: pgf_max
      REAL(dp), DIMENSION(3)                   :: bcell
   END TYPE

! **************************************************************************************************
   TYPE hfx_pgf_list
      TYPE(hfx_pgf_image), DIMENSION(:), POINTER &
         :: image_list
      INTEGER                                  :: nimages
      REAL(dp)                                 :: zetapzetb
      REAL(dp)                                 :: ZetaInv
      REAL(dp)                                 :: zeta, zetb
      INTEGER                                  :: ipgf, jpgf
   END TYPE

! **************************************************************************************************
   TYPE hfx_pgf_product_list
      REAL(dp)                                 :: ra(3), rb(3), rc(3), rd(3)
      REAL(dp)                                 :: ZetapEtaInv
      REAL(dp)                                 :: Rho, RhoInv
      REAL(dp)                                 :: P(3), Q(3), W(3)
      REAL(dp)                                 :: AB(3), CD(3)
      REAL(dp)                                 :: Fm(prim_data_f_size)
   END TYPE

! **************************************************************************************************
   TYPE hfx_block_range_type
      INTEGER        :: istart, iend
      INTEGER(int_8) :: cost
   END TYPE

! **************************************************************************************************
   TYPE hfx_task_list_type
      INTEGER                                  :: thread_id
      INTEGER                                  :: bin_id
      INTEGER(int_8)                           :: cost
   END TYPE

   TYPE :: hfx_compression_type
      TYPE(hfx_container_type), DIMENSION(:), &
         POINTER        :: maxval_container
      TYPE(hfx_cache_type), DIMENSION(:), &
         POINTER            :: maxval_cache
      TYPE(hfx_container_type), DIMENSION(:, :), &
         POINTER        :: integral_containers
      TYPE(hfx_cache_type), DIMENSION(:, :), &
         POINTER            :: integral_caches
      TYPE(hfx_container_type), POINTER :: maxval_container_disk
      TYPE(hfx_cache_type)     :: maxval_cache_disk
      TYPE(hfx_cache_type)      :: integral_caches_disk(64)
      TYPE(hfx_container_type), POINTER, &
         DIMENSION(:)  :: integral_containers_disk
   END TYPE

   TYPE :: block_ind_type
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: ind
   END TYPE

   TYPE hfx_ri_type
      ! input parameters (see input_cp2k_hfx)
      REAL(KIND=dp) :: filter_eps, filter_eps_2c, filter_eps_storage, filter_eps_mo, &
                       eps_lanczos, eps_pgf_orb, eps_eigval
      INTEGER :: t2c_sqrt_order, max_iter_lanczos, flavor, unit_nr_dbcsr, unit_nr, &
                 min_bsize, min_bsize_MO, t2c_method, nelectron_total
      LOGICAL :: check_2c_inv, calc_condnum

      TYPE(libint_potential_type) :: ri_metric

      ! input parameters from hfx
      TYPE(libint_potential_type) :: hfx_pot ! interaction potential
      REAL(KIND=dp) :: eps_schwarz ! integral screening threshold

      LOGICAL :: same_op ! whether RI operator is same as HF potential

      ! default process grid used for 3c tensors
      TYPE(dbcsr_t_pgrid_type), POINTER :: pgrid => NULL()
      TYPE(dbcsr_t_pgrid_type), POINTER :: pgrid_2d => NULL()

      ! distributions for (RI | AO AO) 3c integral tensor (non split)
      TYPE(distribution_3d_type) :: dist_3d
      TYPE(dbcsr_t_distribution_type) :: dist

      ! MO localization
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env => NULL()
      TYPE(section_vals_type), POINTER                   :: loc_subsection => NULL(), print_loc_subsection => NULL()
      LOGICAL                                            :: do_loc

      INTEGER, DIMENSION(:), ALLOCATABLE :: dist1_AO_1, dist1_AO_2, dist1_RI
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist2_AO_1, dist2_AO_2, dist2_RI
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist3_AO_1, dist3_AO_2, dist3_RI
      INTEGER, DIMENSION(:), ALLOCATABLE :: dist4_AO_1, dist4_AO_2, dist4_RI

      ! block sizes for RI and AO tensor dimensions (split)
      INTEGER, DIMENSION(:), ALLOCATABLE :: bsizes_RI, bsizes_AO, bsizes_RI_split, bsizes_AO_split, &
                                            bsizes_RI_fit, bsizes_AO_fit

      ! Note: changed static DIMENSION(1,1) of dbcsr_t_type to allocatables as workaround for gfortran 8.3.0,
      ! with static dimension gfortran gets stuck during compilation

      ! 2c tensor in (RI | RI) format for contraction
      TYPE(dbcsr_t_type), DIMENSION(:, :), ALLOCATABLE    :: t_2c_int

      ! 3c integral tensor in default format
      TYPE(dbcsr_t_type), DIMENSION(:, :), ALLOCATABLE :: t_3c_int

      ! 3c integral tensor in (AO RI | AO) format for contraction
      TYPE(dbcsr_t_type), DIMENSION(:, :), ALLOCATABLE :: t_3c_int_ctr_1
      TYPE(block_ind_type), DIMENSION(:, :), ALLOCATABLE :: blk_indices
      TYPE(dbcsr_t_pgrid_type), POINTER                :: pgrid_1 => NULL()

      ! 3c integral tensor in ( AO | RI AO) (MO) or (AO RI | AO) (RHO) format for contraction
      TYPE(dbcsr_t_type), DIMENSION(:, :), ALLOCATABLE :: t_3c_int_ctr_2
      TYPE(dbcsr_t_pgrid_type), POINTER                :: pgrid_2 => NULL()

      ! 3c integral tensor in ( RI | AO AO ) format for contraction
      TYPE(dbcsr_t_type), DIMENSION(:, :), ALLOCATABLE :: t_3c_int_ctr_3
      TYPE(dbcsr_t_pgrid_type), POINTER                :: pgrid_3 => NULL()

      ! 3c integral tensor in (RI | MO AO ) format for contraction
      TYPE(dbcsr_t_type), DIMENSION(:, :, :), ALLOCATABLE :: t_3c_int_mo
      TYPE(dbcsr_t_type), DIMENSION(:, :, :), ALLOCATABLE :: t_3c_ctr_RI
      TYPE(dbcsr_t_type), DIMENSION(:, :, :), ALLOCATABLE :: t_3c_ctr_KS
      TYPE(dbcsr_t_type), DIMENSION(:, :, :), ALLOCATABLE :: t_3c_ctr_KS_copy
      ! 3c integral tensor in ( RI AO | AO ) format for contraction
      !TYPE(dbcsr_t_type), DIMENSION(:, :), ALLOCATABLE :: t_3c_int_ctr_4
      !TYPE(dbcsr_t_pgrid_type), POINTER                :: pgrid_4 => NULL()

      ! optional: sections for output handling
      ! alternatively set unit_nr_dbcsr (for logging tensor operations) and unit_nr (for general
      ! output) directly
      TYPE(section_vals_type), POINTER :: ri_section => NULL(), hfx_section => NULL()

      ! types of primary and auxiliary basis
      CHARACTER(len=default_string_length) :: orb_basis_type, ri_basis_type

      ! memory reduction factor
      INTEGER :: n_mem

      ! offsets for memory batches
      INTEGER, DIMENSION(:), ALLOCATABLE :: starts_array_mem_block, ends_array_mem_block
      INTEGER, DIMENSION(:), ALLOCATABLE :: starts_array_mem, ends_array_mem

      INTEGER, DIMENSION(:), ALLOCATABLE :: starts_array_RI_mem_block, ends_array_RI_mem_block
      INTEGER, DIMENSION(:), ALLOCATABLE :: starts_array_RI_mem, ends_array_RI_mem

      INTEGER(int_8) :: dbcsr_nflop
      REAL(dp)       :: dbcsr_time
      INTEGER        :: num_pe
      TYPE(hfx_compression_type), DIMENSION(:, :), ALLOCATABLE :: store_3c

   END TYPE

! **************************************************************************************************
!> \brief stores some data used in construction of Kohn-Sham matrix
!> \param potential_parameter stores information on the potential (1/r, erfc(wr)/r
!> \param screening_parameter stores screening infos such as epsilon
!> \param memory_parameter stores infos on memory used for in-core calculations
!> \param periodic_parameter stores information on how to apply pbc
!> \param load_balance_parameter contains infos for Monte Carlo simulated annealing
!> \param general_paramter at the moment stores the fraction of HF amount to be included
!> \param maxval_container stores the maxvals in compressed form
!> \param maxval_cache cache for maxvals in decompressed form
!> \param integral_containers 64 containers for compressed integrals
!> \param integral_caches 64 caches for decompressed integrals
!> \param neighbor_cells manages handling of periodic cells
!> \param distribution_energy stores information on parallelization of energy
!> \param distribution_forces stores information on parallelization of forces
!> \param initial_p stores the initial guess if requested
!> \param is_assoc_atomic_block reflects KS sparsity
!> \param number_of_p_entries Size of P matrix
!> \param n_rep_hf Number of HFX replicas
!> \param b_first_load_balance_x flag to indicate if it is enough just to update
!>        the distribution of the integrals
!> \param full_ks_x full ks matrices
!> \param lib libint type for eris
!> \param basis_info contains information for basis sets
!> \param screen_funct_coeffs_pgf pgf based near field screening coefficients
!> \param pair_dist_radii_pgf pgf based radii coefficients of pair distributions
!> \param screen_funct_coeffs_set set based near field screening coefficients
!> \param screen_funct_coeffs_kind kind based near field screening coefficients
!> \param screen_funct_is_initialized flag that indicates if the coefficients
!>        have already been fitted
!> \par History
!>      11.2006 created [Manuel Guidon]
!>      02.2009 completely rewritten due to new screening
!> \author Manuel Guidon
! **************************************************************************************************
   TYPE hfx_type
      TYPE(hfx_potential_type)                 :: potential_parameter
      TYPE(hfx_screening_type)                 :: screening_parameter
      TYPE(hfx_memory_type)                    :: memory_parameter
      TYPE(hfx_periodic_type)                  :: periodic_parameter
      TYPE(hfx_load_balance_type)              :: load_balance_parameter
      TYPE(hfx_general_type)                   :: general_parameter

      TYPE(hfx_compression_type) :: store_ints
      TYPE(hfx_compression_type) :: store_forces

      TYPE(hfx_cell_type), DIMENSION(:), &
         POINTER                       :: neighbor_cells
      TYPE(hfx_distribution), DIMENSION(:), &
         POINTER         :: distribution_energy
      TYPE(hfx_distribution), DIMENSION(:), &
         POINTER         :: distribution_forces
      INTEGER, DIMENSION(:, :), POINTER         :: is_assoc_atomic_block
      INTEGER                                  :: number_of_p_entries
      TYPE(hfx_basis_type), DIMENSION(:), &
         POINTER           :: basis_parameter
      INTEGER                                  :: n_rep_hf
      LOGICAL                                  :: b_first_load_balance_energy, &
                                                  b_first_load_balance_forces
      REAL(dp), DIMENSION(:, :), POINTER        :: full_ks_alpha
      REAL(dp), DIMENSION(:, :), POINTER        :: full_ks_beta
      TYPE(cp_libint_t)                        :: lib
      TYPE(hfx_basis_info_type)                :: basis_info
      TYPE(hfx_screen_coeff_type), &
         DIMENSION(:, :, :, :, :, :), POINTER     :: screen_funct_coeffs_pgf, &
                                                     pair_dist_radii_pgf
      TYPE(hfx_screen_coeff_type), &
         DIMENSION(:, :, :, :), POINTER         :: screen_funct_coeffs_set
      TYPE(hfx_screen_coeff_type), &
         DIMENSION(:, :), POINTER             :: screen_funct_coeffs_kind
      LOGICAL                                  :: screen_funct_is_initialized
      TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: initial_p
      TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: initial_p_forces
      INTEGER, DIMENSION(:), POINTER           :: map_atom_to_kind_atom
      TYPE(hfx_2D_map), DIMENSION(:), POINTER  :: map_atoms_to_cpus
      INTEGER, DIMENSION(:, :), POINTER         :: atomic_block_offset
      INTEGER, DIMENSION(:, :, :, :), POINTER     :: set_offset
      INTEGER, DIMENSION(:), POINTER           :: block_offset
      TYPE(hfx_block_range_type), DIMENSION(:), &
         POINTER      :: blocks
      TYPE(hfx_task_list_type), DIMENSION(:), &
         POINTER        :: task_list
      REAL(dp), DIMENSION(:, :), POINTER        :: pmax_atom, pmax_atom_forces
      TYPE(cp_libint_t)                         :: lib_deriv
      REAL(dp), DIMENSION(:, :), POINTER        :: pmax_block
      LOGICAL, DIMENSION(:, :), POINTER         :: atomic_pair_list
      LOGICAL, DIMENSION(:, :), POINTER         :: atomic_pair_list_forces
      LOGICAL                                   :: do_hfx_ri
      TYPE(hfx_ri_type), POINTER                :: ri_data
   END TYPE hfx_type

CONTAINS

! **************************************************************************************************
!> \brief - This routine allocates and initializes all types in hfx_data
!> \param x_data contains all relevant data structures for hfx runs
!> \param para_env ...
!> \param hfx_section input section
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param dft_control ...
!> \param cell ...
!> \param do_exx ...
!> \param do_ot ...
!> \param nelectron_total ...
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      - All POINTERS and ALLOCATABLES are allocated, even if their size is
!>        unknown at invocation time
! **************************************************************************************************
   SUBROUTINE hfx_create(x_data, para_env, hfx_section, atomic_kind_set, qs_kind_set, &
                         particle_set, dft_control, cell, do_exx, do_ot, nelectron_total)
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(cp_para_env_type)                             :: para_env
      TYPE(section_vals_type), POINTER                   :: hfx_section
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL, OPTIONAL                                  :: do_exx, do_ot
      INTEGER, OPTIONAL                                  :: nelectron_total

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create', routineP = moduleN//':'//routineN

      CHARACTER(LEN=512)                                 :: error_msg
      CHARACTER(LEN=default_path_length)                 :: char_val
      INTEGER :: handle, i, i_thread, iatom, ikind, int_val, irep, jkind, max_set, n_rep_hf, &
         n_threads, natom, natom_a, natom_b, nkind, nseta, nsetb, pbc_shells, storage_id
      INTEGER, DIMENSION(:), POINTER                     :: atom2kind, kind_of
      LOGICAL                                            :: do_ri, logic_val, my_do_exx
      REAL(dp)                                           :: real_val
      TYPE(hfx_type), POINTER                            :: actual_x_data
      TYPE(section_vals_type), POINTER                   :: hf_pbc_section, hf_sub_section, &
                                                            hfx_ri_section

      CALL timeset(routineN, handle)

      CALL cite_reference(Guidon2008)
      CALL cite_reference(Guidon2009)

      natom = SIZE(particle_set)

      !! There might be 2 hf sections
      CALL section_vals_get(hfx_section, n_repetition=n_rep_hf)
      n_threads = 1
!$    n_threads = omp_get_max_threads()

      CALL section_vals_val_get(hfx_section, "RI%_SECTION_PARAMETERS_", l_val=do_ri)
      IF (do_ri) n_threads = 1 ! RI implementation does not use threads

      ALLOCATE (x_data(n_rep_hf, n_threads))
      DO i_thread = 1, n_threads
         DO irep = 1, n_rep_hf
            actual_x_data => x_data(irep, i_thread)
            !! Get data from input file
            !!
            !! GENERAL params
            CALL section_vals_val_get(hfx_section, "FRACTION", r_val=real_val, i_rep_section=irep)
            actual_x_data%general_parameter%fraction = real_val
            actual_x_data%n_rep_hf = n_rep_hf

            NULLIFY (actual_x_data%map_atoms_to_cpus)

            CALL section_vals_val_get(hfx_section, "TREAT_LSD_IN_CORE", l_val=logic_val, i_rep_section=irep)
            actual_x_data%general_parameter%treat_lsd_in_core = logic_val

            hfx_ri_section => section_vals_get_subs_vals(hfx_section, "RI")
            CALL section_vals_val_get(hfx_ri_section, "_SECTION_PARAMETERS_", l_val=actual_x_data%do_hfx_ri)

            !! MEMORY section
            hf_sub_section => section_vals_get_subs_vals(hfx_section, "MEMORY", i_rep_section=irep)
            CALL parse_memory_section(actual_x_data%memory_parameter, hf_sub_section, storage_id, i_thread, &
                                      n_threads, para_env, irep, skip_disk=.FALSE., skip_in_core_forces=.FALSE.)

            !! PERIODIC section
            hf_sub_section => section_vals_get_subs_vals(hfx_section, "PERIODIC", i_rep_section=irep)
            CALL section_vals_val_get(hf_sub_section, "NUMBER_OF_SHELLS", i_val=int_val)
            actual_x_data%periodic_parameter%number_of_shells = int_val
            actual_x_data%periodic_parameter%mode = int_val
            CALL get_cell(cell=cell, periodic=actual_x_data%periodic_parameter%perd)
            IF (SUM(actual_x_data%periodic_parameter%perd) == 0) THEN
               actual_x_data%periodic_parameter%do_periodic = .FALSE.
            ELSE
               actual_x_data%periodic_parameter%do_periodic = .TRUE.
            END IF

            !! SCREENING section
            hf_sub_section => section_vals_get_subs_vals(hfx_section, "SCREENING", i_rep_section=irep)
            CALL section_vals_val_get(hf_sub_section, "EPS_SCHWARZ", r_val=real_val)
            actual_x_data%screening_parameter%eps_schwarz = real_val
            CALL section_vals_val_get(hf_sub_section, "EPS_SCHWARZ_FORCES", r_val=real_val)
            actual_x_data%screening_parameter%eps_schwarz_forces = real_val
            CALL section_vals_val_get(hf_sub_section, "SCREEN_P_FORCES", l_val=logic_val)
            actual_x_data%screening_parameter%do_p_screening_forces = logic_val
            CALL section_vals_val_get(hf_sub_section, "SCREEN_ON_INITIAL_P", l_val=logic_val)
            actual_x_data%screening_parameter%do_initial_p_screening = logic_val
            actual_x_data%screen_funct_is_initialized = .FALSE.

            !! INTERACTION_POTENTIAL section
            hf_sub_section => section_vals_get_subs_vals(hfx_section, "INTERACTION_POTENTIAL", i_rep_section=irep)
            CALL section_vals_val_get(hf_sub_section, "POTENTIAL_TYPE", i_val=int_val)
            actual_x_data%potential_parameter%potential_type = int_val
            CALL section_vals_val_get(hf_sub_section, "OMEGA", r_val=real_val)
            actual_x_data%potential_parameter%omega = real_val
            CALL section_vals_val_get(hf_sub_section, "SCALE_COULOMB", r_val=real_val)
            actual_x_data%potential_parameter%scale_coulomb = real_val
            CALL section_vals_val_get(hf_sub_section, "SCALE_LONGRANGE", r_val=real_val)
            actual_x_data%potential_parameter%scale_longrange = real_val
            CALL section_vals_val_get(hf_sub_section, "SCALE_GAUSSIAN", r_val=real_val)
            actual_x_data%potential_parameter%scale_gaussian = real_val
            IF (actual_x_data%potential_parameter%potential_type == do_potential_truncated .OR. &
                actual_x_data%potential_parameter%potential_type == do_potential_mix_cl_trunc) THEN
               CALL section_vals_val_get(hf_sub_section, "CUTOFF_RADIUS", r_val=real_val)
               actual_x_data%potential_parameter%cutoff_radius = real_val
               CALL section_vals_val_get(hf_sub_section, "T_C_G_DATA", c_val=char_val)
               CALL compress(char_val, .TRUE.)
               ! ** Check if file is there
               IF (.NOT. file_exists(char_val)) THEN
                  WRITE (error_msg, '(A,A,A)') "Truncated hfx calculation requested. The file containing "// &
                     "the data could not be found at ", TRIM(char_val), " Please check T_C_G_DATA "// &
                     "in the INTERACTION_POTENTIAL section"
                  CPABORT(error_msg)
               ELSE
                  actual_x_data%potential_parameter%filename = char_val
               END IF
            END IF
            IF (actual_x_data%potential_parameter%potential_type == do_potential_short) THEN
               CALL erfc_cutoff(actual_x_data%screening_parameter%eps_schwarz, &
                                actual_x_data%potential_parameter%omega, &
                                actual_x_data%potential_parameter%cutoff_radius)
            END IF

            !! LOAD_BALANCE section
            hf_sub_section => section_vals_get_subs_vals(hfx_section, "LOAD_BALANCE", i_rep_section=irep)
            CALL section_vals_val_get(hf_sub_section, "NBINS", i_val=int_val)
            actual_x_data%load_balance_parameter%nbins = MAX(1, int_val)
            actual_x_data%load_balance_parameter%blocks_initialized = .FALSE.

            CALL section_vals_val_get(hf_sub_section, "RANDOMIZE", l_val=logic_val)
            actual_x_data%load_balance_parameter%do_randomize = logic_val

            actual_x_data%load_balance_parameter%rtp_redistribute = .FALSE.
            IF (ASSOCIATED(dft_control%rtp_control)) &
               actual_x_data%load_balance_parameter%rtp_redistribute = dft_control%rtp_control%hfx_redistribute

            CALL section_vals_val_get(hf_sub_section, "BLOCK_SIZE", i_val=int_val)
            ! negative values ask for a computed default
            IF (int_val <= 0) THEN
               ! this gives a reasonable number of blocks for binning, yet typically results in blocking.
               int_val = CEILING(0.1_dp*natom/ &
                                 REAL(actual_x_data%load_balance_parameter%nbins*n_threads*para_env%num_pe, KIND=dp)**(0.25_dp))
            ENDIF
            ! at least 1 atom per block, and avoid overly large blocks
            actual_x_data%load_balance_parameter%block_size = MIN(max_atom_block, MAX(1, int_val))

            my_do_exx = .FALSE.
            IF (PRESENT(do_exx)) THEN
               my_do_exx = do_exx
            END IF

            ! in the case of RPA exchange (postprocessing after SCF), we do not use ADMM
            IF (my_do_exx) THEN
               CALL hfx_create_basis_types(actual_x_data%basis_parameter, actual_x_data%basis_info, qs_kind_set, &
                                           .FALSE.)
            ELSE
               CALL hfx_create_basis_types(actual_x_data%basis_parameter, actual_x_data%basis_info, qs_kind_set, &
                                           dft_control%do_admm)
            END IF

!!**************************************************************************************************
!! **        !! ** This code writes the contraction routines
!! **        !! ** Very UGLY: BASIS_SET has to be 1 primitive and lmin=lmax=l. For g-functions
!! **        !! **
!! **        !! ** 1  4  4  1  1
!! **        !! **    1.0  1.0
!! **        !! **
!! **        k = max_am - 1
!! **        write(filename,'(A,I0,A)') "sphi",k+1,"a"
!! **        OPEN(UNIT=31415,FILE=filename)
!! **        DO i=ncoset(k)+1,SIZE(sphi_a,1)
!! **          DO j=1,SIZE(sphi_a,2)
!! **            IF( sphi_a(i,j) /= 0.0_dp) THEN
!! **              write(31415,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') "buffer1(i+imax*(",&
!! **                          j,&
!! **                          "-1)) = buffer1(i+imax*(",&
!! **                          j,&
!! **                          "-1)) + work(",&
!! **                          i-ncoset(k),&
!! **                          "+(i-1)*kmax) * sphi_a(",&
!! **                          i-ncoset(k),&
!! **                          ",",&
!! **                          j,&
!! **                          "+s_offset_a1)"
!! **            END IF
!! **          END DO
!! **        END DO
!! **        CLOSE(UNIT=31415)
!! **        write(filename,'(A,I0,A)') "sphi",k+1,"b"
!! **        OPEN(UNIT=31415,FILE=filename)
!! **        DO i=ncoset(k)+1,SIZE(sphi_a,1)
!! **          DO j=1,SIZE(sphi_a,2)
!! **            IF( sphi_a(i,j) /= 0.0_dp) THEN
!! **               write(31415,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') "buffer2(i+imax*(",&
!! **                          j,&
!! **                          "-1)) = buffer2(i+imax*(",&
!! **                          j,&
!! **                          "-1)) + buffer1(",&
!! **                          i-ncoset(k),&
!! **                          "+(i-1)*kmax) * sphi_b(",&
!! **                          i-ncoset(k),&
!! **                          ",",&
!! **                          j,&
!! **                          "+s_offset_b1)"
!! **
!! **            END IF
!! **          END DO
!! **        END DO
!! **        CLOSE(UNIT=31415)
!! **        write(filename,'(A,I0,A)') "sphi",k+1,"c"
!! **        OPEN(UNIT=31415,FILE=filename)
!! **        DO i=ncoset(k)+1,SIZE(sphi_a,1)
!! **          DO j=1,SIZE(sphi_a,2)
!! **            IF( sphi_a(i,j) /= 0.0_dp) THEN
!! **               write(31415,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') "buffer1(i+imax*(",&
!! **                          j,&
!! **                          "-1)) = buffer1(i+imax*(",&
!! **                          j,&
!! **                          "-1)) + buffer2(",&
!! **                          i-ncoset(k),&
!! **                          "+(i-1)*kmax) * sphi_c(",&
!! **                          i-ncoset(k),&
!! **                          ",",&
!! **                          j,&
!! **                          "+s_offset_c1)"
!! **
!! **            END IF
!! **          END DO
!! **        END DO
!! **        CLOSE(UNIT=31415)
!! **        write(filename,'(A,I0,A)') "sphi",k+1,"d"
!! **        OPEN(UNIT=31415,FILE=filename)
!! **        DO i=ncoset(k)+1,SIZE(sphi_a,1)
!! **          DO j=1,SIZE(sphi_a,2)
!! **            IF( sphi_a(i,j) /= 0.0_dp) THEN
!! **
!! **
!! **               write(31415,'(A,I0,A)') "primitives(s_offset_a1+i3, s_offset_b1+i2, s_offset_c1+i1, s_offset_d1+",&
!! **                           j,")= &"
!! **               write(31415,'(A,I0,A)') "primitives(s_offset_a1+i3, s_offset_b1+i2, s_offset_c1+i1, s_offset_d1+",&
!! **                           j,")+ &"
!! **               write(31415,'(A,I0,A,I0,A,I0,A)') "buffer1(",&
!! **                          i-ncoset(k),&
!! **                          "+(i-1)*kmax) * sphi_d(",&
!! **                          i-ncoset(k),&
!! **                          ",",&
!! **                          j,&
!! **                          "+s_offset_d1)"
!! **
!! **
!! **            END IF
!! **          END DO
!! **        END DO
!! **        CLOSE(UNIT=31415)
!! **        stop
!! *************************************************************************************************************************

            IF (actual_x_data%periodic_parameter%do_periodic) THEN
               hf_pbc_section => section_vals_get_subs_vals(hfx_section, "PERIODIC", i_rep_section=irep)
               CALL section_vals_val_get(hf_pbc_section, "NUMBER_OF_SHELLS", i_val=pbc_shells)
               actual_x_data%periodic_parameter%number_of_shells_from_input = pbc_shells
               ALLOCATE (actual_x_data%neighbor_cells(1))
               CALL hfx_create_neighbor_cells(actual_x_data, pbc_shells, cell, i_thread)
            ELSE
               ALLOCATE (actual_x_data%neighbor_cells(1))
               ! ** Initialize this guy to enable non periodic stress regtests
               actual_x_data%periodic_parameter%R_max_stress = 1.0_dp
            END IF

            nkind = SIZE(qs_kind_set, 1)
            max_set = actual_x_data%basis_info%max_set

            !! ** This guy is allocated on the master thread only
            IF (i_thread == 1) THEN
               ALLOCATE (actual_x_data%is_assoc_atomic_block(natom, natom))
               ALLOCATE (actual_x_data%atomic_block_offset(natom, natom))
               ALLOCATE (actual_x_data%set_offset(max_set, max_set, nkind, nkind))
               ALLOCATE (actual_x_data%block_offset(para_env%num_pe + 1))
            END IF

            ALLOCATE (actual_x_data%distribution_forces(1))
            ALLOCATE (actual_x_data%distribution_energy(1))

            actual_x_data%memory_parameter%size_p_screen = 0_int_8
            IF (i_thread == 1) THEN
               ALLOCATE (actual_x_data%atomic_pair_list(natom, natom))
               ALLOCATE (actual_x_data%atomic_pair_list_forces(natom, natom))
            END IF

            IF (actual_x_data%screening_parameter%do_initial_p_screening .OR. &
                actual_x_data%screening_parameter%do_p_screening_forces) THEN
               !! ** This guy is allocated on the master thread only
               IF (i_thread == 1) THEN
                  ALLOCATE (actual_x_data%pmax_atom(natom, natom))
                  ALLOCATE (actual_x_data%initial_p(nkind*(nkind + 1)/2))
                  i = 1
                  DO ikind = 1, nkind
                     CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_a)
                     nseta = actual_x_data%basis_parameter(ikind)%nset
                     DO jkind = ikind, nkind
                        CALL get_atomic_kind(atomic_kind_set(jkind), natom=natom_b)
                        nsetb = actual_x_data%basis_parameter(jkind)%nset
                        ALLOCATE (actual_x_data%initial_p(i)%p_kind(nseta, nsetb, natom_a, natom_b))
                        actual_x_data%memory_parameter%size_p_screen = &
                           actual_x_data%memory_parameter%size_p_screen + nseta*nsetb*natom_a*natom_b
                        i = i + 1
                     END DO
                  END DO
                  IF (actual_x_data%memory_parameter%treat_forces_in_core) THEN
                     ALLOCATE (actual_x_data%pmax_atom_forces(natom, natom))

                     ALLOCATE (actual_x_data%initial_p_forces(nkind*(nkind + 1)/2))
                     i = 1
                     DO ikind = 1, nkind
                        CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_a)
                        nseta = actual_x_data%basis_parameter(ikind)%nset
                        DO jkind = ikind, nkind
                           CALL get_atomic_kind(atomic_kind_set(jkind), natom=natom_b)
                           nsetb = actual_x_data%basis_parameter(jkind)%nset
                           ALLOCATE (actual_x_data%initial_p_forces(i)%p_kind(nseta, nsetb, natom_a, natom_b))
                           actual_x_data%memory_parameter%size_p_screen = &
                              actual_x_data%memory_parameter%size_p_screen + nseta*nsetb*natom_a*natom_b
                           i = i + 1
                        END DO
                     END DO
                  END IF
               END IF
               ALLOCATE (actual_x_data%map_atom_to_kind_atom(natom))
               ALLOCATE (kind_of(natom))
               ALLOCATE (atom2kind(nkind))
               CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

               atom2kind = 0
               DO iatom = 1, natom
                  ikind = kind_of(iatom)
                  atom2kind(ikind) = atom2kind(ikind) + 1
                  actual_x_data%map_atom_to_kind_atom(iatom) = atom2kind(ikind)
               END DO
               DEALLOCATE (kind_of, atom2kind)
            END IF

            ! ** Initialize libint type
            CALL cp_libint_static_init()
            CALL cp_libint_init_eri(actual_x_data%lib, actual_x_data%basis_info%max_am)
            CALL cp_libint_init_eri1(actual_x_data%lib_deriv, actual_x_data%basis_info%max_am)
            CALL cp_libint_set_contrdepth(actual_x_data%lib, 1)
            CALL cp_libint_set_contrdepth(actual_x_data%lib_deriv, 1)

            CALL alloc_containers(actual_x_data%store_ints, 1)
            CALL alloc_containers(actual_x_data%store_forces, 1)

            actual_x_data%store_ints%maxval_cache_disk%element_counter = 1
            ALLOCATE (actual_x_data%store_ints%maxval_container_disk)
            ALLOCATE (actual_x_data%store_ints%maxval_container_disk%first)
            actual_x_data%store_ints%maxval_container_disk%first%prev => NULL()
            actual_x_data%store_ints%maxval_container_disk%first%next => NULL()
            actual_x_data%store_ints%maxval_container_disk%current => actual_x_data%store_ints%maxval_container_disk%first
            actual_x_data%store_ints%maxval_container_disk%current%data = 0
            actual_x_data%store_ints%maxval_container_disk%element_counter = 1
            actual_x_data%store_ints%maxval_container_disk%file_counter = 1
            actual_x_data%store_ints%maxval_container_disk%desc = 'Max_'
            actual_x_data%store_ints%maxval_container_disk%unit = -1
            WRITE (actual_x_data%store_ints%maxval_container_disk%filename, '(A,I0,A,A,A)') &
               TRIM(actual_x_data%memory_parameter%storage_location), &
               storage_id, "_", actual_x_data%store_ints%maxval_container_disk%desc, "6"
            CALL compress(actual_x_data%store_ints%maxval_container_disk%filename, .TRUE.)
            ALLOCATE (actual_x_data%store_ints%integral_containers_disk(64))
            DO i = 1, 64
               actual_x_data%store_ints%integral_caches_disk(i)%element_counter = 1
               actual_x_data%store_ints%integral_caches_disk(i)%data = 0
               ALLOCATE (actual_x_data%store_ints%integral_containers_disk(i)%first)
               actual_x_data%store_ints%integral_containers_disk(i)%first%prev => NULL()
               actual_x_data%store_ints%integral_containers_disk(i)%first%next => NULL()
               actual_x_data%store_ints%integral_containers_disk(i)%current => &
                  actual_x_data%store_ints%integral_containers_disk(i)%first
               actual_x_data%store_ints%integral_containers_disk(i)%current%data = 0
               actual_x_data%store_ints%integral_containers_disk(i)%element_counter = 1
               actual_x_data%store_ints%integral_containers_disk(i)%file_counter = 1
               actual_x_data%store_ints%integral_containers_disk(i)%desc = 'Int_'
               actual_x_data%store_ints%integral_containers_disk(i)%unit = -1
               WRITE (actual_x_data%store_ints%integral_containers_disk(i)%filename, '(A,I0,A,A,I0)') &
                  TRIM(actual_x_data%memory_parameter%storage_location), &
                  storage_id, "_", actual_x_data%store_ints%integral_containers_disk(i)%desc, i
               CALL compress(actual_x_data%store_ints%integral_containers_disk(i)%filename, .TRUE.)
            END DO

            actual_x_data%b_first_load_balance_energy = .TRUE.
            actual_x_data%b_first_load_balance_forces = .TRUE.

            hf_sub_section => section_vals_get_subs_vals(hfx_section, "RI", i_rep_section=irep)
            IF (actual_x_data%do_hfx_ri) THEN
               CPASSERT(PRESENT(do_ot))
               CPASSERT(PRESENT(nelectron_total))
               ALLOCATE (actual_x_data%ri_data)
               CALL hfx_ri_init_read_input_from_hfx(actual_x_data%ri_data, actual_x_data, hfx_section, &
                                                    hf_sub_section, qs_kind_set, &
                                                    particle_set, atomic_kind_set, dft_control, para_env, irep, do_ot, &
                                                    nelectron_total)

            ENDIF
         END DO
      END DO

      DO irep = 1, n_rep_hf
         actual_x_data => x_data(irep, 1)
         CALL hfx_print_info(actual_x_data, hfx_section, irep)
      END DO
      CALL timestop(handle)
   END SUBROUTINE hfx_create

! **************************************************************************************************
!> \brief Read RI input and initialize RI data for use within Hartree-Fock
!> \param ri_data ...
!> \param x_data ...
!> \param hfx_section ...
!> \param ri_section ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param dft_control ...
!> \param para_env ...
!> \param irep ...
!> \param do_ot ...
!> \param nelectron_total ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_init_read_input_from_hfx(ri_data, x_data, hfx_section, ri_section, qs_kind_set, &
                                              particle_set, atomic_kind_set, dft_control, para_env, irep, do_ot, nelectron_total)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(hfx_type), INTENT(INOUT)                      :: x_data
      TYPE(section_vals_type), POINTER                   :: hfx_section, ri_section
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_para_env_type)                             :: para_env
      INTEGER, INTENT(IN)                                :: irep
      LOGICAL, INTENT(IN)                                :: do_ot
      INTEGER, INTENT(IN)                                :: nelectron_total

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_init_read_input_from_hfx', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=512)                                 :: error_msg
      CHARACTER(LEN=default_path_length)                 :: char_val, t_c_filename
      CHARACTER(LEN=default_string_length)               :: orb_basis_type, ri_basis_type
      INTEGER                                            :: handle, unit_nr, unit_nr_dbcsr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: hf_sub_section

      CALL timeset(routineN, handle)

      NULLIFY (hf_sub_section)

      ASSOCIATE (hfx_pot=>ri_data%hfx_pot)
         hfx_pot%potential_type = x_data%potential_parameter%potential_type
         hfx_pot%omega = x_data%potential_parameter%omega
         hfx_pot%cutoff_radius = x_data%potential_parameter%cutoff_radius
      END ASSOCIATE
      ri_data%ri_section => ri_section
      ri_data%hfx_section => hfx_section
      ri_data%eps_schwarz = x_data%screening_parameter%eps_schwarz

      logger => cp_get_default_logger()
      unit_nr_dbcsr = cp_print_key_unit_nr(logger, ri_data%ri_section, "RI_INFO", &
                                           extension=".dbcsrLog")

      unit_nr = cp_print_key_unit_nr(logger, ri_data%hfx_section, "HF_INFO", &
                                     extension=".scfLog")

      hf_sub_section => section_vals_get_subs_vals(hfx_section, "INTERACTION_POTENTIAL", i_rep_section=irep)
      CALL section_vals_val_get(hf_sub_section, "T_C_G_DATA", c_val=char_val)
      CALL compress(char_val, .TRUE.)

      IF (.NOT. file_exists(char_val)) THEN
         WRITE (error_msg, '(A,A,A)') "File not found. Please check T_C_G_DATA "// &
            "in the INTERACTION_POTENTIAL section"
         CPABORT(error_msg)
      ELSE
         t_c_filename = char_val
      END IF

      IF (dft_control%do_admm) THEN
         orb_basis_type = "AUX_FIT"
      ELSE
         orb_basis_type = "ORB"
      ENDIF

      ri_basis_type = "RI_HFX"

      CALL hfx_ri_init_read_input(ri_data, ri_section, qs_kind_set, particle_set, atomic_kind_set, &
                                  orb_basis_type, ri_basis_type, dft_control%auto_basis_ri_hfx, &
                                  para_env, unit_nr, unit_nr_dbcsr, nelectron_total, &
                                  t_c_filename=t_c_filename)

      IF (ri_data%flavor == ri_mo .AND. .NOT. do_ot) THEN
         CPABORT("RI_FLAVOR MO is not consistent with density mixing. Either use OT or RI_FLAVOR RHO.")
      ENDIF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief General routine for reading input of RI section and initializing RI data
!> \param ri_data ...
!> \param ri_section ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param orb_basis_type ...
!> \param ri_basis_type ...
!> \param auto_basis_cntrl ...
!> \param para_env ...
!> \param unit_nr unit number of general output
!> \param unit_nr_dbcsr unit number for logging DBCSR tensor operations
!> \param nelectron_total ...
!> \param t_c_filename ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_init_read_input(ri_data, ri_section, qs_kind_set, &
                                     particle_set, atomic_kind_set, orb_basis_type, ri_basis_type, auto_basis_cntrl, para_env, &
                                     unit_nr, unit_nr_dbcsr, nelectron_total, t_c_filename)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(section_vals_type), POINTER                   :: ri_section
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      CHARACTER(LEN=*), INTENT(IN)                       :: orb_basis_type, ri_basis_type
      INTEGER, INTENT(IN)                                :: auto_basis_cntrl
      TYPE(cp_para_env_type)                             :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr, unit_nr_dbcsr, nelectron_total
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: t_c_filename

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_init_read_input', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: explicit
      REAL(dp)                                           :: eps_storage_scaling

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(ri_section, "EPS_FILTER", r_val=ri_data%filter_eps)
      CALL section_vals_val_get(ri_section, "EPS_FILTER_2C", r_val=ri_data%filter_eps_2c)
      CALL section_vals_val_get(ri_section, "EPS_STORAGE_SCALING", r_val=eps_storage_scaling)
      ri_data%filter_eps_storage = ri_data%filter_eps*eps_storage_scaling
      CALL section_vals_val_get(ri_section, "EPS_FILTER_MO", r_val=ri_data%filter_eps_mo)

      ASSOCIATE (ri_metric=>ri_data%ri_metric, hfx_pot=>ri_data%hfx_pot)
         CALL section_vals_val_get(ri_section, "RI_METRIC", i_val=ri_metric%potential_type, explicit=explicit)
         IF (.NOT. explicit .OR. ri_metric%potential_type == 0) THEN
            ri_metric%potential_type = hfx_pot%potential_type
         ENDIF

         CALL section_vals_val_get(ri_section, "OMEGA", r_val=ri_metric%omega, explicit=explicit)
         IF (.NOT. explicit) THEN
            ri_metric%omega = hfx_pot%omega
         ENDIF

         CALL section_vals_val_get(ri_section, "CUTOFF_RADIUS", r_val=ri_metric%cutoff_radius, explicit=explicit)
         IF (.NOT. explicit) THEN
            ri_metric%cutoff_radius = hfx_pot%cutoff_radius
         ENDIF
      END ASSOCIATE

      CALL section_vals_val_get(ri_section, "2C_MATRIX_FUNCTIONS", i_val=ri_data%t2c_method)
      CALL section_vals_val_get(ri_section, "EPS_EIGVAL", r_val=ri_data%eps_eigval)
      CALL section_vals_val_get(ri_section, "CHECK_2C_MATRIX", l_val=ri_data%check_2c_inv)
      CALL section_vals_val_get(ri_section, "CALC_COND_NUM", l_val=ri_data%calc_condnum)
      CALL section_vals_val_get(ri_section, "SQRT_ORDER", i_val=ri_data%t2c_sqrt_order)
      CALL section_vals_val_get(ri_section, "EPS_LANCZOS", r_val=ri_data%eps_lanczos)
      CALL section_vals_val_get(ri_section, "MAX_ITER_LANCZOS", i_val=ri_data%max_iter_lanczos)
      CALL section_vals_val_get(ri_section, "RI_FLAVOR", i_val=ri_data%flavor)
      CALL section_vals_val_get(ri_section, "EPS_PGF_ORB", r_val=ri_data%eps_pgf_orb)
      CALL section_vals_val_get(ri_section, "MIN_BLOCK_SIZE", i_val=ri_data%min_bsize)
      CALL section_vals_val_get(ri_section, "MIN_BLOCK_SIZE_MO", i_val=ri_data%min_bsize_MO)
      CALL section_vals_val_get(ri_section, "MEMORY_CUT", i_val=ri_data%n_mem)

      IF (ri_data%flavor == ri_pmat) THEN
         ri_data%n_mem = FLOOR(SQRT(REAL(ri_data%n_mem, dp) - 0.1)) + 1
      ENDIF

      ri_data%orb_basis_type = orb_basis_type
      ri_data%ri_basis_type = ri_basis_type
      ri_data%nelectron_total = nelectron_total

      IF (PRESENT(t_c_filename)) THEN
         ri_data%ri_metric%filename = t_c_filename
         ri_data%hfx_pot%filename = t_c_filename
      ENDIF

      ri_data%unit_nr_dbcsr = unit_nr_dbcsr
      ri_data%unit_nr = unit_nr

      ri_data%loc_subsection => section_vals_get_subs_vals(ri_section, "LOCALIZE")
      ri_data%print_loc_subsection => section_vals_get_subs_vals(ri_data%loc_subsection, "PRINT")

      CALL section_vals_val_get(ri_data%loc_subsection, "_SECTION_PARAMETERS_", l_val=ri_data%do_loc)

      CALL hfx_ri_init(ri_data, qs_kind_set, particle_set, atomic_kind_set, para_env, auto_basis_cntrl)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param ri_data ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param para_env ...
!> \param basis_cntrl ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_init(ri_data, qs_kind_set, particle_set, atomic_kind_set, para_env, basis_cntrl)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type)                             :: para_env
      INTEGER, INTENT(IN)                                :: basis_cntrl

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_mem, ikind, iproc, j_mem, &
                                                            MO_dim, mp_comm_3d, natom, nkind, nproc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bsizes_AO_store, bsizes_RI_store, dist1, &
                                                            dist2, dist_AO_1, dist_AO_2, dist_RI
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      LOGICAL                                            :: same_op
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: ri_aux_basis_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      NULLIFY (ri_aux_basis_set, qs_kind)

      natom = SIZE(particle_set)
      nkind = SIZE(qs_kind_set, 1)

      ! initialise RI basis
      DO ikind = 1, nkind
         NULLIFY (ri_aux_basis_set)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=ri_aux_basis_set, &
                          basis_type=ri_data%ri_basis_type)
         IF (.NOT. (ASSOCIATED(ri_aux_basis_set))) THEN
            ! Generate a default basis
            CALL create_ri_aux_basis_set(ri_aux_basis_set, qs_kind, basis_cntrl, ri_data%orb_basis_type)

            CALL add_basis_set_to_container(qs_kind%basis_sets, ri_aux_basis_set, ri_data%ri_basis_type)
         END IF
         CALL init_orb_basis_set(ri_aux_basis_set)
         CALL init_interaction_radii_orb_basis(ri_aux_basis_set, ri_data%eps_pgf_orb)
      END DO

      ASSOCIATE (ri_metric=>ri_data%ri_metric, hfx_pot=>ri_data%hfx_pot)
         IF (ri_metric%potential_type == do_potential_short) THEN
            CALL erfc_cutoff(ri_data%eps_schwarz, ri_metric%omega, ri_metric%cutoff_radius)
         ENDIF

         IF (hfx_pot%potential_type == do_potential_short) THEN
            ! need a more accurate threshold for determining 2-center integral operator range
            ! because stability of matrix inversion/sqrt is sensitive to this
            CALL erfc_cutoff(ri_data%filter_eps_2c, hfx_pot%omega, hfx_pot%cutoff_radius)
         ENDIF
         ! determine whether RI metric is same operator as used in HFX
         same_op = ri_metric%potential_type == hfx_pot%potential_type

         IF (same_op .AND. hfx_pot%potential_type == do_potential_truncated) THEN
            same_op = ABS(ri_metric%cutoff_radius - hfx_pot%cutoff_radius) < 1.0E-16_dp
         ENDIF

         IF (same_op .AND. hfx_pot%potential_type == do_potential_short) THEN
            same_op = ABS(ri_metric%omega - hfx_pot%omega) < 1.0E-16_dp
         ENDIF
      END ASSOCIATE

      ri_data%same_op = same_op

      pdims = 0
      CALL mp_cart_create(para_env%group, 3, pdims, pcoord, mp_comm_3d)

      ALLOCATE (ri_data%bsizes_RI(natom))
      ALLOCATE (ri_data%bsizes_AO(natom))
      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=ri_data%bsizes_RI, basis=basis_set_RI)
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=ri_data%bsizes_AO, basis=basis_set_AO)

      ALLOCATE (dist_RI(natom))
      ALLOCATE (dist_AO_1(natom))
      ALLOCATE (dist_AO_2(natom))
      CALL dbcsr_t_default_distvec(natom, pdims(1), ri_data%bsizes_RI, dist_RI)
      CALL dbcsr_t_default_distvec(natom, pdims(2), ri_data%bsizes_AO, dist_AO_1)
      CALL dbcsr_t_default_distvec(natom, pdims(3), ri_data%bsizes_AO, dist_AO_2)
      CALL distribution_3d_create(dist_3d, dist_RI, dist_ao_1, dist_ao_2, nkind, particle_set, &
                                  mp_comm_3d, own_comm=.TRUE.)

      ALLOCATE (ri_data%pgrid)
      CALL dbcsr_t_pgrid_create(para_env%group, pdims, ri_data%pgrid)

      ALLOCATE (ri_data%pgrid_2d)
      pdims_2d = 0
      CALL dbcsr_t_pgrid_create(para_env%group, pdims_2d, ri_data%pgrid_2d)

      ri_data%dist_3d = dist_3d

      CALL dbcsr_t_distribution_new(ri_data%dist, ri_data%pgrid, &
                                    dist_RI, dist_AO_1, dist_AO_2)

      DEALLOCATE (dist_AO_1, dist_AO_2, dist_RI)

      ! initialize tensors expressed in basis representation

      CALL pgf_block_sizes(atomic_kind_set, basis_set_AO, ri_data%min_bsize, ri_data%bsizes_AO_split)
      CALL pgf_block_sizes(atomic_kind_set, basis_set_RI, ri_data%min_bsize, ri_data%bsizes_RI_split)

      CALL pgf_block_sizes(atomic_kind_set, basis_set_AO, 1, bsizes_AO_store)
      CALL pgf_block_sizes(atomic_kind_set, basis_set_RI, 1, bsizes_RI_store)

      CALL split_block_sizes([SUM(ri_data%bsizes_AO)], ri_data%bsizes_AO_fit, default_block_size)
      CALL split_block_sizes([SUM(ri_data%bsizes_RI)], ri_data%bsizes_RI_fit, default_block_size)

      IF (ri_data%flavor == ri_pmat) THEN

         CALL create_tensor_batches(ri_data%bsizes_AO_split, ri_data%n_mem, ri_data%starts_array_mem, ri_data%ends_array_mem, &
                                    ri_data%starts_array_mem_block, ri_data%ends_array_mem_block)

         CALL create_tensor_batches(ri_data%bsizes_RI_split, ri_data%n_mem, &
                                    ri_data%starts_array_RI_mem, ri_data%ends_array_RI_mem, &
                                    ri_data%starts_array_RI_mem_block, ri_data%ends_array_RI_mem_block)

         ALLOCATE (ri_data%t_3c_int_ctr_1(ri_data%n_mem, ri_data%n_mem))
         CALL create_3c_tensor(ri_data%t_3c_int_ctr_1(1, 1), ri_data%dist1_ao_1, ri_data%dist1_ri, ri_data%dist1_ao_2, &
                               ri_data%pgrid, ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, &
                               [1, 2], [3], &
                               starts_array_block_1=ri_data%starts_array_mem_block, &
                               ends_array_block_1=ri_data%ends_array_mem_block, &
                               starts_array_block_2=ri_data%starts_array_RI_mem_block, &
                               ends_array_block_2=ri_data%ends_array_RI_mem_block, &
                               name="(AO RI | AO)")

         DO i_mem = 1, ri_data%n_mem
            DO j_mem = 1, ri_data%n_mem
               IF (.NOT. (i_mem == 1 .AND. j_mem == 1)) &
                  CALL dbcsr_t_create(ri_data%t_3c_int_ctr_1(1, 1), ri_data%t_3c_int_ctr_1(i_mem, j_mem))
            ENDDO
         ENDDO

         ALLOCATE (ri_data%blk_indices(ri_data%n_mem, ri_data%n_mem))
         ALLOCATE (ri_data%store_3c(ri_data%n_mem, ri_data%n_mem))
         DO i_mem = 1, ri_data%n_mem
         DO j_mem = 1, ri_data%n_mem
            CALL alloc_containers(ri_data%store_3c(i_mem, j_mem), 1)
         ENDDO
         ENDDO

         ALLOCATE (ri_data%t_3c_int_ctr_2(1, 1))
         CALL create_3c_tensor(ri_data%t_3c_int_ctr_2(1, 1), ri_data%dist1_ao_1, ri_data%dist1_ri, ri_data%dist1_ao_2, &
                               ri_data%pgrid, ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, &
                               [1, 2], [3], &
                               starts_array_block_2=ri_data%starts_array_RI_mem_block, &
                               ends_array_block_2=ri_data%ends_array_RI_mem_block, &
                               starts_array_block_3=ri_data%starts_array_mem_block, &
                               ends_array_block_3=ri_data%ends_array_mem_block, &
                               name="(AO RI | AO)")

         ALLOCATE (ri_data%t_3c_int_ctr_3(1, 1))
         CALL create_3c_tensor(ri_data%t_3c_int_ctr_3(1, 1), ri_data%dist3_RI, ri_data%dist3_AO_1, ri_data%dist3_AO_2, &
                               ri_data%pgrid, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                               [1], [2, 3], &
                               starts_array_block_1=ri_data%starts_array_RI_mem_block, &
                               ends_array_block_1=ri_data%ends_array_RI_mem_block, &
                               starts_array_block_2=ri_data%starts_array_mem_block, &
                               ends_array_block_2=ri_data%ends_array_mem_block, &
                               name="(RI | AO AO)")

         ALLOCATE (ri_data%t_2c_int(1, 1))
         CALL create_2c_tensor(ri_data%t_2c_int(1, 1), dist1, dist2, ri_data%pgrid_2d, &
                               ri_data%bsizes_RI_split, ri_data%bsizes_RI_split, &
                               name="(RI | RI)")
         DEALLOCATE (dist1, dist2)

      ELSEIF (ri_data%flavor == ri_mo) THEN
         ALLOCATE (ri_data%t_3c_int(1, 1))
         ALLOCATE (ri_data%t_2c_int(2, 1))

         CALL create_2c_tensor(ri_data%t_2c_int(1, 1), dist1, dist2, ri_data%pgrid_2d, &
                               ri_data%bsizes_RI_fit, ri_data%bsizes_RI_fit, &
                               name="(RI | RI)")

         DEALLOCATE (dist1, dist2)

         CALL mp_environ(nproc, iproc, para_env%group)

         ALLOCATE (ri_data%t_3c_int_ctr_1(1, 1))

         ALLOCATE (ri_data%pgrid_1)
         ALLOCATE (ri_data%pgrid_2)
         pdims = 0

         IF (ri_data%n_mem > ri_data%nelectron_total/2) ri_data%n_mem = MAX(ri_data%nelectron_total/2, 1)
         ! Size of dimension corresponding to MOs is nelectron/2 and divided by the memory factor
         ! we are using ceiling of that division to make sure that no MO dimension (after memory cut)
         ! is larger than this (it is however not a problem for load balancing if actual MO dimension
         ! is slightly smaller)
         MO_dim = MAX((ri_data%nelectron_total/2 - 1)/ri_data%n_mem + 1, 1)
         MO_dim = (MO_dim - 1)/ri_data%min_bsize_MO + 1

         pdims = 0
         CALL dbcsr_t_mp_dims_create(nproc, pdims, [SIZE(ri_data%bsizes_AO_split), SIZE(ri_data%bsizes_RI_split), MO_dim])

         CALL dbcsr_t_pgrid_create(para_env%group, pdims, ri_data%pgrid_1)

         pdims = pdims([3, 2, 1])
         CALL dbcsr_t_pgrid_create(para_env%group, pdims, ri_data%pgrid_2)

         CALL create_3c_tensor(ri_data%t_3c_int_ctr_1(1, 1), ri_data%dist1_ao_1, ri_data%dist1_ri, ri_data%dist1_ao_2, &
                               ri_data%pgrid_1, ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, &
                               [1, 2], [3], name="(AO RI | AO)")

         ALLOCATE (ri_data%t_3c_int_ctr_2(1, 1))
         CALL create_3c_tensor(ri_data%t_3c_int_ctr_2(1, 1), ri_data%dist2_ao_1, ri_data%dist2_ri, ri_data%dist2_ao_2, &
                               ri_data%pgrid_2, ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, &
                               [1], [2, 3], name="(AO | RI AO)")

      ENDIF

      ri_data%dbcsr_nflop = 0
      ri_data%dbcsr_time = 0.0_dp
      ri_data%num_pe = para_env%num_pe

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_write_stats(ri_data)
      TYPE(hfx_ri_type), INTENT(IN)                      :: ri_data

      REAL(dp)                                           :: my_flop_rate

      ASSOCIATE (unit_nr=>ri_data%unit_nr, dbcsr_nflop=>ri_data%dbcsr_nflop, &
                 dbcsr_time=>ri_data%dbcsr_time, num_pe=>ri_data%num_pe)
         my_flop_rate = REAL(dbcsr_nflop, dp)/(1.0E09_dp*ri_data%dbcsr_time)
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(/T2,A,T73,ES8.2)") &
            "RI-HFX PERFORMANCE| DBCSR total number of flops:", REAL(dbcsr_nflop*num_pe, dp)
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T2,A,T66,F15.2)") &
            "RI-HFX PERFORMANCE| DBCSR total execution time:", dbcsr_time
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T2,A,T66,F15.2)") &
            "RI-HFX PERFORMANCE| DBCSR flop rate (Gflops / MPI rank):", my_flop_rate
      END ASSOCIATE
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_release(ri_data)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_mem, ispin, j_mem, unused

      CALL timeset(routineN, handle)

      CALL hfx_ri_write_stats(ri_data)

      IF (ASSOCIATED(ri_data%pgrid)) THEN
         CALL dbcsr_t_pgrid_destroy(ri_data%pgrid)
         DEALLOCATE (ri_data%pgrid)
      ENDIF
      IF (ASSOCIATED(ri_data%pgrid_1)) THEN
         CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_1)
         DEALLOCATE (ri_data%pgrid_1)
      ENDIF
      IF (ASSOCIATED(ri_data%pgrid_2)) THEN
         CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_2)
         DEALLOCATE (ri_data%pgrid_2)
      ENDIF
      IF (ASSOCIATED(ri_data%pgrid_3)) THEN
         CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_3)
         DEALLOCATE (ri_data%pgrid_3)
      ENDIF
      IF (ASSOCIATED(ri_data%pgrid_2d)) THEN
         CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_2d)
         DEALLOCATE (ri_data%pgrid_2d)
      ENDIF

      CALL distribution_3d_destroy(ri_data%dist_3d)
      CALL dbcsr_t_distribution_destroy(ri_data%dist)

      DEALLOCATE (ri_data%bsizes_RI)
      DEALLOCATE (ri_data%bsizes_AO)
      DEALLOCATE (ri_data%bsizes_AO_split)
      DEALLOCATE (ri_data%bsizes_RI_split)
      DEALLOCATE (ri_data%bsizes_AO_fit)
      DEALLOCATE (ri_data%bsizes_RI_fit)

      IF (ri_data%flavor == ri_pmat) THEN
         DO i_mem = 1, ri_data%n_mem
         DO j_mem = 1, ri_data%n_mem
            CALL dealloc_containers(ri_data%store_3c(i_mem, j_mem), unused)
            CALL dbcsr_t_destroy(ri_data%t_3c_int_ctr_1(i_mem, j_mem))
         ENDDO
         ENDDO
         DEALLOCATE (ri_data%t_3c_int_ctr_1)
         CALL dbcsr_t_destroy(ri_data%t_3c_int_ctr_2(1, 1))
         DEALLOCATE (ri_data%t_3c_int_ctr_2)
         DEALLOCATE (ri_data%dist1_RI, ri_data%dist1_AO_1, ri_data%dist1_AO_2)
         DEALLOCATE (ri_data%dist3_RI, ri_data%dist3_AO_1, ri_data%dist3_AO_2)
         CALL dbcsr_t_destroy(ri_data%t_3c_int_ctr_3(1, 1))
         DEALLOCATE (ri_data%t_3c_int_ctr_3)
         CALL dbcsr_t_destroy(ri_data%t_2c_int(1, 1))
         DEALLOCATE (ri_data%t_2c_int)
         DEALLOCATE (ri_data%starts_array_mem_block, ri_data%ends_array_mem_block, &
                     ri_data%starts_array_mem, ri_data%ends_array_mem)
         DEALLOCATE (ri_data%starts_array_RI_mem_block, ri_data%ends_array_RI_mem_block, &
                     ri_data%starts_array_RI_mem, ri_data%ends_array_RI_mem)

         DEALLOCATE (ri_data%blk_indices)
         DEALLOCATE (ri_data%store_3c)
      ELSEIF (ri_data%flavor == ri_mo) THEN
         CALL dbcsr_t_destroy(ri_data%t_3c_int_ctr_1(1, 1))
         CALL dbcsr_t_destroy(ri_data%t_3c_int_ctr_2(1, 1))
         DEALLOCATE (ri_data%t_3c_int_ctr_1)
         DEALLOCATE (ri_data%t_3c_int_ctr_2)

         DO ispin = 1, SIZE(ri_data%t_3c_int_mo, 1)
            CALL dbcsr_t_destroy(ri_data%t_3c_int_mo(ispin, 1, 1))
            CALL dbcsr_t_destroy(ri_data%t_3c_ctr_RI(ispin, 1, 1))
            CALL dbcsr_t_destroy(ri_data%t_3c_ctr_KS(ispin, 1, 1))
            CALL dbcsr_t_destroy(ri_data%t_3c_ctr_KS_copy(ispin, 1, 1))
            CALL dbcsr_t_batched_contract_finalize(ri_data%t_2c_int(ispin, 1))
            CALL dbcsr_t_destroy(ri_data%t_2c_int(ispin, 1))
         ENDDO
         DEALLOCATE (ri_data%t_2c_int)
         DEALLOCATE (ri_data%t_3c_int_mo)
         DEALLOCATE (ri_data%t_3c_ctr_RI)
         DEALLOCATE (ri_data%t_3c_ctr_KS)
         DEALLOCATE (ri_data%t_3c_ctr_KS_copy)
      ENDIF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief - This routine allocates and initializes the basis_info and basis_parameter types
!> \param basis_parameter ...
!> \param basis_info ...
!> \param qs_kind_set ...
!> \param do_admm ...
!> \par History
!>      07.2011 refactored
! **************************************************************************************************
   SUBROUTINE hfx_create_basis_types(basis_parameter, basis_info, qs_kind_set, do_admm)
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      TYPE(hfx_basis_info_type)                          :: basis_info
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      LOGICAL, INTENT(IN)                                :: do_admm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create_basis_types', &
         routineP = moduleN//':'//routineN

      INTEGER :: co_counter, handle, i, ikind, ipgf, iset, j, k, la, max_am_kind, max_coeff, &
         max_nsgfl, max_pgf, max_pgf_kind, max_set, nkind, nl_count, nset, nseta, offset_a, &
         offset_a1, s_offset_nl_a, sgfa, so_counter
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, npgfa, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, nl_a
      REAL(dp), DIMENSION(:, :), POINTER                 :: sphi_a
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_a

      CALL timeset(routineN, handle)

      ! BASIS parameter
      nkind = SIZE(qs_kind_set, 1)
      !
      ALLOCATE (basis_parameter(nkind))
      max_set = 0
      DO ikind = 1, nkind
         IF (.NOT. do_admm) THEN
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_a, basis_type="ORB")
            CALL get_qs_kind_set(qs_kind_set, &
                                 maxsgf=basis_info%max_sgf, &
                                 maxnset=basis_info%max_set, &
                                 maxlgto=basis_info%max_am)

         ELSE
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_a, basis_type="AUX_FIT")
            CALL get_qs_kind_set(qs_kind_set, &
                                 maxsgf=basis_info%max_sgf, &
                                 maxnset=basis_info%max_set, &
                                 maxlgto=basis_info%max_am, &
                                 basis_type="AUX_FIT")
         END IF
         IF (basis_info%max_set < max_set) CPABORT("UNEXPECTED MAX_SET")
         max_set = MAX(max_set, basis_info%max_set)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_a, &
                                lmax=basis_parameter(ikind)%lmax, &
                                lmin=basis_parameter(ikind)%lmin, &
                                npgf=basis_parameter(ikind)%npgf, &
                                nset=basis_parameter(ikind)%nset, &
                                zet=basis_parameter(ikind)%zet, &
                                nsgf_set=basis_parameter(ikind)%nsgf, &
                                first_sgf=basis_parameter(ikind)%first_sgf, &
                                sphi=basis_parameter(ikind)%sphi, &
                                nsgf=basis_parameter(ikind)%nsgf_total, &
                                l=basis_parameter(ikind)%nl, &
                                nshell=basis_parameter(ikind)%nshell, &
                                set_radius=basis_parameter(ikind)%set_radius, &
                                pgf_radius=basis_parameter(ikind)%pgf_radius, &
                                kind_radius=basis_parameter(ikind)%kind_radius)
      END DO
      DO ikind = 1, nkind
         ALLOCATE (basis_parameter(ikind)%nsgfl(0:basis_info%max_am, max_set))
         basis_parameter(ikind)%nsgfl = 0
         nset = basis_parameter(ikind)%nset
         nshell => basis_parameter(ikind)%nshell
         DO iset = 1, nset
            DO i = 0, basis_info%max_am
               nl_count = 0
               DO j = 1, nshell(iset)
                  IF (basis_parameter(ikind)%nl(j, iset) == i) nl_count = nl_count + 1
               END DO
               basis_parameter(ikind)%nsgfl(i, iset) = nl_count
            END DO
         END DO
      END DO

      max_nsgfl = 0
      max_pgf = 0
      DO ikind = 1, nkind
         max_coeff = 0
         max_am_kind = 0
         max_pgf_kind = 0
         npgfa => basis_parameter(ikind)%npgf
         nseta = basis_parameter(ikind)%nset
         nl_a => basis_parameter(ikind)%nsgfl
         la_max => basis_parameter(ikind)%lmax
         la_min => basis_parameter(ikind)%lmin
         DO iset = 1, nseta
            max_pgf_kind = MAX(max_pgf_kind, npgfa(iset))
            max_pgf = MAX(max_pgf, npgfa(iset))
            DO la = la_min(iset), la_max(iset)
               max_nsgfl = MAX(max_nsgfl, nl_a(la, iset))
               max_coeff = MAX(max_coeff, nso(la)*nl_a(la, iset)*nco(la))
               max_am_kind = MAX(max_am_kind, la)
            END DO
         END DO
         ALLOCATE (basis_parameter(ikind)%sphi_ext(max_coeff, 0:max_am_kind, max_pgf_kind, nseta))
         basis_parameter(ikind)%sphi_ext = 0.0_dp
      END DO

      DO ikind = 1, nkind
         sphi_a => basis_parameter(ikind)%sphi
         nseta = basis_parameter(ikind)%nset
         la_max => basis_parameter(ikind)%lmax
         la_min => basis_parameter(ikind)%lmin
         npgfa => basis_parameter(ikind)%npgf
         first_sgfa => basis_parameter(ikind)%first_sgf
         nl_a => basis_parameter(ikind)%nsgfl
         DO iset = 1, nseta
            sgfa = first_sgfa(1, iset)
            DO ipgf = 1, npgfa(iset)
               offset_a1 = (ipgf - 1)*ncoset(la_max(iset))
               s_offset_nl_a = 0
               DO la = la_min(iset), la_max(iset)
                  offset_a = offset_a1 + ncoset(la - 1)
                  co_counter = 0
                  co_counter = co_counter + 1
                  so_counter = 0
                  DO k = sgfa + s_offset_nl_a, sgfa + s_offset_nl_a + nso(la)*nl_a(la, iset) - 1
                     DO i = offset_a + 1, offset_a + nco(la)
                        so_counter = so_counter + 1
                        basis_parameter(ikind)%sphi_ext(so_counter, la, ipgf, iset) = sphi_a(i, k)
                     END DO
                  END DO
                  s_offset_nl_a = s_offset_nl_a + nso(la)*(nl_a(la, iset))
               END DO
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_create_basis_types

! **************************************************************************************************
!> \brief ...
!> \param basis_parameter ...
! **************************************************************************************************
   SUBROUTINE hfx_release_basis_types(basis_parameter)
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_release_basis_types', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      !! BASIS parameter
      DO i = 1, SIZE(basis_parameter)
         DEALLOCATE (basis_parameter(i)%nsgfl)
         DEALLOCATE (basis_parameter(i)%sphi_ext)
      END DO
      DEALLOCATE (basis_parameter)
      CALL timestop(handle)

   END SUBROUTINE hfx_release_basis_types

! **************************************************************************************************
!> \brief - Parses the memory section
!> \param memory_parameter ...
!> \param hf_sub_section ...
!> \param storage_id ...
!> \param i_thread ...
!> \param n_threads ...
!> \param para_env ...
!> \param irep ...
!> \param skip_disk ...
!> \param skip_in_core_forces ...
! **************************************************************************************************
   SUBROUTINE parse_memory_section(memory_parameter, hf_sub_section, storage_id, &
                                   i_thread, n_threads, para_env, irep, skip_disk, skip_in_core_forces)
      TYPE(hfx_memory_type)                              :: memory_parameter
      TYPE(section_vals_type), POINTER                   :: hf_sub_section
      INTEGER, INTENT(OUT), OPTIONAL                     :: storage_id
      INTEGER, INTENT(IN), OPTIONAL                      :: i_thread, n_threads
      TYPE(cp_para_env_type), OPTIONAL                   :: para_env
      INTEGER, INTENT(IN), OPTIONAL                      :: irep
      LOGICAL, INTENT(IN)                                :: skip_disk, skip_in_core_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'parse_memory_section', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=512)                                 :: error_msg
      CHARACTER(LEN=default_path_length)                 :: char_val, filename, orig_wd
      INTEGER                                            :: int_val, stat
      LOGICAL                                            :: check, logic_val
      REAL(dp)                                           :: real_val

      check = (PRESENT(storage_id) .EQV. PRESENT(i_thread)) .AND. &
              (PRESENT(storage_id) .EQV. PRESENT(n_threads)) .AND. &
              (PRESENT(storage_id) .EQV. PRESENT(para_env)) .AND. &
              (PRESENT(storage_id) .EQV. PRESENT(irep))
      CPASSERT(check)

      ! Memory Storage
      CALL section_vals_val_get(hf_sub_section, "MAX_MEMORY", i_val=int_val)
      memory_parameter%max_memory = int_val
      memory_parameter%max_compression_counter = int_val*1024_int_8*128_int_8
      CALL section_vals_val_get(hf_sub_section, "EPS_STORAGE", r_val=real_val)
      memory_parameter%eps_storage_scaling = real_val
      IF (int_val == 0) THEN
         memory_parameter%do_all_on_the_fly = .TRUE.
      ELSE
         memory_parameter%do_all_on_the_fly = .FALSE.
      END IF
      memory_parameter%cache_size = CACHE_SIZE
      memory_parameter%bits_max_val = BITS_MAX_VAL
      memory_parameter%actual_memory_usage = 1
      IF (.NOT. skip_in_core_forces) THEN
         CALL section_vals_val_get(hf_sub_section, "TREAT_FORCES_IN_CORE", l_val=logic_val)
         memory_parameter%treat_forces_in_core = logic_val
      END IF

      ! ** IF MAX_MEM == 0 overwrite this flag to false
      IF (memory_parameter%do_all_on_the_fly) memory_parameter%treat_forces_in_core = .FALSE.

      ! Disk Storage
      IF (.NOT. skip_disk) THEN
         memory_parameter%actual_memory_usage_disk = 1
         CALL section_vals_val_get(hf_sub_section, "MAX_DISK_SPACE", i_val=int_val)
         memory_parameter%max_compression_counter_disk = int_val*1024_int_8*128_int_8
         IF (int_val == 0) THEN
            memory_parameter%do_disk_storage = .FALSE.
         ELSE
            memory_parameter%do_disk_storage = .TRUE.
         END IF
         CALL section_vals_val_get(hf_sub_section, "STORAGE_LOCATION", c_val=char_val)
         CALL compress(char_val, .TRUE.)
         !! Add ending / if necessary

         IF (SCAN(char_val, "/", .TRUE.) /= LEN_TRIM(char_val)) THEN
            WRITE (filename, '(A,A)') TRIM(char_val), "/"
            CALL compress(filename)
         ELSE
            filename = TRIM(char_val)
         END IF
         CALL compress(filename, .TRUE.)

         !! quickly check if we can write on storage_location
         CALL m_getcwd(orig_wd)
         CALL m_chdir(TRIM(filename), stat)
         IF (stat /= 0) THEN
            WRITE (error_msg, '(A,A,A)') "Request for disk storage failed due to unknown error while writing to ", &
               TRIM(filename), ". Please check STORAGE_LOCATION"
            CPABORT(error_msg)
         END IF
         CALL m_chdir(orig_wd, stat)

         memory_parameter%storage_location = filename
         CALL compress(memory_parameter%storage_location, .TRUE.)
      ELSE
         memory_parameter%do_disk_storage = .FALSE.
      END IF
      IF (PRESENT(storage_id)) THEN
         storage_id = (irep - 1)*para_env%num_pe*n_threads + para_env%mepos*n_threads + i_thread - 1
      END IF
   END SUBROUTINE parse_memory_section

! **************************************************************************************************
!> \brief - This routine deallocates all data structures
!> \param x_data contains all relevant data structures for hfx runs
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE hfx_release(x_data)
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, i_thread, irep, n_rep_hf, n_threads
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(hfx_type), POINTER                            :: actual_x_data

!! There might be 2 hf sections

      n_rep_hf = x_data(1, 1)%n_rep_hf
      n_threads = SIZE(x_data, 2)

      IF (x_data(1, 1)%potential_parameter%potential_type == do_potential_truncated .OR. &
          x_data(1, 1)%potential_parameter%potential_type == do_potential_mix_cl_trunc) THEN
         init_t_c_g0_lmax = -1
         CALL free_C0()
      END IF
      DO i_thread = 1, n_threads
         DO irep = 1, n_rep_hf
            actual_x_data => x_data(irep, i_thread)
            DEALLOCATE (actual_x_data%neighbor_cells)
            DEALLOCATE (actual_x_data%distribution_energy)
            DEALLOCATE (actual_x_data%distribution_forces)

            IF (actual_x_data%load_balance_parameter%blocks_initialized) THEN
               DEALLOCATE (actual_x_data%blocks)
               IF (i_thread == 1) THEN
                  DEALLOCATE (actual_x_data%pmax_block)
               END IF
            END IF

            IF (i_thread == 1) THEN
               DEALLOCATE (actual_x_data%atomic_pair_list)
               DEALLOCATE (actual_x_data%atomic_pair_list_forces)
            END IF

            IF (actual_x_data%screening_parameter%do_initial_p_screening .OR. &
                actual_x_data%screening_parameter%do_p_screening_forces) THEN
               IF (i_thread == 1) THEN
                  DEALLOCATE (actual_x_data%pmax_atom)
                  DO i = 1, SIZE(actual_x_data%initial_p)
                     DEALLOCATE (actual_x_data%initial_p(i)%p_kind)
                  END DO
                  DEALLOCATE (actual_x_data%initial_p)
                  IF (actual_x_data%memory_parameter%treat_forces_in_core) THEN
                     DEALLOCATE (actual_x_data%pmax_atom_forces)
                     DO i = 1, SIZE(actual_x_data%initial_p_forces)
                        DEALLOCATE (actual_x_data%initial_p_forces(i)%p_kind)
                     END DO
                     DEALLOCATE (actual_x_data%initial_p_forces)
                  END IF
               END IF
               DEALLOCATE (actual_x_data%map_atom_to_kind_atom)
            END IF
            IF (i_thread == 1) THEN
               DEALLOCATE (actual_x_data%is_assoc_atomic_block)
               DEALLOCATE (actual_x_data%atomic_block_offset)
               DEALLOCATE (actual_x_data%set_offset)
               DEALLOCATE (actual_x_data%block_offset)
            END IF

            !! BASIS parameter
            CALL hfx_release_basis_types(actual_x_data%basis_parameter)

            !MK Release libint and libderiv data structure
            CALL cp_libint_cleanup_eri(actual_x_data%lib)
            CALL cp_libint_cleanup_eri1(actual_x_data%lib_deriv)
            CALL cp_libint_static_cleanup()

            !! Deallocate containers
            CALL dealloc_containers(actual_x_data%store_ints, actual_x_data%memory_parameter%actual_memory_usage)
            CALL dealloc_containers(actual_x_data%store_forces, actual_x_data%memory_parameter%actual_memory_usage)

            !! Deallocate containers
            CALL hfx_init_container(actual_x_data%store_ints%maxval_container_disk, &
                                    actual_x_data%memory_parameter%actual_memory_usage_disk, &
                                    .FALSE.)
            IF (actual_x_data%memory_parameter%do_disk_storage) THEN
               CALL close_file(unit_number=actual_x_data%store_ints%maxval_container_disk%unit, file_status="DELETE")
            END IF
            DEALLOCATE (actual_x_data%store_ints%maxval_container_disk%first)
            DEALLOCATE (actual_x_data%store_ints%maxval_container_disk)

            DO i = 1, 64
               CALL hfx_init_container(actual_x_data%store_ints%integral_containers_disk(i), &
                                       actual_x_data%memory_parameter%actual_memory_usage_disk, &
                                       .FALSE.)
               IF (actual_x_data%memory_parameter%do_disk_storage) THEN
                  CALL close_file(unit_number=actual_x_data%store_ints%integral_containers_disk(i)%unit, file_status="DELETE")
               END IF
               DEALLOCATE (actual_x_data%store_ints%integral_containers_disk(i)%first)
            END DO
            DEALLOCATE (actual_x_data%store_ints%integral_containers_disk)

            ! ** screening functions
            IF (actual_x_data%screen_funct_is_initialized) THEN
               DEALLOCATE (actual_x_data%screen_funct_coeffs_set)
               DEALLOCATE (actual_x_data%screen_funct_coeffs_kind)
               DEALLOCATE (actual_x_data%pair_dist_radii_pgf)
               DEALLOCATE (actual_x_data%screen_funct_coeffs_pgf)
               actual_x_data%screen_funct_is_initialized = .FALSE.
            END IF

            ! ** maps
            IF (ASSOCIATED(actual_x_data%map_atoms_to_cpus)) THEN
               DO i = 1, SIZE(actual_x_data%map_atoms_to_cpus)
                  DEALLOCATE (actual_x_data%map_atoms_to_cpus(i)%iatom_list)
                  DEALLOCATE (actual_x_data%map_atoms_to_cpus(i)%jatom_list)
               END DO
               DEALLOCATE (actual_x_data%map_atoms_to_cpus)
            END IF

            IF (actual_x_data%do_hfx_ri) THEN
               CALL hfx_ri_release(actual_x_data%ri_data)
               IF (ASSOCIATED(actual_x_data%ri_data%ri_section)) THEN
                  logger => cp_get_default_logger()
                  CALL cp_print_key_finished_output(actual_x_data%ri_data%unit_nr_dbcsr, logger, actual_x_data%ri_data%ri_section, &
                                                    "RI_INFO")
               ENDIF
               IF (ASSOCIATED(actual_x_data%ri_data%hfx_section)) THEN
                  logger => cp_get_default_logger()
                  CALL cp_print_key_finished_output(actual_x_data%ri_data%unit_nr, logger, actual_x_data%ri_data%hfx_section, &
                                                    "HF_INFO")
               ENDIF
               DEALLOCATE (actual_x_data%ri_data)
            ENDIF
         END DO

      END DO

      DEALLOCATE (x_data)
   END SUBROUTINE hfx_release

! **************************************************************************************************
!> \brief - This routine computes the neighbor cells that are taken into account
!>        in periodic runs
!> \param x_data contains all relevant data structures for hfx runs
!> \param pbc_shells number of shells taken into account
!> \param cell cell
!> \param i_thread current thread ID
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE hfx_create_neighbor_cells(x_data, pbc_shells, cell, i_thread)
      TYPE(hfx_type), POINTER                            :: x_data
      INTEGER, INTENT(INOUT)                             :: pbc_shells
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: i_thread

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create_neighbor_cells', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=512)                                 :: error_msg
      CHARACTER(LEN=64)                                  :: char_nshells
      INTEGER :: i, idx, ikind, ipgf, iset, ishell, j, jkind, jpgf, jset, jshell, k, kshell, l, &
         m(3), max_shell, nseta, nsetb, total_number_of_cells, ub, ub_max
      INTEGER, DIMENSION(:), POINTER                     :: la_max, lb_max, npgfa, npgfb
      LOGICAL                                            :: image_cell_found, nothing_more_to_add
      REAL(dp) :: cross_product(3), dist_min, distance(14), l_min, normal(3, 6), P(3, 14), &
         plane_vector(3, 2), point_in_plane(3), r(3), R1, R_max, R_max_stress, s(3), x, y, z, Zeta1
      REAL(dp), DIMENSION(:, :), POINTER                 :: zeta, zetb
      TYPE(hfx_cell_type), ALLOCATABLE, DIMENSION(:)     :: tmp_neighbor_cells

      total_number_of_cells = 0

      ! ** Check some settings
      IF (i_thread == 1) THEN
         IF (x_data%potential_parameter%potential_type /= do_potential_truncated .AND. &
             x_data%potential_parameter%potential_type /= do_potential_short .AND. &
             x_data%potential_parameter%potential_type /= do_potential_mix_cl_trunc .AND. &
             x_data%potential_parameter%potential_type /= do_potential_id) THEN
            CALL cp_warn(__LOCATION__, &
                         "Periodic Hartree Fock calculation requested without use "// &
                         "of a truncated or shortrange potential. This may lead to unphysical total energies. "// &
                         "Use a truncated  potential to avoid possible problems.")
         ELSE
            l_min = MIN(plane_distance(1, 0, 0, cell), &
                        plane_distance(0, 1, 0, cell), &
                        plane_distance(0, 0, 1, cell))
            l_min = 0.5_dp*l_min
            IF (x_data%potential_parameter%cutoff_radius >= l_min) THEN
               CALL cp_warn(__LOCATION__, &
                            "Periodic Hartree Fock calculation requested with use "// &
                            "of a truncated or shortrange potential. The cutoff radius is larger than half "// &
                            "the minimal cell dimension. This may lead to unphysical "// &
                            "total energies. Reduce the cutoff radius in order to avoid "// &
                            "possible problems.")
            END IF
         END IF
      END IF

      SELECT CASE (x_data%potential_parameter%potential_type)
      CASE (do_potential_truncated, do_potential_mix_cl_trunc, do_potential_short)
         R_max = 0.0_dp
         DO ikind = 1, SIZE(x_data%basis_parameter)
            la_max => x_data%basis_parameter(ikind)%lmax
            zeta => x_data%basis_parameter(ikind)%zet
            nseta = x_data%basis_parameter(ikind)%nset
            npgfa => x_data%basis_parameter(ikind)%npgf
            DO jkind = 1, SIZE(x_data%basis_parameter)
               lb_max => x_data%basis_parameter(jkind)%lmax
               zetb => x_data%basis_parameter(jkind)%zet
               nsetb = x_data%basis_parameter(jkind)%nset
               npgfb => x_data%basis_parameter(jkind)%npgf
               DO iset = 1, nseta
                  DO jset = 1, nsetb
                     DO ipgf = 1, npgfa(iset)
                        DO jpgf = 1, npgfb(jset)
                           Zeta1 = zeta(ipgf, iset) + zetb(jpgf, jset)
                           R1 = 1.0_dp/SQRT(Zeta1)*mul_fact(la_max(iset) + lb_max(jset))* &
                                SQRT(-LOG(x_data%screening_parameter%eps_schwarz))
                           R_max = MAX(R1, R_max)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

         R_max = 2.0_dp*R_max + x_data%potential_parameter%cutoff_radius
         nothing_more_to_add = .FALSE.
         max_shell = 0
         total_number_of_cells = 0
         ub = 1
         DEALLOCATE (x_data%neighbor_cells)
         ALLOCATE (x_data%neighbor_cells(1))
         x_data%neighbor_cells(1)%cell = 0.0_dp
         x_data%neighbor_cells(1)%cell_r = 0.0_dp

         ! ** What follows is kind of a ray tracing algorithm
         ! ** Given a image cell (ishell, jshell, kshell) we try to figure out the
         ! ** shortest distance of this image cell to the basic unit cell (0,0,0), i.e. the point
         ! ** (0.0, 0.0, 0.0)
         ! ** This is achieved by checking the 8 Corners of the cell, and, in addition, the shortest distance
         ! ** to all 6 faces. The faces are only taken into account if the penetration point of the normal
         ! ** to the plane defined by a face lies within this face.
         ! ** This is very fast, because no trigonometric functions are being used
         ! ** The points are defined as follows
         ! **
         ! **
         ! **               _________________________
         ! **              /P4____________________P8/|
         ! **             / / ___________________/ / |
         ! **            / / /| |               / /  |       z
         ! **           / / / | |              / / . |      /|\  _ y
         ! **          / / /| | |             / / /| |       |   /|
         ! **         / / / | | |            / / / | |       |  /
         ! **        / / /  | | |           / / /| | |       | /
         ! **       / /_/___| | |__________/ / / | | |       |/
         ! **      /P2______| | |_________P6/ /  | | |       ----------> x
         ! **      | _______| | |_________| | |  | | |
         ! **      | | |    | | |________________| | |
         ! **      | | |    |P3___________________P7 |
         ! **      | | |   / / _________________  / /
         ! **      | | |  / / /           | | |/ / /
         ! **      | | | / / /            | | | / /
         ! **      | | |/ / /             | | |/ /
         ! **      | | | / /              | | ' /
         ! **      | | |/_/_______________| |  /
         ! **      | |____________________| | /
         ! **      |P1_____________________P5/
         ! **
         ! **

         DO WHILE (.NOT. nothing_more_to_add)
            ! Calculate distances to the eight points P1 to P8
            image_cell_found = .FALSE.
            ALLOCATE (tmp_neighbor_cells(1:ub))
            DO i = 1, ub - 1
               tmp_neighbor_cells(i) = x_data%neighbor_cells(i)
            END DO
            ub_max = (2*max_shell + 1)**3
            DEALLOCATE (x_data%neighbor_cells)
            ALLOCATE (x_data%neighbor_cells(1:ub_max))
            DO i = 1, ub - 1
               x_data%neighbor_cells(i) = tmp_neighbor_cells(i)
            END DO
            DO i = ub, ub_max
               x_data%neighbor_cells(i)%cell = 0.0_dp
               x_data%neighbor_cells(i)%cell_r = 0.0_dp
            END DO

            DEALLOCATE (tmp_neighbor_cells)

            DO ishell = -max_shell, max_shell
            DO jshell = -max_shell, max_shell
            DO kshell = -max_shell, max_shell
               IF (MAX(ABS(ishell), ABS(jshell), ABS(kshell)) /= max_shell) CYCLE
               idx = 0
               DO j = 0, 1
                  x = -1.0_dp/2.0_dp + j*1.0_dp
                  DO k = 0, 1
                     y = -1.0_dp/2.0_dp + k*1.0_dp
                     DO l = 0, 1
                        z = -1.0_dp/2.0_dp + l*1.0_dp
                        idx = idx + 1
                        P(1, idx) = x + ishell
                        P(2, idx) = y + jshell
                        P(3, idx) = z + kshell
                        CALL scaled_to_real(r, P(:, idx), cell)
                        distance(idx) = SQRT(SUM(r**2))
                        P(1:3, idx) = r
                     END DO
                  END DO
               END DO
               ! Now check distance to Faces and only take them into account if the base point lies within quadrilateral

               ! Face A (1342) 1 is the reference
               idx = idx + 1
               plane_vector(:, 1) = P(:, 3) - P(:, 1)
               plane_vector(:, 2) = P(:, 2) - P(:, 1)
               cross_product(1) = plane_vector(2, 1)*plane_vector(3, 2) - plane_vector(3, 1)*plane_vector(2, 2)
               cross_product(2) = plane_vector(3, 1)*plane_vector(1, 2) - plane_vector(1, 1)*plane_vector(3, 2)
               cross_product(3) = plane_vector(1, 1)*plane_vector(2, 2) - plane_vector(2, 1)*plane_vector(1, 2)
               normal(:, 1) = cross_product/SQRT(SUM(cross_product**2))
               point_in_plane = -normal(:, 1)*(normal(1, 1)*P(1, 1) + normal(2, 1)*P(2, 1) + normal(3, 1)*P(3, 1))

               IF (point_is_in_quadrilateral(P(:, 1), P(:, 3), P(:, 4), P(:, 2), point_in_plane)) THEN
                  distance(idx) = ABS(normal(1, 1)*P(1, 1) + normal(2, 1)*P(2, 1) + normal(3, 1)*P(3, 1))
               ELSE
                  distance(idx) = HUGE(distance(idx))
               END IF

               ! Face B (1562) 1 is the reference
               idx = idx + 1
               plane_vector(:, 1) = P(:, 2) - P(:, 1)
               plane_vector(:, 2) = P(:, 5) - P(:, 1)
               cross_product(1) = plane_vector(2, 1)*plane_vector(3, 2) - plane_vector(3, 1)*plane_vector(2, 2)
               cross_product(2) = plane_vector(3, 1)*plane_vector(1, 2) - plane_vector(1, 1)*plane_vector(3, 2)
               cross_product(3) = plane_vector(1, 1)*plane_vector(2, 2) - plane_vector(2, 1)*plane_vector(1, 2)
               normal(:, 1) = cross_product/SQRT(SUM(cross_product**2))
               point_in_plane = -normal(:, 1)*(normal(1, 1)*P(1, 1) + normal(2, 1)*P(2, 1) + normal(3, 1)*P(3, 1))

               IF (point_is_in_quadrilateral(P(:, 1), P(:, 5), P(:, 6), P(:, 2), point_in_plane)) THEN
                  distance(idx) = ABS(normal(1, 1)*P(1, 1) + normal(2, 1)*P(2, 1) + normal(3, 1)*P(3, 1))
               ELSE
                  distance(idx) = HUGE(distance(idx))
               END IF

               ! Face C (5786) 5 is the reference
               idx = idx + 1
               plane_vector(:, 1) = P(:, 7) - P(:, 5)
               plane_vector(:, 2) = P(:, 6) - P(:, 5)
               cross_product(1) = plane_vector(2, 1)*plane_vector(3, 2) - plane_vector(3, 1)*plane_vector(2, 2)
               cross_product(2) = plane_vector(3, 1)*plane_vector(1, 2) - plane_vector(1, 1)*plane_vector(3, 2)
               cross_product(3) = plane_vector(1, 1)*plane_vector(2, 2) - plane_vector(2, 1)*plane_vector(1, 2)
               normal(:, 1) = cross_product/SQRT(SUM(cross_product**2))
               point_in_plane = -normal(:, 1)*(normal(1, 1)*P(1, 5) + normal(2, 1)*P(2, 5) + normal(3, 1)*P(3, 5))

               IF (point_is_in_quadrilateral(P(:, 5), P(:, 7), P(:, 8), P(:, 6), point_in_plane)) THEN
                  distance(idx) = ABS(normal(1, 1)*P(1, 5) + normal(2, 1)*P(2, 5) + normal(3, 1)*P(3, 5))
               ELSE
                  distance(idx) = HUGE(distance(idx))
               END IF

               ! Face D (3784) 3 is the reference
               idx = idx + 1
               plane_vector(:, 1) = P(:, 7) - P(:, 3)
               plane_vector(:, 2) = P(:, 4) - P(:, 3)
               cross_product(1) = plane_vector(2, 1)*plane_vector(3, 2) - plane_vector(3, 1)*plane_vector(2, 2)
               cross_product(2) = plane_vector(3, 1)*plane_vector(1, 2) - plane_vector(1, 1)*plane_vector(3, 2)
               cross_product(3) = plane_vector(1, 1)*plane_vector(2, 2) - plane_vector(2, 1)*plane_vector(1, 2)
               normal(:, 1) = cross_product/SQRT(SUM(cross_product**2))
               point_in_plane = -normal(:, 1)*(normal(1, 1)*P(1, 3) + normal(2, 1)*P(2, 3) + normal(3, 1)*P(3, 3))

               IF (point_is_in_quadrilateral(P(:, 3), P(:, 7), P(:, 8), P(:, 4), point_in_plane)) THEN
                  distance(idx) = ABS(normal(1, 1)*P(1, 3) + normal(2, 1)*P(2, 3) + normal(3, 1)*P(3, 3))
               ELSE
                  distance(idx) = HUGE(distance(idx))
               END IF

               ! Face E (2684) 2 is the reference
               idx = idx + 1
               plane_vector(:, 1) = P(:, 6) - P(:, 2)
               plane_vector(:, 2) = P(:, 4) - P(:, 2)
               cross_product(1) = plane_vector(2, 1)*plane_vector(3, 2) - plane_vector(3, 1)*plane_vector(2, 2)
               cross_product(2) = plane_vector(3, 1)*plane_vector(1, 2) - plane_vector(1, 1)*plane_vector(3, 2)
               cross_product(3) = plane_vector(1, 1)*plane_vector(2, 2) - plane_vector(2, 1)*plane_vector(1, 2)
               normal(:, 1) = cross_product/SQRT(SUM(cross_product**2))
               point_in_plane = -normal(:, 1)*(normal(1, 1)*P(1, 2) + normal(2, 1)*P(2, 2) + normal(3, 1)*P(3, 2))

               IF (point_is_in_quadrilateral(P(:, 2), P(:, 6), P(:, 8), P(:, 4), point_in_plane)) THEN
                  distance(idx) = ABS(normal(1, 1)*P(1, 2) + normal(2, 1)*P(2, 2) + normal(3, 1)*P(3, 2))
               ELSE
                  distance(idx) = HUGE(distance(idx))
               END IF

               ! Face F (1573) 1 is the reference
               idx = idx + 1
               plane_vector(:, 1) = P(:, 5) - P(:, 1)
               plane_vector(:, 2) = P(:, 3) - P(:, 1)
               cross_product(1) = plane_vector(2, 1)*plane_vector(3, 2) - plane_vector(3, 1)*plane_vector(2, 2)
               cross_product(2) = plane_vector(3, 1)*plane_vector(1, 2) - plane_vector(1, 1)*plane_vector(3, 2)
               cross_product(3) = plane_vector(1, 1)*plane_vector(2, 2) - plane_vector(2, 1)*plane_vector(1, 2)
               normal(:, 1) = cross_product/SQRT(SUM(cross_product**2))
               point_in_plane = -normal(:, 1)*(normal(1, 1)*P(1, 1) + normal(2, 1)*P(2, 1) + normal(3, 1)*P(3, 1))

               IF (point_is_in_quadrilateral(P(:, 1), P(:, 5), P(:, 7), P(:, 3), point_in_plane)) THEN
                  distance(idx) = ABS(normal(1, 1)*P(1, 1) + normal(2, 1)*P(2, 1) + normal(3, 1)*P(3, 1))
               ELSE
                  distance(idx) = HUGE(distance(idx))
               END IF

               dist_min = MINVAL(distance)
               IF (max_shell == 0) THEN
                  image_cell_found = .TRUE.
               END IF
               IF (dist_min < R_max) THEN
                  total_number_of_cells = total_number_of_cells + 1
                  x_data%neighbor_cells(ub)%cell = REAL((/ishell, jshell, kshell/), dp)
                  ub = ub + 1
                  image_cell_found = .TRUE.
               END IF

            END DO
            END DO
            END DO
            IF (image_cell_found) THEN
               max_shell = max_shell + 1
            ELSE
               nothing_more_to_add = .TRUE.
            END IF
         END DO
         ! now remove what is not needed
         ALLOCATE (tmp_neighbor_cells(total_number_of_cells))
         DO i = 1, ub - 1
            tmp_neighbor_cells(i) = x_data%neighbor_cells(i)
         END DO
         DEALLOCATE (x_data%neighbor_cells)
         ! If we only need the supercell, total_number_of_cells is still 0, repair
         IF (total_number_of_cells == 0) THEN
            total_number_of_cells = 1
            ALLOCATE (x_data%neighbor_cells(total_number_of_cells))
            DO i = 1, total_number_of_cells
               x_data%neighbor_cells(i)%cell = 0.0_dp
               x_data%neighbor_cells(i)%cell_r = 0.0_dp
            END DO
         ELSE
            ALLOCATE (x_data%neighbor_cells(total_number_of_cells))
            DO i = 1, total_number_of_cells
               x_data%neighbor_cells(i) = tmp_neighbor_cells(i)
            END DO
         END IF
         DEALLOCATE (tmp_neighbor_cells)

         IF (x_data%periodic_parameter%number_of_shells == do_hfx_auto_shells) THEN
            ! Do nothing
         ELSE
            total_number_of_cells = 0
            DO i = 0, x_data%periodic_parameter%number_of_shells
               total_number_of_cells = total_number_of_cells + count_cells_perd(i, x_data%periodic_parameter%perd)
            END DO
            IF (total_number_of_cells < SIZE(x_data%neighbor_cells)) THEN
               IF (i_thread == 1) THEN

                  WRITE (char_nshells, '(I3)') SIZE(x_data%neighbor_cells)
                  WRITE (error_msg, '(A,A,A)') "Periodic Hartree Fock calculation requested with use "// &
                     "of a truncated potential. The number of shells to be considered "// &
                     "might be too small. CP2K conservatively estimates to need "//TRIM(char_nshells)//" periodic images "// &
                     "Please carefully check if you get converged results."
                  CPWARN(error_msg)
               END IF
            END IF
            total_number_of_cells = 0
            DO i = 0, x_data%periodic_parameter%number_of_shells
               total_number_of_cells = total_number_of_cells + count_cells_perd(i, x_data%periodic_parameter%perd)
            END DO
            DEALLOCATE (x_data%neighbor_cells)

            ALLOCATE (x_data%neighbor_cells(total_number_of_cells))
            m = 0
            i = 1
            DO WHILE (SUM(m**2) <= x_data%periodic_parameter%number_of_shells)
               x_data%neighbor_cells(i)%cell = REAL(m, dp)
               CALL next_image_cell_perd(m, x_data%periodic_parameter%perd)
               i = i + 1
            ENDDO
         END IF
      CASE DEFAULT
         total_number_of_cells = 0
         IF (pbc_shells == -1) pbc_shells = 0
         DO i = 0, pbc_shells
            total_number_of_cells = total_number_of_cells + count_cells_perd(i, x_data%periodic_parameter%perd)
         END DO
         DEALLOCATE (x_data%neighbor_cells)

         ALLOCATE (x_data%neighbor_cells(total_number_of_cells))

         m = 0
         i = 1
         DO WHILE (SUM(m**2) <= pbc_shells)
            x_data%neighbor_cells(i)%cell = REAL(m, dp)
            CALL next_image_cell_perd(m, x_data%periodic_parameter%perd)
            i = i + 1
         ENDDO
      END SELECT

      ! ** Transform into real coord
      DO i = 1, SIZE(x_data%neighbor_cells)
         r = 0.0_dp
         x_data%neighbor_cells(i)%cell_r(:) = 0.0_dp
         s = x_data%neighbor_cells(i)%cell(:)
         CALL scaled_to_real(x_data%neighbor_cells(i)%cell_r, s, cell)
      END DO
      x_data%periodic_parameter%number_of_shells = pbc_shells

      R_max_stress = 0.0_dp
      DO i = 1, SIZE(x_data%neighbor_cells)
         R_max_stress = MAX(R_max_stress, MAXVAL(ABS(x_data%neighbor_cells(i)%cell_r(:))))
      END DO
      R_max_stress = R_max_stress + ABS(MAXVAL(cell%hmat(:, :)))
      x_data%periodic_parameter%R_max_stress = R_max_stress

   END SUBROUTINE hfx_create_neighbor_cells

   ! performs a fuzzy check of being in a quadrilateral
! **************************************************************************************************
!> \brief ...
!> \param A ...
!> \param B ...
!> \param C ...
!> \param D ...
!> \param P ...
!> \return ...
! **************************************************************************************************
   FUNCTION point_is_in_quadrilateral(A, B, C, D, P)
      REAL(dp)                                           :: A(3), B(3), C(3), D(3), P(3)
      LOGICAL                                            :: point_is_in_quadrilateral

      REAL(dp), PARAMETER :: fuzzy = 1000.0_dp*EPSILON(1.0_dp)

      REAL(dp)                                           :: dot00, dot01, dot02, dot11, dot12, &
                                                            invDenom, u, v, v0(3), v1(3), v2(3)

      point_is_in_quadrilateral = .FALSE.

      ! ** Check for both triangles ABC and ACD
      ! **
      ! **     D -------------- C
      ! **    /                /
      ! **   /                /
      ! **  A----------------B
      ! **
      ! **
      ! **

      ! ** ABC

      v0 = D - A
      v1 = C - A
      v2 = P - A

      ! ** Compute dot products
      dot00 = DOT_PRODUCT(v0, v0)
      dot01 = DOT_PRODUCT(v0, v1)
      dot02 = DOT_PRODUCT(v0, v2)
      dot11 = DOT_PRODUCT(v1, v1)
      dot12 = DOT_PRODUCT(v1, v2)

      ! ** Compute barycentric coordinates
      invDenom = 1/(dot00*dot11 - dot01*dot01)
      u = (dot11*dot02 - dot01*dot12)*invDenom
      v = (dot00*dot12 - dot01*dot02)*invDenom
      ! ** Check if point is in triangle
      IF ((u >= 0 - fuzzy) .AND. (v >= 0 - fuzzy) .AND. (u + v <= 1 + fuzzy)) THEN
         point_is_in_quadrilateral = .TRUE.
         RETURN
      END IF
      v0 = C - A
      v1 = B - A
      v2 = P - A

      ! ** Compute dot products
      dot00 = DOT_PRODUCT(v0, v0)
      dot01 = DOT_PRODUCT(v0, v1)
      dot02 = DOT_PRODUCT(v0, v2)
      dot11 = DOT_PRODUCT(v1, v1)
      dot12 = DOT_PRODUCT(v1, v2)

      ! ** Compute barycentric coordinates
      invDenom = 1/(dot00*dot11 - dot01*dot01)
      u = (dot11*dot02 - dot01*dot12)*invDenom
      v = (dot00*dot12 - dot01*dot02)*invDenom

      ! ** Check if point is in triangle
      IF ((u >= 0 - fuzzy) .AND. (v >= 0 - fuzzy) .AND. (u + v <= 1 + fuzzy)) THEN
         point_is_in_quadrilateral = .TRUE.
         RETURN
      END IF

   END FUNCTION point_is_in_quadrilateral

! **************************************************************************************************
!> \brief - This routine deletes all list entries in a container in order to
!>        deallocate the memory.
!> \param container container that contains the compressed elements
!> \param memory_usage ...
!> \param do_disk_storage ...
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE hfx_init_container(container, memory_usage, do_disk_storage)
      TYPE(hfx_container_type)                           :: container
      INTEGER                                            :: memory_usage
      LOGICAL                                            :: do_disk_storage

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_init_container', &
         routineP = moduleN//':'//routineN

      TYPE(hfx_container_node), POINTER                  :: current, next

!! DEALLOCATE memory

      current => container%first
      DO WHILE (ASSOCIATED(current))
         next => current%next
         DEALLOCATE (current)
         current => next
      END DO

      !! Allocate first list entry, init members
      ALLOCATE (container%first)
      container%first%prev => NULL()
      container%first%next => NULL()
      container%current => container%first
      container%current%data = 0
      container%element_counter = 1
      memory_usage = 1

      IF (do_disk_storage) THEN
         !! close the file, if this is no the first time
         IF (container%unit /= -1) THEN
            CALL close_file(unit_number=container%unit)
         END IF
         CALL open_file(file_name=TRIM(container%filename), file_status="UNKNOWN", file_form="UNFORMATTED", file_action="WRITE", &
                        unit_number=container%unit)
      END IF

   END SUBROUTINE hfx_init_container

! **************************************************************************************************
!> \brief - This routine stores the data obtained from the load balance routine
!>        for the energy
!> \param ptr_to_distr contains data to store
!> \param x_data contains all relevant data structures for hfx runs
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE hfx_set_distr_energy(ptr_to_distr, x_data)
      TYPE(hfx_distribution), DIMENSION(:), POINTER      :: ptr_to_distr
      TYPE(hfx_type), POINTER                            :: x_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_set_distr_energy', &
         routineP = moduleN//':'//routineN

      DEALLOCATE (x_data%distribution_energy)

      ALLOCATE (x_data%distribution_energy(SIZE(ptr_to_distr)))
      x_data%distribution_energy = ptr_to_distr

   END SUBROUTINE hfx_set_distr_energy

! **************************************************************************************************
!> \brief - This routine stores the data obtained from the load balance routine
!>        for the forces
!> \param ptr_to_distr contains data to store
!> \param x_data contains all relevant data structures for hfx runs
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE hfx_set_distr_forces(ptr_to_distr, x_data)
      TYPE(hfx_distribution), DIMENSION(:), POINTER      :: ptr_to_distr
      TYPE(hfx_type), POINTER                            :: x_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_set_distr_forces', &
         routineP = moduleN//':'//routineN

      DEALLOCATE (x_data%distribution_forces)

      ALLOCATE (x_data%distribution_forces(SIZE(ptr_to_distr)))
      x_data%distribution_forces = ptr_to_distr

   END SUBROUTINE hfx_set_distr_forces

! **************************************************************************************************
!> \brief - resets the maximum memory usage for a HFX calculation subtracting
!>          all relevant buffers from the input MAX_MEM value and add 10% of
!>          safety margin
!> \param memory_parameter Memory information
!> \param subtr_size_mb size of buffers in MiB
!> \par History
!>      02.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE hfx_reset_memory_usage_counter(memory_parameter, subtr_size_mb)

      TYPE(hfx_memory_type)                              :: memory_parameter
      INTEGER(int_8), INTENT(IN)                         :: subtr_size_mb

      INTEGER(int_8)                                     :: max_memory

      max_memory = memory_parameter%max_memory
      max_memory = max_memory - subtr_size_mb
      IF (max_memory <= 0) THEN
         memory_parameter%do_all_on_the_fly = .TRUE.
         memory_parameter%max_compression_counter = 0
      ELSE
         memory_parameter%do_all_on_the_fly = .FALSE.
         memory_parameter%max_compression_counter = max_memory*1024_int_8*128_int_8
      END IF
   END SUBROUTINE hfx_reset_memory_usage_counter

! **************************************************************************************************
!> \brief - This routine prints some information on HFX
!> \param x_data contains all relevant data structures for hfx runs
!> \param hfx_section HFX input section
!> \par History
!>      03.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! **************************************************************************************************
   SUBROUTINE hfx_print_std_info(x_data, hfx_section)
      TYPE(hfx_type), POINTER                            :: x_data
      TYPE(section_vals_type), POINTER                   :: hfx_section

      INTEGER                                            :: iw
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      iw = cp_print_key_unit_nr(logger, hfx_section, "HF_INFO", &
                                extension=".scfLog")

      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="((T3,A,T73,ES8.1))") &
            "HFX_INFO| EPS_SCHWARZ:     ", x_data%screening_parameter%eps_schwarz
         WRITE (UNIT=iw, FMT="((T3,A,T73,ES8.1))") &
            "HFX_INFO| EPS_SCHWARZ_FORCES     ", x_data%screening_parameter%eps_schwarz_forces
         WRITE (UNIT=iw, FMT="((T3,A,T73,ES8.1))") &
            "HFX_INFO| EPS_STORAGE_SCALING:     ", x_data%memory_parameter%eps_storage_scaling
         WRITE (UNIT=iw, FMT="((T3,A,T61,I20))") &
            "HFX_INFO| NBINS:     ", x_data%load_balance_parameter%nbins
         WRITE (UNIT=iw, FMT="((T3,A,T61,I20))") &
            "HFX_INFO| BLOCK_SIZE:     ", x_data%load_balance_parameter%block_size
         IF (x_data%periodic_parameter%do_periodic) THEN
            IF (x_data%periodic_parameter%mode == -1) THEN
               WRITE (UNIT=iw, FMT="((T3,A,T77,A))") &
                  "HFX_INFO| NUMBER_OF_SHELLS:     ", "AUTO"
            ELSE
               WRITE (UNIT=iw, FMT="((T3,A,T61,I20))") &
                  "HFX_INFO| NUMBER_OF_SHELLS:     ", x_data%periodic_parameter%mode
            END IF
            WRITE (UNIT=iw, FMT="((T3,A,T61,I20))") &
               "HFX_INFO| Number of periodic shells considered:     ", x_data%periodic_parameter%number_of_shells
            WRITE (UNIT=iw, FMT="((T3,A,T61,I20),/)") &
               "HFX_INFO| Number of periodic cells considered:     ", SIZE(x_data%neighbor_cells)
         ELSE
            WRITE (UNIT=iw, FMT="((T3,A,T77,A))") &
               "HFX_INFO| Number of periodic shells considered:     ", "NONE"
            WRITE (UNIT=iw, FMT="((T3,A,T77,A),/)") &
               "HFX_INFO| Number of periodic cells considered:     ", "NONE"
         END IF
      ENDIF
   END SUBROUTINE hfx_print_std_info

! **************************************************************************************************
!> \brief ...
!> \param ri_data ...
!> \param hfx_section ...
! **************************************************************************************************
   SUBROUTINE hfx_print_ri_info(ri_data, hfx_section)
      TYPE(hfx_ri_type), POINTER                         :: ri_data
      TYPE(section_vals_type), POINTER                   :: hfx_section

      INTEGER                                            :: iw
      REAL(dp)                                           :: rc_ang
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: ri_section

      NULLIFY (logger, ri_section)
      logger => cp_get_default_logger()

      ri_section => ri_data%ri_section

      iw = cp_print_key_unit_nr(logger, hfx_section, "HF_INFO", &
                                extension=".scfLog")

      IF (iw > 0) THEN

         ASSOCIATE (ri_metric=>ri_data%ri_metric, hfx_pot=>ri_data%hfx_pot)
            SELECT CASE (ri_metric%potential_type)
            CASE (do_potential_coulomb)
               WRITE (UNIT=iw, FMT="(/T3,A,T74,A)") &
                  "HFX_RI_INFO| RI metric: ", "COULOMB"
            CASE (do_potential_short)
               WRITE (UNIT=iw, FMT="(T3,A,T71,A)") &
                  "HFX_RI_INFO| RI metric: ", "SHORTRANGE"
               WRITE (iw, '(T3,A,T61,F20.10)') &
                  "HFX_RI_INFO| Omega:     ", ri_metric%omega
               rc_ang = cp_unit_from_cp2k(ri_metric%cutoff_radius, "angstrom")
               WRITE (iw, '(T3,A,T61,F20.10)') &
                  "HFX_RI_INFO| Cutoff Radius [angstrom]:     ", rc_ang
            CASE (do_potential_long)
               WRITE (UNIT=iw, FMT="(T3,A,T72,A)") &
                  "HFX_RI_INFO| RI metric: ", "LONGRANGE"
               WRITE (iw, '(T3,A,T61,F20.10)') &
                  "HFX_RI_INFO| Omega:     ", ri_metric%omega
            CASE (do_potential_id)
               WRITE (UNIT=iw, FMT="(T3,A,T73,A)") &
                  "HFX_RI_INFO| RI metric: ", "IDENTITY"
            CASE (do_potential_truncated)
               WRITE (UNIT=iw, FMT="(T3,A,T72,A)") &
                  "HFX_RI_INFO| RI metric: ", "TRUNCATED"
               rc_ang = cp_unit_from_cp2k(ri_metric%cutoff_radius, "angstrom")
               WRITE (iw, '(T3,A,T61,F20.10)') &
                  "HFX_RI_INFO| Cutoff Radius [angstrom]:     ", rc_ang
            END SELECT

         END ASSOCIATE
         SELECT CASE (ri_data%flavor)
         CASE (ri_mo)
            WRITE (UNIT=iw, FMT="(T3, A, T79, A)") &
               "HFX_RI_INFO| RI flavor: ", "MO"
         CASE (ri_pmat)
            WRITE (UNIT=iw, FMT="(T3, A, T78, A)") &
               "HFX_RI_INFO| RI flavor: ", "RHO"
         END SELECT
         SELECT CASE (ri_data%do_loc)
         CASE (.TRUE.)
            WRITE (UNIT=iw, FMT="(T3, A, T80, A)") &
               "HFX_RI_INFO| Localized MOs: ", "T"
         CASE (.FALSE.)
            WRITE (UNIT=iw, FMT="(T3, A, T80, A)") &
               "HFX_RI_INFO| Localized MOs: ", "F"
         END SELECT
         SELECT CASE (ri_data%t2c_method)
         CASE (hfx_ri_do_2c_iter)
            WRITE (UNIT=iw, FMT="(T3, A, T69, A)") &
               "HFX_RI_INFO| Matrix SQRT/INV", "DBCSR / iter"
         CASE (hfx_ri_do_2c_diag)
            WRITE (UNIT=iw, FMT="(T3, A, T65, A)") &
               "HFX_RI_INFO| Matrix SQRT/INV", "Dense / diag"
         END SELECT
         WRITE (UNIT=iw, FMT="(T3, A, T73, ES8.1)") &
            "HFX_RI_INFO| EPS_FILTER", ri_data%filter_eps
         WRITE (UNIT=iw, FMT="(T3, A, T73, ES8.1)") &
            "HFX_RI_INFO| EPS_FILTER 2-center", ri_data%filter_eps_2c
         WRITE (UNIT=iw, FMT="(T3, A, T73, ES8.1)") &
            "HFX_RI_INFO| EPS_FILTER storage", ri_data%filter_eps_storage
         WRITE (UNIT=iw, FMT="(T3, A, T73, ES8.1)") &
            "HFX_RI_INFO| EPS_FILTER MO", ri_data%filter_eps_mo
         WRITE (UNIT=iw, FMT="(T3, A, T73, ES8.1)") &
            "HFX_RI_INFO| EPS_PGF_ORB", ri_data%eps_pgf_orb
         WRITE (UNIT=iw, FMT="((T3, A, T73, ES8.1))") &
            "HFX_RI_INFO| EPS_SCHWARZ:     ", ri_data%eps_schwarz
         WRITE (UNIT=iw, FMT="(T3, A, T78, I3)") &
            "HFX_RI_INFO| Minimum block size", ri_data%min_bsize
         WRITE (UNIT=iw, FMT="(T3, A, T78, I3)") &
            "HFX_RI_INFO| MO block size", ri_data%min_bsize_MO
         SELECT CASE (ri_data%flavor)
         CASE (ri_mo)
            WRITE (UNIT=iw, FMT="(T3, A, T79, I2)") &
               "HFX_RI_INFO| Memory reduction factor", ri_data%n_mem
         CASE (ri_pmat)
            WRITE (UNIT=iw, FMT="(T3, A, T79, I2)") &
               "HFX_RI_INFO| Memory reduction factor", ri_data%n_mem**2
         END SELECT
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param x_data ...
!> \param hfx_section ...
!> \param i_rep ...
! **************************************************************************************************
   SUBROUTINE hfx_print_info(x_data, hfx_section, i_rep)
      TYPE(hfx_type), POINTER                            :: x_data
      TYPE(section_vals_type), POINTER                   :: hfx_section
      INTEGER, INTENT(IN)                                :: i_rep

      INTEGER                                            :: iw
      REAL(dp)                                           :: rc_ang
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      iw = cp_print_key_unit_nr(logger, hfx_section, "HF_INFO", &
                                extension=".scfLog")

      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(/,(T3,A,T61,I20))") &
            "HFX_INFO| Replica ID:     ", i_rep

         WRITE (iw, '(T3,A,T61,F20.10)') &
            "HFX_INFO| FRACTION:     ", x_data%general_parameter%fraction
         SELECT CASE (x_data%potential_parameter%potential_type)
         CASE (do_potential_coulomb)
            WRITE (UNIT=iw, FMT="((T3,A,T74,A))") &
               "HFX_INFO| Interaction Potential:     ", "COULOMB"
         CASE (do_potential_short)
            WRITE (UNIT=iw, FMT="((T3,A,T71,A))") &
               "HFX_INFO| Interaction Potential:    ", "SHORTRANGE"
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
            rc_ang = cp_unit_from_cp2k(x_data%potential_parameter%cutoff_radius, "angstrom")
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| Cutoff Radius [angstrom]:     ", rc_ang
         CASE (do_potential_long)
            WRITE (UNIT=iw, FMT="((T3,A,T72,A))") &
               "HFX_INFO| Interaction Potential:     ", "LONGRANGE"
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
         CASE (do_potential_mix_cl)
            WRITE (UNIT=iw, FMT="((T3,A,T75,A))") &
               "HFX_INFO| Interaction Potential:     ", "MIX_CL"
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| SCALE_COULOMB:     ", x_data%potential_parameter%scale_coulomb
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| SCALE_LONGRANGE:     ", x_data%potential_parameter%scale_longrange
         CASE (do_potential_gaussian)
            WRITE (UNIT=iw, FMT="((T3,A,T73,A))") &
               "HFX_INFO| Interaction Potential:     ", "GAUSSIAN"
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
         CASE (do_potential_mix_lg)
            WRITE (UNIT=iw, FMT="((T3,A,T75,A))") &
               "HFX_INFO| Interaction Potential:    ", "MIX_LG"
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| SCALE_LONGRANGE:     ", x_data%potential_parameter%scale_longrange
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| SCALE_GAUSSIAN:    ", x_data%potential_parameter%scale_gaussian
         CASE (do_potential_id)
            WRITE (UNIT=iw, FMT="((T3,A,T73,A))") &
               "HFX_INFO| Interaction Potential:    ", "IDENTITY"
         CASE (do_potential_truncated)
            WRITE (UNIT=iw, FMT="((T3,A,T72,A))") &
               "HFX_INFO| Interaction Potential:    ", "TRUNCATED"
            rc_ang = cp_unit_from_cp2k(x_data%potential_parameter%cutoff_radius, "angstrom")
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| Cutoff Radius [angstrom]:     ", rc_ang
         CASE (do_potential_mix_cl_trunc)
            WRITE (UNIT=iw, FMT="((T3,A,T65,A))") &
               "HFX_INFO| Interaction Potential:    ", "TRUNCATED MIX_CL"
            rc_ang = cp_unit_from_cp2k(x_data%potential_parameter%cutoff_radius, "angstrom")
            WRITE (iw, '(T3,A,T61,F20.10)') &
               "HFX_INFO| Cutoff Radius [angstrom]:     ", rc_ang
         END SELECT

      END IF
      IF (x_data%do_hfx_ri) THEN
         CALL hfx_print_ri_info(x_data%ri_data, hfx_section)
      ELSE
         CALL hfx_print_std_info(x_data, hfx_section)
      ENDIF

      CALL cp_print_key_finished_output(iw, logger, hfx_section, &
                                        "HF_INFO")
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param DATA ...
!> \param memory_usage ...
! **************************************************************************************************
   SUBROUTINE dealloc_containers(DATA, memory_usage)
      TYPE(hfx_compression_type)                         :: data
      INTEGER                                            :: memory_usage

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dealloc_containers', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bin, i

      DO bin = 1, SIZE(data%maxval_container)
         CALL hfx_init_container(data%maxval_container(bin), memory_usage, &
                                 .FALSE.)
         DEALLOCATE (data%maxval_container(bin)%first)
      END DO
      DEALLOCATE (data%maxval_container)
      DEALLOCATE (data%maxval_cache)

      DO bin = 1, SIZE(data%integral_containers, 2)
         DO i = 1, 64
            CALL hfx_init_container(data%integral_containers(i, bin), memory_usage, &
                                    .FALSE.)
            DEALLOCATE (data%integral_containers(i, bin)%first)
         END DO
      END DO
      DEALLOCATE (data%integral_containers)

      DEALLOCATE (data%integral_caches)

   END SUBROUTINE dealloc_containers

! **************************************************************************************************
!> \brief ...
!> \param DATA ...
!> \param bin_size ...
! **************************************************************************************************
   SUBROUTINE alloc_containers(DATA, bin_size)
      TYPE(hfx_compression_type)                         :: data
      INTEGER, INTENT(IN)                                :: bin_size

      CHARACTER(LEN=*), PARAMETER :: routineN = 'alloc_containers', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bin, i

      ALLOCATE (data%maxval_cache(bin_size))
      DO bin = 1, bin_size
         data%maxval_cache(bin)%element_counter = 1
      END DO
      ALLOCATE (data%maxval_container(bin_size))
      DO bin = 1, bin_size
         ALLOCATE (data%maxval_container(bin)%first)
         data%maxval_container(bin)%first%prev => NULL()
         data%maxval_container(bin)%first%next => NULL()
         data%maxval_container(bin)%current => data%maxval_container(bin)%first
         data%maxval_container(bin)%current%data = 0
         data%maxval_container(bin)%element_counter = 1
      END DO

      ALLOCATE (data%integral_containers(64, bin_size))
      ALLOCATE (data%integral_caches(64, bin_size))

      DO bin = 1, bin_size
         DO i = 1, 64
            data%integral_caches(i, bin)%element_counter = 1
            data%integral_caches(i, bin)%data = 0
            ALLOCATE (data%integral_containers(i, bin)%first)
            data%integral_containers(i, bin)%first%prev => NULL()
            data%integral_containers(i, bin)%first%next => NULL()
            data%integral_containers(i, bin)%current => data%integral_containers(i, bin)%first
            data%integral_containers(i, bin)%current%data = 0
            data%integral_containers(i, bin)%element_counter = 1
         END DO
      END DO

   END SUBROUTINE alloc_containers

END MODULE hfx_types

