!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types and set/get functions for HFX
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE hfx_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
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
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE erf_fn,                          ONLY: erfc
  USE hfx_helpers,                     ONLY: count_cells_perd,&
                                             next_image_cell_perd
  USE hfx_libint_wrapper,              ONLY: initialize_libderiv,&
                                             initialize_libint,&
                                             terminate_libderiv,&
                                             terminate_libint
  USE hfx_libint_wrapper_types,        ONLY: lib_deriv,&
                                             lib_int,&
                                             prim_data_f_size
  USE input_constants,                 ONLY: &
       do_hfx_auto_shells, do_hfx_potential_coulomb, &
       do_hfx_potential_gaussian, do_hfx_potential_id, do_hfx_potential_long, &
       do_hfx_potential_mix_cl, do_hfx_potential_mix_cl_trunc, &
       do_hfx_potential_mix_lg, do_hfx_potential_short, &
       do_hfx_potential_truncated, hfx_do_eval_energy, hfx_do_eval_forces
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp,&
                                             int_8
  USE machine,                         ONLY: m_chdir,&
                                             m_getcwd
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: nco,&
                                             ncoset,&
                                             nso
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE string_utilities,                ONLY: compress
  USE t_c_g0,                          ONLY: free_C0
  USE timings,                         ONLY: timeset,&
                                             timestop

  !$ USE OMP_LIB

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hfx_type, hfx_create, hfx_release,&
            hfx_set_distr_energy,&
            hfx_set_distr_forces,&
            hfx_cell_loop_info_type,&
            hfx_cell_type, hfx_distribution,&
            hfx_potential_type, hfx_screening_type,&
            hfx_memory_type, hfx_load_balance_type, hfx_general_type,&
            hfx_container_type, hfx_cache_type,&
            hfx_basis_type, parse_memory_section,&
            hfx_init_container,&
            hfx_basis_info_type, hfx_screen_coeff_type,&
            hfx_reset_memory_usage_counter, pair_list_type, pair_list_element_type,&
            pair_set_list_type, hfx_p_kind, hfx_2D_map, hfx_pgf_list,&
            hfx_pgf_product_list, hfx_block_range_type,&
            alloc_containers, dealloc_containers, hfx_task_list_type, init_t_c_g0_lmax,&
            hfx_create_neighbor_cells, hfx_create_basis_types, hfx_release_basis_types,&
            set_cell_info, get_cell_index

#define CACHE_SIZE 1024
#define BITS_MAX_VAL 6

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_types'
  INTEGER, PARAMETER, PUBLIC                 :: max_atom_block=32
  INTEGER, PARAMETER, PUBLIC                 :: max_images=27
  REAL(dp), PARAMETER, PUBLIC                :: log_zero = -1000.0_dp
  REAL(dp), PARAMETER, PUBLIC                :: powell_min_log = -20.0_dp
  REAL(KIND=dp), DIMENSION(0:10), &
     PARAMETER, PUBLIC                       :: mul_fact = (/1.0_dp,&
                                                             1.1781_dp,&
                                                             1.3333_dp,&
                                                             1.4726_dp,&
                                                             1.6000_dp,&
                                                             1.7181_dp,&
                                                             1.8286_dp,&
                                                             1.9328_dp,&
                                                             2.0317_dp,&
                                                             2.1261_dp,&
                                                             2.2165_dp/)

  INTEGER, SAVE                                         :: init_t_c_g0_lmax = -1

!***

! *****************************************************************************
  TYPE hfx_potential_type
    INTEGER                                  :: potential_type     !! 1/r/ erfc(wr)/r ...
    REAL(dp)                                 :: omega              !! w
    REAL(dp)                                 :: scale_coulomb      !! scaling factor for mixed potential
    REAL(dp)                                 :: scale_longrange    !! scaling factor for mixed potential
    REAL(dp)                                 :: scale_gaussian     !! scaling factor for mixed potential
    REAL(dp)                                 :: cutoff_radius      !! cutoff radius if cutoff potential in use
    CHARACTER(default_path_length)           :: filename
  END TYPE

! *****************************************************************************
  TYPE hfx_screening_type
    REAL(dp)                                 :: eps_schwarz         !! threshold
    REAL(dp)                                 :: eps_schwarz_forces  !! threshold
    LOGICAL                                  :: do_p_screening_forces  !! screen on P^2 ?
    LOGICAL                                  :: do_initial_p_screening !! screen on inital guess?
  END TYPE
! *****************************************************************************
  TYPE hfx_cell_loop_info_type
    LOGICAL                                  :: dokp
    INTEGER                                  :: image_ks
    INTEGER, DIMENSION(3)                    :: cell_ks
    INTEGER                                  :: image_pm
    INTEGER, DIMENSION(3)                    :: cell_pm
    INTEGER, DIMENSION(:, :), POINTER        :: index_to_cell => NULL()
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index => NULL()
  END TYPE

! *****************************************************************************
  TYPE hfx_memory_type
    INTEGER                                  :: max_memory              !! user def max memory MiB
    INTEGER(int_8)                           :: max_compression_counter !! corresponding number of reals
    INTEGER(int_8)                           :: final_comp_counter_energy
    LOGICAL                                  :: do_all_on_the_fly       !! max mem == 0 ?
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

! *****************************************************************************
  TYPE hfx_periodic_type
    INTEGER                                  :: number_of_shells        !! number of periodic image cells
    LOGICAL                                  :: do_periodic             !! periodic ?
    INTEGER                                  :: perd(3)                 !! x,xy,xyz,...
    INTEGER                                  :: mode
    REAL(dp)                                 :: R_max_stress
    INTEGER                                  :: number_of_shells_from_input
  END TYPE

! *****************************************************************************
  TYPE hfx_load_balance_type
    INTEGER                                  :: nbins
    INTEGER                                  :: block_size
    INTEGER                                  :: nblocks
    LOGICAL                                  :: rtp_redistribute
    LOGICAL                                  :: blocks_initialized
    LOGICAL                                  :: do_randomize
  END TYPE

! *****************************************************************************
  TYPE hfx_general_type
    REAL(dp)                                 :: fraction                 !! for hybrids
    LOGICAL                                  :: treat_lsd_in_core
  END TYPE

! *****************************************************************************
  TYPE hfx_cell_type
    REAL(dp)                                 :: cell(3)
    REAL(dp)                                 :: cell_r(3)
  END TYPE

! *****************************************************************************
  TYPE hfx_distribution
    INTEGER(int_8)                           :: istart
    INTEGER(int_8)                           :: number_of_atom_quartets
    INTEGER(int_8)                           :: cost
    REAL(KIND=dp)                            :: time_first_scf
    REAL(KIND=dp)                            :: time_other_scf
    REAL(KIND=dp)                            :: time_forces
    INTEGER(int_8)                           :: ram_counter
  END TYPE

! *****************************************************************************
  TYPE pair_list_element_type
      INTEGER, DIMENSION(2) :: pair
      INTEGER, DIMENSION(2) :: set_bounds
      INTEGER, DIMENSION(2) :: kind_pair
      REAL(KIND=dp)         :: r1(3),r2(3)
      REAL(KIND=dp)         :: dist2
  END TYPE

  ! *****************************************************************************
  TYPE pair_set_list_type
      INTEGER, DIMENSION(2) :: pair
  END TYPE

! *****************************************************************************
  TYPE pair_list_type
    TYPE(pair_list_element_type), DIMENSION(max_atom_block**2) :: elements
    INTEGER :: n_element
  END TYPE pair_list_type

! *****************************************************************************
  TYPE hfx_compression_type
    INTEGER(int_8), DIMENSION(:), POINTER    :: compressed_data
    INTEGER(int_8)                           :: element_pointer
    INTEGER                                  :: bit_pointer
  END TYPE

! *****************************************************************************
  TYPE hfx_cache_type
    INTEGER(int_8), DIMENSION(CACHE_SIZE)    :: DATA
    INTEGER                                  :: element_counter
  END TYPE

! *****************************************************************************
  TYPE hfx_container_node
    TYPE(hfx_container_node), POINTER        :: next, prev
    INTEGER(int_8), DIMENSION(CACHE_SIZE)    :: DATA
  END TYPE

! *****************************************************************************
  TYPE hfx_container_type
    TYPE(hfx_container_node), POINTER        :: first, current
    INTEGER                                  :: element_counter
    INTEGER(int_8)                           :: file_counter
    CHARACTER(LEN=5)                         :: desc
    INTEGER                                  :: unit
    CHARACTER(default_path_length)           :: filename
  END TYPE

! *****************************************************************************
  TYPE hfx_basis_type
    INTEGER, DIMENSION(:), POINTER           :: lmax
    INTEGER, DIMENSION(:), POINTER           :: lmin
    INTEGER, DIMENSION(:), POINTER           :: npgf
    INTEGER                                  :: nset
    REAL(dp), DIMENSION(:,:), POINTER        :: zet
    INTEGER, DIMENSION(:), POINTER           :: nsgf
    INTEGER, DIMENSION(:,:), POINTER         :: first_sgf
    REAL(dp), DIMENSION(:,:), POINTER        :: sphi
    INTEGER                                  :: nsgf_total
    INTEGER, DIMENSION(:,:), POINTER         :: nl
    INTEGER, DIMENSION(:,:), POINTER         :: nsgfl
    INTEGER, DIMENSION(:), POINTER           :: nshell
    REAL(dp), DIMENSION(:,:,:,:), POINTER &
                                             :: sphi_ext
    REAL(dp), DIMENSION(:), POINTER          :: set_radius
    REAL(dp), DIMENSION(:,:), POINTER        :: pgf_radius
    REAL(dp)                                 :: kind_radius
  END TYPE

! *****************************************************************************
  TYPE hfx_basis_info_type
    INTEGER                                  :: max_set
    INTEGER                                  :: max_sgf
    INTEGER                                  :: max_am
  END TYPE

! *****************************************************************************
  TYPE hfx_screen_coeff_type
    REAL(dp)                                 :: x(2)
  END TYPE

! *****************************************************************************
  TYPE hfx_p_kind
    REAL(dp), DIMENSION(:,:,:,:), POINTER    :: p_kind
  END TYPE

! *****************************************************************************
  TYPE hfx_2D_map
    INTEGER, DIMENSION(:), POINTER           :: iatom_list
    INTEGER, DIMENSION(:), POINTER           :: jatom_list
  END TYPE

! *****************************************************************************
  TYPE hfx_pgf_image
    REAL(dp)                                 :: ra(3), rb(3)
    REAL(dp)                                 :: rab2
    REAL(dp)                                 :: S1234
    REAL(dp)                                 :: P(3)
    REAL(dp)                                 :: R
    REAL(dp)                                 :: pgf_max
    REAL(dp), DIMENSION(3)                   :: bcell
  END TYPE

! *****************************************************************************
  TYPE hfx_pgf_list
    TYPE(hfx_pgf_image), DIMENSION(:), POINTER &
                                             :: image_list
    INTEGER                                  :: nimages
    REAL(dp)                                 :: zetapzetb
    REAL(dp)                                 :: ZetaInv
    REAL(dp)                                 :: zeta, zetb
    INTEGER                                  :: ipgf, jpgf
  END TYPE

! *****************************************************************************
  TYPE hfx_pgf_product_list
    REAL(dp)                                 :: ra(3), rb(3), rc(3), rd(3)
    REAL(dp)                                 :: ZetapEtaInv
    REAL(dp)                                 :: Rho, RhoInv
    REAL(dp)                                 :: P(3), Q(3), W(3)
    REAL(dp)                                 :: AB(3), CD(3)
    REAL(dp)                                 :: Fm(prim_data_f_size)
  END TYPE

! *****************************************************************************
  TYPE hfx_block_range_type
    INTEGER        :: istart,iend
    INTEGER(int_8) :: cost
  END TYPE

! *****************************************************************************
  TYPE hfx_task_list_type
    INTEGER                                  :: thread_id
    INTEGER                                  :: bin_id
    INTEGER(int_8)                           :: cost
  END TYPE

! *****************************************************************************
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
! *****************************************************************************
  TYPE hfx_type
     TYPE(hfx_potential_type)                 :: potential_parameter
     TYPE(hfx_screening_type)                 :: screening_parameter
     TYPE(hfx_memory_type)                    :: memory_parameter
     TYPE(hfx_periodic_type)                  :: periodic_parameter
     TYPE(hfx_load_balance_type)              :: load_balance_parameter
     TYPE(hfx_general_type)                   :: general_parameter
     TYPE(hfx_container_type), DIMENSION(:),&
                               POINTER        :: maxval_container
     TYPE(hfx_cache_type), DIMENSION(:),&
                           POINTER            :: maxval_cache
     TYPE(hfx_container_type), DIMENSION(:,:), &
                               POINTER        :: integral_containers
     TYPE(hfx_cache_type), DIMENSION(:,:),&
                           POINTER            :: integral_caches
     TYPE(hfx_container_type), DIMENSION(:),&
                               POINTER        :: maxval_container_forces
     TYPE(hfx_cache_type), DIMENSION(:),&
                           POINTER            :: maxval_cache_forces
     TYPE(hfx_container_type), DIMENSION(:,:), &
                               POINTER        :: integral_containers_forces
     TYPE(hfx_cache_type), DIMENSION(:,:),&
                           POINTER            :: integral_caches_forces

     TYPE(hfx_container_type), POINTER        :: maxval_container_disk
     TYPE(hfx_cache_type)                     :: maxval_cache_disk
     TYPE(hfx_container_type), POINTER, &
                                DIMENSION(:)  :: integral_containers_disk
     TYPE(hfx_cache_type)                     :: integral_caches_disk(64)

     TYPE(hfx_cell_type), DIMENSION(:), &
                POINTER                       :: neighbor_cells
     TYPE(hfx_distribution), DIMENSION(:),&
                              POINTER         :: distribution_energy
     TYPE(hfx_distribution), DIMENSION(:),&
                              POINTER         :: distribution_forces
     INTEGER, DIMENSION(:,:), POINTER         :: is_assoc_atomic_block
     INTEGER                                  :: number_of_p_entries
     TYPE(hfx_basis_type), DIMENSION(:), &
                            POINTER           :: basis_parameter
     INTEGER                                  :: n_rep_hf
     LOGICAL                                  :: b_first_load_balance_energy, &
                                                 b_first_load_balance_forces
     REAL(dp), DIMENSION(:,:), POINTER        :: full_ks_alpha
     REAL(dp), DIMENSION(:,:), POINTER        :: full_ks_beta
     TYPE(lib_int)                            :: lib
     TYPE(hfx_basis_info_type)                :: basis_info
     TYPE(hfx_screen_coeff_type), &
          DIMENSION(:,:,:,:,:,:), POINTER     :: screen_funct_coeffs_pgf,&
                                                 pair_dist_radii_pgf
     TYPE(hfx_screen_coeff_type), &
          DIMENSION(:,:,:,:), POINTER         :: screen_funct_coeffs_set
     TYPE(hfx_screen_coeff_type), &
          DIMENSION(:,:), POINTER             :: screen_funct_coeffs_kind
     LOGICAL                                  :: screen_funct_is_initialized
     TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: initial_p
     TYPE(hfx_p_kind), DIMENSION(:), POINTER  :: initial_p_forces
     INTEGER, DIMENSION(:), POINTER           :: map_atom_to_kind_atom
     TYPE(hfx_2D_map), DIMENSION(:), POINTER  :: map_atoms_to_cpus
     INTEGER, DIMENSION(:,:), POINTER         :: atomic_block_offset
     INTEGER, DIMENSION(:,:,:,:), POINTER     :: set_offset
     INTEGER, DIMENSION(:), POINTER           :: block_offset
     TYPE(hfx_block_range_type), DIMENSION(:), &
                                 POINTER      :: blocks
     TYPE(hfx_task_list_type), DIMENSION(:), &
                               POINTER        :: task_list
     REAL(dp), DIMENSION(:,:), POINTER        :: pmax_atom, pmax_atom_forces
     TYPE(lib_int)                            :: lib_int
     TYPE(lib_deriv)                          :: lib_deriv
     REAL(dp), DIMENSION(:,:), POINTER        :: pmax_block
     LOGICAL, DIMENSION(:,:), POINTER         :: atomic_pair_list
     LOGICAL, DIMENSION(:,:), POINTER         :: atomic_pair_list_forces
  END TYPE hfx_type

  CONTAINS

! *****************************************************************************
!> \brief - This routine allocates and initializes all types in hfx_data
!> \param x_data contains all relevant data structures for hfx runs
!> \param para_env ...
!> \param hfx_section input section
!> \param natom , max_set: quantities that define some array sizes
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param dft_control ...
!> \param cell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param do_exx ...
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
!> \note
!>      - All POINTERS and ALLOCATABLES are allocated, even if their size is
!>        unknown at invocation time
! *****************************************************************************
  SUBROUTINE hfx_create(x_data,para_env, hfx_section,natom,atomic_kind_set, qs_kind_set, &
                        dft_control, cell, error, do_exx)
    TYPE(hfx_type), DIMENSION(:, :), POINTER :: x_data
    TYPE(cp_para_env_type)                   :: para_env
    TYPE(section_vals_type), POINTER         :: hfx_section
    INTEGER, INTENT(IN)                      :: natom
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, OPTIONAL                        :: do_exx

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=512)                       :: error_msg
    CHARACTER(LEN=default_path_length)       :: char_val
    INTEGER :: handle, i, i_thread, iatom, ikind, int_val, irep, jkind, &
      max_set, n_rep_hf, n_threads, natom_a, natom_b, nkind, nseta, nsetb, &
      pbc_shells, stat, storage_id
    INTEGER, DIMENSION(:), POINTER           :: atom2kind, kind_of
    LOGICAL                                  :: failure, logic_val, my_do_exx
    REAL(dp)                                 :: real_val
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(section_vals_type), POINTER         :: hf_pbc_section, hf_sub_section

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CALL cite_reference(Guidon2008)
    CALL cite_reference(Guidon2009)
    !! There might be 2 hf sections
    CALL section_vals_get(hfx_section,n_repetition=n_rep_hf,error=error)
    n_threads = 1
!$  n_threads = omp_get_max_threads()
    ALLOCATE(x_data(n_rep_hf,n_threads),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i_thread=1,n_threads
      DO irep=1,n_rep_hf
        actual_x_data => x_data(irep,i_thread)
        !! Get data from input file
        !!
        !! GENERAL params
        CALL section_vals_val_get(hfx_section, "FRACTION", r_val=real_val, i_rep_section=irep,error=error)
        actual_x_data%general_parameter%fraction = real_val
        actual_x_data%n_rep_hf = n_rep_hf

        NULLIFY(actual_x_data%map_atoms_to_cpus)

        CALL section_vals_val_get(hfx_section, "TREAT_LSD_IN_CORE", l_val=logic_val, i_rep_section=irep,error=error)
        actual_x_data%general_parameter%treat_lsd_in_core = logic_val

        !! MEMORY section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"MEMORY",i_rep_section=irep, error=error)
        CALL parse_memory_section(actual_x_data%memory_parameter, hf_sub_section, storage_id, i_thread,&
             n_threads, para_env, irep, skip_disk=.FALSE., skip_in_core_forces=.FALSE., error=error)

        !! PERIODIC section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"PERIODIC", i_rep_section=irep, error=error)
        CALL section_vals_val_get(hf_sub_section, "NUMBER_OF_SHELLS", i_val=int_val, error=error)
        actual_x_data%periodic_parameter%number_of_shells = int_val
        actual_x_data%periodic_parameter%mode = int_val
        CALL get_cell(cell=cell,periodic=actual_x_data%periodic_parameter%perd)
        IF( SUM(actual_x_data%periodic_parameter%perd) == 0 ) THEN
          actual_x_data%periodic_parameter%do_periodic = .FALSE.
        ELSE
          actual_x_data%periodic_parameter%do_periodic = .TRUE.
        END IF

        !! SCREENING section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"SCREENING",i_rep_section=irep, error=error)
        CALL section_vals_val_get(hf_sub_section, "EPS_SCHWARZ", r_val=real_val, error=error)
        actual_x_data%screening_parameter%eps_schwarz = real_val
        CALL section_vals_val_get(hf_sub_section, "EPS_SCHWARZ_FORCES", r_val=real_val, error=error)
        actual_x_data%screening_parameter%eps_schwarz_forces = real_val
        CALL section_vals_val_get(hf_sub_section, "SCREEN_P_FORCES", l_val=logic_val, error=error)
        actual_x_data%screening_parameter%do_p_screening_forces = logic_val
        CALL section_vals_val_get(hf_sub_section, "SCREEN_ON_INITIAL_P", l_val=logic_val, error=error)
        actual_x_data%screening_parameter%do_initial_p_screening = logic_val
        actual_x_data%screen_funct_is_initialized = .FALSE.

        !! INTERACTION_POTENTIAL section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"INTERACTION_POTENTIAL",i_rep_section=irep,error=error)
        CALL section_vals_val_get(hf_sub_section, "POTENTIAL_TYPE", i_val=int_val, error=error)
        actual_x_data%potential_parameter%potential_type = int_val
        CALL section_vals_val_get(hf_sub_section, "OMEGA", r_val=real_val, error=error)
        actual_x_data%potential_parameter%omega = real_val
        CALL section_vals_val_get(hf_sub_section, "SCALE_COULOMB", r_val=real_val, error=error)
        actual_x_data%potential_parameter%scale_coulomb = real_val
        CALL section_vals_val_get(hf_sub_section, "SCALE_LONGRANGE", r_val=real_val, error=error)
        actual_x_data%potential_parameter%scale_longrange = real_val
        CALL section_vals_val_get(hf_sub_section, "SCALE_GAUSSIAN", r_val=real_val, error=error)
        actual_x_data%potential_parameter%scale_gaussian = real_val
        IF (actual_x_data%potential_parameter%potential_type == do_hfx_potential_truncated .OR. &
            actual_x_data%potential_parameter%potential_type == do_hfx_potential_mix_cl_trunc) THEN
          CALL section_vals_val_get(hf_sub_section, "CUTOFF_RADIUS", r_val=real_val, error=error)
          actual_x_data%potential_parameter%cutoff_radius = real_val
          CALL section_vals_val_get(hf_sub_section, "T_C_G_DATA", c_val=char_val, error=error)
          CALL compress(char_val,.TRUE.)
          ! ** Check if file is there
          IF ( .NOT. file_exists(char_val) ) THEN
            WRITE(error_msg,'(A,A,A)') "Truncated hfx calculation requested. The file containing "//&
                 "the data could not be found at ", TRIM(char_val), " Please check T_C_G_DATA "//&
                 "in the INTERCATION_POTENTIAL section"
            CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                 error_msg , &
                 error,&
                 only_ionode=.TRUE.)
          ELSE
            actual_x_data%potential_parameter%filename = char_val
          END IF
        END IF
        IF (actual_x_data%potential_parameter%potential_type == do_hfx_potential_short) THEN
            CALL erfc_cutoff(actual_x_data%screening_parameter%eps_schwarz,&
                           actual_x_data%potential_parameter%omega,&
                           actual_x_data%potential_parameter%cutoff_radius, error)
        END IF

        !! LOAD_BALANCE section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"LOAD_BALANCE",i_rep_section=irep, error=error)
        CALL section_vals_val_get(hf_sub_section, "NBINS", i_val=int_val, error=error)
        actual_x_data%load_balance_parameter%nbins = MAX(1,int_val)
        actual_x_data%load_balance_parameter%blocks_initialized = .FALSE.

        CALL section_vals_val_get(hf_sub_section, "RANDOMIZE", l_val=logic_val, error=error)
        actual_x_data%load_balance_parameter%do_randomize = logic_val

        actual_x_data%load_balance_parameter%rtp_redistribute=.FALSE.
        IF(ASSOCIATED(dft_control%rtp_control))&
            actual_x_data%load_balance_parameter%rtp_redistribute=dft_control%rtp_control%hfx_redistribute

        CALL section_vals_val_get(hf_sub_section, "BLOCK_SIZE", i_val=int_val, error=error)
        ! negative values ask for a computed default
        IF (int_val<=0) THEN
           ! this gives a reasonable number of blocks for binning, yet typically results in blocking.
           int_val = CEILING( 0.1_dp * natom / &
                 REAL(actual_x_data%load_balance_parameter%nbins * n_threads * para_env%num_pe ,KIND=dp)**(0.25_dp))
        ENDIF
        ! at least 1 atom per block, and avoid overly large blocks
        actual_x_data%load_balance_parameter%block_size = MIN(max_atom_block, MAX(1, int_val))

        my_do_exx=.FALSE.
        IF(PRESENT(do_exx)) THEN
          my_do_exx=do_exx
        END IF

        ! in the case of RPA exchange (postprocessing after SCF), we do not use ADMM
        IF(my_do_exx) THEN
          CALL hfx_create_basis_types(actual_x_data%basis_parameter,actual_x_data%basis_info, qs_kind_set,&
                          .FALSE., error)
        ELSE
          CALL hfx_create_basis_types(actual_x_data%basis_parameter,actual_x_data%basis_info, qs_kind_set,&
                          dft_control%do_admm, error) 
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
!! ****************************************************************************************************

        IF(actual_x_data%periodic_parameter%do_periodic) THEN
          hf_pbc_section => section_vals_get_subs_vals(hfx_section,"PERIODIC",i_rep_section=irep, error=error)
          CALL section_vals_val_get(hf_pbc_section,"NUMBER_OF_SHELLS",i_val=pbc_shells,error=error)
          actual_x_data%periodic_parameter%number_of_shells_from_input = pbc_shells
          ALLOCATE(actual_x_data%neighbor_cells(1),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL hfx_create_neighbor_cells(actual_x_data,pbc_shells, cell, i_thread, error=error)
        ELSE
          ALLOCATE(actual_x_data%neighbor_cells(1),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! ** Initialize this guy to enable non periodic stress regtests
          actual_x_data%periodic_parameter%R_max_stress = 1.0_dp
        END IF

        nkind = SIZE(qs_kind_set,1)
        max_set=actual_x_data%basis_info%max_set

        !! ** This guy is allocated on the master thread only
        IF (i_thread == 1 ) THEN
          ALLOCATE(actual_x_data%is_assoc_atomic_block(natom,natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(actual_x_data%atomic_block_offset(natom,natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(actual_x_data%set_offset(max_set, max_set, nkind, nkind),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(actual_x_data%block_offset(para_env%num_pe+1),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

        ALLOCATE(actual_x_data%distribution_forces(1),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(actual_x_data%distribution_energy(1),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        actual_x_data%memory_parameter%size_p_screen = 0_int_8
        IF( i_thread == 1 ) THEN
          ALLOCATE(actual_x_data%atomic_pair_list(natom,natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(actual_x_data%atomic_pair_list_forces(natom,natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF(actual_x_data%screening_parameter%do_initial_p_screening .OR. &
           actual_x_data%screening_parameter%do_p_screening_forces) THEN
          !! ** This guy is allocated on the master thread only
          IF (i_thread == 1 ) THEN
            ALLOCATE(actual_x_data%pmax_atom(natom,natom),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            ALLOCATE(actual_x_data%initial_p(nkind*(nkind+1)/2),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            i = 1
            DO ikind=1,nkind
              CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_a)
              nseta=actual_x_data%basis_parameter(ikind)%nset
              DO jkind=ikind,nkind
                CALL get_atomic_kind(atomic_kind_set(jkind), natom=natom_b)
                nsetb=actual_x_data%basis_parameter(jkind)%nset
                ALLOCATE(actual_x_data%initial_p(i)%p_kind(nseta, nsetb, natom_a, natom_b),STAT=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                actual_x_data%memory_parameter%size_p_screen = &
                  actual_x_data%memory_parameter%size_p_screen + nseta*nsetb*natom_a*natom_b
                i = i + 1
              END DO
            END DO
            IF( actual_x_data%memory_parameter%treat_forces_in_core ) THEN
              ALLOCATE(actual_x_data%pmax_atom_forces(natom,natom),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

              ALLOCATE(actual_x_data%initial_p_forces(nkind*(nkind+1)/2),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              i = 1
              DO ikind=1,nkind
                CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_a)
                nseta=actual_x_data%basis_parameter(ikind)%nset
                DO jkind=ikind,nkind
                  CALL get_atomic_kind(atomic_kind_set(jkind), natom=natom_b)
                  nsetb=actual_x_data%basis_parameter(jkind)%nset
                  ALLOCATE(actual_x_data%initial_p_forces(i)%p_kind(nseta, nsetb, natom_a, natom_b),STAT=stat)
                  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                  actual_x_data%memory_parameter%size_p_screen = &
                  actual_x_data%memory_parameter%size_p_screen + nseta*nsetb*natom_a*natom_b
                  i = i + 1
                END DO
              END DO
            END IF
          END IF
          ALLOCATE(actual_x_data%map_atom_to_kind_atom(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(kind_of(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (atom2kind(nkind),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

          atom2kind = 0
          DO iatom=1,natom
            ikind = kind_of(iatom)
            atom2kind(ikind) = atom2kind(ikind) + 1
            actual_x_data%map_atom_to_kind_atom(iatom) = atom2kind(ikind)
          END DO
          DEALLOCATE(kind_of, atom2kind, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

        ! ** Initialize libint type
        CALL initialize_libint(actual_x_data%lib, actual_x_data%basis_info%max_am, error)
        CALL initialize_libderiv(actual_x_data%lib_deriv, actual_x_data%basis_info%max_am, error)

        CALL alloc_containers(actual_x_data, 1, hfx_do_eval_energy, error)
        CALL alloc_containers(actual_x_data, 1, hfx_do_eval_forces, error)

        actual_x_data%maxval_cache_disk%element_counter = 1
        ALLOCATE(actual_x_data%maxval_container_disk,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(actual_x_data%maxval_container_disk%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        actual_x_data%maxval_container_disk%first%prev => NULL()
        actual_x_data%maxval_container_disk%first%next => NULL()
        actual_x_data%maxval_container_disk%current => actual_x_data%maxval_container_disk%first
        actual_x_data%maxval_container_disk%current%data = 0
        actual_x_data%maxval_container_disk%element_counter = 1
        actual_x_data%maxval_container_disk%file_counter = 1
        actual_x_data%maxval_container_disk%desc = 'Max_'
        actual_x_data%maxval_container_disk%unit = -1
        WRITE(actual_x_data%maxval_container_disk%filename,'(A,I0,A,A,A)') TRIM(actual_x_data%memory_parameter%storage_location), &
                                                           storage_id,"_",actual_x_data%maxval_container_disk%desc, "6"
        CALL compress(actual_x_data%maxval_container_disk%filename, .TRUE.)
        ALLOCATE(actual_x_data%integral_containers_disk(64), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO i=1,64
          actual_x_data%integral_caches_disk(i)%element_counter = 1
          actual_x_data%integral_caches_disk(i)%data = 0
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(actual_x_data%integral_containers_disk(i)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          actual_x_data%integral_containers_disk(i)%first%prev => NULL()
          actual_x_data%integral_containers_disk(i)%first%next => NULL()
          actual_x_data%integral_containers_disk(i)%current => actual_x_data%integral_containers_disk(i)%first
          actual_x_data%integral_containers_disk(i)%current%data = 0
          actual_x_data%integral_containers_disk(i)%element_counter = 1
          actual_x_data%integral_containers_disk(i)%file_counter = 1
          actual_x_data%integral_containers_disk(i)%desc = 'Int_'
          actual_x_data%integral_containers_disk(i)%unit = -1
          WRITE(actual_x_data%integral_containers_disk(i)%filename,'(A,I0,A,A,I0)')&
                                             TRIM(actual_x_data%memory_parameter%storage_location), &
                                             storage_id,"_",actual_x_data%integral_containers_disk(i)%desc, i
          CALL compress(actual_x_data%integral_containers_disk(i)%filename, .TRUE.)
        END DO

        actual_x_data%b_first_load_balance_energy = .TRUE.
        actual_x_data%b_first_load_balance_forces = .TRUE.
      END DO
    END DO

    DO irep=1,n_rep_hf
      actual_x_data => x_data(irep,1)
      CALL hfx_print_info(actual_x_data, hfx_section, irep, error=error)
    END DO
    CALL timestop(handle)
  END SUBROUTINE hfx_create

! *****************************************************************************
!> \brief - compute a truncation radius for the shortrange operator
!> \param eps target accuracy
!> \param omg screening parameter
!> \param r_cutoff cutoff radius 
!> \param error ...
!> \par History
!>      10.2012 created [Hossein Banihashemian]
!> \author Hossein Banihashemian
! *****************************************************************************
  SUBROUTINE erfc_cutoff(eps,omg,r_cutoff,error)
    IMPLICIT NONE

    REAL(dp), INTENT(in)  :: eps, omg
    REAL(dp), INTENT(out) :: r_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'erfc_cutoff', &
      routineP = moduleN//':'//routineN

    REAL(dp), PARAMETER :: abstol=1E-10_dp, soltol=1E-16_dp
    REAL(dp) :: r0, f0, fprime0, delta_r
    INTEGER :: iter,handle
    INTEGER, PARAMETER :: iterMAX = 1000
    LOGICAL :: failure

    CALL timeset(routineN,handle)

    ! initial guess assuming that we are in the asymptotic regime of the erf, and the solution is about 10.
    r0 = SQRT(-LOG(eps*omg*10**2))/omg
    CALL eval_transc_func(r0,eps,omg,f0,fprime0)

    DO iter = 1,iterMAX
        delta_r = f0/fprime0
        r0 = r0-delta_r
        CALL eval_transc_func(r0,eps,omg,f0,fprime0)
        IF (ABS(delta_r) .LT. abstol .OR. ABS(f0) .LT. soltol) EXIT
    END DO
    CPPostcondition(iter<=itermax,cp_failure_level,routineP,error,failure)
    r_cutoff = r0

    CALL timestop(handle)
  CONTAINS
! *****************************************************************************
!> \brief ...
!> \param r ...
!> \param eps ...
!> \param omega ...
!> \param fn ...
!> \param df ...
! *****************************************************************************
    SUBROUTINE eval_transc_func(r,eps,omega,fn,df)
    REAL(dp), INTENT(in)                     :: r, eps, omega
    REAL(dp), INTENT(out)                    :: fn, df

      fn = erfc(omega*r) - r*eps
      df = -omega*2*EXP(-omega**2*r**2)/SQRT(pi) - eps
    END SUBROUTINE eval_transc_func
  END SUBROUTINE erfc_cutoff
! *****************************************************************************
!> \brief - This routine allocates and initializes the basis_info and basis_parameter types
!> \param basis_parameter ...
!> \param basis_info ...
!> \param qs_kind_set ...
!> \param do_admm ...
!> \param error ...
!> \par History
!>      07.2011 refactored
! *****************************************************************************
  SUBROUTINE hfx_create_basis_types(basis_parameter,basis_info, qs_kind_set, do_admm, error)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_basis_info_type)                :: basis_info
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    LOGICAL                                  :: do_admm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create_basis_types', &
      routineP = moduleN//':'//routineN

    INTEGER :: co_counter, handle, i, ikind, ipgf, iset, j, k, la, &
      max_am_kind, max_coeff, max_nsgfl, max_pgf, max_pgf_kind, max_set, &
      nkind, nl_count, nset, nseta, offset_a, offset_a1, s_offset_nl_a, sgfa, &
      so_counter, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, la_min, npgfa, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, nl_a
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: sphi_a
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a

    CALL timeset(routineN,handle)

    failure = .FALSE.

    !! BASIS parameter
    nkind = SIZE(qs_kind_set,1)
    ALLOCATE(basis_parameter(nkind), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    max_set = 0
    DO ikind = 1,nkind
        IF( .NOT. do_admm ) THEN
          CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_a, basis_type="ORB", error=error)
          CALL get_qs_kind_set(qs_kind_set,&
                               maxsgf=basis_info%max_sgf,&
                               maxnset=basis_info%max_set,&
                               maxlgto=basis_info%max_am, error=error)

        ELSE
          CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_a, basis_type="AUX_FIT", error=error)
          CALL get_qs_kind_set(qs_kind_set,&
                               maxsgf=basis_info%max_sgf,&
                               maxnset=basis_info%max_set,&
                               maxlgto=basis_info%max_am,&
                               basis_type="AUX_FIT", error=error)
        END IF
        IF (basis_info%max_set<max_set) STOP "UNEXPECTED MAX_SET"
        max_set = MAX(max_set,basis_info%max_set)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_a, &
                               lmax=basis_parameter(ikind)%lmax, &
                               lmin=basis_parameter(ikind)%lmin, &
                               npgf=basis_parameter(ikind)%npgf, &
                               nset=basis_parameter(ikind)%nset, &
                               zet=basis_parameter(ikind)%zet, &
                               nsgf_set=basis_parameter(ikind)%nsgf, &
                               first_sgf=basis_parameter(ikind)%first_sgf, &
                               sphi=basis_parameter(ikind)%sphi, &
                               nsgf=basis_parameter(ikind)%nsgf_total,&
                               l=basis_parameter(ikind)%nl,&
                               nshell=basis_parameter(ikind)%nshell,&
                               set_radius=basis_parameter(ikind)%set_radius,&
                               pgf_radius=basis_parameter(ikind)%pgf_radius,&
                               kind_radius=basis_parameter(ikind)%kind_radius)
    END DO
    DO ikind = 1,nkind
      ALLOCATE(basis_parameter(ikind)%nsgfl(0:basis_info%max_am,max_set))
      basis_parameter(ikind)%nsgfl = 0
      nset = basis_parameter(ikind)%nset
      nshell => basis_parameter(ikind)%nshell
      DO iset = 1,nset
        DO i=0,basis_info%max_am
          nl_count = 0
          DO j=1,nshell(iset)
            IF( basis_parameter(ikind)%nl(j,iset) == i) nl_count = nl_count + 1
          END DO
          basis_parameter(ikind)%nsgfl(i,iset) = nl_count
        END DO
      END DO
    END DO

    max_nsgfl = 0
    max_pgf = 0
    DO ikind = 1,nkind
      max_coeff = 0
      max_am_kind = 0
      max_pgf_kind = 0
      npgfa => basis_parameter(ikind)%npgf
      nseta = basis_parameter(ikind)%nset
      nl_a => basis_parameter(ikind)%nsgfl
      la_max => basis_parameter(ikind)%lmax
      la_min => basis_parameter(ikind)%lmin
      DO iset = 1,nseta
        max_pgf_kind = MAX(max_pgf_kind,npgfa(iset))
        max_pgf = MAX(max_pgf, npgfa(iset))
        DO la = la_min(iset), la_max(iset)
          max_nsgfl = MAX(max_nsgfl, nl_a(la,iset))
          max_coeff = MAX(max_coeff, nso(la)*nl_a(la,iset)*nco(la))
          max_am_kind = MAX(max_am_kind, la)
        END DO
      END DO
      ALLOCATE(basis_parameter(ikind)%sphi_ext(max_coeff,0:max_am_kind, max_pgf_kind, nseta), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      basis_parameter(ikind)%sphi_ext = 0.0_dp
    END DO

    DO ikind = 1,nkind
      sphi_a => basis_parameter(ikind)%sphi
      nseta = basis_parameter(ikind)%nset
      la_max => basis_parameter(ikind)%lmax
      la_min => basis_parameter(ikind)%lmin
      npgfa => basis_parameter(ikind)%npgf
      first_sgfa => basis_parameter(ikind)%first_sgf
      nl_a => basis_parameter(ikind)%nsgfl
      DO iset = 1,nseta
        sgfa = first_sgfa(1,iset)
        DO ipgf = 1,npgfa(iset)
          offset_a1 = (ipgf-1)*ncoset(la_max(iset))
          s_offset_nl_a = 0
          DO la = la_min(iset),la_max(iset)
            offset_a = offset_a1 + ncoset(la-1)
            co_counter = 0
            co_counter = co_counter + 1
            so_counter = 0
            DO k = sgfa+s_offset_nl_a,sgfa+s_offset_nl_a+nso(la)*nl_a(la,iset)-1
              DO i = offset_a+1,offset_a+nco(la)
                so_counter = so_counter + 1
                basis_parameter(ikind)%sphi_ext(so_counter,la,ipgf,iset) = sphi_a(i,k)
              END DO
            END DO
            s_offset_nl_a = s_offset_nl_a + nso(la) * (nl_a(la,iset))
          END DO
        END DO
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE hfx_create_basis_types


! *****************************************************************************
!> \brief ...
!> \param basis_parameter ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE hfx_release_basis_types(basis_parameter,error)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_release_basis_types', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure = .FALSE.

    !! BASIS parameter
    DO i=1,SIZE(basis_parameter)
       DEALLOCATE(basis_parameter(i)%nsgfl,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(basis_parameter(i)%sphi_ext,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END DO
    DEALLOCATE(basis_parameter,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)

  END SUBROUTINE hfx_release_basis_types


! *****************************************************************************
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
!> \param error ...
! *****************************************************************************
  SUBROUTINE parse_memory_section(memory_parameter, hf_sub_section, storage_id,&
       i_thread, n_threads, para_env, irep, skip_disk, skip_in_core_forces, error)
    TYPE(hfx_memory_type)                    :: memory_parameter
    TYPE(section_vals_type), POINTER         :: hf_sub_section
    INTEGER, INTENT(OUT), OPTIONAL           :: storage_id
    INTEGER, INTENT(IN), OPTIONAL            :: i_thread, n_threads
    TYPE(cp_para_env_type), OPTIONAL         :: para_env
    INTEGER, INTENT(IN), OPTIONAL            :: irep
    LOGICAL, INTENT(IN)                      :: skip_disk, skip_in_core_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'parse_memory_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=512)                       :: error_msg
    CHARACTER(LEN=default_path_length)       :: char_val, filename, orig_wd
    INTEGER                                  :: int_val, stat
    LOGICAL                                  :: check, failure, logic_val
    REAL(dp)                                 :: real_val

    failure = .FALSE.
    check = (PRESENT(storage_id).EQV.PRESENT(i_thread)) .AND.&
            (PRESENT(storage_id).EQV.PRESENT(n_threads)).AND.&
            (PRESENT(storage_id).EQV.PRESENT(para_env)) .AND.&
            (PRESENT(storage_id).EQV.PRESENT(irep))
    CPPostcondition(check,cp_failure_level,routineP,error,failure)

    ! Memory Storage
    CALL section_vals_val_get(hf_sub_section, "MAX_MEMORY", i_val=int_val, error=error)
    memory_parameter%max_memory = int_val
    memory_parameter%max_compression_counter = int_val*1024_int_8*128_int_8
    CALL section_vals_val_get(hf_sub_section, "EPS_STORAGE", r_val=real_val, error=error)
    memory_parameter%eps_storage_scaling = real_val
    IF(int_val == 0) THEN
       memory_parameter%do_all_on_the_fly = .TRUE.
    ELSE
       memory_parameter%do_all_on_the_fly = .FALSE.
    END IF
    memory_parameter%cache_size = CACHE_SIZE
    memory_parameter%bits_max_val = BITS_MAX_VAL
    memory_parameter%actual_memory_usage = 1
    IF( .NOT. skip_in_core_forces) THEN
      CALL section_vals_val_get(hf_sub_section, "TREAT_FORCES_IN_CORE", l_val=logic_val, error=error)
      memory_parameter%treat_forces_in_core = logic_val
    END IF

    ! ** IF MAX_MEM == 0 overwrite this flag to false
    IF( memory_parameter%do_all_on_the_fly ) memory_parameter%treat_forces_in_core = .FALSE.

    ! Disk Storage
    IF (.NOT.skip_disk) THEN
       memory_parameter%actual_memory_usage_disk = 1
       CALL section_vals_val_get(hf_sub_section, "MAX_DISK_SPACE", i_val=int_val, error=error)
       memory_parameter%max_compression_counter_disk = int_val*1024_int_8*128_int_8
       IF( int_val == 0 ) THEN
          memory_parameter%do_disk_storage = .FALSE.
       ELSE
          memory_parameter%do_disk_storage = .TRUE.
       END IF
       CALL section_vals_val_get(hf_sub_section, "STORAGE_LOCATION", c_val=char_val, error=error)
       CALL compress(char_val,.TRUE.)
       !! Add ending / if necessary

       IF( SCAN(char_val,"/",.TRUE.) /= LEN_TRIM(char_val)) THEN
          WRITE(filename,'(A,A)') TRIM(char_val),"/"
          CALL compress(filename)
       ELSE
          filename = TRIM(char_val)
       END IF
       CALL compress(filename,.TRUE.)

       !! quickly check if we can write on storage_location
       CALL m_getcwd(orig_wd)
       CALL m_chdir(TRIM(filename),stat)
       IF(stat/=0) THEN
          WRITE(error_msg,'(A,A,A)') "Request for disk storage failed due to unknown error while writing to ",&
               TRIM(filename), ". Please check STORAGE_LOCATION"
          CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
               error_msg , &
               error,failure)
       END IF
       CALL m_chdir(orig_wd,stat)

       memory_parameter%storage_location = filename
       CALL compress(memory_parameter%storage_location, .TRUE.)
    ELSE
       memory_parameter%do_disk_storage = .FALSE.
    END IF
    IF (PRESENT(storage_id)) THEN
       storage_id = (irep-1)*para_env%num_pe*n_threads + para_env%mepos * n_threads + i_thread-1
    END IF
  END SUBROUTINE parse_memory_section

! *****************************************************************************
!> \brief - This routine deallocates all data structures
!> \param x_data contains all relevant data structures for hfx runs
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_release(x_data, error)
    TYPE(hfx_type), DIMENSION(:, :), POINTER :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i_thread, irep, n_rep_hf, &
                                                n_threads, stat
    LOGICAL                                  :: failure
    TYPE(hfx_type), POINTER                  :: actual_x_data

    failure = .FALSE.
    !! There might be 2 hf sections
    n_rep_hf = x_data(1,1)%n_rep_hf
    n_threads = 1
!$  n_threads = omp_get_max_threads()

    IF(x_data(1,1)%potential_parameter%potential_type == do_hfx_potential_truncated .OR. &
       x_data(1,1)%potential_parameter%potential_type == do_hfx_potential_mix_cl_trunc ) THEN
      init_t_c_g0_lmax = -1
      CALL free_C0()
    END IF
    DO i_thread= 1,n_threads
      DO irep=1,n_rep_hf
        actual_x_data => x_data(irep,i_thread)
        DEALLOCATE(actual_x_data%neighbor_cells,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE(actual_x_data%distribution_energy,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE(actual_x_data%distribution_forces,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        IF( actual_x_data%load_balance_parameter%blocks_initialized ) THEN
          DEALLOCATE(actual_x_data%blocks,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF( i_thread == 1 ) THEN
            DEALLOCATE(actual_x_data%pmax_block,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
        END IF

        IF( i_thread == 1 ) THEN
          DEALLOCATE(actual_x_data%atomic_pair_list, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(actual_x_data%atomic_pair_list_forces, STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF(actual_x_data%screening_parameter%do_initial_p_screening .OR. &
           actual_x_data%screening_parameter%do_p_screening_forces) THEN
          IF( i_thread == 1 ) THEN
            DEALLOCATE(actual_x_data%pmax_atom, STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            DO i=1,SIZE(actual_x_data%initial_p)
               DEALLOCATE(actual_x_data%initial_p(i)%p_kind,STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            END DO
            DEALLOCATE(actual_x_data%initial_p,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            IF( actual_x_data%memory_parameter%treat_forces_in_core) THEN
              DEALLOCATE(actual_x_data%pmax_atom_forces, STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              DO i=1,SIZE(actual_x_data%initial_p_forces)
                 DEALLOCATE(actual_x_data%initial_p_forces(i)%p_kind,STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              END DO
              DEALLOCATE(actual_x_data%initial_p_forces,STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            END IF
          END IF
          DEALLOCATE(actual_x_data%map_atom_to_kind_atom,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        IF( i_thread == 1 ) THEN
          DEALLOCATE(actual_x_data%is_assoc_atomic_block,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(actual_x_data%atomic_block_offset,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(actual_x_data%set_offset,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(actual_x_data%block_offset,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

        !! BASIS parameter
        CALL hfx_release_basis_types(actual_x_data%basis_parameter,error)

        !MK Release libint and libderiv data structure
        CALL terminate_libint(actual_x_data%lib)
        CALL terminate_libderiv(actual_x_data%lib_deriv)

        !! Deallocate containers
        CALL dealloc_containers(actual_x_data, hfx_do_eval_energy, error)
        CALL dealloc_containers(actual_x_data, hfx_do_eval_forces, error)

        !! Deallocate containers
        CALL hfx_init_container(actual_x_data%maxval_container_disk, &
                                actual_x_data%memory_parameter%actual_memory_usage_disk, &
                                .FALSE., error)
        IF( actual_x_data%memory_parameter%do_disk_storage) THEN
          CALL close_file(unit_number=actual_x_data%maxval_container_disk%unit,file_status="DELETE")
        END IF
        DEALLOCATE(actual_x_data%maxval_container_disk%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE(actual_x_data%maxval_container_disk,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        DO i=1,64
          CALL hfx_init_container(actual_x_data%integral_containers_disk(i), &
                                  actual_x_data%memory_parameter%actual_memory_usage_disk, &
                                  .FALSE., error)
          IF( actual_x_data%memory_parameter%do_disk_storage) THEN
            CALL close_file(unit_number=actual_x_data%integral_containers_disk(i)%unit,file_status="DELETE")
          END IF
          DEALLOCATE(actual_x_data%integral_containers_disk(i)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END DO
        DEALLOCATE(actual_x_data%integral_containers_disk,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)


        ! ** screening functions
        IF(actual_x_data%screen_funct_is_initialized) THEN
          DEALLOCATE(actual_x_data%screen_funct_coeffs_set,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(actual_x_data%screen_funct_coeffs_kind,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(actual_x_data%pair_dist_radii_pgf,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(actual_x_data%screen_funct_coeffs_pgf,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          actual_x_data%screen_funct_is_initialized = .FALSE.
        END IF

        ! ** maps
        IF( ASSOCIATED(actual_x_data%map_atoms_to_cpus) ) THEN
          DO i=1,SIZE(actual_x_data%map_atoms_to_cpus)
            DEALLOCATE(actual_x_data%map_atoms_to_cpus(i)%iatom_list,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            DEALLOCATE(actual_x_data%map_atoms_to_cpus(i)%jatom_list,STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END DO
          DEALLOCATE(actual_x_data%map_atoms_to_cpus,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
      END DO

    END DO
    DEALLOCATE(x_data,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE hfx_release

! *****************************************************************************
!> \brief - This routine computes the neighbor cells that are taken into account
!>        in periodic runs
!> \param x_data contains all relevant data structures for hfx runs
!> \param pbc_shells number of shells taken into account
!> \param cell cell
!> \param i_thread current thread ID
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_create_neighbor_cells(x_data, pbc_shells, cell, i_thread, &
                                       error)
    TYPE(hfx_type), POINTER                  :: x_data
    INTEGER, INTENT(INOUT)                   :: pbc_shells
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: i_thread
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create_neighbor_cells', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=512)                       :: error_msg
    CHARACTER(LEN=64)                        :: char_nshells
    INTEGER :: i, idx, ikind, ipgf, iset, ishell, j, jkind, jpgf, jset, &
      jshell, k, kshell, l, m(3), max_shell, nseta, nsetb, stat, &
      total_number_of_cells, ub, ub_max
    INTEGER, DIMENSION(:), POINTER           :: la_max, lb_max, npgfa, npgfb
    LOGICAL                                  :: failure, image_cell_found, &
                                                nothing_more_to_add
    REAL(dp) :: cross_product(3), dist_min, distance(14), l_min, normal(3,6), &
      P(3,14), plane_vector(3,2), point_in_plane(3), r(3), R1, R_max, &
      R_max_stress, s(3), x, y, z, Zeta1
    REAL(dp), DIMENSION(:, :), POINTER       :: zeta, zetb
    TYPE(hfx_cell_type), ALLOCATABLE, &
      DIMENSION(:)                           :: tmp_neighbor_cells

    failure = .FALSE.
    total_number_of_cells = 0

    ! ** Check some settings
    IF( i_thread == 1 ) THEN
      IF( x_data%potential_parameter%potential_type /= do_hfx_potential_truncated .AND. &
          x_data%potential_parameter%potential_type /= do_hfx_potential_short .AND. &
          x_data%potential_parameter%potential_type /= do_hfx_potential_mix_cl_trunc) THEN
        CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                   "Periodic Hartree Fock calculation requested without use "//&
                   "of a truncated or shortrange potential. This may lead to unphysical total energies. "//&
                   "Use a truncated  potential to avoid possible problems. "//&
                   CPSourceFileRef,&
                   only_ionode=.TRUE.)
      ELSE
        l_min = MIN( plane_distance(1,0,0,cell), &
                     plane_distance(0,1,0,cell), &
                     plane_distance(0,0,1,cell) )
        l_min = 0.5_dp * l_min
        IF( x_data%potential_parameter%cutoff_radius >= l_min ) THEN
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                   "Periodic Hartree Fock calculation requested with use "//&
                   "of a truncated or shortrange potential. The cutoff radius is larger than half "//&
                   "the minimal cell dimension. This may lead to unphysical "//&
                   "total energies. Reduce the cutoff radius in order to avoid "//&
                   "possible problems. "//&
                   CPSourceFileRef,&
                   only_ionode=.TRUE.)
        END IF
      END IF
    END IF

    SELECT CASE(x_data%potential_parameter%potential_type)
    CASE( do_hfx_potential_truncated, do_hfx_potential_mix_cl_trunc, do_hfx_potential_short )
      R_max = 0.0_dp
      DO ikind=1,SIZE(x_data%basis_parameter)
        la_max => x_data%basis_parameter(ikind)%lmax
        zeta => x_data%basis_parameter(ikind)%zet
        nseta = x_data%basis_parameter(ikind)%nset
        npgfa => x_data%basis_parameter(ikind)%npgf
        DO jkind=1,SIZE(x_data%basis_parameter)
          lb_max => x_data%basis_parameter(jkind)%lmax
          zetb => x_data%basis_parameter(jkind)%zet
          nsetb = x_data%basis_parameter(jkind)%nset
          npgfb => x_data%basis_parameter(jkind)%npgf
          DO iset=1,nseta
            DO jset=1,nsetb
              DO ipgf=1,npgfa(iset)
                DO jpgf=1,npgfb(jset)
                  Zeta1 = zeta(ipgf,iset) + zetb(jpgf,jset)
                  R1 = 1.0_dp/SQRT(Zeta1) * mul_fact(la_max(iset)+lb_max(jset)) * &
                       SQRT(-LOG(x_data%screening_parameter%eps_schwarz))
                  R_max = MAX(R1,R_max)
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
      DEALLOCATE(x_data%neighbor_cells,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(x_data%neighbor_cells(1), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
        ALLOCATE(tmp_neighbor_cells(1:ub),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO i = 1,ub-1
          tmp_neighbor_cells(i) = x_data%neighbor_cells(i)
        END DO
        ub_max = (2*max_shell+1)**3
        DEALLOCATE(x_data%neighbor_cells, STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(x_data%neighbor_cells(1:ub_max), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO i = 1,ub-1
          x_data%neighbor_cells(i) = tmp_neighbor_cells(i)
        END DO
        DO i = ub,ub_max
          x_data%neighbor_cells(i)%cell = 0.0_dp
          x_data%neighbor_cells(i)%cell_r = 0.0_dp
        END DO

        DEALLOCATE(tmp_neighbor_cells, STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        DO ishell=-max_shell,max_shell
        DO jshell=-max_shell,max_shell
        DO kshell=-max_shell,max_shell
          IF( MAX(ABS(ishell),ABS(jshell),ABS(kshell)) /= max_shell) CYCLE
          idx = 0
          DO j=0,1
            x = -1.0_dp/2.0_dp + j*1.0_dp
            DO k=0,1
              y = -1.0_dp/2.0_dp + k*1.0_dp
              DO l=0,1
                z = -1.0_dp/2.0_dp + l*1.0_dp
                idx = idx + 1
                P(1,idx) = x + ishell
                P(2,idx) = y + jshell
                P(3,idx) = z + kshell
                CALL scaled_to_real(r,P(:,idx),cell)
                distance(idx) = SQRT(SUM(r**2))
                P(1:3,idx)=r
              END DO
            END DO
          END DO
          ! Now check distance to Faces and only take them into account if the base point lies within quadrilateral

          ! Face A (1342) 1 is the reference
          idx = idx + 1
          plane_vector(:,1) = P(:,3)-P(:,1)
          plane_vector(:,2) = P(:,2)-P(:,1)
          cross_product(1) = plane_vector(2,1)*plane_vector(3,2) - plane_vector(3,1)*plane_vector(2,2)
          cross_product(2) = plane_vector(3,1)*plane_vector(1,2) - plane_vector(1,1)*plane_vector(3,2)
          cross_product(3) = plane_vector(1,1)*plane_vector(2,2) - plane_vector(2,1)*plane_vector(1,2)
          normal(:,1) = cross_product/SQRT(SUM(cross_product**2))
          point_in_plane = -normal(:,1)*(normal(1,1)*P(1,1)+normal(2,1)*P(2,1)+normal(3,1)*P(3,1))

          IF (point_is_in_quadrilateral(P(:,1),P(:,3),P(:,4),P(:,2),point_in_plane) ) THEN
            distance(idx)=ABS(normal(1,1)*P(1,1)+normal(2,1)*P(2,1)+normal(3,1)*P(3,1))
          ELSE
            distance(idx)=HUGE(distance(idx))
          END IF

          ! Face B (1562) 1 is the reference
          idx = idx + 1
          plane_vector(:,1) = P(:,2)-P(:,1)
          plane_vector(:,2) = P(:,5)-P(:,1)
          cross_product(1) = plane_vector(2,1)*plane_vector(3,2) - plane_vector(3,1)*plane_vector(2,2)
          cross_product(2) = plane_vector(3,1)*plane_vector(1,2) - plane_vector(1,1)*plane_vector(3,2)
          cross_product(3) = plane_vector(1,1)*plane_vector(2,2) - plane_vector(2,1)*plane_vector(1,2)
          normal(:,1) = cross_product/SQRT(SUM(cross_product**2))
          point_in_plane = -normal(:,1)*(normal(1,1)*P(1,1)+normal(2,1)*P(2,1)+normal(3,1)*P(3,1))

          IF (point_is_in_quadrilateral(P(:,1),P(:,5),P(:,6),P(:,2),point_in_plane) ) THEN
            distance(idx)=ABS(normal(1,1)*P(1,1)+normal(2,1)*P(2,1)+normal(3,1)*P(3,1))
          ELSE
            distance(idx)=HUGE(distance(idx))
          END IF

          ! Face C (5786) 5 is the reference
          idx = idx + 1
          plane_vector(:,1) = P(:,7)-P(:,5)
          plane_vector(:,2) = P(:,6)-P(:,5)
          cross_product(1) = plane_vector(2,1)*plane_vector(3,2) - plane_vector(3,1)*plane_vector(2,2)
          cross_product(2) = plane_vector(3,1)*plane_vector(1,2) - plane_vector(1,1)*plane_vector(3,2)
          cross_product(3) = plane_vector(1,1)*plane_vector(2,2) - plane_vector(2,1)*plane_vector(1,2)
          normal(:,1) = cross_product/SQRT(SUM(cross_product**2))
          point_in_plane = -normal(:,1)*(normal(1,1)*P(1,5)+normal(2,1)*P(2,5)+normal(3,1)*P(3,5))

          IF (point_is_in_quadrilateral(P(:,5),P(:,7),P(:,8),P(:,6),point_in_plane) ) THEN
            distance(idx)=ABS(normal(1,1)*P(1,5)+normal(2,1)*P(2,5)+normal(3,1)*P(3,5))
          ELSE
            distance(idx)=HUGE(distance(idx))
          END IF

          ! Face D (3784) 3 is the reference
          idx = idx + 1
          plane_vector(:,1) = P(:,7)-P(:,3)
          plane_vector(:,2) = P(:,4)-P(:,3)
          cross_product(1) = plane_vector(2,1)*plane_vector(3,2) - plane_vector(3,1)*plane_vector(2,2)
          cross_product(2) = plane_vector(3,1)*plane_vector(1,2) - plane_vector(1,1)*plane_vector(3,2)
          cross_product(3) = plane_vector(1,1)*plane_vector(2,2) - plane_vector(2,1)*plane_vector(1,2)
          normal(:,1) = cross_product/SQRT(SUM(cross_product**2))
          point_in_plane = -normal(:,1)*(normal(1,1)*P(1,3)+normal(2,1)*P(2,3)+normal(3,1)*P(3,3))

          IF (point_is_in_quadrilateral(P(:,3),P(:,7),P(:,8),P(:,4),point_in_plane) ) THEN
            distance(idx)=ABS(normal(1,1)*P(1,3)+normal(2,1)*P(2,3)+normal(3,1)*P(3,3))
          ELSE
            distance(idx)=HUGE(distance(idx))
          END IF

          ! Face E (2684) 2 is the reference
          idx = idx + 1
          plane_vector(:,1) = P(:,6)-P(:,2)
          plane_vector(:,2) = P(:,4)-P(:,2)
          cross_product(1) = plane_vector(2,1)*plane_vector(3,2) - plane_vector(3,1)*plane_vector(2,2)
          cross_product(2) = plane_vector(3,1)*plane_vector(1,2) - plane_vector(1,1)*plane_vector(3,2)
          cross_product(3) = plane_vector(1,1)*plane_vector(2,2) - plane_vector(2,1)*plane_vector(1,2)
          normal(:,1) = cross_product/SQRT(SUM(cross_product**2))
          point_in_plane = -normal(:,1)*(normal(1,1)*P(1,2)+normal(2,1)*P(2,2)+normal(3,1)*P(3,2))

          IF (point_is_in_quadrilateral(P(:,2),P(:,6),P(:,8),P(:,4),point_in_plane) ) THEN
            distance(idx)=ABS(normal(1,1)*P(1,2)+normal(2,1)*P(2,2)+normal(3,1)*P(3,2))
          ELSE
            distance(idx)=HUGE(distance(idx))
          END IF

          ! Face F (1573) 1 is the reference
          idx = idx + 1
          plane_vector(:,1) = P(:,5)-P(:,1)
          plane_vector(:,2) = P(:,3)-P(:,1)
          cross_product(1) = plane_vector(2,1)*plane_vector(3,2) - plane_vector(3,1)*plane_vector(2,2)
          cross_product(2) = plane_vector(3,1)*plane_vector(1,2) - plane_vector(1,1)*plane_vector(3,2)
          cross_product(3) = plane_vector(1,1)*plane_vector(2,2) - plane_vector(2,1)*plane_vector(1,2)
          normal(:,1) = cross_product/SQRT(SUM(cross_product**2))
          point_in_plane = -normal(:,1)*(normal(1,1)*P(1,1)+normal(2,1)*P(2,1)+normal(3,1)*P(3,1))

          IF (point_is_in_quadrilateral(P(:,1),P(:,5),P(:,7),P(:,3),point_in_plane) )THEN
            distance(idx)=ABS(normal(1,1)*P(1,1)+normal(2,1)*P(2,1)+normal(3,1)*P(3,1))
          ELSE
            distance(idx)=HUGE(distance(idx))
          END IF

          dist_min=MINVAL(distance)
          IF( max_shell == 0) THEN
            image_cell_found=.TRUE.
          END IF
          IF(dist_min<R_max)THEN
            total_number_of_cells = total_number_of_cells + 1
            x_data%neighbor_cells(ub)%cell = REAL((/ishell,jshell,kshell/),dp)
            ub = ub + 1
            image_cell_found = .TRUE.
          END IF

        END DO
        END DO
        END DO
        IF( image_cell_found) THEN
          max_shell = max_shell + 1
        ELSE
          nothing_more_to_add = .TRUE.
        END IF
      END DO
      ! now remove what is not needed
      ALLOCATE(tmp_neighbor_cells(total_number_of_cells), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO i=1,ub-1
        tmp_neighbor_cells(i) = x_data%neighbor_cells(i)
      END DO
      DEALLOCATE(x_data%neighbor_cells, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ! If we only need the supercell, total_number_of_cells is still 0, repair
      IF( total_number_of_cells == 0 ) THEN
        total_number_of_cells = 1
        ALLOCATE(x_data%neighbor_cells(total_number_of_cells), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO i=1,total_number_of_cells
          x_data%neighbor_cells(i)%cell = 0.0_dp
          x_data%neighbor_cells(i)%cell_r = 0.0_dp
        END DO
      ELSE
        ALLOCATE(x_data%neighbor_cells(total_number_of_cells), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO i=1,total_number_of_cells
          x_data%neighbor_cells(i) = tmp_neighbor_cells(i)
        END DO
      END IF
      DEALLOCATE(tmp_neighbor_cells, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      IF(x_data%periodic_parameter%number_of_shells == do_hfx_auto_shells ) THEN
        ! Do nothing
      ELSE
        total_number_of_cells = 0
        DO i = 0,x_data%periodic_parameter%number_of_shells
          total_number_of_cells = total_number_of_cells + count_cells_perd(i,x_data%periodic_parameter%perd)
        END DO
        IF( total_number_of_cells < SIZE(x_data%neighbor_cells)) THEN
          IF( i_thread == 1) THEN

            WRITE(char_nshells,'(I3)') SIZE(x_data%neighbor_cells)
            WRITE(error_msg,'(A,A,A)') "Periodic Hartree Fock calculation requested with use "//&
                     "of a truncated potential. The number of shells to be considered "//&
                     "might be too small. CP2K conservatively estimates to need "//TRIM(char_nshells)//" periodic images "//&
                     "Please carefully check if you get converged results."

            CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                           error_msg,&
                           error,&
                           only_ionode=.TRUE.)
          END IF
        END IF
        total_number_of_cells = 0
        DO i = 0,x_data%periodic_parameter%number_of_shells
          total_number_of_cells = total_number_of_cells + count_cells_perd(i,x_data%periodic_parameter%perd)
        END DO
        DEALLOCATE(x_data%neighbor_cells, STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        ALLOCATE(x_data%neighbor_cells(total_number_of_cells), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        m = 0
        i = 1
        DO WHILE(SUM(m**2)<=x_data%periodic_parameter%number_of_shells)
           x_data%neighbor_cells(i)%cell = REAL(m,dp)
           CALL next_image_cell_perd(m,x_data%periodic_parameter%perd)
           i=i+1
        ENDDO
      END IF
    CASE DEFAULT
      total_number_of_cells = 0
      IF( pbc_shells == -1 ) pbc_shells = 0
      DO i = 0,pbc_shells
        total_number_of_cells = total_number_of_cells + count_cells_perd(i,x_data%periodic_parameter%perd)
      END DO
      DEALLOCATE(x_data%neighbor_cells,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ALLOCATE(x_data%neighbor_cells(total_number_of_cells), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      m = 0
      i = 1
      DO WHILE(SUM(m**2)<=pbc_shells)
         x_data%neighbor_cells(i)%cell = REAL(m,dp)
         CALL next_image_cell_perd(m,x_data%periodic_parameter%perd)
         i=i+1
      ENDDO
    END SELECT

    ! ** Transform into real coord
    DO i=1,SIZE(x_data%neighbor_cells)
      r = 0.0_dp
      x_data%neighbor_cells(i)%cell_r(:) = 0.0_dp
      s = x_data%neighbor_cells(i)%cell(:)
      CALL scaled_to_real(x_data%neighbor_cells(i)%cell_r,s,cell)
    END DO
    x_data%periodic_parameter%number_of_shells = pbc_shells

    R_max_stress = 0.0_dp
    DO i=1,SIZE(x_data%neighbor_cells)
      R_max_stress = MAX(R_max_stress, MAXVAL(ABS(x_data%neighbor_cells(i)%cell_r(:))))
    END DO
    R_max_stress = R_max_stress + ABS(MAXVAL(cell%hmat(:,:)))
    x_data%periodic_parameter%R_max_stress = R_max_stress

  END SUBROUTINE hfx_create_neighbor_cells

  ! performs a fuzzy check of being in a quadrilateral
! *****************************************************************************
!> \brief ...
!> \param A ...
!> \param B ...
!> \param C ...
!> \param D ...
!> \param P ...
!> \retval point_is_in_quadrilateral ...
! *****************************************************************************
  FUNCTION point_is_in_quadrilateral(A,B,C,D,P)
    REAL(dp)                                 :: A(3), B(3), C(3), D(3), P(3)
    LOGICAL                                  :: point_is_in_quadrilateral

    REAL(dp), PARAMETER :: fuzzy = 1000.0_dp*EPSILON(1.0_dp)

    REAL(dp)                                 :: dot00, dot01, dot02, dot11, &
                                                dot12, invDenom, u, v, v0(3), &
                                                v1(3), v2(3)

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
    invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom
    ! ** Check if point is in triangle
    IF ( (u >= 0-fuzzy) .AND.  (v >= 0-fuzzy) .AND. (u + v <= 1+fuzzy) ) THEN
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
    invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * invDenom
    v = (dot00 * dot12 - dot01 * dot02) * invDenom

    ! ** Check if point is in triangle
    IF ( (u >= 0-fuzzy) .AND.  (v >= 0-fuzzy) .AND. (u + v <= 1+fuzzy) ) THEN
      point_is_in_quadrilateral = .TRUE.
      RETURN
    END IF

  END FUNCTION point_is_in_quadrilateral


! *****************************************************************************
!> \brief - This routine deletes all list entries in a container in order to
!>        deallocate the memory.
!> \param container container that contains the compressed elements
!> \param memory_usage ...
!> \param do_disk_storage ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_init_container(container, memory_usage, do_disk_storage, error)
    TYPE(hfx_container_type)                 :: container
    INTEGER                                  :: memory_usage
    LOGICAL                                  :: do_disk_storage
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_init_container', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(hfx_container_node), POINTER        :: current, next

    failure =.FALSE.

    !! DEALLOCATE memory

    current => container%first
    DO WHILE(ASSOCIATED(current))
      next => current%next
      DEALLOCATE(current)
      current=>next
    END DO

    !! Allocate first list entry, init members
    ALLOCATE(container%first,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    container%first%prev => NULL()
    container%first%next => NULL()
    container%current => container%first
    container%current%data = 0
    container%element_counter = 1
    memory_usage = 1

    IF( do_disk_storage) THEN
      !! close the file, if this is no the first time
      IF(container%unit /= -1 ) THEN
        CALL close_file(unit_number=container%unit)
      END IF
      CALL open_file(file_name=TRIM(container%filename),file_status="UNKNOWN",file_form="UNFORMATTED",file_action="WRITE",&
                     unit_number=container%unit)
    END IF

  END SUBROUTINE hfx_init_container


! *****************************************************************************
!> \brief - This routine stores the data obtained from the load balance routine
!>        for the energy
!> \param ptr_to_distr contains data to store
!> \param x_data contains all relevant data structures for hfx runs
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_set_distr_energy(ptr_to_distr,x_data, error)
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_distr
    TYPE(hfx_type), POINTER                  :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_set_distr_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    DEALLOCATE(x_data%distribution_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%distribution_energy(SIZE(ptr_to_distr)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    x_data%distribution_energy = ptr_to_distr

  END SUBROUTINE hfx_set_distr_energy

! *****************************************************************************
!> \brief - This routine stores the data obtained from the load balance routine
!>        for the forces
!> \param ptr_to_distr contains data to store
!> \param x_data contains all relevant data structures for hfx runs
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_set_distr_forces(ptr_to_distr, x_data, error)
    TYPE(hfx_distribution), DIMENSION(:), &
      POINTER                                :: ptr_to_distr
    TYPE(hfx_type), POINTER                  :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_set_distr_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    DEALLOCATE(x_data%distribution_forces,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%distribution_forces(SIZE(ptr_to_distr)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    x_data%distribution_forces = ptr_to_distr

  END SUBROUTINE hfx_set_distr_forces

! *****************************************************************************
!> \brief - resets the maximum memory usage for a HFX calculation subtracting
!>          all relevant buffers from the input MAX_MEM value and add 10% of
!>          safety margin
!> \param memory_parameter Memory information
!> \param subtr_size_mb size of buffers in MiB
!> \par History
!>      02.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_reset_memory_usage_counter(memory_parameter, subtr_size_mb)

    TYPE(hfx_memory_type)                    :: memory_parameter
    INTEGER(int_8), INTENT(IN)               :: subtr_size_mb

    INTEGER(int_8)                           :: max_memory

    max_memory = memory_parameter%max_memory
    max_memory = max_memory - subtr_size_mb
    IF( max_memory <= 0 ) THEN
      memory_parameter%do_all_on_the_fly = .TRUE.
      memory_parameter%max_compression_counter = 0
    ELSE
      memory_parameter%do_all_on_the_fly = .FALSE.
      memory_parameter%max_compression_counter = max_memory*1024_int_8*128_int_8
    END IF
  END SUBROUTINE hfx_reset_memory_usage_counter

! *****************************************************************************
!> \brief - This routine prints some information on HFX
!> \param x_data contains all relevant data structures for hfx runs
!> \param hfx_section HFX input section
!> \param i_rep current replica ID
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_print_info(x_data, hfx_section, i_rep, error)
    TYPE(hfx_type), POINTER                  :: x_data
    TYPE(section_vals_type), POINTER         :: hfx_section
    INTEGER, INTENT(IN)                      :: i_rep
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: iw
    REAL(dp)                                 :: rc_ang
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    iw = cp_print_key_unit_nr(logger,hfx_section,"HF_INFO",&
         extension=".scfLog",error=error)
    IF (iw>0) THEN
      WRITE (UNIT=iw,FMT="(/,(T3,A,T61,I20))")&
            "HFX_INFO| Replica ID:     ",i_rep
      WRITE (UNIT=iw,FMT="((T3,A,T73,ES8.1))")&
            "HFX_INFO| EPS_SCHWARZ:     ", x_data%screening_parameter%eps_schwarz
      WRITE (UNIT=iw,FMT="((T3,A,T73,ES8.1))")&
            "HFX_INFO| EPS_SCHWARZ_FORCES     ", x_data%screening_parameter%eps_schwarz_forces
      WRITE (UNIT=iw,FMT="((T3,A,T73,ES8.1))")&
            "HFX_INFO| EPS_STORAGE_SCALING:     ", x_data%memory_parameter%eps_storage_scaling
      WRITE (UNIT=iw,FMT="((T3,A,T61,I20))")&
            "HFX_INFO| NBINS:     ", x_data%load_balance_parameter%nbins
      WRITE (UNIT=iw,FMT="((T3,A,T61,I20))")&
            "HFX_INFO| BLOCK_SIZE:     ", x_data%load_balance_parameter%block_size
      WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| FRACTION:     ", x_data%general_parameter%fraction
      SELECT CASE (x_data%potential_parameter%potential_type)
        CASE(do_hfx_potential_coulomb)
          WRITE (UNIT=iw,FMT="((T3,A,T74,A))")&
            "HFX_INFO| Interaction Potential:     ", "COULOMB"
        CASE(do_hfx_potential_short)
          WRITE (UNIT=iw,FMT="((T3,A,T71,A))")&
            "HFX_INFO| Interaction Potential:    ", "SHORTRANGE"
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
          rc_ang = cp_unit_from_cp2k(x_data%potential_parameter%cutoff_radius,"angstrom",error=error)
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| Cutoff Radius [angstrom]:     ", rc_ang
        CASE(do_hfx_potential_long)
          WRITE (UNIT=iw,FMT="((T3,A,T72,A))")&
            "HFX_INFO| Interaction Potential:     ", "LONGRANGE"
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
        CASE(do_hfx_potential_mix_cl)
          WRITE (UNIT=iw,FMT="((T3,A,T75,A))")&
            "HFX_INFO| Interaction Potential:     ", "MIX_CL"
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| SCALE_COULOMB:     ", x_data%potential_parameter%scale_coulomb
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| SCALE_LONGRANGE:     ", x_data%potential_parameter%scale_longrange
        CASE(do_hfx_potential_gaussian)
          WRITE (UNIT=iw,FMT="((T3,A,T73,A))")&
            "HFX_INFO| Interaction Potential:     ", "GAUSSIAN"
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
        CASE(do_hfx_potential_mix_lg)
          WRITE (UNIT=iw,FMT="((T3,A,T75,A))")&
            "HFX_INFO| Interaction Potential:    ", "MIX_LG"
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| Omega:     ", x_data%potential_parameter%omega
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| SCALE_LONGRANGE:     ", x_data%potential_parameter%scale_longrange
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| SCALE_GAUSSIAN:    ", x_data%potential_parameter%scale_gaussian
        CASE(do_hfx_potential_id)
          WRITE (UNIT=iw,FMT="((T3,A,T73,A))")&
            "HFX_INFO| Interaction Potential:    ", "IDENTITY"
        CASE(do_hfx_potential_truncated)
          WRITE (UNIT=iw,FMT="((T3,A,T72,A))")&
            "HFX_INFO| Interaction Potential:    ", "TRUNCATED"
          rc_ang = cp_unit_from_cp2k(x_data%potential_parameter%cutoff_radius,"angstrom",error=error)
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| Cutoff Radius [angstrom]:     ", rc_ang
        CASE(do_hfx_potential_mix_cl_trunc)
          WRITE (UNIT=iw,FMT="((T3,A,T65,A))")&
            "HFX_INFO| Interaction Potential:    ", "TRUNCATED MIX_CL"
          rc_ang = cp_unit_from_cp2k(x_data%potential_parameter%cutoff_radius,"angstrom",error=error)
          WRITE(iw,'(T3,A,T61,F20.10)') &
            "HFX_INFO| Cutoff Radius [angstrom]:     ", rc_ang
      END SELECT
      IF( x_data%periodic_parameter%do_periodic) THEN
        IF( x_data%periodic_parameter%mode == -1 ) THEN
          WRITE (UNIT=iw,FMT="((T3,A,T77,A))")&
              "HFX_INFO| NUMBER_OF_SHELLS:     ", "AUTO"
        ELSE
          WRITE (UNIT=iw,FMT="((T3,A,T61,I20))")&
            "HFX_INFO| NUMBER_OF_SHELLS:     ", x_data%periodic_parameter%mode
        END IF
        WRITE (UNIT=iw,FMT="((T3,A,T61,I20))")&
            "HFX_INFO| Number of periodic shells considered:     ", x_data%periodic_parameter%number_of_shells
        WRITE (UNIT=iw,FMT="((T3,A,T61,I20),/)")&
              "HFX_INFO| Number of periodic cells considered:     ",SIZE(x_data%neighbor_cells)
      ELSE
        WRITE (UNIT=iw,FMT="((T3,A,T77,A))")&
            "HFX_INFO| Number of periodic shells considered:     ", "NONE"
        WRITE (UNIT=iw,FMT="((T3,A,T77,A),/)")&
              "HFX_INFO| Number of periodic cells considered:     ", "NONE"
      END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                     "HF_INFO", error=error)
  END SUBROUTINE hfx_print_info

! *****************************************************************************
!> \brief ...
!> \param x_data ...
!> \param eval_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dealloc_containers(x_data, eval_type, error)
    TYPE(hfx_type), POINTER                  :: x_data
    INTEGER, INTENT(IN)                      :: eval_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dealloc_containers', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bin, i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.


    SELECT CASE(eval_type)
    CASE (hfx_do_eval_energy)
      DO bin=1,SIZE(x_data%maxval_container)
        CALL hfx_init_container(x_data%maxval_container(bin), x_data%memory_parameter%actual_memory_usage, &
                                 .FALSE., error)
        DEALLOCATE(x_data%maxval_container(bin)%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO
      DEALLOCATE(x_data%maxval_container,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(x_data%maxval_cache,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO bin = 1,SIZE(x_data%integral_containers,2)
        DO i=1,64
          CALL hfx_init_container(x_data%integral_containers(i,bin), x_data%memory_parameter%actual_memory_usage, &
                                   .FALSE., error)
          DEALLOCATE(x_data%integral_containers(i,bin)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END DO
      END DO
      DEALLOCATE(x_data%integral_containers,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DEALLOCATE(x_data%integral_caches, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CASE (hfx_do_eval_forces)
      DO bin=1,SIZE(x_data%maxval_container_forces)
        CALL hfx_init_container(x_data%maxval_container_forces(bin), x_data%memory_parameter%actual_memory_usage, &
                                 .FALSE., error)
        DEALLOCATE(x_data%maxval_container_forces(bin)%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END DO
      DEALLOCATE(x_data%maxval_container_forces,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(x_data%maxval_cache_forces,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO bin = 1,SIZE(x_data%integral_containers_forces,2)
        DO i=1,64
          CALL hfx_init_container(x_data%integral_containers_forces(i,bin), x_data%memory_parameter%actual_memory_usage, &
                                   .FALSE., error)
          DEALLOCATE(x_data%integral_containers_forces(i,bin)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END DO
      END DO
      DEALLOCATE(x_data%integral_containers_forces,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DEALLOCATE(x_data%integral_caches_forces, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    END SELECT

  END SUBROUTINE dealloc_containers

! *****************************************************************************
!> \brief ...
!> \param x_data ...
!> \param bin_size ...
!> \param eval_type ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE alloc_containers(x_data, bin_size, eval_type, error)
    TYPE(hfx_type), POINTER                  :: x_data
    INTEGER, INTENT(IN)                      :: bin_size, eval_type
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'alloc_containers', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bin, i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.


    SELECT CASE(eval_type)
    CASE (hfx_do_eval_energy)
      ALLOCATE(x_data%maxval_cache(bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO bin = 1,bin_size
        x_data%maxval_cache(bin)%element_counter = 1
      END DO
      ALLOCATE(x_data%maxval_container(bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO bin = 1,bin_size
        ALLOCATE(x_data%maxval_container(bin)%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        x_data%maxval_container(bin)%first%prev => NULL()
        x_data%maxval_container(bin)%first%next => NULL()
        x_data%maxval_container(bin)%current => x_data%maxval_container(bin)%first
        x_data%maxval_container(bin)%current%data = 0
        x_data%maxval_container(bin)%element_counter = 1
      END DO

      ALLOCATE(x_data%integral_containers(64,bin_size), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(x_data%integral_caches(64,bin_size), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO bin = 1,bin_size
        DO i=1,64
          x_data%integral_caches(i,bin)%element_counter = 1
          x_data%integral_caches(i,bin)%data = 0
          ALLOCATE(x_data%integral_containers(i,bin)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          x_data%integral_containers(i,bin)%first%prev => NULL()
          x_data%integral_containers(i,bin)%first%next => NULL()
          x_data%integral_containers(i,bin)%current => x_data%integral_containers(i,bin)%first
          x_data%integral_containers(i,bin)%current%data = 0
          x_data%integral_containers(i,bin)%element_counter = 1
        END DO
      END DO
    CASE (hfx_do_eval_forces)
      ALLOCATE(x_data%maxval_cache_forces(bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO bin = 1,bin_size
        x_data%maxval_cache_forces(bin)%element_counter = 1
      END DO
      ALLOCATE(x_data%maxval_container_forces(bin_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO bin = 1,bin_size
        ALLOCATE(x_data%maxval_container_forces(bin)%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        x_data%maxval_container_forces(bin)%first%prev => NULL()
        x_data%maxval_container_forces(bin)%first%next => NULL()
        x_data%maxval_container_forces(bin)%current => x_data%maxval_container_forces(bin)%first
        x_data%maxval_container_forces(bin)%current%data = 0
        x_data%maxval_container_forces(bin)%element_counter = 1
      END DO

      ALLOCATE(x_data%integral_containers_forces(64,bin_size), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(x_data%integral_caches_forces(64,bin_size), STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO bin = 1,bin_size
        DO i=1,64
          x_data%integral_caches_forces(i,bin)%element_counter = 1
          x_data%integral_caches_forces(i,bin)%data = 0
          ALLOCATE(x_data%integral_containers_forces(i,bin)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          x_data%integral_containers_forces(i,bin)%first%prev => NULL()
          x_data%integral_containers_forces(i,bin)%first%next => NULL()
          x_data%integral_containers_forces(i,bin)%current => x_data%integral_containers_forces(i,bin)%first
          x_data%integral_containers_forces(i,bin)%current%data = 0
          x_data%integral_containers_forces(i,bin)%element_counter = 1
        END DO
      END DO

    END SELECT

  END SUBROUTINE alloc_containers

! *****************************************************************************
!> \brief ...
!> \param cli ...
! *****************************************************************************
  SUBROUTINE set_cell_info(cli)
    TYPE(hfx_cell_loop_info_type)            :: cli

    INTEGER                                  :: i1, i2, i3, id

    cli%index_to_cell = 0
    DO i3=LBOUND(cli%cell_to_index,3),UBOUND(cli%cell_to_index,3)
      DO i2=LBOUND(cli%cell_to_index,2),UBOUND(cli%cell_to_index,2)
        DO i1=LBOUND(cli%cell_to_index,1),UBOUND(cli%cell_to_index,1)
          id = cli%cell_to_index(i1,i2,i3)
          IF(id > 0 .AND. id <= SIZE(cli%index_to_cell,2)) THEN
            cli%index_to_cell(1,id) = i1
            cli%index_to_cell(2,id) = i2
            cli%index_to_cell(3,id) = i3
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE set_cell_info

! *****************************************************************************
!> \brief ...
!> \param cti ...
!> \param i1 ...
!> \param i2 ...
!> \param i3 ...
!> \retval cindex ...
! *****************************************************************************
  FUNCTION get_cell_index(cti,i1,i2,i3) RESULT(cindex)
    INTEGER, DIMENSION(:, :, :), POINTER     :: cti
    INTEGER                                  :: i1, i2, i3, cindex

    INTEGER                                  :: nl1, nl2, nl3, nu1, nu2, nu3

    cindex = 0

    nl1 = LBOUND(cti,1)
    nl2 = LBOUND(cti,2)
    nl3 = LBOUND(cti,3)
    nu1 = UBOUND(cti,1)
    nu2 = UBOUND(cti,2)
    nu3 = UBOUND(cti,3)

    IF(i1<nl1 .OR. i1>nu1) RETURN
    IF(i2<nl2 .OR. i2>nu2) RETURN
    IF(i3<nl3 .OR. i3>nu3) RETURN
   
    cindex = cti(i1,i2,i3)

  END FUNCTION get_cell_index

END MODULE hfx_types

