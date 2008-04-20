!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types and set/get functions for HFX
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE hfx_helpers,                     ONLY: count_cells,&
                                             next_image_cell
  USE hfx_libint_wrapper,              ONLY: lib_int
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp,&
                                             int_8
  USE machine,                         ONLY: m_chdir,&
                                             m_getcwd
  USE particle_types,                  ONLY: particle_type
  USE string_utilities,                ONLY: compress
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  hfx_type, hfx_create, hfx_release,& 
          hfx_set_distr_energy, &
          hfx_set_distr_forces, &
          hfx_cell_type, hfx_distribution, hfx_box_scheme,&
          hfx_potential_type, hfx_screening_type, hfx_periodic_type,&
          hfx_memory_type, hfx_load_balance_type, hfx_general_type,&
          hfx_compression_type,&
          hfx_container_type, hfx_container_node, hfx_cache_type, &
          hfx_basis_type, &
          hfx_p_type, hfx_create_neighbor_cells, hfx_init_container

#define CACHE_SIZE 1024

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_types'
  INTEGER(int_8), PARAMETER, PRIVATE :: one = 1_int_8
!***

! *****************************************************************************
  TYPE hfx_potential_type
    INTEGER                                  :: potential_type     !! 1/r/ erfc(wr)/r ...
    REAL(dp)                                 :: omega              !! w 
    REAL(dp)                                 :: scale_coulomb      !! scaling factor for mixed potential
    REAL(dp)                                 :: scale_longrange    !! scaling factor for mixed potential
    REAL(dp)                                 :: scale_gaussian     !! scaling factor for mixed potential
  END TYPE

! *****************************************************************************
  TYPE hfx_screening_type
    REAL(dp)                                 :: eps_schwarz         !! threshold
    REAL(dp)                                 :: eps_schwarz_forces  !! threshold
    REAL(dp)                                 :: eps_box             !! threshold
    REAL(dp)                                 :: eps_box_forces      !! threshold
    LOGICAL                                  :: do_p_screening_forces  !! screen on P^2 ?
    LOGICAL                                  :: do_initial_p_screening !! screen on inital guess?
    LOGICAL                                  :: do_init_p_screen_corr  !! correct screen on initial guess?
    REAL(dp)                                 :: initial_p_screening_correction !! corection factor
  END TYPE
  
! *****************************************************************************
  TYPE hfx_memory_type
    INTEGER                                  :: max_memory              !! user def max memory MB
    INTEGER(int_8)                           :: max_compression_counter !! corresponding number of reals 
    LOGICAL                                  :: do_all_on_the_fly       !! max mem == 0 ?
    REAL(dp)                                 :: eps_storage_scaling
    INTEGER                                  :: cache_size
    INTEGER                                  :: actual_memory_usage
    INTEGER                                  :: actual_memory_usage_disk
    INTEGER(int_8)                           :: max_compression_counter_disk
    LOGICAL                                  :: do_disk_storage
    CHARACTER(len=default_path_length)       :: storage_location
    INTEGER(int_8)                           :: ram_counter
  END TYPE

! *****************************************************************************
  TYPE hfx_periodic_type
    INTEGER                                  :: number_of_shells        !! number of periodic image cells 
    LOGICAL                                  :: do_periodic             !! periodic ?
  END TYPE

! *****************************************************************************
  TYPE hfx_load_balance_type
    INTEGER                                  :: nbins                   
    INTEGER                                  :: moves
    INTEGER                                  :: cycles
    INTEGER                                  :: reduction_steps
    REAL(dp)                                 :: reduction_factor
    REAL(dp)                                 :: temperature
    REAL(dp)                                 :: tolerance
    REAL(dp)                                 :: swap_probability
  END TYPE

! *****************************************************************************
  TYPE hfx_general_type
    REAL(dp)                                 :: fraction                 !! for hybrids
  END TYPE

! *****************************************************************************
  TYPE hfx_cell_type
    INTEGER                                  :: cell(3)
  END TYPE

! *****************************************************************************
  TYPE hfx_distribution
    INTEGER                                  :: istart, jstart, kstart, lstart
    INTEGER(int_8)                           :: number_of_integrals
    INTEGER(int_8)                           :: cost
    INTEGER(int_8)                           :: original_cpu_id
    REAL(KIND=dp)                            :: time_first_scf
    REAL(KIND=dp)                            :: time_other_scf
    REAL(KIND=dp)                            :: time_forces
  END TYPE
  
! *****************************************************************************
  TYPE hfx_box_scheme
    REAL(dp)                                 :: center(3)
    REAL(dp)                                 :: K
    REAL(dp)                                 :: R
    REAL(dp)                                 :: alphaInv
    INTEGER                                  :: radius                                 
  END TYPE
  
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
!> \param x_screening_buffer contains screening matrices in non-periodic case
!> \param x_screening_buffer_periodic contains screening matrices in periodic case
!> \param far_field_boxes contains information on far field screening
!> \param neighbor_cells manages handling of periodic cells
!> \param distribution_energy stores information on parallelization of energy
!> \param distribution_forces stores information on parallelization of forces
!> \param initial_p stores the initial guess if requested
!> \par History
!>      11.2006 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  TYPE hfx_type
    TYPE(hfx_potential_type)                 :: potential_parameter
    TYPE(hfx_screening_type)                 :: screening_parameter
    TYPE(hfx_memory_type)                    :: memory_parameter
    TYPE(hfx_periodic_type)                  :: periodic_parameter
    TYPE(hfx_load_balance_type)              :: load_balance_parameter
    TYPE(hfx_general_type)                   :: general_parameter
    TYPE(hfx_container_type), POINTER        :: maxval_container, maxval_container_disk
    TYPE(hfx_cache_type)                     :: maxval_cache, maxval_cache_disk
    TYPE(hfx_container_type), POINTER, &
                               DIMENSION(:)  :: integral_containers, integral_containers_disk
    TYPE(hfx_cache_type)                     :: integral_caches(64), integral_caches_disk(64)

    REAL(dp), DIMENSION(:, :, :, :), POINTER :: x_screening_buffer
    REAL(dp), DIMENSION(:, :, :, :, :), &
              POINTER                        :: x_screening_buffer_periodic
    TYPE(hfx_box_scheme), DIMENSION(:,:),&
                           POINTER           :: far_field_boxes
    TYPE(hfx_cell_type), DIMENSION(:), &
               POINTER                       :: neighbor_cells
    TYPE(hfx_distribution), DIMENSION(:),&
                             POINTER         :: distribution_energy
    TYPE(hfx_distribution), DIMENSION(:),&
                             POINTER         :: distribution_forces
    REAL(dp), DIMENSION(:,:,:,:), POINTER &
                                             :: initial_p
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
  END TYPE

! *****************************************************************************
!> \brief to build arrays of pointers
!> \param hfx_data the pointer to the hfx_type
!> \par History
!>      02.2008 created [mguidon]
!> \author Manuel Guidon
! *****************************************************************************
  TYPE hfx_p_type
     TYPE(hfx_type), POINTER :: hfx_data
  END TYPE hfx_p_type


  CONTAINS

! *****************************************************************************
!> \brief - This routine allocates and initializes all types in hfx_data
!> \param x_data contains all relevant data structures for hfx runs
!> \param hfx_section input section
!> \param natom , max_set: quantities that define some array sizes
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - All POINTERS and ALLOCATABLES are allocated, even if their size is
!>        unknown at invocation time
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_create(x_data,para_env, hfx_section,natom,max_set,atomic_kind_set, particle_set, error)
    TYPE(hfx_type), DIMENSION(:, :), POINTER :: x_data
    TYPE(cp_para_env_type)                   :: para_env
    TYPE(section_vals_type), POINTER         :: hfx_section
    INTEGER, INTENT(IN)                      :: natom, max_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=512)                       :: error_msg
    CHARACTER(LEN=default_path_length)       :: char_val, filename, orig_wd
    INTEGER :: i, i_thread, ikind, int_val, irep, max_memory, n_rep_hf, &
      n_threads, nkind, pbc_shells, stat, storage_id
    LOGICAL                                  :: do_periodic, failure, &
                                                logic_val
    REAL(dp)                                 :: omega, real_val
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(section_vals_type), POINTER         :: hf_pbc_section, hf_sub_section

!$  INTEGER                                  :: omp_get_max_threads

    failure = .FALSE.
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
  
        !! MEMORY section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"MEMORY",i_rep_section=irep, error=error)
        CALL section_vals_val_get(hf_sub_section, "MAX_MEMORY", i_val=int_val, error=error)
        actual_x_data%memory_parameter%max_memory = int_val
        actual_x_data%memory_parameter%max_compression_counter = int_val*1024*128
        CALL section_vals_val_get(hf_sub_section, "EPS_STORAGE_SCALING", r_val=real_val, error=error)
        actual_x_data%memory_parameter%eps_storage_scaling = real_val
        IF(int_val == 0) THEN
          actual_x_data%memory_parameter%do_all_on_the_fly = .TRUE.
        ELSE
          actual_x_data%memory_parameter%do_all_on_the_fly = .FALSE.
        END IF
        actual_x_data%memory_parameter%cache_size = CACHE_SIZE
        actual_x_data%memory_parameter%actual_memory_usage = 1
        actual_x_data%memory_parameter%actual_memory_usage_disk = 1
        CALL section_vals_val_get(hf_sub_section, "MAX_DISK_SPACE", i_val=int_val, error=error)
!!        actual_x_data%memory_parameter%max_memory = int_val
        actual_x_data%memory_parameter%max_compression_counter_disk = int_val*1024_int_8*128_int_8
        IF( int_val == 0 ) THEN
          actual_x_data%memory_parameter%do_disk_storage = .FALSE.
        ELSE
          actual_x_data%memory_parameter%do_disk_storage = .TRUE.
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
          WRITE(error_msg,'(A,A,A)') "HFX_MODULE: Request for disk storage failed due to unknown error while writing to ",&
                                   TRIM(filename), ". Please check STORAGE_LOCATION"
          CALL cp_assert( .FALSE. , cp_failure_level,cp_assertion_failed,routineP,&
                          error_msg , &
                          error,failure)
        END IF
        CALL m_chdir(orig_wd,stat)

        actual_x_data%memory_parameter%storage_location = filename
        
        CALL compress(actual_x_data%memory_parameter%storage_location, .TRUE.)
        storage_id = (irep-1)*para_env%num_pe*n_threads + para_env%mepos * n_threads + i_thread-1

        !! PERIODIC section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"PERIODIC", i_rep_section=irep, error=error)
        CALL section_vals_get(hf_sub_section, explicit=logic_val, error=error)
        actual_x_data%periodic_parameter%do_periodic = logic_val
        CALL section_vals_val_get(hf_sub_section, "NUMBER_OF_SHELLS", i_val=int_val, error=error)
        actual_x_data%periodic_parameter%number_of_shells = int_val
        
        !! SCREENING section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"SCREENING",i_rep_section=irep, error=error)
        CALL section_vals_val_get(hf_sub_section, "EPS_SCHWARZ", r_val=real_val, error=error)
        actual_x_data%screening_parameter%eps_schwarz = real_val
        CALL section_vals_val_get(hf_sub_section, "EPS_BOX", r_val=real_val, error=error)
        actual_x_data%screening_parameter%eps_box = real_val 
        CALL section_vals_val_get(hf_sub_section, "EPS_SCHWARZ_FORCES", r_val=real_val, error=error)
        actual_x_data%screening_parameter%eps_schwarz_forces = real_val 
        CALL section_vals_val_get(hf_sub_section, "EPS_BOX_FORCES", r_val=real_val, error=error)
        actual_x_data%screening_parameter%eps_box_forces = real_val 
        CALL section_vals_val_get(hf_sub_section, "SCREEN_P_FORCES", l_val=logic_val, error=error)
        actual_x_data%screening_parameter%do_p_screening_forces = logic_val 
        CALL section_vals_val_get(hf_sub_section, "SCREEN_ON_INITIAL_P", l_val=logic_val, error=error)
        actual_x_data%screening_parameter%do_initial_p_screening = logic_val 
        CALL section_vals_val_get(hf_sub_section, "P_SCREEN_CORRECTION_FACTOR", r_val=real_val, error=error)
        IF( real_val == 0.0_dp ) THEN
          actual_x_data%screening_parameter%do_init_p_screen_corr = .FALSE. 
        ELSE
          actual_x_data%screening_parameter%do_init_p_screen_corr = .TRUE. 
        END IF
        actual_x_data%screening_parameter%initial_p_screening_correction = real_val 
    
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
    
        !! LOAD_BALANCE section
        hf_sub_section => section_vals_get_subs_vals(hfx_section,"LOAD_BALANCE",i_rep_section=irep, error=error)
        CALL section_vals_val_get(hf_sub_section, "NBINS", i_val=int_val, error=error)
        actual_x_data%load_balance_parameter%nbins = int_val
        CALL section_vals_val_get(hf_sub_section, "MOVES", i_val=int_val, error=error)
        actual_x_data%load_balance_parameter%moves = int_val
        CALL section_vals_val_get(hf_sub_section, "CYCLES", i_val=int_val, error=error)
        actual_x_data%load_balance_parameter%cycles = int_val
        CALL section_vals_val_get(hf_sub_section, "REDUCTION_STEPS", i_val=int_val, error=error)
        actual_x_data%load_balance_parameter%reduction_steps = int_val
        CALL section_vals_val_get(hf_sub_section, "REDUCTION_FACTOR", r_val=real_val, error=error)
        actual_x_data%load_balance_parameter%reduction_factor = real_val
        CALL section_vals_val_get(hf_sub_section, "TEMPERATURE", r_val=real_val, error=error)
        actual_x_data%load_balance_parameter%temperature = real_val
        CALL section_vals_val_get(hf_sub_section, "TOLERANCE", r_val=real_val, error=error)
        actual_x_data%load_balance_parameter%tolerance = real_val
        CALL section_vals_val_get(hf_sub_section, "SWAP_PROBABILITY", r_val=real_val, error=error)
        actual_x_data%load_balance_parameter%swap_probability = real_val
    
        !! BASIS parameter
        nkind = SIZE(atomic_kind_set,1)
        ALLOCATE(actual_x_data%basis_parameter(nkind), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ikind = 1,nkind
            NULLIFY(atom_kind)
            atom_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atom_kind,&
                                 orb_basis_set=orb_basis_a)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_a, &
                                   lmax=actual_x_data%basis_parameter(ikind)%lmax, &
                                   lmin=actual_x_data%basis_parameter(ikind)%lmin, &
                                   npgf=actual_x_data%basis_parameter(ikind)%npgf, &
                                   nset=actual_x_data%basis_parameter(ikind)%nset, &
                                   zet=actual_x_data%basis_parameter(ikind)%zet, &
                                   nsgf_set=actual_x_data%basis_parameter(ikind)%nsgf, &
                                   first_sgf=actual_x_data%basis_parameter(ikind)%first_sgf, &
                                   sphi=actual_x_data%basis_parameter(ikind)%sphi, &
                                   nsgf=actual_x_data%basis_parameter(ikind)%nsgf_total) 
        END DO
    
        IF(actual_x_data%periodic_parameter%do_periodic) THEN
          hf_pbc_section => section_vals_get_subs_vals(hfx_section,"PERIODIC",i_rep_section=irep, error=error)
          CALL section_vals_val_get(hf_pbc_section,"NUMBER_OF_SHELLS",i_val=pbc_shells,error=error)
          ALLOCATE(actual_x_data%neighbor_cells(1),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL hfx_create_neighbor_cells(actual_x_data,pbc_shells, error=error)
          ALLOCATE(actual_x_data%x_screening_buffer_periodic(SIZE(actual_x_data%neighbor_cells),&
                   max_set,max_set,natom,natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ELSE
          ALLOCATE(actual_x_data%x_screening_buffer(max_set,max_set,natom,natom),STAT=stat)
        END IF
        
        ALLOCATE(actual_x_data%is_assoc_atomic_block(natom,natom),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        ALLOCATE(actual_x_data%distribution_forces(1),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(actual_x_data%distribution_energy(1),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        ALLOCATE(actual_x_data%far_field_boxes(1,1),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        
        ALLOCATE(actual_x_data%initial_p(max_set,max_set,natom,natom),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        actual_x_data%maxval_cache%element_counter = 1
        ALLOCATE(actual_x_data%maxval_container,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(actual_x_data%maxval_container%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        actual_x_data%maxval_container%first%prev => NULL()
        actual_x_data%maxval_container%first%next => NULL()
        actual_x_data%maxval_container%current => actual_x_data%maxval_container%first
        actual_x_data%maxval_container%current%data = 0
        actual_x_data%maxval_container%element_counter = 1
    
        ALLOCATE(actual_x_data%integral_containers(64), STAT=stat)
        DO i=1,64
          actual_x_data%integral_caches(i)%element_counter = 1
          actual_x_data%integral_caches(i)%data = 0
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(actual_x_data%integral_containers(i)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          actual_x_data%integral_containers(i)%first%prev => NULL()
          actual_x_data%integral_containers(i)%first%next => NULL()
          actual_x_data%integral_containers(i)%current => actual_x_data%integral_containers(i)%first
          actual_x_data%integral_containers(i)%current%data = 0
          actual_x_data%integral_containers(i)%element_counter = 1
        END DO

        actual_x_data%maxval_cache_disk%element_counter = 1
        ALLOCATE(actual_x_data%maxval_container_disk,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(actual_x_data%maxval_container_disk%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        actual_x_data%maxval_container_disk%first%prev => NULL()
        actual_x_data%maxval_container_disk%first%next => NULL()
        actual_x_data%maxval_container_disk%current => actual_x_data%maxval_container%first
        actual_x_data%maxval_container_disk%current%data = 0
        actual_x_data%maxval_container_disk%element_counter = 1
        actual_x_data%maxval_container_disk%file_counter = 1
        actual_x_data%maxval_container_disk%desc = 'Max_'
        actual_x_data%maxval_container_disk%unit = -1
        WRITE(actual_x_data%maxval_container_disk%filename,'(A,I0,A,A,A)') TRIM(actual_x_data%memory_parameter%storage_location), &
                                                           storage_id,"_",actual_x_data%maxval_container_disk%desc, "6"
        CALL compress(actual_x_data%maxval_container_disk%filename, .TRUE.)
        ALLOCATE(actual_x_data%integral_containers_disk(64), STAT=stat)
        DO i=1,64
          actual_x_data%integral_caches_disk(i)%element_counter = 1
          actual_x_data%integral_caches_disk(i)%data = 0
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(actual_x_data%integral_containers_disk(i)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          actual_x_data%integral_containers_disk(i)%first%prev => NULL()
          actual_x_data%integral_containers_disk(i)%first%next => NULL()
          actual_x_data%integral_containers_disk(i)%current => actual_x_data%integral_containers(i)%first
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
  END SUBROUTINE hfx_create 

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

!$  INTEGER                                  :: omp_get_max_threads


    failure = .FALSE.
    !! There might be 2 hf sections
    n_rep_hf = x_data(1,1)%n_rep_hf
    n_threads = 1
!$  n_threads = omp_get_max_threads()
    DO i_thread= 1,n_threads
      DO irep=1,n_rep_hf
        actual_x_data => x_data(irep,i_thread)
        IF( actual_x_data%periodic_parameter%do_periodic ) THEN
          DEALLOCATE(actual_x_data%x_screening_buffer_periodic,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(actual_x_data%neighbor_cells,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ELSE
          DEALLOCATE(actual_x_data%x_screening_buffer,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF
        DEALLOCATE(actual_x_data%distribution_energy,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE(actual_x_data%distribution_forces,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        DEALLOCATE(actual_x_data%far_field_boxes,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        DEALLOCATE(actual_x_data%initial_p,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        DEALLOCATE(actual_x_data%is_assoc_atomic_block,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        !! BASIS parameter
        DEALLOCATE(actual_x_data%basis_parameter,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        !! Deallocate containers
        CALL hfx_init_container(6, actual_x_data%maxval_container, actual_x_data%memory_parameter%actual_memory_usage, &
                                 .FALSE., error)
        DEALLOCATE(actual_x_data%maxval_container%first,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE(actual_x_data%maxval_container,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
        DO i=1,64   
          CALL hfx_init_container(i, actual_x_data%integral_containers(i), actual_x_data%memory_parameter%actual_memory_usage, &
                                   .FALSE., error)
          DEALLOCATE(actual_x_data%integral_containers(i)%first,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END DO
        DEALLOCATE(actual_x_data%integral_containers,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        !! Deallocate containers
        CALL hfx_init_container(6, actual_x_data%maxval_container_disk, &
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
          CALL hfx_init_container(i, actual_x_data%integral_containers_disk(i), &
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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_create_neighbor_cells(x_data, pbc_shells, error)
    TYPE(hfx_type), POINTER                  :: x_data
    INTEGER, INTENT(IN)                      :: pbc_shells
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create_neighbor_cells', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, m(3), stat, &
                                                total_number_of_cells
    LOGICAL                                  :: failure

    failure = .FALSE.
    total_number_of_cells = 0
    DO i = 0,pbc_shells
      total_number_of_cells = total_number_of_cells + count_cells(i)
    END DO

    DEALLOCATE(x_data%neighbor_cells,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%neighbor_cells(total_number_of_cells), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    m = 0
    i = 1
    DO WHILE(SUM(m**2)<=pbc_shells)
       x_data%neighbor_cells(i)%cell = m
       CALL next_image_cell(m)
       i=i+1
    ENDDO
  END SUBROUTINE hfx_create_neighbor_cells

! *****************************************************************************
!> \brief - This routine deletes all list entries in a container in order to
!>        deallocate the memory.
!> \param container container that contains the compressed elements
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE hfx_init_container(nbits, container, memory_usage, do_disk_storage, error)
    INTEGER                                  :: nbits
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
!> \param ptr_to_distr_energy contains data to store
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

    ALLOCATE(x_data%distribution_energy(SIZE(ptr_to_distr)))
    x_data%distribution_energy = ptr_to_distr

  END SUBROUTINE hfx_set_distr_energy

! *****************************************************************************
!> \brief - This routine stores the data obtained from the load balance routine 
!>        for the forces
!> \param ptr_to_distr_energy contains data to store
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
   
    ALLOCATE(x_data%distribution_forces(SIZE(ptr_to_distr)))
    x_data%distribution_forces = ptr_to_distr

  END SUBROUTINE hfx_set_distr_forces

END MODULE hfx_types

